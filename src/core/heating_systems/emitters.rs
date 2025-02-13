/// This module provides objects to represent radiator and underfloor emitter systems.
use crate::compare_floats::{max_of_2, min_of_2};
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::heating_systems::common::SpaceHeatingService;
use crate::core::heating_systems::heat_pump::{
    BufferTankEmittersData, BufferTankEmittersDataWithResult,
};
use crate::core::material_properties::WATER;
use crate::core::solvers::{fsolve, root};
use crate::core::space_heat_demand::zone::SimpleZone;
use crate::core::units::{
    JOULES_PER_KILOJOULE, KILOJOULES_PER_KILOWATT_HOUR, LITRES_PER_CUBIC_METRE, WATTS_PER_KILOWATT,
};
use crate::corpus::KeyString;
use crate::external_conditions::ExternalConditions;
use crate::input::{
    EcoDesignController, EcoDesignControllerClass, FanSpeedData, FancoilTestData,
    WetEmitter as WetEmitterInput,
};
use crate::simulation_time::SimulationTimeIteration;
use crate::statistics::np_interp;
use crate::StringOrNumber;
use anyhow::{anyhow, bail, Error};
use argmin::core::{CostFunction, Executor};
use argmin::solver::brent::BrentRoot;
use atomic_float::AtomicF64;
use derivative::Derivative;
use itertools::Itertools;
use ode_solvers::{dop_shared::OutputType, Dopri5, System, Vector1};
use ordered_float::OrderedFloat;
use parking_lot::RwLock;
use std::fmt::Debug;
use std::ops::Deref;
use std::sync::atomic::Ordering;
use std::sync::Arc;

type State = Vector1<f64>;
type Time = f64;

/// Convert flow temperature to return temperature using the 6/7th rule.
///
/// Parameters:
///     `flow_temp_celsius` - Flow temperature in degrees Celsius.
///
///     Returns:
///     float: Return temperature in degrees Celsius.
pub fn convert_flow_to_return_temp(flow_temp_celsius: f64) -> f64 {
    (6.0 / 7.0) * flow_temp_celsius
}

#[derive(Derivative)]
#[derivative(Debug)]
pub(crate) struct Emitters {
    pub thermal_mass: f64,
    emitters: Vec<Arc<WetEmitter>>,
    temp_diff_emit_dsgn: f64,
    heat_source: Arc<RwLock<SpaceHeatingService>>,
    #[derivative(Debug = "ignore")]
    zone: Arc<dyn SimpleZone>,
    external_conditions: Arc<ExternalConditions>,
    with_buffer_tank: bool,
    variable_flow_data: VariableFlowData,
    bypass_percentage_recirculated: f64,
    design_flow_temp: f64,
    ecodesign_controller_class: EcoDesignControllerClass,
    min_outdoor_temp: Option<f64>,
    max_outdoor_temp: Option<f64>,
    min_flow_temp: Option<f64>,
    max_flow_temp: Option<f64>,
    temp_emitter_prev: AtomicF64,
    target_flow_temp: AtomicF64, // In Python this is set from inside demand energy and does not exist before then
    output_detailed_results: bool,
    emitters_detailed_results: Option<Arc<RwLock<Vec<EmittersDetailedResult>>>>,
    energy_supply_fan_coil_conn: Option<Arc<EnergySupplyConnection>>,
    min_flow_rate: f64,
    max_flow_rate: f64,
    fancoil: Option<Arc<WetEmitter>>,
}

#[derive(Copy, Clone)]
struct EmittersAndPowerInput<'a> {
    emitters: &'a Emitters,
    power_input: f64,
    temp_diff_max: Option<f64>,
    previous_difference_from_temp_diff_max: Option<f64>,
}

impl EmittersAndPowerInput<'_> {
    pub fn new(
        emitters: &Emitters,
        power_input: f64,
        temp_diff_max: Option<f64>,
        temp_diff_start: f64,
    ) -> EmittersAndPowerInput {
        let previous_difference_from_temp_diff_max = temp_diff_max.map(|max| temp_diff_start - max);

        EmittersAndPowerInput {
            emitters,
            power_input,
            temp_diff_max,
            previous_difference_from_temp_diff_max,
        }
    }

    fn difference_from_temp_diff_max(&self, y: f64) -> f64 {
        y - self.temp_diff_max.unwrap()
    }
}

// Here we're using the ode_solvers crate to replicate
// ODE solving functionality in scipy's solve_ivp
impl System<Time, State> for EmittersAndPowerInput<'_> {
    fn system(&self, _x: Time, y: &State, dy: &mut State) {
        dy[0] = self
            .emitters
            .func_temp_emitter_change_rate(self.power_input)([y[0]]);
    }

    // Stop function called at every successful integration step. The integration is stopped when this function returns true.
    fn solout(&mut self, _x: Time, y: &State, _dy: &State) -> bool {
        if self.temp_diff_max.is_none() {
            // no maximum - keep going
            return false;
        }

        let current_difference = self.difference_from_temp_diff_max(y[0]);

        if self.previous_difference_from_temp_diff_max.is_some() {
            let previous_difference = self.previous_difference_from_temp_diff_max.unwrap();

            // signs are different - we must have passed zero
            if current_difference == 0.
                || signs_are_different(current_difference, previous_difference)
                || previous_difference == 0.
            {
                // passing zero means we hit temp_diff_max, so stop solver
                return true;
            }
        }

        self.previous_difference_from_temp_diff_max = Some(current_difference);
        false
    }
}

fn signs_are_different(a: f64, b: f64) -> bool {
    (b > 0. && a < 0.) || (b < 0. && a > 0.)
}

// Here we're using argmin for root solving on our ode_solver `stepper`
// This is to replicate the `events` feature in scipy's solve_ivp
struct RootProblem<'a> {
    pub stepper: &'a Dopri5<
        f64,
        nalgebra::Matrix<
            f64,
            nalgebra::Const<1>,
            nalgebra::Const<1>,
            nalgebra::ArrayStorage<f64, 1, 1>,
        >,
        EmittersAndPowerInput<'a>,
    >,
    pub max_temp: f64,
}

impl CostFunction for RootProblem<'_> {
    type Param = f64;
    type Output = f64;

    fn cost(&self, x: &Self::Param) -> Result<Self::Output, Error> {
        // Difference between the (interpolated) temperature at time x
        // and the maximum temperature of the emitter
        let cost = self.stepper.dense_output_for_last_step(*x)[0];
        Ok(cost - self.max_temp)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum VariableFlowData {
    Yes,
    No { design_flow_rate: f64 },
}

impl From<VariableFlowData> for bool {
    fn from(v: VariableFlowData) -> Self {
        matches!(v, VariableFlowData::Yes)
    }
}

#[derive(Clone, Debug)]
enum WetEmitter {
    Radiator {
        c: f64,
        n: f64,
        frac_convective: f64,
    },
    Ufh {
        c: f64,
        n: f64,
        equivalent_specific_thermal_mass: f64,
        system_performance_factor: f64,
        emitter_floor_area: f64,
        frac_convective: f64,
    },
    Fancoil {
        n_units: usize,
        frac_convective: f64,
        test_data: FancoilTestData,
        temperature_data: Vec<Vec<f64>>,
        fan_power_data: Vec<StringOrNumber>,
    },
}

impl WetEmitter {
    fn frac_convective(&self) -> f64 {
        *match self {
            WetEmitter::Radiator {
                frac_convective, ..
            } => frac_convective,
            WetEmitter::Ufh {
                frac_convective, ..
            } => frac_convective,
            WetEmitter::Fancoil {
                frac_convective, ..
            } => frac_convective,
        }
    }

    fn c(&self) -> anyhow::Result<f64> {
        match self {
            WetEmitter::Radiator { c, .. } => Ok(*c),
            WetEmitter::Ufh { c, .. } => Ok(*c),
            WetEmitter::Fancoil { .. } => bail!("Fancoil emitters do not have a c value"),
        }
    }

    fn n(&self) -> anyhow::Result<f64> {
        match self {
            WetEmitter::Radiator { n, .. } => Ok(*n),
            WetEmitter::Ufh { n, .. } => Ok(*n),
            WetEmitter::Fancoil { .. } => bail!("Fancoil emitters do not have a n value"),
        }
    }
}

impl Emitters {
    /// Construct an Emitters object
    ///
    /// Arguments:
    /// * `thermal_mass` - thermal mass of emitters, in kWh / K
    /// * `c` - constant from characteristic equation of emitters (e.g. derived from BS EN 442 tests)
    /// * `n` - exponent from characteristic equation of emitters (e.g. derived from BS EN 442 tests)
    /// * `temp_diff_emit_dsgn` - design temperature difference across the emitters, in deg C or K
    /// * `frac_convective` - convective fraction for heating
    /// * `heat_source` - reference to an object representing the system (e.g.
    ///                       boiler or heat pump) providing heat to the emitters
    /// * `zone` - reference to the Zone object representing the zone in which the
    ///             emitters are located
    /// * `simulation_timestep` - timestep length for simulation time being used in this context
    ///
    /// Other variables:
    /// * `temp_emitter_prev` - temperature of the emitters at the end of the
    ///    previous timestep, in deg C
    pub(crate) fn new(
        thermal_mass: Option<f64>,
        emitters: &[WetEmitterInput],
        temp_diff_emit_dsgn: f64,
        variable_flow: bool,
        design_flow_rate: Option<f64>,
        min_flow_rate: Option<f64>,
        max_flow_rate: Option<f64>,
        bypass_percentage_recirculated: Option<f64>,
        heat_source: Arc<RwLock<SpaceHeatingService>>,
        zone: Arc<dyn SimpleZone>,
        external_conditions: Arc<ExternalConditions>,
        ecodesign_controller: EcoDesignController,
        design_flow_temp: f64,
        energy_supply_fan_coil_conn: Option<Arc<EnergySupplyConnection>>,
        output_detailed_results: bool,
        with_buffer_tank: bool,
    ) -> anyhow::Result<Self> {
        let (variable_flow_data, min_flow_rate, max_flow_rate) = if variable_flow {
            if let (Some(min_flow_rate), Some(max_flow_rate)) = (min_flow_rate, max_flow_rate) {
                (VariableFlowData::Yes, min_flow_rate, max_flow_rate)
            } else {
                bail!("Both min_flow_rate and max_flow_rate are required if variable_flow is true")
            }
        } else if let Some(design_flow_rate) = design_flow_rate {
            (
                VariableFlowData::No { design_flow_rate },
                design_flow_rate,
                design_flow_rate,
            )
        } else {
            bail!("design_flow_rate is required if variable_flow is false")
        };
        let bypass_percentage_recirculated = bypass_percentage_recirculated.unwrap_or(0.0);
        let ecodesign_controller_class = ecodesign_controller.ecodesign_control_class;
        let (min_outdoor_temp, max_outdoor_temp, min_flow_temp, max_flow_temp) = if matches!(
            ecodesign_controller_class,
            EcoDesignControllerClass::ClassII
                | EcoDesignControllerClass::ClassIII
                | EcoDesignControllerClass::ClassVI
                | EcoDesignControllerClass::ClassVII
        ) {
            (
                ecodesign_controller.min_outdoor_temp,
                ecodesign_controller.max_outdoor_temp,
                ecodesign_controller.min_flow_temp,
                Some(design_flow_temp),
            )
        } else {
            (None, None, None, None)
        };

        let floor_area = zone.area();
        let mut total_emitter_floor_area: f64 = Default::default();
        let mut fancoil: Option<Arc<WetEmitter>> = Default::default();

        let mut model_emitters: Vec<Option<Arc<WetEmitter>>> = vec![None; emitters.len()];

        // 1. process radiators first
        for (i, emitter) in emitters.iter().enumerate() {
            if let WetEmitterInput::Radiator {
                c,
                n,
                frac_convective,
            } = emitter
            {
                model_emitters[i].replace(Arc::new(WetEmitter::Radiator {
                    c: *c,
                    n: *n,
                    frac_convective: *frac_convective,
                }));
                if thermal_mass.is_none() {
                    // Thermal_mass is a required input for radiators -
                    // not underfloor. This is because the thermal mass of UFH is
                    // included in the UFH-only 'equivalent_specific_thermal_mass' input.
                    // The latter can only be calculated for UFH systems, by definition,
                    // so could not mistakenly be entered for a radiator system.
                    // But for a system containing mix of radiators and UFH, the thermal_mass
                    // input is required - including only the thermal mass of the radiators.
                    bail!("Thermal Mass is required for Radiator type emitters")
                }
            }
        }

        let mut thermal_mass = thermal_mass.unwrap_or(0.0);

        // 2. process UFH emitters next
        for (i, emitter) in emitters.iter().enumerate() {
            if let WetEmitterInput::Ufh {
                equivalent_specific_thermal_mass,
                system_performance_factor,
                emitter_floor_area,
                frac_convective,
            } = emitter
            {
                let n = 1.; // For UFH, BS EN 1264 and 11855 define this as 1 under normal circumstances
                let c = system_performance_factor * emitter_floor_area / WATTS_PER_KILOWATT as f64;
                total_emitter_floor_area += emitter_floor_area;

                // The thermal_mass input from assessor only includes radiators.
                // The equivalent_specific_thermal_mass for UFH (once converted to
                // the same units) needs to be added to this to get the total.
                thermal_mass += equivalent_specific_thermal_mass * emitter_floor_area
                    / KILOJOULES_PER_KILOWATT_HOUR as f64;
                model_emitters[i].replace(Arc::new(WetEmitter::Ufh {
                    c,
                    n,
                    equivalent_specific_thermal_mass: *equivalent_specific_thermal_mass,
                    system_performance_factor: *system_performance_factor,
                    emitter_floor_area: *emitter_floor_area,
                    frac_convective: *frac_convective,
                }));
            }
        }

        // 3. process fancoils last
        for (i, emitter) in emitters.iter().enumerate() {
            if let WetEmitterInput::Fancoil {
                n_units,
                frac_convective,
                fancoil_test_data,
            } = emitter
            {
                let n_units = n_units.unwrap_or(1);
                let test_data = Self::add_temperature_diff_zero(fancoil_test_data);
                let (temperature_data, fan_power_data) =
                    Self::format_fancoil_manufacturer_data(fancoil_test_data)?;

                let emitter = Arc::new(WetEmitter::Fancoil {
                    n_units,
                    frac_convective: *frac_convective,
                    test_data,
                    temperature_data,
                    fan_power_data,
                });
                model_emitters[i].replace(emitter.clone());

                // Only one specification in initial implementation
                fancoil.replace(emitter);
            }
        }

        // TODO (from Python) Calculate convective fraction for UFH from floor surface temperature Tf,
        // and the room air temperature, according the formula below.
        // Ta = self.__zone.temp_internal_air()  # room_air_temp
        // self.__frac_convective = ((8.92 * (Tf - Ta) ** 1.1 / (Tf - Ta)) - 5.5) / (8.92 * ( Tf - Ta ) ** 1.1 / ( Tf - Ta ))
        // Need to come up with a method to calculate floor surface temperature.

        // Final initialisation checks:

        // Ensure total UFH area does not exceed zone area
        if total_emitter_floor_area > floor_area {
            bail!(
                "Total UFH area {total_emitter_floor_area} is bigger than Zone area {floor_area})"
            );
        }

        let model_emitters: Vec<Arc<WetEmitter>> = model_emitters
            .into_iter()
            .map(|emitter| emitter.ok_or_else(|| anyhow!("Emitter type unexpectedly not ingested")))
            .collect::<anyhow::Result<_>>()?;

        // Considering the big differences in the calculation for c, n type emitters (like Radiators and UFH)
        // and fancoils that use a implicit (manufacturer data driven) approach, the initial implementation
        // of several types of emitters does not allow for the mix of fancoils with other systems. This is
        // currently enforced through using flag_fancoil.
        // Ensure only one fancoil is defined and it is the sole emitter type
        if fancoil.is_some() && model_emitters.len() > 1 {
            bail!("Only one fancoil specification can be defined, and it must be the sole emitter type for the zone");
        }

        Ok(Self {
            thermal_mass,
            emitters: model_emitters,
            temp_diff_emit_dsgn,
            heat_source,
            zone,
            external_conditions,
            with_buffer_tank,
            variable_flow_data,
            bypass_percentage_recirculated,
            design_flow_temp,
            ecodesign_controller_class,
            min_outdoor_temp,
            max_outdoor_temp,
            min_flow_temp,
            max_flow_temp,
            temp_emitter_prev: 20.0.into(),
            target_flow_temp: 20.0.into(), // initial value, though expected to be updated before being used
            output_detailed_results,
            emitters_detailed_results: output_detailed_results.then(Default::default),
            energy_supply_fan_coil_conn,
            min_flow_rate,
            max_flow_rate,
            fancoil,
        })
    }

    /// This function appends the product data with a row for delta_T = 0.0, if missing.
    /// At a delta_T = 0, the heat output is assumed to be equal to the fan power only.
    /// NB. The upstream Python version of this function modifies and returns the emitter object.
    /// This method, however, creates a new test data value to be added to the emitter that is part
    /// of the model.
    fn add_temperature_diff_zero(fancoil_test_data: &FancoilTestData) -> FancoilTestData {
        // Extract the fan speed data and fan power for the emitter
        let mut fan_speed_data = fancoil_test_data.fan_speed_data.clone();
        let fan_power_w = fancoil_test_data.fan_power_w.clone();

        // Check if temperature_diff 0.0 exists in the fan_speed_data
        if !fan_speed_data
            .iter()
            .any(|item| item.temperature_diff == 0.0)
        {
            // If missing, create the new entry with power_output as fan_power_W / 1000
            let power_output_zero = fan_power_w
                .iter()
                .map(|fpw| fpw / WATTS_PER_KILOWATT as f64)
                .collect::<Vec<f64>>();
            let new_entry = FanSpeedData {
                temperature_diff: 0.0,
                power_output: power_output_zero,
            };
            fan_speed_data.push(new_entry);
        }

        FancoilTestData {
            fan_speed_data,
            fan_power_w,
        }
    }

    pub fn temp_setpnt(&self, simulation_time_iteration: &SimulationTimeIteration) -> Option<f64> {
        match self.heat_source.read().deref() {
            SpaceHeatingService::HeatPump(heat_pump) => {
                heat_pump.temp_setpnt(simulation_time_iteration)
            }
            SpaceHeatingService::Boiler(boiler) => boiler.temp_setpnt(*simulation_time_iteration),
            SpaceHeatingService::HeatNetwork(heat_network) => {
                heat_network.temperature_setpnt(simulation_time_iteration)
            }
            SpaceHeatingService::HeatBattery(heat_battery) => {
                heat_battery.temp_setpnt(*simulation_time_iteration)
            }
        }
    }

    pub(crate) fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool> {
        match self.heat_source.read().deref() {
            SpaceHeatingService::HeatPump(heat_pump) => {
                heat_pump.in_required_period(simulation_time_iteration)
            }
            SpaceHeatingService::Boiler(boiler) => {
                boiler.in_required_period(*simulation_time_iteration)
            }
            SpaceHeatingService::HeatNetwork(heat_network) => {
                heat_network.in_required_period(simulation_time_iteration)
            }
            SpaceHeatingService::HeatBattery(heat_battery) => {
                heat_battery.in_required_period(*simulation_time_iteration)
            }
        }
    }

    pub(crate) fn frac_convective(&self, simtime: SimulationTimeIteration) -> f64 {
        if let Some(WetEmitter::Fancoil {
            frac_convective, ..
        }) = self.fancoil.as_ref().map(|f| f.as_ref())
        {
            return *frac_convective;
        }

        // For other systems
        let frac_convective: Vec<f64> = self
            .emitters
            .iter()
            .map(|emitter| emitter.frac_convective())
            .collect();

        // weighted average for each emitter
        let power_total_weight = self.power_output_emitter_weight(simtime);

        (0..frac_convective.len())
            .map(|i| power_total_weight[i] * frac_convective[i])
            .sum::<f64>()
    }

    /// Weighted average of emitter power output
    fn power_output_emitter_weight(&self, simtime: SimulationTimeIteration) -> Vec<f64> {
        let t_rm = 20.; // assumed internal air temperature

        // flow and return temperatures
        let (flow_temp, return_temp) = self.temp_flow_return(&simtime);
        let t_e = (flow_temp + return_temp) / 2.;

        let power_total_weight = if t_e > t_rm {
            let (power_total, power_total_list) = self.emitters.iter().fold(
                <(f64, Vec<f64>)>::default(),
                |(mut power_total, mut power_total_list), emitter| {
                    let power_emitter = if let (Ok(c), Ok(n)) = (emitter.c(), emitter.n()) {
                        c * (t_e - t_rm).powf(n)
                    } else {
                        panic!("Fancoil emitter, or an emitter with no c and n values, was logically not expected here");
                    };
                    power_total += power_emitter;
                    power_total_list.push(power_emitter);

                    (power_total, power_total_list)
                },
            );

            power_total_list
                .iter()
                .map(|power_emitter| power_emitter / power_total)
                .collect_vec()
        } else {
            vec![1.0 / self.emitters.len() as f64; self.emitters.len()]
        };

        // If this error is triggered, then there is probably an error in the code above
        // NB. debug_* assertions are elided for production builds
        debug_assert_ne!(
            {
                let six_dec_places_divisor = 1e6;
                (power_total_weight.iter().sum::<f64>() * six_dec_places_divisor).round()
                    / six_dec_places_divisor
            },
            1.0,
            "ERROR: Sum of emitter weightings should equal 1.0, not {}",
            power_total_weight.iter().sum::<f64>()
        );

        power_total_weight
    }

    pub(crate) fn temp_flow_return(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> (f64, f64) {
        let flow_temp = match self.ecodesign_controller_class {
            EcoDesignControllerClass::ClassII
            | EcoDesignControllerClass::ClassIII
            | EcoDesignControllerClass::ClassVI
            | EcoDesignControllerClass::ClassVII => {
                // A heater flow temperature control that varies the flow temperature of
                // water leaving the heat dependant upon prevailing outside temperature
                // and selected weather compensation curve.
                //
                // They feature provision for manual adjustment of the weather
                // compensation curves and therby introduce a technical risk that optimal
                // minimised flow temperatures are not always achieved.

                // use weather temperature at the timestep
                let outside_temp = self.external_conditions.air_temp(simulation_time_iteration);

                let min_flow_temp = self.min_flow_temp.unwrap();
                let max_flow_temp = self.max_flow_temp.unwrap();
                let min_outdoor_temp = self.min_outdoor_temp.unwrap();
                let max_outdoor_temp = self.max_outdoor_temp.unwrap();

                // set outdoor and flow temp limits for weather compensation curve
                if outside_temp < min_outdoor_temp {
                    max_flow_temp
                } else if outside_temp > max_outdoor_temp {
                    min_flow_temp
                } else {
                    // Interpolate
                    // Note: A previous version used numpy interpolate, but this
                    //        seemed to be giving incorrect results, so interpolation
                    //        is implemented manually here.
                    min_flow_temp
                        + (outside_temp - max_outdoor_temp)
                            * ((max_flow_temp - min_flow_temp)
                                / (min_outdoor_temp - max_outdoor_temp))
                }
            }
            _ => self.design_flow_temp,
        };

        let return_temp = if flow_temp >= 70.0 {
            60.0
        } else {
            flow_temp * 6.0 / 7.0
        };

        (flow_temp, return_temp)
    }

    /// Calculate emitter output at given emitter and room temp
    ///
    /// Power output from emitter (eqn from 2020 ASHRAE Handbook p644):
    ///            power_output = c * (T_E - T_rm) ^ n
    ///        where:
    ///            T_E is mean emitter temperature
    ///            T_rm is air temperature in the room/zone
    ///            c and n are characteristic of the emitters (e.g. derived from BS EN 442 tests)
    ///
    /// NB. this method will panic if called where there are fancoils amongst the emitters.
    pub fn power_output_emitter(&self, temp_emitter: f64, temp_rm: f64) -> f64 {
        self.emitters.iter().map(|emitter| {
            if let (Ok(c), Ok(n)) = (emitter.c(), emitter.n()) {
                c * (0_f64).max(temp_emitter - temp_rm).powf(n)
            } else {
                panic!("Fancoil emitter, or an emitter with no c and n values, was logically not expected here");
            }
        }).sum::<f64>()
    }

    /// Calculate emitter temperature that gives required power output at given room temp
    ///
    /// Power output from emitter (eqn from 2020 ASHRAE Handbook p644):
    ///            power_output = c * (T_E - T_rm) ^ n
    ///        where:
    ///            T_E is mean emitter temperature
    ///            T_rm is air temperature in the room/zone
    ///            c and n are characteristic of the emitters (e.g. derived from BS EN 442 tests)
    ///        Rearrange to solve for T_E
    ///
    /// Panics if called when one of the emitters is a fancoil.
    pub fn temp_emitter_req(&self, power_emitter_req: f64, temp_rm: f64) -> anyhow::Result<f64> {
        // extract out c and n values from emitters so we don't need to hold a reference to emitters
        // in following function
        let c_n_pairs = self.extract_c_n_pairs();

        let func_temp_emitter_req = |temp_emitter: f64, _args| {
            Ok(power_emitter_req
                - c_n_pairs
                    .iter()
                    .map(|&(c, n)| c * (temp_emitter - temp_rm).powf(n))
                    .sum::<f64>())
        };

        fsolve(func_temp_emitter_req, temp_rm + 10., [])
    }

    /// Extract out c and n values from emitters so we don't need to hold a reference to emitters elsewhere
    /// Panics if called when one of the emitters is a fancoil.
    fn extract_c_n_pairs(&self) -> Vec<(f64, f64)> {
        self.emitters.iter().map(|emitter| {
            if let (Ok(c), Ok(n)) = (emitter.c(), emitter.n()) {
                (c, n)
            } else {
                panic!("Fancoil emitter, or an emitter with no c and n values, was logically not expected here");
            }
        }).collect()
    }

    pub(crate) fn func_temp_emitter_change_rate(
        &self,
        power_input: f64,
    ) -> impl Fn([f64; 1]) -> f64 {
        /*
            Differential eqn for change rate of emitter temperature, to be solved iteratively

            Derivation:

            Heat balance equation for radiators:
                (T_E(t) - T_E(t-1)) * K_E / timestep = power_input - power_output
            where:
                T_E is mean emitter temperature
                K_E is thermal mass of emitters

            Power output from emitter (eqn from 2020 ASHRAE Handbook p644):
                power_output = c * (T_E(t) - T_rm) ^ n
            where:
                T_rm is air temperature in the room/zone
                c and n are characteristic of the emitters (e.g. derived from BS EN 442 tests)

            Substituting power output eqn into heat balance eqn gives:
                (T_E(t) - T_E(t-1)) * K_E / timestep = power_input - c * (T_E(t) - T_rm) ^ n

            Rearranging gives:
                (T_E(t) - T_E(t-1)) / timestep = (power_input - c * (T_E(t) - T_rm) ^ n) / K_E
            which gives the differential equation as timestep goes to zero:
                d(T_E)/dt = (power_input - c * (T_E - T_rm) ^ n) / K_E

            If T_rm is assumed to be constant over the time period, then the rate of
            change of T_E is the same as the rate of change of deltaT, where:
                deltaT = T_E - T_rm

            Therefore, the differential eqn can be expressed in terms of deltaT:
                d(deltaT)/dt = (power_input - c * deltaT(t) ^ n) / K_E

            This can be solved for deltaT over a specified time period using the
            solve_ivp function from scipy.
        */

        let c_n_pairs = self.extract_c_n_pairs();
        let thermal_mass = self.thermal_mass;

        // (power_input - self.c * max_of_2(0., y).powf(self.n)) / self.thermal_mass
        move |temp_diff: [f64; 1]| -> f64 {
            (power_input
                - c_n_pairs
                    .iter()
                    .map(|&(c, n)| c * 0_f64.max(temp_diff[0]).powf(n))
                    .sum::<f64>())
                / thermal_mass
        }
    }

    /// Calculate emitter temperature after specified time with specified power input
    fn temp_emitter(
        &self,
        time_start: f64,
        time_end: f64,
        temp_emitter_start: f64,
        temp_rm: f64,
        power_input: f64,
        temp_emitter_max: Option<f64>,
    ) -> (f64, Option<f64>) {
        // Calculate emitter temp at start of timestep
        let temp_diff_start = temp_emitter_start - temp_rm;
        let temp_diff_max = temp_emitter_max.map(|emitter_max| emitter_max - temp_rm);

        let emitter_with_power_input =
            EmittersAndPowerInput::new(self, power_input, temp_diff_max, temp_diff_start);

        let f = emitter_with_power_input; // f - Structure implementing the System trait
        let x: Time = time_start; // x - Initial value of the independent variable (usually time)
        let x_end: Time = time_end; // x_end - Final value of the independent variable
        let dx = 0.; // dx - Increment in the dense output. This argument has no effect if the output type is Sparse
        let y0: State = State::new(temp_diff_start); // y - Initial value of the dependent variable(s)

        // scipy implementation for reference:
        // https://github.com/scipy/scipy/blob/6b657ede0c3c4cffef3156229afddf02a2b1d99a/scipy/integrate/_ivp/rk.py#L293
        let rtol = 1e-3; // rtol - set from scipy docs - Relative tolerance used in the computation of the adaptive step size
        let atol = 1e-6; // atol - set from scipy docs - Absolute tolerance used in the computation of the adaptive step size
        let h = 0.; // initial step size - 0
        let safety_factor = 0.9; // matches scipy implementation
        let beta = 0.; // setting this to 0 gives us an alpha of 0.2 and matches scipy's adaptive step size logic (default was 0.04)
        let fac_min = 0.2; // matches scipy implementation
        let fac_max = 10.; // matches scipy implementation
        let h_max = x_end - x;
        let n_max = 100000;
        let n_stiff = 1000;
        let mut stepper = Dopri5::from_param(
            f,
            x,
            x_end,
            dx,
            y0,
            rtol,
            atol,
            safety_factor,
            beta,
            fac_min,
            fac_max,
            h_max,
            h,
            n_max,
            n_stiff,
            OutputType::Sparse,
        );

        let _ = stepper.integrate();

        // similar logic to EmittersAndPowerOutput System
        // in future we could consolidate these
        let temp_diff_max_was_reached = match temp_diff_max {
            Some(temp_diff_max) => {
                let y_count = stepper.y_out().len();
                let current_y = stepper.y_out().last().expect("y_out was empty")[0];
                let previous_y = stepper.y_out().get(y_count - 2).unwrap()[0];

                let current_temp_diff = current_y - temp_diff_max;
                let previous_temp_diff = previous_y - temp_diff_max;

                previous_temp_diff == 0.
                    || current_temp_diff == 0.
                    || signs_are_different(previous_temp_diff, current_temp_diff)
            }
            None => false,
        };

        let temp_emitter;
        let mut time_temp_diff_max_reached: Option<f64> = None;
        if temp_diff_max_was_reached {
            // We stopped early because the temp diff max was passed.
            // The Python code uses a built in feature of scipy's solve_ivp here.
            // when an "event" (in this case, max temp diff) happens a root solver
            // finds the exact x (time) value for that event occuring
            // and sets time_temp_diff_max_reached
            // We use a combination of ode_solvers and argmin to achieve the same.

            // max temp diff was reached, so that should be our result
            temp_emitter = temp_rm + temp_diff_max.unwrap();

            let root_problem = RootProblem {
                stepper: &stepper,
                max_temp: temp_diff_max.unwrap(),
            };

            let previous_step_x = *stepper.x_out().get(stepper.x_out().len() - 2).unwrap();
            let current_step_x = *stepper.x_out().last().unwrap();

            let tol = 1e-3; // From scipy docs (rtol default)
                            // Some time (x) between the previous step and the current step we passed the max temp
                            // Use a root solver to find when that was - i.e. when temp - max = 0
            let solver = BrentRoot::new(previous_step_x, current_step_x, tol);

            let executor = Executor::new(root_problem, solver);
            let res = executor.run();

            if res.is_err() {
                panic!("An error occurred in the root solver for emitters")
            }

            let best_x = res.unwrap().state().best_param;
            time_temp_diff_max_reached = best_x;
        } else {
            let last_y = stepper.y_out().last().expect("y_out was empty")[0];
            let temp_diff_emitter_rm_final = last_y;
            temp_emitter = temp_rm + temp_diff_emitter_rm_final;
        }

        (temp_emitter, time_temp_diff_max_reached)
    }

    fn format_fancoil_manufacturer_data(
        fancoil_test_data: &FancoilTestData,
    ) -> anyhow::Result<(Vec<Vec<f64>>, Vec<StringOrNumber>)> {
        let FancoilTestData {
            fan_speed_data,
            fan_power_w,
        } = fancoil_test_data;

        let temperature_data = fan_speed_data
            .iter()
            .map(|entry| {
                let mut row = vec![entry.temperature_diff];
                row.extend(entry.power_output.iter());
                row
            })
            .collect_vec();

        // Check all the fan speed lists are of the same length
        if let Some(first_row_length) = temperature_data.first().map(|row| row.len()) {
            if !temperature_data
                .iter()
                .all(|row| row.len() == first_row_length)
            {
                bail!("Fan speed lists of fancoil manufacturer data differ in length")
            }
        } else {
            bail!("Fan speed data is unexpectedly empty")
        }

        // Prepare fan power data
        let fan_power_data: Vec<StringOrNumber> = {
            let mut row = vec!["Fan power (W)".into()];
            row.extend(fan_power_w.iter().map(|power| StringOrNumber::from(*power)));
            row
        };

        // Check if the length of fan power matches the number of power outputs
        // (subtract 1 because of the temperature_diff column)
        if fan_power_data.len() != temperature_data[0].len() - 1 {
            bail!("Fan power data length does not match the length of fan speed data",)
        }

        Ok((temperature_data, fan_power_data))
    }

    /// Calculate the power output (kW) from fan coil manufacturer data.
    /// For a given delta T, interpolate values for each fan speed column, and
    /// the maximum is calculated. The actual power output is the minimum between this maximum and
    /// the heat demanded.
    ///
    /// Parameters:
    /// * `delta_T_fancoil` (float): temp dif between primary circuit water temp (average of flow and return temp)
    ///                          and room air temp.
    /// * `temperature_data` (array): product data relating to temperature diff from manufacturer.
    /// * `fan_power_data` (array): product data relating to fan power from manufacturer
    /// * `power_req_from_fan_coil` (float): in kW.
    ///
    ///        Returns:
    ///        Tuple containing the actual power output, fan power and fraction of timestep running.
    fn fancoil_output(
        &self,
        delta_t_fancoil: f64,
        temperature_data: &[Vec<f64>],
        fan_power_data: &[StringOrNumber],
        power_req_from_fan_coil: f64,
    ) -> (f64, f64, f64) {
        let delta_t_values = temperature_data.iter().map(|row| row[0]).collect_vec();
        let _min_delta_t = delta_t_values.iter().min_by(|a, b| a.total_cmp(b));

        let fan_power_values = fan_power_data[1..].iter().map(|&f| {
            <f64>::try_from(f).expect("fan power data is not a number when expected to be")
        });

        // Parsing product data to get outputs and fan speeds
        let interpolated_outputs: Vec<_> = (1..temperature_data[1].len())
            .map(|col| {
                let output_values_for_fan_speed = temperature_data.iter().map(|row| row[col]);
                let delta_t_output_pairs = delta_t_values.iter().zip(output_values_for_fan_speed);
                let unique_delta_t_output_pairs =
                    delta_t_output_pairs.unique_by(|&(&a, b)| (OrderedFloat(a), OrderedFloat(b)));
                let sorted_delta_t_output_pairs = unique_delta_t_output_pairs
                    .map(|(&a, b)| (OrderedFloat(a), OrderedFloat(b)))
                    .sorted()
                    .map(|(a, b)| (a.0, b.0));
                let (sorted_delta_t_values, sorted_output_values_for_fan_speed): (
                    Vec<f64>,
                    Vec<f64>,
                ) = sorted_delta_t_output_pairs.unzip();

                // Find the min and max values from the output
                let _min_output_value = sorted_output_values_for_fan_speed
                    .iter()
                    .min_by(|a, b| a.total_cmp(b))
                    .expect("Expected at least one output value");
                let _max_output_value = sorted_output_values_for_fan_speed
                    .iter()
                    .max_by(|a, b| a.total_cmp(b))
                    .expect("Expected at least one output value");

                // TODO (from Python): Currently interpolation follows a linear equation. We think it can be improved with
                // an equation of the form output = c + a * deltaT ^ b that gives a good fit (where c is the fan power)

                // Interpolate value for the given delta T and fan speed output column
                // NB. fill values are specified here in Python, not yet implemented in Rust (TODO?)
                np_interp(
                    delta_t_fancoil,
                    &sorted_delta_t_values,
                    &sorted_output_values_for_fan_speed,
                )
            })
            .collect_vec();

        let fancoil_max_output = *interpolated_outputs
            .iter()
            .max_by(|a, b| a.total_cmp(b))
            .expect("Expected at least one output value");
        let fancoil_min_output = *interpolated_outputs
            .iter()
            .min_by(|a, b| a.total_cmp(b))
            .expect("Expected at least one output value");

        let mut actual_output = power_req_from_fan_coil.min(fancoil_max_output);
        let fraction_timestep_running = if fancoil_min_output == 0. {
            1.
        } else {
            1f64.min(actual_output / fancoil_min_output)
        };

        let fan_power_value = if actual_output <= 0. {
            actual_output = 0.;
            0.
        } else {
            let interpolated_output_pairs = interpolated_outputs.iter().zip(fan_power_values);
            let unique_interpolated_output_pairs =
                interpolated_output_pairs.unique_by(|&(&a, b)| (OrderedFloat(a), OrderedFloat(b)));
            let sorted_interpolated_output_pairs = unique_interpolated_output_pairs
                .map(|(&a, b)| (OrderedFloat(a), OrderedFloat(b)))
                .sorted()
                .map(|(a, b)| (a.0, b.0));
            let (sorted_interpolated_values, sorted_fan_power_values): (Vec<f64>, Vec<f64>) =
                sorted_interpolated_output_pairs.unzip();

            // Find the min and max fan power values
            let _min_fan_power_value = sorted_fan_power_values
                .iter()
                .min_by(|a, b| a.total_cmp(b))
                .expect("Expected at least one fan power value");
            let _max_fan_power_value = sorted_fan_power_values
                .iter()
                .max_by(|a, b| a.total_cmp(b))
                .expect("Expected at least one fan power value");

            // TODO (from Python): Currently interpolation follows a linear equation. We think it can be improved with
            // an equation of form to be determined that gives a better fit

            // Interpolate fan power without extrapolation
            // NB. in the Python there are defined fill values, unimplemented here
            np_interp(
                actual_output,
                &sorted_interpolated_values,
                &sorted_fan_power_values,
            )
        };

        (actual_output, fan_power_value, fraction_timestep_running)
    }

    fn energy_required_from_heat_source(
        &self,
        energy_demand_heating_period: f64,
        time_heating_start: f64,
        timestep: f64,
        temp_rm_prev: f64,
        temp_emitter_heating_start: f64,
        temp_emitter_req: f64,
        temp_emitter_max: f64,
        temp_return: f64,
        simulation_time: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, bool, Option<BufferTankEmittersDataWithResult>)> {
        // When there is some demand, calculate max. emitter temperature
        // achievable and emitter temperature required, and base calculation
        // on the lower of the two.

        // Calculate extra energy required for emitters to reach temp required
        let energy_req_to_warm_emitters = if let Some(_) = self.fancoil.as_ref() {
            0.0
        } else {
            self.thermal_mass * (temp_emitter_req - temp_emitter_heating_start)
        };

        // Calculate energy input required to meet energy demand
        let energy_req_from_heat_source = max_of_2(
            energy_req_to_warm_emitters + energy_demand_heating_period,
            0.0,
        );
        // potential demand from buffer tank

        let energy_req_from_buffer_tank = energy_req_from_heat_source;

        // === Limit energy to account for maximum emitter temperature ===
        // If emitters are already above max. temp for this timestep,
        // then heat source should provide no energy until emitter temp
        // falls to maximum
        // Otherwise:
        let (energy_provided_by_heat_source_max_min, emitters_data_for_buffer_tank_with_result) =
            if self.fancoil.is_some() || self.temp_emitter_prev() <= temp_emitter_max {
                // If emitters are below max. temp for this timestep, then max energy
                // required from heat source will depend on maximum warm-up rate,
                // which depends on the maximum energy output from the heat source
                let emitters_data_for_buffer_tank = self.with_buffer_tank.then(|| {
                    let power_req_from_buffer_tank = if (timestep - time_heating_start) <= 0.0 {
                        // If there is no time remaining in the timestep, then there
                        // is no power requirement (and we need to avoid div-by-zero)
                        0.0
                    } else {
                        energy_req_from_buffer_tank / (timestep - time_heating_start)
                    };

                    BufferTankEmittersData {
                        temp_emitter_req,
                        power_req_from_buffer_tank,
                        design_flow_temp: self.design_flow_temp,
                        target_flow_temp: self.target_flow_temp.load(Ordering::SeqCst),
                        temp_rm_prev,
                        variable_flow: self.variable_flow_data.into(),
                        temp_diff_emit_dsgn: self.temp_diff_emit_dsgn,
                        min_flow_rate: self.min_flow_rate,
                        max_flow_rate: self.max_flow_rate,
                    }
                });

                self.heat_source.read().energy_output_max(
                    temp_emitter_max,
                    temp_return,
                    emitters_data_for_buffer_tank,
                    simulation_time,
                )?
            } else {
                (Default::default(), None)
            };

        let (energy_req_from_heat_source_max, temp_emitter_max_is_final_temp) =
            if self.fancoil.is_some() {
                (energy_req_from_heat_source, true)
            } else {
                // Radiators and/or UFH
                // Calculate time to reach max. emitter temp at max heat source output
                let power_output_max_min = energy_provided_by_heat_source_max_min / timestep;

                let (temp_emitter, time_temp_emitter_max_reached) = self.temp_emitter(
                    0.0,
                    timestep,
                    self.temp_emitter_prev(),
                    temp_rm_prev,
                    power_output_max_min,
                    Some(temp_emitter_max),
                );

                let (time_in_warmup_cooldown_phase, temp_emitter_max_reached) =
                    match time_temp_emitter_max_reached {
                        None => (timestep, false),
                        Some(time) => (time, true),
                    };

                // Before this time, energy output from heat source is maximum
                let energy_req_from_heat_source_before_temp_emitter_max_reached =
                    power_output_max_min * time_in_warmup_cooldown_phase;

                // After this time, energy output is amount needed to maintain
                // emitter temp (based on emitter output at constant emitter temp)
                let energy_req_from_heat_source_after_temp_emitter_max_reached = self
                    .power_output_emitter(temp_emitter, temp_rm_prev)
                    * (timestep - time_in_warmup_cooldown_phase);

                // Total energy input req from heat source is therefore sum of energy
                // output required before and after max emitter temp reached
                let energy_req_from_heat_source_max =
                    energy_req_from_heat_source_before_temp_emitter_max_reached
                        + energy_req_from_heat_source_after_temp_emitter_max_reached;

                let temp_emitter_max_is_final_temp =
                    temp_emitter_max_reached && temp_emitter_req > temp_emitter_max;

                (
                    energy_req_from_heat_source_max,
                    temp_emitter_max_is_final_temp,
                )
            };

        // Total energy input req from heat source is therefore lower of:
        // - energy output required to meet space heating demand
        // - energy output when emitters reach maximum temperature
        Ok((
            min_of_2(energy_req_from_heat_source, energy_req_from_heat_source_max),
            temp_emitter_max_is_final_temp,
            emitters_data_for_buffer_tank_with_result,
        ))
    }

    fn temp_emitter_prev(&self) -> f64 {
        self.temp_emitter_prev.load(Ordering::SeqCst)
    }

    fn energy_surplus_during_cooldown(
        &self,
        time_cooldown: f64,
        [timestep, energy_demand, temp_rm_prev]: [f64; 3],
    ) -> f64 {
        // Calculate emitter temperature after specified time with no heat input
        let temp_emitter_prev = self.temp_emitter_prev();
        let (time_emitter_no_heat_input, _) = self.temp_emitter(
            0.0,
            time_cooldown,
            temp_emitter_prev,
            temp_rm_prev,
            0.0, // no heat from heat source during initial cool-down
            None,
        );
        let energy_released_from_emitters =
            self.thermal_mass * (temp_emitter_prev - time_emitter_no_heat_input);
        let energy_demand_cooldown = energy_demand * time_cooldown / timestep;

        energy_released_from_emitters - energy_demand_cooldown
    }

    ///
    fn calc_emitter_cooldown(
        &self,
        energy_demand: f64,
        temp_emitter_req: f64,
        temp_rm_prev: f64,
        timestep: f64,
    ) -> anyhow::Result<(f64, f64)> {
        let temp_emitter_prev = self.temp_emitter_prev();
        Ok(if temp_emitter_prev < temp_emitter_req {
            (0.0, temp_emitter_prev)
        } else {
            // Calculate time that emitters are cooling down (accounting for
            // undershoot), during which the heat source does not provide any
            // heat, by iterating to find the end time which leads to the heat
            // output matching the energy demand accumulated so far during the
            // timestep
            // TODO (from Python) Is there a more efficient way to do this than iterating?

            // The starting guess below is the end of the timestep rather
            // than the start because at the start of the timestep the
            // function being solved will effectively be 0 minus 0, which
            // is not the result we are seeking (unless no other exists)
            let time_cooldown = root(
                |time_cooldown, args| self.energy_surplus_during_cooldown(time_cooldown, args),
                timestep,
                [timestep, energy_demand, temp_rm_prev],
                Some(1e-8),
            )?;

            // Limit cooldown time to be within timestep
            let time_heating_start = 0f64.max(time_cooldown.min(timestep));
            // Calculate emitter temperature at heating start time
            let (time_emitter_heating_start, _) = self.temp_emitter(
                0.0,
                time_heating_start,
                temp_emitter_prev,
                temp_rm_prev,
                0.0, // No heat from heat source during initial cool-down
                None,
            );

            (time_heating_start, time_emitter_heating_start)
        })
    }

    /// Demand energy from emitters and calculate how much energy can be provided.
    /// This function is called in a loop inside each timestep, where the return temp
    /// is recalculated in each iteration.
    ///
    /// Arguments:
    /// * `energy_demand` - in kWh
    /// * `temp_flow_target` - flow temp in C
    /// * `temp_return_target` - return temp in C
    /// * `update_heat_source_state` - if False, when heat_source.demand_energy is called,
    ///                                the heat source state does not change.
    /// * `update_temp_emitter_prev` - if False, the emitter temperature is not
    ///                                saved for next timestep.
    /// * `blended_temp_flow` - temp when there is bypass recirculated water.
    ///                         If no recirculated water, the it will be equal to the flow temp.
    fn demand_energy_flow_return(
        &self,
        energy_demand: f64,
        temp_flow_target: f64,
        temp_return_target: f64,
        simtime: SimulationTimeIteration,
        update_heat_source_state: Option<bool>,
        update_temp_emitter_prev: Option<bool>,
        blended_temp_flow: Option<f64>,
    ) -> anyhow::Result<(f64, f64)> {
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);
        let update_temp_emitter_prev = update_temp_emitter_prev.unwrap_or(true);

        let timestep = simtime.timestep;
        let temp_rm_prev = self.zone.temp_internal_air();

        // Calculate target flow and return temperature
        let temp_emitter_max = if let Some(blended_temp_flow) = blended_temp_flow {
            (blended_temp_flow + temp_return_target) / 2.
        } else {
            (temp_flow_target + temp_return_target) / 2.
        };

        // Calculate emitter temperature required
        let power_emitter_req = energy_demand / timestep;
        let temp_emitter_req = if self.fancoil.is_some() {
            temp_emitter_max
        } else {
            self.temp_emitter_req(power_emitter_req, temp_rm_prev)?
        };

        self.target_flow_temp
            .store(temp_flow_target, Ordering::SeqCst);

        let mut emitters_data_for_buffer_tank: Option<BufferTankEmittersDataWithResult> = None;

        let mut temp_emitter_output = StringOrNumber::String(KeyString::from("n/a").unwrap());

        let (
            time_heating_start,
            temp_emitter_heating_start,
            energy_req_from_heat_source,
            temp_emitter_max_is_final_temp,
            fan_energy_kwh,
        ) = if energy_demand <= 0. {
            (0.0, self.temp_emitter_prev(), 0.0, false, 0.0)
        } else {
            let (time_heating_start, temp_emitter_heating_start, fan_energy_kwh) =
                if let Some(WetEmitter::Fancoil {
                    n_units,
                    frac_convective,
                    test_data,
                    temperature_data,
                    fan_power_data,
                }) = self.fancoil.as_ref().map(|f| f.as_ref())
                {
                    let n_units = *n_units as f64;
                    let delta_t_fancoil = temp_emitter_max - temp_rm_prev;
                    let power_req_from_fan_coil = energy_demand / n_units / timestep;
                    let (
                        power_delivered_by_fancoil,
                        fan_power_single_unit,
                        fraction_timestep_running,
                    ) = self.fancoil_output(
                        delta_t_fancoil,
                        temperature_data,
                        fan_power_data,
                        power_req_from_fan_coil,
                    );
                    let power_req_from_heat_source =
                        (power_delivered_by_fancoil - fan_power_single_unit) * n_units;
                    let fan_power = fan_power_single_unit * n_units;
                    let energy_demand = (power_req_from_heat_source + fan_power) * timestep;
                    let fan_energy_kwh = fan_power / WATTS_PER_KILOWATT as f64
                        * timestep
                        * fraction_timestep_running;
                    if update_heat_source_state {
                        if let Some(conn) = self.energy_supply_fan_coil_conn.as_ref() {
                            conn.demand_energy(fan_energy_kwh, simtime.index)?;
                        }
                    }

                    // Emitters (fancoils) don't have a warning up or cooling down period:
                    let time_heating_start = 0.0;
                    let temp_emitter_heating_start = self.temp_emitter_prev();

                    (
                        time_heating_start,
                        temp_emitter_heating_start,
                        fan_energy_kwh,
                    )
                } else {
                    let fan_energy_kwh = 0.;
                    // Emitters (radiators and ufh) warming up or cooling down to a target temperature:
                    // - First we calculate the time taken for the emitters to cool
                    //   before the heating system activates, and the temperature that
                    //   the emitters reach at this time. Note that the emitters will
                    //   cool to below the target temperature so that the total heat
                    //   output in this cooling period matches the demand accumulated so
                    //   far in the timestep (assumed to be proportional to the fraction
                    //   of the timestep that has elapsed)
                    let (time_heating_start, temp_emitter_heating_start) = self
                        .calc_emitter_cooldown(
                            energy_demand,
                            temp_emitter_req,
                            temp_rm_prev,
                            timestep,
                        )?;

                    (
                        time_heating_start,
                        temp_emitter_heating_start,
                        fan_energy_kwh,
                    )
                };

            // Then, we calculate the energy required from the heat source in
            // the remaining part of the timestep - (full timestep for fancoils)
            let (
                energy_req_from_heat_source,
                temp_emitter_max_is_final_temp,
                new_emitters_data_for_buffer_tank,
            ) = self.energy_required_from_heat_source(
                (energy_demand - fan_energy_kwh) * (1.0 - time_heating_start / timestep),
                time_heating_start,
                timestep,
                temp_rm_prev,
                temp_emitter_heating_start,
                temp_emitter_req,
                temp_emitter_max,
                temp_return_target,
                simtime,
            )?;

            if let Some(new_emitters_data_for_buffer_tank) = new_emitters_data_for_buffer_tank {
                emitters_data_for_buffer_tank.replace(new_emitters_data_for_buffer_tank);
            }

            (
                time_heating_start,
                temp_emitter_heating_start,
                energy_req_from_heat_source,
                temp_emitter_max_is_final_temp,
                fan_energy_kwh,
            )
        };

        // Get energy output of heat source (i.e. energy input to emitters)
        // TODO (from Python) Instead of passing temp_flow_req into heating system module,
        // calculate average flow temp achieved across timestep?

        // Catering for the possibility of a BufferTank in the emitters' loop
        let (energy_provided_by_heat_source, _) = self.heat_source.write().demand_energy(
            energy_req_from_heat_source,
            temp_flow_target,
            temp_return_target,
            Some(time_heating_start),
            emitters_data_for_buffer_tank.and_then(|data| self.with_buffer_tank.then_some(data)),
            Some(update_heat_source_state),
            simtime,
        )?;

        let energy_released_from_emitters = if self.fancoil.is_some() {
            energy_provided_by_heat_source + fan_energy_kwh
        } else {
            // Calculate emitter temperature achieved at end of timestep.
            // Do not allow emitter temp to rise above maximum
            // Do not allow emitter temp to fall below room temp
            let temp_emitter = if temp_emitter_max_is_final_temp {
                temp_emitter_max
            } else {
                let power_provided_by_heat_source =
                    energy_provided_by_heat_source / (timestep - time_heating_start);
                let (mut temp_emitter, time_temp_target_reached) = self.temp_emitter(
                    time_heating_start,
                    timestep,
                    temp_emitter_heating_start,
                    temp_rm_prev,
                    power_provided_by_heat_source,
                    Some(temp_emitter_req),
                );
                // If target emitter temperature is reached on warm-up, assume that
                // this is maintained to the end of the timestep. This accounts for
                // overshoot and stabilisation without having to model it explicitly
                if temp_emitter_heating_start < temp_emitter_req
                    && time_temp_target_reached.is_some()
                {
                    temp_emitter = temp_emitter_req;
                }

                temp_emitter
            };

            let temp_emitter = temp_emitter.max(temp_rm_prev);

            temp_emitter_output = StringOrNumber::Float(temp_emitter);

            let energy_released_from_emitters = energy_provided_by_heat_source
                + self.thermal_mass * (self.temp_emitter_prev() - temp_emitter);

            // Save emitter temperature for next timestep
            if update_temp_emitter_prev {
                self.temp_emitter_prev.store(temp_emitter, Ordering::SeqCst);
            }

            energy_released_from_emitters
        };

        // If detailed results flag is set populate dict with values
        if self.output_detailed_results && update_heat_source_state {
            let result = EmittersDetailedResult {
                timestep_index: simtime.index,
                energy_demand,
                temp_emitter_req,
                time_heating_start,
                energy_provided_by_heat_source,
                temp_emitter: temp_emitter_output,
                temp_emitter_max,
                energy_released_from_emitters,
                temp_flow_target,
                temp_return_target,
                temp_emitter_max_is_final_temp,
                energy_req_from_heat_source,
                fan_energy_kwh,
            };
            self.emitters_detailed_results.as_ref().unwrap().write()[simtime.index] = result;
        }

        Ok((energy_released_from_emitters, energy_req_from_heat_source))
    }

    /// Energy released from emitters after doing a previous loop
    /// that updates the return temperature.
    pub(crate) fn demand_energy(
        &mut self,
        energy_demand: f64,
        simulation_time: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // ecodesign controls to determine flow temperature,
        // and 6/7th rule to calculate the initial return temperature
        let (temp_flow_target, temp_return_target) = self.temp_flow_return(&simulation_time);

        // TODO use real method in migration to 0.32 when implemented
        // let (temp_return_target, blended_temp_flow, flow_rate_m3s) =
        //     self.return_temp_from_flow_rate(energy_demand, temp_flow_target, temp_return_target);
        let (temp_return_target, blended_temp_flow, _flow_rate_m3s) =
            (temp_return_target, Default::default(), 0.);

        // Last call to demand_energy_flow_return that updates the heat source state and other internal variables
        // before going to the next timestep.
        let (energy_released_from_emitters, _) = self.demand_energy_flow_return(
            energy_demand,
            temp_flow_target,
            temp_return_target,
            simulation_time,
            Some(true),
            Some(true),
            Some(blended_temp_flow),
        )?;

        Ok(energy_released_from_emitters)
    }

    /// Calculate the return temperature for a given flow temp.
    /// If, for a given design delta T, the corresponding flow rate is in
    /// the allowed range, then the return temp is given directly.
    /// If the flow rate is out of the allowed range or the flow rate is fixed,
    /// (no change with timesteps) then the return temp is calculated by iteration.
    /// If there is bypass recirculated water the blended temperature is returned too.
    ///
    /// Arguments:
    /// * `energy_demand` - in kWh
    /// * `temp_flow_target` - flow temp in C
    /// * `temp_return_target` - return temp in C
    fn return_temp_from_flow_rate(
        &self,
        energy_demand: f64,
        temp_flow_target: f64,
        temp_return_target: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64, f64)> {
        let update_heat_source_state = false;
        let update_temp_emitter_prev = false;
        let specific_heat_capacity = WATER.specific_heat_capacity() / JOULES_PER_KILOJOULE as f64;
        let density = WATER.density() * LITRES_PER_CUBIC_METRE as f64;

        // The heat source can modulate the flow rate
        let flow_rate_m3s = match self.variable_flow_data {
            VariableFlowData::Yes => {
                // The return temperature is calculated from temp_diff_emit_dsgn (not the 6/7th rule).
                let temp_return_target = temp_flow_target - self.temp_diff_emit_dsgn;
                let (energy_released_from_emitters, energy_required_from_heat_source) = self
                    .demand_energy_flow_return(
                        energy_demand,
                        temp_flow_target,
                        temp_return_target,
                        simtime,
                        Some(update_heat_source_state),
                        Some(update_temp_emitter_prev),
                        None,
                    )?;

                if energy_released_from_emitters < 0. || energy_required_from_heat_source <= 0. {
                    return Ok((temp_flow_target, temp_flow_target, 0.0));
                }

                // The flow rate is calculated from energy_released_from_emitters and delta T
                let power_released_from_emitters = energy_released_from_emitters / simtime.timestep;
                let flow_rate_m3s = power_released_from_emitters
                    / (specific_heat_capacity * density * self.temp_diff_emit_dsgn);
                let flow_rate = flow_rate_m3s * LITRES_PER_CUBIC_METRE as f64;

                let (flow_rate, flow_rate_in_range) = if flow_rate < self.min_flow_rate {
                    (self.min_flow_rate, false)
                } else if flow_rate > self.max_flow_rate {
                    (self.max_flow_rate, false)
                } else {
                    (flow_rate, true)
                };

                let flow_rate_m3s = flow_rate / LITRES_PER_CUBIC_METRE as f64;

                if flow_rate_in_range {
                    // The heat source can operate at this flow rate, so no need of loop.

                    // If there is bypass recirculated water, blended temp is calculated and return temp reduced accordingly.
                    let blended_temp_flow_target = self.blended_temp(
                        temp_flow_target,
                        temp_return_target,
                        self.bypass_percentage_recirculated,
                    );
                    let temp_return_target =
                        temp_return_target - (blended_temp_flow_target - temp_flow_target).abs();
                    return Ok((temp_return_target, blended_temp_flow_target, flow_rate_m3s));
                }

                flow_rate_m3s
            }
            VariableFlowData::No { design_flow_rate } => {
                let (emitters_released_from_emitters, energy_required_from_heat_source) = self
                    .demand_energy_flow_return(
                        energy_demand,
                        temp_flow_target,
                        temp_return_target,
                        simtime,
                        update_heat_source_state.into(),
                        update_temp_emitter_prev.into(),
                        Default::default(),
                    )?;

                if energy_required_from_heat_source <= 0. {
                    return Ok((temp_flow_target, temp_flow_target, 0.0));
                } else {
                    self.design_flow_temp / LITRES_PER_CUBIC_METRE as f64
                }
            }
        };

        // Loop when the flow rate is constant (design_flow_rate). The initial return temp is the 6/7th rule.
        // Also, for the case of variable flow rate with flow rate out of the allowed range.
        // In this case the initial return temp is calculated from the temp_diff_emit_dsgn.
        let temp_return_target = self.update_return_temp(
            energy_demand,
            temp_flow_target,
            temp_return_target,
            specific_heat_capacity,
            density,
            flow_rate_m3s,
            simtime,
            update_heat_source_state,
            update_temp_emitter_prev,
        )?;

        // If there is bypass recirculated water, blended temp is calculated and return temp reduced accordingly.
        let blended_temp_flow_target = self.blended_temp(
            temp_flow_target,
            temp_return_target,
            self.bypass_percentage_recirculated,
        );
        let mut temp_return_target =
            temp_return_target - (blended_temp_flow_target - temp_flow_target).abs();

        if self.bypass_percentage_recirculated > 0. {
            // Loop again but this time using blended temp and initial reduced return temp.
            temp_return_target = self.update_return_temp(
                energy_demand,
                blended_temp_flow_target,
                temp_return_target,
                specific_heat_capacity,
                density,
                flow_rate_m3s,
                simtime,
                update_heat_source_state,
                update_temp_emitter_prev,
            )?;
        }

        Ok((temp_return_target, blended_temp_flow_target, flow_rate_m3s))
    }

    fn blended_temp(
        &self,
        temp_flow_target: f64,
        temp_return_target: f64,
        bypass_percentage_recirculated: f64,
    ) -> f64 {
        // When there is bypass recirculated water, the blended temperature is calculated following
        // the formula of final temperature of the water mixture T(final)=(m1*T1+m2*T2)/(m1+m2)
        (temp_flow_target + bypass_percentage_recirculated * temp_return_target)
            / (1. + bypass_percentage_recirculated)
    }

    /// Calculate the return temperature for a given flow temperature using fsolve.
    ///
    /// Arguments:
    /// * `energy_demand` - in kWh
    /// * `temp_flow_target` - flow temp in C
    /// * `temp_return_target` - initial guess for the return temperature.
    /// * `specific_heat_capacity` - water specific heat capacity (kJ/KgC)
    /// * `density` - water density (kg/m3)
    /// * `flow_rate_m3s` - flow rate (m3/s)
    /// * `update_heat_source_state` - if False then heat source state not updated.
    /// * `update_temp_emitter_prev` - if False then emitter temperature is not updated for next time step.
    fn update_return_temp(
        &self,
        energy_demand: f64,
        temp_flow_target: f64,
        temp_return_target: f64,
        specific_heat_capacity: f64,
        density: f64,
        flow_rate_m3s: f64,
        simtime: SimulationTimeIteration,
        update_heat_source_state: bool,
        update_temp_emitter_prev: bool,
    ) -> anyhow::Result<f64> {
        let energy_difference = |temp_return: f64, _: [f64; 0]| -> anyhow::Result<f64> {
            let (energy_released_from_emitters, _) = self.demand_energy_flow_return(
                energy_demand,
                temp_flow_target,
                temp_return,
                simtime,
                Some(update_heat_source_state),
                Some(update_temp_emitter_prev),
                None,
            )?;
            let power_released_from_emitters = energy_released_from_emitters / simtime.timestep;
            let calculated_power =
                specific_heat_capacity * density * flow_rate_m3s * (temp_flow_target - temp_return);
            Ok(power_released_from_emitters - calculated_power)
        };

        // Use fsolve to find the return temperature that makes energy_difference zero
        let initial_guess = temp_return_target;

        // TODO original Python call is:
        // fsolve(energy_difference, initial_guess, xtol=1e-2, maxfev=100)[0]
        // and the extra parameters are not reflected here. this may matter
        let temp_return_target = fsolve(energy_difference, initial_guess, [])?;

        Ok(if temp_return_target > temp_flow_target {
            temp_flow_target
        } else {
            temp_return_target
        })
    }

    // Commented out in Python:
    //
    // /// Return the cumulative running time and throughput factor for the heat source
    // /// Arguments:
    // /// energy_demand -- in kWh
    // /// space_heat_running_time_cumulative
    // ///     -- running time spent on higher-priority space heating services
    // pub(crate) fn running_time_throughput_factor(
    //     &self,
    //     energy_demand: f64,
    //     space_heat_running_time_cumulative: f64,
    //     simulation_time: SimulationTimeIteration,
    // ) -> anyhow::Result<(f64, f64)> {
    //     let timestep = simulation_time.timestep;
    //     let temp_rm_prev = self.zone.temp_internal_air();
    //
    //     // Calculate target flow and return temperature
    //     let (temp_flow_target, temp_return_target) = self.temp_flow_return(&simulation_time);
    //     let temp_emitter_max = (temp_flow_target + temp_return_target) / 2.;
    //
    //     let energy_req_from_heat_source = 0.0;
    //
    //     // commented out while migrating to 0.32
    //     // let energy_req_from_heat_source = if energy_demand > 0. {
    //     //     // Emitters warming up or cooling down to a target temperature
    //     //     self.energy_required_from_heat_source(
    //     //         energy_demand,
    //     //         timestep,
    //     //         temp_rm_prev,
    //     //         temp_emitter_max,
    //     //         temp_return_target,
    //     //         simulation_time,
    //     //     )
    //     //     .0
    //     // } else {
    //     //     // Emitters cooling down or at steady-state with heating off
    //     //     0.
    //     // };
    //
    //     self.heat_source.read().running_time_throughput_factor(
    //         space_heat_running_time_cumulative,
    //         energy_req_from_heat_source,
    //         temp_flow_target,
    //         temp_return_target,
    //         simulation_time,
    //     )
    // }

    pub(crate) fn output_emitter_results(&self) -> Option<Vec<EmittersDetailedResult>> {
        self.emitters_detailed_results
            .as_ref()
            .map(|results| results.read().clone())
    }

    /// Calculate minimum possible energy output
    pub(crate) fn energy_output_min(&self, simtime: SimulationTimeIteration) -> f64 {
        if self.fancoil.is_some() {
            0.
        } else {
            let timestep = simtime.timestep;
            let temp_rm_prev = self.zone.temp_internal_air();

            let (temp_emitter, _) = self.temp_emitter(
                0.0,
                timestep,
                self.temp_emitter_prev(),
                temp_rm_prev,
                0.0,
                None,
            );
            let temp_emitter = temp_emitter.max(temp_rm_prev);

            // Calculate emitter output achieved at end of timestep
            self.thermal_mass * (self.temp_emitter_prev() - temp_emitter)
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct EmittersDetailedResult {
    timestep_index: usize,
    energy_demand: f64,
    temp_emitter_req: f64,
    time_heating_start: f64,
    energy_provided_by_heat_source: f64,
    temp_emitter: StringOrNumber,
    temp_emitter_max: f64,
    energy_released_from_emitters: f64,
    temp_flow_target: f64,
    temp_return_target: f64,
    temp_emitter_max_is_final_temp: bool,
    energy_req_from_heat_source: f64,
    fan_energy_kwh: f64,
}

impl EmittersDetailedResult {
    pub(crate) fn as_string_values(&self) -> Vec<String> {
        vec![
            self.timestep_index.to_string(),
            self.energy_demand.to_string(),
            self.temp_emitter_req.to_string(),
            self.time_heating_start.to_string(),
            self.energy_provided_by_heat_source.to_string(),
            self.temp_emitter.to_string(),
            self.temp_emitter_max.to_string(),
            self.energy_released_from_emitters.to_string(),
            self.temp_flow_target.to_string(),
            self.temp_return_target.to_string(),
            self.temp_emitter_max_is_final_temp.to_string(),
            self.energy_req_from_heat_source.to_string(),
            self.fan_energy_kwh.to_string(),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::controls::time_control::Control;
    use crate::core::controls::time_control::OnOffTimeControl;
    use crate::core::energy_supply::energy_supply::EnergySupplyBuilder;
    use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
    use crate::core::heating_systems::boiler::Boiler;
    use crate::core::heating_systems::boiler::BoilerServiceSpace;
    use crate::external_conditions::DaylightSavingsConfig;
    use crate::external_conditions::ShadingSegment;
    use crate::input::FuelType;
    use crate::input::HeatSourceLocation;
    use crate::input::HeatSourceWetDetails;
    use crate::simulation_time::SimulationTime;
    use crate::simulation_time::SimulationTimeIterator;
    use approx::assert_relative_eq;
    use itertools::Itertools;
    use parking_lot::RwLock;
    use rstest::fixture;
    use rstest::rstest;
    use serde_json::json;

    const EIGHT_DECIMAL_PLACES: f64 = 1e-7;
    const FOUR_DECIMAL_PLACES: f64 = 1e-3;

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 8., 1.)
    }

    #[fixture]
    pub(crate) fn simulation_time_iterator() -> SimulationTimeIterator {
        simulation_time().iter()
    }

    #[fixture]
    pub(crate) fn external_conditions(
        simulation_time_iterator: SimulationTimeIterator,
    ) -> ExternalConditions {
        let wind_speeds = vec![3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4];
        let wind_directions = vec![200., 220., 230., 240., 250., 260., 260., 270.];
        let air_temps = vec![0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0];
        let diffuse_horizontal_radiations = vec![333., 610., 572., 420., 0., 10., 90., 275.];
        let direct_beam_radiations = vec![420., 750., 425., 500., 0., 40., 0., 388.];
        let shading_segments = vec![
            ShadingSegment {
                number: 1,
                start: 180.,
                end: 135.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 2,
                start: 135.,
                end: 90.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 3,
                start: 90.,
                end: 45.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 4,
                start: 45.,
                end: 0.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 5,
                start: 0.,
                end: -45.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 6,
                start: -45.,
                end: -90.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 7,
                start: -90.,
                end: -135.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 8,
                start: -135.,
                end: -180.,
                shading_objects: None,
                ..Default::default()
            },
        ];
        ExternalConditions::new(
            &simulation_time_iterator,
            air_temps,
            wind_speeds,
            wind_directions,
            diffuse_horizontal_radiations,
            direct_beam_radiations,
            vec![0.2; 8760],
            51.42,
            -0.75,
            0,
            0,
            None,
            1.0,
            Some(1),
            Some(DaylightSavingsConfig::NotApplicable),
            false,
            false,
            shading_segments,
        )
    }

    #[fixture]
    pub(crate) fn heat_source(
        simulation_time_iterator: SimulationTimeIterator,
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) -> SpaceHeatingService {
        let boiler_details = HeatSourceWetDetails::Boiler {
            energy_supply: "mains gas".to_string(),
            energy_supply_auxiliary: "mains elec".to_string(),
            rated_power: 2.5, // changed to match Python mocks
            efficiency_full_load: 0.891,
            efficiency_part_load: 0.991,
            boiler_location: HeatSourceLocation::Internal,
            modulation_load: 0.3,
            electricity_circ_pump: 0.06,
            electricity_part_load: 0.0131,
            electricity_full_load: 0.0388,
            electricity_standby: 0.0244,
        };
        let energy_supply = Arc::from(RwLock::from(
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time_iterator.total_steps())
                .build(),
        ));

        let energy_supply_conn_aux =
            EnergySupplyConnection::new(energy_supply.clone(), "end_user_name".into());

        let mut boiler = Boiler::new(
            boiler_details,
            energy_supply,
            energy_supply_conn_aux,
            external_conditions.into(),
            simulation_time.step,
        )
        .unwrap();

        let service_name = "service_name";

        boiler
            .create_service_connection(service_name.into())
            .unwrap();

        let control = Arc::from(Control::OnOffTime(OnOffTimeControl::new(
            vec![true, true, true, true, true, true, true, true]
                .iter()
                .map(|&x| Some(x))
                .collect_vec(),
            0,
            0.25, // to match simulation time
        )));

        let boiler_service_space =
            BoilerServiceSpace::new(Arc::new(RwLock::new(boiler)), service_name.into(), control);

        SpaceHeatingService::Boiler(boiler_service_space)
    }

    struct MockZone;

    impl SimpleZone for MockZone {
        fn temp_internal_air(&self) -> f64 {
            20.0
        }

        fn area(&self) -> f64 {
            80.0
        }
    }

    fn emitters_fixture(
        heat_source: SpaceHeatingService,
        external_conditions: ExternalConditions,
        ecodesign_controller: Option<EcoDesignController>,
    ) -> Emitters {
        let thermal_mass = 0.14;
        let c = 0.08;
        let n = 1.2;
        let temp_diff_emit_dsgn = 10.0;
        let zone = Arc::new(MockZone);
        let canned_value = 20.;
        let ecodesign_controller = ecodesign_controller.unwrap_or_else(|| EcoDesignController {
            ecodesign_control_class: EcoDesignControllerClass::ClassII,
            min_outdoor_temp: Some(-4.),
            max_outdoor_temp: Some(20.),
            min_flow_temp: Some(30.),
        });

        let design_flow_temp = 55.;

        let with_buffer_tank = false;

        let wet_emitters = vec![serde_json::from_value(json!(
         {
          "wet_emitter_type": "radiator",
          "c": 0.08,
          "n": 1.2,
          "frac_convective": 0.4
         }
        ))
        .unwrap()];

        Emitters::new(
            Some(thermal_mass),
            &wet_emitters,
            temp_diff_emit_dsgn,
            true,
            None,
            Some(3.),
            Some(18.),
            Some(0.0),
            Arc::new(RwLock::new(heat_source)),
            zone,
            external_conditions.into(),
            ecodesign_controller,
            design_flow_temp,
            None,
            false,
            with_buffer_tank,
        )
        .unwrap()
    }

    #[fixture]
    fn emitters(
        heat_source: SpaceHeatingService,
        external_conditions: ExternalConditions,
    ) -> Emitters {
        emitters_fixture(heat_source, external_conditions, None)
    }

    #[rstest]
    #[ignore = "blocked by temp_emitters issue"]
    fn test_demand_energy(
        simulation_time_iterator: SimulationTimeIterator,
        mut emitters: Emitters,
    ) {
        let energy_demand_list = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0];
        let mut energy_demand = 0.0;

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            energy_demand += energy_demand_list[t_idx];

            let energy_provided = emitters.demand_energy(energy_demand, t_it).unwrap();
            energy_demand -= energy_provided;

            assert_relative_eq!(
                energy_provided,
                [
                    0.26481930394248643,
                    0.8287480680413242,
                    1.053315069769369,
                    1.053315069769369,
                    0.9604801440326911,
                    0.9419772896929609,
                    0.915353814620655,
                    0.7639281136418886
                ][t_idx]
            );

            assert_relative_eq!(
                emitters.temp_emitter_prev(),
                [
                    35.96557640041081,
                    47.20238095238095,
                    47.20238095238095,
                    47.20238095238095,
                    44.78422619047619,
                    44.78422619047619,
                    43.67306169524251,
                    38.21643231208616
                ][t_idx]
            )
        }
    }

    #[fixture]
    fn emitters_with_different_outdoor_temp(
        heat_source: SpaceHeatingService,
        external_conditions: ExternalConditions,
    ) -> Emitters {
        let ecodesign_controller = EcoDesignController {
            ecodesign_control_class: EcoDesignControllerClass::ClassII,
            min_outdoor_temp: Some(10.),
            max_outdoor_temp: Some(15.),
            min_flow_temp: Some(30.),
        };

        emitters_fixture(heat_source, external_conditions, Some(ecodesign_controller))
    }

    #[fixture]
    fn emitters_with_different_control_class(
        heat_source: SpaceHeatingService,
        external_conditions: ExternalConditions,
    ) -> Emitters {
        let ecodesign_controller = EcoDesignController {
            ecodesign_control_class: EcoDesignControllerClass::ClassIV,
            min_outdoor_temp: Some(-4.),
            max_outdoor_temp: Some(20.),
            min_flow_temp: Some(30.),
        };

        emitters_fixture(heat_source, external_conditions, Some(ecodesign_controller))
    }

    /// Test flow and return temperature based on ecodesign control class
    #[rstest]
    fn test_temp_flow_return(
        emitters: Emitters,
        emitters_with_different_outdoor_temp: Emitters,
        emitters_with_different_control_class: Emitters,
        simulation_time_iterator: SimulationTimeIterator,
        heat_source: SpaceHeatingService,
        external_conditions: ExternalConditions,
    ) {
        let (flow_temp, return_temp) =
            emitters.temp_flow_return(&simulation_time_iterator.current_iteration());

        assert_relative_eq!(flow_temp, 50.8333, max_relative = 1e-2);
        assert_relative_eq!(return_temp, 43.5714, max_relative = 1e-2);

        // Test with different outdoor temp
        let (flow_temp, return_temp) = emitters_with_different_outdoor_temp
            .temp_flow_return(&simulation_time_iterator.current_iteration());

        assert_relative_eq!(flow_temp, 55., max_relative = 1e-2);
        assert_relative_eq!(return_temp, 47.1428, max_relative = 1e-2);

        // Test with different control class
        let (flow_temp, return_temp) = emitters_with_different_control_class
            .temp_flow_return(&simulation_time_iterator.current_iteration());

        assert_relative_eq!(flow_temp, 55.);
        assert_relative_eq!(
            return_temp,
            47.14285714,
            max_relative = EIGHT_DECIMAL_PLACES
        );
    }

    // test_temp_flow_return_invalid_value test not copied across from Python
    // as the compile-time checks cover the invalid value case

    /// Test emitter output at given emitter and room temp
    #[rstest]
    #[ignore = "while emitters module being migrated to 0.32"]
    fn test_power_output_emitter(
        heat_source: SpaceHeatingService,
        external_conditions: ExternalConditions,
    ) {
        let ecodesign_controller = EcoDesignController {
            ecodesign_control_class: EcoDesignControllerClass::ClassII,
            min_outdoor_temp: Some(-4.),
            max_outdoor_temp: Some(20.),
            min_flow_temp: Some(30.),
        };

        let emitters =
            emitters_fixture(heat_source, external_conditions, Some(ecodesign_controller));

        let temp_emitter = 15.;
        let temp_rm = 10.;

        let result = emitters.power_output_emitter(temp_emitter, temp_rm);

        assert_relative_eq!(result, 0.55189186, max_relative = EIGHT_DECIMAL_PLACES);
    }

    /// Test emitter temperature that gives required power output at given room temp
    #[rstest]
    #[ignore = "while emitters module being migrated to 0.32"]
    fn test_temp_emitter_req(emitters: Emitters) {
        let power_emitter_req = 0.22;
        let temp_rm = 2.0;

        let result = emitters
            .temp_emitter_req(power_emitter_req, temp_rm)
            .unwrap();

        assert_relative_eq!(result, 4.32332827, max_relative = EIGHT_DECIMAL_PLACES);
    }

    /// Test Differential eqn is formed for change rate of emitter temperature
    #[rstest]
    #[ignore = "while emitters module being migrated to 0.32"]
    fn test_func_temp_emitter_change_rate(emitters: Emitters) {
        let result = emitters.func_temp_emitter_change_rate(5.)([32.]);

        assert_eq!(result, -0.8571428571428514);
    }

    #[rstest]
    #[ignore = "while emitters module being migrated to 0.32"]
    fn test_temp_emitter_with_no_max(emitters: Emitters) {
        // Test function calculates emitter temperature after specified time with specified power input
        // Check None conditions  are invoked
        let (temp_emitter, time_temp_diff_max_reached) =
            emitters.temp_emitter(0., 2., 5., 10., 0.2, None);

        assert_relative_eq!(
            temp_emitter,
            7.85714285,
            max_relative = EIGHT_DECIMAL_PLACES
        );
        assert!(time_temp_diff_max_reached.is_none());
    }

    #[rstest]
    #[case(0., 2., 70., 10., 0.2, 25., 25., 1.29981138)]
    #[case(0., 5., 25., 3., 0.8, 19., 19., 0.44239778)]
    #[case(5., 25., 6., 14., 0.95, 21., 21., 8.42980041)]
    #[case(0., 1., 45.8, 13., 0., 45.2, 45.2, 0.016104431558688848)]
    #[case(0., 1., 40., 30., 0., 40., 40., 0.)]
    #[ignore = "while emitters module being migrated to 0.32"]
    fn test_temp_emitter_with_max(
        emitters: Emitters,
        #[case] time_start: f64,
        #[case] time_end: f64,
        #[case] temp_emitter_start: f64,
        #[case] temp_rm: f64,
        #[case] power_input: f64,
        #[case] temp_emitter_max: f64,
        #[case] expected_temp: f64,
        #[case] expected_time: f64,
    ) {
        // Check not None conditions are invoked
        // Test when max temp is reached (early exit)
        let (temp_emitter, time_temp_diff_max_reached) = emitters.temp_emitter(
            time_start,
            time_end,
            temp_emitter_start,
            temp_rm,
            power_input,
            Some(temp_emitter_max),
        );

        assert_relative_eq!(
            temp_emitter,
            expected_temp,
            max_relative = EIGHT_DECIMAL_PLACES
        );
        assert_relative_eq!(
            time_temp_diff_max_reached.unwrap(),
            expected_time,
            max_relative = FOUR_DECIMAL_PLACES
        );
    }

    #[rstest]
    #[case(0., 2., 5., 10., 0., 21., 5.)]
    #[ignore = "while emitters module being migrated to 0.32"]
    fn test_temp_emitter_with_max_not_reached(
        emitters: Emitters,
        #[case] time_start: f64,
        #[case] time_end: f64,
        #[case] temp_emitter_start: f64,
        #[case] temp_rm: f64,
        #[case] power_input: f64,
        #[case] temp_emitter_max: f64,
        #[case] expected_temp: f64,
    ) {
        let (temp_emitter, time_temp_diff_max_reached) = emitters.temp_emitter(
            time_start,
            time_end,
            temp_emitter_start,
            temp_rm,
            power_input,
            Some(temp_emitter_max),
        );

        assert_relative_eq!(
            temp_emitter,
            expected_temp,
            max_relative = EIGHT_DECIMAL_PLACES
        );
        assert!(time_temp_diff_max_reached.is_none());
    }

    #[rstest]
    #[ignore = "blocked by temp_emitters issue"]
    fn test_energy_required_from_heat_source(
        simulation_time_iterator: SimulationTimeIterator,
        emitters: Emitters,
    ) {
        // commented out while migrating to 0.32
        // let energy_demand_list = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0];
        // let mut energy_demand = 0.0;
        //
        // for (t_idx, t_it) in simulation_time_iterator.enumerate() {
        //     energy_demand += energy_demand_list[t_idx];
        //
        //     let (energy_req, temp_emitter_max_is_final_temp, _) =
        //         emitters.energy_required_from_heat_source(energy_demand, 1., 10., 25., 30., t_it);
        //
        //     assert_relative_eq!(
        //         energy_req,
        //         [
        //             0.7487346289045738,
        //             2.458950517761754,
        //             2.458950517761754,
        //             2.458950517761754,
        //             2.458950517761754,
        //             2.458950517761754,
        //             2.458950517761754,
        //             2.458950517761754
        //         ][t_idx]
        //     );
        //
        //     assert_eq!(
        //         temp_emitter_max_is_final_temp,
        //         [false, false, true, true, true, true, true, true][t_idx]
        //     );
        // }
    }

    // Python has a test_running_time_throughput_factor which we've not ported for now as it relies
    // on mocking
}
