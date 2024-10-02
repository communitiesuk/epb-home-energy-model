use crate::compare_floats::{max_of_2, min_of_2};
use crate::core::heating_systems::common::SpaceHeatingService;
use crate::core::heating_systems::heat_pump::BufferTankEmittersData;
use crate::corpus::TempInternalAirFn;
use crate::external_conditions::ExternalConditions;
use crate::input::{EcoDesignController, EcoDesignControllerClass};
use crate::simulation_time::SimulationTimeIteration;
use ode_solvers::{dop_shared::OutputType, Dopri5, System, Vector1};
use std::sync::Arc;

type State = Vector1<f64>; // type State = OVector<f32, U3>;
type Time = f64;
// type Result = SolverResult<Time, State>;

/// This module provides objects to represent radiator and underfloor emitter systems.

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

pub struct Emitters {
    pub thermal_mass: f64,
    pub c: f64,
    pub n: f64,
    temp_diff_emit_dsgn: f64,
    frac_convective: f64,
    heat_source: Arc<SpaceHeatingService>,
    temp_internal_air_fn: TempInternalAirFn,
    external_conditions: Arc<ExternalConditions>,
    with_buffer_tank: bool,
    design_flow_temp: f64,
    ecodesign_controller_class: EcoDesignControllerClass,
    min_outdoor_temp: Option<f64>,
    max_outdoor_temp: Option<f64>,
    min_flow_temp: Option<f64>,
    max_flow_temp: Option<f64>,
    simulation_timestep: f64,
    temp_emitter_prev: f64,
    target_flow_temp: Option<f64>, // In Python this is set from inside demand energy and does not exist before then
}

#[derive(Copy, Clone)]
struct EmittersAndPowerInput<'a> {
    pub emitters: &'a Emitters,
    pub power_input: f64,
    pub temp_diff_max: Option<f64>,
}

impl<'a> System<Time, State> for EmittersAndPowerInput<'a> {
    fn system(&self, _x: Time, y: &State, dy: &mut State) {
        dy[0] = self.emitters.func_temp_emitter_change_rate(self.power_input, y[0]);
    }

    // Stop function called at every successful integration step. The integration is stopped when this function returns true.
    fn solout(&mut self, _x: Time, y: &State, _dy: &State) -> bool {
        if self.temp_diff_max.is_none() {
            // no maximum - keep going
            return false
        }

        // we should stop if we cross the max temp difference
        y[0] < self.temp_diff_max.unwrap()
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
        thermal_mass: f64,
        c: f64,
        n: f64,
        temp_diff_emit_dsgn: f64,
        frac_convective: f64,
        heat_source: Arc<SpaceHeatingService>,
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions: Arc<ExternalConditions>,
        ecodesign_controller: EcoDesignController,
        design_flow_temp: f64,
        simulation_timestep: f64,
        with_buffer_tank: bool,
    ) -> Self {
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
        Self {
            thermal_mass,
            c,
            n,
            temp_diff_emit_dsgn,
            frac_convective,
            heat_source,
            temp_internal_air_fn,
            external_conditions,
            with_buffer_tank,
            design_flow_temp,
            ecodesign_controller_class,
            min_outdoor_temp,
            max_outdoor_temp,
            min_flow_temp,
            max_flow_temp,
            simulation_timestep,
            temp_emitter_prev: 20.0,
            target_flow_temp: None,
        }
    }

    pub fn temp_setpnt(&self, simulation_time_iteration: &SimulationTimeIteration) -> Option<f64> {
        match self.heat_source.as_ref() {
            SpaceHeatingService::HeatPump(heat_pump) => {
                heat_pump.temp_setpnt(simulation_time_iteration)
            }
            SpaceHeatingService::Boiler(boiler) => boiler.temp_setpnt(*simulation_time_iteration),
            SpaceHeatingService::HeatNetwork(heat_network) => {
                heat_network.temperature_setpnt(simulation_time_iteration)
            }
            SpaceHeatingService::HeatBattery(_) => unreachable!(),
        }
    }

    pub fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool> {
        match self.heat_source.as_ref() {
            SpaceHeatingService::HeatPump(heat_pump) => {
                heat_pump.in_required_period(simulation_time_iteration)
            }
            SpaceHeatingService::Boiler(boiler) => {
                boiler.in_required_period(*simulation_time_iteration)
            }
            SpaceHeatingService::HeatNetwork(heat_network) => {
                heat_network.in_required_period(simulation_time_iteration)
            }
            SpaceHeatingService::HeatBattery(_) => unreachable!(),
        }
    }

    pub fn frac_convective(&self) -> f64 {
        self.frac_convective
    }

    pub fn temp_flow_return(
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
    pub fn power_output_emitter(&self, temp_emitter: f64, temp_rm: f64) -> f64 {
        self.c * max_of_2(0., temp_emitter - temp_rm).powf(self.n)
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
    pub fn temp_emitter_req(&self, power_emitter_req: f64, temp_rm: f64) -> f64 {
        (power_emitter_req / self.c).powf(1. / self.n) + temp_rm
    }

    pub fn func_temp_emitter_change_rate(
        &self,
        power_input: f64,
        y: f64,
    ) -> f64 {
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

        (power_input - self.c * max_of_2(0., y).powf(self.n)) / self.thermal_mass
    }

    // /// Calculate emitter temperature after specified time with specified power input
    // pub fn temp_emitter(
    //     &self,
    //     time_start: f64,
    //     time_end: f64,
    //     temp_emitter_start: f64,
    //     temp_rm: f64,
    //     power_input: f64,
    //     temp_emitter_max: Option<f64>,
    // ) -> Result<(f64, f64), &'static str> {
    //     // Calculate emitter temp at start of timestep

    //     let temp_diff_start = temp_emitter_start - temp_rm;
    // }

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
        let temp_diff_max = match temp_emitter_max {
            Some(emitter_max) => Some(emitter_max - temp_rm),
            None => None
        };

        let emitter_with_power_input = EmittersAndPowerInput {
            emitters: self,
            power_input,
            temp_diff_max
        };

        let f = emitter_with_power_input.clone(); // f - Structure implementing the System trait
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
        let last_x = *stepper.x_out().last().expect("x_out was empty");

        let mut temp_emitter = 0.;
        let mut time_temp_diff_max_reached: Option<f64> = None;

        let max_temp_diff_was_reached = last_x != x_end;
        if temp_diff_max.is_some() && max_temp_diff_was_reached {
            // We stopped early because the max diff was passed.
            // The Python code uses a built in feature of scipy's solve_ivp here.
            // when an "event" (in this case, max temp diff) happens a root solver
            // finds the exact x (time) value for that event occuring
            // and sets time_temp_diff_max_reached

            // For now we set this to the last time value
            // but Python will have a more exact result
            time_temp_diff_max_reached = Some(last_x);

            // max temp diff was reached, so that should be our result
            temp_emitter = temp_rm + temp_diff_max.unwrap();
        }
        else {
            let temp_diff_emitter_rm_final = stepper.y_out().last().expect("y_out was empty")[0];
            temp_emitter = temp_rm + temp_diff_emitter_rm_final;   
        }

        (temp_emitter, time_temp_diff_max_reached)
    }

    fn energy_required_from_heat_source(
        &self,
        energy_demand: f64,
        timestep: f64,
        temp_rm_prev: f64,
        temp_emitter_max: f64,
        temp_return: f64,
        simulation_time: SimulationTimeIteration,
    ) -> (f64, bool, Option<BufferTankEmittersData>) {
        // When there is some demand, calculate max. emitter temperature
        // achievable and emitter temperature required, and base calculation
        // on the lower of the two.
        // TODO The final calculation of emitter temperature below assumes
        // constant power output from heating source over the timestep.
        // It does not account for:
        // - overshoot/undershoot and then stabilisation of emitter temp.
        // This leads to emitter temp at end of timestep not exactly
        // matching temp_emitter_target.
        // - On warm-up, calculate max. temp achievable, then cap at target
        // and record time this is reached, then assume target temp is
        // maintained? Would there be an overshoot to make up for
        // underheating during warm-up?
        // - On cool-down to lower target temp, calculate lowest temp achieveable,
        // with target temp as floor and record time this is reached, then
        // assume target temp is maintained? Would there be an undershoot
        // to make up for overheating during cool-down?
        // - other services being served by heat source (e.g. DHW) as a
        // higher priority. This may cause intermittent operation or
        // delayed start which could cause the emitters to cool
        // through part of the timestep. We could assume that all the
        // time spent on higher priority services is at the start of
        // the timestep and run solve_ivp for time periods 0 to
        // higher service running time (with no heat input) and higher
        // service running time to timestep end (with heat input). However,
        // this may not be realistic where there are multiple space
        // heating services which in reality would be running at the same
        // time.

        // Calculate emitter temperature required
        let power_emitter_req = energy_demand / timestep;
        let temp_emitter_req = self.temp_emitter_req(power_emitter_req, temp_rm_prev);

        // Calculate extra energy required for emitters to reach temp required
        let energy_req_to_warm_emitters =
            self.thermal_mass * (temp_emitter_req - self.temp_emitter_prev);
        // Calculate energy input required to meet energy demand
        let energy_req_from_heat_source =
            max_of_2(energy_req_to_warm_emitters + energy_demand, 0.0);
        // potential demand from buffer tank

        let energy_req_from_buffer_tank = energy_req_from_heat_source;

        // === Limit energy to account for maximum emitter temperature ===
        let mut emitters_data_for_buffer_tank: Option<BufferTankEmittersData> = None;
        let mut energy_provided_by_heat_source_max_min: f64 = Default::default();

        // If emitters are already above max. temp for this timestep,
        // then heat source should provide no energy until emitter temp
        // falls to maximum
        // Otherwise:
        if self.temp_emitter_prev <= temp_emitter_max {
            // If emitters are below max. temp for this timestep, then max energy
            // required from heat source will depend on maximum warm-up rate,
            // which depends on the maximum energy output from the heat source

            emitters_data_for_buffer_tank = match self.with_buffer_tank {
                true => Some(BufferTankEmittersData {
                    temp_emitter_req,
                    power_req_from_buffer_tank: energy_req_from_buffer_tank / timestep,
                    design_flow_temp: self.design_flow_temp,
                    target_flow_temp: self
                        .target_flow_temp
                        .expect("Expect a target_flow_temp to have been set at this point"),
                    temp_rm_prev,
                }),
                false => None,
            };

            let (
                energy_provided_by_heat_source_max_min_temp,
                emitters_data_for_buffer_tank_with_result,
            ) = self
                .heat_source
                .energy_output_max(
                    temp_emitter_max,
                    temp_return,
                    emitters_data_for_buffer_tank,
                    simulation_time,
                )
                .unwrap();

            energy_provided_by_heat_source_max_min = energy_provided_by_heat_source_max_min_temp;

            emitters_data_for_buffer_tank = match emitters_data_for_buffer_tank_with_result {
                Some(x) => Some(x.data),
                None => None,
            };
        }

        // Calculate time to reach max. emitter temp at max heat source output
        let power_output_max_min = energy_provided_by_heat_source_max_min / timestep;

        let (temp_emitter, time_temp_emitter_max_reached) = self.temp_emitter(
            0.0,
            timestep,
            self.temp_emitter_prev,
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

        // Total energy input req from heat source is therefore lower of:
        // - energy output required to meet space heating demand
        // - energy output when emitters reach maximum temperature
        (
            min_of_2(energy_req_from_heat_source, energy_req_from_heat_source_max),
            temp_emitter_max_is_final_temp,
            emitters_data_for_buffer_tank,
        )
    }

    /// Demand energy from emitters and calculate how much energy can be provided
    /// Arguments:
    /// energy_demand -- in kWh
    fn demand_energy(
        &mut self,
        energy_demand: f64,
        simulation_time: SimulationTimeIteration,
    ) -> f64 {
        let timestep = simulation_time.timestep;
        let temp_rm_prev = &self.temp_internal_air_fn;

        // Calculate target flow and return temperature
        let (temp_flow_target, temp_return_target) = self.temp_flow_return(&simulation_time);
        let temp_emitter_max = (temp_flow_target + temp_return_target) / 2.0;
        self.target_flow_temp = Some(temp_flow_target);

        let mut energy_req_from_heat_source = Default::default();
        let mut temp_emitter_max_is_final_temp = false;
        let mut emitters_data_for_buffer_tank = None;
        if energy_demand <= 0. {
            // Emitters cooling down or at steady-state with heating off
            // energy_req_from_heat_source = 0.0
            // temp_emitter_max_is_final_temp = False
        } else {
            // Emitters warming up or cooling down to a target temperature
            (
                energy_req_from_heat_source,
                temp_emitter_max_is_final_temp,
                emitters_data_for_buffer_tank,
            ) = self.energy_required_from_heat_source(
                energy_demand,
                timestep,
                temp_rm_prev(),
                temp_emitter_max,
                temp_return_target,
                simulation_time,
            )
        }
        todo!("finish this method")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use parking_lot::RwLock;
    use rstest::fixture;
    use rstest::rstest;

    use crate::core::controls::time_control::Control;
    use crate::core::controls::time_control::OnOffTimeControl;
    use crate::core::energy_supply;
    use crate::core::energy_supply::energy_supply::EnergySupply;
    use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
    use crate::core::heating_systems::boiler::Boiler;
    use crate::core::heating_systems::boiler::BoilerServiceSpace;
    use crate::core::heating_systems::boiler::BoilerServiceWaterRegular;
    use crate::core::heating_systems::common::HeatSourceWet;
    use crate::core::space_heat_demand::zone::Zone;
    use crate::corpus::HeatSource;
    use crate::external_conditions::DaylightSavingsConfig;
    use crate::external_conditions::ShadingSegment;
    use crate::input::EnergySupplyType;
    use crate::input::FuelType;
    use crate::input::HeatSourceLocation;
    use crate::input::HeatSourceWetDetails;
    use crate::simulation_time::SimulationTime;
    use crate::simulation_time::SimulationTimeIterator;

    const EIGHT_DECIMAL_PLACES: f64 = 1e-7;

    #[fixture]
    pub(crate) fn simulation_time() -> SimulationTimeIterator {
        SimulationTime::new(0., 2., 0.25).iter()
    }

    #[fixture]
    pub(crate) fn external_conditions(
        simulation_time: SimulationTimeIterator,
    ) -> ExternalConditions {
        let simulation_time_iterator = simulation_time;
        let wind_speeds = vec![3.7, 3.8];
        let wind_directions = vec![200., 220.];
        let air_temps = vec![0.0, 2.5];
        let diffuse_horizontal_radiations = vec![333., 610.];
        let direct_beam_radiations = vec![420., 750.];
        let shading_segments = vec![
            ShadingSegment {
                number: 1,
                start: 180.,
                end: 135.,
                objects: None,
            },
            ShadingSegment {
                number: 2,
                start: 135.,
                end: 90.,
                objects: None,
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
            DaylightSavingsConfig::NotApplicable,
            false,
            false,
            shading_segments,
        )
    }

    #[fixture]
    pub(crate) fn heat_source(
        simulation_time: SimulationTimeIterator,
        external_conditions: ExternalConditions,
    ) -> SpaceHeatingService {
        let boiler_details = HeatSourceWetDetails::Boiler {
            energy_supply: EnergySupplyType::MainsGas,
            energy_supply_auxiliary: EnergySupplyType::Electricity,
            rated_power: 24.,
            efficiency_full_load: 0.891,
            efficiency_part_load: 0.991,
            boiler_location: HeatSourceLocation::Internal,
            modulation_load: 0.3,
            electricity_circ_pump: 0.06,
            electricity_part_load: 0.0131,
            electricity_full_load: 0.0388,
            electricity_standby: 0.0244,
        };
        let energy_supply = Arc::from(RwLock::from(EnergySupply::new(
            FuelType::MainsGas,
            simulation_time.total_steps(),
            None,
            None,
            None,
        )));

        let energy_supply_conn_aux =
            EnergySupplyConnection::new(energy_supply.clone(), "end_user_name".into());

        let boiler = Boiler::new(
            boiler_details,
            energy_supply,
            energy_supply_conn_aux,
            external_conditions.into(),
            1., // TODO is this correct?
        )
        .unwrap();

        let control = Arc::from(Control::OnOffTimeControl(OnOffTimeControl::new(
            vec![true, true, true, true, true, true, true, true],
            0,
            0.,
        )));

        let boiler_service_space = BoilerServiceSpace::new(boiler, "service_name".into(), control);

        SpaceHeatingService::Boiler(boiler_service_space)
    }

    #[fixture]
    pub(crate) fn emitters(
        heat_source: SpaceHeatingService,
        external_conditions: ExternalConditions,
    ) -> Emitters {
        let thermal_mass = 0.14;
        let c = 0.08;
        let n = 1.2;
        let temp_diff_emit_dsgn = 10.0;
        let frac_convective = 0.4;
        let canned_value = 20.;
        let ecodesign_controller = EcoDesignController {
            ecodesign_control_class: EcoDesignControllerClass::ClassII,
            min_outdoor_temp: Some(-4.),
            max_outdoor_temp: Some(20.),
            min_flow_temp: Some(30.),
        };

        let design_flow_temp = 55.;

        // TODO check this is correct
        let simulation_timestep = 1.;
        // TODO check this is correct
        let with_buffer_tank = false;

        Emitters::new(
            thermal_mass,
            c,
            n,
            temp_diff_emit_dsgn,
            frac_convective,
            heat_source.into(),
            Arc::new(move || canned_value),
            external_conditions.into(),
            ecodesign_controller,
            design_flow_temp,
            simulation_timestep,
            with_buffer_tank,
        )
    }

    #[rstest]
    #[ignore = "not yet implemented"]
    fn test_demand_energy(simulation_time: SimulationTimeIterator, mut emitters: Emitters) {
        let energy_demand_list = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0];
        let mut energy_demand = 0.0;

        for (t_idx, t_it) in simulation_time.enumerate() {
            energy_demand += energy_demand_list[t_idx];

            let energy_provided = emitters.demand_energy(energy_demand, t_it);
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
                emitters.temp_emitter_prev,
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

    /// Test flow and return temperature based on ecodesign control class
    #[rstest]
    fn test_temp_flow_return(
        emitters: Emitters,
        simulation_time: SimulationTimeIterator,
        heat_source: SpaceHeatingService,
        external_conditions: ExternalConditions,
    ) {
        let (flow_temp, return_temp) =
            emitters.temp_flow_return(&simulation_time.current_iteration());

        assert_relative_eq!(flow_temp, 50.8333, max_relative = 1e-2);
        assert_relative_eq!(return_temp, 43.5714, max_relative = 1e-2);

        // Test with different outdoor temp
        let ecodesign_controller = EcoDesignController {
            ecodesign_control_class: EcoDesignControllerClass::ClassII,
            min_outdoor_temp: Some(10.),
            max_outdoor_temp: Some(15.),
            min_flow_temp: Some(30.),
        };

        let heat_source = Arc::new(heat_source);

        let emitters = Emitters::new(
            0.14,
            0.08,
            1.2,
            10.,
            0.4,
            heat_source.clone(),
            Arc::new(move || 20.),
            external_conditions.clone().into(),
            ecodesign_controller,
            55.,
            1.,
            false,
        );

        let (flow_temp, return_temp) =
            emitters.temp_flow_return(&simulation_time.current_iteration());

        assert_relative_eq!(flow_temp, 55., max_relative = 1e-2);
        assert_relative_eq!(return_temp, 47.1428, max_relative = 1e-2);

        // Test with different control class
        let ecodesign_controller = EcoDesignController {
            ecodesign_control_class: EcoDesignControllerClass::ClassIV,
            min_outdoor_temp: Some(-4.),
            max_outdoor_temp: Some(20.),
            min_flow_temp: Some(30.),
        };

        let emitters = Emitters::new(
            0.14,
            0.08,
            1.2,
            10.,
            0.4,
            heat_source,
            Arc::new(move || 20.),
            external_conditions.into(),
            ecodesign_controller,
            55.,
            1.,
            false,
        );

        let (flow_temp, return_temp) =
            emitters.temp_flow_return(&simulation_time.current_iteration());

        assert_relative_eq!(flow_temp, 55.);
        assert_relative_eq!(
            return_temp,
            47.14285714,
            max_relative = EIGHT_DECIMAL_PLACES
        );
    }

    // TODO more tests to implement here
    #[rstest]
    // #[ignore = "not yet implemented"]
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
    fn test_temp_emitter_with_max(emitters: Emitters) {
        // Check not None conditions are invoked
        // Test when max temp is reached (early exit)
        let (temp_emitter, _time_temp_diff_max_reached) =
            emitters.temp_emitter(0., 2., 70., 10., 0.2, Some(25.));

        assert_relative_eq!(temp_emitter, 25., max_relative = EIGHT_DECIMAL_PLACES);

        // This assertion has been split into a separate test below
        //assert_relative_eq!(_time_temp_diff_max_reached.unwrap(), 1.29981138);
    }

    #[rstest]
    #[ignore = "we don't currently match the time_temp_diff_max_reached for emitters"]
    fn test_temp_emitter_with_max_reached_time(emitters: Emitters) {
        // this is split out from the test above
        // because of a known issue matching the time_temp_diff_max_reached

        let (_temp_emitter, time_temp_diff_max_reached) =
            emitters.temp_emitter(0., 2., 70., 10., 0.2, Some(25.));

        assert_relative_eq!(time_temp_diff_max_reached.unwrap(), 1.29981138);
    }
}
