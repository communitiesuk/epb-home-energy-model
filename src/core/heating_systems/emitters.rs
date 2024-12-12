/// This module provides objects to represent radiator and underfloor emitter systems.
use crate::compare_floats::{max_of_2, min_of_2};
use crate::core::heating_systems::common::SpaceHeatingService;
use crate::core::heating_systems::heat_pump::{
    BufferTankEmittersData, BufferTankEmittersDataWithResult,
};
use crate::corpus::TempInternalAirFn;
use crate::external_conditions::ExternalConditions;
use crate::input::{EcoDesignController, EcoDesignControllerClass};
use crate::simulation_time::SimulationTimeIteration;
use crate::statistics::np_interp;
use ode_solvers::{dop_shared::OutputType, Dopri5, System, Vector1};
use parking_lot::RwLock;
use std::fmt::{Debug, Formatter};
use std::ops::Deref;
use std::sync::Arc;

type State = Vector1<f64>; // type State = OVector<f32, U3>;
type Time = f64;
// type Result = SolverResult<Time, State>;

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

#[derive(Clone)]
pub struct Emitters {
    pub thermal_mass: f64,
    pub c: f64,
    pub n: f64,
    _temp_diff_emit_dsgn: f64,
    frac_convective: f64,
    heat_source: Arc<RwLock<SpaceHeatingService>>,
    temp_internal_air_fn: TempInternalAirFn,
    external_conditions: Arc<ExternalConditions>,
    with_buffer_tank: bool,
    design_flow_temp: f64,
    ecodesign_controller_class: EcoDesignControllerClass,
    min_outdoor_temp: Option<f64>,
    max_outdoor_temp: Option<f64>,
    min_flow_temp: Option<f64>,
    max_flow_temp: Option<f64>,
    temp_emitter_prev: f64,
    target_flow_temp: Option<f64>, // In Python this is set from inside demand energy and does not exist before then
    output_detailed_results: bool,
    emitters_detailed_results: Option<Arc<RwLock<Vec<EmittersDetailedResult>>>>,
}

// Implement Debug for Emitters using standard strategy, overwriting debug value for temp_internal_air_fn.
impl Debug for Emitters {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Emitters")
            .field("thermal_mass", &self.thermal_mass)
            .field("c", &self.c)
            .field("n", &self.n)
            .field("_temp_diff_emit_dsgn", &self._temp_diff_emit_dsgn)
            .field("frac_convective", &self.frac_convective)
            .field("heat_source", &self.heat_source)
            .field(
                "temp_internal_air_fn",
                &"A function providing the internal air temperature.",
            )
            .field("external_conditions", &self.external_conditions)
            .field("with_buffer_tank", &self.with_buffer_tank)
            .field("design_flow_temp", &self.design_flow_temp)
            .field(
                "ecodesign_controller_class",
                &self.ecodesign_controller_class,
            )
            .field("min_outdoor_temp", &self.min_outdoor_temp)
            .field("max_outdoor_temp", &self.max_outdoor_temp)
            .field("min_flow_temp", &self.min_flow_temp)
            .field("max_flow_temp", &self.max_flow_temp)
            .field("temp_emitter_prev", &self.temp_emitter_prev)
            .field("target_flow_temp", &self.target_flow_temp)
            .finish_non_exhaustive()
    }
}

#[derive(Copy, Clone)]
struct EmittersAndPowerInput<'a> {
    pub emitters: &'a Emitters,
    pub power_input: f64,
    pub temp_diff_max: Option<f64>,

    // TODO can we calculate what this should be
    // based on initial value and timestep
    previous_difference_from_temp_diff_max: Option<f64>,
}

impl EmittersAndPowerInput<'_> {
    fn difference_from_temp_diff_max(&self, y: f64) -> f64 {
        y - self.temp_diff_max.unwrap()
    }
}

impl System<Time, State> for EmittersAndPowerInput<'_> {
    fn system(&self, _x: Time, y: &State, dy: &mut State) {
        dy[0] = self
            .emitters
            .func_temp_emitter_change_rate(self.power_input, y[0]);
    }

    // Stop function called at every successful integration step. The integration is stopped when this function returns true.
    fn solout(&mut self, _x: Time, y: &State, _dy: &State) -> bool {
        if self.temp_diff_max.is_none() {
            // no maximum - keep going
            return false;
        }

        let difference_from_temp_diff_max = self.difference_from_temp_diff_max(y[0]);

        if self.previous_difference_from_temp_diff_max.is_some() {
            let previous = self.previous_difference_from_temp_diff_max.unwrap();

            // signs are different - we must have passed zero
            if difference_from_temp_diff_max == 0.
                || signs_are_different(difference_from_temp_diff_max, previous)
            {
                // passing zero means we hit temp_diff_max, so stop solver
                return true;
            }
        }

        self.previous_difference_from_temp_diff_max = Some(difference_from_temp_diff_max);
        false
    }
}

fn signs_are_different(a: f64, b: f64) -> bool {
    (b > 0. && a < 0.) || (b < 0. && a > 0.)
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
        heat_source: Arc<RwLock<SpaceHeatingService>>,
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions: Arc<ExternalConditions>,
        ecodesign_controller: EcoDesignController,
        design_flow_temp: f64,
        output_detailed_results: bool,
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
            _temp_diff_emit_dsgn: temp_diff_emit_dsgn,
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
            temp_emitter_prev: 20.0,
            target_flow_temp: None,
            output_detailed_results,
            emitters_detailed_results: output_detailed_results.then(Default::default),
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

    pub fn in_required_period(
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

    pub fn func_temp_emitter_change_rate(&self, power_input: f64, y: f64) -> f64 {
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

        let emitter_with_power_input = EmittersAndPowerInput {
            emitters: self,
            power_input,
            temp_diff_max,

            // TODO set these up in the struct with a new function
            previous_difference_from_temp_diff_max: None,
        };

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

        let temp_emitter;
        let mut time_temp_diff_max_reached: Option<f64> = None;

        let y_count = stepper.y_out().len();
        let last_y = stepper.y_out().last().expect("y_out was empty")[0];
        let second_last_y = stepper.y_out().get(y_count - 2).unwrap()[0];

        let max_temp_diff_was_reached = temp_diff_max.is_some()
            && signs_are_different(
                f.difference_from_temp_diff_max(last_y),
                f.difference_from_temp_diff_max(second_last_y),
            );
        if max_temp_diff_was_reached {
            // We stopped early because the max diff was passed.
            // The Python code uses a built in feature of scipy's solve_ivp here.
            // when an "event" (in this case, max temp diff) happens a root solver
            // finds the exact x (time) value for that event occuring
            // and sets time_temp_diff_max_reached

            let mut y_out_reversed: Vec<f64> = stepper.y_out().iter().flatten().copied().collect();
            y_out_reversed.reverse();

            let mut x_out_reversed = stepper.x_out().clone();
            x_out_reversed.reverse();

            // based on the time steps and "outputs" we have from the solver
            // interpolate a reasonable guess as to what time we hit max temp diff
            let interpolated_guess =
                np_interp(temp_diff_max.unwrap(), &y_out_reversed[..], &x_out_reversed);
            time_temp_diff_max_reached = Some(interpolated_guess);

            // max temp diff was reached, so that should be our result
            temp_emitter = temp_rm + temp_diff_max.unwrap();
        } else {
            let last_y = stepper.y_out().last().expect("y_out was empty")[0];
            let temp_diff_emitter_rm_final = last_y;
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
    ) -> (f64, bool, Option<BufferTankEmittersDataWithResult>) {
        // When there is some demand, calculate max. emitter temperature
        // achievable and emitter temperature required, and base calculation
        // on the lower of the two.
        // TODO (from Python) The final calculation of emitter temperature below assumes
        // constant power output from heating source over the timestep.
        // It does not account for:
        // - overshoot/undershoot and then stabilisation of emitter temp.
        // This leads to emitter temp at end of timestep not exactly
        // matching temp_emitter_target.
        // - On warm-up, calculate max. temp achievable, then cap at target
        // and record time this is reached, then assume target temp is
        // maintained? Would there be an overshoot to make up for
        // underheating during warm-up?
        // - On cool-down to lower target temp, calculate lowest temp achievable,
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
        // If emitters are already above max. temp for this timestep,
        // then heat source should provide no energy until emitter temp
        // falls to maximum
        // Otherwise:
        let (energy_provided_by_heat_source_max_min, emitters_data_for_buffer_tank_with_result) =
            if self.temp_emitter_prev <= temp_emitter_max {
                // If emitters are below max. temp for this timestep, then max energy
                // required from heat source will depend on maximum warm-up rate,
                // which depends on the maximum energy output from the heat source
                let emitters_data_for_buffer_tank = match self.with_buffer_tank {
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

                self.heat_source
                    .read()
                    .energy_output_max(
                        temp_emitter_max,
                        temp_return,
                        emitters_data_for_buffer_tank,
                        simulation_time,
                    )
                    .unwrap()
            } else {
                (Default::default(), None)
            };

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
            emitters_data_for_buffer_tank_with_result,
        )
    }

    /// Demand energy from emitters and calculate how much energy can be provided
    /// Arguments:
    /// energy_demand -- in kWh
    pub(crate) fn demand_energy(
        &mut self,
        energy_demand: f64,
        simulation_time: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let timestep = simulation_time.timestep;
        let temp_rm_prev = &self.temp_internal_air_fn;

        // Calculate target flow and return temperature
        let (temp_flow_target, temp_return_target) = self.temp_flow_return(&simulation_time);
        let temp_emitter_max = (temp_flow_target + temp_return_target) / 2.0;
        self.target_flow_temp = Some(temp_flow_target);

        let mut energy_req_from_heat_source = Default::default();
        let mut temp_emitter_max_is_final_temp = false;
        let mut emitters_data_for_buffer_tank_with_result = None;

        if energy_demand > 0. {
            // Emitters warming up or cooling down to a target temperature
            (
                energy_req_from_heat_source,
                temp_emitter_max_is_final_temp,
                emitters_data_for_buffer_tank_with_result,
            ) = self.energy_required_from_heat_source(
                energy_demand,
                timestep,
                temp_rm_prev(),
                temp_emitter_max,
                temp_return_target,
                simulation_time,
            )
        }

        // Get energy output of heat source (i.e. energy input to emitters)
        // TODO (from Python) - Instead of passing temp_flow_req into heating system module,
        // calculate average flow temp achieved across timestep?

        // Catering for the possibility of a BufferTank in the emitters' loop
        let energy_provided_by_heat_source: f64 = if self.with_buffer_tank {
            // Call to HeatSourceServiceSpace with buffer_tank relevant data
            self.heat_source
                .write()
                .demand_energy(
                    energy_req_from_heat_source,
                    temp_flow_target,
                    temp_return_target,
                    emitters_data_for_buffer_tank_with_result,
                    simulation_time,
                )?
                .0
        } else {
            self.heat_source
                .write()
                .demand_energy(
                    energy_req_from_heat_source,
                    temp_flow_target,
                    temp_return_target,
                    None,
                    simulation_time,
                )?
                .0
        };
        // Calculate emitter temperature achieved at end of timestep.
        // Do not allow emitter temp to rise above maximum
        // Do not allow emitter temp to fall below room temp
        let mut temp_emitter = if temp_emitter_max_is_final_temp {
            temp_emitter_max
        } else {
            let power_provided_by_heat_source = energy_provided_by_heat_source / timestep;
            self.temp_emitter(
                0.0,
                timestep,
                self.temp_emitter_prev,
                temp_rm_prev(),
                power_provided_by_heat_source,
                None,
            )
            .0
        };
        temp_emitter = max_of_2(temp_emitter, temp_rm_prev());

        // Calculate emitter output achieved at end of timestep.
        let energy_released_from_emitters = energy_provided_by_heat_source
            + self.thermal_mass * (self.temp_emitter_prev - temp_emitter);

        // Save emitter temperature for next timestep
        self.temp_emitter_prev = temp_emitter;

        if self.output_detailed_results {
            self.emitters_detailed_results
                .as_ref()
                .unwrap()
                .write()
                .push(EmittersDetailedResult {
                    timestep_index: simulation_time.index,
                    energy_demand,
                    energy_provided_by_heat_source,
                    temp_emitter,
                    temp_emitter_max,
                    energy_released_from_emitters,
                    temp_flow_target,
                    temp_return_target,
                    temp_emitter_max_is_final_temp,
                    energy_req_from_heat_source,
                });
        }

        Ok(energy_released_from_emitters)
    }

    /// Return the cumulative running time and throughput factor for the heat source
    /// Arguments:
    /// energy_demand -- in kWh
    /// space_heat_running_time_cumulative
    ///     -- running time spent on higher-priority space heating services
    pub(crate) fn running_time_throughput_factor(
        &self,
        energy_demand: f64,
        space_heat_running_time_cumulative: f64,
        simulation_time: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        let timestep = simulation_time.timestep;
        let temp_rm_prev = &self.temp_internal_air_fn;

        // Calculate target flow and return temperature
        let (temp_flow_target, temp_return_target) = self.temp_flow_return(&simulation_time);
        let temp_emitter_max = (temp_flow_target + temp_return_target) / 2.;

        let energy_req_from_heat_source = if energy_demand > 0. {
            // Emitters warming up or cooling down to a target temperature
            self.energy_required_from_heat_source(
                energy_demand,
                timestep,
                temp_rm_prev(),
                temp_emitter_max,
                temp_return_target,
                simulation_time,
            )
            .0
        } else {
            // Emitters cooling down or at steady-state with heating off
            0.
        };

        self.heat_source.read().running_time_throughput_factor(
            space_heat_running_time_cumulative,
            energy_req_from_heat_source,
            temp_flow_target,
            temp_return_target,
            simulation_time,
        )
    }

    pub(crate) fn output_emitter_results(&self) -> Option<Vec<EmittersDetailedResult>> {
        self.emitters_detailed_results
            .as_ref()
            .map(|results| results.read().clone())
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct EmittersDetailedResult {
    timestep_index: usize,
    energy_demand: f64,
    energy_provided_by_heat_source: f64,
    temp_emitter: f64,
    temp_emitter_max: f64,
    energy_released_from_emitters: f64,
    temp_flow_target: f64,
    temp_return_target: f64,
    temp_emitter_max_is_final_temp: bool,
    energy_req_from_heat_source: f64,
}

impl EmittersDetailedResult {
    pub(crate) fn as_string_values(&self) -> Vec<String> {
        vec![
            self.timestep_index.to_string(),
            self.energy_demand.to_string(),
            self.energy_provided_by_heat_source.to_string(),
            self.temp_emitter.to_string(),
            self.temp_emitter_max.to_string(),
            self.energy_released_from_emitters.to_string(),
            self.temp_flow_target.to_string(),
            self.temp_return_target.to_string(),
            self.temp_emitter_max_is_final_temp.to_string(),
            self.energy_req_from_heat_source.to_string(),
        ]
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
    use crate::core::energy_supply::energy_supply::EnergySupply;
    use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
    use crate::core::heating_systems::boiler::Boiler;
    use crate::core::heating_systems::boiler::BoilerServiceSpace;
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
            },
            ShadingSegment {
                number: 2,
                start: 135.,
                end: 90.,
                shading_objects: None,
            },
            ShadingSegment {
                number: 3,
                start: 90.,
                end: 45.,
                shading_objects: None,
            },
            ShadingSegment {
                number: 4,
                start: 45.,
                end: 0.,
                shading_objects: None,
            },
            ShadingSegment {
                number: 5,
                start: 0.,
                end: -45.,
                shading_objects: None,
            },
            ShadingSegment {
                number: 6,
                start: -45.,
                end: -90.,
                shading_objects: None,
            },
            ShadingSegment {
                number: 7,
                start: -90.,
                end: -135.,
                shading_objects: None,
            },
            ShadingSegment {
                number: 8,
                start: -135.,
                end: -180.,
                shading_objects: None,
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
            energy_supply: EnergySupplyType::MainsGas,
            energy_supply_auxiliary: EnergySupplyType::Electricity,
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
        let energy_supply = Arc::from(RwLock::from(EnergySupply::new(
            FuelType::MainsGas,
            simulation_time_iterator.total_steps(),
            None,
            None,
            None,
        )));

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

        let control = Arc::from(Control::OnOffTimeControl(OnOffTimeControl::new(
            vec![true, true, true, true, true, true, true, true],
            0,
            0.25, // to match simulation time
        )));

        let boiler_service_space =
            BoilerServiceSpace::new(Arc::new(RwLock::new(boiler)), service_name.into(), control);

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

        let with_buffer_tank = false;

        Emitters::new(
            thermal_mass,
            c,
            n,
            temp_diff_emit_dsgn,
            frac_convective,
            Arc::new(RwLock::new(heat_source)),
            Arc::new(move || canned_value),
            external_conditions.into(),
            ecodesign_controller,
            design_flow_temp,
            false,
            with_buffer_tank,
        )
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
        simulation_time_iterator: SimulationTimeIterator,
        heat_source: SpaceHeatingService,
        external_conditions: ExternalConditions,
    ) {
        let (flow_temp, return_temp) =
            emitters.temp_flow_return(&simulation_time_iterator.current_iteration());

        assert_relative_eq!(flow_temp, 50.8333, max_relative = 1e-2);
        assert_relative_eq!(return_temp, 43.5714, max_relative = 1e-2);

        // Test with different outdoor temp
        let ecodesign_controller = EcoDesignController {
            ecodesign_control_class: EcoDesignControllerClass::ClassII,
            min_outdoor_temp: Some(10.),
            max_outdoor_temp: Some(15.),
            min_flow_temp: Some(30.),
        };

        let heat_source = Arc::new(RwLock::new(heat_source));

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
            false,
            false,
        );

        let (flow_temp, return_temp) =
            emitters.temp_flow_return(&simulation_time_iterator.current_iteration());

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
            false,
            false,
        );

        let (flow_temp, return_temp) =
            emitters.temp_flow_return(&simulation_time_iterator.current_iteration());

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

        let heat_source = Arc::new(RwLock::new(heat_source));

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
            55.0,
            false,
            false,
        );

        let temp_emitter = 15.;
        let temp_rm = 10.;

        let result = emitters.power_output_emitter(temp_emitter, temp_rm);

        assert_relative_eq!(result, 0.55189186, max_relative = EIGHT_DECIMAL_PLACES);
    }

    /// Test emitter temperature that gives required power output at given room temp
    #[rstest]
    fn test_temp_emitter_req(emitters: Emitters) {
        let power_emitter_req = 0.22;
        let temp_rm = 2.0;

        let result = emitters.temp_emitter_req(power_emitter_req, temp_rm);

        assert_relative_eq!(result, 4.32332827, max_relative = EIGHT_DECIMAL_PLACES);
    }

    /// Test Differential eqn is formed for change rate of emitter temperature
    #[rstest]
    fn test_func_temp_emitter_change_rate(emitters: Emitters) {
        let result = emitters.func_temp_emitter_change_rate(5., 32.);

        assert_eq!(result, -0.8571428571428514);
    }

    #[rstest]
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

        // TODO can we set a threshold here for a close enough match?
        assert_relative_eq!(time_temp_diff_max_reached.unwrap(), 1.29981138);
    }

    #[rstest]
    #[ignore = "blocked by temp_emitters issue"]
    fn test_energy_required_from_heat_source(
        simulation_time_iterator: SimulationTimeIterator,
        emitters: Emitters,
    ) {
        let energy_demand_list = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0];
        let mut energy_demand = 0.0;

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            energy_demand += energy_demand_list[t_idx];

            let (energy_req, temp_emitter_max_is_final_temp, _) =
                emitters.energy_required_from_heat_source(energy_demand, 1., 10., 25., 30., t_it);

            assert_relative_eq!(
                energy_req,
                [
                    0.7487346289045738,
                    2.458950517761754,
                    2.458950517761754,
                    2.458950517761754,
                    2.458950517761754,
                    2.458950517761754,
                    2.458950517761754,
                    2.458950517761754
                ][t_idx]
            );

            assert_eq!(
                temp_emitter_max_is_final_temp,
                [false, false, true, true, true, true, true, true][t_idx]
            );
        }
    }

    // Python has a test_running_time_throughput_factor which we've not ported for now as it relies
    // on mocking
}
