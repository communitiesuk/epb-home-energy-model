/// This module provides object(s) to model the behaviour of instantaneous electric
/// room heaters.
use crate::core::controls::time_control::SetpointOrCombinationControl;
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::heating_systems::constants::MAX_TEMPERATURE_TOUCHABLE;
use crate::core::solvers::root;
use crate::core::solvers::solve_ivp::{
    solve_ivp, OdeResult, SharedIvpSolveFunction, TerminatingEvent,
};
use crate::corpus::TempInternalAirFn;
use crate::simulation_time::SimulationTimeIteration;
use anyhow::{anyhow, bail};
use approx::relative_eq;
use atomic_float::AtomicF64;
use educe::Educe;
use ndarray::{array, Array1};
#[cfg(test)]
use parking_lot::RwLock;
use std::sync::atomic::Ordering;
use std::sync::Arc;
use tracing::warn;

/// Type to represent instantaneous electric heaters
#[derive(Educe)]
#[educe(Debug)]
pub struct InstantElecHeater {
    rated_power_in_kw: f64,
    frac_convective: f64,
    energy_supply_connection: EnergySupplyConnection,
    #[educe(Debug(ignore))]
    temp_internal_air_fn: TempInternalAirFn,
    simulation_timestep: f64,
    control: Option<SetpointOrCombinationControl>,
    thermal_mass: Option<ThermalMassFields>,
    temp_emitter_prev: AtomicF64,
    max_temperature: f64,
    #[cfg(test)]
    amounts_demanded_from_supply_connection: RwLock<Vec<f64>>,
}

#[derive(Clone, Copy, Debug)]
struct ThermalMassFields {
    thermal_mass: f64,
    c: f64,
    n: f64,
}

impl InstantElecHeater {
    /// Arguments
    /// * `rated_power` - in kW
    /// * `frac_convective` - convective fraction for heating
    /// * `energy_supply_connection` - EnergySupplyConnection value
    /// * `temp_internal_air_fn` - function that can provide the internal air temperature of the zone where the heater is located
    /// * `simulation_timestep` - step in hours for context simulation time
    /// * `control` - reference to a control object which must implement is_on() and setpnt() funcs
    /// * `c` - constant from characteristic equation of emitters (e.g. derived from BS EN 442 style tests)
    /// * `c_per_kw` - constant from characteristic equation of emitters (e.g. derived from BS EN 442 style tests) per kW of rated power
    /// * `n` - exponent from characteristic equation of emitters (e.g. derived from BS EN 442 style tests)
    /// * `thermal_mass` - thermal mass of heater in kWh/K
    /// * `thermal_mass_per_kw` - thermal mass of heater in kWh/K per kW of rated power
    /// * `initial_temperature`
    /// * `max_temperature`
    pub(crate) fn new(
        rated_power_in_kw: f64,
        frac_convective: f64,
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        simulation_timestep: f64,
        control: Option<SetpointOrCombinationControl>,
        c: Option<f64>,
        c_per_kw: Option<f64>,
        n: Option<f64>,
        thermal_mass: Option<f64>,
        thermal_mass_per_kw: Option<f64>,
        initial_temperature: f64,
        max_temperature: Option<f64>,
    ) -> anyhow::Result<Self> {
        let mut c = c;
        if let Some(c_per_kw) = c_per_kw {
            c = Some(c_per_kw * rated_power_in_kw);
        }

        let thermal_mass = thermal_mass_per_kw
            .map(|thermal_mass_per_kw| thermal_mass_per_kw * rated_power_in_kw)
            .or(thermal_mass);

        let thermal_mass = thermal_mass.map(|thermal_mass| {
            let c = if let Some(c) = c {
                c
            } else {
                bail!("Constant from characteristic equation of emitters missing for heater with thermal mass");
            };
            let n = if let Some(n) = n {
                n
            } else {
                bail!("Exponent from characteristic equation of emitters missing for heater with thermal mass");
            };

            Ok(ThermalMassFields { thermal_mass, c, n })
        }).transpose()?;

        Ok(Self {
            rated_power_in_kw,
            frac_convective,
            energy_supply_connection,
            temp_internal_air_fn,
            simulation_timestep,
            control,
            thermal_mass,
            temp_emitter_prev: AtomicF64::new(initial_temperature),
            max_temperature: max_temperature.unwrap_or(MAX_TEMPERATURE_TOUCHABLE),
            #[cfg(test)]
            amounts_demanded_from_supply_connection: Default::default(),
        })
    }

    pub fn temp_setpnt(&self, simtime: &SimulationTimeIteration) -> Option<f64> {
        self.control.as_ref().and_then(|ctrl| ctrl.setpnt(simtime))
    }

    pub fn in_required_period(&self, simtime: &SimulationTimeIteration) -> Option<bool> {
        self.control
            .as_ref()
            .and_then(|ctrl| ctrl.in_required_period(simtime))
    }

    pub fn frac_convective(&self) -> f64 {
        self.frac_convective
    }

    /// Calculate minimum possible energy output
    pub(crate) fn energy_output_min(&self) -> anyhow::Result<f64> {
        let (thermal_mass, thermal_mass_fields) =
            if let Some(thermal_mass_fields) = self.thermal_mass {
                (thermal_mass_fields.thermal_mass, thermal_mass_fields)
            } else {
                return Ok(0.);
            };

        let timestep = self.simulation_timestep;
        let temp_rm_prev = (self.temp_internal_air_fn)();

        let (temp_emitter, _) = self.temp_emitter(
            0.0,
            timestep,
            self.temp_emitter_prev.load(Ordering::SeqCst),
            temp_rm_prev,
            0.0,
            None,
            thermal_mass_fields,
        )?;
        let temp_emitter = temp_emitter.max(temp_rm_prev);

        // Calculate emitter output achieved at end of timestep.
        Ok(thermal_mass * (self.temp_emitter_prev.load(Ordering::SeqCst) - temp_emitter))
    }

    /// Demand energy (in kWh) from the heater
    pub fn demand_energy(
        &self,
        energy_demand: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let (thermal_mass, thermal_mass_fields) =
            if let Some(thermal_mass_fields) = self.thermal_mass {
                (thermal_mass_fields.thermal_mass, thermal_mass_fields)
            } else {
                return self.demand_energy_no_thermal_mass(energy_demand, simtime);
            };

        let timestep = self.simulation_timestep;
        let temp_rm_prev = (self.temp_internal_air_fn)();

        // The emitter coasts (no heat source input) when the requested net output
        // does not exceed the coasting minimum output, i.e. the output the emitter
        // produces with the heat source off. An emitter warmer than the room sheds
        // its stored heat (positive coasting minimum); one below the room absorbs
        // heat as its mass warms (negative coasting minimum). A demand at or below
        // this minimum is already met by coasting, so the heat source stays off; a
        // demand above it engages the heat source to supply the increment the zone
        // requires.
        let energy_coast_threshold = self.energy_output_min()?;

        let (
            time_heating_start,
            temp_emitter_heating_start,
            energy_req_from_heat_source,
            temp_emitter_max_is_final_temp,
            temp_emitter_req,
        ) = if energy_demand < energy_coast_threshold
            || relative_eq!(energy_demand, energy_coast_threshold, epsilon = 1e-10)
        {
            (
                0.0,
                self.temp_emitter_prev.load(Ordering::SeqCst),
                0.0,
                false,
                temp_rm_prev,
            )
        } else {
            // Calculate emitter temperature required
            let power_emitter_req = energy_demand / timestep;
            let temp_emitter_req = if power_emitter_req < 0.0 {
                // A negative net output target cannot be met by an emitter warmer than
                // the room, so the required emitter temperature is the room temperature
                // (giving zero steady-state output). This also avoids evaluating the
                // emitter power law for a negative output, which has no real solution
                // below room temperature. A zero target still solves to the room
                // temperature via temp_emitter_req, so is left to that path.
                temp_rm_prev
            } else {
                self.calculate_emitter_required_temperature(
                    power_emitter_req,
                    temp_rm_prev,
                    thermal_mass_fields,
                )
            };

            // Heater warming up or cooling down to a target temperature:
            // - First we calculate the time taken for the heaters to cool
            //   before the heating system activates, and the temperature that
            //   the heaters reach at this time. Note that the heaters will
            //   cool to below the target temperature so that the total heat
            //   output in this cooling period matches the demand accumulated so
            //   far in the timestep (assumed to be proportional to the fraction
            //   of the timestep that has elapsed)
            let (time_heating_start, temp_emitter_heating_start) = self.calc_emitter_cooldown(
                energy_demand,
                temp_emitter_req,
                temp_rm_prev,
                timestep,
                thermal_mass_fields,
            )?;

            // Then, we calculate the energy required from the heat source in
            // the remaining part of the timestep
            let (energy_req_from_heat_source, temp_emitter_max_is_final_temp) = self
                .energy_required_from_heat_source(
                    energy_demand * (1.0 - time_heating_start / timestep),
                    time_heating_start,
                    timestep,
                    temp_rm_prev,
                    temp_emitter_heating_start,
                    temp_emitter_req,
                    self.max_temperature,
                    thermal_mass_fields,
                )?;

            (
                time_heating_start,
                temp_emitter_heating_start,
                energy_req_from_heat_source,
                temp_emitter_max_is_final_temp,
                temp_emitter_req,
            )
        };

        let energy_provided_by_heat_source = energy_req_from_heat_source
            .min(self.rated_power_in_kw * (timestep - time_heating_start));

        // Calculate heater temperature achieved at end of timestep.
        // Do not allow heater temp to rise above maximum
        // Do not allow heater temp to fall below room temp
        let temp_emitter = if temp_emitter_max_is_final_temp {
            self.max_temperature
        } else {
            let power_provided_by_heat_source =
                energy_provided_by_heat_source / (timestep - time_heating_start);
            let (temp_emitter, time_temp_target_reached) = self.temp_emitter(
                time_heating_start,
                timestep,
                temp_emitter_heating_start,
                temp_rm_prev,
                power_provided_by_heat_source,
                temp_emitter_req.into(),
                thermal_mass_fields,
            )?;

            if temp_emitter_heating_start < temp_emitter_req && time_temp_target_reached.is_some() {
                temp_emitter_req
            } else {
                temp_emitter
            }
        };

        let temp_emitter = temp_emitter.max(temp_rm_prev);

        // Calculate heater output achieved at end of timestep.
        let energy_released_from_emitters = energy_provided_by_heat_source
            + thermal_mass * (self.temp_emitter_prev.load(Ordering::SeqCst) - temp_emitter);

        self.temp_emitter_prev.store(temp_emitter, Ordering::SeqCst);

        // log amount passed to energy supply connection in test environment
        #[cfg(test)]
        {
            self.amounts_demanded_from_supply_connection
                .write()
                .push(energy_provided_by_heat_source);
        }

        self.energy_supply_connection
            .demand_energy(energy_provided_by_heat_source, simtime.index)?;

        Ok(energy_released_from_emitters)
    }

    #[cfg(test)]
    fn passed_energy_demand_values(&self) -> Vec<f64> {
        self.amounts_demanded_from_supply_connection.read().clone()
    }

    /// Demand energy (in kWh) from the heater with no thermal mass
    fn demand_energy_no_thermal_mass(
        &self,
        energy_demand: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let energy_supplied = if self.control.as_ref().is_none_or(|c| c.is_on(&simtime)) {
            energy_demand.min(self.rated_power_in_kw * self.simulation_timestep)
        } else {
            0.0
        };

        // log amount passed to energy supply connection in test environment
        #[cfg(test)]
        {
            self.amounts_demanded_from_supply_connection
                .write()
                .push(energy_supplied);
        }

        self.energy_supply_connection
            .demand_energy(energy_supplied, simtime.index)?;

        Ok(energy_supplied)
    }

    /// Calculate emitter temperature that gives required power output at given room temp
    ///
    ///        Power output from emitter (eqn from 2020 ASHRAE Handbook p644):
    ///           power_output = c * (T_E - T_rm) ^ n
    ///       where:
    ///            T_E is mean emitter temperature
    ///            T_rm is air temperature in the room/zone
    ///            c and n are characteristic of the emitters (e.g. derived from BS EN 442 tests)
    ///        Rearrange to solve for T_E
    fn calculate_emitter_required_temperature(
        &self,
        power_emitter_req: f64,
        temp_rm: f64,
        ThermalMassFields { c, n, .. }: ThermalMassFields,
    ) -> f64 {
        (power_emitter_req / c).powf(1.0 / n) + temp_rm
    }

    /// Differential eqn for change rate of emitter temperature, to be solved iteratively
    fn func_temp_emitter_change_rate(
        power_input: f64,
        ThermalMassFields { c, n, thermal_mass }: ThermalMassFields,
    ) -> SharedIvpSolveFunction {
        // Heat balance equation for radiators:
        //             (T_E(t) - T_E(t-1)) * K_E / timestep = power_input - power_output
        //         where:
        //             T_E is mean emitter temperature
        //             K_E is thermal mass of emitters
        //
        //         Power output from emitter (eqn from 2020 ASHRAE Handbook p644):
        //             power_output = c * (T_E(t) - T_rm) ^ n
        //         where:
        //             T_rm is air temperature in the room/zone
        //             c and n are characteristic of the emitters (e.g. derived from BS EN 442 tests)
        //
        //         Substituting power output eqn into heat balance eqn gives:
        //             (T_E(t) - T_E(t-1)) * K_E / timestep = power_input - c * (T_E(t) - T_rm) ^ n
        //
        //         Rearranging gives:
        //             (T_E(t) - T_E(t-1)) / timestep = (power_input - c * (T_E(t) - T_rm) ^ n) / K_E
        //         which gives the differential equation as timestep goes to zero:
        //             d(T_E)/dt = (power_input - c * (T_E - T_rm) ^ n) / K_E
        //
        //         If T_rm is assumed to be constant over the time period, then the rate of
        //         change of T_E is the same as the rate of change of deltaT, where:
        //             deltaT = T_E - T_rm
        //
        //         Therefore, the differential eqn can be expressed in terms of deltaT:
        //             d(deltaT)/dt = (power_input - c * deltaT(t) ^ n) / K_E
        //
        //         This can be solved for deltaT over a specified time period using the
        //         solve_ivp function from scipy (or equivalent in Rust!).

        // Apply min value of zero to temp_diff because the power law does not
        // work for negative temperature difference
        // solve for temp_diff iteratively
        Arc::new(move |_t, temp_diff| {
            array![(power_input - c * 0.0f64.max(temp_diff[0]).powf(n)) / thermal_mass]
        })
    }

    /// Calculate emitter output at given emitter and room temp
    ///
    ///        Power output from emitter (eqn from 2020 ASHRAE Handbook p644):
    ///            power_output = c * (T_E - T_rm) ^ n
    ///        where:
    ///            T_E is mean emitter temperature
    ///            T_rm is air temperature in the room/zone
    ///            c and n are characteristic of the emitters (e.g. derived from BS EN 442 tests)
    fn power_output_emitter(
        &self,
        temp_emitter: f64,
        temp_rm: f64,
        ThermalMassFields { c, n, .. }: ThermalMassFields,
    ) -> f64 {
        c * 0.0f64.max((temp_emitter - temp_rm).powf(n))
    }

    /// Calculate emitter cooling time and emitter temperature at this time
    fn calc_emitter_cooldown(
        &self,
        energy_demand: f64,
        temp_emitter_req: f64,
        temp_rm_prev: f64,
        timestep: f64,
        thermal_mass_fields: ThermalMassFields,
    ) -> anyhow::Result<(f64, f64)> {
        let temp_emitter_prev = self.temp_emitter_prev.load(Ordering::SeqCst);

        Ok(if temp_emitter_prev < temp_emitter_req {
            (0.0, temp_emitter_prev)
        } else {
            // Calculate time that emitters are cooling down (accounting for
            // undershoot), during which the heat source does not provide any
            // heat, by iterating to find the end time which leads to the heat
            // output matching the energy demand accumulated so far during the
            // timestep
            let energy_surplus_during_cooldown_fn = Box::new(
                |time_cooldown, [timestep, energy_demand, temp_rm_prev]: [f64; 3]| {
                    self.energy_surplus_during_cooldown(
                        time_cooldown,
                        timestep,
                        energy_demand,
                        temp_rm_prev,
                        thermal_mass_fields,
                    )
                },
            );

            let time_cooldown = root(
                energy_surplus_during_cooldown_fn,
                timestep,
                [timestep, energy_demand, temp_rm_prev],
            )?;

            // Limit cooldown time to be within timestep
            let time_heating_start = 0.0f64.max(time_cooldown.min(timestep));

            // Calculate emitter temperature at heating start time
            let (temp_emitter_heating_start, _) = self.temp_emitter(
                0.0,
                time_heating_start,
                temp_emitter_prev,
                temp_rm_prev,
                0.0, // No heat from heat source during initial cool-down
                None,
                thermal_mass_fields,
            )?;

            (time_heating_start, temp_emitter_heating_start)
        })
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
        thermal_mass_fields: ThermalMassFields,
    ) -> anyhow::Result<(f64, Option<f64>)> {
        // Calculate emitter temp at start of timestep
        let temp_diff_start = temp_emitter_start - temp_rm;

        let events = temp_emitter_max.map(|temp_emitter_max| {
            let temp_diff_max = temp_emitter_max - temp_rm;

            [TerminatingEvent::new(
                Arc::new(move |_t: f64, y: &Array1<f64>| -> f64 { y[0] - temp_diff_max }),
                None,
            )]
        });

        // Get function representing change rate equation and solve iteratively
        let func_temp_emitter_change_rate: SharedIvpSolveFunction =
            Self::func_temp_emitter_change_rate(power_input, thermal_mass_fields);

        let temp_diff_emitter_rm_results = solve_ivp(
            &func_temp_emitter_change_rate,
            (time_start, time_end),
            &array![temp_diff_start],
            events.as_ref().map(|events| events.as_slice()),
            None,
            None,
        )?;

        let OdeResult { t_events, y, .. } = temp_diff_emitter_rm_results;

        // Get time at which emitters reach max. temp
        let time_temp_diff_max_reached: Option<f64> = if let Some(ref t_events) = t_events {
            temp_emitter_max.and_then(|_| {
                t_events
                    .first()
                    .and_then(|t_events| t_events.iter().copied().last())
            })
        } else {
            None
        };

        // Get emitter temp at end of timestep
        let temp_diff_emitter_rm_final = y
            .last()
            .ok_or_else(|| anyhow!("A non-empty y array from solve_ivp had no available values"))?
            [0];
        let temp_emitter = temp_rm + temp_diff_emitter_rm_final;

        Ok((temp_emitter, time_temp_diff_max_reached))
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
        thermal_mass_fields: ThermalMassFields,
    ) -> anyhow::Result<(f64, bool)> {
        // When there is some demand, calculate max. emitter temperature
        // achievable and emitter temperature required, and base calculation
        // on the lower of the two.

        let ThermalMassFields { thermal_mass, .. } = thermal_mass_fields;

        // Calculate extra energy required for emitters to reach temp required
        let energy_req_to_warm_emitters =
            thermal_mass * (temp_emitter_req - temp_emitter_heating_start);

        // Calculate energy input required to meet energy demand
        let energy_req_from_heat_source =
            (energy_req_to_warm_emitters + energy_demand_heating_period).max(0.0);

        let energy_provided_by_heat_source_max_min =
            if temp_emitter_heating_start <= temp_emitter_max {
                // If emitters are below max. temp for this timestep, then max energy
                // required from heat source will depend on maximum warm-up rate,
                // which depends on the maximum energy output from the heat source
                energy_req_from_heat_source
                    .min(self.rated_power_in_kw * (timestep - time_heating_start))
            } else {
                // If emitters are already above max. temp for this timestep,
                // then heat source should provide no energy until emitter temp
                // falls to maximum
                0.0
            };

        // Calculate time to reach max. emitter temp at max heat source output
        let power_output_max_min = energy_provided_by_heat_source_max_min / timestep;
        let (temp_emitter, time_temp_emitter_max_reached) = self.temp_emitter(
            time_heating_start,
            timestep,
            temp_emitter_heating_start,
            temp_rm_prev,
            power_output_max_min,
            temp_emitter_max.into(),
            thermal_mass_fields,
        )?;

        let (time_in_warmup_cooldown_phase, temp_emitter_max_reached) =
            if let Some(time_temp_emitter_max_reached) = time_temp_emitter_max_reached {
                (time_temp_emitter_max_reached - time_heating_start, true)
            } else {
                (timestep - time_heating_start, false)
            };

        // Before this time, energy output from heat source is maximum
        let energy_req_from_heat_source_before_temp_emitter_max_reached =
            power_output_max_min * time_in_warmup_cooldown_phase;

        // After this time, energy output is amount needed to maintain
        // emitter temp (based on emitter output at constant emitter temp)
        // Note: the time at steady state in the equation below is the time
        //        remaining after the heating start and warmup/cooldown period
        //        and equals either:
        //        - zero, when time_temp_emitter_max_reached is None
        //        - (timestep - time_temp_emitter_max_reached), for other cases
        let energy_req_from_heat_source_after_temp_emitter_max_reached =
            self.power_output_emitter(temp_emitter, temp_rm_prev, thermal_mass_fields)
                * (timestep - time_heating_start - time_in_warmup_cooldown_phase);

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
        Ok((
            energy_req_from_heat_source.min(energy_req_from_heat_source_max),
            temp_emitter_max_is_final_temp,
        ))
    }

    fn energy_surplus_during_cooldown(
        &self,
        time_cooldown: f64,
        timestep: f64,
        energy_demand: f64,
        temp_rm_prev: f64,
        thermal_mass_fields: ThermalMassFields,
    ) -> f64 {
        // Calculate emitter temperature after specified time with no heat input
        let (temp_emitter_no_heat_input, _) = self
            .temp_emitter(
                0.0,
                time_cooldown,
                self.temp_emitter_prev.load(Ordering::SeqCst),
                temp_rm_prev,
                0.0, // No heat from heat source during initial cool-down
                None,
                thermal_mass_fields,
            )
            .unwrap_or_else(|e| {
                // NB. it may be that warning is not the right thing to do here
                warn!("The temp_emitter function in the instant_elec_heater module errored when not expected: {e}");

                Default::default() // don't panic - reporting no energy surplus seems preferable here
            });

        let ThermalMassFields { thermal_mass, .. } = thermal_mass_fields;

        let energy_released_from_emitters = thermal_mass
            * (self.temp_emitter_prev.load(Ordering::SeqCst) - temp_emitter_no_heat_input);
        let energy_demand_cooldown = energy_demand * time_cooldown / timestep;

        energy_released_from_emitters - energy_demand_cooldown
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::controls::time_control::SetpointTimeControl;
    use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyBuilder};
    use crate::input::FuelType;
    use crate::simulation_time::SimulationTime;
    use approx::assert_relative_eq;
    use parking_lot::RwLock;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use std::sync::Arc;

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 4., 1.)
    }

    fn create_temp_internal_air_fn(canned_value: f64) -> TempInternalAirFn {
        Arc::new(move || canned_value)
    }

    #[fixture]
    fn energy_supply_connection(simulation_time: SimulationTime) -> EnergySupplyConnection {
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.iter().total_steps())
                .build(),
        ));
        EnergySupply::connection(energy_supply, "shower").unwrap()
    }

    #[fixture]
    fn temp_internal_air_fn() -> TempInternalAirFn {
        create_temp_internal_air_fn(20.)
    }

    #[fixture]
    fn control(simulation_time: SimulationTime) -> SetpointOrCombinationControl {
        SetpointOrCombinationControl::SetpointTime(
            SetpointTimeControl::new(
                vec![Some(21.0), Some(21.0), None, Some(21.0)],
                0,
                1.,
                None,
                None,
                simulation_time.step,
            )
            .into(),
        )
    }

    #[fixture]
    fn instant_elec_heater(
        simulation_time: SimulationTime,
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        control: SetpointOrCombinationControl,
    ) -> InstantElecHeater {
        InstantElecHeater::new(
            50.,
            0.4,
            energy_supply_connection,
            temp_internal_air_fn,
            simulation_time.step,
            Some(control),
            1.3.into(),
            None,
            1.2.into(),
            0.14.into(),
            None,
            20.0,
            None,
        )
        .unwrap()
    }

    #[rstest]
    fn test_missing_parameters(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        simulation_time: SimulationTime,
        control: SetpointOrCombinationControl,
    ) {
        assert!(
            InstantElecHeater::new(
                50.,
                0.4,
                energy_supply_connection.clone(),
                temp_internal_air_fn.clone(),
                simulation_time.step,
                Some(control.clone()),
                None,
                None,
                1.2.into(),
                0.14.into(),
                None,
                20.0,
                None
            )
            .is_err(),
            "InstantElecHeater with no c when thermal_mass supplied should error"
        );

        assert!(
            InstantElecHeater::new(
                50.,
                0.4,
                energy_supply_connection,
                temp_internal_air_fn,
                simulation_time.step,
                Some(control),
                Some(1.),
                None,
                None,
                0.14.into(),
                None,
                20.0,
                None
            )
            .is_err(),
            "InstantElecHeater with no n when thermal_mass supplied should error"
        );
    }

    #[rstest]
    /// Test that the energy supplied is the same with setting either c or c_per_kw
    fn test_c_per_kw(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        simulation_time: SimulationTime,
        control: SetpointOrCombinationControl,
    ) {
        let inselecheater1 = InstantElecHeater::new(
            50.,
            0.4,
            energy_supply_connection.clone(),
            temp_internal_air_fn.clone(),
            simulation_time.step,
            Some(control.clone()),
            None,
            Some(1.3 / 50.),
            Some(1.2),
            Some(0.14),
            None,
            20.,
            None,
        )
        .unwrap();

        let inselecheater2 = InstantElecHeater::new(
            50.,
            0.4,
            energy_supply_connection.clone(),
            temp_internal_air_fn.clone(),
            simulation_time.step,
            Some(control.clone()),
            Some(1.3),
            None,
            Some(1.2),
            Some(0.14),
            None,
            20.,
            None,
        )
        .unwrap();

        let simtime = simulation_time.iter().next().unwrap();

        assert_eq!(
            inselecheater1.demand_energy(20., simtime).unwrap(),
            inselecheater2.demand_energy(20., simtime).unwrap()
        );
    }

    #[rstest]
    /// Test that the energy supplied is the same with setting either thermal_mass or thermal_mass_per_kw
    fn test_thermal_mass_per_kw(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        simulation_time: SimulationTime,
        control: SetpointOrCombinationControl,
    ) {
        let inselecheater1 = InstantElecHeater::new(
            50.,
            0.4,
            energy_supply_connection.clone(),
            temp_internal_air_fn.clone(),
            simulation_time.step,
            Some(control.clone()),
            Some(1.3),
            None,
            Some(1.2),
            Some(0.14),
            None,
            20.,
            None,
        )
        .unwrap();

        let inselecheater2 = InstantElecHeater::new(
            50.,
            0.4,
            energy_supply_connection.clone(),
            temp_internal_air_fn.clone(),
            simulation_time.step,
            Some(control.clone()),
            Some(1.3),
            None,
            Some(1.2),
            Some(0.14),
            Some(0.14 / 50.),
            20.,
            None,
        )
        .unwrap();

        let simtime = simulation_time.iter().next().unwrap();

        assert_eq!(
            inselecheater1.demand_energy(20., simtime).unwrap(),
            inselecheater2.demand_energy(20., simtime).unwrap()
        );
    }

    #[rstest]
    /// Test that InstantElecHeater object returns correct energy supplied
    fn test_demand_energy(instant_elec_heater: InstantElecHeater, simulation_time: SimulationTime) {
        let energy_input = [40.0, 100.0, 30.0, 20.0];
        let demand_expected = [
            40.0,
            49.50293629039578,
            28.368450790094787,
            19.166164548593276,
        ];
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                instant_elec_heater
                    .demand_energy(energy_input[t_idx], t_it)
                    .unwrap(),
                demand_expected[t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    /// Test that InstantElecHeater object returns correct energy supplied with a large demand
    fn test_demand_energy_large_demand(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        simulation_time: SimulationTime,
        control: SetpointOrCombinationControl,
    ) {
        let instant_elec_heater = InstantElecHeater::new(
            200.,
            0.4,
            energy_supply_connection,
            temp_internal_air_fn,
            simulation_time.step,
            control.into(),
            Some(1.3),
            None,
            Some(1.2),
            Some(1.),
            None,
            20.,
            None,
        )
        .unwrap();

        let energy_input = [200.0, 50.0, 0.0, 30.0];
        let demand_expected = [
            125.41615257674684,
            49.38851489989834,
            4.473591796719102,
            30.0,
        ];
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                instant_elec_heater
                    .demand_energy(energy_input[t_idx], t_it)
                    .unwrap(),
                demand_expected[t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    /// Test that InstantElecHeater with no thermal mass returns correct energy supplied
    fn test_demand_energy_no_thermal_mass(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        simulation_time: SimulationTime,
        control: SetpointOrCombinationControl,
    ) {
        let instant_elec_heater = InstantElecHeater::new(
            50.,
            0.4,
            energy_supply_connection,
            temp_internal_air_fn,
            simulation_time.step,
            control.into(),
            Some(1.3),
            None,
            Some(1.2),
            None,
            None,
            20.,
            None,
        )
        .unwrap();

        let energy_input = [40.0, 100.0, 30.0, 20.0];
        let demand_expected = [40., 50., 0., 20.];
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                instant_elec_heater
                    .demand_energy(energy_input[t_idx], t_it)
                    .unwrap(),
                demand_expected[t_idx]
            );
        }
    }

    #[rstest]
    /// Test that zero demand draws heat-source energy to maintain emitter
    /// temperature when the emitter ended the previous timestep below the
    /// current room temperature.
    fn test_demand_energy_zero_demand_emitter_below_room(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        simulation_time: SimulationTime,
        control: SetpointOrCombinationControl,
    ) {
        let instant_elec_heater = InstantElecHeater::new(
            10.,
            0.4,
            energy_supply_connection.clone(),
            temp_internal_air_fn,
            simulation_time.step,
            control.into(),
            Some(0.21),
            None,
            Some(1.02),
            Some(0.15),
            None,
            15.,
            None,
        )
        .unwrap();

        let simtime = simulation_time.iter().next().unwrap();

        // Room is 20 (set in setUp); emitter starts at 15 < 20.
        let energy_output = instant_elec_heater.demand_energy(0.0, simtime).unwrap();

        // Heat source provides energy to maintain emitter at room temperature.
        assert_eq!(instant_elec_heater.passed_energy_demand_values()[0], 0.75);

        assert_relative_eq!(energy_output, 0.);
    }

    #[rstest]
    /// Test that minimum demand draws no heat-source energy when the
    /// emitter ended the previous timestep below the current room temperature.
    fn test_demand_energy_min_demand_emitter_below_room(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        simulation_time: SimulationTime,
        control: SetpointOrCombinationControl,
    ) {
        let instant_elec_heater = InstantElecHeater::new(
            10.,
            0.4,
            energy_supply_connection,
            temp_internal_air_fn,
            simulation_time.step,
            control.into(),
            Some(0.21),
            None,
            Some(1.02),
            Some(0.15),
            None,
            15.0,
            None,
        )
        .unwrap();

        let simtime = simulation_time.iter().next().unwrap();

        // Room is 20 (set in setUp); emitter starts at 15 < 20.
        // Set demand to minimum output rather than zero, as negative demand can
        // still trigger heat source operation if it is above the min output.
        let energy_coast_threshold = instant_elec_heater.energy_output_min().unwrap();
        let energy_output = instant_elec_heater
            .demand_energy(energy_coast_threshold, simtime)
            .unwrap();

        // Heat source provides nothing — there is no demand
        assert_eq!(instant_elec_heater.passed_energy_demand_values()[0], 0.);

        // Emitter is colder than the room, so it absorbs heat from the room as
        // it warms to room temperature; the model represents this as a negative
        // heat output to the zone (mirrors Emitters.demand_energy_flow_return).
        // Energy released = thermal_mass * (temp_emitter_prev - temp_emitter_end)
        //                = 0.15 * (15 - 20) = -0.75 kWh.
        assert_relative_eq!(energy_output, -0.75);
    }

    #[rstest]
    fn test_demand_energy_negative_demand_emitter_below_room(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        simulation_time: SimulationTime,
        control: SetpointOrCombinationControl,
    ) {
        let instant_elec_heater = InstantElecHeater::new(
            10.,
            0.4,
            energy_supply_connection,
            temp_internal_air_fn,
            simulation_time.step,
            control.into(),
            Some(0.21),
            None,
            Some(1.02),
            Some(0.15),
            None,
            15.0,
            None,
        )
        .unwrap();

        let simtime = simulation_time.iter().next().unwrap();

        // Room is 20 (set in setUp); emitter starts at 15 < 20.
        // Set demand to half of minimum output to construct a demand value which is less
        // negative than the minimum output
        let energy_coast_threshold = instant_elec_heater.energy_output_min().unwrap();
        let energy_demand = energy_coast_threshold / 2.0;
        let energy_output = instant_elec_heater
            .demand_energy(energy_demand, simtime)
            .unwrap();
        // Heat source provides nothing — there is no demand.
        assert_eq!(instant_elec_heater.passed_energy_demand_values()[0], 0.375);
        // Emitter is colder than the room, so it absorbs heat from the room as
        // it warms to room temperature; the model represents this as a negative
        // heat output to the zone (mirrors Emitters.demand_energy_flow_return).
        // Emitter only warms to 17.5 rather than 20 because heat source activates
        // to meet negative demand above the minimum output (i.e. to reduce the
        // cooling effect of the emitters).
        // Energy released = thermal_mass * (temp_emitter_prev - temp_emitter_end)
        //                 = 0.15 * (15 - 17.5) = -0.75 kWh.
        assert_eq!(energy_output, -0.375);
    }

    /// Test that zero demand draws no heat-source energy when the emitter
    /// ended the previous timestep above the current room temperature.
    ///
    /// This is the standard cooldown case. Prior to the short-circuit, the
    /// root-finder in __calc_emitter_cooldown was invoked with an objective
    /// that is identically zero when energy_demand is zero, leaving its result
    /// undefined. After the short-circuit the ODE simply integrates a passive
    /// cooldown over the timestep.
    #[rstest]
    fn test_demand_energy_zero_demand_emitter_above_room(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        simulation_time: SimulationTime,
        control: SetpointOrCombinationControl,
    ) {
        let instant_elec_heater = InstantElecHeater::new(
            10.,
            0.4,
            energy_supply_connection,
            temp_internal_air_fn,
            simulation_time.step,
            control.into(),
            Some(0.21),
            None,
            Some(1.02),
            Some(0.15),
            None,
            25.0,
            None,
        )
        .unwrap();

        let simtime = simulation_time.iter().next().unwrap();

        // Room is 20 (set in setUp); emitter starts at 25 > 20. With no demand
        // and no input power, solve_ivp integrates the cooldown ODE over the
        // 1 h timestep, yielding temp_emitter_end = 21.2033... °C.
        // Energy released = thermal_mass * (temp_emitter_prev - temp_emitter_end)
        //                 = 0.15 * (25 - 21.2033...) = 0.5694955734921485 kWh.
        let energy_output = instant_elec_heater.demand_energy(0.0, simtime).unwrap();
        assert_eq!(instant_elec_heater.passed_energy_demand_values()[0], 0.);
        assert_relative_eq!(energy_output, 0.5694955734921485)
    }

    #[rstest]
    fn test_demand_energy_negative_demand(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        simulation_time: SimulationTime,
        control: SetpointOrCombinationControl,
    ) {
        let instant_elec_heater = InstantElecHeater::new(
            10.,
            0.4,
            energy_supply_connection,
            temp_internal_air_fn,
            simulation_time.step,
            control.into(),
            Some(0.21),
            None,
            Some(1.02),
            Some(0.15),
            None,
            22.0,
            None,
        )
        .unwrap();

        let simtime = simulation_time.iter().next().unwrap();

        // Room is 20 (set in setUp); emitter starts at 22 > 20. Cooldown ODE
        // over the 1 h timestep gives temp_emitter_end = 20.4938... °C.
        // Energy released = 0.15 * (22 - 20.4938...) = 0.22593164333702534 kWh.
        let energy_output = instant_elec_heater.demand_energy(-1.0e-6, simtime).unwrap();
        assert_eq!(instant_elec_heater.passed_energy_demand_values()[0], 0.);
        assert_relative_eq!(energy_output, 0.22593164333702534);
    }

    #[rstest]
    fn test_temp_setpnt(instant_elec_heater: InstantElecHeater, simulation_time: SimulationTime) {
        let setpoint_expected = [Some(21.0), Some(21.0), None, Some(21.0)];
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                instant_elec_heater.temp_setpnt(&t_it),
                setpoint_expected[t_idx]
            );
        }
    }

    #[rstest]
    fn test_in_required_period(
        instant_elec_heater: InstantElecHeater,
        simulation_time: SimulationTime,
    ) {
        let expected_whether = [true, true, false, true];
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                instant_elec_heater.in_required_period(&t_it),
                Some(expected_whether[t_idx])
            );
        }
    }

    #[rstest]
    fn test_frac_convective(instant_elec_heater: InstantElecHeater) {
        assert_eq!(instant_elec_heater.frac_convective(), 0.4);
    }

    #[rstest]
    fn test_energy_output_min(instant_elec_heater: InstantElecHeater) {
        assert_eq!(instant_elec_heater.energy_output_min().unwrap(), 0.0);
    }

    #[rstest]
    /// Test that the min energy output is above 0 after a previous demand
    fn test_energy_output_min_after_demand(
        instant_elec_heater: InstantElecHeater,
        simulation_time: SimulationTime,
    ) {
        instant_elec_heater
            .demand_energy(100., simulation_time.iter().next().unwrap())
            .unwrap();

        assert_relative_eq!(
            instant_elec_heater.energy_output_min().unwrap(),
            2.9280172879924646
        );
    }

    #[rstest]
    /// Test that the min energy output is 0 with no thermal mass
    fn test_energy_output_min_no_thermal_mass(
        energy_supply_connection: EnergySupplyConnection,
        temp_internal_air_fn: TempInternalAirFn,
        simulation_time: SimulationTime,
        control: SetpointOrCombinationControl,
    ) {
        let instant_elec_heater = InstantElecHeater::new(
            50.,
            0.4,
            energy_supply_connection,
            temp_internal_air_fn,
            simulation_time.step,
            Some(control),
            Some(1.3),
            None,
            None,
            None,
            None,
            20.,
            None,
        )
        .unwrap();

        assert_eq!(instant_elec_heater.energy_output_min().unwrap(), 0.0);
    }

    // the following tests in the upstream Python:
    //
    // test_func_temp_emitter_change_rate_n_is_none
    // test_temp_emitter_req_c_is_none
    // test_temp_emitter_req_n_is_none
    // test_energy_surplus_during_cooldown_n_is_none
    // test_energy_required_from_heat_source_thermal_mass_is_none
    // test_power_output_emitter
    //
    // are all redundant in Rust as the typing used makes case impossible to represent
    // (all tested methods need a ThermalMassFields struct, which is not available when thermal mass is not set)

    #[rstest]
    fn test_energy_required_from_heat_source_above_max(instant_elec_heater: InstantElecHeater) {
        let (energy_input, _) = instant_elec_heater
            .energy_required_from_heat_source(
                10.,
                0.,
                1.,
                20.,
                80.,
                40.,
                75.,
                instant_elec_heater.thermal_mass.unwrap(),
            )
            .unwrap();

        assert_relative_eq!(energy_input, 4.4);
    }

    // test_temp_emitter_exception is low value and not replicated here
}
