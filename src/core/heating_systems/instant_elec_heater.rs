/// This module provides object(s) to model the behaviour of instantaneous electric
/// room heaters.
use crate::core::controls::time_control::{per_control, Control, ControlBehaviour};
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::heating_systems::constants::MAX_TEMPERATURE_TOUCHABLE;
use crate::corpus::TempInternalAirFn;
use crate::simulation_time::SimulationTimeIteration;
use anyhow::bail;
use approx::relative_eq;
use atomic_float::AtomicF64;
use educe::Educe;
use std::sync::atomic::Ordering;

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
    control: Option<Control>,
    thermal_mass: Option<ThermalMassFields>,
    thermal_mass_per_kw: Option<f64>,
    temp_emitter_prev: AtomicF64,
    max_temperature: f64,
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
        control: Option<Control>,
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
            thermal_mass_per_kw,
            temp_emitter_prev: AtomicF64::new(initial_temperature),
            max_temperature: max_temperature.unwrap_or(MAX_TEMPERATURE_TOUCHABLE),
        })
    }

    pub fn temp_setpnt(&self, simulation_time_iteration: &SimulationTimeIteration) -> Option<f64> {
        self.control.as_ref().and_then(
            |ctrl| per_control!(&ctrl, ctrl => { ctrl.setpnt(simulation_time_iteration) }),
        )
    }

    pub fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool> {
        self.control.as_ref().and_then(|ctrl| per_control!(&ctrl, ctrl => { ctrl.in_required_period(simulation_time_iteration) }))
    }

    pub fn frac_convective(&self) -> f64 {
        self.frac_convective
    }

    /// Calculate minimum possible energy output
    pub(crate) fn energy_output_min(&self) -> f64 {
        let thermal_mass = if let Some(ThermalMassFields { thermal_mass, .. }) = self.thermal_mass {
            thermal_mass
        } else {
            return 0.;
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
        );
        let temp_emitter = temp_emitter.max(temp_rm_prev);

        // Calculate emitter output achieved at end of timestep.
        thermal_mass * (self.temp_emitter_prev.load(Ordering::SeqCst) - temp_emitter)
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
        let energy_coast_threshold = self.energy_output_min();

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
            let (time_heating_start, temp_emitter_heating_start) =
                self.calc_emitter_cooldown(energy_demand, temp_emitter_req, temp_rm_prev, timestep);

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
                );

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
            );

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

        self.energy_supply_connection
            .demand_energy(energy_provided_by_heat_source, simtime.index)?;

        Ok(energy_released_from_emitters)
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

    fn calc_emitter_cooldown(
        &self,
        _energy_demand: f64,
        _temp_emitter_req: f64,
        _temp_rm_prev: f64,
        _timestep: f64,
    ) -> (f64, f64) {
        unimplemented!()
    }

    fn temp_emitter(
        &self,
        _time_start: f64,
        _time_end: f64,
        _temp_emitter_start: f64,
        _temp_rm: f64,
        _power_input: f64,
        _temp_emitter_max: Option<f64>,
    ) -> (f64, Option<f64>) {
        unimplemented!()
    }

    fn energy_required_from_heat_source(
        &self,
        _energy_demand_heating_period: f64,
        _time_heating_start: f64,
        _timestep: f64,
        _temp_rm_prev: f64,
        _temp_emitter_heating_start: f64,
        _temp_emitter_req: f64,
        _temp_emitter_max: f64,
    ) -> (f64, bool) {
        unimplemented!()
    }

    fn energy_surplus_during_cooldown(
        &self,
        _time_cooldown: f64,
        _energy_demand: f64,
        _temp_rm_prev: f64,
        ThermalMassFields {
            thermal_mass: _thermal_mass,
            ..
        }: ThermalMassFields,
    ) -> f64 {
        unimplemented!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::controls::time_control::SetpointTimeControl;
    use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyBuilder};
    use crate::input::FuelType;
    use crate::simulation_time::SimulationTime;
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
    fn instant_elec_heater(simulation_time: SimulationTime) -> InstantElecHeater {
        let control = Control::SetpointTime(
            SetpointTimeControl::new(
                vec![Some(21.0), Some(21.0), None, Some(21.0)],
                0,
                1.,
                None,
                None,
                simulation_time.step,
            )
            .into(),
        );
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.iter().total_steps())
                .build(),
        ));
        let energy_supply_conn = EnergySupply::connection(energy_supply, "shower").unwrap();
        let temp_internal_air_fn = create_temp_internal_air_fn(20.);
        InstantElecHeater::new(
            50.,
            0.4,
            energy_supply_conn,
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
    #[ignore = "while migrating to 1.0.0a9"]
    fn test_demand_energy(instant_elec_heater: InstantElecHeater, simulation_time: SimulationTime) {
        let energy_input = [40.0, 100.0, 30.0, 20.0];
        let demand_expected = [40.0, 50.0, 0.0, 20.0];
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
}
