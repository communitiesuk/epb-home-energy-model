/// This module provides object(s) to model the behaviour of instantaneous electric
/// room heaters.
use crate::compare_floats::min_of_2;
use crate::core::controls::time_control::{per_control, Control, ControlBehaviour};
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::simulation_time::SimulationTimeIteration;
use std::sync::Arc;

/// Type to represent instantaneous electric heaters
#[derive(Clone, Debug)]
pub struct InstantElecHeater {
    rated_power_in_kw: f64,
    frac_convective: f64,
    energy_supply_connection: EnergySupplyConnection,
    simulation_timestep: f64,
    control: Option<Arc<Control>>,
}

impl InstantElecHeater {
    /// Arguments
    /// * `rated_power` - in kW
    /// * `frac_convective` - convective fraction for heating
    /// * `energy_supply_connection` - EnergySupplyConnection value
    /// * `simulation_timestep` - step in hours for context simulation time
    /// * `control` - reference to a control object which must implement is_on() and setpnt() funcs
    pub(crate) fn new(
        rated_power_in_kw: f64,
        frac_convective: f64,
        energy_supply_connection: EnergySupplyConnection,
        simulation_timestep: f64,
        control: Option<Arc<Control>>,
    ) -> Self {
        Self {
            rated_power_in_kw,
            frac_convective,
            energy_supply_connection,
            simulation_timestep,
            control,
        }
    }

    pub fn temp_setpnt(&self, simulation_time_iteration: &SimulationTimeIteration) -> Option<f64> {
        self.control.as_ref().and_then(
            |ctrl| per_control!(ctrl.as_ref(), ctrl => { ctrl.setpnt(simulation_time_iteration) }),
        )
    }

    pub fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool> {
        self.control.as_ref().and_then(|ctrl| per_control!(ctrl.as_ref(), ctrl => { ctrl.in_required_period(simulation_time_iteration) }))
    }

    pub fn frac_convective(&self) -> f64 {
        self.frac_convective
    }

    pub(crate) fn energy_output_min(&self) -> f64 {
        0.
    }

    /// Demand energy (in kWh) from the heater
    pub fn demand_energy(&self, energy_demand: f64, simtime: SimulationTimeIteration) -> f64 {
        // Account for time control where present. If no control present, assume
        // system is always active (except for basic thermostatic control, which
        // is implicit in demand calculation).
        let energy_supplied =
            if self.control.as_ref().is_none() || self.control.as_ref().unwrap().is_on(simtime) {
                min_of_2(
                    energy_demand,
                    self.rated_power_in_kw * self.simulation_timestep,
                )
            } else {
                0.
            };

        self.energy_supply_connection
            .demand_energy(energy_supplied, simtime.index)
            .unwrap();

        energy_supplied
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

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 4., 1.)
    }

    #[fixture]
    pub fn instant_elec_heater(simulation_time: SimulationTime) -> InstantElecHeater {
        let control = Control::SetpointTime(
            SetpointTimeControl::new(
                vec![Some(21.0), Some(21.0), None, Some(21.0)],
                0,
                1.,
                None,
                None,
                None,
                Default::default(),
                simulation_time.step,
            )
            .unwrap(),
        );
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
        let energy_supply_conn = EnergySupply::connection(energy_supply, "shower").unwrap();
        InstantElecHeater::new(
            50.,
            0.4,
            energy_supply_conn,
            simulation_time.step,
            Some(Arc::new(control)),
        )
    }

    #[rstest]
    fn test_demand_energy(instant_elec_heater: InstantElecHeater, simulation_time: SimulationTime) {
        let energy_input = [40.0, 100.0, 30.0, 20.0];
        let demand_expected = [40.0, 50.0, 0.0, 20.0];
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                instant_elec_heater.demand_energy(energy_input[t_idx], t_it),
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
