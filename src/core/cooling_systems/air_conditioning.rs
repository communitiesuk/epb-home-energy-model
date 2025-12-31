use super::space_cool_system_base::SpaceCoolSystem;
use crate::compare_floats::max_of_2;
use crate::core::controls::time_control::{per_control, Control, ControlBehaviour};
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::simulation_time::SimulationTimeIteration;
use std::sync::Arc;

/// This module provides objects to air conditioning.

#[derive(Debug)]
pub struct AirConditioning {
    cooling_capacity_in_kw: f64,
    efficiency: f64,
    frac_convective: f64,
    energy_supply_connection: EnergySupplyConnection,
    simulation_timestep: f64,
    control: Arc<Control>,
}

impl AirConditioning {
    /// Construct an air conditioning object
    ///
    /// Arguments:
    /// * `cooling_capacity` - maximum cooling capacity of the system, in kW
    /// * `efficiency` - SEER
    /// * `frac_convective` - convective fraction for cooling
    /// * `energy_supply_connection` = an EnergySupplyConnection value
    /// * `simulation_timestep` - reference to timestep for contextual SimulationTime
    /// * `control` - reference to a control object
    pub(crate) fn new(
        cooling_capacity_in_kw: f64,
        efficiency: f64,
        frac_convective: f64,
        energy_supply_connection: EnergySupplyConnection,
        simulation_timestep: f64,
        control: Arc<Control>,
    ) -> Self {
        Self {
            cooling_capacity_in_kw,
            efficiency,
            frac_convective,
            energy_supply_connection,
            simulation_timestep,
            control,
        }
    }
}

impl SpaceCoolSystem for AirConditioning {
    fn temp_setpnt(&self, simulation_time_iteration: &SimulationTimeIteration) -> Option<f64> {
        per_control!(self.control.as_ref(), ctrl => { ctrl.setpnt(simulation_time_iteration) })
    }

    fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool> {
        per_control!(self.control.as_ref(), ctrl => { ctrl.in_required_period(simulation_time_iteration) })
    }

    fn frac_convective(&self) -> f64 {
        self.frac_convective
    }

    fn energy_output_min(&self) -> f64 {
        0.0
    }

    fn demand_energy(&self, cooling_demand: f64, simtime: SimulationTimeIteration) -> f64 {
        // Account for time control where present. If no control present, assume
        // system is always active (except for basic thermostatic control, which
        // is implicit in demand calculation).
        let cooling_supplied = if self.control.is_on(simtime) {
            max_of_2(
                cooling_demand,
                -self.cooling_capacity_in_kw * self.simulation_timestep,
            )
        } else {
            0.
        };

        self.energy_supply_connection
            .demand_energy(-cooling_supplied / self.efficiency, simtime.index)
            .unwrap();

        cooling_supplied
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
    pub fn aircon(simulation_time: SimulationTime) -> (AirConditioning, Arc<RwLock<EnergySupply>>) {
        let control = SetpointTimeControl::new(
            vec![Some(21.0), Some(21.0), None, Some(21.0)],
            0,
            1.0,
            Default::default(),
            Default::default(),
            1.0,
        );
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
        let energy_supply_conn = EnergySupply::connection(energy_supply.clone(), "aircon").unwrap();
        (
            AirConditioning::new(
                50.,
                2.0,
                0.4,
                energy_supply_conn,
                simulation_time.step,
                Arc::new(Control::SetpointTime(control)),
            ),
            energy_supply,
        )
    }

    #[rstest]
    fn test_demand_energy(
        aircon: (AirConditioning, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
    ) {
        let (aircon, energy_supply) = aircon;
        let inputs = [-40.0, -100.0, -30.0, -20.0];
        let expected_demand = [-40.0, -50.0, 0.0, -20.0];
        let expected_energy_supply_results = [20.0, 25.0, 0.0, 10.0];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                aircon.demand_energy(inputs[t_idx], t_it),
                expected_demand[t_idx],
                "incorrect cooling energy supplied returned"
            );
            assert_eq!(
                energy_supply.read().results_by_end_user()["aircon"][t_idx],
                expected_energy_supply_results[t_idx],
                "incorrect delivered energy demand returned"
            );
        }
    }

    #[rstest]
    fn test_energy_output_min(aircon: (AirConditioning, Arc<RwLock<EnergySupply>>)) {
        let (aircon, _) = aircon;
        assert_eq!(aircon.energy_output_min(), 0.0,);
    }

    #[rstest]
    fn test_temp_setpnt(
        aircon: (AirConditioning, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
    ) {
        let (aircon, _) = aircon;
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                aircon.temp_setpnt(&t_it),
                [Some(21.0), Some(21.0), None, Some(21.0)][t_idx]
            );
        }
    }

    #[rstest]
    fn test_in_required_period(
        aircon: (AirConditioning, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
    ) {
        let (aircon, _) = aircon;
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                aircon.in_required_period(&t_it).unwrap(),
                [true, true, false, true][t_idx]
            );
        }
    }

    #[rstest]
    fn test_frac_convective(aircon: (AirConditioning, Arc<RwLock<EnergySupply>>)) {
        let (aircon, _) = aircon;
        assert_eq!(aircon.frac_convective(), 0.4);
    }
}
