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
    control: Option<Arc<Control>>,
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
    pub fn new(
        cooling_capacity_in_kw: f64,
        efficiency: f64,
        frac_convective: f64,
        energy_supply_connection: EnergySupplyConnection,
        simulation_timestep: f64,
        control: Option<Arc<Control>>,
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

    pub fn temp_setpnt(&self, simulation_time_iteration: &SimulationTimeIteration) -> Option<f64> {
        self.control.as_ref().and_then(|control| per_control!(control.as_ref(), ctrl => { ctrl.setpnt(simulation_time_iteration) }))
    }

    pub fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool> {
        self.control.as_ref().and_then(|control| per_control!(control.as_ref(), ctrl => { ctrl.in_required_period(simulation_time_iteration) }))
    }

    pub fn frac_convective(&self) -> f64 {
        self.frac_convective
    }

    pub fn demand_energy(&self, cooling_demand: f64, simtime: SimulationTimeIteration) -> f64 {
        // Account for time control where present. If no control present, assume
        // system is always active (except for basic thermostatic control, which
        // is implicit in demand calculation).
        let cooling_supplied =
            if self.control.is_none() || self.control.as_ref().unwrap().is_on(simtime) {
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
    use crate::core::controls::time_control::OnOffTimeControl;
    use crate::core::energy_supply::energy_supply::EnergySupply;
    use crate::input::FuelType;
    use crate::simulation_time::SimulationTime;
    use parking_lot::Mutex;
    use rstest::*;

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 4., 1.)
    }

    #[fixture]
    pub fn aircon(simulation_time: SimulationTime) -> (AirConditioning, Arc<Mutex<EnergySupply>>) {
        let control = OnOffTimeControl::new(vec![true, true, false, true], 0, 1.0);
        let energy_supply = Arc::new(Mutex::new(EnergySupply::new(
            FuelType::Electricity,
            simulation_time.total_steps(),
            None,
        )));
        let energy_supply_conn = EnergySupply::connection(energy_supply.clone(), "aircon").unwrap();
        (
            AirConditioning::new(
                50.,
                2.0,
                0.4,
                energy_supply_conn,
                simulation_time.step,
                Some(Arc::new(Control::OnOffTimeControl(control))),
            ),
            energy_supply,
        )
    }

    #[rstest]
    pub fn test_demand_energy(
        aircon: (AirConditioning, Arc<Mutex<EnergySupply>>),
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
                energy_supply.lock().results_by_end_user()["aircon"][t_idx],
                expected_energy_supply_results[t_idx],
                "incorrect delivered energy demand returned"
            );
        }
    }
}
