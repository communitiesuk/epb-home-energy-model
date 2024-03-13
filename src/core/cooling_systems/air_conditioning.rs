use crate::compare_floats::max_of_2;
use crate::core::controls::time_control::{per_control, Control, ControlBehaviour};
use crate::simulation_time::SimulationTimeIteration;
use std::sync::Arc;

/// This module provides objects to air conditioning.

#[derive(Debug)]
pub struct AirConditioning {
    cooling_capacity_in_kw: f64,
    efficiency: f64,
    frac_convective: f64,
    // energy supply
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
    /// * `simulation_timestep` - reference to timestep for contextual SimulationTime
    /// * `control` - reference to a control object
    pub fn new(
        cooling_capacity_in_kw: f64,
        efficiency: f64,
        frac_convective: f64,
        simulation_timestep: f64,
        control: Option<Arc<Control>>,
    ) -> Self {
        Self {
            cooling_capacity_in_kw,
            efficiency,
            frac_convective,
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

    pub fn demand_energy(&self, cooling_demand: f64, timestep_idx: usize) -> f64 {
        // Account for time control where present. If no control present, assume
        // system is always active (except for basic thermostatic control, which
        // is implicit in demand calculation).
        let cooling_supplied =
            if self.control.is_none() || self.control.as_ref().unwrap().is_on(timestep_idx) {
                max_of_2(
                    cooling_demand,
                    -self.cooling_capacity_in_kw * self.simulation_timestep,
                )
            } else {
                0.
            };

        // TODO energy supply reporting

        cooling_supplied
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::controls::time_control::OnOffTimeControl;
    use crate::simulation_time::SimulationTime;
    use rstest::*;

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 4., 1.)
    }

    #[fixture]
    pub fn aircon(simulation_time: SimulationTime) -> AirConditioning {
        let control = OnOffTimeControl::new(vec![true, true, false, true], 0, 1.0);
        AirConditioning::new(
            50.,
            2.0,
            0.4,
            simulation_time.step,
            Some(Arc::new(Control::OnOffTimeControl(control))),
        )
    }

    #[rstest]
    pub fn test_demand_energy(aircon: AirConditioning, simulation_time: SimulationTime) {
        let inputs = [-40.0, -100.0, -30.0, -20.0];
        let expected_demand = [-40.0, -50.0, 0.0, -20.0];

        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_eq!(
                aircon.demand_energy(inputs[t_idx], t_idx),
                expected_demand[t_idx],
                "incorrect cooling energy supplied returned"
            );
            // TODO test for energy supply
        }
    }
}
