use crate::core::controls::time_control::{per_control, Control, ControlBehaviour};
use crate::simulation_time::SimulationTimeIteration;
use std::sync::Arc;

/// This module provides object(s) to model the behaviour of instantaneous electric
/// room heaters.

/// Type to represent instantaneous electric heaters
#[derive(Clone)]
pub struct InstantElecHeater {
    rated_power_in_kw: f64,
    frac_convective: f64,
    simulation_timestep: f64,
    control: Option<Arc<Control>>,
}

impl InstantElecHeater {
    /// Arguments
    /// * `rated_power` - in kW
    /// * `frac_convective` - convective fraction for heating
    /// * `simulation_timestep` - step in hours for context simulation time
    /// * `control` - reference to a control object which must implement is_on() and setpnt() funcs
    pub fn new(
        rated_power_in_kw: f64,
        frac_convective: f64,
        simulation_timestep: f64,
        control: Option<Arc<Control>>,
    ) -> Self {
        Self {
            rated_power_in_kw,
            frac_convective,
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

    /// Demand energy (in kWh) from the heater
    pub fn demand_energy(&self, energy_demand: f64, timestep_idx: usize) -> f64 {
        // Account for time control where present. If no control present, assume
        // system is always active (except for basic thermostatic control, which
        // is implicit in demand calculation).
        let energy_supplied = if self.control.as_ref().is_none()
            || self.control.as_ref().unwrap().is_on(timestep_idx)
        {
            *[
                energy_demand,
                self.rated_power_in_kw * self.simulation_timestep,
            ]
            .iter()
            .max_by(|a, b| a.total_cmp(b).reverse())
            .unwrap()
        } else {
            0.
        };

        // TODO register energy supply demand

        energy_supplied
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
    pub fn instant_elec_heater(simulation_time: SimulationTime) -> InstantElecHeater {
        let control =
            Control::OnOffTimeControl(OnOffTimeControl::new(vec![true, true, false, true], 0, 1.));
        InstantElecHeater::new(50., 0.4, simulation_time.step, Some(Arc::new(control)))
    }

    #[rstest]
    pub fn should_calc_demand_energy(
        instant_elec_heater: InstantElecHeater,
        simulation_time: SimulationTime,
    ) {
        let energy_input = [40.0, 100.0, 30.0, 20.0];
        let demand_expected = [40.0, 50.0, 0.0, 20.0];
        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_eq!(
                instant_elec_heater.demand_energy(energy_input[t_idx], t_idx),
                demand_expected[t_idx]
            );
        }
    }
}
