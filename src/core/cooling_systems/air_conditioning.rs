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
    pub(crate) fn new(
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

    pub(crate) fn energy_output_min(&self) -> f64 {
        0.0
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
