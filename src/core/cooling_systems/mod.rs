use crate::hem_core::simulation_time::SimulationTimeIteration;

pub mod air_conditioning;

pub(crate) trait SpaceCoolSystem {
    /// Return the temperature setpoint for the air conditioning system
    fn temp_setpnt(&self, simulation_time_iteration: &SimulationTimeIteration) -> Option<f64>;

    /// Return true if current time is inside specified time for heating/cooling
    fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool>;

    /// Return the convective fraction for cooling
    fn frac_convective(&self) -> f64;

    /// Return the minimum energy output of the air conditioning system
    fn energy_output_min(&self) -> f64;

    /// Demand energy (in kWh) from the cooling system
    fn demand_energy(&self, cooling_demand: f64, simtime: SimulationTimeIteration) -> f64;
}
