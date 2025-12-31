use crate::hem_core::simulation_time::SimulationTimeIteration;

pub(crate) trait OnSiteGeneration {
    /// Produce energy from the on-site generation system
    fn produce_energy(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)>;

    /// Return whether this unit is considered inside the building or not
    fn inverter_is_inside(&self) -> bool;
}
