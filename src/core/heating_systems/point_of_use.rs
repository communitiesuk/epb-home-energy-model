use crate::compare_floats::min_of_2;
use crate::core::common::WaterSourceWithTemperature;
use crate::core::water_heat_demand::misc::water_demand_to_kwh;
use crate::simulation_time::SimulationTimeIteration;

pub struct PointOfUse {
    power_in_kw: f64,
    efficiency: f64,
    // energy_supply
    cold_feed: WaterSourceWithTemperature,
}

impl PointOfUse {
    pub fn new(
        rated_power_in_kw: f64,
        efficiency: f64,
        cold_feed: WaterSourceWithTemperature,
    ) -> Self {
        Self {
            power_in_kw: rated_power_in_kw,
            efficiency,
            cold_feed,
        }
    }

    pub fn get_cold_water_source(&self) -> &WaterSourceWithTemperature {
        &self.cold_feed
    }

    pub fn demand_hot_water(
        &self,
        volume_demanded: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        let demand_temp = 52.;

        let water_energy_demand = water_demand_to_kwh(
            volume_demanded,
            demand_temp,
            self.cold_feed.temperature(simulation_time_iteration.index),
        );

        // Assumption is that system specified has sufficient capacity to meet any realistic demand
        self.demand_energy(water_energy_demand, simulation_time_iteration)
    }

    /// Demand energy (in kWh) from the heater
    pub fn demand_energy(
        &self,
        energy_demand: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        let fuel_demand = min_of_2(
            energy_demand,
            self.power_in_kw * simulation_time_iteration.timestep * (1. / self.efficiency),
        );

        // TODO call on energy supply

        fuel_demand
    }
}
