use indexmap::IndexMap;

use crate::core::common::WaterSourceWithTemperature;
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::water_heat_demand::dhw_demand::{DemandVolTargetKey, VolumeReference};
use crate::core::water_heat_demand::misc::water_demand_to_kwh;
use crate::simulation_time::SimulationTimeIteration;

pub struct PointOfUse {
    power_in_kw: f64,
    efficiency: f64,
    energy_supply_connection: EnergySupplyConnection,
    cold_feed: WaterSourceWithTemperature,
    temp_hot_water: f64,
}

impl PointOfUse {
    pub fn new(
        efficiency: f64,
        energy_supply_connection: EnergySupplyConnection,
        cold_feed: WaterSourceWithTemperature,
        temp_hot_water: f64,
    ) -> Self {
        Self {
            power_in_kw: 1., // NB. temporary value until this module is migrated to 0.30
            efficiency,
            energy_supply_connection,
            cold_feed,
            temp_hot_water,
        }
    }

    pub fn get_cold_water_source(&self) -> &WaterSourceWithTemperature {
        &self.cold_feed
    }

    pub fn get_temp_hot_water(&self) -> f64 {
        self.temp_hot_water
    }

    pub fn demand_hot_water(
        &self,
        volume_demanded_target: IndexMap<DemandVolTargetKey, VolumeReference>,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        let demand_temp = self.temp_hot_water;

        // (From Python) TODO set required temperature rather than hard coding - also elsewhere in the code
        let volume_demanded = volume_demanded_target
            .get(&DemandVolTargetKey::TempHotWater)
            .map(|volume_reference| volume_reference.warm_vol)
            .unwrap_or(0.0);

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
        // Energy that heater is able to supply is limited by power rating
        let fuel_demand = energy_demand * (1. / self.efficiency);

        self.energy_supply_connection
            .demand_energy(fuel_demand, simulation_time_iteration.index)
            .unwrap();

        fuel_demand
    }
}
