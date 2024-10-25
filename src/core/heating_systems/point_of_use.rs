use crate::core::common::WaterSourceWithTemperature;
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::water_heat_demand::dhw_demand::{DemandVolTargetKey, VolumeReference};
use crate::core::water_heat_demand::misc::water_demand_to_kwh;
use crate::simulation_time::SimulationTimeIteration;
use indexmap::IndexMap;

#[derive(Debug)]
pub struct PointOfUse {
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
            self.cold_feed.temperature(*simulation_time_iteration),
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        core::{
            energy_supply::energy_supply::EnergySupply,
            water_heat_demand::cold_water_source::ColdWaterSource,
        },
        input::FuelType,
        simulation_time::{SimulationTime, SimulationTimeIterator},
    };
    use approx::assert_relative_eq;
    use parking_lot::lock_api::RwLock;
    use rstest::{fixture, rstest};
    use std::sync::Arc;

    #[fixture]
    fn simtime() -> SimulationTime {
        SimulationTime::new(0.0, 2.0, 1.0)
    }

    #[fixture]
    fn simulation_time_iterator(simtime: SimulationTime) -> SimulationTimeIterator {
        simtime.iter()
    }

    #[fixture]
    pub fn energy_supply(simtime: SimulationTime) -> EnergySupply {
        EnergySupply::new(
            FuelType::Electricity,
            simtime.total_steps(),
            None,
            None,
            None,
        )
    }

    #[fixture]
    pub fn point_of_use(energy_supply: EnergySupply, simtime: SimulationTime) -> PointOfUse {
        let efficiency = 1.;
        let energy_supply = Arc::new(RwLock::new(energy_supply));
        let energy_supply_connection =
            EnergySupply::connection(energy_supply.clone(), "electricity").unwrap();
        let cold_water_temps = vec![15., 20., 25.];
        let coldfeed = WaterSourceWithTemperature::ColdWaterSource(Arc::new(ColdWaterSource::new(
            cold_water_temps,
            &simtime,
            0,
            1.,
        )));
        let temp_hot_water = 55.;

        PointOfUse::new(
            efficiency,
            energy_supply_connection,
            coldfeed,
            temp_hot_water,
        )
    }

    #[rstest]
    fn test_demand_hot_water(
        point_of_use: PointOfUse,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        // Test when temp_hot_water is set
        let volume_demanded_target: IndexMap<DemandVolTargetKey, VolumeReference> =
            IndexMap::from([(
                DemandVolTargetKey::TempHotWater,
                VolumeReference {
                    warm_temp: 0.0, // warm_temp not used in this test
                    warm_vol: 60.0,
                },
            )]);

        assert_relative_eq!(
            point_of_use.demand_hot_water(
                volume_demanded_target,
                &simulation_time_iterator.current_iteration()
            ),
            2.7893333333333334
        );

        // Test when temp_hot_water is not set
        let volume_demanded_target: IndexMap<DemandVolTargetKey, VolumeReference> =
            IndexMap::from([]);

        assert_relative_eq!(
            point_of_use.demand_hot_water(
                volume_demanded_target,
                &simulation_time_iterator.current_iteration()
            ),
            0.
        );
    }

    #[rstest]
    fn test_demand_energy(
        point_of_use: PointOfUse,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let energy_demand = 2.0;
        assert_relative_eq!(
            point_of_use
                .demand_energy(energy_demand, &simulation_time_iterator.current_iteration()),
            2.0
        )
    }
}
