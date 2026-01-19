use crate::core::common::WaterSupply;
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::water_heat_demand::misc::{water_demand_to_kwh, WaterEventResult};
use crate::simulation_time::SimulationTimeIteration;

#[derive(Debug, Clone)]
pub(crate) struct PointOfUse {
    efficiency: f64,
    energy_supply_connection: EnergySupplyConnection,
    cold_feed: WaterSupply,
    temp_hot_water: f64,
}

impl PointOfUse {
    pub(crate) fn new(
        efficiency: f64,
        energy_supply_connection: EnergySupplyConnection,
        cold_feed: WaterSupply,
        temp_hot_water: f64,
    ) -> Self {
        Self {
            efficiency,
            energy_supply_connection,
            cold_feed,
            temp_hot_water,
        }
    }

    pub fn get_cold_water_source(&self) -> &WaterSupply {
        &self.cold_feed
    }

    pub fn get_temp_hot_water(
        &self,
        volume_req: f64,
        _volume_req_already: Option<f64>,
    ) -> Vec<(f64, f64)> {
        // Always supplies the whole volume at the same temperature, so list has a single element
        vec![(self.temp_hot_water, volume_req)]
    }

    pub fn demand_hot_water(
        &self,
        usage_events: Vec<WaterEventResult>,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let mut water_energy_demand = 0.;

        for event in usage_events {
            if is_close!(event.volume_hot, 0., rel_tol = 1e-09, abs_tol = 1e-10) {
                continue;
            }
            let list_temp_volume = self
                .cold_feed
                .draw_off_water(event.volume_hot, *simulation_time_iteration)?;
            let sum_t_by_v: f64 = list_temp_volume.iter().map(|(t, v)| t * v).sum();
            let sum_v: f64 = list_temp_volume.iter().map(|(_, v)| v).sum();

            let temp_cold_water = sum_t_by_v / sum_v;

            water_energy_demand +=
                water_demand_to_kwh(event.volume_hot, self.temp_hot_water, temp_cold_water);
        }
        // Assumption is that system specified has sufficient capacity to meet any realistic demand
        Ok(self.demand_energy(water_energy_demand, simulation_time_iteration))
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
    use crate::core::energy_supply::energy_supply::EnergySupplyBuilder;
    use crate::core::water_heat_demand::misc::WaterEventResultType;
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
    fn energy_supply(simtime: SimulationTime) -> EnergySupply {
        EnergySupplyBuilder::new(FuelType::Electricity, simtime.total_steps()).build()
    }

    #[fixture]
    fn point_of_use(energy_supply: EnergySupply) -> PointOfUse {
        let efficiency = 1.;
        let energy_supply = Arc::new(RwLock::new(energy_supply));
        let energy_supply_connection =
            EnergySupply::connection(energy_supply.clone(), "electricity").unwrap();
        let cold_water_temps = vec![15., 20., 25.];
        let coldfeed =
            WaterSupply::ColdWaterSource(Arc::new(ColdWaterSource::new(cold_water_temps, 0, 1.)));
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
        let usage_events = vec![
            WaterEventResult {
                event_result_type: WaterEventResultType::Other,
                volume_hot: 60.,
                volume_warm: 60.,
                temperature_warm: 55.,
            },
            WaterEventResult {
                event_result_type: WaterEventResultType::Other,
                volume_hot: 0.,
                volume_warm: 0.,
                temperature_warm: 55.,
            },
        ];

        assert_relative_eq!(
            point_of_use
                .demand_hot_water(usage_events, &simulation_time_iterator.current_iteration())
                .unwrap(),
            2.7893333333333334
        );

        // Test when events list is empty
        let usage_events = vec![];

        assert_relative_eq!(
            point_of_use
                .demand_hot_water(usage_events, &simulation_time_iterator.current_iteration())
                .unwrap(),
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

    #[rstest]
    fn test_get_cold_water_source(point_of_use: PointOfUse) {
        let actual = point_of_use.get_cold_water_source();
        let expected = &point_of_use.cold_feed;

        match (actual, expected) {
            (WaterSupply::ColdWaterSource(actual), WaterSupply::ColdWaterSource(expected)) => {
                assert_eq!(actual, expected);
            }
            _ => panic!("Expected ColdWaterSource variant"),
        }
    }

    #[rstest]
    fn test_get_temp_hot_water(point_of_use: PointOfUse) {
        assert_eq!(point_of_use.get_temp_hot_water(20., None), vec![(55., 20.)]);
    }
}
