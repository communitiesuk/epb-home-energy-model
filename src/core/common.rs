// location for defining common traits and enums defined across submodules

use crate::core::heating_systems::storage_tank::HotWaterStorageTank;
use crate::core::heating_systems::wwhrs::Wwhrs;
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::simulation_time::SimulationTimeIteration;
use parking_lot::Mutex;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub(crate) enum WaterSupply {
    ColdWaterSource(Arc<ColdWaterSource>),
    Wwhrs(Arc<Mutex<Wwhrs>>),
    Preheated(HotWaterStorageTank),
    #[cfg(test)]
    Mock(MockWaterSupply),
}

impl WaterSupply {
    pub(crate) fn get_temp_cold_water(
        &self,
        volume_needed: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        match self {
            WaterSupply::ColdWaterSource(cold_water_source) => {
                cold_water_source.get_temp_cold_water(volume_needed, simtime)
            }
            WaterSupply::Preheated(storage_tank) => match storage_tank {
                HotWaterStorageTank::StorageTank(rw_lock) => {
                    Ok(rw_lock.read().get_temp_cold_water(volume_needed))
                }
                HotWaterStorageTank::SmartHotWaterTank(rw_lock) => {
                    Ok(rw_lock.read().get_temp_cold_water(volume_needed))
                }
            },
            #[cfg(test)]
            WaterSupply::Mock(mock) => Ok(mock.get_temp_cold_water(volume_needed)),
            _ => unimplemented!("TODO during migration 1.0.0a1"),
        }
    }

    pub(crate) fn draw_off_water(
        &self,
        volume_needed: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        match self {
            WaterSupply::ColdWaterSource(cold_water_source) => {
                cold_water_source.draw_off_water(volume_needed, simtime)
            }
            WaterSupply::Preheated(storage_tank) => match storage_tank {
                HotWaterStorageTank::StorageTank(rw_lock) => {
                    rw_lock.read().draw_off_water(volume_needed, simtime)
                }
                HotWaterStorageTank::SmartHotWaterTank(rw_lock) => {
                    rw_lock.read().draw_off_water(volume_needed, simtime)
                }
            },
            #[cfg(test)]
            WaterSupply::Mock(mock) => Ok(mock.draw_off_water(volume_needed)),
            _ => unimplemented!("TODO during migration 1.0.0a1"),
        }
    }
}

#[cfg(test)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct MockWaterSupply {
    temperature: f64,
}

#[cfg(test)]
impl MockWaterSupply {
    pub(crate) fn new(temperature: f64) -> Self {
        Self { temperature }
    }

    pub(crate) fn get_temp_cold_water(&self, volume_needed: f64) -> Vec<(f64, f64)> {
        vec![(self.temperature, volume_needed)]
    }

    pub(crate) fn draw_off_water(&self, volume_needed: f64) -> Vec<(f64, f64)> {
        self.get_temp_cold_water(volume_needed)
    }
}
