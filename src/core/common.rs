// location for defining common traits and enums defined across submodules

use crate::core::heating_systems::storage_tank::{self, HotWaterStorageTank};
use crate::core::heating_systems::wwhrs::Wwhrs;
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::simulation_time::SimulationTimeIteration;
use parking_lot::Mutex;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub(crate) enum WaterSourceWithTemperature {
    ColdWaterSource(Arc<ColdWaterSource>),
    Wwhrs(Arc<Mutex<Wwhrs>>),
    Preheated(HotWaterStorageTank),
}

impl WaterSourceWithTemperature {
    pub(crate) fn temperature(
        &self,
        simtime: SimulationTimeIteration,
        volume_needed: Option<f64>,
    ) -> f64 {
        match self {
            WaterSourceWithTemperature::ColdWaterSource(cold_water_source) => {
                cold_water_source.temperature(simtime)
            }
            WaterSourceWithTemperature::Wwhrs(w) => w.lock().temperature(),
            WaterSourceWithTemperature::Preheated(source) => match source {
                HotWaterStorageTank::StorageTank(storage_tank) => {
                    todo!("temperature method no longer exists on storage tank")
                }
                HotWaterStorageTank::SmartHotWaterTank(smart_hot_water_tank) => {
                    todo!("temperature method no longer exists on storage tank")
                }
            },
        }
    }

    pub(crate) fn get_temp_cold_water(
        &self,
        volume_needed: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        match self {
            WaterSourceWithTemperature::ColdWaterSource(cold_water_source) => {
                cold_water_source.get_temp_cold_water(volume_needed, simtime)
            }
            WaterSourceWithTemperature::Preheated(storage_tank) => match storage_tank {
                HotWaterStorageTank::StorageTank(rw_lock) => {
                    Ok(rw_lock.read().get_temp_cold_water(volume_needed))
                }
                HotWaterStorageTank::SmartHotWaterTank(rw_lock) => {
                    Ok(rw_lock.read().get_temp_cold_water(volume_needed))
                }
            },
            _ => unimplemented!("TODO during migration 1.0.0a1"),
        }
    }

    pub(crate) fn draw_off_water(
        &self,
        volume_needed: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        match self {
            WaterSourceWithTemperature::ColdWaterSource(cold_water_source) => {
                cold_water_source.draw_off_water(volume_needed, simtime)
            }
            WaterSourceWithTemperature::Preheated(storage_tank) => match storage_tank {
                HotWaterStorageTank::StorageTank(rw_lock) => {
                    rw_lock.read().draw_off_water(volume_needed, simtime)
                }
                HotWaterStorageTank::SmartHotWaterTank(rw_lock) => {
                    rw_lock.read().draw_off_water(volume_needed, simtime)
                }
            },
            _ => unimplemented!("TODO during migration 1.0.0a1"),
        }
    }
}
