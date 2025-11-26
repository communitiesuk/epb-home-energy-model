// location for defining common traits and enums defined across submodules

use crate::core::heating_systems::storage_tank::HotWaterStorageTank;
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

    // NOTE this will be replaced by get_temp_cold_water and draw_off_water
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
                    todo!() // switched to get_temp_cold_water and draw_off_water in latest
                }
                HotWaterStorageTank::SmartHotWaterTank(smart_hot_water_tank) => {
                    smart_hot_water_tank
                        .read()
                        .temperature(volume_needed, simtime)
                }
            },
        }
    }
    
    // NOTE this may need simulation_time in future
    pub(crate) fn get_temp_cold_water(&self, volume_needed: f64, _simulation_time: SimulationTimeIteration) -> Vec<(f64, f64)> {
        match self {
            WaterSourceWithTemperature::ColdWaterSource(_cold_water_source) => todo!(),
            WaterSourceWithTemperature::Wwhrs(_mutex) => todo!(),
            WaterSourceWithTemperature::Preheated(hot_water_storage_tank) => match hot_water_storage_tank {
                HotWaterStorageTank::StorageTank(rw_lock) => rw_lock.read().get_temp_cold_water(volume_needed),
                HotWaterStorageTank::SmartHotWaterTank(_rw_lock) => todo!(),
            },
        }
    }

    pub(crate) fn draw_off_water(&self, volume_needed: f64, simulation_time_iteration: SimulationTimeIteration) -> Vec<(f64, f64)> {
        match self {
            WaterSourceWithTemperature::ColdWaterSource(_cold_water_source) => todo!(),
            WaterSourceWithTemperature::Wwhrs(_mutex) => todo!(),
            WaterSourceWithTemperature::Preheated(hot_water_storage_tank) => match hot_water_storage_tank {
                HotWaterStorageTank::StorageTank(rw_lock) => rw_lock.read().draw_off_water(volume_needed, simulation_time_iteration),
                HotWaterStorageTank::SmartHotWaterTank(_rw_lock) => todo!(),
            },
        }
    }
}
