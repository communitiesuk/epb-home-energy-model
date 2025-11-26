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
                    storage_tank.read().temperature(volume_needed, simtime)
                }
                HotWaterStorageTank::SmartHotWaterTank(smart_hot_water_tank) => {
                    smart_hot_water_tank
                        .read()
                        .temperature(volume_needed, simtime)
                }
            },
        }
    }
    
    // NOTE added for use by storage tank
    // TODO do we need simulation time?
    pub(crate) fn get_temp_cold_water(&self, _volume_needed: f64) -> Vec<(f64, f64)> {
        todo!()
    }

    // NOTE added for use by storage tank
    pub(crate) fn draw_off_water(&self, volume_needed: f64) -> Vec<(f64, f64)> {
        // this seems to just be a proxy for the above method
        self.get_temp_cold_water(volume_needed)
    }
}
