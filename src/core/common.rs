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
                    storage_tank.write().temperature(volume_needed, simtime)
                }
                HotWaterStorageTank::SmartHotWaterTank(_) => {
                    todo!("migration of storage tank module to 0.34")
                }
            },
        }
    }
}
