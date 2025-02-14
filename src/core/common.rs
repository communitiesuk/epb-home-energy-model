// location for defining common traits and enums defined across submodules

use crate::core::heating_systems::storage_tank::StorageTank;
use crate::core::heating_systems::wwhrs::Wwhrs;
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::simulation_time::SimulationTimeIteration;
use anyhow::bail;
use parking_lot::{Mutex, RwLock};
use std::sync::Arc;

#[derive(Clone, Debug)]
pub enum WaterSourceWithTemperature {
    ColdWaterSource(Arc<ColdWaterSource>),
    Wwhrs(Arc<Mutex<Wwhrs>>),
    Preheated(Arc<RwLock<StorageTank>>),
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
            WaterSourceWithTemperature::Preheated(storage_tank) => {
                storage_tank.write().temperature(volume_needed, simtime)
            }
        }
    }

    pub(crate) fn as_cold_water_source(&self) -> anyhow::Result<Arc<ColdWaterSource>> {
        Ok(match self {
            WaterSourceWithTemperature::ColdWaterSource(cold_water_source) => {
                cold_water_source.clone()
            }
            WaterSourceWithTemperature::Wwhrs(_) => {
                bail!("Water source is not a cold water source when it was expected to be.")
            }
            WaterSourceWithTemperature::Preheated(_) => {
                bail!("Storage tank is not a cold water source when it was expected to be.")
            }
        })
    }
}
