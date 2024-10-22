// location for defining common traits and enums defined across submodules

use crate::core::heating_systems::wwhrs::Wwhrs;
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::simulation_time::SimulationTimeIteration;
use anyhow::bail;
use parking_lot::Mutex;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub enum WaterSourceWithTemperature {
    ColdWaterSource(Arc<ColdWaterSource>),
    Wwhrs(Arc<Mutex<Wwhrs>>),
}

impl WaterSourceWithTemperature {
    pub fn temperature(&self, simtime: SimulationTimeIteration) -> f64 {
        match self {
            WaterSourceWithTemperature::ColdWaterSource(cold_water_source) => {
                cold_water_source.temperature(simtime)
            }
            WaterSourceWithTemperature::Wwhrs(w) => w.lock().temperature(),
        }
    }

    pub fn as_cold_water_source(&self) -> anyhow::Result<Arc<ColdWaterSource>> {
        Ok(match self {
            WaterSourceWithTemperature::ColdWaterSource(cold_water_source) => {
                cold_water_source.clone()
            }
            WaterSourceWithTemperature::Wwhrs(_) => {
                bail!("Water source is not a cold water source when it was expected to be.")
            }
        })
    }
}
