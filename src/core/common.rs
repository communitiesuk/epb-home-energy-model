// location for defining common traits and enums defined across submodules

use crate::core::heating_systems::wwhrs::Wwhrs;
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use parking_lot::Mutex;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub enum WaterSourceWithTemperature {
    ColdWaterSource(Arc<ColdWaterSource>),
    Wwhrs(Arc<Mutex<Wwhrs>>),
}

impl WaterSourceWithTemperature {
    pub fn temperature(&self, timestep_idx: usize) -> f64 {
        match self {
            WaterSourceWithTemperature::ColdWaterSource(cold_water_source) => {
                cold_water_source.temperature(timestep_idx)
            }
            WaterSourceWithTemperature::Wwhrs(w) => w.lock().temperature(),
        }
    }
}
