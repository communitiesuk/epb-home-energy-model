// location for defining common traits and enums defined across sub-modules

use crate::core::heating_systems::wwhrs::{WWHRSInstantaneousSystemA, WWHRSInstantaneousSystemC};
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub enum WaterSourceWithTemperature {
    ColdWaterSource(Arc<ColdWaterSource>),
    WwhrsC(Arc<WWHRSInstantaneousSystemC>),
    WwhrsA(Arc<WWHRSInstantaneousSystemA>),
}

impl WaterSourceWithTemperature {
    pub fn temperature(&self, timestep_idx: usize) -> f64 {
        match self {
            WaterSourceWithTemperature::ColdWaterSource(cold_water_source) => {
                cold_water_source.temperature(timestep_idx)
            }
            WaterSourceWithTemperature::WwhrsC(c) => c.temperature(),
            WaterSourceWithTemperature::WwhrsA(a) => a.temperature(),
        }
    }
}
