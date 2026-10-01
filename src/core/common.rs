// location for defining common traits and enums defined across submodules

use crate::core::heating_systems::storage_tank::HotWaterStorageTank;
use crate::core::heating_systems::wwhrs::WwhrsInstantaneous;
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::simulation_time::SimulationTimeIteration;
use parking_lot::Mutex;
#[cfg(test)]
use parking_lot::RwLock;
use std::sync::Arc;

pub trait WaterSupplyBehaviour: Clone {
    fn get_temp_cold_water(
        &self,
        volume_needed: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>>;

    fn draw_off_water(
        &self,
        volume_needed: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>>;

    fn ultimate_cold_water_source(&self) -> Self {
        self.clone()
    }
}

#[derive(Clone, Debug)]
pub enum WaterSupply {
    ColdWaterSource(Arc<ColdWaterSource>),
    Wwhrs(Arc<Mutex<WwhrsInstantaneous>>),
    Preheated(HotWaterStorageTank),
    #[cfg(test)]
    Mock(MockWaterSupply),
    #[cfg(test)]
    VaryingTemp(VaryingTempWaterSupply),
}

impl WaterSupplyBehaviour for WaterSupply {
    fn get_temp_cold_water(
        &self,
        volume_needed: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        match self {
            WaterSupply::ColdWaterSource(cold_water_source) => {
                cold_water_source.get_temp_cold_water(volume_needed, simtime)
            }
            WaterSupply::Wwhrs(wwhrs) => wwhrs.lock().get_temp_cold_water(volume_needed, simtime),
            WaterSupply::Preheated(storage_tank) => match storage_tank {
                HotWaterStorageTank::StorageTank(rw_lock) => {
                    rw_lock.read().get_temp_cold_water(volume_needed, simtime)
                }
                HotWaterStorageTank::SmartHotWaterTank(rw_lock) => {
                    rw_lock.read().get_temp_cold_water(volume_needed, simtime)
                }
                #[cfg(test)]
                HotWaterStorageTank::Mock(_source) => Ok(vec![]),
            },
            #[cfg(test)]
            WaterSupply::Mock(mock) => mock.get_temp_cold_water(volume_needed, simtime),
            #[cfg(test)]
            WaterSupply::VaryingTemp(varying_temp) => {
                varying_temp.get_temp_cold_water(volume_needed, simtime)
            }
        }
    }

    fn draw_off_water(
        &self,
        volume_needed: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        match self {
            WaterSupply::ColdWaterSource(cold_water_source) => {
                cold_water_source.draw_off_water(volume_needed, simtime)
            }
            WaterSupply::Wwhrs(wwhrs) => wwhrs.lock().draw_off_water(volume_needed, simtime),
            WaterSupply::Preheated(storage_tank) => match storage_tank {
                HotWaterStorageTank::StorageTank(rw_lock) => {
                    rw_lock.read().draw_off_water(volume_needed, simtime)
                }
                HotWaterStorageTank::SmartHotWaterTank(rw_lock) => {
                    rw_lock.read().draw_off_water(volume_needed, simtime)
                }
                #[cfg(test)]
                HotWaterStorageTank::Mock(_source) => Ok(vec![]),
            },
            #[cfg(test)]
            WaterSupply::Mock(mock) => mock.draw_off_water(volume_needed, simtime),
            #[cfg(test)]
            WaterSupply::VaryingTemp(varying_temp) => {
                varying_temp.draw_off_water(volume_needed, simtime)
            }
        }
    }

    fn ultimate_cold_water_source(&self) -> Self {
        match self {
            WaterSupply::Preheated(tank) => tank.ultimate_cold_water_source(),
            _ => self.clone(),
        }
    }
}

#[cfg(test)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MockWaterSupply {
    temperature: f64,
}

#[cfg(test)]
impl MockWaterSupply {
    pub(crate) fn new(temperature: f64) -> Self {
        Self { temperature }
    }
}

#[cfg(test)]
impl WaterSupplyBehaviour for MockWaterSupply {
    fn get_temp_cold_water(
        &self,
        volume_needed: f64,
        _simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        Ok(vec![(self.temperature, volume_needed)])
    }

    fn draw_off_water(
        &self,
        volume_needed: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        self.get_temp_cold_water(volume_needed, simtime)
    }
}

#[cfg(test)]
impl Default for MockWaterSupply {
    fn default() -> Self {
        Self { temperature: 10. }
    }
}
#[cfg(test)]
#[derive(Default, Clone, Debug)]
pub struct VaryingTempWaterSupply {
    volumes_passed_to_draw_off_hot_water: Arc<RwLock<Vec<f64>>>,
}

#[cfg(test)]
impl VaryingTempWaterSupply {
    pub fn new(volumes_container: Arc<RwLock<Vec<f64>>>) -> Self {
        Self {
            volumes_passed_to_draw_off_hot_water: volumes_container,
        }
    }

    pub fn register_call_to_draw_off_water(&self, volume: f64) {
        self.volumes_passed_to_draw_off_hot_water
            .write()
            .push(volume);
    }

    pub fn volumes_passed_to_draw_off_water(&self) -> Vec<f64> {
        self.volumes_passed_to_draw_off_hot_water.read().clone()
    }
    pub fn clear_volumes(&self) {
        self.volumes_passed_to_draw_off_hot_water.write().clear();
    }
}

#[cfg(test)]
// Set up cold feed to return different temperatures based on volume
// Simulates drawing from a stratified tank or mixed sources
fn varying_temp_by_volume(volume_needed: f64) -> Vec<(f64, f64)> {
    let volume = volume_needed;

    if volume <= 10. {
        // Small volume - warm water from top of tank
        vec![(15.0, volume)]
    } else if volume <= 30. {
        // Medium volume - mix of warm and cold
        let warm_portion = 10.;
        let cold_portion = volume - 10.;
        vec![(15.0, warm_portion), (8.0, cold_portion)]
    } else {
        // Large volume - mostly cold water
        vec![(15.0, 10.), (8.0, 20.), (5.0, volume - 30.)]
    }
}
#[cfg(test)]
impl WaterSupplyBehaviour for VaryingTempWaterSupply {
    fn get_temp_cold_water(
        &self,
        volume_needed: f64,
        _simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        Ok(varying_temp_by_volume(volume_needed))
    }

    fn draw_off_water(
        &self,
        volume_needed: f64,
        _simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        self.register_call_to_draw_off_water(volume_needed);
        Ok(varying_temp_by_volume(volume_needed))
    }
}
