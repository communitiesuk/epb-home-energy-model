use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::misc::{
    volume_hot_water_required, CallableGetHotWaterTemperature,
};
use crate::input::WaterHeatingEvent;
use crate::simulation_time::SimulationTimeIteration;
use anyhow::bail;
use std::sync::Arc;

#[derive(Debug)]
pub struct Bath {
    size_in_litres: f64,
    cold_water_source: Arc<ColdWaterSource>,
    flowrate: f64,
}

impl Bath {
    pub(crate) fn new(
        size_in_litres: f64,
        cold_water_source: Arc<ColdWaterSource>,
        flowrate: f64,
    ) -> Self {
        Self {
            size_in_litres,
            cold_water_source,
            flowrate,
        }
    }

    pub(crate) fn get_cold_water_source(&self) -> &ColdWaterSource {
        &self.cold_water_source
    }

    pub(crate) fn get_flowrate(&self) -> f64 {
        self.flowrate
    }

    /// Calculate volume of hot water required
    /// (and volume of warm water draining to WWHRS, if applicable)
    pub(crate) fn hot_water_demand(
        &self,
        event: WaterHeatingEvent,
        func_temp_hot_water: &CallableGetHotWaterTemperature,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(Option<f64>, f64)> {
        let peak_flowrate = self.flowrate;

        let (vol_warm_water, bath_duration) = match (event.volume, event.duration) {
            (Some(volume), _) => (volume, volume / peak_flowrate),
            (_, Some(duration)) => (duration * self.flowrate, duration),
            _ => bail!("Invalid bath event {event:?} - must specify either volume or duration"),
        };
        let temp_target = event.temperature;

        let vol_warm_water = vol_warm_water.min(self.size_in_litres);
        let vol_hot_water = volume_hot_water_required(
            vol_warm_water,
            temp_target,
            func_temp_hot_water,
            |x, simtime| self.cold_water_source.get_temp_cold_water(x, simtime),
            simtime,
        )?;
        if let Some(vol_hot_water) = vol_hot_water {
            let vol_cold_water = vol_warm_water - vol_hot_water;
            self.cold_water_source
                .draw_off_water(vol_cold_water, simtime)?;
        }

        Ok((vol_hot_water, vol_warm_water))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation_time::SimulationTime;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    fn bath() -> Bath {
        let cold_water_source = ColdWaterSource::new(vec![2.0, 3.0, 4.0], 0, 1.0);
        Bath::new(100.0, cold_water_source.into(), 4.5)
    }

    #[rstest]
    fn should_give_cold_water_source(bath: Bath) {
        assert_eq!(
            bath.get_cold_water_source(),
            bath.cold_water_source.as_ref(),
            "cold water source not returned"
        );
    }

    fn func_temp_hot_water_fixed(_t: f64) -> f64 {
        52.0
    }

    #[fixture]
    fn func_temp_hot_water() -> CallableGetHotWaterTemperature {
        Box::new(func_temp_hot_water_fixed)
    }

    #[rstest]
    fn test_hot_water_demand(bath: Bath, func_temp_hot_water: CallableGetHotWaterTemperature) {
        let simulation_time = SimulationTime::new(0.0, 3.0, 1.0);

        assert_eq!(
            bath.hot_water_demand(
                WaterHeatingEvent {
                    volume: Some(75.0),
                    temperature: 40.0,
                    start: 0.,
                    duration: None
                },
                &func_temp_hot_water,
                simulation_time.iter().current_iteration()
            )
            .unwrap(),
            (Some(57.0), 75.0),
            "incorrect hot water demand returned"
        );
        assert_eq!(
            bath.hot_water_demand(
                WaterHeatingEvent {
                    volume: Some(200.0),
                    temperature: 40.0,
                    start: 0.,
                    duration: None
                },
                &func_temp_hot_water,
                simulation_time.iter().current_iteration()
            )
            .unwrap(),
            (Some(76.0), 100.0),
            "incorrect hot water demand returned for bath fill volume > bath tub volume"
        );
        assert_eq!(
            bath.hot_water_demand(
                WaterHeatingEvent {
                    duration: Some(16.666666666666667),
                    temperature: 40.0,
                    start: 0.,
                    volume: None
                },
                &func_temp_hot_water,
                simulation_time.iter().current_iteration()
            )
            .unwrap(),
            (Some(57.0), 75.0),
            "incorrect hot water demand returned"
        );
        assert!(bath
            .hot_water_demand(
                WaterHeatingEvent {
                    volume: None,
                    temperature: 40.0,
                    start: 0.,
                    duration: None
                },
                &func_temp_hot_water,
                simulation_time.iter().current_iteration()
            )
            .is_err());
    }
}
