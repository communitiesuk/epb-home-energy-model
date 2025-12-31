use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::misc::volume_hot_water_required;
use crate::input::WaterHeatingEvent;
use crate::simulation_time::SimulationTimeIteration;
use std::sync::Arc;

#[derive(Debug)]
pub struct OtherHotWater {
    flowrate: f64,
    cold_water_source: Arc<ColdWaterSource>,
}

impl OtherHotWater {
    pub(crate) fn new(flowrate: f64, cold_water_source: Arc<ColdWaterSource>) -> Self {
        Self {
            flowrate,
            cold_water_source,
        }
    }

    pub(crate) fn get_cold_water_source(&self) -> &ColdWaterSource {
        &self.cold_water_source
    }

    /// Calculate volume of hot water required
    /// (and volume of warm water draining to WWHRS, if applicable)
    pub(crate) fn hot_water_demand<T: Fn(f64) -> f64>(
        &self,
        event: &WaterHeatingEvent,
        func_temp_hot_water: T,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(Option<f64>, f64)> {
        let total_demand_duration = event.duration;
        let temperature_target = event.temperature;

        // TODO (from Python) Account for behavioural variation factor fbeh (sic)
        let volume_warm_water = self.flowrate * total_demand_duration.unwrap_or(0.0);

        let volume_hot_water = volume_hot_water_required(
            volume_warm_water,
            temperature_target,
            func_temp_hot_water,
            |x, simtime| self.cold_water_source.get_temp_cold_water(x, simtime),
            simtime,
        )?;

        if let Some(volume_hot_water) = volume_hot_water {
            let volume_cold_water = volume_warm_water - volume_hot_water;
            self.cold_water_source
                .draw_off_water(volume_cold_water, simtime)?;
        }

        Ok((volume_hot_water, volume_warm_water))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation_time::SimulationTime;
    use approx::assert_relative_eq;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0.0, 3.0, 1.0)
    }

    #[fixture]
    fn cold_water_source() -> ColdWaterSource {
        ColdWaterSource::new(vec![2.0, 3.0, 4.0], 0, 1.0)
    }

    #[fixture]
    fn other_hot_water(cold_water_source: ColdWaterSource) -> OtherHotWater {
        OtherHotWater::new(5.0, cold_water_source.into())
    }

    #[rstest]
    fn test_get_cold_water_source(
        cold_water_source: ColdWaterSource,
        other_hot_water: OtherHotWater,
    ) {
        let expected_cold_water_source = cold_water_source.clone();
        assert_eq!(
            other_hot_water.get_cold_water_source(),
            &expected_cold_water_source,
            "cold water source not returned"
        );
    }

    fn func_temp_hot_water_fixed(_t: f64) -> f64 {
        52.0
    }

    #[rstest]
    #[ignore = "until completed migration to 1.0.0a1"]
    fn test_hot_water_demand(simulation_time: SimulationTime, other_hot_water: OtherHotWater) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let result = other_hot_water
                .hot_water_demand(
                    &WaterHeatingEvent {
                        start: 0.,
                        temperature: 40.0,
                        duration: Some(4.0),
                        volume: None,
                    },
                    func_temp_hot_water_fixed,
                    t_it,
                )
                .unwrap()
                .0
                .unwrap();
            assert_relative_eq!(result, [15.2, 15.102, 15.0][t_idx], max_relative = 1e-3);
        }
    }
}
