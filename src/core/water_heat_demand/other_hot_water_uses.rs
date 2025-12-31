use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::misc::calc_fraction_hot_water;
use crate::simulation_time::SimulationTimeIteration;
use std::sync::Arc;

#[derive(Debug)]
pub struct OtherHotWater {
    flowrate: f64,
    cold_water_source: Arc<ColdWaterSource>,
}

impl OtherHotWater {
    pub fn new(flowrate: f64, cold_water_source: Arc<ColdWaterSource>) -> Self {
        Self {
            flowrate,
            cold_water_source,
        }
    }

    #[cfg(test)]
    pub fn get_flowrate(&self) -> f64 {
        self.flowrate
    }

    pub fn get_cold_water_source(&self) -> &ColdWaterSource {
        &self.cold_water_source
    }

    /// Calculate volume of hot water required
    /// (and volume of warm water draining to WWHRS, if applicable)
    ///
    /// Arguments:
    /// * `temp_target` - temperature of warm water delivered at tap/outlet head, in Celcius
    /// * `temp_hot_water`
    /// * `total_demand_duration` - cumulative running time of this event during the current
    ///                             timestep, in minutes
    /// * `simtime` - the iteration of the timestep for which we are querying the hot water demand
    pub fn hot_water_demand(
        &self,
        temp_target: f64,
        temp_hot_water: f64,
        total_demand_duration: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        let temp_cold = self.cold_water_source.temperature(simtime);

        // TODO (from Python) Account for behavioural variation factor fbeh (sic)
        let vol_warm_water = self.flowrate * total_demand_duration;
        // ^^^ litres = litres/minute * minutes

        let vol_hot_water =
            vol_warm_water * calc_fraction_hot_water(temp_target, temp_hot_water, temp_cold)?;

        Ok((vol_hot_water, vol_warm_water))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation_time::SimulationTime;
    use approx::assert_relative_eq;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[rstest]
    pub fn should_give_correct_flowrate() {
        let cold_water_source = ColdWaterSource::new(vec![2.0, 3.0, 4.0], 0, 1.0);
        let other_water = OtherHotWater::new(5.0, cold_water_source.into());
        assert_eq!(
            other_water.get_flowrate(),
            5.0,
            "incorrect flow rate returned"
        );
    }

    #[rstest]
    pub fn should_give_cold_water_source() {
        let cold_water_source = ColdWaterSource::new(vec![2.0, 3.0, 4.0], 0, 1.0);
        let expected_cold_water_source = cold_water_source.clone();
        let other_water = OtherHotWater::new(5.0, cold_water_source.into());
        assert_eq!(
            other_water.get_cold_water_source(),
            &expected_cold_water_source,
            "cold water source not returned"
        );
    }

    #[rstest]
    pub fn should_calculate_correct_hot_water_demand() {
        let simulation_time = SimulationTime::new(0.0, 3.0, 1.0);
        let cold_water_source = ColdWaterSource::new(vec![2.0, 3.0, 4.0], 0, 1.0);
        let other_water = OtherHotWater::new(5.0, cold_water_source.into());
        let expected_demands = [15.2, 15.102, 15.0];
        for (idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                other_water
                    .hot_water_demand(40.0, 52.0, 4.0, t_it)
                    .unwrap()
                    .0,
                expected_demands[idx],
                max_relative = 1e-3
            );
        }
    }
}
