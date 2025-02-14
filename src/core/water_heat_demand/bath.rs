use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::misc::frac_hot_water;
use crate::simulation_time::SimulationTimeIteration;
use std::sync::Arc;

#[derive(Debug)]
pub struct Bath {
    size_in_litres: f64,
    cold_water_source: Arc<ColdWaterSource>,
    flowrate: f64,
}

impl Bath {
    pub fn new(
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

    pub fn get_size(&self) -> f64 {
        self.size_in_litres
    }

    pub fn get_cold_water_source(&self) -> &ColdWaterSource {
        &self.cold_water_source
    }

    pub fn get_flowrate(&self) -> f64 {
        self.flowrate
    }

    /// Calculate volume of hot water required
    /// (and volume of warm water draining to WWHRS, if applicable)
    ///
    /// Arguments:
    /// * `temp_target` - temperature of warm water delivered at tap, in Celsius
    /// * `temp_hot_water`
    /// * `simtime` - the iteration of the timestep for which we are querying the hot water demand
    pub fn hot_water_demand(
        &self,
        temp_target: f64,
        temp_hot_water: f64,
        vol_warm_water: f64,
        simtime: SimulationTimeIteration,
    ) -> (f64, f64) {
        let temp_cold = self.cold_water_source.temperature(simtime);

        let vol_warm_water = vol_warm_water.min(self.size_in_litres);

        let vol_hot_water = vol_warm_water * frac_hot_water(temp_target, temp_hot_water, temp_cold);

        (vol_hot_water, vol_warm_water)
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
    pub fn should_get_correct_size(bath: Bath) {
        assert_eq!(bath.get_size(), 100.0, "incorrect size of bath returned");
    }

    #[rstest]
    pub fn should_give_cold_water_source(bath: Bath) {
        assert_eq!(
            bath.get_cold_water_source(),
            bath.cold_water_source.as_ref(),
            "cold water source not returned"
        );
    }

    #[rstest]
    pub fn should_get_correct_flowrate(bath: Bath) {
        assert_eq!(bath.get_flowrate(), 4.5, "incorrect flow rate returned");
    }

    #[rstest]
    pub fn test_hot_water_demand(bath: Bath) {
        let simulation_time = SimulationTime::new(0.0, 3.0, 1.0);

        assert_eq!(
            bath.hot_water_demand(40.0, 52.0, 75.0, simulation_time.iter().next().unwrap()),
            (57.0, 75.0),
            "incorrect hot water demand returned"
        );
        assert_eq!(
            bath.hot_water_demand(40.0, 52.0, 200.0, simulation_time.iter().next().unwrap()),
            (76.0, 100.0),
            "incorrect hot water demand returned for bath fill volume > bath tub volume"
        );
    }
}
