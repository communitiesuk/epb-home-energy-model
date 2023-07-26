use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::misc::frac_hot_water;

pub struct Bath<'a> {
    size_in_litres: f64,
    cold_water_source: &'a ColdWaterSource,
    flowrate: f64,
    temp_hot: f64,
}

impl<'a> Bath<'a> {
    pub fn new(size_in_litres: f64, cold_water_source: &'a ColdWaterSource, flowrate: f64) -> Self {
        Self {
            size_in_litres,
            cold_water_source,
            flowrate,
            temp_hot: 52.0, // TODO Python code has note to get this from somewhere not hard-coded
        }
    }

    pub fn get_size(&self) -> f64 {
        self.size_in_litres
    }

    pub fn get_cold_water_source(&self) -> &ColdWaterSource {
        self.cold_water_source
    }

    pub fn get_flowrate(&self) -> f64 {
        self.flowrate
    }

    pub fn get_temp_hot(&self) -> f64 {
        self.temp_hot
    }

    /// Calculate volume of hot water required
    /// (and volume of warm water draining to WWHRS, if applicable)
    ///
    /// Arguments:
    /// * `temp_target` - temperature of warm water delivered at tap, in Celsius
    /// * `timestep_idx` - the index of the timestep for which we are querying the hot water demand
    pub fn hot_water_demand(&self, temp_target: f64, timestep_idx: usize) -> f64 {
        let temp_cold = self.cold_water_source.temperature(timestep_idx);

        let vol_warm_water = self.size_in_litres; // may wish to modify the volume of water compared to size of bath

        let vol_hot_water = vol_warm_water * frac_hot_water(temp_target, self.temp_hot, temp_cold);

        vol_hot_water
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation_time::SimulationTime;
    use rstest::*;
    use std::rc::Rc;

    #[rstest]
    pub fn should_get_correct_size() {
        let simulation_time = SimulationTime::new(0.0, 3.0, 1.0);
        let cold_water_source = ColdWaterSource::new(vec![2.0, 3.0, 4.0], &simulation_time, 1.0);
        let bath = Bath::new(100.0, &cold_water_source, 4.5);
        assert_eq!(bath.get_size(), 100.0, "incorrect size of bath returned");
    }

    #[rstest]
    pub fn should_give_cold_water_source() {
        let simulation_time = SimulationTime::new(0.0, 3.0, 1.0);
        let cold_water_source = ColdWaterSource::new(vec![2.0, 3.0, 4.0], &simulation_time, 1.0);
        let bath = Bath::new(100.0, &cold_water_source, 4.5);
        assert_eq!(
            bath.get_cold_water_source(),
            &cold_water_source,
            "cold water source not returned"
        );
    }

    #[rstest]
    pub fn should_get_correct_flowrate() {
        let simulation_time = SimulationTime::new(0.0, 3.0, 1.0);
        let cold_water_source = ColdWaterSource::new(vec![2.0, 3.0, 4.0], &simulation_time, 1.0);
        let bath = Bath::new(100.0, &cold_water_source, 4.5);
        assert_eq!(bath.get_flowrate(), 4.5, "incorrect flow rate returned");
    }

    #[rstest]
    pub fn should_get_correct_hot_water_demand() {
        let simulation_time = Rc::new(SimulationTime::new(0.0, 3.0, 1.0));
        let cold_water_source =
            ColdWaterSource::new(vec![2.0, 3.0, 4.0], &simulation_time.clone(), 1.0);
        let bath = Bath::new(100.0, &cold_water_source, 4.5);
        let expected_demands = [76.0, 75.510, 75.0];
        for idx in 0..simulation_time.total_steps() {
            assert_eq!(
                round_by_precision(bath.hot_water_demand(40.0, idx), 1e3),
                round_by_precision(expected_demands[idx as usize], 1e3),
                "incorrect hot water returned"
            );
        }
    }

    fn round_by_precision(src: f64, precision: f64) -> f64 {
        (precision * src).round() / precision
    }
}
