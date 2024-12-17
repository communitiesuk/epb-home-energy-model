use crate::simulation_time::{SimulationTime, SimulationTimeIteration};

/// This module provides objects to represent the source(s) of cold water.

#[derive(Clone, Debug, PartialEq)]
pub struct ColdWaterSource {
    cold_water_temps: Vec<f64>,
    start_day: u32,
    time_series_step: f64,
}

impl ColdWaterSource {
    pub fn new(
        cold_water_temps: Vec<f64>,
        simulation_time: &SimulationTime,
        start_day: u32,
        timestep: f64,
    ) -> Self {
        Self {
            cold_water_temps: cold_water_temps,
            start_day,
            time_series_step: timestep,
        }
    }

    pub fn temperature(&self, simtime: SimulationTimeIteration) -> f64 {
        let time_series_idx = simtime.time_series_idx(self.start_day, self.time_series_step);
        *self.cold_water_temps.get(time_series_idx).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation_time::SimulationTime;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[rstest]
    pub fn should_emit_correct_temperature() {
        let simulation_time = SimulationTime::new(0.0, 8.0, 1.0);
        let water_temps = [0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0];
        let cold_water_source = ColdWaterSource::new(water_temps.into(), &simulation_time, 0, 1.0);
        for (idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                cold_water_source.temperature(t_it),
                water_temps[idx],
                "incorrect water temp returned"
            );
        }
    }
}
