use crate::simulation_time::SimulationTime;

/// This module provides objects to represent the source(s) of cold water.

#[derive(Clone, Debug, PartialEq)]
pub struct ColdWaterSource {
    cold_water_temps: Vec<f64>,
}

impl ColdWaterSource {
    pub fn new(
        cold_water_temps: Vec<f64>,
        simulation_time: &SimulationTime,
        timestep: f64,
    ) -> Self {
        Self {
            cold_water_temps: cold_water_temps
                .iter()
                .map(|t| vec![*t; (timestep / simulation_time.step) as usize])
                .flatten()
                .collect(),
        }
    }

    pub fn temperature(&self, timestep_idx: usize) -> f64 {
        *self.cold_water_temps.get(timestep_idx).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation_time::SimulationTime;
    use rstest::*;

    #[rstest]
    pub fn should_emit_correct_temperature() {
        let simulation_time = SimulationTime::new(0.0, 8.0, 1.0);
        let water_temps = [0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0];
        let cold_water_source = ColdWaterSource::new(water_temps.into(), &simulation_time, 1.0);
        for (idx, _) in simulation_time.iter().enumerate() {
            assert_eq!(
                cold_water_source.temperature(idx),
                water_temps[idx],
                "incorrect water temp returned"
            );
        }
    }
}
