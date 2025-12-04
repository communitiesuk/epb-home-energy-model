use crate::simulation_time::SimulationTimeIteration;
use anyhow::anyhow;

/// This module provides objects to represent the source(s) of cold water.

#[derive(Clone, Debug, PartialEq)]
pub struct ColdWaterSource {
    cold_water_temps: Vec<f64>,
    start_day: u32,
    time_series_step: f64,
}

impl ColdWaterSource {
    /// Construct a ColdWaterSource object
    /// Arguments:
    ///     cold_water_temps: list of cold water temperatures, in deg C (one entry per hour)
    ///     start_day: first day of the time series, day of the year, 0 to 365 (single value)
    ///     time_series_step: timestep of the time series data, in hours
    pub fn new(cold_water_temps: Vec<f64>, start_day: u32, timestep: f64) -> Self {
        Self {
            cold_water_temps,
            start_day,
            time_series_step: timestep,
        }
    }

    // TODO remove temperature function once migration to 1.0.0a1 is complete
    pub fn temperature(&self, simtime: SimulationTimeIteration) -> f64 {
        let time_series_idx = simtime.time_series_idx(self.start_day, self.time_series_step);
        *self.cold_water_temps.get(time_series_idx).unwrap()
    }

    /// Return the cold water temperature and volume for the current timestep
    pub(crate) fn get_temp_cold_water(
        &self,
        volume_needed: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        let time_series_idx = simtime.time_series_idx(self.start_day, self.time_series_step);
        let temperature = *self.cold_water_temps.get(time_series_idx).ok_or_else(|| {
            anyhow!(
                "Cold water temperature not found for index {}",
                time_series_idx
            )
        })?;

        Ok(vec![(temperature, volume_needed)])
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
mod tests {
    use super::*;
    use crate::simulation_time::SimulationTime;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_temperature() {
        let simulation_time = SimulationTime::new(0.0, 8.0, 1.0);
        let water_temps = [0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0];
        let cold_water_source = ColdWaterSource::new(water_temps.into(), 0, 1.0);
        for (idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                cold_water_source.draw_off_water(10., t_it).unwrap(),
                vec!((water_temps[idx], 10.)),
                "incorrect water temp returned"
            );
            assert_eq!(
                cold_water_source.get_temp_cold_water(10., t_it).unwrap(),
                vec!((water_temps[idx], 10.)),
                "incorrect water temp returned"
            );
        }
    }
}
