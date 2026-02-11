use crate::input::WwhrsConfiguration as WwhrsConfigurationInput;
use crate::simulation_time::SimulationTimeIteration;
use crate::{core::water_heat_demand::cold_water_source::ColdWaterSource, statistics::np_interp};
use anyhow::{anyhow, bail};
use atomic_float::AtomicF64;
use fsum::FSum;
use std::sync::atomic::Ordering;
use std::sync::Arc;

/// This module provides types to model waste water heat recovery systems of different kinds.
/// Uses a unified WWHRS class that handles all system types (A, B, C).
// Temperature reduction of water during the shower from temp_target
const DELTA_T_SHOWER: f64 = 6.0;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
/// WWHRS system configuration
/// A - Both shower and water heating system get pre-heated water
/// B - Only shower gets pre-heated water
/// C - Only water heating system gets pre-heated water
pub(crate) enum WwhrsConfiguration {
    AShowerAndWaterHeatingSystem,
    BShower,
    CWaterHeatingSystem,
}

impl From<WwhrsConfigurationInput> for WwhrsConfiguration {
    fn from(value: WwhrsConfigurationInput) -> Self {
        match value {
            WwhrsConfigurationInput::ShowerAndWaterHeatingSystem => {
                Self::AShowerAndWaterHeatingSystem
            }
            WwhrsConfigurationInput::Shower => Self::BShower,
            WwhrsConfigurationInput::WaterHeatingSystem => Self::CWaterHeatingSystem,
        }
    }
}

/// A unified class to represent instantaneous waste water heat recovery systems
///
/// This class can handle all three system configurations (A, B, C) based on the
/// system type specified when calling the calculation methods. Each physical WWHRS
/// unit is defined once and can be connected to multiple showers with different
/// configurations.
#[derive(Debug)]
pub(crate) struct WwhrsInstantaneous {
    cold_water_source: Arc<ColdWaterSource>,
    flow_rates: Vec<f64>,
    system_a_efficiencies: Option<Vec<f64>>,
    system_a_utilisation_factor: Option<f64>,
    system_b_efficiencies: Option<Vec<f64>>,
    system_b_utilisation_factor: Option<f64>,
    system_b_efficiency_factor: Option<f64>,
    system_c_efficiencies: Option<Vec<f64>>,
    system_c_utilisation_factor: Option<f64>,
    system_c_efficiency_factor: Option<f64>,
    preheated_temperature: AtomicF64, // using NaN to represent None here
    preheated_volume_remaining: AtomicF64,
    last_used_time: Option<f64>,
}

pub(crate) struct PerformanceCalculationResult {
    pub(crate) t_cyl_feed: f64,
    pub(crate) flowrate_hot: Option<f64>,
}

impl WwhrsInstantaneous {
    /// Iitialize the WWHRS with efficiency data for all system types.
    ///
    /// Args:
    /// * `flow_rates`: List of test flow rates (e.g., [5, 7, 9, 11, 13])
    /// * `system_a_efficiencies`: Measured efficiencies for System A at test flow rates
    /// * `cold_water_source`: Cold water source object
    /// * `system_a_utilisation_factor`: Utilisation factor for System A
    /// * `system_b_efficiencies`: Efficiencies for System B (optional, will use reduction factor if not provided)
    /// * `system_b_utilisation_factor`: Utilisation factor for System B
    /// * `system_c_efficiencies`: Efficiencies for System C (optional, will use reduction factor if not provided)
    /// * `system_c_utilisation_factor`: Utilisation factor for System C
    /// * `system_b_efficiency_factor`: Reduction factor for System B (default 0.81)
    /// * `system_c_efficiency_factor`: Reduction factor for System C (default 0.88)
    pub(crate) fn new(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        cold_water_source: Arc<ColdWaterSource>,
        system_a_utilisation_factor: Option<f64>,
        system_b_efficiencies: Option<Vec<f64>>,
        system_b_utilisation_factor: Option<f64>,
        system_c_efficiencies: Option<Vec<f64>>,
        system_c_utilisation_factor: Option<f64>,
        system_b_efficiency_factor: Option<f64>,
        system_c_efficiency_factor: Option<f64>,
    ) -> anyhow::Result<Self> {
        let last_used_time = None; // Future expansion - track last use time
        Ok(Self {
            cold_water_source,
            flow_rates,
            system_a_efficiencies,
            system_a_utilisation_factor,
            system_b_efficiencies,
            system_b_utilisation_factor,
            system_b_efficiency_factor,
            system_c_efficiencies,
            system_c_utilisation_factor,
            system_c_efficiency_factor,
            preheated_temperature: AtomicF64::new(f64::NAN), // using NaN to represent None here
            preheated_volume_remaining: AtomicF64::new(0.0),
            last_used_time,
        })
    }

    /// Register pre-heated water available from a shower event.
    ///
    /// Called after a shower event processes through the WWHRS. The pre-heated
    /// water is only available up to the volume of cold water that was drawn
    /// through the heat exchanger during that shower.
    ///
    /// If called multiple times within a timestep (e.g., multiple showers
    /// connected to same WWHRS), volumes are accumulated with weighted average
    /// temperature.
    ///
    /// Args:
    ///   * `temperature` - The pre-heated water temperature (T_cyl_feed)
    ///   * `volume` - Volume of pre-heated water available (litres)
    ///
    /// Returns error if volume is negative
    pub(crate) fn register_preheated_volume(
        &self,
        temperature: f64,
        volume: f64,
    ) -> anyhow::Result<()> {
        if is_close!(volume, 0.0, abs_tol = 1e-10) {
            return Ok(());
        }

        if volume < 0.0 {
            bail!("Volume cannot be negative, got {volume}");
        }

        let current_volume_remaining = self.preheated_volume_remaining.load(Ordering::SeqCst);
        let current_temperature = self.preheated_temperature.load(Ordering::SeqCst);

        if current_volume_remaining > 0. && !current_temperature.is_nan() {
            // Multiple showers in same timestep - weighted average temperature
            let total_volume = self.preheated_volume_remaining.load(Ordering::SeqCst) + volume;
            self.preheated_temperature.store(
                ((current_temperature * current_volume_remaining) + (temperature * volume))
                    / total_volume,
                Ordering::SeqCst,
            );
            self.preheated_volume_remaining
                .store(total_volume, Ordering::SeqCst);
        } else {
            self.preheated_temperature
                .store(temperature, Ordering::SeqCst);
            self.preheated_volume_remaining
                .store(volume, Ordering::SeqCst);
        }

        Ok(())
    }

    /// Reset state at the end of the timestep.
    ///
    /// Pre-heated water is only available while waste water is actively flowing
    /// through the WWHRS heat exchanger. Any pre-heated volume not consumed
    /// during the timestep cannot carry over, as there will be no waste water
    /// flow to maintain the pre-heating.
    ///
    /// This must be called at the end of each timestep to prevent incorrect
    /// energy recovery calculations in subsequent timesteps.
    ///
    /// Note: While draw_off_water() consumes pre-heated volume, it may not
    /// fully deplete it if the total cold water demand is less than the
    /// registered pre-heated volume. This explicit reset ensures correctness
    /// regardless of draw patterns.
    pub(crate) fn timestep_end(&self) {
        self.preheated_temperature.store(f64::NAN, Ordering::SeqCst);
        self.preheated_volume_remaining.store(0.0, Ordering::SeqCst);
    }

    /// Calculate WWHRS performance based on system type.
    ///
    /// Args:
    /// * `system_type`: 'A', 'B', or 'C'
    /// * `temp_target`: Target shower temperature (T_shower)
    /// * `flowrate_waste_water`: Flow rate of waste water
    /// * `temp_hot`: Hot water temperature (required for Systems B and C)
    /// * `volume_cold_water`: Not used in current implementation
    ///
    /// Returns:
    ///
    /// Dictionary with:
    /// * `T_pre`: Pre-heated water temperature
    /// * `T_cyl_feed`: Temperature of water feeding the cylinder
    /// * `m_hot`: Hot water flow rate (if calculable)
    pub(crate) fn calculate_performance(
        &self,
        system_type: WwhrsConfiguration,
        temp_target: f64,
        flowrate_waste_water: f64,
        volume_cold_water: f64,
        temp_hot: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<PerformanceCalculationResult> {
        match system_type {
            WwhrsConfiguration::AShowerAndWaterHeatingSystem => self.calculate_system_a(
                temp_target,
                flowrate_waste_water,
                volume_cold_water,
                temp_hot,
                simulation_time_iteration,
            ),
            WwhrsConfiguration::BShower => self.calculate_system_b(
                temp_target,
                flowrate_waste_water,
                volume_cold_water,
                temp_hot,
                simulation_time_iteration,
            ),
            WwhrsConfiguration::CWaterHeatingSystem => self.calculate_system_c(
                temp_target,
                flowrate_waste_water,
                volume_cold_water,
                temp_hot,
                simulation_time_iteration,
            ),
        }
    }

    /// Calculate performance for System A configuration.
    fn calculate_system_a(
        &self,
        temp_target: f64,
        flowrate_waste_water: f64,
        volume_cold_water: f64,
        temp_hot: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<PerformanceCalculationResult> {
        if self.system_a_utilisation_factor.is_none() {
            anyhow::bail!("system_a_utilisation_factor is required for System A calculation");
        }

        let list_temp_vol = self
            .cold_water_source
            .get_temp_cold_water(volume_cold_water, simulation_time_iteration)?;
        let temp_main = FSum::with_all(list_temp_vol.clone().into_iter().map(|(t, v)| t * v))
            .value()
            / FSum::with_all(list_temp_vol.into_iter().map(|(_, v)| v)).value();

        // Get efficiency for System A
        let efficiency = self.get_efficiency_from_flowrate(
            flowrate_waste_water,
            WwhrsConfiguration::AShowerAndWaterHeatingSystem,
        )? / 100.;
        let eta_uf = efficiency * self.system_a_utilisation_factor.unwrap();

        // Calculate drain temperature
        let temp_drain = temp_target - DELTA_T_SHOWER;

        // For System A: T_pre_A = T_main + η × U_F × (T_drain - T_main)
        let temp_pre = temp_main + eta_uf * (temp_drain - temp_main);

        // For System A: m_hot = flowrate_waste_water * (T_target - T_pre_A) / (temp_hot - T_pre_A)
        let flowrate_hot = if is_close!(temp_hot, temp_pre, abs_tol = 1e-10) {
            None
        } else {
            Some(flowrate_waste_water * (temp_target - temp_pre) / (temp_hot - temp_pre))
        };

        // For System A, both shower and cylinder are fed with pre-heated water
        Ok(PerformanceCalculationResult {
            t_cyl_feed: temp_pre,
            flowrate_hot,
        })
    }

    fn calculate_system_b(
        &self,
        temp_target: f64,
        flowrate_waste_water: f64,
        volume_cold_water: f64,
        temp_hot: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<PerformanceCalculationResult> {
        let list_temp_vol = self
            .cold_water_source
            .get_temp_cold_water(volume_cold_water, simulation_time_iteration)?;
        let temp_main = FSum::with_all(list_temp_vol.clone().into_iter().map(|(t, v)| t * v))
            .value()
            / FSum::with_all(list_temp_vol.into_iter().map(|(_, v)| v)).value();

        // Determine which approach to use based on available data
        let eta_uf = match self.system_b_efficiencies {
            Some(_) => {
                // Approach 1: Pre-corrected System B data
                if self.system_b_utilisation_factor.is_none() {
                    anyhow::bail!(
                        "system_b_utilisation_factor is required when using system_b_efficiencies"
                    );
                }

                let efficiency_adjusted = self.get_efficiency_from_flowrate(
                    flowrate_waste_water,
                    WwhrsConfiguration::BShower,
                )? / 100.0;
                efficiency_adjusted * self.system_b_utilisation_factor.unwrap()
            }
            None => {
                // Approach 2: Convert from System A data
                if self.system_b_utilisation_factor.is_none()
                    || self.system_b_efficiency_factor.is_none()
                {
                    anyhow::bail!("Both system_b_utilisation_factor and system_b_efficiency_factor are required when converting from System A data");
                }

                let base_efficiency = self.get_efficiency_from_flowrate(
                    flowrate_waste_water,
                    WwhrsConfiguration::AShowerAndWaterHeatingSystem,
                )? / 100.0;
                let efficiency_adjusted =
                    base_efficiency * self.system_b_efficiency_factor.unwrap();
                efficiency_adjusted * self.system_b_utilisation_factor.unwrap()
            }
        };

        // Calculate drain temperature
        let temp_drain = temp_target - DELTA_T_SHOWER;

        // Implement algebraic solution from Technical Recommendations
        let temp_pre = if is_close!(temp_hot, temp_target, abs_tol = 1e-10) {
            temp_main
        } else {
            let temp = eta_uf * (temp_drain - temp_main) / (temp_hot - temp_target);
            (temp_main + temp_hot * temp) / (1. + temp)
        };

        // For System B: m_hot = flowrate_waste_water * (T_target - T_pre_B) / (temp_hot - T_pre_B)
        let flowrate_hot = if is_close!(temp_hot, temp_pre, abs_tol = 1e-10) {
            None
        } else {
            Some(flowrate_waste_water * (temp_target - temp_pre) / (temp_hot - temp_pre))
        };

        Ok(PerformanceCalculationResult {
            t_cyl_feed: temp_main,
            flowrate_hot,
        })
    }

    fn calculate_system_c(
        &self,
        temp_target: f64,
        flowrate_waste_water: f64,
        volume_cold_water: f64,
        temp_hot: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<PerformanceCalculationResult> {
        let list_temp_vol = self
            .cold_water_source
            .get_temp_cold_water(volume_cold_water, simulation_time_iteration)?;
        let temp_main = FSum::with_all(list_temp_vol.clone().into_iter().map(|(t, v)| t * v))
            .value()
            / FSum::with_all(list_temp_vol.into_iter().map(|(_, v)| v)).value();

        let eta_uf = match self.system_c_efficiencies {
            Some(_) => {
                // Approach 1: Pre-corrected System C data
                if self.system_c_utilisation_factor.is_none() {
                    anyhow::bail!(
                        "system_c_utilisation_factor is required when using system_c_efficiencies"
                    )
                }

                let efficiency_adjusted = self.get_efficiency_from_flowrate(
                    flowrate_waste_water,
                    WwhrsConfiguration::CWaterHeatingSystem,
                )? / 100.0;
                efficiency_adjusted * self.system_c_utilisation_factor.unwrap()
            }
            None => {
                // Approach 2: Convert from System A data
                if self.system_c_utilisation_factor.is_none()
                    || self.system_c_efficiency_factor.is_none()
                {
                    anyhow::bail!("Both system_c_utilisation_factor and system_c_efficiency_factor are required when converting from System A data")
                }
                let base_efficiency = self.get_efficiency_from_flowrate(
                    flowrate_waste_water,
                    WwhrsConfiguration::AShowerAndWaterHeatingSystem,
                )? / 100.0;
                let efficiency_adjusted =
                    base_efficiency * self.system_c_efficiency_factor.unwrap();
                efficiency_adjusted * self.system_c_utilisation_factor.unwrap()
            }
        };

        // Calculate drain temperature
        let temp_drain = temp_target - DELTA_T_SHOWER;

        // Direct calculation for System C
        let temp_pre = if is_close!(temp_target, temp_main, abs_tol = 1e-10) {
            temp_main
        } else {
            temp_main
                + eta_uf * (temp_drain - temp_main) * (temp_hot - temp_main)
                    / (temp_target - temp_main)
        };

        // For System C: m_hot = flowrate_waste_water * (T_target - T_main) / (temp_hot - T_main)
        let flowrate_hot = if is_close!(temp_hot, temp_main, abs_tol = 1e-10) {
            None
        } else {
            Some(flowrate_waste_water * (temp_target - temp_main) / (temp_hot - temp_main))
        };

        Ok(PerformanceCalculationResult {
            t_cyl_feed: temp_pre,
            flowrate_hot,
        })
    }

    /// Get the interpolated efficiency from the flowrate for specified system type.
    fn get_efficiency_from_flowrate(
        &self,
        flowrate: f64,
        system_type: WwhrsConfiguration,
    ) -> anyhow::Result<f64> {
        let efficiencies = match system_type {
            WwhrsConfiguration::AShowerAndWaterHeatingSystem => {
                self.system_a_efficiencies.clone().ok_or(anyhow!(
                    "System A efficiencies not available - no system_a_efficiencies provided"
                ))?
            }
            WwhrsConfiguration::BShower => self.system_b_efficiencies.clone().ok_or(anyhow!(
                "System B efficiencies not available - no system_b_efficiencies provided"
            ))?,
            WwhrsConfiguration::CWaterHeatingSystem => self.system_c_efficiencies.clone().ok_or(
                anyhow!("System C efficiencies not available - no system_c_efficiencies provided"),
            )?,
        };

        let efficiency = if flowrate <= self.flow_rates[0]
            || is_close!(flowrate, 0., rel_tol = 1e-9, abs_tol = 0.0)
        {
            efficiencies[0]
        } else if flowrate >= self.flow_rates[self.flow_rates.len() - 1]
            || is_close!(
                flowrate,
                self.flow_rates[self.flow_rates.len() - 1],
                rel_tol = 1e-9,
                abs_tol = 0.0
            )
        {
            efficiencies[self.flow_rates.len() - 1]
        } else {
            np_interp(flowrate, &self.flow_rates, &efficiencies)
        };

        Ok(efficiency)
    }

    /// Get the temperature of cold water, accounting for pre-heated availability.
    ///
    /// Returns pre-heated temperature up to the available volume, then mains
    /// temperature for any remainder.
    ///
    /// Args:
    ///   * `volume_needed` - Volume of cold water required (litres)
    ///
    /// Returns:
    ///   List of (temperature, volume) tuples.
    pub(crate) fn get_temp_cold_water(
        &self,
        volume_needed: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        let mut result: Vec<(f64, f64)> = Vec::new();

        if self.preheated_volume_remaining.load(Ordering::SeqCst) > 0.
            && !self.preheated_temperature.load(Ordering::SeqCst).is_nan()
        {
            // Provide re-heated water up to available volume
            let preheated_to_provide =
                volume_needed.min(self.preheated_volume_remaining.load(Ordering::SeqCst));

            if preheated_to_provide > 0. {
                result.push((
                    self.preheated_temperature.load(Ordering::SeqCst),
                    preheated_to_provide,
                ));
            }

            // Remainder from mains
            let remaining = volume_needed - preheated_to_provide;
            if remaining > 0. {
                result.extend(
                    self.cold_water_source
                        .get_temp_cold_water(remaining, simtime)?,
                );
            }
        } else {
            result.extend(
                self.cold_water_source
                    .get_temp_cold_water(volume_needed, simtime)?,
            );
        }

        Ok(result)
    }

    /// Draw water, consuming pre-heated volume if available.
    ///
    /// Unlike get_temp_cold_water(), this method consumes the pre-heated
    /// volume so subsequent draws receive less (or no) pre-heated water.
    ///
    /// Args:
    ///   * `volume_needed` - Volume of water to draw (litres)
    ///
    /// Returns:
    ///   List of (temperature, volume) tuples for the water drawn.
    pub(crate) fn draw_off_water(
        &self,
        volume_needed: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        let mut result: Vec<(f64, f64)> = Vec::new();

        if self.preheated_volume_remaining.load(Ordering::SeqCst) > 0.
            && !self.preheated_temperature.load(Ordering::SeqCst).is_nan()
        {
            // Consume pre-heated water up to available volume
            let preheated_to_use =
                volume_needed.min(self.preheated_volume_remaining.load(Ordering::SeqCst));

            if preheated_to_use > 0. {
                result.push((
                    self.preheated_temperature.load(Ordering::SeqCst),
                    preheated_to_use,
                ));
                self.preheated_volume_remaining
                    .fetch_sub(preheated_to_use, Ordering::SeqCst);
            }

            // Remainder from mains
            let remaining = volume_needed - preheated_to_use;
            if remaining > 0. {
                result.extend(self.cold_water_source.draw_off_water(remaining, simtime)?);
            }
        } else {
            result.extend(
                self.cold_water_source
                    .draw_off_water(volume_needed, simtime)?,
            );
        }

        Ok(result)
    }

    /// Set the time when this WWHRS was last used (for future expansion).
    fn set_last_used_time(&mut self, time: f64) {
        self.last_used_time = Some(time)
    }

    /// Get the time when this WWHRS was last used (for future expansion).
    fn get_last_used_time(&self) -> Option<f64> {
        self.last_used_time
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation_time::SimulationTime;
    use rstest::*;

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 3., 1.)
    }

    #[fixture]
    fn cold_water_source() -> Arc<ColdWaterSource> {
        ColdWaterSource::new(vec![17.0, 17.0, 17.0], 0, 1.0).into()
    }

    #[fixture]
    fn flow_rates() -> Vec<f64> {
        vec![5., 7., 9., 11., 13.]
    }

    #[fixture]
    fn system_a_efficiencies() -> Option<Vec<f64>> {
        vec![44.8, 39.1, 34.8, 31.4, 28.6].into()
    }

    #[fixture]
    fn system_a_utilisation_factor() -> Option<f64> {
        Some(0.7)
    }

    #[rstest]
    fn test_init_with_all_parameters(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        cold_water_source: Arc<ColdWaterSource>,
        system_a_utilisation_factor: Option<f64>,
    ) {
        let system_b_efficiencies = Some(vec![36.3, 31.7, 28.2, 25.4, 23.2]);
        let system_b_utilisation_factor = Some(0.65);
        let system_c_efficiencies = Some(vec![38.9, 34.0, 30.3, 27.3, 24.8]);
        let system_c_utilisation_factor = Some(0.68);
        let system_b_efficiency_factor = Some(0.81);
        let system_c_efficiency_factor = Some(0.88);

        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            system_b_efficiencies,
            system_b_utilisation_factor,
            system_c_efficiencies,
            system_c_utilisation_factor,
            system_b_efficiency_factor,
            system_c_efficiency_factor,
        )
        .unwrap();

        assert_eq!(
            wwhrs
                .get_efficiency_from_flowrate(5., WwhrsConfiguration::AShowerAndWaterHeatingSystem)
                .unwrap(),
            44.8
        );
        assert_eq!(
            wwhrs
                .get_efficiency_from_flowrate(5., WwhrsConfiguration::BShower)
                .unwrap(),
            36.3
        );
        assert_eq!(
            wwhrs
                .get_efficiency_from_flowrate(5., WwhrsConfiguration::CWaterHeatingSystem)
                .unwrap(),
            38.9
        );
    }

    #[rstest]
    fn test_system_a_missing_utilisation_factor(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        assert!(wwhrs
            .calculate_performance(
                WwhrsConfiguration::AShowerAndWaterHeatingSystem,
                35.,
                8.,
                8.,
                55.,
                simulation_time.iter().current_iteration()
            )
            .is_err());
    }

    /// Test that requesting System A efficiency raises error when system_a_efficiencies is None
    #[rstest]
    fn test_get_efficiency_system_a_when_none(
        flow_rates: Vec<f64>,
        cold_water_source: Arc<ColdWaterSource>,
    ) {
        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            None,
            cold_water_source,
            None,
            Some(vec![36.3, 31.7, 28.2, 25.4, 23.2]),
            Some(0.65),
            None,
            None,
            None,
            None,
        )
        .unwrap();

        assert!(wwhrs
            .get_efficiency_from_flowrate(8.0, WwhrsConfiguration::AShowerAndWaterHeatingSystem)
            .is_err_and(|e| {
                e.to_string().contains(
                    "System A efficiencies not available - no system_a_efficiencies provided",
                )
            }));
    }

    //     def test_get_efficiency_system_a_when_none(self):
    //         """Test that requesting System A efficiency raises error when system_a_efficiencies is None"""
    //         wwhrs = WWHRS_Instantaneous(
    //             flow_rates=self.flow_rates,
    //             system_a_efficiencies=None,  # No System A data provided
    //             cold_water_source=self.cold_water_source,
    //             system_b_efficiencies=[36.3, 31.7, 28.2, 25.4, 23.2],
    //             system_b_utilisation_factor=0.65,
    //         )
    //
    //         with self.assertRaises(ValueError) as context:
    //             wwhrs.get_efficiency_from_flowrate(8.0, "A")
    //         self.assertIn(
    //             "System A efficiencies not available - no system_a_efficiencies provided",
    //             str(context.exception),
    //         )

    #[rstest]
    fn test_system_b_conversion_missing_parameters(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            Some(0.7),
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        assert!(wwhrs
            .calculate_performance(
                WwhrsConfiguration::BShower,
                35.,
                8.,
                8.,
                55.,
                simulation_time.iter().current_iteration()
            )
            .is_err());
    }

    #[rstest]
    fn test_system_b_pre_corrected_missing_utilisation_factor(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            Some(0.7),
            Some(vec![36.3, 31.7, 28.2, 25.4, 23.2]),
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        assert!(wwhrs
            .calculate_performance(
                WwhrsConfiguration::BShower,
                35.,
                8.,
                8.,
                55.,
                simulation_time.iter().current_iteration()
            )
            .is_err());
    }

    #[rstest]
    fn test_system_b_conversion_missing_utilisation_factor_only(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            Some(0.7),
            None,
            None,
            None,
            None,
            Some(0.81),
            None,
        )
        .unwrap();

        assert!(wwhrs
            .calculate_performance(
                WwhrsConfiguration::BShower,
                35.,
                8.,
                8.,
                55.,
                simulation_time.iter().current_iteration()
            )
            .is_err());
    }

    #[rstest]
    fn test_system_c_conversion_missing_utilisation_factor_only(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            Some(0.7),
            None,
            None,
            None,
            None,
            None,
            Some(0.88),
        )
        .unwrap();

        assert!(wwhrs
            .calculate_performance(
                WwhrsConfiguration::CWaterHeatingSystem,
                35.,
                8.,
                8.,
                55.,
                simulation_time.iter().current_iteration()
            )
            .is_err());
    }

    #[rstest]
    fn test_system_c_pre_corrected_missing_utilisation_factor(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            Some(0.7),
            None,
            None,
            Some(vec![38.9, 34.0, 30.3, 27.3, 24.8]),
            None,
            None,
            None,
        )
        .unwrap();

        assert!(wwhrs
            .calculate_performance(
                WwhrsConfiguration::CWaterHeatingSystem,
                35.,
                8.,
                8.,
                55.,
                simulation_time.iter().current_iteration()
            )
            .is_err());
    }

    #[rstest]
    fn test_init_with_reduction_factors(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        cold_water_source: Arc<ColdWaterSource>,
        system_a_utilisation_factor: Option<f64>,
    ) {
        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        // System A should work
        assert_eq!(
            wwhrs
                .get_efficiency_from_flowrate(5., WwhrsConfiguration::AShowerAndWaterHeatingSystem)
                .unwrap(),
            44.8
        );

        // Systems B and C should raise errors when no efficiencies provided
        assert!(wwhrs
            .get_efficiency_from_flowrate(5., WwhrsConfiguration::BShower)
            .is_err());
        assert!(wwhrs
            .get_efficiency_from_flowrate(5., WwhrsConfiguration::CWaterHeatingSystem)
            .is_err());
    }

    #[rstest]
    fn test_calculate_performance_system_a(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        let result = wwhrs
            .calculate_performance(
                WwhrsConfiguration::AShowerAndWaterHeatingSystem,
                35.,
                8.,
                8.,
                55.,
                simulation_time.iter().current_iteration(),
            )
            .unwrap();

        // For System A, cylinder is not fed pre-heated water
        assert_eq!(result.t_cyl_feed, 20.1038) // Same as cold water
    }

    #[rstest]
    fn test_calculate_performance_system_b(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            Some(0.65),
            None,
            None,
            Some(0.81),
            None,
        )
        .unwrap();

        let result = wwhrs
            .calculate_performance(
                WwhrsConfiguration::BShower,
                35.,
                8.,
                8.,
                55.,
                simulation_time.iter().current_iteration(),
            )
            .unwrap();

        // For System B, both shower and cylinder get pre-heated water
        assert_eq!(result.t_cyl_feed, 17.); // Same as cold water
        assert!(result.flowrate_hot.is_some());
    }

    #[rstest]
    fn test_calculate_performance_system_b_division_by_zero(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            Some(0.65),
            None,
            None,
            Some(0.81),
            None,
        )
        .unwrap();

        let result = wwhrs
            .calculate_performance(
                WwhrsConfiguration::BShower,
                35.,
                8.,
                8.,
                35.,
                simulation_time.iter().current_iteration(),
            )
            .unwrap();

        assert_eq!(result.flowrate_hot.unwrap(), 8.); // Equal to flowrate_waste_water
    }

    #[rstest]
    fn test_calculate_performance_system_c(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            None,
            None,
            Some(0.68),
            None,
            Some(0.88),
        )
        .unwrap();

        let result = wwhrs
            .calculate_performance(
                WwhrsConfiguration::CWaterHeatingSystem,
                35.,
                8.,
                8.,
                55.,
                simulation_time.iter().current_iteration(),
            )
            .unwrap();

        // For System C, only cylinder gets pre-heated water
        assert_eq!(result.t_cyl_feed, 22.601422933333335);
        assert!(result.flowrate_hot.is_some());
    }

    #[rstest]
    fn test_calculate_performance_system_c_division_by_zero(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        simulation_time: SimulationTime,
    ) {
        // Create cold water source with temp matching target
        let cold_water_source = ColdWaterSource::new(vec![35.0, 35.0, 35.0], 0, 1.0);
        let simulation_time_iteration = simulation_time.iter().next().unwrap();

        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source.into(),
            system_a_utilisation_factor,
            None,
            None,
            None,
            Some(0.68),
            None,
            Some(0.88),
        )
        .unwrap();

        let result = wwhrs
            .calculate_performance(
                WwhrsConfiguration::CWaterHeatingSystem,
                35.,
                8.,
                8.,
                55.,
                simulation_time_iteration,
            )
            .unwrap();

        // Should return temp_main when division by zero would occur
        assert_eq!(result.t_cyl_feed, 35.); // Equal to flowrate_waste_water
    }

    #[rstest]
    fn test_get_efficiency_from_flowrate_interpolation(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
    ) {
        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        // Test exact match
        assert_eq!(
            wwhrs
                .get_efficiency_from_flowrate(5.0, WwhrsConfiguration::AShowerAndWaterHeatingSystem)
                .unwrap(),
            44.8
        );
        assert_eq!(
            wwhrs
                .get_efficiency_from_flowrate(
                    13.0,
                    WwhrsConfiguration::AShowerAndWaterHeatingSystem
                )
                .unwrap(),
            28.6
        );

        // Test interpolation
        let efficiency_8 = wwhrs
            .get_efficiency_from_flowrate(8.0, WwhrsConfiguration::AShowerAndWaterHeatingSystem)
            .unwrap();
        assert!(efficiency_8 > 34.8); // > efficiency at 9
        assert!(efficiency_8 < 39.1); // < efficiency at 7

        // Test below minimum
        assert_eq!(
            wwhrs
                .get_efficiency_from_flowrate(3.0, WwhrsConfiguration::AShowerAndWaterHeatingSystem)
                .unwrap(),
            44.8
        );

        // Test above maximum
        assert_eq!(
            wwhrs
                .get_efficiency_from_flowrate(
                    15.0,
                    WwhrsConfiguration::AShowerAndWaterHeatingSystem
                )
                .unwrap(),
            28.6
        );
    }

    #[rstest]
    fn test_get_efficiency_from_flowrate_all_systems(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
    ) {
        let system_b_efficiencies = vec![36.3, 31.7, 28.2, 25.4, 23.2];
        let system_c_efficiencies = vec![38.9, 34.0, 30.3, 27.3, 24.8];

        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            Some(system_b_efficiencies),
            None,
            Some(system_c_efficiencies),
            None,
            None,
            None,
        )
        .unwrap();

        // Test all systems
        assert_eq!(
            wwhrs
                .get_efficiency_from_flowrate(5.0, WwhrsConfiguration::AShowerAndWaterHeatingSystem)
                .unwrap(),
            44.8
        );
        assert_eq!(
            wwhrs
                .get_efficiency_from_flowrate(5.0, WwhrsConfiguration::BShower)
                .unwrap(),
            36.3
        );
        assert_eq!(
            wwhrs
                .get_efficiency_from_flowrate(5.0, WwhrsConfiguration::CWaterHeatingSystem)
                .unwrap(),
            38.9
        );
    }

    #[rstest]
    fn test_temperature_methods(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let simtime = simulation_time.iter().next().unwrap();

        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        // // Initial temperature should be from cold water source (no pre-heated water)
        assert_eq!(
            wwhrs.get_temp_cold_water(10.0, simtime).unwrap(),
            vec![(17.0, 10.0)]
        );

        // Test registering pre-heated volume
        wwhrs.register_preheated_volume(20.0, 12.0).unwrap();
        assert_eq!(
            wwhrs.get_temp_cold_water(12.0, simtime).unwrap(),
            vec![(20.0, 12.0)]
        );

        // Test that requesting more than registered returns mix of pre-heated and mains
        wwhrs.register_preheated_volume(20.0, 5.0).unwrap();
        let result = wwhrs.get_temp_cold_water(12.0, simtime).unwrap();
        assert_eq!(result[0], (20.0, 12.0));
    }

    #[rstest]
    fn test_last_used_time_tracking(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
    ) {
        let mut wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        // Initially should be None
        assert!(wwhrs.get_last_used_time().is_none());

        // Set last used time
        let test_time = 123.45;
        wwhrs.set_last_used_time(test_time);
        assert_eq!(wwhrs.get_last_used_time().unwrap(), test_time);

        // Update to new time
        let new_time = 456.78;
        wwhrs.set_last_used_time(new_time);
        assert_eq!(wwhrs.get_last_used_time().unwrap(), new_time);
    }

    #[rstest]
    fn test_system_b_with_specific_efficiencies(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let simulation_time_iteration = simulation_time.iter().next().unwrap();
        let system_b_efficiencies = vec![36.3, 31.7, 28.2, 25.4, 23.2];

        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            Some(system_b_efficiencies),
            Some(0.65),
            None,
            None,
            None,
            None,
        )
        .unwrap();

        wwhrs
            .calculate_performance(
                WwhrsConfiguration::BShower,
                35.,
                8.,
                8.,
                55.,
                simulation_time_iteration,
            )
            .unwrap();

        // Verify it uses specific efficiencies, not reduction factor
        // Efficiency at 8 L/min should be interpolated from system_b_efficiencies
        let efficiency_8 = wwhrs
            .get_efficiency_from_flowrate(8.0, WwhrsConfiguration::BShower)
            .unwrap();
        assert_ne!(
            efficiency_8,
            wwhrs
                .get_efficiency_from_flowrate(8.0, WwhrsConfiguration::AShowerAndWaterHeatingSystem)
                .unwrap()
        );
    }

    #[rstest]
    fn test_system_c_with_specific_efficiencies(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let simulation_time_iteration = simulation_time.iter().next().unwrap();
        let system_c_efficiencies = vec![38.9, 34.0, 30.3, 27.3, 24.8];

        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            None,
            Some(system_c_efficiencies),
            Some(0.68),
            None,
            None,
        )
        .unwrap();

        wwhrs
            .calculate_performance(
                WwhrsConfiguration::CWaterHeatingSystem,
                35.,
                8.,
                8.,
                55.,
                simulation_time_iteration,
            )
            .unwrap();

        // Verify it uses specific efficiencies, not reduction factor
        let efficiency_8 = wwhrs
            .get_efficiency_from_flowrate(8.0, WwhrsConfiguration::CWaterHeatingSystem)
            .unwrap();
        assert_ne!(
            efficiency_8,
            wwhrs
                .get_efficiency_from_flowrate(8.0, WwhrsConfiguration::AShowerAndWaterHeatingSystem)
                .unwrap()
        );
    }

    #[rstest]
    fn test_edge_cases_flowrate_cold_water(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let simulation_time_iteration = simulation_time.iter().next().unwrap();

        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            Some(0.65),
            None,
            None,
            Some(0.81),
            None,
        )
        .unwrap();

        // System B uses flowrate_cold_water if provided
        let result = wwhrs
            .calculate_performance(
                WwhrsConfiguration::BShower,
                35.,
                8.,
                6.,
                55.,
                simulation_time_iteration,
            )
            .unwrap();

        assert_eq!(result.flowrate_hot.unwrap(), 3.297999789473684) // Should use provided value
    }

    #[test]
    fn test_module_constant() {
        assert_eq!(DELTA_T_SHOWER, 6.0);
    }

    #[rstest]
    fn test_system_a_temp_hot_equals_temp_pre(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let simulation_time_iteration = simulation_time.iter().next().unwrap();

        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        // Set up conditions where temp_hot will equal temp_pre
        let result = wwhrs
            .calculate_performance(
                WwhrsConfiguration::AShowerAndWaterHeatingSystem,
                35.,
                8.,
                8.,
                20.1038,
                simulation_time_iteration,
            )
            .unwrap();

        assert!(result.flowrate_hot.is_none())
    }

    // skipping test_system_b_missing_temp_hot as we decided to type temp_hot as required (not as an Option)
    // skipping test_system_b_temp_hot_equals_temp_pre as patch difficult to replicate in Rust
    // skipping test_system_c_missing_temp_hot as we decided to type temp_hot as required (not as an Option)

    #[rstest]
    fn test_system_c_temp_hot_equals_temp_main(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        simulation_time: SimulationTime,
    ) {
        let simulation_time_iteration = simulation_time.iter().next().unwrap();
        let cold_water_source = ColdWaterSource::new(vec![55.0, 55.0, 55.0], 0, 1.0).into();

        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            None,
            None,
            Some(0.88),
            None,
            Some(0.68),
        )
        .unwrap();

        let result = wwhrs
            .calculate_performance(
                WwhrsConfiguration::CWaterHeatingSystem,
                35.,
                8.,
                8.,
                55.,
                simulation_time_iteration,
            )
            .unwrap();

        assert!(result.flowrate_hot.is_none())
    }

    #[rstest]
    fn test_draw_off_water_method(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let simtime = simulation_time.iter().next().unwrap();

        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        // With no pre-heated water, should return mains temperature
        let result = wwhrs.draw_off_water(15.0, simtime).unwrap();
        assert_eq!(result, vec![(17.0, 15.0)]);

        // Register pre-heated volume and draw less than available
        wwhrs.register_preheated_volume(25.0, 30.0).unwrap();
        let result = wwhrs.draw_off_water(20.0, simtime).unwrap();
        assert_eq!(result, vec![(25.0, 20.0)]);

        // Draw again - should get remaining 10L pre-heated, then mains
        let result = wwhrs.draw_off_water(15.0, simtime).unwrap();
        assert_eq!(result[0], (25.0, 10.0)); // Remaining pre-heated
        assert_eq!(result[1], (17.0, 5.0)); // Rest from mains

        // Draw again - all pre-heated exhausted, should get mains only
        let result = wwhrs.draw_off_water(10.0, simtime).unwrap();
        assert_eq!(result, vec![(17.0, 10.0)]);
    }

    /// Test that multiple registrations accumulate with weighted average temperature
    #[rstest]
    fn test_register_preheated_volume_accumulation(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        // Register first volume
        wwhrs.register_preheated_volume(20.0, 10.0).unwrap();

        // Register second volume with different size - should accumulate with weighted average
        wwhrs.register_preheated_volume(35.0, 20.0).unwrap();

        // Weighted average: (20*10 + 35*20) / 30 = (200 + 700) / 30 = 30.0
        // Total volume: 30.0
        let result = wwhrs
            .get_temp_cold_water(30.0, simulation_time.iter().current_iteration())
            .unwrap();
        assert_eq!(result, vec![(30.0, 30.0)]);
    }

    /// Test that registering zero volume has no effect
    #[rstest]
    fn test_register_zero_volume(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        // Register zero volume
        wwhrs.register_preheated_volume(25.0, 0.0).unwrap();

        // Should return mains temperature (no pre-heated water registered)
        let result = wwhrs
            .get_temp_cold_water(10.0, simulation_time.iter().current_iteration())
            .unwrap();
        assert_eq!(result, vec![(17.0, 10.0)]);
    }

    /// Test that get_temp_cold_water does not consume pre-heated volume
    #[rstest]
    fn test_preheated_volume_not_consumed_by_get_temp(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        wwhrs.register_preheated_volume(22.0, 15.0).unwrap();

        // Call get_temp_cold_water multiple times
        let result1 = wwhrs
            .get_temp_cold_water(10.0, simulation_time.iter().current_iteration())
            .unwrap();
        let result2 = wwhrs
            .get_temp_cold_water(10.0, simulation_time.iter().current_iteration())
            .unwrap();

        // Both should return pre-heated temperature (volume not consumed)
        assert_eq!(result1, vec![(22.0, 10.0)]);
        assert_eq!(result2, vec![(22.0, 10.0)]);
    }

    /// Test that draw_off_water consumes pre-heated volume
    #[rstest]
    fn test_draw_off_consumes_preheated_volume(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let simtime = simulation_time.iter().current_iteration();

        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        wwhrs.register_preheated_volume(22.0, 15.0).unwrap();

        // First draw consumes pre-heated volume
        let result1 = wwhrs.draw_off_water(10.0, simtime).unwrap();
        assert_eq!(result1, vec![(22.0, 10.0)]);

        // Second draw - only 5L pre-heated remaining
        let result2 = wwhrs.draw_off_water(10.0, simtime).unwrap();
        assert_eq!(result2[0], (22.0, 5.0)); // Remaining pre-heated
        assert_eq!(result2[1], (17.0, 5.0)); // Rest from mains
    }

    /// Test that registering negative volume raises error
    #[rstest]
    fn test_register_negative_volume_raises_error(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
    ) {
        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        assert!(wwhrs
            .register_preheated_volume(25.0, -5.0)
            .is_err_and(|e| e.to_string().contains("cannot be negative")));
    }

    /// Test that timestep_end resets pre-heated water availability
    #[rstest]
    fn test_timestep_end_resets_preheated_volume(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let simtime = simulation_time.iter().current_iteration();

        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        // Register pre-heated volume
        wwhrs.register_preheated_volume(25.0, 20.0).unwrap();

        // Verify pre-heated volume is available
        let result = wwhrs.get_temp_cold_water(10.0, simtime).unwrap();
        assert_eq!(result, vec![(25.0, 10.0)]);

        // Call timestep_end to reset
        wwhrs.timestep_end();

        // After reset, should return mains temperature
        let result = wwhrs.get_temp_cold_water(10.0, simtime).unwrap();
        assert_eq!(result, vec![(17.0, 10.0)]);
    }

    /// Test get_temp_cold_water when requesting more than pre-heated volume available
    #[rstest]
    fn test_get_temp_cold_water_partial_preheated(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Option<Vec<f64>>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time: SimulationTime,
    ) {
        let simtime = simulation_time.iter().current_iteration();

        let wwhrs = WwhrsInstantaneous::new(
            flow_rates,
            system_a_efficiencies,
            cold_water_source,
            system_a_utilisation_factor,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .unwrap();

        // Register limited pre-heated volume
        wwhrs.register_preheated_volume(22.0, 8.0).unwrap();

        // Request more than available - should get mix of pre-heated and mains
        let result = wwhrs.get_temp_cold_water(15.0, simtime).unwrap();

        // First 8L at pre-heated temp, remaining 7L at mains temp (17.0)
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], (22.0, 8.0));
        assert_eq!(result[1], (17.0, 7.0));
    }

    mod test_wwhrs_integration_scenarios {
        use super::*;
        use indexmap::IndexMap;

        #[fixture]
        fn simulation_time() -> SimulationTime {
            SimulationTime::new(0., 3., 1.)
        }

        #[fixture]
        fn cold_water_source() -> Arc<ColdWaterSource> {
            ColdWaterSource::new(vec![5.0, 10.0, 15.0], 0, 1.0).into()
        }

        #[fixture]
        fn flow_rates() -> Vec<f64> {
            vec![5., 7., 9., 11., 13.]
        }

        #[fixture]
        fn efficiencies() -> Option<Vec<f64>> {
            vec![44.8, 39.1, 34.8, 31.4, 28.6].into()
        }

        #[rstest]
        fn test_realistic_shower_scenario(
            flow_rates: Vec<f64>,
            efficiencies: Option<Vec<f64>>,
            cold_water_source: Arc<ColdWaterSource>,
            simulation_time: SimulationTime,
        ) {
            let simulation_time_iteration = simulation_time.iter().next().unwrap();
            let wwhrs = WwhrsInstantaneous::new(
                flow_rates,
                efficiencies,
                cold_water_source,
                Some(0.7),
                None,
                Some(0.65),
                None,
                Some(0.68),
                Some(0.81),
                Some(0.88),
            )
            .unwrap();

            // Typical shower parameters
            let shower_temp = 38.;
            let shower_flow = 10.0;
            let hot_water_temp = 60.0;

            // Test all three systems
            let mut results: IndexMap<WwhrsConfiguration, PerformanceCalculationResult> =
                IndexMap::new();
            for system in [
                WwhrsConfiguration::AShowerAndWaterHeatingSystem,
                WwhrsConfiguration::BShower,
                WwhrsConfiguration::CWaterHeatingSystem,
            ] {
                let result = wwhrs
                    .calculate_performance(
                        system,
                        shower_temp,
                        shower_flow,
                        8.0,
                        hot_water_temp,
                        simulation_time_iteration,
                    )
                    .unwrap();
                results.insert(system, result);
            }

            assert!(results[0].flowrate_hot.is_some());
            assert!(results[1].flowrate_hot.is_some());
            assert!(results[2].flowrate_hot.is_some());
        }

        #[rstest]
        fn test_extreme_conditions(
            flow_rates: Vec<f64>,
            efficiencies: Option<Vec<f64>>,
            simulation_time: SimulationTime,
        ) {
            let simulation_time_iteration = simulation_time.iter().next().unwrap();

            // Very low temperature difference
            let cold_water_source_warm =
                ColdWaterSource::new(vec![35.0, 35.0, 35.0], 0, 1.0).into();
            let wwhrs_warm = WwhrsInstantaneous::new(
                flow_rates,
                efficiencies,
                cold_water_source_warm,
                Some(0.7),
                None,
                None,
                None,
                None,
                None,
                None,
            )
            .unwrap();

            let result_low_diff = wwhrs_warm
                .calculate_performance(
                    WwhrsConfiguration::AShowerAndWaterHeatingSystem,
                    36.0, // Only 1°C above cold water
                    8.0,
                    8.0,
                    55.0,
                    simulation_time_iteration,
                )
                .unwrap();

            // Should still calculate but effect will be minimal
            assert_eq!(result_low_diff.t_cyl_feed, 33.70675);
        }
    }
}
