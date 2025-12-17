use crate::simulation_time::SimulationTimeIteration;
use crate::{core::water_heat_demand::cold_water_source::ColdWaterSource, statistics::np_interp};
use anyhow::anyhow;
use std::sync::Arc;

/// This module provides types to model waste water heat recovery systems of different kinds.
/// Uses a unified WWHRS class that handles all system types (A, B, C).
// Temperature reduction of water during the shower from temp_target
const DELTA_T_SHOWER: f64 = 6.0;

/// A unified class to represent instantaneous waste water heat recovery systems
///
/// This class can handle all three system configurations (A, B, C) based on the
/// system type specified when calling the calculation methods. Each physical WWHRS
/// unit is defined once and can be connected to multiple showers with different
/// configurations.
struct WwhrsInstantaneous {
    cold_water_source: ColdWaterSource,
    flow_rates: Vec<f64>,
    system_a_efficiencies: Vec<f64>,
    system_a_utilisation_factor: Option<f64>,
    system_b_efficiencies: Option<Vec<f64>>,
    system_b_utilisation_factor: Option<f64>,
    system_b_efficiency_factor: Option<f64>,
    system_c_efficiencies: Option<Vec<f64>>,
    system_c_utilisation_factor: Option<f64>,
    system_c_efficiency_factor: Option<f64>,
    stored_temperature: f64,
    last_used_time: Option<f64>,
}

struct PerformanceCalculationResult {
    t_cyl_feed: f64,
    flowrate_hot: Option<f64>,
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
    fn new(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Vec<f64>,
        cold_water_source: ColdWaterSource,
        system_a_utilisation_factor: Option<f64>,
        system_b_efficiencies: Option<Vec<f64>>,
        system_b_utilisation_factor: Option<f64>,
        system_c_efficiencies: Option<Vec<f64>>,
        system_c_utilisation_factor: Option<f64>,
        system_b_efficiency_factor: Option<f64>,
        system_c_efficiency_factor: Option<f64>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<Self> {
        //  Storage for temperature (used by some systems)
        let list_temp_vol =
            cold_water_source.get_temp_cold_water(1.0, simulation_time_iteration)?;
        let stored_temperature = list_temp_vol.iter().map(|(t, v)| t * v).sum::<f64>()
            / list_temp_vol.iter().map(|(_, v)| v).sum::<f64>();
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
            stored_temperature,
            last_used_time,
        })
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
    fn calculate_performance(
        &self,
        system_type: WwhrsType,
        temp_target: f64,
        flowrate_waste_water: f64,
        volume_cold_water: f64,
        temp_hot: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<PerformanceCalculationResult> {
        match system_type {
            WwhrsType::A => self.calculate_system_a(
                temp_target,
                flowrate_waste_water,
                volume_cold_water,
                temp_hot,
                simulation_time_iteration,
            ),
            WwhrsType::B => self.calculate_system_b(
                temp_target,
                flowrate_waste_water,
                volume_cold_water,
                temp_hot,
                simulation_time_iteration,
            ),
            WwhrsType::C => self.calculate_system_c(
                temp_target,
                flowrate_waste_water,
                volume_cold_water,
                temp_hot,
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
        let temp_main = list_temp_vol
            .clone()
            .into_iter()
            .map(|(t, v)| t * v)
            .sum::<f64>()
            / list_temp_vol.into_iter().map(|(_, v)| v).sum::<f64>();

        // Get efficiency for System A
        let efficiency =
            self.get_efficiency_from_flowrate(flowrate_waste_water, WwhrsType::A)? / 100.;
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
            flowrate_hot: flowrate_hot,
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
        let temp_main = list_temp_vol
            .clone()
            .into_iter()
            .map(|(t, v)| t * v)
            .sum::<f64>()
            / list_temp_vol.into_iter().map(|(_, v)| v).sum::<f64>();

        // Determine which approach to use based on available data
        let eta_uf = match self.system_b_efficiencies {
            Some(_) => {
                // Approach 1: Pre-corrected System B data
                if self.system_b_utilisation_factor.is_none() {
                    anyhow::bail!(
                        "system_b_utilisation_factor is required when using system_b_efficiencies"
                    );
                }

                let efficiency_adjusted = self.get_efficiency_from_flowrate(flowrate_waste_water, WwhrsType::B)? / 100.0;
                efficiency_adjusted * self.system_b_utilisation_factor.unwrap()
            }
            None => {
                // Approach 2: Convert from System A data
                if self.system_b_utilisation_factor.is_none()
                    || self.system_b_efficiency_factor.is_none()
                {
                    anyhow::bail!("Both system_b_utilisation_factor and system_b_efficiency_factor are required when converting from System A data");
                }

                let base_efficiency = self.get_efficiency_from_flowrate(flowrate_waste_water, WwhrsType::A)? / 100.0;
                let efficiency_adjusted = base_efficiency * self.system_b_efficiency_factor.unwrap();
                efficiency_adjusted * self.system_b_utilisation_factor.unwrap()
            }
        };

        // Calculate drain temperature
        let temp_drain = temp_target - DELTA_T_SHOWER;

        // Implement algebraic solution from Technical Recommendations
        let temp_pre = if is_close!(temp_hot, temp_target, abs_tol=1e-10) {
            temp_main
        } else {
            let temp = eta_uf * (temp_drain - temp_main) / (temp_hot - temp_target);
            (temp_main + temp_hot * temp) / (1. + temp)
        };

        // For System B: m_hot = flowrate_waste_water * (T_target - T_pre_B) / (temp_hot - T_pre_B)
        let flowrate_hot = if is_close!(temp_hot, temp_pre, abs_tol=1e-10) {
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
        _temp_target: f64,
        _flowrate_waste_water: f64,
        _volume_cold_water: f64,
        _temp_hot: f64,
    ) -> anyhow::Result<PerformanceCalculationResult> {
        match self.system_c_efficiencies {
            Some(_) => {
                // Approach 1: Pre-corrected System C data
                if self.system_c_utilisation_factor.is_none() {
                    anyhow::bail!(
                        "system_c_utilisation_factor is required when using system_c_efficiencies"
                    )
                }
                todo!()
            }
            None => {
                // Approach 2: Convert from System A data
                if self.system_c_utilisation_factor.is_none()
                    || self.system_c_efficiency_factor.is_none()
                {
                    anyhow::bail!("Both system_c_utilisation_factor and system_c_efficiency_factor are required when converting from System A data")
                }
            }
        }
        todo!()
    }

    /// Get the interpolated efficiency from the flowrate for specified system type.
    fn get_efficiency_from_flowrate(
        &self,
        flowrate: f64,
        system_type: WwhrsType,
    ) -> anyhow::Result<f64> {
        let efficiencies = match system_type {
            WwhrsType::A => self.system_a_efficiencies.clone(),
            WwhrsType::B => self
                .system_b_efficiencies
                .clone()
                .ok_or(anyhow!("System B efficiencies expected for WWHR System B"))?,
            WwhrsType::C => self
                .system_c_efficiencies
                .clone()
                .ok_or(anyhow!("System C efficiencies expected for WWHR System C"))?,
        };

        let efficiency = if flowrate <= self.flow_rates[0] {
            efficiencies[0]
        } else if flowrate >= self.flow_rates[self.flow_rates.len() - 1] {
            efficiencies[self.flow_rates.len() - 1]
        } else {
            np_interp(flowrate, &self.flow_rates, &efficiencies)
        };

        Ok(efficiency)
    }
}

pub(crate) enum WwhrsType {
    A,
    B,
    C,
}

#[derive(Clone, Debug)]
pub enum Wwhrs {
    WWHRSInstantaneousSystemB(WWHRSInstantaneousSystemB),
    WWHRSInstantaneousSystemC(WWHRSInstantaneousSystemC),
    WWHRSInstantaneousSystemA(WWHRSInstantaneousSystemA),
}

impl Wwhrs {
    pub fn return_temperature(
        &self,
        temp_target: f64,
        flowrate_waste_water: f64,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        match &self {
            Wwhrs::WWHRSInstantaneousSystemB(system) => {
                system.return_temperature(temp_target, flowrate_waste_water, simtime)
            }
            Wwhrs::WWHRSInstantaneousSystemC(system) => {
                system.return_temperature(temp_target, flowrate_waste_water, simtime)
            }
            Wwhrs::WWHRSInstantaneousSystemA(system) => {
                system.return_temperature(temp_target, flowrate_waste_water, simtime)
            }
        }
    }

    pub fn temperature(&self) -> f64 {
        match self {
            Wwhrs::WWHRSInstantaneousSystemB(_) => {
                unreachable!("A SystemB WWHRS does not expect to have temperature() called on it")
            }
            Wwhrs::WWHRSInstantaneousSystemC(c) => c.temperature(),
            Wwhrs::WWHRSInstantaneousSystemA(a) => a.temperature(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct WWHRSInstantaneousSystemB {
    cold_water_source: Arc<ColdWaterSource>,
    flow_rates: Vec<f64>,
    efficiencies: Vec<f64>,
    utilisation_factor: f64,
}

/// A class to represent instantaneous waste water heat recovery systems with arrangement B
///     
/// For System B WWHRS, output of the heat exchanger is fed to the shower only
impl WWHRSInstantaneousSystemB {
    pub fn new(
        cold_water_source: Arc<ColdWaterSource>,
        flow_rates: Vec<f64>,
        efficiencies: Vec<f64>,
        utilisation_factor: f64,
    ) -> Self {
        Self {
            cold_water_source,
            flow_rates,
            efficiencies,
            utilisation_factor,
        }
    }

    pub fn return_temperature(
        &self,
        temp_target: f64,
        flowrate_waste_water: f64,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        // # TODO (from Python) The cold water flow rate depends on the temperature returned from
        // #      this function, which may create a circularity in the calculation.
        // #      May need to integrate System B into shower module and/or combine
        // #      equations.

        let temp_cold = self.cold_water_source.temperature(simtime);

        // # TODO (from Python) If flowrates have been provided for waste and cold water:
        // #    - Calc heat recovered from waste water. Need to do this per shower
        // #      individually? Need WWHRS_Connection object?
        let wwhrs_efficiency =
            self.get_efficiency_from_flowrate(flowrate_waste_water) * self.utilisation_factor;

        // Calculate temp of pre-heated water based on heat recovered and flow rates
        temp_cold + ((wwhrs_efficiency / 100.0) * (temp_target - temp_cold))
    }

    fn get_efficiency_from_flowrate(&self, flowrate: f64) -> f64 {
        np_interp(flowrate, &self.flow_rates, &self.efficiencies)
    }
}

/// A class to represent instantaneous waste water heat recovery systems with arrangement C
///
/// For System C WWHRS, output of the heat exchanger is fed to the hot water system only
#[derive(Clone, Debug)]
pub struct WWHRSInstantaneousSystemC {
    cold_water_source: Arc<ColdWaterSource>,
    stored_temperature: f64,
    flow_rates: Vec<f64>,
    efficiencies: Vec<f64>,
    utilisation_factor: f64,
}

impl WWHRSInstantaneousSystemC {
    pub(crate) fn new(
        flow_rates: Vec<f64>,
        efficiencies: Vec<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        utilisation_factor: f64,
        initial_simtime: SimulationTimeIteration,
    ) -> Self {
        // assuming the first timestep index is wanted here, but this is unclear!
        let stored_temperature = cold_water_source.temperature(initial_simtime);

        Self {
            cold_water_source,
            stored_temperature,
            flow_rates,
            efficiencies,
            utilisation_factor,
        }
    }

    pub fn set_temperature_for_return(&mut self, water_temperature: f64) {
        self.stored_temperature = water_temperature;
    }

    pub fn return_temperature(
        &self,
        temp_target: f64,
        flowrate_waste_water: f64,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        // # TODO (from Python) The cold water flow rate depends on the temperature returned from
        // #      this function, which may create a circularity in the calculation.
        // #      May need to integrate System B into shower module and/or combine
        // #      equations.

        let temp_cold = self.cold_water_source.temperature(simtime);

        // # TODO (from Python) If flowrates have been provided for waste and cold water:
        // #    - Calc heat recovered from waste water. Need to do this per shower
        // #      individually? Need WWHRS_Connection object?
        let wwhrs_efficiency =
            self.get_efficiency_from_flowrate(flowrate_waste_water) * self.utilisation_factor;

        // Calculate temp of pre-heated water based on heat recovered and flow rates
        temp_cold + ((wwhrs_efficiency / 100.0) * (temp_target - temp_cold))
    }

    fn get_efficiency_from_flowrate(&self, flowrate: f64) -> f64 {
        np_interp(flowrate, &self.flow_rates, &self.efficiencies)
    }

    pub fn temperature(&self) -> f64 {
        self.stored_temperature
    }
}

/// A class to represent instantaneous waste water heat recovery systems with arrangement A
///     
/// For System A WWHRS, output of the heat exchanger is fed to both the shower
/// and the hot water system
#[derive(Clone, Debug)]
pub struct WWHRSInstantaneousSystemA {
    cold_water_source: Arc<ColdWaterSource>,
    stored_temperature: f64,
    flow_rates: Vec<f64>,
    efficiencies: Vec<f64>,
    utilisation_factor: f64,
}

impl WWHRSInstantaneousSystemA {
    pub fn new(
        flow_rates: Vec<f64>,
        efficiencies: Vec<f64>,
        cold_water_source: Arc<ColdWaterSource>,
        utilisation_factor: f64,
        initial_simtime: SimulationTimeIteration,
    ) -> Self {
        // assume that the first timestep is what is wanted here though could be wrong!!
        let stored_temperature = cold_water_source.temperature(initial_simtime);

        Self {
            cold_water_source,
            stored_temperature,
            flow_rates,
            efficiencies,
            utilisation_factor,
        }
    }

    pub fn set_temperature_for_return(&mut self, water_temperature: f64) {
        self.stored_temperature = water_temperature;
    }

    pub fn return_temperature(
        &self,
        temp_target: f64,
        flowrate_waste_water: f64,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        // # TODO (from Python) The cold water flow rate depends on the temperature returned from
        // #      this function, which may create a circularity in the calculation.
        // #      May need to integrate System B into shower module and/or combine
        // #      equations.

        let temp_cold = self.cold_water_source.temperature(simtime);

        // # TODO (from Python) If flowrates have been provided for waste and cold water:
        // #    - Calc heat recovered from waste water. Need to do this per shower
        // #      individually? Need WWHRS_Connection object?
        let wwhrs_efficiency =
            self.get_efficiency_from_flowrate(flowrate_waste_water) * self.utilisation_factor;

        // Calculate temp of pre-heated water based on heat recovered and flow rates
        temp_cold + ((wwhrs_efficiency / 100.0) * (temp_target - temp_cold))
    }

    fn get_efficiency_from_flowrate(&self, flowrate: f64) -> f64 {
        np_interp(flowrate, &self.flow_rates, &self.efficiencies)
    }

    pub fn temperature(&self) -> f64 {
        self.stored_temperature
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation_time::SimulationTime;
    use approx::assert_relative_eq;
    use rstest::*;

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 3., 1.)
    }

    #[fixture]
    fn cold_water_source() -> ColdWaterSource {
        ColdWaterSource::new(vec![17.0, 17.0, 17.0], 0, 1.0)
    }

    #[fixture]
    fn flow_rates() -> Vec<f64> {
        vec![5., 7., 9., 11., 13.]
    }

    #[fixture]
    fn system_a_efficiencies() -> Vec<f64> {
        vec![44.8, 39.1, 34.8, 31.4, 28.6]
    }

    #[fixture]
    fn system_a_utilisation_factor() -> Option<f64> {
        Some(0.7)
    }

    #[rstest]
    fn test_init_with_all_parameters(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Vec<f64>,
        cold_water_source: ColdWaterSource,
        system_a_utilisation_factor: Option<f64>,
        simulation_time: SimulationTime,
    ) {
        let system_b_efficiencies = Some(vec![36.3, 31.7, 28.2, 25.4, 23.2]);
        let system_b_utilisation_factor = Some(0.65);
        let system_c_efficiencies = Some(vec![38.9, 34.0, 30.3, 27.3, 24.8]);
        let system_c_utilisation_factor = Some(0.68);
        let system_b_efficiency_factor = Some(0.81);
        let system_c_efficiency_factor = Some(0.88);
        let simulation_time_iteration = simulation_time.iter().next().unwrap();

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
            simulation_time_iteration,
        )
        .unwrap();

        assert_eq!(
            wwhrs
                .get_efficiency_from_flowrate(5., WwhrsType::A)
                .unwrap(),
            44.8
        );
        assert_eq!(
            wwhrs
                .get_efficiency_from_flowrate(5., WwhrsType::B)
                .unwrap(),
            36.3
        );
        assert_eq!(
            wwhrs
                .get_efficiency_from_flowrate(5., WwhrsType::C)
                .unwrap(),
            38.9
        );
    }

    #[rstest]
    fn test_system_a_missing_utilisation_factor(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Vec<f64>,
        cold_water_source: ColdWaterSource,
        simulation_time: SimulationTime,
    ) {
        let simulation_time_iteration = simulation_time.iter().next().unwrap();

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
            simulation_time_iteration,
        )
        .unwrap();

        assert!(wwhrs
            .calculate_performance(WwhrsType::A, 35., 8., 8., 55., simulation_time_iteration)
            .is_err());
    }

    #[rstest]
    fn test_system_b_conversion_missing_parameters(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Vec<f64>,
        cold_water_source: ColdWaterSource,
        simulation_time: SimulationTime,
    ) {
        let simulation_time_iteration = simulation_time.iter().next().unwrap();

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
            simulation_time_iteration,
        )
            .unwrap();

        assert!(wwhrs
            .calculate_performance(WwhrsType::B, 35., 8., 8., 55., simulation_time_iteration)
            .is_err());
    }

    #[rstest]
    fn test_system_b_pre_corrected_missing_utilisation_factor(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Vec<f64>,
        cold_water_source: ColdWaterSource,
        simulation_time: SimulationTime,
    ) {
        let simulation_time_iteration = simulation_time.iter().next().unwrap();

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
            simulation_time_iteration,
        )
            .unwrap();

        assert!(wwhrs
            .calculate_performance(WwhrsType::B, 35., 8., 8., 55., simulation_time_iteration)
            .is_err());
    }

    #[rstest]
    fn test_system_b_conversion_missing_utilisation_factor_only(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Vec<f64>,
        cold_water_source: ColdWaterSource,
        simulation_time: SimulationTime,
    ) {
        let simulation_time_iteration = simulation_time.iter().next().unwrap();

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
            simulation_time_iteration,
        )
            .unwrap();

        assert!(wwhrs
            .calculate_performance(WwhrsType::B, 35., 8., 8., 55., simulation_time_iteration)
            .is_err());
    }

    #[rstest]
    fn test_system_c_conversion_missing_utilisation_factor_only(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Vec<f64>,
        cold_water_source: ColdWaterSource,
        simulation_time: SimulationTime,
    ) {
        let simulation_time_iteration = simulation_time.iter().next().unwrap();

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
            simulation_time_iteration,
        )
            .unwrap();

        assert!(wwhrs
            .calculate_performance(WwhrsType::C, 35., 8., 8., 55., simulation_time_iteration)
            .is_err());
    }

    #[rstest]
    fn test_system_c_pre_corrected_missing_utilisation_factor(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Vec<f64>,
        cold_water_source: ColdWaterSource,
        simulation_time: SimulationTime,
    ) {
        let simulation_time_iteration = simulation_time.iter().next().unwrap();

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
            simulation_time_iteration,
        )
            .unwrap();

        assert!(wwhrs
            .calculate_performance(WwhrsType::C, 35., 8., 8., 55., simulation_time_iteration)
            .is_err());
    }

    #[rstest]
    fn test_init_with_reduction_factors(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Vec<f64>,
        cold_water_source: ColdWaterSource,
        system_a_utilisation_factor: Option<f64>,
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
            simulation_time_iteration,
        )
            .unwrap();

        // System A should work
        assert_eq!(
            wwhrs
                .get_efficiency_from_flowrate(5., WwhrsType::A)
                .unwrap(),
            44.8
        );

        // Systems B and C should raise errors when no efficiencies provided
        assert!(wwhrs
            .get_efficiency_from_flowrate(5., WwhrsType::B)
            .is_err());
        assert!(wwhrs
            .get_efficiency_from_flowrate(5., WwhrsType::C)
            .is_err());
    }

    #[rstest]
    fn test_calculate_performance_system_a(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Vec<f64>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: ColdWaterSource,
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
            simulation_time_iteration,
        )
            .unwrap();

        let result = wwhrs
            .calculate_performance(WwhrsType::A, 35., 8., 8., 55., simulation_time_iteration)
            .unwrap();

        // For System A, cylinder is not fed pre-heated water
        assert_eq!(result.t_cyl_feed, 20.1038) // Same as cold water
    }

    #[rstest]
    fn test_calculate_performance_system_b(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Vec<f64>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: ColdWaterSource,
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
            simulation_time_iteration,
        )
            .unwrap();

        let result = wwhrs
            .calculate_performance(WwhrsType::B, 35., 8., 8., 55., simulation_time_iteration)
            .unwrap();

        // For System A, cylinder is not fed pre-heated water
        assert_eq!(result.t_cyl_feed, 17.); // Same as cold water
        assert!(result.flowrate_hot.is_some());
    }

    #[rstest]
    fn test_calculate_performance_system_b_division_by_zero(
        flow_rates: Vec<f64>,
        system_a_efficiencies: Vec<f64>,
        system_a_utilisation_factor: Option<f64>,
        cold_water_source: ColdWaterSource,
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
            simulation_time_iteration,
        )
            .unwrap();

        let result = wwhrs
            .calculate_performance(WwhrsType::B, 35., 8., 8., 35., simulation_time_iteration)
            .unwrap();

        assert_eq!(result.flowrate_hot.unwrap(), 8.); // Equal to flowrate_waste_water
    }

    #[fixture]
    fn wwhrs_b() -> WWHRSInstantaneousSystemB {
        let cold_water_source = Arc::from(ColdWaterSource::new(vec![17.0, 17.0, 17.0], 0, 1.0));
        let flow_rates = vec![5., 7., 9., 11., 13.];
        let efficiencies = vec![44.8, 39.1, 34.8, 31.4, 28.6];
        let utilisation_factor = 0.7;

        WWHRSInstantaneousSystemB::new(
            cold_water_source,
            flow_rates,
            efficiencies,
            utilisation_factor,
        )
    }

    #[rstest]
    fn test_return_temperature_for_b(
        wwhrs_b: WWHRSInstantaneousSystemB,
        simulation_time: SimulationTime,
    ) {
        assert_relative_eq!(
            wwhrs_b.return_temperature(35.0, 8.0, simulation_time.iter().current_iteration()),
            21.6557,
            max_relative = 1e-7
        );
    }

    #[rstest]
    fn test_get_efficiency_from_flowrate_for_b(wwhrs_b: WWHRSInstantaneousSystemB) {
        assert_eq!(wwhrs_b.get_efficiency_from_flowrate(5.0), 44.8);
    }

    #[fixture]
    fn wwhrs_c(simulation_time: SimulationTime) -> WWHRSInstantaneousSystemC {
        let cold_water_source = Arc::from(ColdWaterSource::new(vec![17.1, 17.2, 17.3], 0, 1.0));
        let flow_rates = vec![5., 7., 9., 11., 13.];
        let efficiencies = vec![44.8, 39.1, 34.8, 31.4, 28.6];
        let utilisation_factor = 0.7;

        WWHRSInstantaneousSystemC::new(
            flow_rates,
            efficiencies,
            cold_water_source,
            utilisation_factor,
            simulation_time.iter().current_iteration(),
        )
    }
    #[rstest]
    fn test_return_temperature_for_c(
        wwhrs_c: WWHRSInstantaneousSystemC,
        simulation_time: SimulationTime,
    ) {
        assert_relative_eq!(
            wwhrs_c.return_temperature(35.0, 8.0, simulation_time.iter().current_iteration()),
            21.729835,
            max_relative = 1e-7
        );
    }

    #[rstest]
    fn test_get_efficiency_from_flowrate_for_c(wwhrs_c: WWHRSInstantaneousSystemC) {
        assert_relative_eq!(
            wwhrs_c.get_efficiency_from_flowrate(7.0),
            39.1,
            max_relative = 1e-7
        );
    }

    #[rstest]
    fn test_temperature_for_c(
        mut wwhrs_c: WWHRSInstantaneousSystemC,
        simulation_time: SimulationTime,
    ) {
        for (t_idx, _) in simulation_time.iter().enumerate() {
            if t_idx == 2 {
                wwhrs_c.set_temperature_for_return(16.);
            }
            assert_eq!(wwhrs_c.temperature(), [17.1, 17.1, 16.][t_idx]);
        }
    }

    #[fixture]
    fn wwhrs_a(simulation_time: SimulationTime) -> WWHRSInstantaneousSystemA {
        let cold_water_source = Arc::from(ColdWaterSource::new(vec![17.1, 17.2, 17.3], 0, 1.0));
        let flow_rates = vec![5., 7., 9., 11., 13.];
        let efficiencies = vec![44.8, 39.1, 34.8, 31.4, 28.6];
        let utilisation_factor = 0.7;

        WWHRSInstantaneousSystemA::new(
            flow_rates,
            efficiencies,
            cold_water_source,
            utilisation_factor,
            simulation_time.iter().current_iteration(),
        )
    }

    #[rstest]
    fn test_return_temperature_for_a(
        wwhrs_a: WWHRSInstantaneousSystemA,
        simulation_time: SimulationTime,
    ) {
        assert_relative_eq!(
            wwhrs_a.return_temperature(35.0, 8.0, simulation_time.iter().current_iteration()),
            21.729835,
            max_relative = 1e-7
        );
    }

    #[rstest]
    fn test_get_efficiency_from_flowrate_for_a(wwhrs_a: WWHRSInstantaneousSystemA) {
        assert_relative_eq!(
            wwhrs_a.get_efficiency_from_flowrate(8.0),
            36.95,
            max_relative = 1e-7
        );
    }

    #[rstest]
    fn test_temperature_for_a(
        mut wwhrs_a: WWHRSInstantaneousSystemA,
        simulation_time: SimulationTime,
    ) {
        for (t_idx, _) in simulation_time.iter().enumerate() {
            if t_idx == 2 {
                wwhrs_a.set_temperature_for_return(16.);
            }
            assert_eq!(wwhrs_a.temperature(), [17.1, 17.1, 16.][t_idx]);
        }
    }
}
