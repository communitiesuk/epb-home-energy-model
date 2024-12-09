use crate::simulation_time::SimulationTimeIteration;
use crate::{core::water_heat_demand::cold_water_source::ColdWaterSource, statistics::np_interp};

/// This module provides types to model waste water heat recovery systems of different kinds.

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
    cold_water_source: ColdWaterSource,
    flow_rates: Vec<f64>,
    efficiencies: Vec<f64>,
    utilisation_factor: f64,
}

/// A class to represent instantaneous waste water heat recovery systems with arrangement B
///     
/// For System B WWHRS, output of the heat exchanger is fed to the shower only
impl WWHRSInstantaneousSystemB {
    pub fn new(
        cold_water_source: ColdWaterSource,
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
    cold_water_source: ColdWaterSource,
    stored_temperature: f64,
    flow_rates: Vec<f64>,
    efficiencies: Vec<f64>,
    utilisation_factor: f64,
}

impl WWHRSInstantaneousSystemC {
    pub fn new(
        flow_rates: Vec<f64>,
        efficiencies: Vec<f64>,
        cold_water_source: ColdWaterSource,
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
    cold_water_source: ColdWaterSource,
    stored_temperature: f64,
    flow_rates: Vec<f64>,
    efficiencies: Vec<f64>,
    utilisation_factor: f64,
}

impl WWHRSInstantaneousSystemA {
    pub fn new(
        flow_rates: Vec<f64>,
        efficiencies: Vec<f64>,
        cold_water_source: ColdWaterSource,
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
        SimulationTime::new(0., 2., 1.)
    }

    #[fixture]
    fn wwhrs_b(simulation_time: SimulationTime) -> WWHRSInstantaneousSystemB {
        let cold_water_source =
            ColdWaterSource::new(vec![17.0, 17.0, 17.0], &simulation_time, 0, 1.0);
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
        let cold_water_source =
            ColdWaterSource::new(vec![17.1, 17.2, 17.3], &simulation_time, 0, 1.0);
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

    #[fixture]
    fn wwhrs_a(simulation_time: SimulationTime) -> WWHRSInstantaneousSystemA {
        let cold_water_source =
            ColdWaterSource::new(vec![17.1, 17.2, 17.3], &simulation_time, 0, 1.0);
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
}
