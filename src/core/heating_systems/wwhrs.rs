use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use interp::interp;

/// This module provides types to model waste water heat recovery systems of different kinds.

#[derive(Clone)]
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
        timestep_idx: usize,
    ) -> f64 {
        match &self {
            Wwhrs::WWHRSInstantaneousSystemB(system) => {
                system.return_temperature(temp_target, flowrate_waste_water, timestep_idx)
            }
            Wwhrs::WWHRSInstantaneousSystemC(system) => {
                system.return_temperature(temp_target, flowrate_waste_water, timestep_idx)
            }
            Wwhrs::WWHRSInstantaneousSystemA(system) => {
                system.return_temperature(temp_target, flowrate_waste_water, timestep_idx)
            }
        }
    }
}

#[derive(Clone)]
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
        timestep_idx: usize,
    ) -> f64 {
        // # TODO (from Python) The cold water flow rate depends on the temperature returned from
        // #      this function, which may create a circularity in the calculation.
        // #      May need to integrate System B into shower module and/or combine
        // #      equations.

        let temp_cold = self.cold_water_source.temperature(timestep_idx);

        // # TODO If flowrates have been provided for waste and cold water:
        // #    - Calc heat recovered from waste water. Need to do this per shower
        // #      individually? Need WWHRS_Connection object?
        let wwhrs_efficiency =
            self.get_efficiency_from_flowrate(flowrate_waste_water) * self.utilisation_factor;

        // Calculate temp of pre-heated water based on heat recovered and flow rates
        temp_cold + ((wwhrs_efficiency / 100.0) * (temp_target - temp_cold))
    }

    fn get_efficiency_from_flowrate(&self, flowrate: f64) -> f64 {
        interp(&self.flow_rates, &self.efficiencies, flowrate)
    }
}

/// A class to represent instantaneous waste water heat recovery systems with arrangement C
///
/// For System C WWHRS, output of the heat exchanger is fed to the hot water system only
#[derive(Clone)]
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
    ) -> Self {
        // assuming the first timestep index is wanted here, but this is unclear!
        let stored_temperature = cold_water_source.temperature(0usize);

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

    pub fn temperature(&self) -> f64 {
        self.stored_temperature
    }

    pub fn return_temperature(
        &self,
        temp_target: f64,
        flowrate_waste_water: f64,
        timestep_idx: usize,
    ) -> f64 {
        // # TODO (from Python) The cold water flow rate depends on the temperature returned from
        // #      this function, which may create a circularity in the calculation.
        // #      May need to integrate System B into shower module and/or combine
        // #      equations.

        let temp_cold = self.cold_water_source.temperature(timestep_idx);

        // # TODO If flowrates have been provided for waste and cold water:
        // #    - Calc heat recovered from waste water. Need to do this per shower
        // #      individually? Need WWHRS_Connection object?
        let wwhrs_efficiency =
            self.get_efficiency_from_flowrate(flowrate_waste_water) * self.utilisation_factor;

        // Calculate temp of pre-heated water based on heat recovered and flow rates
        temp_cold + ((wwhrs_efficiency / 100.0) * (temp_target - temp_cold))
    }

    fn get_efficiency_from_flowrate(&self, flowrate: f64) -> f64 {
        interp(&self.flow_rates, &self.efficiencies, flowrate)
    }
}

/// A class to represent instantaneous waste water heat recovery systems with arrangement A
///     
/// For System A WWHRS, output of the heat exchanger is fed to both the shower
/// and the hot water system
#[derive(Clone)]
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
    ) -> Self {
        // assume that the first timestep is what is wanted here though could be wrong!!
        let stored_temperature = cold_water_source.temperature(0usize);

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

    pub fn temperature(&self) -> f64 {
        self.stored_temperature
    }

    pub fn return_temperature(
        &self,
        temp_target: f64,
        flowrate_waste_water: f64,
        timestep_idx: usize,
    ) -> f64 {
        // # TODO (from Python) The cold water flow rate depends on the temperature returned from
        // #      this function, which may create a circularity in the calculation.
        // #      May need to integrate System B into shower module and/or combine
        // #      equations.

        let temp_cold = self.cold_water_source.temperature(timestep_idx);

        // # TODO If flowrates have been provided for waste and cold water:
        // #    - Calc heat recovered from waste water. Need to do this per shower
        // #      individually? Need WWHRS_Connection object?
        let wwhrs_efficiency =
            self.get_efficiency_from_flowrate(flowrate_waste_water) * self.utilisation_factor;

        // Calculate temp of pre-heated water based on heat recovered and flow rates
        temp_cold + ((wwhrs_efficiency / 100.0) * (temp_target - temp_cold))
    }

    fn get_efficiency_from_flowrate(&self, flowrate: f64) -> f64 {
        interp(&self.flow_rates, &self.efficiencies, flowrate)
    }
}
