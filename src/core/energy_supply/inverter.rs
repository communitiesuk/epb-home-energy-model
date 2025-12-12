use crate::compare_floats::min_of_2;
use crate::compare_floats::min_of_3;
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::hem_core::simulation_time::SimulationTimeIterator;
use crate::input::InverterType;

/// An object to represent an inverter.  Although primarily for PV systems this may have other uses e.g. wind turbine
#[derive(Debug, Clone)]
pub(crate) struct Inverter {
    energy_supply_connection: EnergySupplyConnection,
    simulation_time: SimulationTimeIterator,
    inverter_peak_power_dc: f64,
    inverter_peak_power_ac: f64,
    inverter_is_inside: bool,
    inverter_type: InverterType,
}
impl Inverter {
    /// Construct an Inverter object
    ///
    /// Arguments:
    /// * `energy_supply_conn`     -- reference to EnergySupplyConnection object
    /// * `simulation_time`        -- reference to SimulationTime object
    /// * `inverter_peak_power_dc` -- Peak power in kW; represents the peak electrical DC power input to the inverter
    /// * `inverter_is_inside`     -- tells us that the inverter is considered inside the building
    /// * `inverter_peak_power_ac` -- Peak power in kW; represents the peak electrical AC power output from the inverter
    /// * `inverter_type`          -- type of inverter to help with calculation of efficiency of inverter when overshading
    pub(crate) fn new(
        energy_supply_connection: EnergySupplyConnection,
        simulation_time: SimulationTimeIterator,
        inverter_peak_power_dc: f64,
        inverter_peak_power_ac: f64,
        inverter_is_inside: bool,
        inverter_type: InverterType,
    ) -> Self {
        Self {
            energy_supply_connection,
            simulation_time,
            inverter_peak_power_dc,
            inverter_peak_power_ac,
            inverter_is_inside,
            inverter_type,
        }
    }

    /// Return the inverter type
    pub(crate) fn r#type(&self) -> InverterType {
        self.inverter_type
    }

    /// Return whether this unit is considered inside the building or not
    pub(crate) fn is_inside(&self) -> bool {
        self.inverter_is_inside
    }

    /// Returns the efficiency of the inverter calculated based on the ratio of input power to the dc capacity of the inverter
    fn efficiency_power_ratio(&self, power_input: f64) -> f64 {
        // TODO (from Python) this should be included in the efficiency_lookup method as an inverter type
        // Calculate the input as a ratio of the dc capacity of the inverter
        // input limited to be within the peak power dc capacity

        let ratio_of_rated_output =
            min_of_2(power_input, self.inverter_peak_power_dc) / self.inverter_peak_power_dc;

        // Using Ratio of Rated Power, calculate Inverter DC to AC efficiency
        // equation was estimated based on graph from
        // https://www.researchgate.net/publication/260286647_Performance_of_PV_inverters figure 9
        
        if ratio_of_rated_output == 0. {
            0.
        } else {
            // Empirical efficiency curve fit primarily based on SMA Sunny Boy inverters (largest market share 2018)
            // assisted with Sungrow and Huawei inverters (largest market share 2019).
            // System of 3 equations to fit efficiency curve for Sunny Boy PV2AC Inverters
            let inverter_dc_ac_efficiency_1 =
                97.2 * (1. - (0.18 / (1. + std::f64::consts::E.powf(21. * ratio_of_rated_output))));
            let inverter_dc_ac_efficiency_2 =
                0.5 * (std::f64::consts::PI * ratio_of_rated_output).cos() + 96.9;
            let inverter_dc_ac_efficiency_3 = 97.2 * (30. * ratio_of_rated_output).tanh();
            min_of_3(
                inverter_dc_ac_efficiency_1,
                inverter_dc_ac_efficiency_2,
                inverter_dc_ac_efficiency_3,
            ) / 100.
        }
    }

    /// Calculate energy from input power, applying efficiency and returning the energy produced
    fn calculate_energy_output(&self, power_input: f64) -> f64 {
        // TODO (from Python) this should use the inverter type to determine the efficiency calculation to be used
        // For now hard coded as the SMA sunny boy for consistency with previous PhotoVoltaicSystem
        // Calculate useful power from ac and dc peak capacity and efficiency
        let mut power = min_of_2(power_input, self.inverter_peak_power_dc);
        power *= self.efficiency_power_ratio(power);
        power = min_of_2(power, self.inverter_peak_power_ac);

        // Convert power to energy
        
        power * self.simulation_time.step_in_hours()
    }

    /// Calculate energy from input power, apply efficiency, and supply to energy supply connection
    pub(crate) fn produce_energy(
        &self,
        power_input: f64,
        timestep_index: usize,
    ) -> anyhow::Result<f64> {
        let energy_produced = self.calculate_energy_output(power_input);
        self.energy_supply_connection
            .supply_energy(energy_produced, timestep_index)?;
        Ok(energy_produced)
    }
}
