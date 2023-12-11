use crate::core::energy_supply::energy_supply::EnergySupply;
// use crate::core::units::WATTS_PER_KILOWATT;
use crate::simulation_time::SimulationTimeIterator;

/// Arguments:
/// * `total_internal_gains` - list of internal gains, in W/m2 (one entry per hour)
/// * `simulation_time_iterator`
/// * `start_day` - first day of time series, day of the year, 0 to 365
/// * `time_series_step` - timestep of the time series data, in hours
pub struct InternalGains<'a> {
    total_internal_gains: Vec<f64>,
    simulation_time_iterator: &'a SimulationTimeIterator,
    start_day: u32,
    time_series_step: f64,
}

impl<'a> InternalGains<'a> {
    /// Return the total internal gain for the current timestep in W
    pub fn total_internal_gain_in_w(&self, zone_area: f64, timestep_idx: usize) -> f64 {
        self.total_internal_gains[timestep_idx] * zone_area
    }
}

/// Arguments:
/// * `total_energy_supply` - list of energy supply from appliances, in W/m2 (one entry per hour)
/// * `connected_energy_supply` - reference to the energy supply attached to the specific appliance
/// * `end_user_name` - name of the energy supply attached to the specific appliance
/// * `gains_fraction` - fraction of energy supply which is counted as an internal gain
/// * `simulation_time_iterator`
/// * `start_day` - first day of time series, day of the year, 0 to 365
/// * `time_series_step` - timestep of the time series data, in hours
pub struct ApplianceGains<'a> {
    total_energy_supply: Vec<f64>,
    connected_energy_supply: &'a EnergySupply,
    end_user_name: &'a str,
    gains_fraction: f64,
    simulation_time_iterator: &'a SimulationTimeIterator,
    start_day: u32,
    time_series_step: f64,
}

impl<'a> ApplianceGains<'a> {
    /// Return the total internal gain for the current timestep in W
    pub fn total_internal_gain_in_w(&self, zone_area: f64, timestep_idx: usize) -> f64 {
        let total_energy_supplied = self.total_energy_supply[timestep_idx];
        let total_energy_supplied_w = total_energy_supplied * zone_area;
        // let total_energy_supplied_kWh =
        //     total_energy_supplied_w / WATTS_PER_KILOWATT as f64 * self.time_series_step;

        // TODO: in Python, this tries to mutate the energy supply object, which is weird
        // need to work out what is happening here and come up with a workaround
        //self.connected_energy_supply.demand_energy(self.end_user_name, total_energy_supplied_kWh);

        total_energy_supplied_w * self.gains_fraction
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::EnergySupplyType;
    use crate::simulation_time::SimulationTime;
    use rstest::*;

    #[fixture]
    pub fn simulation_time_iterator() -> SimulationTimeIterator {
        SimulationTime::new(0.0, 4.0, 1.0).iter()
    }

    #[fixture]
    pub fn total_internal_gains() -> Vec<f64> {
        vec![3.2, 4.6, 7.3, 5.2]
    }

    #[rstest]
    pub fn should_have_correct_total_internal_gain(
        total_internal_gains: Vec<f64>,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let internal_gains = InternalGains {
            total_internal_gains,
            simulation_time_iterator: &(simulation_time_iterator.clone()),
            start_day: 0,
            time_series_step: 1.0,
        };
        let expected = vec![32.0, 46.0, 73.0, 52.0];
        for iteration in simulation_time_iterator {
            assert_eq!(
                internal_gains.total_internal_gain_in_w(
                    10.0,
                    iteration
                        .time_series_idx(internal_gains.start_day, internal_gains.time_series_step)
                ),
                expected[iteration
                    .time_series_idx(internal_gains.start_day, internal_gains.time_series_step)]
            );
        }
    }

    #[rstest]
    pub fn test_total_internal_gain_for_appliance(
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let energy_supply = EnergySupply::new(
            EnergySupplyType::Electricity,
            simulation_time_iterator.clone(),
            None,
        );
        let total_energy_supply = vec![32.0, 46.0, 30.0, 20.0];
        let total_internal_gains = vec![160.0, 230.0, 150.0, 100.0];
        let appliance_gains = ApplianceGains {
            total_energy_supply,
            connected_energy_supply: &energy_supply,
            end_user_name: "lighting",
            gains_fraction: 0.5,
            simulation_time_iterator: &simulation_time_iterator,
            start_day: 0,
            time_series_step: 1.0,
        };
        for iteration in simulation_time_iterator.clone() {
            assert_eq!(
                appliance_gains.total_internal_gain_in_w(10.0, iteration.index),
                total_internal_gains[iteration.index]
            );
            // TODO: Python code has a test here for energy supply results by end user
        }
    }
}
