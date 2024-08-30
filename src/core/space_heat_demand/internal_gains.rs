use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::schedule::TypedScheduleEvent;
use crate::core::units::WATTS_PER_KILOWATT;
use crate::input::{ApplianceGainsDetails, ApplianceGainsDetailsEvent};
// use crate::core::units::WATTS_PER_KILOWATT;
use crate::simulation_time::{SimulationTime, SimulationTimeIteration, SimulationTimeIterator};

/// Arguments:
/// * `total_internal_gains` - list of internal gains, in W/m2 (one entry per hour)
/// * `start_day` - first day of time series, day of the year, 0 to 365
/// * `time_series_step` - timestep of the time series data, in hours
#[derive(Clone, Debug)]
pub struct InternalGains {
    total_internal_gains: Vec<f64>,
    start_day: u32,
    time_series_step: f64,
}

#[derive(Clone, Debug)]
pub enum Gains {
    Internal(InternalGains),
    Appliance(ApplianceGains),
    Event(EventApplianceGains),
}

impl Gains {
    pub fn total_internal_gain_in_w(
        &self,
        zone_area: f64,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        match self {
            Gains::Internal(internal) => internal.total_internal_gain_in_w(zone_area, simtime),
            Gains::Appliance(appliance) => appliance.total_internal_gain_in_w(zone_area, simtime),
            _ => todo!(), // TODO handle Gains::EventAppliance case
        }
    }
}

impl InternalGains {
    pub fn new(total_internal_gains: Vec<f64>, start_day: u32, time_series_step: f64) -> Self {
        InternalGains {
            total_internal_gains,
            start_day,
            time_series_step,
        }
    }

    /// Return the total internal gain for the current timestep in W
    pub fn total_internal_gain_in_w(
        &self,
        zone_area: f64,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        self.total_internal_gains[simtime.time_series_idx(self.start_day, self.time_series_step)]
            * zone_area
    }
}

/// Arguments:
/// * `total_energy_supply` - list of energy supply from appliances, in W/m2 (one entry per hour)
/// * `connected_energy_supply` - reference to the energy supply attached to the specific appliance
/// * `end_user_name` - name of the energy supply attached to the specific appliance
/// * `gains_fraction` - fraction of energy supply which is counted as an internal gain
/// * `start_day` - first day of time series, day of the year, 0 to 365
/// * `time_series_step` - timestep of the time series data, in hours
#[derive(Clone, Debug)]
pub struct ApplianceGains {
    total_energy_supply: Vec<f64>,
    energy_supply_connection: EnergySupplyConnection,
    gains_fraction: f64,
    start_day: u32,
    time_series_step: f64,
}

impl ApplianceGains {
    pub fn new(
        total_energy_supply: Vec<f64>,
        gains_fraction: f64,
        start_day: u32,
        time_series_step: f64,
        energy_supply_connection: EnergySupplyConnection,
    ) -> Self {
        Self {
            total_energy_supply,
            gains_fraction,
            start_day,
            time_series_step,
            energy_supply_connection,
        }
    }

    /// Return the total internal gain for the current timestep in W
    pub fn total_internal_gain_in_w(
        &self,
        zone_area: f64,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        let total_energy_supplied = self.total_energy_supply
            [simtime.time_series_idx(self.start_day, self.time_series_step)];
        let total_energy_supplied_w = total_energy_supplied * zone_area;
        let total_energy_supplied_kwh =
            total_energy_supplied_w / WATTS_PER_KILOWATT as f64 * self.time_series_step;

        self.energy_supply_connection
            .demand_energy(total_energy_supplied_kwh, simtime.index)
            .unwrap();

        total_energy_supplied_w * self.gains_fraction
    }
}

/// An object to represent internal gains and energy consumption from appliances
#[derive(Clone, Debug)]
pub struct EventApplianceGains {
    energy_supply_conn: EnergySupplyConnection,
    gains_fraction: f64,
    start_day: u32,
    time_series_step: f64,
    total_floor_area: f64,
    standby_power: f64,
    usage_events: Vec<ApplianceGainsDetailsEvent>,
    max_shift: f64,
    demand_limit: Option<f64>,
    weight_timeseries: Option<Vec<f64>>,
    otherdemand_timeseries: Option<Vec<f64>>,
    total_power_supply: Vec<f64>,
}

impl EventApplianceGains {
    // Construct a InternalGains object

    //     Arguments:
    //     energy_supply_connection -- reference to EnergySupplyConnection object representing
    //                                 the electricity supply attached to the appliance
    //     simulation_time          -- reference to SimulationTime object
    //     appliance_data           -- dictionary of appliance gains data from project dict, including:
    //                                 gains_fraction    -- proportion of appliance demand turned into heat gains
    //                                 start_day         -- first day of the time series, day of the year, 0 to 365 (single value)
    //                                 time_series_step  -- timestep of the time series data, in hours
    //                                 Standby           -- appliance power consumption when not in use in Watts
    //                                 Events            -- list of appliance usage events, which are dictionaries,
    //                                                      containing demand_W, start and duration,
    //                                                      with start and duration in hours
    //                                 loadshifting      -- (optional) dictionary definining loadshifting parameters
    //                                                      max_shift_hrs
    //                                                      demand_limit_weighted
    //                                                      weight_timeseries
    //                                                      demand_timeseries
    //     TFA                      -- total floor area of dwelling
    pub fn new(
        energy_supply_conn: EnergySupplyConnection,
        simulation_time: &SimulationTimeIterator,
        appliance_data: &ApplianceGainsDetails,
        total_floor_area: f64,
    ) -> Self {
        let standby_power = appliance_data
            .standby
            .expect("standby is expected for EventApplianceGains");
        let usage_events = appliance_data
            .events
            .clone()
            .expect("events are expected for EventApplianceGains");
        let weight_timeseries = match &appliance_data.load_shifting {
            Some(value) => value.weight_timeseries.clone(),
            None => None,
        };
        let demand_limit = match &appliance_data.load_shifting {
            Some(value) => Some(value.demand_limit_weighted),
            None => None,
        };
        let otherdemand_timeseries = match &appliance_data.load_shifting {
            Some(value) => value.demand_timeseries.clone(),
            None => None,
        };
        let time_series_step = appliance_data.time_series_step;
        // TODO remove anything from self that is only used in total_power_supply
        // TODO can any of the clone() calls be passed by reference instead
        Self {
            energy_supply_conn: energy_supply_conn,
            gains_fraction: appliance_data.gains_fraction,
            start_day: appliance_data.start_day,
            time_series_step,
            total_floor_area: total_floor_area,
            standby_power,
            usage_events: usage_events.clone(),
            max_shift: match &appliance_data.load_shifting {
                Some(value) => value.max_shift_hrs / appliance_data.time_series_step,
                None => 0.,
            },
            demand_limit,
            weight_timeseries: weight_timeseries.clone(),
            otherdemand_timeseries: otherdemand_timeseries.clone(),
            total_power_supply: Self::total_power_supply(
                simulation_time,
                time_series_step,
                standby_power,
                usage_events,
                otherdemand_timeseries,
                weight_timeseries,
            ),
        }
    }

    fn total_power_supply(
        simulation_time: &SimulationTimeIterator,
        time_series_step: f64,
        standby_power: f64,
        usage_events: Vec<ApplianceGainsDetailsEvent>,
        otherdemand_timeseries: Option<Vec<f64>>,
        weight_timeseries: Option<Vec<f64>>,
    ) -> Vec<f64> {
        // initialize list with standby power on all timesteps
        let length = (simulation_time.total_steps() as f64 * simulation_time.step_in_hours()
            / time_series_step)
            .ceil();
        let mut total_power_supply = vec![standby_power; length as usize];

        // focus on 2 factors - per appliance shiftability
        // and overall shifting
        for event in usage_events {
            let (s, a) = Self::shift_event(
                event,
                otherdemand_timeseries.clone(),
                weight_timeseries.clone(),
            );
            for (i, x) in a.iter().enumerate() {
                let index = (s + i as f64).floor() as usize % total_power_supply.len();
                total_power_supply[index] += *x as f64;
            }
        }
        total_power_supply
    }

    fn shift_event(
        usage_events: ApplianceGainsDetailsEvent,
        otherdemand_timeseries: Option<Vec<f64>>,
        weight_timeseries: Option<Vec<f64>>,
    ) -> (f64, Vec<u32>) {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::energy_supply::energy_supply::EnergySupply;
    use crate::input::FuelType;
    use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
    use parking_lot::RwLock;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use std::sync::Arc;

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
            start_day: 0,
            time_series_step: 1.0,
        };
        let expected = vec![32.0, 46.0, 73.0, 52.0];
        for iteration in simulation_time_iterator {
            assert_eq!(
                internal_gains.total_internal_gain_in_w(10.0, iteration,),
                expected[iteration
                    .time_series_idx(internal_gains.start_day, internal_gains.time_series_step)]
            );
        }
    }

    #[rstest]
    pub fn test_total_internal_gain_for_appliance(
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let energy_supply = Arc::new(RwLock::new(EnergySupply::new(
            FuelType::Electricity,
            simulation_time_iterator.total_steps(),
            None,
            None,
            None,
        )));
        let energy_supply_connection =
            EnergySupply::connection(energy_supply.clone(), "lighting").unwrap();
        let total_energy_supply = vec![32.0, 46.0, 30.0, 20.0];
        let total_internal_gains = vec![160.0, 230.0, 150.0, 100.0];
        let expected_energy_supply_results = [0.32, 0.46, 0.30, 0.20];
        let appliance_gains = ApplianceGains {
            total_energy_supply,
            energy_supply_connection,
            gains_fraction: 0.5,
            start_day: 0,
            time_series_step: 1.0,
        };
        for iteration in simulation_time_iterator.clone() {
            assert_eq!(
                appliance_gains.total_internal_gain_in_w(10.0, iteration),
                total_internal_gains[iteration.index]
            );
            assert_eq!(
                energy_supply.read().results_by_end_user()["lighting"][iteration.index],
                expected_energy_supply_results[iteration.index],
                "incorrect electricity demand returned"
            );
        }
    }
}
