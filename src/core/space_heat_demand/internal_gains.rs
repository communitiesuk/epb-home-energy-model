use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::units::WATTS_PER_KILOWATT;
use crate::input::{ApplianceGainsDetails, ApplianceGainsDetailsEvent};
use crate::simulation_time::{SimulationTimeIteration, SimulationTimeIterator};
use itertools::Itertools;
use ordered_float::OrderedFloat;

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
    load_shifting_metadata: Option<LoadShiftingMetadata>,
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
        let load_shifting_metadata =
            appliance_data
                .load_shifting
                .as_ref()
                .map(|load_shifting| LoadShiftingMetadata {
                    weight_timeseries: load_shifting
                        .weight_timeseries
                        .as_ref()
                        .cloned()
                        .unwrap_or(vec![]),
                    demand_limit: load_shifting.demand_limit_weighted,
                    otherdemand_timeseries: load_shifting
                        .demand_timeseries
                        .as_ref()
                        .cloned()
                        .unwrap_or(vec![]),
                });
        let time_series_step = appliance_data.time_series_step;
        // TODO remove anything from self that is only used in total_power_supply
        // TODO can any of the clone() calls be passed by reference instead
        let max_shift = match &appliance_data.load_shifting {
            Some(value) => value.max_shift_hrs / appliance_data.time_series_step,
            None => 0.,
        };
        Self {
            energy_supply_conn: energy_supply_conn,
            gains_fraction: appliance_data.gains_fraction,
            start_day: appliance_data.start_day,
            time_series_step,
            total_floor_area: total_floor_area,
            standby_power,
            usage_events: usage_events.clone(),
            max_shift,
            total_power_supply: Self::total_power_supply(
                simulation_time,
                time_series_step,
                standby_power,
                usage_events,
                load_shifting_metadata.as_ref(),
                max_shift,
            ),
            load_shifting_metadata,
        }
    }

    fn total_power_supply(
        simulation_time: &SimulationTimeIterator,
        time_series_step: f64,
        standby_power: f64,
        usage_events: Vec<ApplianceGainsDetailsEvent>,
        load_shifting_metadata: Option<&LoadShiftingMetadata>,
        max_shift: f64,
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
                load_shifting_metadata,
                standby_power,
                time_series_step,
                max_shift,
            );
            for (i, x) in a.iter().enumerate() {
                let index = (s + i as f64).floor() as usize % total_power_supply.len();
                total_power_supply[index] += *x as f64;
            }
        }
        total_power_supply
    }

    fn shift_event(
        usage_event: ApplianceGainsDetailsEvent,
        load_shifting_metadata: Option<&LoadShiftingMetadata>,
        standby_power: f64,
        time_series_step: f64,
        max_shift: f64,
    ) -> (f64, Vec<f64>) {
        // demand limit could also use ie a linear function instead of a hard limit...
        let (s, a) = Self::event_to_schedule(usage_event, standby_power, time_series_step);
        let start_shift = match load_shifting_metadata {
            Some(load_shifting_metadata) if max_shift > 0. => Self::shift_recursive(
                s,
                &a,
                &load_shifting_metadata.otherdemand_timeseries,
                &load_shifting_metadata.weight_timeseries,
                load_shifting_metadata.demand_limit,
                max_shift,
                &mut vec![],
                0,
            ),
            _ => 0,
        };
        (s + start_shift as f64, a)
    }

    /// shifts an event forward in time one timestep at a time,
    /// until either the total weighted demand on that timestep is below demandlimit
    /// or the event has been shifted beyond the maximum allowed number of timesteps
    /// away from its original position. In the latter case, move the event to the
    /// most favourable time within the allowed window
    fn shift_recursive(
        s: f64,
        a: &[f64],
        demand_timeseries: &[f64],
        weight_timeseries: &[f64],
        demandlimit: f64,
        max_shift: f64,
        pos_list: &mut Vec<f64>,
        mut start_shift: usize,
    ) -> usize {
        pos_list.push(0.);
        for (i, x) in a.iter().enumerate() {
            let idx = ((s + i as f64).floor() as usize + start_shift) % demand_timeseries.len();
            let otherdemand = demand_timeseries[idx] * weight_timeseries[idx];
            let newdemand = x * weight_timeseries[idx];
            *pos_list.last_mut().unwrap() += newdemand + otherdemand;
            if newdemand + otherdemand > demandlimit {
                // check if start shift is too high? and if its past limit look up results of
                // each prev shift and choose the best one
                start_shift += 1;
                if start_shift as f64 <= max_shift {
                    start_shift = Self::shift_recursive(
                        s,
                        a,
                        demand_timeseries,
                        weight_timeseries,
                        demandlimit,
                        max_shift,
                        pos_list,
                        start_shift,
                    )
                } else {
                    // choose the timestep within the allowed window with the lowest demand.
                    // also add the entire length of the series - this will be removed again by the modulo operator,
                    // but prevents an infinite loop from occuring
                    let min_index = pos_list
                        .iter()
                        .map(|x| OrderedFloat(*x))
                        .position_min()
                        .unwrap();
                    start_shift = min_index + demand_timeseries.len();
                }
            }
        }
        start_shift
    }

    fn event_to_schedule(
        usage_event: ApplianceGainsDetailsEvent,
        standby_power: f64,
        time_series_step: f64,
    ) -> (f64, Vec<f64>) {
        let ApplianceGainsDetailsEvent {
            start,
            duration,
            demand_w: demand_w_event,
        } = usage_event;
        let start_offset = start % time_series_step;

        // if the event overruns the end of the timestep it starts in,
        // power needs to be allocated to two (or more) timesteps
        // according to the length of time within each timestep the appliance is being used for
        let mut integralx = 0.0;
        let mut res = vec![0.; (duration / time_series_step).ceil() as usize];
        while integralx < duration {
            let segment = (time_series_step - start_offset).min(duration - integralx);
            let idx = (integralx / time_series_step).floor() as usize;
            // subtract standby power from the added event power
            // as it is already accounted for when the list is initialised
            *res.get_mut(idx).unwrap() += (demand_w_event - standby_power) * segment;
            integralx += segment;
        }

        (start, res)
    }
}

#[derive(Clone, Debug)]
struct LoadShiftingMetadata {
    weight_timeseries: Vec<f64>,
    demand_limit: f64,
    otherdemand_timeseries: Vec<f64>,
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
