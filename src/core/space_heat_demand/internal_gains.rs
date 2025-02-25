use crate::compare_floats::min_of_2;
use crate::core::controls::time_control::SmartApplianceControl;
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::units::WATTS_PER_KILOWATT;
use crate::input::{ApplianceGainsDetails, ApplianceGainsDetailsEvent, ApplianceLoadShifting};
use crate::simulation_time::{SimulationTimeIteration, SimulationTimeIterator};
use anyhow::anyhow;
use atomic_float::AtomicF64;
use parking_lot::RwLock;
use std::sync::atomic::Ordering;
use std::sync::Arc;

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

#[derive(Debug)]
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
    ) -> anyhow::Result<f64> {
        match self {
            Gains::Internal(internal) => Ok(internal.total_internal_gain_in_w(zone_area, simtime)),
            Gains::Appliance(appliance) => appliance.total_internal_gain_in_w(zone_area, simtime),
            Gains::Event(event_appliance_gains) => {
                event_appliance_gains.total_internal_gain_in_w(zone_area, simtime)
            }
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
    pub(crate) fn new(
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
    ) -> anyhow::Result<f64> {
        let total_energy_supplied = self.total_energy_supply
            [simtime.time_series_idx(self.start_day, self.time_series_step)];
        let total_energy_supplied_w = total_energy_supplied * zone_area;
        let total_energy_supplied_kwh =
            total_energy_supplied_w / WATTS_PER_KILOWATT as f64 * simtime.timestep;

        self.energy_supply_connection
            .demand_energy(total_energy_supplied_kwh, simtime.index)?;

        Ok(total_energy_supplied_w * self.gains_fraction)
    }
}

/// An object to represent internal gains and energy consumption from appliances
#[derive(Debug)]
pub struct EventApplianceGains {
    energy_supply_conn: EnergySupplyConnection,
    energy_supply_name: String,
    gains_fraction: f64,
    _start_day: u32,
    time_series_step: f64,
    _series_length: usize,
    load_shifting_metadata: Option<LoadShiftingMetadata>,
    max_shift: f64,
    usage_events: Arc<RwLock<Vec<ApplianceGainsDetailsEvent>>>,
    total_floor_area: f64,
    total_power_supply: Vec<AtomicF64>,
    standby_power: f64,
    smart_control: Option<Arc<SmartApplianceControl>>,
    simulation_timestep_count: usize,
    simulation_timestep: f64,
}

impl EventApplianceGains {
    ///  Arguments:
    ///  * `energy_supply_connection` -- reference to EnergySupplyConnection object representing
    ///                                 the electricity supply attached to the appliance
    ///  * `simulation_time`          -- reference to SimulationTime object
    ///  * `appliance_data`           -- dictionary of appliance gains data from project dict, including:
    ///                                 gains_fraction    -- proportion of appliance demand turned into heat gains
    ///                                 start_day         -- first day of the time series, day of the year, 0 to 365 (single value)
    ///                                 time_series_step  -- timestep of the time series data, in hours
    ///                                 Standby           -- appliance power consumption when not in use in Watts
    ///                                 Events            -- list of appliance usage events, which are dictionaries,
    ///                                                      containing demand_W, start and duration,
    ///                                                      with start and duration in hours
    ///                                 loadshifting      -- (optional) dictionary defining loadshifting parameters
    ///                                                          max_shift_hrs         - the maximum time, in hours, that an event
    ///                                                                                  may be shifted away from when it was originally
    ///                                                                                  intended to occur. This may be up to 24 hours.
    ///                                                          weight_timeseries     - this may be, for example, the hourly cost per
    ///                                                                                  kWh of a 7hr tariff, but could be any time series.
    ///                                                                                  The sum of the other demand and the demand of the
    ///                                                                                  appliance in question at any given time is multiplied
    ///                                                                                  by the value of this timeseries to obtain a figure
    ///                                                                                  that is used to determine whether to shift an event,
    ///                                                                                  and when would be the most appropriate time to shift
    ///                                                                                  the event to.
    ///                                                          demand_limit_weighted - value above which the sum demand multiplied by the
    ///                                                                                  weight should not exceed. If a 7hr tariff were used for
    ///                                                                                  the weight timeseries, then this would be a cost. If this value
    ///                                                                                  is 0, then all events will be shifted to the optimal time
    ///                                                                                  within the window, otherwise they will be shifted to the earliest
    ///                                                                                  time at which the weighted demand goes below the limit.
    ///                                                                                  (if there is no time in the window when the weighted demand is below
    ///                                                                                  the limit, then the optimum is chosen.)
    ///  * `total_floor_area`                      -- total floor area of dwelling
    ///  * `smart_control`            -- (optional) reference to a smart control (required for loadshifting, otherwise not needed)
    pub(crate) fn new(
        energy_supply_conn: EnergySupplyConnection,
        simulation_time: &SimulationTimeIterator,
        appliance_data: &ApplianceGainsDetails,
        total_floor_area: f64,
        smart_control: Option<Arc<SmartApplianceControl>>,
    ) -> anyhow::Result<Self> {
        let standby_power = appliance_data
            .standby
            .ok_or_else(|| anyhow!("standby is expected for EventApplianceGains"))?;
        let usage_events = appliance_data
            .events
            .clone()
            .ok_or_else(|| anyhow!("events are expected for EventApplianceGains"))?;
        let load_shifting_metadata = appliance_data
            .load_shifting
            .as_ref()
            .map(|load_shifting| anyhow::Ok::<LoadShiftingMetadata>(load_shifting.try_into()?))
            .transpose()?;
        let time_series_step = appliance_data.time_series_step;
        let series_length = (simulation_time.total_steps() as f64
            / simulation_time.step_in_hours()
            / time_series_step)
            .ceil() as usize;
        let max_shift = match &appliance_data.load_shifting {
            Some(value) => value.max_shift_hrs / simulation_time.step_in_hours(),
            None => -1.,
        };
        Ok(Self {
            energy_supply_conn,
            energy_supply_name: appliance_data.energy_supply.clone(),
            gains_fraction: appliance_data.gains_fraction,
            _start_day: appliance_data.start_day,
            time_series_step,
            _series_length: series_length,
            load_shifting_metadata,
            max_shift,
            usage_events: Arc::new(RwLock::new(usage_events)),
            total_floor_area,
            total_power_supply: [standby_power]
                .into_iter()
                .cycle()
                .take(simulation_time.total_steps())
                .map(AtomicF64::new)
                .collect(),
            standby_power,
            smart_control: match appliance_data.load_shifting.is_some() {
                true => Some(smart_control.ok_or_else(|| {
                    anyhow!("Smart control is required when appliance uses load shifting.")
                })?),
                false => None,
            },
            simulation_timestep_count: simulation_time.total_steps(),
            simulation_timestep: simulation_time.step_in_hours(),
        })
    }

    fn process_events(&self, simtime: SimulationTimeIteration) -> anyhow::Result<()> {
        // adds demand from events up to and including the current timestep
        // to the total annual demand. If there is loadshifting the demand of the event
        // may occur in the future rather than at the time specified
        let usage_event_count = self.usage_events.read().len();
        for _event_idx in 0..usage_event_count {
            let event = if self.usage_events.read()[0].start <= simtime.current_hour() as f64 {
                self.usage_events
                    .write()
                    .pop()
                    .expect("An event was expected to be available to be popped while processing appliance gains events.")
            } else {
                // no events to process yet
                break;
            };
            let (start_idx, power_timesteps) = self.process_event(&event)?;
            for (i, power) in power_timesteps.iter().enumerate() {
                let t_idx = min_of_2(start_idx + i, self.simulation_timestep_count - 1);
                self.total_power_supply[t_idx].fetch_add(*power, Ordering::SeqCst);
                if self.max_shift >= 0. {
                    // a smart control is expected logically to exist here
                    let smart_control = self.smart_control.as_ref().expect("Smart control is expected to exist when max shift is above zero, meaning load shifting is indicated.");
                    smart_control.add_appliance_demand(
                        simtime,
                        power / WATTS_PER_KILOWATT as f64 * self.simulation_timestep,
                        &self.energy_supply_name,
                    )
                }
            }
        }

        Ok(())
    }

    fn process_event(
        &self,
        event: &ApplianceGainsDetailsEvent,
    ) -> anyhow::Result<(usize, Vec<f64>)> {
        let (start_idx, power_list_over_timesteps) = self.event_to_schedule(event);
        Ok((
            start_idx
                + if self.max_shift >= 0. {
                    self.shift_iterative(start_idx, &power_list_over_timesteps, event)?
                } else {
                    0
                },
            power_list_over_timesteps,
        ))
    }

    //     def __process_event(self, eventdict):
    //         start_idx, power_list_over_timesteps = self.__event_to_schedule(eventdict)
    //         if self.__max_shift >= 0:
    //             start_shift = self.__shift_iterative(start_idx, power_list_over_timesteps, eventdict)
    //             return start_idx + start_shift, power_list_over_timesteps
    //         return start_idx, power_list_over_timesteps

    /// shifts an event forward in time one timestep at a time,
    /// until either the total weighted demand on that timestep is below demandlimit
    /// or the event has been shifted beyond the maximum allowed number of timesteps
    /// away from its original position. In the latter case, move the event to the
    /// most favourable time within the allowed window
    fn shift_iterative(
        &self,
        start_idx: usize,
        power_list_over_timesteps: &[f64],
        event: &ApplianceGainsDetailsEvent,
    ) -> anyhow::Result<usize> {
        // pos list will store the total weighted demand (including from the rest of the dwelling)
        // for the usage event happening at the intended time, or 1 timestep into the future, or 2, up to
        // the max shift time.
        // the lowest value in this list will represent the time at which the usage event would result in the
        // lowest demand.

        let ceil_max_shift = self.max_shift.ceil();
        let mut pos_list = vec![
            0.;
            if ceil_max_shift >= 0. {
                ceil_max_shift as usize
            } else {
                0
            }
        ];
        for (start_shift, pos_list_entry) in pos_list.iter_mut().enumerate() {
            for (i, power) in power_list_over_timesteps.iter().enumerate() {
                let t_idx = min_of_2(
                    start_idx + i + start_shift,
                    self.simulation_timestep_count - 1,
                );
                let series_idx = (t_idx as f64 * self.simulation_timestep / self.time_series_step)
                    .floor() as usize;
                let weight_timeseries = &self
                    .load_shifting_metadata
                    .as_ref()
                    .ok_or_else(|| {
                        anyhow!("Internal gains event processing expects load shifting to be set.")
                    })?
                    .weight_timeseries;

                if self.total_power_supply[t_idx].load(Ordering::SeqCst) >= event.duration {
                    // the appliance is already turned on for the entire timestep
                    // cannot put an event here
                    // put arbitrarily high demand at this position so it is not picked
                    *pos_list_entry += 10_000. * weight_timeseries[series_idx];
                    break;
                }

                let other_demand = self
                    .smart_control
                    .as_ref()
                    .ok_or_else(|| {
                        anyhow!("Internal gains event processing expects load shifting to be set.")
                    })?
                    .get_demand(t_idx, &self.energy_supply_name);
                let new_demand = power / WATTS_PER_KILOWATT as f64 * self.simulation_timestep;
                *pos_list_entry += (new_demand + other_demand) * weight_timeseries[series_idx];
            }
            let demand_limit = self
                .load_shifting_metadata
                .as_ref()
                .ok_or_else(|| {
                    anyhow!("Internal gains event processing expects load shifting to be set.")
                })?
                .demand_limit;
            if demand_limit > 0. && *pos_list_entry < demand_limit {
                // demand is below the limit, good enough, no need to look further into the future
                return Ok(start_shift);
            }
        }

        let min_in_pos_list = *pos_list
            .iter()
            .min_by(|&a, &b| a.total_cmp(b))
            .ok_or_else(|| anyhow!("Insufficient max shift size for appliance gains."))?;
        Ok(pos_list
            .iter()
            .position(|&r| r == min_in_pos_list)
            .expect("Position expected to be findable for minimum value in a list."))
    }

    fn event_to_schedule(&self, usage_event: &ApplianceGainsDetailsEvent) -> (usize, Vec<f64>) {
        let ApplianceGainsDetailsEvent {
            start,
            duration,
            demand_w: demand_w_event,
        } = usage_event;
        let start_offset = start % self.simulation_timestep;
        let start_idx = (start / self.simulation_timestep).floor() as usize;

        // if the event overruns the end of the timestep it starts in,
        // power needs to be allocated to two (or more) timesteps
        // according to the length of time within each timestep the appliance is being used for
        let mut integralx = 0.0;
        let mut power_timesteps = vec![0.; (duration / self.simulation_timestep).ceil() as usize];
        while integralx < *duration {
            let segment = (self.simulation_timestep - start_offset).min(duration - integralx);
            let idx = (integralx / self.simulation_timestep).floor() as usize;
            // subtract standby power from the added event power
            // as it is already accounted for when the list is initialised
            *power_timesteps.get_mut(idx).unwrap() +=
                (demand_w_event - self.standby_power) * segment;
            integralx += segment;
        }

        (start_idx, power_timesteps)
    }

    /// Return the total internal gain for the current timestep, in W
    pub(crate) fn total_internal_gain_in_w(
        &self,
        zone_area: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // first process any usage events expected to occur within this timestep
        self.process_events(simtime)?;
        // Forward electricity demand (in kWh) to relevant EnergySupply object
        let total_power_supplied = self.total_power_supply[simtime.index].load(Ordering::SeqCst);
        let total_power_supplied_zone = total_power_supplied * zone_area / self.total_floor_area;
        let total_energy_supplied_kwh =
            total_power_supplied_zone / WATTS_PER_KILOWATT as f64 * simtime.timestep;

        self.energy_supply_conn
            .demand_energy(total_energy_supplied_kwh, simtime.index)?;

        Ok(total_power_supplied_zone * self.gains_fraction)
    }
}

#[derive(Clone, Debug)]
struct LoadShiftingMetadata {
    weight_timeseries: Vec<f64>,
    demand_limit: f64,
}

impl TryFrom<&ApplianceLoadShifting> for LoadShiftingMetadata {
    type Error = anyhow::Error;

    fn try_from(input: &ApplianceLoadShifting) -> Result<Self, Self::Error> {
        Ok(Self {
            weight_timeseries: input.weight_timeseries.as_ref().ok_or_else(|| anyhow!("Expected a weight timeseries to have been available as part of load shifting data for internal gains."))?.clone(),
            demand_limit: input.demand_limit_weighted,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyBuilder};
    use crate::input::{FuelType, SmartApplianceBattery};
    use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
    use indexmap::IndexMap;
    use parking_lot::RwLock;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use std::sync::Arc;

    #[fixture]
    fn simulation_time_iterator() -> SimulationTimeIterator {
        SimulationTime::new(0.0, 4.0, 1.0).iter()
    }

    #[fixture]
    fn total_internal_gains() -> Vec<f64> {
        vec![3.2, 4.6, 7.3, 5.2]
    }

    #[rstest]
    fn should_have_correct_total_internal_gain(
        total_internal_gains: Vec<f64>,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let internal_gains = InternalGains {
            total_internal_gains,
            start_day: 0,
            time_series_step: 1.0,
        };
        let expected = [32.0, 46.0, 73.0, 52.0];
        for iteration in simulation_time_iterator {
            assert_eq!(
                internal_gains.total_internal_gain_in_w(10.0, iteration,),
                expected[iteration
                    .time_series_idx(internal_gains.start_day, internal_gains.time_series_step)]
            );
        }
    }

    #[rstest]
    fn test_total_internal_gain_for_appliance(simulation_time_iterator: SimulationTimeIterator) {
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(
                FuelType::Electricity,
                simulation_time_iterator.total_steps(),
            )
            .build(),
        ));
        let energy_supply_connection =
            EnergySupply::connection(energy_supply.clone(), "lighting").unwrap();
        let total_energy_supply = vec![32.0, 46.0, 30.0, 20.0];
        let total_internal_gains = [160.0, 230.0, 150.0, 100.0];
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
                appliance_gains
                    .total_internal_gain_in_w(10.0, iteration)
                    .unwrap(),
                total_internal_gains[iteration.index]
            );
            assert_eq!(
                energy_supply.read().results_by_end_user()["lighting"][iteration.index],
                expected_energy_supply_results[iteration.index],
                "incorrect electricity demand returned"
            );
        }
    }

    #[fixture]
    fn simulation_time_for_event_appliance_gains() -> SimulationTime {
        SimulationTime::new(0.0, 24.0, 0.5)
    }

    #[fixture]
    fn total_floor_area() -> f64 {
        100.
    }

    #[fixture]
    fn event_appliance_gains(
        simulation_time_for_event_appliance_gains: SimulationTime,
        total_floor_area: f64,
    ) -> EventApplianceGains {
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(
                FuelType::Electricity,
                simulation_time_for_event_appliance_gains.total_steps(),
            )
            .build(),
        ));
        let energy_supply_connection =
            EnergySupply::connection(energy_supply.clone(), "new_connection").unwrap();
        // following is provided in an unrealistic place in the Python test, so providing
        // power timeseries data separately here as data structure matters in the Rust
        let power_timeseries = vec![
            77.70823134533667,
            70.07710122045972,
            66.26153469022015,
            62.445968159980595,
            58.6304045653432,
            66.26153469022015,
            81.52379787557622,
            872.7293737820189,
            448.95976141145366,
            146.38841421163792,
            150.20398074187744,
            246.13290374339178,
            157.8351108667544,
            146.38841421163792,
            150.20398074187744,
            236.48451946096566,
            1235.1378596577206,
            257.03982010376774,
            801.2313187689566,
            207.4374669530622,
            192.17520376770614,
            290.41445627425867,
            138.757284086761,
            100.60162759117186,
            77.70823134533667,
        ];
        let appliance_data: ApplianceGainsDetails = serde_json::from_value(serde_json::json!({
            "type": "Clothes_drying",
            "EnergySupply": "mains elec",
            "start_day": 0,
            "time_series_step": 1,
            "gains_fraction": 0.7,
            "Events": [{"start": 0.1, "duration": 1.75, "demand_W": 900.0},
                      {"start": 5.3, "duration": 1.50, "demand_W": 900.0}
                      ],
            "Standby": 0.5,
            "loadshifting":
                     {"demand_limit_weighted": 0,
                      "max_shift_hrs": 8,
                      "weight": "Tariff",
                      "weight_timeseries": [
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                        1.0,
                      ]
                    },
        }))
        .unwrap();
        let non_appliance_demand_24hr = IndexMap::from([(
            "mains elec".into(),
            vec![
                0.06830825101566576,
                0.060105811973985415,
                0.05305939286989522,
                0.0501236916686892,
                0.011659722362322093,
                0.009841650327447332,
                0.008628179127672608,
                0.00780480411138585,
                0.007342238555297734,
                0.0068683142767418555,
                0.007092339066083696,
                0.0073738789491060025,
                0.00890318157456338,
                0.013196350717229778,
                3.8185258665440713,
                3.686604856404327,
                3.321741503132428,
                2.1009771847288543,
                1.9577604381160962,
                0.982853996628817,
                -0.4690062980892011,
                -0.47106808673125095,
                -0.440083286258449,
                -0.4402744803667663,
                -0.275893537706527,
                -0.27582574287859735,
                -0.02364598193865506,
                -0.022770656131003677,
                0.006386180325667274,
                1.1838622352604393,
                0.016880847238056107,
                0.022939243503258204,
                0.03451315625040322,
                3.6710272517570663,
                3.4183577199797917,
                3.259661970488036,
                2.382744591886866,
                3.347517299485186,
                2.8633293436146374,
                2.194481295688944,
                2.0245859635409174,
                2.160933115928409,
                2.1559541166936143,
                2.078707457615306,
                0.06082872713634447,
                0.05498437799409048,
                0.04615437212925211,
                0.036699730049032445,
            ],
        )]);
        let battery24hr: SmartApplianceBattery = serde_json::from_value(serde_json::json!({
            "battery_state_of_charge": {
                "mains elec": [
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0
                ]
            }
        }))
        .unwrap();
        let smart_control = SmartApplianceControl::new(
            &IndexMap::from([("mains elec".into(), power_timeseries)]),
            &IndexMap::from([(
                "mains elec".into(),
                appliance_data
                    .load_shifting
                    .as_ref()
                    .unwrap()
                    .clone()
                    .weight_timeseries
                    .unwrap(),
            )]),
            appliance_data.time_series_step,
            &simulation_time_for_event_appliance_gains.iter(),
            Some(non_appliance_demand_24hr),
            Some(&battery24hr),
            &IndexMap::from([("mains elec".into(), energy_supply)]),
            vec!["Clothes_drying".into()],
        )
        .unwrap();

        EventApplianceGains::new(
            energy_supply_connection,
            &simulation_time_for_event_appliance_gains.iter(),
            &appliance_data,
            total_floor_area,
            Some(Arc::new(smart_control)),
        )
        .unwrap()
    }

    #[rstest]
    fn test_process_event(event_appliance_gains: EventApplianceGains) {
        let event = ApplianceGainsDetailsEvent {
            start: 3.,
            duration: 1.75,
            demand_w: 900.0,
        };
        assert_eq!(
            event_appliance_gains.process_event(&event).unwrap(),
            (20, vec![449.75, 449.75, 449.75, 224.875])
        );
    }

    #[rstest]
    fn test_event_to_schedule(event_appliance_gains: EventApplianceGains) {
        let event = ApplianceGainsDetailsEvent {
            start: 3.,
            duration: 1.75,
            demand_w: 900.0,
        };
        assert_eq!(
            event_appliance_gains.event_to_schedule(&event),
            (6, vec![449.75, 449.75, 449.75, 224.875])
        );
    }

    #[rstest]
    fn test_total_internal_gain(
        event_appliance_gains: EventApplianceGains,
        simulation_time_for_event_appliance_gains: SimulationTime,
        total_floor_area: f64,
    ) {
        let res = simulation_time_for_event_appliance_gains
            .iter()
            .map(|simtime| {
                event_appliance_gains.total_internal_gain_in_w(total_floor_area, simtime)
            })
            .collect::<Result<Vec<_>, _>>()
            .unwrap();
        assert_eq!(
            res,
            vec![
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                504.07,
                252.20999999999998,
                252.20999999999998,
                94.79749999999994,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                378.1400000000003,
                252.21000000000018,
                315.17499999999944,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35,
                0.35
            ]
        );
    }

    #[rstest]
    fn test_shift_iterative(event_appliance_gains: EventApplianceGains) {
        let event = ApplianceGainsDetailsEvent {
            start: 2.33,
            duration: 1.0,
            demand_w: 900.,
        };
        let (s, a) = (5usize, [600., 300.]);
        assert_eq!(
            event_appliance_gains
                .shift_iterative(s, &a, &event)
                .unwrap(),
            15
        );
    }
}
