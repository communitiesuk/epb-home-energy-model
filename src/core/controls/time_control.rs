// This module provides structs to model time controls

use crate::core::units::{HOURS_PER_DAY, WATTS_PER_KILOWATT};
use crate::external_conditions::ExternalConditions;
use crate::input::{
    ControlCombination, ControlCombinationOperation, ControlCombinations, ControlLogicType,
    ExternalSensor, ExternalSensorCorrelation, HeatSourceControlType, SetpointBoundsInput,
    SmartApplianceBattery, MAIN_REFERENCE,
};
use crate::simulation_time::{SimulationTimeIteration, SimulationTimeIterator, HOURS_IN_DAY};
use anyhow::{anyhow, bail};
use atomic_float::AtomicF64;
use bounded_vec_deque::BoundedVecDeque;
use indexmap::IndexMap;
use itertools::Itertools;
use parking_lot::RwLock;
use smartstring::alias::String;
use std::collections::VecDeque;
use std::fmt::{Debug, Formatter};
use std::iter::repeat;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

#[derive(Debug)]
pub(crate) enum Control {
    OnOffTime(OnOffTimeControl),
    Charge(ChargeControl),
    OnOffMinimisingTime(OnOffCostMinimisingTimeControl),
    SetpointTime(SetpointTimeControl),
    CombinationTime(CombinationTimeControl),
    #[cfg(test)]
    Mock(MockControl),
}

// macro so accessing individual controls through the enum isn't so repetitive
macro_rules! per_control {
    ($val:expr, $pattern:pat => { $res:expr }) => {
        match $val {
            #[allow(noop_method_call)]
            Control::OnOffTime($pattern) => $res,
            #[allow(noop_method_call)]
            Control::Charge($pattern) => $res,
            #[allow(noop_method_call)]
            Control::OnOffMinimisingTime($pattern) => $res,
            #[allow(noop_method_call)]
            Control::SetpointTime($pattern) => $res,
            #[allow(noop_method_call)]
            Control::CombinationTime($pattern) => $res,
            #[cfg(test)]
            #[allow(noop_method_call)]
            Control::Mock($pattern) => $res,
        }
    };
}

use crate::compare_floats::{max_of_2, min_of_2};
use crate::core::energy_supply::energy_supply::EnergySupply;
pub(crate) use per_control;

pub(crate) trait ControlBehaviour: Send + Sync {
    fn in_required_period(
        &self,
        _simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool> {
        None
    }

    fn setpnt(&self, _simulation_time_iteration: &SimulationTimeIteration) -> Option<f64> {
        None
    }

    fn is_on(&self, _simulation_time_iteration: &SimulationTimeIteration) -> bool {
        true
    }
}

impl Debug for dyn ControlBehaviour {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        // if we can downcast self to e.g. Control (if it is one), which we know is Debug, this would be better
        write!(f, "A control object")
    }
}

impl ControlBehaviour for Control {
    fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool> {
        per_control!(self, c => { c.in_required_period(simulation_time_iteration) })
    }

    fn setpnt(&self, simulation_time_iteration: &SimulationTimeIteration) -> Option<f64> {
        per_control!(self, c => { c.setpnt(simulation_time_iteration) })
    }

    fn is_on(&self, simulation_time_iteration: &SimulationTimeIteration) -> bool {
        per_control!(self, c => {c.is_on(&simulation_time_iteration)})
    }
}

#[derive(Debug)]
pub(crate) enum HeatSourceControl {
    HotWaterTimer(Arc<Control>),
    WindowOpening(Arc<Control>),
}

impl HeatSourceControl {
    pub(crate) fn has_type(&self, control_type: HeatSourceControlType) -> bool {
        match control_type {
            HeatSourceControlType::HotWaterTimer => {
                matches!(self, HeatSourceControl::HotWaterTimer(_))
            }
            HeatSourceControlType::WindowOpening => {
                matches!(self, HeatSourceControl::WindowOpening(_))
            }
        }
    }

    pub(crate) fn get(&self) -> Arc<Control> {
        match self {
            HeatSourceControl::HotWaterTimer(control) => control.clone(),
            HeatSourceControl::WindowOpening(control) => control.clone(),
        }
    }
}

/// An object to model a time-only control with on/off (not modulating) operation
#[derive(Clone, Debug)]
pub(crate) struct OnOffTimeControl {
    /// list of boolean values where true means "on" (one entry per hour)
    schedule: Vec<Option<bool>>,
    start_day: u32,
    time_series_step: f64,
}

impl OnOffTimeControl {
    pub(crate) fn new(
        schedule: Vec<Option<bool>>,
        start_day: u32,
        time_series_step: f64,
    ) -> OnOffTimeControl {
        OnOffTimeControl {
            schedule,
            start_day,
            time_series_step,
        }
    }
}

impl ControlBehaviour for OnOffTimeControl {
    fn is_on(&self, timestep: &SimulationTimeIteration) -> bool {
        self.schedule[timestep.time_series_idx(self.start_day, self.time_series_step)]
            .unwrap_or(false)
    }
}

/// An object to model a control that governs electrical charging of a heat storage device
/// that can respond to signals from the grid, for example when carbon intensity is low
#[derive(Debug)]
pub(crate) struct ChargeControl {
    logic_type: ControlLogicType,
    schedule: Vec<bool>,
    start_day: u32,
    time_series_step: f64,
    charge_level: Vec<Option<f64>>,
    temp_charge_cut: Option<f64>,
    temp_charge_cut_delta: Option<Vec<f64>>,
    external_conditions: Option<Arc<ExternalConditions>>,
    external_sensor: Option<ExternalSensor>,
    heat_retention_data: Option<ChargeControlHeatRetentionFields>,
    charge_calc_time: f64,
}

#[derive(Debug)]
pub(crate) struct ChargeControlHeatRetentionFields {
    steps_day: usize,
    demand: Arc<RwLock<BoundedVecDeque<Option<f64>>>>,
    past_ext_temp: Arc<RwLock<BoundedVecDeque<Option<f64>>>>,
    future_ext_temp: Arc<RwLock<BoundedVecDeque<Option<f64>>>>,
    energy_to_store: AtomicF64,
}

impl ChargeControl {
    /// Construct a ChargeControl object
    /// Arguments:
    /// * `logic_type`              - ControlLogicType enum
    /// * `schedule`                - list of boolean values where true means "on" (one entry per hour)
    /// * `simulation_time`         - reference to SimulationTime object
    /// * `start_day`               - first day of the time series, day of the year, 0 to 365 (single value)
    /// * `time_series_step`        - timestep of the time series data, in hours__get_heat_cool_systems_for_zone
    /// * `charge_level`            - Proportion of the charge targeted for each day
    /// * `temp_charge_cut`         - Room temperature at which, if sensed during a charging hour, the control stops charging
    ///                             (Required for AUTOMATIC, CELECT, and HHRSH logic types only)
    /// * `temp_charge_cut_delta`   - array with values for a temperature adjustment which is applied
    ///                                  to the nominal internal air temperature above which the control stops
    ///                                  charging the device with heat.
    ///                                  (Optional for AUTOMATIC, CELECT, and HHRSH logic types)
    /// * `extcond`                 - reference to ExternalConditions object (for HHRSH and HEAT_BATTERY logic)
    /// * `external_sensor`         - external weather sensor that acts as a limiting device to prevent storage
    ///                                  heaters from overcharging (for AUTOMATIC and CELECT logic)
    /// * `charge_calc_time`        - Indicates from which hour of the day the system starts to target the charge level
    ///                             for the next day rather than the current day
    pub(crate) fn new(
        logic_type: ControlLogicType,
        schedule: Vec<bool>,
        simulation_time_iteration: &SimulationTimeIteration,
        start_day: u32,
        time_series_step: f64,
        charge_level: Vec<Option<f64>>,
        temp_charge_cut: Option<f64>,
        temp_charge_cut_delta: Option<Vec<f64>>,
        external_conditions: Option<Arc<ExternalConditions>>,
        external_sensor: Option<ExternalSensor>,
        charge_calc_time: Option<f64>,
    ) -> anyhow::Result<Self> {
        let simulation_timestep = simulation_time_iteration.timestep;
        let charge_calc_time = charge_calc_time.unwrap_or(21.);

        let heat_retention_data: Option<ChargeControlHeatRetentionFields> = match logic_type {
            ControlLogicType::Manual => {
                // do nothing
                None
            }
            ControlLogicType::Automatic => {
                if temp_charge_cut.is_none() {
                    bail!("Automatic ChargeControl definition is missing input parameters.");
                }
                None
            }
            ControlLogicType::Celect => {
                if temp_charge_cut.is_none() {
                    bail!("Celect ChargeControl definition is missing input parameters.");
                }
                None
            }
            ControlLogicType::Hhrsh => {
                if temp_charge_cut.is_none() {
                    bail!("Hhrsh ChargeControl definition is missing input temp_charge_cut parameters.");
                }

                if external_conditions.is_none() {
                    bail!("Hhrsh ChargeControl definition is missing external conditions.");
                }
                let external_conditions = external_conditions.as_ref().unwrap(); // we know it exists at this point

                // Initialize HHRSH-specific attributes
                let steps_day = (HOURS_PER_DAY as f64 / simulation_timestep) as usize;
                let demand = Arc::new(RwLock::new(BoundedVecDeque::from_iter(
                    repeat(None),
                    steps_day,
                )));
                let past_ext_temp = Arc::new(RwLock::new(BoundedVecDeque::from_iter(
                    repeat(None),
                    steps_day,
                )));
                let future_ext_temp = Arc::new(RwLock::new(BoundedVecDeque::from_iter(
                    repeat(Some(0.0)),
                    steps_day,
                )));
                for i in 0..steps_day {
                    future_ext_temp.write().push_back(Some(
                        external_conditions.air_temp_with_offset(simulation_time_iteration, i),
                    ));
                }
                let energy_to_store = AtomicF64::new(0.0);
                Some(ChargeControlHeatRetentionFields {
                    steps_day,
                    demand,
                    past_ext_temp,
                    future_ext_temp,
                    energy_to_store,
                })

                // TODO (from Python) Consider adding solar data for HHRSH logic in addition to heating degree hours.
            }
            ControlLogicType::HeatBattery => {
                // Heat battery doesn't require temp_charge_cut but needs other parameters
                if external_conditions.is_none() {
                    bail!("Heat_battery ChargeControl definition is missing external conditions.");
                }
                let external_conditions = external_conditions.as_ref().unwrap(); // we know it exists at this point

                // Initialize heat battery-specific attributes
                let steps_day = (HOURS_PER_DAY as f64 / simulation_timestep) as usize;
                let demand = Arc::new(RwLock::new(BoundedVecDeque::from_iter(
                    repeat(None),
                    steps_day,
                )));
                let past_ext_temp = Arc::new(RwLock::new(BoundedVecDeque::from_iter(
                    repeat(None),
                    steps_day,
                )));
                let future_ext_temp = Arc::new(RwLock::new(BoundedVecDeque::from_iter(
                    repeat(Some(0.0)),
                    steps_day,
                )));
                for i in 0..steps_day {
                    future_ext_temp.write().push_back(Some(
                        external_conditions.air_temp_with_offset(simulation_time_iteration, i),
                    ));
                }
                let energy_to_store = AtomicF64::new(0.0);
                Some(ChargeControlHeatRetentionFields {
                    steps_day,
                    demand,
                    past_ext_temp,
                    future_ext_temp,
                    energy_to_store,
                })
            }
        };

        Ok(Self {
            logic_type,
            schedule,         // TODO (migration 1.0.0a1) this seems to be optional in Python
            start_day,        // TODO (migration 1.0.0a1) this seems to be optional in Python
            time_series_step, // TODO (migration 1.0.0a1) this seems to be optional in Python
            charge_level,
            temp_charge_cut,
            temp_charge_cut_delta,
            external_conditions,
            external_sensor,
            heat_retention_data,
            charge_calc_time,
        })
    }

    pub(crate) fn logic_type(&self) -> ControlLogicType {
        self.logic_type
    }

    // In Python there is an abstract method target_charge on the ControlCharge class and a concrete implementation here
    /// Return the charge level value from the list given in inputs; one value per day
    pub(crate) fn target_charge(
        &self,
        simtime: SimulationTimeIteration,
        temp_air: Option<f64>,
    ) -> anyhow::Result<f64> {
        // Calculate target charge nominal when unit is on
        let mut target_charge_nominal = if self.is_on(&simtime) {
            self.charge_level
                [simtime.time_series_idx_days(self.start_day, Some(self.charge_calc_time))]
            .unwrap_or_default()
        } else {
            // If unit is off send 0.0 for target charge
            0.
        };

        let target_charge = match self.logic_type {
            ControlLogicType::Manual => target_charge_nominal,
            _ => {
                // automatic, celect and hhrsh control include temperature charge cut logic
                let temp_charge_cut = self.temp_charge_cut_corr(simtime);

                if temp_charge_cut.is_some_and(|temp_charge_cut| {
                    temp_air.is_some_and(|temp_air| temp_air >= temp_charge_cut)
                }) {
                    // Control logic cut when temp_air is over temp_charge cut
                    target_charge_nominal = 0.;
                }

                match self.logic_type {
                    ControlLogicType::Automatic => {
                        // Automatic charge control can be achieved using internal thermostat(s) to
                        // control the extent of charging of the heaters. All or nothing approach

                        // Controls can also be supplemented by an external weather sensor,
                        // which tends to act as a limiting device to prevent the storage heaters from overcharging.
                        if self.external_sensor.is_some() && self.external_conditions.is_some() {
                            let limit = self.get_limit_factor(
                                self.external_conditions
                                    .as_ref()
                                    .unwrap()
                                    .air_temp(&simtime),
                            )?;
                            target_charge_nominal * limit
                        } else {
                            target_charge_nominal
                        }
                    }
                    ControlLogicType::Celect => {
                        // A CELECT-type controller has electronic sensors throughout the dwelling linked
                        // to a central control device. It monitors the individual room sensors and optimises
                        // the charging of all the storage heaters individually (and may select direct acting
                        // heaters in preference to storage heaters).
                        //
                        // Initial CELECT-type logic based on AUTOMATIC until additional literature for
                        // CELECT types is identified
                        //
                        // Controls can also be supplemented by an external weather sensor,
                        // which tends to act as a limiting device to prevent the storage heaters from overcharging.
                        if self.external_sensor.is_some() && self.external_conditions.is_some() {
                            let limit = self.get_limit_factor(
                                self.external_conditions
                                    .as_ref()
                                    .unwrap()
                                    .air_temp(&simtime),
                            )?;
                            target_charge_nominal * limit
                        } else {
                            target_charge_nominal
                        }
                    }
                    ControlLogicType::Hhrsh | ControlLogicType::HeatBattery => {
                        // A ‘high heat retention storage heater’ is one with heat retention not less
                        // than 45% measured according to BS EN 60531. It incorporates a timer, electronic
                        // room thermostat and fan to control the heat output. It is also able to estimate
                        // the next day’s heating demand based on external temperature, room temperature
                        // settings and heat demand periods.
                        if target_charge_nominal != 0. {
                            1.
                        } else {
                            0.
                        }
                    }
                    ControlLogicType::Manual => unreachable!(),
                }
            }
        };

        Ok(target_charge)
    }

    pub(crate) fn energy_to_store(
        &self,
        energy_demand: f64,
        base_temp: f64,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        // ugly, but this method cannot be called when control does not have HHRSH or Heat Battery logic type
        let heat_retention_data = if let Some(heat_retention_data) =
            self.heat_retention_data.as_ref()
        {
            heat_retention_data
        } else {
            unreachable!("energy_to_store() should not be called when control does not have HHRSH or HeatBattery logic type.");
        };
        let ChargeControlHeatRetentionFields {
            steps_day,
            demand,
            past_ext_temp,
            future_ext_temp,
            energy_to_store: energy_to_store_atomic,
        } = heat_retention_data;
        demand.write().push_front(Some(energy_demand));
        if self.external_conditions.is_some() {
            future_ext_temp.write().push_front(Some(
                self.external_conditions
                    .as_ref()
                    .unwrap()
                    .air_temp_with_offset(&simtime, *steps_day),
            ));
            past_ext_temp.write().push_front(Some(
                self.external_conditions
                    .as_ref()
                    .unwrap()
                    .air_temp(&simtime),
            ));
        }

        let future_hdh =
            self.calculate_heating_degree_hours(future_ext_temp.read().as_ref(), base_temp);
        let past_hdh =
            self.calculate_heating_degree_hours(past_ext_temp.read().as_ref(), base_temp);

        let energy_to_store = if !self.is_on(&simtime) {
            0.
        } else {
            match (future_hdh, past_hdh) {
                (None, _) | (_, None) => f64::NAN,
                // Can't calculate tomorrow's demand if no past_hdh, so assume zero to store
                (_, Some(0.)) => 0.,
                (Some(future_hdh), Some(past_hdh)) => {
                    future_hdh / past_hdh * demand.read().iter().flatten().sum::<f64>()
                }
            }
        };

        energy_to_store_atomic.store(energy_to_store, Ordering::SeqCst);

        energy_to_store
    }

    ///Correct nominal/json temp_charge_cut with monthly table
    /// Arguments
    /// returns -- temp_charge_cut (corrected)
    pub(crate) fn temp_charge_cut_corr(&self, simtime: SimulationTimeIteration) -> Option<f64> {
        if self.temp_charge_cut.is_none() {
            // Return None if temp_charge_cut is not set (e.g., for heat batteries)
            None
        } else {
            let temp_charge_cut_delta =
                if let Some(temp_charge_cut_delta) = self.temp_charge_cut_delta.as_ref() {
                    temp_charge_cut_delta
                        [simtime.time_series_idx(self.start_day, self.time_series_step)]
                } else {
                    0.0
                };

            Some(self.temp_charge_cut.unwrap() + temp_charge_cut_delta)
        }
    }

    fn calculate_heating_degree_hours(
        &self,
        temps: &VecDeque<Option<f64>>,
        base_temp: f64,
    ) -> Option<f64> {
        let mut total_hdh = 0.;
        for temp in temps {
            if let Some(temp) = temp {
                let hdh = (base_temp - *temp).max(0.);
                total_hdh += hdh;
            } else {
                return None;
            }
        }

        Some(total_hdh)
    }

    fn get_limit_factor(&self, external_temperature: f64) -> anyhow::Result<f64> {
        let correlation = &self.external_sensor.as_ref().expect("get_limit_factor should not be called on ChargeControl when there is no external sensor").correlation;

        // Edge cases: If temperature is below the first point or above the last point
        let first_correlation = correlation
            .first()
            .expect("External sensor correlation was not expected to be empty.");
        let last_correlation = correlation
            .last()
            .expect("External sensor correlation was not expected to be empty.");
        if external_temperature <= first_correlation.temperature {
            return Ok(first_correlation.max_charge);
        } else if external_temperature >= last_correlation.temperature {
            return Ok(last_correlation.max_charge);
        }

        // Linear interpolation
        for i in 1..correlation.len() {
            let ExternalSensorCorrelation {
                temperature: temp_1,
                max_charge: max_charge_1,
            } = correlation[i - 1];
            let ExternalSensorCorrelation {
                temperature: temp_2,
                max_charge: max_charge_2,
            } = correlation[i];

            if temp_1 <= external_temperature && external_temperature <= temp_2 && temp_1 != temp_2
            {
                // perform linear interpolation
                let slope = (max_charge_2 - max_charge_1) / (temp_2 - temp_1);
                let limit = max_charge_1 + slope * (external_temperature - temp_1);
                return Ok(limit);
            }
        }

        bail!("Calculation of limiting factor linked to external sensor for automatic control failed.")
    }
}

impl ControlBehaviour for ChargeControl {
    // In Python this is inherited from the BoolTimeControl class
    fn is_on(&self, iteration: &SimulationTimeIteration) -> bool {
        self.schedule[iteration.time_series_idx(self.start_day, self.time_series_step)]
    }
}

#[derive(Clone, Debug)]
pub(crate) struct OnOffCostMinimisingTimeControl {
    schedule: Vec<f64>,
    start_day: u32,
    time_series_step: f64,
    time_on_daily: f64,
    on_off_schedule: Vec<bool>,
}

impl OnOffCostMinimisingTimeControl {
    /// Construct an OnOffCostMinimisingControl object
    ///
    /// Arguments:
    /// * `schedule` - list of cost values (one entry per time_series_step)
    /// * `simulation_time` - reference to SimulationTime object
    /// * `start_day` - first day of the time series, day of the year, 0 to 365 (single value)
    /// * `time_series_step` - timestep of the time series data, in hours
    /// * `time_on_daily` - number of "on" hours to be set per day
    pub(crate) fn new(
        schedule: Vec<f64>,
        start_day: u32,
        time_series_step: f64,
        time_on_daily: f64,
    ) -> anyhow::Result<Self> {
        let timesteps_per_day = (HOURS_IN_DAY as f64 / time_series_step) as usize;
        let timesteps_on_daily = (time_on_daily / time_series_step) as usize;
        let time_series_len_days =
            ((schedule.len() as f64 * time_series_step / HOURS_IN_DAY as f64).ceil()) as usize;

        // For each day of schedule, find the specified number of hours with the lowest cost
        let mut on_off_schedule: Vec<bool> = vec![];
        for day in 0..time_series_len_days {
            // Get part of the schedule for current day
            let schedule_day_start = day * timesteps_per_day;
            let schedule_day_end = schedule_day_start + timesteps_per_day;

            // Added below check in the Rust before we try to access `schedule` elements by range
            // just below. This ensures that we handle the case when the end of the range is greater
            // than the length of the schedule (otherwise we'd get a panic in Rust). Python is more
            // lenient and will assume access elements up to the last one and will not error.
            if schedule.len() < schedule_day_end {
                bail!("There is a mismatch between the schedule length and the timesteps per day (hours_per_day / time_series_step)")
            }
            let schedule_day = schedule[schedule_day_start..schedule_day_end].to_vec();

            // Find required number of timesteps with lowest costs
            let mut schedule_day_cost_lowest = schedule_day.clone();
            schedule_day_cost_lowest.sort_by(f64::total_cmp);
            let schedule_day_cost_lowest = schedule_day_cost_lowest[0..timesteps_on_daily]
                .iter()
                .dedup()
                .collect_vec();

            // Initialise boolean schedule for day
            let mut schedule_onoff_day = vec![false; timesteps_per_day];

            // Set lowest cost times to True, then next lowest etc. until required
            // number of timesteps have been set to True
            if schedule_onoff_day.len() != schedule_day.len() {
                bail!("Different lengths for schedule_onoff_day and schedule_day")
            }

            let mut timesteps_to_be_allocated = timesteps_on_daily;
            for cost in schedule_day_cost_lowest.iter() {
                for (idx, entry) in schedule_day.iter().enumerate() {
                    if timesteps_to_be_allocated < 1 {
                        break;
                    }
                    if entry == *cost {
                        schedule_onoff_day[idx] = true;
                        timesteps_to_be_allocated -= 1;
                    }
                }
            }
            // Add day of schedule to overall
            on_off_schedule.extend(schedule_onoff_day);
        }

        Ok(OnOffCostMinimisingTimeControl {
            schedule,
            start_day,
            time_series_step,
            time_on_daily,
            on_off_schedule,
        })
    }
}

impl ControlBehaviour for OnOffCostMinimisingTimeControl {
    /// Return true if control will allow system to run
    fn is_on(&self, timestep: &SimulationTimeIteration) -> bool {
        self.on_off_schedule[timestep.time_series_idx(self.start_day, self.time_series_step)]
    }
}

#[derive(Clone, Debug)]
/// An object to model a control with a setpoint which varies per timestep
pub(crate) struct SetpointTimeControl {
    /// list of float values (one entry per hour)
    schedule: Vec<Option<f64>>,
    /// first day of the time series, day of the year, 0 to 365 (single value)
    start_day: u32,
    /// timestep of the time series data, in hours
    time_series_step: f64,
    setpoint_bounds: SetpointBounds,
    /// how long before heating period the system
    /// should switch on, in hours
    timesteps_advstart: u32,
}

impl SetpointTimeControl {
    /// Construct a SetpointTimeControl object
    ///
    /// Arguments:
    /// * `schedule` - list of float values (one entry per hour)
    /// * `simulation_time` - reference to SimulationTime object
    /// * `start_day` - first day of the time series, day of the year, 0 to 365 (single value)
    /// * `time_series_step` - timestep of the time series data, in hours
    /// * `setpoint_min` - min setpoint allowed
    /// * `setpoint_max` - max setpoint allowed
    /// * `default_to_max` - if both min and max limits are set but setpoint isn't,
    ///                      whether to default to min (False) or max (True)
    /// * `duration_advanced_start` - how long before heating period the system
    ///                               should switch on, in hours
    pub(crate) fn new(
        schedule: Vec<Option<f64>>,
        start_day: u32,
        time_series_step: f64,
        setpoint_bounds: Option<SetpointBoundsInput>,
        duration_advanced_start: Option<f64>,
        timestep: f64,
    ) -> Self {
        let duration_advanced_start = duration_advanced_start.unwrap_or(0.0);
        Self {
            schedule, // in Python now part of initialising base class (FloatOrNoneTimeControl > BaseTimeControl)
            start_day, // in Python now part of initialising base class (FloatOrNoneTimeControl > BaseTimeControl)
            time_series_step, // in Python now part of initialising base class (FloatOrNoneTimeControl > BaseTimeControl)
            setpoint_bounds: setpoint_bounds.into(),
            timesteps_advstart: (duration_advanced_start / timestep).round() as u32,
        }
    }

    // in_required_period method can be found here in Python, in Rust it's part of the
    // implementation block of ControlBehaviour further down

    fn is_on_for_timestep_idx(&self, schedule_idx: usize) -> bool {
        let setpnt = self.schedule[schedule_idx];

        if setpnt.is_none() {
            // Look ahead for duration of warmup period: system is on if setpoint
            // is not None heating period if found
            for timesteps_ahead in 1..(self.timesteps_advstart + 1) {
                let timesteps_ahead = timesteps_ahead as usize;
                if self.schedule.len() <= (schedule_idx + timesteps_ahead) {
                    // Stop looking ahead if we have reached the end of the schedule
                    break;
                }
                if self.schedule[schedule_idx + timesteps_ahead].is_some() {
                    // If heating period starts within duration of warmup period
                    // from now, system is on
                    return true;
                }
            }
        }

        // For this type of control, system is always on if min or max are set
        !(setpnt.is_none() && matches!(self.setpoint_bounds, SetpointBounds::NoSetpoints))
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum SetpointBounds {
    MinAndMax {
        /// Minimum setpoint allowed
        setpoint_min: f64,

        /// Maximum setpoint allowed
        setpoint_max: f64,

        /// If both min and max limits are set but setpoint is not, whether to default to min (false) or max (true)
        default_to_max: bool,
    },
    MinOnly {
        /// Minimum setpoint allowed
        setpoint_min: f64,
    },
    MaxOnly {
        /// Maximum setpoint allowed
        setpoint_max: f64,
    },
    NoSetpoints,
}

impl SetpointBounds {
    fn setpoint_max(&self) -> Option<f64> {
        match self {
            Self::MinAndMax { setpoint_max, .. } => Some(*setpoint_max),
            Self::MaxOnly { setpoint_max } => Some(*setpoint_max),
            _ => None,
        }
    }

    fn setpoint_min(&self) -> Option<f64> {
        match self {
            Self::MinAndMax { setpoint_min, .. } => Some(*setpoint_min),
            Self::MinOnly { setpoint_min } => Some(*setpoint_min),
            _ => None,
        }
    }
}

impl From<Option<SetpointBoundsInput>> for SetpointBounds {
    fn from(value: Option<SetpointBoundsInput>) -> Self {
        match value {
            Some(input) => match input {
                SetpointBoundsInput::MinAndMax {
                    setpoint_min,
                    setpoint_max,
                    default_to_max,
                } => Self::MinAndMax {
                    setpoint_min,
                    setpoint_max,
                    default_to_max,
                },
                SetpointBoundsInput::MinOnly { setpoint_min } => Self::MinOnly { setpoint_min },
                SetpointBoundsInput::MaxOnly { setpoint_max } => Self::MaxOnly { setpoint_max },
            },
            None => Self::NoSetpoints,
        }
    }
}

impl ControlBehaviour for SetpointTimeControl {
    /// Return true if current time is inside specified time for heating/cooling
    ///
    /// (not including timesteps where system is only on due to min or max
    /// setpoint or advanced start)
    fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool> {
        let schedule_idx =
            simulation_time_iteration.time_series_idx(self.start_day, self.time_series_step);
        let setpnt = self.schedule[schedule_idx];
        Some(setpnt.is_some())
    }

    /// Return setpoint for the current timestep
    fn setpnt(&self, simulation_time_iteration: &SimulationTimeIteration) -> Option<f64> {
        let schedule_idx =
            simulation_time_iteration.time_series_idx(self.start_day, self.time_series_step);

        let mut setpnt = self.schedule[schedule_idx];

        if setpnt.is_none() {
            // Look ahead for duration of warmup period and use setpoint from
            // start of heating period if found
            for timesteps_ahead in 1..(self.timesteps_advstart + 1) {
                let timesteps_ahead = timesteps_ahead as usize;
                if self.schedule.len() <= (schedule_idx + timesteps_ahead) {
                    // Stop looking ahead if we have reached the end of the schedule
                    break;
                }
                if let Some(s) = self.schedule[schedule_idx + timesteps_ahead] {
                    // If heating period starts within duration of warmup period
                    // from now, use setpoint from start of heating period
                    setpnt = Some(s);
                }
            }
        }

        match (setpnt, self.setpoint_bounds) {
            // If no setpoint value is in the schedule, use the min/max if set
            (None, SetpointBounds::NoSetpoints) => None, // Use setpnt None
            (None, SetpointBounds::MaxOnly { setpoint_max }) => Some(setpoint_max),
            (None, SetpointBounds::MinOnly { setpoint_min }) => Some(setpoint_min),
            (
                None,
                SetpointBounds::MinAndMax {
                    setpoint_min,
                    setpoint_max,
                    default_to_max,
                },
            ) => match default_to_max {
                true => Some(setpoint_max),
                false => Some(setpoint_min),
            },
            (Some(s), _) => {
                let mut setpnt = s;

                if let Some(setpoint_max) = self.setpoint_bounds.setpoint_max() {
                    // If there is a maximum limit, take the lower of this and the schedule value
                    setpnt = min_of_2(setpoint_max, setpnt);
                }
                if let Some(setpoint_min) = self.setpoint_bounds.setpoint_min() {
                    // If there is a minimum limit, take the higher of this and the schedule value
                    setpnt = max_of_2(setpoint_min, setpnt);
                }
                Some(setpnt)
            }
        }
    }

    /// Return true if control will allow system to run
    fn is_on(&self, timestep: &SimulationTimeIteration) -> bool {
        let schedule_idx = timestep.time_series_idx(self.start_day, self.time_series_step);

        self.is_on_for_timestep_idx(schedule_idx)
    }
}

/// An object for managing loadshifting appliances
#[derive(Debug)]
pub(crate) struct SmartApplianceControl {
    appliance_names: Vec<String>,
    energy_supplies: IndexMap<String, Arc<RwLock<EnergySupply>>>,
    battery_states_of_charge: IndexMap<String, Vec<AtomicF64>>,
    ts_power: IndexMap<String, Vec<AtomicF64>>,
    ts_step: f64,
    simulation_timestep: f64,
    ts_step_ratio: f64,
    non_appliance_demand_24hr: IndexMap<String, Vec<AtomicF64>>,
    buffer_length: usize,
}

impl SmartApplianceControl {
    /// Construct a SmartApplianceControl object
    ///
    /// Arguments:
    /// * `power_timeseries` - dictionary of lists containing expected power for appliances
    ///                        for each energy supply, for the entire length of the simulation
    /// * `timeseries_step` - timestep of the power timeseries
    ///                       (not necessarily equal to simulation_time.timestep())
    /// * `simulation_time` - reference to a SimulationTime object
    /// * `non_appliance_demand_24hr` - dictionary of lists containing 24 hour buffers of
    ///                                 demand per end user for each energy supply
    /// * `battery_24hr` - dictionary of lists containing 24 hour buffers of
    ///                    battery state of charge for each energy supply
    /// * `energysupplies` - dictionary of energysupply objects in the simulation
    /// * `appliances` - list of names of all appliance objects in the simulation
    pub(crate) fn new(
        power_timeseries: &IndexMap<String, Vec<f64>>,
        timeseries_step: f64,
        simulation_time_iterator: &SimulationTimeIterator,
        non_appliance_demand_24hr: IndexMap<String, Vec<f64>>,
        battery_24hr: SmartApplianceBattery,
        energy_supplies: &IndexMap<String, Arc<RwLock<EnergySupply>>>,
        appliance_names: Vec<String>,
    ) -> anyhow::Result<Self> {
        let energy_supplies: IndexMap<String, Arc<RwLock<EnergySupply>>> = energy_supplies
            .iter()
            .filter(|(key, _)| power_timeseries.contains_key(*key))
            .map(|(k, v)| (k.to_owned(), v.clone()))
            .collect();
        let battery_states_of_charge = energy_supplies
            .iter()
            .filter(|&(_name, supply)| supply.read().has_battery())
            .map(|(name, _supply)| {
                (
                    name.to_owned(),
                    battery_24hr.battery_state_of_charge[name]
                        .iter()
                        .map(|x| AtomicF64::new(*x))
                        .collect_vec(),
                )
            })
            .collect();
        for energy_supply in energy_supplies.keys() {
            if power_timeseries[energy_supply].len() as f64 * timeseries_step
                < simulation_time_iterator.total_steps() as f64
                    * simulation_time_iterator.step_in_hours()
            {
                bail!("ERROR: loadshifting power timeseries shorter than simulation length")
            }
        }
        Ok(Self {
            appliance_names,
            energy_supplies,
            battery_states_of_charge,
            ts_power: power_timeseries.iter().map(|(name, series)| (name.to_owned(), series.iter().map(|x| AtomicF64::new(*x)).collect_vec())).collect(),
            ts_step: timeseries_step,
            simulation_timestep: simulation_time_iterator.step_in_hours(),
            ts_step_ratio: simulation_time_iterator.step_in_hours() / timeseries_step,
            buffer_length: non_appliance_demand_24hr.clone().first().as_ref().ok_or_else(|| anyhow!("non_appliance_demand_24hr parameter for SmartApplianceControl cannot be empty."))?.1.len(),
            non_appliance_demand_24hr: non_appliance_demand_24hr.clone()
                .into_iter()
                .map(|(key, value)| (key, value.iter().map(|x| AtomicF64::new(*x)).collect_vec()))
                .collect(),
        })
    }

    fn ts_step(&self, t_idx: usize) -> usize {
        // converts index of simulation time to index of demand or weight timeseries
        (self.ts_step_ratio * t_idx as f64).floor() as usize
    }

    /// returns average energy demand from power timeseries over the current simulation timestep
    fn get_ts_demand(&self, energy_supply: &str, t_idx: usize) -> f64 {
        self.ts_power[energy_supply][self.ts_step(t_idx)].load(Ordering::SeqCst)
            / WATTS_PER_KILOWATT as f64
            * self.simulation_timestep
    }

    /// returns the sum of the anticipated appliance demand,
    /// the demand buffer, and the (negative) battery charge
    pub(crate) fn get_demand(&self, t_idx: usize, energy_supply: &str) -> f64 {
        let idx_24hr = t_idx % self.buffer_length;
        let demand = self.get_ts_demand(energy_supply, t_idx)
            + self.non_appliance_demand_24hr[energy_supply][idx_24hr].load(Ordering::SeqCst);

        if self.battery_states_of_charge.contains_key(energy_supply) {
            demand - self.battery_states_of_charge[energy_supply][idx_24hr].load(Ordering::SeqCst)
        } else {
            demand
        }
    }

    pub(crate) fn add_appliance_demand(
        &self,
        simtime: SimulationTimeIteration,
        demand: f64,
        energy_supply: &str,
    ) {
        // convert demand from appliance usage event to average power over the demand series timestep
        // and add it to the series
        let t_idx = simtime.index;
        self.ts_power[energy_supply][self.ts_step(t_idx)].fetch_add(
            demand * WATTS_PER_KILOWATT as f64 / self.ts_step,
            Ordering::SeqCst,
        );

        // update our prediction of battery charge over the next 24 hours
        if self.battery_states_of_charge.contains_key(energy_supply) {
            // if we expect there will be charge in the battery when this demand occurs, assume
            // the battery supplies as much of it as possible
            let idx_24hr = t_idx % self.buffer_length;
            let max_capacity = self.energy_supplies[energy_supply]
                .read()
                .get_battery_max_capacity()
                .expect("Battery expected to be present and reporting max capacity");
            // max_discharge is a linear function however states it requires input as a 0-1 proportion of total,
            // so divide and then multiply by max capacity in case of future changes
            let max_discharge = -self.energy_supplies[energy_supply]
                .read()
                .get_battery_max_discharge(
                    self.battery_states_of_charge[energy_supply][idx_24hr].load(Ordering::SeqCst)
                        / max_capacity,
                )
                .expect("Battery expected to be present and reporting max capacity")
                * max_capacity;
            // the maths here follows charge_discharge_battery() in ElectricBattery
            let discharge_efficiency = self.energy_supplies[energy_supply]
                .read()
                .get_battery_discharge_efficiency(simtime)
                .expect("Battery expected to be present and reporting max capacity");
            let charge_utilised = min_of_2(
                self.battery_states_of_charge[energy_supply][idx_24hr]
                    .load(Ordering::SeqCst)
                    .max(0.),
                max_discharge.min(demand) * discharge_efficiency,
            );
            // now subtract charge_utilised from the charge stored at every step in the buffer of battery charge.
            // if the battery is already expected to empty at a later time, this will result in the buffer
            // reporting negative charge stored in the battery during the times it is expected to be empty
            // and appliance preferentially not being used at those times
            self.battery_states_of_charge[energy_supply]
                .iter()
                .for_each(|x| {
                    x.fetch_sub(charge_utilised, Ordering::SeqCst);
                });
        }
    }

    pub(crate) fn update_demand_buffer(&self, simtime: SimulationTimeIteration) {
        let t_idx = simtime.index;
        let idx_24hr = t_idx % self.buffer_length;
        for (name, supply) in self.energy_supplies.iter() {
            // total up results for this energy supply but exclude demand from appliances
            // (the demand for which we already know accurately in advance
            // energy generated is negative and if the generation exceeds demand for a given timestep, the total
            // will be negative, and appliance usage events will preferentially be scheduled at that time
            // TODO (from Python) - it is possible to apply a weighting factor to energy generated in the dwelling here
            // (users for whom demand is negative)
            // to make it more or less preferable to use it immediately or export/charge battery
            self.non_appliance_demand_24hr[name][idx_24hr].store(
                supply
                    .read()
                    .results_by_end_user_single_step(t_idx)
                    .iter()
                    .filter_map(|(name, user)| {
                        (!self.appliance_names.contains(name)).then_some(user)
                    })
                    .sum::<f64>(),
                Ordering::SeqCst,
            );

            let supply = supply.read();
            if supply.has_battery() {
                // TODO (from Python) communicate with charge control
                let charge = supply
                    .get_battery_available_charge()
                    .expect("Battery expected to be present and reporting available charge");
                let charge_efficiency = supply
                    .get_battery_charge_efficiency(simtime)
                    .expect("Battery expected to be present and reporting charge efficiency");
                self.battery_states_of_charge[name][idx_24hr]
                    .store(charge * charge_efficiency, Ordering::SeqCst);
            }
        }
    }
}

impl ControlBehaviour for SmartApplianceControl {}

#[derive(Debug)]
/// An object to model a control with nested combinations of other control types
pub(crate) struct CombinationTimeControl {
    combinations: ControlCombinations,
    controls: IndexMap<String, Arc<Control>>,
}

impl CombinationTimeControl {
    /// Construct a CombinationTimeControl object
    ///
    /// Arguments:
    /// * `combination` - mapping of combination names to combination configurations (read-only)
    /// * `controls` - mapping of control names to control instances (read-only)
    /// * `simulation_time` - reference to SimulationTime object
    pub(crate) fn new(
        combinations: ControlCombinations,
        controls: IndexMap<String, Arc<Control>>,
    ) -> anyhow::Result<Self> {
        Self::validate_combinations(&combinations, &controls)?;

        Ok(Self {
            combinations,
            controls,
        })
    }

    // Unlike the upstream Python, we want to validate combinations on the way in so they can't fail
    // during a simulation
    // (Add more conditions if possible)
    fn validate_combinations(
        combinations: &ControlCombinations,
        _controls: &IndexMap<String, Arc<Control>>,
    ) -> anyhow::Result<()> {
        for (name, combination) in [("main", &combinations.main)].into_iter().chain(
            combinations
                .references
                .iter()
                .map(|(name, control)| (name.as_str(), control)),
        ) {
            if combination.controls.len() < 2 {
                bail!("Control combination {name} references fewer than two controls");
            }
        }
        Ok(())
    }

    fn evaluate_boolean_operation_is_on(
        &self,
        operation: ControlCombinationOperation,
        control_results: &mut impl Iterator<Item = bool>,
    ) -> bool {
        match operation {
            ControlCombinationOperation::And => control_results.all(|x| x),
            ControlCombinationOperation::Or => control_results.any(|x| x),
            ControlCombinationOperation::Xor => control_results.filter(|x| *x).count() % 2 == 1,
            ControlCombinationOperation::Not => {
                if control_results.count() != 1 {
                    unreachable!()
                }
                !control_results
                    .next()
                    .expect("Expected control results to have count of 1")
            }
            _ => unreachable!(),
        }
    }

    /// Evaluate a single control
    fn evaluate_control_is_on(&self, control_name: &str, simtime: SimulationTimeIteration) -> bool {
        let control = self.controls[control_name].as_ref();
        control.is_on(&simtime)
    }

    /// Evaluate a combination of controls
    fn evaluate_combination_is_on(
        &self,
        combination_name: &str,
        simtime: SimulationTimeIteration,
    ) -> bool {
        let combination = &self.combinations[combination_name];
        let ControlCombination {
            operation,
            controls,
        } = combination;

        let mut results = controls.iter().map(|control| {
            if self.combinations.contains_key(control) {
                // If the control is a combination, recursively evaluate it
                // Infinite recursion has been avoided by adding checks during control object creation
                self.evaluate_combination_is_on(control, simtime)
            } else {
                // Otherwise, evaluate a single control
                self.evaluate_control_is_on(control, simtime)
            }
        });
        match operation {
            ControlCombinationOperation::And
            | ControlCombinationOperation::Or
            | ControlCombinationOperation::Xor
            | ControlCombinationOperation::Not => {
                self.evaluate_boolean_operation_is_on(*operation, &mut results)
            }
            ControlCombinationOperation::Max
            | ControlCombinationOperation::Min
            | ControlCombinationOperation::Mean => results.any(|x| x),
        }
    }

    /// Evaluate a single control
    fn evaluate_control_in_req_period(
        &self,
        control_name: &str,
        simtime: SimulationTimeIteration,
    ) -> bool {
        let control = self.controls[control_name].as_ref();
        match control {
            c @ Control::OnOffTime(_)
            | c @ Control::Charge(_)
            | c @ Control::OnOffMinimisingTime(_) => c.is_on(&simtime),
            Control::SetpointTime(c) => c
                .in_required_period(&simtime)
                .expect("SetpointTimeControl in_required_period() method will always return Some"),
            _ => unreachable!("CombinationTimeControl only combined OnOffTime, Charge, OnOffMinimisingTime or SetpointTime controls"),
        }
    }

    /// This function processes a combination of control elements applying boolean logic (AND, OR, XOR, etc.) to their evaluation results.
    /// It checks the type of controls , validates allowed combinations and returns the evaluation result based on the specified operation.
    /// Unsupported combinations or operations raise an error.
    fn evaluate_combination_in_req_period(
        &self,
        combination_name: &str,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<bool> {
        let combination = &self.combinations[combination_name];
        let ControlCombination {
            operation,
            controls,
        } = combination;

        let has_on_off = AtomicBool::new(false);
        let has_setpoint = AtomicBool::new(false);

        let results = controls.iter().map(|control| {
            anyhow::Ok(if self.combinations.contains_key(control) {
                self.evaluate_combination_in_req_period(control, simtime)?
            } else {
                // Track the types of controls for logic enforcement
                match self.controls[control].as_ref() {
                    Control::OnOffTime(_)
                    | Control::Charge(_)
                    | Control::OnOffMinimisingTime(_) => {
                        has_on_off.store(true, Ordering::SeqCst);
                    }
                    Control::SetpointTime(_) => {
                        has_setpoint.store(true, Ordering::SeqCst);
                    }
                    _ => unreachable!("CombinationTimeControl only combined OnOffTime, Charge, OnOffMinimisingTime or SetpointTime controls"),
                }
                self.evaluate_control_in_req_period(control, simtime)
            })
        }).collect_vec(); // need to collect the iterator here so has_* variables are populated

        // Ensure valid combinations
        Ok(
            match (
                has_on_off.load(Ordering::SeqCst),
                has_setpoint.load(Ordering::SeqCst),
            ) {
                (true, true) => {
                    if !matches!(operation, ControlCombinationOperation::And) {
                        bail!("OnOff + Setpoint combination in_req_period() only supports the AND operation")
                    }
                    // Combine results using AND for OnOff + Setpoint combination
                    results
                        .into_iter()
                        .process_results(|mut iter| iter.all(|x| x))?
                }
                (true, false) => {
                    bail!(
                        "OnOff + OnOff combination is not applicable for in_req_period() operation"
                    )
                }
                (false, true) => {
                    // Apply operations for Setpoint + Setpoint combinations based on the operation
                    let results = results.into_iter();
                    match operation {
                        ControlCombinationOperation::And => {
                            results.process_results(|mut iter| iter.all(|x| x))?
                        }
                        ControlCombinationOperation::Or => {
                            results.process_results(|mut iter| iter.any(|x| x))?
                        }
                        ControlCombinationOperation::Xor => {
                            // XOR is true if exactly one result is true
                            results.process_results(|iter| iter.filter(|x| *x).count() == 1)?
                        }
                        ControlCombinationOperation::Max => results.process_results(|iter| {
                            iter.max().expect("At least one result was expected")
                        })?,
                        ControlCombinationOperation::Min => results.process_results(|iter| {
                            iter.min().expect("At least one result was expected")
                        })?,
                        ControlCombinationOperation::Mean => {
                            // Mean evaluates to True if average > 0.5
                            let results = results.collect::<anyhow::Result<Vec<_>>>()?;
                            results.iter().cloned().map(f64::from).sum::<f64>()
                                / results.len() as f64
                                > 0.5
                        }
                        _ => {
                            bail!("Unsupported combination operation encountered ('{operation:?}')")
                        }
                    }
                }
                _ => bail!("Invalid combination of controls encountered. No OnOff or Setpoint in in_req_period()"),
            },
        )
    }

    /// Evaluate a single control
    fn evaluate_control_setpnt(
        &self,
        control_name: &str,
        simtime: SimulationTimeIteration,
    ) -> SetpointOrBoolean {
        let control = self.controls[control_name].as_ref();
        match control {
            c @ Control::OnOffTime(_)
            | c @ Control::Charge(_)
            | c @ Control::OnOffMinimisingTime(_) => SetpointOrBoolean::Boolean(c.is_on(&simtime)),
            Control::SetpointTime(c) => SetpointOrBoolean::Setpoint(c.setpnt(&simtime)),
            _ => unreachable!("CombinationTimeControl only combined OnOffTime, Charge, OnOffMinimisingTime or SetpointTime controls"),
        }
    }

    /// Evaluate a combination of controls
    fn evaluate_combination_setpnt(
        &self,
        combination_name: &str,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<SetpointOrBoolean> {
        let combination = &self.combinations[combination_name];
        let ControlCombination {
            operation,
            controls,
        } = combination;

        let has_on_off = AtomicBool::new(false);
        let has_setpoint = AtomicBool::new(false);

        let results = controls
            .iter()
            .map(|control_name| {
                Ok(if self.combinations.contains_key(control_name) {
                    self.evaluate_combination_setpnt(control_name, simtime)?
                } else {
                    // Track the types of controls for logic enforcement
                    match self.controls[control_name].as_ref() {
                        Control::OnOffTime(_)
                        | Control::Charge(_)
                        | Control::OnOffMinimisingTime(_) => {
                            has_on_off.store(true, Ordering::SeqCst);
                        }
                        Control::SetpointTime(_) => {
                            has_setpoint.store(true, Ordering::SeqCst);
                        }
                        _ => unreachable!("CombinationTimeControl only combined OnOffTime, Charge, OnOffMinimisingTime or SetpointTime controls"),
                    }
                    self.evaluate_control_setpnt(control_name, simtime)
                })
            })
            .collect::<anyhow::Result<Vec<_>>>()?;
        // Check a setpnt result is available from previous combination
        if results.iter().any(|x| {
            matches!(
                x,
                SetpointOrBoolean::Setpoint(Some(_)) | SetpointOrBoolean::Boolean(_)
            )
        }) {
            has_setpoint.store(true, Ordering::SeqCst);
        }

        // Ensure valid combinations
        Ok(
            match (
                has_on_off.load(Ordering::SeqCst),
                has_setpoint.load(Ordering::SeqCst),
            ) {
                (true, true) => {
                    if matches!(operation, ControlCombinationOperation::And) {
                        let setpnt_value = results
                            .iter()
                            .filter_map(|x| {
                                if let SetpointOrBoolean::Setpoint(Some(t)) = x {
                                    Some(*t)
                                } else {
                                    None
                                }
                            })
                            .collect_vec();
                        let bool_value = results
                            .iter()
                            .filter_map(|x| {
                                if let SetpointOrBoolean::Boolean(b) = x {
                                    Some(b)
                                } else {
                                    None
                                }
                            })
                            .collect_vec();
                        if setpnt_value.len() > 1 {
                            bail!("Only one numerical value allowed in AND operation")
                        }
                        SetpointOrBoolean::Setpoint(if bool_value.into_iter().all(|x| *x) {
                            Some(*setpnt_value.first().ok_or_else(|| anyhow!("No setpoint values available when performing combination of controls"))?)
                        } else {
                            None
                        })
                    } else {
                        bail!(
                            "OnOff + Setpoint combination setpnt() only supports the AND operation"
                        )
                    }
                }
                (true, false) => {
                    bail!("OnOff + OnOff combination is not applicable for setpnt() operation")
                }
                (false, true) => {
                    // Apply operations for Setpoint + Setpoint combinations based on the operation
                    match operation {
                        ControlCombinationOperation::Max => {
                            *fallible_max_by(results.iter(), |a, b| {
                                Ok(if let (
                                    SetpointOrBoolean::Setpoint(Some(a)),
                                    SetpointOrBoolean::Setpoint(Some(b)),
                                ) = (a, b)
                                {
                                    a.total_cmp(b)
                                } else {
                                    bail!("Setpoint controls that are combined cannot have null values")
                                })
                            })?.expect("Results not expected to be empty")
                        }
                        ControlCombinationOperation::Min => {
                            *fallible_max_by(results.iter(), |a, b| {
                                Ok(if let (
                                    SetpointOrBoolean::Setpoint(Some(a)),
                                    SetpointOrBoolean::Setpoint(Some(b)),
                                ) = (a, b)
                                {
                                    a.total_cmp(b).reverse()
                                } else {
                                    bail!("Setpoint controls that are combined cannot have null values")
                                })
                            })?.expect("Results not expected to be empty")
                        }
                        ControlCombinationOperation::Mean => {
                            let results_sum = results
                                .iter()
                                .filter_map(|x| {
                                    if let SetpointOrBoolean::Setpoint(Some(t)) = x {
                                        Some(*t)
                                    } else {
                                        None
                                    }
                                })
                                .sum::<f64>();
                            SetpointOrBoolean::Boolean(results_sum / results.len() as f64 > 0.5)
                        }
                        _ => {
                            bail!("Unsupported combination operation encountered ('{operation:?}')")
                        }
                    }
                }
                _ => bail!("Invalid combination of controls encountered"),
            },
        )
    }

    #[cfg(test)]
    fn evaluate_control_target_charge(
        &self,
        control_name: &str,
        simtime: SimulationTimeIteration,
        temp_air: Option<f64>,
    ) -> anyhow::Result<Option<f64>> {
        let control = self.controls[control_name].as_ref();
        Ok(if let Control::Charge(c) = control {
            Some(c.target_charge(simtime, temp_air)?)
        } else {
            None
        })
    }

    /// Evaluate the combination for target charge
    #[cfg(test)]
    fn evaluate_combination_target_charge(
        &self,
        combination_name: &str,
        simtime: SimulationTimeIteration,
        temp_air: Option<f64>,
    ) -> anyhow::Result<f64> {
        let combination = &self.combinations[combination_name];
        let ControlCombination { controls, .. } = combination;

        let results = controls
            .iter()
            .map(|control_name| {
                Ok(if self.combinations.contains_key(control_name) {
                    // If the control is a combination, recursively evaluate it
                    // Infinite recursion has been avoided by adding checks during control object creation
                    Some(self.evaluate_combination_target_charge(
                        control_name,
                        simtime,
                        temp_air,
                    )?)
                } else {
                    self.evaluate_control_target_charge(control_name, simtime, temp_air)?
                })
            })
            .collect::<anyhow::Result<Vec<_>>>()?;

        Ok(if results.iter().all(Option::is_none) {
            bail!("Requires atleast one ChargeControl object in combination to determine target charge")
        } else if results.iter().filter_map(|&x| x).sum::<f64>() > 1. {
            bail!("CombinationControl cannot have more than one ChargeControl object to determine target charge")
        } else {
            results
                .iter()
                .filter_map(|&x| x)
                .next()
                .ok_or_else(|| anyhow!("Combination was expected not to be empty"))?
        })
    }

    #[cfg(test)]
    fn target_charge(
        &self,
        simtime: &SimulationTimeIteration,
        temp_air: Option<f64>,
    ) -> anyhow::Result<f64> {
        self.evaluate_combination_target_charge(MAIN_REFERENCE, *simtime, temp_air)
    }
}

impl ControlBehaviour for CombinationTimeControl {
    fn in_required_period(&self, simtime: &SimulationTimeIteration) -> Option<bool> {
        self.evaluate_combination_in_req_period(MAIN_REFERENCE, *simtime)
            .ok()
    }

    fn setpnt(&self, simtime: &SimulationTimeIteration) -> Option<f64> {
        self.evaluate_combination_setpnt("main", *simtime)
            .ok()
            .and_then(|x| {
                if let SetpointOrBoolean::Setpoint(Some(x)) = x {
                    Some(x)
                } else {
                    None
                }
            })
    }

    fn is_on(&self, simtime: &SimulationTimeIteration) -> bool {
        self.evaluate_combination_is_on(MAIN_REFERENCE, *simtime)
    }
}

#[cfg(test)]
#[derive(Copy, Clone, Debug, Default)]
/// A mock control implementation that allows setting canned responses for common control methods.
pub(crate) struct MockControl {
    canned_setpnt: Option<f64>,
    canned_is_on: Option<bool>,
    canned_in_req_period: Option<bool>,
}

#[cfg(test)]
impl MockControl {
    pub(crate) fn new(
        canned_setpnt: Option<f64>,
        canned_is_on: Option<bool>,
        canned_in_req_period: Option<bool>,
    ) -> Self {
        Self {
            canned_setpnt,
            canned_is_on,
            canned_in_req_period,
        }
    }

    pub(crate) fn with_setpnt(setpnt: Option<f64>) -> Self {
        Self {
            canned_setpnt: setpnt,
            ..Default::default()
        }
    }

    pub(crate) fn with_is_on(is_on: bool) -> Self {
        Self {
            canned_is_on: Some(is_on),
            ..Default::default()
        }
    }
}

#[cfg(test)]
impl ControlBehaviour for MockControl {
    fn in_required_period(&self, _simtime: &SimulationTimeIteration) -> Option<bool> {
        self.canned_in_req_period
    }

    fn setpnt(&self, _simtime: &SimulationTimeIteration) -> Option<f64> {
        self.canned_setpnt
    }

    fn is_on(&self, _simtime: &SimulationTimeIteration) -> bool {
        self.canned_is_on.unwrap_or(true)
    }
}

#[derive(Clone, Copy)]
enum SetpointOrBoolean {
    Setpoint(Option<f64>),
    Boolean(bool),
}

// utility functions to assist with comparing setpoint values where the comparisons are fallible
// (because NULLs are implicitly not permitted in combination controls, but may be encountered)

fn fallible_max_by<T, F, E>(iter: impl Iterator<Item = T>, mut compare: F) -> Result<Option<T>, E>
where
    F: FnMut(&T, &T) -> Result<std::cmp::Ordering, E>,
{
    try_reduce(iter, |a, b| match compare(&a, &b)? {
        std::cmp::Ordering::Greater | std::cmp::Ordering::Equal => Ok(a),
        std::cmp::Ordering::Less => Ok(b),
    })
}

fn try_reduce<T, F, E>(mut iter: impl Iterator<Item = T>, mut f: F) -> Result<Option<T>, E>
where
    F: FnMut(T, T) -> Result<T, E>,
{
    let mut acc = match iter.next() {
        Some(value) => value,
        None => return Ok(None),
    };

    for item in iter {
        acc = f(acc, item)?;
    }

    Ok(Some(acc))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external_conditions::DaylightSavingsConfig;
    use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
    use pretty_assertions::assert_eq;
    use rstest::*;
    use serde_json::json;

    #[fixture]
    fn simulation_time() -> SimulationTimeIterator {
        SimulationTime::new(0.0, 8.0, 1.0).iter()
    }

    mod test_on_off_time_control {
        use super::*;
        use pretty_assertions::assert_eq;

        #[test]
        fn test_is_on() {
            let simulation_time_iterator = SimulationTime::new(0.0, 8.0, 1.0).iter();
            let schedule = [true, false, true, true, false, true, false, false];
            let time_control =
                OnOffTimeControl::new(schedule.iter().map(|&v| Some(v)).collect_vec(), 0, 1.0);

            for iteration in simulation_time_iterator {
                assert_eq!(time_control.is_on(&iteration), schedule[iteration.index]);
            }
        }
    }

    mod test_on_off_cost_minimising_time_control {
        use super::*;
        use pretty_assertions::assert_eq;

        #[test]
        fn test_init_invalid_schedule_length() {
            let cost_schedule = [vec![5.0; 7], vec![10.0; 2], vec![7.5; 8], vec![15.0; 6]]
                .to_vec()
                .concat();
            let cost_schedule = [&cost_schedule[..], &cost_schedule[..]].concat();
            let control = OnOffCostMinimisingTimeControl::new(cost_schedule, 0, 1.0, 12.0);
            assert!(control.is_err());
            let error = control.unwrap_err().to_string();
            assert_eq!(
                error,
                "There is a mismatch between the schedule length and the timesteps per day (hours_per_day / time_series_step)"
            );
        }

        #[rstest]
        fn test_is_on() {
            let schedule = [
                vec![5.0; 7],
                vec![10.0; 2],
                vec![7.5; 8],
                vec![15.0; 6],
                vec![5.0],
            ]
            .to_vec()
            .concat();
            let schedule = [&schedule[..], &schedule[..]].concat();
            let cost_minimising_ctrl =
                OnOffCostMinimisingTimeControl::new(schedule, 0, 1.0, 12.0).unwrap();

            let resulting_schedule = [
                vec![true; 7],
                vec![false; 2],
                vec![true; 4],
                vec![false; 4],
                vec![false; 6],
                vec![true],
            ]
            .to_vec()
            .concat();
            let resulting_schedule = [&resulting_schedule[..], &resulting_schedule[..]].concat();
            let simulation_time_iterator = SimulationTime::new(0.0, 48.0, 1.0).iter();
            for iteration in simulation_time_iterator {
                pretty_assertions::assert_eq!(
                    cost_minimising_ctrl.is_on(&iteration),
                    resulting_schedule[iteration.index]
                );
            }
        }
    }

    mod test_setpoint_time_control {
        use super::*;
        use pretty_assertions::assert_eq;

        #[fixture]
        fn simulation_time_iterator() -> SimulationTimeIterator {
            SimulationTime::new(0.0, 8.0, 1.0).iter()
        }
        fn default_schedule() -> Vec<Option<f64>> {
            vec![
                Some(21.0),
                None,
                None,
                Some(21.0),
                None,
                Some(21.0),
                Some(25.0),
                Some(15.0),
            ]
        }

        fn create_time_control(
            schedule: Option<Vec<Option<f64>>>,
            setpoint_bounds: Option<SetpointBoundsInput>,
            duration_advanced_start: Option<f64>,
        ) -> SetpointTimeControl {
            let schedule = schedule.unwrap_or(default_schedule());
            SetpointTimeControl::new(
                schedule,
                0,
                1.,
                setpoint_bounds,
                duration_advanced_start,
                1., // simulation_time.step_in_hours()
            )
        }

        #[fixture]
        fn time_control() -> SetpointTimeControl {
            create_time_control(None, None, None)
        }

        #[fixture]
        fn time_control_min() -> SetpointTimeControl {
            let setpoint_bounds = SetpointBoundsInput::MinOnly { setpoint_min: 16. };
            create_time_control(None, Some(setpoint_bounds), None)
        }

        #[fixture]
        fn time_control_max() -> SetpointTimeControl {
            let setpoint_bounds = SetpointBoundsInput::MaxOnly { setpoint_max: 24. };
            create_time_control(None, Some(setpoint_bounds), None)
        }

        #[fixture]
        fn time_control_min_max() -> SetpointTimeControl {
            let setpoint_bounds = SetpointBoundsInput::MinAndMax {
                setpoint_min: 16.,
                setpoint_max: 24.,
                default_to_max: false,
            };
            create_time_control(None, Some(setpoint_bounds), None)
        }

        #[fixture]
        fn time_control_advstart() -> SetpointTimeControl {
            create_time_control(None, None, Some(1.))
        }

        #[fixture]
        fn time_control_advstart_min_max() -> SetpointTimeControl {
            let setpoint_bounds = SetpointBoundsInput::MinAndMax {
                setpoint_min: 16.,
                setpoint_max: 24.,
                default_to_max: false,
            };
            create_time_control(None, Some(setpoint_bounds), Some(1.))
        }

        #[rstest]
        fn test_in_required_period(
            simulation_time_iterator: SimulationTimeIterator,
            time_control: SetpointTimeControl,
            time_control_min: SetpointTimeControl,
            time_control_max: SetpointTimeControl,
            time_control_min_max: SetpointTimeControl,
            time_control_advstart: SetpointTimeControl,
            time_control_advstart_min_max: SetpointTimeControl,
        ) {
            let expected = [true, false, false, true, false, true, true, true];
            for t_it in simulation_time_iterator {
                assert_eq!(
                    time_control.in_required_period(&t_it).unwrap(),
                    expected[t_it.index],
                    "incorrect in_required_period value returned for control with no min or max set, iteration {}",
                    t_it.index + 1
                );
                assert_eq!(
                    time_control_min.in_required_period(&t_it).unwrap(),
                    expected[t_it.index],
                    "incorrect in_required_period value returned for control with min set, iteration {}",
                    t_it.index + 1
                );
                assert_eq!(
                    time_control_max.in_required_period(&t_it).unwrap(),
                    expected[t_it.index],
                    "incorrect in_required_period value returned for control with max set, iteration {}",
                    t_it.index + 1
                );
                assert_eq!(
                    time_control_min_max.in_required_period(&t_it).unwrap(),
                    expected[t_it.index],
                    "incorrect in_required_period value returned for control with min and max set, iteration {}",
                    t_it.index + 1
                );
                assert_eq!(
                    time_control_advstart.in_required_period(&t_it).unwrap(),
                    expected[t_it.index],
                    "incorrect in_required_period value returned for control with advanced start, iteration {}",
                    t_it.index + 1
                );
                assert_eq!(
                    time_control_advstart_min_max
                        .in_required_period(&t_it)
                        .unwrap(),
                    expected[t_it.index],
                    "incorrect in_required_period value returned for control with advanced start, iteration {}",
                    t_it.index + 1
                );
            }
        }

        #[rstest]
        fn test_is_on(
            simulation_time_iterator: SimulationTimeIterator,
            time_control: SetpointTimeControl,
            time_control_min: SetpointTimeControl,
            time_control_max: SetpointTimeControl,
            time_control_min_max: SetpointTimeControl,
            time_control_advstart: SetpointTimeControl,
            time_control_advstart_min_max: SetpointTimeControl,
        ) {
            for t_it in simulation_time_iterator {
                assert_eq!(time_control.is_on(&t_it),
                           [true, false, false, true, false, true, true, true][t_it.index],
                           "incorrect is_on value returned for control with no min or max set, iteration {}",
                           t_it.index + 1 );
                assert_eq!(
                    time_control_min.is_on(&t_it),
                    true, // Should always be true for this type of control
                    "incorrect is_on value returned for control with min set, iteration {}",
                    t_it.index + 1
                );
                assert_eq!(
                    time_control_max.is_on(&t_it),
                    true, // Should always be true for this type of control
                    "incorrect is_on value returned for control with max set, iteration {}",
                    t_it.index + 1
                );
                assert_eq!(
                    time_control_min_max.is_on(&t_it),
                    true, // Should always be true for this type of control
                    "incorrect is_on value returned for control with min and max set, iteration {}",
                    t_it.index + 1
                );
                assert_eq!(
                    time_control_advstart.is_on(&t_it),
                    [true, false, true, true, true, true, true, true][t_it.index],
                    "incorrect is_on value returned for control with advanced start, iteration {}",
                    t_it.index + 1
                );
                assert_eq!(
                    time_control_advstart_min_max.is_on(&t_it),
                    true,
                    "incorrect is_on value returned for control with advanced start and min/max, iteration {}",
                    t_it.index + 1
                );
            }
        }

        #[rstest]
        fn test_is_on_lookahead(simulation_time_iterator: SimulationTimeIterator) {
            let schedule = vec![None; 24];
            let control = create_time_control(Some(schedule), None, Some(30.));
            assert!(!control.is_on(&simulation_time_iterator.current_iteration()));

            let schedule = vec![Some(20.); 24];
            let control = create_time_control(Some(schedule), None, None);
            assert!(control.is_on(&simulation_time_iterator.current_iteration()));
        }

        #[rstest]
        fn test_setpnt(
            time_control: SetpointTimeControl,
            time_control_min: SetpointTimeControl,
            time_control_max: SetpointTimeControl,
            time_control_min_max: SetpointTimeControl,
            time_control_advstart: SetpointTimeControl,
            time_control_advstart_min_max: SetpointTimeControl,
            simulation_time_iterator: SimulationTimeIterator,
        ) {
            let results_min: [Option<f64>; 8] = [
                Some(21.0),
                Some(16.0),
                Some(16.0),
                Some(21.0),
                Some(16.0),
                Some(21.0),
                Some(25.0),
                Some(16.0),
            ];
            let results_max: [Option<f64>; 8] = [
                Some(21.0),
                Some(24.0),
                Some(24.0),
                Some(21.0),
                Some(24.0),
                Some(21.0),
                Some(24.0),
                Some(15.0),
            ];
            let results_minmax: [Option<f64>; 8] = [
                Some(21.0),
                Some(16.0),
                Some(16.0),
                Some(21.0),
                Some(16.0),
                Some(21.0),
                Some(24.0),
                Some(16.0),
            ];
            let results_advstart: [Option<f64>; 8] = [
                Some(21.0),
                None,
                Some(21.0),
                Some(21.0),
                Some(21.0),
                Some(21.0),
                Some(25.0),
                Some(15.0),
            ];
            let results_advstart_minmax: [Option<f64>; 8] = [
                Some(21.0),
                Some(16.0),
                Some(21.0),
                Some(21.0),
                Some(21.0),
                Some(21.0),
                Some(24.0),
                Some(16.0),
            ];
            for t_it in simulation_time_iterator {
                assert_eq!(
                    time_control.setpnt(&t_it),
                    default_schedule()[t_it.index],
                    "incorrect schedule returned for control with no min or max set, iteration {}",
                    t_it.index + 1
                );
                assert_eq!(
                    time_control_min.setpnt(&t_it),
                    results_min[t_it.index],
                    "incorrect schedule returned for control with min set, iteration {}",
                    t_it.index + 1
                );
                assert_eq!(
                    time_control_max.setpnt(&t_it),
                    results_max[t_it.index],
                    "incorrect schedule returned for control with max set, iteration {}",
                    t_it.index + 1
                );
                assert_eq!(
                    time_control_min_max.setpnt(&t_it),
                    results_minmax[t_it.index],
                    "incorrect schedule returned for control with min and max set, iteration {}",
                    t_it.index + 1
                );
                assert_eq!(
                    time_control_advstart.setpnt(&t_it),
                    results_advstart[t_it.index],
                    "incorrect schedule returned for control with advanced start, iteration {}",
                    t_it.index + 1
                );
                assert_eq!(
                    time_control_advstart_min_max.setpnt(&t_it),
                    results_advstart_minmax[t_it.index],
                    "incorrect schedule returned for control with advanced start and min and max set, iteration {}",
                    t_it.index + 1
                );
            }
        }

        #[rstest]
        fn test_setpnt_lookahead(simulation_time_iterator: SimulationTimeIterator) {
            let schedule = vec![None; 24];
            let control = create_time_control(Some(schedule), None, Some(30.));
            assert!(control
                .setpnt(&simulation_time_iterator.current_iteration())
                .is_none());

            let schedule = vec![Some(20.); 24];
            let control = create_time_control(Some(schedule), None, Some(30.));
            assert_eq!(
                control
                    .setpnt(&simulation_time_iterator.current_iteration())
                    .unwrap(),
                20.
            );
        }

        #[rstest]
        fn test_setpnt_minmax(simulation_time_iterator: SimulationTimeIterator) {
            let schedule = vec![None; 24];
            let control = create_time_control(
                Some(schedule.clone()),
                SetpointBoundsInput::MinOnly { setpoint_min: 10. }.into(),
                None,
            );
            assert_eq!(
                control
                    .setpnt(&simulation_time_iterator.current_iteration())
                    .unwrap(),
                10.
            );

            // skip assertion as not possible to replicate (Rust is more prescriptive here, so not a problem)

            let control = create_time_control(
                Some(schedule.clone()),
                SetpointBoundsInput::MinAndMax {
                    setpoint_min: 10.,
                    setpoint_max: 20.,
                    default_to_max: true,
                }
                .into(),
                Some(30.),
            );
            assert_eq!(
                control
                    .setpnt(&simulation_time_iterator.current_iteration())
                    .unwrap(),
                20.
            );

            let control = create_time_control(
                Some(schedule.clone()),
                SetpointBoundsInput::MinAndMax {
                    setpoint_min: 10.,
                    setpoint_max: 20.,
                    default_to_max: false,
                }
                .into(),
                Some(30.),
            );
            assert_eq!(
                control
                    .setpnt(&simulation_time_iterator.current_iteration())
                    .unwrap(),
                10.
            );
        }
    }

    mod test_smart_appliance_control {
        use super::*;
        use crate::core::energy_supply::elec_battery::ElectricBattery;
        use crate::core::energy_supply::energy_supply::EnergySupplyBuilder;
        use crate::input::{BatteryLocation, FuelType};
        use approx::assert_relative_eq;
        use pretty_assertions::assert_eq;

        #[fixture]
        fn simulation_time_iterator() -> SimulationTimeIterator {
            SimulationTime::new(0., 24., 1.).iter()
        }

        #[fixture]
        // create dummy external conditions
        fn external_conditions(
            simulation_time_iterator: SimulationTimeIterator,
        ) -> ExternalConditions {
            ExternalConditions::new(
                &simulation_time_iterator,
                vec![0.0; 24],
                vec![3.7; 24],
                vec![200.; 24],
                vec![333.; 24],
                vec![0.; 24],
                vec![0.2; 8760],
                51.42,
                -0.75,
                0,
                0,
                None,
                1.,
                None,
                None,
                false,
                false,
                None,
            )
        }

        #[fixture]
        fn energy_supply(
            simulation_time_iterator: SimulationTimeIterator,
            external_conditions: ExternalConditions,
        ) -> Arc<RwLock<EnergySupply>> {
            let electric_battery = ElectricBattery::new(
                100., // significant for test_add_appliance_demand, in Python Magic Mock is used instead
                1., // significant for test_add_appliance_demand, in Python Magic Mock is used instead
                0., // significant for test_add_appliance_demand, in Python Magic Mock is used instead
                0.001,
                1.5,
                -100., // significant for test_add_appliance_demand, in Python Magic Mock is used instead
                BatteryLocation::Inside, // significant for test_add_appliance_demand, in Python Magic Mock is used instead
                false,
                simulation_time_iterator.step_in_hours(),
                Arc::new(external_conditions),
            );

            Arc::new(RwLock::new(
                EnergySupply::new(
                    FuelType::Electricity,
                    simulation_time_iterator.total_steps(),
                    None,
                    Some(electric_battery),
                    None,
                    None,
                )
                .unwrap(),
            ))
        }

        #[fixture]
        fn smart_appliance_control(
            simulation_time_iterator: SimulationTimeIterator,
            energy_supply: Arc<RwLock<EnergySupply>>,
        ) -> SmartApplianceControl {
            let power_timeseries = &IndexMap::from([("mains elec".into(), vec![100.; 12])]);
            let non_appliance_demand_24hr =
                IndexMap::from([("mains elec".into(), vec![[0.1, 0.2]; 6].into_flattened())]);
            let battery_state_of_charge: IndexMap<String, Vec<f64>> =
                IndexMap::from([("mains elec".into(), vec![0.5; 12])]);
            let battery_24hr = SmartApplianceBattery {
                battery_state_of_charge,
                energy_into_battery_from_generation: IndexMap::new(),
                energy_into_battery_from_grid: IndexMap::new(),
                energy_out_of_battery: IndexMap::new(),
            };
            let energy_supplies = &IndexMap::from([("mains elec".into(), energy_supply)]);

            SmartApplianceControl::new(
                power_timeseries,
                2.,
                &simulation_time_iterator,
                non_appliance_demand_24hr,
                battery_24hr,
                energy_supplies,
                vec!["Clothes_drying".into()],
            )
            .unwrap()
        }

        #[rstest]
        fn test_init_invalid_length(
            simulation_time_iterator: SimulationTimeIterator,
            energy_supply: Arc<RwLock<EnergySupply>>,
        ) {
            let battery_state_of_charge: IndexMap<String, Vec<f64>> =
                IndexMap::from([("mains elec".into(), vec![0.; 12])]);
            let battery_24hr = SmartApplianceBattery {
                battery_state_of_charge,
                energy_into_battery_from_generation: IndexMap::new(),
                energy_into_battery_from_grid: IndexMap::new(),
                energy_out_of_battery: IndexMap::new(),
            };
            let smart_appliance_control = SmartApplianceControl::new(
                &IndexMap::from([("mains elec".into(), vec![100.; 11])]),
                2.,
                &simulation_time_iterator,
                IndexMap::from([("mains elec".into(), vec![[0.1, 0.2]; 12].into_flattened())]),
                battery_24hr,
                &IndexMap::from([("mains elec".into(), energy_supply)]),
                vec!["Clothes_drying".into()],
            );

            assert!(smart_appliance_control.is_err());
        }

        #[rstest]
        fn test_ts_step(smart_appliance_control: SmartApplianceControl) {
            assert_eq!(smart_appliance_control.ts_step(0), 0);
            assert_eq!(smart_appliance_control.ts_step(23), 11);
            assert_eq!(smart_appliance_control.ts_step(24), 12);
        }

        #[rstest]
        fn test_add_appliance_demand(
            smart_appliance_control: SmartApplianceControl,
            mut simulation_time_iterator: SimulationTimeIterator,
        ) {
            let iteration = simulation_time_iterator.nth(5).unwrap();
            smart_appliance_control.add_appliance_demand(iteration, 100., "mains elec");
            assert_eq!(smart_appliance_control.get_demand(5, "mains elec"), -9950.2);
        }

        #[rstest]
        fn test_update_demand_buffer(
            smart_appliance_control: SmartApplianceControl,
            mut simulation_time_iterator: SimulationTimeIterator,
        ) {
            for t_it in simulation_time_iterator {
                smart_appliance_control.update_demand_buffer(t_it);
                assert_eq!(
                    smart_appliance_control.get_demand(t_it.index, "mains elec"),
                    0.1
                );
            }
        }

        #[rstest]
        fn test_get_demand(
            smart_appliance_control: SmartApplianceControl,
            mut simulation_time_iterator: SimulationTimeIterator,
        ) {
            assert_relative_eq!(smart_appliance_control.get_demand(0, "mains elec"), -0.3);
            assert_relative_eq!(smart_appliance_control.get_demand(1, "mains elec"), -0.2);
        }

        #[rstest]
        fn test_get_demand_no_battery(mut simulation_time_iterator: SimulationTimeIterator) {
            let smart_appliance_control = SmartApplianceControl::new(
                &IndexMap::from([("mains elec".into(), vec![100.; 12])]),
                2.,
                &simulation_time_iterator,
                IndexMap::from([("mains elec".into(), vec![[0.1, 0.2]; 12].into_flattened())]),
                SmartApplianceBattery {
                    battery_state_of_charge: IndexMap::new(),
                    energy_into_battery_from_generation: IndexMap::new(),
                    energy_into_battery_from_grid: IndexMap::new(),
                    energy_out_of_battery: IndexMap::new(),
                },
                &IndexMap::from([(
                    "mains elec".into(),
                    Arc::new(RwLock::new(
                        EnergySupplyBuilder::new(
                            FuelType::Electricity,
                            simulation_time_iterator.total_steps(),
                        )
                        .build(),
                    )),
                )]),
                vec!["Clothes_drying".into()],
            )
            .unwrap();

            assert_relative_eq!(smart_appliance_control.get_demand(0, "mains elec"), 0.2);
            assert_relative_eq!(smart_appliance_control.get_demand(1, "mains elec"), 0.3);
        }
    }

    #[fixture]
    fn simulation_time_for_charge_control() -> SimulationTime {
        SimulationTime::new(0.0, 24.0, 1.0)
    }

    #[fixture]
    fn schedule_for_charge_control() -> Vec<bool> {
        vec![
            true, true, true, true, true, true, true, true, false, false, false, false, false,
            false, false, false, true, true, true, true, false, false, false, false,
        ]
    }

    mod test_charge_control {
        use super::*;
        use pretty_assertions::assert_eq;

        #[fixture]
        fn charge_control(
            simulation_time_for_charge_control: SimulationTime,
            schedule_for_charge_control: Vec<bool>,
        ) -> ChargeControl {
            let external_conditions = ExternalConditions::new(
                &simulation_time_for_charge_control.iter(),
                vec![
                    19.0, 0.0, 1.0, 2.0, 5.0, 7.0, 6.0, 12.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0,
                    19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0,
                ],
                vec![
                    3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1, 3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1,
                    3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1,
                ],
                vec![
                    300., 250., 220., 180., 150., 120., 100., 80., 60., 40., 20., 10., 50., 100.,
                    140., 190., 200., 320., 330., 340., 350., 355., 315., 5.,
                ],
                vec![
                    0., 0., 0., 0., 35., 73., 139., 244., 320., 361., 369., 348., 318., 249., 225.,
                    198., 121., 68., 19., 0., 0., 0., 0., 0.,
                ],
                vec![
                    0., 0., 0., 0., 0., 0., 7., 53., 63., 164., 339., 242., 315., 577., 385., 285.,
                    332., 126., 7., 0., 0., 0., 0., 0.,
                ],
                vec![
                    0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                    0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                ],
                51.383,
                -0.783,
                0,
                0,
                Some(0),
                1.,
                Some(1),
                Some(DaylightSavingsConfig::NotApplicable),
                false,
                false,
                // following starts/ends are corrected from Python tests which erroneously use previous
                // "start" field instead of "start360" (which has different origin for angle)
                serde_json::from_value(json!([
                    {"start360": 0, "end360": 45},
                    {"start360": 45, "end360": 90},
                    {"start360": 90, "end360": 135},
                    {"start360": 135, "end360": 180,
                        "shading": [
                            {"type": "obstacle", "height": 10.5, "distance": 12}
                        ]
                    },
                    {"start360": 180, "end360": 225},
                    {"start360": 225, "end360": 270},
                    {"start360": 270, "end360": 315},
                    {"start360": 315, "end360": 360}
                ]))
                .unwrap(),
            );
            let external_sensor: ExternalSensor = serde_json::from_value(json!({
                "correlation": [
                    {"temperature": 0.0, "max_charge": 1.0},
                    {"temperature": 10.0, "max_charge": 0.9},
                    {"temperature": 18.0, "max_charge": 0.0}
                ]
            }))
            .unwrap();

            ChargeControl::new(
                ControlLogicType::Automatic,
                schedule_for_charge_control,
                &simulation_time_for_charge_control
                    .iter()
                    .current_iteration(),
                0,
                1.,
                vec![Some(1.0), Some(0.8)],
                Some(15.5),
                None,
                Some(external_conditions.into()),
                Some(external_sensor),
                None,
            )
            .unwrap()
        }

        #[rstest]
        fn test_is_on_for_charge_control(
            charge_control: ChargeControl,
            simulation_time_for_charge_control: SimulationTime,
            schedule_for_charge_control: Vec<bool>,
        ) {
            for (t_idx, t_it) in simulation_time_for_charge_control.iter().enumerate() {
                assert_eq!(
                    charge_control.is_on(&t_it),
                    schedule_for_charge_control[t_idx],
                    "incorrect schedule returned"
                );
            }
        }

        #[rstest]
        fn test_target_charge(
            charge_control: ChargeControl,
            simulation_time_for_charge_control: SimulationTime,
        ) {
            let expected_target_charges = (
                vec![
                    0.0,
                    1.0,
                    0.99,
                    0.98,
                    0.95,
                    0.93,
                    0.9400000000000001,
                    0.675,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                vec![
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            );

            for (t_idx, t_it) in simulation_time_for_charge_control.iter().enumerate() {
                assert_eq!(
                    charge_control.target_charge(t_it, Some(12.5)).unwrap(),
                    expected_target_charges.0[t_idx],
                    "incorrect target charge returned"
                );
            }
            for (t_idx, t_it) in simulation_time_for_charge_control.iter().enumerate() {
                assert_eq!(
                    charge_control.target_charge(t_it, Some(19.5)).unwrap(),
                    expected_target_charges.1[t_idx],
                    "incorrect target charge returned"
                );
            }
        }

        // (from Python) Check correction of nominal/json temp_charge_cut with monthly table.
        //               This function will most likely be superseded when the Electric Storage methodology
        //               is upgraded to consider more realistic manufacturers' controls and corresponding
        //               unit_test will be deprecated.
        #[rstest]
        fn test_temp_charge_cut_corr(
            charge_control: ChargeControl,
            simulation_time_for_charge_control: SimulationTime,
        ) {
            assert_eq!(
                charge_control
                    .temp_charge_cut_corr(simulation_time_for_charge_control.iter().next().unwrap())
                    .unwrap(),
                15.5
            );
        }
    }

    #[fixture]
    fn charge_control_for_combination(
        simulation_time_for_charge_control: SimulationTime,
        schedule_for_charge_control: Vec<bool>,
    ) -> ChargeControl {
        // in the upstream Python tests, a simulation time in injected into the external conditions object
        // used in tests for the combination control that is different than the one iterated on in the test,
        // which means that it is never iterated and reported external conditions are always as per the first
        // timestep. therefore the following is a changed external conditions object that repeats the first value
        // for e.g. air temps, in order to replicate the unrealistic behaviour in the Python.
        let external_conditions = ExternalConditions::new(
            &simulation_time_for_charge_control.iter(),
            vec![
                19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0,
                19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0,
            ],
            vec![
                3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9,
                3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9, 3.9,
            ],
            vec![
                300., 300., 300., 300., 300., 300., 300., 300., 300., 300., 300., 300., 300., 300.,
                300., 300., 300., 300., 300., 300., 300., 300., 300., 300.,
            ],
            vec![
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                0., 0., 0.,
            ],
            vec![
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                0., 0., 0.,
            ],
            vec![
                0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
            ],
            51.383,
            -0.783,
            0,
            0,
            Some(0),
            1.,
            Some(1),
            Some(DaylightSavingsConfig::NotApplicable),
            false,
            false,
            // following starts/ends are corrected from Python tests which erroneously use previous
            // "start" field instead of "start360" (which has different origin for angle)
            serde_json::from_value(json!([
                {"start360": 0, "end360": 45},
                {"start360": 45, "end360": 90},
                {"start360": 90, "end360": 135},
                {"start360": 135, "end360": 180,
                    "shading": [
                        {"type": "obstacle", "height": 10.5, "distance": 12}
                    ]
                },
                {"start360": 180, "end360": 225},
                {"start360": 225, "end360": 270},
                {"start360": 270, "end360": 315},
                {"start360": 315, "end360": 360}
            ]))
            .unwrap(),
        );
        let external_sensor: ExternalSensor = serde_json::from_value(json!({
            "correlation": [
                {"temperature": 0.0, "max_charge": 1.0},
                {"temperature": 10.0, "max_charge": 0.9},
                {"temperature": 18.0, "max_charge": 0.0}
            ]
        }))
        .unwrap();

        ChargeControl::new(
            ControlLogicType::Automatic,
            schedule_for_charge_control,
            &simulation_time_for_charge_control
                .iter()
                .current_iteration(),
            0,
            1.,
            [1.0, 0.8].into_iter().map(Some).collect(),
            Some(15.5),
            None,
            Some(external_conditions.into()),
            Some(external_sensor),
            None,
        )
        .unwrap()
    }

    #[fixture]
    fn controls_for_combination() -> IndexMap<String, Arc<Control>> {
        let cost_schedule = vec![
            5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 10.0, 10.0, 10.0, 10.0,
            10.0, 10.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
        ];
        let cost_minimising_control = Control::OnOffMinimisingTime(
            OnOffCostMinimisingTimeControl::new(
                cost_schedule,
                0,
                1.,
                5.0, // Need 12 "on" hours
            )
            .unwrap(),
        );

        IndexMap::from([
            (
                "ctrl1".into(),
                Control::OnOffTime(OnOffTimeControl::new(
                    [true, true, false, true, true, true, true, true]
                        .into_iter()
                        .map(Some)
                        .collect_vec(),
                    0,
                    1.,
                ))
                .into(),
            ),
            (
                "ctrl2".into(),
                Control::OnOffTime(OnOffTimeControl::new(
                    [false, true, true, false, false, false, true, false]
                        .into_iter()
                        .map(Some)
                        .collect_vec(),
                    0,
                    1.,
                ))
                .into(),
            ),
            (
                "ctrl3".into(),
                Control::OnOffTime(OnOffTimeControl::new(
                    [true, false, true, false, false, false, true, false]
                        .into_iter()
                        .map(Some)
                        .collect_vec(),
                    0,
                    1.,
                ))
                .into(),
            ),
            (
                "ctrl4".into(),
                Control::SetpointTime(SetpointTimeControl::new(
                    [45.0, 47.0, 50.0, 48.0, 48.0, 48.0, 48.0, 48.0]
                        .into_iter()
                        .map(Some)
                        .collect_vec(),
                    0,
                    1.,
                    Default::default(),
                    Default::default(),
                    1.,
                ))
                .into(),
            ),
            (
                "ctrl5".into(),
                Control::SetpointTime(SetpointTimeControl::new(
                    [52.0, 52.0, 52.0, 52.0, 52.0, 52.0, 52.0, 52.0]
                        .into_iter()
                        .map(Some)
                        .collect_vec(),
                    0,
                    1.,
                    Default::default(),
                    Default::default(),
                    1.,
                ))
                .into(),
            ),
            (
                "ctrl6".into(),
                Control::OnOffTime(OnOffTimeControl::new(
                    [true, true, false, true, true, true, true, true]
                        .into_iter()
                        .map(Some)
                        .collect_vec(),
                    0,
                    1.,
                ))
                .into(),
            ),
            (
                "ctrl7".into(),
                Control::OnOffTime(OnOffTimeControl::new(
                    [false, true, false, false, false, false, true, false]
                        .into_iter()
                        .map(Some)
                        .collect_vec(),
                    0,
                    1.,
                ))
                .into(),
            ),
            (
                "ctrl8".into(),
                Control::OnOffTime(OnOffTimeControl::new(
                    [true, false, false, true, true, true, true, true]
                        .into_iter()
                        .map(Some)
                        .collect_vec(),
                    0,
                    1.,
                ))
                .into(),
            ),
            (
                "ctrl9".into(),
                Control::SetpointTime(SetpointTimeControl::new(
                    vec![
                        Some(45.0),
                        None,
                        Some(50.0),
                        Some(48.0),
                        Some(48.0),
                        None,
                        Some(48.0),
                        Some(48.0),
                    ],
                    0,
                    1.,
                    Default::default(),
                    Default::default(),
                    1.,
                ))
                .into(),
            ),
            ("ctrl10".into(), cost_minimising_control.into()),
        ])
    }

    #[fixture]
    fn combination_control_on_off(
        controls_for_combination: IndexMap<String, Arc<Control>>,
    ) -> CombinationTimeControl {
        let combination_on_off: ControlCombinations = serde_json::from_value(json!({
            "main": {"operation": "AND", "controls": ["ctrl1", "ctrl2", "comb1", "comb2"]},
            "comb1": {"operation": "OR", "controls": ["ctrl3", "comb3"]},
            "comb2": {"operation": "MAX", "controls": ["ctrl4", "ctrl5"]},
            "comb3": {"operation": "XOR", "controls": ["ctrl6", "ctrl7", "ctrl8"]}
        }))
        .unwrap();

        CombinationTimeControl::new(combination_on_off, controls_for_combination).unwrap()
    }

    #[fixture]
    fn combination_control_setpoint(
        controls_for_combination: IndexMap<String, Arc<Control>>,
    ) -> CombinationTimeControl {
        let combination_setpoint: ControlCombinations = serde_json::from_value(json!({
            "main": {"operation": "AND", "controls": ["ctrl1", "ctrl2", "comb1"]},
            "comb1": {"operation": "MAX", "controls": ["ctrl4", "ctrl5"]}
        }))
        .unwrap();

        CombinationTimeControl::new(combination_setpoint, controls_for_combination).unwrap()
    }

    #[fixture]
    fn combination_control_req(
        controls_for_combination: IndexMap<String, Arc<Control>>,
    ) -> CombinationTimeControl {
        let combination_req: ControlCombinations = serde_json::from_value(json!({
            "main": {"operation": "AND", "controls": ["ctrl9", "comb1"]},
            "comb1": {"operation": "AND", "controls": ["ctrl4", "ctrl1"]}
        }))
        .unwrap();

        CombinationTimeControl::new(combination_req, controls_for_combination).unwrap()
    }

    #[fixture]
    fn combination_control_on_off_cost(
        controls_for_combination: IndexMap<String, Arc<Control>>,
    ) -> CombinationTimeControl {
        let combination_on_off_cost: ControlCombinations = serde_json::from_value(json!({
            "main": {"operation": "AND", "controls": ["ctrl1", "ctrl2", "comb1"]},
            "comb1": {"operation": "OR", "controls": ["ctrl3", "ctrl10"]}
        }))
        .unwrap();

        CombinationTimeControl::new(combination_on_off_cost, controls_for_combination).unwrap()
    }

    #[fixture]
    fn controls_for_target_charge(
        charge_control_for_combination: ChargeControl,
    ) -> IndexMap<String, Arc<Control>> {
        IndexMap::from([
            (
                "ctrl11".into(),
                Control::OnOffTime(OnOffTimeControl::new(
                    [true, false, false, true, true, true, true, true]
                        .into_iter()
                        .map(Some)
                        .collect_vec(),
                    0,
                    1.,
                ))
                .into(),
            ),
            (
                "ctrl12".into(),
                Control::Charge(charge_control_for_combination).into(),
            ),
            (
                "ctrl13".into(),
                Control::OnOffTime(OnOffTimeControl::new(
                    [true, true, false, false, true, false, true, true]
                        .into_iter()
                        .map(Some)
                        .collect_vec(),
                    0,
                    1.,
                ))
                .into(),
            ),
        ])
    }

    #[fixture]
    fn combination_control_target_charge(
        controls_for_target_charge: IndexMap<String, Arc<Control>>,
    ) -> CombinationTimeControl {
        CombinationTimeControl::new(
            serde_json::from_value(json!({
                "main": {"operation": "AND", "controls": ["ctrl11", "ctrl12"]},
            }))
            .unwrap(),
            controls_for_target_charge,
        )
        .unwrap()
    }

    #[fixture]
    fn combination_control_target_charge1(
        controls_for_target_charge: IndexMap<String, Arc<Control>>,
    ) -> CombinationTimeControl {
        CombinationTimeControl::new(
            serde_json::from_value(json!({
                "main": {"operation": "AND", "controls": ["ctrl11", "ctrl13"]},
            }))
            .unwrap(),
            controls_for_target_charge,
        )
        .unwrap()
    }

    #[fixture]
    fn simulation_time_for_combinations() -> SimulationTime {
        SimulationTime::new(0., 8., 1.)
    }

    #[rstest]
    fn test_is_on_for_combination(
        combination_control_on_off: CombinationTimeControl,
        simulation_time_for_combinations: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time_for_combinations.iter().enumerate() {
            assert_eq!(
                combination_control_on_off.is_on(&t_it),
                [false, false, false, false, false, false, true, false][t_idx]
            );
        }
    }

    #[rstest]
    fn test_setpnt_for_combination(
        combination_control_setpoint: CombinationTimeControl,
        simulation_time_for_combinations: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time_for_combinations.iter().enumerate() {
            assert_eq!(
                combination_control_setpoint.setpnt(&t_it),
                [None, Some(52.0), None, None, None, None, Some(52.0), None][t_idx]
            );
        }
    }

    #[rstest]
    fn test_in_required_period_for_combination(
        combination_control_req: CombinationTimeControl,
        simulation_time_for_combinations: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time_for_combinations.iter().enumerate() {
            assert_eq!(
                combination_control_req.in_required_period(&t_it),
                Some([true, false, false, true, true, false, true, true][t_idx]),
                "incorrect required period returned on iteration {}",
                t_idx + 1
            );
        }
    }

    #[rstest]
    fn test_is_on_cost_for_combination(
        combination_control_on_off_cost: CombinationTimeControl,
        simulation_time_for_combinations: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time_for_combinations.iter().enumerate() {
            assert_eq!(
                combination_control_on_off_cost.is_on(&t_it),
                [false, true, false, false, false, false, true, false][t_idx]
            );
        }
    }

    #[rstest]
    fn test_target_charge_for_combination(
        combination_control_target_charge: CombinationTimeControl,
        combination_control_target_charge1: CombinationTimeControl,
        simulation_time_for_combinations: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time_for_combinations.iter().enumerate() {
            assert_eq!(
                combination_control_target_charge
                    .target_charge(&t_it, None)
                    .unwrap(),
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,][t_idx]
            );
        }

        assert!(combination_control_target_charge1
            .target_charge(
                &simulation_time_for_combinations.iter().next().unwrap(),
                None
            )
            .is_err());
    }

    #[fixture]
    fn controls_for_invalid_combinations(
        charge_control_for_combination: ChargeControl,
    ) -> IndexMap<String, Arc<Control>> {
        IndexMap::from([
            (
                "ctrl14".into(),
                Control::OnOffTime(OnOffTimeControl::new(
                    [true, false, false, true, true, true, true, true]
                        .into_iter()
                        .map(Some)
                        .collect_vec(),
                    0,
                    1.,
                ))
                .into(),
            ),
            (
                "ctrl15".into(),
                Control::Charge(charge_control_for_combination).into(),
            ),
            (
                "ctrl16".into(),
                Control::OnOffTime(OnOffTimeControl::new(
                    [true, true, false, false, true, false, true, true]
                        .into_iter()
                        .map(Some)
                        .collect_vec(),
                    0,
                    1.,
                ))
                .into(),
            ),
            (
                "ctrl17".into(),
                Control::OnOffTime(OnOffTimeControl::new(
                    [true, true, false, false, true, false, true, true]
                        .into_iter()
                        .map(Some)
                        .collect_vec(),
                    0,
                    1.,
                ))
                .into(),
            ),
        ])
    }

    // this test is introduced in the Rust to test up-front validation of combinations
    #[rstest]
    fn test_invalid_combinations_caught_on_instantiation(
        controls_for_invalid_combinations: IndexMap<String, Arc<Control>>,
    ) {
        let invalid_combinations = [
            json!({
                "main": {"operation": "AND", "controls": ["ctrl15"]},
            }), // only one control referenced, regardless of operation
        ];

        for invalid_combination in invalid_combinations {
            assert!(CombinationTimeControl::new(
                serde_json::from_value(invalid_combination).unwrap(),
                controls_for_invalid_combinations.clone(),
            )
            .is_err());
        }
    }
}
