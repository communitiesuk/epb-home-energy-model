// This module provides structs to model time controls

use crate::core::units::{HOURS_PER_DAY, WATTS_PER_KILOWATT};
use crate::external_conditions::ExternalConditions;
use crate::input::{
    ControlCombination, ControlCombinationOperation, ControlCombinations, ControlLogicType,
    ExternalSensor, ExternalSensorCorrelation, HeatSourceControlType, SmartApplianceBattery,
    MAIN_REFERENCE,
};
use crate::simulation_time::{SimulationTimeIteration, SimulationTimeIterator, HOURS_IN_DAY};
use anyhow::{anyhow, bail};
use atomic_float::AtomicF64;
use bounded_vec_deque::BoundedVecDeque;
use indexmap::IndexMap;
use itertools::Itertools;
use parking_lot::RwLock;
use std::collections::VecDeque;
use std::fmt::{Debug, Formatter};
use std::iter::repeat;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

#[derive(Debug)]
pub(crate) enum Control {
    OnOffTime(OnOffTimeControl),
    Charge(ChargeControl),
    OnOffMinimisingTime(OnOffMinimisingTimeControl),
    SetpointTime(SetpointTimeControl),
    SmartAppliance(SmartApplianceControl),
    CombinationTime(CombinationTimeControl),
}

// macro so accessing individual controls through the enum isn't so repetitive
macro_rules! per_control {
    ($val:expr, $pattern:pat => { $res:expr }) => {
        match $val {
            Control::OnOffTime($pattern) => $res,
            Control::Charge($pattern) => $res,
            Control::OnOffMinimisingTime($pattern) => $res,
            Control::SetpointTime($pattern) => $res,
            Control::SmartAppliance($pattern) => $res,
            Control::CombinationTime($pattern) => $res,
        }
    };
}

use crate::compare_floats::min_of_2;
use crate::core::energy_supply::energy_supply::EnergySupply;
pub(crate) use per_control;

pub trait ControlBehaviour: Send + Sync {
    fn in_required_period(
        &self,
        _simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool> {
        None
    }

    fn setpnt(&self, _simulation_time_iteration: &SimulationTimeIteration) -> Option<f64> {
        None
    }
}

impl Debug for dyn ControlBehaviour {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        // if we can downcast self to e.g. Control (if it is one), which we know is Debug, this would be better
        write!(f, "A control object")
    }
}

impl Control {
    pub fn is_on(&self, simulation_time_iteration: SimulationTimeIteration) -> bool {
        per_control!(self, c => {c.is_on(&simulation_time_iteration)})
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
}

#[derive(Debug)]
pub(crate) enum HeatSourceControl {
    HotWaterTimer(Arc<Control>),
    WindowOpening(Arc<Control>),
    WindowOpeningLivingRoom(Arc<Control>),
    WindowOpeningRestOfDwelling(Arc<Control>),
    AlwaysOff(Arc<Control>),
}

impl HeatSourceControl {
    pub fn has_type(&self, control_type: HeatSourceControlType) -> bool {
        match control_type {
            HeatSourceControlType::HotWaterTimer => {
                matches!(self, HeatSourceControl::HotWaterTimer(_))
            }
            HeatSourceControlType::WindowOpening => {
                matches!(self, HeatSourceControl::WindowOpening(_))
            }
            HeatSourceControlType::WindowOpeningLivingRoom => {
                matches!(self, HeatSourceControl::WindowOpeningLivingRoom(_))
            }
            HeatSourceControlType::WindowOpeningRestOfDwelling => {
                matches!(self, HeatSourceControl::WindowOpeningRestOfDwelling(_))
            }
            HeatSourceControlType::AlwaysOff => matches!(self, HeatSourceControl::AlwaysOff(_)),
        }
    }

    pub fn get(&self) -> Arc<Control> {
        match self {
            HeatSourceControl::HotWaterTimer(control) => control.clone(),
            HeatSourceControl::WindowOpening(control) => control.clone(),
            HeatSourceControl::WindowOpeningLivingRoom(control) => control.clone(),
            HeatSourceControl::WindowOpeningRestOfDwelling(control) => control.clone(),
            HeatSourceControl::AlwaysOff(control) => control.clone(),
        }
    }
}

/// An object to model a time-only control with on/off (not modulating) operation
#[derive(Clone, Debug)]
pub struct OnOffTimeControl {
    /// list of boolean values where true means "on" (one entry per hour)
    schedule: Vec<Option<bool>>,
    start_day: u32,
    time_series_step: f64,
}

impl OnOffTimeControl {
    pub fn new(
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

    pub(crate) fn is_on(&self, timestep: &SimulationTimeIteration) -> bool {
        self.schedule[timestep.time_series_idx(self.start_day, self.time_series_step)]
            .unwrap_or(false)
    }
}

impl ControlBehaviour for OnOffTimeControl {}

/// An object to model a control that governs electrical charging of a heat storage device
/// that can respond to signals from the grid, for example when carbon intensity is low
#[derive(Debug)]
pub(crate) struct ChargeControl {
    logic_type: ControlLogicType,
    /// list of boolean values where true means "on" (one entry per hour)
    pub schedule: Vec<bool>,
    pub start_day: u32,
    pub time_series_step: f64,
    /// Proportion of the charge targeted for each day
    pub charge_level: Vec<Option<f64>>,
    temp_charge_cut: Option<f64>,
    temp_charge_cut_delta: Option<Vec<f64>>,
    _min_target_charge_factor: Option<f64>,
    _full_charge_temp_diff: Option<f64>,
    external_conditions: Arc<ExternalConditions>,
    external_sensor: Option<ExternalSensor>,
    hhrsh: Option<ChargeControlHhrshFields>,
}

#[derive(Debug)]
pub(crate) struct ChargeControlHhrshFields {
    steps_day: usize,
    demand: Arc<RwLock<BoundedVecDeque<Option<f64>>>>,
    past_ext_temp: Arc<RwLock<BoundedVecDeque<Option<f64>>>>,
    future_ext_temp: Arc<RwLock<BoundedVecDeque<Option<f64>>>>,
    energy_to_store: AtomicF64,
}

impl ChargeControl {
    pub(crate) fn new(
        logic_type: ControlLogicType,
        schedule: Vec<bool>,
        simulation_timestep: f64,
        start_day: u32,
        time_series_step: f64,
        charge_level: Vec<Option<f64>>,
        temp_charge_cut: Option<f64>,
        temp_charge_cut_delta: Option<Vec<f64>>,
        min_target_charge_factor: Option<f64>,
        full_charge_temp_diff: Option<f64>,
        external_conditions: Arc<ExternalConditions>,
        external_sensor: Option<ExternalSensor>,
    ) -> anyhow::Result<Self> {
        let hhrsh_fields: Option<ChargeControlHhrshFields> = match logic_type {
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
                    bail!("CElect ChargeControl definition is missing input parameters.");
                }
                None
            }
            ControlLogicType::Hhrsh => {
                if temp_charge_cut.is_none() {
                    bail!("Hhrsh ChargeControl definition is missing input parameters.");
                }

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
                let energy_to_store = AtomicF64::new(0.0);
                Some(ChargeControlHhrshFields {
                    steps_day,
                    demand,
                    past_ext_temp,
                    future_ext_temp,
                    energy_to_store,
                })

                // TODO (from Python) Consider adding solar data for HHRSH logic in addition to heating degree hours.
            }
        };

        Ok(Self {
            logic_type,
            schedule,
            start_day,
            time_series_step,
            charge_level,
            temp_charge_cut,
            temp_charge_cut_delta,
            _min_target_charge_factor: min_target_charge_factor,
            _full_charge_temp_diff: full_charge_temp_diff,
            external_conditions,
            external_sensor,
            hhrsh: hhrsh_fields,
        })
    }

    pub(crate) fn logic_type(&self) -> ControlLogicType {
        self.logic_type
    }

    pub(crate) fn is_on(&self, iteration: &SimulationTimeIteration) -> bool {
        self.schedule[iteration.time_series_idx(self.start_day, self.time_series_step)]
    }

    /// Return the charge level value from the list given in inputs; one value per day
    pub(crate) fn target_charge(
        &self,
        simtime: SimulationTimeIteration,
        temp_air: Option<f64>,
    ) -> anyhow::Result<f64> {
        // Calculate target charge nominal when unit is on
        let mut target_charge_nominal = if self.is_on(&simtime) {
            self.charge_level[simtime.time_series_idx_days(self.start_day, self.time_series_step)]
                .unwrap_or_default()
        } else {
            // If unit is off send 0.0 for target charge
            0.
        };

        Ok(match self.logic_type {
            ControlLogicType::Manual => target_charge_nominal,
            _ => {
                // automatic, celect and hhrsh control include temperature charge cut logic
                let temp_charge_cut = self.temp_charge_cut_corr(simtime)?;

                if temp_air.is_some_and(|temp_air| temp_air >= temp_charge_cut) {
                    // Control logic cut when temp_air is over temp_charge cut
                    target_charge_nominal = 0.;
                }
                match self.logic_type {
                    ControlLogicType::Automatic => {
                        // Automatic charge control can be achieved using internal thermostat(s) to
                        // control the extent of charging of the heaters. All or nothing approach

                        // Controls can also be supplemented by an external weather sensor,
                        // which tends to act as a limiting device to prevent the storage heaters from overcharging.
                        if self.external_sensor.is_some() {
                            let limit =
                                self.get_limit_factor(self.external_conditions.air_temp(&simtime))?;
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
                        if self.external_sensor.is_some() {
                            let limit =
                                self.get_limit_factor(self.external_conditions.air_temp(&simtime))?;
                            target_charge_nominal * limit
                        } else {
                            target_charge_nominal
                        }
                    }
                    ControlLogicType::Hhrsh => {
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
        })
    }

    pub(crate) fn energy_to_store(
        &self,
        energy_demand: f64,
        base_temp: f64,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        // ugly, but this method cannot be called when control does not have HHRSH logic type
        if !matches!(self.logic_type, ControlLogicType::Hhrsh) {
            unreachable!("energy_to_store() should not be called when control does not have HHRSH logic type.");
        }
        let ChargeControlHhrshFields {
            steps_day,
            demand,
            past_ext_temp,
            future_ext_temp,
            energy_to_store: energy_to_store_atomic,
        } = self.hhrsh.as_ref().expect("HHRSH fields should be set.");
        demand.write().push_front(Some(energy_demand));
        future_ext_temp.write().push_front(Some(
            self.external_conditions
                .air_temp_with_offset(&simtime, *steps_day),
        ));
        past_ext_temp
            .write()
            .push_front(Some(self.external_conditions.air_temp(&simtime)));

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

    /// Correct nominal/json temp_charge_cut with monthly table
    pub(crate) fn temp_charge_cut_corr(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let temp_charge_cut_delta = if let Some(temp_charge_cut_delta) =
            self.temp_charge_cut_delta.as_ref()
        {
            temp_charge_cut_delta[simtime.time_series_idx(self.start_day, self.time_series_step)]
        } else {
            0.0
        };

        Ok(self.temp_charge_cut.ok_or_else(|| {
            anyhow!("temp_charge_cut in this ChargeControl was expected to be set.")
        })? + temp_charge_cut_delta)
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

impl ControlBehaviour for ChargeControl {}

#[derive(Clone, Debug)]
pub struct OnOffMinimisingTimeControl {
    /// list of boolean values where true means "on" (one entry per hour)
    schedule: Vec<bool>,
    start_day: u32,
    time_series_step: f64,
}

impl OnOffMinimisingTimeControl {
    pub fn new(
        schedule: Vec<f64>,
        start_day: u32,
        time_series_step: f64,
        time_on_daily: f64,
    ) -> Self {
        let timesteps_per_day = (HOURS_IN_DAY as f64 / time_series_step) as usize;
        let timesteps_on_daily = (time_on_daily / time_series_step) as usize;
        let time_series_len_days =
            ((schedule.len() as f64 * time_series_step / HOURS_IN_DAY as f64).ceil()) as usize;

        let mut built_schedule: Vec<bool> = vec![];
        for day in 0..time_series_len_days {
            // Get part of the schedule for current day
            let schedule_day_start = day * timesteps_per_day;
            let schedule_day_end = schedule_day_start + timesteps_per_day;
            let schedule_day = schedule[schedule_day_start..schedule_day_end].to_vec();

            // Find required number of timesteps with lowest costs
            let mut schedule_day_cost_lowest = schedule_day.clone();
            schedule_day_cost_lowest.sort_by(f64::total_cmp);
            let schedule_day_cost_lowest = schedule_day_cost_lowest[0..timesteps_on_daily].to_vec();

            // Initialise boolean schedule for day
            let mut schedule_onoff_day = vec![false; timesteps_per_day];

            assert_eq!(schedule_onoff_day.len(), schedule_day.len());
            let mut timesteps_to_be_allocated = timesteps_on_daily;
            for cost in schedule_day_cost_lowest.iter() {
                for (idx, entry) in schedule_day.iter().enumerate() {
                    if timesteps_to_be_allocated < 1 {
                        break;
                    }
                    if entry == cost {
                        schedule_onoff_day[idx] = true;
                        timesteps_to_be_allocated -= 1;
                    }
                }
            }
            built_schedule.extend(schedule_onoff_day);
        }

        OnOffMinimisingTimeControl {
            schedule: built_schedule,
            start_day,
            time_series_step,
        }
    }

    pub fn is_on(&self, timestep: &SimulationTimeIteration) -> bool {
        self.schedule[timestep.time_series_idx(self.start_day, self.time_series_step)]
    }
}

impl ControlBehaviour for OnOffMinimisingTimeControl {}

#[derive(Clone, Debug)]
pub struct SetpointTimeControl {
    /// list of float values (one entry per hour)
    schedule: Vec<Option<f64>>,
    /// first day of the time series, day of the year, 0 to 365 (single value)
    start_day: u32,
    /// timestep of the time series data, in hours
    time_series_step: f64,
    /// min setpoint allowed
    setpoint_min: Option<f64>,
    /// max setpoint allowed
    setpoint_max: Option<f64>,
    /// if both min and max limits are set but setpoint isn't,
    /// whether to default to min (false) or max (true)
    default_to_max: bool,
    /// how long before heating period the system
    /// should switch on, in hours
    timesteps_advstart: u32,
}

/// Return true if current time is inside specified time for heating/cooling
///
/// (not including timesteps where system is only on due to min or max
/// setpoint or advanced start)
impl SetpointTimeControl {
    pub fn new(
        schedule: Vec<Option<f64>>,
        start_day: u32,
        time_series_step: f64,
        setpoint_min: Option<f64>,
        setpoint_max: Option<f64>,
        default_to_max: Option<bool>,
        duration_advanced_start: f64,
        timestep: f64,
    ) -> Result<Self, &'static str> {
        if setpoint_min.is_some() && setpoint_max.is_some() && default_to_max.is_none() {
            return Err(
                "default_to_max should be set when both setpoint_min and setpoint_max are set",
            );
        }

        Ok(SetpointTimeControl {
            schedule,
            start_day,
            time_series_step,
            setpoint_min,
            setpoint_max,
            default_to_max: default_to_max.unwrap_or(false),
            timesteps_advstart: (duration_advanced_start / timestep).round() as u32,
        })
    }

    pub fn is_on(&self, timestep: &SimulationTimeIteration) -> bool {
        let schedule_idx = timestep.time_series_idx(self.start_day, self.time_series_step);

        self.is_on_for_timestep_idx(schedule_idx)
    }

    fn is_on_for_timestep_idx(&self, schedule_idx: usize) -> bool {
        let setpnt = self.schedule[schedule_idx];

        if setpnt.is_none() {
            for timesteps_ahead in 1..(self.timesteps_advstart + 1) {
                let timesteps_ahead = timesteps_ahead as usize;
                if self.schedule.len() <= (schedule_idx + timesteps_ahead) {
                    break;
                }
                if self.schedule[schedule_idx + timesteps_ahead].is_some() {
                    return true;
                }
            }
        }

        !(setpnt.is_none() && self.setpoint_min.is_none() && self.setpoint_max.is_none())
    }
}

impl ControlBehaviour for SetpointTimeControl {
    fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool> {
        let schedule_idx =
            simulation_time_iteration.time_series_idx(self.start_day, self.time_series_step);

        Some(self.schedule[schedule_idx].is_some())
    }

    fn setpnt(&self, simulation_time_iteration: &SimulationTimeIteration) -> Option<f64> {
        let schedule_idx =
            simulation_time_iteration.time_series_idx(self.start_day, self.time_series_step);

        let mut setpnt = self.schedule[schedule_idx];

        if setpnt.is_none() {
            for timesteps_ahead in 1..(self.timesteps_advstart + 1) {
                let timesteps_ahead = timesteps_ahead as usize;
                if self.schedule.len() <= (schedule_idx + timesteps_ahead) {
                    break;
                }
                if let Some(s) = self.schedule[schedule_idx + timesteps_ahead] {
                    setpnt = Some(s);
                }
            }
        }

        match (setpnt, self.setpoint_min, self.setpoint_max) {
            (None, None, None) => None,
            (None, None, Some(max)) => Some(max),
            (None, Some(min), None) => Some(min),
            (None, Some(min), Some(max)) => match self.default_to_max {
                true => Some(max),
                false => Some(min),
            },
            (Some(_s), _, _) => {
                if self.setpoint_max.is_some() {
                    setpnt = match self.setpoint_max < setpnt {
                        true => self.setpoint_max,
                        false => setpnt,
                    }
                }
                if self.setpoint_min.is_some() {
                    setpnt = match self.setpoint_min > setpnt {
                        true => self.setpoint_min,
                        false => setpnt,
                    }
                }
                setpnt
            }
        }
    }
}

/// An object for managing loadshifting appliances
#[derive(Debug)]
pub(crate) struct SmartApplianceControl {
    appliance_names: Vec<String>,
    energy_supplies: IndexMap<String, Arc<RwLock<EnergySupply>>>,
    battery_states_of_charge: IndexMap<String, Vec<AtomicF64>>,
    ts_power: IndexMap<String, Vec<AtomicF64>>,
    _ts_weight: IndexMap<String, Vec<f64>>,
    ts_step: f64,
    simulation_timestep: f64,
    ts_step_ratio: f64,
    non_appliance_demand_24hr: IndexMap<String, Vec<AtomicF64>>,
    buffer_length: usize,
}

impl SmartApplianceControl {
    /// Arguments:
    /// * `power_timeseries`  - dictionary of lists containing expected power for appliances
    ///                         for each energy supply, for the entire length of the simulation
    /// * `weight_timeseries` - dictionary of lists containing demand weight for each
    ///                         energy supply, for the entire length of the simulation
    /// * `timeseries_step`   - timestep of weight and demand timeseries
    ///                         (not necessarily equal to simulation_time.timestep())
    /// * `simulation_time`   - reference to a SimulationTime object
    /// * `non_appliance_demand_24hr` - dictionary of lists containing 24 hour buffers of
    ///                                 demand per end user for each energy supply
    /// * `battery_24hr`      - dictionary of lists containing 24 hour buffers of
    ///                         battery state of charge for each energy supply
    /// * `energy_supplies`   - dictionary of energy supply objects in the simulation
    /// * `appliance_names`   - list of names of appliances
    pub(crate) fn new(
        power_timeseries: &IndexMap<String, Vec<f64>>,
        weight_timeseries: &IndexMap<String, Vec<f64>>,
        timeseries_step: f64,
        simulation_time_iterator: &SimulationTimeIterator,
        non_appliance_demand_24hr: Option<IndexMap<String, Vec<f64>>>,
        battery_24hr: Option<&SmartApplianceBattery>,
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
                    battery_24hr
                        .as_ref()
                        .expect(
                            "Expected battery 24 h to be set for energy supply that has battery ",
                        )
                        .battery_state_of_charge[name]
                        .iter()
                        .map(|x| AtomicF64::new(*x))
                        .collect_vec(),
                )
            })
            .collect();
        for energy_supply in energy_supplies.keys() {
            if power_timeseries[energy_supply].len() != weight_timeseries[energy_supply].len() {
                bail!("ERROR: loadshifting weight and power timeseries not equal in length")
            }
            if power_timeseries[energy_supply].len() as f64 * timeseries_step
                < simulation_time_iterator.total_steps() as f64
                    * simulation_time_iterator.step_in_hours()
            {
                bail!("ERROR: loadshifting weight+power timeseries shorter than simulation length")
            }
        }
        Ok(Self {
            appliance_names,
            energy_supplies,
            battery_states_of_charge,
            ts_power: power_timeseries.iter().map(|(name, series)| (name.to_owned(), series.iter().map(|x| AtomicF64::new(*x)).collect_vec())).collect(),
            _ts_weight: weight_timeseries.clone(),
            ts_step: timeseries_step,
            simulation_timestep: simulation_time_iterator.step_in_hours(),
            ts_step_ratio: simulation_time_iterator.step_in_hours() / timeseries_step,
            buffer_length: non_appliance_demand_24hr.clone().expect("Expected non_appliance_demand_24hr to be set").first().as_ref().ok_or_else(|| anyhow!("non_appliance_demand_24hr parameter for SmartApplianceControl cannot be empty."))?.1.len(),
            non_appliance_demand_24hr: non_appliance_demand_24hr.clone()
                .expect("Expected non_appliance_demand_24hr to be set")
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

    /// Controls generally have this method, though this control does not in Python implementation.
    /// Defaulting to returning true for a fallback implementation here.
    pub(crate) fn is_on(&self, _simtime: &SimulationTimeIteration) -> bool {
        true
    }
}

impl ControlBehaviour for SmartApplianceControl {}

#[derive(Debug)]
pub(crate) struct CombinationTimeControl {
    combinations: ControlCombinations,
    controls: IndexMap<String, Arc<Control>>,
    start_day: u32,
    time_series_step: f64,
}

impl CombinationTimeControl {
    pub(crate) fn new(
        combinations: ControlCombinations,
        controls: IndexMap<String, Arc<Control>>,
        start_day: u32,
        time_series_step: f64,
    ) -> anyhow::Result<Self> {
        Self::validate_combinations(&combinations, &controls)?;

        Ok(Self {
            combinations,
            controls,
            start_day,
            time_series_step,
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
        control.is_on(simtime)
    }

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
                self.evaluate_combination_is_on(control, simtime)
            } else {
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
            | ControlCombinationOperation::Mean => results.all(|x| x),
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
            | c @ Control::OnOffMinimisingTime(_) => c.is_on(simtime),
            Control::SetpointTime(c) => c
                .in_required_period(&simtime)
                .expect("SetpointTimeControl in_required_period() method will always return Some"),
            _ => unreachable!("CombinationTimeControl only combined OnOffTime, Charge, OnOffMinimisingTime or SetpointTime controls"),
        }
    }

    /// This function processes a combination of control elements,
    /// applying boolean logic (AND, OR, XOR, etc.) to their evaluation results.
    /// It checks the type of controls, validates allowed
    /// combinations and returns the evaluation result based on the specified operation.
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
                _ => bail!("Invalid combination of controls encountered"),
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
            | c @ Control::OnOffMinimisingTime(_) => SetpointOrBoolean::Boolean(c.is_on(simtime)),
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

    pub(crate) fn is_on(&self, simtime: &SimulationTimeIteration) -> bool {
        self.evaluate_combination_is_on(MAIN_REFERENCE, *simtime)
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

    const ON_OFF_SCHEDULE: [bool; 8] = [true, false, true, true, false, true, false, false];

    #[fixture]
    pub fn simulation_time() -> SimulationTimeIterator {
        SimulationTime::new(0.0, 8.0, 1.0).iter()
    }

    #[fixture]
    pub fn on_off_time_control() -> OnOffTimeControl {
        OnOffTimeControl::new(
            ON_OFF_SCHEDULE.iter().map(|&v| Some(v)).collect_vec(),
            0,
            1.0,
        )
    }

    #[rstest]
    pub fn should_be_on_for_on_off_time_control(
        on_off_time_control: OnOffTimeControl,
        simulation_time: SimulationTimeIterator,
    ) {
        for it in simulation_time {
            assert_eq!(on_off_time_control.is_on(&it), ON_OFF_SCHEDULE[it.index]);
        }
    }

    #[fixture]
    pub fn on_off_minimising_control() -> OnOffMinimisingTimeControl {
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
        OnOffMinimisingTimeControl::new(schedule, 0, 1.0, 12.0)
    }

    #[rstest]
    pub fn should_be_on_for_cost_minimising_control(
        on_off_minimising_control: OnOffMinimisingTimeControl,
        simulation_time: SimulationTimeIterator,
    ) {
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
        for it in simulation_time {
            assert_eq!(
                on_off_minimising_control.is_on(&it),
                resulting_schedule[it.index]
            );
        }
    }

    #[fixture]
    pub fn setpoint_schedule() -> Vec<Option<f64>> {
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

    #[fixture]
    pub fn setpoint_time_control(
        setpoint_schedule: Vec<Option<f64>>,
        simulation_time: SimulationTimeIterator,
    ) -> SetpointTimeControl {
        SetpointTimeControl::new(
            setpoint_schedule,
            0,
            1.0,
            None,
            None,
            None,
            Default::default(),
            simulation_time.step_in_hours(),
        )
        .unwrap()
    }

    #[fixture]
    pub fn setpoint_time_control_min(
        setpoint_schedule: Vec<Option<f64>>,
        simulation_time: SimulationTimeIterator,
    ) -> SetpointTimeControl {
        SetpointTimeControl::new(
            setpoint_schedule,
            0,
            1.0,
            Some(16.0),
            None,
            None,
            Default::default(),
            simulation_time.step_in_hours(),
        )
        .unwrap()
    }

    #[fixture]
    pub fn setpoint_time_control_max(
        setpoint_schedule: Vec<Option<f64>>,
        simulation_time: SimulationTimeIterator,
    ) -> SetpointTimeControl {
        SetpointTimeControl::new(
            setpoint_schedule,
            0,
            1.0,
            None,
            Some(24.0),
            None,
            0.0,
            simulation_time.step_in_hours(),
        )
        .unwrap()
    }

    #[fixture]
    pub fn setpoint_time_control_minmax(
        setpoint_schedule: Vec<Option<f64>>,
        simulation_time: SimulationTimeIterator,
    ) -> SetpointTimeControl {
        SetpointTimeControl::new(
            setpoint_schedule,
            0,
            1.0,
            Some(16.0),
            Some(24.0),
            Some(false),
            Default::default(),
            simulation_time.step_in_hours(),
        )
        .unwrap()
    }

    #[fixture]
    pub fn setpoint_time_control_advstart(
        setpoint_schedule: Vec<Option<f64>>,
        simulation_time: SimulationTimeIterator,
    ) -> SetpointTimeControl {
        SetpointTimeControl::new(
            setpoint_schedule,
            0,
            1.0,
            None,
            None,
            Some(false),
            1.0,
            simulation_time.step_in_hours(),
        )
        .unwrap()
    }

    #[fixture]
    pub fn setpoint_time_control_advstart_minmax(
        setpoint_schedule: Vec<Option<f64>>,
        simulation_time: SimulationTimeIterator,
    ) -> SetpointTimeControl {
        SetpointTimeControl::new(
            setpoint_schedule,
            0,
            1.0,
            Some(16.0),
            Some(24.0),
            Some(false),
            1.0,
            simulation_time.step_in_hours(),
        )
        .unwrap()
    }

    #[rstest]
    pub fn should_be_in_required_time_for_setpoint_control(
        setpoint_time_control: SetpointTimeControl,
        setpoint_time_control_min: SetpointTimeControl,
        setpoint_time_control_max: SetpointTimeControl,
        setpoint_time_control_minmax: SetpointTimeControl,
        setpoint_time_control_advstart: SetpointTimeControl,
        setpoint_time_control_advstart_minmax: SetpointTimeControl,
        simulation_time: SimulationTimeIterator,
    ) {
        let results: [bool; 8] = [true, false, false, true, false, true, true, true];
        for it in simulation_time {
            assert_eq!(
                setpoint_time_control.in_required_period(&it).unwrap(),
                results[it.index],
                "incorrect in_required_period value returned for control with no min or max set, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_min.in_required_period(&it).unwrap(),
                results[it.index],
                "incorrect in_required_period value returned for control with min set, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_max.in_required_period(&it).unwrap(),
                results[it.index],
                "incorrect in_required_period value returned for control with max set, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_minmax.in_required_period(&it).unwrap(),
                results[it.index],
                "incorrect in_required_period value returned for control with min and max set, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_advstart.in_required_period(&it).unwrap(),
                results[it.index],
                "incorrect in_required_period value returned for control with advanced start, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_advstart_minmax.in_required_period(&it).unwrap(),
                results[it.index],
                "incorrect in_required_period value returned for control with advanced start and min/max, iteration {}",
                it.index + 1
            );
        }
    }

    #[rstest]
    pub fn should_be_on_for_setpoint_control(
        setpoint_time_control: SetpointTimeControl,
        setpoint_time_control_min: SetpointTimeControl,
        setpoint_time_control_max: SetpointTimeControl,
        setpoint_time_control_minmax: SetpointTimeControl,
        setpoint_time_control_advstart: SetpointTimeControl,
        setpoint_time_control_advstart_minmax: SetpointTimeControl,
        simulation_time: SimulationTimeIterator,
    ) {
        for it in simulation_time {
            assert_eq!(
                setpoint_time_control.is_on(&it),
                [true, false, false, true, false, true, true, true][it.index],
                "incorrect is_on value returned for control with no min or max set, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_min.is_on(&it),
                true, // Should always be true for this type of control
                "incorrect is_on value returned for control with min set, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_max.is_on(&it),
                true, // Should always be true for this type of control
                "incorrect is_on value returned for control with max set, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_minmax.is_on(&it),
                true, // Should always be true for this type of control
                "incorrect is_on value returned for control with min and max set, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_advstart.is_on(&it),
                [true, false, true, true, true, true, true, true][it.index],
                "incorrect is_on value returned for control with advanced start, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_advstart_minmax.is_on(&it),
                true,
                "incorrect is_on value returned for control with advanced start and min/max, iteration {}",
                it.index + 1
            );
        }
    }

    #[rstest]
    pub fn should_have_correct_setpnt(
        setpoint_time_control: SetpointTimeControl,
        setpoint_time_control_min: SetpointTimeControl,
        setpoint_time_control_max: SetpointTimeControl,
        setpoint_time_control_minmax: SetpointTimeControl,
        setpoint_time_control_advstart: SetpointTimeControl,
        setpoint_time_control_advstart_minmax: SetpointTimeControl,
        simulation_time: SimulationTimeIterator,
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
        for it in simulation_time {
            assert_eq!(
                setpoint_time_control.setpnt(&it),
                setpoint_schedule()[it.index],
                "incorrect schedule returned for control with no min or max set, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_min.setpnt(&it),
                results_min[it.index],
                "incorrect schedule returned for control with min set, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_max.setpnt(&it),
                results_max[it.index],
                "incorrect schedule returned for control with max set, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_minmax.setpnt(&it),
                results_minmax[it.index],
                "incorrect schedule returned for control with min and max set, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_advstart.setpnt(&it),
                results_advstart[it.index],
                "incorrect schedule returned for control with advanced start, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_advstart_minmax.setpnt(&it),
                results_advstart_minmax[it.index],
                "incorrect schedule returned for control with advanced start and min and max set, iteration {}",
                it.index + 1
            );
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

    #[fixture]
    fn charge_control(
        simulation_time_for_charge_control: SimulationTime,
        schedule_for_charge_control: Vec<bool>,
    ) -> ChargeControl {
        let external_conditions = ExternalConditions::new(
            &simulation_time_for_charge_control.iter(),
            vec![
                19.0, 0.0, 1.0, 2.0, 5.0, 7.0, 6.0, 12.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0,
                19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0,
            ],
            vec![
                3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1, 3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1,
                3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1,
            ],
            vec![
                300., 250., 220., 180., 150., 120., 100., 80., 60., 40., 20., 10., 50., 100., 140.,
                190., 200., 320., 330., 340., 350., 355., 315., 5.,
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
                {"number": 1, "start360": 0, "end360": 45},
                {"number": 2, "start360": 45, "end360": 90},
                {"number": 3, "start360": 90, "end360": 135},
                {"number": 4, "start360": 135, "end360": 180,
                    "shading": [
                        {"type": "obstacle", "height": 10.5, "distance": 12}
                    ]
                },
                {"number": 5, "start360": 180, "end360": 225},
                {"number": 6, "start360": 225, "end360": 270},
                {"number": 7, "start360": 270, "end360": 315},
                {"number": 8, "start360": 315, "end360": 360}
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
            simulation_time_for_charge_control.step,
            0,
            1.,
            vec![Some(1.0), Some(0.8)],
            Some(15.5),
            None,
            None,
            None,
            external_conditions.into(),
            Some(external_sensor),
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
                {"number": 1, "start360": 0, "end360": 45},
                {"number": 2, "start360": 45, "end360": 90},
                {"number": 3, "start360": 90, "end360": 135},
                {"number": 4, "start360": 135, "end360": 180,
                    "shading": [
                        {"type": "obstacle", "height": 10.5, "distance": 12}
                    ]
                },
                {"number": 5, "start360": 180, "end360": 225},
                {"number": 6, "start360": 225, "end360": 270},
                {"number": 7, "start360": 270, "end360": 315},
                {"number": 8, "start360": 315, "end360": 360}
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
            simulation_time_for_charge_control.step,
            0,
            1.,
            [1.0, 0.8].into_iter().map(Some).collect(),
            Some(15.5),
            None,
            None,
            None,
            external_conditions.into(),
            Some(external_sensor),
        )
        .unwrap()
    }

    #[fixture]
    fn controls_for_combination() -> IndexMap<String, Arc<Control>> {
        let cost_schedule = vec![
            5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 10.0, 10.0, 10.0, 10.0,
            10.0, 10.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
        ];
        let cost_minimising_control =
            Control::OnOffMinimisingTime(OnOffMinimisingTimeControl::new(
                cost_schedule,
                0,
                1.,
                5.0, // Need 12 "on" hours
            ));

        IndexMap::from([
            (
                "ctrl1".to_string(),
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
                "ctrl2".to_string(),
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
                "ctrl3".to_string(),
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
                "ctrl4".to_string(),
                Control::SetpointTime(
                    SetpointTimeControl::new(
                        [45.0, 47.0, 50.0, 48.0, 48.0, 48.0, 48.0, 48.0]
                            .into_iter()
                            .map(Some)
                            .collect_vec(),
                        0,
                        1.,
                        None,
                        None,
                        None,
                        Default::default(),
                        1.,
                    )
                    .unwrap(),
                )
                .into(),
            ),
            (
                "ctrl5".to_string(),
                Control::SetpointTime(
                    SetpointTimeControl::new(
                        [52.0, 52.0, 52.0, 52.0, 52.0, 52.0, 52.0, 52.0]
                            .into_iter()
                            .map(Some)
                            .collect_vec(),
                        0,
                        1.,
                        None,
                        None,
                        None,
                        Default::default(),
                        1.,
                    )
                    .unwrap(),
                )
                .into(),
            ),
            (
                "ctrl6".to_string(),
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
                "ctrl7".to_string(),
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
                "ctrl8".to_string(),
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
                "ctrl9".to_string(),
                Control::SetpointTime(
                    SetpointTimeControl::new(
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
                        None,
                        None,
                        None,
                        Default::default(),
                        1.,
                    )
                    .unwrap(),
                )
                .into(),
            ),
            ("ctrl10".to_string(), cost_minimising_control.into()),
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

        CombinationTimeControl::new(combination_on_off, controls_for_combination, 0, 1.).unwrap()
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

        CombinationTimeControl::new(combination_setpoint, controls_for_combination, 0, 1.).unwrap()
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

        CombinationTimeControl::new(combination_req, controls_for_combination, 0, 1.).unwrap()
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

        CombinationTimeControl::new(combination_on_off_cost, controls_for_combination, 0, 1.)
            .unwrap()
    }

    #[fixture]
    fn controls_for_target_charge(
        charge_control_for_combination: ChargeControl,
    ) -> IndexMap<String, Arc<Control>> {
        IndexMap::from([
            (
                "ctrl11".to_string(),
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
                "ctrl12".to_string(),
                Control::Charge(charge_control_for_combination).into(),
            ),
            (
                "ctrl13".to_string(),
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
            0,
            1.,
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
            0,
            1.,
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
                "ctrl14".to_string(),
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
                "ctrl15".to_string(),
                Control::Charge(charge_control_for_combination).into(),
            ),
            (
                "ctrl16".to_string(),
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
                "ctrl17".to_string(),
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
                0,
                1.,
            )
            .is_err());
        }
    }
}
