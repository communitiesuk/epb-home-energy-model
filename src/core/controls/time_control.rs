// This module provides structs to model time controls

use crate::simulation_time::{SimulationTimeIteration, HOURS_IN_DAY};

pub enum Control {
    OnOffTimeControl(OnOffTimeControl),
    ToUChargeControl(ToUChargeControl),
    OnOffMinimisingTimeControl(OnOffMinimisingTimeControl),
    SetpointTimeControl(SetpointTimeControl),
}

/// An object to model a time-only control with on/off (not modulating) operation
pub struct OnOffTimeControl {
    /// list of boolean values where true means "on" (one entry per hour)
    schedule: Vec<bool>,
    start_day: u32,
    time_series_step: f64,
}

impl OnOffTimeControl {
    pub fn new(schedule: Vec<bool>, start_day: u32, time_series_step: f64) -> OnOffTimeControl {
        OnOffTimeControl {
            schedule,
            start_day,
            time_series_step,
        }
    }
}

impl OnOffTimeControl {
    fn is_on(&self, timestep: &SimulationTimeIteration) -> bool {
        self.schedule[timestep.time_series_idx(self.start_day, self.time_series_step)]
    }
}

/// An object to model a control that governs electrical charging of a heat storage device
/// that can respond to signals from the grid, for example when carbon intensity is low
pub struct ToUChargeControl {
    /// list of boolean values where true means "on" (one entry per hour)
    pub schedule: Vec<bool>,
    pub start_day: u32,
    pub time_series_step: f64,
    /// Proportion of the charge targeted for each day
    pub charge_level: Vec<f64>,
}

impl ToUChargeControl {
    pub fn is_on(&self, iteration: &SimulationTimeIteration) -> bool {
        self.schedule[iteration.time_series_idx(self.start_day, self.time_series_step)]
    }

    pub fn target_charge(&self, timestep: &SimulationTimeIteration) -> f64 {
        self.charge_level[timestep.time_series_idx_days(self.start_day, self.time_series_step)]
    }
}

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
        duration_advanced_start: Option<f64>,
        timestep: f64,
    ) -> Result<Self, &'static str> {
        if setpoint_min.is_some() && setpoint_max.is_some() && default_to_max.is_none() {
            return Err(
                "default_to_max should be set when both setpoint_min and setpoint_max are set",
            );
        }
        let duration_advanced_start = duration_advanced_start.unwrap_or(0.0);

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

    pub fn in_required_period(&self, timestep: &SimulationTimeIteration) -> bool {
        let schedule_idx = timestep.time_series_idx(self.start_day, self.time_series_step);

        self.schedule[schedule_idx].is_some()
    }

    pub fn is_on(&self, timestep: &SimulationTimeIteration) -> bool {
        let schedule_idx = timestep.time_series_idx(self.start_day, self.time_series_step);

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

    pub fn setpnt(&self, timestep: &SimulationTimeIteration) -> Option<f64> {
        let schedule_idx = timestep.time_series_idx(self.start_day, self.time_series_step);

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

#[cfg(test)]
mod test {
    use super::*;
    use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
    use itertools::Itertools;
    use rstest::*;

    const ON_OFF_SCHEDULE: [bool; 8] = [true, false, true, true, false, true, false, false];

    #[fixture]
    pub fn simulation_time() -> SimulationTimeIterator {
        SimulationTime::new(0.0, 8.0, 1.0).iter()
    }

    #[fixture]
    pub fn on_off_time_control() -> OnOffTimeControl {
        OnOffTimeControl::new(ON_OFF_SCHEDULE.to_vec(), 0, 1.0)
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
            None,
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
            None,
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
            None,
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
            None,
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
            Some(1.0),
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
            Some(1.0),
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
                setpoint_time_control.in_required_period(&it),
                results[it.index],
                "incorrect in_required_period value returned for control with no min or max set, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_min.in_required_period(&it),
                results[it.index],
                "incorrect in_required_period value returned for control with min set, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_max.in_required_period(&it),
                results[it.index],
                "incorrect in_required_period value returned for control with max set, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_minmax.in_required_period(&it),
                results[it.index],
                "incorrect in_required_period value returned for control with min and max set, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_advstart.in_required_period(&it),
                results[it.index],
                "incorrect in_required_period value returned for control with advanced start, iteration {}",
                it.index + 1
            );
            assert_eq!(
                setpoint_time_control_advstart_minmax.in_required_period(&it),
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
}
