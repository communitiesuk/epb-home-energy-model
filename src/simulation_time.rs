use serde::Deserialize;

pub const HOURS_IN_DAY: u32 = 24;

// # Define hours that start each month (and end next month). Note there are 13
// # values so that end of final month is handled correctly.
// # E.g. Jan is hours 0-743
const MONTH_START_END_HOURS: [u32; 13] = [
    0, 744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760,
];

#[derive(Clone, Debug, Deserialize)]
pub struct SimulationTime {
    #[serde(rename(deserialize = "start"))]
    start_time: f64,
    #[serde(rename(deserialize = "end"))]
    end_time: f64,
    step: f64,
}

impl SimulationTime {
    pub fn new(start_time: f64, end_time: f64, step: f64) -> Self {
        Self {
            start_time,
            end_time,
            step,
        }
    }

    fn total_steps(&self) -> i32 {
        ((self.end_time - self.start_time) / self.step).ceil() as i32
    }

    pub(crate) fn iter(&self) -> SimulationTimeIterator {
        SimulationTimeIterator::from((*self).clone())
    }
}

#[derive(Clone)]
pub struct SimulationTimeIterator {
    current_index: usize,
    current_time: f64,
    started: bool,
    simulation_time: SimulationTime,
}

impl SimulationTimeIterator {
    fn from(simulation_time: SimulationTime) -> Self {
        SimulationTimeIterator {
            current_index: 0,
            current_time: simulation_time.start_time,
            started: false,
            simulation_time,
        }
    }

    pub fn current_index(&self) -> usize {
        self.current_index
    }

    pub fn time_series_idx(&self, start_day: u32, step: f64) -> usize {
        ((self.current_time - (start_day * HOURS_IN_DAY) as f64) / step) as usize
    }

    pub fn time_series_idx_days(&self, start_day: u32, step: f64) -> u32 {
        let current_day = self.current_time as u32 / HOURS_IN_DAY;
        // # TODO: (Potential) Decide from which hour of the day the system should be targeting next day charge level
        // # Currently 9pm
        if self.current_time.floor() >= 21.0 {
            ((current_day + 1 - start_day) as f64 / step) as u32
        } else {
            ((current_day - start_day) as f64 / step) as u32
        }
    }

    pub fn current_hour(&self) -> u32 {
        self.current_time.floor() as u32
    }

    pub fn current_day(&self) -> u32 {
        self.current_time as u32 / HOURS_IN_DAY
    }

    pub fn current_month(&self) -> Option<u32> {
        let current_hour = self.current_hour();
        for (i, end_hour) in MONTH_START_END_HOURS.iter().enumerate() {
            if current_hour < *end_hour {
                return Some((i - 1) as u32);
            }
        }
        None
    }

    pub fn current_month_start_end_hours(&self) -> (u32, u32) {
        let month_idx = self.current_month().unwrap() as usize;
        (
            MONTH_START_END_HOURS[month_idx],
            MONTH_START_END_HOURS[month_idx + 1],
        )
    }
}

#[derive(Debug)]
pub struct SimulationTimeIteration {
    pub index: usize,
    pub time: f64,
    pub timestep: f64,
}

impl SimulationTimeIteration {
    pub fn current_hour(&self) -> u32 {
        self.time.floor() as u32
    }

    pub fn hour_of_day(&self) -> u32 {
        self.current_hour() % HOURS_IN_DAY
    }

    pub fn current_day(&self) -> u32 {
        self.time as u32 / HOURS_IN_DAY
    }

    pub fn current_month(&self) -> Option<u32> {
        let current_hour = self.current_hour();
        for (i, end_hour) in MONTH_START_END_HOURS.iter().enumerate() {
            if current_hour < *end_hour {
                return Some((i - 1) as u32);
            }
        }
        None
    }

    pub fn current_month_start_end_hours(&self) -> (u32, u32) {
        let month_idx = self.current_month().unwrap() as usize;
        (
            MONTH_START_END_HOURS[month_idx],
            MONTH_START_END_HOURS[month_idx + 1],
        )
    }

    pub fn time_series_idx(&self, start_day: u32, step: f64) -> usize {
        ((self.time - (start_day * HOURS_IN_DAY) as f64) / step) as usize
    }
}

impl Iterator for SimulationTimeIterator {
    type Item = SimulationTimeIteration;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.started && self.simulation_time.start_time != self.simulation_time.end_time {
            self.started = true;
            return Some(SimulationTimeIteration {
                index: 0,
                time: self.simulation_time.start_time,
                timestep: self.simulation_time.step,
            });
        }
        match self.current_time < (self.simulation_time.end_time - self.simulation_time.step) {
            true => {
                self.current_index += 1;
                self.current_time += self.simulation_time.step;
                Some(SimulationTimeIteration {
                    index: self.current_index,
                    time: self.current_time,
                    timestep: self.simulation_time.step,
                })
            }
            false => None,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rstest::*;

    #[fixture]
    pub fn timestep() -> f64 {
        0.5
    }

    #[fixture]
    pub fn simtime() -> SimulationTime {
        SimulationTime {
            start_time: 742.0,
            end_time: 746.0,
            step: timestep(),
        }
    }

    #[rstest]
    fn should_have_correct_total_steps(simtime: SimulationTime) {
        assert_eq!(simtime.total_steps(), 8)
    }

    #[rstest]
    fn should_iterate_correctly(simtime: SimulationTime, timestep: f64) {
        let hours = [742, 742, 743, 743, 744, 744, 745, 745];
        let hours_of_day = [22, 22, 23, 23, 0, 0, 1, 1];
        let current_days = [30, 30, 30, 30, 31, 31, 31, 31];
        let current_months = [0, 0, 0, 0, 1, 1, 1, 1];
        let current_month_start_end_hours = [
            (0, 744),
            (0, 744),
            (0, 744),
            (0, 744),
            (744, 1416),
            (744, 1416),
            (744, 1416),
            (744, 1416),
        ];
        let mut simulation_time_iter = simtime.iter();
        let mut i = 0;
        while let Some(item) = simulation_time_iter.next() {
            assert_eq!(
                item.index, i,
                "current index is {0} with time {1}, but test iterator is {i}",
                item.index, item.time
            );
            assert_eq!(item.time, i as f64 * timestep + 742.0);
            assert_eq!(item.timestep, timestep);
            assert_eq!(item.current_hour(), hours[i]);
            assert_eq!(item.hour_of_day(), hours_of_day[i]);
            assert_eq!(item.current_day(), current_days[i]);
            assert_eq!(
                simulation_time_iter.time_series_idx(0, 1.0),
                hours[i] as usize
            );
            assert_eq!(item.current_month().unwrap(), current_months[i]);
            assert_eq!(
                item.current_month_start_end_hours(),
                current_month_start_end_hours[i]
            );
            i += 1;
        }
    }
}
