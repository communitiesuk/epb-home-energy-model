use crate::input::{Schedule, WaterHeatingEvent};
use serde_json::Value;

pub type BooleanSchedule = Vec<bool>;
pub type NumericSchedule = Vec<Option<f64>>;

const MAIN_SCHEDULE: &str = "main";

pub fn expand_boolean_schedule(schedule: Schedule, nullable: bool) -> BooleanSchedule {
    process_boolean_schedule_entries(schedule.get(MAIN_SCHEDULE).unwrap(), &schedule, nullable)
}

pub fn expand_numeric_schedule(schedule: Schedule, nullable: bool) -> NumericSchedule {
    process_numeric_schedule_entries(schedule.get(MAIN_SCHEDULE).unwrap(), &schedule, nullable)
}

fn process_boolean_schedule_entries(
    entries: &Value,
    schedules: &Schedule,
    nullable: bool,
) -> Vec<bool> {
    let mut expansion = vec![];
    match entries {
        Value::Array(entry_vec) => {
            for entry in entry_vec {
                let extension = process_boolean_schedule_entry(entry, schedules, nullable);
                expansion.extend(extension);
            }
        }
        _ => panic!("An individual schedule was, illegally, not provided as a list."),
    }

    expansion
}

fn process_boolean_schedule_entry(
    entry: &Value,
    subschedules: &Schedule,
    nullable: bool,
) -> Vec<bool> {
    match entry {
        Value::String(subschedule) => process_boolean_schedule_entries(
            subschedules.get(subschedule).unwrap(),
            subschedules,
            nullable,
        ),
        Value::Object(repeated_map) => std::iter::repeat(
            process_boolean_schedule_entry(
                repeated_map.get("value").unwrap(),
                subschedules,
                nullable,
            )
            .into_iter(),
        )
        .take(repeated_map.get("repeat").unwrap().as_u64().unwrap() as usize)
        .flatten()
        .collect(),
        Value::Bool(whether) => vec![*whether],
        Value::Null if nullable => vec![Default::default()],
        _ => panic!("Unexpected value in schedule that was sent."),
    }
}

fn process_numeric_schedule_entries(
    entries: &Value,
    schedules: &Schedule,
    nullable: bool,
) -> NumericSchedule {
    let mut expansion = vec![];
    match entries {
        Value::Array(entry_vec) => {
            for entry in entry_vec {
                let extension = process_numeric_schedule_entry(entry, schedules, nullable);
                expansion.extend(extension);
            }
        }
        _ => panic!("An individual schedule was, illegally, not provided as a list."),
    }

    expansion
}

fn process_numeric_schedule_entry(
    entry: &Value,
    subschedules: &Schedule,
    nullable: bool,
) -> Vec<Option<f64>> {
    match entry {
        Value::String(subschedule) => process_numeric_schedule_entries(
            subschedules.get(subschedule).unwrap(),
            subschedules,
            nullable,
        ),
        Value::Object(repeated_map) => std::iter::repeat(
            process_numeric_schedule_entry(
                repeated_map.get("value").unwrap(),
                subschedules,
                nullable,
            )
            .into_iter(),
        )
        .take(repeated_map.get("repeat").unwrap().as_u64().unwrap() as usize)
        .flatten()
        .collect(),
        Value::Number(number) => vec![number.as_f64()],
        Value::Null if nullable => vec![Default::default()],
        _ => panic!("Unexpected value in schedule that was sent."),
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ScheduleEvent {
    pub start: f64,
    pub duration: Option<f64>,
    pub temperature: Option<f64>,
}

impl From<Value> for ScheduleEvent {
    fn from(value: Value) -> Self {
        match value {
            Value::Object(event_map) => ScheduleEvent {
                start: event_map.get("start").unwrap().as_f64().unwrap(),
                duration: event_map.get("duration").map(|d| d.as_f64().unwrap()),
                temperature: event_map.get("temperature").map(|t| t.as_f64().unwrap()),
            },
            _ => panic!("Expected a JSON object when transforming into a schedule event"),
        }
    }
}

pub fn expand_events(
    events: Vec<Value>,
    simulation_timestep: f64,
    total_timesteps: usize,
) -> Vec<Option<Vec<ScheduleEvent>>> {
    let mut schedule: Vec<Option<Vec<ScheduleEvent>>> = vec![None; total_timesteps];
    for event in events {
        let starting_timestep = (event
            .as_object()
            .unwrap()
            .get("start")
            .unwrap()
            .as_f64()
            .unwrap()
            / simulation_timestep)
            .floor() as usize;
        match schedule.get_mut(starting_timestep).unwrap() {
            Some(events) => events.push(event.into()),
            None => {
                schedule[starting_timestep] = Some(vec![event.into()]);
            }
        }
    }

    schedule
}

impl From<&WaterHeatingEvent> for ScheduleEvent {
    fn from(event: &WaterHeatingEvent) -> Self {
        Self {
            start: event.start,
            duration: event.duration,
            temperature: Some(event.temperature),
        }
    }
}

pub fn expand_water_heating_events(
    events: Vec<&WaterHeatingEvent>,
    simulation_timestep: f64,
    total_timesteps: usize,
) -> Vec<Option<Vec<ScheduleEvent>>> {
    let mut schedule: Vec<Option<Vec<ScheduleEvent>>> = vec![None; total_timesteps];
    for event in events {
        let starting_timestep = (event.start / simulation_timestep).floor() as usize;
        match schedule.get_mut(starting_timestep).unwrap() {
            Some(events) => events.push(event.into()),
            None => schedule[starting_timestep] = Some(vec![event.into()]),
        }
    }

    schedule
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::Schedule;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use serde_json::json;
    use std::collections::HashMap;

    #[fixture]
    pub fn boolean_schedule() -> Schedule {
        HashMap::from([
            (
                "main".to_string(),
                json!([
                    {"value": "weekday", "repeat": 5},
                    "weekend", "weekend",
                ]),
            ),
            (
                "weekday".to_string(),
                json!([
                    {"value": false, "repeat": 7},
                    {"value": true, "repeat": 2},
                    {"value": false, "repeat": 7},
                    {"value": true, "repeat": 7},
                    false,
                ]),
            ),
            (
                "weekend".to_string(),
                json!([
                    {"value": false, "repeat": 7},
                    {"value": true, "repeat": 16},
                    false,
                ]),
            ),
        ])
    }

    #[fixture]
    pub fn boolean_schedule_expanded() -> BooleanSchedule {
        vec![
            // Weekday schedule (Mon)
            false, false, false, false, false, false, false, true, true, false, false, false, false,
            false, false, false, true, true, true, true, true, true, true, false,
            // Weekday schedule (Tue)
            false, false, false, false, false, false, false, true, true, false, false, false, false,
            false, false, false, true, true, true, true, true, true, true, false,
            // Weekday schedule (Wed)
            false, false, false, false, false, false, false, true, true, false, false, false, false,
            false, false, false, true, true, true, true, true, true, true, false,
            // Weekday schedule (Thu)
            false, false, false, false, false, false, false, true, true, false, false, false, false,
            false, false, false, true, true, true, true, true, true, true, false,
            // Weekday schedule (Fri)
            false, false, false, false, false, false, false, true, true, false, false, false, false,
            false, false, false, true, true, true, true, true, true, true, false,
            // Weekend schedule (Sat)
            false, false, false, false, false, false, false, true, true, true, true, true, true,
            true, true, true, true, true, true, true, true, true, true, false,
            // Weekend schedule (Sun)
            false, false, false, false, false, false, false, true, true, true, true, true, true,
            true, true, true, true, true, true, true, true, true, true, false,
        ]
    }

    #[rstest]
    pub fn should_expand_boolean_schedule_correctly(
        boolean_schedule: Schedule,
        boolean_schedule_expanded: BooleanSchedule,
    ) {
        assert_eq!(
            expand_boolean_schedule(boolean_schedule, false),
            boolean_schedule_expanded,
            "Incorrect expansion of Boolean schedule"
        );
    }

    #[fixture]
    pub fn numeric_schedule() -> Schedule {
        HashMap::from([(
            "main".to_string(),
            json!([300.0, 120.0, 220.0, 750.0, 890.0, 150.0, 550.0, 280.0,]),
        )])
    }

    #[fixture]
    pub fn numeric_schedule_expanded() -> NumericSchedule {
        vec![
            Some(300.0),
            Some(120.0),
            Some(220.0),
            Some(750.0),
            Some(890.0),
            Some(150.0),
            Some(550.0),
            Some(280.0),
        ]
    }

    #[rstest]
    pub fn should_expand_numeric_schedule_correctly(
        numeric_schedule: Schedule,
        numeric_schedule_expanded: NumericSchedule,
    ) {
        assert_eq!(
            expand_numeric_schedule(numeric_schedule, false),
            numeric_schedule_expanded,
            "Incorrect expansion of numeric schedule"
        );
    }

    #[fixture]
    pub fn gappy_numeric_schedule() -> Schedule {
        HashMap::from([(
            "main".to_string(),
            json!([300.0, 120.0, null, 750.0, 890.0, null, 550.0, 280.0,]),
        )])
    }

    #[fixture]
    pub fn gappy_numeric_schedule_expanded() -> NumericSchedule {
        vec![
            Some(300.0),
            Some(120.0),
            None,
            Some(750.0),
            Some(890.0),
            None,
            Some(550.0),
            Some(280.0),
        ]
    }

    #[rstest]
    pub fn should_expand_gappy_numeric_schedule_correctly(
        gappy_numeric_schedule: Schedule,
        gappy_numeric_schedule_expanded: NumericSchedule,
    ) {
        assert_eq!(
            expand_numeric_schedule(gappy_numeric_schedule, true),
            gappy_numeric_schedule_expanded,
            "Incorrect expansion of numeric schedule"
        );
    }

    #[fixture]
    pub fn events() -> Vec<Value> {
        json!([
            {"start": 2, "duration": 6},
            {"start": 2.1, "duration": 6},
            {"start": 3, "duration": 6},
        ])
        .as_array()
        .unwrap()
        .clone()
    }

    #[fixture]
    pub fn simulation_timestep() -> f64 {
        0.5
    }

    #[fixture]
    pub fn total_timesteps() -> usize {
        10
    }

    #[fixture]
    pub fn events_schedule() -> Vec<Option<Vec<ScheduleEvent>>> {
        vec![
            None,
            None,
            None,
            None,
            Some(vec![
                ScheduleEvent {
                    start: 2.0,
                    duration: Some(6.0),
                    temperature: None,
                },
                ScheduleEvent {
                    start: 2.1,
                    duration: Some(6.0),
                    temperature: None,
                },
            ]),
            None,
            Some(vec![ScheduleEvent {
                start: 3.0,
                duration: Some(6.0),
                temperature: None,
            }]),
            None,
            None,
            None,
        ]
    }

    #[rstest]
    pub fn should_expand_events_correctly(
        events: Vec<Value>,
        simulation_timestep: f64,
        total_timesteps: usize,
        events_schedule: Vec<Option<Vec<ScheduleEvent>>>,
    ) {
        assert_eq!(
            expand_events(events, simulation_timestep, total_timesteps),
            events_schedule,
            "incorrect expansion of event list to schedule"
        );
    }

    #[fixture]
    pub fn water_heating_events() -> Vec<WaterHeatingEvent> {
        vec![
            WaterHeatingEvent {
                start: 2.0,
                duration: Some(6.0),
                temperature: 41.0,
            },
            WaterHeatingEvent {
                start: 2.1,
                duration: Some(6.0),
                temperature: 42.0,
            },
            WaterHeatingEvent {
                start: 3.0,
                duration: Some(6.0),
                temperature: 43.0,
            },
        ]
    }

    #[fixture]
    pub fn water_events_schedule() -> Vec<Option<Vec<ScheduleEvent>>> {
        vec![
            None,
            None,
            None,
            None,
            Some(vec![
                ScheduleEvent {
                    start: 2.0,
                    duration: Some(6.0),
                    temperature: Some(41.0),
                },
                ScheduleEvent {
                    start: 2.1,
                    duration: Some(6.0),
                    temperature: Some(42.0),
                },
            ]),
            None,
            Some(vec![ScheduleEvent {
                start: 3.0,
                duration: Some(6.0),
                temperature: Some(43.0),
            }]),
            None,
            None,
            None,
        ]
    }

    #[rstest]
    pub fn should_expand_water_events_correctly(
        water_heating_events: Vec<WaterHeatingEvent>,
        simulation_timestep: f64,
        total_timesteps: usize,
        water_events_schedule: Vec<Option<Vec<ScheduleEvent>>>,
    ) {
        assert_eq!(
            expand_water_heating_events(
                water_heating_events.iter().collect(),
                simulation_timestep,
                total_timesteps
            ),
            water_events_schedule,
            "incorrect expansion of event list to schedule"
        );
    }
}
