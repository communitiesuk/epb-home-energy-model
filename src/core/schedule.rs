use crate::input::Schedule;
use serde_json::Value;
use std::collections::HashMap;
use std::ops::Index;

pub type BooleanSchedule = Vec<bool>;
pub type NumericSchedule = Vec<f64>;

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
                let extension = process_boolean_schedule_entry(entry, &schedules, nullable);
                expansion.extend(extension);
            }
        }
        _ => panic!("An individual schedule was, illegally, not provided as a list."),
    }

    println!("expansion: {:?}", expansion);

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
                let extension = process_numeric_schedule_entry(entry, &schedules, nullable);
                expansion.extend(extension);
            }
        }
        _ => panic!("An individual schedule was, illegally, not provided as a list."),
    }

    println!("expansion: {:?}", expansion);

    expansion
}

fn process_numeric_schedule_entry(
    entry: &Value,
    subschedules: &Schedule,
    nullable: bool,
) -> Vec<f64> {
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
        Value::Number(number) => vec![number.as_f64().unwrap()],
        _ => panic!("Unexpected value in schedule that was sent."),
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ScheduleEvent {
    start: f64,
    duration: f64,
}

impl From<Value> for ScheduleEvent {
    fn from(value: Value) -> Self {
        match value {
            Value::Object(event_map) => ScheduleEvent {
                start: event_map.get("start").unwrap().as_f64().unwrap(),
                duration: event_map.get("duration").unwrap().as_f64().unwrap(),
            },
            _ => panic!("Expected a JSON object when transforming into a schedule event"),
        }
    }
}

pub fn expand_events(
    events: Value,
    simulation_timestep: f64,
    total_timesteps: usize,
) -> Vec<Option<Vec<ScheduleEvent>>> {
    let mut schedule: Vec<Option<Vec<ScheduleEvent>>> = vec![None; total_timesteps];
    match events {
        Value::Array(events) => {
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
        }
        _ => panic!("Events list was expected to be a list in the input."),
    }

    schedule
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::input::Schedule;
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
        vec![300.0, 120.0, 220.0, 750.0, 890.0, 150.0, 550.0, 280.0]
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
    pub fn events() -> Value {
        json!([
            {"start": 2, "duration": 6},
            {"start": 2.1, "duration": 6},
            {"start": 3, "duration": 6},
        ])
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
                    duration: 6.0,
                },
                ScheduleEvent {
                    start: 2.1,
                    duration: 6.0,
                },
            ]),
            None,
            Some(vec![ScheduleEvent {
                start: 3.0,
                duration: 6.0,
            }]),
            None,
            None,
            None,
        ]
    }

    #[rstest]
    pub fn should_expand_events_correctly(
        events: Value,
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
}
