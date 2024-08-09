use crate::input::{Schedule, WaterHeatingEvent};
use anyhow::{anyhow, bail};
use serde_json::Value;

pub type BooleanSchedule = Vec<bool>;
pub type NumericSchedule = Vec<Option<f64>>;

const MAIN_SCHEDULE: &str = "main";

pub fn expand_boolean_schedule(schedule: &Schedule, nullable: bool) -> BooleanSchedule {
    process_boolean_schedule_entries(schedule.get(MAIN_SCHEDULE).unwrap(), schedule, nullable)
}

pub fn expand_numeric_schedule(schedule: &Schedule, nullable: bool) -> NumericSchedule {
    process_numeric_schedule_entries(schedule.get(MAIN_SCHEDULE).unwrap(), schedule, nullable)
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

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WaterScheduleEventType {
    Shower,
    Bath,
    Other,
}

#[derive(Clone, Debug, PartialEq)]
pub struct TypedScheduleEvent {
    pub start: f64,
    pub duration: Option<f64>,
    pub temperature: f64,
    pub name: String,
    pub event_type: WaterScheduleEventType,
    pub warm_volume: Option<f64>,
    pub pipework_volume: Option<f64>,
}

impl TypedScheduleEvent {
    pub fn from_simple_event(
        event: ScheduleEvent,
        name: String,
        event_type: WaterScheduleEventType,
    ) -> anyhow::Result<Self> {
        let ScheduleEvent {
            start,
            duration,
            temperature,
        } = event;
        Ok(Self {
            start,
            duration,
            temperature: temperature
                .ok_or_else(|| anyhow!("Temperature was expected to be set on a water event."))?,
            name,
            event_type,
            warm_volume: None,
            pipework_volume: None,
        })
    }
}

impl TryFrom<&Value> for ScheduleEvent {
    type Error = anyhow::Error;

    fn try_from(value: &Value) -> anyhow::Result<Self> {
        match value {
            Value::Object(event_map) => Ok(ScheduleEvent {
                start: event_map
                    .get("start")
                    .ok_or_else(|| anyhow!("Start key was not found in an event."))?
                    .as_f64()
                    .ok_or_else(|| anyhow!("Start value in event was expected to be numeric"))?,
                duration: event_map.get("duration").map(|d| d.as_f64().unwrap()),
                temperature: event_map.get("temperature").map(|t| t.as_f64().unwrap()),
            }),
            _ => bail!("Expected a JSON object when transforming into a schedule event"),
        }
    }
}

/// Construct or update a schedule from a list of events, appending the event type to each event and
/// ensuring events are ordered by 'start' time within each timestep.
///
/// Arguments:
/// * `events` - list of event dictionaries (directly parsed from JSON), where the 'start' element gives
///              the start time of the event, in hours from the start of the simulation
/// * `sim_timestep` - length of simulation timestep, in hours
/// * `tot_timesteps` - total number of timesteps in the simulation
/// * `name`
/// * `event_type` - type of the events being processed (e.g., "Shower", "Bath", "Others")
/// * `schedule` - the existing schedule dictionary to update
pub fn expand_events_from_json_values(
    events: Vec<Value>,
    simulation_timestep: f64,
    total_timesteps: usize,
    name: &str,
    event_type: WaterScheduleEventType,
    schedule: Vec<Option<Vec<TypedScheduleEvent>>>,
) -> anyhow::Result<Vec<Option<Vec<TypedScheduleEvent>>>> {
    expand_events(
        events
            .iter()
            .map(|json| -> anyhow::Result<ScheduleEvent> { ScheduleEvent::try_from(json) })
            .collect::<anyhow::Result<Vec<_>>>()?,
        simulation_timestep,
        total_timesteps,
        name,
        event_type,
        schedule,
    )
}

/// Construct or update a schedule from a list of events, appending the event type to each event and
/// ensuring events are ordered by 'start' time within each timestep.
///
/// Arguments:
/// * `events` - list of schedule events, where the 'start' element gives
///              the start time of the event, in hours from the start of the simulation
/// * `sim_timestep` - length of simulation timestep, in hours
/// * `tot_timesteps` - total number of timesteps in the simulation
/// * `name`
/// * `event_type` - type of the events being processed (e.g., "Shower", "Bath", "Others")
/// * `schedule` - the existing schedule dictionary to update
pub fn expand_events(
    events: Vec<ScheduleEvent>,
    simulation_timestep: f64,
    total_timesteps: usize,
    name: &str,
    event_type: WaterScheduleEventType,
    mut schedule: Vec<Option<Vec<TypedScheduleEvent>>>,
) -> anyhow::Result<Vec<Option<Vec<TypedScheduleEvent>>>> {
    for event in events {
        let starting_timestep = (event.start / simulation_timestep).floor() as usize;

        if starting_timestep < total_timesteps {
            let event_with_type_name =
                TypedScheduleEvent::from_simple_event(event, name.to_string(), event_type)?;

            match schedule.get_mut(starting_timestep).unwrap() {
                Some(events) => {
                    // Insert the event into the correct position to maintain order by 'start' time
                    let mut inserted = false;
                    for (i, existing_event) in events.iter().enumerate() {
                        if existing_event.start > event_with_type_name.start {
                            events.insert(i, event_with_type_name.clone());
                            inserted = true;
                            break;
                        }
                    }
                    if !inserted {
                        events.push(event_with_type_name);
                    }
                }
                None => schedule[starting_timestep] = Some(vec![event_with_type_name]),
            }
        }
    }

    Ok(schedule)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::Schedule;
    use indexmap::IndexMap;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use serde_json::json;

    #[fixture]
    pub fn boolean_schedule() -> Schedule {
        IndexMap::from([
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
            expand_boolean_schedule(&boolean_schedule, false),
            boolean_schedule_expanded,
            "Incorrect expansion of Boolean schedule"
        );
    }

    #[fixture]
    pub fn numeric_schedule() -> Schedule {
        IndexMap::from([(
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
            expand_numeric_schedule(&numeric_schedule, false),
            numeric_schedule_expanded,
            "Incorrect expansion of numeric schedule"
        );
    }

    #[fixture]
    pub fn gappy_numeric_schedule() -> Schedule {
        IndexMap::from([(
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
            expand_numeric_schedule(&gappy_numeric_schedule, true),
            gappy_numeric_schedule_expanded,
            "Incorrect expansion of numeric schedule"
        );
    }

    #[fixture]
    pub fn events() -> Vec<Value> {
        json!([
            {"start": 2, "duration": 6, "temperature": 52},
            {"start": 2.1, "duration": 6, "temperature": 52},
            {"start": 3, "duration": 6, "temperature": 52},
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
    pub fn events_schedule() -> Vec<Option<Vec<TypedScheduleEvent>>> {
        vec![
            None,
            None,
            None,
            None,
            Some(vec![
                TypedScheduleEvent {
                    start: 2.0,
                    duration: Some(6.0),
                    temperature: 52.0,
                    name: "name".to_string(),
                    event_type: WaterScheduleEventType::Shower,
                    warm_volume: None,
                    pipework_volume: None,
                },
                TypedScheduleEvent {
                    start: 2.1,
                    duration: Some(6.0),
                    temperature: 52.0,
                    name: "name".to_string(),
                    event_type: WaterScheduleEventType::Shower,
                    warm_volume: None,
                    pipework_volume: None,
                },
            ]),
            None,
            Some(vec![TypedScheduleEvent {
                start: 3.0,
                duration: Some(6.0),
                temperature: 52.0,
                name: "name".to_string(),
                event_type: WaterScheduleEventType::Shower,
                warm_volume: None,
                pipework_volume: None,
            }]),
            None,
            None,
            None,
        ]
    }

    #[rstest]
    pub fn test_expand_events(
        events: Vec<Value>,
        simulation_timestep: f64,
        total_timesteps: usize,
        events_schedule: Vec<Option<Vec<TypedScheduleEvent>>>,
    ) {
        let schedule_to_test = vec![None; total_timesteps];
        assert_eq!(
            expand_events_from_json_values(
                events,
                simulation_timestep,
                total_timesteps,
                "name",
                WaterScheduleEventType::Shower,
                schedule_to_test
            )
            .unwrap(),
            events_schedule,
            "incorrect expansion of event list to schedule"
        );
    }
}
