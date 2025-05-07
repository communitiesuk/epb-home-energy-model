use crate::input::{WaterHeatingEvent, WaterHeatingEventType};
use anyhow::{anyhow, bail};
#[cfg(test)]
use serde_json::Value;

pub(crate) fn reject_nulls<T>(vec_of_options: Vec<Option<T>>) -> anyhow::Result<Vec<T>> {
    vec_of_options
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .ok_or_else(|| anyhow!("A null was in a schedule when it was not expected."))
}

pub(crate) fn reject_nones<T>(vec_of_options: Vec<Option<T>>) -> anyhow::Result<Vec<Option<T>>> {
    if vec_of_options.iter().any(|o| o.is_none()) {
        bail!("A null value was in a schedule when it was not expected.")
    }
    Ok(vec_of_options)
}

pub(crate) fn expand_boolean_schedule(schedule: &BooleanSchedule) -> Vec<Option<bool>> {
    schedule.expand()
}

pub(crate) fn expand_numeric_schedule(schedule: &NumericSchedule) -> Vec<Option<f64>> {
    schedule.expand()
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ScheduleEvent {
    pub(crate) start: f64,
    pub(crate) duration: Option<f64>,
    pub(crate) volume: Option<f64>,
    pub(crate) temperature: Option<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum WaterScheduleEventType {
    Shower,
    Bath,
    Other,
}

impl From<WaterHeatingEventType> for WaterScheduleEventType {
    fn from(value: WaterHeatingEventType) -> Self {
        match value {
            WaterHeatingEventType::Shower => Self::Shower,
            WaterHeatingEventType::Bath => Self::Bath,
            WaterHeatingEventType::Other => Self::Other,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct TypedScheduleEvent {
    pub(crate) start: f64,
    pub(crate) duration: Option<f64>,
    pub(crate) temperature: f64,
    pub(crate) name: String,
    pub(crate) event_type: WaterScheduleEventType,
    pub(crate) volume: Option<f64>,
    pub(crate) warm_volume: Option<f64>,
    pub(crate) pipework_volume: Option<f64>,
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
            volume,
            temperature,
        } = event;
        Ok(Self {
            start,
            duration,
            temperature: temperature
                .ok_or_else(|| anyhow!("Temperature was expected to be set on a water event."))?,
            name,
            event_type,
            volume,
            warm_volume: None,
            pipework_volume: None,
        })
    }
}

#[cfg(test)]
impl TryFrom<&Value> for ScheduleEvent {
    type Error = anyhow::Error;

    fn try_from(value: &Value) -> anyhow::Result<Self> {
        use anyhow::bail;

        match value {
            Value::Object(event_map) => Ok(ScheduleEvent {
                start: event_map
                    .get("start")
                    .ok_or_else(|| anyhow!("Start key was not found in an event."))?
                    .as_f64()
                    .ok_or_else(|| anyhow!("Start value in event was expected to be numeric"))?,
                duration: event_map.get("duration").map(|d| d.as_f64().unwrap()),
                volume: event_map.get("volume").map(|v| v.as_f64().unwrap()),
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
#[cfg(test)]
pub(crate) fn expand_events_from_json_values(
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
pub(crate) fn expand_events(
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
            volume: event.volume,
        }
    }
}

/// Data structures representing how schedules can be provided as input (in JSON).
pub(crate) mod input {
    use itertools::Itertools;
    use serde::{Deserialize, Serialize};
    use std::collections::HashMap;

    #[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
    #[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
    pub(crate) struct Schedule<T: Copy> {
        pub(crate) main: Vec<ScheduleEntry<T>>,
        #[serde(flatten)]
        pub(crate) references: HashMap<String, ScheduleReferenceEntry<T>>,
    }

    impl<T> Schedule<T>
    where
        T: Copy,
    {
        pub(super) fn expand(&self) -> Vec<Option<T>> {
            self.main
                .iter()
                .flat_map(|entry| self.expand_entry(entry))
                .collect()
        }

        fn expand_entry(&self, entry: &ScheduleEntry<T>) -> Vec<Option<T>> {
            match entry {
                ScheduleEntry::Null(_) => vec![None],
                ScheduleEntry::Value(v) => vec![Some(*v)],
                ScheduleEntry::Repeater(repeater) => std::iter::repeat([match &repeater.value {
                    ScheduleRepeaterValue::Reference(reference) => self.expand_reference(reference),
                    ScheduleRepeaterValue::Entry(ScheduleRepeaterEntry::Null(_)) => {
                        vec![None]
                    }
                    ScheduleRepeaterValue::Entry(ScheduleRepeaterEntry::Value(v)) => {
                        vec![Some(*v)]
                    }
                }])
                .take(repeater.repeat)
                .flatten()
                .flatten()
                .collect_vec(),
                ScheduleEntry::Reference(reference) => self.expand_reference(reference),
            }
        }

        fn expand_reference(&self, reference: &str) -> Vec<Option<T>> {
            match &self.references[reference] {
                ScheduleReferenceEntry::Single(entry) => self.expand_entry(entry),
                ScheduleReferenceEntry::Multi(entries) => entries
                    .iter()
                    .flat_map(|entry| self.expand_entry(entry))
                    .collect_vec(),
            }
        }
    }

    #[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
    #[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
    #[serde(untagged)]
    pub(crate) enum ScheduleEntry<T: Copy> {
        Null(()),
        Value(T),
        Repeater(ScheduleRepeater<T>),
        Reference(String),
    }

    #[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
    #[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
    #[serde(untagged)]
    pub(crate) enum ScheduleReferenceEntry<T: Copy> {
        Multi(Vec<ScheduleEntry<T>>),
        Single(ScheduleEntry<T>),
    }

    #[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
    #[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
    #[serde(untagged)]
    pub(crate) enum ScheduleRepeaterEntry<T> {
        Null(()),
        Value(T),
    }

    #[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
    #[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
    pub(crate) struct ScheduleRepeater<T: Copy> {
        pub(crate) value: ScheduleRepeaterValue<T>,
        pub(crate) repeat: usize,
    }

    #[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
    #[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
    #[serde(untagged)]
    pub(crate) enum ScheduleRepeaterValue<T: Copy> {
        Reference(String),
        Entry(ScheduleRepeaterEntry<T>),
    }

    pub(crate) type BooleanSchedule = Schedule<bool>;
    pub(crate) type NumericSchedule = Schedule<f64>;

    #[cfg(test)]
    mod schedule_test {
        use super::*;
        use rstest::*;
        use serde_json::json;

        #[rstest]
        fn make_boolean_schedule_with_simple_values() {
            let _boolean_schedule: BooleanSchedule = serde_json::from_value(json!({
                "main": [true, false]
            }))
            .unwrap();
        }

        #[rstest]
        fn make_boolean_schedule_with_repeater() {
            let _boolean_schedule: BooleanSchedule = serde_json::from_value(json!({
                "main": [true, false, {"value": true, "repeat": 5}]
            }))
            .unwrap();
        }

        #[rstest]
        fn boolean_schedule_without_main_fails() {
            assert!(serde_json::from_value::<BooleanSchedule>(json!({
                "weekday": [true, false, {"value": true, "repeat": 5}]
            }))
            .is_err());
        }

        #[rstest]
        fn boolean_schedule_allows_nulls() {
            assert!(serde_json::from_value::<BooleanSchedule>(json!({
                "main": [true, null, {"value": true, "repeat": 5}]
            }))
            .is_ok());
        }

        #[rstest]
        fn numeric_schedule_allows_numbers_and_repeats() {
            assert!(serde_json::from_value::<NumericSchedule>(json!({
                "main": [21.0, null, {"value": 22.0, "repeat": 5}]
            }))
            .is_ok());
        }

        #[rstest]
        fn numeric_schedule_with_references_as_entries_allowed() {
            assert!(serde_json::from_value::<BooleanSchedule>(json!({
                "main": [
                    {"value": "weekday", "repeat": 5},
                    "weekend", "weekend",
                ],
                "weekday": [
                    {"value": false, "repeat": 7},
                    {"value": true, "repeat": 2},
                    {"value": false, "repeat": 7},
                    {"value": true, "repeat": 7},
                    false,
                ],
                "weekend": [
                    {"value": false, "repeat": 7},
                    {"value": true, "repeat": 16},
                    false,
                ]
            }))
            .is_ok());
        }

        #[rstest]
        fn schedule_with_repeaters_with_references_allowed() {
            assert!(serde_json::from_value::<NumericSchedule>(json!({
                "main": [{"repeat": 53, "value": "week"}],
                "week": [{"repeat": 5, "value": "weekday"}, {"repeat": 2, "value": "weekend"}],
                "weekday": [22.0],
                "weekend": [23.0],
            }))
            .is_ok());
        }

        #[rstest]
        fn schedule_with_repeaters_with_references_to_single_values_allowed() {
            assert!(serde_json::from_value::<NumericSchedule>(json!({
                "main": [{"repeat": 53, "value": "week"}],
                "week": [{"repeat": 5, "value": "weekday"}, {"repeat": 2, "value": "weekend"}],
                "weekday": 22.0,
                "weekend": 23.0,
            }))
            .is_ok());
        }
    }
}

pub(crate) use input::BooleanSchedule;
pub(crate) use input::NumericSchedule;

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use serde_json::json;

    #[fixture]
    pub fn boolean_schedule() -> BooleanSchedule {
        serde_json::from_value(json!({
            "main": [
                {"value": "weekday", "repeat": 5},
                "weekend", "weekend",
            ],
            "weekday": [
                {"value": false, "repeat": 7},
                {"value": true, "repeat": 2},
                {"value": false, "repeat": 7},
                {"value": true, "repeat": 7},
                false,
            ],
            "weekend": [
                {"value": false, "repeat": 7},
                {"value": true, "repeat": 16},
                false,
            ]
        }))
        .unwrap()
    }

    #[fixture]
    pub fn boolean_schedule_expanded() -> Vec<bool> {
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
        boolean_schedule: BooleanSchedule,
        boolean_schedule_expanded: Vec<bool>,
    ) {
        assert_eq!(
            reject_nulls(expand_boolean_schedule(&boolean_schedule)).unwrap(),
            boolean_schedule_expanded,
            "Incorrect expansion of Boolean schedule"
        );
    }

    #[fixture]
    pub fn numeric_schedule() -> NumericSchedule {
        serde_json::from_value(json!({
            "main": [300.0, 120.0, 220.0, 750.0, 890.0, 150.0, 550.0, 280.0]
        }))
        .unwrap()
    }

    #[fixture]
    pub fn numeric_schedule_expanded() -> Vec<Option<f64>> {
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
        numeric_schedule: NumericSchedule,
        numeric_schedule_expanded: Vec<Option<f64>>,
    ) {
        assert_eq!(
            expand_numeric_schedule(&numeric_schedule),
            numeric_schedule_expanded,
            "Incorrect expansion of numeric schedule"
        );
    }

    #[fixture]
    pub fn gappy_numeric_schedule() -> NumericSchedule {
        serde_json::from_value(json!({
            "main": [300.0, 120.0, null, 750.0, 890.0, null, 550.0, 280.0]
        }))
        .unwrap()
    }

    #[fixture]
    pub fn gappy_numeric_schedule_expanded() -> Vec<Option<f64>> {
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
        gappy_numeric_schedule: NumericSchedule,
        gappy_numeric_schedule_expanded: Vec<Option<f64>>,
    ) {
        assert_eq!(
            expand_numeric_schedule(&gappy_numeric_schedule),
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
                    volume: None,
                    warm_volume: None,
                    pipework_volume: None,
                },
                TypedScheduleEvent {
                    start: 2.1,
                    duration: Some(6.0),
                    temperature: 52.0,
                    name: "name".to_string(),
                    event_type: WaterScheduleEventType::Shower,
                    volume: None,
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
                volume: None,
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
