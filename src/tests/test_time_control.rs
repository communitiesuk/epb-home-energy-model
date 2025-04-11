mod test_time_control {
    use std::sync::Arc;

    use crate::core::controls::time_control::*;
    use crate::external_conditions::{DaylightSavingsConfig, ExternalConditions};
    use crate::input::{ControlCombinations, ControlLogicType, ExternalSensor};
    use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
    use indexmap::IndexMap;
    use itertools::Itertools;
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
