mod elec_battery_tests {
    use std::sync::Arc;

    use crate::core::energy_supply::elec_battery::*;
    use crate::external_conditions::{
        DaylightSavingsConfig, ExternalConditions, ShadingObject, ShadingObjectType, ShadingSegment,
    };
    use crate::input::BatteryLocation;
    use crate::simulation_time::SimulationTime;
    use approx::assert_relative_eq;
    use rstest::*;

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 8., 1.)
    }

    #[fixture]
    pub fn external_conditions(simulation_time: SimulationTime) -> ExternalConditions {
        ExternalConditions::new(
            &simulation_time.iter(),
            vec![0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0],
            vec![3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1],
            vec![0., 20., 40., 60., 0., 20., 40., 60.],
            vec![11., 25., 42., 52., 60., 44., 28., 15.],
            vec![11., 25., 42., 52., 60., 44., 28., 15.],
            vec![0.2; 8],
            51.42,
            -0.75,
            0,
            0,
            Some(0),
            1.,
            Some(1),
            Some(DaylightSavingsConfig::NotApplicable),
            false,
            false,
            vec![
                ShadingSegment {
                    number: 1,
                    start: 180.,
                    end: 135.,
                    shading_objects: None,
                    ..Default::default()
                },
                ShadingSegment {
                    number: 2,
                    start: 135.,
                    end: 90.,
                    shading_objects: Some(vec![ShadingObject {
                        object_type: ShadingObjectType::Overhang,
                        height: 2.2,
                        distance: 6.,
                    }]),
                    ..Default::default()
                },
                ShadingSegment {
                    number: 3,
                    start: 90.,
                    end: 45.,
                    shading_objects: None,
                    ..Default::default()
                },
                ShadingSegment {
                    number: 4,
                    start: 45.,
                    end: 0.,
                    shading_objects: Some(vec![
                        ShadingObject {
                            object_type: ShadingObjectType::Obstacle,
                            height: 40.,
                            distance: 4.,
                        },
                        ShadingObject {
                            object_type: ShadingObjectType::Overhang,
                            height: 3.,
                            distance: 7.,
                        },
                    ]),
                    ..Default::default()
                },
                ShadingSegment {
                    number: 5,
                    start: 0.,
                    end: -45.,
                    shading_objects: Some(vec![ShadingObject {
                        object_type: ShadingObjectType::Obstacle,
                        height: 3.,
                        distance: 8.,
                    }]),
                    ..Default::default()
                },
                ShadingSegment {
                    number: 6,
                    start: -45.,
                    end: -90.,
                    shading_objects: None,
                    ..Default::default()
                },
                ShadingSegment {
                    number: 7,
                    start: -90.,
                    end: -135.,
                    shading_objects: None,
                    ..Default::default()
                },
                ShadingSegment {
                    number: 8,
                    start: -135.,
                    end: -180.,
                    shading_objects: None,
                    ..Default::default()
                },
            ]
            .into(),
        )
    }

    #[fixture]
    pub fn electric_battery(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) -> ElectricBattery {
        ElectricBattery::new(
            2.,
            0.8,
            3.,
            0.001,
            1.5,
            1.5,
            BatteryLocation::Outside,
            false,
            simulation_time.step,
            Arc::new(external_conditions),
        )
    }

    #[rstest]
    pub fn test_charge_discharge_battery(
        electric_battery: ElectricBattery,
        simulation_time: SimulationTime,
    ) {
        let simulation_time = simulation_time.iter().next().unwrap();
        // supply to battery exceeds limit
        assert_relative_eq!(
            electric_battery.charge_discharge_battery(-1_000., false, simulation_time),
            -1.6770509831248424,
            max_relative = 1e-7
        );
        // demand on battery exceeds limit
        assert_relative_eq!(
            electric_battery.charge_discharge_battery(1_000., false, simulation_time),
            1.121472,
            max_relative = 1e-7
        );
        // normal charge
        assert_relative_eq!(
            electric_battery.charge_discharge_battery(-0.2, false, simulation_time),
            0.,
            max_relative = 1e-7
        );
        // normal discharge
        assert_relative_eq!(
            electric_battery.charge_discharge_battery(0.1, false, simulation_time),
            0.0747648,
            max_relative = 1e-7
        );
    }
}
