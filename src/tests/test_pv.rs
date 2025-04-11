mod test_pv {
    use crate::core::energy_supply::pv::*;
    use crate::external_conditions::{ExternalConditions, WindowShadingObject};
    use crate::input::{InverterType, OnSiteGenerationVentilationStrategy};
    use std::sync::Arc;

    use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyBuilder};
    use crate::external_conditions::{
        DaylightSavingsConfig, ShadingObject, ShadingObjectType, ShadingSegment,
    };
    use crate::input::FuelType;
    use crate::simulation_time::SimulationTime;
    use approx::assert_relative_eq;
    use parking_lot::RwLock;
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
            vec![220., 230., 240., 250., 260., 270., 270., 280.],
            vec![11., 25., 42., 52., 60., 44., 28., 15.],
            vec![11., 25., 42., 52., 60., 44., 28., 15.],
            vec![0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2],
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
                    ..Default::default()
                },
                ShadingSegment {
                    number: 7,
                    start: -90.,
                    end: -135.,
                    ..Default::default()
                },
                ShadingSegment {
                    number: 8,
                    start: -135.,
                    end: -180.,
                    ..Default::default()
                },
            ]
            .into(),
        )
    }

    #[fixture]
    pub fn pv(
        simulation_time: SimulationTime,
        external_conditions: ExternalConditions,
    ) -> (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>) {
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
        let energy_supply_conn =
            EnergySupply::connection(energy_supply.clone(), "pv generation without shading")
                .unwrap();
        let pv = PhotovoltaicSystem::new(
            2.5,
            OnSiteGenerationVentilationStrategy::ModeratelyVentilated,
            30.,
            0.,
            10.,
            2.,
            3.,
            Arc::new(external_conditions),
            energy_supply_conn,
            simulation_time.step,
            vec![],
            2.5,
            0.05,
            false,
            InverterType::OptimisedInverter,
        );
        (pv, energy_supply)
    }

    #[fixture]
    fn pv_with_shading(
        simulation_time: SimulationTime,
        external_conditions: ExternalConditions,
    ) -> (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>) {
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
        let energy_supply_conn =
            EnergySupply::connection(energy_supply.clone(), "pv generation with shading").unwrap();
        let pv = PhotovoltaicSystem::new(
            2.5,
            OnSiteGenerationVentilationStrategy::ModeratelyVentilated,
            30.,
            0.,
            10.,
            2.,
            3.,
            Arc::new(external_conditions),
            energy_supply_conn,
            simulation_time.step,
            vec![
                WindowShadingObject::Overhang {
                    depth: 0.5,
                    distance: 0.5,
                },
                WindowShadingObject::SideFinLeft {
                    depth: 0.25,
                    distance: 0.1,
                },
                WindowShadingObject::SideFinRight {
                    depth: 0.25,
                    distance: 0.1,
                },
            ],
            2.5,
            0.02,
            true,
            InverterType::StringInverter,
        );
        (pv, energy_supply)
    }

    #[rstest]
    fn test_is_inside(
        pv: (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>),
        pv_with_shading: (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>),
    ) {
        let (pv, _) = pv;
        let (pv_with_shading, _) = pv_with_shading;
        assert!(!pv.inverter_is_inside());
        assert!(pv_with_shading.inverter_is_inside());
    }

    #[rstest]
    fn test_produce_energy(
        pv: (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
    ) {
        let (pv, energy_supply) = pv;
        let expected_generation_results = [
            -0.002911179810082315,
            -0.01585915973389526,
            -0.03631681332778666,
            -0.0462218185635626,
            -0.05,
            -0.03841528069730012,
            -0.019985927280177524,
            -0.014819433057321862,
        ];
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            pv.produce_energy(t_it);
            assert_relative_eq!(
                energy_supply.read().results_by_end_user()["pv generation without shading"][t_idx],
                expected_generation_results[t_idx],
                max_relative = 1e-6
            );
        }
    }

    #[rstest]
    fn test_produce_energy_with_shading(
        pv_with_shading: (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
    ) {
        let (pv, energy_supply) = pv_with_shading;
        let expected_generation_results = [
            -0.0015507260447403823,
            -0.008732593505451556,
            -0.02,
            -0.02,
            -0.02,
            -0.013971934370722397,
            -0.006779589823217942,
            -0.007020372160065822,
        ];
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            pv.produce_energy(t_it);
            assert_relative_eq!(
                energy_supply.read().results_by_end_user()["pv generation with shading"][t_idx],
                expected_generation_results[t_idx],
                max_relative = 1e-6
            );
        }
    }

    #[rstest]
    fn test_energy_produced_and_energy_lost(
        pv: (PhotovoltaicSystem, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
    ) {
        let (pv, _) = pv;
        let expected_energy_produced = [
            0.002911179810082315,
            0.01585915973389526,
            0.03631681332778666,
            0.0462218185635626,
            0.05,
            0.03841528069730012,
            0.019985927280177524,
            0.014819433057321862,
        ];
        let expected_energy_lost = [
            0.012982394066666526,
            0.022261865114937766,
            0.024010573568482414,
            0.023376882776544372,
            0.02560556391283768,
            0.02392305675244094,
            0.023189444399645136,
            0.021949092776458276,
        ];
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let (energy_produced, energy_lost) = pv.produce_energy(t_it);
            assert_relative_eq!(
                energy_produced,
                expected_energy_produced[t_idx],
                max_relative = 1e-6
            );
            assert_relative_eq!(
                energy_lost,
                expected_energy_lost[t_idx],
                max_relative = 1e-6
            );
        }
    }
}
