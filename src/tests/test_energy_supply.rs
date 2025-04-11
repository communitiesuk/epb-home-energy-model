mod test_energy_supply {
    use std::sync::atomic::Ordering;
    use std::sync::Arc;

    use crate::core::energy_supply::elec_battery::ElectricBattery;
    use crate::core::energy_supply::energy_supply::*;
    use crate::core::heating_systems::storage_tank::SurplusDiverting;
    use crate::external_conditions::{DaylightSavingsConfig, ExternalConditions};
    use crate::input::{BatteryLocation, FuelType, SecondarySupplyType};
    use crate::simulation_time::{SimulationTime, SimulationTimeIteration};
    use atomic_float::AtomicF64;
    use itertools::Itertools;
    use parking_lot::RwLock;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use serde_json::json;

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0.0, 8.0, 1.0)
    }

    #[fixture]
    pub fn energy_supply<'a>(simulation_time: SimulationTime) -> EnergySupply {
        let mut energy_supply =
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time.total_steps()).build();
        energy_supply.register_end_user_name("shower".to_string());
        energy_supply.register_end_user_name("bath".to_string());

        energy_supply
    }

    #[fixture]
    pub fn energy_supply_connections(
        energy_supply: EnergySupply,
    ) -> (
        EnergySupplyConnection,
        EnergySupplyConnection,
        Arc<RwLock<EnergySupply>>,
    ) {
        let shared_supply = Arc::new(RwLock::new(energy_supply));
        let energy_connection_1 = EnergySupplyConnection {
            energy_supply: shared_supply.clone(),
            end_user_name: "shower".to_string(),
        };
        let energy_connection_2 = EnergySupplyConnection {
            energy_supply: shared_supply.clone(),
            end_user_name: "bath".to_string(),
        };
        (energy_connection_1, energy_connection_2, shared_supply)
    }

    #[fixture]
    pub fn energy_supply_connection_1<'a>(energy_supply: EnergySupply) -> EnergySupplyConnection {
        EnergySupplyConnection {
            energy_supply: Arc::new(RwLock::new(energy_supply)),
            end_user_name: "shower".to_string(),
        }
    }

    #[fixture]
    pub fn energy_supply_connection_2<'a>(energy_supply: EnergySupply) -> EnergySupplyConnection {
        EnergySupplyConnection {
            energy_supply: Arc::new(RwLock::new(energy_supply)),
            end_user_name: "bath".to_string(),
        }
    }

    #[rstest]
    pub fn test_init_demand_list(simulation_time: SimulationTime) {
        assert_eq!(
            init_demand_list(simulation_time.total_steps()),
            [0.; 8].into_iter().map(AtomicF64::new).collect::<Vec<_>>()
        );
    }

    #[rstest]
    fn test_energy_out(energy_supply: EnergySupply, simulation_time: SimulationTime) {
        // Check with existing end user name
        let amount_demand = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0];

        for (t_idx, _) in simulation_time.iter().enumerate() {
            energy_supply
                .energy_out("shower", amount_demand[t_idx], t_idx)
                .unwrap();
            assert_eq!(
                energy_supply.energy_out_by_end_user["shower"][t_idx].load(Ordering::SeqCst),
                amount_demand[t_idx]
            );
        }
        // Check an error is raised with new end user name
        assert!(energy_supply.energy_out("electricshower", 10., 0).is_err());
    }

    #[fixture]
    fn pv_diverter() -> Arc<RwLock<dyn SurplusDiverting>> {
        struct NullDiverter;

        impl SurplusDiverting for NullDiverter {
            fn divert_surplus(
                &self,
                _surplus: f64,
                _simtime: SimulationTimeIteration,
            ) -> anyhow::Result<f64> {
                Ok(0.)
            }
        }

        Arc::new(RwLock::new(NullDiverter))
    }

    #[rstest]
    fn test_connect_diverter(
        mut energy_supply: EnergySupply,
        pv_diverter: Arc<RwLock<dyn SurplusDiverting>>,
    ) {
        assert!(energy_supply.diverter.is_none());
        energy_supply.connect_diverter(pv_diverter.clone()).unwrap();
        assert!(energy_supply.diverter.is_some());
        assert!(energy_supply.connect_diverter(pv_diverter.clone()).is_err());
    }

    #[rstest]
    fn test_demand_energy(energy_supply: EnergySupply, simulation_time: SimulationTime) {
        let amount_demanded = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0];
        for (t_idx, _) in simulation_time.iter().enumerate() {
            energy_supply
                .demand_energy("shower", amount_demanded[t_idx], t_idx)
                .unwrap();
            assert_eq!(
                energy_supply.demand_total[t_idx].load(Ordering::SeqCst),
                amount_demanded[t_idx]
            );
            assert_eq!(
                energy_supply.demand_by_end_user["shower"][t_idx].load(Ordering::SeqCst),
                amount_demanded[t_idx]
            );
        }
    }

    #[rstest]
    fn test_supply_energy(energy_supply: EnergySupply, simulation_time: SimulationTime) {
        let amount_produced = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0];
        for (t_idx, _) in simulation_time.iter().enumerate() {
            energy_supply
                .supply_energy("shower", amount_produced[t_idx], t_idx)
                .unwrap();
            assert_eq!(
                energy_supply.demand_total[t_idx].load(Ordering::SeqCst),
                [-10.0, -20.0, -30.0, -40.0, -50.0, -60.0, -70.0, -80.0][t_idx]
            );
            assert_eq!(
                energy_supply.demand_by_end_user["shower"][t_idx].load(Ordering::SeqCst),
                [-10.0, -20.0, -30.0, -40.0, -50.0, -60.0, -70.0, -80.0][t_idx]
            );
        }
    }

    const EXPECTED_TOTAL_DEMANDS: [f64; 8] =
        [50.0, 120.0, 190.0, 260.0, 330.0, 400.0, 470.0, 540.0];

    #[rstest]
    pub fn test_results_total(energy_supply: EnergySupply, simulation_time: SimulationTime) {
        for simtime in simulation_time.iter() {
            let _ = energy_supply.demand_energy(
                "shower",
                (simtime.index as f64 + 1.0) * 50.0,
                simtime.index,
            );
            let _ = energy_supply.demand_energy("bath", simtime.index as f64 * 20.0, simtime.index);
            assert_eq!(
                energy_supply.results_total()[simtime.index],
                EXPECTED_TOTAL_DEMANDS[simtime.index],
                "incorrect total demand energy returned on iteration {} (1-indexed)",
                simtime.index + 1
            )
        }
    }

    const EXPECTED_TOTAL_DEMANDS_BY_END_USER: [[f64; 8]; 2] = [
        [50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0],
        [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0],
    ];

    #[rstest]
    pub fn test_results_by_end_user(
        energy_supply_connections: (
            EnergySupplyConnection,
            EnergySupplyConnection,
            Arc<RwLock<EnergySupply>>,
        ),
        simulation_time: SimulationTime,
    ) {
        let (energy_connection_1, energy_connection_2, energy_supply) = energy_supply_connections;
        for simtime in simulation_time.iter() {
            let _ = energy_connection_1
                .demand_energy((simtime.index as f64 + 1.0) * 50.0, simtime.index);
            let _ = energy_connection_2.demand_energy(simtime.index as f64 * 20.0, simtime.index);
            assert_eq!(
                energy_supply.read().results_by_end_user()["shower"][simtime.index],
                EXPECTED_TOTAL_DEMANDS_BY_END_USER[0][simtime.index]
            );
            assert_eq!(
                energy_supply.read().results_by_end_user()["bath"][simtime.index],
                EXPECTED_TOTAL_DEMANDS_BY_END_USER[1][simtime.index]
            );
        }
    }

    const EXPECTED_BETA_FACTORS: [f64; 8] = [
        1.0,
        0.8973610789278808,
        0.4677549807236648,
        0.3297589507351858,
        0.2578125,
        0.2,
        0.16319444444444445,
        0.1377551020408163,
    ];
    const EXPECTED_SURPLUSES: [f64; 8] = [
        0.0,
        -8.21111368576954,
        -170.3184061684273,
        -482.57355547066624,
        -950.0,
        -1600.0,
        -2410.0,
        -3380.0,
    ];
    const EXPECTED_DEMANDS_NOT_MET: [f64; 8] = [
        50.0,
        48.21111368576953,
        40.31840616842726,
        22.573555470666236,
        0.0,
        0.0,
        0.0,
        0.0,
    ];

    #[rstest]
    pub fn test_beta_factor(
        energy_supply_connections: (
            EnergySupplyConnection,
            EnergySupplyConnection,
            Arc<RwLock<EnergySupply>>,
        ),
        simulation_time: SimulationTime,
    ) {
        let (energy_connection_1, energy_connection_2, energy_supply) = energy_supply_connections;
        let energy_connection_3 = EnergySupply::connection(energy_supply.clone(), "PV").unwrap();
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            energy_connection_1
                .demand_energy((t_idx as f64 + 1.) * 50., t_idx)
                .unwrap();
            energy_connection_2
                .demand_energy(t_idx as f64 * 20., t_idx)
                .unwrap();
            energy_connection_3
                .supply_energy(t_idx as f64 * t_idx as f64 * 80., t_idx)
                .unwrap();

            let energy_supply = energy_supply.read();
            energy_supply
                .calc_energy_import_export_betafactor(t_it)
                .unwrap();

            assert_eq!(
                energy_supply.get_beta_factor()[t_idx],
                EXPECTED_BETA_FACTORS[t_idx],
                "incorrect beta factor returned"
            );
            assert_eq!(
                energy_supply.get_energy_export()[t_idx],
                EXPECTED_SURPLUSES[t_idx],
                "incorrect energy export returned"
            );
            assert_eq!(
                energy_supply.get_energy_import()[t_idx],
                EXPECTED_DEMANDS_NOT_MET[t_idx],
                "incorrect energy import returned"
            );
        }
    }

    #[fixture]
    fn external_conditions(simulation_time: SimulationTime) -> ExternalConditions {
        ExternalConditions::new(
            &simulation_time.iter(),
            vec![0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0],
            vec![3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1],
            vec![0., 20., 40., 60., 0., 20., 40., 60.],
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
            serde_json::from_value(json!([
                // upstream Python gives old 'start' fields, but we need 'start360' here
                {"number": 1, "start360": 0, "end360": 45},
                {"number": 2, "start360": 45, "end360": 90,
                 "shading": [
                     {"type": "overhang", "height": 2.2, "distance": 6}
                     ]
                 },
                {"number": 3, "start360": 90, "end360": 135},
                {"number": 4, "start360": 135, "end360": 180,
                 "shading": [
                     {"type": "obstacle", "height": 40, "distance": 4},
                     {"type": "overhang", "height": 3, "distance": 7}
                     ]
                 },
                {"number": 5, "start360": 180, "end360": 225,
                 "shading": [
                     {"type": "obstacle", "height": 3, "distance": 8},
                     ]
                 },
                {"number": 6, "start360": 225, "end360": 270},
                {"number": 7, "start360": 270, "end360": 315},
                {"number": 8, "start360": 315, "end360": 360}
            ]))
            .unwrap(),
        )
    }

    #[rstest]
    fn test_calc_energy_import_export_betafactor(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        let amount_demanded = [50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0];
        let amount_produced = [50.0, 90.0, 130.0, 210.0, 2300.0, 290.0, 300.0, 350.0];

        let elec_battery = ElectricBattery::new(
            2.,
            0.8,
            3.,
            0.001,
            1.5,
            1.5,
            BatteryLocation::Outside,
            false,
            simulation_time.step,
            Arc::new(external_conditions.clone()),
        );

        let builder =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps());
        let energy_supply = builder.with_electric_battery(elec_battery).build();

        let energy_supply = Arc::new(RwLock::new(energy_supply));

        // test with elec battery
        let _shower_connection = EnergySupply::connection(energy_supply.clone(), "shower").unwrap();
        let _bath_connection = EnergySupply::connection(energy_supply.clone(), "bath").unwrap();

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let energy_supply = energy_supply.read();
            energy_supply
                .demand_energy("shower", amount_demanded[t_idx], t_idx)
                .unwrap();
            energy_supply
                .supply_energy("bath", amount_produced[t_idx], t_idx)
                .unwrap();
            energy_supply
                .calc_energy_import_export_betafactor(t_it)
                .unwrap();
        }

        {
            let energy_supply = energy_supply.read();

            assert_eq!(
                energy_supply
                    .demand_total
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![0., 10., 20., -10., -2050., 10., 50., 50.]
            );

            assert_eq!(
                energy_supply
                    .demand_not_met
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    15.138528000000004,
                    34.37872423953049,
                    52.991809292243516,
                    63.07009986960006,
                    -2.842170943040401e-14,
                    99.5880926222444,
                    124.3891811736907,
                    140.57522007095577
                ]
            );

            assert_eq!(
                energy_supply
                    .supply_surplus
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    -14.582949016875158,
                    -24.59889302603037,
                    -32.991809292243516,
                    -73.07009986960004,
                    -2050.,
                    -89.5880926222444,
                    -74.38918117369072,
                    -90.5752200709558
                ]
            );

            assert_eq!(
                energy_supply
                    .energy_generated_consumed
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    33.739999999999995,
                    65.40110697396963,
                    97.00819070775648,
                    136.92990013039994,
                    250.00000000000003,
                    200.4119073777556,
                    225.6108188263093,
                    259.42477992904423
                ]
            );

            assert_eq!(
                energy_supply
                    .energy_into_battery_from_generation
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![1.6770509831248424, -0., -0., -0., -0., -0., -0., -0.]
            );

            assert_eq!(
                energy_supply
                    .energy_out_of_battery
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![-1.121472, -0.2201687864998738, -0., -0., 0., -0., -0., -0.]
            );

            assert_eq!(
                energy_supply
                    .energy_diverted
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![0., 0., 0., 0., 0., 0., 0., 0.]
            );

            assert_eq!(
                energy_supply
                    .beta_factor
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    0.6748,
                    0.7266789663774403,
                    0.7462168515981268,
                    0.652047143478095,
                    0.10869565217391305,
                    0.6910755426819158,
                    0.7520360627543643,
                    0.7412136569401263
                ]
            );
        }

        // Test with PV diverter
        struct MockDiverter;

        impl SurplusDiverting for MockDiverter {
            fn divert_surplus(
                &self,
                _supply_surplus: f64,
                _simulation_time_iteration: SimulationTimeIteration,
            ) -> anyhow::Result<f64> {
                Ok(10.)
            }
        }

        let diverter = Arc::new(RwLock::new(MockDiverter));
        energy_supply.write().connect_diverter(diverter).unwrap();

        for (t_idx, simtime) in simulation_time.iter().enumerate() {
            let energy_supply = energy_supply.read();

            energy_supply
                .demand_energy("shower", amount_demanded[t_idx], t_idx)
                .unwrap();
            energy_supply
                .supply_energy("bath", amount_produced[t_idx], t_idx)
                .unwrap();
            energy_supply
                .calc_energy_import_export_betafactor(simtime)
                .unwrap();
        }

        {
            let energy_supply = energy_supply.read();

            assert_eq!(
                energy_supply
                    .demand_total
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![0., 20., 40., -20., -4100., 20., 100., 100.]
            );

            assert_eq!(
                energy_supply
                    .demand_not_met
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    47.65852800000002,
                    103.57651029159123,
                    158.97542787673055,
                    189.21029960880017,
                    -8.526512829121202e-14,
                    298.7642778667332,
                    373.1675435210721,
                    421.7256602128673
                ]
            );

            assert_eq!(
                energy_supply
                    .supply_surplus
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    -37.10294901687516,
                    -63.79667907809112,
                    -88.97542787673055,
                    -209.21029960880014,
                    -6140.0,
                    -258.7642778667332,
                    -213.16754352107216,
                    -261.7256602128674
                ]
            );

            assert_eq!(
                energy_supply
                    .energy_generated_consumed
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    101.21999999999998,
                    196.2033209219089,
                    291.0245721232694,
                    410.78970039119986,
                    750.0000000000001,
                    601.2357221332668,
                    676.832456478928,
                    778.2743397871327
                ]
            );

            assert_eq!(
                energy_supply
                    .energy_into_battery_from_generation
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![-0., -0., -0., -0., -0., -0., -0., -0.]
            );

            assert_eq!(
                energy_supply
                    .energy_out_of_battery
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![-0., -0., -0., -0., 0., -0., -0., -0.]
            );

            assert_eq!(
                energy_supply
                    .energy_diverted
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![10., 10., 10., 10., 10., 10., 10., 10.]
            );

            assert_eq!(
                energy_supply
                    .beta_factor
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    0.6748,
                    0.7266789663774403,
                    0.7462168515981268,
                    0.652047143478095,
                    0.10869565217391305,
                    0.6910755426819158,
                    0.7520360627543643,
                    0.7412136569401263
                ]
            );
        }

        // important so the energy supply Arc only has one strong reference to the energy supply
        // below where Arc::into_inner is called, otherwise that call would fail
        drop(_shower_connection);
        drop(_bath_connection);

        // LOOK AWAY 👀👀👀👀👀
        // (the upstream Python shared the same electric battery across energy supplies in this test,
        // so its internal state is not isolated - therefore we need to cannibalise the previous energy
        // supply here for scraps (the electric battery) for use in the next set of assertions)
        let elec_battery = Arc::into_inner(energy_supply)
            .unwrap()
            .into_inner()
            .electric_battery
            .unwrap();

        // Set priority
        let priority = vec![
            SecondarySupplyType::Diverter,
            SecondarySupplyType::ElectricBattery,
        ];

        let mut builder =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps());
        builder = builder
            .with_electric_battery(elec_battery)
            .with_priority(priority);

        let energy_supply = Arc::new(RwLock::new(builder.build()));

        let _shower_connection = EnergySupply::connection(energy_supply.clone(), "shower").unwrap();
        let _bath_connection = EnergySupply::connection(energy_supply.clone(), "bath").unwrap();

        let diverter = Arc::new(RwLock::new(MockDiverter));
        energy_supply.write().connect_diverter(diverter).unwrap();

        for (t_idx, simtime) in simulation_time.iter().enumerate() {
            let energy_supply = energy_supply.read();

            energy_supply
                .demand_energy("shower", amount_demanded[t_idx], t_idx)
                .unwrap();
            energy_supply
                .supply_energy("bath", amount_produced[t_idx], t_idx)
                .unwrap();
            energy_supply
                .calc_energy_import_export_betafactor(simtime)
                .unwrap();
        }

        {
            let energy_supply = energy_supply.read();

            assert_eq!(
                energy_supply
                    .demand_total
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![0., 10., 20., -10., -2050., 10., 50., 50.]
            );

            assert_eq!(
                energy_supply
                    .demand_not_met
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    16.260000000000005,
                    34.59889302603037,
                    52.991809292243516,
                    63.07009986960006,
                    -2.842170943040401e-14,
                    99.5880926222444,
                    124.3891811736907,
                    140.57522007095577
                ]
            );

            assert_eq!(
                energy_supply
                    .supply_surplus
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    -6.260000000000002,
                    -14.598893026030371,
                    -22.991809292243516,
                    -63.07009986960004,
                    -2040.0,
                    -79.5880926222444,
                    -64.38918117369072,
                    -80.5752200709558
                ]
            );

            assert_eq!(
                energy_supply
                    .energy_generated_consumed
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    33.739999999999995,
                    65.40110697396963,
                    97.00819070775648,
                    136.92990013039994,
                    250.00000000000003,
                    200.4119073777556,
                    225.6108188263093,
                    259.42477992904423
                ]
            );

            assert_eq!(
                energy_supply
                    .energy_into_battery_from_generation
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![-0., -0., -0., -0., -0., -0., -0., -0.]
            );

            assert_eq!(
                energy_supply
                    .energy_out_of_battery
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![-0., -0., -0., -0., -0., -0., -0., -0.]
            );

            assert_eq!(
                energy_supply
                    .energy_diverted
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![10., 10., 10., 10., 10., 10., 10., 10.]
            );

            assert_eq!(
                energy_supply
                    .beta_factor
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    0.6748,
                    0.7266789663774403,
                    0.7462168515981268,
                    0.652047143478095,
                    0.10869565217391305,
                    0.6910755426819158,
                    0.7520360627543643,
                    0.7412136569401263
                ]
            );
        }
    }

    #[rstest]
    pub fn test_energy_supply_without_export(simulation_time: SimulationTime) {
        let mut builder =
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time.total_steps());
        builder = builder.with_export_capable(false);
        let energy_supply = builder.build();
        let shared_supply = Arc::new(RwLock::new(energy_supply));
        let energy_connection_1 =
            EnergySupply::connection(shared_supply.clone(), "shower").unwrap();
        let energy_connection_2 = EnergySupply::connection(shared_supply.clone(), "bath").unwrap();
        for t_it in simulation_time.iter() {
            let t_idx = t_it.index;
            energy_connection_1
                .demand_energy(((t_idx + 1) * 50) as f64, t_idx)
                .unwrap();
            energy_connection_2
                .demand_energy((t_idx * 20) as f64, t_idx)
                .unwrap();
            assert_eq!(
                shared_supply.read().get_energy_export()[t_idx],
                0.,
                "incorrect energy export returned"
            );
        }
    }
}
