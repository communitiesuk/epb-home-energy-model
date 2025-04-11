mod test_air_conditioning {
    use std::sync::Arc;

    use crate::core::controls::time_control::{Control, SetpointTimeControl};
    use crate::core::cooling_systems::air_conditioning::*;
    use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyBuilder};
    use crate::input::FuelType;
    use crate::simulation_time::SimulationTime;
    use parking_lot::RwLock;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 4., 1.)
    }

    #[fixture]
    pub fn aircon(simulation_time: SimulationTime) -> (AirConditioning, Arc<RwLock<EnergySupply>>) {
        let control = SetpointTimeControl::new(
            vec![Some(21.0), Some(21.0), None, Some(21.0)],
            0,
            1.0,
            None,
            None,
            None,
            Default::default(),
            1.0,
        )
        .unwrap();
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
        let energy_supply_conn = EnergySupply::connection(energy_supply.clone(), "aircon").unwrap();
        (
            AirConditioning::new(
                50.,
                2.0,
                0.4,
                energy_supply_conn,
                simulation_time.step,
                Some(Arc::new(Control::SetpointTime(control))),
            ),
            energy_supply,
        )
    }

    #[rstest]
    fn test_demand_energy(
        aircon: (AirConditioning, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
    ) {
        let (aircon, energy_supply) = aircon;
        let inputs = [-40.0, -100.0, -30.0, -20.0];
        let expected_demand = [-40.0, -50.0, 0.0, -20.0];
        let expected_energy_supply_results = [20.0, 25.0, 0.0, 10.0];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                aircon.demand_energy(inputs[t_idx], t_it),
                expected_demand[t_idx],
                "incorrect cooling energy supplied returned"
            );
            assert_eq!(
                energy_supply.read().results_by_end_user()["aircon"][t_idx],
                expected_energy_supply_results[t_idx],
                "incorrect delivered energy demand returned"
            );
        }
    }

    #[rstest]
    fn test_energy_output_min(aircon: (AirConditioning, Arc<RwLock<EnergySupply>>)) {
        let (aircon, _) = aircon;
        assert_eq!(aircon.energy_output_min(), 0.0,);
    }

    #[rstest]
    fn test_temp_setpnt(
        aircon: (AirConditioning, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
    ) {
        let (aircon, _) = aircon;
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                aircon.temp_setpnt(&t_it),
                [Some(21.0), Some(21.0), None, Some(21.0)][t_idx]
            );
        }
    }

    #[rstest]
    fn test_in_required_period(
        aircon: (AirConditioning, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
    ) {
        let (aircon, _) = aircon;
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                aircon.in_required_period(&t_it).unwrap(),
                [true, true, false, true][t_idx]
            );
        }
    }

    #[rstest]
    fn test_frac_convective(aircon: (AirConditioning, Arc<RwLock<EnergySupply>>)) {
        let (aircon, _) = aircon;
        assert_eq!(aircon.frac_convective(), 0.4);
    }
}
