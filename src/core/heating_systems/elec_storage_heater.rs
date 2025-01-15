use itertools::Itertools;

use crate::{
    core::{
        controls::time_control::{ChargeControl, SetpointTimeControl},
        energy_supply::energy_supply::EnergySupplyConnection,
        space_heat_demand::zone::Zone,
    },
    external_conditions::ExternalConditions,
    input::ElectricStorageHeaterAirFlowType,
    simulation_time::{SimulationTime, SimulationTimeIteration},
};

#[derive(Debug)]
pub struct ElecStorageHeater {
    pwr_in: f64,
    pwr_instant: f64,
    storage_capacity: f64,
    air_flow_type: ElectricStorageHeaterAirFlowType,
    frac_convective: f64,
    n_units: i32,
    zone: Zone,
    energy_supply_conn: EnergySupplyConnection,
    simulation_time: SimulationTime,
    control: SetpointTimeControl,
    charge_control: ChargeControl,
    fan_pwr: f64,
    external_conditions: ExternalConditions,
    temp_air: f64,
    state_of_charge: f64,
    esh_min_output: Vec<(f64, f64)>,
    esh_max_output: Vec<(f64, f64)>,
    demand_met: f64,
    demand_unmet: f64,
    zone_setpoint_init: f64,
    soc_max_array: Vec<f64>,
    power_max_array: Vec<f64>,
    soc_min_array: Vec<f64>,
    power_min_array: Vec<f64>,
    // TODO review - do we need to keep these as public properties?
    pub energy_for_fan: f64,
    pub energy_instant: f64,
    pub energy_charged: f64,
    pub energy_delivered: f64,
}

impl ElecStorageHeater {
    pub fn new(
        pwr_in: f64,
        rated_power_instant: f64,
        storage_capacity: f64,
        air_flow_type: ElectricStorageHeaterAirFlowType,
        frac_convective: f64,
        fan_pwr: f64,
        n_units: i32,
        zone: Zone,
        energy_supply_conn: EnergySupplyConnection,
        simulation_time: SimulationTime,
        control: SetpointTimeControl,
        charge_control: ChargeControl,
        esh_min_output: Vec<(f64, f64)>,
        esh_max_output: Vec<(f64, f64)>,
        ext_cond: ExternalConditions,
    ) -> Self {
        // Arguments:
        // pwr_in               -- in kW (Charging)
        // rated_power_instant  -- in kW (Instant backup)
        // storage_capacity     -- in kWh
        // air_flow_type        -- string specifying type of Electric Storage Heater:
        //                      -- fan-assisted
        //                      -- damper-only
        // frac_convective      -- convective fraction for heating (TODO: Check if necessary)
        // fan_pwr              -- Fan power [W]
        // n_units              -- number of units install in zone
        // zone                 -- zone where the unit(s) is/are installed
        // energy_supply_conn   -- reference to EnergySupplyConnection object
        // simulation_time      -- reference to SimulationTime object
        // control              -- reference to a control object which must implement is_on() and setpnt() funcs
        // charge_control       -- reference to a ChargeControl object which must implement different logic types
        //                         for charging the Electric Storage Heaters.
        // ESH_min_output       -- Data from test showing the output from the storage heater when not actively
        //                         outputting heat, i.e. case losses only (with units kW)
        // ESH_max_output       -- Data from test showing the output from the storage heater when it is actively
        //                         outputting heat, e.g. damper open / fan running (with units kW)
        // extcond              -- reference to ExternalConditions object

        let zone_setpoint_init = zone.setpnt_init();
        let temp_air = zone.temp_internal_air();

        // Convert ESH_max_output to NumPy arrays without sorting
        let soc_max_array = esh_max_output.iter().map(|f| f.0).collect_vec();
        let power_max_array = esh_max_output.iter().map(|f| f.1).collect_vec();

        // Convert ESH_min_output to NumPy arrays without sorting
        let soc_min_array = esh_min_output.iter().map(|f| f.0).collect_vec();
        let power_min_array = esh_min_output.iter().map(|f| f.1).collect_vec();

        // Validate that both SOC arrays start at 0.0 and end at 1.0
        // TODO Result or more specific panic
        if !is_close!(*soc_max_array.first().unwrap(), 0.) {
            panic!("The first SOC value in esh_max_output must be 0.0 (fully discharged).");
        }

        if !is_close!(*soc_max_array.last().unwrap(), 1.) {
            panic!("The last SOC value in esh_max_output must be 1.0 (fully charged).");
        }

        if !is_close!(*soc_min_array.first().unwrap(), 0.) {
            panic!("The first SOC value in esh_min_output must be 0.0 (fully discharged).");
        }

        if !is_close!(*soc_min_array.last().unwrap(), 1.) {
            panic!("The last SOC value in esh_min_output must be 1.0 (fully charged).");
        }
        
        Self {
            pwr_in,
            pwr_instant: rated_power_instant,
            storage_capacity,
            air_flow_type,
            frac_convective,
            n_units,
            zone,
            energy_supply_conn,
            simulation_time,
            control,
            charge_control,
            fan_pwr,
            external_conditions: ext_cond,
            temp_air,
            state_of_charge: 0.,
            esh_min_output,
            esh_max_output,
            demand_met: 0.,
            demand_unmet: 0.,
            zone_setpoint_init,
            soc_max_array,
            power_max_array,
            soc_min_array,
            power_min_array,
            // TODO ...
            energy_for_fan: 0.,
            energy_instant: 0.,
            energy_charged: 0.,
            energy_delivered: 0.
        }
    }

    pub fn energy_output_min(&self, simulation_time_iteration: &SimulationTimeIteration) -> f64 {
        todo!()
    }

    pub fn energy_output_max(&self, simulation_time_iteration: &SimulationTimeIteration) -> (f64, f64, f64, f64) {
        todo!()
    }
    
    pub fn demand_energy(&self, energy_demand: f64, simulation_time_iteration: &SimulationTimeIteration) -> f64 {
        todo!()
    }

    pub fn target_electric_charge(&self, simulation_time_iteration: &SimulationTimeIteration) -> f64 {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        core::{
            controls::time_control::SetpointTimeControl,
            space_heat_demand::zone::Zone,
        },
        external_conditions::{DaylightSavingsConfig, ExternalConditions},
        input::{ControlLogicType, ExternalSensor},
        simulation_time::{SimulationTime, SimulationTimeIteration, SimulationTimeIterator},
    };
    use approx::assert_relative_eq;
    use rstest::{fixture, rstest};
    use serde_json::json;
    use super::*;

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 2., 1.)
    }

    #[fixture]
    pub fn simulation_time_iterator(
        simulation_time: SimulationTime,
    ) -> SimulationTimeIterator {
        simulation_time.iter()
    }

    #[fixture]
    pub fn simulation_time_iteration(
        simulation_time_iterator: SimulationTimeIterator,
    ) -> SimulationTimeIteration {
        simulation_time_iterator.current_iteration()
    }

    #[fixture]
    fn external_conditions(simulation_time: SimulationTime) -> ExternalConditions {
        ExternalConditions::new(
            &simulation_time.iter(),
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
            // TODO check these are correct
            serde_json::from_value(json!(
                [
                    {"number": 1, "start360": 180, "end360": 135},
                    {"number": 2, "start360": 135, "end360": 90},
                    {"number": 3, "start360": 90, "end360": 45},
                    {"number": 4, "start360": 45, "end360": 0,
                        "shading": [
                            {"type": "obstacle", "height": 10.5, "distance": 12}
                        ]
                    },
                    {"number": 5, "start360": 0, "end360": -45},
                    {"number": 6, "start360": -45, "end360": -90},
                    {"number": 7, "start360": -90, "end360": -135},
                    {"number": 8, "start360": -135, "end360": -180}
                ]
            ))
            .unwrap(),
        )
    }

    #[fixture]
    fn external_sensor() -> ExternalSensor {
        serde_json::from_value(json!({
            "correlation": [
                {"temperature": 0.0, "max_charge": 1.0},
                {"temperature": 10.0, "max_charge": 0.9},
                {"temperature": 18.0, "max_charge": 0.5}
            ]
        }))
        .unwrap()
    }

    #[fixture]
    fn charge_control(
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
    ) -> ChargeControl {
        ChargeControl::new(
            ControlLogicType::Automatic,
            vec![
                true, true, true, true, true, true, true, true, false, false, false, false, false,
                false, false, false, true, true, true, true, false, false, false, false,
            ],
            1.,
            0,
            1.,
            vec![1.0, 0.8],
            Some(22.),
            None,
            None,
            None,
            external_conditions.into(),
            Some(external_sensor),
        )
        .unwrap()
    }

    #[fixture]
    fn elec_storage_heater(
        simulation_time: SimulationTime,
        charge_control: ChargeControl,
        external_conditions: ExternalConditions,
    ) -> ElecStorageHeater {
        let zone: Zone = todo!();
        let energy_supply: EnergySupplyConnection = todo!();
        let control: SetpointTimeControl = todo!();

        let esh_min_output = vec![(0.0, 0.0), (0.5, 0.02), (1.0, 0.05)];
        let esh_max_output = vec![(0.0, 0.0), (0.5, 1.5), (1.0, 3.0)];

        ElecStorageHeater::new(
            3.5,
            2.5,
            10.0,
            ElectricStorageHeaterAirFlowType::FanAssisted,
            0.7,
            11.,
            1,
            zone,
            energy_supply,
            simulation_time,
            control,
            charge_control,
            esh_min_output,
            esh_max_output,
            external_conditions, // TODO this is None in Python
        )
    }

    #[rstest]
    #[ignore = "not yet implemented"]
    fn test_energy_output_min(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: ElecStorageHeater
    ) {
        // Test minimum energy output calculation across all timesteps.

        let expected_min_energy_output = [
            0.019999999999999997,
            0.030419151282454364,
            0.0406826518009459,
            0.046014105108884346,
            0.04674323706081289,
            0.045800000000000014,
            0.046400000000000004,
            0.038,
            0.046886897763894056,
            0.03215061095009046,
            0.021233700713726503,
            0.01542628227371725,
            0.011428072222071864,
            0.008466125322130943,
            0.006271860545638567,
            0.004646308882570838,
            0.010432746398675731,
            0.01897277720547578,
            0.019999999999999997,
            0.019999999999999997,
            0.020056174317410178,
            0.014844117518306485,
            0.010996794239857175,
            0.008146626650951819,
        ]; // Actual minimum energy output for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let min_energy_output = elec_storage_heater.energy_output_min(&t_it);

            // TODO is this line needed?
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);

            assert_relative_eq!(min_energy_output, expected_min_energy_output[t_idx]);
        }
        assert!(false);
    }

    #[rstest]
    #[ignore = "not yet implemented"]
    fn test_energy_output_max(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: ElecStorageHeater
    ) {
        // Test maximum energy output calculation across all timesteps.
        let expected_max_energy_output = [
            1.5, 1.772121660521405, 2.2199562136927717, 2.5517202117781994,
            2.7913851590672585, 2.7899999999999996, 2.8200000000000003, 2.4000000000000004,
            2.463423313846487, 1.8249489529640162, 1.3519554011630448, 1.0015529506734968,
            0.7419686857505708, 0.5496640579327374, 0.40720123344887615, 0.30166213975932143,
            0.6996897293886958, 1.3814284569589004, 1.5, 1.5,
            1.3009346098448467, 0.9637557636923015, 0.713967931076402, 0.5289205810615784
        ]; // Expected max energy output for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let max_energy_output = elec_storage_heater.energy_output_max(&t_it);

            // TODO is this line needed
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);
            let (energy, _, _, _) = max_energy_output;

            assert_relative_eq!(energy, expected_max_energy_output[t_idx]);
        }
    }

    #[rstest]
    #[ignore = "not yet implemented"]
    fn test_electric_charge(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: ElecStorageHeater) {
            // Test electric charge calculation across all timesteps.
            let expected_target_elec_charge = [
            0.5, 1.0, 0.99, 0.98,
            0.95, 0.93, 0.9400000000000001, 0.8,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.5, 0.5, 0.5, 0.5,
            0.0, 0.0, 0.0, 0.0 
        ]; // Expected target charge for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let target_elec_charge = elec_storage_heater.target_electric_charge(&t_it);
            assert_relative_eq!(target_elec_charge, expected_target_elec_charge[t_idx]);
        }
    }

    #[rstest]
    #[ignore = "not yet implemented"]
    pub fn test_demand_energy(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: ElecStorageHeater
    ) {
        let expected_energy = [
            4.0, 4.272121660521405, 4.719956213692772, 5.0,
            5.0, 5.0, 5.0, 4.9,
            4.963423313846487, 4.324948952964016, 3.851955401163045, 3.5015529506734966,
            3.241968685750571, 3.0496640579327376, 2.907201233448876, 2.8016621397593213,
            3.199689729388696, 3.8814284569589006, 4.0, 4.0,
            3.8009346098448464, 3.4637557636923013, 3.213967931076402, 3.0289205810615782
        ]; // Expected energy for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let energy_out = elec_storage_heater.demand_energy(5.0, &t_it);
            assert_relative_eq!(energy_out, expected_energy[t_idx]);
        }
    }

    #[rstest]
    #[ignore = "not yet implemented"]
    pub fn test_energy_for_fan(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: ElecStorageHeater
    ) {
        let expected_energy_for_fan = [
            0.003666666666666666, 0.002707621094790285, 0.0019580976410457644, 0.0018707482993197276,
            0.0019298245614035089, 0.0019713261648745518, 0.001950354609929078, 0.0022916666666666662,
            0.0021220628632204496, 0.002233750362976654, 0.0023578475867206904, 0.0024965444862427343,
            0.0026525785014320812, 0.0028294170564121318, 0.003031518268419362, 0.0032647119841350417,
            0.003666666666666666, 0.003666666666666666, 0.003666666666666666, 0.003666666666666666,
            0.0021220628632204505, 0.0022337503629766544, 0.0023578475867206913, 0.002496544486242736
        ]; // Expected energy for fan for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);
            let energy_for_fan = elec_storage_heater.energy_for_fan;
            assert_relative_eq!(energy_for_fan, expected_energy_for_fan[t_idx]);
        }
    }

    #[rstest]
    #[ignore = "not yet implemented"]
    pub fn test_energy_instant(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: ElecStorageHeater
    ) {
        let expected_energy_instant = [
            2.5, 2.5, 2.5, 2.4482797882218006,
            2.208614840932741, 2.2100000000000004, 2.18, 2.5,
            2.5, 2.5, 2.5, 2.5,
            2.5, 2.5, 2.5, 2.5,
            2.5, 2.5, 2.5, 2.5,
            2.5, 2.5, 2.5, 2.5
        ]; // Expected backup energy instant for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);
            let energy_instant = elec_storage_heater.energy_instant;
            assert_relative_eq!(energy_instant, expected_energy_instant[t_idx]);
        }
    }

    #[rstest]
    #[ignore = "not yet implemented"]
    pub fn test_energy_charged(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: ElecStorageHeater
    ) {
        let expected_energy_charged = [
            1.5, 3.500000000000001, 3.5, 3.5,
            3.3397997646920756, 2.7899999999999996, 2.82, 2.4000000000000004,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            3.500000000000001, 2.738269398472182, 1.5, 1.5,
            0.0, 0.0, 0.0, 0.0
        ]; // Expected energy charged for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);
            let energy_charged = elec_storage_heater.energy_charged;
            assert_relative_eq!(energy_charged, expected_energy_charged[t_idx]);
        }
    }

    #[rstest]
    #[ignore = "not yet implemented"]
    pub fn test_energy_stored_delivered(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: ElecStorageHeater
    )
    {
        let expected_energy_delivered = [
            1.5, 1.772121660521405, 2.219956213692772, 2.5517202117781994,
            2.791385159067259, 2.7899999999999996, 2.82, 2.4000000000000004,
            2.4634233138464876, 1.8249489529640166, 1.3519554011630448, 1.0015529506734968,
            0.7419686857505706, 0.5496640579327375, 0.4072012334488761, 0.30166213975932155,
            0.6996897293886956, 1.3814284569589008, 1.5, 1.5,
            1.300934609844847, 0.9637557636923016, 0.713967931076402, 0.5289205810615786
        ]; // Expected energy stored delivered for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);
            let energy_delivered = elec_storage_heater.energy_delivered;
            assert_relative_eq!(energy_delivered, expected_energy_delivered[t_idx]);
        }
    }
}
