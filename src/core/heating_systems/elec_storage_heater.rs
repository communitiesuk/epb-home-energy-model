use crate::{
    core::{
        controls::time_control::{ChargeControl, Control, SetpointTimeControl},
        energy_supply::energy_supply::EnergySupplyConnection,
        space_heat_demand::zone::Zone,
    },
    external_conditions::ExternalConditions,
    input::ElectricStorageHeaterAirFlowType,
    simulation_time::SimulationTime,
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
    state_of_charge: f64,
    esh_min_output: Vec<(f64, f64)>,
    esh_max_output: Vec<(f64, f64)>,
    demand_met: f64,
    demand_unmet: f64,
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
            state_of_charge: 0.,
            esh_min_output,
            esh_max_output,
            demand_met: 0.,
            demand_unmet: 0.,
            // zone_setpnt_init: // TODO
            // soc_max_array
            // ower_max_array
            // soc_min_array
            // power_min_array
            // TODO ...
        }
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use crate::{
        core::{
            controls::time_control::SetpointTimeControl,
            space_heat_demand::zone::Zone,
        },
        external_conditions::{DaylightSavingsConfig, ExternalConditions},
        input::{ControlLogicType, ExternalSensor},
        simulation_time::{SimulationTime, SimulationTimeIteration, SimulationTimeIterator},
    };
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
    fn test_energy_output_min(
        simulation_time_iterator: SimulationTimeIterator,
        external_conditions: ExternalConditions,
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

        for t_idx in simulation_time_iterator.enumerate() {
            // TODO
        }
        assert!(false);
    }
}
