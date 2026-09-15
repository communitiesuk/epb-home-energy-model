use crate::compare_floats::min_of_2;
use crate::core::common::WaterSupply;
use crate::core::controls::time_control::{Control, RangeTimeControl};
use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::core::heating_systems::boiler::{
    BoilerForBoilerService, BoilerServiceSpace, BoilerServiceWaterCombi, BoilerServiceWaterRegular,
    IncorrectBoilerDataType,
};
use crate::external_conditions::ExternalConditions;
use crate::hem_core::simulation_time::SimulationTimeIteration;
use crate::input::{FuelType, HeatSourceWetDetails, HotWaterSourceDetails};
use indexmap::IndexMap;
use parking_lot::RwLock;
use std::sync::Arc;

/// An object to represent a direct electric boiler
#[derive(Debug)]
pub struct DirectElectricBoiler {
    energy_supply: Arc<RwLock<EnergySupply>>,
    simulation_timestep: f64,
    external_conditions: Arc<ExternalConditions>,
    energy_supply_connections: IndexMap<smartstring::alias::String, EnergySupplyConnection>,
    energy_supply_connection_aux: EnergySupplyConnection,
    service_results: RwLock<Vec<crate::core::heating_systems::boiler::ServiceResult>>,
    boiler_power: f64,
    power_circ_pump: f64,
    power_standby: f64,
    total_time_running_current_timestep: f64,
}
impl DirectElectricBoiler {
    /// Construct a Boiler object
    fn new(
        boiler_data: HeatSourceWetDetails,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_conn_name_auxiliary: &str,
        simulation_timestep: f64,
        external_conditions: Arc<ExternalConditions>,
    ) -> anyhow::Result<Self> {
        let fuel_type = energy_supply.read().fuel_type();
        if !(matches!(fuel_type, FuelType::Electricity)) {
            anyhow::bail!("DirectElectricBoiler requires an electricity energy supply. Got {fuel_type:?} instead.")
        }
        let energy_supply_connection_aux =
            EnergySupply::connection(energy_supply.clone(), energy_supply_conn_name_auxiliary)?;

        match boiler_data {
            HeatSourceWetDetails::DirectElectricBoiler {
                // boiler properties
                rated_power,
                // electricity properties
                electricity_circ_pump,
                electricity_standby,
                ..
            } => Ok(DirectElectricBoiler {
                energy_supply,
                simulation_timestep,
                external_conditions,
                energy_supply_connections: Default::default(),
                energy_supply_connection_aux,
                service_results: Default::default(),
                boiler_power: rated_power,
                power_circ_pump: electricity_circ_pump,
                power_standby: electricity_standby,
                total_time_running_current_timestep: 0.,
            }),
            _ => Err(IncorrectBoilerDataType)?,
        }
    }

    /// Create an EnergySupplyConnection for the service name given
    fn create_service_connection(&mut self, service_name: &str) -> anyhow::Result<()> {
        // Check that service_name is not already registered
        if self.energy_supply_connections.contains_key(service_name) {
            anyhow::bail!("Service name already used: {service_name}");
            // TODO (from Python) Exit just the current case instead of whole program entirely?
        }

        // Set up EnergySupplyConnection for this service
        self.energy_supply_connections.insert(
            service_name.into(),
            EnergySupply::connection(self.energy_supply.clone(), service_name).unwrap(),
        );
        Ok(())
    }

    /// Return a BoilerServiceWater object and create an EnergySupplyConnection for it
    ///
    /// Arguments:
    /// * `boiler_data` - boiler hot water heating properties
    /// * `service_name` - name of the service demanding energy from the boiler
    /// * `temp_hot_water` - temperature of the hot water to be provided, in deg C
    /// * `cold_feed` - reference to ColdWaterSource object
    pub(crate) fn create_service_hot_water_combi(
        boiler: Arc<RwLock<Self>>,
        service_name: &str,
        boiler_data: HotWaterSourceDetails,
        temp_hot_water: f64,
        cold_feed: WaterSupply,
    ) -> Result<BoilerServiceWaterCombi, IncorrectBoilerDataType> {
        boiler
            .write()
            .create_service_connection(service_name)
            .unwrap();
        BoilerServiceWaterCombi::new(
            BoilerForBoilerService::DirectElectricBoiler(boiler.clone()),
            boiler_data,
            service_name.into(),
            temp_hot_water,
            cold_feed,
            boiler.read().simulation_timestep,
        )
    }

    /// Return a BoilerServiceWaterRegular object and create an EnergySupplyConnection for it.
    ///
    /// Arguments:
    /// `service_name` - name of the service demanding energy from the boiler
    /// `controlmin` - reference to a control object which must select current
    ///                the minimum timestep temperature
    /// `controlmax` - reference to a control object which must select current
    ///                the maximum timestep temperature
    /// `control` - reference to a RangeTimeControl object, combining controlmax and controlmin.
    ///             Takes precedence if set.
    fn create_service_hot_water_regular(
        boiler: Arc<RwLock<Self>>,
        service_name: &str,
        control_min: Arc<Control>,
        control_max: Arc<Control>,
        control: Option<Arc<RangeTimeControl>>,
    ) -> anyhow::Result<BoilerServiceWaterRegular> {
        boiler.write().create_service_connection(service_name)?;
        BoilerServiceWaterRegular::new(
            BoilerForBoilerService::DirectElectricBoiler(boiler.clone()),
            service_name.into(),
            control_min,
            control_max,
            control,
        )
    }

    /// Return a BoilerServiceSpace object and create an EnergySupplyConnection for it
    ///
    /// Arguments:
    /// * `service_name` - name of the service demanding energy from the boiler
    /// * `control` - reference to a control object which must implement is_on() and setpnt() funcs
    fn create_service_space_heating(
        boiler: Arc<RwLock<Self>>,
        service_name: &str,
        control: Arc<Control>, // TODO 1.0.0a9 this is a ControlSetPoint in Python
    ) -> anyhow::Result<BoilerServiceSpace> {
        boiler.write().create_service_connection(service_name)?;
        Ok(BoilerServiceSpace::new(
            BoilerForBoilerService::DirectElectricBoiler(boiler.clone()),
            service_name.into(),
            control,
        ))
    }

    /// Return the rated power so that running time reflects actual on-time.
    ///
    /// Direct electric boilers switch on at rated power and off again — they
    /// do not modulate down to spread a small delivery over the entire
    /// available window.
    fn calc_current_boiler_power(&self, _energy_output_provided: f64, time_available: f64) -> f64 {
        if time_available <= 0. {
            return 0.;
        }
        self.boiler_power
    }

    fn calc_energy_output_provided(&self, energy_output_required: f64, time_available: f64) -> f64 {
        let energy_output_max_power = self.boiler_power * time_available;
        let energy_output_provided = min_of_2(energy_output_required, energy_output_max_power);

        energy_output_provided
    }

    /// Calculate time available for the current service
    // Assumes that time spent on other services is evenly spread throughout
    // the timestep so the adjustment for start time below is a proportional
    // reduction of the overall time available, not simply a subtraction
    fn time_available(&self, time_start: f64, time_elapsed_hp: Option<f64>) -> f64 {
        let timestep = self.simulation_timestep;
        let total_time_running_current_timestep = if let Some(time_elapsed_hp) = time_elapsed_hp {
            time_elapsed_hp
        } else {
            self.total_time_running_current_timestep
        };
        let time_available =
            (timestep - total_time_running_current_timestep) * (1. - time_start / timestep);
        time_available
    }

    fn time_running(&self) {
        todo!()
    }

    pub(crate) fn demand_energy(&self) {
        todo!()
    }

    fn electrical_energy_demand(&self) {
        todo!()
    }

    fn calc_auxiliary_energy(&self) {
        todo!()
    }

    pub(crate) fn timestep_end(&self, _simtime: SimulationTimeIteration) -> anyhow::Result<()> {
        todo!()
    }

    fn energy_output_max(&self) {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::controls::time_control::SetpointTimeControl;
    use crate::core::energy_supply::energy_supply::EnergySupplyBuilder;
    use crate::core::units::Orientation360;
    use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
    use crate::hem_core::external_conditions::{DaylightSavingsConfig, ShadingSegment};
    use crate::hem_core::simulation_time::SimulationTime;
    use approx::assert_relative_eq;
    use rstest::{fixture, rstest};
    use serde_json::json;

    #[fixture]
    fn boiler_data() -> HeatSourceWetDetails {
        HeatSourceWetDetails::DirectElectricBoiler {
            energy_supply: "mains_gas".into(),
            electricity_circ_pump: 0.0600,
            electricity_standby: 0.0244,
            rated_power: 24.0,
        }
    }

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 2., 1.)
    }

    #[fixture]
    fn external_conditions(simulation_time: SimulationTime) -> ExternalConditions {
        ExternalConditions::new(
            &simulation_time.iter(),
            vec![0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0],
            vec![3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4],
            vec![200.0, 220.0, 230.0, 240.0, 250.0, 260.0, 260.0, 270.0]
                .into_iter()
                .map(Into::into)
                .collect(),
            vec![333.0, 610.0, 572.0, 420.0, 0.0, 10.0, 90.0, 275.0],
            vec![420.0, 750.0, 425.0, 500.0, 0.0, 40.0, 0.0, 388.0],
            vec![0.2; 8760],
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
                    start360: Orientation360::create_from_180(180.).unwrap(),
                    end360: Orientation360::create_from_180(135.).unwrap(),
                    ..Default::default()
                },
                ShadingSegment {
                    start360: Orientation360::create_from_180(135.).unwrap(),
                    end360: Orientation360::create_from_180(90.).unwrap(),
                    ..Default::default()
                },
                ShadingSegment {
                    start360: Orientation360::create_from_180(90.).unwrap(),
                    end360: Orientation360::create_from_180(45.).unwrap(),
                    ..Default::default()
                },
                ShadingSegment {
                    start360: Orientation360::create_from_180(45.).unwrap(),
                    end360: Orientation360::create_from_180(0.).unwrap(),
                    ..Default::default()
                },
                ShadingSegment {
                    start360: Orientation360::create_from_180(0.).unwrap(),
                    end360: Orientation360::create_from_180(-45.).unwrap(),
                    ..Default::default()
                },
                ShadingSegment {
                    start360: Orientation360::create_from_180(-45.).unwrap(),
                    end360: Orientation360::create_from_180(-90.).unwrap(),
                    ..Default::default()
                },
                ShadingSegment {
                    start360: Orientation360::create_from_180(-90.).unwrap(),
                    end360: Orientation360::create_from_180(-135.).unwrap(),
                    ..Default::default()
                },
                ShadingSegment {
                    start360: Orientation360::create_from_180(-135.).unwrap(),
                    end360: Orientation360::create_from_180(-180.).unwrap(),
                    ..Default::default()
                },
            ]
            .into(),
        )
    }

    #[fixture]
    fn boiler(
        boiler_data: HeatSourceWetDetails,
        simulation_time: SimulationTime,
        external_conditions: ExternalConditions,
    ) -> DirectElectricBoiler {
        DirectElectricBoiler::new(
            boiler_data,
            Arc::new(RwLock::new(
                EnergySupplyBuilder::new(FuelType::Electricity, &simulation_time.iter()).build(),
            )),
            "boiler aux",
            simulation_time.step,
            external_conditions.into(),
        )
        .unwrap()
    }

    /// Test creation of EnergySupplyConnection for the service name given
    #[rstest]
    fn test_create_service_connection(mut boiler: DirectElectricBoiler) {
        let service_name = "new_service";
        // Ensure the service name does not exist in __energy_supply_connections
        assert!(!boiler.energy_supply_connections.contains_key(service_name));

        // Call the method under test
        boiler.create_service_connection(service_name).unwrap();

        // Check that the service name was added to __energy_supply_connections
        assert!(boiler.energy_supply_connections.contains_key(service_name));

        // Check system exit when connection is created with existing service name
        assert!(boiler.create_service_connection(service_name).is_err());
    }

    /// Check BoilerServiceWaterCombi object is created correctly
    #[rstest]
    fn test_create_service_hot_water_combi(
        boiler: DirectElectricBoiler,
        simulation_time: SimulationTime,
    ) {
        let service_name = "service_hot_water_combi";
        let coldfeed =
            ColdWaterSource::new(vec![1.0, 1.2], simulation_time.iter().current_day(), 1.);
        let temp_hot_water = 50.;
        let boiler_data: HotWaterSourceDetails = serde_json::from_value(json!({
            "type": "CombiBoiler",
            "combi_boiler_type": "KeepHot",
            "ColdWaterSource": "mains water",
            "HeatSourceWet": "hp",
            "separate_DHW_tests": "M&L",
            "combi_keep_hot_fuel": "Mixed",
            "keep_hot_test_hours": 24,
            "rejected_energy_1": 0.0004,
            "storage_loss_factor_2": 0.91574,
            "rejected_factor_3": 0,
            "daily_HW_usage": 120,
            "setpoint_temp": 60.0,
        }))
        .unwrap();
        let boiler_service_result = DirectElectricBoiler::create_service_hot_water_combi(
            Arc::new(RwLock::new(boiler)),
            service_name,
            boiler_data,
            temp_hot_water,
            WaterSupply::ColdWaterSource(Arc::new(coldfeed)),
        );
        assert!(boiler_service_result.is_ok());
    }

    /// Check the function returns BoilerServiceWaterRegular object
    #[rstest]
    fn test_create_service_hot_water_regular(
        boiler: DirectElectricBoiler,
        simulation_time: SimulationTime,
    ) {
        let service_name = "service_hot_water_regular";
        let control_min =
            SetpointTimeControl::new(vec![None, None], 0, 1., None, None, simulation_time.step);
        let control_max =
            SetpointTimeControl::new(vec![None, None], 0, 1., None, None, simulation_time.step);

        let boiler_service_result = DirectElectricBoiler::create_service_hot_water_regular(
            Arc::new(RwLock::new(boiler)),
            service_name,
            Arc::new(Control::SetpointTime(control_min)),
            Arc::new(Control::SetpointTime(control_max)),
            None,
        );
        assert!(boiler_service_result.is_ok());
    }

    /// Check the function returns BoilerServiceSpace object
    #[rstest]
    fn test_create_service_space_heating(
        boiler: DirectElectricBoiler,
        simulation_time: SimulationTime,
    ) {
        let service_name = "BoilerServiceSpace";
        let control =
            SetpointTimeControl::new(vec![None, None], 0, 1., None, None, simulation_time.step);

        let boiler_service_result = DirectElectricBoiler::create_service_space_heating(
            Arc::new(RwLock::new(boiler)),
            service_name,
            Arc::new(Control::SetpointTime(control)),
        );
        assert!(boiler_service_result.is_ok());
    }

    #[rstest]
    fn test_calc_current_boiler_power(boiler: DirectElectricBoiler) {
        assert_relative_eq!(boiler.calc_current_boiler_power(10., 0.), 0.);

        // Returns rated power (24 kW) regardless of energy delivered, because
        // electric boilers operate at rated power and switch off — they do not
        // modulate down to spread delivery across the available window.
        assert_relative_eq!(boiler.calc_current_boiler_power(10., 3.), 24.);
    }

    #[rstest]
    fn test_calc_energy_output_provided(boiler: DirectElectricBoiler) {
        assert_relative_eq!(boiler.calc_energy_output_provided(5., 1.), 5.);
        assert_relative_eq!(boiler.calc_energy_output_provided(25., 1.), 24.);
    }

    #[rstest]
    fn test_time_available(boiler: DirectElectricBoiler, simulation_time: SimulationTime) {
        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_relative_eq!(boiler.time_available(0., None), &[1., 1.][t_idx]);
        }

        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_relative_eq!(boiler.time_available(0.2, Some(0.5)), &[0.4, 0.4][t_idx]);
        }
    }
}
