use crate::core::common::WaterSupply;
use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::core::heating_systems::boiler::IncorrectBoilerDataType;
use crate::external_conditions::ExternalConditions;
use crate::hem_core::simulation_time::SimulationTimeIteration;
use crate::input::{FuelType, HeatSourceWetDetails};
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
    fn create_service_hot_water_combi(
        &mut self,
        service_name: &str,
        _boiler_data: HeatSourceWetDetails,
        _temp_hot_water: f64,
        _cold_feed: WaterSupply,
    ) -> Result<(), anyhow::Error> {
        self.create_service_connection(service_name)?;
        todo!()
        // BoilerServiceWaterCombi::new(
        //     self,
        //     boiler_data,
        //     service_name.parse()?,
        //     temp_hot_water,
        //     cold_feed,
        //     self.simulation_timestep,
        // )?;
    }

    fn create_service_hot_water_regular(&self) {
        todo!()
    }

    fn create_service_space_heating(&self) {
        todo!()
    }

    fn calc_current_boiler_power(&self) {
        todo!()
    }

    fn calc_energy_output_provided(&self) {
        todo!()
    }

    fn time_available(&self) {
        todo!()
    }

    fn time_running(&self) {
        todo!()
    }

    fn demand_energy(&self) {
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
    use crate::core::energy_supply::energy_supply::EnergySupplyBuilder;
    use crate::core::units::Orientation360;
    use crate::hem_core::external_conditions::{DaylightSavingsConfig, ShadingSegment};
    use crate::hem_core::simulation_time::SimulationTime;
    use rstest::{fixture, rstest};

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
}
