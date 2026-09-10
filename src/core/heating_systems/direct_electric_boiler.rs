use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::core::heating_systems::boiler::IncorrectBoilerDataType;
use crate::external_conditions::ExternalConditions;
use crate::input::{FuelType, HeatSourceWetDetails};
use indexmap::IndexMap;
use parking_lot::RwLock;
use std::sync::Arc;

/// An object to represent a direct electric boiler
struct DirectElectricBoiler {
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
}
