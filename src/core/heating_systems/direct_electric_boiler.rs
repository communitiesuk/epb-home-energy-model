use crate::compare_floats::min_of_2;
use crate::core::common::WaterSupply;
use crate::core::controls::time_control::{
    Control, RangeTimeControl, SetpointOrCombinationControl,
};
use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::core::heating_systems::boiler::{
    BoilerForBoilerService, BoilerServiceSpace, BoilerServiceWaterCombi, BoilerServiceWaterRegular,
    CombiBoilerConfig, IncorrectBoilerDataType, KeepHotCombiBoilerConfig, ServiceResult,
    ServiceType,
};
use crate::external_conditions::ExternalConditions;
use crate::hem_core::simulation_time::SimulationTimeIteration;
use crate::input::{
    CombiBoilerType, CombiKeepHotFuel, FuelType, HeatSourceWetDetails, HotWaterSourceDetails,
};
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
            EnergySupply::connection(self.energy_supply.clone(), service_name)?,
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
        // TODO: look at improving error handling here and in boiler.rs
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
        control_min: Option<SetpointOrCombinationControl>,
        control_max: Option<SetpointOrCombinationControl>,
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
        control: Control, // TODO 1.0.0a9 this is a ControlSetPoint in Python
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
        min_of_2(energy_output_required, energy_output_max_power)
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

        (timestep - total_time_running_current_timestep) * (1. - time_start / timestep)
    }

    fn time_running(&self, energy_output_provided: f64, time_available: f64) -> f64 {
        // Calculate running time of Boiler
        let current_boiler_power =
            self.calc_current_boiler_power(energy_output_provided, time_available);

        if current_boiler_power <= 0.0 {
            0.0
        } else {
            min_of_2(
                energy_output_provided / current_boiler_power,
                time_available,
            )
        }
    }

    /// Calculate energy required by boiler to satisfy demand for the service indicated.
    pub(crate) fn demand_energy(
        &mut self,
        service_name: &str,
        service_type: ServiceType,
        energy_output_required: f64,
        temp_flow: f64,
        temp_return_feed: Option<f64>,
        hybrid_service: Option<bool>,
        time_start: Option<f64>,
        time_elapsed_hp: Option<f64>,
        update_heat_source_state: Option<bool>,
        combi_boiler_config: Option<CombiBoilerConfig>,
    ) -> anyhow::Result<(f64, Option<f64>)> {
        // Account for time control where present. If no control present, assume
        // system is always active (except for basic thermostatic control, which
        // is implicit in demand calculation).
        // if self.__control is None or self.__control.is_on():
        //     // Energy that heater is able to supply is limited by power rating
        //     energy_output_provided = min(energy_output_required, self.__boiler_power * self.__simulation_time.timestep())
        // else:
        //

        let time_start = time_start.unwrap_or(0.0);
        let hybrid_service_bool = hybrid_service.unwrap_or(false);
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);
        let combi_boiler_config = combi_boiler_config.unwrap_or(CombiBoilerConfig {
            combi_loss: 0.,
            combi_type: Default::default(),
            keep_hot_config: None,
        });

        let time_available = self.time_available(time_start, time_elapsed_hp);
        let energy_output_provided =
            self.calc_energy_output_provided(energy_output_required, time_available);

        // TODO (from Python) Ideally, the boiler power used for the running time calculation
        //      would account for space heating demand for all zones, but the
        //      calculation flow does not allow for this without circularity.
        //      Therefore, the value for time running returned from this function
        //      (used in the hybrid HP calculation) will be slightly inaccurate.

        let time_running_current_service =
            self.time_running(energy_output_provided, time_available);

        if update_heat_source_state {
            self.total_time_running_current_timestep += time_running_current_service;

            let combi_boiler_config = match service_type {
                ServiceType::WaterCombi => {
                    let keep_hot_config = match combi_boiler_config.combi_type {
                        CombiBoilerType::KeepHot => {
                            Some(combi_boiler_config.keep_hot_config.unwrap_or(
                                KeepHotCombiBoilerConfig {
                                    keep_hot_on: true,
                                    keep_hot_fuel: CombiKeepHotFuel::MainBoilerFuel,
                                },
                            ))
                        }
                        _ => None,
                    };

                    Some(CombiBoilerConfig {
                        combi_loss: combi_boiler_config.combi_loss,
                        combi_type: combi_boiler_config.combi_type,
                        keep_hot_config,
                    })
                }
                _ => None,
            };

            // Save results that are needed later (in the timestep_end function)
            let service_result = ServiceResult {
                service_name: service_name.into(),
                service_type,
                temp_flow,
                temp_return_feed,
                energy_output_required,
                energy_output_provided,
                time_available,
                _time_start: time_start,
                _time_elapsed_hp: time_elapsed_hp,
                combi_boiler_config,
            };

            self.service_results.write().push(service_result);
        }

        Ok(if hybrid_service_bool {
            (energy_output_provided, Some(time_running_current_service))
        } else {
            (energy_output_provided, None)
        })
    }

    /// Calculate boiler electrical demand for all services (excl. auxiliary),
    /// and request this from relevant EnergySupplyConnection
    fn electrical_energy_demand(&self, timestep_idx: usize) -> anyhow::Result<()> {
        for service_data in self.service_results.read().iter() {
            let ServiceResult {
                service_name,
                temp_return_feed,
                energy_output_provided,
                ..
            } = service_data;
            let demand = if temp_return_feed.is_some() {
                *energy_output_provided
            } else {
                0.
            };

            self
                .energy_supply_connections
                .get(service_name)
                .ok_or_else(|| anyhow::anyhow!("Expected direct electric boiler energy_supply_connection with service name {service_name}"))?
                .demand_energy(demand, timestep_idx)?;
        }
        Ok(())
    }

    /// Calculation of boiler electrical consumption
    fn calc_auxiliary_energy(
        &self,
        timestep_idx: usize,
        time_remaining_current_timestep: f64,
    ) -> anyhow::Result<()> {
        // Energy used by circulation pump
        let mut energy_aux = self.total_time_running_current_timestep * self.power_circ_pump;

        // Energy used in standby mode
        energy_aux += self.power_standby * time_remaining_current_timestep;
        self.energy_supply_connection_aux
            .demand_energy(energy_aux, timestep_idx)
    }

    /// Calculations to be done at the end of each timestep
    pub(crate) fn timestep_end(&mut self, simtime: SimulationTimeIteration) -> anyhow::Result<()> {
        self.electrical_energy_demand(simtime.index)?;

        let time_remaining_current_step =
            simtime.timestep - self.total_time_running_current_timestep;

        self.calc_auxiliary_energy(simtime.index, time_remaining_current_step)?;

        // Variables below need to be reset at the end of each timestep
        self.total_time_running_current_timestep = Default::default();
        self.service_results = Default::default();

        Ok(())
    }

    pub(crate) fn energy_output_max(
        &self,
        time_start: Option<f64>,
        time_elapsed_hp: Option<f64>,
    ) -> f64 {
        let time_start = time_start.unwrap_or(0.);
        let time_available = self.time_available(time_start, time_elapsed_hp);

        self.boiler_power * time_available
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
                EnergySupplyBuilder::new(
                    FuelType::Electricity,
                    simulation_time.iter().total_steps(),
                )
                .build(),
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
            Some(SetpointOrCombinationControl::SetpointTime(control_min.into())),
            Some(SetpointOrCombinationControl::SetpointTime(control_max.into())),
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
            Control::SetpointTime(control.into()),
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

    #[rstest]
    fn test_demand_energy(mut boiler: DirectElectricBoiler, simulation_time: SimulationTime) {
        // Test with different values of  service types and hybrid_service_bool
        boiler
            .create_service_connection("boiler_demand_energy")
            .unwrap();

        // Both timesteps deliver the full 10 kWh because the boiler's running
        // time (10 / 24 h) is well within the 1 h timestep, leaving capacity
        // for the second call.
        for (t_idx, _) in simulation_time.iter().enumerate() {
            let result = boiler
                .demand_energy(
                    "boiler_demand_energy",
                    ServiceType::WaterCombi,
                    10.,
                    45.,
                    Some(37.),
                    Some(false),
                    None,
                    None,
                    None,
                    Some(CombiBoilerConfig {
                        combi_loss: 1.2,
                        combi_type: CombiBoilerType::KeepHot,
                        keep_hot_config: None,
                    }),
                )
                .unwrap();

            assert_eq!(result, [(10., None), (10., None)][t_idx]);
        }

        boiler
            .create_service_connection("boiler_demand_energy_with_hybrid")
            .unwrap();

        for (t_idx, _) in simulation_time.iter().enumerate() {
            let result = boiler
                .demand_energy(
                    "boiler_demand_energy_with_hybrid",
                    ServiceType::Space,
                    100.,
                    45.,
                    Some(37.),
                    Some(true),
                    None,
                    Some(0.),
                    None,
                    None,
                )
                .unwrap();

            assert_eq!(result, [(24., Some(1.)), (24., Some(1.))][t_idx]);
        }

        // Test with time_elapsed_hp
        boiler
            .create_service_connection("boiler_demand_energy_hybrid_time_elapsed")
            .unwrap();

        for (t_idx, _) in simulation_time.iter().enumerate() {
            let result = boiler
                .demand_energy(
                    "boiler_demand_energy_hybrid_time_elapsed",
                    ServiceType::Space,
                    100.,
                    45.,
                    Some(37.),
                    Some(true),
                    None,
                    Some(0.5),
                    None,
                    None,
                )
                .unwrap();

            assert_eq!(result, [(12., Some(0.5)), (12., Some(0.5))][t_idx]);
        }

        let sercive_results = &boiler.service_results.read();
        let service_names = sercive_results
            .iter()
            .map(|service_result| service_result.service_name.as_str())
            .collect::<Vec<&str>>();

        assert!(service_names.contains(&"boiler_demand_energy"));
        assert!(service_names.contains(&"boiler_demand_energy_with_hybrid"));
        assert!(service_names.contains(&"boiler_demand_energy_hybrid_time_elapsed"));
    }

    #[rstest]
    /// Two services in the same timestep share the boiler's capacity.
    ///
    /// The 24 kW boiler has 1 h available. The first service draws 10 kWh
    /// (running time 10/24 h), leaving 14/24 h for the second. Requesting
    /// 20 kWh on the second call is therefore capped at 14 kWh.
    fn test_demand_energy_two_calls_same_timestep(mut boiler: DirectElectricBoiler) {
        boiler.create_service_connection("service_a").unwrap();

        boiler.create_service_connection("service_b").unwrap();

        // First call: 10 kWh easily within the 24 kWh single-timestep capacity
        let result_a = boiler
            .demand_energy(
                "service_a",
                ServiceType::Space,
                10.,
                45.,
                Some(37.),
                Some(false),
                None,
                None,
                None,
                None,
            )
            .unwrap();

        // 10 kWh < 24 kWh capacity, so delivered in full
        assert_eq!(result_a, (10., None));

        // Second call in the same timestep (no timestep_end between):
        // remaining time = 1 - 10/24 = 14/24 h, max output = 24 × 14/24 = 14 kWh
        let result_b = boiler
            .demand_energy(
                "service_b",
                ServiceType::Space,
                20.,
                45.,
                Some(37.),
                Some(false),
                None,
                None,
                None,
                None,
            )
            .unwrap();

        assert_relative_eq!(result_b.0, 14.);
        assert_eq!(result_b.1, None);
    }

    #[ignore = "usage of mocks in Python, won't replicate for now"]
    #[rstest]
    fn test_fuel_demand() {
        todo!()
    }

    #[ignore = "usage of mocks in Python, won't replicate for now"]
    #[rstest]
    fn test_fuel_demand_with_no_return_feed() {
        todo!()
    }

    #[rstest]
    /// Check boiler electrical consumption
    fn test_calc_auxiliary_energy(boiler: DirectElectricBoiler) {
        // Check the function runs without throwing errors
        assert!(boiler.calc_auxiliary_energy(1, 0.).is_ok())

        // Mocks used for second assertion, won't replicate for now
    }

    #[ignore = "usage of mocks in Python, won't replicate for now"]
    #[rstest]
    fn test_calc_auxiliary_energy_with_space_heating() {
        todo!()
    }

    #[rstest]
    fn test_timestep_end(mut boiler: DirectElectricBoiler, simulation_time: SimulationTime) {
        boiler
            .create_service_connection("boiler_demand_energy")
            .unwrap();

        boiler
            .demand_energy(
                "boiler_demand_energy",
                ServiceType::WaterCombi,
                10.,
                45.,
                Some(60.),
                Some(false),
                None,
                None,
                None,
                None,
            )
            .unwrap();

        // Running time = energy_delivered / rated_power = 10 / 24
        assert_relative_eq!(
            boiler.total_time_running_current_timestep,
            0.4166666666666667
        );

        assert_eq!(
            boiler.service_results.read()[0].service_name.as_str(),
            "boiler_demand_energy"
        );

        // Call the method under test
        boiler
            .timestep_end(simulation_time.iter().current_iteration())
            .unwrap();

        // Assertions to check if the internal state was updated correctly
        assert_eq!(boiler.total_time_running_current_timestep, 0.);
        assert!(boiler.service_results.read().is_empty());
    }

    #[rstest]
    fn test_energy_output_max(boiler: DirectElectricBoiler) {
        assert_relative_eq!(boiler.energy_output_max(Some(0.), None), 24.);
        assert_relative_eq!(boiler.energy_output_max(Some(0.5), None), 12.);
    }

    #[rstest]
    /// DirectElectricBoiler must reject a non-electricity energy supply.
    fn test_init_raises_for_non_electricity_supply(
        boiler_data: HeatSourceWetDetails,
        simulation_time: SimulationTime,
        external_conditions: ExternalConditions,
    ) {
        let result = DirectElectricBoiler::new(
            boiler_data,
            Arc::new(RwLock::new(
                EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time.iter().total_steps())
                    .build(),
            )),
            "boiler aux",
            simulation_time.step,
            external_conditions.into(),
        );

        let error = result.err().unwrap();
        assert!(error
            .to_string()
            .contains("DirectElectricBoiler requires an electricity energy supply."));
    }
}
