use crate::compare_floats::{max_of_2, min_of_2};
use crate::core::common::WaterSourceWithTemperature;
use crate::core::controls::time_control::{per_control, Control, ControlBehaviour};
use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::core::material_properties::WATER;
use crate::core::units::{HOURS_PER_DAY, WATTS_PER_KILOWATT};
use crate::core::water_heat_demand::dhw_demand::{DemandVolTargetKey, VolumeReference};
use crate::simulation_time::SimulationTimeIteration;
use anyhow::bail;
use indexmap::IndexMap;
use parking_lot::{Mutex, RwLock};
use std::collections::HashMap;
use std::sync::Arc;

#[derive(Debug)]
pub struct HeatNetworkServiceWaterDirect {
    heat_network: Arc<Mutex<HeatNetwork>>,
    service_name: String,
    temperature_hot_water: f64, // in C
    cold_feed: WaterSourceWithTemperature,
}

/// An object to represent a water heating service provided by a heat network.
///
/// This object contains the parts of the heat network calculation that are
/// specific to providing hot water directly to the dwelling.
impl HeatNetworkServiceWaterDirect {
    /// Arguments:
    /// * `heat_network` - reference to the HeatNetwork object providing the service
    /// * `service_name` - name of the service demanding energy from the heat network
    /// * `temp_hot_water` - temperature of the hot water to be provided, in deg C
    /// * `cold_feed` - reference to ColdWaterSource object
    pub(crate) fn new(
        heat_network: Arc<Mutex<HeatNetwork>>,
        service_name: String,
        temperature_hot_water: f64,
        cold_feed: WaterSourceWithTemperature,
    ) -> Self {
        Self {
            heat_network,
            service_name,
            temperature_hot_water,
            cold_feed,
        }
    }

    pub fn get_cold_water_source(&self) -> &WaterSourceWithTemperature {
        &self.cold_feed
    }

    pub fn temp_hot_water(&self) -> f64 {
        self.temperature_hot_water
    }

    /// Demand energy for hot water (in kWh) from the heat network
    pub fn demand_hot_water(
        &self,
        volume_demanded_target: IndexMap<DemandVolTargetKey, VolumeReference>,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        // Calculate energy needed to meet hot water demand
        let volume_demanded = volume_demanded_target
            .get(&DemandVolTargetKey::TempHotWater)
            .map(|volume_reference| volume_reference.warm_vol)
            .unwrap_or(0.0);
        let energy_content_kwh_per_litre = WATER.volumetric_energy_content_kwh_per_litre(
            self.temperature_hot_water,
            self.cold_feed.temperature(simtime, None),
        );
        let energy_demand = volume_demanded * energy_content_kwh_per_litre;

        self.heat_network.lock().demand_energy(
            &self.service_name,
            energy_demand,
            None,
            None,
            simtime.index,
        )
    }
}

/// An object to represent a water heating service provided by a heat network.
///
/// This object contains the parts of the heat network calculation that are
/// specific to providing hot water to the dwelling via a hot water cylinder.
#[derive(Clone, Debug)]
pub struct HeatNetworkServiceWaterStorage {
    heat_network: Arc<Mutex<HeatNetwork>>,
    service_name: String,
    control: Arc<Control>,
    control_min: Arc<Control>,
    control_max: Arc<Control>,
}

impl HeatNetworkServiceWaterStorage {
    /// Arguments:
    /// * `heat_network` - reference to the HeatNetwork object providing the service
    /// * `service_name` - name of the service demanding energy from the heat network
    /// * `control_min` - reference to a control object which must select current the minimum timestep temperature
    /// * `control_max` - reference to a control object which must select current the maximum timestep temperature
    pub(crate) fn new(
        heat_network: Arc<Mutex<HeatNetwork>>,
        service_name: String,
        control_min: Arc<Control>,
        control_max: Arc<Control>,
    ) -> Self {
        let control = control_min.clone();

        Self {
            heat_network,
            service_name,
            control,
            control_min,
            control_max,
        }
    }

    pub(crate) fn _temp_setpnt(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> (Option<f64>, Option<f64>) {
        (
            self.control_min.setpnt(simulation_time_iteration),
            self.control_max.setpnt(simulation_time_iteration),
        )
    }

    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        _temp_return: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        if !self.is_on(simulation_time_iteration) {
            return 0.;
        }

        self.heat_network.lock().demand_energy(
            &self.service_name,
            energy_demand,
            None,
            None,
            simulation_time_iteration.index,
        )
    }

    pub fn energy_output_max(
        &self,
        _temp_return: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        if !self.is_on(simulation_time_iteration) {
            return 0.;
        }

        let control_max_setpnt = self.control_max.setpnt(simulation_time_iteration);

        self.heat_network
            .lock()
            .energy_output_max(control_max_setpnt, None)
    }

    fn is_on(&self, simulation_time_iteration: &SimulationTimeIteration) -> bool {
        per_control!(self.control.as_ref(), ctrl => { ctrl.is_on(simulation_time_iteration) })
    }
}

#[derive(Debug)]
pub struct HeatNetworkServiceSpace {
    heat_network: Arc<Mutex<HeatNetwork>>,
    service_name: String,
    control: Arc<Control>,
}

impl HeatNetworkServiceSpace {
    pub(crate) fn new(
        heat_network: Arc<Mutex<HeatNetwork>>,
        service_name: String,
        control: Arc<Control>,
    ) -> Self {
        Self {
            heat_network,
            service_name,
            control,
        }
    }

    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        _temp_flow: f64,
        _temp_return: f64,
        time_start: Option<f64>,
        update_heat_source_state: Option<bool>,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        let time_start = time_start.unwrap_or(0.);
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);

        if !self.is_on(simulation_time_iteration) {
            return 0.;
        }

        self.heat_network.lock().demand_energy(
            &self.service_name,
            energy_demand,
            Some(time_start),
            Some(update_heat_source_state),
            simulation_time_iteration.index,
        )
    }

    pub fn energy_output_max(
        &self,
        temp_output: f64,
        _temp_return_feed: f64,
        time_start: Option<f64>,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        let time_start = time_start.unwrap_or(0.);

        if !self.is_on(simulation_time_iteration) {
            return 0.;
        }

        self.heat_network
            .lock()
            .energy_output_max(Some(temp_output), Some(time_start))
    }

    pub fn temperature_setpnt(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<f64> {
        per_control!(&self.control.as_ref(), ctrl => { ctrl.setpnt(simulation_time_iteration) })
    }

    pub fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool> {
        per_control!(&self.control.as_ref(), ctrl => { <_ as ControlBehaviour>::in_required_period(ctrl, simulation_time_iteration) })
    }

    fn is_on(&self, simulation_time_iteration: &SimulationTimeIteration) -> bool {
        per_control!(&self.control.as_ref(), ctrl => { ctrl.is_on(simulation_time_iteration) })
    }
}

#[derive(Clone, Debug)]
pub(crate) struct HeatNetwork {
    power_max_in_kw: f64,
    daily_loss: f64,                         // in kWh/day
    building_level_distribution_losses: f64, // in watts
    energy_supply: Arc<RwLock<EnergySupply>>,
    energy_supply_connections: HashMap<String, EnergySupplyConnection>,
    energy_supply_connection_aux: EnergySupplyConnection,
    energy_supply_connection_building_level_distribution_losses: EnergySupplyConnection,
    total_time_running_current_timestep: f64,
    simulation_timestep: f64,
}

impl HeatNetwork {
    pub(crate) fn new(
        power_max_in_kw: f64,
        daily_loss: f64,
        building_level_distribution_losses: f64,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_conn_name_auxiliary: String,
        energy_supply_conn_name_building_level_distribution_losses: String,
        simulation_timestep: f64,
    ) -> Self {
        Self {
            power_max_in_kw,
            daily_loss,
            building_level_distribution_losses,
            energy_supply: energy_supply.clone(),
            energy_supply_connections: Default::default(),
            energy_supply_connection_aux: EnergySupply::connection(
                energy_supply.clone(),
                energy_supply_conn_name_auxiliary.as_str(),
            )
            .unwrap(),
            energy_supply_connection_building_level_distribution_losses: EnergySupply::connection(
                energy_supply,
                energy_supply_conn_name_building_level_distribution_losses.as_str(),
            )
            .unwrap(),
            total_time_running_current_timestep: Default::default(),
            simulation_timestep,
        }
    }

    /// Create an EnergySupplyConnection for the service name given
    pub fn create_service_connection(
        heat_network: Arc<Mutex<Self>>,
        service_name: &str,
    ) -> anyhow::Result<()> {
        if heat_network
            .lock()
            .energy_supply_connections
            .contains_key(service_name)
        {
            bail!("Error: Service name already used: {service_name}");
        }
        let energy_supply = heat_network.lock().energy_supply.clone();

        // Set up EnergySupplyConnection for this service
        heat_network.lock().energy_supply_connections.insert(
            service_name.to_string(),
            EnergySupply::connection(energy_supply, service_name).unwrap(),
        );

        Ok(())
    }

    pub fn create_service_hot_water_direct(
        heat_network: Arc<Mutex<Self>>,
        service_name: String,
        temperature_hot_water: f64,
        cold_feed: WaterSourceWithTemperature,
    ) -> HeatNetworkServiceWaterDirect {
        Self::create_service_connection(heat_network.clone(), service_name.as_str()).unwrap();

        HeatNetworkServiceWaterDirect::new(
            heat_network,
            service_name,
            temperature_hot_water,
            cold_feed,
        )
    }

    pub(crate) fn create_service_hot_water_storage(
        heat_network: Arc<Mutex<Self>>,
        service_name: String,
        control_min: Arc<Control>,
        control_max: Arc<Control>,
    ) -> HeatNetworkServiceWaterStorage {
        Self::create_service_connection(heat_network.clone(), service_name.as_str()).unwrap();

        HeatNetworkServiceWaterStorage::new(heat_network, service_name, control_min, control_max)
    }

    pub(crate) fn create_service_space_heating(
        heat_network: Arc<Mutex<Self>>,
        service_name: String,
        control: Arc<Control>,
    ) -> HeatNetworkServiceSpace {
        Self::create_service_connection(heat_network.clone(), service_name.as_str()).unwrap();

        HeatNetworkServiceSpace::new(heat_network, service_name, control)
    }

    /// Calculate the maximum energy output of the heat network, accounting
    /// for time spent on higher-priority services.
    fn energy_output_max(&self, _temp_output: Option<f64>, time_start: Option<f64>) -> f64 {
        let time_start = time_start.unwrap_or(0.);

        let time_available = self.time_available(time_start, self.simulation_timestep);

        self.power_max_in_kw * time_available
    }

    /// Calculate time available for the current service
    fn time_available(&self, time_start: f64, timestep: f64) -> f64 {
        // Assumes that time spent on other services is evenly spread throughout
        // the timestep so the adjustment for start time below is a proportional
        // reduction of the overall time available, not simply a subtraction
        (timestep - self.total_time_running_current_timestep) * (1. - time_start / timestep)
    }

    /// Calculate energy required by heat network to satisfy demand for the service indicated.
    pub fn demand_energy(
        &mut self,
        service_name: &str,
        energy_output_required: f64,
        time_start: Option<f64>,
        update_heat_source_state: Option<bool>,
        timestep_idx: usize,
    ) -> f64 {
        let time_start = time_start.unwrap_or(0.);
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);
        let energy_output_max = self.energy_output_max(None, None);
        if energy_output_max == 0. {
            return energy_output_max;
        }
        let energy_output_provided =
            max_of_2(0., min_of_2(energy_output_required, energy_output_max));

        if update_heat_source_state {
            self.energy_supply_connections[service_name]
                .demand_energy(energy_output_provided, timestep_idx)
                .unwrap();
        }

        let time_available = self.time_available(time_start, self.simulation_timestep);

        if update_heat_source_state {
            self.total_time_running_current_timestep +=
                (energy_output_provided / energy_output_max) * time_available;
        }

        energy_output_provided
    }

    /// Calculations to be done at the end of each timestep
    pub fn timestep_end(&mut self, timestep_idx: usize) {
        // Energy required to overcome losses
        self.energy_supply_connection_aux
            .demand_energy(self.hiu_loss(), timestep_idx)
            .unwrap();
        self.energy_supply_connection_building_level_distribution_losses
            .demand_energy(self.building_level_loss(), timestep_idx)
            .unwrap();

        // Variables below need to be reset at the end of each timestep
        self.total_time_running_current_timestep = Default::default();
    }

    /// Standing heat loss from the HIU (heat interface unit) in kWh
    pub fn hiu_loss(&self) -> f64 {
        // daily_loss to be sourced from the PCDB, in kWh/day
        self.daily_loss / HOURS_PER_DAY as f64 * self.simulation_timestep
    }

    /// Converts building level distribution loss from watts to kWh
    pub fn building_level_loss(&self) -> f64 {
        self.building_level_distribution_losses / WATTS_PER_KILOWATT as f64
            * self.simulation_timestep
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::controls::time_control::SetpointTimeControl;
    use crate::core::energy_supply::energy_supply::EnergySupplyBuilder;
    use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
    use crate::input::FuelType;
    use crate::simulation_time::SimulationTime;
    use approx::assert_relative_eq;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    pub fn two_len_simulation_time() -> SimulationTime {
        SimulationTime::new(0., 2., 1.)
    }

    const SERVICE_NAME: &str = "heat_network_service";

    #[fixture]
    fn dummy_heat_network(two_len_simulation_time: SimulationTime) -> Arc<Mutex<HeatNetwork>> {
        Arc::new(Mutex::new(HeatNetwork::new(
            0.0,
            0.0,
            0.0,
            Arc::new(RwLock::new(
                EnergySupplyBuilder::new(
                    FuelType::Electricity,
                    two_len_simulation_time.total_steps(),
                )
                .build(),
            )),
            "aux".to_string(),
            "distro_losses".to_string(),
            1.0,
        )))
    }

    #[rstest]
    fn test_service_is_on(
        two_len_simulation_time: SimulationTime,
        dummy_heat_network: Arc<Mutex<HeatNetwork>>,
    ) {
        let control = SetpointTimeControl::new(
            vec![Some(21.0), Some(21.0), None],
            0,
            1.0,
            None,
            None,
            None,
            Default::default(),
            two_len_simulation_time.step,
        )
        .unwrap();

        // there is no base call in Rust so use one of the concrete implementations
        let heat_network_service = HeatNetworkServiceWaterStorage::new(
            dummy_heat_network.clone(),
            SERVICE_NAME.to_owned(),
            Arc::new(Control::SetpointTime(control.clone())),
            Arc::new(Control::SetpointTime(control.clone())),
        );
        assert!(heat_network_service.is_on(&two_len_simulation_time.iter().next().unwrap()));

        let heat_network_service_no_control = HeatNetworkServiceWaterStorage::new(
            dummy_heat_network,
            SERVICE_NAME.to_owned(),
            Arc::new(Control::SetpointTime(control.clone())),
            Arc::new(Control::SetpointTime(control.clone())),
        );
        assert!(
            heat_network_service_no_control.is_on(&two_len_simulation_time.iter().next().unwrap())
        );
    }

    #[fixture]
    fn heat_network_for_water_direct(
        two_len_simulation_time: SimulationTime,
    ) -> Arc<Mutex<HeatNetwork>> {
        let energy_supply =
            EnergySupplyBuilder::new(FuelType::Custom, two_len_simulation_time.total_steps())
                .build();
        let energy_supply_conn_name_auxiliary = "heat_network_auxiliary";
        let energy_supply_conn_name_building_level_distribution_losses =
            "HeatNetwork_building_level_distribution_losses";
        Arc::new(Mutex::new(HeatNetwork::new(
            18.0,
            1.0,
            0.8,
            Arc::new(RwLock::new(energy_supply)),
            energy_supply_conn_name_auxiliary.to_owned(),
            energy_supply_conn_name_building_level_distribution_losses.to_owned(),
            two_len_simulation_time.step,
        )))
    }

    #[fixture]
    fn heat_network_service_water_direct(
        two_len_simulation_time: SimulationTime,
        heat_network_for_water_direct: Arc<Mutex<HeatNetwork>>,
    ) -> HeatNetworkServiceWaterDirect {
        //

        let heat_network = heat_network_for_water_direct;

        HeatNetwork::create_service_connection(heat_network.clone(), "heat_network_test").unwrap();

        let cold_water_temps = vec![1.0, 1.2];
        let cold_feed = ColdWaterSource::new(cold_water_temps, 0, two_len_simulation_time.step);
        let return_temp = 60.;

        HeatNetworkServiceWaterDirect::new(
            heat_network,
            "heat_network_test".to_owned(),
            return_temp,
            WaterSourceWithTemperature::ColdWaterSource(Arc::new(cold_feed)),
        )
    }

    #[rstest]
    fn test_heat_network_service_water(
        two_len_simulation_time: SimulationTime,
        heat_network_service_water_direct: HeatNetworkServiceWaterDirect,
        heat_network_for_water_direct: Arc<Mutex<HeatNetwork>>,
    ) {
        let volume_demanded = [
            IndexMap::from([
                (
                    41.0.into(),
                    VolumeReference {
                        warm_temp: 41.0,
                        warm_vol: 48.0,
                    },
                ),
                (
                    43.0.into(),
                    VolumeReference {
                        warm_temp: 43.0,
                        warm_vol: 100.0,
                    },
                ),
                (
                    40.0.into(),
                    VolumeReference {
                        warm_temp: 40.0,
                        warm_vol: 0.0,
                    },
                ),
                (
                    DemandVolTargetKey::TempHotWater,
                    VolumeReference {
                        warm_temp: 55.0,
                        warm_vol: 110.59194954841298,
                    },
                ),
            ]),
            IndexMap::from([
                (
                    41.0.into(),
                    VolumeReference {
                        warm_temp: 41.0,
                        warm_vol: 48.0,
                    },
                ),
                (
                    DemandVolTargetKey::TempHotWater,
                    VolumeReference {
                        warm_temp: 55.0,
                        warm_vol: 32.60190808710678,
                    },
                ),
            ]),
        ];

        for (t_idx, t_it) in two_len_simulation_time.iter().enumerate() {
            assert_relative_eq!(
                heat_network_service_water_direct
                    .demand_hot_water(volume_demanded[t_idx].clone(), t_it),
                [7.5834, 2.2279][t_idx],
                max_relative = 1e-4
            );
            heat_network_for_water_direct.lock().timestep_end(t_idx);
        }
    }

    #[fixture]
    #[once]
    fn heat_network_for_water_storage(
        two_len_simulation_time: SimulationTime,
    ) -> Arc<Mutex<HeatNetwork>> {
        let energy_supply =
            EnergySupplyBuilder::new(FuelType::Custom, two_len_simulation_time.total_steps())
                .build();
        let energy_supply_conn_name_auxiliary = "heat_network_auxiliary";
        let energy_supply_conn_name_building_level_distribution_losses =
            "HeatNetwork_building_level_distribution_losses";

        Arc::new(Mutex::new(HeatNetwork::new(
            7.0,
            1.0,
            0.8,
            Arc::new(RwLock::new(energy_supply)),
            energy_supply_conn_name_auxiliary.to_owned(),
            energy_supply_conn_name_building_level_distribution_losses.to_owned(),
            two_len_simulation_time.step,
        )))
    }

    #[fixture]
    fn heat_network_service_water_storage(
        heat_network_for_water_storage: &Arc<Mutex<HeatNetwork>>,
        two_len_simulation_time: SimulationTime,
    ) -> HeatNetworkServiceWaterStorage {
        let heat_network = heat_network_for_water_storage;

        let _ = HeatNetwork::create_service_connection(heat_network.clone(), "heat_network_test");

        let control_min = SetpointTimeControl::new(
            vec![Some(52.), Some(52.), None],
            0,
            1.0,
            None,
            None,
            None,
            Default::default(),
            two_len_simulation_time.step,
        )
        .unwrap();

        let control_max = SetpointTimeControl::new(
            vec![Some(60.), Some(60.), None],
            0,
            1.0,
            None,
            None,
            None,
            Default::default(),
            two_len_simulation_time.step,
        )
        .unwrap();

        HeatNetworkServiceWaterStorage::new(
            heat_network.clone(),
            "heat_network_test".to_owned(),
            Arc::new(Control::SetpointTime(control_min)),
            Arc::new(Control::SetpointTime(control_max)),
        )
    }

    #[rstest]
    fn test_heat_network_service_water_storage(
        two_len_simulation_time: SimulationTime,
        mut heat_network_service_water_storage: HeatNetworkServiceWaterStorage,
        heat_network_for_water_storage: &Arc<Mutex<HeatNetwork>>,
    ) {
        let heat_network = heat_network_for_water_storage.clone();
        let energy_demanded = [10.0, 2.0];
        for (t_idx, t_it) in two_len_simulation_time.iter().enumerate() {
            assert_eq!(
                heat_network_service_water_storage.demand_energy(
                    energy_demanded[t_idx],
                    60.,
                    &t_it
                ),
                [7.0, 2.0][t_idx]
            );
            heat_network.lock().timestep_end(t_idx);
        }
    }

    #[rstest]
    fn test_energy_output_max_for_water_storage(
        two_len_simulation_time: SimulationTime,
        heat_network_service_water_storage: HeatNetworkServiceWaterStorage,
        heat_network_for_water_storage: &Arc<Mutex<HeatNetwork>>,
    ) {
        let heat_network = heat_network_for_water_storage.clone();
        let temp_return = [10.0, 15.0];
        for (t_idx, t_it) in two_len_simulation_time.iter().enumerate() {
            assert_eq!(
                heat_network_service_water_storage.energy_output_max(temp_return[t_idx], &t_it),
                [7.0, 7.0][t_idx]
            );
            heat_network.lock().timestep_end(t_idx);
        }
    }

    #[fixture]
    pub fn three_len_simulation_time() -> SimulationTime {
        SimulationTime::new(0., 3., 1.)
    }

    #[fixture]
    #[once]
    fn heat_network_for_service_space(
        three_len_simulation_time: SimulationTime,
    ) -> Arc<Mutex<HeatNetwork>> {
        let energy_supply =
            EnergySupplyBuilder::new(FuelType::MainsGas, three_len_simulation_time.total_steps())
                .build();
        let energy_supply_conn_name_auxiliary = "Boiler_auxiliary";
        let energy_supply_conn_name_building_level_distribution_losses =
            "HeatNetwork_building_level_distribution_losses";

        Arc::new(Mutex::new(HeatNetwork::new(
            5.0,
            1.0,
            0.8,
            Arc::new(RwLock::new(energy_supply)),
            energy_supply_conn_name_auxiliary.to_owned(),
            energy_supply_conn_name_building_level_distribution_losses.to_owned(),
            three_len_simulation_time.step,
        )))
    }

    #[fixture]
    fn heat_network_service_space(
        three_len_simulation_time: SimulationTime,
        heat_network_for_service_space: &Arc<Mutex<HeatNetwork>>,
    ) -> HeatNetworkServiceSpace {
        let heat_network = heat_network_for_service_space;

        let _ = HeatNetwork::create_service_connection(heat_network.clone(), "heat_network_test");

        let control = Control::SetpointTime(
            SetpointTimeControl::new(
                vec![Some(21.0), Some(21.0), None],
                0,
                1.0,
                None,
                None,
                None,
                Default::default(),
                three_len_simulation_time.step,
            )
            .unwrap(),
        );

        HeatNetworkServiceSpace::new(
            heat_network.clone(),
            "heat_network_test".to_owned(),
            Arc::new(control),
        )
    }

    #[rstest]
    fn test_heat_network_service_space(
        three_len_simulation_time: SimulationTime,
        heat_network_for_service_space: &Arc<Mutex<HeatNetwork>>,
        mut heat_network_service_space: HeatNetworkServiceSpace,
    ) {
        let heat_network = heat_network_for_service_space.clone();
        let energy_demanded = [10.0, 2.0, 2.0];
        let temp_flow = [55.0, 65.0, 65.0];
        let temp_return = [50.0, 60.0, 60.0];
        for (t_idx, t_it) in three_len_simulation_time.iter().enumerate() {
            assert_eq!(
                heat_network_service_space.demand_energy(
                    energy_demanded[t_idx],
                    temp_flow[t_idx],
                    temp_return[t_idx],
                    None,
                    None,
                    &t_it
                ),
                [5.0, 2.0, 0.0][t_idx]
            );
            heat_network.lock().timestep_end(t_idx);
        }
    }

    #[rstest]
    fn test_energy_output_max_for_service_space(
        three_len_simulation_time: SimulationTime,
        heat_network_for_service_space: &Arc<Mutex<HeatNetwork>>,
        heat_network_service_space: HeatNetworkServiceSpace,
    ) {
        let temp_output = [55.0, 65.0, 65.0];
        let temp_return_feed = [10.0, 15.0, 20.0];
        for (t_idx, t_it) in three_len_simulation_time.iter().enumerate() {
            assert_eq!(
                heat_network_service_space.energy_output_max(
                    temp_output[t_idx],
                    temp_return_feed[t_idx],
                    None,
                    &t_it
                ),
                [5.0, 5.0, 0.0][t_idx]
            );
            heat_network_for_service_space.lock().timestep_end(t_idx);
        }
    }

    #[fixture]
    #[once]
    fn energy_supply_for_heat_network(
        two_len_simulation_time: SimulationTime,
    ) -> Arc<RwLock<EnergySupply>> {
        Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Custom, two_len_simulation_time.total_steps())
                .build(),
        ))
    }

    #[fixture]
    #[once]
    fn heat_network(
        two_len_simulation_time: SimulationTime,
        energy_supply_for_heat_network: &Arc<RwLock<EnergySupply>>,
    ) -> Arc<Mutex<HeatNetwork>> {
        let energy_supply_conn_name_auxiliary = "heat_network_auxiliary";
        let energy_supply_conn_name_building_level_distribution_losses =
            "HeatNetwork_building_level_distribution_losses";

        let heat_network = Arc::new(Mutex::new(HeatNetwork::new(
            6.0,
            0.24,
            0.8,
            energy_supply_for_heat_network.clone(),
            energy_supply_conn_name_auxiliary.to_string(),
            energy_supply_conn_name_building_level_distribution_losses.to_string(),
            two_len_simulation_time.step,
        )));

        let _ = HeatNetwork::create_service_connection(heat_network.clone(), "heat_network_test");

        heat_network
    }

    #[rstest]
    fn test_energy_output_provided(
        two_len_simulation_time: SimulationTime,
        heat_network: &Arc<Mutex<HeatNetwork>>,
        energy_supply_for_heat_network: &Arc<RwLock<EnergySupply>>,
    ) {
        let heat_network = heat_network.clone();
        let energy_supply = energy_supply_for_heat_network.clone();
        let energy_output_required = [2.0, 10.0];
        for (t_idx, _) in two_len_simulation_time.iter().enumerate() {
            assert_eq!(
                heat_network.lock().demand_energy(
                    "heat_network_test",
                    energy_output_required[t_idx],
                    None,
                    None,
                    t_idx
                ),
                [2.0, 6.0][t_idx]
            );
            heat_network.lock().timestep_end(t_idx);
            assert_eq!(
                energy_supply.read().results_by_end_user()["heat_network_test"][t_idx],
                [2.0, 6.0][t_idx],
            );
            assert_eq!(
                energy_supply.read().results_by_end_user()["heat_network_auxiliary"][t_idx],
                [0.01, 0.01][t_idx]
            );
        }
    }

    #[rstest]
    fn test_hiu_loss(
        two_len_simulation_time: SimulationTime,
        heat_network: &Arc<Mutex<HeatNetwork>>,
    ) {
        let heat_network = heat_network.clone();
        for _ in two_len_simulation_time.iter() {
            assert_eq!(heat_network.lock().hiu_loss(), 0.01);
        }
    }

    #[rstest]
    fn test_building_level_distribution_losses(
        two_len_simulation_time: SimulationTime,
        heat_network: &Arc<Mutex<HeatNetwork>>,
    ) {
        let heat_network = heat_network.clone();
        for _ in two_len_simulation_time.iter() {
            assert_eq!(heat_network.lock().building_level_loss(), 0.0008);
        }
    }

    #[rstest]
    fn test_create_service_connection(heat_network: &Arc<Mutex<HeatNetwork>>) {
        let heat_network = heat_network.clone();
        let service_name = "new_service";
        assert!(!heat_network
            .lock()
            .energy_supply_connections
            .contains_key(service_name));

        HeatNetwork::create_service_connection(heat_network.clone(), service_name).unwrap();

        assert!(heat_network
            .lock()
            .energy_supply_connections
            .contains_key(service_name));

        assert!(
            HeatNetwork::create_service_connection(heat_network.clone(), service_name).is_err()
        );
    }

    #[rstest]
    fn test_create_service_hot_water_direct(
        two_len_simulation_time: SimulationTime,
        heat_network: &Arc<Mutex<HeatNetwork>>,
    ) {
        let heat_network = heat_network.clone();
        let service_name = "hot_water_direct";
        let temp_hot_water = 50.;
        let cold_feed = WaterSourceWithTemperature::ColdWaterSource(Arc::new(
            ColdWaterSource::new(vec![1.0, 1.2], 0, two_len_simulation_time.step),
        ));

        HeatNetwork::create_service_hot_water_direct(
            heat_network.clone(),
            service_name.to_owned(),
            temp_hot_water,
            cold_feed,
        );

        assert!(heat_network
            .lock()
            .energy_supply_connections
            .contains_key(service_name));
    }

    #[rstest]
    fn test_create_service_space_heating(
        three_len_simulation_time: SimulationTime,
        heat_network: &Arc<Mutex<HeatNetwork>>,
    ) {
        let heat_network = heat_network.clone();
        let service_name = "hot_water_space";

        HeatNetwork::create_service_space_heating(
            heat_network.clone(),
            service_name.to_owned(),
            Arc::new(Control::SetpointTime(
                SetpointTimeControl::new(
                    vec![Some(21.0), Some(21.0), None],
                    0,
                    1.0,
                    None,
                    None,
                    None,
                    Default::default(),
                    three_len_simulation_time.step,
                )
                .unwrap(),
            )),
        );

        assert!(heat_network
            .lock()
            .energy_supply_connections
            .contains_key(service_name));
    }

    #[rstest]
    #[ignore = "test unimplemented as would need mocked methods on an energy supply connection"]
    fn test_timestep_end() {
        // complete me (uses mocks)
    }
}
