use crate::compare_floats::{max_of_2, min_of_2};
use crate::core::common::{WaterSupply, WaterSupplyBehaviour};
use crate::core::controls::time_control::{per_control, Control, ControlBehaviour};
use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::core::units::{HOURS_PER_DAY, WATTS_PER_KILOWATT};
use crate::core::water_heat_demand::misc::{water_demand_to_kwh, WaterEventResult};
use crate::simulation_time::SimulationTimeIteration;
use anyhow::bail;
use indexmap::IndexMap;
use parking_lot::{Mutex, RwLock};
use smartstring::alias::String;
use std::ops::Deref;
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct HeatNetworkServiceWaterDirect {
    heat_network: Arc<Mutex<HeatNetwork>>,
    service_name: String,
    temperature_hot_water: f64, // in C
    cold_feed: WaterSupply,
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
        cold_feed: WaterSupply,
    ) -> Self {
        Self {
            heat_network,
            service_name,
            temperature_hot_water,
            cold_feed,
        }
    }

    pub(crate) fn get_cold_water_source(&self) -> &WaterSupply {
        &self.cold_feed
    }

    pub fn get_temp_hot_water(
        &self,
        volume_req: f64,
        _volume_req_already: Option<f64>,
    ) -> Vec<(f64, f64)> {
        // Always supplies the whole volume at the same temperature, so list has a single element
        vec![(self.temperature_hot_water, volume_req)]
    }

    /// Demand energy for hot water (in kWh) from the heat network
    pub(crate) fn demand_hot_water(
        &self,
        usage_events: Vec<WaterEventResult>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let mut energy_demand = 0.;

        for event in usage_events {
            if is_close!(event.volume_hot, 0., rel_tol = 1e-09, abs_tol = 1e-10) {
                continue;
            }
            let list_temp_volume = self.cold_feed.draw_off_water(event.volume_hot, simtime)?;
            let sum_t_by_v: f64 = list_temp_volume.iter().map(|(t, v)| t * v).sum();
            let sum_v: f64 = list_temp_volume.iter().map(|(_, v)| v).sum();

            let temp_cold_water = sum_t_by_v / sum_v;

            energy_demand += water_demand_to_kwh(
                event.volume_hot,
                self.temperature_hot_water,
                temp_cold_water,
            );
        }

        Ok(self.heat_network.lock().demand_energy(
            &self.service_name,
            energy_demand,
            None,
            None,
            simtime.index,
        ))
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
    _control_min: Arc<Control>,
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
        control_min: Arc<Control>, // in Python this is ControlSetPoint
        control_max: Arc<Control>, // in Python this is ControlSetPoint
    ) -> Self {
        let control = control_min.clone();

        Self {
            heat_network,
            service_name,
            control,
            _control_min: control_min,
            control_max,
        }
    }

    /// Return setpoint (not necessarily temperature)
    pub(crate) fn setpnt(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (Option<f64>, Option<f64>) {
        (
            self._control_min.setpnt(&simulation_time_iteration),
            self.control_max.setpnt(&simulation_time_iteration),
        )
    }

    pub fn demand_energy(
        &self,
        energy_demand: f64,
        _temp_flow: f64,
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
        _temp_flow: f64,
        _temp_return: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        if !self.is_on(simulation_time_iteration) {
            return 0.;
        }

        self.heat_network.lock().energy_output_max(None)
    }

    fn is_on(&self, simulation_time_iteration: &SimulationTimeIteration) -> bool {
        per_control!(self.control.as_ref(), ctrl => { ctrl.is_on(simulation_time_iteration) })
    }
}

#[derive(Clone, Debug)]
pub struct HeatNetworkServiceSpace {
    heat_network: Arc<Mutex<HeatNetwork>>,
    service_name: String,
    control: Arc<Control>,
}

impl HeatNetworkServiceSpace {
    pub(crate) fn new(
        heat_network: Arc<Mutex<HeatNetwork>>,
        service_name: String,
        control: Arc<Control>, // in Python this is ControlSetPoint
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
        _temp_output: f64,
        _temp_return_feed: f64,
        time_start: Option<f64>,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        let time_start = time_start.unwrap_or(0.);

        if !self.is_on(simulation_time_iteration) {
            return 0.;
        }

        self.heat_network.lock().energy_output_max(Some(time_start))
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
        per_control!(&self.control.as_ref(), ctrl => { <_ as ControlBehaviour>::in_required_period(ctrl.deref(), simulation_time_iteration) })
    }

    fn is_on(&self, simulation_time_iteration: &SimulationTimeIteration) -> bool {
        per_control!(&self.control.as_ref(), ctrl => { ctrl.is_on(simulation_time_iteration) })
    }
}

#[derive(Clone, Debug)]
pub(crate) struct HeatNetwork {
    power_max_in_kw: f64,
    daily_loss: f64, // in kWh/day
    power_circ_pump: f64,
    power_aux: f64,
    building_level_distribution_losses: f64, // in watts
    energy_supply: Arc<RwLock<EnergySupply>>,
    energy_supply_connections: IndexMap<String, EnergySupplyConnection>,
    energy_supply_connection_aux: EnergySupplyConnection,
    energy_supply_connection_building_level_distribution_losses: EnergySupplyConnection,
    total_time_running_current_timestep: f64,
    simulation_timestep: f64,
}

impl HeatNetwork {
    pub(crate) fn new(
        power_max_in_kw: f64,
        daily_loss: f64,
        power_circ_pump: f64,
        power_aux: f64,
        building_level_distribution_losses: f64,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_conn_name_auxiliary: String,
        energy_supply_conn_name_building_level_distribution_losses: String,
        simulation_timestep: f64,
    ) -> Self {
        Self {
            power_max_in_kw,
            daily_loss,
            power_circ_pump,
            power_aux,
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
            service_name.into(),
            EnergySupply::connection(energy_supply, service_name)?,
        );

        Ok(())
    }

    pub fn create_service_hot_water_direct(
        heat_network: Arc<Mutex<Self>>,
        service_name: &str,
        temperature_hot_water: f64,
        cold_feed: WaterSupply,
    ) -> HeatNetworkServiceWaterDirect {
        Self::create_service_connection(heat_network.clone(), service_name).unwrap();

        HeatNetworkServiceWaterDirect::new(
            heat_network,
            service_name.into(),
            temperature_hot_water,
            cold_feed,
        )
    }

    pub(crate) fn create_service_hot_water_storage(
        heat_network: Arc<Mutex<Self>>,
        service_name: &str,
        control_min: Arc<Control>, // in Python this is ControlSetPoint
        control_max: Arc<Control>, // in Python this is ControlSetPoint
    ) -> HeatNetworkServiceWaterStorage {
        Self::create_service_connection(heat_network.clone(), service_name).unwrap();

        HeatNetworkServiceWaterStorage::new(
            heat_network,
            service_name.into(),
            control_min,
            control_max,
        )
    }

    pub(crate) fn create_service_space_heating(
        heat_network: Arc<Mutex<Self>>,
        service_name: &str,
        control: Arc<Control>, // in Python this is ControlSetPoint
    ) -> HeatNetworkServiceSpace {
        Self::create_service_connection(heat_network.clone(), service_name).unwrap();

        HeatNetworkServiceSpace::new(heat_network, service_name.into(), control)
    }

    /// Calculate the maximum energy output of the heat network, accounting
    /// for time spent on higher-priority services.
    fn energy_output_max(&self, time_start: Option<f64>) -> f64 {
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
        let energy_output_max = self.energy_output_max(None);
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
    pub fn timestep_end(&mut self, timestep_idx: usize) -> anyhow::Result<()> {
        // Energy required to overcome losses
        self.energy_supply_connection_aux
            .demand_energy(self.hiu_loss(), timestep_idx)?;
        self.energy_supply_connection_building_level_distribution_losses
            .demand_energy(self.building_level_loss(), timestep_idx)?;

        let timestep = self.simulation_timestep;
        let time_remaining_current_timestep = timestep - self.total_time_running_current_timestep;

        self.calc_auxiliary_energy(timestep, time_remaining_current_timestep, timestep_idx)?;

        // Variables below need to be reset at the end of each timestep
        self.total_time_running_current_timestep = Default::default();

        Ok(())
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

    /// Calculation of energy from pump
    fn calc_auxiliary_energy(
        &self,
        timestep: f64,
        _time_remaining_current_timestep: f64,
        timestep_idx: usize,
    ) -> anyhow::Result<()> {
        let mut energy_aux = self.total_time_running_current_timestep * self.power_circ_pump;
        energy_aux += timestep * self.power_aux;
        self.energy_supply_connection_aux
            .demand_energy(energy_aux, timestep_idx)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::controls::time_control::{OnOffTimeControl, SetpointTimeControl};
    use crate::core::energy_supply::energy_supply::EnergySupplyBuilder;
    use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
    use crate::core::water_heat_demand::misc::WaterEventResultType;
    use crate::input::FuelType;
    use crate::simulation_time::SimulationTime;
    use approx::assert_relative_eq;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    fn two_len_simulation_time() -> SimulationTime {
        SimulationTime::new(0., 2., 1.)
    }

    const SERVICE_NAME: &str = "heat_network_service";

    #[fixture]
    fn dummy_heat_network(two_len_simulation_time: SimulationTime) -> Arc<Mutex<HeatNetwork>> {
        Arc::new(Mutex::new(HeatNetwork::new(
            0.0,
            0.0,
            0.0,
            0.,
            0.8,
            Arc::new(RwLock::new(
                EnergySupplyBuilder::new(
                    FuelType::Electricity,
                    two_len_simulation_time.total_steps(),
                )
                .build(),
            )),
            "aux".into(),
            "distro_losses".into(),
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
            Default::default(),
            Default::default(),
            two_len_simulation_time.step,
        );

        // there is no base call in Rust so use one of the concrete implementations
        let heat_network_service = HeatNetworkServiceWaterStorage::new(
            dummy_heat_network.clone(),
            SERVICE_NAME.into(),
            Arc::new(Control::SetpointTime(control.clone())),
            Arc::new(Control::SetpointTime(control.clone())),
        );
        assert!(heat_network_service.is_on(&two_len_simulation_time.iter().next().unwrap()));

        let heat_network_service_no_control = HeatNetworkServiceWaterStorage::new(
            dummy_heat_network,
            SERVICE_NAME.into(),
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
            0.,
            0.,
            0.8,
            Arc::new(RwLock::new(energy_supply)),
            energy_supply_conn_name_auxiliary.into(),
            energy_supply_conn_name_building_level_distribution_losses.into(),
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
            "heat_network_test".into(),
            return_temp,
            WaterSupply::ColdWaterSource(Arc::new(cold_feed)),
        )
    }

    #[rstest]
    fn test_heat_network_service_water(
        two_len_simulation_time: SimulationTime,
        heat_network_service_water_direct: HeatNetworkServiceWaterDirect,
        heat_network_for_water_direct: Arc<Mutex<HeatNetwork>>,
    ) {
        let expected_demand = [7.5834, 2.2279];
        let usage_events = [
            vec![
                WaterEventResult {
                    event_result_type: WaterEventResultType::Other,
                    temperature_warm: 60.,
                    #[allow(clippy::excessive_precision)]
                    volume_warm: 34.93868988826640,
                    #[allow(clippy::excessive_precision)]
                    volume_hot: 34.93868988826640,
                    event_duration: 5.0, // temporary - to be updated in 1.0.06a migration
                },
                WaterEventResult {
                    event_result_type: WaterEventResultType::Other,
                    temperature_warm: 60.,
                    #[allow(clippy::excessive_precision)]
                    volume_warm: 75.65325966014560,
                    #[allow(clippy::excessive_precision)]
                    volume_hot: 75.65325966014560,
                    event_duration: 5.0, // temporary - to be updated in 1.0.06a migration
                },
                WaterEventResult {
                    event_result_type: WaterEventResultType::Other,
                    temperature_warm: 60.,
                    volume_warm: 0.,
                    volume_hot: 0.,
                    event_duration: 5.0, // temporary - to be updated in 1.0.06a migration
                },
            ],
            vec![WaterEventResult {
                event_result_type: WaterEventResultType::Other,
                temperature_warm: 60.,
                volume_warm: 32.60190808710678,
                volume_hot: 32.60190808710678,
                event_duration: 5.0, // temporary - to be updated in 1.0.06a migration
            }],
        ];

        for (t_idx, t_it) in two_len_simulation_time.iter().enumerate() {
            assert_relative_eq!(
                heat_network_service_water_direct
                    .demand_hot_water(usage_events[t_idx].clone(), t_it)
                    .unwrap(),
                expected_demand[t_idx],
                max_relative = 1e-3
            );
            heat_network_for_water_direct
                .lock()
                .timestep_end(t_idx)
                .unwrap();
        }
    }

    #[rstest]
    fn test_get_cold_water_source(
        heat_network_service_water_direct: HeatNetworkServiceWaterDirect,
    ) {
        let expected = &heat_network_service_water_direct.cold_feed;
        let actual = heat_network_service_water_direct.get_cold_water_source();

        match (actual, expected) {
            (WaterSupply::ColdWaterSource(actual), WaterSupply::ColdWaterSource(expected)) => {
                assert_eq!(actual, expected);
            }
            _ => panic!("Expected ColdWaterSource variant"),
        }
    }

    #[rstest]
    fn test_get_temp_hot_water(heat_network_service_water_direct: HeatNetworkServiceWaterDirect) {
        assert_eq!(
            heat_network_service_water_direct.get_temp_hot_water(12., None),
            [(60., 12.)]
        );
    }

    #[rstest]
    fn test_demand_hot_water_empty_volume_demanded_target(
        heat_network_service_water_direct: HeatNetworkServiceWaterDirect,
        two_len_simulation_time: SimulationTime,
    ) {
        assert_eq!(
            heat_network_service_water_direct
                .demand_hot_water(vec![], two_len_simulation_time.iter().current_iteration())
                .unwrap(),
            0.
        );
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
            0.,
            0.,
            0.8,
            Arc::new(RwLock::new(energy_supply)),
            energy_supply_conn_name_auxiliary.into(),
            energy_supply_conn_name_building_level_distribution_losses.into(),
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
            Default::default(),
            Default::default(),
            two_len_simulation_time.step,
        );

        let control_max = SetpointTimeControl::new(
            vec![Some(60.), Some(60.), None],
            0,
            1.0,
            Default::default(),
            Default::default(),
            two_len_simulation_time.step,
        );

        HeatNetworkServiceWaterStorage::new(
            heat_network.clone(),
            "heat_network_test".into(),
            Arc::new(Control::SetpointTime(control_min)),
            Arc::new(Control::SetpointTime(control_max)),
        )
    }

    #[rstest]
    fn test_heat_network_service_water_storage(
        two_len_simulation_time: SimulationTime,
        heat_network_service_water_storage: HeatNetworkServiceWaterStorage,
        heat_network_for_water_storage: &Arc<Mutex<HeatNetwork>>,
    ) {
        let heat_network = heat_network_for_water_storage.clone();
        let energy_demanded = [10.0, 2.0];
        for (t_idx, t_it) in two_len_simulation_time.iter().enumerate() {
            assert_eq!(
                heat_network_service_water_storage.demand_energy(
                    energy_demanded[t_idx],
                    60.,
                    60.,
                    &t_it
                ),
                [7.0, 2.0][t_idx]
            );
            heat_network.lock().timestep_end(t_idx).unwrap();
        }
    }

    #[rstest]
    fn test_energy_output_max_for_water_storage(
        two_len_simulation_time: SimulationTime,
        heat_network_service_water_storage: HeatNetworkServiceWaterStorage,
        heat_network_for_water_storage: &Arc<Mutex<HeatNetwork>>,
    ) {
        let heat_network = heat_network_for_water_storage.clone();
        let temp_flow = [15.0, 20.0];
        let temp_return = [10.0, 15.0];
        for (t_idx, t_it) in two_len_simulation_time.iter().enumerate() {
            assert_eq!(
                heat_network_service_water_storage.energy_output_max(
                    temp_flow[t_idx],
                    temp_return[t_idx],
                    &t_it
                ),
                [7.0, 7.0][t_idx]
            );
            heat_network.lock().timestep_end(t_idx).unwrap();
        }
    }

    #[rstest]
    fn test_demand_energy_if_off_for_water_storage(
        heat_network_for_water_storage: &Arc<Mutex<HeatNetwork>>,
        two_len_simulation_time: SimulationTime,
    ) {
        let control = Arc::new(Control::OnOffTime(OnOffTimeControl::new(
            vec![Some(false), Some(false)],
            0,
            1.,
        )));
        let control_max = Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(60.), Some(60.), None],
            0,
            1.0,
            Default::default(),
            Default::default(),
            two_len_simulation_time.step,
        )));
        let heat_network_service_water_storage = HeatNetworkServiceWaterStorage::new(
            heat_network_for_water_storage.clone(),
            "heat_network_test".into(),
            control,
            control_max,
        );

        assert_eq!(
            heat_network_service_water_storage.demand_energy(
                10.,
                0., // Python passes None here but the parameter is unused
                0., // Python passes None here but the parameter is unused
                &two_len_simulation_time.iter().current_iteration()
            ),
            0.
        );
    }

    #[rstest]
    fn test_energy_output_max_if_off_for_water_storage(
        heat_network_for_water_storage: &Arc<Mutex<HeatNetwork>>,
        two_len_simulation_time: SimulationTime,
    ) {
        let control = Arc::new(Control::OnOffTime(OnOffTimeControl::new(
            vec![Some(false), Some(false)],
            0,
            1.,
        )));
        let control_max = Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(60.), Some(60.), None],
            0,
            1.0,
            Default::default(),
            Default::default(),
            two_len_simulation_time.step,
        )));
        let heat_network_service_water_storage = HeatNetworkServiceWaterStorage::new(
            heat_network_for_water_storage.clone(),
            "heat_network_test".into(),
            control,
            control_max,
        );

        assert_eq!(
            heat_network_service_water_storage.energy_output_max(
                10.,
                10.,
                &two_len_simulation_time.iter().current_iteration()
            ),
            0.
        );
    }

    #[fixture]
    fn three_len_simulation_time() -> SimulationTime {
        SimulationTime::new(0., 3., 1.)
    }

    // no longer using #[once] or Arc<Mutex<>> for this fixture
    // as they caused a race condition between tests
    #[fixture]
    fn heat_network_for_service_space(three_len_simulation_time: SimulationTime) -> HeatNetwork {
        let energy_supply =
            EnergySupplyBuilder::new(FuelType::MainsGas, three_len_simulation_time.total_steps())
                .build();
        let energy_supply_conn_name_auxiliary = "Boiler_auxiliary";
        let energy_supply_conn_name_building_level_distribution_losses =
            "HeatNetwork_building_level_distribution_losses";

        HeatNetwork::new(
            5.0,
            1.0,
            0.,
            0.,
            0.8,
            Arc::new(RwLock::new(energy_supply)),
            energy_supply_conn_name_auxiliary.into(),
            energy_supply_conn_name_building_level_distribution_losses.into(),
            three_len_simulation_time.step,
        )
    }

    #[fixture]
    fn heat_network_service_space(
        three_len_simulation_time: SimulationTime,
        heat_network_for_service_space: HeatNetwork,
    ) -> HeatNetworkServiceSpace {
        let heat_network = Arc::new(Mutex::new(heat_network_for_service_space.clone()));

        let _ = HeatNetwork::create_service_connection(heat_network.clone(), "heat_network_test");

        let control = Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(21.0), Some(21.0), None],
            0,
            1.0,
            Default::default(),
            Default::default(),
            three_len_simulation_time.step,
        ));

        HeatNetworkServiceSpace::new(
            heat_network.clone(),
            "heat_network_test".into(),
            Arc::new(control),
        )
    }

    #[rstest]
    fn test_heat_network_service_space(
        three_len_simulation_time: SimulationTime,
        mut heat_network_service_space: HeatNetworkServiceSpace,
    ) {
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
                [5.0, 0.0, 0.0][t_idx]
            );
        }
    }

    #[rstest]
    fn test_energy_output_max_for_service_space(
        three_len_simulation_time: SimulationTime,
        heat_network_for_service_space: HeatNetwork,
        heat_network_service_space: HeatNetworkServiceSpace,
    ) {
        let mut heat_network_for_service_space = heat_network_for_service_space.clone();
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
            heat_network_for_service_space.timestep_end(t_idx).unwrap();
        }
    }

    #[rstest]
    fn test_temp_setpnt_for_service_space(
        heat_network_service_space: HeatNetworkServiceSpace,
        three_len_simulation_time: SimulationTime,
    ) {
        for (t_idx, t_it) in three_len_simulation_time.iter().enumerate() {
            assert_eq!(
                heat_network_service_space.temperature_setpnt(&t_it),
                [Some(21.), Some(21.), None][t_idx]
            );
        }
    }

    #[rstest]
    fn test_in_required_period_for_service_space(
        heat_network_service_space: HeatNetworkServiceSpace,
        three_len_simulation_time: SimulationTime,
    ) {
        for (t_idx, t_it) in three_len_simulation_time.iter().enumerate() {
            assert_eq!(
                heat_network_service_space.in_required_period(&t_it),
                [Some(true), Some(true), Some(false)][t_idx]
            );
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
            0.08,
            0.04,
            0.8,
            energy_supply_for_heat_network.clone(),
            energy_supply_conn_name_auxiliary.into(),
            energy_supply_conn_name_building_level_distribution_losses.into(),
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
            heat_network.lock().timestep_end(t_idx).unwrap();
            assert_eq!(
                energy_supply.read().results_by_end_user()["heat_network_test"][t_idx],
                [2.0, 6.0][t_idx],
            );
            assert_relative_eq!(
                energy_supply.read().results_by_end_user()["heat_network_auxiliary"][t_idx],
                [0.07666, 0.13][t_idx],
                max_relative = 1e-4
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
        let cold_feed = WaterSupply::ColdWaterSource(Arc::new(ColdWaterSource::new(
            vec![1.0, 1.2],
            0,
            two_len_simulation_time.step,
        )));

        HeatNetwork::create_service_hot_water_direct(
            heat_network.clone(),
            service_name,
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
            service_name,
            Arc::new(Control::SetpointTime(SetpointTimeControl::new(
                vec![Some(21.0), Some(21.0), None],
                0,
                1.0,
                Default::default(),
                Default::default(),
                three_len_simulation_time.step,
            ))),
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

    #[rstest]
    fn test_create_service_hot_water_storage(
        three_len_simulation_time: SimulationTime,
        heat_network: &Arc<Mutex<HeatNetwork>>,
    ) {
        let control_min = Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(52.), Some(52.), None],
            0,
            1.0,
            Default::default(),
            Default::default(),
            three_len_simulation_time.step,
        )));

        let control_max = Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(60.), Some(60.), None],
            0,
            1.0,
            Default::default(),
            Default::default(),
            three_len_simulation_time.step,
        )));

        let water_storage = HeatNetwork::create_service_hot_water_storage(
            heat_network.clone(),
            "name",
            control_min,
            control_max,
        );
        assert_eq!(
            water_storage.setpnt(three_len_simulation_time.iter().current_iteration()),
            (Some(52.), Some(60.))
        );
    }

    #[rstest]
    fn test_demand_energy_with_zero_power(
        two_len_simulation_time: SimulationTime,
        energy_supply_for_heat_network: &Arc<RwLock<EnergySupply>>,
    ) {
        let energy_supply_conn_name_auxiliary = "heat_network_auxiliary_new";
        let energy_supply_conn_name_building_level_distribution_losses =
            "HeatNetwork_building_level_distribution_losses_new";

        let heat_network = Arc::new(Mutex::new(HeatNetwork::new(
            0.,
            0.24,
            0.,
            0.,
            0.8,
            energy_supply_for_heat_network.clone(),
            energy_supply_conn_name_auxiliary.into(),
            energy_supply_conn_name_building_level_distribution_losses.into(),
            two_len_simulation_time.step,
        )));

        let _ =
            HeatNetwork::create_service_connection(heat_network.clone(), "heat_network_test_new");
        assert_eq!(
            heat_network
                .lock()
                .demand_energy("name", 10., None, None, 1),
            0.
        );
    }
}
