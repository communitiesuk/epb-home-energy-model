use crate::core::common::WaterSourceWithTemperature;
use crate::core::controls::time_control::{per_control, Control, ControlBehaviour};
use crate::core::material_properties::WATER;
use crate::core::units::{HOURS_PER_DAY, WATTS_PER_KILOWATT};
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::simulation_time::SimulationTimeIteration;
use std::sync::{Arc, Mutex};

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
    pub fn new(
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

    /// Demand energy for hot water (in kWh) from the heat network
    pub fn demand_hot_water(&mut self, volume_demanded: f64, timestep_idx: usize) -> f64 {
        // Calculate energy needed to meet hot water demand
        let energy_content_kwh_per_litre = WATER.volumetric_energy_content_kWh_per_litre(
            self.temperature_hot_water,
            self.cold_feed.temperature(timestep_idx),
        );
        let energy_demand = volume_demanded * energy_content_kwh_per_litre;

        self.heat_network
            .lock()
            .expect("could not get mutex lock on heat network")
            .demand_energy(&self.service_name, energy_demand)
    }
}

/// An object to represent a water heating service provided by a heat network.
///
/// This object contains the parts of the heat network calculation that are
/// specific to providing hot water to the dwelling via a hot water cylinder.
#[derive(Clone)]
pub struct HeatNetworkServiceWaterStorage {
    heat_network: Arc<Mutex<HeatNetwork>>,
    service_name: String,
    temperature_hot_water: f64,
    control: Option<Arc<Control>>,
}

impl HeatNetworkServiceWaterStorage {
    /// Arguments:
    /// * `heat_network` - reference to the HeatNetwork object providing the service
    /// * `service_name` - name of the service demanding energy from the heat network
    /// * `temperature_hot_water` - temperature of the hot water to be provided, in deg C
    /// * `control` - optional reference to a control object
    pub fn new(
        heat_network: Arc<Mutex<HeatNetwork>>,
        service_name: String,
        temperature_hot_water: f64,
        control: Option<Arc<Control>>,
    ) -> Self {
        Self {
            heat_network,
            service_name,
            temperature_hot_water,
            control,
        }
    }

    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        if !self.is_on(simulation_time_iteration) {
            return 0.;
        }

        self.heat_network
            .lock()
            .expect("could not get mutex lock on heat network")
            .demand_energy(&self.service_name, energy_demand)
    }

    pub fn energy_output_max(&self, simulation_time_iteration: &SimulationTimeIteration) -> f64 {
        if !self.is_on(simulation_time_iteration) {
            return 0.;
        }

        self.heat_network
            .lock()
            .expect("could not get mutex lock on heat network")
            .energy_output_max(Some(self.temperature_hot_water))
    }

    fn is_on(&self, simulation_time_iteration: &SimulationTimeIteration) -> bool {
        match &self.control {
            Some(ctrl) => {
                per_control!(ctrl.as_ref(), ctrl => { ctrl.is_on(simulation_time_iteration) })
            }
            None => true,
        }
    }
}

pub struct HeatNetworkServiceSpace {
    heat_network: Arc<Mutex<HeatNetwork>>,
    service_name: String,
    control: Arc<Control>,
}

impl HeatNetworkServiceSpace {
    pub fn new(
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
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        if !self.is_on(simulation_time_iteration) {
            return 0.;
        }

        self.heat_network
            .lock()
            .expect("could not get mutex lock on heat network")
            .demand_energy(&self.service_name, energy_demand)
    }

    pub fn energy_output_max(
        &self,
        temp_output: f64,
        _temp_return_feed: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        if !self.is_on(simulation_time_iteration) {
            return 0.;
        }

        self.heat_network
            .lock()
            .expect("could not get mutex lock on heat network")
            .energy_output_max(Some(temp_output))
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

#[derive(Clone)]
pub struct HeatNetwork {
    power_max_in_kw: f64,
    daily_loss: f64,                         // in kWh/day
    building_level_distribution_losses: f64, // in watts
    // energy_supply
    total_time_running_current_timestep: f64,
    simulation_timestep: f64,
}

impl HeatNetwork {
    pub fn new(
        power_max_in_kw: f64,
        daily_loss: f64,
        building_level_distribution_losses: f64,
        simulation_timestep: f64,
    ) -> Self {
        Self {
            power_max_in_kw,
            daily_loss,
            building_level_distribution_losses,
            total_time_running_current_timestep: Default::default(),
            simulation_timestep,
        }
    }

    pub fn create_service_hot_water_direct(
        heat_network: Arc<Mutex<Self>>,
        service_name: String,
        temperature_hot_water: f64,
        cold_feed: WaterSourceWithTemperature,
    ) -> HeatNetworkServiceWaterDirect {
        // TODO create energy supply service connection for this new object

        HeatNetworkServiceWaterDirect::new(
            heat_network,
            service_name,
            temperature_hot_water,
            cold_feed,
        )
    }

    pub fn create_service_hot_water_storage(
        heat_network: Arc<Mutex<Self>>,
        service_name: String,
        temperature_hot_water: f64,
        control: Option<Arc<Control>>,
    ) -> HeatNetworkServiceWaterStorage {
        // TODO create energy supply service connection for this new object

        HeatNetworkServiceWaterStorage::new(
            heat_network,
            service_name,
            temperature_hot_water,
            control,
        )
    }

    pub fn create_service_space_heating(
        heat_network: Arc<Mutex<Self>>,
        service_name: String,
        control: Arc<Control>,
    ) -> HeatNetworkServiceSpace {
        // TODO create energy supply service connection for this new object

        HeatNetworkServiceSpace::new(heat_network, service_name, control)
    }

    /// Calculate the maximum energy output of the heat network, accounting
    /// for time spent on higher-priority services.
    fn energy_output_max(&self, _temp_output: Option<f64>) -> f64 {
        let time_available = self.simulation_timestep - self.total_time_running_current_timestep;

        self.power_max_in_kw * time_available
    }

    /// Calculate energy required by heat network to satisfy demand for the service indicated.
    pub fn demand_energy(&mut self, _service_name: &str, energy_output_required: f64) -> f64 {
        let energy_output_max = self.energy_output_max(None);
        if energy_output_max == 0. {
            return energy_output_max;
        }
        let energy_output_provided = *[
            0.,
            *[energy_output_required, energy_output_max]
                .iter()
                .max_by(|a, b| a.total_cmp(b).reverse())
                .unwrap(),
        ]
        .iter()
        .max_by(|a, b| a.total_cmp(b))
        .unwrap();

        // TODO demand energy from energy supply

        let time_available = self.simulation_timestep - self.total_time_running_current_timestep;
        self.total_time_running_current_timestep +=
            (energy_output_provided / energy_output_max) * time_available;

        energy_output_provided
    }

    pub fn timestep_end(&mut self) {
        // TODO record energy supply stuff

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
    use crate::simulation_time::SimulationTime;
    use rstest::*;

    #[fixture]
    pub fn two_len_simulation_time() -> SimulationTime {
        SimulationTime::new(0., 2., 1.)
    }

    #[fixture]
    pub fn heat_network(two_len_simulation_time: SimulationTime) -> HeatNetwork {
        HeatNetwork::new(6.0, 0.24, 0.8, two_len_simulation_time.step)
    }

    #[rstest]
    pub fn should_calc_heat_network_energy_output_provider(
        mut heat_network: HeatNetwork,
        two_len_simulation_time: SimulationTime,
    ) {
        let energy_output_required = [2.0, 10.0];
        let expected_provided = [2.0, 6.0];
        for (t_idx, _) in two_len_simulation_time.iter().enumerate() {
            assert_eq!(
                heat_network.demand_energy("heat_network_test", energy_output_required[t_idx]),
                expected_provided[t_idx]
            );
            heat_network.timestep_end();
            // TODO test energy supply stuff
        }
    }

    #[rstest]
    pub fn should_calc_correct_hiu_loss(
        heat_network: HeatNetwork,
        two_len_simulation_time: SimulationTime,
    ) {
        for _ in two_len_simulation_time.iter() {
            assert_eq!(heat_network.hiu_loss(), 0.01, "incorrect HIU loss returned");
        }
    }

    #[rstest]
    pub fn should_calc_building_level_distribution_losses(
        heat_network: HeatNetwork,
        two_len_simulation_time: SimulationTime,
    ) {
        for _ in two_len_simulation_time.iter() {
            assert_eq!(
                heat_network.building_level_loss(),
                0.0008,
                "incorrect building level distribution losses returned"
            );
        }
    }

    #[fixture]
    pub fn heat_network_for_water_direct(
        two_len_simulation_time: SimulationTime,
    ) -> Arc<Mutex<HeatNetwork>> {
        Arc::new(Mutex::new(HeatNetwork::new(
            18.0,
            1.0,
            0.8,
            two_len_simulation_time.step,
        )))
    }

    #[fixture]
    pub fn heat_network_water_direct(
        heat_network_for_water_direct: Arc<Mutex<HeatNetwork>>,
        two_len_simulation_time: SimulationTime,
    ) -> HeatNetworkServiceWaterDirect {
        // TODO register energy supply

        let cold_water_temps = [1.0, 1.2];
        let cold_feed = ColdWaterSource::new(cold_water_temps.into(), &two_len_simulation_time, 1.);
        let return_temp = 60.;
        HeatNetworkServiceWaterDirect::new(
            heat_network_for_water_direct,
            "heat_network_test".to_string(),
            return_temp,
            WaterSourceWithTemperature::ColdWaterSource(Arc::new(cold_feed)),
        )
    }

    #[rstest]
    pub fn should_calc_demand_hot_water_for_water_direct(
        mut heat_network_water_direct: HeatNetworkServiceWaterDirect,
        mut heat_network_for_water_direct: Arc<Mutex<HeatNetwork>>,
        two_len_simulation_time: SimulationTime,
    ) {
        let volume_demanded = [50.0, 100.0];
        let expected_demand = [3.429, 6.834];
        for (t_idx, _) in two_len_simulation_time.iter().enumerate() {
            assert_eq!(
                round_by_precision(
                    heat_network_water_direct.demand_hot_water(volume_demanded[t_idx], t_idx),
                    1e3
                ),
                round_by_precision(expected_demand[t_idx], 1e3),
            );
            heat_network_for_water_direct.lock().unwrap().timestep_end();
        }
    }

    #[fixture]
    pub fn heat_network_for_water_storage(
        two_len_simulation_time: SimulationTime,
    ) -> Arc<Mutex<HeatNetwork>> {
        Arc::new(Mutex::new(HeatNetwork::new(
            7.0,
            1.0,
            0.8,
            two_len_simulation_time.step,
        )))
    }

    #[fixture]
    pub fn heat_network_water_storage(
        heat_network_for_water_storage: Arc<Mutex<HeatNetwork>>,
    ) -> HeatNetworkServiceWaterStorage {
        // TODO register supply

        HeatNetworkServiceWaterStorage::new(
            heat_network_for_water_storage,
            "heat_network_test".to_string(),
            60.,
            None,
        )
    }

    #[rstest]
    pub fn should_calc_demand_energy_for_water_storage(
        mut heat_network_water_storage: HeatNetworkServiceWaterStorage,
        two_len_simulation_time: SimulationTime,
    ) {
        let energy_demanded = [10.0, 2.0];
        let expected = [7.0, 2.0];
        for (t_idx, t_it) in two_len_simulation_time.iter().enumerate() {
            assert_eq!(
                heat_network_water_storage.demand_energy(energy_demanded[t_idx], &t_it),
                expected[t_idx]
            );
            heat_network_water_storage
                .heat_network
                .lock()
                .unwrap()
                .timestep_end();
        }
    }

    #[fixture]
    pub fn three_len_simulation_time() -> SimulationTime {
        SimulationTime::new(0., 3., 1.)
    }

    #[fixture]
    pub fn heat_network_for_service_space(
        three_len_simulation_time: SimulationTime,
    ) -> Arc<Mutex<HeatNetwork>> {
        Arc::new(Mutex::new(HeatNetwork::new(
            5.0,
            1.0,
            0.8,
            three_len_simulation_time.step,
        )))
    }

    #[fixture]
    pub fn heat_network_service_space(
        heat_network_for_service_space: Arc<Mutex<HeatNetwork>>,
        three_len_simulation_time: SimulationTime,
    ) -> HeatNetworkServiceSpace {
        // TODO register service supply

        HeatNetworkServiceSpace::new(
            heat_network_for_service_space,
            "heat_network_test".to_string(),
            Arc::new(Control::SetpointTimeControl(
                SetpointTimeControl::new(
                    vec![Some(21.0), Some(21.0), None],
                    0,
                    1.0,
                    None,
                    None,
                    None,
                    None,
                    three_len_simulation_time.step,
                )
                .unwrap(),
            )),
        )
    }

    #[rstest]
    pub fn should_calc_demand_energy_for_service_space(
        mut heat_network_service_space: HeatNetworkServiceSpace,
        three_len_simulation_time: SimulationTime,
    ) {
        let energy_demanded = [10., 2., 2.];
        let temp_flow = [55., 65., 65.];
        let temp_return = [50., 60., 60.];
        let expected_demand = [5.0, 2.0, 0.0];
        for (t_idx, t_it) in three_len_simulation_time.iter().enumerate() {
            assert_eq!(
                heat_network_service_space.demand_energy(
                    energy_demanded[t_idx],
                    temp_flow[t_idx],
                    temp_return[t_idx],
                    &t_it
                ),
                expected_demand[t_idx]
            );
            heat_network_service_space
                .heat_network
                .lock()
                .unwrap()
                .timestep_end();
        }
    }

    fn round_by_precision(src: f64, precision: f64) -> f64 {
        (precision * src).round() / precision
    }
}
