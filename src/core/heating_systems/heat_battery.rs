/// This module provides object(s) to model the behaviour of heat batteries.
use crate::core::common::WaterSourceWithTemperature;
use crate::core::controls::time_control::{per_control, Control, ControlBehaviour};
use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::core::heating_systems::heat_pump::ServiceResult;
use crate::core::material_properties::WATER;
use crate::core::schedule::TypedScheduleEvent;
use crate::core::units::{
    KILOJOULES_PER_KILOWATT_HOUR, SECONDS_PER_HOUR, SECONDS_PER_MINUTE, WATTS_PER_KILOWATT,
};
use crate::input::HeatSourceWetDetails;
use crate::simulation_time::{SimulationTimeIteration, SimulationTimeIterator};
use anyhow::bail;
use indexmap::IndexMap;
use itertools::Itertools;
use parking_lot::{Mutex, RwLock};
use std::collections::HashMap;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub enum ServiceType {
    WaterRegular,
    Space,
}

pub enum OperationMode {
    Normal,
    OnlyCharging,
    Losses,
}

/// An object to represent a water heating service provided by a regular heat battery.
///
/// This object contains the parts of the heat battery calculation that are
/// specific to providing hot water.
#[derive(Clone, Debug)]
pub struct HeatBatteryServiceWaterRegular {
    heat_battery: Arc<Mutex<HeatBattery>>,
    service_name: String,
    cold_feed: WaterSourceWithTemperature,
    control: Arc<Control>,
    _control_min: Arc<Control>,
    control_max: Arc<Control>,
}

impl HeatBatteryServiceWaterRegular {
    /// Arguments:
    /// * `heat_battery` - reference to the Heat Battery object providing the service
    /// * `service_name` - name of the service demanding energy
    /// * `control_min` - reference to a control object which must select current the minimum timestep temperature
    /// * `control_max` - reference to a control object which must select current the maximum timestep temperature
    pub(crate) fn new(
        heat_battery: Arc<Mutex<HeatBattery>>,
        service_name: String,
        cold_feed: WaterSourceWithTemperature,
        control_min: Arc<Control>,
        control_max: Arc<Control>,
    ) -> Self {
        let control = control_min.clone();

        Self {
            heat_battery,
            service_name,
            cold_feed,
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

    fn get_cold_water_source(&self) -> &WaterSourceWithTemperature {
        &self.cold_feed
    }

    fn get_temp_hot_water(&self, simulation_time_iteration: SimulationTimeIteration) -> f64 {
        let volume = 20.; // Nominal volumen to calculate water temperature from battery
        let inlet_temp = self.cold_feed.temperature(simulation_time_iteration, None);

        self.heat_battery
            .lock()
            .get_temp_hot_water(inlet_temp, volume)
    }

    fn demand_hot_water(
        &self,
        usage_events: Option<Vec<TypedScheduleEvent>>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let mut energy_demand = 0.;
        let mut temp_hot_water = 0.;

        // Filtering out IES events that don't get added a 'warm_volume' when processing
        // the dhw_demand calculation
        let filtered_events = usage_events
            .clone()
            .into_iter()
            .flatten()
            .filter(|e| e.warm_volume.is_some())
            .collect_vec();

        for event in filtered_events {
            let warm_temp = event.temperature;
            let warm_volume = event.warm_volume;

            if warm_temp > temp_hot_water {
                temp_hot_water = warm_temp;
            }

            let energy_content_kwh_per_litre = WATER.volumetric_energy_content_kwh_per_litre(
                warm_temp,
                self.cold_feed.temperature(simulation_time_iteration, None),
            );

            energy_demand += warm_volume.unwrap() * energy_content_kwh_per_litre;
        }

        let service_on = self.is_on(simulation_time_iteration);

        if !service_on {
            energy_demand = 0.;
        }

        self.heat_battery.lock().demand_energy(
            &*self.service_name,
            ServiceType::WaterRegular,
            energy_demand,
            self.cold_feed.temperature(simulation_time_iteration, None),
            Some(temp_hot_water),
            service_on,
            None,
            Some(true),
            simulation_time_iteration.index,
        )
    }

    /// Demand energy (in kWh) from the heat_battery
    pub fn demand_energy(
        &self,
        energy_demand: f64,
        temp_return: f64,
        update_heat_source_state: Option<bool>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let service_on = self.is_on(simulation_time_iteration);
        let energy_demand = if !service_on { 0.0 } else { energy_demand };
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);

        self.heat_battery.lock().demand_energy(
            &self.service_name,
            ServiceType::WaterRegular,
            energy_demand,
            temp_return,
            self.control_max.setpnt(&simulation_time_iteration),
            service_on,
            None,
            Some(update_heat_source_state),
            simulation_time_iteration.index,
        )
    }

    pub fn energy_output_max(
        &self,
        temp_flow: f64,
        _temp_return: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        let service_on = self.is_on(simulation_time_iteration);
        if !service_on {
            return 0.;
        }

        self.heat_battery.lock().energy_output_max(temp_flow, None)
    }

    fn is_on(&self, simulation_time_iteration: SimulationTimeIteration) -> bool {
        per_control!(self.control.as_ref(), ctrl => { ctrl.is_on(&simulation_time_iteration) })
    }
}

#[derive(Clone, Debug)]
pub struct HeatBatteryServiceSpace {
    heat_battery: Arc<Mutex<HeatBattery>>,
    service_name: String,
    control: Arc<Control>,
}

/// An object to represent a space heating service provided by a heat_battery to e.g. radiators.
///
/// This object contains the parts of the heat battery calculation that are
/// specific to providing space heating.
impl HeatBatteryServiceSpace {
    pub(crate) fn new(
        heat_battery: Arc<Mutex<HeatBattery>>,
        service_name: String,
        control: Arc<Control>,
    ) -> Self {
        Self {
            heat_battery,
            service_name,
            control,
        }
    }

    pub fn temp_setpnt(&self, simulation_time_iteration: SimulationTimeIteration) -> Option<f64> {
        per_control!(self.control.as_ref(), ctrl => { ctrl.setpnt(&simulation_time_iteration) })
    }

    pub fn in_required_period(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> Option<bool> {
        per_control!(self.control.as_ref(), ctrl => { ctrl.in_required_period(&simulation_time_iteration) })
    }

    /// Demand energy (in kWh) from the heat battery
    pub fn demand_energy(
        &self,
        energy_demand: f64,
        temp_flow: f64,
        temp_return: f64,
        time_start: Option<f64>,
        update_heat_source_state: Option<bool>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let _time_start = time_start.unwrap_or(0.);
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);
        let service_on = self.is_on(simulation_time_iteration);

        let energy_demand = if !service_on { 0.0 } else { energy_demand };

        self.heat_battery.lock().demand_energy(
            &self.service_name,
            ServiceType::Space,
            energy_demand,
            temp_return,
            Some(temp_flow),
            service_on,
            None,
            Some(update_heat_source_state),
            simulation_time_iteration.index,
        )
    }

    fn is_on(&self, simulation_time_iteration: SimulationTimeIteration) -> bool {
        per_control!(self.control.as_ref(), ctrl => { ctrl.is_on(&simulation_time_iteration) })
    }

    pub(crate) fn energy_output_max(
        &self,
        temp_output: f64,
        _temp_return_feed: f64,
        time_start: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        let time_start = time_start.unwrap_or(0.);

        if !self.is_on(simtime) {
            return 0.;
        }

        self.heat_battery
            .lock()
            .energy_output_max(temp_output, Some(time_start))
    }
}

const HEAT_BATTERY_TIME_UNIT: u32 = SECONDS_PER_HOUR;
const N_ZONES: usize = 8;
const HB_TIME_STEP: f64 = 20.;
const MINIMUM_TIME_REQUIRED_TO_RUN: f64 = 120.;
const INITIAL_INLET_TEMP: f64 = 10.;
const ESTIMATED_OUTLET_TEMP: f64 = 53.;

// nothing seems to read this - check upstream whether service_results field is necessary
#[derive(Clone, Debug)]
#[allow(dead_code)]
struct HeatBatteryResult {
    service_name: String,
    service_type: ServiceType,
    service_on: bool,
    energy_output_required: f64,
    temp_output: Option<f64>,
    temp_inlet: f64,
    time_running: f64,
    energy_left_in_pipe: f64,
    temperature_left_in_pipe: f64,
    energy_delivered_hb: f64,
    energy_delivered_backup: f64,
    energy_delivered_total: f64,
    energy_delivered_low_temp: f64,
    energy_charged_during_service: f64,
    hb_zone_temperatures: Vec<f64>,
    current_hb_power: f64,
}

#[derive(Clone, Debug)]
pub struct HeatBattery {
    simulation_time: Arc<SimulationTimeIterator>,
    energy_supply: Arc<RwLock<EnergySupply>>,
    energy_supply_connection: EnergySupplyConnection,
    energy_supply_connections: HashMap<String, EnergySupplyConnection>,
    pwr_in: f64,
    max_rated_losses: f64,
    power_circ_pump: f64,
    power_standby: f64,
    n_units: usize,
    charge_control: Arc<Control>, // ChargeControl variant expected
    // nothing external seems to read this - check upstream whether service_results field is necessary
    service_results: Vec<HeatBatteryResult>,
    total_time_running_current_timestep: f64,
    flag_first_call: bool,
    charge_level: f64,
    pipe_energy: IndexMap<String, IndexMap<String, f64>>,
    energy_charged: f64,
    simultaneous_charging_and_discharging: bool,
    max_temp_of_charge: f64,
    zone_temp_c_dist_initial: Vec<f64>,
    heat_storage_kj_per_k_above: f64,
    heat_storage_kj_per_k_below: f64,
    heat_storage_kj_per_k_during: f64,
    phase_transition_temperature_upper: f64,
    phase_transition_temperature_lower: f64,
    velocity_in_hex_tube: f64,
    capillary_diameter_m: f64,
    a: f64,
    b: f64,
    heat_exchanger_surface_area_m2: f64,
    flow_rate_l_per_min: f64,
    detailed_results: Option<Vec<Arc<RwLock<Vec<ServiceResult>>>>>,
}

impl HeatBattery {
    pub(crate) fn new(
        heat_battery_details: &HeatSourceWetDetails,
        charge_control: Arc<Control>,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_connection: EnergySupplyConnection,
        simulation_time: Arc<SimulationTimeIterator>,
        output_detailed_results: Option<bool>,
    ) -> Self {
        let output_detailed_results = output_detailed_results.unwrap_or(false);

        let (
            pwr_in,
            max_rated_losses,
            power_circ_pump,
            power_standby,
            n_units,
            simultaneous_charging_and_discharging,
            max_temp_of_charge,
            heat_storage_kj_per_k_above,
            heat_storage_kj_per_k_below,
            heat_storage_kj_per_k_during,
            phase_transition_temperature_upper,
            phase_transition_temperature_lower,
            velocity_in_hex_tube,
            capillary_diameter_m,
            a,
            b,
            heat_exchanger_surface_area_m2,
            flow_rate_l_per_min,
            ..,
        ) = if let HeatSourceWetDetails::HeatBattery {
            rated_charge_power: pwr_in,
            max_rated_losses,
            electricity_circ_pump: power_circ_pump,
            electricity_standby: power_standby,
            number_of_units: n_units,
            simultaneous_charging_and_discharging,
            max_temperature,
            heat_storage_kj_per_k_above,
            heat_storage_kj_per_k_below,
            heat_storage_kj_per_k_during,
            phase_transition_temperature_upper,
            phase_transition_temperature_lower,
            velocity_in_hex_tube,
            capillary_diameter_m,
            a,
            b,
            heat_exchanger_surface_area_m2,
            flow_rate_l_per_min,
            ..
        } = heat_battery_details
        {
            (
                *pwr_in,
                *max_rated_losses,
                *power_circ_pump,
                *power_standby,
                *n_units,
                *simultaneous_charging_and_discharging,
                *max_temperature,
                *heat_storage_kj_per_k_above,
                *heat_storage_kj_per_k_below,
                *heat_storage_kj_per_k_during,
                *phase_transition_temperature_upper,
                *phase_transition_temperature_lower,
                *velocity_in_hex_tube,
                *capillary_diameter_m,
                *a,
                *b,
                *heat_exchanger_surface_area_m2,
                *flow_rate_l_per_min,
            )
        } else {
            unreachable!()
        };

        let detailed_results = if output_detailed_results {
            Some(vec![])
        } else {
            None
        };

        Self {
            simulation_time,
            energy_supply,
            energy_supply_connection,
            energy_supply_connections: Default::default(),
            pwr_in,
            max_rated_losses,
            power_circ_pump,
            power_standby,
            n_units,
            charge_control,
            service_results: Default::default(),
            total_time_running_current_timestep: Default::default(),
            flag_first_call: true,
            charge_level: Default::default(),
            pipe_energy: Default::default(),
            energy_charged: Default::default(),
            simultaneous_charging_and_discharging,
            max_temp_of_charge,
            zone_temp_c_dist_initial: vec![max_temp_of_charge; N_ZONES],
            heat_storage_kj_per_k_above,
            heat_storage_kj_per_k_below,
            heat_storage_kj_per_k_during,
            phase_transition_temperature_upper,
            phase_transition_temperature_lower,
            velocity_in_hex_tube,
            capillary_diameter_m,
            a,
            b,
            heat_exchanger_surface_area_m2,
            flow_rate_l_per_min,
            detailed_results,
        }
    }

    pub fn create_service_connection(
        heat_battery: Arc<Mutex<Self>>,
        service_name: &str,
    ) -> anyhow::Result<()> {
        if heat_battery
            .lock()
            .energy_supply_connections
            .contains_key(service_name)
        {
            bail!("Error: Service name already used: {service_name}");
        }
        let energy_supply = heat_battery.lock().energy_supply.clone();

        // Set up EnergySupplyConnection for this service
        heat_battery.lock().energy_supply_connections.insert(
            service_name.to_string(),
            EnergySupply::connection(energy_supply, service_name).unwrap(),
        );

        Ok(())
    }

    /// Return a HeatBatteryServiceWaterRegular object and create an EnergySupplyConnection for it
    ///
    /// Arguments:
    /// * `heat_battery` - reference to heat battery
    /// * `service_name` - name of the service demanding energy from the heat battery
    /// * `cold_feed` - reference to ColdWaterSource object
    /// * `control_min` - reference to a control object which must select current the minimum timestep temperature
    /// * `control_max` - reference to a control object which must select current the maximum timestep temperature
    pub(crate) fn create_service_hot_water_regular(
        heat_battery: Arc<Mutex<Self>>,
        service_name: &str,
        cold_feed: WaterSourceWithTemperature,
        control_min: Arc<Control>,
        control_max: Arc<Control>,
    ) -> HeatBatteryServiceWaterRegular {
        Self::create_service_connection(heat_battery.clone(), service_name).unwrap();
        HeatBatteryServiceWaterRegular::new(
            heat_battery,
            service_name.to_string(),
            cold_feed,
            control_min,
            control_max,
        )
    }

    pub(crate) fn create_service_space_heating(
        heat_battery: Arc<Mutex<Self>>,
        service_name: &str,
        control: Arc<Control>,
    ) -> HeatBatteryServiceSpace {
        Self::create_service_connection(heat_battery.clone(), service_name).unwrap();
        HeatBatteryServiceSpace::new(heat_battery, service_name.to_string(), control)
    }

    /// Calculates power required for unit
    ///
    /// Arguments
    /// * `time` - current time period that we are looking at
    /// * `simulation_time_iteration` - an iteration of the contextual simulation time
    ///
    /// returns -- Power required in watts
    fn electric_charge(&self) -> f64 {
        if per_control!(self.charge_control.as_ref(), ctrl => { ctrl.is_on(&self.simulation_time.current_iteration()) })
        {
            self.pwr_in
        } else {
            0.0
        }
    }

    /// Calculate time available for the current service
    fn time_available(&self, time_start: f64, timestep: f64) -> f64 {
        // Assumes that time spent on other services is evenly spread throughout
        // the timestep so the adjustment for start time below is a proportional
        // reduction of the overall time available, not simply a subtraction
        (timestep - self.total_time_running_current_timestep) * (1. - time_start / timestep)
    }

    fn calculate_heat_transfer_coeff(
        a: f64,
        b: f64,
        flow_rate_l_per_min: f64,
        reynold_number_at_1_l_per_min: f64,
    ) -> f64 {
        // Equations parameters A and B are based on test data
        // Consider adding further documentation and evidence for this in future updates
        a * (reynold_number_at_1_l_per_min * flow_rate_l_per_min).ln() + b
    }

    fn calculate_heat_transfer_kw_per_k(heat_transfer_coeff: f64, surface_area_m2: f64) -> f64 {
        (heat_transfer_coeff * surface_area_m2) / WATTS_PER_KILOWATT as f64
    }

    /// Heat transfer from heat battery zone to water flowing through it.
    ///     UAZ(n) = UA1Z(n) ------- (a) When the heat battery is discharging e.g. hot water heating mode.
    ///     UAZ(n) = UA2Z(n) ------- (b) When the heat battery is charging via external heat source
    ///     Q3Z(n) = mWCW(twoZ(n) – twiZ(n) )= UAZ(n)(TZ(n) – (twiZ(n) + twoZ(n) )/2) ----- (1)
    ///     Q3Z(n) = Heat transfer rate between PCM and the water flowing through it, (W)
    ///     mW = water mass flow rate, (kg/s)
    ///     CW = Specific heat of water, (J/(kg.K)
    ///     twoZ(n) = Water outlet temperature from zone, n, (oC)
    ///     twiZ(n) = Water inlet temperature from zone, n, (oC)
    ///     UAZ(n) = Overall heat transfer coefficient of heat exchanger in zone, n, (W/k)
    ///     TZ(n) = Heat battery zone temperature, (oC)
    /// Outlet temperature twoZ is calculated by resolving the equation (1)
    fn calculate_outlet_temp_c(
        heat_transfer_kw_per_k: f64,
        zone_temp_c: f64,
        inlet_temp_c: f64,
        flow_rate_kg_per_s: f64,
    ) -> f64 {
        (2. * heat_transfer_kw_per_k * zone_temp_c - heat_transfer_kw_per_k * inlet_temp_c
            + 2. * flow_rate_kg_per_s
                * WATER.specific_heat_capacity_kwh()
                * KILOJOULES_PER_KILOWATT_HOUR as f64
                * inlet_temp_c)
            / (2.
                * flow_rate_kg_per_s
                * WATER.specific_heat_capacity_kwh()
                * KILOJOULES_PER_KILOWATT_HOUR as f64
                + heat_transfer_kw_per_k)
    }

    /// Calculate the kinematic viscosity of water (m²/s) based on average circuit temperature.
    /// This method uses a quadratic approximation to estimate the kinematic viscosity
    /// of water as a function of the average temperature of the secondary circuit.
    /// The equation used is:
    ///     ν = a * T_avg² + b * T_avg + c
    /// where:
    ///     - ν is the kinematic viscosity in m²/s
    ///     - T_avg is the average of the inlet and outlet temperatures in °C
    ///     - a, b, c are experimentally determined coefficients:
    ///         a = 0.000000000145238
    ///         b = -0.0000000248238
    ///         c = 0.000001432
    /// These coefficients are likely derived from experimental test data or
    /// thermodynamic property tables for water within a specific temperature range
    /// relevant to secondary circuit operation.
    /// Parameters:
    ///     inlet_temp_C (float): The inlet temperature of the circuit in °C.
    ///     outlet_temp_C (float): The outlet temperature of the circuit in °C.
    /// Returns:
    ///     float: The kinematic viscosity of water in m²/s.
    /// Notes:
    ///     - This approximation is valid for the expected operating range of secondary
    ///       circuits (e.g., HVAC or hydronic systems) and may lose accuracy outside
    ///       typical temperature ranges (e.g., 0–100 °C).
    ///     - The coefficients are fixed constants based on empirical data and are not
    ///       variables in this implementation.
    fn calculate_water_kinematic_viscosity_m2_per_s(inlet_temp_c: f64, outlet_temp_c: f64) -> f64 {
        let average_temp = (inlet_temp_c + outlet_temp_c) / 2.;

        0.000000000145238 * average_temp.powf(2.) - 0.0000000248238 * average_temp + 0.000001432
    }

    fn calculate_reynold_number_at_1_l_per_min(
        water_kinematic_viscosity_m2_per_s: f64,
        velocity_in_hex_tube: f64,
        diameter_m: f64,
    ) -> f64 {
        (velocity_in_hex_tube * diameter_m) / water_kinematic_viscosity_m2_per_s
    }

    fn get_zone_properties(
        &self,
        index: usize,
        mode: OperationMode,
        zone_temp_c_dist: Vec<f64>,
        inlet_temp_c: f64,
        inlet_temp_c_zone: f64,
        q_max_kj: f64,
        reynold_number_at_1_l_per_min: f64,
        flow_rate_kg_per_s: f64,
        time_step_s: f64,
    ) -> (f64, usize, f64, f64) {
        match mode {
            OperationMode::OnlyCharging => {
                let zone_index = zone_temp_c_dist.iter().len() - index - 1;
                let zone_temp_c_start = zone_temp_c_dist[zone_index];
                (0., zone_index, zone_temp_c_start, 0.)
            }
            OperationMode::Losses => {
                let zone_temp_c_start = zone_temp_c_dist[index];
                let energy_transf = if zone_temp_c_start > inlet_temp_c {
                    q_max_kj / zone_temp_c_dist.len() as f64
                } else {
                    0.
                };
                (energy_transf, index, zone_temp_c_start, 0.)
            }
            OperationMode::Normal => {
                // NORMAL mode include battery primarily hydraulic charging or discharging with or without simultaneous electric charging
                let zone_temp_c_start = zone_temp_c_dist[index];
                let heat_transfer_coeff = Self::calculate_heat_transfer_coeff(
                    self.a,
                    self.b,
                    self.flow_rate_l_per_min,
                    reynold_number_at_1_l_per_min,
                );
                let heat_transfer_kw_per_k = Self::calculate_heat_transfer_kw_per_k(
                    heat_transfer_coeff,
                    self.heat_exchanger_surface_area_m2,
                );

                // Calculate outlet temperature and heat exchange for this zone
                let outlet_temp_c = Self::calculate_outlet_temp_c(
                    heat_transfer_kw_per_k,
                    zone_temp_c_start,
                    inlet_temp_c_zone,
                    flow_rate_kg_per_s,
                );
                let energy_transf = WATER.specific_heat_capacity_kwh()
                    * KILOJOULES_PER_KILOWATT_HOUR as f64
                    * flow_rate_kg_per_s
                    * (outlet_temp_c - inlet_temp_c_zone)
                    * time_step_s;

                (energy_transf, index, zone_temp_c_start, outlet_temp_c)
            }
        }
    }

    fn calculate_zone_energy_required(&self, zone_temp_c_start: f64, target_temp: f64) -> f64 {
        if zone_temp_c_start >= self.phase_transition_temperature_upper {
            self.heat_storage_kj_per_k_above * (zone_temp_c_start - target_temp)
        } else if zone_temp_c_start >= self.phase_transition_temperature_lower {
            if target_temp > self.phase_transition_temperature_upper {
                self.heat_storage_kj_per_k_above
                    * (self.phase_transition_temperature_upper - target_temp)
                    + self.heat_storage_kj_per_k_during
                        * (zone_temp_c_start - self.phase_transition_temperature_upper)
            } else {
                self.heat_storage_kj_per_k_during * (zone_temp_c_start - target_temp)
            }
        } else {
            if target_temp > self.phase_transition_temperature_upper {
                self.heat_storage_kj_per_k_above
                    * (self.phase_transition_temperature_upper - target_temp)
                    + self.heat_storage_kj_per_k_during
                        * (self.phase_transition_temperature_lower
                            - self.phase_transition_temperature_upper)
                    + self.heat_storage_kj_per_k_below
                        * (zone_temp_c_start - self.phase_transition_temperature_lower)
            } else if target_temp > self.phase_transition_temperature_lower {
                self.heat_storage_kj_per_k_during
                    * (self.phase_transition_temperature_lower - target_temp)
                    + self.heat_storage_kj_per_k_below
                        * (zone_temp_c_start - self.phase_transition_temperature_lower)
            } else {
                self.heat_storage_kj_per_k_below * (zone_temp_c_start - target_temp)
            }
        }
    }

    fn process_zone_simultaneous_charging(
        &self,
        zone_temp_c_start: f64,
        target_temp: f64,
        q_max_kj: f64,
        energy_transf: f64,
        energy_charged: f64,
    ) -> (f64, f64, f64) {
        let mut q_max_kj = q_max_kj;
        let mut energy_charged = energy_charged;
        let mut energy_transf = energy_transf;

        if zone_temp_c_start < target_temp {
            // zone initially below full charge
            let mut q_required =
                self.calculate_zone_energy_required(zone_temp_c_start, target_temp);

            if energy_transf >= 0. {
                // inlet water withdraws energy from battery
                if -q_max_kj >= energy_transf {
                    // Charging is enough to recover energy withdrawn and possibly more
                    q_max_kj += energy_transf;
                    energy_charged += energy_transf / KILOJOULES_PER_KILOWATT_HOUR as f64;
                    energy_transf = 0.;

                    if q_max_kj > q_required {
                        // Charging is not enough to push zone temperature to target
                        q_required = q_max_kj;
                        energy_charged += -q_max_kj / KILOJOULES_PER_KILOWATT_HOUR as f64;
                        q_max_kj = 0.;
                    } else {
                        // Charging is enough to push zone temperature to target temperature
                        q_max_kj -= q_required;
                        energy_charged += -q_required / KILOJOULES_PER_KILOWATT_HOUR as f64;
                    }
                    // Update zone temperature with energy from charging
                    energy_transf += q_required;
                } else {
                    // Charging can only recover partially the energy withdrawn
                    energy_transf += q_max_kj;
                    energy_charged += -q_max_kj / KILOJOULES_PER_KILOWATT_HOUR as f64;
                    q_max_kj = 0.;
                }
            } else {
                // inlet water adds energy to battery
                if q_max_kj + energy_transf > q_required {
                    // inlet water + charging is not enough to push zone temperature to target
                    q_required = q_max_kj + energy_transf;
                    energy_charged += -q_max_kj / KILOJOULES_PER_KILOWATT_HOUR as f64;
                    q_max_kj = 0.;

                    energy_transf = q_required;
                } else {
                    // inlet temperature + charging can take zone temperature to target temp
                    if !(energy_transf < q_required) {
                        //
                        if q_max_kj > q_required - energy_transf {
                            // Charging cannot take zone temperature to target after zone warmed by inlet water...
                            q_required = q_max_kj + energy_transf;
                            energy_charged += -q_max_kj / KILOJOULES_PER_KILOWATT_HOUR as f64;
                            q_max_kj = 0.;

                            energy_transf += q_required;
                        } else {
                            // There is plenty of charging after taking zone temperature to target
                            q_max_kj -= q_required - energy_transf;
                            energy_charged +=
                                -(q_required - energy_transf) / KILOJOULES_PER_KILOWATT_HOUR as f64;
                            energy_transf = q_required;
                        }
                    }
                }
            }
        } else {
            // zone initially fully charged
            if energy_transf >= 0. {
                // inlet water withdraws energy from battery
                if -q_max_kj > energy_transf {
                    // Charging is enough to recover energy withdrawn
                    q_max_kj += energy_transf;
                    energy_charged += energy_transf / KILOJOULES_PER_KILOWATT_HOUR as f64;
                    energy_transf = 0.;
                } else {
                    // Charging can only recover partially the energy withdrawn
                    energy_transf += q_max_kj;
                    energy_charged += -q_max_kj / KILOJOULES_PER_KILOWATT_HOUR as f64;
                    q_max_kj = 0.;
                }
            }
        }

        (q_max_kj, energy_charged, energy_transf)
    }

    fn calculate_new_zone_temperature(&self) {
        todo!("0.34")
    }

    fn process_heat_battery_zones(
        &self,
        inlet_temp_c: f64,
        zone_temp_c_dist: &Vec<f64>,
        flow_rate_kg_per_s: f64,
        time_step_s: f64,
        reynold_number: f64,
        pwr_in: Option<f64>,
    ) -> (f64, f64, Vec<f64>, f64) {
        todo!()
    }

    fn charge_battery_hydraulic(&self) {
        todo!("0.34")
    }

    fn charge_battery(&self) {
        todo!("0.34")
    }

    fn battery_heat_loss(&self) {
        todo!("0.34")
    }

    fn get_temp_hot_water(&self, inlet_temp: f64, volume: f64) -> f64 {
        todo!("0.34 migration")
    }

    /// Calculate the maximum energy output of the heat battery, accounting
    /// for time spent on higher-priority services.
    pub fn energy_output_max(&mut self, temp_output: f64, time_start: Option<f64>) -> f64 {
        // Return the energy the battery can provide assuming the HB temperature inlet
        // is constant during HEM time step equal to the required emitter temperature (temp_output)
        // Maximum energy for a given HB zones temperature distribution and inlet temperature.
        // The calculation methodology is the same as described in the demand_energy function.
        let time_start = time_start.unwrap_or(0.);
        let timestep = self.simulation_time.step_in_hours();
        let total_time_s = timestep * SECONDS_PER_HOUR as f64;

        // time_step_s for HB calculation is a sensitive inputs for the process as, the longer it is, the
        // lower the accuracy due to maintaining Reynolds number working in intervals where the properties
        // of the fluid have changed sufficiently to degrade the accuracy of the calculation.
        // This is critical for the demand_energy function but less so for the energy_output_max as this
        // only provides an estimation of the heat capacity of the battery and can be slightly overestitmated
        // with a longer time step that reduces the calculation time, which is a critical factor for HEM.
        // However, from current testing, time_step_s longer than 100 might cause instabilities in the calculation
        // leading to failure to complete. Thus, we are capping the max timestep to 100 seconds.
        let time_step_s = (HB_TIME_STEP * 5.).min(100.);

        let pwr_in = self.electric_charge();

        // Initial Reynold number
        let water_kinematic_viscosity_m2_per_s = Self::calculate_water_kinematic_viscosity_m2_per_s(
            INITIAL_INLET_TEMP,
            ESTIMATED_OUTLET_TEMP,
        );
        let reynold_number_at_1_l_per_min = Self::calculate_reynold_number_at_1_l_per_min(
            water_kinematic_viscosity_m2_per_s,
            self.velocity_in_hex_tube,
            self.capillary_diameter_m,
        );

        let flow_rate_kg_per_s =
            (self.flow_rate_l_per_min / SECONDS_PER_MINUTE as f64) * WATER.density();

        let zone_temp_c_dist = self.zone_temp_c_dist_initial.clone();
        let mut energy_delivered_hb = 0.;
        let inlet_temp_c = temp_output;
        let n_time_steps = (total_time_s / time_step_s) as usize;

        for _ in 0..n_time_steps {
            // Processing HB zones
            let (outlet_temp_c, zone_temp_c_dist, energy_transf_delivered, _) = self
                .process_heat_battery_zones(
                    inlet_temp_c,
                    &zone_temp_c_dist,
                    flow_rate_kg_per_s,
                    time_step_s,
                    reynold_number_at_1_l_per_min,
                    None,
                );

            // RN for next time step
            let water_kinematic_viscosity_m2_per_s =
                Self::calculate_water_kinematic_viscosity_m2_per_s(inlet_temp_c, outlet_temp_c);
            let reynold_number_at_1_l_per_min = Self::calculate_reynold_number_at_1_l_per_min(
                water_kinematic_viscosity_m2_per_s,
                self.velocity_in_hex_tube,
                self.capillary_diameter_m,
            );

            let energy_delivered_ts: f64 = energy_transf_delivered.iter().sum();

            if outlet_temp_c > temp_output {
                // In this new method, adjust total energy to make more real with the 6 ts we have configured
                energy_delivered_hb += energy_delivered_ts
            } else {
                break;
            }
        }

        if energy_delivered_hb < 0. {
            energy_delivered_hb = 0.
        }

        energy_delivered_hb * self.n_units as f64
    }

    fn first_call(&self) {
        todo!("0.34")
    }

    pub fn demand_energy(
        &mut self,
        service_name: &str,
        service_type: ServiceType,
        mut energy_output_required: f64,
        temp_return_feed: f64,
        temp_output: Option<f64>,
        service_on: bool,
        time_start: Option<f64>,
        update_heat_source_state: Option<bool>,
        timestep_idx: usize,
    ) -> anyhow::Result<f64> {
        // Return the energy provided by the HB during a HEM time step (assuming
        // an inlet temperature constant) and update the HB state (zones distribution temperatures)
        // The HEM time step is divided into sub-timesteps. For each sub-timestep the zones temperature are
        // calculated (loop through zones).
        let time_start = time_start.unwrap_or(0.);
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);
        let timestep = self.simulation_time.step_in_hours();
        let time_available = self.time_available(time_start, timestep);

        // demand_energy is called for each service in each timestep
        // Some calculations are only required once per timestep
        // Perform these calculations here
        if self.flag_first_call {
            self.first_call();
        }

        let pwr_in = if self.simultaneous_charging_and_discharging {
            self.electric_charge()
        } else {
            0.
        };

        // TODO 0.34 do we need a new type for pipe_energy?
        if temp_output.is_none()
            || temp_output.unwrap() <= self.pipe_energy[service_name]["temperature"]
        {
            if energy_output_required > self.pipe_energy[service_name]["energy"] {
                energy_output_required -= self.pipe_energy[service_name]["energy"];
                self.pipe_energy[service_name]["energy"] = 0.;
                self.pipe_energy[service_name]["temperature"] = 0.;
            } else {
                self.pipe_energy[service_name]["energy"] -= energy_output_required;
                energy_output_required = 0.;
            }
        }

        // Distributing energy demand through all units
        let energy_demand = energy_output_required / self.n_units as f64;

        // Initial Reynold number
        let water_kinematic_viscosity_m2_per_s = Self::calculate_water_kinematic_viscosity_m2_per_s(
            INITIAL_INLET_TEMP,
            ESTIMATED_OUTLET_TEMP,
        );
        let reynold_number_at_1_l_per_min = Self::calculate_reynold_number_at_1_l_per_min(
            water_kinematic_viscosity_m2_per_s,
            self.velocity_in_hex_tube,
            self.capillary_diameter_m,
        );

        let flow_rate_kg_per_s =
            (self.flow_rate_l_per_min / SECONDS_PER_MINUTE as f64) * WATER.density();

        let mut energy_delivered_hb = 0.;
        let mut total_energy_low_temp = 0.;
        let inlet_temp_c = temp_return_feed;
        let zone_temp_c_dist = self.zone_temp_c_dist_initial.clone();

        if energy_output_required <= 0. {
            if update_heat_source_state {
                self.service_results.push(HeatBatteryResult {
                    service_name: service_name.to_string(),
                    service_type,
                    service_on,
                    energy_output_required,
                    temp_output,
                    temp_inlet: temp_return_feed,
                    time_running: 0.,
                    energy_left_in_pipe: self.pipe_energy[service_name]["energy"],
                    temperature_left_in_pipe: self.pipe_energy[service_name]["temperature"],
                    energy_delivered_hb: 0.,
                    energy_delivered_backup: 0.,
                    energy_delivered_total: 0.,
                    energy_delivered_low_temp: 0.,
                    energy_charged_during_service: 0.,
                    hb_zone_temperatures: zone_temp_c_dist,
                    current_hb_power: Default::default(),
                });
            }
            return Ok(energy_delivered_hb);
        }

        let mut time_step_s = 1.;
        let mut time_running_current_service = 0.;

        let mut flag_minimum_run = false; // False: supply energy to emitter; True: running water to complete loop
        let mut energy_charged = 0.;
        let mut time_extra = Default::default();

        while time_step_s > 0. {
            // Processing HB zones
            let (
                outlet_temp_c,
                zone_temp_c_dist,
                energy_transf_delivered,
                energy_charged_during_battery_time_step,
            ) = self.process_heat_battery_zones(
                inlet_temp_c,
                &zone_temp_c_dist,
                flow_rate_kg_per_s,
                time_step_s,
                reynold_number_at_1_l_per_min,
                Some(pwr_in),
            );

            if update_heat_source_state {
                self.energy_charged += energy_charged_during_battery_time_step;
            }
            energy_charged += energy_charged_during_battery_time_step;

            time_running_current_service += time_step_s;

            // RN for next time step
            let water_kinematic_viscosity_m2_per_s =
                Self::calculate_water_kinematic_viscosity_m2_per_s(temp_return_feed, outlet_temp_c);
            let reynold_number_at_1_l_per_min = Self::calculate_reynold_number_at_1_l_per_min(
                water_kinematic_viscosity_m2_per_s,
                self.velocity_in_hex_tube,
                self.capillary_diameter_m,
            );

            let energy_delivered_ts =
                energy_transf_delivered.iter().sum::<f64>() / KILOJOULES_PER_KILOWATT_HOUR as f64;

            if temp_output.is_none() || outlet_temp_c > temp_output.unwrap() {
                if !flag_minimum_run {
                    energy_delivered_hb += energy_delivered_ts; // demand_per_time_step_kwh
                    let max_instant_power = energy_delivered_ts / time_step_s;

                    if max_instant_power > 0. {
                        time_step_s = (energy_demand - energy_delivered_hb) / max_instant_power;
                    }
                } else {
                    if energy_delivered_ts != 0. {
                        let current_energy = self.pipe_energy[service_name]["energy"];
                        let current_temperature = self.pipe_energy[service_name]["temperature"];
                        let new_temperature = ((current_temperature * current_energy)
                            + (outlet_temp_c * energy_delivered_ts))
                            / (current_energy + energy_delivered_ts);

                        self.pipe_energy[service_name]["energy"] += energy_delivered_ts;
                        self.pipe_energy[service_name]["temperature"] = new_temperature;
                    }
                }

                if time_step_s > HB_TIME_STEP {
                    time_step_s = HB_TIME_STEP;
                }

                if (energy_demand - energy_delivered_hb) < 0.0001 {
                    // Energy supplied, run to complete water loop
                    if time_running_current_service > MINIMUM_TIME_REQUIRED_TO_RUN {
                        break;
                    }

                    if !flag_minimum_run {
                        time_extra = MINIMUM_TIME_REQUIRED_TO_RUN - time_running_current_service;
                        flag_minimum_run = true;
                    } else {
                        time_extra -= time_step_s;
                    }

                    if time_extra > HB_TIME_STEP {
                        time_step_s = HB_TIME_STEP
                    } else {
                        time_step_s = time_extra
                    }
                }

                if time_running_current_service + time_step_s
                    > time_available * SECONDS_PER_HOUR as f64
                {
                    time_step_s =
                        time_available * SECONDS_PER_HOUR as f64 - time_running_current_service
                }
            } else {
                // outlet_temp_c is below required temperature
                total_energy_low_temp += energy_delivered_ts;

                if energy_delivered_hb > 0. {
                    if time_running_current_service > MINIMUM_TIME_REQUIRED_TO_RUN {
                        break;
                    }

                    if !flag_minimum_run {
                        time_extra = MINIMUM_TIME_REQUIRED_TO_RUN - time_running_current_service;
                        flag_minimum_run = true;
                    } else {
                        time_extra -= time_step_s;
                    }

                    if time_extra > HB_TIME_STEP {
                        time_step_s = HB_TIME_STEP
                    } else {
                        time_step_s = time_extra
                    }
                } else {
                    break;
                }
            }
        }

        if update_heat_source_state {
            self.zone_temp_c_dist_initial = zone_temp_c_dist.clone();

            self.total_time_running_current_timestep +=
                time_running_current_service / SECONDS_PER_HOUR as f64;

            let current_hb_power = if time_running_current_service > 0. {
                energy_delivered_hb * SECONDS_PER_HOUR as f64 / time_running_current_service
            } else {
                Default::default()
            };
            // TODO (from Python) Clarify whether Heat Batteries can have direct electric backup if depleted
            self.service_results.push(HeatBatteryResult {
                service_name: service_name.to_string(),
                service_type,
                service_on,
                energy_output_required,
                temp_output,
                temp_inlet: temp_return_feed,
                time_running: time_running_current_service,
                energy_left_in_pipe: self.pipe_energy[service_name]["energy"],
                temperature_left_in_pipe: self.pipe_energy[service_name]["temperature"],
                energy_delivered_hb: energy_delivered_hb * self.n_units as f64,
                energy_delivered_backup: 0.,
                energy_delivered_total: energy_delivered_hb * self.n_units as f64 + 0.,
                energy_delivered_low_temp: total_energy_low_temp * self.n_units as f64,
                energy_charged_during_service: energy_charged * self.n_units as f64,
                hb_zone_temperatures: zone_temp_c_dist,
                current_hb_power,
            });
        }

        Ok(energy_delivered_hb * self.n_units as f64)
    }

    /// Calculation of heat battery auxilary energy consumption
    fn calc_auxiliary_energy(
        &self,
        _timestep: f64,
        time_remaining_current_timestep: f64,
        timestep_idx: usize,
    ) {
        // Energy used by circulation pump
        let mut energy_aux = self.total_time_running_current_timestep * self.power_circ_pump;

        // Energy used in standby mode
        energy_aux += self.power_standby * time_remaining_current_timestep;

        self.energy_supply_connection
            .demand_energy(energy_aux, timestep_idx)
            .unwrap();
    }

    /// Calculations to be done at the end of each timestep
    pub fn timestep_end(&mut self, timestep_idx: usize) -> anyhow::Result<()> {
        // TODO 0.34 - commented out/deleted code in this method when updating methods above, this method still needs to be updated
        let timestep = self.simulation_time.step_in_hours();
        let time_remaining_current_timestep = timestep - self.total_time_running_current_timestep;

        if self.flag_first_call {
            self.first_call();
        }

        // Calculating auxiliary energy to provide services during timestep
        self.calc_auxiliary_energy(timestep, time_remaining_current_timestep, timestep_idx);

        // Completing any charging left in the timestep and removing all losses from the charge level
        // Calculating heat battery losses in timestep to correct charge level
        // Currently assumed all losses are to the exterior independently of the
        // heat battery location

        let charge_level = self.charge_level;
        let target_charge = self.target_charge()?;

        // Calculate new charge level after accounting for energy in and out and cap at target_charge

        // charge_level += delta_charge_level;
        // if charge_level > target_charge {
        //     e_in -= (charge_level - target_charge) * self.heat_storage_capacity;
        //     if e_in < 0.0 {
        //         e_in = 0.;
        //         charge_level -= delta_charge_level;
        //         delta_charge_level = -e_loss / self.heat_storage_capacity;
        //         charge_level += delta_charge_level;
        //     } else {
        //         charge_level = target_charge;
        //     }
        // }

        self.charge_level = charge_level;

        let current_hour = self.simulation_time.current_hour();

        // Preparing Heat battery for next time step
        // Variables below need to be reset at the end of each timestep
        // Picking target charge level from control
        let time_range = (current_hour + 1) * HEAT_BATTERY_TIME_UNIT;

        let target_charge = self.target_charge()?;
        let charge_level_qin = self.charge_level;

        // Calculate max charge level possible in next timestep
        // if charge_level_qin < target_charge {
        //     delta_charge_level = q_in_ts * timestep / self.heat_storage_capacity;
        //     charge_level_qin += delta_charge_level;
        //     if charge_level_qin > target_charge {
        //         charge_level_qin = target_charge;
        //     }
        // }

        // Estimating output rate at average of capacity in timestep
        // let delta_charge_level = max_output * timestep / self.heat_storage_capacity;
        // self.q_out_ts =
        //     Some(self.lab_test_rated_output(charge_level_qin - delta_charge_level / 2.));
        // self.q_loss_ts = Some(self.lab_test_losses(charge_level_qin - delta_charge_level / 2.));

        self.total_time_running_current_timestep = Default::default();
        self.service_results = Default::default();

        Ok(())
    }

    fn output_detailed_results() {
        todo!("0.34")
    }

    fn target_charge(&self) -> anyhow::Result<f64> {
        match self.charge_control.as_ref() {
            Control::Charge(ctrl) => {
                ctrl.target_charge(self.simulation_time.current_iteration(), None)
            }
            _ => unreachable!(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::core::common::WaterSourceWithTemperature;
    use crate::core::controls::time_control::SetpointTimeControl;
    use crate::core::controls::time_control::{ChargeControl, Control};
    use crate::core::energy_supply::energy_supply::{
        EnergySupply, EnergySupplyBuilder, EnergySupplyConnection,
    };
    use crate::core::heating_systems::heat_battery::HeatBattery;
    use crate::core::heating_systems::heat_battery::HeatBatteryServiceSpace;
    use crate::core::heating_systems::heat_battery::HeatBatteryServiceWaterRegular;
    use crate::core::heating_systems::heat_battery::ServiceType;
    use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
    use crate::external_conditions::{DaylightSavingsConfig, ExternalConditions};
    use crate::input::HeatSourceWetDetails;
    use crate::input::{ControlLogicType, HeatSourceLocation};
    use crate::input::{ExternalSensor, FuelType};
    use crate::simulation_time::SimulationTimeIteration;
    use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
    use approx::assert_relative_eq;
    use itertools::Itertools;
    use parking_lot::{Mutex, RwLock};
    use rstest::fixture;
    use rstest::rstest;
    use serde_json::json;
    use std::sync::Arc;

    const SERVICE_NAME: &str = "TestService";

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 2., 1.)
    }

    #[fixture]
    pub fn simulation_time_iterator(
        simulation_time: SimulationTime,
    ) -> Arc<SimulationTimeIterator> {
        simulation_time.iter().into()
    }

    #[fixture]
    pub fn simulation_time_iteration(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) -> SimulationTimeIteration {
        simulation_time_iterator.current_iteration()
    }

    #[fixture]
    fn external_sensor() -> ExternalSensor {
        serde_json::from_value(json!({
            "correlation": [
                {"temperature": 0.0, "max_charge": 1.0},
                {"temperature": 10.0, "max_charge": 0.9},
                {"temperature": 18.0, "max_charge": 0.0}
            ]
        }))
        .unwrap()
    }

    #[fixture]
    fn external_conditions(simulation_time: SimulationTime) -> ExternalConditions {
        ExternalConditions::new(
            &simulation_time.iter(),
            vec![0.0, 2.5],
            vec![3.7, 3.8],
            vec![200., 220.],
            vec![333., 610.],
            vec![420., 750.],
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
            // following shading segments are corrected from upstream Python, which uses angles measured from wrong origin
            serde_json::from_value(json!(
                [
                    {"number": 1, "start360": 0, "end360": 45},
                    {"number": 2, "start360": 45, "end360": 90},
                ]
            ))
            .unwrap(),
        )
    }

    #[fixture]
    fn battery_control_off(
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
    ) -> Control {
        create_control_with_value(false, external_conditions, external_sensor)
    }

    #[fixture]
    fn battery_control_on(
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
    ) -> Control {
        create_control_with_value(true, external_conditions, external_sensor)
    }

    fn create_control_with_value(
        boolean: bool,
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
    ) -> Control {
        Control::Charge(
            ChargeControl::new(
                ControlLogicType::Manual,
                vec![boolean],
                1.,
                0,
                1.,
                vec![Some(0.2)],
                None,
                None,
                None,
                None,
                external_conditions.into(),
                Some(external_sensor),
            )
            .unwrap(),
        )
    }

    fn create_heat_battery(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        control: Control,
    ) -> Arc<Mutex<HeatBattery>> {
        let _labs_tests_rated_output = vec![
            (0.0, 0.0),
            (0.08, 0.00),
            (0.17, 0.05),
            (0.19, 0.10),
            (0.21, 0.15),
            (0.23, 0.21),
            (0.25, 0.23),
            (0.28, 0.26),
            (0.31, 0.29),
            (0.34, 0.32),
            (0.38, 0.36),
            (0.42, 0.41),
            (0.47, 0.45),
            (0.52, 0.51),
            (0.58, 0.57),
            (0.64, 0.64),
            (0.72, 0.71),
            (0.8, 0.8),
            (0.89, 0.89),
            (1.0, 1.0),
        ];
        let _labs_tests_rated_output_enhanced = vec![
            (0.0, 0.0),
            (0.101, 0.0),
            (0.12, 0.18),
            (0.144, 0.235),
            (0.175, 0.313),
            (0.215, 0.391),
            (0.266, 0.486),
            (0.328, 0.607),
            (0.406, 0.728),
            (0.494, 0.795),
            (0.587, 0.825),
            (0.683, 0.875),
            (0.781, 0.906),
            (0.891, 0.953),
            (0.981, 0.992),
            (1.0, 1.0),
        ];
        let _labs_tests_losses = vec![
            (0.0, 0.),
            (0.16, 0.13),
            (0.17, 0.15),
            (0.19, 0.17),
            (0.21, 0.18),
            (0.23, 0.21),
            (0.25, 0.23),
            (0.28, 0.26),
            (0.31, 0.29),
            (0.34, 0.32),
            (0.38, 0.36),
            (0.42, 0.41),
            (0.47, 0.45),
            (0.52, 0.51),
            (0.58, 0.57),
            (0.64, 0.64),
            (0.72, 0.71),
            (0.8, 0.8),
            (0.89, 0.89),
            (1.0, 1.0),
        ];
        let heat_battery_details: &HeatSourceWetDetails = &HeatSourceWetDetails::HeatBattery {
            energy_supply: "mains elec".to_string(),
            heat_battery_location: Some(HeatSourceLocation::Internal),
            electricity_circ_pump: 0.06,
            electricity_standby: 0.0244,
            rated_charge_power: 20.0,
            max_rated_losses: 0.22,
            number_of_units: 1,
            control_charge: "hb_charge_control".into(),
            simultaneous_charging_and_discharging: false,
            heat_storage_kj_per_k_above: 47.6875,
            heat_storage_kj_per_k_below: 38.15,
            heat_storage_kj_per_k_during: 1539.625,
            phase_transition_temperature_upper: 59.,
            phase_transition_temperature_lower: 57.,
            max_temperature: 80.,
            velocity_in_hex_tube: 0.035,
            capillary_diameter_m: 0.0065,
            a: 19.744,
            b: -105.5,
            heat_exchanger_surface_area_m2: 8.83,
            flow_rate_l_per_min: 10.,
        };

        let energy_supply: Arc<RwLock<EnergySupply>> = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time_iterator.total_steps())
                .build(),
        ));

        let energy_supply_connection: EnergySupplyConnection =
            EnergySupply::connection(energy_supply.clone(), "WaterHeating").unwrap();

        let heat_battery = Arc::new(Mutex::new(HeatBattery::new(
            heat_battery_details,
            control.into(),
            energy_supply,
            energy_supply_connection,
            simulation_time_iterator,
            None,
        )));

        HeatBattery::create_service_connection(heat_battery.clone(), SERVICE_NAME).unwrap();

        heat_battery
    }

    fn create_setpoint_time_control(schedule: Vec<Option<f64>>) -> Control {
        Control::SetpointTime(
            SetpointTimeControl::new(schedule, 0, 1., None, None, None, Default::default(), 1.)
                .unwrap(),
        )
    }

    fn get_service_names_from_results(heat_battery: Arc<Mutex<HeatBattery>>) -> Vec<String> {
        heat_battery
            .lock()
            .service_results
            .iter()
            .map(|result| result.service_name.clone())
            .collect_vec()
    }

    // in Python this test is called test_service_is_on_with_control
    #[rstest]
    fn test_service_is_on_when_service_control_is_on(
        simulation_time_iteration: SimulationTimeIteration,
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        // Test when controlvent is provided and returns True
        let service_control_on: Control = create_setpoint_time_control(vec![Some(21.0)]);

        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off);

        let heat_battery_service = HeatBatteryServiceSpace::new(
            heat_battery.clone(),
            SERVICE_NAME.into(),
            service_control_on.into(),
        );

        assert!(heat_battery_service.is_on(simulation_time_iteration));

        let service_control_off: Control = create_setpoint_time_control(vec![None]);

        let heat_battery_service: HeatBatteryServiceSpace = HeatBatteryServiceSpace::new(
            heat_battery,
            SERVICE_NAME.into(),
            service_control_off.into(),
        );

        assert!(!heat_battery_service.is_on(simulation_time_iteration));
    }

    // test_service_is_on_without_control
    #[rstest]
    fn test_service_with_no_service_control_is_always_on_for_water_regular(
        simulation_time_iteration: SimulationTimeIteration,
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off);
        let control_min = create_setpoint_time_control(vec![
            Some(52.),
            None,
            None,
            None,
            Some(52.),
            Some(52.),
            Some(52.),
            Some(52.),
        ]);
        let control_max = create_setpoint_time_control(vec![
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
        ]);
        let cold_water_source = ColdWaterSource::new(vec![1.0, 1.2], 0, 1.);
        let heat_battery_service: HeatBatteryServiceWaterRegular =
            HeatBatteryServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                WaterSourceWithTemperature::ColdWaterSource(Arc::new(cold_water_source)),
                Arc::new(control_min),
                Arc::new(control_max),
            );

        assert!(heat_battery_service.is_on(simulation_time_iteration));
    }

    #[rstest]
    #[ignore = "while migrating to 0.34"]
    fn test_demand_energy_when_service_control_on_for_water_regular(
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
        let energy_demand = 10.;
        let temp_return = 40.;

        let control_min = create_setpoint_time_control(vec![
            Some(52.),
            None,
            None,
            None,
            Some(52.),
            Some(52.),
            Some(52.),
            Some(52.),
        ]);
        let control_max = create_setpoint_time_control(vec![
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
        ]);

        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        let cold_water_source = ColdWaterSource::new(vec![1.0, 1.2], 0, 1.);
        let heat_battery_service: HeatBatteryServiceWaterRegular =
            HeatBatteryServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                WaterSourceWithTemperature::ColdWaterSource(Arc::new(cold_water_source)),
                Arc::new(control_min),
                Arc::new(control_max),
            );

        let result = heat_battery_service
            .demand_energy(energy_demand, temp_return, None, simulation_time_iteration)
            .unwrap();

        assert_relative_eq!(result, 4.358566028225806);
    }

    // In Python this is test_demand_energy_service_off
    #[rstest]
    #[ignore = "while migrating to 0.34"]
    fn test_demand_energy_returns_zero_when_service_control_is_off_for_water_regular(
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
        let energy_demand = 10.;
        let temp_return = 40.;

        let service_control_off = Arc::new(create_setpoint_time_control(vec![None]));

        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        let cold_water_source = ColdWaterSource::new(vec![1.0, 1.2], 0, 1.);
        let heat_battery_service: HeatBatteryServiceWaterRegular =
            HeatBatteryServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                WaterSourceWithTemperature::ColdWaterSource(Arc::new(cold_water_source)),
                service_control_off.clone(),
                service_control_off,
            );

        let result = heat_battery_service
            .demand_energy(energy_demand, temp_return, None, simulation_time_iteration)
            .unwrap();

        assert_eq!(result, 0.);
    }

    // In Python this is test_energy_output_max_service_on
    #[rstest]
    #[ignore = "while migrating to 0.34"]
    fn test_energy_output_max_when_service_control_on_for_water_regular(
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let temp_return = 40.;
        let battery_control_on: Control = Control::Charge(
            ChargeControl::new(
                ControlLogicType::Manual,
                vec![true, true, true],
                1.,
                0,
                1.,
                [1.5, 1.6].into_iter().map(Into::into).collect(), // these values change the result
                None,
                None,
                None,
                None,
                external_conditions.into(),
                Some(external_sensor),
            )
            .unwrap(),
        );
        let control_min = create_setpoint_time_control(vec![
            Some(52.),
            None,
            None,
            None,
            Some(52.),
            Some(52.),
            Some(52.),
            Some(52.),
        ]);
        let control_max = create_setpoint_time_control(vec![
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
        ]);
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        let cold_water_source = ColdWaterSource::new(vec![1.0, 1.2], 0, 1.);
        let heat_battery_service: HeatBatteryServiceWaterRegular =
            HeatBatteryServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                WaterSourceWithTemperature::ColdWaterSource(Arc::new(cold_water_source)),
                Arc::new(control_min),
                Arc::new(control_max),
            );

        let result = heat_battery_service.energy_output_max(
            Default::default(),
            temp_return,
            simulation_time_iteration,
        );

        assert_relative_eq!(result, 5.637774816176471);
    }

    #[rstest]
    #[ignore = "while migrating to 0.34"]
    fn test_energy_output_max_service_off_for_water_regular(
        // In Python this is test_energy_output_max_service_off
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_off: Control,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off);
        let control_min = create_setpoint_time_control(vec![
            Some(52.),
            None,
            None,
            None,
            Some(52.),
            Some(52.),
            Some(52.),
            Some(52.),
        ]);
        let control_max = create_setpoint_time_control(vec![
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
        ]);

        let cold_water_source = ColdWaterSource::new(vec![1.0, 1.2], 0, 1.);
        let heat_battery_service: HeatBatteryServiceWaterRegular =
            HeatBatteryServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                WaterSourceWithTemperature::ColdWaterSource(Arc::new(cold_water_source)),
                Arc::new(control_min),
                Arc::new(control_max),
            );

        let temp_return = 40.;
        let result = heat_battery_service.energy_output_max(
            Default::default(),
            temp_return,
            simulation_time_iteration,
        );

        assert_eq!(result, 0.);
    }

    #[rstest]
    fn test_temp_setpnt_for_space(
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_off: Control,
    ) {
        let first_scheduled_temp = Some(21.);
        let ctrl: Control = create_setpoint_time_control(vec![first_scheduled_temp]);
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off);
        let heat_battery_space =
            HeatBatteryServiceSpace::new(heat_battery, SERVICE_NAME.into(), ctrl.into());

        assert_eq!(
            heat_battery_space.temp_setpnt(simulation_time_iteration),
            first_scheduled_temp
        );
    }

    #[rstest]
    fn test_in_required_period_for_space(
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_off: Control,
    ) {
        let ctrl: Control = create_setpoint_time_control(vec![Some(21.)]);
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off);
        let heat_battery_space =
            HeatBatteryServiceSpace::new(heat_battery, SERVICE_NAME.into(), ctrl.into());

        assert_eq!(
            heat_battery_space.in_required_period(simulation_time_iteration),
            Some(true)
        );
    }

    // TODO test_demand_energy_for_space

    #[rstest]
    #[ignore = "while migrating to 0.34"]
    fn test_demand_energy_service_off_for_space(
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_off: Control,
    ) {
        let energy_demand = 10.;
        let temp_return = 40.;
        let temp_flow = 1.;
        let ctrl: Control = create_setpoint_time_control(vec![None]);
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off);
        let heat_battery_space =
            HeatBatteryServiceSpace::new(heat_battery, SERVICE_NAME.into(), ctrl.into());
        let result = heat_battery_space
            .demand_energy(
                energy_demand,
                temp_return,
                temp_flow,
                None,
                None,
                simulation_time_iteration,
            )
            .unwrap();
        assert_eq!(result, 0.);
    }

    // in Python this test is called test_energy_output_max_service_on
    #[rstest]
    #[ignore = "while migrating to 0.34"]
    fn test_energy_output_max_service_on_for_space(
        battery_control_on: Control,
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let temp_output = 70.;
        let temp_return = 40.;
        let time_start = 0.1;
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        let service_control_on: Control =
            create_setpoint_time_control(vec![Some(21.0), Some(21.0)]);

        let heat_battery_service: HeatBatteryServiceSpace = HeatBatteryServiceSpace::new(
            heat_battery,
            SERVICE_NAME.into(),
            service_control_on.into(),
        );

        let result = heat_battery_service.energy_output_max(
            temp_output,
            temp_return,
            Some(time_start),
            simulation_time_iteration,
        );

        assert_relative_eq!(result, 4.037907837701614);
    }

    // in Python this test is called test_energy_output_max_service_off
    #[rstest]
    fn test_energy_output_max_service_off_for_space(
        battery_control_on: Control,
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let temp_output = 70.;
        let temp_return = 40.;
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        let service_control_off: Control = create_setpoint_time_control(vec![None]);

        let heat_battery_service: HeatBatteryServiceSpace = HeatBatteryServiceSpace::new(
            heat_battery,
            SERVICE_NAME.into(),
            service_control_off.into(),
        );

        let result = heat_battery_service.energy_output_max(
            temp_output,
            temp_return,
            None,
            simulation_time_iteration,
        );

        assert_relative_eq!(result, 0.);
    }

    #[rstest]
    fn test_create_service_connection(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        let create_connection_result =
            HeatBattery::create_service_connection(heat_battery.clone(), "new service");
        assert!(create_connection_result.is_ok());
        let create_connection_result =
            HeatBattery::create_service_connection(heat_battery, "new service");
        assert!(create_connection_result.is_err()) // second attempt to create a service connection with same name should error
    }

    #[rstest]
    fn test_electric_charge(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_off: Control,
        battery_control_on: Control,
    ) {
        // electric charge should be 0 when battery control is off
        let heat_battery =
            create_heat_battery(simulation_time_iterator.clone(), battery_control_off);
        assert_relative_eq!(heat_battery.lock().electric_charge(), 0.0);

        // electric charge should be calculated when battery control is on
        let heat_battery =
            create_heat_battery(simulation_time_iterator.clone(), battery_control_on);
        assert_relative_eq!(heat_battery.lock().electric_charge(), 20.0);
    }

    #[rstest]
    #[ignore = "while migrating to 0.34"]
    fn test_lab_test_rated_output(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_off: Control,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off);

        // assert_relative_eq!(heat_battery.lock().lab_test_rated_output(5.), 15.);
    }

    #[rstest]
    #[ignore = "while migrating to 0.34"]
    fn test_first_call(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
        simulation_time: SimulationTime,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        for (t_idx, _) in simulation_time.iter().enumerate() {
            heat_battery.lock().first_call();

            assert!(!heat_battery.lock().flag_first_call);
            // assert_relative_eq!(heat_battery.lock().q_in_ts.unwrap(), 20.);
            // assert_relative_eq!(heat_battery.lock().q_out_ts.unwrap(), 4.358566028225806);
            // assert_relative_eq!(heat_battery.lock().q_loss_ts.unwrap(), 0.031277812499999995);

            heat_battery.lock().timestep_end(t_idx).unwrap();
        }
    }

    #[rstest]
    #[ignore = "while migrating to 0.34"]
    fn test_demand_energy(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        simulation_time: SimulationTime,
        battery_control_on: Control,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);

        let service_name = "new_service";
        HeatBattery::create_service_connection(heat_battery.clone(), service_name).unwrap();

        for (t_idx, _) in simulation_time.iter().enumerate() {
            let demand_energy_actual = heat_battery
                .clone()
                .lock()
                .demand_energy(
                    service_name,
                    ServiceType::WaterRegular,
                    5.,
                    40.,
                    None,
                    false, // TODO 0.34 passed in to get it compiling, check if value is correct
                    None,
                    None,
                    t_idx,
                )
                .unwrap();

            assert_relative_eq!(demand_energy_actual, 4.358566028225806);

            let service_names_in_results = get_service_names_from_results(heat_battery.clone());

            assert!(service_names_in_results.contains(&service_name.into()));

            assert_relative_eq!(
                heat_battery.lock().charge_level,
                [0.19512695199092742, 0.2][t_idx]
            );

            assert_relative_eq!(heat_battery.lock().total_time_running_current_timestep, 1.);

            heat_battery.lock().timestep_end(t_idx).unwrap();
        }
    }

    #[rstest]
    /// Check heat battery auxilary energy consumption
    fn test_calc_auxiliary_energy(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
        let heat_battery =
            create_heat_battery(simulation_time_iterator.clone(), battery_control_on);

        heat_battery.lock().calc_auxiliary_energy(
            1.0,
            0.5,
            simulation_time_iterator.current_index(),
        );

        let results_by_end_user = heat_battery
            .lock()
            .energy_supply
            .read()
            .results_by_end_user();

        let end_user_name = heat_battery
            .lock()
            .energy_supply_connection
            .end_user_name
            .clone();

        let results_by_end_user = results_by_end_user.get(&end_user_name).unwrap();

        assert_eq!(*results_by_end_user, vec![0.0122, 0.]);
    }

    #[rstest]
    #[ignore = "while migrating to 0.34"]
    fn test_timestep_end(
        external_sensor: ExternalSensor,
        external_conditions: ExternalConditions,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        // not using the fixture here
        // because we need to set different charge_levels
        let battery_control_on: Control = Control::Charge(
            ChargeControl::new(
                ControlLogicType::Manual,
                vec![true, true, true],
                1.,
                0,
                1.,
                [1.0, 1.5].into_iter().map(Into::into).collect(),
                None,
                None,
                None,
                None,
                external_conditions.into(),
                Some(external_sensor),
            )
            .unwrap(),
        );

        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        let service_name = "new_timestep_end_service";
        HeatBattery::create_service_connection(heat_battery.clone(), service_name).unwrap();

        let t_idx = 0;
        heat_battery
            .lock()
            .demand_energy(
                service_name,
                ServiceType::WaterRegular,
                5.0,
                40.,
                None,
                false, // TODO 0.34 passed in to get it compiling, check if value is correct
                None,
                None,
                t_idx,
            )
            .unwrap();

        // assert_relative_eq!(heat_battery.lock().q_in_ts.unwrap(), 20.);
        // assert_relative_eq!(heat_battery.lock().q_out_ts.unwrap(), 5.637774816176471);
        // assert_relative_eq!(heat_battery.lock().q_loss_ts.unwrap(), 0.03929547794117647);
        assert_relative_eq!(
            heat_battery.lock().total_time_running_current_timestep,
            0.8868747268254661
        );

        let service_names_in_results = get_service_names_from_results(heat_battery.clone());

        assert!(service_names_in_results.contains(&service_name.into()));

        heat_battery.lock().timestep_end(t_idx).unwrap();

        assert!(!heat_battery.lock().flag_first_call);
        // assert_relative_eq!(heat_battery.lock().q_in_ts.unwrap(), 20.);
        // assert_relative_eq!(heat_battery.lock().q_out_ts.unwrap(), 10.001923317091928);
        // assert_relative_eq!(heat_battery.lock().q_loss_ts.unwrap(), 0.07624000227068732);
        assert_relative_eq!(heat_battery.lock().total_time_running_current_timestep, 0.0);
        assert_eq!(heat_battery.lock().service_results.len(), 0);
    }

    #[rstest]
    #[ignore = "while migrating to 0.34"]
    fn test_energy_output_max(
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        simulation_time: SimulationTime,
    ) {
        // not using the fixture here
        // because we need to set different charge_levels
        let battery_control_on: Control = Control::Charge(
            ChargeControl::new(
                ControlLogicType::Manual,
                vec![true, true, true],
                1.,
                0,
                1.,
                [1.5, 1.6].into_iter().map(Into::into).collect(), // these values change the result
                None,
                None,
                None,
                None,
                external_conditions.into(),
                Some(external_sensor),
            )
            .unwrap(),
        );

        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);

        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                heat_battery.lock().energy_output_max(0., None),
                [5.637774816176471, 11.13482970854502][t_idx]
            );

            heat_battery.lock().timestep_end(t_idx).unwrap();
        }
    }
}
