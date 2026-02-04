use crate::compare_floats::min_of_2;
/// This module provides object(s) to model the behaviour of heat batteries.
use crate::core::common::{WaterSupply, WaterSupplyBehaviour};
use crate::core::controls::time_control::{per_control, Control, ControlBehaviour};
use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::core::heating_systems::common::HeatingServiceType;
use crate::core::material_properties::WATER;
use crate::core::units::{
    KILOJOULES_PER_KILOWATT_HOUR, MILLIMETRES_IN_METRE, SECONDS_PER_HOUR, SECONDS_PER_MINUTE,
    WATTS_PER_KILOWATT,
};
use crate::core::water_heat_demand::misc::{
    calculate_volume_weighted_average_temperature, water_demand_to_kwh, WaterEventResult,
};
use crate::corpus::{ResultParamValue, ResultsAnnual, ResultsPerTimestep};
use crate::input::{HeatBattery as HeatBatteryInput, HeatSourceWetDetails};
use crate::simulation_time::{SimulationTimeIteration, SimulationTimeIterator};
use anyhow::bail;
use atomic_float::AtomicF64;
use fsum::FSum;
use indexmap::IndexMap;
use itertools::Itertools;
use parking_lot::RwLock;
use smartstring::alias::String;
use std::ops::Deref;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use thiserror::Error;

#[derive(Clone, Copy, Debug)]
pub(crate) enum HeatBatteryPcmOperationMode {
    Normal,
    OnlyCharging,
    Losses,
}

/// An object to represent a water heating service provided by a regular heat battery.
///
/// This object contains the parts of the heat battery calculation that are
/// specific to providing hot water.
#[derive(Clone, Debug)]
pub struct HeatBatteryPcmServiceWaterRegular {
    heat_battery: Arc<RwLock<HeatBatteryPcm>>,
    service_name: String,
    cold_feed: WaterSupply,
    control: Arc<Control>,
    control_min: Arc<Control>,
    control_max: Arc<Control>,
}

impl HeatBatteryPcmServiceWaterRegular {
    /// Arguments:
    /// * `heat_battery` - reference to the Heat Battery object providing the service
    /// * `service_name` - name of the service demanding energy
    /// * `cold_feed` - reference to ColdWaterSource object
    /// * `control_min` - reference to a control object which must select current the minimum timestep temperature
    /// * `control_max` - reference to a control object which must select current the maximum timestep temperature
    pub(crate) fn new(
        heat_battery: Arc<RwLock<HeatBatteryPcm>>,
        service_name: String,
        cold_feed: WaterSupply,
        control_min: Arc<Control>,
        control_max: Arc<Control>,
    ) -> Self {
        let control = control_min.clone();

        Self {
            heat_battery,
            service_name,
            cold_feed,
            control,
            control_min,
            control_max,
        }
    }

    /// Return setpoint (not necessarily temperature)
    pub(crate) fn setpnt(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (Option<f64>, Option<f64>) {
        (
            self.control_min.setpnt(&simulation_time_iteration),
            self.control_max.setpnt(&simulation_time_iteration),
        )
    }

    /// Demand energy (in kWh) from the heat_battery
    pub(crate) fn demand_energy(
        &self,
        energy_demand: f64,
        temp_flow: Option<f64>,
        temp_return: f64,
        update_heat_source_state: Option<bool>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let service_on = self.is_on(simulation_time_iteration);
        let energy_demand = if !service_on { 0.0 } else { energy_demand };
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);

        self.heat_battery.read().demand_energy(
            &self.service_name,
            HeatingServiceType::DomesticHotWaterRegular,
            energy_demand,
            temp_return,
            temp_flow,
            service_on,
            None,
            Some(update_heat_source_state),
        )
    }

    pub(crate) fn energy_output_max(
        &self,
        temp_flow: f64,
        _temp_return: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let service_on = self.is_on(simulation_time_iteration);
        if !service_on {
            return Ok(0.);
        }

        self.heat_battery.read().energy_output_max(temp_flow, None)
    }

    fn is_on(&self, simtime: SimulationTimeIteration) -> bool {
        self.control.is_on(&simtime)
    }
}

/// An object to represent a direct water heating service provided by a heat battery.
///
/// This is similar to a combi boiler or HIU providing hot water on demand.
#[derive(Debug, Clone)]
pub(crate) struct HeatBatteryPcmServiceWaterDirect {
    heat_battery: Arc<RwLock<HeatBatteryPcm>>,
    service_name: String,
    setpoint_temp: f64,
    cold_feed: WaterSupply,
}

impl HeatBatteryPcmServiceWaterDirect {
    /// Arguments:
    /// * `heat_battery` - reference to the HeatBatteryPCM object providing the service
    /// * `service_name` - name of the service demanding energy from the heat battery
    /// * `setpoint_temp` - temperature of hot water to be provided, in deg C
    /// * `cold_feed` - reference to ColdWaterSource object
    fn new(
        heat_battery: Arc<RwLock<HeatBatteryPcm>>,
        service_name: String,
        setpoint_temp: f64,
        cold_feed: WaterSupply,
    ) -> Self {
        Self {
            heat_battery,
            service_name,
            setpoint_temp,
            cold_feed,
        }
    }

    pub(crate) fn get_cold_water_source(&self) -> &WaterSupply {
        &self.cold_feed
    }

    fn temp_hot_water(
        &self,
        vol: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let list_temp_vol = self
            .cold_feed
            .get_temp_cold_water(vol, simulation_time_iteration)?;
        let inlet_temp =
            calculate_volume_weighted_average_temperature(list_temp_vol, Some(vol), None)?;

        self.heat_battery
            .read()
            .get_temp_hot_water(inlet_temp, vol, self.setpoint_temp)
    }

    pub(crate) fn get_temp_hot_water(
        &self,
        volume_req: f64,
        volume_req_already: Option<f64>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        let volume_req_already = volume_req_already.unwrap_or(0.);

        if is_close!(volume_req, 0., rel_tol = 1e-09, abs_tol = 1e-10) {
            return Ok(vec![]);
        }

        let volume_req_cumulative = volume_req + volume_req_already;
        let temp_hot_water_cumulative =
            self.temp_hot_water(volume_req_cumulative, simulation_time_iteration)?;

        // Base temperature on the part of the draw-off for volume_req, and
        // ignore any volume previously considered
        let temp_hot_water_req =
            if is_close!(volume_req_already, 0., rel_tol = 1e-09, abs_tol = 1e-10) {
                temp_hot_water_cumulative
            } else {
                let temp_hot_water_req_already =
                    self.temp_hot_water(volume_req_already, simulation_time_iteration)?;

                (temp_hot_water_cumulative * volume_req_cumulative
                    - temp_hot_water_req_already * volume_req_already)
                    / volume_req
            };

        Ok(vec![(temp_hot_water_req, volume_req)])
    }

    /// Process hot water demand directly from dry core heat battery
    pub(crate) fn demand_hot_water(
        &self,
        usage_events: Option<Vec<WaterEventResult>>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let mut energy_demand = 0.;
        let mut total_volume = 0.;
        let mut weighted_cold_temp_sum = 0.;

        if let Some(events) = usage_events {
            for event in events {
                if is_close!(event.volume_hot, 0., rel_tol = 1e-09, abs_tol = 1e-10) {
                    continue;
                }
                // Skip this event if no temperature available
                if let Some(hot_temp) = self
                    .get_temp_hot_water(event.volume_hot, None, simtime)?
                    .first()
                    .map(|(t, _v)| t)
                {
                    let list_temp_vol = self.cold_feed.draw_off_water(event.volume_hot, simtime)?;
                    let cold_temp = calculate_volume_weighted_average_temperature(
                        list_temp_vol,
                        Some(event.volume_hot), // This validates the volume
                        None,
                    )?;

                    // Calculate energy needed to heat water
                    energy_demand += water_demand_to_kwh(event.volume_hot, *hot_temp, cold_temp);

                    // Accumulate for weighted average cold water temperature
                    total_volume += event.volume_hot;
                    weighted_cold_temp_sum += cold_temp * event.volume_hot;
                }
            }
        }

        // Calculate weighted average cold water temperature
        let cold_water_temp = if total_volume > 0. {
            weighted_cold_temp_sum / total_volume
        } else {
            // Fallback to sampling method if no events processed
            let cold_water_temp_vol = self.cold_feed.get_temp_cold_water(1., simtime)?;

            calculate_volume_weighted_average_temperature(cold_water_temp_vol, Some(1.), None)?
        };

        // Demand energy from heat battery
        self.heat_battery.read().demand_energy(
            self.service_name.as_str(),
            HeatingServiceType::DomesticHotWaterDirect,
            energy_demand,
            cold_water_temp, // return temperature (cold water inlet)
            None,            // flow temperature (hot water outlet)
            true,
            None,
            Some(true),
        )
    }
}

#[derive(Clone, Debug)]
pub struct HeatBatteryPcmServiceSpace {
    heat_battery: Arc<RwLock<HeatBatteryPcm>>,
    service_name: String,
    control: Arc<Control>,
}

/// An object to represent a space heating service provided by a heat_battery to e.g. radiators.
///
/// This object contains the parts of the heat battery calculation that are
/// specific to providing space heating.
impl HeatBatteryPcmServiceSpace {
    pub(crate) fn new(
        heat_battery: Arc<RwLock<HeatBatteryPcm>>,
        service_name: String,
        control: Arc<Control>, // in Python this is SetpointTimeControl | CombinationTimeControl
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

        self.heat_battery.read().demand_energy(
            &self.service_name,
            HeatingServiceType::Space,
            energy_demand,
            temp_return,
            Some(temp_flow),
            service_on,
            None,
            Some(update_heat_source_state),
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
    ) -> anyhow::Result<f64> {
        let time_start = time_start.unwrap_or(0.);

        if !self.is_on(simtime) {
            return Ok(0.);
        }

        self.heat_battery
            .read()
            .energy_output_max(temp_output, Some(time_start))
    }
}

const DEFAULT_N_LAYERS: usize = 8; // Number of calculation layers in heat battery
const DEFAULT_TIME_STEP_SECONDS: f64 = 20.; // Time step for iterative calculations (seconds)
const DEFAULT_INLET_TEMP_CELSIUS: f64 = 10.; // Initial inlet temperature for Reynolds number calculation (°C)
const DEFAULT_OUTLET_TEMP_CELSIUS: f64 = 53.; // Estimated outlet temperature for Reynolds number calculation (°C)

// nothing seems to read this - check upstream whether service_results field is necessary
#[derive(Clone, Debug)]
#[allow(dead_code)]
struct HeatBatteryResult {
    service_name: String,
    service_type: Option<HeatingServiceType>,
    service_on: bool,
    energy_output_required: f64,
    temp_output: Option<f64>,
    temp_inlet: f64,
    time_running: f64,
    energy_delivered_hb: f64,
    energy_delivered_backup: f64,
    energy_delivered_total: f64,
    energy_charged_during_service: f64,
    hb_zone_temperatures: Vec<f64>,
    current_hb_power: f64,
}

impl HeatBatteryResult {
    fn param(&self, param: &str) -> ResultParamValue {
        match param {
            "service_name" => ResultParamValue::from(self.service_name.clone()),
            "service_type" => self
                .service_type
                .as_ref()
                .map(|service_type| ResultParamValue::from(String::from(service_type.to_string())))
                .unwrap_or(ResultParamValue::Empty),
            "service_on" => self.service_on.into(),
            "energy_output_required" => self.energy_output_required.into(),
            "temp_output" => self.temp_output.into(),
            "temp_inlet" => self.temp_inlet.into(),
            "time_running" => self.time_running.into(),
            "energy_delivered_HB" => self.energy_delivered_hb.into(),
            "energy_delivered_backup" => self.energy_delivered_backup.into(),
            "energy_delivered_total" => self.energy_delivered_total.into(),
            "energy_charged_during_service" => self.energy_charged_during_service.into(),
            "current_hb_power" => self.current_hb_power.into(),
            _ => panic!("Unknown parameter: {}", param),
        }
    }
}

const OUTPUT_PARAMETERS: [(&str, Option<&str>, bool); 13] = [
    ("service_name", None, false),
    ("service_type", None, false),
    ("service_on", None, false),
    ("energy_output_required", Some("kWh"), true),
    ("temp_output", Some("degC"), false),
    ("temp_inlet", Some("degC"), false),
    ("time_running", Some("secs"), true),
    ("energy_delivered_HB", Some("kWh"), true),
    ("energy_delivered_backup", Some("kWh"), true),
    ("energy_delivered_total", Some("kWh"), true),
    ("energy_charged_during_service", Some("kWh"), true),
    ("hb_zone_temperatures", Some("degC"), false),
    ("current_hb_power", Some("kW"), false),
];
const AUX_PARAMETERS: [(&str, Option<&str>, bool); 6] = [
    ("energy_aux", Some("kWh"), true),
    ("battery_losses", Some("kWh"), true),
    ("Temps_after_losses", Some("degC"), false),
    ("total_charge", Some("kWh"), true),
    ("end_of_timestep_charge", Some("kWh"), true),
    ("hb_after_only_charge_zone_temp", Some("degC"), false),
];

#[derive(Debug)]
struct HeatBatteryTimestepSummary {
    energy_aux: f64,
    battery_losses: f64,
    temps_after_losses: Vec<f64>,
    total_charge: f64,
    end_of_timestep_charge: f64,
    hb_after_only_charge_zone_temp: Vec<f64>,
}

impl HeatBatteryTimestepSummary {
    fn param(&self, param: &str) -> ResultParamValue {
        match param {
            "energy_aux" => self.energy_aux.into(),
            "battery_losses" => self.battery_losses.into(),
            "total_charge" => self.total_charge.into(),
            "end_of_timestep_charge" => self.end_of_timestep_charge.into(),
            _ => panic!("Parameter {param} not recognised"),
        }
    }
}

#[derive(Debug)]
struct HeatBatteryTimestepResult {
    results: Vec<HeatBatteryResult>,
    summary: HeatBatteryTimestepSummary,
}

#[derive(Clone, Debug)]
struct PipeEnergy {
    energy: f64,
    temperature: f64,
}

#[derive(Debug)]
pub struct HeatBatteryPcm {
    simulation_time: Arc<SimulationTimeIterator>,
    energy_supply: Arc<RwLock<EnergySupply>>,
    energy_supply_connection: EnergySupplyConnection,
    energy_supply_connections: IndexMap<String, EnergySupplyConnection>,
    pwr_in: f64,
    max_rated_losses: f64,
    power_circ_pump: f64,
    power_standby: f64,
    n_units: usize,
    charge_control: Arc<Control>, // ChargeControl variant expected
    // nothing external seems to read this - check upstream whether service_results field is necessary
    service_results: Arc<RwLock<Vec<HeatBatteryResult>>>,
    total_time_running_current_timestep: AtomicF64,
    pump_running_time_current_timestep: AtomicF64,
    flag_first_call: AtomicBool,
    #[allow(dead_code)]
    charge_level: f64,
    battery_losses: AtomicF64,
    n_layers: usize,
    hb_time_step: f64,
    initial_inlet_temp: f64,
    estimated_outlet_temp: f64,
    energy_charged: AtomicF64,
    simultaneous_charging_and_discharging: bool,
    max_temp_of_charge: f64,
    zone_temp_c_dist_initial: Arc<RwLock<Vec<f64>>>,
    heat_storage_kj_per_k_above: f64,
    heat_storage_kj_per_k_below: f64,
    heat_storage_kj_per_k_during: f64,
    phase_transition_temperature_upper: f64,
    phase_transition_temperature_lower: f64,
    velocity_in_hex_tube: f64,
    capillary_diameter_m: f64,
    a: f64,
    b: f64,
    flow_rate_l_per_min: f64,
    detailed_results: Option<Arc<RwLock<Vec<HeatBatteryTimestepResult>>>>,
}

impl HeatBatteryPcm {
    pub(crate) fn new(
        heat_battery_details: &HeatSourceWetDetails,
        charge_control: Arc<Control>,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_connection: EnergySupplyConnection,
        simulation_time: Arc<SimulationTimeIterator>,
        n_layers: Option<usize>,
        hb_time_step: Option<f64>,
        initial_inlet_temp: Option<f64>,
        estimated_outlet_temp: Option<f64>,
        output_detailed_results: Option<bool>,
    ) -> Self {
        let n_layers = n_layers.unwrap_or(DEFAULT_N_LAYERS);
        let hb_time_step = hb_time_step.unwrap_or(DEFAULT_TIME_STEP_SECONDS);

        let initial_inlet_temp = initial_inlet_temp.unwrap_or(DEFAULT_INLET_TEMP_CELSIUS);
        let estimated_outlet_temp = estimated_outlet_temp.unwrap_or(DEFAULT_OUTLET_TEMP_CELSIUS);
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
            inlet_diameter_mm,
            a,
            b,
            flow_rate_l_per_min,
            temp_init,
            ..,
        ) = if let HeatSourceWetDetails::HeatBattery {
            battery:
                HeatBatteryInput::Pcm {
                    rated_charge_power: pwr_in,
                    max_rated_losses,
                    electricity_circ_pump: power_circ_pump,
                    electricity_standby: power_standby,
                    number_of_units: n_units,
                    simultaneous_charging_and_discharging,
                    max_temperature,
                    heat_storage_k_j_per_k_above_phase_transition: heat_storage_kj_per_k_above,
                    heat_storage_k_j_per_k_below_phase_transition: heat_storage_kj_per_k_below,
                    heat_storage_k_j_per_k_during_phase_transition: heat_storage_kj_per_k_during,
                    phase_transition_temperature_upper,
                    phase_transition_temperature_lower,
                    velocity_in_hex_tube_at_1_l_per_min_m_per_s: velocity_in_hex_tube,
                    inlet_diameter_mm,
                    a,
                    b,
                    flow_rate_l_per_min,
                    temp_init,
                    ..
                },
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
                *inlet_diameter_mm,
                *a,
                *b,
                *flow_rate_l_per_min,
                *temp_init,
            )
        } else {
            unreachable!()
        };

        let detailed_results = if output_detailed_results {
            Some(Default::default())
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
            pump_running_time_current_timestep: Default::default(),
            flag_first_call: true.into(),
            charge_level: Default::default(),
            battery_losses: Default::default(),
            n_layers,
            hb_time_step,
            initial_inlet_temp,
            estimated_outlet_temp,
            energy_charged: Default::default(),
            simultaneous_charging_and_discharging,
            max_temp_of_charge,
            zone_temp_c_dist_initial: Arc::new(RwLock::new(vec![temp_init; n_layers])),
            heat_storage_kj_per_k_above: heat_storage_kj_per_k_above / n_layers as f64,
            heat_storage_kj_per_k_below: heat_storage_kj_per_k_below / n_layers as f64,
            heat_storage_kj_per_k_during: heat_storage_kj_per_k_during / n_layers as f64,
            phase_transition_temperature_upper,
            phase_transition_temperature_lower,
            velocity_in_hex_tube,
            capillary_diameter_m: inlet_diameter_mm / MILLIMETRES_IN_METRE as f64,
            a,
            b,
            flow_rate_l_per_min,
            detailed_results,
        }
    }

    fn create_service_connection(
        heat_battery: Arc<RwLock<Self>>,
        service_name: &str,
    ) -> anyhow::Result<()> {
        if heat_battery
            .read()
            .energy_supply_connections
            .contains_key(service_name)
        {
            bail!("Error: Service name already used: {service_name}");
        }
        let energy_supply = heat_battery.read().energy_supply.clone();

        // Set up EnergySupplyConnection for this service
        heat_battery.write().energy_supply_connections.insert(
            service_name.into(),
            EnergySupply::connection(energy_supply, service_name)?,
        );

        Ok(())
    }

    /// Return a HeatBatteryPcmServiceWaterRegular object and create an EnergySupplyConnection for it
    ///
    /// Arguments:
    /// * `heat_battery` - reference to heat battery
    /// * `service_name` - name of the service demanding energy from the heat battery
    /// * `cold_feed` - reference to ColdWaterSource object
    /// * `control_min` - reference to a control object which must select current the minimum timestep temperature
    /// * `control_max` - reference to a control object which must select current the maximum timestep temperature
    pub(crate) fn create_service_hot_water_regular(
        heat_battery: Arc<RwLock<Self>>,
        service_name: &str,
        cold_feed: WaterSupply,
        control_min: Arc<Control>,
        control_max: Arc<Control>,
    ) -> anyhow::Result<HeatBatteryPcmServiceWaterRegular> {
        Self::create_service_connection(heat_battery.clone(), service_name)?;
        Ok(HeatBatteryPcmServiceWaterRegular::new(
            heat_battery,
            service_name.into(),
            cold_feed,
            control_min,
            control_max,
        ))
    }

    /// Return a HeatBatteryPCMServiceWaterDirect object and create an EnergySupplyConnection for it
    ///
    /// Arguments:
    /// * `heat_battery` - reference to heat battery
    /// * `service_name` - name of the service demanding energy from the heat battery
    /// * `setpoint_temp` - temperature of hot water to be provided, in deg C
    /// * `cold_feed` - reference to ColdWaterSource object
    pub(crate) fn create_service_hot_water_direct(
        heat_battery: Arc<RwLock<Self>>,
        service_name: &str,
        setpoint_temp: f64,
        cold_feed: WaterSupply,
    ) -> anyhow::Result<HeatBatteryPcmServiceWaterDirect> {
        Self::create_service_connection(heat_battery.clone(), service_name)?;
        Ok(HeatBatteryPcmServiceWaterDirect::new(
            heat_battery,
            service_name.into(),
            setpoint_temp,
            cold_feed,
        ))
    }

    /// Return a HeatBatteryPCMServiceSpace object and create an EnergySupplyConnection for it
    ///
    /// Arguments:
    /// * `heat_battery` - reference to heat battery
    /// * `service_name` - name of the service demanding energy from the heat battery
    /// * `control` - reference to a control object which must implement is_on() and setpnt() funcs
    pub(crate) fn create_service_space_heating(
        heat_battery: Arc<RwLock<Self>>,
        service_name: &str,
        control: Arc<Control>, // in Python this is SetpointTimeControl | CombinationTimeControl
    ) -> anyhow::Result<HeatBatteryPcmServiceSpace> {
        Self::create_service_connection(heat_battery.clone(), service_name)?;
        Ok(HeatBatteryPcmServiceSpace::new(
            heat_battery,
            service_name.into(),
            control,
        ))
    }

    /// Return battery losses
    pub(crate) fn get_battery_losses(&self) -> f64 {
        let battery_losses = self.battery_losses.load(Ordering::SeqCst) * self.n_units as f64;
        self.battery_losses.store(0., Ordering::SeqCst);
        battery_losses
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
        (timestep
            - self
                .total_time_running_current_timestep
                .load(Ordering::SeqCst))
            * (1. - time_start / timestep)
    }

    fn calculate_heat_transfer_kw_per_k(
        a: f64,
        b: f64,
        flow_rate_l_per_min: f64,
        reynold_number_at_1_l_per_min: f64,
    ) -> f64 {
        (a * (reynold_number_at_1_l_per_min * flow_rate_l_per_min).ln() + b)
            / WATTS_PER_KILOWATT as f64
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

        0.000000000145238 * average_temp.powi(2) - 0.0000000248238 * average_temp + 0.000001432
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
        mode: &HeatBatteryPcmOperationMode,
        zone_temp_c_dist: &[f64],
        inlet_temp_c: f64,
        inlet_temp_c_zone: f64,
        q_max_kj: f64,
        reynold_number_at_1_l_per_min: f64,
        flow_rate_kg_per_s: f64,
        time_step_s: f64,
    ) -> (f64, usize, f64, f64) {
        match mode {
            HeatBatteryPcmOperationMode::OnlyCharging => {
                let zone_index = zone_temp_c_dist.iter().len() - index - 1;
                let zone_temp_c_start = zone_temp_c_dist[zone_index];
                (0., zone_index, zone_temp_c_start, 0.)
            }
            HeatBatteryPcmOperationMode::Losses => {
                let zone_temp_c_start = zone_temp_c_dist[index];
                let energy_transf = if zone_temp_c_start > inlet_temp_c {
                    q_max_kj / zone_temp_c_dist.len() as f64
                } else {
                    0.
                };
                (energy_transf, index, zone_temp_c_start, 0.)
            }
            HeatBatteryPcmOperationMode::Normal => {
                // NORMAL mode include battery primarily hydraulic charging or discharging with or without simultaneous electric charging
                let zone_temp_c_start = zone_temp_c_dist[index];
                let heat_transfer_kw_per_k = Self::calculate_heat_transfer_kw_per_k(
                    self.a,
                    self.b,
                    self.flow_rate_l_per_min,
                    reynold_number_at_1_l_per_min,
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
        } else if target_temp > self.phase_transition_temperature_upper {
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
                    if energy_transf >= q_required {
                        // There is plenty of charging after taking zone temperature to target
                        q_max_kj -= q_required - energy_transf;
                        energy_charged +=
                            -(q_required - energy_transf) / KILOJOULES_PER_KILOWATT_HOUR as f64;
                        energy_transf = q_required;
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

    /// ranges _1, _2, and _3 refer to:
    /// _1: temperature of PCM above transition phase
    /// _2: temperature of PCM within transition phase
    /// _3: temperature of PCM below transition phase
    fn calculate_new_zone_temperature(
        &self,
        zone_temp_c_start: f64,
        mut energy_transf: f64,
    ) -> f64 {
        let mut delta_temp_1 = 0.;
        let mut delta_temp_2 = 0.;
        let mut delta_temp_3 = 0.;

        if energy_transf > 0. {
            // zone delivering energy to water
            if zone_temp_c_start >= self.phase_transition_temperature_upper {
                let heat_range_1 = (zone_temp_c_start - self.phase_transition_temperature_upper)
                    * self.heat_storage_kj_per_k_above;
                let heat_range_2 = (self.phase_transition_temperature_upper
                    - self.phase_transition_temperature_lower)
                    * self.heat_storage_kj_per_k_during;

                if energy_transf <= heat_range_1 {
                    delta_temp_1 = energy_transf / self.heat_storage_kj_per_k_above;
                } else {
                    delta_temp_1 = zone_temp_c_start - self.phase_transition_temperature_upper;

                    energy_transf -= heat_range_1;
                    if energy_transf <= heat_range_2 {
                        delta_temp_2 = energy_transf / self.heat_storage_kj_per_k_during;
                    } else {
                        delta_temp_2 = self.phase_transition_temperature_upper
                            - self.phase_transition_temperature_lower;
                        energy_transf -= heat_range_2;
                        delta_temp_3 = energy_transf / self.heat_storage_kj_per_k_below;
                    }
                }
            } else if self.phase_transition_temperature_lower <= zone_temp_c_start
                && zone_temp_c_start < self.phase_transition_temperature_upper
            {
                let heat_range_2 = (zone_temp_c_start - self.phase_transition_temperature_lower)
                    * self.heat_storage_kj_per_k_during;

                if energy_transf <= heat_range_2 {
                    delta_temp_2 = energy_transf / self.heat_storage_kj_per_k_during;
                } else {
                    delta_temp_2 = zone_temp_c_start - self.phase_transition_temperature_lower;
                    energy_transf -= heat_range_2;
                    delta_temp_3 = energy_transf / self.heat_storage_kj_per_k_below
                }
            } else {
                delta_temp_3 = energy_transf / self.heat_storage_kj_per_k_below;
            }
        } else if energy_transf < 0. {
            // zone retrieving energy from water
            if zone_temp_c_start <= self.phase_transition_temperature_lower {
                let heat_range_3 = (zone_temp_c_start - self.phase_transition_temperature_lower)
                    * self.heat_storage_kj_per_k_below;
                let heat_range_2 = (self.phase_transition_temperature_lower
                    - self.phase_transition_temperature_upper)
                    * self.heat_storage_kj_per_k_during;

                if energy_transf >= heat_range_3 {
                    delta_temp_3 = energy_transf / self.heat_storage_kj_per_k_below;
                } else {
                    delta_temp_3 = zone_temp_c_start - self.phase_transition_temperature_lower;

                    energy_transf -= heat_range_3;
                    if energy_transf >= heat_range_2 {
                        delta_temp_2 = energy_transf / self.heat_storage_kj_per_k_during;
                    } else {
                        delta_temp_2 = self.phase_transition_temperature_lower
                            - self.phase_transition_temperature_upper;

                        energy_transf -= heat_range_2;
                        delta_temp_1 = energy_transf / self.heat_storage_kj_per_k_above;
                    }
                }
            } else if self.phase_transition_temperature_lower < zone_temp_c_start
                && zone_temp_c_start <= self.phase_transition_temperature_upper
            {
                let heat_range_2 = (zone_temp_c_start - self.phase_transition_temperature_upper)
                    * self.heat_storage_kj_per_k_during;

                if energy_transf >= heat_range_2 {
                    delta_temp_2 = energy_transf / self.heat_storage_kj_per_k_during;
                } else {
                    delta_temp_2 = zone_temp_c_start - self.phase_transition_temperature_upper;
                    energy_transf -= heat_range_2;
                    delta_temp_1 = energy_transf / self.heat_storage_kj_per_k_above;
                }
            } else {
                delta_temp_1 = energy_transf / self.heat_storage_kj_per_k_above;
            }
        }
        zone_temp_c_start - (delta_temp_1 + delta_temp_2 + delta_temp_3)
    }

    fn process_heat_battery_zones(
        &self,
        inlet_temp_c: f64,
        zone_temp_c_dist: &mut [f64],
        flow_rate_kg_per_s: f64,
        time_step_s: f64,
        reynold_number_at_1_l_per_min: f64,
        pwr_in: Option<f64>,
        mode: Option<HeatBatteryPcmOperationMode>,
    ) -> anyhow::Result<(f64, Vec<f64>, Vec<f64>, f64)> {
        let pwr_in = pwr_in.unwrap_or(0.);
        let mode = mode.unwrap_or(HeatBatteryPcmOperationMode::Normal);
        let target_temp = self.max_temp_of_charge * self.target_charge()?;
        let mut energy_charged = 0.;

        let mut q_max_kj =
            -pwr_in * time_step_s / SECONDS_PER_HOUR as f64 * KILOJOULES_PER_KILOWATT_HOUR as f64;

        let mut energy_transf_delivered = vec![0.; self.n_layers];
        let mut inlet_temp_c_zone = inlet_temp_c;
        let mut energy_transf;
        let mut zone_index;
        let mut zone_temp_c_start;
        let mut outlet_temp_c = Default::default();

        for index in 0..zone_temp_c_dist.iter().len() {
            // Get zone index, starting temperature, outlet temperature and energy_transfer based on operation mode
            (energy_transf, zone_index, zone_temp_c_start, outlet_temp_c) = self
                .get_zone_properties(
                    index,
                    &mode,
                    zone_temp_c_dist,
                    inlet_temp_c,
                    inlet_temp_c_zone,
                    q_max_kj,
                    reynold_number_at_1_l_per_min,
                    flow_rate_kg_per_s,
                    time_step_s,
                );

            energy_transf_delivered[zone_index] += energy_transf;

            // Process energy transfer in zone with simultaneous charging
            if q_max_kj < 0. {
                (q_max_kj, energy_charged, energy_transf) = self
                    .process_zone_simultaneous_charging(
                        zone_temp_c_start,
                        target_temp,
                        q_max_kj,
                        energy_transf,
                        energy_charged,
                    );
            };

            // Recalculate zone temperatures after energy transfer
            zone_temp_c_dist[zone_index] =
                self.calculate_new_zone_temperature(zone_temp_c_start, energy_transf);

            // Update values for the next iteration
            inlet_temp_c_zone = outlet_temp_c
        }

        Ok((
            outlet_temp_c,
            zone_temp_c_dist.to_owned(),
            energy_transf_delivered,
            energy_charged,
        ))
    }

    /// Charge the battery (update the zones temperature)
    /// It follows the same methodology as energy_demand function
    fn _charge_battery_hydraulic(&mut self, inlet_temp_c: f64) -> anyhow::Result<f64> {
        let total_time_s = self.simulation_time.step_in_hours() * SECONDS_PER_HOUR as f64;
        let time_step_s = self.hb_time_step;

        // Initial Reynold number
        let mut water_kinematic_viscosity_m2_per_s =
            Self::calculate_water_kinematic_viscosity_m2_per_s(
                self.initial_inlet_temp,
                self.estimated_outlet_temp,
            );
        let mut reynold_number_at_1_l_per_min = Self::calculate_reynold_number_at_1_l_per_min(
            water_kinematic_viscosity_m2_per_s,
            self.velocity_in_hex_tube,
            self.capillary_diameter_m,
        );

        let flow_rate_kg_per_s =
            (self.flow_rate_l_per_min / SECONDS_PER_MINUTE as f64) * WATER.density();
        let n_time_steps = (total_time_s / time_step_s) as usize;
        let mut zone_temp_c_dist = self.zone_temp_c_dist_initial.read().clone();
        let mut total_charge = 0.;

        for _ in 0..n_time_steps {
            // Processing HB zones
            let (outlet_temp_c, zone_temp_c_dist_new, energy_transf_charged, _) = self
                .process_heat_battery_zones(
                    inlet_temp_c,
                    &mut zone_temp_c_dist,
                    flow_rate_kg_per_s,
                    time_step_s,
                    reynold_number_at_1_l_per_min,
                    Some(0.),
                    None,
                )?;

            zone_temp_c_dist = zone_temp_c_dist_new;

            // RN for next time step
            water_kinematic_viscosity_m2_per_s =
                Self::calculate_water_kinematic_viscosity_m2_per_s(inlet_temp_c, outlet_temp_c);
            reynold_number_at_1_l_per_min = Self::calculate_reynold_number_at_1_l_per_min(
                water_kinematic_viscosity_m2_per_s,
                self.velocity_in_hex_tube,
                self.capillary_diameter_m,
            );

            // Equivalent of using Python's math.fsum instead of sum() for better numerical accuracy with floating point arithmetic
            let energy_charged_during_battery_time_step =
                FSum::with_all(&energy_transf_charged).value();

            if outlet_temp_c < inlet_temp_c {
                total_charge += energy_charged_during_battery_time_step;
            } else {
                break;
            }
        }

        *self.zone_temp_c_dist_initial.write() = zone_temp_c_dist;

        Ok(total_charge)
    }

    /// Charge the battery (update the zones temperature)
    /// It follows the same methodology as energy_demand function
    fn charge_battery(&self) -> anyhow::Result<(f64, Vec<f64>)> {
        let timestep = self.simulation_time.step_in_hours();
        let time_available = self.time_available(0., timestep);

        let pwr_in = self.electric_charge();
        let time_step_s = time_available * SECONDS_PER_HOUR as f64;
        let mut zone_temp_c_dist = self.zone_temp_c_dist_initial.read().clone();

        // Processing HB zones
        let (_, zone_temp_c_dist, _, energy_charged_during_battery_time_step) = self
            .process_heat_battery_zones(
                0.,
                &mut zone_temp_c_dist,
                0.,
                time_step_s,
                0.,
                Some(pwr_in),
                Some(HeatBatteryPcmOperationMode::OnlyCharging),
            )?;

        self.energy_charged
            .fetch_add(energy_charged_during_battery_time_step, Ordering::SeqCst);
        *self.zone_temp_c_dist_initial.write() = zone_temp_c_dist.clone();

        Ok((energy_charged_during_battery_time_step, zone_temp_c_dist))
    }

    fn battery_heat_loss(&self) -> anyhow::Result<(f64, Vec<f64>)> {
        // Battery losses
        let timestep = self.simulation_time.step_in_hours();
        let time_step_s = timestep * SECONDS_PER_HOUR as f64;

        let mut zone_temp_c_dist = self.zone_temp_c_dist_initial.read().clone();

        // Processing HB zones
        let (_, zone_temp_c_dist, energy_loss, _) = self.process_heat_battery_zones(
            22.,
            &mut zone_temp_c_dist,
            0.,
            time_step_s,
            0.,
            Some(-self.max_rated_losses),
            Some(HeatBatteryPcmOperationMode::Losses),
        )?;

        *self.zone_temp_c_dist_initial.write() = zone_temp_c_dist.clone();

        // Equivalent of using Python's math.fsum instead of sum() for better numerical accuracy with floating point arithmetic
        Ok((
            FSum::with_all(&energy_loss).value() / KILOJOULES_PER_KILOWATT_HOUR as f64,
            zone_temp_c_dist,
        ))
    }

    fn get_temp_hot_water(
        &self,
        inlet_temp: f64,
        volume: f64,
        setpoint_temp: f64,
    ) -> anyhow::Result<f64> {
        let total_time_s = volume / self.flow_rate_l_per_min * SECONDS_PER_MINUTE as f64;

        let time_step_s = self.hb_time_step;

        let pwr_in = self.electric_charge();

        // Initial Reynold number
        let mut water_kinematic_viscosity_m2_per_s =
            Self::calculate_water_kinematic_viscosity_m2_per_s(
                self.initial_inlet_temp,
                self.estimated_outlet_temp,
            );
        let mut reynold_number_at_1_l_per_min = Self::calculate_reynold_number_at_1_l_per_min(
            water_kinematic_viscosity_m2_per_s,
            self.velocity_in_hex_tube,
            self.capillary_diameter_m,
        );

        let flow_rate_kg_per_s =
            (self.flow_rate_l_per_min / SECONDS_PER_MINUTE as f64) * WATER.density();

        let mut zone_temp_c_dist = self.zone_temp_c_dist_initial.read().clone();
        let mut inlet_temp_c = inlet_temp;
        let mut outlet_temp_c = inlet_temp_c; // initialise, though expectation is this will be overridden in loop

        let n_time_steps = if total_time_s > time_step_s {
            (total_time_s / time_step_s) as usize
        } else {
            1
        };

        for _ in 0..n_time_steps {
            (outlet_temp_c, zone_temp_c_dist, _, _) = self.process_heat_battery_zones(
                inlet_temp_c,
                &mut zone_temp_c_dist,
                flow_rate_kg_per_s,
                time_step_s,
                reynold_number_at_1_l_per_min,
                pwr_in.into(),
                None,
            )?;

            // RN for next time step
            water_kinematic_viscosity_m2_per_s =
                Self::calculate_water_kinematic_viscosity_m2_per_s(inlet_temp_c, outlet_temp_c);
            reynold_number_at_1_l_per_min = Self::calculate_reynold_number_at_1_l_per_min(
                water_kinematic_viscosity_m2_per_s,
                self.velocity_in_hex_tube,
                self.capillary_diameter_m,
            );

            inlet_temp_c = outlet_temp_c;
        }

        Ok(min_of_2(outlet_temp_c, setpoint_temp))
    }

    /// Calculate the maximum energy output of the heat battery, accounting
    /// for time spent on higher-priority services.
    fn energy_output_max(&self, temp_output: f64, time_start: Option<f64>) -> anyhow::Result<f64> {
        // Return the energy the battery can provide assuming the HB temperature inlet
        // is constant during HEM time step equal to the required emitter temperature (temp_output)
        // Maximum energy for a given HB zones temperature distribution and inlet temperature.
        // The calculation methodology is the same as described in the demand_energy function.
        let time_start = time_start.unwrap_or(0.);
        let timestep = self.simulation_time.step_in_hours();
        let time_available = self.time_available(time_start, timestep);
        let total_time_s = time_available * SECONDS_PER_HOUR as f64;

        // time_step_s for HB calculation is a sensitive inputs for the process as, the longer it is, the
        // lower the accuracy due to maintaining Reynolds number working in intervals where the properties
        // of the fluid have changed sufficiently to degrade the accuracy of the calculation.
        // This is critical for the demand_energy function but less so for the energy_output_max as this
        // only provides an estimation of the heat capacity of the battery and can be slightly overestitmated
        // with a longer time step that reduces the calculation time, which is a critical factor for HEM.
        // However, from current testing, time_step_s longer than 100 might cause instabilities in the calculation
        // leading to failure to complete. Thus, we are capping the max timestep to 100 seconds.
        let time_step_s = (self.hb_time_step * 5.).min(100.);

        let pwr_in = self.electric_charge();

        // Initial Reynold number
        let mut water_kinematic_viscosity_m2_per_s =
            Self::calculate_water_kinematic_viscosity_m2_per_s(
                self.initial_inlet_temp,
                self.estimated_outlet_temp,
            );
        let mut reynold_number_at_1_l_per_min = Self::calculate_reynold_number_at_1_l_per_min(
            water_kinematic_viscosity_m2_per_s,
            self.velocity_in_hex_tube,
            self.capillary_diameter_m,
        );

        let flow_rate_kg_per_s =
            (self.flow_rate_l_per_min / SECONDS_PER_MINUTE as f64) * WATER.density();

        let mut zone_temp_c_dist = self.zone_temp_c_dist_initial.read().deref().clone();
        let mut energy_delivered_hb = 0.;
        let inlet_temp_c = temp_output;
        let n_time_steps = (total_time_s / time_step_s) as usize;

        for _ in 0..n_time_steps {
            // Processing HB zones
            let (outlet_temp_c, zone_temp_c_dist_new, energy_transf_delivered, _) = self
                .process_heat_battery_zones(
                    inlet_temp_c,
                    &mut zone_temp_c_dist,
                    flow_rate_kg_per_s,
                    time_step_s,
                    reynold_number_at_1_l_per_min,
                    pwr_in.into(),
                    None,
                )?;

            zone_temp_c_dist = zone_temp_c_dist_new;

            // RN for next time step
            water_kinematic_viscosity_m2_per_s =
                Self::calculate_water_kinematic_viscosity_m2_per_s(inlet_temp_c, outlet_temp_c);
            reynold_number_at_1_l_per_min = Self::calculate_reynold_number_at_1_l_per_min(
                water_kinematic_viscosity_m2_per_s,
                self.velocity_in_hex_tube,
                self.capillary_diameter_m,
            );

            // Equivalent of using Python's math.fsum instead of sum() for better numerical accuracy with floating point arithmetic
            let energy_delivered_ts = FSum::with_all(&energy_transf_delivered).value();

            if outlet_temp_c > temp_output {
                // In this new method, adjust total energy to make more real with the 6 ts we have configured
                energy_delivered_hb += energy_delivered_ts;
            } else {
                break;
            }
        }

        if energy_delivered_hb < 0. {
            energy_delivered_hb = 0.
        }

        Ok(energy_delivered_hb * self.n_units as f64)
    }

    fn first_call(&self) {
        self.flag_first_call.store(false, Ordering::SeqCst);
    }

    fn demand_energy(
        &self,
        service_name: &str,
        service_type: HeatingServiceType,
        energy_output_required: f64,
        temp_return_feed: f64,
        temp_output: Option<f64>,
        service_on: bool,
        time_start: Option<f64>,
        update_heat_source_state: Option<bool>,
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
        if self.flag_first_call.load(Ordering::SeqCst) {
            self.first_call();
        }

        let pwr_in = if self.simultaneous_charging_and_discharging {
            self.electric_charge()
        } else {
            0.
        };

        // Distributing energy demand through all units
        let energy_demand = energy_output_required / self.n_units as f64;

        // Initial Reynold number
        let mut water_kinematic_viscosity_m2_per_s =
            Self::calculate_water_kinematic_viscosity_m2_per_s(
                self.initial_inlet_temp,
                self.estimated_outlet_temp,
            );
        let mut reynold_number_at_1_l_per_min = Self::calculate_reynold_number_at_1_l_per_min(
            water_kinematic_viscosity_m2_per_s,
            self.velocity_in_hex_tube,
            self.capillary_diameter_m,
        );

        let flow_rate_kg_per_s =
            (self.flow_rate_l_per_min / SECONDS_PER_MINUTE as f64) * WATER.density();

        let mut energy_delivered_hb = 0.;
        let inlet_temp_c = temp_return_feed;
        let mut zone_temp_c_dist = self.zone_temp_c_dist_initial.read().clone();

        if energy_output_required < 0.
            || is_close!(energy_output_required, 0., rel_tol = 1e-09, abs_tol = 1e-10)
        {
            if update_heat_source_state {
                self.service_results.write().push(HeatBatteryResult {
                    service_name: service_name.into(),
                    service_type: service_type.into(),
                    service_on,
                    energy_output_required,
                    temp_output,
                    temp_inlet: temp_return_feed,
                    time_running: 0.,
                    energy_delivered_hb: 0.,
                    energy_delivered_backup: 0.,
                    energy_delivered_total: 0.,
                    energy_charged_during_service: 0.,
                    hb_zone_temperatures: zone_temp_c_dist,
                    current_hb_power: Default::default(),
                });
            }
            return Ok(0.);
        }

        let mut time_step_s = 1.;
        let mut time_running_current_service = 0.;

        let mut energy_charged = 0.;

        let mut outlet_temp_c = None;

        while time_step_s > 0. {
            // Processing HB zones
            let (
                outlet_temp_c_new,
                zone_temp_c_dist_new,
                energy_transf_delivered,
                energy_charged_during_battery_time_step,
            ) = self.process_heat_battery_zones(
                inlet_temp_c,
                &mut zone_temp_c_dist,
                flow_rate_kg_per_s,
                time_step_s,
                reynold_number_at_1_l_per_min,
                Some(pwr_in),
                None,
            )?;

            zone_temp_c_dist = zone_temp_c_dist_new;
            outlet_temp_c = Some(outlet_temp_c_new);

            if update_heat_source_state {
                self.energy_charged
                    .fetch_add(energy_charged_during_battery_time_step, Ordering::SeqCst);
            }
            energy_charged += energy_charged_during_battery_time_step;

            time_running_current_service += time_step_s;

            // RN for next time step
            water_kinematic_viscosity_m2_per_s = Self::calculate_water_kinematic_viscosity_m2_per_s(
                temp_return_feed,
                outlet_temp_c_new,
            );
            reynold_number_at_1_l_per_min = Self::calculate_reynold_number_at_1_l_per_min(
                water_kinematic_viscosity_m2_per_s,
                self.velocity_in_hex_tube,
                self.capillary_diameter_m,
            );

            let energy_transf_delivered_sum = FSum::with_all(&energy_transf_delivered).value();

            //  We assume the heat battery controls would stop the pump if the HB is absorbing energy from the emitters loop instead of contributing to it
            if energy_transf_delivered_sum < 0. {
                // Break prevents negative energy output by stopping before the current sub-timestep's result
                // is added to energy_delivered_HB. This occurs when the heat battery zones have cooled to the
                // point where they would absorb heat from the inlet flow rather than deliver it. By breaking
                // here, energy_delivered_HB retains only the positive contributions from previous sub-timesteps
                // where the battery was actively heating the flow, ensuring the function never returns negative
                // energy delivery values.
                break;
            }

            // Equivalent of using Python's math.fsum instead of sum() for better numerical accuracy with floating point arithmetic
            let energy_delivered_ts: f64 =
                energy_transf_delivered_sum / KILOJOULES_PER_KILOWATT_HOUR as f64;
            energy_delivered_hb += energy_delivered_ts; // demand_per_time_step_kwh
                                                        // balance = total_energy - energy_charged
            let max_instant_power = energy_delivered_ts / time_step_s;

            if max_instant_power > 0. {
                time_step_s = (energy_demand - energy_delivered_hb) / max_instant_power;
            }

            if time_step_s > self.hb_time_step {
                time_step_s = self.hb_time_step;
            }

            if is_close!(
                energy_demand,
                energy_delivered_hb,
                rel_tol = 1e-09,
                abs_tol = 1e-10
            ) || energy_delivered_hb > energy_demand
            {
                break;
            }

            if time_running_current_service + time_step_s > time_available * SECONDS_PER_HOUR as f64
            {
                time_step_s =
                    time_available * SECONDS_PER_HOUR as f64 - time_running_current_service;
            }
        }

        if update_heat_source_state {
            *self.zone_temp_c_dist_initial.write() = zone_temp_c_dist.clone();

            self.total_time_running_current_timestep.fetch_add(
                time_running_current_service / SECONDS_PER_HOUR as f64,
                Ordering::SeqCst,
            );

            // Track pump running time (only for regular DHW and space heating)
            // Direct DHW services don't use circulation pumps
            match service_type {
                HeatingServiceType::DomesticHotWaterRegular | HeatingServiceType::Space => {
                    self.pump_running_time_current_timestep.fetch_add(
                        time_running_current_service / SECONDS_PER_HOUR as f64,
                        Ordering::SeqCst,
                    );
                }
                HeatingServiceType::DomesticHotWaterDirect => (), // Direct DHW doesn't use circulation pump
                _ => bail!("Unexpected service type: {service_type}"),
            }

            let current_hb_power = if time_running_current_service > 0. {
                energy_delivered_hb * SECONDS_PER_HOUR as f64 / time_running_current_service
            } else {
                Default::default()
            };
            // TODO (from Python) Clarify whether Heat Batteries can have direct electric backup if depleted
            self.service_results.write().push(HeatBatteryResult {
                service_name: service_name.into(),
                service_type: service_type.into(),
                service_on,
                energy_output_required,
                temp_output: outlet_temp_c,
                temp_inlet: temp_return_feed,
                time_running: time_running_current_service,
                energy_delivered_hb: energy_delivered_hb * self.n_units as f64,
                energy_delivered_backup: 0.,
                energy_delivered_total: energy_delivered_hb * self.n_units as f64 + 0.,
                energy_charged_during_service: energy_charged * self.n_units as f64,
                hb_zone_temperatures: zone_temp_c_dist,
                current_hb_power: current_hb_power * self.n_units as f64,
            });
        }

        Ok(energy_delivered_hb * self.n_units as f64)
    }

    /// Calculation of heat battery auxiliary energy consumption
    fn calc_auxiliary_energy(
        &self,
        _timestep: f64,
        time_remaining_current_timestep: f64,
        timestep_idx: usize,
    ) -> anyhow::Result<f64> {
        // Energy used by circulation pump (for regular hot water and space heating services)
        let mut energy_aux = self
            .pump_running_time_current_timestep
            .load(Ordering::SeqCst)
            * self.power_circ_pump;

        // Energy used in standby mode
        energy_aux += self.power_standby * time_remaining_current_timestep;

        self.energy_supply_connection
            .demand_energy(energy_aux, timestep_idx)?;

        Ok(energy_aux)
    }

    /// Calculations to be done at the end of each timestep
    pub(crate) fn timestep_end(&self, timestep_idx: usize) -> anyhow::Result<()> {
        let timestep = self.simulation_time.step_in_hours();
        let time_remaining_current_timestep = timestep
            - self
                .total_time_running_current_timestep
                .load(Ordering::SeqCst);

        if self.flag_first_call.load(Ordering::SeqCst) {
            self.first_call();
        }
        self.flag_first_call.store(true, Ordering::SeqCst);

        // Calculating auxiliary energy to provide services during timestep
        let energy_aux =
            self.calc_auxiliary_energy(timestep, time_remaining_current_timestep, timestep_idx)?;

        let (battery_losses, zone_temp_c_after_losses) = self.battery_heat_loss()?;
        self.battery_losses.store(battery_losses, Ordering::SeqCst);

        // Charging battery for the remainder of the timestep
        let (end_of_ts_charge, zone_temp_c_after_charging) = if self
            .charge_control
            .is_on(&self.simulation_time.current_iteration())
        {
            self.charge_battery()?
        } else {
            (0., self.zone_temp_c_dist_initial.read().clone())
        };

        self.energy_supply_connection.demand_energy(
            self.energy_charged.load(Ordering::SeqCst) * self.n_units as f64,
            timestep_idx,
        )?;

        // If detailed results are to be output, save the results from the current timestep
        if let Some(detailed_results) = self.detailed_results.as_ref() {
            let service_results = self.service_results.read();
            let services_called: IndexMap<&String, &HeatBatteryResult> = service_results
                .iter()
                .map(|result| (&result.service_name, result))
                .collect();

            // Ensure all registered services have an entry in the results
            let mut ordered_service_results =
                Vec::with_capacity(self.energy_supply_connections.len());
            let initial_temps = self.zone_temp_c_dist_initial.read().clone();

            for service_name in self.energy_supply_connections.keys() {
                if let Some(result) = services_called.get(service_name) {
                    // Service was called, use its results
                    ordered_service_results.push((*result).clone());
                } else {
                    // Service was not called, create a placeholder entry
                    ordered_service_results.push(HeatBatteryResult {
                        service_name: service_name.clone(),
                        service_type: None, // Unknown since service wasn't called
                        service_on: false,
                        energy_output_required: 0.,
                        temp_output: None,
                        temp_inlet: 0.,
                        time_running: 0.,
                        energy_delivered_hb: 0.,
                        energy_delivered_backup: 0.,
                        energy_delivered_total: 0.,
                        energy_charged_during_service: 0.,
                        hb_zone_temperatures: initial_temps.clone(),
                        current_hb_power: 0.,
                    });
                }
            }

            // Add auxiliary results at the end
            let n_units = self.n_units as f64;
            let battery_losses = self.battery_losses.load(Ordering::SeqCst) * n_units;
            let total_charge = self.energy_charged.load(Ordering::SeqCst) * n_units;

            detailed_results.write().push(HeatBatteryTimestepResult {
                results: ordered_service_results,
                summary: HeatBatteryTimestepSummary {
                    energy_aux: energy_aux * n_units,
                    battery_losses,
                    temps_after_losses: zone_temp_c_after_losses,
                    total_charge,
                    end_of_timestep_charge: end_of_ts_charge * n_units,
                    hb_after_only_charge_zone_temp: zone_temp_c_after_charging,
                },
            });
        }

        self.total_time_running_current_timestep
            .store(Default::default(), Ordering::SeqCst);
        self.pump_running_time_current_timestep
            .store(Default::default(), Ordering::SeqCst);
        *self.service_results.write() = Default::default();
        self.energy_charged
            .store(Default::default(), Ordering::SeqCst);

        Ok(())
    }

    /// Output detailed results of heat battery calculation
    pub(crate) fn output_detailed_results(
        &self,
        hot_water_energy_output: &IndexMap<Arc<str>, Vec<ResultParamValue>>,
        hot_water_source_name_for_heat_battery_service: &IndexMap<Arc<str>, Arc<str>>,
    ) -> Result<(ResultsPerTimestep, ResultsAnnual), OutputDetailedResultsNotEnabledError> {
        let detailed_results = self
            .detailed_results
            .as_ref()
            .ok_or(OutputDetailedResultsNotEnabledError)?;

        let mut results_per_timestep: ResultsPerTimestep =
            [("auxiliary".into(), Default::default())].into();

        // Report auxiliary parameters (not specific to a service)
        for (parameter, param_unit, _) in AUX_PARAMETERS {
            if ["Temps_after_losses", "hb_after_only_charge_zone_temp"].contains(&parameter) {
                let mut labels: Option<Vec<Arc<str>>> = Default::default();
                for service_results in detailed_results.read().iter() {
                    let summary = &service_results.summary;
                    let param_values = match parameter {
                        "Temps_after_losses" => &summary.temps_after_losses,
                        "hb_after_only_charge_zone_temp" => &summary.hb_after_only_charge_zone_temp,
                        _ => unreachable!(),
                    };
                    // Determine the number of elements in the list for this parameter
                    if labels.is_none() {
                        labels = Some(
                            param_values
                                .iter()
                                .enumerate()
                                .map(|(i, _)| format!("{parameter}{i}").into())
                                .collect(),
                        );
                    }
                    for (label, result) in labels.as_ref().unwrap().iter().zip(param_values) {
                        results_per_timestep["auxiliary"]
                            .entry((label.clone(), param_unit.map(Into::into)))
                            .or_default()
                            .push(result.into());
                    }
                }
            } else {
                // Default behaviour for scalar parameters
                let mut param_results = vec![];
                for service_results in detailed_results.read().iter() {
                    let result = &service_results.summary.param(parameter);

                    param_results.push(result.clone());
                }
                results_per_timestep["auxiliary"].insert(
                    (parameter.into(), param_unit.map(Into::into)),
                    param_results,
                );
            }
        }

        // For each service, report required output parameters
        for (service_idx, service_name) in self.energy_supply_connections.keys().enumerate() {
            let service_name: Arc<str> = service_name.as_str().into();
            let mut current_results: ResultPerTimestep = Default::default();

            // Look up each required parameter
            for (parameter, param_unit, _) in OUTPUT_PARAMETERS {
                // Look up value of required parameter in each timestep
                for service_results in detailed_results.read().iter() {
                    let current_result = &service_results.results[service_idx];
                    if parameter == "hb_zone_temperatures" {
                        let labels: Vec<Arc<str>> = (0..current_result.hb_zone_temperatures.len())
                            .map(|i| format!("{parameter}{i}").into())
                            .collect_vec();
                        for (label, result) in labels
                            .into_iter()
                            .zip(current_result.hb_zone_temperatures.iter())
                        {
                            current_results
                                .entry((label.clone(), param_unit.map(|x| x.to_string().into())))
                                .or_default()
                                .push(result.into());
                        }
                    } else {
                        let result = current_result.param(parameter);
                        current_results
                            .entry((parameter.into(), param_unit.map(|x| x.to_string().into())))
                            .or_default()
                            .push(result);
                    }
                }
            }
            // For water heating service, record hot water energy delivered from tank
            current_results.insert(("energy_delivered_H4".into(), Some("kWh".into())), {
                if detailed_results.read().first().unwrap().results[service_idx].service_type
                    == Some(HeatingServiceType::DomesticHotWaterRegular)
                {
                    let energy_delivered_total_len = current_results
                        .get(&("energy_delivered_total".into(), Some("kWh".into())))
                        .map(|vec| vec.len())
                        .unwrap_or(0);
                    // For DHW, need to include storage and primary circuit losses.
                    // Can do this by replacing H5 numerator with total energy
                    // draw-off from hot water cylinder.
                    let hws_name = &hot_water_source_name_for_heat_battery_service[&service_name];

                    if !hot_water_energy_output.contains_key(hws_name) {
                        vec![ResultParamValue::Empty; energy_delivered_total_len]
                    } else {
                        hot_water_energy_output[hws_name].clone()
                    }
                } else {
                    current_results[&("energy_delivered_total".into(), Some("kWh".into()))].clone()
                }
            });

            results_per_timestep.insert(service_name.clone(), current_results);
        }

        let mut results_annual: ResultsAnnual = [
            (
                "Overall".into(),
                OUTPUT_PARAMETERS
                    .iter()
                    .filter_map(|(parameter, param_units, incl_in_manual)| {
                        incl_in_manual.then_some((
                            (
                                (*parameter).into(),
                                param_units.map(|x| x.to_string().into()),
                            ),
                            0.0f64.into(),
                        ))
                    })
                    .collect(),
            ),
            ("auxiliary".into(), Default::default()),
        ]
        .into();
        results_annual["Overall"].insert(
            ("energy_delivered_H4".into(), Some("kWh".into())),
            0.0f64.into(),
        );
        // Report auxiliary parameters (not specific to a service)
        for (parameter, param_unit, incl_in_annual) in AUX_PARAMETERS.iter() {
            if *incl_in_annual {
                results_annual["auxiliary"].insert(
                    ((*parameter).into(), param_unit.map(Into::into)),
                    ResultParamValue::from(
                        FSum::with_all(
                            results_per_timestep["auxiliary"][&(
                                (*parameter).into(),
                                param_unit.map(|x| x.to_string().into()),
                            )]
                                .iter()
                                .map(ResultParamValue::as_f64),
                        )
                        .value(),
                    ),
                );
            }
        }
        // For each service, report required output parameters
        for service_name in self.energy_supply_connections.keys() {
            let service_name: Arc<str> = service_name.as_str().into();
            results_annual.insert(service_name.clone(), Default::default());
            for (parameter, param_unit, incl_in_annual) in OUTPUT_PARAMETERS {
                if incl_in_annual {
                    let parameter_annual_total = ResultParamValue::from(
                        FSum::with_all(
                            results_per_timestep[&service_name]
                                [&(parameter.into(), param_unit.map(Into::into))]
                                .iter()
                                .map(ResultParamValue::as_f64),
                        )
                        .value(),
                    );
                    results_annual[&service_name].insert(
                        (parameter.into(), param_unit.map(Into::into)),
                        parameter_annual_total.clone(),
                    );
                    *results_annual["Overall"]
                        .entry((parameter.into(), param_unit.map(Into::into)))
                        .or_insert(ResultParamValue::Number(0.)) += parameter_annual_total;
                }
            }
            if results_per_timestep[&service_name]
                [&("energy_delivered_H4".into(), Some("kWh".into()))]
                .contains(&ResultParamValue::Empty)
            {
                results_annual.get_mut(&service_name).unwrap().insert(
                    ("energy_delivered_H4".into(), Some("kWh".into())),
                    ResultParamValue::Empty,
                );
                results_annual.get_mut("Overall").unwrap().insert(
                    ("energy_delivered_H4".into(), Some("kWh".into())),
                    ResultParamValue::Empty,
                );
            } else {
                results_annual.get_mut(&service_name).unwrap().insert(
                    ("energy_delivered_H4".into(), Some("kWh".into())),
                    ResultParamValue::from(
                        FSum::with_all(
                            results_per_timestep[&service_name]
                                [&("energy_delivered_H4".into(), Some("kWh".into()))]
                                .iter()
                                .map(ResultParamValue::as_f64),
                        )
                        .value(),
                    ),
                );

                if results_annual["Overall"][&("energy_delivered_H4".into(), Some("kWh".into()))]
                    != ResultParamValue::Empty
                {
                    let service_energy_delivered = results_annual[&service_name]
                        [&("energy_delivered_H4".into(), Some("kWh".into()))]
                        .clone();
                    results_annual["Overall"]
                        [&("energy_delivered_H4".into(), Some("kWh".into()))] +=
                        service_energy_delivered;
                }
            }
        }

        Ok((results_per_timestep, results_annual))
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

#[derive(Debug, Error)]
#[error("Tried to call output_detailed_results when option to collect detailed results was not selected")]
pub(crate) struct OutputDetailedResultsNotEnabledError;

type ResultPerTimestep = IndexMap<(Arc<str>, Option<Arc<str>>), Vec<ResultParamValue>>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::common::{MockWaterSupply, WaterSupply};
    use crate::core::controls::time_control::{ChargeControl, Control};
    use crate::core::controls::time_control::{MockControl, SetpointTimeControl};
    use crate::core::energy_supply::energy_supply::{
        EnergySupply, EnergySupplyBuilder, EnergySupplyConnection,
    };
    use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
    use crate::external_conditions::{DaylightSavingsConfig, ExternalConditions};
    use crate::input::{
        ControlLogicType, ExternalSensor, FuelType, HeatBattery as HeatBatteryInput,
        HeatSourceWetDetails,
    };
    use crate::simulation_time::{SimulationTime, SimulationTimeIteration, SimulationTimeIterator};
    use approx::assert_relative_eq;
    use indexmap::indexmap;
    use itertools::Itertools;
    use parking_lot::RwLock;
    use rstest::*;
    use serde_json::json;
    use smartstring::alias::String;
    use std::sync::atomic::Ordering;
    use std::sync::Arc;

    const SERVICE_NAME: &str = "TestService";

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 2., 1.)
    }

    #[fixture]
    fn simulation_time_iterator(simulation_time: SimulationTime) -> Arc<SimulationTimeIterator> {
        simulation_time.iter().into()
    }

    #[fixture]
    fn simulation_time_iteration(
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
                    {"start360": 0, "end360": 45},
                    {"start360": 45, "end360": 90},
                ]
            ))
            .unwrap(),
        )
    }

    #[fixture]
    fn battery_control_off(
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> Control {
        create_control_with_value(
            false,
            external_conditions,
            external_sensor,
            simulation_time_iteration,
        )
    }

    #[fixture]
    fn battery_control_on(
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> Control {
        create_control_with_value(
            true,
            external_conditions,
            external_sensor,
            simulation_time_iteration,
        )
    }

    fn create_control_with_value(
        boolean: bool,
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> Control {
        Control::Charge(
            ChargeControl::new(
                ControlLogicType::Manual,
                vec![boolean],
                &simulation_time_iteration,
                0,
                1.,
                vec![Some(0.2)],
                None,
                None,
                Some(external_conditions.into()),
                Some(external_sensor),
                None,
            )
            .unwrap(),
        )
    }

    fn create_heat_battery(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        control: Control,
        output_detailed_results: Option<bool>,
    ) -> Arc<RwLock<HeatBatteryPcm>> {
        let heat_battery_details: &HeatSourceWetDetails = &HeatSourceWetDetails::HeatBattery {
            battery: HeatBatteryInput::Pcm {
                energy_supply: "mains elec".into(),
                electricity_circ_pump: 0.06,
                electricity_standby: 0.0244,
                rated_charge_power: 20.0,
                max_rated_losses: 0.1,
                number_of_units: 1,
                control_charge: "hb_charge_control".into(),
                simultaneous_charging_and_discharging: false,
                heat_storage_k_j_per_k_above_phase_transition: 381.5,
                heat_storage_k_j_per_k_below_phase_transition: 305.2,
                heat_storage_k_j_per_k_during_phase_transition: 12317.,
                phase_transition_temperature_upper: 59.,
                phase_transition_temperature_lower: 57.,
                max_temperature: 80.,
                temp_init: 80.,
                velocity_in_hex_tube_at_1_l_per_min_m_per_s: 0.035,
                inlet_diameter_mm: 0.0065,
                a: 19.744,
                b: -105.5,
                flow_rate_l_per_min: 10.,
            },
        };

        let energy_supply: Arc<RwLock<EnergySupply>> = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time_iterator.total_steps())
                .build(),
        ));

        let energy_supply_connection: EnergySupplyConnection =
            EnergySupply::connection(energy_supply.clone(), "WaterHeating").unwrap();

        let heat_battery = Arc::new(RwLock::new(HeatBatteryPcm::new(
            heat_battery_details,
            control.into(),
            energy_supply,
            energy_supply_connection,
            simulation_time_iterator,
            Some(8),
            Some(20.),
            None,
            None,
            output_detailed_results,
        )));

        HeatBatteryPcm::create_service_connection(heat_battery.clone(), SERVICE_NAME).unwrap();

        heat_battery
    }

    fn create_setpoint_time_control(schedule: Vec<Option<f64>>) -> Control {
        Control::SetpointTime(SetpointTimeControl::new(
            schedule,
            0,
            1.,
            Default::default(),
            Default::default(),
            1.,
        ))
    }

    fn get_service_names_from_results(heat_battery: Arc<RwLock<HeatBatteryPcm>>) -> Vec<String> {
        heat_battery
            .read()
            .service_results
            .read()
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

        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);

        let heat_battery_service = HeatBatteryPcmServiceSpace::new(
            heat_battery.clone(),
            SERVICE_NAME.into(),
            service_control_on.into(),
        );

        assert!(heat_battery_service.is_on(simulation_time_iteration));

        let service_control_off: Control = create_setpoint_time_control(vec![None]);

        let heat_battery_service: HeatBatteryPcmServiceSpace = HeatBatteryPcmServiceSpace::new(
            heat_battery,
            SERVICE_NAME.into(),
            service_control_off.into(),
        );

        assert!(!heat_battery_service.is_on(simulation_time_iteration));
    }

    #[fixture]
    fn heat_battery_service_water_direct(
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) -> HeatBatteryPcmServiceWaterDirect {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);
        let cold_feed =
            WaterSupply::ColdWaterSource(Arc::new(ColdWaterSource::new(vec![1.0, 1.2], 0, 1.)));
        let service_name = "WaterHeating".into();

        HeatBatteryPcmServiceWaterDirect::new(heat_battery, service_name, 60., cold_feed.clone())
    }

    #[rstest]
    fn test_get_cold_water_source_for_water_direct(
        heat_battery_service_water_direct: HeatBatteryPcmServiceWaterDirect,
    ) {
        let expected =
            WaterSupply::ColdWaterSource(Arc::new(ColdWaterSource::new(vec![1.0, 1.2], 0, 1.)));
        let actual = heat_battery_service_water_direct.get_cold_water_source();

        match (actual, expected) {
            (WaterSupply::ColdWaterSource(actual), WaterSupply::ColdWaterSource(expected)) => {
                assert_eq!(actual, &expected);
            }
            _ => panic!("Expected ColdWaterSource variant"),
        }
    }

    #[rstest]
    #[ignore = "test to be updated during migration to 1.0.0a6"]
    fn test_get_temp_hot_water_for_water_direct(
        mut heat_battery_service_water_direct: HeatBatteryPcmServiceWaterDirect,
        simulation_time_iteration: SimulationTimeIteration,
    ) {
        heat_battery_service_water_direct.cold_feed = WaterSupply::Mock(MockWaterSupply::new(25.));

        let expected = vec![(60., 20.)];
        let actual = heat_battery_service_water_direct
            .get_temp_hot_water(20., None, simulation_time_iteration)
            .unwrap();

        assert_eq!(actual, expected)
    }

    // skipping following python tests due to mocking:
    // test_demand_hot_water, test_demand_hot_water_fallback_path

    fn create_service_water_regular_with_controls(
        battery_control: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) -> HeatBatteryPcmServiceWaterRegular {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control, None);
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
        let cold_water_source =
            WaterSupply::ColdWaterSource(Arc::new(ColdWaterSource::new(vec![1.0, 1.2], 0, 1.)));

        HeatBatteryPcmServiceWaterRegular::new(
            heat_battery,
            SERVICE_NAME.into(),
            cold_water_source,
            Arc::new(control_min),
            Arc::new(control_max),
        )
    }

    // test_service_is_on_without_control
    #[rstest]
    fn test_service_with_no_service_control_is_always_on_for_water_regular(
        simulation_time_iteration: SimulationTimeIteration,
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let heat_battery_service = create_service_water_regular_with_controls(
            battery_control_off,
            simulation_time_iterator,
        );

        assert!(heat_battery_service.is_on(simulation_time_iteration));
    }

    #[rstest]
    fn test_setpnt_for_water_regular(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        simulation_time: SimulationTime,
        battery_control_off: Control,
    ) {
        let service = create_service_water_regular_with_controls(
            battery_control_off,
            simulation_time_iterator,
        );

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let (control_min, control_max) = service.setpnt(t_it);

            assert_eq!(
                control_min,
                [
                    Some(52.),
                    None,
                    None,
                    None,
                    Some(52.),
                    Some(52.),
                    Some(52.),
                    Some(52.)
                ][t_idx]
            );
            assert_eq!(control_max, Some(55.));
        }
    }

    // In Python this is test_demand_energy_service_off
    #[rstest]
    fn test_demand_energy_returns_zero_when_service_control_is_off_for_water_regular(
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
        let energy_demand = 10.;
        let temp_flow = 55.;
        let temp_return = 40.;

        let service_control_off = Arc::new(create_setpoint_time_control(vec![None]));

        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on, None);
        let cold_water_source = ColdWaterSource::new(vec![1.0, 1.2], 0, 1.);
        let heat_battery_service: HeatBatteryPcmServiceWaterRegular =
            HeatBatteryPcmServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                WaterSupply::ColdWaterSource(Arc::new(cold_water_source)),
                service_control_off.clone(),
                service_control_off,
            );

        let result = heat_battery_service
            .demand_energy(
                energy_demand,
                Some(temp_flow),
                temp_return,
                None,
                simulation_time_iteration,
            )
            .unwrap();

        assert_eq!(result, 0.);
    }

    // In Python this is test_energy_output_max_service_on
    #[rstest]
    #[ignore = "test to be updated during migration to 1.0.0a6"]
    fn test_energy_output_max_when_service_control_on_for_water_regular(
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
        let heat_battery_service = create_service_water_regular_with_controls(
            battery_control_on,
            simulation_time_iterator,
        );

        let temp_flow = 50.0;
        let temp_return = 40.0;
        let result = heat_battery_service
            .energy_output_max(temp_flow, temp_return, simulation_time_iteration)
            .unwrap();

        assert_relative_eq!(result, 72279.10023958197);
    }

    #[rstest]
    #[ignore = "test to be updated during migration to 1.0.0a6"]
    fn test_energy_output_max_service_off_for_water_regular(
        // In Python this is test_energy_output_max_service_off
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_off: Control,
    ) {
        let heat_battery_service = create_service_water_regular_with_controls(
            battery_control_off,
            simulation_time_iterator,
        );

        let temp_flow = 50.0;
        let temp_return = 40.0;
        let result = heat_battery_service
            .energy_output_max(temp_flow, temp_return, simulation_time_iteration)
            .unwrap();

        assert_relative_eq!(result, 28882.5139822234, epsilon = 1e-7);
    }

    #[rstest]
    fn test_temp_setpnt_for_space(
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_off: Control,
    ) {
        let first_scheduled_temp = Some(21.);
        let ctrl: Control = create_setpoint_time_control(vec![first_scheduled_temp]);
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);
        let heat_battery_space =
            HeatBatteryPcmServiceSpace::new(heat_battery, SERVICE_NAME.into(), ctrl.into());

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
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);
        let heat_battery_space =
            HeatBatteryPcmServiceSpace::new(heat_battery, SERVICE_NAME.into(), ctrl.into());

        assert_eq!(
            heat_battery_space.in_required_period(simulation_time_iteration),
            Some(true)
        );
    }

    #[rstest]
    fn test_demand_energy_service_off_for_space(
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_off: Control,
    ) {
        let energy_demand = 10.;
        let temp_return = 40.;
        let temp_flow = 1.;
        let time_start = 0.2;
        let ctrl: Control = create_setpoint_time_control(vec![None]);
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);
        let heat_battery_space =
            HeatBatteryPcmServiceSpace::new(heat_battery, SERVICE_NAME.into(), ctrl.into());
        let result = heat_battery_space
            .demand_energy(
                energy_demand,
                temp_flow,
                temp_return,
                Some(time_start),
                None,
                simulation_time_iteration,
            )
            .unwrap();
        assert_eq!(result, 0.);
    }

    // skipping python's test_energy_output_max_service_on due to mocking

    // in Python this test is called test_energy_output_max_service_off
    #[rstest]
    fn test_energy_output_max_service_off_for_space(
        battery_control_on: Control,
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let temp_output = 70.;
        let temp_return = 40.;
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on, None);
        let service_control_off: Control = create_setpoint_time_control(vec![None]);

        let heat_battery_service: HeatBatteryPcmServiceSpace = HeatBatteryPcmServiceSpace::new(
            heat_battery,
            SERVICE_NAME.into(),
            service_control_off.into(),
        );

        let result = heat_battery_service
            .energy_output_max(temp_output, temp_return, None, simulation_time_iteration)
            .unwrap();

        assert_relative_eq!(result, 0.);
    }

    #[rstest]
    fn test_create_service_connection(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on, None);
        let create_connection_result =
            HeatBatteryPcm::create_service_connection(heat_battery.clone(), "new service");
        assert!(create_connection_result.is_ok());
        assert!(heat_battery
            .read()
            .energy_supply_connections
            .contains_key("new service"));
        let create_connection_result =
            HeatBatteryPcm::create_service_connection(heat_battery, "new service");
        assert!(create_connection_result.is_err()) // second attempt to create a service connection with same name should error
    }

    #[rstest]
    fn test_create_service_hot_water_direct(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on, None);
        let cold_feed =
            WaterSupply::ColdWaterSource(Arc::new(ColdWaterSource::new(vec![1.0, 1.2], 0, 1.)));
        let service = HeatBatteryPcm::create_service_hot_water_direct(
            heat_battery.clone(),
            "new_service",
            60.,
            cold_feed.clone(),
        )
        .unwrap();

        let actual = service.get_cold_water_source();

        match (actual, cold_feed) {
            (WaterSupply::ColdWaterSource(actual), WaterSupply::ColdWaterSource(cold_feed)) => {
                assert_eq!(actual, &cold_feed);
            }
            _ => panic!("Expected ColdWaterSource variant"),
        }

        assert!(heat_battery
            .read()
            .energy_supply_connections
            .contains_key("new_service"));
    }

    #[rstest]
    fn test_create_service_space_heating(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        simulation_time_iteration: SimulationTimeIteration,
        battery_control_off: Control,
    ) {
        let control = Arc::new(create_setpoint_time_control(vec![Some(21.0)]));
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);
        let service = HeatBatteryPcm::create_service_space_heating(
            heat_battery.clone(),
            "new_service",
            control,
        )
        .unwrap();

        assert!(service.is_on(simulation_time_iteration));
        assert!(heat_battery
            .read()
            .energy_supply_connections
            .contains_key("new_service"));
    }

    #[rstest]
    fn test_electric_charge(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_off: Control,
        battery_control_on: Control,
    ) {
        // electric charge should be 0 when battery control is off
        let heat_battery =
            create_heat_battery(simulation_time_iterator.clone(), battery_control_off, None);
        assert_relative_eq!(heat_battery.read().electric_charge(), 0.0);

        // electric charge should be calculated when battery control is on
        let heat_battery =
            create_heat_battery(simulation_time_iterator.clone(), battery_control_on, None);
        assert_relative_eq!(heat_battery.read().electric_charge(), 20.0);
    }

    #[rstest]
    fn test_first_call(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
        simulation_time: SimulationTime,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on, None);
        for (t_idx, _) in simulation_time.iter().enumerate() {
            heat_battery.read().first_call();

            assert!(!heat_battery.read().flag_first_call.load(Ordering::SeqCst));

            heat_battery.read().timestep_end(t_idx).unwrap();
        }
    }

    #[rstest]
    #[ignore = "test to be updated during migration to 1.0.0a6"]
    fn test_demand_energy(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        simulation_time: SimulationTime,
        battery_control_on: Control,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on, None);

        let expected_zone_temp_c_dist = [
            vec![
                79.71165314809511,
                79.85379912318692,
                79.92587158056449,
                79.96241457173316,
                79.98094301175232,
                79.99033751063061,
                79.99510081553287,
                79.99751596017077,
            ], // First timestep
            vec![
                78.48854379731785,
                78.76743300209962,
                78.90934369283018,
                78.9815529739174,
                79.01829519996613,
                79.03699050188325,
                79.04650298972031,
                79.05134304583224,
            ], // Second timestep
        ];

        let service_name = "new_service";
        HeatBatteryPcm::create_service_connection(heat_battery.clone(), service_name).unwrap();

        for (t_idx, _) in simulation_time.iter().enumerate() {
            let demand_energy_actual = heat_battery
                .clone()
                .read()
                .demand_energy(
                    service_name,
                    HeatingServiceType::DomesticHotWaterRegular,
                    5.,
                    40.,
                    Some(52.5),
                    true,
                    Some(1.), // the Python here erroneously uses too many arguments to demand_energy so this is to fake the equivalent in the Rust, for example the Python True is understood as the number 1
                    None,
                )
                .unwrap();

            assert_relative_eq!(
                demand_energy_actual,
                [0.007714304589733515, 0.007530418147738887][t_idx]
            );

            let service_names_in_results = get_service_names_from_results(heat_battery.clone());

            assert!(service_names_in_results.contains(&service_name.into()));

            assert_eq!(heat_battery.read().charge_level, [0.0, 0.0][t_idx]);

            assert_relative_eq!(
                heat_battery
                    .read()
                    .total_time_running_current_timestep
                    .load(Ordering::SeqCst),
                [0.0002777777777777778, 0.0002777777777777778][t_idx]
            );

            assert_eq!(
                heat_battery.read().zone_temp_c_dist_initial.read().clone(),
                expected_zone_temp_c_dist[t_idx]
            );

            heat_battery.read().timestep_end(t_idx).unwrap();
        }
    }

    fn create_heat_battery_pcm(
        external_sensor: ExternalSensor,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        external_conditions: ExternalConditions,
    ) -> Arc<RwLock<HeatBatteryPcm>> {
        let control = Control::Charge(
            ChargeControl::new(
                ControlLogicType::Manual,
                vec![false],
                &simulation_time_iterator.current_iteration(),
                0,
                1.,
                vec![Some(0.2), Some(0.3)],
                None,
                None,
                Some(external_conditions.into()),
                Some(external_sensor),
                None,
            )
            .unwrap(),
        );
        let heat_battery = create_heat_battery(simulation_time_iterator, control, None);
        HeatBatteryPcm::create_service_connection(heat_battery.clone(), "new_service").unwrap();

        heat_battery
    }

    // skipping python's test_demand_energy_simultaneous_charging_and_discharging due to mocking

    #[rstest]
    #[ignore = "test to be updated during migration to 1.0.0a6"]
    fn test_demand_energy_simultaneous_no_temp_output(
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);
        assert_relative_eq!(
            heat_battery
                .read()
                .demand_energy(
                    SERVICE_NAME,
                    HeatingServiceType::DomesticHotWaterRegular,
                    0.08,
                    40.,
                    None,
                    true,
                    None,
                    None
                )
                .unwrap(),
            0.08021138263537801
        );

        assert_relative_eq!(
            heat_battery
                .read()
                .demand_energy(
                    SERVICE_NAME,
                    HeatingServiceType::DomesticHotWaterRegular,
                    0.06,
                    40.,
                    None,
                    true,
                    None,
                    None
                )
                .unwrap(),
            0.06018673551977593
        );

        // Battery losses
        assert_eq!(heat_battery.read().get_battery_losses(), 0.);
    }

    #[rstest]
    #[ignore = "test to be updated during migration to 1.0.0a6"]
    fn test_demand_energy_other(
        external_sensor: ExternalSensor,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        external_conditions: ExternalConditions,
    ) {
        let heat_battery = create_heat_battery_pcm(
            external_sensor.clone(),
            simulation_time_iterator.clone(),
            external_conditions.clone(),
        );
        assert_eq!(
            heat_battery
                .read()
                .demand_energy(
                    "new_service",
                    HeatingServiceType::DomesticHotWaterRegular,
                    0.08,
                    40.,
                    Some(40.),
                    true,
                    None,
                    None
                )
                .unwrap(),
            0.08021138263537801
        );

        let heat_battery = create_heat_battery_pcm(
            external_sensor.clone(),
            simulation_time_iterator.clone(),
            external_conditions.clone(),
        );
        heat_battery.write().hb_time_step = 119.;

        assert_eq!(
            heat_battery
                .read()
                .demand_energy(
                    "new_service",
                    HeatingServiceType::DomesticHotWaterRegular,
                    0.08,
                    40.,
                    Some(40.),
                    true,
                    None,
                    None
                )
                .unwrap(),
            0.08021138263537801
        );

        let heat_battery = create_heat_battery_pcm(
            external_sensor.clone(),
            simulation_time_iterator.clone(),
            external_conditions.clone(),
        );
        heat_battery.write().hb_time_step = 20.;

        assert_eq!(
            heat_battery
                .read()
                .demand_energy(
                    "new_service",
                    HeatingServiceType::DomesticHotWaterRegular,
                    0.08,
                    40.,
                    Some(80.),
                    true,
                    None,
                    None
                )
                .unwrap(),
            0.08021138263537801
        );

        let heat_battery = create_heat_battery_pcm(
            external_sensor,
            simulation_time_iterator,
            external_conditions,
        );

        assert_eq!(
            heat_battery
                .read()
                .demand_energy(
                    "new_service",
                    HeatingServiceType::DomesticHotWaterRegular,
                    0.08,
                    40.,
                    Some(79.),
                    true,
                    None,
                    None
                )
                .unwrap(),
            0.08021138263537801
        );
    }

    #[rstest]
    fn test_dhw_service_demand_hot_water(
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        simulation_time_iteration: SimulationTimeIteration,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);
        let cold_feed =
            WaterSupply::ColdWaterSource(Arc::new(ColdWaterSource::new(vec![1.0, 1.2], 0, 1.)));
        let service = HeatBatteryPcm::create_service_hot_water_direct(
            heat_battery,
            "dhw_complex",
            65., // High setpoint
            cold_feed.clone(),
        )
        .unwrap();

        let actual = service.get_cold_water_source();

        match (actual, cold_feed) {
            (WaterSupply::ColdWaterSource(actual), WaterSupply::ColdWaterSource(cold_feed)) => {
                assert_eq!(actual, &cold_feed);
            }
            _ => panic!("Expected ColdWaterSource variant"),
        }
        // Python tests with usage events here using mocking that is not easy to replicate

        // Test with no usage events
        let energy_no_usage = service
            .demand_hot_water(None, simulation_time_iteration)
            .unwrap();

        assert_eq!(energy_no_usage, 0.);
    }

    #[rstest]
    /// Check heat battery auxilary energy consumption
    fn test_calc_auxiliary_energy(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
        let heat_battery =
            create_heat_battery(simulation_time_iterator.clone(), battery_control_on, None);

        heat_battery
            .read()
            .calc_auxiliary_energy(1.0, 0.5, simulation_time_iterator.current_index())
            .unwrap();

        let results_by_end_user = heat_battery
            .read()
            .energy_supply
            .read()
            .results_by_end_user();

        let end_user_name: Arc<str> = heat_battery
            .read()
            .energy_supply_connection
            .end_user_name
            .to_string()
            .into();

        let results_by_end_user = results_by_end_user.get(&end_user_name).unwrap();

        assert_eq!(*results_by_end_user, vec![0.0122, 0.]);
    }

    #[rstest]
    #[ignore = "test to be updated during migration to 1.0.0a6"]
    fn test_timestep_end(
        external_sensor: ExternalSensor,
        external_conditions: ExternalConditions,
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        // not using the fixture here
        // because we need to set different charge_levels
        let battery_control_on: Control = Control::Charge(
            ChargeControl::new(
                ControlLogicType::Manual,
                vec![true, true, true],
                &simulation_time_iteration,
                0,
                1.,
                [1.0, 1.5].into_iter().map(Into::into).collect(),
                None,
                None,
                Some(external_conditions.into()),
                Some(external_sensor),
                None,
            )
            .unwrap(),
        );

        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on, None);
        let service_name = "new_timestep_end_service";
        HeatBatteryPcm::create_service_connection(heat_battery.clone(), service_name).unwrap();

        let t_idx = 0;
        heat_battery
            .read()
            .demand_energy(
                service_name,
                HeatingServiceType::DomesticHotWaterRegular,
                5.0,
                40.,
                Some(55.),
                true,
                None,
                None,
            )
            .unwrap();

        assert_relative_eq!(
            heat_battery
                .read()
                .total_time_running_current_timestep
                .load(Ordering::SeqCst),
            0.25690463025906096
        );

        let service_names_in_results = get_service_names_from_results(heat_battery.clone());

        assert!(service_names_in_results.contains(&service_name.into()));

        heat_battery.read().timestep_end(t_idx).unwrap();

        // Assertions to check if the internal state was updated correctly
        assert!(heat_battery.read().flag_first_call.load(Ordering::SeqCst)); // Python has double negative here

        assert_relative_eq!(
            heat_battery
                .read()
                .total_time_running_current_timestep
                .load(Ordering::SeqCst),
            0.0
        );
        assert_eq!(heat_battery.read().service_results.read().len(), 0);
    }

    #[rstest]
    #[ignore = "test to be updated during migration to 1.0.0a6"]
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
                &simulation_time_iterator.current_iteration(),
                0,
                1.,
                [1.5, 1.6].into_iter().map(Into::into).collect(), // these values change the result
                None,
                None,
                Some(external_conditions.clone().into()),
                Some(external_sensor.clone()),
                None,
            )
            .unwrap(),
        );

        let heat_battery =
            create_heat_battery(simulation_time_iterator.clone(), battery_control_on, None);

        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                heat_battery.read().energy_output_max(0., None).unwrap(),
                [108864.87597021714, 124118.95144251334][t_idx]
            );

            heat_battery.read().timestep_end(t_idx).unwrap();
        }

        let battery_control_on: Control = Control::Charge(
            ChargeControl::new(
                ControlLogicType::Manual,
                vec![true, true, true],
                &simulation_time_iterator.current_iteration(),
                0,
                1.,
                [1.5, 1.6].into_iter().map(Into::into).collect(), // these values change the result
                None,
                None,
                Some(external_conditions.into()),
                Some(external_sensor),
                None,
            )
            .unwrap(),
        );
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on, None);

        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                heat_battery.read().energy_output_max(90., None).unwrap(),
                [0., 72281.56558957469][t_idx]
            );

            heat_battery.read().timestep_end(t_idx).unwrap();
        }
    }

    #[rstest]
    fn test_get_zone_properties_losses(
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        // Test that get_zone_properties returns the correct energy_transf with losses model
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);
        let (energy_transf, _, _, _) = heat_battery.read().get_zone_properties(
            0,
            &HeatBatteryPcmOperationMode::Losses,
            &[42., 57., 58., 58., 59., 59., 60., 61.],
            40.,
            40.,
            5.,
            414.,
            0.16,
            20.,
        );

        assert_eq!(energy_transf, 0.625);
    }

    #[rstest]
    fn test_get_zone_properties_no_energy_transf(
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        // Test that get_zone_properties returns energy_transf as 0 with losses model and higher zone_temp_c_start than inlet_temp_c
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);
        let (energy_transf, _, _, _) = heat_battery.read().get_zone_properties(
            0,
            &HeatBatteryPcmOperationMode::Losses,
            &[42., 57., 58., 58., 59., 59., 60., 61.],
            45.,
            45.,
            5.,
            414.,
            0.16,
            20.,
        );

        assert_eq!(energy_transf, 0.);
    }

    // skipping python's test_get_zone_properties_invalid_mode as mode can't be invalid in rust

    #[rstest]
    fn test_calculate_zone_energy_required(
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);

        let required = heat_battery.read().calculate_zone_energy_required(50., 80.);

        assert_relative_eq!(required, -4347.7375);

        let required = heat_battery.read().calculate_zone_energy_required(58., 80.);

        assert_relative_eq!(required, -2541.0625);

        let required = heat_battery.read().calculate_zone_energy_required(60., 80.);

        assert_relative_eq!(required, -953.75);

        let required = heat_battery
            .read()
            .calculate_zone_energy_required(58., 58.5);

        assert_relative_eq!(required, -769.8125);

        let required = heat_battery
            .read()
            .calculate_zone_energy_required(55., 58.5);

        assert_relative_eq!(required, -2385.7375);

        let required = heat_battery
            .read()
            .calculate_zone_energy_required(60., 58.5);

        assert_relative_eq!(required, 71.53125);

        let required = heat_battery.read().calculate_zone_energy_required(50., 55.);

        assert_relative_eq!(required, -190.75);

        let required = heat_battery.read().calculate_zone_energy_required(58., 55.);

        assert_relative_eq!(required, 4618.875);

        let required = heat_battery.read().calculate_zone_energy_required(60., 55.);

        assert_relative_eq!(required, 238.4375);
    }

    #[rstest]
    fn test_process_zone_simultaneous_charging(
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);

        let (q_max_kj, energy_charged, energy_transf) = heat_battery
            .read()
            .process_zone_simultaneous_charging(58., 120., -2000., 1900., 0.);

        assert_relative_eq!(q_max_kj, 0.);
        assert_relative_eq!(energy_charged, 0.5555555555555556);
        assert_relative_eq!(energy_transf, -100.);

        let (q_max_kj, energy_charged, energy_transf) = heat_battery
            .read()
            .process_zone_simultaneous_charging(58., 120., -2000., -1900., 0.);

        assert_relative_eq!(q_max_kj, 0.);
        assert_relative_eq!(energy_charged, 0.5555555555555556);
        assert_relative_eq!(energy_transf, -3900.);

        let (q_max_kj, energy_charged, energy_transf) = heat_battery
            .read()
            .process_zone_simultaneous_charging(58., 120., -3000., -1900., 0.);

        assert_relative_eq!(q_max_kj, -451.4375);
        assert_relative_eq!(energy_charged, 0.7079340277777778);
        assert_relative_eq!(energy_transf, -4448.5625);
    }

    // skipping python's test_process_zone_simultaneous_charging_warning1 and
    // test_process_zone_simultaneous_charging_warning2 as we haven't incorporated these warnings

    #[rstest]
    #[case(55., 3000., -23.63695937090432)]
    #[case(58., 3000., 18.720183486238533)]
    #[case(60., 3000., 57.08244702443777)]
    #[case(55., 4000., -49.84927916120577)]
    #[case(58., 4000., -7.49213630406291)]
    #[case(60., 4000., 34.11500655307995)]
    #[case(60., 10., 59.79030144167759)]
    #[case(50., -4000., 72.70799475753604)]
    #[case(50., -3000., 58.77507509945603)]
    #[case(50., -100., 52.62123197903014)]
    #[case(58., -2000., 68.65399737876803)]
    #[case(58., -1000., 58.64950880896322)]
    #[case(60., -1000., 80.96985583224115)]
    fn test_calculate_new_zone_temperature(
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        #[case] zone_temp_c_start: f64,
        #[case] energy_transf: f64,
        #[case] expected: f64,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);
        let result = heat_battery
            .read()
            .calculate_new_zone_temperature(zone_temp_c_start, energy_transf);

        assert_relative_eq!(result, expected);
    }

    #[rstest]
    #[ignore = "test to be updated during migration to 1.0.0a6"]
    fn test_charge_battery_hydraulic(
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);

        assert_relative_eq!(
            heat_battery.write()._charge_battery_hydraulic(70.).unwrap(),
            0.
        );
        assert_relative_eq!(
            heat_battery.write()._charge_battery_hydraulic(80.).unwrap(),
            -138.85748246864733
        );
        assert_relative_eq!(
            heat_battery.write()._charge_battery_hydraulic(90.).unwrap(),
            -3814.99999900312
        );
    }

    #[rstest]
    #[ignore = "test to be updated during migration to 1.0.0a6"]
    fn test_get_temp_hot_water(
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);

        assert_relative_eq!(
            heat_battery
                .read()
                .get_temp_hot_water(50., 20., 80.)
                .unwrap(),
            79.70798180572169
        );
        assert_relative_eq!(
            heat_battery
                .read()
                .get_temp_hot_water(50., 10., 80.)
                .unwrap(),
            79.8652529090689
        );
        assert_relative_eq!(
            heat_battery
                .read()
                .get_temp_hot_water(40., 10., 80.)
                .unwrap(),
            79.81947841211459
        );
        assert_relative_eq!(
            heat_battery
                .read()
                .get_temp_hot_water(60., 1., 65.)
                .unwrap(),
            65.
        );
    }

    // skipping python's test_energy_output_max_negative as unable to replicate patch object

    #[fixture]
    fn heat_battery_no_service_connection(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) -> Arc<RwLock<HeatBatteryPcm>> {
        let heat_battery_details: &HeatSourceWetDetails = &HeatSourceWetDetails::HeatBattery {
            battery: HeatBatteryInput::Pcm {
                energy_supply: "mains elec".into(),
                electricity_circ_pump: 0.06,
                electricity_standby: 0.0244,
                rated_charge_power: 20.0,
                max_rated_losses: 0.1,
                number_of_units: 1,
                control_charge: "hb_charge_control".into(),
                simultaneous_charging_and_discharging: false,
                heat_storage_k_j_per_k_above_phase_transition: 381.5,
                heat_storage_k_j_per_k_below_phase_transition: 305.2,
                heat_storage_k_j_per_k_during_phase_transition: 12317.,
                phase_transition_temperature_upper: 59.,
                phase_transition_temperature_lower: 57.,
                max_temperature: 80.,
                temp_init: 80.,
                velocity_in_hex_tube_at_1_l_per_min_m_per_s: 0.035,
                inlet_diameter_mm: 0.0065,
                a: 19.744,
                b: -105.5,
                flow_rate_l_per_min: 10.,
            },
        };

        let energy_supply: Arc<RwLock<EnergySupply>> = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time_iterator.total_steps())
                .build(),
        ));

        let energy_supply_connection: EnergySupplyConnection =
            EnergySupply::connection(energy_supply.clone(), "WaterHeating").unwrap();

        Arc::new(RwLock::new(HeatBatteryPcm::new(
            heat_battery_details,
            battery_control_on.into(),
            energy_supply,
            energy_supply_connection,
            simulation_time_iterator,
            Some(8),
            Some(20.),
            None,
            None,
            Some(true),
        )))
    }

    #[rstest]
    #[ignore = "test to be updated during migration to 1.0.0a6"]
    fn test_output_detailed_results_water_regular(
        simulation_time: SimulationTime,
        heat_battery_no_service_connection: Arc<RwLock<HeatBatteryPcm>>,
    ) {
        let heat_battery = heat_battery_no_service_connection;
        let cold_feed =
            WaterSupply::ColdWaterSource(Arc::new(ColdWaterSource::new(vec![1.0, 1.2], 0, 1.)));
        let service_name = "new_service";

        HeatBatteryPcm::create_service_hot_water_regular(
            heat_battery.clone(),
            service_name,
            cold_feed.clone(),
            Arc::new(Control::Mock(MockControl::default())),
            Arc::new(Control::Mock(MockControl::default())),
        )
        .unwrap();

        let mut expected_results_per_timestep: ResultsPerTimestep = indexmap! {
            "auxiliary".into() => indexmap! {
                ("energy_aux".into(), Some("kWh".into())) => vec![0.06.into(), 0.06.into()],
                ("battery_losses".into(), Some("kWh".into())) => vec![0.1.into(), 0.1.into()],
                ("Temps_after_losses0".into(), Some("degC".into())) => vec![38.82044560943649.into(), 38.82044560943607.into()],
                ("Temps_after_losses1".into(), Some("degC".into())) => vec![38.820445609439716.into(), 38.82044560943589.into()],
                ("Temps_after_losses2".into(), Some("degC".into())) => vec![38.82044560954896.into(), 38.8204456094357.into()],
                ("Temps_after_losses3".into(), Some("degC".into())) => vec![38.82044561251537.into(), 38.820445609433904.into()],
                ("Temps_after_losses4".into(), Some("degC".into())) => vec![38.82044568377406.into(), 38.82044560942411.into()],
                ("Temps_after_losses5".into(), Some("degC".into())) => vec![38.82044727208915.into(), 38.820445609382055.into()],
                ("Temps_after_losses6".into(), Some("degC".into())) => vec![38.82048124427147.into(), 38.8204456092257.into()],
                ("Temps_after_losses7".into(), Some("degC".into())) => vec![38.82118099093947.into(), 38.820445608708134.into()],
                ("total_charge".into(), Some("kWh".into())) => vec![0.0.into(); 2],
                ("end_of_timestep_charge".into(), Some("kWh".into())) => vec![0.0.into(); 2],
                ("hb_after_only_charge_zone_temp0".into(), Some("degC".into())) => vec![38.82044560943649.into(), 38.82044560943607.into()],
                ("hb_after_only_charge_zone_temp1".into(), Some("degC".into())) => vec![38.820445609439716.into(), 38.82044560943589.into()],
                ("hb_after_only_charge_zone_temp2".into(), Some("degC".into())) => vec![38.82044560954896.into(), 38.8204456094357.into()],
                ("hb_after_only_charge_zone_temp3".into(), Some("degC".into())) => vec![38.82044561251537.into(), 38.820445609433904.into()],
                ("hb_after_only_charge_zone_temp4".into(), Some("degC".into())) => vec![38.82044568377406.into(), 38.82044560942411.into()],
                ("hb_after_only_charge_zone_temp5".into(), Some("degC".into())) => vec![38.82044727208915.into(), 38.820445609382055.into()],
                ("hb_after_only_charge_zone_temp6".into(), Some("degC".into())) => vec![38.82048124427147.into(), 38.8204456092257.into()],
                ("hb_after_only_charge_zone_temp7".into(), Some("degC".into())) => vec![38.82118099093947.into(), 38.820445608708134.into()],
            },
            "new_service".into() => indexmap! {
                ("service_name".into(), None) => vec![ResultParamValue::String("new_service".into()); 2],
                ("service_type".into(), None) => vec![ResultParamValue::String(HeatingServiceType::DomesticHotWaterRegular.to_string().into()); 2],
                ("service_on".into(), None) => vec![ResultParamValue::Boolean(true); 2],
                ("energy_output_required".into(), Some("kWh".into())) => vec![100.0.into(); 2],
                ("temp_output".into(), Some("degC".into())) => vec![40.000471231805946.into(), 40.0.into()],
                ("temp_inlet".into(), Some("degC".into())) => vec![40.0.into(); 2],
                ("time_running".into(), Some("secs".into())) => vec![3600.0.into(); 2],
                ("energy_delivered_HB".into(), Some("kWh".into())) => vec![10.509408477594047.into(), (-0.09999181091672799).into()],
                ("energy_delivered_backup".into(), Some("kWh".into())) => vec![0.0.into(); 2],
                ("energy_delivered_total".into(), Some("kWh".into())) => vec![10.509408477594047.into(), (-0.09999181091672799).into()],
                ("energy_charged_during_service".into(), Some("kWh".into())) => vec![0.0.into(); 2],
                ("hb_zone_temperatures0".into(), Some("degC".into())) => vec![
                    40.00000000000006.into(),
                    39.99999999999964.into(),
                ],
                ("hb_zone_temperatures1".into(), Some("degC".into())) => vec![40.00000000000328.into(), 40.0.into()],
                ("hb_zone_temperatures2".into(), Some("degC".into())) => vec![40.00000000011253.into(), 40.0.into()],
                ("hb_zone_temperatures3".into(), Some("degC".into())) => vec![40.00000000307894.into(), 40.0.into()],
                ("hb_zone_temperatures4".into(), Some("degC".into())) => vec![40.00000007433763.into(), 40.0.into()],
                ("hb_zone_temperatures5".into(), Some("degC".into())) => vec![40.000001662652714.into(), 40.0.into()],
                ("hb_zone_temperatures6".into(), Some("degC".into())) => vec![40.00003563483504.into(), 40.0.into()],
                ("hb_zone_temperatures7".into(), Some("degC".into())) => vec![40.000735381503034.into(), 40.0.into()],
                ("current_hb_power".into(), Some("kW".into())) => vec![10.509408477594047.into(), (-0.09999181091672799).into()],
                ("energy_delivered_H4".into(), Some("kWh".into())) => vec![100.0.into()],
            },
        };

        let mut expected_results_annual: ResultsAnnual = indexmap! {
            "Overall".into() => indexmap! {
                ("energy_output_required".into(), Some("kWh".into())) => 200.0.into(),
                ("time_running".into(), Some("secs".into())) => 7200.0.into(),
                ("energy_delivered_HB".into(), Some("kWh".into())) => 10.409416666677318.into(),
                ("energy_delivered_backup".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_total".into(), Some("kWh".into())) => 10.409416666677318.into(),
                ("energy_charged_during_service".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_H4".into(), Some("kWh".into())) => 100.0.into(),
            },
            "auxiliary".into() => indexmap! {
                ("energy_aux".into(), Some("kWh".into())) => 0.12.into(),
                ("battery_losses".into(), Some("kWh".into())) => 0.2.into(),
                ("total_charge".into(), Some("kWh".into())) => 0.0.into(),
                ("end_of_timestep_charge".into(), Some("kWh".into())) => 0.0.into(),
            },
            "new_service".into() => indexmap! {
                ("energy_output_required".into(), Some("kWh".into())) => 200.0.into(),
                ("time_running".into(), Some("secs".into())) => 7200.0.into(),
                ("energy_delivered_HB".into(), Some("kWh".into())) => 10.409416666677318.into(),
                ("energy_delivered_backup".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_total".into(), Some("kWh".into())) => 10.409416666677318.into(),
                ("energy_charged_during_service".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_H4".into(), Some("kWh".into())) => 100.0.into(),
            },
        };

        for (t_idx, _) in simulation_time.iter().enumerate() {
            heat_battery
                .read()
                .demand_energy(
                    service_name,
                    HeatingServiceType::DomesticHotWaterRegular,
                    100.,
                    40.,
                    Some(50.),
                    true,
                    None,
                    Some(true),
                )
                .unwrap();

            heat_battery.read().timestep_end(t_idx).unwrap();
        }

        let (results_per_timestep, results_annual) = heat_battery
            .read()
            .output_detailed_results(
                &indexmap! { "hwsname".into() => vec![100.0.into()] },
                &indexmap! { service_name.into() => "hwsname".into()},
            )
            .unwrap();

        assert_eq!(
            results_per_timestep.keys().collect_vec(),
            expected_results_per_timestep.keys().collect_vec()
        );
        assert_eq!(
            results_annual.keys().collect_vec(),
            expected_results_annual.keys().collect_vec()
        );

        let assert_value =
            |actual: &ResultParamValue, expected: &ResultParamValue| match (actual, expected) {
                (ResultParamValue::Number(actual_num), ResultParamValue::Number(expected_num)) => {
                    assert_relative_eq!(actual_num, expected_num, max_relative = 1e-7);
                }
                _ => assert_eq!(actual, expected,),
            };

        for (key, expected_results) in &expected_results_per_timestep {
            let actual_results = &results_per_timestep[key];

            assert_eq!(
                actual_results.keys().collect_vec(),
                expected_results.keys().collect_vec()
            );

            for (inner_key, expected_vec) in expected_results {
                for (actual, expected) in actual_results[inner_key].iter().zip(expected_vec) {
                    assert_value(actual, expected);
                }
            }
        }

        for (key, expected_results) in &expected_results_annual {
            let actual_results = &results_annual[key];

            assert_eq!(
                actual_results.keys().collect_vec(),
                expected_results.keys().collect_vec()
            );

            for (inner_key, value) in expected_results {
                assert_value(&actual_results[inner_key], value);
            }
        }

        // Test case where hot water source is not in hot water energy source data
        expected_results_per_timestep[service_name].insert(
            ("energy_delivered_H4".into(), Some("kWh".into())),
            vec![ResultParamValue::Empty; 2],
        );
        expected_results_annual[service_name].insert(
            ("energy_delivered_H4".into(), Some("kWh".into())),
            ResultParamValue::Empty,
        );
        expected_results_annual["Overall"].insert(
            ("energy_delivered_H4".into(), Some("kWh".into())),
            ResultParamValue::Empty,
        );

        let (results_per_timestep, results_annual) = heat_battery
            .read()
            .output_detailed_results(
                &indexmap! { "hwsname".into() => vec![100.0.into()] },
                &indexmap! { service_name.into() => "hwsname_other".into()},
            )
            .unwrap();

        assert_eq!(
            results_per_timestep.keys().collect_vec(),
            expected_results_per_timestep.keys().collect_vec()
        );
        assert_eq!(
            results_annual.keys().collect_vec(),
            expected_results_annual.keys().collect_vec()
        );

        for (key, expected_results) in &expected_results_per_timestep {
            let actual_results = &results_per_timestep[key];

            assert_eq!(
                actual_results.keys().collect_vec(),
                expected_results.keys().collect_vec()
            );

            for (inner_key, expected_vec) in expected_results {
                for (actual, expected) in actual_results[inner_key].iter().zip(expected_vec) {
                    assert_value(actual, expected);
                }
            }
        }

        for (key, expected_results) in &expected_results_annual {
            let actual_results = &results_annual[key];

            assert_eq!(
                actual_results.keys().collect_vec(),
                expected_results.keys().collect_vec()
            );

            for (inner_key, value) in expected_results {
                assert_value(&actual_results[inner_key], value);
            }
        }
    }

    #[rstest]
    #[ignore = "test to be updated during migration to 1.0.0a6"]
    fn test_output_detailed_results_space(
        simulation_time: SimulationTime,
        heat_battery_no_service_connection: Arc<RwLock<HeatBatteryPcm>>,
    ) {
        let heat_battery = heat_battery_no_service_connection;
        let service_name = "new_service";
        let control = create_setpoint_time_control(vec![]);

        HeatBatteryPcm::create_service_space_heating(
            heat_battery.clone(),
            service_name,
            Arc::new(control),
        )
        .unwrap();

        let expected_results_per_timestep: ResultsPerTimestep = indexmap! {
            "auxiliary".into() => indexmap! {
                ("energy_aux".into(), Some("kWh".into())) => vec![0.06.into(), 0.06.into()],
                ("battery_losses".into(), Some("kWh".into())) => vec![0.1.into(), 0.1.into()],
                ("Temps_after_losses0".into(), Some("degC".into())) => vec![38.82044560943649.into(), 38.82044560943607.into()],
                ("Temps_after_losses1".into(), Some("degC".into())) => vec![38.820445609439716.into(), 38.82044560943589.into()],
                ("Temps_after_losses2".into(), Some("degC".into())) => vec![38.82044560954896.into(), 38.8204456094357.into()],
                ("Temps_after_losses3".into(), Some("degC".into())) => vec![38.82044561251537.into(), 38.820445609433904.into()],
                ("Temps_after_losses4".into(), Some("degC".into())) => vec![38.82044568377406.into(), 38.82044560942411.into()],
                ("Temps_after_losses5".into(), Some("degC".into())) => vec![38.82044727208915.into(), 38.820445609382055.into()],
                ("Temps_after_losses6".into(), Some("degC".into())) => vec![38.82048124427147.into(), 38.8204456092257.into()],
                ("Temps_after_losses7".into(), Some("degC".into())) => vec![38.82118099093947.into(), 38.820445608708134.into()],
                ("total_charge".into(), Some("kWh".into())) => vec![0.0.into(), 0.0.into()],
                ("end_of_timestep_charge".into(), Some("kWh".into())) => vec![0.0.into(), 0.0.into()],
                ("hb_after_only_charge_zone_temp0".into(), Some("degC".into())) => vec![38.82044560943649.into(), 38.82044560943607.into()],
                ("hb_after_only_charge_zone_temp1".into(), Some("degC".into())) => vec![38.820445609439716.into(), 38.82044560943589.into()],
                ("hb_after_only_charge_zone_temp2".into(), Some("degC".into())) => vec![38.82044560954896.into(), 38.8204456094357.into()],
                ("hb_after_only_charge_zone_temp3".into(), Some("degC".into())) => vec![38.82044561251537.into(), 38.820445609433904.into()],
                ("hb_after_only_charge_zone_temp4".into(), Some("degC".into())) => vec![38.82044568377406.into(), 38.82044560942411.into()],
                ("hb_after_only_charge_zone_temp5".into(), Some("degC".into())) => vec![38.82044727208915.into(), 38.820445609382055.into()],
                ("hb_after_only_charge_zone_temp6".into(), Some("degC".into())) => vec![38.82048124427147.into(), 38.8204456092257.into()],
                ("hb_after_only_charge_zone_temp7".into(), Some("degC".into())) => vec![38.82118099093947.into(), 38.820445608708134.into()],
            },
            "new_service".into() => indexmap! {
                ("service_name".into(), None) => vec![ResultParamValue::String("new_service".into()), ResultParamValue::String("new_service".into())],
                ("service_type".into(), None) => vec![ResultParamValue::String(HeatingServiceType::Space.to_string().into()), ResultParamValue::String(HeatingServiceType::Space.to_string().into())],
                ("service_on".into(), None) => vec![ResultParamValue::Boolean(true), ResultParamValue::Boolean(true)],
                ("energy_output_required".into(), Some("kWh".into())) => vec![100.0.into(), 100.0.into()],
                ("temp_output".into(), Some("degC".into())) => vec![40.000471231805946.into(), 39.99999999956042.into()],
                ("temp_inlet".into(), Some("degC".into())) => vec![40.0.into(), 40.0.into()],
                ("time_running".into(), Some("secs".into())) => vec![3600.0.into(), 3600.0.into()],
                ("energy_delivered_HB".into(), Some("kWh".into())) => vec![10.509408477594047.into(), (-0.09999181091672799).into()],
                ("energy_delivered_backup".into(), Some("kWh".into())) => vec![0.0.into(), 0.0.into()],
                ("energy_delivered_total".into(), Some("kWh".into())) => vec![10.509408477594047.into(), (-0.09999181091672799).into()],
                ("energy_charged_during_service".into(), Some("kWh".into())) => vec![0.0.into(), 0.0.into()],
                ("hb_zone_temperatures0".into(), Some("degC".into())) => vec![40.00000000000006.into(), 39.99999999999964.into()],
                ("hb_zone_temperatures1".into(), Some("degC".into())) => vec![40.00000000000328.into(), 39.99999999999946.into()],
                ("hb_zone_temperatures2".into(), Some("degC".into())) => vec![40.00000000011253.into(), 39.99999999999927.into()],
                ("hb_zone_temperatures3".into(), Some("degC".into())) => vec![40.00000000307894.into(), 39.99999999999747.into()],
                ("hb_zone_temperatures4".into(), Some("degC".into())) => vec![40.00000007433763.into(), 39.99999999998768.into()],
                ("hb_zone_temperatures5".into(), Some("degC".into())) => vec![40.000001662652714.into(), 39.99999999994562.into()],
                ("hb_zone_temperatures6".into(), Some("degC".into())) => vec![40.00003563483504.into(), 39.99999999978927.into()],
                ("hb_zone_temperatures7".into(), Some("degC".into())) => vec![40.000735381503034.into(), 39.9999999992717.into()],
                ("current_hb_power".into(), Some("kW".into())) => vec![10.509408477594047.into(), (-0.09999181091672799).into()],
                ("energy_delivered_H4".into(), Some("kWh".into())) => vec![10.509408477594047.into(), (-0.09999181091672799).into()],
            },
        };

        let expected_results_annual: ResultsAnnual = indexmap! {
            "Overall".into() => indexmap! {
                ("energy_output_required".into(), Some("kWh".into())) => 200.0.into(),
                ("time_running".into(), Some("secs".into())) => 7200.0.into(),
                ("energy_delivered_HB".into(), Some("kWh".into())) => 10.409416666677318.into(),
                ("energy_delivered_backup".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_total".into(), Some("kWh".into())) => 10.409416666677318.into(),
                ("energy_charged_during_service".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_H4".into(), Some("kWh".into())) => 10.409416666677318.into(),
            },
            "auxiliary".into() => indexmap! {
                ("energy_aux".into(), Some("kWh".into())) => 0.12.into(),
                ("battery_losses".into(), Some("kWh".into())) => 0.2.into(),
                ("total_charge".into(), Some("kWh".into())) => 0.0.into(),
                ("end_of_timestep_charge".into(), Some("kWh".into())) => 0.0.into(),
            },
            "new_service".into() => indexmap! {
                ("energy_output_required".into(), Some("kWh".into())) => 200.0.into(),
                ("time_running".into(), Some("secs".into())) => 7200.0.into(),
                ("energy_delivered_HB".into(), Some("kWh".into())) => 10.409416666677318.into(),
                ("energy_delivered_backup".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_total".into(), Some("kWh".into())) => 10.409416666677318.into(),
                ("energy_charged_during_service".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_H4".into(), Some("kWh".into())) => 10.409416666677318.into(),
            },
        };

        for (t_idx, _) in simulation_time.iter().enumerate() {
            heat_battery
                .read()
                .demand_energy(
                    service_name,
                    HeatingServiceType::Space,
                    100.,
                    40.,
                    Some(55.),
                    true,
                    None,
                    Some(true),
                )
                .unwrap();

            heat_battery.read().timestep_end(t_idx).unwrap();
        }

        let (results_per_timestep, results_annual) = heat_battery
            .read()
            .output_detailed_results(&indexmap! {}, &indexmap! {})
            .unwrap();

        assert_eq!(
            results_per_timestep.keys().collect_vec(),
            expected_results_per_timestep.keys().collect_vec()
        );
        assert_eq!(
            results_annual.keys().collect_vec(),
            expected_results_annual.keys().collect_vec()
        );

        let assert_value =
            |actual: &ResultParamValue, expected: &ResultParamValue| match (actual, expected) {
                (ResultParamValue::Number(actual_num), ResultParamValue::Number(expected_num)) => {
                    assert_relative_eq!(actual_num, expected_num, max_relative = 1e-7);
                }
                _ => assert_eq!(actual, expected,),
            };

        for (key, expected_results) in &expected_results_per_timestep {
            let actual_results = &results_per_timestep[key];

            assert_eq!(
                actual_results.keys().collect_vec(),
                expected_results.keys().collect_vec()
            );

            for (inner_key, expected_vec) in expected_results {
                for (actual, expected) in actual_results[inner_key].iter().zip(expected_vec) {
                    assert_value(actual, expected);
                }
            }
        }

        for (key, expected_results) in &expected_results_annual {
            let actual_results = &results_annual[key];

            assert_eq!(
                actual_results.keys().collect_vec(),
                expected_results.keys().collect_vec()
            );

            for (inner_key, value) in expected_results {
                assert_value(&actual_results[inner_key], value);
            }
        }
    }

    #[rstest]
    fn test_output_detailed_results_none(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
        // Test that calling output_detailed_results errors when output_detailed_results on heat_battery is false
        let heat_battery =
            create_heat_battery(simulation_time_iterator, battery_control_on, Some(false));

        assert!(heat_battery
            .read()
            .output_detailed_results(&indexmap! {}, &indexmap! {})
            .is_err());
    }

    #[rstest]
    #[ignore = "test to be updated during migration to 1.0.0a6"]
    fn test_demand_energy_low_temp_minimum_run_coverage(
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);
        heat_battery
            .write()
            .energy_supply
            .write()
            .set_fuel_type(FuelType::MainsGas);
        heat_battery.write().hb_time_step = 5.; // Small time step
                                                // Set all zones to high temperature
        heat_battery.write().zone_temp_c_dist_initial = Arc::new(RwLock::new(vec![50.2; 8]));
        HeatBatteryPcm::create_service_connection(heat_battery.clone(), "test_service").unwrap();
        // Very small energy demand that will be satisfied in first loop iteration
        // But will need to continue running to meet minimum time
        heat_battery
            .read()
            .demand_energy(
                "test_service",
                HeatingServiceType::DomesticHotWaterRegular,
                0.1,
                40.,
                Some(50.),
                true,
                Some(0.),
                Some(true),
            )
            .unwrap();

        //Check that minimum time was enforced
        let service_result = heat_battery
            .read()
            .service_results
            .read()
            .last()
            .unwrap()
            .clone();

        assert_relative_eq!(service_result.time_running, 51.08689856959955);
    }

    #[rstest]
    fn test_timestep_end_with_uncalled_services(
        heat_battery_no_service_connection: Arc<RwLock<HeatBatteryPcm>>,
    ) {
        let heat_battery = heat_battery_no_service_connection;

        // Create three services
        let service1 = "water_heating";
        let service2 = "space_heating_zone1";
        let service3 = "space_heating_zone2";

        HeatBatteryPcm::create_service_connection(heat_battery.clone(), service1).unwrap();
        HeatBatteryPcm::create_service_connection(heat_battery.clone(), service2).unwrap();
        HeatBatteryPcm::create_service_connection(heat_battery.clone(), service3).unwrap();

        // In timestep 1: Call only service1 and service3 (skip service2)
        heat_battery
            .read()
            .demand_energy(
                service1,
                HeatingServiceType::DomesticHotWaterRegular,
                5.,
                40.,
                Some(55.),
                true,
                Some(0.),
                Some(true),
            )
            .unwrap();

        heat_battery
            .read()
            .demand_energy(
                service3,
                HeatingServiceType::Space,
                3.,
                35.,
                Some(50.),
                true,
                Some(0.),
                Some(true),
            )
            .unwrap();

        heat_battery.read().timestep_end(0).unwrap();

        {
            let hb_guard = heat_battery.read();
            let detailed_results_guard = hb_guard.detailed_results.as_ref().unwrap().read();

            // Check that detailed results were created
            assert_eq!(detailed_results_guard.len(), 1);

            let timestep_results = &detailed_results_guard[0].results;

            assert_eq!(timestep_results.len(), 3); // In Python this is 4 (Should have 3 service results + 1 auxiliary result = 4 total)

            // Check service1 (was called)
            assert_eq!(timestep_results[0].service_name, service1);
            assert_eq!(
                timestep_results[0].service_type.unwrap(),
                HeatingServiceType::DomesticHotWaterRegular
            );
            assert!(timestep_results[0].service_on);
            assert!(timestep_results[0].time_running > 0.);

            // Check service2 (was NOT called - should have placeholder values)
            assert_eq!(timestep_results[1].service_name, service2);
            assert!(timestep_results[1].service_type.is_none());
            assert!(!timestep_results[1].service_on);
            assert_eq!(timestep_results[1].energy_output_required, 0.);
            assert_eq!(timestep_results[1].time_running, 0.);
            assert_eq!(timestep_results[1].energy_delivered_hb, 0.);
            assert_eq!(timestep_results[1].current_hb_power, 0.);

            // Check service3 (was called)
            assert_eq!(timestep_results[2].service_name, service3);
            assert_eq!(
                timestep_results[2].service_type.unwrap(),
                HeatingServiceType::Space
            );
            assert!(timestep_results[2].service_on);
            assert!(timestep_results[2].time_running > 0.);

            // Check auxiliary results
            let summary = &detailed_results_guard[0].summary;
            assert!(summary.energy_aux >= 0.);
            assert!(summary.battery_losses >= 0.);
            assert!(!summary.temps_after_losses.is_empty());
            assert!(summary.total_charge >= 0.);
            assert!(summary.end_of_timestep_charge >= 0.);
            assert!(!summary.hb_after_only_charge_zone_temp.is_empty());
        }

        // In timestep 2: Call only service2 (skip service1 and service3)
        heat_battery
            .read()
            .demand_energy(
                service2,
                HeatingServiceType::Space,
                4.,
                38.,
                Some(52.),
                true,
                Some(0.),
                Some(true),
            )
            .unwrap();

        heat_battery.read().timestep_end(1).unwrap();

        {
            let hb_guard = heat_battery.read();
            let detailed_results_guard = hb_guard.detailed_results.as_ref().unwrap().read();

            // Check second timestep results
            assert_eq!(detailed_results_guard.len(), 2);

            let timestep2_results = &detailed_results_guard[1].results;

            // service1 should have placeholder values this time
            assert_eq!(timestep2_results[0].service_name, service1);
            assert!(timestep2_results[0].service_type.is_none());
            assert!(!timestep2_results[0].service_on);
            assert_eq!(timestep2_results[0].time_running, 0.);

            // service2 should have actual values
            assert_eq!(timestep2_results[1].service_name, service2);
            assert_eq!(
                timestep2_results[1].service_type.unwrap(),
                HeatingServiceType::Space
            );
            assert!(timestep2_results[1].service_on);
            assert!(timestep2_results[1].time_running > 0.);

            // service3 should have placeholder values
            assert_eq!(timestep2_results[2].service_name, service3);
            assert!(timestep2_results[2].service_type.is_none());
            assert!(!timestep2_results[2].service_on);
            assert_eq!(timestep2_results[2].time_running, 0.);
        }
    }

    #[rstest]
    fn test_timestep_end_no_services_called(
        heat_battery_no_service_connection: Arc<RwLock<HeatBatteryPcm>>,
    ) {
        // Test timestep_end when no services are called but services are registered
        let heat_battery = heat_battery_no_service_connection;

        // Create services but don't call them
        let service1 = "water_heating";
        let service2 = "space_heating";

        HeatBatteryPcm::create_service_connection(heat_battery.clone(), service1).unwrap();
        HeatBatteryPcm::create_service_connection(heat_battery.clone(), service2).unwrap();

        // Call timestep_end without calling any services
        heat_battery.read().timestep_end(0).unwrap();

        let hb_guard = heat_battery.read();
        let detailed_results_guard = hb_guard.detailed_results.as_ref().unwrap().read();

        //  Check that detailed results were created with placeholder entries
        assert_eq!(detailed_results_guard.len(), 1);

        let timestep_results = &detailed_results_guard[0].results;

        assert_eq!(timestep_results.len(), 2); // In Python this is 3 (Should have 2 service results + 1 auxiliary result = 3 total)

        for (i, result) in timestep_results.iter().enumerate() {
            assert_eq!(result.service_name, [service1, service2][i]);
            assert!(result.service_type.is_none());
            assert!(!result.service_on);
            assert_eq!(result.energy_output_required, 0.);
            assert_eq!(result.time_running, 0.);
            assert_eq!(result.energy_delivered_hb, 0.);
            assert_eq!(result.current_hb_power, 0.);
        }
    }

    #[rstest]
    fn test_heat_battery_create_service_connection_already_exists(
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);
        let service_name = "test_service";
        let cold_feed =
            WaterSupply::ColdWaterSource(Arc::new(ColdWaterSource::new(vec![1.0, 1.2], 0, 1.)));

        let result = HeatBatteryPcm::create_service_hot_water_regular(
            heat_battery.clone(),
            service_name,
            cold_feed.clone(),
            Arc::new(Control::Mock(MockControl::default())),
            Arc::new(Control::Mock(MockControl::default())),
        );

        assert!(result.is_ok());

        let result = HeatBatteryPcm::create_service_hot_water_regular(
            heat_battery,
            service_name,
            cold_feed,
            Arc::new(Control::Mock(MockControl::default())),
            Arc::new(Control::Mock(MockControl::default())),
        );

        assert!(result.is_err())
    }

    // skipping python's test_heat_battery_edge_case_zero_timestep as function can't return None

    // skipping python's test_heat_battery_process_zone_edge_cases as function can't return None and does return 4 values

    // skipping python's test_heat_battery_charge_battery_hydraulic_edge_cases as function can't return None

    #[rstest]
    fn test_heat_battery_energy_output_max_boundary_conditions(
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);

        // Test with very low output temperature
        let result = heat_battery
            .read()
            .energy_output_max(10., Some(0.))
            .unwrap();

        // The method returns energy based on zone temps, not necessarily 0
        assert!(result >= 0.);

        //Test with temperature at threshold
        let result = heat_battery
            .read()
            .energy_output_max(45., Some(0.))
            .unwrap();

        assert!(result >= 0.);
    }

    // skipping python's test_heat_battery_service_cold_water_source_not_set as cold feed not optional

    #[rstest]
    fn test_heat_battery_zero_volume_zones(
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);

        // Test with zero energy transfer
        let result = heat_battery.read().calculate_new_zone_temperature(50., 0.);

        assert_eq!(result, 50.);
    }

    #[rstest]
    fn test_heat_battery_all_zones_below_threshold(
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        // Request high output temperature that no zone can provide
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off, None);

        let result = heat_battery
            .read()
            .energy_output_max(80., Some(0.))
            .unwrap();

        assert_eq!(result, 0.);
    }

    // skipping python's test_demand_hot_water_with_varying_cold_temperatures and
    // test_demand_hot_water_zero_volume_continue due to mocking
}
