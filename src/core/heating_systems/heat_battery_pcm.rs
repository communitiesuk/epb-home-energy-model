/// This module provides object(s) to model the behaviour of heat batteries.
use crate::core::common::WaterSourceWithTemperature;
use crate::core::controls::time_control::{per_control, Control, ControlBehaviour};
use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::core::heating_systems::common::HeatingServiceType;
use crate::core::material_properties::WATER;
use crate::core::schedule::TypedScheduleEvent;
use crate::core::units::{
    KILOJOULES_PER_KILOWATT_HOUR, SECONDS_PER_HOUR, SECONDS_PER_MINUTE, WATTS_PER_KILOWATT,
};
use crate::corpus::{
    ResultParamValue, ResultsAnnual as CorpusResultsAnnual,
    ResultsPerTimestep as CorpusResultsPerTimestep,
};
use crate::input::{HeatBattery as HeatBatteryInput, HeatSourceWetDetails};
use crate::simulation_time::{SimulationTimeIteration, SimulationTimeIterator};
use crate::StringOrNumber;
use anyhow::anyhow;
use anyhow::bail;
use atomic_float::AtomicF64;
use indexmap::IndexMap;
use itertools::Itertools;
use parking_lot::RwLock;
use smartstring::alias::String;
use std::collections::HashMap;
use std::ops::{Deref, DerefMut};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

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
    cold_feed: WaterSourceWithTemperature,
    control: Option<Arc<Control>>,
    _control_min: Option<Arc<Control>>,
    control_max: Option<Arc<Control>>,
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
        cold_feed: WaterSourceWithTemperature,
        control_min: Option<Arc<Control>>,
        control_max: Option<Arc<Control>>,
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
    ) -> anyhow::Result<(Option<f64>, Option<f64>)> {
        Ok((
            self._control_min
                .as_ref()
                .ok_or_else(|| {
                    anyhow!("control min expected on heat battery when setpnt is called")
                })?
                .setpnt(&simulation_time_iteration),
            self.control_max
                .as_ref()
                .ok_or_else(|| {
                    anyhow!("control max expected on heat battery when setpnt is called")
                })?
                .setpnt(&simulation_time_iteration),
        ))
    }

    // TODO remove below 3 methods 1.0.0a1
    pub(crate) fn get_cold_water_source(&self) -> &WaterSourceWithTemperature {
        &self.cold_feed
    }

    pub(crate) fn get_temp_hot_water(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let volume = 20.; // Nominal volumen to calculate water temperature from battery
        let inlet_temp = self.cold_feed.temperature(simulation_time_iteration, None);

        self.heat_battery
            .read()
            .get_temp_hot_water(inlet_temp, volume)
    }

    pub(crate) fn demand_hot_water(
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

        let service_on = self.is_on(simulation_time_iteration)?;

        if !service_on {
            energy_demand = 0.;
        }

        self.heat_battery.read().demand_energy(
            &self.service_name,
            HeatingServiceType::DomesticHotWaterRegular,
            energy_demand,
            self.cold_feed.temperature(simulation_time_iteration, None),
            Some(temp_hot_water),
            service_on,
            None,
            Some(true),
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
        let service_on = self.is_on(simulation_time_iteration)?;
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
        let service_on = self.is_on(simulation_time_iteration)?;
        if !service_on {
            return Ok(0.);
        }

        self.heat_battery.read().energy_output_max(temp_flow, None)
    }

    fn is_on(&self, simulation_time_iteration: SimulationTimeIteration) -> anyhow::Result<bool> {
        if let Some(control) = &self.control {
            Ok(control.is_on(simulation_time_iteration))
        } else {
            Ok(true)
        }
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

const DEFAULT_N_ZONES: usize = 8;
const DEFAULT_HB_TIME_STEP: f64 = 20.;
const DEFAULT_MINIMUM_TIME_REQUIRED_TO_RUN: f64 = 120.;
const DEFAULT_INITIAL_INLET_TEMP: f64 = 10.;
const DEFAULT_ESTIMATED_OUTLET_TEMP: f64 = 53.;

// nothing seems to read this - check upstream whether service_results field is necessary
#[derive(Clone, Debug)]
#[allow(dead_code)]
struct HeatBatteryResult {
    service_name: String,
    service_type: HeatingServiceType,
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

impl HeatBatteryResult {
    fn param_value_as_string(&self, param: &str) -> String {
        match param {
            "service_name" => self.service_name.clone(),
            "service_type" => format!("{:?}", self.service_type).into(),
            "service_on" => self.service_on.to_string().into(),
            "energy_output_required" => self.energy_output_required.to_string().into(),
            "temp_output" => match self.temp_output {
                Some(temp) => temp.to_string().into(),
                None => "".into(),
            },
            "temp_inlet" => self.temp_inlet.to_string().into(),
            "time_running" => self.time_running.to_string().into(),
            "energy_left_in_pipe" => self.energy_left_in_pipe.to_string().into(),
            "temperature_left_in_pipe" => self.temperature_left_in_pipe.to_string().into(),
            "energy_delivered_hb" => self.energy_delivered_hb.to_string().into(),
            "energy_delivered_backup" => self.energy_delivered_backup.to_string().into(),
            "energy_delivered_total" => self.energy_delivered_total.to_string().into(),
            "energy_delivered_low_temp" => self.energy_delivered_low_temp.to_string().into(),
            "energy_charged_during_service" => {
                self.energy_charged_during_service.to_string().into()
            }
            "hb_zone_temperatures" => self.hb_zone_temperatures.iter().join(",").into(),
            "current_hb_power" => self.current_hb_power.to_string().into(),
            _ => panic!("Unknown parameter: {}", param),
        }
    }
}

#[derive(Debug)]
struct HeatBatteryTimestepSummary {
    energy_aux: f64,
    battery_losses: f64,
    temps_after_losses: Vec<f64>,
    total_charge: f64,
    end_of_timestep_charge: f64,
    hb_after_only_charge_zone_temp: Vec<f64>,
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
    energy_supply_connections: HashMap<String, EnergySupplyConnection>,
    pwr_in: f64,
    max_rated_losses: f64,
    power_circ_pump: f64,
    power_standby: f64,
    n_units: usize,
    charge_control: Arc<Control>, // ChargeControl variant expected
    // nothing external seems to read this - check upstream whether service_results field is necessary
    service_results: Arc<RwLock<Vec<HeatBatteryResult>>>,
    total_time_running_current_timestep: AtomicF64,
    flag_first_call: AtomicBool,
    #[allow(dead_code)]
    charge_level: f64,
    n_zones: usize,
    hb_time_step: f64,
    minimum_time_required_to_run: f64,
    initial_inlet_temp: f64,
    estimated_outlet_temp: f64,
    pipe_energy: Arc<RwLock<IndexMap<String, PipeEnergy>>>,
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
    heat_exchanger_surface_area_m2: f64,
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
        n_zones: Option<usize>,
        hb_time_step: Option<f64>,
        minimum_time_required_to_run: Option<f64>,
        initial_inlet_temp: Option<f64>,
        estimated_outlet_temp: Option<f64>,
        output_detailed_results: Option<bool>,
    ) -> Self {
        let n_zones = n_zones.unwrap_or(DEFAULT_N_ZONES);
        let hb_time_step = hb_time_step.unwrap_or(DEFAULT_HB_TIME_STEP);
        let minimum_time_required_to_run =
            minimum_time_required_to_run.unwrap_or(DEFAULT_MINIMUM_TIME_REQUIRED_TO_RUN);
        let initial_inlet_temp = initial_inlet_temp.unwrap_or(DEFAULT_INITIAL_INLET_TEMP);
        let estimated_outlet_temp = estimated_outlet_temp.unwrap_or(DEFAULT_ESTIMATED_OUTLET_TEMP);
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
                    capillary_diameter_m,
                    a,
                    b,
                    heat_exchanger_surface_area_m2,
                    flow_rate_l_per_min,
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
            flag_first_call: true.into(),
            charge_level: Default::default(),
            n_zones,
            hb_time_step,
            minimum_time_required_to_run,
            initial_inlet_temp,
            estimated_outlet_temp,
            pipe_energy: Default::default(),
            energy_charged: Default::default(),
            simultaneous_charging_and_discharging,
            max_temp_of_charge,
            zone_temp_c_dist_initial: Arc::new(RwLock::new(vec![max_temp_of_charge; n_zones])),
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

        // Set up PipeEnergy for this service to store extra
        // energy pushed into the pipe to run the battery and temperature
        heat_battery.read().pipe_energy.write().insert(
            service_name.into(),
            PipeEnergy {
                energy: 0.,
                temperature: 0.,
            },
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
        cold_feed: WaterSourceWithTemperature,
        control_min: Option<Arc<Control>>,
        control_max: Option<Arc<Control>>,
    ) -> HeatBatteryPcmServiceWaterRegular {
        Self::create_service_connection(heat_battery.clone(), service_name).unwrap();
        HeatBatteryPcmServiceWaterRegular::new(
            heat_battery,
            service_name.into(),
            cold_feed,
            control_min,
            control_max,
        )
    }

    pub(crate) fn create_service_space_heating(
        heat_battery: Arc<RwLock<Self>>,
        service_name: &str,
        control: Arc<Control>,
    ) -> HeatBatteryPcmServiceSpace {
        Self::create_service_connection(heat_battery.clone(), service_name).unwrap();
        HeatBatteryPcmServiceSpace::new(heat_battery, service_name.into(), control)
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

        let mut energy_transf_delivered = vec![0.; self.n_zones];
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

            let energy_charged_during_battery_time_step = energy_transf_charged.iter().sum::<f64>();

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

        Ok((
            energy_loss.iter().sum::<f64>() / KILOJOULES_PER_KILOWATT_HOUR as f64,
            zone_temp_c_dist,
        ))
    }

    fn get_temp_hot_water(&self, inlet_temp: f64, volume: f64) -> anyhow::Result<f64> {
        let total_time_s = volume / self.flow_rate_l_per_min * SECONDS_PER_MINUTE as f64;

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

        let mut zone_temp_c_dist = self.zone_temp_c_dist_initial.read().clone();
        let mut inlet_temp_c = inlet_temp;

        let n_time_steps = if total_time_s > time_step_s {
            (total_time_s / time_step_s) as usize
        } else {
            1
        };

        let mut outlet_temp_c = inlet_temp_c; // initialise, though expectation is this will be overridden in loop

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

        Ok(outlet_temp_c)
    }

    /// Calculate the maximum energy output of the heat battery, accounting
    /// for time spent on higher-priority services.
    pub(crate) fn energy_output_max(
        &self,
        temp_output: f64,
        _time_start: Option<f64>,
    ) -> anyhow::Result<f64> {
        // Return the energy the battery can provide assuming the HB temperature inlet
        // is constant during HEM time step equal to the required emitter temperature (temp_output)
        // Maximum energy for a given HB zones temperature distribution and inlet temperature.
        // The calculation methodology is the same as described in the demand_energy function.
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

            let energy_delivered_ts: f64 = energy_transf_delivered.iter().sum();

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

    pub(crate) fn demand_energy(
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
        let mut energy_output_required = energy_output_required;

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

        {
            let mut pipe_energy = self.pipe_energy.write();
            let service_pipe_energy = &mut pipe_energy.deref_mut()[service_name];
            if temp_output.is_none() || temp_output.unwrap() <= service_pipe_energy.temperature {
                if energy_output_required > service_pipe_energy.energy {
                    energy_output_required -= service_pipe_energy.energy;
                    service_pipe_energy.energy = 0.;
                    service_pipe_energy.temperature = 0.;
                } else {
                    service_pipe_energy.energy -= energy_output_required;
                    energy_output_required = 0.;
                }
            }
        }

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
        let mut total_energy_low_temp = 0.;
        let inlet_temp_c = temp_return_feed;
        let mut zone_temp_c_dist = self.zone_temp_c_dist_initial.read().clone();

        if energy_output_required <= 0. {
            if update_heat_source_state {
                let pipe_energy = self.pipe_energy.read();
                let service_pipe_energy = &pipe_energy.deref()[service_name];
                self.service_results.write().push(HeatBatteryResult {
                    service_name: service_name.into(),
                    service_type,
                    service_on,
                    energy_output_required,
                    temp_output,
                    temp_inlet: temp_return_feed,
                    time_running: 0.,
                    energy_left_in_pipe: service_pipe_energy.energy,
                    temperature_left_in_pipe: service_pipe_energy.temperature,
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

            if update_heat_source_state {
                self.energy_charged
                    .fetch_add(energy_charged_during_battery_time_step, Ordering::SeqCst);
            }
            energy_charged += energy_charged_during_battery_time_step;

            time_running_current_service += time_step_s;

            // RN for next time step
            water_kinematic_viscosity_m2_per_s =
                Self::calculate_water_kinematic_viscosity_m2_per_s(temp_return_feed, outlet_temp_c);
            reynold_number_at_1_l_per_min = Self::calculate_reynold_number_at_1_l_per_min(
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
                } else if energy_delivered_ts != 0. {
                    let mut pipe_energy = self.pipe_energy.write();
                    let service_pipe_energy = &mut pipe_energy.deref_mut()[service_name];
                    let current_energy = service_pipe_energy.energy;
                    let current_temperature = service_pipe_energy.temperature;
                    let new_temperature = ((current_temperature * current_energy)
                        + (outlet_temp_c * energy_delivered_ts))
                        / (current_energy + energy_delivered_ts);

                    service_pipe_energy.energy += energy_delivered_ts;
                    service_pipe_energy.temperature = new_temperature;
                }

                if time_step_s > self.hb_time_step {
                    time_step_s = self.hb_time_step;
                }

                if (energy_demand - energy_delivered_hb) < 0.0001 {
                    // Energy supplied, run to complete water loop
                    if time_running_current_service > self.minimum_time_required_to_run {
                        break;
                    }

                    if !flag_minimum_run {
                        time_extra =
                            self.minimum_time_required_to_run - time_running_current_service;
                        flag_minimum_run = true;
                    } else {
                        time_extra -= time_step_s;
                    }

                    if time_extra > self.hb_time_step {
                        time_step_s = self.hb_time_step
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
                    if time_running_current_service > self.minimum_time_required_to_run {
                        break;
                    }

                    if !flag_minimum_run {
                        time_extra =
                            self.minimum_time_required_to_run - time_running_current_service;
                        flag_minimum_run = true;
                    } else {
                        time_extra -= time_step_s;
                    }

                    if time_extra > self.hb_time_step {
                        time_step_s = self.hb_time_step
                    } else {
                        time_step_s = time_extra
                    }
                } else {
                    break;
                }
            }
        }

        if update_heat_source_state {
            *self.zone_temp_c_dist_initial.write() = zone_temp_c_dist.clone();

            self.total_time_running_current_timestep.fetch_add(
                time_running_current_service / SECONDS_PER_HOUR as f64,
                Ordering::SeqCst,
            );

            let current_hb_power = if time_running_current_service > 0. {
                energy_delivered_hb * SECONDS_PER_HOUR as f64 / time_running_current_service
            } else {
                Default::default()
            };
            // TODO (from Python) Clarify whether Heat Batteries can have direct electric backup if depleted
            let pipe_energy = self.pipe_energy.read();
            let service_pipe_energy = &pipe_energy.deref()[service_name];
            self.service_results.write().push(HeatBatteryResult {
                service_name: service_name.into(),
                service_type,
                service_on,
                energy_output_required,
                temp_output,
                temp_inlet: temp_return_feed,
                time_running: time_running_current_service,
                energy_left_in_pipe: service_pipe_energy.energy,
                temperature_left_in_pipe: service_pipe_energy.temperature,
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
    ) -> anyhow::Result<f64> {
        // Energy used by circulation pump
        let mut energy_aux = self
            .total_time_running_current_timestep
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

        // Calculating heat battery losses in timestep to correct charge level
        // Currently assumed all losses are to the exterior independently of the
        // heat battery location
        // TODO (from Python): Assign thermal losses to relevant zone if heat battery is not outdoors.
        let (battery_losses, zone_temp_c_after_losses) = self.battery_heat_loss()?;

        // Charging battery for the remaining of the timestep
        let (end_of_ts_charge, zone_temp_c_after_charging) = if self
            .charge_control
            .is_on(self.simulation_time.current_iteration())
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
        if self.detailed_results.is_some() {
            let results = self.service_results.read().clone();
            self.detailed_results
                .as_ref()
                .unwrap()
                .write()
                .push(HeatBatteryTimestepResult {
                    results,
                    summary: HeatBatteryTimestepSummary {
                        energy_aux: energy_aux * self.n_units as f64,
                        battery_losses: battery_losses * self.n_units as f64,
                        temps_after_losses: zone_temp_c_after_losses,
                        total_charge: self.energy_charged.load(Ordering::SeqCst)
                            * self.n_units as f64,
                        end_of_timestep_charge: end_of_ts_charge * self.n_units as f64,
                        hb_after_only_charge_zone_temp: zone_temp_c_after_charging,
                    },
                });
        }

        self.total_time_running_current_timestep
            .store(Default::default(), Ordering::SeqCst);
        *self.service_results.write() = Default::default();
        self.energy_charged
            .store(Default::default(), Ordering::SeqCst);

        Ok(())
    }

    /// Output detailed results of heat battery calculation
    pub(crate) fn output_detailed_results(
        &self,
        hot_water_energy_output: &[f64],
    ) -> anyhow::Result<(ResultsPerTimestep, ResultsAnnual)> {
        // Define parameters to output
        // Last element of each tuple controls whether item is summed for annual total

        let output_parameters: [(&'static str, Option<&'static str>, bool); 16] = [
            ("service_name", None, false),
            ("service_type", None, false),
            ("service_on", None, false),
            ("energy_output_required", "kWh".into(), true),
            ("temp_output", "degC".into(), false),
            ("temp_inlet", "degC".into(), false),
            ("time_running", "secs".into(), true),
            ("energy_left_in_pipe", "kWh".into(), true),
            ("temperature_left_in_pipe", "degC".into(), false),
            ("energy_delivered_HB", "kWh".into(), true),
            ("energy_delivered_backup", "kWh".into(), true),
            ("energy_delivered_total", "kWh".into(), true),
            ("energy_delivered_low_temp", "kWh".into(), true),
            ("energy_charged_during_service", "kWh".into(), true),
            ("hb_zone_temperatures", "degC".into(), false),
            ("current_hb_power", "kW".into(), false),
        ];
        let aux_parameters = [
            ("energy_aux", "kWh", true),
            ("battery_losses", "kWh", true),
            ("Temps_after_losses", "degC", false),
            ("total_charge", "kWh", true),
            ("end_of_timestep_charge", "kWh", true),
            ("hb_after_only_charge_zone_temp", "degC", false),
        ];

        let mut results_per_timestep: ResultsPerTimestep =
            [("auxiliary".into(), Default::default())].into();

        // Report auxiliary parameters (not specific to a service)
        for (parameter, param_unit, _) in aux_parameters.iter() {
            if ["Temps_after_losses", "hb_after_only_charge_zone_temp"].contains(parameter) {
                let mut labels: Option<Vec<String>> = Default::default();
                for service_results in self
                    .detailed_results
                    .as_ref()
                    .expect("Detailed results accessed on heat battery when none collected.")
                    .read()
                    .iter()
                {
                    let summary = &service_results.summary;
                    let param_values = match *parameter {
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
                        results_per_timestep["auxiliary"][&(
                            String::from(label.as_str()),
                            String::from(*param_unit).into(),
                        )]
                            .push(result.into());
                    }
                }
            } else {
                // Default behaviour for scalar parameters
                let mut param_results = vec![];
                for service_results in self.detailed_results.as_ref().unwrap().read().iter() {
                    let summary = &service_results.summary;
                    let result = match *parameter {
                        "energy_aux" => summary.energy_aux,
                        "battery_losses" => summary.battery_losses,
                        "total_charge" => summary.total_charge,
                        "end_of_timestep_charge" => summary.end_of_timestep_charge,
                        _ => unreachable!(),
                    };
                    param_results.push(result.into());
                }
                results_per_timestep["auxiliary"]
                    [&(String::from(*parameter), String::from(*param_unit).into())] = param_results;
            }
        }

        // For each service, report required output parameters
        for (service_idx, service_name) in self.energy_supply_connections.keys().enumerate() {
            let mut current_results: IndexMap<(String, Option<String>), Vec<StringOrNumber>> =
                Default::default();

            // Look up each required parameter
            for (parameter, param_unit, _) in output_parameters.iter() {
                // Look up value of required parameter in each timestep
                for service_results in self
                    .detailed_results
                    .as_ref()
                    .expect("Detailed results accessed on heat battery when none collected.")
                    .read()
                    .iter()
                {
                    let current_result = &service_results.results[service_idx];
                    if *parameter == "hb_zone_temperatures" {
                        let labels = (0..current_result.hb_zone_temperatures.len())
                            .map(|i| format!("{parameter}{i}"))
                            .collect_vec();
                        for (label, result) in labels
                            .into_iter()
                            .zip(current_result.hb_zone_temperatures.iter())
                        {
                            current_results
                                .entry((String::from(label), param_unit.map(String::from)))
                                .or_default()
                                .push(result.into());
                        }
                    } else {
                        let result = current_result.param_value_as_string(parameter);
                        current_results
                            .entry((String::from(*parameter), param_unit.map(String::from)))
                            .or_default()
                            .push(result.into());
                    }
                }
            }
            // For water heating service, record hot water energy delivered from tank
            current_results[&("energy_delivered_H4".into(), Some("kWh".into()))] = if self
                .detailed_results
                .as_ref()
                .unwrap()
                .read()
                .first()
                .unwrap()
                .results[service_idx]
                .service_type
                == HeatingServiceType::DomesticHotWaterRegular
            {
                // For DHW, need to include storage and primary circuit losses.
                // Can do this by replacing H4 numerator with total energy
                // draw-off from hot water cylinder.
                // TODO (from Python) Note that the below assumes that there is only one water
                //       heating service and therefore that all hot water energy
                //       output is assigned to that service. If the model changes in
                //       future to allow more than one hot water system, this code may
                //       need to be revised to handle that scenario.
                hot_water_energy_output.iter().map(|x| x.into()).collect()
            } else {
                // TODO (from Python) Note that the below assumes there is no buffer tank for
                //       space heating, which is not currently included in the
                //       model. If this is included in future, this code will need
                //       to be revised.
                current_results[&("energy_delivered_total".into(), Some("kWh".into()))].clone()
            };

            results_per_timestep.insert(service_name.to_owned(), current_results);
        }

        let mut results_annual: ResultsAnnual = [
            (
                "Overall".into(),
                output_parameters
                    .iter()
                    .filter_map(|(parameter, param_units, incl_in_manual)| {
                        incl_in_manual.then_some((
                            (String::from(*parameter), param_units.map(String::from)),
                            0.0,
                        ))
                    })
                    .collect(),
            ),
            ("auxiliary".into(), Default::default()),
        ]
        .into();
        results_annual["Overall"].insert(("energy_delivered_H4".into(), Some("kWh".into())), 0.0);
        // Report auxiliary parameters (not specific to a service)
        for (parameter, param_unit, incl_in_annual) in aux_parameters.iter() {
            if *incl_in_annual {
                results_annual["auxiliary"].insert(
                    (String::from(*parameter), String::from(*param_unit).into()),
                    results_per_timestep["auxiliary"]
                        [&(String::from(*parameter), String::from(*param_unit).into())]
                        .iter()
                        .cloned()
                        .map(f64::from)
                        .sum::<f64>(),
                );
            }
        }
        // For each service, report required output parameters
        for service_name in self.energy_supply_connections.keys() {
            let mut current_annual_results: IndexMap<(String, Option<String>), f64> =
                Default::default();
            for (parameter, param_unit, incl_in_annual) in output_parameters.iter() {
                if *incl_in_annual {
                    let parameter_annual_total = results_per_timestep[service_name]
                        [&(String::from(*parameter), param_unit.map(String::from))]
                        .iter()
                        .cloned()
                        .map(f64::from)
                        .sum::<f64>();
                    current_annual_results.insert(
                        (String::from(*parameter), param_unit.map(String::from)),
                        parameter_annual_total,
                    );
                    results_annual["Overall"]
                        [&(String::from(*parameter), param_unit.map(String::from))] +=
                        parameter_annual_total;
                }
            }
            current_annual_results.insert(
                ("energy_delivered_H4".into(), Some("kWh".into())),
                results_per_timestep[service_name]
                    [&("energy_delivered_H4".into(), Some("kWh".into()))]
                    .iter()
                    .cloned()
                    .map(f64::from)
                    .sum::<f64>(),
            );
            results_annual["Overall"][&("energy_delivered_H4".into(), Some("kWh".into()))] +=
                results_annual[service_name][&("energy_delivered_H4".into(), Some("kWh".into()))];

            results_annual.insert(service_name.to_owned(), current_annual_results);
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

pub(crate) type ResultsPerTimestep =
    IndexMap<String, IndexMap<(String, Option<String>), Vec<StringOrNumber>>>;
pub(crate) type ResultsAnnual = IndexMap<String, IndexMap<(String, Option<String>), f64>>;

pub(crate) fn to_corpus_results_per_timestep(
    results: ResultsPerTimestep,
) -> CorpusResultsPerTimestep {
    results
        .into_iter()
        .map(|(key, value)| {
            (
                key,
                value
                    .into_iter()
                    .map(|((key1, key2), value)| {
                        (
                            (key1, key2),
                            value.into_iter().map(ResultParamValue::from).collect(),
                        )
                    })
                    .collect(),
            )
        })
        .collect()
}

pub(crate) fn to_corpus_results_annual(results: ResultsAnnual) -> CorpusResultsAnnual {
    results
        .into_iter()
        .map(|(key, value)| {
            (
                key,
                value
                    .into_iter()
                    .map(|((key1, key2), value)| {
                        ((key1, key2), vec![ResultParamValue::from(value)])
                    })
                    .collect(),
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use crate::core::common::WaterSourceWithTemperature;
    use crate::core::controls::time_control::SetpointTimeControl;
    use crate::core::controls::time_control::{ChargeControl, Control};
    use crate::core::energy_supply::energy_supply::{
        EnergySupply, EnergySupplyBuilder, EnergySupplyConnection,
    };
    use crate::core::heating_systems::common::HeatingServiceType;
    use crate::core::heating_systems::heat_battery_pcm::HeatBatteryPcm;
    use crate::core::heating_systems::heat_battery_pcm::HeatBatteryPcmServiceSpace;
    use crate::core::heating_systems::heat_battery_pcm::HeatBatteryPcmServiceWaterRegular;
    use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
    use crate::external_conditions::{DaylightSavingsConfig, ExternalConditions};
    use crate::input::{
        ControlLogicType, ExternalSensor, FuelType, HeatBattery as HeatBatteryInput,
        HeatSourceWetDetails,
    };
    use crate::simulation_time::SimulationTimeIteration;
    use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
    use approx::assert_relative_eq;
    use itertools::Itertools;
    use parking_lot::RwLock;
    use rstest::fixture;
    use rstest::rstest;
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
                heat_storage_k_j_per_k_above_phase_transition: 47.6875,
                heat_storage_k_j_per_k_below_phase_transition: 38.15,
                heat_storage_k_j_per_k_during_phase_transition: 1539.625,
                phase_transition_temperature_upper: 59.,
                phase_transition_temperature_lower: 57.,
                max_temperature: 80.,
                temp_init: 80.,
                velocity_in_hex_tube_at_1_l_per_min_m_per_s: 0.035,
                capillary_diameter_m: 0.0065,
                a: 19.744,
                b: -105.5,
                heat_exchanger_surface_area_m2: 8.83,
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
            Some(120.),
            None,
            None,
            None,
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

        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off);

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
        let heat_battery_service: HeatBatteryPcmServiceWaterRegular =
            HeatBatteryPcmServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                WaterSourceWithTemperature::ColdWaterSource(Arc::new(cold_water_source)),
                Some(Arc::new(control_min)),
                Some(Arc::new(control_max)),
            );

        assert!(heat_battery_service
            .is_on(simulation_time_iteration)
            .unwrap());
    }

    #[rstest]
    fn test_demand_energy_when_service_control_on_for_water_regular(
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
        let energy_demand = 10.;
        let temp_flow = 55.;
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
        let heat_battery_service: HeatBatteryPcmServiceWaterRegular =
            HeatBatteryPcmServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                WaterSourceWithTemperature::ColdWaterSource(Arc::new(cold_water_source)),
                Some(Arc::new(control_min)),
                Some(Arc::new(control_max)),
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

        assert_relative_eq!(result, 9.198558500698649);
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

        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        let cold_water_source = ColdWaterSource::new(vec![1.0, 1.2], 0, 1.);
        let heat_battery_service: HeatBatteryPcmServiceWaterRegular =
            HeatBatteryPcmServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                WaterSourceWithTemperature::ColdWaterSource(Arc::new(cold_water_source)),
                Some(service_control_off.clone()),
                Some(service_control_off),
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
    fn test_energy_output_max_when_service_control_on_for_water_regular(
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
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
        let heat_battery_service: HeatBatteryPcmServiceWaterRegular =
            HeatBatteryPcmServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                WaterSourceWithTemperature::ColdWaterSource(Arc::new(cold_water_source)),
                Some(Arc::new(control_min)),
                Some(Arc::new(control_max)),
            );

        let temp_flow = 50.0;
        let temp_return = 40.0;
        let result = heat_battery_service
            .energy_output_max(temp_flow, temp_return, simulation_time_iteration)
            .unwrap();

        assert_relative_eq!(result, 72279.10023958197);
    }

    #[rstest]
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
        let heat_battery_service: HeatBatteryPcmServiceWaterRegular =
            HeatBatteryPcmServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                WaterSourceWithTemperature::ColdWaterSource(Arc::new(cold_water_source)),
                Some(Arc::new(control_min)),
                Some(Arc::new(control_max)),
            );

        let temp_flow = 50.0;
        let temp_return = 40.0;
        let result = heat_battery_service
            .energy_output_max(temp_flow, temp_return, simulation_time_iteration)
            .unwrap();

        assert_eq!(result, 28882.5139822234);
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
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off);
        let heat_battery_space =
            HeatBatteryPcmServiceSpace::new(heat_battery, SERVICE_NAME.into(), ctrl.into());

        assert_eq!(
            heat_battery_space.in_required_period(simulation_time_iteration),
            Some(true)
        );
    }

    // TODO test_demand_energy_for_space

    #[rstest]
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
            HeatBatteryPcmServiceSpace::new(heat_battery, SERVICE_NAME.into(), ctrl.into());
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

        let heat_battery_service: HeatBatteryPcmServiceSpace = HeatBatteryPcmServiceSpace::new(
            heat_battery,
            SERVICE_NAME.into(),
            service_control_on.into(),
        );

        let result = heat_battery_service
            .energy_output_max(
                temp_output,
                temp_return,
                Some(time_start),
                simulation_time_iteration,
            )
            .unwrap();

        assert_relative_eq!(result, 25080.879624795467);
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
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        let create_connection_result =
            HeatBatteryPcm::create_service_connection(heat_battery.clone(), "new service");
        assert!(create_connection_result.is_ok());
        let create_connection_result =
            HeatBatteryPcm::create_service_connection(heat_battery, "new service");
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
        assert_relative_eq!(heat_battery.read().electric_charge(), 0.0);

        // electric charge should be calculated when battery control is on
        let heat_battery =
            create_heat_battery(simulation_time_iterator.clone(), battery_control_on);
        assert_relative_eq!(heat_battery.read().electric_charge(), 20.0);
    }

    #[rstest]
    fn test_first_call(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
        simulation_time: SimulationTime,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        for (t_idx, _) in simulation_time.iter().enumerate() {
            heat_battery.read().first_call();

            assert!(!heat_battery.read().flag_first_call.load(Ordering::SeqCst));

            heat_battery.read().timestep_end(t_idx).unwrap();
        }
    }

    #[rstest]
    fn test_demand_energy(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        simulation_time: SimulationTime,
        battery_control_on: Control,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);

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

    #[rstest]
    /// Check heat battery auxilary energy consumption
    fn test_calc_auxiliary_energy(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
        let heat_battery =
            create_heat_battery(simulation_time_iterator.clone(), battery_control_on);

        heat_battery
            .read()
            .calc_auxiliary_energy(1.0, 0.5, simulation_time_iterator.current_index())
            .unwrap();

        let results_by_end_user = heat_battery
            .read()
            .energy_supply
            .read()
            .results_by_end_user();

        let end_user_name = heat_battery
            .read()
            .energy_supply_connection
            .end_user_name
            .clone();

        let results_by_end_user = results_by_end_user.get(&end_user_name).unwrap();

        assert_eq!(*results_by_end_user, vec![0.0122, 0.]);
    }

    #[rstest]
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
            0.25690083365212835
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
                heat_battery.read().energy_output_max(0., None).unwrap(),
                [108864.87597021714, 124118.95144251334][t_idx]
            );

            heat_battery.read().timestep_end(t_idx).unwrap();
        }
    }
}
