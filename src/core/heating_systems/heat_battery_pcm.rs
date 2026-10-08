use crate::compare_floats::min_of_2;
/// This module provides object(s) to model the behaviour of heat batteries.
use crate::core::common::{WaterSupply, WaterSupplyBehaviour};
use crate::core::controls::time_control::{
    per_control, ChargeControl, Control, ControlBehaviour, RangeTimeControl,
};
use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::core::heating_systems::boiler::BoilerServiceWaterRegular;
use crate::core::heating_systems::common::HeatingServiceType;
use crate::core::heating_systems::heat_battery_drycore::HeatBatteryDryCoreServiceWaterRegular;
use crate::core::heating_systems::heat_network::HeatNetworkServiceWaterStorage;
use crate::core::heating_systems::heat_pump::HeatPumpServiceWater;
use crate::core::heating_systems::primary_pipework_losses_mixin::PrimaryPipeworkLossesMixin;
use crate::core::heating_systems::storage_tank::{
    DEFAULT_AMBIENT_TEMPERATURE, THERMAL_CONSTANTS_F_STO_M,
};
use crate::core::material_properties::WATER;
use crate::core::pipework::Pipework;
use crate::core::units::{
    self, KILOJOULES_PER_KILOWATT_HOUR, MILLIMETRES_IN_METRE, SECONDS_PER_HOUR, SECONDS_PER_MINUTE,
    WATTS_PER_KILOWATT,
};
use crate::core::water_heat_demand::misc::{
    calculate_volume_weighted_average_temperature, water_demand_to_kwh, WaterEventResult,
};
use crate::corpus::{ResultParamValue, ResultsAnnual, ResultsPerTimestep, TempInternalAirFn};
use crate::hem_core::external_conditions::ExternalConditions;
use crate::hem_core::simulation_time::SimulationTimeIterator;
use crate::input::{
    HeatBattery as HeatBatteryInput, HeatSourceWetDetails, PcmBatteryChargingConfiguration,
    ScheduleUnit,
};
use crate::simulation_time::SimulationTimeIteration;
use anyhow::{anyhow, bail, Ok};
use approx::relative_eq;
use arcstr::ArcStr;
use atomic_float::AtomicF64;
use fsum::FSum;
use indexmap::IndexMap;
use itertools::Itertools;
use parking_lot::RwLock;
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

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum ChargingSourceType {
    DirectElectric,
    HeatSourceWet,
}

/// Configuration for a single charging source in the HeatSource input format.
///
/// Each source has its own RangeTimeControl for hysteresis scheduling
/// and type-specific parameters (rated power for electric, heat source service
/// reference and flow temperature limit for hydronic).
///
/// Attributes:
///     source_type: ChargingSourceType.DIRECT_ELECTRIC for direct electric element,
///         ChargingSourceType.HEAT_SOURCE_WET for hydronic charging from a wet heat source.
///     control: RangeTimeControl providing hysteresis thresholds —
///         lower setpoint triggers charging start, upper triggers stop.
///         Setpoint units determined by schedule_unit.
///     rated_charge_power: Rated charging power in kW (electric sources only).
///     heat_source_service: Reference to the heat source service object that
///         provides hot water for charging (hydronic sources only). Set by
///         project.py during post-construction linking; None until then.
///     temp_flow_max: Maximum flow temperature (°C) the heat source
///         should provide when charging (hydronic sources only). Used to
///         cap the heat source flow temperature request during battery charging
///         and to estimate return temperature for energy demand calculations.
///     flow_rate_charging_l_per_min: Flow rate through the heat exchanger
///         during hydronic charging (litre/minute). May differ from the
///         battery's main flow_rate_l_per_min if charging and discharging
///         circuits have different pipework (hydronic sources only).
///     schedule_unit: How to interpret the control schedule setpoints.
///         "soc" (default): values are state-of-charge fractions (0–1).
///         "temperature": values are temperatures (°C), converted to SOC
///         internally for comparison against the battery's current state.
/// Enum representing the union of all heat source service types that can provide hydronic charging.
/// Defined here (not in _base.py) to avoid circular imports — the concrete
/// service types are defined across multiple modules that import from _base.py.

#[derive(Debug, Clone)]
pub(crate) enum HeatSourceWetService {
    HeatPumpServiceWater(HeatPumpServiceWater),
    BoilerServiceWaterRegular(BoilerServiceWaterRegular),
    HeatBatteryPCMServiceWaterRegular(HeatBatteryPcmServiceWaterRegular),
    HeatBatteryDryCoreServiceWaterRegular(HeatBatteryDryCoreServiceWaterRegular),
    HeatNetworkServiceWaterStorage(HeatNetworkServiceWaterStorage),
}
// Set to be ergonomic for current use cases by providing a unified interface for all heat source services.
// may need to be more flexible in the future
impl HeatSourceWetService {
    fn energy_output_max(
        &self,
        temp_flow: f64,
        temp_return: f64,
        simtime: &SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        match self {
            Self::HeatPumpServiceWater(service) => Ok(service
                .energy_output_max(temp_flow, temp_return, *simtime)?
                .0),
            Self::BoilerServiceWaterRegular(service) => {
                Ok(service.energy_output_max(temp_flow, temp_return, None, *simtime))
            }
            Self::HeatBatteryPCMServiceWaterRegular(service) => {
                service.energy_output_max(temp_flow, temp_return, *simtime, false)
            }
            Self::HeatBatteryDryCoreServiceWaterRegular(service) => {
                service.energy_output_max(temp_flow, temp_return, *simtime)
            }
            Self::HeatNetworkServiceWaterStorage(service) => {
                Ok(service.energy_output_max(temp_flow, temp_return, simtime))
            }
        }
    }

    fn demand_energy(
        &self,
        energy_demand: f64,
        temp_flow: f64,
        temp_return: f64,
        simtime: &SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        match self {
            Self::HeatPumpServiceWater(service) => {
                service.demand_energy(energy_demand, Some(temp_flow), Some(temp_return), *simtime)
            }
            Self::BoilerServiceWaterRegular(service) => Ok(service
                .demand_energy(
                    energy_demand,
                    temp_flow,
                    Some(temp_return),
                    None,
                    None,
                    Some(true),
                    *simtime,
                )?
                .0),
            Self::HeatBatteryPCMServiceWaterRegular(service) => service.demand_energy(
                energy_demand,
                Some(temp_flow),
                Some(temp_return),
                Some(true),
                *simtime,
                false,
            ),
            Self::HeatBatteryDryCoreServiceWaterRegular(service) => service.demand_energy(
                energy_demand,
                Some(temp_flow),
                temp_return,
                Some(true),
                *simtime,
            ),
            Self::HeatNetworkServiceWaterStorage(service) => {
                service.demand_energy(energy_demand, temp_flow, Some(temp_return), simtime)
            }
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct HeatBatteryChargingSource {
    source_type: ChargingSourceType,
    control: Arc<RangeTimeControl>,
    rated_charge_power: Option<f64>,
    heat_source_service: Option<HeatSourceWetService>,
    temp_flow_max: f64,
    flow_rate_charging_l_per_min: Option<f64>,
    hex_a: Option<f64>,
    hex_b: Option<f64>,
    hex_velocity_at_1_l_per_min: Option<f64>,
    hex_capillary_diameter_m: Option<f64>,
    schedule_unit: ScheduleUnit,
}

///    Check that no two charging sources have overlapping active periods.
///
///    For each timestep in the simulation, at most one source may have a
///    non-null schedule_upper (i.e. be in an active or transition period).
///    Raises ValueError if any overlap is found.
///
///    Iterates the shared SimulationTime through all timesteps so that
///    each source's RangeTimeControl.setpnt() returns the correct
///    per-timestep value, then resets to initial state. Must be called
///    during construction before the simulation loop begins.
///
///    Args:
///        heat_source_data: Dict of charging sources with their controls.
///        battery_name: Name of the heat battery, for error messages.
///        simtime: Shared SimulationTime iterator used to advance through
///            timesteps for per-step schedule evaluation.
fn validate_no_schedule_overlap(
    heat_source_data: IndexMap<ArcStr, HeatBatteryChargingSource>,
    battery_name: &str,
    simtime_iterator: &SimulationTimeIterator,
) -> anyhow::Result<()> {
    if simtime_iterator.current_index() != 0 {
        bail!(
            "HeatBattery '{}': validate_no_schedule_overlap must be called before the simulation starts (current timestep index: {}).",
            battery_name,
            simtime_iterator.current_index()
        );
    }

    for (t_idx, t_it) in simtime_iterator.clone().enumerate() {
        let mut active_sources: Vec<ArcStr> = Vec::new();
        for (src_name, charging_source) in &heat_source_data {
            if let (_, Some(_)) = charging_source.control.setpnt_range_time_control(&t_it) {
                active_sources.push(src_name.clone());
            }
        }
        if active_sources.len() > 1 {
            bail!(
                "HeatBattery '{}': charging sources {:?} have overlapping active schedules at timestep index {} (first overlapping timestep). Each source must have non-overlapping RangeTimeControl schedules.",
                battery_name,
                active_sources,
                t_idx
            );
        }
    }
    Ok(())
}

/// An object to represent a water heating service provided by a regular heat battery.
///
/// This object contains the parts of the heat battery calculation that are
/// specific to providing hot water.
#[derive(Debug, Clone)]
pub(crate) struct HeatBatteryPcmServiceWaterRegular {
    heat_battery: Arc<RwLock<HeatBatteryPcm>>,
    service_name: ArcStr,
    cold_feed: WaterSupply,
    control: Arc<RangeTimeControl>,
}

impl HeatBatteryPcmServiceWaterRegular {
    /// Arguments:
    /// * `heat_battery` - reference to the Heat Battery object providing the service
    /// * `service_name` - name of the service demanding energy
    /// * `cold_feed` - reference to ColdWaterSource object
    /// * `control` - reference to a RangeTimeControl object, combining controlmax and controlmin.
    ///    Takes precedence if set in python, one is constructed from 2 control objects otherwise.
    ///
    /// * From python, used to create a range time control before this step
    /// * `control_min` - reference to a control object which must select current the minimum timestep temperature
    /// * `control_max` - reference to a control object which must select current the maximum timestep temperature
    pub(crate) fn new(
        heat_battery: Arc<RwLock<HeatBatteryPcm>>,
        service_name: ArcStr,
        cold_feed: WaterSupply,
        control: Arc<RangeTimeControl>,
    ) -> Self {
        Self {
            heat_battery,
            service_name,
            cold_feed,
            control,
        }
    }

    /// Return setpoint (not necessarily temperature)
    pub(crate) fn setpnt(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (Option<f64>, Option<f64>) {
        self.control
            .setpnt_range_time_control(&simulation_time_iteration)
    }

    /// Demand energy (in kWh) from the heat_battery
    pub(crate) fn demand_energy(
        &self,
        energy_demand: f64,
        temp_flow: Option<f64>,
        temp_return: Option<f64>,
        update_heat_source_state: Option<bool>,
        simtime: SimulationTimeIteration,
        ignore_standard_ctrl: bool,
    ) -> anyhow::Result<f64> {
        let service_on = self.is_on(simtime) || ignore_standard_ctrl;
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
            simtime,
        )
    }

    pub(crate) fn energy_output_max(
        &self,
        temp_flow: f64,
        temp_return: f64,
        simtime: SimulationTimeIteration,
        ignore_standard_ctrl: bool,
    ) -> anyhow::Result<f64> {
        let service_on = self.is_on(simtime) || ignore_standard_ctrl;
        if !service_on {
            return Ok(0.);
        }

        self.heat_battery
            .read()
            .energy_output_max(temp_flow, temp_return, Some(0.), simtime)
    }

    fn is_on(&self, simtime: SimulationTimeIteration) -> bool {
        self.control.is_on(&simtime)
    }
}

/// An object to represent a direct water heating service provided by a heat battery.
///
/// This is similar to a combi boiler or HIU providing hot water on demand.
#[derive(Debug)]
pub struct HeatBatteryPcmServiceWaterDirect {
    heat_battery: Arc<RwLock<HeatBatteryPcm>>,
    service_name: ArcStr,
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
        service_name: ArcStr,
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

    fn temp_hot_water(&self, vol: f64, simtime: SimulationTimeIteration) -> anyhow::Result<f64> {
        let list_temp_vol = self.cold_feed.get_temp_cold_water(vol, simtime)?;
        let inlet_temp =
            calculate_volume_weighted_average_temperature(list_temp_vol, Some(vol), None)?;

        self.heat_battery
            .read()
            .get_temp_hot_water(inlet_temp, vol, self.setpoint_temp, simtime)
    }

    pub(crate) fn get_temp_hot_water(
        &self,
        volume_req: f64,
        volume_req_already: Option<f64>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        let volume_req_already = volume_req_already.unwrap_or(0.);

        if relative_eq!(volume_req, 0., max_relative = 1e-09, epsilon = 1e-10) {
            bail!("volume_req must be non-zero");
        }

        let volume_req_cumulative = volume_req + volume_req_already;
        let temp_hot_water_cumulative =
            self.temp_hot_water(volume_req_cumulative, simulation_time_iteration)?;

        // Base temperature on the part of the draw-off for volume_req, and
        // ignore any volume previously considered
        let temp_hot_water_req = if relative_eq!(
            volume_req_already,
            0.,
            max_relative = 1e-09,
            epsilon = 1e-10
        ) {
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
                if relative_eq!(event.volume_hot, 0., max_relative = 1e-09, epsilon = 1e-10) {
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
                } else {
                    bail!("No hot water temperatures were available for an event, unexpectedly");
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
            Some(cold_water_temp), // return temperature (cold water inlet)
            None,                  // flow temperature (hot water outlet)
            true,
            None,
            Some(true),
            simtime,
        )
    }
}

#[derive(Clone, Debug)]
pub struct HeatBatteryPcmServiceSpace {
    heat_battery: Arc<RwLock<HeatBatteryPcm>>,
    service_name: ArcStr,
    control: Control,
}

/// An object to represent a space heating service provided by a heat_battery to e.g. radiators.
///
/// This object contains the parts of the heat battery calculation that are
/// specific to providing space heating.
impl HeatBatteryPcmServiceSpace {
    pub(crate) fn new(
        heat_battery: Arc<RwLock<HeatBatteryPcm>>,
        service_name: ArcStr,
        control: Control, // in Python this is SetpointTimeControl | CombinationTimeControl
    ) -> Self {
        Self {
            heat_battery,
            service_name,
            control,
        }
    }

    pub fn temp_setpnt(&self, simulation_time_iteration: SimulationTimeIteration) -> Option<f64> {
        per_control!(&self.control, ctrl => { ctrl.setpnt(&simulation_time_iteration) })
    }

    pub fn in_required_period(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> Option<bool> {
        per_control!(&self.control, ctrl => { ctrl.in_required_period(&simulation_time_iteration) })
    }

    /// Demand energy (in kWh) from the heat battery
    pub fn demand_energy(
        &self,
        energy_demand: f64,
        temp_flow: f64,
        temp_return: f64,
        time_start: Option<f64>,
        update_heat_source_state: Option<bool>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);
        let service_on = self.is_on(simtime);

        let energy_demand = if !service_on { 0.0 } else { energy_demand };

        self.heat_battery.read().demand_energy(
            &self.service_name,
            HeatingServiceType::Space,
            energy_demand,
            Some(temp_return),
            Some(temp_flow),
            service_on,
            time_start,
            Some(update_heat_source_state),
            simtime,
        )
    }

    fn is_on(&self, simulation_time_iteration: SimulationTimeIteration) -> bool {
        per_control!(&self.control, ctrl => { ctrl.is_on(&simulation_time_iteration) })
    }

    pub(crate) fn energy_output_max(
        &self,
        temp_output: f64,
        temp_return_feed: f64,
        time_start: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let time_start = time_start.unwrap_or(0.);

        if !self.is_on(simtime) {
            return Ok(0.);
        }

        self.heat_battery.read().energy_output_max(
            temp_output,
            temp_return_feed,
            Some(time_start),
            simtime,
        )
    }
}

const DEFAULT_N_LAYERS: usize = 8; // Number of calculation layers in heat battery
const DEFAULT_TIME_STEP_SECONDS: f64 = 20.; // Time step for iterative calculations (seconds)
const DEFAULT_INLET_TEMP_CELSIUS: f64 = 10.; // Initial inlet temperature for Reynolds number calculation (°C)
const DEFAULT_OUTLET_TEMP_CELSIUS: f64 = 53.; // Estimated outlet temperature for Reynolds number calculation (°C)

// Surrounding air temperature assumed during the standing-loss
// characterisation test. max_rated_losses is the loss measured with the
// battery fully charged (at max_temperature) against this ambient, so the
// rated temperature difference is (max_temperature - this value). The value
// mirrors the reference ambient used for hot water cylinder standby losses
// (BS EN 12897:2016).
// nothing seems to read this - check upstream whether service_results field is necessary
const TEMP_AMBIENT_RATED_LOSSES_C: f64 = 20.0;

// Near-equilibrium per-sub-timestep energy transfers are small differences of
// larger flows, so their floating-point cancellation floor sits well above the
// 1e-10 used for one-off energy comparisons elsewhere. Testing such a value
// against zero with too tight a tolerance lets it take a different branch on
// different platforms, which accumulates into a divergent result (seen as a
// Windows-vs-Linux e2e difference). 1e-6 is the tightest tolerance that holds
// the decision stable across platforms, and is still ~7 orders of magnitude
// below any meaningful transfer (a delivering zone moves tens of kJ per step).
// Temperature comparisons do not suffer this cancellation, so they keep 1e-10.
const NEGLIGIBLE_ENERGY_KJ: f64 = 1e-6;
const NEGLIGIBLE_TEMP_DIFF_C: f64 = 1e-10;

// Charging flow temperature is set above the target uniform PCM temperature
// by this heat-exchanger approach difference. A heat source must run its
// flow above the store's target temperature to drive heat across the
// exchanger; charging at exactly the target would collapse the driving
// temperature difference to zero as the layers approach it, so the store
// would only ever approach the target asymptotically and never reach it
// within a timestep. 5 °C is a typical charging approach difference. The
// flow temperature is capped at the source's maximum, and the energy demand
// is capped at the target SOC, so the store charges towards its target
// temperature and does not exceed it.
const CHARGE_APPROACH_TEMP_DIFF_C: f64 = 5.0;
#[derive(Clone, Debug)]
#[allow(dead_code)]
struct HeatBatteryResult {
    service_name: ArcStr,
    service_type: Option<HeatingServiceType>,
    service_on: bool,
    energy_output_required: f64,
    temp_output: Option<f64>,
    temp_inlet: Option<f64>,
    time_running: f64,
    energy_delivered_hb: f64,
    energy_delivered_backup: f64,
    energy_delivered_total: f64,
    energy_charged_during_service: f64,
    hb_zone_temperatures: Vec<f64>,
    current_hb_power: ResultParamValue,
}

impl HeatBatteryResult {
    fn param(&self, param: &str) -> ResultParamValue {
        match param {
            "service_name" => ResultParamValue::String(self.service_name.clone()),
            "service_type" => self
                .service_type
                .as_ref()
                .map(|service_type| ResultParamValue::String((*service_type).into()))
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
            "current_hb_power" => self.current_hb_power.clone(),
            _ => panic!("Unknown parameter: {}", param),
        }
    }
}

// Define parameters to output. The third tuple element controls whether each
// parameter is summed for the annual total.
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
const AUX_PARAMETERS: [(&str, Option<&str>, bool); 7] = [
    ("energy_aux", Some("kWh"), true),
    ("battery_losses", Some("kWh"), true),
    ("Temps_after_losses", Some("degC"), false),
    ("total_charge", Some("kWh"), true),
    ("end_of_timestep_charge", Some("kWh"), true),
    ("hb_after_only_charge_zone_temp", Some("degC"), false),
    ("state_of_charge", Some("ratio"), false),
];

#[derive(Debug)]
struct HeatBatteryTimestepSummary {
    energy_aux: f64,
    battery_losses: f64,
    temps_after_losses: Vec<f64>,
    total_charge: f64,
    end_of_timestep_charge: f64,
    hb_after_only_charge_zone_temp: Vec<f64>,
    state_of_charge: f64,
}

impl HeatBatteryTimestepSummary {
    fn param(&self, param: &str) -> ResultParamValue {
        match param {
            "energy_aux" => self.energy_aux.into(),
            "battery_losses" => self.battery_losses.into(),
            "total_charge" => self.total_charge.into(),
            "end_of_timestep_charge" => self.end_of_timestep_charge.into(),
            "state_of_charge" => self.state_of_charge.into(),
            _ => panic!("Parameter {param} not recognised"),
        }
    }
}

#[derive(Debug)]
struct HeatBatteryTimestepResult {
    results: Vec<HeatBatteryResult>,
    summary: HeatBatteryTimestepSummary,
}

fn timestep_record_for_service<'a>(
    timestep_service_results: &'a HeatBatteryTimestepResult,
    service_name: &ArcStr,
) -> &'a HeatBatteryResult {
    timestep_service_results
        .results
        .iter()
        .find(|record| &record.service_name == service_name)
        .expect("each timestep must contain a record for every registered service")
}

#[derive(Clone, Debug)]
struct PipeEnergy {
    energy: f64,
    temperature: f64,
}

#[derive(Debug)]
pub struct HeatBatteryPcm {
    hb_time_step: f64,
    n_layers: usize,
    initial_inlet_temp: f64,
    estimated_outlet_temp: f64,
    energy_supply: Arc<RwLock<EnergySupply>>,
    energy_supply_connection: EnergySupplyConnection,
    simulation_time_step: f64,
    external_conditions: Arc<ExternalConditions>,
    energy_supply_connections: IndexMap<ArcStr, EnergySupplyConnection>,
    use_heatsource_data: bool,
    heat_source_data: Option<IndexMap<ArcStr, HeatBatteryChargingSource>>,
    charging_active: Arc<RwLock<IndexMap<ArcStr, bool>>>,
    charge_control: Option<Arc<ChargeControl>>,
    pwr_in: f64,
    max_rated_losses: f64,
    power_circ_pump: f64,
    power_standby: f64,
    n_units: usize,
    time_unit: u32,
    total_time_running_current_timestep: AtomicF64,
    pump_running_time_current_timestep: AtomicF64,
    time_running_direct_current_timestep: AtomicF64,
    flag_first_call: AtomicBool,
    charge_level: f64,
    energy_charged_total: AtomicF64,
    energy_charged_electric: AtomicF64,
    battery_losses: AtomicF64,
    pipework_primary_gains_kwh: AtomicF64,
    simultaneous_charging_and_discharging: bool,
    max_temp_of_charge: f64,
    energy_stored_max: f64,
    temp_diff_rated_losses: f64,
    zone_temp_c_dist_initial: Arc<RwLock<Vec<f64>>>,
    heat_storage_kj_per_k_below: f64,
    heat_storage_kj_per_k_during: f64,
    heat_storage_kj_per_k_above: f64,
    phase_transition_temperature_upper: f64,
    phase_transition_temperature_lower: f64,
    velocity_in_hex_tube: f64,
    capillary_diameter_m: f64,
    a: f64,
    b: f64,
    flow_rate_l_per_min: f64,
    service_results: Arc<RwLock<Vec<HeatBatteryResult>>>,
    energy_delivered_by_service: Arc<RwLock<IndexMap<ArcStr, f64>>>,
    detailed_results: Option<Arc<RwLock<Vec<HeatBatteryTimestepResult>>>>,
    flag_1_warning: [bool; 2],
    temp_ref: f64,
    pipework: PrimaryPipeworkLossesMixin,
}
/// PCM heat battery that can be charged electrically, hydronically, or both.
/// Models a phase-change-material heat battery that can be charged by electric
/// elements, hydronic heat sources (heat pumps, solar thermal, etc.), or a
/// combination of both with non-overlapping schedules. Provides space heating
/// and hot water (regular and direct).
impl HeatBatteryPcm {
    pub(crate) fn new(
        heat_battery_details: &HeatSourceWetDetails,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_connection: EnergySupplyConnection,
        simulation_time: SimulationTimeIterator,
        external_conditions: Arc<ExternalConditions>,
        temp_min_useful: Option<f64>,
        temp_internal_air_callback: TempInternalAirFn,
        charge_control: Option<Arc<ChargeControl>>,
        heat_source_data: Option<IndexMap<ArcStr, HeatBatteryChargingSource>>,
        n_layers: Option<usize>,
        hb_time_step: Option<f64>,
        initial_inlet_temp: Option<f64>,
        estimated_outlet_temp: Option<f64>,
        output_detailed_results: Option<bool>,
        primary_pipework: Option<Vec<Pipework>>,
    ) -> anyhow::Result<Self> {
        // 20secs is the current preferred timestep to run the heat battery iterative calculations
        let hb_time_step = hb_time_step.unwrap_or(DEFAULT_TIME_STEP_SECONDS);
        let n_layers = n_layers.unwrap_or(DEFAULT_N_LAYERS);
        let initial_inlet_temp = initial_inlet_temp.unwrap_or(DEFAULT_INLET_TEMP_CELSIUS);
        let estimated_outlet_temp = estimated_outlet_temp.unwrap_or(DEFAULT_OUTLET_TEMP_CELSIUS);
        let output_detailed_results = output_detailed_results.unwrap_or(false);

        // Determine charging configuration mode: HeatSource dict with
        // per-source RangeTimeControl vs single ChargeControl + rated_charge_power
        let use_heatsource_data = heat_source_data.is_some();
        let pipework = PrimaryPipeworkLossesMixin::new(
            primary_pipework.unwrap_or_default(),
            Arc::new(|_| 0.0), // TODO: THis needs thinking about
            temp_internal_air_callback,
            &simulation_time.current_iteration(),
        );
        // Per-source hysteresis state: tracks whether each source is currently
        // in its active charging band (SOC below upper setpoint after being
        // triggered by SOC falling below lower setpoint)
        let charging_active = if let Some(heat_source_data) = &heat_source_data {
            let charging_active: IndexMap<ArcStr, bool> = heat_source_data
                .iter()
                .to_owned()
                .map(|(k, _)| (k.clone(), false))
                .collect();
            Arc::new(RwLock::new(charging_active))
        } else {
            Arc::new(RwLock::new(IndexMap::new()))
        };

        // Warn when temperature-based schedules exceed the heat source's flow
        // temperature limit — the battery cannot reach those targets, so
        // charging will cycle inefficiently and may degrade performance.
        for (_, src) in &heat_source_data.clone().unwrap_or_default() {
            if src.schedule_unit == ScheduleUnit::Temperature {
                HeatBatteryPcm::warn_if_schedule_exceeds_flow_temp(src, &simulation_time)?;
            }
        }

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
                    charging_config,
                    max_rated_losses,
                    electricity_circ_pump: power_circ_pump,
                    electricity_standby: power_standby,
                    number_of_units: n_units,
                    simultaneous_charging_and_discharging,
                    max_temperature,
                    heat_storage_kj_per_k_above_phase_transition: heat_storage_kj_per_k_above,
                    heat_storage_kj_per_k_below_phase_transition: heat_storage_kj_per_k_below,
                    heat_storage_kj_per_k_during_phase_transition: heat_storage_kj_per_k_during,
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
            let pwr_in = match charging_config {
                PcmBatteryChargingConfiguration::ChargeControl {
                    rated_charge_power, ..
                } => rated_charge_power,
                PcmBatteryChargingConfiguration::RangeControl { .. } => {
                    todo!("as part of migration to alpha9")
                }
            };

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

        let detailed_results: Option<Arc<RwLock<Vec<HeatBatteryTimestepResult>>>> =
            if output_detailed_results {
                Some(Default::default())
            } else {
                None
            };

        let temp_diff_rated_losses = max_temp_of_charge - TEMP_AMBIENT_RATED_LOSSES_C;
        if temp_diff_rated_losses <= 0.0 {
            bail!(
                "Heat battery max_temperature must exceed the rated-loss reference ambient temperature ({TEMP_AMBIENT_RATED_LOSSES_C} °C)."
            );
        }
        // Cumulative energy delivered to each service over the whole calculation, in kWh
        // (across all units). Unlike __service_results this is never reset at timestep end; it
        // is read after the simulation to apportion the battery's charging energy between the
        // services it feeds in proportion to the output each received.
        let energy_delivered_by_service: Arc<RwLock<IndexMap<ArcStr, f64>>> =
            Arc::new(RwLock::new(IndexMap::new()));

        // Minimum useful temperature for SOC calculation — the temperature at which
        // the battery is considered fully discharged (SOC=0).
        let temp_ref = temp_min_useful.unwrap_or(0.);

        if max_temp_of_charge < temp_ref
            || relative_eq!(max_temp_of_charge, temp_ref, epsilon = 1e-10)
        {
            bail!(
                "HeatBatteryPCM: max_temperature ({max_temp_of_charge} °C) must be 
                greater than the SOC reference temperature ({temp_ref} °C). 
                A battery with zero storable energy cannot function."
            );
        }
        let energy_stored_max = HeatBatteryPcm::calculate_layer_energy_stored(
            max_temp_of_charge,
            temp_ref,
            phase_transition_temperature_lower,
            phase_transition_temperature_upper,
            heat_storage_kj_per_k_below / n_layers as f64,
            heat_storage_kj_per_k_during / n_layers as f64,
            heat_storage_kj_per_k_above / n_layers as f64,
        ) * n_layers as f64;
        Ok(Self {
            hb_time_step,
            n_layers,
            initial_inlet_temp,
            estimated_outlet_temp,
            energy_supply,
            energy_supply_connection,
            simulation_time_step: simulation_time.step_in_hours(),
            external_conditions,
            energy_supply_connections: Default::default(),
            use_heatsource_data,
            heat_source_data,
            charging_active,
            charge_control,
            pwr_in,
            max_rated_losses,
            power_circ_pump,
            power_standby,
            n_units,
            time_unit: units::SECONDS_PER_HOUR,
            total_time_running_current_timestep: Default::default(),
            pump_running_time_current_timestep: Default::default(),
            time_running_direct_current_timestep: Default::default(),
            flag_first_call: true.into(),
            charge_level: Default::default(),
            energy_charged_total: Default::default(),
            energy_charged_electric: Default::default(),
            battery_losses: Default::default(),
            pipework_primary_gains_kwh: Default::default(),
            simultaneous_charging_and_discharging,
            max_temp_of_charge,
            energy_stored_max,
            temp_diff_rated_losses,
            zone_temp_c_dist_initial: Arc::new(RwLock::new(vec![temp_init; n_layers])),
            heat_storage_kj_per_k_below: heat_storage_kj_per_k_below / n_layers as f64,
            heat_storage_kj_per_k_during: heat_storage_kj_per_k_during / n_layers as f64,
            heat_storage_kj_per_k_above: heat_storage_kj_per_k_above / n_layers as f64,
            phase_transition_temperature_upper,
            phase_transition_temperature_lower,
            velocity_in_hex_tube,
            capillary_diameter_m: inlet_diameter_mm / MILLIMETRES_IN_METRE as f64,
            a,
            b,
            flow_rate_l_per_min,
            service_results: Default::default(),
            energy_delivered_by_service,
            detailed_results,
            flag_1_warning: [true; 2],
            temp_ref: temp_min_useful.unwrap_or(0.),
            pipework,
        })
    }

    fn warn_if_schedule_exceeds_flow_temp(
        source: &HeatBatteryChargingSource,
        simtime_iterator: &SimulationTimeIterator,
    ) -> anyhow::Result<()> {
        if simtime_iterator.current_index() != 0 {
            bail!(
                "warn_if_schedule_exceeds_flow_temp must not be called after the simulation has started, as it resets simulation_time."
            )
        }
        let mut temp_observed_max = source.temp_flow_max;
        for _ in simtime_iterator.clone().enumerate() {
            let (lower, upper) = source
                .control
                .setpnt_range_time_control(&simtime_iterator.current_iteration());
            if let Some(lower) = lower {
                temp_observed_max = lower.max(temp_observed_max);
            }
            if let Some(upper) = upper {
                temp_observed_max = upper.max(temp_observed_max);
            }
        }
        Ok(())
    }
    /// Calculate energy stored in a single layer above a reference temperature.
    ///
    /// Accounts for the three thermal regimes of PCM:
    /// - Below phase transition: sensible heat with below-transition heat capacity
    /// - During phase transition: latent + sensible with during-transition heat capacity
    /// - Above phase transition: sensible heat with above-transition heat capacity
    ///
    /// The six mutually exclusive cases are determined by where temp_ref and
    /// temp_layer sit relative to the phase transition band (temp_lower, temp_upper).
    ///
    /// Args:
    ///     temp_layer: Current temperature of the layer (°C).
    ///     temp_ref: Reference temperature (°C) — energy below this is not counte
    ///  Returns:
    ///      Energy stored in the layer above temp_ref (kJ). Returns 0 if
    ///     temp_layer <= temp_ref.
    fn calculate_layer_energy_stored(
        temp_layer: f64,
        temp_ref: f64,
        temp_lower: f64,
        temp_upper: f64,
        cap_below: f64,
        cap_during: f64,
        cap_above: f64,
    ) -> f64 {
        if temp_layer <= temp_ref {
            return 0.0;
        }

        if temp_layer <= temp_lower {
            // Both below phase transition
            (temp_layer - temp_ref) * cap_below
        } else if temp_ref >= temp_upper {
            // Both above phase transition
            (temp_layer - temp_ref) * cap_above
        } else if temp_ref >= temp_lower {
            // temp_ref in transition band
            if temp_layer <= temp_upper {
                // Both in transition band
                (temp_layer - temp_ref) * cap_during
            } else {
                // temp_ref in transition, temp_layer above
                (temp_upper - temp_ref) * cap_during + (temp_layer - temp_upper) * cap_above
            }
        } else {
            // temp_ref below transition band
            if temp_layer <= temp_upper {
                // temp_ref below, temp_layer in transition
                (temp_lower - temp_ref) * cap_below + (temp_layer - temp_lower) * cap_during
            } else {
                // Spans all three regions
                (temp_lower - temp_ref) * cap_below
                    + (temp_upper - temp_lower) * cap_during
                    + (temp_layer - temp_upper) * cap_above
            }
        }
    }

    /// Calculate energy-based state of charge for the PCM heat battery.
    ///
    /// For each layer, calculates stored energy above temp_ref across the
    /// three PCM thermal regimes (below, during, and above phase transition).
    /// SOC is the ratio of total stored energy to maximum storable energy (all
    /// layers at temp_charge_max).
    ///
    /// Args:
    ///     zone_temps: List of current temperatures for each battery layer (°C).
    ///
    /// Returns:
    ///     State of charge as a float between 0.0 (fully discharged to
    ///     temp_ref) and 1.0 (all layers at temp_charge_max).
    fn calc_state_of_charge(&self, zone_temps: Vec<f64>) -> anyhow::Result<f64> {
        let temp_ref = self.temp_ref;

        let energy_stored_total: f64 = zone_temps
            .iter()
            .map(|&t| {
                HeatBatteryPcm::calculate_layer_energy_stored(
                    t,
                    temp_ref,
                    self.phase_transition_temperature_lower,
                    self.phase_transition_temperature_upper,
                    self.heat_storage_kj_per_k_below,
                    self.heat_storage_kj_per_k_during,
                    self.heat_storage_kj_per_k_above,
                )
            })
            .sum();

        let mut soc = energy_stored_total / self.energy_stored_max;

        if soc > 1.0 + 1e-10 {
            bail!(
                "State of charge exceeds 1.0 {} Zone temperatures exceed temp_charge_max ({} °C).",
                soc,
                self.max_temp_of_charge
            );
        }

        soc = soc.min(1.0);

        Ok(soc)
    }

    /// Convert a temperature setpoint to equivalent state of charge.
    ///
    /// Assumes all layers are at the given temperature — the simplifying
    /// assumption appropriate for a schedule setpoint (uniform target).
    ///
    /// Args:
    ///     temp: Target temperature in °C.
    ///
    /// Returns:
    ///      SOC value (0–1) corresponding to the given temperature.
    ///
    fn temp_to_soc(&self, temp: f64) -> anyhow::Result<f64> {
        let temp_ref = self.temp_ref;
        let energy_stored = HeatBatteryPcm::calculate_layer_energy_stored(
            temp,
            temp_ref,
            self.phase_transition_temperature_lower,
            self.phase_transition_temperature_upper,
            self.heat_storage_kj_per_k_below,
            self.heat_storage_kj_per_k_during,
            self.heat_storage_kj_per_k_above,
        ) * self.n_layers as f64;

        let mut soc = energy_stored / self.energy_stored_max;

        if soc > 1.0 + 1e-10 {
            bail!(
                "State of charge exceeds 1.0 ({soc}). Target temperature ({temp} °C) exceeds temp_charge_max ({} °C).",
                self.max_temp_of_charge
            );
        }

        soc = soc.min(1.0);

        Ok(soc)
    }

    /// Convert a state of charge to the equivalent uniform temperature.
    ///
    ///  Inverts the piecewise-linear energy function used by __temp_to_soc.
    ///  Assumes all layers are at the same temperature (uniform target), which
    ///  is appropriate for converting a schedule SOC target to a temperature
    ///  setpoint
    ///
    ///  The energy function has three thermal regimes separated by the phase
    ///  transition band. This method computes the target energy from the SOC,
    ///  determines which regime it falls in, and inverts the corresponding
    ///  linear segment to recover the temperature
    ///
    ///  Args:
    ///      soc: State of charge (0.0–1.0)
    ///  Returns:
    ///      Temperature in °C corresponding to the given SOC
    fn soc_to_temp(&self, soc: f64) -> f64 {
        if soc <= 0.0 {
            return self.temp_ref;
        }
        if soc >= 1.0 {
            return self.max_temp_of_charge;
        }
        let temp_ref = self.temp_ref;
        let temp_lower = self.phase_transition_temperature_lower;
        let temp_upper = self.phase_transition_temperature_upper;
        let cap_below = self.heat_storage_kj_per_k_below;
        let cap_during = self.heat_storage_kj_per_k_during;
        let cap_above = self.heat_storage_kj_per_k_above;

        let energy_max = HeatBatteryPcm::calculate_layer_energy_stored(
            self.max_temp_of_charge,
            temp_ref,
            temp_lower,
            temp_upper,
            cap_below,
            cap_during,
            cap_above,
        );
        let energy_target = soc * energy_max;

        // Energy at the boundaries of the phase transition band, accumulated
        // from temp_ref upward. These define the regime boundaries in energy
        // space.
        let energy_at_lower = HeatBatteryPcm::calculate_layer_energy_stored(
            temp_lower, temp_ref, temp_lower, temp_upper, cap_below, cap_during, cap_above,
        );
        let energy_at_upper = HeatBatteryPcm::calculate_layer_energy_stored(
            temp_upper, temp_ref, temp_lower, temp_upper, cap_below, cap_during, cap_above,
        );
        if energy_target <= energy_at_lower {
            // Target falls in the below-transition regime
            temp_ref + (energy_target / cap_below)
        } else if energy_target <= energy_at_upper {
            // Target falls in the phase-transition regime
            temp_lower + ((energy_target - energy_at_lower) / cap_during)
        } else {
            // Target falls in the above-transition regime
            temp_upper + ((energy_target - energy_at_upper) / cap_above)
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
    /// * `control` - reference to a RangeTimeControl
    pub(crate) fn create_service_hot_water_regular(
        heat_battery: Arc<RwLock<Self>>,
        service_name: &str,
        cold_feed: WaterSupply,
        control: Arc<RangeTimeControl>,
    ) -> anyhow::Result<HeatBatteryPcmServiceWaterRegular> {
        Self::create_service_connection(heat_battery.clone(), service_name)?;
        Ok(HeatBatteryPcmServiceWaterRegular::new(
            heat_battery,
            service_name.into(),
            cold_feed,
            control,
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
        control: Control, // in Python this is SetpointTimeControl | CombinationTimeControl
    ) -> anyhow::Result<HeatBatteryPcmServiceSpace> {
        Self::create_service_connection(heat_battery.clone(), service_name)?;
        Ok(HeatBatteryPcmServiceSpace::new(
            heat_battery,
            service_name.into(),
            control,
        ))
    }

    /// Return recoverable standing losses and accumulated pipework gains for the timestep.
    ///
    ///        Includes both standing heat losses from the battery casing and any
    ///        internal pipework gains from hydronic charging. Both are recoverable
    ///        as dwelling internal gains.
    ///
    ///        Only the share of standing losses that reaches the heated space is
    ///        returned, given by the thermal loss recovery factor f_sto_m
    ///        (BS EN 15316-5:2017 Table B.3). This matches StorageTank, whose
    ///        recoverable storage losses already carry f_sto_m, so centralised
    ///        storage technologies are compared on a consistent basis. The remaining
    ///        losses escape the dwelling.
    ///
    ///        Returns:
    ///            Total recoverable losses across all units, in kWh.
    pub(crate) fn get_battery_losses(&self) -> f64 {
        let battery_losses = self.battery_losses.load(Ordering::SeqCst)
            * self.n_units as f64
            * THERMAL_CONSTANTS_F_STO_M;
        let pipework_gains = self.pipework_primary_gains_kwh.load(Ordering::SeqCst);
        self.pipework_primary_gains_kwh.store(0., Ordering::SeqCst);
        self.battery_losses.store(0., Ordering::SeqCst);
        battery_losses + pipework_gains
    }

    /// Return cumulative energy delivered to each service over the calculation.
    ///
    /// A running total, accumulated across timesteps and never reset, of the energy the
    /// battery delivers to each service (keyed by energy supply connection name). It is read
    /// after the simulation to apportion the battery's charging energy — metered on a single
    /// connection rather than per service — between the services in proportion to the output
    /// each received.
    ///
    /// Returns:
    ///     Mapping of service connection name to cumulative delivered energy (kWh), across all
    ///     units.
    pub(crate) fn energy_delivered_by_service(&self) -> IndexMap<ArcStr, f64> {
        self.energy_delivered_by_service.read().clone()
    }

    /// Calculate electric charging power for the current timestep (ChargeControl mode).
    ///
    /// Arguments
    /// * `simtime` - an iteration of the contextual simulation time
    ///
    /// In ChargeControl mode, returns the rated charge power when the
    /// ChargeControl is on, otherwise 0. In RangeTimeControl mode, always
    /// returns 0 because charging dispatch is handled per-source in timestep_end().
    ///   Returns:
    ///      Charging power in kW.
    fn electric_charge(&self, simtime: SimulationTimeIteration) -> f64 {
        if let Some(control) = &self.charge_control {
            if control.is_on(&simtime) {
                self.pwr_in
            } else {
                0.0
            }
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
    #[allow(clippy::unreadable_literal)]
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

    /// Return the energy transfer and starting temperature for a single zone.
    ///
    ///  The behaviour depends on the operation mode: charging fills zones from
    ///  the top down, standing losses remove heat in proportion to each zone's
    ///  temperature above the surroundings, and normal operation exchanges heat
    ///  with the water flowing through the heat exchanger.
    ///
    ///  Args:
    ///      index: Iteration index over the zones.
    ///      mode: Operation mode selecting the zone calculation.
    ///      zone_temp_c_dist: Current temperature of each zone, in °C.
    ///      inlet_temp_c: Temperature entering the battery (the surrounding air
    ///    temperature in the losses mode), in °C.
    ///      inlet_temp_c_zone: Temperature entering this zone, in °C.
    ///      Q_max_kJ: Energy available for transfer over the timestep, in kJ.
    ///      reynold_number_at_1_l_per_min: Reynolds number at 1 litre/minute.
    ///      flow_rate_kg_per_s: Heat exchanger flow rate, in kg/s.
    ///      time_step_s: Length of the calculation step, in seconds.
    ///
    ///  Returns:
    ///      A tuple of the energy transferred for this zone (in kJ), the zone
    ///      index, the zone's starting temperature (in °C), and the outlet
    ///      temperature (in °C).
    fn get_zone_properties(
        &self,
        index: usize,
        mode: &HeatBatteryPcmOperationMode,
        zone_temp_c_dist: &[f64],
        inlet_temp_c: f64,
        inlet_temp_c_zone: f64,
        q_max_kj: f64,
        reynold_number_at_1_l_per_min: f64,
        time_step_s: f64,
        flow_rate_l_per_min: f64,
        hex_a: Option<f64>,
        hex_b: Option<f64>,
    ) -> (f64, usize, f64, f64) {
        match mode {
            HeatBatteryPcmOperationMode::OnlyCharging => {
                let zone_index = zone_temp_c_dist.iter().len() - index - 1;
                let zone_temp_c_start = zone_temp_c_dist[zone_index];
                (0., zone_index, zone_temp_c_start, 0.)
            }
            HeatBatteryPcmOperationMode::Losses => {
                let zone_temp_c_start = zone_temp_c_dist[index];
                // Standing loss for this zone scales with its temperature above the
                // surroundings, relative to the temperature difference at which the
                // rated loss was characterised (BS EN 15316-5:2017, as for
                // StorageTank). Q_max_kJ carries the full rated-loss energy for the
                // timestep; share it equally across zones and weight by the per-zone
                // temperature ratio. Clamp at zero so a zone at or below the
                // surrounding temperature does not gain heat.
                let temp_diff_zone = zone_temp_c_start - inlet_temp_c;
                let energy_transf = (q_max_kj / zone_temp_c_dist.len() as f64 * temp_diff_zone
                    / self.temp_diff_rated_losses)
                    .max(0.0);
                (energy_transf, index, zone_temp_c_start, 0.)
            }
            HeatBatteryPcmOperationMode::Normal => {
                // NORMAL mode include battery primarily hydraulic charging or discharging with or without simultaneous electric charging
                let zone_temp_c_start = zone_temp_c_dist[index];
                let flow_rate_kg_per_s =
                    (flow_rate_l_per_min / SECONDS_PER_MINUTE as f64) * WATER.density();
                let effective_a = hex_a.unwrap_or(self.a);
                let effective_b = hex_b.unwrap_or(self.b);
                // The A * ln(Re * V) + B correlation returns the UA value for the
                // whole heat battery heat exchanger [W/K]. UA is an extensive
                // quantity (UA = U * area), so when the exchanger is discretised
                // into n_layers zones in series, each zone owns 1 / n_layers of the
                // surface area and therefore a conductance of UA / n_layers. Dividing
                // here keeps the total NTU (and hence the modelled heat transfer)
                // invariant to the chosen layer count. The zone storage capacities
                // are divided per layer in the same way (see __init__).
                let heat_transfer_kw_per_k = Self::calculate_heat_transfer_kw_per_k(
                    effective_a,
                    effective_b,
                    flow_rate_l_per_min,
                    reynold_number_at_1_l_per_min,
                ) / self.n_layers as f64;

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

    fn calculate_zone_energy_required(
        &self,
        zone_temp_c_start: f64,
        temp_charge_target: f64,
    ) -> f64 {
        if zone_temp_c_start >= self.phase_transition_temperature_upper {
            self.heat_storage_kj_per_k_above * (zone_temp_c_start - temp_charge_target)
        } else if zone_temp_c_start >= self.phase_transition_temperature_lower {
            if temp_charge_target > self.phase_transition_temperature_upper {
                self.heat_storage_kj_per_k_above
                    * (self.phase_transition_temperature_upper - temp_charge_target)
                    + self.heat_storage_kj_per_k_during
                        * (zone_temp_c_start - self.phase_transition_temperature_upper)
            } else {
                self.heat_storage_kj_per_k_during * (zone_temp_c_start - temp_charge_target)
            }
        } else if temp_charge_target > self.phase_transition_temperature_upper {
            self.heat_storage_kj_per_k_above
                * (self.phase_transition_temperature_upper - temp_charge_target)
                + self.heat_storage_kj_per_k_during
                    * (self.phase_transition_temperature_lower
                        - self.phase_transition_temperature_upper)
                + self.heat_storage_kj_per_k_below
                    * (zone_temp_c_start - self.phase_transition_temperature_lower)
        } else if temp_charge_target > self.phase_transition_temperature_lower {
            self.heat_storage_kj_per_k_during
                * (self.phase_transition_temperature_lower - temp_charge_target)
                + self.heat_storage_kj_per_k_below
                    * (zone_temp_c_start - self.phase_transition_temperature_lower)
        } else {
            self.heat_storage_kj_per_k_below * (zone_temp_c_start - temp_charge_target)
        }
    }

    fn process_zone_simultaneous_charging(
        &self,
        zone_temp_c_start: f64,
        temp_charge_target: f64,
        q_max_kj: f64,
        energy_transf: f64,
        energy_charged_electric: f64,
    ) -> (f64, f64, f64) {
        let mut q_max_kj = q_max_kj;
        let mut energy_charged_electric = energy_charged_electric;
        let mut energy_transf = energy_transf;

        if zone_temp_c_start < temp_charge_target {
            // zone initially below full charge
            let mut q_required =
                self.calculate_zone_energy_required(zone_temp_c_start, temp_charge_target);

            if energy_transf >= 0. {
                // inlet water withdraws energy from battery
                if -q_max_kj >= energy_transf {
                    // Charging is enough to recover energy withdrawn and possibly more
                    q_max_kj += energy_transf;
                    energy_charged_electric += energy_transf / KILOJOULES_PER_KILOWATT_HOUR as f64;
                    energy_transf = 0.;

                    if q_max_kj > q_required {
                        // Charging is not enough to push zone temperature to target
                        q_required = q_max_kj;
                        energy_charged_electric += -q_max_kj / KILOJOULES_PER_KILOWATT_HOUR as f64;
                        q_max_kj = 0.;
                    } else {
                        // Charging is enough to push zone temperature to target temperature
                        q_max_kj -= q_required;
                        energy_charged_electric +=
                            -q_required / KILOJOULES_PER_KILOWATT_HOUR as f64;
                    }
                    // Update zone temperature with energy from charging
                    energy_transf += q_required;
                } else {
                    // Charging can only recover partially the energy withdrawn
                    energy_transf += q_max_kj;
                    energy_charged_electric += -q_max_kj / KILOJOULES_PER_KILOWATT_HOUR as f64;
                    q_max_kj = 0.;
                }
            } else {
                // inlet water adds energy to battery
                if q_max_kj + energy_transf > q_required {
                    // inlet water + charging is not enough to push zone temperature to target
                    q_required = q_max_kj + energy_transf;
                    energy_charged_electric += -q_max_kj / KILOJOULES_PER_KILOWATT_HOUR as f64;
                    q_max_kj = 0.;

                    energy_transf = q_required;
                } else {
                    // inlet temperature + charging can take zone temperature to target temp
                    if energy_transf >= q_required {
                        // There is plenty of charging after taking zone temperature to target
                        q_max_kj -= q_required - energy_transf;
                        energy_charged_electric +=
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
                    energy_charged_electric += energy_transf / KILOJOULES_PER_KILOWATT_HOUR as f64;
                    energy_transf = 0.;
                } else {
                    // Charging can only recover partially the energy withdrawn
                    energy_transf += q_max_kj;
                    energy_charged_electric += -q_max_kj / KILOJOULES_PER_KILOWATT_HOUR as f64;
                    q_max_kj = 0.;
                }
            }
        }

        (q_max_kj, energy_charged_electric, energy_transf)
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
        if relative_eq!(energy_transf, 0.0, epsilon = NEGLIGIBLE_ENERGY_KJ) {
            // Negligible transfer: zone temperature is unchanged. Guarding the sign
            // test stops a noise-floor value taking the delivering vs retrieving
            // branch differently across platforms.
            return zone_temp_c_start;
        }
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
        time_step_s: f64,
        reynold_number_at_1_l_per_min: f64,
        flow_rate_l_per_min: f64,
        pwr_in: Option<f64>,
        mode: Option<HeatBatteryPcmOperationMode>,
        target_charge_fraction: f64,
        hex_a: Option<f64>,
        hex_b: Option<f64>,
    ) -> anyhow::Result<(f64, Vec<f64>, f64)> {
        let pwr_in = pwr_in.unwrap_or(0.);
        let mode = mode.unwrap_or(HeatBatteryPcmOperationMode::Normal);
        // target_charge_fraction is the SOC target (0–1) for charging, passed
        // by the caller. In ChargeControl mode this comes from ChargeControl.target_charge();
        // in RangeTimeControl mode from the source's upper setpoint.
        // Non-charging callers (discharge, losses) use the default 1.0.
        let target_temp = self.soc_to_temp(target_charge_fraction);
        let mut energy_charged_electric = 0.;

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
                    time_step_s,
                    flow_rate_l_per_min,
                    hex_a,
                    hex_b,
                );

            energy_transf_delivered[zone_index] += energy_transf;

            // Process energy transfer in zone with simultaneous charging
            if q_max_kj < 0. {
                (q_max_kj, energy_charged_electric, energy_transf) = self
                    .process_zone_simultaneous_charging(
                        zone_temp_c_start,
                        target_temp,
                        q_max_kj,
                        energy_transf,
                        energy_charged_electric,
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
            energy_transf_delivered,
            energy_charged_electric,
        ))
    }

    /// Charge the battery via hot water flow (hydronic charging).
    ///
    /// Simulates heat exchange between hot water flowing through the heat
    /// exchanger and the PCM zones. Iterates in sub-timesteps of
    /// __hb_time_step seconds, recalculating Reynolds number each iteration.
    /// Stops when:
    ///     - The outlet temperature reaches the inlet temperature (no further
    ///   heat transfer possible)
    ///     - The accumulated energy reaches energy_limit_kWh (conservation cap)
    ///     - Time runs out
    ///
    /// Args:
    ///     inlet_temp_C: Temperature of incoming hot water from the heat
    ///         source (°C). Typically the heat source flow temperature.
    ///     time_available_hrs: Time available for charging in this timestep
    ///         (hours). Accounts for time already spent on other services.
    ///     flow_rate_l_per_min: Flow rate through the charging heat exchanger
    ///         (litre/minute). From the hydronic source configuration.
    ///     target_charge_fraction: SOC target (0–1) from the source's
    ///         RangeTimeControl upper setpoint. Limits the temperature
    ///         target during zone heat exchange.
    ///     energy_limit_kWh: Maximum energy the battery can absorb (kWh),
    ///         based on the heat source capacity minus pipework losses.
    ///         Prevents the battery from absorbing more energy than the
    ///         heat source delivers.
    ///
    /// Returns:
    ///     Tuple of (energy_charged_kWh, updated_zone_temperatures).
    ///     energy_charged_kWh is positive when the battery absorbs energy.
    ///
    fn charge_battery_hydronic(
        &self,
        inlet_temp_c: f64,
        time_available_hrs: f64,
        flow_rate_l_per_min: f64,
        target_charge_fraction: f64,
        energy_limit_kwh: f64,
        hex_a: f64,
        hex_b: f64,
        hex_velocity_at_1_l_per_min: f64,
        hex_capillary_diameter_m: f64,
    ) -> anyhow::Result<(f64, Vec<f64>)> {
        let total_time_s = time_available_hrs * SECONDS_PER_HOUR as f64;
        let time_step_s = self.hb_time_step;

        // Initial Reynolds number
        let mut water_kinematic_viscosity_m2_per_s =
            Self::calculate_water_kinematic_viscosity_m2_per_s(
                self.initial_inlet_temp,
                self.estimated_outlet_temp,
            );
        let mut reynold_number_at_1_l_per_min = Self::calculate_reynold_number_at_1_l_per_min(
            water_kinematic_viscosity_m2_per_s,
            hex_velocity_at_1_l_per_min,
            hex_capillary_diameter_m,
        );
        let n_time_steps = (total_time_s / time_step_s) as usize;
        let mut zone_temp_c_dist = self.zone_temp_c_dist_initial.read().clone();
        let mut energy_absorbed_kj = 0.0;
        // Iterate through sub-timesteps, accumulating energy transferred to
        // the battery. energy_transf from __process_heat_battery_zones is
        // negative when the battery absorbs heat (water loses energy), so
        // we negate the sum to get positive energy_charged_kWh.
        let energy_limit_kj = energy_limit_kwh * units::KILOJOULES_PER_KILOWATT_HOUR as f64;
        for _ in 0..n_time_steps {
            let mut zone_temp_c_dist_prev = zone_temp_c_dist.clone();
            let (outlet_temp_c, energy_transf_charged, _) = self.process_heat_battery_zones(
                inlet_temp_c,
                &mut zone_temp_c_dist,
                time_step_s,
                reynold_number_at_1_l_per_min,
                flow_rate_l_per_min,
                Some(0.),
                Some(HeatBatteryPcmOperationMode::Normal),
                target_charge_fraction,
                Some(hex_a),
                Some(hex_b),
            )?;

            //  Recalculate Reynolds number for next sub-timestep
            water_kinematic_viscosity_m2_per_s =
                Self::calculate_water_kinematic_viscosity_m2_per_s(inlet_temp_c, outlet_temp_c);
            reynold_number_at_1_l_per_min = Self::calculate_reynold_number_at_1_l_per_min(
                water_kinematic_viscosity_m2_per_s,
                hex_velocity_at_1_l_per_min,
                hex_capillary_diameter_m,
            );

            // energy_transf_per_zone values are negative when the battery
            // absorbs heat from hot water (outlet cooler than inlet). Continue
            // charging only while the outlet is meaningfully below the inlet;
            // treating a near-equal outlet and inlet as "stop" keeps the
            // decision platform-independent at the noise floor.
            if outlet_temp_c < inlet_temp_c
                && !relative_eq!(
                    outlet_temp_c,
                    inlet_temp_c,
                    epsilon = NEGLIGIBLE_TEMP_DIFF_C
                )
            {
                let energy_this_step_kj = -FSum::with_all(&energy_transf_charged).value();
                energy_absorbed_kj += energy_this_step_kj;
                if energy_absorbed_kj > energy_limit_kj
                    || relative_eq!(energy_absorbed_kj, energy_limit_kj, epsilon = 1e-10)
                {
                    if energy_this_step_kj > 0. {
                        let energy_allowed_kj =
                            energy_limit_kj - (energy_absorbed_kj - energy_this_step_kj);
                        let fraction = energy_allowed_kj / energy_this_step_kj;
                        let (_, _, _) = self.process_heat_battery_zones(
                            inlet_temp_c,
                            zone_temp_c_dist_prev.as_mut_slice(),
                            time_step_s * fraction,
                            reynold_number_at_1_l_per_min,
                            flow_rate_l_per_min,
                            Some(0.),
                            Some(HeatBatteryPcmOperationMode::Normal),
                            target_charge_fraction,
                            Some(hex_a),
                            Some(hex_b),
                        )?;
                        zone_temp_c_dist = zone_temp_c_dist_prev;
                        energy_absorbed_kj = energy_limit_kj;
                    }
                    break;
                }
            } else {
                break;
            }
        }

        // Convert kJ to kWh (positive = energy absorbed by battery)
        let energy_charged_kwh = energy_absorbed_kj / units::KILOJOULES_PER_KILOWATT_HOUR as f64;
        self.energy_charged_total
            .fetch_add(energy_charged_kwh, Ordering::SeqCst);
        *self.zone_temp_c_dist_initial.write() = zone_temp_c_dist.clone();

        Ok((energy_charged_kwh, zone_temp_c_dist))
    }

    /// Calculate the energy demand to charge the battery from current SOC.
    ///
    /// Uses the SOC deficit and the battery's energy capacity to estimate
    /// the energy needed in kWh. When temp_flow_max is provided (hydronic
    /// charging), the target temperature is capped to avoid requesting more
    /// energy than the heat source can thermodynamically deliver. The
    /// target_soc parameter limits how far to charge — preventing
    /// overcharging beyond the RangeTimeControl upper setpoint within a
    /// single timestep.
    ///
    /// Args:
    ///     temp_flow_max: Maximum flow temperature the heat source can
    ///         deliver (°C). When provided, caps the charge target at
    ///         min(temp_flow_max, temp_charge_max). Without this cap, a heat
    ///         source with a lower flow temperature than the battery's max
    ///         would be billed for energy it cannot deliver, breaking
    ///         conservation of energy.
    ///     target_soc: SOC target (0–1) from the source's RangeTimeControl
    ///         upper setpoint. Energy demand is capped at the deficit to
    ///         reach this SOC rather than full charge.
    ///
    /// Returns:
    ///     Energy demand in kWh (positive). Returns 0.0 if battery is at or
    ///     above the target SOC.
    fn calculate_charge_energy_demand(
        &self,
        temp_flow_max: Option<f64>,
        target_soc: f64,
    ) -> anyhow::Result<f64> {
        let current_soc =
            self.calc_state_of_charge(self.zone_temp_c_dist_initial.read().clone())?;
        if current_soc > target_soc || relative_eq!(current_soc, target_soc, epsilon = 1e-10) {
            return Ok(0.0);
        }

        let temp_ref = self.temp_ref;
        let energy_current_kj = self.energy_stored_max * current_soc;

        // Demand based on SOC target (relative to full battery capacity at
        // temp_charge_max). This is how much energy is needed to reach the
        // target SOC regardless of flow temperature limitations.
        let energy_at_target_soc_kj = self.energy_stored_max * target_soc;
        let demand_soc_kj = energy_at_target_soc_kj - energy_current_kj;

        // Demand based on flow temperature limit (maximum energy the heat
        // source can thermodynamically deliver). When temp_flow_max <
        // temp_charge_max, the battery cannot be heated beyond temp_flow_max
        // even if the SOC target is higher.
        let energy_deficit_kj = if let Some(temp_flow_max) = temp_flow_max {
            let energy_at_flow_max_kj = HeatBatteryPcm::calculate_layer_energy_stored(
                temp_flow_max,
                temp_ref,
                self.phase_transition_temperature_lower,
                self.phase_transition_temperature_upper,
                self.heat_storage_kj_per_k_below,
                self.heat_storage_kj_per_k_during,
                self.heat_storage_kj_per_k_above,
            ) * self.n_layers as f64;
            let demand_flow_kj = (energy_at_flow_max_kj - energy_current_kj).max(0.0);
            demand_soc_kj.min(demand_flow_kj)
        } else {
            demand_soc_kj
        };

        Ok(energy_deficit_kj / units::KILOJOULES_PER_KILOWATT_HOUR as f64)
    }

    /// It follows the same methodology as energy_demand function
    /// Charge the battery electrically at the given rated power.
    ///
    /// Applies electric heating to all battery zones using the
    /// `OnlyCharging` operation mode. The charge target limits the zone
    /// temperatures so the battery does not charge beyond the active control
    /// setpoint during this timestep.
    ///
    /// # Arguments
    ///
    /// * `rated_power` - Electric charging power in kW.
    /// * `target_charge_fraction` - SOC target in the range 0–1. This limits
    ///   the temperature target during zone heat exchange.
    ///
    /// # Returns
    ///
    /// A tuple containing:
    ///
    /// * The energy charged during this timestep in kWh.
    /// * The updated zone temperatures.
    ///
    /// # Errors
    ///
    /// Returns an error if processing the battery zones fails.
    fn charge_battery_electric(
        &self,
        rated_power: f64,
        target_charge_fraction: f64,
    ) -> anyhow::Result<(f64, Vec<f64>)> {
        let timestep = self.simulation_time_step;
        let time_available = self.time_available(0., timestep);
        let time_step_s = time_available * SECONDS_PER_HOUR as f64;

        let mut zone_temp_c_dist = self.zone_temp_c_dist_initial.read().clone();

        // Apply the electric charging power across all battery zones.
        let (_, _, energy_charged_electric_substep) = self.process_heat_battery_zones(
            0.,
            &mut zone_temp_c_dist,
            time_step_s,
            0.,
            self.flow_rate_l_per_min,
            Some(rated_power),
            Some(HeatBatteryPcmOperationMode::OnlyCharging),
            target_charge_fraction,
            None,
            None,
        )?;

        // The returned energy is already in kWh.
        self.energy_charged_total
            .fetch_add(energy_charged_electric_substep, Ordering::SeqCst);

        self.energy_charged_electric
            .fetch_add(energy_charged_electric_substep, Ordering::SeqCst);

        *self.zone_temp_c_dist_initial.write() = zone_temp_c_dist.clone();

        Ok((energy_charged_electric_substep, zone_temp_c_dist))
    }

    /// Unified charging entry point for both control modes.
    ///
    /// Dispatches to the appropriate charging path based on the battery's
    /// configuration:
    /// - RangeTimeControl mode (HeatSource dict): per-source SOC-based
    ///   charging with hysteresis, supporting both electric and hydronic sources.
    /// - ChargeControl mode: single electric element controlled by
    ///   ChargeControl.is_on() and target_charge().
    ///
    /// Args:
    ///     time_remaining_current_timestep: Time remaining in the current
    ///         timestep (hours), after services have consumed their share.
    ///
    /// Returns:
    ///      Tuple of (energy_charged_kWh, updated_zone_temperatures).
    fn charge_battery(
        &self,
        time_remaining_current_timestep: f64,
        simtime: &SimulationTimeIteration,
    ) -> anyhow::Result<(f64, Vec<f64>)> {
        let zone_temp_c_after_charging = self.zone_temp_c_dist_initial.read().clone();

        if self.use_heatsource_data {
            // RangeTimeControl: per-source charging dispatch (SOC or temperature-based).
            // The order of processing heat sources shouldn't matter because
            // their schedules should not overlap.
            for (src, source) in self.heat_source_data.as_ref().unwrap_or(&IndexMap::new()) {
                self.determine_heat_source_switch_on(src, source, simtime)?;
                self.determine_heat_source_switch_off(src, source, simtime)?;
                if *self.charging_active.read().get(src).unwrap_or(&false) {
                    // Use the upper setpoint from the source's RangeTimeControl
                    // as the SOC target for charging. This limits both the
                    // temperature target in zone heat exchange and the energy
                    // demand calculation, preventing overcharging beyond the
                    // hysteresis upper threshold within a single timestep.
                    let (_, setpnt_upper) = self.resolve_setpoints(source, simtime)?;
                    let target_soc = if let Some(setpnt_upper) = setpnt_upper {
                        setpnt_upper
                    } else {
                        bail!("Upper setpoint should be guaranteed by switch_off logic");
                    };
                    match source.source_type {
                        ChargingSourceType::DirectElectric => {
                            if let Some(rated_charge_power) = source.rated_charge_power {
                                return self
                                    .charge_battery_electric(rated_charge_power, target_soc);
                            } else {
                                // maybe this is already enforced elsewhere
                                bail!("Rated charge power must not be None for DirectElectric charging source");
                            }
                        }
                        ChargingSourceType::HeatSourceWet => {
                            if source.heat_source_service.is_some() {
                                return self.charge_from_heat_source(
                                    source,
                                    time_remaining_current_timestep,
                                    target_soc,
                                    simtime,
                                );
                            } else {
                                bail!("Heat source service must not be None for HeatSourceWet charging source");
                            }
                        }
                    }
                } else {
                    if source.source_type == ChargingSourceType::HeatSourceWet {
                        // Hydronic source not active — call demand_energy(0) so
                        // the heat source records a zero-demand service call
                        // every timestep (ensures consistent detailed results
                        // without placeholders)
                        if let HeatBatteryChargingSource {
                            heat_source_service: Some(heat_source_service),
                            temp_flow_max,
                            flow_rate_charging_l_per_min: Some(flow_rate_charging_l_per_min),
                            hex_a: Some(hex_a),
                            hex_b: Some(hex_b),
                            hex_velocity_at_1_l_per_min: Some(hex_velocity_at_1_l_per_min),
                            hex_capillary_diameter_m: Some(hex_capillary_diameter_m),
                            ..
                        } = source
                        {
                            let temp_return = self.estimate_return_temp(
                                *temp_flow_max,
                                *flow_rate_charging_l_per_min,
                                *hex_a,
                                *hex_b,
                                *hex_velocity_at_1_l_per_min,
                                *hex_capillary_diameter_m,
                            )?;

                            heat_source_service.demand_energy(
                                0.0,
                                *temp_flow_max,
                                temp_return,
                                simtime,
                            )?;
                            // Report zero input to pipework loss tracker for event
                            // boundary detection
                            self.pipework.calculate_primary_pipework_losses(
                                0.0,
                                Some(source.temp_flow_max),
                                Some(true),
                                simtime,
                            )?;
                        }
                    }
                }
            }
            Ok((0., zone_temp_c_after_charging))
        } else if let Some(charge_control) = &self.charge_control {
            if charge_control.is_on(simtime) {
                // ChargeControl: single electric element with temperature-based proxy
                let pwr_in = self.electric_charge(*simtime);
                let target = charge_control.target_charge(*simtime, None)?;
                return self.charge_battery_electric(pwr_in, target);
            }
            Ok((0., zone_temp_c_after_charging))
        } else {
            bail!("No charging configuration: neither HeatSource dict nor ChargeControl is configured.")
        }
    }

    /// Charge the battery from a wet heat source (e.g. heat pump).
    ///
    /// Iterative temperature refinement algorithm:
    /// 1. Estimate return temperature via one-pass heat exchange simulation
    ///    on a copy of zone state (non-mutating)
    /// 2. Query energy_output_max() (side-effect-free) to get heat source capacity
    ///    at estimated flow/return temperatures
    /// 3. Calculate battery energy demand, cap to heat source capacity
    /// 4. Make a single demand_energy() call with refined temperatures
    ///
    /// Args:
    ///     source: HeatBatteryChargingSource with heat_source_service,
    ///         temp_flow_max, and control.
    ///     time_available_hrs: Remaining timestep time in hours.
    ///     target_charge_fraction: SOC target (0–1) from the source's
    ///         RangeTimeControl upper setpoint.
    ///
    /// Returns:
    ///     Tuple of (energy_charged_kWh, updated_zone_temperatures).
    ///
    fn charge_from_heat_source(
        &self,
        source: &HeatBatteryChargingSource,
        time_available_hrs: f64,
        target_charge_fraction: f64,
        simtime: &SimulationTimeIteration,
    ) -> anyhow::Result<(f64, Vec<f64>)> {
        if let HeatBatteryChargingSource {
            heat_source_service: Some(heat_source_service),
            flow_rate_charging_l_per_min: Some(flow_rate_charging_l_per_min),
            hex_a: Some(hex_a),
            hex_b: Some(hex_b),
            hex_velocity_at_1_l_per_min: Some(hex_velocity_at_1_l_per_min),
            hex_capillary_diameter_m: Some(hex_capillary_diameter_m),
            ..
        } = source
        {
            // Charge at the flow temperature needed to reach the SOC target plus a
            // heat-exchanger approach difference, not the source's peak.
            // __soc_to_temp gives the uniform PCM temperature for the target charge
            // fraction; the approach difference above it drives heat across the
            // exchanger so the store can actually reach the target rather than only
            // approach it asymptotically. Querying a wet source (e.g. a heat pump)
            // at temp_flow_max would understate its COP and capacity. The source
            // cannot deliver above its own maximum, so bound by that.

            //           let   temp_flow = min(
            // +            source.temp_flow_max,
            // +            self.__soc_to_temp(target_charge_fraction) + self.CHARGE_APPROACH_TEMP_DIFF_C,
            // +        )
            let temp_flow = (self.soc_to_temp(target_charge_fraction)
                + CHARGE_APPROACH_TEMP_DIFF_C)
                .min(source.temp_flow_max);
            let temp_return = self.estimate_return_temp(
                temp_flow,
                *flow_rate_charging_l_per_min,
                *hex_a,
                *hex_b,
                *hex_velocity_at_1_l_per_min,
                *hex_capillary_diameter_m,
            )?;
            // Calculate how much energy one unit needs, capped to what the
            // heat source can deliver at its max flow temperature and limited
            // to the SOC target from the source's RangeTimeControl.
            let energy_demand_per_unit = self.calculate_charge_energy_demand(
                Some(source.temp_flow_max),
                target_charge_fraction,
            )?;

            if energy_demand_per_unit <= 0.
                || relative_eq!(energy_demand_per_unit, 0., epsilon = 1e-10)
            {
                // No charging needed — still call demand_energy(0) so the heat
                // source records a zero-demand service call, and run pipework loss
                // calc to detect end-of-event transition.
                heat_source_service.demand_energy(
                    0.,
                    source.temp_flow_max,
                    temp_return,
                    simtime,
                )?;
                self.pipework.calculate_primary_pipework_losses(
                    0.0,
                    Some(source.temp_flow_max),
                    Some(true),
                    simtime,
                )?;
                return Ok((0., self.zone_temp_c_dist_initial.read().clone()));
            }
            // Scale to total demand across all units for heat source interaction
            let energy_demand_total = energy_demand_per_unit * self.n_units as f64;

            // Use energy_output_max() to check heat source capacity at the estimated
            // temperatures. This is side-effect-free (no energy consumption
            // recorded), giving the maximum energy the heat source can deliver.
            let energy_max =
                heat_source_service.energy_output_max(temp_flow, temp_return, simtime)?;

            if relative_eq!(energy_max, 0.0, epsilon = 1e-10) {
                // Heat source has no capacity — call demand_energy(0) for reporting
                heat_source_service.demand_energy(0., temp_flow, temp_return, simtime)?;
                self.pipework.calculate_primary_pipework_losses(
                    0.,
                    Some(temp_flow),
                    Some(true),
                    simtime,
                )?;
                return Ok((0., self.zone_temp_c_dist_initial.read().clone()));
            } else if energy_max < 0. {
                bail!(
                    "energy_output_max returned negative value ({energy_max}). This indicates a bug in the heat source implementation."
                );
            } else {
                // Cap demand to heat source capacity before calculating pipework losses,
                // so losses are based on the achievable energy transfer, not the
                // uncapped demand.
                let energy_demand_total = energy_demand_total.min(energy_max);

                // Account for primary pipework losses: the heat source must deliver
                // extra energy to compensate for losses in the pipes between it and
                // the battery (following StorageTank's heat_source_output pattern)
                let (pipework_losses_kwh, primary_gains_w) =
                    self.pipework.calculate_primary_pipework_losses(
                        energy_demand_total,
                        Some(temp_flow),
                        Some(true),
                        simtime,
                    )?;

                // Energy available for the battery is what the heat source can
                // deliver minus what's lost in the pipework.
                let energy_for_battery =
                    (energy_max - pipework_losses_kwh).min(energy_demand_total);
                let energy_for_battery_per_unit = energy_for_battery / self.n_units as f64;

                // Simulate the physical heat exchange with energy conservation cap.
                // The simulation stops when the battery has absorbed
                // energy_for_battery_per_unit, preventing it from absorbing more
                // than the heat source delivers minus pipework losses.
                let (energy_charged, zone_temps) = self.charge_battery_hydronic(
                    temp_flow,
                    time_available_hrs,
                    *flow_rate_charging_l_per_min,
                    target_charge_fraction,
                    energy_for_battery_per_unit,
                    *hex_a,
                    *hex_b,
                    *hex_velocity_at_1_l_per_min,
                    *hex_capillary_diameter_m,
                )?;

                // Bill the heat source for the actual energy absorbed by the battery
                // plus pipework losses. This ensures the HP is only charged for what
                // was actually used, not the pre-estimated demand.
                let energy_hp_demand = energy_charged * self.n_units as f64 + pipework_losses_kwh;
                heat_source_service.demand_energy(
                    energy_hp_demand,
                    temp_flow,
                    temp_return,
                    simtime,
                )?;

                // Accumulate internal pipework gains for return via get_battery_losses().
                // Convert W to kWh for consistency with __battery_losses (both are in kWh).
                if primary_gains_w > 0. {
                    self.pipework_primary_gains_kwh.fetch_add(
                        primary_gains_w / WATTS_PER_KILOWATT as f64 * self.simulation_time_step,
                        Ordering::SeqCst,
                    );
                }
                Ok((energy_charged, zone_temps))
            }
        } else {
            bail!("Incomplete heat source configuration, missing required fields. heat_source_service, temp_flow_max, flow_rate_charging_l_per_min, hex_a, hex_b, hex_velocity_at_1_l_per_min, hex_capillary_diameter_m are required.")
        }
    }
    /// Activate charging for a source when SOC falls to/below the lower setpoint.
    ///
    /// Follows StorageTank._determine_heat_source_switch_on pattern. The
    /// RangeTimeControl provides (lower, upper) setpoints for hysteresis,
    /// resolved to SOC via __resolve_setpoints (which converts temperature
    /// setpoints to SOC when schedule_unit="temperature"):
    /// - Both non-None → active period: start charging when SOC <= lower
    /// - lower is None, upper non-None → transition period: don't start new
    ///   charging (existing active state persists until switch-off)
    /// - Both None → off period: don't start charging
    ///
    /// Args:
    ///     source_name: Key into __charging_active for this source.
    ///     source: HeatBatteryChargingSource with control and config.
    fn determine_heat_source_switch_on(
        &self,
        source_name: &str,
        source: &HeatBatteryChargingSource,
        simtime: &SimulationTimeIteration,
    ) -> anyhow::Result<()> {
        let (setpnt_lower, setpnt_upper) = self.resolve_setpoints(source, simtime)?;
        if setpnt_upper.is_some() {
            let soc = self.calc_state_of_charge(self.zone_temp_c_dist_initial.read().clone())?;
            if let Some(setpnt_lower) = setpnt_lower {
                if soc <= setpnt_lower || relative_eq!(soc, setpnt_lower, epsilon = 1e-10) {
                    // Active period: start charging when SOC <= lower
                    self.charging_active
                        .write()
                        .insert(source_name.into(), true);
                    return Ok(());
                }
            }
        }
        // Off period — control is entirely off, don't start charging
        Ok(())
    }
    /// Deactivate charging for a source when SOC reaches the upper setpoint.
    ///
    /// Follows StorageTank._determine_heat_source_switch_off pattern.
    /// Setpoints are resolved to SOC via __resolve_setpoints (converting
    /// from temperature when schedule_unit="temperature"). Turns off when:
    /// - setpnt_upper is None (off period), or
    /// - SOC >= upper setpoint (battery sufficiently charged)
    ///
    /// Args:
    ///     source_name: Key into __charging_active for this source.
    ///     source: HeatBatteryChargingSource with control and config.
    ///
    fn determine_heat_source_switch_off(
        &self,
        source_name: &str,
        source: &HeatBatteryChargingSource,
        simtime: &SimulationTimeIteration,
    ) -> anyhow::Result<()> {
        let (_, setpnt_upper) = self.resolve_setpoints(source, simtime)?;
        let soc = self.calc_state_of_charge(self.zone_temp_c_dist_initial.read().clone())?;
        if setpnt_upper.is_none_or(|upper| soc > upper || relative_eq!(soc, upper, epsilon = 1e-10))
        {
            self.charging_active
                .write()
                .insert(source_name.into(), false);
        }
        Ok(())
    }

    /// Estimate the return temperature from a single heat exchange pass.
    ///
    /// Runs one sub-timestep of zone-by-zone heat exchange on a copy of
    /// the current zone temperatures (non-mutating). The outlet temperature
    /// from this pass reflects the heat exchanger effectiveness, flow rate,
    /// and current zone thermal state — giving a physics-based estimate of
    /// what temperature the water returns to the heat source.
    ///
    /// Args:
    ///     temp_flow: Inlet temperature from the heat source (°C).
    ///     flow_rate_l_per_min: Flow rate through the charging heat exchanger
    ///         (litre/minute).
    ///
    /// Returns:
    ///     Estimated outlet (return) temperature in °C.
    fn estimate_return_temp(
        &self,
        temp_flow: f64,
        flow_rate_charging_l_per_min: f64,
        hex_a: f64,
        hex_b: f64,
        hex_velocity_at_1_l_per_min: f64,
        hex_capillary_diameter_m: f64,
    ) -> anyhow::Result<f64> {
        let water_kinematic_viscosity_m2_per_s =
            HeatBatteryPcm::calculate_water_kinematic_viscosity_m2_per_s(
                self.initial_inlet_temp,
                self.estimated_outlet_temp,
            );
        let reynold_number_at_1_l_per_min = HeatBatteryPcm::calculate_reynold_number_at_1_l_per_min(
            water_kinematic_viscosity_m2_per_s,
            hex_velocity_at_1_l_per_min,
            hex_capillary_diameter_m,
        );

        // Run one sub-timestep on a copy of zone temps (non-mutating)
        let (temp_outlet_c, _, _) = self.process_heat_battery_zones(
            temp_flow,
            self.zone_temp_c_dist_initial.read().clone().as_mut_slice(),
            self.hb_time_step,
            reynold_number_at_1_l_per_min,
            flow_rate_charging_l_per_min,
            Some(0.0),
            Some(HeatBatteryPcmOperationMode::Normal),
            0.,
            Some(hex_a),
            Some(hex_b),
        )?;
        Ok(temp_outlet_c)
    }

    /// Calculate the standing heat loss over the timestep and update zone temperatures.
    ///
    /// Losses are applied per zone in proportion to each zone's temperature
    /// above the assumed surrounding air temperature, scaled by the temperature
    /// difference at which the rated loss was characterised.
    ///
    /// Returns:
    ///     A tuple of the total standing loss over the timestep (in kWh) and
    ///     the updated zone temperatures (in °C).
    ///
    fn battery_heat_loss(&self) -> anyhow::Result<(f64, Vec<f64>)> {
        // Battery losses
        let timestep = self.simulation_time_step;
        let time_step_s = timestep * SECONDS_PER_HOUR as f64; // time_available * SECONDS_PER_HOUR;

        let mut zone_temp_c_dist = self.zone_temp_c_dist_initial.read().clone();

        // Processing HB zones. The surrounding air temperature is assumed fixed,
        // consistent with the hot water cylinder standby-loss calculation.
        let (_, energy_loss, _) = self.process_heat_battery_zones(
            DEFAULT_AMBIENT_TEMPERATURE,
            &mut zone_temp_c_dist,
            time_step_s,
            time_step_s,
            0.,
            Some(-self.max_rated_losses),
            Some(HeatBatteryPcmOperationMode::Losses),
            0.,
            None,
            None,
        )?;

        *self.zone_temp_c_dist_initial.write() = zone_temp_c_dist.clone();

        //Equivalent of using Python's math.fsum instead of sum() for better numerical accuracy with floating point arithmetic
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
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let total_time_s = volume / self.flow_rate_l_per_min * SECONDS_PER_MINUTE as f64;

        let time_step_s = self.hb_time_step;

        let pwr_in = self.electric_charge(simtime);

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

        let mut zone_temp_c_dist = self.zone_temp_c_dist_initial.read().clone();
        let mut inlet_temp_c = inlet_temp;
        let mut outlet_temp_c = inlet_temp_c; // initialise, though expectation is this will be overridden in loop

        let n_time_steps = if total_time_s > time_step_s {
            (total_time_s / time_step_s) as usize
        } else {
            1
        };

        for _ in 0..n_time_steps {
            (outlet_temp_c, _, _) = self.process_heat_battery_zones(
                inlet_temp_c,
                &mut zone_temp_c_dist,
                time_step_s,
                reynold_number_at_1_l_per_min,
                self.flow_rate_l_per_min,
                Some(pwr_in),
                None,
                0.,
                None, // Todo - temp values as part of 1.0.0a9
                None, // Todo - temp values as part of 1.0.0a9
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

    /// Return the maximum energy the battery can deliver over the timestep.
    ///
    /// The heat-exchanger inlet is held constant at the return-feed temperature for the
    /// HEM timestep and the full positive heat transfer to the water is summed, stopping
    /// only once the core has cooled to where it would absorb heat from the flow rather
    /// than deliver it. As in demand_energy, the required flow temperature does not gate
    /// delivery: the heat transferred to the loop is driven by the return-feed inlet and
    /// the core state, so this returns the same energy demand_energy delivers when given
    /// unlimited demand.
    ///
    /// Args:
    ///     temp_output: required emitter flow temperature, in °C. Accepted for call-site
    ///         symmetry with demand_energy; it does not gate the deliverable maximum.
    ///     temp_return_feed: heat-exchanger inlet temperature, in °C
    ///     time_start: start time within the timestep, in hours
    ///
    /// Returns:
    ///     Maximum deliverable energy across all units, in kWh.
    fn energy_output_max(
        &self,
        _temp_output: f64,
        temp_return_feed: f64,
        time_start: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let time_start = time_start.unwrap_or(0.);
        let timestep = self.simulation_time_step;
        let time_available = self.time_available(time_start, timestep);

        let total_time_s = time_available * SECONDS_PER_HOUR as f64;
        // Integrate on the same sub-timestep demand_energy caps at, so the maximum equals
        // what demand_energy delivers. A coarser step degrades accuracy because the
        // Reynolds number is held over intervals where the fluid properties have changed
        // enough to matter, leaving the maximum offset from the deliverable energy.
        let time_step_s = self.hb_time_step;
        // Only charge while discharging if the battery supports it, matching
        // demand_energy. Otherwise charging is deferred to timestep_end, so
        // assuming it here would overstate the deliverable ceiling that
        // demand_energy can reach.

        let pwr_in = if self.simultaneous_charging_and_discharging {
            self.electric_charge(simtime)
        } else {
            0.
        };

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

        let mut zone_temp_c_dist = self.zone_temp_c_dist_initial.read().deref().clone();
        let mut energy_delivered_hb = 0.;
        let mut inlet_temp_c = temp_return_feed;
        let n_time_steps = (total_time_s / time_step_s) as usize;

        for _ in 0..n_time_steps {
            // Processing HB zones
            let (outlet_temp_c, energy_transf_delivered, _) = self.process_heat_battery_zones(
                inlet_temp_c,
                &mut zone_temp_c_dist,
                time_step_s,
                reynold_number_at_1_l_per_min,
                self.flow_rate_l_per_min,
                Some(pwr_in),
                None,
                0.0,
                None,
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

            // Equivalent of using Python's math.fsum instead of sum() for better numerical accuracy with floating point arithmetic
            let energy_delivered_kj = FSum::with_all(&energy_transf_delivered).value();
            let energy_delivered_ts = energy_delivered_kj / KILOJOULES_PER_KILOWATT_HOUR as f64;

            // Stop before counting a sub-step in which the battery would absorb heat from
            // the inlet flow rather than deliver it, matching demand_energy. The required
            // flow temperature does not gate delivery, so the maximum is the full positive
            // heat transfer the core can drive from the return-feed inlet. A negligibly-
            // negative result is floating-point noise and is treated as zero so the stop
            // decision is identical across platforms.
            if energy_delivered_kj < 0.
                && !relative_eq!(energy_delivered_kj, 0.0, epsilon = NEGLIGIBLE_ENERGY_KJ)
            {
                break;
            }

            // In this new method, adjust total energy to make more real with the 6 ts we have configured
            energy_delivered_hb += energy_delivered_ts;

            inlet_temp_c = temp_return_feed
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
        temp_return_feed: Option<f64>,
        temp_output: Option<f64>,
        service_on: bool,
        time_start: Option<f64>,
        update_heat_source_state: Option<bool>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // Return the energy provided by the HB during a HEM time step (assuming
        // an inlet temperature constant) and update the HB state (zones distribution temperatures)
        // The HEM time step is divided into sub-timesteps. For each sub-timestep the zones temperature are
        // calculated (loop through zones).
        let time_start = time_start.unwrap_or(0.);
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);
        let timestep = self.simulation_time_step;
        let time_available = self.time_available(time_start, timestep);

        // demand_energy is called for each service in each timestep
        // Some calculations are only required once per timestep
        // Perform these calculations here
        if self.flag_first_call.load(Ordering::SeqCst) {
            self.first_call();
        }

        let pwr_in = if self.simultaneous_charging_and_discharging {
            self.electric_charge(simtime)
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

        let mut energy_delivered_hb = 0.;
        // inlet_temp_c assignment moved down in comparison with Python, as we have to deal with None case
        let mut zone_temp_c_dist = self.zone_temp_c_dist_initial.read().clone();

        if energy_output_required < 0.
            || relative_eq!(
                energy_output_required,
                0.,
                max_relative = 1e-09,
                epsilon = 1e-10
            )
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
                    current_hb_power: ResultParamValue::Empty,
                });
            }
            return Ok(0.);
        }

        let temp_return_feed = temp_return_feed.ok_or_else(|| anyhow!("temp_return_feed value was expected to be set for demand_energy method on HeatBatteryPcm when energy_output_required > 0"))?;
        let inlet_temp_c = temp_return_feed;

        let mut time_step_s = 1.;
        let mut time_running_current_service = 0.;

        let mut energy_charged_electric_service = 0.;

        let mut outlet_temp_c = None;

        while time_step_s > 0. {
            // Processing HB zones
            let (outlet_temp_c_new, energy_transf_delivered, energy_charged_electric_substep) =
                self.process_heat_battery_zones(
                    inlet_temp_c,
                    &mut zone_temp_c_dist,
                    time_step_s,
                    reynold_number_at_1_l_per_min,
                    self.flow_rate_l_per_min,
                    Some(pwr_in),
                    None,
                    0.0,
                    None,
                    None,
                )?;

            outlet_temp_c = Some(outlet_temp_c_new);

            if update_heat_source_state {
                self.energy_charged_electric
                    .fetch_add(energy_charged_electric_substep, Ordering::SeqCst);
            }
            energy_charged_electric_service += energy_charged_electric_substep;

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

            let energy_delivered_kj = FSum::with_all(&energy_transf_delivered).value();

            //  A negligibly-negative result is floating-point noise, not real
            //  absorption, so it is treated as zero to keep the stop decision
            //  identical across platforms.
            if energy_delivered_kj < 0. {
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
                energy_delivered_kj / KILOJOULES_PER_KILOWATT_HOUR as f64;
            energy_delivered_hb += energy_delivered_ts; // demand_per_time_step_kwh
                                                        // balance = total_energy - energy_charged
            let max_instant_power = energy_delivered_ts / time_step_s;

            if max_instant_power > 0. {
                time_step_s = (energy_demand - energy_delivered_hb) / max_instant_power;
            }

            if time_step_s > self.hb_time_step {
                time_step_s = self.hb_time_step;
            }

            if relative_eq!(
                energy_demand,
                energy_delivered_hb,
                max_relative = 1e-09,
                epsilon = 1e-10
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

            // Track pump running time (only for regular DHW and space heating services)
            // Direct DHW services don't use circulation pumps
            match service_type {
                HeatingServiceType::DomesticHotWaterRegular | HeatingServiceType::Space => {
                    self.pump_running_time_current_timestep.fetch_add(
                        time_running_current_service / SECONDS_PER_HOUR as f64,
                        Ordering::SeqCst,
                    );
                }
                HeatingServiceType::DomesticHotWaterDirect => {
                    // Direct DHW doesn't use circulation pump but track time
                    // separately — it uses a different heat exchanger and can
                    // happen simultaneously with hydronic charging.
                    self.time_running_direct_current_timestep.fetch_add(
                        time_running_current_service / SECONDS_PER_HOUR as f64,
                        Ordering::SeqCst,
                    );
                } // Direct DHW doesn't use circulation pump
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
                temp_inlet: Some(temp_return_feed),
                time_running: time_running_current_service,
                energy_delivered_hb: energy_delivered_hb * self.n_units as f64,
                energy_delivered_backup: 0.,
                energy_delivered_total: energy_delivered_hb * self.n_units as f64 + 0.,
                energy_charged_during_service: energy_charged_electric_service
                    * self.n_units as f64,
                hb_zone_temperatures: zone_temp_c_dist,
                current_hb_power: ResultParamValue::Number(current_hb_power * self.n_units as f64),
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
    pub(crate) fn timestep_end(&self, simtime: SimulationTimeIteration) -> anyhow::Result<()> {
        let timestep = self.simulation_time_step;
        let time_remaining_current_timestep = timestep
            - self
                .total_time_running_current_timestep
                .load(Ordering::SeqCst);
        // Direct DHW uses a separate heat exchanger from hydronic charging,
        // so time spent on direct DHW is still available for charging.
        let time_remaining_for_charging = self
            .time_running_direct_current_timestep
            .fetch_add(time_remaining_current_timestep, Ordering::SeqCst);

        if self.flag_first_call.load(Ordering::SeqCst) {
            self.first_call();
        }
        self.flag_first_call.store(true, Ordering::SeqCst);

        // Calculating auxiliary energy to provide services during timestep
        let energy_aux =
            self.calc_auxiliary_energy(timestep, time_remaining_current_timestep, simtime.index)?;

        let (battery_losses, zone_temp_c_after_losses) = self.battery_heat_loss()?;
        self.battery_losses.store(battery_losses, Ordering::SeqCst);

        // Charging battery for the remaining of the timestep. Direct DHW
        // time is added back because it uses a separate heat exchanger and
        // can happen simultaneously with hydronic charging.
        let (end_of_ts_charge, zone_temp_c_after_charging) = if self
            .charge_control
            .clone()
            .is_some_and(|cc| cc.is_on(&simtime))
        {
            self.charge_battery(time_remaining_for_charging, &simtime)?
        } else {
            (0., self.zone_temp_c_dist_initial.read().clone())
        };

        self.energy_supply_connection.demand_energy(
            self.energy_charged_electric.load(Ordering::SeqCst) * self.n_units as f64,
            simtime.index,
        )?;

        // Calculate SOC at end of timestep (after charging and losses)
        let soc_end_of_timestep =
            self.calc_state_of_charge(self.zone_temp_c_dist_initial.read().clone())?;

        // If detailed results are to be output, save the results from the current timestep
        if let Some(detailed_results) = self.detailed_results.as_ref() {
            let service_results = self.service_results.read();
            let services_called: IndexMap<ArcStr, &HeatBatteryResult> = service_results
                .iter()
                .map(|result| (result.service_name.clone(), result))
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
                        temp_inlet: None,
                        time_running: 0.,
                        energy_delivered_hb: 0.,
                        energy_delivered_backup: 0.,
                        energy_delivered_total: 0.,
                        energy_charged_during_service: 0.,
                        hb_zone_temperatures: initial_temps.clone(),
                        current_hb_power: 0.0.into(),
                    });
                }
            }

            // Add auxiliary results at the end
            let n_units = self.n_units as f64;
            let battery_losses = self.battery_losses.load(Ordering::SeqCst) * n_units;
            let total_charge = self.energy_charged_total.load(Ordering::SeqCst) * n_units;

            detailed_results.write().push(HeatBatteryTimestepResult {
                results: ordered_service_results,
                summary: HeatBatteryTimestepSummary {
                    energy_aux: energy_aux * n_units,
                    battery_losses,
                    temps_after_losses: zone_temp_c_after_losses,
                    total_charge,
                    end_of_timestep_charge: end_of_ts_charge * n_units,
                    hb_after_only_charge_zone_temp: zone_temp_c_after_charging,
                    state_of_charge: soc_end_of_timestep,
                },
            });
        }

        // Accumulate the energy delivered to each service this timestep into the running total
        // before __service_results is reset, so the per-service delivery survives for the
        // post-simulation charging-energy apportionment.
        let mut energy_delivered_by_service = self.energy_delivered_by_service.write();
        for result in self.service_results.read().iter() {
            *energy_delivered_by_service
                .entry(result.service_name.clone())
                .or_default() += result.energy_delivered_total;
        }
        self.total_time_running_current_timestep
            .store(Default::default(), Ordering::SeqCst);
        self.pump_running_time_current_timestep
            .store(Default::default(), Ordering::SeqCst);
        self.time_running_direct_current_timestep
            .store(Default::default(), Ordering::SeqCst);
        *self.service_results.write() = Default::default();
        self.energy_charged_total
            .store(Default::default(), Ordering::SeqCst);
        self.energy_charged_electric
            .store(Default::default(), Ordering::SeqCst);

        Ok(())
    }

    /// Output detailed results of heat battery calculation.
    pub(crate) fn output_detailed_results(
        &self,
        _hot_water_energy_output: &IndexMap<ArcStr, Vec<ResultParamValue>>,
        _hot_water_source_name_for_heat_battery_service: &IndexMap<ArcStr, ArcStr>,
    ) -> Result<(ResultsPerTimestep, ResultsAnnual), OutputDetailedResultsNotEnabledError> {
        let detailed_results = self
            .detailed_results
            .as_ref()
            .ok_or(OutputDetailedResultsNotEnabledError)?;
        let detailed_results = detailed_results.read();

        let mut results_per_timestep: ResultsPerTimestep =
            [("auxiliary".into(), Default::default())].into();

        // Report auxiliary parameters (not specific to a service)
        for (parameter, param_unit, _) in AUX_PARAMETERS {
            // Check if the parameter is a list to handle individual elements.
            if ["Temps_after_losses", "hb_after_only_charge_zone_temp"].contains(&parameter) {
                let mut labels: Option<Vec<ArcStr>> = Default::default();
                for service_results in detailed_results.iter() {
                    let summary = &service_results.summary;
                    let param_values = match parameter {
                        "Temps_after_losses" => &summary.temps_after_losses,
                        "hb_after_only_charge_zone_temp" => &summary.hb_after_only_charge_zone_temp,
                        _ => unreachable!(),
                    };
                    // Create labels for each element in the list.
                    // Determine the number of elements in the list for this parameter.
                    if labels.is_none() {
                        labels = Some(
                            param_values
                                .iter()
                                .enumerate()
                                .map(|(i, _)| format!("{parameter}{i}").into())
                                .collect(),
                        );
                    }
                    // Append each element to the respective label.
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
                for service_results in detailed_results.iter() {
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
        for service_name in self.energy_supply_connections.keys() {
            let service_name: ArcStr = service_name.into();
            let mut results_for_service: ResultPerTimestep = Default::default();

            // Look up each required parameter. Match the service by name rather than by
            // position (see _timestep_record_for_service), because services are recorded in
            // the order they are called, which need not match registration order.
            for (parameter, param_unit, _) in OUTPUT_PARAMETERS {
                // Create labels for each element in the list using the last timestep.
                if parameter == "hb_zone_temperatures" {
                    let last_record = timestep_record_for_service(
                        detailed_results
                            .last()
                            .expect("detailed results must contain a timestep"),
                        &service_name,
                    );
                    let labels = (0..last_record.hb_zone_temperatures.len())
                        .map(|i| format!("{parameter}{i}").into())
                        .collect_vec();
                    for label in labels {
                        results_for_service.insert(
                            (label, param_unit.map(|unit| unit.to_string().into())),
                            Vec::new(),
                        );
                    }
                } else {
                    results_for_service.insert(
                        (
                            parameter.into(),
                            param_unit.map(|unit| unit.to_string().into()),
                        ),
                        Vec::new(),
                    );
                }

                // Look up value of required parameter in each timestep
                for service_results in detailed_results.iter() {
                    let record = timestep_record_for_service(service_results, &service_name);
                    if parameter == "hb_zone_temperatures" {
                        let labels: Vec<ArcStr> = (0..record.hb_zone_temperatures.len())
                            .map(|i| format!("{parameter}{i}").into())
                            .collect_vec();
                        for (label, result) in labels.iter().zip(&record.hb_zone_temperatures) {
                            let key = (
                                label.clone(),
                                param_unit.map(|unit| unit.to_string().into()),
                            );
                            results_for_service
                                .get_mut(&key)
                                .expect("zone count must match the final timestep")
                                .push(result.into());
                        }
                    } else {
                        let result = record.param(parameter);
                        results_for_service
                            .get_mut(&(
                                parameter.into(),
                                param_unit.map(|unit| unit.to_string().into()),
                            ))
                            .expect("output parameter series must be initialized")
                            .push(result);
                    }
                }
            }

            results_per_timestep.insert(service_name.clone(), results_for_service);
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
            let service_name: ArcStr = service_name.into();
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
        }
        std::result::Result::Ok((results_per_timestep, results_annual))
    }

    fn target_charge(&self, simtime: SimulationTimeIteration) -> anyhow::Result<f64> {
        if let Some(charge) = &self.charge_control {
            charge.target_charge(simtime, None)
        } else {
            unreachable!()
        }
    }

    /// Get the current setpoints for a charging source, converting if needed.
    ///
    ///        When schedule_unit is "temperature", converts the raw temperature
    ///        setpoints from the RangeTimeControl to SOC values using the battery's
    ///        energy calculation. When schedule_unit is "soc", returns unchanged.
    ///
    ///        Args:
    ///            source: Charging source with control and schedule_unit.
    ///
    ///        Returns:
    ///            Tuple of (setpnt_lower, setpnt_upper) as SOC values (0–1).
    ///
    fn resolve_setpoints(
        &self,
        source: &HeatBatteryChargingSource,
        simtime: &SimulationTimeIteration,
    ) -> anyhow::Result<(Option<f64>, Option<f64>)> {
        let (mut setpnt_lower, mut setpnt_upper) =
            source.control.setpnt_range_time_control(simtime);
        if source.schedule_unit == ScheduleUnit::Temperature {
            if let Some(lower) = setpnt_lower {
                setpnt_lower = Some(self.temp_to_soc(lower)?);
            }
            if let Some(upper) = setpnt_upper {
                setpnt_upper = Some(self.temp_to_soc(upper)?);
            }
        }
        // Reject invalid combination: lower set but upper not set.
        // Follows StorageTank._retrieve_setpnt() convention.
        if setpnt_upper.is_none() && setpnt_lower.is_some() {
            bail!("schedule_lower must be None when schedule_upper is None");
        }

        Ok((setpnt_lower, setpnt_upper))
    }
}

#[derive(Debug, Error)]
#[error("Tried to call output_detailed_results when option to collect detailed results was not selected")]
pub(crate) struct OutputDetailedResultsNotEnabledError;

type ResultPerTimestep = IndexMap<(ArcStr, Option<ArcStr>), Vec<ResultParamValue>>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::common::{MockWaterSupply, VaryingTempWaterSupply};
    use crate::core::controls::time_control::{
        ChargeControl, Control, ScheduleOrControl, SetpointTimeControl,
    };
    use crate::core::energy_supply::energy_supply::{
        EnergySupply, EnergySupplyBuilder, EnergySupplyConnection,
    };
    use crate::core::water_heat_demand::misc::WaterEventResultType;
    use crate::external_conditions::{DaylightSavingsConfig, ExternalConditions};
    use crate::input::{
        ControlLogicType, ExternalSensor, FuelType, HeatBattery as HeatBatteryInput,
        HeatSourceWetDetails, PcmBatteryChargingConfiguration,
    };
    use crate::simulation_time::{SimulationTime, SimulationTimeIteration, SimulationTimeIterator};
    use approx::{assert_relative_eq, assert_relative_ne};
    use indexmap::indexmap;
    use itertools::Itertools;
    use parking_lot::RwLock;
    use rstest::*;
    use serde_json::json;
    use std::sync::atomic::Ordering;
    use std::sync::Arc;

    const SERVICE_NAME: &str = "TestService";

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 2., 1.)
    }

    #[fixture]
    fn simulation_time_iterator(simulation_time: SimulationTime) -> SimulationTimeIterator {
        simulation_time.iter()
    }

    #[fixture]
    fn simulation_time_iteration(
        simulation_time_iterator: SimulationTimeIterator,
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
            vec![200., 220.].into_iter().map(Into::into).collect(),
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
        simulation_time_iterator: SimulationTimeIterator,
    ) -> Arc<ChargeControl> {
        create_control_with_value(
            false,
            external_conditions,
            external_sensor,
            simulation_time_iterator,
        )
    }

    #[fixture]
    fn battery_control_on(
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
        simulation_time_iterator: SimulationTimeIterator,
    ) -> Arc<ChargeControl> {
        create_control_with_value(
            true,
            external_conditions,
            external_sensor,
            simulation_time_iterator,
        )
    }
    fn create_heat_battery_details(
        simultaneuous_charging_and_discharging: Option<bool>,
        max_temperature: Option<f64>,
    ) -> HeatSourceWetDetails {
        HeatSourceWetDetails::HeatBattery {
            battery: HeatBatteryInput::Pcm {
                energy_supply: "mains elec".into(),
                electricity_circ_pump: 0.06,
                electricity_standby: 0.0244,
                max_rated_losses: 0.1,
                number_of_units: 1,
                charging_config: PcmBatteryChargingConfiguration::ChargeControl {
                    control_charge: "hb_charge_control".into(),
                    rated_charge_power: 20.0,
                },
                simultaneous_charging_and_discharging: simultaneuous_charging_and_discharging
                    .unwrap_or(false),
                heat_storage_kj_per_k_above_phase_transition: 381.5,
                heat_storage_kj_per_k_below_phase_transition: 305.2,
                heat_storage_kj_per_k_during_phase_transition: 12317.,
                phase_transition_temperature_upper: 59.,
                phase_transition_temperature_lower: 57.,
                max_temperature: max_temperature.unwrap_or(80.),
                temp_init: 80.,
                velocity_in_hex_tube_at_1_l_per_min_m_per_s: 0.035,
                inlet_diameter_mm: 6.5,
                a: 174.33952,
                b: -931.565,
                flow_rate_l_per_min: 10.,
            },
        }
    }

    fn create_control_with_value(
        boolean: bool,
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
        simulation_time_iterator: SimulationTimeIterator,
    ) -> Arc<ChargeControl> {
        ChargeControl::new(
            ControlLogicType::Manual,
            ScheduleOrControl::Schedule(vec![boolean, boolean]),
            &simulation_time_iterator,
            0,
            1.,
            vec![Some(0.2)],
            None,
            None,
            Some(external_conditions.into()),
            Some(external_sensor),
            None,
        )
        .unwrap()
        .into()
    }

    fn temp_air_int_callback() -> TempInternalAirFn {
        Arc::new(|| 22.)
    }

    #[derive(Default)]
    struct HeatBatteryPcmTestOverrides {
        n_layers: Option<usize>,
        temp_min_useful: Option<f64>,
        energy_supply_end_usr_name: Option<String>,
        simultaneuous_charging_and_discharging: Option<bool>,
    }

    fn create_heat_battery(
        control: Arc<ChargeControl>,
        _output_detailed_results: Option<bool>,
    ) -> Arc<RwLock<HeatBatteryPcm>> {
        create_default_heat_battery_with_overides(
            control,
            _output_detailed_results,
            HeatBatteryPcmTestOverrides::default(),
        )
    }
    /// Creates a default heat battery with the given overrides.
    /// # Arguments
    ///
    /// * `control` - The charge control for the heat battery.
    /// * `_output_detailed_results` - Whether to output detailed results.
    /// * `overrides` - The overrides for the heat battery configuration.
    ///     - use Default::default() for any fields not explicitly overridden.
    fn create_default_heat_battery_with_overides(
        control: Arc<ChargeControl>,
        _output_detailed_results: Option<bool>,
        overrides: HeatBatteryPcmTestOverrides,
    ) -> Arc<RwLock<HeatBatteryPcm>> {
        let n_layers = overrides.n_layers.unwrap_or(8);
        let hb_time_step = 20.;
        let energy_supply_end_usr_name = overrides
            .energy_supply_end_usr_name
            .unwrap_or_else(|| "WaterHeating".to_string());
        let simulation_time = simulation_time();
        let energy_supply: Arc<RwLock<EnergySupply>> = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time.total_steps()).build(),
        ));
        let external_conditions = external_conditions(simulation_time);
        let energy_supply_connection: EnergySupplyConnection =
            EnergySupply::connection(energy_supply.clone(), &energy_supply_end_usr_name).unwrap();
        let heat_battery_details =
            create_heat_battery_details(overrides.simultaneuous_charging_and_discharging, None);

        let heat_battery = Arc::new(RwLock::new(
            HeatBatteryPcm::new(
                &heat_battery_details,
                energy_supply.into(),
                energy_supply_connection,
                simulation_time.iter(),
                external_conditions.into(),
                overrides.temp_min_useful,
                temp_air_int_callback(),
                Some(control),
                None,
                Some(n_layers),
                Some(hb_time_step),
                None,
                None,
                None,
                None,
            )
            .unwrap(),
        ));

        HeatBatteryPcm::create_service_connection(heat_battery.clone(), SERVICE_NAME).unwrap();

        heat_battery
    }

    fn create_setpoint_time_control(schedule: Vec<Option<f64>>) -> Control {
        Control::SetpointTime(
            SetpointTimeControl::new(schedule, 0, 1., Default::default(), Default::default(), 1.)
                .into(),
        )
    }

    fn get_service_names_from_results(heat_battery: Arc<RwLock<HeatBatteryPcm>>) -> Vec<ArcStr> {
        heat_battery
            .read()
            .service_results
            .read()
            .iter()
            .map(|result| result.service_name.clone())
            .collect_vec()
    }

    #[rstest]
    fn test_per_zone_conductance_is_total_ua_divided_by_n_layers(
        battery_control_on: Arc<ChargeControl>,
    ) {
        let n_layers = 4; // Example value for the number of layers
        let heat_battery_lock = create_default_heat_battery_with_overides(
            battery_control_on.clone(),
            None,
            HeatBatteryPcmTestOverrides {
                n_layers: Some(n_layers),
                ..Default::default()
            },
        );
        let heat_battery = heat_battery_lock.read();

        let pcm_temp_c = 70.;
        let inlet_temp_c = 40.;

        // get_zone_properties now derives the mass flow internally from the
        // volumetric flow, so drive the comparison from the same l/min the model
        // uses and convert it the same way for the manual outlet calculations.
        let flow_rate_l_per_min = heat_battery.flow_rate_l_per_min;
        let flow_rate_kg_per_s = flow_rate_l_per_min / 60.0 * WATER.density();
        let reynold_number_at_1_l_per_min = 1234.0;

        // Whole-exchanger conductance from the correlation, and its per-zone share.
        let ua_total_kw_per_k = HeatBatteryPcm::calculate_heat_transfer_kw_per_k(
            heat_battery.a,
            heat_battery.b,
            flow_rate_l_per_min,
            reynold_number_at_1_l_per_min,
        );
        let outlet_divided = HeatBatteryPcm::calculate_outlet_temp_c(
            ua_total_kw_per_k / n_layers as f64,
            pcm_temp_c,
            inlet_temp_c,
            flow_rate_kg_per_s,
        );
        let outlet_undivided = HeatBatteryPcm::calculate_outlet_temp_c(
            ua_total_kw_per_k,
            pcm_temp_c,
            inlet_temp_c,
            flow_rate_kg_per_s,
        );

        // The split must actually change the answer for n_layers > 1, otherwise
        // the assertion below would hold even if the division were removed.
        assert_relative_ne!(outlet_divided, outlet_undivided, epsilon = 1e-3);
        let mut zone_temp_dist = vec![pcm_temp_c; n_layers];
        let (_, zone_index, zone_temp_c_start, outlet_temp_c) = heat_battery.get_zone_properties(
            0,
            &HeatBatteryPcmOperationMode::Normal,
            &mut zone_temp_dist,
            inlet_temp_c,
            inlet_temp_c,
            0.,
            reynold_number_at_1_l_per_min,
            60.,
            flow_rate_l_per_min,
            None,
            None,
        );
        assert_eq!(zone_index, 0);
        assert_eq!(zone_temp_c_start, pcm_temp_c);
        assert_relative_eq!(outlet_temp_c, outlet_divided, epsilon = 1e-9);
    }

    fn discharge(
        n_layers: usize,
        inlet_temp_c: f64,
        pcm_temp_c: f64,
        battery_control: Arc<ChargeControl>,
    ) -> (f64, f64) {
        let hb_lock = create_default_heat_battery_with_overides(
            battery_control.clone(),
            None,
            HeatBatteryPcmTestOverrides {
                n_layers: Some(n_layers),
                ..Default::default()
            },
        );
        let mut zone_temp_c_dist = vec![pcm_temp_c; n_layers];
        let hb = hb_lock.read();
        // Physical flow and Reynolds number, derived the same way the model
        // does, so the only thing differing between passes is n_layers.
        let viscosity =
            HeatBatteryPcm::calculate_water_kinematic_viscosity_m2_per_s(inlet_temp_c, pcm_temp_c);
        let reynold_number_at_1_l_per_min = HeatBatteryPcm::calculate_reynold_number_at_1_l_per_min(
            viscosity,
            hb.velocity_in_hex_tube,
            hb.capillary_diameter_m,
        );
        let (outlet_temp_c, _, energy_transfer_delivered) = hb
            .process_heat_battery_zones(
                inlet_temp_c,
                &mut zone_temp_c_dist,
                60.,
                reynold_number_at_1_l_per_min,
                hb.flow_rate_l_per_min,
                None,
                Some(HeatBatteryPcmOperationMode::Normal),
                0.,
                None,
                None,
            )
            .unwrap();
        (outlet_temp_c, energy_transfer_delivered)
    }
    #[rstest]
    fn test_discharge_outputs_invariant_to_n_layers(battery_control_on: Arc<ChargeControl>) {
        // Discharge outlet temperature and energy must not depend on n_layers.
        //
        // n_layers is a numerical discretisation parameter, so a single discharge
        // pass through a uniformly charged battery must give the same physical
        // result whatever value it takes. Splitting the whole-exchanger UA per
        // layer keeps the total NTU invariant; applying the full UA per layer
        // inflates NTU by the layer count and makes the result diverge.
        let inlet_temp_c = 40.0;
        let pcm_temp_c = 80.0;

        let (outlet_4, energy_4) =
            discharge(4, inlet_temp_c, pcm_temp_c, battery_control_on.clone());
        let (outlet_8, energy_8) =
            discharge(8, inlet_temp_c, pcm_temp_c, battery_control_on.clone());
        let (outlet_16, energy_16) =
            discharge(16, inlet_temp_c, pcm_temp_c, battery_control_on.clone());

        // Tolerances sit well below the multi-degree / hundreds-of-kJ spread the
        // unfixed full-UA-per-layer code produces, but absorb the residual
        // arithmetic-mean discretisation error between layer counts.
        assert_relative_eq!(outlet_8, outlet_4, epsilon = 0.1);
        assert_relative_eq!(outlet_16, outlet_8, epsilon = 0.1);
        assert_relative_eq!(energy_8, energy_4, epsilon = 10.0);
        assert_relative_eq!(energy_16, energy_8, epsilon = 10.0);
    }

    #[rstest]
    fn test_demand_energy_inlet_hotter_than_battery_delivers_zero(
        battery_control_on: Arc<ChargeControl>,
        simulation_time_iteration: SimulationTimeIteration,
    ) {
        // demand_energy delivers nothing when the return-feed inlet is hotter than the core.
        //
        // If the inlet flow is hotter than every battery zone, the sub-timestep would
        // transfer heat into the battery rather than out of it; demand_energy breaks before
        // counting any negative contribution, so no energy is delivered.
        let hb_lock = create_heat_battery(battery_control_on, None);
        HeatBatteryPcm::create_service_connection(hb_lock.clone(), "hot_inlet_service").unwrap();
        let hb = hb_lock.read();
        let delivered = hb
            .demand_energy(
                "hot_inlet_service".into(),
                HeatingServiceType::DomesticHotWaterRegular,
                5.0,
                Some(90.0),
                Some(95.0),
                true,
                None,
                None,
                simulation_time_iteration,
            )
            .unwrap();
        assert_eq!(delivered, 0.0);
    }
    // in Python this test is called test_service_is_on_with_control
    #[rstest]
    fn test_service_is_on_when_service_control_is_on(
        simulation_time_iteration: SimulationTimeIteration,
        battery_control_off: Arc<ChargeControl>,
    ) {
        // Test when controlvent is provided and returns True
        let service_control_on: Control =
            create_setpoint_time_control(vec![Some(21.0), Some(21.0)]);

        let heat_battery = create_heat_battery(battery_control_off, None);

        let heat_battery_service = HeatBatteryPcmServiceSpace::new(
            heat_battery.clone(),
            SERVICE_NAME.into(),
            service_control_on,
        );

        assert!(heat_battery_service.is_on(simulation_time_iteration));

        let service_control_off: Control = create_setpoint_time_control(vec![None, None]);

        let heat_battery_service: HeatBatteryPcmServiceSpace =
            HeatBatteryPcmServiceSpace::new(heat_battery, SERVICE_NAME.into(), service_control_off);

        assert!(!heat_battery_service.is_on(simulation_time_iteration));
    }

    #[fixture]
    fn heat_battery_service_water_direct(
        battery_control_off: Arc<ChargeControl>,
    ) -> HeatBatteryPcmServiceWaterDirect {
        let heat_battery = create_heat_battery(battery_control_off, None);
        let mock_cold_feed = WaterSupply::Mock(MockWaterSupply::new(10.));
        let service_name = "WaterHeating".into();

        HeatBatteryPcmServiceWaterDirect::new(heat_battery, service_name, 60., mock_cold_feed)
    }

    #[rstest]
    fn test_get_cold_water_source_for_water_direct(
        heat_battery_service_water_direct: HeatBatteryPcmServiceWaterDirect,
    ) {
        let _expected = &WaterSupply::Mock(MockWaterSupply::new(10.));

        let actual = heat_battery_service_water_direct.get_cold_water_source();

        if let WaterSupply::Mock(mock) = actual {
            assert_eq!(mock, &MockWaterSupply::new(10.));
        } else {
            panic!("Expected a MockWaterSupply");
        }
    }

    #[rstest]
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
        battery_control: Arc<ChargeControl>,
        force_is_on: Option<bool>,
    ) -> HeatBatteryPcmServiceWaterRegular {
        let heat_battery = create_heat_battery(battery_control, None);
        let simulation_time = SimulationTime::new(0., 8., 1.);
        let simulation_time_iterator = simulation_time.iter();

        let shedule_lower_defaults = vec![
            Some(52.),
            None,
            None,
            None,
            Some(52.),
            Some(52.),
            Some(52.),
            Some(52.),
        ];
        let shedule_upper_defaults = vec![
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
        ];
        let schedule_lower_always_on = vec![Some(52.); 8];
        let schedule_upper_always_on = vec![Some(55.); 8];
        let schedule_always_off = vec![None; 8];

        let (schedule_lower, schedule_upper) = match force_is_on {
            Some(true) => (
                schedule_lower_always_on.clone(),
                schedule_upper_always_on.clone(),
            ),
            Some(false) => (schedule_always_off.clone(), schedule_always_off.clone()),
            None => (
                shedule_lower_defaults.clone(),
                shedule_upper_defaults.clone(),
            ),
        };
        let range_time_control = Arc::new(
            RangeTimeControl::new(
                ScheduleOrControl::Schedule(schedule_lower.clone()),
                ScheduleOrControl::Schedule(schedule_upper.clone()),
                simulation_time_iterator,
                0,
                1.,
                None,
            )
            .unwrap(),
        );

        let mock_cold_feed = WaterSupply::Mock(MockWaterSupply::new(10.));

        HeatBatteryPcmServiceWaterRegular::new(
            heat_battery,
            SERVICE_NAME.into(),
            mock_cold_feed,
            range_time_control,
        )
    }

    // test_service_is_on_without_control
    #[rstest]
    fn test_service_with_no_service_control_is_always_on_for_water_regular(
        simulation_time_iteration: SimulationTimeIteration,
        battery_control_off: Arc<ChargeControl>,
    ) {
        let heat_battery_service =
            create_service_water_regular_with_controls(battery_control_off, Some(true));

        assert!(heat_battery_service.is_on(simulation_time_iteration));
    }

    // this test accounts for python test test_setpnt_with_range_time_control as we only deal with RTC in Rust
    #[rstest]
    fn test_setpnt_for_water_regular(
        simulation_time_iterator: SimulationTimeIterator,
        battery_control_off: Arc<ChargeControl>,
    ) {
        let service = create_service_water_regular_with_controls(battery_control_off, None);

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
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
        simulation_time_iterator: SimulationTimeIterator,
        battery_control_on: Arc<ChargeControl>,
    ) {
        let energy_demand = 10.;
        let temp_flow = 55.;
        let temp_return = 40.;

        let range_time_control = Arc::new(
            RangeTimeControl::new(
                ScheduleOrControl::Schedule(vec![None]),
                ScheduleOrControl::Schedule(vec![None]),
                simulation_time_iterator.clone(),
                0,
                1.,
                None,
            )
            .unwrap(),
        );

        let heat_battery = create_heat_battery(battery_control_on, None);
        let mock_cold_feed = WaterSupply::Mock(MockWaterSupply::new(10.));
        let heat_battery_service: HeatBatteryPcmServiceWaterRegular =
            HeatBatteryPcmServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                mock_cold_feed,
                range_time_control,
            );

        let result = heat_battery_service
            .demand_energy(
                energy_demand,
                Some(temp_flow),
                Some(temp_return),
                None,
                simulation_time_iteration,
                false,
            )
            .unwrap();

        assert_eq!(result, 0.);
    }

    // skipped test_control_off_bypassed_by_ignore_standard_ctrl due to mocking and the minimal complexity of the change
    // skipped test_energy_output_max_when_service_control_on_for_water_regular as the python just tests the
    // method passes the correct arguments to the next call

    #[rstest]
    fn test_energy_output_max_service_off_for_water_regular(
        // In Python this is test_energy_output_max_service_off
        simulation_time_iteration: SimulationTimeIteration,
        battery_control_off: Arc<ChargeControl>,
    ) {
        let heat_battery_service =
            create_service_water_regular_with_controls(battery_control_off, Some(false));

        let temp_flow = 50.0;
        let temp_return = 40.0;
        let result = heat_battery_service
            // added false to match signature not yet ported for 1.0.0a9
            .energy_output_max(temp_flow, temp_return, simulation_time_iteration, false)
            .unwrap();

        assert_relative_eq!(result, 0., epsilon = 1e-7);
    }

    #[rstest]
    fn test_temp_setpnt_for_space(
        simulation_time_iteration: SimulationTimeIteration,
        battery_control_off: Arc<ChargeControl>,
    ) {
        let first_scheduled_temp = Some(21.);
        let ctrl: Control = create_setpoint_time_control(vec![first_scheduled_temp]);
        let heat_battery = create_heat_battery(battery_control_off, None);
        let heat_battery_space =
            HeatBatteryPcmServiceSpace::new(heat_battery, SERVICE_NAME.into(), ctrl);

        assert_eq!(
            heat_battery_space.temp_setpnt(simulation_time_iteration),
            first_scheduled_temp
        );
    }

    #[rstest]
    fn test_in_required_period_for_space(
        simulation_time_iteration: SimulationTimeIteration,
        battery_control_off: Arc<ChargeControl>,
    ) {
        let ctrl: Control = create_setpoint_time_control(vec![Some(21.)]);
        let heat_battery = create_heat_battery(battery_control_off, None);
        let heat_battery_space =
            HeatBatteryPcmServiceSpace::new(heat_battery, SERVICE_NAME.into(), ctrl);

        assert_eq!(
            heat_battery_space.in_required_period(simulation_time_iteration),
            Some(true)
        );
    }

    #[rstest]
    fn test_demand_energy_service_off_for_space(
        simulation_time_iteration: SimulationTimeIteration,

        battery_control_off: Arc<ChargeControl>,
    ) {
        let energy_demand = 10.;
        let temp_return = 40.;
        let temp_flow = 1.;
        let time_start = 0.2;
        let ctrl: Control = create_setpoint_time_control(vec![None]);
        let heat_battery = create_heat_battery(battery_control_off, None);
        let heat_battery_space =
            HeatBatteryPcmServiceSpace::new(heat_battery, SERVICE_NAME.into(), ctrl);
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
        battery_control_on: Arc<ChargeControl>,
        simulation_time_iteration: SimulationTimeIteration,
    ) {
        let temp_output = 70.;
        let temp_return = 40.;
        let heat_battery = create_heat_battery(battery_control_on, None);
        let service_control_off: Control = create_setpoint_time_control(vec![None]);

        let heat_battery_service: HeatBatteryPcmServiceSpace =
            HeatBatteryPcmServiceSpace::new(heat_battery, SERVICE_NAME.into(), service_control_off);

        let result = heat_battery_service
            .energy_output_max(temp_output, temp_return, None, simulation_time_iteration)
            .unwrap();

        assert_relative_eq!(result, 0.);
    }

    #[rstest]
    fn test_create_service_connection(battery_control_on: Arc<ChargeControl>) {
        let heat_battery = create_heat_battery(battery_control_on, None);
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
    fn test_create_service_hot_water_direct(battery_control_on: Arc<ChargeControl>) {
        let heat_battery = create_heat_battery(battery_control_on, None);
        let mock_cold_feed = WaterSupply::Mock(MockWaterSupply::new(10.));
        let service = HeatBatteryPcm::create_service_hot_water_direct(
            heat_battery.clone(),
            "new_service",
            60.,
            mock_cold_feed,
        )
        .unwrap();

        let actual = service.get_cold_water_source();

        if let WaterSupply::Mock(mock) = actual {
            assert_eq!(mock, &MockWaterSupply::new(10.));
        } else {
            panic!("Expected a MockWaterSupply");
        }

        assert!(heat_battery
            .read()
            .energy_supply_connections
            .contains_key("new_service"));
    }

    #[rstest]
    fn test_create_service_space_heating(
        simulation_time_iteration: SimulationTimeIteration,
        battery_control_off: Arc<ChargeControl>,
    ) {
        let control = create_setpoint_time_control(vec![Some(21.0)]);
        let heat_battery = create_heat_battery(battery_control_off, None);
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
        simulation_time_iteration: SimulationTimeIteration,
        battery_control_off: Arc<ChargeControl>,
        battery_control_on: Arc<ChargeControl>,
    ) {
        // electric charge should be 0 when battery control is off
        let heat_battery = create_heat_battery(battery_control_off, None);

        assert_relative_eq!(
            heat_battery
                .read()
                .electric_charge(simulation_time_iteration),
            0.0
        );

        // electric charge should be calculated when battery control is on
        let heat_battery = create_heat_battery(battery_control_on, None);
        assert_relative_eq!(
            heat_battery
                .read()
                .electric_charge(simulation_time_iteration),
            20.0
        );
    }

    #[rstest]
    fn test_first_call(battery_control_on: Arc<ChargeControl>, simulation_time: SimulationTime) {
        let heat_battery = create_heat_battery(battery_control_on, None);
        for t_it in simulation_time.iter() {
            heat_battery.read().first_call();

            assert!(!heat_battery.read().flag_first_call.load(Ordering::SeqCst));

            heat_battery.read().timestep_end(t_it).unwrap();
        }
    }

    #[rstest]
    fn test_demand_energy(simulation_time: SimulationTime, battery_control_on: Arc<ChargeControl>) {
        let heat_battery = create_heat_battery(battery_control_on, None);

        let expected_zone_temp_c_dist = [
            vec![
                79.95404545080883,
                79.95765589934315,
                79.96098269068034,
                79.96404811052943,
                79.96687269370845,
                79.96947536170435,
                79.97187354942538,
                79.97408332199502,
            ], // First timestep
            vec![
                78.9034689243783,
                78.91053402045067,
                78.91704436349183,
                78.92304351282773,
                78.92857160749791,
                78.93366563481702,
                78.9383596778487,
                78.94268514344776,
            ], // Second timestep
        ];

        let service_name = "new_service";
        HeatBatteryPcm::create_service_connection(heat_battery.clone(), service_name).unwrap();

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let demand_energy_actual = heat_battery
                .clone()
                .read()
                .demand_energy(
                    service_name,
                    HeatingServiceType::DomesticHotWaterRegular,
                    5.,
                    Some(40.),
                    Some(52.5),
                    true,
                    Some(1.), // Python passes True positionally as time_start.
                    None,
                    t_it,
                )
                .unwrap();

            assert_relative_eq!(
                demand_energy_actual,
                [0.0037217831482159714, 0.0036250231710519888][t_idx]
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

            heat_battery.read().timestep_end(t_it).unwrap();
        }
    }

    fn create_heat_battery_pcm(
        external_sensor: ExternalSensor,
        simulation_time_iterator: SimulationTimeIterator,
        external_conditions: ExternalConditions,
    ) -> Arc<RwLock<HeatBatteryPcm>> {
        let control = Arc::new(
            ChargeControl::new(
                ControlLogicType::Manual,
                ScheduleOrControl::Schedule(vec![false]),
                &simulation_time_iterator,
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
        let heat_battery = create_heat_battery(control, None);
        HeatBatteryPcm::create_service_connection(heat_battery.clone(), "new_service").unwrap();

        heat_battery
    }

    // skipping python's test_demand_energy_simultaneous_charging_and_discharging due to mocking

    #[rstest]
    fn test_demand_energy_simultaneous_no_temp_output(
        battery_control_off: Arc<ChargeControl>,
        simulation_time_iteration: SimulationTimeIteration,
    ) {
        let heat_battery = create_heat_battery(battery_control_off, None);
        assert_relative_eq!(
            heat_battery
                .read()
                .demand_energy(
                    SERVICE_NAME,
                    HeatingServiceType::DomesticHotWaterRegular,
                    0.08,
                    Some(40.),
                    None,
                    true,
                    None,
                    None,
                    simulation_time_iteration
                )
                .unwrap(),
            0.0849030119531927
        );

        assert_relative_eq!(
            heat_battery
                .read()
                .demand_energy(
                    SERVICE_NAME,
                    HeatingServiceType::DomesticHotWaterRegular,
                    0.06,
                    Some(40.),
                    None,
                    true,
                    None,
                    None,
                    simulation_time_iteration
                )
                .unwrap(),
            0.06506455710564082,
        );

        // With no accumulated losses, nothing is creditable to internal gains
        // self.assertEqual(self.heatbattery.get_battery_losses(), 0.0)
        //
        // Only the recoverable fraction (f_sto_m = 0.75, BS EN 15316-5:2017
        // Table B.3) of accumulated standing losses is returned for crediting
        // to internal gains; the remainder escapes the dwelling.
        // self.heatbattery._HeatBatteryPCM__battery_losses = 1.0  # type: ignore[reportAttributeAccessIssue]
        // 1.0 kWh raw loss * 1 unit (number_of_units) * 0.75 recovery factor = 0.75 kWh
        // self.assertEqual(self.heatbattery.get_battery_losses(), 0.75)
        // The accumulator is reset once read
        assert_eq!(heat_battery.read().get_battery_losses(), 0.);
    }

    #[rstest]
    fn test_demand_energy_other(
        external_sensor: ExternalSensor,
        simulation_time_iterator: SimulationTimeIterator,
        external_conditions: ExternalConditions,
    ) {
        let heat_battery = create_heat_battery_pcm(
            external_sensor.clone(),
            simulation_time_iterator.clone(),
            external_conditions.clone(),
        );
        let simtime = simulation_time_iterator.current_iteration();
        assert_relative_eq!(
            heat_battery
                .read()
                .demand_energy(
                    "new_service",
                    HeatingServiceType::DomesticHotWaterRegular,
                    0.08,
                    Some(40.),
                    Some(40.),
                    true,
                    None,
                    None,
                    simtime
                )
                .unwrap(),
            0.0849030119531927
        );

        let heat_battery = create_heat_battery_pcm(
            external_sensor.clone(),
            simulation_time_iterator.clone(),
            external_conditions.clone(),
        );
        heat_battery.write().hb_time_step = 119.;

        assert_relative_eq!(
            heat_battery
                .read()
                .demand_energy(
                    "new_service",
                    HeatingServiceType::DomesticHotWaterRegular,
                    0.08,
                    Some(40.),
                    Some(40.),
                    true,
                    None,
                    None,
                    simtime
                )
                .unwrap(),
            0.08691254317474525
        );

        let heat_battery = create_heat_battery_pcm(
            external_sensor.clone(),
            simulation_time_iterator.clone(),
            external_conditions.clone(),
        );
        heat_battery.write().hb_time_step = 20.;

        assert_relative_eq!(
            heat_battery
                .read()
                .demand_energy(
                    "new_service",
                    HeatingServiceType::DomesticHotWaterRegular,
                    0.08,
                    Some(40.),
                    Some(80.),
                    true,
                    None,
                    None,
                    simtime
                )
                .unwrap(),
            0.0849030119531927
        );

        let heat_battery = create_heat_battery_pcm(
            external_sensor,
            simulation_time_iterator,
            external_conditions,
        );

        assert_relative_eq!(
            heat_battery
                .read()
                .demand_energy(
                    "new_service",
                    HeatingServiceType::DomesticHotWaterRegular,
                    0.08,
                    Some(40.),
                    Some(79.),
                    true,
                    None,
                    None,
                    simtime
                )
                .unwrap(),
            0.0849030119531927
        );
    }

    #[rstest]
    fn test_dhw_service_demand_hot_water(
        battery_control_off: Arc<ChargeControl>,
        simulation_time_iteration: SimulationTimeIteration,
    ) {
        let heat_battery = create_heat_battery(battery_control_off, None);
        let mock_cold_feed = WaterSupply::Mock(MockWaterSupply::new(10.));
        let service = HeatBatteryPcm::create_service_hot_water_direct(
            heat_battery,
            "dhw_complex",
            65., // High setpoint
            mock_cold_feed,
        )
        .unwrap();

        let actual = service.get_cold_water_source();

        if let WaterSupply::Mock(mock) = actual {
            assert_eq!(mock, &MockWaterSupply::new(10.));
        } else {
            panic!("Expected a MockWaterSupply");
        }

        // Test with usage events
        let usage_events = vec![
            WaterEventResult {
                event_result_type: WaterEventResultType::Other,
                temperature_warm: 40.0,
                volume_warm: 50.0,
                volume_hot: 8.0,
                event_duration: 0.0,
            },
            WaterEventResult {
                event_result_type: WaterEventResultType::Other,
                temperature_warm: 35.0,
                volume_warm: 0.0,
                volume_hot: 0.0,
                event_duration: 0.0,
            },
        ];

        let energy = service
            .demand_hot_water(Some(usage_events), simulation_time_iteration)
            .unwrap();

        assert_eq!(energy, 0.46581279161760664);

        // Test with no usage events
        let energy_no_usage = service
            .demand_hot_water(None, simulation_time_iteration)
            .unwrap();

        assert_eq!(energy_no_usage, 0.);
    }

    // Skipping Python's test_calc_auxiliary_energy due to mocking (only assertion uses assert_called_once_with)
    // Skipping Python's test_calc_auxiliary_energy_space_heating due to mocking (only assertion uses assert_called_once_with)

    /// Check heat battery auxiliary energy includes standby power when no services are called
    #[rstest]
    fn test_calc_auxiliary_energy_no_services(
        simulation_time_iterator: SimulationTimeIterator,
        battery_control_on: Arc<ChargeControl>,
    ) {
        // Don't create any services
        let heat_battery = create_heat_battery(battery_control_on, None);

        let result = heat_battery
            .read()
            .calc_auxiliary_energy(1.0, 0.5, simulation_time_iterator.current_index())
            .unwrap();

        // Should only have standby power (no pump power since no services were called)
        let expected_energy_aux = heat_battery.read().power_standby * 0.5;

        assert_relative_eq!(result, expected_energy_aux);
    }

    /// Check that direct DHW service doesn't contribute to pump running time
    #[rstest]
    fn test_calc_auxiliary_energy_direct_dhw_no_pump_contribution(
        simulation_time_iterator: SimulationTimeIterator,
        battery_control_on: Arc<ChargeControl>,
    ) {
        let heat_battery = create_heat_battery(battery_control_on, None);
        let mock_cold_feed = WaterSupply::Mock(MockWaterSupply::new(10.));
        // Create only a direct hot water service
        HeatBatteryPcm::create_service_hot_water_direct(
            heat_battery.clone(),
            "dhw_direct",
            60.,
            mock_cold_feed,
        )
        .unwrap();

        // Simulate demand that sets total_time_running but should not affect pump time
        heat_battery
            .write()
            .total_time_running_current_timestep
            .store(0.5, Ordering::SeqCst);
        heat_battery
            .write()
            .pump_running_time_current_timestep
            .store(0., Ordering::SeqCst);

        let result = heat_battery
            .read()
            .calc_auxiliary_energy(1.0, 0.5, simulation_time_iterator.current_index())
            .unwrap();

        // Only standby power, no pump power since pump time is 0
        let expected_energy_aux = heat_battery.read().power_standby * 0.5;

        assert_relative_eq!(result, expected_energy_aux);
    }

    #[rstest]
    fn test_timestep_end(
        external_sensor: ExternalSensor,
        external_conditions: ExternalConditions,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        // not using the fixture here
        // because we need to set different charge_levels
        let battery_control_on: Arc<ChargeControl> = ChargeControl::new(
            ControlLogicType::Manual,
            ScheduleOrControl::Schedule(vec![true, true, true]),
            &simulation_time_iterator,
            0,
            1.,
            [1.0, 1.5].into_iter().map(Into::into).collect(),
            None,
            None,
            Some(external_conditions.into()),
            Some(external_sensor),
            None,
        )
        .unwrap()
        .into();

        let heat_battery = create_heat_battery(battery_control_on, None);
        let service_name = "new_timestep_end_service";
        HeatBatteryPcm::create_service_connection(heat_battery.clone(), service_name).unwrap();

        let simtime = simulation_time_iterator.current_iteration();
        heat_battery
            .read()
            .demand_energy(
                service_name,
                HeatingServiceType::DomesticHotWaterRegular,
                5.0,
                Some(40.),
                Some(55.),
                true,
                None,
                None,
                simtime,
            )
            .unwrap();

        assert_relative_eq!(
            heat_battery
                .read()
                .total_time_running_current_timestep
                .load(Ordering::SeqCst),
            0.6256433826572602
        );

        let service_names_in_results = get_service_names_from_results(heat_battery.clone());

        assert!(service_names_in_results.contains(&service_name.into()));

        heat_battery.read().timestep_end(simtime).unwrap();

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
        simulation_time_iterator: SimulationTimeIterator,
        simulation_time: SimulationTime,
    ) {
        // not using the fixture here
        // because we need to set different charge_levels
        let battery_control_on: Arc<ChargeControl> = ChargeControl::new(
            ControlLogicType::Manual,
            ScheduleOrControl::Schedule(vec![true, true]),
            &simulation_time_iterator,
            0,
            1.,
            [0.95, 0.99].into_iter().map(Into::into).collect(), // these values change the result
            None,
            None,
            Some(external_conditions.clone().into()),
            Some(external_sensor.clone()),
            None,
        )
        .unwrap()
        .into();

        let heat_battery = create_heat_battery(battery_control_on.clone(), None);

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            // Inlet to the heat exchanger is the return-feed temperature, not the
            // required flow temperature, so the available output reflects the true
            // temperature difference driving heat transfer. The battery does not
            // charge while discharging (simultaneous_charging_and_discharging is
            // False here), matching demand_energy, so the ceiling reflects only the
            // stored energy.
            assert_relative_eq!(
                heat_battery
                    .read()
                    .energy_output_max(40., 30., None, t_it)
                    .unwrap(),
                [9.850915334607903, 9.78071474696652][t_idx],
                epsilon = 5e-8
            );

            heat_battery.read().timestep_end(t_it).unwrap();
        }

        let heat_battery = create_heat_battery(battery_control_on, None);

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            // The required flow temperature does not gate delivery: even when 90 °C
            // cannot be reached, the battery still transfers heat to the 40 °C return
            // feed, and energy_output_max counts that full positive delivery just as
            // demand_energy does.
            assert_relative_eq!(
                heat_battery
                    .read()
                    .energy_output_max(90., 40., None, t_it)
                    .unwrap(),
                [7.390489238946976, 7.326317707345256][t_idx],
                epsilon = 5e-8
            );

            heat_battery.read().timestep_end(t_it).unwrap();
        }
    }

    #[rstest]
    fn test_get_zone_properties_losses(battery_control_off: Arc<ChargeControl>) {
        // Test that the losses model scales a zone's loss by its temperature difference.
        //
        // The zone's share of the rated loss energy (Q_max_kJ split across the
        // eight zones) is scaled by the ratio of the zone-to-surroundings
        // temperature difference to the difference at which the rated loss was
        // characterised (max_temperature 80 °C minus the 20 °C reference ambient,
        // i.e. 60 K).
        let heat_battery = create_heat_battery(battery_control_off, None);
        let (energy_transf, _, _, _) = heat_battery.read().get_zone_properties(
            0,
            &HeatBatteryPcmOperationMode::Losses,
            &[42., 57., 58., 58., 59., 59., 60., 61.],
            40.,
            40.,
            5.,
            414.,
            20.,
            10.,
            None,
            None,
        );
        // 5 / 8 zones * (42 - 40) K / 60 K = 0.0208333...
        assert_eq!(energy_transf, 0.020833333333333332);
    }

    #[rstest]
    fn test_get_zone_properties_no_energy_transf(battery_control_off: Arc<ChargeControl>) {
        // Test that the losses model clamps to zero when a zone is at or below the surroundings
        let heat_battery = create_heat_battery(battery_control_off, None);
        let (energy_transf, _, _, _) = heat_battery.read().get_zone_properties(
            0,
            &HeatBatteryPcmOperationMode::Losses,
            &[42., 57., 58., 58., 59., 59., 60., 61.],
            45.,
            45.,
            5.,
            414.,
            20.,
            10.,
            None,
            None,
        );

        assert_eq!(energy_transf, 0.);
    }

    // skipping python's test_get_zone_properties_invalid_mode as mode can't be invalid in rust

    #[rstest]
    fn test_get_zone_properties_losses_at_rated_condition(battery_control_off: Arc<ChargeControl>) {
        // Test that a zone at the rated condition loses its full equal share of the rated loss.
        //
        // When a zone sits at max_temperature (80 °C) and the surroundings are at
        // the rated reference ambient (20 °C), the temperature-difference ratio is
        // one, so the zone loses an equal split of the rated loss energy
        // (Q_max_kJ / number of zones) — the same as the unscaled model.
        let heat_battery = create_heat_battery(battery_control_off, None);

        let (energy_transf, _, _, _) = heat_battery.read().get_zone_properties(
            0,
            &HeatBatteryPcmOperationMode::Losses,
            &[80., 80., 80., 80., 80., 80., 80., 80.],
            20.,
            20.,
            5.,
            414.,
            20.,
            12.9847456,
            None,
            None,
        );

        // 5 / 8 zones * (80 - 20) K / 60 K = 0.625
        assert_relative_eq!(energy_transf, 0.625);
    }

    #[rstest]
    fn test_calculate_zone_energy_required(battery_control_off: Arc<ChargeControl>) {
        let heat_battery = create_heat_battery(battery_control_off, None);

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
    fn test_process_zone_simultaneous_charging(battery_control_off: Arc<ChargeControl>) {
        let heat_battery = create_heat_battery(battery_control_off, None);

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

        let (q_max_kj, energy_charged, energy_transf) = heat_battery
            .read()
            .process_zone_simultaneous_charging(58., 120., -1000., 1900., 0.);

        assert_relative_eq!(q_max_kj, 0.);
        assert_relative_eq!(energy_charged, 0.2777777777777778);
        assert_relative_eq!(energy_transf, 900.);
    }

    #[rstest]
    fn test_process_zone_simultaneous_charging_partial_recovery(
        battery_control_off: Arc<ChargeControl>,
    ) {
        // Test simultaneous charging when charging power is insufficient to recover all withdrawn energy.

        // energy_transf=1900 (withdrawn), q_max_kj=-500 (charging power)
        // -q_max_kj=500 < 1900=energy_transf -> partial recovery branch
        let heat_battery = create_heat_battery(battery_control_off, None);
        let (q_max_kj, energy_charged, energy_transf) = heat_battery
            .read()
            .process_zone_simultaneous_charging(58., 120., -500., 1900., 0.);

        // Charging partially recovers: energy_transf += q_max_kj = 1900 + (-500) = 1400
        assert_relative_eq!(energy_transf, 1400.);
        // All charging power used: q_max_kj = 0
        assert_relative_eq!(q_max_kj, 0.);
        // energy_charged = -q_max_kj / kJ_per_kWh = 500 / 3600
        assert_relative_eq!(energy_charged, 0.13888888888888888);
    }

    #[rstest]
    fn test_process_zone_simultaneous_charging_fully_charged_partial_recovery(
        battery_control_off: Arc<ChargeControl>,
    ) {
        // Test simultaneous charging when the zone starts fully charged and the charging power is insufficient to recover all the energy withdrawn by the inlet water.

        // zone_temp_c_start=58 >= temp_charge_target=20 -> zone starts fully charged,
        // routing through the fully-charged branch (distinct from the below-target
        // branch exercised by test_process_zone_simultaneous_charging_partial_recovery).
        // energy_transf=1900 (withdrawn), q_max_kj=-500 (charging power):
        // -q_max_kj=500 < 1900=energy_transf -> partial recovery branch.
        let heat_battery = create_heat_battery(battery_control_off, None);
        let (q_max_kj, energy_charged, energy_transf) = heat_battery
            .read()
            .process_zone_simultaneous_charging(58., 20., -500., 1900., 0.);

        // Charging partially recovers: energy_transf += q_max_kj = 1900 + (-500) = 1400
        assert_relative_eq!(energy_transf, 1400.);
        // All charging power used: q_max_kj = 0
        assert_relative_eq!(q_max_kj, 0.);
        // energy_charged = -q_max_kj / kJ_per_kWh = 500 / 3600
        assert_relative_eq!(energy_charged, 0.13888888888888888);
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
        battery_control_off: Arc<ChargeControl>,
        #[case] zone_temp_c_start: f64,
        #[case] energy_transf: f64,
        #[case] expected: f64,
    ) {
        let heat_battery = create_heat_battery(battery_control_off, None);
        let result = heat_battery
            .read()
            .calculate_new_zone_temperature(zone_temp_c_start, energy_transf);

        assert_relative_eq!(result, expected);
    }

    #[rstest]
    fn test_charge_battery_hydronic(battery_control_off: Arc<ChargeControl>) {
        let heat_battery = create_heat_battery(battery_control_off, None);
        let charge = |inlet_temp_c| {
            let heat_battery_guard = heat_battery.read();
            let (energy, _) = heat_battery_guard
                .charge_battery_hydronic(
                    inlet_temp_c,
                    1.0,
                    10.0,
                    1.0,
                    f64::INFINITY,
                    174.33952,
                    -931.565,
                    0.035,
                    6.5 / 1000.0,
                )
                .unwrap();
            energy
        };

        assert_relative_eq!(charge(70.0), 0.0);
        assert_relative_eq!(charge(80.0), 65.72932929968711 / 3600.0);
        assert_relative_eq!(charge(90.0), 3751.226341395129 / 3600.0);
    }

    #[rstest]
    fn test_charge_battery_hydronic_energy_limit(battery_control_off: Arc<ChargeControl>) {
        let heat_battery = create_heat_battery(battery_control_off, None);
        let (energy, zone_temps) = heat_battery
            .read()
            .charge_battery_hydronic(
                90.0,
                1.0,
                10.0,
                1.0,
                0.5,
                174.33952,
                -931.565,
                0.035,
                6.5 / 1000.0,
            )
            .unwrap();

        assert_relative_eq!(energy, 0.5);
        let expected_zone_temps = [
            85.87275090774999,
            85.5034629005933,
            85.15503953386835,
            84.82654437538658,
            84.51705558542992,
            84.22566956638788,
            83.95150395049762,
            83.69370000554323,
        ];
        assert_eq!(zone_temps.len(), expected_zone_temps.len());
        for (actual, expected) in zone_temps.iter().zip(expected_zone_temps) {
            assert_relative_eq!(*actual, expected, epsilon = 1e-10);
        }
    }

    #[rstest]
    fn test_heat_battery_charge_battery_hydronic_edge_cases(
        battery_control_off: Arc<ChargeControl>,
    ) {
        // Test hydronic charging with boundary temperature conditions.
        let heat_battery = create_heat_battery(battery_control_off, None);

        // Very high inlet temperature -> significant charging
        let (energy, zones) = heat_battery
            .read()
            .charge_battery_hydronic(
                95.0,
                1.0,
                10.0,
                1.0,
                f64::INFINITY,
                174.33952,
                -931.565,
                0.035,
                6.5 / 1000.0,
            )
            .unwrap();
        assert_relative_eq!(energy, 1.5619144113112322);
        assert_eq!(zones.len(), 8);

        // Inlet at max storage temp - no energy transfer
        let (energy, _) = heat_battery
            .read()
            .charge_battery_hydronic(
                80.0,
                1.0,
                10.0,
                1.0,
                f64::INFINITY,
                174.33952,
                -931.565,
                0.035,
                6.5 / 1000.0,
            )
            .unwrap();
        assert_relative_eq!(energy, 0.0);
    }

    #[rstest]
    fn test_get_temp_hot_water(
        battery_control_off: Arc<ChargeControl>,
        simulation_time_iteration: SimulationTimeIteration,
    ) {
        let heat_battery = create_heat_battery(battery_control_off, None);
        let simtime = simulation_time_iteration;

        assert_relative_eq!(
            heat_battery
                .read()
                .get_temp_hot_water(50., 20., 80., simtime)
                .unwrap(),
            78.79649236439744
        );
        assert_relative_eq!(
            heat_battery
                .read()
                .get_temp_hot_water(50., 10., 80., simtime)
                .unwrap(),
            76.27701993478482
        );
        assert_relative_eq!(
            heat_battery
                .read()
                .get_temp_hot_water(40., 10., 80., simtime)
                .unwrap(),
            74.84690041657751,
        );
        assert_relative_eq!(
            heat_battery
                .read()
                .get_temp_hot_water(60., 1., 65., simtime)
                .unwrap(),
            65.
        );
    }

    // skipping python's test_energy_output_max_negative as unable to replicate patch object

    #[fixture]
    fn heat_battery_no_service_connection(
        simulation_time_iterator: SimulationTimeIterator,
        battery_control_on: Arc<ChargeControl>,
        external_conditions: ExternalConditions,
    ) -> Arc<RwLock<HeatBatteryPcm>> {
        let heat_battery_details: &HeatSourceWetDetails = &HeatSourceWetDetails::HeatBattery {
            battery: HeatBatteryInput::Pcm {
                energy_supply: "mains elec".into(),
                electricity_circ_pump: 0.06,
                electricity_standby: 0.0244,
                max_rated_losses: 0.1,
                number_of_units: 1,
                charging_config: PcmBatteryChargingConfiguration::ChargeControl {
                    control_charge: "hb_charge_control".into(),
                    rated_charge_power: 20.0,
                },
                simultaneous_charging_and_discharging: false,
                heat_storage_kj_per_k_above_phase_transition: 381.5,
                heat_storage_kj_per_k_below_phase_transition: 305.2,
                heat_storage_kj_per_k_during_phase_transition: 12317.,
                phase_transition_temperature_upper: 59.,
                phase_transition_temperature_lower: 57.,
                max_temperature: 80.,
                temp_init: 80.,
                velocity_in_hex_tube_at_1_l_per_min_m_per_s: 0.035,
                inlet_diameter_mm: 6.5,
                a: 174.33952,
                b: -931.565,
                flow_rate_l_per_min: 10.,
            },
        };

        let energy_supply: Arc<RwLock<EnergySupply>> = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time_iterator.total_steps())
                .build(),
        ));

        let energy_supply_connection: EnergySupplyConnection =
            EnergySupply::connection(energy_supply.clone(), "WaterHeating").unwrap();

        Arc::new(RwLock::new(
            HeatBatteryPcm::new(
                &heat_battery_details,
                energy_supply.into(),
                energy_supply_connection,
                simulation_time_iterator,
                external_conditions.into(),
                None,
                temp_air_int_callback(),
                battery_control_on.into(),
                None,
                Some(8),
                Some(20.),
                None,
                None,
                Some(true),
                None,
            )
            .unwrap(),
        ))
    }

    #[rstest]
    fn test_output_detailed_results_water_regular(
        simulation_time: SimulationTime,
        heat_battery_no_service_connection: Arc<RwLock<HeatBatteryPcm>>,
    ) {
        let heat_battery = heat_battery_no_service_connection;
        let mock_cold_feed = WaterSupply::Mock(MockWaterSupply::new(10.));
        let service_name = "new_service";

        let range_time_control = RangeTimeControl::new(
            ScheduleOrControl::Schedule(vec![]),
            ScheduleOrControl::Schedule(vec![]),
            simulation_time.iter(),
            0,
            1.,
            None,
        )
        .unwrap();

        HeatBatteryPcm::create_service_hot_water_regular(
            heat_battery.clone(),
            service_name,
            mock_cold_feed,
            range_time_control.into(),
        )
        .unwrap();

        let expected_results_per_timestep: ResultsPerTimestep = indexmap! {
            "auxiliary".into() => indexmap! {
                ("energy_aux".into(), Some("kWh".into())) => vec![0.06.into(), 0.06.into()],
                ("battery_losses".into(), Some("kWh".into())) => vec![0.06830005082775674.into(), 0.0411152773366954.into()],
                ("Temps_after_losses0".into(), Some("degC".into())) => vec![52.09705202087361.into(), 39.56966446334649.into()],
                ("Temps_after_losses1".into(), Some("degC".into())) => vec![57.06703887771538.into(), 39.616910213222845.into()],
                ("Temps_after_losses2".into(), Some("degC".into())) => vec![57.254259772056514.into(), 39.69421041287997.into()],
                ("Temps_after_losses3".into(), Some("degC".into())) => vec![57.42527650005248.into(), 39.81901697860602.into()],
                ("Temps_after_losses4".into(), Some("degC".into())) => vec![57.581399089081124.into(), 40.018095866602636.into()],
                ("Temps_after_losses5".into(), Some("degC".into())) => vec![57.72384825236103.into(), 40.3327159263897.into()],
                ("Temps_after_losses6".into(), Some("degC".into())) => vec![57.8537148110231.into(), 40.82647693395094.into()],
                ("Temps_after_losses7".into(), Some("degC".into())) => vec![57.97201857678866.into(), 41.59642394900029.into()],
                ("total_charge".into(), Some("kWh".into())) => vec![0.0.into(); 2],
                ("end_of_timestep_charge".into(), Some("kWh".into())) => vec![0.0.into(); 2],
                ("hb_after_only_charge_zone_temp0".into(), Some("degC".into())) => vec![52.09705202087361.into(), 39.56966446334649.into()],
                ("hb_after_only_charge_zone_temp1".into(), Some("degC".into())) => vec![57.06703887771538.into(), 39.616910213222845.into()],
                ("hb_after_only_charge_zone_temp2".into(), Some("degC".into())) => vec![57.254259772056514.into(), 39.69421041287997.into()],
                ("hb_after_only_charge_zone_temp3".into(), Some("degC".into())) => vec![57.42527650005248.into(), 39.81901697860602.into()],
                ("hb_after_only_charge_zone_temp4".into(), Some("degC".into())) => vec![57.581399089081124.into(), 40.018095866602636.into()],
                ("hb_after_only_charge_zone_temp5".into(), Some("degC".into())) => vec![57.72384825236103.into(), 40.3327159263897.into()],
                ("hb_after_only_charge_zone_temp6".into(), Some("degC".into())) => vec![57.8537148110231.into(), 40.82647693395094.into()],
                ("hb_after_only_charge_zone_temp7".into(), Some("degC".into())) => vec![57.97201857678866.into(), 41.59642394900029.into()],
                ("state_of_charge".into(), Some("ratio".into())) => vec![0.463198529730531.into(), 0.24507891561838296.into()],
            },
            "new_service".into() => indexmap! {
                ("service_name".into(), None) => vec![ResultParamValue::String(arcstr::literal!("new_service")); 2],
                ("service_type".into(), None) => vec![ResultParamValue::String(HeatingServiceType::DomesticHotWaterRegular.to_string().into()); 2],
                ("service_on".into(), None) => vec![ResultParamValue::Boolean(true); 2],
                ("energy_output_required".into(), Some("kWh".into())) => vec![100.0.into(); 2],
                ("temp_output".into(), Some("degC".into())) => vec![48.83495256437514.into(), 40.41092501564586.into()],
                ("temp_inlet".into(), Some("degC".into())) => vec![40.0.into(); 2],
                ("time_running".into(), Some("secs".into())) => vec![3600.0.into(), 3600.0.into()],
                ("energy_delivered_HB".into(), Some("kWh".into())) => vec![7.393523697804974.into(), 2.9908624775073656.into()],
                ("energy_delivered_backup".into(), Some("kWh".into())) => vec![0.0.into(); 2],
                ("energy_delivered_total".into(), Some("kWh".into())) => vec![7.393523697804974.into(), 2.9908624775073656.into()],
                ("energy_charged_during_service".into(), Some("kWh".into())) => vec![0.0.into(); 2],
                ("hb_zone_temperatures0".into(), Some("degC".into())) => vec![
                    52.82092338492856.into(),
                    40.04231816247777.into(),
                ],
                ("hb_zone_temperatures1".into(), Some("degC".into())) => vec![57.0870536801901.into(), 40.09051135386234.into()],
                ("hb_zone_temperatures2".into(), Some("degC".into())) => vec![57.27436582019495.into(), 40.169361691213126.into()],
                ("hb_zone_temperatures3".into(), Some("degC".into())) => vec![57.44546589644598.into(), 40.2966710624016.into()],
                ("hb_zone_temperatures4".into(), Some("degC".into())) => vec![57.60166457478777.into(), 40.499742174087984.into()],
                ("hb_zone_temperatures5".into(), Some("degC".into())) => vec![57.74418316337672.into(), 40.820671459672916.into()],
                ("hb_zone_temperatures6".into(), Some("degC".into())) => vec![57.87411301497617.into(), 41.32433409171734.into()],
                ("hb_zone_temperatures7".into(), Some("degC".into())) => vec![57.99247443833206.into(), 42.10972122070484.into()],
                ("current_hb_power".into(), Some("kW".into())) => vec![7.393523697804974.into(), 2.9908624775073656.into()],
            },
        };

        let expected_results_annual: ResultsAnnual = indexmap! {
            "Overall".into() => indexmap! {
                ("energy_output_required".into(), Some("kWh".into())) => 200.0.into(),
                ("time_running".into(), Some("secs".into())) => 7200.0.into(),
                ("energy_delivered_HB".into(), Some("kWh".into())) => 10.38438617531234.into(),
                ("energy_delivered_backup".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_total".into(), Some("kWh".into())) => 10.38438617531234.into(),
                ("energy_charged_during_service".into(), Some("kWh".into())) => 0.0.into(),
            },
            "auxiliary".into() => indexmap! {
                ("energy_aux".into(), Some("kWh".into())) => 0.12.into(),
                ("battery_losses".into(), Some("kWh".into())) => 0.10941532816445214.into(),
                ("total_charge".into(), Some("kWh".into())) => 0.0.into(),
                ("end_of_timestep_charge".into(), Some("kWh".into())) => 0.0.into(),
            },
            "new_service".into() => indexmap! {
                ("energy_output_required".into(), Some("kWh".into())) => 200.0.into(),
                ("time_running".into(), Some("secs".into())) => 7200.0.into(),
                ("energy_delivered_HB".into(), Some("kWh".into())) => 10.38438617531234.into(),
                ("energy_delivered_backup".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_total".into(), Some("kWh".into())) => 10.38438617531234.into(),
                ("energy_charged_during_service".into(), Some("kWh".into())) => 0.0.into(),
            },
        };

        for t_it in simulation_time.iter() {
            heat_battery
                .read()
                .demand_energy(
                    service_name,
                    HeatingServiceType::DomesticHotWaterRegular,
                    100.,
                    Some(40.),
                    Some(55.),
                    true,
                    None,
                    Some(true),
                    t_it,
                )
                .unwrap();

            heat_battery.read().timestep_end(t_it).unwrap();
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
    fn test_output_detailed_results_space(
        simulation_time: SimulationTime,
        heat_battery_no_service_connection: Arc<RwLock<HeatBatteryPcm>>,
    ) {
        let heat_battery = heat_battery_no_service_connection;
        let service_name = "new_service";
        let control = create_setpoint_time_control(vec![]);

        HeatBatteryPcm::create_service_space_heating(heat_battery.clone(), service_name, control)
            .unwrap();

        let expected_results_per_timestep: ResultsPerTimestep = indexmap! {
            "auxiliary".into() => indexmap! {
                ("energy_aux".into(), Some("kWh".into())) => vec![0.06.into(), 0.06.into()],
                ("battery_losses".into(), Some("kWh".into())) => vec![0.06830005082775674.into(), 0.0411152773366954.into()],
                ("Temps_after_losses0".into(), Some("degC".into())) => vec![52.09705202087361.into(), 39.56966446334649.into()],
                ("Temps_after_losses1".into(), Some("degC".into())) => vec![57.06703887771538.into(), 39.616910213222845.into()],
                ("Temps_after_losses2".into(), Some("degC".into())) => vec![57.254259772056514.into(), 39.69421041287997.into()],
                ("Temps_after_losses3".into(), Some("degC".into())) => vec![57.42527650005248.into(), 39.81901697860602.into()],
                ("Temps_after_losses4".into(), Some("degC".into())) => vec![57.581399089081124.into(), 40.018095866602636.into()],
                ("Temps_after_losses5".into(), Some("degC".into())) => vec![57.72384825236103.into(), 40.3327159263897.into()],
                ("Temps_after_losses6".into(), Some("degC".into())) => vec![57.8537148110231.into(), 40.82647693395094.into()],
                ("Temps_after_losses7".into(), Some("degC".into())) => vec![57.97201857678866.into(), 41.59642394900029.into()],
                ("total_charge".into(), Some("kWh".into())) => vec![0.0.into(), 0.0.into()],
                ("end_of_timestep_charge".into(), Some("kWh".into())) => vec![0.0.into(), 0.0.into()],
                ("hb_after_only_charge_zone_temp0".into(), Some("degC".into())) => vec![52.09705202087361.into(), 39.56966446334649.into()],
                ("hb_after_only_charge_zone_temp1".into(), Some("degC".into())) => vec![57.06703887771538.into(), 39.616910213222845.into()],
                ("hb_after_only_charge_zone_temp2".into(), Some("degC".into())) => vec![57.254259772056514.into(), 39.69421041287997.into()],
                ("hb_after_only_charge_zone_temp3".into(), Some("degC".into())) => vec![57.42527650005248.into(), 39.81901697860602.into()],
                ("hb_after_only_charge_zone_temp4".into(), Some("degC".into())) => vec![57.581399089081124.into(), 40.018095866602636.into()],
                ("hb_after_only_charge_zone_temp5".into(), Some("degC".into())) => vec![57.72384825236103.into(), 40.3327159263897.into()],
                ("hb_after_only_charge_zone_temp6".into(), Some("degC".into())) => vec![57.8537148110231.into(), 40.82647693395094.into()],
                ("hb_after_only_charge_zone_temp7".into(), Some("degC".into())) => vec![57.97201857678866.into(), 41.59642394900029.into()],
                ("state_of_charge".into(), Some("ratio".into())) => vec![0.463198529730531.into(), 0.24507891561838296.into()],
            },
            "new_service".into() => indexmap! {
                ("service_name".into(), None) => vec![ResultParamValue::String("new_service".into()), ResultParamValue::String("new_service".into())],
                ("service_type".into(), None) => vec![ResultParamValue::String(HeatingServiceType::Space.to_string().into()), ResultParamValue::String(HeatingServiceType::Space.to_string().into())],
                ("service_on".into(), None) => vec![ResultParamValue::Boolean(true), ResultParamValue::Boolean(true)],
                ("energy_output_required".into(), Some("kWh".into())) => vec![100.0.into(), 100.0.into()],
                ("temp_output".into(), Some("degC".into())) => vec![48.83495256437514.into(), 40.41092501564586.into()],
                ("temp_inlet".into(), Some("degC".into())) => vec![40.0.into(), 40.0.into()],
                ("time_running".into(), Some("secs".into())) => vec![3600.0.into(), 3600.0.into()],
                ("energy_delivered_HB".into(), Some("kWh".into())) => vec![7.393523697804974.into(), 2.9908624775073656.into()],
                ("energy_delivered_backup".into(), Some("kWh".into())) => vec![0.0.into(), 0.0.into()],
                ("energy_delivered_total".into(), Some("kWh".into())) => vec![7.393523697804974.into(), 2.9908624775073656.into()],
                ("energy_charged_during_service".into(), Some("kWh".into())) => vec![0.0.into(), 0.0.into()],
                ("hb_zone_temperatures0".into(), Some("degC".into())) => vec![52.82092338492856.into(), 40.04231816247777.into()],
                ("hb_zone_temperatures1".into(), Some("degC".into())) => vec![57.0870536801901.into(), 40.09051135386234.into()],
                ("hb_zone_temperatures2".into(), Some("degC".into())) => vec![57.27436582019495.into(), 40.169361691213126.into()],
                ("hb_zone_temperatures3".into(), Some("degC".into())) => vec![57.44546589644598.into(), 40.2966710624016.into()],
                ("hb_zone_temperatures4".into(), Some("degC".into())) => vec![57.60166457478777.into(), 40.499742174087984.into()],
                ("hb_zone_temperatures5".into(), Some("degC".into())) => vec![57.74418316337672.into(), 40.820671459672916.into()],
                ("hb_zone_temperatures6".into(), Some("degC".into())) => vec![57.87411301497617.into(), 41.32433409171734.into()],
                ("hb_zone_temperatures7".into(), Some("degC".into())) => vec![57.99247443833206.into(), 42.10972122070484.into()],
                ("current_hb_power".into(), Some("kW".into())) => vec![7.393523697804974.into(), 2.9908624775073656.into()],
            },
        };

        let expected_results_annual: ResultsAnnual = indexmap! {
            "Overall".into() => indexmap! {
                ("energy_output_required".into(), Some("kWh".into())) => 200.0.into(),
                ("time_running".into(), Some("secs".into())) => 7200.0.into(),
                ("energy_delivered_HB".into(), Some("kWh".into())) => 10.38438617531234.into(),
                ("energy_delivered_backup".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_total".into(), Some("kWh".into())) => 10.38438617531234.into(),
                ("energy_charged_during_service".into(), Some("kWh".into())) => 0.0.into(),
            },
            "auxiliary".into() => indexmap! {
                ("energy_aux".into(), Some("kWh".into())) => 0.12.into(),
                ("battery_losses".into(), Some("kWh".into())) => 0.10941532816445214.into(),
                ("total_charge".into(), Some("kWh".into())) => 0.0.into(),
                ("end_of_timestep_charge".into(), Some("kWh".into())) => 0.0.into(),
            },
            "new_service".into() => indexmap! {
                ("energy_output_required".into(), Some("kWh".into())) => 200.0.into(),
                ("time_running".into(), Some("secs".into())) => 7200.0.into(),
                ("energy_delivered_HB".into(), Some("kWh".into())) => 10.38438617531234.into(),
                ("energy_delivered_backup".into(), Some("kWh".into())) => 0.0.into(),
                ("energy_delivered_total".into(), Some("kWh".into())) => 10.38438617531234.into(),
                ("energy_charged_during_service".into(), Some("kWh".into())) => 0.0.into(),
            },
        };

        for t_it in simulation_time.iter() {
            heat_battery
                .read()
                .demand_energy(
                    service_name,
                    HeatingServiceType::Space,
                    100.,
                    Some(40.),
                    Some(55.),
                    true,
                    None,
                    Some(true),
                    t_it,
                )
                .unwrap();

            heat_battery.read().timestep_end(t_it).unwrap();
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
    fn test_output_detailed_results_none(battery_control_on: Arc<ChargeControl>) {
        // Test that calling output_detailed_results errors when output_detailed_results on heat_battery is false
        let heat_battery = create_heat_battery(battery_control_on, Some(false));

        assert!(heat_battery
            .read()
            .output_detailed_results(&indexmap! {}, &indexmap! {})
            .is_err());
    }

    #[rstest]
    fn test_demand_energy_low_temp_minimum_run_coverage(
        battery_control_off: Arc<ChargeControl>,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let heat_battery = create_heat_battery(battery_control_off, None);
        let simtime = simulation_time_iterator.current_iteration();
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
                Some(40.),
                Some(50.),
                true,
                Some(0.),
                Some(true),
                simtime,
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

        assert_relative_eq!(
            service_result.time_running,
            105.4442943033978,
            max_relative = 1e-7
        );
    }

    #[rstest]
    fn test_timestep_end_with_uncalled_services(
        heat_battery_no_service_connection: Arc<RwLock<HeatBatteryPcm>>,
        mut simulation_time_iterator: SimulationTimeIterator,
    ) {
        let heat_battery = heat_battery_no_service_connection;
        let simtime = simulation_time_iterator.current_iteration();

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
                Some(40.),
                Some(55.),
                true,
                Some(0.),
                Some(true),
                simtime,
            )
            .unwrap();

        heat_battery
            .read()
            .demand_energy(
                service3,
                HeatingServiceType::Space,
                3.,
                Some(35.),
                Some(50.),
                true,
                Some(0.),
                Some(true),
                simtime,
            )
            .unwrap();

        heat_battery.read().timestep_end(simtime).unwrap();

        simulation_time_iterator.next();
        let simtime = simulation_time_iterator.current_iteration();

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
                Some(38.),
                Some(52.),
                true,
                Some(0.),
                Some(true),
                simtime,
            )
            .unwrap();

        heat_battery.read().timestep_end(simtime).unwrap();

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
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        // Test timestep_end when no services are called but services are registered
        let heat_battery = heat_battery_no_service_connection;
        let simtime = simulation_time_iterator.current_iteration();

        // Create services but don't call them
        let service1 = "water_heating";
        let service2 = "space_heating";

        HeatBatteryPcm::create_service_connection(heat_battery.clone(), service1).unwrap();
        HeatBatteryPcm::create_service_connection(heat_battery.clone(), service2).unwrap();

        // Call timestep_end without calling any services
        heat_battery.read().timestep_end(simtime).unwrap();

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
        battery_control_off: Arc<ChargeControl>,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let heat_battery = create_heat_battery(battery_control_off, None);
        let service_name = "test_service";
        let mock_cold_feed = WaterSupply::Mock(MockWaterSupply::new(10.));

        let range_time_control = Arc::new(
            RangeTimeControl::new(
                ScheduleOrControl::Schedule(vec![]),
                ScheduleOrControl::Schedule(vec![]),
                simulation_time_iterator,
                0,
                1.,
                None,
            )
            .unwrap(),
        );

        let result = HeatBatteryPcm::create_service_hot_water_regular(
            heat_battery.clone(),
            service_name,
            mock_cold_feed.clone(),
            range_time_control.clone(),
        );

        assert!(result.is_ok());

        let result = HeatBatteryPcm::create_service_hot_water_regular(
            heat_battery,
            service_name,
            mock_cold_feed.clone(),
            range_time_control,
        );

        assert!(result.is_err())
    }

    // skipping python's test_heat_battery_edge_case_zero_timestep as function can't return None

    // skipping python's test_heat_battery_process_zone_edge_cases as function can't return None and does return 4 values

    // skipping python's test_heat_battery_charge_battery_hydraulic_edge_cases as function can't return None

    #[rstest]
    fn test_heat_battery_energy_output_max_boundary_conditions(
        battery_control_off: Arc<ChargeControl>,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        // Test energy_output_max with boundary temperature conditions
        let simtime = simulation_time_iterator.current_iteration();
        let heat_battery = create_heat_battery(battery_control_off, None);

        // Test with very low output temperature
        let result = heat_battery
            .read()
            .energy_output_max(10., 5., Some(0.), simtime)
            .unwrap();

        // The method returns energy based on zone temps, not necessarily 0
        assert_relative_eq!(result, 12.689231024029406, epsilon = 1e-7);

        //Test with temperature at threshold
        let result = heat_battery
            .read()
            .energy_output_max(45., 40., Some(0.), simtime)
            .unwrap();

        assert_relative_eq!(result, 7.390489238946976, epsilon = 1e-7);
    }

    // skipping python's test_heat_battery_service_cold_water_source_not_set as cold feed not optional

    #[rstest]
    fn test_heat_battery_zero_volume_zones(battery_control_off: Arc<ChargeControl>) {
        let heat_battery = create_heat_battery(battery_control_off, None);

        // Test with zero energy transfer
        let result = heat_battery.read().calculate_new_zone_temperature(50., 0.);

        assert_eq!(result, 50.);
    }

    #[rstest]
    fn test_heat_battery_all_zones_below_threshold(
        battery_control_off: Arc<ChargeControl>,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let simtime = simulation_time_iterator.current_iteration();
        // Request high output temperature that no zone can provide
        let heat_battery = create_heat_battery(battery_control_off, None);

        let result = heat_battery
            .read()
            .energy_output_max(80., 80., Some(0.), simtime)
            .unwrap();

        assert_relative_eq!(result, 0., epsilon = 1e-7);
    }

    /// Test DHW service with cold water temperature that varies with volume demanded
    #[rstest]
    fn test_demand_hot_water_with_varying_cold_temperatures(
        battery_control_off: Arc<ChargeControl>,
        simulation_time: SimulationTime,
    ) {
        let simtime = simulation_time.iter().current_iteration();
        let heat_battery = create_heat_battery(battery_control_off, None);

        // Set up cold feed to return different temperatures based on volume
        // Simulates drawing from a stratified tank or mixed sources
        fn varying_temp_by_volume(volume_needed: f64) -> Vec<(f64, f64)> {
            let volume = volume_needed;

            if volume <= 10. {
                // Small volume - warm water from top of tank
                vec![(15.0, volume)]
            } else if volume <= 30. {
                // Medium volume - mix of warm and cold
                let warm_portion = 10.;
                let cold_portion = volume - 10.;
                vec![(15.0, warm_portion), (8.0, cold_portion)]
            } else {
                // Large volume - mostly cold water
                vec![(15.0, 10.), (8.0, 20.), (5.0, volume - 30.)]
            }
        }

        let volumes_container: Arc<RwLock<Vec<f64>>> = Default::default();

        let mock_cold_feed =
            WaterSupply::VaryingTemp(VaryingTempWaterSupply::new(volumes_container.clone()));

        let service = HeatBatteryPcm::create_service_hot_water_direct(
            heat_battery.clone(),
            "dhw_varying_temp",
            65.0,
            mock_cold_feed,
        )
        .unwrap();

        // Test with different volume events
        let usage_events = vec![
            WaterEventResult {
                event_result_type: WaterEventResultType::Other, // the Python uses a nonexistent "HandWash" type here - this is the best equivalent
                temperature_warm: 35.0,
                volume_warm: 5.0,
                volume_hot: 5.0, // Small - should get 15°C
                event_duration: 0.,
            },
            WaterEventResult {
                event_result_type: WaterEventResultType::Shower,
                temperature_warm: 38.0,
                volume_warm: 40.0,
                volume_hot: 25.0, // Medium - should get mix (15°C and 8°C)
                event_duration: 0.,
            },
            WaterEventResult {
                event_result_type: WaterEventResultType::Bath,
                temperature_warm: 40.0,
                volume_warm: 80.0,
                volume_hot: 50.0, // Large - should get mix of all three temps
                event_duration: 0.,
            },
        ];

        // Execute
        let energy = service
            .demand_hot_water(usage_events.into(), simtime)
            .unwrap();

        // Varify draw_off_water was called with correct volumes
        let draw_volumes = volumes_container.read().clone();
        assert_eq!(draw_volumes.len(), 3);

        assert_eq!(draw_volumes[0], 5.0); // First event volume
        assert_eq!(draw_volumes[1], 25.0); // Second event volume
        assert_eq!(draw_volumes[2], 50.0); // Third event volume

        // Energy should be calculated based on varying temperatures
        assert_relative_eq!(energy, 5.0569591506976606, epsilon = 1e-7);

        // Test that different volumes give different inlet temperatures
        // Reset and test with single large volume
        volumes_container.write().clear();

        let single_large_event = vec![WaterEventResult {
            event_result_type: WaterEventResultType::Bath,
            temperature_warm: 40.0,
            volume_warm: 80.0,
            volume_hot: 40.0,
            event_duration: 0.,
        }];

        let energy_large = service
            .demand_hot_water(single_large_event.into(), simtime)
            .unwrap();

        // For 40L: 10L@15°C + 20L@8°C + 10L@5°C
        // Average = (150 + 160 + 50) / 40 = 9°C

        // Now test with equivalent volume but as small draws
        volumes_container.write().clear();

        let multiple_small_events = vec![
            WaterEventResult {
                event_result_type: WaterEventResultType::Other, // Python uses nonexistent type "Small" here
                temperature_warm: 40.0,
                volume_warm: 10.0,
                volume_hot: 8.0,
                event_duration: 0.,
            },
            WaterEventResult {
                event_result_type: WaterEventResultType::Other, // Python uses nonexistent type "Small" here
                temperature_warm: 40.0,
                volume_warm: 10.0,
                volume_hot: 8.0,
                event_duration: 0.,
            },
            WaterEventResult {
                event_result_type: WaterEventResultType::Other, // Python uses nonexistent type "Small" here
                temperature_warm: 40.0,
                volume_warm: 10.0,
                volume_hot: 8.0,
                event_duration: 0.,
            },
            WaterEventResult {
                event_result_type: WaterEventResultType::Other, // Python uses nonexistent type "Small" here
                temperature_warm: 40.0,
                volume_warm: 10.0,
                volume_hot: 8.0,
                event_duration: 0.,
            },
            WaterEventResult {
                event_result_type: WaterEventResultType::Other, // Python uses nonexistent type "Small" here
                temperature_warm: 40.0,
                volume_warm: 10.0,
                volume_hot: 8.0,
                event_duration: 0.,
            },
        ];

        let energy_small_batches = service
            .demand_hot_water(multiple_small_events.into(), simtime)
            .unwrap();

        // Small batches all get 15°C water, so should need less energy than large draw
        // (less heating required when inlet is 15°C vs 9°C average)
        assert!(energy_small_batches < energy_large);
        assert_relative_eq!(energy_small_batches, 1.4186641524169479, epsilon = 1e-7);
        assert_relative_eq!(energy_large, 2.2854792470827485, epsilon = 1e-7);
    }

    #[rstest]
    fn test_energy_delivered_by_service(
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
        simulation_time: SimulationTime,
    ) {
        // Two services share the battery, as hot water and space heating do in production; the
        // accessor must keep their delivered energy separate for the per-service apportionment.
        let simtime_iterator = simulation_time.iter();
        let control = ChargeControl::new(
            ControlLogicType::Manual,
            ScheduleOrControl::Schedule(vec![true; simulation_time.total_steps()]),
            &simtime_iterator,
            0,
            1.,
            vec![Some(1.0), Some(1.5)],
            None,
            None,
            Some(external_conditions.into()),
            Some(external_sensor),
            None,
        )
        .unwrap()
        .into();
        let heat_battery = create_default_heat_battery_with_overides(
            control,
            None,
            HeatBatteryPcmTestOverrides {
                temp_min_useful: Some(30.0),
                ..Default::default()
            },
        );
        // No timestep has ended yet, so nothing has been accumulated.

        assert!(heat_battery.read().energy_delivered_by_service().is_empty());

        HeatBatteryPcm::create_service_connection(heat_battery.clone(), "hw_service").unwrap();
        HeatBatteryPcm::create_service_connection(heat_battery.clone(), "sh_service").unwrap();
        // Demand 3.0 kWh for hot water and 2.0 kWh for space heating in each of two timesteps.
        for simtime in simulation_time.iter() {
            heat_battery
                .read()
                .demand_energy(
                    "hw_service",
                    HeatingServiceType::DomesticHotWaterRegular,
                    3.0,
                    Some(40.0),
                    Some(55.0),
                    true,
                    None,
                    Some(true),
                    simtime,
                )
                .unwrap();
            heat_battery
                .read()
                .demand_energy(
                    "sh_service",
                    HeatingServiceType::Space,
                    2.0,
                    Some(40.0),
                    Some(55.0),
                    true,
                    None,
                    Some(true),
                    simtime,
                )
                .unwrap();

            heat_battery.read().timestep_end(simtime).unwrap();
            assert!(heat_battery.read().service_results.read().is_empty());
        }

        let delivered = heat_battery.read().energy_delivered_by_service();
        let mut service_names: Vec<_> = delivered.keys().map(|name| name.as_str()).collect();
        service_names.sort_unstable();

        // Both demands are fully met each timestep, so the running totals are the per-service
        // demand summed over the two timesteps: 3.0 x 2 and 2.0 x 2, kept separate by service.
        assert_eq!(service_names, vec!["hw_service", "sh_service"]);
        assert_relative_eq!(delivered["hw_service"], 6.0, epsilon = 1e-7);
        assert_relative_eq!(delivered["sh_service"], 4.0, epsilon = 1e-7);
    }

    #[rstest]
    fn test_energy_output_max_simultaneous_charging(
        external_sensor: ExternalSensor,
        external_conditions: ExternalConditions,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        // With simultaneous charging, the ceiling includes charge replenished while discharging.
        //
        // A battery that charges while discharging is held warmer over the timestep, so
        // energy_output_max draws on the charging power in addition to stored energy.
        let battery_control: Arc<ChargeControl> = ChargeControl::new(
            ControlLogicType::Manual,
            ScheduleOrControl::Schedule(vec![true, true]),
            &simulation_time_iterator,
            0,
            1.,
            [1.5, 1.6].into_iter().map(Into::into).collect(),
            None,
            None,
            Some(external_conditions.into()),
            Some(external_sensor),
            None,
        )
        .unwrap()
        .into();
        let heat_battery_lock = create_default_heat_battery_with_overides(
            battery_control,
            None,
            HeatBatteryPcmTestOverrides {
                temp_min_useful: Some(30.),
                energy_supply_end_usr_name: Some("sim_charge".to_string()),
                simultaneuous_charging_and_discharging: Some(true),
                ..Default::default()
            },
        );
        let heat_battery = heat_battery_lock.read();
        assert_relative_eq!(
            heat_battery
                .energy_output_max(
                    40.0,
                    30.0,
                    None,
                    simulation_time_iterator.current_iteration()
                )
                .unwrap(),
            17.717635418087347,
            epsilon = 1e-7
        );
    }

    #[rstest]
    fn test_energy_output_max_matches_demand_energy_ceiling(
        external_sensor: ExternalSensor,
        external_conditions: ExternalConditions,
        simulation_time_iterator: SimulationTimeIterator,
        simulation_time_iteration: SimulationTimeIteration,
    ) {
        // energy_output_max equals the energy demand_energy delivers at unlimited demand.
        //
        // Both integrate the same heat transfer to the water, neither gates on the required
        // flow temperature, and both run at the same sub-timestep. So giving demand_energy
        // effectively unlimited demand reaches the energy_output_max ceiling. Under the
        // per-layer conductance the discharge no longer saturates within the timestep, so
        // the two integration schemes agree closely (well within 0.1%) but no longer to
        // floating-point precision.
        let battery_control: Arc<ChargeControl> = ChargeControl::new(
            ControlLogicType::Manual,
            ScheduleOrControl::Schedule(vec![true; simulation_time_iterator.total_steps()]),
            &simulation_time_iterator,
            0,
            1.,
            [1.5, 1.6].into_iter().map(Into::into).collect(),
            None,
            None,
            Some(external_conditions.into()),
            Some(external_sensor),
            None,
        )
        .unwrap()
        .into();

        // Build a fresh battery for each call so both see the identical initial state
        let energy_max = create_default_heat_battery_with_overides(
            battery_control.clone(),
            None,
            HeatBatteryPcmTestOverrides {
                temp_min_useful: Some(30.0),
                energy_supply_end_usr_name: Some("ceiling_max".to_string()),
                ..Default::default()
            },
        )
        .read()
        .energy_output_max(40.0, 30.0, None, simulation_time_iteration)
        .unwrap();
        let energy_delivered = create_default_heat_battery_with_overides(
            battery_control,
            None,
            HeatBatteryPcmTestOverrides {
                temp_min_useful: Some(30.0),
                energy_supply_end_usr_name: Some("ceiling_demand".to_string()),
                ..Default::default()
            },
        )
        .read()
        .demand_energy(
            "ceiling_demand",
            HeatingServiceType::DomesticHotWaterRegular,
            1.0e6,
            Some(30.0),
            Some(40.0),
            true,
            None,
            Some(false),
            simulation_time_iteration,
        )
        .unwrap();

        assert_relative_eq!(energy_max, 9.850915334607903, epsilon = 1e-7);
        assert_relative_eq!(energy_delivered, 9.852626311604713, epsilon = 1e-7);
        // The two integration schemes agree to well within 0.1%; they no longer
        // coincide to floating-point precision because the discharge does not
        // saturate the timestep under the per-layer conductance.
        assert!((energy_max - energy_delivered).abs() <= 0.01);
    }

    #[rstest]
    fn test_init_rejects_max_temperature_at_or_below_rated_ambient(
        external_conditions: ExternalConditions,
        simulation_time_iterator: SimulationTimeIterator,
        battery_control_on: Arc<ChargeControl>,
    ) {
        // Test that construction fails if max_temperature does not exceed the rated-loss ambient.
        //
        // The rated temperature difference (max_temperature minus the reference
        // ambient) is the denominator when scaling standing losses, so a
        // non-positive difference is a degenerate input that must be rejected.

        let energy_supply: Arc<RwLock<EnergySupply>> = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time_iterator.total_steps())
                .build(),
        ));

        let energy_supply_connection: EnergySupplyConnection =
            EnergySupply::connection(energy_supply.clone(), "WaterHeating").unwrap();
        let heat_data = create_heat_battery_details(None, Some(TEMP_AMBIENT_RATED_LOSSES_C));
        let heat_battery = HeatBatteryPcm::new(
            &heat_data,
            energy_supply.clone().into(),
            energy_supply_connection,
            simulation_time_iterator,
            external_conditions.into(),
            Some(30.),
            temp_air_int_callback(),
            Some(battery_control_on),
            None,
            Some(8),
            None,
            None,
            None,
            None,
            None,
        );
        assert!(heat_battery.is_err());
        assert_eq!(
            heat_battery.unwrap_err().to_string(),
            format!("Heat battery max_temperature must exceed the rated-loss reference ambient temperature ({TEMP_AMBIENT_RATED_LOSSES_C} °C).")
        );
    }

    // skipping python's test_demand_hot_water_zero_volume_continue due to mocking

    /// Tests for validate_no_schedule_overlap (Deviation 3 fix).
    /// Uses real RangeTimeControl objects to verify that overlapping active
    /// schedules are rejected and non-overlapping schedules are accepted.
    mod test_schedule_overlap_validation {
        use crate::core::controls::time_control::RangeTimeControl;

        use super::*;

        #[derive(Debug, Clone)]
        struct MockWaterSupply;

        //mock all as they don't matter
        impl WaterSupplyBehaviour for MockWaterSupply {
            fn draw_off_water(
                &self,
                _: f64,
                _: SimulationTimeIteration,
            ) -> anyhow::Result<Vec<(f64, f64)>> {
                Ok(vec![])
            }
            fn get_temp_cold_water(
                &self,
                _: f64,
                _: SimulationTimeIteration,
            ) -> anyhow::Result<Vec<(f64, f64)>> {
                Ok(vec![])
            }
            fn ultimate_cold_water_source(&self) -> Self {
                Self {}
            }
        }

        #[fixture]
        fn simtime() -> SimulationTime {
            SimulationTime::new(0., 4., 1.)
        }

        /// Create a RangeTimeControl with given schedule lists.
        fn make_range_time_control(
            schedule_lower: Vec<Option<f64>>,
            schedule_upper: Vec<Option<f64>>,
            simtime: SimulationTime,
        ) -> Arc<RangeTimeControl> {
            RangeTimeControl::new(
                ScheduleOrControl::Schedule(schedule_lower),
                ScheduleOrControl::Schedule(schedule_upper),
                simtime.iter(),
                0,
                1.0,
                None,
            )
            .unwrap()
            .into()
        }

        #[rstest]
        fn test_non_overlapping_schedules_pass(simtime: SimulationTime) {
            // specify concrete type that satisfies WaterSupplyBehaviour
            let ctrl_a = make_range_time_control(
                vec![Some(0.2), Some(0.2), None, None],
                vec![Some(0.8), Some(0.8), None, None],
                simtime,
            );
            let ctrl_b = make_range_time_control(
                vec![None, None, Some(0.2), Some(0.2)],
                vec![None, None, Some(0.8), Some(0.8)],
                simtime,
            );
            let sources: IndexMap<ArcStr, HeatBatteryChargingSource> = {
                let mut m = IndexMap::new();
                m.insert(
                    "electric".into(),
                    HeatBatteryChargingSource {
                        source_type: ChargingSourceType::DirectElectric,
                        control: ctrl_a,
                        rated_charge_power: Some(5.0),
                        flow_rate_charging_l_per_min: None,
                        temp_flow_max: 0.,
                        hex_a: None,
                        hex_b: None,
                        schedule_unit: ScheduleUnit::StateOfCharge,
                        hex_velocity_at_1_l_per_min: None,
                        hex_capillary_diameter_m: None,
                        heat_source_service: Option::<HeatSourceWetService>::None,
                    },
                );
                m.insert(
                    "hydronic".into(),
                    HeatBatteryChargingSource {
                        source_type: ChargingSourceType::HeatSourceWet,
                        control: ctrl_b,
                        temp_flow_max: 65.0,
                        flow_rate_charging_l_per_min: Some(10.0),
                        hex_a: Some(174.33952),
                        hex_b: Some(-931.565),
                        hex_velocity_at_1_l_per_min: Some(0.035),
                        hex_capillary_diameter_m: Some(6.5 / 1000.0),
                        heat_source_service: Option::<HeatSourceWetService>::None,
                        schedule_unit: Default::default(),
                        rated_charge_power: None,
                    },
                );
                m
            };
            // Should not raise
            validate_no_schedule_overlap(sources, "test_battery", &simtime.iter()).unwrap();
        }

        /// Overlapping schedules (both active at t1) should return an error.
        #[rstest]
        fn test_overlapping_schedule_raises(simtime: SimulationTime) {
            let ctrl_a = make_range_time_control(
                vec![Some(0.2), Some(0.2), None, None],
                vec![Some(0.8), Some(0.8), None, None],
                simtime,
            );
            let ctrl_b = make_range_time_control(
                vec![None, Some(0.2), Some(0.2), None],
                vec![None, Some(0.8), Some(0.8), None],
                simtime,
            );
            let sources: IndexMap<ArcStr, HeatBatteryChargingSource> = {
                let mut m = IndexMap::new();
                m.insert(
                    "electric".into(),
                    HeatBatteryChargingSource {
                        source_type: ChargingSourceType::DirectElectric,
                        control: ctrl_a,
                        rated_charge_power: Some(5.0),
                        flow_rate_charging_l_per_min: None,
                        heat_source_service: Option::<HeatSourceWetService>::None,
                        hex_a: None,
                        temp_flow_max: 0.,
                        hex_b: None,
                        hex_velocity_at_1_l_per_min: None,
                        hex_capillary_diameter_m: None,
                        schedule_unit: Default::default(),
                    },
                );
                m.insert(
                    "hydronic".into(),
                    HeatBatteryChargingSource {
                        source_type: ChargingSourceType::HeatSourceWet,
                        control: ctrl_b,
                        rated_charge_power: Some(3.0),
                        flow_rate_charging_l_per_min: None,
                        hex_a: Some(174.33952),
                        temp_flow_max: 0.,
                        hex_b: Some(-931.565),
                        hex_velocity_at_1_l_per_min: Some(0.035),
                        hex_capillary_diameter_m: Some(6.5 / 1000.0),
                        heat_source_service: Option::<HeatSourceWetService>::None,
                        schedule_unit: Default::default(),
                    },
                );
                m
            };
            assert!(validate_no_schedule_overlap(sources, "test", &simtime.iter()).is_err());
        }
        /// Single source can never overlap — validation accepts it.
        #[rstest]
        fn single_source_always_passes(simtime: SimulationTime) {
            let ctrl_a = make_range_time_control(
                vec![Some(0.2), Some(0.2), Some(0.2), Some(0.2)],
                vec![Some(0.8), Some(0.8), Some(0.8), Some(0.8)],
                simtime,
            );
            let sources: IndexMap<ArcStr, HeatBatteryChargingSource> = {
                let mut m = IndexMap::new();
                m.insert(
                    "a".into(),
                    HeatBatteryChargingSource {
                        source_type: ChargingSourceType::DirectElectric,
                        control: ctrl_a,
                        rated_charge_power: Some(5.0),
                        flow_rate_charging_l_per_min: None,
                        hex_a: None,
                        temp_flow_max: 0.,
                        hex_b: None,
                        hex_velocity_at_1_l_per_min: None,
                        hex_capillary_diameter_m: None,
                        heat_source_service: Option::<HeatSourceWetService>::None,
                        schedule_unit: Default::default(),
                    },
                );
                m
            };
            assert!(validate_no_schedule_overlap(sources, "test", &simtime.iter()).is_ok());
        }

        // skipped test_simtime_reset_after_validation and test_simtime_reset_on_error
        //from Python as we don't have a mutable reference to SimulationTime in Rust.

        ///Transition period (lower=None, upper=non-null) counts as active for overlap.
        /// A transition period means existing charging may continue, so it's an
        /// active period from the overlap perspective.
        #[rstest]
        fn test_transition_period_counts_as_active(simtime: SimulationTime) {
            // Source A: fully active all timesteps
            let ctrl_a = make_range_time_control(
                vec![Some(0.2), Some(0.2), Some(0.2), Some(0.2)],
                vec![Some(0.8), Some(0.8), Some(0.8), Some(0.8)],
                simtime,
            );
            // Source B: transition at t2 (lower=None, upper=Some(0.8))
            let ctrl_b = make_range_time_control(
                vec![None, None, None, None],
                vec![None, None, Some(0.8), None],
                simtime,
            );
            let sources: IndexMap<ArcStr, HeatBatteryChargingSource> = {
                let mut m = IndexMap::new();
                m.insert(
                    "a".into(),
                    HeatBatteryChargingSource {
                        source_type: ChargingSourceType::DirectElectric,
                        control: ctrl_a,
                        rated_charge_power: Some(5.0),
                        flow_rate_charging_l_per_min: None,
                        hex_a: None,
                        temp_flow_max: 0.,
                        hex_b: None,
                        hex_velocity_at_1_l_per_min: None,
                        hex_capillary_diameter_m: None,
                        heat_source_service: Option::<HeatSourceWetService>::None,
                        schedule_unit: Default::default(),
                    },
                );
                m.insert(
                    "b".into(),
                    HeatBatteryChargingSource {
                        source_type: ChargingSourceType::DirectElectric,
                        control: ctrl_b,
                        rated_charge_power: Some(3.0),
                        flow_rate_charging_l_per_min: None,
                        hex_a: None,
                        temp_flow_max: 0.,
                        hex_b: None,
                        hex_velocity_at_1_l_per_min: None,
                        hex_capillary_diameter_m: None,
                        heat_source_service: Option::<HeatSourceWetService>::None,
                        schedule_unit: Default::default(),
                    },
                );
                m
            };
            assert!(validate_no_schedule_overlap(sources, "test", &simtime.iter()).is_err());
        }
    }
}
