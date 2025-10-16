#![allow(unused_variables)]

use crate::core::heating_systems::heat_pump::TestLetter;
use crate::core::schedule::{BooleanSchedule, NumericSchedule};
use crate::corpus::Corpus;
use crate::external_conditions::{DaylightSavingsConfig, ShadingSegment, WindowShadingObject};
use crate::simulation_time::SimulationTime;
use anyhow::{anyhow, bail};
use indexmap::IndexMap;
use itertools::Itertools;
use jsonschema::{BasicOutput, Validator};
use serde::de::Error;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use serde_enum_str::{Deserialize_enum_str, Serialize_enum_str};
use serde_json::{json, Map, Value as JsonValue};
use serde_repr::{Deserialize_repr, Serialize_repr};
use serde_valid::json::ToJsonString;
use serde_valid::validation::error::{Format, Message};
use serde_valid::{MinimumError, Validate};
use smartstring::alias::String;
use std::collections::HashSet;
use std::fmt::{Display, Formatter};
use std::io::{BufReader, Read};
use std::ops::Index;
use std::sync::Arc;
use std::sync::LazyLock;
use thiserror::Error;

pub fn ingest_for_processing(json: impl Read) -> Result<InputForProcessing, anyhow::Error> {
    InputForProcessing::init_with_json(json)
}

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(rename_all = "PascalCase", deny_unknown_fields)]
pub struct Input {
    #[serde(rename = "temp_internal_air_static_calcs")]
    pub(crate) temp_internal_air_static_calcs: f64,
    pub(crate) simulation_time: SimulationTime,
    pub(crate) external_conditions: Arc<ExternalConditionsInput>,
    pub(crate) internal_gains: InternalGains,
    #[serde(default)]
    pub(crate) appliance_gains: ApplianceGains,
    pub(crate) cold_water_source: ColdWaterSourceInput,
    #[serde(default)]
    #[validate(custom = validate_only_storage_tanks)]
    pub(crate) pre_heated_water_source: IndexMap<String, HotWaterSourceDetails>,
    pub(crate) energy_supply: EnergySupplyInput,
    pub(crate) control: Control,
    #[serde(default, skip_serializing_if = "IndexMap::is_empty")]
    pub(crate) smart_appliance_controls: IndexMap<String, SmartApplianceControlDetails>,
    pub(crate) hot_water_source: HotWaterSource,
    pub(crate) hot_water_demand: HotWaterDemand,
    #[serde(rename = "Events")]
    pub(crate) water_heating_events: WaterHeatingEvents,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) space_heat_system: Option<SpaceHeatSystem>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) space_cool_system: Option<SpaceCoolSystem>,
    pub(crate) zone: ZoneDictionary,
    // following fields marked as possibly dead code are likely to be used by wrappers, but worth checking when compiling input schema
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) heat_source_wet: Option<HeatSourceWet>,
    #[serde(rename = "WWHRS", skip_serializing_if = "Option::is_none")]
    pub(crate) waste_water_heat_recovery: Option<WasteWaterHeatRecovery>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) on_site_generation: Option<OnSiteGeneration>,
    pub(crate) infiltration_ventilation: InfiltrationVentilation,
}

fn validate_only_storage_tanks(
    sources: &IndexMap<String, HotWaterSourceDetails>,
) -> Result<(), serde_valid::validation::Error> {
    sources.values()
        .all(|details| matches!(details, HotWaterSourceDetails::StorageTank { .. }))
        .then_some(())
        .ok_or_else(|| serde_valid::validation::Error::Custom("PreHeatedWaterSource input can only contain HotWaterSource data of the type StorageTank".to_owned()))
}

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub(crate) struct ExternalConditionsInput {
    /// List of external air temperatures, one entry per hour (unit: ˚C)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) air_temperatures: Option<Vec<f64>>,
    /// List of wind speeds, one entry per hour (unit: m/s)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) wind_speeds: Option<Vec<f64>>,
    /// List of wind directions in degrees where North=0, East=90, South=180, West=270. Values range: 0 to 360. Wind direction is reported by the direction from which it originates, e.g. a southerly (180 degree) wind blows from the south to the north. (unit: ˚)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) wind_directions: Option<Vec<f64>>,
    /// List of diffuse horizontal radiation values, one entry per hour (unit: W/m²)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) diffuse_horizontal_radiation: Option<Vec<f64>>,
    /// List of direct beam radiation values, one entry per hour (unit: W/m²)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) direct_beam_radiation: Option<Vec<f64>>,
    /// List of ground reflectivity values, 0 to 1, one entry per hour
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) solar_reflectivity_of_ground: Option<Vec<f64>>,
    /// Latitude of weather station, angle from south (unit: ˚)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) latitude: Option<f64>,
    /// Longitude of weather station, easterly +ve westerly -ve (unit: ˚)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) longitude: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) timezone: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) start_day: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) end_day: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) time_series_step: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) january_first: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) daylight_savings: Option<DaylightSavingsConfig>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) leap_day_included: Option<bool>,
    /// A flag to indicate whether direct beam radiation from climate data needs to be converted from horizontal to normal incidence; if normal direct beam radiation values are provided then no conversion is needed
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) direct_beam_conversion_needed: Option<bool>,
    /// Data splitting the ground plane into segments (8-36) and giving height and distance to shading objects surrounding the building
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) shading_segments: Option<Vec<ShadingSegment>>,
}

#[derive(Clone, Debug, Default, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub(crate) struct InternalGains {
    #[serde(
        alias = "total internal gains",
        skip_serializing_if = "Option::is_none"
    )]
    #[validate]
    pub total_internal_gains: Option<InternalGainsDetails>,
    #[serde(rename = "metabolic gains", skip_serializing_if = "Option::is_none")]
    pub metabolic_gains: Option<InternalGainsDetails>,
    #[serde(rename = "EvaporativeLosses", skip_serializing_if = "Option::is_none")]
    pub evaporative_losses: Option<InternalGainsDetails>,
    #[serde(rename = "ColdWaterLosses", skip_serializing_if = "Option::is_none")]
    pub cold_water_losses: Option<InternalGainsDetails>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub other: Option<InternalGainsDetails>,
}

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub(crate) struct InternalGainsDetails {
    #[validate(minimum = 0)]
    #[validate(maximum = 365)]
    pub(crate) start_day: u32,
    pub(crate) time_series_step: f64,
    pub(crate) schedule: NumericSchedule,
}

pub(crate) type ApplianceGains = IndexMap<String, ApplianceGainsDetails>;

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub(crate) struct ApplianceGainsDetails {
    /// First day of the time series, day of the year, 0 to 365
    #[validate(minimum = 0)]
    #[validate(maximum = 365)]
    pub(crate) start_day: u32,
    /// Timestep of the time series data (unit: hours)
    pub(crate) time_series_step: f64,
    /// Proportion of appliance demand turned into heat gains (no unit)
    pub(crate) gains_fraction: f64,
    #[serde(rename = "EnergySupply")]
    pub(crate) energy_supply: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) schedule: Option<NumericSchedule>,
    /// Appliance power consumption when not in use (unit: W)
    #[serde(rename = "Standby", skip_serializing_if = "Option::is_none")]
    pub(crate) standby: Option<f64>,
    #[serde(rename = "Events", skip_serializing_if = "Option::is_none")]
    pub(crate) events: Option<Vec<ApplianceGainsEvent>>,
    #[serde(rename = "loadshifting", skip_serializing_if = "Option::is_none")]
    pub(crate) load_shifting: Option<ApplianceLoadShifting>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) priority: Option<isize>,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct ApplianceGainsEvent {
    pub start: f64,
    pub duration: f64,
    #[serde(rename = "demand_W")]
    pub demand_w: f64,
}

pub(crate) type EnergySupplyInput = IndexMap<String, EnergySupplyDetails>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields, rename_all = "snake_case")]
pub(crate) struct EnergySupplyDetails {
    pub(crate) fuel: FuelType,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) diverter: Option<EnergyDiverter>,
    #[serde(rename = "ElectricBattery", skip_serializing_if = "Option::is_none")]
    pub(crate) electric_battery: Option<ElectricBattery>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) factor: Option<CustomEnergySourceFactor>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) priority: Option<Vec<EnergySupplyPriorityEntry>>,
    /// Denotes that this energy supply can export its surplus supply
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) is_export_capable: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) threshold_charges: Option<[f64; 12]>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) threshold_prices: Option<[f64; 12]>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) tariff: Option<EnergySupplyTariff>,
}

impl EnergySupplyDetails {
    pub(crate) fn with_fuel(fuel_type: FuelType) -> Self {
        Self {
            fuel: fuel_type,
            diverter: None,
            electric_battery: None,
            factor: None,
            priority: None,
            is_export_capable: None,
            threshold_charges: None,
            threshold_prices: None,
            tariff: None,
        }
    }
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum EnergySupplyTariff {
    #[serde(rename = "Standard Tariff")]
    Standard,
    #[serde(rename = "7-Hour Off Peak Tariff")]
    SevenHourOffPeak,
    #[serde(rename = "10-Hour Off Peak Tariff")]
    TenHourOffPeak,
    #[serde(rename = "Variable Time of Day Tariff")]
    VariableTimeOfDay,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum EnergySupplyPriorityEntry {
    ElectricBattery,
    #[serde(rename = "diverter")]
    Diverter,
}

// It's not completely clear at the moment what the difference between fuel type and energy supply type is,
// but electricity and gas each seem to be indicated using different strings between fuel and energy supply
// in the input examples, so keeping them separate for the time being
// (It's also hard to see some of these as types of fuel)
// We expect to reach more clarity when the input specification is finalised.
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "snake_case")]
pub(crate) enum FuelType {
    Electricity,
    MainsGas,
    Custom,
    #[serde(rename = "LPG_bulk")]
    LpgBulk,
    #[serde(rename = "LPG_bottled")]
    LpgBottled,
    #[serde(rename = "LPG_condition_11F")]
    LpgCondition11F,
    UnmetDemand,
    EnergyFromEnvironment,
}

impl From<&EnergySupplyDetails> for FuelType {
    fn from(value: &EnergySupplyDetails) -> Self {
        value.fuel
    }
}

impl Display for FuelType {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            serde_json::to_value(self).unwrap().as_str().unwrap()
        )
    }
}

impl From<FuelType> for String {
    fn from(value: FuelType) -> Self {
        serde_json::to_value(value)
            .unwrap()
            .as_str()
            .unwrap()
            .into()
    }
}

impl TryFrom<EnergySupplyType> for FuelType {
    type Error = anyhow::Error;

    fn try_from(value: EnergySupplyType) -> Result<Self, Self::Error> {
        Ok(match value {
            EnergySupplyType::Electricity => FuelType::Electricity,
            EnergySupplyType::MainsGas => FuelType::MainsGas,
            EnergySupplyType::Custom => FuelType::Custom,
            EnergySupplyType::LpgBulk => FuelType::LpgBulk,
            EnergySupplyType::UnmetDemand => FuelType::UnmetDemand,
            _ => {
                bail!("No fuel type defined to map the energy supply type {value:?}")
            }
        })
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum EnergySupplyType {
    #[serde(rename = "mains elec")]
    Electricity,
    #[serde(rename = "mains gas", alias = "mains_gas")]
    MainsGas,
    #[serde(rename = "_unmet_demand")]
    UnmetDemand,
    #[serde(rename = "custom")]
    Custom,
    #[serde(rename = "bulk LPG")]
    LpgBulk,
    #[serde(rename = "LPG_bottled")]
    LpgBottled,
    #[serde(rename = "LPG_condition_11F")]
    LpgCondition11F,
    #[serde(rename = "heat network")]
    HeatNetwork,
    #[cfg(feature = "fhs")]
    #[serde(rename = "_notional_heat_network")]
    NotionalHeatNetwork,
}

impl Display for EnergySupplyType {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            serde_json::to_value(self).unwrap().as_str().unwrap()
        )
    }
}

impl From<EnergySupplyType> for String {
    fn from(value: EnergySupplyType) -> Self {
        serde_json::to_value(value)
            .unwrap()
            .as_str()
            .unwrap()
            .into()
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "PascalCase")]
#[serde(deny_unknown_fields)]
pub(crate) struct EnergyDiverter {
    pub(crate) heat_source: DiverterHeatSourceType,
    #[serde(rename = "Controlmax", skip_serializing_if = "Option::is_none")]
    pub(crate) control_max: Option<String>,
}

impl Default for EnergyDiverter {
    fn default() -> Self {
        Self {
            heat_source: DiverterHeatSourceType::Immersion,
            control_max: None,
        }
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum StorageTankType {
    #[serde(rename = "hw cylinder")]
    HotWaterCylinder,
}

impl StorageTankType {
    // implementation here could be derived via serde stuff, but keeping simple/ duplicated for now
    pub fn matches(&self, type_string: &str) -> bool {
        match self {
            StorageTankType::HotWaterCylinder => type_string == "hw cylinder",
        }
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum DiverterHeatSourceType {
    #[serde(rename = "immersion")]
    Immersion,
}

impl DiverterHeatSourceType {
    // implementation here could be derived via serde stuff, but keeping simple/ duplicated for now
    pub fn matches(&self, type_string: &str) -> bool {
        match self {
            DiverterHeatSourceType::Immersion => type_string == "immersion",
        }
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct ElectricBattery {
    pub capacity: f64,
    pub charge_discharge_efficiency_round_trip: f64,
    pub battery_age: f64,
    pub minimum_charge_rate_one_way_trip: f64,
    pub maximum_charge_rate_one_way_trip: f64,
    pub maximum_discharge_rate_one_way_trip: f64,
    pub battery_location: BatteryLocation,
    // #[serde(skip_serializing_if = "Option::is_none")]
    pub grid_charging_possible: bool,
}

#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "lowercase")]
pub enum BatteryLocation {
    Inside,
    Outside,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct CustomEnergySourceFactor {
    #[serde(rename = "Emissions Factor kgCO2e/kWh")]
    pub emissions_factor_kg_co2e_k_wh: f64,
    #[serde(rename = "Emissions Factor kgCO2e/kWh including out-of-scope emissions")]
    pub emissions_factor_kg_co2e_k_wh_including_out_of_scope_emissions: f64,
    #[serde(rename = "Primary Energy Factor kWh/kWh delivered")]
    pub primary_energy_factor_k_wh_k_wh_delivered: f64,
}

pub(crate) type ColdWaterSourceInput = IndexMap<ColdWaterSourceType, ColdWaterSourceDetails>;

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub struct ColdWaterSourceDetails {
    /// First day of the time series, day of the year, 0 to 365
    #[validate(minimum = 0)]
    #[validate(maximum = 365)]
    pub(crate) start_day: u32,
    /// List of cold water temperatures, one entry per hour (unit: ˚C)
    pub(crate) temperatures: Vec<f64>,
    /// Timestep of the time series data (unit: hours)
    pub(crate) time_series_step: f64,
}

pub(crate) type ExtraControls = IndexMap<String, ControlDetails>;

#[derive(Clone, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
pub(crate) struct Control {
    #[serde(skip_serializing_if = "Option::is_none", rename = "hw timer")]
    pub(crate) hot_water_timer: Option<ControlDetails>,
    #[serde(skip_serializing_if = "Option::is_none", rename = "window opening")]
    pub(crate) window_opening: Option<ControlDetails>,
    #[serde(flatten)]
    pub(crate) extra: ExtraControls,
}

impl Control {
    pub(crate) fn get(&self, key: &str) -> Option<&ControlDetails> {
        match key {
            "hw timer" => self.hot_water_timer.as_ref(),
            "window opening" => self.window_opening.as_ref(),
            reference => self.extra.get(reference),
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(tag = "type", deny_unknown_fields)]
pub(crate) enum ControlDetails {
    #[serde(rename = "OnOffTimeControl")]
    OnOffTimer {
        #[serde(skip_serializing_if = "Option::is_none")]
        allow_null: Option<bool>,
        /// First day of the time series, day of the year, 0 to 365
        #[validate(minimum = 0)]
        #[validate(maximum = 365)]
        start_day: u32,
        time_series_step: f64,
        schedule: BooleanSchedule,
    },
    #[serde(rename = "OnOffCostMinimisingTimeControl")]
    OnOffCostMinimising {
        /// First day of the time series, day of the year, 0 to 365
        #[validate(minimum = 0)]
        #[validate(maximum = 365)]
        start_day: u32,
        /// Timestep of the time series data (unit: hours)
        time_series_step: f64,
        /// Number of 'on' hours to be set per day
        time_on_daily: f64,
        schedule: NumericSchedule,
    },
    #[serde(rename = "SetpointTimeControl")]
    SetpointTimer {
        /// First day of the time series, day of the year, 0 to 365
        #[validate(minimum = 0)]
        #[validate(maximum = 365)]
        start_day: u32,
        /// Timestep of the time series data (unit: hours)
        time_series_step: f64,
        /// How long before heating period the system should switch on (unit: hours)
        #[serde(skip_serializing_if = "Option::is_none")]
        advanced_start: Option<f64>,
        /// Minimum setpoint allowed
        #[serde(skip_serializing_if = "Option::is_none")]
        setpoint_min: Option<f64>,
        /// Maximum setpoint allowed
        #[serde(skip_serializing_if = "Option::is_none")]
        setpoint_max: Option<f64>,
        /// If both min and max limits are set but setpoint is not, whether to default to min (false) or max (true)
        #[serde(skip_serializing_if = "Option::is_none")]
        default_to_max: Option<bool>,
        /// list of float values (one entry per hour)
        schedule: NumericSchedule,
    },
    #[serde(rename = "ChargeControl")]
    ChargeTarget {
        /// First day of the time series, day of the year, 0 to 365
        #[validate(minimum = 0)]
        #[validate(maximum = 365)]
        start_day: u32,
        time_series_step: f64,
        #[serde(skip_serializing_if = "Option::is_none")]
        logic_type: Option<ControlLogicType>,
        /// Proportion of the charge targeted for each day
        #[serde(skip_serializing_if = "Option::is_none")]
        charge_level: Option<ChargeLevel>,
        #[serde(skip_serializing_if = "Option::is_none")]
        external_sensor: Option<ExternalSensor>,
        #[serde(skip_serializing_if = "Option::is_none")]
        temp_charge_cut: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        temp_charge_cut_delta: Option<NumericSchedule>,
        /// List of boolean values where true means 'on' (one entry per hour)
        schedule: BooleanSchedule,
    },
    #[serde(rename = "CombinationTimeControl")]
    CombinationTime { combination: ControlCombinations },
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(untagged)]
pub(crate) enum ChargeLevel {
    Single(f64),
    List(Vec<f64>),
    Schedule(NumericSchedule),
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct ExternalSensor {
    pub(crate) correlation: Vec<ExternalSensorCorrelation>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct ExternalSensorCorrelation {
    pub(crate) temperature: f64,
    pub(crate) max_charge: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) struct ControlCombinations {
    pub(crate) main: ControlCombination,
    #[serde(flatten)]
    pub(crate) references: IndexMap<String, ControlCombination>,
}

pub(crate) const MAIN_REFERENCE: &str = "main";

impl ControlCombinations {
    pub(crate) fn contains_key(&self, key: &str) -> bool {
        if key == MAIN_REFERENCE {
            return true;
        }
        self.references.contains_key(key)
    }
}

impl Index<&str> for ControlCombinations {
    type Output = ControlCombination;

    fn index(&self, index: &str) -> &Self::Output {
        if index == MAIN_REFERENCE {
            return &self.main;
        }
        &self.references[index]
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct ControlCombination {
    pub(crate) operation: ControlCombinationOperation,
    // unable currently to use serde_valid built-in validations to validate length of string from
    // smartstring crate, so using custom validation instead
    #[validate(custom = validate_length_minimum_one)]
    pub(crate) controls: Vec<String>,
}

fn validate_length_minimum_one(sources: &[String]) -> Result<(), serde_valid::validation::Error> {
    validate_length_minimum::<1>(sources)
}

fn validate_length_minimum<const T: usize>(
    sources: &[String],
) -> Result<(), serde_valid::validation::Error> {
    sources
        .iter()
        .all(|string| string.len() >= T)
        .then_some(())
        .ok_or_else(|| {
            serde_valid::validation::Error::Minimum(Message::new(
                MinimumError::new(2),
                Format::Default,
            ))
        })
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields, rename_all = "UPPERCASE")]
pub(crate) enum ControlCombinationOperation {
    And,
    Or,
    Xor,
    Not,
    Max,
    Min,
    Mean,
}

#[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct SmartApplianceBattery {
    pub(crate) battery_state_of_charge: IndexMap<String, Vec<f64>>,
    #[serde(default)]
    pub(crate) energy_into_battery_from_generation: IndexMap<String, Vec<f64>>,
    #[serde(default)]
    pub(crate) energy_into_battery_from_grid: IndexMap<String, Vec<f64>>,
    #[serde(default)]
    pub(crate) energy_out_of_battery: IndexMap<String, Vec<f64>>,
}

impl ControlDetails {
    pub(crate) fn start_day(&self) -> anyhow::Result<u32> {
        match self {
            ControlDetails::OnOffTimer { start_day, .. } => Ok(*start_day),
            ControlDetails::OnOffCostMinimising { start_day, .. } => Ok(*start_day),
            ControlDetails::SetpointTimer { start_day, .. } => Ok(*start_day),
            ControlDetails::ChargeTarget { start_day, .. } => Ok(*start_day),
            _ => Err(anyhow!("Start day was not available")),
        }
    }

    pub(crate) fn time_series_step(&self) -> anyhow::Result<f64> {
        match self {
            ControlDetails::OnOffTimer {
                time_series_step, ..
            } => Ok(*time_series_step),
            ControlDetails::OnOffCostMinimising {
                time_series_step, ..
            } => Ok(*time_series_step),
            ControlDetails::SetpointTimer {
                time_series_step, ..
            } => Ok(*time_series_step),
            ControlDetails::ChargeTarget {
                time_series_step, ..
            } => Ok(*time_series_step),
            _ => Err(anyhow!("Time series step was not available")),
        }
    }

    pub(crate) fn numeric_schedule(&self) -> anyhow::Result<&NumericSchedule> {
        match self {
            ControlDetails::OnOffCostMinimising { schedule, .. } => Ok(schedule),
            ControlDetails::SetpointTimer { schedule, .. } => Ok(schedule),
            _ => Err(anyhow!("Numeric schedule was not available")),
        }
    }
}

#[derive(Clone, Copy, Debug, Default, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "snake_case")]
pub(crate) enum ControlLogicType {
    #[default]
    Manual,
    Automatic,
    Celect,
    // high heat retention storage heater
    Hhrsh,
    HeatBattery,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct SmartApplianceControlDetails {
    #[serde(rename = "Appliances")]
    pub(crate) appliances: Vec<ApplianceKey>,
    #[serde(rename = "battery24hr")]
    pub(crate) battery_24hr: SmartApplianceBattery,
    pub(crate) non_appliance_demand_24hr: IndexMap<String, Vec<f64>>,
    pub(crate) power_timeseries: IndexMap<String, Vec<f64>>,
    pub(crate) time_series_step: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct HotWaterSource {
    #[serde(rename = "hw cylinder")]
    pub(crate) hot_water_cylinder: HotWaterSourceDetails,
}

impl HotWaterSource {
    pub(crate) fn as_index_map(&self) -> IndexMap<String, HotWaterSourceDetails> {
        IndexMap::from([("hw cylinder".into(), self.hot_water_cylinder.clone())])
    }
}

#[derive(Clone, Deserialize_enum_str, PartialEq, Debug, Serialize_enum_str)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum BoilerHotWaterTest {
    #[serde(rename = "M&L")]
    ML,
    #[serde(rename = "M&S")]
    MS,
    #[serde(rename = "M_only")]
    MOnly,
    #[serde(rename = "No_additional_tests")]
    NoAdditionalTests,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "type", deny_unknown_fields)]
pub(crate) enum HotWaterSourceDetails {
    StorageTank {
        /// Total volume of tank (unit: litre)
        #[validate(minimum = 0.0)]
        volume: f64,
        /// Measured standby losses due to cylinder insulation at standardised conditions (unit: kWh/24h)
        daily_losses: f64,
        /// Surface area of the heat exchanger within the storage tank (unit: m²)
        #[serde(skip_serializing_if = "Option::is_none")]
        heat_exchanger_surface_area: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        min_temp: Option<f64>,
        /// Initial temperature of the storage tank at the start of simulation (unit: ˚C)
        init_temp: f64,
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceReference,
        /// Map of heating systems connected to the storage tank
        #[serde(rename = "HeatSource")]
        heat_source: IndexMap<String, HeatSource>,
        /// List of primary pipework components connected to the storage tank
        #[serde(skip_serializing_if = "Option::is_none")]
        primary_pipework: Option<Vec<WaterPipework>>,
    },
    CombiBoiler {
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceType,
        #[serde(rename = "HeatSourceWet")]
        heat_source_wet: HeatSourceWetType,
        #[serde(rename = "separate_DHW_tests")]
        separate_dhw_tests: BoilerHotWaterTest,
        #[serde(skip_serializing_if = "Option::is_none")]
        rejected_energy_1: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        storage_loss_factor_2: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        rejected_factor_3: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        setpoint_temp: Option<f64>,
        #[serde(rename = "daily_HW_usage")]
        daily_hw_usage: f64,
    },
    #[serde(rename = "HIU")]
    Hiu {
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceType,
        #[serde(rename = "HeatSourceWet")]
        heat_source_wet: HeatSourceWetType,
        #[serde(skip_serializing_if = "Option::is_none")]
        setpoint_temp: Option<f64>,
    },
    PointOfUse {
        efficiency: f64,
        #[serde(rename = "EnergySupply")]
        energy_supply: String,
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceType,
        setpoint_temp: f64,
    },
    SmartHotWaterTank {
        volume: f64,
        #[serde(rename = "power_pump_kW")]
        power_pump_kw: f64,
        max_flow_rate_pump_l_per_min: f64,
        temp_usable: f64,
        /// Reference to a control schedule of maximum state of charge values
        temp_setpnt_max: String,
        daily_losses: f64,
        init_temp: f64,
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceType,
        #[serde(rename = "EnergySupply_pump")]
        energy_supply_pump: String,
        #[serde(rename = "HeatSource")]
        heat_source: IndexMap<String, HeatSource>,
        #[serde(skip_serializing_if = "Option::is_none")]
        primary_pipework: Option<Vec<WaterPipeworkSimple>>,
    },
    HeatBattery {
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceType,
        #[serde(rename = "HeatSourceWet")]
        heat_source_wet: String,
    },
}

impl HotWaterSourceDetails {
    pub fn volume(&self) -> Option<f64> {
        if let Self::StorageTank { volume, .. } = self {
            Some(*volume)
        } else {
            None
        }
    }
}

pub(crate) trait HotWaterSourceDetailsForProcessing {
    fn is_storage_tank(&self) -> bool;
    fn is_combi_boiler(&self) -> bool;
    fn is_hiu(&self) -> bool;
    fn is_point_of_use(&self) -> bool;
    fn is_smart_hot_water_tank(&self) -> bool;
    fn set_control_min_name_for_storage_tank_heat_sources(
        &mut self,
        control_name: &str,
    ) -> anyhow::Result<()>;
    fn set_control_max_name_for_storage_tank_heat_sources(
        &mut self,
        control_name: &str,
    ) -> anyhow::Result<()>;
    fn set_control_min_name_for_smart_hot_water_tank_heat_sources(
        &mut self,
        control_name: &str,
    ) -> anyhow::Result<()>;
    fn set_control_max_name_for_smart_hot_water_tank_heat_sources(
        &mut self,
        control_name: &str,
    ) -> anyhow::Result<()>;
    fn set_temp_setpoint_max_for_smart_hot_water_tank_heat_sources(
        &mut self,
        temp_setpoint_max_name: &str,
    ) -> anyhow::Result<()>;
}

pub(crate) struct HotWaterSourceDetailsJsonMap<'a>(
    pub(crate) &'a mut Map<std::string::String, JsonValue>,
);

impl HotWaterSourceDetailsForProcessing for HotWaterSourceDetailsJsonMap<'_> {
    fn is_storage_tank(&self) -> bool {
        self.0
            .get("type")
            .and_then(|source_type| source_type.as_str())
            .is_some_and(|source_type| source_type == "StorageTank")
    }

    fn is_combi_boiler(&self) -> bool {
        self.0
            .get("type")
            .and_then(|source_type| source_type.as_str())
            .is_some_and(|source_type| source_type == "CombiBoiler")
    }

    fn is_hiu(&self) -> bool {
        self.0
            .get("type")
            .and_then(|source_type| source_type.as_str())
            .is_some_and(|source_type| source_type == "HIU")
    }

    fn is_point_of_use(&self) -> bool {
        self.0
            .get("type")
            .and_then(|source_type| source_type.as_str())
            .is_some_and(|source_type| source_type == "PointOfUse")
    }

    fn is_smart_hot_water_tank(&self) -> bool {
        self.0
            .get("type")
            .and_then(|source_type| source_type.as_str())
            .is_some_and(|source_type| source_type == "SmartHotWaterTank")
    }

    fn set_control_min_name_for_storage_tank_heat_sources(
        &mut self,
        control_name: &str,
    ) -> anyhow::Result<()> {
        if !self.is_storage_tank() {
            return Ok(());
        }

        if let Some(heat_sources) = self.0.get_mut("HeatSource").and_then(|v| v.as_object_mut()) {
            for heat_source in heat_sources.values_mut().flat_map(|v| v.as_object_mut()) {
                heat_source.insert("Controlmin".into(), json!(control_name));
            }
        }

        Ok(())
    }

    fn set_control_max_name_for_storage_tank_heat_sources(
        &mut self,
        control_name: &str,
    ) -> anyhow::Result<()> {
        if !self.is_storage_tank() {
            return Ok(());
        }

        if let Some(heat_sources) = self.0.get_mut("HeatSource").and_then(|v| v.as_object_mut()) {
            for heat_source in heat_sources.values_mut().flat_map(|v| v.as_object_mut()) {
                heat_source.insert("Controlmax".into(), json!(control_name));
            }
        }

        Ok(())
    }

    fn set_control_min_name_for_smart_hot_water_tank_heat_sources(
        &mut self,
        control_name: &str,
    ) -> anyhow::Result<()> {
        if !self.is_smart_hot_water_tank() {
            return Ok(());
        }

        if let Some(heat_sources) = self.0.get_mut("HeatSource").and_then(|v| v.as_object_mut()) {
            for heat_source in heat_sources.values_mut().flat_map(|v| v.as_object_mut()) {
                heat_source.insert("Controlmin".into(), json!(control_name));
            }
        }

        Ok(())
    }

    fn set_control_max_name_for_smart_hot_water_tank_heat_sources(
        &mut self,
        control_name: &str,
    ) -> anyhow::Result<()> {
        if !self.is_smart_hot_water_tank() {
            return Ok(());
        }

        if let Some(heat_sources) = self.0.get_mut("HeatSource").and_then(|v| v.as_object_mut()) {
            for heat_source in heat_sources.values_mut().flat_map(|v| v.as_object_mut()) {
                heat_source.insert("Controlmax".into(), json!(control_name));
            }
        }

        Ok(())
    }

    fn set_temp_setpoint_max_for_smart_hot_water_tank_heat_sources(
        &mut self,
        temp_setpoint_max_name: &str,
    ) -> anyhow::Result<()> {
        if !self.is_smart_hot_water_tank() {
            return Ok(());
        }

        if let Some(heat_sources) = self.0.get_mut("HeatSource").and_then(|v| v.as_object_mut()) {
            for heat_source in heat_sources.values_mut().flat_map(|v| v.as_object_mut()) {
                heat_source.insert("temp_setpnt_max".into(), json!(temp_setpoint_max_name));
            }
        }

        Ok(())
    }
}

#[derive(Copy, Clone, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum ColdWaterSourceType {
    #[serde(rename = "mains water")]
    MainsWater,
    #[serde(rename = "header tank")]
    HeaderTank,
}

#[derive(Clone, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(untagged)]
pub enum ColdWaterSourceReference {
    Type(ColdWaterSourceType),
    Reference(String),
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum HeatSourceWetType {
    #[serde(rename = "boiler")]
    Boiler,
    HeatNetwork,
    #[serde(rename = "hp")]
    HeatPump,
    #[serde(untagged)]
    Other(String),
}

impl HeatSourceWetType {
    /// Convert the type to a canonical string based on the input format to be used in e.g. energy supply names
    pub fn to_canonical_string(&self) -> String {
        serde_json::to_value(self).unwrap().as_str().unwrap().into()
    }
}

#[derive(Clone, Copy, Debug, Deserialize_enum_str, PartialEq, Serialize_enum_str)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum HeatSourceControlType {
    #[serde(rename = "hw timer")]
    HotWaterTimer,
    #[serde(rename = "window opening")]
    WindowOpening,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "type", deny_unknown_fields)]
pub(crate) enum HeatSource {
    ImmersionHeater {
        power: f64,
        #[serde(rename = "EnergySupply")]
        energy_supply: String,
        /// Reference to a control schedule of minimum temperature setpoints
        #[serde(rename = "Controlmin", skip_serializing_if = "Option::is_none")]
        control_min: Option<String>,
        /// Reference to a control schedule of maximum temperature setpoints
        #[serde(rename = "Controlmax", skip_serializing_if = "Option::is_none")]
        control_max: Option<String>,
        heater_position: f64,
        #[serde(skip_serializing_if = "Option::is_none")]
        thermostat_position: Option<f64>,
    },
    SolarThermalSystem {
        #[serde(rename = "sol_loc")]
        solar_cell_location: SolarCollectorLoopLocation,
        area_module: f64,
        #[validate(minimum = 1)]
        modules: usize,
        peak_collector_efficiency: f64,
        incidence_angle_modifier: f64,
        first_order_hlc: f64,
        second_order_hlc: f64,
        collector_mass_flow_rate: f64,
        power_pump: f64,
        power_pump_control: f64,
        #[serde(rename = "EnergySupply")]
        energy_supply: String,
        tilt: f64,
        #[serde(
            rename = "orientation360",
            deserialize_with = "deserialize_orientation",
            serialize_with = "serialize_orientation"
        )]
        orientation: f64,
        solar_loop_piping_hlc: f64,
        heater_position: f64,
        /// Required for StorageTank but not for SmartHotWaterTank
        #[serde(skip_serializing_if = "Option::is_none")]
        thermostat_position: Option<f64>,
        /// Reference to a control schedule of maximum temperature setpoints
        #[serde(rename = "Controlmax")]
        control_max: String,
    },
    #[serde(rename = "HeatSourceWet")]
    ServiceWaterRegular {
        name: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        temp_flow_limit_upper: Option<f64>,
        #[serde(rename = "EnergySupply")]
        energy_supply: String,
        /// Reference to a control schedule of minimum temperature setpoints
        #[serde(rename = "Controlmin", skip_serializing_if = "Option::is_none")]
        control_min: Option<String>,
        /// Reference to a control schedule of maximum temperature setpoints
        #[serde(rename = "Controlmax", skip_serializing_if = "Option::is_none")]
        control_max: Option<String>,
        heater_position: f64,
        /// Required for StorageTank but not for SmartHotWaterTank
        #[serde(skip_serializing_if = "Option::is_none")]
        thermostat_position: Option<f64>,
    },
    #[serde(rename = "HeatPump_HWOnly")]
    HeatPumpHotWaterOnly {
        power_max: f64,
        vol_hw_daily_average: f64,
        tank_volume_declared: f64,
        heat_exchanger_surface_area_declared: f64,
        daily_losses_declared: f64,
        in_use_factor_mismatch: f64,
        test_data: HeatPumpHotWaterTestData,
        #[serde(rename = "EnergySupply")]
        energy_supply: String,
        /// Reference to a control schedule of minimum temperature setpoints
        #[serde(rename = "Controlmin")]
        control_min: String,
        /// Reference to a control schedule of maximum temperature setpoints
        #[serde(rename = "Controlmax")]
        control_max: String,
        heater_position: f64,
        thermostat_position: f64,
    },
}

impl HeatSource {
    pub(crate) fn heater_position(&self) -> f64 {
        match self {
            HeatSource::ImmersionHeater {
                heater_position, ..
            } => *heater_position,
            HeatSource::SolarThermalSystem {
                heater_position, ..
            } => *heater_position,
            HeatSource::ServiceWaterRegular {
                heater_position, ..
            } => *heater_position,
            HeatSource::HeatPumpHotWaterOnly {
                heater_position, ..
            } => *heater_position,
        }
    }

    pub(crate) fn thermostat_position(&self) -> Option<f64> {
        match self {
            HeatSource::ImmersionHeater {
                thermostat_position,
                ..
            } => *thermostat_position,
            HeatSource::SolarThermalSystem {
                thermostat_position,
                ..
            } => *thermostat_position,
            HeatSource::ServiceWaterRegular {
                thermostat_position,
                ..
            } => *thermostat_position,
            HeatSource::HeatPumpHotWaterOnly {
                thermostat_position,
                ..
            } => Some(*thermostat_position),
        }
    }

    pub(crate) fn energy_supply_name(&self) -> &str {
        match self {
            HeatSource::ImmersionHeater { energy_supply, .. } => energy_supply,
            HeatSource::SolarThermalSystem { energy_supply, .. } => energy_supply,
            HeatSource::ServiceWaterRegular { energy_supply, .. } => energy_supply,
            HeatSource::HeatPumpHotWaterOnly { energy_supply, .. } => energy_supply,
        }
    }
}

/// Location of the main part of the solar thermal collector loop piping.
///
/// This affects the ambient temperature used for heat loss calculations
/// in the collector loop piping.
#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum SolarCollectorLoopLocation {
    #[serde(rename = "OUT")]
    Out,
    #[serde(rename = "HS")]
    Hs,
    #[serde(rename = "NHS")]
    Nhs,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct HeatPumpHotWaterTestData {
    #[serde(rename = "L", skip_serializing_if = "Option::is_none")]
    pub(crate) l: Option<HeatPumpHotWaterOnlyTestDatum>,
    #[serde(rename = "M")]
    pub(crate) m: HeatPumpHotWaterOnlyTestDatum,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct HeatPumpHotWaterOnlyTestDatum {
    /// CoP measured during EN 16147 test
    pub(crate) cop_dhw: f64,
    /// daily energy requirement (kWh/day) for tapping profile used for test
    pub(crate) hw_tapping_prof_daily_total: f64,
    /// electrical input energy (kWh) measured in EN 16147 test over 24 hrs
    pub(crate) energy_input_measured: f64,
    /// standby power (W) measured in EN 16147 test
    pub(crate) power_standby: f64,
    /// daily hot water vessel heat loss
    /// (kWh/day) for a 45 K temperature difference between vessel
    /// and surroundings, tested in accordance with BS 1566 or
    /// EN 12897 or any equivalent standard. Vessel must be same
    /// as that used during EN 16147 test
    pub(crate) hw_vessel_loss_daily: f64,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct WaterPipeworkSimple {
    pub location: WaterPipeworkLocation,
    pub internal_diameter_mm: f64,
    pub length: f64,
    // remainder of fields are not in the input definition for Python HEM 0.36 but storage tank logic seems to need them for now
    #[serde(skip_serializing_if = "Option::is_none")]
    pub external_diameter_mm: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub insulation_thermal_conductivity: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub insulation_thickness_mm: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub surface_reflectivity: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pipe_contents: Option<WaterPipeContentsType>,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct WaterPipework {
    pub location: WaterPipeworkLocation,
    pub internal_diameter_mm: f64,
    pub external_diameter_mm: f64,
    pub length: f64,
    pub insulation_thermal_conductivity: f64,
    pub insulation_thickness_mm: f64,
    pub surface_reflectivity: bool,
    pub pipe_contents: WaterPipeContentsType,
}

impl TryFrom<WaterPipeworkSimple> for WaterPipework {
    type Error = anyhow::Error;

    fn try_from(loose: WaterPipeworkSimple) -> Result<Self, Self::Error> {
        Ok(Self {
            location: loose.location,
            internal_diameter_mm: loose.internal_diameter_mm,
            external_diameter_mm: loose
                .external_diameter_mm
                .ok_or_else(|| anyhow!("Missing external_diameter_mm value in water pipework."))?,
            length: loose.length,
            insulation_thermal_conductivity: loose.insulation_thermal_conductivity.ok_or_else(
                || anyhow!("Missing insulation_thermal_conductivity value in water pipework."),
            )?,
            insulation_thickness_mm: loose.insulation_thickness_mm.ok_or_else(|| {
                anyhow!("Missing insulation_thickness_mm value in water pipework.")
            })?,
            surface_reflectivity: loose
                .surface_reflectivity
                .ok_or_else(|| anyhow!("Missing surface_reflectivity value in water pipework."))?,
            pipe_contents: loose
                .pipe_contents
                .ok_or_else(|| anyhow!("Missing pipe_contents value in water pipework."))?,
        })
    }
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum WaterPipeworkLocation {
    #[serde(rename = "internal")]
    Internal,
    #[serde(rename = "external")]
    External,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum WaterPipeContentsType {
    #[serde(rename = "water")]
    Water,
    #[serde(rename = "glycol25")]
    Glycol25,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub(crate) struct HotWaterDemand {
    #[serde(rename = "Shower", skip_serializing_if = "Option::is_none")]
    pub(crate) shower: Option<Showers>,
    #[serde(rename = "Bath", skip_serializing_if = "Option::is_none")]
    pub(crate) bath: Option<Baths>,
    #[serde(rename = "Other", skip_serializing_if = "Option::is_none")]
    pub(crate) other_water_use: Option<OtherWaterUses>,
    #[serde(rename = "Distribution", skip_serializing_if = "Option::is_none")]
    pub(crate) water_distribution: Option<WaterDistribution>,
}

#[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub struct Showers(pub IndexMap<String, Shower>);

impl Showers {
    /// Provide shower field names as strings.
    pub fn keys(&self) -> Vec<String> {
        self.0.keys().cloned().collect()
    }

    pub fn name_refers_to_instant_electric_shower(&self, name: &str) -> bool {
        self.0
            .get(name)
            .is_some_and(|shower| matches!(shower, Shower::InstantElectricShower { .. }))
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields, tag = "type")]
pub enum Shower {
    MixerShower {
        /// Shower flow rate (unit: litre/minute)
        flowrate: f64,
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceType,
        /// Reference to a key in Input.WWHRS
        #[serde(rename = "WWHRS", skip_serializing_if = "Option::is_none")]
        waste_water_heat_recovery_system: Option<String>,
    },
    #[serde(rename = "InstantElecShower")]
    InstantElectricShower {
        rated_power: f64,
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceType,
        #[serde(rename = "EnergySupply")]
        energy_supply: String,
    },
}

#[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct Baths(pub IndexMap<String, BathDetails>);

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct BathDetails {
    /// Volume held by bath (unit: litre)
    pub(crate) size: f64,
    #[serde(rename = "ColdWaterSource")]
    pub(crate) cold_water_source: ColdWaterSourceType,
    /// Tap/outlet flow rate (unit: litre/minute)
    pub(crate) flowrate: f64,
}

#[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct OtherWaterUses(pub IndexMap<String, OtherWaterUse>);

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct OtherWaterUse {
    /// Tap/outlet flow rate (unit: litre/minute)
    pub(crate) flowrate: f64,
    #[serde(rename = "ColdWaterSource")]
    pub(crate) cold_water_source: ColdWaterSourceType,
}

pub(crate) type WaterDistribution = Vec<WaterPipeworkSimple>;

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(rename_all = "PascalCase")]
pub(crate) struct WaterHeatingEvents {
    #[serde(default)]
    pub(crate) shower: IndexMap<String, Vec<WaterHeatingEvent>>,
    #[serde(default)]
    pub(crate) bath: IndexMap<String, Vec<WaterHeatingEvent>>,
    #[serde(default)]
    pub(crate) other: IndexMap<String, Vec<WaterHeatingEvent>>,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub struct WaterHeatingEvent {
    pub start: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub duration: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub volume: Option<f64>,
    pub temperature: f64,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum WaterHeatingEventType {
    Shower,
    Bath,
    Other,
}

pub(crate) type SpaceHeatSystem = IndexMap<String, SpaceHeatSystemDetails>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields, tag = "type")]
pub(crate) enum SpaceHeatSystemDetails {
    #[serde(rename = "InstantElecHeater")]
    InstantElectricHeater {
        rated_power: f64,
        #[serde(rename = "EnergySupply")]
        energy_supply: String,
        #[serde(rename = "Control")]
        control: String,
        /// Convective fraction for heating
        frac_convective: f64,
    },
    #[serde(rename = "ElecStorageHeater")]
    ElectricStorageHeater {
        #[serde(skip_serializing_if = "Option::is_none")]
        advanced_start: Option<f64>,
        pwr_in: f64,
        /// (instant backup) (unit: kW)
        rated_power_instant: f64,
        storage_capacity: f64,
        air_flow_type: ElectricStorageHeaterAirFlowType,
        frac_convective: f64,
        /// Fan power (unit: W)
        fan_pwr: f64,
        n_units: u32,
        #[serde(rename = "EnergySupply")]
        energy_supply: String,
        #[serde(rename = "Control")]
        control: String,
        #[serde(rename = "ControlCharger")]
        control_charger: String,
        /// The zone where the unit(s) is/are installed
        #[serde(rename = "Zone")]
        zone: String,
        #[serde(skip_serializing_if = "Option::is_none")]
        temp_setback: Option<f64>,
        #[serde(rename = "ESH_min_output")]
        esh_min_output: Vec<(f64, f64)>,
        #[serde(rename = "ESH_max_output")]
        esh_max_output: Vec<(f64, f64)>,
    },
    WetDistribution {
        #[serde(skip_serializing_if = "Option::is_none")]
        thermal_mass: Option<f64>,
        #[serde(default)]
        emitters: Vec<WetEmitter>,
        #[serde(rename = "EnergySupply", skip_serializing_if = "Option::is_none")]
        energy_supply: Option<String>,
        temp_diff_emit_dsgn: f64,
        variable_flow: bool,
        #[serde(skip_serializing_if = "Option::is_none")]
        design_flow_rate: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        min_flow_rate: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        max_flow_rate: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        bypass_percentage_recirculated: Option<f64>,
        #[serde(rename = "HeatSource")]
        heat_source: SpaceHeatSystemHeatSource,
        #[serde(rename = "Control")]
        control: String,
        ecodesign_controller: EcoDesignController,
        design_flow_temp: i32,
        #[serde(rename = "Zone")]
        zone: String,
    },
    WarmAir {
        frac_convective: f64,
        #[serde(rename = "HeatSource")]
        heat_source: SpaceHeatSystemHeatSource,
        #[serde(rename = "Control")]
        control: String,
    },
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(
    deny_unknown_fields,
    tag = "wet_emitter_type",
    rename_all = "lowercase"
)]
pub(crate) enum WetEmitter {
    Radiator {
        c: f64,
        n: f64,
        frac_convective: f64,
    },
    Ufh {
        equivalent_specific_thermal_mass: f64,
        system_performance_factor: f64,
        emitter_floor_area: f64,
        frac_convective: f64,
    },
    Fancoil {
        #[serde(default = "default_n_units")]
        n_units: usize,
        frac_convective: f64,
        fancoil_test_data: FancoilTestData,
    },
}

const fn default_n_units() -> usize {
    1
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) struct FancoilTestData {
    pub(crate) fan_speed_data: Vec<FanSpeedData>,
    #[serde(rename = "fan_power_W")]
    pub(crate) fan_power_w: Vec<f64>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) struct FanSpeedData {
    pub(crate) temperature_diff: f64,
    pub(crate) power_output: Vec<f64>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct SpaceHeatSystemHeatSource {
    pub name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub temp_flow_limit_upper: Option<f64>,
}

// it is unclear whether this struct should be used - see reference to the struct above
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[allow(dead_code)]
#[serde(deny_unknown_fields)]
pub(crate) struct EcoDesignController {
    pub(crate) ecodesign_control_class: EcoDesignControllerClass,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) min_outdoor_temp: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) max_outdoor_temp: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) min_flow_temp: Option<f64>,
}

#[derive(Clone, Copy, Debug, Deserialize_repr, PartialEq, Serialize_repr)]
#[repr(u8)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum EcoDesignControllerClass {
    /// On/off room thermostat
    ClassI = 1,
    /// Weather compensator with modulating heaters
    ClassII = 2,
    /// Weather compensator with on/off heaters
    ClassIII = 3,
    /// TPI room thermostat with on/off heaters
    ClassIV = 4,
    /// Modulating room thermostat with modulating heaters
    ClassV = 5,
    /// Weather compensator with room sensor for modulating heaters
    ClassVI = 6,
    /// Weather compensator with room sensor for on/off heaters
    ClassVII = 7,
    /// Multi room temperature control with modulating heaters
    ClassVIII = 8,
}

#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "lowercase")]
pub enum MVHRLocation {
    Inside,
    Outside,
}

#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "kebab-case")]
pub(crate) enum ElectricStorageHeaterAirFlowType {
    FanAssisted,
    DamperOnly,
}

pub(crate) type ZoneDictionary = IndexMap<String, ZoneInput>;

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub(crate) struct ZoneInput {
    /// Heating system details of the zone. References a key in $.SpaceHeatSystem
    #[serde(
        rename = "SpaceHeatSystem",
        skip_serializing_if = "SystemReference::is_none",
        default
    )]
    pub(crate) space_heat_system: SystemReference,
    /// Cooling system details of the zone. References a key in $.SpaceCoolSystem
    #[serde(
        rename = "SpaceCoolSystem",
        skip_serializing_if = "SystemReference::is_none",
        default
    )]
    pub(crate) space_cool_system: SystemReference,
    /// Useful floor area of the zone. (Unit: m²)
    #[validate(minimum = 0.)]
    pub(crate) area: f64,
    /// Total volume of the zone. (Unit: m³)
    #[validate(minimum = 0.)]
    pub(crate) volume: f64,
    /// Basis for zone temperature control.
    #[serde(rename = "temp_setpnt_basis", skip_serializing_if = "Option::is_none")]
    pub(crate) temp_setpnt_basis: Option<ZoneTemperatureControlBasis>,
    /// Setpoint temperature to use during initialisation (unit: ˚C)
    pub(crate) temp_setpnt_init: f64,
    /// Map of building elements present in the zone (e.g. walls, floors, windows, etc.).
    #[serde(rename = "BuildingElement")]
    pub(crate) building_elements: IndexMap<String, BuildingElement>,
    /// Overall heat transfer coefficient of the thermal bridge (in W/K), or dictionary of linear thermal transmittance details of the thermal bridges in the zone.
    #[serde(rename = "ThermalBridging")]
    pub(crate) thermal_bridging: ThermalBridging,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(untagged)]
pub(crate) enum SystemReference {
    None(()),
    Single(String),
    #[validate(unique_items)]
    Multiple(Vec<String>),
}

impl SystemReference {
    fn is_none(&self) -> bool {
        matches!(self, Self::None(_))
    }
}

impl Default for SystemReference {
    fn default() -> Self {
        Self::None(())
    }
}

impl IntoIterator for SystemReference {
    type Item = String;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        match self {
            Self::None(_) => vec![].into_iter(),
            Self::Single(s) => vec![s].into_iter(),
            Self::Multiple(v) => v.into_iter(),
        }
    }
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum SpaceHeatControlType {
    #[serde(rename = "livingroom")]
    LivingRoom,
    #[serde(rename = "restofdwelling")]
    RestOfDwelling,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct ZoneLighting {
    #[serde(rename = "efficacy")]
    efficacy: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    bulbs: Option<IndexMap<String, ZoneLightingBulbs>>,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct ZoneLightingBulbs {
    pub(crate) count: usize,
    pub(crate) power: f64,
    pub(crate) efficacy: f64,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "lowercase")]
pub enum ZoneTemperatureControlBasis {
    // for dry-build temperature
    Air,
    // for operative temperature
    Operative,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(tag = "type")]
pub(crate) enum BuildingElement {
    #[serde(rename = "BuildingElementOpaque")]
    Opaque {
        #[serde(skip_serializing_if = "Option::is_none")]
        is_unheated_pitched_roof: Option<bool>,
        /// Solar absorption coefficient at the external surface (dimensionless)
        solar_absorption_coeff: f64,
        #[serde(skip_serializing_if = "Option::is_none")]
        u_value: Option<f64>,
        /// Thermal resistance (unit: m².K/W)
        #[serde(skip_serializing_if = "Option::is_none")]
        thermal_resistance_construction: Option<f64>,
        /// Areal heat capacity (unit: J/m².K)
        areal_heat_capacity: f64,
        /// Mass distribution class of the building element, one of: evenly distributed (D); concentrated on external side (E); concentrated on internal side (I); concentrated on internal and external sides (IE); concentrated in middle (M)
        mass_distribution_class: MassDistributionClass,
        is_external_door: Option<bool>,
        /// Tilt angle of the surface from horizontal, between 0 and 180, where 0 means the external surface is facing up, 90 means the external surface is vertical and 180 means the external surface is facing down (unit: ˚)
        pitch: f64,
        #[serde(
            rename = "orientation360",
            deserialize_with = "deserialize_orientation",
            serialize_with = "serialize_orientation"
        )]
        orientation: f64,
        /// The distance between the ground and the lowest edge of the element (unit: m)
        base_height: f64,
        /// The height of the building element (unit: m)
        height: f64,
        /// The width of the building element (unit: m)
        width: f64,
        /// Net area of the opaque building element (i.e. minus any windows / doors / etc.) (unit: m²)
        area: f64,
    },
    #[serde(rename = "BuildingElementTransparent")]
    Transparent {
        #[serde(skip_serializing_if = "Option::is_none")]
        u_value: Option<f64>,
        #[serde(
            rename = "Control_WindowOpenable",
            skip_serializing_if = "Option::is_none"
        )]
        window_openable_control: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        thermal_resistance_construction: Option<f64>,
        /// Tilt angle of the surface from horizontal, between 0 and 180, where 0 means the external surface is facing up, 90 means the external surface is vertical and 180 means the external surface is facing down (unit: ˚
        pitch: f64,
        #[serde(
            rename = "orientation360",
            deserialize_with = "deserialize_orientation",
            serialize_with = "serialize_orientation"
        )]
        orientation: f64,
        /// Total solar energy transmittance of the transparent part of the window
        g_value: f64,
        /// The frame area fraction of window, ratio of the projected frame area to the overall projected area of the glazed element of the window
        frame_area_fraction: f64,
        /// The distance between the ground and the lowest edge of the element (unit: m)
        base_height: f64,
        /// The height of the building element (unit: m)
        height: f64,
        /// The width of the building element (unit: m)
        width: f64,
        free_area_height: f64,
        mid_height: f64,
        max_window_open_area: f64,
        window_part_list: Vec<WindowPart>,
        shading: Vec<WindowShadingObject>,
        #[serde(default)]
        treatment: Option<Vec<WindowTreatment>>,
    },
    #[serde(rename = "BuildingElementGround")]
    Ground {
        /// Area of this building element within the zone (unit: m²)
        area: f64,
        /// Total area of the building element across entire dwelling; if the Floor is divided among several zones, this is the total area across all zones (unit: m²)
        total_area: f64,
        /// Tilt angle of the surface from horizontal, between 0 and 180, where 0 means the external surface is facing up, 90 means the external surface is vertical and 180 means the external surface is facing down (unit: ˚)
        pitch: f64,
        /// Steady-state thermal transmittance of floor, including the effect of the ground (calculated for the entire ground floor, even if it is distributed among several zones) (unit: W/m2.K)
        u_value: f64,
        /// Total thermal resistance of all layers in the floor construction (unit: m².K/W)
        thermal_resistance_floor_construction: f64,
        /// Areal heat capacity of the ground floor element (unit: J/m2.K)
        areal_heat_capacity: f64,
        /// Mass distribution class of the building element, one of: evenly distributed (D); concentrated on external side (E); concentrated on internal side (I); concentrated on internal and external sides (IE); concentrated in middle (M)
        mass_distribution_class: MassDistributionClass,
        /// Perimeter of the floor; calculated for the entire ground floor, even if it is distributed among several zones (unit: m)
        perimeter: f64,
        /// Linear thermal transmittance of the junction between the floor and the walls (unit: W/m.K)
        psi_wall_floor_junc: f64,
        /// Thickness of the walls (unit: m)
        thickness_walls: f64,
        #[serde(flatten)]
        floor_data: FloorData,
    },
    #[serde(rename = "BuildingElementAdjacentConditionedSpace")]
    AdjacentConditionedSpace {
        area: f64,
        /// Tilt angle of the surface from horizontal, between 0 and 180, where 0 means the external surface is facing up, 90 means the external surface is vertical and 180 means the external surface is facing down (unit: ˚)
        pitch: f64,
        #[serde(skip_serializing_if = "Option::is_none")]
        u_value: Option<f64>,
        /// Thermal resistance (unit: m².K/W)
        #[serde(skip_serializing_if = "Option::is_none")]
        thermal_resistance_construction: Option<f64>,
        areal_heat_capacity: f64,
        mass_distribution_class: MassDistributionClass,
    },
    #[serde(rename = "BuildingElementAdjacentUnconditionedSpace_Simple")]
    AdjacentUnconditionedSpace {
        /// Area of this building element (unit: m²)
        area: f64,
        /// Tilt angle of the surface from horizontal, between 0 and 180, where 0 means the external surface is facing up, 90 means the external surface is vertical and 180 means the external surface is facing down (unit: ˚)
        pitch: f64,
        #[serde(skip_serializing_if = "Option::is_none")]
        u_value: Option<f64>,
        /// Thermal resistance (unit: m2.K/W)
        #[serde(skip_serializing_if = "Option::is_none")]
        thermal_resistance_construction: Option<f64>,
        /// Effective thermal resistance of unheated space (unit: m².K/W)
        thermal_resistance_unconditioned_space: f64,
        /// Areal heat capacity (unit: J/m2.K)
        areal_heat_capacity: f64,
        mass_distribution_class: MassDistributionClass,
    },
}

impl BuildingElement {
    pub(crate) fn pitch(&self) -> f64 {
        *match self {
            BuildingElement::Opaque { pitch, .. } => pitch,
            BuildingElement::Transparent { pitch, .. } => pitch,
            BuildingElement::Ground { pitch, .. } => pitch,
            BuildingElement::AdjacentConditionedSpace { pitch, .. } => pitch,
            BuildingElement::AdjacentUnconditionedSpace { pitch, .. } => pitch,
        }
    }

    pub(crate) fn height(&self) -> Option<f64> {
        match self {
            BuildingElement::Opaque { height, .. } => Some(*height),
            BuildingElement::Transparent { height, .. } => Some(*height),
            _ => None,
        }
    }

    pub(crate) fn width(&self) -> Option<f64> {
        match self {
            BuildingElement::Opaque { width, .. } => Some(*width),
            BuildingElement::Transparent { width, .. } => Some(*width),
            _ => None,
        }
    }

    pub(crate) fn u_value(&self) -> Option<f64> {
        match self {
            BuildingElement::Opaque { u_value, .. } => *u_value,
            BuildingElement::Transparent { u_value, .. } => *u_value,
            BuildingElement::Ground { u_value, .. } => Some(*u_value),
            BuildingElement::AdjacentConditionedSpace { u_value, .. } => *u_value,
            BuildingElement::AdjacentUnconditionedSpace { u_value, .. } => *u_value,
        }
    }

    pub(crate) fn orientation(&self) -> Option<f64> {
        match self {
            BuildingElement::Opaque { orientation, .. } => Some(*orientation),
            BuildingElement::Transparent { orientation, .. } => Some(*orientation),
            _ => None,
        }
    }
}

pub(crate) trait TransparentBuildingElement {
    fn set_window_openable_control(&mut self, control: &str);
    fn is_security_risk(&self) -> bool;
    fn treatment(&mut self) -> Option<Vec<&mut Map<std::string::String, JsonValue>>>;
}

pub(crate) struct TransparentBuildingElementJsonValue<'a>(
    pub(crate) &'a mut Map<std::string::String, JsonValue>,
);

impl TransparentBuildingElement for TransparentBuildingElementJsonValue<'_> {
    fn set_window_openable_control(&mut self, control: &str) {
        self.0
            .insert("Control_WindowOpenable".to_string(), json!(control));
    }

    fn is_security_risk(&self) -> bool {
        self.0
            .get("security_risk")
            .and_then(|v| v.as_bool())
            .unwrap_or(false)
    }

    fn treatment(&mut self) -> Option<Vec<&mut Map<std::string::String, JsonValue>>> {
        self.0
            .get_mut("treatment")
            .and_then(|v| v.as_array_mut())
            .map(|v| v.iter_mut().flat_map(|v| v.as_object_mut()).collect())
    }
}

pub(crate) trait GroundBuildingElement {
    fn set_u_value(&mut self, new_u_value: f64);
    fn set_thermal_resistance_floor_construction(&mut self, new_r_f: f64);
    fn set_psi_wall_floor_junc(&mut self, new_psi_wall_floor_junc: f64);
}

pub(crate) struct GroundBuildingElementJsonValue<'a>(
    pub(crate) &'a mut Map<std::string::String, JsonValue>,
);

impl GroundBuildingElement for GroundBuildingElementJsonValue<'_> {
    fn set_u_value(&mut self, new_u_value: f64) {
        self.0.insert("u_value".to_string(), json!(new_u_value));
    }

    fn set_thermal_resistance_floor_construction(
        &mut self,
        new_thermal_resistance_floor_construction: f64,
    ) {
        self.0.insert(
            "thermal_resistance_floor_construction".to_string(),
            json!(new_thermal_resistance_floor_construction),
        );
    }

    fn set_psi_wall_floor_junc(&mut self, new_psi_wall_floor_junc: f64) {
        self.0.insert(
            "psi_wall_floor_junc".to_string(),
            json!(new_psi_wall_floor_junc),
        );
    }
}

pub(crate) trait UValueEditableBuildingElement {
    fn set_u_value(&mut self, new_u_value: f64);
    fn pitch(&self) -> JsonAccessResult<f64>;
    fn is_opaque(&self) -> bool;
    fn is_external_door(&self) -> Option<bool>;
    fn remove_thermal_resistance_construction(&mut self);
    fn height(&self) -> Option<f64>;
    fn width(&self) -> Option<f64>;
    fn u_value(&self) -> Option<f64>;
}

pub(crate) struct UValueEditableBuildingElementJsonValue<'a>(
    pub(crate) &'a mut Map<std::string::String, JsonValue>,
);

impl UValueEditableBuildingElement for UValueEditableBuildingElementJsonValue<'_> {
    fn set_u_value(&mut self, new_u_value: f64) {
        self.0.insert("u_value".to_string(), json!(new_u_value));
    }

    fn pitch(&self) -> JsonAccessResult<f64> {
        self.0
            .get("pitch")
            .ok_or(json_error("Pitch field not provided"))?
            .as_f64()
            .ok_or(json_error("Pitch field did not provide number"))
    }

    fn is_opaque(&self) -> bool {
        self.0
            .get("type")
            .and_then(|v| v.as_str())
            .is_some_and(|building_type| building_type == "BuildingElementOpaque")
    }

    fn is_external_door(&self) -> Option<bool> {
        self.0.get("is_external_door").and_then(|v| v.as_bool())
    }

    fn remove_thermal_resistance_construction(&mut self) {
        self.0.remove("thermal_resistance_construction");
    }

    fn height(&self) -> Option<f64> {
        self.0.get("height").and_then(|v| v.as_f64())
    }

    fn width(&self) -> Option<f64> {
        self.0.get("width").and_then(|v| v.as_f64())
    }

    fn u_value(&self) -> Option<f64> {
        self.0.get("u_value").and_then(|v| v.as_f64())
    }
}

// special deserialization logic so that orientations are normalized correctly on the way in
pub(crate) fn deserialize_orientation<'de, D>(deserializer: D) -> Result<f64, D::Error>
where
    D: Deserializer<'de>,
{
    let orientation360_value: f64 = Deserialize::deserialize(deserializer)?;
    Ok(init_orientation(orientation360_value))
}

// special serialization logic so that orientations are un-normalized correctly on way out!
pub(crate) fn serialize_orientation<S>(normalized: &f64, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let orientation360_value = init_orientation(*normalized);
    orientation360_value.serialize(serializer)
}

pub(crate) fn init_orientation(value: f64) -> f64 {
    // Convert orientation from 0-360 (clockwise) to -180 to +180 (anticlockwise)
    180. - value
}

#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum MassDistributionClass {
    D,
    E,
    I,
    IE,
    M,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct WindowPart {
    pub(crate) mid_height_air_flow_path: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct WindowTreatment {
    #[serde(rename = "type")]
    pub(crate) treatment_type: WindowTreatmentType,
    pub(crate) controls: WindowTreatmentControl,
    pub(crate) delta_r: f64,
    pub(crate) trans_red: f64,
    #[serde(
        rename = "Control_closing_irrad",
        skip_serializing_if = "Option::is_none"
    )]
    pub(crate) closing_irradiance_control: Option<String>,
    #[serde(
        rename = "Control_opening_irrad",
        skip_serializing_if = "Option::is_none"
    )]
    pub(crate) opening_irradiance_control: Option<String>,
    #[serde(rename = "Control_open", skip_serializing_if = "Option::is_none")]
    pub(crate) open_control: Option<String>,
    #[serde(
        skip_serializing_if = "Option::is_none",
        deserialize_with = "deserialize_possible_string_for_boolean",
        default
    )]
    pub(crate) is_open: Option<bool>,
    #[serde(default)]
    pub(crate) opening_delay_hrs: f64,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields, rename_all = "lowercase")]
pub enum WindowTreatmentType {
    Curtains,
    Blinds,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields, rename_all = "snake_case")]
pub enum WindowTreatmentControl {
    Manual,
    ManualMotorised,
    AutoMotorised,
    #[serde(rename = "combined_light_blind_HVAC")]
    CombinedLightBlindHvac,
}

pub(crate) fn deserialize_possible_string_for_boolean<'de, D>(
    deserializer: D,
) -> Result<Option<bool>, D::Error>
where
    D: Deserializer<'de>,
{
    let json_value: JsonValue = Deserialize::deserialize(deserializer)?;
    Ok(match json_value {
        JsonValue::Null => None,
        JsonValue::Bool(x) => Some(x),
        JsonValue::String(_) => None,
        _ => {
            return Err(Error::custom(
                "expected boolean (expected) or string (ignored)",
            ))
        }
    })
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "floor_type")]
pub(crate) enum FloorData {
    #[serde(rename = "Slab_no_edge_insulation")]
    SlabNoEdgeInsulation,
    #[serde(rename = "Slab_edge_insulation")]
    SlabEdgeInsulation {
        #[validate(min_items = 1)]
        edge_insulation: Vec<EdgeInsulation>,
    },
    // (non-optional fields in the schema for SuspendedFloor, but values do seem expected)
    #[serde(rename = "Suspended_floor")]
    SuspendedFloor {
        height_upper_surface: f64,
        /// Thermal transmittance of walls above ground (unit: W/m².K)
        #[serde(rename = "thermal_transm_walls")]
        thermal_transmission_walls: f64,
        /// Area of ventilation openings per perimeter (unit: m²/m)
        area_per_perimeter_vent: f64,
        /// Wind shielding factor
        shield_fact_location: WindShieldLocation,
        /// Thermal resistance of insulation on base of underfloor space (unit: m².K/W)
        #[serde(rename = "thermal_resist_insul")]
        thermal_resistance_of_insulation: f64,
    },
    #[serde(rename = "Heated_basement")]
    HeatedBasement {
        /// Depth of basement floor below ground level (unit: m)
        depth_basement_floor: f64,
        /// Thermal resistance of walls of the basement (unit: m².K/W)
        #[serde(rename = "thermal_resist_walls_base")]
        thermal_resistance_of_basement_walls: f64,
    },
    #[serde(rename = "Unheated_basement")]
    UnheatedBasement {
        /// Thermal transmittance of floor above basement (unit: W/m².K)
        #[serde(rename = "thermal_transm_envi_base")]
        thermal_transmittance_of_floor_above_basement: f64,
        /// Thermal transmittance of walls above ground (unit: W/m².K)
        #[serde(rename = "thermal_transm_walls")]
        thermal_transmission_walls: f64,
        /// Depth of basement floor below ground level (unit: m)
        depth_basement_floor: f64,
        /// Height of the basement walls above ground level (unit: m)
        height_basement_walls: f64,
        /// Thermal resistance of walls of the basement (unit: m².K/W)
        #[serde(rename = "thermal_resist_walls_base")]
        thermal_resistance_of_basement_walls: f64,
    },
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum WindShieldLocation {
    Sheltered,
    Average,
    Exposed,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "type")]
pub enum EdgeInsulation {
    #[serde(rename = "horizontal")]
    Horizontal {
        width: f64,
        edge_thermal_resistance: f64,
    },
    #[serde(rename = "vertical")]
    Vertical {
        depth: f64,
        edge_thermal_resistance: f64,
    },
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(untagged)]
pub enum ThermalBridging {
    Elements(IndexMap<String, ThermalBridgingDetails>),
    Number(f64),
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(tag = "type", deny_unknown_fields)]
pub enum ThermalBridgingDetails {
    #[serde(rename = "ThermalBridgeLinear")]
    Linear {
        linear_thermal_transmittance: f64,
        length: f64,
    },
    #[serde(rename = "ThermalBridgePoint")]
    Point {
        #[serde(rename = "heat_transfer_coeff")]
        heat_transfer_coefficient: f64,
    },
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum HeatingControlType {
    #[serde(rename = "SeparateTimeAndTempControl")]
    SeparateTimeAndTemperatureControl,
    #[serde(rename = "SeparateTempControl")]
    SeparateTemperatureControl,
}

pub(crate) type SpaceCoolSystem = IndexMap<String, SpaceCoolSystemDetails>;

#[derive(Clone, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(tag = "type", deny_unknown_fields)]
pub(crate) enum SpaceCoolSystemDetails {
    AirConditioning {
        cooling_capacity: f64,
        efficiency: f64,
        frac_convective: f64,
        #[serde(rename = "EnergySupply")]
        energy_supply: String,
        #[serde(rename = "Control")]
        control: String,
    },
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum WaterHeatingSchedule {
    AllDay,
    HeatingHours,
}

pub(crate) type HeatSourceWet = IndexMap<String, HeatSourceWetDetails>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[allow(clippy::large_enum_variant)]
#[serde(tag = "type", deny_unknown_fields)]
pub(crate) enum HeatSourceWetDetails {
    HeatPump {
        #[serde(rename = "EnergySupply")]
        energy_supply: String,
        source_type: HeatPumpSourceType,
        #[serde(
            rename = "EnergySupply_heat_network",
            skip_serializing_if = "Option::is_none"
        )]
        energy_supply_heat_network: Option<String>,
        #[serde(skip_serializing_if = "Option::is_none")]
        temp_distribution_heat_network: Option<f64>,
        sink_type: HeatPumpSinkType,
        #[serde(rename = "backup_ctrl_type")]
        backup_control_type: HeatPumpBackupControlType,
        #[serde(skip_serializing_if = "Option::is_none")]
        time_delay_backup: Option<f64>,
        modulating_control: bool,
        #[serde(skip_serializing_if = "Option::is_none")]
        min_modulation_rate_20: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        min_modulation_rate_35: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        min_modulation_rate_55: Option<f64>,
        time_constant_onoff_operation: f64,
        #[serde(skip_serializing_if = "Option::is_none")]
        temp_return_feed_max: Option<f64>,
        temp_lower_operating_limit: f64,
        min_temp_diff_flow_return_for_hp_to_operate: f64,
        var_flow_temp_ctrl_during_test: bool,
        #[serde(skip_serializing_if = "Option::is_none")]
        power_heating_warm_air_fan: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        power_heating_circ_pump: Option<f64>,
        power_source_circ_pump: f64,
        power_standby: f64,
        power_crankcase_heater: f64,
        power_off: f64,
        #[serde(skip_serializing_if = "Option::is_none")]
        power_max_backup: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        eahp_mixed_max_temp: Option<f64>,
        #[serde(skip_serializing_if = "Option::is_none")]
        eahp_mixed_min_temp: Option<f64>,
        #[serde(
            rename = "MechanicalVentilation",
            skip_serializing_if = "Option::is_none"
        )]
        mechanical_ventilation: Option<String>,
        #[serde(rename = "BufferTank", skip_serializing_if = "Option::is_none")]
        buffer_tank: Option<Box<HeatPumpBufferTank>>,
        #[serde(rename = "test_data_EN14825")]
        #[validate]
        test_data: Vec<HeatPumpTestDatum>,
        boiler: Option<Box<HeatPumpBoiler>>,
    },
    Boiler {
        #[serde(rename = "EnergySupply")]
        energy_supply: String,
        #[serde(rename = "EnergySupply_aux")]
        energy_supply_aux: String,
        rated_power: f64,
        efficiency_full_load: f64,
        efficiency_part_load: f64,
        boiler_location: HeatSourceLocation,
        #[validate(maximum = 1.)]
        modulation_load: f64,
        electricity_circ_pump: f64,
        electricity_part_load: f64,
        electricity_full_load: f64,
        electricity_standby: f64,
    },
    HeatBattery {
        #[serde(rename = "EnergySupply")]
        energy_supply: String,
        heat_battery_location: Option<HeatSourceLocation>,
        electricity_circ_pump: f64,
        electricity_standby: f64,
        // in kW (Charging)
        rated_charge_power: f64,
        // in kW (Losses to internal or external)
        max_rated_losses: f64,
        // number of units installed in zone
        number_of_units: usize,
        simultaneous_charging_and_discharging: bool,
        #[serde(rename = "ControlCharge")]
        control_charge: String,
        #[serde(rename = "heat_storage_zone_material_kJ_per_K_above_Phase_transition")]
        heat_storage_zone_material_k_j_per_k_above_phase_transition: f64,
        #[serde(rename = "heat_storage_zone_material_kJ_per_K_below_Phase_transition")]
        heat_storage_zone_material_k_j_per_k_below_phase_transition: f64,
        #[serde(rename = "heat_storage_zone_material_kJ_per_K_during_Phase_transition")]
        heat_storage_zone_material_k_j_per_k_during_phase_transition: f64,
        phase_transition_temperature_upper: f64,
        phase_transition_temperature_lower: f64,
        max_temperature: f64,
        #[serde(rename = "velocity_in_HEX_tube_at_1_l_per_min_m_per_s")]
        velocity_in_hex_tube_at_1_l_per_min_m_per_s: f64,
        capillary_diameter_m: f64,
        #[serde(rename = "A")]
        a: f64,
        #[serde(rename = "B")]
        b: f64,
        heat_exchanger_surface_area_m2: f64,
        flow_rate_l_per_min: f64,
    },
    #[serde(rename = "HIU")]
    Hiu {
        #[serde(rename = "EnergySupply")]
        energy_supply: String,
        power_max: f64,
        #[serde(rename = "HIU_daily_loss")]
        hiu_daily_loss: f64,
        building_level_distribution_losses: f64,
    },
}

impl From<&HeatPumpBoiler> for HeatSourceWetDetails {
    fn from(value: &HeatPumpBoiler) -> Self {
        HeatSourceWetDetails::Boiler {
            energy_supply: value.energy_supply.clone(),
            energy_supply_aux: value.energy_supply_aux.clone(),
            rated_power: value.rated_power,
            efficiency_full_load: value.efficiency_full_load,
            efficiency_part_load: value.efficiency_part_load,
            boiler_location: value.boiler_location,
            modulation_load: value.modulation_load,
            electricity_circ_pump: value.electricity_circ_pump,
            electricity_part_load: value.electricity_part_load,
            electricity_full_load: value.electricity_full_load,
            electricity_standby: value.electricity_standby,
        }
    }
}

#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum HeatPumpSourceType {
    Ground,
    OutsideAir,
    ExhaustAirMEV,
    ExhaustAirMVHR,
    ExhaustAirMixed,
    WaterGround,
    WaterSurface,
    HeatNetwork,
}

#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum HeatPumpSinkType {
    Water,
    Air,
    Glycol25,
}

#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum HeatPumpBackupControlType {
    None,
    TopUp,
    Substitute,
}

#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct HeatPumpBufferTank {
    pub daily_losses: f64,
    pub volume: f64,
    pub pump_fixed_flow_rate: f64,
    pub pump_power_at_flow_rate: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Validate, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct HeatPumpTestDatum {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) air_flow_rate: Option<f64>,
    pub(crate) test_letter: TestLetter,
    pub(crate) capacity: f64,
    pub(crate) cop: f64,
    #[serde(rename = "degradation_coeff")]
    pub(crate) degradation_coefficient: f64,
    pub(crate) design_flow_temp: f64,
    #[validate(minimum = -273.15)]
    pub(crate) temp_outlet: f64,
    #[validate(minimum = -273.15)]
    pub(crate) temp_source: f64,
    pub(crate) temp_test: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) eahp_mixed_ext_air_ratio: Option<f64>,
    // // following is not defined in Python input, but
    // #[serde(skip)]
    // pub(crate) ext_air_ratio: Option<f64>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct HeatPumpBoiler {
    #[serde(rename = "EnergySupply")]
    pub(crate) energy_supply: String,
    #[serde(rename = "EnergySupply_aux")]
    pub(crate) energy_supply_aux: String,
    rated_power: f64,
    efficiency_full_load: f64,
    efficiency_part_load: f64,
    boiler_location: HeatSourceLocation,
    modulation_load: f64,
    electricity_circ_pump: f64,
    electricity_part_load: f64,
    electricity_full_load: f64,
    electricity_standby: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate]
    pub(crate) cost_schedule_hybrid: Option<BoilerCostScheduleHybrid>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct BoilerCostScheduleHybrid {
    #[validate(minimum = 0)]
    #[validate(maximum = 365)]
    pub cost_schedule_start_day: u32,
    pub cost_schedule_time_series_step: f64,
    pub cost_schedule_hp: NumericSchedule,
    pub cost_schedule_boiler: NumericSchedule,
}

#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum HeatSourceLocation {
    #[serde(rename = "internal")]
    Internal,
    #[serde(rename = "external")]
    External,
}

pub(crate) type WasteWaterHeatRecovery = IndexMap<String, WasteWaterHeatRecoveryDetails>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct WasteWaterHeatRecoveryDetails {
    #[serde(rename = "type")]
    pub(crate) system_type: WasteWaterHeatRecoverySystemType,
    #[serde(rename = "ColdWaterSource")]
    pub(crate) cold_water_source: ColdWaterSourceType,
    pub(crate) flow_rates: Vec<f64>,
    pub(crate) efficiencies: Vec<f64>,
    pub(crate) utilisation_factor: f64,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum WasteWaterHeatRecoverySystemType {
    #[serde(rename = "WWHRS_InstantaneousSystemA")]
    SystemA,
    #[serde(rename = "WWHRS_InstantaneousSystemB")]
    SystemB,
    #[serde(rename = "WWHRS_InstantaneousSystemC")]
    SystemC,
}

pub(crate) type OnSiteGeneration = IndexMap<String, OnSiteGenerationDetails>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "type", deny_unknown_fields)]
pub(crate) enum OnSiteGenerationDetails {
    PhotovoltaicSystem {
        /// Peak power; represents the electrical power of a photovoltaic system with a given area for a solar irradiance of 1 kW/m² on this surface (at 25 degrees) (unit: kW)
        peak_power: f64,
        ventilation_strategy: PhotovoltaicVentilationStrategy,
        /// The tilt angle (inclination) of the PV panel from horizontal, measured upwards facing, 0 to 90 (unit: ˚)
        pitch: f64,
        #[serde(
            rename = "orientation360",
            deserialize_with = "deserialize_orientation",
            serialize_with = "serialize_orientation"
        )]
        orientation: f64,
        /// The distance between the ground and the lowest edge of the PV array (unit: m)
        base_height: f64,
        /// Height of the PV array (unit: m)
        height: f64,
        /// Width of the PV panel (unit: m)
        width: f64,
        #[serde(rename = "EnergySupply")]
        energy_supply: String,
        shading: Vec<WindowShadingObject>,
        inverter_peak_power_ac: f64,
        inverter_peak_power_dc: f64,
        /// Whether the inverter is considered inside the building
        inverter_is_inside: bool,
        inverter_type: InverterType,
    },
}

#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum PhotovoltaicVentilationStrategy {
    #[serde(rename = "unventilated")]
    Unventilated,
    #[serde(rename = "moderately_ventilated")]
    ModeratelyVentilated,
    #[serde(rename = "strongly_or_forced_ventilated")]
    StronglyOrForcedVentilated,
    #[serde(rename = "rear_surface_free")]
    RearSurfaceFree,
}

#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "snake_case")]
pub(crate) enum InverterType {
    OptimisedInverter,
    StringInverter,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct WindowOpeningForCooling {
    pub equivalent_area: f64,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub struct General {
    pub storeys_in_building: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub storey_of_dwelling: Option<usize>,
    pub build_type: BuildType,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum BuildType {
    #[serde(rename = "house")]
    House,
    #[serde(rename = "flat")]
    Flat,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct InfiltrationVentilation {
    pub(crate) cross_vent_possible: bool,
    /// Indicates the exposure to wind of an air flow path on a facade (can can be open, normal and shielded)
    pub(crate) shield_class: VentilationShieldClass,
    pub(crate) terrain_class: TerrainClass,
    /// Altitude of dwelling above sea level (unit: m)
    pub(crate) altitude: f64,
    /// Minimum ACH (Air Changes per Hour) limit
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = 0.)]
    pub(crate) ach_min_static_calcs: Option<f64>,
    /// Maximum ACH (Air Changes per Hour) limit
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = 0.)]
    pub(crate) ach_max_static_calcs: Option<f64>,
    #[serde(
        rename = "Control_VentAdjustMin",
        skip_serializing_if = "Option::is_none"
    )]
    pub(crate) control_vent_adjust_min: Option<String>,
    #[serde(
        rename = "Control_VentAdjustMax",
        skip_serializing_if = "Option::is_none"
    )]
    pub(crate) control_vent_adjust_max: Option<String>,
    /// Initial vent position, 0 = vents closed and 1 = vents fully open
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) vent_opening_ratio_init: Option<f64>,
    #[serde(
        rename = "Control_WindowAdjust",
        skip_serializing_if = "Option::is_none"
    )]
    pub(crate) control_window_adjust: Option<String>,
    /// Provides details about available non-mechanical ventilation systems
    #[serde(rename = "Vents")]
    pub(crate) vents: IndexMap<String, Vent>,
    /// List of the required inputs for Leaks
    #[serde(rename = "Leaks")]
    pub(crate) leaks: VentilationLeaks,
    /// Provides details about available mechanical ventilation systems
    #[serde(rename = "MechanicalVentilation", default)]
    pub(crate) mechanical_ventilation: IndexMap<String, MechanicalVentilation>,
    /// Base height of the ventilation zone relative to ground (m)
    pub(crate) ventilation_zone_base_height: f64,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum VentilationShieldClass {
    Open,
    Normal,
    Shielded,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum TerrainClass {
    OpenWater,
    OpenField,
    Suburban,
    Urban,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct Vent {
    pub(crate) mid_height_air_flow_path: f64,
    pub(crate) area_cm2: f64,
    /// Reference pressure difference for an air terminal device (unit: Pa)
    pub(crate) pressure_difference_ref: f64,
    #[serde(
        rename = "orientation360",
        deserialize_with = "deserialize_orientation",
        serialize_with = "serialize_orientation"
    )]
    pub(crate) orientation: f64,
    /// Tilt angle of the surface from horizontal, between 0 and 180, where 0 means the external surface is facing up, 90 means the external surface is vertical and 180 means the external surface is facing down (unit: ˚
    pub(crate) pitch: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct VentilationLeaks {
    pub ventilation_zone_height: f64,
    /// Reference pressure difference (unit: Pa)
    pub test_pressure: f64,
    /// Flow rate through
    pub test_result: f64,
    /// Reference area of the envelope airtightness index
    pub env_area: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct MechanicalVentilation {
    #[serde(rename = "sup_air_flw_ctrl")]
    pub(crate) supply_air_flow_rate_control: SupplyAirFlowRateControlType,
    #[serde(rename = "sup_air_temp_ctrl")]
    pub(crate) supply_air_temperature_control_type: AcceptedSupplyAirTemperatureControlType,
    pub(crate) vent_type: MechVentType,
    #[serde(rename = "mvhr_eff", skip_serializing_if = "Option::is_none")]
    pub(crate) mvhr_efficiency: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) mvhr_location: Option<MVHRLocation>,
    #[serde(rename = "Control", skip_serializing_if = "Option::is_none")]
    pub(crate) control: Option<String>,
    /// Specific fan power, inclusive of any in use factors (unit: W/l/s)
    #[serde(rename = "SFP")]
    pub(crate) sfp: f64,
    #[serde(rename = "EnergySupply")]
    pub(crate) energy_supply: String,
    /// unit: m³/hour
    pub(crate) design_outdoor_air_flow_rate: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) ductwork: Option<Vec<MechanicalVentilationDuctwork>>,
}

pub(crate) trait MechanicalVentilationForProcessing {
    fn vent_is_type(&self, vent_type: &str) -> bool;
    fn measured_fan_power(&self) -> Option<f64>;
    fn measured_air_flow_rate(&self) -> Option<f64>;
    fn set_sfp(&mut self, sfp: f64);
    fn set_control(&mut self, control: &str);
}

pub(crate) struct MechanicalVentilationJsonValue<'a>(
    pub(crate) &'a mut Map<std::string::String, JsonValue>,
);

impl MechanicalVentilationForProcessing for MechanicalVentilationJsonValue<'_> {
    fn vent_is_type(&self, vent_type: &str) -> bool {
        self.0
            .get("vent_type")
            .and_then(|value_type| value_type.as_str())
            .is_some_and(|existing_type| existing_type == vent_type)
    }

    fn measured_fan_power(&self) -> Option<f64> {
        self.0.get("measured_fan_power").and_then(|m| m.as_f64())
    }

    fn measured_air_flow_rate(&self) -> Option<f64> {
        self.0
            .get("measured_air_flow_rate")
            .and_then(|m| m.as_f64())
    }

    fn set_sfp(&mut self, sfp: f64) {
        self.0.insert("SFP".into(), json!(sfp));
    }

    fn set_control(&mut self, control: &str) {
        self.0.insert("Control".into(), json!(control));
    }
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum SupplyAirFlowRateControlType {
    #[serde(rename = "ODA")]
    Oda,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum SupplyAirTemperatureControlType {
    #[serde(rename = "CONST")]
    Constant,
    #[serde(rename = "NO_CTRL")]
    NoControl,
    #[serde(rename = "LOAD_COM")]
    LoadCom,
}

// The SupplyAirTemperatureControlType values that are currently accepted/ implemented.
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum AcceptedSupplyAirTemperatureControlType {
    #[serde(rename = "NO_CTRL")]
    NoControl,
}

impl From<AcceptedSupplyAirTemperatureControlType> for SupplyAirTemperatureControlType {
    fn from(value: AcceptedSupplyAirTemperatureControlType) -> Self {
        match value {
            AcceptedSupplyAirTemperatureControlType::NoControl => {
                SupplyAirTemperatureControlType::NoControl
            }
        }
    }
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum MechVentType {
    #[serde(rename = "Intermittent MEV")]
    IntermittentMev,
    #[serde(rename = "Centralised continuous MEV")]
    CentralisedContinuousMev,
    #[serde(rename = "Decentralised continuous MEV")]
    DecentralisedContinuousMev,
    #[serde(rename = "MVHR")]
    Mvhr,
    #[serde(rename = "PIV")]
    Piv,
}

impl Display for MechVentType {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            serde_json::to_value(self).unwrap().as_str().unwrap()
        )
    }
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct MechanicalVentilationDuctwork {
    pub(crate) cross_section_shape: DuctShape,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) duct_perimeter_mm: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) internal_diameter_mm: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) external_diameter_mm: Option<f64>,
    pub(crate) length: f64,
    pub(crate) insulation_thermal_conductivity: f64,
    pub(crate) insulation_thickness_mm: f64,
    pub(crate) reflective: bool,
    pub(crate) duct_type: DuctType,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "snake_case")]
pub enum DuctShape {
    Circular,
    Rectangular,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "snake_case")]
pub enum DuctType {
    Intake,
    Supply,
    Extract,
    Exhaust,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct AirTerminalDevice {
    area_cm2: f64,
    pressure_difference_ref: f64,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct CombustionAppliance {
    pub(crate) supply_situation: CombustionAirSupplySituation,
    pub(crate) exhaust_situation: FlueGasExhaustSituation,
    pub(crate) fuel_type: CombustionFuelType,
    pub(crate) appliance_type: CombustionApplianceType,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum CombustionAirSupplySituation {
    #[serde(rename = "room_air")]
    RoomAir,
    #[serde(rename = "outside")]
    Outside,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum FlueGasExhaustSituation {
    #[serde(rename = "into_room")]
    IntoRoom,
    #[serde(rename = "into_separate_duct")]
    IntoSeparateDuct,
    #[serde(rename = "into_mech_vent")]
    IntoMechVent,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "snake_case")]
pub enum CombustionFuelType {
    Wood,
    Gas,
    Oil,
    Coal,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "snake_case")]
pub enum CombustionApplianceType {
    OpenFireplace,
    ClosedWithFan,
    OpenGasFlueBalancer,
    OpenGasKitchenStove,
    OpenGasFire,
    ClosedFire,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum ApplianceKey {
    Fridge,
    Freezer,
    #[serde(rename = "Fridge-Freezer")]
    FridgeFreezer,
    #[serde(rename = "Otherdevices")]
    OtherDevices,
    Dishwasher,
    #[serde(rename = "Clothes_washing")]
    ClothesWashing,
    #[serde(rename = "Clothes_drying")]
    ClothesDrying,
    Oven,
    Hobs,
    Microwave,
    Kettle,
    #[serde(rename = "lighting")]
    Lighting,
}

impl Display for ApplianceKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            serde_json::to_value(self).unwrap().as_str().unwrap()
        )
    }
}

impl From<ApplianceKey> for String {
    fn from(appliance_key: ApplianceKey) -> Self {
        String::from(
            serde_json::to_value(appliance_key)
                .unwrap()
                .as_str()
                .unwrap(),
        )
    }
}

impl From<&ApplianceKey> for String {
    fn from(appliance_key: &ApplianceKey) -> Self {
        (*appliance_key).into()
    }
}

impl TryFrom<&str> for ApplianceKey {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        serde_json::from_str(format!("\"{}\"", value).as_str()).map_err(|err| anyhow!(err))
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(untagged)]
pub enum ApplianceEntry {
    Object(Appliance),
    Reference(ApplianceReference),
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub struct Appliance {
    #[serde(rename = "kWh_per_100cycle", skip_serializing_if = "Option::is_none")]
    pub(crate) kwh_per_100_cycle: Option<f64>,
    #[serde(rename = "loadshifting", skip_serializing_if = "Option::is_none")]
    pub(crate) load_shifting: Option<ApplianceLoadShifting>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) kg_load: Option<f64>,
    #[serde(rename = "kWh_per_annum", skip_serializing_if = "Option::is_none")]
    pub(crate) kwh_per_annum: Option<f64>,
    #[serde(rename = "Energysupply", skip_serializing_if = "Option::is_none")]
    pub(crate) energy_supply: Option<EnergySupplyType>,
    #[serde(rename = "kWh_per_cycle", skip_serializing_if = "Option::is_none")]
    pub(crate) kwh_per_cycle: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) standard_use: Option<f64>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub(crate) struct ApplianceLoadShifting {
    #[serde(rename = "Control", skip_serializing_if = "Option::is_none")]
    pub(crate) control: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) priority: Option<isize>,
    pub(crate) max_shift_hrs: f64,
    pub(crate) demand_limit_weighted: f64,
    pub(crate) weight_timeseries: Vec<f64>,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum ApplianceReference {
    Default,
    #[serde(rename = "Not Installed")]
    NotInstalled,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub struct Tariff {
    #[serde(rename = "schedule")]
    schedule: NumericSchedule,
}

static FHS_SCHEMA_VALIDATOR: LazyLock<Validator> = LazyLock::new(|| {
    let schema = serde_json::from_str(include_str!("../schemas/input.schema.json")).unwrap();
    jsonschema::validator_for(&schema).unwrap()
});

#[derive(Debug, Error)]
#[error("Error accessing JSON during FHS preprocessing: {0}")]
pub(crate) struct JsonAccessError(String);

pub(crate) fn json_error<T: Into<String>>(message: T) -> JsonAccessError {
    JsonAccessError(message.into())
}

pub(crate) type JsonAccessResult<T> = Result<T, JsonAccessError>;

#[derive(Clone, Debug)]
pub struct InputForProcessing {
    input: JsonValue,
}

/// This type makes methods available for restricted access by wrappers,
/// in order to work towards a reasonable API for wrappers to interact with inputs rather than
/// the more brittle approach of allowing full access to the input data structure.
/// If the full access is encapsulated within methods here, it becomes possible to update the
/// underlying structure without breaking wrappers.
impl InputForProcessing {
    pub fn init_with_json(json: impl Read) -> Result<Self, anyhow::Error> {
        let input_for_processing = Self::init_with_json_skip_validation(json)?;

        if let BasicOutput::Invalid(errors) = FHS_SCHEMA_VALIDATOR
            .apply(&input_for_processing.input)
            .basic()
        {
            bail!(
                "Invalid JSON against the FHS schema: {}",
                serde_json::to_value(errors)?.to_json_string_pretty()?
            ); // TODO build this handling logic out
        }

        Ok(input_for_processing)
    }

    pub(crate) fn init_with_json_skip_validation(json: impl Read) -> Result<Self, anyhow::Error> {
        let reader = BufReader::new(json);

        let input: JsonValue = serde_json::from_reader(reader)?;

        Ok(Self { input })
    }

    pub(crate) fn as_input(&self) -> Input {
        serde_json::from_value(self.input.to_owned()).unwrap()
    }

    pub(crate) fn finalize(self) -> Result<Input, serde_json::Error> {
        // NB. this _might_ in time be a good point to perform a validation against the core schema - or it might not
        serde_json::from_value(self.input)
    }

    fn root(&self) -> JsonAccessResult<&Map<std::string::String, JsonValue>> {
        self.input
            .as_object()
            .ok_or(json_error("Root document is not an object"))
    }

    fn root_mut(&mut self) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.input
            .as_object_mut()
            .ok_or(json_error("Root document is not an object"))
    }

    fn set_on_root_key(&mut self, root_key: &str, value: JsonValue) -> JsonAccessResult<&mut Self> {
        self.root_mut()?.insert(root_key.into(), value);

        Ok(self)
    }

    fn remove_root_key(&mut self, root_key: &str) -> JsonAccessResult<&mut Self> {
        self.root_mut()?.remove(root_key);

        Ok(self)
    }

    fn root_object(
        &self,
        root_key: &str,
    ) -> JsonAccessResult<&Map<std::string::String, JsonValue>> {
        self.root()?
            .get(root_key)
            .ok_or(json_error(format!("No {root_key} node found")))?
            .as_object()
            .ok_or(json_error(format!("{root_key} node was not an object")))
    }

    fn root_object_mut(
        &mut self,
        root_key: &str,
    ) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.root_mut()?
            .get_mut(root_key)
            .ok_or(json_error(format!("No {root_key} node found")))?
            .as_object_mut()
            .ok_or(json_error(format!("{root_key} node was not an object")))
    }

    /// Uses entry API to ensure that root key is created if it does not already exist.
    fn root_object_entry_mut(
        &mut self,
        root_key: &str,
    ) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.root_mut()?
            .entry(root_key)
            .or_insert(json!({}))
            .as_object_mut()
            .ok_or(json_error(format!("{root_key} node was not an object")))
    }

    fn optional_root_object(
        &self,
        root_key: &str,
    ) -> JsonAccessResult<Option<&Map<std::string::String, JsonValue>>> {
        Ok(self.root()?.get(root_key).and_then(|v| v.as_object()))
    }

    pub(crate) fn set_simulation_time(
        &mut self,
        simulation_time: SimulationTime,
    ) -> anyhow::Result<&mut Self> {
        self.set_on_root_key("SimulationTime", serde_json::to_value(simulation_time)?)
            .map_err(Into::into)
    }

    pub(crate) fn set_temp_internal_air_static_calcs(
        &mut self,
        temp_internal_air_static_calcs: Option<f64>,
    ) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key(
            "temp_internal_air_static_calcs",
            temp_internal_air_static_calcs.into(),
        )
    }

    pub(crate) fn merge_external_conditions_data(
        &mut self,
        external_conditions_data: Option<ExternalConditionsInput>,
    ) -> anyhow::Result<()> {
        if let Some(external_conditions) = external_conditions_data {
            let shading_segments = self
                .root_object("ExternalConditions")?
                .get("shading_segments")
                .cloned()
                .unwrap_or(json!([]));
            let mut new_external_conditions = serde_json::to_value(external_conditions)?;
            let new_external_conditions_map = new_external_conditions.as_object_mut().ok_or(json_error("External conditions was not a JSON object when it was expected to be provided as one"))?;
            new_external_conditions_map.insert("shading_segments".into(), shading_segments);
            self.set_on_root_key("ExternalConditions", new_external_conditions)?;
        }

        Ok(())
    }

    pub(crate) fn reset_internal_gains(&mut self) -> JsonAccessResult<&Self> {
        *(self.root_object_mut("InternalGains")?) = Map::new();

        Ok(self)
    }

    fn zone_node(&self) -> JsonAccessResult<&serde_json::Map<std::string::String, JsonValue>> {
        self.root_object("Zone")
    }

    fn zone_node_mut(
        &mut self,
    ) -> JsonAccessResult<&mut serde_json::Map<std::string::String, JsonValue>> {
        self.root_object_mut("Zone")
    }

    fn specific_zone(
        &self,
        zone_key: &str,
    ) -> JsonAccessResult<&serde_json::Map<std::string::String, JsonValue>> {
        self.zone_node()?
            .get(zone_key)
            .ok_or(json_error(format!("Zone key {zone_key} did not exist")))?
            .as_object()
            .ok_or(json_error("Zone node was not an object"))
    }

    fn specific_zone_mut(
        &mut self,
        zone_key: &str,
    ) -> JsonAccessResult<&mut serde_json::Map<std::string::String, JsonValue>> {
        self.zone_node_mut()?
            .get_mut(zone_key)
            .ok_or(json_error(format!("Zone key {zone_key} did not exist")))?
            .as_object_mut()
            .ok_or(json_error("Zone node was not an object"))
    }

    pub(crate) fn total_zone_area(&self) -> JsonAccessResult<f64> {
        self.zone_node()?
            .values()
            .map(|z| {
                z.get("area")
                    .ok_or(json_error("Area field not found on zone"))?
                    .as_f64()
                    .ok_or(json_error("Area field not a number"))
            })
            .sum::<JsonAccessResult<f64>>()
    }

    pub(crate) fn total_zone_volume(&self) -> JsonAccessResult<f64> {
        Ok(self
            .zone_node()?
            .values()
            .map(|z| {
                z.get("volume")
                    .ok_or(json_error("Volume field not found on zone"))?
                    .as_number()
                    .ok_or(json_error("Volume field not a number"))?
                    .as_f64()
                    .ok_or(json_error("Volume field not a number"))
            })
            .collect::<JsonAccessResult<Vec<_>>>()?
            .into_iter()
            .sum::<f64>())
    }

    pub(crate) fn area_for_zone(&self, zone: &str) -> anyhow::Result<f64> {
        Ok(self
            .zone_node()?
            .get(zone)
            .ok_or(anyhow!("Used zone key for a zone that does not exist"))?
            .get("area")
            .ok_or(json_error("Area not found on zone"))?
            .as_number()
            .ok_or(json_error("Area on zone was not a number"))?
            .as_f64()
            .ok_or(json_error("Area number could not be read as a number"))?)
    }

    #[cfg(test)]
    pub(crate) fn all_thermal_bridgings(&self) -> JsonAccessResult<Vec<&JsonValue>> {
        Ok(self
            .zone_node()?
            .values()
            .flat_map(|z| z.get("ThermalBridging"))
            .collect::<Vec<_>>())
    }

    pub(crate) fn all_thermal_bridging_elements(
        &mut self,
    ) -> JsonAccessResult<Vec<&mut Map<std::string::String, JsonValue>>> {
        let zones = self.zone_node_mut()?;
        let mut result = Vec::new();
        for zone in zones.values_mut() {
            if let Some(thermal_bridging) = zone.get_mut("ThermalBridging") {
                if let Some(obj) = thermal_bridging.as_object_mut() {
                    result.push(obj);
                }
            }
        }

        Ok(result)
    }

    pub(crate) fn number_of_bedrooms(&self) -> JsonAccessResult<Option<usize>> {
        match self.input.get("NumberOfBedrooms") {
            None => Ok(None),
            Some(JsonValue::Number(n)) => Ok(Some(
                n.as_u64()
                    .ok_or(json_error("NumberOfBedrooms not a positive integer"))?
                    as usize,
            )),
            Some(_) => Err(json_error("NumberOfBedrooms not a number")),
        }
    }

    pub(crate) fn number_of_wet_rooms(&self) -> JsonAccessResult<Option<usize>> {
        match self.input.get("NumberOfWetRooms") {
            None => Ok(None),
            Some(JsonValue::Number(n)) => Ok(Some(
                n.as_u64()
                    .ok_or(json_error("NumberOfWetRooms not a positive integer"))?
                    as usize,
            )),
            Some(_) => Err(json_error("NumberOfWetRooms not a number")),
        }
    }

    fn internal_gains_mut(&mut self) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.root_object_entry_mut("InternalGains")
    }

    pub(crate) fn set_metabolic_gains(
        &mut self,
        start_day: u32,
        time_series_step: f64,
        schedule_json: JsonValue,
    ) -> anyhow::Result<&Self> {
        self.internal_gains_mut()?.insert(
            "metabolic gains".into(),
            json!({
                "start_day": start_day,
                "time_series_step": time_series_step,
                "schedule": schedule_json,
            }),
        );

        Ok(self)
    }

    pub(crate) fn set_evaporative_losses(
        &mut self,
        start_day: u32,
        time_series_step: f64,
        schedule_json: JsonValue,
    ) -> anyhow::Result<&Self> {
        self.internal_gains_mut()?.insert(
            "EvaporativeLosses".into(),
            json!({
                "start_day": start_day,
                "time_series_step": time_series_step,
                "schedule": schedule_json,
            }),
        );

        Ok(self)
    }

    pub(crate) fn set_cold_water_losses(
        &mut self,
        start_day: u32,
        time_series_step: f64,
        schedule_json: JsonValue,
    ) -> anyhow::Result<&Self> {
        self.internal_gains_mut()?.insert(
            "ColdWaterLosses".into(),
            json!({
                "start_day": start_day,
                "time_series_step": time_series_step,
                "schedule": schedule_json,
            }),
        );

        Ok(self)
    }

    pub(crate) fn heating_control_type(&self) -> JsonAccessResult<Option<HeatingControlType>> {
        self.root()?
            .get("HeatingControlType")
            .map(
                |node| match serde_json::from_value::<HeatingControlType>(node.to_owned()) {
                    Ok(t) => Ok(t),
                    Err(_) => Err(json_error(
                        "Could not parse HeatingControlType into a known value",
                    )),
                },
            )
            .transpose()
    }

    pub(crate) fn set_heating_control_type(
        &mut self,
        heating_control_type_value: JsonValue,
    ) -> anyhow::Result<&mut Self> {
        self.set_on_root_key("HeatingControlType", heating_control_type_value)
            .map_err(Into::into)
    }

    pub(crate) fn add_control(
        &mut self,
        control_key: &str,
        control_json: JsonValue,
    ) -> JsonAccessResult<&Self> {
        self.root_object_entry_mut("Control")?
            .insert(control_key.into(), control_json);

        Ok(self)
    }

    pub(crate) fn remove_all_smart_appliance_controls(&mut self) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key("SmartApplianceControls", json!({}))
    }

    fn smart_appliance_controls_mut(
        &mut self,
    ) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.root_object_entry_mut("SmartApplianceControls")
    }

    pub(crate) fn add_smart_appliance_control(
        &mut self,
        smart_control_name: &str,
        control: JsonValue,
    ) -> JsonAccessResult<&Self> {
        self.smart_appliance_controls_mut()?
            .insert(smart_control_name.into(), control);

        Ok(self)
    }

    pub(super) fn set_non_appliance_demand_24hr_on_smart_appliance_control(
        &mut self,
        smart_control_name: &str,
        non_appliance_demand_24hr_input: IndexMap<String, Vec<f64>>,
    ) -> JsonAccessResult<&Self> {
        if let Some(ref mut control) = self
            .smart_appliance_controls_mut()?
            .get_mut(smart_control_name)
            .and_then(|v| v.as_object_mut())
        {
            control.insert(
                "non_appliance_demand_24hr".into(),
                json!(non_appliance_demand_24hr_input),
            );
        }

        Ok(self)
    }

    pub(super) fn set_battery24hr_on_smart_appliance_control(
        &mut self,
        smart_control_name: &str,
        battery24hr_input: SmartApplianceBattery,
    ) -> JsonAccessResult<&Self> {
        if let Some(ref mut control) = self
            .smart_appliance_controls_mut()?
            .get_mut(smart_control_name)
            .and_then(|v| v.as_object_mut())
        {
            control.insert("battery24hr".into(), json!(battery24hr_input));
        }

        Ok(self)
    }

    pub(crate) fn zone_keys(&self) -> JsonAccessResult<Vec<String>> {
        Ok(self.zone_node()?.keys().map(String::from).collect())
    }

    #[cfg(test)]
    pub(crate) fn all_init_temp_setpoints(&self) -> JsonAccessResult<Vec<Option<f64>>> {
        Ok(self
            .zone_node()?
            .values()
            .map(|zone| zone.get("temp_setpnt_init").and_then(|t| t.as_f64()))
            .collect())
    }

    pub(crate) fn set_init_temp_setpoint_for_zone(
        &mut self,
        zone: &str,
        temperature: f64,
    ) -> JsonAccessResult<&Self> {
        self.specific_zone_mut(zone)?
            .insert("temp_setpnt_init".into(), json!(temperature));
        Ok(self)
    }

    pub(crate) fn space_heat_control_for_zone(&self, zone: &str) -> anyhow::Result<Option<String>> {
        Ok(self
            .specific_zone(zone)?
            .get("SpaceHeatControl")
            .and_then(|field| field.as_str())
            .map(String::from))
    }

    pub(crate) fn space_heat_system_for_zone(&self, zone: &str) -> JsonAccessResult<Vec<String>> {
        Ok(match self.specific_zone(zone)?.get("SpaceHeatSystem") {
            Some(JsonValue::String(system)) => vec![String::from(system)],
            Some(JsonValue::Array(systems)) => systems
                .iter()
                .map(|system| {
                    Ok(String::from(system.as_str().ok_or(json_error(
                        "Space heat system list contained a non-string",
                    ))?))
                })
                .collect::<Result<Vec<_>, _>>()?,
            _ => vec![],
        })
    }

    pub(crate) fn set_space_heat_system_for_zone(
        &mut self,
        zone: &str,
        system_name: &str,
    ) -> anyhow::Result<&Self> {
        let zone = self.specific_zone_mut(zone)?;
        zone.insert("SpaceHeatSystem".into(), system_name.into());

        Ok(self)
    }

    pub(crate) fn space_cool_system_for_zone(&self, zone: &str) -> JsonAccessResult<Vec<String>> {
        Ok(match self.specific_zone(zone)?.get("SpaceCoolSystem") {
            Some(JsonValue::String(system)) => vec![String::from(system)],
            Some(JsonValue::Array(systems)) => systems
                .iter()
                .map(|s| {
                    Ok(String::from(s.as_str().ok_or(json_error(
                        "SpaceCoolSystem list contained non-strings",
                    ))?))
                })
                .collect::<Result<_, _>>()?,
            _ => vec![],
        })
    }

    pub(crate) fn set_space_cool_system_for_zone(
        &mut self,
        zone: &str,
        system_name: &str,
    ) -> anyhow::Result<&Self> {
        let zone = self.specific_zone_mut(zone)?;
        zone.insert("SpaceCoolSystem".into(), system_name.into());

        Ok(self)
    }

    pub(crate) fn lighting_efficacy_for_zone(&self, zone: &str) -> JsonAccessResult<Option<f64>> {
        Ok(self
            .specific_zone(zone)?
            .get("Lighting")
            .and_then(|v| v.as_object())
            .and_then(|lighting| lighting.get("efficacy"))
            .and_then(|efficacy| efficacy.as_f64()))
    }

    pub(crate) fn set_lighting_efficacy_for_all_zones(
        &mut self,
        efficacy: f64,
    ) -> JsonAccessResult<&Self> {
        for lighting in self
            .zone_node_mut()?
            .values_mut()
            .filter_map(|zone| zone.get_mut("Lighting").and_then(|l| l.as_object_mut()))
        {
            lighting.insert("efficacy".into(), efficacy.into());
        }

        Ok(self)
    }

    pub(crate) fn all_zones_have_bulbs(&self) -> JsonAccessResult<bool> {
        Ok(self.zone_node()?.values().all(|zone| {
            zone.get("Lighting")
                .and_then(|l| l.as_object())
                .and_then(|l| l.get("bulbs"))
                .is_some_and(|bulbs| bulbs.is_object())
        }))
    }

    pub(crate) fn light_bulbs_for_each_zone(
        &self,
    ) -> JsonAccessResult<IndexMap<String, Map<std::string::String, JsonValue>>> {
        Ok(self
            .zone_node()?
            .iter()
            .map(|(zone_name, zone)| {
                let bulbs = zone
                    .get("Lighting")
                    .and_then(|lighting| lighting.get("bulbs"))
                    .and_then(|bulbs| bulbs.as_object());
                (
                    String::from(zone_name),
                    bulbs.map(ToOwned::to_owned).unwrap_or_default(),
                )
            })
            .collect())
    }

    pub fn set_control_window_opening_for_zone(
        &mut self,
        zone: &str,
        opening_type: Option<&str>,
    ) -> anyhow::Result<&Self> {
        self.specific_zone_mut(zone)?
            .insert("Control_WindowOpening".into(), json!(opening_type));

        Ok(self)
    }

    pub(crate) fn set_control_string_for_space_heat_system(
        &mut self,
        space_heat_system: &str,
        control_string: &str,
    ) -> anyhow::Result<&Self> {
        self.root_object_mut("SpaceHeatSystem")?
            .get_mut(space_heat_system)
            .ok_or(anyhow!(
                "There is no provided space heat system with the name '{space_heat_system}'"
            ))?
            .as_object_mut()
            .ok_or(json_error("Space heat system was not an object"))?
            .insert("Control".into(), json!(control_string));

        Ok(self)
    }

    pub(crate) fn set_control_string_for_space_cool_system(
        &mut self,
        space_cool_system: &str,
        control_string: &str,
    ) -> anyhow::Result<&Self> {
        self.root_object_mut("SpaceCoolSystem")?
            .get_mut(space_cool_system)
            .ok_or(anyhow!(
                "There is no provided space cool system with the name '{space_cool_system}'"
            ))?
            .as_object_mut()
            .ok_or(json_error("Space cool system was not an object"))?
            .insert("Control".into(), json!(control_string));

        Ok(self)
    }

    pub(crate) fn has_named_smart_appliance_control(&self, name: &str) -> JsonAccessResult<bool> {
        Ok(self
            .optional_root_object("SmartApplianceControls")?
            .is_some_and(|controls| controls.contains_key(name)))
    }

    pub(crate) fn set_efficiency_for_all_space_cool_systems(
        &mut self,
        efficiency: f64,
    ) -> JsonAccessResult<()> {
        let systems = self.root_object_entry_mut("SpaceCoolSystem")?;
        for system in systems.values_mut().flat_map(|s| s.as_object_mut()) {
            system.insert("efficiency".into(), json!(efficiency));
        }

        Ok(())
    }

    pub(crate) fn set_frac_convective_for_all_space_cool_systems(
        &mut self,
        frac_convective: f64,
    ) -> JsonAccessResult<()> {
        let systems = self.root_object_entry_mut("SpaceCoolSystem")?;
        for system in systems.values_mut().flat_map(|s| s.as_object_mut()) {
            system.insert("frac_convective".into(), json!(frac_convective));
        }

        Ok(())
    }

    pub(crate) fn set_energy_supply_for_all_space_cool_systems(
        &mut self,
        energy_supply_name: &str,
    ) -> JsonAccessResult<()> {
        let systems = self.root_object_entry_mut("SpaceCoolSystem")?;
        for system in systems.values_mut().flat_map(|s| s.as_object_mut()) {
            system.insert("EnergySupply".into(), json!(energy_supply_name));
        }

        Ok(())
    }

    #[cfg(test)]
    pub(crate) fn space_cool_system(
        &self,
    ) -> JsonAccessResult<Option<&Map<std::string::String, JsonValue>>> {
        self.optional_root_object("SpaceCoolSystem")
    }

    pub(crate) fn space_heat_system_keys(&self) -> JsonAccessResult<Vec<String>> {
        Ok(match self.optional_root_object("SpaceHeatSystem")? {
            Some(space_heat_system) => space_heat_system.keys().map(String::from).collect(),
            None => vec![],
        })
    }

    pub(crate) fn temperature_setback_for_space_heat_system(
        &self,
        space_heat_system: &str,
    ) -> JsonAccessResult<Option<f64>> {
        let space_heat_systems = self.optional_root_object("SpaceHeatSystem")?;
        let space_heat_systems = match space_heat_systems {
            Some(ref space_heat_systems) => space_heat_systems,
            None => return Ok(None),
        };
        let space_heat_system = space_heat_systems.get(space_heat_system);

        Ok(space_heat_system.and_then(|space_heat_system| {
            space_heat_system.as_object().and_then(|space_heat_system| {
                space_heat_system
                    .get("temp_setback")
                    .and_then(|temp_setback| temp_setback.as_f64())
            })
        }))
    }

    pub(crate) fn temperature_setback_for_space_cool_system(
        &self,
        space_cool_system: &str,
    ) -> JsonAccessResult<Option<f64>> {
        let space_cool_systems = self.optional_root_object("SpaceCoolSystem")?;
        let space_cool_systems = match space_cool_systems {
            Some(ref space_heat_systems) => space_heat_systems,
            None => return Ok(None),
        };
        let space_cool_system = space_cool_systems.get(space_cool_system);

        Ok(space_cool_system.and_then(|space_cool_system| {
            space_cool_system.as_object().and_then(|space_cool_system| {
                space_cool_system
                    .get("temp_setback")
                    .and_then(|temp_setback| temp_setback.as_f64())
            })
        }))
    }

    pub(crate) fn advanced_start_for_space_cool_system(
        &self,
        space_cool_system: &str,
    ) -> JsonAccessResult<Option<f64>> {
        let space_cool_systems = self.optional_root_object("SpaceCoolSystem")?;
        let space_cool_systems = match space_cool_systems {
            Some(ref space_heat_systems) => space_heat_systems,
            None => return Ok(None),
        };
        let space_cool_system = space_cool_systems.get(space_cool_system);

        Ok(space_cool_system.and_then(|space_cool_system| {
            space_cool_system.as_object().and_then(|space_cool_system| {
                space_cool_system
                    .get("advanced_start")
                    .and_then(|temp_setback| temp_setback.as_f64())
            })
        }))
    }

    #[cfg(feature = "fhs")]
    pub(crate) fn advanced_start_for_space_heat_system(
        &self,
        space_heat_system: &str,
    ) -> JsonAccessResult<Option<f64>> {
        let space_heat_systems = self.optional_root_object("SpaceHeatSystem")?;
        let space_heat_systems = match space_heat_systems {
            Some(ref space_heat_systems) => space_heat_systems,
            None => return Ok(None),
        };
        let space_heat_system = space_heat_systems.get(space_heat_system);

        Ok(space_heat_system.and_then(|space_heat_system| {
            space_heat_system.as_object().and_then(|space_heat_system| {
                space_heat_system
                    .get("advanced_start")
                    .and_then(|temp_setback| temp_setback.as_f64())
            })
        }))
    }

    pub(crate) fn set_advance_start_for_space_heat_system(
        &mut self,
        space_heat_system: &str,
        new_advanced_start: f64,
    ) -> anyhow::Result<&Self> {
        self.root_object_entry_mut("SpaceHeatSystem")?
            .get_mut(space_heat_system)
            .ok_or(anyhow!(
                "There is no provided space heat system with the name '{space_heat_system}'"
            ))?
            .as_object_mut()
            .ok_or(json_error(
                "The indicated space heat system was not an object",
            ))?
            .insert("advanced_start".into(), new_advanced_start.into());
        Ok(self)
    }

    pub(crate) fn set_temperature_setback_for_space_heat_systems(
        &mut self,
        new_temperature_setback: Option<f64>,
    ) -> anyhow::Result<()> {
        self.root_object_entry_mut("SpaceHeatSystem")?
            .values_mut()
            .flat_map(|system| system.as_object_mut())
            .for_each(|system_details| {
                system_details.insert("temp_setback".into(), new_temperature_setback.into());
            });
        Ok(())
    }

    #[cfg(test)]
    pub(crate) fn heat_source_for_space_heat_system(
        &self,
        space_heat_system: &str,
    ) -> JsonAccessResult<Option<&JsonValue>> {
        let space_heat_systems = self.optional_root_object("SpaceHeatSystem")?;
        let space_heat_systems = match space_heat_systems {
            Some(ref space_heat_systems) => space_heat_systems,
            None => return Ok(None),
        };
        let space_heat_system = space_heat_systems.get(space_heat_system);

        Ok(space_heat_system.and_then(|space_heat_system| {
            space_heat_system
                .as_object()
                .and_then(|space_heat_system| space_heat_system.get("HeatSource"))
        }))
    }

    pub(crate) fn set_heat_source_for_all_space_heat_systems(
        &mut self,
        heat_source: SpaceHeatSystemHeatSource,
    ) -> anyhow::Result<()> {
        self.root_object_entry_mut("SpaceHeatSystem")?
            .values_mut()
            .flat_map(|system| system.as_object_mut())
            .for_each(|system_details| {
                system_details.insert("HeatSource".into(), json!(heat_source));
            });
        Ok(())
    }

    pub(crate) fn set_hot_water_source(
        &mut self,
        hot_water_source: JsonValue,
    ) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key("HotWaterSource", hot_water_source)
    }

    pub(crate) fn hot_water_source(
        &self,
    ) -> JsonAccessResult<&Map<std::string::String, JsonValue>> {
        self.root_object("HotWaterSource")
    }

    pub(crate) fn hot_water_source_mut(
        &mut self,
    ) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.root_object_mut("HotWaterSource")
    }

    pub(crate) fn names_of_energy_supplies_with_diverters(&self) -> JsonAccessResult<Vec<String>> {
        Ok(self
            .root_object("EnergySupply")?
            .iter()
            .filter_map(|(energy_supply_name, energy_supply)| {
                energy_supply.as_object().and_then(|energy_supply| {
                    energy_supply
                        .get("diverter")
                        .map(|_| String::from(energy_supply_name))
                })
            })
            .collect_vec())
    }

    pub(crate) fn set_control_max_name_for_energy_supply_diverter(
        &mut self,
        energy_supply_name: &str,
        control_max_name: &str,
    ) -> JsonAccessResult<&Self> {
        self.root_object_mut("EnergySupply")?
            .get_mut(energy_supply_name)
            .ok_or(json_error(format!(
                "There is no provided energy supply with the name '{energy_supply_name}'"
            )))?
            .as_object_mut()
            .ok_or(json_error("The indicated energy supply was not an object"))?
            .get_mut(energy_supply_name)
            .map(|energy_supply| {
                energy_supply.get("diverter").as_mut().map(|diverter| {
                    diverter.get("Controlmax").replace(&json!(control_max_name));
                })
            });
        Ok(self)
    }

    pub(crate) fn set_lighting_gains(
        &mut self,
        gains_details: JsonValue,
    ) -> JsonAccessResult<&Self> {
        self.set_gains_for_field("lighting", gains_details)
    }

    pub(crate) fn set_topup_gains(&mut self, gains_details: JsonValue) -> JsonAccessResult<&Self> {
        self.set_gains_for_field("topup", gains_details)
    }

    pub(crate) fn set_gains_for_field(
        &mut self,
        field: impl Into<std::string::String>,
        gains_details: JsonValue,
    ) -> JsonAccessResult<&Self> {
        self.root_object_entry_mut("ApplianceGains")?
            .insert(field.into(), gains_details);

        Ok(self)
    }

    pub(crate) fn energy_supply_type_for_appliance_gains_field(
        &self,
        field: &str,
    ) -> Option<String> {
        self.root_object("ApplianceGains")
            .ok()
            .and_then(|appliance_gains| appliance_gains.get(field))
            .and_then(|details| {
                details
                    .get("EnergySupply")
                    .and_then(|energy_supply| energy_supply.as_str())
                    .map(String::from)
            })
    }

    pub(crate) fn clear_appliance_gains(&mut self) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key("ApplianceGains", json!({}))
    }

    pub(crate) fn set_priority_for_gains_appliance(
        &mut self,
        priority: isize,
        appliance: &str,
    ) -> anyhow::Result<()> {
        self.root_object_entry_mut("ApplianceGains")?
            .get_mut(appliance)
            .ok_or_else(|| anyhow!("Encountered bad appliance gains reference {appliance:?}"))?
            .as_object_mut()
            .ok_or_else(|| anyhow!("Appliance gains reference was not a JSON object"))?
            .insert("priority".into(), json!(priority));

        Ok(())
    }

    pub(crate) fn fuel_type_for_energy_supply_reference(
        &self,
        reference: &str,
    ) -> anyhow::Result<String> {
        Ok(self
            .root_object("EnergySupply")?
            .get(reference)
            .ok_or(anyhow!(
                "Energy supply with reference '{reference}' could not be found"
            ))?
            .get("fuel")
            .ok_or(json_error("Energy supply object did not have a fuel field"))?
            .as_str()
            .ok_or(json_error(
                "Energy supply fuel field expected to have a string value",
            ))?
            .into())
    }

    pub(crate) fn shower_flowrates(&self) -> JsonAccessResult<IndexMap<String, f64>> {
        let showers = match self
            .hot_water_demand()?
            .get("Shower")
            .and_then(|s| s.as_object())
        {
            None => return Ok(Default::default()),
            Some(showers) => showers,
        };

        Ok(showers
            .iter()
            .filter_map(|(name, shower)| {
                shower
                    .get("flowrate")
                    .and_then(|s| s.as_f64())
                    .map(|flow_rate| (String::from(name), flow_rate))
            })
            .collect())
    }

    pub(crate) fn reset_water_heating_events(&mut self) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key("Events", json!({}))
    }

    pub(crate) fn showers(&self) -> JsonAccessResult<Option<&Map<std::string::String, JsonValue>>> {
        self.hot_water_demand()?
            .get("Shower")
            .map(|showers| {
                showers
                    .as_object()
                    .ok_or(json_error("Shower was not an object"))
            })
            .transpose()
    }

    pub(crate) fn shower_keys(&self) -> JsonAccessResult<Vec<String>> {
        Ok(self
            .hot_water_demand()?
            .get("Shower")
            .map_or(vec![], |showers| match showers.as_object() {
                Some(showers) => showers.keys().map(String::from).collect(),
                None => vec![],
            }))
    }

    pub(crate) fn shower_name_refers_to_instant_electric(&self, name: &str) -> bool {
        self.hot_water_demand()
            .ok()
            .and_then(|demand| demand.get("Shower"))
            .and_then(|showers| showers.get(name))
            .and_then(|shower| shower.get("type"))
            .and_then(|shower_type| shower_type.as_str())
            .is_some_and(|shower_type| shower_type == "InstantElecShower")
    }

    pub(crate) fn register_wwhrs_name_on_mixer_shower(
        &mut self,
        wwhrs: &str,
    ) -> anyhow::Result<()> {
        let mixer_shower = self
            .hot_water_demand_mut()?
            .get_mut("Shower")
            .ok_or(json_error("Shower node not set on HotWaterDemand"))?
            .get_mut("mixer")
            .ok_or(json_error(
                "A mixer shower with the reference 'mixer' was expected to have been set",
            ))?
            .as_object_mut()
            .ok_or(json_error("Mixer shower was not a JSON object"))?;
        mixer_shower.insert("WWHRS".into(), json!(wwhrs));

        Ok(())
    }

    pub(crate) fn baths(&self) -> JsonAccessResult<Option<&Map<std::string::String, JsonValue>>> {
        Ok(self
            .hot_water_demand()?
            .get("Bath")
            .and_then(|baths| baths.as_object()))
    }

    pub(crate) fn bath_keys(&self) -> JsonAccessResult<Vec<String>> {
        Ok(self
            .hot_water_demand()?
            .get("Bath")
            .map_or(vec![], |baths| match baths.as_object() {
                Some(baths) => baths.keys().map(String::from).collect(),
                None => vec![],
            }))
    }

    pub(crate) fn size_for_bath_field(&self, field: &str) -> JsonAccessResult<Option<f64>> {
        Ok(self
            .hot_water_demand()?
            .get("Bath")
            .and_then(|baths| baths.as_object())
            .and_then(|bath| bath.get(field))
            .and_then(|bath| bath.get("size"))
            .and_then(|size| size.as_f64()))
    }

    pub(crate) fn flowrate_for_bath_field(&self, field: &str) -> JsonAccessResult<Option<f64>> {
        Ok(self
            .hot_water_demand()?
            .get("Bath")
            .and_then(|baths| baths.as_object())
            .and_then(|bath| bath.get(field))
            .and_then(|bath| bath.get("flowrate"))
            .and_then(|flowrate| flowrate.as_f64()))
    }

    pub(crate) fn other_water_uses(
        &self,
    ) -> JsonAccessResult<Option<&Map<std::string::String, JsonValue>>> {
        Ok(self
            .hot_water_demand()?
            .get("Other")
            .and_then(|other| other.as_object()))
    }

    pub(crate) fn other_water_use_keys(&self) -> JsonAccessResult<Vec<String>> {
        Ok(self
            .other_water_uses()?
            .map(|other| other.keys().map(String::from).collect())
            .unwrap_or(vec![]))
    }

    pub(crate) fn flow_rate_for_other_water_use_field(
        &self,
        field: &str,
    ) -> JsonAccessResult<Option<f64>> {
        Ok(self
            .hot_water_demand()?
            .get("Other")
            .and_then(|others| others.as_object())
            .and_then(|other| other.get(field))
            .and_then(|other| other.get("flowrate"))
            .and_then(|flowrate| flowrate.as_f64()))
    }

    pub(crate) fn set_other_water_use_details(
        &mut self,
        cold_water_source_type: &str,
        flowrate: f64,
    ) -> JsonAccessResult<()> {
        let other_details = json!({
            "flowrate": flowrate,
            "ColdWaterSource": cold_water_source_type,
        });

        let other_water_uses = self
            .hot_water_demand_mut()?
            .entry("Other")
            .or_insert(json!({}))
            .as_object_mut()
            .ok_or(json_error("Other water uses not provided as an object"))?;
        other_water_uses.insert("other".into(), other_details);

        Ok(())
    }

    pub(crate) fn water_distribution(&self) -> anyhow::Result<Option<WaterDistribution>> {
        Ok(self
            .root_object("HotWaterDemand")?
            .get("Distribution")
            .map(|node| serde_json::from_value::<WaterDistribution>(node.to_owned()))
            .transpose()?)
    }

    /// Override all the vol_hw_daily_average values on the heat pump hot water only heat sources.
    pub(crate) fn override_vol_hw_daily_average_on_heat_pumps(
        &mut self,
        vol_hw_daily_average: f64,
    ) {
        let heat_sources = match self
            .root_object_mut("HotWaterSource")
            .ok()
            .and_then(|hot_water_source| hot_water_source.get_mut("hw cylinder"))
            .and_then(|cylinder| cylinder.as_object_mut())
            .and_then(|cylinder| {
                cylinder
                    .get("type")
                    .and_then(|source_type| source_type.as_str())
                    .is_some_and(|source_type| source_type == "StorageTank")
                    .then_some(cylinder)
            }) {
            None => return,
            Some(heat_sources) => heat_sources,
        };

        for heat_source in heat_sources.values_mut().flat_map(|hs| hs.as_object_mut()) {
            if heat_source
                .get("type")
                .and_then(|heat_source_type| heat_source_type.as_str())
                .is_some_and(|heat_source_type| heat_source_type == "HeatPump_HWOnly")
            {
                heat_source.insert("vol_hw_daily_average".into(), json!(vol_hw_daily_average));
            }
        }
    }

    pub(crate) fn part_g_compliance(&self) -> JsonAccessResult<Option<bool>> {
        self.root()?
            .get("PartGcompliance")
            .map(|node| {
                node.as_bool()
                    .ok_or(json_error("Part G compliance was not passed as a boolean"))
            })
            .transpose()
    }

    pub(crate) fn set_part_g_compliance(&mut self, is_compliant: bool) -> JsonAccessResult<&Self> {
        self.root_mut()?
            .insert("PartGcompliance".into(), is_compliant.into());

        Ok(self)
    }

    pub(crate) fn add_water_heating_event(
        &mut self,
        event_type: &str,
        subtype_name: &str,
        event: JsonValue,
    ) -> JsonAccessResult<&Self> {
        let node_for_type = self
            .root_object_entry_mut("Events")?
            .entry(event_type)
            .or_insert(json!({}))
            .as_object_mut()
            .ok_or(json_error("Events node was not an object"))?;
        let node_for_subtype = node_for_type.entry(subtype_name).or_insert(json!([]));
        if let Some(events) = node_for_subtype.as_array_mut() {
            events.push(event);
        } else {
            return Err(json_error(format!(
                "Events node at '{event_type}' -> '{subtype_name}' was not an array"
            )));
        }

        Ok(self)
    }

    pub(crate) fn water_heating_events_of_types(
        &self,
        event_types: &[&str],
    ) -> JsonAccessResult<Vec<JsonValue>> {
        Ok(self
            .root_object("Events")?
            .iter()
            .filter(|(event_type, _)| event_types.contains(&&***event_type))
            .flat_map(|(_, events)| {
                events.as_object().map(|events| {
                    events
                        .values()
                        .into_iter()
                        .map(JsonValue::as_array)
                        .flatten()
                })
            })
            .flatten()
            .flatten()
            .cloned()
            .collect_vec())
    }

    pub(crate) fn cold_water_source_has_header_tank(&self) -> JsonAccessResult<bool> {
        Ok(self
            .root_object("ColdWaterSource")?
            .contains_key("header tank"))
    }

    pub(crate) fn set_cold_water_source_by_key(
        &mut self,
        key: &str,
        source_details: JsonValue,
    ) -> JsonAccessResult<&Self> {
        self.root_object_mut("ColdWaterSource")?
            .insert(key.into(), source_details);

        Ok(self)
    }

    pub(crate) fn set_hot_water_cylinder(
        &mut self,
        source_value: JsonValue,
    ) -> JsonAccessResult<&Self> {
        let hot_water_source = self.root_object_mut("HotWaterSource")?;
        hot_water_source.insert("hw cylinder".into(), source_value);

        Ok(self)
    }

    fn hot_water_demand(&self) -> JsonAccessResult<&Map<std::string::String, JsonValue>> {
        self.root_object("HotWaterDemand")
    }

    fn hot_water_demand_mut(
        &mut self,
    ) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.root_object_mut("HotWaterDemand")
    }

    pub(crate) fn set_water_distribution(
        &mut self,
        distribution_value: JsonValue,
    ) -> JsonAccessResult<&Self> {
        self.hot_water_demand_mut()?
            .insert("Distribution".into(), distribution_value);

        Ok(self)
    }

    pub(crate) fn set_shower(&mut self, shower_value: JsonValue) -> anyhow::Result<&Self> {
        self.hot_water_demand_mut()?
            .insert("Shower".into(), shower_value);

        Ok(self)
    }

    pub(crate) fn set_bath(&mut self, bath_value: JsonValue) -> anyhow::Result<&Self> {
        self.hot_water_demand_mut()?
            .insert("Bath".into(), bath_value);

        Ok(self)
    }

    pub(crate) fn set_other_water_use(
        &mut self,
        other_water_use_value: JsonValue,
    ) -> anyhow::Result<&Self> {
        self.hot_water_demand_mut()?
            .insert("Other".into(), other_water_use_value);

        Ok(self)
    }

    pub(crate) fn remove_wwhrs(&mut self) -> JsonAccessResult<&mut Self> {
        self.remove_root_key("WWHRS")
    }

    pub(crate) fn wwhrs(&self) -> anyhow::Result<Option<WasteWaterHeatRecovery>> {
        Ok(self
            .root()?
            .get("WWHRS")
            .map(|wwhrs| match serde_json::from_value(wwhrs.to_owned()) {
                Ok(wwhrs) => Ok(wwhrs),
                Err(err) => Err(err),
            })
            .transpose()?)
    }

    pub(crate) fn set_wwhrs(&mut self, wwhrs: JsonValue) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key("WWHRS", wwhrs)
    }

    pub(crate) fn remove_space_heat_systems(&mut self) -> JsonAccessResult<&mut Self> {
        self.remove_root_key("SpaceHeatSystem")
    }

    #[cfg(test)]
    pub(crate) fn space_heat_system_for_key(
        &self,
        key: &str,
    ) -> JsonAccessResult<Option<&JsonValue>> {
        Ok(self.root_object("SpaceHeatSystem")?.get(key))
    }

    pub(crate) fn set_space_heat_system_for_key(
        &mut self,
        key: &str,
        space_heat_system_value: JsonValue,
    ) -> JsonAccessResult<&Self> {
        self.root_object_entry_mut("SpaceHeatSystem")?
            .insert(key.into(), space_heat_system_value);
        Ok(self)
    }

    pub(crate) fn remove_space_cool_systems(&mut self) -> JsonAccessResult<&mut Self> {
        self.remove_root_key("SpaceCoolSystem")
    }

    pub(crate) fn set_space_cool_system_for_key(
        &mut self,
        key: &str,
        space_cool_system_value: JsonValue,
    ) -> JsonAccessResult<&Self> {
        self.root_object_entry_mut("SpaceCoolSystem")?
            .insert(key.into(), space_cool_system_value);
        Ok(self)
    }

    pub(crate) fn set_on_site_generation(
        &mut self,
        on_site_generation: JsonValue,
    ) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key("OnSiteGeneration", on_site_generation)
    }

    pub(crate) fn remove_on_site_generation(&mut self) -> JsonAccessResult<&mut Self> {
        self.remove_root_key("OnSiteGeneration")
    }

    pub(crate) fn on_site_generation(
        &self,
    ) -> JsonAccessResult<Option<&Map<std::string::String, JsonValue>>> {
        self.optional_root_object("OnSiteGeneration")
    }

    pub(crate) fn remove_all_diverters_from_energy_supplies(
        &mut self,
    ) -> JsonAccessResult<&mut Self> {
        self.root_object_entry_mut("EnergySupply")?
            .values_mut()
            .filter_map(|value| value.as_object_mut())
            .for_each(|energy_supply| {
                energy_supply.remove("diverter");
            });
        Ok(self)
    }

    pub(crate) fn add_energy_supply_for_key(
        &mut self,
        energy_supply_key: &str,
        energy_supply_details: JsonValue,
    ) -> JsonAccessResult<()> {
        self.root_object_entry_mut("EnergySupply")?
            .insert(energy_supply_key.into(), energy_supply_details);

        Ok(())
    }

    #[cfg(test)]
    pub(crate) fn energy_supply_by_key(
        &self,
        energy_supply_key: &str,
    ) -> JsonAccessResult<Option<&Map<std::string::String, JsonValue>>> {
        Ok(self
            .root_object("EnergySupply")?
            .get(energy_supply_key)
            .and_then(|energy_supply| energy_supply.as_object()))
    }

    #[cfg(test)]
    pub(crate) fn add_diverter_to_energy_supply(
        &mut self,
        energy_supply_key: &str,
        diverter: JsonValue,
    ) -> JsonAccessResult<()> {
        if let Some(energy_supply) = self
            .root_object_mut("EnergySupply")?
            .get_mut(energy_supply_key)
            .and_then(|energy_supply| energy_supply.as_object_mut())
        {
            energy_supply.insert("diverter".into(), diverter);
        }
        Ok(())
    }

    #[cfg(test)]
    pub(crate) fn add_electric_battery_to_energy_supply(
        &mut self,
        energy_supply_key: &str,
        electric_battery: JsonValue,
    ) -> JsonAccessResult<()> {
        if let Some(energy_supply) = self
            .root_object_mut("EnergySupply")?
            .get_mut(energy_supply_key)
            .and_then(|energy_supply| energy_supply.as_object_mut())
        {
            energy_supply.insert("ElectricBattery".into(), electric_battery);
        }

        Ok(())
    }

    pub(crate) fn remove_all_batteries_from_energy_supplies(
        &mut self,
    ) -> JsonAccessResult<&mut Self> {
        for energy_supply in self
            .root_object_mut("EnergySupply")?
            .values_mut()
            .flat_map(|v| v.as_object_mut())
        {
            energy_supply.remove("ElectricBattery");
        }
        Ok(self)
    }

    pub(crate) fn external_conditions(&self) -> anyhow::Result<ExternalConditionsInput> {
        serde_json::from_value(
            self.root()?
                .get("ExternalConditions")
                .ok_or(json_error("ExternalConditions not found"))?
                .to_owned(),
        )
        .map_err(Into::into)
    }

    fn all_building_elements_mut_of_types(
        &mut self,
        types: &[&str],
    ) -> JsonAccessResult<Vec<&mut Map<std::string::String, JsonValue>>> {
        Ok(self
            .zone_node_mut()?
            .values_mut()
            .filter_map(|zone| {
                zone.get_mut("BuildingElement")
                    .and_then(|building_element_node| building_element_node.as_object_mut())
            })
            .flat_map(|building_elements| building_elements.values_mut())
            .filter(|building_element| {
                building_element
                    .get("type")
                    .and_then(|building_element_type| building_element_type.as_str())
                    .is_some_and(|building_element_type| types.contains(&building_element_type))
            })
            .filter_map(|element| element.as_object_mut())
            .collect())
    }

    pub(crate) fn all_transparent_building_elements_mut(
        &mut self,
    ) -> JsonAccessResult<Vec<&mut Map<std::string::String, JsonValue>>> {
        self.all_building_elements_mut_of_types(&["BuildingElementTransparent"])
    }

    pub(crate) fn all_ground_building_elements_mut(
        &mut self,
    ) -> JsonAccessResult<Vec<&mut Map<std::string::String, JsonValue>>> {
        self.all_building_elements_mut_of_types(&["BuildingElementGround"])
    }

    pub(crate) fn all_opaque_and_adjztu_building_elements_mut_u_values(
        &mut self,
    ) -> JsonAccessResult<Vec<&mut Map<std::string::String, JsonValue>>> {
        self.all_building_elements_mut_of_types(&[
            "BuildingElementOpaque",
            "BuildingElementAdjacentUnconditionedSpace_Simple",
        ])
    }

    pub(crate) fn max_base_height_from_building_elements(&self) -> JsonAccessResult<Option<f64>> {
        Ok(self
            .zone_node()?
            .values()
            .filter_map(|zone| zone.get("BuildingElement"))
            .filter_map(|building_element| building_element.as_object())
            .flat_map(|building_element_node| building_element_node.values())
            .filter_map(|building_element| {
                building_element.get("base_height").and_then(|h| h.as_f64())
            })
            .max_by(|a, b| a.total_cmp(b)))
    }

    pub(crate) fn set_numeric_field_for_building_element(
        &mut self,
        building_element_reference: &str,
        field: &str,
        value: f64,
    ) -> anyhow::Result<()> {
        *self.zone_node_mut()?
            .values_mut()
            .filter_map(|zone| zone.get_mut("BuildingElement").and_then(|el| el.as_object_mut()))
            .flatten()
            .find(|(name, _value)| *name == building_element_reference)
            .ok_or(anyhow!("Could not find building element with reference '{building_element_reference}'"))?
            .1
            .get_mut(field)
            .ok_or(anyhow!("Could not find field '{field}' on building element with reference '{building_element_reference}'"))? = json!(value);

        Ok(())
    }

    pub(crate) fn all_building_elements(
        &self,
    ) -> anyhow::Result<IndexMap<String, BuildingElement>> {
        self.zone_node()?
            .values()
            .filter_map(|zone| zone.get("BuildingElement").and_then(|el| el.as_object()))
            .flatten()
            .map(|(name, el)| Ok((String::from(name), serde_json::from_value(el.to_owned())?)))
            .collect()
    }

    #[cfg(test)]
    pub(crate) fn building_element_by_key(
        &self,
        zone_key: &str,
        key: &str,
    ) -> JsonAccessResult<&Map<std::string::String, JsonValue>> {
        self.specific_zone(zone_key)?
            .get("BuildingElement")
            .ok_or(json_error("BuildingElement node not present"))?
            .as_object()
            .ok_or(json_error("BuildingElement node was not an object"))?
            .get(key)
            .ok_or(json_error(format!(
                "BuildingElement with name {key} was not present"
            )))?
            .as_object()
            .ok_or(json_error(
                "Building element with name {key} not provided as an object",
            ))
    }

    pub(crate) fn all_energy_supply_fuel_types(&self) -> JsonAccessResult<HashSet<String>> {
        let mut fuel_types = HashSet::new();
        for fuel in self
            .root_object("EnergySupply")?
            .values()
            .flat_map(|supply| supply.get("fuel").map(|fuel| fuel.as_str()))
            .flatten()
        {
            fuel_types.insert(String::from(fuel));
        }

        Ok(fuel_types)
    }

    pub(crate) fn has_appliances(&self) -> JsonAccessResult<bool> {
        Ok(self.root()?.contains_key("Appliances"))
    }

    pub(crate) fn merge_in_appliances(
        &mut self,
        appliances: &IndexMap<&str, JsonValue>,
    ) -> anyhow::Result<()> {
        let mut appliances_value = serde_json::to_value(appliances.to_owned())?;
        let appliances = appliances_value
            .as_object_mut()
            .ok_or(anyhow!("Appliances were not an object when expected to be"))?;
        let existing_appliances = self.root_object_entry_mut("Appliances")?;
        existing_appliances.append(appliances);

        Ok(())
    }

    pub(crate) fn remove_appliance(&mut self, appliance_key: &str) -> JsonAccessResult<&Self> {
        self.root_object_entry_mut("Appliances")?
            .remove(appliance_key);

        Ok(self)
    }

    pub(crate) fn appliances_contain_key(&self, name: &str) -> bool {
        self.root_object("Appliances")
            .ok()
            .is_some_and(|appliances| appliances.contains_key(name))
    }

    pub(crate) fn appliance_key_has_reference(
        &self,
        key: &str,
        reference: &str,
    ) -> JsonAccessResult<bool> {
        let empty_map = Map::new();
        Ok(self
            .root_object("Appliances")
            .unwrap_or(&empty_map)
            .get(key)
            .and_then(|value| value.as_str())
            .is_some_and(|appliance_reference| appliance_reference == reference))
    }

    pub(crate) fn appliance_keys(&self) -> JsonAccessResult<Vec<String>> {
        let empty_map = Map::new();
        Ok(self
            .root_object("Appliances")
            .unwrap_or(&empty_map)
            .keys()
            .map(String::from)
            .collect())
    }

    pub(crate) fn appliance_with_key(&self, key: &str) -> JsonAccessResult<Option<&JsonValue>> {
        Ok(match self.root_object("Appliances") {
            Err(_) => return Ok(None),
            Ok(appliances) => appliances.get(key),
        })
    }

    pub(crate) fn appliance_with_key_mut(
        &mut self,
        key: &str,
    ) -> JsonAccessResult<Option<&mut JsonValue>> {
        Ok(self.root_object_entry_mut("Appliances")?.get_mut(key))
    }

    pub(crate) fn clone_appliances(&self) -> Map<std::string::String, JsonValue> {
        self.root_object("Appliances")
            .cloned()
            .unwrap_or(Map::new())
    }

    pub(crate) fn tariff_schedule(&self) -> anyhow::Result<Option<NumericSchedule>> {
        self.root_object("Tariff")
            .ok()
            .cloned()
            .and_then(|tariff| tariff.get("schedule").cloned())
            .map(|schedule| serde_json::from_value(schedule.clone()))
            .transpose()
            .map_err(|err| anyhow!(err))
    }

    pub(crate) fn energy_supply_for_appliance(&self, key: &str) -> anyhow::Result<&str> {
        let appliances = self.root_object("Appliances")?;

        appliances
            .get(key)
            .and_then(|appliance| appliance.as_object())
            .ok_or_else(|| anyhow!("No {key} object in appliances input"))?
            .get("Energysupply")
            .and_then(|supply| supply.as_str())
            .ok_or_else(|| anyhow!("No energy supply for appliance '{key}'"))
    }

    pub(crate) fn loadshifting_for_appliance(
        &self,
        appliance_key: &str,
    ) -> JsonAccessResult<Option<Map<std::string::String, JsonValue>>> {
        let appliance = self.appliance_with_key(appliance_key)?;

        Ok(appliance
            .and_then(|appliance| appliance.get("loadshifting"))
            .and_then(|load_shifting| load_shifting.as_object())
            .cloned())
    }

    pub(crate) fn set_loadshifting_for_appliance(
        &mut self,
        appliance_key: &str,
        new_load_shifting: JsonValue,
    ) -> JsonAccessResult<()> {
        let mut appliance = self.appliance_with_key_mut(appliance_key)?;
        if let Some(appliance) = appliance
            .as_mut()
            .and_then(|appliance| appliance.as_object_mut())
        {
            appliance.insert("loadshifting".into(), new_load_shifting);
        }

        Ok(())
    }

    fn infiltration_ventilation_node_mut(
        &mut self,
    ) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.root_object_mut("InfiltrationVentilation")
    }

    pub(crate) fn mechanical_ventilations_for_processing(
        &mut self,
    ) -> JsonAccessResult<Vec<&mut Map<std::string::String, JsonValue>>> {
        let mech_vents = match self
            .infiltration_ventilation_node_mut()?
            .get_mut("MechanicalVentilation")
            .and_then(|v| v.as_object_mut())
        {
            None => return Ok(Vec::new()),
            Some(mech_vents) => mech_vents,
        };
        Ok(mech_vents
            .values_mut()
            .filter_map(|v| v.as_object_mut())
            .collect())
    }

    pub(crate) fn keyed_mechanical_ventilations_for_processing(
        &mut self,
    ) -> JsonAccessResult<IndexMap<String, &mut Map<std::string::String, JsonValue>>> {
        let mech_vents = match self
            .infiltration_ventilation_node_mut()?
            .get_mut("MechanicalVentilation")
            .and_then(|v| v.as_object_mut())
        {
            None => return Ok(Default::default()),
            Some(mech_vents) => mech_vents,
        };
        Ok(mech_vents
            .iter_mut()
            .filter_map(|(name, v)| {
                let mech_vent = match v.as_object_mut() {
                    None => return None,
                    Some(mech_vent) => mech_vent,
                };
                Some((String::from(name), mech_vent))
            })
            .collect())
    }

    pub(crate) fn has_mechanical_ventilation(&self) -> bool {
        self.root_object("InfiltrationVentilation")
            .ok()
            .is_some_and(|node| node.contains_key("MechanicalVentilation"))
    }

    pub(crate) fn reset_mechanical_ventilation(&mut self) -> JsonAccessResult<&Self> {
        self.root_object_entry_mut("InfiltrationVentilation")?
            .remove("MechanicalVentilation");
        Ok(self)
    }

    pub(crate) fn add_mechanical_ventilation(
        &mut self,
        vent_name: &str,
        mech_vent: JsonValue,
    ) -> anyhow::Result<()> {
        let infiltration_ventilation_node = self
            .input
            .get_mut("InfiltrationVentilation")
            .ok_or(json_error("InfiltrationVentilation node not found"))?
            .as_object_mut()
            .ok_or(json_error("InfiltrationVentilation node is not an object"))?;
        let mech_vent_map = infiltration_ventilation_node
            .entry("MechanicalVentilation")
            .or_insert(json!({}))
            .as_object_mut()
            .ok_or(json_error("MechanicalVentilation node is not an object"))?;
        mech_vent_map.insert(vent_name.into(), mech_vent);

        Ok(())
    }

    pub(crate) fn appliance_gains_events(
        &self,
    ) -> anyhow::Result<IndexMap<String, Vec<ApplianceGainsEvent>>> {
        let appliance_gains = match self.root_object("ApplianceGains") {
            Ok(appliance_gains) => appliance_gains,
            Err(_) => return Ok(IndexMap::new()),
        };
        appliance_gains
            .iter()
            .map(
                |(name, gain)| -> Result<(String, Vec<ApplianceGainsEvent>), _> {
                    Ok((
                        String::from(name),
                        serde_json::from_value(
                            gain.get("Events")
                                .and_then(|events| events.is_array().then_some(events))
                                .cloned()
                                .unwrap_or(json!([])),
                        )?,
                    ))
                },
            )
            .collect::<anyhow::Result<_>>()
    }

    pub(crate) fn set_window_adjust_control_for_infiltration_ventilation(
        &mut self,
        control: &str,
    ) -> JsonAccessResult<&Self> {
        self.infiltration_ventilation_node_mut()?
            .insert("Control_WindowAdjust".into(), control.into());
        Ok(self)
    }

    pub(crate) fn set_vent_adjust_min_control_for_infiltration_ventilation(
        &mut self,
        control: &str,
    ) -> JsonAccessResult<&Self> {
        self.infiltration_ventilation_node_mut()?
            .insert("Control_VentAdjustMin".into(), control.into());
        Ok(self)
    }

    pub(crate) fn set_vent_adjust_max_control_for_infiltration_ventilation(
        &mut self,
        control: &str,
    ) -> JsonAccessResult<&Self> {
        self.infiltration_ventilation_node_mut()?
            .insert("Control_VentAdjustMax".into(), control.into());
        Ok(self)
    }

    pub(crate) fn infiltration_ventilation_is_noise_nuisance(&self) -> bool {
        self.root_object("InfiltrationVentilation")
            .ok()
            .and_then(|infiltration| infiltration.get("noise_nuisance"))
            .and_then(|nuisance| nuisance.as_bool())
            .unwrap_or(false)
    }

    pub(crate) fn infiltration_ventilation_mut(
        &mut self,
    ) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.root_object_entry_mut("InfiltrationVentilation")
    }

    pub(crate) fn set_heat_source_wet(
        &mut self,
        heat_source_wet: JsonValue,
    ) -> JsonAccessResult<()> {
        self.root_mut()?
            .insert("HeatSourceWet".into(), heat_source_wet);
        Ok(())
    }

    pub(crate) fn heat_source_wet(&self) -> anyhow::Result<IndexMap<String, HeatSourceWetDetails>> {
        self.root()?
            .get("HeatSourceWet")
            .and_then(|value| value.as_object())
            .into_iter()
            .flatten()
            .map(|(name, source)| Ok((String::from(name), serde_json::from_value(source.clone())?)))
            .collect::<anyhow::Result<_, _>>()
    }

    pub(crate) fn cold_water_source(&self) -> anyhow::Result<ColdWaterSourceInput> {
        Ok(serde_json::from_value(
            self.root()?
                .get("ColdWaterSource")
                .cloned()
                .ok_or(json_error("ColdWaterSource was not present"))?,
        )?)
    }

    #[cfg(test)]
    pub(crate) fn set_storeys_in_building(&mut self, storeys: usize) -> JsonAccessResult<&Self> {
        self.root_object_mut("General")?
            .insert("storeys_in_building".into(), json!(storeys));

        Ok(self)
    }

    pub(crate) fn storeys_in_building(&self) -> JsonAccessResult<usize> {
        Ok(self
            .input
            .get("General")
            .ok_or(json_error("General node not found"))?
            .get("storeys_in_building")
            .ok_or(json_error("storeys_in_building field not found"))?
            .as_u64()
            .ok_or(json_error(
                "storeys_in_building field is not a positive integer",
            ))? as usize)
    }

    pub(crate) fn build_type(&self) -> JsonAccessResult<String> {
        Ok(self
            .root_object("General")?
            .get("build_type")
            .ok_or(json_error(
                "There was no build_type field on the General input object",
            ))?
            .as_str()
            .ok_or(json_error("The build_type field was not a string"))?
            .into())
    }

    pub(crate) fn hot_water_cylinder_volume(&self) -> JsonAccessResult<Option<f64>> {
        Ok(self
            .root_object("HotWaterSource")?
            .get("hw cylinder")
            .and_then(|cylinder| cylinder.get("volume"))
            .and_then(|v| v.as_f64()))
    }

    pub(crate) fn ground_floor_area(&self) -> JsonAccessResult<Option<f64>> {
        Ok(self
            .root()?
            .get("GroundFloorArea")
            .and_then(|area| area.as_f64()))
    }

    pub(crate) fn primary_pipework_clone(&self) -> anyhow::Result<Option<Vec<WaterPipework>>> {
        Ok(self
            .hot_water_source()?
            .get("hw cylinder")
            .and_then(|cylinder| cylinder.get("primary_pipework"))
            .and_then(|primary_pipework| primary_pipework.is_array().then_some(primary_pipework))
            .map(|primary_pipework| serde_json::from_value(primary_pipework.to_owned()))
            .transpose()?)
    }

    pub(crate) fn water_heating_event_by_type_and_name(
        &self,
        event_type: &str,
        event_name: &str,
    ) -> anyhow::Result<Option<Vec<WaterHeatingEvent>>> {
        let result = Ok(self
            .root_object("Events")?
            .get(event_type)
            .and_then(|event_group| event_group.get(event_name))
            .map(|events| serde_json::from_value(events.to_owned()))
            .transpose()?);

        if let Ok(None) = result {
            println!(
                "No events found for event type {} and name {}",
                event_type, event_name
            );
            println!("Events: {:?}", self.root_object("Events")?.keys());
        }

        result
    }

    pub(crate) fn part_o_active_cooling_required(&self) -> JsonAccessResult<Option<bool>> {
        Ok(match self.input.get("PartO_active_cooling_required") {
            None => None,
            Some(JsonValue::Bool(whether)) => Some(*whether),
            Some(_) => {
                return Err(json_error(
                    "PartO_active_cooling_required field not a boolean",
                ))
            }
        })
    }

    #[cfg(test)]
    pub(crate) fn set_part_o_active_cooling_required(
        &mut self,
        required: bool,
    ) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key("PartO_active_cooling_required", json!(required))
    }

    #[cfg(test)]
    pub(crate) fn set_zone(&mut self, zone: JsonValue) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key("Zone", zone)
    }
}

impl TryFrom<&InputForProcessing> for Corpus {
    type Error = anyhow::Error;

    fn try_from(input: &InputForProcessing) -> Result<Self, Self::Error> {
        Corpus::from_inputs(
            &serde_json::from_value(input.input.to_owned())?,
            None,
            None,
            &Default::default(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;
    use std::fs::File;
    use walkdir::{DirEntry, WalkDir};

    #[fixture]
    fn core_files() -> Vec<DirEntry> {
        WalkDir::new("./examples/input/core")
            .into_iter()
            .filter_map(Result::ok)
            .filter(|e| {
                !e.file_type().is_dir()
                    && e.file_name().to_str().unwrap().ends_with("json")
                    && !e
                        .path()
                        .parent()
                        .unwrap()
                        .to_str()
                        .unwrap()
                        .ends_with("results") // don't test against files in results output directories
            })
            .collect_vec()
    }

    #[rstest]
    fn should_successfully_parse_all_demo_files(core_files: Vec<DirEntry>) {
        for entry in core_files {
            let parsed = ingest_for_processing(File::open(entry.path()).unwrap());
            assert!(
                parsed.is_ok(),
                "error was {:?} when parsing file {}",
                parsed.err().unwrap(),
                entry.file_name().to_str().unwrap()
            );
        }
    }

    /// this test should become redundant once there is an official JSON Schema definition for HEM.
    /// until that is the case, this checks each example input JSON against the schema generated by the schema-gen package.
    #[rstest]
    fn test_all_demo_files_pass_input_schema(core_files: Vec<DirEntry>) {
        let mut erroring_files: usize = Default::default();
        let mut error_outputs: Vec<std::string::String> = Default::default();

        for entry in core_files {
            let json_to_validate =
                serde_json::from_reader(BufReader::new(File::open(entry.path()).unwrap())).unwrap();
            let errors = FHS_SCHEMA_VALIDATOR.iter_errors(&json_to_validate);
            let mut at_least_one_error = false;
            for error in errors {
                at_least_one_error = true;
                error_outputs.push(format!(
                    "{} at path \"{}\":\n{error:?}",
                    entry.file_name().to_str().unwrap().to_owned(),
                    error.instance_path
                ));
            }
            if at_least_one_error {
                erroring_files += 1;
            }
        }

        assert!(
            error_outputs.is_empty(),
            "Schema created {} individual validation errors across {} files.\n\n{}",
            error_outputs.len(),
            erroring_files,
            error_outputs.join("\n\n")
        );
    }

    #[rstest]
    fn test_all_demo_files_deserialize_and_serialize(core_files: Vec<DirEntry>) {
        for entry in core_files {
            let input: Input =
                serde_json::from_reader(BufReader::new(File::open(entry.path()).unwrap())).expect(
                    format!(
                        "Failed deserializing {}",
                        entry.file_name().to_str().unwrap()
                    )
                    .as_str(),
                );
            let json = serde_json::to_string_pretty(&input.clone()).unwrap();
            let recreated_input: Input = serde_json::from_str(&json).unwrap();
            assert_eq!(input, recreated_input,);
        }
    }
}

#[cfg(test)]
mod accessors_tests {
    use super::*;
    use rstest::*;

    #[fixture]
    fn events_input() -> InputForProcessing {
        let events_input_json = json!({
            "Events": {
                "Shower": {
                  "IES": [
                    {
                      "start": 4.1,
                      "duration": 6,
                      "temperature": 41.0
                    },
                    {
                      "start": 4.5,
                      "duration": 6,
                      "temperature": 41.0
                    },
                    {
                      "start": 6,
                      "duration": 6,
                      "temperature": 41.0
                    }
                  ],
                  "mixer": [
                    {
                      "start": 7,
                      "duration": 6,
                      "temperature": 41.0
                    }
                  ]
                }
          }
        });

        InputForProcessing {
            input: events_input_json,
        }
    }

    #[rstest]
    fn test_water_heating_event_by_type_and_name_when_exists(events_input: InputForProcessing) {
        assert_eq!(
            events_input
                .water_heating_event_by_type_and_name("Shower", "mixer")
                .unwrap(),
            Some(vec![WaterHeatingEvent {
                start: 7.,
                duration: Some(6.),
                volume: None,
                temperature: 41.0
            }])
        );
    }

    #[fixture]
    fn hot_water_cylinder_input() -> InputForProcessing {
        let hot_water_source_json = json!({
            "HotWaterSource": {
                "hw cylinder": {
                  "type": "StorageTank",
                  "volume": 80.0,
                  "daily_losses": 1.68,
                  "min_temp": 52.0,
                  "setpoint_temp": 55.0,
                  "ColdWaterSource": "mains water",
                  "HeatSource": {
                    "hp": {
                      "type": "HeatSourceWet",
                      "name": "hp",
                      "temp_flow_limit_upper": 65,
                      "ColdWaterSource": "mains water",
                      "EnergySupply": "mains elec",
                      "Control": "hw timer",
                      "heater_position": 0.1,
                      "thermostat_position": 0.33
                    }
                  },
                  "primary_pipework": [
                    {
                      "location": "external",
                      "internal_diameter_mm": 26.,
                      "external_diameter_mm": 28.,
                      "length": 3.0,
                      "insulation_thermal_conductivity": 0.037,
                      "insulation_thickness_mm": 25.,
                      "surface_reflectivity": false,
                      "pipe_contents": "water"
                    }
                  ]
                }
              }
        });

        InputForProcessing {
            input: hot_water_source_json,
        }
    }

    #[rstest]
    fn test_hot_water_cylinder_volume(hot_water_cylinder_input: InputForProcessing) {
        assert_eq!(
            hot_water_cylinder_input
                .hot_water_cylinder_volume()
                .unwrap(),
            Some(80.0)
        )
    }

    #[rstest]
    fn test_set_gains_for_field() {
        let mut input = InputForProcessing {
            input: json!({
                "ApplianceGains": {
                    "Clothes_washing": 2,
                }
            }),
        };

        let expected_appliance_gains = json!({
            "Clothes_washing": 2,
            "Clothes_drying": 42,
        });

        input
            .set_gains_for_field("Clothes_drying", json!(42))
            .unwrap();

        assert_eq!(
            json!(input.root_object("ApplianceGains").unwrap()),
            expected_appliance_gains
        );
    }

    #[rstest]
    fn test_water_heating_events_of_types(events_input: InputForProcessing) {
        let actual = events_input
            .water_heating_events_of_types(&["Shower"])
            .unwrap();
        let expected = vec![
            json!({
              "start": 4.1,
              "duration": 6,
              "temperature": 41.0
            }),
            json!({
              "start": 4.5,
              "duration": 6,
              "temperature": 41.0
            }),
            json!({
              "start": 6,
              "duration": 6,
              "temperature": 41.0
            }),
            json!({
              "start": 7,
              "duration": 6,
              "temperature": 41.0
            }),
        ];
        assert_eq!(actual, expected);
    }
}
