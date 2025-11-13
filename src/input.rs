#![allow(unused_variables)]

use crate::core::heating_systems::heat_pump::TestLetter;
use crate::core::schedule::{BooleanSchedule, NumericSchedule};
use crate::external_conditions::{DaylightSavingsConfig, ShadingSegment, WindowShadingObject};
use crate::simulation_time::SimulationTime;
use anyhow::{anyhow, bail};
use indexmap::IndexMap;
use jsonschema::Validator;
use serde::de::Error;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use serde_enum_str::{Deserialize_enum_str, Serialize_enum_str};
use serde_json::{json, Map, Value as JsonValue};
use serde_repr::{Deserialize_repr, Serialize_repr};
use serde_valid::validation::error::{Format, Message};
use serde_valid::{MinimumError, Validate};
use smartstring::alias::String;
use std::fmt::{Display, Formatter};
use std::ops::Index;
use std::sync::Arc;
use std::sync::LazyLock;

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(rename_all = "PascalCase")] // TODO: add `deny_unknown_fields` declaration back in for versions newer than 0.36
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
    pub energy_supply: EnergySupplyInput,
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
    pub zone: ZoneDictionary,
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
pub struct ExternalConditionsInput {
    /// List of external air temperatures, one entry per hour (unit: ˚C)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub air_temperatures: Option<Vec<f64>>,
    /// List of wind speeds, one entry per hour (unit: m/s)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub wind_speeds: Option<Vec<f64>>,
    /// List of wind directions in degrees where North=0, East=90, South=180, West=270. Values range: 0 to 360. Wind direction is reported by the direction from which it originates, e.g. a southerly (180 degree) wind blows from the south to the north. (unit: ˚)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub wind_directions: Option<Vec<f64>>,
    /// List of diffuse horizontal radiation values, one entry per hour (unit: W/m²)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub diffuse_horizontal_radiation: Option<Vec<f64>>,
    /// List of direct beam radiation values, one entry per hour (unit: W/m²)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub direct_beam_radiation: Option<Vec<f64>>,
    /// List of ground reflectivity values, 0 to 1, one entry per hour
    #[serde(skip_serializing_if = "Option::is_none")]
    pub solar_reflectivity_of_ground: Option<Vec<f64>>,
    /// Latitude of weather station, angle from south (unit: ˚)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub latitude: Option<f64>,
    /// Longitude of weather station, easterly +ve westerly -ve (unit: ˚)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub longitude: Option<f64>,
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
    pub direct_beam_conversion_needed: Option<bool>,
    /// Data splitting the ground plane into segments (8-36) and giving height and distance to shading objects surrounding the building
    #[serde(skip_serializing_if = "Option::is_none")]
    pub shading_segments: Option<Vec<ShadingSegment>>,
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
// #[serde(deny_unknown_fields)] // TODO: add back in for versions after 0.36
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
pub struct EnergySupplyDetails {
    pub fuel: FuelType,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) diverter: Option<EnergyDiverter>,
    #[serde(rename = "ElectricBattery", skip_serializing_if = "Option::is_none")]
    pub(crate) electric_battery: Option<ElectricBattery>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub factor: Option<CustomEnergySourceFactor>,
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
    #[cfg(feature = "fhs")] // TODO review after migration as we expect it may be removed
    pub fn with_fuel(fuel_type: FuelType) -> Self {
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
pub enum FuelType {
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
pub enum EnergySupplyType {
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

pub type ColdWaterSourceInput = IndexMap<ColdWaterSourceType, ColdWaterSourceDetails>;

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub struct ColdWaterSourceDetails {
    /// First day of the time series, day of the year, 0 to 365
    #[validate(minimum = 0)]
    #[validate(maximum = 365)]
    pub start_day: u32,
    /// List of cold water temperatures, one entry per hour (unit: ˚C)
    pub temperatures: Vec<f64>,
    /// Timestep of the time series data (unit: hours)
    pub time_series_step: f64,
}

pub(crate) type ExtraControls = IndexMap<String, ControlDetails>;

#[derive(Clone, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
pub struct Control {
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
pub struct SmartApplianceBattery {
    pub battery_state_of_charge: IndexMap<String, Vec<f64>>,
    #[serde(default)]
    pub energy_into_battery_from_generation: IndexMap<String, Vec<f64>>,
    #[serde(default)]
    pub energy_into_battery_from_grid: IndexMap<String, Vec<f64>>,
    #[serde(default)]
    pub energy_out_of_battery: IndexMap<String, Vec<f64>>,
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
pub struct HotWaterSource {
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

pub trait HotWaterSourceDetailsForProcessing {
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

pub struct HotWaterSourceDetailsJsonMap<'a>(pub &'a mut Map<std::string::String, JsonValue>);

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
#[serde(tag = "type")] // TODO: possibly restore `deny_unknown_fields` annotation after 0.36
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
pub struct Baths(pub IndexMap<String, BathDetails>);

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct BathDetails {
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
pub struct OtherWaterUses(pub IndexMap<String, OtherWaterUse>);

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct OtherWaterUse {
    /// Tap/outlet flow rate (unit: litre/minute)
    pub(crate) flowrate: f64,
    #[serde(rename = "ColdWaterSource")]
    pub(crate) cold_water_source: ColdWaterSourceType,
}

pub type WaterDistribution = Vec<WaterPipeworkSimple>;

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

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
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
#[serde(tag = "type")] // TODO: possibly restore `deny_unknown_fields` annotation after 0.36
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
        variable_flow: Option<bool>, // TODO: restore as non-Option after 0.36 if possible
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
// #[serde(deny_unknown_fields)] // TODO: restore after 0.36 if possible
pub struct ZoneInput {
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
    pub area: f64,
    /// Total volume of the zone. (Unit: m³)
    #[validate(minimum = 0.)]
    pub(crate) volume: f64,
    /// Basis for zone temperature control.
    #[serde(rename = "temp_setpnt_basis", skip_serializing_if = "Option::is_none")]
    pub(crate) temp_setpnt_basis: Option<ZoneTemperatureControlBasis>,
    /// Setpoint temperature to use during initialisation (unit: ˚C)
    pub(crate) temp_setpnt_init: Option<f64>, // restore as non `Option` after 0.36
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
pub enum BuildingElement {
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
    pub fn pitch(&self) -> f64 {
        *match self {
            BuildingElement::Opaque { pitch, .. } => pitch,
            BuildingElement::Transparent { pitch, .. } => pitch,
            BuildingElement::Ground { pitch, .. } => pitch,
            BuildingElement::AdjacentConditionedSpace { pitch, .. } => pitch,
            BuildingElement::AdjacentUnconditionedSpace { pitch, .. } => pitch,
        }
    }

    pub fn height(&self) -> Option<f64> {
        match self {
            BuildingElement::Opaque { height, .. } => Some(*height),
            BuildingElement::Transparent { height, .. } => Some(*height),
            _ => None,
        }
    }

    pub fn width(&self) -> Option<f64> {
        match self {
            BuildingElement::Opaque { width, .. } => Some(*width),
            BuildingElement::Transparent { width, .. } => Some(*width),
            _ => None,
        }
    }

    pub fn u_value(&self) -> Option<f64> {
        match self {
            BuildingElement::Opaque { u_value, .. } => *u_value,
            BuildingElement::Transparent { u_value, .. } => *u_value,
            BuildingElement::Ground { u_value, .. } => Some(*u_value),
            BuildingElement::AdjacentConditionedSpace { u_value, .. } => *u_value,
            BuildingElement::AdjacentUnconditionedSpace { u_value, .. } => *u_value,
        }
    }

    pub fn orientation(&self) -> Option<f64> {
        match self {
            BuildingElement::Opaque { orientation, .. } => Some(*orientation),
            BuildingElement::Transparent { orientation, .. } => Some(*orientation),
            _ => None,
        }
    }
}

pub trait TransparentBuildingElement {
    fn set_window_openable_control(&mut self, control: &str);
    fn is_security_risk(&self) -> bool;
    fn treatment(&mut self) -> Option<Vec<&mut Map<std::string::String, JsonValue>>>;
}

pub struct TransparentBuildingElementJsonValue<'a>(pub &'a mut Map<std::string::String, JsonValue>);

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

pub trait GroundBuildingElement {
    fn set_u_value(&mut self, new_u_value: f64);
    fn set_thermal_resistance_floor_construction(&mut self, new_r_f: f64);
    fn set_psi_wall_floor_junc(&mut self, new_psi_wall_floor_junc: f64);
}

pub struct GroundBuildingElementJsonValue<'a>(pub &'a mut Map<std::string::String, JsonValue>);

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
pub enum MassDistributionClass {
    D,
    E,
    I,
    IE,
    M,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct WindowPart {
    pub(crate) mid_height_air_flow_path: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct WindowTreatment {
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
pub enum FloorData {
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
#[serde(tag = "type")] // TODO: possibly restore `deny_unknown_fields` serde annotation after 0.36 (once FHS is extracted)
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
#[serde(tag = "type")] // TODO: possibly restore `deny_unknown_fields` serde annotation after 0.36 (once FHS is extracted)
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

pub type HeatSourceWet = IndexMap<String, HeatSourceWetDetails>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[allow(clippy::large_enum_variant)]
#[serde(tag = "type", deny_unknown_fields)]
pub enum HeatSourceWetDetails {
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
pub enum HeatPumpSourceType {
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
pub enum HeatPumpSinkType {
    Water,
    Air,
    Glycol25,
}

#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum HeatPumpBackupControlType {
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
pub struct HeatPumpTestDatum {
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

pub type WasteWaterHeatRecovery = IndexMap<String, WasteWaterHeatRecoveryDetails>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct WasteWaterHeatRecoveryDetails {
    #[serde(rename = "type")]
    pub(crate) system_type: WasteWaterHeatRecoverySystemType,
    #[serde(rename = "ColdWaterSource")]
    pub cold_water_source: ColdWaterSourceType,
    pub flow_rates: Vec<f64>,
    pub efficiencies: Vec<f64>,
    pub utilisation_factor: f64,
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
// #[serde(deny_unknown_fields)] // TODO: restore `deny_unknown_fields` declaration after 0.36
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
// #[serde(deny_unknown_fields)] // TODO: possibly restore `deny_unknown_fields` declaration after 0.36
pub struct MechanicalVentilation {
    #[serde(rename = "sup_air_flw_ctrl")]
    pub(crate) supply_air_flow_rate_control: SupplyAirFlowRateControlType,
    #[serde(rename = "sup_air_temp_ctrl")]
    pub(crate) supply_air_temperature_control_type: SupplyAirTemperatureControlType,
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

pub trait MechanicalVentilationForProcessing {
    fn vent_is_type(&self, vent_type: &str) -> bool;
    fn measured_fan_power(&self) -> Option<f64>;
    fn measured_air_flow_rate(&self) -> Option<f64>;
    fn set_sfp(&mut self, sfp: f64);
    fn set_control(&mut self, control: &str);
}

pub struct MechanicalVentilationJsonValue<'a>(pub &'a mut Map<std::string::String, JsonValue>);

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
// #[serde(deny_unknown_fields)] // TODO: possibly restore `deny_unknown_fields` declaration after 0.36
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

// The calc_htc_hlp function in the corpus module needs reduced access to an input
// though this may be in the context of a wrapper which cannot guarantee that other data is in the
// right shape. Abstracting this to a trait allows us to define an input that allows partial deserialisation
// from the underlying JSON, and so can ignore areas of the input that it does not refer to.
pub trait InputForCalcHtcHlp {
    fn simulation_time(&self) -> &SimulationTime;
    fn energy_supply(&self) -> &EnergySupplyInput;
    fn external_conditions(&self) -> &ExternalConditionsInput;
    fn control(&self) -> &Control;
    fn infiltration_ventilation(&self) -> &InfiltrationVentilation;
    fn zone(&self) -> &ZoneDictionary;
    fn temp_internal_air_static_calcs(&self) -> f64;
}

impl InputForCalcHtcHlp for Input {
    fn simulation_time(&self) -> &SimulationTime {
        &self.simulation_time
    }

    fn energy_supply(&self) -> &EnergySupplyInput {
        &self.energy_supply
    }

    fn external_conditions(&self) -> &ExternalConditionsInput {
        self.external_conditions.as_ref()
    }

    fn control(&self) -> &Control {
        &self.control
    }

    fn infiltration_ventilation(&self) -> &InfiltrationVentilation {
        &self.infiltration_ventilation
    }

    fn zone(&self) -> &ZoneDictionary {
        &self.zone
    }

    fn temp_internal_air_static_calcs(&self) -> f64 {
        self.temp_internal_air_static_calcs
    }
}

// The purpose of this struct is to allow deserialisation of an input just containing the data needed for the
// calc_htc_hlp function in the corpus module, so that we can ignore other areas of the data that may not be in the
// expected shape for a core input.
#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(rename_all = "PascalCase")]
pub struct ReducedInputForCalcHtcHlp {
    #[serde(rename = "temp_internal_air_static_calcs")]
    pub(crate) temp_internal_air_static_calcs: f64,
    pub(crate) simulation_time: SimulationTime,
    pub(crate) external_conditions: Arc<ExternalConditionsInput>,
    pub(crate) energy_supply: EnergySupplyInput,
    pub(crate) control: Control,
    pub(crate) zone: ZoneDictionary,
    pub(crate) infiltration_ventilation: InfiltrationVentilation,
}

impl InputForCalcHtcHlp for ReducedInputForCalcHtcHlp {
    fn simulation_time(&self) -> &SimulationTime {
        &self.simulation_time
    }

    fn energy_supply(&self) -> &EnergySupplyInput {
        &self.energy_supply
    }

    fn external_conditions(&self) -> &ExternalConditionsInput {
        self.external_conditions.as_ref()
    }

    fn control(&self) -> &Control {
        &self.control
    }

    fn infiltration_ventilation(&self) -> &InfiltrationVentilation {
        &self.infiltration_ventilation
    }

    fn zone(&self) -> &ZoneDictionary {
        &self.zone
    }

    fn temp_internal_air_static_calcs(&self) -> f64 {
        self.temp_internal_air_static_calcs
    }
}

#[expect(unused)]
static CORE_SCHEMA_VALIDATOR: LazyLock<Validator> = LazyLock::new(|| {
    let schema = serde_json::from_str(include_str!("../schemas/input_core.schema.json")).unwrap();
    jsonschema::validator_for(&schema).unwrap()
});

#[expect(unused)]
static CORE_INCLUDING_FHS_VALIDATOR: LazyLock<Validator> = LazyLock::new(|| {
    let schema = serde_json::from_str(include_str!(
        "../schemas/input_core_allowing_fhs.schema.json"
    ))
    .unwrap();
    jsonschema::validator_for(&schema).unwrap()
});

#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;
    use rstest::*;
    use std::fs::File;
    use std::io::BufReader;
    use walkdir::{DirEntry, WalkDir};

    #[fixture]
    fn core_files() -> Vec<DirEntry> {
        files_with_root("./examples/input/core")
    }

    fn files_with_root(root: &str) -> Vec<DirEntry> {
        WalkDir::new(root)
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
    fn test_all_demo_files_deserialize_and_serialize(core_files: Vec<DirEntry>) {
        for entry in core_files {
            let input: Input =
                serde_json::from_reader(BufReader::new(File::open(entry.path()).unwrap()))
                    .unwrap_or_else(|_| {
                        panic!(
                            "Failed deserializing {}",
                            entry.file_name().to_str().unwrap()
                        )
                    });
            let json = serde_json::to_string_pretty(&input.clone()).unwrap();
            let recreated_input: Input = serde_json::from_str(&json).unwrap();
            assert_eq!(input, recreated_input,);
        }
    }
}
