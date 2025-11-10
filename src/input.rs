#![allow(unused_variables)]

use crate::core::heating_systems::heat_pump::TestLetter;
use crate::core::schedule::{BooleanSchedule, NumericSchedule};
use crate::core::units::calculate_thermal_resistance_of_virtual_layer;
use crate::corpus::Corpus;
use crate::external_conditions::{ShadingSegment, WindowShadingObject};
use crate::simulation_time::SimulationTime;
use crate::{HEM_VERSION, HOURS_TO_END_DEC};
use anyhow::{anyhow, bail};
use indexmap::IndexMap;
use itertools::Itertools;
use jsonschema::{BasicOutput, Validator};
use monostate::MustBe;
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

pub fn ingest_for_processing(
    json: impl Read,
    schema_reference: &SchemaReference,
) -> Result<InputForProcessing, anyhow::Error> {
    InputForProcessing::init_with_json(json, schema_reference)
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum SchemaReference {
    Core,
    #[cfg(feature = "fhs")]
    Fhs,
}

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(rename_all = "PascalCase")] // TODO: add `deny_unknown_fields` declaration back in for versions newer than 0.36
#[validate(custom = validate_shower_waste_water_heat_recovery_systems)]
#[validate(custom = validate_exhaust_air_heat_pump_ventilation_compatibility)]
pub struct Input {
    /// Metadata for the input file
    #[serde(rename = "metadata")]
    #[validate]
    metadata: InputMetadata,

    #[serde(default)]
    #[validate]
    pub(crate) appliance_gains: ApplianceGains,

    #[validate]
    pub(crate) cold_water_source: ColdWaterSourceInput,

    #[validate]
    pub(crate) control: Control,

    #[validate]
    pub(crate) energy_supply: EnergySupplyInput,

    #[serde(rename = "Events")]
    #[validate]
    pub(crate) water_heating_events: WaterHeatingEvents,

    #[validate]
    pub(crate) external_conditions: Arc<ExternalConditionsInput>,

    /// Dictionary of available wet heat sources, keyed by user-defined names (e.g., 'boiler', 'hp', 'HeatNetwork', 'hb1'). Other models reference these keys via their heat_source_wet fields.
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate]
    pub(crate) heat_source_wet: Option<HeatSourceWet>,

    #[validate]
    pub(crate) hot_water_demand: HotWaterDemand,

    #[validate]
    pub(crate) hot_water_source: HotWaterSource,

    #[validate]
    pub(crate) infiltration_ventilation: InfiltrationVentilation,

    #[validate]
    pub(crate) internal_gains: InternalGains,

    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate]
    pub(crate) on_site_generation: Option<OnSiteGeneration>,

    #[serde(default)]
    #[validate(custom = validate_only_storage_tanks)]
    #[validate]
    pub(crate) pre_heated_water_source: IndexMap<std::string::String, HotWaterSourceDetails>,

    #[validate]
    pub(crate) simulation_time: SimulationTime,

    #[serde(default, skip_serializing_if = "IndexMap::is_empty")]
    #[validate]
    pub(crate) smart_appliance_controls:
        IndexMap<std::string::String, SmartApplianceControlDetails>,

    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate]
    pub(crate) space_cool_system: Option<SpaceCoolSystem>,

    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate]
    pub(crate) space_heat_system: Option<SpaceHeatSystem>,

    #[serde(rename = "WWHRS", skip_serializing_if = "Option::is_none")]
    #[validate]
    pub(crate) waste_water_heat_recovery: Option<WasteWaterHeatRecovery>,

    #[validate]
    pub(crate) zone: ZoneDictionary,

    #[serde(rename = "temp_internal_air_static_calcs")]
    pub(crate) temp_internal_air_static_calcs: f64,
}

#[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize, Validate)]
pub(crate) struct InputMetadata {
    #[validate(custom = validate_hem_core_version)]
    hem_core_version: String,
}

fn validate_hem_core_version(hem_core_version: &str) -> Result<(), serde_valid::validation::Error> {
    (hem_core_version == HEM_VERSION)
        .then_some(())
        .ok_or_else(|| serde_valid::validation::Error::Custom(
            format!("hem_core_version was expected to be '{HEM_VERSION}' but was provided as '{hem_core_version}'")
        ))
}

fn validate_shower_waste_water_heat_recovery_systems(
    input: &Input,
) -> Result<(), serde_valid::validation::Error> {
    for (shower_name, shower) in input.hot_water_demand.shower.0.iter() {
        if let (
            Shower::MixerShower {
                wwhrs_config:
                    Some(MixerShowerWwhrsConfiguration {
                        waste_water_heat_recovery_system: wwhrs,
                        ..
                    }),
                ..
            },
            Some(wwhrs_systems),
        ) = (shower, input.waste_water_heat_recovery.as_ref())
        {
            if !wwhrs_systems.contains_key(wwhrs.as_str()) {
                return custom_validation_error(
                    format!("WWHRS value '{wwhrs}' not found in Input.WWHRS (from Input.HotWaterDemand.Shower['{shower_name}'].WWHRS)")
                );
            }
        }
    }

    Ok(())
}

/// Validate that exhaust air heat pumps are compatible with ventilation systems.
fn validate_exhaust_air_heat_pump_ventilation_compatibility(
    input: &Input,
) -> Result<(), serde_valid::validation::Error> {
    if let Some(heat_source_wet) = input.heat_source_wet.as_ref() {
        fn vent_is_incompatible(vent: &MechVentData) -> bool {
            matches!(
                vent,
                MechVentData::IntermittentMev { .. }
                    | MechVentData::DecentralisedContinuousMev { .. }
            )
        }

        let exhaust_air_source_types = [
            HeatPumpSourceType::ExhaustAirMEV,
            HeatPumpSourceType::ExhaustAirMVHR,
            HeatPumpSourceType::ExhaustAirMixed,
        ];

        let exhaust_air_heat_pumps: IndexMap<String, HeatPumpSourceType> = heat_source_wet
            .iter()
            .filter_map(|(name, source)| {
                if let HeatSourceWetDetails::HeatPump { source_type, .. } = source {
                    if exhaust_air_source_types.contains(&source_type) {
                        Some((name.into(), source_type.clone()))
                    } else {
                        None
                    }
                } else {
                    None
                }
            })
            .collect();

        if exhaust_air_heat_pumps.is_empty() {
            return Ok(());
        }

        let incompatible_vents: IndexMap<String, &'static str> = input
            .infiltration_ventilation
            .mechanical_ventilation
            .iter()
            .filter_map(|(vent_name, vent_data)| {
                vent_is_incompatible(&vent_data.vent_data)
                    .then_some((vent_name.into(), vent_data.vent_data.vent_type()))
            })
            .collect();

        if !incompatible_vents.is_empty() {
            let mut incompatibilities: Vec<String> = Default::default();
            for (vent_name, vent_type) in incompatible_vents.iter() {
                for (heat_source_name, heat_source_type) in exhaust_air_heat_pumps.iter() {
                    incompatibilities
                        .push(format!("Exhaust air heat pump '{heat_source_name}' is incompatible with ventilation system '{vent_name}' (vent type: {vent_type})").into());
                }
            }

            return custom_validation_error(
                format!(
                    "System incompatibilities found: {}. Exhaust air heat pumps do not work with Intermittent MEV or Decentralised continuous MEV.",
                    incompatibilities.join("; ")
                )
            );
        }
    }

    Ok(())
}

fn validate_only_storage_tanks(
    sources: &IndexMap<std::string::String, HotWaterSourceDetails>,
) -> Result<(), serde_valid::validation::Error> {
    sources.values()
        .all(|details| matches!(details, HotWaterSourceDetails::StorageTank { .. }))
        .then_some(())
        .ok_or_else(|| serde_valid::validation::Error::Custom("PreHeatedWaterSource input can only contain HotWaterSource data of the type StorageTank".to_owned()))
}

#[derive(Clone, Debug, Default, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub struct ExternalConditionsInput {
    /// List of external air temperatures, one entry per hour (unit: ˚C)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(custom = validate_air_temperatures)]
    pub(crate) air_temperatures: Option<Vec<f64>>,

    /// List of diffuse horizontal radiation values, one entry per hour (unit: W/m²)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(custom = validate_all_items_in_option_non_negative)]
    pub(crate) diffuse_horizontal_radiation: Option<Vec<f64>>,

    /// A flag to indicate whether direct beam radiation from climate data needs to be converted from horizontal to normal incidence; if normal direct beam radiation values are provided then no conversion is needed
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) direct_beam_conversion_needed: Option<bool>,

    /// List of direct beam radiation values, one entry per hour (unit: W/m²)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(custom = validate_all_items_in_option_non_negative)]
    pub(crate) direct_beam_radiation: Option<Vec<f64>>,

    /// Latitude of weather station, angle from south (unit: ˚)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = -90.)]
    #[validate(maximum = 90.)]
    pub(crate) latitude: Option<f64>,

    /// Longitude of weather station, easterly +ve westerly -ve (unit: ˚)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = -180.)]
    #[validate(maximum = 180.)]
    pub(crate) longitude: Option<f64>,

    /// Data splitting the ground plane into segments (8-36) and giving height and distance to shading objects surrounding the building
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate]
    pub(crate) shading_segments: Option<Vec<ShadingSegment>>,

    /// List of ground reflectivity values, 0 to 1, one entry per hour
    /// Each item represents the fraction of solar radiation incident on the ground that is reflected. Also called the albedo.
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(custom = validate_all_items_in_option_non_negative)]
    #[validate(custom = |v| validate_all_items_in_option_at_most_n(v, 1.))]
    pub(crate) solar_reflectivity_of_ground: Option<Vec<f64>>,

    /// List of wind directions in degrees where North=0, East=90, South=180, West=270. Values range: 0 to 360. Wind direction is reported by the direction from which it originates, e.g. a southerly (180 degree) wind blows from the south to the north. (unit: ˚)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(custom = validate_all_items_in_option_non_negative)]
    #[validate(custom = |v| validate_all_items_in_option_at_most_n(v, 360.))]
    pub(crate) wind_directions: Option<Vec<f64>>,

    /// List of wind speeds, one entry per hour (unit: m/s)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(custom = validate_all_items_in_option_non_negative)]
    pub(crate) wind_speeds: Option<Vec<f64>>,
}

impl ExternalConditionsInput {
    pub(crate) fn are_all_fields_set(&self) -> bool {
        [
            self.air_temperatures.is_some(),
            self.diffuse_horizontal_radiation.is_some(),
            self.direct_beam_conversion_needed.is_some(),
            self.direct_beam_radiation.is_some(),
            self.latitude.is_some(),
            self.longitude.is_some(),
            self.shading_segments.is_some(),
            self.solar_reflectivity_of_ground.is_some(),
            self.wind_directions.is_some(),
            self.wind_speeds.is_some(),
        ]
        .iter()
        .all(|whether| *whether)
            && if let Self {
                air_temperatures: Some(air_temperatures),
                diffuse_horizontal_radiation: Some(diffuse_horizontal_radiation),
                direct_beam_radiation: Some(direct_beam_radiation),
                solar_reflectivity_of_ground: Some(solar_reflectivity_of_ground),
                wind_directions: Some(wind_directions),
                wind_speeds: Some(wind_speeds),
                ..
            } = self
            {
                [
                    air_temperatures,
                    diffuse_horizontal_radiation,
                    direct_beam_radiation,
                    solar_reflectivity_of_ground,
                    wind_directions,
                    wind_speeds,
                ]
                .iter()
                .all(|items| items.len() >= HOURS_TO_END_DEC as usize)
            } else {
                false
            }
    }
}

fn validate_air_temperatures(
    air_temps: &Option<Vec<f64>>,
) -> Result<(), serde_valid::validation::Error> {
    if air_temps.iter().flatten().all(|&v| v >= -273.15) {
        Ok(())
    } else {
        custom_validation_error(
            "Some air temperatures contained values that were below -273.15˚C.".to_string(),
        )
    }
}

#[derive(Clone, Debug, Default, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
pub(crate) struct InternalGains {
    #[serde(
        alias = "total internal gains",
        skip_serializing_if = "Option::is_none"
    )]
    #[validate]
    pub(crate) total_internal_gains: Option<InternalGainsDetails>,

    #[serde(rename = "metabolic gains", skip_serializing_if = "Option::is_none")]
    #[validate]
    pub(crate) metabolic_gains: Option<InternalGainsDetails>,

    #[serde(rename = "EvaporativeLosses", skip_serializing_if = "Option::is_none")]
    #[validate]
    pub(crate) evaporative_losses: Option<InternalGainsDetails>,

    #[serde(rename = "ColdWaterLosses", skip_serializing_if = "Option::is_none")]
    #[validate]
    pub(crate) cold_water_losses: Option<InternalGainsDetails>,

    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate]
    pub(crate) other: Option<InternalGainsDetails>,
}

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub(crate) struct InternalGainsDetails {
    #[validate(minimum = 0)]
    #[validate(maximum = 365)]
    pub(crate) start_day: u32,
    #[validate(minimum = 0.)]
    #[validate(maximum = 24.)]
    pub(crate) time_series_step: f64,
    pub(crate) schedule: NumericSchedule,
}

pub(crate) type ApplianceGains = IndexMap<std::string::String, ApplianceGainsDetails>;

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub(crate) struct ApplianceGainsDetails {
    /// List of appliance usage events
    #[serde(rename = "Events", skip_serializing_if = "Option::is_none")]
    #[validate]
    pub(crate) events: Option<Vec<ApplianceGainsEvent>>,

    /// Appliance power consumption when not in use (unit: W)
    #[serde(rename = "Standby", skip_serializing_if = "Option::is_none")]
    #[validate(minimum = 0.)]
    pub(crate) standby: Option<f64>,

    #[serde(rename = "EnergySupply")]
    pub(crate) energy_supply: String,

    /// Proportion of appliance demand turned into heat gains (dimensionless, 0-1)
    #[validate(minimum = 0.)]
    #[validate(maximum = 1.)]
    pub(crate) gains_fraction: f64,

    /// Load shifting configuration for smart appliance control
    #[serde(rename = "loadshifting", skip_serializing_if = "Option::is_none")]
    #[validate]
    pub(crate) load_shifting: Option<ApplianceLoadShifting>,

    /// Priority level for load shifting (lower numbers = higher priority)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) priority: Option<isize>,

    /// Power consumption schedule (one entry per hour)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) schedule: Option<NumericSchedule>,

    /// First day of the time series, day of the year, 0 to 365
    #[validate(minimum = 0)]
    #[validate(maximum = 365)]
    pub(crate) start_day: u32,

    /// Timestep of the time series data (unit: hours)
    #[validate(minimum = 0.)]
    #[validate(maximum = 24.)]
    pub(crate) time_series_step: f64,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct ApplianceGainsEvent {
    /// Hours from start of simulation to when the appliance event begins (unit: hours)
    #[validate(minimum = 0.)]
    pub start: f64,

    /// Duration of the appliance event (unit: hours)
    #[validate(minimum = 0.)]
    pub duration: f64,

    /// Electrical power consumption during the appliance event (unit: W)
    #[serde(rename = "demand_W")]
    #[validate(minimum = 0.)]
    pub demand_w: f64,
}

pub(crate) type EnergySupplyInput = IndexMap<std::string::String, EnergySupplyDetails>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields, rename_all = "snake_case")]
pub struct EnergySupplyDetails {
    /// Type of combustion fuel
    pub(crate) fuel: FuelType,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) diverter: Option<EnergyDiverter>,

    /// Indicates that an electric battery is present
    #[serde(rename = "ElectricBattery", skip_serializing_if = "Option::is_none")]
    pub(crate) electric_battery: Option<ElectricBattery>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) factor: Option<CustomEnergySourceFactor>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) priority: Option<Vec<EnergySupplyPriorityEntry>>,

    /// Denotes that this energy supply can export its surplus supply
    pub(crate) is_export_capable: bool,

    /// Level of battery charge above which grid prohibited from charging battery (monthly values) (0 - 1)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(custom = validate_threshold_value_fractions)]
    pub(crate) threshold_charges: Option<[f64; 12]>,

    /// Grid price below which battery is permitted to charge from grid (monthly values) (unit: p/kWh)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(custom = validate_threshold_value_fractions)]
    pub(crate) threshold_prices: Option<[f64; 12]>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) tariff: Option<EnergySupplyTariff>,
}

fn validate_threshold_value_fractions(
    values: &Option<[f64; 12]>,
) -> Result<(), serde_valid::validation::Error> {
    if let Some(values) = values {
        if !values.iter().all(|&v| v >= 0. && v <= 1.) {
            return custom_validation_error("Some threshold values for an energy supply contained numbers that were not fractions between 0 and 1 inclusive.".to_string());
        }
    }

    Ok(())
}

impl EnergySupplyDetails {
    pub(crate) fn with_fuel(fuel_type: FuelType) -> Self {
        Self {
            fuel: fuel_type,
            diverter: None,
            electric_battery: None,
            factor: None,
            priority: None,
            is_export_capable: false,
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
    /// Reference to a control schedule of maximum temperature setpoints. References a key in $.Control.
    #[serde(rename = "Controlmax")]
    pub(crate) control_max: String,
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

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct ElectricBattery {
    /// The maximum capacity of the battery (unit: kWh)
    #[validate(minimum = 0.)]
    pub capacity: f64,

    /// Charge/discharge round trip efficiency of battery system (greater than 0, up to 1)
    #[validate(exclusive_minimum = 0.)]
    #[validate(maximum = 1.)]
    pub charge_discharge_efficiency_round_trip: f64,

    /// The starting age of the battery (in years)
    #[validate(minimum = 0.)]
    pub battery_age: f64,

    /// The maximum charge rate one way trip the battery allows (unit: kW)
    #[validate(minimum = 0.)]
    pub minimum_charge_rate_one_way_trip: f64,

    /// The maximum discharge rate one way trip the battery allows (unit: kW)
    #[validate(minimum = 0.)]
    pub maximum_charge_rate_one_way_trip: f64,

    /// The minimum charge rate one way trip the battery allows (unit: kW)
    #[validate(minimum = 0.)]
    pub maximum_discharge_rate_one_way_trip: f64,

    /// The location of the battery (inside/outside)
    pub battery_location: BatteryLocation,

    /// Is charging from the grid possible?
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

pub(crate) type ColdWaterSourceInput = IndexMap<std::string::String, ColdWaterSourceDetails>;

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
    #[validate(custom = validate_running_water_temperatures)]
    pub(crate) temperatures: Vec<f64>,

    /// Duration in hours (must be within 24-hour period)
    #[validate(minimum = 0.)]
    #[validate(maximum = 24.)]
    pub(crate) time_series_step: f64,
}

fn validate_running_water_temperatures(
    temps: &[f64],
) -> Result<(), serde_valid::validation::Error> {
    if !temps.iter().all(|&t| t >= 0. && t <= 100.) {
        custom_validation_error("Some cold water temperatures contained numbers that were not within the range 0 to 100 inclusive.".into())
    } else {
        Ok(())
    }
}

pub(crate) type ExtraControls = IndexMap<std::string::String, ControlDetails>;

/// Control schedule configuration for heating and energy systems.
#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
pub struct Control {
    #[serde(skip_serializing_if = "Option::is_none", rename = "hw timer")]
    #[validate]
    pub(crate) hot_water_timer: Option<ControlDetails>,

    #[serde(skip_serializing_if = "Option::is_none", rename = "window opening")]
    #[validate]
    pub(crate) window_opening: Option<ControlDetails>,

    #[serde(flatten)]
    #[validate]
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
#[serde(tag = "type")]
#[validate(custom = validate_logic_type_for_charge_target)]
pub(crate) enum ControlDetails {
    #[serde(rename = "OnOffTimeControl")]
    OnOffTimer {
        #[serde(skip_serializing_if = "Option::is_none")]
        allow_null: Option<bool>,

        /// First day of the time series, day of the year, 0 to 365
        #[validate(minimum = 0)]
        #[validate(maximum = 365)]
        start_day: u32,

        #[validate(minimum = 0.)]
        #[validate(maximum = 24.)]
        time_series_step: f64,

        /// List of boolean values where true means on, one entry per hour
        schedule: BooleanSchedule,
    },
    #[serde(rename = "OnOffCostMinimisingTimeControl")]
    OnOffCostMinimising {
        /// First day of the time series, day of the year, 0 to 365
        #[validate(minimum = 0)]
        #[validate(maximum = 365)]
        start_day: u32,

        /// Timestep of the time series data (unit: hours)
        #[validate(minimum = 0.)]
        #[validate(maximum = 24.)]
        time_series_step: f64,

        /// Number of 'on' hours to be set per day
        #[validate(exclusive_minimum = 0.)]
        #[validate(maximum = 24.)]
        time_on_daily: f64,

        /// List of cost values (one entry per time_series_step)
        schedule: NumericSchedule,
    },
    #[serde(rename = "SetpointTimeControl")]
    SetpointTimer {
        /// First day of the time series, day of the year, 0 to 365
        #[validate(minimum = 0)]
        #[validate(maximum = 365)]
        start_day: u32,

        /// Timestep of the time series data (unit: hours)
        #[validate(minimum = 0.)]
        #[validate(maximum = 24.)]
        time_series_step: f64,

        /// How long before heating period the system should switch on (unit: hours)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        advanced_start: Option<f64>,

        /// list of float values (one entry per hour)
        schedule: NumericSchedule,

        #[serde(flatten)]
        #[validate]
        setpoint_bounds: Option<SetpointBoundsInput>,
    },
    #[serde(rename = "ChargeControl")]
    ChargeTarget {
        /// Proportion of the charge targeted for each day
        #[serde(skip_serializing_if = "Option::is_none")]
        charge_level: Option<ChargeLevel>,

        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate]
        external_sensor: Option<ExternalSensor>,

        #[serde(skip_serializing_if = "Option::is_none")]
        logic_type: Option<ControlLogicType>,

        /// List of boolean values where true means 'on' (one entry per hour)
        schedule: BooleanSchedule,

        /// Temperature at which charging should stop
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = -273.15)]
        temp_charge_cut: Option<f64>,

        /// Temperature delta schedule for charge cut-off adjustment (unit: ˚C)
        #[serde(skip_serializing_if = "Option::is_none")]
        temp_charge_cut_delta: Option<NumericSchedule>,

        /// Indicates from which hour of the day the system starts to target the charge level for the next day rather than the current day
        #[serde(default = "default_charge_calc_time")]
        #[validate(minimum = 0.)]
        #[validate(exclusive_maximum = 24.)]
        charge_calc_time: f64,

        /// First day of the time series, day of the year, 0 to 365
        #[validate(minimum = 0)]
        #[validate(maximum = 365)]
        start_day: u32,

        #[validate(minimum = 0.)]
        #[validate(maximum = 24.)]
        time_series_step: f64,
    },
    #[serde(rename = "CombinationTimeControl")]
    CombinationTime { combination: ControlCombinations },
}

const fn default_charge_calc_time() -> f64 {
    21.
}

fn validate_logic_type_for_charge_target(
    data: &ControlDetails,
) -> Result<(), serde_valid::validation::Error> {
    if let ControlDetails::ChargeTarget {
        logic_type:
            Some(ControlLogicType::Automatic | ControlLogicType::Celect | ControlLogicType::Hhrsh),
        temp_charge_cut: None,
        ..
    } = data
    {
        custom_validation_error(
            "logic types 'automatic', 'celect' and 'hhrsh' requires temp_charge_cut to be set"
                .into(),
        )
    } else {
        Ok(())
    }
}

/// Enum to encapsulate possible combinations of setpoint data
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(untagged)]
#[validate(custom = validate_setpoint_bounds)]
pub(crate) enum SetpointBoundsInput {
    MinAndMax {
        /// Minimum setpoint allowed
        setpoint_min: f64,

        /// Maximum setpoint allowed
        setpoint_max: f64,

        /// If both min and max limits are set but setpoint is not, whether to default to min (false) or max (true)
        default_to_max: bool,
    },
    MinOnly {
        /// Minimum setpoint allowed
        setpoint_min: f64,
    },
    MaxOnly {
        /// Maximum setpoint allowed
        setpoint_max: f64,
    },
}

fn validate_setpoint_bounds(
    bounds: &SetpointBoundsInput,
) -> Result<(), serde_valid::validation::Error> {
    match bounds {
        SetpointBoundsInput::MinAndMax {
            setpoint_min,
            setpoint_max,
            ..
        } if setpoint_max <= setpoint_min => {
            custom_validation_error("setpoint_max must be greater than setpoint_min".into())
        }
        _ => Ok(()),
    }
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

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct ExternalSensor {
    #[validate]
    pub(crate) correlation: Vec<ExternalSensorCorrelation>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct ExternalSensorCorrelation {
    /// External temperature data point for the corresponding maximum charge level (unit: Celsius)
    #[validate(minimum = -273.15)]
    pub(crate) temperature: f64,

    /// Maximum charge level permitted by the control for a given external temperature
    #[validate(minimum = 0.)]
    #[validate(maximum = 1.)]
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

#[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct SmartApplianceBattery {
    /// Dictionary of lists containing the battery state of charge for each timestep for each energy supply
    #[validate(custom = validate_battery_state_fractions)]
    #[validate(custom = |v| validate_all_sublists_non_empty(v, "SmartApplianceBattery"))]
    #[validate(custom = |v| validate_map_non_empty(v, "SmartApplianceBattery"))]
    pub(crate) battery_state_of_charge: IndexMap<String, Vec<f64>>,

    /// Dictionary of lists containing energy sent to the battery from generation for each timestep for each energy supply (unit: kWh)
    #[serde(default)]
    #[validate(custom = |v| validate_all_sublists_non_empty(v, "SmartApplianceBattery"))]
    #[validate(custom = |v| validate_map_non_empty(v, "SmartApplianceBattery"))]
    pub(crate) energy_into_battery_from_generation: IndexMap<String, Vec<f64>>,

    /// Dictionary of lists containing energy sent to the battery from the grid for each timestep for each energy supply (unit: kWh)
    #[serde(default)]
    #[validate(custom = |v| validate_all_sublists_non_empty(v, "SmartApplianceBattery"))]
    #[validate(custom = |v| validate_map_non_empty(v, "SmartApplianceBattery"))]
    pub(crate) energy_into_battery_from_grid: IndexMap<String, Vec<f64>>,

    /// Dictionary of lists containing energy drawn from the battery for each timestep for each energy supply (unit: kWh)
    #[serde(default)]
    #[validate(custom = |v| validate_all_sublists_non_empty(v, "SmartApplianceBattery"))]
    #[validate(custom = |v| validate_map_non_empty(v, "SmartApplianceBattery"))]
    pub(crate) energy_out_of_battery: IndexMap<String, Vec<f64>>,
}

fn validate_battery_state_fractions(
    data: &IndexMap<String, Vec<f64>>,
) -> Result<(), serde_valid::validation::Error> {
    data
        .values()
        .flatten()
        .all(|&fraction| fraction >= 0. && fraction <= 1.)
        .then_some(())
        .ok_or_else(|| serde_valid::validation::Error::Custom("battery_state_of_charge of SmartApplianceBattery was provided with numbers that were invalid fractions.".to_string()))
}

fn validate_all_sublists_non_empty(
    data: &IndexMap<String, Vec<f64>>,
    struct_type: &str,
) -> Result<(), serde_valid::validation::Error> {
    data
        .values()
        .all(|sublist| !sublist.is_empty())
        .then_some(())
        .ok_or_else(|| serde_valid::validation::Error::Custom(format!("A {struct_type} value had an empty list of numbers provided. Lists must not be empty.")))
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

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct SmartApplianceControlDetails {
    /// List of names of all appliance objects in the simulation
    #[serde(rename = "Appliances")]
    pub(crate) appliances: Vec<ApplianceKey>,

    #[serde(rename = "battery24hr")]
    #[validate]
    pub(crate) battery_24hr: SmartApplianceBattery,

    /// Dictionary of lists containing demand per end user for each timestep for each energy supply (unit: W)
    #[validate(custom = |v| validate_all_sublists_non_empty(v, "SmartApplianceControlDetails"))]
    pub(crate) non_appliance_demand_24hr: IndexMap<String, Vec<f64>>,

    /// Dictionary of lists containing expected power for appliances for each energy supply, for the entire length of the simulation (unit: W)
    #[validate(custom = |v| validate_all_sublists_non_empty(v, "SmartApplianceControlDetails"))]
    #[validate(custom = |v| validate_map_non_empty(v, "SmartApplianceControlDetails"))]
    pub(crate) power_timeseries: IndexMap<String, Vec<f64>>,

    /// Timestep of the power time series (unit: hours)
    #[validate(minimum = 0.)]
    #[validate(maximum = 24.)]
    pub(crate) time_series_step: f64,
}

fn validate_map_non_empty<T>(
    data: &IndexMap<String, T>,
    struct_type: &str,
) -> Result<(), serde_valid::validation::Error> {
    (!data.is_empty()).then_some(()).ok_or_else(|| {
        serde_valid::validation::Error::Custom(format!(
            "A {struct_type} value contained a field of a map (object) type that was empty."
        ))
    })
}

pub(crate) type HotWaterSource = IndexMap<std::string::String, HotWaterSourceDetails>;

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
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: String,

        /// Map of heating systems connected to the storage tank
        #[serde(rename = "HeatSource")]
        #[validate]
        heat_source: IndexMap<std::string::String, HeatSource>,

        /// Measured standby losses due to cylinder insulation at standardised conditions (unit: kWh/24h)
        #[validate(minimum = 0.)]
        daily_losses: f64,

        /// Surface area of the heat exchanger within the storage tank (unit: m²)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        heat_exchanger_surface_area: Option<f64>,

        /// Initial temperature of the storage tank at the start of simulation (unit: ˚C)
        #[validate(minimum = -273.15)]
        init_temp: f64,

        /// List of primary pipework components connected to the storage tank
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate]
        primary_pipework: Option<Vec<WaterPipework>>,

        /// Total volume of tank (unit: litre)
        #[validate(minimum = 0.)]
        volume: f64,
    },
    CombiBoiler {
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: String,

        #[serde(rename = "HeatSourceWet")]
        heat_source_wet: String,

        /// Type of separate domestic hot water test performed on the combi boiler (M&L, M&S, M_only, or No_additional_tests)
        #[serde(rename = "separate_DHW_tests")]
        separate_dhw_tests: BoilerHotWaterTest,

        /// Rejected energy factor 1 for combi boiler efficiency calculations (unit: kWh)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        rejected_energy_1: Option<f64>,

        /// Storage loss factor 2 for combi boiler efficiency calculations (dimensionless)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        storage_loss_factor_2: Option<f64>,

        /// Rejected energy factor 3 for combi boiler efficiency calculations (dimensionless)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        rejected_factor_3: Option<f64>,

        /// Temperature setpoint for the combi boiler hot water output (unit: ˚C)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = -273.15)]
        setpoint_temp: Option<f64>,

        /// Daily hot water usage for the combi boiler system (unit: litre/day)
        #[serde(rename = "daily_HW_usage")]
        #[validate(minimum = 0.)]
        daily_hw_usage: f64,
    },
    #[serde(rename = "HIU")]
    Hiu {
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: String,

        #[serde(rename = "HeatSourceWet")]
        heat_source_wet: String,

        /// Temperature setpoint for the HIU hot water output (unit: ˚C)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = -273.15)]
        setpoint_temp: Option<f64>,
    },
    PointOfUse {
        /// Thermal efficiency of the point-of-use water heater (dimensionless, 0-1)
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        efficiency: Option<f64>,

        #[serde(rename = "EnergySupply")]
        energy_supply: String,

        #[serde(rename = "ColdWaterSource")]
        cold_water_source: String,

        /// Temperature setpoint for the point-of-use water heater output (unit: ˚C)
        #[validate(minimum = -273.15)]
        setpoint_temp: f64,
    },
    SmartHotWaterTank {
        /// Total volume of tank (unit: litre)
        #[validate(minimum = 0.)]
        volume: f64,

        /// Electrical power consumption of the pump (unit: kW)
        #[serde(rename = "power_pump_kW")]
        #[validate(minimum = 0.)]
        power_pump_kw: f64,

        /// Maximum flow rate of the pump (unit: litre/minute)
        #[validate(minimum = 0.)]
        max_flow_rate_pump_l_per_min: f64,

        /// Temperature below which water is considered unusable (unit: ˚C)
        #[validate(minimum = -273.15)]
        temp_usable: f64,

        /// Reference to a control schedule of maximum state of charge values
        temp_setpnt_max: String,

        /// Daily standby losses due to tank insulation at standardised conditions (unit: kWh/24h)
        #[validate(minimum = 0.)]
        daily_losses: f64,

        /// Initial temperature of the smart hot water tank at the start of simulation (unit: ˚C)
        #[validate(minimum = -273.15)]
        init_temp: f64,

        #[serde(rename = "ColdWaterSource")]
        cold_water_source: String,

        #[serde(rename = "EnergySupply_pump")]
        energy_supply_pump: String,

        /// Dictionary of heating systems connected to the smart hot water tank
        #[serde(rename = "HeatSource")]
        #[validate]
        heat_source: IndexMap<std::string::String, HeatSource>,

        /// List of primary pipework components connected to the smart hot water tank
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate]
        primary_pipework: Option<Vec<WaterPipeworkSimple>>,
    },
    HeatBattery {
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: String,

        #[serde(rename = "HeatSourceWet")]
        heat_source_wet: String,

        /// Temperature setpoint for the heat battery hot water output (unit: ˚C)
        #[validate(minimum = -273.15)]
        setpoint_temp: f64,
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
        /// (unit: kW)
        #[validate(minimum = 0.)]
        power: f64,

        #[serde(rename = "EnergySupply")]
        energy_supply: String,

        /// Reference to a control schedule of minimum temperature setpoints
        #[serde(rename = "Controlmin", skip_serializing_if = "Option::is_none")]
        control_min: Option<String>,

        /// Reference to a control schedule of maximum temperature setpoints
        #[serde(rename = "Controlmax", skip_serializing_if = "Option::is_none")]
        control_max: Option<String>,

        /// Vertical position of the heater within the tank, as a fraction of the tank height (0 = bottom, 1 = top). Dimensionless.
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        heater_position: f64,

        /// Vertical position of the thermostat within the tank, as a fraction of the tank height (0 = bottom, 1 = top). Dimensionless. Required for StorageTank but not for SmartHotWaterTank.
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        thermostat_position: Option<f64>,
    },
    SolarThermalSystem {
        /// Location of the main part of the collector loop piping
        #[serde(rename = "sol_loc")]
        solar_cell_location: SolarCollectorLoopLocation,

        /// Collector module reference area (unit: m2)
        #[validate(minimum = 0.)]
        area_module: f64,

        /// Number of collector modules installed
        #[validate(minimum = 1)]
        modules: usize,

        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        peak_collector_efficiency: f64,

        /// Hemispherical incidence angle modifier
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        incidence_angle_modifier: f64,

        /// First order heat loss coefficient
        #[validate(minimum = 0.)]
        first_order_hlc: f64,

        /// Second order heat loss coefficient
        #[validate(minimum = 0.)]
        second_order_hlc: f64,

        /// Mass flow rate solar loop (unit: kg/s)
        #[validate(minimum = 0.)]
        collector_mass_flow_rate: f64,

        /// Power of collector pump (unit: W)
        #[validate(minimum = 0.)]
        power_pump: f64,

        /// Power of collector pump controller (unit: W)
        #[validate(minimum = 0.)]
        power_pump_control: f64,

        #[serde(rename = "EnergySupply")]
        energy_supply: String,

        /// Tilt angle (inclination) of the PV panel from horizontal,
        /// measured upwards facing, 0 to 90, in degrees.
        /// 0=horizontal surface, 90=vertical surface.
        /// Needed to calculate solar irradiation at the panel surface.
        #[validate(minimum = 0.)]
        #[validate(maximum = 90.)]
        tilt: f64,

        #[serde(
            rename = "orientation360",
            deserialize_with = "deserialize_orientation",
            serialize_with = "serialize_orientation"
        )]
        #[validate(minimum = -180.)]
        #[validate(maximum = 180.)]
        orientation: f64,

        /// Heat loss coefficient of the collector loop piping (unit: W/K)
        #[validate(minimum = 0.)]
        solar_loop_piping_hlc: f64,

        /// Vertical position of the heater within the tank, as a fraction of the tank height (0 = bottom, 1 = top). Dimensionless.
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        heater_position: f64,

        /// Vertical position of the thermostat within the tank, as a fraction of the tank height (0 = bottom, 1 = top). Dimensionless. Required for StorageTank but not for SmartHotWaterTank.
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        #[serde(skip_serializing_if = "Option::is_none")]
        thermostat_position: Option<f64>,

        /// Reference to a control schedule of maximum temperature setpoints. References a key in $.Control.
        #[serde(rename = "Controlmax")]
        control_max: String,
    },
    #[serde(rename = "HeatSourceWet")]
    ServiceWaterRegular {
        /// User-defined name for this heat source.
        name: String,

        /// Upper operating limit for flow temperature (unit: °C). Optional.
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = -273.15)]
        temp_flow_limit_upper: Option<f64>,

        #[serde(rename = "EnergySupply")]
        energy_supply: String,

        /// Reference to a control schedule of minimum temperature setpoints
        #[serde(rename = "Controlmin", skip_serializing_if = "Option::is_none")]
        control_min: Option<String>,

        /// Reference to a control schedule of maximum temperature setpoints
        #[serde(rename = "Controlmax", skip_serializing_if = "Option::is_none")]
        control_max: Option<String>,

        /// Vertical position of the heater within the tank, as a fraction of the tank height (0 = bottom, 1 = top). Dimensionless
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        heater_position: f64,

        /// Vertical position of the thermostat within the tank, as a fraction of the tank height (0 = bottom, 1 = top). Dimensionless. Required for StorageTank but not for SmartHotWaterTank.
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        thermostat_position: Option<f64>,
    },
    #[serde(rename = "HeatPump_HWOnly")]
    HeatPumpHotWaterOnly {
        /// (unit: kW)
        #[validate(minimum = 0.)]
        power_max: f64,

        /// Annual average hot water use for the dwelling (unit: litres/day)
        #[validate(exclusive_minimum = 0.)]
        vol_hw_daily_average: f64,

        /// Tank volume stored in the database (unit: litres)
        #[validate(exclusive_minimum = 0.)]
        tank_volume_declared: f64,

        /// Surface area of heat exchanger stored in the database (unit: m2)
        #[validate(minimum = 0.)]
        heat_exchanger_surface_area_declared: f64,

        /// Standing heat loss (unit: kWh/day)
        #[validate(minimum = 0.)]
        daily_losses_declared: f64,

        /// In use factor to be applied to heat pump efficiency
        #[validate(minimum = 0.)]
        in_use_factor_mismatch: f64,

        /// Dictionary with keys denoting tapping profile letter (M or L)
        test_data: HeatPumpHotWaterTestData,

        #[serde(rename = "EnergySupply")]
        energy_supply: String,

        /// Reference to a control schedule of minimum temperature setpoints. References a key in $.Control.
        #[serde(rename = "Controlmin")]
        control_min: String,

        /// Reference to a control schedule of maximum temperature setpoints. References a key in $.Control.
        #[serde(rename = "Controlmax")]
        control_max: String,

        /// Vertical position of the heater within the tank, as a fraction of the tank height (0 = bottom, 1 = top). Dimensionless.
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        heater_position: f64,

        /// Vertical position of the thermostat within the tank, as a fraction of the tank height (0 = bottom, 1 = top). Dimensionless. Required for StorageTank but not for SmartHotWaterTank.
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        thermostat_position: Option<f64>,
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
            } => *thermostat_position,
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

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct HeatPumpHotWaterOnlyTestDatum {
    /// CoP measured during EN 16147 test
    #[validate(minimum = 0.)]
    pub(crate) cop_dhw: f64,

    /// daily energy requirement (kWh/day) for tapping profile used for test
    #[validate(minimum = 0.)]
    pub(crate) hw_tapping_prof_daily_total: f64,

    /// electrical input energy (kWh) measured in EN 16147 test over 24 hrs
    #[validate(minimum = 0.)]
    pub(crate) energy_input_measured: f64,

    /// standby power (W) measured in EN 16147 test
    #[validate(minimum = 0.)]
    pub(crate) power_standby: f64,

    /// daily hot water vessel heat loss
    /// (kWh/day) for a 45 K temperature difference between vessel
    /// and surroundings, tested in accordance with BS 1566 or
    /// EN 12897 or any equivalent standard. Vessel must be same
    /// as that used during EN 16147 test (unit: kWh/day)
    #[validate(minimum = 0.)]
    pub(crate) hw_vessel_loss_daily: f64,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct WaterPipeworkSimple {
    pub location: WaterPipeworkLocation,

    /// (unit: mm)
    #[validate(exclusive_minimum = 0.)]
    pub internal_diameter_mm: f64,

    /// (unit: m)
    #[validate(exclusive_minimum = 0.)]
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
    pub pipe_contents: Option<PipeworkContents>,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct WaterPipework {
    /// Location of the pipework (internal or external)
    pub location: WaterPipeworkLocation,

    /// (unit: mm)
    #[validate(exclusive_minimum = 0.)]
    pub internal_diameter_mm: f64,

    /// (unit: mm)
    #[validate(exclusive_minimum = 0.)]
    pub external_diameter_mm: f64,

    /// (unit: m)
    #[validate(minimum = 0.)]
    pub length: f64,

    /// Thermal conductivity of the insulation (unit: W / m K)
    #[validate(exclusive_minimum = 0.)]
    pub insulation_thermal_conductivity: f64,

    /// (unit: mm)
    #[validate(minimum = 0.)]
    pub insulation_thickness_mm: f64,

    pub surface_reflectivity: bool,

    /// Contents of the pipework (water or glycol25)
    pub pipe_contents: PipeworkContents,
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
pub enum PipeworkContents {
    #[serde(rename = "water")]
    Water,

    #[serde(rename = "glycol25")]
    Glycol25,
}

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub(crate) struct HotWaterDemand {
    #[serde(default, rename = "Shower")]
    #[validate]
    pub(crate) shower: Showers,

    #[serde(default, rename = "Bath")]
    #[validate]
    pub(crate) bath: Baths,

    #[serde(default, rename = "Other")]
    #[validate]
    pub(crate) other_water_use: OtherWaterUses,

    #[serde(default, rename = "Distribution")]
    #[validate]
    pub(crate) water_distribution: WaterDistribution,
}

#[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[validate]
pub struct Showers(#[validate] pub IndexMap<std::string::String, Shower>);

impl Showers {
    /// Provide shower field names as strings.
    pub fn keys(&self) -> Vec<String> {
        self.0.keys().map(Into::into).collect()
    }

    pub fn name_refers_to_instant_electric_shower(&self, name: &str) -> bool {
        self.0
            .get(name)
            .is_some_and(|shower| matches!(shower, Shower::InstantElectricShower { .. }))
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields, tag = "type")]
pub enum Shower {
    MixerShower {
        /// Shower flow rate (unit: litre/minute)
        #[validate(minimum = 0.)]
        flowrate: f64,

        #[serde(rename = "ColdWaterSource")]
        cold_water_source: String,

        /// Reference to HotWaterSource object that provides hot water to this shower. If only one HotWaterSource is defined, then this will be assumed by default
        #[serde(rename = "HotWaterSource", skip_serializing_if = "Option::is_none")]
        hot_water_source: Option<String>,

        #[serde(flatten, skip_serializing_if = "Option::is_none")]
        wwhrs_config: Option<MixerShowerWwhrsConfiguration>,
    },
    #[serde(rename = "InstantElecShower")]
    InstantElectricShower {
        /// Shower's rated electrical power (unit: kW)
        #[validate(minimum = 0.)]
        rated_power: f64,

        #[serde(rename = "ColdWaterSource")]
        cold_water_source: String,

        #[serde(rename = "EnergySupply")]
        energy_supply: String,
    },
}

#[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub struct MixerShowerWwhrsConfiguration {
    /// Reference to a key in Input.WWHRS
    #[serde(rename = "WWHRS")]
    pub(crate) waste_water_heat_recovery_system: String,

    /// WWHRS system configuration for this shower connection
    #[serde(
        rename = "WWHRS_configuration",
        skip_serializing_if = "Option::is_none"
    )]
    pub(crate) wwhrs_configuration: Option<WwhrsConfiguration>,
}

/// WWHRS system configuration
///
///    A - Both shower and water heating system get pre-heated water
///    B - Only shower gets pre-heated water
///    C - Only water heating system gets pre-heated water
#[derive(Clone, Copy, Debug, Default, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum WwhrsConfiguration {
    #[serde(rename = "A")]
    #[default]
    ShowerAndWaterHeatingSystem,

    #[serde(rename = "B")]
    Shower,

    #[serde(rename = "C")]
    WaterHeatingSystem,
}

#[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) struct Baths(#[validate] pub IndexMap<std::string::String, BathDetails>);

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct BathDetails {
    /// Volume held by bath (unit: litre)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) size: f64,

    #[serde(rename = "ColdWaterSource")]
    pub(crate) cold_water_source: String,

    /// Reference to HotWaterSource object that provides hot water to this bath. If only one HotWaterSource is defined, then this will be assumed by default
    #[serde(rename = "HotWaterSource", skip_serializing_if = "Option::is_none")]
    pub(crate) hot_water_source: Option<String>,

    /// Tap/outlet flow rate (unit: litre/minute)
    #[validate(minimum = 0.)]
    pub(crate) flowrate: f64,
}

#[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct OtherWaterUses(#[validate] pub IndexMap<std::string::String, OtherWaterUse>);

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct OtherWaterUse {
    /// Tap/outlet flow rate (unit: litre/minute)
    #[validate(minimum = 0.)]
    pub(crate) flowrate: f64,

    #[serde(rename = "ColdWaterSource")]
    pub(crate) cold_water_source: String,

    /// Reference to HotWaterSource object that provides hot water to this tapping point. If only one HotWaterSource is defined, then this will be assumed by default
    #[serde(rename = "HotWaterSource", skip_serializing_if = "Option::is_none")]
    pub(crate) hot_water_source: Option<String>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(untagged)]
pub(crate) enum WaterDistribution {
    List(#[validate] Vec<WaterPipeworkSimple>),
    Map(#[validate] IndexMap<std::string::String, Vec<WaterPipeworkSimple>>),
}

impl Default for WaterDistribution {
    fn default() -> Self {
        Self::Map(Default::default())
    }
}

#[derive(Clone, Debug, Default, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(rename_all = "PascalCase")]
pub(crate) struct WaterHeatingEvents {
    /// Dictionary of shower water heating events, where keys are shower names and values are lists of events
    #[serde(default)]
    #[validate]
    pub(crate) shower: IndexMap<std::string::String, Vec<WaterHeatingEvent>>,

    /// Dictionary of bath water heating events, where keys are bath names and values are lists of events
    #[serde(default)]
    #[validate]
    pub(crate) bath: IndexMap<std::string::String, Vec<WaterHeatingEvent>>,

    /// Dictionary of other water heating events (e.g., taps, sinks), where keys are event names and values are lists of events
    #[serde(default)]
    #[validate]
    pub(crate) other: IndexMap<std::string::String, Vec<WaterHeatingEvent>>,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub struct WaterHeatingEvent {
    /// Hours from start of simulation to when the water heating event begins (unit: hours)
    #[validate(minimum = 0.)]
    pub start: f64,

    /// Duration of the water heating event (unit: minutes)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = 0.)]
    pub duration: Option<f64>,

    /// Volume of water for the event (unit: litre)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = 0.)]
    pub volume: Option<f64>,

    /// Target temperature for the water heating event (unit: ˚C)
    #[validate(minimum = -273.15)]
    pub temperature: f64,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum WaterHeatingEventType {
    Shower,
    Bath,
    Other,
}

pub(crate) type SpaceHeatSystem = IndexMap<std::string::String, SpaceHeatSystemDetails>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "type")]
pub(crate) enum SpaceHeatSystemDetails {
    #[serde(rename = "InstantElecHeater")]
    InstantElectricHeater {
        //// Rated power of the instant electric heater. (Unit: kW)
        #[validate(minimum = 0.)]
        rated_power: f64,

        #[serde(rename = "EnergySupply")]
        energy_supply: String,

        #[serde(rename = "Control")]
        control: String,

        /// Convective fraction for heating
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        frac_convective: f64,
    },
    #[serde(rename = "ElecStorageHeater")]
    ElectricStorageHeater {
        #[serde(rename = "ControlCharger")]
        control_charger: String,

        /// Maximum output of the electric storage heater. (Unit: kW)
        /// Data from test showing the output from the storage heater when it is actively
        /// outputting heat, e.g. damper open / fan running.
        #[validate(custom = |v| validate_dry_core_output(v, "dry_core_max_output"))]
        dry_core_max_output: Vec<[f64; 2]>,

        /// Minimum output of the electric storage heater. (Unit: kW)
        /// Data from test showing the output from the storage heater when not actively
        /// outputting heat, i.e. case losses only
        #[validate(custom = |v| validate_dry_core_output(v, "dry_core_min_output"))]
        dry_core_min_output: Vec<[f64; 2]>,

        #[serde(rename = "EnergySupply")]
        energy_supply: String,

        air_flow_type: ElectricStorageHeaterAirFlowType,

        #[serde(rename = "Control")]
        control: String,

        /// Fan power (unit: W)
        #[validate(minimum = 0.)]
        fan_pwr: f64,

        /// Convective fraction for heating
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        frac_convective: f64,

        /// Number of units installed in the zone.
        #[validate(exclusive_minimum = 0)]
        n_units: u32,

        /// The rated power of the heating element which charges the storage medium with heat (unit: kW)
        #[validate(minimum = 0.)]
        pwr_in: f64,

        /// State of charge at initialisation of dry core heat storage (ratio)
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        state_of_charge_init: f64,

        /// The rated power output of the instantaneous backup heater (unit: kW)
        #[validate(minimum = 0.)]
        rated_power_instant: f64,

        /// Storage capacity of the electric storage heater. (Unit: kWh)
        #[validate(exclusive_minimum = 0.)]
        storage_capacity: f64,

        /// The zone where the unit(s) is/are installed
        #[serde(rename = "Zone")]
        zone: String,
    },
    WetDistribution {
        #[serde(rename = "HeatSource")]
        #[validate]
        heat_source: SpaceHeatSystemHeatSource,

        /// Fraction of return back into flow water
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        bypass_fraction_recirculated: Option<f64>,

        /// Design flow temperature. (Unit: ˚C)
        #[validate(exclusive_minimum = 0)]
        design_flow_temp: i32,

        /// Wet emitter details of the heating system.
        #[validate(min_items = 1)]
        #[validate]
        emitters: Vec<WetEmitter>,

        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        #[validate]
        pipework: Vec<WaterPipework>,

        ecodesign_controller: EcoDesignController,

        /// Design temperature difference across the emitters. (Unit: deg C or K)
        #[validate(minimum = 0.)]
        temp_diff_emit_dsgn: f64,

        #[serde(skip_serializing_if = "Option::is_none")]
        /// Thermal mass of the emitters. (Unit: kWh/K)
        #[validate(minimum = 0.)]
        thermal_mass: Option<f64>,

        #[serde(rename = "Control")]
        control: String,

        #[serde(rename = "EnergySupply", skip_serializing_if = "Option::is_none")]
        energy_supply: Option<String>,

        /// Zone in which the emitters are located. References a key in $.Zone
        #[serde(rename = "Zone")]
        zone: String,

        #[serde(flatten)]
        #[validate]
        flow_data: FlowData,
    },
    WarmAir {
        /// Convective fraction for heating
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        frac_convective: f64,

        #[serde(rename = "HeatSource")]
        #[validate]
        heat_source: SpaceHeatSystemHeatSource,

        #[serde(rename = "Control")]
        control: String,
    },
}

fn validate_dry_core_output(
    output_data: &[[f64; 2]],
    field: &str,
) -> Result<(), serde_valid::validation::Error> {
    //ensure body of data has at least 2 pairs
    if output_data.len() < 2 {
        return custom_validation_error(format!(
            "The field {field} for an electric storage heater must have at least 2 pairs of data."
        ));
    }

    // Convert ESH_***_output to NumPy arrays without sorting
    let soc_values = output_data.iter().map(|f| f[0]).collect_vec();

    // Validate that SOC array is in strictly increasing order
    if !soc_values.iter().tuple_windows().all(|(a, b)| a <= b) {
        return custom_validation_error(format!(
            "{field} SOC values must be in increasing order (from 0.0 to 1.0)."
        ));
    }

    // Validate that both SOC arrays start at 0.0 and end at 1.0
    if !is_close!(*soc_values.first().unwrap(), 0.) {
        return custom_validation_error(format!(
            "The first SOC value in {field} must be 0.0 (fully discharged)."
        ));
    }

    if !is_close!(*soc_values.last().unwrap(), 1.) {
        return custom_validation_error(format!(
            "The last SOC value in {field} must be 1.0 (fully charged)."
        ));
    }

    Ok(())
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(untagged)]
pub(crate) enum FlowData {
    Variable {
        #[serde(rename = "variable_flow")]
        _variable_flow: MustBe!(true),

        /// Maximum flow rate allowed (unit: litres/min)
        #[validate(minimum = 0.)]
        max_flow_rate: f64,

        /// Minimum flow rate allowed (unit: litres/min)
        #[validate(minimum = 0.)]
        min_flow_rate: f64,
    },
    Design {
        #[serde(rename = "variable_flow")]
        _variable_flow: MustBe!(false),

        /// Constant flow rate if the heat source can't modulate flow rate (unit: l/s)
        #[validate(minimum = 0.)]
        design_flow_rate: f64,
    },
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "wet_emitter_type", rename_all = "lowercase")]
pub(crate) enum WetEmitter {
    Radiator {
        /// Exponent from characteristic equation of emitters (e.g. derived from BS EN 442 tests)
        #[serde(rename = "n")]
        #[validate(exclusive_minimum = 0.)]
        exponent: f64,

        /// Convective fraction for heating
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        frac_convective: f64,

        #[serde(flatten)]
        #[validate]
        constant_data: RadiatorConstantData,
    },
    Ufh {
        /// Equivalent thermal mass per m² of floor area for under-floor heating systems (unit: kJ/m²K)
        #[validate(exclusive_minimum = 0.)]
        equivalent_specific_thermal_mass: f64,

        /// Heat output per m² of floor area for under-floor heating systems (unit: W/m²K)
        #[validate(exclusive_minimum = 0.)]
        system_performance_factor: f64,

        /// (unit: m²)
        #[validate(minimum = 0.)]
        emitter_floor_area: f64,

        /// Convective fraction for heating
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        frac_convective: f64,
    },
    Fancoil {
        /// Number of units of this specification of fancoil in Zone
        #[serde(default = "default_n_units")]
        #[validate(exclusive_minimum = 0)]
        n_units: usize,

        /// Convective fraction for heating
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        frac_convective: f64,

        /// Manufacturer's data for fancoil unit
        #[validate]
        fancoil_test_data: FancoilTestData,
    },
}

/// Enum encapsulating rule for radiators that either the `c` field should be provided, or `c_per_m` and `length` should be
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[serde(untagged)]
pub(crate) enum RadiatorConstantData {
    Constant {
        /// Constant from characteristic equation of emitters (e.g. derived from BS EN 442 tests)
        #[serde(rename = "c")]
        #[validate(exclusive_minimum = 0.)]
        constant: f64,
    },
    ConstantUsingLength {
        /// Constant from characteristic equation of emitters (e.g. derived from BS EN 442 tests) per the length of the emitter
        #[serde(rename = "c_per_m")]
        #[validate(exclusive_minimum = 0.)]
        constant_per_m: f64,

        /// The length of the emitter (unit: m)
        #[validate(exclusive_minimum = 0.)]
        length: f64,
    },
}

impl RadiatorConstantData {
    pub(crate) fn constant(&self) -> f64 {
        match self {
            Self::Constant { constant } => *constant,
            Self::ConstantUsingLength {
                constant_per_m,
                length,
            } => *constant_per_m * *length,
        }
    }
}

const fn default_n_units() -> usize {
    1
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[validate(custom = validate_fancoil_test_data)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) struct FancoilTestData {
    #[validate]
    pub(crate) fan_speed_data: Vec<FanSpeedData>,

    /// A list of fan powers for which heat output data is provided (unit: W)
    #[serde(rename = "fan_power_W")]
    #[validate(custom = validate_all_items_non_negative)]
    pub(crate) fan_power_w: Vec<f64>,
}

fn validate_fancoil_test_data(
    data: &FancoilTestData,
) -> Result<(), serde_valid::validation::Error> {
    // Check all the fan speed lists are of the same length
    if !data
        .fan_speed_data
        .iter()
        .map(|entry| entry.power_output.len())
        .all_equal()
    {
        return custom_validation_error(
            "Fan speed lists of fancoil manufacturer data differ in length".to_string(),
        );
    }

    if data
        .fan_speed_data
        .first()
        .is_none_or(|entry| entry.power_output.len() != data.fan_power_w.len())
    {
        return custom_validation_error(
            "Fan power data length does not match the length of fan speed data".to_string(),
        );
    }

    Ok(())
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) struct FanSpeedData {
    /// Difference in temperature between the hot water supplied to the fan coil and the air in the room (unit: Kelvin)
    #[validate(minimum = 0.)]
    pub(crate) temperature_diff: f64,

    /// Heat output for a specific test temperature difference (unit: kW)
    #[validate(custom = validate_all_items_non_negative)]
    pub(crate) power_output: Vec<f64>,
}

fn validate_all_items_non_negative(items: &[f64]) -> Result<(), serde_valid::validation::Error> {
    items
        .iter()
        .all(|item| item >= &0.)
        .then(|| Ok(()))
        .unwrap_or_else(|| custom_validation_error("All items must be non-negative".to_string()))
}

fn validate_all_items_in_option_non_negative(
    items: &Option<Vec<f64>>,
) -> Result<(), serde_valid::validation::Error> {
    items
        .iter()
        .flatten()
        .all(|item| item >= &0.)
        .then(|| Ok(()))
        .unwrap_or_else(|| custom_validation_error("All items must be non-negative".to_string()))
}

fn validate_all_items_in_option_at_most_n(
    items: &Option<Vec<f64>>,
    n: f64,
) -> Result<(), serde_valid::validation::Error> {
    items
        .iter()
        .flatten()
        .all(|item| item <= &n)
        .then(|| Ok(()))
        .unwrap_or_else(|| custom_validation_error("All items must be non-negative".to_string()))
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct SpaceHeatSystemHeatSource {
    pub(crate) name: String,

    /// Upper operating limit for temperature (unit: deg C)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = 0.)]
    pub(crate) temp_flow_limit_upper: Option<f64>,
}

// it is unclear whether this struct should be used - see reference to the struct above
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[allow(dead_code)]
#[serde(deny_unknown_fields)]
pub(crate) struct EcoDesignController {
    pub(crate) ecodesign_control_class: EcoDesignControllerClass,

    /// Minimum outdoor temperature (unit: Celsius)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) min_outdoor_temp: Option<f64>,

    /// Maximum outdoor temperature (unit: Celsius)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) max_outdoor_temp: Option<f64>,

    /// Minimum flow temperature (unit: Celsius)
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

pub(crate) type ZoneDictionary = IndexMap<std::string::String, ZoneInput>;

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
#[validate(custom = validate_system_list_no_duplicates)]
pub struct ZoneInput {
    /// Map of building elements present in the zone (e.g. walls, floors, windows, etc.).
    #[serde(rename = "BuildingElement")]
    #[validate]
    pub(crate) building_elements: IndexMap<std::string::String, BuildingElement>,

    /// Cooling system details of the zone. References a key in $.SpaceCoolSystem
    #[serde(
        rename = "SpaceCoolSystem",
        skip_serializing_if = "SystemReference::is_none",
        default
    )]
    #[validate]
    pub(crate) space_cool_system: SystemReference,

    /// Heating system details of the zone. References a key in $.SpaceHeatSystem
    #[serde(
        rename = "SpaceHeatSystem",
        skip_serializing_if = "SystemReference::is_none",
        default
    )]
    #[validate]
    pub(crate) space_heat_system: SystemReference,

    /// Overall heat transfer coefficient of the thermal bridge (in W/K), or dictionary of linear thermal transmittance details of the thermal bridges in the zone.
    #[serde(rename = "ThermalBridging")]
    #[validate]
    pub(crate) thermal_bridging: ThermalBridging,

    /// Useful floor area of the zone. (Unit: m²)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) area: f64,

    /// Basis for zone temperature control.
    #[serde(rename = "temp_setpnt_basis", skip_serializing_if = "Option::is_none")]
    pub(crate) temp_setpnt_basis: Option<ZoneTemperatureControlBasis>,

    /// Setpoint temperature to use during initialisation (unit: ˚C)
    #[validate(minimum = -273.15)]
    pub(crate) temp_setpnt_init: f64,

    /// Total volume of the zone. (Unit: m³)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) volume: f64,
}

fn validate_system_list_no_duplicates(
    zone: &ZoneInput,
) -> Result<(), serde_valid::validation::Error> {
    if zone
        .space_heat_system
        .has_intersecting_entries_with(&zone.space_cool_system)
    {
        custom_validation_error(
            "Space heat and space cool systems cannot have overlapping entries".into(),
        )
    } else {
        Ok(())
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(untagged)]
pub(crate) enum SystemReference {
    None(()),
    Single(String),
    Multiple(#[validate(unique_items)] Vec<String>),
}

impl SystemReference {
    fn is_none(&self) -> bool {
        matches!(self, Self::None(_))
    }

    fn has_intersecting_entries_with(&self, other: &Self) -> bool {
        if let (Self::Multiple(a), Self::Multiple(b)) = (self, other) {
            a.iter().any(|a_entry| b.contains(a_entry))
        } else {
            false
        }
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

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(tag = "type")]
#[validate(custom = validate_u_value_and_thermal_resistance_floor_construction)]
pub(crate) enum BuildingElement {
    #[serde(rename = "BuildingElementOpaque")]
    Opaque {
        #[serde(skip_serializing_if = "Option::is_none")]
        is_unheated_pitched_roof: Option<bool>,

        /// Solar absorption coefficient at the external surface (dimensionless)
        solar_absorption_coeff: f64,

        #[serde(flatten)]
        #[validate]
        u_value_input: UValueInput,

        /// Areal heat capacity (unit: J/m².K)
        #[validate(exclusive_minimum = 0.)]
        areal_heat_capacity: f64,

        /// Mass distribution class of the building element, one of: evenly distributed (D); concentrated on external side (E); concentrated on internal side (I); concentrated on internal and external sides (IE); concentrated in middle (M)
        mass_distribution_class: MassDistributionClass,

        is_external_door: Option<bool>,

        /// Tilt angle of the surface from horizontal, between 0 and 180, where 0 means the external surface is facing up, 90 means the external surface is vertical and 180 means the external surface is facing down (unit: ˚)
        #[validate(minimum = 0.)]
        #[validate(maximum = 180.)]
        pitch: f64,

        #[serde(
            rename = "orientation360",
            skip_serializing_if = "Option::is_none",
            deserialize_with = "deserialize_orientation_option",
            serialize_with = "serialize_orientation_option"
        )]
        #[validate(minimum = -180.)]
        #[validate(maximum = 180.)]
        orientation: Option<f64>,

        /// The distance between the ground and the lowest edge of the element (unit: m)
        #[validate(minimum = 0.)]
        base_height: f64,

        #[serde(flatten)]
        #[validate]
        area_input: BuildingElementAreaOrHeightWidthInput,
    },

    #[serde(rename = "BuildingElementTransparent")]
    Transparent {
        #[serde(flatten)]
        #[validate]
        u_value_input: UValueInput,

        /// Areal heat capacity (J/m².K)
        #[validate(exclusive_minimum = 0.)]
        #[serde(default = "default_areal_heat_capacity_for_windows")]
        areal_heat_capacity: f64,

        #[serde(
            rename = "Control_WindowOpenable",
            skip_serializing_if = "Option::is_none"
        )]
        window_openable_control: Option<String>,

        /// Tilt angle of the surface from horizontal, between 0 and 180, where 0 means the external surface is facing up, 90 means the external surface is vertical and 180 means the external surface is facing down (unit: ˚
        #[validate(minimum = 0.)]
        #[validate(maximum = 180.)]
        pitch: f64,

        #[serde(
            rename = "orientation360",
            skip_serializing_if = "Option::is_none",
            deserialize_with = "deserialize_orientation_option",
            serialize_with = "serialize_orientation_option"
        )]
        #[validate(minimum = -180.)]
        #[validate(maximum = 180.)]
        orientation: Option<f64>,

        /// Total solar energy transmittance of the transparent part of the window
        #[validate(minimum = 0.)]
        g_value: f64,

        /// The frame area fraction of window, ratio of the projected frame area to the overall projected area of the glazed element of the window
        frame_area_fraction: f64,

        /// The distance between the ground and the lowest edge of the element (unit: m)
        #[validate(minimum = 0.)]
        base_height: f64,

        #[serde(flatten)]
        #[validate]
        #[validate(custom = validate_area_height_width)]
        area_input: BuildingElementAreaOrHeightWidthInput,

        /// Height of the openable area, corrected for obstruction due to the window frame of the openable section (unit: m)
        #[validate(minimum = 0.)]
        free_area_height: f64,

        /// Height of the mid-point of the window, relative to its base (unit: m)
        #[validate(exclusive_minimum = 0.)]
        mid_height: f64,

        /// Openable area of the window ignoring the obstructing affect of the frame of the openable part
        #[validate(minimum = 0.)]
        max_window_open_area: f64,

        #[validate]
        window_part_list: Vec<WindowPart>,

        #[validate]
        shading: Vec<WindowShadingObject>,

        #[serde(default)]
        #[validate]
        treatment: Vec<WindowTreatment>,
    },

    #[serde(rename = "BuildingElementGround")]
    Ground {
        /// Area of this building element within the zone (unit: m²)
        #[validate(exclusive_minimum = 0.)]
        area: f64,

        /// Total area of the building element across entire dwelling; if the Floor is divided among several zones, this is the total area across all zones (unit: m²)
        #[validate(exclusive_minimum = 0.)]
        total_area: f64,

        /// Tilt angle of the surface from horizontal, between 0 and 180, where 0 means the external surface is facing up, 90 means the external surface is vertical and 180 means the external surface is facing down (unit: ˚)
        #[validate(minimum = 0.)]
        #[validate(maximum = 180.)]
        pitch: f64,

        /// Steady-state thermal transmittance of floor, including the effect of the ground (calculated for the entire ground floor, even if it is distributed among several zones) (unit: W/m2.K)
        #[validate(exclusive_minimum = 0.)]
        u_value: f64,

        /// Total thermal resistance of all layers in the floor construction (unit: m².K/W)
        #[validate(exclusive_minimum = 0.)]
        thermal_resistance_floor_construction: f64,

        #[validate(exclusive_minimum = 0.)]
        /// Areal heat capacity of the ground floor element (unit: J/m2.K)
        #[validate(exclusive_minimum = 0.)]
        areal_heat_capacity: f64,

        /// Mass distribution class of the building element, one of: evenly distributed (D); concentrated on external side (E); concentrated on internal side (I); concentrated on internal and external sides (IE); concentrated in middle (M)
        mass_distribution_class: MassDistributionClass,

        /// Perimeter of the floor; calculated for the entire ground floor, even if it is distributed among several zones (unit: m)
        #[validate(exclusive_minimum = 0.)]
        perimeter: f64,

        /// Linear thermal transmittance of the junction between the floor and the walls (unit: W/m.K)
        psi_wall_floor_junc: f64,

        /// Thickness of the walls (unit: m)
        #[validate(exclusive_minimum = 0.)]
        thickness_walls: f64,

        /// Data specific to types of floors
        #[serde(flatten)]
        #[validate]
        floor_data: FloorData,
    },

    #[serde(rename = "BuildingElementAdjacentConditionedSpace")]
    AdjacentConditionedSpace {
        #[serde(skip_serializing_if = "Option::is_none", default)]
        #[validate(exclusive_minimum = 0.)]
        area: Option<f64>,

        /// Tilt angle of the surface from horizontal, between 0 and 180, where 0 means the external surface is facing up, 90 means the external surface is vertical and 180 means the external surface is facing down (unit: ˚)
        #[validate(minimum = 0.)]
        #[validate(maximum = 180.)]
        pitch: f64,

        #[serde(flatten)]
        #[validate]
        u_value_input: UValueInput,

        /// Areal heat capacity (J/m².K)
        #[validate(exclusive_minimum = 0.)]
        areal_heat_capacity: f64,

        mass_distribution_class: MassDistributionClass,
    },

    #[serde(rename = "BuildingElementAdjacentUnconditionedSpace_Simple")]
    AdjacentUnconditionedSpace {
        /// Area of this building element (unit: m²)
        #[serde(skip_serializing_if = "Option::is_none", default)]
        #[validate(exclusive_minimum = 0.)]
        area: Option<f64>,

        /// Tilt angle of the surface from horizontal, between 0 and 180, where 0 means the external surface is facing up, 90 means the external surface is vertical and 180 means the external surface is facing down (unit: ˚)
        #[validate(minimum = 0.)]
        #[validate(maximum = 180.)]
        pitch: f64,

        #[serde(flatten)]
        #[validate]
        u_value_input: UValueInput,

        /// Effective thermal resistance of unheated space (unit: m².K/W)
        #[validate(exclusive_minimum = 0.)]
        thermal_resistance_unconditioned_space: f64,

        /// Areal heat capacity (unit: J/m2.K)
        #[validate(exclusive_minimum = 0.)]
        areal_heat_capacity: f64,

        mass_distribution_class: MassDistributionClass,
    },
}

const fn default_areal_heat_capacity_for_windows() -> f64 {
    // Windows have much lower thermal mass than walls
    1_000. // J/m².K for typical glazing
}

/// Enum to encapsulate when either u_value or thermal_resistance_construction should be provided, but not both
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(untagged)]
pub(crate) enum UValueInput {
    UValue {
        /// U-value (W/m²·K)
        #[validate(exclusive_minimum = 0.)]
        u_value: f64,
    },
    ThermalResistanceConstruction {
        /// Thermal resistance (m².K/W)
        #[validate(exclusive_minimum = 0.)]
        thermal_resistance_construction: f64,
    },
}

/// Enum to encapsulate when height and width can be provided, or area (or both though they need to match)
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) struct BuildingElementAreaOrHeightWidthInput {
    /// Area of the building element (m²)
    #[validate(exclusive_minimum = 0.)]
    area: Option<f64>,

    #[serde(flatten)]
    pub(crate) height_and_width: Option<BuildingElementHeightWidthInput>,
}

impl BuildingElementAreaOrHeightWidthInput {
    pub(crate) fn area(&self) -> f64 {
        self.area.unwrap_or(self.height_and_width.unwrap_or_else(|| panic!("Building element area/ height and width code has a logic error and could not proceed.")).area())
    }
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) struct BuildingElementHeightWidthInput {
    /// Height of the building element (m)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) height: f64,

    /// Width of the building element (m)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) width: f64,
}

impl BuildingElementHeightWidthInput {
    fn area(&self) -> f64 {
        self.height * self.width
    }
}

fn validate_area_height_width(
    data: &BuildingElementAreaOrHeightWidthInput,
) -> Result<(), serde_valid::validation::Error> {
    match data {
        BuildingElementAreaOrHeightWidthInput {
            area: None,
            height_and_width: None,
        } => custom_validation_error("Building element input needed to specify an area value if a height and width pair of values was not provided.".into()),
        BuildingElementAreaOrHeightWidthInput {
            area: Some(area),
            height_and_width: Some(height_and_width),
        } if !is_close!(*area, height_and_width.area()) => custom_validation_error(
            format!(
                "Building element specified an area of {area} but a height and width pair of {} and {} making an area {} was also provided. These areas need to align.",
                height_and_width.height,
                height_and_width.width,
                height_and_width.area()
            )
        ),
        _ => Ok(())
    }
}

fn validate_u_value_and_thermal_resistance_floor_construction(
    data: &BuildingElement,
) -> Result<(), serde_valid::validation::Error> {
    // this validation only applies for ground building elements
    if let BuildingElement::Ground {
        u_value,
        thermal_resistance_floor_construction,
        ..
    } = data
    {
        return match calculate_thermal_resistance_of_virtual_layer(*u_value, *thermal_resistance_floor_construction) {
            Ok(_) => Ok(()),
            Err(_) => custom_validation_error(format!("Ground building element provided with u_value {u_value} and thermal_resistance_floor_construction {thermal_resistance_floor_construction} values that were not compatible with each other."))
        };
    }

    Ok(())
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

    // pub(crate) fn height(&self) -> Option<f64> {
    //     match self {
    //         BuildingElement::Opaque { height, .. } => Some(*height),
    //         BuildingElement::Transparent { height, .. } => Some(*height),
    //         _ => None,
    //     }
    // }
    //
    // pub(crate) fn width(&self) -> Option<f64> {
    //     match self {
    //         BuildingElement::Opaque { width, .. } => Some(*width),
    //         BuildingElement::Transparent { width, .. } => Some(*width),
    //         _ => None,
    //     }
    // }

    #[cfg(test)]
    pub(crate) fn u_value(&self) -> Option<f64> {
        match self {
            BuildingElement::Opaque { u_value_input, .. } => {
                if let UValueInput::UValue { u_value } = u_value_input {
                    Some(*u_value)
                } else {
                    None
                }
            }
            BuildingElement::Transparent { u_value_input, .. } => {
                if let UValueInput::UValue { u_value } = u_value_input {
                    Some(*u_value)
                } else {
                    None
                }
            }
            BuildingElement::Ground { u_value, .. } => Some(*u_value),
            BuildingElement::AdjacentConditionedSpace { u_value_input, .. } => {
                if let UValueInput::UValue { u_value } = u_value_input {
                    Some(*u_value)
                } else {
                    None
                }
            }
            BuildingElement::AdjacentUnconditionedSpace { u_value_input, .. } => {
                if let UValueInput::UValue { u_value } = u_value_input {
                    Some(*u_value)
                } else {
                    None
                }
            }
        }
    }

    pub(crate) fn orientation(&self) -> Option<f64> {
        match self {
            BuildingElement::Opaque { orientation, .. } => *orientation,
            BuildingElement::Transparent { orientation, .. } => *orientation,
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
        self.0.shift_remove("thermal_resistance_construction");
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

// special deserialization logic so that orientations are normalized correctly on the way in
pub(crate) fn deserialize_orientation_option<'de, D>(
    deserializer: D,
) -> Result<Option<f64>, D::Error>
where
    D: Deserializer<'de>,
{
    let orientation360_value: Option<f64> = Deserialize::deserialize(deserializer)?;
    Ok(orientation360_value.map(init_orientation))
}

// special serialization logic so that orientations are un-normalized correctly on way out!
pub(crate) fn serialize_orientation_option<S>(
    normalized: &Option<f64>,
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let orientation360_value = normalized.map(init_orientation);
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

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct WindowPart {
    /// (unit: m)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) mid_height_air_flow_path: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct WindowTreatment {
    #[serde(rename = "type")]
    pub(crate) treatment_type: WindowTreatmentType,

    pub(crate) controls: WindowTreatmentControl,

    /// Additional thermal resistance provided by a window treatment (unit: m²K/W)
    pub(crate) delta_r: f64,

    /// Dimensionless factor describing the reduction in the amount of transmitted radiation due to a window treatment
    #[validate(minimum = 0.)]
    #[validate(maximum = 1.)]
    pub(crate) trans_red: f64,

    /// Irradiation level above which the window treatment is assumed to be closed (unit: W/m²). References a key in $.Control.
    #[serde(
        rename = "Control_closing_irrad",
        skip_serializing_if = "Option::is_none"
    )]
    pub(crate) closing_irradiance_control: Option<String>,

    /// Irradiation level below which a window treatment is assumed to be open (unit: W/m²). References a key in $.Control.
    #[serde(
        rename = "Control_opening_irrad",
        skip_serializing_if = "Option::is_none"
    )]
    pub(crate) opening_irradiance_control: Option<String>,

    /// Reference to a time control object containing a schedule of booleans describing when a window treatment is open. References a key in $.Control.
    #[serde(rename = "Control_open", skip_serializing_if = "Option::is_none")]
    pub(crate) open_control: Option<String>,

    /// A boolean describing the state of the window treatment
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub(crate) is_open: Option<bool>,

    /// Time delay enforced before a window treatment may be opened after the conditions for its opening are met (unit: hours)
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

impl WindowTreatmentControl {
    fn is_manual(&self) -> bool {
        matches!(self, Self::Manual | Self::ManualMotorised)
    }

    fn is_automatic(&self) -> bool {
        matches!(self, Self::AutoMotorised | Self::CombinedLightBlindHvac)
    }
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
        #[validate]
        edge_insulation: Vec<EdgeInsulation>,
    },

    // (non-optional fields in the schema for SuspendedFloor, but values do seem expected)
    #[serde(rename = "Suspended_floor")]
    SuspendedFloor {
        #[validate(exclusive_minimum = 0.)]
        height_upper_surface: f64,

        /// Thermal transmittance of walls above ground (unit: W/m².K)
        #[serde(rename = "thermal_transm_walls")]
        thermal_transmission_walls: f64,

        /// Area of ventilation openings per perimeter (unit: m²/m)
        #[validate(exclusive_minimum = 0.)]
        area_per_perimeter_vent: f64,

        /// Wind shielding factor
        shield_fact_location: WindShieldLocation,

        /// Thermal resistance of insulation on base of underfloor space (unit: m².K/W)
        #[serde(rename = "thermal_resist_insul")]
        #[validate(exclusive_minimum = 0.)]
        thermal_resistance_of_insulation: f64,
    },

    #[serde(rename = "Heated_basement")]
    HeatedBasement {
        /// Depth of basement floor below ground level (unit: m)
        #[validate(minimum = 0.)]
        depth_basement_floor: f64,

        /// Thermal resistance of walls of the basement (unit: m².K/W)
        #[serde(rename = "thermal_resist_walls_base")]
        #[validate(exclusive_minimum = 0.)]
        thermal_resistance_of_basement_walls: f64,
    },

    #[serde(rename = "Unheated_basement")]
    UnheatedBasement {
        /// Thermal transmittance of floor above basement (unit: W/m².K)
        #[serde(rename = "thermal_transm_envi_base")]
        #[validate(exclusive_minimum = 0.)]
        thermal_transmittance_of_floor_above_basement: f64,

        /// Thermal transmittance of walls above ground (unit: W/m².K)
        #[serde(rename = "thermal_transm_walls")]
        #[validate(exclusive_minimum = 0.)]
        thermal_transmission_walls: f64,

        /// Depth of basement floor below ground level (unit: m)
        #[validate(minimum = 0.)]
        depth_basement_floor: f64,

        /// Height of the basement walls above ground level (unit: m)
        #[validate(exclusive_minimum = 0.)]
        height_basement_walls: f64,

        /// Thermal resistance of walls of the basement (unit: m².K/W)
        #[serde(rename = "thermal_resist_walls_base")]
        #[validate(exclusive_minimum = 0.)]
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

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "type")]
pub enum EdgeInsulation {
    #[serde(rename = "horizontal")]
    Horizontal {
        /// (unit: m)
        #[validate(minimum = 0.)]
        width: f64,

        /// Thermal resistance of floor edge insulation (unit: m²K/W)
        #[validate(exclusive_minimum = 0.)]
        edge_thermal_resistance: f64,
    },
    #[serde(rename = "vertical")]
    Vertical {
        /// (unit: m)
        #[validate(minimum = 0.)]
        depth: f64,

        /// Thermal resistance of floor edge insulation (unit: m²K/W)
        #[validate(exclusive_minimum = 0.)]
        edge_thermal_resistance: f64,
    },
}

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(untagged)]
pub enum ThermalBridging {
    Elements(#[validate] IndexMap<std::string::String, ThermalBridgingDetails>),
    Number(f64),
}

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(tag = "type")] // TODO: possibly restore `deny_unknown_fields` serde annotation after 0.36 (once FHS is extracted)
pub enum ThermalBridgingDetails {
    #[serde(rename = "ThermalBridgeLinear")]
    Linear {
        /// Linear thermal transmittance of the thermal bridge. (Unit: W/m.K)
        // NB. this can be positive or negative, so no validation of zero bounds necessary.
        linear_thermal_transmittance: f64,

        /// Length of the thermal bridge over which the linear thermal transmittance applies. (Unit: m)
        #[validate(exclusive_minimum = 0.)]
        length: f64,
    },
    #[serde(rename = "ThermalBridgePoint")]
    Point {
        /// Heat transfer coefficient of the thermal bridge. (Unit: W/K)
        // NB. this can be positive or negative. SAP 2012, Appendix K has negative point thermal bridges.
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

pub(crate) type SpaceCoolSystem = IndexMap<std::string::String, SpaceCoolSystemDetails>;

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(tag = "type")] // TODO: possibly restore `deny_unknown_fields` serde annotation after 0.36 (once FHS is extracted)
pub(crate) enum SpaceCoolSystemDetails {
    AirConditioning {
        /// Maximum cooling capacity of the system (unit: kW)
        #[validate(minimum = 0.)]
        cooling_capacity: f64,

        /// Efficiency of the air conditioning system. SEER (Seasonal energy efficiency ratio)
        #[validate(minimum = 0.)]
        efficiency: f64,

        /// Convective fraction for cooling
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
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

pub(crate) type HeatSourceWet = IndexMap<std::string::String, HeatSourceWetDetails>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[allow(clippy::large_enum_variant)]
#[serde(tag = "type")]
pub(crate) enum HeatSourceWetDetails {
    HeatPump {
        /// Optional buffer tank configuration for the heat pump system
        #[serde(rename = "BufferTank", skip_serializing_if = "Option::is_none")]
        #[validate(custom = validate_boxed_in_option)]
        buffer_tank: Option<Box<HeatPumpBufferTank>>,

        #[serde(rename = "EnergySupply")]
        energy_supply: String,

        /// References a key in $.EnergySupply for heat network energy supply
        #[serde(
            rename = "EnergySupply_heat_network",
            skip_serializing_if = "Option::is_none"
        )]
        energy_supply_heat_network: Option<String>,

        /// References a key in $.MechanicalVentilation
        #[serde(
            rename = "MechanicalVentilation",
            skip_serializing_if = "Option::is_none"
        )]
        mechanical_ventilation: Option<String>,

        /// Type of backup control for the heat pump system
        #[serde(rename = "backup_ctrl_type")]
        backup_control_type: HeatPumpBackupControlType,

        /// Optional boiler configuration used as backup for the heat pump
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(custom = validate_boxed_in_option)]
        boiler: Option<Box<HeatPumpBoiler>>,

        /// Maximum temperature for exhaust air heat pump mixed operation (unit: ˚C)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = -273.15)]
        eahp_mixed_max_temp: Option<f64>,

        /// Minimum temperature for exhaust air heat pump mixed operation (unit: ˚C)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = -273.15)]
        eahp_mixed_min_temp: Option<f64>,

        /// Minimum modulation rate at 20°C flow temperature (dimensionless, 0-1)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        min_modulation_rate_20: Option<f64>,

        /// Minimum modulation rate at 35°C flow temperature (dimensionless, 0-1)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        min_modulation_rate_35: Option<f64>,

        /// Minimum modulation rate at 55°C flow temperature (dimensionless, 0-1)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        min_modulation_rate_55: Option<f64>,

        /// Minimum temperature difference between flow and return for heat pump operation (unit: K)
        min_temp_diff_flow_return_for_hp_to_operate: f64,

        /// Whether the heat pump uses modulating control
        modulating_control: bool,

        /// Power consumption of crankcase heater (unit: kW)
        #[validate(minimum = 0.)]
        power_crankcase_heater: f64,

        /// Power consumption of heating circuit pump (unit: kW)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        power_heating_circ_pump: Option<f64>,

        /// Power consumption of warm air fan (unit: kW)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        power_heating_warm_air_fan: Option<f64>,

        /// Maximum backup power (unit: kW)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        power_max_backup: Option<f64>,

        /// Power consumption when heat pump is off (unit: kW)
        #[validate(minimum = 0.)]
        power_off: f64,

        /// Power consumption of source circuit pump (unit: kW)
        #[validate(minimum = 0.)]
        power_source_circ_pump: f64,

        /// Power consumption in standby mode (unit: kW)
        #[validate(minimum = 0.)]
        power_standby: f64,

        /// Type of heat sink for the heat pump
        sink_type: HeatPumpSinkType,

        /// Type of heat source for the heat pump
        source_type: HeatPumpSourceType,

        /// Distribution temperature for heat network (unit: ˚C)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = -273.15)]
        temp_distribution_heat_network: Option<f64>,

        /// Lower temperature limit for heat pump operation (unit: ˚C)
        #[validate(minimum = -273.15)]
        temp_lower_operating_limit: f64,

        /// Maximum return feed temperature (unit: ˚C)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = -273.15)]
        temp_return_feed_max: Option<f64>,

        #[serde(rename = "test_data_EN14825")]
        #[validate]
        test_data_en14825: Vec<HeatPumpTestDatum>,

        /// Time constant for on/off operation (unit: hours)
        #[validate(minimum = 0.)]
        time_constant_onoff_operation: f64,

        /// Time delay before backup operation (unit: hours)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        time_delay_backup: Option<f64>,

        /// Whether variable flow temperature control was used during testing
        var_flow_temp_ctrl_during_test: bool,
    },
    Boiler {
        #[serde(rename = "EnergySupply")]
        energy_supply: String,

        /// References a key in $.EnergySupply for auxiliary electrical power
        #[serde(rename = "EnergySupply_aux")]
        energy_supply_aux: String,

        #[validate(exclusive_minimum = 0.)]
        rated_power: f64,

        /// Boiler efficiency at full load (dimensionless, 0-1)
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        efficiency_full_load: f64,

        /// Boiler efficiency at part load (dimensionless, 0-1)
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        efficiency_part_load: f64,

        /// Location of the boiler (internal or external to the building)
        boiler_location: HeatSourceLocation,

        /// Modulation load ratio (dimensionless, 0-1)
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        modulation_load: f64,

        /// Electrical power consumption of circulation pump (unit: kW)
        #[validate(minimum = 0.)]
        electricity_circ_pump: f64,

        /// Electrical power consumption at part load (unit: kW)
        #[validate(minimum = 0.)]
        electricity_part_load: f64,

        /// Electrical power consumption at full load (unit: kW)
        #[validate(minimum = 0.)]
        electricity_full_load: f64,

        /// Electrical power consumption in standby mode (unit: kW)
        #[validate(minimum = 0.)]
        electricity_standby: f64,
    },
    HeatBattery {
        #[serde(flatten)]
        #[validate]
        battery: HeatBattery,
    },
    #[serde(rename = "HIU")]
    Hiu {
        #[serde(rename = "EnergySupply")]
        energy_supply: String,

        /// Maximum power output of the HIU (unit: kW)
        #[validate(exclusive_minimum = 0.)]
        power_max: f64,

        /// Daily heat losses from the HIU (unit: kWh/day)
        #[serde(rename = "HIU_daily_loss")]
        #[validate(minimum = 0.)]
        hiu_daily_loss: f64,

        /// Heat losses from building-level distribution pipework (unit: W)
        #[validate(minimum = 0.)]
        building_level_distribution_losses: f64,

        /// Power consumption of heating circulation pump (unit: kW)
        #[validate(minimum = 0.)]
        power_circ_pump: Option<f64>,

        /// Power consumption of auxiliary electrical usage (unit: kW)
        #[validate(minimum = 0.)]
        power_aux: Option<f64>,
    },
}

fn validate_boxed_in_option<T: Validate>(
    value: &Option<Box<T>>,
) -> Result<(), serde_valid::validation::Error> {
    if let Some(boxed) = value.as_ref().map(|v| v.as_ref()) {
        boxed
            .validate()
            .map_err(|e| serde_valid::validation::Error::Custom(format!("Validation error: {e}")))
    } else {
        Ok(())
    }
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

#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct HeatPumpBufferTank {
    /// Standing heat loss (unit: kWh/day)
    #[validate(minimum = 0.)]
    pub daily_losses: f64,

    /// Volume of the buffer tank (unit: litre)
    #[validate(minimum = 0.)]
    pub volume: f64,

    /// Flow rate of the buffer tank - emitters loop (unit: l/min)
    #[validate(minimum = 0.)]
    pub pump_fixed_flow_rate: f64,

    /// Pump power of the buffer tank - emitters loop (unit: W)
    #[validate(minimum = 0.)]
    pub pump_power_at_flow_rate: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Validate, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct HeatPumpTestDatum {
    /// Air flow rate through the heat pump for the test condition (unit: m³/h)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = 0.)]
    pub(crate) air_flow_rate: Option<f64>,

    /// Heat output capacity at this test condition (unit: kW)
    #[validate(minimum = 0.)]
    pub(crate) capacity: f64,

    /// Coefficient of performance at this test condition (dimensionless)
    #[validate(minimum = 0.)]
    pub(crate) cop: f64,

    /// Design flow temperature for the heating system (unit: Celsius)
    pub(crate) design_flow_temp: f64,

    /// Ratio of external air to recirculated air for exhaust air heat pumps (dimensionless)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = 0.)]
    #[validate(maximum = 1.)]
    pub(crate) eahp_mixed_ext_air_ratio: Option<f64>,

    /// Heat pump outlet temperature for the test condition (unit: Celsius)
    #[validate(minimum = -273.15)]
    pub(crate) temp_outlet: f64,

    /// Heat pump source temperature for the test condition (unit: Celsius)
    #[validate(minimum = -273.15)]
    pub(crate) temp_source: f64,

    /// Ambient air temperature for the test condition (unit: Celsius)
    #[validate(minimum = -273.15)]
    pub(crate) temp_test: f64,

    pub(crate) test_letter: TestLetter,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct HeatPumpBoiler {
    #[serde(rename = "EnergySupply")]
    pub(crate) energy_supply: String,

    /// References a key in $.EnergySupply for auxiliary electrical power
    #[serde(rename = "EnergySupply_aux")]
    pub(crate) energy_supply_aux: String,

    /// Rated power output of the boiler (unit: kW)
    #[validate(exclusive_minimum = 0.)]
    rated_power: f64,

    /// Boiler efficiency at full load (dimensionless, 0-1)
    #[validate(minimum = 0.)]
    #[validate(maximum = 1.)]
    efficiency_full_load: f64,

    /// Boiler efficiency at part load (dimensionless, 0-1)
    #[validate(minimum = 0.)]
    #[validate(maximum = 1.)]
    efficiency_part_load: f64,

    /// Location of the boiler (internal or external to the building)
    boiler_location: HeatSourceLocation,

    /// Modulation load ratio (dimensionless, 0-1)
    #[validate(minimum = 0.)]
    #[validate(maximum = 1.)]
    modulation_load: f64,

    /// Electrical power consumption of circulation pump (unit: kW)
    #[validate(minimum = 0.)]
    electricity_circ_pump: f64,

    /// Electrical power consumption at part load (unit: kW)
    #[validate(minimum = 0.)]
    electricity_part_load: f64,

    /// Electrical power consumption at full load (unit: kW)
    #[validate(minimum = 0.)]
    electricity_full_load: f64,

    /// Electrical power consumption in standby mode (unit: kW)
    #[validate(minimum = 0.)]
    electricity_standby: f64,

    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate]
    pub(crate) cost_schedule_hybrid: Option<BoilerCostScheduleHybrid>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct BoilerCostScheduleHybrid {
    /// Cost data for the fuel used by the hybrid's heat pump (can be any units, typically p/kWh)
    pub(crate) cost_schedule_boiler: NumericSchedule,

    /// Cost data for the fuel used by the hybrid's boiler (can be any units, typically p/kWh)
    pub(crate) cost_schedule_hp: NumericSchedule,

    /// Day on which the cost data series begins
    #[validate(minimum = 0)]
    #[validate(maximum = 365)]
    pub(crate) cost_schedule_start_day: u32,

    /// Time step of the cost data series
    pub(crate) cost_schedule_time_series_step: f64,
}

#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum HeatSourceLocation {
    #[serde(rename = "internal")]
    Internal,
    #[serde(rename = "external")]
    External,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "battery_type")]
pub(crate) enum HeatBattery {
    #[serde(rename = "pcm")]
    Pcm {
        #[serde(rename = "EnergySupply")]
        energy_supply: String,

        /// Initial temperature of the PCM heat battery at the start of simulation (unit: ˚C)
        #[validate(minimum = -273.15)]
        temp_init: f64,

        /// Electrical power consumption of circulation pump (unit: kW)
        #[validate(minimum = 0.)]
        electricity_circ_pump: f64,

        /// Electrical power consumption in standby mode (unit: kW)
        #[validate(minimum = 0.)]
        electricity_standby: f64,

        /// Rated charging power (unit: kW)
        #[validate(exclusive_minimum = 0.)]
        rated_charge_power: f64,

        /// Maximum rated heat losses (unit: kW)
        #[validate(minimum = 0.)]
        max_rated_losses: f64,

        /// Number of heat battery units
        number_of_units: usize,

        /// Whether the heat battery can charge and discharge simultaneously
        simultaneous_charging_and_discharging: bool,

        #[serde(rename = "ControlCharge")]
        control_charge: String,

        /// Heat capacity of storage above phase transition (unit: kJ/K)
        #[serde(rename = "heat_storage_kJ_per_K_above_Phase_transition")]
        #[validate(exclusive_minimum = 0.)]
        heat_storage_k_j_per_k_above_phase_transition: f64,

        /// Heat capacity of storage below phase transition (unit: kJ/K)
        #[serde(rename = "heat_storage_kJ_per_K_below_Phase_transition")]
        #[validate(exclusive_minimum = 0.)]
        heat_storage_k_j_per_k_below_phase_transition: f64,

        /// Heat capacity of storage during phase transition (unit: kJ/K)
        #[serde(rename = "heat_storage_kJ_per_K_during_Phase_transition")]
        #[validate(exclusive_minimum = 0.)]
        heat_storage_k_j_per_k_during_phase_transition: f64,

        /// Upper temperature limit for phase transition (unit: ˚C)
        #[validate(minimum = -273.15)]
        phase_transition_temperature_upper: f64,

        /// Lower temperature limit for phase transition (unit: ˚C)
        #[validate(minimum = -273.15)]
        phase_transition_temperature_lower: f64,

        /// Maximum operating temperature (unit: ˚C)
        #[validate(minimum = -273.15)]
        max_temperature: f64,

        /// Velocity in heat exchanger tube at 1 litre/minute flow rate (unit: m/s)
        #[serde(rename = "velocity_in_HEX_tube_at_1_l_per_min_m_per_s")]
        #[validate(exclusive_minimum = 0.)]
        velocity_in_hex_tube_at_1_l_per_min_m_per_s: f64,

        /// Diameter of capillary tubes (unit: m)
        #[validate(exclusive_minimum = 0.)]
        capillary_diameter_m: f64,

        /// Heat battery parameter A (dimensionless)
        #[serde(rename = "A")]
        a: f64,

        /// Heat battery parameter B (dimensionless)
        #[serde(rename = "B")]
        b: f64,

        /// Surface area of heat exchanger (unit: m²)
        #[validate(exclusive_minimum = 0.)]
        heat_exchanger_surface_area_m2: f64,

        /// Flow rate through the heat battery (unit: litre/minute)
        #[validate(exclusive_minimum = 0.)]
        flow_rate_l_per_min: f64,
    },
    #[serde(rename = "dry_core")]
    DryCore {
        #[serde(rename = "ControlCharge")]
        control_charge: String,

        #[serde(rename = "EnergySupply")]
        energy_supply: String,

        electricity_circ_pump: f64,

        electricity_standby: f64,

        /// Charging power (kW)
        pwr_in: f64,

        /// State of charge at initialisation of dry core heat storage (ratio)
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        state_of_charge_init: f64,

        /// Rated instantaneous power output (kW)
        rated_power_instant: f64,

        heat_storage_capacity: f64,

        #[validate(minimum = 1)]
        number_of_units: usize,

        /// Lookup table for minimum output based on charge level
        dry_core_min_output: Vec<[f64; 2]>,

        /// Lookup table for maximum output based on charge level
        dry_core_max_output: Vec<[f64; 2]>,

        /// Fan power (W)
        #[serde(rename = "fan_pwr")]
        fan_power: f64,
    },
}

pub(crate) type WasteWaterHeatRecovery =
    IndexMap<std::string::String, WasteWaterHeatRecoveryDetails>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
#[validate(custom = validate_flow_rates_and_efficiencies_length)]
pub(crate) struct WasteWaterHeatRecoveryDetails {
    #[serde(rename = "type")]
    _type: MustBe!("WWHRS_Instantaneous"),

    #[serde(rename = "ColdWaterSource")]
    pub(crate) cold_water_source: String,

    /// Test flow rates in litres per minute (e.g., [5., 7., 9., 11., 13.])
    #[validate(custom = validate_all_items_non_negative)]
    pub(crate) flow_rates: Vec<f64>,

    /// Measured efficiencies for System A at the test flow rates
    #[validate(custom = validate_all_items_non_negative)]
    pub(crate) system_a_efficiencies: Vec<f64>,

    /// Utilisation factor for System A
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = 0.)]
    #[validate(maximum = 1.)]
    pub(crate) system_a_utilisation_factor: Option<f64>,

    /// Measured efficiencies for System B (optional, uses system_b_efficiency_factor if not provided)
    #[validate(custom = validate_all_items_in_option_non_negative)]
    pub(crate) system_b_efficiencies: Option<Vec<f64>>,

    /// Utilisation factor for System B. Required when using either system_b_efficiencies (pre-corrected data) or when converting system_a_efficiencies to System B (used with system_b_efficiency_factor).
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = 0.)]
    #[validate(maximum = 1.)]
    pub(crate) system_b_utilisation_factor: Option<f64>,

    /// Measured efficiencies for System C (optional, uses system_c_efficiency_factor if not provided)
    #[validate(custom = validate_all_items_in_option_non_negative)]
    pub(crate) system_c_efficiencies: Option<Vec<f64>>,

    /// Utilisation factor for System C. Required when using either system_c_efficiencies (pre-corrected data) or when converting system_a_efficiencies to System C (used with system_c_efficiency_factor).
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = 0.)]
    #[validate(maximum = 1.)]
    pub(crate) system_c_utilisation_factor: Option<f64>,

    /// Reduction factor for converting System A efficiency data to System B (default 0.81). Only used when system_b_efficiencies is not provided.
    #[serde(default = "default_system_b_efficiency_factor")]
    #[validate(minimum = 0.)]
    #[validate(maximum = 1.)]
    pub(crate) system_b_efficiency_factor: f64,

    /// Reduction factor for converting System A efficiency data to System C (default 0.88). Only used when system_c_efficiencies is not provided.
    #[serde(default = "default_system_c_efficiency_factor")]
    #[validate(minimum = 0.)]
    #[validate(maximum = 1.)]
    pub(crate) system_c_efficiency_factor: f64,
}

const fn default_system_b_efficiency_factor() -> f64 {
    0.81
}

const fn default_system_c_efficiency_factor() -> f64 {
    0.88
}

fn validate_flow_rates_and_efficiencies_length(
    data: &WasteWaterHeatRecoveryDetails,
) -> Result<(), serde_valid::validation::Error> {
    // Validate system A (always required)
    if data.flow_rates.len() != data.system_a_efficiencies.len() {
        return custom_validation_error(
            "flow_rates and system_a_efficiencies must have the same length".into(),
        );
    }

    // Validate system B if provided
    if data
        .system_b_efficiencies
        .as_ref()
        .is_some_and(|efficiencies| data.flow_rates.len() != efficiencies.len())
    {
        return custom_validation_error(
            "flow_rates and system_b_efficiencies must have the same length".into(),
        );
    }

    // Validate system C if provided
    if data
        .system_c_efficiencies
        .as_ref()
        .is_some_and(|efficiencies| data.flow_rates.len() != efficiencies.len())
    {
        return custom_validation_error(
            "flow_rates and system_c_efficiencies must have the same length".into(),
        );
    }

    Ok(())
}

fn custom_validation_error(
    message: std::string::String,
) -> Result<(), serde_valid::validation::Error> {
    Err(serde_valid::validation::Error::Custom(message))
}

pub(crate) type OnSiteGeneration = IndexMap<std::string::String, PhotovoltaicInputs>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(untagged)]
pub(crate) enum PhotovoltaicInputs {
    DeprecatedStyle(#[validate] PhotovoltaicSystem),
    WithPanels(#[validate] PhotovoltaicSystemWithPanels),
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) struct PhotovoltaicPanel {
    /// Peak power; represents the electrical power of a photovoltaic system
    /// with a given area for a solar irradiance of 1 kW/m² on this surface (at 25 degrees)
    /// (unit: kW)
    #[validate(minimum = 0.)]
    peak_power: f64,

    ventilation_strategy: PhotovoltaicVentilationStrategy,

    /// The tilt angle (inclination) of the PV panel from horizontal, measured upwards facing, 0 to 90 (unit: ˚)
    #[validate(minimum = 0.)]
    #[validate(maximum = 90.)]
    pitch: f64,

    #[serde(
        rename = "orientation360",
        deserialize_with = "deserialize_orientation",
        serialize_with = "serialize_orientation"
    )]
    #[validate(minimum = -180.)]
    #[validate(maximum = 180.)]
    orientation: f64,

    /// The distance between the ground and the lowest edge of the PV array (unit: m)
    #[validate(minimum = 0.)]
    base_height: f64,

    /// Height of the PV array (unit: m)
    height: f64,

    /// Width of the PV panel (unit: m)
    #[validate(exclusive_minimum = 0.)]
    width: f64,

    #[validate]
    shading: Vec<WindowShadingObject>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) struct PhotovoltaicSystemWithPanels {
    #[serde(rename = "type")]
    _type: MustBe!("PhotovoltaicSystem"),

    #[serde(rename = "EnergySupply")]
    energy_supply: String,

    /// Whether the inverter is considered inside the building
    inverter_is_inside: bool,

    /// Peak power; represents the peak electrical AC power output from the inverter (unit: kW)
    #[validate(minimum = 0.)]
    inverter_peak_power_ac: f64,

    /// Peak power; represents the peak electrical DC power input to the inverter (unit: kW)
    #[validate(minimum = 0.)]
    inverter_peak_power_dc: f64,

    inverter_type: InverterType,

    #[validate(min_items = 1)]
    panels: Vec<PhotovoltaicPanel>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct PhotovoltaicSystem {
    #[serde(rename = "type")]
    _type: MustBe!("PhotovoltaicSystem"),

    /// Peak power; represents the electrical power of a photovoltaic system with a given area for a solar irradiance of 1 kW/m² on this surface (at 25 degrees) (unit: kW)
    #[validate(minimum = 0.)]
    pub(crate) peak_power: f64,

    pub(crate) ventilation_strategy: PhotovoltaicVentilationStrategy,

    /// The tilt angle (inclination) of the PV panel from horizontal, measured upwards facing, 0 to 90 (unit: ˚)
    #[validate(minimum = 0.)]
    #[validate(maximum = 90.)]
    pub(crate) pitch: f64,

    #[serde(
        rename = "orientation360",
        deserialize_with = "deserialize_orientation",
        serialize_with = "serialize_orientation"
    )]
    #[validate(minimum = -180.)]
    #[validate(maximum = 180.)]
    pub(crate) orientation: f64,

    /// The distance between the ground and the lowest edge of the PV array (unit: m)
    #[validate(minimum = 0.)]
    pub(crate) base_height: f64,

    /// Height of the PV array (unit: m)
    #[validate(minimum = 0.)]
    pub(crate) height: f64,

    /// Width of the PV panel (unit: m)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) width: f64,

    #[serde(rename = "EnergySupply")]
    pub(crate) energy_supply: String,

    pub(crate) shading: Vec<WindowShadingObject>,

    /// Peak power; represents the peak electrical AC power output from the inverter (unit: kW)
    #[validate(minimum = 0.)]
    pub(crate) inverter_peak_power_ac: f64,

    /// Peak power; represents the peak electrical DC power input to the inverter (unit: kW)
    #[validate(minimum = 0.)]
    pub(crate) inverter_peak_power_dc: f64,

    /// Whether the inverter is considered inside the building
    pub(crate) inverter_is_inside: bool,

    pub(crate) inverter_type: InverterType,
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
    #[serde(
        rename = "Control_VentAdjustMax",
        skip_serializing_if = "Option::is_none"
    )]
    pub(crate) control_vent_adjust_max: Option<String>,

    #[serde(
        rename = "Control_VentAdjustMin",
        skip_serializing_if = "Option::is_none"
    )]
    pub(crate) control_vent_adjust_min: Option<String>,

    #[serde(
        rename = "Control_WindowAdjust",
        skip_serializing_if = "Option::is_none"
    )]
    pub(crate) control_window_adjust: Option<String>,

    /// List of the required inputs for Leaks
    #[serde(rename = "Leaks")]
    pub(crate) leaks: VentilationLeaks,

    /// Provides details about available mechanical ventilation systems
    #[serde(rename = "MechanicalVentilation", default)]
    #[validate]
    pub(crate) mechanical_ventilation: IndexMap<std::string::String, MechanicalVentilation>,

    /// Provides details about available non-mechanical ventilation systems
    #[serde(rename = "Vents")]
    #[validate]
    pub(crate) vents: IndexMap<std::string::String, Vent>,

    /// Maximum ACH (Air Changes per Hour) limit
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = 0.)]
    pub(crate) ach_max_static_calcs: Option<f64>,

    /// Minimum ACH (Air Changes per Hour) limit
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = 0.)]
    pub(crate) ach_min_static_calcs: Option<f64>,

    /// Altitude of dwelling above sea level (unit: m)
    pub(crate) altitude: f64,

    pub(crate) cross_vent_possible: bool,

    /// Indicates the exposure to wind of an air flow path on a facade (can can be open, normal and shielded)
    pub(crate) shield_class: VentilationShieldClass,

    pub(crate) terrain_class: TerrainClass,

    /// Base height of the ventilation zone relative to ground (m)
    #[validate(minimum = 0.)]
    pub(crate) ventilation_zone_base_height: f64,

    /// Initial vent position, 0 = vents closed and 1 = vents fully open
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = 0.)]
    #[validate(maximum = 1.)]
    pub(crate) vent_opening_ratio_init: Option<f64>,
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

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct Vent {
    /// Mid height of air flow path relative to ventilation zone (unit: m)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) mid_height_air_flow_path: f64,

    /// Equivalent area of a vent (unit: cm2)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) area_cm2: f64,

    /// Reference pressure difference for an air terminal device (unit: Pa)
    #[validate(minimum = 0.)]
    pub(crate) pressure_difference_ref: f64,

    #[serde(
        rename = "orientation360",
        deserialize_with = "deserialize_orientation",
        serialize_with = "serialize_orientation"
    )]
    #[validate(minimum = -180.)]
    #[validate(maximum = 180.)]
    pub(crate) orientation: f64,

    /// Tilt angle of the surface from horizontal, between 0 and 180, where 0 means the external surface is facing up, 90 means the external surface is vertical and 180 means the external surface is facing down (unit: ˚
    #[validate(minimum = 0.)]
    #[validate(maximum = 180.)]
    pub(crate) pitch: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct VentilationLeaks {
    /// Height of ventilation zone (unit: m)
    #[validate(minimum = 0.)]
    pub(crate) ventilation_zone_height: f64,

    /// Reference pressure difference (unit: Pa)
    #[validate(minimum = 0.)]
    pub(crate) test_pressure: f64,

    /// Flow rate through (unit: m³/h.m²)
    #[validate(minimum = 0.)]
    pub(crate) test_result: f64,

    /// Reference area of the envelope airtightness index
    #[validate(minimum = 0.)]
    pub(crate) env_area: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[validate(custom = validate_supply_air_temp_ctrl)]
pub struct MechanicalVentilation {
    /// Supply air flow rate control
    #[serde(rename = "sup_air_flw_ctrl")]
    pub(crate) supply_air_flow_rate_control: SupplyAirFlowRateControlType,

    /// Supply air temperature control
    #[serde(rename = "sup_air_temp_ctrl")]
    pub(crate) supply_air_temperature_control_type: SupplyAirTemperatureControlType,

    /// MVHR efficiency
    #[serde(rename = "mvhr_eff", skip_serializing_if = "Option::is_none")]
    #[validate(minimum = 0.)]
    pub(crate) mvhr_efficiency: Option<f64>,

    /// Location of the MVHR unit (inside or outside the thermal envelope)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) mvhr_location: Option<MVHRLocation>,

    #[serde(rename = "Control", skip_serializing_if = "Option::is_none")]
    pub(crate) control: Option<String>,

    /// Specific fan power, inclusive of any in use factors (unit: W/l/s)
    #[serde(rename = "SFP")]
    #[validate(minimum = 0.)]
    pub(crate) sfp: f64,

    /// Adjustment factor to be applied to SFP to account for e.g. type of ducting. Typical range 1 - 2.5
    #[serde(default = "default_sfp_in_use_factor", rename = "SFP_in_use_factor")]
    #[validate(minimum = 1.)]
    pub(crate) sfp_in_use_factor: f64,

    #[serde(rename = "EnergySupply")]
    pub(crate) energy_supply: String,

    /// (unit: m³/hour)
    #[validate(minimum = 0.)]
    pub(crate) design_outdoor_air_flow_rate: f64,

    /// List of ductworks installed in this ventilation system
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub(crate) ductwork: Vec<MechanicalVentilationDuctwork>,

    #[serde(flatten)]
    #[validate]
    pub(crate) vent_data: MechVentData,
}

const fn default_sfp_in_use_factor() -> f64 {
    1.
}

/// Validate that only implemented supply air temperature control types are used.
fn validate_supply_air_temp_ctrl(
    data: &MechanicalVentilation,
) -> Result<(), serde_valid::validation::Error> {
    match data.supply_air_temperature_control_type {
        SupplyAirTemperatureControlType::NoControl => Ok(()),
        _ => custom_validation_error(
            format!(
                "Supply air temperature control type {} is not currently implemented. Only NO_CTRL is supported. Other values would be silently overwritten by the ventilation engine.",
                serde_json::to_value(data.supply_air_temperature_control_type).unwrap().as_str().unwrap()
            )
        ),
    }
}

/// Enum to encapsulate the data provided with the different types of vent.
// NB. This type replaces the upstream MechVentType, which is reflected in the variants here.
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "vent_type")]
pub enum MechVentData {
    #[serde(rename = "MVHR")]
    Mvhr {
        position_intake: MechanicalVentilationPosition,

        position_exhaust: MechanicalVentilationPosition,
    },
    #[serde(rename = "Intermittent MEV")]
    IntermittentMev {
        #[serde(flatten, deserialize_with = "deserialize_maybe_nested_position")]
        position_exhaust: MechanicalVentilationPosition,
    },
    #[serde(rename = "Centralised continuous MEV")]
    CentralisedContinuousMev {
        #[serde(flatten, deserialize_with = "deserialize_maybe_nested_position")]
        position_exhaust: MechanicalVentilationPosition,
    },
    #[serde(rename = "Decentralised continuous MEV")]
    DecentralisedContinuousMev {
        #[serde(flatten, deserialize_with = "deserialize_maybe_nested_position")]
        position_exhaust: MechanicalVentilationPosition,
    },
    #[serde(rename = "Positive input ventilation")]
    PositiveInputVentilation {
        #[serde(flatten, deserialize_with = "deserialize_maybe_nested_position")]
        position_exhaust: MechanicalVentilationPosition,
    },
}

// this macro allows for position_exhaust fields to either be nested under that field in the JSON,
// or to be found at the root level, and to deserialise to an identical representation either way
serde_with::flattened_maybe!(deserialize_maybe_nested_position, "position_exhaust");

impl MechVentData {
    /// Check if this ventilation type is a balanced system (both supply and extract).
    fn is_balanced(&self) -> bool {
        matches!(self, Self::Mvhr { .. })
    }

    /// Check if this ventilation type is an extract-only system.
    fn is_extract_only(&self) -> bool {
        matches!(
            self,
            Self::IntermittentMev { .. }
                | Self::CentralisedContinuousMev { .. }
                | Self::DecentralisedContinuousMev { .. }
        )
    }

    /// Check if this ventilation type is a supply-only system.
    fn is_supply_only(&self) -> bool {
        matches!(self, Self::PositiveInputVentilation { .. })
    }

    /// Check if this ventilation type includes supply air.
    fn has_supply(&self) -> bool {
        self.is_balanced() || self.is_supply_only()
    }

    /// Check if this ventilation type includes extract air.
    fn has_extract(&self) -> bool {
        self.is_balanced() || self.is_extract_only()
    }

    /// Check if this ventilation type operates continuously.
    fn is_continuous(&self) -> bool {
        matches!(
            self,
            Self::CentralisedContinuousMev { .. }
                | Self::DecentralisedContinuousMev { .. }
                | Self::Mvhr { .. }
        )
    }

    /// Check if this ventilation type operates intermittently.
    fn is_intermittent(&self) -> bool {
        matches!(self, Self::IntermittentMev { .. })
    }

    fn vent_type(&self) -> &'static str {
        match self {
            Self::Mvhr { .. } => "MVHR",
            Self::IntermittentMev { .. } => "Intermittent MEV",
            Self::CentralisedContinuousMev { .. } => "Centralised continuous MEV",
            Self::DecentralisedContinuousMev { .. } => "Decentralised continuous MEV",
            Self::PositiveInputVentilation { .. } => "Positive input ventilation",
        }
    }
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub struct MechanicalVentilationPosition {
    #[serde(
        rename = "orientation360",
        deserialize_with = "deserialize_orientation",
        serialize_with = "serialize_orientation"
    )]
    #[validate(minimum = -180.)]
    #[validate(maximum = 180.)]
    orientation: f64,

    /// Tilt angle of the surface from horizontal, between 0 and 180, where 0 means the external surface is facing up, 90 means the external surface is vertical and 180 means the external surface is facing down (unit: ˚
    #[validate(minimum = 0.)]
    #[validate(maximum = 180.)]
    pub(crate) pitch: f64,

    /// Mid height of air flow path relative to ventilation zone (unit: m)
    #[validate(exclusive_minimum = 0.)]
    mid_height_air_flow_path: f64,
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

    #[serde(rename = "LOAD")]
    Load,
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

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct MechanicalVentilationDuctwork {
    /// Whether the cross-section of duct is circular or rectangular (square)
    pub(crate) cross_section_shape: DuctShape,

    #[serde(skip_serializing_if = "Option::is_none")]
    /// Cross-sectional perimeter length of rectangular ductwork (unit: mm)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) duct_perimeter_mm: Option<f64>,

    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(exclusive_minimum = 0.)]
    pub(crate) internal_diameter_mm: Option<f64>,

    /// (unit: mm)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(exclusive_minimum = 0.)]
    pub(crate) external_diameter_mm: Option<f64>,

    /// (unit: m)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) length: f64,

    /// Thermal conductivity of the insulation (unit: W / m K)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) insulation_thermal_conductivity: f64,

    /// (unit: mm)
    #[validate(minimum = 0.)]
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

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct AirTerminalDevice {
    /// Equivalent area of the air terminal device (unit: cm2)
    #[validate(exclusive_minimum = 0.)]
    area_cm2: f64,

    /// Reference pressure difference for an air terminal device (unit: Pa)
    #[validate(minimum = 0.)]
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

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
// #[serde(deny_unknown_fields)] // TODO: possibly restore `deny_unknown_fields` declaration after 0.36
pub(crate) struct ApplianceLoadShifting {
    #[serde(rename = "Control", skip_serializing_if = "Option::is_none")]
    pub(crate) control: Option<String>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) priority: Option<isize>,

    /// Maximum time that an event may be shifted away from when it was originally intended to occur. This may be up to 24 hours. (unit: hours)
    #[validate(minimum = 0.)]
    #[validate(maximum = 24.)]
    pub(crate) max_shift_hrs: f64,

    /// Value above which the sum demand multiplied by the
    /// weight should not exceed. If a 7hr tariff were used for
    /// the weight timeseries, then this would be a cost. If this value
    /// is 0, then all events will be shifted to the optimal time
    /// within the window, otherwise they will be shifted to the earliest
    /// time at which the weighted demand goes below the limit.
    /// (if there is no time in the window when the weighted demand is below
    /// the limit, then the optimum is chosen.)
    #[validate(minimum = 0.)]
    pub(crate) demand_limit_weighted: f64,

    /// This may be, for example, the hourly cost per
    /// kWh of a 7hr tariff, but could be any time series.
    /// The sum of the other demand and the demand of the
    /// appliance in question at any given time is multiplied
    /// by the value of this timeseries to obtain a figure
    /// that is used to determine whether to shift an event,
    /// and when would be the most appropriate time to shift
    /// the event to.
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

static FHS_SCHEMA_VALIDATOR: LazyLock<Validator> = LazyLock::new(|| {
    let schema = serde_json::from_str(include_str!("../schemas/input_fhs.schema.json")).unwrap();
    jsonschema::validator_for(&schema).unwrap()
});

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

#[derive(Debug, Error)]
#[error("Error accessing JSON during FHS preprocessing: {0}")]
pub(crate) struct JsonAccessError(String);

pub(crate) fn json_error<T: Into<String>>(message: T) -> JsonAccessError {
    JsonAccessError(message.into())
}

pub(crate) type JsonAccessResult<T> = Result<T, JsonAccessError>;

#[derive(Clone, Debug)]
pub struct InputForProcessing {
    pub(crate) input: JsonValue,
}

/// This type makes methods available for restricted access by wrappers,
/// in order to work towards a reasonable API for wrappers to interact with inputs rather than
/// the more brittle approach of allowing full access to the input data structure.
/// If the full access is encapsulated within methods here, it becomes possible to update the
/// underlying structure without breaking wrappers.
impl InputForProcessing {
    pub fn init_with_json(
        json: impl Read,
        schema_reference: &SchemaReference,
    ) -> Result<Self, anyhow::Error> {
        let input_for_processing = Self::init_with_json_skip_validation(json)?;

        let validator = match schema_reference {
            SchemaReference::Core => &CORE_SCHEMA_VALIDATOR,
            #[cfg(feature = "fhs")]
            SchemaReference::Fhs => &FHS_SCHEMA_VALIDATOR,
        };

        if let BasicOutput::Invalid(errors) = validator.apply(&input_for_processing.input).basic() {
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

    pub(crate) fn as_input(&self) -> anyhow::Result<Input> {
        serde_json::from_value(self.input.to_owned()).map_err(|err| anyhow!(err))
    }

    pub(crate) fn as_input_for_calc_htc_hlp(&self) -> anyhow::Result<ReducedInputForCalcHtcHlp> {
        serde_json::from_value(self.input.to_owned()).map_err(|err| anyhow!(err))
    }

    pub(crate) fn finalize(self) -> anyhow::Result<Input> {
        // NB. this _might_ in time be a good point to perform a validation against the core schema - or it might not
        // if let BasicOutput::Invalid(errors) =
        //     CORE_INCLUDING_FHS_VALIDATOR.apply(&self.input).basic()
        // {
        //     bail!(
        //         "Wrapper formed invalid JSON for the core schema: {}",
        //         serde_json::to_value(errors)?.to_json_string_pretty()?
        //     );
        // }
        serde_json::from_value(self.input).map_err(|err| anyhow!(err))
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
        self.root_mut()?.shift_remove(root_key);

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
        self.root_mut()?.insert("InternalGains".into(), json!({}));

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

    #[cfg(test)]
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

    pub(crate) fn water_distribution(&self) -> anyhow::Result<WaterDistribution> {
        Ok(self
            .root_object("HotWaterDemand")?
            .get("Distribution")
            .map(|node| serde_json::from_value::<WaterDistribution>(node.to_owned()))
            .transpose()?
            .unwrap_or_default())
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
                events
                    .as_object()
                    .map(|events| events.values().filter_map(JsonValue::as_array))
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

    #[cfg(test)]
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
                energy_supply.shift_remove("diverter");
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
            energy_supply.shift_remove("ElectricBattery");
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
        // we use .shift_remove instead of remove here
        // to preserve the relative order of the appliances
        self.root_object_entry_mut("Appliances")?
            .shift_remove(appliance_key);

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
            .shift_remove("MechanicalVentilation");
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

    pub(crate) fn heat_source_wet(
        &self,
    ) -> anyhow::Result<IndexMap<std::string::String, HeatSourceWetDetails>> {
        self.root()?
            .get("HeatSourceWet")
            .and_then(|value| value.as_object())
            .into_iter()
            .flatten()
            .map(|(name, source)| Ok((name.to_owned(), serde_json::from_value(source.clone())?)))
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
    use crate::core::space_heat_demand::ventilation::MechVentType;
    use rstest::*;
    use serde::de::DeserializeOwned;
    use std::fs::File;
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
    fn should_successfully_parse_all_core_demo_files(core_files: Vec<DirEntry>) {
        for entry in core_files {
            let parsed =
                ingest_for_processing(File::open(entry.path()).unwrap(), &SchemaReference::Core);
            assert!(
                parsed.is_ok(),
                "error was {:?} when parsing file {}",
                parsed.err().unwrap(),
                entry.file_name().to_str().unwrap()
            );
        }
    }

    #[rstest]
    fn should_successfully_deserialise_all_core_demo_files(core_files: Vec<DirEntry>) {
        for entry in core_files {
            let input: Result<Input, _> =
                serde_json::from_reader(File::open(entry.path()).unwrap());
            assert!(
                input.is_ok(),
                "error was {:?} when parsing file {}",
                input.err().unwrap(),
                entry.file_name().to_str().unwrap()
            );
        }
    }

    #[rstest]
    fn test_all_demo_files_deserialize_and_serialize(core_files: Vec<DirEntry>) {
        for entry in core_files {
            let input: Input =
                serde_json::from_reader(BufReader::new(File::open(entry.path()).unwrap()))
                    .unwrap_or_else(|e| {
                        panic!(
                            "Failed deserializing {} - {e:?}",
                            entry.file_name().to_str().unwrap()
                        )
                    });
            let json = serde_json::to_string_pretty(&input.clone()).unwrap();
            let recreated_input: Input = serde_json::from_str(&json).unwrap();
            assert_eq!(input, recreated_input,);
        }
    }

    #[fixture]
    fn baseline_demo_file_json() -> JsonValue {
        let file = File::open("./examples/input/core/demo.json").unwrap();
        serde_json::from_reader(file).unwrap()
    }

    fn merge_json_onto_base(base: JsonValue, add: JsonValue) -> JsonValue {
        let mut merged = base;
        for (key, value) in add.as_object().unwrap() {
            merged
                .as_object_mut()
                .unwrap()
                .insert(key.to_string(), value.clone());
        }
        merged
    }

    /// Test WWHRS validation with empty configuration.
    #[rstest]
    fn test_validate_shower_waste_water_heat_recovery_systems_empty(
        baseline_demo_file_json: JsonValue,
    ) {
        let modified_input = merge_json_onto_base(
            baseline_demo_file_json,
            json!({
                "WWHRS": null,
                "HotWaterDemand": {}
            }),
        );
        let input: Result<Input, _> = serde_json::from_value(modified_input);
        assert!(input.is_ok());
    }

    /// Test WWHRS validation with valid mixer shower configuration.
    #[rstest]
    fn test_validate_shower_waste_water_heat_recovery_systems_valid_mixer_shower(
        baseline_demo_file_json: JsonValue,
    ) {
        let modified_input = merge_json_onto_base(
            baseline_demo_file_json,
            json!({
                "WWHRS": {
                    "WWHRS1": {
                        "ColdWaterSource": "mains water",
                        "flow_rates": [5, 7, 9, 11, 13],
                        "system_a_efficiencies": [44.8, 39.1, 34.8, 31.4, 28.6],
                        "system_a_utilisation_factor": 0.7,
                        "system_b_efficiency_factor": 0.81,
                        "type": "WWHRS_Instantaneous",
                    }
                },
                "HotWaterDemand": {
                    "Shower": {
                        "Shower1": {
                            "type": "MixerShower",
                            "ColdWaterSource": "mains water",
                            "flowrate": 15,
                            "WWHRS": "WWHRS1",
                            "WWHRS_configuration": "B",
                        },
                    },
                },
            }),
        );
        let input: Result<Input, _> = serde_json::from_value(modified_input);
        assert!(input.is_ok());
    }

    // there is a test in the Python here that checks for invalid references - not bringing this over because
    // we will detect bad references while building the corpus

    static MINIMAL_TEST_DATA: LazyLock<JsonValue> = LazyLock::new(|| {
        json!({
            "air_flow_rate": 100.0,
            "test_letter": "A",
            "capacity": 8.4,
            "cop": 4.6,
            "design_flow_temp": 35,
            "temp_outlet": 34,
            "temp_source": 0,
            "temp_test": -7,
        })
    });

    fn create_heat_pump_config(source_type: &HeatPumpSourceType) -> JsonValue {
        json!({
            "hp": {
                "type": "HeatPump",
                "EnergySupply": "mains elec",
                "source_type": source_type,
                "sink_type": "Water",
                "backup_ctrl_type": "TopUp",
                "time_delay_backup": 1.0,
                "modulating_control": true,
                "min_modulation_rate_35": 0.35,
                "min_modulation_rate_55": 0.4,
                "time_constant_onoff_operation": 140,
                "temp_return_feed_max": 70.0,
                "temp_lower_operating_limit": -5.0,
                "min_temp_diff_flow_return_for_hp_to_operate": 0.0,
                "var_flow_temp_ctrl_during_test": true,
                "power_heating_circ_pump": 0.015,
                "power_source_circ_pump": 0.010,
                "power_standby": 0.015,
                "power_crankcase_heater": 0.01,
                "power_off": 0.015,
                "power_max_backup": 3.0,
                "MechanicalVentilation": "mechvent1",
                "test_data_EN14825": [*MINIMAL_TEST_DATA],
            }
        })
    }

    fn create_ventilation_config(vent_type: &MechVentType) -> JsonValue {
        if *vent_type == MechVentType::Mvhr {
            json!({
                "sup_air_flw_ctrl": "ODA",
                "sup_air_temp_ctrl": "NO_CTRL",
                "vent_type": vent_type,
                "SFP": 1.5,
                "EnergySupply": "mains elec",
                "design_outdoor_air_flow_rate": 0.5,
                "position_intake": {
                    "orientation360": 180,
                    "pitch": 90,
                    "mid_height_air_flow_path": 3.0,
                },
                "position_exhaust": {
                    "orientation360": 0,
                    "pitch": 90,
                    "mid_height_air_flow_path": 2.0,
                },
                "mvhr_eff": 0.80,
                "mvhr_location": "outside",
                "ductwork": [
                    {
                        "cross_section_shape": "circular",
                        "internal_diameter_mm": 200,
                        "external_diameter_mm": 300,
                        "length": 10.0,
                        "insulation_thermal_conductivity": 0.023,
                        "insulation_thickness_mm": 100,
                        "reflective": false,
                        "duct_type": "supply",
                    },
                    {
                        "cross_section_shape": "circular",
                        "internal_diameter_mm": 200,
                        "external_diameter_mm": 300,
                        "length": 10.0,
                        "insulation_thermal_conductivity": 0.023,
                        "insulation_thickness_mm": 100,
                        "reflective": false,
                        "duct_type": "extract",
                    },
                    {
                        "cross_section_shape": "circular",
                        "internal_diameter_mm": 200,
                        "external_diameter_mm": 300,
                        "length": 10.0,
                        "insulation_thermal_conductivity": 0.023,
                        "insulation_thickness_mm": 100,
                        "reflective": false,
                        "duct_type": "intake",
                    },
                    {
                        "cross_section_shape": "circular",
                        "internal_diameter_mm": 200,
                        "external_diameter_mm": 300,
                        "length": 10.0,
                        "insulation_thermal_conductivity": 0.023,
                        "insulation_thickness_mm": 100,
                        "reflective": false,
                        "duct_type": "exhaust",
                    },
                ],
            })
        } else {
            json!({
                "sup_air_flw_ctrl": "ODA",
                "sup_air_temp_ctrl": "NO_CTRL",
                "vent_type": vent_type,
                "Control": format!("{}_Control", vent_type.to_string().replace(" ", "_")),
                "SFP": 1.5,
                "EnergySupply": "mains elec",
                "design_outdoor_air_flow_rate": 0.5,
                "orientation360": 180,
                "pitch": 90,
                "mid_height_air_flow_path": 2,
            })
        }
    }

    fn create_control_config() -> JsonValue {
        json!({
            "type": "SetpointTimeControl",
            "start_day": 0,
            "time_series_step": 0.5,
            "schedule": {"main": [{"repeat": 7, "value": "day"}]},
        })
    }

    /// Helper method to create complete exhaust air heat pump configuration.
    fn create_exhaust_air_heat_pump_config(
        base_input: JsonValue,
        source_type: HeatPumpSourceType,
        vent_type: MechVentType,
    ) -> JsonValue {
        let mut modified_input = base_input.clone();
        let modified_input_obj = modified_input.as_object_mut().unwrap();

        modified_input_obj.insert(
            "HeatSourceWet".into(),
            create_heat_pump_config(&source_type),
        );
        let mut cloned_infiltration_ventilation = base_input["InfiltrationVentilation"].clone();
        let new_infiltration_ventilation = cloned_infiltration_ventilation.as_object_mut().unwrap();
        new_infiltration_ventilation.insert(
            "MechanicalVentilation".into(),
            json!({"mechvent1": create_ventilation_config(&vent_type)}),
        );
        modified_input_obj.insert(
            "InfiltrationVentilation".into(),
            json!(new_infiltration_ventilation),
        );

        modified_input
    }

    /// Test that compatible exhaust air heat pump and ventilation combinations pass validation.
    #[rstest]
    fn test_validate_exhaust_air_heat_pump_ventilation_compatibility_valid_combinations(
        baseline_demo_file_json: JsonValue,
    ) {
        let valid_combinations = [
            (HeatPumpSourceType::ExhaustAirMEV, MechVentType::Mvhr),
            (
                HeatPumpSourceType::ExhaustAirMEV,
                MechVentType::CentralisedContinuousMev,
            ),
            (HeatPumpSourceType::ExhaustAirMVHR, MechVentType::Mvhr),
            (HeatPumpSourceType::ExhaustAirMixed, MechVentType::Mvhr),
        ];

        for (source_type, vent_type) in valid_combinations {
            let mut modified_input = create_exhaust_air_heat_pump_config(
                baseline_demo_file_json.clone(),
                source_type,
                vent_type,
            );
            let modified_input_control_node = modified_input
                .get_mut("Control")
                .unwrap()
                .as_object_mut()
                .unwrap();
            modified_input_control_node.insert("MVHR_Control".into(), create_control_config());
            modified_input_control_node.insert(
                "Centralised_continuous_MEV_Control".into(),
                create_control_config(),
            );
            let input: Result<Input, _> = serde_json::from_value(modified_input);
            assert!(input.is_ok());
            let input = input.unwrap();
            assert!(
                input.validate().is_ok(),
                "compatibility expected: source_type: {:?}, vent_type: {:?}",
                source_type,
                vent_type
            );
        }
    }

    #[rstest]
    fn test_validate_exhaust_air_heat_pump_ventilation_compatibility_invalid_combinations(
        baseline_demo_file_json: JsonValue,
    ) {
        let invalid_combinations = [
            (
                HeatPumpSourceType::ExhaustAirMEV,
                MechVentType::IntermittentMev,
            ),
            (
                HeatPumpSourceType::ExhaustAirMEV,
                MechVentType::DecentralisedContinuousMev,
            ),
            (
                HeatPumpSourceType::ExhaustAirMVHR,
                MechVentType::IntermittentMev,
            ),
            (
                HeatPumpSourceType::ExhaustAirMixed,
                MechVentType::DecentralisedContinuousMev,
            ),
        ];

        for (source_type, vent_type) in invalid_combinations {
            let mut modified_input = create_exhaust_air_heat_pump_config(
                baseline_demo_file_json.clone(),
                source_type,
                vent_type,
            );
            let modified_input_control_node = modified_input
                .get_mut("Control")
                .unwrap()
                .as_object_mut()
                .unwrap();
            modified_input_control_node
                .insert("Intermittent_MEV_Control".into(), create_control_config());
            modified_input_control_node.insert(
                "Decentralised_continuous_MEV_Control".into(),
                create_control_config(),
            );
            let input: Input = serde_json::from_value(modified_input).unwrap();
            assert!(input.validate().is_err());
        }
    }

    /// Test edge cases where validation should pass regardless of configuration.
    #[rstest]
    fn test_validate_exhaust_air_heat_pump_ventilation_compatibility_edge_cases(
        baseline_demo_file_json: JsonValue,
    ) {
        // Test case 1: No heat pumps
        let mut modified_input = baseline_demo_file_json.clone();
        let modified_input_infiltration_ventilation = modified_input
            .get_mut("InfiltrationVentilation")
            .unwrap()
            .as_object_mut()
            .unwrap();
        modified_input_infiltration_ventilation.insert(
            "MechanicalVentilation".into(),
            json!({
                "mechvent1": create_ventilation_config(&MechVentType::IntermittentMev)
            }),
        );
        let modified_input_control_node = modified_input
            .get_mut("Control")
            .unwrap()
            .as_object_mut()
            .unwrap();
        modified_input_control_node
            .insert("Intermittent_MEV_Control".into(), create_control_config());
        // should be ok
        let input: Result<Input, _> = serde_json::from_value(modified_input);
        let input = input.unwrap();
        assert!(input.validate().is_ok());

        // Test case 2: No mechanical ventilation
        let mut modified_input = baseline_demo_file_json.clone();
        let modified_input_obj = modified_input.as_object_mut().unwrap();
        modified_input_obj.insert(
            "HeatSourceWet".into(),
            create_heat_pump_config(&HeatPumpSourceType::ExhaustAirMEV),
        );
        let infiltration_ventilation_node = modified_input
            .get_mut("InfiltrationVentilation")
            .unwrap()
            .as_object_mut()
            .unwrap();
        infiltration_ventilation_node.insert(
            "MechanicalVentilation".into(),
            json!({
                "mechvent1": create_ventilation_config(&MechVentType::CentralisedContinuousMev)
            }),
        );
        let modified_input_control_node = modified_input
            .get_mut("Control")
            .unwrap()
            .as_object_mut()
            .unwrap();
        modified_input_control_node.insert(
            "Centralised_continuous_MEV_Control".into(),
            create_control_config(),
        );
        // should be ok
        let input: Result<Input, _> = serde_json::from_value(modified_input);
        let input = input.unwrap();
        assert!(input.validate().is_ok());
    }

    fn assert_range_constraints<T>(valid_example: JsonValue, inputs: JsonValue)
    where
        T: DeserializeOwned + Validate,
    {
        let input_under_test = merge_json_onto_base(valid_example, inputs);
        let input: Result<T, _> = serde_json::from_value(input_under_test);
        assert!(match input {
            Ok(input) => input.validate().is_err(),
            Err(_) => true,
        });
    }

    mod time_series {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(ApplianceGainsDetails {
                events: None,
                standby: None,
                energy_supply: Default::default(),
                gains_fraction: 0.0,
                load_shifting: None,
                priority: None,
                schedule: None,
                start_day: 12,
                time_series_step: 12.0,
            })
            .unwrap()
        }

        #[rstest(inputs,
            // don't need to check start_day is greater than zero as this is enforced by u32 type
            case::at_most_365(json!({"start_day": 366})),
            case::at_least_zero_time_series_step(json!({"time_series_step": -1.0})),
            case::at_most_24_time_series_step(json!({"time_series_step": 25.0})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<ApplianceGainsDetails>(valid_example, inputs);
        }
    }

    mod metadata {
        use super::*;

        fn test_validate_hem_core_version() {
            let metadata = InputMetadata {
                hem_core_version: "0.0".into(),
            };
            assert!(metadata.validate().is_err());
        }
    }

    mod appliance_gains_event {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(ApplianceGainsEvent {
                duration: 7.,
                start: 6.,
                demand_w: 2000.,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::demand_at_least_zero(json!({"demand_W": -1})),
            case::duration_at_least_zero(json!({"duration": -1})),
            case::start_at_least_zero(json!({"start": -1})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<ApplianceGainsEvent>(valid_example, inputs);
        }
    }

    mod appliance_load_shifting {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(ApplianceLoadShifting {
                control: None,
                priority: None,
                max_shift_hrs: 12.,
                demand_limit_weighted: 14.,
                weight_timeseries: vec![],
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::demand_limit_at_least_zero(json!({"demand_limit_weighted": -1})),
            case::max_shift_hrs_at_least_zero(json!({"max_shift_hrs": -1})),
            case::max_shift_hrs_at_most_24(json!({"max_shift_hrs": 25})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<ApplianceLoadShifting>(valid_example, inputs);
        }
    }

    mod cold_water_source {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(ColdWaterSourceDetails {
                start_day: 1,
                temperatures: vec![14.],
                time_series_step: 1.,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::temperatures_at_least_zero(json!({"temperatures": [-1]})),
            case::temperatures_at_most_100(json!({"temperatures": [101]})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<ColdWaterSourceDetails>(valid_example, inputs);
        }
    }

    mod electric_battery {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(ElectricBattery {
                capacity: 1.,
                charge_discharge_efficiency_round_trip: 0.8,
                battery_age: 1.,
                minimum_charge_rate_one_way_trip: 0.001,
                maximum_charge_rate_one_way_trip: 1.5,
                maximum_discharge_rate_one_way_trip: 1.25,
                battery_location: BatteryLocation::Outside,
                grid_charging_possible: true,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::age_at_least_zero(json!({"battery_age": -1})),
            case::capacity_at_least_zero(json!({"capacity": -1})),
            case::maximum_charge_rate_at_least_zero(json!({"maximum_charge_rate_one_way_trip": -1})
            ),
            case::maximum_discharge_rate_at_least_zero(json!({"maximum_discharge_rate_one_way_trip": -1})
            ),
            case::minimum_charge_rate_than_or_equal_to_zero(json!({"minimum_charge_rate_one_way_trip": -1})
            ),
            case::charge_discharge_efficiency_at_least_zero(json!({"charge_discharge_efficiency_round_trip": -1})
            ),
            case::charge_discharge_efficiency_at_most_one(json!({"charge_discharge_efficiency_round_trip": 2})
            ),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<ElectricBattery>(valid_example, inputs);
        }
    }

    mod external_sensor_correlation {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(ExternalSensorCorrelation {
                temperature: 15.,
                max_charge: 0.4,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::temperature_at_least_absolute_zero(json!({"temperature": -274})),
            case::max_charge_at_least_zero(json!({"max_charge": -1})),
            case::max_charge_at_most_one(json!({"max_charge": 2})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<ExternalSensorCorrelation>(valid_example, inputs);
        }
    }

    mod fan_speed_data {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(FanSpeedData {
                temperature_diff: 10.,
                power_output: vec![],
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::power_output_at_least_zero(json!({"power_output": [-1]})),
            case::temperature_diff_at_least_zero(json!({"temperature_diff": -1})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<FanSpeedData>(valid_example, inputs);
        }
    }

    mod boiler {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(HeatPumpBoiler {
                energy_supply: "mains elec".into(),
                energy_supply_aux: "mains elec".into(),
                rated_power: 10.,
                efficiency_full_load: 1.,
                efficiency_part_load: 1.,
                boiler_location: HeatSourceLocation::Internal,
                modulation_load: 1.,
                electricity_circ_pump: 10.,
                electricity_part_load: 1.,
                electricity_full_load: 1.,
                electricity_standby: 10.,
                cost_schedule_hybrid: None,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::efficiency_full_load_at_least_zero(json!({"efficiency_full_load": -1})
            ),
            case::efficiency_full_load_at_most_one(json!({"efficiency_full_load": 2})
            ),
            case::efficiency_part_load_at_least_zero(json!({"efficiency_part_load": -1})
            ),
            case::efficiency_part_load_at_most_one(json!({"efficiency_part_load": 2})
            ),
            case::modulation_load_at_least_zero(json!({"modulation_load": -1})),
            case::modulation_load_at_most_one(json!({"modulation_load": 2})),
            case::electricity_circ_pump_at_least_zero(json!({"electricity_circ_pump": -1})
            ),
            case::electricity_full_load_at_least_zero(json!({"electricity_full_load": -1})
            ),
            case::electricity_part_load_at_least_zero(json!({"electricity_part_load": -1})
            ),
            case::electricity_standby_at_least_zero(json!({"electricity_standby": -1})
            ),
            case::rated_power_greater_than_zero(json!({"rated_power": 0})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<HeatPumpBoiler>(valid_example, inputs);
        }
    }

    mod heat_source_wet {
        use super::*;

        mod heat_battery {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(HeatBattery::Pcm {
                    energy_supply: "mains elec".into(),
                    temp_init: 25.,
                    electricity_circ_pump: 0.0600,
                    electricity_standby: 0.0244,
                    rated_charge_power: 10.,
                    max_rated_losses: 0.22,
                    number_of_units: 1,
                    simultaneous_charging_and_discharging: true,
                    control_charge: "control".into(),
                    heat_storage_k_j_per_k_above_phase_transition: 381.5,
                    heat_storage_k_j_per_k_below_phase_transition: 305.2,
                    heat_storage_k_j_per_k_during_phase_transition: 12317.,
                    phase_transition_temperature_upper: 59.,
                    phase_transition_temperature_lower: 57.,
                    max_temperature: 25.,
                    velocity_in_hex_tube_at_1_l_per_min_m_per_s: 0.035,
                    capillary_diameter_m: 0.0065,
                    a: 0.4,
                    b: 0.5,
                    heat_exchanger_surface_area_m2: 8.83,
                    flow_rate_l_per_min: 10.,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::capillary_diameter_m_greater_than_zero(json!({"capillary_diameter_m": 0})),
                case::flow_rate_l_per_min_greater_than_zero(json!({"flow_rate_l_per_min": 0})),
                case::heat_exchanger_surface_area_m2_greater_than_zero(json!({"heat_exchanger_surface_area_m2": 0})
                ),
                case::heat_storage_kj_per_k_above_phase_transition_greater_than_zero(json!({"heat_storage_kJ_per_K_above_Phase_transition": 0})
                ),
                case::heat_storage_kj_per_k_below_phase_transition_greater_than_zero(json!({"heat_storage_kJ_per_K_below_Phase_transition": 0})
                ),
                case::heat_storage_kj_per_k_during_phase_transition_greater_than_zero(json!({"heat_storage_kJ_per_K_during_Phase_transition": 0})
                ),
                case::rated_charge_power_greater_than_zero(json!({"rated_charge_power": 0})),
                case::velocity_in_hex_tube_at_1_l_per_min_m_per_s_greater_than_zero(json!({"velocity_in_HEX_tube_at_1_l_per_min_m_per_s": 0})
                ),
                case::electricity_circ_pump_at_least_zero(json!({"electricity_circ_pump": -1})
                ),
                case::electricity_standby_at_least_zero(json!({"electricity_standby": -1})
                ),
                case::max_rated_losses_at_least_zero(json!({"max_rated_losses": -1})
                ),
                // case::number_of_units_at_least_zero(json!({"number_of_units": -1})), // not needed as enforced by u32 type
                case::max_temperature_at_least_absolute_zero(json!({"max_temperature": -9999})
                ),
                case::phase_transition_temperature_upper_at_least_absolute_zero(json!({"phase_transition_temperature_upper": -9999})
                ),
                case::phase_transition_temperature_lower_at_least_absolute_zero(json!({"phase_transition_temperature_lower": -9999})
                ),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<HeatBattery>(valid_example, inputs);
            }
        }

        mod hiu {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(HeatSourceWetDetails::Hiu {
                    energy_supply: "mains elec".into(),
                    power_max: 3.0,
                    hiu_daily_loss: 0.8,
                    building_level_distribution_losses: 62.,
                    power_circ_pump: None,
                    power_aux: None,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::power_max_greater_than_zero(json!({"power_max": 0})),
                case::hiu_daily_loss_at_least_zero(json!({"HIU_daily_loss": -1})),
                case::building_level_distribution_losses_at_least_zero(json!({"building_level_distribution_losses": -1})
                ),
                case::power_circ_pump_at_least_zero(json!({"power_circ_pump": -1})),
                case::power_aux_at_least_zero(json!({"power_aux": -1})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<HeatSourceWetDetails>(valid_example, inputs);
            }
        }

        mod heat_pump {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(HeatSourceWetDetails::HeatPump {
                    buffer_tank: None,
                    energy_supply: "mains elec".into(),
                    energy_supply_heat_network: None,
                    mechanical_ventilation: None,
                    backup_control_type: HeatPumpBackupControlType::TopUp,
                    boiler: None,
                    eahp_mixed_max_temp: Some(10.),
                    eahp_mixed_min_temp: Some(10.),
                    min_modulation_rate_20: None,
                    min_modulation_rate_35: None,
                    min_modulation_rate_55: None,
                    min_temp_diff_flow_return_for_hp_to_operate: 1.,
                    modulating_control: true,
                    power_crankcase_heater: 0.01,
                    power_heating_circ_pump: Some(0.015),
                    power_heating_warm_air_fan: Some(0.015),
                    power_max_backup: Some(3.0),
                    power_off: 0.015,
                    power_source_circ_pump: 0.010,
                    power_standby: 0.015,
                    sink_type: HeatPumpSinkType::Water,
                    source_type: HeatPumpSourceType::ExhaustAirMixed,
                    temp_distribution_heat_network: Some(10.),
                    temp_lower_operating_limit: -5.,
                    temp_return_feed_max: Some(70.),
                    test_data_en14825: vec![],
                    time_constant_onoff_operation: 140.,
                    time_delay_backup: None,
                    var_flow_temp_ctrl_during_test: true,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::eahp_mixed_max_temp_at_least_absolute_zero(json!({"eahp_mixed_max_temp": -274})),
                case::eahp_mixed_min_temp_at_least_absolute_zero(json!({"eahp_mixed_min_temp": -274})),
                case::temp_distribution_heat_network_at_least_absolute_zero(json!({"temp_distribution_heat_network": -274})),
                case::temp_lower_operating_limit_at_least_absolute_zero(json!({"temp_lower_operating_limit": -274})),
                case::temp_return_feed_max_at_least_absolute_zero(json!({"temp_return_feed_max": -274})),
                case::min_modulation_rate_20_at_least_zero(json!({"min_modulation_rate_20": -1})),
                case::min_modulation_rate_20_at_most_one(json!({"min_modulation_rate_20": 2})),
                case::min_modulation_rate_35_at_least_zero(json!({"min_modulation_rate_35": -1})),
                case::min_modulation_rate_35_at_most_one(json!({"min_modulation_rate_35": 2})),
                case::min_modulation_rate_55_at_least_zero(json!({"min_modulation_rate_55": -1})),
                case::min_modulation_rate_55_at_most_one(json!({"min_modulation_rate_55": 2})),
                case::power_crankcase_heater_at_least_zero(json!({"power_crankcase_heater": -1})),
                case::power_heating_circ_pump_at_least_zero(json!({"power_heating_circ_pump": -1})),
                case::power_heating_warm_air_fan_at_least_zero(json!({"power_heating_warm_air_fan": -1})),
                case::power_max_backup_at_least_zero(json!({"power_max_backup": -1})),
                case::power_off_at_least_zero(json!({"power_off": -1})),
                case::power_source_circ_pump_at_least_zero(json!({"power_source_circ_pump": -1})),
                case::power_standby_at_least_zero(json!({"power_standby": -1})),
                case::time_constant_onoff_operation_at_least_zero(json!({"time_constant_onoff_operation": -1})),
                case::time_delay_backup_at_least_zero(json!({"time_delay_backup": -1})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<HeatSourceWetDetails>(valid_example, inputs);
            }
        }
    }

    mod hot_water_source {
        use super::*;

        mod combi_boiler {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(HotWaterSourceDetails::CombiBoiler {
                    cold_water_source: "cold water source".into(),
                    heat_source_wet: "heat source wet".into(),
                    separate_dhw_tests: BoilerHotWaterTest::MOnly,
                    rejected_energy_1: None,
                    storage_loss_factor_2: None,
                    rejected_factor_3: None,
                    setpoint_temp: Some(10.),
                    daily_hw_usage: 120.,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::daily_hw_usage_at_least_zero(json!({"daily_HW_usage": -1})),
                case::rejected_energy_1_at_least_zero(json!({"rejected_energy_1": -1})
                ),
                case::rejected_factor_3_at_least_zero(json!({"rejected_factor_3": -1})
                ),
                case::storage_loss_factor_2_at_least_zero(json!({"storage_loss_factor_2": -1})
                ),
                case::setpoint_temp_at_least_absolute_zero(json!({"setpoint_temp": -9999})
                ),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<HotWaterSourceDetails>(valid_example, inputs);
            }
        }

        mod hui {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(HotWaterSourceDetails::Hiu {
                    cold_water_source: "cold water source".into(),
                    heat_source_wet: "heat source wet".into(),
                    setpoint_temp: None,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::setpoint_temp_at_least_absolute_zero(json!({"setpoint_temp": -9999})
                ),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<HotWaterSourceDetails>(valid_example, inputs);
            }
        }

        mod point_of_use {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(HotWaterSourceDetails::PointOfUse {
                    efficiency: Some(0.7),
                    energy_supply: "mains elec".into(),
                    cold_water_source: "cold water source".into(),
                    setpoint_temp: 25.,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::efficiency_at_least_zero(json!({"efficiency": -1})),
                case::efficiency_at_most_one(json!({"efficiency": 2})),
                case::setpoint_temp_at_least_absolute_zero(json!({"setpoint_temp": -9999})
                ),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<HotWaterSourceDetails>(valid_example, inputs);
            }
        }

        mod storage_tank {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(HotWaterSourceDetails::StorageTank {
                    cold_water_source: "cold water source".into(),
                    heat_source: IndexMap::from([(
                        "hp".into(),
                        HeatSource::ServiceWaterRegular {
                            name: "hp".into(),
                            temp_flow_limit_upper: Some(65.),
                            energy_supply: "mains elec".into(),
                            control_min: Some("min_temp".into()),
                            control_max: Some("setpoint_temp_max".into()),
                            heater_position: 0.1,
                            thermostat_position: Some(0.33),
                        },
                    )]),
                    daily_losses: 10.,
                    heat_exchanger_surface_area: None,
                    init_temp: 10.,
                    primary_pipework: None,
                    volume: 100.,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::init_temp_at_least_absolute_zero(json!({"init_temp": -274})),
                case::daily_losses_at_least_zero(json!({"daily_losses": -1})),
                case::heat_exchanger_surface_area_at_least_zero(json!({"heat_exchanger_surface_area": -1})),
                case::volume_at_least_zero(json!({"volume": -1})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<HotWaterSourceDetails>(valid_example, inputs);
            }
        }
    }

    mod mechanical_ventilation_ductwork {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(MechanicalVentilationDuctwork {
                cross_section_shape: DuctShape::Circular,
                duct_perimeter_mm: None,
                internal_diameter_mm: None,
                external_diameter_mm: None,
                length: 100.,
                insulation_thermal_conductivity: 0.8,
                insulation_thickness_mm: 20.,
                reflective: true,
                duct_type: DuctType::Exhaust,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::duct_perimeter_mm_greater_than_zero(json!({"duct_perimeter_mm": 0})),
            case::external_diameter_mm_greater_than_zero(json!({"external_diameter_mm": 0})),
            case::insulation_thermal_conductivity_greater_than_zero(json!({"insulation_thermal_conductivity": 0})
            ),
            case::internal_diameter_mm_greater_than_zero(json!({"internal_diameter_mm": 0})),
            case::length_greater_than_zero(json!({"length": 0})),
            case::insulation_thickness_mm_greater_than_zero(json!({"insulation_thickness_mm": -1})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<MechanicalVentilationDuctwork>(valid_example, inputs);
        }
    }

    mod other_water_use {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(OtherWaterUse {
                flowrate: 10.,
                cold_water_source: "cold water source".into(),
                hot_water_source: None,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::flowrate_at_least_zero(json!({"flowrate": -1})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<OtherWaterUse>(valid_example, inputs);
        }
    }

    // upstream Python has tests that schedule repeater `repeat` fields are non-negative - unnecessary here as field uses unsigned integer type

    mod shower {
        use super::*;

        mod mixer {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(Shower::MixerShower {
                    flowrate: 10.,
                    cold_water_source: "cold water source".into(),
                    hot_water_source: None,
                    wwhrs_config: Default::default(),
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::flowrate_at_least_zero(json!({"flowrate": -1})),
                case::wwhrs_configuration_not_provided_when_wwhrs_not_provided(json!({"WWHRS": null, "WWHRS_configuration": "A"})
                )
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<Shower>(valid_example, inputs);
            }
        }

        mod instant_electric {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(Shower::InstantElectricShower {
                    rated_power: 5.,
                    cold_water_source: "cold water source".into(),
                    energy_supply: "mains elec".into(),
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::rated_power_at_least_zero(json!({"rated_power": -1})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<Shower>(valid_example, inputs);
            }
        }
    }

    mod smart_appliance_battery {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            let base_map = IndexMap::from([("test".into(), vec![0.1, 0.2, 0.3, 0.4])]);

            serde_json::to_value(SmartApplianceBattery {
                battery_state_of_charge: base_map.clone(),
                energy_into_battery_from_generation: base_map.clone(),
                energy_into_battery_from_grid: base_map.clone(),
                energy_out_of_battery: base_map.clone(),
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::battery_state_of_change_at_least_one_map_item(json!({"battery_state_of_charge": {}})),
            case::battery_state_of_change_at_least_one_list_item(json!({"battery_state_of_charge": {"_unmet_demand": []}})),
            case::battery_state_of_change_each_at_most_one(json!({"battery_state_of_charge": {"_unmet_demand": [0, 1, 2, -1]}})),
            case::battery_state_of_change_each_at_least_zero(json!({"battery_state_of_charge": {"_unmet_demand": [0, 0.1, 0.2, -1]}})),
            case::energy_into_battery_from_generation_at_least_one_map_item(json!({"energy_into_battery_from_generation": {}})),
            case::energy_into_battery_from_generation_at_least_one_list_item(json!({"energy_into_battery_from_generation": {"_unmet_demand": []}})),
            case::energy_into_battery_from_grid_at_least_one_map_item(json!({"energy_into_battery_from_grid": {}})),
            case::energy_into_battery_from_grid_at_least_one_list_item(json!({"energy_into_battery_from_grid": {"_unmet_demand": []}})),
            case::energy_out_of_battery_at_least_one_map_item(json!({"energy_out_of_battery": {}})),
            case::energy_out_of_battery_at_least_one_list_item(json!({"energy_out_of_battery": {"_unmet_demand": []}})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<SmartApplianceBattery>(valid_example, inputs);
        }
    }

    mod smart_appliance_control {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            let base_map = IndexMap::from([("test".into(), vec![0.1, 0.2, 0.3, 0.4])]);

            serde_json::to_value(SmartApplianceControlDetails {
                appliances: vec![],
                battery_24hr: SmartApplianceBattery {
                    battery_state_of_charge: base_map.clone(),
                    energy_into_battery_from_generation: base_map.clone(),
                    energy_into_battery_from_grid: base_map.clone(),
                    energy_out_of_battery: base_map.clone(),
                },
                non_appliance_demand_24hr: IndexMap::from([("test".into(), vec![1., 2., 3., 4.])]),
                power_timeseries: IndexMap::from([("test".into(), vec![1., 2., 3., 4.])]),
                time_series_step: 1.,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::non_appliance_demand_24hr_at_least_one_list_item(json!({"non_appliance_demand_24hr": {"_unmet_demand": []}})),
            case::power_timeseries_at_least_one_list_item(json!({"power_timeseries": {"_unmet_demand": []}})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<SmartApplianceControlDetails>(valid_example, inputs);
        }
    }

    mod space_heat_system_heat_source {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(SpaceHeatSystemHeatSource {
                name: "test".into(),
                temp_flow_limit_upper: None,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::temp_flow_limit_upper_at_least_zero(json!({"temp_flow_limit_upper": -1})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<SpaceHeatSystemHeatSource>(valid_example, inputs);
        }
    }

    mod thermal_bridging_linear {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(ThermalBridgingDetails::Linear {
                linear_thermal_transmittance: 0.9,
                length: 10.,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::length_greater_than_zero(json!({"length": 0})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<ThermalBridgingDetails>(valid_example, inputs);
        }
    }

    mod vent {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(Vent {
                mid_height_air_flow_path: 1.3,
                area_cm2: 120.,
                pressure_difference_ref: 3.4,
                orientation: 180. - 123.,
                pitch: 45.,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::area_cm2_greater_than_zero(json!({"area_cm2": 0})),
            case::mid_height_air_flow_path_greater_than_zero(json!({"mid_height_air_flow_path": 0})),
            case::orientation360_at_least_zero(json!({"orientation360": -1})),
            case::orientation360_at_most_360(json!({"orientation360": 361})),
            case::pitch_at_least_zero(json!({"pitch": -1})),
            case::pitch_at_most_one(json!({"pitch": 181})),
            case::pressure_difference_ref_at_least_zero(json!({"pressure_difference_ref": -1})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<Vent>(valid_example, inputs);
        }
    }

    mod ventilation_leaks {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(VentilationLeaks {
                ventilation_zone_height: 6.,
                test_pressure: 50.,
                test_result: 1.2,
                env_area: 220.,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::env_area_at_least_zero(json!({"env_area": -1})),
            case::test_pressure_at_least_zero(json!({"test_pressure": -1})),
            case::test_result_at_least_zero(json!({"test_result": -1})),
            case::ventilation_zone_height_at_least_zero(json!({"ventilation_zone_height": -1})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<VentilationLeaks>(valid_example, inputs);
        }
    }

    mod water_heating_event {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(WaterHeatingEvent {
                start: 1.,
                duration: None,
                volume: None,
                temperature: 12.,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::duration_at_least_zero(json!({"duration": -1})),
            case::start_at_least_zero(json!({"start": -1})),
            case::volume_at_least_zero(json!({"volume": -1})),
            case::temperature_at_least_absolute_zero(json!({"temperature": -274})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<WaterHeatingEvent>(valid_example, inputs);
        }
    }

    mod water_pipework {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(WaterPipework {
                location: WaterPipeworkLocation::Internal,
                internal_diameter_mm: 25.,
                external_diameter_mm: 27.,
                length: 10.,
                insulation_thermal_conductivity: 1.3,
                insulation_thickness_mm: 34.,
                surface_reflectivity: false,
                pipe_contents: PipeworkContents::Water,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::external_diameter_mm_greater_than_zero(json!({"external_diameter_mm": 0})),
            case::insulation_thermal_conductivity_greater_than_zero(json!({"insulation_thermal_conductivity": 0})),
            case::internal_diameter_mm_greater_than_zero(json!({"internal_diameter_mm": 0})),
            case::insulation_thickness_mm_at_least_zero(json!({"insulation_thickness_mm": -1})),
            case::length_at_least_zero(json!({"length": -1})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<WaterPipework>(valid_example, inputs);
        }
    }

    mod water_pipework_simple {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(WaterPipeworkSimple {
                location: WaterPipeworkLocation::External,
                internal_diameter_mm: 25.,
                length: 10.,
                external_diameter_mm: None,
                insulation_thermal_conductivity: None,
                insulation_thickness_mm: None,
                surface_reflectivity: None,
                pipe_contents: None,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::internal_diameter_mm_greater_than_zero(json!({"internal_diameter_mm": 0})),
            case::length_greater_than_zero(json!({"length": 0})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<WaterPipeworkSimple>(valid_example, inputs);
        }
    }

    mod wet_emitter {
        use super::*;

        mod radiator {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(WetEmitter::Radiator {
                    exponent: 1.2,
                    frac_convective: 0.9,
                    constant_data: RadiatorConstantData::Constant { constant: 1.2 },
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::c_greater_than_zero(json!({"c": 0})),
                case::n_greater_than_zero(json!({"n": 0})),
                case::frac_convective_at_least_zero(json!({"frac_convective": -1})),
                case::frac_convective_at_most_one(json!({"frac_convective": 2})),
                case::length_greater_than_zero(json!({"c": null, "length": 0, "c_per_m": 1})),
                case::c_per_m_greater_than_zero(json!({"c": null, "c_per_m": 0, "length": 3})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<WetEmitter>(valid_example, inputs);
            }
        }

        mod ufh {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(WetEmitter::Ufh {
                    equivalent_specific_thermal_mass: 80.,
                    system_performance_factor: 5.,
                    emitter_floor_area: 40.,
                    frac_convective: 0.9,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::equivalent_specific_thermal_mass_greater_than_zero(json!({"equivalent_specific_thermal_mass": 0})),
                case::system_performance_factor_greater_than_zero(json!({"system_performance_factor": 0})),
                case::emitter_floor_area_at_least_zero(json!({"emitter_floor_area": -1})),
                case::frac_convective_at_least_zero(json!({"frac_convective": -1})),
                case::frac_convective_at_most_one(json!({"frac_convective": 2})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<WetEmitter>(valid_example, inputs);
            }
        }

        mod fancoil {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(WetEmitter::Fancoil {
                    n_units: 1,
                    frac_convective: 0.8,
                    fancoil_test_data: FancoilTestData {
                        fan_speed_data: vec![
                            FanSpeedData {
                                temperature_diff: 80.,
                                power_output: vec![2.7, 3.6, 5., 5.3, 6.2, 7.4],
                            },
                            FanSpeedData {
                                temperature_diff: 70.,
                                power_output: vec![2.3, 3.1, 4.2, 4.5, 5.3, 6.3],
                            },
                            FanSpeedData {
                                temperature_diff: 60.,
                                power_output: vec![1.9, 2.6, 3.5, 3.8, 4.4, 5.3],
                            },
                            FanSpeedData {
                                temperature_diff: 50.,
                                power_output: vec![1.5, 2., 2.8, 3., 3.5, 4.2],
                            },
                        ],
                        fan_power_w: vec![15., 19., 25., 33., 43., 56.],
                    },
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::n_units_greater_than_zero(json!({"n_units": 0})),
                case::frac_convective_at_least_zero(json!({"frac_convective": -1})),
                case::frac_convective_at_most_one(json!({"frac_convective": 2})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<WetEmitter>(valid_example, inputs);
            }
        }
    }

    mod window_shading {
        use super::*;

        mod obstacle {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(WindowShadingObject::Obstacle {
                    height: 10.,
                    distance: 10.,
                    transparency: 0.6,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::height_at_least_zero(json!({"height": -1})),
                case::distance_at_least_zero(json!({"distance": -1})),
                case::transparency_at_least_zero(json!({"transparency": -1})),
                case::transparency_at_most_one(json!({"transparency": 2})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<WindowShadingObject>(valid_example, inputs);
            }
        }

        mod object {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(WindowShadingObject::SideFinLeft {
                    depth: 10.,
                    distance: 1.5,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::depth_at_least_zero(json!({"depth": -1})),
                case::distance_at_least_zero(json!({"distance": -1})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<WindowShadingObject>(valid_example, inputs);
            }
        }
    }

    mod bath {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(BathDetails {
                size: 10.,
                cold_water_source: "cold water source".into(),
                hot_water_source: None,
                flowrate: 2.,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::flowrate_at_least_zero(json!({"flowrate": -1})),
            case::size_greater_than_zero(json!({"size": 0})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<BathDetails>(valid_example, inputs);
        }
    }

    mod energy_supply {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(EnergySupplyDetails {
                fuel: FuelType::EnergyFromEnvironment,
                diverter: None,
                electric_battery: None,
                factor: None,
                priority: None,
                is_export_capable: true,
                threshold_charges: None,
                threshold_prices: None,
                tariff: None,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::threshold_charges_at_least_12_items(json!({"threshold_charges": [0, 1, 1]})),
            case::threshold_charges_at_most_12_items(json!({"threshold_charges": [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]})),
            case::threshold_charges_item_at_most_one(json!({"threshold_charges": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2]})),
            case::threshold_charges_item_at_least_zero(json!({"threshold_charges": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1]})),
            case::threshold_prices_at_least_12_items(json!({"threshold_prices": [0, 1, 1]})),
            case::threshold_prices_at_most_12_items(json!({"threshold_prices": [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<EnergySupplyDetails>(valid_example, inputs);
        }
    }

    mod heat_source {
        use super::*;

        mod immersion_heater {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(HeatSource::ImmersionHeater {
                    power: 12.,
                    energy_supply: "mains elec".into(),
                    control_min: Some("control min".into()),
                    control_max: Some("control max".into()),
                    heater_position: 0.7,
                    thermostat_position: None,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::power_at_least_zero(json!({"power": -1})),
                case::heater_position_at_least_zero(json!({"heater_position": -1})),
                case::heater_position_at_most_one(json!({"heater_position": 2})),
                case::thermostat_position_at_least_zero(json!({"thermostat_position": -1})),
                case::thermostat_position_at_most_one(json!({"thermostat_position": 2})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<HeatSource>(valid_example, inputs);
            }
        }

        mod solar_thermal_system {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(HeatSource::SolarThermalSystem {
                    solar_cell_location: SolarCollectorLoopLocation::Out,
                    area_module: 10.,
                    modules: 1,
                    peak_collector_efficiency: 0.8,
                    incidence_angle_modifier: 0.9,
                    first_order_hlc: 3.5,
                    second_order_hlc: 0.,
                    collector_mass_flow_rate: 1.,
                    power_pump: 100.,
                    power_pump_control: 10.,
                    energy_supply: "mains elec".into(),
                    tilt: 56.,
                    orientation: -54.,
                    solar_loop_piping_hlc: 0.5,
                    heater_position: 0.8,
                    thermostat_position: None,
                    control_max: "control max".into(),
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::area_module_at_least_zero(json!({"area_module": -1})),
                case::collector_mass_flow_rate_at_least_zero(json!({"collector_mass_flow_rate": -1})),
                case::first_order_hlc_at_least_zero(json!({"first_order_hlc": -1})),
                case::incidence_angle_modifier_at_least_zero(json!({"incidence_angle_modifier": -1})),
                case::incidence_angle_modifier_at_least_zero(json!({"incidence_angle_modifier": 2})),
                case::peak_collector_efficiency_at_least_zero(json!({"peak_collector_efficiency": -1})),
                case::peak_collector_efficiency_at_most_one(json!({"peak_collector_efficiency": 2})),
                case::power_pump_at_least_zero(json!({"power_pump": -1})),
                case::power_pump_control_at_least_zero(json!({"power_pump_control": -1})),
                case::second_order_hlc_at_least_zero(json!({"second_order_hlc": -1})),
                case::solar_loop_piping_hlc_at_least_zero(json!({"solar_loop_piping_hlc": -1})),
                case::modules_greater_than_zero(json!({"modules": 0})),
                case::orientation_at_most_180(json!({"orientation360": -1})),
                case::orientation_at_least_minus_180(json!({"orientation360": 361})),
                case::heater_position_at_least_zero(json!({"heater_position": -1})),
                case::heater_position_at_most_one(json!({"heater_position": 2})),
                case::thermostat_position_at_least_zero(json!({"thermostat_position": -1})),
                case::thermostat_position_at_most_one(json!({"thermostat_position": 2})),
                case::tilt_at_least_zero(json!({"tilt": -1})),
                case::tilt_at_most_90(json!({"tilt": 91})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<HeatSource>(valid_example, inputs);
            }
        }

        mod wet_service {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(HeatSource::ServiceWaterRegular {
                    name: Default::default(),
                    temp_flow_limit_upper: None,
                    energy_supply: Default::default(),
                    heater_position: 0.8,
                    thermostat_position: None,
                    control_max: Some("control max".into()),
                    control_min: None,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::heater_position_at_least_zero(json!({"heater_position": -1})),
                case::heater_position_at_most_one(json!({"heater_position": 2})),
                case::thermostat_position_at_least_zero(json!({"thermostat_position": -1})),
                case::thermostat_position_at_most_one(json!({"thermostat_position": 2})),
                case::temp_flow_limit_upper_at_least_absolute_zero(json!({"temp_flow_limit_upper": -274})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<HeatSource>(valid_example, inputs);
            }
        }

        mod heat_pump_hot_water_only {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(HeatSource::HeatPumpHotWaterOnly {
                    power_max: 10.,
                    vol_hw_daily_average: 10.,
                    tank_volume_declared: 10.,
                    heat_exchanger_surface_area_declared: 1.3,
                    daily_losses_declared: 2.3,
                    in_use_factor_mismatch: 0.4,
                    test_data: HeatPumpHotWaterTestData {
                        l: None,
                        m: HeatPumpHotWaterOnlyTestDatum {
                            cop_dhw: 2.5,
                            hw_tapping_prof_daily_total: 5.845,
                            energy_input_measured: 2.338,
                            power_standby: 0.02,
                            hw_vessel_loss_daily: 2.0,
                        },
                    },
                    energy_supply: "mains elec".into(),
                    control_min: "control min".into(),
                    control_max: "control max".into(),
                    heater_position: 1.,
                    thermostat_position: None,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::heater_position_at_least_zero(json!({"heater_position": -1})),
                case::heater_position_at_most_one(json!({"heater_position": 2})),
                case::thermostat_position_at_least_zero(json!({"thermostat_position": -1})),
                case::thermostat_position_at_most_one(json!({"thermostat_position": 2})),
                case::daily_losses_declared_at_least_zero(json!({"daily_losses_declared": -1})),
                case::heat_exchanger_surface_area_declared_at_least_zero(json!({"heat_exchanger_surface_area_declared": -1})),
                case::in_use_factor_mismatch_at_least_zero(json!({"in_use_factor_mismatch": -1})),
                case::power_max_at_least_zero(json!({"power_max": -1})),
                case::vol_hw_daily_average_greater_than_zero(json!({"vol_hw_daily_average": 0})),
                case::tank_volume_declared_greater_than_zero(json!({"tank_volume_declared": 0})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<HeatSource>(valid_example, inputs);
            }
        }
    }

    mod smart_hot_water_tank {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(HotWaterSourceDetails::SmartHotWaterTank {
                volume: 10.,
                power_pump_kw: 5.,
                max_flow_rate_pump_l_per_min: 10.,
                temp_usable: 40.,
                temp_setpnt_max: "test".into(),
                daily_losses: 2.3,
                init_temp: 15.,
                cold_water_source: "cold water source".into(),
                energy_supply_pump: "mains elec".into(),
                heat_source: Default::default(),
                primary_pipework: None,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::daily_losses_at_least_zero(json!({"daily_losses": -1})),
            case::max_flow_rate_pump_l_per_min_at_least_zero(json!({"max_flow_rate_pump_l_per_min": -1})),
            case::power_pump_kw_at_least_zero(json!({"power_pump_kW": -1})),
            case::volume_at_least_zero(json!({"volume": -1})),
            case::init_temp_at_least_absolute_zero(json!({"init_temp": -274})),
            case::temp_usable_at_least_absolute_zero(json!({"temp_usable": -274})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<HotWaterSourceDetails>(valid_example, inputs);
        }
    }

    mod mechanical_ventilation {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(MechanicalVentilation {
                supply_air_flow_rate_control: SupplyAirFlowRateControlType::Oda,
                supply_air_temperature_control_type: SupplyAirTemperatureControlType::NoControl,
                mvhr_efficiency: None,
                mvhr_location: None,
                control: None,
                sfp: 1.2,
                sfp_in_use_factor: 0.0,
                energy_supply: "mains elec".into(),
                design_outdoor_air_flow_rate: 10.,
                ductwork: vec![],
                vent_data: MechVentData::Mvhr {
                    position_intake: MechanicalVentilationPosition {
                        orientation: 0.,
                        pitch: 90.,
                        mid_height_air_flow_path: 3.0,
                    },
                    position_exhaust: MechanicalVentilationPosition {
                        orientation: 180.,
                        pitch: 90.,
                        mid_height_air_flow_path: 2.0,
                    },
                },
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::sfp_at_least_zero(json!({"SFP": -1})),
            case::design_outdoor_air_flow_rate_at_least_zero(json!({"design_outdoor_air_flow_rate": -1})),
            case::mvhr_eff_at_least_zero(json!({"mvhr_eff": -1})),
            case::supply_air_temperature_control_type_must_be_no_control(json!({"sup_air_temp_ctrl": SupplyAirTemperatureControlType::Constant})),
            case::mech_vent_mvhr_needs_both_position_intake_and_position_exhaust(json!({"vent_type": "MVHR", "position_intake": null})),
            case::mech_vent_mvhr_needs_both_position_intake_and_position_exhaust(json!({"vent_type": "MVHR", "position_exhaust": null})),
            case::mech_vent_mvhr_uses_position_intake_and_position_exhaust_not_legacy_fields(json!({"orientation360": 234})),
            case::mech_vent_mvhr_uses_position_intake_and_position_exhaust_not_legacy_fields(json!({"pitch": 12})),
            case::mech_vent_mvhr_uses_position_intake_and_position_exhaust_not_legacy_fields(json!({"mid_height_air_flow_path": 4})),
            case::mech_vent_piv_uses_position_exhaust_xor_legacy_fields(json!({
                    "vent_type": "Positive input ventilation",
                    "position_exhaust": {
                        "orientation360": 0,
                        "pitch": 90,
                        "mid_height_air_flow_path": 2.0,
                    },
                    "orientation360": 345,
                })),
            case::mech_vent_piv_uses_position_exhaust_xor_legacy_fields(json!({
                    "vent_type": "Positive input ventilation",
                    "position_exhaust": {
                        "orientation360": 0,
                        "pitch": 90,
                        "mid_height_air_flow_path": 2.0,
                    },
                    "pitch": 34,
                })),
            case::mech_vent_piv_uses_position_exhaust_xor_legacy_fields(json!({
                    "vent_type": "Positive input ventilation",
                    "position_exhaust": {
                        "orientation360": 0,
                        "pitch": 90,
                        "mid_height_air_flow_path": 2.0,
                    },
                    "mid_height_air_flow_path": 3,
                })),
            case::mech_vent_piv_uses_position_exhaust_or_legacy_fields(json!({
                    "vent_type": "Positive input ventilation",
                    "position_exhaust": null,
                    "orientation360": null,
                    "pitch": null,
                    "mid_height_air_flow_path": 3,
                })),
            case::mech_vent_piv_uses_position_exhaust_or_legacy_fields(json!({
                    "vent_type": "Positive input ventilation",
                    "position_exhaust": null,
                    "orientation360": 234,
                    "pitch": null,
                    "mid_height_air_flow_path": null,
                })),
            case::mech_vent_piv_uses_position_exhaust_or_legacy_fields(json!({
                    "vent_type": "Positive input ventilation",
                    "position_exhaust": null,
                    "orientation360": null,
                    "pitch": 23,
                    "mid_height_air_flow_path": null,
                })),
            case::mech_vent_piv_does_not_have_position_intake(json!({
                    "vent_type": "Positive input ventilation",
                    "position_exhaust": {
                        "orientation360": 0,
                        "pitch": 90,
                        "mid_height_air_flow_path": 2.0,
                    },
                    "position_intake": {
                        "orientation360": 180,
                        "pitch": 90,
                        "mid_height_air_flow_path": 3.0,
                    },
                })),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<MechanicalVentilation>(valid_example, inputs);
        }
    }

    mod photovoltaic_system {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(PhotovoltaicSystem {
                _type: Default::default(),
                peak_power: 10.,
                ventilation_strategy: PhotovoltaicVentilationStrategy::Unventilated,
                pitch: 90.,
                orientation: 150.,
                base_height: 10.,
                height: 10.,
                width: 10.,
                energy_supply: "mains elec".into(),
                shading: vec![],
                inverter_peak_power_ac: 1.4,
                inverter_peak_power_dc: 3.5,
                inverter_is_inside: true,
                inverter_type: InverterType::StringInverter,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::base_height_at_least_zero(json!({"base_height": -1})),
            case::height_at_least_zero(json!({"height": -1})),
            case::inverter_peak_power_ac_at_least_zero(json!({"inverter_peak_power_ac": -1})),
            case::inverter_peak_power_dc_at_least_zero(json!({"inverter_peak_power_dc": -1})),
            case::peak_power_at_least_zero(json!({"peak_power": -1})),
            case::width_greater_than_zero(json!({"width": 0})),
            case::orientation_at_most_180(json!({"orientation360": -1})),
            case::orientation_at_least_minus_180(json!({"orientation360": 361})),
            case::pitch_at_least_zero(json!({"pitch": -1})),
            case::pitch_at_most_90(json!({"pitch": 91})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<PhotovoltaicSystem>(valid_example, inputs);
        }
    }

    mod photovoltaic_system_with_panels {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(PhotovoltaicSystemWithPanels {
                _type: Default::default(),
                energy_supply: "mains elec".into(),
                inverter_peak_power_ac: 1.4,
                inverter_peak_power_dc: 3.5,
                inverter_is_inside: true,
                panels: vec![PhotovoltaicPanel {
                    peak_power: 2.5,
                    ventilation_strategy: PhotovoltaicVentilationStrategy::ModeratelyVentilated,
                    pitch: 30.,
                    orientation: 0.,
                    base_height: 10.,
                    height: 1.,
                    width: 1.,
                    shading: vec![],
                }],
                inverter_type: InverterType::StringInverter,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::inverter_peak_power_ac_at_least_zero(json!({"inverter_peak_power_ac": -1})),
            case::inverter_peak_power_dc_at_least_zero(json!({"inverter_peak_power_dc": -1})),
            case::panels_at_least_one_item(json!({"panels": []})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<PhotovoltaicSystemWithPanels>(valid_example, inputs);
        }
    }

    mod photovoltaic_panel {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(PhotovoltaicPanel {
                peak_power: 1.3,
                ventilation_strategy: PhotovoltaicVentilationStrategy::Unventilated,
                pitch: 34.,
                orientation: -65.,
                base_height: 1.2,
                height: 2.3,
                width: 4.3,
                shading: vec![],
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::peak_power_at_least_zero(json!({"peak_power": -1})),
            case::base_height_at_least_zero(json!({"base_height": -1})),
            case::width_greater_than_zero(json!({"width": 0})),
            case::orientation_at_most_180(json!({"orientation360": -1})),
            case::orientation_at_least_minus_180(json!({"orientation360": 361})),
            case::pitch_at_least_zero(json!({"pitch": -1})),
            case::pitch_at_most_90(json!({"pitch": 91})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<PhotovoltaicPanel>(valid_example, inputs);
        }
    }

    mod shading_object {
        use super::*;
        use crate::hem_core::external_conditions::{ShadingObject, ShadingObjectType};

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(ShadingObject {
                object_type: ShadingObjectType::Obstacle,
                height: 10.,
                distance: 10.,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::distance_greater_than_zero(json!({"distance": 0})),
            case::height_greater_than_zero(json!({"height": 0})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<ShadingObject>(valid_example, inputs);
        }
    }

    mod shading_segment {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(ShadingSegment {
                start: 180.,
                end: -180.,
                shading_objects: vec![],
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::end_at_most_180(json!({"end360": -1})),
            case::end_at_least_minus_180(json!({"end360": 361})),
            case::start_at_most_180(json!({"start360": -1})),
            case::start_at_least_minus_180(json!({"start360": 361})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<ShadingSegment>(valid_example, inputs);
        }
    }

    mod air_conditioning {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(SpaceCoolSystemDetails::AirConditioning {
                cooling_capacity: 1.2,
                efficiency: 4.3,
                frac_convective: 0.2,
                energy_supply: "mains elec".into(),
                control: "control".into(),
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::frac_convective_at_least_zero(json!({"frac_convective": -1})),
            case::frac_convective_at_most_one(json!({"frac_convective": 2})),
            case::efficiency_at_least_zero(json!({"efficiency": -1})),
            case::cooling_capacity_at_least_zero(json!({"cooling_capacity": -1})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<SpaceCoolSystemDetails>(valid_example, inputs);
        }
    }

    mod space_heat_system {
        use super::*;

        mod instant_electric_heater {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(SpaceHeatSystemDetails::InstantElectricHeater {
                    rated_power: 6.0,
                    frac_convective: 0.3,
                    energy_supply: "mains elec".into(),
                    control: "control".into(),
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::frac_convective_at_least_zero(json!({"frac_convective": -1})),
                case::frac_convective_at_most_one(json!({"frac_convective": 2})),
                case::rated_power_at_least_zero(json!({"rated_power": -1})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<SpaceHeatSystemDetails>(valid_example, inputs);
            }
        }

        mod wet_distribution {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(SpaceHeatSystemDetails::WetDistribution {
                    heat_source: SpaceHeatSystemHeatSource {
                        name: "immersion_system_heating".into(),
                        temp_flow_limit_upper: None,
                    },
                    bypass_fraction_recirculated: Some(0.0),
                    design_flow_temp: 55,
                    emitters: vec![WetEmitter::Radiator {
                        exponent: 1.2,
                        frac_convective: 0.4,
                        constant_data: RadiatorConstantData::Constant { constant: 0.08 },
                    }],
                    pipework: vec![],
                    ecodesign_controller: EcoDesignController {
                        ecodesign_control_class: EcoDesignControllerClass::ClassVI,
                        min_outdoor_temp: None,
                        max_outdoor_temp: None,
                        min_flow_temp: None,
                    },
                    temp_diff_emit_dsgn: 10.0,
                    energy_supply: None,
                    zone: "zone_test".into(),
                    control: "control_test".into(),
                    thermal_mass: None,
                    flow_data: FlowData::Variable {
                        _variable_flow: Default::default(),
                        max_flow_rate: 18.,
                        min_flow_rate: 3.,
                    },
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::min_and_max_flow_rate_needed_if_variable_flow_true(json!({"variable_flow": true, "min_flow_rate": null, "max_flow_rate": null})),
                case::min_and_max_flow_rate_needed_if_variable_flow_true(json!({"variable_flow": true, "min_flow_rate": 4, "max_flow_rate": null})),
                case::min_and_max_flow_rate_needed_if_variable_flow_true(json!({"variable_flow": true, "min_flow_rate": null, "max_flow_rate": 4})),
                case::design_flow_rate_needed_if_variable_flow_false(json!({"variable_flow": false, "design_flow_rate": null})),
                case::bypass_fraction_recirculated_at_most_one(json!({"bypass_fraction_recirculated": 2})),
                case::bypass_fraction_recirculated_at_least_zero(json!({"bypass_fraction_recirculated": -1})),
                case::design_flow_rate_at_least_zero(json!({"variable_flow": false, "design_flow_rate": -1})),
                case::max_flow_rate_at_least_zero(json!({"max_flow_rate": -1})),
                case::min_flow_rate_at_least_zero(json!({"min_flow_rate": -1})),
                case::temp_diff_emit_dsgn_at_least_zero(json!({"temp_diff_emit_dsgn": -1})),
                case::thermal_mass_at_least_zero(json!({"thermal_mass": -1})),
                case::design_flow_temp_greater_than_zero(json!({"design_flow_temp": 0})),
                case::emitters_at_least_one_item(json!({"emitters": []})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<SpaceHeatSystemDetails>(valid_example, inputs);
            }
        }

        mod warm_air {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(SpaceHeatSystemDetails::WarmAir {
                    frac_convective: 0.8,
                    heat_source: SpaceHeatSystemHeatSource {
                        name: "immersion_system_heating".into(),
                        temp_flow_limit_upper: None,
                    },
                    control: "control".into(),
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::frac_convective_at_least_zero(json!({"frac_convective": -1})),
                case::frac_convective_at_most_one(json!({"frac_convective": 2})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<SpaceHeatSystemDetails>(valid_example, inputs);
            }
        }

        mod electric_storage_heater {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(SpaceHeatSystemDetails::ElectricStorageHeater {
                    control_charger: "control_charger_test".into(),
                    dry_core_max_output: vec![[0.0, 0.0], [0.5, 1.5], [1.0, 3.0]],
                    dry_core_min_output: vec![[0.0, 0.0], [0.5, 0.02], [1.0, 0.05]],
                    energy_supply: "energy_supply_test".into(),
                    frac_convective: 0.2,
                    n_units: 1,
                    pwr_in: 3.7,
                    state_of_charge_init: 0.,
                    rated_power_instant: 2.5,
                    storage_capacity: 20.,
                    control: "control_test".into(),
                    air_flow_type: ElectricStorageHeaterAirFlowType::FanAssisted,
                    fan_pwr: 11.0,
                    zone: "zone_test".into(),
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::first_soc_value_in_dry_core_max_output_must_be_zero(json!({"dry_core_max_output": [[0.1, 0.0], [0.5, 1.5], [1.0, 3.0]]})),
                case::last_soc_value_in_dry_core_max_output_must_be_one(json!({"dry_core_max_output": [[0.0, 0.0], [0.5, 1.5], [0.9, 3.0]]})),
                case::dry_core_max_output_values_must_be_increasing(json!({"dry_core_max_output": [[0.0, 0.0], [0.7, 1.5], [0.5, 3.0], [1.0, 3.5]]})),
                case::first_soc_value_in_dry_core_min_output_must_be_zero(json!({"dry_core_min_output": [[0.1, 0.0], [0.5, 0.02], [1.0, 0.05]]})),
                case::last_soc_value_in_dry_core_min_output_must_be_one(json!({"dry_core_min_output": [[0.0, 0.0], [0.5, 0.02], [0.9, 0.05]]})),
                case::dry_core_min_output_values_must_be_increasing(json!({"dry_core_min_output": [[0.0, 0.0], [0.7, 0.02], [0.5, 0.03], [1.0, 0.05]]})),
                case::dry_core_max_output_must_have_at_least_two_items(json!({"dry_core_max_output": [[]]})),
                case::dry_core_max_output_item_at_least_zero(json!({"dry_core_max_output": [[-1], [-1]]})),
                case::dry_core_min_output_must_have_at_least_two_items(json!({"dry_core_min_output": [[]]})),
                case::dry_core_min_output_item_at_least_zero(json!({"dry_core_min_output": [[-1], [-1]]})),
                case::fan_pwr_at_least_zero(json!({"fan_pwr": -1})),
                case::pwr_in_at_least_zero(json!({"pwr_in": -1})),
                case::rated_power_instant_at_least_zero(json!({"rated_power_instant": -1})),
                case::frac_convective_at_least_zero(json!({"frac_convective": -1})),
                case::frac_convective_at_most_one(json!({"frac_convective": 2})),
                case::n_units_greater_than_zer(json!({"n_units": 0})),
                case::storage_capacity_greater_than_zero(json!({"storage_capacity": 0})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<SpaceHeatSystemDetails>(valid_example, inputs);
            }

            #[rstest]
            fn test_validate_dry_core_output_when_no_input(valid_example: JsonValue) {
                let modified_input =
                    merge_json_onto_base(valid_example.clone(), json!({"dry_core_max_output": []}));
                let space_heat_system: SpaceHeatSystemDetails =
                    serde_json::from_value(modified_input).unwrap();
                if let SpaceHeatSystemDetails::ElectricStorageHeater {
                    dry_core_max_output,
                    ..
                } = space_heat_system
                {
                    assert!(dry_core_max_output.is_empty());
                } else {
                    panic!("ElectricStorageHeater variant expected");
                }

                let modified_input =
                    merge_json_onto_base(valid_example, json!({"dry_core_min_output": []}));
                let space_heat_system: SpaceHeatSystemDetails =
                    serde_json::from_value(modified_input).unwrap();
                if let SpaceHeatSystemDetails::ElectricStorageHeater {
                    dry_core_min_output,
                    ..
                } = space_heat_system
                {
                    assert!(dry_core_min_output.is_empty());
                } else {
                    panic!("ElectricStorageHeater variant expected");
                }
            }
        }
    }

    mod control {
        use super::*;

        mod charge_target {
            use super::*;
            use crate::core::schedule::input::Schedule;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(ControlDetails::ChargeTarget {
                    charge_level: None,
                    external_sensor: None,
                    logic_type: None,
                    schedule: Schedule {
                        main: vec![],
                        references: Default::default(),
                    },
                    temp_charge_cut: None,
                    temp_charge_cut_delta: None,
                    charge_calc_time: default_charge_calc_time(),
                    start_day: 0,
                    time_series_step: 0.1,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::logic_type_automatic_requires_temp_charge_cut(json!({"logic_type": "automatic"})),
                case::logic_type_celect_requires_temp_charge_cut(json!({"logic_type": "celect"})),
                case::logic_type_hhrsh_requires_temp_charge_cut(json!({"logic_type": "hhrsh"})),
                case::temp_charge_cut_at_least_absolute_zero(json!({"temp_charge_cut": -274})),
                case::charge_calc_time_at_least_zero(json!({"charge_calc_time": -1})),
                case::charge_calc_time_less_than_24(json!({"charge_calc_time": 24})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<ControlDetails>(valid_example, inputs);
            }
        }

        mod setpoint_timer {
            use super::*;

            #[fixture]
            fn base_data() -> JsonValue {
                json!({
                    "type": "SetpointTimeControl",
                    "start_day": 0,
                    "time_series_step": 1.0,
                    "schedule": {"main": [{"value": 21.0, "repeat": 24}]},
                })
            }

            /// Test that valid setpoint bounds pass validation.
            #[rstest]
            fn test_valid_setpoint_bounds(base_data: JsonValue) {
                let data = merge_json_onto_base(
                    base_data,
                    json!({"setpoint_min": 18.0, "setpoint_max": 25.0, "default_to_max": false}),
                );

                let control: ControlDetails = serde_json::from_value(data).unwrap();
                assert!(matches!(control, ControlDetails::SetpointTimer { .. }));
                if let ControlDetails::SetpointTimer {
                    setpoint_bounds:
                        Some(
                            SetpointBoundsInput::MinAndMax {
                                setpoint_min,
                                setpoint_max,
                                ..
                            },
                            ..,
                        ),
                    ..
                } = control
                {
                    assert_eq!(setpoint_min, 18.0);
                    assert_eq!(setpoint_max, 25.0);
                } else {
                    panic!("control was expected to be a setpoint timer type");
                }
            }

            #[rstest]
            fn test_invalid_setpoint_bounds(base_data: JsonValue) {
                let data = merge_json_onto_base(
                    base_data,
                    json!({
                        "setpoint_min": 18.0,
                        "setpoint_max": 18.0,  // Equal to min - should fail
                        "default_to_max": false, // need to specify default_to_max when min and max present
                    }),
                );

                let control: Result<ControlDetails, _> = serde_json::from_value(data);

                assert!(match control {
                    Ok(control) => control.validate().is_err(),
                    Err(_) => true,
                });
            }

            #[rstest(inputs,
                case::advanced_start_at_least_zero(json!({"advanced_start": -1})),
            )]
            fn test_validate_range_constraints(base_data: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<ControlDetails>(base_data, inputs);
            }
        }

        mod on_off_cost_minimising {
            use super::*;
            use crate::core::schedule::input::Schedule;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(ControlDetails::OnOffCostMinimising {
                    schedule: Schedule {
                        main: vec![],
                        references: Default::default(),
                    },
                    start_day: 0,
                    time_series_step: 1.,
                    time_on_daily: 1.,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::time_on_daily_greater_than_zero(json!({"time_on_daily": 0})),
                case::time_on_daily_at_most_24(json!({"time_on_daily": 25})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<ControlDetails>(valid_example, inputs);
            }
        }
    }

    mod waste_water_heat_recovery_system {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(WasteWaterHeatRecoveryDetails {
                _type: Default::default(),
                cold_water_source: "header tank".into(),
                flow_rates: vec![],
                system_a_efficiencies: vec![],
                system_a_utilisation_factor: Some(1.),
                system_b_efficiencies: None,
                system_b_utilisation_factor: None,
                system_c_efficiencies: None,
                system_c_utilisation_factor: None,
                system_b_efficiency_factor: default_system_b_efficiency_factor(),
                system_c_efficiency_factor: default_system_c_efficiency_factor(),
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::flow_rates_and_system_a_efficiencies_must_have_same_length(json!({
                    "system_a_efficiencies": [1, 2, 3],
                    "system_b_efficiencies": [1, 2, 3, 4],
                    "system_c_efficiencies": [1, 2, 3, 4],
                    "flow_rates": [1, 2, 3, 4],
                })),
            case::flow_rates_and_system_b_efficiencies_must_have_same_length(json!({
                    "system_a_efficiencies": [1, 2, 3, 4],
                    "system_b_efficiencies": [1, 2, 3],
                    "system_c_efficiencies": [1, 2, 3, 4],
                    "flow_rates": [1, 2, 3, 4],
                })),
            case::flow_rates_and_system_c_efficiencies_must_have_same_length(json!({
                    "system_a_efficiencies": [1, 2, 3, 4],
                    "system_b_efficiencies": [1, 2, 3, 4],
                    "system_c_efficiencies": [1, 2, 3],
                    "flow_rates": [1, 2, 3, 4],
                })),
            case::system_a_efficiencies_item_at_least_zero(json!({"system_a_efficiencies": [-1]})),
            case::flow_rates_item_at_least_zero(json!({"flow_rates": [-1]})),
            case::system_a_utilisation_factor_at_least_zero(json!({"system_a_utilisation_factor": -1})),
            case::system_a_utilisation_factor_at_most_one(json!({"system_a_utilisation_factor": 2})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<WasteWaterHeatRecoveryDetails>(valid_example, inputs);
        }

        /// Test that system B and C fields are truly optional.
        #[rstest]
        fn test_wwhrs_system_b_and_c_fields_optional() {
            let cold_water_source = "test_cold_water_source".into();

            let wwhrs = WasteWaterHeatRecoveryDetails {
                _type: Default::default(),
                cold_water_source,
                flow_rates: vec![5., 7., 9., 11., 13.],
                system_a_efficiencies: vec![50., 55., 60., 65., 70.],
                system_a_utilisation_factor: Some(0.9),
                system_b_efficiencies: None,
                system_b_utilisation_factor: None,
                system_c_efficiencies: None,
                system_c_utilisation_factor: None,
                system_b_efficiency_factor: default_system_b_efficiency_factor(),
                system_c_efficiency_factor: default_system_c_efficiency_factor(),
            };

            assert!(wwhrs.system_b_efficiencies.is_none());
            assert!(wwhrs.system_c_efficiencies.is_none());
            assert!(wwhrs.system_b_utilisation_factor.is_none());
            assert!(wwhrs.system_c_utilisation_factor.is_none());
        }

        /// Test WWHRS with all three systems configured.
        #[rstest]
        fn test_wwhrs_with_all_systems() {
            let cold_water_source = "test_cold_water_source".into();

            let wwhrs = WasteWaterHeatRecoveryDetails {
                _type: Default::default(),
                cold_water_source,
                flow_rates: vec![5., 7., 9., 11., 13.],
                system_a_efficiencies: vec![50., 55., 60., 65., 70.],
                system_a_utilisation_factor: Some(0.9),
                system_b_efficiencies: Some(vec![40., 44., 48., 52., 56.]),
                system_b_utilisation_factor: Some(0.85),
                system_c_efficiencies: Some(vec![45., 49., 53., 57., 61.]),
                system_c_utilisation_factor: Some(0.88),
                system_b_efficiency_factor: 0.81,
                system_c_efficiency_factor: 0.88,
            };

            // Test that all values are stored correctly
            assert_eq!(
                wwhrs.system_b_efficiencies,
                Some(vec![40., 44., 48., 52., 56.])
            );
            assert_eq!(
                wwhrs.system_c_efficiencies,
                Some(vec![45., 49., 53., 57., 61.])
            );
            assert_eq!(wwhrs.system_b_utilisation_factor, Some(0.85));
            assert_eq!(wwhrs.system_c_utilisation_factor, Some(0.88));
        }
    }

    mod building_element {
        use super::*;

        mod geometric_base {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(BuildingElement::AdjacentConditionedSpace {
                    area: None,
                    pitch: 90.,
                    u_value_input: UValueInput::UValue { u_value: 0.0 },
                    areal_heat_capacity: 0.0,
                    mass_distribution_class: MassDistributionClass::D,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::area_greater_than_zero(json!({"area": 0})),
                case::pitch_at_least_zero(json!({"pitch": -1})),
                case::pitch_at_most_zero(json!({"pitch": 181})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<BuildingElement>(valid_example, inputs);
            }
        }

        mod thermal_properties_base {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(BuildingElement::AdjacentConditionedSpace {
                    area: None,
                    pitch: 80.,
                    u_value_input: UValueInput::UValue { u_value: 0.3 },
                    areal_heat_capacity: 0.0,
                    mass_distribution_class: MassDistributionClass::D,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::height_greater_than_zero(json!({"height": 0})),
                case::width_greater_than_zero(json!({"width": 0})),
                case::thermal_resistance_construction_greater_than_zero(json!({"thermal_resistance_construction": 0})),
                case::u_value_greater_than_zero(json!({"u_value": 0})),
                case::must_specify_thermal_resistance_construction_or_u_value(json!({"thermal_resistance_construction": null, "u_value": null})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<BuildingElement>(valid_example, inputs);
            }
        }

        mod spatial_properties_base {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(BuildingElement::Opaque {
                    is_unheated_pitched_roof: None,
                    pitch: 80.,
                    orientation: None,
                    base_height: 0.0,
                    u_value_input: UValueInput::UValue { u_value: 0.3 },
                    areal_heat_capacity: 0.0,
                    mass_distribution_class: MassDistributionClass::D,
                    solar_absorption_coeff: 0.0,
                    is_external_door: None,
                    area_input: BuildingElementAreaOrHeightWidthInput {
                        area: None,
                        height_and_width: None,
                    },
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::areal_heat_capacity_greater_than_zero(json!({"areal_heat_capacity": 0})),
                case::thermal_resistance_construction_greater_than_zero(json!({"thermal_resistance_construction": 0})),
                case::orientation_at_most_180(json!({"orientation360": -1})),
                case::orientation_at_least_minus_180(json!({"orientation360": 361})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<BuildingElement>(valid_example, inputs);
            }
        }

        mod ground_base {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(BuildingElement::Ground {
                    area: 20.,
                    total_area: 15.,
                    pitch: 90.,
                    u_value: 1.4,
                    thermal_resistance_floor_construction: 0.2,
                    areal_heat_capacity: 19200.,
                    mass_distribution_class: MassDistributionClass::I,
                    perimeter: 16.,
                    psi_wall_floor_junc: 0.5,
                    thickness_walls: 0.2,
                    floor_data: FloorData::SlabNoEdgeInsulation,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::calculated_r_vi_greater_than_zero(json!({"u_value": 1, "thermal_resistance_floor_construction": 1})),
                case::area_greater_than_zero(json!({"area": 0})),
                case::thickness_walls_greater_than_zero(json!({"thickness_walls": 0})),
                case::total_area_greater_than_zero(json!({"total_area": 0})),
                case::perimeter_greater_than_zero(json!({"perimeter": 0})),
                case::areal_heat_capacity_greater_than_zero(json!({"areal_heat_capacity": 0})),
                case::thermal_resistance_floor_construction_greater_than_zero(json!({"thermal_resistance_floor_construction": 0})),
                case::u_value_greater_than_zero(json!({"u_value": 0})),
                // height_upper_surface only present in suspended floor, so moved that test there
                case::pitch_at_least_zero(json!({"pitch": -1})),
                case::pitch_at_most_zero(json!({"pitch": 181})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<BuildingElement>(valid_example, inputs);
            }
        }

        mod transparent {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(BuildingElement::Transparent {
                    u_value_input: UValueInput::ThermalResistanceConstruction {
                        thermal_resistance_construction: 0.74,
                    },
                    pitch: 45.,
                    orientation: None,
                    g_value: 10.,
                    frame_area_fraction: 0.25,
                    base_height: 10.,
                    area_input: BuildingElementAreaOrHeightWidthInput {
                        area: None,
                        height_and_width: Some(BuildingElementHeightWidthInput {
                            height: 10.,
                            width: 5.,
                        }),
                    },
                    free_area_height: 1.6,
                    mid_height: 40.,
                    max_window_open_area: 3.,
                    window_part_list: vec![],
                    shading: vec![],
                    areal_heat_capacity: default_areal_heat_capacity_for_windows(),
                    window_openable_control: None,
                    treatment: vec![],
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::max_window_open_area_at_most_area(json!({"max_window_open_area": 9999, "width": 5, "height": 10})),
                case::area_or_height_and_width_provided(json!({"area": null, "height": null, "width": 12})),
                case::area_or_height_and_width_provided(json!({"area": null, "height": 12, "width": null})),
                case::base_height_at_least_zero(json!({"base_height": -1})),
                case::free_area_height_at_least_zero(json!({"free_area_height": -1})),
                case::g_value_at_least_zero(json!({"g_value": -1})),
                case::max_window_open_area_at_least_zero(json!({"max_window_open_area": -1})),
                case::mid_height_greater_than_zero(json!({"mid_height": 0})),
                case::thermal_resistance_construction_greater_than_zero(json!({"thermal_resistance_construction": 0})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<BuildingElement>(valid_example, inputs);
            }
        }

        mod opaque {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(BuildingElement::Opaque {
                    is_unheated_pitched_roof: None,
                    solar_absorption_coeff: 0.8,
                    u_value_input: UValueInput::UValue { u_value: 1.2 },
                    pitch: 1.2,
                    orientation: None,
                    base_height: 10.,
                    area_input: BuildingElementAreaOrHeightWidthInput {
                        area: None,
                        height_and_width: None,
                    },
                    areal_heat_capacity: default_areal_heat_capacity_for_windows(),
                    mass_distribution_class: MassDistributionClass::I,
                    is_external_door: None,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::base_height_at_least_zero(json!({"base_height": -1})),
                case::solar_absorption_coeff_at_least_zero(json!({"solar_absorption_coeff": -1})),
                case::solar_absorption_coeff_at_most_one(json!({"solar_absorption_coeff": 2})),
                case::thermal_resistance_construction_greater_than_zero(json!({"thermal_resistance_construction": 0})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<BuildingElement>(valid_example, inputs);
            }
        }

        mod adjacent_conditioned {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(BuildingElement::AdjacentConditionedSpace {
                    area: None,
                    u_value_input: UValueInput::UValue { u_value: 1.2 },
                    pitch: 1.2,
                    areal_heat_capacity: default_areal_heat_capacity_for_windows(),
                    mass_distribution_class: MassDistributionClass::I,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::thermal_resistance_construction_greater_than_zero(json!({"thermal_resistance_construction": 0, "u_value": null})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<BuildingElement>(valid_example, inputs);
            }
        }

        mod adjacent_unconditioned {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(BuildingElement::AdjacentUnconditionedSpace {
                    area: None,
                    u_value_input: UValueInput::UValue { u_value: 0.7 },
                    pitch: 1.2,
                    areal_heat_capacity: default_areal_heat_capacity_for_windows(),
                    mass_distribution_class: MassDistributionClass::I,
                    thermal_resistance_unconditioned_space: 1.3,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::thermal_resistance_construction_greater_than_zero(json!({"thermal_resistance_construction": 0, "u_value": null})),
                case::thermal_resistance_unconditioned_space_greater_than_zero(json!({"thermal_resistance_unconditioned_space": 0})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<BuildingElement>(valid_example, inputs);
            }
        }

        mod ground {
            use super::*;

            mod suspended_floor {
                use super::*;

                #[fixture]
                fn valid_example() -> JsonValue {
                    serde_json::to_value(BuildingElement::Ground {
                        area: 200.,
                        total_area: 210.,
                        pitch: 1.2,
                        u_value: 0.7,
                        thermal_resistance_floor_construction: 0.2,
                        areal_heat_capacity: 1.2,
                        mass_distribution_class: MassDistributionClass::I,
                        perimeter: 100.,
                        psi_wall_floor_junc: 0.4,
                        thickness_walls: 20.,
                        floor_data: FloorData::SuspendedFloor {
                            height_upper_surface: 2.3, // missing from upstream Python definition, but giving a valid input here
                            thermal_transmission_walls: 0.8,
                            area_per_perimeter_vent: 4.2,
                            shield_fact_location: WindShieldLocation::Sheltered,
                            thermal_resistance_of_insulation: 0.7,
                        },
                    })
                    .unwrap()
                }

                #[rstest(inputs,
                    case::area_per_perimeter_vent_greater_than_zero(json!({"area_per_perimeter_vent": 0, "u_value": null})),
                    case::thermal_transm_walls_greater_than_zero(json!({"thermal_transm_walls": 0, "u_value": null})),
                    case::thermal_resist_insul_greater_than_zero(json!({"thermal_resist_insul": 0})),
                    case::height_upper_surface_greater_than_zero(json!({"height_upper_surface": 0})),
                )]
                fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                    assert_range_constraints::<BuildingElement>(valid_example, inputs);
                }
            }

            mod heated_basement {
                use super::*;

                #[fixture]
                fn valid_example() -> JsonValue {
                    serde_json::to_value(BuildingElement::Ground {
                        area: 200.,
                        total_area: 210.,
                        pitch: 1.2,
                        u_value: 0.7,
                        thermal_resistance_floor_construction: 0.2,
                        areal_heat_capacity: 1.2,
                        mass_distribution_class: MassDistributionClass::I,
                        perimeter: 100.,
                        psi_wall_floor_junc: 0.4,
                        thickness_walls: 20.,
                        floor_data: FloorData::HeatedBasement {
                            depth_basement_floor: 10.,
                            thermal_resistance_of_basement_walls: 0.15,
                        },
                    })
                    .unwrap()
                }

                #[rstest(inputs,
                    case::height_basement_walls_greater_than_zero(json!({"height_basement_walls": 0, "u_value": null})),
                    case::thermal_resist_walls_base_greater_than_zero(json!({"thermal_resist_walls_base": 0, "u_value": null})),
                    case::thermal_transm_envi_base_greater_than_zero(json!({"thermal_transm_envi_base": 0, "u_value": null})),
                    case::thermal_transm_walls_greater_than_zero(json!({"thermal_transm_walls": 0, "u_value": null})),
                    case::depth_basement_floor_at_least_zero(json!({"depth_basement_floor": -1, "u_value": null})),
                )]
                fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                    assert_range_constraints::<BuildingElement>(valid_example, inputs);
                }
            }

            mod unheated_basement {
                use super::*;

                #[fixture]
                fn valid_example() -> JsonValue {
                    serde_json::to_value(BuildingElement::Ground {
                        area: 200.,
                        total_area: 210.,
                        pitch: 1.2,
                        u_value: 0.7,
                        thermal_resistance_floor_construction: 0.2,
                        areal_heat_capacity: 1.2,
                        mass_distribution_class: MassDistributionClass::I,
                        perimeter: 100.,
                        psi_wall_floor_junc: 0.4,
                        thickness_walls: 20.,
                        floor_data: FloorData::UnheatedBasement {
                            thermal_transmittance_of_floor_above_basement: 1.2,
                            thermal_transmission_walls: 0.15,
                            depth_basement_floor: 10.,
                            height_basement_walls: 10.,
                            thermal_resistance_of_basement_walls: 0.15,
                        },
                    })
                    .unwrap()
                }

                #[rstest(inputs,
                    case::thermal_resist_walls_base_greater_than_zero(json!({"thermal_resist_walls_base": 0, "u_value": null})),
                    case::thermal_transm_envi_base_greater_than_zero(json!({"thermal_transm_envi_base": 0, "u_value": null})),
                    case::thermal_transm_walls_greater_than_zero(json!({"thermal_transm_walls": 0, "u_value": null})),
                    case::depth_basement_floor_at_least_zero(json!({"depth_basement_floor": -1, "u_value": null})),
                    case::height_basement_walls_greater_than_zero(json!({"height_basement_walls": 0, "u_value": null})),
                )]
                fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                    assert_range_constraints::<BuildingElement>(valid_example, inputs);
                }
            }
        }
    }

    mod fancoil_test_data {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(FancoilTestData {
                fan_speed_data: vec![
                    FanSpeedData {
                        temperature_diff: 80.0,
                        power_output: vec![2.7, 3.6, 5., 5.3, 6.2, 7.4],
                    },
                    FanSpeedData {
                        temperature_diff: 70.0,
                        power_output: vec![2.3, 3.1, 4.2, 4.5, 5.3, 6.3],
                    },
                    FanSpeedData {
                        temperature_diff: 60.0,
                        power_output: vec![1.9, 2.6, 3.5, 3.8, 4.4, 5.3],
                    },
                    FanSpeedData {
                        temperature_diff: 50.0,
                        power_output: vec![1.5, 2., 2.8, 3., 3.5, 4.2],
                    },
                    FanSpeedData {
                        temperature_diff: 40.0,
                        power_output: vec![1.1, 1.5, 2.05, 2.25, 2.6, 3.15],
                    },
                    FanSpeedData {
                        temperature_diff: 30.0,
                        power_output: vec![0.7, 0.97, 1.32, 1.49, 1.7, 2.09],
                    },
                    FanSpeedData {
                        temperature_diff: 20.0,
                        power_output: vec![0.3, 0.44, 0.59, 0.73, 0.8, 1.03],
                    },
                    FanSpeedData {
                        temperature_diff: 10.0,
                        power_output: vec![0., 0., 0., 0., 0., 0.],
                    },
                ],
                fan_power_w: vec![15., 19., 25., 33., 43., 56.],
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::fan_speed_lists_have_same_length(json!({
                    "fan_speed_data": [
                        {"temperature_diff": 80.0, "power_output": [2, 4, 6, 8, 9, 11]},
                        {"temperature_diff": 70.0, "power_output": [2, 4, 6, 8, 9]},
                    ]
                })),
            case::fan_power_list_same_length_as_fan_speed_lists(json!({
                    "fan_speed_data": [
                        {"temperature_diff": 80.0, "power_output": [2, 4, 6, 8, 9]},
                        {"temperature_diff": 70.0, "power_output": [2, 4, 6, 8, 9]},
                    ],
                    "fan_power_W": [15, 19, 25, 33, 43, 56],
                })),
            case::fan_power_at_least_zero(json!({"fan_power_W": [-1]})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<FancoilTestData>(valid_example, inputs);
        }
    }

    mod appliance_gains {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(ApplianceGainsDetails {
                events: None,
                standby: None,
                energy_supply: "mains elec".into(),
                gains_fraction: 0.8,
                load_shifting: None,
                priority: None,
                schedule: None,
                start_day: 0,
                time_series_step: 1.0,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::standby_at_least_zero(json!({"Standby": -1})),
            case::gains_fraction_at_most_one(json!({"gains_fraction": 2})),
            case::gains_fraction_at_least_zero(json!({"gains_fraction": -1})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<ApplianceGainsDetails>(valid_example, inputs);
        }
    }

    mod boiler_cost_schedule_hybrid {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(BoilerCostScheduleHybrid {
                cost_schedule_boiler: NumericSchedule {
                    main: vec![],
                    references: Default::default(),
                },
                cost_schedule_hp: NumericSchedule {
                    main: vec![],
                    references: Default::default(),
                },
                cost_schedule_start_day: 0,
                cost_schedule_time_series_step: 1.,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::cost_schedule_start_day_at_least_zero(json!({"cost_schedule_start_day": -1})),
            case::cost_schedule_start_day_at_most_365(json!({"cost_schedule_start_day": 366})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<BoilerCostScheduleHybrid>(valid_example, inputs);
        }
    }

    mod external_conditions {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(ExternalConditionsInput {
                air_temperatures: None,
                diffuse_horizontal_radiation: None,
                direct_beam_conversion_needed: None,
                direct_beam_radiation: None,
                latitude: None,
                longitude: None,
                shading_segments: None,
                solar_reflectivity_of_ground: None,
                wind_directions: None,
                wind_speeds: None,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::air_temperatures_should_be_at_least_absolute_zero(json!({"air_temperatures": [-274]})),
            case::latitude_should_be_at_least_minus_90(json!({"latitude": -91})),
            case::latitude_should_be_at_most_90(json!({"latitude": 91})),
            case::longitude_should_be_at_least_minus_180(json!({"latitude": -181})),
            case::longitude_should_be_at_most_180(json!({"latitude": 181})),
            case::solar_reflectivity_of_ground_should_be_at_least_zero(json!({"solar_reflectivity_of_ground": [-1]})),
            case::solar_reflectivity_of_ground_should_be_at_most_one(json!({"solar_reflectivity_of_ground": [2]})),
            case::wind_speeds_should_be_at_least_zero(json!({"wind_speeds": [-1]})),
            case::diffuse_horizontal_radiation_should_be_at_least_zero(json!({"diffuse_horizontal_radiation": [-1]})),
            case::direct_beam_radiation_should_be_at_least_zero(json!({"direct_beam_radiation": [-1]})),
            case::wind_directions_should_be_at_least_zero(json!({"wind_directions": [-1]})),
            case::wind_directions_should_be_at_most_360(json!({"wind_directions": [361]})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<ExternalConditionsInput>(valid_example, inputs);
        }

        #[rstest]
        fn test_are_all_fields_set() {
            let valid_example = ExternalConditionsInput {
                air_temperatures: Some(vec![0.; 8760]),
                diffuse_horizontal_radiation: Some(vec![0.; 8760]),
                direct_beam_conversion_needed: false.into(),
                direct_beam_radiation: Some(vec![0.; 8760]),
                latitude: Some(13.9),
                longitude: Some(34.2),
                shading_segments: Some(vec![]),
                solar_reflectivity_of_ground: Some(vec![0.; 8760]),
                wind_directions: Some(vec![0.; 8760]),
                wind_speeds: Some(vec![0.; 8760]),
            };
            assert!(valid_example.are_all_fields_set());

            let invalid_example = ExternalConditionsInput {
                air_temperatures: None,
                diffuse_horizontal_radiation: None,
                direct_beam_conversion_needed: None,
                direct_beam_radiation: None,
                latitude: None,
                longitude: None,
                shading_segments: None,
                solar_reflectivity_of_ground: None,
                wind_directions: None,
                wind_speeds: None,
            };
            assert!(!invalid_example.are_all_fields_set());
        }

        #[fixture]
        fn base_for_length_test() -> JsonValue {
            json!({
                "air_temperatures": vec![0.0; 8760],
                "wind_speeds": vec![0.0; 8760],
                "wind_directions": vec![0.0; 8760],
                "diffuse_horizontal_radiation": vec![0.0; 8760],
                "direct_beam_conversion_needed": false,
                "direct_beam_radiation": vec![0.0; 8760],
                "latitude": 13.9,
                "longitude": 34.2,
                "shading_segments": [],
                "solar_reflectivity_of_ground": vec![0.0; 8760],
            })
        }

        #[rstest(inputs,
            case::fail_when_air_temperatures_wrong_length(json!({"air_temperatures": [0]})),
            case::fail_when_wind_speeds_wrong_length(json!({"wind_speeds": [0]})),
            case::fail_when_wind_directions_wrong_length(json!({"wind_directions": [0]})),
            case::fail_when_diffuse_horizontal_radiation_wrong_length(json!({"diffuse_horizontal_radiation": [0]})),
            case::fail_when_direct_beam_radiation_wrong_length(json!({"direct_beam_radiation": [0]})),
            case::fail_when_solar_reflectivity_of_ground_wrong_length(json!({"solar_reflectivity_of_ground": [0]})),
        )]
        fn test_are_all_fields_set_invalid_lengths(
            base_for_length_test: JsonValue,
            inputs: JsonValue,
        ) {
            let modified_input = merge_json_onto_base(base_for_length_test, inputs);
            let external_conditions: ExternalConditionsInput =
                serde_json::from_value(modified_input).unwrap();
            assert!(!external_conditions.are_all_fields_set());
        }
    }

    mod heat_pump_buffer_tank {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(HeatPumpBufferTank {
                daily_losses: 1.68,
                volume: 40.,
                pump_fixed_flow_rate: 15.,
                pump_power_at_flow_rate: 0.04,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::daily_losses_at_least_zero(json!({"daily_losses": -1})),
            case::pump_fixed_flow_rate_at_least_zero(json!({"pump_fixed_flow_rate": -1})),
            case::pump_power_at_flow_rate_at_least_zero(json!({"pump_power_at_flow_rate": -1})),
            case::volume_at_least_zero(json!({"volume": -1})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<HeatPumpBufferTank>(valid_example, inputs);
        }
    }

    mod heat_pump_hot_water_only_test_datum {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(HeatPumpHotWaterOnlyTestDatum {
                cop_dhw: 2.5,
                hw_tapping_prof_daily_total: 5.8,
                energy_input_measured: 2.4,
                power_standby: 0.012,
                hw_vessel_loss_daily: 2.0,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::cop_dhw_at_least_zero(json!({"cop_dhw": -1})),
            case::energy_input_measured_at_least_zero(json!({"energy_input_measured": -1})),
            case::hw_tapping_prof_daily_total_at_least_zero(json!({"hw_tapping_prof_daily_total": -1})),
            case::hw_vessel_loss_daily_at_least_zero(json!({"hw_vessel_loss_daily": -1})),
            case::power_standby_at_least_zero(json!({"power_standby": -1})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<HeatPumpHotWaterOnlyTestDatum>(valid_example, inputs);
        }
    }

    mod infiltration_ventilation {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(InfiltrationVentilation {
                control_vent_adjust_max: None,
                control_vent_adjust_min: None,
                control_window_adjust: None,
                leaks: VentilationLeaks {
                    ventilation_zone_height: 6.,
                    test_pressure: 50.,
                    test_result: 1.2,
                    env_area: 220.,
                },
                mechanical_ventilation: Default::default(),
                vents: Default::default(),
                ach_max_static_calcs: None,
                ach_min_static_calcs: None,
                altitude: 10.,
                cross_vent_possible: true,
                shield_class: VentilationShieldClass::Shielded,
                terrain_class: TerrainClass::Urban,
                ventilation_zone_base_height: 10.,
                vent_opening_ratio_init: None,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::ach_max_static_calcs_at_least_zero(json!({"ach_max_static_calcs": -1})),
            case::ventilation_zone_base_height_at_least_zero(json!({"ventilation_zone_base_height": -1})),
            case::ach_min_static_calcs_at_least_zero(json!({"ach_min_static_calcs": -1})),
            case::vent_opening_ratio_init_at_least_zero(json!({"vent_opening_ratio_init": -1})),
            case::vent_opening_ratio_init_at_most_zero(json!({"vent_opening_ratio_init": 2})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<InfiltrationVentilation>(valid_example, inputs);
        }
    }

    mod zone {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(ZoneInput {
                building_elements: Default::default(),
                space_cool_system: Default::default(),
                space_heat_system: Default::default(),
                thermal_bridging: ThermalBridging::Number(1.9),
                area: 200.,
                temp_setpnt_basis: None,
                temp_setpnt_init: 10.,
                volume: 100.,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::area_greater_than_zero(json!({"area": 0})),
            case::volume_greater_than_zero(json!({"volume": 0})),
            case::temp_setpnt_init_at_least_absolute_zero(json!({"temp_setpnt_init": -274})),
            case::space_heat_system_list_does_not_allow_duplicates(json!({"SpaceHeatSystem": ["mains", "mains", "other"]})),
            case::space_cool_system_list_does_not_allow_duplicates(json!({"SpaceCoolSystem": ["mains", "mains", "other"]})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<ZoneInput>(valid_example, inputs);
        }
    }

    mod window_treatment {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(WindowTreatment {
                treatment_type: WindowTreatmentType::Curtains,
                controls: WindowTreatmentControl::AutoMotorised,
                delta_r: 10.,
                trans_red: 0.3,
                closing_irradiance_control: None,
                opening_irradiance_control: None,
                open_control: None,
                is_open: None,
                opening_delay_hrs: 0.0, // valid value for field we are marking as required
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::trans_red_at_least_zero(json!({"trans_red": -1})),
            case::trans_red_at_most_one(json!({"trans_red": 2})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<WindowTreatment>(valid_example, inputs);
        }
    }

    mod air_terminal_device {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(AirTerminalDevice {
                area_cm2: 100.,
                pressure_difference_ref: 1.2,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::area_cm2_greater_than_zero(json!({"area_cm2": 0})),
            case::pressure_difference_ref_at_least_zero(json!({"pressure_difference_ref": -1})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<AirTerminalDevice>(valid_example, inputs);
        }
    }

    mod edge_insulation {
        use super::*;

        mod horizontal {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(EdgeInsulation::Horizontal {
                    width: 10.,
                    edge_thermal_resistance: 1.2,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::edge_thermal_resistance_greater_than_zero(json!({"edge_thermal_resistance": 0})),
                case::width_at_least_zero(json!({"width": -1})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<EdgeInsulation>(valid_example, inputs);
            }
        }

        mod vertical {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(EdgeInsulation::Vertical {
                    depth: 10.,
                    edge_thermal_resistance: 1.2,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::edge_thermal_resistance_greater_than_zero(json!({"edge_thermal_resistance": 0})),
                case::depth_at_least_zero(json!({"depth": -1})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<EdgeInsulation>(valid_example, inputs);
            }
        }
    }

    mod window_part {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(WindowPart {
                mid_height_air_flow_path: 1.5,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::mid_height_air_flow_path_greater_than_zero(json!({"mid_height_air_flow_path": 0})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<WindowPart>(valid_example, inputs);
        }
    }

    mod heat_pump_test_datum {
        use super::*;

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(HeatPumpTestDatum {
                air_flow_rate: None,
                capacity: 10.,
                cop: 1.2,
                design_flow_temp: 12.,
                eahp_mixed_ext_air_ratio: None,
                temp_outlet: 12.,
                temp_source: 13.,
                temp_test: 14.,
                test_letter: TestLetter::A,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::air_flow_rate_at_least_zero(json!({"air_flow_rate": -1})),
            case::capacity_at_least_zero(json!({"capacity": -1})),
            case::cop_at_least_zero(json!({"cop": -1})),
            case::eahp_mixed_ext_air_ratio_at_least_zero(json!({"eahp_mixed_ext_air_ratio": -1})),
            case::eahp_mixed_ext_air_ratio_at_most_one(json!({"eahp_mixed_ext_air_ratio": 2})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<HeatPumpTestDatum>(valid_example, inputs);
        }
    }

    // test for ScheduleRepeater unnecessary as `repeat` field uses usize type which cannot be negative

    // test for ControlCombinations unnecessary as asserted properties are enforced at type level
    // same applies for Schedule

    mod internal_gains {
        use super::*;

        #[rstest]
        fn test_cold_water_losses_access() {
            let internal_gains: InternalGains = serde_json::from_value(json!({
                "ColdWaterLosses": {"start_day": 1, "time_series_step": 1, "schedule": {"main": []}}
            }))
            .unwrap();

            assert!(internal_gains.cold_water_losses.is_some());
        }

        #[rstest]
        fn test_evaporative_losses_access() {
            let internal_gains: InternalGains = serde_json::from_value(json!({
                "EvaporativeLosses": {
                    "start_day": 1,
                    "time_series_step": 1,
                    "schedule": {"main": []},
                }
            }))
            .unwrap();

            assert!(internal_gains.evaporative_losses.is_some());
        }

        #[rstest]
        fn test_metabolic_gains_access() {
            let internal_gains: InternalGains = serde_json::from_value(json!({
                "metabolic gains": {"start_day": 1, "time_series_step": 1, "schedule": {"main": []}}
            }))
            .unwrap();

            assert!(internal_gains.metabolic_gains.is_some());
        }

        #[rstest]
        fn test_other_access() {
            let internal_gains: InternalGains = serde_json::from_value(
                json!({"other": {"start_day": 1, "time_series_step": 1, "schedule": {"main": []}}}),
            )
            .unwrap();

            assert!(internal_gains.other.is_some());
        }

        #[rstest]
        fn test_total_internal_gains_access() {
            let internal_gains: InternalGains = serde_json::from_value(json!({
                "total internal gains": {
                    "start_day": 1,
                    "time_series_step": 1,
                    "schedule": {"main": []},
                }
            }))
            .unwrap();
            assert!(internal_gains.total_internal_gains.is_some());

            let internal_gains: InternalGains = serde_json::from_value(json!({
                "total_internal_gains": {
                    "start_day": 1,
                    "time_series_step": 1,
                    "schedule": {"main": []},
                }
            }))
            .unwrap();
            assert!(internal_gains.total_internal_gains.is_some());
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

    #[rstest]
    fn test_reset_internal_gains() {
        let base_input = json!({
            "InternalGains": {
                "metabolic gains": {
                    "start_day": 0,
                    "time_series_step": 1,
                    "schedule": {
                        "main": [1305.6, 1876.8, 2978.4, 2121.6, 3631.2, 2284.8, 4161.6, 3304.8]
                    }
                }
            }
        });
        let mut input = InputForProcessing { input: base_input };
        input.reset_internal_gains().unwrap();
        assert_eq!(input.input, json!({"InternalGains": {}}));
    }
}
