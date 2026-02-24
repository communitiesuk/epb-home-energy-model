#![allow(unused_variables)]

use crate::core::heating_systems::heat_pump::TestLetter;
use crate::core::schedule::{BooleanSchedule, NumericSchedule};
use crate::core::units::{calculate_thermal_resistance_of_virtual_layer, Orientation360};
use crate::external_conditions::{ShadingSegment, WindowShadingObject};
use crate::read_weather_file::ExternalConditions as ExternalConditionsFromFile;
use crate::simulation_time::SimulationTime;
use crate::HEM_VERSION;
use anyhow::{anyhow, bail};
use indexmap::IndexMap;
use itertools::Itertools;
use jsonschema::Validator;
use monostate::MustBe;
use serde::{Deserialize, Serialize};
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

const HOURS_IN_YEAR: usize = 8760;

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(rename_all = "PascalCase", deny_unknown_fields)]
#[validate(custom = validate_shower_waste_water_heat_recovery_systems)]
#[validate(custom = validate_exhaust_air_heat_pump_ventilation_compatibility)]
#[validate(custom = validate_time_series)]
pub struct Input {
    /// Metadata for the input file
    #[serde(rename = "metadata", skip_serializing_if = "Option::is_none")]
    #[validate]
    metadata: Option<InputMetadata>,

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

impl Input {
    pub fn space_heat_system(&self) -> Option<&SpaceHeatSystem> {
        self.space_heat_system.as_ref()
    }

    pub fn hot_water_source(&self) -> &HotWaterSource {
        &self.hot_water_source
    }
}

#[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
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
                    if exhaust_air_source_types.contains(source_type) {
                        Some((name.into(), *source_type))
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

fn validate_time_series(input: &Input) -> Result<(), serde_valid::validation::Error> {
    let Input {
        cold_water_source,
        simulation_time,
        ..
    } = input;

    for cold_water_source in cold_water_source.values() {
        let total_steps = (simulation_time.end_time() - simulation_time.start_time()).ceil()
            / cold_water_source.time_series_step;
        if (cold_water_source.temperatures.len() as f64) < total_steps {
            return custom_validation_error("ColdWaterSource.temperatures does not contain enough values to cover the simulation.".to_string());
        }
    }

    Ok(())
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
    #[validate(custom = validate_all_items_in_option)]
    pub(crate) wind_directions: Option<Vec<Orientation360>>,

    /// List of wind speeds, one entry per hour (unit: m/s)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(custom = validate_all_items_in_option_non_negative)]
    pub(crate) wind_speeds: Option<Vec<f64>>,
}

impl ExternalConditionsInput {
    /// Assert that all required fields are set and valid, allowing the model to run without a weather file.
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
                    air_temperatures.len(),
                    diffuse_horizontal_radiation.len(),
                    direct_beam_radiation.len(),
                    solar_reflectivity_of_ground.len(),
                    wind_directions.len(),
                    wind_speeds.len(),
                ]
                .iter()
                .all(|&items| items >= HOURS_IN_YEAR)
            } else {
                false
            }
    }
}

impl From<ExternalConditionsFromFile> for ExternalConditionsInput {
    fn from(weather_file_conditions: ExternalConditionsFromFile) -> Self {
        ExternalConditionsInput {
            air_temperatures: weather_file_conditions.air_temperatures.into(),
            diffuse_horizontal_radiation: weather_file_conditions
                .diffuse_horizontal_radiation
                .into(),
            direct_beam_conversion_needed: None,
            direct_beam_radiation: weather_file_conditions.direct_beam_radiation.into(),
            latitude: None,
            longitude: None,
            shading_segments: None,
            solar_reflectivity_of_ground: weather_file_conditions
                .solar_reflectivity_of_ground
                .into(),
            wind_directions: weather_file_conditions.wind_directions.into(),
            wind_speeds: weather_file_conditions.wind_speeds.into(),
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
    #[validate(exclusive_minimum = 0.)]
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
    #[validate(exclusive_minimum = 0.)]
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
    #[validate(exclusive_minimum = 0.)]
    pub duration: f64,

    /// Electrical power consumption during the appliance event (unit: W)
    #[serde(rename = "demand_W")]
    #[validate(exclusive_minimum = 0.)]
    pub demand_w: f64,
}

pub(crate) type EnergySupplyInput = IndexMap<std::string::String, EnergySupplyDetails>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields, rename_all = "snake_case")]
pub struct EnergySupplyDetails {
    /// Type of fuel
    pub fuel: FuelType,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) diverter: Option<EnergyDiverter>,

    /// Indicates that an electric battery is present
    #[serde(rename = "ElectricBattery", skip_serializing_if = "Option::is_none")]
    pub(crate) electric_battery: Option<ElectricBattery>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub factor: Option<CustomEnergySourceFactor>,

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
        if !values.iter().all(|&v| (0. ..=1.).contains(&v)) {
            return custom_validation_error("Some threshold values for an energy supply contained numbers that were not fractions between 0 and 1 inclusive.".to_string());
        }
    }

    Ok(())
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

#[derive(Clone, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
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

    #[serde(untagged)]
    Other(String),
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
    #[validate(exclusive_minimum = 0.)]
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
    #[validate(exclusive_minimum = 0.)]
    pub maximum_charge_rate_one_way_trip: f64,

    /// The minimum charge rate one way trip the battery allows (unit: kW)
    #[validate(exclusive_minimum = 0.)]
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

pub type ColdWaterSourceInput = IndexMap<std::string::String, ColdWaterSourceDetails>;

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
    #[validate(custom = validate_running_water_temperatures)]
    pub temperatures: Vec<f64>,

    /// Duration in hours (must be within 24-hour period)
    #[validate(minimum = 0.)]
    #[validate(maximum = 24.)]
    pub time_series_step: f64,
}

fn validate_running_water_temperatures(
    temps: &[f64],
) -> Result<(), serde_valid::validation::Error> {
    if !temps.iter().all(|&t| (0. ..=100.).contains(&t)) {
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
pub struct SmartApplianceBattery {
    /// Dictionary of lists containing the battery state of charge for each timestep for each energy supply
    #[validate(custom = validate_battery_state_fractions)]
    #[validate(custom = |v| validate_all_sublists_non_empty(v, "SmartApplianceBattery"))]
    #[validate(custom = |v| validate_map_non_empty(v, "SmartApplianceBattery"))]
    pub battery_state_of_charge: IndexMap<String, Vec<f64>>,

    /// Dictionary of lists containing energy sent to the battery from generation for each timestep for each energy supply (unit: kWh)
    #[serde(default)]
    #[validate(custom = |v| validate_all_sublists_non_empty(v, "SmartApplianceBattery"))]
    #[validate(custom = |v| validate_map_non_empty(v, "SmartApplianceBattery"))]
    pub energy_into_battery_from_generation: IndexMap<String, Vec<f64>>,

    /// Dictionary of lists containing energy sent to the battery from the grid for each timestep for each energy supply (unit: kWh)
    #[serde(default)]
    #[validate(custom = |v| validate_all_sublists_non_empty(v, "SmartApplianceBattery"))]
    #[validate(custom = |v| validate_map_non_empty(v, "SmartApplianceBattery"))]
    pub energy_into_battery_from_grid: IndexMap<String, Vec<f64>>,

    /// Dictionary of lists containing energy drawn from the battery for each timestep for each energy supply (unit: kWh)
    #[serde(default)]
    #[validate(custom = |v| validate_all_sublists_non_empty(v, "SmartApplianceBattery"))]
    #[validate(custom = |v| validate_map_non_empty(v, "SmartApplianceBattery"))]
    pub energy_out_of_battery: IndexMap<String, Vec<f64>>,
}

fn validate_battery_state_fractions(
    data: &IndexMap<String, Vec<f64>>,
) -> Result<(), serde_valid::validation::Error> {
    data
        .values()
        .flatten()
        .all(|&fraction| (0. ..=1.).contains(&fraction))
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

pub type HotWaterSource = IndexMap<std::string::String, HotWaterSourceDetails>;

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
#[validate(custom = validate_dhw_tests_inputs)]
pub enum HotWaterSourceDetails {
    StorageTank {
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: String,

        /// Map of heating systems connected to the storage tank
        #[serde(rename = "HeatSource")]
        #[validate]
        heat_source: IndexMap<std::string::String, HeatSource>,

        /// Measured standby losses due to cylinder insulation at standardised conditions (unit: kWh/24h)
        #[validate(exclusive_minimum = 0.)]
        daily_losses: f64,

        /// Surface area of the heat exchanger within the storage tank (unit: m²)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(exclusive_minimum = 0.)]
        heat_exchanger_surface_area: Option<f64>,

        /// Initial temperature of the storage tank at the start of simulation (unit: ˚C)
        #[validate(minimum = 0.)]
        #[validate(maximum = 100.)]
        init_temp: f64,

        /// List of primary pipework components connected to the storage tank
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate]
        primary_pipework: Option<Vec<WaterPipework>>,

        /// Total volume of tank (unit: litre)
        #[validate(exclusive_minimum = 0.)]
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

        /// Storage loss factor 1 for combi boiler efficiency calculations (unit: kWh/day)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        storage_loss_factor_1: Option<f64>,

        /// Storage loss factor 2 for combi boiler efficiency calculations (unit: kWh/day)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        storage_loss_factor_2: Option<f64>,

        /// Rejected energy factor 3 for combi boiler efficiency calculations (dimensionless)
        #[serde(skip_serializing_if = "Option::is_none")]
        rejected_factor_3: Option<f64>,

        /// Temperature setpoint for the combi boiler hot water output (unit: ˚C)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        #[validate(maximum = 100.)]
        setpoint_temp: Option<f64>,

        /// Daily hot water usage for the combi boiler system (unit: litre/day)
        #[serde(rename = "daily_HW_usage")]
        #[validate(exclusive_minimum = 0.)]
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
        #[validate(minimum = 0.)]
        #[validate(maximum = 100.)]
        setpoint_temp: Option<f64>,
    },
    PointOfUse {
        /// Thermal efficiency of the point-of-use water heater (dimensionless, 0-1)
        #[validate(exclusive_minimum = 0.)]
        #[validate(maximum = 1.)]
        efficiency: Option<f64>,

        #[serde(rename = "EnergySupply")]
        energy_supply: String,

        #[serde(rename = "ColdWaterSource")]
        cold_water_source: String,

        /// Temperature setpoint for the point-of-use water heater output (unit: ˚C)
        #[validate(minimum = 0.)]
        #[validate(maximum = 100.)]
        setpoint_temp: f64,
    },
    SmartHotWaterTank {
        /// Total volume of tank (unit: litre)
        #[validate(exclusive_minimum = 0.)]
        volume: f64,

        /// Electrical power consumption of the pump (unit: kW)
        #[serde(rename = "power_pump_kW")]
        #[validate(exclusive_minimum = 0.)]
        power_pump_kw: f64,

        /// Maximum flow rate of the pump (unit: litre/minute)
        #[validate(exclusive_minimum = 0.)]
        max_flow_rate_pump_l_per_min: f64,

        /// Temperature below which water is considered unusable (unit: ˚C)
        #[validate(minimum = 0.)]
        #[validate(maximum = 100.)]
        temp_usable: f64,

        /// Reference to a control schedule of maximum state of charge values
        temp_setpnt_max: String,

        /// Daily standby losses due to tank insulation at standardised conditions (unit: kWh/24h)
        #[validate(exclusive_minimum = 0.)]
        daily_losses: f64,

        /// Initial temperature of the smart hot water tank at the start of simulation (unit: ˚C)
        #[validate(minimum = 0.)]
        #[validate(maximum = 100.)]
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
        #[validate(minimum = 0.)]
        #[validate(maximum = 100.)]
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

    pub(crate) fn cold_water_source(&self) -> &str {
        match self {
            HotWaterSourceDetails::StorageTank {
                cold_water_source, ..
            }
            | HotWaterSourceDetails::CombiBoiler {
                cold_water_source, ..
            }
            | HotWaterSourceDetails::Hiu {
                cold_water_source, ..
            }
            | HotWaterSourceDetails::PointOfUse {
                cold_water_source, ..
            }
            | HotWaterSourceDetails::SmartHotWaterTank {
                cold_water_source, ..
            }
            | HotWaterSourceDetails::HeatBattery {
                cold_water_source, ..
            } => cold_water_source,
        }
    }

    pub fn contains_heat_source(&self, source_key: &str) -> bool {
        match self {
            HotWaterSourceDetails::StorageTank { heat_source, .. } => {
                heat_source.contains_key(source_key)
            }
            HotWaterSourceDetails::SmartHotWaterTank { heat_source, .. } => {
                heat_source.contains_key(source_key)
            }
            _ => false,
        }
    }

    /// Provides access to check if provided string can be found within the name of a heat source wet, if applicable
    pub fn contains_heat_source_wet_reference(&self, source_key_part: &str) -> bool {
        match self {
            HotWaterSourceDetails::CombiBoiler {
                heat_source_wet, ..
            } => heat_source_wet.contains(source_key_part),
            HotWaterSourceDetails::Hiu {
                heat_source_wet, ..
            } => heat_source_wet.contains(source_key_part),
            HotWaterSourceDetails::HeatBattery {
                heat_source_wet, ..
            } => heat_source_wet.contains(source_key_part),
            _ => false,
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

fn validate_dhw_tests_inputs(
    hw_source_details: &HotWaterSourceDetails,
) -> Result<(), serde_valid::validation::Error> {
    let (
        separate_dhw_tests,
        rejected_energy_1,
        rejected_factor_3,
        storage_loss_factor_1,
        storage_loss_factor_2,
    ) = if let HotWaterSourceDetails::CombiBoiler {
        separate_dhw_tests,
        rejected_energy_1,
        rejected_factor_3,
        storage_loss_factor_1,
        storage_loss_factor_2,
        ..
    } = hw_source_details
    {
        (
            separate_dhw_tests,
            rejected_energy_1,
            rejected_factor_3,
            storage_loss_factor_1,
            storage_loss_factor_2,
        )
    } else {
        return Ok(());
    };
    if matches!(
        separate_dhw_tests,
        BoilerHotWaterTest::ML | BoilerHotWaterTest::MS
    ) {
        if rejected_energy_1.is_none()
            || rejected_factor_3.is_none()
            || storage_loss_factor_2.is_none()
        {
            return custom_validation_error("Loss factors r1, F2, and F3 are required when a combi boiler is tested to two profiles.".into());
        } else if storage_loss_factor_1.is_some() {
            return custom_validation_error(
                "storage_loss_factor_1 invalid input for combis tested to two profiles.".into(),
            );
        }
    } else if matches!(
        separate_dhw_tests,
        BoilerHotWaterTest::MOnly | BoilerHotWaterTest::NoAdditionalTests
    ) {
        if rejected_energy_1.is_none() || storage_loss_factor_1.is_none() {
            return custom_validation_error(
                "Loss factors r1, and F1, are required when a combi boiler is tested to profile M, or not tested."
                    .into(),
            );
        } else if storage_loss_factor_2.is_some() {
            return custom_validation_error(
                "storage_loss_factor_2 invalid input for combis tested to one profile, or not tested."
                    .into(),
            );
        } else if rejected_factor_3.is_some() {
            return custom_validation_error(
                "rejected_factor_3 invalid input for combis tested to one profile, or not tested."
                    .into(),
            );
        }
    }

    Ok(())
}

// fn validate_dry_core_output(
//     output_data: &[[f64; 2]],
//     field: &str,
// ) -> Result<(), serde_valid::validation::Error> {
//     //ensure body of data has at least 2 pairs
//     if output_data.len() < 2 {
//         return custom_validation_error(format!(
//             "The field {field} for an electric storage heater must have at least 2 pairs of data."
//         ));
//     }
//
//     // Convert ESH_***_output to NumPy arrays without sorting
//     let soc_values = output_data.iter().map(|f| f[0]).collect_vec();
//
//     // Validate that SOC array is in strictly increasing order
//     if !soc_values.iter().tuple_windows().all(|(a, b)| a <= b) {
//         return custom_validation_error(format!(
//             "{field} SOC values must be in increasing order (from 0.0 to 1.0)."
//         ));
//     }
//
//     // Validate that both SOC arrays start at 0.0 and end at 1.0
//     if !is_close!(*soc_values.first().unwrap(), 0.) {
//         return custom_validation_error(format!(
//             "The first SOC value in {field} must be 0.0 (fully discharged)."
//         ));
//     }
//
//     if !is_close!(*soc_values.last().unwrap(), 1.) {
//         return custom_validation_error(format!(
//             "The last SOC value in {field} must be 1.0 (fully charged)."
//         ));
//     }
//
//     Ok(())
// }

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
pub enum HeatSource {
    ImmersionHeater {
        /// (unit: kW)
        #[validate(exclusive_minimum = 0.)]
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
        #[validate(exclusive_minimum = 0.)]
        area_module: f64,

        /// Number of collector modules installed
        #[validate(minimum = 1)]
        modules: usize,

        #[validate(exclusive_minimum = 0.)]
        #[validate(maximum = 1.)]
        peak_collector_efficiency: f64,

        /// Hemispherical incidence angle modifier
        #[validate(exclusive_minimum = 0.)]
        #[validate(maximum = 1.)]
        incidence_angle_modifier: f64,

        /// First order heat loss coefficient
        #[validate(exclusive_minimum = 0.)]
        first_order_hlc: f64,

        /// Second order heat loss coefficient
        #[validate(minimum = 0.)]
        second_order_hlc: f64,

        /// Mass flow rate solar loop (unit: kg/s)
        #[validate(exclusive_minimum = 0.)]
        collector_mass_flow_rate: f64,

        /// Power of collector pump (unit: kW)
        #[validate(minimum = 0.)]
        power_pump: f64,

        /// Power of collector pump controller (unit: kW)
        #[validate(minimum = 0.)]
        power_pump_control: f64,

        #[serde(rename = "EnergySupply")]
        energy_supply: String,

        /// Tilt angle (inclination) of the solar thermal panel from horizontal,
        /// measured upwards facing, 0 to 90, in degrees.
        /// 0=horizontal surface, 90=vertical surface.
        /// Needed to calculate solar irradiation at the panel surface.
        #[validate(minimum = 0.)]
        #[validate(maximum = 90.)]
        tilt: f64,

        #[validate]
        orientation360: Orientation360,

        /// Heat loss coefficient of the collector loop piping (unit: W/K)
        #[validate(exclusive_minimum = 0.)]
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
        #[validate(exclusive_minimum = 0.)]
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
        #[validate(exclusive_minimum = 0.)]
        power_max: f64,

        /// Annual average hot water use for the dwelling (unit: litres/day)
        #[validate(exclusive_minimum = 0.)]
        vol_hw_daily_average: f64,

        /// Tank volume stored in the database (unit: litres)
        #[validate(exclusive_minimum = 0.)]
        tank_volume_declared: f64,

        /// Surface area of heat exchanger stored in the database (unit: m2)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(exclusive_minimum = 0.)]
        heat_exchanger_surface_area_declared: Option<f64>,

        /// Standing heat loss (unit: kWh/day)
        #[validate(exclusive_minimum = 0.)]
        daily_losses_declared: f64,

        /// In use factor to be applied to heat pump efficiency
        #[validate(exclusive_minimum = 0.)]
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
pub enum SolarCollectorLoopLocation {
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
pub struct HeatPumpHotWaterTestData {
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
    #[validate(exclusive_minimum = 0.)]
    pub(crate) cop_dhw: f64,

    /// daily energy requirement (kWh/day) for tapping profile used for test
    #[validate(exclusive_minimum = 0.)]
    pub(crate) hw_tapping_prof_daily_total: f64,

    /// electrical input energy (kWh) measured in EN 16147 test over 24 hrs
    #[validate(exclusive_minimum = 0.)]
    pub(crate) energy_input_measured: f64,

    /// standby power (W) measured in EN 16147 test
    #[validate(exclusive_minimum = 0.)]
    pub(crate) power_standby: f64,

    /// daily hot water vessel heat loss
    /// (kWh/day) for a 45 K temperature difference between vessel
    /// and surroundings, tested in accordance with BS 1566 or
    /// EN 12897 or any equivalent standard. Vessel must be same
    /// as that used during EN 16147 test (unit: kWh/day)
    #[validate(exclusive_minimum = 0.)]
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
    #[validate(exclusive_minimum = 0.)]
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
        #[validate(exclusive_minimum = 0.)]
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
        #[validate(exclusive_minimum = 0.)]
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
pub struct Baths(#[validate] pub IndexMap<std::string::String, BathDetails>);

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct BathDetails {
    /// Volume held by bath (unit: litre)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) size: f64,

    #[serde(rename = "ColdWaterSource")]
    pub(crate) cold_water_source: String,

    /// Reference to HotWaterSource object that provides hot water to this bath. If only one HotWaterSource is defined, then this will be assumed by default
    #[serde(rename = "HotWaterSource", skip_serializing_if = "Option::is_none")]
    pub(crate) hot_water_source: Option<String>,

    /// Tap/outlet flow rate (unit: litre/minute)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) flowrate: f64,
}

#[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct OtherWaterUses(#[validate] pub IndexMap<std::string::String, OtherWaterUse>);

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct OtherWaterUse {
    /// Tap/outlet flow rate (unit: litre/minute)
    #[validate(exclusive_minimum = 0.)]
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
pub enum WaterDistribution {
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

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct WaterHeatingEvent {
    /// Hours from start of simulation to when the water heating event begins (unit: hours)
    #[validate(minimum = 0.)]
    pub start: f64,

    /// Duration of the water heating event (unit: minutes)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(exclusive_minimum = 0.)]
    pub duration: Option<f64>,

    /// Volume of water for the event (unit: litre)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(exclusive_minimum = 0.)]
    pub volume: Option<f64>,

    /// Target temperature for the water heating event (unit: ˚C)
    #[validate(minimum = 0.)]
    #[validate(maximum = 100.)]
    pub temperature: f64,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum WaterHeatingEventType {
    Shower,
    Bath,
    Other,
}

pub type SpaceHeatSystem = IndexMap<std::string::String, SpaceHeatSystemDetails>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "type")]
pub enum SpaceHeatSystemDetails {
    #[serde(rename = "InstantElecHeater")]
    InstantElectricHeater {
        //// Rated power of the instant electric heater. (Unit: kW)
        #[validate(exclusive_minimum = 0.)]
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
        #[validate(exclusive_minimum = 0.)]
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
        #[validate(exclusive_maximum = 1.)]
        bypass_fraction_recirculated: Option<f64>,

        /// Design flow temperature. (Unit: ˚C)
        #[validate(exclusive_minimum = 0.)]
        design_flow_temp: f64,

        /// Wet emitter details of the heating system.
        #[validate(min_items = 1)]
        #[validate]
        emitters: Vec<WetEmitter>,

        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        #[validate]
        pipework: Vec<WaterPipework>,

        #[validate]
        ecodesign_controller: EcoDesignController,

        /// Design temperature difference across the emitters. (Unit: deg C or K)
        #[validate(exclusive_minimum = 0.)]
        temp_diff_emit_dsgn: f64,

        #[serde(skip_serializing_if = "Option::is_none")]
        /// Thermal mass of the emitters. (Unit: kWh/K)
        #[validate(exclusive_minimum = 0.)]
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
pub enum FlowData {
    Variable {
        #[serde(rename = "variable_flow")]
        #[cfg_attr(feature = "arbitrary", arbitrary(value = MustBe!(true)))]
        _variable_flow: MustBe!(true),

        /// Maximum flow rate allowed (unit: litres/min)
        #[validate(exclusive_minimum = 0.)]
        max_flow_rate: f64,

        /// Minimum flow rate allowed (unit: litres/min)
        #[validate(exclusive_minimum = 0.)]
        min_flow_rate: f64,
    },
    Design {
        #[serde(rename = "variable_flow")]
        #[cfg_attr(feature = "arbitrary", arbitrary(value = MustBe!(false)))]
        _variable_flow: MustBe!(false),

        /// Constant flow rate if the heat source can't modulate flow rate (unit: l/s)
        #[validate(exclusive_minimum = 0.)]
        design_flow_rate: f64,
    },
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "wet_emitter_type", rename_all = "lowercase")]
#[validate(custom = validate_radiator_required_fields)]
pub enum WetEmitter {
    Radiator {
        /// Exponent from characteristic equation of emitters (e.g. derived from BS EN 442 tests)
        #[serde(rename = "n")]
        #[validate(exclusive_minimum = 0.)]
        exponent: f64,

        /// Convective fraction for heating
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        frac_convective: f64,

        // NB. thermal_mass* field should be logically exclusive but can coexist in upstream schema,
        // so collecting as separate fields and applying validation logic for them
        /// Thermal mass of the radiator (unit: kWh/K)
        #[validate(exclusive_minimum = 0.)]
        thermal_mass: Option<f64>,

        /// Thermal mass per meter length of the radiator (unit: kWh/K/m)
        #[validate(exclusive_minimum = 0.)]
        thermal_mass_per_m: Option<f64>,

        #[serde(rename = "c")]
        #[validate(exclusive_minimum = 0.)]
        constant: Option<f64>,

        #[serde(rename = "c_per_m")]
        #[validate(exclusive_minimum = 0.)]
        constant_per_m: Option<f64>,

        /// The length of the emitter (unit: m)
        #[validate(exclusive_minimum = 0.)]
        length: Option<f64>,
    },
    Ufh {
        /// Equivalent thermal mass per m² of floor area for under-floor heating systems (unit: kJ/m²K)
        #[validate(exclusive_minimum = 0.)]
        equivalent_specific_thermal_mass: f64,

        /// Heat output per m² of floor area for under-floor heating systems (unit: W/m²K)
        #[validate(exclusive_minimum = 0.)]
        system_performance_factor: f64,

        /// (unit: m²)
        #[validate(exclusive_minimum = 0.)]
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

fn validate_radiator_required_fields(
    data: &WetEmitter,
) -> Result<(), serde_valid::validation::Error> {
    match data {
        WetEmitter::Radiator {
            exponent,
            frac_convective,
            thermal_mass,
            thermal_mass_per_m,
            constant,
            constant_per_m,
            length,
        } => {
            if constant.is_none() && constant_per_m.is_none() {
                return custom_validation_error("Must provide 'c' or 'c_per_m'".to_string());
            }

            if constant_per_m.is_some() && length.is_none() {
                return custom_validation_error(
                    "Must specify 'length' when 'constant_per_m'  is provided".to_string(),
                );
            }

            if thermal_mass_per_m.is_some() && length.is_none() {
                return custom_validation_error(
                    "Must specify 'length' when 'thermal_mass_per_m'  is provided".to_string(),
                );
            }

            Ok(())
        }
        _ => Ok(()),
    }
}

const fn default_n_units() -> usize {
    1
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[validate(custom = validate_fancoil_test_data)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub struct FancoilTestData {
    #[validate]
    pub(crate) fan_speed_data: Vec<FanSpeedData>,

    /// A list of fan powers for which heat output data is provided (unit: W)
    #[serde(rename = "fan_power_W")]
    #[validate(custom = validate_all_items_at_least_zero)]
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
    #[validate(exclusive_minimum = 0.)]
    pub(crate) temperature_diff: f64,

    /// Heat output for a specific test temperature difference (unit: kW)
    #[validate(custom = validate_all_items_non_negative)]
    pub(crate) power_output: Vec<f64>,
}

fn validate_all_items_non_negative(items: &[f64]) -> Result<(), serde_valid::validation::Error> {
    if items.iter().all(|item| item >= &0.) {
        Ok(())
    } else {
        custom_validation_error("All items must be non-negative".to_string())
    }
}

fn validate_all_items_at_least_zero(items: &[f64]) -> Result<(), serde_valid::validation::Error> {
    if items.iter().all(|item| item > &0.) {
        Ok(())
    } else {
        custom_validation_error("All items must be at least zero".to_string())
    }
}

fn validate_all_items_in_option_non_negative(
    items: &Option<Vec<f64>>,
) -> Result<(), serde_valid::validation::Error> {
    if items.iter().flatten().all(|item| item >= &0.) {
        Ok(())
    } else {
        custom_validation_error("All items must be non-negative".to_string())
    }
}

fn validate_all_items_in_option_at_most_n(
    items: &Option<Vec<f64>>,
    n: f64,
) -> Result<(), serde_valid::validation::Error> {
    if items.iter().flatten().all(|item| item <= &n) {
        Ok(())
    } else {
        custom_validation_error("All items must be non-negative".to_string())
    }
}

fn validate_all_items_in_option<T: Validate>(
    items: &Option<Vec<T>>,
) -> Result<(), serde_valid::validation::Error> {
    if items.iter().flatten().all(|item| item.validate().is_ok()) {
        Ok(())
    } else {
        custom_validation_error("All items must be valid".to_string())
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct SpaceHeatSystemHeatSource {
    pub(crate) name: String,

    /// Upper operating limit for temperature (unit: deg C)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(exclusive_minimum = 0.)]
    pub(crate) temp_flow_limit_upper: Option<f64>,
}

// it is unclear whether this struct should be used - see reference to the struct above
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[allow(dead_code)]
#[serde(deny_unknown_fields)]
pub struct EcoDesignController {
    pub(crate) ecodesign_control_class: EcoDesignControllerClass,

    /// Minimum outdoor temperature (unit: Celsius)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = -273.15)]
    pub(crate) min_outdoor_temp: Option<f64>,

    /// Maximum outdoor temperature (unit: Celsius)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = -273.15)]
    pub(crate) max_outdoor_temp: Option<f64>,

    /// Minimum flow temperature (unit: Celsius)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(exclusive_minimum = 0.)]
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
pub enum ElectricStorageHeaterAirFlowType {
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
    pub area: f64,

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

// From BR 443: The values under "horizontal" apply to heat flow
// directions +/- 30 degrees from horizontal plane.
pub(crate) const PITCH_LIMIT_HORIZ_CEILING: f64 = 60.0;
pub(crate) const PITCH_LIMIT_HORIZ_FLOOR: f64 = 120.0;

#[derive(Clone, Debug, Deserialize, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(tag = "type")]
#[validate(custom = validate_u_value_and_thermal_resistance_floor_construction)]
#[validate(custom = validate_max_window_open_area_for_transparent)]
pub enum BuildingElement {
    #[serde(rename = "BuildingElementOpaque")]
    Opaque {
        #[serde(skip_serializing_if = "Option::is_none")]
        is_unheated_pitched_roof: Option<bool>,

        /// Solar absorption coefficient at the external surface (dimensionless)
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
        solar_absorption_coeff: f64,

        #[serde(flatten)]
        #[validate]
        u_value_input: UValueInput,

        /// Areal heat capacity (unit: J/m².K)
        #[validate(exclusive_minimum = 0.)]
        areal_heat_capacity: f64,

        /// Mass distribution class of the building element, one of: evenly distributed (D); concentrated on external side (E); concentrated on internal side (I); concentrated on internal and external sides (IE); concentrated in middle (M)
        mass_distribution_class: MassDistributionClass,

        /// Tilt angle of the surface from horizontal, between 0 and 180, where 0 means the external surface is facing up, 90 means the external surface is vertical and 180 means the external surface is facing down (unit: ˚)
        #[validate(minimum = 0.)]
        #[validate(maximum = 180.)]
        pitch: f64,

        #[validate]
        orientation360: Option<Orientation360>,

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
        control_window_openable: Option<String>,

        /// Tilt angle of the surface from horizontal, between 0 and 180, where 0 means the external surface is facing up, 90 means the external surface is vertical and 180 means the external surface is facing down (unit: ˚
        #[validate(minimum = 0.)]
        #[validate(maximum = 180.)]
        pitch: f64,

        #[validate]
        orientation360: Orientation360,

        /// Total solar energy transmittance of the transparent part of the window
        #[validate(minimum = 0.)]
        g_value: f64,

        /// The frame area fraction of window, ratio of the projected frame area to the overall projected area of the glazed element of the window
        #[validate(minimum = 0.)]
        #[validate(maximum = 1.)]
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
        #[validate(exclusive_minimum = 0.)]
        area: f64,

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
        #[validate(exclusive_minimum = 0.)]
        area: f64,

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

    /// Party wall element for all party wall types, with specific handling for cavity air movement heat loss.
    ///
    /// For cavity party walls, this accounts for non-negligible heat loss due to air movement
    /// within the cavity, as identified by research. The thermal resistance of the
    /// unconditioned space (cavity) is derived from the party wall cavity type and lining type.
    /// Five cavity types are supported:
    ///   solid, unfilled_unsealed, unfilled_sealed, filled_sealed, and defined_resistance.
    /// Two lining types are supported:
    ///   wet_plaster, and dry_lined
    #[serde(rename = "BuildingElementPartyWall")]
    PartyWall {
        /// Tilt angle of the surface from horizontal, between 60 and 120 degrees (wall range), where 90 means vertical (unit: °)
        #[validate(minimum = PITCH_LIMIT_HORIZ_CEILING)]
        #[validate(maximum = PITCH_LIMIT_HORIZ_FLOOR)]
        pitch: f64,

        /// Area of this building element (unit: m²)
        #[validate(exclusive_minimum = 0.)]
        area: f64,

        /// Areal heat capacity (unit: J/m2.K)
        #[validate(exclusive_minimum = 0.)]
        areal_heat_capacity: f64,

        mass_distribution_class: MassDistributionClass,

        #[serde(flatten)]
        #[validate]
        party_wall_cavity_data: PartyWallCavityData,

        #[serde(flatten)]
        #[validate]
        u_value_input: UValueInput,
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
pub enum UValueInput {
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
pub struct BuildingElementAreaOrHeightWidthInput {
    /// Area of the building element (m²)
    #[validate(exclusive_minimum = 0.)]
    area: Option<f64>,

    #[serde(flatten)]
    #[validate]
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

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "party_wall_cavity_type", rename_all = "snake_case")]
pub enum PartyWallCavityData {
    /// Solid wall or structurally insulated panel
    Solid,
    /// Unfilled cavity with no effective edge sealing
    UnfilledUnsealed {
        party_wall_lining_type: PartyWallLiningType,
    },
    /// Unfilled cavity with effective sealing
    UnfilledSealed {
        party_wall_lining_type: PartyWallLiningType,
    },
    /// Fully filled cavity with effective sealing
    FilledSealed {
        party_wall_lining_type: PartyWallLiningType,
    },
    /// Fully filled cavity with no effective edge sealing
    FilledUnsealed,
    /// User-defined thermal resistance
    DefinedResistance {
        #[validate(exclusive_minimum = 0.)]
        thermal_resistance_cavity: Option<f64>,
    },
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

fn validate_max_window_open_area_for_transparent(
    element: &BuildingElement,
) -> Result<(), serde_valid::validation::Error> {
    let (area, max_window_open_area) = if let BuildingElement::Transparent {
        area_input,
        max_window_open_area,
        ..
    } = element
    {
        (area_input.area(), *max_window_open_area)
    } else {
        return Ok(());
    };

    if max_window_open_area > area {
        return custom_validation_error(
            "max_window_open_area must be less than or equal to the area".to_string(),
        );
    }

    Ok(())
}

impl BuildingElement {
    pub fn pitch(&self) -> f64 {
        *match self {
            BuildingElement::Opaque { pitch, .. } => pitch,
            BuildingElement::Transparent { pitch, .. } => pitch,
            BuildingElement::Ground { pitch, .. } => pitch,
            BuildingElement::AdjacentConditionedSpace { pitch, .. } => pitch,
            BuildingElement::AdjacentUnconditionedSpace { pitch, .. } => pitch,
            BuildingElement::PartyWall { pitch, .. } => pitch,
        }
    }

    pub fn u_value(&self) -> Option<f64> {
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
            BuildingElement::PartyWall { u_value_input, .. } => {
                if let UValueInput::UValue { u_value } = u_value_input {
                    Some(*u_value)
                } else {
                    None
                }
            }
        }
    }

    pub fn orientation(&self) -> Option<Orientation360> {
        match self {
            BuildingElement::Opaque {
                orientation360: orientation,
                ..
            } => *orientation,
            BuildingElement::Transparent {
                orientation360: orientation,
                ..
            } => Some(*orientation),
            _ => None,
        }
    }

    #[cfg(test)]
    pub fn remove_window_openable_control(&mut self) {
        if let BuildingElement::Transparent {
            control_window_openable: window_openable_control,
            ..
        } = self
        {
            *window_openable_control = None
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

/// Types of party wall cavity configurations
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "snake_case")]
pub(crate) enum PartyWallCavityType {
    /// Solid wall or structurally insulated panel
    Solid,
    /// Unfilled cavity with no effective edge sealing
    UnfilledUnsealed,
    /// Unfilled cavity with effective sealing
    UnfilledSealed,
    /// Fully filled cavity with effective sealing
    FilledSealed,
    /// Fully filled cavity with no effective edge sealing
    FilledUnsealed,
    /// User-defined thermal resistance
    DefinedResistance,
}

impl From<PartyWallCavityData> for PartyWallCavityType {
    fn from(value: PartyWallCavityData) -> Self {
        match value {
            PartyWallCavityData::Solid => Self::Solid,
            PartyWallCavityData::UnfilledUnsealed { .. } => Self::UnfilledUnsealed,
            PartyWallCavityData::UnfilledSealed { .. } => Self::UnfilledSealed,
            PartyWallCavityData::FilledSealed { .. } => Self::FilledSealed,
            PartyWallCavityData::FilledUnsealed => Self::FilledUnsealed,
            PartyWallCavityData::DefinedResistance { .. } => Self::DefinedResistance,
        }
    }
}

impl PartyWallCavityData {
    pub(crate) fn party_wall_lining_type(&self) -> Option<PartyWallLiningType> {
        match self {
            Self::UnfilledUnsealed {
                party_wall_lining_type,
            }
            | Self::UnfilledSealed {
                party_wall_lining_type,
            }
            | Self::FilledSealed {
                party_wall_lining_type,
            } => Some(*party_wall_lining_type),
            _ => None,
        }
    }

    pub(crate) fn thermal_resistance_cavity(&self) -> Option<f64> {
        match self {
            Self::DefinedResistance {
                thermal_resistance_cavity,
            } => *thermal_resistance_cavity,
            _ => None,
        }
    }
}

/// Types of party wall lining
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "snake_case")]
pub enum PartyWallLiningType {
    WetPlaster,
    DryLined,
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

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct WindowPart {
    /// (unit: m)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) mid_height_air_flow_path: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct WindowTreatment {
    #[serde(rename = "type")]
    pub(crate) treatment_type: WindowTreatmentType,

    pub(crate) controls: WindowTreatmentControl,

    /// Additional thermal resistance provided by a window treatment (unit: m²K/W)
    #[validate(exclusive_minimum = 0.)]
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
    pub(crate) control_closing_irrad: Option<String>,

    /// Irradiation level below which a window treatment is assumed to be open (unit: W/m²). References a key in $.Control.
    #[serde(
        rename = "Control_opening_irrad",
        skip_serializing_if = "Option::is_none"
    )]
    pub(crate) control_opening_irrad: Option<String>,

    /// Reference to a time control object containing a schedule of booleans describing when a window treatment is open. References a key in $.Control.
    #[serde(rename = "Control_open", skip_serializing_if = "Option::is_none")]
    pub(crate) control_open: Option<String>,

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
pub enum FloorData {
    #[serde(rename = "Slab_no_edge_insulation")]
    SlabNoEdgeInsulation,

    #[serde(rename = "Slab_edge_insulation")]
    SlabEdgeInsulation {
        #[serde(default)]
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
        #[validate(exclusive_minimum = 0.)]
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

        // Edge insulation not typically used for suspended floors
        #[serde(default)]
        #[validate]
        edge_insulation: Vec<EdgeInsulation>,
    },

    #[serde(rename = "Heated_basement")]
    HeatedBasement {
        /// Depth of basement floor below ground level (unit: m)
        #[validate(exclusive_minimum = 0.)]
        depth_basement_floor: f64,

        /// Thermal resistance of walls of the basement (unit: m².K/W)
        #[serde(rename = "thermal_resist_walls_base")]
        #[validate(exclusive_minimum = 0.)]
        thermal_resistance_of_basement_walls: f64,

        // Optional - edge insulation can be used with basements
        #[serde(default)]
        #[validate]
        edge_insulation: Vec<EdgeInsulation>,
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
        #[validate(exclusive_minimum = 0.)]
        depth_basement_floor: f64,

        /// Height of the basement walls above ground level (unit: m)
        #[validate(exclusive_minimum = 0.)]
        height_basement_walls: f64,

        /// Thermal resistance of walls of the basement (unit: m².K/W)
        #[serde(rename = "thermal_resist_walls_base")]
        #[validate(exclusive_minimum = 0.)]
        thermal_resistance_of_basement_walls: f64,

        // Optional - edge insulation can be used with basements
        #[serde(default)]
        #[validate]
        edge_insulation: Vec<EdgeInsulation>,
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
        #[validate(exclusive_minimum = 0.)]
        width: f64,

        /// Thermal resistance of floor edge insulation (unit: m²K/W)
        #[validate(exclusive_minimum = 0.)]
        edge_thermal_resistance: f64,
    },
    #[serde(rename = "vertical")]
    Vertical {
        /// (unit: m)
        #[validate(exclusive_minimum = 0.)]
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
#[serde(tag = "type", deny_unknown_fields)]
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
#[serde(tag = "type", deny_unknown_fields)]
pub(crate) enum SpaceCoolSystemDetails {
    AirConditioning {
        /// Maximum cooling capacity of the system (unit: kW)
        #[validate(exclusive_minimum = 0.)]
        cooling_capacity: f64,

        /// Efficiency of the air conditioning system. SEER (Seasonal energy efficiency ratio)
        #[validate(exclusive_minimum = 0.)]
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

pub type HeatSourceWet = IndexMap<std::string::String, HeatSourceWetDetails>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[allow(clippy::large_enum_variant)]
#[serde(tag = "type")]
#[validate(custom = validate_backup_configuration)]
pub enum HeatSourceWetDetails {
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
        #[validate(minimum = 0.)]
        min_temp_diff_flow_return_for_hp_to_operate: f64,

        /// Whether the heat pump uses modulating control
        modulating_control: bool,

        /// Power consumption of crankcase heater (unit: kW)
        #[validate(minimum = 0.)]
        power_crankcase_heater: f64,

        /// Power consumption of heating circuit pump (unit: kW)
        #[validate(minimum = 0.)]
        #[serde(default = "default_power_heating_circ_pump")]
        power_heating_circ_pump: f64,

        /// Power consumption of warm air fan (unit: kW)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(minimum = 0.)]
        power_heating_warm_air_fan: Option<f64>,

        /// Maximum backup power (unit: kW)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(exclusive_minimum = 0.)]
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
        #[validate(exclusive_minimum = 0.)]
        temp_distribution_heat_network: Option<f64>,

        /// Lower temperature limit for heat pump operation (unit: ˚C)
        #[validate(minimum = -273.15)]
        temp_lower_operating_limit: f64,

        /// Maximum return feed temperature (unit: ˚C)
        #[serde(skip_serializing_if = "Option::is_none")]
        #[validate(exclusive_minimum = 0.)]
        temp_return_feed_max: Option<f64>,

        #[serde(rename = "test_data_EN14825")]
        #[validate]
        test_data_en14825: Vec<HeatPumpTestDatum>,

        /// Time constant for on/off operation (unit: hours)
        #[validate(exclusive_minimum = 0.)]
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
        #[validate(exclusive_minimum = 0.)]
        #[validate(maximum = 1.)]
        efficiency_full_load: f64,

        /// Boiler efficiency at part load (dimensionless, 0-1)
        #[validate(exclusive_minimum = 0.)]
        #[validate(maximum = 1.12)]
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
        #[validate(exclusive_minimum = 0.)]
        electricity_part_load: f64,

        /// Electrical power consumption at full load (unit: kW)
        #[validate(exclusive_minimum = 0.)]
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
        #[validate(exclusive_minimum = 0.)]
        hiu_daily_loss: f64,

        /// Heat losses from building-level distribution pipework (unit: W)
        #[validate(minimum = 0.)]
        building_level_distribution_losses: f64,

        /// Power consumption of heating circulation pump (unit: kW)
        #[validate(minimum = 0.)]
        power_circ_pump: Option<f64>,

        /// Power consumption of auxiliary electrical usage (unit: kW)
        #[validate(exclusive_minimum = 0.)]
        power_aux: Option<f64>,
    },
}

const fn default_power_heating_circ_pump() -> f64 {
    0.
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

#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct HeatPumpBufferTank {
    /// Standing heat loss (unit: kWh/day)
    #[validate(exclusive_minimum = 0.)]
    pub daily_losses: f64,

    /// Volume of the buffer tank (unit: litre)
    #[validate(exclusive_minimum = 0.)]
    pub volume: f64,

    /// Flow rate of the buffer tank - emitters loop (unit: l/min)
    #[validate(exclusive_minimum = 0.)]
    pub pump_fixed_flow_rate: f64,

    /// Pump power of the buffer tank - emitters loop (unit: W)
    #[validate(exclusive_minimum = 0.)]
    pub pump_power_at_flow_rate: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Validate, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct HeatPumpTestDatum {
    /// Air flow rate through the heat pump for the test condition (unit: m³/h)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(exclusive_minimum = 0.)]
    pub(crate) air_flow_rate: Option<f64>,

    /// Heat output capacity at this test condition (unit: kW)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) capacity: f64,

    /// Coefficient of performance at this test condition (dimensionless)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) cop: f64,

    /// Design flow temperature for the heating system (unit: Celsius)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) design_flow_temp: f64,

    /// Ratio of external air to recirculated air for exhaust air heat pumps (dimensionless)
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = 0.)]
    #[validate(maximum = 1.)]
    pub(crate) eahp_mixed_ext_air_ratio: Option<f64>,

    /// Heat pump outlet temperature for the test condition (unit: Celsius)
    #[validate(exclusive_minimum = 0.)]
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
    #[validate(exclusive_minimum = 0.)]
    #[validate(maximum = 1.)]
    efficiency_full_load: f64,

    /// Boiler efficiency at part load (dimensionless, 0-1)
    #[validate(exclusive_minimum = 0.)]
    #[validate(maximum = 1.12)]
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
    #[validate(exclusive_minimum = 0.)]
    electricity_part_load: f64,

    /// Electrical power consumption at full load (unit: kW)
    #[validate(exclusive_minimum = 0.)]
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
    /// Cost data for the fuel used by the hybrid system's boiler (can be any units, typically p/kWh,
    /// as long as they are consistent across input fields
    pub(crate) cost_schedule_boiler: NumericSchedule,

    /// Cost data for the fuel used by the hybrid system's
    /// heat pump (can be any units, typically p/kWh, as long as they are consistent
    /// across input fields)"
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
pub enum HeatBattery {
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
        #[validate(exclusive_minimum = 0.)]
        electricity_standby: f64,

        /// Rated charging power (unit: kW)
        #[validate(exclusive_minimum = 0.)]
        rated_charge_power: f64,

        /// Maximum rated heat losses (unit: kW)
        #[validate(exclusive_minimum = 0.)]
        max_rated_losses: f64,

        /// Number of heat battery units
        #[validate(minimum = 1)]
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

        /// Inlet diameter of capillary tubes (unit: mm)
        #[validate(exclusive_minimum = 0.)]
        inlet_diameter_mm: f64,

        /// Heat battery parameter A (dimensionless)
        #[serde(rename = "A")]
        a: f64,

        /// Heat battery parameter B (dimensionless)
        #[serde(rename = "B")]
        b: f64,

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

fn validate_backup_configuration(
    data: &HeatSourceWetDetails,
) -> Result<(), serde_valid::validation::Error> {
    let (boiler, power_max_backup, backup_ctrl_type, time_delay_backup) =
        if let HeatSourceWetDetails::HeatPump {
            boiler,
            power_max_backup,
            backup_control_type,
            time_delay_backup,
            ..
        } = data
        {
            (
                boiler,
                power_max_backup,
                backup_control_type,
                time_delay_backup,
            )
        } else {
            return Ok(());
        };

    if boiler.is_some() && power_max_backup.is_some() {
        return custom_validation_error(
            "power_max_backup and boiler can not both be set.".to_string(),
        );
    }
    if matches!(backup_ctrl_type, HeatPumpBackupControlType::None) {
        if boiler.is_some() {
            return custom_validation_error(
                "boiler can not be set if backup_ctrl_type is 'None'.".to_string(),
            );
        }
        if power_max_backup.is_some() {
            return custom_validation_error(
                "power_max_backup can not be set if backup_ctrl_type is 'None'.".to_string(),
            );
        }
    } else {
        if time_delay_backup.is_none() {
            return custom_validation_error(
                "time_delay_backup is required if backup_ctrl_type is set.".to_string(),
            );
        }
        if boiler.is_none() && power_max_backup.is_none() {
            return custom_validation_error(
                "Either power_max_backup or boiler is required if backup_ctrl_type is set."
                    .to_string(),
            );
        }
    }

    Ok(())
}

pub type WasteWaterHeatRecovery = IndexMap<std::string::String, WasteWaterHeatRecoveryDetails>;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
#[validate(custom = validate_at_least_one_efficiency_set)]
#[validate(custom = validate_flow_rates_and_efficiencies_length)]
pub struct WasteWaterHeatRecoveryDetails {
    #[serde(rename = "type")]
    #[cfg_attr(feature = "arbitrary", arbitrary(value = MustBe!("WWHRS_Instantaneous")))]
    _type: MustBe!("WWHRS_Instantaneous"),

    #[serde(rename = "ColdWaterSource")]
    pub(crate) cold_water_source: String,

    /// Test flow rates in litres per minute (e.g., [5., 7., 9., 11., 13.])
    #[validate(custom = validate_all_items_at_least_zero)]
    pub(crate) flow_rates: Vec<f64>,

    /// Measured efficiencies for System A at the test flow rates
    #[validate(custom = validate_all_items_in_option_valid_system_efficiency)]
    pub(crate) system_a_efficiencies: Option<Vec<f64>>,

    /// Utilisation factor for System A
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(exclusive_minimum = 0.)]
    #[validate(maximum = 1.)]
    pub(crate) system_a_utilisation_factor: Option<f64>,

    /// Measured efficiencies for System B (optional, uses system_b_efficiency_factor if not provided)
    #[validate(custom = validate_all_items_in_option_valid_system_efficiency)]
    pub(crate) system_b_efficiencies: Option<Vec<f64>>,

    /// Utilisation factor for System B. Required when using either system_b_efficiencies (pre-corrected data) or when converting system_a_efficiencies to System B (used with system_b_efficiency_factor).
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(exclusive_minimum = 0.)]
    #[validate(maximum = 1.)]
    pub(crate) system_b_utilisation_factor: Option<f64>,

    /// Measured efficiencies for System C (optional, uses system_c_efficiency_factor if not provided)
    #[validate(custom = validate_all_items_in_option_valid_system_efficiency)]
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

fn validate_at_least_one_efficiency_set(
    data: &WasteWaterHeatRecoveryDetails,
) -> Result<(), serde_valid::validation::Error> {
    if let WasteWaterHeatRecoveryDetails {
        system_a_efficiencies: None,
        system_b_efficiencies: None,
        system_c_efficiencies: None,
        ..
    } = data
    {
        return custom_validation_error("At least one efficiency dataset must be provided: system_a_efficiencies, system_b_efficiencies, or system_c_efficiencies".to_string());
    }

    Ok(())
}

fn validate_flow_rates_and_efficiencies_length(
    data: &WasteWaterHeatRecoveryDetails,
) -> Result<(), serde_valid::validation::Error> {
    // Validate system A (always required)
    if data.flow_rates.len() != data.system_a_efficiencies.iter().flatten().count() {
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

fn validate_all_items_in_option_valid_system_efficiency(
    items: &Option<Vec<f64>>,
) -> Result<(), serde_valid::validation::Error> {
    if items
        .iter()
        .flatten()
        .all(|&item| item > 0. && item <= 100.)
    {
        Ok(())
    } else {
        custom_validation_error(
            "All items must be a number representing a valid efficiency for a WWHRS system"
                .to_string(),
        )
    }
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
    #[validate(exclusive_minimum = 0.)]
    pub(crate) peak_power: f64,

    pub(crate) ventilation_strategy: PhotovoltaicVentilationStrategy,

    /// The tilt angle (inclination) of the PV panel from horizontal, measured upwards facing, 0 to 90 (unit: ˚)
    #[validate(minimum = 0.)]
    #[validate(maximum = 90.)]
    pub(crate) pitch: f64,

    #[validate]
    pub(crate) orientation360: Orientation360,

    /// The distance between the ground and the lowest edge of the PV array (unit: m)
    #[validate(minimum = 0.)]
    pub(crate) base_height: f64,

    /// Height of the PV array (unit: m)
    pub(crate) height: f64,

    /// Width of the PV panel (unit: m)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) width: f64,

    #[validate]
    pub(crate) shading: Vec<WindowShadingObject>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) struct PhotovoltaicSystemWithPanels {
    #[serde(rename = "type")]
    #[cfg_attr(feature = "arbitrary", arbitrary(value = MustBe!("PhotovoltaicSystem")))]
    _type: MustBe!("PhotovoltaicSystem"),

    #[serde(rename = "EnergySupply")]
    pub(crate) energy_supply: String,

    /// Whether the inverter is considered inside the building
    pub(crate) inverter_is_inside: bool,

    /// Peak power; represents the peak electrical AC power output from the inverter (unit: kW)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) inverter_peak_power_ac: f64,

    /// Peak power; represents the peak electrical DC power input to the inverter (unit: kW)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) inverter_peak_power_dc: f64,

    pub(crate) inverter_type: InverterType,

    #[validate(min_items = 1)]
    pub(crate) panels: Vec<PhotovoltaicPanel>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct PhotovoltaicSystem {
    #[serde(rename = "type")]
    #[cfg_attr(feature = "arbitrary", arbitrary(value = MustBe!("PhotovoltaicSystem")))]
    _type: MustBe!("PhotovoltaicSystem"),

    /// Peak power; represents the electrical power of a photovoltaic system with a given area for a solar irradiance of 1 kW/m² on this surface (at 25 degrees) (unit: kW)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) peak_power: f64,

    pub(crate) ventilation_strategy: PhotovoltaicVentilationStrategy,

    /// The tilt angle (inclination) of the PV panel from horizontal, measured upwards facing, 0 to 90 (unit: ˚)
    #[validate(minimum = 0.)]
    #[validate(maximum = 90.)]
    pub(crate) pitch: f64,

    #[validate]
    pub(crate) orientation360: Orientation360,

    /// The distance between the ground and the lowest edge of the PV array (unit: m)
    #[validate(minimum = 0.)]
    pub(crate) base_height: f64,

    /// Height of the PV array (unit: m)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) height: f64,

    /// Width of the PV panel (unit: m)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) width: f64,

    #[serde(rename = "EnergySupply")]
    pub(crate) energy_supply: String,

    pub(crate) shading: Vec<WindowShadingObject>,

    /// Peak power; represents the peak electrical AC power output from the inverter (unit: kW)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) inverter_peak_power_ac: f64,

    /// Peak power; represents the peak electrical DC power input to the inverter (unit: kW)
    #[validate(exclusive_minimum = 0.)]
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
#[serde(deny_unknown_fields)]
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
    pub(crate) ventilation_zone_base_height: f64,

    /// Initial vent position, 0 = vents closed and 1 = vents fully open
    #[serde(skip_serializing_if = "Option::is_none")]
    #[validate(minimum = 0.)]
    #[validate(maximum = 1.)]
    pub(crate) vent_opening_ratio_init: Option<f64>,
}

impl InfiltrationVentilation {
    pub fn mechanical_ventilation(&self) -> &IndexMap<std::string::String, MechanicalVentilation> {
        &self.mechanical_ventilation
    }
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
    #[validate(exclusive_minimum = 0.)]
    pub(crate) pressure_difference_ref: f64,

    #[validate]
    pub(crate) orientation360: Orientation360,

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
    #[validate(exclusive_minimum = 0.)]
    pub(crate) ventilation_zone_height: f64,

    /// Reference pressure difference (unit: Pa)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) test_pressure: f64,

    /// Flow rate through (unit: m³/h.m²)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) test_result: f64,

    /// Reference area of the envelope airtightness index
    #[validate(exclusive_minimum = 0.)]
    pub(crate) env_area: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
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
    #[validate(maximum = 1.)]
    pub(crate) mvhr_efficiency: Option<f64>,

    /// Location of the MVHR unit (inside or outside the thermal envelope)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) mvhr_location: Option<MVHRLocation>,

    #[serde(rename = "Control", skip_serializing_if = "Option::is_none")]
    pub(crate) control: Option<String>,

    /// Specific fan power, inclusive of any in use factors (unit: W/l/s)
    #[serde(rename = "SFP")]
    #[validate(exclusive_minimum = 0.)]
    pub(crate) sfp: f64,

    /// Adjustment factor to be applied to SFP to account for e.g. type of ducting. Typical range 1 - 2.5
    #[serde(default = "default_sfp_in_use_factor", rename = "SFP_in_use_factor")]
    #[validate(minimum = 1.)]
    pub(crate) sfp_in_use_factor: f64,

    #[serde(rename = "EnergySupply")]
    pub(crate) energy_supply: String,

    /// (unit: m³/hour)
    #[validate(exclusive_minimum = 0.)]
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

    pub(crate) fn position_exhaust(&self) -> (Orientation360, f64, f64) {
        match self {
            Self::Mvhr {
                position_exhaust, ..
            }
            | Self::IntermittentMev {
                position_exhaust, ..
            }
            | Self::CentralisedContinuousMev {
                position_exhaust, ..
            }
            | Self::DecentralisedContinuousMev {
                position_exhaust, ..
            }
            | Self::PositiveInputVentilation {
                position_exhaust, ..
            } => {
                let MechanicalVentilationPosition {
                    orientation360: orientation,
                    pitch,
                    mid_height_air_flow_path,
                } = *position_exhaust;

                (orientation, pitch, mid_height_air_flow_path)
            }
        }
    }

    pub(crate) fn position_intake(&self) -> Option<(Orientation360, f64, f64)> {
        match self {
            Self::Mvhr {
                position_intake, ..
            } => Some((
                position_intake.orientation360,
                position_intake.pitch,
                position_intake.mid_height_air_flow_path,
            )),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub struct MechanicalVentilationPosition {
    #[validate]
    pub(crate) orientation360: Orientation360,

    /// Tilt angle of the surface from horizontal, between 0 and 180, where 0 means the external surface is facing up, 90 means the external surface is vertical and 180 means the external surface is facing down (unit: ˚
    #[validate(minimum = 0.)]
    #[validate(maximum = 180.)]
    pub(crate) pitch: f64,

    /// Mid height of air flow path relative to ventilation zone (unit: m)
    #[validate(exclusive_minimum = 0.)]
    pub(crate) mid_height_air_flow_path: f64,
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
    // commented out from Python
    // #[serde(rename = "LOAD")]
    // Load,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum SupplyAirTemperatureControlType {
    // commented out from Python
    // #[serde(rename = "CONST")]
    // Constant,
    #[serde(rename = "NO_CTRL")]
    NoControl,
    // commented out from Python
    // #[serde(rename = "LOAD_COM")]
    // LoadCom,
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
    #[validate(exclusive_minimum = 0.)]
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
#[serde(deny_unknown_fields)]
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

#[expect(unused)]
static CORE_SCHEMA_VALIDATOR: LazyLock<Validator> = LazyLock::new(|| {
    let schema = serde_json::from_str(include_str!("../schemas/core-input.schema.json")).unwrap();
    jsonschema::validator_for(&schema).unwrap()
});

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::space_heat_demand::ventilation::MechVentType;
    use rstest::*;
    use serde::de::DeserializeOwned;
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
                    // exclude FHS input files from unit test runs as they are very long
                    // and rest of files give a good spread
                    && !e.file_name().to_str().unwrap().contains("FHS")
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
            let parsed: Result<Input, _> =
                serde_json::from_reader(BufReader::new(File::open(entry.path()).unwrap()));
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
        let file = File::open("./examples/input/core/short/demo.json").unwrap();
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

    #[rstest]
    fn test_validate_time_series(baseline_demo_file_json: JsonValue) {
        let mut modified_input = baseline_demo_file_json.clone();
        modified_input["ColdWaterSource"]["mains water"]["temperatures"] = json!([0.0]);

        let input = serde_json::from_value::<Input>(modified_input).unwrap();
        if let Err(e) = input.validate() {
            assert!(
                e.to_string().contains("ColdWaterSource.temperatures does not contain enough values to cover the simulation."),
            );
        } else {
            panic!("Expected validation error");
        }
        assert!(input.validate().is_err());
    }

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
            case::at_least_zero_time_series_step(json!({"time_series_step": 0})),
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
            case::demand_greater_than_zero(json!({"demand_W": 0})),
            case::duration_greater_than_zero(json!({"duration": 0})),
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
            case::capacity_greater_than_zero(json!({"capacity": 0})),
            case::maximum_charge_rate_greater_than_zero(json!({"maximum_charge_rate_one_way_trip": 0})
            ),
            case::maximum_discharge_rate_greater_than_zero(json!({"maximum_discharge_rate_one_way_trip": 0})
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
            case::temperature_diff_greater_than_zero(json!({"temperature_diff": 0})),
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
            case::efficiency_full_load_greater_than_zero(json!({"efficiency_full_load": 0})
            ),
            case::efficiency_full_load_at_most_one(json!({"efficiency_full_load": 2})
            ),
            case::efficiency_part_load_greater_than_zero(json!({"efficiency_part_load": 0})
            ),
            case::efficiency_part_load_at_most_one_point_one_two(json!({"efficiency_part_load": 2})
            ),
            case::modulation_load_at_least_zero(json!({"modulation_load": -1})),
            case::modulation_load_at_most_one(json!({"modulation_load": 2})),
            case::electricity_circ_pump_at_least_zero(json!({"electricity_circ_pump": -1})
            ),
            case::electricity_full_load_greater_than_zero(json!({"electricity_full_load": 0})
            ),
            case::electricity_part_load_greater_than_zero(json!({"electricity_part_load": 0})
            ),
            case::electricity_standby_at_least_zero(json!({"electricity_standby": -0.1})
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
                    inlet_diameter_mm: 6.5,
                    a: 3.532,
                    b: 4.415,
                    flow_rate_l_per_min: 10.,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::inlet_diameter_mm_greater_than_zero(json!({"inlet_diameter_mm": 0})),
                case::flow_rate_l_per_min_greater_than_zero(json!({"flow_rate_l_per_min": 0})),
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
                case::electricity_standby_greater_than_zero(json!({"electricity_standby": 0})
                ),
                case::max_rated_losses_greater_than_zero(json!({"max_rated_losses": 0})
                ),
                case::number_of_units_greater_than_zero(json!({"number_of_units": 0})),
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
                case::hiu_daily_loss_greater_than_zero(json!({"HIU_daily_loss": 0})),
                case::building_level_distribution_losses_at_least_zero(json!({"building_level_distribution_losses": -1})
                ),
                case::power_circ_pump_at_least_zero(json!({"power_circ_pump": -1})),
                case::power_aux_greater_than_zero(json!({"power_aux": 0})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<HeatSourceWetDetails>(valid_example, inputs);
            }
        }

        mod heat_pump {
            use super::*;

            #[fixture]
            fn valid_heat_pump() -> JsonValue {
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
                    power_heating_circ_pump: 0.015,
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
                    time_delay_backup: Some(0.5),
                    var_flow_temp_ctrl_during_test: true,
                })
                .unwrap()
            }

            #[fixture]
            fn valid_boiler() -> HeatPumpBoiler {
                HeatPumpBoiler {
                    energy_supply: "mains_gas".into(),
                    energy_supply_aux: "mains_elec".into(),
                    rated_power: 6.0,
                    efficiency_full_load: 0.9,
                    efficiency_part_load: 0.7,
                    boiler_location: HeatSourceLocation::Internal,
                    modulation_load: 0.0,
                    electricity_circ_pump: 0.0,
                    electricity_part_load: 0.1,
                    electricity_full_load: 0.2,
                    electricity_standby: 0.01,
                    cost_schedule_hybrid: None,
                }
            }

            #[rstest(inputs,
                case::eahp_mixed_max_temp_at_least_absolute_zero(json!({"eahp_mixed_max_temp": -274})
                ),
                case::eahp_mixed_min_temp_at_least_absolute_zero(json!({"eahp_mixed_min_temp": -274})
                ),
                case::temp_distribution_heat_network_greater_than_zero(json!({"temp_distribution_heat_network": 0})
                ),
                case::temp_lower_operating_limit_at_least_absolute_zero(json!({"temp_lower_operating_limit": -274})
                ),
                case::temp_return_feed_max_greater_than_zero(json!({"temp_return_feed_max": 0})
                ),
                case::min_modulation_rate_20_at_least_zero(json!({"min_modulation_rate_20": -1})),
                case::min_modulation_rate_20_at_most_one(json!({"min_modulation_rate_20": 2})),
                case::min_modulation_rate_35_at_least_zero(json!({"min_modulation_rate_35": -1})),
                case::min_modulation_rate_35_at_most_one(json!({"min_modulation_rate_35": 2})),
                case::min_modulation_rate_55_at_least_zero(json!({"min_modulation_rate_55": -1})),
                case::min_modulation_rate_55_at_most_one(json!({"min_modulation_rate_55": 2})),
                case::power_crankcase_heater_at_least_zero(json!({"power_crankcase_heater": -1})),
                case::power_heating_circ_pump_at_least_zero(json!({"power_heating_circ_pump": -1})),
                case::power_heating_warm_air_fan_at_least_zero(json!({"power_heating_warm_air_fan": -1})
                ),
                case::power_max_backup_greater_than_zero(json!({"power_max_backup": 0})),
                case::power_off_at_least_zero(json!({"power_off": -1})),
                case::power_source_circ_pump_at_least_zero(json!({"power_source_circ_pump": -1})),
                case::power_standby_at_least_zero(json!({"power_standby": -0.1})),
                case::time_constant_onoff_operation_greater_than_zero(json!({"time_constant_onoff_operation": 0})
                ),
                case::time_delay_backup_at_least_zero(json!({"time_delay_backup": -1})),
            )]
            fn test_validate_range_constraints(valid_heat_pump: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<HeatSourceWetDetails>(valid_heat_pump, inputs);
            }

            #[rstest]
            #[case(json!({
                "boiler": null,
                "backup_ctrl_type": "TopUp",
                "power_max_backup": 3.0,
                "time_delay_backup": 0.5,
            }), true, None)]
            #[case(json!({
                "boiler": null,
                "backup_ctrl_type": "Substitute",
                "power_max_backup": 3.0,
                "time_delay_backup": 0.5,
            }), true, None)]
            #[case(json!({
                "boiler": true,
                "backup_ctrl_type": "TopUp",
                "power_max_backup": null,
                "time_delay_backup": 0.5,
            }), true, None)]
            #[case(json!({
                "boiler": true,
                "backup_ctrl_type": "Substitute",
                "power_max_backup": null,
                "time_delay_backup": 0.5,
            }), true, None)]
            #[case(json!({
                "boiler": null,
                "backup_ctrl_type": "None",
                "power_max_backup": 3.0,
                "time_delay_backup": 0.5,
            }), false, "power_max_backup can not be set if backup_ctrl_type is 'None'.".into())]
            #[case(json!({
                "boiler": true,
                "backup_ctrl_type": "None",
                "power_max_backup": null,
                "time_delay_backup": 0.5,
            }), false, "boiler can not be set if backup_ctrl_type is 'None'.".into())]
            #[case(json!({
                "boiler": null,
                "backup_ctrl_type": "TopUp",
                "power_max_backup": 3.0,
                "time_delay_backup": null,
            }), false, "time_delay_backup is required if backup_ctrl_type is set.".into())]
            #[case(json!({
                "boiler": true,
                "backup_ctrl_type": "None",
                "power_max_backup": 3.0,
                "time_delay_backup": 0.5,
            }), false, "power_max_backup and boiler can not both be set.".into())]
            #[case(json!({
                "boiler": null,
                "backup_ctrl_type": "TopUp",
                "power_max_backup": null,
                "time_delay_backup": 0.5,
            }), false, "Either power_max_backup or boiler is required if backup_ctrl_type is set.".into())]
            fn test_validate_backup_configuration(
                #[case] mut inputs: JsonValue,
                #[case] valid: bool,
                #[case] exception_match: Option<&str>,
                mut valid_heat_pump: JsonValue,
                valid_boiler: HeatPumpBoiler,
            ) {
                if let JsonValue::Bool(true) = inputs["boiler"] {
                    inputs["boiler"] = serde_json::to_value(valid_boiler).unwrap();
                }
                valid_heat_pump
                    .as_object_mut()
                    .unwrap()
                    .extend(inputs.as_object().unwrap().to_owned());
                if valid {
                    let heat_pump = serde_json::from_value::<HeatSourceWetDetails>(valid_heat_pump);
                    assert!(heat_pump.is_ok());
                    assert!(heat_pump.unwrap().validate().is_ok());
                } else {
                    assert!(
                        serde_json::from_value::<HeatSourceWetDetails>(valid_heat_pump).is_ok_and(
                            |heat_pump| {
                                if let Err(e) = heat_pump.validate() {
                                    e.to_string().contains(exception_match.unwrap())
                                } else {
                                    panic!("Validation of heat pump was expected to fail.")
                                }
                            }
                        )
                    );
                }
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
                    rejected_energy_1: Some(0.0004),
                    storage_loss_factor_1: Some(1.35),
                    storage_loss_factor_2: None,
                    rejected_factor_3: None,
                    setpoint_temp: Some(10.),
                    daily_hw_usage: 120.,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::daily_hw_usage_greater_than_zero(json!({"daily_HW_usage": 0})),
                case::rejected_energy_1_at_least_zero(json!({"rejected_energy_1": -0.1})),
                case::storage_loss_factor_1_at_least_zero(json!({"storage_loss_factor_1": -0.1})),
                case::storage_loss_factor_2_at_least_zero(json!({"storage_loss_factor_2": -0.1})),
                case::setpoint_temp_at_least_zero(json!({"setpoint_temp": -9999})),
                case::setpoint_temp_at_most_a_hundred(json!({"setpoint_temp": 101})),
                case::storage_loss_factor_2_invalid_input_for_combis(json!({"storage_loss_factor_2": 2.3})
                ),
                case::rejected_factor_3_invalid_input_for_combis(json!({"rejected_factor_3": 0.0001})
                ),
                case::factor_1_required_when_combi_boiler_is_profile_m(json!({"storage_loss_factor_1": null})
                ),
                case::input_at_least_zero(json!({
                    "separate_DHW_tests": "M&L",
                    "rejected_factor_3": 0.0002,
                    "storage_loss_factor_1": null,
                    "storage_loss_factor_2": -0.1,
                })), // this test is likely misnamed, but follows upstream message in Python
                case::storage_loss_factor_1_invalid_input_for_combis(json!({
                    "separate_DHW_tests": "M&L",
                    "rejected_factor_3": 0.0002,
                    "storage_loss_factor_2": 1.67,
                })),
                case::loss_factors_r1_f2_and_f3_required_when_combi_tested_to_two_profiles(json!({
                    "separate_DHW_tests": "M&L",
                    "rejected_factor_3": 0.0002,
                    "storage_loss_factor_1": null,
                })),
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
                case::setpoint_temp_at_least_zero(json!({"setpoint_temp": -9999})),
                case::setpoint_temp_at_most_a_hundred(json!({"setpoint_temp": 101})),
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
                case::efficiency_greater_than_zero(json!({"efficiency": 0})),
                case::efficiency_at_most_one(json!({"efficiency": 2})),
                case::setpoint_temp_at_least_zero(json!({"setpoint_temp": -9999})),
                case::setpoint_temp_at_most_a_hundred(json!({"setpoint_temp": 101})),
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
                case::init_temp_at_least_zero(json!({"init_temp": -2})),
                case::init_temp_at_most_a_hundred(json!({"init_temp": 101})),
                case::daily_losses_greater_than_zero(json!({"daily_losses": 0})),
                case::heat_exchanger_surface_area_greater_than_zero(json!({"heat_exchanger_surface_area": 0})
                ),
                case::volume_greater_than_zero(json!({"volume": 0})),
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
            case::flowrate_greater_than_zero(json!({"flowrate": 0})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<OtherWaterUse>(valid_example, inputs);
        }
    }

    mod schedule_repeater {
        use super::*;
        use crate::core::schedule::input::{
            ScheduleRepeater, ScheduleRepeaterEntry, ScheduleRepeaterValue,
        };

        #[fixture]
        fn valid_example() -> JsonValue {
            serde_json::to_value(ScheduleRepeater {
                repeat: 10,
                value: ScheduleRepeaterValue::Entry(ScheduleRepeaterEntry::<()>::Null(())),
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::repeat_at_least_one(json!({"repeat": -1})),
        )]
        fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            assert_range_constraints::<ScheduleRepeater<()>>(valid_example, inputs);
        }
    }

    // separate tests for different ScheduleRepeater variants are not needed

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
                case::flowrate_greater_than_zero(json!({"flowrate": 0})),
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
                case::rated_power_greater_than_zero(json!({"rated_power": 0})),
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
            case::battery_state_of_change_at_least_one_map_item(json!({"battery_state_of_charge": {}})
            ),
            case::battery_state_of_change_at_least_one_list_item(json!({"battery_state_of_charge": {"_unmet_demand": []}})
            ),
            case::battery_state_of_change_each_at_most_one(json!({"battery_state_of_charge": {"_unmet_demand": [0, 1, 2, -1]}})
            ),
            case::battery_state_of_change_each_at_least_zero(json!({"battery_state_of_charge": {"_unmet_demand": [0, 0.1, 0.2, -1]}})
            ),
            case::energy_into_battery_from_generation_at_least_one_map_item(json!({"energy_into_battery_from_generation": {}})
            ),
            case::energy_into_battery_from_generation_at_least_one_list_item(json!({"energy_into_battery_from_generation": {"_unmet_demand": []}})
            ),
            case::energy_into_battery_from_grid_at_least_one_map_item(json!({"energy_into_battery_from_grid": {}})
            ),
            case::energy_into_battery_from_grid_at_least_one_list_item(json!({"energy_into_battery_from_grid": {"_unmet_demand": []}})
            ),
            case::energy_out_of_battery_at_least_one_map_item(json!({"energy_out_of_battery": {}})),
            case::energy_out_of_battery_at_least_one_list_item(json!({"energy_out_of_battery": {"_unmet_demand": []}})
            ),
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
            case::non_appliance_demand_24hr_at_least_one_list_item(json!({"non_appliance_demand_24hr": {"_unmet_demand": []}})
            ),
            case::power_timeseries_at_least_one_list_item(json!({"power_timeseries": {"_unmet_demand": []}})
            ),
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
            case::temp_flow_limit_upper_greater_than_zero(json!({"temp_flow_limit_upper": 0})),
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
                orientation360: 123.0.into(),
                pitch: 45.,
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::area_cm2_greater_than_zero(json!({"area_cm2": 0})),
            case::mid_height_air_flow_path_greater_than_zero(json!({"mid_height_air_flow_path": 0})
            ),
            case::orientation360_at_least_zero(json!({"orientation360": -1})),
            case::orientation360_at_most_360(json!({"orientation360": 361})),
            case::pitch_at_least_zero(json!({"pitch": -1})),
            case::pitch_at_most_one(json!({"pitch": 181})),
            case::pressure_difference_ref_greater_than_zero(json!({"pressure_difference_ref": 0})),
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
            case::env_area_greater_than_zero(json!({"env_area": 0})),
            case::test_pressure_greater_than_zero(json!({"test_pressure": 0})),
            case::test_result_greater_than_zero(json!({"test_result": 0})),
            case::ventilation_zone_height_greater_than_zero(json!({"ventilation_zone_height": 0})),
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
            case::duration_greater_than_zero(json!({"duration": 0})),
            case::start_at_least_zero(json!({"start": -1})),
            case::volume_greater_than_zero(json!({"volume": 0})),
            case::temperature_at_least_zero(json!({"temperature": -1})),
            case::temperature_at_most_a_hundred(json!({"temperature": 101})),
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
            case::insulation_thermal_conductivity_greater_than_zero(json!({"insulation_thermal_conductivity": 0})
            ),
            case::internal_diameter_mm_greater_than_zero(json!({"internal_diameter_mm": 0})),
            case::insulation_thickness_mm_at_least_zero(json!({"insulation_thickness_mm": -1})),
            case::length_greater_than_zero(json!({"length": 0})),
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
                    thermal_mass: None,
                    thermal_mass_per_m: None,
                    constant_per_m: None,
                    length: None,
                    constant: Some(1.2),
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
                case::must_specify_length_when_thermal_mass_per_m_provided(json!({
                    "thermal_mass_per_m": 5,
                    "length": null,
                })),
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
                case::equivalent_specific_thermal_mass_greater_than_zero(json!({"equivalent_specific_thermal_mass": 0})
                ),
                case::system_performance_factor_greater_than_zero(json!({"system_performance_factor": 0})
                ),
                case::emitter_floor_area_greater_than_zero(json!({"emitter_floor_area": 0})),
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
                case::height_greater_than_zero(json!({"height": 0})),
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
                case::depth_greater_than_zero(json!({"depth": 0})),
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
            case::flowrate_greater_than_zero(json!({"flowrate": 0})),
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
            case::threshold_charges_at_most_12_items(json!({"threshold_charges": [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]})
            ),
            case::threshold_charges_item_at_most_one(json!({"threshold_charges": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2]})
            ),
            case::threshold_charges_item_at_least_zero(json!({"threshold_charges": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1]})
            ),
            case::threshold_prices_at_least_12_items(json!({"threshold_prices": [0, 1, 1]})),
            case::threshold_prices_at_most_12_items(json!({"threshold_prices": [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]})
            ),
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
                case::power_greater_than_zero(json!({"power": 0})),
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
                    orientation360: 234.0.into(),
                    solar_loop_piping_hlc: 0.5,
                    heater_position: 0.8,
                    thermostat_position: None,
                    control_max: "control max".into(),
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::area_module_greater_than_zero(json!({"area_module": 0})),
                case::collector_mass_flow_rate_greater_than_zero(json!({"collector_mass_flow_rate": 0})
                ),
                case::first_order_hlc_greater_than_zero(json!({"first_order_hlc": 0})),
                case::incidence_angle_modifier_greater_than_zero(json!({"incidence_angle_modifier": 0})
                ),
                case::incidence_angle_modifier_at_most_one(json!({"incidence_angle_modifier": 2})
                ),
                case::peak_collector_efficiency_greater_than_zero(json!({"peak_collector_efficiency": 0})
                ),
                case::peak_collector_efficiency_at_most_one(json!({"peak_collector_efficiency": 2})
                ),
                case::power_pump_at_least_zero(json!({"power_pump": -1})),
                case::power_pump_control_at_least_zero(json!({"power_pump_control": -1})),
                case::second_order_hlc_at_least_zero(json!({"second_order_hlc": -1})),
                case::solar_loop_piping_hlc_greater_than_zero(json!({"solar_loop_piping_hlc": 0})),
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
                case::temp_flow_limit_upper_greater_than_zero(json!({"temp_flow_limit_upper": 0})
                ),
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
                    heat_exchanger_surface_area_declared: Some(1.3),
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
                case::daily_losses_declared_greater_than_zero(json!({"daily_losses_declared": 0})),
                case::heat_exchanger_surface_area_declared_greater_than_zero(json!({"heat_exchanger_surface_area_declared": 0})
                ),
                case::in_use_factor_mismatch_greater_than_zero(json!({"in_use_factor_mismatch": 0})),
                case::power_max_greater_than_zero(json!({"power_max": 0})),
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
            case::daily_losses_greater_than_zero(json!({"daily_losses": 0})),
            case::max_flow_rate_pump_l_per_min_greater_than_zero(json!({"max_flow_rate_pump_l_per_min": 0})
            ),
            case::power_pump_kw_greater_than_zero(json!({"power_pump_kW": 0})),
            case::volume_greater_than_zero(json!({"volume": 0})),
            case::init_temp_at_least_zero(json!({"init_temp": -2})),
            case::temp_usable_at_least_zero(json!({"temp_usable": -2})),
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
                        orientation360: 180.0.into(),
                        pitch: 90.,
                        mid_height_air_flow_path: 3.0,
                    },
                    position_exhaust: MechanicalVentilationPosition {
                        orientation360: 0.0.into(),
                        pitch: 90.,
                        mid_height_air_flow_path: 2.0,
                    },
                },
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::sfp_greater_than_zero(json!({"sfp": -1})),
            case::design_outdoor_air_flow_rate_greater_than_zero(json!({"design_outdoor_air_flow_rate": -1})),
            case::mvhr_eff_at_least_zero(json!({"mvhr_eff": -1})),
            case::mech_vent_mvhr_needs_both_position_intake_and_position_exhaust(json!({"vent_type": "MVHR", "position_intake": null})
            ),
            case::mech_vent_mvhr_needs_both_position_intake_and_position_exhaust(json!({"vent_type": "MVHR", "position_exhaust": null})
            ),
            case::mech_vent_mvhr_uses_position_intake_and_position_exhaust_not_legacy_fields(json!({"orientation360": 234})
            ),
            case::mech_vent_mvhr_uses_position_intake_and_position_exhaust_not_legacy_fields(json!({"pitch": 12})
            ),
            case::mech_vent_mvhr_uses_position_intake_and_position_exhaust_not_legacy_fields(json!({"mid_height_air_flow_path": 4})
            ),
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
                orientation360: 30.0.into(),
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
            case::height_greater_than_zero(json!({"height": 0})),
            case::inverter_peak_power_ac_greater_than_zero(json!({"inverter_peak_power_ac": 0})),
            case::inverter_peak_power_dc_greater_than_zero(json!({"inverter_peak_power_dc": 0})),
            case::peak_power_greater_than_zero(json!({"peak_power": 0})),
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
                    orientation360: 180.0.into(),
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
            case::inverter_peak_power_ac_greater_than_zero(json!({"inverter_peak_power_ac": 0})),
            case::inverter_peak_power_dc_greater_than_zero(json!({"inverter_peak_power_dc": 0})),
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
                orientation360: Orientation360::from(245.),
                base_height: 1.2,
                height: 2.3,
                width: 4.3,
                shading: vec![],
            })
            .unwrap()
        }

        #[rstest(inputs,
            case::peak_power_greater_than_zero(json!({"peak_power": 0})),
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
                start360: 0.0.into(),
                end360: 360.0.into(),
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
            case::efficiency_greater_than_zero(json!({"efficiency": 0})),
            case::cooling_capacity_greater_than_zero(json!({"cooling_capacity": 0})),
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
                case::rated_power_greater_than_zero(json!({"rated_power": 0})),
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
                    design_flow_temp: 55.0,
                    emitters: vec![WetEmitter::Radiator {
                        exponent: 1.2,
                        frac_convective: 0.4,
                        thermal_mass: None,
                        thermal_mass_per_m: None,
                        constant_per_m: None,
                        length: None,
                        constant: Some(0.08),
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
                case::min_and_max_flow_rate_needed_if_variable_flow_true(json!({"variable_flow": true, "min_flow_rate": null, "max_flow_rate": null})
                ),
                case::min_and_max_flow_rate_needed_if_variable_flow_true(json!({"variable_flow": true, "min_flow_rate": 4, "max_flow_rate": null})
                ),
                case::min_and_max_flow_rate_needed_if_variable_flow_true(json!({"variable_flow": true, "min_flow_rate": null, "max_flow_rate": 4})
                ),
                case::design_flow_rate_needed_if_variable_flow_false(json!({"variable_flow": false, "design_flow_rate": null})
                ),
                case::bypass_fraction_recirculated_less_than_one(json!({"bypass_fraction_recirculated": 1})
                ),
                case::bypass_fraction_recirculated_at_least_zero(json!({"bypass_fraction_recirculated": -1})
                ),
                case::design_flow_rate_greater_than_zero(json!({"variable_flow": false, "design_flow_rate": 0})
                ),
                case::max_flow_rate_greater_than_zero(json!({"max_flow_rate": 0})),
                case::min_flow_rate_greater_than_zero(json!({"min_flow_rate": 0})),
                case::temp_diff_emit_dsgn_greater_than_zero(json!({"temp_diff_emit_dsgn": 0})),
                case::thermal_mass_greater_than_zero(json!({"thermal_mass": 0})),
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
                case::first_soc_value_in_dry_core_max_output_must_be_zero(json!({"dry_core_max_output": [[0.1, 0.0], [0.5, 1.5], [1.0, 3.0]]})
                ),
                case::last_soc_value_in_dry_core_max_output_must_be_one(json!({"dry_core_max_output": [[0.0, 0.0], [0.5, 1.5], [0.9, 3.0]]})
                ),
                case::dry_core_max_output_values_must_be_increasing(json!({"dry_core_max_output": [[0.0, 0.0], [0.7, 1.5], [0.5, 3.0], [1.0, 3.5]]})
                ),
                case::first_soc_value_in_dry_core_min_output_must_be_zero(json!({"dry_core_min_output": [[0.1, 0.0], [0.5, 0.02], [1.0, 0.05]]})
                ),
                case::last_soc_value_in_dry_core_min_output_must_be_one(json!({"dry_core_min_output": [[0.0, 0.0], [0.5, 0.02], [0.9, 0.05]]})
                ),
                case::dry_core_min_output_values_must_be_increasing(json!({"dry_core_min_output": [[0.0, 0.0], [0.7, 0.02], [0.5, 0.03], [1.0, 0.05]]})
                ),
                case::dry_core_max_output_must_have_at_least_two_items(json!({"dry_core_max_output": [[]]})
                ),
                case::dry_core_max_output_item_at_least_zero(json!({"dry_core_max_output": [[-1], [-1]]})
                ),
                case::dry_core_min_output_must_have_at_least_two_items(json!({"dry_core_min_output": [[]]})
                ),
                case::dry_core_min_output_item_at_least_zero(json!({"dry_core_min_output": [[-1], [-1]]})
                ),
                case::fan_pwr_at_least_zero(json!({"fan_pwr": -1})),
                case::pwr_in_greate_than_zero(json!({"pwr_in": 0})),
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
                case::logic_type_automatic_requires_temp_charge_cut(json!({"logic_type": "automatic"})
                ),
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
                system_a_efficiencies: Some(vec![]),
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
            case::system_a_efficiencies_item_greater_than_zero(json!({"system_a_efficiencies": [0]})),
            case::system_a_efficiencies_item_at_most_a_hundred(json!({"system_a_efficiencies": [101]})),
            case::flow_rates_item_greater_than_zero(json!({"flow_rates": [0]})),
            case::system_a_utilisation_factor_greater_than_zero(json!({"system_a_utilisation_factor": 0})
            ),
            case::system_a_utilisation_factor_at_most_one(json!({"system_a_utilisation_factor": 2})
            ),
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
                system_a_efficiencies: Some(vec![50., 55., 60., 65., 70.]),
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
                system_a_efficiencies: Some(vec![50., 55., 60., 65., 70.]),
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

        fn test_no_efficiencies_raises_error() {
            let wwhrs = WasteWaterHeatRecoveryDetails {
                _type: MustBe!("WWHRS_Instantaneous"),
                cold_water_source: "header tank".into(),
                flow_rates: vec![5., 7., 9., 11., 13.],
                system_a_efficiencies: None,
                system_a_utilisation_factor: None,
                system_b_efficiencies: None,
                system_b_utilisation_factor: None,
                system_c_efficiencies: None,
                system_c_utilisation_factor: None,
                system_b_efficiency_factor: Default::default(),
                system_c_efficiency_factor: Default::default(),
            };

            assert!(wwhrs.validate().is_err());
        }
    }

    mod building_element {
        use super::*;

        mod common {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(BuildingElement::Opaque {
                    is_unheated_pitched_roof: None,
                    pitch: 90.,
                    orientation360: None,
                    base_height: 0.0,
                    u_value_input: UValueInput::UValue { u_value: 0.3 },
                    areal_heat_capacity: 0.0,
                    mass_distribution_class: MassDistributionClass::D,
                    solar_absorption_coeff: 0.0,
                    area_input: BuildingElementAreaOrHeightWidthInput {
                        area: None,
                        height_and_width: None,
                    },
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::pitch_at_least_zero(json!({"pitch": -1})),
                case::pitch_at_most_zero(json!({"pitch": 181})),
                case::u_value_greater_than_zero(json!({"u_value": 0})),
                case::u_value_greater_than_zero(json!({"u_value": -0.5})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<BuildingElement>(valid_example, inputs);
            }

            // skipping test_u_value_optional as testing something a bit redundant
            // esp as cannot instantiate CommonBase object in Rust
        }

        mod not_ground {
            use super::*;

            #[fixture]
            fn valid_example_with_thermal_resistance() -> JsonValue {
                serde_json::to_value(BuildingElement::Opaque {
                    is_unheated_pitched_roof: None,
                    pitch: 60.,
                    orientation360: None,
                    base_height: 0.0,
                    u_value_input: UValueInput::ThermalResistanceConstruction {
                        thermal_resistance_construction: 4.0,
                    },
                    areal_heat_capacity: 0.0,
                    mass_distribution_class: MassDistributionClass::D,
                    solar_absorption_coeff: 0.0,
                    area_input: BuildingElementAreaOrHeightWidthInput {
                        area: None,
                        height_and_width: None,
                    },
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::u_value_or_thermal_resistance_construction_required(json!({"u_value": null, "thermal_resistance_construction": null})),
                case::thermal_resistance_construction_greater_than_zero(json!({"thermal_resistance_construction": 0})),
                case::thermal_resistance_construction_greater_than_zero(json!({"thermal_resistance_construction": -1.0})),
            )]
            fn test_validate(valid_example_with_thermal_resistance: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<BuildingElement>(
                    valid_example_with_thermal_resistance,
                    inputs,
                );
            }
        }

        mod not_transparent {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(BuildingElement::Opaque {
                    is_unheated_pitched_roof: None,
                    pitch: 60.,
                    orientation360: None,
                    base_height: 0.0,
                    u_value_input: UValueInput::ThermalResistanceConstruction {
                        thermal_resistance_construction: 4.0,
                    },
                    areal_heat_capacity: 15000.0,
                    mass_distribution_class: MassDistributionClass::I,
                    solar_absorption_coeff: 0.0,
                    area_input: BuildingElementAreaOrHeightWidthInput {
                        area: Some(25.0),
                        height_and_width: None,
                    },
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::area_greater_than_zero(json!({"area": 0})),
                case::area_greater_than_zero(json!({"area": -5})),
                case::areal_heat_capacity_greater_than_zero(json!({"areal_heat_capacity": 0})),
                case::areal_heat_capacity_greater_than_zero(json!({"areal_heat_capacity": -1000})),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<BuildingElement>(valid_example, inputs);
            }

            // skipping other python tests here as redundant
        }

        mod exposed_to_solar_radiation {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(BuildingElement::Opaque {
                    is_unheated_pitched_roof: None,
                    pitch: 60.,
                    orientation360: Some(0.0.into()),
                    base_height: 0.5,
                    u_value_input: UValueInput::ThermalResistanceConstruction {
                        thermal_resistance_construction: 4.0,
                    },
                    areal_heat_capacity: 15000.0,
                    mass_distribution_class: MassDistributionClass::I,
                    solar_absorption_coeff: 0.0,
                    area_input: BuildingElementAreaOrHeightWidthInput {
                        area: None,
                        height_and_width: Some(BuildingElementHeightWidthInput {
                            height: 3.0,
                            width: 8.0,
                        }),
                    },
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::orientation360_at_least_zero(json!({"orientation360": -1})),
                case::orientation360_at_most_360(json!({"orientation360": 361})),
                case::base_height_at_least_zero(json!({"base_height": -0.5})),
                case::height_greater_than_zero(json!({"height": 0})),
                case::height_greater_than_zero(json!({"height": -1.0})),
                case::width_greater_than_zero(json!({"width": 0})),
                case::width_greater_than_zero(json!({"width": -2.0})),
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
                    orientation360: 180.0.into(),
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
                    control_window_openable: None,
                    treatment: vec![],
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::max_window_open_area_at_most_area(json!({"max_window_open_area": 9999, "width": 5, "height": 10})
                ),
                case::base_height_at_least_zero(json!({"base_height": -1})),
                case::free_area_height_at_least_zero(json!({"free_area_height": -1})),
                case::g_value_at_least_zero(json!({"g_value": -1})),
                case::max_window_open_area_at_least_zero(json!({"max_window_open_area": -1})),
                case::mid_height_greater_than_zero(json!({"mid_height": 0})),
                case::thermal_resistance_construction_greater_than_zero(json!({"thermal_resistance_construction": 0})
                ),
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
                    orientation360: Some(0.0.into()),
                    base_height: 10.,
                    area_input: BuildingElementAreaOrHeightWidthInput {
                        area: Some(24.0),
                        height_and_width: BuildingElementHeightWidthInput {
                            height: 3.0,
                            width: 8.0,
                        }
                        .into(),
                    },
                    areal_heat_capacity: default_areal_heat_capacity_for_windows(),
                    mass_distribution_class: MassDistributionClass::I,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::base_height_at_least_zero(json!({"base_height": -1})),
                case::solar_absorption_coeff_at_least_zero(json!({"solar_absorption_coeff": -1})),
                case::solar_absorption_coeff_at_most_one(json!({"solar_absorption_coeff": 2})),
                case::thermal_resistance_construction_greater_than_zero(json!({"thermal_resistance_construction": 0, "u_value": null})
                ),
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
                    area: 20.0,
                    u_value_input: UValueInput::UValue { u_value: 1.2 },
                    pitch: 1.2,
                    areal_heat_capacity: 15000.,
                    mass_distribution_class: MassDistributionClass::I,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::thermal_resistance_construction_greater_than_zero(json!({"thermal_resistance_construction": 0, "u_value": null})
                ),
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
                    area: 20.0,
                    u_value_input: UValueInput::UValue { u_value: 0.7 },
                    pitch: 1.2,
                    areal_heat_capacity: 15000.,
                    mass_distribution_class: MassDistributionClass::I,
                    thermal_resistance_unconditioned_space: 1.3,
                })
                .unwrap()
            }

            #[rstest(inputs,
                case::thermal_resistance_construction_greater_than_zero(json!({"thermal_resistance_construction": 0, "u_value": null})
                ),
                case::thermal_resistance_unconditioned_space_greater_than_zero(json!({"thermal_resistance_unconditioned_space": 0})
                ),
            )]
            fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
                assert_range_constraints::<BuildingElement>(valid_example, inputs);
            }
        }

        mod party_wall {
            use super::*;

            #[fixture]
            fn valid_example() -> JsonValue {
                serde_json::to_value(BuildingElement::PartyWall {
                    area: 15.0,
                    u_value_input: UValueInput::UValue { u_value: 0.3 },
                    pitch: 90.,
                    areal_heat_capacity: 120000.,
                    mass_distribution_class: MassDistributionClass::I,
                    party_wall_cavity_data: PartyWallCavityData::UnfilledSealed {
                        party_wall_lining_type: PartyWallLiningType::DryLined,
                    },
                })
                .unwrap()
            }

            // test_invalid_party_wall_lining_type_raises_error is expressed in type system

            #[rstest]
            fn test_party_wall_lining_type_provided(mut valid_example: JsonValue) {
                valid_example["party_wall_cavity_type"] = "unfilled_unsealed".into();
                valid_example["party_wall_lining_type"] = ().into();

                assert!(serde_json::from_value::<BuildingElement>(valid_example).is_err());
            }

            // test_defined_resistance_without_thermal_resistance_raises_error error is expressed in type system
        }

        mod ground {
            use super::*;

            // for deletion during 1.0.0a6
            // mod base {
            //     use super::*;
            //
            //     #[fixture]
            //     fn valid_example() -> JsonValue {
            //         serde_json::to_value(BuildingElement::Ground {
            //             area: 20.,
            //             total_area: 15.,
            //             pitch: 90.,
            //             u_value: 1.4,
            //             thermal_resistance_floor_construction: 0.2,
            //             areal_heat_capacity: 19200.,
            //             mass_distribution_class: MassDistributionClass::I,
            //             perimeter: 16.,
            //             psi_wall_floor_junc: 0.5,
            //             thickness_walls: 0.2,
            //             floor_data: FloorData::SlabNoEdgeInsulation,
            //         })
            //         .unwrap()
            //     }
            //
            //     #[rstest(inputs,
            //         case::calculated_r_vi_greater_than_zero(json!({"u_value": 1, "thermal_resistance_floor_construction": 1})
            //         ),
            //         case::area_greater_than_zero(json!({"area": 0})),
            //         case::thickness_walls_greater_than_zero(json!({"thickness_walls": 0})),
            //         case::total_area_greater_than_zero(json!({"total_area": 0})),
            //         case::perimeter_greater_than_zero(json!({"perimeter": 0})),
            //         case::areal_heat_capacity_greater_than_zero(json!({"areal_heat_capacity": 0})),
            //         case::thermal_resistance_floor_construction_greater_than_zero(json!({"thermal_resistance_floor_construction": 0})
            //         ),
            //         case::u_value_greater_than_zero(json!({"u_value": 0})),
            //         // height_upper_surface only present in suspended floor, so moved that test there
            //         case::pitch_at_least_zero(json!({"pitch": -1})),
            //         case::pitch_at_most_zero(json!({"pitch": 181})),
            //     )]
            //     fn test_validate_range_constraints(valid_example: JsonValue, inputs: JsonValue) {
            //         assert_range_constraints::<BuildingElement>(valid_example, inputs);
            //     }
            // }

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
                            height_upper_surface: 0.5,
                            thermal_transmission_walls: 1.5,
                            area_per_perimeter_vent: 0.0015,
                            shield_fact_location: WindShieldLocation::Average,
                            thermal_resistance_of_insulation: 0.5,
                            edge_insulation: Default::default(),
                        },
                    })
                    .unwrap()
                }

                #[rstest(inputs,
                    case::area_per_perimeter_vent_greater_than_zero(json!({"area_per_perimeter_vent": 0, "u_value": null})
                    ),
                    case::thermal_transm_walls_greater_than_zero(json!({"thermal_transm_walls": 0, "u_value": null})
                    ),
                    case::thermal_resist_insul_greater_than_zero(json!({"thermal_resist_insul": 0})
                    ),
                    case::area_greater_than_zero(json!({"area": 0})),
                    case::thickness_walls_greater_than_zero(json!({"thickness_walls": 0})),
                    case::total_area_greater_than_zero(json!({"total_area": 0})),
                    case::perimeter_greater_than_zero(json!({"perimeter": 0})),
                    case::areal_heat_capacity_greater_than_zero(json!({"areal_heat_capacity": 0})),
                    case::thermal_resistance_floor_construction_greater_than_zero(json!({"thermal_resistance_floor_construction": 0})
                    ),
                    case::u_value_greater_than_zero(json!({"u_value": 0})),
                    case::pitch_at_least_zero(json!({"pitch": -1})),
                    case::pitch_at_most_zero(json!({"pitch": 181})),
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
                            edge_insulation: Default::default(),
                        },
                    })
                    .unwrap()
                }

                #[rstest(inputs,
                    case::thermal_resist_walls_base_greater_than_zero(json!({"thermal_resist_walls_base": 0, "u_value": null})
                    ),
                    case::depth_basement_floor_greater_than_zero(json!({"depth_basement_floor": 0, "u_value": null})
                    ),
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
                            edge_insulation: Default::default(),
                        },
                    })
                    .unwrap()
                }

                #[rstest(inputs,
                    case::thermal_resist_walls_base_greater_than_zero(json!({"thermal_resist_walls_base": 0, "u_value": null})
                    ),
                    case::thermal_transm_envi_base_greater_than_zero(json!({"thermal_transm_envi_base": 0, "u_value": null})
                    ),
                    case::thermal_transm_walls_greater_than_zero(json!({"thermal_transm_walls": 0, "u_value": null})
                    ),
                    case::depth_basement_floor_greater_than_zero(json!({"depth_basement_floor": 0, "u_value": null})
                    ),
                    case::height_basement_walls_greater_than_zero(json!({"height_basement_walls": 0, "u_value": null})
                    ),
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
            case::fan_power_greater_than_zero(json!({"fan_power_W": [0]})),
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
            case::air_temperatures_should_be_at_least_absolute_zero(json!({"air_temperatures": [-274]})
            ),
            case::latitude_should_be_at_least_minus_90(json!({"latitude": -91})),
            case::latitude_should_be_at_most_90(json!({"latitude": 91})),
            case::longitude_should_be_at_least_minus_180(json!({"latitude": -181})),
            case::longitude_should_be_at_most_180(json!({"latitude": 181})),
            case::solar_reflectivity_of_ground_should_be_at_least_zero(json!({"solar_reflectivity_of_ground": [-1]})
            ),
            case::solar_reflectivity_of_ground_should_be_at_most_one(json!({"solar_reflectivity_of_ground": [2]})
            ),
            case::wind_speeds_should_be_at_least_zero(json!({"wind_speeds": [-1]})),
            case::diffuse_horizontal_radiation_should_be_at_least_zero(json!({"diffuse_horizontal_radiation": [-1]})
            ),
            case::direct_beam_radiation_should_be_at_least_zero(json!({"direct_beam_radiation": [-1]})
            ),
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
                wind_directions: Some(vec![0.0.into(); 8760]),
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

        fn test_are_all_fields_set_invalid_lengths1() {
            let valid_example = ExternalConditionsInput {
                air_temperatures: Some(vec![0.; 8760]),
                diffuse_horizontal_radiation: Some(vec![0.; 8760]),
                direct_beam_conversion_needed: false.into(),
                direct_beam_radiation: Some(vec![0.; 8760]),
                latitude: Some(13.9),
                longitude: Some(34.2),
                shading_segments: Some(vec![]),
                solar_reflectivity_of_ground: Some(vec![0.; 8760]),
                wind_directions: Some(vec![0.0.into(); 8760]),
                wind_speeds: Some(vec![0.; 8760]),
            };

            for field in [
                "air_temperatures",
                "wind_speeds",
                "wind_directions",
                "diffuse_horizontal_radiation",
                "direct_beam_radiation",
                "solar_reflectivity_of_ground",
            ] {
                let mut input = serde_json::to_value(valid_example.clone()).unwrap();
                input[field] = [0.0].into();
                let new_input = serde_json::from_value::<ExternalConditionsInput>(input).unwrap();
                assert!(!new_input.are_all_fields_set());
            }
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
            case::daily_losses_greater_than_zero(json!({"daily_losses": 0})),
            case::pump_fixed_flow_rate_greater_than_zero(json!({"pump_fixed_flow_rate": 0})),
            case::pump_power_at_flow_rate_greater_than_zero(json!({"pump_power_at_flow_rate": 0})),
            case::volume_greater_than_zero(json!({"volume": 0})),
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
            case::cop_dhw_greater_than_zero(json!({"cop_dhw": 0})),
            case::energy_input_measured_greater_than_zero(json!({"energy_input_measured": 0})),
            case::hw_tapping_prof_daily_total_greater_than_zero(json!({"hw_tapping_prof_daily_total": 0})
            ),
            case::hw_vessel_loss_daily_greater_than_zero(json!({"hw_vessel_loss_daily": 0})),
            case::power_standby_greater_than_zero(json!({"power_standby": 0})),
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
            case::space_heat_system_list_does_not_allow_duplicates(json!({"SpaceHeatSystem": ["mains", "mains", "other"]})
            ),
            case::space_cool_system_list_does_not_allow_duplicates(json!({"SpaceCoolSystem": ["mains", "mains", "other"]})
            ),
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
                control_closing_irrad: None,
                control_opening_irrad: None,
                control_open: None,
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
            case::pressure_difference_ref_greater_than_zero(json!({"pressure_difference_ref": 0})),
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
                case::edge_thermal_resistance_greater_than_zero(json!({"edge_thermal_resistance": 0})
                ),
                case::width_greater_than_zero(json!({"width": 0})),
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
                case::edge_thermal_resistance_greater_than_zero(json!({"edge_thermal_resistance": 0})
                ),
                case::depth_greater_than_zero(json!({"depth": 0})),
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
            case::mid_height_air_flow_path_greater_than_zero(json!({"mid_height_air_flow_path": 0})
            ),
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
            case::air_flow_rate_greater_than_zero(json!({"air_flow_rate": 0})),
            case::capacity_greater_than_zero(json!({"capacity": 0})),
            case::cop_greater_than_zero(json!({"cop": 0})),
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

    // no need to cover TestUniqueStringList as this is covered by unique_items constraint in serde_valid
}
