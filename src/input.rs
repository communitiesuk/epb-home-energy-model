use crate::core::heating_systems::heat_pump::TestLetter;
use crate::core::schedule::{BooleanSchedule, NumericSchedule};
use crate::external_conditions::{DaylightSavingsConfig, ShadingSegment, WindowShadingObject};
use crate::simulation_time::SimulationTime;
use anyhow::{anyhow, bail};
use indexmap::{Equivalent, IndexMap};
use itertools::Itertools;
use serde::{Deserialize, Deserializer, Serialize};
use serde_enum_str::{Deserialize_enum_str, Serialize_enum_str};
use serde_json::{json, Value};
use serde_repr::Deserialize_repr;
use serde_valid::Validate;
use std::borrow::Borrow;
use std::collections::HashSet;
use std::fmt::{Display, Formatter};
use std::io::{BufReader, Read};
use std::sync::Arc;

pub fn ingest_for_processing(json: impl Read) -> Result<InputForProcessing, anyhow::Error> {
    InputForProcessing::init_with_json(json)
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "PascalCase", deny_unknown_fields)]
pub struct Input {
    pub simulation_time: SimulationTime,
    pub external_conditions: Arc<ExternalConditionsInput>,
    pub(crate) internal_gains: InternalGains,
    pub appliance_gains: ApplianceGains,
    pub cold_water_source: ColdWaterSourceInput,
    pub energy_supply: EnergySupplyInput,
    #[serde(deserialize_with = "deserialize_control")]
    pub control: Control,
    pub hot_water_source: HotWaterSource,
    pub hot_water_demand: HotWaterDemand,
    #[serde(rename = "Events")]
    pub(crate) water_heating_events: WaterHeatingEvents,
    pub(crate) space_heat_system: Option<SpaceHeatSystem>,
    pub space_cool_system: Option<SpaceCoolSystem>,
    pub ventilation: Option<Ventilation>,
    pub zone: ZoneDictionary,
    // following fields marked as possibly dead code are likely to be used by wrappers, but worth checking when compiling input schema
    #[allow(dead_code)]
    #[serde(rename = "PartGcompliance")]
    part_g_compliance: Option<bool>,
    #[allow(dead_code)]
    #[serde(rename = "PartO_active_cooling_required")]
    part_o_active_cooling_required: Option<bool>,
    #[allow(dead_code)]
    ground_floor_area: Option<f64>,
    number_of_bedrooms: Option<usize>,
    #[allow(dead_code)]
    number_of_wet_rooms: Option<usize>,
    heating_control_type: Option<HeatingControlType>,
    #[allow(dead_code)]
    #[serde(rename = "WaterHeatSchedDefault")]
    default_water_heating_schedule: Option<WaterHeatingSchedule>,
    pub heat_source_wet: Option<HeatSourceWet>,
    #[serde(rename = "WWHRS")]
    pub waste_water_heat_recovery: Option<WasteWaterHeatRecovery>,
    pub on_site_generation: Option<OnSiteGeneration>,
    #[serde(rename = "Window_Opening_For_Cooling")]
    pub window_opening_for_cooling: Option<WindowOpeningForCooling>,
    pub general: General,
    pub infiltration_ventilation: InfiltrationVentilation,
    #[serde(rename = "Appliances")]
    pub(crate) appliances: Option<IndexMap<ApplianceKey, ApplianceEntry>>,
    pub tariff: Option<Tariff>,
}

#[derive(Debug, Default, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct ExternalConditionsInput {
    pub air_temperatures: Option<Vec<f64>>,
    pub wind_speeds: Option<Vec<f64>>,
    pub wind_directions: Option<Vec<f64>>,
    // check upstream whether anything uses this
    #[serde(rename = "ground_temperatures")]
    pub(crate) _ground_temperatures: Option<Vec<f64>>,
    pub diffuse_horizontal_radiation: Option<Vec<f64>>,
    pub direct_beam_radiation: Option<Vec<f64>>,
    pub solar_reflectivity_of_ground: Option<Vec<f64>>,
    pub latitude: Option<f64>,
    pub longitude: Option<f64>,
    pub timezone: Option<u32>,
    pub start_day: Option<u32>,
    pub end_day: Option<u32>,
    pub time_series_step: Option<f64>,
    pub january_first: Option<u32>,
    pub daylight_savings: Option<DaylightSavingsConfig>,
    pub leap_day_included: Option<bool>,
    pub direct_beam_conversion_needed: Option<bool>,
    pub shading_segments: Vec<ShadingSegment>,
}

#[derive(Debug, Default, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct InternalGains {
    #[serde(alias = "total internal gains")]
    pub total_internal_gains: Option<InternalGainsDetails>,
    #[serde(rename = "metabolic gains")]
    pub metabolic_gains: Option<InternalGainsDetails>,
    pub evaporative_losses: Option<InternalGainsDetails>,
    #[serde(rename = "ColdWaterLosses")]
    pub cold_water_losses: Option<InternalGainsDetails>,
    pub other: Option<InternalGainsDetails>,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct InternalGainsDetails {
    pub start_day: u32,
    pub time_series_step: f64,
    pub(crate) schedule: NumericSchedule,
}

pub type ApplianceGains = IndexMap<String, ApplianceGainsDetails>;

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct ApplianceGainsDetails {
    #[serde(rename(deserialize = "type"))]
    _gain_type: Option<String>,
    pub start_day: u32,
    pub time_series_step: f64,
    pub gains_fraction: f64,
    #[serde(rename = "EnergySupply")]
    pub energy_supply: EnergySupplyType,
    pub(crate) schedule: Option<NumericSchedule>,
    // In the Python code these fields are
    // set in the FHS wrapper
    #[serde(rename = "Standby")]
    pub standby: Option<f64>,
    #[serde(rename = "Events")]
    pub events: Option<Vec<ApplianceGainsDetailsEvent>>,
    #[serde(rename = "loadshifting")]
    pub load_shifting: Option<ApplianceLoadShifting>,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct ApplianceGainsDetailsEvent {
    pub start: f64,
    pub duration: f64,
    pub demand_w: f64,
}

#[derive(Clone, Copy, Debug, Deserialize_enum_str, Eq, Hash, PartialEq, Serialize_enum_str)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum EnergySupplyKey {
    #[serde(rename = "mains elec")]
    MainsElectricity,
    #[serde(rename = "mains gas")]
    MainsGas,
    #[serde(rename = "bulk LPG")]
    BulkLpg,
    #[serde(rename = "heat network")]
    HeatNetwork,
}

impl Equivalent<str> for EnergySupplyKey {
    fn equivalent(&self, key: &str) -> bool {
        if let Ok(key) = serde_json::from_value::<EnergySupplyKey>(json!(key)) {
            key == *self
        } else {
            false
        }
    }
}

impl Borrow<str> for EnergySupplyKey {
    fn borrow(&self) -> &str {
        match serde_json::to_value(self)
            .expect("Expected an EnergySupplyKey to be turned into a string")
        {
            Value::String(string) => match string.as_str() {
                "mains elec" => "mains elec",
                "mains gas" => "mains gas",
                "bulk LPG" => "bulk LPG",
                "heat network" => "heat network",
                _ => unreachable!(),
            },
            _ => unreachable!(),
        }
    }
}

impl From<&EnergySupplyKey> for String {
    fn from(value: &EnergySupplyKey) -> Self {
        serde_json::to_value(value)
            .unwrap()
            .as_str()
            .unwrap()
            .to_string()
    }
}

#[cfg(test)]
mod energy_supply_key_tests {
    use crate::input::{EnergySupplyDetails, EnergySupplyKey, FuelType};
    use indexmap::IndexMap;
    use rstest::*;

    #[rstest]
    fn deserialize_energy_supply_key() {
        assert_eq!(
            "bulk LPG".parse::<EnergySupplyKey>().unwrap(),
            EnergySupplyKey::BulkLpg
        );
        assert_eq!(
            "mains elec".parse::<EnergySupplyKey>().unwrap(),
            EnergySupplyKey::MainsElectricity
        );
    }

    #[rstest]
    fn use_energy_supply_key_to_get_supply_details() {
        let supplies = IndexMap::from([
            (
                EnergySupplyKey::MainsElectricity,
                EnergySupplyDetails {
                    fuel: FuelType::Electricity,
                    diverter: None,
                    electric_battery: None,
                    factor: None,
                    priority: None,
                    is_export_capable: None,
                },
            ),
            (
                EnergySupplyKey::MainsGas,
                EnergySupplyDetails {
                    fuel: FuelType::Electricity,
                    diverter: None,
                    electric_battery: None,
                    factor: None,
                    priority: None,
                    is_export_capable: None,
                },
            ),
        ]);
        assert!(supplies.get(&EnergySupplyKey::MainsElectricity).is_some());
        assert!(supplies
            .get(&("mains gas".parse::<EnergySupplyKey>().unwrap()))
            .is_some());
    }
}

pub type EnergySupplyInput = IndexMap<EnergySupplyKey, EnergySupplyDetails>;

#[derive(Clone, Debug, Deserialize, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "fuel", deny_unknown_fields, rename_all = "snake_case")]
pub struct EnergySupplyDetails {
    pub fuel: FuelType,
    pub diverter: Option<EnergyDiverter>,
    #[serde(rename = "ElectricBattery")]
    pub electric_battery: Option<ElectricBattery>,
    pub factor: Option<CustomEnergySourceFactor>,
    pub priority: Option<Vec<SecondarySupplyType>>,
    pub is_export_capable: Option<bool>,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum SecondarySupplyType {
    ElectricBattery,
    #[serde(rename = "diverter")]
    Diverter,
}

// TODO clarify further
// It's not completely clear at the moment what the difference between fuel type and energy supply type is,
// but electricity and gas each seem to be indicated using different strings between fuel and energy supply
// in the input examples, so keeping them separate for the time being
// (It's also hard to see some of these as types of fuel)
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "snake_case")]
pub enum FuelType {
    Electricity,
    MainsGas,
    Custom,
    #[serde(rename = "LPG_bulk")]
    LpgBulk,
    UnmetDemand,
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

#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum EnergySupplyType {
    #[serde(rename = "mains elec")]
    Electricity,
    #[serde(rename = "mains gas", alias = "mains_gas")]
    MainsGas,
    #[serde(rename = "unmet_demand")]
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

#[derive(Clone, Debug, Deserialize, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "PascalCase")]
#[serde(deny_unknown_fields)]
pub struct EnergyDiverter {
    pub storage_tank: StorageTankType,
    pub heat_source: DiverterHeatSourceType,
}

impl Default for EnergyDiverter {
    fn default() -> Self {
        Self {
            storage_tank: StorageTankType::HotWaterCylinder,
            heat_source: DiverterHeatSourceType::Immersion,
        }
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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

#[derive(Clone, Debug, Deserialize, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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

#[derive(Clone, Debug, Deserialize, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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
}

#[derive(Copy, Clone, Debug, Deserialize, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "lowercase")]
pub enum BatteryLocation {
    Inside,
    Outside,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct CustomEnergySourceFactor {
    #[serde(rename = "Emissions Factor kgCO2e/kWh")]
    pub emissions: f64,
    #[serde(rename = "Emissions Factor kgCO2e/kWh including out-of-scope emissions")]
    pub emissions_including_out_of_scope: f64,
    #[serde(rename = "Primary Energy Factor kWh/kWh delivered")]
    pub primary_energy_factor: f64,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct ColdWaterSourceInput {
    #[serde(rename = "mains water")]
    pub mains_water: Option<ColdWaterSourceDetails>,
    #[serde(rename = "header tank")]
    pub header_tank: Option<ColdWaterSourceDetails>,
}

impl ColdWaterSourceInput {
    fn has_header_tank(&self) -> bool {
        self.header_tank.is_some()
    }

    fn set_cold_water_source_details_by_key(&mut self, key: &str, source: ColdWaterSourceDetails) {
        match key {
            "mains water" => {
                self.mains_water = Some(source);
            }
            "header tank" => {
                self.header_tank = Some(source);
            }
            _ => unreachable!(),
        }
    }
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct ColdWaterSourceDetails {
    pub start_day: u32,
    pub temperatures: Vec<f64>,
    pub time_series_step: f64,
}

pub(crate) type CoreControls = Vec<HeatSourceControl>;

pub(crate) type ExtraControls = IndexMap<String, ControlDetails>;

#[derive(Debug)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub struct Control {
    pub(crate) core: CoreControls,
    pub(crate) extra: ExtraControls,
}

// specialised deserialisation logic for converting a map of controls into a list of HeatSourceControl structs
fn deserialize_control<'de, D>(deserializer: D) -> Result<Control, D::Error>
where
    D: Deserializer<'de>,
{
    let map: IndexMap<String, ControlDetails> = Deserialize::deserialize(deserializer)?;
    let mut core: CoreControls = Default::default();
    let mut extra: ExtraControls = Default::default();
    for (control_type, control_details) in map {
        match control_type.as_str() {
            // following strings need to be in sync with HeatSourceControlType known values
            "hw timer" => {
                core.push(HeatSourceControl::HotWaterTimer(control_details));
            }
            "window opening" => {
                core.push(HeatSourceControl::WindowOpening(control_details));
            }
            // there are some extra control definitions from time to time called things like
            // "hw timer 2" and "zone 1 radiators timer" - can only presume now to store these keys as-is
            // and perhaps match on other references
            other => {
                extra.insert(other.to_string(), control_details);
            }
        }
    }
    Ok(Control { core, extra })
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "type", deny_unknown_fields)]
pub(crate) enum ControlDetails {
    #[serde(rename = "OnOffTimeControl")]
    OnOffTime {
        start_day: u32,
        time_series_step: f64,
        #[allow(dead_code)]
        logic_type: Option<ControlLogicType>,
        schedule: BooleanSchedule,
    },
    #[serde(rename = "OnOffCostMinimisingTimeControl")]
    OnOffCostMinimisingTime {
        start_day: u32,
        time_series_step: f64,
        #[allow(dead_code)]
        logic_type: Option<ControlLogicType>,
        time_on_daily: Option<f64>,
        schedule: NumericSchedule,
    },
    #[serde(rename = "SetpointTimeControl")]
    SetpointTime {
        start_day: u32,
        time_series_step: f64,
        advanced_start: Option<f64>,
        #[allow(dead_code)]
        logic_type: Option<ControlLogicType>,
        setpoint_min: Option<f64>,
        setpoint_max: Option<f64>,
        default_to_max: Option<bool>,
        schedule: NumericSchedule,
    },
    #[serde(rename = "ToUChargeControl")]
    ToUCharge {
        start_day: u32,
        time_series_step: f64,
        #[allow(dead_code)]
        logic_type: Option<ControlLogicType>,
        charge_level: Option<ChargeLevel>,
        #[serde(rename = "target_charge")]
        _target_charge: Option<f64>,
        schedule: BooleanSchedule,
    },
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(untagged)]
pub(crate) enum ChargeLevel {
    Single(f64),
    List(Vec<f64>),
}

impl ControlDetails {
    pub(crate) fn start_day(&self) -> u32 {
        match self {
            ControlDetails::OnOffTime { start_day, .. } => *start_day,
            ControlDetails::OnOffCostMinimisingTime { start_day, .. } => *start_day,
            ControlDetails::SetpointTime { start_day, .. } => *start_day,
            ControlDetails::ToUCharge { start_day, .. } => *start_day,
        }
    }

    pub(crate) fn time_series_step(&self) -> f64 {
        match self {
            ControlDetails::OnOffTime {
                time_series_step, ..
            } => *time_series_step,
            ControlDetails::OnOffCostMinimisingTime {
                time_series_step, ..
            } => *time_series_step,
            ControlDetails::SetpointTime {
                time_series_step, ..
            } => *time_series_step,
            ControlDetails::ToUCharge {
                time_series_step, ..
            } => *time_series_step,
        }
    }

    pub(crate) fn numeric_schedule(&self) -> anyhow::Result<&NumericSchedule> {
        match self {
            ControlDetails::OnOffCostMinimisingTime { schedule, .. } => Ok(schedule),
            ControlDetails::SetpointTime { schedule, .. } => Ok(schedule),
            _ => Err(anyhow::anyhow!("Numeric schedule was not available")),
        }
    }
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum ControlLogicType {
    Manual,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct HotWaterSource {
    #[serde(rename = "hw cylinder")]
    pub hot_water_cylinder: HotWaterSourceDetails,
}

impl HotWaterSource {
    fn source_keys(&self) -> Vec<String> {
        vec!["hw cylinder".to_string()]
    }

    fn hot_water_source_for_processing(
        &mut self,
        source_key: &str,
    ) -> &mut impl HotWaterSourceDetailsForProcessing {
        match source_key {
            "hw cylinder" => &mut self.hot_water_cylinder,
            _ => unreachable!("requested a source key {source_key} that is not known"),
        }
    }

    pub fn as_index_map(&self) -> IndexMap<String, HotWaterSourceDetails> {
        IndexMap::from([("hw cylinder".to_string(), self.hot_water_cylinder.clone())])
    }
}

#[derive(Clone, Deserialize_enum_str, PartialEq, Debug)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "type", deny_unknown_fields)]
pub enum HotWaterSourceDetails {
    StorageTank {
        volume: f64,
        daily_losses: f64,
        heat_exchanger_surface_area: Option<f64>,
        min_temp: f64,
        setpoint_temp: f64,
        #[serde(rename = "Control_hold_at_setpnt")]
        control_hold_at_setpoint: Option<String>,
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceType,
        #[serde(rename = "HeatSource")]
        heat_source: IndexMap<String, HeatSource>,
        primary_pipework: Option<Vec<WaterPipework>>,
    },
    CombiBoiler {
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceType,
        #[serde(rename = "HeatSourceWet")]
        heat_source_wet: HeatSourceWetType,
        #[serde(rename = "Control")]
        control: HeatSourceControlType,
        #[serde(rename = "separate_DHW_tests")]
        separate_dhw_tests: BoilerHotWaterTest,
        rejected_energy_1: f64,
        fuel_energy_2: f64,
        rejected_energy_2: f64,
        storage_loss_factor_2: f64,
        rejected_factor_3: f64,
        setpoint_temp: Option<f64>,
        #[serde(rename = "daily_HW_usage")]
        daily_hot_water_usage: f64,
    },
    #[serde(rename(deserialize = "HIU"))]
    Hiu {
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceType,
        #[serde(rename = "HeatSourceWet")]
        heat_source_wet: HeatSourceWetType,
        #[serde(rename = "Control")]
        control: HeatSourceControlType,
        setpoint_temp: f64,
    },
    PointOfUse {
        efficiency: f64,
        #[serde(rename = "EnergySupply")]
        energy_supply: EnergySupplyType,
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceType,
        setpoint_temp: Option<f64>,
    },
    HeatBattery {
        // tbc
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
    #[allow(dead_code)]
    fn is_heat_battery(&self) -> bool;
    fn set_control_hold_at_setpoint(&mut self, control_name: impl Into<String>);
    fn set_control_name_for_heat_sources(
        &mut self,
        control_name: impl Into<String>,
    ) -> anyhow::Result<()>;
    fn set_min_temp_and_setpoint_temp_if_storage_tank(&mut self, min_temp: f64, setpoint_temp: f64);
    fn set_setpoint_temp(&mut self, setpoint_temp: f64);
}

impl HotWaterSourceDetailsForProcessing for HotWaterSourceDetails {
    fn is_storage_tank(&self) -> bool {
        matches!(self, Self::StorageTank { .. })
    }

    fn is_combi_boiler(&self) -> bool {
        matches!(self, Self::CombiBoiler { .. })
    }

    fn is_hiu(&self) -> bool {
        matches!(self, Self::Hiu { .. })
    }

    fn is_point_of_use(&self) -> bool {
        matches!(self, Self::PointOfUse { .. })
    }

    #[allow(dead_code)]
    fn is_heat_battery(&self) -> bool {
        matches!(self, Self::HeatBattery { .. })
    }

    fn set_control_hold_at_setpoint(&mut self, control_name: impl Into<String>) {
        if let Self::StorageTank {
            ref mut control_hold_at_setpoint,
            ..
        } = self
        {
            *control_hold_at_setpoint = Some(control_name.into());
        }
    }

    fn set_control_name_for_heat_sources(
        &mut self,
        control_name: impl Into<String>,
    ) -> anyhow::Result<()> {
        if let Self::StorageTank {
            ref mut heat_source,
            ..
        } = self
        {
            let control_type_for_heat_sources = control_name.into().try_into()?;
            for heat_source in heat_source.values_mut() {
                heat_source.set_control(control_type_for_heat_sources);
            }
        }

        Ok(())
    }

    fn set_min_temp_and_setpoint_temp_if_storage_tank(
        &mut self,
        min_temp: f64,
        setpoint_temp: f64,
    ) {
        if let HotWaterSourceDetails::StorageTank {
            min_temp: ref mut min_temp_store,
            setpoint_temp: ref mut setpoint_temp_store,
            ..
        } = self
        {
            *min_temp_store = min_temp;
            *setpoint_temp_store = setpoint_temp;
        }
    }

    fn set_setpoint_temp(&mut self, setpoint_temp: f64) {
        match self {
            HotWaterSourceDetails::StorageTank {
                setpoint_temp: ref mut source_setpoint_temp,
                ..
            } => {
                *source_setpoint_temp = setpoint_temp;
            }
            HotWaterSourceDetails::CombiBoiler {
                setpoint_temp: ref mut source_setpoint_temp,
                ..
            } => {
                *source_setpoint_temp = Some(setpoint_temp);
            }
            HotWaterSourceDetails::Hiu {
                setpoint_temp: ref mut source_setpoint_temp,
                ..
            } => {
                *source_setpoint_temp = setpoint_temp;
            }
            HotWaterSourceDetails::PointOfUse {
                setpoint_temp: ref mut source_setpoint_temp,
                ..
            } => {
                *source_setpoint_temp = Some(setpoint_temp);
            }
            HotWaterSourceDetails::HeatBattery {} => {}
        }
    }
}

#[derive(Copy, Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum ColdWaterSourceType {
    #[serde(rename = "mains water")]
    MainsWater,
    #[serde(rename = "header tank")]
    HeaderTank,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum HeatSourceWetType {
    #[serde(rename = "boiler")]
    Boiler,
    HeatNetwork,
    #[serde(rename = "hp")]
    HeatPump,
}

impl HeatSourceWetType {
    /// Convert the type to a canonical string based on the input format to be used in e.g. energy supply names
    pub fn to_canonical_string(&self) -> String {
        serde_json::to_value(self)
            .unwrap()
            .as_str()
            .unwrap()
            .to_owned()
    }
}

#[derive(Clone, Copy, Debug, Deserialize_enum_str)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum HeatSourceControlType {
    #[serde(rename = "hw timer")]
    HotWaterTimer,
    #[serde(rename = "window opening")]
    WindowOpening,
    #[serde(rename = "WindowOpening_LivingRoom")]
    WindowOpeningLivingRoom,
    #[serde(rename = "WindowOpening_RestOfDwelling")]
    WindowOpeningRestOfDwelling,
    #[serde(rename = "always off")]
    AlwaysOff,
}

#[derive(Debug)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum HeatSourceControl {
    HotWaterTimer(ControlDetails),
    WindowOpening(ControlDetails),
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "type", deny_unknown_fields)]
pub enum HeatSource {
    ImmersionHeater {
        power: f64,
        #[serde(rename = "EnergySupply")]
        energy_supply: EnergySupplyType,
        #[serde(rename = "Control")]
        control: Option<HeatSourceControlType>,
        heater_position: f64,
        thermostat_position: f64,
    },
    SolarThermalSystem {
        #[serde(rename(deserialize = "sol_loc"))]
        solar_cell_location: SolarCellLocation,
        area_module: f64,
        modules: usize,
        peak_collector_efficiency: f64,
        incidence_angle_modifier: f64,
        first_order_hlc: f64,
        second_order_hlc: f64,
        collector_mass_flow_rate: f64,
        power_pump: f64,
        power_pump_control: f64,
        #[serde(rename = "EnergySupply")]
        energy_supply: EnergySupplyType,
        tilt: f64,
        #[serde(rename = "orientation360")]
        #[serde(deserialize_with = "deserialize_orientation")]
        orientation: f64,
        solar_loop_piping_hlc: f64,
        heater_position: f64,
        thermostat_position: f64,
        control: Option<HeatSourceControlType>,
    },
    #[serde(rename(deserialize = "HeatSourceWet"))]
    Wet {
        name: String,
        temp_flow_limit_upper: Option<f64>,
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: Option<ColdWaterSourceType>,
        #[serde(rename = "EnergySupply")]
        energy_supply: EnergySupplyType,
        #[serde(rename = "Control")]
        control: Option<HeatSourceControlType>,
        heater_position: f64,
        thermostat_position: f64,
        #[serde(rename = "temp_return")]
        temperature_return: Option<f64>,
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
        energy_supply: EnergySupplyType,
        #[serde(rename = "Control")]
        control: HeatSourceControlType,
        heater_position: f64,
        thermostat_position: f64,
    },
}

impl HeatSource {
    pub fn heater_position(&self) -> f64 {
        match self {
            HeatSource::ImmersionHeater {
                heater_position, ..
            } => *heater_position,
            HeatSource::SolarThermalSystem {
                heater_position, ..
            } => *heater_position,
            HeatSource::Wet {
                heater_position, ..
            } => *heater_position,
            HeatSource::HeatPumpHotWaterOnly {
                heater_position, ..
            } => *heater_position,
        }
    }

    pub fn thermostat_position(&self) -> f64 {
        match self {
            HeatSource::ImmersionHeater {
                thermostat_position,
                ..
            } => *thermostat_position,
            HeatSource::SolarThermalSystem {
                thermostat_position,
                ..
            } => *thermostat_position,
            HeatSource::Wet {
                thermostat_position,
                ..
            } => *thermostat_position,
            HeatSource::HeatPumpHotWaterOnly {
                thermostat_position,
                ..
            } => *thermostat_position,
        }
    }

    pub fn energy_supply_type(&self) -> EnergySupplyType {
        match self {
            HeatSource::ImmersionHeater { energy_supply, .. } => *energy_supply,
            HeatSource::SolarThermalSystem { energy_supply, .. } => *energy_supply,
            HeatSource::Wet { energy_supply, .. } => *energy_supply,
            HeatSource::HeatPumpHotWaterOnly { energy_supply, .. } => *energy_supply,
        }
    }

    pub fn set_control(&mut self, control_type: HeatSourceControlType) {
        match self {
            HeatSource::ImmersionHeater {
                ref mut control, ..
            } => {
                *control = Some(control_type);
            }
            HeatSource::SolarThermalSystem {
                ref mut control, ..
            } => {
                *control = Some(control_type);
            }
            HeatSource::Wet {
                ref mut control, ..
            } => {
                *control = Some(control_type);
            }
            HeatSource::HeatPumpHotWaterOnly {
                ref mut control, ..
            } => {
                *control = control_type;
            }
        }
    }
}

#[derive(Copy, Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum SolarCellLocation {
    #[serde(rename = "OUT")]
    Out,
    #[serde(rename = "HS")]
    Hs,
    #[serde(rename = "NHS")]
    Nhs,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct HeatPumpHotWaterTestData {
    #[serde(rename = "L")]
    pub l: Option<HeatPumpHotWaterOnlyTestDatum>,
    #[serde(rename = "M")]
    pub m: HeatPumpHotWaterOnlyTestDatum,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct HeatPumpHotWaterOnlyTestDatum {
    // CoP measured during EN 16147 test
    pub cop_dhw: f64,
    // daily energy requirement (kWh/day) for tapping profile used for test
    #[serde(rename = "hw_tapping_prof_daily_total")]
    pub hw_tapping_prof_daily: f64,
    // electrical input energy (kWh) measured in EN 16147 test over 24 hrs
    pub energy_input_measured: f64,
    // standby power (W) measured in EN 16147 test
    pub power_standby: f64,
    // daily hot water vessel heat loss
    // (kWh/day) for a 45 K temperature difference between vessel
    // and surroundings, tested in accordance with BS 1566 or
    // EN 12897 or any equivalent standard. Vessel must be same
    // as that used during EN 16147 test
    pub hw_vessel_loss_daily: f64,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct WaterPipeworkSimple {
    pub location: WaterPipeworkLocation,
    pub internal_diameter_mm: f64,
    pub length: f64,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum WaterPipeworkLocation {
    #[serde(rename = "internal")]
    Internal,
    #[serde(rename = "external")]
    External,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum WaterPipeContentsType {
    #[serde(rename = "water")]
    Water,
    #[serde(rename = "air")]
    Air,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct HotWaterDemand {
    #[serde(rename = "Shower")]
    pub shower: Option<Showers>,
    #[serde(rename = "Bath")]
    pub bath: Option<Baths>,
    #[serde(rename = "Other")]
    pub other_water_use: Option<OtherWaterUses>,
    #[serde(rename = "Distribution")]
    pub water_distribution: Option<WaterDistribution>,
}

#[derive(Clone, Debug, Default, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub struct Showers(pub IndexMap<String, Shower>);

impl Showers {
    /// Provide shower field names as strings.
    pub fn keys(&self) -> Vec<String> {
        self.0.keys().cloned().collect()
    }

    pub fn name_refers_to_instant_electric_shower(&self, name: &str) -> bool {
        self.0.get(name).map_or(false, |shower| {
            matches!(shower, Shower::InstantElectricShower { .. })
        })
    }
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields, tag = "type")]
pub enum Shower {
    MixerShower {
        flowrate: f64,
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceType,
        #[serde(rename = "WWHRS")]
        waste_water_heat_recovery: Option<String>,
    },
    #[serde(rename = "InstantElecShower")]
    InstantElectricShower {
        rated_power: f64,
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceType,
        #[serde(rename = "EnergySupply")]
        energy_supply: EnergySupplyType,
    },
}

#[derive(Clone, Debug, Default, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct Baths(pub IndexMap<String, BathDetails>);

impl Baths {
    /// Provide bath field names as strings.
    pub fn keys(&self) -> Vec<String> {
        self.0.keys().cloned().collect()
    }

    pub fn size_for_field(&self, field: &str) -> Option<f64> {
        self.0.get(field).map(|bath| bath.size)
    }

    pub fn flowrate_for_field(&self, field: &str) -> Option<f64> {
        self.0.get(field).map(|bath| bath.flowrate)
    }
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct BathDetails {
    pub size: f64,
    #[serde(rename = "ColdWaterSource")]
    pub cold_water_source: ColdWaterSourceType,
    pub flowrate: f64,
}

#[derive(Clone, Debug, Default, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct OtherWaterUses(pub IndexMap<String, OtherWaterUseDetails>);

impl OtherWaterUses {
    /// Provide other water use field names as strings.
    pub fn keys(&self) -> Vec<String> {
        self.0.keys().cloned().collect()
    }

    pub fn flowrate_for_field(&self, field: &str) -> Option<f64> {
        self.0.get(field).map(|other| other.flowrate)
    }
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct OtherWaterUseDetails {
    pub flowrate: f64,
    #[serde(rename = "ColdWaterSource")]
    pub cold_water_source: ColdWaterSourceType,
}

pub type WaterDistribution = Vec<WaterPipeworkSimple>;

#[derive(Debug, Default, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) struct WaterHeatingEvents(
    pub(crate) IndexMap<WaterHeatingEventType, IndexMap<String, Vec<WaterHeatingEvent>>>,
);

impl WaterHeatingEvents {
    fn add_event_for_type_and_name(
        &mut self,
        event_type: WaterHeatingEventType,
        name: &str,
        event: WaterHeatingEvent,
    ) {
        self.0
            .entry(event_type)
            .or_default()
            .entry(name.to_string())
            .or_default()
            .push(event);
    }
}

pub(crate) trait WaterHeatingEventsForProcessing {
    fn water_heating_events_of_types(
        &self,
        event_types: &[WaterHeatingEventType],
    ) -> Vec<WaterHeatingEvent>;
}

impl WaterHeatingEventsForProcessing for WaterHeatingEvents {
    fn water_heating_events_of_types(
        &self,
        event_types: &[WaterHeatingEventType],
    ) -> Vec<WaterHeatingEvent> {
        self.0
            .iter()
            .filter(|(event_type, _)| event_types.contains(event_type))
            .flat_map(|(_, events)| events.values())
            .flatten()
            .copied()
            .collect_vec()
    }
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct WaterHeatingEvent {
    pub start: f64,
    pub duration: Option<f64>,
    pub temperature: f64,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum WaterHeatingEventType {
    Shower,
    Bath,
    Other,
}

pub(crate) type SpaceHeatSystem = IndexMap<String, SpaceHeatSystemDetails>;

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields, tag = "type")]
pub(crate) enum SpaceHeatSystemDetails {
    #[serde(rename = "InstantElecHeater")]
    InstantElectricHeater {
        temp_setback: Option<f64>,
        rated_power: f64,
        #[serde(rename = "EnergySupply")]
        energy_supply: EnergySupplyType,
        #[serde(rename = "Control")]
        control: Option<String>,
        // not sure what the possible options are here yet
        frac_convective: f64,
        #[serde(rename = "Zone")]
        _zone: Option<String>,
    },
    #[serde(rename = "ElecStorageHeater")]
    #[allow(dead_code)]
    ElectricStorageHeater {
        temp_charge_cut: f64,
        rated_power: f64,
        rated_power_instant: f64,
        air_flow_type: String,
        // don't know what the possible values are here yet
        temp_dis_safe: f64,
        thermal_mass: f64,
        frac_convective: f64,
        #[serde(rename = "U_ins")]
        u_ins: f64,
        mass_core: f64,
        c_pcore: f64,
        temp_core_target: f64,
        #[serde(rename = "A_core")]
        a_core: f64,
        c_wall: f64,
        n_wall: f64,
        thermal_mass_wall: f64,
        fan_pwr: f64,
        n_units: u32,
        #[serde(rename = "EnergySupply")]
        energy_supply: EnergySupplyType,
        #[serde(rename = "Control")]
        control: Option<String>,
        // don't know possible options here
        #[serde(rename = "ControlCharger")]
        control_charger: String,
        // don't know possible options here
        #[serde(rename = "Zone")]
        zone: String, // think these are just arbitrary names?
    },
    #[allow(dead_code)]
    WetDistribution {
        advanced_start: Option<f64>,
        thermal_mass: f64,
        c: f64,
        n: f64,
        temp_diff_emit_dsgn: f64,
        frac_convective: f64,
        // unclear which values are possible here
        #[serde(rename = "HeatSource")]
        heat_source: SpaceHeatSystemHeatSource,
        #[serde(rename = "Control")]
        control: Option<String>,
        // check upstream if this is used
        ecodesign_controller: EcoDesignController,
        design_flow_temp: i32,
        #[serde(rename = "Zone")]
        zone: String, // as above, these are likely arbitrary names
    },
    #[allow(dead_code)]
    WarmAir {
        temp_diff_emit_dsgn: f64,
        frac_convective: f64,
        #[serde(rename = "HeatSource")]
        heat_source: SpaceHeatSystemHeatSource,
        #[serde(rename = "Control")]
        control: Option<String>,
    },
}

impl SpaceHeatSystemDetails {
    pub fn set_control(&mut self, control_string: impl Into<String>) -> anyhow::Result<&Self> {
        match self {
            SpaceHeatSystemDetails::InstantElectricHeater {
                ref mut control, ..
            } => {
                *control = Some(control_string.into());
            }
            SpaceHeatSystemDetails::ElectricStorageHeater {
                ref mut control, ..
            } => {
                *control = Some(control_string.into());
            }
            SpaceHeatSystemDetails::WetDistribution {
                ref mut control, ..
            } => {
                *control = Some(control_string.into());
            }
            SpaceHeatSystemDetails::WarmAir {
                ref mut control, ..
            } => {
                *control = Some(control_string.into());
            }
        }
        Ok(self)
    }

    pub fn temp_setback(&self) -> Option<f64> {
        match self {
            SpaceHeatSystemDetails::InstantElectricHeater { temp_setback, .. } => *temp_setback,
            _ => None,
        }
    }

    pub fn advanced_start(&self) -> Option<f64> {
        match self {
            SpaceHeatSystemDetails::WetDistribution { advanced_start, .. } => *advanced_start,
            _ => None,
        }
    }
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct SpaceHeatSystemHeatSource {
    pub name: String,
    pub temp_flow_limit_upper: Option<f64>,
}

// it is unclear whether this struct should be used - see reference to the struct above
#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[allow(dead_code)]
#[serde(deny_unknown_fields)]
pub(crate) struct EcoDesignController {
    pub(crate) ecodesign_control_class: EcoDesignControllerClass,
    pub(crate) min_outdoor_temp: Option<f64>,
    pub(crate) max_outdoor_temp: Option<f64>,
    pub(crate) min_flow_temp: Option<f64>,
}

#[derive(Clone, Copy, Debug, Deserialize_repr, PartialEq)]
#[repr(u8)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum EcoDesignControllerClass {
    // exact possible classes tbc
    ClassI = 1,
    ClassII = 2,
    ClassIII = 3,
    ClassIV = 4,
    ClassV = 5,
    ClassVI = 6,
    ClassVII = 7,
    ClassVIII = 8,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "type", deny_unknown_fields)]
pub enum Ventilation {
    #[serde(rename = "NatVent")]
    Natural { req_ach: f64 },
    #[serde(rename = "WHEV")]
    Whev {
        req_ach: f64,
        #[serde(rename = "SFP")]
        sfp: f64,
        #[serde(rename = "EnergySupply")]
        energy_supply: EnergySupplyType,
    },
}

impl Ventilation {
    pub fn req_ach(&self) -> f64 {
        match self {
            Ventilation::Natural { req_ach, .. } => *req_ach,
            Ventilation::Whev { req_ach, .. } => *req_ach,
        }
    }
}

#[derive(Copy, Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "lowercase")]
pub enum MVHRLocation {
    Inside,
    Outside,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct Infiltration {
    pub storeys_in_building: u32,
    pub storey_of_dwelling: Option<u32>,
    pub shelter: InfiltrationShelterType,
    pub build_type: InfiltrationBuildType,
    pub test_result: f64,
    pub test_type: InfiltrationTestType,
    pub env_area: f64,
    pub volume: f64,
    pub sheltered_sides: u32,
    pub open_chimneys: u32,
    pub open_flues: u32,
    pub closed_fire: u32,
    pub flues_d: u32,
    pub flues_e: u32,
    pub blocked_chimneys: u32,
    pub extract_fans: u32,
    pub passive_vents: u32,
    pub gas_fires: u32,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "lowercase")]
pub enum InfiltrationShelterType {
    VerySheltered,
    Sheltered,
    Normal,
    Exposed,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "lowercase")]
pub enum InfiltrationBuildType {
    House,
    Flat,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum InfiltrationTestType {
    #[serde(rename = "50Pa")]
    FiftyPascals,
    #[serde(rename = "4Pa")]
    FourPascals,
}

pub type ZoneDictionary = IndexMap<String, ZoneInput>;

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct ZoneInput {
    #[serde(rename = "SpaceHeatSystem")]
    pub space_heat_system: Option<String>,
    #[serde(rename = "SpaceCoolSystem")]
    pub space_cool_system: Option<String>,
    #[serde(rename = "SpaceHeatControl")]
    pub space_heat_control: Option<SpaceHeatControlType>,
    // don't know what the options are yet
    #[serde(rename = "Control_WindowOpening")]
    pub control_window_opening: Option<HeatSourceControlType>,
    pub area: f64,
    pub volume: f64,
    // check upstream whether this is used
    #[serde(rename = "Lighting")]
    pub lighting: Option<ZoneLighting>,
    // check upstream whether these two are used
    #[serde(rename = "temp_setpnt_heat")]
    _temp_setpnt_heat: Option<f64>,
    #[serde(rename = "temp_setpnt_cool")]
    _temp_setpnt_cool: Option<f64>,
    pub temp_setpnt_init: Option<f64>,
    #[serde(rename = "BuildingElement")]
    pub building_elements: IndexMap<String, BuildingElement>,
    #[serde(rename = "ThermalBridging")]
    pub thermal_bridging: ThermalBridging,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum SpaceHeatControlType {
    #[serde(rename = "livingroom")]
    LivingRoom,
    #[serde(rename = "restofdwelling")]
    RestOfDwelling,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct ZoneLighting {
    #[serde(rename = "efficacy")]
    efficacy: f64,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "type", deny_unknown_fields)]
pub enum BuildingElement {
    #[serde(rename = "BuildingElementOpaque")]
    Opaque {
        is_unheated_pitched_roof: Option<bool>,
        a_sol: f64,
        u_value: Option<f64>,
        r_c: Option<f64>,
        k_m: f64,
        mass_distribution_class: MassDistributionClass,
        is_external_door: Option<bool>,
        pitch: f64,
        #[serde(rename = "orientation360")]
        #[serde(deserialize_with = "deserialize_orientation")]
        orientation: f64,
        base_height: f64,
        height: f64,
        width: f64,
        area: f64,
        h_ci: Option<f64>,
        h_ri: Option<f64>,
        h_ce: Option<f64>,
        h_re: Option<f64>,
    },
    #[serde(rename = "BuildingElementTransparent")]
    Transparent {
        u_value: Option<f64>,
        #[serde(rename = "area")]
        // area is sometimes present but not expected to be used
        _area: Option<f64>,
        #[serde(rename = "Control_WindowOpenable")]
        window_openable_control: Option<String>, // unclear how this might be used
        r_c: Option<f64>,
        pitch: f64,
        #[serde(rename = "orientation360")]
        #[serde(deserialize_with = "deserialize_orientation")]
        orientation: f64,
        g_value: f64,
        frame_area_fraction: f64,
        base_height: f64,
        height: f64,
        width: f64,
        // following attributes may be FHS only
        free_area_height: Option<f64>,
        mid_height: Option<f64>,
        max_window_open_area: Option<f64>,
        security_risk: Option<bool>,
        window_part_list: Option<Vec<WindowPart>>,
        // end possible FHS only attributes
        shading: Vec<WindowShadingObject>,
    },
    #[serde(rename = "BuildingElementGround")]
    Ground {
        area: f64,
        total_area: f64,
        pitch: f64,
        u_value: f64,
        r_f: f64,
        k_m: f64,
        mass_distribution_class: MassDistributionClass,
        floor_type: FloorType,
        height_upper_surface: Option<f64>,
        #[serde(rename = "thermal_transm_walls")]
        thermal_transmission_walls: Option<f64>,
        #[serde(rename = "thermal_resist_insul")]
        thermal_resistance_of_insulation: Option<f64>,
        area_per_perimeter_vent: Option<f64>,
        shield_fact_location: Option<WindShieldLocation>,
        thickness_walls: f64,
        depth_basement_floor: Option<f64>,
        #[serde(rename = "thermal_resist_walls_base")]
        thermal_resistance_of_basement_walls: Option<f64>,
        #[serde(rename = "thermal_transm_envi_base")]
        thermal_transmittance_of_floor_above_basement: Option<f64>,
        height_basement_walls: Option<f64>,
        perimeter: f64,
        psi_wall_floor_junc: f64,
        edge_insulation: Option<Vec<EdgeInsulation>>,
    },
    #[serde(rename = "BuildingElementAdjacentZTC")]
    AdjacentZTC {
        area: f64,
        pitch: f64,
        u_value: Option<f64>,
        r_c: Option<f64>,
        k_m: f64,
        mass_distribution_class: MassDistributionClass,
    },
    #[serde(rename = "BuildingElementAdjacentZTU_Simple")]
    AdjacentZTUSimple {
        area: f64,
        pitch: f64,
        u_value: Option<f64>,
        r_c: Option<f64>,
        r_u: f64,
        k_m: f64,
        mass_distribution_class: MassDistributionClass,
    },
}

impl BuildingElement {
    pub(crate) fn pitch(&self) -> f64 {
        *match self {
            BuildingElement::Opaque { pitch, .. } => pitch,
            BuildingElement::Transparent { pitch, .. } => pitch,
            BuildingElement::Ground { pitch, .. } => pitch,
            BuildingElement::AdjacentZTC { pitch, .. } => pitch,
            BuildingElement::AdjacentZTUSimple { pitch, .. } => pitch,
        }
    }
}

pub(crate) trait TransparentBuildingElement {
    fn set_window_openable_control(&mut self, control: &str);
    fn is_security_risk(&self) -> bool;
}

impl TransparentBuildingElement for BuildingElement {
    fn set_window_openable_control(&mut self, control: &str) {
        match self {
            BuildingElement::Transparent {
                ref mut window_openable_control,
                ..
            } => {
                *window_openable_control = Some(control.to_owned());
            }
            _ => unreachable!(),
        }
    }

    fn is_security_risk(&self) -> bool {
        match self {
            BuildingElement::Transparent { security_risk, .. } => security_risk.unwrap_or(false),
            _ => unreachable!(),
        }
    }
}

// special deserialization logic so that orientations are normalized correctly on the way in
pub fn deserialize_orientation<'de, D>(deserializer: D) -> Result<f64, D::Error>
where
    D: Deserializer<'de>,
{
    let orientation360_value: f64 = Deserialize::deserialize(deserializer)?;
    Ok(init_orientation(orientation360_value))
}

pub(crate) fn init_orientation(value: f64) -> f64 {
    // Convert orientation from 0-360 (clockwise) to -180 to +180 (anticlockwise)
    180. - value
}

#[derive(Copy, Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum MassDistributionClass {
    D,
    E,
    I,
    IE,
    M,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct WindowPart {
    pub mid_height_air_flow_path: f64,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum FloorType {
    #[serde(rename = "Slab_no_edge_insulation")]
    SlabNoEdgeInsulation,
    #[serde(rename = "Slab_edge_insulation")]
    SlabEdgeInsulation,
    #[serde(rename = "Suspended_floor")]
    SuspendedFloor,
    #[serde(rename = "Heated_basement")]
    HeatedBasement,
    #[serde(rename = "Unheated_basement")]
    UnheatedBasement,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum WindShieldLocation {
    Sheltered,
    Average,
    Exposed,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(untagged)]
pub enum ThermalBridging {
    ThermalBridgingElements(IndexMap<String, ThermalBridgingDetails>),
    ThermalBridgingNumber(f64),
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "type", deny_unknown_fields)]
pub enum ThermalBridgingDetails {
    #[serde(rename(deserialize = "ThermalBridgeLinear"))]
    Linear {
        linear_thermal_transmittance: f64,
        length: f64,
        #[allow(dead_code)]
        junction_type: Option<String>,
    },
    #[serde(rename = "ThermalBridgePoint")]
    Point {
        #[serde(rename = "heat_transfer_coeff")]
        heat_transfer_coefficient: f64,
    },
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum HeatingControlType {
    #[serde(rename = "SeparateTimeAndTempControl")]
    SeparateTimeAndTemperatureControl,
    #[serde(rename = "SeparateTempControl")]
    SeparateTemperatureControl,
}

pub type SpaceCoolSystem = IndexMap<String, SpaceCoolSystemDetails>;

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct SpaceCoolSystemDetails {
    #[serde(rename(deserialize = "type"))]
    pub system_type: SpaceCoolSystemType,
    temp_setback: Option<f64>,
    pub cooling_capacity: f64,
    pub efficiency: f64,
    pub frac_convective: f64,
    #[serde(rename = "EnergySupply")]
    pub energy_supply: EnergySupplyType,
    #[serde(rename = "Control")]
    pub control: Option<String>,
}

impl SpaceCoolSystemDetails {
    pub fn set_control(&mut self, control_string: impl Into<String>) -> anyhow::Result<&Self> {
        self.control = Some(control_string.into());
        Ok(self)
    }

    pub fn temp_setback(&self) -> Option<f64> {
        self.temp_setback
    }
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum SpaceCoolSystemType {
    AirConditioning,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum WaterHeatingSchedule {
    AllDay,
    HeatingHours,
}

pub type HeatSourceWet = IndexMap<String, HeatSourceWetDetails>;

#[derive(Clone, Debug, Deserialize, Validate)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[allow(clippy::large_enum_variant)]
#[serde(tag = "type", deny_unknown_fields)]
pub enum HeatSourceWetDetails {
    HeatPump {
        #[serde(rename = "EnergySupply")]
        energy_supply: EnergySupplyType,
        source_type: HeatPumpSourceType,
        #[serde(rename = "EnergySupply_heat_network")]
        energy_supply_heat_network: Option<String>,
        // unclear what this is
        temp_distribution_heat_network: Option<f64>,
        sink_type: HeatPumpSinkType,
        #[serde(rename = "backup_ctrl_type")]
        backup_control_type: HeatPumpBackupControlType,
        time_delay_backup: f64,
        modulating_control: bool,
        min_modulation_rate_20: Option<f64>,
        min_modulation_rate_35: Option<f64>,
        min_modulation_rate_55: Option<f64>,
        time_constant_onoff_operation: f64,
        temp_return_feed_max: f64,
        temp_lower_operating_limit: f64,
        min_temp_diff_flow_return_for_hp_to_operate: f64,
        var_flow_temp_ctrl_during_test: bool,
        power_heating_warm_air_fan: Option<f64>,
        power_heating_circ_pump: f64,
        power_source_circ_pump: f64,
        power_standby: f64,
        power_crankcase_heater: f64,
        power_off: f64,
        power_max_backup: Option<f64>,
        eahp_mixed_max_temp: Option<f64>,
        eahp_mixed_min_temp: Option<f64>,
        #[serde(rename = "MechanicalVentilation")]
        mechanical_ventilation: Option<String>,
        #[serde(rename = "BufferTank")]
        buffer_tank: Option<Box<HeatPumpBufferTank>>,
        test_data: Vec<HeatPumpTestDatum>,
        boiler: Option<Box<HeatPumpBoiler>>,
    },
    Boiler {
        #[serde(rename = "EnergySupply")]
        energy_supply: EnergySupplyType,
        #[serde(rename = "EnergySupply_aux")]
        energy_supply_auxiliary: EnergySupplyType,
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
        energy_supply: EnergySupplyType,
        heat_battery_location: HeatSourceLocation,
        electricity_circ_pump: f64,
        electricity_standby: f64,
        // in kW (Charging)
        rated_charge_power: f64,
        // in kWh
        heat_storage_capacity: f64,
        // in kW (Output to hot water and space heat services)
        max_rated_heat_output: f64,
        // in kW (Losses to internal or external)
        max_rated_losses: f64,
        // number of units installed in zone
        number_of_units: usize,
        #[serde(rename = "ControlCharge")]
        control_charge: String,
    },
    #[serde(rename = "HIU")]
    Hiu {
        #[serde(rename = "EnergySupply")]
        energy_supply: EnergySupplyType,
        power_max: f64,
        #[serde(rename = "HIU_daily_loss")]
        hiu_daily_loss: f64,
        building_level_distribution_losses: f64,
    },
}

impl From<&HeatPumpBoiler> for HeatSourceWetDetails {
    fn from(value: &HeatPumpBoiler) -> Self {
        HeatSourceWetDetails::Boiler {
            energy_supply: value.energy_supply,
            energy_supply_auxiliary: value.energy_supply_auxiliary,
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

#[derive(Copy, Clone, Debug, Deserialize, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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

#[derive(Copy, Clone, Debug, Deserialize, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum HeatPumpSinkType {
    Water,
    Air,
}

#[derive(Copy, Clone, Debug, Deserialize, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum HeatPumpBackupControlType {
    None,
    TopUp,
    Substitute,
}

#[derive(Copy, Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct HeatPumpBufferTank {
    pub daily_losses: f64,
    pub volume: f64,
    pub pump_fixed_flow_rate: f64,
    pub pump_power_at_flow_rate: f64,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct HeatPumpTestDatum {
    pub air_flow_rate: Option<f64>,
    pub(crate) test_letter: TestLetter,
    pub capacity: f64,
    pub cop: f64,
    #[serde(rename = "degradation_coeff")]
    pub degradation_coefficient: f64,
    pub design_flow_temp: f64,
    pub temp_outlet: f64,
    pub temp_source: f64,
    pub temp_test: f64,
    pub eahp_mixed_ext_air_ratio: Option<f64>,
    #[serde(skip)]
    pub ext_air_ratio: Option<f64>,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct HeatPumpBoiler {
    #[serde(rename = "EnergySupply")]
    pub(crate) energy_supply: EnergySupplyType,
    #[serde(rename = "EnergySupply_aux")]
    pub(crate) energy_supply_auxiliary: EnergySupplyType,
    rated_power: f64,
    efficiency_full_load: f64,
    efficiency_part_load: f64,
    boiler_location: HeatSourceLocation,
    modulation_load: f64,
    electricity_circ_pump: f64,
    electricity_part_load: f64,
    electricity_full_load: f64,
    electricity_standby: f64,
    pub(crate) cost_schedule_hybrid: Option<BoilerCostScheduleHybrid>,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub(crate) struct BoilerCostScheduleHybrid {
    pub cost_schedule_start_day: u32,
    pub cost_schedule_time_series_step: f64,
    pub cost_schedule_hp: NumericSchedule,
    pub cost_schedule_boiler: NumericSchedule,
}

#[derive(Copy, Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum HeatSourceLocation {
    #[serde(rename = "internal")]
    Internal,
    #[serde(rename = "external")]
    External,
}

pub type WasteWaterHeatRecovery = IndexMap<String, WasteWaterHeatRecoveryDetails>;

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct WasteWaterHeatRecoveryDetails {
    #[serde(rename = "type")]
    pub system_type: WwhrsType,
    #[serde(rename = "ColdWaterSource")]
    pub cold_water_source: ColdWaterSourceType,
    pub flow_rates: Vec<f64>,
    pub efficiencies: Vec<f64>,
    pub utilisation_factor: f64,
    pub electrical_consumption: Option<f64>,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum WwhrsType {
    #[serde(rename = "WWHRS_InstantaneousSystemA")]
    SystemA,
    #[serde(rename = "WWHRS_InstantaneousSystemB")]
    SystemB,
    #[serde(rename = "WWHRS_InstantaneousSystemC")]
    SystemC,
}

pub type OnSiteGeneration = IndexMap<String, OnSiteGenerationDetails>;

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(tag = "type", deny_unknown_fields)]
pub enum OnSiteGenerationDetails {
    PhotovoltaicSystem {
        peak_power: f64,
        ventilation_strategy: OnSiteGenerationVentilationStrategy,
        pitch: f64,
        #[serde(rename = "orientation360")]
        #[serde(deserialize_with = "deserialize_orientation")]
        orientation: f64,
        base_height: f64,
        height: f64,
        width: f64,
        #[serde(rename = "EnergySupply")]
        energy_supply: EnergySupplyType,
        shading: Vec<WindowShadingObject>,
        inverter_peak_power: f64,
        inverter_is_inside: bool,
    },
}

#[derive(Copy, Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum OnSiteGenerationVentilationStrategy {
    #[serde(rename = "unventilated")]
    Unventilated,
    #[serde(rename = "moderately_ventilated")]
    ModeratelyVentilated,
    #[serde(rename = "strongly_or_forced_ventilated")]
    StronglyOrForcedVentilated,
    #[serde(rename = "rear_surface_free")]
    RearSurfaceFree,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct WindowOpeningForCooling {
    pub equivalent_area: f64,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct General {
    pub storeys_in_building: usize,
    pub storey_of_dwelling: Option<usize>,
    pub build_type: BuildType,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum BuildType {
    #[serde(rename = "house")]
    House,
    #[serde(rename = "flat")]
    Flat,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct InfiltrationVentilation {
    pub(crate) cross_vent_factor: bool,
    pub(crate) shield_class: VentilationShieldClass,
    pub(crate) terrain_class: TerrainClass,
    pub(crate) altitude: f64,
    #[serde(rename = "Control_WindowAdjust")]
    pub(crate) window_adjust_control: Option<String>, // don't know what this can be
    noise_nuisance: Option<bool>,
    #[serde(rename = "Vents")]
    pub(crate) vents: IndexMap<String, Vent>,
    #[serde(rename = "Leaks")]
    pub(crate) leaks: VentilationLeaks,
    #[serde(rename = "MechanicalVentilation", default)]
    pub(crate) mechanical_ventilation: IndexMap<String, MechanicalVentilation>,
    #[serde(rename = "AirTerminalDevices")]
    _air_terminal_devices: Option<IndexMap<String, AirTerminalDevice>>,
    #[serde(rename = "PDUs")]
    _pdus: IndexMap<String, ()>, // don't know what this looks like yet
    #[serde(rename = "Cowls")]
    _cowls: IndexMap<String, ()>, // don't know what this looks like yet
    #[serde(rename = "CombustionAppliances")]
    pub(crate) combustion_appliances: IndexMap<String, CombustionAppliance>,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum VentilationShieldClass {
    Open,
    Normal,
    Shielded,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum TerrainClass {
    #[serde(rename = "Open terrain")]
    OpenTerrain,
    Country,
    Urban,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct Vent {
    pub(crate) mid_height_air_flow_path: f64,
    pub(crate) area_cm2: f64,
    pub(crate) pressure_difference_ref: f64,
    #[serde(rename = "orientation360")]
    #[serde(deserialize_with = "deserialize_orientation")]
    pub(crate) orientation: f64,
    pub(crate) pitch: f64,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct VentilationLeaks {
    pub ventilation_zone_height: f64,
    pub test_pressure: f64,
    pub test_result: f64,
    pub env_area: f64,
    // following values appear to be usually overridden
    #[serde(rename = "area_roof")]
    _area_roof: Option<f64>,
    #[serde(rename = "area_facades")]
    _area_facades: Option<f64>,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct MechanicalVentilation {
    #[serde(rename = "vent_sys_op")]
    _vent_sys_op: Option<String>, // this seems useless/unreferenced
    #[serde(rename = "sup_air_flw_ctrl")]
    pub(crate) supply_air_flow_rate_control: SupplyAirFlowRateControlType,
    #[serde(rename = "sup_air_temp_ctrl")]
    pub(crate) supply_air_temperature_control_type: SupplyAirTemperatureControlType,
    #[serde(rename = "design_zone_cooling_covered_by_mech_vent")]
    _design_zone_cooling_covered_by_mechanical_vent: Option<f64>,
    #[serde(rename = "design_zone_heating_covered_by_mech_vent")]
    _design_zone_heating_covered_by_mechanical_vent: Option<f64>,
    pub(crate) vent_type: VentType,
    #[serde(rename = "mvhr_eff")]
    pub(crate) mvhr_efficiency: Option<f64>,
    pub(crate) mvhr_location: Option<MVHRLocation>,
    #[serde(rename = "Control")]
    pub(crate) control: Option<String>,
    #[serde(rename = "SFP")]
    pub(crate) sfp: Option<f64>,
    #[serde(rename = "measured_fan_power")]
    pub(crate) measured_fan_power: Option<f64>,
    #[serde(rename = "measured_air_flow_rate")]
    pub(crate) measured_air_flow_rate: Option<f64>,
    #[serde(rename = "EnergySupply")]
    pub(crate) energy_supply: EnergySupplyType,
    pub(crate) design_outdoor_air_flow_rate: f64,
    pub(crate) ductwork: Option<Vec<MechanicalVentilationDuctwork>>,
}

pub(crate) trait MechanicalVentilationForProcessing {
    fn vent_type(&self) -> VentType;
    fn measured_fan_power(&self) -> Option<f64>;
    fn measured_air_flow_rate(&self) -> Option<f64>;
    fn set_sfp(&mut self, sfp: f64);
    fn set_control(&mut self, control: &str);
}

impl MechanicalVentilationForProcessing for MechanicalVentilation {
    fn set_sfp(&mut self, sfp: f64) {
        self.sfp = Some(sfp);
    }

    fn vent_type(&self) -> VentType {
        self.vent_type
    }

    fn measured_fan_power(&self) -> Option<f64> {
        self.measured_fan_power
    }

    fn measured_air_flow_rate(&self) -> Option<f64> {
        self.measured_air_flow_rate
    }

    fn set_control(&mut self, control: &str) {
        self.control = Some(control.to_owned());
    }
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum SupplyAirFlowRateControlType {
    ODA,
    #[serde(rename = "LOAD")]
    Load,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum SupplyAirTemperatureControlType {
    #[serde(rename = "CONST")]
    Constant,
    #[serde(rename = "NO_CTRL")]
    NoControl,
    #[serde(rename = "LOAD_COM")]
    LoadCom,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum VentType {
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

impl Display for VentType {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            serde_json::to_value(self).unwrap().as_str().unwrap()
        )
    }
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct MechanicalVentilationDuctwork {
    pub(crate) cross_section_shape: DuctShape,
    pub(crate) duct_perimeter_mm: Option<f64>,
    pub(crate) internal_diameter_mm: Option<f64>,
    pub(crate) external_diameter_mm: Option<f64>,
    pub(crate) length: f64,
    pub(crate) insulation_thermal_conductivity: f64,
    pub(crate) insulation_thickness_mm: f64,
    pub(crate) reflective: bool,
    pub(crate) duct_type: DuctType,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "snake_case")]
pub enum DuctShape {
    Circular,
    Rectangular,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "snake_case")]
pub enum DuctType {
    Intake,
    Supply,
    Extract,
    Exhaust,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
// type is not read in 0.30 version of code - please remove following once data is used
#[allow(dead_code)]
pub struct AirTerminalDevice {
    area_cm2: f64,
    pressure_difference_ref: f64,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct CombustionAppliance {
    pub(crate) supply_situation: CombustionAirSupplySituation,
    pub(crate) exhaust_situation: FlueGasExhaustSituation,
    pub(crate) fuel_type: CombustionFuelType,
    pub(crate) appliance_type: CombustionApplianceType,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum CombustionAirSupplySituation {
    #[serde(rename = "room_air")]
    RoomAir,
    #[serde(rename = "outside")]
    Outside,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum FlueGasExhaustSituation {
    #[serde(rename = "into_room")]
    IntoRoom,
    #[serde(rename = "into_separate_duct")]
    IntoSeparateDuct,
    #[serde(rename = "into_mech_vent")]
    IntoMechVent,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "snake_case")]
pub enum CombustionFuelType {
    Wood,
    Gas,
    Oil,
    Coal,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub(crate) enum ApplianceKey {
    Fridge,
    Freezer,
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
}

impl ApplianceKey {
    pub(crate) fn is_clothes_appliance(&self) -> bool {
        [Self::ClothesWashing, Self::ClothesDrying].contains(self)
    }
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
        appliance_key.to_string()
    }
}

impl TryFrom<&str> for ApplianceKey {
    type Error = anyhow::Error;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        serde_json::from_str(format!("\"{}\"", value).as_str()).map_err(|err| anyhow!(err))
    }
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(untagged)]
pub enum ApplianceEntry {
    Object(Appliance),
    Reference(ApplianceReference),
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct Appliance {
    #[serde(rename = "kWh_per_100cycle")]
    pub(crate) kwh_per_100_cycle: Option<f64>,
    #[serde(rename = "loadshifting")]
    pub(crate) load_shifting: Option<ApplianceLoadShifting>,
    pub(crate) kg_load: Option<f64>,
    #[serde(rename = "kWh_per_annum")]
    pub(crate) kwh_per_annum: Option<f64>,
    #[serde(rename = "Energysupply")]
    pub(crate) energy_supply: Option<EnergySupplyType>,
    #[serde(rename = "kWh_per_cycle")]
    pub(crate) kwh_per_cycle: Option<f64>,
    pub(crate) standard_use: Option<f64>,
}

impl Appliance {
    pub(crate) fn with_energy_supply(energy_supply: EnergySupplyType, kwh_per_cycle: f64) -> Self {
        Self {
            kwh_per_100_cycle: None,
            load_shifting: None,
            kg_load: None,
            kwh_per_annum: None,
            energy_supply: Some(energy_supply),
            kwh_per_cycle: Some(kwh_per_cycle),
            standard_use: None,
        }
    }

    pub(crate) fn with_kwh_per_cycle(kwh_per_cycle: f64) -> Self {
        Self {
            kwh_per_100_cycle: None,
            load_shifting: None,
            kg_load: None,
            kwh_per_annum: None,
            energy_supply: None,
            kwh_per_cycle: Some(kwh_per_cycle),
            standard_use: None,
        }
    }

    pub(crate) fn with_kwh_per_annum(kwh_per_annum: f64) -> Self {
        Self {
            kwh_per_100_cycle: None,
            load_shifting: None,
            kg_load: None,
            kwh_per_annum: Some(kwh_per_annum),
            energy_supply: None,
            kwh_per_cycle: None,
            standard_use: None,
        }
    }

    pub(crate) fn with_kwh_per_100_cycle(kwh_per_100_cycle: f64, kg_load: Option<f64>) -> Self {
        Self {
            kwh_per_100_cycle: Some(kwh_per_100_cycle),
            load_shifting: None,
            kg_load,
            kwh_per_annum: None,
            energy_supply: None,
            kwh_per_cycle: None,
            standard_use: None,
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct ApplianceLoadShifting {
    pub max_shift_hrs: f64,
    pub demand_limit_weighted: f64,
    #[serde(rename = "weight")]
    _weight: String, // not sure yet what these can be
    // In Python these are set from the FHS wrapper
    pub weight_timeseries: Option<Vec<f64>>,
    pub demand_timeseries: Option<Vec<f64>>,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum ApplianceReference {
    Default,
    #[serde(rename = "Not Installed")]
    NotInstalled,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields)]
pub struct Tariff {
    #[serde(rename = "schedule")]
    schedule: NumericSchedule,
}

#[derive(Debug)]
pub struct InputForProcessing {
    input: Input,
}

/// This type makes methods available for restricted access by wrappers,
/// in order to work towards a reasonable API for wrappers to interact with inputs rather than
/// the more brittle approach of allowing full access to the input data structure.
/// If the full access is encapsulated within methods here, it becomes possible to update the
/// underlying structure without breaking wrappers.
impl InputForProcessing {
    pub fn init_with_json(json: impl Read) -> Result<Self, anyhow::Error> {
        let reader = BufReader::new(json);

        let input: Input = serde_json::from_reader(reader)?;

        Ok(Self { input })
    }

    pub fn finalize(self) -> Input {
        self.input
    }

    pub fn set_simulation_time(&mut self, simulation_time: SimulationTime) -> &Self {
        self.input.simulation_time = simulation_time;
        self
    }

    pub(crate) fn merge_external_conditions_data(
        &mut self,
        external_conditions_data: Option<ExternalConditionsInput>,
    ) {
        if let Some(external_conditions) = external_conditions_data {
            let shading_segments = self.input.external_conditions.shading_segments.clone();
            let mut new_external_conditions: ExternalConditionsInput = external_conditions;
            new_external_conditions.shading_segments = shading_segments;
            self.input.external_conditions = Arc::new(new_external_conditions);
        }
    }

    pub fn reset_internal_gains(&mut self) -> &Self {
        self.input.internal_gains = Default::default();
        self
    }

    pub fn total_zone_area(&self) -> f64 {
        self.input.zone.values().map(|z| z.area).sum::<f64>()
    }

    pub(crate) fn total_zone_volume(&self) -> f64 {
        self.input.zone.values().map(|z| z.volume).sum::<f64>()
    }

    pub fn area_for_zone(&self, zone: &str) -> anyhow::Result<f64> {
        Ok(self
            .input
            .zone
            .get(zone)
            .ok_or(anyhow!("Used zone key for a zone that does not exist"))?
            .area)
    }

    pub fn number_of_bedrooms(&self) -> Option<usize> {
        self.input.number_of_bedrooms
    }

    pub fn set_metabolic_gains(
        &mut self,
        start_day: u32,
        time_series_step: f64,
        schedule_json: Value,
    ) -> anyhow::Result<&Self> {
        self.input.internal_gains.metabolic_gains = Some(InternalGainsDetails {
            start_day,
            time_series_step,
            schedule: serde_json::from_value(schedule_json)?,
        });

        Ok(self)
    }

    pub fn set_evaporative_losses(
        &mut self,
        start_day: u32,
        time_series_step: f64,
        schedule_json: Value,
    ) -> anyhow::Result<&Self> {
        self.input.internal_gains.evaporative_losses = Some(InternalGainsDetails {
            start_day,
            time_series_step,
            schedule: serde_json::from_value(schedule_json)?,
        });

        Ok(self)
    }

    pub(crate) fn set_cold_water_losses(
        &mut self,
        start_day: u32,
        time_series_step: f64,
        schedule_json: Value,
    ) -> anyhow::Result<&Self> {
        self.input.internal_gains.cold_water_losses = Some(InternalGainsDetails {
            start_day,
            time_series_step,
            schedule: serde_json::from_value(schedule_json)?,
        });

        Ok(self)
    }

    pub fn heating_control_type(&self) -> Option<HeatingControlType> {
        self.input.heating_control_type
    }

    pub fn set_heating_control_type(
        &mut self,
        heating_control_type_value: Value,
    ) -> anyhow::Result<&Self> {
        self.input.heating_control_type = Some(serde_json::from_value(heating_control_type_value)?);
        Ok(self)
    }

    pub fn add_control(
        &mut self,
        control_key: impl Into<String>,
        control_json: Value,
    ) -> anyhow::Result<&Self> {
        self.input
            .control
            .extra
            .insert(control_key.into(), serde_json::from_value(control_json)?);

        Ok(self)
    }

    pub fn zone_keys(&self) -> Vec<String> {
        self.input.zone.keys().cloned().collect()
    }

    pub fn set_init_temp_setpoint_for_zone(
        &mut self,
        zone: &str,
        temperature: f64,
    ) -> anyhow::Result<&Self> {
        self.input
            .zone
            .get_mut(zone)
            .ok_or(anyhow!("Used zone key for a zone that does not exist"))?
            .temp_setpnt_init = Some(temperature);
        Ok(self)
    }

    pub fn space_heat_control_for_zone(
        &self,
        zone: &str,
    ) -> anyhow::Result<Option<SpaceHeatControlType>> {
        Ok(self
            .input
            .zone
            .get(zone)
            .ok_or(anyhow!("Used zone key for a zone that does not exist"))?
            .space_heat_control)
    }

    pub fn space_heat_system_for_zone(&self, zone: &str) -> anyhow::Result<Option<String>> {
        Ok(self
            .input
            .zone
            .get(zone)
            .ok_or(anyhow!("Used zone key for a zone that does not exist"))?
            .space_heat_system
            .clone())
    }

    pub fn set_space_heat_system_for_zone(
        &mut self,
        zone: &str,
        system_name: &str,
    ) -> anyhow::Result<&Self> {
        let zone = self
            .input
            .zone
            .get_mut(zone)
            .ok_or(anyhow!("Used zone key for a zone that does not exist"))?;
        zone.space_heat_system = Some(system_name.to_string());
        Ok(self)
    }

    pub fn space_cool_system_for_zone(&self, zone: &str) -> anyhow::Result<Option<String>> {
        Ok(self
            .input
            .zone
            .get(zone)
            .ok_or(anyhow!("Used zone key for a zone that does not exist"))?
            .space_cool_system
            .clone())
    }

    pub fn set_space_cool_system_for_zone(
        &mut self,
        zone: &str,
        system_name: &str,
    ) -> anyhow::Result<&Self> {
        let zone = self
            .input
            .zone
            .get_mut(zone)
            .ok_or(anyhow!("Used zone key for a zone that does not exist"))?;
        zone.space_cool_system = Some(system_name.to_string());
        Ok(self)
    }

    pub fn lighting_efficacy_for_zone(&self, zone: &str) -> anyhow::Result<Option<f64>> {
        Ok(self
            .input
            .zone
            .get(zone)
            .ok_or(anyhow!("Used zone key for a zone that does not exist"))?
            .lighting
            .as_ref()
            .map(|lighting| lighting.efficacy))
    }

    pub fn set_lighting_efficacy_for_all_zones(&mut self, efficacy: f64) -> &Self {
        for lighting in self
            .input
            .zone
            .values_mut()
            .filter_map(|zone| zone.lighting.as_mut())
        {
            lighting.efficacy = efficacy;
        }
        self
    }

    pub fn set_control_window_opening_for_zone(
        &mut self,
        zone: &str,
        opening_type: Option<HeatSourceControlType>,
    ) -> anyhow::Result<&Self> {
        self.input
            .zone
            .get_mut(zone)
            .ok_or(anyhow!("Used zone key for a zone that does not exist"))?
            .control_window_opening = opening_type;
        Ok(self)
    }

    pub fn set_control_string_for_space_heat_system(
        &mut self,
        space_heat_system: &str,
        control_string: &str,
    ) -> anyhow::Result<&Self> {
        self.input
            .space_heat_system
            .as_mut()
            .ok_or(anyhow!(
                "There is no space heat system provided at the root of the input"
            ))?
            .get_mut(space_heat_system)
            .ok_or(anyhow!(
                "There is no provided space heat system with the name '{space_heat_system}'"
            ))?
            .set_control(control_string)?;
        Ok(self)
    }

    pub fn set_control_string_for_space_cool_system(
        &mut self,
        space_cool_system: &str,
        control_string: &str,
    ) -> anyhow::Result<&Self> {
        self.input
            .space_cool_system
            .as_mut()
            .ok_or(anyhow!(
                "There is no space cool system provided at the root of the input"
            ))?
            .get_mut(space_cool_system)
            .ok_or(anyhow!(
                "There is no provided space cool system with the name '{space_cool_system}'"
            ))?
            .set_control(control_string)?;
        Ok(self)
    }

    pub fn temperature_setback_for_space_heat_system(
        &self,
        space_heat_system: &str,
    ) -> anyhow::Result<Option<f64>> {
        Ok(self
            .input
            .space_heat_system
            .as_ref()
            .ok_or(anyhow!(
                "There is no space heat system provided at the root of the input"
            ))?
            .get(space_heat_system)
            .ok_or(anyhow!(
                "There is no provided space heat system with the name '{space_heat_system}'"
            ))?
            .temp_setback())
    }

    pub fn temperature_setback_for_space_cool_system(
        &self,
        space_cool_system: &str,
    ) -> anyhow::Result<Option<f64>> {
        Ok(self
            .input
            .space_cool_system
            .as_ref()
            .ok_or(anyhow!(
                "There is no space cool system provided at the root of the input"
            ))?
            .get(space_cool_system)
            .ok_or(anyhow!(
                "There is no provided space heat system with the name '{space_cool_system}'"
            ))?
            .temp_setback())
    }

    pub fn advanced_start_for_space_heat_system(
        &self,
        space_heat_system: &str,
    ) -> anyhow::Result<Option<f64>> {
        Ok(self
            .input
            .space_heat_system
            .as_ref()
            .ok_or(anyhow!(
                "There is no space heat system provided at the root of the input"
            ))?
            .get(space_heat_system)
            .ok_or(anyhow!(
                "There is no provided space heat system with the name '{space_heat_system}'"
            ))?
            .advanced_start())
    }

    pub fn hot_water_source_keys(&self) -> Vec<String> {
        self.input.hot_water_source.source_keys()
    }

    pub fn hot_water_source_details_for_key(
        &mut self,
        source_key: &str,
    ) -> &mut impl HotWaterSourceDetailsForProcessing {
        self.input
            .hot_water_source
            .hot_water_source_for_processing(source_key)
    }

    pub fn set_lighting_gains(&mut self, gains_details: Value) -> anyhow::Result<&Self> {
        self.set_gains_for_field("lighting", gains_details)
    }

    pub fn set_gains_for_field(
        &mut self,
        field: impl Into<String>,
        gains_details: Value,
    ) -> anyhow::Result<&Self> {
        let gains_details: ApplianceGainsDetails = serde_json::from_value(gains_details)?;
        self.input
            .appliance_gains
            .insert(field.into(), gains_details);

        Ok(self)
    }

    pub fn energy_supply_type_for_appliance_gains_field(&self, field: &str) -> Option<String> {
        self.input
            .appliance_gains
            .get(field)
            .map(|details| details.energy_supply.to_string())
    }

    pub fn reset_appliance_gains_field(&mut self, field: &str) -> anyhow::Result<()> {
        self.input.appliance_gains.shift_remove(field);

        Ok(())
    }

    pub fn clear_appliance_gains(&mut self) {
        self.input.appliance_gains.clear();
    }

    pub fn fuel_type_for_energy_supply_field(&self, field: &str) -> anyhow::Result<String> {
        Ok(self
            .input
            .energy_supply
            .get(&(field.parse::<EnergySupplyKey>()?))
            .ok_or(anyhow!(
                "Fuel type not provided for energy supply field '{field}'"
            ))?
            .fuel
            .to_string())
    }

    pub fn shower_flowrates(&self) -> IndexMap<String, f64> {
        self.input
            .hot_water_demand
            .shower
            .as_ref()
            .map(|showers| {
                showers
                    .0
                    .iter()
                    .filter_map(|(name, shower)| match shower {
                        Shower::MixerShower { flowrate, .. } => Some((name.clone(), *flowrate)),
                        _ => None,
                    })
                    .collect()
            })
            .unwrap_or_default()
    }

    pub fn reset_water_heating_events(&mut self) {
        self.input.water_heating_events = Default::default();
    }

    pub fn shower_keys(&self) -> Vec<String> {
        match self.input.hot_water_demand.shower.as_ref() {
            Some(shower) => shower.keys(),
            None => Default::default(),
        }
    }

    pub fn shower_name_refers_to_instant_electric(&self, name: &str) -> bool {
        match &self.input.hot_water_demand.shower {
            Some(shower) => shower.name_refers_to_instant_electric_shower(name),
            None => false,
        }
    }

    pub fn bath_keys(&self) -> Vec<String> {
        match self.input.hot_water_demand.bath.as_ref() {
            Some(bath) => bath.keys(),
            None => Default::default(),
        }
    }

    pub fn size_for_bath_field(&self, field: &str) -> Option<f64> {
        self.input
            .hot_water_demand
            .bath
            .as_ref()
            .and_then(|bath| bath.size_for_field(field))
    }

    pub fn flowrate_for_bath_field(&self, field: &str) -> Option<f64> {
        self.input
            .hot_water_demand
            .bath
            .as_ref()
            .and_then(|bath| bath.flowrate_for_field(field))
    }

    pub fn other_water_use_keys(&self) -> Vec<String> {
        match self.input.hot_water_demand.other_water_use.as_ref() {
            Some(other) => other.keys(),
            None => Default::default(),
        }
    }

    pub fn flow_rate_for_other_water_use_field(&self, field: &str) -> Option<f64> {
        self.input
            .hot_water_demand
            .other_water_use
            .as_ref()
            .and_then(|other| other.flowrate_for_field(field))
    }

    pub fn set_other_water_use_details(
        &mut self,
        cold_water_source_type: impl Into<String>,
        flowrate: f64,
    ) -> anyhow::Result<()> {
        let other_details = OtherWaterUseDetails {
            flowrate,
            cold_water_source: serde_json::from_str(cold_water_source_type.into().as_str())?,
        };

        match self.input.hot_water_demand.other_water_use {
            Some(ref mut other_water_use) => {
                other_water_use.0.insert("other".to_string(), other_details);
            }
            None => {
                self.input.hot_water_demand.other_water_use = Some(OtherWaterUses(IndexMap::from(
                    [("other".to_string(), other_details)],
                )));
            }
        }

        Ok(())
    }

    /// Override all the vol_hw_daily_average values on the heat pump hot water only heat sources.
    pub fn override_vol_hw_daily_average_on_heat_pumps(&mut self, vol_hw_daily_average: f64) {
        if let HotWaterSourceDetails::StorageTank {
            ref mut heat_source,
            ..
        } = self.input.hot_water_source.hot_water_cylinder
        {
            let heat_source_keys: Vec<String> = heat_source.keys().cloned().collect();
            for heat_source_key in heat_source_keys {
                if let HeatSource::HeatPumpHotWaterOnly {
                    vol_hw_daily_average: ref mut input_hw_daily_average,
                    ..
                } = heat_source.get_mut(&heat_source_key).unwrap()
                {
                    *input_hw_daily_average = vol_hw_daily_average;
                }
            }
        }
    }

    pub fn part_g_compliance(&self) -> Option<bool> {
        self.input.part_g_compliance
    }

    pub fn set_part_g_compliance(&mut self, is_compliant: bool) -> &Self {
        self.input.part_g_compliance = Some(is_compliant);
        self
    }

    pub fn add_water_heating_event(
        &mut self,
        event_type: WaterHeatingEventType,
        subtype_name: &str,
        event: WaterHeatingEvent,
    ) {
        self.input.water_heating_events.add_event_for_type_and_name(
            event_type,
            subtype_name,
            event,
        );
    }

    pub(crate) fn water_heating_events_of_types(
        &self,
        event_types: &[WaterHeatingEventType],
    ) -> Vec<WaterHeatingEvent> {
        self.input
            .water_heating_events
            .water_heating_events_of_types(event_types)
    }

    pub fn defines_window_opening_for_cooling(&self) -> bool {
        self.input.window_opening_for_cooling.is_some()
    }

    pub fn cold_water_source_has_header_tank(&self) -> bool {
        self.input.cold_water_source.has_header_tank()
    }

    pub fn set_cold_water_source_by_key(
        &mut self,
        key: &str,
        source_details: Value,
    ) -> anyhow::Result<&Self> {
        self.input
            .cold_water_source
            .set_cold_water_source_details_by_key(key, serde_json::from_value(source_details)?);
        Ok(self)
    }

    pub fn set_hot_water_cylinder(&mut self, source_value: Value) -> anyhow::Result<&Self> {
        self.input.hot_water_source.hot_water_cylinder = serde_json::from_value(source_value)?;
        Ok(self)
    }

    pub fn set_water_distribution(&mut self, distribution_value: Value) -> anyhow::Result<&Self> {
        self.input.hot_water_demand.water_distribution =
            Some(serde_json::from_value(distribution_value)?);
        Ok(self)
    }

    pub fn set_shower(&mut self, shower_value: Value) -> anyhow::Result<&Self> {
        self.input.hot_water_demand.shower = Some(serde_json::from_value(shower_value)?);
        Ok(self)
    }

    pub fn set_bath(&mut self, bath_value: Value) -> anyhow::Result<&Self> {
        self.input.hot_water_demand.bath = Some(serde_json::from_value(bath_value)?);
        Ok(self)
    }

    pub fn set_other_water_use(&mut self, other_water_use_value: Value) -> anyhow::Result<&Self> {
        self.input.hot_water_demand.other_water_use =
            Some(serde_json::from_value(other_water_use_value)?);
        Ok(self)
    }

    pub fn remove_wwhrs(&mut self) -> &Self {
        self.input.waste_water_heat_recovery = None;
        self
    }

    pub fn remove_space_heat_systems(&mut self) -> &Self {
        self.input.space_heat_system = None;
        self
    }

    pub fn set_space_heat_system_for_key(
        &mut self,
        key: &str,
        space_heat_system_value: Value,
    ) -> anyhow::Result<&Self> {
        let empty_details: IndexMap<String, SpaceHeatSystemDetails> = Default::default();
        let system_details: SpaceHeatSystemDetails =
            serde_json::from_value(space_heat_system_value)?;
        let systems = self.input.space_heat_system.get_or_insert(empty_details);
        systems.insert(key.to_string(), system_details);
        Ok(self)
    }

    pub fn remove_space_cool_systems(&mut self) -> &Self {
        self.input.space_cool_system = None;
        self
    }

    pub fn set_space_cool_system_for_key(
        &mut self,
        key: &str,
        space_cool_system_value: Value,
    ) -> anyhow::Result<&Self> {
        let empty_details: IndexMap<String, SpaceCoolSystemDetails> = Default::default();
        let system_details: SpaceCoolSystemDetails =
            serde_json::from_value(space_cool_system_value)?;
        let systems = self.input.space_cool_system.get_or_insert(empty_details);
        systems.insert(key.to_string(), system_details);
        Ok(self)
    }

    pub fn remove_on_site_generation(&mut self) -> &mut Self {
        self.input.on_site_generation = None;
        self
    }

    pub fn remove_all_diverters_from_energy_supplies(&mut self) -> &mut Self {
        for energy_supply in self.input.energy_supply.values_mut() {
            energy_supply.diverter = None;
        }
        self
    }

    pub fn remove_all_batteries_from_energy_supplies(&mut self) -> &mut Self {
        for energy_supply in self.input.energy_supply.values_mut() {
            energy_supply.electric_battery = None;
        }
        self
    }

    pub(crate) fn external_conditions(&self) -> Arc<ExternalConditionsInput> {
        self.input.external_conditions.clone()
    }

    pub(crate) fn all_transparent_building_elements(&self) -> Vec<&BuildingElement> {
        self.input
            .zone
            .values()
            .flat_map(|zone| zone.building_elements.values())
            .filter(|el| matches!(el, BuildingElement::Transparent { .. }))
            .collect()
    }

    pub(crate) fn all_transparent_building_elements_mut(
        &mut self,
    ) -> Vec<&mut impl TransparentBuildingElement> {
        self.input
            .zone
            .values_mut()
            .flat_map(|zone| zone.building_elements.values_mut())
            .filter(|el| matches!(el, BuildingElement::Transparent { .. }))
            .collect()
    }

    pub(crate) fn all_energy_supply_fuel_types(&self) -> HashSet<FuelType> {
        let mut fuel_types: HashSet<FuelType> = Default::default();
        for energy_supply in self.input.energy_supply.values() {
            fuel_types.insert(energy_supply.fuel);
        }

        fuel_types
    }

    pub(crate) fn has_appliances(&self) -> bool {
        self.input.appliances.is_some()
    }

    pub(crate) fn merge_in_appliances(&mut self, appliances: &IndexMap<ApplianceKey, Appliance>) {
        let mut appliances: IndexMap<ApplianceKey, ApplianceEntry> = appliances
            .iter()
            .map(|(k, v)| (*k, ApplianceEntry::Object(v.clone())))
            .collect();
        if let Some(ref mut existing_appliances) = self.input.appliances.as_mut() {
            existing_appliances.append(&mut appliances);
        } else {
            self.input.appliances = Some(appliances);
        }
    }

    pub(crate) fn remove_appliance(&mut self, appliance_key: &ApplianceKey) {
        if let Some(ref mut appliances) = self.input.appliances.as_mut() {
            appliances.shift_remove_entry(appliance_key);
        }
    }

    pub(crate) fn appliances_contain_key(&self, name: &ApplianceKey) -> bool {
        self.input
            .appliances
            .as_ref()
            .is_some_and(|appliances| appliances.contains_key(name))
    }

    pub(crate) fn appliance_key_has_reference(
        &self,
        key: &ApplianceKey,
        reference: &ApplianceReference,
    ) -> bool {
        self.input.appliances.as_ref().is_some_and(|appliances| {
            appliances.get(key).is_some_and(|appliance| {
                if let ApplianceEntry::Reference(appliance_reference) = appliance {
                    appliance_reference == reference
                } else {
                    false
                }
            })
        })
    }

    pub(crate) fn appliance_with_key(&self, key: &ApplianceKey) -> Option<&ApplianceEntry> {
        self.input
            .appliances
            .as_ref()
            .and_then(|appliances| appliances.get(key))
    }

    pub(crate) fn clone_appliances(&self) -> IndexMap<ApplianceKey, ApplianceEntry> {
        self.input.appliances.clone().unwrap_or_default()
    }

    pub(crate) fn tariff_schedule(&self) -> Option<&NumericSchedule> {
        self.input.tariff.as_ref().map(|tariff| &tariff.schedule)
    }

    pub(crate) fn keys_for_appliance_gains_with_load_shifting(&self) -> Vec<String> {
        self.input
            .appliance_gains
            .iter()
            .filter(|&(_, gain)| gain.load_shifting.is_some())
            .map(|(key, _)| key.clone())
            .collect::<Vec<_>>()
    }

    pub(crate) fn set_load_shifting_demand_timeseries_for_appliance(
        &mut self,
        appliance_name: &str,
        timeseries: Vec<f64>,
    ) {
        if let Some(ApplianceEntry::Object(Appliance {
            load_shifting: Some(ref mut load_shifting),
            ..
        })) = ApplianceKey::try_from(appliance_name)
            .ok()
            .and_then(|appliance_key| {
                self.input
                    .appliances
                    .as_mut()
                    .map(|gains| gains.get_mut(&appliance_key))
            })
            .flatten()
        {
            load_shifting.demand_timeseries = Some(timeseries);
        }
    }

    pub(crate) fn mechanical_ventilations_for_processing(
        &mut self,
    ) -> Vec<&mut impl MechanicalVentilationForProcessing> {
        self.input
            .infiltration_ventilation
            .mechanical_ventilation
            .values_mut()
            .collect::<Vec<_>>()
    }

    pub(crate) fn keyed_mechanical_ventilations_for_processing(
        &mut self,
    ) -> &mut IndexMap<String, impl MechanicalVentilationForProcessing> {
        &mut self.input.infiltration_ventilation.mechanical_ventilation
    }

    pub(crate) fn has_mechanical_ventilation(&self) -> bool {
        !self
            .input
            .infiltration_ventilation
            .mechanical_ventilation
            .is_empty()
    }

    pub(crate) fn reset_mechanical_ventilation(&mut self) {
        self.input.infiltration_ventilation.mechanical_ventilation = Default::default();
    }

    pub(crate) fn add_mechanical_ventilation(
        &mut self,
        vent_name: &str,
        mech_vent: Value,
    ) -> anyhow::Result<()> {
        self.input
            .infiltration_ventilation
            .mechanical_ventilation
            .insert(vent_name.to_owned(), serde_json::from_value(mech_vent)?);

        Ok(())
    }

    pub(crate) fn appliance_gains_events(
        &self,
    ) -> IndexMap<String, Vec<ApplianceGainsDetailsEvent>> {
        self.input
            .appliance_gains
            .iter()
            .map(|(name, gain)| {
                (
                    name.clone(),
                    gain.events.as_ref().unwrap_or(&vec![]).clone(),
                )
            })
            .collect()
    }

    pub(crate) fn set_window_adjust_control_for_infiltration_ventilation(&mut self, control: &str) {
        self.input.infiltration_ventilation.window_adjust_control = Some(control.to_owned());
    }

    pub(crate) fn infiltration_ventilation_is_noise_nuisance(&self) -> bool {
        self.input
            .infiltration_ventilation
            .noise_nuisance
            .unwrap_or(false)
    }

    pub(crate) fn heat_source_wet(&self) -> Option<&IndexMap<String, HeatSourceWetDetails>> {
        self.input.heat_source_wet.as_ref()
    }

    pub(crate) fn cold_water_source(&self) -> ColdWaterSourceInput {
        self.input.cold_water_source.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;
    use std::fs::File;
    use walkdir::{DirEntry, WalkDir};

    #[fixture]
    fn files() -> Vec<DirEntry> {
        WalkDir::new("./examples/input")
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
    fn should_successfully_parse_all_demo_files(files: Vec<DirEntry>) {
        for entry in files {
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
    fn test_all_demo_files_pass_input_schema(files: Vec<DirEntry>) {
        let schema = serde_json::from_str(include_str!("../schemas/input.schema.json")).unwrap();
        let validator = jsonschema::validator_for(&schema).unwrap();

        let mut erroring_files: usize = Default::default();
        let mut error_outputs: Vec<String> = Default::default();

        for entry in files {
            let json_to_validate =
                serde_json::from_reader(BufReader::new(File::open(entry.path()).unwrap())).unwrap();
            if let Err(errors) = validator.validate(&json_to_validate) {
                erroring_files += 1;
                for error in errors {
                    error_outputs.push(format!(
                        "{} at path \"{}\":\n{error:?}",
                        entry.file_name().to_str().unwrap().to_owned(),
                        error.instance_path.to_string()
                    ));
                }
            };
        }

        assert!(
            error_outputs.is_empty(),
            "Schema created {} individual validation errors across {} files.\n\n{}",
            error_outputs.len(),
            erroring_files,
            error_outputs.join("\n\n")
        );
    }
}
