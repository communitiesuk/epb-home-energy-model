use crate::external_conditions::{DaylightSavingsConfig, ShadingSegment, WindowShadingObject};
use crate::simulation_time::SimulationTime;
use anyhow::{anyhow, bail};
use arrayvec::ArrayString;
use indexmap::{Equivalent, IndexMap};
use serde::{Deserialize, Deserializer, Serialize};
use serde_enum_str::Deserialize_enum_str;
use serde_json::{json, Value};
use std::borrow::Borrow;
use std::fmt::{Display, Formatter};
use std::io::{BufReader, Read};
use std::sync::Arc;

pub fn ingest_for_processing(json: impl Read) -> Result<InputForProcessing, anyhow::Error> {
    InputForProcessing::init_with_json(json)
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(rename_all = "PascalCase", deny_unknown_fields)]
pub struct Input {
    pub simulation_time: SimulationTime,
    pub external_conditions: Arc<ExternalConditionsInput>,
    pub internal_gains: InternalGains,
    pub appliance_gains: ApplianceGains,
    pub cold_water_source: ColdWaterSourceInput,
    pub energy_supply: EnergySupplyInput,
    #[serde(deserialize_with = "deserialize_control")]
    pub control: Control,
    pub hot_water_source: HotWaterSource,
    pub shower: Option<Shower>,
    pub bath: Option<Bath>,
    #[serde(rename = "Other")]
    pub other_water_use: Option<OtherWaterUse>,
    #[serde(rename = "Distribution")]
    pub water_distribution: Option<WaterDistribution>,
    #[serde(rename = "Events")]
    pub water_heating_events: WaterHeatingEvents,
    pub space_heat_system: Option<SpaceHeatSystem>,
    pub space_cool_system: Option<SpaceCoolSystem>,
    pub ventilation: Option<Ventilation>,
    pub infiltration: Infiltration,
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
    #[allow(dead_code)]
    number_of_bedrooms: Option<usize>,
    #[allow(dead_code)]
    number_of_wet_rooms: Option<usize>,
    #[allow(dead_code)]
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
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ExternalConditionsInput {
    pub air_temperatures: Option<Vec<f64>>,
    pub wind_speeds: Option<Vec<f64>>,
    // check upstream whether anything uses this
    #[serde(rename = "ground_temperatures")]
    _ground_temperatures: Option<Vec<f64>>,
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
#[serde(deny_unknown_fields)]
pub struct InternalGains {
    #[serde(alias = "total internal gains")]
    pub total_internal_gains: Option<InternalGainsDetails>,
    #[serde(alias = "metabolic gains")]
    pub metabolic_gains: Option<InternalGainsDetails>,
    pub evaporative_losses: Option<InternalGainsDetails>,
    pub other: Option<InternalGainsDetails>,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub struct InternalGainsDetails {
    pub start_day: u32,
    pub time_series_step: f64,
    pub schedule: InternalGainsSchedule,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct InternalGainsSchedule {
    pub main: Value,
    pub day: Option<Value>,
    pub week: Option<Value>,
    pub weekday: Option<Value>,
    pub weekend: Option<Value>,
}

pub type ApplianceGains = IndexMap<String, ApplianceGainsDetails>;

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ApplianceGainsDetails {
    // check upstream whether type is used here
    #[serde(rename(deserialize = "type"))]
    gain_type: Option<ApplianceGainType>,
    pub start_day: u32,
    pub time_series_step: f64,
    pub gains_fraction: f64,
    #[serde(alias = "EnergySupply")]
    pub energy_supply: EnergySupplyType,
    pub schedule: Schedule,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields, rename_all = "snake_case")]
pub enum ApplianceGainType {
    Lighting,
    Cooking,
}

#[derive(Clone, Copy, Debug, Deserialize_enum_str, Eq, Hash, PartialEq, Serialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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
                "mains gas" => "mains elec",
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
        value.try_into().unwrap()
    }
}

pub type EnergySupplyInput = IndexMap<EnergySupplyKey, EnergySupplyDetails>;

#[derive(Clone, Debug, Deserialize, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(tag = "fuel", deny_unknown_fields, rename_all = "snake_case")]
pub struct EnergySupplyDetails {
    pub fuel: FuelType,
    pub diverter: Option<EnergyDiverter>,
    #[serde(rename = "ElectricBattery")]
    pub electric_battery: Option<ElectricBattery>,
    pub factor: Option<CustomEnergySourceFactor>,
}

/// TODO clarify further
/// It's not completely clear at the moment what the difference between fuel type and energy supply type is,
/// but electricity and gas each seem to be indicated using different strings between fuel and energy supply
/// in the input examples, so keeping them separate for the time being
/// (It's also hard to see some of these as types of fuel)
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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
        let json_string = serde_json::to_string(self).unwrap();
        write!(f, "{}", json_string)
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
pub enum EnergySupplyType {
    #[serde(rename = "mains elec")]
    Electricity,
    #[serde(rename = "mains gas")]
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
        let json_string = serde_json::to_string(self).unwrap();
        write!(f, "{}", json_string)
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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
#[serde(deny_unknown_fields)]
pub struct ElectricBattery {
    pub capacity: f64,
    pub charge_discharge_efficiency: f64,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct CustomEnergySourceFactor {
    #[serde(rename = "Emissions Factor kgCO2e/kWh")]
    pub emissions: f64,
    #[serde(rename = "Emissions Factor kgCO2e/kWh including out-of-scope emissions")]
    pub emissions_including_out_of_scope: f64,
    #[serde(rename = "Primary Energy Factor kWh/kWh delivered")]
    pub primary_energy_factor: f64,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ColdWaterSourceDetails {
    pub start_day: u32,
    pub temperatures: Vec<f64>,
    pub time_series_step: f64,
}

pub type Schedule = IndexMap<String, Value>; // TODO: possible values are too undefined and unpredictable to reverse-engineer at time of writing! (2023-07-06)

pub type CoreControls = Vec<HeatSourceControl>;

pub type ExtraControls = IndexMap<String, ControlDetails>;

#[derive(Debug)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub struct Control {
    pub core: CoreControls,
    pub extra: ExtraControls,
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
#[serde(tag = "type", deny_unknown_fields)]
pub enum ControlDetails {
    #[serde(rename = "OnOffTimeControl")]
    OnOffTime {
        start_day: u32,
        time_series_step: f64,
        advanced_start: Option<u32>,
        logic_type: Option<ControlLogicType>,
        schedule: Schedule,
    },
    #[serde(rename = "OnOffCostMinimisingTimeControl")]
    OnOffCostMinimisingTime {
        start_day: u32,
        time_series_step: f64,
        advanced_start: Option<u32>,
        logic_type: Option<ControlLogicType>,
        time_on_daily: Option<f64>,
        schedule: Schedule,
    },
    #[serde(rename = "SetpointTimeControl")]
    SetpointTime {
        start_day: u32,
        time_series_step: f64,
        advanced_start: Option<f64>,
        logic_type: Option<ControlLogicType>,
        setpoint_min: Option<f64>,
        setpoint_max: Option<f64>,
        default_to_max: Option<bool>,
        schedule: Schedule,
    },
    #[serde(rename = "ToUChargeControl")]
    ToUCharge {
        start_day: u32,
        time_series_step: f64,
        advanced_start: Option<u32>,
        logic_type: Option<ControlLogicType>,
        charge_level: Option<Value>,
        target_charge: Option<f64>,
        schedule: Schedule,
    },
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub enum ControlLogicType {
    Manual,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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
#[serde(tag = "type", deny_unknown_fields)]
pub enum HotWaterSourceDetails {
    StorageTank {
        volume: f64,
        daily_losses: f64,
        min_temp: f64,
        setpoint_temp: f64,
        #[serde(rename = "Control_hold_at_setpnt")]
        control_hold_at_setpoint: Option<String>,
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceType,
        #[serde(rename = "HeatSource")]
        heat_source: IndexMap<String, HeatSource>,
        primary_pipework: Option<WaterPipework>,
    },
    CombiBoiler {
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceType,
        #[serde(rename = "HeatSourceWet")]
        heat_source_wet: HeatSourceWetType,
        #[serde(alias = "Control")]
        control: HeatSourceControlType,
        #[serde(rename = "separate_DHW_tests")]
        separate_dhw_tests: BoilerHotWaterTest,
        rejected_energy_1: f64,
        fuel_energy_2: f64,
        rejected_energy_2: f64,
        storage_loss_factor_2: f64,
        rejected_factor_3: f64,
        #[serde(rename = "daily_HW_usage")]
        daily_hot_water_usage: f64,
    },
    #[serde(rename(deserialize = "HIU"))]
    Hiu {
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceType,
        #[serde(rename = "HeatSourceWet")]
        heat_source_wet: HeatSourceWetType,
        #[serde(alias = "Control")]
        control: HeatSourceControlType,
    },
    PointOfUse {
        power: f64,
        efficiency: f64,
        #[serde(rename = "EnergySupply")]
        energy_supply: EnergySupplyType,
        #[serde(rename = "ColdWaterSource")]
        cold_water_source: ColdWaterSourceType,
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
}

#[derive(Copy, Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub enum ColdWaterSourceType {
    #[serde(rename = "mains water")]
    MainsWater,
    #[serde(rename = "header tank")]
    HeaderTank,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub enum HeatSourceWetType {
    #[serde(alias = "boiler")]
    Boiler,
    HeatNetwork,
    HeatPump,
}

impl HeatSourceWetType {
    /// Convert the type to a canonical string based on the input format to be used in e.g. energy supply names
    pub fn to_canonical_string(&self) -> String {
        serde_json::to_value(self).unwrap().to_string()
    }
}

#[derive(Clone, Copy, Debug, Deserialize_enum_str)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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
pub enum HeatSourceControl {
    HotWaterTimer(ControlDetails),
    WindowOpening(ControlDetails),
    WindowOpeningLivingRoom(ControlDetails),
    WindowOpeningRestOfDwelling(ControlDetails),
    AlwaysOff(ControlDetails),
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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
        cold_water_source: ColdWaterSourceType,
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
pub struct HeatPumpHotWaterTestData {
    #[serde(rename = "L")]
    pub l: Option<HeatPumpHotWaterOnlyTestDatum>,
    #[serde(rename = "M")]
    pub m: HeatPumpHotWaterOnlyTestDatum,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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
#[serde(deny_unknown_fields)]
pub struct WaterPipework {
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
pub enum WaterPipeContentsType {
    #[serde(rename = "water")]
    Water,
    #[serde(rename = "air")]
    Air,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub struct Shower {
    pub mixer: MixerShower,
    #[serde(rename = "IES")]
    pub ies: Option<InstantElectricShower>,
}

impl Shower {
    /// Provide shower field names as strings.
    pub fn keys(&self) -> Vec<String> {
        let mut keys = vec!["mixer".to_string()];
        if self.ies.is_some() {
            keys.push("IES".to_string());
        }

        keys
    }

    pub fn name_refers_to_instant_electric_shower(&self, name: &str) -> bool {
        name == "ies" && self.ies.is_some()
    }
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct MixerShower {
    #[serde(rename = "type")]
    pub shower_type: ShowerType,
    pub flowrate: f64,
    #[serde(rename = "ColdWaterSource")]
    pub cold_water_source: ColdWaterSourceType,
    #[serde(rename = "WWHRS")]
    pub waste_water_heat_recovery: Option<Value>, // unclear what these can be yet
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct InstantElectricShower {
    #[serde(rename = "type")]
    pub shower_type: ShowerType,
    // somewhat of a redundant value possibly
    pub rated_power: f64,
    #[serde(rename = "ColdWaterSource")]
    pub cold_water_source: ColdWaterSourceType,
    #[serde(rename = "EnergySupply")]
    pub energy_supply: EnergySupplyType,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub enum ShowerType {
    MixerShower,
    #[serde(alias = "InstantElecShower")]
    InstantElectricShower,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Bath {
    pub medium: Option<BathDetails>,
}

impl Bath {
    /// Provide bath field names as strings.
    pub fn keys(&self) -> Vec<String> {
        let mut keys = vec![];
        if self.medium.is_some() {
            keys.push("medium".to_string());
        }

        keys
    }

    pub fn size_for_field(&self, field: &str) -> Option<f64> {
        (field == "medium")
            .then(|| self.medium.as_ref().map(|details| details.size))
            .unwrap_or_default()
    }

    pub fn flowrate_for_field(&self, field: &str) -> Option<f64> {
        (field == "medium")
            .then(|| self.medium.as_ref().map(|details| details.flowrate))
            .unwrap_or_default()
    }
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub struct BathDetails {
    pub size: f64,
    #[serde(rename = "ColdWaterSource")]
    pub cold_water_source: ColdWaterSourceType,
    pub flowrate: f64,
}

#[derive(Clone, Debug, Default, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct OtherWaterUse {
    pub other: Option<OtherWaterUseDetails>,
}

impl OtherWaterUse {
    /// Provide other water use field names as strings.
    pub fn keys(&self) -> Vec<String> {
        let mut keys = vec![];
        if self.other.is_some() {
            keys.push("other".to_string());
        }

        keys
    }

    pub fn flowrate_for_field(&self, field: &str) -> Option<f64> {
        (field == "other")
            .then(|| self.other.as_ref().map(|details| details.flowrate))
            .unwrap_or_default()
    }
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct OtherWaterUseDetails {
    pub flowrate: f64,
    #[serde(rename = "ColdWaterSource")]
    pub cold_water_source: ColdWaterSourceType,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct WaterDistribution {
    pub internal: WaterPipework,
    pub external: WaterPipework,
}

#[derive(Debug, Default, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(rename_all = "PascalCase", deny_unknown_fields)]
pub struct WaterHeatingEvents {
    #[serde(default)]
    pub shower: ShowerEvents,
    #[serde(default)]
    pub bath: BathEvents,
    #[serde(default)]
    pub other: OtherWaterHeatingEvents,
}

impl WaterHeatingEvents {
    fn add_event_for_type_and_name(
        &mut self,
        event_type: WaterHeatingEventType,
        name: &str,
        event: WaterHeatingEvent,
    ) {
        match event_type {
            WaterHeatingEventType::Shower => {
                self.shower.add_event_for_name(name, event);
            }
            WaterHeatingEventType::Bath => {
                self.bath.add_event_for_name(name, event);
            }
            WaterHeatingEventType::Other => {
                self.other.add_event_for_name(name, event);
            }
        }
    }
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct WaterHeatingEvent {
    pub start: f64,
    pub duration: Option<f64>,
    pub temperature: f64,
}

#[derive(Clone, Copy, PartialEq)]
pub enum WaterHeatingEventType {
    Shower,
    Bath,
    Other,
}

#[derive(Debug, Deserialize, Default)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ShowerEvents {
    #[serde(alias = "IES")]
    pub ies: Vec<WaterHeatingEvent>,
    pub mixer: Vec<WaterHeatingEvent>,
}

impl ShowerEvents {
    fn add_event_for_name(&mut self, name: &str, event: WaterHeatingEvent) {
        match name {
            "ies" => {
                self.ies.push(event);
            }
            "mixer" => {
                self.mixer.push(event);
            }
            _ => {}
        }
    }
}

#[derive(Debug, Deserialize, Default)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct BathEvents {
    pub medium: Vec<WaterHeatingEvent>,
}

impl BathEvents {
    fn add_event_for_name(&mut self, name: &str, event: WaterHeatingEvent) {
        if name != "medium" {
            return;
        }
        self.medium.push(event);
    }
}

#[derive(Debug, Deserialize, Default)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct OtherWaterHeatingEvents {
    pub other: Vec<WaterHeatingEvent>,
}

impl OtherWaterHeatingEvents {
    fn add_event_for_name(&mut self, name: &str, event: WaterHeatingEvent) {
        if name != "other" {
            return;
        }
        self.other.push(event);
    }
}

pub type SpaceHeatSystem = IndexMap<String, SpaceHeatSystemDetails>;

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields, tag = "type")]
pub enum SpaceHeatSystemDetails {
    #[serde(alias = "InstantElecHeater")]
    InstantElectricHeater {
        temp_setback: Option<f64>,
        rated_power: f64,
        #[serde(rename = "EnergySupply")]
        energy_supply: EnergySupplyType,
        #[serde(alias = "Control")]
        control: Option<String>,
        // not sure what the possible options are here yet
        frac_convective: f64,
        #[serde(alias = "Zone")]
        zone: Option<String>,
    },
    #[serde(alias = "ElecStorageHeater")]
    ElectricStorageHeater {
        temp_charge_cut: f64,
        rated_power: f64,
        rated_power_instant: f64,
        air_flow_type: String,
        // don't know what the possible values are here yet
        temp_dis_safe: f64,
        thermal_mass: f64,
        frac_convective: f64,
        #[serde(alias = "U_ins")]
        u_ins: f64,
        mass_core: f64,
        c_pcore: f64,
        temp_core_target: f64,
        #[serde(alias = "A_core")]
        a_core: f64,
        c_wall: f64,
        n_wall: f64,
        thermal_mass_wall: f64,
        fan_pwr: f64,
        n_units: u32,
        #[serde(rename = "EnergySupply")]
        energy_supply: EnergySupplyType,
        #[serde(alias = "Control")]
        control: Option<String>,
        // don't know possible options here
        #[serde(rename = "ControlCharger")]
        control_charger: String,
        // don't know possible options here
        #[serde(alias = "Zone")]
        zone: String, // think these are just arbitrary names?
    },
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
        #[serde(alias = "Control")]
        control: Option<String>,
        // check upstream if this is used
        ecodesign_controller: EcoDesignController,
        design_flow_temp: i32,
        #[serde(alias = "Zone")]
        zone: String, // as above, these are likely arbitrary names
    },
    WarmAir {
        temp_diff_emit_dsgn: f64,
        frac_convective: f64,
        #[serde(rename = "HeatSource")]
        heat_source: SpaceHeatSystemHeatSource,
        #[serde(alias = "Control")]
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
#[serde(deny_unknown_fields)]
pub struct SpaceHeatSystemHeatSource {
    pub name: String,
    pub temp_flow_limit_upper: Option<f64>,
}

// it is unclear whether this struct should be used - see reference to the struct above
#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[allow(dead_code)]
#[serde(deny_unknown_fields)]
pub struct EcoDesignController {
    ecodesign_control_class: u32,
    min_outdoor_temp: Option<i32>,
    max_outdoor_temp: Option<i32>,
    min_flow_temp: Option<i32>,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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
    #[serde(rename = "MVHR")]
    Mvhr {
        req_ach: f64,
        #[serde(rename = "SFP")]
        sfp: f64,
        efficiency: f64,
        #[serde(rename = "EnergySupply")]
        energy_supply: EnergySupplyType,
        ductwork: VentilationDuctwork,
    },
}

impl Ventilation {
    pub fn req_ach(&self) -> f64 {
        match self {
            Ventilation::Natural { req_ach, .. } => *req_ach,
            Ventilation::Whev { req_ach, .. } => *req_ach,
            Ventilation::Mvhr { req_ach, .. } => *req_ach,
        }
    }
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct VentilationDuctwork {
    pub internal_diameter_mm: f64,
    pub external_diameter_mm: f64,
    pub length_in: f64,
    pub length_out: f64,
    pub insulation_thermal_conductivity: f64,
    pub insulation_thickness_mm: f64,
    pub reflective: bool,
    #[serde(rename = "MVHR_location")]
    pub mvhr_location: MVHRLocation,
}

#[derive(Copy, Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(rename_all = "lowercase")]
pub enum MVHRLocation {
    Inside,
    Outside,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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
#[serde(rename_all = "lowercase")]
pub enum InfiltrationShelterType {
    VerySheltered,
    Sheltered,
    Normal,
    Exposed,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(rename_all = "lowercase")]
pub enum InfiltrationBuildType {
    House,
    Flat,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub enum InfiltrationTestType {
    #[serde(rename = "50Pa")]
    FiftyPascals,
    #[serde(rename = "4Pa")]
    FourPascals,
}

pub type ZoneDictionary = IndexMap<String, ZoneInput>;

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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
pub enum SpaceHeatControlType {
    #[serde(rename = "livingroom")]
    LivingRoom,
    #[serde(rename = "restofdwelling")]
    RestOfDwelling,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ZoneLighting {
    #[serde(rename = "efficacy")]
    efficacy: f64,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(tag = "type", deny_unknown_fields)]
pub enum BuildingElement {
    #[serde(rename = "BuildingElementOpaque")]
    Opaque {
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
        shading: Vec<WindowShadingObject>,
    },
    #[serde(rename = "BuildingElementGround")]
    Ground {
        area: f64,
        pitch: f64,
        u_value: f64,
        r_f: f64,
        k_m: f64,
        mass_distribution_class: MassDistributionClass,
        h_pi: f64,
        h_pe: f64,
        perimeter: f64,
        psi_wall_floor_junc: f64,
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

// special deserialization logic so that orientations are normalized correctly on the way in
fn deserialize_orientation<'de, D>(deserializer: D) -> Result<f64, D::Error>
where
    D: Deserializer<'de>,
{
    let orientation360_value: f64 = Deserialize::deserialize(deserializer)?;
    Ok(init_orientation(orientation360_value))
}

fn init_orientation(value: f64) -> f64 {
    // Convert orientation from 0-360 (clockwise) to -180 to +180 (anticlockwise)
    180. - value
}

#[derive(Copy, Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub enum MassDistributionClass {
    D,
    E,
    I,
    IE,
    M,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(untagged)]
pub enum ThermalBridging {
    ThermalBridgingElements(IndexMap<String, ThermalBridgingDetails>),
    ThermalBridgingNumber(f64),
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(tag = "type")]
pub enum ThermalBridgingDetails {
    #[serde(rename(deserialize = "ThermalBridgeLinear"))]
    Linear {
        linear_thermal_transmittance: f64,
        length: f64,
    },
    #[serde(rename = "ThermalBridgePoint")]
    Point {
        #[serde(alias = "heat_transfer_coeff")]
        heat_transfer_coefficient: f64,
    },
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub enum HeatingControlType {
    #[serde(rename = "SeparateTimeAndTempControl")]
    SeparateTimeAndTemperatureControl,
    #[serde(rename = "SeparateTempControl")]
    SeparateTemperatureControl,
}

pub type SpaceCoolSystem = IndexMap<String, SpaceCoolSystemDetails>;

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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
pub enum SpaceCoolSystemType {
    AirConditioning,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub enum WaterHeatingSchedule {
    AllDay,
    HeatingHours,
}

pub type HeatSourceWet = IndexMap<String, HeatSourceWetDetails>;

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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
        power_heating_circ_pump: f64,
        power_source_circ_pump: f64,
        power_standby: f64,
        power_crankcase_heater: f64,
        power_off: f64,
        power_max_backup: f64,
        test_data: Vec<HeatPumpTestDatum>,
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

#[derive(Copy, Clone, Debug, Deserialize, PartialEq)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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
pub enum HeatPumpSinkType {
    Water,
    Air,
}

#[derive(Copy, Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub enum HeatPumpBackupControlType {
    None,
    TopUp,
    Substitute,
}

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct HeatPumpTestDatum {
    pub air_flow_rate: Option<f64>,
    pub test_letter: TestLetter,
    pub capacity: f64,
    pub cop: f64,
    #[serde(rename = "degradation_coeff")]
    pub degradation_coefficient: f64,
    pub design_flow_temp: f64,
    pub temp_outlet: f64,
    pub temp_source: f64,
    pub temp_test: f64,
}

pub type TestLetter = ArrayString<2>;

#[derive(Copy, Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub enum HeatSourceLocation {
    #[serde(alias = "internal")]
    Internal,
    #[serde(alias = "external")]
    External,
}

pub type WasteWaterHeatRecovery = IndexMap<String, WasteWaterHeatRecoveryDetails>;

#[derive(Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub struct WasteWaterHeatRecoveryDetails {
    #[serde(rename = "type")]
    pub system_type: WwhrsType,
    #[serde(rename = "ColdWaterSource")]
    pub cold_water_source: ColdWaterSourceType,
    pub flow_rates: Vec<f64>,
    pub efficiencies: Vec<f64>,
    pub utilisation_factor: f64,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub enum WwhrsType {
    #[serde(alias = "WWHRS_InstantaneousSystemA")]
    SystemA,
    #[serde(alias = "WWHRS_InstantaneousSystemB")]
    SystemB,
    #[serde(alias = "WWHRS_InstantaneousSystemC")]
    SystemC,
}

pub type OnSiteGeneration = IndexMap<String, OnSiteGenerationDetails>;

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
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
    },
}

#[derive(Copy, Clone, Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub enum OnSiteGenerationVentilationStrategy {
    #[serde(alias = "unventilated")]
    Unventilated,
    #[serde(alias = "moderately_ventilated")]
    ModeratelyVentilated,
    #[serde(alias = "strongly_or_forced_ventilated")]
    StronglyOrForcedVentilated,
    #[serde(alias = "rear_surface_free")]
    RearSurfaceFree,
}

#[derive(Debug, Deserialize)]
#[cfg_attr(feature = "schemars", derive(schemars::JsonSchema))]
pub struct WindowOpeningForCooling {
    pub equivalent_area: f64,
    // control: String,
}

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

    pub fn reset_internal_gains(&mut self) -> &Self {
        self.input.internal_gains = Default::default();
        self
    }

    pub fn total_zone_area(&self) -> f64 {
        self.input.zone.values().map(|z| z.area).sum::<f64>()
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

    pub fn heating_control_type(&self) -> Option<HeatingControlType> {
        self.input.heating_control_type
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

    pub fn space_cool_system_for_zone(&self, zone: &str) -> anyhow::Result<Option<String>> {
        Ok(self
            .input
            .zone
            .get(zone)
            .ok_or(anyhow!("Used zone key for a zone that does not exist"))?
            .space_cool_system
            .clone())
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
        let lighting_entry = self.input.appliance_gains.entry(field.into());
        lighting_entry.and_modify(|entry| *entry = gains_details);

        Ok(self)
    }

    pub fn appliance_gains_fields_for_cooking(&self) -> Vec<String> {
        self.input
            .appliance_gains
            .iter()
            .filter(|(_field, gains)| matches!(gains.gain_type, Some(ApplianceGainType::Cooking)))
            .map(|(field, _gains)| field.clone())
            .collect()
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

    pub fn fuel_type_for_energy_supply_field(&self, field: &str) -> anyhow::Result<String> {
        Ok(self
            .input
            .energy_supply
            .get(field)
            .ok_or(anyhow!(
                "Fuel type not provided for energy supply field '{field}'"
            ))?
            .fuel
            .to_string())
    }

    pub fn shower_flowrate(&self) -> Option<f64> {
        self.input
            .shower
            .as_ref()
            .map(|shower| shower.mixer.flowrate)
    }

    pub fn reset_water_heating_events(&mut self) {
        self.input.water_heating_events = Default::default();
    }

    pub fn shower_keys(&self) -> Vec<String> {
        match self.input.shower.as_ref() {
            Some(shower) => shower.keys(),
            None => Default::default(),
        }
    }

    pub fn shower_name_refers_to_instant_electric(&self, name: &str) -> bool {
        match &self.input.shower {
            Some(shower) => shower.name_refers_to_instant_electric_shower(name),
            None => false,
        }
    }

    pub fn bath_keys(&self) -> Vec<String> {
        match self.input.bath.as_ref() {
            Some(bath) => bath.keys(),
            None => Default::default(),
        }
    }

    pub fn size_for_bath_field(&self, field: &str) -> Option<f64> {
        self.input
            .bath
            .as_ref()
            .and_then(|bath| bath.size_for_field(field))
    }

    pub fn flowrate_for_bath_field(&self, field: &str) -> Option<f64> {
        self.input
            .bath
            .as_ref()
            .and_then(|bath| bath.flowrate_for_field(field))
    }

    pub fn other_water_use_keys(&self) -> Vec<String> {
        match self.input.other_water_use.as_ref() {
            Some(other) => other.keys(),
            None => Default::default(),
        }
    }

    pub fn flow_rate_for_other_water_use_field(&self, field: &str) -> Option<f64> {
        self.input
            .other_water_use
            .as_ref()
            .and_then(|other| other.flowrate_for_field(field))
    }

    pub fn set_other_water_use_details(
        &mut self,
        cold_water_source_type: impl Into<String>,
        flowrate: f64,
    ) -> anyhow::Result<()> {
        let other_details = Some(OtherWaterUseDetails {
            flowrate,
            cold_water_source: serde_json::from_str(cold_water_source_type.into().as_str())?,
        });
        match self.input.other_water_use {
            Some(ref mut other_water_use) => {
                other_water_use.other = other_details;
            }
            None => {
                self.input.other_water_use = Some(OtherWaterUse {
                    other: other_details,
                });
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;
    use std::fs::File;
    use walkdir::WalkDir;

    #[rstest]
    fn should_successfully_parse_all_demo_files() {
        for entry in WalkDir::new("./examples/input")
            .into_iter()
            .filter_map(Result::ok)
            .filter(|e| {
                !e.file_type().is_dir() && e.file_name().to_str().unwrap().ends_with("json")
            })
        {
            let parsed = ingest_for_processing(File::open(entry.path()).unwrap());
            assert!(
                parsed.is_ok(),
                "error was {:?} when parsing file {}",
                parsed.err().unwrap(),
                entry.file_name().to_str().unwrap()
            );
        }
    }
}
