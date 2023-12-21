use crate::external_conditions::{DaylightSavingsConfig, ShadingSegment, WindowShadingObject};
use crate::simulation_time::SimulationTime;
use indexmap::IndexMap;
use serde::Deserialize;
use serde_json::Value;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;

pub fn parse_input_file(file: &Path) -> Result<Input, Box<dyn Error>> {
    let file = File::open(file)?;
    let reader = BufReader::new(file);

    let project_data: Input = serde_json::from_reader(reader)?;

    Ok(project_data)
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "PascalCase", deny_unknown_fields)]
pub struct Input {
    simulation_time: Arc<SimulationTime>,
    external_conditions: Arc<ExternalConditionsInput>,
    internal_gains: InternalGains,
    appliance_gains: ApplianceGains,
    cold_water_source: ColdWaterSource,
    energy_supply: EnergySupplyInput,
    control: Control,
    hot_water_source: HotWaterSource,
    shower: Option<Shower>,
    bath: Option<Bath>,
    #[serde(rename(deserialize = "Other"))]
    other_water_use: Option<OtherWaterUse>,
    #[serde(rename(deserialize = "Distribution"))]
    water_distribution: Option<WaterDistribution>,
    #[serde(rename(deserialize = "Events"))]
    water_heating_events: WaterHeatingEvents,
    space_heat_system: Option<SpaceHeatSystem>,
    space_cool_system: Option<SpaceCoolSystem>,
    ventilation: Option<Ventilation>,
    infiltration: Infiltration,
    zone: ZoneDictionary,
    #[serde(rename(deserialize = "PartGcompliance"))]
    part_g_compliance: Option<bool>,
    #[serde(rename(deserialize = "PartO_active_cooling_required"))]
    part_o_active_cooling_required: Option<bool>,
    ground_floor_area: Option<f64>,
    number_of_bedrooms: Option<u32>,
    number_of_wet_rooms: Option<u32>,
    heating_control_type: Option<HeatingControlType>,
    #[serde(rename(deserialize = "WaterHeatSchedDefault"))]
    default_water_heating_schedule: Option<WaterHeatingSchedule>,
    heat_source_wet: Option<HeatSourceWet>,
    #[serde(rename(deserialize = "WWHRS"))]
    waste_water_heat_recovery: Option<WasteWaterHeatRecovery>,
    on_site_generation: Option<OnSiteGeneration>,
    #[serde(rename(deserialize = "Window_Opening_For_Cooling"))]
    window_opening_for_cooling: Option<WindowOpeningForCooling>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExternalConditionsInput {
    air_temperatures: Option<Vec<f64>>,
    wind_speeds: Option<Vec<f64>>,
    ground_temperatures: Option<Vec<f64>>,
    diffuse_horizontal_radiation: Option<Vec<f64>>,
    direct_beam_radiation: Option<Vec<f64>>,
    solar_reflectivity_of_ground: Option<Vec<f64>>,
    latitude: Option<f64>,
    longitude: Option<f64>,
    timezone: Option<u32>,
    start_day: Option<u32>,
    end_day: Option<u32>,
    time_series_step: Option<f64>,
    january_first: Option<u32>,
    daylight_savings: Option<DaylightSavingsConfig>,
    leap_day_included: Option<bool>,
    direct_beam_conversion_needed: Option<bool>,
    shading_segments: Vec<ShadingSegment>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InternalGains {
    #[serde(alias = "total internal gains")]
    total_internal_gains: Option<InternalGainsDetails>,
    #[serde(alias = "metabolic gains")]
    metabolic_gains: Option<InternalGainsDetails>,
    other: Option<InternalGainsDetails>,
}

#[derive(Debug, Deserialize)]
pub struct InternalGainsDetails {
    start_day: u32,
    time_series_step: f64,
    schedule: InternalGainsSchedule,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InternalGainsSchedule {
    main: Value, // TODO: possible values are too undefined and unpredictable to reverse-engineer at time of writing! (2023-07-06)
    day: Option<Vec<f64>>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ApplianceGains {
    lighting: Option<ApplianceGainsDetails>,
    cooking: Option<ApplianceGainsDetails>,
    cooking1: Option<ApplianceGainsDetails>, // TODO not sure how stable these numbered keys are but we'll go with this for now
    cooking2: Option<ApplianceGainsDetails>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ApplianceGainsDetails {
    #[serde(rename(deserialize = "type"))]
    gain_type: Option<ApplianceGainType>,
    start_day: u32,
    time_series_step: f64,
    gains_fraction: f64,
    #[serde(alias = "EnergySupply")]
    energy_supply: EnergySupplyType,
    schedule: Schedule,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "snake_case")]
pub enum ApplianceGainType {
    Lighting,
    Cooking,
}

// NB. assuming for now the fields in this struct map to the FuelCode enum
// there may be a may to map this a priori but keeping them manually in sync for now
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EnergySupplyInput {
    #[serde(alias = "mains elec", alias = "mains_elec")]
    mains_electricity: Option<EnergySupplyDetails>,
    #[serde(alias = "mains gas")]
    mains_gas: Option<EnergySupplyDetails>,
    #[serde(rename(deserialize = "bulk LPG"))]
    bulk_lpg: Option<EnergySupplyDetails>,
    #[serde(alias = "heat network")]
    heat_network: Option<HeatNetwork>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EnergySupplyDetails {
    fuel: EnergySupplyType,
    diverter: Option<EnergyDiverter>,
    #[serde(rename(deserialize = "ElectricBattery"))]
    electric_battery: Option<ElectricBattery>,
}

#[derive(Clone, Copy, Debug, Deserialize)]
pub enum EnergySupplyType {
    #[serde(alias = "mains elec", alias = "mains_elec", alias = "electricity")]
    Electricity,
    #[serde(alias = "mains_gas", alias = "mains gas")]
    MainsGas,
    #[serde(rename(deserialize = "unmet_demand"))]
    UnmetDemand,
    #[serde(rename(deserialize = "custom"))]
    Custom,
    #[serde(alias = "bulk LPG", alias = "LPG_bulk")]
    LpgBulk,
    #[serde(rename(deserialize = "LPG_bottled"))]
    LpgBottled,
    #[serde(rename(deserialize = "LPG_condition_11F"))]
    LpgCondition11F,
    #[serde(rename(deserialize = "heat network"))]
    HeatNetwork,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "PascalCase")]
#[serde(deny_unknown_fields)]
pub struct EnergyDiverter {
    storage_tank: StorageTankType,
    heat_source: DiverterHeatSourceType,
}

#[derive(Debug, Deserialize)]
pub enum StorageTankType {
    #[serde(rename(deserialize = "hw cylinder"))]
    HotWaterCylinder,
}

#[derive(Debug, Deserialize)]
pub enum DiverterHeatSourceType {
    #[serde(rename(deserialize = "immersion"))]
    Immersion,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ElectricBattery {
    capacity: f64,
    charge_discharge_efficiency: f64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HeatNetwork {
    fuel: EnergySupplyType,
    factor: HeatNetworkFactor,
}

pub type HeatNetworkFactor = HashMap<String, f64>; // don't really know what these values can be yet

#[derive(Debug, Deserialize)]
pub struct ColdWaterSource {
    #[serde(alias = "mains water")]
    mains_water: Option<ColdWaterSourceDetails>,
    #[serde(alias = "header tank")]
    header_tank: Option<ColdWaterSourceDetails>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ColdWaterSourceDetails {
    start_day: u32,
    temperatures: Vec<f64>,
    time_series_step: f64,
}

pub type Schedule = HashMap<String, Value>; // TODO: possible values are too undefined and unpredictable to reverse-engineer at time of writing! (2023-07-06)

pub type Control = HashMap<String, ControlDetails>;

#[derive(Debug, Deserialize)]
#[serde(tag = "type", deny_unknown_fields)]
pub enum ControlDetails {
    OnOffTimeControl {
        start_day: u32,
        time_series_step: f64,
        advanced_start: Option<u32>,
        logic_type: Option<ControlLogicType>,
        schedule: Schedule,
    },
    OnOffCostMinimisingTimeControl {
        start_day: u32,
        time_series_step: f64,
        advanced_start: Option<u32>,
        logic_type: Option<ControlLogicType>,
        time_on_daily: Option<f64>,
        schedule: Schedule,
    },
    SetpointTimeControl {
        start_day: u32,
        time_series_step: f64,
        advanced_start: Option<u32>,
        logic_type: Option<ControlLogicType>,
        setpoint_min: Option<f64>,
        schedule: Schedule,
    },
    ToUChargeControl {
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
pub enum ControlLogicType {
    Manual,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HotWaterSource {
    #[serde[rename(deserialize = "hw cylinder")]]
    hot_water_cylinder: HotWaterSourceDetails,
}

#[derive(Debug, Deserialize)]
#[serde(tag = "type", deny_unknown_fields)]
pub enum HotWaterSourceDetails {
    StorageTank {
        volume: f64,
        daily_losses: f64,
        min_temp: f64,
        setpoint_temp: f64,
        #[serde(rename(deserialize = "Control_hold_at_setpnt"))]
        control_hold_at_setpoint: Option<Value>, // randomly appears in one of the demo files, seems anomalous
        #[serde(rename(deserialize = "ColdWaterSource"))]
        cold_water_source: ColdWaterSourceType,
        #[serde(rename(deserialize = "HeatSource"))]
        heat_source: Value, // to be finalised as it is not quite clear what form this node can take
        primary_pipework: Option<WaterPipework>,
    },
    CombiBoiler {
        #[serde(rename(deserialize = "ColdWaterSource"))]
        cold_water_source: ColdWaterSourceType,
        #[serde(rename(deserialize = "HeatSourceWet"))]
        heat_source_wet: HeatSourceWetType,
        #[serde(alias = "Control")]
        control: HeatSourceControlType,
        #[allow(non_snake_case)]
        separate_DHW_tests: Value, // only known value here is "M&L" so looks too early to say this is an enum
        rejected_energy_1: f64,
        fuel_energy_2: f64,
        rejected_energy_2: f64,
        storage_loss_factor_2: f64,
        rejected_factor_3: f64,
        #[serde(rename(deserialize = "daily_HW_usage"))]
        daily_hot_water_usage: f64,
    },
    HIU {
        #[serde(rename(deserialize = "ColdWaterSource"))]
        cold_water_source: ColdWaterSourceType,
        #[serde(rename(deserialize = "HeatSourceWet"))]
        heat_source_wet: HeatSourceWetType,
        #[serde(alias = "Control")]
        control: HeatSourceControlType,
    },
    PointOfUse {
        power: f64,
        efficiency: f64,
        #[serde[rename(deserialize = "EnergySupply")]]
        energy_supply: EnergySupplyType,
        #[serde(rename(deserialize = "ColdWaterSource"))]
        cold_water_source: ColdWaterSourceType,
    },
}

#[derive(Debug, Deserialize)]
pub enum ColdWaterSourceType {
    #[serde(rename(deserialize = "mains water"))]
    MainsWater,
    #[serde(rename(deserialize = "header tank"))]
    HeaderTank,
}

#[derive(Debug, Deserialize)]
pub enum HeatSourceWetType {
    #[serde(alias = "boiler")]
    Boiler,
    HeatNetwork,
    HeatPump,
}

#[derive(Debug, Deserialize)]
pub enum HeatSourceControlType {
    #[serde(rename(deserialize = "hw timer"))]
    HotWaterTimer,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WaterPipework {
    internal_diameter_mm: f64,
    external_diameter_mm: f64,
    length: f64,
    insulation_thermal_conductivity: f64,
    insulation_thickness_mm: f64,
    surface_reflectivity: bool,
    pipe_contents: WaterPipeContentsType,
}

#[derive(Debug, Deserialize)]
pub enum WaterPipeContentsType {
    #[serde(alias = "water")]
    Water,
}

#[derive(Debug, Deserialize)]
pub struct Shower {
    mixer: MixerShower,
    #[serde(rename(deserialize = "IES"))]
    ies: Option<InstantElectricShower>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MixerShower {
    #[serde(rename(deserialize = "type"))]
    shower_type: ShowerType,
    flowrate: f64,
    #[serde(rename(deserialize = "ColdWaterSource"))]
    cold_water_source: ColdWaterSourceType,
    #[serde(rename(deserialize = "WWHRS"))]
    waste_water_heat_recovery: Option<Value>, // unclear what these can be yet
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InstantElectricShower {
    #[serde(rename(deserialize = "type"))]
    shower_type: ShowerType, // somewhat of a redundant value possibly
    rated_power: f64,
    #[serde(rename(deserialize = "ColdWaterSource"))]
    cold_water_source: ColdWaterSourceType,
    #[serde(rename(deserialize = "EnergySupply"))]
    energy_supply: EnergySupplyType,
}

#[derive(Debug, Deserialize)]
pub enum ShowerType {
    MixerShower,
    #[serde(alias = "InstantElecShower")]
    InstantElectricShower,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Bath {
    medium: Option<BathDetails>,
}

#[derive(Debug, Deserialize)]
pub struct BathDetails {
    size: u32,
    #[serde(rename(deserialize = "ColdWaterSource"))]
    cold_water_source: ColdWaterSourceType,
    flowrate: f64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OtherWaterUse {
    other: Option<OtherWaterUseDetails>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OtherWaterUseDetails {
    flowrate: f64,
    #[serde(rename(deserialize = "ColdWaterSource"))]
    cold_water_source: ColdWaterSourceType,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WaterDistribution {
    internal: WaterPipework,
    external: WaterPipework,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "PascalCase", deny_unknown_fields)]
pub struct WaterHeatingEvents {
    shower: Option<ShowerEvents>,
    bath: Option<BathEvents>,
    other: Option<OtherWaterHeatingEvents>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WaterHeatingEvent {
    start: f64,
    duration: Option<f64>,
    temperature: f64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ShowerEvents {
    #[serde(alias = "IES")]
    ies: Vec<WaterHeatingEvent>,
    mixer: Vec<WaterHeatingEvent>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BathEvents {
    medium: Vec<WaterHeatingEvent>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OtherWaterHeatingEvents {
    other: Vec<WaterHeatingEvent>,
}

pub type SpaceHeatSystem = HashMap<String, SpaceHeatSystemDetails>;

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields, tag = "type")]
pub enum SpaceHeatSystemDetails {
    #[serde(alias = "InstantElecHeater")]
    InstantElectricHeater {
        temp_setback: Option<f64>,
        rated_power: f64,
        #[serde(rename(deserialize = "EnergySupply"))]
        energy_supply: EnergySupplyType,
        #[serde(alias = "Control")]
        control: Option<String>, // not sure what the possible options are here yet
        frac_convective: f64,
        #[serde(alias = "Zone")]
        zone: Option<String>,
    },
    #[serde(alias = "ElecStorageHeater")]
    ElectricStorageHeater {
        temp_charge_cut: f64,
        rated_power: f64,
        rated_power_instant: f64,
        air_flow_type: String, // don't know what the possible values are here yet
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
        #[serde(rename(deserialize = "EnergySupply"))]
        energy_supply: EnergySupplyType,
        #[serde(alias = "Control")]
        control: Option<String>, // don't know possible options here
        #[serde(rename(deserialize = "ControlCharger"))]
        control_charger: String, // don't know possible options here
        #[serde(alias = "Zone")]
        zone: String, // think these are just arbitrary names?
    },
    WetDistribution {
        advanced_start: Option<u32>,
        thermal_mass: f64,
        c: f64,
        n: f64,
        temp_diff_emit_dsgn: f64,
        frac_convective: f64,
        #[serde(rename(deserialize = "HeatSource"))]
        heat_source: Value, // unclear which values are possible here
        #[serde(alias = "Control")]
        control: Option<String>, // don't know possible values
        ecodesign_controller: EcoDesignController,
        design_flow_temp: i32,
        #[serde(alias = "Zone")]
        zone: String, // as above, these are likely arbitrary names
    },
    WarmAir {
        temp_diff_emit_dsgn: f64,
        frac_convective: f64,
        #[serde(rename(deserialize = "HeatSource"))]
        heat_source: Value,
        #[serde(alias = "Control")]
        control: Option<String>,
    },
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EcoDesignController {
    ecodesign_control_class: u32,
    min_outdoor_temp: Option<i32>,
    max_outdoor_temp: Option<i32>,
    min_flow_temp: Option<i32>,
}

#[derive(Debug, Deserialize)]
#[serde(tag = "type", deny_unknown_fields)]
pub enum Ventilation {
    #[serde(rename(deserialize = "NatVent"))]
    NaturalVentilation { req_ach: f64 },
    WHEV {
        req_ach: f64,
        #[serde(rename(deserialize = "SFP"))]
        sfp: f64,
        #[serde(rename(deserialize = "EnergySupply"))]
        energy_supply: EnergySupplyType,
    },
    MVHR {
        req_ach: f64,
        #[serde(rename(deserialize = "SFP"))]
        sfp: f64,
        efficiency: f64,
        #[serde(rename(deserialize = "EnergySupply"))]
        energy_supply: EnergySupplyType,
        ductwork: VentilationDuctwork,
    },
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VentilationDuctwork {
    internal_diameter_mm: f64,
    external_diameter_mm: f64,
    length_in: f64,
    length_out: f64,
    insulation_thermal_conductivity: f64,
    insulation_thickness_mm: f64,
    reflective: bool,
    #[serde(rename(deserialize = "MVHR_location"))]
    mvhr_location: MVHRLocation,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum MVHRLocation {
    Inside,
    Outside,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Infiltration {
    storeys_in_building: u32,
    storey_of_dwelling: Option<u32>,
    shelter: InfiltrationShelterType,
    build_type: InfiltrationBuildType,
    test_result: f64,
    test_type: InfiltrationTestType,
    env_area: f64,
    volume: f64,
    sheltered_sides: u32,
    open_chimneys: u32,
    open_flues: u32,
    closed_fire: u32,
    flues_d: u32,
    flues_e: u32,
    blocked_chimneys: u32,
    extract_fans: u32,
    passive_vents: u32,
    gas_fires: u32,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum InfiltrationShelterType {
    VerySheltered,
    Sheltered,
    Normal,
    Exposed,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum InfiltrationBuildType {
    House,
    Flat,
}

#[derive(Debug, Deserialize)]
pub enum InfiltrationTestType {
    #[serde(rename(deserialize = "50Pa"))]
    FiftyPascals,
    #[serde(rename(deserialize = "4Pa"))]
    FourPascals,
}

pub type ZoneDictionary = HashMap<String, ZoneInput>;

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ZoneInput {
    #[serde(rename(deserialize = "SpaceHeatSystem"))]
    space_heat_system: String,
    #[serde(rename(deserialize = "SpaceCoolSystem"))]
    space_cool_system: Option<String>,
    #[serde(rename(deserialize = "SpaceHeatControl"))]
    space_heat_control: Option<String>, // don't know what the options are yet
    #[serde(rename(deserialize = "Control_WindowOpening"))]
    control_window_opening: Option<String>, // don't know what the options are yet
    area: f64,
    volume: f64,
    #[serde(rename(deserialize = "Lighting"))]
    lighting: Option<ZoneLighting>,
    temp_setpnt_heat: Option<f64>,
    temp_setpnt_cool: Option<f64>,
    temp_setpnt_init: Option<f64>,
    #[serde(rename(deserialize = "BuildingElement"))]
    building_elements: IndexMap<String, BuildingElement>,
    #[serde(rename(deserialize = "ThermalBridging"))]
    thermal_bridging: Value, // this can be either a float or a hashmap of thermal bridging details - see commented out structs below
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ZoneLighting {
    efficacy: f64,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(tag = "type", deny_unknown_fields)]
pub enum BuildingElement {
    #[serde(rename(deserialize = "BuildingElementOpaque"))]
    Opaque {
        a_sol: f64,
        u_value: Option<f64>,
        r_c: Option<f64>,
        k_m: f64,
        mass_distribution_class: MassDistributionClass,
        is_external_door: Option<bool>,
        pitch: f64,
        orientation360: f64,
        base_height: f64,
        height: f64,
        width: f64,
        area: f64,
        h_ci: Option<f64>,
        h_ri: Option<f64>,
        h_ce: Option<f64>,
        h_re: Option<f64>,
    },
    #[serde(rename(deserialize = "BuildingElementTransparent"))]
    Transparent {
        u_value: Option<f64>,
        area: Option<f64>,
        r_c: Option<f64>,
        pitch: f64,
        orientation360: f64,
        g_value: f64,
        frame_area_fraction: f64,
        base_height: f64,
        height: f64,
        width: f64,
        shading: Vec<WindowShadingObject>,
    },
    #[serde(rename(deserialize = "BuildingElementGround"))]
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
    #[serde(rename(deserialize = "BuildingElementAdjacentZTC"))]
    AdjacentZTC {
        area: f64,
        pitch: f64,
        u_value: Option<f64>,
        r_c: Option<f64>,
        k_m: f64,
        mass_distribution_class: MassDistributionClass,
    },
    #[serde(rename(deserialize = "BuildingElementAdjacentZTU_Simple"))]
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

#[derive(Copy, Clone, Debug, Deserialize)]
pub enum MassDistributionClass {
    D,
    E,
    I,
    IE,
    M,
}

// following _should_ work in theory for deserializing a Zone's ThermalBridging value which is currently
// either a map of thermal bridging details or a single float - going to use serde_json::Value for now
// as this ambiguity will hopefully go away of its own accord as
// #[derive(Debug, Deserialize)]
// #[serde(untagged)]
// pub enum ThermalBridging {
//     ThermalBridgingElements(IndexMap<String, ThermalBridgingDetails>),
//     ThermalBridgingNumber(f64),
// }
//
// #[derive(Debug, Deserialize)]
// #[serde(tag = "type")]
// pub enum ThermalBridgingDetails {
//     #[serde(rename(deserialize = "ThermalBridgeLinear"))]
//     Linear {
//         linear_thermal_transmittance: f64,
//         length: f64,
//     },
//     #[serde(rename(deserialize = "ThermalBridgePoint"))]
//     Point {
//         #[serde(alias = "heat_transfer_coeff")]
//         heat_transfer_coefficient: f64,
//     },
// }

#[derive(Debug, Deserialize)]
pub enum HeatingControlType {
    #[serde(rename(deserialize = "SeparateTempControl"))]
    SeparateTemperatureControl,
}

pub type SpaceCoolSystem = HashMap<String, SpaceCoolSystemDetails>;

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SpaceCoolSystemDetails {
    #[serde(rename(deserialize = "type"))]
    system_type: SpaceCoolSystemType,
    temp_setback: Option<f64>,
    cooling_capacity: f64,
    efficiency: f64,
    frac_convective: f64,
    #[serde(rename(deserialize = "EnergySupply"))]
    energy_supply: EnergySupplyType,
    #[serde(rename(deserialize = "Control"))]
    control: Option<String>,
}

#[derive(Debug, Deserialize)]
pub enum SpaceCoolSystemType {
    AirConditioning,
}

#[derive(Debug, Deserialize)]
pub enum WaterHeatingSchedule {
    AllDay,
    HeatingHours,
}

pub type HeatSourceWet = HashMap<String, HeatSourceWetDetails>;

#[derive(Debug, Deserialize)]
#[serde(tag = "type", deny_unknown_fields)]
pub enum HeatSourceWetDetails {
    HeatPump {
        #[serde(rename(deserialize = "EnergySupply"))]
        energy_supply: EnergySupplyType,
        source_type: HeatPumpSourceType,
        #[serde(rename(deserialize = "EnergySupply_heat_network"))]
        energy_supply_heat_network: Option<String>, // unclear what this is
        temp_distribution_heat_network: Option<f64>,
        sink_type: HeatPumpSinkType,
        #[serde(rename(deserialize = "backup_ctrl_type"))]
        backup_control_type: HeatPumpBackupControlType,
        time_delay_backup: f64,
        modulating_control: bool,
        min_modulation_rate_20: Option<f64>,
        min_modulation_rate_35: Option<f64>,
        min_modulation_rate_55: Option<f64>,
        time_constant_onoff_operation: u32,
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
        test_data: Vec<HeatSourceTestDatum>,
    },
    Boiler {
        #[serde(rename(deserialize = "EnergySupply"))]
        energy_supply: EnergySupplyType,
        #[serde(rename(deserialize = "EnergySupply_aux"))]
        energy_supply_auxiliary: Option<EnergySupplyType>,
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
        #[serde(rename(deserialize = "EnergySupply"))]
        energy_supply: EnergySupplyType,
        heat_battery_location: HeatSourceLocation,
        electricity_circ_pump: f64,
        electricity_standby: f64,
        rated_charge_power: f64,
        heat_storage_capacity: f64,
        max_rated_heat_output: f64,
        max_rated_losses: f64,
        number_of_units: u32,
        #[serde(rename(deserialize = "ControlCharge"))]
        control_charge: String,
    },
    HIU {
        #[serde(rename(deserialize = "EnergySupply"))]
        energy_supply: EnergySupplyType,
        power_max: f64,
        #[serde(rename(deserialize = "HIU_daily_loss"))]
        hiu_daily_loss: f64,
        building_level_distribution_losses: f64,
    },
}

#[derive(Debug, Deserialize)]
pub enum HeatPumpSourceType {
    OutsideAir,
    ExhaustAirMEV,
    HeatNetwork,
}

#[derive(Debug, Deserialize)]
pub enum HeatPumpSinkType {
    Water,
    Air,
}

#[derive(Debug, Deserialize)]
pub enum HeatPumpBackupControlType {
    TopUp,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HeatSourceTestDatum {
    air_flow_rate: Option<f64>,
    test_letter: char,
    capacity: f64,
    cop: f64,
    #[serde(rename(deserialize = "degradation_coeff"))]
    degradation_coefficient: f64,
    design_flow_temp: i32,
    temp_outlet: i32,
    temp_source: i32,
    temp_test: i32,
}

#[derive(Debug, Deserialize)]
pub enum HeatSourceLocation {
    #[serde(alias = "internal")]
    Internal,
}

pub type WasteWaterHeatRecovery = HashMap<String, WasteWaterHeatRecoveryDetails>;

#[derive(Debug, Deserialize)]
pub struct WasteWaterHeatRecoveryDetails {
    #[serde(rename(deserialize = "type"))]
    system_type: String, // unclear what these can be yet
    #[serde(rename(deserialize = "ColdWaterSource"))]
    cold_water_source: ColdWaterSourceType,
    flow_rates: Vec<f64>,
    efficiencies: Vec<f64>,
    utilisation_factor: f64,
}

pub type OnSiteGeneration = HashMap<String, OnSiteGenerationDetails>;

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OnSiteGenerationDetails {
    #[serde(rename(deserialize = "type"))]
    generation_type: OnSiteGenerationType,
    peak_power: f64,
    ventilation_strategy: OnSiteGenerationVentilationStrategy,
    pitch: f64,
    orientation360: f64,
    base_height: f64,
    height: f64,
    width: f64,
    #[serde(rename(deserialize = "EnergySupply"))]
    energy_supply: EnergySupplyType,
}

#[derive(Debug, Deserialize)]
pub enum OnSiteGenerationType {
    PhotovoltaicSystem,
}

#[derive(Debug, Deserialize)]
pub enum OnSiteGenerationVentilationStrategy {
    #[serde(alias = "moderately_ventilated")]
    ModeratelyVentilated,
}

#[derive(Debug, Deserialize)]
pub struct WindowOpeningForCooling {
    equivalent_area: f64,
    // control: String,
}

#[cfg(test)]
mod test {
    use super::*;
    use rstest::*;
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
            let parsed = parse_input_file(entry.path());
            assert!(
                parsed.is_ok(),
                "error was {:?} when parsing file {}",
                parsed.err().unwrap(),
                entry.file_name().to_str().unwrap()
            );
        }
    }
}
