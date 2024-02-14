use crate::core::common::WaterSourceWithTemperature;
use crate::core::controls::time_control::{
    Control, OnOffMinimisingTimeControl, OnOffTimeControl, SetpointTimeControl, ToUChargeControl,
};
use crate::core::ductwork::Ductwork;
use crate::core::energy_supply::energy_supply::{EnergySupplies, EnergySupply};
use crate::core::heating_systems::boiler::{Boiler, BoilerServiceWaterCombi};
use crate::core::heating_systems::common::HeatSourceWet;
use crate::core::heating_systems::heat_battery::HeatBattery;
use crate::core::heating_systems::heat_network::{HeatNetwork, HeatNetworkServiceWaterDirect};
use crate::core::heating_systems::heat_pump::{HeatPump, HeatPumpHotWaterOnly};
use crate::core::heating_systems::point_of_use::PointOfUse;
use crate::core::heating_systems::storage_tank::{
    HeatSourceWithStorageTank, ImmersionHeater, SolarThermalSystem, StorageTank,
};
use crate::core::heating_systems::wwhrs::{
    WWHRSInstantaneousSystemA, WWHRSInstantaneousSystemB, WWHRSInstantaneousSystemC, Wwhrs,
};
use crate::core::material_properties::WATER;
use crate::core::schedule::{
    expand_boolean_schedule, expand_numeric_schedule, expand_water_heating_events, ScheduleEvent,
};
use crate::core::space_heat_demand::building_element::area_for_building_element_input;
use crate::core::space_heat_demand::internal_gains::{ApplianceGains, InternalGains};
use crate::core::space_heat_demand::thermal_bridge::{ThermalBridge, ThermalBridging};
use crate::core::space_heat_demand::ventilation_element::{
    air_change_rate_to_flow_rate, MechanicalVentilationHeatRecovery, NaturalVentilation,
    VentilationElement, VentilationElementInfiltration, WholeHouseExtractVentilation,
    WindowOpeningForCooling,
};
use crate::core::space_heat_demand::zone::{NamedBuildingElement, Zone};
use crate::core::units::{LITRES_PER_CUBIC_METRE, MILLIMETRES_IN_METRE};
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::dhw_demand::DomesticHotWaterDemand;
use crate::external_conditions::ExternalConditions;
use crate::input::HeatSourceWetType::HeatPump as HeatPumpInput;
use crate::input::{
    ApplianceGains as ApplianceGainsInput, ApplianceGainsDetails, BuildingElement,
    ColdWaterSourceDetails, ColdWaterSourceInput, ColdWaterSourceType, Control as ControlInput,
    ControlDetails, EnergyDiverter, EnergySupplyDetails, EnergySupplyInput, EnergySupplyType,
    ExternalConditionsInput, HeatSource as HeatSourceInput, HeatSourceControl,
    HeatSourceControlType, HeatSourceWetDetails, HeatSourceWetType, HotWaterSourceDetails,
    Infiltration, Input, InternalGains as InternalGainsInput, InternalGainsDetails,
    ThermalBridging as ThermalBridgingInput, ThermalBridgingDetails, Ventilation,
    WasteWaterHeatRecovery, WasteWaterHeatRecoveryDetails, WaterHeatingEvent, WaterHeatingEvents,
    WindowOpeningForCooling as WindowOpeningForCoolingInput, WwhrsType, ZoneDictionary, ZoneInput,
};
use crate::simulation_time::{SimulationTime, SimulationTimeIteration, SimulationTimeIterator};
use indexmap::IndexMap;
use serde_json::Value;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::sync::{Arc, Mutex};

// TODO make this a runtime parameter?
const DETAILED_OUTPUT_HEATING_COOLING: bool = true;

pub struct Corpus {
    pub external_conditions: Arc<ExternalConditions>,
    pub infiltration: VentilationElementInfiltration,
    pub cold_water_sources: ColdWaterSources,
    pub energy_supplies: EnergySupplies,
    pub internal_gains: InternalGainsCollection,
    pub controls: Controls,
    pub wwhrs: HashMap<String, Wwhrs>,
    pub event_schedules: HotWaterEventSchedules,
    pub domestic_hot_water_demand: DomesticHotWaterDemand,
    pub ventilation: Option<VentilationElement>,
    pub space_heating_ductwork: Option<Ductwork>,
    pub zones: HashMap<String, Zone>,
    pub heat_system_name_for_zone: HashMap<String, Option<String>>,
    pub cool_system_name_for_zone: HashMap<String, Option<String>>,
    pub total_floor_area: f64,
    pub total_volume: f64,
    pub wet_heat_sources: HashMap<String, Arc<WetHeatSource>>,
    pub hot_water_sources: HashMap<String, HotWaterSource>,
}

impl TryFrom<Input> for Corpus {
    type Error = ();

    fn try_from(input: Input) -> Result<Self, Self::Error> {
        let simulation_time_iterator = Arc::new(input.simulation_time.iter());

        let external_conditions = Arc::new(external_conditions_from_input(
            input.external_conditions.clone(),
            simulation_time_iterator.clone().as_ref(),
        ));

        let diverters: Diverters = (&input.energy_supply).into();

        let cold_water_sources =
            cold_water_sources_from_input(input.cold_water_source, &input.simulation_time);
        let wwhrs = wwhrs_from_input(input.waste_water_heat_recovery, &cold_water_sources);

        let mut energy_supplies = energy_supplies_from_input(
            input.energy_supply,
            simulation_time_iterator.clone().as_ref(),
        );

        let controls = control_from_input(input.control, simulation_time_iterator.clone().as_ref());

        let event_schedules = event_schedules_from_input(
            input.water_heating_events,
            simulation_time_iterator.as_ref(),
        );

        let domestic_hot_water_demand = DomesticHotWaterDemand::new(
            input.shower,
            input.bath,
            input.other_water_use,
            match &input.hot_water_source.hot_water_cylinder {
                HotWaterSourceDetails::PointOfUse { .. } => None,
                _ => input.water_distribution,
            },
            &cold_water_sources,
            &wwhrs,
            &energy_supplies,
            event_schedules.clone(),
        );

        let infiltration = infiltration_from_input(input.infiltration);

        let space_heating_ductwork = ductwork_from_ventilation_input(&input.ventilation);

        let ventilation = input.ventilation.as_ref().map(|v| {
            ventilation_from_input(&v, &infiltration, simulation_time_iterator.clone().as_ref())
        });

        let opening_area_total_from_zones = opening_area_total_from_zones(&input.zone);

        let mut heat_system_name_for_zone: HashMap<String, Option<String>> = Default::default();
        let mut cool_system_name_for_zone: HashMap<String, Option<String>> = Default::default();

        let zones: HashMap<String, Zone> = input
            .zone
            .iter()
            .map(|(i, zone)| {
                ((*i).clone(), {
                    let (zone_for_corpus, heat_system_name, cool_system_name) = zone_from_input(
                        zone,
                        opening_area_total_from_zones,
                        &input.window_opening_for_cooling,
                        &controls,
                        ventilation.as_ref(),
                        external_conditions.clone(),
                        &infiltration,
                        simulation_time_iterator.clone().as_ref(),
                    );
                    heat_system_name_for_zone.insert((*i).clone(), heat_system_name);
                    cool_system_name_for_zone.insert((*i).clone(), cool_system_name);

                    zone_for_corpus
                })
            })
            .collect();

        if !has_unique_some_values(&heat_system_name_for_zone)
            || !has_unique_some_values(&cool_system_name_for_zone)
        {
            return Err(());
        }

        // TODO: there needs to be some equivalent here of the Python code that builds the dict __energy_supply_conn_unmet_demand_zone

        let (total_floor_area, total_volume) = zones.values().fold((0., 0.), |acc, zone| {
            (zone.area() + acc.0, zone.volume() + acc.1)
        });

        let mut internal_gains = internal_gains_from_input(input.internal_gains);

        apply_appliance_gains_from_input(
            &mut internal_gains,
            input.appliance_gains,
            total_floor_area,
        );

        let wet_heat_sources: HashMap<String, Arc<WetHeatSource>> = input
            .heat_source_wet
            .unwrap_or_default()
            .iter()
            .map(|(name, heat_source_wet_details)| {
                (
                    (*name).clone(),
                    Arc::new(heat_source_wet_from_input(
                        (*heat_source_wet_details).clone(),
                        external_conditions.clone(),
                        simulation_time_iterator.clone(),
                        ventilation.as_ref(),
                        input.ventilation.as_ref().map(|v| v.req_ach()),
                        total_volume,
                        &controls,
                    )),
                )
            })
            .collect();

        let mut hot_water_sources: HashMap<String, HotWaterSource> = Default::default();
        hot_water_sources.insert(
            "hw cylinder".to_string(),
            hot_water_source_from_input(
                "hw cylinder".to_string(),
                input.hot_water_source.hot_water_cylinder,
                &cold_water_sources,
                &wet_heat_sources,
                &wwhrs,
                &controls,
                simulation_time_iterator.clone().as_ref(),
                external_conditions.clone(),
            ),
        );

        Ok(Self {
            external_conditions,
            infiltration,
            cold_water_sources,
            energy_supplies,
            internal_gains,
            controls,
            wwhrs,
            event_schedules,
            domestic_hot_water_demand,
            ventilation,
            space_heating_ductwork,
            zones,
            heat_system_name_for_zone,
            cool_system_name_for_zone,
            total_floor_area,
            total_volume,
            wet_heat_sources,
            hot_water_sources,
        })
    }
}

fn has_unique_some_values<K, V: Eq + Hash>(map: &HashMap<K, Option<V>>) -> bool {
    let some_values: Vec<&V> = map.values().flat_map(|v| v.iter()).collect();
    let value_set: HashSet<&&V> = some_values.iter().collect();
    some_values.len() == value_set.len()
}

fn external_conditions_from_input(
    input: Arc<ExternalConditionsInput>,
    simulation_time: &SimulationTimeIterator,
) -> ExternalConditions {
    ExternalConditions::new(
        simulation_time,
        input.air_temperatures.clone().unwrap_or_default(),
        input.wind_speeds.clone().unwrap_or_default(),
        input
            .diffuse_horizontal_radiation
            .clone()
            .unwrap_or_default(),
        input.direct_beam_radiation.clone().unwrap_or_default(),
        input
            .solar_reflectivity_of_ground
            .clone()
            .unwrap_or_default(),
        input.latitude.unwrap(),
        input.longitude.unwrap(),
        input.timezone.unwrap(),
        input.start_day.unwrap_or(0),
        input.end_day,
        input.time_series_step.unwrap_or(1.0),
        input.january_first,
        input.daylight_savings.clone().unwrap(),
        input.leap_day_included.unwrap_or(false),
        input.direct_beam_conversion_needed.unwrap_or(false),
        input.shading_segments.clone(),
    )
}

fn infiltration_from_input(input: Infiltration) -> VentilationElementInfiltration {
    VentilationElementInfiltration::new(
        input.storeys_in_building,
        input.shelter,
        input.build_type,
        input.test_result,
        input.test_type,
        input.env_area,
        input.volume,
        input.sheltered_sides,
        input.open_chimneys,
        input.open_flues,
        input.closed_fire,
        input.flues_d,
        input.flues_e,
        input.blocked_chimneys,
        input.extract_fans,
        input.passive_vents,
        input.gas_fires,
        input.storey_of_dwelling,
    )
}

pub struct ColdWaterSources {
    mains_water: Option<ColdWaterSource>,
    header_tank: Option<ColdWaterSource>,
}

impl ColdWaterSources {
    pub fn ref_for_mains_water(&self) -> Option<ColdWaterSource> {
        self.mains_water.clone()
    }

    pub fn ref_for_header_tank(&self) -> Option<ColdWaterSource> {
        self.header_tank.clone()
    }

    pub fn ref_for_type(&self, source_type: ColdWaterSourceType) -> Option<ColdWaterSource> {
        match source_type {
            ColdWaterSourceType::MainsWater => self.ref_for_mains_water(),
            ColdWaterSourceType::HeaderTank => self.ref_for_header_tank(),
        }
    }
}

fn cold_water_sources_from_input(
    input: ColdWaterSourceInput,
    simulation_time: &SimulationTime,
) -> ColdWaterSources {
    ColdWaterSources {
        mains_water: input
            .mains_water
            .map(|details| cold_water_source_from_input_details(details, simulation_time)),
        header_tank: input
            .header_tank
            .map(|details| cold_water_source_from_input_details(details, simulation_time)),
    }
}

fn cold_water_source_from_input_details(
    details: ColdWaterSourceDetails,
    simulation_time: &SimulationTime,
) -> ColdWaterSource {
    ColdWaterSource::new(
        details.temperatures,
        simulation_time,
        details.time_series_step,
    )
}

fn energy_supplies_from_input(
    input: EnergySupplyInput,
    simulation_time_iterator: &SimulationTimeIterator,
) -> EnergySupplies {
    EnergySupplies {
        mains_electricity: energy_supply_from_input(
            input.mains_electricity,
            simulation_time_iterator,
        ),
        mains_gas: energy_supply_from_input(input.mains_gas, simulation_time_iterator),
        bulk_lpg: energy_supply_from_input(input.bulk_lpg, simulation_time_iterator),
        heat_network: input.heat_network,
        unmet_demand: EnergySupply::new(
            EnergySupplyType::UnmetDemand,
            simulation_time_iterator.total_steps(),
            Default::default(),
        ),
    }
}

fn energy_supply_from_input(
    input: Option<EnergySupplyDetails>,
    simulation_time_iterator: &SimulationTimeIterator,
) -> Option<EnergySupply> {
    match input {
        Some(details) => Some(EnergySupply::new(
            details.fuel,
            simulation_time_iterator.total_steps(),
            details.electric_battery,
        )),
        None => None,
    }
}

struct Diverters {
    pub mains_electricity: Option<EnergyDiverter>,
    pub mains_gas: Option<EnergyDiverter>,
    pub bulk_lpg: Option<EnergyDiverter>,
}

impl From<&EnergySupplyInput> for Diverters {
    fn from(input: &EnergySupplyInput) -> Self {
        Self {
            mains_electricity: diverter_from_energy_supply(&input.mains_electricity),
            mains_gas: diverter_from_energy_supply(&input.mains_gas),
            bulk_lpg: diverter_from_energy_supply(&input.bulk_lpg),
        }
    }
}

fn diverter_from_energy_supply(supply: &Option<EnergySupplyDetails>) -> Option<EnergyDiverter> {
    supply
        .as_ref()
        .map(|supply| supply.diverter.clone().unwrap_or_default())
}

pub struct InternalGainsCollection {
    total_internal_gains: Option<InternalGains>,
    metabolic_gains: Option<InternalGains>,
    lighting: Option<ApplianceGains>,
    cooking: Option<ApplianceGains>,
    cooking1: Option<ApplianceGains>,
    cooking2: Option<ApplianceGains>,
    other: Option<InternalGains>,
}

fn internal_gains_from_input(input: InternalGainsInput) -> InternalGainsCollection {
    InternalGainsCollection {
        total_internal_gains: input.total_internal_gains.map(internal_gains_from_details),
        metabolic_gains: input.metabolic_gains.map(internal_gains_from_details),
        lighting: None,
        cooking: None,
        cooking1: None,
        cooking2: None,
        other: input.other.map(internal_gains_from_details),
    }
}

fn internal_gains_from_details(details: InternalGainsDetails) -> InternalGains {
    InternalGains::new(
        expand_numeric_schedule(
            HashMap::from([("main".to_string(), details.schedule.main)]),
            false,
        ),
        details.start_day,
        details.time_series_step,
    )
}

pub struct Controls {
    core: Vec<HeatSourceControl<Option<Arc<Control>>>>,
    extra: HashMap<String, Arc<Control>>,
}

impl Controls {
    pub fn new(
        core: Vec<HeatSourceControl<Option<Arc<Control>>>>,
        extra: HashMap<String, Arc<Control>>,
    ) -> Self {
        Self { core, extra }
    }

    pub fn get(&self, control_type: &HeatSourceControlType) -> Option<&Arc<Control>> {
        self.core
            .iter()
            .find(|heat_source_control| heat_source_control.get(control_type).is_some())
            .and_then(|heat_source_control| heat_source_control.get(control_type).unwrap().as_ref())
    }

    // access a control using a string, possibly because it is one of the "extra" controls
    pub fn get_with_string(&self, control_name: &str) -> Option<&Arc<Control>> {
        match control_name {
            // hard-code ways of resolving to core control types (for now)
            "hw timer" => self.get(&HeatSourceControlType::HotWaterTimer),
            "window opening" => self.get(&HeatSourceControlType::WindowOpening),
            other => self.extra.get(other),
        }
    }
}

fn control_from_input(
    control_input: ControlInput,
    simulation_time_iterator: &SimulationTimeIterator,
) -> Controls {
    let mut core: Vec<HeatSourceControl<Option<Arc<Control>>>> = Default::default();
    let mut extra: HashMap<String, Arc<Control>> = Default::default();

    // this is very ugly(!) but is just a reflection of the lack of clarity in the schema
    // and the way the variants-struct crate works;
    // we should be able to improve it in time
    for control in control_input.core {
        match control {
            HeatSourceControl {
                hot_water_timer: Some(control),
                ..
            } => {
                core.push(HeatSourceControl::new(
                    Some(Arc::new(single_control_from_details(
                        control,
                        simulation_time_iterator,
                    ))),
                    None,
                ));
            }
            HeatSourceControl {
                window_opening: Some(control),
                ..
            } => {
                core.push(HeatSourceControl::new(
                    None,
                    Some(Arc::new(single_control_from_details(
                        control,
                        simulation_time_iterator,
                    ))),
                ));
            }
            unknown => panic!(
                "incorrectly formed HeatSourceControl struct encountered: {:?}",
                unknown
            ),
        }
    }
    for (name, control) in control_input.extra {
        extra.insert(
            name,
            Arc::new(single_control_from_details(
                control,
                simulation_time_iterator,
            )),
        );
    }

    Controls { core, extra }
}

fn single_control_from_details(
    details: ControlDetails,
    simulation_time_iterator: &SimulationTimeIterator,
) -> Control {
    match details {
        ControlDetails::OnOffTime {
            start_day,
            time_series_step,
            schedule,
            ..
        } => Control::OnOffTimeControl(OnOffTimeControl::new(
            expand_boolean_schedule(schedule, false),
            start_day,
            time_series_step,
        )),
        ControlDetails::OnOffCostMinimisingTime {
            start_day,
            time_series_step,
            time_on_daily,
            schedule,
            ..
        } => Control::OnOffMinimisingTimeControl(OnOffMinimisingTimeControl::new(
            expand_numeric_schedule(schedule, false),
            start_day,
            time_series_step,
            time_on_daily.unwrap_or_default(),
        )),
        ControlDetails::SetpointTime {
            start_day,
            time_series_step,
            advanced_start,
            setpoint_min,
            setpoint_max,
            default_to_max,
            schedule,
            ..
        } => Control::SetpointTimeControl(
            SetpointTimeControl::new(
                expand_numeric_schedule(schedule, true)
                    .iter()
                    .map(|s| Some(*s))
                    .collect(),
                start_day,
                time_series_step,
                setpoint_min,
                setpoint_max,
                default_to_max,
                advanced_start,
                simulation_time_iterator.step_in_hours(),
            )
            .unwrap(),
        ),
        ControlDetails::ToUCharge {
            start_day,
            time_series_step,
            charge_level,
            schedule,
            ..
        } => {
            // Simulation manual charge control
            // Set charge level to 1.0 (max) for each day of simulation (plus 1)
            let vec_size = ((simulation_time_iterator.total_steps() as f64
                * simulation_time_iterator.step_in_hours()
                / 24.0)
                + 1.0)
                .ceil() as usize;
            let mut charge_level_vec: Vec<f64> = vec![1.0; vec_size];
            // if charge_level is present in input, overwrite initial vector
            // user can specify a vector with all days (plus 1), or as a single float value to be used for each day
            if let Some(charge) = charge_level {
                match charge {
                    Value::Array(charge_vec) => {
                        charge_level_vec = charge_vec.iter().map(|v| v.as_f64().unwrap()).collect();
                    }
                    Value::Number(charge) => {
                        charge_level_vec = vec![charge.as_f64().unwrap(); vec_size];
                    }
                    _ => {
                        panic!("Control charge value must be either a number of a list of numbers")
                    }
                }
            }

            Control::ToUChargeControl(ToUChargeControl {
                schedule: expand_boolean_schedule(schedule, false),
                start_day,
                time_series_step,
                charge_level: charge_level_vec,
            })
        }
    }
}

fn wwhrs_from_input(
    wwhrs: Option<WasteWaterHeatRecovery>,
    cold_water_sources: &ColdWaterSources,
) -> HashMap<String, Wwhrs> {
    let mut wwhr_systems: HashMap<String, Wwhrs> = HashMap::from([]);
    if let Some(systems) = wwhrs {
        for (name, system) in systems {
            wwhr_systems
                .entry(name)
                .or_insert(wwhr_system_from_details(system, cold_water_sources));
        }
    }

    wwhr_systems
}

fn wwhr_system_from_details(
    system: WasteWaterHeatRecoveryDetails,
    cold_water_sources: &ColdWaterSources,
) -> Wwhrs {
    match system.system_type {
        WwhrsType::SystemA => Wwhrs::WWHRSInstantaneousSystemA(WWHRSInstantaneousSystemA::new(
            system.flow_rates,
            system.efficiencies,
            get_cold_water_source_ref_for_type(system.cold_water_source, cold_water_sources)
                .unwrap(),
            system.utilisation_factor,
        )),
        WwhrsType::SystemB => Wwhrs::WWHRSInstantaneousSystemB(WWHRSInstantaneousSystemB::new(
            get_cold_water_source_ref_for_type(system.cold_water_source, cold_water_sources)
                .unwrap(),
            system.flow_rates,
            system.efficiencies,
            system.utilisation_factor,
        )),
        WwhrsType::SystemC => Wwhrs::WWHRSInstantaneousSystemC(WWHRSInstantaneousSystemC::new(
            system.flow_rates,
            system.efficiencies,
            get_cold_water_source_ref_for_type(system.cold_water_source, cold_water_sources)
                .unwrap(),
            system.utilisation_factor,
        )),
    }
}

fn get_cold_water_source_ref_for_type(
    source_type: ColdWaterSourceType,
    cold_water_sources: &ColdWaterSources,
) -> Option<ColdWaterSource> {
    match source_type {
        ColdWaterSourceType::MainsWater => cold_water_sources.ref_for_mains_water(),
        ColdWaterSourceType::HeaderTank => cold_water_sources.ref_for_header_tank(),
    }
}

pub type EventSchedule = Vec<Option<Vec<ScheduleEvent>>>;

#[derive(Clone)]
pub struct HotWaterEventSchedules {
    pub shower: HashMap<String, EventSchedule>,
    pub bath: HashMap<String, EventSchedule>,
    pub other: HashMap<String, EventSchedule>,
}

fn event_schedules_from_input(
    events: WaterHeatingEvents,
    simulation_time_iterator: &SimulationTimeIterator,
) -> HotWaterEventSchedules {
    let mut shower_schedules: HashMap<String, EventSchedule> = Default::default();
    if let Some(shower_events) = events.shower {
        shower_schedules.insert(
            "ies".to_string(),
            schedule_event_from_input(shower_events.ies.iter().collect(), simulation_time_iterator),
        );
        shower_schedules.insert(
            "mixer".to_string(),
            schedule_event_from_input(
                shower_events.mixer.iter().collect(),
                simulation_time_iterator,
            ),
        );
    }

    let mut bath_schedules: HashMap<String, EventSchedule> = Default::default();
    if let Some(bath_events) = events.bath {
        bath_schedules.insert(
            "medium".to_string(),
            schedule_event_from_input(
                bath_events.medium.iter().collect(),
                simulation_time_iterator,
            ),
        );
    }

    let mut other_schedules: HashMap<String, EventSchedule> = Default::default();
    if let Some(other_events) = events.other {
        other_schedules.insert(
            "other".to_string(),
            schedule_event_from_input(
                other_events.other.iter().collect(),
                simulation_time_iterator,
            ),
        );
    }

    HotWaterEventSchedules {
        shower: shower_schedules,
        bath: bath_schedules,
        other: other_schedules,
    }
}

fn schedule_event_from_input(
    events_input: Vec<&WaterHeatingEvent>,
    simulation_time_iterator: &SimulationTimeIterator,
) -> EventSchedule {
    let sim_timestep = simulation_time_iterator.step_in_hours();
    let total_timesteps = simulation_time_iterator.total_steps();
    expand_water_heating_events(events_input, sim_timestep, total_timesteps)
}

fn ductwork_from_ventilation_input(ventilation: &Option<Ventilation>) -> Option<Ductwork> {
    ventilation.as_ref().and_then(|v| match v {
        Ventilation::Mvhr { ductwork, .. } => Some(Ductwork::new(
            ductwork.internal_diameter_mm / MILLIMETRES_IN_METRE as f64,
            ductwork.external_diameter_mm / MILLIMETRES_IN_METRE as f64,
            ductwork.length_in,
            ductwork.length_out,
            ductwork.insulation_thermal_conductivity,
            ductwork.insulation_thickness_mm / MILLIMETRES_IN_METRE as f64,
            ductwork.reflective,
            ductwork.mvhr_location,
        )),
        _ => None,
    })
}

fn ventilation_from_input<'a>(
    ventilation: &'a Ventilation,
    infiltration: &VentilationElementInfiltration,
    simulation_time: &'a SimulationTimeIterator,
) -> VentilationElement {
    match ventilation {
        Ventilation::Whev { req_ach, sfp, .. } => {
            VentilationElement::Whev(WholeHouseExtractVentilation::new(
                *req_ach,
                *sfp,
                infiltration.infiltration_rate(),
                // energy_supply_from_type_for_ventilation(energy_supply, energy_supplies),
                // "Ventilation system".to_string(),
                simulation_time.step_in_hours(),
            ))
        }
        Ventilation::Mvhr {
            req_ach,
            sfp,
            efficiency,
            ..
        } => VentilationElement::Mvhr(MechanicalVentilationHeatRecovery::new(
            *req_ach,
            *sfp,
            *efficiency,
            // energy_supply_from_type_for_ventilation(energy_supply, energy_supplies),
            // "Ventilation system".to_string(),
            simulation_time.step_in_hours(),
        )),
        Ventilation::Natural { req_ach } => VentilationElement::Natural(NaturalVentilation::new(
            *req_ach,
            infiltration.infiltration_rate(),
        )),
    }
}

fn energy_supply_from_type_for_ventilation<'a>(
    energy_supply_type: &'a EnergySupplyType,
    energy_supplies: &'a mut EnergySupplies,
) -> &'a mut EnergySupply {
    match energy_supply_type {
        EnergySupplyType::Electricity => energy_supplies.mains_electricity.as_mut().unwrap(),
        EnergySupplyType::MainsGas => energy_supplies.mains_gas.as_mut().unwrap(),
        EnergySupplyType::UnmetDemand => &mut energy_supplies.unmet_demand,
        EnergySupplyType::LpgBulk => energy_supplies.bulk_lpg.as_mut().unwrap(),
        // commenting out for now as this is a different type and might not be used in this context
        // EnergySupplyType::HeatNetwork => energy_supplies.heat_network.unwrap(),
        _ => panic!("Unexpected energy supply type listed for ventilation."),
    }
}

fn opening_area_total_from_zones(zones: &ZoneDictionary) -> f64 {
    zones
        .iter()
        .flat_map(|(_, zone)| {
            zone.building_elements
                .iter()
                .map(|(_, building_element)| match building_element {
                    BuildingElement::Transparent { height, width, .. } => height * width,
                    _ => 0.,
                })
        })
        .sum()
}

fn check_space_heat_systems_unique_to_zones(zones: &ZoneDictionary) -> Result<(), &'static str> {
    let res = {
        let mut name_set = HashSet::new();
        zones.iter().all(|(_, zone)| match &zone.space_heat_system {
            Some(system_name) => name_set.insert(system_name),
            None => true,
        })
    };
    if res {
        Ok(())
    } else {
        Err("A space heat system was declared as used in more than one zone.")
    }
}

fn check_space_cool_systems_unique_to_zones(zones: &ZoneDictionary) -> Result<(), &'static str> {
    let res = {
        let mut name_set = HashSet::new();
        zones.iter().all(|(_, zone)| match &zone.space_cool_system {
            Some(system_name) => name_set.insert(system_name),
            None => true,
        })
    };
    if res {
        Ok(())
    } else {
        Err("A space cool system was declared as used in more than one zone.")
    }
}

fn thermal_bridging_from_input(input: &ThermalBridgingInput) -> ThermalBridging {
    match input {
        ThermalBridgingInput::ThermalBridgingElements(input_bridges) => ThermalBridging::Bridges({
            let mut bridges = IndexMap::new();
            bridges.extend(input_bridges.iter().map(|(name, details)| {
                (
                    name.clone(),
                    match details {
                        ThermalBridgingDetails::Linear {
                            linear_thermal_transmittance,
                            length,
                        } => ThermalBridge::Linear {
                            linear_thermal_transmittance: *linear_thermal_transmittance,
                            length: *length,
                        },
                        ThermalBridgingDetails::Point {
                            heat_transfer_coefficient,
                        } => ThermalBridge::Point {
                            heat_transfer_coefficient: *heat_transfer_coefficient,
                        },
                    },
                )
            }));
            bridges
        }),
        ThermalBridgingInput::ThermalBridgingNumber(num) => ThermalBridging::Number(*num),
    }
}

fn zone_from_input<'a>(
    input: &ZoneInput,
    opening_area_total: f64,
    window_opening_for_cooling: &Option<WindowOpeningForCoolingInput>,
    controls: &'a Controls,
    ventilation: Option<&'a VentilationElement>,
    external_conditions: Arc<ExternalConditions>,
    infiltration: &'a VentilationElementInfiltration,
    simulation_time_iterator: &'a SimulationTimeIterator,
) -> (Zone, Option<String>, Option<String>) {
    let heat_system_name = input.space_heat_system.clone();
    let cool_system_name = input.space_cool_system.clone();

    let vent_cool_extra = window_opening_for_cooling.as_ref().map(|opening| {
        let openings: HashMap<&String, &BuildingElement> = input
            .building_elements
            .iter()
            .filter(|(_, el)| matches!(el, BuildingElement::Transparent { .. }))
            .collect();
        let opening_area_zone: f64 = openings
            .values()
            .map(|op| area_for_building_element_input(op))
            .sum();
        let opening_area_equivalent =
            opening.equivalent_area * opening_area_zone / opening_area_total;
        let control = input
            .control_window_opening
            .as_ref()
            .and_then(|opening_control| {
                controls
                    .get(opening_control)
                    .and_then(|c| match c.as_ref() {
                        Control::SetpointTimeControl(ctrl) => Some((*ctrl).clone()),
                        _ => None,
                    })
            });
        let natvent = ventilation.and_then(|v| match v {
            VentilationElement::Natural(natural_ventilation) => {
                Some((*natural_ventilation).clone())
            }
            _ => None,
        });
        let named_openings: Vec<NamedBuildingElement> = openings
            .iter()
            .map(|(name, element)| NamedBuildingElement {
                name: name.to_string(),
                element: (*element).clone(),
            })
            .collect::<Vec<_>>();
        WindowOpeningForCooling::new(
            opening_area_equivalent,
            external_conditions.clone(),
            named_openings,
            control,
            natvent,
        )
    });

    let infiltration_ventilation = VentilationElement::Infiltration((*infiltration).clone());
    let mut vent_elements: Vec<VentilationElement> = vec![infiltration_ventilation];

    if let Some(v) = ventilation {
        vent_elements.push((*v).clone());
    }

    (
        Zone::new(
            input.area,
            input.volume,
            input.building_elements.clone(),
            thermal_bridging_from_input(&input.thermal_bridging),
            vent_elements,
            vent_cool_extra,
            external_conditions.air_temp_for_timestep_idx(simulation_time_iterator.current_index()),
            input.temp_setpnt_init.unwrap(),
            external_conditions.clone(),
            simulation_time_iterator,
        ),
        heat_system_name,
        cool_system_name,
    )
}

fn apply_appliance_gains_from_input(
    internal_gains_collection: &mut InternalGainsCollection,
    input: ApplianceGainsInput,
    total_floor_area: f64,
) {
    if let Some(details) = input.lighting {
        internal_gains_collection.lighting = Some(appliance_gains_from_single_input(
            details,
            "lighting".to_string(),
            total_floor_area,
        ));
    }
    if let Some(details) = input.cooking {
        internal_gains_collection.cooking = Some(appliance_gains_from_single_input(
            details,
            "cooking".to_string(),
            total_floor_area,
        ));
    }
    if let Some(details) = input.cooking1 {
        internal_gains_collection.cooking1 = Some(appliance_gains_from_single_input(
            details,
            "cooking1".to_string(),
            total_floor_area,
        ));
    }
    if let Some(details) = input.cooking2 {
        internal_gains_collection.cooking2 = Some(appliance_gains_from_single_input(
            details,
            "cooking2".to_string(),
            total_floor_area,
        ));
    }
}

fn appliance_gains_from_single_input(
    input: ApplianceGainsDetails,
    supply_end_user_name: String,
    total_floor_area: f64,
) -> ApplianceGains {
    let total_energy_supply = expand_numeric_schedule(input.schedule, false)
        .iter()
        .map(|energy_data| energy_data / total_floor_area)
        .collect();

    ApplianceGains::new(
        total_energy_supply,
        supply_end_user_name,
        input.gains_fraction,
        input.start_day,
        input.time_series_step,
    )
}

#[derive(Clone)]
pub enum HeatSource {
    Storage(HeatSourceWithStorageTank),
    Wet(Box<HeatSourceWet>),
}

impl HeatSource {
    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        match self {
            HeatSource::Storage(ref mut storage) => match storage {
                HeatSourceWithStorageTank::Immersion(imm) => imm
                    .lock()
                    .expect("Expected to obtain a mutex lock on an immersion heater value")
                    .demand_energy(energy_demand, simulation_time_iteration.index),
                HeatSourceWithStorageTank::Solar(ref mut solar) => {
                    solar.demand_energy(energy_demand)
                }
            },
            HeatSource::Wet(ref mut wet) => match wet.as_mut() {
                HeatSourceWet::WaterCombi(_) => {
                    panic!("not expected? this value does not have a demand_energy method")
                    // the Python uses duck-typing here but there is no method for this type
                }
                HeatSourceWet::WaterRegular(ref mut r) => {
                    r.demand_energy(energy_demand, simulation_time_iteration.index)
                }
                HeatSourceWet::Space(_) => {
                    panic!("not expected? this value does not have a demand_energy method")
                    // the Python uses duck-typing here but there is no method for this type
                }
                HeatSourceWet::HeatNetworkWaterStorage(ref mut h) => {
                    h.demand_energy(energy_demand, simulation_time_iteration)
                }
                HeatSourceWet::HeatBatteryHotWater(ref mut h) => {
                    h.demand_energy(energy_demand, *simulation_time_iteration)
                }
                HeatSourceWet::HeatPumpWater(ref mut h) => {
                    h.demand_energy(energy_demand, simulation_time_iteration)
                }
                HeatSourceWet::HeatPumpWaterOnly(h) => {
                    h.demand_energy(energy_demand, simulation_time_iteration.index)
                }
            },
        }
    }
}

#[derive(Clone)]
pub struct PositionedHeatSource {
    pub heat_source: HeatSource,
    pub heater_position: f64,
    pub thermostat_position: f64,
}

pub enum WetHeatSource {
    HeatPump(HeatPump),
    Boiler(Boiler),
    Hiu(HeatNetwork),
    HeatBattery(HeatBattery), // TODO to be implemented
}

fn heat_source_wet_from_input(
    input: HeatSourceWetDetails,
    external_conditions: Arc<ExternalConditions>,
    simulation_time: Arc<SimulationTimeIterator>,
    ventilation: Option<&VentilationElement>,
    ventilation_req_ach: Option<f64>,
    total_volume: f64,
    controls: &Controls,
) -> WetHeatSource {
    match &input {
        HeatSourceWetDetails::HeatPump { source_type, .. } => {
            let throughput_exhaust_air = if source_type.is_exhaust_air() {
                // Check that ventilation system is compatible with exhaust air HP
                if ventilation.is_none()
                    || !matches!(
                        ventilation.unwrap(),
                        VentilationElement::Mvhr(_) | VentilationElement::Whev(_)
                    )
                {
                    panic!("Exhaust air heat pump requires ventilation to be MVHR or WHEV.")
                }
                Some(
                    air_change_rate_to_flow_rate(ventilation_req_ach.unwrap(), total_volume)
                        * LITRES_PER_CUBIC_METRE as f64,
                )
            } else {
                None
            };

            WetHeatSource::HeatPump(
                HeatPump::new(
                    &input,
                    simulation_time.step_in_hours(),
                    external_conditions.clone(),
                    throughput_exhaust_air,
                    DETAILED_OUTPUT_HEATING_COOLING,
                )
                .unwrap(),
            )
        }
        HeatSourceWetDetails::Boiler { .. } => WetHeatSource::Boiler(
            Boiler::new(
                input,
                external_conditions.clone(),
                simulation_time.step_in_hours(),
            )
            .expect("could not construct boiler value from provided data"),
        ),
        HeatSourceWetDetails::Hiu {
            power_max,
            hiu_daily_loss,
            building_level_distribution_losses,
            ..
        } => WetHeatSource::Hiu(HeatNetwork::new(
            *power_max,
            *hiu_daily_loss,
            *building_level_distribution_losses,
            simulation_time.step_in_hours(),
        )),
        HeatSourceWetDetails::HeatBattery { control_charge, .. } => {
            let heat_source = WetHeatSource::HeatBattery(HeatBattery::new(
                &input,
                controls
                    .get_with_string(control_charge)
                    .unwrap_or_else(|| {
                        panic!(
                            "expected a control to be registered with the name '{control_charge}'"
                        )
                    })
                    .clone(),
                simulation_time,
                external_conditions.clone(),
            ));
            heat_source
        }
    }
}

fn heat_source_from_input(
    input: HeatSourceInput,
    temp_setpoint: f64,
    wet_heat_sources: &HashMap<String, Arc<WetHeatSource>>,
    simulation_time: &SimulationTimeIterator,
    controls: &Controls,
    cold_water_sources: &ColdWaterSources,
    external_conditions: Arc<ExternalConditions>,
) -> HeatSource {
    // TODO add in all the stuff to do with energy supply

    match input {
        HeatSourceInput::ImmersionHeater { power, control, .. } => HeatSource::Storage(
            HeatSourceWithStorageTank::Immersion(Arc::new(Mutex::new(ImmersionHeater::new(
                power,
                simulation_time.step_in_hours(),
                control.and_then(|ctrl| controls.get(&ctrl).map(|c| (*c).clone())),
            )))),
        ),
        HeatSourceInput::SolarThermalSystem {
            solar_cell_location,
            area_module,
            modules,
            peak_collector_efficiency,
            incidence_angle_modifier,
            first_order_hlc,
            second_order_hlc,
            collector_mass_flow_rate,
            power_pump,
            power_pump_control,
            tilt,
            orientation,
            solar_loop_piping_hlc,
            ..
        } => HeatSource::Storage(HeatSourceWithStorageTank::Solar(SolarThermalSystem::new(
            solar_cell_location,
            area_module,
            modules,
            peak_collector_efficiency,
            incidence_angle_modifier,
            first_order_hlc,
            second_order_hlc,
            collector_mass_flow_rate,
            power_pump,
            power_pump_control,
            tilt,
            orientation,
            solar_loop_piping_hlc,
            external_conditions.clone(),
            simulation_time.step_in_hours(),
            WATER.clone(),
        ))),
        HeatSourceInput::Wet {
            name,
            cold_water_source: cold_water_source_type,
            control,
            temp_flow_limit_upper,
            ..
        } => {
            let cold_water_source = cold_water_sources
                .ref_for_type(cold_water_source_type)
                .expect("Expected a cold water source to be available to a boiler heat source.");
            let energy_supply_conn_name = format!("{name}_water_heating");
            let heat_source_wet = wet_heat_sources
                .get(&name)
                .unwrap_or_else(|| {
                    panic!("Expected a wet heat source registered with the name '{name}'.")
                })
                .clone();
            let source_control = control.and_then(|ctrl| controls.get(&ctrl).map(|c| (*c).clone()));

            match heat_source_wet.as_ref() {
                WetHeatSource::HeatPump(heat_pump) => HeatSource::Wet(Box::new(
                    HeatSourceWet::HeatPumpWater(HeatPump::create_service_hot_water(
                        Arc::new(Mutex::new((*heat_pump).clone())),
                        energy_supply_conn_name,
                        temp_setpoint,
                        55.,
                        temp_flow_limit_upper
                            .expect("temp_flow_limit_upper field was expected to be set"),
                        Arc::new(cold_water_source),
                        source_control,
                    )),
                )),
                WetHeatSource::Boiler(boiler) => HeatSource::Wet(Box::new(
                    HeatSourceWet::WaterRegular(boiler.create_service_hot_water_regular(
                        energy_supply_conn_name,
                        temp_setpoint,
                        cold_water_source,
                        55.,
                        source_control,
                    )),
                )),
                WetHeatSource::Hiu(heat_network) => {
                    HeatSource::Wet(Box::new(HeatSourceWet::HeatNetworkWaterStorage(
                        HeatNetwork::create_service_hot_water_storage(
                            Arc::new(Mutex::new((*heat_network).clone())),
                            energy_supply_conn_name,
                            temp_setpoint,
                            source_control,
                        ),
                    )))
                }
                WetHeatSource::HeatBattery(battery) => {
                    HeatSource::Wet(Box::new(HeatSourceWet::HeatBatteryHotWater(
                        HeatBattery::create_service_hot_water_regular(
                            Arc::new(Mutex::new((*battery).clone())),
                            &energy_supply_conn_name,
                            temp_setpoint,
                            Arc::new(cold_water_source),
                            55.,
                            source_control,
                        ),
                    )))
                }
            }
        }
        HeatSourceInput::HeatPumpHotWaterOnly {
            power_max,
            vol_hw_daily_average,
            ref test_data,
            energy_supply,
            control,
            heater_position,
            thermostat_position,
        } => HeatSource::Wet(Box::new(HeatSourceWet::HeatPumpWaterOnly(
            HeatPumpHotWaterOnly::new(
                power_max,
                &test_data,
                vol_hw_daily_average,
                simulation_time.step_in_hours(),
                controls.get(&control).map(|c| (*c).clone()),
            ),
        ))),
    }
}

enum HotWaterSource {
    StorageTank(StorageTank),
    CombiBoiler(BoilerServiceWaterCombi),
    PointOfUse(PointOfUse),
    HeatNetwork(HeatNetworkServiceWaterDirect),
    HeatBattery(()),
}

fn hot_water_source_from_input(
    source_name: String,
    input: HotWaterSourceDetails,
    cold_water_sources: &ColdWaterSources,
    wet_heat_sources: &HashMap<String, Arc<WetHeatSource>>,
    wwhrs: &HashMap<String, Wwhrs>,
    controls: &Controls,
    simulation_time: &SimulationTimeIterator,
    external_conditions: Arc<ExternalConditions>,
) -> HotWaterSource {
    let cloned_input = input.clone();
    match input {
        HotWaterSourceDetails::StorageTank {
            volume,
            daily_losses,
            min_temp,
            setpoint_temp,
            control_hold_at_setpoint,
            cold_water_source: cold_water_source_type,
            primary_pipework,
            heat_source,
        } => {
            let mut cold_water_source: WaterSourceWithTemperature =
                cold_water_source_for_type(cold_water_source_type, cold_water_sources);
            if !wwhrs.is_empty() {
                for heat_recovery_system in wwhrs.values() {
                    match heat_recovery_system {
                        Wwhrs::WWHRSInstantaneousSystemC(c) => {
                            cold_water_source =
                                WaterSourceWithTemperature::WwhrsC(Arc::new((*c).clone()));
                        }
                        Wwhrs::WWHRSInstantaneousSystemA(a) => {
                            cold_water_source =
                                WaterSourceWithTemperature::WwhrsA(Arc::new((*a).clone()));
                        }
                        _ => {}
                    }
                }
            }
            let pipework = primary_pipework.and_then(|p| p.into());
            let mut heat_sources: IndexMap<String, PositionedHeatSource> = Default::default();
            for (name, hs) in heat_source {
                let heater_position = hs.heater_position();
                let thermostat_position = hs.thermostat_position();
                heat_sources.insert(
                    name,
                    PositionedHeatSource {
                        heat_source: heat_source_from_input(
                            hs,
                            setpoint_temp,
                            wet_heat_sources,
                            simulation_time,
                            controls,
                            cold_water_sources,
                            external_conditions.clone(),
                        ),
                        heater_position,
                        thermostat_position,
                    },
                );
            }
            let ctrl_hold_at_setpoint = control_hold_at_setpoint
                .and_then(|ctrl| controls.get_with_string(ctrl.as_str()).cloned());
            HotWaterSource::StorageTank(StorageTank::new(
                volume,
                daily_losses,
                min_temp,
                setpoint_temp,
                cold_water_source,
                simulation_time.step_in_hours(),
                heat_sources,
                pipework,
                ctrl_hold_at_setpoint,
                WATER.clone(),
            ))
            // TODO add diverters stuff
        }
        HotWaterSourceDetails::CombiBoiler {
            cold_water_source: cold_water_source_type,
            heat_source_wet: heat_source_wet_type,
            ..
        } => {
            let cold_water_source =
                cold_water_source_for_type(cold_water_source_type, cold_water_sources);
            let energy_supply_conn_name = "boiler_water_heating".to_string(); // making assumption wet heat source is boiler, as this is only one allowable
            let heat_source_wet = match heat_source_wet_type {
                HeatSourceWetType::Boiler => {
                    match wet_heat_sources
                        .get("boiler")
                        .expect("Expected a boiler as wet heat source")
                        .as_ref()
                    {
                        WetHeatSource::Boiler(boiler) => boiler,
                        _ => panic!("Did not expect a heat source type that was not a boiler"),
                    }
                }
                _ => panic!("Did not expect a heat source type that was not a boiler"),
            };
            HotWaterSource::CombiBoiler(
                heat_source_wet
                    .create_service_hot_water_combi(
                        cloned_input,
                        energy_supply_conn_name,
                        60.,
                        cold_water_source,
                    )
                    .expect("expected to be able to instantiate a combi boiler object"),
            )
        }
        HotWaterSourceDetails::PointOfUse {
            power,
            efficiency,
            cold_water_source: cold_water_source_type,
            ..
        } => {
            let _energy_supply_conn_name = source_name;
            // TODO energy supply stuff
            let cold_water_source =
                cold_water_source_for_type(cold_water_source_type, cold_water_sources);
            HotWaterSource::PointOfUse(PointOfUse::new(power, efficiency, cold_water_source))
        }
        HotWaterSourceDetails::Hiu {
            cold_water_source: cold_water_source_type,
            heat_source_wet: heat_source_wet_type,
            control,
        } => {
            let energy_supply_conn_name = source_name;
            let cold_water_source =
                cold_water_source_for_type(cold_water_source_type, cold_water_sources);
            let heat_source_wet = match heat_source_wet_type {
                HeatSourceWetType::HeatNetwork => {
                    match wet_heat_sources
                        .get("HeatNetwork")
                        .expect("expected a heat network in this context")
                        .as_ref()
                    {
                        WetHeatSource::Hiu(heat_network) => heat_network,
                        _ => panic!("expected a heat network in this context"),
                    }
                }
                _ => panic!("expected a heat network in this context"),
            };
            HotWaterSource::HeatNetwork(HeatNetwork::create_service_hot_water_direct(
                Arc::new(Mutex::new((*heat_source_wet).clone())),
                energy_supply_conn_name,
                60.,
                cold_water_source,
            ))
        }
        HotWaterSourceDetails::HeatBattery { .. } => todo!(), // TODO is from Python
    }
}

fn cold_water_source_for_type(
    cold_water_source_type: ColdWaterSourceType,
    cold_water_sources: &ColdWaterSources,
) -> WaterSourceWithTemperature {
    WaterSourceWithTemperature::ColdWaterSource(Arc::new(match cold_water_source_type {
        ColdWaterSourceType::MainsWater => cold_water_sources
            .mains_water
            .as_ref()
            .expect("referenced cold water source was expected to exist")
            .clone(),
        ColdWaterSourceType::HeaderTank => cold_water_sources
            .header_tank
            .as_ref()
            .expect("referenced cold water source was expected to exist")
            .clone(),
    }))
}
