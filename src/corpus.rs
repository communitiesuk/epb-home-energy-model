use crate::core::controls::time_control::{
    Control, OnOffMinimisingTimeControl, OnOffTimeControl, SetpointTimeControl, ToUChargeControl,
};
use crate::core::energy_supply::energy_supply::{EnergySupplies, EnergySupply};
use crate::core::heating_systems::wwhrs::{
    WWHRSInstantaneousSystemA, WWHRSInstantaneousSystemB, WWHRSInstantaneousSystemC, Wwhrs,
};
use crate::core::schedule::{
    expand_boolean_schedule, expand_events, expand_numeric_schedule, expand_water_heating_events,
    ScheduleEvent,
};
use crate::core::space_heat_demand::internal_gains::InternalGains;
use crate::core::space_heat_demand::ventilation_element::VentilationElementInfiltration;
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::external_conditions::ExternalConditions;
use crate::input::{
    ColdWaterSourceDetails, ColdWaterSourceInput, ColdWaterSourceType, Control as ControlInput,
    ControlDetails, EnergyDiverter, EnergySupplyDetails, EnergySupplyInput, EnergySupplyType,
    ExternalConditionsInput, Infiltration, Input, InternalGains as InternalGainsInput,
    InternalGainsDetails, WasteWaterHeatRecovery, WasteWaterHeatRecoveryDetails, WaterHeatingEvent,
    WaterHeatingEvents, WwhrsType,
};
use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
use serde_json::Value;
use std::collections::HashMap;
use std::sync::Arc;

pub struct Corpus {
    pub external_conditions: ExternalConditions,
    pub infiltration: VentilationElementInfiltration,
    pub cold_water_sources: ColdWaterSources,
    pub energy_supplies: EnergySupplies,
    pub internal_gains: InternalGainsCollection,
    pub control: HashMap<String, Control>,
    pub wwhrs: HashMap<String, Wwhrs>,
    pub event_schedules: HotWaterEventSchedules,
}

impl From<Input> for Corpus {
    fn from(input: Input) -> Self {
        let simulation_time_iterator = Arc::new(input.simulation_time.iter());

        let external_conditions = external_conditions_from_input(
            input.external_conditions.clone(),
            simulation_time_iterator.clone().as_ref(),
        );

        let diverters: Diverters = (&input.energy_supply).into();

        let cold_water_sources =
            cold_water_sources_from_input(input.cold_water_source, &input.simulation_time);
        let wwhrs = wwhrs_from_input(input.waste_water_heat_recovery, &cold_water_sources);

        let energy_supplies = energy_supplies_from_input(
            input.energy_supply,
            simulation_time_iterator.clone().as_ref(),
        );

        Self {
            external_conditions: external_conditions,
            infiltration: infiltration_from_input(input.infiltration),
            cold_water_sources,
            energy_supplies,
            internal_gains: internal_gains_from_input(input.internal_gains),
            control: control_from_input(input.control, simulation_time_iterator.clone().as_ref()),
            wwhrs,
            event_schedules: event_schedules_from_input(
                input.water_heating_events,
                simulation_time_iterator.as_ref(),
            ),
        }
    }
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
}

fn cold_water_sources_from_input(
    input: ColdWaterSourceInput,
    simulation_time: &SimulationTime,
) -> ColdWaterSources {
    ColdWaterSources {
        mains_water: match input.mains_water {
            Some(details) => Some(cold_water_source_from_input_details(
                details,
                simulation_time,
            )),
            None => None,
        },
        header_tank: match input.header_tank {
            Some(details) => Some(cold_water_source_from_input_details(
                details,
                simulation_time,
            )),
            None => None,
        },
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
    match supply {
        Some(supply) => Some(supply.diverter.clone().unwrap_or_default()),
        None => None,
    }
}

pub struct InternalGainsCollection {
    total_internal_gains: Option<InternalGains>,
    metabolic_gains: Option<InternalGains>,
    other: Option<InternalGains>,
}

fn internal_gains_from_input(input: InternalGainsInput) -> InternalGainsCollection {
    InternalGainsCollection {
        total_internal_gains: match input.total_internal_gains {
            Some(details) => Some(internal_gains_from_details(details)),
            None => None,
        },
        metabolic_gains: match input.metabolic_gains {
            Some(details) => Some(internal_gains_from_details(details)),
            None => None,
        },
        other: match input.other {
            Some(details) => Some(internal_gains_from_details(details)),
            None => None,
        },
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

fn control_from_input(
    control: ControlInput,
    simulation_time_iterator: &SimulationTimeIterator,
) -> HashMap<String, Control> {
    let mut controls: Vec<(String, Control)> = vec![];

    for (name, control) in control {
        controls.push((
            name,
            match control {
                ControlDetails::OnOffTimeControl {
                    start_day,
                    time_series_step,
                    schedule,
                    ..
                } => Control::OnOffTimeControl(OnOffTimeControl::new(
                    expand_boolean_schedule(schedule, false),
                    start_day,
                    time_series_step,
                )),
                ControlDetails::OnOffCostMinimisingTimeControl {
                    start_day,
                    time_series_step,
                    advanced_start,
                    logic_type,
                    time_on_daily,
                    schedule,
                } => Control::OnOffMinimisingTimeControl(OnOffMinimisingTimeControl::new(
                    expand_numeric_schedule(schedule, false),
                    start_day,
                    time_series_step,
                    time_on_daily.unwrap_or_default(),
                )),
                ControlDetails::SetpointTimeControl {
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
                ControlDetails::ToUChargeControl {
                    start_day,
                    time_series_step,
                    advanced_start,
                    logic_type,
                    charge_level,
                    target_charge,
                    schedule,
                } => {
                    // Simulation manual charge control
                    // Set charge level to 1.0 (max) for each day of simulation (plus 1)
                    let vec_size = ((simulation_time_iterator.total_steps() as f64
                        * simulation_time_iterator.step_in_hours()
                        / 24.0)
                        + 1.0)
                        .ceil() as usize;
                    let charge_level_vec: Vec<f64> = vec![1.0; vec_size];
                    // if charge_level is present in input, overwrite initial vector
                    // user can specify a vector with all days (plus 1), or as a single float value to be used for each day
                    if let Some(charge) = charge_level {
                        match charge {
                            Value::Array(charge_vec) => {
                                let charge_level_vec: Vec<f64> =
                                    charge_vec.iter().map(|v| v.as_f64().unwrap()).collect();
                            }
                            Value::Number(charge) => {
                                let charge_level_vec = vec![charge.as_f64().unwrap(); vec_size];
                            }
                            _ => panic!(
                                "Control charge value must be either a number of a list of numbers"
                            ),
                        }
                    }

                    Control::ToUChargeControl(ToUChargeControl {
                        schedule: expand_boolean_schedule(schedule, false),
                        start_day,
                        time_series_step,
                        charge_level: charge_level_vec,
                    })
                }
            },
        ));
    }

    controls.into_iter().collect::<HashMap<_, _>>()
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

fn get_cold_water_source_ref_for_type<'a>(
    source_type: ColdWaterSourceType,
    cold_water_sources: &ColdWaterSources,
) -> Option<ColdWaterSource> {
    match source_type {
        ColdWaterSourceType::MainsWater => cold_water_sources.ref_for_mains_water(),
        ColdWaterSourceType::HeaderTank => cold_water_sources.ref_for_header_tank(),
    }
}

pub type EventSchedule = Vec<Option<Vec<ScheduleEvent>>>;

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
