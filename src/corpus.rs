use crate::core::common::WaterSourceWithTemperature;
use crate::core::controls::time_control::{
    ChargeControl, CombinationTimeControl, Control, ControlBehaviour, HeatSourceControl,
    OnOffMinimisingTimeControl, OnOffTimeControl, SetpointTimeControl, SmartApplianceControl,
};
use crate::core::cooling_systems::air_conditioning::AirConditioning;
use crate::core::energy_supply::elec_battery::ElectricBattery;
use crate::core::energy_supply::energy_supply::{
    EnergySupply, EnergySupplyBuilder, EnergySupplyConnection, EnergySupplyTariffInput,
    ENERGY_FROM_ENVIRONMENT_SUPPLY_NAME, UNMET_DEMAND_SUPPLY_NAME,
};
use crate::core::energy_supply::inverter::Inverter;
use crate::core::energy_supply::pv::{PhotovoltaicPanel, PhotovoltaicSystem};
use crate::core::heating_systems::boiler::{Boiler, BoilerServiceWaterCombi};
use crate::core::heating_systems::common::{HeatSourceWet, SpaceHeatSystem, SpaceHeatingService};
use crate::core::heating_systems::elec_storage_heater::{
    ElecStorageHeater, StorageHeaterDetailedResult,
};
use crate::core::heating_systems::emitters::{Emitters, EmittersDetailedResult};
use crate::core::heating_systems::heat_battery_pcm::{
    HeatBatteryPcm, HeatBatteryPcmServiceWaterDirect,
};
use crate::core::heating_systems::heat_network::{HeatNetwork, HeatNetworkServiceWaterDirect};
use crate::core::heating_systems::heat_pump::{HeatPump, HeatPumpHotWaterOnly};
use crate::core::heating_systems::instant_elec_heater::InstantElecHeater;
use crate::core::heating_systems::point_of_use::PointOfUse;
use crate::core::heating_systems::storage_tank::{
    HeatSourceWithStorageTank, HotWaterStorageTank, ImmersionHeater, PVDiverter,
    PositionedHeatSource, SmartHotWaterTank, SolarThermalSystem, StorageTank,
    StorageTankDetailedResult,
};
use crate::core::heating_systems::wwhrs::Wwhrs;
use crate::core::material_properties::WATER;
use crate::core::schedule::{
    expand_boolean_schedule, expand_events, expand_numeric_schedule, reject_nones, reject_nulls,
    ScheduleEvent, TypedScheduleEvent, WaterScheduleEventType,
};
use crate::core::space_heat_demand::building_element::{
    convert_uvalue_to_resistance, BuildingElement, BuildingElementAdjacentConditionedSpace,
    BuildingElementAdjacentUnconditionedSpaceSimple, BuildingElementGround, H_CE, H_RE,
    PITCH_LIMIT_HORIZ_CEILING, PITCH_LIMIT_HORIZ_FLOOR, R_SI_DOWNWARDS, R_SI_HORIZONTAL,
    R_SI_UPWARDS,
};
use crate::core::space_heat_demand::internal_gains::{
    ApplianceGains, EventApplianceGains, Gains, InternalGains,
};
use crate::core::space_heat_demand::thermal_bridge::{ThermalBridge, ThermalBridging};
use crate::core::space_heat_demand::ventilation::{
    InfiltrationVentilation, MechVentType, MechanicalVentilation, VentilationDetailedResult,
};
use crate::core::space_heat_demand::zone::{
    calc_vent_heat_transfer_coeff, AirChangesPerHourArgument, HeatBalance, HeatBalanceFieldName,
    Zone, ZoneTempInternalAir,
};
use crate::core::units::{kelvin_to_celsius, SECONDS_PER_HOUR, WATTS_PER_KILOWATT};
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::dhw_demand::{
    DemandVolTargetKey, DomesticHotWaterDemand, DomesticHotWaterDemandData, VolumeReference,
};
use crate::core::water_heat_demand::misc::water_demand_to_kwh;
use crate::external_conditions::{create_external_conditions, ExternalConditions};
use crate::input::{
    ApplianceGains as ApplianceGainsInput, ApplianceGainsDetails,
    BuildingElement as BuildingElementInput, ChargeLevel, ColdWaterSourceDetails,
    ColdWaterSourceInput, Control as ControlInput, ControlCombinations, ControlDetails, DuctType,
    EnergyDiverter, EnergySupplyDetails, EnergySupplyInput, FlowData, FuelType,
    HeatBattery as HeatBatteryInput, HeatPumpSourceType, HeatSource as HeatSourceInput,
    HeatSourceControlType, HeatSourceWetDetails, HotWaterSourceDetails,
    InfiltrationVentilation as InfiltrationVentilationInput, Input, InputForCalcHtcHlp,
    InternalGains as InternalGainsInput, InternalGainsDetails, OnSiteGeneration,
    PhotovoltaicInputs, PhotovoltaicSystem as PhotovoltaicSystemInput,
    SpaceCoolSystem as SpaceCoolSystemInput, SpaceCoolSystemDetails,
    SpaceHeatSystem as SpaceHeatSystemInput, SpaceHeatSystemDetails, SystemReference,
    ThermalBridging as ThermalBridgingInput, ThermalBridgingDetails, UValueInput, VentilationLeaks,
    WasteWaterHeatRecovery, WasteWaterHeatRecoveryDetails, WaterHeatingEvent, WaterHeatingEvents,
    WaterPipework, ZoneDictionary, ZoneInput, ZoneTemperatureControlBasis, MAIN_REFERENCE,
};
use crate::simulation_time::{SimulationTimeIteration, SimulationTimeIterator};
use crate::StringOrNumber;
use anyhow::{anyhow, bail};
use atomic_float::AtomicF64;
use indexmap::IndexMap;
#[cfg(feature = "indicatif")]
use indicatif::ProgressIterator;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use parking_lot::{Mutex, RwLock};
use serde_enum_str::{Deserialize_enum_str, Serialize_enum_str};
use smartstring::alias::String;
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::hash::Hash;
use std::io::{BufReader, Cursor, Read};
use std::iter::Sum;
use std::ops::{Add, AddAssign, Div};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

/// As of adopting Rust 1.82 as an MSRV we'll be able to declare this using constants as it supports floating-point arithmetic at compile time
fn temp_setpnt_heat_none() -> f64 {
    kelvin_to_celsius(0.0).expect("Not below absolute zero")
}

fn temp_setpnt_cool_none() -> f64 {
    kelvin_to_celsius(1.4e32).expect("Not below absolute zero")
}

// used for calculations re internal gains from pipework
const FRAC_DHW_ENERGY_INTERNAL_GAINS: f64 = 0.25;

fn control_from_input(
    control_input: &ControlInput,
    external_conditions: Arc<ExternalConditions>,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<Controls> {
    let mut core: Vec<HeatSourceControl> = Default::default();
    let mut extra: HashMap<String, Arc<Control>> = Default::default();

    // this is very ugly(!) but is just a reflection of the lack of clarity in the schema
    // and the way the variants-struct crate works;
    // we should be able to improve it in time
    if let Some(control) = &control_input.hot_water_timer.as_ref() {
        if let Some(ctrl) = single_control_from_details(
            control,
            external_conditions.clone(),
            simulation_time_iterator,
            control_input,
        )? {
            core.push(HeatSourceControl::HotWaterTimer(Arc::new(ctrl)));
        }
    }
    if let Some(control) = &control_input.window_opening.as_ref() {
        if let Some(ctrl) = single_control_from_details(
            control,
            external_conditions.clone(),
            simulation_time_iterator,
            control_input,
        )? {
            core.push(HeatSourceControl::WindowOpening(Arc::new(ctrl)));
        }
    }
    for (name, control) in &control_input.extra {
        if let Some(ctrl) = single_control_from_details(
            control,
            external_conditions.clone(),
            simulation_time_iterator,
            control_input,
        )? {
            extra.insert(name.to_string().into(), Arc::new(ctrl));
        }
    }

    Ok(Controls::new(core, extra))
}

fn single_control_from_details(
    details: &ControlDetails,
    external_conditions: Arc<ExternalConditions>,
    simulation_time_iterator: &SimulationTimeIterator,
    control_input: &ControlInput,
) -> anyhow::Result<Option<Control>> {
    Ok(match details {
        ControlDetails::OnOffTimer {
            start_day,
            time_series_step,
            schedule,
            allow_null,
            ..
        } => {
            let nullable = allow_null.unwrap_or(false);
            let schedule = if nullable {
                expand_boolean_schedule(schedule)
            } else {
                reject_nones(expand_boolean_schedule(schedule))?
            };
            Control::OnOffTime(OnOffTimeControl::new(
                schedule,
                *start_day,
                *time_series_step,
            ))
            .into()
        }
        ControlDetails::SetpointTimer {
            start_day,
            time_series_step,
            advanced_start,
            setpoint_bounds,
            schedule,
            ..
        } => Control::SetpointTime(SetpointTimeControl::new(
            expand_numeric_schedule(schedule),
            *start_day,
            *time_series_step,
            *setpoint_bounds,
            advanced_start.unwrap_or(0.),
            simulation_time_iterator.step_in_hours(),
        ))
        .into(),
        ControlDetails::ChargeTarget {
            charge_level,
            external_sensor,
            schedule,
            start_day,
            time_series_step,
            temp_charge_cut,
            temp_charge_cut_delta,
            logic_type,
            ..
        } => {
            let schedule = reject_nulls(expand_boolean_schedule(schedule))?;
            let temp_charge_cut_delta = temp_charge_cut_delta
                .as_ref()
                .map(|schedule| anyhow::Ok(reject_nulls(expand_numeric_schedule(schedule))?))
                .transpose()?;
            let logic_type = logic_type.unwrap_or_default();

            // Simulation manual charge control
            // Set charge level to 1.0 (max) for each day of simulation (plus 1)
            let vec_size = ((simulation_time_iterator.total_steps() as f64
                * simulation_time_iterator.step_in_hours()
                / 24.0)
                + 1.0)
                .ceil() as usize;
            // if charge_level is present in input, overwrite initial vector
            // user can specify a vector with all days (plus 1), or as a single float value to be used for each day
            let charge_level_vec = if let Some(charge) = charge_level {
                let required_length = vec_size;
                let mut charge_level_vec = match charge {
                    ChargeLevel::List(charge_vec) => charge_vec.iter().map(|x| Some(*x)).collect(),
                    ChargeLevel::Single(charge) => {
                        vec![Some(*charge); vec_size]
                    }
                    ChargeLevel::Schedule(schedule) => expand_numeric_schedule(schedule),
                };
                // Ensure charge_level has the required length by appending the last value as many times as necessary
                if charge_level_vec.len() < required_length {
                    let extension = vec![
                        charge_level_vec.iter().copied().last().ok_or_else(
                            || anyhow!("Provided charge level data was empty.")
                        )?;
                        required_length - charge_level_vec.len()
                    ];
                    charge_level_vec.extend(extension);
                }

                charge_level_vec
            } else {
                vec![Some(1.0); vec_size]
            };

            Control::Charge(ChargeControl::new(
                logic_type,
                schedule,
                simulation_time_iterator.step_in_hours(),
                *start_day,
                *time_series_step,
                charge_level_vec,
                *temp_charge_cut,
                temp_charge_cut_delta,
                None,
                None,
                external_conditions.clone(),
                external_sensor.clone(),
            )?)
            .into()
        }
        ControlDetails::OnOffCostMinimising {
            start_day,
            time_series_step,
            time_on_daily,
            schedule,
            ..
        } => Control::OnOffMinimisingTime(OnOffMinimisingTimeControl::new(
            reject_nulls(expand_numeric_schedule(schedule))?,
            *start_day,
            *time_series_step,
            *time_on_daily,
        ))
        .into(),
        ControlDetails::CombinationTime { combination } => {
            // resolved controls needs to be: IndexMap<String, Arc<Control>>

            /// Recursively collects all unique controls from the combination control dictionary.
            ///
            /// Args:
            ///    * `control_combination` : The combination control dictionary.
            ///    * `current_key` (str): The current key to process in the dictionary.
            ///    * `visited` (set): A set to track visited keys to detect circular references.
            fn collect_controls(
                control_combinations: &ControlCombinations,
                current_key: &str,
                visited: Option<&mut HashSet<String>>,
                control_input: &ControlInput,
                external_conditions: Arc<ExternalConditions>,
                simulation_time_iterator: &SimulationTimeIterator,
            ) -> anyhow::Result<IndexMap<String, Arc<Control>>> {
                let mut empty_set: HashSet<String> = Default::default();
                let visited = visited.unwrap_or(&mut empty_set);
                let mut controls: HashSet<String> = Default::default();

                if visited.contains(current_key) {
                    bail!("Error: Circular reference detected involving '{current_key}' in CombinationTimeControl. Exiting program.")
                }
                visited.insert(current_key.into());

                if control_combinations.contains_key(current_key) {
                    let control = &control_combinations[current_key];
                    for ctrl in control.controls.iter() {
                        if control_combinations.contains_key(ctrl) {
                            // If the control is another combination
                            controls.extend(
                                collect_controls(
                                    control_combinations,
                                    ctrl,
                                    Some(visited),
                                    control_input,
                                    external_conditions.clone(),
                                    simulation_time_iterator,
                                )?
                                .keys()
                                .cloned(),
                            );
                        } else {
                            controls.insert(ctrl.clone());
                        }
                    }
                }

                let control_instances = controls
                    .iter()
                    .filter_map(|control_name| {
                        let control_details = match control_input.get(control_name) {
                            Some(control_details) => control_details,
                            None => return Some(Err(anyhow!("There was a reference to a control with name '{control_name}' that was not provided.")))
                        };
                        let control = match single_control_from_details(
                            control_details,
                            external_conditions.clone(),
                            simulation_time_iterator,
                            control_input,
                        ) {
                            Ok(c) => c,
                            Err(_) => return Some(Err(anyhow!("The control name '{control_name}' refers to a control that cannot be included as part of a combination.")))
                        };
                        control.map(|control| {
                            Ok((
                                control_name.clone(),
                                Arc::new(control),
                            ))
                        })
                    })
                    .try_collect()?;

                Ok(control_instances)
            }

            let resolved_controls = collect_controls(
                combination,
                MAIN_REFERENCE,
                None,
                control_input,
                external_conditions,
                simulation_time_iterator,
            )?;

            Control::CombinationTime(CombinationTimeControl::new(
                combination.clone(),
                resolved_controls,
            )?)
            .into()
        }
    })
}

fn init_resistance_or_uvalue_from_data(
    thermal_resistance_construction: Option<f64>,
    u_value: Option<f64>,
    pitch: f64,
) -> anyhow::Result<f64> {
    Ok(
        if let Some(thermal_resistance_construction) = thermal_resistance_construction {
            thermal_resistance_construction
        } else {
            convert_uvalue_to_resistance(
            u_value.ok_or_else(|| {
                anyhow!(
                    "Neither thermal_resistance_construction nor u_value were provided for one of the building element inputs."
                )
            })?,
            pitch,
        )
        },
    )
}

fn init_resistance_or_uvalue_from_input_struct(
    data: &UValueInput,
    pitch: f64,
) -> anyhow::Result<f64> {
    let (thermal_resistance_construction, u_value) = match data {
        UValueInput::UValue { u_value, .. } => (None, Some(*u_value)),
        UValueInput::ThermalResistanceConstruction {
            thermal_resistance_construction,
        } => (Some(*thermal_resistance_construction), None),
    };
    init_resistance_or_uvalue_from_data(thermal_resistance_construction, u_value, pitch)
}

/// Return thermal resistance of construction (thermal_resistance_construction) based on alternative inputs
///
/// User will either provide thermal_resistance_construction directly or provide u_value which needs to be converted
fn init_resistance_or_uvalue(element: &BuildingElementInput) -> anyhow::Result<f64> {
    fn params_from_u_value_input_and_pitch(
        u_value_input: &UValueInput,
        pitch: f64,
    ) -> (Option<f64>, Option<f64>, f64) {
        let (thermal_resistance_construction, u_value) = match u_value_input {
            UValueInput::UValue { u_value, .. } => (None, Some(*u_value)),
            UValueInput::ThermalResistanceConstruction {
                thermal_resistance_construction,
            } => (Some(*thermal_resistance_construction), None),
        };
        (thermal_resistance_construction, u_value, pitch)
    }

    let (thermal_resistance_construction, u_value, pitch) = match element {
        BuildingElementInput::Opaque {
            u_value_input,
            pitch,
            ..
        } => params_from_u_value_input_and_pitch(u_value_input, *pitch),
        BuildingElementInput::Transparent {
            u_value_input,
            pitch,
            ..
        } => params_from_u_value_input_and_pitch(u_value_input, *pitch),
        BuildingElementInput::Ground { pitch, u_value, .. } => (None, Some(*u_value), *pitch),
        BuildingElementInput::AdjacentConditionedSpace {
            pitch,
            u_value_input,
            ..
        } => params_from_u_value_input_and_pitch(u_value_input, *pitch),
        BuildingElementInput::AdjacentUnconditionedSpace {
            pitch,
            u_value_input,
            ..
        } => params_from_u_value_input_and_pitch(u_value_input, *pitch),
    };
    init_resistance_or_uvalue_from_data(thermal_resistance_construction, u_value, pitch)
}

/// Calculate heat transfer coefficient (HTC) and heat loss parameter (HLP)
/// according to the SAP10.2 specification
pub fn calc_htc_hlp<T: InputForCalcHtcHlp>(input: &T) -> anyhow::Result<HtcHlpCalculation> {
    let simtime = input.simulation_time();
    let external_conditions = Arc::from(create_external_conditions(
        (*input.external_conditions()).clone(),
        &simtime.iter(),
    )?);
    let energy_supply_unmet_demand =
        EnergySupplyBuilder::new(FuelType::UnmetDemand, simtime.total_steps()).build();
    let mut energy_supplies: IndexMap<String, Arc<RwLock<EnergySupply>>> = [(
        "_unmet_demand".into(),
        Arc::new(RwLock::new(energy_supply_unmet_demand)),
    )]
    .into();
    for (name, data) in input.energy_supply().iter() {
        energy_supplies.insert(
            name.into(),
            Arc::new(RwLock::new(
                EnergySupplyBuilder::new(data.fuel, simtime.total_steps()).build(),
            )),
        );
    }

    let controls = control_from_input(
        input.control(),
        external_conditions.clone(),
        &simtime.iter(),
    )?;

    let ventilation = InfiltrationVentilation::create(
        input.infiltration_ventilation(),
        input.zone(),
        false,
        &energy_supplies,
        &controls,
    )?;

    fn calc_heat_loss(data: &BuildingElementInput) -> anyhow::Result<f64> {
        // Calculate thermal_resistance_construction from u_value if only the latter has been provided
        let thermal_resistance_construction = init_resistance_or_uvalue(data)?;
        let r_si = match data.pitch() {
            PITCH_LIMIT_HORIZ_CEILING..=PITCH_LIMIT_HORIZ_FLOOR => R_SI_HORIZONTAL,
            ..PITCH_LIMIT_HORIZ_CEILING => R_SI_UPWARDS,
            PITCH_LIMIT_HORIZ_FLOOR.. => R_SI_DOWNWARDS,
            _ => unreachable!("Rust cannot tell that above is exhaustive"),
        };
        let r_se = 1. / (H_CE + H_RE);

        Ok(match data {
            BuildingElementInput::Opaque { area_input, .. } => {
                let u_value = 1.0 / (thermal_resistance_construction + r_se + r_si);
                area_input.area() * u_value
            }
            BuildingElementInput::Transparent { area_input, .. } => {
                let r_curtains_blinds = 0.04;
                let u_value =
                    1.0 / (thermal_resistance_construction + r_se + r_si + r_curtains_blinds);
                let area = area_input.area();
                area * u_value
            }
            BuildingElementInput::Ground { area, u_value, .. } => *area * *u_value,
            BuildingElementInput::AdjacentConditionedSpace { .. } => 0.,
            BuildingElementInput::AdjacentUnconditionedSpace { area, .. } => {
                let u_value = 1.0 / (thermal_resistance_construction + r_se + r_si);
                area.ok_or_else(|| {
                    anyhow!("AdjacentConditionedSpace building element is expected have an area")
                })? * u_value
            }
        })
    }

    fn calc_heat_transfer_coeff(data: &ThermalBridgingInput) -> f64 {
        // If data is for individual thermal bridges, initialise the relevant
        // objects and return a list of them. Otherwise, just use the overall
        // figure given.
        match data {
            ThermalBridgingInput::Elements(bridges) => bridges
                .values()
                .map(|bridge| match bridge {
                    ThermalBridgingDetails::Linear {
                        linear_thermal_transmittance,
                        length,
                        ..
                    } => *linear_thermal_transmittance * *length,
                    ThermalBridgingDetails::Point {
                        heat_transfer_coefficient,
                    } => *heat_transfer_coefficient,
                })
                .sum::<f64>(),
            ThermalBridgingInput::Number(num) => *num,
        }
    }

    let calc_htc = |zone: &ZoneInput| -> anyhow::Result<(f64, f64, f64)> {
        let wind_speed = external_conditions.wind_speed_annual().ok_or_else(|| {
            anyhow!("Expected external conditions to contain data for entire year")
        })?;
        let wind_direction = external_conditions.wind_direction_annual();
        let temp_int_air = input.temp_internal_air_static_calcs();
        let temp_ext_air = external_conditions.air_temp_annual_daily_average_min();
        let ach_min = input.infiltration_ventilation().ach_min_static_calcs;
        let ach_max = input.infiltration_ventilation().ach_max_static_calcs;
        let initial_r_v_arg = input
            .infiltration_ventilation()
            .vent_opening_ratio_init
            .unwrap_or(1.);

        // Adjust vent position if required to attempt to meet min or max ach
        let r_v_arg = ventilation.find_r_v_arg_within_bounds(
            ach_min,
            ach_max,
            initial_r_v_arg,
            wind_speed,
            wind_direction,
            temp_int_air,
            temp_ext_air,
            Some(0.),
            0.,
            None,
            simtime.iter().current_iteration(),
        )?;
        let air_changes_per_hour = ventilation.calc_air_changes_per_hour(
            wind_speed,
            wind_direction,
            temp_int_air,
            temp_ext_air,
            r_v_arg,
            Some(0.),
            0.,
            Some(ReportingFlag::Min),
            simtime.iter().current_iteration(),
        )?;
        let total_vent_heat_loss = calc_vent_heat_transfer_coeff(zone.volume, air_changes_per_hour);

        // Calculate fabric heat loss and total floor area
        let total_fabric_heat_loss = zone
            .building_elements
            .values()
            .map(calc_heat_loss)
            .try_collect::<f64, Vec<f64>, anyhow::Error>()?
            .iter()
            .sum::<f64>();

        // Read in thermal bridging data
        let tb_heat_trans_coeff = calc_heat_transfer_coeff(&zone.thermal_bridging);

        Ok((
            total_fabric_heat_loss,
            tb_heat_trans_coeff,
            total_vent_heat_loss,
        ))
    };

    // Calculate the total fabric heat loss, total heat capacity, total ventilation heat
    // loss and total heat transfer coeffient for thermal bridges across all zones

    let mut htc_map: IndexMap<String, f64> = Default::default();
    let mut hlp_map: IndexMap<String, f64> = Default::default();
    let mut zone_area: IndexMap<String, f64> = Default::default();

    for (z_name, zone) in input.zone().iter() {
        let (fabric_heat_loss, thermal_bridges, vent_heat_loss) = calc_htc(zone)?;
        // Calculate the heat transfer coefficent (HTC), in W / K
        // TODO (from Python) check ventilation losses are correct
        let htc = fabric_heat_loss + thermal_bridges + vent_heat_loss;
        let hlp = htc / zone.area;
        htc_map.insert(z_name.into(), htc);
        hlp_map.insert(z_name.into(), hlp);
        zone_area.insert(z_name.into(), zone.area);
    }

    let total_htc = htc_map.values().sum::<f64>();
    let total_floor_area = zone_area.values().sum::<f64>();
    let total_hlp = total_htc / total_floor_area;

    Ok(HtcHlpCalculation {
        total_htc,
        total_hlp,
        _htc_map: htc_map,
        _hlp_map: hlp_map,
    })
}

pub struct HtcHlpCalculation {
    pub(crate) total_htc: f64,
    pub(crate) total_hlp: f64,
    pub _htc_map: IndexMap<String, f64>,
    pub(crate) _hlp_map: IndexMap<String, f64>,
}

#[derive(Debug)]
pub struct Corpus {
    pub(crate) simulation_time: Arc<SimulationTimeIterator>,
    pub(crate) external_conditions: Arc<ExternalConditions>,
    pre_heated_water_sources: IndexMap<String, HotWaterStorageTank>,
    pub(crate) energy_supplies: IndexMap<String, Arc<RwLock<EnergySupply>>>,
    pub(crate) internal_gains: InternalGainsCollection,
    pub(crate) domestic_hot_water_demand: DomesticHotWaterDemand,
    r_v_arg: AtomicF64,
    pub(crate) ventilation: Arc<InfiltrationVentilation>,
    pub(crate) zones: IndexMap<String, Arc<Zone>>,
    pub(crate) energy_supply_conn_unmet_demand_zone: IndexMap<String, Arc<EnergySupplyConnection>>,
    pub(crate) heat_system_name_for_zone: IndexMap<String, Vec<String>>,
    pub(crate) cool_system_name_for_zone: IndexMap<String, Vec<String>>,
    pub total_floor_area: f64,
    pub(crate) total_volume: f64,
    pub(crate) wet_heat_sources: IndexMap<String, WetHeatSource>,
    pub(crate) hot_water_sources: IndexMap<String, HotWaterSource>,
    pub(crate) heat_sources_wet_with_buffer_tank: Vec<String>,
    pub(crate) space_heat_systems: IndexMap<String, Arc<Mutex<SpaceHeatSystem>>>,
    pub(crate) space_cool_systems: IndexMap<String, AirConditioning>,
    pub(crate) on_site_generation: IndexMap<String, PhotovoltaicSystem>,
    pub(crate) diverters: Vec<Arc<RwLock<PVDiverter>>>,
    required_vent_data: Option<RequiredVentData>,
    energy_supply_conn_names_for_hot_water_source: IndexMap<String, Vec<String>>,
    energy_supply_conn_names_for_heat_systems: IndexMap<String, String>,
    timestep_end_calcs: Arc<RwLock<Vec<WetHeatSource>>>,
    initial_loop: AtomicBool,
    detailed_output_heating_cooling: bool,
    vent_adjust_min_control: Option<Arc<Control>>,
    vent_adjust_max_control: Option<Arc<Control>>,
    temp_internal_air_prev: Arc<RwLock<f64>>,
    smart_appliance_controls: IndexMap<String, Arc<SmartApplianceControl>>,
}

impl Corpus {
    pub fn from_inputs(
        input: &Input,
        external_conditions: Option<&ExternalConditions>,
        tariff_file_path: Option<&str>,
        output_options: &OutputOptions,
    ) -> anyhow::Result<Self> {
        let simulation_time_iterator = Arc::new(input.simulation_time.iter());

        let external_conditions = Arc::new(match external_conditions {
            Some(external_conditions) => external_conditions.clone(),
            None => create_external_conditions(
                input.external_conditions.as_ref().to_owned(),
                &simulation_time_iterator,
            )?,
        });

        let diverter_types: DiverterTypes = input
            .energy_supply
            .iter()
            .flat_map(|(name, supply)| {
                diverter_from_energy_supply(supply).map(|diverter| (name.into(), diverter))
            })
            .collect();
        let mut diverters: Vec<Arc<RwLock<PVDiverter>>> = Default::default();

        let cold_water_sources = cold_water_sources_from_input(&input.cold_water_source);
        let wwhrs = wwhrs_from_input(
            input.waste_water_heat_recovery.as_ref(),
            &cold_water_sources,
            simulation_time_iterator.current_iteration(),
        )?;

        let mut energy_supplies = energy_supplies_from_input(
            &input.energy_supply,
            simulation_time_iterator.clone().as_ref(),
            tariff_file_path,
            external_conditions.clone(),
        )?;

        let controls = control_from_input(
            &input.control,
            external_conditions.clone(),
            simulation_time_iterator.clone().as_ref(),
        )?;

        let event_schedules = event_schedules_from_input(
            &input.water_heating_events,
            simulation_time_iterator.as_ref(),
        )?;

        let domestic_hot_water_demand = DomesticHotWaterDemand::new(
            input.hot_water_demand.shower.clone(),
            input.hot_water_demand.bath.clone(),
            input.hot_water_demand.other_water_use.clone(),
            // match &input.hot_water_source.hot_water_cylinder {
            //     HotWaterSourceDetails::PointOfUse { .. } => None,
            //     _ => Some(input.hot_water_demand.water_distribution.clone()),
            // },
            Default::default(), // TODO: migrate properly for 1.0.0a1
            &cold_water_sources,
            &wwhrs,
            &energy_supplies,
            event_schedules,
        )?;

        let total_volume = input.zone.values().map(|zone| zone.volume).sum::<f64>();

        let mut heat_system_name_for_zone: IndexMap<String, Vec<String>> = Default::default();
        let mut cool_system_name_for_zone: IndexMap<String, Vec<String>> = Default::default();

        // infiltration ventilation
        let (
            infiltration_ventilation,
            window_adjust_control,
            vent_adjust_min_control,
            vent_adjust_max_control,
        ) = infiltration_ventilation_from_input(
            &input.zone,
            &input.infiltration_ventilation,
            &controls,
            &mut energy_supplies,
            output_options.detailed_output_heating_cooling,
        )?;

        let infiltration_ventilation = Arc::from(infiltration_ventilation);

        let required_vent_data = required_vent_data_from_input(&input.control)?;

        let zones: IndexMap<String, Arc<Zone>> = input
            .zone
            .iter()
            .map(|(zone_name, zone)| -> anyhow::Result<(String, Arc<Zone>)> {
                Ok((zone_name.into(), {
                    let zone_for_corpus = zone_from_input(
                        zone,
                        zone_name,
                        &mut heat_system_name_for_zone,
                        &mut cool_system_name_for_zone,
                        external_conditions.clone(),
                        infiltration_ventilation.clone(),
                        window_adjust_control.clone(),
                        &controls,
                        output_options.print_heat_balance,
                        simulation_time_iterator.clone().as_ref(),
                    )?;

                    Arc::new(zone_for_corpus)
                }))
            })
            .collect::<anyhow::Result<_>>()?;

        let energy_supply_conn_unmet_demand_zone = set_up_energy_supply_unmet_demand_zones(
            energy_supplies[UNMET_DEMAND_SUPPLY_NAME].clone(),
            &input.zone,
        );

        let total_floor_area = zones.values().fold(0., |acc, zone| zone.area() + acc);

        // Internal gains is an ordered IndexMap. This is because load shifting behaviours
        // of appliance gains depend on other energy demand in the dwelling at any given time,
        // so depend on the order in which gains are considered by the engine.
        // See check_priority() in apply_appliance_gains_from_input()
        let mut internal_gains =
            internal_gains_from_input(&input.internal_gains, total_floor_area)?;

        // setup smart control for loadshifting
        let mut smart_appliance_controls: IndexMap<String, Arc<SmartApplianceControl>> =
            Default::default();
        for (smart_appliance_name, smart_appliance_data) in &input.smart_appliance_controls {
            // TODO (from Python) - power_timeseries is a redundant input,
            // we could obtain power_timeseries here from the list smartappctrldata['Appliances'],
            // by looking up each string in the list in proj_dict['ApplianceGains'], and
            // summing together the demand schedules of listed items
            // this will require the use of EventApplianceGains.__event_to_schedule()
            // for event based appliance use
            smart_appliance_controls.insert(
                smart_appliance_name.into(),
                SmartApplianceControl::new(
                    &smart_appliance_data.power_timeseries,
                    smart_appliance_data.time_series_step,
                    &simulation_time_iterator,
                    smart_appliance_data.non_appliance_demand_24hr.clone(),
                    smart_appliance_data.battery_24hr.clone(),
                    &energy_supplies,
                    smart_appliance_data.appliances.iter().map_into().collect(),
                )?
                .into(),
            );
        }

        //  Add internal gains from applicances to the internal gains dictionary and
        //  create an energy supply connection for appliances
        // work out order in which to process loadshifting appliances
        apply_appliance_gains_from_input(
            &mut internal_gains,
            &input.appliance_gains,
            &mut energy_supplies,
            total_floor_area,
            &smart_appliance_controls,
            simulation_time_iterator.as_ref(),
        )?;

        let timestep_end_calcs: Arc<RwLock<Vec<WetHeatSource>>> = Default::default();
        let mut heat_sources_wet_with_buffer_tank: Vec<String> = vec![];
        let mechanical_ventilations = infiltration_ventilation.mech_vents();

        let temp_internal_air_prev: Arc<RwLock<f64>> = Default::default();

        let mut wet_heat_sources: IndexMap<String, WetHeatSource> = input
            .heat_source_wet
            .clone()
            .unwrap_or_default()
            .iter()
            .map(|(name, heat_source_wet_details)| {
                let heat_source = heat_source_wet_from_input(
                    name,
                    (*heat_source_wet_details).clone(),
                    external_conditions.clone(),
                    simulation_time_iterator.clone(),
                    mechanical_ventilations,
                    zones.len(),
                    shareable_fn(&temp_internal_air_prev),
                    &controls,
                    &mut energy_supplies,
                    output_options.detailed_output_heating_cooling,
                )?;
                timestep_end_calcs.write().push(heat_source.clone());
                if let HeatSourceWetDetails::HeatPump {
                    buffer_tank: Some(_),
                    ..
                } = heat_source_wet_details
                {
                    heat_sources_wet_with_buffer_tank.push(name.into());
                }
                anyhow::Ok((name.into(), heat_source))
            })
            .collect::<anyhow::Result<IndexMap<String, WetHeatSource>>>()?;

        let mut energy_supply_conn_names_for_hot_water_source: IndexMap<String, Vec<String>> =
            Default::default();

        // processing pre-heated sources
        let mut pre_heated_water_sources: IndexMap<String, HotWaterStorageTank> =
            Default::default();

        for (source_name, source_details) in &input.pre_heated_water_source {
            let (heat_source, energy_conn_names) = hot_water_source_from_input(
                source_name.into(),
                source_details,
                &cold_water_sources,
                &pre_heated_water_sources,
                &mut wet_heat_sources,
                &wwhrs,
                &controls,
                &mut energy_supplies,
                &diverter_types,
                &mut diverters,
                shareable_fn(&temp_internal_air_prev),
                simulation_time_iterator.clone().as_ref(),
                external_conditions.clone(),
                output_options.detailed_output_heating_cooling,
            )?;
            energy_supply_conn_names_for_hot_water_source
                .insert(source_name.into(), energy_conn_names);
            if let HotWaterSource::PreHeated(source) = heat_source {
                pre_heated_water_sources.insert(source_name.into(), source);
            } else {
                bail!("Pre-heated water sources must be storage tanks");
            }
        }

        let mut hot_water_sources: IndexMap<String, HotWaterSource> = Default::default();
        for (name, data) in input.hot_water_source.iter() {
            let (hot_water_source, hw_cylinder_conn_names) = hot_water_source_from_input(
                name.into(),
                data,
                &cold_water_sources,
                &pre_heated_water_sources,
                &mut wet_heat_sources,
                &wwhrs,
                &controls,
                &mut energy_supplies,
                &diverter_types,
                &mut diverters,
                shareable_fn(&temp_internal_air_prev),
                simulation_time_iterator.clone().as_ref(),
                external_conditions.clone(),
                output_options.detailed_output_heating_cooling,
            )?;
            hot_water_sources.insert(name.into(), hot_water_source);
            energy_supply_conn_names_for_hot_water_source
                .insert(name.into(), hw_cylinder_conn_names);
        }

        let mut heat_system_names_requiring_overvent: Vec<String> = Default::default();

        let (space_heat_systems, energy_supply_conn_names_for_heat_systems) = input
            .space_heat_system
            .as_ref()
            .map(|system| {
                anyhow::Ok(space_heat_systems_from_input(
                    system,
                    &controls,
                    &mut energy_supplies,
                    simulation_time_iterator.as_ref(),
                    &wet_heat_sources,
                    &mut heat_system_names_requiring_overvent,
                    &heat_system_name_for_zone,
                    &zones,
                    &heat_sources_wet_with_buffer_tank
                        .iter()
                        .cloned()
                        .collect_vec(),
                    external_conditions.clone(),
                    output_options.detailed_output_heating_cooling,
                )?)
            })
            .transpose()?
            .unwrap_or_default();

        let space_cool_systems = input
            .space_cool_system
            .as_ref()
            .map(|system| {
                anyhow::Ok(space_cool_systems_from_input(
                    system,
                    cool_system_name_for_zone
                        .values()
                        .flatten()
                        .map(|s| s.as_str())
                        .unique()
                        .collect::<Vec<_>>(),
                    &controls,
                    &mut energy_supplies,
                    &simulation_time_iterator,
                )?)
            })
            .transpose()?
            .unwrap_or_default();

        let on_site_generation = input
            .on_site_generation
            .as_ref()
            .map(|on_site_generation| {
                anyhow::Ok(on_site_generation_from_input(
                    on_site_generation,
                    &mut energy_supplies,
                    external_conditions.clone(),
                    &simulation_time_iterator,
                )?)
            })
            .transpose()?
            .unwrap_or_default();

        Ok(Self {
            simulation_time: simulation_time_iterator,
            external_conditions,
            pre_heated_water_sources,
            energy_supplies,
            internal_gains,
            domestic_hot_water_demand,
            r_v_arg: AtomicF64::new(
                input
                    .infiltration_ventilation
                    .vent_opening_ratio_init
                    .unwrap_or(1.), // default to 1 if unspecified
            ),
            ventilation: infiltration_ventilation,
            zones,
            energy_supply_conn_unmet_demand_zone,
            heat_system_name_for_zone,
            cool_system_name_for_zone,
            total_floor_area,
            total_volume,
            wet_heat_sources,
            hot_water_sources,
            heat_sources_wet_with_buffer_tank,
            space_heat_systems,
            space_cool_systems,
            on_site_generation,
            diverters,
            required_vent_data,
            energy_supply_conn_names_for_hot_water_source,
            energy_supply_conn_names_for_heat_systems,
            timestep_end_calcs,
            initial_loop: AtomicBool::new(false),
            detailed_output_heating_cooling: output_options.detailed_output_heating_cooling,
            vent_adjust_min_control,
            vent_adjust_max_control,
            temp_internal_air_prev,
            smart_appliance_controls,
        })
    }

    pub fn total_floor_area(&self) -> f64 {
        self.total_floor_area
    }

    /// Calculate the total heat capacity normalised for floor area
    pub fn calc_hcp(&self) -> f64 {
        // TODO (from Python) party walls and solid doors should be exluded according to SAP spec - if party walls are
        // assumed to be ZTU building elements this could be set to zero?
        let total_heat_capacity = self
            .zones
            .values()
            .map(|zone| zone.total_heat_capacity())
            .sum::<f64>();

        total_heat_capacity / self.total_floor_area
    }

    /// Calculate the heat loss form factor, defined as exposed area / floor area
    pub fn calc_hlff(&self) -> f64 {
        let total_heat_loss_area = self
            .zones
            .values()
            .map(|zone| zone.total_heat_loss_area())
            .sum::<f64>();

        total_heat_loss_area / self.total_floor_area
    }

    pub fn update_temp_internal_air(&self) {
        *self.temp_internal_air_prev.write() =
            temp_internal_air_for_zones(&self.zones, self.total_volume);
    }

    /// Return the volume-weighted average internal air temperature from the previous timestep
    /// Some parts of the calculation rely on the whole-dwelling internal air
    /// temperature before it has been calculated for the current timestep, so
    /// we use the air temperature calculated in the previous timestep as an
    /// approximation. Note that this returns a stored value rather than
    /// calculating from the internal air temperature of each Zone object,
    /// because this function may be called after the temperatures of some Zone
    /// objects have been updated for the current timestep but before the
    /// temperatures of other Zone objects have been updated, which would be
    /// inconsistent.
    fn temp_internal_air_prev_timestep(&self) -> f64 {
        *self.temp_internal_air_prev.read()
    }

    /// Return:
    ///  # - losses from internal distribution pipework (kWh)
    ///  # - losses from external distribution pipework (kWh)
    ///  # - internal gains due to hot water use (kWh)
    fn pipework_losses_and_internal_gains_from_hw(
        &self,
        delta_t_h: f64,
        vol_hot_water_at_tapping_point: f64,
        hw_duration: f64,
        no_of_hw_events: usize,
        temp_hot_water: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (f64, f64, f64) {
        let (pw_losses_internal, pw_losses_external) = self.calc_pipework_losses(
            delta_t_h,
            hw_duration,
            no_of_hw_events,
            temp_hot_water,
            simulation_time_iteration,
        );

        let gains_internal_dhw_use = FRAC_DHW_ENERGY_INTERNAL_GAINS
            * water_demand_to_kwh(
                vol_hot_water_at_tapping_point,
                temp_hot_water,
                self.temp_internal_air_prev_timestep(),
            );

        (
            pw_losses_internal,
            pw_losses_external,
            gains_internal_dhw_use,
        )
    }

    fn pipework_losses_and_internal_gains_from_hw_storage_tank(
        &self,
        delta_t_h: f64,
        volume_water_remove_from_tank: f64,
        hw_duration: f64,
        no_of_hw_events: usize,
        temp_final_drawoff: f64,
        temp_average_drawoff: f64,
        temp_hot_water: f64,
        vol_hot_water_equiv_elec_shower: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (f64, f64, f64) {
        let (pw_losses_internal, pw_losses_external) = self.calc_pipework_losses(
            delta_t_h,
            hw_duration,
            no_of_hw_events,
            temp_final_drawoff,
            simulation_time_iteration,
        );
        let gains_internal_dhw_use_storagetank = FRAC_DHW_ENERGY_INTERNAL_GAINS
            * water_demand_to_kwh(
                volume_water_remove_from_tank,
                temp_average_drawoff,
                self.temp_internal_air_prev_timestep(),
            );

        let gains_internal_dhw_use_ies = FRAC_DHW_ENERGY_INTERNAL_GAINS
            * water_demand_to_kwh(
                vol_hot_water_equiv_elec_shower,
                temp_hot_water,
                self.temp_internal_air_prev_timestep(),
            );

        let gains_internal_dhw_use =
            gains_internal_dhw_use_storagetank + gains_internal_dhw_use_ies;

        // Return:
        // losses from internal distribution pipework (kWh)
        // losses from external distribution pipework (kWh)
        // internal gains due to hot water use (kWh)
        (
            pw_losses_internal,
            pw_losses_external,
            gains_internal_dhw_use,
        )
    }

    fn calc_pipework_losses(
        &self,
        delta_t_h: f64,
        hw_duration: f64,
        no_of_hw_events: usize,
        temp_hot_water: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (f64, f64) {
        let demand_water_temperature = temp_hot_water;
        let internal_air_temperature = self.temp_internal_air_prev_timestep();
        let external_air_temperature = self
            .external_conditions
            .air_temp(&simulation_time_iteration);

        self.domestic_hot_water_demand.calc_pipework_losses(
            delta_t_h,
            hw_duration,
            no_of_hw_events,
            demand_water_temperature,
            internal_air_temperature,
            external_air_temperature,
        )
    }

    /// Calculate the losses in the buffer tank
    fn calc_internal_gains_buffer_tank(&self) -> f64 {
        self.heat_sources_wet_with_buffer_tank
            .iter()
            .map(
                |heat_source_name| match self.wet_heat_sources.get(heat_source_name) {
                    Some(heat_source) => match heat_source {
                        WetHeatSource::HeatPump(heat_pump) => heat_pump.lock().buffer_int_gains(),
                        _ => unreachable!(),
                    },
                    None => 0.,
                },
            )
            .sum::<f64>()
    }

    /// Calculate the losses/gains in the MVHR ductwork
    fn calc_internal_gains_ductwork(&self, simulation_time: SimulationTimeIteration) -> f64 {
        let mut internal_gains_ductwork_watts = 0.0;
        let space_heating_ductwork = self.ventilation.space_heating_ductworks();
        for mvhr_ductwork in space_heating_ductwork.values() {
            // assume MVHR unit is running 100% of the time
            for duct in mvhr_ductwork {
                match duct.duct_type() {
                    DuctType::Intake | DuctType::Exhaust => {
                        // Heat loss from intake or exhaust ducts is to zone, so add
                        // to internal gains (may be negative gains)
                        internal_gains_ductwork_watts += duct.total_duct_heat_loss(
                            self.temp_internal_air_prev_timestep(),
                            self.external_conditions.air_temp(&simulation_time),
                        );
                    }
                    DuctType::Supply | DuctType::Extract => {
                        // Heat loss from supply and extract ducts is to outside, so
                        // subtract from internal gains
                        internal_gains_ductwork_watts -= duct.total_duct_heat_loss(
                            self.temp_internal_air_prev_timestep(),
                            self.external_conditions.air_temp(&simulation_time),
                        );
                    }
                }
            }
        }
        internal_gains_ductwork_watts
    }

    fn space_heat_internal_gains_for_zone(
        &self,
        zone: &Zone,
        gains_internal_dhw: f64,
        internal_gains_ductwork_per_m3: f64,
        gains_internal_buffer_tank: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // Initialise to dhw internal gains split proportionally to zone floor area
        let mut gains_internal_zone =
            (gains_internal_buffer_tank + gains_internal_dhw) * zone.area() / self.total_floor_area;

        for internal_gains in self.internal_gains.values() {
            gains_internal_zone += internal_gains.total_internal_gain_in_w(zone.area(), simtime)?;
        }

        // Add gains from ventilation fans (also calculates elec demand from fans)
        // TODO (from Python) Remove the branch on the type of ventilation (find a better way)
        let mech_vents = self.ventilation.mech_vents();
        for mech_vent in mech_vents.iter() {
            gains_internal_zone +=
                mech_vent.fans(zone.volume(), self.total_volume, None, &simtime)?;
            gains_internal_zone += internal_gains_ductwork_per_m3 * zone.volume();
        }

        Ok(gains_internal_zone)
    }

    /// Look up relevant heating and cooling systems for the specified zone
    fn heat_cool_systems_for_zone(
        &self,
        z_name: &str,
        simtime: SimulationTimeIteration,
    ) -> (Vec<String>, Vec<String>, SetpointsAndConvectiveFractions) {
        let SetpointsAndConvectiveFractions {
            temp_setpnt_heat: temp_setpnt_heat_system,
            temp_setpnt_cool: temp_setpnt_cool_system,
            frac_convective_heat: frac_convective_heat_system,
            frac_convective_cool: frac_convective_cool_system,
        } = self.setpoints_and_convective_fractions(
            &self.heat_system_name_for_zone[z_name],
            &self.cool_system_name_for_zone[z_name],
            simtime,
        );

        // Sort heating and cooling systems by setpoint (highest first for
        // heating, lowest first for cooling)
        // In the event of two systems having the same setpoint, the one
        // listed first by the user takes priority
        let h_name_list_sorted: Vec<String> = temp_setpnt_heat_system
            .iter()
            .sorted_by(|a, b| OrderedFloat(*a.1).cmp(&OrderedFloat(*b.1)))
            .rev()
            .map(|x| x.0.to_owned())
            .collect();
        let c_name_list_sorted: Vec<String> = temp_setpnt_cool_system
            .iter()
            .sorted_by(|a, b| OrderedFloat(*a.1).cmp(&OrderedFloat(*b.1)))
            .map(|x| x.0.to_owned())
            .collect();

        (
            h_name_list_sorted,
            c_name_list_sorted,
            SetpointsAndConvectiveFractions {
                temp_setpnt_heat: temp_setpnt_heat_system,
                temp_setpnt_cool: temp_setpnt_cool_system,
                frac_convective_heat: frac_convective_heat_system,
                frac_convective_cool: frac_convective_cool_system,
            },
        )
    }

    fn setpoints_and_convective_fractions(
        &self,
        h_name_list: &Vec<String>,
        c_name_list: &Vec<String>,
        simtime: SimulationTimeIteration,
    ) -> SetpointsAndConvectiveFractions {
        let mut frac_convective_heat: IndexMap<String, f64> = Default::default();
        let mut frac_convective_cool: IndexMap<String, f64> = Default::default();
        let mut temp_setpnt_heat: IndexMap<String, f64> = Default::default();
        let mut temp_setpnt_cool: IndexMap<String, f64> = Default::default();

        for h_name in h_name_list {
            match h_name.as_str() {
                h_name @ "" => {
                    frac_convective_heat.insert((*h_name).into(), 1.0);
                    temp_setpnt_heat.insert(h_name.into(), temp_setpnt_heat_none());
                }
                h_name => {
                    let space_heat_system = self.space_heat_systems.get(h_name).unwrap().lock();
                    frac_convective_heat
                        .insert(h_name.into(), space_heat_system.frac_convective(simtime));
                    temp_setpnt_heat.insert(
                        (*h_name).into(),
                        space_heat_system
                            .temp_setpnt(simtime)
                            .unwrap_or_else(temp_setpnt_heat_none),
                    );
                }
            }
        }

        for c_name in c_name_list {
            match c_name.as_str() {
                c_name @ "" => {
                    frac_convective_cool.insert((*c_name).into(), 1.0);
                    temp_setpnt_cool.insert(c_name.into(), temp_setpnt_cool_none());
                }
                c_name => {
                    let space_cool_system = self.space_cool_systems.get(c_name).unwrap();
                    frac_convective_cool.insert(c_name.into(), space_cool_system.frac_convective());
                    temp_setpnt_cool.insert(
                        c_name.into(),
                        space_cool_system
                            .temp_setpnt(&simtime)
                            .unwrap_or_else(temp_setpnt_cool_none),
                    );
                }
            }
        }

        SetpointsAndConvectiveFractions {
            temp_setpnt_heat,
            temp_setpnt_cool,
            frac_convective_heat,
            frac_convective_cool,
        }
    }

    fn gains_heat_cool(
        &self,
        delta_t_h: f64,
        hc_output_convective: &IndexMap<String, f64>,
        hc_output_radiative: &IndexMap<String, f64>,
    ) -> (f64, f64) {
        let gains_heat_cool_convective =
            hc_output_convective.values().sum::<f64>() * WATTS_PER_KILOWATT as f64 / delta_t_h;
        let gains_heat_cool_radiative =
            hc_output_radiative.values().sum::<f64>() * WATTS_PER_KILOWATT as f64 / delta_t_h;

        (gains_heat_cool_convective, gains_heat_cool_radiative)
    }

    /// Calculate the incoming air changes per hour
    /// initial_p_z_ref_guess is used for calculation in first timestep.
    /// Later timesteps use the previous timesteps p_z_ref of max and min ACH ,respective to calc.
    fn calc_air_changes_per_hour(
        &self,
        wind_speed: f64,
        wind_direction: f64,
        temp_int_air: f64,
        temp_ext_air: f64,
        r_v_arg: f64,
        r_w_arg: Option<f64>,
        initial_p_z_ref_guess: f64,
        reporting_flag: ReportingFlag,
        internal_pressure_window: &mut HashMap<ReportingFlag, f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let current_internal_pressure_window = if self.initial_loop.load(Ordering::SeqCst) {
            self.ventilation.calculate_internal_reference_pressure(
                initial_p_z_ref_guess,
                wind_speed,
                wind_direction,
                temp_int_air,
                temp_ext_air,
                r_v_arg,
                r_w_arg,
                simtime,
            )?
        } else {
            self.ventilation.calculate_internal_reference_pressure(
                internal_pressure_window[&reporting_flag],
                wind_speed,
                wind_direction,
                temp_int_air,
                temp_ext_air,
                r_v_arg,
                r_w_arg,
                simtime,
            )?
        };

        internal_pressure_window.insert(reporting_flag, current_internal_pressure_window);

        let incoming_air_flow = self.ventilation.incoming_air_flow(
            current_internal_pressure_window,
            wind_speed,
            wind_direction,
            temp_int_air,
            temp_ext_air,
            r_v_arg,
            r_w_arg,
            Some(reporting_flag),
            Some(true),
            simtime,
        )?;

        Ok(incoming_air_flow / self.total_volume)
    }

    /// Get minimum output for each heating/cooling system in the specified zone
    fn heat_cool_system_output_min(
        &self,
        h_name_list_sorted_zone: &HashMap<&str, Vec<String>>,
        c_name_list_sorted_zone: &HashMap<&str, Vec<String>>,
        frac_convective_heat_zone_system: &HashMap<&str, IndexMap<String, f64>>,
        frac_convective_cool_zone_system: &HashMap<&str, IndexMap<String, f64>>,
        z_name: &str,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<HeatCoolOutputs> {
        let h_output_min: IndexMap<String, f64> = h_name_list_sorted_zone[z_name]
            .iter()
            .filter(|h_name| h_name.as_str() != "") // we need to exclude the empty string as it stands for None (yes, we're stringly typing this)
            .map(|h_name| -> anyhow::Result<(String, f64)> {
                Ok((
                    h_name.clone(),
                    self.space_heat_systems
                        .get(h_name)
                        .unwrap()
                        .lock()
                        .energy_output_min(simulation_time_iteration)?,
                ))
            })
            .try_collect()?;
        let c_output_min = c_name_list_sorted_zone[z_name]
            .iter()
            .filter(|c_name| c_name.as_str() != "") // we need to exclude the empty string as it stands for None (yes, we're stringly typing this)
            .map(|c_name| {
                (
                    c_name.clone(),
                    self.space_cool_systems
                        .get(c_name)
                        .unwrap()
                        .energy_output_min(),
                )
            })
            .collect::<IndexMap<_, _>>();

        let mut hc_output_min = h_output_min.clone();
        hc_output_min.extend(c_output_min);
        hc_output_min.extend(IndexMap::from([("".into(), 0.0)])); // empty string used here as equivalent of None in Python

        let mut frac_convective_system = frac_convective_heat_zone_system[z_name].clone();
        frac_convective_system.extend(frac_convective_cool_zone_system[z_name].clone());

        let hc_output_convective = h_name_list_sorted_zone[z_name]
            .iter()
            .chain(c_name_list_sorted_zone[z_name].iter())
            .map(|hc_name| {
                (
                    hc_name.clone(),
                    hc_output_min[hc_name.as_str()] * frac_convective_system[hc_name],
                )
            })
            .collect::<IndexMap<_, _>>();
        let hc_output_radiative = h_name_list_sorted_zone[z_name]
            .iter()
            .chain(c_name_list_sorted_zone[z_name].iter())
            .map(|hc_name| {
                (
                    hc_name.clone(),
                    hc_output_min[hc_name] - hc_output_convective[hc_name],
                )
            })
            .collect::<IndexMap<_, _>>();

        Ok(HeatCoolOutputs {
            hc_output_convective,
            hc_output_radiative,
            hc_output_min,
        })
    }

    /// Calculate space heating demand, heating system output and temperatures
    ///
    /// Arguments:
    /// * `delta_t_h` - calculation timestep, in hours
    /// * `gains_internal_dhw` - internal gains from hot water system for this timestep, in W
    fn calc_space_heating(
        &self,
        delta_t_h: f64,
        gains_internal_dhw: f64,
        internal_pressure_window: &mut HashMap<ReportingFlag, f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<SpaceHeatingCalculation> {
        let wind_speed = self.external_conditions.wind_speed(&simtime);
        let wind_direction = self.external_conditions.wind_direction(simtime);
        let temp_ext_air = self.external_conditions.air_temp(&simtime);
        let temp_int_air = self.temp_internal_air_prev_timestep();
        let ach_min = self
            .vent_adjust_min_control
            .as_ref()
            .and_then(|c| c.setpnt(&simtime));
        let ach_max = self
            .vent_adjust_max_control
            .as_ref()
            .and_then(|c| c.setpnt(&simtime));
        // Calculate timestep in seconds
        let delta_t = delta_t_h * SECONDS_PER_HOUR as f64;

        let internal_gains_ductwork = self.calc_internal_gains_ductwork(simtime);
        let internal_gains_ductwork_per_m3 = internal_gains_ductwork / self.total_volume;

        let internal_gains_buffer_tank = self.calc_internal_gains_buffer_tank();

        // Adjust vent position if required to attempt to meet min or max ach
        self.r_v_arg.store(
            self.ventilation.find_r_v_arg_within_bounds(
                ach_min,
                ach_max,
                self.r_v_arg.load(Ordering::SeqCst),
                wind_speed,
                wind_direction,
                temp_int_air,
                temp_ext_air,
                Some(0.),
                0.,
                None,
                simtime,
            )?,
            Ordering::SeqCst,
        );

        // Windows shut
        let ach_windows_shut = self.calc_air_changes_per_hour(
            wind_speed,
            wind_direction,
            temp_int_air,
            temp_ext_air,
            self.r_v_arg.load(Ordering::SeqCst),
            Some(0.),
            0.,
            ReportingFlag::Min,
            internal_pressure_window,
            simtime,
        )?;

        // Windows fully open
        let ach_windows_open = self.calc_air_changes_per_hour(
            wind_speed,
            wind_direction,
            temp_int_air,
            temp_ext_air,
            self.r_v_arg.load(Ordering::SeqCst),
            Some(1.),
            0.,
            ReportingFlag::Max,
            internal_pressure_window,
            simtime,
        )?;

        // To indicate the future loop should involve the p_Z_ref from previous calc
        self.initial_loop.store(false, Ordering::SeqCst);

        let ach_target = if let Some(required_vent_data) = self.required_vent_data.as_ref() {
            let ach_target = required_vent_data.schedule[simtime.time_series_idx(
                required_vent_data.start_day,
                required_vent_data.time_series_step,
            )];

            ach_windows_shut.max(ach_target.unwrap_or(ach_windows_open).min(ach_windows_open))
        } else {
            ach_windows_shut
        };

        let mut gains_internal_zone: HashMap<String, f64> = Default::default();
        let mut gains_solar_zone: HashMap<String, f64> = Default::default();
        let mut h_name_list_sorted_zone: HashMap<&str, Vec<String>> = Default::default();
        let mut c_name_list_sorted_zone: HashMap<&str, Vec<String>> = Default::default();
        let mut temp_setpnt_heat_zone_system: HashMap<&str, IndexMap<String, f64>> =
            Default::default();
        let mut temp_setpnt_cool_zone_system: HashMap<&str, IndexMap<String, f64>> =
            Default::default();
        let mut frac_convective_heat_zone_system: HashMap<&str, IndexMap<String, f64>> =
            Default::default();
        let mut frac_convective_cool_zone_system: HashMap<&str, IndexMap<String, f64>> =
            Default::default();
        let mut ach_cooling_zone: HashMap<&str, f64> = Default::default();
        let mut ach_to_trigger_heating_zone: HashMap<&str, Option<f64>> = Default::default();
        let mut internal_air_temp: HashMap<String, f64> = Default::default();
        let mut operative_temp: HashMap<String, f64> = Default::default();
        let mut space_heat_demand_zone: HashMap<String, f64> = Default::default();
        let mut space_cool_demand_zone: HashMap<String, f64> = Default::default();
        let mut space_heat_provided_system: HashMap<String, f64> = Default::default();
        let mut space_cool_provided_system: HashMap<String, f64> = Default::default();
        let mut heat_balance_map: HashMap<String, Option<HeatBalance>> = Default::default();

        // Average supply temperature
        let avg_air_supply_temp = self.external_conditions.air_temp(&simtime);

        for (z_name, zone) in self.zones.iter() {
            let z_name = z_name.as_str();
            // Calculate internal and solar gains
            gains_internal_zone.insert(
                z_name.into(),
                self.space_heat_internal_gains_for_zone(
                    zone,
                    gains_internal_dhw,
                    internal_gains_ductwork_per_m3,
                    internal_gains_buffer_tank,
                    simtime,
                )?,
            );
            gains_solar_zone.insert(z_name.into(), zone.gains_solar(simtime));

            // Get heating and cooling characteristics for the current zone
            let (
                h_name_list_sorted_zone_current,
                c_name_list_sorted_zone_current,
                SetpointsAndConvectiveFractions {
                    temp_setpnt_heat: temp_setpnt_heat_zone_system_current,
                    temp_setpnt_cool: temp_setpnt_cool_zone_system_current,
                    frac_convective_heat: frac_convective_heat_zone_system_current,
                    frac_convective_cool: frac_convective_cool_zone_system_current,
                },
            ) = self.heat_cool_systems_for_zone(z_name, simtime);

            h_name_list_sorted_zone.insert(z_name, h_name_list_sorted_zone_current);
            c_name_list_sorted_zone.insert(z_name, c_name_list_sorted_zone_current);
            temp_setpnt_heat_zone_system.insert(z_name, temp_setpnt_heat_zone_system_current);
            temp_setpnt_cool_zone_system.insert(z_name, temp_setpnt_cool_zone_system_current);
            frac_convective_heat_zone_system
                .insert(z_name, frac_convective_heat_zone_system_current);
            frac_convective_cool_zone_system
                .insert(z_name, frac_convective_cool_zone_system_current);

            // Calculate space heating demand based on highest-priority systems,
            // assuming no output from any other systems

            let (
                space_heat_demand_zone_current,
                space_cool_demand_zone_current,
                ach_cooling_zone_current,
                ach_to_trigger_heating_zone_current,
            ) = zone.space_heat_cool_demand(
                delta_t_h,
                temp_ext_air,
                gains_internal_zone[z_name],
                gains_solar_zone[z_name],
                frac_convective_heat_zone_system[z_name][&h_name_list_sorted_zone[z_name][0]],
                frac_convective_cool_zone_system[z_name][&c_name_list_sorted_zone[z_name][0]],
                temp_setpnt_heat_zone_system[z_name][&h_name_list_sorted_zone[z_name][0]],
                temp_setpnt_cool_zone_system[z_name][&c_name_list_sorted_zone[z_name][0]],
                avg_air_supply_temp,
                None,
                None,
                AirChangesPerHourArgument::TargetAndWindowsOpen {
                    ach_target,
                    ach_windows_open,
                },
                simtime,
            )?;

            space_heat_demand_zone.insert(z_name.into(), space_heat_demand_zone_current);
            space_cool_demand_zone.insert(z_name.into(), space_cool_demand_zone_current);
            ach_cooling_zone.insert(z_name, ach_cooling_zone_current);
            ach_to_trigger_heating_zone.insert(z_name, ach_to_trigger_heating_zone_current);
        }

        // Ventilation required, including for cooling
        let is_heating_demand = space_heat_demand_zone.values().any(|&demand| demand > 0.0);
        let is_cooling_demand = space_cool_demand_zone.values().any(|&demand| demand < 0.0);

        let ach_cooling = if is_heating_demand {
            // Do not open windows any further than required for ventilation
            // requirement if there is any heating demand in any zone
            ach_target
        } else if is_cooling_demand {
            let ach_cooling = ach_target;

            // In this case, will need to recalculate space cooling demand for
            // each zone, this time assuming no window opening for all zones
            // TODO (from Python) There might be a way to make this more efficient and reduce
            //      the number of times the heat balance solver has to run, but
            //      this would require a wider refactoring of the zone module's
            //      space_heat_cool_demand function
            for (z_name, zone) in self.zones.iter() {
                let z_name = z_name.as_str();
                let (space_heat_demand_zone_current, space_cool_demand_zone_current, _, _) = zone
                    .space_heat_cool_demand(
                    delta_t_h,
                    temp_ext_air,
                    gains_internal_zone[z_name],
                    gains_solar_zone[z_name],
                    frac_convective_heat_zone_system[z_name][&h_name_list_sorted_zone[z_name][0]],
                    frac_convective_cool_zone_system[z_name][&c_name_list_sorted_zone[z_name][0]],
                    temp_setpnt_heat_zone_system[z_name][&h_name_list_sorted_zone[z_name][0]],
                    temp_setpnt_cool_zone_system[z_name][&c_name_list_sorted_zone[z_name][0]],
                    avg_air_supply_temp,
                    None,
                    None,
                    AirChangesPerHourArgument::Cooling { ach_cooling },
                    simtime,
                )?;
                space_heat_demand_zone.insert(z_name.into(), space_heat_demand_zone_current);
                space_cool_demand_zone.insert(z_name.into(), space_cool_demand_zone_current);
            }

            ach_cooling
        } else {
            // Subject to the above/below limits, take the maximum required window
            // opening from across all the zones
            let mut ach_cooling = *ach_cooling_zone
                .values()
                .max_by(|a, b| a.total_cmp(b))
                .unwrap();

            // Do not open windows to an extent where it would cause any zone
            // temperature to fall below the heating setpoint for that zone
            let ach_to_trigger_heating_list: Vec<f64> = ach_to_trigger_heating_zone
                .values()
                .filter_map(|x| *x)
                .collect();
            if !ach_to_trigger_heating_list.is_empty() {
                ach_cooling = ach_to_trigger_heating_list
                    .iter()
                    .max_by(|a, b| a.total_cmp(b).reverse())
                    .unwrap()
                    .min(ach_cooling);
            }

            // Do not reduce air change rate below ventilation requirement even
            // if it would help with temperature regulation

            ach_cooling.max(ach_target)
        };

        // Calculate heating/cooling system response and temperature achieved in each zone
        for (z_name, zone) in self.zones.iter() {
            let z_name = z_name.as_str();
            let c_name_list_hashset = c_name_list_sorted_zone[z_name]
                .iter()
                .filter(|name| !name.is_empty()) // Ignore the "" placeholder to ensure zones with no system don't trigger the uniqueness error
                .collect::<HashSet<_>>();
            let h_name_list_hashset = h_name_list_sorted_zone[z_name]
                .iter()
                .filter(|name| !name.is_empty()) // Ignore the "" placeholder to ensure zones with no system don't trigger the uniqueness error
                .collect::<HashSet<_>>();
            let intersection = h_name_list_hashset
                .intersection(&c_name_list_hashset)
                .collect::<Vec<_>>();
            if !intersection.is_empty() {
                bail!("All heating and cooling systems must have unique names")
            };
            // drop temporary values for working out intersection from scope
            drop(intersection);
            drop(c_name_list_hashset);
            drop(h_name_list_hashset);

            // Initialise system outputs to minimum output for each heating and cooling system
            let HeatCoolOutputs {
                mut hc_output_convective,
                mut hc_output_radiative,
                hc_output_min,
            } = self.heat_cool_system_output_min(
                &h_name_list_sorted_zone,
                &c_name_list_sorted_zone,
                &frac_convective_heat_zone_system,
                &frac_convective_cool_zone_system,
                z_name,
                simtime,
            )?;
            let mut space_heat_demand_zone_system = h_name_list_sorted_zone[z_name]
                .iter()
                .map(|h_name| (h_name.as_str(), 0.0))
                .collect::<IndexMap<_, _>>();
            let mut space_cool_demand_zone_system = c_name_list_sorted_zone[z_name]
                .iter()
                .map(|c_name| (c_name.as_str(), 0.0))
                .collect::<IndexMap<_, _>>();
            let mut space_heat_provided_zone_system = h_name_list_sorted_zone[z_name]
                .iter()
                .map(|h_name| (h_name.as_str(), 0.0))
                .collect::<IndexMap<_, _>>();
            let mut space_cool_provided_zone_system = c_name_list_sorted_zone[z_name]
                .iter()
                .map(|c_name| (c_name.as_str(), 0.0))
                .collect::<IndexMap<_, _>>();

            let mut h_idx: usize = 0;
            let mut c_idx: usize = 0;
            let _space_heat_running_time_cumulative = 0.0;

            // these following variables need to be referenced outside the while loop with their last value retained

            let mut frac_convective_heat = 1.0;
            let mut frac_convective_cool = 1.0;

            while h_idx < h_name_list_sorted_zone[z_name].len()
                && c_idx < c_name_list_sorted_zone[z_name].len()
            {
                let h_name = &h_name_list_sorted_zone[z_name][h_idx].as_str();
                let c_name = &c_name_list_sorted_zone[z_name][c_idx].as_str();
                frac_convective_heat = frac_convective_heat_zone_system[z_name][h_name.to_owned()];
                frac_convective_cool = frac_convective_cool_zone_system[z_name][c_name.to_owned()];
                let temp_setpnt_heat = temp_setpnt_heat_zone_system[z_name][h_name.to_owned()];
                let temp_setpnt_cool = temp_setpnt_cool_zone_system[z_name][c_name.to_owned()];

                // Calculate space heating/cooling demand, accounting for any
                // output from systems (either output already calculated for
                // higher-priority systems, or min output of current or
                // lower-priority systems).
                let (gains_heat_cool_convective, gains_heat_cool_radiative) =
                    self.gains_heat_cool(delta_t_h, &hc_output_convective, &hc_output_radiative);
                if gains_heat_cool_convective == 0.0 && gains_heat_cool_radiative == 0.0 {
                    // If there is no output from any systems, then don't need to
                    // calculate demand again
                    space_heat_demand_zone_system.insert(h_name, space_heat_demand_zone[z_name]);
                    space_cool_demand_zone_system.insert(c_name, space_cool_demand_zone[z_name]);
                } else {
                    let (
                        space_heat_demand_zone_system_current,
                        space_cool_demand_zone_system_current,
                        ach_cooling_zone_current,
                        _,
                    ) = zone.space_heat_cool_demand(
                        delta_t_h,
                        temp_ext_air,
                        gains_internal_zone[z_name],
                        gains_solar_zone[z_name],
                        frac_convective_heat,
                        frac_convective_cool,
                        temp_setpnt_heat,
                        temp_setpnt_cool,
                        avg_air_supply_temp,
                        Some(gains_heat_cool_convective),
                        Some(gains_heat_cool_radiative),
                        AirChangesPerHourArgument::Cooling { ach_cooling },
                        simtime,
                    )?;
                    space_heat_demand_zone_system
                        .insert(h_name, space_heat_demand_zone_system_current);
                    space_cool_demand_zone_system
                        .insert(c_name, space_cool_demand_zone_system_current);
                    ach_cooling_zone.insert(z_name, ach_cooling_zone_current);

                    // Space heating/cooling demand calculated above already assumes
                    // minimum output from all systems, so we need to add this on
                    // for the current system. If minimum output is enough to meet
                    // demand, then demand calculated above will be zero. This is
                    // okay because calling the demand_energy function for the
                    // heating/cooling system later with an input of zero will still
                    // result in the minimum output being provided.
                    if space_heat_demand_zone_system[h_name] > 0. {
                        space_heat_demand_zone_system[h_name] += hc_output_min[h_name.to_owned()]
                    } else if space_cool_demand_zone_system[c_name] < 0. {
                        space_cool_demand_zone_system[c_name] += hc_output_min[c_name.to_owned()]
                    }
                }

                // If any heating systems potentially require overventilation,
                // calculate running time and throughput factor for current service
                // based on space heating demand assuming only overventilation
                // required for DHW

                // NB. there is a large chunk of Python in the upstream code here that is commented out there currently, so not brought over

                // Calculate heating/cooling provided
                if space_heat_demand_zone_system[h_name] > 0.0 {
                    space_heat_provided_zone_system.insert(
                        h_name,
                        self.space_heat_systems[h_name.to_owned()]
                            .lock()
                            .demand_energy(space_heat_demand_zone_system[h_name], simtime)?,
                    );
                    hc_output_convective.insert(
                        (*h_name).into(),
                        space_heat_provided_zone_system[h_name] * frac_convective_heat,
                    );
                    hc_output_radiative.insert(
                        (*h_name).into(),
                        space_heat_provided_zone_system[h_name] * (1.0 - frac_convective_heat),
                    );
                    // If heating has been provided, then next iteration of loop
                    // should use next-priority heating system
                    h_idx += 1;
                }
                if space_cool_demand_zone_system[c_name] < 0.0 {
                    space_cool_provided_zone_system.insert(
                        c_name,
                        self.space_cool_systems[c_name.to_owned()]
                            .demand_energy(space_cool_demand_zone_system[c_name], simtime),
                    );
                    hc_output_convective.insert(
                        (*c_name).into(),
                        space_cool_provided_zone_system[c_name] * frac_convective_cool,
                    );
                    hc_output_radiative.insert(
                        (*c_name).into(),
                        space_cool_provided_zone_system[c_name] * (1.0 - frac_convective_cool),
                    );
                    // If cooling has been provided, then next iteration of loop
                    // should use next-priority cooling system
                    c_idx += 1;
                }

                // Terminate loop if there is no more demand
                if space_heat_demand_zone_system[h_name] <= 0.0
                    && space_cool_demand_zone_system[c_name] >= 0.0
                {
                    break;
                }
            }

            // Call any remaining heating and cooling systems with zero demand
            for h_name in h_name_list_sorted_zone[z_name][h_idx..].iter() {
                let h_name = h_name.as_str();
                space_heat_provided_zone_system.insert(
                    h_name,
                    if !h_name.is_empty() {
                        self.space_heat_systems[h_name]
                            .lock()
                            .demand_energy(0.0, simtime)?
                    } else {
                        0.0
                    },
                );
                hc_output_convective.insert(
                    (*h_name).into(),
                    space_heat_provided_zone_system[h_name] * frac_convective_heat,
                );
                hc_output_radiative.insert(
                    (*h_name).into(),
                    space_heat_provided_zone_system[h_name] * (1.0 - frac_convective_heat),
                );
            }
            for c_name in c_name_list_sorted_zone[z_name][c_idx..].iter() {
                let c_name = c_name.as_str();
                space_cool_provided_zone_system.insert(
                    c_name,
                    if !c_name.is_empty() {
                        self.space_cool_systems[c_name].demand_energy(0.0, simtime)
                    } else {
                        0.0
                    },
                );
                hc_output_convective.insert(
                    (*c_name).into(),
                    space_cool_provided_zone_system[c_name] * frac_convective_cool,
                );
                hc_output_radiative.insert(
                    (*c_name).into(),
                    space_cool_provided_zone_system[c_name] * (1.0 - frac_convective_cool),
                );
            }
            // Calculate unmet demand
            self.unmet_demand(
                delta_t_h,
                temp_ext_air,
                z_name,
                zone,
                gains_internal_zone[z_name],
                gains_solar_zone[z_name],
                &temp_setpnt_heat_zone_system[z_name],
                &temp_setpnt_cool_zone_system[z_name],
                &frac_convective_heat_zone_system[z_name],
                &frac_convective_cool_zone_system[z_name],
                &h_name_list_sorted_zone[z_name],
                &c_name_list_sorted_zone[z_name],
                space_heat_demand_zone[z_name],
                space_cool_demand_zone[z_name],
                &hc_output_convective,
                &hc_output_radiative,
                ach_windows_open,
                ach_target,
                avg_air_supply_temp,
                simtime,
            )?;

            // Sum heating gains (+ve) and cooling gains (-ve) and convert from kWh to W
            let hc_output_convective_total = hc_output_convective.values().sum::<f64>();
            let hc_output_radiative_total = hc_output_radiative.values().sum::<f64>();

            let gains_heat_cool = (hc_output_convective_total + hc_output_radiative_total)
                * WATTS_PER_KILOWATT as f64
                / delta_t_h;
            let frac_convective = if gains_heat_cool != 0.0 {
                hc_output_convective_total
                    / (hc_output_convective_total + hc_output_radiative_total)
            } else {
                1.0
            };

            // Calculate final temperatures achieved
            heat_balance_map.insert(
                z_name.into(),
                zone.update_temperatures(
                    delta_t,
                    temp_ext_air,
                    gains_internal_zone[z_name],
                    gains_solar_zone[z_name],
                    gains_heat_cool,
                    frac_convective,
                    ach_cooling,
                    avg_air_supply_temp,
                    simtime,
                )?,
            );
            internal_air_temp.insert(z_name.into(), zone.temp_internal_air());
            operative_temp.insert(z_name.into(), zone.temp_operative());

            for h_name in h_name_list_sorted_zone[z_name].iter() {
                *space_heat_provided_system
                    .entry(h_name.to_owned())
                    .or_insert(0.0) += space_heat_provided_zone_system[h_name.as_str()];
            }
            for c_name in c_name_list_sorted_zone[z_name].iter() {
                *space_cool_provided_system
                    .entry(c_name.to_owned())
                    .or_insert(0.0) += space_cool_provided_zone_system[c_name.as_str()];
            }
        }

        Ok(SpaceHeatingCalculation {
            gains_internal_zone,
            gains_solar_zone,
            operative_temp,
            internal_air_temp,
            space_heat_demand_zone,
            space_cool_demand_zone,
            space_heat_provided_system,
            space_cool_provided_system,
            internal_gains_ductwork,
            heat_balance_map,
        })
    }

    /// Determine highest-priority system that is in its required heating
    /// or cooling period (and not just the setback period)
    ///
    /// Arguments:
    /// * `hc_name_list_sorted` - list of heating or cooling systems (not combined
    ///                                list), sorted in order of priority
    /// * `space_heat_cool_systems` - dict of space heating or cooling system objects
    ///                                   (not combined list)
    fn highest_priority_required_system(
        &self,
        hc_name_list_sorted: &[String],
        space_heat_cool_systems: SpaceHeatCoolSystems,
        simtime: SimulationTimeIteration,
    ) -> Option<String> {
        let mut hc_name_highest_req = Default::default();
        for hc_name in hc_name_list_sorted {
            if !hc_name.is_empty()
                && space_heat_cool_systems
                    .in_required_period_for_name(hc_name, simtime)
                    .unwrap_or(false)
            {
                hc_name_highest_req = Some(hc_name.to_owned());
                break;
            }
        }

        hc_name_highest_req
    }

    /// Calculate how much space heating / cooling demand is unmet
    fn unmet_demand(
        &self,
        delta_t_h: f64,
        temp_ext_air: f64,
        z_name: &str,
        zone: &Zone,
        gains_internal: f64,
        gains_solar: f64,
        temp_setpnt_heat_system: &IndexMap<String, f64>,
        temp_setpnt_cool_system: &IndexMap<String, f64>,
        frac_convective_heat_system: &IndexMap<String, f64>,
        frac_convective_cool_system: &IndexMap<String, f64>,
        h_name_list_sorted: &[String],
        c_name_list_sorted: &[String],
        space_heat_demand: f64,
        space_cool_demand: f64,
        hc_output_convective: &IndexMap<String, f64>,
        hc_output_radiative: &IndexMap<String, f64>,
        ach_max: f64,
        ach_target: f64,
        avg_air_supply_temp: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<()> {
        // Note: Use demand calculated based on highest-priority systems
        // Note: Demand is not considered unmet if it is outside the
        //       required heating/cooling period (which does not include
        //       times when the system is on due to setback or advanced
        //       start). If different systems have different required
        //       heating/cooling periods, unmet demand will be based on the
        //       system with the highest setpoint, ignoring any systems that
        //       are not in required periods (e.g. systems that are in
        //       setback or advanced start periods).
        // Note: Need to check that demand is non-zero, to avoid
        //       reporting unmet demand when heating system is absorbing
        //       energy from zone or cooling system is releasing energy
        //       to zone, which may be the case in some timesteps for
        //       systems with significant thermal mass.

        // Determine highest-priority system that is in its required heating
        // or cooling period (and not just the setback period)
        let h_name_highest_req = self.highest_priority_required_system(
            h_name_list_sorted,
            SpaceHeatCoolSystems::Heat(&self.space_heat_systems),
            simtime,
        );
        let c_name_highest_req = self.highest_priority_required_system(
            c_name_list_sorted,
            SpaceHeatCoolSystems::Cool(&self.space_cool_systems),
            simtime,
        );

        let gains_heat = h_name_list_sorted
            .iter()
            .map(|h_name| {
                hc_output_convective[h_name.as_str()] + hc_output_radiative[h_name.as_str()]
            })
            .sum::<f64>();
        let gains_cool = c_name_list_sorted
            .iter()
            .map(|c_name| {
                hc_output_convective[c_name.as_str()] + hc_output_radiative[c_name.as_str()]
            })
            .sum::<f64>();
        let energy_shortfall_heat = 0.0f64.max(space_heat_demand - gains_heat);
        let energy_shortfall_cool = 0.0f64.max(-(space_cool_demand - gains_cool));

        if (h_name_highest_req.is_some() && space_heat_demand > 0.0 && energy_shortfall_heat > 0.0)
            || (c_name_highest_req.is_some()
                && space_cool_demand < 0.0
                && energy_shortfall_cool > 0.0)
        {
            let (unmet_demand_heat, unmet_demand_cool) = if (energy_shortfall_heat > 0.0
                && h_name_highest_req
                    .as_ref()
                    .is_some_and(|h_name| h_name != &h_name_list_sorted[0]))
                || (energy_shortfall_cool > 0.0
                    && c_name_highest_req
                        .as_ref()
                        .is_some_and(|c_name| c_name != &c_name_list_sorted[0]))
            {
                // If the highest-priority system is not in required heating
                // period, but a lower-priority system is, calculate demand
                // based on the highest-priority system that is in required
                // heating period
                //
                // Handle case where no heating/cooling system is in required
                // period. In this case, there will be no heat output anyway so
                // the convective fraction doesn't matter
                let (frac_convective_heat, temp_setpnt_heat) =
                    if let Some(h_name) = h_name_highest_req.as_ref() {
                        (
                            frac_convective_heat_system[h_name],
                            temp_setpnt_heat_system[h_name],
                        )
                    } else {
                        (1.0, temp_setpnt_heat_none())
                    };
                let (frac_convective_cool, temp_setpnt_cool) =
                    if let Some(c_name) = c_name_highest_req.as_ref() {
                        (
                            frac_convective_cool_system[c_name],
                            temp_setpnt_cool_system[c_name],
                        )
                    } else {
                        (1.0, temp_setpnt_cool_none())
                    };

                let (space_heat_demand_req, space_cool_demand_req, _, _) = zone
                    .space_heat_cool_demand(
                        delta_t_h,
                        temp_ext_air,
                        gains_internal,
                        gains_solar,
                        frac_convective_heat,
                        frac_convective_cool,
                        temp_setpnt_heat,
                        temp_setpnt_cool,
                        avg_air_supply_temp,
                        None,
                        None,
                        AirChangesPerHourArgument::TargetAndWindowsOpen {
                            ach_target,
                            ach_windows_open: ach_max,
                        },
                        simtime,
                    )?;

                (
                    0.0f64.max(space_heat_demand_req - gains_heat),
                    0.0f64.max(-(space_cool_demand_req - gains_cool)),
                )
            } else {
                // If highest-priority system is in required heating period,
                // use the demand already calculated for the zone
                (energy_shortfall_heat, energy_shortfall_cool)
            };

            self.energy_supply_conn_unmet_demand_zone[z_name]
                .demand_energy(unmet_demand_heat + unmet_demand_cool, simtime.index)?;
        }

        Ok(())
    }

    pub fn run(&self) -> anyhow::Result<RunResults> {
        let simulation_time = self.simulation_time.as_ref().to_owned();
        let vec_capacity = || Vec::with_capacity(simulation_time.total_steps());

        let mut timestep_array = vec_capacity();
        let mut gains_internal_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut gains_solar_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut operative_temp_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut internal_air_temp_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut space_heat_demand_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut space_cool_demand_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut space_heat_provided_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut space_cool_provided_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut zone_list: Vec<String> = Default::default();
        let mut hot_water_demand_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut hot_water_energy_demand_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut hot_water_energy_demand_dict_incl_pipework: IndexMap<String, Vec<f64>> =
            Default::default();
        let mut hot_water_energy_output_dict: IndexMap<&str, Vec<f64>> = Default::default();
        let mut hot_water_duration_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut hot_water_no_events_dict: IndexMap<String, Vec<usize>> = Default::default();
        let mut hot_water_pipework_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut ductwork_gains_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut hot_water_primary_pipework_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut hot_water_storage_losses_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut heat_balance_all_dict: HeatBalanceAllResults = IndexMap::from([
            (HeatBalanceFieldName::AirNode, Default::default()),
            (HeatBalanceFieldName::InternalBoundary, Default::default()),
            (HeatBalanceFieldName::ExternalBoundary, Default::default()),
        ]);
        let mut heat_source_wet_results_dict: IndexMap<String, ResultsPerTimestep> =
            Default::default();
        let mut heat_source_wet_results_annual_dict: IndexMap<String, ResultsAnnual> =
            Default::default();
        let mut emitters_output_dict: IndexMap<String, Vec<EmittersDetailedResult>> =
            Default::default();
        let mut vent_output_list: Vec<VentilationDetailedResult> = Default::default();
        let mut esh_output_dict: IndexMap<String, Vec<StorageHeaterDetailedResult>> =
            Default::default();
        let mut hot_water_source_results_dict: IndexMap<String, Vec<StorageTankDetailedResult>> =
            Default::default();

        for z_name in self.zones.keys() {
            gains_internal_dict.insert(z_name.clone(), vec_capacity());
            gains_solar_dict.insert(z_name.clone(), vec_capacity());
            operative_temp_dict.insert(z_name.clone(), vec_capacity());
            internal_air_temp_dict.insert(z_name.clone(), vec_capacity());
            space_heat_demand_dict.insert(z_name.clone(), vec_capacity());
            space_cool_demand_dict.insert(z_name.clone(), vec_capacity());
            zone_list.push(z_name.clone());
            for heat_balance_value in heat_balance_all_dict.values_mut() {
                heat_balance_value.insert(z_name.clone(), Default::default());
            }
        }

        for z_h_names in self.heat_system_name_for_zone.values() {
            for h_name in z_h_names {
                space_heat_provided_dict.insert(h_name.to_owned(), vec_capacity());
            }
        }

        for z_c_names in self.cool_system_name_for_zone.values() {
            for c_name in z_c_names {
                space_cool_provided_dict.insert(c_name.to_owned(), vec_capacity());
            }
        }

        hot_water_demand_dict.insert("demand".into(), vec_capacity());
        hot_water_energy_demand_dict.insert("energy_demand".into(), vec_capacity());
        hot_water_energy_demand_dict_incl_pipework
            .insert("energy_demand_incl_pipework_loss".into(), vec_capacity());
        hot_water_energy_output_dict.insert("energy_output", vec_capacity());
        hot_water_duration_dict.insert("duration".into(), vec_capacity());
        hot_water_no_events_dict.insert(
            "no_events".into(),
            Vec::with_capacity(simulation_time.total_steps()),
        );
        hot_water_pipework_dict.insert("pw_losses".into(), vec_capacity());
        ductwork_gains_dict.insert("ductwork_gains".into(), vec_capacity());
        hot_water_primary_pipework_dict.insert("primary_pw_losses".into(), vec_capacity());
        hot_water_storage_losses_dict.insert("storage_losses".into(), vec_capacity());
        self.initial_loop.store(true, Ordering::SeqCst);
        let mut internal_pressure_window: HashMap<ReportingFlag, f64> = Default::default();

        let delta_t_h = simulation_time.step_in_hours();

        #[cfg(feature = "indicatif")]
        let simulation_time_iter = simulation_time.progress();
        #[cfg(not(feature = "indicatif"))]
        let simulation_time_iter = simulation_time;

        for t_it in simulation_time_iter {
            timestep_array.push(t_it.time);
            self.update_temp_internal_air();
            let temp_hot_water = self.hot_water_sources["hw cylinder"].temp_hot_water()?;
            let _temp_final_drawoff = temp_hot_water;
            let _temp_average_drawoff = temp_hot_water;

            let DomesticHotWaterDemandData {
                hw_demand_vol,
                hw_demand_vol_target,
                hw_vol_at_tapping_points,
                hw_duration,
                all_events: no_events,
                hw_energy_demand,
                usage_events: _,
                vol_hot_water_equiv_elec_shower: _,
            } = self
                .domestic_hot_water_demand
                .hot_water_demand(t_it, temp_hot_water)?;

            // Running heat sources of pre-heated tanks and updating thermal losses, etc.
            for source in self.pre_heated_water_sources.values() {
                match source {
                    HotWaterStorageTank::StorageTank(storage_tank) => {
                        storage_tank.read().demand_hot_water(None, t_it)?;
                    }
                    HotWaterStorageTank::SmartHotWaterTank(smart_storage_tank) => {
                        smart_storage_tank.read().demand_hot_water(None, t_it)?;
                    }
                }
            }
            // TODO (Python) Remove hard-coding of hot water source name
            // TODO (Python) Reporting of the hot water energy output assumes that there
            //      is only one water heating system. If the model changes in
            //      future to allow more than one hot water system, this code may
            //      need to be revised to handle that scenario.

            let (hw_energy_output, pw_losses_internal, pw_losses_external, gains_internal_dhw_use) =
                if let HotWaterSource::PreHeated(_source) = &self.hot_water_sources["hw cylinder"] {
                    unimplemented!("To be implemented as part of migration to 1_0_a1")
                } else if let HotWaterSource::HeatBattery(source) =
                    &self.hot_water_sources["hw cylinder"]
                {
                    let hw_energy_output = source.demand_hot_water(None, t_it)?; // TODO 1.0.0a1 update to use common demand_hot_water method once its signature has been updated and pass in usage_events
                    let (pw_losses_internal, pw_losses_external, gains_internal_dhw_use) = self
                        .pipework_losses_and_internal_gains_from_hw(
                            delta_t_h,
                            hw_vol_at_tapping_points,
                            hw_duration,
                            no_events,
                            temp_hot_water,
                            t_it,
                        );
                    (
                        hw_energy_output,
                        pw_losses_internal,
                        pw_losses_external,
                        gains_internal_dhw_use,
                    )
                } else {
                    let hw_energy_output = self.hot_water_sources["hw cylinder"]
                        .demand_hot_water(hw_demand_vol_target, t_it)?;

                    let (pw_losses_internal, pw_losses_external, gains_internal_dhw_use) = self
                        .pipework_losses_and_internal_gains_from_hw(
                            delta_t_h,
                            hw_vol_at_tapping_points,
                            hw_duration,
                            no_events,
                            temp_hot_water,
                            t_it,
                        );

                    (
                        hw_energy_output,
                        pw_losses_internal,
                        pw_losses_external,
                        gains_internal_dhw_use,
                    )
                };

            // Convert from litres to kWh
            let cold_water_source = self.hot_water_sources["hw cylinder"].get_cold_water_source();
            let cold_water_temperature = cold_water_source.temperature(t_it, None);
            let hw_energy_demand_incl_pipework_loss =
                water_demand_to_kwh(hw_demand_vol, temp_hot_water, cold_water_temperature);
            let mut gains_internal_dhw = (pw_losses_internal + gains_internal_dhw_use)
                * WATTS_PER_KILOWATT as f64
                / t_it.timestep;
            match self.hot_water_sources.get("hw cylinder").unwrap() {
                HotWaterSource::PreHeated(ref source) => match source {
                    HotWaterStorageTank::StorageTank(storage_tank) => {
                        gains_internal_dhw += storage_tank.read().internal_gains();
                    }
                    HotWaterStorageTank::SmartHotWaterTank(smart_hot_water_tank) => {
                        gains_internal_dhw += smart_hot_water_tank.read().internal_gains();
                    }
                },
                HotWaterSource::CombiBoiler(ref source) => {
                    gains_internal_dhw += source.internal_gains();
                }
                _ => {}
            }

            // loop through on-site energy generation
            for pv in self.on_site_generation.values() {
                // Get energy produced for the current timestep
                let (_energy_produced, energy_lost) = pv.produce_energy(t_it)?;
                // Add the energy lost figure to the internal gains if it is considered inside the building
                if pv.inverter_is_inside() {
                    gains_internal_dhw += energy_lost * WATTS_PER_KILOWATT as f64 / delta_t_h;
                }
            }

            // Addition of primary_pipework_losses_kWh for reporting as part of investigation of (upstream BRE) issue #31225: FDEV A082
            let (primary_pw_losses, storage_losses) =
                if let HotWaterSource::PreHeated(source) = &self.hot_water_sources["hw cylinder"] {
                    match source {
                        HotWaterStorageTank::StorageTank(_storage_tank) => {
                            unimplemented!(".to_report() no longer exists on storage tank");
                            //storage_tank.read().to_report()
                        }
                        HotWaterStorageTank::SmartHotWaterTank(_smart_hot_water_tank) => {
                            unimplemented!(".to_report() no longer exists on storage tank");
                            //smart_hot_water_tank.read().to_report()
                        }
                    }
                } else {
                    (0.0, 0.0)
                };

            let SpaceHeatingCalculation {
                gains_internal_zone,
                gains_solar_zone,
                operative_temp,
                internal_air_temp,
                space_heat_demand_zone,
                space_cool_demand_zone,
                space_heat_provided_system: space_heat_provided,
                space_cool_provided_system: space_cool_provided,
                internal_gains_ductwork: ductwork_gains,
                heat_balance_map: heat_balance_dict,
            } = self.calc_space_heating(
                t_it.timestep,
                gains_internal_dhw,
                &mut internal_pressure_window,
                t_it,
            )?;

            // Perform calculations that can only be done after all heating
            // services have been calculated
            for system in self.timestep_end_calcs.read().iter() {
                system.timestep_end(t_it)?;
            }

            for (z_name, gains_internal) in gains_internal_zone {
                gains_internal_dict
                    .get_mut(z_name.as_str())
                    .unwrap()
                    .push(gains_internal);
            }

            for (z_name, gains_solar) in gains_solar_zone {
                gains_solar_dict
                    .get_mut(z_name.as_str())
                    .unwrap()
                    .push(gains_solar);
            }

            for (z_name, temp) in operative_temp {
                operative_temp_dict
                    .get_mut(z_name.as_str())
                    .unwrap()
                    .push(temp);
            }

            for (z_name, temp) in internal_air_temp {
                internal_air_temp_dict
                    .get_mut(z_name.as_str())
                    .unwrap()
                    .push(temp);
            }

            for (z_name, demand) in space_heat_demand_zone {
                space_heat_demand_dict
                    .get_mut(z_name.as_str())
                    .unwrap()
                    .push(demand);
            }

            for (z_name, demand) in space_cool_demand_zone {
                space_cool_demand_dict
                    .get_mut(z_name.as_str())
                    .unwrap()
                    .push(demand);
            }

            for (h_name, output) in space_heat_provided {
                space_heat_provided_dict
                    .get_mut(&h_name)
                    .unwrap()
                    .push(output);
            }

            for (c_name, output) in space_cool_provided {
                space_cool_provided_dict
                    .get_mut(&c_name)
                    .unwrap()
                    .push(output);
            }

            for (z_name, hb_dict) in heat_balance_dict {
                if let Some(hb_dict) = hb_dict {
                    for (hb_name, gains_losses) in hb_dict.as_index_map() {
                        for (heat_gains_losses_name, heat_gains_losses_value) in gains_losses {
                            heat_balance_all_dict
                                .get_mut(&hb_name)
                                .unwrap()
                                .get_mut(&z_name)
                                .unwrap()
                                .entry(heat_gains_losses_name)
                                .or_default()
                                .push(heat_gains_losses_value);
                        }
                    }
                }
            }

            hot_water_demand_dict
                .get_mut("demand")
                .unwrap()
                .push(hw_demand_vol);
            hot_water_energy_demand_dict
                .get_mut("energy_demand")
                .unwrap()
                .push(hw_energy_demand);
            hot_water_energy_demand_dict_incl_pipework
                .get_mut("energy_demand_incl_pipework_loss")
                .unwrap()
                .push(hw_energy_demand_incl_pipework_loss);
            hot_water_energy_output_dict
                .get_mut("energy_output")
                .unwrap()
                .push(hw_energy_output);
            hot_water_duration_dict
                .get_mut("duration")
                .unwrap()
                .push(hw_duration);
            hot_water_no_events_dict
                .get_mut("no_events")
                .unwrap()
                .push(no_events);
            hot_water_pipework_dict
                .get_mut("pw_losses")
                .unwrap()
                .push(pw_losses_internal + pw_losses_external);
            ductwork_gains_dict
                .get_mut("ductwork_gains")
                .unwrap()
                .push(ductwork_gains);
            hot_water_primary_pipework_dict
                .get_mut("primary_pw_losses")
                .unwrap()
                .push(primary_pw_losses);
            hot_water_storage_losses_dict
                .get_mut("storage_losses")
                .unwrap()
                .push(storage_losses);

            for supply in self.energy_supplies.values() {
                anyhow::Ok(supply.read().calc_energy_import_export_betafactor(t_it)?)?;
                anyhow::Ok(
                    supply
                        .read()
                        .calc_energy_import_from_grid_to_battery(t_it)?,
                )?;
            }

            for diverter in &self.diverters {
                diverter.write().timestep_end();
            }

            for (_, control) in &self.smart_appliance_controls {
                control.update_demand_buffer(t_it);
            }
        }

        // Return results from all energy supplies
        let mut results_totals: IndexMap<String, Vec<f64>> = Default::default();
        let mut results_end_user: IndexMap<String, IndexMap<String, Vec<f64>>> = Default::default();
        let mut energy_import: IndexMap<String, Vec<f64>> = Default::default();
        let mut energy_export: IndexMap<String, Vec<f64>> = Default::default();
        let mut energy_generated_consumed: IndexMap<String, Vec<f64>> = Default::default();
        let mut energy_to_storage: IndexMap<String, Vec<f64>> = Default::default();
        let mut energy_from_storage: IndexMap<String, Vec<f64>> = Default::default();
        let mut storage_from_grid: IndexMap<String, Vec<f64>> = Default::default();
        let mut battery_state_of_charge: IndexMap<String, Vec<f64>> = Default::default();
        let mut energy_diverted: IndexMap<String, Vec<f64>> = Default::default();
        let mut betafactor: IndexMap<String, Vec<f64>> = Default::default();
        for (name, supply) in self
            .energy_supplies
            .iter()
            .map(|(name, supply)| (name.to_owned(), Arc::clone(supply)))
        {
            let supply = supply.read();
            results_totals.insert(name.clone(), supply.results_total());
            results_end_user.insert(name.clone(), supply.results_by_end_user().to_owned());
            energy_import.insert(name.clone(), supply.get_energy_import().to_owned());
            energy_export.insert(name.clone(), supply.get_energy_export().to_owned());
            energy_generated_consumed.insert(
                name.clone(),
                supply.get_energy_generated_consumed().to_owned(),
            );
            let (energy_to, energy_from, storage_from, state_of_charge) =
                supply.get_energy_to_from_battery();
            energy_to_storage.insert(name.clone(), energy_to.to_owned());
            energy_from_storage.insert(name.clone(), energy_from.to_owned());
            storage_from_grid.insert(name.clone(), storage_from.to_owned());
            battery_state_of_charge.insert(name.clone(), state_of_charge.to_owned());
            energy_diverted.insert(name.clone(), supply.get_energy_diverted().to_owned());
            betafactor.insert(name, supply.get_beta_factor().to_owned());
        }

        let hot_water_energy_out: IndexMap<String, Vec<f64>> = IndexMap::from([(
            "hw cylinder".into(),
            hot_water_energy_output_dict
                .get("energy_output")
                .unwrap()
                .to_owned(),
        )]);
        let dhw_cop_dict = self.heat_cool_cop(
            &hot_water_energy_out,
            &results_end_user,
            self.energy_supply_conn_names_for_hot_water_source.clone(),
        );
        let heat_cop_dict = self.heat_cool_cop(
            &space_heat_provided_dict,
            &results_end_user,
            self.energy_supply_conn_names_for_heat_systems
                .iter()
                .map(|(k, v)| (k.clone(), vec![v.clone()]))
                .collect(),
        );
        let cool_cop_dict = self.heat_cool_cop(
            &space_cool_provided_dict,
            &results_end_user,
            self.space_cool_systems
                .keys()
                .map(|system_name| (system_name.clone(), vec![system_name.clone()]))
                .collect(),
        );

        let zone_dict = IndexMap::from([
            (ZoneResultKey::InternalGains, gains_internal_dict),
            (ZoneResultKey::SolarGains, gains_solar_dict),
            (ZoneResultKey::OperativeTemp, operative_temp_dict),
            (ZoneResultKey::InternalAirTemp, internal_air_temp_dict),
            (ZoneResultKey::SpaceHeatDemand, space_heat_demand_dict),
            (ZoneResultKey::SpaceCoolDemand, space_cool_demand_dict),
        ]);

        let hc_system_dict = IndexMap::from([
            (
                HeatingCoolingSystemResultKey::HeatingSystemOutput,
                space_heat_provided_dict,
            ),
            (
                HeatingCoolingSystemResultKey::CoolingSystemOutput,
                space_cool_provided_dict,
            ),
        ]);

        let hot_water_dict = IndexMap::from([
            (
                HotWaterResultKey::HotWaterDemand,
                HotWaterResultMap::Float(hot_water_demand_dict),
            ),
            (
                HotWaterResultKey::HotWaterEnergyDemandIncludingPipeworkLoss,
                HotWaterResultMap::Float(hot_water_energy_demand_dict_incl_pipework),
            ),
            (
                HotWaterResultKey::HotWaterEnergyDemand,
                HotWaterResultMap::Float(hot_water_energy_demand_dict),
            ),
            (
                HotWaterResultKey::HotWaterDuration,
                HotWaterResultMap::Float(hot_water_duration_dict),
            ),
            (
                HotWaterResultKey::HotWaterEvents,
                HotWaterResultMap::Int(hot_water_no_events_dict),
            ),
            (
                HotWaterResultKey::PipeworkLosses,
                HotWaterResultMap::Float(hot_water_pipework_dict),
            ),
            (
                HotWaterResultKey::PrimaryPipeworkLosses,
                HotWaterResultMap::Float(hot_water_primary_pipework_dict),
            ),
            (
                HotWaterResultKey::StorageLosses,
                HotWaterResultMap::Float(hot_water_storage_losses_dict),
            ),
        ]);

        // Report detailed outputs from heat source wet objects, if requested and available
        // TODO (from Python) Note that the below assumes that there is only one water
        //      heating service and therefore that all hot water energy
        //      output is assigned to that service. If the model changes in
        //      future to allow more than one hot water system, this code may
        //      need to be revised to handle that scenario.
        if self.detailed_output_heating_cooling {
            for (name, heat_source_wet) in self.wet_heat_sources.iter() {
                if let Some((results, results_annual)) = heat_source_wet
                    .output_detailed_results(&hot_water_energy_output_dict["energy_output"])
                {
                    heat_source_wet_results_dict.insert(name.clone(), results);
                    heat_source_wet_results_annual_dict.insert(name.clone(), results_annual);
                }
            }

            // Emitter detailed output results are stored with respect to heat_system_name
            for (heat_system_name, heat_system) in self.space_heat_systems.iter() {
                if let Some(emitters_output) = heat_system.lock().output_emitter_results() {
                    emitters_output_dict.insert(heat_system_name.to_owned(), emitters_output);
                }
            }

            // ESH detailed output results are stored with respect to heat_system_name
            for (heat_system_name, heat_system) in self.space_heat_systems.iter() {
                if let Some(esh_output) = heat_system.lock().output_esh_results() {
                    esh_output_dict.insert(heat_system_name.to_owned(), esh_output);
                }
            }

            vent_output_list = self.ventilation.output_vent_results().read().clone();

            // Detailed output results collected from storage tank class function
            for (name, hot_water_source) in self.hot_water_sources.iter() {
                if let HotWaterSource::PreHeated(source) = hot_water_source {
                    match source {
                        HotWaterStorageTank::StorageTank(storage_tank) => {
                            if let Some(hot_water_source_output) =
                                storage_tank.read().output_results()
                            {
                                hot_water_source_results_dict
                                    .insert(name.to_owned(), hot_water_source_output);
                            }
                        }
                        HotWaterStorageTank::SmartHotWaterTank(smart_storage_tank) => {
                            if let Some(hot_water_source_output) =
                                smart_storage_tank.read().output_results()
                            {
                                hot_water_source_results_dict
                                    .insert(name.to_owned(), hot_water_source_output);
                            }
                        }
                    }
                }
            }
        }

        Ok(RunResults {
            timestep_array,
            results_totals,
            results_end_user,
            energy_import,
            energy_export,
            energy_generated_consumed,
            energy_to_storage,
            energy_from_storage,
            storage_from_grid,
            battery_state_of_charge,
            energy_diverted,
            betafactor,
            zone_dict,
            zone_list,
            hc_system_dict,
            hot_water_dict,
            heat_cop_dict,
            cool_cop_dict,
            dhw_cop_dict,
            ductwork_gains: ductwork_gains_dict,
            heat_balance_dict: heat_balance_all_dict,
            heat_source_wet_results_dict,
            heat_source_wet_results_annual_dict,
            emitters_output_dict,
            esh_output_dict,
            vent_output_list,
            hot_water_source_results_dict,
        })
    }

    /// Calculate overall CoP over calculation period for each heating and cooling system
    fn heat_cool_cop(
        &self,
        energy_provided: &IndexMap<String, Vec<f64>>,
        results_end_user: &IndexMap<String, IndexMap<String, Vec<f64>>>,
        energy_supply_conn_name_for_space_hc_system: IndexMap<String, Vec<String>>,
    ) -> IndexMap<String, NumberOrDivisionByZero> {
        let mut hc_output_overall: IndexMap<String, f64> = Default::default();
        let mut hc_input_overall: IndexMap<String, f64> = Default::default();
        let mut cop_dict: IndexMap<String, NumberOrDivisionByZero> = Default::default();
        for (hc_name, hc_output) in energy_provided {
            hc_output_overall.insert(hc_name.to_owned(), hc_output.iter().sum::<f64>().abs());
            hc_input_overall.insert(hc_name.to_owned(), 0.);
            let energy_supply_conn_names =
                match energy_supply_conn_name_for_space_hc_system.get(hc_name.as_str()) {
                    Some(hc_name) => hc_name.clone(),
                    None => vec![],
                };
            for (fuel_name, fuel_summary) in results_end_user {
                if fuel_name == UNMET_DEMAND_SUPPLY_NAME
                    || fuel_name == ENERGY_FROM_ENVIRONMENT_SUPPLY_NAME
                {
                    continue;
                }
                for (conn_name, energy_cons) in fuel_summary {
                    if energy_supply_conn_names.contains(conn_name) {
                        *hc_input_overall.get_mut(hc_name).unwrap() +=
                            energy_cons.iter().sum::<f64>();
                    }
                }
            }

            cop_dict.insert(
                hc_name.to_owned(),
                if hc_input_overall[hc_name] > 0. {
                    NumberOrDivisionByZero::Number(
                        hc_output_overall[hc_name] / hc_input_overall[hc_name],
                    )
                } else {
                    NumberOrDivisionByZero::DivisionByZero
                },
            );
        }

        cop_dict
    }
}

struct HeatCoolOutputs {
    hc_output_convective: IndexMap<String, f64>,
    hc_output_radiative: IndexMap<String, f64>,
    hc_output_min: IndexMap<String, f64>,
}

// let (mut hc_output_convective, mut hc_output_radiative, hc_output_min) = self
//                 .heat_cool_system_output_min(
//                     &h_name_list_sorted_zone,
//                     &c_name_list_sorted_zone,
//                     &frac_convective_heat_zone_system,
//                     &frac_convective_cool_zone_system,
//                     z_name,
//                     simtime,
//                 )?;

#[derive(Debug, Default)]
pub struct OutputOptions {
    pub print_heat_balance: bool,
    pub detailed_output_heating_cooling: bool,
}

struct SetpointsAndConvectiveFractions {
    temp_setpnt_heat: IndexMap<String, f64>,
    temp_setpnt_cool: IndexMap<String, f64>,
    frac_convective_heat: IndexMap<String, f64>,
    frac_convective_cool: IndexMap<String, f64>,
}

#[derive(Clone, Debug)]
pub enum HotWaterResultMap {
    Float(IndexMap<String, Vec<f64>>),
    Int(IndexMap<String, Vec<usize>>),
}

#[derive(Clone, Copy, Debug)]
pub enum NumberOrDivisionByZero {
    Number(f64),
    DivisionByZero,
}

impl Display for NumberOrDivisionByZero {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                NumberOrDivisionByZero::Number(number) => Cow::Owned(format!("{}", number)),
                NumberOrDivisionByZero::DivisionByZero => Cow::Borrowed("DIV/0"),
            }
        )
    }
}

pub type ResultsEndUser = IndexMap<String, IndexMap<String, Vec<f64>>>;

enum SpaceHeatCoolSystems<'a> {
    Heat(&'a IndexMap<String, Arc<Mutex<SpaceHeatSystem>>>),
    Cool(&'a IndexMap<String, AirConditioning>),
}

impl SpaceHeatCoolSystems<'_> {
    fn in_required_period_for_name(
        &self,
        system_name: &str,
        simtime: SimulationTimeIteration,
    ) -> Option<bool> {
        match self {
            SpaceHeatCoolSystems::Heat(heat) => {
                heat[system_name].lock().in_required_period(simtime)
            }
            SpaceHeatCoolSystems::Cool(cool) => cool[system_name].in_required_period(&simtime),
        }
    }
}

pub type ColdWaterSources = IndexMap<String, Arc<ColdWaterSource>>;

fn cold_water_sources_from_input(input: &ColdWaterSourceInput) -> ColdWaterSources {
    input
        .iter()
        .map(|(source_type, source_details)| {
            (
                source_type.into(),
                Arc::from(cold_water_source_from_input_details(source_details)),
            )
        })
        .collect()
}

fn cold_water_source_from_input_details(details: &ColdWaterSourceDetails) -> ColdWaterSource {
    ColdWaterSource::new(
        details.temperatures.clone(),
        details.start_day,
        details.time_series_step,
    )
}

fn energy_supplies_from_input(
    input: &EnergySupplyInput,
    simulation_time_iterator: &SimulationTimeIterator,
    tariff_data_file: Option<&str>,
    external_conditions: Arc<ExternalConditions>,
) -> anyhow::Result<IndexMap<String, Arc<RwLock<EnergySupply>>>> {
    let mut supplies: IndexMap<String, Arc<RwLock<EnergySupply>>> = IndexMap::new();

    // set up supply representing unmet demand
    supplies.insert(
        UNMET_DEMAND_SUPPLY_NAME.into(),
        Arc::new(RwLock::new(
            EnergySupplyBuilder::new(
                FuelType::UnmetDemand,
                simulation_time_iterator.total_steps(),
            )
            .build(),
        )),
    );
    supplies.insert(
        ENERGY_FROM_ENVIRONMENT_SUPPLY_NAME.into(),
        Arc::new(RwLock::new(
            EnergySupplyBuilder::new(
                FuelType::EnergyFromEnvironment,
                simulation_time_iterator.total_steps(),
            )
            .build(),
        )),
    );
    for (name, supply) in input {
        let energy_supply = energy_supply_from_input(
            supply,
            simulation_time_iterator,
            tariff_data_file,
            external_conditions.clone(),
        )?;
        supplies.insert(name.into(), energy_supply);
    }

    Ok(supplies)
}

fn energy_supply_from_input(
    input: &EnergySupplyDetails,
    simulation_time_iterator: &SimulationTimeIterator,
    tariff_file_path: Option<&str>,
    external_conditions: Arc<ExternalConditions>,
) -> anyhow::Result<Arc<RwLock<EnergySupply>>> {
    Ok(Arc::new(RwLock::new({
        let mut builder =
            EnergySupplyBuilder::new(input.fuel, simulation_time_iterator.total_steps());
        if let Some(battery) = input.electric_battery.as_ref() {
            builder = builder.with_electric_battery(ElectricBattery::from_input(
                battery,
                simulation_time_iterator.step_in_hours(),
                external_conditions,
            ))
        }
        if let Some(priority) = input.priority.as_ref() {
            builder = builder.with_priority(priority.clone());
        }
        builder = builder.with_export_capable(input.is_export_capable);

        if input
            .electric_battery
            .as_ref()
            .is_some_and(|battery| battery.grid_charging_possible)
        {
            let tariff_data: Box<dyn Read> = match tariff_file_path {
                // fall back to using tariff data for entire year for now
                None => Box::new(Cursor::new(include_str!(
                    "../examples/tariff_data/tariff_data_25-06-2024.csv"
                ))),
                Some(tariff_file_path) => Box::new(BufReader::new(
                    File::open(tariff_file_path)
                        .expect("Provided tariff file at provided path was not found."),
                )),
            };
            builder = builder.with_tariff_input(EnergySupplyTariffInput::new(
                input.tariff.ok_or_else(|| anyhow!("Energy supply with electric battery that allows grid charging expected tariff to be indicated"))?,
                tariff_data,
                input.threshold_charges.map(|threshold_charges| threshold_charges.to_vec()),
                input.threshold_prices.map(|threshold_prices| threshold_prices.to_vec()),
            ))?;
        }

        builder.build()
    })))
}

type DiverterTypes = IndexMap<String, EnergyDiverter>;

// struct DiverterTypes {
//     pub mains_electricity: Option<EnergyDiverter>,
//     pub mains_gas: Option<EnergyDiverter>,
//     pub bulk_lpg: Option<EnergyDiverter>,
// }

// impl From<&EnergySupplyInput> for DiverterTypes {
//     fn from(input: &EnergySupplyInput) -> Self {
//         Self {
//             mains_electricity: diverter_from_energy_supply(
//                 input.get(&EnergySupplyKey::MainsElectricity),
//             ),
//             mains_gas: diverter_from_energy_supply(input.get(&EnergySupplyKey::MainsGas)),
//             bulk_lpg: diverter_from_energy_supply(input.get(&EnergySupplyKey::BulkLpg)),
//         }
//     }
// }
//
// impl DiverterTypes {
//     pub fn get_for_supply_type(
//         &self,
//         energy_supply_type: EnergySupplyType,
//     ) -> Option<&EnergyDiverter> {
//         match energy_supply_type {
//             EnergySupplyType::Electricity => self.mains_electricity.as_ref(),
//             EnergySupplyType::MainsGas => self.mains_gas.as_ref(),
//             EnergySupplyType::LpgBulk => self.bulk_lpg.as_ref(),
//             _ => None,
//         }
//     }
// }

fn diverter_from_energy_supply(supply: &EnergySupplyDetails) -> Option<EnergyDiverter> {
    supply.diverter.clone()
}

pub(crate) type InternalGainsCollection = IndexMap<Arc<str>, Gains>;

fn internal_gains_from_input(
    input: &InternalGainsInput,
    total_floor_area: f64,
) -> anyhow::Result<InternalGainsCollection> {
    let mut gains_collection = InternalGainsCollection::from([]);
    if let Some(internal_gains) = input.total_internal_gains.as_ref() {
        gains_collection.insert(
            "total_internal_gains".into(),
            Gains::Internal(internal_gains_from_details(
                internal_gains,
                total_floor_area,
            )?),
        );
    }
    if let Some(internal_gains) = input.metabolic_gains.as_ref() {
        gains_collection.insert(
            "metabolic_gains".into(),
            Gains::Internal(internal_gains_from_details(
                internal_gains,
                total_floor_area,
            )?),
        );
    }
    if let Some(internal_gains) = input.evaporative_losses.as_ref() {
        gains_collection.insert(
            "evaporative_losses".into(),
            Gains::Internal(internal_gains_from_details(
                internal_gains,
                total_floor_area,
            )?),
        );
    }
    if let Some(internal_gains) = input.cold_water_losses.as_ref() {
        gains_collection.insert(
            "coldwaterlosses".into(),
            Gains::Internal(internal_gains_from_details(
                internal_gains,
                total_floor_area,
            )?),
        );
    }
    if let Some(internal_gains) = input.other.as_ref() {
        gains_collection.insert(
            "other".into(),
            Gains::Internal(internal_gains_from_details(
                internal_gains,
                total_floor_area,
            )?),
        );
    }

    Ok(gains_collection)
}

fn internal_gains_from_details(
    details: &InternalGainsDetails,
    total_floor_area: f64,
) -> anyhow::Result<InternalGains> {
    Ok(InternalGains::new(
        convert_energy_to_wm2(details, total_floor_area)?,
        details.start_day,
        details.time_series_step,
    ))
}

fn convert_energy_to_wm2(
    internal_gains_details: &InternalGainsDetails,
    total_floor_area: f64,
) -> anyhow::Result<Vec<f64>> {
    let schedule = &internal_gains_details.schedule;
    Ok(reject_nulls(expand_numeric_schedule(schedule))?
        .iter()
        .map(|energy_data| energy_data / total_floor_area)
        .collect())
}

#[derive(Debug)]
pub struct Controls {
    core: Vec<HeatSourceControl>,
    extra: HashMap<String, Arc<Control>>,
}

impl Controls {
    pub(crate) fn new(core: Vec<HeatSourceControl>, extra: HashMap<String, Arc<Control>>) -> Self {
        Self { core, extra }
    }

    pub(crate) fn get(&self, control_type: &HeatSourceControlType) -> Option<Arc<Control>> {
        self.core
            .iter()
            .find(|heat_source_control| heat_source_control.has_type(*control_type))
            .map(|heat_source_control| heat_source_control.get())
    }

    pub(crate) fn get_with_string(&self, control_name: &str) -> Option<Arc<Control>> {
        match control_name {
            // hard-code ways of resolving to core control types (for now)
            "hw timer" => self.get(&HeatSourceControlType::HotWaterTimer),
            "window opening" => self.get(&HeatSourceControlType::WindowOpening),
            other => self.extra.get(other).cloned(),
        }
    }
}

fn wwhrs_from_input(
    wwhrs: Option<&WasteWaterHeatRecovery>,
    cold_water_sources: &ColdWaterSources,
    initial_simtime: SimulationTimeIteration,
) -> anyhow::Result<IndexMap<String, Arc<Mutex<Wwhrs>>>> {
    let mut wwhr_systems: IndexMap<String, Arc<Mutex<Wwhrs>>> = IndexMap::from([]);
    if let Some(systems) = wwhrs {
        for (name, system) in systems {
            wwhr_systems
                .entry(name.into())
                .or_insert(Arc::new(Mutex::new(wwhr_system_from_details(
                    system.clone(),
                    cold_water_sources,
                    initial_simtime,
                )?)));
        }
    }

    Ok(wwhr_systems)
}

fn wwhr_system_from_details(
    _system: WasteWaterHeatRecoveryDetails,
    _cold_water_sources: &ColdWaterSources,
    _initial_simtime: SimulationTimeIteration,
) -> anyhow::Result<Wwhrs> {
    todo!();
    // Ok(match system.system_type {
    //     WasteWaterHeatRecoverySystemType::SystemA => {
    //         Wwhrs::WWHRSInstantaneousSystemA(WWHRSInstantaneousSystemA::new(
    //             system.flow_rates,
    //             system.efficiencies,
    //             get_cold_water_source_ref_for_type(system.cold_water_source, cold_water_sources)
    //                 .ok_or_else(|| {
    //                     anyhow!(
    //                         "Could not find cold water source '{:?}'",
    //                         system.cold_water_source
    //                     )
    //                 })?,
    //             system.utilisation_factor,
    //             initial_simtime,
    //         ))
    //     }
    //     WasteWaterHeatRecoverySystemType::SystemB => {
    //         Wwhrs::WWHRSInstantaneousSystemB(WWHRSInstantaneousSystemB::new(
    //             get_cold_water_source_ref_for_type(system.cold_water_source, cold_water_sources)
    //                 .ok_or_else(|| {
    //                     anyhow!(
    //                         "Could not find cold water source '{:?}'",
    //                         system.cold_water_source
    //                     )
    //                 })?,
    //             system.flow_rates,
    //             system.efficiencies,
    //             system.utilisation_factor,
    //         ))
    //     }
    //     WasteWaterHeatRecoverySystemType::SystemC => {
    //         Wwhrs::WWHRSInstantaneousSystemC(WWHRSInstantaneousSystemC::new(
    //             system.flow_rates,
    //             system.efficiencies,
    //             get_cold_water_source_ref_for_type(system.cold_water_source, cold_water_sources)
    //                 .ok_or_else(|| {
    //                     anyhow!(
    //                         "Could not find cold water source '{:?}'",
    //                         system.cold_water_source
    //                     )
    //                 })?,
    //             system.utilisation_factor,
    //             initial_simtime,
    //         ))
    //     }
    // })
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) enum ReportingFlag {
    Min,
    Max,
}

impl Display for ReportingFlag {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                ReportingFlag::Min => "min",
                ReportingFlag::Max => "max",
            }
        )
    }
}

fn temp_internal_air_for_zones(zones: &IndexMap<String, Arc<Zone>>, total_volume: f64) -> f64 {
    let internal_air_temperature = zones
        .values()
        .map(|zone| zone.temp_internal_air() * zone.volume())
        .sum::<f64>();

    internal_air_temperature / total_volume
}

pub(crate) type TempInternalAirFn = Arc<dyn Fn() -> f64 + Send + Sync>;

fn shareable_fn(num: &Arc<RwLock<f64>>) -> TempInternalAirFn {
    let clone = num.clone();
    Arc::from(move || *clone.read())
}

/// A struct definition to encapsulate results from a corpus run.
#[derive(Debug)]
pub struct RunResults {
    pub timestep_array: Vec<f64>,
    pub results_totals: IndexMap<String, Vec<f64>>,
    pub results_end_user: ResultsEndUser,
    pub energy_import: IndexMap<String, Vec<f64>>,
    pub energy_export: IndexMap<String, Vec<f64>>,
    pub(crate) energy_generated_consumed: IndexMap<String, Vec<f64>>,
    pub energy_to_storage: IndexMap<String, Vec<f64>>,
    pub energy_from_storage: IndexMap<String, Vec<f64>>,
    pub storage_from_grid: IndexMap<String, Vec<f64>>,
    pub battery_state_of_charge: IndexMap<String, Vec<f64>>,
    pub(crate) energy_diverted: IndexMap<String, Vec<f64>>,
    pub(crate) betafactor: IndexMap<String, Vec<f64>>,
    pub(crate) zone_dict: IndexMap<ZoneResultKey, IndexMap<String, Vec<f64>>>,
    pub(crate) zone_list: Vec<String>,
    pub(crate) hc_system_dict: IndexMap<HeatingCoolingSystemResultKey, IndexMap<String, Vec<f64>>>,
    pub(crate) hot_water_dict: IndexMap<HotWaterResultKey, HotWaterResultMap>,
    pub(crate) heat_cop_dict: IndexMap<String, NumberOrDivisionByZero>,
    pub(crate) cool_cop_dict: IndexMap<String, NumberOrDivisionByZero>,
    pub(crate) dhw_cop_dict: IndexMap<String, NumberOrDivisionByZero>,
    pub(crate) ductwork_gains: IndexMap<String, Vec<f64>>,
    pub(crate) heat_balance_dict: HeatBalanceAllResults,
    pub(crate) heat_source_wet_results_dict: IndexMap<String, ResultsPerTimestep>,
    pub(crate) heat_source_wet_results_annual_dict: IndexMap<String, ResultsAnnual>,
    pub(crate) emitters_output_dict: IndexMap<String, Vec<EmittersDetailedResult>>,
    pub(crate) esh_output_dict: IndexMap<String, Vec<StorageHeaterDetailedResult>>,
    pub(crate) vent_output_list: Vec<VentilationDetailedResult>,
    pub(crate) hot_water_source_results_dict: IndexMap<String, Vec<StorageTankDetailedResult>>,
}

impl RunResults {
    pub fn space_heat_demand_total(&self) -> f64 {
        self.zone_dict[&ZoneResultKey::SpaceHeatDemand]
            .values()
            .map(|v| v.iter().sum::<f64>())
            .sum::<f64>()
    }

    pub fn space_cool_demand_total(&self) -> f64 {
        self.zone_dict[&ZoneResultKey::SpaceCoolDemand]
            .values()
            .map(|v| v.iter().sum::<f64>())
            .sum::<f64>()
    }
}

pub(crate) type HeatBalanceAllResults =
    IndexMap<HeatBalanceFieldName, IndexMap<String, IndexMap<String, Vec<f64>>>>;

struct SpaceHeatingCalculation {
    gains_internal_zone: HashMap<String, f64>,
    gains_solar_zone: HashMap<String, f64>,
    operative_temp: HashMap<String, f64>,
    internal_air_temp: HashMap<String, f64>,
    space_heat_demand_zone: HashMap<String, f64>,
    space_cool_demand_zone: HashMap<String, f64>,
    space_heat_provided_system: HashMap<String, f64>,
    space_cool_provided_system: HashMap<String, f64>,
    internal_gains_ductwork: f64,
    heat_balance_map: HashMap<String, Option<HeatBalance>>,
}

#[derive(Clone, Copy, Deserialize_enum_str, Debug, Eq, Hash, PartialEq, Serialize_enum_str)]
pub(crate) enum ZoneResultKey {
    #[serde(rename = "internal gains")]
    InternalGains,
    #[serde(rename = "solar gains")]
    SolarGains,
    #[serde(rename = "operative temp")]
    OperativeTemp,
    #[serde(rename = "internal air temp")]
    InternalAirTemp,
    #[serde(rename = "space heat demand")]
    SpaceHeatDemand,
    #[serde(rename = "space cool demand")]
    SpaceCoolDemand,
}

impl ZoneResultKey {
    pub(crate) fn as_str(&self) -> &'static str {
        // if there's a clean way of reusing the serde deserialization map above, this would be preferable
        // but for now, this minor duplication will do
        match self {
            ZoneResultKey::InternalGains => "internal gains",
            ZoneResultKey::SolarGains => "solar gains",
            ZoneResultKey::OperativeTemp => "operative temp",
            ZoneResultKey::InternalAirTemp => "internal air temp",
            ZoneResultKey::SpaceHeatDemand => "space heat demand",
            ZoneResultKey::SpaceCoolDemand => "space cool demand",
        }
    }
}

#[derive(Clone, Copy, Deserialize_enum_str, Debug, Eq, Hash, PartialEq, Serialize_enum_str)]
pub(crate) enum HeatingCoolingSystemResultKey {
    #[serde(rename = "Heating system")]
    HeatingSystem,
    #[serde(rename = "Heating system output")]
    HeatingSystemOutput,
    #[serde(rename = "Cooling system")]
    CoolingSystem,
    #[serde(rename = "Cooling system output")]
    CoolingSystemOutput,
}

#[derive(Clone, Copy, Deserialize_enum_str, Debug, Eq, Hash, PartialEq, Serialize_enum_str)]
pub(crate) enum HotWaterResultKey {
    #[serde(rename = "Hot water demand")]
    HotWaterDemand,
    #[serde(rename = "Hot water energy demand incl pipework_loss")]
    HotWaterEnergyDemandIncludingPipeworkLoss,
    #[serde(rename = "Hot water energy demand")]
    HotWaterEnergyDemand,
    #[serde(rename = "Hot water duration")]
    HotWaterDuration,
    #[serde(rename = "Hot Water Events")]
    HotWaterEvents,
    #[serde(rename = "Pipework losses")]
    PipeworkLosses,
    #[serde(rename = "Primary pipework losses")]
    PrimaryPipeworkLosses,
    #[serde(rename = "Storage losses")]
    StorageLosses,
}

fn get_cold_water_source_ref_for_type(
    source_type: &str,
    cold_water_sources: &ColdWaterSources,
) -> Option<Arc<ColdWaterSource>> {
    cold_water_sources.get(source_type).cloned()
}

pub(crate) type EventSchedule = Vec<Option<Vec<TypedScheduleEvent>>>;

fn event_schedules_from_input(
    events: &WaterHeatingEvents,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<EventSchedule> {
    let mut schedule: EventSchedule = vec![None; simulation_time_iterator.total_steps()];
    for (name, events) in &events.shower {
        schedule = schedule_event_from_input(
            events.iter().collect(),
            name.as_str(),
            WaterScheduleEventType::Shower,
            schedule,
            simulation_time_iterator,
        )?;
    }

    for (name, events) in &events.bath {
        schedule = schedule_event_from_input(
            events.iter().collect(),
            name.as_str(),
            WaterScheduleEventType::Bath,
            schedule,
            simulation_time_iterator,
        )?;
    }

    for (name, events) in &events.other {
        schedule = schedule_event_from_input(
            events.iter().collect(),
            name.as_str(),
            WaterScheduleEventType::Other,
            schedule,
            simulation_time_iterator,
        )?;
    }

    Ok(schedule)
}

fn schedule_event_from_input(
    events_input: Vec<&WaterHeatingEvent>,
    name: &str,
    event_type: WaterScheduleEventType,
    existing_schedule: EventSchedule,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<EventSchedule> {
    let sim_timestep = simulation_time_iterator.step_in_hours();
    let total_timesteps = simulation_time_iterator.total_steps();

    expand_events(
        events_input
            .iter()
            .map(|event| ScheduleEvent::from(*event))
            .collect::<Vec<_>>(),
        sim_timestep,
        total_timesteps,
        name,
        event_type,
        existing_schedule,
    )
}

fn thermal_bridging_from_input(input: &ThermalBridgingInput) -> ThermalBridging {
    match input {
        ThermalBridgingInput::Elements(input_bridges) => ThermalBridging::Bridges({
            let mut bridges = IndexMap::new();
            bridges.extend(input_bridges.iter().map(|(name, details)| {
                (
                    name.into(),
                    match details {
                        ThermalBridgingDetails::Linear {
                            linear_thermal_transmittance,
                            length,
                            ..
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
        ThermalBridgingInput::Number(num) => ThermalBridging::Number(*num),
    }
}

fn zone_from_input(
    input: &ZoneInput,
    zone_name: &str,
    heat_system_name_for_zone: &mut IndexMap<String, Vec<String>>,
    cool_system_name_for_zone: &mut IndexMap<String, Vec<String>>,
    external_conditions: Arc<ExternalConditions>,
    infiltration_ventilation: Arc<InfiltrationVentilation>,
    window_adjust_control: Option<Arc<dyn ControlBehaviour>>,
    controls: &Controls,
    print_heat_balance: bool,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<Zone> {
    let heat_system_name = input.space_heat_system.clone();
    let cool_system_name = input.space_cool_system.clone();

    let heat_system_names = match heat_system_name {
        SystemReference::None(_) => vec!["".into()], // equivalent of [None] in Python - we are using empty string to denote absence rather than using Option<String> everywhere
        SystemReference::Single(name) => vec![name.clone()],
        SystemReference::Multiple(names) => names.clone(),
    };

    for zone_h_name in heat_system_name_for_zone.values() {
        let zone_h_name_set: HashSet<String> = HashSet::from_iter(zone_h_name.iter().cloned());
        let h_overassigned: Vec<String> = HashSet::from_iter(heat_system_names.clone().into_iter())
            .intersection(&zone_h_name_set)
            .filter(|&name| !name.is_empty())
            .cloned()
            .collect_vec();
        if !h_overassigned.is_empty() {
            bail!(
                "Invalid input: SpaceHeatSystem ({}) has been assigned to more than one Zone",
                h_overassigned.into_iter().join(", ")
            )
        }
    }

    heat_system_name_for_zone.insert(zone_name.into(), heat_system_names);

    let cool_system_names = match cool_system_name {
        SystemReference::None(_) => vec!["".into()], // equivalent of [None] in Python - we are using empty string to denote absence rather than using Option<String> everywhere
        SystemReference::Single(name) => vec![name.clone()],
        SystemReference::Multiple(names) => names.clone(),
    };

    for zone_c_name in cool_system_name_for_zone.values() {
        let zone_c_name_set: HashSet<String> = HashSet::from_iter(zone_c_name.iter().cloned());
        let c_overassigned: Vec<String> = HashSet::from_iter(cool_system_names.clone().into_iter())
            .intersection(&zone_c_name_set)
            .filter(|&name| !name.is_empty())
            .cloned()
            .collect_vec();
        if !c_overassigned.is_empty() {
            bail!(
                "Invalid input: SpaceCoolSystem ({}) has been assigned to more than one Zone",
                c_overassigned.into_iter().join(", ")
            )
        }
    }

    cool_system_name_for_zone.insert(zone_name.into(), cool_system_names);

    // Default setpoint basis to 'operative' if not provided in input
    let temp_setpnt_basis = input
        .temp_setpnt_basis
        .unwrap_or(ZoneTemperatureControlBasis::Operative);

    Zone::new(
        input.area,
        input.volume,
        input
            .building_elements
            .iter()
            .map(|(element_name, el)| {
                Ok((
                    element_name.into(),
                    building_element_from_input(
                        el,
                        external_conditions.clone(),
                        controls,
                        simulation_time_iterator,
                    )?,
                ))
            })
            .collect::<anyhow::Result<IndexMap<String, Arc<BuildingElement>>>>()?,
        thermal_bridging_from_input(&input.thermal_bridging),
        infiltration_ventilation,
        external_conditions.air_temp(&simulation_time_iterator.current_iteration()),
        input.temp_setpnt_init,
        temp_setpnt_basis,
        window_adjust_control,
        print_heat_balance,
        simulation_time_iterator,
    )
}

#[allow(clippy::type_complexity)]
fn infiltration_ventilation_from_input(
    zones: &ZoneDictionary,
    input: &InfiltrationVentilationInput,
    controls: &Controls,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    detailed_output_heating_cooling: bool,
) -> anyhow::Result<(
    InfiltrationVentilation,
    Option<Arc<dyn ControlBehaviour>>,
    Option<Arc<Control>>,
    Option<Arc<Control>>,
)> {
    let window_adjust_control = input
        .control_window_adjust
        .as_ref()
        .and_then(|ctrl_name| controls.get_with_string(ctrl_name))
        .map(|ctrl| ctrl.clone() as Arc<dyn ControlBehaviour>);
    let vent_adjust_min_control = input
        .control_vent_adjust_min
        .as_ref()
        .and_then(|ctrl_name| controls.get_with_string(ctrl_name));
    let vent_adjust_max_control = input
        .control_vent_adjust_max
        .as_ref()
        .and_then(|ctrl_name| controls.get_with_string(ctrl_name));

    let ventilation = InfiltrationVentilation::create(
        input,
        zones,
        detailed_output_heating_cooling,
        energy_supplies,
        controls,
    )?;

    Ok((
        ventilation,
        window_adjust_control,
        vent_adjust_min_control,
        vent_adjust_max_control,
    ))
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct CompletedVentilationLeaks {
    pub(crate) ventilation_zone_height: f64,
    pub(crate) test_pressure: f64,
    pub(crate) test_result: f64,
    pub(crate) area_roof: f64,
    pub(crate) area_facades: f64,
    pub(crate) env_area: f64,
    pub(crate) altitude: f64,
}

impl CompletedVentilationLeaks {
    pub(crate) fn complete_input(
        input: &InfiltrationVentilationInput,
        area_facades: f64,
        area_roof: f64,
    ) -> Self {
        let VentilationLeaks {
            ventilation_zone_height,
            test_pressure,
            test_result,
            env_area,
            ..
        } = input.leaks;
        Self {
            ventilation_zone_height,
            test_pressure,
            test_result,
            area_roof,
            area_facades,
            env_area,
            altitude: input.altitude,
        }
    }
}

fn building_element_from_input(
    input: &BuildingElementInput,
    external_conditions: Arc<ExternalConditions>,
    _controls: &Controls,
    _simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<Arc<BuildingElement>> {
    Ok(Arc::from(match input {
        BuildingElementInput::Opaque {
            // is_unheated_pitched_roof,
            // area_input,
            // pitch,
            // solar_absorption_coeff,
            // u_value_input,
            // areal_heat_capacity,
            // mass_distribution_class,
            // orientation,
            // base_height,
            ..
        } => {
            todo!("BuildingElementInput::Opaque cannot yet be built from input in 0.40");
            // let is_unheated_pitched_roof = if *pitch < PITCH_LIMIT_HORIZ_CEILING {
            //     is_unheated_pitched_roof
            //         .ok_or_else(|| anyhow!("Pitch of opaque building element was {pitch} degrees, so it is necessary for this element to indicate whether this is an unheated pitched roof."))?
            // } else {
            //     false
            // };
            //
            // BuildingElement::Opaque(BuildingElementOpaque::new(
            //     *area,
            //     is_unheated_pitched_roof,
            //     *pitch,
            //     *solar_absorption_coeff,
            //     init_resistance_or_uvalue_from_input_struct(u_value_input, *pitch)?,
            //     *areal_heat_capacity,
            //     *mass_distribution_class,
            //     *orientation,
            //     *base_height,
            //     *height,
            //     *width,
            //     external_conditions,
            // ))
        }
        BuildingElementInput::Transparent {
            // u_value_input,
            // pitch,
            // orientation,
            // g_value,
            // frame_area_fraction,
            // base_height,
            // shading,
            // treatment,
            ..
        } => {
            todo!("Transparent building elements are not yet supported in this 0.40 version");
            // BuildingElement::Transparent(BuildingElementTransparent::new(
            //     *pitch,
            //     init_resistance_or_uvalue_from_input_struct(u_value_input, *pitch)?,
            //     *orientation,
            //     *g_value,
            //     *frame_area_fraction,
            //     *base_height,
            //     *height,
            //     *width,
            //     Some(shading.clone()),
            //     treatment
            //         .iter()
            //         .map(|t| {
            //             WindowTreatment::from_input(
            //                 t,
            //                 controls,
            //                 simulation_time_iterator.current_hour(),
            //             )
            //         })
            //         .collect_vec(),
            //     external_conditions,
            // ))
        }
        BuildingElementInput::Ground {
            area,
            total_area,
            pitch,
            u_value,
            thermal_resistance_floor_construction,
            areal_heat_capacity,
            mass_distribution_class,
            floor_data,
            thickness_walls,
            perimeter,
            psi_wall_floor_junc,
            ..
        } => BuildingElement::Ground(BuildingElementGround::new(
            *total_area,
            *area,
            *pitch,
            *u_value,
            *thermal_resistance_floor_construction,
            *areal_heat_capacity,
            *mass_distribution_class,
            floor_data,
            *thickness_walls,
            *perimeter,
            *psi_wall_floor_junc,
            external_conditions,
        )?),
        BuildingElementInput::AdjacentConditionedSpace {
            area,
            pitch,
            u_value_input,
            areal_heat_capacity,
            mass_distribution_class,
            ..
        } => {
            BuildingElement::AdjacentConditionedSpace(BuildingElementAdjacentConditionedSpace::new(
                area.ok_or_else(|| anyhow!("AdjacentConditionedSpace building element is expected have an area"))?,
                *pitch,
                init_resistance_or_uvalue_from_input_struct(u_value_input, *pitch)?,
                *areal_heat_capacity,
                *mass_distribution_class,
                external_conditions,
            ))
        }
        BuildingElementInput::AdjacentUnconditionedSpace {
            area,
            pitch,
            u_value_input,
            thermal_resistance_unconditioned_space,
            areal_heat_capacity,
            mass_distribution_class,
        } => BuildingElement::AdjacentUnconditionedSpaceSimple(
            BuildingElementAdjacentUnconditionedSpaceSimple::new(
                area.ok_or_else(|| anyhow!("AdjacentUnconditionedSpace building element is expected have an area"))?,
                *pitch,
                init_resistance_or_uvalue_from_input_struct(u_value_input, *pitch)?,
                *thermal_resistance_unconditioned_space,
                *areal_heat_capacity,
                *mass_distribution_class,
                external_conditions,
            ),
        ),
    }))
}

fn set_up_energy_supply_unmet_demand_zones(
    unmet_demand_supply: Arc<RwLock<EnergySupply>>,
    zones: &ZoneDictionary,
) -> IndexMap<String, Arc<EnergySupplyConnection>> {
    let mut energy_supplies: IndexMap<String, Arc<EnergySupplyConnection>> = Default::default();

    for name in zones.keys() {
        energy_supplies.insert(
            name.into(),
            Arc::new(EnergySupply::connection(unmet_demand_supply.clone(), name.as_str()).unwrap()),
        );
    }

    energy_supplies
}

fn apply_appliance_gains_from_input(
    internal_gains_collection: &mut InternalGainsCollection,
    input: &ApplianceGainsInput,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    total_floor_area: f64,
    smart_appliance_controls: &IndexMap<String, Arc<SmartApplianceControl>>,
    simulation_time: &SimulationTimeIterator,
) -> anyhow::Result<()> {
    fn check_priority(appliance_gains: &ApplianceGainsInput) -> ApplianceGainsInput {
        let (defined_priority, priorities): (Vec<std::string::String>, Vec<isize>) =
            appliance_gains
                .iter()
                .filter_map(|(app, x)| x.priority.as_ref().map(|priority| (app.clone(), priority)))
                .unzip();
        let mut lowest_priority = priorities.into_iter().max().unwrap_or(0);
        let mut new_details = appliance_gains.clone();
        for appliance in appliance_gains.keys() {
            if defined_priority.contains(appliance) {
                continue;
            }
            lowest_priority += 1;
            new_details
                .get_mut(appliance)
                .unwrap()
                .priority
                .replace(lowest_priority);
        }

        new_details
    }

    let sorted_input = check_priority(input)
        .sorted_by(|_, details1, _, details2| {
            details1.priority.expect(
                "All details in the output from check_priority were expected to have a priority.",
            ).cmp(&details2.priority.expect("All details in the output from check_priority were expected to have a priority."))
        })
        .rev()
        .collect::<IndexMap<_, _>>();

    for (name, gains_details) in sorted_input {
        let energy_supply_conn = EnergySupply::connection(
            energy_supplies
                .get(&gains_details.energy_supply)
                .ok_or_else(|| {
                    anyhow!(
                        "Appliance gains reference an undeclared energy supply '{}'.",
                        gains_details.energy_supply
                    )
                })?
                .clone(),
            &name,
        )?;

        let smart_control = gains_details
            .load_shifting
            .as_ref()
            .and_then(|load_shifting| {
                load_shifting
                    .control
                    .as_ref()
                    .and_then(|load_shifting_control| {
                        smart_appliance_controls.get(load_shifting_control)
                    })
            });

        let gains = if gains_details.events.is_some() && gains_details.standby.is_some() {
            Gains::Event(EventApplianceGains::new(
                energy_supply_conn,
                simulation_time,
                &gains_details,
                total_floor_area,
                smart_control.cloned(),
            )?)
        } else {
            Gains::Appliance(appliance_gains_from_single_input(
                &gains_details,
                energy_supply_conn,
                total_floor_area,
            )?)
        };

        internal_gains_collection.insert(name.into(), gains);
    }

    Ok(())
}

fn appliance_gains_from_single_input(
    input: &ApplianceGainsDetails,
    energy_supply_connection: EnergySupplyConnection,
    total_floor_area: f64,
) -> anyhow::Result<ApplianceGains> {
    let total_energy_supply = reject_nulls(expand_numeric_schedule(
        input
            .schedule
            .as_ref()
            .ok_or_else(|| anyhow!("Appliance gains did not have schedule when expected."))?,
    ))?
    .iter()
    .map(|energy_data| energy_data / total_floor_area)
    .collect();

    Ok(ApplianceGains::new(
        total_energy_supply,
        input.gains_fraction,
        input.start_day,
        input.time_series_step,
        energy_supply_connection,
    ))
}

#[derive(Debug)]
pub(crate) enum HeatSource {
    Storage(HeatSourceWithStorageTank),
    Wet(Box<HeatSourceWet>),
}

impl HeatSource {
    pub(crate) fn setpnt(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(Option<f64>, Option<f64>)> {
        match self {
            HeatSource::Storage(ref storage) => match storage {
                HeatSourceWithStorageTank::Immersion(imm) => Ok(imm.lock().setpnt(simtime)),
                HeatSourceWithStorageTank::Solar(ref solar) => Ok(solar.lock().setpnt(&simtime)),
            },
            HeatSource::Wet(ref heat_source) => heat_source.setpnt(simtime),
        }
    }
    pub(crate) fn _demand_energy(
        &mut self,
        energy_demand: f64,
        temp_return: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        match self {
            HeatSource::Storage(ref mut storage) => match storage {
                HeatSourceWithStorageTank::Immersion(imm) => imm
                    .lock()
                    .demand_energy(energy_demand, simulation_time_iteration),
                HeatSourceWithStorageTank::Solar(ref solar) => Ok(solar
                    .lock()
                    .demand_energy(energy_demand, simulation_time_iteration.index)),
            },
            HeatSource::Wet(ref mut wet) => {
                wet.demand_energy(energy_demand, None, temp_return, simulation_time_iteration)
            }
        }
    }

    pub fn as_immersion_heater(&self) -> Option<Arc<Mutex<ImmersionHeater>>> {
        match self {
            HeatSource::Storage(HeatSourceWithStorageTank::Immersion(immersion)) => {
                Some(immersion.clone())
            }
            _ => None,
        }
    }
}
#[derive(Clone, Debug)]
pub(crate) enum WetHeatSource {
    HeatPump(Arc<Mutex<HeatPump>>),
    Boiler(Arc<RwLock<Boiler>>),
    Hiu(Arc<Mutex<HeatNetwork>>),
    HeatBattery(Arc<RwLock<HeatBatteryPcm>>),
}

impl WetHeatSource {
    pub fn timestep_end(&self, simtime: SimulationTimeIteration) -> anyhow::Result<()> {
        match self {
            WetHeatSource::HeatPump(heat_pump) => heat_pump.lock().timestep_end(simtime.index)?,
            WetHeatSource::Boiler(boiler) => boiler.write().timestep_end(simtime)?,
            WetHeatSource::Hiu(heat_network) => heat_network.lock().timestep_end(simtime.index),
            WetHeatSource::HeatBattery(heat_battery) => {
                heat_battery.read().timestep_end(simtime.index)?
            }
        }

        Ok(())
    }

    pub(crate) fn create_service_hot_water_combi(
        &mut self,
        boiler_data: HotWaterSourceDetails,
        service_name: &str,
        temp_hot_water: f64,
        cold_feed: WaterSourceWithTemperature,
    ) -> anyhow::Result<BoilerServiceWaterCombi> {
        match self {
            WetHeatSource::HeatPump(heat_pump) => heat_pump.lock().create_service_hot_water_combi(
                boiler_data,
                service_name,
                temp_hot_water,
                cold_feed,
            ),
            WetHeatSource::Boiler(ref mut boiler) => Boiler::create_service_hot_water_combi(
                boiler.clone(),
                boiler_data,
                service_name,
                temp_hot_water,
                cold_feed,
            )
            .map_err(|err| anyhow!(format!("{err}"))),
            _ => {
                bail!("Expect to only be able to create a hot water combi service for boilers and heat pumps.")
            }
        }
    }

    fn output_detailed_results(
        &self,
        hot_water_energy_output: &[f64],
    ) -> Option<(ResultsPerTimestep, ResultsAnnual)> {
        match self {
            WetHeatSource::HeatPump(heat_pump) => {
                let (results_per_timestep, results_annual) =
                    heat_pump.lock().clone().output_detailed_results(
                        hot_water_energy_output
                            .iter()
                            .map(|&x| x.into())
                            .collect_vec(),
                    );
                Some((
                    crate::core::heating_systems::heat_pump::to_corpus_results_per_timestep(
                        results_per_timestep,
                    ),
                    crate::core::heating_systems::heat_pump::to_corpus_results_annual(
                        results_annual,
                    ),
                ))
            }
            WetHeatSource::HeatBattery(heat_battery) => heat_battery
                .read()
                .output_detailed_results(hot_water_energy_output)
                .ok()
                .map(|(results_per_timestep, results_annual)| {
                    (
                        crate::core::heating_systems::heat_battery_pcm::to_corpus_results_per_timestep(
                            results_per_timestep,
                        ),
                        crate::core::heating_systems::heat_battery_pcm::to_corpus_results_annual(
                            results_annual,
                        ),
                    )
                }),
            _ => None,
        }
    }
}

pub type ResultsPerTimestep =
    IndexMap<String, IndexMap<(String, Option<String>), Vec<ResultParamValue>>>;
pub type ResultsAnnual =
    IndexMap<String, IndexMap<(String, Option<String>), Vec<ResultParamValue>>>;

#[derive(Clone, Debug)]
pub enum ResultParamValue {
    String(String),
    Number(f64),
    Boolean(bool),
}

impl ResultParamValue {
    pub fn as_f64(&self) -> f64 {
        match self {
            ResultParamValue::Number(num) => *num,
            _ => 0.,
        }
    }
}

impl Sum for ResultParamValue {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        // expect everything to be a number (or infer zero), return a number variant
        Self::Number(
            iter.map(|value| match value {
                ResultParamValue::Number(num) => num,
                _ => 0.,
            })
            .sum::<f64>(),
        )
    }
}

impl Add for ResultParamValue {
    type Output = ResultParamValue;

    fn add(self, rhs: Self) -> Self::Output {
        ResultParamValue::Number(self.as_f64() + rhs.as_f64())
    }
}

impl Add for &ResultParamValue {
    type Output = ResultParamValue;

    fn add(self, rhs: Self) -> Self::Output {
        ResultParamValue::Number(self.as_f64() + rhs.as_f64())
    }
}

impl AddAssign for ResultParamValue {
    fn add_assign(&mut self, rhs: Self) {
        *self = Self::Number(self.as_f64() + rhs.as_f64());
    }
}

impl Div for ResultParamValue {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self::Number(self.as_f64() / rhs.as_f64())
    }
}

impl PartialEq<f64> for ResultParamValue {
    fn eq(&self, other: &f64) -> bool {
        self.as_f64() == *other
    }
}

impl From<f64> for ResultParamValue {
    fn from(value: f64) -> Self {
        Self::Number(value)
    }
}

impl From<ResultParamValue> for String {
    fn from(value: ResultParamValue) -> Self {
        value.to_string().into()
    }
}

impl From<&ResultParamValue> for String {
    fn from(value: &ResultParamValue) -> Self {
        value.to_string().into()
    }
}

impl From<String> for ResultParamValue {
    fn from(value: String) -> Self {
        Self::String(value)
    }
}

impl From<StringOrNumber> for ResultParamValue {
    fn from(value: StringOrNumber) -> Self {
        match value {
            StringOrNumber::String(string) => Self::String(string),
            StringOrNumber::Float(number) => Self::Number(number),
            StringOrNumber::Integer(number) => Self::Number(number as f64),
        }
    }
}

impl Display for ResultParamValue {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                ResultParamValue::String(string) => string.to_string(),
                ResultParamValue::Number(number) => number.to_string(),
                ResultParamValue::Boolean(boolean) => boolean.to_string(),
            }
        )
    }
}

fn heat_source_wet_from_input(
    energy_supply_name: &str,
    input: HeatSourceWetDetails,
    external_conditions: Arc<ExternalConditions>,
    simulation_time: Arc<SimulationTimeIterator>,
    mechanical_ventilations: &[Arc<MechanicalVentilation>],
    number_of_zones: usize,
    temp_internal_air_fn: TempInternalAirFn,
    controls: &Controls,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    detailed_output_heating_cooling: bool,
) -> anyhow::Result<WetHeatSource> {
    match &input {
        HeatSourceWetDetails::HeatPump {
            source_type,
            energy_supply_heat_network,
            energy_supply,
            boiler,
            ..
        } => {
            let throughput_exhaust_air = if source_type.is_exhaust_air() {
                // Check that ventilation system is compatible with exhaust air HP
                // the Python code here will assign to throughput_exhaust_air on the last iteration, though mech vents collection is not ordered - this may be erroneous
                // reported to BRE as possible bug https://dev.azure.com/BreGroup/SAP%2011/_workitems/edit/45508
                let mut throughput_exhaust_air: Option<f64> = Default::default();
                for mech_vent in mechanical_ventilations.iter() {
                    match mech_vent.vent_type() {
                        MechVentType::IntermittentMev
                        | MechVentType::DecentralisedContinuousMev => {
                            bail!("Exhaust air heat pump does not work with Intermittent MEV or Decentralised continuous MEV.")
                        }
                        _ => {
                            throughput_exhaust_air =
                                Some(mech_vent.as_ref().design_outdoor_air_flow_rate_m3_h);
                            // In Python code data is appended here to a project instance variable called __heat_source_wet_names_requiring_overvent
                            // but this is never read, so we are leaving out the implementation here
                        }
                    }
                }
                throughput_exhaust_air
            } else {
                None
            };

            let energy_supply_heat_source = if matches!(
                source_type,
                HeatPumpSourceType::HeatNetwork
            ) {
                let energy_supply_heat_network = energy_supply_heat_network.as_ref().ok_or_else(|| anyhow!("A heat pump with a heat network source is expected to reference an energy supply for the heat network."))?;
                Some(energy_supplies.get(energy_supply_heat_network).ok_or_else(|| anyhow!("A heat network with a heat network source references an undeclared energy supply '{energy_supply_heat_network}'."))?.clone())
            } else {
                Some(energy_supplies.get(ENERGY_FROM_ENVIRONMENT_SUPPLY_NAME).ok_or_else(|| anyhow!("A heat pump with a '{source_type:?}' source is expected to have an energy supply representing the environment set up."))?.clone())
            };

            let (boiler, cost_schedule_hybrid_hp) = if let Some(boiler) = boiler {
                let energy_supply_boiler = energy_supplies
                    .get(&boiler.energy_supply)
                    .ok_or_else(|| {
                        anyhow!(
                            "A boiler references an undeclared energy supply '{}'.",
                            boiler.energy_supply
                        )
                    })?
                    .clone();
                let energy_supply_aux_boiler = energy_supplies
                    .get(&boiler.energy_supply_aux)
                    .ok_or_else(|| {
                        anyhow!(
                            "A boiler references an undeclared energy supply '{}'.",
                            boiler.energy_supply_aux
                        )
                    })?
                    .clone();
                let energy_supply_conn_aux_boiler = EnergySupply::connection(
                    energy_supply_aux_boiler,
                    format!("Boiler_auxiliary: {energy_supply_name}").as_str(),
                )?;

                let cost_schedule_hybrid_hp = boiler.cost_schedule_hybrid.clone();

                let boiler = Boiler::new(
                    boiler.as_ref().into(),
                    energy_supply_boiler,
                    energy_supply_conn_aux_boiler,
                    external_conditions.clone(),
                    simulation_time.step_in_hours(),
                )?;

                (Some(boiler), cost_schedule_hybrid_hp)
            } else {
                Default::default()
            };

            let energy_supply = energy_supplies
                .get(energy_supply)
                .ok_or_else(|| {
                    anyhow!("Heat pump references an undeclared energy supply '{energy_supply}'.")
                })?
                .clone();
            let energy_supply_conn_name_auxiliary =
                format!("HeatPump_auxiliary: {energy_supply_name}");

            Ok(WetHeatSource::HeatPump(Arc::new(Mutex::new(
                HeatPump::new(
                    &input,
                    energy_supply,
                    &energy_supply_conn_name_auxiliary,
                    simulation_time.step_in_hours(),
                    external_conditions.clone(),
                    number_of_zones,
                    throughput_exhaust_air,
                    energy_supply_heat_source,
                    detailed_output_heating_cooling,
                    boiler.map(|boiler: Boiler| Arc::new(RwLock::new(boiler))),
                    cost_schedule_hybrid_hp,
                    temp_internal_air_fn,
                )?,
            ))))
        }
        HeatSourceWetDetails::Boiler {
            energy_supply,
            energy_supply_aux: energy_supply_auxiliary,
            ..
        } => {
            let energy_supply = energy_supplies
                .get(energy_supply)
                .ok_or_else(|| {
                    anyhow!("Boiler references undeclared energy supply '{energy_supply}'.")
                })?
                .clone();
            let energy_supply_aux = energy_supplies.get(energy_supply_auxiliary).ok_or_else(|| anyhow!("Boiler references undeclared auxiliary energy supply '{energy_supply_auxiliary}'."))?.clone();
            let aux_supply_name = format!("Boiler_auxiliary: {energy_supply_name}");
            let energy_supply_conn_aux =
                EnergySupply::connection(energy_supply_aux.clone(), aux_supply_name.as_str())?;

            Ok(WetHeatSource::Boiler(Arc::new(RwLock::new(
                Boiler::new(
                    input,
                    energy_supply,
                    energy_supply_conn_aux,
                    external_conditions.clone(),
                    simulation_time.step_in_hours(),
                )
                .expect("could not construct boiler value from provided data"),
            ))))
        }
        HeatSourceWetDetails::Hiu {
            power_max,
            hiu_daily_loss,
            building_level_distribution_losses,
            energy_supply,
            ..
        } => {
            let energy_supply = energy_supplies
                .get(energy_supply)
                .ok_or_else(|| {
                    anyhow!("HIU references an undeclared energy supply '{energy_supply}'.")
                })?
                .clone();
            let energy_supply_conn_name_auxiliary =
                String::from(["HeatNetwork_auxiliary: ", energy_supply_name].concat());
            let energy_supply_conn_name_building_level_distribution_losses = String::from(
                [
                    "HeatNetwork_building_level_distribution_losses: ",
                    energy_supply_name,
                ]
                .concat(),
            );

            Ok(WetHeatSource::Hiu(Arc::new(Mutex::new(HeatNetwork::new(
                *power_max,
                *hiu_daily_loss,
                0.06, // TODO update during 1.0.0a1 migration
                0.,   // TODO update during 1.0.0a1 migration
                *building_level_distribution_losses,
                energy_supply,
                energy_supply_conn_name_auxiliary,
                energy_supply_conn_name_building_level_distribution_losses,
                simulation_time.step_in_hours(),
            )))))
        }
        HeatSourceWetDetails::HeatBattery {
            battery:
                HeatBatteryInput::Pcm {
                    control_charge,
                    energy_supply,
                    ..
                },
        } => {
            let energy_supply = energy_supplies
                .get(energy_supply)
                .ok_or_else(|| {
                    anyhow!(
                        "Heat battery references an undeclared energy supply '{energy_supply}'."
                    )
                })?
                .clone();
            let energy_supply_conn =
                EnergySupply::connection(energy_supply.clone(), energy_supply_name)?;

            let heat_source =
                WetHeatSource::HeatBattery(Arc::new(RwLock::new(HeatBatteryPcm::new(
                    &input,
                    controls
                        .get_with_string(control_charge)
                        .unwrap_or_else(|| {
                            panic!(
                            "expected a control to be registered with the name '{control_charge}'"
                        )
                        })
                        .clone(),
                    energy_supply,
                    energy_supply_conn,
                    simulation_time,
                    Some(8),
                    Some(20.),
                    Some(10.),
                    Some(53.),
                    Some(detailed_output_heating_cooling),
                ))));
            Ok(heat_source)
        }
        HeatSourceWetDetails::HeatBattery {
            battery: HeatBatteryInput::DryCore { .. },
        } => unimplemented!(), // TODO: presumably dry core batteries will be mapped during migration to 1.0.0a1
    }
}

fn heat_source_from_input(
    name: &str,
    input: &HeatSourceInput,
    cold_water_source: &WaterSourceWithTemperature,
    volume: f64,
    daily_losses: f64,
    heat_exchanger_surface_area: Option<f64>,
    wet_heat_sources: &IndexMap<String, WetHeatSource>,
    simulation_time: &SimulationTimeIterator,
    controls: &Controls,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    temp_internal_air_fn: TempInternalAirFn,
    external_conditions: Arc<ExternalConditions>,
) -> anyhow::Result<(HeatSource, String)> {
    match input {
        HeatSourceInput::ImmersionHeater {
            power,
            control_min,
            control_max,
            energy_supply,
            ..
        } => {
            let energy_supply = energy_supplies.get(energy_supply).ok_or_else(|| anyhow!("Immersion heater references an undeclared energy supply '{energy_supply}'."))?.clone();
            let energy_supply_conn = EnergySupply::connection(energy_supply.clone(), name)?;
            let control_min = control_min
                .as_ref()
                .and_then(|ctrl| controls.get_with_string(ctrl)).ok_or_else(|| anyhow!("A control indicated by `control_min` is needed for an ImmersionHeater object."))?;
            let control_max = control_max
                .as_ref()
                .and_then(|ctrl| controls.get_with_string(ctrl)).ok_or_else(|| anyhow!("A control indicated by `control_max` is needed for an ImmersionHeater object."))?;

            Ok((
                HeatSource::Storage(HeatSourceWithStorageTank::Immersion(Arc::new(Mutex::new(
                    ImmersionHeater::new(
                        *power,
                        energy_supply_conn,
                        simulation_time.step_in_hours(),
                        Some(control_min),
                        Some(control_max),
                    ),
                )))),
                name.into(),
            ))
        }
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
            energy_supply,
            control_max,
            ..
        } => {
            let energy_supply = energy_supplies.get(energy_supply).ok_or_else(|| anyhow!("Solar thermal system references an undeclared energy supply '{energy_supply}'."))?.clone();
            let energy_supply_conn = EnergySupply::connection(energy_supply.clone(), name)?;
            let control_max = controls
                .get_with_string(control_max)
                .ok_or_else(|| anyhow!("A control indicated by `control_max` is needed for a SolarThermalSystem object."))?;

            let energy_supply_from_environment = energy_supplies
                .get(ENERGY_FROM_ENVIRONMENT_SUPPLY_NAME)
                .ok_or_else(|| anyhow!("An energy supply representing energy from the environment is expected to have been set up."))?;
            let energy_supply_from_environment_conn =
                EnergySupply::connection(energy_supply_from_environment.clone(), name)?;
            let contents = &WATER;

            Ok((
                HeatSource::Storage(HeatSourceWithStorageTank::Solar(Arc::new(Mutex::new(
                    SolarThermalSystem::new(
                        *solar_cell_location,
                        *area_module,
                        *modules,
                        *peak_collector_efficiency,
                        *incidence_angle_modifier,
                        *first_order_hlc,
                        *second_order_hlc,
                        *collector_mass_flow_rate,
                        *power_pump,
                        *power_pump_control,
                        energy_supply_conn,
                        *tilt,
                        *orientation,
                        *solar_loop_piping_hlc,
                        external_conditions.clone(),
                        temp_internal_air_fn,
                        simulation_time.step_in_hours(),
                        control_max,
                        **contents,
                        Some(energy_supply_from_environment_conn),
                    ),
                )))),
                name.into(),
            ))
        }
        HeatSourceInput::ServiceWaterRegular {
            name,
            control_min,
            control_max,
            temp_flow_limit_upper,
            ..
        } => {
            let energy_supply_conn_name: String = format!("{name}_water_heating").into();
            let heat_source_wet = wet_heat_sources
                .get(name)
                .ok_or_else(|| {
                    anyhow!("Expected a wet heat source registered with the name '{name}'.")
                })?
                .clone();
            let control_min = control_min
                .as_ref()
                .and_then(|ctrl| controls.get_with_string(ctrl)).ok_or_else(|| anyhow!("A control indicated by `control_min` is needed for wet heat source with the name '{name}'"))?;
            let control_max = control_max
                .as_ref()
                .and_then(|ctrl| controls.get_with_string(ctrl)).ok_or_else(|| anyhow!("A control indicated by `control_max` is needed for wet heat source with the name '{name}'"))?;
            let temp_flow_limit_upper = temp_flow_limit_upper.ok_or_else(|| {
                anyhow!(
                    "A temp_flow_limit_upper is needed for wet heat source with the name '{name}'"
                )
            })?;
            let mut heat_source_wet_clone = heat_source_wet.clone();

            Ok((
                match heat_source_wet_clone {
                    WetHeatSource::HeatPump(heat_pump) => HeatSource::Wet(Box::new(
                        HeatSourceWet::HeatPumpWater(HeatPump::create_service_hot_water(
                            heat_pump.clone(),
                            &energy_supply_conn_name,
                            temp_flow_limit_upper,
                            Arc::new(cold_water_source.clone()),
                            control_min,
                            control_max,
                        )?),
                    )),
                    WetHeatSource::Boiler(ref mut boiler) => HeatSource::Wet(Box::new(
                        HeatSourceWet::WaterRegular(Boiler::create_service_hot_water_regular(
                            boiler.clone(),
                            energy_supply_conn_name.as_str(),
                            control_min,
                            control_max,
                        )?),
                    )),
                    WetHeatSource::Hiu(heat_network) => {
                        HeatSource::Wet(Box::new(HeatSourceWet::HeatNetworkWaterStorage(
                            HeatNetwork::create_service_hot_water_storage(
                                heat_network,
                                &energy_supply_conn_name,
                                control_min,
                                control_max,
                            ),
                        )))
                    }
                    WetHeatSource::HeatBattery(battery) => {
                        HeatSource::Wet(Box::new(HeatSourceWet::HeatBatteryHotWater(
                            HeatBatteryPcm::create_service_hot_water_regular(
                                battery,
                                &energy_supply_conn_name,
                                cold_water_source.clone(),
                                Some(control_min),
                                Some(control_max),
                            )?,
                        )))
                    }
                },
                energy_supply_conn_name,
            ))
        }
        HeatSourceInput::HeatPumpHotWaterOnly {
            power_max,
            vol_hw_daily_average,
            tank_volume_declared,
            heat_exchanger_surface_area_declared,
            daily_losses_declared,
            ref test_data,
            energy_supply,
            control_min,
            control_max,
            in_use_factor_mismatch,
            ..
        } => {
            let energy_supply = energy_supplies.get(energy_supply).ok_or_else(|| anyhow!("Heat pump (hot water only) references an undeclared energy supply '{energy_supply}'."))?.clone();
            let energy_supply_conn_name = name;
            let energy_supply_connection =
                EnergySupply::connection(energy_supply.clone(), energy_supply_conn_name)?;
            let control_min = controls
                .get_with_string(control_min)
                .ok_or_else(|| anyhow!("A control indicated by `control_min` is needed for a HeatPumpHotWaterOnly object."))?;
            let control_max = controls
                .get_with_string(control_max)
                .ok_or_else(|| anyhow!("A control indicated by `control_max` is needed for a HeatPumpHotWaterOnly object."))?;

            Ok((
                HeatSource::Wet(Box::new(HeatSourceWet::HeatPumpWaterOnly(
                    HeatPumpHotWaterOnly::new(
                        *power_max,
                        energy_supply_connection,
                        test_data,
                        *vol_hw_daily_average,
                        volume,
                        daily_losses,
                        heat_exchanger_surface_area.ok_or(anyhow!("A heat exchanger surface_area is expected for a HeatPumpHotWaterOnly heat source"))?,
                        *in_use_factor_mismatch,
                        *tank_volume_declared,
                        *heat_exchanger_surface_area_declared,
                        *daily_losses_declared,
                        simulation_time.step_in_hours(),
                        control_min,
                        control_max,
                    ),
                ))),
                energy_supply_conn_name.into(),
            ))
        }
    }
}

#[derive(Debug)]
pub(crate) enum HotWaterSource {
    PreHeated(HotWaterStorageTank),
    CombiBoiler(BoilerServiceWaterCombi),
    PointOfUse(PointOfUse),
    HeatNetwork(HeatNetworkServiceWaterDirect),
    HeatBattery(HeatBatteryPcmServiceWaterDirect),
}

impl HotWaterSource {
    pub fn get_cold_water_source(&self) -> WaterSourceWithTemperature {
        match self {
            HotWaterSource::PreHeated(source) => match source {
                HotWaterStorageTank::StorageTank(storage_tank) => {
                    storage_tank.read().get_cold_water_source().clone()
                }
                HotWaterStorageTank::SmartHotWaterTank(smart_storage_tank) => {
                    smart_storage_tank.read().get_cold_water_source().clone()
                }
            },
            HotWaterSource::CombiBoiler(source) => source.get_cold_water_source().clone(),
            HotWaterSource::PointOfUse(source) => source.get_cold_water_source().clone(),
            HotWaterSource::HeatNetwork(source) => source.get_cold_water_source().clone(),
            HotWaterSource::HeatBattery(source) => source.get_cold_water_source().clone(),
        }
    }

    pub(crate) fn temp_hot_water(&self) -> anyhow::Result<f64> {
        Ok(match self {
            HotWaterSource::PreHeated(source) => match source {
                HotWaterStorageTank::StorageTank(_storage_tank) => {
                    unimplemented!("WIP - storage tank migration")
                    // storage_tank.read().get_temp_hot_water()
                }
                HotWaterStorageTank::SmartHotWaterTank(_smart_storage_tank) => {
                    unimplemented!("WIP - storage tank migration")
                    // smart_storage_tank.read().get_temp_hot_water()
                }
            },
            HotWaterSource::CombiBoiler(_combi) => {
                todo!("Probably gets removed/moved as part of migration to 1.0.01a")
            } // combi.get_temp_hot_water(),
            HotWaterSource::PointOfUse(_point_of_use) => {
                todo!("Probably gets removed/moved as part of migration to 1.0.01a")
            }
            HotWaterSource::HeatNetwork(_heat_network) => {
                todo!("Probably gets removed/moved as part of migration to 1.0.01a")
            }
            HotWaterSource::HeatBattery(_source) => {
                todo!("Probably gets removed/moved as part of migration to 1.0.01a")
            }
        })
    }

    pub fn demand_hot_water(
        &self,
        _vol_demand_target: IndexMap<DemandVolTargetKey, VolumeReference>,
        _simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        Ok(match self {
            HotWaterSource::PreHeated(_) => {
                // StorageTank does not match the same method signature or return type as all other Hot Water sources
                panic!("demand_hot_water for HotWaterSource::StorageTank should be called directly on the HotWaterSource::StorageTank");
            }
            HotWaterSource::CombiBoiler(ref _source) => {
                todo!("To do, this probably gets removed as part of migration to 1.0.0a1");
            }
            // source
            //     .demand_hot_water(vol_demand_target, simulation_time_iteration)
            //     .expect("Combi boiler could not calc demand hot water."),
            HotWaterSource::PointOfUse(ref _source) => {
                todo!("To do, this probably gets removed as part of migration to 1.0.0a1");
            }
            HotWaterSource::HeatNetwork(ref _source) => {
                todo!("To do, this probably gets removed as part of migration to 1.0.0a1");
            }
            HotWaterSource::HeatBattery(_source) => {
                todo!("To do, this probably gets removed as part of migration to 1.0.0a1");
            }
        })
    }
}

fn hot_water_source_from_input(
    source_name: String,
    input: &HotWaterSourceDetails,
    cold_water_sources: &ColdWaterSources,
    pre_heated_water_sources: &IndexMap<String, HotWaterStorageTank>,
    wet_heat_sources: &mut IndexMap<String, WetHeatSource>,
    wwhrs: &IndexMap<String, Arc<Mutex<Wwhrs>>>,
    controls: &Controls,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    diverter_types: &DiverterTypes,
    diverters: &mut Vec<Arc<RwLock<PVDiverter>>>,
    temp_internal_air_fn: TempInternalAirFn,
    simulation_time: &SimulationTimeIterator,
    external_conditions: Arc<ExternalConditions>,
    detailed_output_heating_cooling: bool,
) -> anyhow::Result<(HotWaterSource, Vec<String>)> {
    let mut energy_supply_conn_names = vec![];
    let cloned_input = input.clone();

    let cold_water_source_for_hot_water_tank =
        |cold_water_source_type: &str| -> anyhow::Result<WaterSourceWithTemperature> {
            pre_heated_water_sources
                .get(cold_water_source_type)
                .map(|source| WaterSourceWithTemperature::Preheated(source.clone()))
                .or(wwhrs.get(cold_water_source_type).map(|source| WaterSourceWithTemperature::Wwhrs(source.clone())))
                .ok_or_else(|| anyhow!("Could not find pre-heated or WWHRS water source for name '{cold_water_source_type}'"))
        };

    let mut heat_sources_for_hot_water_tank =
        |cold_water_source: WaterSourceWithTemperature,
         heat_exchanger_surface_area: &Option<f64>,
         heat_source: &IndexMap<std::string::String, HeatSourceInput>,
         volume: &f64,
         daily_losses: &f64|
         -> anyhow::Result<IndexMap<String, PositionedHeatSource>> {
            let mut heat_sources: IndexMap<String, PositionedHeatSource> = Default::default();

            let heat_exchanger_surface_area =
                heat_exchanger_surface_area.and_then(|surface_area| {
                    heat_source
                        .values()
                        .any(|source| {
                            matches!(source, HeatSourceInput::HeatPumpHotWaterOnly { .. })
                        })
                        .then_some(surface_area)
                });

            // With pre-heated tanks we allow now tanks not to have a heat source as the 'cold' feed
            // could be a pre-heated source or wwhr that might be enough
            let mut used_heat_source_names: Vec<String> = Default::default();
            for (name, hs) in heat_source {
                let name = String::from(name);
                if used_heat_source_names.contains(&name) {
                    return Err(anyhow!("Duplicate heat source name detected: {name}"));
                }
                used_heat_source_names.push(name.clone());

                let heater_position = hs.heater_position();
                let thermostat_position = match input {
                    HotWaterSourceDetails::StorageTank { .. } => hs.thermostat_position(),
                    HotWaterSourceDetails::SmartHotWaterTank { .. } => None,
                    _ => {
                        unreachable!()
                    }
                };

                let (heat_source, energy_supply_conn_name) = heat_source_from_input(
                    name.as_str(),
                    hs,
                    &cold_water_source,
                    *volume,
                    *daily_losses,
                    heat_exchanger_surface_area,
                    wet_heat_sources,
                    simulation_time,
                    controls,
                    energy_supplies,
                    temp_internal_air_fn.clone(),
                    external_conditions.clone(),
                )?;
                let heat_source = Arc::new(Mutex::new(heat_source));

                heat_sources.insert(
                    name,
                    PositionedHeatSource {
                        heat_source: heat_source.clone(),
                        heater_position,
                        thermostat_position,
                    },
                );
                energy_supply_conn_names.push(energy_supply_conn_name);
            }

            Ok(heat_sources)
        };

    let mut connect_diverter_for_hot_water_tank = |energy_supplies: IndexMap<
        String,
        Arc<RwLock<EnergySupply>>,
    >,
                                                   heat_source: &IndexMap<
        std::string::String,
        HeatSourceInput,
    >,
                                                   heat_sources: &IndexMap<
        String,
        PositionedHeatSource,
    >,
                                                   pre_heated_tank: HotWaterStorageTank|
     -> anyhow::Result<()> {
        for (heat_source_name, hs) in heat_source {
            let heat_source_name = String::from(heat_source_name);
            let energy_supply_name = hs.energy_supply_name();
            if let Some(diverter) = diverter_types.get(energy_supply_name) {
                if diverter.heat_source.matches(&heat_source_name) {
                    let energy_supply = energy_supplies.get(energy_supply_name).ok_or_else(|| anyhow!("Heat source references an undeclared energy supply '{energy_supply_name}'."))?.clone();

                    let positioned_heat_source = &heat_sources.get(&heat_source_name);
                    let immersion_heater = positioned_heat_source
                        .unwrap()
                        .heat_source
                        .lock()
                        .as_immersion_heater();

                    if let Some(im) = immersion_heater {
                        let control_max = controls.get_with_string(&diverter.control_max);
                        let pv_diverter = PVDiverter::new(
                            &pre_heated_tank,
                            im,
                            heat_source_name.clone(),
                            control_max,
                        );
                        energy_supply
                            .write()
                            .connect_diverter(pv_diverter.clone())
                            .unwrap();
                        diverters.push(pv_diverter);
                    }
                }
            }
        }
        Ok(())
    };

    let hot_water_source = match input {
        HotWaterSourceDetails::StorageTank {
            volume,
            daily_losses,
            heat_exchanger_surface_area,
            init_temp,
            cold_water_source: cold_water_source_type,
            primary_pipework,
            heat_source,
            ..
        } => {
            let cold_water_source = cold_water_source_for_hot_water_tank(cold_water_source_type)?;
            // At this point in the Python, the internal_diameter and external_diameter fields on
            // primary_pipework are updated, this is done in Pipework.rs in the Rust
            let primary_pipework_lst = primary_pipework.as_ref();
            let heat_sources = heat_sources_for_hot_water_tank(
                cold_water_source.clone(),
                heat_exchanger_surface_area,
                heat_source,
                volume,
                daily_losses,
            )?;

            let storage_tank = Arc::new(RwLock::new(StorageTank::new(
                *volume,
                *daily_losses,
                *init_temp,
                cold_water_source,
                simulation_time.step_in_hours(),
                heat_sources.clone(),
                temp_internal_air_fn,
                external_conditions,
                Some(24),
                primary_pipework_lst,
                *WATER,
                None,
                None,
                None,
                detailed_output_heating_cooling,
            )?));

            connect_diverter_for_hot_water_tank(
                energy_supplies.clone(),
                heat_source,
                &heat_sources,
                HotWaterStorageTank::StorageTank(storage_tank.clone()),
            )?;

            HotWaterSource::PreHeated(HotWaterStorageTank::StorageTank(storage_tank))
        }
        HotWaterSourceDetails::SmartHotWaterTank {
            volume,
            daily_losses,
            init_temp,
            power_pump_kw,
            max_flow_rate_pump_l_per_min,
            temp_usable,
            temp_setpnt_max,
            cold_water_source,
            primary_pipework,
            heat_source,
            energy_supply_pump,
            ..
        } => {
            let cold_water_source_type = cold_water_source;
            let cold_water_source = cold_water_source_for_hot_water_tank(cold_water_source_type)?;

            // At this point in the Python, the internal_diameter and external_diameter fields on
            // primary_pipework are updated, this is done in Pipework.rs in the Rust
            let primary_pipework_lst = primary_pipework.as_ref();
            let heat_sources = heat_sources_for_hot_water_tank(
                cold_water_source.clone(),
                &None,
                heat_source,
                volume,
                daily_losses,
            )?;

            if !heat_sources
                .values()
                .map(|source| source.heater_position)
                .all_equal()
            {
                return Err(anyhow!(
                    "For SmartHotWaterTank, heater position must be the same for all heat sources"
                ));
            }

            let energy_supply_pump = energy_supplies[energy_supply_pump].clone();
            let pump_source_name = format!("Smart_hot_water_tank_pump: {source_name}");
            let energy_supply_conn_pump =
                EnergySupply::connection(energy_supply_pump, &pump_source_name)?;
            let smart_hot_water_tank = Arc::new(RwLock::new(SmartHotWaterTank::new(
                *volume,
                *daily_losses,
                *init_temp,
                *power_pump_kw,
                *max_flow_rate_pump_l_per_min,
                *temp_usable,
                controls.get_with_string(temp_setpnt_max).ok_or_else(|| anyhow!("A control indicated by `temp_setpnt_max` is needed for a SmartHotWaterTank object."))?,
                cold_water_source.clone(),
                simulation_time.step_in_hours(),
                heat_sources.clone(),
                temp_internal_air_fn,
                external_conditions,
                detailed_output_heating_cooling.into(),
                100.into(),
                primary_pipework_lst
                    .map(|primary_pipework| {
                        anyhow::Ok(
                            primary_pipework
                                .iter()
                                .copied()
                                .map(WaterPipework::try_from)
                                .collect::<Result<Vec<_>, _>>()?,
                        )
                    })
                    .transpose()?
                    .as_ref(),
                energy_supply_conn_pump,
                None,
            )?));

            connect_diverter_for_hot_water_tank(
                energy_supplies.clone(),
                heat_source,
                &heat_sources,
                HotWaterStorageTank::SmartHotWaterTank(smart_hot_water_tank.clone()),
            )?;

            HotWaterSource::PreHeated(HotWaterStorageTank::SmartHotWaterTank(smart_hot_water_tank))
        }
        HotWaterSourceDetails::CombiBoiler {
            cold_water_source: cold_water_source_type,
            heat_source_wet: heat_source_wet_type,
            setpoint_temp,
            ..
        } => {
            let cold_water_source =
                cold_water_source_for_type(cold_water_source_type, cold_water_sources)?;
            let energy_supply_conn_name =
                String::from([heat_source_wet_type, "_water_heating"].concat());
            energy_supply_conn_names.push(energy_supply_conn_name.clone());
            let heat_source_wet = wet_heat_sources.get_mut(heat_source_wet_type).ok_or_else(|| {
                anyhow!("Expected '{heat_source_wet_type}' to have been defined as a wet heat source")
            })?;
            HotWaterSource::CombiBoiler(
                heat_source_wet
                    .create_service_hot_water_combi(
                        cloned_input,
                        energy_supply_conn_name.as_str(),
                        setpoint_temp.ok_or_else(|| {
                            anyhow!("A setpoint temp was expected on a combi boiler input.")
                        })?,
                        cold_water_source,
                    )
                    .expect("expected to be able to instantiate a combi boiler object"),
            )
        }
        HotWaterSourceDetails::PointOfUse {
            efficiency,
            cold_water_source: cold_water_source_type,
            energy_supply,
            setpoint_temp,
            ..
        } => {
            let energy_supply = energy_supplies.get(energy_supply).ok_or_else(|| anyhow!("Point of use hot water source references an undeclared energy supply '{energy_supply}'."))?.clone();
            let energy_supply_conn_name = source_name;
            energy_supply_conn_names.push(energy_supply_conn_name.clone());
            let energy_supply_conn =
                EnergySupply::connection(energy_supply.clone(), &energy_supply_conn_name)?;
            let cold_water_source =
                cold_water_source_for_type(cold_water_source_type, cold_water_sources)?;
            HotWaterSource::PointOfUse(PointOfUse::new(
                efficiency.ok_or_else(|| anyhow!("An efficiency value was expected on a point of use hot water source input."))?, // TODO: review as part of migration to 1.0.0a1 as efficiency may now be optional
                energy_supply_conn,
                cold_water_source,
                *setpoint_temp,
            ))
        }
        HotWaterSourceDetails::Hiu {
            cold_water_source: cold_water_source_type,
            heat_source_wet: heat_source_wet_type,
            setpoint_temp,
            ..
        } => {
            let energy_supply_conn_name =
                String::from([heat_source_wet_type, "_water_heating"].concat());
            energy_supply_conn_names.push(energy_supply_conn_name.clone());
            let cold_water_source =
                cold_water_source_for_type(cold_water_source_type, cold_water_sources)?;
            let heat_source_wet = match heat_source_wet_type.as_str() {
                "HeatNetwork" => {
                    match wet_heat_sources
                        .get("HeatNetwork")
                        .expect("expected a heat network in this context")
                    {
                        WetHeatSource::Hiu(heat_network) => heat_network.clone(),
                        _ => panic!("expected a heat network in this context"),
                    }
                }
                _ => panic!("expected a heat network in this context"),
            };
            HotWaterSource::HeatNetwork(HeatNetwork::create_service_hot_water_direct(
                heat_source_wet.clone(),
                &energy_supply_conn_name,
                (*setpoint_temp).ok_or_else(|| {
                    anyhow!(
                        "A setpoint_temp value was expected on a point of use hot water source."
                    )
                })?,
                cold_water_source,
            ))
        }
        HotWaterSourceDetails::HeatBattery {
            cold_water_source,
            heat_source_wet: heat_source_wet_name,
            setpoint_temp,
            ..
        } => {
            let energy_supply_conn_name: String =
                format!("{}_water_heating", heat_source_wet_name).into();
            energy_supply_conn_names.push(energy_supply_conn_name.clone());
            let cold_water_source =
                cold_water_source_for_type(cold_water_source, cold_water_sources)?;
            let heat_source_wet = wet_heat_sources.get(heat_source_wet_name).ok_or_else(|| {
                anyhow!(
                    "Expected a wet heat source registered with the name '{heat_source_wet_name}'."
                )
            })?;
            let heat_battery = match heat_source_wet {
                WetHeatSource::HeatBattery(heat_battery) => heat_battery.clone(),
                _ => unreachable!("heat source wet was expected to be a heat battery"),
            };

            HotWaterSource::HeatBattery(HeatBatteryPcm::create_service_hot_water_direct(
                heat_battery,
                &energy_supply_conn_name,
                *setpoint_temp,
                cold_water_source,
            )?)
        }
    };

    Ok((hot_water_source, energy_supply_conn_names))
}

fn cold_water_source_for_type(
    cold_water_source_type: &str,
    cold_water_sources: &ColdWaterSources,
) -> anyhow::Result<WaterSourceWithTemperature> {
    Ok(WaterSourceWithTemperature::ColdWaterSource(
        cold_water_sources
            .get(cold_water_source_type)
            .ok_or_else(|| anyhow!("referenced cold water source was expected to exist"))?
            .clone(),
    ))
}

fn space_heat_systems_from_input(
    input: &SpaceHeatSystemInput,
    controls: &Controls,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    simulation_time: &SimulationTimeIterator,
    heat_sources_wet: &IndexMap<String, WetHeatSource>,
    heat_system_names_requiring_overvent: &mut Vec<String>,
    heat_system_name_for_zone: &IndexMap<String, Vec<String>>,
    zones: &IndexMap<String, Arc<Zone>>,
    heat_sources_wet_with_buffer_tank: &[String],
    external_conditions: Arc<ExternalConditions>,
    detailed_output_heating_cooling: bool,
) -> anyhow::Result<SpaceHeatSystemsWithEnergyConnections> {
    let mut energy_conn_names_for_systems: IndexMap<String, String> = Default::default();
    let space_heat_systems = input
        .iter()
        .filter(|(system_name, _)| heat_system_name_for_zone.values().flatten().any(|heat_system_name| heat_system_name == system_name.as_str()))
        .map(|(system_name, space_heat_system_details)| {
            let system_name = String::from(system_name);
            Ok((
                system_name.clone(),
                Arc::new(Mutex::new(match space_heat_system_details {
                    SpaceHeatSystemDetails::InstantElectricHeater {
                        rated_power,
                        control,
                        frac_convective,
                        energy_supply,
                        ..
                    } => {
                        let energy_supply = energy_supplies.get(energy_supply).ok_or_else(|| anyhow!("Space heat system references an undeclared energy supply '{energy_supply}'."))?.clone();
                        let energy_supply_conn_name = system_name.clone();
                        energy_conn_names_for_systems.insert(system_name.clone(), energy_supply_conn_name.clone());
                        let energy_supply_conn = EnergySupply::connection(energy_supply, energy_supply_conn_name.as_str()).unwrap();
                        SpaceHeatSystem::Instant(InstantElecHeater::new(
                            *rated_power,
                            *frac_convective,
                            energy_supply_conn,
                            simulation_time.step_in_hours(),
                            controls.get_with_string(control),
                        ))
                    }
                    SpaceHeatSystemDetails::ElectricStorageHeater { pwr_in, rated_power_instant, storage_capacity, air_flow_type, frac_convective, fan_pwr, n_units, energy_supply, zone, control, control_charger, dry_core_min_output, dry_core_max_output, .. } => {
                        let energy_supply = energy_supplies.get(energy_supply).ok_or_else(|| anyhow!("Space heat system references an undeclared energy supply '{energy_supply}'."))?.clone();
                        let energy_supply_conn_name = system_name.clone();
                        energy_conn_names_for_systems.insert(system_name.clone(), energy_supply_conn_name.clone());
                        let energy_supply_conn = EnergySupply::connection(energy_supply, energy_supply_conn_name.as_str()).unwrap();

                        let zone = zones.get(zone).ok_or_else(|| anyhow!("Space heat system references an undeclared zone '{zone}'."))?.clone();
                        let zone_setpoint_init = zone.setpnt_init();
                        let control = controls.get_with_string(control).ok_or_else(|| anyhow!("A control object was expected for an electric storage heater"))?;
                        let charge_control = controls.get_with_string(control_charger).ok_or_else(|| anyhow!("Space heat system references an invalid charge control name '{control_charger}'"))?;
                        SpaceHeatSystem::ElecStorage(ElecStorageHeater::new(*pwr_in, *rated_power_instant, *storage_capacity, *air_flow_type, *frac_convective, *fan_pwr, *n_units, zone_setpoint_init, ZoneTempInternalAir(zone).as_fn(), energy_supply_conn, simulation_time, control, charge_control, dry_core_min_output.clone(), dry_core_max_output.clone(), external_conditions.clone(), Some(detailed_output_heating_cooling))?)
                    }
                    SpaceHeatSystemDetails::WetDistribution { emitters, energy_supply, flow_data, bypass_fraction_recirculated, heat_source, temp_diff_emit_dsgn, control, thermal_mass, ecodesign_controller, design_flow_temp, zone, .. } => {
                        let heat_source_name = &heat_source.name;
                        let temp_flow_limit_upper = &heat_source.temp_flow_limit_upper;

                        let energy_supply_conn_name = String::from([heat_source_name, "_space_heating: ", &system_name].concat());
                        energy_conn_names_for_systems.insert(system_name.clone(), energy_supply_conn_name.clone());

                        let heat_source = heat_sources_wet.get(&heat_source.name).ok_or_else(|| anyhow!("A heat source name provided under the name '{heat_source_name}' was expected when setting up space heat systems in the calculation corpus."))?;
                        let mut with_buffer_tank = false;

                        let control = controls.get_with_string(control).ok_or_else(|| anyhow!("A control object was expected for wet heat source: '{heat_source_name}'"))?;

                        let heat_source_service: SpaceHeatingService =
                            match heat_source {
                                WetHeatSource::HeatPump(heat_pump) => {
                                    // TODO (from Python) If EAHP, feed zone volume into function below

                                    // For HPs, checking if there's a buffer tank to inform both the service space heating
                                    // and the emitters of its presence.
                                    if heat_sources_wet_with_buffer_tank.contains(heat_source_name) {
                                        with_buffer_tank = true;
                                    }

                                    let volume_heated = total_volume_heated_by_system(zones, heat_system_name_for_zone, &system_name);

                                    let heat_source_service = HeatPump::create_service_space_heating(
                                        heat_pump.clone(),
                                        &energy_supply_conn_name,
                                        temp_flow_limit_upper.expect("Expected a temp_flow_limit_upper to be present for a heat pump"),
                                        *temp_diff_emit_dsgn, control,
                                        volume_heated);

                                    if heat_pump.lock().source_is_exhaust_air() {
                                        // Record heating system as potentially requiring overventilation
                                        heat_system_names_requiring_overvent.push((system_name).clone());
                                    }
                                    SpaceHeatingService::HeatPump(heat_source_service)
                                }
                                WetHeatSource::Boiler(boiler) => {
                                    let heat_source_service = Boiler::create_service_space_heating(boiler.clone(), &energy_supply_conn_name, control);
                                    SpaceHeatingService::Boiler(heat_source_service)
                                }
                                WetHeatSource::Hiu(heat_network) => {
                                    let heat_source_service = HeatNetwork::create_service_space_heating(heat_network.clone(), &energy_supply_conn_name, control);
                                    SpaceHeatingService::HeatNetwork(heat_source_service)
                                }
                                WetHeatSource::HeatBattery(heat_battery) => {
                                    let heat_source_service = HeatBatteryPcm::create_service_space_heating(heat_battery.clone(), &energy_supply_conn_name, control)?;
                                    SpaceHeatingService::HeatBattery(heat_source_service)
                                }
                            };

                        let energy_supply_fc_conn  = if energy_supply.is_none() {
                            None
                        } else {
                            let energy_supply_name = energy_supply.clone().unwrap();
                            let energy_supply = energy_supplies.get(&energy_supply_name).ok_or_else(|| anyhow!("Space heat system references an undeclared energy supply '{energy_supply_name}'."))?.clone();
                            let energy_supply_fc_conn_name: String = format!("FC_fan {system_name}").into();
                            energy_conn_names_for_systems.insert(system_name.clone(), energy_supply_fc_conn_name.clone());
                            Some(Arc::new(EnergySupply::connection(energy_supply, energy_supply_fc_conn_name.as_str()).unwrap()))
                        };

                        let space_heater = Emitters::new(
                            *thermal_mass,
                            emitters,
                            &[], // <---- pipework goes here!!
                            *temp_diff_emit_dsgn,
                            matches!(flow_data, FlowData::Variable {..}),
                            if let FlowData::Design {design_flow_rate, ..} = flow_data {
                                Some(*design_flow_rate)
                            } else {
                                None
                            },
                            if let FlowData::Variable {min_flow_rate, ..} = flow_data {
                                Some(*min_flow_rate)
                            } else {
                                None
                            },
                            if let FlowData::Variable {max_flow_rate, ..} = flow_data {
                                Some(*max_flow_rate)
                            } else {
                                None
                            },
                            *bypass_fraction_recirculated,
                            Arc::new(RwLock::new(heat_source_service)),
                            zones.get(zone).ok_or_else(|| anyhow!("Space heat system wet distribution had reference to undeclared zone with name '{zone}'"))?.clone(),
                            // zone area
                            external_conditions.clone(),
                            *ecodesign_controller,
                            *design_flow_temp as f64,
                            20., // replace!!!!! this uses a function to provide an initial value
                            energy_supply_fc_conn,
                            detailed_output_heating_cooling.into(),
                            with_buffer_tank.into(),
                        )?;
                        SpaceHeatSystem::WetDistribution(space_heater)
                    }
                    SpaceHeatSystemDetails::WarmAir {
                        frac_convective,
                        heat_source,
                        control,
                        ..
                    } => {
                        let heat_source_name = &heat_source.name;
                        let energy_supply_conn_name = String::from([heat_source_name, "_space_heating: ", &system_name].concat());
                        energy_conn_names_for_systems.insert(system_name.clone(), energy_supply_conn_name.clone());
                        let heat_source = heat_sources_wet.get(&heat_source.name).ok_or_else(|| anyhow!("A heat source name provided under the name '{heat_source_name}' was expected when setting up space heat systems in the calculation corpus."))?;
                        let control = controls.get_with_string(control).ok_or_else(|| anyhow!("Unknown control object reference '{control}' encountered"))?;

                        match heat_source {
                            WetHeatSource::HeatPump(heat_pump) => {
                                if heat_pump.lock().source_is_exhaust_air() {
                                    heat_system_names_requiring_overvent.push(system_name.clone());
                                }
                                let volume_heated = total_volume_heated_by_system(zones, heat_system_name_for_zone, &system_name);
                                SpaceHeatSystem::WarmAir(HeatPump::create_service_space_heating_warm_air(heat_pump.clone(), &energy_supply_conn_name, control, *frac_convective, volume_heated).unwrap())
                            }
                            _ => panic!("The heat source referenced by details about warm air space heating with the name '{heat_source_name}' was expected to be a heat pump."),
                        }
                    }
                })),
            ))
        })
        .collect::<anyhow::Result<IndexMap<_, _>>>()?;
    Ok((space_heat_systems, energy_conn_names_for_systems))
}

type SpaceHeatSystemsWithEnergyConnections = (
    IndexMap<String, Arc<Mutex<SpaceHeatSystem>>>,
    IndexMap<String, String>,
);

fn space_cool_systems_from_input(
    input: &SpaceCoolSystemInput,
    cool_system_names_for_zone: Vec<&str>,
    controls: &Controls,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<IndexMap<String, AirConditioning>> {
    input
        .iter()
        .filter(|(system_name, _)| cool_system_names_for_zone.contains(&system_name.as_str()))
        .map(|(system_name, space_cool_system_details)| {
            if !matches!(
                space_cool_system_details,
                SpaceCoolSystemDetails::AirConditioning { .. }
            ) {
                unreachable!(
                    "There are no known space cool system types other than air conditioning."
                )
            }
            let SpaceCoolSystemDetails::AirConditioning {
                cooling_capacity,
                efficiency,
                frac_convective,
                control,
                energy_supply,
                ..
            } = space_cool_system_details;
            let energy_supply = energy_supplies.get(energy_supply).ok_or_else(|| anyhow!("Space cool system references an undeclared energy supply '{energy_supply}'."))?.clone();
            let energy_supply_conn_name = system_name;
            let energy_supply_conn =
                EnergySupply::connection(energy_supply, energy_supply_conn_name).unwrap();
            let control = controls.get_with_string(control).ok_or_else(|| anyhow!("The control reference '{control}' was expected to refer to a known control."))?;

            Ok((
                energy_supply_conn_name.into(),
                AirConditioning::new(
                    *cooling_capacity,
                    *efficiency,
                    *frac_convective,
                    energy_supply_conn,
                    simulation_time_iterator.step_in_hours(),
                    control,
                ),
            ))
        })
        .collect::<anyhow::Result<IndexMap<_, _>>>()
}

fn on_site_generation_from_input(
    input: &OnSiteGeneration,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    external_conditions: Arc<ExternalConditions>,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<IndexMap<String, PhotovoltaicSystem>> {
    input
        .iter()
        .map(|(name, generation_details)| {
            Ok((name.into(), {
                let generation_details = match generation_details {
                    PhotovoltaicInputs::DeprecatedStyle(generation_details) => generation_details,
                    PhotovoltaicInputs::WithPanels(_) => todo!("New variant not yet implemented for 1.0.0a1"),
                };
                let PhotovoltaicSystemInput {
                    peak_power,
                    ventilation_strategy,
                    pitch,
                    orientation,
                    base_height,
                    height,
                    width,
                    energy_supply,
                    shading,
                    inverter_peak_power_dc,
                    inverter_peak_power_ac,
                    inverter_is_inside,
                    inverter_type,
                    ..
                } = generation_details;
                let energy_supply = energy_supplies.get(energy_supply).ok_or_else(|| anyhow!("On site generation (photovoltaic) references an undeclared energy supply '{energy_supply}'."))?.clone();
                let energy_supply_conn = EnergySupply::connection(energy_supply, name).unwrap();
                let panels = vec![PhotovoltaicPanel::new(*peak_power, *ventilation_strategy, *pitch, *orientation, *base_height, *height, *width, simulation_time_iterator.step_in_hours(), shading.to_vec())]; // TODO review migration alpha1
                let inverter = Inverter::new(
                    energy_supply_conn,
                    simulation_time_iterator.clone(),
                    *inverter_peak_power_dc,
                    *inverter_peak_power_ac,
                    *inverter_is_inside,
                    *inverter_type,
                );
                PhotovoltaicSystem::new(
                    external_conditions.clone(),
                    panels,
                    inverter
                )
            }))
        })
        .collect::<anyhow::Result<IndexMap<_, _>>>()
}

fn total_volume_heated_by_system(
    zones: &IndexMap<String, Arc<Zone>>,
    heat_system_name_for_zone: &IndexMap<String, Vec<String>>,
    heat_system_name: &str,
) -> f64 {
    zones
        .iter()
        .filter_map(|(z_name, zone)| {
            if let Some(system_names) = heat_system_name_for_zone.get(z_name) {
                (system_names.iter().any(|name| heat_system_name == name)).then(|| zone.volume())
            } else {
                None
            }
        })
        .sum::<f64>()
}

fn required_vent_data_from_input(input: &ControlInput) -> anyhow::Result<Option<RequiredVentData>> {
    input
        .extra
        .get("required_vent")
        .map(|ctrl| {
            anyhow::Ok(RequiredVentData {
                schedule: expand_numeric_schedule(ctrl.numeric_schedule()?),
                start_day: ctrl.start_day()?,
                time_series_step: ctrl.time_series_step()?,
            })
        })
        .transpose()
}

#[derive(Clone, Debug)]
struct RequiredVentData {
    schedule: Vec<Option<f64>>,
    start_day: u32,
    time_series_step: f64,
}
