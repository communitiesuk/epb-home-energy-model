use crate::core::common::WaterSupply;
use crate::core::controls::time_control::{
    ChargeControl, CombinationTimeControl, Control, ControlBehaviour, HeatSourceControl,
    OnOffCostMinimisingTimeControl, OnOffTimeControl, SetpointTimeControl, SmartApplianceControl,
};
use crate::core::cooling_systems::air_conditioning::AirConditioning;
use crate::core::cooling_systems::space_cool_system_base::SpaceCoolSystem;
use crate::core::energy_supply::elec_battery::ElectricBattery;
use crate::core::energy_supply::energy_supply::{
    EnergySupply, EnergySupplyBuilder, EnergySupplyConnection, EnergySupplyTariffInput,
    ENERGY_FROM_ENVIRONMENT_SUPPLY_NAME, UNMET_DEMAND_SUPPLY_NAME,
};
use crate::core::energy_supply::inverter::Inverter;
use crate::core::energy_supply::on_site_generation_base::OnSiteGeneration;
use crate::core::energy_supply::pv::{PhotovoltaicPanel, PhotovoltaicSystem};
use crate::core::heating_systems::boiler::{Boiler, BoilerServiceWaterCombi};
use crate::core::heating_systems::common::{
    HeatBatteryServiceSpace, HeatBatteryWaterService, HeatSourceWet, SpaceHeatSystem,
    SpaceHeatingService,
};
use crate::core::heating_systems::elec_storage_heater::{
    ElecStorageHeater, StorageHeaterDetailedResult,
};
use crate::core::heating_systems::emitters::Emitters;
use crate::core::heating_systems::heat_battery_drycore::{
    HeatBatteryDryCore, HeatBatteryDryCoreServiceWaterDirect,
};
use crate::core::heating_systems::heat_battery_pcm::{
    HeatBatteryPcm, HeatBatteryPcmServiceWaterDirect,
};
use crate::core::heating_systems::heat_network::{HeatNetwork, HeatNetworkServiceWaterDirect};
use crate::core::heating_systems::heat_pump::{
    HeatPump, HeatPumpEmitterType, HeatPumpHotWaterOnly,
};
use crate::core::heating_systems::instant_elec_heater::InstantElecHeater;
use crate::core::heating_systems::point_of_use::PointOfUse;
use crate::core::heating_systems::storage_tank::{
    HeatSourceWithStorageTank, HotWaterStorageTank, ImmersionHeater, PVDiverter,
    PositionedHeatSource, SmartHotWaterTank, SolarThermalSystem, StorageTank,
};
use crate::core::heating_systems::wwhrs::WwhrsInstantaneous;
use crate::core::material_properties::WATER;
use crate::core::schedule::{
    expand_boolean_schedule, expand_events, expand_numeric_schedule, reject_nones, reject_nulls,
    ScheduleEvent, TypedScheduleEvent, WaterScheduleEventType,
};
use crate::core::space_heat_demand::building_element::{
    calculate_cavity_resistance, convert_uvalue_to_resistance, BuildingElement,
    BuildingElementAdjacentConditionedSpace, BuildingElementAdjacentUnconditionedSpaceSimple,
    BuildingElementGround, BuildingElementOpaque, BuildingElementPartyWall,
    BuildingElementTransparent, WindowTreatment, H_CE, H_RE, R_SI_DOWNWARDS, R_SI_HORIZONTAL,
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
use crate::core::units::{kelvin_to_celsius, Orientation360, SECONDS_PER_HOUR, WATTS_PER_KILOWATT};
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::dhw_demand::ELECTRIC_SHOWERS_HWS_NAME;
use crate::core::water_heat_demand::dhw_demand::{DomesticHotWaterDemand, WaterHeatingCalculation};
use crate::core::water_heat_demand::misc::WaterEventResult;
use crate::external_conditions::{create_external_conditions, ExternalConditions};
use crate::hem_core::simulation_time::SimulationTime;
use crate::input::{
    ApplianceGains as ApplianceGainsInput, ApplianceGainsDetails,
    BuildingElement as BuildingElementInput, BuildingElementHeightWidthInput, ChargeLevel,
    ColdWaterSourceDetails, ColdWaterSourceInput, Control as ControlInput, ControlCombinations,
    ControlDetails, EnergyDiverter, EnergySupplyDetails, EnergySupplyInput, FlowData, FuelType,
    HeatBattery as HeatBatteryInput, HeatPumpSourceType, HeatSource as HeatSourceInput,
    HeatSourceControlType, HeatSourceWetDetails, HotWaterSourceDetails,
    InfiltrationVentilation as InfiltrationVentilationInput, Input, InputForCalcHtcHlp,
    InternalGains as InternalGainsInput, InternalGainsDetails,
    OnSiteGeneration as OnSiteGenerationInput, PartyWallCavityType, PhotovoltaicInputs,
    PhotovoltaicSystem as PhotovoltaicSystemInput,
    PhotovoltaicSystemWithPanels as PhotovoltaicSystemWithPanelsInput,
    SpaceCoolSystem as SpaceCoolSystemInput, SpaceCoolSystemDetails,
    SpaceHeatSystem as SpaceHeatSystemInput, SpaceHeatSystemDetails, SystemReference,
    ThermalBridging as ThermalBridgingInput, ThermalBridgingDetails, UValueInput, VentilationLeaks,
    WasteWaterHeatRecovery, WasteWaterHeatRecoveryDetails, WaterHeatingEvent, WaterHeatingEvents,
    WaterPipework, WetEmitter, ZoneDictionary, ZoneInput, ZoneTemperatureControlBasis,
    MAIN_REFERENCE, PITCH_LIMIT_HORIZ_CEILING, PITCH_LIMIT_HORIZ_FLOOR,
};
use crate::input_dependency_resolvers::{
    build_preheated_water_source_dependency_graph, topological_sort_preheated_water_sources,
};
use crate::output::{
    Output, OutputCop, OutputCore, OutputHeatingCoolingSystem, OutputHotWaterSystems, OutputStatic,
    OutputSummary, OutputSummaryEnergySupply, OutputSummaryPeakElectricityConsumption,
    OutputZoneData,
};
use crate::output::{OutputEmitters, OutputMetadata};
use crate::simulation_time::{SimulationTimeIteration, SimulationTimeIterator};
use crate::statistics::percentile;
use crate::StringOrNumber;
use crate::{convert_profile_to_daily, HEM_VERSION};
use anyhow::{anyhow, bail};
use atomic_float::AtomicF64;
use chrono::{prelude::*, TimeDelta};
use erased_serde::__private::serde::Serializer;
use fsum::FSum;
use indexmap::IndexMap;
#[cfg(feature = "indicatif")]
use indicatif::ProgressIterator;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use parking_lot::{Mutex, RwLock};
use serde::{Deserialize, Serialize};
use serde_enum_str::{Deserialize_enum_str, Serialize_enum_str};
use smartstring::alias::String;
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::default::Default;
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
            *advanced_start,
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
            charge_calc_time,
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
                simulation_time_iterator,
                *start_day,
                *time_series_step,
                charge_level_vec,
                *temp_charge_cut,
                temp_charge_cut_delta,
                Some(external_conditions.clone()),
                external_sensor.clone(),
                Some(*charge_calc_time),
            )?)
            .into()
        }
        ControlDetails::OnOffCostMinimising {
            start_day,
            time_series_step,
            time_on_daily,
            schedule,
            ..
        } => Control::OnOffMinimisingTime(OnOffCostMinimisingTimeControl::new(
            reject_nulls(expand_numeric_schedule(schedule))?,
            simulation_time_iterator,
            *start_day,
            *time_series_step,
            *time_on_daily,
        )?)
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
        BuildingElementInput::PartyWall {
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
                area * u_value
            }
            BuildingElementInput::PartyWall {
                area,
                party_wall_cavity_data,
                ..
            } => {
                //Party wall calculations follow same approach as adjacent unconditioned space
                let r_se = calculate_cavity_resistance(
                    &PartyWallCavityType::from(*party_wall_cavity_data),
                    &party_wall_cavity_data.party_wall_lining_type(),
                    Some(thermal_resistance_construction), // Reported to DESNZ as suspect thermal_resistance_cavity should be passed in here instead
                )?;
                let u_value = 1.0 / (thermal_resistance_construction + r_se + r_si);
                area * u_value
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
        let wind_speed = external_conditions.wind_speed_annual()?;
        let wind_direction = external_conditions.wind_direction_annual()?;
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
            true.into(),
            simtime.iter().current_iteration(),
        )?;
        let total_vent_heat_loss = calc_vent_heat_transfer_coeff(zone.volume, air_changes_per_hour);

        // Calculate fabric heat loss and total floor area
        let total_fabric_heat_loss = FSum::with_all(
            zone.building_elements
                .values()
                .map(calc_heat_loss)
                .try_collect::<f64, Vec<f64>, anyhow::Error>()?,
        )
        .value();

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

    let total_htc = FSum::with_all(htc_map.values()).value();
    let total_floor_area = FSum::with_all(zone_area.values()).value();
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
    pub(crate) domestic_hot_water_demand:
        DomesticHotWaterDemand<HotWaterSource, HotWaterStorageTank>,
    r_v_arg: AtomicF64,
    pub(crate) ventilation: Arc<InfiltrationVentilation>,
    pub(crate) zones: IndexMap<Arc<str>, Arc<Zone>>,
    pub(crate) energy_supply_conn_unmet_demand_zone: IndexMap<String, Arc<EnergySupplyConnection>>,
    pub(crate) heat_system_name_for_zone: IndexMap<Arc<str>, Vec<Arc<str>>>,
    pub(crate) cool_system_name_for_zone: IndexMap<Arc<str>, Vec<Arc<str>>>,
    pub total_floor_area: f64,
    pub(crate) total_volume: f64,
    pub(crate) heat_sources_wet: IndexMap<String, WetHeatSource>,
    pub(crate) hot_water_sources: IndexMap<Arc<str>, HotWaterSource>,
    pub(crate) heat_sources_wet_with_buffer_tank: Vec<String>,
    pub(crate) space_heat_systems: IndexMap<Arc<str>, Arc<Mutex<SpaceHeatSystem>>>,
    pub(crate) space_cool_systems: IndexMap<Arc<str>, AirConditioning>,
    pub(crate) on_site_generation: IndexMap<String, PhotovoltaicSystem>,
    pub(crate) diverters: Vec<Arc<RwLock<PVDiverter>>>,
    required_vent_data: Option<RequiredVentData>,
    energy_supply_conn_names_for_hot_water_source: IndexMap<String, Vec<String>>,
    energy_supply_conn_names_for_heat_systems: IndexMap<Arc<str>, Arc<str>>,
    hotwatersource_name_for_heatsourcewet_service: IndexMap<Arc<str>, Arc<str>>,
    timestep_end_calcs: Arc<RwLock<Vec<HeatSystem>>>,
    initial_loop: AtomicBool,
    detailed_output_heating_cooling: bool,
    vent_adjust_min_control: Option<Arc<Control>>,
    vent_adjust_max_control: Option<Arc<Control>>,
    temp_internal_air_prev: Arc<AtomicF64>,
    smart_appliance_controls: IndexMap<String, Arc<SmartApplianceControl>>,
    input: Arc<Input>,
}

impl Corpus {
    pub fn from_inputs(
        input: Arc<Input>,
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

        let total_volume = FSum::with_all(input.zone.values().map(|zone| zone.volume)).value();

        let mut heat_system_name_for_zone: IndexMap<Arc<str>, Vec<Arc<str>>> = Default::default();
        let mut cool_system_name_for_zone: IndexMap<Arc<str>, Vec<Arc<str>>> = Default::default();

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

        let zones: IndexMap<Arc<str>, Arc<Zone>> = input
            .zone
            .iter()
            .map(
                |(zone_name, zone)| -> anyhow::Result<(Arc<str>, Arc<Zone>)> {
                    Ok((zone_name.to_string().into(), {
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
                },
            )
            .collect::<anyhow::Result<_>>()?;

        let energy_supply_conn_unmet_demand_zone = set_up_energy_supply_unmet_demand_zones(
            energy_supplies[UNMET_DEMAND_SUPPLY_NAME].clone(),
            &input.zone,
        );
        let total_floor_area = FSum::with_all(zones.values().map(|zone| zone.area())).value();

        // Internal gains is an ordered IndexMap. This is because load shifting behaviours
        // of appliance gains depend on other energy demand in the dwelling at any given time,
        // so depend on the order in which gains are considered by the engine.
        // See check_priority() in apply_appliance_gains_from_input()
        let mut internal_gains = internal_gains_from_input(
            &input.internal_gains,
            total_floor_area,
            simulation_time_iterator.as_ref(),
        )?;

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

        let timestep_end_calcs: Arc<RwLock<Vec<HeatSystem>>> = Default::default();

        // Register WWHRS objects for timestep_end calls
        for wwhrs in wwhrs.values() {
            timestep_end_calcs
                .write()
                .push(HeatSystem::WwhrsSystem(wwhrs.clone()));
        } // TODO review along with HeatSystem enum

        let mut heat_sources_wet_with_buffer_tank: Vec<String> = vec![];
        let mechanical_ventilations = infiltration_ventilation.mech_vents();

        let temp_internal_air_prev: Arc<AtomicF64> = Arc::new(AtomicF64::new(
            temp_internal_air_for_zones(&zones, total_volume),
        ));

        let mut heat_sources_wet: IndexMap<String, WetHeatSource> = input
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
                timestep_end_calcs
                    .write()
                    .push(HeatSystem::WetSystem(heat_source.clone()));
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
        let mut hot_water_source_name_for_service: IndexMap<Arc<str>, Arc<str>> =
            Default::default();
        let mut used_heat_source_names: HashSet<String> = Default::default();

        // Track pre-heat sources and WWHRS allocated to ensure single allocation only
        // This is required to avoid double-counting of the saving without significant additional
        // book-keeping code or a complete re-conceptualisation of the water heating calculation,
        // to handle an arrangement which is unlikely to occur in practice
        let mut cold_water_sources_already_allocated: HashSet<String> = Default::default();

        // processing pre-heated sources
        let mut pre_heated_water_sources: IndexMap<String, HotWaterStorageTank> =
            Default::default();

        let init_order = {
            let preheated_water_source_dependency_graph =
                build_preheated_water_source_dependency_graph(&input.pre_heated_water_source);
            topological_sort_preheated_water_sources(&preheated_water_source_dependency_graph)?
        };
        for source_name in init_order.iter() {
            let source_details = &input.pre_heated_water_source[source_name];
            let (heat_source, energy_conn_names, preheated_source_names_for_service) =
                hot_water_source_from_input(
                    source_name,
                    source_details,
                    &cold_water_sources,
                    &pre_heated_water_sources,
                    &mut heat_sources_wet,
                    &wwhrs,
                    &controls,
                    &mut energy_supplies,
                    &diverter_types,
                    &mut diverters,
                    shareable_fn(&temp_internal_air_prev),
                    simulation_time_iterator.clone().as_ref(),
                    external_conditions.clone(),
                    output_options.detailed_output_heating_cooling,
                    &mut used_heat_source_names,
                    &mut cold_water_sources_already_allocated,
                )?;
            energy_supply_conn_names_for_hot_water_source
                .insert(source_name.to_string().into(), energy_conn_names);
            if let HotWaterSource::PreHeated(source) = heat_source {
                pre_heated_water_sources.insert(source_name.to_string().into(), source);
            } else {
                bail!("Pre-heated water sources must be storage tanks");
            }
            hot_water_source_name_for_service.extend(
                preheated_source_names_for_service
                    .into_iter()
                    .map(|(x, y)| (x.as_str().into(), y.as_str().into())),
            );
        }

        let mut hot_water_sources: IndexMap<Arc<str>, HotWaterSource> = Default::default();
        for (name, data) in input.hot_water_source.iter() {
            let (hot_water_source, hw_cylinder_conn_names, source_names_for_service) =
                hot_water_source_from_input(
                    name,
                    data,
                    &cold_water_sources,
                    &pre_heated_water_sources,
                    &mut heat_sources_wet,
                    &wwhrs,
                    &controls,
                    &mut energy_supplies,
                    &diverter_types,
                    &mut diverters,
                    shareable_fn(&temp_internal_air_prev),
                    simulation_time_iterator.clone().as_ref(),
                    external_conditions.clone(),
                    output_options.detailed_output_heating_cooling,
                    &mut used_heat_source_names,
                    &mut cold_water_sources_already_allocated,
                )?;
            hot_water_sources.insert(name.to_string().into(), hot_water_source);
            energy_supply_conn_names_for_hot_water_source
                .insert(name.to_string().into(), hw_cylinder_conn_names);
            hot_water_source_name_for_service.extend(
                source_names_for_service
                    .into_iter()
                    .map(|(x, y)| (x.as_str().into(), y.as_str().into())),
            );
        }

        let domestic_hot_water_demand = DomesticHotWaterDemand::new(
            &input.hot_water_demand.shower,
            &input.hot_water_demand.bath,
            &input.hot_water_demand.other_water_use,
            &input.hot_water_demand.water_distribution,
            &cold_water_sources,
            &wwhrs,
            &energy_supplies,
            event_schedules,
            hot_water_sources
                .iter()
                .map(|(k, v)| (k.clone(), v.clone()))
                .collect(),
            pre_heated_water_sources.clone(),
        )?;

        let mut heat_system_names_requiring_overvent: Vec<Arc<str>> = Default::default();

        let (space_heat_systems, energy_supply_conn_names_for_heat_systems) = input
            .space_heat_system
            .as_ref()
            .map(|system| {
                anyhow::Ok(space_heat_systems_from_input(
                    system,
                    &controls,
                    &mut energy_supplies,
                    simulation_time_iterator.as_ref(),
                    &heat_sources_wet,
                    &mut heat_system_names_requiring_overvent,
                    &heat_system_name_for_zone,
                    &zones,
                    &heat_sources_wet_with_buffer_tank
                        .iter()
                        .cloned()
                        .collect_vec(),
                    external_conditions.clone(),
                    output_options.detailed_output_heating_cooling,
                    temp_internal_air_prev.load(Ordering::SeqCst),
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
                        .map(|s| s.as_ref())
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
            heat_sources_wet,
            hot_water_sources,
            heat_sources_wet_with_buffer_tank,
            space_heat_systems,
            space_cool_systems,
            on_site_generation,
            diverters,
            required_vent_data,
            energy_supply_conn_names_for_hot_water_source,
            energy_supply_conn_names_for_heat_systems,
            hotwatersource_name_for_heatsourcewet_service: hot_water_source_name_for_service,
            timestep_end_calcs,
            initial_loop: AtomicBool::new(false),
            detailed_output_heating_cooling: output_options.detailed_output_heating_cooling,
            vent_adjust_min_control,
            vent_adjust_max_control,
            temp_internal_air_prev,
            smart_appliance_controls,
            input,
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
    pub(crate) fn calc_hlff(&self) -> f64 {
        let total_heat_loss_area = self
            .zones
            .values()
            .map(|zone| zone.total_heat_loss_area())
            .sum::<f64>();

        total_heat_loss_area / self.total_floor_area
    }

    fn update_temp_internal_air(&self) {
        self.temp_internal_air_prev.store(
            temp_internal_air_for_zones(&self.zones, self.total_volume),
            Ordering::SeqCst,
        );
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
        self.temp_internal_air_prev.load(Ordering::SeqCst)
    }

    /// Calculate the losses in the buffer tank
    fn calc_internal_gains_buffer_tank(&self) -> f64 {
        self.heat_sources_wet_with_buffer_tank
            .iter()
            .map(
                |heat_source_name| match self.heat_sources_wet.get(heat_source_name) {
                    Some(heat_source) => match heat_source {
                        WetHeatSource::HeatPump(heat_pump) => heat_pump.lock().buffer_int_gains(),
                        _ => unreachable!(),
                    },
                    None => 0.,
                },
            )
            .sum::<f64>()
    }

    /// Calculate the losses/gains in the MVHR ductwork, in Watts
    fn calc_internal_gains_ductwork(&self, simulation_time: SimulationTimeIteration) -> f64 {
        let temp_outdoor_air = self.external_conditions.air_temp(&simulation_time);
        let temp_indoor_air = self.temp_internal_air_prev_timestep();

        self.ventilation
            .calc_internal_gains_ductwork(temp_outdoor_air, temp_indoor_air)
    }

    fn space_heat_internal_gains_for_zone(
        &self,
        zone: &Zone,
        gains_internal_dhw_on_site_generation: f64,
        gains_internal_hb: f64,
        internal_gains_ductwork_per_m3: f64,
        gains_internal_buffer_tank: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // Initialise to dhw internal gains split proportionally to zone floor area
        let mut gains_internal_zone = (gains_internal_buffer_tank
            + gains_internal_dhw_on_site_generation
            + gains_internal_hb)
            * zone.area()
            / self.total_floor_area;

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
    ) -> anyhow::Result<HeatCoolSystemsForZone> {
        let SetpointsAndConvectiveFractions {
            temp_setpnt_heat: temp_setpnt_heat_system,
            temp_setpnt_cool: temp_setpnt_cool_system,
            frac_convective_heat: frac_convective_heat_system,
            frac_convective_cool: frac_convective_cool_system,
        } = self.setpoints_and_convective_fractions(
            &self.heat_system_name_for_zone[z_name],
            &self.cool_system_name_for_zone[z_name],
            simtime,
        )?;

        // Sort heating and cooling systems by setpoint (highest first for
        // heating, lowest first for cooling)
        // In the event of two systems having the same setpoint, the one
        // listed first by the user takes priority
        let h_name_list_sorted: Vec<Arc<str>> = temp_setpnt_heat_system
            .iter()
            .sorted_by(|a, b| OrderedFloat(*a.1).cmp(&OrderedFloat(*b.1)))
            .rev()
            .map(|x| x.0.to_owned())
            .collect();
        let c_name_list_sorted: Vec<Arc<str>> = temp_setpnt_cool_system
            .iter()
            .sorted_by(|a, b| OrderedFloat(*a.1).cmp(&OrderedFloat(*b.1)))
            .map(|x| x.0.to_owned())
            .collect();

        Ok(HeatCoolSystemsForZone {
            h_sorted_names: h_name_list_sorted,
            c_sorted_names: c_name_list_sorted,
            setpoints_and_convective_fractions: SetpointsAndConvectiveFractions {
                temp_setpnt_heat: temp_setpnt_heat_system,
                temp_setpnt_cool: temp_setpnt_cool_system,
                frac_convective_heat: frac_convective_heat_system,
                frac_convective_cool: frac_convective_cool_system,
            },
        })
    }

    fn setpoints_and_convective_fractions(
        &self,
        h_name_list: &Vec<Arc<str>>,
        c_name_list: &Vec<Arc<str>>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<SetpointsAndConvectiveFractions> {
        let mut frac_convective_heat: IndexMap<Arc<str>, f64> = Default::default();
        let mut frac_convective_cool: IndexMap<Arc<str>, f64> = Default::default();
        let mut temp_setpnt_heat: IndexMap<Arc<str>, f64> = Default::default();
        let mut temp_setpnt_cool: IndexMap<Arc<str>, f64> = Default::default();

        for h_name in h_name_list {
            match h_name {
                h_name if h_name.as_ref() == "" => {
                    frac_convective_heat.insert(h_name.clone(), 1.0);
                    temp_setpnt_heat.insert(h_name.clone(), temp_setpnt_heat_none());
                }
                h_name => {
                    let space_heat_system = self.space_heat_systems.get(h_name).unwrap().lock();
                    frac_convective_heat
                        .insert(h_name.clone(), space_heat_system.frac_convective(simtime));
                    temp_setpnt_heat.insert(
                        h_name.clone(),
                        space_heat_system
                            .temp_setpnt(simtime)
                            .unwrap_or_else(temp_setpnt_heat_none),
                    );
                }
            }
        }

        for c_name in c_name_list {
            match c_name {
                c_name if c_name.as_ref() == "" => {
                    frac_convective_cool.insert(c_name.clone(), 1.0);
                    temp_setpnt_cool.insert(c_name.clone(), temp_setpnt_cool_none());
                }
                c_name => {
                    let space_cool_system = self.space_cool_systems.get(c_name).unwrap();
                    frac_convective_cool
                        .insert(c_name.clone(), space_cool_system.frac_convective());
                    temp_setpnt_cool.insert(
                        c_name.clone(),
                        space_cool_system
                            .temp_setpnt(&simtime)
                            .unwrap_or_else(temp_setpnt_cool_none),
                    );
                }
            }
        }

        Ok(SetpointsAndConvectiveFractions {
            temp_setpnt_heat,
            temp_setpnt_cool,
            frac_convective_heat,
            frac_convective_cool,
        })
    }

    fn gains_heat_cool(
        &self,
        delta_t_h: f64,
        hc_output_convective: &IndexMap<Arc<str>, f64>,
        hc_output_radiative: &IndexMap<Arc<str>, f64>,
    ) -> (f64, f64) {
        let gains_heat_cool_convective = FSum::with_all(hc_output_convective.values()).value()
            * WATTS_PER_KILOWATT as f64
            / delta_t_h;
        let gains_heat_cool_radiative = FSum::with_all(hc_output_radiative.values()).value()
            * WATTS_PER_KILOWATT as f64
            / delta_t_h;

        (gains_heat_cool_convective, gains_heat_cool_radiative)
    }

    /// Calculate the incoming air changes per hour
    /// initial_p_z_ref_guess is used for calculation in first timestep.
    /// Later timesteps use the previous timesteps p_z_ref of max and min ACH ,respective to calc.
    fn calc_air_changes_per_hour(
        &self,
        wind_speed: f64,
        wind_direction: Orientation360,
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
        h_name_list_sorted_zone: &HashMap<&str, Vec<Arc<str>>>,
        c_name_list_sorted_zone: &HashMap<&str, Vec<Arc<str>>>,
        frac_convective_heat_zone_system: &HashMap<&str, IndexMap<Arc<str>, f64>>,
        frac_convective_cool_zone_system: &HashMap<&str, IndexMap<Arc<str>, f64>>,
        z_name: &str,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<HeatCoolOutputs> {
        let h_output_min: IndexMap<Arc<str>, f64> = h_name_list_sorted_zone[z_name]
            .iter()
            .filter(|h_name| h_name.as_ref() != "") // we need to exclude the empty string as it stands for None (yes, we're stringly typing this)
            .map(|h_name| -> anyhow::Result<(Arc<str>, f64)> {
                Ok((
                    h_name.clone(),
                    self.space_heat_systems
                        .get(h_name.as_ref())
                        .unwrap()
                        .lock()
                        .energy_output_min(simulation_time_iteration)?,
                ))
            })
            .try_collect()?;
        let c_output_min = c_name_list_sorted_zone[z_name]
            .iter()
            .filter(|c_name| c_name.as_ref() != "") // we need to exclude the empty string as it stands for None (yes, we're stringly typing this)
            .map(|c_name| {
                (
                    c_name.clone(),
                    self.space_cool_systems
                        .get(c_name.as_ref())
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
                    hc_output_min[hc_name.as_ref()] * frac_convective_system[hc_name],
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
    /// * `gains_internal_dhw_on_site_generation` - internal gains from hot water system and on site generation for this timestep, in W
    /// * `gains_internal_hb` - internal gains from central heat battery systems for this timestep, in W
    fn calc_space_heating(
        &self,
        delta_t_h: f64,
        gains_internal_dhw_on_site_generation: f64,
        gains_internal_hb: f64,
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

        let mut gains_internal_zone: HashMap<Arc<str>, f64> = Default::default();
        let mut gains_solar_zone: HashMap<Arc<str>, f64> = Default::default();
        let mut h_name_list_sorted_zone: HashMap<&str, Vec<Arc<str>>> = Default::default();
        let mut c_name_list_sorted_zone: HashMap<&str, Vec<Arc<str>>> = Default::default();
        let mut temp_setpnt_heat_zone_system: HashMap<&str, IndexMap<Arc<str>, f64>> =
            Default::default();
        let mut temp_setpnt_cool_zone_system: HashMap<&str, IndexMap<Arc<str>, f64>> =
            Default::default();
        let mut frac_convective_heat_zone_system: HashMap<&str, IndexMap<Arc<str>, f64>> =
            Default::default();
        let mut frac_convective_cool_zone_system: HashMap<&str, IndexMap<Arc<str>, f64>> =
            Default::default();
        let mut ach_cooling_zone: HashMap<&str, f64> = Default::default();
        let mut ach_to_trigger_heating_zone: HashMap<&str, Option<f64>> = Default::default();
        let mut internal_air_temp: HashMap<Arc<str>, f64> = Default::default();
        let mut operative_temp: HashMap<Arc<str>, f64> = Default::default();
        let mut space_heat_demand_zone: HashMap<Arc<str>, f64> = Default::default();
        let mut space_cool_demand_zone: HashMap<Arc<str>, f64> = Default::default();
        let mut space_heat_provided_system: HashMap<Arc<str>, f64> = Default::default();
        let mut space_cool_provided_system: HashMap<Arc<str>, f64> = Default::default();
        let mut heat_balance_map: HashMap<Arc<str>, Option<HeatBalance>> = Default::default();

        // Average supply temperature
        let avg_air_supply_temp = self.external_conditions.air_temp(&simtime);

        for (z_name, zone) in self.zones.iter() {
            let z_name = z_name.as_ref();
            // Calculate internal and solar gains
            gains_internal_zone.insert(
                z_name.into(),
                self.space_heat_internal_gains_for_zone(
                    zone,
                    gains_internal_dhw_on_site_generation,
                    gains_internal_hb,
                    internal_gains_ductwork_per_m3,
                    internal_gains_buffer_tank,
                    simtime,
                )?,
            );
            gains_solar_zone.insert(z_name.into(), zone.gains_solar(simtime));

            // Get heating and cooling characteristics for the current zone
            let HeatCoolSystemsForZone {
                h_sorted_names: h_name_list_sorted_zone_current,
                c_sorted_names: c_name_list_sorted_zone_current,
                setpoints_and_convective_fractions:
                    SetpointsAndConvectiveFractions {
                        temp_setpnt_heat: temp_setpnt_heat_zone_system_current,
                        temp_setpnt_cool: temp_setpnt_cool_zone_system_current,
                        frac_convective_heat: frac_convective_heat_zone_system_current,
                        frac_convective_cool: frac_convective_cool_zone_system_current,
                    },
            } = self.heat_cool_systems_for_zone(z_name, simtime)?;

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
                let z_name = z_name.as_ref();
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
            let z_name = z_name.as_ref();
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
                .map(|h_name| (h_name.as_ref(), 0.0))
                .collect::<IndexMap<_, _>>();
            let mut space_cool_demand_zone_system = c_name_list_sorted_zone[z_name]
                .iter()
                .map(|c_name| (c_name.as_ref(), 0.0))
                .collect::<IndexMap<_, _>>();
            let mut space_heat_provided_zone_system = h_name_list_sorted_zone[z_name]
                .iter()
                .map(|h_name| (h_name.as_ref(), 0.0))
                .collect::<IndexMap<_, _>>();
            let mut space_cool_provided_zone_system = c_name_list_sorted_zone[z_name]
                .iter()
                .map(|c_name| (c_name.as_ref(), 0.0))
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
                let h_name = &h_name_list_sorted_zone[z_name][h_idx].as_ref();
                let c_name = &c_name_list_sorted_zone[z_name][c_idx].as_ref();
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
                let h_name = h_name.as_ref();
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
                let c_name = c_name.as_ref();
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
            let hc_output_convective_total = FSum::with_all(hc_output_convective.values()).value();
            let hc_output_radiative_total = FSum::with_all(hc_output_radiative.values()).value();

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
                    .or_insert(0.0) += space_heat_provided_zone_system[h_name.as_ref()];
            }
            for c_name in c_name_list_sorted_zone[z_name].iter() {
                *space_cool_provided_system
                    .entry(c_name.to_owned())
                    .or_insert(0.0) += space_cool_provided_zone_system[c_name.as_ref()];
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
        hc_name_list_sorted: &[Arc<str>],
        space_heat_cool_systems: SpaceHeatCoolSystems,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Option<Arc<str>>> {
        let mut hc_name_highest_req = Default::default();
        for hc_name in hc_name_list_sorted {
            if !hc_name.is_empty()
                && space_heat_cool_systems
                    .in_required_period_for_name(hc_name, simtime)
                    .unwrap_or(false)
            {
                hc_name_highest_req = Some(hc_name.clone());
                break;
            }
        }

        Ok(hc_name_highest_req)
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
        temp_setpnt_heat_system: &IndexMap<Arc<str>, f64>,
        temp_setpnt_cool_system: &IndexMap<Arc<str>, f64>,
        frac_convective_heat_system: &IndexMap<Arc<str>, f64>,
        frac_convective_cool_system: &IndexMap<Arc<str>, f64>,
        h_name_list_sorted: &[Arc<str>],
        c_name_list_sorted: &[Arc<str>],
        space_heat_demand: f64,
        space_cool_demand: f64,
        hc_output_convective: &IndexMap<Arc<str>, f64>,
        hc_output_radiative: &IndexMap<Arc<str>, f64>,
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
        )?;
        let c_name_highest_req = self.highest_priority_required_system(
            c_name_list_sorted,
            SpaceHeatCoolSystems::Cool(&self.space_cool_systems),
            simtime,
        )?;

        let gains_heat = FSum::with_all(h_name_list_sorted.iter().map(|h_name| {
            hc_output_convective[h_name.as_ref()] + hc_output_radiative[h_name.as_ref()]
        }))
        .value();
        let gains_cool = FSum::with_all(c_name_list_sorted.iter().map(|c_name| {
            hc_output_convective[c_name.as_ref()] + hc_output_radiative[c_name.as_ref()]
        }))
        .value();
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

    pub fn run(&self) -> anyhow::Result<Output> {
        let simulation_time = self.simulation_time.as_ref().to_owned();
        let vec_capacity = || Vec::with_capacity(simulation_time.total_steps());

        let mut timestep_array = vec_capacity();
        let mut gains_internal_dict: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        let mut gains_solar_dict: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        let mut operative_temp_dict: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        let mut internal_air_temp_dict: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        let mut space_heat_demand_dict: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        let mut space_cool_demand_dict: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        let mut space_heat_provided_dict: IndexMap<Option<Arc<str>>, Vec<f64>> = Default::default();
        let mut space_cool_provided_dict: IndexMap<Option<Arc<str>>, Vec<f64>> = Default::default();
        let mut zone_list: Vec<Arc<str>> = Default::default();
        let mut heat_balance_all_dict: HeatBalanceAllResults = IndexMap::from([
            (HeatBalanceFieldName::AirNode, Default::default()),
            (HeatBalanceFieldName::InternalBoundary, Default::default()),
            (HeatBalanceFieldName::ExternalBoundary, Default::default()),
        ]);
        let mut heat_source_wet_results_dict: IndexMap<Arc<str>, ResultsPerTimestep> =
            Default::default();
        let mut heat_source_wet_results_annual_dict: IndexMap<Arc<str>, ResultsAnnual> =
            Default::default();
        let mut emitters_output_dict: IndexMap<Arc<str>, Vec<OutputEmitters>> = Default::default();
        let mut esh_output_dict: IndexMap<Arc<str>, Vec<StorageHeaterDetailedResult>> =
            Default::default();
        let mut vent_output_list: Vec<VentilationDetailedResult> = Default::default();
        let mut hot_water_source_results_summary: IndexMap<Arc<str>, Vec<Vec<StringOrNumber>>> =
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
                let h_name = match h_name.as_ref() {
                    "" => None,
                    _ => Some(h_name.clone()),
                };
                space_heat_provided_dict.insert(h_name, vec_capacity());
            }
        }

        for z_c_names in self.cool_system_name_for_zone.values() {
            for c_name in z_c_names {
                let c_name = match c_name.as_ref() {
                    "" => None,
                    _ => Some(c_name.clone()),
                };
                space_cool_provided_dict.insert(c_name, vec_capacity());
            }
        }

        let mut list_hot_water_source_names_incl_electric_showers: Vec<Arc<str>> = self
            .hot_water_sources
            .keys()
            .map(|x| x.to_string().into())
            .collect();
        list_hot_water_source_names_incl_electric_showers
            .push(ELECTRIC_SHOWERS_HWS_NAME.to_string().into());

        let mut hot_water_demand: IndexMap<Arc<str>, Vec<f64>> =
            list_hot_water_source_names_incl_electric_showers
                .iter()
                .map(|hws_name| (hws_name.clone(), vec_capacity()))
                .collect();
        let mut hot_water_energy_demand_at_tapping_points: IndexMap<Arc<str>, Vec<f64>> =
            list_hot_water_source_names_incl_electric_showers
                .iter()
                .map(|hws_name| (hws_name.clone(), vec_capacity()))
                .collect();
        let mut hot_water_energy_demand_at_hot_water_source: IndexMap<Arc<str>, Vec<f64>> = self
            .hot_water_sources
            .keys()
            .map(|hws_name| (hws_name.clone(), vec_capacity()))
            .collect();
        let mut hot_water_energy_output: IndexMap<Arc<str>, Vec<f64>> = self
            .hot_water_sources
            .keys()
            .map(|hws_name| (hws_name.clone(), vec_capacity()))
            .collect();
        let mut hot_water_duration: IndexMap<Arc<str>, Vec<f64>> =
            list_hot_water_source_names_incl_electric_showers
                .iter()
                .map(|hws_name| (hws_name.clone(), vec_capacity()))
                .collect();
        let mut hot_water_no_events: IndexMap<Arc<str>, Vec<f64>> =
            list_hot_water_source_names_incl_electric_showers
                .iter()
                .map(|hws_name| (hws_name.clone(), vec_capacity()))
                .collect();
        let mut hot_water_pipework: IndexMap<Arc<str>, Vec<f64>> =
            list_hot_water_source_names_incl_electric_showers
                .iter()
                .map(|hws_name| (hws_name.clone(), vec_capacity()))
                .collect();
        let mut hot_water_primary_pipework: IndexMap<Arc<str>, Vec<f64>> = self
            .hot_water_sources
            .keys()
            .map(|hws_name| (hws_name.clone(), vec_capacity()))
            .collect();
        let mut hot_water_storage_losses: IndexMap<Arc<str>, Vec<f64>> = self
            .hot_water_sources
            .keys()
            .map(|hws_name| (hws_name.clone(), vec_capacity()))
            .collect();
        let mut ductwork_gains_list = vec_capacity();
        self.initial_loop.store(true, Ordering::SeqCst);
        let mut internal_pressure_window: HashMap<ReportingFlag, f64> = Default::default();

        let delta_t_h = simulation_time.step_in_hours();

        // Loop over each timestep
        #[cfg(feature = "indicatif")]
        let simulation_time_iter = simulation_time.progress();
        #[cfg(not(feature = "indicatif"))]
        let simulation_time_iter = simulation_time;

        for t_it in simulation_time_iter {
            timestep_array.push(t_it.time);
            self.update_temp_internal_air();

            let WaterHeatingCalculation {
                hw_demand_vol,
                hw_duration,
                no_events,
                hw_energy_demand_at_tapping_points,
                hw_energy_demand_at_hot_water_source,
                hw_energy_output,
                pw_losses_total,
                primary_pw_losses,
                storage_losses,
                gains_internal_dhw,
            } = self.domestic_hot_water_demand.calc_water_heating(
                t_it,
                self.temp_internal_air_prev_timestep(),
                self.external_conditions.air_temp(&t_it),
            )?;

            let mut gains_internal_dhw_on_site_generation =
                FSum::with_all(gains_internal_dhw.values()).value();

            let mut gains_internal_hb = 0.;
            // Adding heat battery losses to internal gains
            for heat_source_wet in self.heat_sources_wet.values() {
                if let WetHeatSource::HeatBattery(battery) = heat_source_wet {
                    gains_internal_hb +=
                        battery.get_battery_losses() * WATTS_PER_KILOWATT as f64 / t_it.timestep
                }
            }

            // loop through on-site energy generation
            for pv in self.on_site_generation.values() {
                // Get energy produced for the current timestep
                let (_energy_produced, energy_lost) = pv.produce_energy(t_it)?;
                // Add the energy lost figure to the internal gains if it is considered inside the building
                if pv.inverter_is_inside() {
                    gains_internal_dhw_on_site_generation +=
                        energy_lost * WATTS_PER_KILOWATT as f64 / t_it.timestep;
                }
            }

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
                delta_t_h,
                gains_internal_dhw_on_site_generation,
                gains_internal_hb,
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
                    .get_mut(&z_name)
                    .unwrap()
                    .push(gains_internal);
            }

            for (z_name, gains_solar) in gains_solar_zone {
                gains_solar_dict.get_mut(&z_name).unwrap().push(gains_solar);
            }

            for (z_name, temp) in operative_temp {
                operative_temp_dict.get_mut(&z_name).unwrap().push(temp);
            }

            for (z_name, temp) in internal_air_temp {
                internal_air_temp_dict.get_mut(&z_name).unwrap().push(temp);
            }

            for (z_name, demand) in space_heat_demand_zone {
                space_heat_demand_dict
                    .get_mut(&z_name)
                    .unwrap()
                    .push(demand);
            }

            for (z_name, demand) in space_cool_demand_zone {
                space_cool_demand_dict
                    .get_mut(&z_name)
                    .unwrap()
                    .push(demand);
            }

            for (h_name, output) in space_heat_provided {
                space_heat_provided_dict
                    .get_mut(&Some(h_name.clone()))
                    .unwrap()
                    .push(output);
            }

            for (c_name, output) in space_cool_provided {
                space_cool_provided_dict
                    .get_mut(&Some(c_name.clone()))
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

            fn append_dict_vals_to_dict_lists(
                destination: &mut IndexMap<Arc<str>, Vec<f64>>,
                origin: &IndexMap<Arc<str>, f64>,
            ) {
                for (key, value) in origin {
                    destination[key].push(*value);
                }
            }

            append_dict_vals_to_dict_lists(&mut hot_water_demand, &hw_demand_vol);
            append_dict_vals_to_dict_lists(
                &mut hot_water_energy_demand_at_tapping_points,
                &hw_energy_demand_at_tapping_points,
            );
            append_dict_vals_to_dict_lists(
                &mut hot_water_energy_demand_at_hot_water_source,
                &hw_energy_demand_at_hot_water_source,
            );
            append_dict_vals_to_dict_lists(&mut hot_water_energy_output, &hw_energy_output);
            append_dict_vals_to_dict_lists(&mut hot_water_duration, &hw_duration);
            let no_events = no_events
                .iter()
                .map(|(k, v)| (k.clone(), *v as f64))
                .collect();
            append_dict_vals_to_dict_lists(&mut hot_water_no_events, &no_events);
            append_dict_vals_to_dict_lists(&mut hot_water_pipework, &pw_losses_total);
            append_dict_vals_to_dict_lists(&mut hot_water_primary_pipework, &primary_pw_losses);
            append_dict_vals_to_dict_lists(&mut hot_water_storage_losses, &storage_losses);

            ductwork_gains_list.push(ductwork_gains);

            for supply in self.energy_supplies.values() {
                supply.read().calc_energy_import_export_betafactor(t_it)?;
                supply
                    .read()
                    .calc_energy_import_from_grid_to_battery(t_it)?;
                supply.read().timestep_end();
            }

            for diverter in &self.diverters {
                diverter.write().timestep_end();
            }

            for (_, control) in &self.smart_appliance_controls {
                control.update_demand_buffer(t_it);
            }
        }

        // Report detailed outputs from heat source wet objects, if requested and available
        // to take same type (IndexMap<String, Vec<ResultParamValue> or IndexMap<String, Vec<f64>)
        let hot_water_energy_output_as_result_param_value: IndexMap<
            Arc<str>,
            Vec<ResultParamValue>,
        > = hot_water_energy_output
            .iter()
            .map(|(key, value)| {
                (
                    key.clone(),
                    value
                        .iter()
                        .map(|element: &f64| element.into())
                        .collect_vec(),
                )
            })
            .collect();
        if self.detailed_output_heating_cooling {
            for (name, heat_source_wet) in self.heat_sources_wet.iter() {
                let name: Arc<str> = name.as_str().into();
                if let Some((results, results_annual)) = heat_source_wet.output_detailed_results(
                    &hot_water_energy_output_as_result_param_value,
                    &self.hotwatersource_name_for_heatsourcewet_service,
                ) {
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
                                hot_water_source_results_summary
                                    .insert(name.to_string().into(), hot_water_source_output);
                            }
                        }
                        HotWaterStorageTank::SmartHotWaterTank(smart_storage_tank) => {
                            if let Some(hot_water_source_output) =
                                smart_storage_tank.read().output_results()
                            {
                                hot_water_source_results_summary
                                    .insert(name.to_string().into(), hot_water_source_output);
                            }
                        }
                    }
                }
            }
        }

        // Return results from all energy supplies
        let mut results_totals: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        let mut results_end_user: IndexMap<Arc<str>, IndexMap<Arc<str>, Vec<f64>>> =
            Default::default();
        let mut energy_import: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        let mut energy_export: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        let mut grid_to_consumption: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        let mut generation_to_grid: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        let mut energy_generated_consumed: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        let mut energy_to_storage: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        let mut energy_from_storage: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        let mut storage_from_grid: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        let mut battery_state_of_charge: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        let mut energy_diverted: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        let mut beta_factor: IndexMap<Arc<str>, Vec<f64>> = Default::default();
        for (name, supply) in self
            .energy_supplies
            .iter()
            .map(|(name, supply)| (name.to_owned(), Arc::clone(supply)))
        {
            let name: Arc<str> = name.to_string().into();
            let supply = supply.read();
            results_totals.insert(name.clone(), supply.results_total());
            results_end_user.insert(name.clone(), supply.results_by_end_user().to_owned());
            energy_import.insert(name.clone(), supply.get_energy_import().to_owned());
            energy_export.insert(name.clone(), supply.get_energy_export().to_owned());
            grid_to_consumption.insert(name.clone(), supply.get_grid_to_consumption().to_owned());
            generation_to_grid.insert(name.clone(), supply.get_energy_export().to_owned());
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
            beta_factor.insert(name, supply.get_beta_factor().to_owned());
        }

        let dhw_cop_dict = self.heat_cool_cop(
            &hot_water_energy_output
                .into_iter()
                .map(|(k, v)| (Some(k.clone()), v))
                .collect(),
            &results_end_user,
            self.energy_supply_conn_names_for_hot_water_source
                .iter()
                .map(|(k, v)| {
                    (
                        k.to_string().into(),
                        v.iter().map(|x| x.as_str().into()).collect(),
                    )
                })
                .collect(),
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

        let output_zone_data = OutputZoneData {
            internal_gains: gains_internal_dict,
            solar_gains: gains_solar_dict,
            operative_temp: operative_temp_dict,
            internal_air_temp: internal_air_temp_dict,
            space_heat_demand: space_heat_demand_dict,
            space_cool_demand: space_cool_demand_dict,
        };

        let output_core = OutputCore {
            timestep_array,
            results_totals,
            results_end_user,
            energy_import,
            energy_export,
            grid_to_consumption,
            generation_to_grid,
            energy_generated_consumed,
            energy_to_storage,
            energy_from_storage,
            storage_from_grid,
            battery_state_of_charge,
            energy_diverted,
            beta_factor,
            zone_data: output_zone_data,
            zone_list,
            heating_cooling_system: OutputHeatingCoolingSystem {
                heating_system_output: space_heat_provided_dict,
                cooling_system_output: space_cool_provided_dict,
            },
            hot_water_systems: OutputHotWaterSystems {
                demand: hot_water_demand,
                energy_demand_at_hot_water_source: hot_water_energy_demand_at_hot_water_source,
                energy_demand_at_tapping_points: hot_water_energy_demand_at_tapping_points,
                duration: hot_water_duration,
                events_count: hot_water_no_events,
                losses_pipework: hot_water_pipework,
                losses_primary_pipework: hot_water_primary_pipework,
                losses_storage: hot_water_storage_losses,
            },
            cop: OutputCop {
                space_heating_system: heat_cop_dict,
                space_cooling_system: cool_cop_dict,
                hot_water_system: dhw_cop_dict,
            },
            ductwork_gains: ductwork_gains_list,
            heat_balance_all: heat_balance_all_dict
                .into_iter()
                .map(|(k, v)| (Arc::<str>::from(k), v))
                .collect(), // TODO (from Python) could be output object too fixed keys...
            heat_source_wet_results: heat_source_wet_results_dict,
            heat_source_wet_results_annual: heat_source_wet_results_annual_dict,
            hot_water_source_results_summary,
            emitters: emitters_output_dict
                .into_iter()
                .map(|(k, v)| (k, v.into_iter().enumerate().collect()))
                .collect(),
            electric_storage_heaters: esh_output_dict
                .into_iter()
                .map(|(k, v)| {
                    (
                        k,
                        v.into_iter()
                            .enumerate()
                            .map(|(k, v)| (k, v.into()))
                            .collect(),
                    )
                })
                .collect(),
            ventilation: vent_output_list.into_iter().map(Into::into).collect(),
        };

        Ok(Output {
            static_: self.calculate_output_static()?,
            summary: self.calculate_output_summary(&output_core),
            core: output_core,
            metadata: OutputMetadata {
                hem_core_version: HEM_VERSION.into(),
            },
        })
    }

    /// Calculate overall CoP over calculation period for each heating and cooling system
    fn heat_cool_cop(
        &self,
        energy_provided: &IndexMap<Option<Arc<str>>, Vec<f64>>,
        results_end_user: &IndexMap<Arc<str>, IndexMap<Arc<str>, Vec<f64>>>,
        energy_supply_conn_name_for_space_hc_system: IndexMap<Arc<str>, Vec<Arc<str>>>,
    ) -> IndexMap<Arc<str>, NumberOrDivisionByZero> {
        let mut hc_output_overall: IndexMap<Arc<str>, f64> = Default::default();
        let mut hc_input_overall: IndexMap<Arc<str>, f64> = Default::default();
        let mut cop_dict: IndexMap<Arc<str>, NumberOrDivisionByZero> = Default::default();
        // TODO review hc_name type, we're using "" instead of None
        for (hc_name, hc_output) in energy_provided {
            let hc_name = match hc_name.as_ref() {
                Some(hc_name) if !hc_name.is_empty() => hc_name,
                _ => continue,
            };
            hc_output_overall.insert(hc_name.clone(), FSum::with_all(hc_output).value().abs());
            hc_input_overall.insert(hc_name.clone(), 0.);
            let energy_supply_conn_names =
                energy_supply_conn_name_for_space_hc_system[hc_name].clone();
            for (fuel_name, fuel_summary) in results_end_user {
                if fuel_name.as_ref() == UNMET_DEMAND_SUPPLY_NAME
                    || fuel_name.as_ref() == ENERGY_FROM_ENVIRONMENT_SUPPLY_NAME
                {
                    continue;
                }
                for (conn_name, energy_cons) in fuel_summary {
                    if energy_supply_conn_names.contains(conn_name) {
                        *hc_input_overall.get_mut(hc_name).unwrap() +=
                            FSum::with_all(energy_cons).value();
                    }
                }
            }

            cop_dict.insert(
                hc_name.clone(),
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

    fn calculate_output_static(&self) -> anyhow::Result<OutputStatic> {
        let HtcHlpCalculation {
            total_htc: heat_trans_coeff,
            total_hlp: heat_loss_param,
            ..
        } = calc_htc_hlp(self.input.as_ref())?;

        let heat_capacity_param = self.calc_hcp();
        let heat_loss_form_factor = self.calc_hlff();

        Ok(OutputStatic {
            heat_transfer_coefficient: heat_trans_coeff,
            heat_loss_param,
            heat_capacity_param,
            heat_loss_form_factor,
            temperature_air_internal: self.input.temp_internal_air_static_calcs(),
            temperature_air_external: self.external_conditions.air_temp_annual_daily_average_min(),
        })
    }

    fn calculate_output_summary(&self, output_core: &OutputCore) -> OutputSummary {
        let mut energy_supply_stats: IndexMap<Arc<str>, OutputSummaryEnergySupply> =
            Default::default();
        for (key, result) in &output_core.results_end_user {
            let mut total_generated = 0.;
            let mut total_consumed = 0.;
            for result_value in result.values() {
                let sum_arr = FSum::with_all(result_value).value();
                if sum_arr < 0. {
                    total_generated += sum_arr.abs();
                } else {
                    total_consumed += sum_arr;
                }
            }

            let grid_to_consumption = FSum::with_all(&output_core.grid_to_consumption[key]).value();
            let generation_to_grid = FSum::with_all(&output_core.generation_to_grid[key])
                .value()
                .abs();
            let gen_to_storage = FSum::with_all(&output_core.energy_to_storage[key]).value();
            let storage_to_consumption = FSum::with_all(&output_core.energy_from_storage[key])
                .value()
                .abs();
            let gen_to_diverter = FSum::with_all(&output_core.energy_diverted[key]).value();
            let total_gross_import = FSum::with_all(&output_core.energy_import[key]).value();
            let total_gross_export = FSum::with_all(&output_core.energy_export[key]).value();

            let total_energy_into_battery =
                gen_to_storage + FSum::with_all(&output_core.storage_from_grid[key]).value();
            let storage_eff = if total_energy_into_battery > 0. {
                storage_to_consumption / total_energy_into_battery
            } else {
                f64::NAN
            };

            energy_supply_stats.insert(
                key.clone(),
                OutputSummaryEnergySupply {
                    generation: total_generated,
                    consumption: total_consumed,
                    generation_to_consumption: output_core.energy_generated_consumed[key]
                        .iter()
                        .sum(),
                    generation_to_grid,
                    grid_to_consumption,
                    net_import: total_gross_import + total_gross_export,
                    generation_to_storage: gen_to_storage,
                    storage_to_consumption,
                    grid_to_storage: FSum::with_all(&output_core.storage_from_grid[key])
                        .value()
                        .abs(),
                    generation_to_diverter: gen_to_diverter,
                    storage_efficiency: storage_eff,
                    total_gross_import,
                    total_gross_export,
                },
            );
        }

        // Delivered energy by end-use and by fuel
        // TODO (from Python): Ensure end_uses not consuming fuel directly are filtered out on this report
        let mut delivered_energy_dict: IndexMap<Arc<str>, IndexMap<Arc<str>, f64>> =
            [("total".into(), IndexMap::from([("total".into(), 0.)]))].into();
        for (fuel, end_uses) in &output_core.results_end_user {
            // TODO (from Python) are these keys EnergySupplyType ? Why hot water source names too?
            let fuel_found_in_hot_water_sources =
                self.hot_water_sources.keys().collect_vec().contains(&fuel);

            if !fuel_found_in_hot_water_sources && fuel.as_ref() != "_unmet_demand" {
                delivered_energy_dict.insert(fuel.clone(), IndexMap::from([("total".into(), 0.)]));

                for (end_use, delivered_energy) in end_uses {
                    let sum_delivered_energy = FSum::with_all(delivered_energy).value();
                    if sum_delivered_energy > 0.
                        || is_close!(sum_delivered_energy, 0., rel_tol = 1e-09, abs_tol = 1e-10)
                    {
                        delivered_energy_dict[fuel].insert(end_use.clone(), sum_delivered_energy);
                        delivered_energy_dict[fuel]["total"] += sum_delivered_energy;
                        delivered_energy_dict["total"]
                            .entry(end_use.clone())
                            .and_modify(|v| *v += sum_delivered_energy)
                            .or_insert(sum_delivered_energy);
                        delivered_energy_dict["total"]["total"] += sum_delivered_energy
                    };
                }
            }
        }

        let mut hot_water_demand_daily_75th_percentile_dict: IndexMap<Arc<str>, f64> =
            IndexMap::new();
        let simulation_time: SimulationTime = self.simulation_time.as_ref().into();
        for hws_name in self.hot_water_sources.keys() {
            let daily_hw_demand = convert_profile_to_daily(
                &output_core
                    .hot_water_systems
                    .energy_demand_at_hot_water_source[hws_name.as_ref()],
                simulation_time.step,
            );
            hot_water_demand_daily_75th_percentile_dict
                .insert(hws_name.clone(), percentile(&daily_hw_demand, 75));
        }

        let space_heat_demand_total =
            FSum::with_all(output_core.zone_data.space_heat_demand.values().flatten()).value();

        let space_cool_demand_total =
            FSum::with_all(output_core.zone_data.space_cool_demand.values().flatten()).value();

        OutputSummary {
            total_floor_area: self.total_floor_area,
            space_heat_demand_total,
            space_cool_demand_total,
            electricity_peak_consumption: self.calculate_peak_electricity_consumption(output_core),
            energy_supply: energy_supply_stats,
            delivered_energy: delivered_energy_dict,
            hot_water_demand_daily_75th_percentile: hot_water_demand_daily_75th_percentile_dict,
        }
    }

    /// Finds the peak electricity consumption from the simulation output.
    fn calculate_peak_electricity_consumption(
        &self,
        output_core: &OutputCore,
    ) -> OutputSummaryPeakElectricityConsumption {
        // Initialize the SimulationTime object
        let mut timestep_to_date: IndexMap<OrderedFloat<f64>, DateTime<Utc>> = Default::default();
        // Set the base for any non-leap year
        let base_time = Utc.with_ymd_and_hms(2023, 1, 1, 0, 0, 0).unwrap();
        let simulation_time: SimulationTime = self.simulation_time.as_ref().into();
        // The step must reflect hour or half hour in the year (hour 0 to hour 8759)
        // Starts on the start timestep.
        let mut step = simulation_time.start_time();
        for hours in output_core.timestep_array.iter() {
            timestep_to_date.insert(
                OrderedFloat(step),
                base_time + TimeDelta::hours(*hours as i64),
            );
            step += 1.;
        }

        // Get peak electricity consumption, and when it happens.
        // Initialize the SimulationTime object
        let start_timestep = simulation_time.start_time();
        let stepping = simulation_time.step;

        // Get Energy Supply objects with fuel type 'electricity'.
        let electricity_keys: Vec<Arc<str>> = self
            .input
            .energy_supply
            .iter()
            .filter_map(|(key, value)| {
                (value.fuel == FuelType::Electricity).then_some(Arc::<str>::from(key.to_string()))
            })
            .collect_vec();

        // Calculate net import per timestep by adding gross import and export figures.
        // Add because export figures already negative.
        let net_import_per_timestep = (0..output_core.timestep_array.len())
            .map(|i| {
                FSum::with_all(electricity_keys.iter().map(|key| {
                    output_core.energy_import[key][i] + output_core.energy_export[key][i]
                }))
                .value()
            })
            .collect_vec();

        // Find peak electricity consumption
        let peak_elec_consumption = *net_import_per_timestep
            .iter()
            .max_by(|a, b| a.total_cmp(b))
            .expect("Expected to be able to find a max for a non-empty list of net imports");
        let index_peak_elec_consumption = net_import_per_timestep
            .iter()
            .position(|v| *v == peak_elec_consumption)
            .expect("Expected to be able to find the index for the peak that we just found");

        // must reflect hour or half hour in the year (hour 0 to hour 8759)
        // to work with the dictionary below timestep_to_date
        // hence + start_timestep
        let peak_step = index_peak_elec_consumption as f64 + start_timestep;
        let peak_datetime = timestep_to_date[&OrderedFloat(peak_step)];

        OutputSummaryPeakElectricityConsumption {
            peak: peak_elec_consumption,
            index: index_peak_elec_consumption,
            month: peak_datetime.month() as u8,
            day: peak_datetime.day() as u8,
            hour: ((peak_step % (24. / stepping)) * stepping) as u8,
        }
    }
}

struct HeatCoolOutputs {
    hc_output_convective: IndexMap<Arc<str>, f64>,
    hc_output_radiative: IndexMap<Arc<str>, f64>,
    hc_output_min: IndexMap<Arc<str>, f64>,
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
    temp_setpnt_heat: IndexMap<Arc<str>, f64>,
    temp_setpnt_cool: IndexMap<Arc<str>, f64>,
    frac_convective_heat: IndexMap<Arc<str>, f64>,
    frac_convective_cool: IndexMap<Arc<str>, f64>,
}

#[derive(Clone, Debug)]
pub enum HotWaterResultMap {
    Float(IndexMap<Arc<str>, Vec<f64>>),
    Int(IndexMap<Arc<str>, Vec<usize>>),
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq)]
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

impl From<Option<f64>> for NumberOrDivisionByZero {
    fn from(value: Option<f64>) -> Self {
        match value {
            Some(value) => NumberOrDivisionByZero::Number(value),
            None => NumberOrDivisionByZero::DivisionByZero,
        }
    }
}

impl Serialize for NumberOrDivisionByZero {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        match self {
            NumberOrDivisionByZero::Number(number) => serializer.serialize_f64(*number),
            NumberOrDivisionByZero::DivisionByZero => serializer.serialize_str("DIV/0"),
        }
    }
}

pub type ResultsEndUser = IndexMap<String, IndexMap<String, Vec<f64>>>;

enum SpaceHeatCoolSystems<'a> {
    Heat(&'a IndexMap<Arc<str>, Arc<Mutex<SpaceHeatSystem>>>),
    Cool(&'a IndexMap<Arc<str>, AirConditioning>),
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

struct HeatCoolSystemsForZone {
    h_sorted_names: Vec<Arc<str>>,
    c_sorted_names: Vec<Arc<str>>,
    setpoints_and_convective_fractions: SetpointsAndConvectiveFractions,
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
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<InternalGainsCollection> {
    let mut gains_collection = InternalGainsCollection::from([]);
    if let Some(internal_gains) = input.total_internal_gains.as_ref() {
        gains_collection.insert(
            "total_internal_gains".into(),
            Gains::Internal(internal_gains_from_details(
                internal_gains,
                total_floor_area,
                simulation_time_iterator,
            )?),
        );
    }
    if let Some(internal_gains) = input.metabolic_gains.as_ref() {
        gains_collection.insert(
            "metabolic_gains".into(),
            Gains::Internal(internal_gains_from_details(
                internal_gains,
                total_floor_area,
                simulation_time_iterator,
            )?),
        );
    }
    if let Some(internal_gains) = input.evaporative_losses.as_ref() {
        gains_collection.insert(
            "evaporative_losses".into(),
            Gains::Internal(internal_gains_from_details(
                internal_gains,
                total_floor_area,
                simulation_time_iterator,
            )?),
        );
    }
    if let Some(internal_gains) = input.cold_water_losses.as_ref() {
        gains_collection.insert(
            "coldwaterlosses".into(),
            Gains::Internal(internal_gains_from_details(
                internal_gains,
                total_floor_area,
                simulation_time_iterator,
            )?),
        );
    }
    if let Some(internal_gains) = input.other.as_ref() {
        gains_collection.insert(
            "other".into(),
            Gains::Internal(internal_gains_from_details(
                internal_gains,
                total_floor_area,
                simulation_time_iterator,
            )?),
        );
    }

    Ok(gains_collection)
}

fn internal_gains_from_details(
    details: &InternalGainsDetails,
    total_floor_area: f64,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<InternalGains> {
    InternalGains::new(
        convert_energy_to_wm2(details, total_floor_area)?,
        details.start_day,
        details.time_series_step,
        simulation_time_iterator,
    )
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
) -> anyhow::Result<IndexMap<String, Arc<Mutex<WwhrsInstantaneous>>>> {
    let mut wwhr_systems: IndexMap<String, Arc<Mutex<WwhrsInstantaneous>>> = IndexMap::from([]);
    if let Some(systems) = wwhrs {
        for (name, system) in systems {
            wwhr_systems
                .entry(name.into())
                .or_insert(Arc::new(Mutex::new(wwhr_system_from_details(
                    system.clone(),
                    cold_water_sources,
                )?)));
        }
    }

    Ok(wwhr_systems)
}

fn wwhr_system_from_details(
    system: WasteWaterHeatRecoveryDetails,
    cold_water_sources: &ColdWaterSources,
) -> anyhow::Result<WwhrsInstantaneous> {
    let cold_water_source = cold_water_sources.get(&system.cold_water_source).ok_or_else(|| anyhow!("Cold water source '{}' referenced by WWHRS input not found in cold water sources list.", system.cold_water_source))?;

    // Get efficiency data for all systems if provided
    let WasteWaterHeatRecoveryDetails {
        flow_rates,
        system_a_efficiencies,
        system_a_utilisation_factor,
        system_b_efficiencies,
        system_b_utilisation_factor,
        system_c_efficiencies,
        system_c_utilisation_factor,
        system_b_efficiency_factor,
        system_c_efficiency_factor,
        ..
    } = system;

    WwhrsInstantaneous::new(
        flow_rates,
        system_a_efficiencies,
        cold_water_source.clone(),
        system_a_utilisation_factor,
        system_b_efficiencies,
        system_b_utilisation_factor,
        system_c_efficiencies,
        system_c_utilisation_factor,
        system_b_efficiency_factor.into(),
        system_c_efficiency_factor.into(),
    )
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

fn temp_internal_air_for_zones(zones: &IndexMap<Arc<str>, Arc<Zone>>, total_volume: f64) -> f64 {
    let internal_air_temperature = zones
        .values()
        .map(|zone| zone.temp_internal_air() * zone.volume())
        .sum::<f64>();

    internal_air_temperature / total_volume
}

pub(crate) type TempInternalAirFn = Arc<dyn Fn() -> f64 + Send + Sync>;

fn shareable_fn(num: &Arc<AtomicF64>) -> TempInternalAirFn {
    let clone = num.clone();
    Arc::from(move || clone.load(Ordering::SeqCst))
}

pub(crate) type HeatBalanceAllResults =
    IndexMap<HeatBalanceFieldName, IndexMap<Arc<str>, IndexMap<Arc<str>, Vec<f64>>>>;

struct SpaceHeatingCalculation {
    gains_internal_zone: HashMap<Arc<str>, f64>,
    gains_solar_zone: HashMap<Arc<str>, f64>,
    operative_temp: HashMap<Arc<str>, f64>,
    internal_air_temp: HashMap<Arc<str>, f64>,
    space_heat_demand_zone: HashMap<Arc<str>, f64>,
    space_cool_demand_zone: HashMap<Arc<str>, f64>,
    space_heat_provided_system: HashMap<Arc<str>, f64>,
    space_cool_provided_system: HashMap<Arc<str>, f64>,
    internal_gains_ductwork: f64,
    heat_balance_map: HashMap<Arc<str>, Option<HeatBalance>>,
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
    heat_system_name_for_zone: &mut IndexMap<Arc<str>, Vec<Arc<str>>>,
    cool_system_name_for_zone: &mut IndexMap<Arc<str>, Vec<Arc<str>>>,
    external_conditions: Arc<ExternalConditions>,
    infiltration_ventilation: Arc<InfiltrationVentilation>,
    window_adjust_control: Option<Arc<dyn ControlBehaviour>>,
    controls: &Controls,
    print_heat_balance: bool,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<Zone> {
    let heat_system_name = input.space_heat_system.clone();
    let cool_system_name = input.space_cool_system.clone();

    let heat_system_names: Vec<Arc<str>> = match heat_system_name {
        SystemReference::None(_) => vec!["".into()], // equivalent of [None] in Python - we are using empty string to denote absence rather than using Option<String> everywhere
        SystemReference::Single(name) => vec![name.clone()],
        SystemReference::Multiple(names) => names.clone(),
    }
    .into_iter()
    .map(|name| name.to_string().into())
    .collect();

    for zone_h_name in heat_system_name_for_zone.values() {
        let zone_h_name_set: HashSet<Arc<str>> = HashSet::from_iter(zone_h_name.iter().cloned());
        let h_overassigned: Vec<Arc<str>> =
            HashSet::from_iter(heat_system_names.clone().into_iter())
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

    let cool_system_names: Vec<Arc<str>> = match cool_system_name {
        SystemReference::None(_) => vec!["".into()], // equivalent of [None] in Python - we are using empty string to denote absence rather than using Option<String> everywhere
        SystemReference::Single(name) => vec![name.clone()],
        SystemReference::Multiple(names) => names.clone(),
    }
    .into_iter()
    .map(|name| name.to_string().into())
    .collect();

    for zone_c_name in cool_system_name_for_zone.values() {
        let zone_c_name_set: HashSet<Arc<str>> = HashSet::from_iter(zone_c_name.iter().cloned());
        let c_overassigned: Vec<Arc<str>> =
            HashSet::from_iter(cool_system_names.clone().into_iter())
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
    controls: &Controls,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<Arc<BuildingElement>> {
    Ok(Arc::from(match input {
        BuildingElementInput::Opaque {
            is_unheated_pitched_roof,
            area_input,
            pitch,
            solar_absorption_coeff,
            u_value_input,
            areal_heat_capacity,
            mass_distribution_class,
            orientation360: orientation,
            base_height,
            ..
        } => {
            let is_unheated_pitched_roof = if *pitch < PITCH_LIMIT_HORIZ_CEILING {
                is_unheated_pitched_roof
                    .ok_or_else(|| anyhow!("Pitch of opaque building element was {pitch} degrees, so it is necessary for this element to indicate whether this is an unheated pitched roof."))?
            } else {
                false
            };

            let BuildingElementHeightWidthInput { height, width } =
                area_input.height_and_width.ok_or_else(|| {
                    anyhow!("Height and width of opaque building element must be provided.")
                })?;

            BuildingElement::Opaque(BuildingElementOpaque::new(
                area_input.area(),
                is_unheated_pitched_roof,
                *pitch,
                *solar_absorption_coeff,
                init_resistance_or_uvalue_from_input_struct(u_value_input, *pitch)?,
                *areal_heat_capacity,
                *mass_distribution_class,
                *orientation,
                *base_height,
                height,
                width,
                external_conditions,
            ))
        }
        BuildingElementInput::Transparent {
            area_input,
            u_value_input,
            pitch,
            orientation360: orientation,
            g_value,
            frame_area_fraction,
            base_height,
            shading,
            treatment,
            ..
        } => {
            let BuildingElementHeightWidthInput { height, width } =
                area_input.height_and_width.ok_or_else(|| {
                    anyhow!("Height and width of transparent building element must be provided.")
                })?;

            BuildingElement::Transparent(BuildingElementTransparent::new(
                *pitch,
                init_resistance_or_uvalue_from_input_struct(u_value_input, *pitch)?,
                Some(*orientation),
                *g_value,
                *frame_area_fraction,
                *base_height,
                height,
                width,
                Some(shading.clone()),
                treatment
                    .iter()
                    .map(|t| {
                        WindowTreatment::from_input(
                            t,
                            controls,
                            simulation_time_iterator.current_hour(),
                        )
                    })
                    .collect_vec(),
                external_conditions,
            ))
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
                *area,
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
                *area,
                *pitch,
                init_resistance_or_uvalue_from_input_struct(u_value_input, *pitch)?,
                *thermal_resistance_unconditioned_space,
                *areal_heat_capacity,
                *mass_distribution_class,
                external_conditions,
            ),
        ),
        BuildingElementInput::PartyWall {
            area,
            pitch,
            u_value_input,
            party_wall_cavity_data,
            areal_heat_capacity,
            mass_distribution_class,
        } =>
        // Handle party wall with automatic or manual cavity resistance
        {
            BuildingElement::PartyWall(BuildingElementPartyWall::new(
                *area,
                *pitch,
                init_resistance_or_uvalue_from_input_struct(u_value_input, *pitch)?,
                PartyWallCavityType::from(*party_wall_cavity_data),
                party_wall_cavity_data.party_wall_lining_type(),
                party_wall_cavity_data.thermal_resistance_cavity(),
                *areal_heat_capacity,
                *mass_distribution_class,
                external_conditions,
            )?)
        }
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
                simulation_time,
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
    simulation_time_iterator: &SimulationTimeIterator,
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

    ApplianceGains::new(
        total_energy_supply,
        input.gains_fraction,
        input.start_day,
        input.time_series_step,
        simulation_time_iterator,
        energy_supply_connection,
    )
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
    HeatBattery(HeatBattery),
}

#[derive(Clone, Debug)]
pub(crate) enum HeatBattery {
    DryCore(Arc<HeatBatteryDryCore>),
    Pcm(Arc<RwLock<HeatBatteryPcm>>),
}

impl HeatBattery {
    pub(crate) fn get_battery_losses(&self) -> f64 {
        match self {
            HeatBattery::DryCore(heat_battery) => heat_battery.get_battery_losses(),
            HeatBattery::Pcm(heat_battery) => heat_battery.read().get_battery_losses(),
        }
    }

    pub(crate) fn timestep_end(&self, simtime: SimulationTimeIteration) -> anyhow::Result<()> {
        match self {
            HeatBattery::DryCore(heat_battery) => heat_battery.timestep_end(simtime)?,
            HeatBattery::Pcm(heat_battery) => heat_battery.read().timestep_end(simtime.index)?,
        }

        Ok(())
    }

    pub(crate) fn output_detailed_results(
        &self,
        hot_water_energy_output: &IndexMap<Arc<str>, Vec<ResultParamValue>>,
        hot_water_source_name_for_heat_battery_service: &IndexMap<Arc<str>, Arc<str>>,
    ) -> Option<(ResultsPerTimestep, ResultsAnnual)> {
        match self {
            HeatBattery::DryCore(drycore) => drycore.output_detailed_results().into(),
            HeatBattery::Pcm(pcm) => pcm
                .read()
                .output_detailed_results(
                    hot_water_energy_output,
                    hot_water_source_name_for_heat_battery_service,
                )
                .unwrap_or_else(|e| panic!("{e}"))
                .into(),
        }
    }
}

impl WetHeatSource {
    pub(crate) fn timestep_end(&self, simtime: SimulationTimeIteration) -> anyhow::Result<()> {
        match self {
            WetHeatSource::HeatPump(heat_pump) => heat_pump.lock().timestep_end(simtime.index)?,
            WetHeatSource::Boiler(boiler) => boiler.write().timestep_end(simtime)?,
            WetHeatSource::Hiu(heat_network) => heat_network.lock().timestep_end(simtime.index)?,
            WetHeatSource::HeatBattery(heat_battery) => heat_battery.timestep_end(simtime)?,
        }

        Ok(())
    }

    pub(crate) fn create_service_hot_water_combi(
        &mut self,
        boiler_data: HotWaterSourceDetails,
        service_name: &str,
        temp_hot_water: f64,
        cold_feed: WaterSupply,
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
        hot_water_energy_output: &IndexMap<Arc<str>, Vec<ResultParamValue>>,
        hot_water_source_name_for_heat_battery_service: &IndexMap<Arc<str>, Arc<str>>,
    ) -> Option<(ResultsPerTimestep, ResultsAnnual)> {
        match self {
            WetHeatSource::HeatPump(heat_pump) => heat_pump
                .lock()
                .clone()
                .output_detailed_results(
                    hot_water_energy_output,
                    hot_water_source_name_for_heat_battery_service,
                )
                .ok(),
            WetHeatSource::HeatBattery(heat_battery) => heat_battery.output_detailed_results(
                hot_water_energy_output,
                hot_water_source_name_for_heat_battery_service,
            ),
            _ => None,
        }
    }
}

// TODO - this enum is a placeholder, to review later and potentially implement a trait instead
#[derive(Debug)]
pub(crate) enum HeatSystem {
    WetSystem(WetHeatSource),
    WwhrsSystem(Arc<Mutex<WwhrsInstantaneous>>),
}

impl HeatSystem {
    pub(crate) fn timestep_end(&self, simtime: SimulationTimeIteration) -> anyhow::Result<()> {
        match self {
            HeatSystem::WetSystem(wet_heat_source) => wet_heat_source.timestep_end(simtime),
            HeatSystem::WwhrsSystem(wwhrs) => {
                wwhrs.lock().timestep_end();
                Ok(())
            }
        }
    }
}

pub type ResultsPerTimestep =
    IndexMap<Arc<str>, IndexMap<(Arc<str>, Option<Arc<str>>), Vec<ResultParamValue>>>;
pub type ResultsAnnual =
    IndexMap<Arc<str>, IndexMap<(Arc<str>, Option<Arc<str>>), ResultParamValue>>;

#[derive(Clone, Debug, Deserialize, PartialEq)]
pub enum ResultParamValue {
    String(String),
    Number(f64),
    Boolean(bool),
    Empty,
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

impl From<&f64> for ResultParamValue {
    fn from(value: &f64) -> Self {
        Self::Number(*value)
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

impl From<bool> for ResultParamValue {
    fn from(value: bool) -> Self {
        Self::Boolean(value)
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

impl<T: Into<ResultParamValue>> From<Option<T>> for ResultParamValue {
    fn from(value: Option<T>) -> Self {
        match value {
            Some(value) => value.into(),
            None => Self::Empty,
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
                ResultParamValue::Empty => "".to_string(),
            }
        )
    }
}

impl Serialize for ResultParamValue {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match self {
            ResultParamValue::String(string) => serializer.serialize_str(string),
            ResultParamValue::Number(number) => serializer.serialize_f64(*number),
            ResultParamValue::Boolean(boolean) => serializer.serialize_bool(*boolean),
            ResultParamValue::Empty => serializer.serialize_none(),
        }
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
                    simulation_time.as_ref(),
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
            power_circ_pump,
            power_aux,
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
                power_circ_pump.unwrap_or(0.06),
                power_aux.unwrap_or(0.),
                *building_level_distribution_losses,
                energy_supply,
                energy_supply_conn_name_auxiliary,
                energy_supply_conn_name_building_level_distribution_losses,
                simulation_time.step_in_hours(),
            )))))
        }
        HeatSourceWetDetails::HeatBattery { battery } => match battery {
            HeatBatteryInput::Pcm {
                control_charge,
                energy_supply,
                ..
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
                    WetHeatSource::HeatBattery(HeatBattery::Pcm(Arc::new(RwLock::new(HeatBatteryPcm::new(
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
                        None,
                        None,
                        None,
                        None,
                        Some(detailed_output_heating_cooling),
                    )))));
                Ok(heat_source)
            }
            HeatBatteryInput::DryCore {
                control_charge,
                energy_supply,
                number_of_units,
                ..
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

                Ok(WetHeatSource::HeatBattery(HeatBattery::DryCore(HeatBatteryDryCore::new(
                    battery,
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
                    Some(*number_of_units as u32),
                    simulation_time.step_in_hours(),
                    detailed_output_heating_cooling.into(),
                )?)))
            }
        },
    }
}

struct HeatSourceFromInput {
    heat_source: HeatSource,
    energy_supply_conn_name: String,
    heat_source_name_pair: Option<(String, String)>,
}

fn heat_source_from_input(
    name: &str,
    hot_water_source_name: &str,
    input: &HeatSourceInput,
    cold_water_source: &WaterSupply,
    volume: f64,
    daily_losses: f64,
    heat_exchanger_surface_area: Option<f64>,
    wet_heat_sources: &IndexMap<String, WetHeatSource>,
    simulation_time: &SimulationTimeIterator,
    controls: &Controls,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    temp_internal_air_fn: TempInternalAirFn,
    external_conditions: Arc<ExternalConditions>,
) -> anyhow::Result<HeatSourceFromInput> {
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

            Ok(HeatSourceFromInput {
                heat_source: HeatSource::Storage(HeatSourceWithStorageTank::Immersion(Arc::new(
                    Mutex::new(ImmersionHeater::new(
                        *power,
                        energy_supply_conn,
                        simulation_time.step_in_hours(),
                        Some(control_min),
                        Some(control_max),
                    )),
                ))),
                energy_supply_conn_name: name.into(),
                heat_source_name_pair: None,
            })
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
            orientation360: orientation,
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

            Ok(HeatSourceFromInput {
                heat_source: HeatSource::Storage(HeatSourceWithStorageTank::Solar(Arc::new(
                    Mutex::new(SolarThermalSystem::new(
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
                    )),
                ))),
                energy_supply_conn_name: name.into(),
                heat_source_name_pair: None,
            })
        }
        HeatSourceInput::ServiceWaterRegular {
            name,
            control_min,
            control_max,
            temp_flow_limit_upper,
            ..
        } => {
            let energy_supply_conn_name: String =
                format!("{name}_water_heating: {hot_water_source_name}").into();
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
            let mut heat_source_wet_clone = heat_source_wet.clone();

            Ok(HeatSourceFromInput {
                heat_source: match heat_source_wet_clone {
                    WetHeatSource::HeatPump(heat_pump) => HeatSource::Wet(Box::new(
                        HeatSourceWet::HeatPumpWater(HeatPump::create_service_hot_water(
                            heat_pump.clone(),
                            &energy_supply_conn_name,
                            temp_flow_limit_upper.ok_or_else(|| anyhow!("A temp_flow_limit_upper is needed for heat pump with the name '{name}'"))?,
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
                    WetHeatSource::HeatBattery(battery) => HeatSource::Wet(Box::new(
                        HeatSourceWet::HeatBatteryHotWater(match battery {
                            HeatBattery::DryCore(dry_core) => HeatBatteryWaterService::DryCore(
                                HeatBatteryDryCore::create_service_hot_water_regular(
                                    dry_core,
                                    &energy_supply_conn_name,
                                    cold_water_source.clone(),
                                    control_min,
                                    control_max,
                                )?,
                            ),
                            HeatBattery::Pcm(pcm) => HeatBatteryWaterService::Pcm(
                                HeatBatteryPcm::create_service_hot_water_regular(
                                    pcm,
                                    &energy_supply_conn_name,
                                    cold_water_source.clone(),
                                    control_min,
                                    control_max,
                                )?,
                            ),
                        }),
                    )),
                },
                energy_supply_conn_name: energy_supply_conn_name.clone(),
                heat_source_name_pair: (energy_supply_conn_name, hot_water_source_name.into())
                    .into(),
            })
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

            Ok(HeatSourceFromInput {
                heat_source: HeatSource::Wet(Box::new(HeatSourceWet::HeatPumpWaterOnly(
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
                energy_supply_conn_name: energy_supply_conn_name.into(),
                heat_source_name_pair: None,
            })
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) enum HotWaterSource {
    PreHeated(HotWaterStorageTank),
    CombiBoiler(Arc<BoilerServiceWaterCombi>),
    PointOfUse(Arc<PointOfUse>),
    HeatNetwork(Arc<HeatNetworkServiceWaterDirect>),
    HeatBattery(HeatBatteryHotWaterSource),
}

pub(crate) trait HotWaterSourceBehaviour: std::fmt::Debug + Clone {
    fn get_cold_water_source(&self) -> WaterSupply;
    fn demand_hot_water(
        &self,
        usage_events: Vec<WaterEventResult>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64>;
    fn get_temp_hot_water(
        &self,
        volume_required: f64,
        volume_required_already: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>>;
    fn internal_gains(&self) -> Option<f64> {
        None
    }
    fn get_losses_from_primary_pipework_and_storage(&self) -> (f64, f64) {
        (0., 0.)
    }
    fn is_point_of_use(&self) -> bool {
        false
    }
}

impl HotWaterSourceBehaviour for HotWaterSource {
    fn get_cold_water_source(&self) -> WaterSupply {
        match self {
            HotWaterSource::PreHeated(source) => source.get_cold_water_source(),
            HotWaterSource::CombiBoiler(source) => source.get_cold_water_source().clone(),
            HotWaterSource::PointOfUse(source) => source.get_cold_water_source().clone(),
            HotWaterSource::HeatNetwork(source) => source.get_cold_water_source().clone(),
            HotWaterSource::HeatBattery(source) => source.get_cold_water_source().clone(),
        }
    }

    fn demand_hot_water(
        &self,
        usage_events: Vec<WaterEventResult>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        Ok(match self {
            HotWaterSource::PreHeated(hot_water_storage_tank) => {
                hot_water_storage_tank.demand_hot_water(usage_events, simtime)?
            }
            HotWaterSource::CombiBoiler(boiler_service_water_combi) => {
                boiler_service_water_combi.demand_hot_water(usage_events, simtime)?
            }
            HotWaterSource::PointOfUse(point_of_use) => {
                point_of_use.demand_hot_water(usage_events, &simtime)?
            }
            HotWaterSource::HeatNetwork(heat_network_service_water_direct) => {
                heat_network_service_water_direct.demand_hot_water(usage_events, simtime)?
            }
            HotWaterSource::HeatBattery(heat_battery_hot_water_source) => {
                heat_battery_hot_water_source.demand_hot_water(usage_events, simtime)?
            }
        })
    }

    fn get_temp_hot_water(
        &self,
        volume_required: f64,
        volume_required_already: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        match self {
            HotWaterSource::PreHeated(hot_water_storage_tank) => hot_water_storage_tank
                .get_temp_hot_water(volume_required, volume_required_already, simtime),
            HotWaterSource::CombiBoiler(boiler_service_water_combi) => {
                Ok(boiler_service_water_combi
                    .get_temp_hot_water(volume_required, Some(volume_required_already)))
            }
            HotWaterSource::PointOfUse(point_of_use) => {
                Ok(point_of_use.get_temp_hot_water(volume_required, Some(volume_required_already)))
            }
            HotWaterSource::HeatNetwork(heat_network_service_water_direct) => {
                Ok(heat_network_service_water_direct
                    .get_temp_hot_water(volume_required, Some(volume_required_already)))
            }
            HotWaterSource::HeatBattery(heat_battery_hot_water_source) => {
                heat_battery_hot_water_source.get_temp_hot_water(
                    volume_required,
                    volume_required_already,
                    simtime,
                )
            }
        }
    }

    // Calls internal_gains on hot water source where available
    fn internal_gains(&self) -> Option<f64> {
        match &self {
            HotWaterSource::PreHeated(hot_water_storage_tank) => {
                hot_water_storage_tank.internal_gains()
            }
            HotWaterSource::CombiBoiler(boiler_service_water_combi) => {
                Some(boiler_service_water_combi.internal_gains())
            }
            HotWaterSource::PointOfUse(_) => None,
            HotWaterSource::HeatNetwork(_) => None,
            HotWaterSource::HeatBattery(_) => None,
        }
    }

    // Calls get_losses_from_primary_pipework_and_storage on hot water source where available, otherwise returns 0s.
    fn get_losses_from_primary_pipework_and_storage(&self) -> (f64, f64) {
        match &self {
            HotWaterSource::PreHeated(hot_water_storage_tank) => {
                hot_water_storage_tank.get_losses_from_primary_pipework_and_storage()
            }
            _ => (0., 0.),
        }
    }

    fn is_point_of_use(&self) -> bool {
        matches!(&self, HotWaterSource::PointOfUse(_))
    }
}

#[derive(Clone, Debug)]
pub(crate) enum HeatBatteryHotWaterSource {
    Pcm(Arc<HeatBatteryPcmServiceWaterDirect<WaterSupply>>),
    DryCore(Arc<HeatBatteryDryCoreServiceWaterDirect<WaterSupply>>),
}

impl HotWaterSourceBehaviour for HeatBatteryHotWaterSource {
    fn get_cold_water_source(&self) -> WaterSupply {
        match self {
            HeatBatteryHotWaterSource::Pcm(pcm) => pcm.get_cold_water_source().clone(),
            HeatBatteryHotWaterSource::DryCore(dry_core) => {
                dry_core.get_cold_water_source().clone()
            }
        }
    }

    fn demand_hot_water(
        &self,
        usage_events: Vec<WaterEventResult>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        match self {
            HeatBatteryHotWaterSource::Pcm(pcm) => {
                pcm.demand_hot_water(usage_events.into(), simtime)
            }
            HeatBatteryHotWaterSource::DryCore(dry_core) => {
                dry_core.demand_hot_water(usage_events.into(), simtime)
            }
        }
    }

    fn get_temp_hot_water(
        &self,
        volume_required: f64,
        volume_required_already: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        match self {
            HeatBatteryHotWaterSource::Pcm(pcm) => {
                pcm.get_temp_hot_water(volume_required, volume_required_already.into(), simtime)
            }
            HeatBatteryHotWaterSource::DryCore(dry_core) => dry_core.get_temp_hot_water(
                volume_required,
                volume_required_already.into(),
                simtime,
            ),
        }
    }
}

fn hot_water_source_from_input(
    name: &str,
    input: &HotWaterSourceDetails,
    cold_water_sources: &ColdWaterSources,
    pre_heated_water_sources: &IndexMap<String, HotWaterStorageTank>,
    wet_heat_sources: &mut IndexMap<String, WetHeatSource>,
    wwhrs: &IndexMap<String, Arc<Mutex<WwhrsInstantaneous>>>,
    controls: &Controls,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    diverter_types: &DiverterTypes,
    diverters: &mut Vec<Arc<RwLock<PVDiverter>>>,
    temp_internal_air_fn: TempInternalAirFn,
    simulation_time: &SimulationTimeIterator,
    external_conditions: Arc<ExternalConditions>,
    detailed_output_heating_cooling: bool,
    used_heat_source_names: &mut HashSet<String>,
    cold_water_sources_already_allocated: &mut HashSet<String>,
) -> anyhow::Result<(HotWaterSource, Vec<String>, IndexMap<String, String>)> {
    let mut energy_supply_conn_names = vec![];
    let mut hot_water_source_name_for_service: IndexMap<String, String> = Default::default();
    let cloned_input = input.clone();

    let cold_water_source_for_hot_water_tank =
        |cold_water_source_type: &str,
         cold_water_sources_already_allocated: &mut HashSet<String>|
         -> anyhow::Result<WaterSupply> {
            cold_water_sources.get(cold_water_source_type).map(|source| Ok(WaterSupply::ColdWaterSource(source.clone()))).or_else(|| {
                let source = pre_heated_water_sources
                    .get(cold_water_source_type)
                    .map(|source| WaterSupply::Preheated(source.clone()))
                    .or(wwhrs.get(cold_water_source_type).map(|source| WaterSupply::Wwhrs(source.clone())));
                if source.is_some() {
                    if cold_water_sources_already_allocated.contains(cold_water_source_type) {
                        return Some(Err(anyhow!("Cannot allocate PreHeatedWaterSource or WWHRS to more than one HotWaterSource: {cold_water_source_type}")));
                    } else {
                        cold_water_sources_already_allocated.insert(cold_water_source_type.into());
                    }
                }

                Ok(source).transpose()
            })
                .transpose()?
                .ok_or_else(|| anyhow!("Could not find pre-heated or WWHRS water source for name '{cold_water_source_type}'"))
        };

    let mut heat_sources_for_hot_water_tank = |cold_water_source: WaterSupply,
                                               heat_exchanger_surface_area: &Option<f64>,
                                               heat_source: &IndexMap<
        std::string::String,
        HeatSourceInput,
    >,
                                               volume: &f64,
                                               daily_losses: &f64|
     -> anyhow::Result<
        IndexMap<String, PositionedHeatSource>,
    > {
        let mut heat_sources: IndexMap<String, PositionedHeatSource> = Default::default();

        let heat_exchanger_surface_area = heat_exchanger_surface_area.and_then(|surface_area| {
            heat_source
                .values()
                .any(|source| matches!(source, HeatSourceInput::HeatPumpHotWaterOnly { .. }))
                .then_some(surface_area)
        });

        // With pre-heated tanks we allow now tanks not to have a heat source as the 'cold' feed
        // could be a pre-heated source or wwhr that might be enough
        for (heat_source_name, heat_source_data) in heat_source {
            let heat_source_name = String::from(heat_source_name);
            if used_heat_source_names.contains(&heat_source_name) {
                return Err(anyhow!(
                    "Duplicate heat source name detected: {heat_source_name}"
                ));
            }
            used_heat_source_names.insert(heat_source_name.clone());

            let heater_position = heat_source_data.heater_position();
            let thermostat_position = match input {
                HotWaterSourceDetails::StorageTank { .. } => heat_source_data.thermostat_position(),
                HotWaterSourceDetails::SmartHotWaterTank { .. } => None,
                _ => {
                    unreachable!()
                }
            };

            let HeatSourceFromInput {
                heat_source,
                energy_supply_conn_name,
                heat_source_name_pair,
            } = heat_source_from_input(
                heat_source_name.as_str(),
                name,
                heat_source_data,
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
                heat_source_name,
                PositionedHeatSource {
                    heat_source: heat_source.clone(),
                    heater_position,
                    thermostat_position,
                },
            );
            energy_supply_conn_names.push(energy_supply_conn_name);
            if let Some((energy_supply_conn_name, hot_water_source_name)) = heat_source_name_pair {
                hot_water_source_name_for_service
                    .insert(energy_supply_conn_name, hot_water_source_name);
            }
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
            let cold_water_source = cold_water_source_for_hot_water_tank(
                cold_water_source_type,
                cold_water_sources_already_allocated,
            )?;
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
                &simulation_time.current_iteration(),
                heat_sources.clone(),
                temp_internal_air_fn,
                external_conditions,
                detailed_output_heating_cooling,
                Some(24),
                primary_pipework_lst,
                *WATER,
                None,
                None,
                None,
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
            let cold_water_source = cold_water_source_for_hot_water_tank(
                cold_water_source_type,
                cold_water_sources_already_allocated,
            )?;

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
            let pump_source_name = format!("Smart_hot_water_tank_pump: {name}");
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
                &simulation_time.current_iteration(),
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
                    .expect("expected to be able to instantiate a combi boiler object")
                    .into(),
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
            let energy_supply_conn_name = name;
            energy_supply_conn_names.push(energy_supply_conn_name.to_string().into());
            let energy_supply_conn =
                EnergySupply::connection(energy_supply.clone(), energy_supply_conn_name)?;
            let cold_water_source =
                cold_water_source_for_type(cold_water_source_type, cold_water_sources)?;
            HotWaterSource::PointOfUse(PointOfUse::new(
                efficiency.ok_or_else(|| anyhow!("An efficiency value was expected on a point of use hot water source input."))?,
                energy_supply_conn,
                cold_water_source,
                *setpoint_temp,
            ).into())
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
            HotWaterSource::HeatNetwork(
                HeatNetwork::create_service_hot_water_direct(
                    heat_source_wet.clone(),
                    &energy_supply_conn_name,
                    (*setpoint_temp).ok_or_else(|| {
                        anyhow!(
                        "A setpoint_temp value was expected on a point of use hot water source."
                    )
                    })?,
                    cold_water_source,
                )
                .into(),
            )
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

            HotWaterSource::HeatBattery(match heat_battery {
                HeatBattery::DryCore(dry_core) => HeatBatteryHotWaterSource::DryCore(Arc::new(
                    HeatBatteryDryCore::create_service_hot_water_direct(
                        dry_core,
                        &energy_supply_conn_name,
                        *setpoint_temp,
                        cold_water_source,
                    )?,
                )),
                HeatBattery::Pcm(pcm) => HeatBatteryHotWaterSource::Pcm(Arc::new(
                    HeatBatteryPcm::create_service_hot_water_direct(
                        pcm,
                        &energy_supply_conn_name,
                        *setpoint_temp,
                        cold_water_source,
                    )?,
                )),
            })
        }
    };

    Ok((
        hot_water_source,
        energy_supply_conn_names,
        hot_water_source_name_for_service,
    ))
}

fn cold_water_source_for_type(
    cold_water_source_type: &str,
    cold_water_sources: &ColdWaterSources,
) -> anyhow::Result<WaterSupply> {
    Ok(WaterSupply::ColdWaterSource(
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
    heat_system_names_requiring_overvent: &mut Vec<Arc<str>>,
    heat_system_name_for_zone: &IndexMap<Arc<str>, Vec<Arc<str>>>,
    zones: &IndexMap<Arc<str>, Arc<Zone>>,
    heat_sources_wet_with_buffer_tank: &[String],
    external_conditions: Arc<ExternalConditions>,
    detailed_output_heating_cooling: bool,
    initial_temp: f64,
) -> anyhow::Result<SpaceHeatSystemsWithEnergyConnections> {
    let mut energy_conn_names_for_systems: IndexMap<Arc<str>, Arc<str>> = Default::default();
    let space_heat_systems = input
        .iter()
        .filter(|(system_name, _)| heat_system_name_for_zone.values().flatten().any(|heat_system_name| heat_system_name.as_ref() == system_name.as_str()))
        .map(|(system_name, space_heat_system_details)| {
            let system_name: Arc<str> = system_name.to_string().into();
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
                        let energy_supply_conn = EnergySupply::connection(energy_supply, energy_supply_conn_name.as_ref()).unwrap();
                        SpaceHeatSystem::Instant(InstantElecHeater::new(
                            *rated_power,
                            *frac_convective,
                            energy_supply_conn,
                            simulation_time.step_in_hours(),
                            controls.get_with_string(control),
                        ))
                    }
                    SpaceHeatSystemDetails::ElectricStorageHeater { pwr_in, rated_power_instant, storage_capacity, air_flow_type, frac_convective, fan_pwr, n_units, energy_supply, zone, control, control_charger, dry_core_min_output, dry_core_max_output, state_of_charge_init, .. } => {
                        let zone: Arc<str> = zone.as_str().into();
                        let energy_supply = energy_supplies.get(energy_supply).ok_or_else(|| anyhow!("Space heat system references an undeclared energy supply '{energy_supply}'."))?.clone();
                        let energy_supply_conn_name = system_name.clone();
                        energy_conn_names_for_systems.insert(system_name.clone(), energy_supply_conn_name.clone());
                        let energy_supply_conn = EnergySupply::connection(energy_supply, energy_supply_conn_name.as_ref()).unwrap();

                        let zone = zones.get(&zone).ok_or_else(|| anyhow!("Space heat system references an undeclared zone '{zone}'."))?.clone();
                        let zone_setpoint_init = zone.setpnt_init();
                        let control = controls.get_with_string(control).ok_or_else(|| anyhow!("A control object was expected for an electric storage heater"))?;
                        let charge_control = controls.get_with_string(control_charger).ok_or_else(|| anyhow!("Space heat system references an invalid charge control name '{control_charger}'"))?;
                        SpaceHeatSystem::ElecStorage(ElecStorageHeater::new(*pwr_in, *rated_power_instant, *storage_capacity, *air_flow_type, *frac_convective, *fan_pwr, *n_units, zone_setpoint_init, ZoneTempInternalAir(zone).as_fn(), energy_supply_conn, simulation_time, control, charge_control, dry_core_min_output.clone(), dry_core_max_output.clone(), external_conditions.clone(), *state_of_charge_init, Some(detailed_output_heating_cooling))?)
                    }
                    SpaceHeatSystemDetails::WetDistribution { emitters, energy_supply, flow_data, bypass_fraction_recirculated, heat_source, temp_diff_emit_dsgn, control, thermal_mass, ecodesign_controller, design_flow_temp, zone, pipework, .. } => {
                        let zone: Arc<str> = zone.as_str().into();
                        let heat_source_name = &heat_source.name;
                        let temp_flow_limit_upper = &heat_source.temp_flow_limit_upper;

                        let energy_supply_conn_name = String::from([heat_source_name, "_space_heating: ", &system_name].concat());
                        energy_conn_names_for_systems.insert(system_name.clone(), energy_supply_conn_name.to_string().into());

                        let heat_source = heat_sources_wet.get(&heat_source.name).ok_or_else(|| anyhow!("A heat source name provided under the name '{heat_source_name}' was expected when setting up space heat systems in the calculation corpus."))?;
                        let mut with_buffer_tank = false;

                        let control = controls.get_with_string(control).ok_or_else(|| anyhow!("A control object was expected for wet heat source: '{heat_source_name}'"))?;

                        let heat_source_service: SpaceHeatingService =
                            match heat_source {
                                WetHeatSource::HeatPump(heat_pump) => {
                                    // For HPs, need to know emitter type for inertia calculations
                                    let emitter_type = if emitters.iter().any(|e| matches!(e, WetEmitter::Fancoil { .. })) {
                                        HeatPumpEmitterType::FanCoils
                                    } else {
                                        HeatPumpEmitterType::RadiatorsUfh
                                    };
                                    // For HPs, checking if there's a buffer tank to inform both the service space heating
                                    // and the emitters of its presence.
                                    if heat_sources_wet_with_buffer_tank.contains(heat_source_name) {
                                        with_buffer_tank = true;
                                    }

                                    let volume_heated = total_volume_heated_by_system(zones, heat_system_name_for_zone, &system_name);

                                    let heat_source_service = HeatPump::create_service_space_heating(
                                        heat_pump.clone(),
                                        &energy_supply_conn_name,
                                        emitter_type,
                                        temp_flow_limit_upper.expect("Expected a temp_flow_limit_upper to be present for a heat pump"),
                                        *temp_diff_emit_dsgn,
                                        *design_flow_temp,
                                        control,
                                        volume_heated);

                                    if heat_pump.lock().source_is_exhaust_air() {
                                        // Record heating system as potentially requiring overventilation
                                        heat_system_names_requiring_overvent.push((system_name).clone());
                                    }
                                    SpaceHeatingService::HeatPump(heat_source_service?)
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
                                    SpaceHeatingService::HeatBattery(match heat_battery {
                                        HeatBattery::DryCore(dry_core) => HeatBatteryServiceSpace::DryCore(HeatBatteryDryCore::create_service_space_heating(dry_core.clone(), &energy_supply_conn_name, control.into())?),
                                        HeatBattery::Pcm(pcm) =>
                                            HeatBatteryServiceSpace::Pcm(HeatBatteryPcm::create_service_space_heating(pcm.clone(), &energy_supply_conn_name, control)?)
                                    })
                                }
                            };

                        let energy_supply_fc_conn = if energy_supply.is_none() {
                            None
                        } else {
                            let energy_supply_name = energy_supply.clone().unwrap();
                            let energy_supply = energy_supplies.get(&energy_supply_name).ok_or_else(|| anyhow!("Space heat system references an undeclared energy supply '{energy_supply_name}'."))?.clone();
                            let energy_supply_fc_conn_name: Arc<str> = format!("FC_fan {system_name}").into();
                            energy_conn_names_for_systems.insert(system_name.clone(), energy_supply_fc_conn_name.clone());
                            Some(Arc::new(EnergySupply::connection(energy_supply, energy_supply_fc_conn_name.as_ref()).unwrap()))
                        };

                        let space_heater = Emitters::new(
                            *thermal_mass,
                            emitters,
                            pipework,
                            *temp_diff_emit_dsgn,
                            matches!(flow_data, FlowData::Variable {..}),
                            if let FlowData::Design { design_flow_rate, .. } = flow_data {
                                Some(*design_flow_rate)
                            } else {
                                None
                            },
                            if let FlowData::Variable { min_flow_rate, .. } = flow_data {
                                Some(*min_flow_rate)
                            } else {
                                None
                            },
                            if let FlowData::Variable { max_flow_rate, .. } = flow_data {
                                Some(*max_flow_rate)
                            } else {
                                None
                            },
                            *bypass_fraction_recirculated,
                            Arc::new(RwLock::new(heat_source_service)),
                            zones.get(&zone).ok_or_else(|| anyhow!("Space heat system wet distribution had reference to undeclared zone with name '{zone}'"))?.clone(),
                            // zone area
                            external_conditions.clone(),
                            *ecodesign_controller,
                            *design_flow_temp,
                            initial_temp,
                            simulation_time.total_steps(),
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
                        let energy_supply_conn_name: Arc<str> = [heat_source_name, "_space_heating: ", &system_name].concat().into();
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
    IndexMap<Arc<str>, Arc<Mutex<SpaceHeatSystem>>>,
    IndexMap<Arc<str>, Arc<str>>,
);

fn space_cool_systems_from_input(
    input: &SpaceCoolSystemInput,
    cool_system_names_for_zone: Vec<&str>,
    controls: &Controls,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<IndexMap<Arc<str>, AirConditioning>> {
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
            let energy_supply_conn_name: Arc<str> = system_name.to_string().into();
            let energy_supply_conn =
                EnergySupply::connection(energy_supply, energy_supply_conn_name.as_ref()).unwrap();
            let control = controls.get_with_string(control).ok_or_else(|| anyhow!("The control reference '{control}' was expected to refer to a known control."))?;

            Ok((
                energy_supply_conn_name,
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
    input: &OnSiteGenerationInput,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    external_conditions: Arc<ExternalConditions>,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<IndexMap<String, PhotovoltaicSystem>> {
    input
        .iter()
        .map(|(name, generation_details)| {
            Ok((name.into(), {
                let (
                    panels,
                    inverter_peak_power_dc,
                    inverter_peak_power_ac,
                    inverter_is_inside,
                    inverter_type,
                    energy_supply
                ) = match generation_details {
                    PhotovoltaicInputs::DeprecatedStyle(PhotovoltaicSystemInput {
                        peak_power,
                        ventilation_strategy,
                        pitch,
                                                            orientation360: orientation,
                        base_height,
                        height,
                        width,
                        energy_supply,
                        shading,
                        inverter_peak_power_dc,
                        inverter_peak_power_ac,
                        inverter_is_inside,
                        inverter_type, ..
                    }) => {
                        (
                            vec![PhotovoltaicPanel::new(*peak_power, *ventilation_strategy, *pitch, *orientation, *base_height, *height, *width, simulation_time_iterator.step_in_hours(), shading.to_vec())],
                            *inverter_peak_power_dc,
                            *inverter_peak_power_ac,
                            *inverter_is_inside,
                            *inverter_type,
                            energy_supply
                        )
                    }
                    PhotovoltaicInputs::WithPanels(PhotovoltaicSystemWithPanelsInput {
                       energy_supply, inverter_is_inside, inverter_peak_power_ac, inverter_peak_power_dc, inverter_type, panels, ..
                    }) => {
                        (
                            panels
                                .iter()
                                .map(|panel| PhotovoltaicPanel::new(
                                    panel.peak_power,
                                    panel.ventilation_strategy,
                                    panel.pitch,
                                    panel.orientation360,
                                    panel.base_height,
                                    panel.height,
                                    panel.width,
                                    simulation_time_iterator.step_in_hours(),
                                    panel.shading.to_vec())
                                )
                                .collect(),
                            *inverter_peak_power_dc,
                            *inverter_peak_power_ac,
                            *inverter_is_inside,
                            *inverter_type,
                            energy_supply
                        )
                    }
                };

                let energy_supply = energy_supplies.get(energy_supply).ok_or_else(|| anyhow!("On site generation (photovoltaic) references an undeclared energy supply '{energy_supply}'."))?.clone();
                let energy_supply_conn = EnergySupply::connection(energy_supply, name).unwrap();
                let inverter = Inverter::new(
                    energy_supply_conn,
                    simulation_time_iterator.step_in_hours(),
                    inverter_peak_power_dc,
                    inverter_peak_power_ac,
                    inverter_is_inside,
                    inverter_type,
                );
                PhotovoltaicSystem::new(
                    external_conditions.clone(),
                    panels,
                    inverter,
                )
            }))
        })
        .collect::<anyhow::Result<IndexMap<_, _>>>()
}

fn total_volume_heated_by_system(
    zones: &IndexMap<Arc<str>, Arc<Zone>>,
    heat_system_name_for_zone: &IndexMap<Arc<str>, Vec<Arc<str>>>,
    heat_system_name: &str,
) -> f64 {
    FSum::with_all(zones.iter().filter_map(|(z_name, zone)| {
        if let Some(system_names) = heat_system_name_for_zone.get(z_name) {
            (system_names
                .iter()
                .any(|name| heat_system_name == name.as_ref()))
            .then(|| zone.volume())
        } else {
            None
        }
    }))
    .value()
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

#[cfg(test)]
mod tests {
    use crate::corpus::Corpus;
    use crate::input::{HotWaterSourceDetails, Input};
    use rstest::{fixture, rstest};
    use serde_json::json;
    use std::sync::Arc;

    #[fixture]
    fn minimal_input() -> Input {
        const NUM_TIMESTEPS: usize = 8;
        let air_temperatures = [10.; NUM_TIMESTEPS];
        let wind_speeds = [4.; NUM_TIMESTEPS];
        let wind_directions = [180.; NUM_TIMESTEPS];
        let diffuse_horizontal_radiation = [100.; NUM_TIMESTEPS];
        let direct_beam_radiation = [200.; NUM_TIMESTEPS];
        let solar_reflectivity_of_ground = [0.2; NUM_TIMESTEPS];
        let schedule_main = [0.; NUM_TIMESTEPS];
        let temperatures = [10.; NUM_TIMESTEPS];

        serde_json::from_value(json!({
            "temp_internal_air_static_calcs": 20.0,
            "SimulationTime": {
                "start": 0,
                "end": NUM_TIMESTEPS,
                "step": 1,
            },
            "ExternalConditions": {
                "air_temperatures": air_temperatures,
                "wind_speeds": wind_speeds,
                "wind_directions": wind_directions,
                "diffuse_horizontal_radiation": diffuse_horizontal_radiation,
                "direct_beam_radiation": direct_beam_radiation,
                "solar_reflectivity_of_ground": solar_reflectivity_of_ground,
                "latitude": 51.5,
                "longitude": -0.1,
                "direct_beam_conversion_needed": false,
                "shading_segments": [
                    {"start360": 0, "end360": 360},
                ],
            },
            "InternalGains": {
                "total_internal_gains": {
                    "start_day": 0,
                    "time_series_step": 1,
                    "schedule": {
                        "main": schedule_main,
                    },
                }
            },
            "ApplianceGains": {},
            "Zone": {
                "zone1": {
                    "area": 50.0,
                    "volume": 100.0,
                    "temp_setpnt_init": 21.0,
                    "BuildingElement": {
                        "wall1": {
                            "type": "BuildingElementOpaque",
                            "width": 5.0,
                            "height": 2.5,
                            "area": 12.5,
                            "base_height": 0,
                            "orientation360": 0,
                            "pitch": 90,
                            "solar_absorption_coeff": 0.6,
                            "u_value": 0.18,
                            "areal_heat_capacity": 110000,
                            "mass_distribution_class": "IE",
                        }
                    },
                    "ThermalBridging": {},
                }
            },
            "ColdWaterSource": {
                "mains water": {
                    "start_day": 0,
                    "time_series_step": 1,
                    "temperatures": temperatures,
                }
            },
            "EnergySupply": {
                "mains elec": {
                    "fuel": "electricity",
                    "is_export_capable": false,
                }
            },
            "Control": {
                "min_temp": {
                    "type": "SetpointTimeControl",
                    "start_day": 0,
                    "time_series_step": 1,
                    "schedule": {
                        "main": [
                            {
                                "value": 50.0,
                                "repeat": NUM_TIMESTEPS,
                            }
                        ]
                    },
                },
                "setpoint_temp_max": {
                    "type": "SetpointTimeControl",
                    "start_day": 0,
                    "time_series_step": 1,
                    "schedule": {
                        "main": [
                            {
                                "value": 60.0,
                                "repeat": NUM_TIMESTEPS,
                            }
                        ]
                    },
                },
            },
            "Events": {
                "Shower": {},
                "Bath": {},
                "Other": {},
            },
            "InfiltrationVentilation": {
                "cross_vent_possible": false,
                "shield_class": "Normal",
                "terrain_class": "OpenField",
                "ventilation_zone_base_height": 2.5,
                "altitude": 30,
                "Vents": {},
                "Leaks": {
                    "ventilation_zone_height": 6,
                    "test_pressure": 50,
                    "test_result": 1.2,
                    "env_area": 220,
                },
            },
            "HotWaterDemand": {
                "Shower": {},
                "Bath": {},
                "Other": {},
                "Distribution": {},
            },
            "HotWaterSource": {},
            "WWHRS": {},
        }))
        .unwrap()
    }

    fn create_preheated_water_source_storage_tank(
        name: &str,
        cold_water_source_name: &str,
    ) -> HotWaterSourceDetails {
        serde_json::from_value(json!(
        {"type": "StorageTank",
        "volume": 24.0,
        "daily_losses": 1.55,
        "init_temp": 48.0,
        "ColdWaterSource": cold_water_source_name,
        "HeatSource": {
            format!("{name}_immersion"): {
                "type": "ImmersionHeater",
                "power": 3.0,
                "EnergySupply": "mains elec",
                "Controlmin": "min_temp",
                "Controlmax": "setpoint_temp_max",
                "heater_position": 0.3,
                "thermostat_position": 0.33}}
            }))
        .unwrap()
    }

    /// Test that PreHeatedWaterSource objects can be initialized in any order
    #[rstest]
    fn test_preheated_water_source_initialization_order_independent(mut minimal_input: Input) {
        let tank1 = create_preheated_water_source_storage_tank("tank1", "tank2");
        let tank2 = create_preheated_water_source_storage_tank("tank2", "mains water");

        minimal_input.pre_heated_water_source =
            serde_json::from_value(json!({"tank1": tank1, "tank2": tank2})).unwrap();

        let corpus =
            Corpus::from_inputs(Arc::new(minimal_input), None, None, &Default::default()).unwrap();

        let pre_heated_sources = corpus.pre_heated_water_sources;

        assert!(pre_heated_sources.contains_key("tank1"));
        assert!(pre_heated_sources.contains_key("tank2"));
    }

    ///  Test that circular references between PreHeatedWaterSource objects are detected and raise an appropriate error.
    #[rstest]
    fn test_preheated_water_source_circular_reference_detection(mut minimal_input: Input) {
        let tank1 = create_preheated_water_source_storage_tank("tank1", "tank2");
        let tank2 = create_preheated_water_source_storage_tank("tank2", "tank1");

        minimal_input.pre_heated_water_source =
            serde_json::from_value(json!({"tank1": tank1, "tank2": tank2})).unwrap();

        let result = Corpus::from_inputs(Arc::new(minimal_input), None, None, &Default::default());

        assert_eq!(
            result.unwrap_err().to_string(),
            "A circular dependency was found between defined preheated water sources."
        )
    }

    /// Test that chains of dependencies work regardless of initialization order
    #[rstest]
    fn test_preheated_water_source_chain_of_dependencies(mut minimal_input: Input) {
        let tank1 = create_preheated_water_source_storage_tank("tank1", "tank2");
        let tank2 = create_preheated_water_source_storage_tank("tank2", "tank3");
        let tank3 = create_preheated_water_source_storage_tank("tank3", "mains water");

        minimal_input.pre_heated_water_source =
            serde_json::from_value(json!({"tank1": tank1, "tank2": tank2, "tank3": tank3}))
                .unwrap();

        let corpus =
            Corpus::from_inputs(Arc::new(minimal_input), None, None, &Default::default()).unwrap();

        let pre_heated_sources = corpus.pre_heated_water_sources;

        assert!(pre_heated_sources.contains_key("tank1"));
        assert!(pre_heated_sources.contains_key("tank2"));
        assert!(pre_heated_sources.contains_key("tank3"));
    }
}
