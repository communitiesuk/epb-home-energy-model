#![allow(dead_code)]

use crate::compare_floats::min_of_2;
use crate::core::energy_supply::energy_supply::EnergySupplies;
use crate::core::heating_systems::wwhrs::{WWHRSInstantaneousSystemB, Wwhrs};
use crate::core::schedule::{expand_events, TypedScheduleEvent};
use crate::core::units::convert_profile_to_daily;
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::dhw_demand::{
    DomesticHotWaterDemand, DomesticHotWaterDemandData,
};
use crate::core::water_heat_demand::misc::water_demand_to_kwh;
use crate::corpus::ColdWaterSources;
use crate::input::{
    BuildType, ColdWaterSourceType, EnergySupplyType, WaterHeatingEventType, WaterPipeContentsType,
    WaterPipework, WaterPipeworkLocation,
};
use crate::simulation_time::SimulationTime;
use crate::statistics::{np_interp, percentile};
use crate::wrappers::future_homes_standard::future_homes_standard::{
    calc_n_occupants, calc_nbeds, create_cold_water_feed_temps, create_hot_water_use_pattern,
    HW_TEMPERATURE, SIMTIME_END, SIMTIME_START, SIMTIME_STEP,
};
use crate::{
    compare_floats::max_of_2,
    core::{
        space_heat_demand::building_element::{pitch_class, HeatFlowDirection},
        units::{LITRES_PER_CUBIC_METRE, SECONDS_PER_HOUR},
    },
    input::{
        BuildingElement, GroundBuildingElement, HeatPumpSourceType,
        HeatSourceWetDetails::{HeatPump, Hiu},
        InputForProcessing, ThermalBridgingDetails, UValueEditableBuildingElement,
    },
    wrappers::future_homes_standard::future_homes_standard::calc_tfa,
};
use anyhow::{anyhow, bail};
use indexmap::IndexMap;
use itertools::Itertools;
use parking_lot::Mutex;
use serde_json::json;
use std::collections::HashMap;
use std::sync::{Arc, LazyLock};

const NOTIONAL_WWHRS: &str = "Notional_Inst_WWHRS";
const NOTIONAL_HP: &str = "notional_HP";
const NOTIONAL_BATH_NAME: &str = "medium";
const NOTIONAL_SHOWER_NAME: &str = "mixer";
const NOTIONAL_OTHER_HW_NAME: &str = "other";

/// Apply assumptions and pre-processing steps for the Future Homes Standard Notional building
pub(crate) fn apply_fhs_notional_preprocessing(
    input: &mut InputForProcessing,
    fhs_notional_a_assumptions: bool,
    _fhs_notional_b_assumptions: bool,
    fhs_fee_notional_a_assumptions: bool,
    fhs_fee_notional_b_assumptions: bool,
) -> anyhow::Result<()> {
    let is_notional_a = fhs_notional_a_assumptions || fhs_fee_notional_a_assumptions;
    let is_fee = fhs_fee_notional_a_assumptions || fhs_fee_notional_b_assumptions;
    // Check if a heat network is present
    let _is_heat_network = check_heatnetwork_present(input);

    // Determine cold water source
    let cold_water_source = input.cold_water_source();

    if cold_water_source.len() != 1 {
        bail!("The FHS Notional wrapper expects exactly one cold water type to be set.");
    }

    let cold_water_source = {
        let first_cold_water_source = cold_water_source.first();
        let (cold_water_source, _) = first_cold_water_source.as_ref().unwrap();
        **cold_water_source
    };

    // Retrieve the number of bedrooms and total volume
    let bedroom_number = input.number_of_bedrooms().ok_or_else(|| {
        anyhow!("FHS Notional wrapper expects number of bedrooms to be provided.")
    })?;
    let total_volume = input.total_zone_volume();

    // Determine the TFA
    let total_floor_area = calc_tfa(input);

    edit_lighting_efficacy(input);
    edit_opaque_adjztu_elements(input)?;
    edit_transparent_element(input)?;
    edit_glazing_for_glazing_limit(input, total_floor_area)?;
    edit_ground_floors(input)?;
    edit_thermal_bridging(input)?;

    // modify bath, shower and other dhw characteristics
    edit_bath_shower_other(input, cold_water_source)?;

    // add WWHRS if needed (and remove any existing systems)
    // TODO enable once function is implemented
    // remove_wwhrs_if_present(input);
    add_wwhrs(input, cold_water_source, is_notional_a, is_fee)?;

    // modify hot water distribution
    edit_hot_water_distribution(input, total_floor_area)?;

    // remove on-site generation, pv diverter or electric battery if present
    // TODO enable following calls once functions implemented
    remove_onsite_generation_if_present(input);
    remove_pv_diverter_if_present(input);
    remove_electric_battery_if_present(input);

    // modify ventilation
    let minimum_ach =
        minimum_air_change_rate(input, total_floor_area, total_volume, bedroom_number);
    // convert to m3/h
    let _minimum_air_flow_rate = minimum_ach * total_volume;
    // TODO enable following call once function implemented
    // edit_infiltration_ventilation(input, is_notional_a, _minimum_air_flow_rate);

    // edit space heating system
    // TODO enable following call once function implemented
    // edit_space_heating_system(input, cold_water_source, total_floor_area, _is_heat_network, is_fee);

    // modify air-conditioning
    edit_space_cool_system(input)?;

    // add solar pv
    add_solar_pv(input, is_notional_a, is_fee, total_floor_area)?;

    Ok(())
}

fn check_heatnetwork_present(input: &InputForProcessing) -> bool {
    let mut is_heat_network = false;

    let heat_source_wet = input.heat_source_wet();

    if heat_source_wet.is_some() {
        let sources = heat_source_wet
            .unwrap()
            .iter()
            .map(|(_, source)| source)
            .collect_vec();

        for source in sources {
            match source {
                Hiu { .. } => is_heat_network = true,
                HeatPump {
                    source_type: HeatPumpSourceType::HeatNetwork,
                    ..
                } => is_heat_network = true,
                _ => {}
            };
        }
    }
    is_heat_network
}

/// Apply notional lighting efficacy
/// efficacy = 120 lm/W
fn edit_lighting_efficacy(input: &mut InputForProcessing) {
    let lighting_efficacy = 120.0;
    input.set_lighting_efficacy_for_all_zones(lighting_efficacy);
}

fn edit_infiltration_ventilation() {
    todo!()
}

/// Apply notional u-value (W/m2K) to:
///
/// external elements: walls (0.18), doors (1.0), roofs (0.11), exposed floors (0.13)
/// elements adjacent to unheated space: walls (0.18), ceilings (0.11), floors (0.13)
/// to differenciate external doors from walls, user input: is_external_door
fn edit_opaque_adjztu_elements(input: &mut InputForProcessing) -> anyhow::Result<()> {
    let mut opaque_adjztu_building_elements =
        input.all_opaque_and_adjztu_building_elements_mut_u_values();

    for building_element in opaque_adjztu_building_elements.iter_mut() {
        let pitch_class = pitch_class(building_element.pitch());
        match pitch_class {
            HeatFlowDirection::Downwards => {
                building_element.set_u_value(0.13);
            }
            HeatFlowDirection::Upwards => {
                building_element.set_u_value(0.11);
            }
            HeatFlowDirection::Horizontal => {
                building_element.set_u_value(0.18);
                // exception if external door
                if building_element.is_opaque() {
                    match building_element.is_external_door() {
                        None => {
                            bail!("Opaque building element input needed value for is_external_door field.");
                        }
                        Some(true) => {
                            building_element.set_u_value(1.0);
                        }
                        _ => {}
                    }
                }
            }
        }
        // remove the r_c input if it was there, as engine would prioritise r_c over u_value
        building_element.remove_r_c();
    }

    Ok(())
}

/// Apply notional u-value to windows & glazed doors and rooflights
/// for windows and glazed doors
/// u-value is 1.2
/// for rooflights
/// u-value is 1.7
/// the max rooflight area is exactly defined as:
/// Max area of glazing if rooflight, as a % of TFA = 25% of TFA - % reduction
/// where % reduction = area of actual rooflight as a % of TFA * ((actual u-value of rooflight - 1.2)/1.2)
/// interpret the instruction for max rooflight area as:
/// max_area_reduction_factor = total_rooflight_area / TFA * ((average_uvalue - 1.2)/1.2)
/// where
///     total_rooflight_area = total area of all rooflights combined
///     average_uvalue = area weighted average actual rooflight u-value
/// max_rooflight_area = maximum allowed total area of all rooflights combined
/// max_rooflight_area = TFA*0.25*max_area_reduction_factor
/// TODO (from Python) - awaiting confirmation from DLUHC/DESNZ that interpretation is correct
fn edit_transparent_element(input: &mut InputForProcessing) -> anyhow::Result<()> {
    let mut _total_rooflight_area = 0.;
    let mut _sum_uval_times_area = 0.;

    let mut building_elements = input.all_transparent_building_elements_mut_u_values();

    for building_element in building_elements.iter_mut() {
        let pitch_class = pitch_class(building_element.pitch());
        match pitch_class {
            HeatFlowDirection::Upwards => {
                // rooflight
                let height = building_element.height().ok_or_else(|| {
                    anyhow!(
                        "FHS notional wrapper needs transparent building elements to have a height set."
                    )
                })?;
                let width = building_element.width().ok_or_else(|| {
                    anyhow!(
                        "FHS notional wrapper needs transparent building elements to have a width set."
                    )
                })?;
                let rooflight_area = height * width;
                _total_rooflight_area += rooflight_area;

                let current_u_value = building_element.u_value().ok_or_else(|| {
                    anyhow!(
                        "FHS notional wrapper needs transparent building elements to have u values set."
                    )
                })?;
                _sum_uval_times_area += current_u_value * rooflight_area;
                building_element.set_u_value(1.7);
                building_element.remove_r_c();
            }
            _ => {
                // if it is not a roof light, it is a glazed door or window
                building_element.set_u_value(1.2);
                building_element.remove_r_c();
            }
        }
    }

    Ok(())
}

///Split windows/rooflights and walls/roofs into dictionaries.
fn split_glazing_and_walls(
    input: &mut InputForProcessing,
) -> (Vec<&mut BuildingElement>, Vec<&mut BuildingElement>) {
    let mut windows_rooflight = vec![];
    let mut walls_roofs = vec![];

    let building_elements = input.all_building_elements_mut();

    for building_element in building_elements {
        match building_element {
            BuildingElement::Transparent { .. } => windows_rooflight.push(building_element),
            BuildingElement::Opaque { .. } => walls_roofs.push(building_element),
            _ => continue,
        }
    }

    (windows_rooflight, walls_roofs)
}

///Calculate difference between old  and new glazing area and adjust the glazing areas
fn calculate_area_diff_and_adjust_glazing_area(
    linear_reduction_factor: f64,
    window_rooflight_element: &mut BuildingElement,
) -> f64 {
    if let BuildingElement::Transparent {
        ref mut height,
        ref mut width,
        ..
    } = window_rooflight_element
    {
        let old_area = *height * *width;
        *height *= linear_reduction_factor;
        *width *= linear_reduction_factor;

        let new_area = *height * *width;

        old_area - new_area
    } else {
        panic!("This function expects to be called for a Transparent BuildingElement only.")
    }
}

///Find all walls/roofs with same orientation and pitch as this window/rooflight.
fn find_walls_roofs_with_same_orientation_and_pitch(
    wall_roofs: &[&mut BuildingElement],
    window_rooflight_element: &BuildingElement,
) -> anyhow::Result<Vec<usize>> {
    let window_rooflight_pitch = window_rooflight_element.pitch();
    let window_rooflight_orientation = window_rooflight_element.orientation();

    let mut indices: Vec<usize> = Default::default();

    for (i, el) in wall_roofs.iter().enumerate() {
        if match el {
            BuildingElement::Opaque {
                pitch, orientation, ..
            } => {
                window_rooflight_orientation.is_some_and(|window_rooflight_orientation| {
                    *orientation == window_rooflight_orientation
                }) && *pitch == window_rooflight_pitch
            }
            _ => false,
        } {
            indices.push(i);
        }
    }

    if indices.is_empty() {
        bail!(
            "There are no walls/roofs with the same orientation and pitch as the window/rooflight"
        );
    }

    Ok(indices)
}

/// Calculate max glazing area fraction for notional building, adjusted for rooflights
fn calc_max_glazing_area_fraction(
    input: &InputForProcessing,
    total_floor_area: f64,
) -> anyhow::Result<f64> {
    let mut total_rooflight_area = 0.0;
    let mut sum_uval_times_area = 0.0;

    for element in input.all_building_elements() {
        if pitch_class(element.pitch()) != HeatFlowDirection::Upwards {
            continue;
        }

        if let BuildingElement::Transparent {
            height,
            width,
            u_value,
            ..
        } = element
        {
            let u_value = u_value.ok_or_else(|| anyhow!("FHS notional wrapper needs transparent building elements to have u values set."))?;
            let rooflight_area = height * width;

            total_rooflight_area += rooflight_area;
            sum_uval_times_area += rooflight_area * u_value;
        }
    }

    let rooflight_correction_factor = if total_rooflight_area == 0.0 {
        0.0
    } else {
        let average_rooflight_uval = sum_uval_times_area / total_rooflight_area;
        let rooflight_proportion = total_rooflight_area / total_floor_area;
        (rooflight_proportion * (average_rooflight_uval - 1.2) / 1.2).max(0.0)
    };

    Ok(0.25 - rooflight_correction_factor)
}

/// Resize window/rooflight and wall/roofs to meet glazing limits
fn edit_glazing_for_glazing_limit(
    input: &mut InputForProcessing,
    total_floor_area: f64,
) -> anyhow::Result<()> {
    let total_glazing_area: f64 = input
        .all_building_elements()
        .iter()
        .filter_map(|el| match el {
            BuildingElement::Transparent { .. } => Some(el.height().unwrap() * el.width().unwrap()),
            _ => None,
        })
        .sum();

    let max_glazing_area_fraction = calc_max_glazing_area_fraction(input, total_floor_area)?;
    let max_glazing_area = max_glazing_area_fraction * total_floor_area;

    let (windows_rooflight, mut walls_roofs) = split_glazing_and_walls(input);

    if total_glazing_area > max_glazing_area {
        let linear_reduction_factor = (max_glazing_area / total_glazing_area).sqrt();

        for window_rooflight_element in windows_rooflight {
            let area_diff = calculate_area_diff_and_adjust_glazing_area(
                linear_reduction_factor,
                window_rooflight_element,
            );

            let same_orientation_indices = find_walls_roofs_with_same_orientation_and_pitch(
                &walls_roofs,
                window_rooflight_element,
            )?;

            let wall_roof_area_total = same_orientation_indices
                .iter()
                .filter_map(|i| match walls_roofs[*i] {
                    BuildingElement::Opaque { ref area, .. } => Some(*area),
                    _ => None,
                })
                .sum::<f64>();

            for i in same_orientation_indices.iter() {
                let wall_roof = walls_roofs.get_mut(*i).unwrap();
                if let BuildingElement::Opaque { ref mut area, .. } = wall_roof {
                    let wall_roof_prop = *area / wall_roof_area_total;

                    *area += area_diff * wall_roof_prop;
                }
            }
        }
    }
    Ok(())
}

/// Apply notional building ground specifications
///
///     u-value = 0.13 W/m2.K
///     thermal resistance of the floor construction,excluding the ground, r_f = 6.12 m2.K/W
///     linear thermal transmittance, psi_wall_floor_junc = 0.16 W/m.K
pub(crate) fn edit_ground_floors(input: &mut InputForProcessing) -> anyhow::Result<()> {
    // TODO (from Python) waiting from MHCLG/DESNZ for clarification if basement floors and basement walls are treated the same

    for building_element in input.all_ground_building_elements_mut() {
        building_element.set_u_value(0.13);
        building_element.set_r_f(6.12);
        building_element.set_psi_wall_floor_junc(0.16);
    }

    Ok(())
}

/// The notional building must follow the same thermal bridges as specified in
/// SAP10.2 Table R2
///
/// TODO (from Python) - how to deal with ThermalBridging when lengths are not specified?
pub(crate) fn edit_thermal_bridging(input: &mut InputForProcessing) -> anyhow::Result<()> {
    let mut thermal_bridging_elements = input.all_thermal_bridging_elements();

    for element in thermal_bridging_elements
        .iter_mut()
        .flat_map(|group| group.values_mut())
    {
        match element {
            ThermalBridgingDetails::Point {
                heat_transfer_coefficient,
                ..
            } => {
                *heat_transfer_coefficient = 0.;
            }
            ThermalBridgingDetails::Linear {
                junction_type,
                linear_thermal_transmittance,
                ..
            } => {
                let junction_type = junction_type.as_ref().and_then(|junc| if TABLE_R2.contains_key(junc.as_str()) {Some(junc.clone())} else {None}).ok_or_else(|| anyhow!("Thermal bridging junction type was expected to be set and one of the values in SAP10.2 Table R2."))?;
                *linear_thermal_transmittance = TABLE_R2[junction_type.as_str()];
            }
        }
    }

    Ok(())
}

/// Table R2 from SAP10.2
static TABLE_R2: LazyLock<HashMap<&'static str, f64>> = LazyLock::new(|| {
    HashMap::from([
        ("E1", 0.05),
        ("E2", 0.05),
        ("E3", 0.05),
        ("E4", 0.05),
        ("E5", 0.16),
        ("E19", 0.07),
        ("E20", 0.32),
        ("E21", 0.32),
        ("E22", 0.07),
        ("E6", 0.),
        ("E7", 0.07),
        ("E8", 0.),
        ("E9", 0.02),
        ("E23", 0.02),
        ("E10", 0.06),
        ("E24", 0.24),
        ("E11", 0.04),
        ("E12", 0.06),
        ("E13", 0.08),
        ("E14", 0.08),
        ("E15", 0.56),
        ("E16", 0.09),
        ("E17", -0.09),
        ("E18", 0.06),
        ("E25", 0.06),
        ("P1", 0.08),
        ("P6", 0.07),
        ("P2", 0.),
        ("P3", 0.),
        ("P7", 0.16),
        ("P8", 0.24),
        ("P4", 0.12),
        ("P5", 0.08),
        ("R1", 0.08),
        ("R2", 0.06),
        ("R3", 0.08),
        ("R4", 0.08),
        ("R5", 0.04),
        ("R6", 0.06),
        ("R7", 0.04),
        ("R8", 0.06),
        ("R9", 0.04),
        ("R10", 0.08),
        ("R11", 0.08),
    ])
});

fn edit_add_heatnetwork_heating() {
    todo!()
}

fn edit_add_default_space_heating_system() {
    todo!()
}

fn edit_default_space_heating_distribution_system() {
    todo!()
}

fn edit_heatnetwork_space_heating_distribution_system() {
    todo!()
}
fn edit_bath_shower_other(
    input: &mut InputForProcessing,
    cold_water_source_type: ColdWaterSourceType,
) -> anyhow::Result<()> {
    // Define Bath, Shower, and Other DHW outlet
    let notional_bath = json!({ NOTIONAL_BATH_NAME: {
            "ColdWaterSource": cold_water_source_type,
            "flowrate": 12,
            "size": 73
        }
    });
    input.set_bath(notional_bath)?;

    let notional_shower = json!({ NOTIONAL_SHOWER_NAME: {
            "ColdWaterSource": cold_water_source_type,
            "flowrate": 8,
            "type": "MixerShower"
        }
    });
    input.set_shower(notional_shower)?;

    let notional_other_hw = json!({ NOTIONAL_OTHER_HW_NAME: {
            "ColdWaterSource": cold_water_source_type,
            "flowrate": 6,
        }
    });
    input.set_other_water_use(notional_other_hw)?;

    Ok(())
}

fn remove_wwhrs_if_present(input: &mut InputForProcessing) {
    if input.wwhrs().is_some() {
        input.remove_wwhrs();
    }
}

fn add_wwhrs(
    input: &mut InputForProcessing,
    cold_water_source_type: ColdWaterSourceType,
    is_notional_a: bool,
    is_fee: bool,
) -> anyhow::Result<()> {
    // TODO (from Python) Storeys in dwelling is not currently collected as an input, so use
    //      storeys in building for houses and assume 1 for flats. Note that this
    //      means that maisonettes cannot be handled at present.

    let storeys_in_building = match input.build_type() {
        BuildType::House => input.storeys_in_building(),
        BuildType::Flat => 1,
    };

    // add WWHRS if more than 1 storeys in dwelling, notional A and not FEE
    if storeys_in_building > 1 && is_notional_a && !is_fee {
        input.register_wwhrs_name_on_mixer_shower(NOTIONAL_WWHRS)?;
        input.set_wwhrs(json!({
            NOTIONAL_WWHRS: {
                "ColdWaterSource": cold_water_source_type,
                "efficiencies": [50, 50],
                "flow_rates": [0, 100],
                "type": "WWHRS_InstantaneousSystemB",
                "utilisation_factor": 0.98
            }
        }))?;
    }

    Ok(())
}

fn calculate_daily_losses(cylinder_vol: f64) -> f64 {
    const CYLINDER_LOSS: f64 = 0.005;
    const FACTORY_INSULATED_THICKNESS_COEFF: f64 = 0.55;
    const THICKNESS: f64 = 120.; // mm

    // calculate cylinder factor insulated factor
    let cylinder_heat_loss_factor =
        CYLINDER_LOSS + FACTORY_INSULATED_THICKNESS_COEFF / (THICKNESS + 4.0);

    // calculate volume factor
    let vol_factor = (120. / cylinder_vol).powf(1. / 3.);

    // Temperature factor
    let temp_factor = 0.6 * 0.9;

    // Calculate daily losses
    cylinder_heat_loss_factor * vol_factor * temp_factor * cylinder_vol
}

fn calc_daily_hw_demand(
    input: &mut InputForProcessing,
    total_floor_area: f64,
    cold_water_source_type: ColdWaterSourceType,
) -> anyhow::Result<Vec<f64>> {
    // create SimulationTime
    let simtime = SimulationTime::new(SIMTIME_START, SIMTIME_END, SIMTIME_STEP);

    // create ColdWaterSource
    let cold_water_feed_temps = create_cold_water_feed_temps(input)?;
    let cold_water_sources: ColdWaterSources = input
        .cold_water_source()
        .iter()
        .map(|(key, source)| {
            (
                *key,
                ColdWaterSource::new(
                    source.temperatures.clone(),
                    &simtime,
                    source.start_day,
                    source.time_series_step,
                ),
            )
        })
        .collect();

    let wwhrs: IndexMap<String, Arc<Mutex<Wwhrs>>> = if let Some(waste_water_heat_recovery) =
        input.wwhrs()
    {
        let notional_wwhrs = waste_water_heat_recovery.get(NOTIONAL_WWHRS).ok_or_else(|| anyhow!("A {} entry for WWHRS was expected to have been set in the FHS Notional wrapper.", NOTIONAL_WWHRS))?;
        [(
            NOTIONAL_WWHRS.to_owned(),
            Arc::new(Mutex::new(Wwhrs::WWHRSInstantaneousSystemB(
                WWHRSInstantaneousSystemB::new(
                    cold_water_sources
                        .get(&notional_wwhrs.cold_water_source)
                        .ok_or_else(|| {
                            anyhow!(
                                "A cold water source could not be found with the type '{:?}'.",
                                notional_wwhrs.cold_water_source
                            )
                        })?
                        .clone(),
                    notional_wwhrs.flow_rates.clone(),
                    notional_wwhrs.efficiencies.clone(),
                    notional_wwhrs.utilisation_factor,
                ),
            ))),
        )]
        .into()
    } else {
        Default::default()
    };

    let nbeds = calc_nbeds(input)?;
    let number_of_occupants = calc_n_occupants(total_floor_area, nbeds)?;
    create_hot_water_use_pattern(input, number_of_occupants, &cold_water_feed_temps)?;
    let sim_timestep = simtime.step;
    let total_timesteps = simtime.total_steps();
    let event_types_names_list = [
        (WaterHeatingEventType::Shower, NOTIONAL_SHOWER_NAME),
        (WaterHeatingEventType::Bath, NOTIONAL_BATH_NAME),
        (WaterHeatingEventType::Other, NOTIONAL_OTHER_HW_NAME),
    ];

    // Initialize a single schedule dictionary
    let mut event_schedules: Vec<Option<Vec<TypedScheduleEvent>>> = vec![None; total_timesteps];

    // Populate the event_schedules dictionary using the modified expand_events function
    for (event_type, event_name) in event_types_names_list {
        let event_data = input.water_heating_event_by_type_and_name(event_type, event_name).ok_or_else(|| anyhow!("FHS Notional wrapper expected water heating events with type '{event_type:?}' and name '{event_name}'"))?.iter().map(Into::into).collect::<Vec<_>>();
        event_schedules = expand_events(
            event_data,
            sim_timestep,
            total_timesteps,
            event_name,
            event_type.into(),
            event_schedules,
        )?;
    }

    let dhw_demand = DomesticHotWaterDemand::new(
        input.showers().cloned().unwrap_or_default(),
        input.baths().cloned().unwrap_or_default(),
        input.other_water_uses().cloned().unwrap_or_default(),
        input.water_distribution().cloned(),
        &cold_water_sources,
        &wwhrs,
        &EnergySupplies::default(),
        event_schedules,
    )?;

    // For each timestep, calculate HW draw
    let total_steps = simtime.total_steps();
    let mut hw_energy_demand = vec![0.0; total_steps];
    for (t_idx, t_it) in simtime.iter().enumerate() {
        let DomesticHotWaterDemandData { hw_demand_vol, .. } =
            dhw_demand.hot_water_demand(t_it, HW_TEMPERATURE);

        // Convert from litres to kWh
        let cold_water_temperature = cold_water_sources[&cold_water_source_type].temperature(t_it);
        hw_energy_demand[t_idx] =
            water_demand_to_kwh(hw_demand_vol, HW_TEMPERATURE, cold_water_temperature);
    }

    Ok(convert_profile_to_daily(&hw_energy_demand, simtime.step))
}

fn edit_storagetank(
    input: &mut InputForProcessing,
    cold_water_source_type: ColdWaterSourceType,
    total_floor_area: f64,
) -> anyhow::Result<()> {
    let cylinder_vol = match input.hot_water_cylinder_volume() {
        Some(volume) => volume,
        None => {
            let daily_hwd = calc_daily_hw_demand(input, total_floor_area, cold_water_source_type)?;
            calculate_cylinder_volume(&daily_hwd)
        }
    };

    // Calculate daily losses
    let daily_losses = calculate_daily_losses(cylinder_vol);

    // Modify primary pipework characteristics
    let primary_pipework = edit_primary_pipework(input, total_floor_area)?;

    // Modify cylinder characteristics
    input.set_hot_water_cylinder(json!({
        "ColdWaterSource": cold_water_source_type,
            "HeatSource": {
                NOTIONAL_HP: {
                    "ColdWaterSource": cold_water_source_type,
                    "EnergySupply": "mains elec",
                    "heater_position": 0.1,
                    "name": NOTIONAL_HP,
                    "temp_flow_limit_upper": 60,
                    "thermostat_position": 0.1,
                    "type": "HeatSourceWet"
                }
            },
            "daily_losses": daily_losses,
            "type": "StorageTank",
            "volume": cylinder_vol,
            "primary_pipework": primary_pipework
    }))?;

    Ok(())
}

fn edit_primary_pipework(
    input: &InputForProcessing,
    total_floor_area: f64,
) -> anyhow::Result<Vec<WaterPipework>> {
    // Define minimum values
    let internal_diameter_mm_min = 20.;
    let external_diameter_mm_min = 22.;
    let insulation_thickness_mm_min = 25.;
    let surface_reflectivity = false;
    let pipe_contents = WaterPipeContentsType::Water;
    let insulation_thermal_conductivity = 0.035;

    let length_max = match input.build_type() {
        BuildType::Flat => 0.05 * total_floor_area,
        BuildType::House => {
            0.05 * input.ground_floor_area().ok_or_else(|| {
                anyhow!("FHS Notional wrapper expected ground floor area to be set for a house.")
            })?
        }
    };

    let mut primary_pipework = input.primary_pipework_clone();

    if primary_pipework.is_none() {
        primary_pipework = Some(vec![serde_json::from_value::<WaterPipework>(json!({
            "location": "internal",
            "internal_diameter_mm": internal_diameter_mm_min,
            "external_diameter_mm": external_diameter_mm_min,
            "length": length_max,
            "insulation_thermal_conductivity": insulation_thermal_conductivity,
            "insulation_thickness_mm": insulation_thickness_mm_min,
            "surface_reflectivity": surface_reflectivity,
            "pipe_contents": pipe_contents
        }))?]);
    } else {
        for pipework in primary_pipework.as_mut().unwrap().iter_mut() {
            let length = pipework.length;
            let internal_diameter_mm = pipework.internal_diameter_mm.max(internal_diameter_mm_min);
            let external_diameter_mm = pipework.external_diameter_mm.max(external_diameter_mm_min);

            // Update insulation thickness based on internal diameter
            let adjusted_insulation_thickness_mm_min = if internal_diameter_mm > 25. {
                35.
            } else {
                insulation_thickness_mm_min
            };

            // Primary pipework should not be greater than maximum length
            let length = length.min(length_max);

            // Update pipework
            *pipework = serde_json::from_value(json!({
                "location": "internal",
                "internal_diameter_mm": internal_diameter_mm,
                "external_diameter_mm": external_diameter_mm,
                "length": length,
                "insulation_thermal_conductivity": insulation_thermal_conductivity,
                "insulation_thickness_mm": adjusted_insulation_thickness_mm_min,
                "surface_reflectivity": surface_reflectivity,
                "pipe_contents": pipe_contents
            }))?;
        }
    }

    Ok(primary_pipework.unwrap())
}

fn edit_hot_water_distribution(
    input: &mut InputForProcessing,
    total_floor_area: f64,
) -> anyhow::Result<()> {
    // hot water dictionary
    let mut hot_water_distribution_inner_list = vec![];

    for item in input.water_distribution().into_iter().flatten() {
        // only include internal pipework in notional buildings
        if item.location == WaterPipeworkLocation::Internal {
            hot_water_distribution_inner_list.push(item);
        }
    }

    // Create an empty list to store updated dictionaries
    let mut updated_hot_water_distribution_inner_list = vec![];

    // Defaults
    let internal_diameter_mm_min = 13.;
    let external_diameter_mm_min = 15.;
    let insulation_thickness_mm = 20.;

    let length_max = match input.build_type() {
        BuildType::Flat => 0.2 * total_floor_area,
        BuildType::House => {
            0.2 * input
                .ground_floor_area()
                .ok_or_else(|| anyhow!("Notional wrapper expected a ground floor area to be set"))?
        }
    };

    // Iterate over hot_water_distribution_inner_list
    for hot_water_distribution_inner in hot_water_distribution_inner_list {
        // hot water distribution (inner) length should not be greater than maximum length

        let length_actual = hot_water_distribution_inner.length;
        let length = min_of_2(length_actual, length_max);

        // Update internal diameter to minimum if not present and should not be lower than the minimum
        let internal_diameter_mm = hot_water_distribution_inner
            .internal_diameter_mm
            .unwrap_or(internal_diameter_mm_min);
        let internal_diameter_mm = internal_diameter_mm.max(internal_diameter_mm_min);

        // Update external diameter to minimum if not present and should not be lower than the minimum
        let external_diameter_mm = hot_water_distribution_inner
            .external_diameter_mm
            .unwrap_or(external_diameter_mm_min);
        let external_diameter_mm = external_diameter_mm.max(external_diameter_mm_min);

        // Update insulation thickness based on internal diameter
        let adjusted_insulation_thickness_mm = if internal_diameter_mm > 25. {
            24.
        } else {
            insulation_thickness_mm
        };

        let pipework_to_update = json!({
            "location": "internal",
            "external_diameter_mm": external_diameter_mm,
            "insulation_thermal_conductivity": 0.035,
            "insulation_thickness_mm": adjusted_insulation_thickness_mm,
            "internal_diameter_mm": internal_diameter_mm,
            "length": length,
            "pipe_contents": "water",
            "surface_reflectivity": false
        });

        updated_hot_water_distribution_inner_list.push(pipework_to_update);
    }

    input.set_water_distribution(serde_json::Value::Array(
        updated_hot_water_distribution_inner_list,
    ))?;

    Ok(())
}

fn remove_pv_diverter_if_present(input: &mut InputForProcessing) {
    input.remove_all_diverters_from_energy_supplies();
}

fn remove_electric_battery_if_present(input: &mut InputForProcessing) {
    input.remove_all_batteries_from_energy_supplies();
}

/// Calculate effective air change rate according to Part F 1.24 a
pub fn minimum_air_change_rate(
    _input: &InputForProcessing,
    total_floor_area: f64,
    total_volume: f64,
    bedroom_number: usize,
) -> f64 {
    // minimum ventilation rates method B
    let min_ventilation_rates_b = [19, 25, 31, 37, 43];

    // Calculate minimum whole dwelling ventilation rate l/s method A
    let min_ventilation_rate_a = total_floor_area * 0.3;

    // Calculate minimum whole dwelling ventilation rate l/s method B
    let min_ventilation_rate_b = if bedroom_number <= 5 {
        min_ventilation_rates_b[bedroom_number - 1]
    } else {
        min_ventilation_rates_b.last().unwrap() + (bedroom_number - 5) * 6
    };

    // Calculate air change rate ACH
    (max_of_2(min_ventilation_rate_a, min_ventilation_rate_b as f64) / total_volume)
        * SECONDS_PER_HOUR as f64
        / LITRES_PER_CUBIC_METRE as f64
}

fn edit_space_heating_system() {
    todo!()
}

fn edit_space_cool_system(input: &mut InputForProcessing) -> anyhow::Result<()> {
    let part_o_active_cooling_required = input.part_o_active_cooling_required().unwrap_or(false);

    if part_o_active_cooling_required {
        input.set_efficiency_for_all_space_cool_systems(5.1)?;
        input.set_frac_convective_for_all_space_cool_systems(0.95)?;
        input.set_energy_supply_for_all_space_cool_systems(EnergySupplyType::Electricity)?;
    }

    Ok(())
}

fn calc_design_capacity() {
    todo!()
}

fn initialise_temperature_setpoints() {
    todo!()
}

fn remove_onsite_generation_if_present(input: &mut InputForProcessing) {
    if input.on_site_generation().is_some() {
        input.remove_on_site_generation();
    }
}

fn add_solar_pv(
    input: &mut InputForProcessing,
    is_notional_a: bool,
    is_fee: bool,
    total_floor_area: f64,
) -> anyhow::Result<()> {
    let number_of_storeys = input.storeys_in_building();

    // PV is included in the notional if the building contains 15 stories or
    // less that contain dwellings.
    if number_of_storeys <= 15 && is_notional_a && !is_fee {
        let ground_floor_area = input
            .ground_floor_area()
            .ok_or_else(|| anyhow!("Notional wrapped expected ground floor area to be set"))?;
        let (peak_kw, base_height_pv) = match input.build_type() {
            BuildType::House => {
                let peak_kw = ground_floor_area * 0.4 / 4.5;
                let base_height_pv = input.max_base_height_from_building_elements().ok_or_else(|| anyhow!("Notional wrapper expected at least one building element with a base height"))?;

                (peak_kw, base_height_pv)
            }
            BuildType::Flat => {
                let peak_kw = total_floor_area * 0.4 / (4.5 * number_of_storeys as f64);
                let zone_total_volume = input.total_zone_volume();
                let zone_total_area = input.total_zone_area();
                let base_height_pv =
                    (zone_total_volume / zone_total_area + 0.3) * number_of_storeys as f64;

                (peak_kw, base_height_pv)
            }
        };

        // PV array area
        let pv_area = 4.5 * peak_kw;

        // PV width and height based on 2:1 aspect ratio
        let pv_height = (pv_area / 2.).powf(0.5);
        let pv_width = 2. * pv_height;

        let solar_pv = json!({
            "PV1": {
                "EnergySupply": "mains elec",
                "orientation360": 180.,
                "peak_power": peak_kw,
                "inverter_peak_power": peak_kw,
                "inverter_is_inside": false,
                "pitch": 45.,
                "type": "PhotovoltaicSystem",
                "ventilation_strategy": "moderately_ventilated",
                "base_height": base_height_pv,
                "height":pv_height,
                "width":pv_width,
                "shading": [],
                }
        });

        input.set_on_site_generation(solar_pv)?;
    }

    Ok(())
}

fn calculate_cylinder_volume(daily_hwd: &[f64]) -> f64 {
    // Data from the table
    let percentiles_kwh = [3.7, 4.4, 5.2, 5.9, 6.7, 7.4, 8.1, 8.9, 9.6, 10.3, 11.1];
    let vessel_sizes_litres = [
        165., 190., 215., 240., 265., 290., 315., 340., 365., 390., 415.,
    ];

    // Calculate the 75th percentile of daily hot water demand
    let percentile_75_kwh = percentile(daily_hwd, 75);

    // Use linear interpolation to find the appropriate vessel size
    let interpolated_size_litres =
        np_interp(percentile_75_kwh, &percentiles_kwh, &vessel_sizes_litres);
    let mut interpolated_size_litres = interpolated_size_litres.round();

    // If the size of the hot water storage vessel is unavailable, the next
    // largest size available should be selected
    if !vessel_sizes_litres.contains(&interpolated_size_litres) {
        for size in vessel_sizes_litres {
            if size > interpolated_size_litres {
                interpolated_size_litres = size;
                break;
            }
        }
    }

    interpolated_size_litres
}

#[cfg(test)]
mod tests {
    use crate::core::space_heat_demand::building_element::{pitch_class, HeatFlowDirection};

    use super::*;
    use crate::input::{
        self, EnergySupplyKey, EnergySupplyType, OnSiteGeneration, WaterPipeworkSimple,
    };
    use crate::input::{
        Baths, HotWaterSource, OtherWaterUses, Shower, Showers, ThermalBridging,
        ThermalBridgingDetails, WasteWaterHeatRecovery, ZoneDictionary,
    };
    use approx::assert_relative_eq;
    use rstest::{fixture, rstest};
    use serde_json::json;
    use std::borrow::BorrowMut;
    use std::io::{BufReader, Cursor};

    #[fixture]
    fn test_input() -> InputForProcessing {
        let reader = BufReader::new(Cursor::new(include_str!(
            "./test_future_homes_standard_notional_input_data.json"
        )));
        InputForProcessing::init_with_json(reader).expect(
            "expected valid test_future_homes_standard_notional_input_data.json to be present",
        )
    }

    #[ignore = "still to complete"]
    #[rstest]
    // this test does not exist in Python HEM
    fn test_apply_fhs_notional_preprocessing(mut test_input: InputForProcessing) {
        let fhs_notional_a_assumptions = true;
        let _fhs_notional_b_assumptions = false;
        let fhs_fee_notional_a_assumptions = false;
        let fhs_fee_notional_b_assumptions = false;

        let _actual = apply_fhs_notional_preprocessing(
            &mut test_input,
            fhs_notional_a_assumptions,
            _fhs_notional_b_assumptions,
            fhs_fee_notional_a_assumptions,
            fhs_fee_notional_b_assumptions,
        );
        todo!()
    }

    #[rstest]
    // this test does not exist in Python HEM
    fn test_check_heatnetwork_present(test_input: InputForProcessing) {
        assert!(!check_heatnetwork_present(&test_input));
    }

    #[rstest]
    fn test_edit_lighting_efficacy(mut test_input: InputForProcessing) {
        let test_input = test_input.borrow_mut();
        edit_lighting_efficacy(test_input);

        for zone in test_input.zone_keys() {
            let lighting_efficacy = test_input.lighting_efficacy_for_zone(&zone);
            assert_eq!(
                lighting_efficacy
                    .unwrap()
                    .expect("expected lighting in zone and efficacy in lighting"),
                120.
            )
        }
    }

    #[rstest]
    fn test_edit_opaque_ajdztu_elements(mut test_input: InputForProcessing) {
        edit_opaque_adjztu_elements(&mut test_input).unwrap();

        // not using the building_element_by_key method here to closly match the Python test

        for building_element in test_input.all_building_elements() {
            if let input::BuildingElement::Opaque { .. }
            | input::BuildingElement::AdjacentZTUSimple { .. } = building_element
            {
                if let Some(u_value) = building_element.u_value() {
                    match pitch_class(building_element.pitch()) {
                        HeatFlowDirection::Downwards => {
                            // this assertion is not reached with the current test data
                            assert_eq!(u_value, 0.13);
                        }
                        HeatFlowDirection::Upwards => {
                            assert_eq!(u_value, 0.11);
                        }
                        HeatFlowDirection::Horizontal => {
                            let expected_u_value = if let input::BuildingElement::Opaque {
                                is_external_door: Some(true),
                                ..
                            } = building_element
                            {
                                1.0
                            } else {
                                0.18
                            };
                            assert_eq!(u_value, expected_u_value);
                        }
                    }
                }
            }
        }
    }

    // this test does not exist in Python HEM
    #[rstest]
    fn test_edit_transparent_element(mut test_input: InputForProcessing) {
        edit_transparent_element(&mut test_input).unwrap();

        let BuildingElement::Transparent { u_value, r_c, .. } =
            test_input.building_element_by_key("zone 2", "skylight 0")
        else {
            panic!("Skylight 0 in Zone 2 should be set up as a transparent building element")
        };

        assert_eq!(*u_value, Some(1.7));
        assert_eq!(*r_c, None);

        let BuildingElement::Transparent { u_value, r_c, .. } =
            test_input.building_element_by_key("zone 1", "window 0")
        else {
            panic!("Window 0 in Zone 1 should be set up as a transparent building element")
        };

        assert_eq!(*u_value, Some(1.2));
        assert_eq!(*r_c, None);

        let BuildingElement::Transparent { u_value, r_c, .. } =
            test_input.building_element_by_key("zone 2", "window 0")
        else {
            panic!("Window 0 in Zone 2 should be set up as a transparent building element")
        };

        assert_eq!(*u_value, Some(1.2));
        assert_eq!(*r_c, None);
    }

    #[rstest]
    fn test_edit_ground_floors(mut test_input: InputForProcessing) {
        let test_input = test_input.borrow_mut();

        edit_ground_floors(test_input).unwrap();

        let BuildingElement::Ground {
            u_value,
            r_f,
            psi_wall_floor_junc,
            ..
        } = test_input.building_element_by_key("zone 1", "ground")
        else {
            panic!("ground in Zone 1 should be set up as a ground building element")
        };
        assert_eq!(*u_value, 0.13);
        assert_eq!(*r_f, 6.12);
        assert_eq!(*psi_wall_floor_junc, 0.16);

        let BuildingElement::Ground {
            u_value,
            r_f,
            psi_wall_floor_junc,
            ..
        } = test_input.building_element_by_key("zone 2", "ground")
        else {
            panic!("ground in Zone 2 should be set up as a ground building element")
        };
        assert_eq!(*u_value, 0.13);
        assert_eq!(*r_f, 6.12);
        assert_eq!(*psi_wall_floor_junc, 0.16);
    }

    #[rstest]
    fn test_edit_thermal_bridgings(mut test_input: InputForProcessing) {
        let test_input = test_input.borrow_mut();

        edit_thermal_bridging(test_input).unwrap();

        for thermal_bridging in test_input.all_thermal_bridgings() {
            match thermal_bridging {
                ThermalBridging::ThermalBridgingNumber(_) => {}
                ThermalBridging::ThermalBridgingElements(elements) => {
                    for bridging in elements.values() {
                        match bridging {
                            ThermalBridgingDetails::Point {
                                heat_transfer_coefficient,
                                ..
                            } => {
                                assert_eq!(*heat_transfer_coefficient, 0.0);
                            }
                            ThermalBridgingDetails::Linear {
                                junction_type,
                                linear_thermal_transmittance,
                                ..
                            } => {
                                let junction_type = junction_type.as_ref().unwrap().as_str();
                                assert!(TABLE_R2.contains_key(junction_type));
                                assert_eq!(*linear_thermal_transmittance, TABLE_R2[junction_type]);
                            }
                        }
                    }
                }
            }
        }
    }

    #[rstest]
    fn test_calc_max_glazing_area_fraction(mut test_input: InputForProcessing) {
        test_input.set_zone(zone_input_for_max_glazing_area_test(1.5, None));
        assert_eq!(
            calc_max_glazing_area_fraction(&test_input, 80.0).unwrap(),
            0.24375,
            "incorrect max glazing area fraction"
        );
        test_input.set_zone(zone_input_for_max_glazing_area_test(1.0, None));
        assert_eq!(
            calc_max_glazing_area_fraction(&test_input, 80.0).unwrap(),
            0.25,
            "incorrect max glazing area fraction"
        );
        test_input.set_zone(zone_input_for_max_glazing_area_test(1.5, Some(90.)));
        assert_eq!(
            calc_max_glazing_area_fraction(&test_input, 80.0).unwrap(),
            0.25,
            "incorrect max glazing area fraction when there are no rooflights"
        );
    }

    // this test does not exist in Python HEM
    #[rstest]
    fn test_edit_glazing_for_glazing_limit(mut test_input: InputForProcessing) {
        let total_floor_area = 1000.;
        edit_glazing_for_glazing_limit(&mut test_input, total_floor_area).unwrap();

        let skylight = test_input.building_element_by_key("zone 2", "skylight 0");
        let roof = test_input.building_element_by_key("zone 2", "roof 0");

        assert_relative_eq!(skylight.width().unwrap(), 280.056016805602);
        assert_relative_eq!(skylight.height().unwrap(), 0.8751750525175062);

        if let BuildingElement::Opaque { area, .. } = roof {
            assert_relative_eq!(*area, 269.9019607843137);
        } else {
            unreachable!()
        }
    }

    // this test does not exist in Python HEM
    #[rstest]
    fn test_calculate_area_diff_and_adjust_glazing_area(mut test_input: InputForProcessing) {
        let linear_reduction_factor: f64 = 0.7001400420140049;

        let window = test_input.building_element_by_key_mut("zone 1", "window 0");

        let area_diff =
            calculate_area_diff_and_adjust_glazing_area(linear_reduction_factor, window);

        assert_relative_eq!(area_diff, 2.549019607843137);
    }

    fn zone_input_for_max_glazing_area_test(
        u_value: f64,
        pitch_override: Option<f64>,
    ) -> ZoneDictionary {
        serde_json::from_value(json!({"test_zone": {
            "BuildingElement": {
                "test_rooflight": {
                    "type": "BuildingElementTransparent",
                    "pitch": pitch_override.unwrap_or(0.0),
                    "height": 2.0,
                    "width": 1.0,
                    "u_value": u_value,
                    "orientation360": 0.0,
                    "g_value": 0.0,
                    "frame_area_fraction": 0.0,
                    "base_height": 0.0,
                    "shading": [],
                }
            },
            "ThermalBridging": 1.0,
            "area": 10.0,
            "volume": 20.0
        }}))
        .unwrap()
    }

    #[rstest]
    fn test_edit_bath_shower_other(mut test_input: InputForProcessing) {
        // this is the only cold water source type in the test input JSON file
        let cold_water_source_type = ColdWaterSourceType::MainsWater;
        let cold_water_source_type_string = "mains water";

        edit_bath_shower_other(&mut test_input, cold_water_source_type).unwrap();

        let expected_baths: Baths = serde_json::from_value(json!({ "medium": {
            "ColdWaterSource": cold_water_source_type_string,
            "flowrate": 12,
            "size": 73
        }}))
        .unwrap();

        let expected_showers: Showers = serde_json::from_value(json!({"mixer": {
            "ColdWaterSource": cold_water_source_type_string,
            "flowrate": 8,
            "type": "MixerShower"
        }}))
        .unwrap();

        let expected_other: OtherWaterUses = serde_json::from_value(json!({"other": {
            "ColdWaterSource": cold_water_source_type_string,
            "flowrate": 6,
        }}))
        .unwrap();

        assert_eq!(test_input.baths().unwrap().clone(), expected_baths);
        assert_eq!(test_input.showers().unwrap().clone(), expected_showers);
        assert_eq!(
            test_input.other_water_uses().unwrap().clone(),
            expected_other
        );
    }

    // this test does not exist in Python HEM
    #[rstest]
    fn test_remove_wwhrs_if_present(mut test_input: InputForProcessing) {
        test_input
            .register_wwhrs_name_on_mixer_shower(NOTIONAL_WWHRS)
            .unwrap();
        test_input
            .set_wwhrs(json!({
                NOTIONAL_WWHRS: {
                    "ColdWaterSource": ColdWaterSourceType::MainsWater,
                    "efficiencies": [50, 50],
                    "flow_rates": [0, 100],
                    "type": "WWHRS_InstantaneousSystemB",
                    "utilisation_factor": 0.98
                }
            }))
            .unwrap();
        remove_wwhrs_if_present(&mut test_input);
        assert_eq!(test_input.wwhrs(), None);
    }

    #[rstest]
    fn test_add_wwhrs(mut test_input: InputForProcessing) {
        test_input.set_storeys_in_building(2);

        let cold_water_source_type = ColdWaterSourceType::MainsWater;

        add_wwhrs(&mut test_input, cold_water_source_type, true, false).unwrap();

        let expected_wwhrs: WasteWaterHeatRecovery = serde_json::from_value(json!({
            "Notional_Inst_WWHRS": {
               "ColdWaterSource": cold_water_source_type,
               "efficiencies": [50, 50],
               "flow_rates": [0, 100],
               "type": "WWHRS_InstantaneousSystemB",
               "utilisation_factor": 0.98
           }
        }))
        .unwrap();

        assert_eq!(test_input.wwhrs().unwrap().clone(), expected_wwhrs);

        match &test_input.showers().as_ref().unwrap().0["mixer"] {
            Shower::MixerShower {
                waste_water_heat_recovery,
                ..
            } => {
                assert_eq!(
                    waste_water_heat_recovery.as_ref().unwrap(),
                    "Notional_Inst_WWHRS"
                )
            }
            _ => {
                panic!()
            }
        }
    }

    #[rstest]
    fn test_add_no_wwhrs_for_one_storey_buildings(mut test_input: InputForProcessing) {
        let cold_water_source_type = ColdWaterSourceType::MainsWater;

        add_wwhrs(&mut test_input, cold_water_source_type, true, false).unwrap();

        assert!(test_input.wwhrs().is_none());

        match &test_input.showers().as_ref().unwrap().0["mixer"] {
            Shower::MixerShower {
                waste_water_heat_recovery,
                ..
            } => assert!(waste_water_heat_recovery.is_none()),
            _ => {
                panic!()
            }
        }
    }

    #[rstest]
    fn test_calculate_daily_losses() {
        let cylinder_vol = 265.;
        let actual_daily_losses = calculate_daily_losses(cylinder_vol);
        let expected_daily_losses = 1.03685;
        assert_relative_eq!(
            actual_daily_losses,
            expected_daily_losses,
            max_relative = 1E-6
        );
    }

    #[rstest]
    fn test_edit_storagetank(mut test_input: InputForProcessing) {
        let cold_water_source_type = ColdWaterSourceType::MainsWater;
        let total_floor_area = calc_tfa(&test_input);

        edit_storagetank(&mut test_input, cold_water_source_type, total_floor_area).unwrap();

        let expected_primary_pipework = json!([{
            "location": "internal",
            "internal_diameter_mm": 26,
                "external_diameter_mm": 28,
                "length": 2.5,
                "insulation_thermal_conductivity": 0.035,
                "insulation_thickness_mm": 35,
                "surface_reflectivity": false,
                "pipe_contents": "water"
        }]);

        let expected_hotwater_source: HotWaterSource = serde_json::from_value(json!({
            "hw cylinder": {
                "ColdWaterSource": cold_water_source_type,
                "HeatSource": {
                    "notional_HP": {
                        "ColdWaterSource": cold_water_source_type,
                        "EnergySupply": "mains elec",
                        "heater_position": 0.1,
                        "name": "notional_HP",
                        "temp_flow_limit_upper": 60,
                        "thermostat_position": 0.1,
                        "type": "HeatSourceWet"
                    }
                },
                "daily_losses": 0.46660029577109363,
                "type": "StorageTank",
                "volume": 80.0,
                "primary_pipework": expected_primary_pipework,
            }
        }))
        .unwrap();

        assert_eq!(test_input.hot_water_source(), &expected_hotwater_source);
    }

    // Below test fails due to randomisation issue
    #[ignore]
    #[rstest]
    fn test_calc_daily_hw_demand(mut test_input: InputForProcessing) {
        let cold_water_source_type = ColdWaterSourceType::MainsWater;
        let total_floor_area = calc_tfa(&test_input);

        // Add notional objects that affect HW demand calc
        edit_bath_shower_other(&mut test_input, cold_water_source_type).unwrap();

        let daily_hwd =
            calc_daily_hw_demand(&mut test_input, total_floor_area, cold_water_source_type);

        let expected_result = [
            1.9247705151905075,
            5.576203038843596,
            7.6007157714945075,
            7.449862571901721,
            7.461412530591756,
            3.7286348666793794,
            2.3556348486180076,
            9.244086081999619,
            7.809417885183313,
            6.720229563250378,
            9.309311605930931,
            6.4626685134725,
            4.926337875707286,
            2.955998159025295,
            4.378691240445834,
            4.183481647765038,
            2.8501758203228684,
            3.9789394569708225,
            3.6457378380926326,
            3.5542994822205465,
            7.5509033652933395,
            8.588631512469588,
            3.3982177377487055,
            9.837603310451735,
            1.7405253832752825,
            5.647133414069438,
            5.267627347117556,
            4.695838488918805,
            9.155438749768443,
            4.288055751730484,
            5.621785072585681,
            1.8740633510531937,
            4.331598555562393,
            1.826166638255688,
            3.68962861221586,
            2.219500171925045,
            1.7967277069158922,
            4.317531786893869,
            2.57706747992738,
            12.306282500041988,
            2.8753642751239674,
            3.7600221662740894,
            7.360079797196356,
            2.6296124207791345,
            7.807057152573838,
            4.43609522704038,
            1.9156220967756887,
            5.284452049284063,
            2.4949323051085392,
            6.497533306250856,
            3.824180666742241,
            9.728040912664783,
            1.814920720580427,
            5.55633979352759,
            2.6953349766591557,
            5.095593369748756,
            8.561784374656341,
            8.156292707910048,
            8.475773955460514,
            5.0216984064446315,
            3.3270019319900848,
            9.907481952712754,
            4.543546158764719,
            2.9445325844941674,
            4.288447997961151,
            2.1088537397078175,
            5.750511597879198,
            2.5039708929181157,
            8.375000524872696,
            2.7868751749184786,
            3.511742128969035,
            5.85703895590205,
            2.105618084694723,
            8.036785305151378,
            3.0368154032538444,
            1.6594446006198589,
            4.613839563966503,
            2.3055377201064493,
            5.658198924950215,
            5.407084330954116,
            4.542787190556694,
            3.270221582578504,
            5.987669970564481,
            4.42558556578354,
            8.714265525023105,
            4.5049508141313686,
            5.483469087886154,
            3.9045131886808133,
            6.794033200950009,
            6.047933970327377,
            5.795256117004792,
            4.3744770522502625,
            1.7876728415727776,
            2.9060549429238627,
            3.7227364706001675,
            3.009488393003571,
            2.759132853869388,
            6.335308821993266,
            3.9259359337662816,
            4.383460811361665,
            6.6585522587615475,
            2.081858784992064,
            2.1247837102426796,
            6.227994536424524,
            2.433778062146155,
            1.908817497769811,
            3.5887148755723746,
            7.636320821287575,
            3.618316428089196,
            5.355530952160414,
            2.9681786058728816,
            2.1597750260596595,
            3.040417418105342,
            4.609059851140694,
            3.583962229108492,
            3.3040593630692907,
            6.3482511235842916,
            2.613016804330107,
            6.909558537262041,
            5.800486476904296,
            3.8753800442292112,
            1.2409240762417555,
            4.029465630663848,
            1.6107646096574633,
            4.631302847695539,
            8.855234198347574,
            4.197124682287838,
            2.5210598417535834,
            1.6581406681030086,
            5.92361889002503,
            5.8780307729859045,
            3.9856748937058715,
            4.271843630442291,
            4.253576957532686,
            7.591420868012066,
            2.3114540617002106,
            3.017444562575851,
            5.904036521258879,
            2.204938284391645,
            2.2754858559645026,
            3.6840009029293195,
            3.089890246704239,
            4.253465238337842,
            2.8324597660832893,
            3.1612981304237384,
            8.325633143994937,
            1.8553836461789126,
            2.2969172668319158,
            2.397706693159479,
            3.609517691408761,
            10.08182682035413,
            1.8217819704354354,
            3.3592832481688077,
            2.2211629069161405,
            5.567710458036769,
            1.3306568484728694,
            5.901484368346956,
            3.561281306410421,
            3.8332381578251282,
            4.94069039508088,
            4.08203202380226,
            1.9418831608216365,
            1.2464840645384074,
            2.6492862159406174,
            1.7845473877485638,
            5.166469854690617,
            3.1823256992387385,
            5.092541331537495,
            3.6610555908624582,
            3.789039816649784,
            4.723335268842766,
            1.6700781436939853,
            1.2551420298064355,
            7.732197294565827,
            4.260643127341771,
            6.152444120828289,
            2.227502511042676,
            1.3340698420794308,
            1.6510608368657897,
            4.39552333267791,
            1.6351642634229449,
            6.720559087600739,
            7.567106844111134,
            1.2519921711264994,
            1.4799154955400764,
            1.8700439605296832,
            1.093617514945719,
            9.01275905028318,
            4.061754231796428,
            3.164922807079911,
            3.143045427063748,
            1.5041370739547428,
            3.2724657427650397,
            3.417354138344965,
            6.542082791701956,
            5.041036314516693,
            1.6045703502330106,
            3.3414484076218662,
            2.7331540049373424,
            3.79330061680888,
            8.63779421882033,
            3.1801753567323807,
            4.895499582941207,
            1.2561920524543828,
            3.4964415152459996,
            4.782566213187739,
            4.551034876923523,
            2.73614901845206,
            1.7168281520084105,
            3.9217030138268765,
            7.248845868769475,
            1.8022793537823816,
            5.481676355430747,
            8.653978722908901,
            5.951251247784173,
            2.888960086348835,
            3.750271258530173,
            4.517376304949411,
            6.051232026772094,
            3.5886647912475036,
            2.037234182343205,
            5.271392422028121,
            2.8572620643850284,
            4.773005975750628,
            6.66144901036031,
            3.6212058838221375,
            2.5257762803929107,
            2.5906765222752455,
            2.1090611314869174,
            6.929197380913999,
            1.4098980958771328,
            7.496508082680751,
            5.608986478326592,
            2.175678455800954,
            6.6351554715417445,
            3.805222198144343,
            4.523453262751246,
            1.8416115511475466,
            5.044645652087459,
            1.4989714285040214,
            5.844474409322508,
            4.536862694476934,
            1.3965850769957864,
            4.149045526448823,
            2.7871770590564218,
            3.052621019556693,
            3.870741903589138,
            1.6479621875927128,
            6.834061150286603,
            3.6489191796482485,
            3.4859686803754033,
            5.409698524229855,
            1.7565170501496026,
            4.283453014699891,
            5.128253620366694,
            4.869904892953897,
            2.2324390356156365,
            4.614479831383094,
            5.666982084230242,
            3.340756906512815,
            3.5150237548404157,
            4.993388314090618,
            1.7142049244624797,
            1.9049902623443256,
            4.402662535986304,
            8.771161770936283,
            4.042976328289106,
            3.8639348269469123,
            2.7236767359560314,
            1.7110462015095553,
            1.6382477505296515,
            4.569157052294179,
            2.054799833316604,
            4.0237481102445996,
            6.56031064983421,
            2.691363891269495,
            2.588759457744727,
            2.9779783115754013,
            3.126843169996894,
            3.885405142958953,
            8.053033414022122,
            2.77506564016161,
            3.0381062855040297,
            9.90450274144717,
            4.886805876724099,
            3.330121759629845,
            2.1624814665627743,
            3.531203591624277,
            1.3779531011260544,
            1.9168434526044662,
            3.0949139923816262,
            3.657899611662933,
            2.2883012081393046,
            5.967044820769308,
            3.881705509404706,
            3.5073762276178146,
            3.3449301687443573,
            7.45422464851573,
            5.932213161460171,
            5.21131900484489,
            1.5239650565327352,
            3.536571967919709,
            1.932913987992797,
            3.1286191235700156,
            9.524030341529306,
            5.428105346824928,
            3.2860976763418686,
            6.331638947122046,
            7.701823489325471,
            2.6230074721579166,
            2.149101843913856,
            3.3708235823553254,
            3.2545080383982756,
            3.4500463414113263,
            2.3467076117668233,
            10.675412547771414,
            2.5992125777205843,
            3.3360493969934764,
            3.840959539506068,
            9.850869691653068,
            3.060016235742392,
            2.253805838484651,
            3.712146561928004,
            5.439778162375586,
            10.03866779257671,
            1.94369044273424,
            4.191128810936307,
            2.3592056494636005,
            4.499657272276414,
            12.942247872384844,
            6.581192605543767,
            7.549078732731062,
            3.49403131303363,
            8.340052577722927,
            5.494971722668838,
            13.511708835252188,
            4.693266365762429,
            1.7489836873383677,
            4.9237907492514035,
            4.014679006692068,
            3.32372447820178,
            1.7182659660860025,
            16.998877577184835,
            3.467478966807464,
            4.433296318304915,
            6.181473028863488,
            3.524968879202839,
            2.892857484960585,
            2.8087293669247653,
            5.579701659245514,
            2.4146151787212347,
            4.091392881686932,
            6.685613403077936,
            7.926823659636203,
            4.431672187411348,
            8.485557188901405,
            3.751914981465912,
            2.5251253885007907,
            6.250233950684245,
            2.848781401781672,
            4.347786769133044,
            4.066268433487575,
            9.530055175582689,
            5.614737683967288,
            2.3168528411417855,
        ];

        assert_eq!(daily_hwd.as_ref().unwrap().len(), expected_result.len());

        for (x, y) in daily_hwd.unwrap().iter().zip(expected_result.iter()) {
            assert_relative_eq!(x, y, max_relative = 1e-4)
        }
    }

    #[rstest]
    fn test_edit_hot_water_distribution(mut test_input: InputForProcessing) {
        let tfa = calc_tfa(&test_input);
        edit_hot_water_distribution(&mut test_input, tfa).unwrap();

        let expected_hot_water_distribution_inner: WaterPipeworkSimple =
            serde_json::from_value(json!(
                    {
                        "location": "internal",
                        "external_diameter_mm": 27,
                        "insulation_thermal_conductivity": 0.035,
                        "insulation_thickness_mm": 20,
                        "internal_diameter_mm": 25,
                        "length": 8.0,
                        "pipe_contents": "water",
                        "surface_reflectivity": false
                    }
            ))
            .unwrap();

        let actual_hot_water_distribution_inner =
            test_input.water_distribution().unwrap().first().unwrap();

        assert_eq!(
            *actual_hot_water_distribution_inner,
            expected_hot_water_distribution_inner
        );
    }

    // this test does not exist in Python HEM
    #[rstest]
    fn test_remove_pv_diverter_if_present(mut test_input: InputForProcessing) {
        let diverter = json!({
            "StorageTank": "hw cylinder",
            "HeatSource": "immersion"
        });
        let energy_supply_key = EnergySupplyKey::MainsElectricity;
        test_input.add_diverter_to_energy_supply(energy_supply_key, diverter);

        remove_pv_diverter_if_present(&mut test_input);
        let energy_supply = test_input.energy_supply_by_key(energy_supply_key).unwrap();
        assert!(energy_supply.diverter.is_none())
    }

    // this test does not exist in Python HEM
    #[rstest]
    fn test_remove_electric_battery_if_present(mut test_input: InputForProcessing) {
        let electric_battery = json!({
            "capacity": 5,
            "charge_discharge_efficiency_round_trip": 10,
            "battery_age": 2,
            "minimum_charge_rate_one_way_trip": 42,
            "maximum_charge_rate_one_way_trip": 43,
            "maximum_discharge_rate_one_way_trip": 44,
            "battery_location": "inside"
        });
        let energy_supply_key = EnergySupplyKey::MainsElectricity;
        test_input.add_electric_battery_to_energy_supply(energy_supply_key, electric_battery);

        remove_electric_battery_if_present(&mut test_input);
        let energy_supply = test_input.energy_supply_by_key(energy_supply_key).unwrap();
        assert!(energy_supply.electric_battery.is_none());
    }

    #[rstest]
    fn test_edit_space_cool_system(mut test_input: InputForProcessing) {
        test_input.set_part_o_active_cooling_required(true);
        let _ = edit_space_cool_system(&mut test_input);
        let space_cool_system = test_input.space_cool_system().unwrap();

        for system in space_cool_system.values() {
            assert_eq!(system.efficiency, 5.1);
            assert_eq!(system.frac_convective, 0.95);
            assert_eq!(system.energy_supply, EnergySupplyType::Electricity);
        }
    }

    // this test does not exist in Python HEM
    #[rstest]
    fn test_remove_onsite_generation_if_present(mut test_input: InputForProcessing) {
        test_input
            .set_on_site_generation(json!({
                    "PV 1": {
                        "type": "PhotovoltaicSystem",
                        "peak_power": 2.5,
                        "ventilation_strategy": "moderately_ventilated",
                        "pitch": 30,
                        "orientation360": 180,
                        "base_height":10,
                        "height":1,
                        "width":1,
                        "EnergySupply": "mains elec",
                        "shading": [],
                        "inverter_peak_power": 2,
                        "inverter_is_inside": false
                    }
            }))
            .unwrap();

        remove_onsite_generation_if_present(&mut test_input);
        assert_eq!(test_input.on_site_generation(), None);
    }

    #[rstest]
    fn test_add_solar_pv_house_only(mut test_input: InputForProcessing) {
        let expected_result: OnSiteGeneration = serde_json::from_value(json!({"PV1": {
                "EnergySupply": "mains elec",
                "orientation360": 180,
                "peak_power": 4.444444444444445,
                "inverter_peak_power": 4.444444444444445,
                "inverter_is_inside": false,
                "pitch": 45,
                "type": "PhotovoltaicSystem",
                "ventilation_strategy": "moderately_ventilated",
                "shading": [],
                "base_height": 1,
                "width": 6.324555320336759,
                "height": 3.1622776601683795
                }
        }))
        .unwrap();

        let is_notional_a = true;
        let is_fee = false;
        let total_floor_area = calc_tfa(&test_input);

        add_solar_pv(&mut test_input, is_notional_a, is_fee, total_floor_area).unwrap();

        assert_eq!(*test_input.on_site_generation().unwrap(), expected_result);
    }
}
