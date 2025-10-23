use super::future_homes_standard::{
    ENERGY_SUPPLY_NAME_ELECTRICITY, HW_TEMPERATURE, LIVING_ROOM_SETPOINT_FHS,
    REST_OF_DWELLING_SETPOINT_FHS, SIMTIME_END, SIMTIME_START, SIMTIME_STEP, calc_n_occupants,
    calc_nbeds, calc_tfa, create_cold_water_feed_temps, create_hot_water_use_pattern,
    minimum_air_change_rate, set_temp_internal_static_calcs,
};
use crate::future_homes_standard::fhs_hw_events::STANDARD_BATH_SIZE;
use anyhow::{anyhow, bail};
use hem::compare_floats::min_of_2;
use hem::core::heating_systems::wwhrs::{WWHRSInstantaneousSystemB, Wwhrs};
use hem::core::schedule::{TypedScheduleEvent, expand_events};
use hem::core::units::{
    JOULES_PER_KILOJOULE, JOULES_PER_KILOWATT_HOUR, WATTS_PER_KILOWATT, convert_profile_to_daily,
};
use hem::core::water_heat_demand::cold_water_source::ColdWaterSource;
use hem::core::water_heat_demand::dhw_demand::{
    DomesticHotWaterDemand, DomesticHotWaterDemandData,
};
use hem::core::water_heat_demand::misc::water_demand_to_kwh;
use hem::corpus::{ColdWaterSources, HtcHlpCalculation, calc_htc_hlp};
use hem::input::{
    BuildingElement, ColdWaterSourceType, GroundBuildingElement, GroundBuildingElementJsonValue,
    HeatPumpSourceType, HeatSourceWetDetails, InputForProcessing, JsonAccessResult,
    SpaceHeatSystemHeatSource, UValueEditableBuildingElement,
    UValueEditableBuildingElementJsonValue, WaterPipeContentsType, WaterPipework,
    WaterPipeworkLocation,
};
use hem::simulation_time::SimulationTime;
use hem::statistics::{np_interp, percentile};
use hem::{
    compare_floats::max_of_2,
    core::space_heat_demand::building_element::{HeatFlowDirection, pitch_class},
};
use indexmap::IndexMap;
use parking_lot::Mutex;
use serde_json::{Value, json};
use smartstring::alias::String;
use std::collections::HashMap;
use std::sync::{Arc, LazyLock};
use tracing::instrument;

const NOTIONAL_WWHRS: &str = "Notional_Inst_WWHRS";
const NOTIONAL_HIU: &str = "notionalHIU";
const NOTIONAL_HP: &str = "notional_HP";
const NOTIONAL_BATH_NAME: &str = "medium";
const NOTIONAL_SHOWER_NAME: &str = "mixer";
const NOTIONAL_OTHER_HW_NAME: &str = "other";
const HEATING_PATTERN: &str = "HeatingPattern_Null";

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
    let is_heat_network = check_heatnetwork_present(input)?;

    // Determine cold water source
    let cold_water_source = input.cold_water_source()?;

    if cold_water_source.len() != 1 {
        bail!("The FHS Notional wrapper expects exactly one cold water type to be set.");
    }

    let cold_water_source = {
        let first_cold_water_source = cold_water_source.first();
        let (cold_water_source, _) = first_cold_water_source.as_ref().unwrap();
        **cold_water_source
    };

    // Retrieve the number of bedrooms and total volume
    let bedroom_number = input.number_of_bedrooms()?.ok_or_else(|| {
        anyhow!("FHS Notional wrapper expects number of bedrooms to be provided.")
    })?;
    let total_volume = input.total_zone_volume()?;

    // Determine the TFA
    let total_floor_area = calc_tfa(input)?;

    edit_lighting_efficacy(input)?;
    edit_opaque_adjztu_elements(input)?;
    edit_transparent_element(input)?;
    edit_glazing_for_glazing_limit(input, total_floor_area)?;
    edit_ground_floors(input)?;
    edit_thermal_bridging(input)?;

    // modify bath, shower and other dhw characteristics
    edit_bath_shower_other(input, cold_water_source)?;

    // add WWHRS if needed (and remove any existing systems)
    remove_wwhrs_if_present(input)?;
    add_wwhrs(input, cold_water_source, is_notional_a, is_fee)?;

    // modify hot water distribution
    edit_hot_water_distribution(input, total_floor_area)?;

    // remove on-site generation, pv diverter or electric battery if present
    remove_onsite_generation_if_present(input)?;
    remove_pv_diverter_if_present(input)?;
    remove_electric_battery_if_present(input)?;

    // modify ventilation
    let minimum_air_change_rate =
        minimum_air_change_rate(input, total_floor_area, total_volume, bedroom_number);
    // convert to m3/h
    let minimum_air_flow_rate = minimum_air_change_rate * total_volume;
    edit_infiltration_ventilation(input, is_notional_a, minimum_air_flow_rate)?;

    // edit space heating system
    edit_space_heating_system(
        input,
        cold_water_source,
        total_floor_area,
        is_heat_network,
        is_fee,
    )?;

    // modify air-conditioning
    edit_space_cool_system(input)?;

    // add solar pv
    add_solar_pv(input, is_notional_a, is_fee, total_floor_area)?;

    Ok(())
}

fn check_heatnetwork_present(input: &InputForProcessing) -> anyhow::Result<bool> {
    Ok(input.heat_source_wet()?.values().any(|source| {
        matches!(
            source,
            HeatSourceWetDetails::Hiu { .. }
                | HeatSourceWetDetails::HeatPump {
                    source_type: HeatPumpSourceType::HeatNetwork,
                    ..
                }
        )
    }))
}

/// Apply notional lighting efficacy
/// efficacy = 120 lm/W
fn edit_lighting_efficacy(input: &mut InputForProcessing) -> JsonAccessResult<()> {
    let lighting_efficacy = 120.0;
    input.set_lighting_efficacy_for_all_zones(lighting_efficacy)?;

    Ok(())
}

/// Apply Notional infiltration specifications
/// Notional option A pressure test result at 50Pa = 4 m3/h.m2
/// Notional option B pressure test result at 50Pa = 5 m3/h.m2
/// All passive openings count are set to zero
/// Mechanical extract fans count follows the Actual dwelling,
/// with the exception that there must be at least one per wet room
fn edit_infiltration_ventilation(
    input: &mut InputForProcessing,
    is_notional_a: bool,
    minimum_air_flow_rate: f64,
) -> anyhow::Result<()> {
    // pressure test results dependent on Notional option A or B
    let test_result = if is_notional_a { 4. } else { 5. };

    let number_of_wet_rooms = input.number_of_wet_rooms()?;

    let infiltration_ventilation = input.infiltration_ventilation_mut()?;

    {
        let leaks = infiltration_ventilation
            .entry("Leaks")
            .or_insert(json!({}))
            .as_object_mut()
            .ok_or(anyhow::anyhow!("Leaks was expected to be an object"))?;
        leaks.insert("test_pressure".into(), json!(50.));
        leaks.insert("test_result".into(), json!(test_result));
    }

    // all openings set to 0
    // delete all combustion appliances Cowls and PDUs.
    infiltration_ventilation.insert("PDUs".into(), json!({}));
    infiltration_ventilation.insert("Cowls".into(), json!({}));
    infiltration_ventilation.insert("CombustionAppliances".into(), json!({}));

    if is_notional_a {
        // Notional option A uses continuous extract, so no intermittent extract fans
        // Continuous decentralised mechanical extract ventilation

        infiltration_ventilation.insert(
            "MechanicalVentilation".into(),
            json!({
            "Decentralised_Continuous_MEV_for_notional":{
                "sup_air_flw_ctrl": "ODA",
                "sup_air_temp_ctrl": "CONST",
                "vent_type": "Decentralised continuous MEV",
                "SFP":0.15,
                "EnergySupply": "mains elec",
                "design_outdoor_air_flow_rate": minimum_air_flow_rate
            }}),
        );
    } else {
        // extract_fans follow the same as the actual dwelling
        // but there must be a minimum of one extract fan
        // per wet room, as per ADF guidance
        let wet_rooms_count = number_of_wet_rooms.ok_or_else(|| {
            anyhow!("missing NumberOfWetRooms - required for FHS notional building")
        })?;
        if wet_rooms_count <= 1 {
            bail!("invalid/missing NumberOfWetRooms ({wet_rooms_count})");
        }
        let mut mech_vents: IndexMap<String, Value> = Default::default();
        for i in 0..wet_rooms_count {
            let mech_vent = json!(
                    {
                    "sup_air_flw_ctrl": "ODA",
                    "sup_air_temp_ctrl": "CONST",
                    "vent_type": "Intermittent MEV",
                    "SFP": 0.15,
                    "EnergySupply": "mains elec",
                    "design_outdoor_air_flow_rate": 80
                }
            );
            mech_vents.insert(i.to_string().into(), mech_vent);
        }
        infiltration_ventilation.insert("MechanicalVentilation".into(), json!(mech_vents));
    }

    Ok(())
}

/// Apply notional u-value (W/m2K) to:
///
/// external elements: walls (0.18), doors (1.0), roofs (0.11), exposed floors (0.13)
/// elements adjacent to unheated space: walls (0.18), ceilings (0.11), floors (0.13)
/// to differentiate external doors from walls, user input: is_external_door
fn edit_opaque_adjztu_elements(input: &mut InputForProcessing) -> anyhow::Result<()> {
    let mut opaque_adjztu_building_elements =
        input.all_opaque_and_adjztu_building_elements_mut_u_values()?;

    for mut building_element in opaque_adjztu_building_elements
        .iter_mut()
        .map(|json_map| UValueEditableBuildingElementJsonValue(json_map))
    {
        let pitch_class = pitch_class(building_element.pitch()?);
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
                            bail!(
                                "Opaque building element input needed value for is_external_door field."
                            );
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
        building_element.remove_thermal_resistance_construction();
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

    let mut building_elements = input.all_transparent_building_elements_mut()?;

    for mut building_element in building_elements
        .iter_mut()
        .map(|json_map| UValueEditableBuildingElementJsonValue(json_map))
    {
        let pitch_class = pitch_class(building_element.pitch()?);
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
                building_element.remove_thermal_resistance_construction();
            }
            _ => {
                // if it is not a roof light, it is a glazed door or window
                building_element.set_u_value(1.2);
                building_element.remove_thermal_resistance_construction();
            }
        }
    }

    Ok(())
}

///Split windows/rooflights and walls/roofs into dictionaries.
fn split_glazing_and_walls(
    input: &mut InputForProcessing,
) -> anyhow::Result<(
    IndexMap<String, BuildingElement>,
    IndexMap<String, BuildingElement>,
)> {
    let mut windows_rooflight: IndexMap<String, BuildingElement> = Default::default();
    let mut walls_roofs: IndexMap<String, BuildingElement> = Default::default();

    let building_elements = input.all_building_elements()?;

    for (name, building_element) in building_elements {
        match building_element {
            BuildingElement::Transparent { .. } => {
                windows_rooflight.insert(String::from(name), building_element);
            }
            BuildingElement::Opaque { .. } => {
                walls_roofs.insert(String::from(name), building_element);
            }
            _ => continue,
        }
    }

    Ok((windows_rooflight, walls_roofs))
}

///Calculate difference between old  and new glazing area and adjust the glazing areas
fn calculate_area_diff_and_adjust_glazing_area(
    input: &mut InputForProcessing,
    linear_reduction_factor: f64,
    window_rooflight_element: &BuildingElement,
    building_element_reference: &str,
) -> anyhow::Result<f64> {
    if let BuildingElement::Transparent { height, width, .. } = window_rooflight_element {
        let old_area = height * width;
        let new_height = height * linear_reduction_factor;
        let new_width = width * linear_reduction_factor;

        input.set_numeric_field_for_building_element(
            building_element_reference,
            "height",
            new_height,
        )?;
        input.set_numeric_field_for_building_element(
            building_element_reference,
            "width",
            new_width,
        )?;

        let new_area = new_height * new_width;

        Ok(old_area - new_area)
    } else {
        panic!("This function expects to be called for a Transparent BuildingElement only.")
    }
}

///Find all walls/roofs with same orientation and pitch as this window/rooflight.
fn find_walls_roofs_with_same_orientation_and_pitch(
    wall_roofs: &[&BuildingElement],
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

    for element in input.all_building_elements()?.values() {
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
            let u_value = u_value.ok_or_else(|| {
                anyhow!(
                    "FHS notional wrapper needs transparent building elements to have u values set."
                )
            })?;
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
        .all_building_elements()?
        .values()
        .filter_map(|el| match el {
            BuildingElement::Transparent { .. } => Some(el.height().unwrap() * el.width().unwrap()),
            _ => None,
        })
        .sum();

    let max_glazing_area_fraction = calc_max_glazing_area_fraction(input, total_floor_area)?;
    let max_glazing_area = max_glazing_area_fraction * total_floor_area;

    let (windows_rooflight, walls_roofs) = split_glazing_and_walls(input)?;

    if total_glazing_area > max_glazing_area {
        let linear_reduction_factor = (max_glazing_area / total_glazing_area).sqrt();
        // TODO: deal with case where linear_reduction_factor is NaN (sqrt() is NaN if called on a
        //       negative number, max_glazing_area could come back as a negative number from calc_max_glazing_area_fraction
        //       To do this, we may need to capture a sample input that induces this to happen in the Python, and request
        //       upstream for how to deal with this.

        for (building_element_reference, window_rooflight_element) in windows_rooflight {
            let area_diff = calculate_area_diff_and_adjust_glazing_area(
                input,
                linear_reduction_factor,
                &window_rooflight_element,
                &building_element_reference,
            )?;

            let same_orientation_indices = find_walls_roofs_with_same_orientation_and_pitch(
                &walls_roofs.values().collect::<Vec<_>>(),
                &window_rooflight_element,
            )?;

            let wall_roof_area_total = same_orientation_indices
                .iter()
                .filter_map(|i| match walls_roofs.values().nth(*i).unwrap() {
                    BuildingElement::Opaque { ref area, .. } => Some(*area),
                    _ => None,
                })
                .sum::<f64>();

            for i in same_orientation_indices.iter() {
                let wall_roof = walls_roofs.values().nth(*i).unwrap();
                if let BuildingElement::Opaque { ref area, .. } = wall_roof {
                    let wall_roof_prop = *area / wall_roof_area_total;

                    let new_area = area + area_diff * wall_roof_prop;

                    input.set_numeric_field_for_building_element(
                        &building_element_reference,
                        "area",
                        new_area,
                    )?;
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

    for mut building_element in input
        .all_ground_building_elements_mut()?
        .into_iter()
        .map(GroundBuildingElementJsonValue)
    {
        building_element.set_u_value(0.13);
        building_element.set_thermal_resistance_floor_construction(6.12);
        building_element.set_psi_wall_floor_junc(0.16);
    }

    Ok(())
}

/// The notional building must follow the same thermal bridges as specified in
/// SAP10.2 Table R2
///
/// TODO (from Python) - how to deal with ThermalBridging when lengths are not specified?
pub(crate) fn edit_thermal_bridging(input: &mut InputForProcessing) -> anyhow::Result<()> {
    let mut thermal_bridging_elements = input.all_thermal_bridging_elements()?;

    for element in thermal_bridging_elements
        .iter_mut()
        .flat_map(|group| group.values_mut())
        .filter_map(|bridging| bridging.as_object_mut())
    {
        let bridge_type: String = element
            .get("type")
            .and_then(|bridge_type| bridge_type.as_str())
            .ok_or_else(|| anyhow!("Thermal bridging type was expected to be set."))?
            .into();
        match bridge_type.as_str() {
            "ThermalBridgePoint" => {
                element.insert("heat_transfer_coeff".into(), json!(0.));
            }
            "ThermalBridgeLinear" => {
                let junction_type = element.get("junction_type").and_then(|junc| junc.as_str()).and_then(|junc| if TABLE_R2.contains_key(junc) {Some(junc)} else {None}).ok_or_else(|| anyhow!("Thermal bridging junction type was expected to be set and one of the values in SAP10.2 Table R2."))?;
                element.insert(
                    "linear_thermal_transmittance".into(),
                    json!(TABLE_R2[junction_type]),
                );
            }
            unknown_type => bail!(
                "Thermal bridging type was expected to be set and either ThermalBridgePoint or ThermalBridgeLinear was expected, but {unknown_type} was found."
            ),
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

///  Apply heat network settings to notional building calculation in project_dict.
fn edit_add_heatnetwork_heating(
    input: &mut InputForProcessing,
    cold_water_source: ColdWaterSourceType,
) -> anyhow::Result<()> {
    let heat_network_name = "_notional_heat_network";

    let notional_heat_network = json!(
     {
        NOTIONAL_HIU: {
            "type": "HIU",
            "EnergySupply": heat_network_name,
            "power_max": 45,
            "HIU_daily_loss": 0.8,
            "building_level_distribution_losses": 62,
        }
    });

    let notional_hot_water_source = json!({
        "hw cylinder": {
            "type": "HIU",
            "ColdWaterSource": cold_water_source,
            "HeatSourceWet": NOTIONAL_HIU,
            }
    });

    let heat_network_fuel_data = json!({
        "fuel": "custom",
        "factor":{
            "Emissions Factor kgCO2e/kWh": 0.033,
            "Emissions Factor kgCO2e/kWh including out-of-scope emissions": 0.033,
            "Primary Energy Factor kWh/kWh delivered": 0.75
            }
    });

    input.set_heat_source_wet(notional_heat_network)?;
    input.set_hot_water_source(notional_hot_water_source)?;
    input.add_energy_supply_for_key(heat_network_name, heat_network_fuel_data)?;

    Ok(())
}

fn edit_add_default_space_heating_system(
    input: &mut InputForProcessing,
    design_capacity_overall: f64,
) -> anyhow::Result<()> {
    let factors_35 = IndexMap::from([
        ("A", 1.00),
        ("B", 0.62),
        ("C", 0.55),
        ("D", 0.47),
        ("F", 1.05),
    ]);
    let factors_55 = IndexMap::from([
        ("A", 0.99),
        ("B", 0.60),
        ("C", 0.49),
        ("D", 0.51),
        ("F", 1.03),
    ]);

    let mut capacity_results_dict_35: IndexMap<&str, f64> = Default::default();
    for (record, factor) in factors_35 {
        let value = round_by_precision(factor * design_capacity_overall, 1e3);
        capacity_results_dict_35.insert(record, value);
    }

    let mut capacity_results_dict_55: IndexMap<&str, f64> = Default::default();
    for (record, factor) in factors_55 {
        let value = round_by_precision(factor * design_capacity_overall, 1e3);
        capacity_results_dict_55.insert(record, value);
    }

    let notional_hp = serde_json::from_value(json!(
     {
        NOTIONAL_HP: {
            "EnergySupply": "mains elec",
            "backup_ctrl_type": "TopUp",
            "min_modulation_rate_35": 0.4,
            "min_modulation_rate_55": 0.4,
            "min_temp_diff_flow_return_for_hp_to_operate": 0,
            "modulating_control": true,
            "power_crankcase_heater": 0.01,
            "power_heating_circ_pump": capacity_results_dict_55["F"] * 0.003,
            "power_max_backup": 3,
            "power_off": 0,
            "power_source_circ_pump": 0.01,
            "power_standby": 0.01,
            "sink_type": "Water",
            "source_type": "OutsideAir",
            "temp_lower_operating_limit": -10,
            "temp_return_feed_max": 60,
            "test_data_EN14825": [
                {
                    "capacity": capacity_results_dict_35["A"],
                    "cop": 2.79,
                    "degradation_coeff": 0.9,
                    "design_flow_temp": 35,
                    "temp_outlet": 34,
                    "temp_source": -7,
                    "temp_test": -7,
                    "test_letter": "A"
                },
                {
                    "capacity": capacity_results_dict_35["B"],
                    "cop": 4.29,
                    "degradation_coeff": 0.9,
                    "design_flow_temp": 35,
                    "temp_outlet": 30,
                    "temp_source": 2,
                    "temp_test": 2,
                    "test_letter": "B"
                },
                {
                    "capacity": capacity_results_dict_35["C"],
                    "cop": 5.91,
                    "degradation_coeff": 0.9,
                    "design_flow_temp": 35,
                    "temp_outlet": 27,
                    "temp_source": 7,
                    "temp_test": 7,
                    "test_letter": "C"
                },
                {
                    "capacity": capacity_results_dict_35["D"],
                    "cop": 8.02,
                    "degradation_coeff": 0.9,
                    "design_flow_temp": 35,
                    "temp_outlet": 24,
                    "temp_source": 12,
                    "temp_test": 12,
                    "test_letter": "D"
                },
                {
                    "capacity": capacity_results_dict_35["F"],
                    "cop": 2.49,
                    "degradation_coeff": 0.9,
                    "design_flow_temp": 35,
                    "temp_outlet": 35,
                    "temp_source": -10,
                    "temp_test": -10,
                    "test_letter": "F"
                },
                {
                    "capacity": capacity_results_dict_55["A"],
                    "cop": 2.03,
                    "degradation_coeff": 0.9,
                    "design_flow_temp": 55,
                    "temp_outlet": 52,
                    "temp_source": -7,
                    "temp_test": -7,
                    "test_letter": "A"
                },
                {
                    "capacity": capacity_results_dict_55["B"],
                    "cop": 3.12,
                    "degradation_coeff": 0.9,
                    "design_flow_temp": 55,
                    "temp_outlet": 42,
                    "temp_source": 2,
                    "temp_test": 2,
                    "test_letter": "B"
                },
                {
                    "capacity": capacity_results_dict_55["C"],
                    "cop": 4.41,
                    "degradation_coeff": 0.9,
                    "design_flow_temp": 55,
                    "temp_outlet": 36,
                    "temp_source": 7,
                    "temp_test": 7,
                    "test_letter": "C"
                },
                {
                    "capacity": capacity_results_dict_55["D"],
                    "cop": 6.30,
                    "degradation_coeff": 0.9,
                    "design_flow_temp": 55,
                    "temp_outlet": 30,
                    "temp_source": 12,
                    "temp_test": 12,
                    "test_letter": "D"
                },
                {
                    "capacity": capacity_results_dict_55["F"],
                    "cop": 1.87,
                    "degradation_coeff": 0.9,
                    "design_flow_temp": 55,
                    "temp_outlet": 55,
                    "temp_source": -10,
                    "temp_test": -10,
                    "test_letter": "F"
                }
            ],
            "time_constant_onoff_operation": 120,
            "time_delay_backup": 1,
            "type": "HeatPump",
            "var_flow_temp_ctrl_during_test": true
        }
    }))?;

    input.set_heat_source_wet(notional_hp)?;
    Ok(())
}

/// Apply distribution system details to notional building calculation
fn edit_default_space_heating_distribution_system(
    input: &mut InputForProcessing,
    design_capacity: &IndexMap<String, f64>,
) -> anyhow::Result<()> {
    let setpoint_for_sizing = max_of_2(LIVING_ROOM_SETPOINT_FHS, REST_OF_DWELLING_SETPOINT_FHS);

    let design_flow_temp = 45.;
    let design_flow_rate = 12.; // TODO (from Python) what value should this be?
    let n: f64 = 1.34;
    let c_per_rad = 1.89 / (50_f64).powf(n);
    let power_output_per_rad = c_per_rad * (design_flow_temp - setpoint_for_sizing).powf(n);

    // thermal mass specified in kJ/K but required in kWh/K
    let thermal_mass_per_rad = 51.8 * JOULES_PER_KILOJOULE as f64 / JOULES_PER_KILOWATT_HOUR as f64;

    // Initialise space heating system in project dict
    input.remove_space_heat_systems()?;

    for zone_name in input.zone_keys()? {
        let system_name = format!("{zone_name}_SpaceHeatSystem_Notional");
        input.set_space_heat_system_for_zone(&zone_name, &system_name)?;
        let heatsourcewet_name = input
            .heat_source_wet()?
            .first()
            .ok_or_else(|| {
                anyhow!("FHS Notional wrapper expected at least one heat source wet to be defined.")
            })?
            .0
            .clone();

        // Calculate number of radiators
        let emitter_cap = design_capacity.get(&zone_name).ok_or_else(|| {
            anyhow!("FHS Notional wrapper expected a design capacity with name: {zone_name}.")
        })?;
        let number_of_rads = (emitter_cap / power_output_per_rad).ceil();

        // Calculate c and thermal mass
        let c = number_of_rads * c_per_rad;
        let thermal_mass = number_of_rads * thermal_mass_per_rad;

        let space_heat_system_value = json!({
            "type": "WetDistribution",
            "advanced_start": 1,
            "thermal_mass": thermal_mass,
            "emitters": [{"wet_emitter_type": "radiator",
                          "frac_convective": 0.7,
                          "c": c,
                          "n": n}],
            "temp_diff_emit_dsgn": 5,
            "HeatSource": {
                "name": heatsourcewet_name,
                "temp_flow_limit_upper": 65.0
            },
            "ecodesign_controller": {
                    "ecodesign_control_class": 2,
                    "max_outdoor_temp": 20,
                    "min_flow_temp": 21,
                    "min_outdoor_temp": 0
                    },
            "Control": HEATING_PATTERN,
            "design_flow_temp": design_flow_temp as i32,
            "design_flow_rate": design_flow_rate,
            "Zone": zone_name,
            "temp_setback" : 18
        });

        // Create radiator dict for zone
        input.set_space_heat_system_for_key(&system_name, space_heat_system_value)?;
    }

    Ok(())
}

/// Edit distribution system details to notional building heat network
fn edit_heatnetwork_space_heating_distribution_system(
    input: &mut InputForProcessing,
) -> anyhow::Result<()> {
    let space_heat_system_keys: Vec<String> = input.space_heat_system_keys()?;

    for system in space_heat_system_keys {
        input.set_advance_start_for_space_heat_system(&system, 1.)?;
    }

    input.set_temperature_setback_for_space_heat_systems(None)?;

    let notional_heat_source: SpaceHeatSystemHeatSource =
        serde_json::from_value(json!({"name": NOTIONAL_HIU}))?;
    input.set_heat_source_for_all_space_heat_systems(notional_heat_source)?;

    Ok(())
}
fn edit_bath_shower_other(
    input: &mut InputForProcessing,
    cold_water_source_type: ColdWaterSourceType,
) -> anyhow::Result<()> {
    // Define Bath, Shower, and Other DHW outlet
    let notional_bath = json!({ NOTIONAL_BATH_NAME: {
            "ColdWaterSource": cold_water_source_type,
            "flowrate": 12,
            "size": STANDARD_BATH_SIZE
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

fn remove_wwhrs_if_present(input: &mut InputForProcessing) -> anyhow::Result<()> {
    if input.wwhrs()?.is_some() {
        input.remove_wwhrs()?;
    }

    Ok(())
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

    let storeys_in_building = match input.build_type()?.as_str() {
        "house" => input.storeys_in_building()?,
        "flat" => 1,
        unknown_type => bail!("Encountered unexpected building type '{unknown_type}'"),
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
        .cold_water_source()?
        .iter()
        .map(|(key, source)| {
            (
                *key,
                Arc::from(ColdWaterSource::new(
                    source.temperatures.clone(),
                    source.start_day,
                    source.time_series_step,
                )),
            )
        })
        .collect();

    let wwhrs: IndexMap<String, Arc<Mutex<Wwhrs>>> = if let Some(waste_water_heat_recovery) =
        input.wwhrs()?
    {
        let notional_wwhrs = waste_water_heat_recovery.get(NOTIONAL_WWHRS).ok_or_else(|| anyhow!("A {} entry for WWHRS was expected to have been set in the FHS Notional wrapper.", NOTIONAL_WWHRS))?;
        [(
            String::from(NOTIONAL_WWHRS),
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
        ("Shower", NOTIONAL_SHOWER_NAME),
        ("Bath", NOTIONAL_BATH_NAME),
        ("Other", NOTIONAL_OTHER_HW_NAME),
    ];

    // Initialize a single schedule dictionary
    let mut event_schedules: Vec<Option<Vec<TypedScheduleEvent>>> = vec![None; total_timesteps];

    // Populate the event_schedules dictionary using the modified expand_events function
    for (event_type, event_name) in event_types_names_list {
        let event_data = input.water_heating_event_by_type_and_name(event_type, event_name)?.ok_or_else(|| anyhow!("FHS Notional wrapper expected water heating events with type '{event_type}' and name '{event_name}'"))?.iter().map(Into::into).collect::<Vec<_>>();
        event_schedules = expand_events(
            event_data,
            sim_timestep,
            total_timesteps,
            event_name,
            serde_json::from_value(json!(event_type))?,
            event_schedules,
        )?;
    }

    let dhw_demand = DomesticHotWaterDemand::new(
        input
            .showers()?
            .map(|showers| serde_json::from_value(json!(showers)))
            .transpose()?
            .unwrap_or_default(),
        input
            .baths()?
            .map(|baths| serde_json::from_value(json!(baths)))
            .transpose()?
            .unwrap_or_default(),
        input
            .other_water_uses()?
            .map(|other_water_uses| serde_json::from_value(json!(other_water_uses)))
            .transpose()?
            .unwrap_or_default(),
        input.water_distribution()?.clone(),
        &cold_water_sources,
        &wwhrs,
        &Default::default(),
        event_schedules,
    )?;

    // For each timestep, calculate HW draw
    let total_steps = simtime.total_steps();
    let mut hw_energy_demand = vec![0.0; total_steps];
    for (t_idx, t_it) in simtime.iter().enumerate() {
        let DomesticHotWaterDemandData { hw_demand_vol, .. } =
            dhw_demand.hot_water_demand(t_it, HW_TEMPERATURE)?;

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
    let cylinder_vol = match input.hot_water_cylinder_volume()? {
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

    let length_max = match input.build_type()?.as_str() {
        "flat" => 0.05 * total_floor_area,
        "house" => {
            0.05 * input.ground_floor_area()?.ok_or_else(|| {
                anyhow!("FHS Notional wrapper expected ground floor area to be set for a house.")
            })?
        }
        unknown_type => bail!("Encountered unexpected building type '{unknown_type}'"),
    };

    let mut primary_pipework = input.primary_pipework_clone()?;

    match primary_pipework {
        None => {
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
        }
        Some(ref mut primary_pipework) => {
            for pipework in primary_pipework.iter_mut() {
                let length = pipework.length;
                let internal_diameter_mm =
                    pipework.internal_diameter_mm.max(internal_diameter_mm_min);
                let external_diameter_mm =
                    pipework.external_diameter_mm.max(external_diameter_mm_min);

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
    }

    Ok(primary_pipework.unwrap())
}

fn edit_hot_water_distribution(
    input: &mut InputForProcessing,
    total_floor_area: f64,
) -> anyhow::Result<()> {
    // hot water dictionary
    let mut hot_water_distribution_inner_list = vec![];

    for item in input.water_distribution()?.into_iter().flatten() {
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

    let length_max = match input.build_type()?.as_str() {
        "flat" => 0.2 * total_floor_area,
        "house" => {
            0.2 * input
                .ground_floor_area()?
                .ok_or_else(|| anyhow!("Notional wrapper expected a ground floor area to be set"))?
        }
        unknown_type => bail!("Encountered unexpected building type '{unknown_type}'"),
    };

    // Iterate over hot_water_distribution_inner_list
    for hot_water_distribution_inner in hot_water_distribution_inner_list {
        // hot water distribution (inner) length should not be greater than maximum length

        let length_actual = hot_water_distribution_inner.length;
        let length = min_of_2(length_actual, length_max);

        // Update internal diameter to minimum if not present and should not be lower than the minimum
        let internal_diameter_mm = hot_water_distribution_inner.internal_diameter_mm;
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

fn remove_pv_diverter_if_present(
    input: &mut InputForProcessing,
) -> JsonAccessResult<&mut InputForProcessing> {
    input.remove_all_diverters_from_energy_supplies()
}

fn remove_electric_battery_if_present(
    input: &mut InputForProcessing,
) -> JsonAccessResult<&mut InputForProcessing> {
    input.remove_all_batteries_from_energy_supplies()
}

fn edit_space_heating_system(
    input: &mut InputForProcessing,
    cold_water_source: ColdWaterSourceType,
    total_floor_area: f64,
    is_heat_network: bool,
    is_fee: bool,
) -> anyhow::Result<()> {
    // FEE calculation which doesn't need the space heating system at this stage.
    if !is_fee {
        // If Actual dwelling is heated with heat networks - Notional heated with HIU.
        // Otherwise, notional heated with an air to water heat pump
        if is_heat_network {
            edit_add_heatnetwork_heating(input, cold_water_source)?;
            edit_heatnetwork_space_heating_distribution_system(input)?;
        } else {
            let (design_capacity_map, design_capacity_overall) = calc_design_capacity(input)?;
            edit_add_default_space_heating_system(input, design_capacity_overall)?;
            edit_default_space_heating_distribution_system(input, &design_capacity_map)?;
            edit_storagetank(input, cold_water_source, total_floor_area)?;
        }
    }

    Ok(())
}

fn edit_space_cool_system(input: &mut InputForProcessing) -> anyhow::Result<()> {
    let part_o_active_cooling_required = input.part_o_active_cooling_required()?.unwrap_or(false);

    if part_o_active_cooling_required {
        input.set_efficiency_for_all_space_cool_systems(5.1)?;
        input.set_frac_convective_for_all_space_cool_systems(0.95)?;
        input.set_energy_supply_for_all_space_cool_systems(ENERGY_SUPPLY_NAME_ELECTRICITY)?;
    }

    Ok(())
}

#[instrument(skip(input))]
fn calc_design_capacity(
    input: &InputForProcessing,
) -> anyhow::Result<(IndexMap<String, f64>, f64)> {
    // Create a deep copy as init_resistance_or_uvalue() will add u_value & r_c
    // which will raise warning when called second time
    let mut clone = input.clone();

    // Calculate heat transfer coefficients and heat loss parameters
    set_temp_internal_static_calcs(&mut clone)?;
    let HtcHlpCalculation {
        htc_map: htc_dict, ..
    } = calc_htc_hlp(&clone.as_input_for_calc_htc_hlp()?)?;

    // Calculate design capacity
    let min_air_temp = *input.external_conditions()?.air_temperatures.as_ref().ok_or_else(|| anyhow!("FHS Notional wrapper expected to have air temperatures merged onto the input structure."))?.iter().min_by(|a, b| a.total_cmp(b)).ok_or_else(|| anyhow!("FHS Notional wrapper expects air temperature list set on input structure not to be empty."))?;
    let set_point = LIVING_ROOM_SETPOINT_FHS.max(REST_OF_DWELLING_SETPOINT_FHS);
    let temperature_difference = set_point - min_air_temp;
    let design_capacity_map: IndexMap<String, f64> = input
        .zone_keys()?
        .into_iter()
        .map(|key| {
            (key.to_owned(), {
                let design_heat_loss = htc_dict[&key] * temperature_difference;
                let design_capacity = 2. * design_heat_loss;
                design_capacity / WATTS_PER_KILOWATT as f64
            })
        })
        .collect();

    let design_capacity_overall = design_capacity_map.values().sum::<f64>();

    Ok((design_capacity_map, design_capacity_overall))
}

/// Initialise temperature setpoints for all zones.
/// The initial set point is needed to call the Project class.
/// Set as 18C for now. The FHS wrapper will overwrite temp_setpnt_init '''
#[cfg(test)]
fn initialise_temperature_setpoints(input: &mut InputForProcessing) -> anyhow::Result<()> {
    for zone_key in input.zone_keys()? {
        input.set_init_temp_setpoint_for_zone(zone_key.as_str(), 18.)?;
    }
    Ok(())
}

fn remove_onsite_generation_if_present(input: &mut InputForProcessing) -> JsonAccessResult<()> {
    input.remove_on_site_generation()?;
    Ok(())
}

fn add_solar_pv(
    input: &mut InputForProcessing,
    is_notional_a: bool,
    is_fee: bool,
    total_floor_area: f64,
) -> anyhow::Result<()> {
    let number_of_storeys = input.storeys_in_building()?;

    // PV is included in the notional if the building contains 15 stories or
    // less that contain dwellings.
    if number_of_storeys <= 15 && is_notional_a && !is_fee {
        let ground_floor_area = input
            .ground_floor_area()?
            .ok_or_else(|| anyhow!("Notional wrapped expected ground floor area to be set"))?;
        let (peak_kw, base_height_pv) = match input.build_type()?.as_str() {
            "house" => {
                let peak_kw = ground_floor_area * 0.4 / 4.5;
                let base_height_pv = input.max_base_height_from_building_elements()?.ok_or_else(|| anyhow!("Notional wrapper expected at least one building element with a base height"))?;

                (peak_kw, base_height_pv)
            }
            "flat" => {
                let peak_kw = total_floor_area * 0.4 / (4.5 * number_of_storeys as f64);
                let zone_total_volume = input.total_zone_volume()?;
                let zone_total_area = input.total_zone_area()?;
                let base_height_pv =
                    (zone_total_volume / zone_total_area + 0.3) * number_of_storeys as f64;

                (peak_kw, base_height_pv)
            }
            unknown_type => bail!("Unexpected building type '{unknown_type}' encountered"),
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
                "inverter_peak_power_ac": peak_kw,
                "inverter_peak_power_dc": peak_kw,
                "inverter_is_inside": false,
                "inverter_type": "optimised_inverter",
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

fn round_by_precision(src: f64, precision: f64) -> f64 {
    (precision * src).round() / precision
}

#[cfg(test)]
mod tests {
    use crate::core::space_heat_demand::building_element::{HeatFlowDirection, pitch_class};

    use super::*;
    use crate::input::{
        self, EnergySupplyDetails, HeatSourceWet, HeatSourceWetDetails, WaterPipeworkSimple,
    };
    use crate::input::{HotWaterSource, WasteWaterHeatRecovery};
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
        InputForProcessing::init_with_json_skip_validation(reader).expect(
            "expected valid test_future_homes_standard_notional_input_data.json to be present",
        )
    }

    #[ignore = "currently failing as calc_design_capacity is failing (our test data is missing some expected fields on external conditions)"]
    #[rstest]
    // this test does not exist in Python HEM
    fn test_apply_fhs_notional_preprocessing(mut test_input: InputForProcessing) {
        let fhs_notional_a_assumptions = true;
        let _fhs_notional_b_assumptions = false;
        let fhs_fee_notional_a_assumptions = false;
        let fhs_fee_notional_b_assumptions = false;

        let actual = apply_fhs_notional_preprocessing(
            &mut test_input,
            fhs_notional_a_assumptions,
            _fhs_notional_b_assumptions,
            fhs_fee_notional_a_assumptions,
            fhs_fee_notional_b_assumptions,
        );
        assert!(actual.is_ok())
    }

    #[rstest]
    // this test does not exist in Python HEM
    fn test_check_heatnetwork_present(test_input: InputForProcessing) {
        assert!(!check_heatnetwork_present(&test_input).unwrap());
    }

    #[rstest]
    fn test_edit_lighting_efficacy(mut test_input: InputForProcessing) {
        let test_input = test_input.borrow_mut();
        edit_lighting_efficacy(test_input).unwrap();

        for zone in test_input.zone_keys().unwrap() {
            let lighting_efficacy = test_input.lighting_efficacy_for_zone(&zone).unwrap();
            assert_eq!(
                lighting_efficacy.expect("expected lighting in zone and efficacy in lighting"),
                120.
            )
        }
    }

    #[rstest]
    fn test_edit_opaque_ajdztu_elements(mut test_input: InputForProcessing) {
        edit_opaque_adjztu_elements(&mut test_input).unwrap();

        // not using the building_element_by_key method here to closly match the Python test

        for building_element in test_input.all_building_elements().unwrap().values() {
            if let input::BuildingElement::Opaque { .. }
            | input::BuildingElement::AdjacentUnconditionedSpace { .. } = building_element
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

        let zone_1_window_0_element = test_input
            .building_element_by_key("zone 1", "window 0")
            .unwrap();

        assert_eq!(
            zone_1_window_0_element
                .get("u_value")
                .and_then(|v| v.as_f64()),
            Some(1.2)
        );
        assert!(
            zone_1_window_0_element
                .get("thermal_resistance_construction")
                .is_none()
        );

        let zone_2_window_0_element = test_input
            .building_element_by_key("zone 2", "window 0")
            .unwrap();

        assert_eq!(
            zone_2_window_0_element
                .get("u_value")
                .and_then(|v| v.as_f64()),
            Some(1.2)
        );
        assert!(
            zone_2_window_0_element
                .get("thermal_resistance_construction")
                .is_none()
        );
    }

    #[rstest]
    fn test_edit_ground_floors(mut test_input: InputForProcessing) {
        let test_input = test_input.borrow_mut();

        edit_ground_floors(test_input).unwrap();

        let zone_1_ground_element = test_input
            .building_element_by_key("zone 1", "ground")
            .unwrap();

        assert_eq!(
            zone_1_ground_element
                .get("u_value")
                .and_then(|v| v.as_f64()),
            Some(0.13)
        );
        assert_eq!(
            zone_1_ground_element
                .get("thermal_resistance_floor_construction")
                .and_then(|v| v.as_f64()),
            Some(6.12)
        );
        assert_eq!(
            zone_1_ground_element
                .get("psi_wall_floor_junc")
                .and_then(|v| v.as_f64()),
            Some(0.16)
        );

        let zone_2_ground_element = test_input
            .building_element_by_key("zone 2", "ground")
            .unwrap();

        assert_eq!(
            zone_2_ground_element
                .get("u_value")
                .and_then(|v| v.as_f64()),
            Some(0.13)
        );
        assert_eq!(
            zone_2_ground_element
                .get("thermal_resistance_floor_construction")
                .and_then(|v| v.as_f64()),
            Some(6.12)
        );
        assert_eq!(
            zone_2_ground_element
                .get("psi_wall_floor_junc")
                .and_then(|v| v.as_f64()),
            Some(0.16)
        );
    }

    #[rstest]
    fn test_edit_thermal_bridgings(mut test_input: InputForProcessing) {
        let test_input = test_input.borrow_mut();

        edit_thermal_bridging(test_input).unwrap();

        for thermal_bridging in test_input
            .all_thermal_bridgings()
            .unwrap()
            .into_iter()
            .filter_map(|t| t.as_object())
            .flat_map(|v| v.values())
            .filter_map(|v| v.as_object())
        {
            let bridging_type = thermal_bridging.get("type").unwrap().as_str().unwrap();
            if bridging_type == "ThermalBridgePoint" {
                assert_eq!(
                    thermal_bridging
                        .get("heat_transfer_coeff")
                        .unwrap()
                        .as_f64()
                        .unwrap(),
                    0.0
                );
            } else if bridging_type == "ThermalBridgeLinear" {
                let junction_type = thermal_bridging
                    .get("junction_type")
                    .unwrap()
                    .as_str()
                    .unwrap();
                assert!(TABLE_R2.contains_key(junction_type));
                assert_eq!(
                    thermal_bridging
                        .get("linear_thermal_transmittance")
                        .unwrap()
                        .as_f64()
                        .unwrap(),
                    TABLE_R2[junction_type]
                );
            } else {
                panic!("Unknown thermal bridging type '{bridging_type}' encountered");
            }
        }
    }

    #[rstest]
    fn test_calc_max_glazing_area_fraction(mut test_input: InputForProcessing) {
        test_input
            .set_zone(zone_input_for_max_glazing_area_test(1.5, None))
            .unwrap();
        assert_eq!(
            calc_max_glazing_area_fraction(&test_input, 80.0).unwrap(),
            0.24375,
            "incorrect max glazing area fraction"
        );
        test_input
            .set_zone(zone_input_for_max_glazing_area_test(1.0, None))
            .unwrap();
        assert_eq!(
            calc_max_glazing_area_fraction(&test_input, 80.0).unwrap(),
            0.25,
            "incorrect max glazing area fraction"
        );
        test_input
            .set_zone(zone_input_for_max_glazing_area_test(1.5, Some(90.)))
            .unwrap();
        assert_eq!(
            calc_max_glazing_area_fraction(&test_input, 80.0).unwrap(),
            0.25,
            "incorrect max glazing area fraction when there are no rooflights"
        );
    }

    // this test does not exist in Python HEM
    #[rstest]
    fn test_calculate_area_diff_and_adjust_glazing_area(mut test_input: InputForProcessing) {
        let linear_reduction_factor: f64 = 0.7001400420140049;

        let window: BuildingElement = serde_json::from_value(json!(
            test_input
                .building_element_by_key("zone 1", "window 0")
                .unwrap()
        ))
        .unwrap();

        let area_diff = calculate_area_diff_and_adjust_glazing_area(
            &mut test_input,
            linear_reduction_factor,
            &window,
            "window 0",
        )
        .unwrap();

        assert_relative_eq!(area_diff, 2.549019607843137);
    }

    fn zone_input_for_max_glazing_area_test(u_value: f64, pitch_override: Option<f64>) -> Value {
        json!({"test_zone": {
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
                    "free_area_height": 0.0,
                    "base_height": 0.0,
                    "max_window_open_area": 3,
                    "mid_height": 1.5,
                    "shading": [],
                    "window_part_list": []
                }
            },
            "ThermalBridging": 1.0,
            "area": 10.0,
            "volume": 20.0
        }})
    }

    // this test does not exist in Python HEM
    #[rstest]
    fn test_edit_add_heatnetwork_heating(mut test_input: InputForProcessing) {
        let heat_network_name = "_notional_heat_network";

        let expected_heat_source_wet: HeatSourceWet = serde_json::from_value(json!({
            NOTIONAL_HIU: {
                "type": "HIU",
                "EnergySupply": heat_network_name,
                "power_max": 45,
                "HIU_daily_loss": 0.8,
                "building_level_distribution_losses": 62,
            }
        }))
        .unwrap();

        let expected_hot_water_source: HotWaterSource = serde_json::from_value(json!({
        "hw cylinder": {
            "type": "HIU",
            "ColdWaterSource": "mains water",
            "HeatSourceWet": NOTIONAL_HIU,
            }
        }))
        .unwrap();

        let heat_network_name = "_notional_heat_network";
        let expected_heat_network_fuel_data: EnergySupplyDetails = serde_json::from_value(json!({
            "fuel": "custom",
            "factor":{
                "Emissions Factor kgCO2e/kWh": 0.033,
                "Emissions Factor kgCO2e/kWh including out-of-scope emissions": 0.033,
                "Primary Energy Factor kWh/kWh delivered": 0.75
                }
        }))
        .unwrap();

        edit_add_heatnetwork_heating(&mut test_input, ColdWaterSourceType::MainsWater).unwrap();

        assert_eq!(
            test_input.heat_source_wet().unwrap(),
            expected_heat_source_wet
        );

        assert_eq!(
            json!(test_input.hot_water_source().unwrap()),
            json!(expected_hot_water_source)
        );

        assert_eq!(
            json!(test_input.energy_supply_by_key(heat_network_name).unwrap()),
            json!(expected_heat_network_fuel_data)
        );
    }

    // this test does not exist in Python HEM
    #[rstest]
    fn test_edit_heatnetwork_space_heating_distribution_system(mut test_input: InputForProcessing) {
        edit_heatnetwork_space_heating_distribution_system(&mut test_input).unwrap();

        for system in test_input.space_heat_system_keys().unwrap() {
            let advanced_start = test_input
                .advanced_start_for_space_heat_system(&system)
                .unwrap();
            assert_eq!(advanced_start, Some(1.));

            let temp_setback = test_input
                .temperature_setback_for_space_heat_system(&system)
                .unwrap();
            assert_eq!(temp_setback, None);

            let heat_source = test_input
                .heat_source_for_space_heat_system(&system)
                .unwrap()
                .unwrap();

            let expected_heat_source = json!({"name": NOTIONAL_HIU});

            assert_eq!(*heat_source, expected_heat_source)
        }
    }

    #[rstest]
    fn test_edit_bath_shower_other(mut test_input: InputForProcessing) {
        // this is the only cold water source type in the test input JSON file
        let cold_water_source_type = ColdWaterSourceType::MainsWater;
        let cold_water_source_type_string = "mains water";

        edit_bath_shower_other(&mut test_input, cold_water_source_type).unwrap();

        let expected_baths = json!({ "medium": {
            "ColdWaterSource": cold_water_source_type_string,
            "flowrate": 12,
            "size": 180.
        }})
        .as_object()
        .unwrap()
        .clone();

        let expected_showers = json!({"mixer": {
            "ColdWaterSource": cold_water_source_type_string,
            "flowrate": 8,
            "type": "MixerShower"
        }})
        .as_object()
        .unwrap()
        .clone();

        let expected_other = json!({"other": {
            "ColdWaterSource": cold_water_source_type_string,
            "flowrate": 6,
        }})
        .as_object()
        .unwrap()
        .clone();

        assert_eq!(test_input.baths().unwrap().unwrap().clone(), expected_baths);
        assert_eq!(
            test_input.showers().unwrap().unwrap().clone(),
            expected_showers
        );
        assert_eq!(
            test_input.other_water_uses().unwrap().unwrap().clone(),
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
        remove_wwhrs_if_present(&mut test_input).unwrap();
        assert!(test_input.wwhrs().unwrap().is_none());
    }

    #[rstest]
    fn test_add_wwhrs(mut test_input: InputForProcessing) {
        test_input.set_storeys_in_building(2).unwrap();

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

        assert_eq!(test_input.wwhrs().unwrap().unwrap().clone(), expected_wwhrs);

        let mixer_shower = test_input
            .showers()
            .unwrap()
            .unwrap()
            .get("mixer")
            .and_then(|mixer| mixer.as_object())
            .unwrap();

        assert_eq!(
            mixer_shower.get("WWHRS").and_then(|wwhrs| wwhrs.as_str()),
            Some("Notional_Inst_WWHRS")
        );
    }

    #[rstest]
    fn test_add_no_wwhrs_for_one_storey_buildings(mut test_input: InputForProcessing) {
        let cold_water_source_type = ColdWaterSourceType::MainsWater;

        add_wwhrs(&mut test_input, cold_water_source_type, true, false).unwrap();

        assert!(test_input.wwhrs().unwrap().is_none());

        let mixer_shower = test_input
            .showers()
            .unwrap()
            .unwrap()
            .get("mixer")
            .and_then(|mixer| mixer.as_object())
            .unwrap();

        assert!(!mixer_shower.contains_key("WWHRS"));
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
        let total_floor_area = calc_tfa(&test_input).unwrap();

        edit_storagetank(&mut test_input, cold_water_source_type, total_floor_area).unwrap();

        let expected_primary_pipework = json!([{
            "location": "internal",
            "internal_diameter_mm": 26.,
                "external_diameter_mm": 28.,
                "length": 2.5,
                "insulation_thermal_conductivity": 0.035,
                "insulation_thickness_mm": 35.,
                "surface_reflectivity": false,
                "pipe_contents": "water"
        }]);

        let expected_hotwater_source = json!({
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
        })
        .as_object()
        .unwrap()
        .clone();

        assert_eq!(
            test_input.hot_water_source().unwrap().clone(),
            expected_hotwater_source
        );
    }

    #[ignore = "Below test used to fail due to randomisation issue, now numbers are even further apart"]
    #[rstest]
    fn test_calc_daily_hw_demand(mut test_input: InputForProcessing) {
        let cold_water_source_type = ColdWaterSourceType::MainsWater;
        let total_floor_area = calc_tfa(&test_input).unwrap();

        // Add notional objects that affect HW demand calc
        edit_bath_shower_other(&mut test_input, cold_water_source_type).unwrap();

        let daily_hwd =
            calc_daily_hw_demand(&mut test_input, total_floor_area, cold_water_source_type);

        let expected_result = [
            4.494866624219392,
            3.6262929661140406,
            2.4792292219363525,
            7.396905949799487,
            2.334290833000372,
            8.938831222369114,
            4.218109245384848,
            3.2095837032274623,
            1.562913391249543,
            8.846662366836481,
            2.573896298797947,
            3.8158857819823955,
            1.8643761342051466,
            1.456499804780102,
            7.426850861522029,
            2.9833503486722512,
            4.217424343066319,
            8.086072907696455,
            4.14341306475027,
            5.363210797769194,
            4.51254160486555,
            4.535867190456099,
            2.7857687977141605,
            1.7560127175725972,
            10.333211623720878,
            2.0533256568949536,
            10.5846961515653,
            2.9116693757992294,
            6.246398935042146,
            1.6696701053184573,
            13.851788339322859,
            3.492087111231953,
            7.271874351886643,
            3.4529488454587005,
            12.374228697052182,
            4.392672154556575,
            1.6028771659496917,
            2.5058963927074895,
            2.075293843606148,
            2.949279475580221,
            8.392209203216268,
            7.314951072027724,
            3.8238937613049946,
            6.812712493984371,
            4.537127728764957,
            6.858233389853893,
            3.994128161632102,
            6.612136728233484,
            10.073004625703325,
            7.1389148991972755,
            1.5377879632668527,
            2.5192092256423533,
            3.5974699592273436,
            2.677722497971631,
            6.641600203287786,
            2.108183420048377,
            2.0324151238352606,
            4.5025717651118,
            1.6287189927962715,
            5.184027724364109,
            2.19733248287286,
            5.538684171666794,
            1.6808281306652284,
            5.413255340871867,
            2.5025554646129446,
            6.9674570352988034,
            4.018149967311069,
            3.598667197534451,
            2.2197290514730836,
            6.818451857176455,
            5.796189225222955,
            7.768257372957939,
            1.9635622695280748,
            4.990639078067053,
            11.326110200617702,
            7.793122367331328,
            4.50364508936643,
            8.379833734745256,
            3.308002750963755,
            5.125036944678628,
            5.130855988470622,
            6.976241946425853,
            9.280525199389762,
            6.879123493336726,
            3.7542978536142275,
            1.0782890932651108,
            6.663709522738787,
            5.7746999922120015,
            3.9743519683699224,
            9.8172995461514,
            2.9545593496596627,
            5.318321839987381,
            3.3213819919472374,
            7.238785487773112,
            0.975438526348773,
            6.4976602476390495,
            3.7068430878370746,
            3.6626428865004845,
            3.448702673630278,
            2.836638476910602,
            4.504302459687092,
            4.594078645382666,
            1.280400852785038,
            9.635660153417774,
            2.3614201923456397,
            5.013620232473153,
            5.308539972801421,
            3.238066845393417,
            6.598303437213936,
            2.659608088311913,
            3.4249044366870596,
            1.7403514603158758,
            9.599864960640643,
            4.369113075109336,
            3.8042018874949823,
            4.28862554376783,
            8.382412615836817,
            3.6774875962929796,
            9.929521288784244,
            5.062516904173654,
            1.295233711655901,
            3.821798499477692,
            2.7132360178922594,
            4.1887507596892,
            3.0863014076672695,
            5.419182763235196,
            2.4073147138874753,
            4.213814051467208,
            1.4251763271057125,
            4.63864991810561,
            1.8216774464333805,
            7.563505390005953,
            3.555241721557862,
            4.493983266747359,
            3.6876604931200268,
            2.454316031896153,
            6.607387606094413,
            9.661565220506915,
            4.386150963483148,
            8.17494299730526,
            6.3198003788420865,
            11.467482765051136,
            4.4065433247605155,
            4.988557708764494,
            9.059519529597852,
            3.3755152852188606,
            7.698743531592137,
            6.618981677978387,
            4.175477039072181,
            6.806814349660336,
            7.431397132641212,
            3.2631899144907384,
            5.0699911933702815,
            2.544651729021151,
            6.080829912290311,
            3.258481663966277,
            3.2938506150971927,
            1.8260826000310022,
            2.7299288206743455,
            6.385202523388037,
            7.33598893338676,
            6.840098800550557,
            2.4629260392399206,
            2.822974223355313,
            4.03397696765668,
            7.1488374756688975,
            3.5278223212553437,
            1.6660138380568987,
            3.458555531243357,
            4.197013547018917,
            4.16975870859787,
            5.562332837781567,
            5.725468422378183,
            6.185468819167943,
            1.6837093971173656,
            4.104228396783036,
            2.1451522407332986,
            5.200237362413139,
            2.084669219978378,
            2.143187834435002,
            2.3140330225844843,
            2.5126535024521788,
            4.292203829906608,
            6.2948386960261375,
            2.707084807447703,
            1.430079063200245,
            3.1398877317179585,
            4.624382313637398,
            2.1098013499095423,
            3.9693834315834158,
            2.918849367120602,
            1.7223877188894419,
            4.2052444383809435,
            1.7027379387189652,
            1.7409058342821224,
            9.73885109912091,
            3.8320864374919834,
            2.8551701405557335,
            3.3985085530029173,
            2.6417203955118143,
            3.1546804210730217,
            3.4473648609972964,
            4.3394904975655955,
            1.8618334598784554,
            1.187119680626635,
            3.3999014822138482,
            2.4231306986854895,
            4.0855931662787555,
            1.1110467622212452,
            3.0836645479450717,
            5.041340866799629,
            2.8584761392149347,
            3.7750089454569724,
            3.98317291735132,
            1.795715129859829,
            2.627605288805403,
            5.475886512190622,
            3.4271418225056154,
            2.891347603259713,
            5.552587133534232,
            5.9809436633734885,
            2.767071076874223,
            4.710760448075293,
            2.376717698170292,
            4.942802828102577,
            4.892581657820345,
            2.6791926869893503,
            2.4743683782040664,
            1.7379083377994877,
            5.778130144567433,
            7.089071071443536,
            2.0388630666174756,
            3.7782560363505686,
            2.0730543304536373,
            2.1948457120426,
            3.361267582386128,
            4.2629464057701245,
            5.958719480369529,
            3.843708413984395,
            3.4815545720306273,
            8.026655051714712,
            6.732224042552772,
            1.4786422506253278,
            3.359516228956052,
            5.051731271772764,
            10.37713698283845,
            1.5329362087999223,
            6.88186935703012,
            2.2867563460747355,
            3.5695696056108077,
            3.930672254899223,
            3.0399623738750345,
            3.209364172407534,
            3.038123333541644,
            1.7884030890335403,
            6.617270158127451,
            4.7892852624442686,
            7.98246376739204,
            8.787993619430916,
            5.4370125075956715,
            5.631570198072755,
            3.7085584331780943,
            1.5882618579464969,
            6.8903268532947735,
            7.892748258525165,
            4.307146684097607,
            5.325141240693638,
            5.615893606452018,
            3.382501768861289,
            2.341364633783292,
            5.952009533729282,
            3.9511068446824225,
            1.878750671506351,
            8.770931877395236,
            7.139918473168391,
            2.968787917613602,
            6.133155615519703,
            4.665146442297377,
            3.5212090189137006,
            4.272327030053521,
            1.8181271956714553,
            1.2111424719202177,
            3.8362418637305393,
            1.6897828694837667,
            4.081067466294491,
            4.733604465939571,
            2.796815803783686,
            3.542465234414504,
            0.9548600743010305,
            2.270717485512143,
            3.850180844042854,
            3.7517662603259643,
            2.9551810686059867,
            5.147502087772008,
            2.4467917467578144,
            5.105007513097308,
            8.408655228616226,
            4.452315949253354,
            1.6886214201468253,
            3.2675270705264667,
            5.874045148360757,
            3.0273104135176405,
            3.7648099268073,
            8.321357616729175,
            6.922623016214074,
            1.678742522381662,
            2.631473336660425,
            1.9769260252425107,
            8.54513049934888,
            5.606712381312642,
            9.985899928272206,
            7.7550045545739374,
            2.5269986968302973,
            8.642277130729743,
            3.817375807234058,
            5.305481369727255,
            8.292051764966633,
            4.4453842092352,
            4.349003461844681,
            6.000704722101477,
            4.543551953871819,
            5.833324443935714,
            3.3688153740004076,
            1.1431546305228522,
            5.498587072359388,
            2.703070385560106,
            3.9334169401183137,
            9.643230396608962,
            3.4603187156827,
            7.691852031027734,
            7.22790045250162,
            4.393578726180066,
            5.243417451059749,
            8.13349302370389,
            5.2583811234088245,
            3.546269300522664,
            3.506822851905734,
            8.301287815488369,
            4.300791041878211,
            6.6587500730057325,
            3.9709462505155106,
            3.362464817847712,
            11.37915331454732,
            8.068138598813995,
            7.916480638467263,
            5.12202392506206,
            8.685405827800933,
            3.6092106401749424,
            5.91911663192843,
            9.953524458486692,
            4.472235413408162,
            5.318791897610933,
            5.403579710683875,
            10.195682092743064,
            14.299048571158615,
            12.38411860673397,
            2.1620234107802583,
            8.521664591626005,
            12.330847080589812,
            3.136419959075777,
            5.542427971288237,
            4.941578224663928,
            4.725295110261525,
            5.165968526411404,
            6.110105031454435,
        ];

        for (actual, expected) in daily_hwd.unwrap().iter().zip(expected_result.iter()) {
            assert_relative_eq!(actual, expected, max_relative = 1e-4)
        }
    }

    #[rstest]
    fn test_edit_hot_water_distribution(mut test_input: InputForProcessing) {
        let tfa = calc_tfa(&test_input).unwrap();
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

        let actual_hot_water_distribution_inner = test_input
            .water_distribution()
            .unwrap()
            .unwrap()
            .first()
            .cloned()
            .unwrap();

        assert_eq!(
            actual_hot_water_distribution_inner,
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
        let energy_supply_key = ENERGY_SUPPLY_NAME_ELECTRICITY;
        let _ = test_input.add_diverter_to_energy_supply(energy_supply_key, diverter);

        remove_pv_diverter_if_present(&mut test_input).unwrap();
        let energy_supply = test_input
            .energy_supply_by_key(energy_supply_key)
            .unwrap()
            .unwrap();
        assert!(!energy_supply.contains_key("diverter"))
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
            "battery_location": "inside",
            "grid_charging_possible": false
        });
        let energy_supply_key = ENERGY_SUPPLY_NAME_ELECTRICITY;
        let _ =
            test_input.add_electric_battery_to_energy_supply(energy_supply_key, electric_battery);

        let _ = remove_electric_battery_if_present(&mut test_input);
        let energy_supply = test_input
            .energy_supply_by_key(energy_supply_key)
            .unwrap()
            .unwrap();
        assert!(!energy_supply.contains_key("ElectricBattery"));
    }

    #[rstest]
    fn test_edit_space_cool_system(mut test_input: InputForProcessing) {
        test_input.set_part_o_active_cooling_required(true).unwrap();
        let _ = edit_space_cool_system(&mut test_input);
        let space_cool_system = test_input.space_cool_system().unwrap().unwrap();

        for system in space_cool_system.values() {
            assert_eq!(system.get("efficiency").and_then(|v| v.as_f64()), Some(5.1));
            assert_eq!(
                system.get("frac_convective").and_then(|v| v.as_f64()),
                Some(0.95)
            );
            assert_eq!(
                system.get("EnergySupply").and_then(|v| v.as_str()),
                Some(ENERGY_SUPPLY_NAME_ELECTRICITY)
            );
        }
    }

    // this test does not exist in Python HEM
    #[rstest]
    #[ignore = "This currently fails because our test data does not have all expected fields on ExternalConditions."]
    fn test_design_capacity(test_input: InputForProcessing) {
        let actual_design_capacity = calc_design_capacity(&test_input).unwrap();
        assert_eq!(
            actual_design_capacity.0.get("zone 1").unwrap(),
            &5.356813765662826
        );
        assert_eq!(
            actual_design_capacity.0.get("zone 2").unwrap(),
            &5.356813765662826
        );
        assert_eq!(actual_design_capacity.1, 10.713627531325653);
    }

    // this test does not exist in Python HEM
    #[rstest]
    fn test_initialise_temperature_setpoints(mut test_input: InputForProcessing) {
        initialise_temperature_setpoints(&mut test_input).unwrap();

        let temp_setpoints = test_input.all_init_temp_setpoints().unwrap();

        for temp_setpoint in temp_setpoints {
            assert_eq!(temp_setpoint, Some(18.));
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

        remove_onsite_generation_if_present(&mut test_input).unwrap();
        assert!(test_input.on_site_generation().unwrap().is_none());
    }

    #[rstest]
    fn test_add_solar_pv_house_only(mut test_input: InputForProcessing) {
        let expected_result = json!({"PV1": {
                "EnergySupply": "mains elec",
                "orientation360": 180.,
                "peak_power": 4.444444444444445,
                "inverter_peak_power_ac": 4.444444444444445,
                "inverter_peak_power_dc": 4.444444444444445,
                "inverter_is_inside": false,
                "inverter_type": "optimised_inverter",
                "pitch": 45.,
                "type": "PhotovoltaicSystem",
                "ventilation_strategy": "moderately_ventilated",
                "shading": [],
                "base_height": 1.,
                "width": 6.324555320336759,
                "height": 3.1622776601683795
                }
        })
        .as_object()
        .cloned()
        .unwrap();

        let is_notional_a = true;
        let is_fee = false;
        let total_floor_area = calc_tfa(&test_input).unwrap();

        add_solar_pv(&mut test_input, is_notional_a, is_fee, total_floor_area).unwrap();

        assert_eq!(
            test_input.on_site_generation().unwrap().unwrap().to_owned(),
            expected_result
        );
    }

    #[rstest]
    fn test_edit_add_default_space_heating_system(mut test_input: InputForProcessing) {
        let expected: IndexMap<String, HeatSourceWetDetails> = serde_json::from_value(json!(
         {
            "notional_HP": {
                "EnergySupply": "mains elec",
                "backup_ctrl_type": "TopUp",
                "min_modulation_rate_35": 0.4,
                "min_modulation_rate_55": 0.4,
                "min_temp_diff_flow_return_for_hp_to_operate": 0,
                "modulating_control": true,
                "power_crankcase_heater": 0.01,
                "power_heating_circ_pump": 0.022866,
                "power_max_backup": 3,
                "power_off": 0,
                "power_source_circ_pump": 0.01,
                "power_standby": 0.01,
                "sink_type": "Water",
                "source_type": "OutsideAir",
                "temp_lower_operating_limit": -10,
                "temp_return_feed_max": 60,
                "test_data_EN14825": [
                 {
                     "capacity": 7.4,
                     "cop": 2.79,
                     "degradation_coeff": 0.9,
                     "design_flow_temp": 35,
                     "temp_outlet": 34,
                     "temp_source": -7,
                     "temp_test": -7,
                     "test_letter": "A"
                 },
                 {
                     "capacity": 4.588,
                     "cop": 4.29,
                     "degradation_coeff": 0.9,
                     "design_flow_temp": 35,
                     "temp_outlet": 30,
                     "temp_source": 2,
                     "temp_test": 2,
                     "test_letter": "B"
                 },
                 {
                     "capacity": 4.07,
                     "cop": 5.91,
                     "degradation_coeff": 0.9,
                     "design_flow_temp": 35,
                     "temp_outlet": 27,
                     "temp_source": 7,
                     "temp_test": 7,
                     "test_letter": "C"
                 },
                 {
                     "capacity": 3.478,
                     "cop": 8.02,
                     "degradation_coeff": 0.9,
                     "design_flow_temp": 35,
                     "temp_outlet": 24,
                     "temp_source": 12,
                     "temp_test": 12,
                     "test_letter": "D"
                 },
                 {
                     "capacity": 7.77,
                     "cop": 2.49,
                     "degradation_coeff": 0.9,
                     "design_flow_temp": 35,
                     "temp_outlet": 35,
                     "temp_source": -10,
                     "temp_test": -10,
                     "test_letter": "F"
                 },
                 {
                     "capacity": 7.326,
                     "cop": 2.03,
                     "degradation_coeff": 0.9,
                     "design_flow_temp": 55,
                     "temp_outlet": 52,
                     "temp_source": -7,
                     "temp_test": -7,
                     "test_letter": "A"
                 },
                 {
                     "capacity": 4.44,
                     "cop": 3.12,
                     "degradation_coeff": 0.9,
                     "design_flow_temp": 55,
                     "temp_outlet": 42,
                     "temp_source": 2,
                     "temp_test": 2,
                     "test_letter": "B"
                 },
                 {
                     "capacity": 3.626,
                     "cop": 4.41,
                     "degradation_coeff": 0.9,
                     "design_flow_temp": 55,
                     "temp_outlet": 36,
                     "temp_source": 7,
                     "temp_test": 7,
                     "test_letter": "C"
                 },
                 {
                     "capacity": 3.774,
                     "cop": 6.30,
                     "degradation_coeff": 0.9,
                     "design_flow_temp": 55,
                     "temp_outlet": 30,
                     "temp_source": 12,
                     "temp_test": 12,
                     "test_letter": "D"
                 },
                 {
                     "capacity": 7.622,
                     "cop": 1.87,
                     "degradation_coeff": 0.9,
                     "design_flow_temp": 55,
                     "temp_outlet": 55,
                     "temp_source": -10,
                     "temp_test": -10,
                     "test_letter": "F"
                 }
             ],
                "time_constant_onoff_operation": 120,
                "time_delay_backup": 1,
                "type": "HeatPump",
                "var_flow_temp_ctrl_during_test": true
            }
        }))
        .unwrap();

        let design_capacity_overall = 7.4;
        edit_add_default_space_heating_system(&mut test_input, design_capacity_overall).unwrap();

        // TODO currently because of a custom impl PartialEq for HeatPumpTestDatum
        // this test passes when test data does not match exactly

        assert_eq!(test_input.heat_source_wet().unwrap(), expected);
    }

    // this test does not exist in Python HEM
    #[rstest]
    fn test_edit_default_space_heating_distribution_system(mut test_input: InputForProcessing) {
        let design_capacity: IndexMap<String, f64> =
            serde_json::from_value(json!({"zone 1": 0., "zone 2": 0})).unwrap();

        edit_default_space_heating_distribution_system(&mut test_input, &design_capacity).unwrap();

        for zone_key in test_input.zone_keys().unwrap() {
            let expected_space_heat_system_name = zone_key.clone() + "_SpaceHeatSystem_Notional";

            let actual_space_heat_system_name_in_zone = test_input
                .space_heat_system_for_zone(&zone_key)
                .unwrap()
                .first()
                .unwrap()
                .to_owned();

            let actual_space_heat_system = test_input
                .space_heat_system_for_key(&expected_space_heat_system_name)
                .unwrap();

            let expected_space_heat_systems: Value = json!({
                "zone 1_SpaceHeatSystem_Notional":
                {
                    "Control": "HeatingPattern_Null",
                    "HeatSource": {"name": "hp", "temp_flow_limit_upper": 65.0},
                    "Zone": "zone 1",
                    "advanced_start": 1,
                    "design_flow_temp": 45,
                    "design_flow_rate": 12.,
                    "ecodesign_controller": {
                        "ecodesign_control_class": 2,
                        "max_outdoor_temp": 20,
                        "min_flow_temp": 21,
                        "min_outdoor_temp": 0},
                    "temp_diff_emit_dsgn": 5,
                    "temp_setback": 18,
                    "thermal_mass": 0.0,
                    "emitters": [{"wet_emitter_type": "radiator",
                          "frac_convective": 0.7,
                          "c": 0.,
                          "n": 1.34}],
                    "type": "WetDistribution"
                },
                "zone 2_SpaceHeatSystem_Notional":
                {
                    "Control": "HeatingPattern_Null",
                    "HeatSource": {"name": "hp", "temp_flow_limit_upper": 65.0},
                    "Zone": "zone 2",
                    "advanced_start": 1,
                    "design_flow_temp": 45,
                    "design_flow_rate": 12.,
                    "ecodesign_controller": {
                        "ecodesign_control_class": 2,
                        "max_outdoor_temp": 20,
                        "min_flow_temp": 21,
                        "min_outdoor_temp": 0},
                    "temp_diff_emit_dsgn": 5,
                    "temp_setback": 18,
                    "thermal_mass": 0.0,
                    "emitters": [{"wet_emitter_type": "radiator",
                          "frac_convective": 0.7,
                          "c": 0.,
                          "n": 1.34}],
                    "type": "WetDistribution"
                }
            }
            );

            let expected_space_heat_system = expected_space_heat_systems
                .get::<&std::string::String>(&expected_space_heat_system_name.clone().into())
                .unwrap();

            assert_eq!(
                actual_space_heat_system_name_in_zone,
                expected_space_heat_system_name
            );
            assert_eq!(actual_space_heat_system, Some(expected_space_heat_system))
        }
    }
}
