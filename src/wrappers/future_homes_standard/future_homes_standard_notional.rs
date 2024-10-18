#![allow(dead_code)]

use crate::{
    compare_floats::max_of_2,
    core::{
        space_heat_demand::building_element::{pitch_class, HeatFlowDirection},
        units::{LITRES_PER_CUBIC_METRE, SECONDS_PER_HOUR},
    },
    input::{
        GroundBuildingElement, HeatPumpSourceType,
        HeatSourceWetDetails::{HeatPump, Hiu},
        InputForProcessing, ThermalBridgingDetails, UValueEditableBuildingElement,
    },
    wrappers::future_homes_standard::future_homes_standard::calc_tfa,
};
use anyhow::{anyhow, bail};
use itertools::Itertools;
use std::collections::HashMap;
use std::sync::LazyLock;

/// Apply assumptions and pre-processing steps for the Future Homes Standard Notional building
fn apply_fhs_not_preprocessing(
    mut input: InputForProcessing,
    fhs_not_a_assumptions: bool,
    _fhs_not_b_assumptions: bool,
    fhs_fee_not_a_assumptions: bool,
    fhs_fee_not_b_assumptions: bool,
) {
    let _is_not_a = fhs_not_a_assumptions || fhs_fee_not_a_assumptions;
    let _is_fee = fhs_fee_not_a_assumptions || fhs_fee_not_b_assumptions;
    // Check if a heat network is present
    let _is_heat_network = check_heatnetwork_present(&input);

    // Determine cold water source
    let cold_water_type = input.cold_water_source();

    let _cold_water_source = match (cold_water_type.mains_water, cold_water_type.header_tank) {
        (Some(source), None) => source,
        (None, Some(source)) => source,
        _ => panic!("Error: There should be exactly one cold water type"),
    };

    // Retrieve the number of bedrooms and total volume
    let _bedroom_number = input.number_of_bedrooms();
    // Loop through zones to sum up volume.
    let _total_volume = input.total_zone_area();

    // Determine the TFA
    let _tfa = calc_tfa(&input);
    edit_lighting_efficacy(&mut input);

    todo!()
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
/// TODO - how to deal with ThermalBridging when lengths are not specified?
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

/// Calculate effective air change rate accoring to according to Part F 1.24 a
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

#[cfg(test)]
mod tests {
    use crate::core::space_heat_demand::building_element::{pitch_class, HeatFlowDirection};

    use super::*;
    use crate::input;
    use crate::input::{ThermalBridging, ThermalBridgingDetails};
    use rstest::{fixture, rstest};
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
    // test written in Rust, not present in Python
    fn test_apply_fhs_not_preprocessing(test_input: InputForProcessing) {
        let fhs_not_a_assumptions = true;
        let _fhs_not_b_assumptions = false;
        let fhs_fee_not_a_assumptions = false;
        let fhs_fee_not_b_assumptions = false;

        let _actual = apply_fhs_not_preprocessing(
            test_input,
            fhs_not_a_assumptions,
            _fhs_not_b_assumptions,
            fhs_fee_not_a_assumptions,
            fhs_fee_not_b_assumptions,
        );
        todo!()
    }

    #[rstest]
    /// test written in Rust, not present in Python
    fn test_check_heatnetwork_present(test_input: InputForProcessing) {
        assert_eq!(check_heatnetwork_present(&test_input), false);
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

        for building_element in test_input.all_building_elements() {
            if let input::BuildingElement::Opaque { .. }
            | input::BuildingElement::AdjacentZTUSimple { .. } = building_element
            {
                if let Some(u_value) = building_element.u_value() {
                    match pitch_class(building_element.pitch()) {
                        HeatFlowDirection::Downwards => {
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

    #[rstest]
    fn test_edit_ground_floors(mut test_input: InputForProcessing) {
        let test_input = test_input.borrow_mut();

        edit_ground_floors(test_input).unwrap();

        for building_element in test_input.all_building_elements() {
            if let input::BuildingElement::Ground {
                u_value,
                r_f,
                psi_wall_floor_junc,
                ..
            } = building_element
            {
                assert_eq!(*u_value, 0.13);
                assert_eq!(*r_f, 6.12);
                assert_eq!(*psi_wall_floor_junc, 0.16);
            }
        }
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
}
