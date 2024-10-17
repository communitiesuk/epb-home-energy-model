#![allow(dead_code)]
use core::panic;

use itertools::Itertools;

use crate::{
    compare_floats::max_of_2,
    core::units::{LITRES_PER_CUBIC_METRE, SECONDS_PER_HOUR},
    input::{
        HeatPumpSourceType,
        HeatSourceWetDetails::{HeatPump, Hiu},
        InputForProcessing,
    },
    wrappers::future_homes_standard::future_homes_standard::calc_tfa,
};

/// Apply assumptions and pre-processing steps for the Future Homes Standard Notional building
fn apply_fhs_not_preprocessing(
    mut input: InputForProcessing,
    fhs_not_a_assumptions: bool,
    _fhs_not_b_assumptions: bool,
    fhs_fee_not_a_assumptions: bool,
    fhs_fee_not_b_assumptions: bool,
) -> () {
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
fn edit_lighting_efficacy(input: &mut InputForProcessing) -> () {
    let lighting_efficacy = 120.0;
    input.set_lighting_efficacy_for_all_zones(lighting_efficacy);
}

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
    use super::*;
    use rstest::{fixture, rstest};
    use std::borrow::BorrowMut;
    use std::fs::File;
    use std::io::BufReader;
    use std::path::Path;

    #[fixture]
    fn test_input() -> InputForProcessing {
        let reader = BufReader::new(File::open(Path::new("src/wrappers/future_homes_standard/test_future_homes_standard_notional_input_data.json")).unwrap());
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
}
