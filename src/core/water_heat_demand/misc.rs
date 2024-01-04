use crate::core::material_properties::WATER;

/// Calculate the fraction of hot water required when mixing hot and cold
/// water to achieve a target temperature
///
/// Arguments:
/// * `temp_target` -- temperature to be achieved, in any units
/// * `temp_hot`    -- temperature of hot water to be mixed, in same units as temp_target
/// * `temp_cold`   -- temperature of cold water to be mixed, in same units as temp_target
pub fn frac_hot_water(temp_target: f64, temp_hot: f64, temp_cold: f64) -> f64 {
    (temp_target - temp_cold) / (temp_hot - temp_cold)
}

/// Calculates the kWh energy content of the hot water demand.
///          
/// Arguments:
/// * `litres_demand`  -- hot water demand in litres
/// * `demand_temp`    -- temperature of hot water inside the pipe, in degrees C
/// * `cold_temp`     -- temperature outside the pipe, in degrees C
pub fn water_demand_to_kWh(litres_demand: f64, demand_temp: f64, cold_temp: f64) -> f64 {
    &WATER.volumetric_energy_content_kWh_per_litre(demand_temp, cold_temp) * litres_demand
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

    #[rstest]
    pub fn should_calculate_correct_frac_hot_water() {
        assert_eq!(
            frac_hot_water(40.0, 55.0, 5.0),
            0.7,
            "incorrect fraction of hot water returned"
        );
    }

    #[rstest]
    pub fn should_calculate_correct_water_demand_to_kWh() {
        let litres_demand = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0];
        let demand_temp = [40.0, 35.0, 37.0, 39.0, 40.0, 38.0, 39.0, 40.0];
        let cold_temp = [5.0, 4.0, 5.0, 6.0, 5.0, 4.0, 3.0, 4.0];
        for i in 0..8 {
            assert_eq!(
                round_by_precision(
                    water_demand_to_kWh(litres_demand[i], demand_temp[i], cold_temp[i]),
                    1e5
                ),
                round_by_precision(
                    [0.20339, 0.36029, 0.55787, 0.76707, 1.01694, 1.18547, 1.46440, 1.6736][i],
                    1e5
                ),
                "incorrect water demand to kWh returned"
            );
        }
    }

    fn round_by_precision(src: f64, precision: f64) -> f64 {
        (precision * src).round() / precision
    }
}