use crate::compare_floats::min_of_2;
use thiserror::Error;

pub const JOULES_PER_KILOWATT_HOUR: u32 = 3_600_000;
pub const KILOJOULES_PER_KILOWATT_HOUR: u32 = 3_600;
pub const JOULES_PER_KILOJOULE: u32 = 1_000;
pub const WATTS_PER_KILOWATT: u32 = 1_000;
pub const LITRES_PER_CUBIC_METRE: u32 = 1_000;
pub(crate) const _M3_PER_S_TO_L_PER_MIN: f64 = 60_000.;
pub const MINUTES_PER_HOUR: u32 = 60;
pub const SECONDS_PER_MINUTE: u32 = 60;
pub const SECONDS_PER_HOUR: u32 = 3_600;
pub const HOURS_PER_DAY: u32 = 24;
pub const DAYS_PER_YEAR: u32 = 365;
pub const DAYS_IN_MONTH: [u32; 12] = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
pub const MILLIMETRES_IN_METRE: u32 = 1_000;

pub fn average_monthly_to_annual(list_monthly_averages: [f64; 12]) -> f64 {
    list_monthly_averages
        .iter()
        .enumerate()
        .map(|(month_idx, month_ave)| month_ave * DAYS_IN_MONTH[month_idx] as f64)
        .sum::<f64>()
        / DAYS_IN_MONTH.iter().sum::<u32>() as f64
}

pub fn convert_profile_to_daily(timestep_totals: &[f64], timestep: f64) -> Vec<f64> {
    let total_steps = timestep_totals.len();
    let steps_per_day = (HOURS_PER_DAY as f64 / timestep).floor() as usize;
    (0..total_steps)
        .step_by(steps_per_day)
        .map(|y| {
            timestep_totals[y..min_of_2(y + steps_per_day, timestep_totals.len())]
                .iter()
                .sum::<f64>()
        })
        .collect()
}

// Python equivalent of this has an allow_none parameter. We know when we have Options in Rust so we
// can check for Some value in calling code - no need to have defensive code/ extra params here.
pub(crate) fn kelvin_to_celsius(temp_k: f64) -> Result<f64, BelowAbsoluteZeroError> {
    if temp_k < 0.0 {
        Err(BelowAbsoluteZeroError::from_k(temp_k))
    } else {
        Ok(temp_k - 273.15)
    }
}

pub(crate) fn celsius_to_kelvin(temp_c: f64) -> Result<f64, BelowAbsoluteZeroError> {
    if temp_c < -273.15 {
        Err(BelowAbsoluteZeroError::from_c(temp_c))
    } else {
        Ok(temp_c + 273.15)
    }
}

#[derive(Debug, Error)]
#[error("A temperature of {k}ºK/{}ºC was encountered, which is less than absolute zero", k - 273.15)]
pub(crate) struct BelowAbsoluteZeroError {
    k: f64,
}

impl BelowAbsoluteZeroError {
    fn from_k(k: f64) -> Self {
        Self { k }
    }

    fn from_c(c: f64) -> Self {
        Self { k: c + 273.15 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[rstest]
    pub fn should_do_correct_temperature_conversions() {
        assert_eq!(
            celsius_to_kelvin(20.0).unwrap(),
            293.15,
            "incorrect conversion of Celsius to Kelvin"
        );
        assert_eq!(
            kelvin_to_celsius(5.0).unwrap(),
            -268.15,
            "incorrect conversion to Kelvin to Celsius"
        );
        for i in -10..80 {
            assert_eq!(
                kelvin_to_celsius(celsius_to_kelvin(i as f64).unwrap()).unwrap(),
                i as f64,
                "round trip temperature conversion (C to K to C) failed to return orig value"
            );
        }
    }

    #[rstest]
    pub fn should_convert_average_monthly_to_annual() {
        let list_monthly_averages = [
            4.3, 4.9, 6.5, 8.9, 11.7, 14.6, 16.6, 16.4, 14.1, 10.6, 7.1, 4.2,
        ];
        assert_eq!(
            average_monthly_to_annual(list_monthly_averages),
            10.020547945205479,
            "incorrect conversion of monthly averages to annual average"
        );
    }

    #[rstest]
    pub fn should_convert_profile_to_daily() {
        let mut list_timestep_totals = vec![1.0; 48];
        list_timestep_totals.extend((0..48).map(|x| x as f64 / 2.0));
        assert_eq!(
            convert_profile_to_daily(&list_timestep_totals, 0.5),
            vec![48.0, 564.0],
            "incorrect conversion of per-timestep profile to daily profile"
        );
    }
}
