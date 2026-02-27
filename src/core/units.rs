use crate::compare_floats::min_of_2;
use fsum::FSum;
use serde::{Deserialize, Serialize};
use serde_valid::Validate;
use std::fmt::Display;
use std::str::FromStr;
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
pub(crate) const KNOTS_PER_METRES_PER_SECOND: f64 = 1. / (1852. / 3600.);
pub const MILLIMETRES_IN_METRE: u32 = 1_000;

pub(crate) fn average_monthly_to_annual(list_monthly_averages: [f64; 12]) -> f64 {
    FSum::with_all(
        list_monthly_averages
            .iter()
            .enumerate()
            .map(|(month_idx, month_ave)| month_ave * DAYS_IN_MONTH[month_idx] as f64),
    )
    .value()
        / DAYS_IN_MONTH.iter().sum::<u32>() as f64
}

// Python equivalent of this has an allow_none parameter. We know when we have Options in Rust so we
// can check for Some value in calling code - no need to have defensive code/ extra params here.
pub(crate) fn celsius_to_kelvin(temp_c: f64) -> Result<f64, BelowAbsoluteZeroError> {
    if temp_c < -273.15 {
        Err(BelowAbsoluteZeroError::from_c(temp_c))
    } else {
        Ok(temp_c + 273.15)
    }
}

pub(crate) fn kelvin_to_celsius(temp_k: f64) -> Result<f64, BelowAbsoluteZeroError> {
    if temp_k < 0.0 {
        Err(BelowAbsoluteZeroError::from_k(temp_k))
    } else {
        Ok(temp_k - 273.15)
    }
}

pub fn convert_profile_to_daily(timestep_totals: &[f64], timestep: f64) -> Vec<f64> {
    let total_steps = timestep_totals.len();
    let steps_per_day = (HOURS_PER_DAY as f64 / timestep).floor() as usize;
    (0..total_steps)
        .step_by(steps_per_day)
        .map(|y| {
            FSum::with_all(
                timestep_totals[y..min_of_2(y + steps_per_day, timestep_totals.len())].iter(),
            )
            .value()
        })
        .collect()
}

pub(crate) fn calculate_thermal_resistance_of_virtual_layer(
    u_value: f64,
    thermal_resistance_floor_construction: f64,
) -> Result<f64, BadThermalResistanceCalculationError> {
    // Thermal properties of ground from BS EN ISO 13370:2017 Table 7
    // Use values for clay or silt (same as BR 443 and SAP 10)
    let thermal_conductivity = 1.5; // in W/(m.K)

    // Calculate thermal resistance and heat capacity of fixed ground layer
    // using BS EN ISO 13370:2017
    let thickness_ground_layer = 0.5; // in m. Specified in BS EN ISO 52016-1:2017 section 6.5.8.2

    // thermal resistance in (m2.K)/W
    let r_gr = thickness_ground_layer / thermal_conductivity;

    // Calculate thermal resistance of virtual layer using BS EN ISO 13370:2017 Equation (F1)
    let r_si = 0.17; // ISO 6946 - internal surface resistance
    let r_vi = (1.0 / u_value) - r_si - thermal_resistance_floor_construction - r_gr; // in m2.K/W

    // BS EN ISO 13370:2017 Table 2 validity interval r_vi > 0
    match r_vi {
        r_vi if r_vi <= 0. => Err(BadThermalResistanceCalculationError),
        _ => Ok(r_vi),
    }
}

#[derive(Debug, Error)]
#[error("A temperature of {k}ºK/{}ºC was encountered, which is less than absolute zero", k - 273.15)]
pub(crate) struct BelowAbsoluteZeroError {
    k: f64,
}

#[derive(Debug, Error)]
#[error("Thermal resistance (r_vi) calculation produced a negative value. Check u-value and thermal_resistance_floor_construction inputs for floors.")]
pub(crate) struct BadThermalResistanceCalculationError;

impl BelowAbsoluteZeroError {
    fn from_k(k: f64) -> Self {
        Self { k }
    }

    fn from_c(c: f64) -> Self {
        Self { k: c + 273.15 }
    }
}

#[derive(Clone, Copy, Debug, Default, Deserialize, PartialEq, PartialOrd, Serialize, Validate)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(transparent)]
#[repr(transparent)]
pub struct Orientation360(
    #[validate(minimum = 0.)]
    #[validate(maximum = 360.)]
    f64,
);

impl Orientation360 {
    pub(crate) fn new(angle: f64) -> Result<Self, Orientation360Error> {
        if !(0. ..=360.).contains(&angle) {
            return Err(Orientation360Error::InvalidAngle);
        }

        Ok(Self(angle))
    }

    pub(crate) fn angle(&self) -> f64 {
        self.0
    }

    pub(crate) fn transform_to_180(&self) -> f64 {
        180. - self.0
    }

    pub(crate) fn create_from_180(angle180: f64) -> Result<Self, Orientation360Error> {
        if !(-180. ..=180.).contains(&angle180) {
            return Err(Orientation360Error::InvalidAngleFrom180);
        }
        Ok(Self(180. - angle180))
    }

    // Python __str__ method is covered equivalently by Display trait impl

    pub(crate) fn orientation_difference(
        orientation1: Orientation360,
        orientation2: Orientation360,
    ) -> f64 {
        let op_rel_orientation = (orientation1.angle() - orientation2.angle()).abs();

        if op_rel_orientation > 180. {
            360. - op_rel_orientation
        } else {
            op_rel_orientation
        }
    }
}

impl Display for Orientation360 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl From<f64> for Orientation360 {
    fn from(angle: f64) -> Self {
        Self::new(angle).expect("Angle must be between 0 and 360 degrees inclusive")
    }
}

impl FromStr for Orientation360 {
    type Err = Orientation360Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let angle = s
            .parse::<f64>()
            .map_err(|_| Orientation360Error::InvalidAngle)?;
        Self::new(angle).map_err(|_| Orientation360Error::InvalidAngle)
    }
}

#[derive(Clone, Copy, Debug, Error)]
pub enum Orientation360Error {
    #[error("Angle must be between 0 and 360 degrees inclusive")]
    InvalidAngle,
    #[error("Angle for create_from_180 must be between -180 and 180 degrees inclusive")]
    InvalidAngleFrom180,
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[rstest]
    fn should_do_correct_temperature_conversions() {
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
    fn should_convert_average_monthly_to_annual() {
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
    fn should_convert_profile_to_daily() {
        let mut list_timestep_totals = vec![1.0; 48];
        list_timestep_totals.extend((0..48).map(|x| x as f64 / 2.0));
        assert_eq!(
            convert_profile_to_daily(&list_timestep_totals, 0.5),
            vec![48.0, 564.0],
            "incorrect conversion of per-timestep profile to daily profile"
        );
    }

    mod orientation360 {
        use super::*;
        use pretty_assertions::assert_eq;

        #[rstest]
        fn test_orientation360_angle() {
            assert_eq!(Orientation360::new(180.).unwrap().angle(), 180.);
        }

        #[rstest]
        fn test_orientation360_invalid_angle() {
            assert!(Orientation360::new(-10.).is_err());
            assert!(Orientation360::new(380.).is_err());
        }

        #[rstest]
        fn test_orientation360_str() {
            assert_eq!(format!("{}", Orientation360(180.)), "180");
        }

        #[rstest]
        #[case(0., 180.)]
        #[case(90., 90.)]
        #[case(180., 0.)]
        #[case(270., -90.)]
        #[case(360., -180.)]
        fn test_orientation360_transform_to_180(#[case] value: f64, #[case] expected_result: f64) {
            assert_eq!(Orientation360(value).transform_to_180(), expected_result);
        }

        #[rstest]
        #[case(-90., 270.)]
        #[case(90., 90.)]
        #[case(-180., 360.)]
        #[case(0., 180.)]
        #[case(180., 0.)]
        fn test_orientation360_create_from_180(#[case] value: f64, #[case] expected_result: f64) {
            let orientation360 = Orientation360::create_from_180(value).unwrap();
            assert_eq!(orientation360.angle(), expected_result);
        }

        #[rstest]
        fn test_orientation360_create_from_180_invalid_angle() {
            assert!(Orientation360::create_from_180(190.).is_err());
            assert!(Orientation360::create_from_180(-190.).is_err());
        }
    }
}
