use crate::core::material_properties::WATER;
use crate::core::schedule::WaterScheduleEventType;
use crate::hem_core::simulation_time::SimulationTimeIteration;
use anyhow::bail;
use format_num::format_num;
use fsum::FSum;
use itertools::Itertools;
use smartstring::alias::String;

// Fraction of domestic hot water energy that becomes internal gains
// This applies to both hot water usage and combi boiler losses
pub(crate) const FRAC_DHW_ENERGY_INTERNAL_GAINS: f64 = 0.25;

pub(crate) type CallableGetHotWaterTemperature = Box<dyn Fn(f64) -> anyhow::Result<f64>>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum WaterEventResultType {
    Shower,
    Bath,
    Other,
    PipeFlush,
}

impl WaterEventResultType {
    fn abbreviation(&self) -> &'static str {
        match self {
            Self::Shower => "S",
            Self::Bath => "B",
            Self::Other => "O",
            Self::PipeFlush => "P",
        }
    }
}

impl From<WaterScheduleEventType> for WaterEventResultType {
    fn from(value: WaterScheduleEventType) -> Self {
        match value {
            WaterScheduleEventType::Shower => WaterEventResultType::Shower,
            WaterScheduleEventType::Bath => WaterEventResultType::Bath,
            WaterScheduleEventType::Other => WaterEventResultType::Other,
        }
    }
}

#[derive(Clone, Copy, Debug)]
/// Result of processing a single water use event
pub(crate) struct WaterEventResult {
    pub(crate) event_result_type: WaterEventResultType,
    pub(crate) temperature_warm: f64, // Temperature of water at outlet (Celsius)
    pub(crate) volume_warm: f64,      // Volume of water at outlet (litres)
    pub(crate) volume_hot: f64,       // Hot water demand volume (litres)
    pub(crate) event_duration: f64,   // Duration of hot water event (minutes)
}

impl WaterEventResult {
    pub(crate) fn format_event(&self) -> String {
        let hot_volume = self.volume_hot;
        let warm_volume = self.volume_warm;
        let temperature = self.temperature_warm;
        let abbrev = self.event_result_type.abbreviation();

        format!(
            "{abbrev}: {} ({} @ {})",
            format_num!(".10g", hot_volume),
            format_num!(".10g", warm_volume),
            format_num!(".10g", temperature)
        )
        .into()
    }
}

pub(crate) fn summarise_events(events: &[WaterEventResult]) -> String {
    events
        .iter()
        .map(|event| event.format_event())
        .join(" | ")
        .into()
}

#[derive(Debug, thiserror::Error)]
#[error("Cannot achieve temperature {temperature_target} by mixing water at {temperature_cold} and {temperature_hot}")]
pub(crate) struct UnachievableTemperatureByMixingError {
    temperature_target: f64,
    temperature_hot: f64,
    temperature_cold: f64,
}

/// Calculate the fraction of hot water required when mixing hot and cold
/// water to achieve a target temperature
///
/// Arguments:
/// * `temp_target` -- temperature to be achieved, in any units
/// * `temp_hot`    -- temperature of hot water to be mixed, in same units as temp_target
/// * `temp_cold`   -- temperature of cold water to be mixed, in same units as temp_target
pub(crate) fn calc_fraction_hot_water(
    temperature_target: f64,
    temperature_hot: f64,
    temperature_cold: f64,
) -> Result<f64, UnachievableTemperatureByMixingError> {
    let fraction = (temperature_target - temperature_cold) / (temperature_hot - temperature_cold);
    if (fraction < 0.0 && !is_close!(fraction, 0.0, abs_tol = 1e-10))
        || (fraction > 1.0 && !is_close!(fraction, 1.0, abs_tol = 1e-10))
    {
        return Err(UnachievableTemperatureByMixingError {
            temperature_target,
            temperature_hot,
            temperature_cold,
        });
    }

    Ok(fraction)
}

/// Calculates the kWh energy content of the hot water demand.
///          
/// Arguments:
/// * `litres_demand`  -- hot water demand in litres
/// * `demand_temp`    -- temperature of hot water inside the pipe, in degrees C
/// * `cold_temp`     -- temperature outside the pipe, in degrees C
pub fn water_demand_to_kwh(
    litres_demand: f64,
    demand_temperature: f64,
    cold_temperature: f64,
) -> f64 {
    WATER.volumetric_energy_content_kwh_per_litre(demand_temperature, cold_temperature)
        * litres_demand
}

pub(crate) fn volume_hot_water_required<
    T: Fn(f64) -> anyhow::Result<f64>,
    U: Fn(f64, SimulationTimeIteration) -> anyhow::Result<Vec<(f64, f64)>>,
>(
    volume_warm_water: f64,
    temperature_target: f64,
    func_temperature_hot_water: T,
    func_temperature_cold_water: U,
    simtime: SimulationTimeIteration,
) -> anyhow::Result<Option<f64>> {
    let mut temperature_hot_water = func_temperature_hot_water(volume_warm_water * 0.5)?;
    let mut list_temperature_volume =
        func_temperature_cold_water(volume_warm_water * 0.5, simtime)?;
    let mut temperature_cold_water =
        FSum::with_all(list_temperature_volume.iter().map(|(t, v)| t * v)).value()
            / FSum::with_all(list_temperature_volume.iter().map(|(_, v)| v)).value();
    let mut temperature_warm_water = temperature_hot_water;
    let mut volume_hot_water: f64 = Default::default();
    let mut was_in_loop = false;

    while !is_close!(temperature_warm_water, temperature_target, abs_tol = 1e-10) {
        was_in_loop = true;
        // Calculate the volume of hot/cold water needed if heating from cold water source
        if temperature_target > temperature_hot_water {
            return Ok(None);
        }
        volume_hot_water = volume_warm_water
            * calc_fraction_hot_water(
                temperature_target,
                temperature_hot_water,
                temperature_cold_water,
            )?;
        let volume_cold_water = volume_warm_water - volume_hot_water;
        // Calculate the temperature of hot and warm (mixed) water given the volumes calculated
        temperature_hot_water = func_temperature_hot_water(volume_hot_water)?;
        let _ = std::mem::replace(
            &mut list_temperature_volume,
            func_temperature_cold_water(volume_cold_water, simtime)?,
        );
        temperature_cold_water = FSum::with_all(list_temperature_volume.iter().map(|(t, v)| t * v))
            .value()
            / FSum::with_all(list_temperature_volume.iter().map(|(_, v)| v)).value();
        temperature_warm_water = (volume_hot_water * temperature_hot_water
            + volume_cold_water * temperature_cold_water)
            / volume_warm_water;
    }

    Ok(if was_in_loop {
        Some(volume_hot_water)
    } else {
        None
    })
}

/// Calculate volume-weighted average temperature from list of (temperature, volume) pairs.
///
/// Arguments:
/// * `temp_volume_pairs`  -- List of (temperature, volume) tuples
/// * `expected_volume`    -- Expected total volume. If provided, validates that actual total matches.
/// * `tolerance`     -- Tolerance for volume validation (absolute tolerance)
///
/// Returns:
///     Volume-weighted average temperature
pub(crate) fn calculate_volume_weighted_average_temperature(
    temp_volume_pairs: Vec<(f64, f64)>,
    expected_volume: Option<f64>,
    tolerance: Option<f64>,
) -> anyhow::Result<f64> {
    if temp_volume_pairs.is_empty() {
        bail!("Cannot calculate weighted average: temp_volume_pairs is empty")
    }

    let tolerance = tolerance.unwrap_or(1e-10);
    let mut temp_volume_products = vec![];
    let mut volumes = vec![];

    for (temp, volume) in temp_volume_pairs {
        temp_volume_products.push(temp * volume);
        volumes.push(volume);
    }

    // Equivalent of using Python's math.fsum instead of sum() for better numerical accuracy with floating point arithmetic
    let weighted_temp_sum = FSum::with_all(&temp_volume_products).value();
    let total_volume = FSum::with_all(&volumes).value();

    if total_volume == 0. {
        bail!("Cannot calculate weighted average: total volume is zero");
    }

    // Validate expected volume if provided
    if let Some(expected_volume) = expected_volume {
        if !is_close!(total_volume, expected_volume, rel_tol = tolerance) {
            bail!("Volume mismatch: expected {expected_volume}, got {total_volume}");
        }
    }

    Ok(weighted_temp_sum / total_volume)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use pretty_assertions::assert_eq;

    #[test]
    fn test_frac_hot_water() {
        assert_eq!(
            calc_fraction_hot_water(40.0, 55.0, 5.0).unwrap(),
            0.7,
            "incorrect fraction of hot water returned"
        );

        assert!(calc_fraction_hot_water(60.0, 55.0, 5.0).is_err());
        assert!(calc_fraction_hot_water(0.0, 55.0, 5.0).is_err());
    }

    #[test]
    fn test_water_demand_to_kwh() {
        let litres_demand = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0];
        let demand_temp = [40.0, 35.0, 37.0, 39.0, 40.0, 38.0, 39.0, 40.0];
        let cold_temp = [5.0, 4.0, 5.0, 6.0, 5.0, 4.0, 3.0, 4.0];
        for i in 0..8 {
            assert_relative_eq!(
                water_demand_to_kwh(litres_demand[i], demand_temp[i], cold_temp[i]),
                [0.20339, 0.36029, 0.55787, 0.76707, 1.01694, 1.18547, 1.46440, 1.6736][i],
                max_relative = 1e-5
            );
        }
    }

    #[test]
    fn test_valid_calculation() {
        // Test normal volume-weighted average calculation
        let temp_vol_pairs = vec![(10., 5.), (20., 3.), (30., 2.)];
        let expected_temp = 17.;

        let result =
            calculate_volume_weighted_average_temperature(temp_vol_pairs, None, None).unwrap();

        assert!(is_close!(result, expected_temp, rel_tol = 1e-10));
    }

    #[test]
    fn test_volume_validation_success() {
        // Test that volume validation passes when volumes match
        let temp_vol_pairs = vec![(15., 3.), (25., 7.)];
        let expected_volume = 10.;

        let result = calculate_volume_weighted_average_temperature(
            temp_vol_pairs,
            Some(expected_volume),
            None,
        )
        .unwrap();

        assert!(is_close!(result, 22., rel_tol = 1e-10));
    }

    #[test]
    #[should_panic = "Volume mismatch: expected 5, got 10"]
    fn test_volume_validation_failure() {
        // Test that volume validation fails when volumes don't match
        let temp_vol_pairs = vec![(15., 3.), (25., 7.)]; // Total = 10.0
        let expected_volume = 5.;

        calculate_volume_weighted_average_temperature(temp_vol_pairs, Some(expected_volume), None)
            .unwrap();
    }
    #[test]
    #[should_panic = "temp_volume_pairs is empty"]
    fn test_empty_list_error() {
        // Test that empty list raises appropriate error
        calculate_volume_weighted_average_temperature(vec![], None, None).unwrap();
    }

    #[test]
    #[should_panic = "total volume is zero"]
    fn test_zero_volume_error() {
        let temp_vol_pairs = vec![(15., 0.), (25., 0.)];

        calculate_volume_weighted_average_temperature(temp_vol_pairs, None, None).unwrap();
    }
}
