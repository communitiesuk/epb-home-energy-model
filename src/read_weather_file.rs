use crate::core::units::{Orientation360, KNOTS_PER_METRES_PER_SECOND};
use csv::ReaderBuilder as CsvReaderBuilder;
use std::io::Read;
use thiserror::Error;

const EPW_COLUMN_LONGITUDE: usize = 7;
const EPW_COLUMN_LATITUDE: usize = 6;
const EPW_COLUMN_TIMEZONE: usize = 8; // time zone in hours relative to GMT (LOCATION record)
const EPW_COLUMN_AIR_TEMP: usize = 6; // dry bulb temp in degrees
const EPW_COLUMN_WIND_SPEED: usize = 21; // wind speed in m/sec
const EPW_COLUMN_WIND_DIRECTION: usize = 20; // wind direction in degrees
const EPW_COLUMN_DNI_RAD: usize = 14; // direct beam normal irradiation in Wh/m2
const EPW_COLUMN_DIF_RAD: usize = 15; // diffuse irradiation (horizontal plane) in Wh/m2

const SOLAR_REFLECTIVITY_OF_GROUND: f64 = 0.2;

// supported range of time zones, in hours ahead of UTC (as in EnergyPlus)
const TIMEZONE_MIN: f64 = -12.;
const TIMEZONE_MAX: f64 = 14.;

#[derive(Clone, Debug)]
pub struct ExternalConditions {
    pub air_temperatures: Vec<f64>,
    pub wind_speeds: Vec<f64>,
    pub wind_directions: Vec<Orientation360>,
    pub diffuse_horizontal_radiation: Vec<f64>,
    pub direct_beam_radiation: Vec<f64>,
    pub solar_reflectivity_of_ground: Vec<f64>,
    pub longitude: f64,
    pub latitude: f64,
    pub timezone: f64,
    pub direct_beam_conversion_needed: bool,
}

const LIKELY_STEP_COUNT: usize = 8760; // hours in non-leap year

pub fn epw_weather_data_to_external_conditions(
    file: impl Read,
) -> ReadWeatherFileResult<ExternalConditions> {
    let mut reader = CsvReaderBuilder::new()
        .flexible(true)
        .has_headers(false)
        .from_reader(file);

    let mut air_temperatures = Vec::with_capacity(LIKELY_STEP_COUNT);
    let mut wind_speeds = Vec::with_capacity(LIKELY_STEP_COUNT);
    let mut wind_directions = Vec::with_capacity(LIKELY_STEP_COUNT);
    let mut diff_hor_rad = Vec::with_capacity(LIKELY_STEP_COUNT);
    let mut dir_beam_rad = Vec::with_capacity(LIKELY_STEP_COUNT);
    let mut ground_solar_reflc = Vec::with_capacity(LIKELY_STEP_COUNT);
    let mut latitude: Option<f64> = None;
    let mut longitude: Option<f64> = None;
    let mut timezone: Option<f64> = None;

    for (i, result) in reader.records().enumerate() {
        let record: csv::StringRecord = result.unwrap();
        if i == 0 {
            latitude.replace(record.get(EPW_COLUMN_LATITUDE).unwrap().parse().unwrap());
            longitude.replace(record.get(EPW_COLUMN_LONGITUDE).unwrap().parse().unwrap());
            timezone.replace(parse_timezone(record.get(EPW_COLUMN_TIMEZONE))?);
        } else if i >= 8 {
            air_temperatures.push(record.get(EPW_COLUMN_AIR_TEMP).unwrap().parse().unwrap());
            wind_speeds.push(record.get(EPW_COLUMN_WIND_SPEED).unwrap().parse().unwrap());
            wind_directions.push(
                record
                    .get(EPW_COLUMN_WIND_DIRECTION)
                    .unwrap()
                    .parse()
                    .unwrap(),
            );
            dir_beam_rad.push(record.get(EPW_COLUMN_DNI_RAD).unwrap().parse().unwrap());
            diff_hor_rad.push(record.get(EPW_COLUMN_DIF_RAD).unwrap().parse().unwrap());
            ground_solar_reflc.push(SOLAR_REFLECTIVITY_OF_GROUND);
        }
    }

    let external_conditions = ExternalConditions {
        air_temperatures,
        wind_speeds,
        wind_directions,
        diffuse_horizontal_radiation: diff_hor_rad,
        direct_beam_radiation: dir_beam_rad,
        solar_reflectivity_of_ground: ground_solar_reflc,
        latitude: latitude.unwrap(),
        longitude: longitude.unwrap(),
        timezone: timezone.unwrap(),
        direct_beam_conversion_needed: false,
    };

    validate_weather_data(&external_conditions)?;
    Ok(external_conditions)
}

const CIBSE_COLUMN_LONGITUDE: usize = 3;
const CIBSE_COLUMN_LATITUDE: usize = 1;
const CIBSE_COLUMN_AIR_TEMP: usize = 6; // dry bulb temp in degrees
const CIBSE_COLUMN_WIND_SPEED: usize = 11; // wind speed in knots
const CIBSE_COLUMN_WIND_DIRECTION: usize = 10; // wind direction in degrees
const CIBSE_COLUMN_GHI_RAD: usize = 12; // global irradiation (horizontal plane) in Wh/m2
const CIBSE_COLUMN_DIF_RAD: usize = 13; // diffuse irradiation (horizontal plane) in Wh/m2

pub fn cibse_weather_data_to_external_conditions(
    file: impl Read,
) -> ReadWeatherFileResult<ExternalConditions> {
    let mut reader = CsvReaderBuilder::new()
        .flexible(true)
        .has_headers(false)
        .from_reader(file);

    let mut air_temperatures = Vec::with_capacity(LIKELY_STEP_COUNT);
    let mut wind_speeds = Vec::with_capacity(LIKELY_STEP_COUNT);
    let mut wind_directions = Vec::with_capacity(LIKELY_STEP_COUNT);
    let mut diffuse_horizontal_radiation = Vec::with_capacity(LIKELY_STEP_COUNT);
    let mut direct_beam_radiation = Vec::with_capacity(LIKELY_STEP_COUNT);
    let mut ground_solar_reflc = Vec::with_capacity(LIKELY_STEP_COUNT);
    let mut latitude: Option<f64> = None;
    let mut longitude: Option<f64> = None;

    for (i, result) in reader.records().enumerate() {
        let record: csv::StringRecord = result.unwrap();
        if i == 5 {
            longitude.replace(record[CIBSE_COLUMN_LONGITUDE].parse().unwrap());
            latitude.replace(record[CIBSE_COLUMN_LATITUDE].parse().unwrap());
        } else if i >= 32 {
            air_temperatures.push(record[CIBSE_COLUMN_AIR_TEMP].parse().unwrap());
            wind_speeds.push(
                record[CIBSE_COLUMN_WIND_SPEED].parse::<f64>().unwrap()
                    / KNOTS_PER_METRES_PER_SECOND,
            );
            wind_directions.push(record[CIBSE_COLUMN_WIND_DIRECTION].parse().unwrap());
            // no DNI direct irradiation in file need to extract from global and diffuse values
            let global_horiz_irr: f64 = record[CIBSE_COLUMN_GHI_RAD].parse().unwrap();
            let diffuse_horiz_irr: f64 = record[CIBSE_COLUMN_DIF_RAD].parse().unwrap();
            direct_beam_radiation.push(global_horiz_irr - diffuse_horiz_irr);
            diffuse_horizontal_radiation.push(record[CIBSE_COLUMN_DIF_RAD].parse().unwrap());
            ground_solar_reflc.push(SOLAR_REFLECTIVITY_OF_GROUND);
        }
    }

    let external_conditions = ExternalConditions {
        air_temperatures,
        wind_speeds,
        wind_directions,
        diffuse_horizontal_radiation,
        direct_beam_radiation,
        solar_reflectivity_of_ground: ground_solar_reflc,
        latitude: latitude.unwrap(),
        longitude: longitude.unwrap(),
        // CIBSE weather files are for UK locations and use GMT
        timezone: 0.,
        // CIBSE format provides direct beam as horizontal irradiance; conversion to normal plane required
        direct_beam_conversion_needed: true,
    };

    validate_weather_data(&external_conditions)?;
    Ok(external_conditions)
}

/// Parse a time zone in hours ahead of UTC, which may be fractional (e.g. 9.5 for UTC+9:30).
fn parse_timezone(field: Option<&str>) -> ReadWeatherFileResult<f64> {
    let field = field.unwrap_or_default();
    match field.trim().parse::<f64>() {
        // the range check also rejects NaN and infinite values
        Ok(timezone) if (TIMEZONE_MIN..=TIMEZONE_MAX).contains(&timezone) => Ok(timezone),
        _ => Err(ReadWeatherFileError::InvalidTimezone(field.to_string())),
    }
}

fn validate_weather_data(external_conditions: &ExternalConditions) -> ReadWeatherFileResult<()> {
    if [
        external_conditions.air_temperatures.len(),
        external_conditions.wind_speeds.len(),
        external_conditions.wind_directions.len(),
        external_conditions.diffuse_horizontal_radiation.len(),
        external_conditions.direct_beam_radiation.len(),
        external_conditions.solar_reflectivity_of_ground.len(),
    ]
    .iter()
    .any(|&i| i != 8760)
    {
        Err(ReadWeatherFileError::InvalidLength)
    } else {
        Ok(())
    }
}

#[derive(Debug, Error, PartialEq)]
pub enum ReadWeatherFileError {
    #[error("Weather data should contain at least 8760 entries")]
    InvalidLength,
    #[error(
        "Time zone in weather file should be a number of hours ahead of UTC from {min} to {max}, but was {0:?}",
        min = TIMEZONE_MIN,
        max = TIMEZONE_MAX
    )]
    InvalidTimezone(String),
}

pub type ReadWeatherFileResult<T> = Result<T, ReadWeatherFileError>;

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

    #[fixture]
    fn cibse_weather_file() -> &'static [u8] {
        include_bytes!("../examples/weather_data/London_weather_CIBSE_format.csv")
    }

    #[fixture]
    fn epw_weather_file() -> &'static [u8] {
        include_bytes!("../examples/weather_data/London_weather_EnergyPlus_format.epw")
    }

    #[rstest]
    fn test_cibse_weather_data_to_external_conditions(cibse_weather_file: impl Read) {
        let external_conditions =
            cibse_weather_data_to_external_conditions(cibse_weather_file).unwrap();
        assert!([
            external_conditions.air_temperatures.len(),
            external_conditions.wind_speeds.len(),
            external_conditions.wind_directions.len(),
            external_conditions.diffuse_horizontal_radiation.len(),
            external_conditions.direct_beam_radiation.len(),
            external_conditions.solar_reflectivity_of_ground.len()
        ]
        .iter()
        .all(|&v| v == 8760));
    }

    /// The EPW weather file with its LOCATION record replaced by one for a site at 48.1 N, 11.6 E
    /// with the given time zone field (or none)
    fn epw_weather_file_with_timezone(
        epw_weather_file: &[u8],
        timezone_field: Option<&str>,
    ) -> String {
        let (_, data_after_location) = std::str::from_utf8(epw_weather_file)
            .unwrap()
            .split_once('\n')
            .unwrap();
        let location = "LOCATION,unknown,-,unknown,unknown,unknown,48.1,11.6";
        match timezone_field {
            Some(timezone_field) => {
                format!("{location},{timezone_field},520.0\r\n{data_after_location}")
            }
            None => format!("{location}\r\n{data_after_location}"),
        }
    }

    #[rstest]
    #[case("1.0", 1.)]
    #[case("-5", -5.)]
    #[case("9.5", 9.5)]
    #[case("5.75", 5.75)]
    #[case("-3.5", -3.5)]
    #[case("-12", -12.)]
    #[case("13", 13.)]
    #[case("14.0", 14.)]
    fn test_epw_weather_data_reads_timezone(
        epw_weather_file: &'static [u8],
        #[case] timezone_field: &str,
        #[case] expected_timezone: f64,
    ) {
        let external_conditions =
            epw_weather_data_to_external_conditions(epw_weather_file).unwrap();
        assert_eq!(external_conditions.timezone, 0.);

        let weather_file = epw_weather_file_with_timezone(epw_weather_file, Some(timezone_field));
        let external_conditions =
            epw_weather_data_to_external_conditions(weather_file.as_bytes()).unwrap();
        assert_eq!(external_conditions.timezone, expected_timezone);
        assert_eq!(external_conditions.latitude, 48.1);
        assert_eq!(external_conditions.longitude, 11.6);
    }

    #[rstest]
    #[case::missing(None)]
    #[case::empty(Some(""))]
    #[case::not_a_number(Some("UTC+1"))]
    #[case::nan(Some("NaN"))]
    #[case::infinite(Some("inf"))]
    #[case::negative_infinite(Some("-inf"))]
    #[case::above_range(Some("15"))]
    #[case::above_range_fractional(Some("14.5"))]
    #[case::below_range(Some("-13"))]
    fn test_epw_weather_data_rejects_invalid_timezone(
        epw_weather_file: &'static [u8],
        #[case] timezone_field: Option<&str>,
    ) {
        let weather_file = epw_weather_file_with_timezone(epw_weather_file, timezone_field);
        let error = epw_weather_data_to_external_conditions(weather_file.as_bytes()).unwrap_err();
        assert_eq!(
            error,
            ReadWeatherFileError::InvalidTimezone(timezone_field.unwrap_or_default().to_string())
        );
        assert!(error.to_string().contains("from -12 to 14"));
    }

    #[rstest]
    fn test_weather_data_to_external_conditions(epw_weather_file: impl Read) {
        let external_conditions =
            epw_weather_data_to_external_conditions(epw_weather_file).unwrap();
        assert!([
            external_conditions.air_temperatures.len(),
            external_conditions.wind_speeds.len(),
            external_conditions.wind_directions.len(),
            external_conditions.diffuse_horizontal_radiation.len(),
            external_conditions.direct_beam_radiation.len(),
            external_conditions.solar_reflectivity_of_ground.len()
        ]
        .iter()
        .all(|&v| v == 8760));
    }
}

// pub air_temperatures: Vec<f64>,
//     pub wind_speeds: Vec<f64>,
//     pub wind_directions: Vec<f64>,
//     pub diffuse_horizontal_radiation: Vec<f64>,
//     pub direct_beam_radiation: Vec<f64>,
//     pub solar_reflectivity_of_ground: Vec<f64>,
