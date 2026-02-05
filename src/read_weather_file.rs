use crate::core::units::{Orientation360, KNOTS_PER_METRES_PER_SECOND};
use csv::ReaderBuilder as CsvReaderBuilder;
use std::io::Read;

const EPW_COLUMN_LONGITUDE: usize = 7;
const EPW_COLUMN_LATITUDE: usize = 6;
const EPW_COLUMN_AIR_TEMP: usize = 6; // dry bulb temp in degrees
const EPW_COLUMN_WIND_SPEED: usize = 21; // wind speed in m/sec
const EPW_COLUMN_WIND_DIRECTION: usize = 20; // wind direction in degrees
const EPW_COLUMN_DNI_RAD: usize = 14; // direct beam normal irradiation in Wh/m2
const EPW_COLUMN_DIF_RAD: usize = 15; // diffuse irradiation (horizontal plane) in Wh/m2

const SOLAR_REFLECTIVITY_OF_GROUND: f64 = 0.2;

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
    pub direct_beam_conversion_needed: bool,
}

const LIKELY_STEP_COUNT: usize = 8760; // hours in non-leap year

pub fn epw_weather_data_to_external_conditions(
    file: impl Read,
) -> Result<ExternalConditions, &'static str> {
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

    for (i, result) in reader.records().enumerate() {
        let record: csv::StringRecord = result.unwrap();
        if i == 0 {
            latitude.replace(record.get(EPW_COLUMN_LATITUDE).unwrap().parse().unwrap());
            longitude.replace(record.get(EPW_COLUMN_LONGITUDE).unwrap().parse().unwrap());
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

    Ok(ExternalConditions {
        air_temperatures,
        wind_speeds,
        wind_directions,
        diffuse_horizontal_radiation: diff_hor_rad,
        direct_beam_radiation: dir_beam_rad,
        solar_reflectivity_of_ground: ground_solar_reflc,
        latitude: latitude.unwrap(),
        longitude: longitude.unwrap(),
        direct_beam_conversion_needed: false,
    })
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
) -> Result<ExternalConditions, &'static str> {
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

    Ok(ExternalConditions {
        air_temperatures,
        wind_speeds,
        wind_directions,
        diffuse_horizontal_radiation,
        direct_beam_radiation,
        solar_reflectivity_of_ground: ground_solar_reflc,
        latitude: latitude.unwrap(),
        longitude: longitude.unwrap(),
        // Conversion is not needed as direct irradiation will be normal plane from this file
        direct_beam_conversion_needed: false,
    })
}

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
