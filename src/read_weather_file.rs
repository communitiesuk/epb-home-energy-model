use csv::ReaderBuilder as CsvReaderBuilder;
use std::io::Read;

const COLUMN_LONGITUDE: usize = 7;
const COLUMN_LATITUDE: usize = 6;
const COLUMN_AIR_TEMP: usize = 6; // dry bulb temp in degrees
const COLUMN_WIND_SPEED: usize = 21; // wind speed in m/sec
const COLUMN_WIND_DIRECTION: usize = 20; // wind direction in degrees
const COLUMN_DNI_RAD: usize = 14; // direct beam normal irradiation in Wh/m2
const COLUMN_DIF_RAD: usize = 15; // diffuse irradiation (horizontal plane) in Wh/m2
const _COLUMN_GROUND_REFLECT: usize = 32;

#[derive(Clone, Debug)]
pub struct ExternalConditions {
    pub air_temperatures: Vec<f64>,
    pub wind_speeds: Vec<f64>,
    pub wind_directions: Vec<f64>,
    pub diffuse_horizontal_radiation: Vec<f64>,
    pub direct_beam_radiation: Vec<f64>,
    pub solar_reflectivity_of_ground: Vec<f64>,
    pub longitude: f64,
    pub latitude: f64,
    pub direct_beam_conversion_needed: bool,
}

pub fn weather_data_to_vec(file: impl Read) -> Result<ExternalConditions, &'static str> {
    let mut reader = CsvReaderBuilder::new()
        .flexible(true)
        .has_headers(false)
        .from_reader(file);

    let mut air_temperatures = vec![];
    let mut wind_speeds = vec![];
    let mut wind_directions = vec![];
    let mut diff_hor_rad = vec![];
    let mut dir_beam_rad = vec![];
    let mut ground_solar_reflc = vec![];
    let mut latitude: Option<f64> = None;
    let mut longitude: Option<f64> = None;

    for (i, result) in reader.records().enumerate() {
        let record: csv::StringRecord = result.unwrap();
        if i == 0 {
            latitude = Some(record.get(COLUMN_LATITUDE).unwrap().parse().unwrap());
            longitude = Some(record.get(COLUMN_LONGITUDE).unwrap().parse().unwrap());
        } else if i >= 8 {
            air_temperatures.push(record.get(COLUMN_AIR_TEMP).unwrap().parse().unwrap());
            wind_speeds.push(record.get(COLUMN_WIND_SPEED).unwrap().parse().unwrap());
            wind_directions.push(record.get(COLUMN_WIND_DIRECTION).unwrap().parse().unwrap());
            dir_beam_rad.push(record.get(COLUMN_DNI_RAD).unwrap().parse().unwrap());
            diff_hor_rad.push(record.get(COLUMN_DIF_RAD).unwrap().parse().unwrap());
            ground_solar_reflc.push(0.2);
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
