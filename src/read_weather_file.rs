use csv::ReaderBuilder as CsvReaderBuilder;
use std::fs::File;

const COLUMN_LONGITUDE: usize = 7;
const COLUMN_LATITUDE: usize = 6;
const COLUMN_AIR_TEMP: usize = 6; // dry bulb temp in degrees
const COLUMN_WIND_SPEED: usize = 21; // wind speed in m/sec
const COLUMN_DNI_RAD: usize = 14; // direct beam normal irradiation in Wh/m2
const COLUMN_DIF_RAD: usize = 15; // diffuse irradiation (horizontal plane) in Wh/m2
const COLUMN_GROUND_REFLECT: usize = 32;

#[derive(Debug)]
pub struct ExternalConditions {
    air_temperatures: Vec<f64>,
    wind_speeds: Vec<f64>,
    diffuse_horizontal_radiation: Vec<f64>,
    direct_beam_radiation: Vec<f64>,
    solar_reflectivity_of_ground: Vec<f64>,
    longitude: f64,
    latitude: f64,
    direct_beam_conversion_needed: bool,
}

struct ExternalConditionDatum {
    air_temperature: f64,
    wind_speed: f64,
    diffuse_horizontal_radiation: f64,
    direct_beam_radiation: f64,
    solar_reflectivity_of_ground: f64,
}

pub fn weather_data_to_vec(file: &str) -> Result<ExternalConditions, &'static str> {
    let file = match File::open(file) {
        Ok(f) => f,
        Err(_) => {
            return Err("The weather file provided did not exist!");
        }
    };
    let mut reader = CsvReaderBuilder::new()
        .flexible(true)
        .has_headers(false)
        .from_reader(file);

    let mut air_temperatures = vec![];
    let mut wind_speeds = vec![];
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
            dir_beam_rad.push(record.get(COLUMN_DIF_RAD).unwrap().parse().unwrap());
            diff_hor_rad.push(record.get(COLUMN_DNI_RAD).unwrap().parse().unwrap());
            ground_solar_reflc.push(0.2); // this could be an upstream bug as is not reading from the file?
        }
    }

    Ok(ExternalConditions {
        air_temperatures,
        wind_speeds,
        diffuse_horizontal_radiation: diff_hor_rad,
        direct_beam_radiation: dir_beam_rad,
        solar_reflectivity_of_ground: ground_solar_reflc,
        latitude: latitude.unwrap(),
        longitude: longitude.unwrap(),
        direct_beam_conversion_needed: false,
    })
}
