#![allow(clippy::too_many_arguments)]

mod compare_floats;
pub mod core;
pub mod corpus;
mod external_conditions;
mod input;
pub mod read_weather_file;
mod simulation_time;

#[macro_use]
extern crate is_close;
extern crate lazy_static;

use crate::input::{parse_input, ExternalConditionsInput};
use crate::read_weather_file::ExternalConditions as ExternalConditionsFromFile;
use std::collections::HashMap;

pub use crate::corpus::RunResults;
use crate::corpus::{Corpus, KeyString};
use crate::external_conditions::{DaylightSavingsConfig, ExternalConditions};
use crate::simulation_time::SimulationTime;
use itertools::Itertools;
use lazy_static::lazy_static;
use std::error::Error;
use std::ffi::OsStr;
use std::io::Read;
use std::ops::Deref;
use std::path::Path;
use std::sync::Arc;

pub fn run_project<'a, T>(
    input: T,
    external_conditions_data: Option<ExternalConditionsFromFile>,
    _preprocess_only: bool,
    _fhs_assumptions: bool,
    _fhs_fee_assumptions: bool,
    _fhs_not_a_assumptions: bool,
    _fhs_not_b_assumptions: bool,
    _heat_balance: bool,
) -> Result<RunResults, anyhow::Error>
where
    T: Read,
{
    let project_data = parse_input(input);

    // println!("{:?}", project_data);

    let input = project_data.unwrap();

    let external_conditions = external_conditions_from_input(
        input.external_conditions.clone(),
        external_conditions_data,
        input.simulation_time.clone(),
    );

    let mut corpus: Corpus = Corpus::from_inputs(input, external_conditions)?;

    Ok(corpus.run().clone())
}

fn external_conditions_from_input(
    input: Arc<ExternalConditionsInput>,
    external_conditions_data: Option<ExternalConditionsFromFile>,
    simulation_time: Arc<SimulationTime>,
) -> ExternalConditions {
    match external_conditions_data {
        Some(ec) => ExternalConditions::new(
            &simulation_time.iter(),
            ec.air_temperatures,
            ec.wind_speeds,
            ec.diffuse_horizontal_radiation,
            ec.direct_beam_radiation,
            ec.solar_reflectivity_of_ground,
            ec.latitude,
            ec.longitude,
            0,
            0,
            Some(365),
            1.0,
            None,
            DaylightSavingsConfig::NotApplicable,
            false,
            ec.direct_beam_conversion_needed,
            input.shading_segments.clone(),
        ),
        None => ExternalConditions::new(
            &simulation_time.iter(),
            input.air_temperatures.clone().unwrap_or_default(),
            input.wind_speeds.clone().unwrap_or_default(),
            input
                .diffuse_horizontal_radiation
                .clone()
                .unwrap_or_default(),
            input.direct_beam_radiation.clone().unwrap_or_default(),
            input
                .solar_reflectivity_of_ground
                .clone()
                .unwrap_or_default(),
            input.latitude.unwrap_or(55.0),
            input.longitude.unwrap_or(0.0),
            0,
            0,
            Some(365),
            1.0,
            None,
            DaylightSavingsConfig::NotApplicable,
            false,
            input.direct_beam_conversion_needed.unwrap_or(false),
            input.shading_segments.clone(), //imperfect but this should be quite small...
        ),
    }
}

lazy_static! {
    pub static ref UNITS_MAP: HashMap<KeyString, &'static str> = HashMap::from([
        ("Internal gains".try_into().unwrap(), "[W]"),
        ("Solar gains".try_into().unwrap(), "[W]"),
        ("Operative temp".try_into().unwrap(), "[deg C]"),
        ("Internal air temp".try_into().unwrap(), "[deg C]"),
        ("Space heat demand".try_into().unwrap(), "[kWh]"),
        ("Space cool demand".try_into().unwrap(), "[kWh]"),
        ("Hot water demand".try_into().unwrap(), "[litres]"),
        ("Hot water energy demand".try_into().unwrap(), "[kWh]"),
        ("Hot water duration".try_into().unwrap(), "[mins]"),
        ("Hot Water Events".try_into().unwrap(), "[count]"),
        ("Pipework losses".try_into().unwrap(), "[kWh]")
    ]);
}
