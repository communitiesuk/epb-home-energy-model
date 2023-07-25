pub mod core;
mod external_conditions;
mod input;
mod project;
pub mod read_weather_file;
mod simulation_time;

#[macro_use]
extern crate is_close;

use crate::input::parse_input_file;
use crate::read_weather_file::ExternalConditions;


use std::error::Error;
use std::ffi::OsStr;
use std::path::Path;

pub fn run_project(
    input_file: &str,
    _external_conditions_data: Option<ExternalConditions>,
    _preprocess_only: bool,
    _fhs_assumptions: bool,
    _fhs_fee_assumptions: bool,
    _heat_balance: bool,
) -> Result<(), Box<dyn Error>> {
    let input_file_ext = Path::new(input_file).extension().and_then(OsStr::to_str);
    let input_file_stem = match input_file_ext {
        Some(ext) => &input_file[..(input_file.len() - ext.len() - 1)],
        None => input_file,
    };
    let _output_file = format_args!("{input_file_stem}_results.csv");
    let _output_file_static = format_args!("{input_file_stem}_results_static.csv");
    let _output_file_summary = format_args!("{input_file_stem}_results_summary.csv");

    println!("about to try and open {}", input_file);

    let project_data = parse_input_file(Path::new(input_file));

    println!("{:?}", project_data);

    Ok(())
}
