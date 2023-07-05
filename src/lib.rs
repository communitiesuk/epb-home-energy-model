mod external_conditions;
mod input;
pub mod read_weather_file;
mod simulation_time;

use crate::input::parse_input_file;
use crate::read_weather_file::ExternalConditions;
use serde_json;
use serde_json::Value;
use std::error::Error;
use std::ffi::OsStr;
use std::path::Path;

pub fn run_project(
    input_file: &str,
    external_conditions_data: Option<ExternalConditions>,
    preprocess_only: bool,
    fhs_assumptions: bool,
    fhs_fee_assumptions: bool,
    heat_balance: bool,
) -> Result<(), Box<dyn Error>> {
    let input_file_ext = Path::new(input_file).extension().and_then(OsStr::to_str);
    let input_file_stem = match input_file_ext {
        Some(ext) => &input_file[..(input_file.len() - ext.len() - 1)],
        None => input_file,
    };
    let output_file = format_args!("{input_file_stem}_results.csv");
    let output_file_static = format_args!("{input_file_stem}_results_static.csv");
    let output_file_summary = format_args!("{input_file_stem}_results_summary.csv");

    println!("about to try and open {}", input_file);

    let project_data = parse_input_file(Path::new(input_file));

    println!("{:?}", project_data);

    Ok(())
}
