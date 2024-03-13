#![allow(clippy::too_many_arguments)]

mod compare_floats;
pub mod core;
mod corpus;
mod external_conditions;
mod input;
pub mod read_weather_file;
mod simulation_time;

#[macro_use]
extern crate is_close;
extern crate lazy_static;

use crate::input::{parse_input_file, ExternalConditionsInput, Input, Ventilation};
use crate::read_weather_file::ExternalConditions as ExternalConditionsFromFile;
use std::borrow::Cow;
use std::collections::HashMap;

use crate::core::units::{SECONDS_PER_HOUR, WATTS_PER_KILOWATT};
use crate::corpus::{Corpus, HotWaterResultMap};
use crate::external_conditions::{DaylightSavingsConfig, ExternalConditions};
use crate::simulation_time::{SimulationTime, SimulationTimeIteration, SimulationTimeIterator};
use csv::Writer;
use itertools::Itertools;
use lazy_static::lazy_static;
use std::error::Error;
use std::ffi::OsStr;
use std::io::Read;
use std::ops::Deref;
use std::path::Path;
use std::sync::Arc;

pub fn run_project(
    input_file: &str,
    external_conditions_data: Option<ExternalConditionsFromFile>,
    _preprocess_only: bool,
    _fhs_assumptions: bool,
    _fhs_fee_assumptions: bool,
    _fhs_not_a_assumptions: bool,
    _fhs_not_b_assumptions: bool,
    _heat_balance: bool,
) -> Result<(), Box<dyn Error>> {
    let input_file_ext = Path::new(input_file).extension().and_then(OsStr::to_str);
    let input_file_stem = match input_file_ext {
        Some(ext) => &input_file[..(input_file.len() - ext.len() - 1)],
        None => input_file,
    };
    let output_file_detailed = format!("{input_file_stem}_results.csv");
    let _output_file_static = format!("{input_file_stem}_results_static.csv");
    let _output_file_summary = format!("{input_file_stem}_results_summary.csv");

    println!("about to try and open {}", input_file);

    let project_data = parse_input_file(Path::new(input_file));

    // println!("{:?}", project_data);

    let input = project_data.unwrap();

    let external_conditions = external_conditions_from_input(
        input.external_conditions.clone(),
        external_conditions_data,
        input.simulation_time.clone(),
    );

    let mut corpus: Corpus = Corpus::from_inputs(input, external_conditions)?;

    let (
        timestep_array,
        results_totals,
        results_end_user,
        energy_import,
        energy_export,
        energy_generated_consumed,
        energy_to_storage,
        energy_from_storage,
        energy_diverted,
        betafactor,
        zone_dict,
        zone_list,
        hc_system_dict,
        hot_water_dict,
        heat_cop_dict,
        cool_cop_dict,
        dhw_cop_dict,
        ductwork_gains,
        heat_balance_dict,
        heat_source_wet_results_dict,
        heat_source_wet_results_annual_dict,
    ) = corpus.run();

    let _ = write_core_output_file(
        output_file_detailed,
        timestep_array,
        results_totals,
        results_end_user,
        energy_import,
        energy_export,
        energy_generated_consumed,
        energy_to_storage,
        energy_from_storage,
        energy_diverted,
        betafactor,
        zone_dict,
        zone_list,
        hc_system_dict,
        hot_water_dict,
        ductwork_gains,
    );

    Ok(())
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

fn write_core_output_file(
    output_file: String,
    timestep_array: Vec<f64>,
    results_totals: HashMap<&str, Vec<f64>>,
    _results_end_user: HashMap<&str, Vec<f64>>,
    energy_import: HashMap<&str, Vec<f64>>,
    energy_export: HashMap<&str, Vec<f64>>,
    energy_generated_consumed: HashMap<&str, Vec<f64>>,
    energy_to_storage: HashMap<&str, Vec<f64>>,
    energy_from_storage: HashMap<&str, Vec<f64>>,
    energy_diverted: HashMap<&str, Vec<f64>>,
    betafactor: HashMap<&str, Vec<f64>>,
    zone_dict: HashMap<&str, HashMap<String, Vec<f64>>>,
    zone_list: Vec<&str>,
    hc_system_dict: HashMap<&str, HashMap<String, Vec<f64>>>,
    hot_water_dict: HashMap<&str, HotWaterResultMap>,
    ductwork_gains: HashMap<&str, Vec<f64>>,
) -> Result<(), anyhow::Error> {
    let mut writer = Writer::from_path(output_file)?;

    let mut headings: Vec<Cow<'static, str>> = vec!["Timestep".into()];
    let mut units_row = vec!["[count]"];
    for totals_key in results_totals.keys() {
        let totals_header = format!("{totals_key} total");
        headings.push(totals_header.into());
        units_row.push("[kWh]");
        // TODO enable when energy supplies implemented
        // for end_user_key in results_end_user[totals_key].keys() {
        //     headings.push(end_user_key.into());
        //     units_row.push("[kWh]");
        // }
        headings.push(format!("{totals_key} import").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key} export").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key} generated and consumed").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key} beta factor").into());
        units_row.push("[ratio]");
        headings.push(format!("{totals_key} to storage").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key} from storage").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key} diverted").into());
        units_row.push("[kWh]");

        for zone in &zone_list {
            for zone_outputs in zone_dict.keys() {
                let zone_headings = format!("{zone_outputs} {zone}");
                headings.push(zone_headings.into());
                if UNITS_MAP.keys().contains(zone_outputs) {
                    units_row.push(UNITS_MAP[zone_outputs]);
                } else {
                    units_row.push("Unit not defined");
                }
            }
        }

        for system in hc_system_dict.keys() {
            for hc_name in hc_system_dict[system].keys() {
                let hc_system_headings = format!("{system} {hc_name}");
                headings.push(hc_system_headings.into());
                units_row.push("[kWh]");
            }
        }

        headings.push("Ductwork gains".into());
        units_row.push("[kWh]");

        // Write headings and units to output file
        writer.write_record(headings.iter().map(|heading| heading.as_ref()))?;
        writer.write_record(&units_row)?;

        for (t_idx, timestep) in timestep_array.iter().enumerate() {
            let mut energy_use_row = vec![];
            let mut zone_row = vec![];
            let mut hc_system_row = vec![];
            let mut hw_system_row = vec![];
            let mut hw_system_row_energy = vec![];
            let mut hw_system_row_duration = vec![];
            let mut hw_system_row_events = vec![];
            let mut pw_losses_row = vec![];
            let mut ductwork_row = vec![];
            let mut _energy_shortfall: Vec<f64> = vec![];
            for totals_key in results_totals.keys() {
                energy_use_row.push(results_totals[totals_key][t_idx]);
                // TODO uncomment when energy supplies implemented
                // for end_user_key in results_end_user[totals_key] {
                //     energy_use_row.push(results_end_user[totals_key][end_user_key][t_idx]);
                // }
                energy_use_row.push(energy_import[totals_key][t_idx]);
                energy_use_row.push(energy_export[totals_key][t_idx]);
                energy_use_row.push(energy_generated_consumed[totals_key][t_idx]);
                energy_use_row.push(betafactor[totals_key][t_idx]);
                energy_use_row.push(energy_to_storage[totals_key][t_idx]);
                energy_use_row.push(energy_from_storage[totals_key][t_idx]);
                energy_use_row.push(energy_diverted[totals_key][t_idx]);
            }

            // Loop over results separated by zone
            for zone in &zone_list {
                for zone_outputs in zone_dict.keys() {
                    zone_row.push(zone_dict[zone_outputs][zone.to_owned()][t_idx]);
                }
            }

            // Loop over heading and cooling system demand
            for system in hc_system_dict.keys() {
                for hc_name in hc_system_dict[system].keys() {
                    hc_system_row.push(hc_system_dict[system][hc_name][t_idx]);
                }
            }

            // Loop over hot water demand
            if let HotWaterResultMap::Float(map) = &hot_water_dict["Hot water demand"] {
                hw_system_row.push(map["demand"][t_idx]);
            }
            if let HotWaterResultMap::Float(map) = &hot_water_dict["Hot water energy demand"] {
                hw_system_row_energy.push(map["energy_demand"][t_idx]);
            }
            if let HotWaterResultMap::Float(map) = &hot_water_dict["Hot water duration"] {
                hw_system_row_duration.push(map["duration"][t_idx]);
            }
            if let HotWaterResultMap::Float(map) = &hot_water_dict["Pipework losses"] {
                pw_losses_row.push(map["pw_losses"][t_idx]);
            }
            if let HotWaterResultMap::Int(map) = &hot_water_dict["Hot Water Events"] {
                hw_system_row_events.push(map["no_events"][t_idx]);
            }
            ductwork_row.push(ductwork_gains["ductwork_gains"][t_idx]);

            // create row of outputs and write to output file
            // TODO work out how to compile a row of results from the following
            // let mut row = vec![];
            // row.append(&mut vec![t_idx.to_string().as_str()]);
            // row.append(
            //     energy_use_row
            //         .into_iter()
            //         .map(|val| val.to_string().as_mut())
            //         .collect(),
            // );
            // row.append(
            //     zone_row
            //         .into_iter()
            //         .map(|val| val.to_string().as_mut())
            //         .collect(),
            // );
            // row.append(
            //     hc_system_row
            //         .into_iter()
            //         .map(|val| val.to_string().as_mut())
            //         .collect(),
            // );
            // row.append(
            //     hw_system_row
            //         .into_iter()
            //         .map(|val| val.to_string().as_mut())
            //         .collect(),
            // );
            // row.append(
            //     hw_system_row_energy
            //         .into_iter()
            //         .map(|val| val.to_string().as_mut())
            //         .collect(),
            // );
            // row.append(
            //     hw_system_row_duration
            //         .into_iter()
            //         .map(|val| val.to_string().as_mut())
            //         .collect(),
            // );
            // row.append(
            //     hw_system_row_events
            //         .into_iter()
            //         .map(|val| val.to_string().as_mut())
            //         .collect(),
            // );
            // row.append(
            //     pw_losses_row
            //         .into_iter()
            //         .map(|val| val.to_string().as_mut())
            //         .collect(),
            // );
            // row.append(
            //     ductwork_row
            //         .into_iter()
            //         .map(|val| val.to_string().as_mut())
            //         .collect(),
            // );
            // row.append(
            //     energy_shortfall
            //         .into_iter()
            //         .map(|val| val.to_string().as_mut())
            //         .collect(),
            // );
            //
            // writer.write_record(&row)?;
        }
    }

    writer.flush()?;

    Ok(())
}

lazy_static! {
    pub static ref UNITS_MAP: HashMap<&'static str, &'static str> = HashMap::from([
        ("Internal gains", "[W]"),
        ("Solar gains", "[W]"),
        ("Operative temp", "[deg C]"),
        ("Internal air temp", "[deg C]"),
        ("Space heat demand", "[kWh]"),
        ("Space cool demand", "[kWh]"),
        ("Hot water demand", "[litres]"),
        ("Hot water energy demand", "[kWh]"),
        ("Hot water duration", "[mins]"),
        ("Hot Water Events", "[count]"),
        ("Pipework losses", "[kWh]")
    ]);
}
