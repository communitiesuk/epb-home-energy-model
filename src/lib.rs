#![allow(clippy::too_many_arguments)]

mod compare_floats;
pub mod core;
pub mod corpus;
mod external_conditions;
pub mod input;
pub mod output;
pub mod read_weather_file;
mod simulation_time;
mod statistics;
mod wrappers;

#[macro_use]
extern crate is_close;
extern crate lazy_static;

pub use crate::corpus::RunResults;
use crate::corpus::{Corpus, HotWaterResultMap, KeyString};
use crate::external_conditions::{DaylightSavingsConfig, ExternalConditions};
use crate::input::{ingest_for_processing, ExternalConditionsInput};
use crate::output::Output;
use crate::read_weather_file::ExternalConditions as ExternalConditionsFromFile;
use crate::simulation_time::SimulationTime;
use crate::wrappers::future_homes_standard::future_homes_standard::apply_fhs_preprocessing;
use csv::WriterBuilder;
use indexmap::IndexMap;
use lazy_static::lazy_static;
use std::borrow::Cow;
use std::io::Read;
use std::sync::Arc;

pub fn run_project(
    input: impl Read,
    output: impl Output,
    external_conditions_data: Option<ExternalConditionsFromFile>,
    _preprocess_only: bool,
    fhs_assumptions: bool,
    _fhs_fee_assumptions: bool,
    fhs_not_a_assumptions: bool,
    fhs_not_b_assumptions: bool,
    _heat_balance: bool,
) -> Result<(), anyhow::Error> {
    let mut input_for_processing = ingest_for_processing(input)?;

    // do wrapper pre-processing here
    if fhs_assumptions || fhs_not_a_assumptions || fhs_not_b_assumptions {
        apply_fhs_preprocessing(&mut input_for_processing)?;
    }

    let input = input_for_processing.finalize();

    let external_conditions = external_conditions_from_input(
        input.external_conditions.clone(),
        external_conditions_data,
        input.simulation_time,
    );

    let mut corpus: Corpus = Corpus::from_inputs(input, Some(external_conditions))?;

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
        _heat_cop_dict,
        _cool_cop_dict,
        _dhw_cop_dict,
        ductwork_gains,
        _heat_balance_dict,
        _heat_source_wet_results_dict,
        _heat_source_wet_results_annual_dict,
    ) = corpus.run();

    write_core_output_file(
        output,
        OutputFileArgs {
            output_key: "results".to_string(),
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
        },
    )?;

    Ok(())
}

fn external_conditions_from_input(
    input: Arc<ExternalConditionsInput>,
    external_conditions_data: Option<ExternalConditionsFromFile>,
    simulation_time: SimulationTime,
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
    pub static ref UNITS_MAP: IndexMap<KeyString, &'static str> = IndexMap::from([
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

struct OutputFileArgs {
    output_key: String,
    timestep_array: Vec<f64>,
    results_totals: IndexMap<KeyString, Vec<f64>>,
    results_end_user: IndexMap<KeyString, IndexMap<String, Vec<f64>>>,
    energy_import: IndexMap<KeyString, Vec<f64>>,
    energy_export: IndexMap<KeyString, Vec<f64>>,
    energy_generated_consumed: IndexMap<KeyString, Vec<f64>>,
    energy_to_storage: IndexMap<KeyString, Vec<f64>>,
    energy_from_storage: IndexMap<KeyString, Vec<f64>>,
    energy_diverted: IndexMap<KeyString, Vec<f64>>,
    betafactor: IndexMap<KeyString, Vec<f64>>,
    zone_dict: IndexMap<KeyString, IndexMap<KeyString, Vec<f64>>>,
    zone_list: Vec<KeyString>,
    hc_system_dict: IndexMap<KeyString, IndexMap<KeyString, Vec<f64>>>,
    hot_water_dict: IndexMap<KeyString, HotWaterResultMap>,
    ductwork_gains: IndexMap<KeyString, Vec<f64>>,
}

fn write_core_output_file(output: impl Output, args: OutputFileArgs) -> Result<(), anyhow::Error> {
    let OutputFileArgs {
        output_key,
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
    } = args;
    println!("writing out to {output_key}");
    let writer = output.writer_for_location_key(&output_key)?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    let mut headings: Vec<Cow<'static, str>> = vec!["Timestep".into()];
    let mut units_row = vec!["[count]"];

    for totals_key in results_totals.keys() {
        let totals_header = format!("{totals_key} total");
        headings.push(totals_header.into());
        units_row.push("[kWh]");
        for end_user_key in results_end_user[totals_key].keys() {
            headings.push((*end_user_key).clone().into());
            units_row.push("[kWh]");
        }
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
    }

    for zone in &zone_list {
        for zone_outputs in zone_dict.keys() {
            let zone_headings = format!("{zone_outputs} {zone}");
            headings.push(zone_headings.into());
            if UNITS_MAP.contains_key(zone_outputs) {
                units_row.push(UNITS_MAP[zone_outputs]);
            } else {
                units_row.push("Unit not defined");
            }
        }
    }

    for system in hc_system_dict.keys() {
        if hc_system_dict[system].keys().len() == 0 {
            // edge case - if we don't have any hc systems keys use None
            let none_heading = format!("{system} None");
            headings.push(none_heading.into());
            units_row.push("[kWh]");
            continue;
        }

        for hc_name in hc_system_dict[system].keys() {
            let hc_system_headings = format!("{system} {hc_name}");
            headings.push(hc_system_headings.into());
            units_row.push("[kWh]");
        }
    }

    for system in hot_water_dict.keys() {
        headings.push(system.to_string().into());

        if UNITS_MAP.contains_key(system) {
            units_row.push(UNITS_MAP[system]);
        } else {
            units_row.push("Unit not defined");
        }
    }

    headings.push("Ductwork gains".into());
    units_row.push("[kWh]");

    // Write headings and units to output file
    writer.write_record(headings.iter().map(|heading| heading.as_ref()))?;
    writer.write_record(&units_row)?;

    for (t_idx, _timestep) in timestep_array.iter().enumerate() {
        let mut energy_use_row = vec![];
        let mut zone_row = vec![];
        let mut hc_system_row = vec![];
        let mut hw_system_row = vec![];
        let mut hw_system_row_energy = vec![];
        let mut hw_system_row_duration = vec![];
        let mut hw_system_row_events = vec![];
        let mut pw_losses_row = vec![];
        let mut ductwork_row = vec![];
        let energy_shortfall: Vec<f64> = vec![];
        for totals_key in results_totals.keys() {
            energy_use_row.push(results_totals[totals_key][t_idx]);
            for (end_user_key, _) in results_end_user[totals_key].iter().enumerate() {
                energy_use_row.push(results_end_user[totals_key][end_user_key][t_idx]);
            }
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
                zone_row.push(zone_dict[zone_outputs][zone][t_idx]);
            }
        }

        // Loop over heating and cooling system demand
        for system in hc_system_dict.keys() {
            // if we don't have any data for this system
            // we still need to output 0s for the "{system} None" column
            // Expect this to move into the calculation logic in future
            if hc_system_dict[system].keys().len() == 0 {
                hc_system_row.push(0.);
            }

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
        let mut row: Vec<String> = vec![];
        row.append(&mut vec![t_idx.to_string()]);
        row.append(
            &mut energy_use_row
                .into_iter()
                .map(|val| val.to_string())
                .collect(),
        );
        row.append(&mut zone_row.into_iter().map(|val| val.to_string()).collect());
        row.append(
            &mut hc_system_row
                .into_iter()
                .map(|val| val.to_string())
                .collect(),
        );
        row.append(
            &mut hw_system_row
                .into_iter()
                .map(|val| val.to_string())
                .collect(),
        );
        row.append(
            &mut hw_system_row_energy
                .into_iter()
                .map(|val| val.to_string())
                .collect(),
        );
        row.append(
            &mut hw_system_row_duration
                .into_iter()
                .map(|val| val.to_string())
                .collect(),
        );
        row.append(
            &mut hw_system_row_events
                .into_iter()
                .map(|val| val.to_string())
                .collect(),
        );
        row.append(
            &mut pw_losses_row
                .into_iter()
                .map(|val| val.to_string())
                .collect(),
        );
        row.append(
            &mut ductwork_row
                .into_iter()
                .map(|val| val.to_string())
                .collect(),
        );
        row.append(
            &mut energy_shortfall
                .into_iter()
                .map(|val| val.to_string())
                .collect(),
        );

        writer.write_record(&row)?;
    }

    println!("flushing out CSV");
    writer.flush()?;

    Ok(())
}
