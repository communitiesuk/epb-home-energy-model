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

use crate::core::units::convert_profile_to_daily;
pub use crate::corpus::RunResults;
use crate::corpus::{Corpus, HotWaterResultMap, KeyString, NumberOrDivisionByZero};
use crate::external_conditions::{DaylightSavingsConfig, ExternalConditions};
use crate::input::{ingest_for_processing, ExternalConditionsInput, HotWaterSourceDetails, Input};
use crate::output::Output;
use crate::read_weather_file::ExternalConditions as ExternalConditionsFromFile;
use crate::simulation_time::SimulationTime;
use crate::statistics::percentile;
use crate::wrappers::future_homes_standard::future_homes_standard::apply_fhs_preprocessing;
use anyhow::anyhow;
use csv::WriterBuilder;
use indexmap::IndexMap;
use lazy_static::lazy_static;
use std::borrow::Cow;
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::io::Read;
use std::sync::Arc;
use wrappers::future_homes_standard::future_homes_standard::apply_fhs_postprocessing;
use wrappers::future_homes_standard::future_homes_standard_fee::{
    apply_fhs_fee_postprocessing, apply_fhs_fee_preprocessing,
};

pub fn run_project(
    input: impl Read,
    output: impl Output,
    external_conditions_data: Option<ExternalConditionsFromFile>,
    _preprocess_only: bool,
    fhs_assumptions: bool,
    fhs_fee_assumptions: bool,
    fhs_not_a_assumptions: bool,
    fhs_not_b_assumptions: bool,
    fhs_fee_not_a_assumptions: bool,
    fhs_fee_not_b_assumptions: bool,
    heat_balance: bool,
    detailed_output_heating_cooling: bool,
) -> Result<(), anyhow::Error> {
    let mut input_for_processing = ingest_for_processing(input)?;

    // do wrapper pre-processing here
    if fhs_assumptions || fhs_not_a_assumptions || fhs_not_b_assumptions {
        apply_fhs_preprocessing(&mut input_for_processing)?;
    } else if fhs_fee_assumptions || fhs_fee_not_a_assumptions || fhs_fee_not_b_assumptions {
        apply_fhs_fee_preprocessing(&mut input_for_processing)?;
    }

    let input = input_for_processing.finalize();

    let external_conditions = external_conditions_from_input(
        input.external_conditions.clone(),
        external_conditions_data,
        input.simulation_time,
    );

    let summary_input_digest: SummaryInputDigest = (&input).into();

    let mut corpus: Corpus = Corpus::from_inputs(&input, Some(external_conditions))?;

    let RunResults {
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
        _heat_source_wet_results_dict,
        _heat_source_wet_results_annual_dict,
        emitters_output_dict: _emitters_output_dict,
        vent_output_list: _vent_output_list,
    } = corpus.run()?;

    write_core_output_file(
        &output,
        OutputFileArgs {
            output_key: "results".to_string(),
            timestep_array: &timestep_array,
            results_totals: &results_totals,
            results_end_user: &results_end_user,
            energy_import: &energy_import,
            energy_export: &energy_export,
            energy_generated_consumed: &energy_generated_consumed,
            energy_to_storage: &energy_to_storage,
            energy_from_storage: &energy_from_storage,
            energy_diverted: &energy_diverted,
            betafactor: &betafactor,
            zone_dict: &zone_dict,
            zone_list: &zone_list,
            hc_system_dict,
            hot_water_dict: &hot_water_dict,
            ductwork_gains,
        },
    )?;

    if heat_balance {
        for (_hb_name, _hb_map) in heat_balance_dict {
            // TODO: write out heat balance files
        }
    }

    if detailed_output_heating_cooling {
        // TODO: write out heat source wet outputs
    }

    // Sum per-timestep figures as needed
    let space_heat_demand_total = zone_dict["space heat demand"]
        .values()
        .map(|v| v.iter().sum::<f64>())
        .sum::<f64>();
    let space_cool_demand_total = zone_dict["space cool demand"]
        .values()
        .map(|v| v.iter().sum::<f64>())
        .sum::<f64>();
    let total_floor_area = corpus.total_floor_area;
    let daily_hw_demand = convert_profile_to_daily(
        match &hot_water_dict["Hot water energy demand incl pipework_loss"] {
            HotWaterResultMap::Float(results) => results
                .get(&KeyString::from("energy_demand_incl_pipework_loss")?)
                .ok_or(anyhow!(
                    "Hot water energy demand incl pipework_loss field not set in hot water output"
                ))?,
            HotWaterResultMap::Int(_) => unreachable!(
                "Hot water energy demand incl pipework_loss is not expected to be an integer"
            ),
        },
        corpus.simulation_time.step_in_hours(),
    );
    let daily_hw_demand_75th_percentile = percentile(&daily_hw_demand, 75);

    write_core_output_file_summary(
        &output,
        SummaryOutputFileArgs {
            output_key: "results_summary".to_string(),
            input: summary_input_digest,
            timestep_array: &timestep_array,
            results_end_user: &results_end_user,
            energy_generated_consumed: &energy_generated_consumed,
            energy_to_storage: &energy_to_storage,
            energy_from_storage: &energy_from_storage,
            energy_diverted: &energy_diverted,
            energy_import: &energy_import,
            energy_export: &energy_export,
            space_heat_demand_total,
            space_cool_demand_total,
            total_floor_area,
            heat_cop_dict,
            cool_cop_dict,
            dhw_cop_dict,
            daily_hw_demand_75th_percentile,
        },
    )?;

    let (heat_transfer_coefficient, heat_loss_parameter, _, _) = corpus.calc_htc_hlp();
    let heat_capacity_parameter = corpus.calc_hcp();
    let heat_loss_form_factor = corpus.calc_hlff();

    write_core_output_file_static(&output, StaticOutputFileArgs {
        output_key: "results_static".to_string(),
        heat_transfer_coefficient,
        heat_loss_parameter,
        heat_capacity_parameter,
        heat_loss_form_factor
    })?;

    if fhs_assumptions || fhs_not_a_assumptions || fhs_not_b_assumptions {
        let notional = fhs_not_a_assumptions || fhs_not_b_assumptions;
        apply_fhs_postprocessing(
            &input,
            &output,
            &energy_import,
            &energy_export,
            &results_end_user,
            &timestep_array,
            notional,
        )?;
    } else if fhs_fee_assumptions || fhs_fee_not_a_assumptions || fhs_fee_not_b_assumptions {
        apply_fhs_fee_postprocessing(
            &output,
            total_floor_area,
            space_heat_demand_total,
            space_cool_demand_total,
        )?;
    }

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
            ec.wind_directions,
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
            input.wind_directions.clone().unwrap_or_default(),
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
    pub static ref UNITS_MAP: IndexMap<&'static str, &'static str> = IndexMap::from([
        ("internal gains", "[W]"),
        ("solar gains", "[W]"),
        ("operative temp", "[deg C]"),
        ("internal air temp", "[deg C]"),
        ("space heat demand", "[kWh]"),
        ("space cool demand", "[kWh]"),
        (
            "DHW: demand volume (including distribution pipework losses)",
            "[litres]"
        ),
        (
            "DHW: demand energy (including distribution pipework losses)",
            "[kWh]"
        ),
        (
            "DHW: demand energy (excluding distribution pipework losses)",
            "[kWh]"
        ),
        ("DHW: total event duration", "[mins]"),
        ("DHW: number of events", "[count]"),
        ("DHW: distribution pipework losses", "[kWh]"),
        ("DHW: primary pipework losses", "[kWh]"),
        ("DHW: storage losses", "[kWh]"),
    ]);
}

struct OutputFileArgs<'a> {
    output_key: String,
    timestep_array: &'a [f64],
    results_totals: &'a IndexMap<KeyString, Vec<f64>>,
    results_end_user: &'a IndexMap<KeyString, IndexMap<String, Vec<f64>>>,
    energy_import: &'a IndexMap<KeyString, Vec<f64>>,
    energy_export: &'a IndexMap<KeyString, Vec<f64>>,
    energy_generated_consumed: &'a IndexMap<KeyString, Vec<f64>>,
    energy_to_storage: &'a IndexMap<KeyString, Vec<f64>>,
    energy_from_storage: &'a IndexMap<KeyString, Vec<f64>>,
    energy_diverted: &'a IndexMap<KeyString, Vec<f64>>,
    betafactor: &'a IndexMap<KeyString, Vec<f64>>,
    zone_dict: &'a IndexMap<&'static str, IndexMap<KeyString, Vec<f64>>>,
    zone_list: &'a [KeyString],
    hc_system_dict: IndexMap<&'static str, IndexMap<String, Vec<f64>>>,
    hot_water_dict: &'a IndexMap<&'static str, HotWaterResultMap>,
    ductwork_gains: IndexMap<KeyString, Vec<f64>>,
}

fn write_core_output_file(output: &impl Output, args: OutputFileArgs) -> Result<(), anyhow::Error> {
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

    // hot_water_dict headings
    for system in hot_water_dict.keys() {
        let mut system = *system;
        if system == "Hot water demand" {
            system = "DHW: demand volume (including distribution pipework losses)";
        }
        if system == "Hot water energy demand" {
            system = "DHW: demand energy (excluding distribution pipework losses)";
        }
        if system == "Hot water energy demand incl pipework_loss" {
            system = "DHW: demand energy (including distribution pipework losses)";
        }
        if system == "Hot water duration" {
            system = "DHW: total event duration";
        }
        if system == "Hot Water Events" {
            system = "DHW: number of events";
        }
        if system == "Pipework losses" {
            system = "DHW: distribution pipework losses";
        }
        if system == "Primary pipework losses" {
            system = "DHW: primary pipework losses";
        }
        if system == "Storage losses" {
            system = "DHW: storage losses";
        }
        headings.push(system.into());
        if UNITS_MAP.contains_key(system) {
            units_row.push(UNITS_MAP[system]);
        } else {
            units_row.push("Unit not defined");
        }
    }

    headings.push("Ventilation: Ductwork gains".into());
    units_row.push("[kWh]");

    for zone in zone_list.iter() {
        for zone_outputs in zone_dict.keys() {
            let zone_headings = format!("{zone}: {zone_outputs}");
            headings.push(zone_headings.into());
            if UNITS_MAP.contains_key(zone_outputs) {
                units_row.push(UNITS_MAP[zone_outputs]);
            } else {
                units_row.push("Unit not defined");
            }
        }
    }

    // hc_system_dict holds heating demand and output as first level keys
    // and the system name as second level keys.
    // Reorganising this dictionary so system names can be grouped together

    // Initialize the reorganized dictionary for grouping systems in hc_system_dict
    let mut reorganised_dict: IndexMap<String, IndexMap<&'static str, Vec<f64>>> =
        Default::default();

    // Iterate over the original map
    for (key, value) in hc_system_dict {
        // Iterate over the nested map
        for (nested_key, nested_value) in value {
            // Add the nested_value to the corresponding entry in reorganized_dict
            reorganised_dict
                .entry(nested_key)
                .or_default()
                .insert(key, nested_value);
        }
    }

    // Loop over reorganised dictionary to add  column and unit headers
    // Check if the system name is set ,else add a designated empty 'None' string
    for system in reorganised_dict.keys() {
        let system_label = if !system.is_empty() {
            system.as_str()
        } else {
            "None"
        };
        for hc_name in reorganised_dict[system].keys() {
            let hc_system = if ["Heating system", "Cooling system"].contains(hc_name) {
                let alternate_name = "energy demand";
                format!("{system_label}: {alternate_name}")
            } else if ["Heating system output", "Cooling system output"].contains(hc_name) {
                let alternate_name = "energy output";
                format!("{system_label}: {alternate_name}")
            } else {
                format!("{system_label}: {hc_name}")
            };
            headings.push(hc_system.into());
            units_row.push("[kWh]");
        }
    }

    for totals_key in results_totals.keys() {
        let totals_header = format!("{totals_key}: total");
        headings.push(totals_header.into());
        units_row.push("[kWh]");
        for end_user_key in results_end_user[totals_key].keys() {
            headings.push(format!("{totals_key}: {end_user_key}").into());
            units_row.push("[kWh]");
        }
        headings.push(format!("{totals_key}: import").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key}: export").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key}: generated and consumed").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key}: beta factor").into());
        units_row.push("[ratio]");
        headings.push(format!("{totals_key}: to storage").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key}: from storage").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key}: diverted").into());
        units_row.push("[kWh]");
    }

    // Write headings and units to output file
    writer.write_record(headings.iter().map(|heading| heading.as_ref()))?;
    writer.write_record(&units_row)?;

    for (t_idx, _timestep) in timestep_array.iter().enumerate() {
        let mut energy_use_row = vec![];
        let mut zone_row = vec![];
        let mut hc_system_row = vec![];
        let mut hw_system_row = vec![];
        let mut hw_system_row_energy = vec![];
        let mut hw_system_row_energy_with_pipework_losses = vec![];
        let mut hw_system_row_duration = vec![];
        let mut hw_system_row_events = vec![];
        let mut pw_losses_row = vec![];
        let mut primary_pw_losses_row = vec![];
        let mut storage_losses_row = vec![];
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
        for zone in zone_list.iter() {
            for zone_outputs in zone_dict.keys() {
                zone_row.push(zone_dict[zone_outputs][zone][t_idx]);
            }
        }

        // Loop over heating and cooling system demand
        for hc_system in reorganised_dict.values() {
            if hc_system.keys().len() == 0 {
                hc_system_row.push(0.);
            }

            for hc_name in hc_system.keys() {
                hc_system_row.push(hc_system[hc_name][t_idx]);
            }
        }

        // Loop over hot water demand
        if let HotWaterResultMap::Float(map) = &hot_water_dict["Hot water demand"] {
            hw_system_row.push(map["demand"][t_idx]);
        }
        if let HotWaterResultMap::Float(map) = &hot_water_dict["Hot water energy demand"] {
            hw_system_row_energy.push(map["energy_demand"][t_idx]);
        }
        if let HotWaterResultMap::Float(map) =
            &hot_water_dict["Hot water energy demand incl pipework_loss"]
        {
            hw_system_row_energy_with_pipework_losses
                .push(map["energy_demand_incl_pipework_loss"][t_idx]);
        }
        if let HotWaterResultMap::Float(map) = &hot_water_dict["Hot water duration"] {
            hw_system_row_duration.push(map["duration"][t_idx]);
        }
        if let HotWaterResultMap::Float(map) = &hot_water_dict["Pipework losses"] {
            pw_losses_row.push(map["pw_losses"][t_idx]);
        }
        if let HotWaterResultMap::Float(map) = &hot_water_dict["Primary pipework losses"] {
            primary_pw_losses_row.push(map["primary_pw_losses"][t_idx]);
        }
        if let HotWaterResultMap::Float(map) = &hot_water_dict["Storage losses"] {
            storage_losses_row.push(map["storage_losses"][t_idx]);
        }
        if let HotWaterResultMap::Int(map) = &hot_water_dict["Hot Water Events"] {
            hw_system_row_events.push(map["no_events"][t_idx]);
        }
        ductwork_row.push(ductwork_gains["ductwork_gains"][t_idx]);

        // create row of outputs and write to output file
        let mut row: Vec<String> = vec![];
        row.append(&mut vec![t_idx.to_string()]);
        row.append(
            &mut hw_system_row
                .into_iter()
                .map(|val| val.to_string())
                .collect(),
        );
        row.append(
            &mut hw_system_row_energy_with_pipework_losses
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
            &mut primary_pw_losses_row
                .into_iter()
                .map(|val| val.to_string())
                .collect(),
        );
        row.append(
            &mut storage_losses_row
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
        row.append(&mut zone_row.into_iter().map(|val| val.to_string()).collect());
        row.append(
            &mut hc_system_row
                .into_iter()
                .map(|val| val.to_string())
                .collect(),
        );
        row.append(
            &mut energy_use_row
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

struct SummaryOutputFileArgs<'a> {
    output_key: String,
    input: SummaryInputDigest,
    timestep_array: &'a [f64],
    results_end_user: &'a IndexMap<KeyString, IndexMap<String, Vec<f64>>>,
    energy_generated_consumed: &'a IndexMap<KeyString, Vec<f64>>,
    energy_to_storage: &'a IndexMap<KeyString, Vec<f64>>,
    energy_from_storage: &'a IndexMap<KeyString, Vec<f64>>,
    energy_diverted: &'a IndexMap<KeyString, Vec<f64>>,
    energy_import: &'a IndexMap<KeyString, Vec<f64>>,
    energy_export: &'a IndexMap<KeyString, Vec<f64>>,
    space_heat_demand_total: f64,
    space_cool_demand_total: f64,
    total_floor_area: f64,
    heat_cop_dict: IndexMap<String, NumberOrDivisionByZero>,
    cool_cop_dict: IndexMap<String, NumberOrDivisionByZero>,
    dhw_cop_dict: IndexMap<String, NumberOrDivisionByZero>,
    daily_hw_demand_75th_percentile: f64,
}

/// A digest of data from the input that is used when generating a summary file.
struct SummaryInputDigest {
    simulation_time: SimulationTime,
    hot_water_source_digests: IndexMap<String, SummaryInputHotWaterSourceDigest>,
}

impl From<&Input> for SummaryInputDigest {
    fn from(input: &Input) -> Self {
        Self {
            simulation_time: input.simulation_time,
            hot_water_source_digests: input
                .hot_water_source
                .as_index_map()
                .iter()
                .map(|(key, details)| (key.clone(), details.into()))
                .collect::<IndexMap<_, _>>(),
        }
    }
}

/// A digest specifically of hot water source data from the input that is used for the summary file.
#[derive(Clone, Copy)]
struct SummaryInputHotWaterSourceDigest {
    source_is_storage_tank: bool,
    source_volume: Option<f64>,
}

impl From<&HotWaterSourceDetails> for SummaryInputHotWaterSourceDigest {
    fn from(value: &HotWaterSourceDetails) -> Self {
        Self {
            source_is_storage_tank: matches!(value, HotWaterSourceDetails::StorageTank { .. }),
            source_volume: value.volume(),
        }
    }
}

fn write_core_output_file_summary(
    output: &impl Output,
    args: SummaryOutputFileArgs,
) -> Result<(), anyhow::Error> {
    let SummaryOutputFileArgs {
        output_key,
        input,
        timestep_array,
        results_end_user,
        energy_generated_consumed,
        energy_to_storage,
        energy_from_storage,
        energy_diverted,
        energy_import,
        energy_export,
        space_heat_demand_total,
        space_cool_demand_total,
        total_floor_area,
        heat_cop_dict,
        cool_cop_dict,
        dhw_cop_dict,
        daily_hw_demand_75th_percentile,
    } = args;

    println!("writing out to {output_key}");

    // Electricity breakdown
    let (elec_generated, elec_consumed) = results_end_user["mains elec"].iter().fold(
        (0.0, 0.0),
        |(elec_generated_acc, elec_consumed_acc), (_end_use, values)| {
            let values_sum = values.iter().sum::<f64>();
            if values_sum < 0.0 {
                (elec_generated_acc + values_sum.abs(), elec_consumed_acc)
            } else {
                (elec_generated_acc, elec_consumed_acc + values_sum)
            }
        },
    );
    let gen_to_consumption = energy_generated_consumed["mains elec"].iter().sum::<f64>();
    let grid_to_consumption = energy_import["mains elec"].iter().sum::<f64>();
    let generation_to_grid = energy_export["mains elec"].iter().sum::<f64>().abs();
    let net_import = grid_to_consumption - generation_to_grid;
    let gen_to_storage = energy_to_storage["mains elec"].iter().sum::<f64>();
    let storage_to_consumption = energy_from_storage["mains elec"].iter().sum::<f64>().abs();
    let gen_to_diverter = energy_diverted["mains elec"].iter().sum::<f64>();
    let storage_eff = if gen_to_storage > 0.0 {
        NumberOrDivisionByZero::Number(storage_to_consumption / gen_to_storage)
    } else {
        NumberOrDivisionByZero::DivisionByZero
    };

    // get peak electricity consumption, and when it happens
    let start_timestep = input.simulation_time.start_time();
    let stepping = input.simulation_time.step;
    // Calculate net import by adding gross import and export figures. Add
    // because export figures are already negative
    let net_import_per_timestep = {
        let mains_energy_import = &energy_import["mains elec"];
        let mains_energy_export = &energy_export["mains elec"];
        (0..timestep_array.len())
            .map(|i| mains_energy_import[i] + mains_energy_export[i])
            .collect::<Vec<_>>()
    };
    let peak_elec_consumption = net_import_per_timestep
        .iter()
        .max_by(|a, b| a.total_cmp(b))
        .unwrap();
    let index_peak_elec_consumption = net_import_per_timestep
        .iter()
        .position(|&x| x == *peak_elec_consumption)
        .unwrap();

    // must reflect hour or half hour in the year (hour 0 to hour 8759)
    // to work with the dictionary below timestep_to_date
    // hence + start_timestep
    let step_peak_elec_consumption = index_peak_elec_consumption as f64 + start_timestep;

    let months_start_end_timesteps = IndexMap::from([
        ("JAN", (0., HOURS_TO_END_JAN / stepping - 1.)),
        (
            "FEB",
            (
                HOURS_TO_END_JAN / stepping,
                HOURS_TO_END_FEB / stepping - 1.,
            ),
        ),
        (
            "MAR",
            (
                HOURS_TO_END_FEB / stepping,
                HOURS_TO_END_MAR / stepping - 1.,
            ),
        ),
        (
            "APR",
            (
                HOURS_TO_END_MAR / stepping,
                HOURS_TO_END_APR / stepping - 1.,
            ),
        ),
        (
            "MAY",
            (
                HOURS_TO_END_APR / stepping,
                HOURS_TO_END_MAY / stepping - 1.,
            ),
        ),
        (
            "JUN",
            (
                HOURS_TO_END_MAY / stepping,
                HOURS_TO_END_JUN / stepping - 1.,
            ),
        ),
        (
            "JUL",
            (
                HOURS_TO_END_JUN / stepping,
                HOURS_TO_END_JUL / stepping - 1.,
            ),
        ),
        (
            "AUG",
            (
                HOURS_TO_END_JUL / stepping,
                HOURS_TO_END_AUG / stepping - 1.,
            ),
        ),
        (
            "SEP",
            (
                HOURS_TO_END_AUG / stepping,
                HOURS_TO_END_SEP / stepping - 1.,
            ),
        ),
        (
            "OCT",
            (
                HOURS_TO_END_SEP / stepping,
                HOURS_TO_END_OCT / stepping - 1.,
            ),
        ),
        (
            "NOV",
            (
                HOURS_TO_END_OCT / stepping,
                HOURS_TO_END_NOV / stepping - 1.,
            ),
        ),
        (
            "DEC",
            (
                HOURS_TO_END_NOV / stepping,
                HOURS_TO_END_DEC / stepping - 1.,
            ),
        ),
    ]);
    let mut timestep_to_date: HashMap<usize, HourForTimestep> = Default::default();
    let mut step = start_timestep;
    for _ in timestep_array {
        for (month, (start, end)) in months_start_end_timesteps.iter() {
            if step <= end.floor() && step >= start.floor() {
                let hour_of_year = step * stepping;
                let hour_start_month = start * stepping;
                let hour_of_month = hour_of_year - hour_start_month;
                // add +1 to day_of_month for first day to be day 1 (not day 0)
                let day_of_month = (hour_of_month / 24.).floor() as usize + 1;
                // add +1 to hour_of_month for first hour to be hour 1 (not hour 0)
                let hour_of_day = ((step % (24. / stepping)) * stepping) + 1.;
                timestep_to_date.insert(
                    step as usize,
                    HourForTimestep {
                        month,
                        day: day_of_month,
                        hour: hour_of_day,
                    },
                );
            }
        }
        step += 1.;
    }

    // Delivered energy by end-use and by fuel
    let mut delivered_energy_map: IndexMap<KeyString, IndexMap<KeyString, f64>> =
        IndexMap::from([(KeyString::from("total").unwrap(), Default::default())]);
    for (fuel, end_uses) in results_end_user {
        if ["_unmet_demand", "hw cylinder"].contains(&fuel.as_str()) {
            continue;
        }
        let mut fuel_results: IndexMap<KeyString, f64> =
            IndexMap::from([(KeyString::from("total").unwrap(), Default::default())]);
        for (end_use, delivered_energy) in end_uses {
            let delivered_energy_sum = delivered_energy.iter().sum::<f64>();
            if delivered_energy_sum >= 0. {
                fuel_results.insert(
                    KeyString::from(end_use).expect("End use was too long to fit in a KeyString."),
                    delivered_energy_sum,
                );
                *fuel_results.get_mut("total").unwrap() += delivered_energy_sum;
                *delivered_energy_map["total"]
                    .entry(
                        KeyString::from(end_use)
                            .expect("End use was too long to fit in a KeyString."),
                    )
                    .or_default() += delivered_energy_sum;
                *delivered_energy_map["total"]
                    .entry(KeyString::from("total").unwrap())
                    .or_default() += delivered_energy_sum;
            }
        }
        delivered_energy_map.insert(*fuel, fuel_results);
    }

    let mut delivered_energy_rows_title =
        vec![
            KeyString::from("Delivered energy by end-use (below) and fuel (right) [kWh/m2]")
                .unwrap(),
        ];
    let mut delivered_energy_rows = vec![vec![StringOrNumber::String(
        KeyString::from("total").unwrap(),
    )]];
    for (fuel, end_uses) in delivered_energy_map {
        delivered_energy_rows_title.push(fuel);
        for row in delivered_energy_rows.iter_mut() {
            row.push(StringOrNumber::Number(0.));
        }
        for (end_use, value) in end_uses {
            let mut end_use_found = false;
            for row in delivered_energy_rows.iter_mut() {
                if row.contains(&StringOrNumber::String(end_use)) {
                    end_use_found = true;
                    *row.get_mut(
                        delivered_energy_rows_title
                            .iter()
                            .position(|&x| x == fuel)
                            .unwrap(),
                    )
                    .unwrap() = StringOrNumber::Number(value / total_floor_area);
                }
            }
            if !end_use_found {
                let mut new_row =
                    vec![StringOrNumber::Number(0.); delivered_energy_rows_title.len()];
                *new_row.get_mut(0).unwrap() = StringOrNumber::String(end_use);
                *new_row
                    .get_mut(
                        delivered_energy_rows_title
                            .iter()
                            .position(|&x| x == fuel)
                            .unwrap(),
                    )
                    .unwrap() = StringOrNumber::Number(value / total_floor_area);
                delivered_energy_rows.push(new_row);
            }
        }
    }

    let heat_cop_rows = heat_cop_dict
        .iter()
        .map(|(h_name, h_cop)| vec![StringOrNumber::from(h_name.as_str()), (*h_cop).into()])
        .collect::<Vec<_>>();
    let cool_cop_rows = cool_cop_dict
        .iter()
        .map(|(c_name, c_cop)| vec![StringOrNumber::from(c_name.as_str()), (*c_cop).into()])
        .collect::<Vec<_>>();
    let mut dhw_cop_rows = dhw_cop_dict
        .iter()
        .map(|(hw_name, hw_cop)| vec![hw_name.as_str().into(), (*hw_cop).into()])
        .collect::<Vec<_>>();

    let writer = output.writer_for_location_key(&output_key)?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    let blank_line: Vec<&'static str> = vec![];

    writer.write_record(["Energy Demand Summary"])?;
    writer.write_record(["", "", "Total"])?;
    writer.write_record([
        "Space heat demand".to_string(),
        "kWh/m2".to_string(),
        (space_heat_demand_total / total_floor_area).to_string(),
    ])?;
    writer.write_record([
        "Space cool demand".to_string(),
        "kWh/m2".to_string(),
        (space_cool_demand_total / total_floor_area).to_string(),
    ])?;
    writer.write_record(&blank_line)?;
    writer.write_record(["Electricity Summary"])?;
    writer.write_record(["", "kWh", "timestep", "month", "day", "hour of day"])?;
    writer.write_record([
        "Peak half-hour consumption".to_string(),
        peak_elec_consumption.to_string(),
        index_peak_elec_consumption.to_string(),
        timestep_to_date[&(step_peak_elec_consumption as usize)]
            .month
            .to_string(),
        timestep_to_date[&(step_peak_elec_consumption as usize)]
            .day
            .to_string(),
        timestep_to_date[&(step_peak_elec_consumption as usize)]
            .hour
            .to_string(),
    ])?;
    writer.write_record(["", "", "Total"])?;
    writer.write_record([
        "Consumption".to_string(),
        "kWh".to_string(),
        elec_consumed.to_string(),
    ])?;
    writer.write_record([
        "Generation".to_string(),
        "kWh".to_string(),
        elec_generated.to_string(),
    ])?;
    writer.write_record([
        "Generation to consumption (immediate excl. diverter)".to_string(),
        "kWh".to_string(),
        gen_to_consumption.to_string(),
    ])?;
    writer.write_record([
        "Generation to storage".to_string(),
        "kWh".to_string(),
        gen_to_storage.to_string(),
    ])?;
    writer.write_record([
        "Generation to diverter".to_string(),
        "kWh".to_string(),
        gen_to_diverter.to_string(),
    ])?;
    writer.write_record([
        "Generation to grid (export)".to_string(),
        "kWh".to_string(),
        generation_to_grid.to_string(),
    ])?;
    writer.write_record([
        "Storage to consumption".to_string(),
        "kWh".to_string(),
        storage_to_consumption.to_string(),
    ])?;
    writer.write_record([
        "Grid to consumption (import)".to_string(),
        "kWh".to_string(),
        grid_to_consumption.to_string(),
    ])?;
    writer.write_record([
        String::from("Net import").into_bytes(),
        String::from("kWh").into_bytes(),
        format!("{}", net_import).into_bytes(),
    ])?;
    writer.write_record([
        "Storage round-trip efficiency".to_string(),
        "ratio".to_string(),
        storage_eff.to_string(),
    ])?;
    writer.write_record(&blank_line)?;
    writer.write_record(["Delivered Energy Summary"])?;
    writer.write_record(
        delivered_energy_rows_title
            .iter()
            .map(|x| format!("{}", x).into_bytes()),
    )?;
    for row in delivered_energy_rows {
        writer.write_record(row.iter().map(|x| format!("{}", x).into_bytes()))?;
    }

    if !dhw_cop_rows.is_empty() {
        writer.write_record(&blank_line)?;
        writer.write_record([
            "Hot water system",
            "Overall CoP",
            "Daily HW demand ([kWh] 75th percentile)",
            "HW cylinder volume (litres)",
        ])?;
        for row in dhw_cop_rows.iter_mut() {
            row.push(StringOrNumber::Number(daily_hw_demand_75th_percentile));
            row.push({
                let hot_water_source: SummaryInputHotWaterSourceDigest = *input
                    .hot_water_source_digests
                    .get(&row[0].to_string())
                    .ok_or_else(|| {
                        anyhow!(
                            "Could not find hot water source digest for row '{}'",
                            row[0]
                        )
                    })?;
                if hot_water_source.source_is_storage_tank {
                    StringOrNumber::Number(hot_water_source.source_volume.unwrap())
                } else {
                    StringOrNumber::String(KeyString::from("N/A").unwrap())
                }
            });
        }
        for row in dhw_cop_rows {
            writer.write_record(row.iter().map(|x| format!("{}", x).into_bytes()))?;
        }
    }
    if !heat_cop_rows.is_empty() {
        writer.write_record(&blank_line)?;
        writer.write_record(["Space heating system", "Overall CoP"])?;
        for row in heat_cop_rows {
            writer.write_record(row.iter().map(|x| format!("{}", x).into_bytes()))?;
        }
    }
    if !cool_cop_rows.is_empty() {
        writer.write_record(&blank_line)?;
        writer.write_record(["Space cooling system", "Overall CoP"])?;
        for row in cool_cop_rows {
            writer.write_record(row.iter().map(|x| format!("{}", x).into_bytes()))?;
        }
    }

    println!("flushing out summary CSV");
    writer.flush()?;

    Ok(())
}

struct StaticOutputFileArgs {
    output_key: String,
    heat_transfer_coefficient: f64,
    heat_loss_parameter: f64,
    heat_capacity_parameter: f64,
    heat_loss_form_factor: f64
}

fn write_core_output_file_static(output: &impl Output, args: StaticOutputFileArgs) -> Result<(), anyhow::Error> {
    let StaticOutputFileArgs {
        output_key,
        heat_transfer_coefficient,
        heat_loss_parameter,
        heat_capacity_parameter,
        heat_loss_form_factor
    } = args;

    println!("writing out to {output_key}");

    let writer = output.writer_for_location_key(&output_key)?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    writer.write_record(["Heat transfer coefficient".to_owned(), "W / K".to_owned(), heat_transfer_coefficient.to_string() ])?;
    writer.write_record(["Heat loss parameter".to_owned(), "W / m2.K".to_owned(), heat_loss_parameter.to_string() ])?;
    writer.write_record(["Heat capacity parameter".to_owned(), "kJ / m2.K".to_owned(), heat_capacity_parameter.to_string() ])?;
    writer.write_record(["Heat loss form factor".to_owned(), "".to_owned(), heat_loss_form_factor.to_string() ])?;

    println!("flushing out static CSV");
    writer.flush()?;

    Ok(())
}

const HOURS_TO_END_JAN: f64 = 744.;
const HOURS_TO_END_FEB: f64 = 1416.;
const HOURS_TO_END_MAR: f64 = 2160.;
const HOURS_TO_END_APR: f64 = 2880.;
const HOURS_TO_END_MAY: f64 = 3624.;
const HOURS_TO_END_JUN: f64 = 4344.;
const HOURS_TO_END_JUL: f64 = 5088.;
const HOURS_TO_END_AUG: f64 = 5832.;
const HOURS_TO_END_SEP: f64 = 6552.;
const HOURS_TO_END_OCT: f64 = 7296.;
const HOURS_TO_END_NOV: f64 = 8016.;
const HOURS_TO_END_DEC: f64 = 8760.;

struct HourForTimestep {
    month: &'static str,
    day: usize,
    hour: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum StringOrNumber {
    String(KeyString),
    Number(f64),
}

impl From<NumberOrDivisionByZero> for StringOrNumber {
    fn from(value: NumberOrDivisionByZero) -> Self {
        match value {
            NumberOrDivisionByZero::Number(number) => StringOrNumber::Number(number),
            NumberOrDivisionByZero::DivisionByZero => {
                StringOrNumber::String(KeyString::from("DIV/0").unwrap())
            }
        }
    }
}

impl Display for StringOrNumber {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                StringOrNumber::String(string) => Cow::Borrowed(string.as_str()),
                StringOrNumber::Number(number) => Cow::Owned(number.to_string()),
            }
        )
    }
}

impl From<StringOrNumber> for Vec<u8> {
    fn from(value: StringOrNumber) -> Self {
        format!("{}", value).into_bytes()
    }
}

impl From<StringOrNumber> for KeyString {
    fn from(value: StringOrNumber) -> Self {
        KeyString::from(format!("{}", value).as_str()).unwrap()
    }
}

impl From<&str> for StringOrNumber {
    fn from(value: &str) -> Self {
        StringOrNumber::String(KeyString::from(value).unwrap())
    }
}

impl From<f64> for StringOrNumber {
    fn from(value: f64) -> Self {
        StringOrNumber::Number(value)
    }
}
