#![allow(
    clippy::too_many_arguments,
    clippy::doc_overindented_list_items,
    dead_code
)]

#[macro_use]
extern crate is_close;
pub mod compare_floats;
pub mod core;
pub mod corpus;
pub mod errors;

mod hem_core;

pub mod input;
mod output;
pub mod output_writer;
pub mod read_weather_file;
pub mod statistics;

use crate::core::units::{convert_profile_to_daily, WATTS_PER_KILOWATT};
use crate::corpus::{
    Corpus, NumberOrDivisionByZero, OutputOptions, ResultsAnnual, ResultsPerTimestep,
};
use crate::errors::{HemCoreError, HemError, NotImplementedError};
use crate::external_conditions::ExternalConditions;
use crate::input::{ExternalConditionsInput, HotWaterSourceDetails, Input};
use crate::output::{Output, OutputEmitters, OutputStatic, OUTPUT_ZONE_DATA_FIELD_HEADINGS};
use crate::output_writer::OutputWriter;
use crate::read_weather_file::ExternalConditions as ExternalConditionsFromFile;
use crate::simulation_time::SimulationTime;
use anyhow::{anyhow, bail};
use convert_case::{Case, Casing};
use csv::WriterBuilder;
use erased_serde::Serialize as ErasedSerialize;
use hem_core::external_conditions;
use hem_core::simulation_time;
use indexmap::IndexMap;
use serde::{Serialize, Serializer};
use smartstring::alias::String;
use std::borrow::Cow;
use std::fmt::{Debug, Display, Formatter};
use std::io::Read;
use std::ops::AddAssign;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::{Arc, LazyLock};
use tracing::{debug, instrument};
use zerocopy::IntoBytes;

pub const HEM_VERSION: &str = "1.0.0a1";
pub const HEM_VERSION_DATE: &str = "2025-10-02";

#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
pub enum OutputFormat {
    Json,
    Csv,
}

#[derive(Serialize)]
pub struct HemResponse {
    #[serde(flatten)]
    payload: Box<dyn ErasedSerialize + Send + Sync + 'static>,
}

impl HemResponse {
    pub fn new(payload: impl ErasedSerialize + Send + Sync + 'static) -> Self {
        Self {
            payload: Box::new(payload),
        }
    }

    pub fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.payload.serialize(serializer)
    }
}

#[instrument(skip_all)]
pub fn run_project_from_input_file(
    input: impl Read,
    output_writer: &impl OutputWriter,
    external_conditions_data: Option<ExternalConditionsFromFile>,
    output_formats: Option<&Vec<OutputFormat>>,
    tariff_data_file: Option<&str>,
    heat_balance: bool,
    detailed_output_heating_cooling: bool,
) -> Result<CalculationResult, HemError> {
    #[instrument(skip_all)]
    fn finalize(input: impl Read) -> anyhow::Result<Input> {
        let input = serde_json::from_reader(input)?;
        // NB. this _might_ in time be a good point to perform a validation against the core schema - or it might not
        // if let BasicOutput::Invalid(errors) =
        //     CORE_INCLUDING_FHS_VALIDATOR.apply(&self.input).basic()
        // {
        //     bail!(
        //         "Wrapper formed invalid JSON for the core schema: {}",
        //         serde_json::to_value(errors)?.to_json_string_pretty()?
        //     );
        // }

        serde_json::from_value(input).map_err(|err| anyhow!(err))
    }
    let input = finalize(input)?;

    let results = run_project(
        input,
        external_conditions_data,
        tariff_data_file,
        heat_balance,
        detailed_output_heating_cooling,
    )?;

    if let Some(output_formats) = output_formats {
        let steps_in_hours = results.input.simulation_time.step;

        write_core_output_files(
            &results.output,
            results.input.as_ref(),
            output_writer,
            output_formats,
            steps_in_hours,
            heat_balance,
            detailed_output_heating_cooling,
        )?;
    }

    Ok(results)
}

#[instrument(skip_all)]
pub fn run_project(
    input: Input,
    external_conditions_data: Option<ExternalConditionsFromFile>,
    tariff_data_file: Option<&str>,
    heat_balance: bool,
    detailed_output_heating_cooling: bool,
) -> Result<CalculationResult, HemError> {
    catch_unwind(AssertUnwindSafe(|| {
        #[instrument(skip_all)]
        fn merge_external_conditions_data(
            mut input: Input,
            external_conditions_data: Option<ExternalConditionsInput>,
        ) -> anyhow::Result<Input> {
            if external_conditions_data.is_none()
                && !input.external_conditions.are_all_fields_set() {
                bail!("No weather data found. Please provide a weather file or complete weather data in the input file.");
            }
            if let Some(external_conditions) = external_conditions_data {
                let shading_segments = &input.external_conditions.shading_segments;
                let mut new_external_conditions = external_conditions.clone();
                new_external_conditions.shading_segments = shading_segments.clone();
                input.external_conditions = Arc::from(new_external_conditions);
            }

            Ok(input)
        }

        let input: Arc<Input> = merge_external_conditions_data(input, external_conditions_data
            .as_ref()
            .map(ExternalConditionsInput::from))?
            .into();

        // 1. Determine external conditions to use for calculations.
        #[instrument(skip_all)]
        fn resolve_external_conditions(
            input: &Input,
        ) -> ExternalConditions {
            external_conditions_from_input(
                input.external_conditions.clone(),
                input.simulation_time,
            )
        }

        let corpus = {
            let external_conditions = resolve_external_conditions(&input);

            // 2. Build corpus from input and external conditions.
            #[instrument(skip_all)]
            fn build_corpus(
                input: Arc<Input>,
                external_conditions: &ExternalConditions,
                tariff_data_file: Option<&str>,
                print_heat_balance: bool,
                detailed_output_heating_cooling: bool,
            ) -> anyhow::Result<Corpus> {
                let output_options = OutputOptions { print_heat_balance, detailed_output_heating_cooling };
                Corpus::from_inputs(input, Some(external_conditions), tariff_data_file, &output_options)
            }

            build_corpus(input.clone(), &external_conditions, tariff_data_file, heat_balance, detailed_output_heating_cooling).map_err(|e| {
                capture_specific_error_case(&e).unwrap_or_else(|| HemError::InvalidRequest(e))
            })?
        };

        // 3. Run HEM calculation.
        #[instrument(skip_all)]
        fn run_hem_calculation(
            corpus: &Corpus,
        ) -> anyhow::Result<Output> {
            corpus.run()
        }

        // catch_unwind here catches any downstream panics so we can at least map to the right HemError variant
        let output = match catch_unwind(AssertUnwindSafe(|| {
            run_hem_calculation(&corpus).map_err(|e| {
                capture_specific_error_case(&e)
                    .unwrap_or_else(|| HemError::FailureInCalculation(HemCoreError::new(e)))
            })
        })) {
            Ok(results) => results?,
            Err(panic) => {
                return Err(HemError::PanicInCalculation(
                    panic
                        .downcast_ref::<&str>()
                        .map_or("Error not captured", |v| v)
                        .to_owned(),
                ))
            }
        };

        Ok(CalculationResult::new(input, output))
    }))
        .map_err(|e| {
            HemError::GeneralPanic(
                e.downcast_ref::<&str>()
                    .map_or("Uncaught panic - could not capture more information; to debug further, try replicating in debug mode.", |v| v)
                    .to_owned(),
            )
        })?
}

#[instrument(skip_all)]
fn write_core_output_files(
    output: &Output,
    primary_input: &Input,
    output_writer: &impl OutputWriter,
    output_formats: &[OutputFormat],
    hour_per_step: f64,
    heat_balance: bool,
    detailed_output_heating_cooling: bool,
) -> anyhow::Result<()> {
    if output_writer.is_noop() {
        return Ok(());
    }

    if output_formats.contains(&OutputFormat::Json) {
        write_output_json_file("output", output, output_writer)?;
    }

    if output_formats.contains(&OutputFormat::Csv) {
        let input = primary_input;

        write_core_output_file_static(&output.static_, "results_static", output_writer)?;

        write_core_output_file(output, "results", output_writer)?;

        write_core_output_file_summary(output, "results_summary", output_writer, input)?;

        if heat_balance {
            for (hb_name, hb_map) in &output.core.heat_balance_all {
                let output_key = format!("results_heat_balance_{}", hb_name.to_case(Case::Snake));
                write_core_output_file_heat_balance(
                    output_key.as_str(),
                    &output.core.timestep_array,
                    hour_per_step,
                    hb_map,
                    output_writer,
                )?;
            }
        }

        if detailed_output_heating_cooling {
            for (heat_source_wet_name, heat_source_wet_results) in
                &output.core.heat_source_wet_results
            {
                let output_key = format!("results_heat_source_wet__{heat_source_wet_name}");
                write_core_output_file_heat_source_wet(
                    output_key.as_str(),
                    &output.core.timestep_array,
                    heat_source_wet_results,
                    output_writer,
                )?;
            }

            for (heat_source_wet_name, heat_source_wet_results_annual) in
                &output.core.heat_source_wet_results_annual
            {
                let output_key = format!("results_heat_source_wet__{heat_source_wet_name}");
                write_core_output_file_heat_source_wet_summary(
                    output_key.as_str(),
                    heat_source_wet_results_annual,
                    output_writer,
                )?;
            }

            // Function call to write detailed ventilation results
            let vent_output_key = "ventilation_results";
            write_core_output_file_ventilation_detailed(
                vent_output_key,
                &output.core.ventilation,
                output_writer,
            )?;

            for (hot_water_source_name, hot_water_source_results) in
                &output.core.hot_water_source_results_summary
            {
                let hot_water_source_file = format!(
                    "results_hot_water_source_summary__{}",
                    hot_water_source_name.replace(" ", "_")
                );
                write_core_output_file_hot_water_source_summary(
                    hot_water_source_file.as_str(),
                    hot_water_source_results,
                    output_writer,
                )?;
            }

            // Create a file for emitters detailed output and write
            let emitters_output_prefix = "results_emitters_";
            write_core_output_file_emitters_detailed(
                emitters_output_prefix,
                &output.core.emitters,
                output_writer,
            )?;

            // Create a file for esh detailed output and write
            let esh_output_prefix = "results_esh_";
            write_core_output_file_esh_detailed(
                esh_output_prefix,
                &output.core.electric_storage_heaters,
                output_writer,
            )?;
        }
    }

    Ok(())
}

fn capture_specific_error_case(e: &anyhow::Error) -> Option<HemError> {
    if let Some(e) = e.downcast_ref::<NotImplementedError>() {
        return Some(HemError::NotImplemented(e.clone()));
    }

    None
}

pub struct CalculationResult {
    pub output: Output,
    pub input: Arc<Input>,
}

impl CalculationResult {
    fn new(input: Arc<Input>, output: Output) -> CalculationResult {
        Self { output, input }
    }
}

pub fn external_conditions_from_input(
    input: Arc<ExternalConditionsInput>,
    simulation_time: SimulationTime,
) -> ExternalConditions {
    ExternalConditions::new(
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
        None,
        false,
        input.direct_beam_conversion_needed.unwrap_or(false),
        input.shading_segments.clone(),
    )
}

pub static UNITS_MAP: LazyLock<IndexMap<&'static str, &'static str>> = LazyLock::new(|| {
    IndexMap::from([
        ("internal gains", "[W]"),
        ("solar gains", "[W]"),
        ("operative temp", "[deg C]"),
        ("internal air temp", "[deg C]"),
        ("space heat demand", "[kWh]"),
        ("space cool demand", "[kWh]"),
        (
            "hot water volume required from hot water source",
            "[litres]",
        ),
        ("hot water energy demand at hot water source", "[kWh]"),
        (
            "hot water energy demand at connected tapping points",
            "[kWh]",
        ),
        ("total event duration", "[mins]"),
        ("number of events", "[count]"),
        ("distribution pipework losses", "[kWh]"),
        ("primary pipework losses", "[kWh]"),
        ("storage losses", "[kWh]"),
    ])
});

type ReorganisedMapForOutput = IndexMap<Option<Arc<str>>, IndexMap<Arc<str>, Vec<f64>>>;

fn write_core_output_file(
    output: &Output,
    output_key: &str,
    output_writer: &impl OutputWriter,
) -> anyhow::Result<()> {
    let writer = output_writer.writer_for_location_key(output_key, "csv")?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    let mut headings: Vec<Cow<'static, str>> = vec!["Timestep".into()];
    let mut units_row = vec!["[count]"];

    // Headings mapping for each hot water system.
    for (field_alias, field_keys) in output.core.hot_water_systems.fields() {
        for hws_name in field_keys {
            headings.push(format!("{hws_name}: {field_alias}").into());
            units_row.push(
                UNITS_MAP
                    .get(field_alias.to_string().as_str())
                    .unwrap_or(&"Unit not defined (add to UNITS_MAP)"),
            );
        }
    }

    headings.push("Ventilation: Ductwork gains".into());
    units_row.push("[kWh]");

    for zone in output.core.zone_list.iter() {
        for field_name in OUTPUT_ZONE_DATA_FIELD_HEADINGS {
            let zone_headings = format!("{zone}: {field_name}").into();
            headings.push(zone_headings);
            units_row.push(
                UNITS_MAP
                    .get(field_name)
                    .unwrap_or(&"Unit not defined (add to UNITS_MAP)"),
            );
        }
    }

    // OutputHeatingCoolingSystem holds heating demand and output as first level keys
    // and the system name as second level keys.
    // Reorganising this dictionary so system names can be grouped together

    // Initialize the reorganized dictionary for grouping systems from OutputHeatingCoolingSystem
    let mut reorganised_dict: ReorganisedMapForOutput = Default::default();

    // Iterate over the original structures
    for (key, value) in [
        (
            "heating_system_output",
            &output.core.heating_cooling_system.heating_system_output,
        ),
        (
            "cooling_system_output",
            &output.core.heating_cooling_system.cooling_system_output,
        ),
    ] {
        // Iterate over the nested map
        for (nested_key, nested_value) in value {
            // Check if the nested_key already exists in reorganized_dict
            reorganised_dict
                .entry(nested_key.clone())
                .or_default()
                .insert(key.into(), nested_value.clone());
        }
    }

    // Loop over reorganised dictionary to add column and unit headers
    // Check if the system name is set, else add a designated empty 'None' string
    for system in reorganised_dict.keys() {
        if let Some(system) = system {
            for hc_name in reorganised_dict[&Some(system.clone())].keys() {
                let hc_system = if ["heating_system_output", "cooling_system_output"]
                    .contains(&hc_name.as_ref())
                {
                    let alternate_name = "energy output";
                    format!("{system}: {alternate_name}")
                } else {
                    format!("{system}: {hc_name}")
                };
                headings.push(hc_system.into());
                units_row.push("[kWh]");
            }
        } else {
            for hc_name in reorganised_dict[&None].keys() {
                let hc_system = if ["heating_system_output", "cooling_system_output"]
                    .contains(&hc_name.as_ref())
                {
                    let alternate_name = "energy output";
                    format!("None: {alternate_name}")
                } else {
                    format!("None: {hc_name}")
                };
                headings.push(hc_system.into());
                units_row.push("[kWh]");
            }
        }
    }

    for totals_key in output.core.results_totals.keys() {
        let totals_header = format!("{totals_key}: total");
        headings.push(totals_header.into());
        units_row.push("[kWh]");
        for end_user_key in output.core.results_end_user[totals_key].keys() {
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
        headings.push(format!("{totals_key}: generation to storage").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key}: from storage").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key}: grid to storage").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key}: battery charge level").into());
        units_row.push("[ratio]");
        headings.push(format!("{totals_key}: diverted").into());
        units_row.push("[kWh]");
    }

    // Write headings and units to output file
    writer.write_record(headings.iter().map(|heading| heading.as_ref()))?;
    writer.write_record(&units_row)?;

    for t_idx in 0..output.core.timestep_array.len() {
        let mut energy_use_row = vec![];
        let mut zone_row = vec![];
        let mut hc_system_row: Vec<f64> = vec![];
        let mut dhw_row: Vec<f64> = vec![];
        let energy_shortfall: Vec<f64> = vec![]; // TODO (from Python) this seem to never be populated ???

        for totals_key in output.core.results_totals.keys() {
            energy_use_row.push(output.core.results_totals[totals_key][t_idx]);
            for end_user_key in output.core.results_end_user[totals_key].keys() {
                energy_use_row.push(output.core.results_end_user[totals_key][end_user_key][t_idx]);
            }
            energy_use_row.push(output.core.energy_import[totals_key][t_idx]);
            energy_use_row.push(output.core.energy_export[totals_key][t_idx]);
            energy_use_row.push(output.core.generation_to_grid[totals_key][t_idx]);
            energy_use_row.push(output.core.energy_generated_consumed[totals_key][t_idx]);
            energy_use_row.push(output.core.beta_factor[totals_key][t_idx]);
            energy_use_row.push(output.core.energy_to_storage[totals_key][t_idx]);
            energy_use_row.push(output.core.energy_from_storage[totals_key][t_idx]);
            energy_use_row.push(output.core.storage_from_grid[totals_key][t_idx]);
            energy_use_row.push(output.core.battery_state_of_charge[totals_key][t_idx]);
            energy_use_row.push(output.core.energy_diverted[totals_key][t_idx]);
        }

        // Loop over results separated by zone
        for zone in output.core.zone_list.iter() {
            for zone_outputs in output.core.zone_data.ordered_values_for_zone(zone) {
                zone_row.push(zone_outputs[t_idx]);
            }
        }

        // Loop over system names and print the heating and cooling energy demand and output
        for system in reorganised_dict.keys() {
            let system_key = &system.as_ref().map(|x| x.clone());
            for hc_name in reorganised_dict[system_key].keys() {
                hc_system_row.push(reorganised_dict[system_key][hc_name][t_idx]);
            }
        }

        // Column order determined by field order on OutputHotWaterSystems
        for field in output.core.hot_water_systems.ordered_values() {
            for hws_name in field.keys() {
                dhw_row.push(field[hws_name][t_idx]);
            }
        }

        let mut row: Vec<StringOrNumber> = vec![t_idx.into()];
        row.append(&mut dhw_row.into_iter().map(StringOrNumber::from).collect());
        row.append(&mut vec![output.core.ductwork_gains[t_idx].into()]);
        row.append(
            &mut energy_shortfall
                .into_iter()
                .map(StringOrNumber::from)
                .collect(),
        );
        row.append(&mut zone_row.into_iter().map(StringOrNumber::from).collect());
        row.append(
            &mut hc_system_row
                .into_iter()
                .map(StringOrNumber::from)
                .collect(),
        );
        row.append(
            &mut energy_use_row
                .into_iter()
                .map(StringOrNumber::from)
                .collect(),
        );

        writer.write_record(row.iter().map(StringOrNumber::as_bytes))?;
    }

    Ok(())
}

fn write_core_output_file_summary(
    output: &Output,
    output_key: &str,
    output_writer: &impl OutputWriter,
    input: &Input,
) -> Result<(), anyhow::Error> {
    if output_writer.is_noop() {
        return Ok(());
    }
    debug!("writing out to {output_key}");

    let mut delivered_energy_rows_title =
        vec!["Delivered energy by end-use (below) and fuel (right) [kWh/m2]".into()];
    let mut delivered_energy_rows = vec![vec![StringOrNumber::String("total".into())]];

    for (fuel, end_uses) in output.summary.delivered_energy_by_floor_area() {
        delivered_energy_rows_title.push(fuel.clone());
        let index = delivered_energy_rows_title.len() - 1;
        for row in delivered_energy_rows.iter_mut() {
            row.push(StringOrNumber::Float(0.));
        }
        for (end_use, value) in end_uses {
            let mut end_use_found = false;
            for row in delivered_energy_rows.iter_mut() {
                if row.contains(&StringOrNumber::String(end_use.to_string().into())) {
                    end_use_found = true;
                    row[index] = StringOrNumber::Float(value);
                }
            }
            if !end_use_found {
                let mut new_row =
                    vec![StringOrNumber::Float(0.); delivered_energy_rows_title.len()];
                *new_row.get_mut(0).unwrap() = StringOrNumber::String(end_use.to_string().into());
                new_row[index] = StringOrNumber::Float(value);
                delivered_energy_rows.push(new_row);
            }
        }
    }

    // Output "DIV/0" in place of None
    let heat_cop_rows = output
        .core
        .cop
        .space_heating_system
        .iter()
        .map(|(h_name, h_cop)| {
            vec![
                StringOrNumber::from(h_name.to_string()),
                match h_cop {
                    NumberOrDivisionByZero::Number(number) => number.into(),
                    NumberOrDivisionByZero::DivisionByZero => "DIV/0".into(),
                },
            ]
        })
        .collect::<Vec<_>>();
    let cool_cop_rows = output
        .core
        .cop
        .space_cooling_system
        .iter()
        .map(|(c_name, c_cop)| {
            vec![
                StringOrNumber::from(c_name.to_string()),
                match c_cop {
                    NumberOrDivisionByZero::Number(number) => number.into(),
                    NumberOrDivisionByZero::DivisionByZero => "DIV/0".into(),
                },
            ]
        })
        .collect::<Vec<_>>();
    let mut dhw_cop_rows = output
        .core
        .cop
        .hot_water_system
        .iter()
        .map(|(hw_name, hw_cop)| {
            vec![
                StringOrNumber::from(hw_name.to_string()),
                match hw_cop {
                    NumberOrDivisionByZero::Number(number) => number.into(),
                    NumberOrDivisionByZero::DivisionByZero => "DIV/0".into(),
                },
            ]
        })
        .collect::<Vec<_>>();

    let writer = output_writer.writer_for_location_key(output_key, "csv")?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    let blank_line: Vec<&'static str> = vec![];

    writer.write_record(["Energy Demand Summary"])?;
    writer.write_record(["", "", "Total"])?;
    writer.write_record([
        "Space heat demand".to_string(),
        "kWh/m2".to_string(),
        output.summary.space_heat_demand_by_floor_area().to_string(),
    ])?;
    writer.write_record([
        "Space cool demand".to_string(),
        "kWh/m2".to_string(),
        output.summary.space_cool_demand_by_floor_area().to_string(),
    ])?;
    writer.write_record(&blank_line)?;
    writer.write_record(["Energy Supply Summary"])?;
    writer.write_record(["", "kWh", "timestep", "month", "day", "hour of day"])?;

    let peak_consumption = &output.summary.electricity_peak_consumption;
    writer.write_record([
        "Peak half-hour consumption (electricity)".to_string(), // TODO (from Python) technically per-step, not half-hour
        peak_consumption.peak.to_string(),
        peak_consumption.index.to_string(),
        peak_consumption.index.to_string(),
        MONTH_NAMES[peak_consumption.month as usize].to_string(),
        peak_consumption.day.to_string(),
        peak_consumption.hour.to_string(),
    ])?;
    writer.write_record(&blank_line)?;

    let mut header_row = vec!["".to_owned(), "Total".to_owned()];
    header_row.extend(output.summary.energy_supply.keys().map(|x| x.to_string()));
    writer.write_record(&header_row)?;
    let fields = [
        // Label, unit, OutputSummaryEnergySupply field
        (
            "Consumption",
            "kWh",
            EnergySupplyStatKey::ElectricityConsumed,
        ),
        (
            "Generation",
            "kWh",
            EnergySupplyStatKey::ElectricityGenerated,
        ),
        (
            "Generation to consumption (immediate excl. diverter)",
            "kWh",
            EnergySupplyStatKey::GenerationToConsumption,
        ),
        (
            "Generation to storage",
            "kWh",
            EnergySupplyStatKey::GenerationToStorage,
        ),
        (
            "Generation to diverter",
            "kWh",
            EnergySupplyStatKey::GenerationToDiverter,
        ),
        (
            "Generation to grid (export)",
            "kWh",
            EnergySupplyStatKey::GenerationToGrid,
        ),
        (
            "Storage to consumption",
            "kWh",
            EnergySupplyStatKey::StorageToConsumption,
        ),
        (
            "Grid to storage",
            "kWh",
            EnergySupplyStatKey::StorageFromGrid,
        ),
        (
            "Grid to consumption (import)",
            "kWh",
            EnergySupplyStatKey::GridToConsumption,
        ),
        ("Net import", "kWh", EnergySupplyStatKey::NetImport),
        (
            "Storage round-trip efficiency",
            "ratio",
            EnergySupplyStatKey::StorageEfficiency,
        ),
    ];
    for (label, unit, field) in fields {
        let mut row: Vec<std::string::String> = vec![label.into(), unit.into()];
        for stat in output.summary.energy_supply.values() {
            let value = stat.field(&field);
            let value = if field == EnergySupplyStatKey::StorageEfficiency && value.is_nan() {
                "DIV/0".into()
            } else {
                value.to_string()
            };
            row.push(value);
        }
        writer.write_record(&row)?;
    }
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
            let hws_name: Arc<str> = row[0].clone().into();
            row.push(StringOrNumber::Float(
                output.summary.hot_water_demand_daily_75th_percentile[hws_name.as_ref()],
            ));

            row.push({
                if let HotWaterSourceDetails::StorageTank { volume, .. } = input
                    .hot_water_source
                    .get(&hws_name.to_string())
                    .ok_or_else(|| {
                        anyhow!("Could not find hot water source with name '{hws_name}'")
                    })?
                {
                    volume.into()
                } else {
                    StringOrNumber::String("N/A".into())
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

    debug!("flushing out summary CSV");
    writer.flush()?;

    Ok(())
}

const MONTH_NAMES: [&str; 12] = [
    "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC",
];

fn write_output_json_file(
    output_key: &str,
    output: &Output,
    output_writer: &impl OutputWriter,
) -> anyhow::Result<()> {
    if output_writer.is_noop() {
        return Ok(());
    }

    debug!("writing JSON output file");

    let writer = output_writer.writer_for_location_key(output_key, "json")?;
    serde_json::to_writer_pretty(writer, &output)?;

    Ok(())
}

#[derive(Clone, Copy)]
struct EnergySupplyStat {
    elec_generated: f64,
    elec_consumed: f64,
    gen_to_consumption: f64,
    grid_to_consumption: f64,
    generation_to_grid: f64,
    net_import: f64,
    gen_to_storage: f64,
    storage_to_consumption: f64,
    storage_from_grid: f64,
    gen_to_diverter: f64,
    storage_eff: NumberOrDivisionByZero,
}

impl EnergySupplyStat {
    fn display_for_key(&self, key: &EnergySupplyStatKey) -> String {
        match key {
            EnergySupplyStatKey::ElectricityGenerated => self.elec_generated.to_string().into(),
            EnergySupplyStatKey::ElectricityConsumed => self.elec_consumed.to_string().into(),
            EnergySupplyStatKey::GenerationToConsumption => {
                self.gen_to_consumption.to_string().into()
            }
            EnergySupplyStatKey::GridToConsumption => self.grid_to_consumption.to_string().into(),
            EnergySupplyStatKey::GenerationToGrid => self.generation_to_grid.to_string().into(),
            EnergySupplyStatKey::NetImport => self.net_import.to_string().into(),
            EnergySupplyStatKey::GenerationToStorage => self.gen_to_storage.to_string().into(),
            EnergySupplyStatKey::StorageToConsumption => {
                self.storage_to_consumption.to_string().into()
            }
            EnergySupplyStatKey::StorageFromGrid => self.storage_from_grid.to_string().into(),
            EnergySupplyStatKey::GenerationToDiverter => self.gen_to_diverter.to_string().into(),
            EnergySupplyStatKey::StorageEfficiency => self.storage_eff.to_string().into(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum EnergySupplyStatKey {
    ElectricityGenerated,
    ElectricityConsumed,
    GenerationToConsumption,
    GridToConsumption,
    GenerationToGrid,
    NetImport,
    GenerationToStorage,
    StorageToConsumption,
    StorageFromGrid,
    GenerationToDiverter,
    StorageEfficiency,
}

fn write_core_output_file_static(
    output: &OutputStatic,
    output_key: &str,
    output_writer: &impl OutputWriter,
) -> Result<(), anyhow::Error> {
    debug!("writing out to {output_key}");

    let writer = output_writer.writer_for_location_key(output_key, "csv")?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    writer.write_record([
        "Heat transfer coefficient".to_owned(),
        "W / K".to_owned(),
        output.heat_transfer_coefficient.to_string(),
    ])?;
    writer.write_record([
        "Heat loss parameter".to_owned(),
        "W / m2.K".to_owned(),
        output.heat_loss_param.to_string(),
    ])?;
    writer.write_record([
        "Heat capacity parameter".to_owned(),
        "kJ / m2.K".to_owned(),
        output.heat_capacity_param.to_string(),
    ])?;
    writer.write_record([
        "Heat loss form factor".to_owned(),
        "".to_owned(),
        output.heat_loss_form_factor.to_string(),
    ])?;
    writer.write_record(["Assumptions used for HTC/HLP calculation:"])?;
    writer.write_record([
        "Internal air temperature".to_owned(),
        "Celsius".to_owned(),
        output.temperature_air_internal.to_string(),
    ])?;
    writer.write_record([
        "External air temperature".to_owned(),
        "Celsius".to_owned(),
        output.temperature_air_external.to_string(),
    ])?;

    debug!("flushing out static CSV");
    writer.flush()?;

    Ok(())
}

fn write_core_output_file_heat_balance(
    output_key: &str,
    timestep_array: &[f64],
    hour_per_step: f64,
    heat_balance_map: &IndexMap<Arc<str>, IndexMap<Arc<str>, Vec<f64>>>,
    output_writer: &impl OutputWriter,
) -> Result<(), anyhow::Error> {
    let writer = output_writer.writer_for_location_key(output_key, "csv")?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    let mut headings = vec!["Timestep".to_string()];
    let mut units_row = vec!["index".to_string()];
    let mut rows = vec![vec![StringOrNumber::String("".into())]];

    let mut headings_annual = vec!["".to_string()];
    let mut units_annual = vec!["".to_string()];

    let mut number_of_zones = 0usize;

    for (z_name, heat_loss_gain_map) in heat_balance_map {
        for heat_loss_gain_name in heat_loss_gain_map.keys() {
            headings.push(format!("{}: {}", z_name, heat_loss_gain_name));
            units_row.push("[W]".to_string());
        }
        number_of_zones += 1;
    }

    let mut annual_totals: Vec<StringOrNumber> = vec![
        StringOrNumber::Float(0.);
        heat_balance_map
            .values()
            .map(|x| x.len())
            .max()
            .ok_or_else(|| anyhow!(
                "Could not write heat balance file as heat balance data was empty."
            ))?
            * number_of_zones
    ];
    annual_totals.insert(0, StringOrNumber::String("".into()));

    for (z_name, heat_loss_gain_map) in heat_balance_map {
        for heat_loss_gain_name in heat_loss_gain_map.keys() {
            headings_annual.push(format!("{}: total {}", z_name, heat_loss_gain_name));
            units_annual.push("[kWh]".to_string());
        }
    }

    for t_idx in 0..timestep_array.len() {
        let mut row = vec![StringOrNumber::Integer(t_idx)];
        let mut annual_totals_index = 1;
        for heat_loss_gain_map in heat_balance_map.values() {
            for heat_loss_gain_value in heat_loss_gain_map.values() {
                row.push(StringOrNumber::Float(heat_loss_gain_value[t_idx]));
                annual_totals[annual_totals_index] += StringOrNumber::Float(
                    heat_loss_gain_value[t_idx] * hour_per_step / WATTS_PER_KILOWATT as f64,
                );
                annual_totals_index += 1;
            }
        }
        rows.push(row);
    }

    writer.write_record(&headings_annual)?;
    writer.write_record(&units_annual)?;
    writer.write_record(annual_totals.iter().map(|x| format!("{}", x).into_bytes()))?;
    writer.write_record([""])?;
    writer.write_record(&headings)?;
    writer.write_record(&units_row)?;
    for row in rows {
        writer.write_record(row.iter().map(|x| format!("{}", x).into_bytes()))?;
    }

    Ok(())
}

type HeatSourceWetServiceResultColumn = Vec<(Arc<str>, Option<Arc<str>>)>;

fn write_core_output_file_heat_source_wet(
    output_key: &str,
    timestep_array: &[f64],
    heat_source_wet_results: &ResultsPerTimestep,
    output_writer: &impl OutputWriter,
) -> Result<(), anyhow::Error> {
    // Repeat column headings for each service
    let mut col_headings: Vec<Arc<str>> = vec!["Timestep".into()];
    let mut col_units_row: Vec<Arc<str>> = vec!["count".into()];
    let mut columns: IndexMap<Arc<str>, HeatSourceWetServiceResultColumn> = Default::default();

    for (service_name, service_results) in heat_source_wet_results.iter() {
        columns.insert(
            service_name.clone(),
            service_results.keys().cloned().collect(),
        );
        col_headings.extend(
            service_results
                .keys()
                .cloned()
                .collect::<IndexMap<_, _>>()
                .values()
                .map(|col_heading| match col_heading {
                    None => service_name.clone(),
                    Some(col_heading) => format!("{service_name}: {col_heading}").into(),
                })
                .collect::<Vec<Arc<str>>>(),
        );
        col_units_row.extend(
            service_results
                .keys()
                .cloned()
                .collect::<IndexMap<_, _>>()
                .keys()
                .cloned()
                .collect::<Vec<Arc<str>>>(),
        );
    }

    let writer = output_writer.writer_for_location_key(output_key, "csv")?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    // Write column headings and units
    writer.write_record(
        col_headings
            .iter()
            .map(|x| x.as_bytes())
            .collect::<Vec<_>>(),
    )?;
    writer.write_record(
        col_units_row
            .iter()
            .map(|x| x.as_bytes())
            .collect::<Vec<_>>(),
    )?;

    // Write rows
    for t_idx in 0..timestep_array.len() {
        let mut row: Vec<String> = vec![t_idx.to_string().into()];
        for (service_name, service_results) in heat_source_wet_results {
            row.extend(
                columns[service_name]
                    .iter()
                    .map(|col| service_results[col][t_idx].clone().into()),
            );
        }
        writer.write_record(row.iter().map(|x| x.to_string().into_bytes()))?;
    }

    Ok(())
}

fn write_core_output_file_heat_source_wet_summary(
    output_key: &str,
    heat_source_wet_results_annual: &ResultsAnnual,
    output_writer: &impl OutputWriter,
) -> Result<(), anyhow::Error> {
    let writer = output_writer.writer_for_location_key(output_key, "csv")?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    for (service_name, service_results) in heat_source_wet_results_annual {
        writer.write_record([service_name.to_string()])?;
        for (name, value) in service_results.iter() {
            writer.write_record([
                name.0.as_bytes(),
                name.1.as_ref().map(|x| x.as_bytes()).unwrap_or_default(),
                String::from(value).as_bytes(),
            ])?;
        }
        writer.write_record([""])?;
    }

    Ok(())
}

fn write_core_output_file_emitters_detailed(
    output_prefix: &str,
    emitters_output_map: &IndexMap<Arc<str>, IndexMap<usize, OutputEmitters>>,
    output_writer: &impl OutputWriter,
) -> Result<(), anyhow::Error> {
    for (emitter_name, emitters_detailed_results) in emitters_output_map {
        let output_key = format!("{}{}", output_prefix, emitter_name);
        let writer = output_writer.writer_for_location_key(&output_key, "csv")?;
        let mut writer = WriterBuilder::new()
            .flexible(true)
            .has_headers(false)
            .from_writer(writer);

        writer.write_record([
            "timestep",
            "demand_energy",
            "temp_emitter_req",
            "time_before_heating_start",
            "energy_provided_by_heat_source",
            "temp_emitter",
            "temp_emitter_max",
            "energy_released_from_emitters",
            "temp_flow_target",
            "temp_return_target",
            "temp_emitter_max_is_final_temp",
            "energy_req_from_heat_source",
            "fan_energy_kWh",
        ])?;
        writer.write_record([
            "[count]",
            "[kWh]",
            "[Celsius]",
            "[hours]",
            "[kWh]",
            "[Celsius]",
            "[Celsius]",
            "[kWh]",
            "[Celsius]",
            "[Celsius]",
            "[Boolean]",
            "[kWh]",
            "[kWh]",
        ])?;
        for emitters_detailed_result in emitters_detailed_results.values() {
            writer.serialize(emitters_detailed_result)?;
        }
    }

    Ok(())
}

fn write_core_output_file_esh_detailed(
    output_prefix: &str,
    esh_output: &IndexMap<Arc<str>, IndexMap<usize, Vec<f64>>>,
    output_writer: &impl OutputWriter,
) -> Result<(), anyhow::Error> {
    let headings = [
        "timestep",
        "n_units",
        "demand_energy",
        "energy_delivered",
        "energy_instant",
        "energy_charged",
        "energy_for_fan",
        "state_of_charge",
        "final_soc_ivp",
        "time_used_max",
    ];
    let units_row = [
        "[count]", "[count]", "[kWh]", "[kWh]", "[kWh]", "[kWh]", "[kWh]", "[ratio]", "[ratio]",
        "[hours]",
    ];

    for (esh, esh_output) in esh_output {
        let output_key = format!("{output_prefix}_{esh}");
        let writer = output_writer.writer_for_location_key(&output_key, "csv")?;
        let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);
        writer.write_record(headings)?;
        writer.write_record(units_row)?;
        for esh_results in esh_output.values() {
            writer.serialize(esh_results)?;
        }

        writer.flush()?;
    }

    Ok(())
}

fn write_core_output_file_ventilation_detailed(
    output_key: &str,
    vent_output_list: &[Vec<StringOrNumber>],
    output_writer: &impl OutputWriter,
) -> Result<(), anyhow::Error> {
    let writer = output_writer.writer_for_location_key(output_key, "csv")?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    writer.write_record([
        "Timestep",
        "Incoming air changes per hour",
        "Vent opening ratio",
        "incoming air flow",
        "total_volume",
        "air changes per hour",
        "Internal temperature",
        "Internal reference pressure",
        "Air mass flow rate entering through window opening",
        "Air mass flow rate leaving through window opening",
        "Air mass flow rate entering through vents (openings in the external envelope)",
        "Air mass flow rate leaving through vents (openings in the external envelope)",
        "Air mass flow rate entering through envelope leakage",
        "Air mass flow rate leaving through envelope leakage",
        "Air mass flow rate entering through combustion appliances",
        "Air mass flow rate leaving through combustion appliances",
        "Air mass flow rate entering through passive or hybrid duct",
        "Air mass flow rate leaving through passive or hybrid duct",
        "Supply air mass flow rate going to ventilation zone",
        "Extract air mass flow rate from a ventilation zone",
        "Extract air mass flow rate from heat recovery",
        "Total air mass flow rate entering the zone",
        "Total air mass flow rate leaving the zone",
    ])?;
    writer.write_record([
        "[count]",
        "[indicator]",
        "[ratio]",
        "[m3/h]",
        "[m3]",
        "[ACH]",
        "[Celsius]",
        "[Pa]",
        "[kg/h]",
        "[kg/h]",
        "[kg/h]",
        "[kg/h]",
        "[kg/h]",
        "[kg/h]",
        "[kg/h]",
        "[kg/h]",
        "[kg/h]",
        "[kg/h]",
        "[kg/h]",
        "[kg/h]",
        "[kg/h]",
        "[kg/h]",
        "[kg/h]",
    ])?;

    for ventilation_results in vent_output_list.iter() {
        writer.write_record(ventilation_results.iter().map(StringOrNumber::as_bytes))?;
    }

    Ok(())
}

fn write_core_output_file_hot_water_source_summary(
    output_key: &str,
    hot_water_source_results: &[Vec<StringOrNumber>],
    output_writer: &impl OutputWriter,
) -> anyhow::Result<()> {
    let writer = output_writer.writer_for_location_key(output_key, "csv")?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    for hot_water_source_row in hot_water_source_results {
        writer.write_record(hot_water_source_row.iter().map(StringOrNumber::as_bytes))?;
    }

    Ok(())
}

#[derive(Clone, Debug, PartialEq, Serialize)]
enum StringOrNumber {
    String(String),
    Float(f64),
    Integer(usize),
}

impl StringOrNumber {
    pub(crate) fn as_bytes(&self) -> &[u8] {
        match self {
            StringOrNumber::String(s) => s.as_bytes(),
            StringOrNumber::Float(f) => f.as_bytes(),
            StringOrNumber::Integer(i) => i.as_bytes(),
        }
    }
}

impl From<NumberOrDivisionByZero> for StringOrNumber {
    fn from(value: NumberOrDivisionByZero) -> Self {
        match value {
            NumberOrDivisionByZero::Number(number) => StringOrNumber::Float(number),
            NumberOrDivisionByZero::DivisionByZero => StringOrNumber::String("DIV/0".into()),
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
                StringOrNumber::Float(number) => Cow::Owned(number.to_string()),
                StringOrNumber::Integer(number) => Cow::Owned(number.to_string()),
            }
        )
    }
}

impl From<StringOrNumber> for Vec<u8> {
    fn from(value: StringOrNumber) -> Self {
        format!("{}", value).into_bytes()
    }
}

impl From<StringOrNumber> for String {
    fn from(value: StringOrNumber) -> Self {
        format!("{}", value).into()
    }
}

impl From<StringOrNumber> for Arc<str> {
    fn from(value: StringOrNumber) -> Self {
        format!("{}", value).into()
    }
}

impl From<&str> for StringOrNumber {
    fn from(value: &str) -> Self {
        StringOrNumber::String(value.into())
    }
}

impl From<String> for StringOrNumber {
    fn from(value: String) -> Self {
        StringOrNumber::String(value.clone())
    }
}

impl From<std::string::String> for StringOrNumber {
    fn from(value: std::string::String) -> Self {
        StringOrNumber::String(value.into())
    }
}

impl From<&f64> for StringOrNumber {
    fn from(value: &f64) -> Self {
        StringOrNumber::Float(*value)
    }
}

impl From<f64> for StringOrNumber {
    fn from(value: f64) -> Self {
        StringOrNumber::Float(value)
    }
}

impl From<usize> for StringOrNumber {
    fn from(value: usize) -> Self {
        StringOrNumber::Integer(value)
    }
}

impl From<u32> for StringOrNumber {
    fn from(value: u32) -> Self {
        StringOrNumber::Integer(value as usize)
    }
}

impl From<StringOrNumber> for f64 {
    fn from(value: StringOrNumber) -> Self {
        match value {
            StringOrNumber::Float(number) => number,
            StringOrNumber::Integer(number) => number as f64,
            StringOrNumber::String(_) => 0.,
        }
    }
}

impl From<&StringOrNumber> for f64 {
    fn from(value: &StringOrNumber) -> Self {
        match value {
            StringOrNumber::Float(number) => *number,
            StringOrNumber::Integer(number) => *number as f64,
            StringOrNumber::String(_) => 0.,
        }
    }
}

impl AddAssign for StringOrNumber {
    fn add_assign(&mut self, rhs: Self) {
        match (self, rhs) {
            (StringOrNumber::Float(a), StringOrNumber::Float(b)) => *a += b,
            (StringOrNumber::Float(a), StringOrNumber::Integer(b)) => *a += b as f64,
            _ => panic!(
                "Cannot add to a non-float StringOrNumber, and cannot add strings to floats."
            ),
        }
    }
}

impl From<&ExternalConditionsFromFile> for ExternalConditionsInput {
    fn from(value: &ExternalConditionsFromFile) -> Self {
        Self {
            air_temperatures: Some(value.air_temperatures.clone()),
            wind_speeds: Some(value.wind_speeds.clone()),
            wind_directions: Some(value.wind_directions.clone()),
            diffuse_horizontal_radiation: Some(value.diffuse_horizontal_radiation.clone()),
            direct_beam_radiation: Some(value.direct_beam_radiation.clone()),
            solar_reflectivity_of_ground: Some(value.solar_reflectivity_of_ground.clone()),
            latitude: Some(value.latitude),
            longitude: Some(value.longitude),
            direct_beam_conversion_needed: Some(value.direct_beam_conversion_needed),
            ..Default::default()
        }
    }
}
