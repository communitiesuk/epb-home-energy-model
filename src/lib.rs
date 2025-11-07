#![allow(clippy::too_many_arguments, clippy::doc_overindented_list_items)]

mod compare_floats;
pub mod core;
pub mod corpus;
pub mod errors;
mod hem_core;
pub mod input;
pub mod output;
pub mod read_weather_file;
mod statistics;
mod wrappers;

#[macro_use]
extern crate is_close;

use crate::core::heating_systems::elec_storage_heater::StorageHeaterDetailedResult;
use crate::core::heating_systems::emitters::EmittersDetailedResult;
use crate::core::heating_systems::storage_tank::StorageTankDetailedResult;
use crate::core::space_heat_demand::ventilation::VentilationDetailedResult;
use crate::core::units::{convert_profile_to_daily, WATTS_PER_KILOWATT};
pub use crate::corpus::RunResults;
use crate::corpus::{
    calc_htc_hlp, Corpus, HeatingCoolingSystemResultKey, HotWaterResultKey, HotWaterResultMap,
    HtcHlpCalculation, NumberOrDivisionByZero, ResultsAnnual, ResultsEndUser, ResultsPerTimestep,
    ZoneResultKey,
};
use crate::errors::{HemCoreError, HemError, NotImplementedError, PostprocessingError};
use crate::external_conditions::ExternalConditions;
use crate::input::{
    ingest_for_processing, ExternalConditionsInput, FuelType, HotWaterSourceDetails, Input,
    InputForProcessing, SchemaReference,
};
use crate::output::Output;
use crate::read_weather_file::ExternalConditions as ExternalConditionsFromFile;
use crate::simulation_time::SimulationTime;
use crate::statistics::percentile;
#[cfg(feature = "fhs")]
use crate::wrappers::future_homes_standard::{FhsComplianceWrapper, FhsSingleCalcWrapper};
pub use crate::wrappers::HemResponse;
use crate::wrappers::{ChosenWrapper, HemWrapper, PassthroughHemWrapper};
use anyhow::anyhow;
use bitflags::bitflags;
use chrono::prelude::*;
use chrono::{TimeDelta, Utc};
use convert_case::{Case, Casing};
use csv::WriterBuilder;
use hem_core::external_conditions;
use hem_core::simulation_time;
use indexmap::IndexMap;
use itertools::Itertools;
use rayon::prelude::*;
use smartstring::alias::String;
use std::borrow::Cow;
use std::collections::HashMap;
use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::io::Read;
use std::ops::AddAssign;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::{Arc, LazyLock};
use tracing::{debug, error, instrument};

pub const HEM_VERSION: &str = "1.0.0a1";
pub const HEM_VERSION_DATE: &str = "2025-10-02";

#[instrument(skip_all)]
pub fn run_project(
    input: impl Read,
    output: impl Output,
    external_conditions_data: Option<ExternalConditionsFromFile>,
    tariff_data_file: Option<&str>,
    flags: &ProjectFlags,
) -> Result<Option<HemResponse>, HemError> {
    catch_unwind(AssertUnwindSafe(|| {
        // 1. ingest/ parse input and enter preprocessing stage
        #[instrument(skip_all)]
        fn ingest_input_and_start_preprocessing(
            input: impl Read,
            external_conditions_data: Option<&ExternalConditionsFromFile>,
            schema_reference: &SchemaReference,
        ) -> anyhow::Result<InputForProcessing> {
            let mut input_for_processing = ingest_for_processing(input, schema_reference)?;

            input_for_processing
                .merge_external_conditions_data(external_conditions_data.map(|x| x.into()))?;
            Ok(input_for_processing)
        }

        fn choose_schema_reference(flags: &ProjectFlags) -> SchemaReference {
            let mut schema_reference = SchemaReference::Core;
            #[cfg(feature = "fhs")]
            if flags.intersects(ProjectFlags::FHS_ASSUMPTIONS
                | ProjectFlags::FHS_FEE_ASSUMPTIONS
                | ProjectFlags::FHS_NOT_A_ASSUMPTIONS
                | ProjectFlags::FHS_NOT_B_ASSUMPTIONS
                | ProjectFlags::FHS_FEE_NOT_A_ASSUMPTIONS
                | ProjectFlags::FHS_FEE_NOT_B_ASSUMPTIONS | ProjectFlags::FHS_COMPLIANCE) {
                schema_reference = SchemaReference::Fhs;
            }

            schema_reference
        }

        let input_for_processing =
            ingest_input_and_start_preprocessing(input, external_conditions_data.as_ref(), &choose_schema_reference(flags))?;

        fn choose_wrapper(flags: &ProjectFlags) -> ChosenWrapper {
            #[cfg(feature = "fhs")]
            {
                if flags.contains(ProjectFlags::FHS_COMPLIANCE) {
                    ChosenWrapper::FhsCompliance(FhsComplianceWrapper::new())
                } else if flags.intersects(
                    ProjectFlags::FHS_ASSUMPTIONS
                        | ProjectFlags::FHS_FEE_ASSUMPTIONS
                        | ProjectFlags::FHS_NOT_A_ASSUMPTIONS
                        | ProjectFlags::FHS_NOT_B_ASSUMPTIONS
                        | ProjectFlags::FHS_FEE_NOT_A_ASSUMPTIONS
                        | ProjectFlags::FHS_FEE_NOT_B_ASSUMPTIONS,
                ) {
                    ChosenWrapper::FhsSingleCalc(FhsSingleCalcWrapper::new())
                } else {
                    ChosenWrapper::Passthrough(PassthroughHemWrapper::new())
                }
            }
            #[cfg(not(feature = "fhs"))]
            {
                ChosenWrapper::Passthrough(PassthroughHemWrapper::new())
            }
        }

        let wrapper = choose_wrapper(flags);

        // 2. apply preprocessing from wrappers
        #[instrument(skip_all)]
        fn apply_preprocessing_from_wrappers(
            input_for_processing: InputForProcessing,
            wrapper: &impl HemWrapper,
            flags: &ProjectFlags,
        ) -> anyhow::Result<HashMap<CalculationKey, Input>> {
            wrapper.apply_preprocessing(input_for_processing, flags)
        }

        let input = match catch_unwind(AssertUnwindSafe(|| {
            apply_preprocessing_from_wrappers(input_for_processing, &wrapper, flags)
                .map_err(HemError::InvalidRequest)
        })) {
            Ok(result) => result?,
            Err(panic) => {
                return Err(HemError::PanicInWrapper(
                    panic
                        .downcast_ref::<&str>()
                        .map_or("Error not captured", |v| v)
                        .to_owned(),
                ))
            }
        };

        let cloned_input = input.get(&CalculationKey::Primary).cloned();

        // 2b.(!) If preprocess-only flag is present and there is a primary calculation key, write out preprocess file
        if flags.contains(ProjectFlags::PRE_PROCESS_ONLY) {
            if let Some(input) = input.get(&CalculationKey::Primary) {
                write_preproc_file(input, &output, "preproc", "json")?;
            } else {
                error!("Preprocess-only flag only set up to work with a calculation using a primary calculation key (i.e. not FHS compliance)");
            }

            return Ok(None);
        }

        #[instrument(skip_all)]
        fn write_preproc_file(input: &Input, output: &impl Output, location_key: &str, file_extension: &str) -> anyhow::Result<()> {
            let writer = output.writer_for_location_key(location_key, file_extension)?;
            if let Err(e) = serde_json::to_writer_pretty(writer, input) {
                error!("Could not write out preprocess file: {}", e);
            }

            Ok(())
        }


        // 3. Determine external conditions to use for calculations.
        #[instrument(skip_all)]
        fn resolve_external_conditions(
            input: &HashMap<CalculationKey, Input>,
            external_conditions_data: Option<ExternalConditionsFromFile>,
        ) -> HashMap<CalculationKey, ExternalConditions> {
            input
                .par_iter()
                .map(|(key, input)| {
                    (
                        *key,
                        external_conditions_from_input(
                            input.external_conditions.clone(),
                            external_conditions_data.clone(),
                            input.simulation_time,
                        ),
                    )
                })
                .collect()
        }

        let corpora = {
            let external_conditions = resolve_external_conditions(&input, external_conditions_data);

            // 4. Build corpus from input and external conditions.
            #[instrument(skip_all)]
            fn build_corpus(
                input: &HashMap<CalculationKey, Input>,
                external_conditions: &HashMap<CalculationKey, ExternalConditions>,
                tariff_data_file: Option<&str>,
                flags: &ProjectFlags,
            ) -> anyhow::Result<HashMap<CalculationKey, Corpus>> {
                let output_options = flags.into();
                iterate_maps(input, external_conditions)
                    .map(|(key, input, external_conditions)| {
                        anyhow::Ok((*key, Corpus::from_inputs(input, Some(external_conditions), tariff_data_file, &output_options)?))
                    })
                    .collect()
            }

            build_corpus(&input, &external_conditions, tariff_data_file, flags).map_err(|e| {
                capture_specific_error_case(&e).unwrap_or_else(|| HemError::InvalidRequest(e))
            })?
        };

        // 5. Run HEM calculation(s).
        #[instrument(skip_all)]
        fn run_hem_calculation(
            corpora: &HashMap<CalculationKey, Corpus>,
        ) -> anyhow::Result<HashMap<CalculationKey, RunResults>> {
            corpora
                .par_iter()
                .map(|(key, corpus)| anyhow::Ok((*key, corpus.run()?)))
                .collect()
        }

        // catch_unwind here catches any downstream panics so we can at least map to the right HemError variant
        let run_results = match catch_unwind(AssertUnwindSafe(|| {
            run_hem_calculation(&corpora).map_err(|e| {
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

        let contextualised_results: HashMap<_, _> = run_results
            .iter()
            .map(|(key, results)| {
                (
                    *key,
                    CalculationResultsWithContext::new(&input[key], &corpora[key], results),
                )
            })
            .collect();

        // 6. Write out to core output files.
        #[instrument(skip_all)]
        fn write_core_output_files(
            primary_input: Option<&Input>,
            output: &impl Output,
            results: &HashMap<CalculationKey, CalculationResultsWithContext>,
            corpora: &HashMap<CalculationKey, Corpus>,
            flags: &ProjectFlags,
        ) -> anyhow::Result<()> {
            if let Some(results) = results.get(&CalculationKey::Primary) {
                if output.is_noop() {
                    return Ok(());
                }
                let CalculationResultsWithContext {
                    results:
                    RunResults {
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
                        heat_balance_dict,
                        heat_source_wet_results_dict,
                        heat_source_wet_results_annual_dict,
                        vent_output_list,
                        emitters_output_dict,
                        storage_from_grid,
                        battery_state_of_charge,
                        esh_output_dict,
                        hot_water_source_results_dict,
                        ..
                    },
                    ..
                } = results;
                write_core_output_file(
                    output,
                    OutputFileArgs {
                        output_key: "results".into(),
                        timestep_array,
                        results_totals,
                        results_end_user,
                        energy_import,
                        energy_export,
                        energy_generated_consumed,
                        energy_to_storage,
                        energy_from_storage,
                        storage_from_grid,
                        battery_state_of_charge,
                        energy_diverted,
                        betafactor,
                        zone_dict,
                        zone_list,
                        hc_system_dict,
                        hot_water_dict,
                        ductwork_gains,
                    },
                )?;

                if flags.contains(ProjectFlags::HEAT_BALANCE) {
                    let hour_per_step = corpora[&CalculationKey::Primary].simulation_time.step_in_hours();
                    for (hb_name, hb_map) in heat_balance_dict.iter() {
                        let output_key = format!("results_heat_balance_{}", hb_name.to_string().to_case(Case::Snake)).into();
                        write_core_output_file_heat_balance(output, HeatBalanceOutputFileArgs {
                            output_key,
                            timestep_array,
                            hour_per_step,
                            heat_balance_map: hb_map,
                        })?;
                    }
                }

                if flags.contains(ProjectFlags::DETAILED_OUTPUT_HEATING_COOLING) {
                    for (heat_source_wet_name, heat_source_wet_results) in heat_source_wet_results_dict.iter() {
                        let output_key = format!("results_heat_source_wet__{heat_source_wet_name}");
                        write_core_output_file_heat_source_wet(output, &output_key, timestep_array, heat_source_wet_results)?;
                    }
                    for (heat_source_wet_name, heat_source_wet_results_annual) in heat_source_wet_results_annual_dict.iter() {
                        let output_key = format!("results_heat_source_wet_summary__{heat_source_wet_name}");
                        write_core_output_file_heat_source_wet_summary(output, &output_key, heat_source_wet_results_annual)?;
                    }
                    // Function call to write detailed ventilation results
                    let vent_output_file = "ventilation_results";
                    write_core_output_file_ventilation_detailed(output, vent_output_file, vent_output_list)?;
                    for (hot_water_source_name, hot_water_source_results) in hot_water_source_results_dict.iter() {
                        let hot_water_source_file = format!("results_hot_water_source_summary__{}", hot_water_source_name.replace(" ", "_"));
                        write_core_output_file_hot_water_source_summary(output, &hot_water_source_file, hot_water_source_results);
                    }
                }

                write_core_output_file_summary(output, results.try_into()?)?;

                let corpus = results.context.corpus;

                let primary_input = primary_input.ok_or_else(|| anyhow!("Primary input should be available as there is a primary calculation."))?;

                let HtcHlpCalculation { total_htc: heat_transfer_coefficient, total_hlp: heat_loss_parameter, .. } = calc_htc_hlp(primary_input)?;
                let heat_capacity_parameter = corpus.calc_hcp();
                let heat_loss_form_factor = corpus.calc_hlff();

                write_core_output_file_static(
                    output,
                    StaticOutputFileArgs {
                        output_key: "results_static".into(),
                        heat_transfer_coefficient,
                        heat_loss_parameter,
                        heat_capacity_parameter,
                        heat_loss_form_factor,
                        temp_internal_air: primary_input.temp_internal_air_static_calcs,
                        temp_external_air: corpus.external_conditions.air_temp_annual_daily_average_min(),
                    },
                )?;

                if flags.contains(ProjectFlags::DETAILED_OUTPUT_HEATING_COOLING) {
                    let output_prefix = "results_emitters_";
                    write_core_output_file_emitters_detailed(output, output_prefix, emitters_output_dict)?;

                    let esh_output_prefix = "results_esh_";
                    write_core_output_file_esh_detailed(output, esh_output_prefix, esh_output_dict)?;
                }
            }

            Ok(())
        }

        write_core_output_files(cloned_input.as_ref(), &output, &contextualised_results, &corpora, flags)?;

        // 7. Run wrapper post-processing and capture any output.
        #[instrument(skip_all)]
        fn run_wrapper_postprocessing(
            output: &impl Output,
            results: &HashMap<CalculationKey, CalculationResultsWithContext>,
            wrapper: &impl HemWrapper,
            flags: &ProjectFlags,
        ) -> anyhow::Result<Option<HemResponse>> {
            wrapper.apply_postprocessing(output, results, flags)
        }

        run_wrapper_postprocessing(&output, &contextualised_results, &wrapper, flags)
            .map_err(|e| HemError::ErrorInPostprocessing(PostprocessingError::new(e)))
    }))
        .map_err(|e| {
            HemError::GeneralPanic(
                e.downcast_ref::<&str>()
                    .map_or("Uncaught panic - could not capture more information; to debug further, try replicating in debug mode.", |v| v)
                    .to_owned(),
            )
        })?
}

fn capture_specific_error_case(e: &anyhow::Error) -> Option<HemError> {
    if let Some(e) = e.downcast_ref::<NotImplementedError>() {
        return Some(HemError::NotImplemented(e.clone()));
    }

    None
}

#[derive(Clone, Copy)]
pub(crate) struct CalculationContext<'a> {
    input: &'a Input,
    corpus: &'a Corpus,
}

#[derive(Clone, Copy)]
pub(crate) struct CalculationResultsWithContext<'a> {
    results: &'a RunResults,
    context: CalculationContext<'a>,
}

impl<'a> CalculationResultsWithContext<'a> {
    fn new(
        input: &'a Input,
        corpus: &'a Corpus,
        results: &'a RunResults,
    ) -> CalculationResultsWithContext<'a> {
        Self {
            results,
            context: CalculationContext { input, corpus },
        }
    }
}

impl CalculationResultsWithContext<'_> {
    fn daily_hw_demand_percentile(&self, percentage: usize) -> anyhow::Result<f64> {
        Ok(percentile(
            &convert_profile_to_daily(
                match &self.results.hot_water_dict
                    [&HotWaterResultKey::HotWaterEnergyDemandIncludingPipeworkLoss]
                {
                    HotWaterResultMap::Float(results) => results
                        .get("energy_demand_incl_pipework_loss")
                        .ok_or(anyhow!(
                    "Hot water energy demand incl pipework_loss field not set in hot water output"
                ))?,
                    HotWaterResultMap::Int(_) => unreachable!(
                        "Hot water energy demand incl pipework_loss is not expected to be an integer"
                    ),
                },
                self.context.corpus.simulation_time.step_in_hours(),
            ),
            percentage,
        ))
    }
}

bitflags! {
    pub struct ProjectFlags: u32 {
        const PRE_PROCESS_ONLY = 0b1;
        const HEAT_BALANCE = 0b10;
        const DETAILED_OUTPUT_HEATING_COOLING = 0b100;
        // start FHS flags from 2^8
        #[cfg(feature = "fhs")]
        const FHS_ASSUMPTIONS = 0b100000000;
        #[cfg(feature = "fhs")]
        const FHS_FEE_ASSUMPTIONS = 0b1000000000;
        #[cfg(feature = "fhs")]
        const FHS_NOT_A_ASSUMPTIONS = 0b10000000000;
        #[cfg(feature = "fhs")]
        const FHS_NOT_B_ASSUMPTIONS = 0b100000000000;
        #[cfg(feature = "fhs")]
        const FHS_FEE_NOT_A_ASSUMPTIONS = 0b1000000000000;
        #[cfg(feature = "fhs")]
        const FHS_FEE_NOT_B_ASSUMPTIONS = 0b10000000000000;
        #[cfg(feature = "fhs")]
        const FHS_COMPLIANCE = 0b100000000000000;
    }
}

#[derive(Copy, Clone, Debug, Eq, Hash, PartialEq)]
pub(crate) enum CalculationKey {
    Primary,
    #[cfg(feature = "fhs")]
    Fhs,
    #[cfg(feature = "fhs")]
    FhsFee,
    #[cfg(feature = "fhs")]
    FhsNotional,
    #[cfg(feature = "fhs")]
    FhsNotionalFee,
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
            None,
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
            None,
            false,
            input.direct_beam_conversion_needed.unwrap_or(false),
            input.shading_segments.clone(),
        ),
    }
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
            "DHW: demand volume (including distribution pipework losses)",
            "[litres]",
        ),
        (
            "DHW: demand energy (including distribution pipework losses)",
            "[kWh]",
        ),
        (
            "DHW: demand energy (excluding distribution pipework losses)",
            "[kWh]",
        ),
        ("DHW: total event duration", "[mins]"),
        ("DHW: number of events", "[count]"),
        ("DHW: distribution pipework losses", "[kWh]"),
        ("DHW: primary pipework losses", "[kWh]"),
        ("DHW: storage losses", "[kWh]"),
    ])
});

struct OutputFileArgs<'a> {
    output_key: String,
    timestep_array: &'a [f64],
    results_totals: &'a IndexMap<String, Vec<f64>>,
    results_end_user: &'a IndexMap<String, IndexMap<String, Vec<f64>>>,
    energy_import: &'a IndexMap<String, Vec<f64>>,
    energy_export: &'a IndexMap<String, Vec<f64>>,
    energy_generated_consumed: &'a IndexMap<String, Vec<f64>>,
    energy_to_storage: &'a IndexMap<String, Vec<f64>>,
    energy_from_storage: &'a IndexMap<String, Vec<f64>>,
    storage_from_grid: &'a IndexMap<String, Vec<f64>>,
    battery_state_of_charge: &'a IndexMap<String, Vec<f64>>,
    energy_diverted: &'a IndexMap<String, Vec<f64>>,
    betafactor: &'a IndexMap<String, Vec<f64>>,
    zone_dict: &'a IndexMap<ZoneResultKey, IndexMap<String, Vec<f64>>>,
    zone_list: &'a [String],
    hc_system_dict: &'a IndexMap<HeatingCoolingSystemResultKey, IndexMap<String, Vec<f64>>>,
    hot_water_dict: &'a IndexMap<HotWaterResultKey, HotWaterResultMap>,
    ductwork_gains: &'a IndexMap<String, Vec<f64>>,
}

fn write_core_output_file(output: &impl Output, args: OutputFileArgs) -> anyhow::Result<()> {
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
        storage_from_grid,
        battery_state_of_charge,
        energy_diverted,
        betafactor,
        zone_dict,
        zone_list,
        hc_system_dict,
        hot_water_dict,
        ductwork_gains,
    } = args;
    debug!("writing out to {output_key}");
    let writer = output.writer_for_location_key(&output_key, "csv")?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    let mut headings: Vec<Cow<'static, str>> = vec!["Timestep".into()];
    let mut units_row = vec!["[count]"];

    // hot_water_dict headings
    for system in hot_water_dict.keys() {
        let system = match system {
            HotWaterResultKey::HotWaterDemand => {
                "DHW: demand volume (including distribution pipework losses)"
            }
            HotWaterResultKey::HotWaterEnergyDemand => {
                "DHW: demand energy (excluding distribution pipework losses)"
            }
            HotWaterResultKey::HotWaterEnergyDemandIncludingPipeworkLoss => {
                "DHW: demand energy (including distribution pipework losses)"
            }
            HotWaterResultKey::HotWaterDuration => "DHW: total event duration",
            HotWaterResultKey::HotWaterEvents => "DHW: number of events",
            HotWaterResultKey::PipeworkLosses => "DHW: distribution pipework losses",
            HotWaterResultKey::PrimaryPipeworkLosses => "DHW: primary pipework losses",
            HotWaterResultKey::StorageLosses => "DHW: storage losses",
        };
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
            if UNITS_MAP.contains_key(zone_outputs.as_str()) {
                units_row.push(UNITS_MAP[zone_outputs.as_str()]);
            } else {
                units_row.push("Unit not defined");
            }
        }
    }

    // hc_system_dict holds heating demand and output as first level keys
    // and the system name as second level keys.
    // Reorganising this dictionary so system names can be grouped together

    // Initialize the reorganized dictionary for grouping systems in hc_system_dict
    let mut reorganised_dict: IndexMap<String, IndexMap<HeatingCoolingSystemResultKey, Vec<f64>>> =
        Default::default();

    // Iterate over the original map
    for (key, value) in hc_system_dict {
        // Iterate over the nested map
        for (nested_key, nested_value) in value {
            // Add the nested_value to the corresponding entry in reorganized_dict
            reorganised_dict
                .entry(nested_key.clone())
                .or_default()
                .insert(*key, nested_value.clone());
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
            let hc_system = if matches!(
                hc_name,
                HeatingCoolingSystemResultKey::HeatingSystem
                    | HeatingCoolingSystemResultKey::CoolingSystem
            ) {
                let alternate_name = "energy demand";
                format!("{system_label}: {alternate_name}")
            } else if matches!(
                hc_name,
                HeatingCoolingSystemResultKey::HeatingSystemOutput
                    | HeatingCoolingSystemResultKey::CoolingSystemOutput
            ) {
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
            energy_use_row.push(storage_from_grid[totals_key][t_idx]);
            energy_use_row.push(battery_state_of_charge[totals_key][t_idx]);
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
        if let HotWaterResultMap::Float(map) = &hot_water_dict[&HotWaterResultKey::HotWaterDemand] {
            hw_system_row.push(map["demand"][t_idx]);
        }
        if let HotWaterResultMap::Float(map) =
            &hot_water_dict[&HotWaterResultKey::HotWaterEnergyDemand]
        {
            hw_system_row_energy.push(map["energy_demand"][t_idx]);
        }
        if let HotWaterResultMap::Float(map) =
            &hot_water_dict[&HotWaterResultKey::HotWaterEnergyDemandIncludingPipeworkLoss]
        {
            hw_system_row_energy_with_pipework_losses
                .push(map["energy_demand_incl_pipework_loss"][t_idx]);
        }
        if let HotWaterResultMap::Float(map) = &hot_water_dict[&HotWaterResultKey::HotWaterDuration]
        {
            hw_system_row_duration.push(map["duration"][t_idx]);
        }
        if let HotWaterResultMap::Float(map) = &hot_water_dict[&HotWaterResultKey::PipeworkLosses] {
            pw_losses_row.push(map["pw_losses"][t_idx]);
        }
        if let HotWaterResultMap::Float(map) =
            &hot_water_dict[&HotWaterResultKey::PrimaryPipeworkLosses]
        {
            primary_pw_losses_row.push(map["primary_pw_losses"][t_idx]);
        }
        if let HotWaterResultMap::Float(map) = &hot_water_dict[&HotWaterResultKey::StorageLosses] {
            storage_losses_row.push(map["storage_losses"][t_idx]);
        }
        if let HotWaterResultMap::Int(map) = &hot_water_dict[&HotWaterResultKey::HotWaterEvents] {
            hw_system_row_events.push(map["no_events"][t_idx]);
        }
        ductwork_row.push(ductwork_gains["ductwork_gains"][t_idx]);

        // create row of outputs and write to output file
        let mut row: Vec<String> = vec![];
        row.append(&mut vec![t_idx.to_string().into()]);
        row.append(
            &mut hw_system_row
                .into_iter()
                .map(|x| x.to_string().into())
                .collect(),
        );
        row.append(
            &mut hw_system_row_energy_with_pipework_losses
                .into_iter()
                .map(|x| x.to_string().into())
                .collect(),
        );
        row.append(
            &mut hw_system_row_energy
                .into_iter()
                .map(|x| x.to_string().into())
                .collect(),
        );
        row.append(
            &mut hw_system_row_duration
                .into_iter()
                .map(|x| x.to_string().into())
                .collect(),
        );
        row.append(
            &mut hw_system_row_events
                .into_iter()
                .map(|x| x.to_string().into())
                .collect(),
        );
        row.append(
            &mut pw_losses_row
                .into_iter()
                .map(|x| x.to_string().into())
                .collect(),
        );
        row.append(
            &mut primary_pw_losses_row
                .into_iter()
                .map(|x| x.to_string().into())
                .collect(),
        );
        row.append(
            &mut storage_losses_row
                .into_iter()
                .map(|x| x.to_string().into())
                .collect(),
        );
        row.append(
            &mut ductwork_row
                .into_iter()
                .map(|x| x.to_string().into())
                .collect(),
        );
        row.append(
            &mut energy_shortfall
                .into_iter()
                .map(|x| x.to_string().into())
                .collect(),
        );
        row.append(&mut zone_row.into_iter().map(|x| x.to_string().into()).collect());
        row.append(
            &mut hc_system_row
                .into_iter()
                .map(|x| x.to_string().into())
                .collect(),
        );
        row.append(
            &mut energy_use_row
                .into_iter()
                .map(|x| x.to_string().into())
                .collect(),
        );

        writer.write_record(&row)?;
    }

    debug!("flushing out CSV");
    writer.flush()?;

    Ok(())
}

struct SummaryOutputFileArgs<'a> {
    output_key: String,
    input: SummaryInputDigest,
    timestep_array: &'a [f64],
    results_end_user: &'a IndexMap<String, IndexMap<String, Vec<f64>>>,
    energy_generated_consumed: &'a IndexMap<String, Vec<f64>>,
    energy_to_storage: &'a IndexMap<String, Vec<f64>>,
    energy_from_storage: &'a IndexMap<String, Vec<f64>>,
    energy_diverted: &'a IndexMap<String, Vec<f64>>,
    energy_import: &'a IndexMap<String, Vec<f64>>,
    energy_export: &'a IndexMap<String, Vec<f64>>,
    storage_from_grid: &'a IndexMap<String, Vec<f64>>,
    space_heat_demand_total: f64,
    space_cool_demand_total: f64,
    total_floor_area: f64,
    heat_cop_dict: &'a IndexMap<String, NumberOrDivisionByZero>,
    cool_cop_dict: &'a IndexMap<String, NumberOrDivisionByZero>,
    dhw_cop_dict: &'a IndexMap<String, NumberOrDivisionByZero>,
    daily_hw_demand_75th_percentile: f64,
}

/// A digest of data from the input that is used when generating a summary file.
#[derive(Clone)]
struct SummaryInputDigest {
    simulation_time: SimulationTime,
    hot_water_source_digests: IndexMap<String, SummaryInputHotWaterSourceDigest>,
    electricity_keys: Vec<String>,
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
            electricity_keys: input
                .energy_supply
                .iter()
                .filter(|&(_, energy_supply_details)| {
                    energy_supply_details.fuel == FuelType::Electricity
                })
                .map(|(key, _)| key.into())
                .collect(),
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

impl<'a> TryFrom<&CalculationResultsWithContext<'a>> for SummaryOutputFileArgs<'a> {
    type Error = anyhow::Error;

    fn try_from(value: &CalculationResultsWithContext<'a>) -> Result<Self, Self::Error> {
        let RunResults {
            timestep_array,
            results_end_user,
            energy_import,
            energy_export,
            energy_generated_consumed,
            energy_to_storage,
            energy_from_storage,
            storage_from_grid,
            energy_diverted,
            heat_cop_dict,
            cool_cop_dict,
            dhw_cop_dict,
            ..
        } = value.results;
        Ok(SummaryOutputFileArgs {
            output_key: "results_summary".into(),
            input: value.context.input.into(),
            timestep_array,
            results_end_user,
            energy_generated_consumed,
            energy_to_storage,
            energy_from_storage,
            energy_diverted,
            energy_import,
            energy_export,
            storage_from_grid,
            space_heat_demand_total: value.results.space_heat_demand_total(),
            space_cool_demand_total: value.results.space_cool_demand_total(),
            total_floor_area: value.context.corpus.total_floor_area,
            heat_cop_dict,
            cool_cop_dict,
            dhw_cop_dict,
            daily_hw_demand_75th_percentile: value.daily_hw_demand_percentile(75)?,
        })
    }
}

impl<'a> TryFrom<&CalculationResultsWithContext<'a>> for SummaryDataArgs<'a> {
    type Error = anyhow::Error;

    fn try_from(results: &CalculationResultsWithContext<'a>) -> Result<Self, Self::Error> {
        let RunResults {
            timestep_array,
            results_end_user,
            energy_import,
            energy_export,
            energy_generated_consumed,
            energy_to_storage,
            energy_from_storage,
            storage_from_grid,
            energy_diverted,
            ..
        } = results.results;
        Ok(SummaryDataArgs {
            timestep_array,
            input: results.context.input.into(),
            results_end_user,
            energy_generated_consumed,
            energy_to_storage,
            energy_from_storage,
            storage_from_grid,
            energy_diverted,
            energy_import,
            energy_export,
        })
    }
}

fn write_core_output_file_summary(
    output: &impl Output,
    args: SummaryOutputFileArgs,
) -> Result<(), anyhow::Error> {
    if output.is_noop() {
        return Ok(());
    }
    let SummaryOutputFileArgs {
        output_key,
        input,
        timestep_array,
        results_end_user,
        energy_generated_consumed,
        energy_to_storage,
        energy_from_storage,
        storage_from_grid,
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

    debug!("writing out to {output_key}");

    let SummaryData {
        delivered_energy_map,
        stats,
        peak_elec_consumption,
        index_peak_elec_consumption,
        step_peak_elec_consumption,
        timestep_to_date,
    } = build_summary_data(SummaryDataArgs {
        timestep_array,
        input: input.clone(),
        results_end_user,
        energy_generated_consumed,
        energy_to_storage,
        energy_from_storage,
        storage_from_grid,
        energy_diverted,
        energy_import,
        energy_export,
    });

    let mut delivered_energy_rows_title =
        vec!["Delivered energy by end-use (below) and fuel (right) [kWh/m2]".into()];
    let mut delivered_energy_rows = vec![vec![StringOrNumber::String("total".into())]];
    for (fuel, end_uses) in delivered_energy_map {
        delivered_energy_rows_title.push(fuel.clone());
        for row in delivered_energy_rows.iter_mut() {
            row.push(StringOrNumber::Float(0.));
        }
        for (end_use, value) in end_uses {
            let mut end_use_found = false;
            for row in delivered_energy_rows.iter_mut() {
                if row.contains(&StringOrNumber::String(end_use.clone())) {
                    end_use_found = true;
                    *row.get_mut(
                        delivered_energy_rows_title
                            .iter()
                            .position(|x| x == &fuel)
                            .unwrap(),
                    )
                    .unwrap() = StringOrNumber::Float(value / total_floor_area);
                }
            }
            if !end_use_found {
                let mut new_row =
                    vec![StringOrNumber::Float(0.); delivered_energy_rows_title.len()];
                *new_row.get_mut(0).unwrap() = StringOrNumber::String(end_use);
                *new_row
                    .get_mut(
                        delivered_energy_rows_title
                            .iter()
                            .position(|x| x == &fuel)
                            .unwrap(),
                    )
                    .unwrap() = StringOrNumber::Float(value / total_floor_area);
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
        .map(|(hw_name, hw_cop)| vec![StringOrNumber::from(hw_name.as_str()), (*hw_cop).into()])
        .collect::<Vec<_>>();

    let writer = output.writer_for_location_key(&output_key, "csv")?;
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
    writer.write_record(["Energy Supply Summary"])?;
    writer.write_record(["", "kWh", "timestep", "month", "day", "hour of day"])?;
    writer.write_record([
        "Peak half-hour consumption (electricity)".to_string(),
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
    writer.write_record(&blank_line)?;
    let mut header_row = vec!["".to_owned(), "Total".to_owned()];
    header_row.extend(stats.keys().map(|x| x.to_string()));
    writer.write_record(&header_row)?;
    let fields = [
        ("Consumption", "kWh", EnergySupplyStatKey::ElecConsumed),
        ("Generation", "kWh", EnergySupplyStatKey::ElecGenerated),
        (
            "Generation to consumption (immediate excl. diverter)",
            "kWh",
            EnergySupplyStatKey::GenToConsumption,
        ),
        (
            "Generation to storage",
            "kWh",
            EnergySupplyStatKey::GenToStorage,
        ),
        (
            "Generation to diverter",
            "kWh",
            EnergySupplyStatKey::GenToDiverter,
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
            EnergySupplyStatKey::StorageEff,
        ),
    ];
    for field in fields {
        let mut row: Vec<String> = vec![field.0.into(), field.1.into()];
        for stat in stats.values() {
            row.push(stat.display_for_key(&field.2));
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
            let first_item_as_string = String::from(row[0].clone());
            row.push(StringOrNumber::Float(daily_hw_demand_75th_percentile));
            row.push({
                let hot_water_source: SummaryInputHotWaterSourceDigest = *input
                    .hot_water_source_digests
                    .get(&first_item_as_string)
                    .ok_or_else(|| {
                        anyhow!(
                            "Could not find hot water source digest for row '{}'",
                            row[0]
                        )
                    })?;
                if hot_water_source.source_is_storage_tank {
                    StringOrNumber::Float(hot_water_source.source_volume.unwrap())
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

struct SummaryDataArgs<'a> {
    timestep_array: &'a [f64],
    input: SummaryInputDigest,
    results_end_user: &'a ResultsEndUser,
    energy_generated_consumed: &'a IndexMap<String, Vec<f64>>,
    energy_to_storage: &'a IndexMap<String, Vec<f64>>,
    energy_from_storage: &'a IndexMap<String, Vec<f64>>,
    storage_from_grid: &'a IndexMap<String, Vec<f64>>,
    energy_diverted: &'a IndexMap<String, Vec<f64>>,
    energy_import: &'a IndexMap<String, Vec<f64>>,
    energy_export: &'a IndexMap<String, Vec<f64>>,
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
            EnergySupplyStatKey::ElecGenerated => self.elec_generated.to_string().into(),
            EnergySupplyStatKey::ElecConsumed => self.elec_consumed.to_string().into(),
            EnergySupplyStatKey::GenToConsumption => self.gen_to_consumption.to_string().into(),
            EnergySupplyStatKey::GridToConsumption => self.grid_to_consumption.to_string().into(),
            EnergySupplyStatKey::GenerationToGrid => self.generation_to_grid.to_string().into(),
            EnergySupplyStatKey::NetImport => self.net_import.to_string().into(),
            EnergySupplyStatKey::GenToStorage => self.gen_to_storage.to_string().into(),
            EnergySupplyStatKey::StorageToConsumption => {
                self.storage_to_consumption.to_string().into()
            }
            EnergySupplyStatKey::StorageFromGrid => self.storage_from_grid.to_string().into(),
            EnergySupplyStatKey::GenToDiverter => self.gen_to_diverter.to_string().into(),
            EnergySupplyStatKey::StorageEff => self.storage_eff.to_string().into(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum EnergySupplyStatKey {
    ElecGenerated,
    ElecConsumed,
    GenToConsumption,
    GridToConsumption,
    GenerationToGrid,
    NetImport,
    GenToStorage,
    StorageToConsumption,
    StorageFromGrid,
    GenToDiverter,
    StorageEff,
}

fn build_summary_data(args: SummaryDataArgs) -> SummaryData {
    let SummaryDataArgs {
        timestep_array,
        input,
        results_end_user,
        energy_generated_consumed,
        energy_to_storage,
        energy_from_storage,
        storage_from_grid,
        energy_diverted,
        energy_import,
        energy_export,
    } = args;

    // Energy Supply breakdown for all EnergySupply objects
    let stats: IndexMap<String, EnergySupplyStat> = results_end_user
        .iter()
        .map(|(key, value)| {
            (key.clone(), {
                let (elec_generated, elec_consumed) = value.iter().fold(
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

                let grid_to_consumption = energy_import[key].iter().sum::<f64>();
                let generation_to_grid = energy_export[key].iter().sum::<f64>().abs();
                let gen_to_storage = energy_to_storage[key].iter().sum::<f64>();
                let storage_to_consumption = energy_from_storage[key].iter().sum::<f64>().abs();
                let gen_to_diverter = energy_diverted[key].iter().sum::<f64>();
                let storage_eff = if gen_to_storage > 0.0 {
                    NumberOrDivisionByZero::Number(
                        storage_to_consumption
                            / (gen_to_storage + storage_from_grid[key].iter().sum::<f64>()),
                    )
                } else {
                    NumberOrDivisionByZero::DivisionByZero
                };

                EnergySupplyStat {
                    elec_generated,
                    elec_consumed,
                    gen_to_consumption: energy_generated_consumed[key].iter().sum::<f64>(),
                    grid_to_consumption,
                    generation_to_grid,
                    net_import: grid_to_consumption - generation_to_grid,
                    gen_to_storage,
                    storage_to_consumption,
                    storage_from_grid: storage_from_grid[key].iter().sum::<f64>().abs(),
                    gen_to_diverter,
                    storage_eff,
                }
            })
        })
        .collect();

    // get peak electricity consumption, and when it happens
    let start_timestep = input.simulation_time.start_time();
    let stepping = input.simulation_time.step;

    // Get energy supply objects with fuel type 'electricity'
    let electricity_keys = &input.electricity_keys;

    // Calculate net import by adding gross import and export figures. Add
    // because export figures are already negative
    let net_import_per_timestep = (0..timestep_array.len())
        .map(|i| {
            electricity_keys
                .iter()
                .map(|k| energy_import[k][i] + energy_export[k][i])
                .sum::<f64>()
        })
        .collect::<Vec<_>>();

    // Find peak electricity consumption
    let peak_elec_consumption = *net_import_per_timestep
        .iter()
        .max_by(|a, b| a.total_cmp(b))
        .expect("Net import per timestep collection was empty.");
    let index_peak_elec_consumption = net_import_per_timestep
        .iter()
        .position(|&x| x == peak_elec_consumption)
        .expect("Could not find index for peak electricity consumption.");

    // must reflect hour or half hour in the year (hour 0 to hour 8759)
    // to work with the dictionary below timestep_to_date
    // hence + start_timestep
    let step_peak_elec_consumption = index_peak_elec_consumption as f64 + start_timestep;

    let mut timestep_to_date: HashMap<usize, HourForTimestep> = Default::default();

    // Set the base for any non-leap year
    let base_time = Utc.with_ymd_and_hms(2023, 1, 1, 0, 0, 0).unwrap();

    let mut step = start_timestep as usize;
    for t in timestep_array.iter() {
        let current_time = base_time + TimeDelta::minutes((*t * 60.).round() as i64);
        timestep_to_date.insert(
            step,
            HourForTimestep {
                month: current_time
                    .format("%b")
                    .to_string()
                    .to_ascii_uppercase()
                    .into(),
                day: current_time.day() as usize,
                hour: (step as f64 % (24. / stepping)) * stepping + 1.,
            },
        );
        step += 1;
    }

    // Delivered energy by end-use and by fuel
    // TODO (from Python) Ensure end_uses not consuming fuel directly are filtered out on this report
    let mut delivered_energy_map: IndexMap<String, IndexMap<String, f64>> =
        IndexMap::from([("total".into(), Default::default())]);
    for (fuel, end_uses) in results_end_user {
        if ["_unmet_demand", "hw cylinder"].contains(&fuel.as_str()) {
            continue;
        }
        let mut fuel_results: IndexMap<String, f64> =
            IndexMap::from([("total".into(), Default::default())]);
        for (end_use, delivered_energy) in end_uses {
            let delivered_energy_sum = delivered_energy.iter().sum::<f64>();
            if delivered_energy_sum >= 0. {
                fuel_results.insert(end_use.to_owned(), delivered_energy_sum);
                *fuel_results
                    .get_mut("total")
                    .expect("Total key was not present in fuel results.") += delivered_energy_sum;
                *delivered_energy_map["total"]
                    .entry(end_use.to_owned())
                    .or_default() += delivered_energy_sum;
                *delivered_energy_map["total"]
                    .entry("total".into())
                    .or_default() += delivered_energy_sum;
            }
        }
        delivered_energy_map.insert(fuel.clone(), fuel_results);
    }

    SummaryData {
        delivered_energy_map,
        stats,
        peak_elec_consumption,
        index_peak_elec_consumption,
        step_peak_elec_consumption,
        timestep_to_date,
    }
}

struct SummaryData {
    delivered_energy_map: IndexMap<String, IndexMap<String, f64>>,
    stats: IndexMap<String, EnergySupplyStat>,
    peak_elec_consumption: f64,
    index_peak_elec_consumption: usize,
    step_peak_elec_consumption: f64,
    timestep_to_date: HashMap<usize, HourForTimestep>,
}

struct StaticOutputFileArgs {
    output_key: String,
    heat_transfer_coefficient: f64,
    heat_loss_parameter: f64,
    heat_capacity_parameter: f64,
    heat_loss_form_factor: f64,
    temp_internal_air: f64,
    temp_external_air: f64,
}

fn write_core_output_file_static(
    output: &impl Output,
    args: StaticOutputFileArgs,
) -> Result<(), anyhow::Error> {
    let StaticOutputFileArgs {
        output_key,
        heat_transfer_coefficient,
        heat_loss_parameter,
        heat_capacity_parameter,
        heat_loss_form_factor,
        temp_internal_air,
        temp_external_air,
    } = args;

    debug!("writing out to {output_key}");

    let writer = output.writer_for_location_key(&output_key, "csv")?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    writer.write_record([
        "Heat transfer coefficient".to_owned(),
        "W / K".to_owned(),
        heat_transfer_coefficient.to_string(),
    ])?;
    writer.write_record([
        "Heat loss parameter".to_owned(),
        "W / m2.K".to_owned(),
        heat_loss_parameter.to_string(),
    ])?;
    writer.write_record([
        "Heat capacity parameter".to_owned(),
        "kJ / m2.K".to_owned(),
        heat_capacity_parameter.to_string(),
    ])?;
    writer.write_record([
        "Heat loss form factor".to_owned(),
        "".to_owned(),
        heat_loss_form_factor.to_string(),
    ])?;
    writer.write_record(["Assumptions used for HTC/HLP calculation:"])?;
    writer.write_record([
        "Internal air temperature".to_owned(),
        "Celsius".to_owned(),
        temp_internal_air.to_string(),
    ])?;
    writer.write_record([
        "External air temperature".to_owned(),
        "Celsius".to_owned(),
        temp_external_air.to_string(),
    ])?;

    debug!("flushing out static CSV");
    writer.flush()?;

    Ok(())
}

struct HeatBalanceOutputFileArgs<'a> {
    output_key: String,
    timestep_array: &'a [f64],
    hour_per_step: f64,
    heat_balance_map: &'a IndexMap<String, IndexMap<String, Vec<f64>>>,
}

fn write_core_output_file_heat_balance(
    output: &impl Output,
    args: HeatBalanceOutputFileArgs,
) -> Result<(), anyhow::Error> {
    let HeatBalanceOutputFileArgs {
        output_key,
        timestep_array,
        hour_per_step,
        heat_balance_map,
    } = args;

    let writer = output.writer_for_location_key(&output_key, "csv")?;
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

    for (t_idx, _) in timestep_array.iter().enumerate() {
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

fn write_core_output_file_heat_source_wet(
    output: &impl Output,
    output_key: &str,
    timestep_array: &[f64],
    heat_source_wet_results: &ResultsPerTimestep,
) -> Result<(), anyhow::Error> {
    // Repeat column headings for each service
    let mut col_headings = vec![String::from("Timestep")];
    let mut col_units_row = vec![String::from("count")];
    let mut columns: IndexMap<String, Vec<(String, Option<String>)>> = Default::default();

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
                .collect::<Vec<String>>(),
        );
        col_units_row.extend(
            service_results
                .keys()
                .cloned()
                .collect::<IndexMap<_, _>>()
                .keys()
                .cloned()
                .collect::<Vec<String>>(),
        );
    }

    let writer = output.writer_for_location_key(output_key, "csv")?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    // Write column headings and units
    writer.write_record(&col_headings)?;
    writer.write_record(&col_units_row)?;

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
    output: &impl Output,
    output_key: &str,
    heat_source_wet_results_annual: &ResultsAnnual,
) -> Result<(), anyhow::Error> {
    let writer = output.writer_for_location_key(output_key, "csv")?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    for (service_name, service_results) in heat_source_wet_results_annual.iter() {
        writer.write_record([service_name.to_string()])?;
        for (name, value) in service_results.iter() {
            writer.write_record([
                name.0.as_bytes(),
                name.1.as_ref().map(|x| x.as_bytes()).unwrap_or_default(),
                value
                    .iter()
                    .map(String::from)
                    .collect_vec()
                    .join(", ")
                    .as_bytes(),
            ])?;
        }
        writer.write_record([""])?;
    }

    Ok(())
}

fn write_core_output_file_emitters_detailed(
    output: &impl Output,
    output_prefix: &str,
    emitters_output_dict: &IndexMap<String, Vec<EmittersDetailedResult>>,
) -> Result<(), anyhow::Error> {
    for (emitter_name, emitters_detailed_results) in emitters_output_dict {
        let output_key = format!("{}{}", output_prefix, emitter_name);
        let writer = output.writer_for_location_key(&output_key, "csv")?;
        let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

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
        ])?;
        for emitters_detailed_result in emitters_detailed_results {
            writer.write_record(emitters_detailed_result.as_string_values())?;
        }
    }

    Ok(())
}

fn write_core_output_file_esh_detailed(
    output: &impl Output,
    output_prefix: &str,
    esh_output: &IndexMap<String, Vec<StorageHeaterDetailedResult>>,
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
        let writer = output.writer_for_location_key(&output_key, "csv")?;
        let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);
        writer.write_record(headings)?;
        writer.write_record(units_row)?;
        for esh_results in esh_output {
            writer.write_record(esh_results.as_string_values())?;
        }
    }

    Ok(())
}

fn write_core_output_file_ventilation_detailed(
    output: &impl Output,
    output_key: &str,
    vent_output_list: &[VentilationDetailedResult],
) -> Result<(), anyhow::Error> {
    let writer = output.writer_for_location_key(output_key, "csv")?;
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
        writer.write_record(ventilation_results.as_string_values())?;
    }

    Ok(())
}

fn write_core_output_file_hot_water_source_summary(
    _output: &impl Output,
    _output_file: &str,
    _hot_water_source_results: &[StorageTankDetailedResult],
) {
    // TODO complete when hot water source results defined
}

const HOURS_TO_END_DEC: f64 = 8760.;

struct HourForTimestep {
    month: String,
    day: usize,
    hour: f64,
}

#[derive(Clone, Debug, PartialEq)]
enum StringOrNumber {
    String(String),
    Float(f64),
    Integer(usize),
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

impl From<f64> for StringOrNumber {
    fn from(value: f64) -> Self {
        StringOrNumber::Float(value)
    }
}

impl From<&f64> for StringOrNumber {
    fn from(value: &f64) -> Self {
        StringOrNumber::Float(*value)
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

/// Utility function for iterating multiple hashmaps with same keys.
fn iterate_maps<'a: 'b, 'b, K: Eq + Hash, V, W>(
    m1: &'a HashMap<K, V>,
    m2: &'b HashMap<K, W>,
) -> impl Iterator<Item = (&'a K, &'a V, &'b W)> {
    m1.iter().map(move |(k, v1)| (k, v1, m2.get(k).unwrap()))
}
