#![allow(clippy::too_many_arguments)]

mod compare_floats;
pub mod core;
pub mod corpus;
pub mod errors;
mod external_conditions;
pub mod input;
pub mod output;
pub mod read_weather_file;
mod simulation_time;
mod statistics;
mod wrappers;

#[macro_use]
extern crate is_close;

use crate::core::units::convert_profile_to_daily;
pub use crate::corpus::RunResults;
use crate::corpus::{
    Corpus, HeatingCoolingSystemResultKey, HotWaterResultKey, HotWaterResultMap, KeyString,
    NumberOrDivisionByZero, ResultsEndUser, ZoneResultKey,
};
use crate::errors::{HemCoreError, HemError, PostprocessingError, NotImplementedError};
use crate::external_conditions::ExternalConditions;
use crate::input::{
    ingest_for_processing, ExternalConditionsInput, HotWaterSourceDetails, Input,
    InputForProcessing,
};
use crate::output::Output;
use crate::read_weather_file::ExternalConditions as ExternalConditionsFromFile;
use crate::simulation_time::SimulationTime;
use crate::statistics::percentile;
#[cfg(feature = "fhs")]
use crate::wrappers::future_homes_standard::{FhsComplianceWrapper, FhsSingleCalcWrapper};
use crate::wrappers::{ChosenWrapper, HemResponse, HemWrapper, PassthroughHemWrapper};
use anyhow::anyhow;
use bitflags::bitflags;
use csv::WriterBuilder;
use indexmap::IndexMap;
use rayon::prelude::*;
use std::borrow::Cow;
use std::collections::HashMap;
use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::io::Read;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::{Arc, LazyLock};
use tracing::{debug, instrument};

pub const HEM_VERSION: &str = "0.30";
pub const HEM_VERSION_DATE: &str = "2024-06-25";
#[cfg(feature = "fhs")]
pub const FHS_VERSION: &str = "0.21";
#[cfg(feature = "fhs")]
pub const FHS_VERSION_DATE: &str = "2024-06-25";

#[instrument(skip_all)]
pub fn run_project(
    input: impl Read + Debug,
    output: impl Output,
    external_conditions_data: Option<ExternalConditionsFromFile>,
    flags: &ProjectFlags,
) -> Result<Option<HemResponse>, HemError> {
    catch_unwind(AssertUnwindSafe(|| {
        // 1. ingest/ parse input and enter preprocessing stage
        #[instrument(skip_all)]
        fn ingest_input_and_start_preprocessing(
            input: impl Read + Debug,
            external_conditions_data: Option<&ExternalConditionsFromFile>,
        ) -> anyhow::Result<InputForProcessing> {
            let mut input_for_processing = ingest_for_processing(input)?;

            input_for_processing
                .merge_external_conditions_data(external_conditions_data.map(|x| x.into()));
            Ok(input_for_processing)
        }

        let input_for_processing =
            ingest_input_and_start_preprocessing(input, external_conditions_data.as_ref())?;

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
            ) -> anyhow::Result<HashMap<CalculationKey, Corpus>> {
                // TODO: parallel iterate this
                iterate_maps(input, external_conditions)
                    .map(|(key, input, external_conditions)| {
                        anyhow::Ok((*key, Corpus::from_inputs(input, Some(external_conditions))?))
                    })
                    .collect()
            }

            build_corpus(&input, &external_conditions).map_err(|e| match e.downcast_ref::<NotImplementedError>() {
                Some(e) => HemError::NotImplemented(e.clone()),
                None => HemError::InvalidRequest(e),
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
            run_hem_calculation(&corpora)
                .map_err(|e| HemError::FailureInCalculation(HemCoreError::new(e)))
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
            output: &impl Output,
            results: &HashMap<CalculationKey, CalculationResultsWithContext>,
            flags: &ProjectFlags,
        ) -> anyhow::Result<()> {
            if let Some(results) = results.get(&CalculationKey::Primary) {
                {
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
                                ..
                            },
                        ..
                    } = results;
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
                }

                if flags.contains(ProjectFlags::HEAT_BALANCE) {
                    for (_hb_name, _hb_map) in results.results.heat_balance_dict.iter() {
                        // TODO: write out heat balance files
                    }
                }

                if flags.contains(ProjectFlags::DETAILED_OUTPUT_HEATING_COOLING) {
                    // TODO: write out heat source wet outputs
                }

                write_core_output_file_summary(output, results.try_into()?)?;

                let corpus = results.context.corpus;

                let (heat_transfer_coefficient, heat_loss_parameter, _, _) = corpus.calc_htc_hlp();
                let heat_capacity_parameter = corpus.calc_hcp();
                let heat_loss_form_factor = corpus.calc_hlff();

                write_core_output_file_static(
                    output,
                    StaticOutputFileArgs {
                        output_key: "results_static".to_string(),
                        heat_transfer_coefficient,
                        heat_loss_parameter,
                        heat_capacity_parameter,
                        heat_loss_form_factor,
                    },
                )?;
            }

            Ok(())
        }

        write_core_output_files(&output, &contextualised_results, flags)?;

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
                .map_or("Error not captured", |v| v)
                .to_owned(),
        )
    })?
}

#[derive(Clone, Copy)]
pub(crate) struct CalculationContext<'a> {
    input: &'a Input,
    corpus: &'a Corpus,
}

#[derive(Clone, Copy)]
pub(crate) struct CalculationResultsWithContext<'a> {
    results: &'a RunResults<'a>,
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
                        .get(&KeyString::from("energy_demand_incl_pipework_loss")?)
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
            input.shading_segments.clone(), //imperfect but this should be quite small...
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
    results_totals: &'a IndexMap<KeyString, Vec<f64>>,
    results_end_user: &'a IndexMap<KeyString, IndexMap<String, Vec<f64>>>,
    energy_import: &'a IndexMap<KeyString, Vec<f64>>,
    energy_export: &'a IndexMap<KeyString, Vec<f64>>,
    energy_generated_consumed: &'a IndexMap<KeyString, Vec<f64>>,
    energy_to_storage: &'a IndexMap<KeyString, Vec<f64>>,
    energy_from_storage: &'a IndexMap<KeyString, Vec<f64>>,
    energy_diverted: &'a IndexMap<KeyString, Vec<f64>>,
    betafactor: &'a IndexMap<KeyString, Vec<f64>>,
    zone_dict: &'a IndexMap<ZoneResultKey, IndexMap<KeyString, Vec<f64>>>,
    zone_list: &'a [KeyString],
    hc_system_dict: &'a IndexMap<HeatingCoolingSystemResultKey, IndexMap<String, Vec<f64>>>,
    hot_water_dict: &'a IndexMap<HotWaterResultKey, HotWaterResultMap>,
    ductwork_gains: &'a IndexMap<KeyString, Vec<f64>>,
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
        energy_diverted,
        betafactor,
        zone_dict,
        zone_list,
        hc_system_dict,
        hot_water_dict,
        ductwork_gains,
    } = args;
    debug!("writing out to {output_key}");
    let writer = output.writer_for_location_key(&output_key)?;
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

    debug!("flushing out CSV");
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
            energy_diverted,
            heat_cop_dict,
            cool_cop_dict,
            dhw_cop_dict,
            ..
        } = value.results;
        Ok(SummaryOutputFileArgs {
            output_key: "results_summary".to_string(),
            input: value.context.input.into(),
            timestep_array,
            results_end_user,
            energy_generated_consumed,
            energy_to_storage,
            energy_from_storage,
            energy_diverted,
            energy_import,
            energy_export,
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
        elec_generated,
        elec_consumed,
        gen_to_consumption,
        gen_to_diverter,
        gen_to_storage,
        grid_to_consumption,
        generation_to_grid,
        storage_to_consumption,
        storage_eff,
        net_import,
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
        energy_diverted,
        energy_import,
        energy_export,
    });

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

    debug!("flushing out summary CSV");
    writer.flush()?;

    Ok(())
}

struct SummaryDataArgs<'a> {
    timestep_array: &'a [f64],
    input: SummaryInputDigest,
    results_end_user: &'a ResultsEndUser,
    energy_generated_consumed: &'a IndexMap<KeyString, Vec<f64>>,
    energy_to_storage: &'a IndexMap<KeyString, Vec<f64>>,
    energy_from_storage: &'a IndexMap<KeyString, Vec<f64>>,
    energy_diverted: &'a IndexMap<KeyString, Vec<f64>>,
    energy_import: &'a IndexMap<KeyString, Vec<f64>>,
    energy_export: &'a IndexMap<KeyString, Vec<f64>>,
}

fn build_summary_data(args: SummaryDataArgs) -> SummaryData {
    let SummaryDataArgs {
        timestep_array,
        input,
        results_end_user,
        energy_generated_consumed,
        energy_to_storage,
        energy_from_storage,
        energy_diverted,
        energy_import,
        energy_export,
    } = args;

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
    let peak_elec_consumption = *net_import_per_timestep
        .iter()
        .max_by(|a, b| a.total_cmp(b))
        .unwrap();
    let index_peak_elec_consumption = net_import_per_timestep
        .iter()
        .position(|&x| x == peak_elec_consumption)
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

    SummaryData {
        delivered_energy_map,
        elec_generated,
        elec_consumed,
        gen_to_consumption,
        gen_to_diverter,
        gen_to_storage,
        grid_to_consumption,
        generation_to_grid,
        storage_to_consumption,
        storage_eff,
        net_import,
        peak_elec_consumption,
        index_peak_elec_consumption,
        step_peak_elec_consumption,
        timestep_to_date,
    }
}

struct SummaryData {
    delivered_energy_map: IndexMap<KeyString, IndexMap<KeyString, f64>>,
    elec_generated: f64,
    elec_consumed: f64,
    gen_to_consumption: f64,
    gen_to_diverter: f64,
    gen_to_storage: f64,
    grid_to_consumption: f64,
    generation_to_grid: f64,
    storage_to_consumption: f64,
    storage_eff: NumberOrDivisionByZero,
    net_import: f64,
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
    } = args;

    debug!("writing out to {output_key}");

    let writer = output.writer_for_location_key(&output_key)?;
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

    debug!("flushing out static CSV");
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

impl From<&ExternalConditionsFromFile> for ExternalConditionsInput {
    fn from(value: &ExternalConditionsFromFile) -> Self {
        Self {
            air_temperatures: Some(value.air_temperatures.clone()),
            wind_speeds: Some(value.wind_speeds.clone()),
            wind_directions: Some(value.wind_directions.clone()),
            diffuse_horizontal_radiation: Some(value.wind_directions.clone()),
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
