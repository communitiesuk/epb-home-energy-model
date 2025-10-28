use crate::future_homes_standard::input::{InputForProcessing, ingest_for_processing};
#[cfg(feature = "fhs")]
pub use crate::future_homes_standard::{FhsComplianceWrapper, FhsSingleCalcWrapper};

use crate::future_homes_standard::{FhsComplianceWrapper, FhsSingleCalcWrapper};
use anyhow::anyhow;
use hem::corpus::{Corpus, HtcHlpCalculation, calc_htc_hlp};
use hem::errors::{HemCoreError, HemError, PostprocessingError};
use hem::external_conditions::ExternalConditions;
use hem::input::{Input, SchemaReference};
use hem::output::Output;
use hem::read_weather_file::ExternalConditions as ExternalConditionsFromFile;
use hem::{
    CalculationKey, CalculationResultsWithContext, HemResponse, ProjectFlags, RunResults,
    external_conditions_from_input,
};
use rayon::iter::IntoParallelRefIterator;
use std::collections::HashMap;
use std::io::Read;
use std::panic::{AssertUnwindSafe, catch_unwind};
use tracing::{error, instrument};

pub mod future_homes_standard;

/// Common trait for a wrapper of the HEM methodology, which in its preprocessing stage is able to
/// produce inputs for the HEM core for a particular purpose, and in its postprocessing stage is
/// able to output data as a side effect as well as returning a JSON-serializable response that can
/// be consumed or sent out from an API.
pub trait HemWrapper {
    fn apply_preprocessing(
        &self,
        input: InputForProcessing,
        flags: &ProjectFlags,
    ) -> anyhow::Result<HashMap<CalculationKey, Input>>;
    fn apply_postprocessing(
        &self,
        output: &impl Output,
        results: &HashMap<CalculationKey, CalculationResultsWithContext>,
        flags: &ProjectFlags,
    ) -> anyhow::Result<Option<HemResponse>>;
}

/// A HEM wrapper that does nothing, so can be used in cases when the input for core HEM
/// should be passed directly without mutation.
pub struct PassthroughHemWrapper;

impl PassthroughHemWrapper {
    pub fn new() -> Self {
        Self {}
    }
}

impl HemWrapper for PassthroughHemWrapper {
    fn apply_preprocessing(
        &self,
        input: InputForProcessing,
        _flags: &ProjectFlags,
    ) -> anyhow::Result<HashMap<CalculationKey, Input>> {
        Ok(HashMap::from([(
            CalculationKey::Primary,
            input.finalize()?,
        )]))
    }

    fn apply_postprocessing(
        &self,
        _output: &impl Output,
        _results: &HashMap<CalculationKey, CalculationResultsWithContext>,
        _flags: &ProjectFlags,
    ) -> anyhow::Result<Option<HemResponse>> {
        Ok(None)
    }
}

/// An enum to wrap the known wrappers that could be chosen for a given invocation.
pub enum ChosenWrapper {
    Passthrough(PassthroughHemWrapper),
    #[cfg(feature = "fhs")]
    FhsSingleCalc(FhsSingleCalcWrapper),
    #[cfg(feature = "fhs")]
    FhsCompliance(FhsComplianceWrapper),
}

impl HemWrapper for ChosenWrapper {
    fn apply_preprocessing(
        &self,
        input: InputForProcessing,
        flags: &ProjectFlags,
    ) -> anyhow::Result<HashMap<CalculationKey, Input>> {
        match self {
            ChosenWrapper::Passthrough(wrapper) => {
                <PassthroughHemWrapper as HemWrapper>::apply_preprocessing(wrapper, input, flags)
            }
            #[cfg(feature = "fhs")]
            ChosenWrapper::FhsSingleCalc(wrapper) => {
                <FhsSingleCalcWrapper as HemWrapper>::apply_preprocessing(wrapper, input, flags)
            }
            #[cfg(feature = "fhs")]
            ChosenWrapper::FhsCompliance(wrapper) => {
                <FhsComplianceWrapper as HemWrapper>::apply_preprocessing(wrapper, input, flags)
            }
        }
    }

    fn apply_postprocessing(
        &self,
        output: &impl Output,
        results: &HashMap<CalculationKey, CalculationResultsWithContext>,
        flags: &ProjectFlags,
    ) -> anyhow::Result<Option<HemResponse>> {
        match self {
            ChosenWrapper::Passthrough(wrapper) => {
                wrapper.apply_postprocessing(output, results, flags)
            }
            #[cfg(feature = "fhs")]
            ChosenWrapper::FhsSingleCalc(wrapper) => {
                wrapper.apply_postprocessing(output, results, flags)
            }
            #[cfg(feature = "fhs")]
            ChosenWrapper::FhsCompliance(wrapper) => {
                wrapper.apply_postprocessing(output, results, flags)
            }
        }
    }
}

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
        ) -> anyhow::Result<InputForProcessing> {
            let mut input_for_processing = ingest_for_processing(input)?;

            input_for_processing
                .merge_external_conditions_data(external_conditions_data.map(|x| x.into()))?;
            Ok(input_for_processing)
        }

        let input_for_processing =
            ingest_input_and_start_preprocessing(input, external_conditions_data.as_ref())?;

        fn choose_wrapper(flags: &ProjectFlags) -> ChosenWrapper {
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
                    ChosenWrapper::Passthrough(PassthroughHemWrapper::new()) // TODO review
                }
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

        let external_conditions = resolve_external_conditions(&input, external_conditions_data);

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
