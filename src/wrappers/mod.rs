use crate::input::{Input, InputForProcessing};
use crate::output::Output;
#[cfg(feature = "fhs")]
use crate::wrappers::future_homes_standard::{FhsComplianceWrapper, FhsSingleCalcWrapper};
use crate::{CalculationKey, CalculationResultsWithContext, ProjectFlags};
use std::collections::HashMap;

#[cfg(feature = "fhs")]
pub mod future_homes_standard;

/// Common trait for a wrapper of the HEM methodology, which in its preprocessing stage is able to
/// produce inputs for the HEM core for a particular purpose, and in its postprocessing stage is
/// able to output data as a side effect as well as returning a JSON-serializable response that can
/// be consumed or sent out from an API.
pub(crate) trait HemWrapper {
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
    ) -> anyhow::Result<()>;
}

/// A HEM wrapper that does nothing, so can be used in cases when the input for core HEM
/// should be passed directly without mutation.
pub(crate) struct PassthroughHemWrapper;

impl PassthroughHemWrapper {
    pub(crate) fn new() -> Self {
        Self {}
    }
}

impl HemWrapper for PassthroughHemWrapper {
    fn apply_preprocessing(
        &self,
        input: InputForProcessing,
        _flags: &ProjectFlags,
    ) -> anyhow::Result<HashMap<CalculationKey, Input>> {
        Ok(HashMap::from([(CalculationKey::Primary, input.finalize())]))
    }

    fn apply_postprocessing(
        &self,
        _output: &impl Output,
        _results: &HashMap<CalculationKey, CalculationResultsWithContext>,
        _flags: &ProjectFlags,
    ) -> anyhow::Result<()> {
        Ok(())
    }
}

/// An enum to wrap the known wrappers that could be chosen for a given invocation.
pub(crate) enum ChosenWrapper {
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
            ChosenWrapper::Passthrough(wrapper) => wrapper.apply_preprocessing(input, flags),
            #[cfg(feature = "fhs")]
            ChosenWrapper::FhsSingleCalc(wrapper) => wrapper.apply_preprocessing(input, flags),
            #[cfg(feature = "fhs")]
            ChosenWrapper::FhsCompliance(wrapper) => wrapper.apply_preprocessing(input, flags),
        }
    }

    fn apply_postprocessing(
        &self,
        output: &impl Output,
        results: &HashMap<CalculationKey, CalculationResultsWithContext>,
        flags: &ProjectFlags,
    ) -> anyhow::Result<()> {
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
