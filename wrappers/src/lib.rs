use crate::future_homes_standard::input::InputForProcessing;
#[cfg(feature = "fhs")]
pub use crate::future_homes_standard::{FhsComplianceWrapper, FhsSingleCalcWrapper};

use crate::future_homes_standard::{FhsComplianceWrapper, FhsSingleCalcWrapper};
use hem::input::Input;
use hem::output::Output;
use hem::{CalculationKey, CalculationResultsWithContext, HemResponse, ProjectFlags};
use std::collections::HashMap;

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
            ChosenWrapper::Passthrough(PassthroughHemWrapper::new())
        }
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
