use crate::CalculationKey;
use crate::HemWrapper;
use crate::future_homes_standard::fhs_compliance_response::{
    CalculatedComplianceResult, FhsComplianceResponse,
};
use crate::future_homes_standard::future_homes_standard_notional::apply_fhs_notional_preprocessing;
use crate::future_homes_standard::input::InputForProcessing;

use future_homes_standard::{apply_fhs_postprocessing, apply_fhs_preprocessing};
use future_homes_standard_fee::{apply_fhs_fee_postprocessing, apply_fhs_fee_preprocessing};
use hem::output::Output;
use hem::{
    CalculationContext, CalculationResultsWithContext, HemResponse, ProjectFlags, RunResults,
};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::LazyLock;

mod fhs_appliance;
mod fhs_compliance_response;
mod fhs_hw_events;
#[allow(clippy::module_inception)]
pub mod future_homes_standard;
pub(crate) mod future_homes_standard_fee;
pub(crate) mod future_homes_standard_notional;
pub mod input;

/// A HEM wrapper for all single calculations using the FHS wrapper.
pub struct FhsSingleCalcWrapper;

impl FhsSingleCalcWrapper {
    pub fn new() -> Self {
        Self {}
    }
}

impl HemWrapper for FhsSingleCalcWrapper {
    fn apply_preprocessing(
        &self,
        mut input: InputForProcessing,
        flags: &ProjectFlags,
    ) -> anyhow::Result<HashMap<CalculationKey, InputForProcessing>> {
        do_fhs_preprocessing(&mut input, flags)?;
        Ok(HashMap::from([(CalculationKey::Primary, input)]))
    }

    fn apply_postprocessing(
        &self,
        output: &impl Output,
        results: &HashMap<CalculationKey, CalculationResultsWithContext>,
        flags: &ProjectFlags,
    ) -> anyhow::Result<Option<HemResponse>> {
        let results = results
            .get(&CalculationKey::Primary)
            .expect("A primary calculation was expected in the FHS single calc wrapper");
        do_fhs_postprocessing(output, results, flags)
    }
}

/// A HEM wrapper for full FHS compliance calculations.
pub struct FhsComplianceWrapper;

impl FhsComplianceWrapper {
    pub fn new() -> Self {
        Self {}
    }
}

impl HemWrapper for FhsComplianceWrapper {
    fn apply_preprocessing(
        &self,
        input: InputForProcessing,
        _flags: &ProjectFlags,
    ) -> anyhow::Result<HashMap<CalculationKey, InputForProcessing>> {
        vec![input; FHS_COMPLIANCE_CALCULATIONS.len()]
            .into_par_iter()
            .enumerate()
            .map(|(i, mut input)| {
                let (key, flags) = &FHS_COMPLIANCE_CALCULATIONS[i];
                do_fhs_preprocessing(&mut input, flags)?;
                Ok((*key, input))
            })
            .collect::<anyhow::Result<HashMap<CalculationKey, InputForProcessing>>>()
    }

    fn apply_postprocessing(
        &self,
        output: &impl Output,
        results: &HashMap<CalculationKey, CalculationResultsWithContext>,
        _flags: &ProjectFlags,
    ) -> anyhow::Result<Option<HemResponse>> {
        FHS_COMPLIANCE_CALCULATIONS
            .par_iter()
            .map(|(key, flags)| {
                do_fhs_postprocessing(output, &results[key], flags)?;
                Ok(())
            })
            .collect::<anyhow::Result<()>>()?;

        let compliance_result = CalculatedComplianceResult::try_from(results)?;
        let compliance_response = FhsComplianceResponse::build_from(&compliance_result)?;

        Ok(Some(HemResponse::new(compliance_response)))
    }
}

static FHS_COMPLIANCE_CALCULATIONS: LazyLock<[(CalculationKey, ProjectFlags); 4]> =
    LazyLock::new(|| {
        [
            (CalculationKey::Fhs, ProjectFlags::FHS_ASSUMPTIONS),
            (CalculationKey::FhsFee, ProjectFlags::FHS_FEE_ASSUMPTIONS),
            (
                CalculationKey::FhsNotional,
                ProjectFlags::FHS_NOT_A_ASSUMPTIONS | ProjectFlags::FHS_NOT_B_ASSUMPTIONS,
            ),
            (
                CalculationKey::FhsNotionalFee,
                ProjectFlags::FHS_FEE_NOT_A_ASSUMPTIONS | ProjectFlags::FHS_FEE_NOT_B_ASSUMPTIONS,
            ),
        ]
    });

fn do_fhs_preprocessing(
    input_for_processing: &mut InputForProcessing,
    flags: &ProjectFlags,
) -> anyhow::Result<()> {
    // Apply required preprocessing steps, if any
    // TODO (from Python) Implement notional runs (the below treats them the same as the equivalent non-notional runs)
    if flags.intersects(
        ProjectFlags::FHS_NOT_A_ASSUMPTIONS
            | ProjectFlags::FHS_NOT_B_ASSUMPTIONS
            | ProjectFlags::FHS_FEE_NOT_A_ASSUMPTIONS
            | ProjectFlags::FHS_FEE_NOT_B_ASSUMPTIONS,
    ) {
        apply_fhs_notional_preprocessing(
            input_for_processing,
            flags.contains(ProjectFlags::FHS_NOT_A_ASSUMPTIONS),
            flags.contains(ProjectFlags::FHS_NOT_B_ASSUMPTIONS),
            flags.contains(ProjectFlags::FHS_FEE_NOT_A_ASSUMPTIONS),
            flags.contains(ProjectFlags::FHS_FEE_NOT_B_ASSUMPTIONS),
        )?;
    }
    if flags.intersects(
        ProjectFlags::FHS_ASSUMPTIONS
            | ProjectFlags::FHS_NOT_A_ASSUMPTIONS
            | ProjectFlags::FHS_NOT_B_ASSUMPTIONS,
    ) {
        apply_fhs_preprocessing(input_for_processing, Some(false), None)?;
    } else if flags.intersects(
        ProjectFlags::FHS_FEE_ASSUMPTIONS
            | ProjectFlags::FHS_FEE_NOT_A_ASSUMPTIONS
            | ProjectFlags::FHS_FEE_NOT_B_ASSUMPTIONS,
    ) {
        apply_fhs_fee_preprocessing(input_for_processing)?;
    }

    Ok(())
}

fn do_fhs_postprocessing(
    output: &impl Output,
    results: &CalculationResultsWithContext,
    flags: &ProjectFlags,
) -> anyhow::Result<Option<HemResponse>> {
    let input = &results.context.input;
    let RunResults {
        timestep_array,
        results_end_user,
        energy_import,
        energy_export,
        ..
    } = &results.results;

    if flags.intersects(
        ProjectFlags::FHS_ASSUMPTIONS
            | ProjectFlags::FHS_NOT_A_ASSUMPTIONS
            | ProjectFlags::FHS_NOT_B_ASSUMPTIONS,
    ) {
        let notional = flags
            .intersects(ProjectFlags::FHS_NOT_A_ASSUMPTIONS | ProjectFlags::FHS_NOT_B_ASSUMPTIONS);
        apply_fhs_postprocessing(
            input,
            output,
            energy_import,
            energy_export,
            results_end_user,
            timestep_array,
            notional,
        )?;
    } else if flags.intersects(
        ProjectFlags::FHS_FEE_ASSUMPTIONS
            | ProjectFlags::FHS_FEE_NOT_A_ASSUMPTIONS
            | ProjectFlags::FHS_FEE_NOT_B_ASSUMPTIONS,
    ) {
        let CalculationResultsWithContext {
            results,
            context: CalculationContext { corpus, .. },
        } = results;
        apply_fhs_fee_postprocessing(
            output,
            corpus.total_floor_area,
            results.space_heat_demand_total(),
            results.space_cool_demand_total(),
        )?;
    }

    Ok(None)
}
