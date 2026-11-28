use crate::{CalculationResult, RunInput};
use std::collections::HashMap;
use std::future::Future;
use std::hash::Hash;

/// Public trait defining how multiple HEM calculations can be run together, whether serially, in parallel on the same machine, or fanned out to different computing nodes.
pub trait RunsBatchedHem {
    type Key: Hash + Eq;

    fn run_batched_hem(
        &self,
        inputs: &HashMap<Self::Key, RunInput>,
    ) -> impl Future<Output = Result<HashMap<Self::Key, CalculationResult>, BatchError>>;
}

pub trait RunsBatchedHemAsync: Send {
    type Key: Hash + Eq + Send + Sync;

    fn run_batched_hem(
        &self,
        inputs: &HashMap<Self::Key, RunInput>,
    ) -> impl Future<Output = Result<HashMap<Self::Key, CalculationResult>, BatchError>>;
}

pub struct BatchError;

#[cfg(feature = "use-aws-step-functions")]
pub mod aws_step_functions {
    use crate::output::Output;
    use crate::runner::{BatchError, RunsBatchedHemAsync};
    use crate::{CalculationResult, RunInput};
    use aws_sdk_sfn::Client;
    use serde::de::Error;
    use serde_json::value::RawValue;
    use std::collections::HashMap;
    use std::hash::Hash;
    use std::marker::PhantomData;
    use std::sync::Arc;
    use tracing::info;

    pub struct AwsStepFunctionHemRunner<'a, T> {
        state_machine_arn: String,
        client: &'a Client,
        _key: PhantomData<T>,
    }

    impl<'a, T> AwsStepFunctionHemRunner<'a, T> {
        pub fn new(state_machine_arn: String, client: &'a Client) -> Self {
            Self {
                state_machine_arn,
                client,
                _key: PhantomData,
            }
        }
    }

    impl<'a, T: Hash + Eq + Send + Sync + Clone> RunsBatchedHemAsync
        for AwsStepFunctionHemRunner<'a, T>
    {
        type Key = T;

        async fn run_batched_hem(
            &self,
            inputs: &HashMap<Self::Key, RunInput<'_>>,
        ) -> Result<HashMap<Self::Key, CalculationResult>, BatchError> {
            info!("in run_batched_hem in AwsStepFunctionHemRunner!");
            let (keys, payloads): (Vec<Self::Key>, Vec<_>) = inputs
                .iter()
                .map(|(key, payload)| (key.clone(), payload))
                .unzip();
            let payloads: Vec<Box<RawValue>> = payloads
                .into_iter()
                .map(|payload| {
                    let json = payload
                        .as_json_string()
                        .map_err(|_| serde_json::Error::custom("Failed to parse JSON"))?;
                    serde_json::from_str::<Box<RawValue>>(&json)
                })
                .collect::<Result<Vec<_>, _>>()
                .map_err(|_| BatchError)?;

            let out = self
                .client
                .start_sync_execution()
                .state_machine_arn(self.state_machine_arn.clone())
                .input(serde_json::to_string(&payloads).map_err(|_| BatchError)?)
                .send()
                .await
                .map_err(|_| BatchError)?;

            if let Some(output) = out.output() {
                let outputs: Vec<Output> = serde_json::from_str(output).map_err(|_| BatchError)?;
                let calculation_results: HashMap<Self::Key, CalculationResult> = outputs
                    .into_iter()
                    .enumerate()
                    .map(|(i, output)| {
                        Ok::<(Self::Key, CalculationResult), BatchError>((
                            keys[i].clone(),
                            CalculationResult::new(
                                Arc::new(
                                    serde_json::from_str(payloads[i].get())
                                        .map_err(|_| BatchError)?,
                                ),
                                output,
                            ),
                        ))
                    })
                    .collect::<Result<_, _>>()?;
                Ok(calculation_results)
            } else {
                Err(BatchError)
            }
        }
    }
}

#[cfg(feature = "rayon")]
pub mod rayon {
    use crate::output_writer::OutputWriter;
    use crate::runner::{BatchError, RunsBatchedHemAsync};
    use crate::{run_project, CalculationResult, RunInput};
    use rayon::prelude::*;
    use std::collections::HashMap;
    use std::hash::Hash;
    use std::marker::PhantomData;
    use std::sync::Arc;

    struct RayonRunner<T, U: OutputWriter> {
        output_writer: U,
        heat_balance: bool,
        detailed_output_heating_cooling: bool,
        _key: PhantomData<T>,
    }

    impl<T, U> RayonRunner<T, U>
    where
        U: OutputWriter,
    {
        pub fn new(
            output_writer: U,
            heat_balance: Option<bool>,
            detailed_output_heating_cooling: Option<bool>,
        ) -> Self {
            let heat_balance = heat_balance.unwrap_or(false);
            let detailed_output_heating_cooling = detailed_output_heating_cooling.unwrap_or(false);

            Self {
                output_writer,
                heat_balance,
                detailed_output_heating_cooling,
                _key: PhantomData,
            }
        }
    }

    impl<T: Hash + Eq + Send + Sync + Clone, U: OutputWriter> RunsBatchedHemAsync
        for RayonRunner<T, U>
    {
        type Key = T;

        async fn run_batched_hem(
            &self,
            inputs: &HashMap<Self::Key, RunInput<'_>>,
        ) -> Result<HashMap<Self::Key, CalculationResult>, BatchError> {
            inputs
                .par_iter()
                .map(|(key, input)| {
                    let finalized = match input {
                        RunInput::Json(json) => json.clone().finalize(),
                        RunInput::Read(read) => {
                            let json_value = serde_json::from_reader(read.read())?;
                            json_value.finalize()
                        }
                    }
                    .map_err(|err| BatchError)?;

                    run_project(
                        finalized,
                        None,
                        None,
                        self.heat_balance,
                        self.detailed_output_heating_cooling,
                    )
                    .map(|result_value| (*key, result_value))
                })
                .collect::<Result<_, _>>()
        }
    }
}

//         let contextualised_results: Result<HashMap<CalculationKey, CalculationResult>, HemError> =
//             inputs_by_key.par_iter()
//                 .map(|(key, input_value)| {
//                     let finalized = input_value.0.clone().finalize()?; // TODO avoid cloning here!
//                     home_energy_model::run_project(finalized, external_conditions_data.clone(), tariff_data_file, heat_balance, detailed_output_heating_cooling)
//                         .map(|result_value| (*key, result_value))
//                     }).collect();
