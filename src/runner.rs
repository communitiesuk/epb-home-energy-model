use crate::{CalculationResult, RunInput};
#[cfg(feature = "rayon")]
use rayon::prelude::*;
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

// pub async fn run_map_express<T, R>(
//     client: &Client,
//     state_machine_arn: &str,
//     items: Vec<T>,
// ) -> anyhow::Result<Vec<R>>
// where
//     T: Serialize,
//     R: for<'de> Deserialize<'de>,
// {
//     let out = client
//         .start_sync_execution()
//         .state_machine_arn(state_machine_arn)
//         .input(serde_json::to_string(&BatchInput { items })?)
//         .send()
//         .await?;
//
//     if out.status().as_str() != "SUCCEEDED" {
//         return Err(anyhow::anyhow!(
//             "execution ended with {}: {:?}",
//             out.status().as_str(),
//             out.cause()
//         ));
//     }
//
//     let output = out.output().ok_or_else(|| anyhow::anyhow!("missing output"))?;
//     Ok(serde_json::from_str(output)?)
// }
