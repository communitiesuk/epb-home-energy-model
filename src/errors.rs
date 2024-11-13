use thiserror::Error;

#[derive(Debug, Error)]
pub enum HemError {
    #[error("Request was considered invalid due to error: {0}")]
    InvalidRequest(#[from] anyhow::Error),
    #[error("Uncaught error during wrapper preprocessing: {0}")]
    PanicInWrapper(String),
    #[error("Error identified during HEM calculation: {0}")]
    FailureInCalculation(#[from] HemCoreError),
    #[error("Uncaught error during HEM calculation: {0}")]
    PanicInCalculation(String),
    #[error("Error during wrapper postprocessing: {0}")]
    ErrorInPostprocessing(PostprocessingError),
    #[error("General uncaught error: {0}")]
    GeneralPanic(String),
}

#[derive(Debug, Error)]
#[error(transparent)]
pub struct HemCoreError {
    error: anyhow::Error,
}

impl HemCoreError {
    pub(crate) fn new(error: anyhow::Error) -> Self {
        Self { error }
    }
}

#[derive(Debug, Error)]
#[error(transparent)]
pub struct PostprocessingError {
    error: anyhow::Error,
}

impl PostprocessingError {
    pub(crate) fn new(error: anyhow::Error) -> Self {
        Self { error }
    }
}
