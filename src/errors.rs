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
    #[error("{0}")]
    NotImplemented(NotImplementedError),
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
    pub fn new(error: anyhow::Error) -> Self {
        Self { error }
    }
}

/// An error representing that an area of functionality has not been implemented.
#[derive(Clone, Debug, Error)]
#[error("Not implemented: {0}")]
pub struct NotImplementedError(String);

impl NotImplementedError {
    pub(crate) fn new(message: &str) -> Self {
        NotImplementedError(message.to_string())
    }
}
