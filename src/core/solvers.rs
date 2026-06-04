use eqsolver::single_variable::FDNewton;
use numpy::ndarray::ArrayD;
use numpy::PyReadonlyArrayDyn;
use parking_lot::Mutex;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::IntoPyDict;
use std::sync::Arc;

pub(crate) fn fsolve(func: impl Fn(f64) -> f64 + Copy, x0: f64) -> anyhow::Result<f64> {
    let solver = FDNewton::new(func);

    solver.solve(x0).map_err(|e| anyhow::anyhow!(e))
}

#[pyclass]
struct FallibleRustCallback {
    // We use Box<dyn Fn> to store the closure.
    // Note: It needs to be Send + Sync because it will need to be passed to a Python thread.
    inner: Box<dyn Fn(f64) -> anyhow::Result<f64> + Send + Sync>,
    last_error: Arc<Mutex<Option<anyhow::Error>>>,
}

#[pymethods]
impl FallibleRustCallback {
    // The `__call__` method makes the object act like a function in Python.
    #[pyo3(signature = (arg, /))]
    fn __call__(&self, arg: f64) -> PyResult<f64> {
        // Execute the boxed closure
        match (self.inner)(arg) {
            Ok(result) => Ok(result),
            Err(e) => {
                let msg = format!("{e:#}");
                *self.last_error.lock() = Some(e);
                Err(PyRuntimeError::new_err(msg))
            }
        }
    }
}

// TODO this is from scipy
// Find equivalent function in a Rust library or implement
pub(crate) fn bisect<'a>(
    func: Box<dyn Fn(f64) -> anyhow::Result<f64> + Send + Sync>,
    a: f64,
    b: f64,
    xtol: f64,
) -> anyhow::Result<f64> {
    let rust_callback = FallibleRustCallback {
        inner: func,
        last_error: Arc::new(Mutex::new(None)),
    };

    Python::attach(|py| {
        let zeroes_module = py.import("scipy.optimize._zeros_py")?;
        let bisect = zeroes_module.getattr("bisect")?;
        let kwargs = [("xtol", xtol)].into_py_dict(py)?;
        let result = bisect.call((rust_callback, a, b), Some(&kwargs))?;
        Ok(result.extract()?)
    })
    .map_err(move |e: PyErr| anyhow::anyhow!(e))
}

#[pyclass]
struct FallibleRootRustCallback {
    // We use Box<dyn Fn> to store the closure.
    // Note: It needs to be Send + Sync because it will need to be passed to a Python thread.
    inner: Box<dyn Fn(f64, [f64; 3]) -> anyhow::Result<f64> + Send + Sync>,
    last_error: Arc<Mutex<Option<anyhow::Error>>>,
}

#[pymethods]
impl FallibleRootRustCallback {
    // The `__call__` method makes the object act like a function in Python.
    #[pyo3(signature = (x0, args, /))]
    fn __call__(&self, x0: f64, args: [f64; 3]) -> PyResult<f64> {
        // Execute the boxed closure
        match (self.inner)(x0, args) {
            Ok(result) => Ok(result),
            Err(e) => {
                let msg = format!("{e:#}");
                *self.last_error.lock() = Some(e);
                Err(PyRuntimeError::new_err(msg))
            }
        }
    }
}

// An viable equivalent of scipy.optimize.root
pub(crate) fn root(
    fun: Box<dyn Fn(f64, [f64; 3]) -> anyhow::Result<f64> + Send + Sync>,
    x0: f64,
    args: [f64; 3],
    tol: Option<f64>,
) -> anyhow::Result<f64> {
    let rust_callback = FallibleRootRustCallback {
        inner: fun,
        last_error: Arc::new(Mutex::new(None)),
    };

    let result = Python::attach(move |py| {
        let optimize = py.import("scipy.optimize")?;
        let root = optimize.getattr("root")?;

        let kwargs = [("tol", tol)].into_py_dict(py)?;

        Ok(root
            .call((rust_callback, x0, args), Some(&kwargs))?
            .extract::<RootResult>()?)
    })
    .map_err(move |e: PyErr| anyhow::anyhow!(e))?;

    let RootResult { x } = result;

    Ok(x[0])
}

#[derive(Debug)]
struct RootResult {
    x: ArrayD<f64>,
}

impl<'a, 'py> FromPyObject<'a, 'py> for RootResult {
    type Error = PyErr;

    fn extract(ob: Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        let x: PyReadonlyArrayDyn<'py, f64> = ob.getattr("x")?.extract()?;
        let x = x.as_array().to_owned();

        Ok(Self { x })
    }
}
