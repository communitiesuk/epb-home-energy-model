use eqsolver::single_variable::FDNewton;
use parking_lot::Mutex;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::IntoPyDict;
use roots::{find_root_brent, SimpleConvergency};
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

    // let mut convergency = SimpleConvergency {
    //     eps: xtol,
    //     max_iter: 100, // default for bisect in scipy
    // };
    //
    // // For the time being we use the brent root solver
    // find_root_brent::<f64, _>(
    //     a.max(0.), // ensure first bracket is at least zero
    //     b,
    //     func_modified,
    //     &mut convergency,
    // )
    // .map_err(|e| anyhow::anyhow!(e))
}

// An viable equivalent of scipy.optimize.root
pub(crate) fn root<const ARGCOUNT: usize>(
    fun: impl Fn(f64, [f64; ARGCOUNT]) -> f64,
    x0: f64,
    args: [f64; ARGCOUNT],
    tol: Option<f64>,
) -> anyhow::Result<f64> {
    let mut convergency = SimpleConvergency {
        eps: tol.unwrap_or(1e-8),
        max_iter: 40, // picked intended as reasonably conservative default
    };
    let guess_interval = 5.; // initial guess for guess interval

    // TODO can we use BrentRoot in argmin for this?
    find_root_brent::<f64, _>(
        (x0 - guess_interval).max(0.), // ensure first bracket is at least zero
        x0 + guess_interval,
        |f| fun(f, args),
        &mut convergency,
    )
    .map_err(|e| anyhow::anyhow!(e))
}
