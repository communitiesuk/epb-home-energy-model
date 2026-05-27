use anyhow::anyhow;
use roots::{find_root_brent, SimpleConvergency};
use pyo3::prelude::*;

#[pyclass]
struct RustCallback {
    // We use Box<dyn Fn> to store the closure.
    // Note: It needs to be Send + Sync because it will need to be passed to a Python thread.
    inner: Box<dyn Fn(f64) -> f64 + Send + Sync>,
}

#[pymethods]
impl RustCallback {
    // The `__call__` method makes the object act like a function in Python.
    #[pyo3(signature = (arg, /))]
    fn __call__(&self, arg: f64) -> PyResult<f64> {
        // Execute the boxed closure
        let result = (self.inner)(arg);
        Ok(result)
    }
}

pub(crate) fn fsolve(
    func: Box<dyn Fn(f64) -> f64 + Send + Sync>,
    x0: f64,
) -> anyhow::Result<f64> {
    let callback = RustCallback { inner: func };

    Python::attach(|py| -> PyResult<f64> {
        let optimize = py.import("scipy.optimize")?;
        let fsolve = optimize.getattr("fsolve")?;
        let result = fsolve.call1((callback, x0))?.extract::<f64>()?;
        Ok(result)
    }).map_err( |e| anyhow!(e))
}

// TODO this is from scipy
// Find equivalent function in a Rust library or implement
pub(crate) fn bisect(
    func: impl Fn(f64) -> anyhow::Result<f64>,
    a: f64,
    b: f64,
    xtol: f64,
) -> anyhow::Result<f64> {
    let func_modified = |x| {
        // TODO handle errors here
        func(x).unwrap()
    };

    let mut convergency = SimpleConvergency {
        eps: xtol,
        max_iter: 100, // default for bisect in scipy
    };

    // For the time being we use the brent root solver
    find_root_brent::<f64, _>(
        a.max(0.), // ensure first bracket is at least zero
        b,
        func_modified,
        &mut convergency,
    )
    .map_err(|e| anyhow::anyhow!(e))
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
