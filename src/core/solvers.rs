use anyhow::anyhow;
use eqsolver::single_variable::FDNewton;
use itertools::Itertools;
use numpy::ndarray::ArrayD;
use numpy::PyReadonlyArrayDyn;
use parking_lot::Mutex;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::{IntoPyDict, PyDict};
use std::sync::Arc;

pub fn fsolve(func: impl Fn(f64) -> f64 + Copy, x0: f64) -> anyhow::Result<f64> {
    let solver = FDNewton::new(func);

    solver.solve(x0).map_err(|e| anyhow::anyhow!(e))
}

pub mod bisect {
    use std::fmt;

    pub fn bisect(
        func: impl Fn(f64) -> anyhow::Result<f64>,
        a: f64,
        b: f64,
        xtol: f64,
    ) -> Result<(f64, RootResults), BisectError> {
        let mut a = a;
        let mut b = b;

        let rtol = 8.881784197001252e-16; // 4 * f64::EPSILON
        let maxiter = 100;

        let mut func_calls = 0;
        let mut iterations = 0;

        let mut fa = func(a)?;
        let mut fb = func(b)?;
        func_calls += 2;

        // 2. Initial bracket validation
        if fa * fb > 0.0 {
            return Err(BisectError::SignError(
                "f(a) and f(b) must have different signs".to_string(),
            ));
        }

        // Quick check if endpoints are already perfectly zero
        if fa == 0.0 {
            return Ok((
                a,
                RootResults {
                    root: a,
                    iterations,
                    function_calls: func_calls,
                    converged: true,
                    flag: "converged".into(),
                },
            ));
        }
        if fb == 0.0 {
            return Ok((
                b,
                RootResults {
                    root: b,
                    iterations,
                    function_calls: func_calls,
                    converged: true,
                    flag: "converged".into(),
                },
            ));
        }

        // Orient interval boundaries such that f(a) < 0
        if fa > 0.0 {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut fa, &mut fb);
        }

        let mut mid = a;

        // 3. Main Bisection Loop
        while iterations < maxiter {
            iterations += 1;

            // Midpoint calculation designed to minimize floating-point roundoff
            mid = a + (b - a) * 0.5;
            let fmid = func(mid)?;
            func_calls += 1;

            if fmid == 0.0 {
                break;
            } else if fmid < 0.0 {
                a = mid;
            } else {
                b = mid;
            }

            // SciPy's convergence criterion formula
            let delta = (b - a).abs();
            let threshold = xtol + rtol * mid.abs();
            if delta <= threshold {
                break;
            }
        }

        let converged = iterations <= maxiter;
        let flag = if converged {
            "converged".to_string()
        } else {
            format!("Failed to converge after {} iterations.", maxiter)
        };

        let results = RootResults {
            root: mid,
            iterations,
            function_calls: func_calls,
            converged,
            flag,
        };

        if !converged {
            Err(BisectError::ConvergenceError(results))
        } else {
            Ok((mid, results))
        }
    }

    #[derive(Debug, Clone)]
    pub struct RootResults {
        pub root: f64,
        pub iterations: usize,
        pub function_calls: usize,
        pub converged: bool,
        pub flag: String,
    }

    #[derive(Debug)]
    pub enum BisectError {
        SignError(String),
        ConvergenceError(RootResults),
        FuncError(anyhow::Error),
    }

    impl From<anyhow::Error> for BisectError {
        fn from(err: anyhow::Error) -> Self {
            BisectError::FuncError(err)
        }
    }

    impl fmt::Display for BisectError {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            match self {
                BisectError::SignError(msg) => write!(f, "ValueError: {}", msg),
                BisectError::ConvergenceError(res) => write!(f, "RuntimeError: {}", res.flag),
                BisectError::FuncError(err) => write!(f, "Error in passed-in function: {}", err),
            }
        }
    }

    impl std::error::Error for BisectError {}
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
pub fn root(
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

pub fn solve_ivp<const ARGCOUNT: usize>(
    func: Box<dyn Fn(f64, &[f64]) -> PyResult<Vec<f64>> + Send + Sync>,
    t_span: (f64, f64),
    y0: [f64; ARGCOUNT],
    events: Option<Vec<TerminalFunction>>,
    method: Option<&'static str>,
    rtol: Option<f64>,
    atol: Option<f64>,
) -> anyhow::Result<OdeResult> {
    Python::attach(move |py| -> PyResult<OdeResult> {
        let callback = OdeSolverCallback { inner: func };
        let integrate = py.import("scipy.integrate")?;
        let solve_ivp = integrate.getattr("solve_ivp")?;
        let kwargs = PyDict::new(py);
        kwargs.set_item("events", events)?;
        kwargs.set_item("method", method.unwrap_or("RK45"))?;
        kwargs.set_item("rtol", rtol.unwrap_or(1e-3))?;
        kwargs.set_item("atol", atol.unwrap_or(1e-6))?;
        solve_ivp
            .call((callback, t_span, y0), Some(&kwargs))?
            .extract::<OdeResult>()
    })
    .map_err(|e| anyhow!(e))
}

#[pyclass]
pub struct TerminalFunction {
    pub inner: Box<dyn Fn(f64, &[f64]) -> f64 + Send + Sync>,
    is_terminal: Option<bool>,
    direction: Option<f64>,
}

impl TerminalFunction {
    pub fn new(
        inner: Box<dyn Fn(f64, &[f64]) -> f64 + Send + Sync>,
        is_terminal: Option<bool>,
        direction: Option<f64>,
    ) -> Self {
        Self {
            inner,
            is_terminal,
            direction,
        }
    }

    fn is_terminal(&self) -> bool {
        self.is_terminal.unwrap_or(false)
    }

    fn direction(&self) -> f64 {
        self.direction.unwrap_or(0.)
    }
}

#[pymethods]
impl TerminalFunction {
    #[pyo3(signature = (t, y, /))]
    fn __call__(&self, t: f64, y: Vec<f64>) -> PyResult<f64> {
        // Execute the boxed closure
        let result = (self.inner)(t, &y);
        Ok(result)
    }

    fn __getattribute__(&self, name: String) -> PyResult<Option<f64>> {
        Ok(match name.as_str() {
            "terminal" => Some(if self.is_terminal() { 1.0 } else { 0.0 }), // a float value that Python would consider truthy
            "direction" => Some(self.direction()),
            _ => None,
        })
    }
}

#[pyclass]
struct OdeSolverCallback {
    // We use Box<dyn Fn> to store the closure.
    // Note: It needs to be Send + Sync because it will need to be passed to a Python thread.
    inner: Box<dyn Fn(f64, &[f64]) -> PyResult<Vec<f64>> + Send + Sync>,
}

#[pymethods]
impl OdeSolverCallback {
    // The `__call__` method makes the object act like a function in Python.
    #[pyo3(signature = (arg1, arg2, /))]
    fn __call__(&self, arg1: f64, arg2: Vec<f64>) -> PyResult<Vec<f64>> {
        // Execute the boxed closure
        (self.inner)(arg1, &arg2)
    }
}

#[derive(Debug)]
pub struct OdeResult {
    pub y: Vec<ArrayD<f64>>,
    pub t_events: Vec<ArrayD<f64>>,
    pub t: ArrayD<f64>,
}

impl<'a, 'py> FromPyObject<'a, 'py> for OdeResult {
    type Error = PyErr;

    fn extract(ob: Borrowed<'a, 'py, PyAny>) -> Result<Self, Self::Error> {
        let y: Vec<PyReadonlyArrayDyn<'py, f64>> = ob.getattr("y")?.extract()?;
        let t_events: Option<Vec<PyReadonlyArrayDyn<'py, f64>>> =
            ob.getattr("t_events")?.extract()?;

        let y = y.into_iter().map(|y| y.as_array().to_owned()).collect_vec();
        let t_events = t_events.map(|t_events| {
            t_events
                .into_iter()
                .map(|row| row.as_array().to_owned())
                .collect()
        });

        let t: PyReadonlyArrayDyn<'py, f64> = ob.getattr("t")?.extract()?;
        let t = t.as_array().to_owned();

        Ok(Self {
            y,
            t_events: t_events.unwrap_or_default(),
            t,
        })
    }
}

pub fn interp1d(x: &[f64], y: &[f64], fill_value: Interp1dFillValue) -> Py<PyAny> {
    Python::attach(|py| {
        let interpolate = py.import("scipy.interpolate")?;
        let interp1d = interpolate.getattr("interp1d")?;

        let kwargs = PyDict::new(py);
        kwargs.set_item("kind", "linear")?;
        match fill_value {
            Interp1dFillValue::Extrapolate => {
                kwargs.set_item("fill_value", "extrapolate")?;
            }
            Interp1dFillValue::FillValues(values) => {
                kwargs.set_item("fill_value", values)?;
            }
        }
        kwargs.set_item("bounds_error", false)?;
        kwargs.set_item("assume_sorted", true)?;

        interp1d
            .call((x, y), Some(&kwargs))
            .map(|func| func.unbind())
    })
    .unwrap()
}

pub enum Interp1dFillValue {
    Extrapolate,
    FillValues((f64, f64)),
}
