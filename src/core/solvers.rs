use anyhow::anyhow;
use eqsolver::single_variable::FDNewton;
use interp::{interp_slice, InterpMode};
use ivp::prelude::*;
use roots::{find_root_secant, SimpleConvergency};
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

// A viable equivalent of scipy.optimize.root
pub(crate) fn root<const ARGCOUNT: usize>(
    fun: impl Fn(f64, [f64; ARGCOUNT]) -> f64,
    x0: f64,
    args: [f64; ARGCOUNT],
    tol: Option<f64>,
) -> anyhow::Result<f64> {
    let tol = tol.unwrap_or(1e-9);

    let mut convergency = SimpleConvergency {
        eps: tol,
        max_iter: 100,
    };

    let x1 = x0 + 5.0;

    let mut wrapper = |x: f64| -> f64 { fun(x, args) };

    let found_root = find_root_secant(x0, x1, &mut wrapper, &mut convergency)
        .map_err(|e| anyhow::anyhow!("Root optimisation failed: {:?}", e))?;

    Ok(found_root)
}

pub fn solve_ivp(
    func: Arc<dyn Fn(f64, &[f64]) -> Vec<f64> + Send + Sync>,
    t_span: (f64, f64),
    y0: &[f64],
    events: Option<TerminatingEvents>,
    rtol: Option<f64>,
    atol: Option<f64>,
) -> anyhow::Result<OdeResult> {
    let rtol = rtol.unwrap_or(1e-3);
    let atol = atol.unwrap_or(1e-6);

    let system = IvpConfig {
        func: func.clone(),
        events,
    };

    let (t0, tf) = t_span;

    let ivp_solver = Ivp::first_order(&system, t0, tf, y0)
        .method(Method::DOPRI5)
        .rtol(rtol)
        .atol(atol);

    let solution = ivp_solver
        .solve()
        .map_err(|e| anyhow!("IVP solver failed: {:?}", e))?;

    Ok(OdeResult {
        y: transpose(solution.y),
        t_events: transpose(solution.t_events),
        t: solution.t,
    })
}

fn transpose(matrix: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    if matrix.is_empty() || matrix[0].is_empty() {
        return Vec::new();
    }

    let num_columns = matrix[0].len();

    (0..num_columns)
        .map(|col_idx| matrix.iter().map(|row| row[col_idx]).collect::<Vec<f64>>())
        .collect()
}

struct IvpConfig {
    func: Arc<dyn Fn(f64, &[f64]) -> Vec<f64> + Send + Sync>,
    events: Option<TerminatingEvents>,
}

impl FirstOrderSystem for IvpConfig {
    fn derivative(&self, x: f64, y: &[f64], dydx: &mut [f64]) {
        let derived = (self.func)(x, y);
        dydx.clone_from_slice(&derived);
    }

    fn events(&self, x: f64, y: &[f64], out: &mut [f64]) {
        if let Some(events) = &self.events {
            for (i, event_func) in events.funcs.iter().enumerate() {
                let result = event_func(x, y);
                out[i] = result;
            }
        }
    }

    fn n_events(&self) -> usize {
        self.events
            .as_ref()
            .map(|events| events.funcs.len())
            .unwrap_or(0)
    }

    fn event_config(&self, _event_index: usize) -> EventConfig {
        let events = self.events.as_ref().expect(
            "Expected event_config only to have been passed if there are events to look up.",
        );

        let mut event_config = EventConfig::new();
        event_config.terminal();
        event_config.direction(events.direction.unwrap_or_default().into());

        event_config
    }
}

pub struct OdeResult {
    pub y: Vec<Vec<f64>>,
    pub t_events: Vec<Vec<f64>>,
    pub t: Vec<f64>,
}

pub struct TerminatingEvents {
    funcs: Vec<Arc<dyn Fn(f64, &[f64]) -> f64 + Send + Sync>>,
    direction: Option<TerminateDirection>,
}

impl TerminatingEvents {
    pub fn new(
        funcs: Vec<Arc<dyn Fn(f64, &[f64]) -> f64 + Send + Sync>>,
        direction: Option<TerminateDirection>,
    ) -> Self {
        Self { funcs, direction }
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub enum TerminateDirection {
    #[default]
    Both,
    Positive,
    Negative,
}

impl From<TerminateDirection> for Direction {
    fn from(direction: TerminateDirection) -> Self {
        match direction {
            TerminateDirection::Both => Direction::All,
            TerminateDirection::Positive => Direction::Positive,
            TerminateDirection::Negative => Direction::Negative,
        }
    }
}

pub fn interp1d<'a>(
    x: Vec<f64>,
    y: Vec<f64>,
    fill_value: Interp1dFillValue,
) -> Arc<dyn Fn(&[f64]) -> Vec<f64> + Send + Sync> {
    let interp_mode = match fill_value {
        Interp1dFillValue::Extrapolate => InterpMode::Extrapolate,
        Interp1dFillValue::FillValues(values) => InterpMode::Constant(values.0),
    };

    let func =
        move |x_new: &[f64]| interp_slice(&x, &y, x_new.try_into().unwrap(), &interp_mode).to_vec();

    Arc::new(func)
}

pub enum Interp1dFillValue {
    Extrapolate,
    FillValues((f64, f64)),
}
