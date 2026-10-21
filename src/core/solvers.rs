use argmin::core::{CostFunction, Executor};
use argmin::solver::brent::BrentRoot;
use eqsolver::single_variable::FDNewton;
use interp::{interp_slice, InterpMode};
use std::sync::Arc;

pub fn fsolve(func: impl Fn(f64) -> f64 + Copy, x0: f64) -> anyhow::Result<f64> {
    let f0 = func(x0);
    if f0.abs() < 1e-12 {
        return Ok(x0);
    }

    if let Ok(result) = FDNewton::new(func).solve(x0) {
        return Ok(result);
    }

    fn bracket_root(func: impl Fn(f64) -> f64, x0: f64, f0: f64) -> anyhow::Result<(f64, f64)> {
        let sign0 = f0.signum();

        // Start with a fine initial step to catch roots very close to x0, then
        // expand geometrically. Starting at 1e-4 rather than 1.0 avoids jumping
        // over narrow roots near the initial guess.
        let mut step = 1e-4_f64;
        let max_step = 1e10_f64;

        while step <= max_step {
            let lo = x0 - step;
            let hi = x0 + step;

            let f_lo = func(lo);
            if f_lo == 0.0 || f_lo.signum() != sign0 {
                return Ok((lo, x0));
            }

            let f_hi = func(hi);
            if f_hi == 0.0 || f_hi.signum() != sign0 {
                return Ok((x0, hi));
            }

            // Also check if lo and hi bracket each other directly (catches cases
            // where the root is symmetric around x0 but x0 itself is not near it).
            if f_lo.signum() != f_hi.signum() {
                return Ok((lo, hi));
            }

            step *= 2.0;
        }

        anyhow::bail!(
            "fsolve: could not bracket a root near x0={x0} after Newton's method failed \
         (f(x0)={f0:.3e}, searched radius up to {max_step:.3e})"
        )
    }

    // Fallback as FDNewton solution can, rarely, fail to converge:
    // expand outward from x0 until we find a sign change bracket,
    // then use Brent's method (superlinear convergence, bisection reliability).
    let (a, b) = bracket_root(func, x0, f0)?;

    struct F<G>(G);
    impl<G: Fn(f64) -> f64> CostFunction for F<G> {
        type Param = f64;
        type Output = f64;
        fn cost(&self, x: &f64) -> Result<f64, argmin::core::Error> {
            Ok(self.0(*x))
        }
    }

    let solver = BrentRoot::new(a, b, 1e-8);
    let executor = Executor::new(F(func), solver);
    executor
        .run()
        .map_err(|e| anyhow::anyhow!("fsolve Brent fallback failed: {e}"))?
        .state()
        .best_param
        .ok_or_else(|| anyhow::anyhow!("fsolve Brent fallback: no root found"))
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
            return Err(BisectError::Sign(
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
            Err(BisectError::Convergence(results))
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
        Sign(String),
        Convergence(RootResults),
        Function(anyhow::Error),
    }

    impl From<anyhow::Error> for BisectError {
        fn from(err: anyhow::Error) -> Self {
            BisectError::Function(err)
        }
    }

    impl fmt::Display for BisectError {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            match self {
                BisectError::Sign(msg) => write!(f, "ValueError: {}", msg),
                BisectError::Convergence(res) => write!(f, "RuntimeError: {}", res.flag),
                BisectError::Function(err) => write!(f, "Error in passed-in function: {}", err),
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
) -> anyhow::Result<f64> {
    let f = |x: f64| -> f64 { fun(x, args) };

    fsolve(f, x0)
}

/// Direct port of scipy.optimize.root_scalar (as of scipy 1.14.1), including ports of the
/// underlying solvers: bisect, brentq, brenth and ridder (implemented in C in scipy), toms748
/// and newton (covering the newton, secant and halley methods).
///
/// Differences from the Python interface:
/// - extra arguments (`args`) are not supported, as closures can capture these instead
/// - the only solver-specific option (`k` for toms748) is a field on the options struct
/// - as `full_output=True` and `disp=False` are always passed by root_scalar to the underlying
///   solvers, non-convergence is reported in the returned `RootResults` rather than as an error,
///   and Python's RuntimeWarnings are logged as warnings
pub mod root_scalar {
    use std::cell::{Cell, RefCell};
    use std::fmt;
    use thiserror::Error;
    use tracing::warn;

    const DEFAULT_ITER: usize = 100;
    const DEFAULT_XTOL: f64 = 2e-12;
    const DEFAULT_RTOL: f64 = 4. * f64::EPSILON;

    const NEWTON_DEFAULT_TOL: f64 = 1.48e-8;
    const NEWTON_DEFAULT_MAXITER: usize = 50;
    const NEWTON_DEFAULT_RTOL: f64 = 0.0;

    const TOMS748_DEFAULT_K: usize = 1;

    const CONVERGED: &str = "converged";
    const CONVERR: &str = "convergence error";

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub enum Method {
        Bisect,
        Brentq,
        Brenth,
        Ridder,
        Toms748,
        Newton,
        Secant,
        Halley,
    }

    impl fmt::Display for Method {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            write!(
                f,
                "{}",
                match self {
                    Method::Bisect => "bisect",
                    Method::Brentq => "brentq",
                    Method::Brenth => "brenth",
                    Method::Ridder => "ridder",
                    Method::Toms748 => "toms748",
                    Method::Newton => "newton",
                    Method::Secant => "secant",
                    Method::Halley => "halley",
                }
            )
        }
    }

    /// The value returned by the objective function. In Python the objective function returns
    /// either a number, or (when `fprime=True` and/or `fprime2=True`) a sequence containing the
    /// value followed by the first and (optionally) second derivatives - here these are
    /// represented by `f64`, `(f64, f64)` and `(f64, f64, f64)` respectively.
    pub trait ObjectiveOutput: Copy {
        fn value(&self) -> f64;
        fn first_derivative(&self) -> Option<f64> {
            None
        }
        fn second_derivative(&self) -> Option<f64> {
            None
        }
    }

    impl ObjectiveOutput for f64 {
        fn value(&self) -> f64 {
            *self
        }
    }

    impl ObjectiveOutput for (f64, f64) {
        fn value(&self) -> f64 {
            self.0
        }
        fn first_derivative(&self) -> Option<f64> {
            Some(self.1)
        }
    }

    impl ObjectiveOutput for (f64, f64, f64) {
        fn value(&self) -> f64 {
            self.0
        }
        fn first_derivative(&self) -> Option<f64> {
            Some(self.1)
        }
        fn second_derivative(&self) -> Option<f64> {
            Some(self.2)
        }
    }

    pub type DerivativeFn<'a> = Box<dyn Fn(f64) -> anyhow::Result<f64> + 'a>;

    /// A scalar function as used internally by the solvers
    type ScalarFn<'a> = dyn Fn(f64) -> Result<f64, RootScalarError> + 'a;

    /// Equivalent of the values that `fprime` and `fprime2` can take in Python.
    pub enum Derivative<'a> {
        /// Equivalent of passing a callable returning the derivative.
        Function(DerivativeFn<'a>),
        /// Equivalent of passing `True`: the objective function returns the derivative(s)
        /// along with its value.
        FromObjective,
    }

    impl<'a> Derivative<'a> {
        pub fn function(func: impl Fn(f64) -> anyhow::Result<f64> + 'a) -> Self {
            Derivative::Function(Box::new(func))
        }
    }

    /// Optional arguments to `root_scalar`, mirroring the keyword arguments of the Python
    /// function. Use `..Default::default()` for any arguments that are not needed.
    #[derive(Default)]
    pub struct RootScalarOptions<'a> {
        pub method: Option<Method>,
        pub bracket: Option<(f64, f64)>,
        pub fprime: Option<Derivative<'a>>,
        pub fprime2: Option<Derivative<'a>>,
        pub x0: Option<f64>,
        pub x1: Option<f64>,
        pub xtol: Option<f64>,
        pub rtol: Option<f64>,
        pub maxiter: Option<usize>,
        /// Corresponds to the `k` solver option for toms748 (passed in `options` in Python)
        pub k: Option<usize>,
    }

    #[derive(Clone, Debug, PartialEq)]
    pub struct RootResults {
        pub root: f64,
        /// This is `None` where scipy reports NaN iterations (i.e. when the solver was stopped
        /// due to the function returning NaN)
        pub iterations: Option<usize>,
        pub function_calls: usize,
        pub converged: bool,
        pub flag: String,
        pub method: Method,
    }

    #[derive(Debug, Error)]
    pub enum RootScalarError {
        #[error("ValueError: {0}")]
        Value(String),
        #[error("TypeError: {0}")]
        Type(String),
        #[error("Error in passed-in function: {0}")]
        Function(#[from] anyhow::Error),
    }

    /// Find a root of a scalar function.
    ///
    /// Port of scipy.optimize.root_scalar - see its documentation for details of the available
    /// methods and the arguments each of them requires.
    pub fn root_scalar<T: ObjectiveOutput>(
        f: impl Fn(f64) -> anyhow::Result<T>,
        options: RootScalarOptions,
    ) -> Result<RootResults, RootScalarError> {
        let RootScalarOptions {
            method,
            bracket,
            fprime,
            fprime2,
            x0,
            x1,
            xtol,
            rtol,
            maxiter,
            k,
        } = options;

        // fun also returns the derivative(s)
        let (fprime, fprime2) = match (fprime, fprime2) {
            (_, Some(Derivative::FromObjective)) => (
                Some(Derivative::FromObjective),
                Some(Derivative::FromObjective),
            ),
            (fprime, fprime2) => (fprime, fprime2),
        };
        let is_memoized = matches!(fprime, Some(Derivative::FromObjective));
        let memoized = MemoizeDer::new(&f);

        let func = |x: f64| -> Result<f64, RootScalarError> {
            if is_memoized {
                memoized.call(x)
            } else {
                Ok(f(x)?.value())
            }
        };
        let fprime: Option<Box<ScalarFn>> = match fprime {
            Some(Derivative::Function(fprime)) => Some(Box::new(move |x| Ok(fprime(x)?))),
            Some(Derivative::FromObjective) => Some(Box::new(|x| memoized.fprime(x))),
            None => None,
        };
        let fprime2: Option<Box<ScalarFn>> = match fprime2 {
            Some(Derivative::Function(fprime2)) => Some(Box::new(move |x| Ok(fprime2(x)?))),
            Some(Derivative::FromObjective) => Some(Box::new(|x| memoized.fprime2(x))),
            None => None,
        };

        // Pick a method if not specified.
        // Use the "best" method available for the situation.
        let method = match method {
            Some(method) => method,
            None => {
                if bracket.is_some() {
                    Method::Brentq
                } else if x0.is_some() {
                    if fprime.is_some() {
                        if fprime2.is_some() {
                            Method::Halley
                        } else {
                            Method::Newton
                        }
                    } else if x1.is_some() {
                        Method::Secant
                    } else {
                        Method::Newton
                    }
                } else {
                    return Err(RootScalarError::Value(
                        "Unable to select a solver as neither bracket nor starting point provided."
                            .into(),
                    ));
                }
            }
        };

        // toms748 is the only underlying solver that accepts `k`
        if k.is_some() && method != Method::Toms748 {
            let underlying = match method {
                Method::Secant | Method::Halley => Method::Newton,
                method => method,
            };
            return Err(RootScalarError::Type(format!(
                "{underlying}() got an unexpected keyword argument 'k'"
            )));
        }

        let mut sol = match method {
            Method::Bisect | Method::Ridder | Method::Brentq | Method::Brenth | Method::Toms748 => {
                let (a, b) = bracket.ok_or_else(|| {
                    RootScalarError::Value(format!("Bracket needed for {method}"))
                })?;
                let xtol = xtol.unwrap_or(DEFAULT_XTOL);
                let rtol = rtol.unwrap_or(DEFAULT_RTOL);
                let maxiter = maxiter.unwrap_or(DEFAULT_ITER);

                let result = match method {
                    Method::Toms748 => toms748(
                        &func,
                        a,
                        b,
                        k.unwrap_or(TOMS748_DEFAULT_K),
                        xtol,
                        rtol,
                        maxiter,
                    ),
                    _ => {
                        let solver = match method {
                            Method::Bisect => zeros::bisect,
                            Method::Ridder => zeros::ridder,
                            Method::Brentq => zeros::brentq,
                            Method::Brenth => zeros::brenth,
                            _ => unreachable!(),
                        };
                        call_solver(solver, &func, a, b, xtol, rtol, maxiter)
                    }
                };

                match result {
                    Ok(raw) => raw.into_results(method),
                    // scipy wraps the function so that it raises an error (rather than returning
                    // incorrect results) when the callable returns a NaN, so the iteration count
                    // is not available.
                    Err(Interrupt::NanValue { x, function_calls }) => RootResults {
                        root: x,
                        iterations: None,
                        function_calls,
                        converged: false,
                        flag: nan_value_message(x),
                        method,
                    },
                    Err(Interrupt::Error(e)) => return Err(e),
                }
            }
            Method::Secant => {
                let x0 = x0.ok_or_else(|| {
                    RootScalarError::Value(format!("x0 must not be None for {method}"))
                })?;
                newton(
                    &func,
                    x0,
                    None,
                    xtol.unwrap_or(NEWTON_DEFAULT_TOL),
                    maxiter.unwrap_or(NEWTON_DEFAULT_MAXITER),
                    None,
                    x1,
                    rtol.unwrap_or(NEWTON_DEFAULT_RTOL),
                )?
            }
            Method::Newton => {
                let x0 = x0.ok_or_else(|| {
                    RootScalarError::Value(format!("x0 must not be None for {method}"))
                })?;
                // approximate fprime with finite differences if not provided
                let approx_fprime = |x: f64| -> Result<f64, RootScalarError> {
                    approx_derivative_2_point(|x| Ok(f(x)?.value()), x)
                };
                let fprime: &ScalarFn = match fprime.as_deref() {
                    Some(fprime) => fprime,
                    None => &approx_fprime,
                };
                newton(
                    &func,
                    x0,
                    Some(fprime),
                    xtol.unwrap_or(NEWTON_DEFAULT_TOL),
                    maxiter.unwrap_or(NEWTON_DEFAULT_MAXITER),
                    None,
                    None,
                    rtol.unwrap_or(NEWTON_DEFAULT_RTOL),
                )?
            }
            Method::Halley => {
                let x0 = x0.ok_or_else(|| {
                    RootScalarError::Value(format!("x0 must not be None for {method}"))
                })?;
                let fprime = fprime.as_deref().ok_or_else(|| {
                    RootScalarError::Value(format!("fprime must be specified for {method}"))
                })?;
                let fprime2 = fprime2.as_deref().ok_or_else(|| {
                    RootScalarError::Value(format!("fprime2 must be specified for {method}"))
                })?;
                newton(
                    &func,
                    x0,
                    Some(fprime),
                    xtol.unwrap_or(NEWTON_DEFAULT_TOL),
                    maxiter.unwrap_or(NEWTON_DEFAULT_MAXITER),
                    Some(fprime2),
                    None,
                    rtol.unwrap_or(NEWTON_DEFAULT_RTOL),
                )?
            }
        };

        if is_memoized {
            // Replace the function_calls count with the memoized count.
            // Avoids double and triple-counting.
            sol.function_calls = memoized.n_calls.get();
        }

        Ok(sol)
    }

    /// Caches the value and derivative(s) of the function each time it is called.
    ///
    /// This is a simplistic memoizer that calls and caches a single value of `f(x)`.
    struct MemoizeDer<'a, T, F> {
        fun: &'a F,
        vals: RefCell<Option<(f64, T)>>,
        n_calls: Cell<usize>,
    }

    impl<'a, T: ObjectiveOutput, F: Fn(f64) -> anyhow::Result<T>> MemoizeDer<'a, T, F> {
        fn new(fun: &'a F) -> Self {
            Self {
                fun,
                vals: RefCell::new(None),
                n_calls: Cell::new(0),
            }
        }

        fn vals(&self, x: f64) -> Result<T, RootScalarError> {
            // Derivative may be requested before the function itself, always check
            if let Some((cached_x, vals)) = *self.vals.borrow() {
                if x == cached_x {
                    return Ok(vals);
                }
            }
            let fg = (self.fun)(x)?;
            *self.vals.borrow_mut() = Some((x, fg));
            self.n_calls.set(self.n_calls.get() + 1);
            Ok(fg)
        }

        fn call(&self, x: f64) -> Result<f64, RootScalarError> {
            Ok(self.vals(x)?.value())
        }

        fn fprime(&self, x: f64) -> Result<f64, RootScalarError> {
            self.vals(x)?.first_derivative().ok_or_else(|| {
                RootScalarError::Type(
                    "fprime=True requires the function to return its first derivative".into(),
                )
            })
        }

        fn fprime2(&self, x: f64) -> Result<f64, RootScalarError> {
            self.vals(x)?.second_derivative().ok_or_else(|| {
                RootScalarError::Type(
                    "fprime2=True requires the function to return its second derivative".into(),
                )
            })
        }
    }

    /// Equivalent of scipy.optimize._numdiff.approx_derivative(f, x0, method='2-point') for a
    /// scalar function of a scalar, with no bounds.
    fn approx_derivative_2_point(
        fun: impl Fn(f64) -> Result<f64, RootScalarError>,
        x0: f64,
    ) -> Result<f64, RootScalarError> {
        let f0 = fun(x0)?;
        // sign_x0 needs to be 1 when x0 == 0
        let sign_x0 = if x0 >= 0. { 1. } else { -1. };
        let rstep = f64::EPSILON.sqrt();
        let h = rstep * sign_x0 * np_maximum(1.0, x0.abs());
        let x1 = x0 + h;
        let dx = x1 - x0; // Recompute dx as exactly representable number.
        let df = fun(x1)? - f0;
        Ok(df / dx)
    }

    /// Raw result tuple from a solver: (root, function calls, iterations, converged)
    struct RawResult {
        root: f64,
        function_calls: usize,
        iterations: usize,
        converged: bool,
    }

    impl RawResult {
        fn into_results(self, method: Method) -> RootResults {
            RootResults {
                root: self.root,
                iterations: Some(self.iterations),
                function_calls: self.function_calls,
                converged: self.converged,
                flag: if self.converged { CONVERGED } else { CONVERR }.into(),
                method,
            }
        }
    }

    /// Ways in which a bracketing solver can be interrupted
    enum Interrupt {
        /// The function returned NaN (equivalent of the error raised by scipy's _wrap_nan_raise)
        NanValue {
            x: f64,
            function_calls: usize,
        },
        Error(RootScalarError),
    }

    impl From<RootScalarError> for Interrupt {
        fn from(err: RootScalarError) -> Self {
            Interrupt::Error(err)
        }
    }

    fn value_error(msg: String) -> Interrupt {
        Interrupt::Error(RootScalarError::Value(msg))
    }

    fn nan_value_message(x: f64) -> String {
        format!(
            "The function value at x={} is NaN; solver cannot continue.",
            py_repr(x)
        )
    }

    /// Equivalent of scipy's _wrap_nan_raise
    struct NanRaise<'a> {
        f: &'a ScalarFn<'a>,
        function_calls: Cell<usize>,
    }

    impl<'a> NanRaise<'a> {
        fn new(f: &'a ScalarFn) -> Self {
            Self {
                f,
                function_calls: Cell::new(0),
            }
        }

        fn call(&self, x: f64) -> Result<f64, Interrupt> {
            let fx = (self.f)(x)?;
            self.function_calls.set(self.function_calls.get() + 1);
            if fx.is_nan() {
                return Err(Interrupt::NanValue {
                    x,
                    function_calls: self.function_calls.get(),
                });
            }
            Ok(fx)
        }
    }

    type ZerosSolver = fn(
        &dyn Fn(f64) -> Result<f64, Interrupt>,
        f64,
        f64,
        f64,
        f64,
        usize,
        &mut zeros::SolverStats,
    ) -> Result<f64, Interrupt>;

    /// Equivalent of the Python wrappers bisect/ridder/brentq/brenth plus the call_solver
    /// function in scipy's zeros.c
    fn call_solver(
        solver: ZerosSolver,
        f: &ScalarFn,
        a: f64,
        b: f64,
        xtol: f64,
        rtol: f64,
        maxiter: usize,
    ) -> Result<RawResult, Interrupt> {
        if xtol <= 0. {
            return Err(value_error(format!(
                "xtol too small ({} <= 0)",
                py_format_g(xtol)
            )));
        }
        if rtol < DEFAULT_RTOL {
            return Err(value_error(format!(
                "rtol too small ({} < {})",
                py_format_g(rtol),
                py_format_g(DEFAULT_RTOL)
            )));
        }
        let f = NanRaise::new(f);

        let mut solver_stats = zeros::SolverStats::default();
        let zero = solver(&|x| f.call(x), a, b, xtol, rtol, maxiter, &mut solver_stats)?;

        let converged = match solver_stats.error_num {
            zeros::ErrorNum::Converged => true,
            zeros::ErrorNum::SignErr => {
                return Err(value_error(
                    "f(a) and f(b) must have different signs".into(),
                ))
            }
            zeros::ErrorNum::ConvErr => false,
            zeros::ErrorNum::InProgress => unreachable!("solvers always set a final status"),
        };

        Ok(RawResult {
            root: zero,
            function_calls: solver_stats.funcalls,
            iterations: solver_stats.iterations,
            converged,
        })
    }

    /// Ports of the C solvers in scipy/optimize/Zeros
    mod zeros {
        use super::Interrupt;

        #[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
        pub(super) enum ErrorNum {
            Converged,
            SignErr,
            ConvErr,
            #[default]
            InProgress,
        }

        // NB. in the C code iterations is left uninitialised when the solver returns early
        // (i.e. if a bracket endpoint is a root) - here it is reported as 0
        #[derive(Debug, Default)]
        pub(super) struct SolverStats {
            pub(super) funcalls: usize,
            pub(super) iterations: usize,
            pub(super) error_num: ErrorNum,
        }

        fn signbit(x: f64) -> bool {
            x.is_sign_negative()
        }

        pub(super) fn bisect(
            f: &dyn Fn(f64) -> Result<f64, Interrupt>,
            mut xa: f64,
            xb: f64,
            xtol: f64,
            rtol: f64,
            iter: usize,
            solver_stats: &mut SolverStats,
        ) -> Result<f64, Interrupt> {
            solver_stats.error_num = ErrorNum::InProgress;

            let fa = f(xa)?;
            let fb = f(xb)?;
            solver_stats.funcalls = 2;
            if fa == 0. {
                solver_stats.error_num = ErrorNum::Converged;
                return Ok(xa);
            }
            if fb == 0. {
                solver_stats.error_num = ErrorNum::Converged;
                return Ok(xb);
            }
            if signbit(fa) == signbit(fb) {
                solver_stats.error_num = ErrorNum::SignErr;
                return Ok(0.);
            }
            let mut dm = xb - xa;
            solver_stats.iterations = 0;
            for _ in 0..iter {
                solver_stats.iterations += 1;
                dm *= 0.5;
                let xm = xa + dm;
                let fm = f(xm)?;
                solver_stats.funcalls += 1;
                if signbit(fm) == signbit(fa) {
                    xa = xm;
                }
                if fm == 0. || dm.abs() < xtol + rtol * xm.abs() {
                    solver_stats.error_num = ErrorNum::Converged;
                    return Ok(xm);
                }
            }
            solver_stats.error_num = ErrorNum::ConvErr;
            Ok(xa)
        }

        pub(super) fn ridder(
            f: &dyn Fn(f64) -> Result<f64, Interrupt>,
            mut xa: f64,
            mut xb: f64,
            xtol: f64,
            rtol: f64,
            iter: usize,
            solver_stats: &mut SolverStats,
        ) -> Result<f64, Interrupt> {
            fn sign(a: f64) -> f64 {
                if a > 0. {
                    1.
                } else {
                    -1.
                }
            }

            let mut xn = 0.0;

            solver_stats.error_num = ErrorNum::InProgress;

            let mut tol = xtol + rtol * c_min(xa.abs(), xb.abs());
            let mut fa = f(xa)?;
            let mut fb = f(xb)?;
            solver_stats.funcalls = 2;
            if fa == 0. {
                solver_stats.error_num = ErrorNum::Converged;
                return Ok(xa);
            }
            if fb == 0. {
                solver_stats.error_num = ErrorNum::Converged;
                return Ok(xb);
            }
            if signbit(fa) == signbit(fb) {
                solver_stats.error_num = ErrorNum::SignErr;
                return Ok(0.);
            }

            solver_stats.iterations = 0;
            for _ in 0..iter {
                solver_stats.iterations += 1;
                let dm = 0.5 * (xb - xa);
                let xm = xa + dm;
                let fm = f(xm)?;
                let dn = sign(fb - fa) * dm * fm / (fm * fm - fa * fb).sqrt();
                xn = xm - sign(dn) * c_min(dn.abs(), dm.abs() - 0.5 * tol);
                let fn_ = f(xn)?;
                solver_stats.funcalls += 2;
                if signbit(fn_) != signbit(fm) {
                    xa = xn;
                    fa = fn_;
                    xb = xm;
                    fb = fm;
                } else if signbit(fn_) != signbit(fa) {
                    xb = xn;
                    fb = fn_;
                } else {
                    xa = xn;
                    fa = fn_;
                }
                tol = xtol + rtol * xn;
                if fn_ == 0.0 || (xb - xa).abs() < tol {
                    solver_stats.error_num = ErrorNum::Converged;
                    return Ok(xn);
                }
            }
            solver_stats.error_num = ErrorNum::ConvErr;
            Ok(xn)
        }

        pub(super) fn brentq(
            f: &dyn Fn(f64) -> Result<f64, Interrupt>,
            xa: f64,
            xb: f64,
            xtol: f64,
            rtol: f64,
            iter: usize,
            solver_stats: &mut SolverStats,
        ) -> Result<f64, Interrupt> {
            brent_impl(
                f,
                xa,
                xb,
                xtol,
                rtol,
                iter,
                solver_stats,
                BrentVariant::Quadratic,
            )
        }

        pub(super) fn brenth(
            f: &dyn Fn(f64) -> Result<f64, Interrupt>,
            xa: f64,
            xb: f64,
            xtol: f64,
            rtol: f64,
            iter: usize,
            solver_stats: &mut SolverStats,
        ) -> Result<f64, Interrupt> {
            brent_impl(
                f,
                xa,
                xb,
                xtol,
                rtol,
                iter,
                solver_stats,
                BrentVariant::Hyperbolic,
            )
        }

        #[derive(Clone, Copy)]
        enum BrentVariant {
            /// brentq: inverse quadratic extrapolation
            Quadratic,
            /// brenth: hyperbolic extrapolation
            Hyperbolic,
        }

        // brentq.c and brenth.c are identical other than the extrapolation step
        //
        // At the top of the loop the situation is the following:
        //
        //  1. the root is bracketed between xa and xb
        //  2. xa is the most recent estimate
        //  3. xp is the previous estimate
        //  4. |fp| < |fb|
        //
        // The order of xa and xp doesn't matter, but assume xp < xb. Then xa lies to
        // the right of xp and the assumption is that xa is increasing towards the root.
        // In this situation we will attempt quadratic extrapolation as long as the
        // condition
        //
        //  *  |fa| < |fp| < |fb|
        //
        // is satisfied. That is, the function value is decreasing as we go along.
        // Note the 4 above implies that the right inequlity already holds.
        //
        // The first check is that xa is still to the left of the root. If not, xb is
        // replaced by xp and the interval reverses, with xb < xa. In this situation
        // we will try linear interpolation. That this has happened is signaled by the
        // equality xb == xp;
        //
        // The second check is that |fa| < |fb|. If this is not the case, we swap
        // xa and xb and resort to bisection.
        #[allow(clippy::too_many_arguments)]
        fn brent_impl(
            f: &dyn Fn(f64) -> Result<f64, Interrupt>,
            xa: f64,
            xb: f64,
            xtol: f64,
            rtol: f64,
            iter: usize,
            solver_stats: &mut SolverStats,
            variant: BrentVariant,
        ) -> Result<f64, Interrupt> {
            let mut xpre = xa;
            let mut xcur = xb;
            let mut xblk = 0.;
            let mut fblk = 0.;
            let mut spre = 0.;
            let mut scur = 0.;
            // the tolerance is 2*delta
            let mut delta;

            solver_stats.error_num = ErrorNum::InProgress;

            let mut fpre = f(xpre)?;
            let mut fcur = f(xcur)?;
            solver_stats.funcalls = 2;
            if fpre == 0. {
                solver_stats.error_num = ErrorNum::Converged;
                return Ok(xpre);
            }
            if fcur == 0. {
                solver_stats.error_num = ErrorNum::Converged;
                return Ok(xcur);
            }
            if signbit(fpre) == signbit(fcur) {
                solver_stats.error_num = ErrorNum::SignErr;
                return Ok(0.);
            }
            solver_stats.iterations = 0;
            for _ in 0..iter {
                solver_stats.iterations += 1;
                if fpre != 0. && fcur != 0. && (signbit(fpre) != signbit(fcur)) {
                    xblk = xpre;
                    fblk = fpre;
                    scur = xcur - xpre;
                    spre = scur;
                }
                if fblk.abs() < fcur.abs() {
                    xpre = xcur;
                    xcur = xblk;
                    xblk = xpre;

                    fpre = fcur;
                    fcur = fblk;
                    fblk = fpre;
                }

                delta = (xtol + rtol * xcur.abs()) / 2.;
                let sbis = (xblk - xcur) / 2.;
                if fcur == 0. || sbis.abs() < delta {
                    solver_stats.error_num = ErrorNum::Converged;
                    return Ok(xcur);
                }

                if spre.abs() > delta && fcur.abs() < fpre.abs() {
                    let stry = if xpre == xblk {
                        // interpolate
                        -fcur * (xcur - xpre) / (fcur - fpre)
                    } else {
                        // extrapolate
                        let dpre = (fpre - fcur) / (xpre - xcur);
                        let dblk = (fblk - fcur) / (xblk - xcur);
                        match variant {
                            BrentVariant::Quadratic => {
                                -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre))
                            }
                            BrentVariant::Hyperbolic => {
                                -fcur * (fblk - fpre) / (fblk * dpre - fpre * dblk)
                            }
                        }
                    };
                    if 2. * stry.abs() < c_min(spre.abs(), 3. * sbis.abs() - delta) {
                        // good short step
                        spre = scur;
                        scur = stry;
                    } else {
                        // bisect
                        spre = sbis;
                        scur = sbis;
                    }
                } else {
                    // bisect
                    spre = sbis;
                    scur = sbis;
                }

                xpre = xcur;
                fpre = fcur;
                if scur.abs() > delta {
                    xcur += scur;
                } else {
                    xcur += if sbis > 0. { delta } else { -delta };
                }

                fcur = f(xcur)?;
                solver_stats.funcalls += 1;
            }
            solver_stats.error_num = ErrorNum::ConvErr;
            Ok(xcur)
        }

        /// Equivalent of the C macro MIN(a, b) ((a) < (b) ? (a) : (b))
        fn c_min(a: f64, b: f64) -> f64 {
            if a < b {
                a
            } else {
                b
            }
        }
    }

    /// Port of scipy.optimize.newton for a scalar x0 (with full_output=True, disp=False)
    #[allow(clippy::too_many_arguments)]
    fn newton(
        func: &ScalarFn,
        x0: f64,
        fprime: Option<&ScalarFn>,
        tol: f64,
        maxiter: usize,
        fprime2: Option<&ScalarFn>,
        x1: Option<f64>,
        rtol: f64,
    ) -> Result<RootResults, RootScalarError> {
        if tol <= 0. {
            return Err(RootScalarError::Value(format!(
                "tol too small ({} <= 0)",
                py_format_g(tol)
            )));
        }
        if maxiter < 1 {
            return Err(RootScalarError::Value(
                "maxiter must be greater than 0".into(),
            ));
        }

        let results = |root: f64,
                       function_calls: usize,
                       iterations: usize,
                       converged: bool,
                       method: Method| {
            Ok(RawResult {
                root,
                function_calls,
                iterations,
                converged,
            }
            .into_results(method))
        };

        let mut p0 = x0;
        let mut funcalls = 0;
        // maxiter >= 1, so these will always be set by the loops below
        let mut p = f64::NAN;
        let mut itr = 0;
        let method;

        if let Some(fprime) = fprime {
            // Newton-Raphson method
            let mut method_ = Method::Newton;
            for itr_ in 0..maxiter {
                itr = itr_;
                // first evaluate fval
                let fval = func(p0)?;
                funcalls += 1;
                // If fval is 0, a root has been found, then terminate
                if fval == 0. {
                    return results(p0, funcalls, itr, true, method_);
                }
                let fder = fprime(p0)?;
                funcalls += 1;
                if fder == 0. {
                    warn!("Derivative was zero.");
                    return results(p0, funcalls, itr + 1, false, method_);
                }
                let mut newton_step = fval / fder;
                if let Some(fprime2) = fprime2 {
                    let fder2 = fprime2(p0)?;
                    funcalls += 1;
                    method_ = Method::Halley;
                    // Halley's method:
                    //   newton_step /= (1.0 - 0.5 * newton_step * fder2 / fder)
                    // Only do it if denominator stays close enough to 1
                    // Rationale: If 1-adj < 0, then Halley sends x in the
                    // opposite direction to Newton. Doesn't happen if x is close
                    // enough to root.
                    let adj = newton_step * fder2 / fder / 2.;
                    if adj.abs() < 1. {
                        newton_step /= 1.0 - adj;
                    }
                }
                p = p0 - newton_step;
                if np_isclose(p, p0, rtol, tol) {
                    return results(p, funcalls, itr + 1, true, method_);
                }
                p0 = p;
            }
            method = method_;
        } else {
            // Secant method
            method = Method::Secant;
            let mut p1 = match x1 {
                Some(x1) => {
                    if x1 == x0 {
                        return Err(RootScalarError::Value("x1 and x0 must be different".into()));
                    }
                    x1
                }
                None => {
                    let eps = 1e-4;
                    let p1 = x0 * (1. + eps);
                    p1 + if p1 >= 0. { eps } else { -eps }
                }
            };
            let mut q0 = func(p0)?;
            funcalls += 1;
            let mut q1 = func(p1)?;
            funcalls += 1;
            if q1.abs() < q0.abs() {
                std::mem::swap(&mut p0, &mut p1);
                std::mem::swap(&mut q0, &mut q1);
            }
            for itr_ in 0..maxiter {
                itr = itr_;
                if q1 == q0 {
                    if p1 != p0 {
                        warn!("Tolerance of {} reached.", py_repr(p1 - p0));
                    }
                    p = (p1 + p0) / 2.0;
                    return results(p, funcalls, itr + 1, false, method);
                } else if q1.abs() > q0.abs() {
                    p = (-q0 / q1 * p1 + p0) / (1. - q0 / q1);
                } else {
                    p = (-q1 / q0 * p0 + p1) / (1. - q1 / q0);
                }
                if np_isclose(p, p1, rtol, tol) {
                    return results(p, funcalls, itr + 1, true, method);
                }
                p0 = p1;
                q0 = q1;
                p1 = p;
                q1 = func(p1)?;
                funcalls += 1;
            }
        }

        results(p, funcalls, itr + 1, false, method)
    }

    // TOMS "Algorithm 748: Enclosing Zeros of Continuous Functions", by
    //  Alefeld, G. E. and Potra, F. A. and Shi, Yixun

    /// Port of scipy.optimize.toms748 (with full_output=True, disp=False)
    fn toms748(
        f: &ScalarFn,
        a: f64,
        b: f64,
        k: usize,
        xtol: f64,
        rtol: f64,
        maxiter: usize,
    ) -> Result<RawResult, Interrupt> {
        if xtol <= 0. {
            return Err(value_error(format!(
                "xtol too small ({} <= 0)",
                py_format_g(xtol)
            )));
        }
        if rtol < DEFAULT_RTOL / 4. {
            return Err(value_error(format!(
                "rtol too small ({} < {})",
                py_format_g(rtol),
                py_format_g(DEFAULT_RTOL / 4.)
            )));
        }
        if maxiter < 1 {
            return Err(value_error("maxiter must be greater than 0".into()));
        }
        if !a.is_finite() {
            return Err(value_error(format!("a is not finite {}", py_repr(a))));
        }
        if !b.is_finite() {
            return Err(value_error(format!("b is not finite {}", py_repr(b))));
        }
        if a >= b {
            return Err(value_error(format!(
                "a and b are not an interval [{}, {}]",
                py_repr(a),
                py_repr(b)
            )));
        }
        if k < 1 {
            return Err(value_error(format!("k too small ({k} < 1)")));
        }

        let f = NanRaise::new(f);
        let mut solver = Toms748Solver::new(&f, xtol, rtol, maxiter, k);
        solver.solve(a, b)
    }

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum Toms748Status {
        Converged,
        ConvErr,
        InProgress,
    }

    /// Solve f(x) == 0 using Algorithm748 of Alefeld, Potro & Shi.
    struct Toms748Solver<'a, 'f> {
        f: &'a NanRaise<'f>,
        function_calls: usize,
        iterations: usize,
        k: usize,
        // ab=[a,b] is a global interval containing a root
        ab: [f64; 2],
        // fab is function values at a, b
        fab: [f64; 2],
        d: Option<f64>,
        fd: Option<f64>,
        e: Option<f64>,
        fe: Option<f64>,
        xtol: f64,
        rtol: f64,
        maxiter: usize,
    }

    impl<'a, 'f> Toms748Solver<'a, 'f> {
        const MU: f64 = 0.5;
        const K_MIN: usize = 1;
        const K_MAX: usize = 100; // A very high value for real usage. Expect 1, 2, maybe 3.

        fn new(f: &'a NanRaise<'f>, xtol: f64, rtol: f64, maxiter: usize, k: usize) -> Self {
            // Silently replace a low value of k with 1
            let mut k = k.max(Self::K_MIN);
            // Noisily replace a high value of k with K_MAX
            if k > Self::K_MAX {
                warn!("toms748: Overriding k: ->{}", Self::K_MAX);
                k = Self::K_MAX;
            }

            Self {
                f,
                function_calls: 0,
                iterations: 0,
                k,
                ab: [f64::NAN, f64::NAN],
                fab: [f64::NAN, f64::NAN],
                d: None,
                fd: None,
                e: None,
                fe: None,
                xtol,
                rtol,
                maxiter,
            }
        }

        /// Call the user-supplied function, update book-keeping
        fn callf(&mut self, x: f64) -> Result<f64, Interrupt> {
            let fx = self.f.call(x)?;
            self.function_calls += 1;
            if !fx.is_finite() {
                return Err(value_error(format!(
                    "Invalid function value: f({x:.6}) -> {} ",
                    py_repr(fx)
                )));
            }
            Ok(fx)
        }

        /// Package the result and statistics
        fn get_result(&self, x: f64, converged: bool) -> RawResult {
            RawResult {
                root: x,
                function_calls: self.function_calls,
                iterations: self.iterations,
                converged,
            }
        }

        fn update_bracket(&mut self, c: f64, fc: f64) -> (f64, f64) {
            update_bracket(&mut self.ab, &mut self.fab, c, fc)
        }

        fn ab_midpoint(&self) -> f64 {
            (self.ab[0] + self.ab[1]) / 2.0
        }

        /// Prepare for the iterations.
        fn start(&mut self, a: f64, b: f64) -> Result<(Toms748Status, f64), Interrupt> {
            self.function_calls = 0;
            self.iterations = 0;

            self.ab = [a, b];
            if !a.is_finite() {
                return Err(value_error(format!("Invalid x value: {} ", py_repr(a))));
            }
            if !b.is_finite() {
                return Err(value_error(format!("Invalid x value: {} ", py_repr(b))));
            }

            let fa = self.callf(a)?;
            if fa == 0. {
                return Ok((Toms748Status::Converged, a));
            }
            let fb = self.callf(b)?;
            if fb == 0. {
                return Ok((Toms748Status::Converged, b));
            }

            if np_sign(fb) * np_sign(fa) > 0. {
                return Err(value_error(format!(
                    "f(a) and f(b) must have different signs, but f({})={}, f({})={} ",
                    py_format_e(a),
                    py_format_e(fa),
                    py_format_e(b),
                    py_format_e(fb)
                )));
            }
            self.fab = [fa, fb];

            Ok((Toms748Status::InProgress, self.ab_midpoint()))
        }

        /// Determine the current status.
        fn get_status(&self) -> (Toms748Status, f64) {
            let [a, b] = self.ab;
            if np_isclose(a, b, self.rtol, self.xtol) {
                return (Toms748Status::Converged, self.ab_midpoint());
            }
            if self.iterations >= self.maxiter {
                return (Toms748Status::ConvErr, self.ab_midpoint());
            }
            (Toms748Status::InProgress, self.ab_midpoint())
        }

        /// Perform one step in the algorithm.
        ///
        /// Implements Algorithm 4.1(k=1) or 4.2(k=2) in [APS1995]
        fn iterate(&mut self) -> Result<(Toms748Status, f64), Interrupt> {
            self.iterations += 1;
            let eps = f64::EPSILON;
            let (mut d, mut fd, mut e, mut fe) = (self.d, self.fd, self.e, self.fe);
            let ab_width = self.ab[1] - self.ab[0]; // Need the start width below
            let mut c = None;

            for nsteps in 2..(self.k + 2) {
                // If the f-values are sufficiently separated, perform an inverse
                // polynomial interpolation step. Otherwise, nsteps repeats of
                // an approximate Newton-Raphson step.
                if notclose(
                    &[Some(self.fab[0]), Some(self.fab[1]), fd, fe],
                    0.,
                    32. * eps,
                ) {
                    let c0 = inverse_poly_zero(
                        self.ab[0],
                        self.ab[1],
                        d.expect("d is set when f-values are not close"),
                        e.expect("e is set when f-values are not close"),
                        self.fab[0],
                        self.fab[1],
                        fd.expect("fd is set when f-values are not close"),
                        fe.expect("fe is set when f-values are not close"),
                    );
                    if self.ab[0] < c0 && c0 < self.ab[1] {
                        c = Some(c0);
                    }
                }
                // NB. as in scipy, c is not reset between steps
                let c_val = match c {
                    Some(c) => c,
                    None => {
                        let c_new = newton_quadratic(
                            self.ab,
                            self.fab,
                            d.expect("d is always set after the first step"),
                            fd.expect("fd is always set after the first step"),
                            nsteps,
                        );
                        c = Some(c_new);
                        c_new
                    }
                };

                let fc = self.callf(c_val)?;
                if fc == 0. {
                    return Ok((Toms748Status::Converged, c_val));
                }

                // re-bracket
                e = d;
                fe = fd;
                let (new_d, new_fd) = self.update_bracket(c_val, fc);
                d = Some(new_d);
                fd = Some(new_fd);
            }

            // u is the endpoint with the smallest f-value
            let uix = if self.fab[0].abs() < self.fab[1].abs() {
                0
            } else {
                1
            };
            let (u, fu) = (self.ab[uix], self.fab[uix]);

            let a_dd = compute_divided_differences(&self.ab, &self.fab, uix == 0)[1];
            let mut c = u - 2. * fu / a_dd;
            if (c - u).abs() > 0.5 * (self.ab[1] - self.ab[0]) {
                c = self.ab_midpoint();
            } else if np_isclose(c, u, eps, 0.) {
                // c didn't change (much).
                // Either because the f-values at the endpoints have vastly
                // differing magnitudes, or because the root is very close to
                // that endpoint
                let frs = [frexp_exponent(self.fab[0]), frexp_exponent(self.fab[1])];
                if frs[uix] < frs[1 - uix] - 50 {
                    // Differ by more than 2**50
                    c = (31. * self.ab[uix] + self.ab[1 - uix]) / 32.;
                } else {
                    // Make a bigger adjustment, about the
                    // size of the requested tolerance.
                    let mm = if uix == 0 { 1. } else { -1. };
                    let adj = mm * c.abs() * self.rtol + mm * self.xtol;
                    c = u + adj;
                }
                if !(self.ab[0] < c && c < self.ab[1]) {
                    c = self.ab_midpoint();
                }
            }

            let fc = self.callf(c)?;
            if fc == 0. {
                return Ok((Toms748Status::Converged, c));
            }

            e = d;
            fe = fd;
            let (new_d, new_fd) = self.update_bracket(c, fc);
            d = Some(new_d);
            fd = Some(new_fd);

            // If the width of the new interval did not decrease enough, bisect
            if self.ab[1] - self.ab[0] > Self::MU * ab_width {
                e = d;
                fe = fd;
                let z = self.ab_midpoint();
                let fz = self.callf(z)?;
                if fz == 0. {
                    return Ok((Toms748Status::Converged, z));
                }
                let (new_d, new_fd) = self.update_bracket(z, fz);
                d = Some(new_d);
                fd = Some(new_fd);
            }

            // Record d and e for next iteration
            self.d = d;
            self.fd = fd;
            self.e = e;
            self.fe = fe;

            Ok(self.get_status())
        }

        /// Solve f(x) = 0 given an interval containing a root.
        fn solve(&mut self, a: f64, b: f64) -> Result<RawResult, Interrupt> {
            let (status, xn) = self.start(a, b)?;
            if status == Toms748Status::Converged {
                return Ok(self.get_result(xn, true));
            }

            // The first step only has two x-values.
            let mut c = secant(self.ab, self.fab);
            if !(self.ab[0] < c && c < self.ab[1]) {
                c = self.ab_midpoint();
            }
            let fc = self.callf(c)?;
            if fc == 0. {
                return Ok(self.get_result(c, true));
            }

            let (d, fd) = self.update_bracket(c, fc);
            self.d = Some(d);
            self.fd = Some(fd);
            self.e = None;
            self.fe = None;
            self.iterations += 1;

            loop {
                let (status, xn) = self.iterate()?;
                match status {
                    Toms748Status::Converged => return Ok(self.get_result(xn, true)),
                    Toms748Status::ConvErr => return Ok(self.get_result(xn, false)),
                    Toms748Status::InProgress => {}
                }
            }
        }
    }

    /// Ensure not None, not 0, all finite, and not very close to each other
    fn notclose(fs: &[Option<f64>], rtol: f64, atol: f64) -> bool {
        let Some(fs) = fs.iter().copied().collect::<Option<Vec<f64>>>() else {
            return false;
        };
        fs.iter().all(|&f| f != 0.)
            && fs.iter().all(|f| f.is_finite())
            && !(0..fs.len().saturating_sub(1)).any(|i| {
                fs[i + 1..]
                    .iter()
                    .any(|&f| np_isclose(fs[i], f, rtol, atol))
            })
    }

    /// Perform a secant step, taking a little care
    fn secant(xvals: [f64; 2], fvals: [f64; 2]) -> f64 {
        // Secant has many "mathematically" equivalent formulations
        // x2 = x0 - (x1 - x0)/(f1 - f0) * f0
        //    = x1 - (x1 - x0)/(f1 - f0) * f1
        //    = (-x1 * f0 + x0 * f1) / (f1 - f0)
        //    = (-f0 / f1 * x1 + x0) / (1 - f0 / f1)
        //    = (-f1 / f0 * x0 + x1) / (1 - f1 / f0)
        let [x0, x1] = xvals;
        let [f0, f1] = fvals;
        if f0 == f1 {
            return f64::NAN;
        }
        if f1.abs() > f0.abs() {
            (-f0 / f1 * x1 + x0) / (1. - f0 / f1)
        } else {
            (-f1 / f0 * x0 + x1) / (1. - f1 / f0)
        }
    }

    /// Update a bracket given (c, fc), return the discarded endpoints.
    fn update_bracket(ab: &mut [f64; 2], fab: &mut [f64; 2], c: f64, fc: f64) -> (f64, f64) {
        let [fa, _fb] = *fab;
        let idx = if np_sign(fa) * np_sign(fc) > 0. { 0 } else { 1 };
        let (rx, rfx) = (ab[idx], fab[idx]);
        fab[idx] = fc;
        ab[idx] = c;
        (rx, rfx)
    }

    /// Return the main diagonal (or last row) of the matrix of divided differences for the
    /// xvals, fvals pairs: f[a], f[a, b] and f[a, b, c].
    /// If forward is false, return f[c], f[b, c], f[a, b, c].
    ///
    /// (Equivalent of scipy's _compute_divided_differences with full=False.)
    fn compute_divided_differences(xvals: &[f64], fvals: &[f64], forward: bool) -> Vec<f64> {
        let mut dd = fvals.to_vec();
        let mut row = fvals.to_vec();
        dd[0] = if forward {
            fvals[0]
        } else {
            fvals[fvals.len() - 1]
        };
        for i in 1..xvals.len() {
            let len_row = row.len();
            row = (0..len_row - 1)
                .map(|j| (row[j + 1] - row[j]) / (xvals[i + j] - xvals[j]))
                .collect();
            dd[i] = if forward { row[0] } else { row[row.len() - 1] };
        }
        dd
    }

    /// Compute p(x) for the polynomial passing through the specified locations.
    ///
    /// Use Neville's algorithm to compute p(x) where p is the minimal degree
    /// polynomial passing through the points xvals, fvals
    fn interpolated_poly(xvals: &[f64], fvals: &[f64], x: f64) -> f64 {
        let n = xvals.len();
        let mut q = vec![vec![0.; n]; n];
        let mut d = vec![vec![0.; n]; n];
        for i in 0..n {
            q[i][0] = fvals[i];
            d[i][0] = fvals[i];
        }
        for k in 1..n {
            for j in 0..(n - k) {
                let alpha = d[k + j][k - 1] - q[k - 1 + j][k - 1];
                let diffik = xvals[j] - xvals[k + j];
                q[k + j][k] = (xvals[k + j] - x) / diffik * alpha;
                d[k + j][k] = (xvals[j] - x) / diffik * alpha;
            }
        }
        // Expect Q[-1, 1:] to be small relative to Q[-1, 0] as x approaches a root
        q[n - 1][1..].iter().sum::<f64>() + q[n - 1][0]
    }

    /// Inverse cubic interpolation f-values -> x-values
    ///
    /// Given four points (fa, a), (fb, b), (fc, c), (fd, d) with
    /// fa, fb, fc, fd all distinct, find poly IP(y) through the 4 points
    /// and compute x=IP(0).
    #[allow(clippy::too_many_arguments)]
    fn inverse_poly_zero(
        a: f64,
        b: f64,
        c: f64,
        d: f64,
        fa: f64,
        fb: f64,
        fc: f64,
        fd: f64,
    ) -> f64 {
        interpolated_poly(&[fa, fb, fc, fd], &[a, b, c, d], 0.)
    }

    /// Apply Newton-Raphson like steps, using divided differences to approximate f'
    ///
    /// ab is a real interval [a, b] containing a root,
    /// fab holds the real values of f(a), f(b)
    /// d is a real number outside [ab, b]
    /// k is the number of steps to apply
    fn newton_quadratic(ab: [f64; 2], fab: [f64; 2], d: f64, fd: f64, k: usize) -> f64 {
        let [a, b] = ab;
        let [fa, fb] = fab;
        let dd = compute_divided_differences(&[a, b, d], &[fa, fb, fd], true);
        let (big_b, big_a) = (dd[1], dd[2]);

        // p is the quadratic polynomial through the 3 points
        let p = |x: f64| -> f64 {
            // Horner evaluation of fa + B * (x - a) + A * (x - a) * (x - b)
            (big_a * (x - b) + big_b) * (x - a) + fa
        };

        let mut r;
        if big_a == 0. {
            r = a - fa / big_b;
        } else {
            r = if np_sign(big_a) * np_sign(fa) > 0. {
                a
            } else {
                b
            };
            // Apply k Newton-Raphson steps to p(x), starting from x=r
            for _ in 0..k {
                let r1 = r - p(r) / (big_b + big_a * (2. * r - a - b));
                if !(ab[0] < r1 && r1 < ab[1]) {
                    if ab[0] < r && r < ab[1] {
                        return r;
                    }
                    r = (ab[0] + ab[1]) / 2.0;
                    break;
                }
                r = r1;
            }
        }

        r
    }

    /// Equivalent of np.isclose for scalars (with equal_nan=False)
    fn np_isclose(a: f64, b: f64, rtol: f64, atol: f64) -> bool {
        ((a - b).abs() <= atol + rtol * b.abs() && b.is_finite()) || a == b
    }

    /// Equivalent of np.sign for a scalar
    fn np_sign(x: f64) -> f64 {
        if x > 0. {
            1.
        } else if x < 0. {
            -1.
        } else if x == 0. {
            0.
        } else {
            f64::NAN
        }
    }

    /// Equivalent of np.maximum for scalars (which propagates NaNs)
    fn np_maximum(a: f64, b: f64) -> f64 {
        if a.is_nan() || b.is_nan() {
            f64::NAN
        } else {
            a.max(b)
        }
    }

    /// The exponent from np.frexp, i.e. e such that x = m * 2**e with 0.5 <= |m| < 1
    fn frexp_exponent(x: f64) -> i32 {
        if x == 0. || !x.is_finite() {
            return 0;
        }
        let bits = x.to_bits();
        let biased_exponent = ((bits >> 52) & 0x7ff) as i32;
        if biased_exponent == 0 {
            // subnormal - normalise first
            return frexp_exponent(x * 2f64.powi(54)) - 54;
        }
        biased_exponent - 1022
    }

    /// Format a float as Python's repr() would
    fn py_repr(x: f64) -> String {
        if x.is_nan() {
            return "nan".into();
        }
        if x.is_infinite() {
            return if x > 0. { "inf" } else { "-inf" }.into();
        }
        // Rust's {:e} formatting gives the shortest round-trip digits, as does Python's repr
        let sci = format!("{:e}", x.abs());
        let (mantissa, exponent) = sci.split_once('e').unwrap();
        let exponent: i32 = exponent.parse().unwrap();
        let digits = mantissa.replace('.', "");
        let sign = if x.is_sign_negative() { "-" } else { "" };

        if (-4..16).contains(&exponent) {
            let formatted = if exponent >= 0 {
                let int_len = (exponent + 1) as usize;
                if digits.len() <= int_len {
                    format!("{digits}{}.0", "0".repeat(int_len - digits.len()))
                } else {
                    format!("{}.{}", &digits[..int_len], &digits[int_len..])
                }
            } else {
                format!("0.{}{digits}", "0".repeat((-exponent - 1) as usize))
            };
            format!("{sign}{formatted}")
        } else {
            let exp_sign = if exponent < 0 { "-" } else { "+" };
            format!("{sign}{mantissa}e{exp_sign}{:02}", exponent.abs())
        }
    }

    /// Format a float as Python's "%g" would
    fn py_format_g(x: f64) -> String {
        if !x.is_finite() {
            return py_repr(x);
        }
        let precision = 6;
        let sci = format!("{:.*e}", precision - 1, x);
        let (mantissa, exponent) = sci.split_once('e').unwrap();
        let exponent: i32 = exponent.parse().unwrap();
        let strip = |s: &str| -> String {
            if s.contains('.') {
                s.trim_end_matches('0').trim_end_matches('.').to_string()
            } else {
                s.to_string()
            }
        };
        if (-4..precision as i32).contains(&exponent) {
            strip(&format!(
                "{:.*}",
                (precision as i32 - 1 - exponent) as usize,
                x
            ))
        } else {
            let exp_sign = if exponent < 0 { "-" } else { "+" };
            format!("{}e{exp_sign}{:02}", strip(mantissa), exponent.abs())
        }
    }

    /// Format a float as Python's "{:e}" would
    fn py_format_e(x: f64) -> String {
        if !x.is_finite() {
            return py_repr(x);
        }
        let sci = format!("{:.6e}", x);
        let (mantissa, exponent) = sci.split_once('e').unwrap();
        let exponent: i32 = exponent.parse().unwrap();
        let exp_sign = if exponent < 0 { "-" } else { "+" };
        format!("{mantissa}e{exp_sign}{:02}", exponent.abs())
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        enum Expected {
            Ok {
                root: f64,
                iterations: Option<usize>,
                // scipy's C solvers leave iterations uninitialised when an endpoint is a root
                check_iterations: bool,
                function_calls: usize,
                converged: bool,
                flag: &'static str,
                method: Method,
            },
            ValueErr(&'static str),
            TypeErr(&'static str),
        }

        fn check(
            description: &str,
            actual: Result<RootResults, RootScalarError>,
            expected: Expected,
        ) -> Option<String> {
            let mismatch = |details: String| Some(format!("{description}: {details}"));
            match (actual, expected) {
                (
                    Ok(actual),
                    Expected::Ok {
                        root,
                        iterations,
                        check_iterations,
                        function_calls,
                        converged,
                        flag,
                        method,
                    },
                ) => {
                    let root_matches = actual.root.to_bits() == root.to_bits()
                        || (actual.root.is_nan() && root.is_nan());
                    if !root_matches
                        || (check_iterations && actual.iterations != iterations)
                        || actual.function_calls != function_calls
                        || actual.converged != converged
                        || actual.flag != flag
                        || actual.method != method
                    {
                        return mismatch(format!(
                            "expected root={root:?} iterations={iterations:?} function_calls={function_calls} converged={converged} flag={flag:?} method={method}, got {actual:?}"
                        ));
                    }
                    None
                }
                (Err(RootScalarError::Value(actual)), Expected::ValueErr(expected))
                | (Err(RootScalarError::Type(actual)), Expected::TypeErr(expected))
                    if actual == expected =>
                {
                    None
                }
                (Ok(actual), _) => mismatch(format!("expected an error, got {actual:?}")),
                (Err(actual), _) => mismatch(format!("unexpected error {actual}")),
            }
        }

        #[test]
        fn test_root_scalar_against_scipy() {
            // Expected values generated by running scipy.optimize.root_scalar from scipy 1.14.1
            // on equivalent Python functions
            let cases: Vec<(&str, Result<RootResults, RootScalarError>, Expected)> = vec![
                (
                    "cubic bisect bracket=(0, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Bisect),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0000000000004547_f64,
                        iterations: Some(41),
                        check_iterations: true,
                        function_calls: 43,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Bisect,
                    },
                ),
                (
                    "cubic brentq bracket=(0, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Brentq),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(10),
                        check_iterations: true,
                        function_calls: 11,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Brentq,
                    },
                ),
                (
                    "cubic brenth bracket=(0, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Brenth),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(10),
                        check_iterations: true,
                        function_calls: 11,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Brenth,
                    },
                ),
                (
                    "cubic ridder bracket=(0, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Ridder),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.000000000000369_f64,
                        iterations: Some(6),
                        check_iterations: true,
                        function_calls: 14,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Ridder,
                    },
                ),
                (
                    "cubic toms748 bracket=(0, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(5),
                        check_iterations: true,
                        function_calls: 11,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "cubic toms748 bracket=(0, 3) k=2",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            k: Some(2),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 12,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "cubic toms748 bracket=(0, 3) k=3",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            k: Some(3),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 13,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "cubic bisect bracket=(0, 3) xtol=1e-06 rtol=1e-10",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Bisect),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            xtol: Some(1.0e-06_f64),
                            rtol: Some(1.0e-10_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.9999997615814209_f64,
                        iterations: Some(22),
                        check_iterations: true,
                        function_calls: 24,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Bisect,
                    },
                ),
                (
                    "cubic brentq bracket=(0, 3) xtol=1e-06 rtol=1e-10",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Brentq),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            xtol: Some(1.0e-06_f64),
                            rtol: Some(1.0e-10_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.9999998220370444_f64,
                        iterations: Some(9),
                        check_iterations: true,
                        function_calls: 10,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Brentq,
                    },
                ),
                (
                    "cubic brenth bracket=(0, 3) xtol=1e-06 rtol=1e-10",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Brenth),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            xtol: Some(1.0e-06_f64),
                            rtol: Some(1.0e-10_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.000000129284141_f64,
                        iterations: Some(9),
                        check_iterations: true,
                        function_calls: 10,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Brenth,
                    },
                ),
                (
                    "cubic ridder bracket=(0, 3) xtol=1e-06 rtol=1e-10",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Ridder),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            xtol: Some(1.0e-06_f64),
                            rtol: Some(1.0e-10_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0000004963614826_f64,
                        iterations: Some(5),
                        check_iterations: true,
                        function_calls: 12,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Ridder,
                    },
                ),
                (
                    "cubic toms748 bracket=(0, 3) xtol=1e-06 rtol=1e-10",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            xtol: Some(1.0e-06_f64),
                            rtol: Some(1.0e-10_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.9999999989165773_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 10,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "cosx bisect bracket=(0, 1)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Bisect),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332147577_f64,
                        iterations: Some(39),
                        check_iterations: true,
                        function_calls: 41,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Bisect,
                    },
                ),
                (
                    "cosx brentq bracket=(0, 1)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Brentq),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151559_f64,
                        iterations: Some(7),
                        check_iterations: true,
                        function_calls: 8,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Brentq,
                    },
                ),
                (
                    "cosx brenth bracket=(0, 1)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Brenth),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.739085133215158_f64,
                        iterations: Some(7),
                        check_iterations: true,
                        function_calls: 8,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Brenth,
                    },
                ),
                (
                    "cosx ridder bracket=(0, 1)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Ridder),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151607_f64,
                        iterations: Some(5),
                        check_iterations: true,
                        function_calls: 12,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Ridder,
                    },
                ),
                (
                    "cosx toms748 bracket=(0, 1)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151607_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 8,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "cosx toms748 bracket=(0, 1) k=2",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            k: Some(2),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151607_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 7,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "cosx toms748 bracket=(0, 1) k=3",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            k: Some(3),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151607_f64,
                        iterations: Some(2),
                        check_iterations: true,
                        function_calls: 7,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "cosx bisect bracket=(0, 1) xtol=1e-06 rtol=1e-10",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Bisect),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            xtol: Some(1.0e-06_f64),
                            rtol: Some(1.0e-10_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851974487305_f64,
                        iterations: Some(20),
                        check_iterations: true,
                        function_calls: 22,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Bisect,
                    },
                ),
                (
                    "cosx brentq bracket=(0, 1) xtol=1e-06 rtol=1e-10",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Brentq),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            xtol: Some(1.0e-06_f64),
                            rtol: Some(1.0e-10_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851295647751_f64,
                        iterations: Some(6),
                        check_iterations: true,
                        function_calls: 7,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Brentq,
                    },
                ),
                (
                    "cosx brenth bracket=(0, 1) xtol=1e-06 rtol=1e-10",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Brenth),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            xtol: Some(1.0e-06_f64),
                            rtol: Some(1.0e-10_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851305146882_f64,
                        iterations: Some(6),
                        check_iterations: true,
                        function_calls: 7,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Brenth,
                    },
                ),
                (
                    "cosx ridder bracket=(0, 1) xtol=1e-06 rtol=1e-10",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Ridder),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            xtol: Some(1.0e-06_f64),
                            rtol: Some(1.0e-10_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390846399335982_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 10,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Ridder,
                    },
                ),
                (
                    "cosx toms748 bracket=(0, 1) xtol=1e-06 rtol=1e-10",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            xtol: Some(1.0e-06_f64),
                            rtol: Some(1.0e-10_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332150754_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 7,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "triple bisect bracket=(0, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok((x - 1.) * (x - 1.) * (x - 1.)) },
                        RootScalarOptions {
                            method: Some(Method::Bisect),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0000000000004547_f64,
                        iterations: Some(41),
                        check_iterations: true,
                        function_calls: 43,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Bisect,
                    },
                ),
                (
                    "triple brentq bracket=(0, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok((x - 1.) * (x - 1.) * (x - 1.)) },
                        RootScalarOptions {
                            method: Some(Method::Brentq),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0000000002456189_f64,
                        iterations: Some(100),
                        check_iterations: true,
                        function_calls: 102,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Brentq,
                    },
                ),
                (
                    "triple brenth bracket=(0, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok((x - 1.) * (x - 1.) * (x - 1.)) },
                        RootScalarOptions {
                            method: Some(Method::Brenth),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0000000002139096_f64,
                        iterations: Some(100),
                        check_iterations: true,
                        function_calls: 102,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Brenth,
                    },
                ),
                (
                    "triple ridder bracket=(0, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok((x - 1.) * (x - 1.) * (x - 1.)) },
                        RootScalarOptions {
                            method: Some(Method::Ridder),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.9999999999995511_f64,
                        iterations: Some(34),
                        check_iterations: true,
                        function_calls: 70,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Ridder,
                    },
                ),
                (
                    "triple toms748 bracket=(0, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok((x - 1.) * (x - 1.) * (x - 1.)) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0000000000001854_f64,
                        iterations: Some(31),
                        check_iterations: true,
                        function_calls: 89,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "triple toms748 bracket=(0, 3) k=2",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok((x - 1.) * (x - 1.) * (x - 1.)) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            k: Some(2),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0000000000004496_f64,
                        iterations: Some(32),
                        check_iterations: true,
                        function_calls: 121,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "triple toms748 bracket=(0, 3) k=3",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok((x - 1.) * (x - 1.) * (x - 1.)) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            k: Some(3),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0000000000001434_f64,
                        iterations: Some(31),
                        check_iterations: true,
                        function_calls: 151,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "kepler bisect bracket=(0, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Bisect),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745907_f64,
                        iterations: Some(41),
                        check_iterations: true,
                        function_calls: 43,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Bisect,
                    },
                ),
                (
                    "kepler brentq bracket=(0, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Brentq),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745114_f64,
                        iterations: Some(9),
                        check_iterations: true,
                        function_calls: 10,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Brentq,
                    },
                ),
                (
                    "kepler brenth bracket=(0, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Brenth),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745323_f64,
                        iterations: Some(9),
                        check_iterations: true,
                        function_calls: 10,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Brenth,
                    },
                ),
                (
                    "kepler ridder bracket=(0, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Ridder),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868735368_f64,
                        iterations: Some(6),
                        check_iterations: true,
                        function_calls: 14,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Ridder,
                    },
                ),
                (
                    "kepler toms748 bracket=(0, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745323_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 8,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "kepler toms748 bracket=(0, 3) k=2",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            k: Some(2),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745323_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 9,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "kepler toms748 bracket=(0, 3) k=3",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            k: Some(3),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745323_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 10,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "kepler bisect bracket=(0, 3) xtol=1e-06 rtol=1e-10",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Bisect),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            xtol: Some(1.0e-06_f64),
                            rtol: Some(1.0e-10_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620860576629639_f64,
                        iterations: Some(22),
                        check_iterations: true,
                        function_calls: 24,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Bisect,
                    },
                ),
                (
                    "kepler brentq bracket=(0, 3) xtol=1e-06 rtol=1e-10",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Brentq),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            xtol: Some(1.0e-06_f64),
                            rtol: Some(1.0e-10_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866788187573_f64,
                        iterations: Some(8),
                        check_iterations: true,
                        function_calls: 9,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Brentq,
                    },
                ),
                (
                    "kepler brenth bracket=(0, 3) xtol=1e-06 rtol=1e-10",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Brenth),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            xtol: Some(1.0e-06_f64),
                            rtol: Some(1.0e-10_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620868679023552_f64,
                        iterations: Some(8),
                        check_iterations: true,
                        function_calls: 9,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Brenth,
                    },
                ),
                (
                    "kepler ridder bracket=(0, 3) xtol=1e-06 rtol=1e-10",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Ridder),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            xtol: Some(1.0e-06_f64),
                            rtol: Some(1.0e-10_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620864897324472_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 10,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Ridder,
                    },
                ),
                (
                    "kepler toms748 bracket=(0, 3) xtol=1e-06 rtol=1e-10",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            xtol: Some(1.0e-06_f64),
                            rtol: Some(1.0e-10_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745323_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 8,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "cosx bisect bracket=(0, 1) maxiter=3",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Bisect),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            maxiter: Some(3),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.625_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 5,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Bisect,
                    },
                ),
                (
                    "nanf bisect bracket=(0, 1)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> {
                            Ok(if x > 0.5 { f64::NAN } else { x - 0.7 })
                        },
                        RootScalarOptions {
                            method: Some(Method::Bisect),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: None,
                        check_iterations: true,
                        function_calls: 2,
                        converged: false,
                        flag: r#"The function value at x=1.0 is NaN; solver cannot continue."#,
                        method: Method::Bisect,
                    },
                ),
                (
                    "inff bisect bracket=(0, 1)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> {
                            Ok(if x > 0.5 { f64::INFINITY } else { x - 0.7 })
                        },
                        RootScalarOptions {
                            method: Some(Method::Bisect),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.500000000001819_f64,
                        iterations: Some(39),
                        check_iterations: true,
                        function_calls: 41,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Bisect,
                    },
                ),
                (
                    "cubic fp bisect bracket=(0, 3) fprime=True",
                    root_scalar(
                        |x: f64| -> anyhow::Result<(f64, f64)> { Ok((x * x * x - 1., 3. * x * x)) },
                        RootScalarOptions {
                            method: Some(Method::Bisect),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            fprime: Some(Derivative::FromObjective),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0000000000004547_f64,
                        iterations: Some(41),
                        check_iterations: true,
                        function_calls: 43,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Bisect,
                    },
                ),
                (
                    "cosx brentq bracket=(0, 1) maxiter=3",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Brentq),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            maxiter: Some(3),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390910582853698_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 5,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Brentq,
                    },
                ),
                (
                    "nanf brentq bracket=(0, 1)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> {
                            Ok(if x > 0.5 { f64::NAN } else { x - 0.7 })
                        },
                        RootScalarOptions {
                            method: Some(Method::Brentq),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: None,
                        check_iterations: true,
                        function_calls: 2,
                        converged: false,
                        flag: r#"The function value at x=1.0 is NaN; solver cannot continue."#,
                        method: Method::Brentq,
                    },
                ),
                (
                    "inff brentq bracket=(0, 1)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> {
                            Ok(if x > 0.5 { f64::INFINITY } else { x - 0.7 })
                        },
                        RootScalarOptions {
                            method: Some(Method::Brentq),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.49999999999868106_f64,
                        iterations: Some(42),
                        check_iterations: true,
                        function_calls: 43,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Brentq,
                    },
                ),
                (
                    "cubic fp brentq bracket=(0, 3) fprime=True",
                    root_scalar(
                        |x: f64| -> anyhow::Result<(f64, f64)> { Ok((x * x * x - 1., 3. * x * x)) },
                        RootScalarOptions {
                            method: Some(Method::Brentq),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            fprime: Some(Derivative::FromObjective),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(10),
                        check_iterations: true,
                        function_calls: 11,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Brentq,
                    },
                ),
                (
                    "cosx brenth bracket=(0, 1) maxiter=3",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Brenth),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            maxiter: Some(3),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390895164441823_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 5,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Brenth,
                    },
                ),
                (
                    "nanf brenth bracket=(0, 1)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> {
                            Ok(if x > 0.5 { f64::NAN } else { x - 0.7 })
                        },
                        RootScalarOptions {
                            method: Some(Method::Brenth),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: None,
                        check_iterations: true,
                        function_calls: 2,
                        converged: false,
                        flag: r#"The function value at x=1.0 is NaN; solver cannot continue."#,
                        method: Method::Brenth,
                    },
                ),
                (
                    "inff brenth bracket=(0, 1)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> {
                            Ok(if x > 0.5 { f64::INFINITY } else { x - 0.7 })
                        },
                        RootScalarOptions {
                            method: Some(Method::Brenth),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.49999999999868106_f64,
                        iterations: Some(42),
                        check_iterations: true,
                        function_calls: 43,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Brenth,
                    },
                ),
                (
                    "cubic fp brenth bracket=(0, 3) fprime=True",
                    root_scalar(
                        |x: f64| -> anyhow::Result<(f64, f64)> { Ok((x * x * x - 1., 3. * x * x)) },
                        RootScalarOptions {
                            method: Some(Method::Brenth),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            fprime: Some(Derivative::FromObjective),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(10),
                        check_iterations: true,
                        function_calls: 11,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Brenth,
                    },
                ),
                (
                    "cosx ridder bracket=(0, 1) maxiter=3",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Ridder),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            maxiter: Some(3),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851399705524_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 8,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Ridder,
                    },
                ),
                (
                    "nanf ridder bracket=(0, 1)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> {
                            Ok(if x > 0.5 { f64::NAN } else { x - 0.7 })
                        },
                        RootScalarOptions {
                            method: Some(Method::Ridder),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: None,
                        check_iterations: true,
                        function_calls: 2,
                        converged: false,
                        flag: r#"The function value at x=1.0 is NaN; solver cannot continue."#,
                        method: Method::Ridder,
                    },
                ),
                (
                    "inff ridder bracket=(0, 1)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> {
                            Ok(if x > 0.5 { f64::INFINITY } else { x - 0.7 })
                        },
                        RootScalarOptions {
                            method: Some(Method::Ridder),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.9999999999009478_f64,
                        iterations: Some(100),
                        check_iterations: true,
                        function_calls: 202,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Ridder,
                    },
                ),
                (
                    "cubic fp ridder bracket=(0, 3) fprime=True",
                    root_scalar(
                        |x: f64| -> anyhow::Result<(f64, f64)> { Ok((x * x * x - 1., 3. * x * x)) },
                        RootScalarOptions {
                            method: Some(Method::Ridder),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            fprime: Some(Derivative::FromObjective),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.000000000000369_f64,
                        iterations: Some(6),
                        check_iterations: true,
                        function_calls: 14,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Ridder,
                    },
                ),
                (
                    "cosx toms748 bracket=(0, 1) maxiter=3",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            maxiter: Some(3),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332150754_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 7,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "nanf toms748 bracket=(0, 1)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> {
                            Ok(if x > 0.5 { f64::NAN } else { x - 0.7 })
                        },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: None,
                        check_iterations: true,
                        function_calls: 2,
                        converged: false,
                        flag: r#"The function value at x=1.0 is NaN; solver cannot continue."#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "inff toms748 bracket=(0, 1)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> {
                            Ok(if x > 0.5 { f64::INFINITY } else { x - 0.7 })
                        },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"Invalid function value: f(1.000000) -> inf "#),
                ),
                (
                    "cubic fp toms748 bracket=(0, 3) fprime=True",
                    root_scalar(
                        |x: f64| -> anyhow::Result<(f64, f64)> { Ok((x * x * x - 1., 3. * x * x)) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            fprime: Some(Derivative::FromObjective),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(5),
                        check_iterations: true,
                        function_calls: 11,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "cosx bisect bracket=(0, 1) maxiter=0",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Bisect),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            maxiter: Some(0),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.0_f64,
                        iterations: Some(0),
                        check_iterations: true,
                        function_calls: 2,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Bisect,
                    },
                ),
                (
                    "cosx toms748 bracket=(0, 1) k=150",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            k: Some(150),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151607_f64,
                        iterations: Some(2),
                        check_iterations: true,
                        function_calls: 7,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "cubic bracket=(0, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            bracket: Some((0.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(10),
                        check_iterations: true,
                        function_calls: 11,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Brentq,
                    },
                ),
                (
                    "ident toms748 bracket=(0, 1)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.0_f64,
                        iterations: Some(0),
                        check_iterations: true,
                        function_calls: 1,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "shifted toms748 bracket=(0, 1)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 1.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(0),
                        check_iterations: true,
                        function_calls: 2,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Toms748,
                    },
                ),
                (
                    "cubic newton x0=0.2 fprime=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(0.2_f64),
                            fprime: Some(Derivative::function(|x| Ok(3. * x * x))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(11),
                        check_iterations: true,
                        function_calls: 22,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "cubic newton x0=0.2",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(0.2_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(11),
                        check_iterations: true,
                        function_calls: 22,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "cubic secant x0=0.2",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Secant),
                            x0: Some(0.2_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0000000000000002_f64,
                        iterations: Some(22),
                        check_iterations: true,
                        function_calls: 23,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Secant,
                    },
                ),
                (
                    "cubic secant x0=0.2 x1=0.7",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Secant),
                            x0: Some(0.2_f64),
                            x1: Some(0.7_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(11),
                        check_iterations: true,
                        function_calls: 12,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Secant,
                    },
                ),
                (
                    "cubic halley x0=0.2 fprime=fn fprime2=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(0.2_f64),
                            fprime: Some(Derivative::function(|x| Ok(3. * x * x))),
                            fprime2: Some(Derivative::function(|x| Ok(6. * x))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(7),
                        check_iterations: true,
                        function_calls: 22,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "cubic fp newton x0=0.2 fprime=True",
                    root_scalar(
                        |x: f64| -> anyhow::Result<(f64, f64)> { Ok((x * x * x - 1., 3. * x * x)) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(0.2_f64),
                            fprime: Some(Derivative::FromObjective),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(11),
                        check_iterations: true,
                        function_calls: 11,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "cubic fpp halley x0=0.2 fprime=True fprime2=True",
                    root_scalar(
                        |x: f64| -> anyhow::Result<(f64, f64, f64)> {
                            Ok((x * x * x - 1., 3. * x * x, 6. * x))
                        },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(0.2_f64),
                            fprime: Some(Derivative::FromObjective),
                            fprime2: Some(Derivative::FromObjective),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(7),
                        check_iterations: true,
                        function_calls: 8,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "cubic fpp halley x0=0.2 fprime2=True",
                    root_scalar(
                        |x: f64| -> anyhow::Result<(f64, f64, f64)> {
                            Ok((x * x * x - 1., 3. * x * x, 6. * x))
                        },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(0.2_f64),
                            fprime2: Some(Derivative::FromObjective),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(7),
                        check_iterations: true,
                        function_calls: 8,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "cubic fp halley x0=0.2 fprime=True fprime2=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<(f64, f64)> { Ok((x * x * x - 1., 3. * x * x)) },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(0.2_f64),
                            fprime: Some(Derivative::FromObjective),
                            fprime2: Some(Derivative::function(|x| Ok(6. * x))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(7),
                        check_iterations: true,
                        function_calls: 8,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "cubic x0=0.2 fprime=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            x0: Some(0.2_f64),
                            fprime: Some(Derivative::function(|x| Ok(3. * x * x))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(11),
                        check_iterations: true,
                        function_calls: 22,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "cubic x0=0.2 fprime=fn fprime2=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            x0: Some(0.2_f64),
                            fprime: Some(Derivative::function(|x| Ok(3. * x * x))),
                            fprime2: Some(Derivative::function(|x| Ok(6. * x))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(7),
                        check_iterations: true,
                        function_calls: 22,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "cubic x0=0.2",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            x0: Some(0.2_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(11),
                        check_iterations: true,
                        function_calls: 22,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "cubic x0=0.2 x1=0.22000000000000003",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            x0: Some(0.2_f64),
                            x1: Some(0.22000000000000003_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(22),
                        check_iterations: true,
                        function_calls: 23,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Secant,
                    },
                ),
                (
                    "cubic newton x0=0.2 fprime=fn maxiter=3",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(0.2_f64),
                            maxiter: Some(3),
                            fprime: Some(Derivative::function(|x| Ok(3. * x * x))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 3.7765082728873756_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 6,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "cubic newton x0=0.2 fprime=fn rtol=0.001",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(0.2_f64),
                            rtol: Some(0.001_f64),
                            fprime: Some(Derivative::function(|x| Ok(3. * x * x))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0000000000575857_f64,
                        iterations: Some(10),
                        check_iterations: true,
                        function_calls: 20,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "cubic secant x0=0.2 xtol=0.001 maxiter=5",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Secant),
                            x0: Some(0.2_f64),
                            xtol: Some(0.001_f64),
                            maxiter: Some(5),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.24644402344342206_f64,
                        iterations: Some(5),
                        check_iterations: true,
                        function_calls: 7,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Secant,
                    },
                ),
                (
                    "cubic halley x0=0.2 fprime=fn fprime2=fn xtol=1e-12 rtol=1e-12",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(0.2_f64),
                            xtol: Some(1.0e-12_f64),
                            rtol: Some(1.0e-12_f64),
                            fprime: Some(Derivative::function(|x| Ok(3. * x * x))),
                            fprime2: Some(Derivative::function(|x| Ok(6. * x))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(7),
                        check_iterations: true,
                        function_calls: 22,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "cosx newton x0=1.0 fprime=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(1.0_f64),
                            fprime: Some(Derivative::function(|x| Ok(-x.sin() - 1.))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151607_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 8,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "cosx newton x0=1.0",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(1.0_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151607_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 8,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "cosx secant x0=1.0",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Secant),
                            x0: Some(1.0_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151608_f64,
                        iterations: Some(5),
                        check_iterations: true,
                        function_calls: 6,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Secant,
                    },
                ),
                (
                    "cosx secant x0=1.0 x1=1.5",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Secant),
                            x0: Some(1.0_f64),
                            x1: Some(1.5_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151827_f64,
                        iterations: Some(5),
                        check_iterations: true,
                        function_calls: 6,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Secant,
                    },
                ),
                (
                    "cosx halley x0=1.0 fprime=fn fprime2=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(1.0_f64),
                            fprime: Some(Derivative::function(|x| Ok(-x.sin() - 1.))),
                            fprime2: Some(Derivative::function(|x| Ok(-x.cos()))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151607_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 9,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "cosx fp newton x0=1.0 fprime=True",
                    root_scalar(
                        |x: f64| -> anyhow::Result<(f64, f64)> { Ok((x.cos() - x, -x.sin() - 1.)) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(1.0_f64),
                            fprime: Some(Derivative::FromObjective),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151607_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 4,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "cosx fpp halley x0=1.0 fprime=True fprime2=True",
                    root_scalar(
                        |x: f64| -> anyhow::Result<(f64, f64, f64)> {
                            Ok((x.cos() - x, -x.sin() - 1., -x.cos()))
                        },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(1.0_f64),
                            fprime: Some(Derivative::FromObjective),
                            fprime2: Some(Derivative::FromObjective),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151607_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 3,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "cosx fpp halley x0=1.0 fprime2=True",
                    root_scalar(
                        |x: f64| -> anyhow::Result<(f64, f64, f64)> {
                            Ok((x.cos() - x, -x.sin() - 1., -x.cos()))
                        },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(1.0_f64),
                            fprime2: Some(Derivative::FromObjective),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151607_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 3,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "cosx fp halley x0=1.0 fprime=True fprime2=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<(f64, f64)> { Ok((x.cos() - x, -x.sin() - 1.)) },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(1.0_f64),
                            fprime: Some(Derivative::FromObjective),
                            fprime2: Some(Derivative::function(|x| Ok(-x.cos()))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151607_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 3,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "cosx x0=1.0 fprime=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            x0: Some(1.0_f64),
                            fprime: Some(Derivative::function(|x| Ok(-x.sin() - 1.))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151607_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 8,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "cosx x0=1.0 fprime=fn fprime2=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            x0: Some(1.0_f64),
                            fprime: Some(Derivative::function(|x| Ok(-x.sin() - 1.))),
                            fprime2: Some(Derivative::function(|x| Ok(-x.cos()))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151607_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 9,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "cosx x0=1.0",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            x0: Some(1.0_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151607_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 8,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "cosx x0=1.0 x1=1.1",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            x0: Some(1.0_f64),
                            x1: Some(1.1_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151611_f64,
                        iterations: Some(5),
                        check_iterations: true,
                        function_calls: 6,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Secant,
                    },
                ),
                (
                    "cosx newton x0=1.0 fprime=fn maxiter=3",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(1.0_f64),
                            maxiter: Some(3),
                            fprime: Some(Derivative::function(|x| Ok(-x.sin() - 1.))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.739085133385284_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 6,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "cosx newton x0=1.0 fprime=fn rtol=0.001",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(1.0_f64),
                            rtol: Some(0.001_f64),
                            fprime: Some(Derivative::function(|x| Ok(-x.sin() - 1.))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.739085133385284_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 6,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "cosx secant x0=1.0 xtol=0.001 maxiter=5",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Secant),
                            x0: Some(1.0_f64),
                            xtol: Some(0.001_f64),
                            maxiter: Some(5),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390865256303388_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 4,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Secant,
                    },
                ),
                (
                    "cosx halley x0=1.0 fprime=fn fprime2=fn xtol=1e-12 rtol=1e-12",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x.cos() - x) },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(1.0_f64),
                            xtol: Some(1.0e-12_f64),
                            rtol: Some(1.0e-12_f64),
                            fprime: Some(Derivative::function(|x| Ok(-x.sin() - 1.))),
                            fprime2: Some(Derivative::function(|x| Ok(-x.cos()))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.7390851332151607_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 10,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "kepler newton x0=1.0 fprime=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(1.0_f64),
                            fprime: Some(Derivative::function(|x| Ok(1. - 0.9 * x.cos()))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745323_f64,
                        iterations: Some(6),
                        check_iterations: true,
                        function_calls: 12,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "kepler newton x0=1.0",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(1.0_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745323_f64,
                        iterations: Some(6),
                        check_iterations: true,
                        function_calls: 12,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "kepler secant x0=1.0",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Secant),
                            x0: Some(1.0_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745325_f64,
                        iterations: Some(8),
                        check_iterations: true,
                        function_calls: 9,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Secant,
                    },
                ),
                (
                    "kepler secant x0=1.0 x1=1.5",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Secant),
                            x0: Some(1.0_f64),
                            x1: Some(1.5_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.862086686874532_f64,
                        iterations: Some(7),
                        check_iterations: true,
                        function_calls: 8,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Secant,
                    },
                ),
                (
                    "kepler halley x0=1.0 fprime=fn fprime2=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(1.0_f64),
                            fprime: Some(Derivative::function(|x| Ok(1. - 0.9 * x.cos()))),
                            fprime2: Some(Derivative::function(|x| Ok(0.9 * x.sin()))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745323_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 13,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "kepler fp newton x0=1.0 fprime=True",
                    root_scalar(
                        |x: f64| -> anyhow::Result<(f64, f64)> {
                            Ok((x - 0.9 * x.sin() - 1., 1. - 0.9 * x.cos()))
                        },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(1.0_f64),
                            fprime: Some(Derivative::FromObjective),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745323_f64,
                        iterations: Some(6),
                        check_iterations: true,
                        function_calls: 6,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "kepler fpp halley x0=1.0 fprime=True fprime2=True",
                    root_scalar(
                        |x: f64| -> anyhow::Result<(f64, f64, f64)> {
                            Ok((x - 0.9 * x.sin() - 1., 1. - 0.9 * x.cos(), 0.9 * x.sin()))
                        },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(1.0_f64),
                            fprime: Some(Derivative::FromObjective),
                            fprime2: Some(Derivative::FromObjective),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745323_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 5,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "kepler fpp halley x0=1.0 fprime2=True",
                    root_scalar(
                        |x: f64| -> anyhow::Result<(f64, f64, f64)> {
                            Ok((x - 0.9 * x.sin() - 1., 1. - 0.9 * x.cos(), 0.9 * x.sin()))
                        },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(1.0_f64),
                            fprime2: Some(Derivative::FromObjective),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745323_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 5,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "kepler fp halley x0=1.0 fprime=True fprime2=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<(f64, f64)> {
                            Ok((x - 0.9 * x.sin() - 1., 1. - 0.9 * x.cos()))
                        },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(1.0_f64),
                            fprime: Some(Derivative::FromObjective),
                            fprime2: Some(Derivative::function(|x| Ok(0.9 * x.sin()))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745323_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 5,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "kepler x0=1.0 fprime=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            x0: Some(1.0_f64),
                            fprime: Some(Derivative::function(|x| Ok(1. - 0.9 * x.cos()))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745323_f64,
                        iterations: Some(6),
                        check_iterations: true,
                        function_calls: 12,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "kepler x0=1.0 fprime=fn fprime2=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            x0: Some(1.0_f64),
                            fprime: Some(Derivative::function(|x| Ok(1. - 0.9 * x.cos()))),
                            fprime2: Some(Derivative::function(|x| Ok(0.9 * x.sin()))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745323_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 13,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "kepler x0=1.0",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            x0: Some(1.0_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745323_f64,
                        iterations: Some(6),
                        check_iterations: true,
                        function_calls: 12,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "kepler x0=1.0 x1=1.1",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            x0: Some(1.0_f64),
                            x1: Some(1.1_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745716_f64,
                        iterations: Some(7),
                        check_iterations: true,
                        function_calls: 8,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Secant,
                    },
                ),
                (
                    "kepler newton x0=1.0 fprime=fn maxiter=3",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(1.0_f64),
                            maxiter: Some(3),
                            fprime: Some(Derivative::function(|x| Ok(1. - 0.9 * x.cos()))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8638802540484245_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 6,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "kepler newton x0=1.0 fprime=fn rtol=0.001",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(1.0_f64),
                            rtol: Some(0.001_f64),
                            fprime: Some(Derivative::function(|x| Ok(1. - 0.9 * x.cos()))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620877869546315_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 8,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "kepler secant x0=1.0 xtol=0.001 maxiter=5",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Secant),
                            x0: Some(1.0_f64),
                            xtol: Some(0.001_f64),
                            maxiter: Some(5),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620570720048764_f64,
                        iterations: Some(5),
                        check_iterations: true,
                        function_calls: 7,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Secant,
                    },
                ),
                (
                    "kepler halley x0=1.0 fprime=fn fprime2=fn xtol=1e-12 rtol=1e-12",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 0.9 * x.sin() - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(1.0_f64),
                            xtol: Some(1.0e-12_f64),
                            rtol: Some(1.0e-12_f64),
                            fprime: Some(Derivative::function(|x| Ok(1. - 0.9 * x.cos()))),
                            fprime2: Some(Derivative::function(|x| Ok(0.9 * x.sin()))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.8620866868745323_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 13,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Halley,
                    },
                ),
                (
                    "steep newton x0=0.29",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok((1000. * x - 300.).atan()) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(0.29_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1403199.9271292472_f64,
                        iterations: Some(4),
                        check_iterations: true,
                        function_calls: 8,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "steep secant x0=0.29",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok((1000. * x - 300.).atan()) },
                        RootScalarOptions {
                            method: Some(Method::Secant),
                            x0: Some(0.29_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: -9935490040481.354_f64,
                        iterations: Some(9),
                        check_iterations: true,
                        function_calls: 10,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Secant,
                    },
                ),
                (
                    "cubic newton x0=-0.0",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(-0.0_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: -0.0_f64,
                        iterations: Some(1),
                        check_iterations: true,
                        function_calls: 2,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "cubic newton x0=-5.0",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(-5.0_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(20),
                        check_iterations: true,
                        function_calls: 40,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "cubic secant x0=-0.0",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Secant),
                            x0: Some(-0.0_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.9999115140982547e-16_f64,
                        iterations: Some(3),
                        check_iterations: true,
                        function_calls: 4,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Secant,
                    },
                ),
                (
                    "sqplus1 newton x0=0.0 fprime=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x + 1.) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(0.0_f64),
                            fprime: Some(Derivative::function(|x| Ok(2. * x))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.0_f64,
                        iterations: Some(1),
                        check_iterations: true,
                        function_calls: 2,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "sqplus1 halley x0=0.0 fprime=fn fprime2=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x + 1.) },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(0.0_f64),
                            fprime: Some(Derivative::function(|x| Ok(2. * x))),
                            fprime2: Some(Derivative::function(|_x| Ok(2.))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 0.0_f64,
                        iterations: Some(1),
                        check_iterations: true,
                        function_calls: 2,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "sqplus1 newton x0=0.5 fprime=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x + 1.) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(0.5_f64),
                            fprime: Some(Derivative::function(|x| Ok(2. * x))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: -0.8707527744354187_f64,
                        iterations: Some(50),
                        check_iterations: true,
                        function_calls: 100,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "constant secant x0=1.0",
                    root_scalar(
                        |_x: f64| -> anyhow::Result<f64> { Ok(1.) },
                        RootScalarOptions {
                            method: Some(Method::Secant),
                            x0: Some(1.0_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0001_f64,
                        iterations: Some(1),
                        check_iterations: true,
                        function_calls: 2,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Secant,
                    },
                ),
                (
                    "constant newton x0=1.0",
                    root_scalar(
                        |_x: f64| -> anyhow::Result<f64> { Ok(1.) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(1.0_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(1),
                        check_iterations: true,
                        function_calls: 2,
                        converged: false,
                        flag: r#"convergence error"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "shifted halley x0=1.0 fprime=fn fprime2=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(1.0_f64),
                            fprime: Some(Derivative::function(|_x| Ok(1.))),
                            fprime2: Some(Derivative::function(|_x| Ok(0.))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(0),
                        check_iterations: true,
                        function_calls: 1,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "shifted newton x0=1.0 fprime=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(1.0_f64),
                            fprime: Some(Derivative::function(|_x| Ok(1.))),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(0),
                        check_iterations: true,
                        function_calls: 1,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Newton,
                    },
                ),
                (
                    "shifted secant x0=1.0",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Secant),
                            x0: Some(1.0_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::Ok {
                        root: 1.0_f64,
                        iterations: Some(2),
                        check_iterations: true,
                        function_calls: 3,
                        converged: true,
                        flag: r#"converged"#,
                        method: Method::Secant,
                    },
                ),
                (
                    "cubic brentq bracket=(2, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Brentq),
                            bracket: Some((2.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"f(a) and f(b) must have different signs"#),
                ),
                (
                    "cubic toms748 bracket=(2, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((2.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(
                        r#"f(a) and f(b) must have different signs, but f(2.000000e+00)=7.000000e+00, f(3.000000e+00)=2.600000e+01 "#,
                    ),
                ),
                (
                    "cubic toms748 bracket=(3, 0)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((3.0_f64, 0.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"a and b are not an interval [3.0, 0.0]"#),
                ),
                (
                    "cubic brentq bracket=(0, 3) xtol=0.0",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Brentq),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            xtol: Some(0.0_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"xtol too small (0 <= 0)"#),
                ),
                (
                    "cubic bisect bracket=(0, 3) rtol=1e-17",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Bisect),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            rtol: Some(1.0e-17_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"rtol too small (1e-17 < 8.88178e-16)"#),
                ),
                (
                    "cubic toms748 bracket=(0, 3) rtol=1e-17",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            rtol: Some(1.0e-17_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"rtol too small (1e-17 < 2.22045e-16)"#),
                ),
                (
                    "cubic toms748 bracket=(0, 3) xtol=-0.001",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            xtol: Some(-0.001_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"xtol too small (-0.001 <= 0)"#),
                ),
                (
                    "cubic toms748 bracket=(0, 3) maxiter=0",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            maxiter: Some(0),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"maxiter must be greater than 0"#),
                ),
                (
                    "cubic toms748 bracket=(0, 3) k=0",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            k: Some(0),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"k too small (0 < 1)"#),
                ),
                (
                    "cubic toms748 bracket=(0, inf)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Toms748),
                            bracket: Some((0.0_f64, f64::INFINITY)),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"b is not finite inf"#),
                ),
                (
                    "cubic newton x0=0.2 xtol=0.0",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(0.2_f64),
                            xtol: Some(0.0_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"tol too small (0 <= 0)"#),
                ),
                (
                    "cubic newton x0=0.2 maxiter=0",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            x0: Some(0.2_f64),
                            maxiter: Some(0),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"maxiter must be greater than 0"#),
                ),
                (
                    "cubic secant x0=0.2 x1=0.2",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Secant),
                            x0: Some(0.2_f64),
                            x1: Some(0.2_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"x1 and x0 must be different"#),
                ),
                (
                    "cubic halley x0=0.2 fprime=fn",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(0.2_f64),
                            fprime: Some(Derivative::function(|x| Ok(3. * x * x))),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"fprime2 must be specified for halley"#),
                ),
                (
                    "cubic halley x0=0.2",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Halley),
                            x0: Some(0.2_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"fprime must be specified for halley"#),
                ),
                (
                    "cubic brentq x0=0.2",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Brentq),
                            x0: Some(0.2_f64),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"Bracket needed for brentq"#),
                ),
                (
                    "cubic newton bracket=(0, 3)",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Newton),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"x0 must not be None for newton"#),
                ),
                (
                    "cubic secant",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Secant),
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(r#"x0 must not be None for secant"#),
                ),
                (
                    "cubic",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            ..Default::default()
                        },
                    ),
                    Expected::ValueErr(
                        r#"Unable to select a solver as neither bracket nor starting point provided."#,
                    ),
                ),
                (
                    "cubic brentq bracket=(0, 3) k=2",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Brentq),
                            bracket: Some((0.0_f64, 3.0_f64)),
                            k: Some(2),
                            ..Default::default()
                        },
                    ),
                    Expected::TypeErr(r#"brentq() got an unexpected keyword argument 'k'"#),
                ),
                (
                    "cubic secant x0=0.2 k=2",
                    root_scalar(
                        |x: f64| -> anyhow::Result<f64> { Ok(x * x * x - 1.) },
                        RootScalarOptions {
                            method: Some(Method::Secant),
                            x0: Some(0.2_f64),
                            k: Some(2),
                            ..Default::default()
                        },
                    ),
                    Expected::TypeErr(r#"newton() got an unexpected keyword argument 'k'"#),
                ),
            ];
            let total = cases.len();
            let failures: Vec<String> = cases
                .into_iter()
                .filter_map(|(description, actual, expected)| check(description, actual, expected))
                .collect();
            assert!(
                failures.is_empty(),
                "{} of {total} cases failed:\n{}",
                failures.len(),
                failures.join("\n")
            );
        }
    }
}

/// direct port of scipy.integrate.solve_ivp but only including areas used in this library
pub mod solve_ivp {
    use crate::core::solvers::solve_ivp::base_solver::{DenseOutput, Status};
    use anyhow::Error;
    use argmin::core::{CostFunction, Executor};
    use argmin::solver::brent::BrentRoot;
    use ndarray::{s, Array, Array1, Axis, Dimension, Zip};
    use std::cmp::Ordering;
    use std::sync::Arc;
    use thiserror::Error;

    pub type SharedIvpSolveFunction = Arc<dyn Fn(f64, &Array1<f64>) -> Array1<f64> + Send + Sync>;
    pub type SharedIvpEventFunction = Arc<dyn Fn(f64, &Array1<f64>) -> f64 + Send + Sync>;

    pub fn solve_ivp(
        func: &SharedIvpSolveFunction,
        t_span: (f64, f64),
        y0: &Array1<f64>,
        events: Option<&[TerminatingEvent]>,
        rtol: Option<f64>,
        atol: Option<f64>,
    ) -> Result<OdeResult, SolverError> {
        let (t0, tf) = t_span;

        let mut solver = rk45::Rk45Solver::new(func, t0, y0, tf, rtol, atol);

        let mut ts = vec![t0];
        let mut ys = vec![y0.to_owned()];

        let mut event_data: Option<EventData> = events.map(|events| {
            let (max_events, event_dir) = prepare_events(events);

            EventData {
                t_events: vec![vec![]; events.len()],
                y_events: vec![vec![]; events.len()],
                events,
                max_events,
                event_dir,
                event_count: Array1::zeros(events.len()),
                g: events.iter().map(|event| (event.func)(t0, y0)).collect(),
            }
        });

        let mut _status: Option<i8> = None;

        while _status.is_none() {
            solver
                .step()
                .map_err(|e| SolverError::MainSolver(e.into()))?;

            if solver.status == Status::Finished {
                _status = Some(0);
            } else if solver.status == Status::Failed {
                _status = Some(-1);
                break;
            }

            let t_old = solver
                .t_old
                .expect("t_old is expected to have been initialised after the first step");
            let mut t = solver.t;
            let mut y = solver.y.clone();

            let mut sol: Option<Arc<dyn DenseOutput>> = None;

            if let Some(EventData {
                events,
                ref event_dir,
                ref mut g,
                ref mut event_count,
                ref max_events,
                ref mut t_events,
                ref mut y_events,
                ..
            }) = event_data
            {
                let g_new: Array1<f64> = events.iter().map(|event| (event.func)(t, &y)).collect();
                let active_events = find_active_events(g, &g_new, event_dir);
                if !active_events.is_empty() {
                    if sol.is_none() {
                        sol = solver.dense_output().into();
                    }

                    let sol = &sol.unwrap();

                    for active_event in active_events.iter() {
                        event_count[*active_event] += 1;
                    }
                    let (root_indices, roots, terminate) = handle_events(
                        sol,
                        events,
                        active_events,
                        event_count,
                        max_events,
                        t_old,
                        t,
                    )?;

                    for (&e, &te) in root_indices.iter().zip(roots.iter()) {
                        t_events[e].push(te);
                        y_events[e].push(sol.call_impl(te));
                    }

                    if terminate {
                        _status = 1.into();
                        t = roots.last().copied().unwrap();
                        y = sol.call_impl(t);
                    }
                }
                *g = g_new;
            }

            ts.push(t);
            ys.push(y);
        }

        let (t_events, _y_events): (Option<TEvents>, Option<YEvents>) = if let Some(EventData {
            ref t_events,
            ref y_events,
            ..
        }) = event_data
        {
            (
                Some(
                    t_events
                        .iter()
                        .map(|te| te.iter().copied().collect())
                        .collect(),
                ),
                Some(
                    y_events
                        .iter()
                        .map(|ye| ye.iter().map(|ye| ye.to_owned()).collect())
                        .collect(),
                ),
            )
        } else {
            (None, None)
        };

        let ts = Array1::from_vec(ts);
        let (ys, _) = Array1::from_vec(ys)
            .insert_axis(Axis(1))
            .into_raw_vec_and_offset();

        Ok(OdeResult {
            t: ts,
            y: ys,
            t_events,
            success: true,
        })
    }

    struct EventData<'a> {
        t_events: Vec<Vec<f64>>,
        y_events: YEvents,
        events: &'a [TerminatingEvent],
        max_events: Array1<f64>,
        event_dir: Array1<TerminateDirection>,
        event_count: Array1<usize>,
        g: Array1<f64>,
    }

    type TEvents = Vec<Array1<f64>>;
    type YEvents = Vec<Vec<Array1<f64>>>;

    fn prepare_events(events: &[TerminatingEvent]) -> (Array1<f64>, Array1<TerminateDirection>) {
        let mut max_events = Array1::zeros(events.len());
        let mut event_dir = Array1::from_elem(events.len(), TerminateDirection::Both);

        for (i, event) in events.iter().enumerate() {
            max_events[i] = if event.is_terminal() {
                1.0
            } else {
                f64::INFINITY
            };
            event_dir[i] = event.direction;
        }

        (max_events, event_dir)
    }

    fn handle_events(
        sol: &Arc<dyn DenseOutput>,
        events: &[TerminatingEvent],
        mut active_events: Array1<usize>,
        event_count: &Array1<usize>,
        max_events: &Array1<f64>,
        t_old: f64,
        t: f64,
    ) -> Result<(Array1<usize>, Array1<f64>, bool), SolverError> {
        let mut roots: Array1<f64> = events
            .iter()
            .enumerate()
            .filter(|(i, _event)| {
                active_events
                    .as_slice()
                    .is_some_and(|slice| slice.contains(i))
            })
            .map(|(_i, event)| solve_event_equation(&event.func, sol, t_old, t))
            .collect::<Result<_, _>>()?;

        let terminate = active_events
            .iter()
            .any(|active_event| event_count[*active_event] >= max_events[*active_event] as usize);

        if terminate {
            let order = argsort(
                &roots,
                t.partial_cmp(&t_old).ok_or_else(|| {
                    SolverError::Event(NotANumberEncounteredDuringEventHandling.into())
                })? != Ordering::Greater,
            );
            active_events = order.iter().map(|&i| active_events[i]).collect();
            roots = order.iter().map(|&i| roots[i]).collect();

            let count_view = event_count.select(Axis(0), active_events.as_slice().unwrap());
            let max_view = max_events.select(Axis(0), active_events.as_slice().unwrap());
            let t = Zip::indexed(&count_view)
                .and(&max_view)
                .fold(None, |acc, idx, &count, &max| {
                    if acc.is_some() {
                        acc
                    } else if count as f64 >= max {
                        Some(idx)
                    } else {
                        None
                    }
                })
                .unwrap();

            active_events = active_events.slice(s![..t + 1]).to_owned();
            roots = roots.slice(s![..t + 1]).to_owned();
        }

        Ok((active_events, roots, terminate))
    }

    #[derive(Debug, Error)]
    #[error("A NaN value was encountered during event handling for either t or t_old")]
    struct NotANumberEncounteredDuringEventHandling;

    /// Perform equivalent of `np.argsort`. NB. will panic if passed an empty array.
    fn argsort<T>(arr: &Array1<T>, reverse: bool) -> Array1<usize>
    where
        T: PartialOrd,
    {
        let mut indices: Vec<usize> = (0..arr.len()).collect();

        indices.sort_by(|&a, &b| arr[a].partial_cmp(&arr[b]).unwrap());

        if reverse {
            indices.reverse();
        }

        Array1::from_vec(indices)
    }

    fn find_active_events(
        g: &Array1<f64>,
        g_new: &Array1<f64>,
        direction: &Array1<TerminateDirection>,
    ) -> Array1<usize> {
        let mut indices = Vec::new();

        Zip::indexed(g)
            .and(g_new)
            .and(direction)
            .for_each(|idx, &g_val, &g_new_val, &dir_val| {
                let up = g_val <= 0.0 && g_new_val >= 0.0;
                let down = g_val >= 0.0 && g_new_val <= 0.0;
                let either = up || down;

                let mask = (up && dir_val > TerminateDirection::Both)
                    || (down && dir_val < TerminateDirection::Both)
                    || (either && dir_val == TerminateDirection::Both);

                if mask {
                    indices.push(idx);
                }
            });

        Array1::from(indices)
    }

    fn solve_event_equation(
        event: &SharedIvpEventFunction,
        sol: &Arc<dyn DenseOutput>,
        t_old: f64,
        t: f64,
    ) -> Result<f64, SolverError> {
        let problem = EventRootFunction {
            event: Arc::clone(event),
            sol: Arc::clone(sol),
        };

        let solver = BrentRoot::new(t_old, t, 1e-8);

        let executor = Executor::new(problem, solver);

        let result = executor.run();

        match result.map(|r| r.state().best_param) {
            Ok(Some(t)) => Ok(t),
            Ok(None) => Err(SolverError::Event(Box::new(NoRootFoundError))),
            Err(err) => Err(SolverError::Event(err.into())),
        }
    }

    struct EventRootFunction {
        event: SharedIvpEventFunction,
        sol: Arc<dyn DenseOutput>,
    }

    #[derive(Debug, Error)]
    #[error("No root found for IVP solving event equation")]
    struct NoRootFoundError;

    impl CostFunction for EventRootFunction {
        type Param = f64;
        type Output = f64;

        fn cost(&self, t: &Self::Param) -> Result<Self::Output, Error> {
            Ok((self.event)(*t, &(self.sol).call_impl(*t)))
        }
    }

    #[derive(Debug, Error)]
    pub enum SolverError {
        #[error("Error when solving events: {0}")]
        Event(Box<dyn std::error::Error + Send + Sync>),
        #[error("Error in main IVP solver: {0}")]
        MainSolver(Box<dyn std::error::Error + Send + Sync>),
        #[error("Error happened during event termination")]
        Termination,
    }

    pub struct TerminatingEvent {
        func: SharedIvpEventFunction,
        direction: TerminateDirection,
    }

    impl TerminatingEvent {
        pub fn new(func: SharedIvpEventFunction, direction: Option<TerminateDirection>) -> Self {
            Self {
                func,
                direction: direction.unwrap_or_default(),
            }
        }

        pub fn is_terminal(&self) -> bool {
            true
        }
    }

    #[derive(Clone, Copy, Debug, Default, Eq, PartialEq, PartialOrd)]
    pub enum TerminateDirection {
        #[default]
        Both = 0,
        Positive = 1,
        Negative = -1,
    }

    pub struct OdeResult {
        pub y: Vec<Array1<f64>>,
        pub t_events: Option<Vec<Array1<f64>>>,
        pub t: Array1<f64>,
        pub success: bool,
    }

    /// module to do with the RK45 method for scipy
    mod rk45 {
        #![allow(clippy::unreadable_literal)]

        use crate::core::solvers::solve_ivp::base_solver::{DefaultDenseOutput, DenseOutput};
        use crate::core::solvers::solve_ivp::{select_initial_step, SharedIvpSolveFunction};
        use ndarray::{array, concatenate, s, Array1, Array2, Axis, Zip};
        use std::sync::{Arc, LazyLock};
        use thiserror::Error;

        #[allow(clippy::many_single_char_names)]
        fn rk_step(
            fun: &SharedIvpSolveFunction,
            t: f64,
            y: &Array1<f64>,
            f: &Array1<f64>,
            h: f64,
            a: &Array2<f64>,
            b: &Array1<f64>,
            c: &Array1<f64>,
            k: &mut Array2<f64>,
        ) -> (Array1<f64>, Array1<f64>) {
            // 1. K[0] = f
            // Assigning a 1D vector to a specific row of our 2D matrix
            k.row_mut(0).assign(f);

            let n_stages = b.len();

            // 2. Loop over the intermediate stages
            for s in 1..n_stages {
                let c = c[s];

                // Extract a slice of the coefficients up to index s: a[:s]
                let a_slice = a.slice(s![s, ..s]);

                // Extract a slice of previous K stages up to row s: K[:s]
                // In Python: np.dot(K[:s].T, a[:s])
                // In Rust: We match the matrix multiplication mechanics:
                // Transposing rows into columns with .t(), then computing a matrix-vector dot product
                let k_slice = k.slice(s![..s, ..]);
                let dy = k_slice.t().dot(&a_slice) * h;

                // K[s] = fun(t + c * h, y + dy)
                let next_stage_y = y + &dy;
                let f_stage = fun(t + c * h, &next_stage_y);
                k.row_mut(s).assign(&f_stage);
            }

            // 3. Compute final y prediction
            // Python: y_new = y + h * np.dot(K[:-1].T, B)
            let k_all_stages = k.slice(s![..n_stages, ..]);
            let y_new = y + (k_all_stages.t().dot(b) * h);

            // 4. Compute final f derivative evaluation
            // f_new = fun(t + h, y_new)
            let f_new = fun(t + h, &y_new);

            // 5. K[-1] = f_new
            // Assigning the last row of the K matrix
            let last_row_idx = k.nrows() - 1;
            k.row_mut(last_row_idx).assign(&f_new);

            (y_new, f_new)
        }

        #[derive(Clone)]
        pub(super) struct Rk45Solver {
            n: usize,
            pub(crate) status: super::base_solver::Status,
            t_bound: f64,
            direction: i8,
            pub(crate) t: f64,
            pub(crate) y: Array1<f64>,
            pub(crate) t_old: Option<f64>,
            step_size: Option<f64>,
            nfev: usize,
            fun: SharedIvpSolveFunction,
            y_old: Option<Array1<f64>>,
            max_step: f64,
            rtol: f64,
            atol: f64,
            f: Array1<f64>,
            h_abs: f64,
            k: Array2<f64>,
            h_previous: Option<f64>,
        }

        const ORDER: usize = 5;
        const ERROR_ESTIMATOR_ORDER: usize = 4;
        const N_STAGES: usize = 6;
        static C: LazyLock<Array1<f64>> =
            LazyLock::new(|| array![0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0]);
        static A: LazyLock<Array2<f64>> = LazyLock::new(|| {
            array!(
                [0., 0., 0., 0., 0.],
                [1. / 5., 0., 0., 0., 0.],
                [3. / 40., 9. / 40., 0., 0., 0.],
                [44. / 45., -56. / 15., 32. / 9., 0., 0.],
                [
                    19372. / 6561.,
                    -25360. / 2187.,
                    64448. / 6561.,
                    -212. / 729.,
                    0.
                ],
                [
                    9017. / 3168.,
                    -355. / 33.,
                    46732. / 5247.,
                    49. / 176.,
                    -5103. / 18656.
                ]
            )
        });
        static B: LazyLock<Array1<f64>> = LazyLock::new(|| {
            array![
                35. / 384.,
                0.,
                500. / 1113.,
                125. / 192.,
                -2187. / 6784.,
                11. / 84.
            ]
        });
        static E: LazyLock<Array1<f64>> = LazyLock::new(|| {
            array![
                -71. / 57600.,
                0.,
                71. / 16695.,
                -71. / 1920.,
                17253. / 339200.,
                -22. / 525.,
                1. / 40.
            ]
        });
        static P: LazyLock<Array2<f64>> = LazyLock::new(|| {
            array![
                [
                    1.,
                    -8048581381. / 2820520608.,
                    8663915743. / 2820520608.,
                    -12715105075. / 11282082432.
                ],
                [0., 0., 0., 0.],
                [
                    0.,
                    131558114200. / 32700410799.,
                    -68118460800. / 10900136933.,
                    87487479700. / 32700410799.
                ],
                [
                    0.,
                    -1754552775. / 470086768.,
                    14199869525. / 1410260304.,
                    -10690763975. / 1880347072.
                ],
                [
                    0.,
                    127303824393. / 49829197408.,
                    -318862633887. / 49829197408.,
                    701980252875. / 199316789632.
                ],
                [
                    0.,
                    -282668133. / 205662961.,
                    2019193451. / 616988883.,
                    -1453857185. / 822651844.
                ],
                [
                    0.,
                    40617522. / 29380423.,
                    -110615467. / 29380423.,
                    69997945. / 29380423.
                ]
            ]
        });
        const ERROR_EXPONENT: f64 = -1. / (ERROR_ESTIMATOR_ORDER as f64 + 1.);

        // Multiply steps computed from asymptotic behaviour of errors by this.
        const SAFETY: f64 = 0.9;

        // Minimum allowed decrease in a step size.
        const MIN_FACTOR: f64 = 0.2;

        // Maximum allowed increase in a step size.
        const MAX_FACTOR: f64 = 10.;

        impl Rk45Solver {
            pub(super) fn new(
                fun: &SharedIvpSolveFunction,
                t0: f64,
                y0: &Array1<f64>,
                t_bound: f64,
                rtol: Option<f64>,
                atol: Option<f64>,
            ) -> Self {
                let f = fun(t0, y0);
                let t = t0;
                let y = y0.to_owned();
                let max_step = f64::INFINITY;
                let direction = if t_bound != t0 { sign(t_bound - t0) } else { 1 };
                let rtol = rtol.unwrap_or(1e-3);
                let atol = atol.unwrap_or(1e-6);
                let h_abs = select_initial_step(
                    fun,
                    t0,
                    y0,
                    t_bound,
                    max_step,
                    &f,
                    direction as f64,
                    ERROR_ESTIMATOR_ORDER,
                    rtol,
                    atol,
                );
                let n = y.len();

                Self {
                    n,
                    status: Default::default(),
                    t_bound,
                    direction,
                    t,
                    y,
                    t_old: None,
                    step_size: None,
                    nfev: 0,
                    fun: fun.clone(),
                    y_old: None,
                    max_step,
                    rtol,
                    atol,
                    f,
                    h_abs,
                    k: Array2::zeros((N_STAGES + 1, n)),
                    h_previous: None,
                }
            }

            fn step_size(&self) -> Option<f64> {
                self.t_old.map(|t_old| (self.t - t_old).abs())
            }

            pub(crate) fn step(&mut self) -> Result<(), Rk45SolverError> {
                if self.status != super::base_solver::Status::Running {
                    panic!("Attempt to step on a failed or finished solver.");
                }

                if self.n == 0 || self.t == self.t_bound {
                    // Handle corner cases of empty solver or no integration.
                    self.t_old = self.t.into();
                    self.t = self.t_bound;
                    self.status = super::base_solver::Status::Finished;

                    Ok(())
                } else {
                    let t = self.t;
                    let result = self.step_impl();
                    if result.is_err() {
                        self.status = super::base_solver::Status::Failed;
                    } else {
                        self.t_old = t.into();
                        if self.direction as f64 * (self.t - self.t_bound) >= 0. {
                            self.status = super::base_solver::Status::Finished;
                        }
                    }

                    result
                }
            }

            fn step_impl(&mut self) -> Result<(), Rk45SolverError> {
                let t = self.t;
                let y = &self.y;

                let max_step = self.max_step;
                let rtol = self.rtol;
                let atol = self.atol;

                let min_step =
                    10. * f64::abs(nextafter(t, self.direction as f64 * f64::INFINITY) - t);

                let mut h_abs = if self.h_abs > max_step {
                    max_step
                } else if self.h_abs < min_step {
                    min_step
                } else {
                    self.h_abs
                };

                let mut step_accepted = false;
                let mut step_rejected = false;

                let mut h = f64::default();
                let mut f_new = Array1::zeros(y.dim());
                let mut y_new = Array1::zeros(y.dim());

                let mut t_new = f64::default();

                while !step_accepted {
                    if h_abs < min_step {
                        return Err(Rk45SolverError::TooSmallStep);
                    }

                    h = h_abs * self.direction as f64;
                    t_new = t + h;

                    if self.direction as f64 * (t_new - self.t_bound) > 0. {
                        t_new = self.t_bound;
                    }

                    h = t_new - t;
                    h_abs = f64::abs(h);

                    let (y_new_step, f_new_step) =
                        rk_step(&self.fun, t, y, &self.f, h, &A, &B, &C, &mut self.k);
                    (y_new, f_new) = (y_new_step, f_new_step);

                    let mut scale = Array1::<f64>::zeros(y_new.dim());
                    Zip::from(&mut scale).and(y).and(&y_new).for_each(
                        |scale_val, y_val, y_new_val| {
                            *scale_val = atol + y_val.abs().max(y_new_val.abs()) * rtol;
                        },
                    );
                    let error_norm = self.estimate_error_norm(&self.k, h, &scale);

                    if error_norm < 1.0 {
                        let mut factor = if error_norm == 0.0 {
                            MAX_FACTOR
                        } else {
                            MAX_FACTOR.min(SAFETY * error_norm.powf(ERROR_EXPONENT))
                        };

                        if step_rejected {
                            factor = factor.min(1.);
                        }

                        h_abs *= factor;

                        step_accepted = true;
                    } else {
                        h_abs *= MIN_FACTOR.max(SAFETY * error_norm.powf(ERROR_EXPONENT));
                        step_rejected = true;
                    }
                }

                self.h_previous = h.into();
                self.y_old = y.to_owned().into();

                self.t = t_new;
                self.y = y_new;

                self.h_abs = h_abs;
                self.f = f_new;

                Ok(())
            }

            pub(crate) fn dense_output(&self) -> Arc<dyn DenseOutput> {
                if self.n == 0 || self.t_old.is_some_and(|t_old| self.t == t_old) {
                    Arc::new(
                        crate::core::solvers::solve_ivp::base_solver::ConstantDenseOutput::new(
                            self.t_old
                                .expect("Expected to be past the first step with t_old set"),
                            self.t,
                            self.y.clone(),
                        ),
                    )
                } else {
                    Arc::new(self.dense_output_impl())
                }
            }

            fn estimate_error(&self, k: &Array2<f64>, h: f64) -> Array1<f64> {
                let e: &Array1<f64> = &E;

                k.t().dot(e) * h
            }

            fn estimate_error_norm(&self, k: &Array2<f64>, h: f64, scale: &Array1<f64>) -> f64 {
                super::norm(&(self.estimate_error(k, h) / scale))
            }

            fn dense_output_impl(&self) -> impl super::base_solver::DenseOutput {
                let p: &Array2<f64> = &P;
                let q = self.k.t().dot(p);

                RkDenseOutput::new(
                    self.t_old.expect(
                        "t_old expects to have been initialised before emitting a dense output",
                    ),
                    self.t,
                    self.y_old.clone().expect(
                        "y_old expects to have been initialised before emitting a dense output",
                    ),
                    q,
                )
            }
        }

        #[derive(Debug, Error)]
        pub enum Rk45SolverError {
            #[error("Required step size is less than spacing between numbers.")]
            TooSmallStep,
        }

        // equivalent of np.nextafter
        fn nextafter(x: f64, y: f64) -> f64 {
            if y > x {
                x.next_up()
            } else if y < x {
                x.next_down()
            } else {
                x
            }
        }

        struct RkDenseOutput {
            _default: DefaultDenseOutput,
            h: f64,
            q: Array2<f64>,
            order: usize,
            y_old: Array1<f64>,
        }

        impl RkDenseOutput {
            fn new(t_old: f64, t: f64, y_old: Array1<f64>, q: Array2<f64>) -> Self {
                let order = q.shape()[1] - 1;

                Self {
                    _default: super::base_solver::dense_output_init(t_old, t),
                    h: t - t_old,
                    q,
                    order,
                    y_old,
                }
            }

            fn t_old(&self) -> f64 {
                self._default.t_old
            }
        }

        impl DenseOutput for RkDenseOutput {
            fn call_impl(&self, t: f64) -> Array1<f64> {
                let t = array![t];
                let x = (&t - self.t_old()) / self.h;

                let p = if t.ndim() == 0 {
                    let p = x.broadcast((self.order + 1, x.len())).unwrap();
                    p.flatten().cumprod(Axis(0))
                } else {
                    let views = vec![x.view(); self.order + 1];
                    let p = concatenate(Axis(0), &views).unwrap();
                    p.flatten().cumprod(Axis(0))
                };

                let mut y = self.h * self.q.dot(&p);
                if y.ndim() == 2 {
                    let new_line = self.y_old.clone().insert_axis(Axis(1));
                    y += &new_line;
                } else {
                    y += &self.y_old;
                }

                y
            }
        }

        fn sign(value: f64) -> i8 {
            if value < 0. {
                -1
            } else if value == 0. {
                0
            } else {
                1
            }
        }
    }

    fn select_initial_step(
        fun: &SharedIvpSolveFunction,
        t0: f64,
        y0: &Array1<f64>,
        t_bound: f64,
        max_step: f64,
        f0: &Array1<f64>,
        direction: f64,
        order: usize,
        rtol: f64,
        atol: f64,
    ) -> f64 {
        if y0.is_empty() {
            return f64::INFINITY;
        }

        let interval_length = (t_bound - t0).abs();
        if interval_length == 0.0 {
            return 0.0;
        }

        let scale: Array1<f64> = y0.iter().map(|&y| (y.abs() * rtol) + atol).collect();

        let d0 = norm(&(y0 / &scale));
        let d1 = norm(&(f0 / &scale));

        let h0 = if d0 < 1e-5 || d1 < 1e-5 {
            1e-6
        } else {
            0.01 * d0 / d1
        };

        // Check t0+h0*direction doesn't take us beyond t_bound
        let h0 = h0.min(interval_length);
        let y1 = y0 + h0 * direction * f0;
        let f1 = fun(t0 + h0 * direction, &y1);
        let d2 = norm(&((f1 - f0) / &scale)) / h0;

        let h1 = if d1 <= 1e-15 && d2 <= 1e-15 {
            1e-6_f64.max(h0 * 1e-3)
        } else {
            (0.01 / d1.max(d2)).powf(1.0 / (order + 1) as f64)
        };

        [100. * h0, h1, interval_length, max_step]
            .into_iter()
            .min_by(|a, b| a.total_cmp(b))
            .unwrap()
    }

    // functional equivalent of the common norm function inside scipy.integrate ivp library
    fn norm<D: Dimension>(arr: &Array<f64, D>) -> f64 {
        numpy_linalg_norm_equivalent(arr) / (arr.len() as f64).powf(0.5)
    }

    fn numpy_linalg_norm_equivalent<D: Dimension>(arr: &Array<f64, D>) -> f64 {
        arr.iter().map(|&x| x * x).sum::<f64>().sqrt()
    }

    mod base_solver {
        use ndarray::Array1;

        #[derive(Clone, Copy, Debug, Default, PartialEq)]
        pub enum Status {
            #[default]
            Running,
            Finished,
            Failed,
        }

        pub trait DenseOutput {
            fn call_impl(&self, t: f64) -> Array1<f64>;
        }

        pub struct DefaultDenseOutput {
            pub(super) t_old: f64,
            t: f64,
            t_min: f64,
            t_max: f64,
        }

        pub struct ConstantDenseOutput {
            _default: DefaultDenseOutput,
            value: Array1<f64>,
        }

        impl ConstantDenseOutput {
            pub fn new(t_old: f64, t: f64, value: Array1<f64>) -> Self {
                Self {
                    _default: dense_output_init(t_old, t),
                    value,
                }
            }
        }

        pub(super) fn dense_output_init(t_old: f64, t: f64) -> DefaultDenseOutput {
            DefaultDenseOutput {
                t_old,
                t,
                t_min: t.min(t_old),
                t_max: t.max(t_old),
            }
        }

        impl DenseOutput for ConstantDenseOutput {
            fn call_impl(&self, _t: f64) -> Array1<f64> {
                self.value.clone()
            }
        }
    }
}

pub type SharedInterpolationFunction = Arc<dyn Fn(&[f64]) -> Vec<f64> + Send + Sync>;

pub fn interp1d(
    x: Vec<f64>,
    y: Vec<f64>,
    fill_value: Interp1dFillValue,
) -> SharedInterpolationFunction {
    let interp_mode = match fill_value {
        Interp1dFillValue::Extrapolate => InterpMode::Extrapolate,
        Interp1dFillValue::FillValues(values) => InterpMode::Constant(values.0),
    };

    let func = move |x_new: &[f64]| interp_slice(&x, &y, x_new, &interp_mode).to_vec();

    Arc::new(func)
}

pub enum Interp1dFillValue {
    Extrapolate,
    FillValues((f64, f64)),
}
