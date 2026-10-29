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

/// direct port of scipy.integrate.solve_ivp but only including areas used in this library
///
/// The state is a fixed-size array so that the integration loop never touches the heap.
pub mod solve_ivp {
    use argmin::core::{CostFunction, Executor};
    use argmin::solver::brent::BrentRoot;
    use thiserror::Error;

    /// Integrate `fun` over `t_span` from `y0` using the RK45 method.
    ///
    /// All events are terminal, so integration stops at the first event that occurs, and only
    /// the final state is returned (which is all that callers in this library use).
    pub fn solve_ivp<const N: usize>(
        mut fun: impl FnMut(f64, &[f64; N]) -> [f64; N],
        t_span: (f64, f64),
        y0: [f64; N],
        events: &[TerminatingEvent<'_, N>],
        rtol: Option<f64>,
        atol: Option<f64>,
    ) -> Result<OdeResult<N>, SolverError> {
        let (t0, tf) = t_span;

        let mut solver = rk45::Rk45Solver::new(&mut fun, t0, y0, tf, rtol, atol);

        // allocated once per solve and reused on every step
        let mut g: Vec<f64> = events.iter().map(|event| (event.func)(t0, &y0)).collect();
        let mut g_new = g.clone();

        loop {
            let finished = solver.step(&mut fun)?;

            if !events.is_empty() {
                let (t_old, t) = (solver.t_old, solver.t);

                for (g_new_val, event) in g_new.iter_mut().zip(events) {
                    *g_new_val = (event.func)(t, &solver.y);
                }

                // equivalent of np.argsort over the roots of the active events (reversed when
                // integrating backwards) and then taking the first, as all events are terminal
                let reverse = t <= t_old;
                let sol = solver.dense_output();
                let mut first_event: Option<(usize, f64)> = None;

                for (i, event) in events.iter().enumerate() {
                    if !is_active_event(g[i], g_new[i], event.direction) {
                        continue;
                    }

                    let root = solve_event_equation(event.func, &sol, t_old, t)?;

                    let is_first = match first_event {
                        None => true,
                        Some((_, first_root)) if reverse => root >= first_root,
                        Some((_, first_root)) => root < first_root,
                    };
                    if is_first {
                        first_event = Some((i, root));
                    }
                }

                if let Some((i, t_event)) = first_event {
                    return Ok(OdeResult {
                        t: t_event,
                        y: sol.call(t_event),
                        event: Some((i, t_event)),
                    });
                }

                std::mem::swap(&mut g, &mut g_new);
            }

            if finished {
                return Ok(OdeResult {
                    t: solver.t,
                    y: solver.y,
                    event: None,
                });
            }
        }
    }

    fn is_active_event(g: f64, g_new: f64, direction: TerminateDirection) -> bool {
        let up = g <= 0.0 && g_new >= 0.0;
        let down = g >= 0.0 && g_new <= 0.0;

        match direction {
            TerminateDirection::Both => up || down,
            TerminateDirection::Positive => up,
            TerminateDirection::Negative => down,
        }
    }

    fn solve_event_equation<const N: usize>(
        event: &dyn Fn(f64, &[f64; N]) -> f64,
        sol: &rk45::RkDenseOutput<N>,
        t_old: f64,
        t: f64,
    ) -> Result<f64, SolverError> {
        let problem = EventRootFunction { event, sol };

        let solver = BrentRoot::new(t_old, t, 1e-8);

        let executor = Executor::new(problem, solver);

        let result = executor.run();

        match result.map(|r| r.state().best_param) {
            Ok(Some(t)) => Ok(t),
            Ok(None) => Err(SolverError::Event(Box::new(NoRootFoundError))),
            Err(err) => Err(SolverError::Event(err.into())),
        }
    }

    struct EventRootFunction<'a, const N: usize> {
        event: &'a dyn Fn(f64, &[f64; N]) -> f64,
        sol: &'a rk45::RkDenseOutput<N>,
    }

    #[derive(Debug, Error)]
    #[error("No root found for IVP solving event equation")]
    struct NoRootFoundError;

    impl<const N: usize> CostFunction for EventRootFunction<'_, N> {
        type Param = f64;
        type Output = f64;

        fn cost(&self, t: &Self::Param) -> Result<Self::Output, argmin::core::Error> {
            Ok((self.event)(*t, &self.sol.call(*t)))
        }
    }

    #[derive(Debug, Error)]
    pub enum SolverError {
        #[error("Error when solving events: {0}")]
        Event(Box<dyn std::error::Error + Send + Sync>),
        #[error("Error in main IVP solver: {0}")]
        MainSolver(#[from] rk45::Rk45SolverError),
    }

    pub struct TerminatingEvent<'a, const N: usize> {
        func: &'a dyn Fn(f64, &[f64; N]) -> f64,
        direction: TerminateDirection,
    }

    impl<'a, const N: usize> TerminatingEvent<'a, N> {
        pub fn new(
            func: &'a dyn Fn(f64, &[f64; N]) -> f64,
            direction: Option<TerminateDirection>,
        ) -> Self {
            Self {
                func,
                direction: direction.unwrap_or_default(),
            }
        }
    }

    #[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
    pub enum TerminateDirection {
        #[default]
        Both,
        Positive,
        Negative,
    }

    pub struct OdeResult<const N: usize> {
        /// final time reached (the time of the terminating event, if one occurred)
        pub t: f64,
        /// final state
        pub y: [f64; N],
        /// index of the event that terminated integration, and the time at which it occurred
        pub event: Option<(usize, f64)>,
    }

    /// module to do with the RK45 method for scipy
    mod rk45 {
        #![allow(clippy::unreadable_literal)]

        use super::select_initial_step;
        use thiserror::Error;

        const ERROR_ESTIMATOR_ORDER: usize = 4;
        const N_STAGES: usize = 6;
        const C: [f64; N_STAGES] = [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0];
        const A: [[f64; N_STAGES - 1]; N_STAGES] = [
            [0., 0., 0., 0., 0.],
            [1. / 5., 0., 0., 0., 0.],
            [3. / 40., 9. / 40., 0., 0., 0.],
            [44. / 45., -56. / 15., 32. / 9., 0., 0.],
            [
                19372. / 6561.,
                -25360. / 2187.,
                64448. / 6561.,
                -212. / 729.,
                0.,
            ],
            [
                9017. / 3168.,
                -355. / 33.,
                46732. / 5247.,
                49. / 176.,
                -5103. / 18656.,
            ],
        ];
        const B: [f64; N_STAGES] = [
            35. / 384.,
            0.,
            500. / 1113.,
            125. / 192.,
            -2187. / 6784.,
            11. / 84.,
        ];
        const E: [f64; N_STAGES + 1] = [
            -71. / 57600.,
            0.,
            71. / 16695.,
            -71. / 1920.,
            17253. / 339200.,
            -22. / 525.,
            1. / 40.,
        ];
        const P: [[f64; 4]; N_STAGES + 1] = [
            [
                1.,
                -8048581381. / 2820520608.,
                8663915743. / 2820520608.,
                -12715105075. / 11282082432.,
            ],
            [0., 0., 0., 0.],
            [
                0.,
                131558114200. / 32700410799.,
                -68118460800. / 10900136933.,
                87487479700. / 32700410799.,
            ],
            [
                0.,
                -1754552775. / 470086768.,
                14199869525. / 1410260304.,
                -10690763975. / 1880347072.,
            ],
            [
                0.,
                127303824393. / 49829197408.,
                -318862633887. / 49829197408.,
                701980252875. / 199316789632.,
            ],
            [
                0.,
                -282668133. / 205662961.,
                2019193451. / 616988883.,
                -1453857185. / 822651844.,
            ],
            [
                0.,
                40617522. / 29380423.,
                -110615467. / 29380423.,
                69997945. / 29380423.,
            ],
        ];
        const ERROR_EXPONENT: f64 = -1. / (ERROR_ESTIMATOR_ORDER as f64 + 1.);

        // Multiply steps computed from asymptotic behaviour of errors by this.
        const SAFETY: f64 = 0.9;

        // Minimum allowed decrease in a step size.
        const MIN_FACTOR: f64 = 0.2;

        // Maximum allowed increase in a step size.
        const MAX_FACTOR: f64 = 10.;

        #[allow(clippy::many_single_char_names)]
        fn rk_step<const N: usize>(
            fun: &mut impl FnMut(f64, &[f64; N]) -> [f64; N],
            t: f64,
            y: &[f64; N],
            f: &[f64; N],
            h: f64,
            k: &mut [[f64; N]; N_STAGES + 1],
        ) -> ([f64; N], [f64; N]) {
            // K[0] = f
            k[0] = *f;

            for s in 1..N_STAGES {
                // dy = np.dot(K[:s].T, a[:s]) * h
                // K[s] = fun(t + c * h, y + dy)
                let mut y_stage = *y;
                for i in 0..N {
                    let mut dy = 0.;
                    for j in 0..s {
                        dy += k[j][i] * A[s][j];
                    }
                    y_stage[i] += dy * h;
                }
                k[s] = fun(t + C[s] * h, &y_stage);
            }

            // y_new = y + h * np.dot(K[:-1].T, B)
            let mut y_new = *y;
            for i in 0..N {
                let mut dy = 0.;
                for j in 0..N_STAGES {
                    dy += k[j][i] * B[j];
                }
                y_new[i] += dy * h;
            }

            // f_new = fun(t + h, y_new)
            let f_new = fun(t + h, &y_new);

            // K[-1] = f_new
            k[N_STAGES] = f_new;

            (y_new, f_new)
        }

        pub(super) struct Rk45Solver<const N: usize> {
            t_bound: f64,
            direction: f64,
            pub(super) t: f64,
            pub(super) y: [f64; N],
            pub(super) t_old: f64,
            y_old: [f64; N],
            max_step: f64,
            rtol: f64,
            atol: f64,
            f: [f64; N],
            h_abs: f64,
            k: [[f64; N]; N_STAGES + 1],
        }

        impl<const N: usize> Rk45Solver<N> {
            pub(super) fn new(
                fun: &mut impl FnMut(f64, &[f64; N]) -> [f64; N],
                t0: f64,
                y0: [f64; N],
                t_bound: f64,
                rtol: Option<f64>,
                atol: Option<f64>,
            ) -> Self {
                let f = fun(t0, &y0);
                let max_step = f64::INFINITY;
                let direction = if t_bound != t0 {
                    sign(t_bound - t0)
                } else {
                    1.
                };
                let rtol = rtol.unwrap_or(1e-3);
                let atol = atol.unwrap_or(1e-6);
                let h_abs = select_initial_step(
                    fun,
                    t0,
                    &y0,
                    t_bound,
                    max_step,
                    &f,
                    direction,
                    ERROR_ESTIMATOR_ORDER,
                    rtol,
                    atol,
                );

                Self {
                    t_bound,
                    direction,
                    t: t0,
                    y: y0,
                    t_old: t0,
                    y_old: y0,
                    max_step,
                    rtol,
                    atol,
                    f,
                    h_abs,
                    k: [[0.; N]; N_STAGES + 1],
                }
            }

            /// Advance one step, returning whether the solver has finished (reached t_bound).
            pub(super) fn step(
                &mut self,
                fun: &mut impl FnMut(f64, &[f64; N]) -> [f64; N],
            ) -> Result<bool, Rk45SolverError> {
                if N == 0 || self.t == self.t_bound {
                    // Handle corner cases of empty solver or no integration.
                    self.t_old = self.t;
                    self.t = self.t_bound;

                    return Ok(true);
                }

                let t = self.t;
                self.step_impl(fun)?;
                self.t_old = t;

                Ok(self.direction * (self.t - self.t_bound) >= 0.)
            }

            fn step_impl(
                &mut self,
                fun: &mut impl FnMut(f64, &[f64; N]) -> [f64; N],
            ) -> Result<(), Rk45SolverError> {
                let t = self.t;
                let y = self.y;

                let min_step = 10. * f64::abs(nextafter(t, self.direction * f64::INFINITY) - t);

                let mut h_abs = if self.h_abs > self.max_step {
                    self.max_step
                } else if self.h_abs < min_step {
                    min_step
                } else {
                    self.h_abs
                };

                let mut step_rejected = false;

                let (t_new, y_new, f_new) = loop {
                    if h_abs < min_step {
                        return Err(Rk45SolverError::TooSmallStep);
                    }

                    let mut t_new = t + h_abs * self.direction;

                    if self.direction * (t_new - self.t_bound) > 0. {
                        t_new = self.t_bound;
                    }

                    let h = t_new - t;
                    h_abs = f64::abs(h);

                    let (y_new, f_new) = rk_step(fun, t, &y, &self.f, h, &mut self.k);

                    let error_norm = self.estimate_error_norm(h, &y, &y_new);

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

                        break (t_new, y_new, f_new);
                    } else {
                        h_abs *= MIN_FACTOR.max(SAFETY * error_norm.powf(ERROR_EXPONENT));
                        step_rejected = true;
                    }
                };

                self.y_old = y;

                self.t = t_new;
                self.y = y_new;

                self.h_abs = h_abs;
                self.f = f_new;

                Ok(())
            }

            /// norm(np.dot(K.T, E) * h / scale)
            fn estimate_error_norm(&self, h: f64, y: &[f64; N], y_new: &[f64; N]) -> f64 {
                super::norm((0..N).map(|i| {
                    let scale = self.atol + y[i].abs().max(y_new[i].abs()) * self.rtol;
                    let mut error = 0.;
                    for (k_row, e) in self.k.iter().zip(E) {
                        error += k_row[i] * e;
                    }
                    error * h / scale
                }))
            }

            pub(super) fn dense_output(&self) -> RkDenseOutput<N> {
                // Q = K.T.dot(P)
                let mut q = [[0.; 4]; N];
                for (i, q_row) in q.iter_mut().enumerate() {
                    for (col, q_val) in q_row.iter_mut().enumerate() {
                        for (k_row, p_row) in self.k.iter().zip(P) {
                            *q_val += k_row[i] * p_row[col];
                        }
                    }
                }

                RkDenseOutput {
                    t_old: self.t_old,
                    h: self.t - self.t_old,
                    q,
                    y_old: self.y_old,
                    y: self.y,
                }
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

        pub(super) struct RkDenseOutput<const N: usize> {
            t_old: f64,
            h: f64,
            q: [[f64; 4]; N],
            y_old: [f64; N],
            y: [f64; N],
        }

        impl<const N: usize> RkDenseOutput<N> {
            pub(super) fn call(&self, t: f64) -> [f64; N] {
                if self.h == 0. {
                    // equivalent of ConstantDenseOutput, used when no step has been taken
                    return self.y;
                }

                let x = (t - self.t_old) / self.h;

                // p = np.cumprod(np.tile(x, order + 1))
                let x2 = x * x;
                let x3 = x2 * x;
                let p = [x, x2, x3, x3 * x];

                std::array::from_fn(|i| {
                    let q = &self.q[i];
                    self.h * (q[0] * p[0] + q[1] * p[1] + q[2] * p[2] + q[3] * p[3]) + self.y_old[i]
                })
            }
        }

        fn sign(value: f64) -> f64 {
            if value < 0. {
                -1.
            } else if value == 0. {
                0.
            } else {
                1.
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn select_initial_step<const N: usize>(
        fun: &mut impl FnMut(f64, &[f64; N]) -> [f64; N],
        t0: f64,
        y0: &[f64; N],
        t_bound: f64,
        max_step: f64,
        f0: &[f64; N],
        direction: f64,
        order: usize,
        rtol: f64,
        atol: f64,
    ) -> f64 {
        if N == 0 {
            return f64::INFINITY;
        }

        let interval_length = (t_bound - t0).abs();
        if interval_length == 0.0 {
            return 0.0;
        }

        let scale: [f64; N] = std::array::from_fn(|i| (y0[i].abs() * rtol) + atol);

        let d0 = norm((0..N).map(|i| y0[i] / scale[i]));
        let d1 = norm((0..N).map(|i| f0[i] / scale[i]));

        let h0 = if d0 < 1e-5 || d1 < 1e-5 {
            1e-6
        } else {
            0.01 * d0 / d1
        };

        // Check t0+h0*direction doesn't take us beyond t_bound
        let h0 = h0.min(interval_length);
        let y1: [f64; N] = std::array::from_fn(|i| y0[i] + h0 * direction * f0[i]);
        let f1 = fun(t0 + h0 * direction, &y1);
        let d2 = norm((0..N).map(|i| (f1[i] - f0[i]) / scale[i])) / h0;

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
    fn norm(values: impl ExactSizeIterator<Item = f64>) -> f64 {
        let len = values.len();
        values.map(|x| x * x).sum::<f64>().sqrt() / (len as f64).powf(0.5)
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

/// Equivalent of `interp::interp` with `InterpMode::Extrapolate`, but with the slopes and
/// intercepts precomputed so that repeated evaluation (e.g. inside an ODE function) does not
/// allocate.
#[derive(Clone, Debug)]
pub struct LinearInterpolator {
    x: Vec<f64>,
    y_limits: (f64, f64),
    slopes: Vec<f64>,
    intercepts: Vec<f64>,
}

impl LinearInterpolator {
    pub fn new(x: &[f64], y: &[f64]) -> Self {
        // as with `interp`, ignore any excess entries in either x or y
        let min_len = x.len().min(y.len());

        let slopes: Vec<f64> = (1..min_len)
            .map(|i| (y[i] - y[i - 1]) / (x[i] - x[i - 1]))
            .collect();
        let intercepts = slopes
            .iter()
            .enumerate()
            .map(|(i, m)| y[i] - x[i] * m)
            .collect();
        let y_limits = match min_len {
            0 => (0., 0.),
            _ => (y[0], y[min_len - 1]),
        };

        Self {
            x: x[..min_len].to_vec(),
            y_limits,
            slopes,
            intercepts,
        }
    }

    pub fn interp(&self, xp: f64) -> f64 {
        match self.x.len() {
            0 => 0.,
            1 => self.y_limits.0,
            len => {
                if xp == self.x[0] {
                    return self.y_limits.0;
                }
                if xp == self.x[len - 1] {
                    return self.y_limits.1;
                }

                let i = self
                    .x
                    .partition_point(|&probe| probe < xp)
                    .saturating_sub(1)
                    .min(len - 2);

                self.slopes[i] * xp + self.intercepts[i]
            }
        }
    }
}
