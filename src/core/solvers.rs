use anyhow::anyhow;
use argmin::core::CostFunction;
use differential_equations::prelude::*;
use eqsolver::single_variable::FDNewton;
use interp::{interp_slice, InterpMode};
use itertools::Itertools;
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

/// direct port of scipy.integrate.solve_ivp but only including areas used in this library
mod solve_ivp {
    use ndarray::Array1;
    use ndarray_linalg::Norm;
    use std::sync::Arc;

    pub fn solve_ivp(
        func: &Arc<dyn Fn(f64, &Array1<f64>) -> Array1<f64> + Send + Sync>,
        t_span: (f64, f64),
        y0: &Array1<f64>,
        events: Option<&[TerminatingEvent]>,
        rtol: Option<f64>,
        atol: Option<f64>,
    ) -> anyhow::Result<OdeResult> {
        let (t0, tf) = t_span;

        let solver = rk45::RK45Solver::new(func, t0, y0, tf);

        todo!()
    }

    //
    //     if method in METHODS:
    //         method = METHODS[method]
    //
    //     solver = method(fun, t0, y0, tf, vectorized=vectorized, **options)
    //
    //     if t_eval is None:
    //         ts = [t0]
    //         ys = [y0]
    //     elif t_eval is not None and dense_output:
    //         ts = []
    //         ti = [t0]
    //         ys = []
    //     else:
    //         ts = []
    //         ys = []
    //
    //     interpolants = []
    //
    //     if events is not None:
    //         events, max_events, event_dir = prepare_events(events)
    //         event_count = np.zeros(len(events))
    //         if args is not None:
    //             # Wrap user functions in lambdas to hide the additional parameters.
    //             # The original event function is passed as a keyword argument to the
    //             # lambda to keep the original function in scope (i.e., avoid the
    //             # late binding closure "gotcha").
    //             events = [lambda t, x, event=event: event(t, x, *args)
    //                       for event in events]
    //         g = [event(t0, y0) for event in events]
    //         t_events = [[] for _ in range(len(events))]
    //         y_events = [[] for _ in range(len(events))]
    //     else:
    //         t_events = None
    //         y_events = None
    //
    //     status = None
    //     while status is None:
    //         message = solver.step()
    //
    //         if solver.status == 'finished':
    //             status = 0
    //         elif solver.status == 'failed':
    //             status = -1
    //             break
    //
    //         t_old = solver.t_old
    //         t = solver.t
    //         y = solver.y
    //
    //         if dense_output:
    //             sol = solver.dense_output()
    //             interpolants.append(sol)
    //         else:
    //             sol = None
    //
    //         if events is not None:
    //             g_new = [event(t, y) for event in events]
    //             active_events = find_active_events(g, g_new, event_dir)
    //             if active_events.size > 0:
    //                 if sol is None:
    //                     sol = solver.dense_output()
    //
    //                 event_count[active_events] += 1
    //                 root_indices, roots, terminate = handle_events(
    //                     sol, events, active_events, event_count, max_events,
    //                     t_old, t)
    //
    //                 for e, te in zip(root_indices, roots):
    //                     t_events[e].append(te)
    //                     y_events[e].append(sol(te))
    //
    //                 if terminate:
    //                     status = 1
    //                     t = roots[-1]
    //                     y = sol(t)
    //
    //             g = g_new
    //
    //         if t_eval is None:
    //             ts.append(t)
    //             ys.append(y)
    //         else:
    //             # The value in t_eval equal to t will be included.
    //             if solver.direction > 0:
    //                 t_eval_i_new = np.searchsorted(t_eval, t, side='right')
    //                 t_eval_step = t_eval[t_eval_i:t_eval_i_new]
    //             else:
    //                 t_eval_i_new = np.searchsorted(t_eval, t, side='left')
    //                 # It has to be done with two slice operations, because
    //                 # you can't slice to 0th element inclusive using backward
    //                 # slicing.
    //                 t_eval_step = t_eval[t_eval_i_new:t_eval_i][::-1]
    //
    //             if t_eval_step.size > 0:
    //                 if sol is None:
    //                     sol = solver.dense_output()
    //                 ts.append(t_eval_step)
    //                 ys.append(sol(t_eval_step))
    //                 t_eval_i = t_eval_i_new
    //
    //         if t_eval is not None and dense_output:
    //             ti.append(t)
    //
    //     message = MESSAGES.get(status, message)
    //
    //     if t_events is not None:
    //         t_events = [np.asarray(te) for te in t_events]
    //         y_events = [np.asarray(ye) for ye in y_events]
    //
    //     if t_eval is None:
    //         ts = np.array(ts)
    //         ys = np.vstack(ys).T
    //     elif ts:
    //         ts = np.hstack(ts)
    //         ys = np.hstack(ys)
    //
    //     if dense_output:
    //         if t_eval is None:
    //             sol = OdeSolution(
    //                 ts, interpolants, alt_segment=True if method in [BDF, LSODA] else False
    //             )
    //         else:
    //             sol = OdeSolution(
    //                 ti, interpolants, alt_segment=True if method in [BDF, LSODA] else False
    //             )
    //     else:
    //         sol = None
    //
    //     return OdeResult(t=ts, y=ys, sol=sol, t_events=t_events, y_events=y_events,
    //                      nfev=solver.nfev, njev=solver.njev, nlu=solver.nlu,
    //                      status=status, message=message, success=status >= 0)

    pub struct TerminatingEvent {
        func: Arc<dyn Fn(f64, &[f64]) -> f64 + Send + Sync>,
        terminal: bool,
        direction: TerminateDirection,
    }

    impl TerminatingEvent {
        pub fn new(
            func: Arc<dyn Fn(f64, &[f64]) -> f64 + Send + Sync>,
            terminal: bool,
            direction: Option<TerminateDirection>,
        ) -> Self {
            Self {
                func,
                terminal,
                direction: direction.unwrap_or_default(),
            }
        }
    }

    #[derive(Clone, Copy, Debug, Default)]
    pub enum TerminateDirection {
        #[default]
        Both,
        Positive,
        Negative,
    }

    pub struct OdeResult {
        pub y: Vec<Vec<f64>>,
        pub t_events: Vec<Vec<f64>>,
        pub t: Vec<f64>,
    }

    /// module to do with the RK45 method for scipy
    mod rk45 {
        use crate::core::solvers::solve_ivp::select_initial_step;
        use ndarray::{array, Array1, Array2};
        use std::sync::{Arc, LazyLock};
        use thiserror::Error;

        fn rk_step(
            fun: &Arc<dyn Fn(f64, &Array1<f64>) -> Array1<f64> + Send + Sync>,
            t: f64,
            y: &Array1<f64>,
            f: &Array1<f64>,
            h: f64,
            a: &Array2<f64>,
            b: &Array1<f64>,
            c: &Array1<f64>,
            k: &mut Array2<f64>,
        ) -> (Array2<f64>, Array2<f64>) {
            // what's the Rust equiv here???
            // k[0] = f.clone();

            //     K[0] = f
            //     for s, (a, c) in enumerate(zip(A[1:], C[1:]), start=1):
            //         dy = np.dot(K[:s].T, a[:s]) * h
            //         K[s] = fun(t + c * h, y + dy)
            //
            //     y_new = y + h * np.dot(K[:-1].T, B)
            //     f_new = fun(t + h, y_new)
            //
            //     K[-1] = f_new

            todo!()
        }

        #[derive(Clone)]
        pub(super) struct RK45Solver {
            n: usize,
            status: super::base_solver::Status,
            t_bound: f64,
            direction: i8,
            t: f64,
            y: Array1<f64>,
            t_old: Option<f64>,
            step_size: Option<f64>,
            nfev: usize,
            fun: Arc<dyn Fn(f64, &Array1<f64>) -> Array1<f64> + Send + Sync>,
            y_old: Option<()>,
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

        impl RK45Solver {
            pub(super) fn new(
                fun: &Arc<dyn Fn(f64, &Array1<f64>) -> Array1<f64> + Send + Sync>,
                t0: f64,
                y0: &Array1<f64>,
                t_bound: f64,
            ) -> Self {
                let f = fun(t0, y0);
                let t = t0;
                let y = y0.to_owned();
                let max_step = f64::INFINITY;
                let direction = if t_bound != t0 { sign(t_bound - t0) } else { 1 };
                let rtol = 1e-3;
                let atol = 1e-6;
                let h_abs = select_initial_step(
                    fun,
                    t0,
                    y0,
                    t_bound,
                    max_step,
                    &f,
                    direction as f64,
                    ORDER,
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
                if self.t_old.is_none() {
                    None
                } else {
                    Some((self.t - self.t_old.unwrap()).abs())
                }
            }

            fn step(&mut self) -> Result<(), Rk45SolverError> {
                if self.status != crate::core::solvers::solve_ivp::base_solver::Status::Running {
                    panic!("Attempt to step on a failed or finished solver.");
                }

                if self.n == 0 || self.t == self.t_bound {
                    // Handle corner cases of empty solver or no integration.
                    self.t_old = self.t.into();
                    self.t = self.t_bound;
                    self.status = crate::core::solvers::solve_ivp::base_solver::Status::Finished;

                    Ok(())
                } else {
                    let t = self.t;
                    let result = self.step_impl();
                    if result.is_err() {
                        self.status = crate::core::solvers::solve_ivp::base_solver::Status::Failed;
                    } else {
                        self.t_old = t.into();
                        if self.direction as f64 * (self.t - self.t_bound) >= 0. {
                            self.status =
                                crate::core::solvers::solve_ivp::base_solver::Status::Finished;
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

                while !step_accepted {
                    if h_abs < min_step {
                        return Err(Rk45SolverError::TooSmallStep);
                    }

                    let h = h_abs * self.direction as f64;
                    let mut t_new = t + h;

                    if self.direction as f64 * (t_new - self.t_bound) > 0. {
                        t_new = self.t_bound;
                    }

                    let h = t_new - t;
                    h_abs = f64::abs(h);

                    let (y_new, f_new) =
                        rk_step(&self.fun, t, y, &self.f, h, &A, &B, &C, &mut self.k);

                    // implement rk_step
                }

                todo!()

                //         //
                //         //         while not step_accepted:
                //         //             if h_abs < min_step:
                //         //                 return False, self.TOO_SMALL_STEP
                //         //
                //         //             h = h_abs * self.direction
                //         //             t_new = t + h
                //         //
                //         //             if self.direction * (t_new - self.t_bound) > 0:
                //         //                 t_new = self.t_bound
                //         //
                //         //             h = t_new - t
                //         //             h_abs = np.abs(h)
                //         //
                //         //             y_new, f_new = rk_step(self.fun, t, y, self.f, h, self.A,
                //         //                                    self.B, self.C, self.K)
                //         //             scale = atol + np.maximum(np.abs(y), np.abs(y_new)) * rtol
                //         //             error_norm = self._estimate_error_norm(self.K, h, scale)
                //         //
                //         //             if error_norm < 1:
                //         //                 if error_norm == 0:
                //         //                     factor = MAX_FACTOR
                //         //                 else:
                //         //                     factor = min(MAX_FACTOR,
                //         //                                  SAFETY * error_norm ** self.error_exponent)
                //         //
                //         //                 if step_rejected:
                //         //                     factor = min(1, factor)
                //         //
                //         //                 h_abs *= factor
                //         //
                //         //                 step_accepted = True
                //         //             else:
                //         //                 h_abs *= max(MIN_FACTOR,
                //         //                              SAFETY * error_norm ** self.error_exponent)
                //         //                 step_rejected = True
                //         //
                //         //         self.h_previous = h
                //         //         self.y_old = y
                //         //
                //         //         self.t = t_new
                //         //         self.y = y_new
                //         //
                //         //         self.h_abs = h_abs
                //         //         self.f = f_new
                //         //
                //         //         return True, None
            }

            fn dense_output(&self) -> impl super::base_solver::DenseOutput {
                if self.n == 0 || self.t_old.is_some_and(|t_old| self.t == t_old) {
                    crate::core::solvers::solve_ivp::base_solver::ConstantDenseOutput::new(
                        self.t_old
                            .expect("Expected to be past the first step with t_old set"),
                        self.t,
                        self.y.clone(),
                    )
                } else {
                    unimplemented!()
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

        // class RungeKutta(OdeSolver):

        //
        //     def _estimate_error(self, K, h):
        //         return np.dot(K.T, self.E) * h
        //
        //     def _estimate_error_norm(self, K, h, scale):
        //         return norm(self._estimate_error(K, h) / scale)
        //
        //     def _step_impl(self):
        //         t = self.t
        //         y = self.y
        //
        //         max_step = self.max_step
        //         rtol = self.rtol
        //         atol = self.atol
        //
        //         min_step = 10 * np.abs(np.nextafter(t, self.direction * np.inf) - t)
        //
        //         if self.h_abs > max_step:
        //             h_abs = max_step
        //         elif self.h_abs < min_step:
        //             h_abs = min_step
        //         else:
        //             h_abs = self.h_abs
        //
        //         step_accepted = False
        //         step_rejected = False
        //
        //         while not step_accepted:
        //             if h_abs < min_step:
        //                 return False, self.TOO_SMALL_STEP
        //
        //             h = h_abs * self.direction
        //             t_new = t + h
        //
        //             if self.direction * (t_new - self.t_bound) > 0:
        //                 t_new = self.t_bound
        //
        //             h = t_new - t
        //             h_abs = np.abs(h)
        //
        //             y_new, f_new = rk_step(self.fun, t, y, self.f, h, self.A,
        //                                    self.B, self.C, self.K)
        //             scale = atol + np.maximum(np.abs(y), np.abs(y_new)) * rtol
        //             error_norm = self._estimate_error_norm(self.K, h, scale)
        //
        //             if error_norm < 1:
        //                 if error_norm == 0:
        //                     factor = MAX_FACTOR
        //                 else:
        //                     factor = min(MAX_FACTOR,
        //                                  SAFETY * error_norm ** self.error_exponent)
        //
        //                 if step_rejected:
        //                     factor = min(1, factor)
        //
        //                 h_abs *= factor
        //
        //                 step_accepted = True
        //             else:
        //                 h_abs *= max(MIN_FACTOR,
        //                              SAFETY * error_norm ** self.error_exponent)
        //                 step_rejected = True
        //
        //         self.h_previous = h
        //         self.y_old = y
        //
        //         self.t = t_new
        //         self.y = y_new
        //
        //         self.h_abs = h_abs
        //         self.f = f_new
        //
        //         return True, None
        //
        //     def _dense_output_impl(self):
        //         Q = self.K.T.dot(self.P)
        //         return RkDenseOutput(self.t_old, self.t, self.y_old, Q)
        //

        // class RkDenseOutput(DenseOutput):
        //     def __init__(self, t_old, t, y_old, Q):
        //         super().__init__(t_old, t)
        //         self.h = t - t_old
        //         self.Q = Q
        //         self.order = Q.shape[1] - 1
        //         self.y_old = y_old
        //
        //     def _call_impl(self, t):
        //         x = (t - self.t_old) / self.h
        //         if t.ndim == 0:
        //             p = np.tile(x, self.order + 1)
        //             p = np.cumprod(p)
        //         else:
        //             p = np.tile(x, (self.order + 1, 1))
        //             p = np.cumprod(p, axis=0)
        //         y = self.h * np.dot(self.Q, p)
        //         if y.ndim == 2:
        //             y += self.y_old[:, None]
        //         else:
        //             y += self.y_old
        //
        //         return y

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
        fun: &Arc<dyn Fn(f64, &Array1<f64>) -> Array1<f64> + Send + Sync>,
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
        if y0.len() == 0 {
            return f64::INFINITY;
        }

        let interval_length = (t_bound - t0).abs();
        if interval_length == 0.0 {
            return 0.0;
        }

        let scale: Array1<f64> = y0.iter().map(|&y| (y.abs() * rtol) + atol).collect();

        let d0 = (y0 / &scale).norm();
        let d1 = (f0 / &scale).norm();

        let h0 = if d0 < 1e-5 || d1 < 1e-5 {
            1e-6
        } else {
            0.01 * d0 / d1
        };

        // Check t0+h0*direction doesn't take us beyond t_bound
        let h0 = h0.min(interval_length);
        let y1 = y0 + h0 * direction * f0;
        let f1 = fun(t0 + h0 * direction, &y1);
        let d2 = ((f1 - f0) / &scale).norm() / h0;

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
            fn call_impl(&self, t: f64) -> f64;
        }

        pub struct DefaultDenseOutput {
            t_old: f64,
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

        fn dense_output_init(t_old: f64, t: f64) -> DefaultDenseOutput {
            DefaultDenseOutput {
                t_old,
                t,
                t_min: t.min(t_old),
                t_max: t.max(t_old),
            }
        }

        impl DenseOutput for ConstantDenseOutput {
            fn call_impl(&self, t: f64) -> f64 {
                //         if t.ndim == 0:
                //             return self.value
                //         else:
                //             ret = np.empty((self.value.shape[0], t.shape[0]))
                //             ret[:] = self.value[:, None]
                //             return ret
                todo!()
            }
        }

        //
        //     def __call__(self, t):
        //         """Evaluate the interpolant.
        //
        //         Parameters
        //         ----------
        //         t : float or array_like with shape (n_points,)
        //             Points to evaluate the solution at.
        //
        //         Returns
        //         -------
        //         y : ndarray, shape (n,) or (n, n_points)
        //             Computed values. Shape depends on whether `t` was a scalar or a
        //             1-D array.
        //         """
        //         t = np.asarray(t)
        //         if t.ndim > 1:
        //             raise ValueError("`t` must be a float or a 1-D array.")
        //         return self._call_impl(t)
        //
        //     def _call_impl(self, t):
        //         raise NotImplementedError
        //
        //
        // class ConstantDenseOutput(DenseOutput):
        //     """Constant value interpolator.
        //
        //     This class used for degenerate integration cases: equal integration limits
        //     or a system with 0 equations.
        //     """
        //     def __init__(self, t_old, t, value):
        //         super().__init__(t_old, t)
        //         self.value = value
        //
        //     def _call_impl(self, t):
        //         if t.ndim == 0:
        //             return self.value
        //         else:
        //             ret = np.empty((self.value.shape[0], t.shape[0]))
        //             ret[:] = self.value[:, None]
        //             return ret
    }
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

    let system = OdeConfiguration { func: func.clone() };

    let solver = ExplicitRungeKutta::dopri5().rtol(rtol).atol(atol);

    let (t0, tf) = t_span;

    let ode_events: Vec<OdeEvent> = events
        .map(|events| {
            events
                .funcs
                .iter()
                .map(|func| OdeEvent {
                    func: func.clone(),
                    direction: events.direction.unwrap_or_default().into(),
                })
                .collect_vec()
        })
        .unwrap_or_default();

    let solution = match &ode_events.as_slice() {
        &[] => IVP::ode(&system, t0, tf, y0.to_vec())
            .method(solver)
            .solve()
            .map_err(|e| anyhow!("IVP solver failed: {:?}", e))?,
        &[event1, event2] => IVP::ode(&system, t0, tf, y0.to_vec())
            .method(solver)
            .event(event1)
            .event(event2)
            .solve()
            .map_err(|e| anyhow!("IVP solver failed: {:?}", e))?,
        &[event] => IVP::ode(&system, t0, tf, y0.to_vec())
            .method(solver)
            .event(event)
            .solve()
            .map_err(|e| anyhow!("IVP solver failed: {:?}", e))?,
        _ => unimplemented!("The solve_ivp function does not support more than two events at the same time (though could be extended if needed)"),
    };

    // t_event is equivalent digest of the scipy.optimize.t_events array but with just the final value, and only provided if there was a termination i.e. the solution was interrupted
    let t_event: Option<f64> = (solution.status == Status::Interrupted && solution.t.len() > 0)
        .then(|| solution.t.last().copied().unwrap());

    Ok(OdeResult {
        y: transpose(solution.y),
        t_event,
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

struct OdeConfiguration {
    func: Arc<dyn Fn(f64, &[f64]) -> Vec<f64> + Send + Sync>,
}

type OdeRealNumber = f64;
type RealGrouping = Vec<OdeRealNumber>;

impl ODE<OdeRealNumber, RealGrouping> for OdeConfiguration {
    fn diff(&self, t: f64, y: &Vec<f64>, dydt: &mut Vec<f64>) {
        let _ = std::mem::replace(dydt, (self.func)(t, y));
    }
}

struct OdeEvent {
    direction: CrossingDirection,
    func: Arc<dyn Fn(f64, &[f64]) -> f64 + Send + Sync>,
}

impl Event<OdeRealNumber, RealGrouping> for OdeEvent {
    fn config(&self) -> EventConfig {
        EventConfig::new(self.direction, None).terminal()
    }

    fn event(&self, t: OdeRealNumber, y: &RealGrouping) -> OdeRealNumber {
        // calculate time event may have occurred
        //

        (self.func)(t, y)
    }
}

struct EventTimeProblem(Arc<dyn Fn(f64, &[f64]) -> f64 + Send + Sync>);

impl CostFunction for EventTimeProblem {
    type Param = f64;
    type Output = f64;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, anyhow::Error> {
        Ok((self.0)(*param, &[]))
    }
}

pub struct OdeResult {
    pub y: Vec<Vec<f64>>,
    pub t_event: Option<f64>,
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

impl From<TerminateDirection> for CrossingDirection {
    fn from(direction: TerminateDirection) -> Self {
        match direction {
            TerminateDirection::Both => CrossingDirection::Both,
            TerminateDirection::Positive => CrossingDirection::Positive,
            TerminateDirection::Negative => CrossingDirection::Negative,
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
