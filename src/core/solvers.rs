use eqsolver::single_variable::FDNewton;
use interp::{interp_slice, InterpMode};
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
pub mod solve_ivp {
    use crate::core::solvers::solve_ivp::base_solver::{DenseOutput, Status};
    use ndarray::{s, Array, Array1, Axis, Dimension, Zip};
    use roots::{find_root_brent, SimpleConvergency};
    use std::sync::Arc;
    use thiserror::Error;

    pub fn solve_ivp(
        func: &Arc<dyn Fn(f64, &Array1<f64>) -> Array1<f64> + Send + Sync>,
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
            let _result = solver
                .step()
                .map_err(|e| SolverError::MainSolverError(e.into()))?;

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
                if active_events.len() > 0 {
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

        let (t_events, _y_events): (Option<Vec<Array1<f64>>>, Option<Vec<Array1<Array1<f64>>>>) =
            if let Some(EventData {
                ref t_events,
                ref y_events,
                ..
            }) = event_data
            {
                (
                    Some(
                        t_events
                            .iter()
                            .map(|te| te.iter().map(|&te| te.into()).collect())
                            .collect(),
                    ),
                    Some(
                        y_events
                            .iter()
                            .map(|ye| ye.iter().map(|ye| ye.to_owned().into()).collect())
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

    struct EventData<'a> {
        t_events: Vec<Vec<f64>>,
        y_events: Vec<Vec<Array1<f64>>>,
        events: &'a [TerminatingEvent],
        max_events: Array1<f64>,
        event_dir: Array1<TerminateDirection>,
        event_count: Array1<usize>,
        g: Array1<f64>,
    }

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
            let order = argsort(&roots, !(t > t_old));
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
        event: &Arc<dyn Fn(f64, &Array1<f64>) -> f64 + Send + Sync>,
        sol: &Arc<dyn DenseOutput>,
        t_old: f64,
        t: f64,
    ) -> Result<f64, SolverError> {
        let f = move |t| event(t, &sol.call_impl(t));

        let mut convergency = SimpleConvergency {
            eps: f64::EPSILON,
            max_iter: 1000,
        };

        find_root_brent(t_old, t, f, &mut convergency)
            .map_err(|err| SolverError::EventError(err.into()))
    }

    #[derive(Debug, Error)]
    pub enum SolverError {
        #[error("Error when solving events: {0}")]
        EventError(Box<dyn std::error::Error + Send + Sync>),
        #[error("Error in main IVP solver: {0}")]
        MainSolverError(Box<dyn std::error::Error + Send + Sync>),
        #[error("Error happened during event termination")]
        TerminationError,
    }

    pub struct TerminatingEvent {
        func: Arc<dyn Fn(f64, &Array1<f64>) -> f64 + Send + Sync>,
        direction: TerminateDirection,
    }

    impl TerminatingEvent {
        pub fn new(
            func: Arc<dyn Fn(f64, &Array1<f64>) -> f64 + Send + Sync>,
            direction: Option<TerminateDirection>,
        ) -> Self {
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
        use crate::core::solvers::solve_ivp::base_solver::{DefaultDenseOutput, DenseOutput};
        use crate::core::solvers::solve_ivp::select_initial_step;
        use ndarray::{array, concatenate, s, Array1, Array2, Axis, Zip};
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
            fun: Arc<dyn Fn(f64, &Array1<f64>) -> Array1<f64> + Send + Sync>,
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
                fun: &Arc<dyn Fn(f64, &Array1<f64>) -> Array1<f64> + Send + Sync>,
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
                if self.t_old.is_none() {
                    None
                } else {
                    Some((self.t - self.t_old.unwrap()).abs())
                }
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
        use ndarray::{array, Array1};

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
            fn call_impl(&self, t: f64) -> Array1<f64> {
                array![t]
            }
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
