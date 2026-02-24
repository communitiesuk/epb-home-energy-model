use roots::{find_root_brent, SimpleConvergency};

// TODO this is from scipy
// Find equivalent function in a Rust library or implement
pub(crate) fn fsolve<const ARGCOUNT: usize>(
    _func: impl Fn(f64, [f64; ARGCOUNT]) -> anyhow::Result<f64>,
    x0: f64,
    _args: [f64; ARGCOUNT],
) -> anyhow::Result<f64> {
    // For the time being we use the other root solver we have
    // TODO implement fsolve
    Ok(x0) // stub response
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
