use roots::{find_root_brent, SimpleConvergency};

// TODO this is from scipy
// Find equivalent function in a Rust library or implement
pub(crate) fn fsolve<const ARGCOUNT: usize>(
    _func: impl FnOnce(f64, [f64; ARGCOUNT]) -> anyhow::Result<f64>,
    x0: f64,
    _args: [f64; ARGCOUNT],
) -> anyhow::Result<f64> {
    // Stub implementation for the time being
    Ok(x0)
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

    find_root_brent::<f64, _>(
        (x0 - guess_interval).max(0.), // ensure first bracket is at least zero
        x0 + guess_interval,
        |f| fun(f, args),
        &mut convergency,
    )
    .map_err(|e| anyhow::anyhow!(e))
}
