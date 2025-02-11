// TODO this is from scipy
// Find equivalent function in a Rust library or implement
pub(crate) fn fsolve<const ARGCOUNT: usize>(
    _func: impl FnOnce(f64, [f64; ARGCOUNT]) -> f64,
    x0: f64,
    _args: [f64; ARGCOUNT],
) -> f64 {
    // Stub implementation for the time being
    x0
}
