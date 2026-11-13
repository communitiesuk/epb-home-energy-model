/// A simple statistics module with some utility functions such as calculation of percentiles.
use interp::{interp, InterpMode};

use std::sync::atomic::AtomicUsize;
#[cfg(feature = "cache-np-interp")]
use std::{collections::HashMap, sync::LazyLock};
#[cfg(feature = "cache-np-interp")]
use itertools::Itertools;
#[cfg(feature = "cache-np-interp")]
use ordered_float::OrderedFloat;
#[cfg(feature = "cache-np-interp")]
use parking_lot::RwLock;

#[allow(dead_code)]
pub fn percentile(mut numbers: Vec<f64>, percentile: usize) -> f64 {
    // use implementation from rust test::stats (currently unstable)
    // https://github.com/rust-lang/rust/blob/main/library/test/src/stats.rs#L258-L280
    numbers.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let sorted_samples = numbers;
    if sorted_samples.len() == 1 {
        return sorted_samples[0];
    }
    let zero: f64 = 0.0;
    let pct = percentile as f64;
    assert!(zero <= pct);
    let hundred = 100_f64;
    assert!(pct <= hundred);
    if pct == hundred {
        return sorted_samples[sorted_samples.len() - 1];
    }
    let length = (sorted_samples.len() - 1) as f64;
    let rank = (pct / hundred) * length;
    let lrank = rank.floor();
    let d = rank - lrank;
    let n = lrank as usize;
    let lo = sorted_samples[n];
    let hi = sorted_samples[n + 1];
    lo + (hi - lo) * d
}

#[cfg(feature = "cache-np-interp")]
#[derive(Eq, PartialEq, Clone, Debug, Hash)]
struct CacheKey {
    input: OrderedFloat<f64>,
    x: Vec<OrderedFloat<f64>>,
    y: Vec<OrderedFloat<f64>>,
}
pub static CACHE_HITS: AtomicUsize = AtomicUsize::new(0);
pub static CACHE_MISSES: AtomicUsize = AtomicUsize::new(0);

#[cfg(feature = "cache-np-interp")]
static CACHE: LazyLock<RwLock<HashMap<CacheKey, f64>>> = LazyLock::new(|| {
    RwLock::new(HashMap::new())
});

pub fn np_interp(input: f64, x: &[f64], y: &[f64]) -> f64 {
    // dbg!("np_interp called with input:", input, "x:", x, "y:", y);

    // TODO rounding floats could be more performant to make cache keys
    // can we just make this into a string 

    #[cfg(feature = "cache-np-interp")]
    let cache_key = CacheKey {
        input: OrderedFloat(input),
        x: x.iter().map(|&f| OrderedFloat(f)).collect_vec(),
        y: y.iter().map(|&f| OrderedFloat(f)).collect_vec(),
    };

    #[cfg(feature = "cache-np-interp")]
    {
        // scope this to avoid a deadlock
        let cache = CACHE.read();
        if let Some(&result) = cache.get(&cache_key) {
            CACHE_HITS.fetch_add(1usize, std::sync::atomic::Ordering::Relaxed);
            return result;
        }
    }

    let result = np_interp_impl(input, x, y);

    #[cfg(feature = "cache-np-interp")]
    CACHE_MISSES.fetch_add(1usize, std::sync::atomic::Ordering::Relaxed);

    #[cfg(feature = "cache-np-interp")]
    CACHE.write().insert(cache_key, result);
    result
}

/// This function matches the behaviour Numpy interp
/// https://numpy.org/doc/stable/reference/generated/numpy.interp.html
pub fn np_interp_impl(input: f64, x: &[f64], y: &[f64]) -> f64 {

    if x.is_empty() {
        panic!("x cannot be empty");
    }

    if y.is_empty() {
        panic!("y cannot be empty");
    }

    if x.len() != y.len() {
        panic!("x and y must be of equal length");
    }

    interp(x, y, input, &InterpMode::FirstLast)
}

/// This function matches the behaviour Numpy interp
/// https://numpy.org/doc/stable/reference/generated/numpy.interp.html
pub fn np_interp_with_extrapolate(input: f64, x: &[f64], y: &[f64]) -> f64 {

    if x.is_empty() {
        panic!("x cannot be empty");
    }

    if y.is_empty() {
        panic!("y cannot be empty");
    }

    if x.len() != y.len() {
        panic!("x and y must be of equal length");
    }

    interp(x, y, input, &InterpMode::Extrapolate)
}

pub fn linspace(start: f64, end: f64, num: i32) -> Vec<f64> {
    let step = (end - start) / f64::from(num - 1);
    (0..num).map(|n| start + (f64::from(n) * step)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::{assert_relative_eq, assert_ulps_eq};
    use interp::interp;
    use rstest::*;

    #[fixture]
    fn numbers() -> Vec<f64> {
        vec![9.0, 3.0, 3.0, 4.0, 5.0, 4.9, 8.0, 3.3, 2.0, 0.1]
    }

    #[fixture]
    fn other_numbers() -> Vec<f64> {
        vec![14.2, 16.5, 3.4, 7.8, 18.4, 10.5, 7.4, 2.9, 22.3, 16.0]
    }

    #[rstest]
    fn test_percentile(numbers: Vec<f64>) {
        assert_relative_eq!(percentile(numbers.clone(), 70), 4.95, max_relative = 1e-2);
        assert_relative_eq!(percentile(numbers, 50), 3.65, max_relative = 1e-2);
    }

    #[rstest]
    fn test_percentile_with_other_cases(other_numbers: Vec<f64>) {
        assert_ulps_eq!(percentile(other_numbers, 80), 16.88);
    }

    #[rstest]
    #[case(0.)]
    #[case(1.)]
    #[case(1.5)]
    #[case(5.)]
    fn test_np_interp_given_number_in_expected_range(#[case] input: f64) {
        let x = [0., 1., 2., 3., 4., 5.];
        let y = [0., 10., 20., 30., 40., 50.];

        let actual = np_interp(input, &x, &y);
        let expected = interp(&x, &y, input, &InterpMode::Extrapolate);

        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case(5.1, 50.)]
    #[case(10., 50.)]
    fn test_np_interp_given_number_above_expected_range_returns_last_element_of_y(
        #[case] input: f64,
        #[case] expected: f64,
    ) {
        let x = [0., 1., 2., 3., 4., 5.];
        let y = [0., 10., 20., 30., 40., 50.];

        let actual = np_interp(input, &x, &y);

        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case(5.1, 0.)]
    #[case(10., 0.)]
    fn test_np_interp_given_number_above_expected_range_returns_last_element_of_y_for_inverse(
        #[case] input: f64,
        #[case] expected: f64,
    ) {
        let x = [0., 1., 2., 3., 4., 5.];
        let y = [50., 40., 30., 20., 10., 0.];

        let actual = np_interp(input, &x, &y);

        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case(-0.1, 0.)]
    #[case(-25., 0.)]
    fn test_np_interp_given_number_below_expected_range_returns_first_element_of_y(
        #[case] input: f64,
        #[case] expected: f64,
    ) {
        let x = [0., 1., 2., 3., 4., 5.];
        let y = [0., 10., 20., 30., 40., 50.];

        let actual = np_interp(input, &x, &y);

        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case(-0.1, 50.)]
    #[case(-25., 50.)]
    fn test_np_interp_given_number_below_expected_range_returns_first_element_of_y_for_inverse(
        #[case] input: f64,
        #[case] expected: f64,
    ) {
        let x = [0., 1., 2., 3., 4., 5.];
        let y = [50., 40., 30., 20., 10., 0.];

        let actual = np_interp(input, &x, &y);
        assert_eq!(actual, expected);
    }

    #[test]
    #[should_panic(expected = "x cannot be empty")]
    fn test_np_interp_empty_x_panics() {
        let x = [];
        let y = [1., 2., 3., 4., 5.];

        np_interp(3., &x, &y);
    }

    #[test]
    #[should_panic(expected = "y cannot be empty")]
    fn test_np_interp_empty_y_panics() {
        let x = [1., 2., 3., 4., 5.];
        let y = [];

        np_interp(3., &x, &y);
    }

    #[test]
    #[should_panic(expected = "x and y must be of equal length")]
    fn test_np_interp_different_length_x_and_y_panics() {
        let x = [1., 2., 3., 4., 5.];
        let y = [10., 20.];

        np_interp(3., &x, &y);
    }
}
