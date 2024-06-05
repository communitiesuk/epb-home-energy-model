/// A simple statistics module with some utility functions such as calculation of percentiles.
use statrs::statistics::{Data, OrderStatistics};

#[allow(dead_code)]
pub fn percentile(numbers: &[f64], percentile: usize) -> f64 {
    let numbers = numbers.iter().cloned().collect::<Vec<_>>();
    let mut data = Data::new(numbers);

    data.percentile(percentile)
}

#[cfg(test)]
mod tests {
    use super::*;
    use assert_approx_eq::assert_approx_eq;
    use rstest::*;

    #[fixture]
    fn numbers() -> [f64; 10] {
        [9.0, 3.0, 3.0, 4.0, 5.0, 4.9, 8.0, 3.3, 2.0, 0.1]
    }

    #[fixture]
    fn other_numbers() -> [f64; 10] {
        [14.2, 16.5, 3.4, 7.8, 18.4, 10.5, 7.4, 2.9, 22.3, 16.0]
    }

    #[rstest]
    fn test_percentile(numbers: [f64; 10]) {
        assert_approx_eq!(percentile(&numbers, 70), 4.95, 1e-2);
        assert_approx_eq!(percentile(&numbers, 50), 3.65, 1e-2);
    }

    #[ignore = "there is some divergence with the linear interpolation in numpy's percentile method"]
    #[rstest]
    fn test_percentile_with_other_cases(other_numbers: [f64; 10]) {
        assert_approx_eq!(percentile(&other_numbers, 80), 16.88, 1e-2);
    }
}
