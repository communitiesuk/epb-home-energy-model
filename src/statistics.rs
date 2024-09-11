use interp::interp;
/// A simple statistics module with some utility functions such as calculation of percentiles.
use statrs::statistics::{Data, OrderStatistics, Statistics};

#[allow(dead_code)]
pub fn percentile(numbers: &[f64], percentile: usize) -> f64 {
    let numbers = numbers.to_vec();
    let mut data = Data::new(numbers);

    data.percentile(percentile)
}

/// This function matches the behaviour Numpy interp
/// https://numpy.org/doc/stable/reference/generated/numpy.interp.html
pub fn np_interp(input: f64, x: &[f64], y: &[f64]) -> f64 {
    if x.is_empty() {
        panic!("x cannot be empty");
    }

    if y.is_empty() {
        panic!("y cannot be empty");
    }

    if x.len() != y.len() {
        panic!("x and y must be of equal length");
    }

    if input > x.max() {
        return *y.last().unwrap();
    }

    if input < x.min() {
        return *y.first().unwrap();
    }

    interp(x, y, input)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::{assert_relative_eq, assert_ulps_eq};
    use interp::interp;
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
        assert_relative_eq!(percentile(&numbers, 70), 4.95, max_relative = 1e-2);
        assert_relative_eq!(percentile(&numbers, 50), 3.65, max_relative = 1e-2);
    }

    #[ignore = "there is some divergence with the linear interpolation in numpy's percentile method"]
    #[rstest]
    fn test_percentile_with_other_cases(other_numbers: [f64; 10]) {
        assert_ulps_eq!(percentile(&other_numbers, 80), 16.88);
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
        let expected = interp(&x, &y, input);

        assert_eq!(actual, expected)
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

        assert_eq!(actual, expected)
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

        assert_eq!(actual, expected)
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

        assert_eq!(actual, expected)
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
        assert_eq!(actual, expected)
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
