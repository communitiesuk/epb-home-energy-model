pub fn min_of_2<T: PartialOrd + Copy>(first: T, second: T) -> T {
    if first < second {
        first
    } else {
        second
    }
}

pub(crate) fn min_of_3<T: PartialOrd + Copy>(first: T, second: T, third: T) -> T {
    let first_min = min_of_2(first, second);
    min_of_2(first_min, third)
}

pub fn max_of_2<T: PartialOrd + Copy>(first: T, second: T) -> T {
    if first > second {
        first
    } else {
        second
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[rstest]
    pub fn should_calc_2_as_min_of_2_and_4_ints() {
        assert_eq!(min_of_2(2, 4), 2);
    }

    #[rstest]
    pub fn should_calc_2_as_min_of_4_and_2_ints() {
        assert_eq!(min_of_2(4, 2), 2);
    }

    #[rstest]
    pub fn should_calc_2_as_min_of_4_and_2_floats() {
        assert_eq!(min_of_2(4., 2.), 2.);
    }

    #[rstest]
    pub fn should_calc_4_as_max_of_2_and_4_ints() {
        assert_eq!(max_of_2(2, 4), 4);
    }

    #[rstest]
    pub fn should_calc_4_as_max_of_4_and_2_ints() {
        assert_eq!(max_of_2(4, 2), 4);
    }

    #[rstest]
    pub fn should_calc_4_as_max_of_4_and_2_floats() {
        assert_eq!(max_of_2(4., 2.), 4.);
    }

    #[rstest]
    pub fn should_calc_a_as_min_of_b_and_a_chars() {
        assert_eq!(min_of_2('b', 'a'), 'a');
    }

    #[test]
    fn should_calc_2_as_min_of_2_3_4() {
        assert_eq!(min_of_3(2., 3., 4.), 2.);
    }

    #[test]
    fn should_calc_2_as_min_of_3_4_2() {
        assert_eq!(min_of_3(3., 4., 2.), 2.);
    }

    #[test]
    fn should_calc_0_as_min_of_1_0_2() {
        assert_eq!(min_of_3(1., 0., 2.), 0.);
    }

    #[test]
    fn should_calc_minus_10_as_min_of_three_numbers() {
        assert_eq!(min_of_3(33., -10., 1.55), -10.);
    }
}
