pub fn min_of_2<T: PartialOrd + Copy>(first: T, second: T) -> T {
    if first < second {
        first
    } else {
        second
    }
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
}
