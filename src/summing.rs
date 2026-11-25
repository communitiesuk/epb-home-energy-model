use accurate::{sum::Sum2, traits::*};
use std::borrow::Borrow;

/// Function for performing fast precise summing, especially for a large series of numbers.
pub fn precise_sum<I>(numbers: I) -> f64
where
    I: IntoIterator,
    I::Item: Borrow<f64>,
{
    sum2_simd::<8, I>(numbers)
}

#[inline]
pub fn sum2_simd<const LANES: usize, I>(values: I) -> f64
where
    I: IntoIterator,
    I::Item: Borrow<f64>,
{
    let mut lanes = [Sum2::<f64>::zero(); LANES];
    let mut iter = values.into_iter();

    'outer: loop {
        let mut chunk = [0.0; LANES];
        for (idx, slot) in chunk.iter_mut().enumerate() {
            match iter.next() {
                Some(v) => *slot = *v.borrow(),
                None => {
                    for i in 0..idx {
                        lanes[i] += chunk[i];
                    }
                    break 'outer;
                }
            }
        }

        for i in 0..LANES {
            lanes[i] += chunk[i];
        }
    }

    // Evaluate each lane to `f64` via `.sum()`, then combine with `Sum2`
    lanes
        .into_iter()
        .map(|lane| lane.sum())
        .sum_with_accumulator::<Sum2<f64>>()
}
