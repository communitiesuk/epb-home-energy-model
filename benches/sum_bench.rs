use accurate::{
    sum::{Kahan, Klein, NaiveSum, Neumaier, OnlineExactSum, Sum2, Sum3},
    traits::*,
};
use divan::{black_box, Bencher};

const LEN: usize = 17_520;

fn main() {
    divan::main();
}

fn sample_data() -> Vec<f64> {
    (0..LEN)
        .map(|i| {
            let base = (i as f64) * 0.137;
            if i % 97 == 0 {
                base * 1e6
            } else {
                base.sin().abs() * 12.3456789
            }
        })
        .collect()
}

#[divan::bench(types = [
    NaiveSum<f64>,
    Kahan<f64>,
    Neumaier<f64>,
    Klein<f64>,
    Sum2<f64>,
    Sum3<f64>,
    OnlineExactSum<f64>,
])]
fn sum_17520<Acc>(bencher: Bencher)
where
    Acc: SumAccumulator<f64>,
{
    bencher.with_inputs(sample_data).bench_refs(|data| {
        black_box(data.as_slice())
            .iter()
            .copied()
            .sum_with_accumulator::<Acc>()
    });
}

#[divan::bench(consts = [2, 4, 8])]
fn sum2_simd_lanes<const LANES: usize>(bencher: Bencher) {
    bencher
        .with_inputs(sample_data)
        .bench_refs(|data| sum2_simd::<LANES>(black_box(data.as_slice())));
}

#[inline]
pub fn sum2_simd<const LANES: usize>(values: &[f64]) -> f64 {
    let mut lanes = [Sum2::<f64>::zero(); LANES];
    let mut chunks = values.chunks_exact(LANES);

    for chunk in &mut chunks {
        for i in 0..LANES {
            lanes[i] += chunk[i];
        }
    }

    for (i, &rem) in chunks.remainder().iter().enumerate() {
        lanes[i] += rem;
    }

    // Evaluate each lane to `f64` via `.sum()`, then combine with `Sum2`
    lanes
        .into_iter()
        .map(|lane| lane.sum())
        .sum_with_accumulator::<Sum2<f64>>()
}
