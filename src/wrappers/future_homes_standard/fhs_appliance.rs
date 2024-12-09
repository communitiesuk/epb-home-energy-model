use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, Normal, Poisson};
use rand_pcg::Pcg64;

use crate::{
    core::units::{DAYS_PER_YEAR, HOURS_PER_DAY, WATTS_PER_KILOWATT},
    input::ApplianceGainsDetailsEvent,
};

pub(super) struct FhsAppliance {
    pub(super) standby_w: f64,
    pub(super) gains_frac: f64,
    pub(super) event_list: Vec<ApplianceGainsDetailsEvent>,
    pub(super) flat_schedule: Vec<f64>,
}

impl FhsAppliance {
    pub(super) fn new(
        util_unit: f64,
        annual_use_per_unit: f64,
        op_kwh: f64,
        event_duration: f64,
        standby_w: f64,
        gains_frac: f64,
        flat_profile: &[f64],
        seed: Option<usize>,
        duration_std_dev: Option<f64>,
    ) -> anyhow::Result<Self> {
        let annual_expected_uses = util_unit * annual_use_per_unit;
        let seed = seed.unwrap_or(DEFAULT_SEED);
        let duration_std_dev = duration_std_dev.unwrap_or(DEFAULT_DURATION_STD_DEV);
        let annual_expected_demand = annual_expected_uses * op_kwh
            + standby_w
                * ((HOURS_PER_DAY * DAYS_PER_YEAR) as f64 - annual_expected_uses * event_duration)
                / WATTS_PER_KILOWATT as f64;
        let (event_list, flat_schedule, standby_w) = Self::build_sched(
            flat_profile,
            seed,
            annual_expected_uses,
            annual_expected_demand,
            op_kwh,
            standby_w,
            event_duration,
            duration_std_dev,
        )?;
        Ok(Self {
            standby_w,
            gains_frac,
            event_list,
            flat_schedule,
        })
    }

    fn build_sched(
        flat_profile: &[f64],
        seed: usize,
        annual_expected_uses: f64,
        annual_expected_demand: f64,
        op_kwh: f64,
        standby_w: f64,
        event_duration: f64,
        duration_std_dev: f64,
    ) -> anyhow::Result<(Vec<ApplianceGainsDetailsEvent>, Vec<f64>, f64)> {
        // upstream Python here constructs a seed sequence from consecutive numbers - instead, here we sum the series
        let mut appliance_rng = Pcg64::seed_from_u64(
            (0..(flat_profile.len() + annual_expected_uses.ceil() as usize))
                .map(|x| (x + seed) as u64)
                .sum::<u64>(),
        );
        let events = flat_profile
            .iter()
            .map(|x| {
                let lambda =
                    (x * annual_expected_uses / DAYS_PER_YEAR as f64).max(f64::MIN_POSITIVE); // ensure lambda > 0. for poisson distribution
                let poisson = Poisson::new(lambda)?;
                Ok(poisson.sample(&mut appliance_rng))
            })
            .collect::<anyhow::Result<Vec<_>>>()?;
        let num_events = events.iter().copied().map(|x| x as usize).sum::<usize>();

        let normal_distribution = Normal::new(0., duration_std_dev)?;
        let mut event_size_deviations = normal_distribution
            .sample_iter(&mut appliance_rng)
            .take(num_events)
            .collect::<Vec<_>>();

        for deviation in event_size_deviations.iter_mut() {
            if *deviation < -1. {
                *deviation = normal_distribution.sample(&mut appliance_rng).max(-1.);
            }
        }

        let norm_events = num_events as f64 + event_size_deviations.iter().sum::<f64>();
        // adjustment in overall event size for random variation
        let f_appliance = (norm_events * op_kwh
            + standby_w * ((HOURS_PER_DAY * DAYS_PER_YEAR) as f64 - norm_events * event_duration)
                / WATTS_PER_KILOWATT as f64)
            / annual_expected_demand;

        // TODO (from Python) - this could be modifying the duration instead as with HW
        let expected_demand_w_event =
            op_kwh * WATTS_PER_KILOWATT as f64 / event_duration / f_appliance;
        let standby_adjusted = standby_w / f_appliance;
        let standby_w = standby_adjusted;
        let mut eventlist: Vec<ApplianceGainsDetailsEvent> = vec![];
        let mut sched = vec![standby_adjusted; flat_profile.len()];
        let flat_profile_len = flat_profile.len();

        let mut event_count: usize = Default::default();
        for (step, num_events_in_step) in events.into_iter().enumerate() {
            let mut start_offset = appliance_rng.gen::<f64>();
            for e in 0..(num_events_in_step.floor() as usize) {
                let demand_w_event = expected_demand_w_event;
                let duration = event_duration * (1. + event_size_deviations[event_count]);
                event_count += 1;
                // step will depend on timestep of flatprofile, always hourly so no adjustment
                eventlist.push(ApplianceGainsDetailsEvent {
                    start: step as f64 + start_offset,
                    duration,
                    demand_w: demand_w_event,
                });

                // build the flattened profile for use with loadshifting
                let mut integralx: f64 = Default::default();
                while integralx < duration {
                    let segment = (start_offset.ceil() - start_offset).min(duration - integralx);
                    sched[(step + (start_offset + integralx).floor() as usize)
                        % flat_profile_len] += (demand_w_event - standby_adjusted) * segment;
                    integralx += segment;
                }
                start_offset += e as f64 * duration;
            }
        }

        Ok((eventlist, sched, standby_w))
    }
}

const DEFAULT_SEED: usize = 37;
const DEFAULT_DURATION_STD_DEV: f64 = 0.;
