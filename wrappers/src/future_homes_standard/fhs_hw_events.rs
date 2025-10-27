use crate::future_homes_standard::future_homes_standard::HourlyHotWaterEvent;
use crate::future_homes_standard::input::InputForProcessing;
use anyhow::{anyhow, bail};
use csv::Reader;
use hem::core::water_heat_demand::misc::frac_hot_water;
use hem::input::WaterHeatingEventType;
use indexmap::IndexMap;
use parking_lot::Mutex;
use partial_application::partial;
use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, Poisson};
use rand_mt::Mt64;
use rand_pcg::Pcg64;
use serde::Deserialize;
use smartstring::alias::String;
use std::fmt::{Debug, Formatter};
use std::io::{BufReader, Cursor};
use std::iter::Iterator;
use std::rc::Rc;
use strum::EnumCount;
use strum_macros::EnumCount as EnumCountMacro;

// utility for applying the SAP10.2 monthly factors
const MONTH_HOUR_STARTS: [usize; 12] = [
    744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760,
];
// from SAP10.2 J5
const BEHAVIOURAL_HW_FACTOR_M: [f64; 12] = [
    1.035, 1.021, 1.007, 0.993, 0.979, 0.965, 0.965, 0.979, 0.993, 1.007, 1.021, 1.035,
];
// from SAP10.2 j2
const OTHER_HW_FACTOR_M: [f64; 13] = [
    1.10, 1.06, 1.02, 0.98, 0.94, 0.90, 0.90, 0.94, 0.98, 1.02, 1.06, 1.10, 1.00,
];
const STANDARD_FILL: f64 = 73.;
pub(super) const STANDARD_BATH_SIZE: f64 = 180.;

fn bath_size_displaced(n_occupants: f64, bath_size: f64) -> anyhow::Result<f64> {
    // number of adults and children derived from Metabolic gains BSA calc
    let n_adults =
        (2.0001 * n_occupants.powf(0.8492) - 1.07451 * n_occupants) / (1.888074 - 1.07451);
    let n_children = n_occupants - n_adults;
    // average occupant weight, same assumptions as metabolic gain
    let w_kg = (78.6 * n_adults + 33.01 * n_children) / n_occupants;
    // assume density equal to water (in reality varies)
    let occupant_displacement_l = w_kg * 1.;
    let fill_vol_l = (occupant_displacement_l + STANDARD_FILL) * (bath_size / STANDARD_BATH_SIZE)
        - occupant_displacement_l;
    if fill_vol_l <= 0. {
        bail!("Bath too small, fill_vol_l: {}", fill_vol_l)
    }
    Ok(fill_vol_l)
}

#[derive(Debug)]
pub struct DrawoffGenerator {
    showers: Vec<Drawoff>,
    baths: Vec<Drawoff>,
    other: Vec<Drawoff>,
    which_shower: isize,
    which_bath: isize,
    which_other: isize,
}

pub fn reset_events_and_provide_drawoff_generator(
    n_occupants: f64,
    input: &mut InputForProcessing,
    fhw: f64,
    event_temperature: f64,
    hw_temperature: f64,
    cold_water_feed_temps: &[f64],
    part_g_bonus: f64,
) -> anyhow::Result<DrawoffGenerator> {
    let mut showers: Vec<Drawoff> = Default::default();
    let mut baths: Vec<Drawoff> = Default::default();
    let mut other: Vec<Drawoff> = Default::default();

    let shower_duration_func = move |event: DrawoffEvent| -> f64 {
        let shower_month_index = MONTH_HOUR_STARTS
            .iter()
            .position(|&value| value as f64 > event.time)
            .unwrap();
        event.duration * fhw * BEHAVIOURAL_HW_FACTOR_M[shower_month_index]
    };

    let bath_duration_func = move |bath_size: f64, flowrate: f64, event: DrawoffEvent| -> f64 {
        let bath_month_index = MONTH_HOUR_STARTS
            .iter()
            .position(|&value| value as f64 > event.time)
            .unwrap();
        let vol = bath_size * fhw * BEHAVIOURAL_HW_FACTOR_M[bath_month_index];
        // bath size is already a volume of warm water (not hot water)
        // so application frac_hw is unnecessary here
        vol / flowrate
    };

    let other_duration_func_gen = || {
        let cold_water_feed_temps: Vec<f64> = cold_water_feed_temps.to_vec();
        move |flow_rate: f64, event: DrawoffEvent| -> f64 {
            let other_month_index = MONTH_HOUR_STARTS
                .iter()
                .position(|&value| value as f64 > event.time)
                .unwrap();
            let frac_hw = frac_hot_water(
                event_temperature,
                hw_temperature,
                cold_water_feed_temps[event.time.floor() as usize],
            );
            (event.volume / frac_hw / flow_rate)
                * fhw
                * OTHER_HW_FACTOR_M[other_month_index]
                * part_g_bonus
        }
    };

    // set up events section of input
    // check if showers/baths are present
    // if multiple showers/baths are present, we need to cycle through them
    // if either is missing replace with the one that is present,
    // if neither is present, "other" events with same consumption as a bath should be used
    input.reset_water_heating_events()?;

    for shower in input.shower_keys()? {
        showers.push(Drawoff {
            event_type: "Shower".into(),
            name: shower,
            duration_fn: Rc::new(Mutex::new(Box::new(shower_duration_func))),
        });
    }

    for bath in input.bath_keys()? {
        let bath_size = input.size_for_bath_field(bath.as_str())?.ok_or(anyhow!(
            "Tried to access a bath input with a nonexistent key '{bath}'"
        ))?;
        // displacement of average occupant subtracted from volume of bath tub to work out fill volume
        let bath_size = bath_size_displaced(n_occupants, bath_size)?;
        let bath_flowrate = input
            .flowrate_for_bath_field(bath.as_str())?
            .ok_or(anyhow!(
                "Tried to access bath input with a nonexistent key '{bath}'"
            ))?;
        baths.push(Drawoff {
            event_type: "Bath".into(),
            duration_fn: Rc::new(Mutex::new(Box::new(partial!(move bath_duration_func =>
                bath_size,
                bath_flowrate,
                _
            )))),
            name: bath,
        });
    }

    for other_name in input.other_water_use_keys()? {
        let other_flow_rate = input
            .flow_rate_for_other_water_use_field(other_name.as_str())?
            .ok_or(anyhow!(
                "Tried to access an input for other water use with a nonexistent key '{other_name}'"
            ))?;
        let other_duration_func = other_duration_func_gen();
        other.push(Drawoff {
            event_type: "Other".into(),
            duration_fn: Rc::new(Mutex::new(Box::new(partial!(move other_duration_func =>
                other_flow_rate,
                _
            )))),
            name: other_name,
        })
    }

    // if there are no other events we need to add them
    // using a default of 8.0l/min flowrate -
    // event duration is calculated such as to deliver a fixed volume of water (otherdurationfunc above)
    // so this choice only affects how sharp peaks in HW demand can be.
    if other.is_empty() {
        let feed_type = if input.cold_water_source_has_header_tank()? {
            "header tank"
        } else {
            "mains water"
        };

        input.set_other_water_use_details(feed_type, 8.0)?;
        let other_flow_rate = input
            .flow_rate_for_other_water_use_field("other")?
            .unwrap_or_else(|| {
                panic!(
                    "Tried to access an input for other water use with a nonexistent key 'other'"
                )
            });
        let other_duration_func = other_duration_func_gen();
        other.push(Drawoff {
            event_type: "Other".into(),
            name: "other".into(),
            duration_fn: Rc::new(Mutex::new(Box::new(partial!(move other_duration_func =>
                other_flow_rate,
                _
            )))),
        })
    }

    // if no shower present, baths should be taken and vice versa.
    // If neither is present then bath sized drawoff
    match (showers.is_empty(), baths.is_empty()) {
        (true, false) => {
            showers.clone_from(&baths);
        }
        (false, true) => {
            baths.clone_from(&showers);
        }
        (true, true) => {
            // bath sized events occur whenever a shower or bath would
            // if there are no shower or bath facilities in the dwelling
            // using a default of 180L tub and 8.0l/min flowrate
            let bath_size = bath_size_displaced(n_occupants, STANDARD_BATH_SIZE)?;
            baths.push(Drawoff {
                event_type: "Other".into(),
                name: "other".into(),
                duration_fn: Rc::new(Mutex::new(Box::new(
                    partial!(move bath_duration_func => bath_size, 8.0, _),
                ))),
            });
            showers.push(Drawoff {
                event_type: "Other".into(),
                name: "other".into(),
                duration_fn: Rc::new(Mutex::new(Box::new(
                    partial!(move bath_duration_func => bath_size, 8.0, _),
                ))),
            });
        }
        _ => {}
    }

    Ok(DrawoffGenerator {
        showers,
        baths,
        other,
        which_shower: -1,
        which_bath: -1,
        which_other: -1,
    })
}

impl DrawoffGenerator {
    pub fn get_shower(&mut self) -> &mut Drawoff {
        self.which_shower = ((self.which_shower + 1) as usize % self.showers.len()) as isize;
        self.showers.get_mut(self.which_shower as usize).unwrap()
    }

    pub fn get_bath(&mut self) -> &mut Drawoff {
        self.which_bath = ((self.which_bath + 1) as usize % self.baths.len()) as isize;
        self.baths.get_mut(self.which_bath as usize).unwrap()
    }

    pub fn get_other(&mut self) -> &mut Drawoff {
        self.which_other = ((self.which_other + 1) as usize % self.other.len()) as isize;
        self.other.get_mut(self.which_other as usize).unwrap()
    }
}

type DurationFn = Rc<Mutex<Box<dyn FnMut(DrawoffEvent) -> f64>>>;

#[derive(Clone)]
pub struct Drawoff {
    pub event_type: String,
    pub name: String,
    pub duration_fn: DurationFn,
}

impl Drawoff {
    pub fn call_duration_fn(&mut self, event: HourEvent) -> f64 {
        self.duration_fn.lock()(event.into())
    }
}

impl Debug for Drawoff {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Drawoff")
            .field("event_type", &self.event_type)
            .field("name", &self.name)
            .field("duration_fn", &"duration_fn content hidden")
            .finish_non_exhaustive()
    }
}

pub struct DrawoffEvent {
    time: f64,
    _event_type: WaterHeatingEventType,
    volume: f64,
    duration: f64,
}

impl From<HourEvent> for DrawoffEvent {
    fn from(value: HourEvent) -> Self {
        Self {
            time: value.time,
            _event_type: value.event_type.into(),
            volume: value.volume,
            duration: value.duration,
        }
    }
}

impl From<SimpleLabelBasedOn900KSample> for WaterHeatingEventType {
    fn from(value: SimpleLabelBasedOn900KSample) -> Self {
        if value.is_shower_type() {
            Self::Shower
        } else if value.is_bath_type() {
            Self::Bath
        } else {
            Self::Other
        }
    }
}

#[derive(Debug)]
pub struct HotWaterEventGenerator {
    week: IndexMap<DayOfWeek, IndexMap<SimpleLabelBasedOn900KSample, DayDigest>>,
    rng: Mt64,
}

impl HotWaterEventGenerator {
    pub fn new(
        daily_dhw_vol: f64,
        hw_seed: Option<u64>,
        banding: Option<BandingOrigin>,
    ) -> anyhow::Result<Self> {
        let banding = banding.unwrap_or(BandingOrigin::Correct);
        let mut week: IndexMap<DayOfWeek, IndexMap<SimpleLabelBasedOn900KSample, DayDigest>> =
            IndexMap::from([
                (DayOfWeek::Monday, Default::default()),
                (DayOfWeek::Tuesday, Default::default()),
                (DayOfWeek::Wednesday, Default::default()),
                (DayOfWeek::Thursday, Default::default()),
                (DayOfWeek::Friday, Default::default()),
                (DayOfWeek::Saturday, Default::default()),
                (DayOfWeek::Sunday, Default::default()),
            ]);

        let rng = SeedableRng::seed_from_u64(hw_seed.unwrap_or(RNG_SEED));
        let mut rng_poisson = Pcg64::seed_from_u64(hw_seed.unwrap_or(RNG_SEED));

        let mut decile: i8 = -1;
        let mut banding_correction = 1.0;

        let _target_dhw_vol = daily_dhw_vol;

        let mut bands_file_data = Vec::with_capacity(DECILE_BANDING_FILE.lines().count() - 1);
        {
            let mut decile_bands_reader =
                Reader::from_reader(BufReader::new(Cursor::new(DECILE_BANDING_FILE)));
            for band in decile_bands_reader.deserialize() {
                let band: DecileBand =
                    band.expect("There was a problem reading the static decile banding file.");
                bands_file_data.push(band);
                if daily_dhw_vol >= band.min_daily_dhw_vol && daily_dhw_vol < band.max_daily_dhw_vol
                {
                    decile = band.decile - 1;
                    banding_correction = daily_dhw_vol / band.calibration_daily_dhw_vol;
                }
            }
            if decile == -1 {
                if daily_dhw_vol < bands_file_data[0].min_daily_dhw_vol {
                    decile = 0;
                    banding_correction = daily_dhw_vol
                        / bands_file_data
                            .first()
                            .expect("No decile bands were read from the file")
                            .calibration_daily_dhw_vol;
                } else if daily_dhw_vol > bands_file_data[9].min_daily_dhw_vol {
                    decile = 9;
                    banding_correction = daily_dhw_vol
                        / bands_file_data
                            .last()
                            .expect("No decile bands were read from the file")
                            .calibration_daily_dhw_vol;
                }
            }
            if decile == -1 {
                bail!("HW decile error, exiting");
            }
        }

        if let BandingOrigin::NotCorrect = banding {
            banding_correction = 1.0;
        }

        {
            let mut decile_events_reader =
                Reader::from_reader(BufReader::new(Cursor::new(DECILE_EVENTS_FILE)));
            for event in decile_events_reader.deserialize() {
                let event: DecileEvent =
                    event.expect("There was a problem reading the static decile events file");
                if event.decile - 1 == decile {
                    week.entry(event.day_name).or_default().insert(
                        event.simple_labels2_based_on_900k_sample,
                        DayDigest {
                            event_count: event.event_count,
                            median_event_volume: event.median_event_volume,
                            mean_event_volume: event.mean_event_volume,
                            median_duration: event.median_duration / 60.,
                            mean_duration: event.mean_duration / 60.,
                            hourly_event_counts: [0usize; 24],
                            hourly_event_distribution: Default::default(),
                        },
                    );
                }
            }
        }

        {
            let mut decile_event_times_reader =
                Reader::from_reader(BufReader::new(Cursor::new(DECILE_EVENT_TIMES_FILE)));
            for times in decile_event_times_reader.deserialize() {
                let times: DecileEventTimes =
                    times.expect("There was a problem reading the static decile event times file");
                let current_day_digest = week
                    .get_mut(&times.day_name)
                    .unwrap()
                    .get_mut(&times.simple_labels2_based_on_900k_sample)
                    .unwrap();
                current_day_digest.hourly_event_counts[times.hour] = times.event_count;
            }
        }

        let week_keys = week.keys().cloned().collect::<Vec<_>>();
        for day in week_keys {
            let week_day_keys = week[&day].keys().cloned().collect::<Vec<_>>();
            for event_type in week_day_keys {
                let current_digest = week.get_mut(&day).unwrap().get_mut(&event_type).unwrap();
                let hourly_event_counts = current_digest.hourly_event_counts;
                let sum_event_count = hourly_event_counts.iter().sum::<usize>();

                // successive calls to a poisson distribution with fixed seed
                // will yield the same answer - have to ask rng to generate an
                // array of poisson samples and draw from it.
                // generate array of size 53 as each hour is unique per week of the year
                current_digest.hourly_event_distribution =
                    Some(hourly_event_counts.map(|x| PoissonDistribution {
                        poisson_arr: {
                            let lambda = banding_correction * x as f64 * current_digest.event_count
                                / sum_event_count as f64;
                            let poisson = Poisson::new(lambda).unwrap_or_else(|_| {
                                panic!(
                                    "Unable to create a poisson generator with the lambda {lambda}"
                                )
                            });
                            (0..POISSON_DISTRIBUTION_SIZE)
                                .map(|_| poisson.sample(&mut rng_poisson))
                                .collect::<Vec<f64>>()
                        },
                        poisson_arr_idx: 0,
                    }));
            }
        }

        Ok(Self { week, rng })
    }

    fn events_in_hour(
        &mut self,
        time: usize,
        event_type: SimpleLabelBasedOn900KSample,
        digest: &mut DayDigest,
    ) -> Vec<HourEvent> {
        let mut out: Vec<HourEvent> = Default::default();
        let hourly_event_distribution = digest.hourly_event_distribution.as_mut().expect(
            "Found a DayDigest value that did not have its hourly event distribution set on it",
        );
        let count = hourly_event_distribution[time % 24].poisson_arr
            [hourly_event_distribution[time % 24].poisson_arr_idx];
        hourly_event_distribution[time % 24].poisson_arr_idx += 1;
        for _ in 0..(count as usize) {
            out.push(HourEvent {
                time: time as f64 + self.rng.random::<f64>(),
                event_type,
                volume: digest.mean_event_volume,
                duration: digest.mean_duration,
            });
        }

        out
    }

    pub fn overlap_check(
        &mut self,
        hourly_events: &mut Vec<Vec<HourlyHotWaterEvent>>,
        matching_types: &[WaterHeatingEventType],
        event_start: f64,
        duration: f64,
    ) {
        let mut event_start = event_start;
        let hourly_events_idx = event_start.floor() as usize;
        let event_start_events: Vec<HourlyHotWaterEvent> =
            hourly_events[hourly_events_idx].to_vec();
        for existing_event in event_start_events {
            if matching_types.contains(&existing_event.event_type)
                && (event_start >= existing_event.start && event_start < existing_event.end)
                || (event_start + duration / 60. >= existing_event.start
                    && event_start + duration / 60. < existing_event.end)
            {
                // events are overlapping, and we need to reroll the time until they aren't
                event_start = self.reroll_event_time(event_start);
                self.overlap_check(hourly_events, matching_types, event_start, duration);
            }
        }
    }

    /// Sometimes events will overlap, and we need to change the time so they don't
    /// do this by adding random value between 0-30 mins to current time
    /// until it does not overlap with anything
    fn reroll_event_time(&mut self, time: f64) -> f64 {
        (time + self.rng.random::<f64>() / 2.) % 8760.
    }

    pub fn build_annual_hw_events(&mut self, start_day: usize) -> anyhow::Result<Vec<HourEvent>> {
        let mut list_days: Vec<IndexMap<SimpleLabelBasedOn900KSample, DayDigest>> =
            self.week.values().cloned().collect::<Vec<_>>();
        let mut annual_hw_events: Vec<HourEvent> =
            Vec::with_capacity(8760 * SimpleLabelBasedOn900KSample::COUNT);
        for day in 0..365 {
            for hour in 0..24 {
                let day_keys = list_days[(day + start_day) % 7]
                    .keys()
                    .cloned()
                    .collect::<Vec<_>>();
                for event_type in day_keys {
                    let mut events_to_append = self.events_in_hour(
                        hour + (day * 24),
                        event_type,
                        list_days[day % 7].get_mut(&event_type).unwrap(),
                    );
                    annual_hw_events.append(&mut events_to_append);
                }
            }
        }

        Ok(annual_hw_events)
    }
}

const RNG_SEED: u64 = 37;

#[derive(Clone, Debug, Deserialize)]
struct DayDigest {
    event_count: f64,
    #[allow(dead_code)]
    median_event_volume: f64,
    mean_event_volume: f64,
    #[allow(dead_code)]
    #[serde(rename = "median_dur")]
    median_duration: f64,
    #[serde(rename = "mean_dur")]
    mean_duration: f64,
    hourly_event_counts: [usize; 24],
    hourly_event_distribution: Option<[PoissonDistribution; 24]>,
}

#[derive(Clone, Debug, Default, Deserialize)]
struct PoissonDistribution {
    poisson_arr: Vec<f64>,
    poisson_arr_idx: usize,
}

const POISSON_DISTRIBUTION_SIZE: usize = 53;

#[derive(Deserialize)]
struct DecileEvent {
    event_count: f64,
    decile: i8,
    day_name: DayOfWeek,
    #[allow(dead_code)]
    day_wk_mon1: usize,
    simple_labels2_based_on_900k_sample: SimpleLabelBasedOn900KSample,
    median_event_volume: f64,
    mean_event_volume: f64,
    #[serde(rename = "median_dur")]
    median_duration: f64,
    #[serde(rename = "mean_dur")]
    mean_duration: f64,
}

#[derive(Deserialize)]
struct DecileEventTimes {
    hour: usize,
    event_count: usize,
    #[allow(dead_code)]
    no_hholds: Option<usize>,
    day_name: DayOfWeek,
    #[allow(dead_code)]
    day_wk_mon1: Option<usize>,
    simple_labels2_based_on_900k_sample: SimpleLabelBasedOn900KSample,
}

#[derive(Clone, Copy, Debug)]
pub struct HourEvent {
    pub time: f64,
    pub event_type: SimpleLabelBasedOn900KSample,
    pub volume: f64,
    duration: f64,
}

pub enum BandingOrigin {
    Correct,
    #[allow(dead_code)]
    NotCorrect,
}

#[derive(Clone, Copy, Deserialize, Eq, Hash, PartialEq, Debug)]
enum DayOfWeek {
    Monday,
    Tuesday,
    Wednesday,
    Thursday,
    Friday,
    Saturday,
    Sunday,
}

#[derive(Clone, Copy, Debug, Deserialize, EnumCountMacro, Eq, Hash, PartialEq)]
pub enum SimpleLabelBasedOn900KSample {
    #[serde(rename = "0_small_tap")]
    SmallTap,
    #[serde(rename = "1_long_tap")]
    LongTap,
    #[serde(rename = "2_shower")]
    Shower,
    #[serde(rename = "3_bath")]
    Bath,
    #[serde(rename = "2b_shower_big")]
    ShowerBig,
}

impl SimpleLabelBasedOn900KSample {
    pub fn is_shower_type(&self) -> bool {
        matches!(
            self,
            SimpleLabelBasedOn900KSample::Shower | SimpleLabelBasedOn900KSample::ShowerBig
        )
    }

    pub fn is_bath_type(&self) -> bool {
        matches!(self, Self::Bath)
    }
}

// bundle the data CSVs into the binary
static DECILE_BANDING_FILE: &str = include_str!("decile_banding.csv");
static DECILE_EVENTS_FILE: &str = include_str!("day_of_week_events_by_decile.csv");
static DECILE_EVENT_TIMES_FILE: &str = include_str!("day_of_week_events_by_decile_event_times.csv");

// represents a record from the decile bands file
#[derive(Clone, Copy, Debug, Deserialize)]
struct DecileBand {
    #[allow(dead_code)]
    no_hholds: usize,
    decile: i8,
    #[allow(dead_code)]
    median_daily_dhw_vol: f64,
    min_daily_dhw_vol: f64,
    max_daily_dhw_vol: f64,
    calibration_daily_dhw_vol: f64,
    #[allow(dead_code)]
    #[serde(rename = "calibration_DHW_variance")]
    calibration_dhw_variance: f64,
}
