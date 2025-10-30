#![allow(non_snake_case)]

use crate::compare_floats::{max_of_2, min_of_2};
use crate::core::space_heat_demand::building_element::sky_view_factor;
use crate::core::units::HOURS_PER_DAY;
use crate::input::{deserialize_orientation, serialize_orientation, ExternalConditionsInput};
use crate::simulation_time::{SimulationTimeIteration, SimulationTimeIterator, HOURS_IN_DAY};
use anyhow::{anyhow, bail};
#[cfg(test)]
use approx::{AbsDiffEq, RelativeEq};
use itertools::Itertools;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
pub enum DaylightSavingsConfig {
    #[serde(rename = "applicable and taken into account")]
    ApplicableAndTakenIntoAccount,
    #[serde(rename = "applicable but not taken into account")]
    ApplicableButNotTakenIntoAccount,
    #[serde(rename = "not applicable")]
    NotApplicable,
}

#[derive(Clone, Debug, Default, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
// #[serde(deny_unknown_fields)] // TODO: restore this for all versions after 0.36
pub struct ShadingSegment {
    #[serde(rename = "start360")]
    #[serde(
        deserialize_with = "deserialize_orientation",
        serialize_with = "serialize_orientation"
    )]
    pub(crate) start: f64,
    #[serde(rename = "end360")]
    #[serde(
        deserialize_with = "deserialize_orientation",
        serialize_with = "serialize_orientation"
    )]
    pub(crate) end: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(rename = "shading")]
    pub(crate) shading_objects: Option<Vec<ShadingObject>>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[cfg_attr(test, derive(PartialEq))]
#[serde(deny_unknown_fields)]
pub(crate) struct ShadingObject {
    #[serde(rename = "type")]
    pub(crate) object_type: ShadingObjectType,
    pub(crate) height: f64,
    pub(crate) distance: f64,
}

#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(deny_unknown_fields, rename_all = "lowercase", tag = "type")]
pub enum WindowShadingObject {
    Obstacle {
        height: f64,
        distance: f64,
        transparency: f64,
    },
    Overhang {
        depth: f64,
        distance: f64,
    },
    SideFinRight {
        depth: f64,
        distance: f64,
    },
    SideFinLeft {
        depth: f64,
        distance: f64,
    },
    Reveal {
        depth: f64,
        distance: f64,
    },
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "arbitrary", derive(arbitrary::Arbitrary))]
#[serde(rename_all = "lowercase")]
pub(crate) enum ShadingObjectType {
    Obstacle,
    Overhang,
}

#[derive(Clone, Debug)]
pub struct ExternalConditions {
    air_temps: Vec<f64>,
    wind_speeds: Vec<f64>,
    wind_directions: Vec<f64>,
    diffuse_horizontal_radiations: Vec<f64>,
    direct_beam_radiations: Vec<f64>,
    solar_reflectivity_of_ground: Vec<f64>,
    pub(crate) latitude: f64,
    #[allow(dead_code)]
    pub(crate) longitude: f64,
    #[allow(dead_code)]
    pub(crate) timezone: i32,
    pub(crate) start_day: u32,
    time_series_step: f64,
    shading_segments: Option<Vec<ShadingSegment>>,
    solar_declinations: Vec<f64>,
    solar_hour_angles: Vec<f64>,
    solar_altitudes: Vec<f64>,
    solar_zenith_angles: Vec<f64>,
    solar_azimuth_angles: Vec<f64>,
    f1_circumsolar_brightness_coefficients: Vec<f64>,
    f2_horizontal_brightness_coefficients: Vec<f64>,
}

/// Arguments:
/// * `simulation_time` - reference to SimulationTime iteration
/// * `air_temps` - list of external air temperatures, in deg C (one entry per hour)
/// * `wind_speeds` - list of wind speeds, in m/s (one entry per hour)
/// * `wind_directions` - list of wind directions in degrees where North=0, East=90,
///                       South=180, West=270. Values range: 0 to 360.
///                       Wind direction is reported by the direction from which it originates.
///                       E.g, a southerly (180 degree) wind blows from the south to the north.
/// * `diffuse_horizontal_radiation` - list of diffuse horizontal radiation values, in W/m2 (one entry per hour)
/// * `direct_beam_radiation` - list of direct beam radiation values, in W/m2 (one entry per hour)
/// * `solar_reflectivity_of_ground` - list of ground reflectivity values, 0 to 1 (one entry per hour)
/// * `latitude` - latitude of weather station, angle from south, in degrees (single value)
/// * `longitude` - longitude of weather station, easterly +ve westerly -ve, in degrees (single value)
/// * `timezone` - timezone of weather station, -12 to 12 (single value)
/// * `start_day` - first day of the time series, day of the year, 0 to 365 (single value)
/// * `end_day` - last day of the time series, day of the year, 0 to 365 (single value)
/// * `time_series_step` - timestep of the time series data, in hours
/// * `january_first` - day of the week for January 1st, monday to sunday, 1 to 7 (single value)
/// * `daylight_savings` - handling of daylight savings time, (single value)
///                        e.g. applicable and taken into account,
///                        applicable but not taken into account,
///                        not applicable
/// * `leap_day_included` - whether climate data includes a leap day, true or false (single value)
/// * `direct_beam_conversion_needed` - A flag to indicate whether direct beam radiation from climate data needs to be
///                                     converted from horizontal to normal incidence. If normal direct beam radiation
///                                     values are provided then no conversion is needed.
/// * `shading_segments` - data splitting the ground plane into segments (8-36) and giving height
///                        and distance to shading objects surrounding the building
impl ExternalConditions {
    pub fn new(
        simulation_time: &SimulationTimeIterator,
        air_temps: Vec<f64>,
        wind_speeds: Vec<f64>,
        wind_directions: Vec<f64>,
        diffuse_horizontal_radiations: Vec<f64>,
        direct_beam_radiations: Vec<f64>,
        solar_reflectivity_of_ground: Vec<f64>,
        latitude: f64,
        longitude: f64,
        timezone: i32,
        start_day: u32,
        _end_day: Option<u32>,
        time_series_step: f64,
        _january_first: Option<u32>,
        _daylight_savings: Option<DaylightSavingsConfig>,
        leap_day_included: bool,
        direct_beam_conversion_needed: bool,
        shading_segments: Option<Vec<ShadingSegment>>,
    ) -> Self {
        let days_in_year = if leap_day_included { 366 } else { 365 };
        let hours_in_year = days_in_year * HOURS_IN_DAY;
        let time_shift = init_time_shift(timezone, longitude);

        // # Calculate earth orbit deviation for each day of year
        let earth_orbit_deviations = (0..days_in_year)
            .map(init_earth_orbit_deviation)
            .collect::<Vec<f64>>();

        let extra_terrestrial_radiation = (0..days_in_year)
            .map(|day| init_extra_terrestrial_radiation(earth_orbit_deviations[day as usize]))
            .collect::<Vec<f64>>();

        let solar_declinations = (0..days_in_year)
            .map(|day| init_solar_declination(earth_orbit_deviations[day as usize]))
            .collect::<Vec<f64>>();

        let equations_of_time = (0..days_in_year)
            .map(init_equation_of_time)
            .collect::<Vec<f64>>();

        let solar_times = (0..hours_in_year)
            .map(|hour| {
                init_solar_time(
                    hour % 24,
                    equations_of_time[hour.div_euclid(24) as usize],
                    time_shift,
                    hour,
                )
            })
            .collect::<Vec<f64>>();

        let solar_hour_angles = (0..hours_in_year)
            .map(|hour| init_solar_hour_angle(solar_times[hour as usize]))
            .collect::<Vec<f64>>();

        let solar_altitudes = (0..hours_in_year)
            .map(|hour| {
                init_solar_altitude(
                    latitude,
                    solar_declinations[hour.div_euclid(24) as usize],
                    solar_hour_angles[hour as usize],
                )
            })
            .collect::<Vec<f64>>();

        let solar_zenith_angles = (0..hours_in_year)
            .map(|hour| init_solar_zenith_angle(solar_altitudes[hour as usize]))
            .collect::<Vec<f64>>();

        let solar_azimuth_angles = (0..hours_in_year)
            .map(|hour| {
                init_solar_azimuth_angle(
                    latitude,
                    solar_declinations[hour.div_euclid(24) as usize],
                    solar_hour_angles[hour as usize],
                    solar_altitudes[hour as usize],
                )
            })
            .collect::<Vec<f64>>();

        let air_masses = (0..hours_in_year)
            .map(|hour| init_air_mass(solar_altitudes[hour as usize]))
            .collect::<Vec<f64>>();

        let simtime = simulation_time.clone();
        let direct_beam_radiations = simtime
            .map(|it| {
                init_direct_beam_radiation(
                    direct_beam_conversion_needed,
                    direct_beam_radiations[it.time_series_idx(start_day, time_series_step)],
                    solar_altitudes[it.current_hour() as usize],
                )
            })
            .collect::<Vec<f64>>();

        let simtime = simulation_time.clone();
        let diffuse_horizontal_radiations = simtime
            .map(|it| {
                diffuse_horizontal_radiations[it.time_series_idx(start_day, time_series_step)]
            })
            .collect::<Vec<f64>>();

        let simtime = simulation_time.clone();
        let dimensionless_clearness_parameters = simtime
            .map(|it| {
                init_dimensionless_clearness_parameter(
                    diffuse_horizontal_radiations[it.index],
                    direct_beam_radiations[it.index],
                    solar_altitudes[it.current_hour() as usize],
                )
            })
            .collect::<Vec<f64>>();

        let simtime = simulation_time.clone();
        let dimensionless_sky_brightness_parameters = simtime
            .map(|it| {
                init_dimensionless_sky_brightness_parameter(
                    air_masses[it.current_hour() as usize],
                    diffuse_horizontal_radiations[it.index],
                    extra_terrestrial_radiation[it.current_day() as usize],
                )
            })
            .collect::<Vec<f64>>();

        let simtime = simulation_time.clone();
        // # Calculate circumsolar brightness coefficient, F1 for each timestep

        let f1_circumsolar_brightness_coefficients = simtime
            .map(|it| {
                init_f1_circumsolar_brightness_coefficient(
                    dimensionless_clearness_parameters[it.index],
                    dimensionless_sky_brightness_parameters[it.index],
                    solar_zenith_angles[it.current_hour() as usize],
                )
            })
            .collect::<Vec<f64>>();

        let simtime = simulation_time.clone();
        // # Calculate horizontal brightness coefficient, F2 for each timestep
        let f2_horizontal_brightness_coefficients = simtime
            .map(|it| {
                init_f2_horizontal_brightness_coefficient(
                    dimensionless_clearness_parameters[it.index],
                    dimensionless_sky_brightness_parameters[it.index],
                    solar_zenith_angles[it.current_hour() as usize],
                )
            })
            .collect::<Vec<f64>>();

        Self {
            air_temps,
            wind_speeds,
            wind_directions,
            diffuse_horizontal_radiations,
            direct_beam_radiations,
            solar_reflectivity_of_ground,
            latitude,
            longitude,
            timezone,
            start_day,
            time_series_step,
            shading_segments,
            solar_declinations,
            solar_hour_angles,
            solar_altitudes,
            solar_zenith_angles,
            solar_azimuth_angles,
            f1_circumsolar_brightness_coefficients,
            f2_horizontal_brightness_coefficients,
        }
    }

    pub(crate) fn air_temp(&self, simtime: &SimulationTimeIteration) -> f64 {
        self.air_temp_with_offset(simtime, 0)
    }

    // This method represents using the Python air_temp method with an offset parameter provided.
    pub(crate) fn air_temp_with_offset(
        &self,
        simtime: &SimulationTimeIteration,
        offset: usize,
    ) -> f64 {
        let mut idx = simtime.time_series_idx(self.start_day, self.time_series_step) + offset;
        if idx >= self.air_temps.len() {
            idx -= self.air_temps.len();
        }
        self.air_temps[idx]
    }

    fn _air_temp_for_timestep_idx(&self, timestep_idx: usize) -> f64 {
        self.air_temps[timestep_idx]
    }

    pub(crate) fn air_temp_annual(&self) -> Option<f64> {
        if self.air_temps.len() != 8760 {
            return None;
        }
        let sum: f64 = self.air_temps.iter().sum();
        Some(sum / self.air_temps.len() as f64)
    }

    pub(crate) fn air_temp_monthly(&self, current_month_start_end_hours: (u32, u32)) -> f64 {
        let (idx_start, idx_end) = current_month_start_end_hours;
        let (idx_start, idx_end) = (idx_start as usize, idx_end as usize);
        let air_temps_month = &self.air_temps[idx_start..idx_end];
        let sum: f64 = air_temps_month.iter().sum();
        sum / air_temps_month.len() as f64
    }

    pub(crate) fn air_temp_annual_daily_average_min(&self) -> f64 {
        // only works if data for a whole year has been provided
        debug_assert!(self.air_temps.len() == 8760);
        // determine the air temperatures for each day
        let no_of_days = self.air_temps.len() / HOURS_PER_DAY as usize;
        let daily_averages = (0..no_of_days)
            .map(|i| {
                self.air_temps[(i * HOURS_PER_DAY as usize)..((i + 1) * HOURS_PER_DAY as usize)]
                    .iter()
                    .sum::<f64>()
                    / HOURS_PER_DAY as f64
            })
            .collect_vec();
        daily_averages
            .into_iter()
            .min_by(|a, b| a.total_cmp(b))
            .unwrap()
    }

    pub fn wind_speed(&self, simulation_time: &SimulationTimeIteration) -> f64 {
        self.wind_speed_for_timestep_idx(
            simulation_time.time_series_idx(self.start_day, self.time_series_step),
        )
    }

    fn wind_speed_for_timestep_idx(&self, timestep_idx: usize) -> f64 {
        self.wind_speeds[timestep_idx]
    }

    pub fn wind_speed_annual(&self) -> Option<f64> {
        if self.wind_speeds.len() != (8760.0 / self.time_series_step) as usize {
            return None;
        }
        let sum: f64 = self.wind_speeds.iter().sum();
        Some(sum / self.wind_speeds.len() as f64)
    }

    pub fn wind_direction(&self, simulation_time: SimulationTimeIteration) -> f64 {
        self.wind_directions[simulation_time.time_series_idx(self.start_day, self.time_series_step)]
    }

    /// Return the average wind direction for the whole year
    pub fn wind_direction_annual(&self) -> f64 {
        // only works if data for whole year has been provided
        debug_assert!(self.wind_speeds.len() == 8760);
        debug_assert!(self.wind_directions.len() == 8760);
        let (x_total, y_total) = self
            .wind_speeds
            .iter()
            .zip(self.wind_directions.iter())
            .fold(
                (0., 0.),
                |(x_total, y_total), (wind_speed, wind_direction)| {
                    (
                        x_total + wind_speed * wind_direction.to_radians().cos(),
                        y_total + wind_speed * wind_direction.to_radians().sin(),
                    )
                },
            );
        // Take average of x and y for each timestep and then convert back to angle
        let x_average = x_total / self.wind_directions.len() as f64;
        let y_average = y_total / self.wind_directions.len() as f64;
        // NB. the atan2 implementation in Rust currently is consistent with the implementation in libc,
        // but this could potentially change in the future - the precision is marked as "unspecified".
        let wind_direction_average = y_average.atan2(x_average).to_degrees();

        wind_direction_average.rem_euclid(360.) // cannot use % operator here as we need lowest non-negative remainder, which % does not give us
    }

    pub fn diffuse_horizontal_radiation(&self, timestep_idx: usize) -> f64 {
        // self.diffuse_horizontal_radiations[self.simulation_time.current_index()]
        self.diffuse_horizontal_radiations[timestep_idx]
    }

    pub fn direct_beam_radiation(&self, timestep_idx: usize) -> f64 {
        // self.direct_beam_radiations[self.simulation_time.current_index()]
        self.direct_beam_radiations[timestep_idx]
    }

    /// Return clockwise 0/360 orientation angle from anti-clockwise -180/+180 basis
    fn orientation360(&self, orientation: f64) -> f64 {
        180. - orientation
    }

    pub fn solar_reflectivity_of_ground(&self, simulation_time: &SimulationTimeIteration) -> f64 {
        self.solar_reflectivity_of_ground
            [simulation_time.time_series_idx(self.start_day, self.time_series_step)]
    }

    fn solar_angle_of_incidence(
        &self,
        tilt: f64,
        orientation: f64,
        simulation_time: &SimulationTimeIteration,
    ) -> f64 {
        // """  calculates the solar angle of incidence, which is the angle of incidence of the
        // solar beam on an inclined surface and is determined as function of the solar hour angle
        // and solar declination
        //
        // Arguments:
        // tilt           -- is the tilt angle of the inclined surface from horizontal, measured
        //                   upwards facing, 0 to 180, in degrees;
        // orientation    -- is the orientation angle of the inclined surface, expressed as the
        //                   geographical azimuth angle of the horizontal projection of the inclined
        //                   surface normal, -180 to 180, in degrees;
        // simulation_time - an iteration of the current simulation time
        // """

        //set up/ shadow some vars as radians for trig stuff
        let tilt = tilt.to_radians();
        let orientation = orientation.to_radians();
        let latitude = self.latitude.to_radians();
        let solar_declination: f64 =
            self.solar_declinations[simulation_time.current_day() as usize].to_radians();
        let solar_hour_angle: f64 =
            self.solar_hour_angles[simulation_time.current_hour() as usize].to_radians();

        (solar_declination.sin() * latitude.sin() * tilt.cos()
            - solar_declination.sin() * latitude.cos() * tilt.sin() * orientation.cos()
            + solar_declination.cos() * latitude.cos() * tilt.cos() * solar_hour_angle.cos()
            + solar_declination.cos()
                * latitude.sin()
                * tilt.sin()
                * orientation.cos()
                * solar_hour_angle.cos()
            + solar_declination.cos() * tilt.sin() * orientation.sin() * solar_hour_angle.sin())
        .acos()
        .to_degrees()
    }

    /// calculates the azimuth angle between sun and the inclined surface,
    /// needed as input for the calculation of the irradiance in case of solar shading by objects
    ///
    /// Arguments:
    ///
    /// * `orientation` - is the orientation angle of the inclined surface, expressed as the
    ///                   geographical azimuth angle of the horizontal projection of the inclined
    ///                   surface normal, -180 to 180, in degrees;
    #[cfg(test)]
    fn sun_surface_azimuth(&self, orientation: f64, simtime: SimulationTimeIteration) -> f64 {
        let current_hour = simtime.current_hour();
        let test_angle = self.solar_hour_angles[current_hour as usize] - orientation;

        if test_angle > 180. {
            -360. + test_angle
        } else if test_angle < -180. {
            360. + test_angle
        } else {
            test_angle
        }
    }

    /// calculates the tilt angle between sun and the inclined surface,
    /// needed as input for the calculation of the irradiance in case of solar shading by objects
    ///
    /// Arguments:
    ///
    /// * `tilt` - is the tilt angle of the inclined surface from horizontal, measured
    ///            upwards facing, 0 to 180, in degrees;
    #[cfg(test)]
    fn sun_surface_tilt(&self, tilt: f64, simtime: SimulationTimeIteration) -> f64 {
        let current_hour = simtime.current_hour();
        let test_angle = tilt - self.solar_zenith_angles[current_hour as usize];

        if test_angle > 180. {
            -360. + test_angle
        } else if test_angle < -180. {
            360. + test_angle
        } else {
            test_angle
        }
    }

    fn direct_irradiance(
        &self,
        tilt: f64,
        orientation: f64,
        simulation_time: &SimulationTimeIteration,
    ) -> f64 {
        // """  calculates the direct irradiance on the inclined surface, determined as function
        // of cosine of the solar angle of incidence and the direct normal (beam) solar irradiance
        // NOTE The solar beam irradiance is defined as falling on an surface normal to the solar beam.
        // This is not the same as direct horizontal radiation.
        //
        // Arguments:
        // tilt           -- is the tilt angle of the inclined surface from horizontal, measured
        //                   upwards facing, 0 to 180, in degrees;
        // orientation    -- is the orientation angle of the inclined surface, expressed as the
        //                   geographical azimuth angle of the horizontal projection of the inclined
        //                   surface normal, -180 to 180, in degrees;
        //
        // """
        let direct_irradiance = self.direct_beam_radiation(simulation_time.index)
            * self
                .solar_angle_of_incidence(tilt, orientation, simulation_time)
                .to_radians()
                .cos();
        if direct_irradiance < 0.0 {
            0.0
        } else {
            direct_irradiance
        }
    }

    fn a_over_b(
        &self,
        tilt: f64,
        orientation: f64,
        simulation_time: &SimulationTimeIteration,
    ) -> f64 {
        // """  calculates the ratio of the parameters a and b
        //
        // Arguments:
        // tilt           -- is the tilt angle of the inclined surface from horizontal, measured
        //                   upwards facing, 0 to 180, in degrees;
        // orientation    -- is the orientation angle of the inclined surface, expressed as the
        //                   geographical azimuth angle of the horizontal projection of the inclined
        //                   surface normal, -180 to 180, in degrees;
        // """

        // #dimensionless parameters a & b
        // #describing the incidence-weighted solid angle sustained by the circumsolar region as seen
        // #respectively by the tilted surface and the horizontal.
        let a = max_of_2(
            0.,
            self.solar_angle_of_incidence(tilt, orientation, simulation_time)
                .to_radians()
                .cos(),
        );
        let b = max_of_2(
            85.0f64.to_radians().cos(),
            self.solar_zenith_angles[simulation_time.current_hour() as usize]
                .to_radians()
                .cos(),
        );

        a / b
    }

    fn diffuse_irradiance(
        &self,
        tilt: f64,
        orientation: f64,
        simulation_time: &SimulationTimeIteration,
    ) -> DiffuseIrradiance {
        // """  calculates the diffuse part of the irradiance on the surface (without ground reflection)
        //
        // Arguments:
        // tilt           -- is the tilt angle of the inclined surface from horizontal, measured
        //                   upwards facing, 0 to 180, in degrees;
        // orientation    -- is the orientation angle of the inclined surface, expressed as the
        //                   geographical azimuth angle of the horizontal projection of the inclined
        //                   surface normal, -180 to 180, in degrees;
        // """

        // #first set up parameters needed for the calculation
        let gsol_d = self.diffuse_horizontal_radiation(simulation_time.index);
        let f1 = self.f1_circumsolar_brightness_coefficients[simulation_time.index];
        let f2 = self.f2_horizontal_brightness_coefficients[simulation_time.index];

        // # Calculate components of diffuse radiation
        let diffuse_irr_sky = gsol_d * (1.0 - f1) * ((1.0 + tilt.to_radians().cos()) / 2.0);
        let diffuse_irr_circumsolar =
            self.circumsolar_irradiance(tilt, orientation, simulation_time);
        let diffuse_irr_horiz = gsol_d * f2 * tilt.to_radians().sin();

        DiffuseIrradiance(
            diffuse_irr_sky + diffuse_irr_circumsolar + diffuse_irr_horiz,
            diffuse_irr_sky,
            diffuse_irr_circumsolar,
            diffuse_irr_horiz,
        )
    }

    fn ground_reflection_irradiance(
        &self,
        tilt: f64,
        simulation_time: &SimulationTimeIteration,
    ) -> f64 {
        // """  calculates the contribution of the ground reflection to the irradiance on the inclined surface,
        // determined as function of global horizontal irradiance, which in this case is calculated from the solar
        // altitude, diffuse and beam solar irradiance and the solar reflectivity of the ground
        //
        // Arguments:
        // tilt           -- is the tilt angle of the inclined surface from horizontal, measured
        //                   upwards facing, 0 to 180, in degrees;
        // """

        (self.diffuse_horizontal_radiation(simulation_time.index)
            + self.direct_beam_radiation(simulation_time.index)
                * self.solar_altitudes[simulation_time.current_hour() as usize]
                    .to_radians()
                    .sin())
            * self.solar_reflectivity_of_ground(simulation_time)
            * ((1.0 - tilt.to_radians().cos()) / 2.0)
    }

    fn circumsolar_irradiance(
        &self,
        tilt: f64,
        orientation: f64,
        simulation_time: &SimulationTimeIteration,
    ) -> f64 {
        // """  calculates the circumsolar_irradiance
        //
        // Arguments:
        // tilt           -- is the tilt angle of the inclined surface from horizontal, measured
        //                   upwards facing, 0 to 180, in degrees;
        // orientation    -- is the orientation angle of the inclined surface, expressed as the
        //                   geographical azimuth angle of the horizontal projection of the inclined
        //                   surface normal, -180 to 180, in degrees;
        // """

        self.diffuse_horizontal_radiation(simulation_time.index)
            * self.f1_circumsolar_brightness_coefficients[simulation_time.index]
            * self.a_over_b(tilt, orientation, simulation_time)
    }

    fn calculated_direct_irradiance(
        &self,
        tilt: f64,
        orientation: f64,
        simulation_time: &SimulationTimeIteration,
    ) -> f64 {
        // """  calculates the total direct irradiance on an inclined surface including circumsolar
        //
        // Arguments:
        // tilt           -- is the tilt angle of the inclined surface from horizontal, measured
        //                   upwards facing, 0 to 180, in degrees;
        // orientation    -- is the orientation angle of the inclined surface, expressed as the
        //                   geographical azimuth angle of the horizontal projection of the inclined
        //                   surface normal, -180 to 180, in degrees;
        // """

        self.direct_irradiance(tilt, orientation, simulation_time)
            + self.circumsolar_irradiance(tilt, orientation, simulation_time)
    }

    fn calculated_diffuse_irradiance(
        &self,
        tilt: f64,
        orientation: f64,
        simulation_time: &SimulationTimeIteration,
    ) -> f64 {
        // """  calculates the total diffuse irradiance on an inclined surface excluding circumsolar
        // and including ground reflected irradiance
        //
        // Arguments:
        // tilt           -- is the tilt angle of the inclined surface from horizontal, measured
        //                   upwards facing, 0 to 180, in degrees;
        // orientation    -- is the orientation angle of the inclined surface, expressed as the
        //                   geographical azimuth angle of the horizontal projection of the inclined
        //                   surface normal, -180 to 180, in degrees;
        // """

        let DiffuseIrradiance(diffuse_irr_total, _, diffuse_irr_circumsolar, _) =
            self.diffuse_irradiance(tilt, orientation, simulation_time);

        diffuse_irr_total - diffuse_irr_circumsolar
            + self.ground_reflection_irradiance(tilt, simulation_time)
    }

    pub fn calculated_total_solar_irradiance(
        &self,
        tilt: f64,
        orientation: f64,
        simulation_time: &SimulationTimeIteration,
    ) -> f64 {
        // """  calculates the hemispherical or total solar irradiance on the inclined surface
        // without the effect of shading
        //
        // Arguments:
        // tilt           -- is the tilt angle of the inclined surface from horizontal, measured
        //                   upwards facing, 0 to 180, in degrees;
        // orientation    -- is the orientation angle of the inclined surface, expressed as the
        //                   geographical azimuth angle of the horizontal projection of the inclined
        //                   surface normal, -180 to 180, in degrees;
        //
        // """

        self.calculated_direct_irradiance(tilt, orientation, simulation_time)
            + self.calculated_diffuse_irradiance(tilt, orientation, simulation_time)
    }

    pub(crate) fn calculated_direct_diffuse_total_irradiance(
        &self,
        tilt: f64,
        orientation: f64,
        diffuse_breakdown: bool,
        simulation_time: &SimulationTimeIteration,
    ) -> CalculatedDirectDiffuseTotalIrradiance {
        // NB. the original Python implementation uses an internal cache for each timestep here
        // it only retains one set of results at a time so it may be of limited use, and i've skipped reimplementing it
        // in a first pass of conversion into Rust

        let DiffuseIrradiance(
            diffuse_irr_total,
            diffuse_irr_sky,
            diffuse_irr_circumsolar,
            diffuse_irr_horiz,
        ) = self.diffuse_irradiance(tilt, orientation, simulation_time);

        let ground_reflection_irradiance = self.ground_reflection_irradiance(tilt, simulation_time);

        let calculated_direct =
            self.direct_irradiance(tilt, orientation, simulation_time) + diffuse_irr_circumsolar;
        let calculated_diffuse =
            diffuse_irr_total - diffuse_irr_circumsolar + ground_reflection_irradiance;
        let total_irradiance = calculated_direct + calculated_diffuse;

        CalculatedDirectDiffuseTotalIrradiance(
            calculated_direct,
            calculated_diffuse,
            total_irradiance,
            diffuse_breakdown.then_some(DiffuseBreakdown {
                sky: diffuse_irr_sky,
                circumsolar: diffuse_irr_circumsolar,
                horiz: diffuse_irr_horiz,
                ground_refl: ground_reflection_irradiance,
            }),
        )
    }

    fn outside_solar_beam(
        &self,
        tilt: f64,
        orientation: f64,
        simulation_time: &SimulationTimeIteration,
    ) -> bool {
        // """ checks if the shaded surface is in the view of the solar beam.
        // if not, then shading is complete, total direct rad = 0 and no further
        // shading calculation needed for this object for this time step. returns
        // a flag for whether the surface is outside solar beam
        //
        // Arguments:
        // tilt           -- is the tilt angle of the inclined surface from horizontal, measured
        //                   upwards facing, 0 to 180, in degrees;
        // orientation    -- is the orientation angle of the inclined surface, expressed as the
        //                   geographical azimuth angle of the horizontal projection of the
        //                   inclined surface normal, -180 to 180, in degrees;
        //
        // """

        let current_hour_idx = simulation_time.current_hour() as usize;

        let test1 = orientation - self.solar_azimuth_angles[current_hour_idx];
        let test1 = if test1 > 180. {
            test1 - 360.
        } else if test1 < -180. {
            test1 + 360.
        } else {
            test1
        };
        let test2 = tilt - self.solar_altitudes[current_hour_idx];

        !(-90.0..=90.0).contains(&test1) || !(-90.0..=90.0).contains(&test2)
    }

    fn get_segment(
        &self,
        simulation_time: &SimulationTimeIteration,
    ) -> anyhow::Result<ShadingSegment> {
        // """ for complex (environment) shading objects, we need to know which
        // segment the azimuth of the sun occupies at each timestep
        //
        // """

        let current_hour_idx = simulation_time.current_hour() as usize;
        let azimuth = self.solar_azimuth_angles[current_hour_idx];

        let mut previous_segment_end: Option<f64> = None;

        if let Some(shading_segments) = self.shading_segments.as_ref() {
            for segment in shading_segments {
                if let Some(previous_segment_end) = previous_segment_end {
                    if previous_segment_end != segment.start {
                        return Err(anyhow!("No gaps between shading segments allowed"));
                    }
                }
                previous_segment_end = Some(segment.end);
                if segment.end > segment.start {
                    return Err(anyhow!(
                        "End orientation is less than the start orientation. Check shading inputs.",
                    ));
                }
                if azimuth < segment.start && azimuth > segment.end {
                    return Ok(segment.clone());
                }
            }
        }

        Err(anyhow!(
            "Solar segment was not found - this is an unexpected error"
        ))
    }

    fn obstacle_shading_height(
        &self,
        base_height_of_k: f64,
        height_of_obstacle: f64,
        horiz_distance_from_surface_to_obstacle: f64,
        simulation_time: &SimulationTimeIteration,
    ) -> f64 {
        // """ calculates the height of the shading on the shaded surface (k),
        // from the shading obstacle in segment i at time t. Note that "obstacle"
        // has a specific meaning in ISO 52016 Annex F
        //
        // Arguments:
        // base_height_of_k - Hkbase        -- is the base height of the shaded surface k, in m
        // height_of_obstacle - Hobst         -- is the height of the shading obstacle, p, in segment i, in m
        // horiz_distance_from_surface_to_obstacle - Lkobst        -- is the horizontal distance between the shaded surface k, in m
        //                  and the shading obstacle p in segment i, in m
        // """
        let hshade = height_of_obstacle
            - base_height_of_k
            - horiz_distance_from_surface_to_obstacle
                * self.solar_altitudes[simulation_time.current_hour() as usize]
                    .to_radians()
                    .tan();
        if hshade < 0.0 {
            0.0
        } else {
            hshade
        }
    }

    // following references nonexistent solar_altitude method in the original python, so leaving incomplete
    /// calculates the height of the shading on the shaded surface (k),
    /// from the shading overhang in segment i at time t. Note that "overhang"
    /// has a specific meaning in ISO 52016 Annex F
    ///
    /// Arguments:
    /// * `shaded_surface_height` - is the height of the shaded surface, k, in m (referred to as Hk in the Python)
    /// * `base_shaded_surface_height` - is the base height of the shaded surface k, in m (referred to as Hkbase in the Python)
    /// * `lowest_height_of_overhang` - is the lowest height of the overhang q, in segment i, in m (referred to as Hovh in the Python)
    /// * `horiz_distance_from_surface_to_overhang` - is the horizontal distance between the shaded surface k
    ///                                               and the shading overhang, q, in segment i, in m (referred to as Lkovh in the Python)
    fn overhang_shading_height(
        &self,
        shaded_surface_height: f64,
        base_shaded_surface_height: f64,
        lowest_height_of_overhang: f64,
        horiz_distance_from_surface_to_overhang: f64,
        simulation_time: SimulationTimeIteration,
    ) -> f64 {
        let current_hour = simulation_time.current_hour();
        0.0f64.max(
            shaded_surface_height + base_shaded_surface_height - lowest_height_of_overhang
                + horiz_distance_from_surface_to_overhang
                    * self.solar_altitudes[current_hour as usize]
                        .to_radians()
                        .tan(),
        )
    }

    /// calculates the shading factor of direct radiation due to external
    /// shading objects
    ///
    /// Arguments:
    /// * `height` - is the height of the shaded surface (if surface is tilted then
    ///                   this must be the vertical projection of the height), in m
    /// * `base_height` - is the base height of the shaded surface k, in m
    /// * `width` - is the width of the shaded surface, in m
    /// * `orientation` - is the orientation angle of the inclined surface, expressed as the
    ///                   geographical azimuth angle of the horizontal projection of the
    ///                   inclined surface normal, -180 to 180, in degrees;
    /// * `window_shading` - data on overhangs and side fins associated to this building element
    ///                   includes the shading object type, depth, anf distance from element
    pub fn direct_shading_reduction_factor(
        &self,
        base_height: f64,
        height: f64,
        width: f64,
        orientation: f64,
        window_shading: Option<&[WindowShadingObject]>,
        simulation_time: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // start with default assumption of no shading
        let mut hshade_obst = 0.0;
        let mut hshade_ovh = 0.0;
        let mut wfinr = 0.0;
        let mut wfinl = 0.0;

        // #first process the distant (environment) shading for this building element
        let segment = self.get_segment(&simulation_time)?;

        if let Some(shading_objects) = segment.shading_objects {
            for shading_object in shading_objects {
                match shading_object.object_type {
                    ShadingObjectType::Obstacle => {
                        let new_shade_height = self.obstacle_shading_height(
                            base_height,
                            shading_object.height,
                            shading_object.distance,
                            &simulation_time,
                        );
                        hshade_obst = max_of_2(hshade_obst, new_shade_height);
                    }
                    ShadingObjectType::Overhang => {
                        let new_shade_height = self.overhang_shading_height(
                            height,
                            base_height,
                            shading_object.height,
                            shading_object.distance,
                            simulation_time,
                        );

                        hshade_ovh = max_of_2(hshade_ovh, new_shade_height);
                    }
                }
            }
        }

        // then check if there is any simple shading on this building element
        // (note only applicable to transparent building elements so window_shading
        // will always be False for other elements)
        if let Some(window_shading) = window_shading {
            let current_hour = simulation_time.current_hour();
            let altitude = self.solar_altitudes[current_hour as usize];
            let azimuth = self.solar_azimuth_angles[current_hour as usize];
            for shading_object in window_shading {
                match shading_object {
                    WindowShadingObject::Obstacle { .. } => {
                        // For nearby obstacles, skip this loop. These will be dealt with later
                        continue;
                    }
                    WindowShadingObject::Overhang { depth, distance } => {
                        let new_shade_height = (depth * altitude.to_radians().tan()
                            / (azimuth - orientation).to_radians().cos())
                            - distance;
                        hshade_ovh = max_of_2(hshade_ovh, new_shade_height);
                    }
                    WindowShadingObject::SideFinRight { depth, distance } => {
                        // check if the sun is in the opposite direction
                        let check = azimuth - orientation;
                        let new_finrshade = if check > 0. {
                            0.
                        } else {
                            depth * (azimuth - orientation).to_radians().tan() - distance
                        };
                        wfinr = max_of_2(wfinr, new_finrshade);
                    }
                    WindowShadingObject::SideFinLeft { depth, distance } => {
                        // check if the sun is in the opposite direction
                        let check = azimuth - orientation;
                        let new_finlshade = if check < 0. {
                            0.
                        } else {
                            depth * (azimuth - orientation).to_radians().tan() - distance
                        };
                        wfinl = max_of_2(wfinl, new_finlshade);
                    }
                    _ => return Err(anyhow!("unexpected window shading object type encountered")),
                }
            }
        }

        // The height of the shade on the shaded surface from all obstacles is the
        // largest of all, with as maximum value the height of the shaded object
        let hk_obst = min_of_2(height, hshade_obst);

        // The height of the shade on the shaded surface from all overhangs is the
        // largest of all, with as maximum value the height of the shaded object
        let hk_ovh = min_of_2(height, hshade_ovh);

        // The height of the remaining sunlit area on the shaded surface from
        // all obstacles and all overhangs
        let hk_sun = max_of_2(0., height - (hk_obst + hk_ovh));

        // The width of the shade on the shaded surface from all right side fins
        // is the largest of all, with as maximum value the width of the shaded object
        let wk_finr = min_of_2(width, wfinr);

        // The width of the shade on the shaded surface from all left side fins
        // is the largest of all, with as maximum value the width of the shaded object
        let wk_finl = min_of_2(width, wfinl);

        // The width of the remaining sunlit area on the shaded surface from all
        // right hand side fins and all left hand side fins
        let wk_sun = max_of_2(0., width - (wk_finr + wk_finl));

        // And then the direct shading reduction factor of the shaded surface for
        // obstacles, overhangs and side fins
        let mut fdir = (hk_sun * wk_sun) / (height * width);

        if let Some(window_shading) = window_shading {
            for shade_obj in window_shading {
                if let WindowShadingObject::Obstacle {
                    height: shading_height,
                    distance,
                    transparency,
                } = shade_obj
                {
                    let new_shade_height = self.obstacle_shading_height(
                        base_height,
                        *shading_height,
                        *distance,
                        &simulation_time,
                    );

                    let new_shade_trans = *transparency;

                    // repeat Fdir assessment for each near obstacle to find largest shading effect
                    let hk_obst = height.min(new_shade_height);
                    let hk_sun = 0.0f64.max(height - (hk_obst + hk_ovh))
                        + (hk_obst.min(height - hk_ovh) * new_shade_trans);

                    fdir = fdir.min((hk_sun * wk_sun) / (height * width));
                }
            }
        }

        Ok(fdir)
    }

    /// calculates the shading factor of diffuse radiation due to external shading objects
    ///
    /// Arguments:
    /// * `height` - is the height of the shaded surface (if surface is tilted then
    ///              this must be the vertical projection of the height), in m
    /// * `base_height` - is the base height of the shaded surface k, in m
    /// * `width` - is the width of the shaded surface, in m
    /// * `orientation` - is the orientation angle of the inclined surface, expressed as the
    ///                   geographical azimuth angle of the horizontal projection of the
    ///                   inclined surface normal, -180 to 180, in degrees;
    /// * `window_shading` - data on overhangs and side fins associated to this building element
    ///                      includes the shading object type, depth, and distance from element
    fn diffuse_shading_reduction_factor(
        &self,
        diffuse_breakdown: DiffuseBreakdown,
        tilt: f64,
        height: f64,
        base_height: f64,
        width: f64,
        orientation: f64,
        window_shading: Option<&Vec<WindowShadingObject>>,
        f_sky: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // Note: Shading factor for circumsolar radiation is same as for direct.
        //       As circumsolar radiation will be subtracted from diffuse and
        //       added to direct later on, we don't need to do anything for
        //       circumsolar radiation here and it is excluded.
        let diffuse_irr_sky = diffuse_breakdown.sky;
        let diffuse_irr_hor = diffuse_breakdown.horiz;
        let diffuse_irr_ref = diffuse_breakdown.ground_refl;
        let diffuse_irr_total = diffuse_irr_sky + diffuse_irr_hor + diffuse_irr_ref;

        // PD CEN ISO/TR 52016-2:2017 Section F.6.2 is not clearly defined and has been ignored.
        // Calculation for remote obstacles uses a similar method to the simple facade object
        // and self-shading corrections in F6.3. Equation F.7 is now assumed to include a
        // further term of (FvW - FvRM), where FvRM is the view factor between the element
        // and a remote obstacle.
        //
        // Assumes for all elements that the angle between element and baseline horizon is 0.
        // Any significant height variation between element and unobstructed horizon should
        // be considered as an obstacle to be included.
        //
        // Any element height that projects above obstruction is consider unobstructed. The
        // shading angle to the sky is taken from the midpoint of the section of element
        // below the obstruction.

        /// Returns intersection between element shaded arc and shading segment
        fn interval_intersect(a: (f64, f64), b: (f64, f64)) -> f64 {
            0.0f64.max(a.1.min(b.1) - a.0.max(b.0))
        }

        /// Returns angle ranges included in element shaded arc with 0/360 crossover
        /// split if required plus total angle of arc
        fn arc_angle(arc_srt: f64, arc_fsh: f64) -> ArcAngle {
            let (arc1, arc2, deg_arc, rarc1, rarc2) = if arc_srt < arc_fsh {
                // Define front arc as single arc and split rear arc either side of 0/360 boundary
                (
                    (arc_srt, arc_fsh),
                    (0., 0.),
                    arc_fsh - arc_srt,
                    (arc_fsh, 360.),
                    (0., arc_srt),
                )
            } else {
                // Define rear arc as single arc and split front arc either side of 0/360 boundary
                (
                    (arc_srt, 360.),
                    (0., arc_fsh),
                    (360. - arc_srt) + arc_fsh,
                    (arc_fsh, arc_srt),
                    (0., 0.),
                )
            };

            let arc_ang = (arc1, arc2);
            let rarc_ang = (rarc1, rarc2);

            (arc_ang, rarc_ang, deg_arc)
        }

        /// Returns angle ranges included in shading segment with 0/360 crossover
        /// split if required plus total angle of segment
        fn seg_angle(seg_srt: f64, seg_fsh: f64) -> SegAngle {
            let (seg1, seg2, deg_seg) = if seg_srt < seg_fsh {
                ((seg_srt, seg_fsh), (0., 0.), seg_fsh - seg_srt)
            } else {
                // Treat as separate segments either side of 0/360 boundary
                ((seg_srt, 360.), (0., seg_fsh), (360. - seg_srt) + seg_fsh)
            };

            let seg_ang = (seg1, seg2);

            (seg_ang, deg_seg)
        }

        // Determine start and end orientations for potential forward shading arc by remote obstacles
        // 180deg arc assumed unless horizontal

        let (arc_srt, arc_fsh) = if tilt > 0. {
            let orient360 = self.orientation360(orientation);
            if (90. ..=270.).contains(&orient360) {
                (orient360 - 90., orient360 + 90.)
            } else if orient360 < 90. {
                (orient360 + 270., orient360 + 90.)
            } else {
                (orient360 - 90., orient360 - 270.)
            }
        } else {
            (0., 360.)
        };

        // Define arcs to the front and rear of element
        let (arc_ang, rarc_ang, deg_arc) = arc_angle(arc_srt, arc_fsh);

        let mut f_sky_new = 0.;

        if let Some(shading_segments) = self.shading_segments.as_ref() {
            for segment in shading_segments {
                let seg_srt = 180. - segment.start; // Segment start angle (clockwise)
                let seg_fsh = 180. - segment.end; // Segment end angle (clockwise)

                // Define segment
                let (seg_ang, deg_seg) = seg_angle(seg_srt, seg_fsh);

                // Compare sub-arcs and sub-segments for overlap - front
                let ap00 = interval_intersect(arc_ang.0, seg_ang.0);
                let ap11 = interval_intersect(arc_ang.1, seg_ang.1);
                let ap01 = interval_intersect(arc_ang.0, seg_ang.1);
                let ap10 = interval_intersect(arc_ang.1, seg_ang.0);

                // Proportion of front arc shaded by segment
                let arc_prop = (ap00 + ap11 + ap01 + ap10) / deg_arc;

                // Compare sub-arcs and sub-segments for overlap - rear
                let rap00 = interval_intersect(rarc_ang.0, seg_ang.0);
                let rap11 = interval_intersect(rarc_ang.1, seg_ang.1);
                let rap01 = interval_intersect(rarc_ang.0, seg_ang.1);
                let rap10 = interval_intersect(rarc_ang.1, seg_ang.0);

                // Proportion of rearward arc within segment
                let rarc_prop = (rap00 + rap11 + rap01 + rap10) / deg_arc;

                // Segment f_sky contribution in forward direction
                let f_sky_seg_front = if tilt == 0. {
                    f_sky * (deg_seg / 360.)
                } else {
                    arc_prop * 0.5f64.min(f_sky)
                };

                // For tilted surface, segment f_sky contribution in rear direction
                let f_sky_seg_rear = if tilt > 0. && tilt < 90. {
                    rarc_prop * 0.0f64.max(f_sky - 0.5)
                } else {
                    0.
                };

                let mut f_sky_ft = f_sky_seg_front;
                let mut f_sky_rr = f_sky_seg_rear;

                if let Some(shading) = segment.shading_objects.as_ref() {
                    for shade_obj in shading {
                        match shade_obj {
                            ShadingObject {
                                object_type: ShadingObjectType::Obstacle,
                                height: shading_height,
                                distance,
                            } => {
                                // added below bail to match the (erroneous) Python behaviour which
                                // has been reported upstream (schema allows 0 but will error here)
                                if *distance <= 0. {
                                    bail!("distance for shading objects with type 'obstacle' should be greater than zero");
                                }
                                let h_shade = 0.0f64.max(*shading_height - base_height);

                                if f_sky == 1. {
                                    let alpha_obst = (h_shade / *distance).atan().to_degrees();
                                    f_sky_ft = f_sky_ft
                                        .min(f_sky_seg_front * alpha_obst.to_radians().cos());
                                } else if f_sky > 0. {
                                    if f_sky_seg_front > 0. {
                                        // height element is above obstacle (zero if not)
                                        let h_above = 0.0f64.max(height - h_shade);
                                        // proportion of element above obstacle
                                        let p_above = h_above / height;
                                        // angle between midpoint of shaded section and obstacle
                                        let alpha_obst = ((h_shade - (height.min(h_shade) / 2.))
                                            / *distance)
                                            .atan()
                                            .to_degrees();
                                        // Determine if obstacle gives largest reduction to the segment f_sky contribution
                                        f_sky_ft = f_sky_ft.min(0.0f64.max(
                                            f_sky_seg_front
                                                - 0.5
                                                    * arc_prop
                                                    * (1. - alpha_obst.to_radians().cos())
                                                    * (1. - p_above),
                                        ));
                                    }

                                    if f_sky_seg_rear > 0. {
                                        // Determine if potential for shading from obstacles to rear of tilted surface exist
                                        // Projected height of element at obstacle distance
                                        let h_eff = height + (*distance * tilt.to_radians().tan());
                                        if h_eff < h_shade {
                                            // angle from element midpoint to top of obstacle
                                            let alpha_obst =
                                                (h_shade / *distance).atan().to_degrees();
                                            // Determine new rear sky view factor directly from shading angle (alpha_obst) using standard 0.5 x cos(angle) method
                                            // as shading is now determined by this angle not the tilt angle which is smaller
                                            f_sky_rr = f_sky_rr.min(
                                                rarc_prop * 0.5 * alpha_obst.to_radians().cos(),
                                            );
                                        }
                                    }
                                }
                            }
                            ShadingObject {
                                object_type: ShadingObjectType::Overhang,
                                height: shading_height,
                                distance,
                            } => {
                                let h_shade = 0.0f64.max(*shading_height - base_height);

                                if f_sky == 1. {
                                    let alpha_ovh = (h_shade / *distance).atan().to_degrees();
                                    f_sky_ft = f_sky_ft
                                        .min(f_sky_seg_front * (1. - alpha_ovh.to_radians().cos()));
                                } else if f_sky > 0. {
                                    if f_sky_seg_front > 0. {
                                        // height element is below overhang (zero if not)
                                        let h_below = height.min(h_shade);
                                        // proportion of element below overhang
                                        let p_below = h_below / height;
                                        let alpha_ovh = ((h_shade - (height.min(h_shade) / 2.))
                                            / *distance)
                                            .atan()
                                            .to_degrees();
                                        // determine if overhang gives largest reduction to the segment f_sky contribution
                                        f_sky_ft = f_sky_ft.min(
                                            0.5 * arc_prop
                                                * (1. - alpha_ovh.to_radians().cos())
                                                * p_below,
                                        );
                                    }

                                    if f_sky_seg_rear > 0. {
                                        // Determine if potential for shading from overhangs to rear of tilted surface exist
                                        // Projected height of element at overhang distance
                                        let h_eff = height + (*distance * tilt.to_radians().tan());
                                        if h_eff < h_shade {
                                            // angle from element midpoint to top of obstacle
                                            let alpha_ovh =
                                                (h_shade / *distance).atan().to_degrees();
                                            f_sky_rr = f_sky_rr.min(
                                                rarc_prop * 0.5 * tilt.to_radians().cos()
                                                    - alpha_ovh.to_radians().cos(),
                                            );
                                        } else {
                                            f_sky_rr = 0.;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                f_sky_new += f_sky_ft + f_sky_rr;
            }
        }

        // Calculate Fdiff for remote obstacles
        // Allow for tilt = 180deg (i.e. f_sky = 0) case
        let f_sh_dif_rem = if f_sky > 0. {
            1. - ((f_sky - f_sky_new) / f_sky)
        } else {
            1.
        };

        // Assumes for remote objects that the reduction in sky view factor is
        // matched by an equivalent increase in ground reflected irradiance.
        let fdiff_ro = if f_sky != 1. {
            let f_sh_ref_rem = (1. - f_sky_new) / (1. - f_sky);
            (f_sh_dif_rem * (diffuse_irr_sky + diffuse_irr_hor) + f_sh_ref_rem * diffuse_irr_ref)
                / diffuse_irr_total
        } else {
            // Effective tilt angle equivalent of shading
            let angle_eff = ((2. * f_sky_new) - 1.).acos().to_degrees();
            let diffuse_irr_ref_new = self.ground_reflection_irradiance(angle_eff, &simtime);
            (f_sh_dif_rem * (diffuse_irr_sky + diffuse_irr_hor) + diffuse_irr_ref_new)
                / diffuse_irr_total
        };

        //
        // Calculate shading from fins and overhangs (PD CEN ISO/TR 52016-2:2017 Section F.6.3)
        // TODO (from Python) Alpha is not defined in this standard but is possibly the angular
        //      height of the horizon. This would make some sense as raising the
        //      horizon angle would have essentially the same effect as increasing
        //      the tilt of the building element so it would make sense to add it
        //      to beta when calculating the sky view factor. The overall effect
        //      of changing alpha when there are no fins or overhangs seems to be
        //      to decrease the shading factor (meaning more shading) for diffuse
        //      radiation from sky and horizon, and to increase the shading factor
        //      for radiation reflected from the ground (sometimes to above 1),
        //      which would also seem to make sense as in this case there is less
        //      sky and more ground in view than in the basic assumption of
        //      perfectly flat surroundings.
        // TODO (from Python) Should angular height of horizon be user input or derived from
        //      calculation of distant shading objects above? Set to zero for now
        let angular_height_of_horizon = 0.0f64;
        let alpha = angular_height_of_horizon.to_radians();

        // TODO (from Python) Beta is not defined in this standard but is used for tilt of the
        //      building element in BS EN ISO 52016-1:2017, so assuming the same.
        //      This seems to give sensible numbers for the sky view factor
        //      F_w_sky when alpha = 0
        let beta = tilt.to_radians();

        let fdiff_list = if window_shading.is_some() {
            // create lists of diffuse shading factors to keep the largest one
            // in case there are multiple shading objects
            let mut fdiff_list: Vec<f64> = vec![];

            // # Unpack window shading details
            let mut ovh_D_L_ls = vec![vec![0.0, 1.0]]; // # [D,L] - L cannot be zero as this leads to divide-by-zero later on
            let mut finR_D_L_ls = vec![vec![0.0, 1.0]]; //# [D,L] - L cannot be zero as this leads to divide-by-zero later on
            let mut finL_D_L_ls = vec![vec![0.0, 1.0]]; // # [D,L] - L cannot be zero as this leads to divide-by-zero later on
            let mut obs_H_L_ls = vec![vec![0.0, 1.0, 0.0]]; // [H,L,Trans]

            if let Some(shading_objects) = window_shading {
                for shading_object in shading_objects.iter() {
                    match shading_object {
                        WindowShadingObject::Overhang { depth, distance } => {
                            ovh_D_L_ls.push(vec![*depth, *distance]);
                        }
                        WindowShadingObject::SideFinLeft { depth, distance } => {
                            finR_D_L_ls.push(vec![*depth, *distance]);
                        }
                        WindowShadingObject::SideFinRight { depth, distance } => {
                            finL_D_L_ls.push(vec![*depth, *distance]);
                        }
                        WindowShadingObject::Obstacle {
                            height,
                            distance,
                            transparency,
                        } => {
                            obs_H_L_ls.push(vec![*height, *distance, *transparency]);
                        }
                        WindowShadingObject::Reveal { .. } => {
                            bail!("shading object type 'reveal' not allowed in context of calculating diffuse shading reduction factor")
                        }
                    }
                }
            }

            // #the default values should not be used if shading is specified
            if ovh_D_L_ls.len() >= 2 {
                ovh_D_L_ls.remove(0);
            }
            if finR_D_L_ls.len() >= 2 {
                finR_D_L_ls.remove(0);
            }
            if finL_D_L_ls.len() >= 2 {
                finL_D_L_ls.remove(0);
            }
            if obs_H_L_ls.len() >= 2 {
                obs_H_L_ls.remove(0);
            }

            let mut f_sh_dif: f64 = Default::default();
            let mut f_sh_ref: f64 = Default::default();

            // #perform the diff shading calculation for each combination of overhangs and fins
            for iteration_vec in [ovh_D_L_ls, finR_D_L_ls, finL_D_L_ls, obs_H_L_ls]
                .iter()
                .multi_cartesian_product()
            {
                let (ovh_D_L, finR_D_L, finL_D_L, obs_H_L) = (
                    iteration_vec[0],
                    iteration_vec[1],
                    iteration_vec[2],
                    iteration_vec[3],
                );
                let d_ovh = ovh_D_L[0];
                let l_ovh = ovh_D_L[1];
                let d_finL = finL_D_L[0];
                let l_finL = finL_D_L[1];
                let d_finR = finR_D_L[0];
                let l_finR = finR_D_L[1];
                let h_obs = obs_H_L[0];
                let l_obs = obs_H_L[1];
                let t_obs = obs_H_L[2];
                // # Calculate required geometric ratios
                // # Note: PD CEN ISO/TR 52016-2:2017 Section F.6.3 refers to ISO 52016-1:2017
                // #       Section F.5.5.1.6 for the definition of P1 and P2. However, this
                // #       section does not exist. Therefore, these definitions have been
                // #       taken from Section F.3.5.1.2 instead, also supported by Table F.6
                // #       in PD CEN ISO/TR 52016-2:2017. These sources define P1 and P2
                // #       differently for fins and for overhangs so it is assumed that
                // #       should also apply here.
                let p1_ovh = d_ovh / height;
                let p2_ovh = l_ovh / height;
                let p1_finL = d_finL / width;
                let p2_finL = l_finL / width;
                let p1_finR = d_finR / width;
                let p2_finR = l_finR / width;

                // # Calculate view factors (eqns F.15 to F.18) required for eqns F.9 to F.14
                // # Note: The equations in the standard refer to P1 and P2, but as per the
                // #       comment above, there are different definitions of these for fins
                // #       and for overhangs. The decision on which ones to use for each of
                // #       the equations below has been made depending on which of the
                // #       subsequent equations the resulting variables are used in (e.g.
                // #       F_w_s is used to calculate F_sh_dif_fins so we use P1 and P2 for
                // #       fins).
                // # Note: For F_w_r, we could set P1 equal to P1 for fins and P2 equal to
                // #       P1 (not P2) for overhangs, as this appears to be consistent with
                // #       example in Table F.6
                // # F_w_r = 1 - exp(-0.8632 * (P1_fin + P1_ovh))
                // # Note: Formula in standard for view factor to fins seems to assume that
                // #       fins are the same on each side. Therefore, here we take the
                // #       average of this view factor calculated with the dimensions of
                // #       each fin.
                let f_w_s = (0.6514
                    * (1.0 - (p2_finL / (p1_finL.powi(2) + p2_finL.powi(2)).sqrt()))
                    + 0.6514 * (1.0 - (p2_finR / (p1_finR.powi(2) + p2_finR.powi(2)).sqrt())))
                    / 2.0;
                let f_w_o = 0.3282 * (1.0 - (p2_ovh / (p1_ovh.powi(2) + p2_ovh.powi(2)).sqrt()));
                let f_w_sky = (1.0 - (alpha + beta - 90.0f64.to_radians()).sin()) / 2.0;

                // # Calculate denominators of eqns F.9 to F.14
                let view_factor_sky_no_obstacles = (1.0 + beta.cos()) / 2.0;
                let view_factor_ground_no_obstacles = (1.0 - beta.cos()) / 2.0;

                // # Setback and remote obstacles (eqns F.9 and F.10): Top half of each eqn
                // # is view factor to sky (F.9) or ground (F.10) with setback and distant
                // # obstacles
                // # TODO (from Python) Uncomment these lines when definitions of P1 and P2 in formula
                // #      for F_w_r have been confirmed.
                // # if view_factor_sky_no_obstacles == 0:
                // #     # Shading makes no difference if sky not visible (avoid divide-by-zero)
                // #     F_sh_dif_setback = 1.0
                // # else:
                // #     F_sh_dif_setback = (1 - F_w_r) * F_w_sky \
                // #                      / view_factor_sky_no_obstacles
                // # if view_factor_ground_no_obstacles == 0:
                // #     # Shading makes no difference if ground not visible (avoid divide-by-zero)
                // #     F_sh_ref_setback = 1.0
                // # else:
                // #     F_sh_ref_setback = (1 - F_w_r) * (1 - F_w_sky) \
                // #                      / view_factor_ground_no_obstacles
                //
                // # Fins and remote obstacles (eqns F.11 and F.12): Top half of each eqn
                // # is view factor to sky (F.11) or ground (F.12) with fins and distant
                // # obstacles
                let f_sh_dif_fins = if view_factor_sky_no_obstacles == 0.0 {
                    1.0
                } else {
                    (1.0 - f_w_s) * f_w_sky / view_factor_sky_no_obstacles
                };
                let f_sh_ref_fins = if view_factor_ground_no_obstacles == 0.0 {
                    1.0
                } else {
                    (1.0 - f_w_s) * (1.0 - f_w_sky) / view_factor_ground_no_obstacles
                };

                // Overhangs and remote obstacles (eqns F.13 and F.14)
                // Top half of eqn F.13 is view factor to sky with overhangs
                let f_sh_dif_overhangs = if view_factor_sky_no_obstacles == 0.0 {
                    1.0
                } else {
                    (f_w_sky - f_w_o) / view_factor_sky_no_obstacles
                };

                // Top half of eqn F.14 is view factor to ground with distant obstacles,
                // but does not account for overhangs blocking any part of the view of
                // the ground, presumably because this will not happen in the vast
                // majority of cases
                let f_sh_ref_overhangs = if view_factor_ground_no_obstacles == 0.0 {
                    // Shading makes no difference if ground not visible (avoid divide-by-zero)
                    1.0
                } else {
                    (1.0 - f_w_sky) / view_factor_ground_no_obstacles
                };

                // Obstacles adjacent to surface (e.g. balcony rails, garden walls). Not explicitly
                // covered by 52016-1 or -2, therefore derived from first principles and general
                // method basis.

                let net_shade_height = h_obs - base_height;
                let f_sh_dif_obs = if view_factor_sky_no_obstacles == 0. || net_shade_height <= 0. {
                    // Shading makes no difference if sky not visible (avoid divide-by-zero)
                    1.0
                } else {
                    // height of element above obstacle
                    let height_above_obstacle = 0.0f64.max(height - net_shade_height);
                    // proportion of element above obstacle
                    let prop_above_obstacle = height_above_obstacle / height;
                    // angle between midpoint of shaded section and top of obstacle
                    let angle_obst = ((net_shade_height / 2.) / l_obs).atan().to_degrees();
                    // Sky view factor reduction
                    let f_w_ob = view_factor_sky_no_obstacles
                        .min((1. - (90. - angle_obst).to_radians().sin()) * 0.5)
                        * (1. - prop_above_obstacle)
                        * (1. - t_obs);
                    (view_factor_sky_no_obstacles - f_w_ob) / view_factor_sky_no_obstacles
                };

                // The impact of obstacles on ground reflected irradiance is difficult to calculate
                // as there is no defined reference ground distance to determine shading impact.
                // A reflected shading reduction factor of 1 is therefore assumed for obstacles
                // on the assumption that any ground reflected irradiance lost will be offset by
                // diffuse reflectance from the obstacle.
                let f_sh_ref_obs = 1.0;

                // Keep the smallest of the three shading reduction factors as the
                // diffuse or reflected shading factor. Also enforce that these cannot be
                // negative (which may happen with some extreme tilt values)
                // TODO (from Python) Add setback shading factors to the arguments to min function when
                //      definitions of P1 and P2 in formula for F_w_r have been confirmed.
                // F_sh_dif = max(0.0, min(F_sh_dif_setback, F_sh_dif_fins, F_sh_dif_overhangs))
                // F_sh_ref = max(0.0, min(F_sh_ref_setback, F_sh_ref_fins, F_sh_ref_overhangs))
                f_sh_dif = max_of_2(
                    0.,
                    [f_sh_dif_fins, f_sh_dif_overhangs, f_sh_dif_obs]
                        .into_iter()
                        .min_by(|a, b| a.total_cmp(b))
                        .unwrap(),
                );
                f_sh_ref = max_of_2(
                    0.,
                    [f_sh_ref_fins, f_sh_ref_overhangs, f_sh_ref_obs]
                        .into_iter()
                        .min_by(|a, b| a.total_cmp(b))
                        .unwrap(),
                );
            }

            if diffuse_irr_total == 0. {
                bail!("Zero diffuse radiation with non-zero direct radiation.");
            }
            let fdiff = (f_sh_dif * (diffuse_irr_sky + diffuse_irr_hor)
                + f_sh_ref * diffuse_irr_ref)
                / diffuse_irr_total;
            fdiff_list.push(fdiff);

            fdiff_list
        } else {
            vec![1.]
        };

        let fdiff = fdiff_list.iter().min_by(|a, b| a.total_cmp(b)).ok_or_else(|| anyhow!("Diffuse shading reduction factor could not be calculated as fdiff_list was empty."))?;
        Ok(fdiff.min(fdiff_ro))
    }

    pub(crate) fn shading_reduction_factor_direct_diffuse(
        &self,
        base_height: f64,
        height: f64,
        width: f64,
        tilt: f64,
        orientation: f64,
        window_shading: &[WindowShadingObject],
        simulation_time: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        // """ calculates the direct and diffuse shading factors due to external
        // shading objects
        //
        // Arguments:
        // height         -- is the height of the shaded surface (if surface is tilted then
        //                   this must be the vertical projection of the height), in m
        // base_height    -- is the base height of the shaded surface k, in m
        // width          -- is the width of the shaded surface, in m
        // orientation    -- is the orientation angle of the inclined surface, expressed as the
        //                   geographical azimuth angle of the horizontal projection of the
        //                   inclined surface normal, -180 to 180, in degrees;
        // tilt           -- is the tilt angle of the inclined surface from horizontal, measured
        //                   upwards facing, 0 to 180, in degrees;
        // window_shading -- data on overhangs and side fins associated to this building element
        //                   includes the shading object type, depth, anf distance from element
        // """

        // # first check if there is any radiation. This is needed to prevent a potential
        // # divide by zero error in the final step, but also, if there is no radiation
        // # then shading is irrelevant and we can skip the whole calculation
        let CalculatedDirectDiffuseTotalIrradiance(direct, diffuse, _, diffuse_breakdown) = self
            .calculated_direct_diffuse_total_irradiance(tilt, orientation, true, &simulation_time);
        if direct + diffuse == 0.0 {
            return Ok((0.0, 0.0));
        }

        let mut window_shading_expanded: Vec<WindowShadingObject> = vec![];
        for shading in window_shading {
            if let WindowShadingObject::Reveal { depth, distance } = shading {
                window_shading_expanded.push(WindowShadingObject::Overhang {
                    depth: *depth,
                    distance: *distance,
                });
                window_shading_expanded.push(WindowShadingObject::SideFinLeft {
                    depth: *depth,
                    distance: *distance,
                });
                window_shading_expanded.push(WindowShadingObject::SideFinRight {
                    depth: *depth,
                    distance: *distance,
                });
            } else {
                window_shading_expanded.push(*shading);
            }
        }

        // # first check if the surface is outside the solar beam
        // # if so then direct shading is complete and we don't need to
        // # calculate shading from objects
        // # TODO (from Python): The outside solar beam condition is based on a vertical projection
        //                       of the surface and does not account for the condition where a
        //                       surface that is only slightly pitched is exposed to direct solar
        //                       radiation when the sun is high (e.g. a surface pitched slightly
        //                       to the north will be exposed to direct solar radiation when the
        //                       sun is high in the southern sky). As the solar radiation
        //                       calculation already accounts for the situation where the sun is
        //                       actually behind the surface (accounting for the combination of
        //                       orientation and pitch), there is no need to zero it using the
        //                       shading factor. For now, we set the shading factor to 1 and ignore
        //                       shading from objects on the other side of the building (which if
        //                       significantly pitched would have to be relatively tall and/or
        //                       very close to cast a shadow on the surface in question anyway),
        //                       so that results in the unshaded case will be correct.
        let fdir = if self.outside_solar_beam(tilt, orientation, &simulation_time) {
            1.0
        } else {
            self.direct_shading_reduction_factor(
                base_height,
                height,
                width,
                orientation,
                Some(&window_shading_expanded),
                simulation_time,
            )?
        };

        let f_sky = sky_view_factor(&tilt);
        let fdiff = self.diffuse_shading_reduction_factor(
            diffuse_breakdown.expect("expected diffuse breakdown to be available"),
            tilt,
            height,
            base_height,
            width,
            orientation,
            Some(&window_shading_expanded),
            f_sky,
            simulation_time,
        )?;

        Ok((fdir, fdiff))
    }

    pub(crate) fn surface_irradiance(
        &self,
        base_height: f64,
        projected_height: f64,
        width: f64,
        tilt: f64,
        orientation: f64,
        window_shading: &[WindowShadingObject],
        simulation_time: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let CalculatedDirectDiffuseTotalIrradiance(i_sol_dir, i_sol_dif, _, _) = self
            .calculated_direct_diffuse_total_irradiance(tilt, orientation, false, &simulation_time);
        let (f_sh_dir, f_sh_dif) = self.shading_reduction_factor_direct_diffuse(
            base_height,
            projected_height,
            width,
            tilt,
            orientation,
            window_shading,
            simulation_time,
        )?;

        Ok(i_sol_dif * f_sh_dif + i_sol_dir * f_sh_dir)
    }

    #[cfg(feature = "fhs")]
    pub fn sun_above_horizon(&self, simtime: SimulationTimeIteration) -> bool {
        let solar_angle = self.solar_angle_of_incidence(0., 0., &simtime);
        solar_angle < 90.
    }
}

type ArcAngle = (((f64, f64), (f64, f64)), ((f64, f64), (f64, f64)), f64);
type SegAngle = (((f64, f64), (f64, f64)), f64);

#[derive(PartialEq, Debug)]
pub(crate) struct CalculatedDirectDiffuseTotalIrradiance(
    pub(crate) f64,
    pub(crate) f64,
    pub(crate) f64,
    pub(crate) Option<DiffuseBreakdown>,
);

// implement traits for comparing floats within CalculatedDirectDiffuseTotalIrradiance in tests
#[cfg(test)]
impl AbsDiffEq<Self> for CalculatedDirectDiffuseTotalIrradiance {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-8
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
            && self.1.abs_diff_eq(&other.1, epsilon)
            && self.2.abs_diff_eq(&other.2, epsilon)
            && match (&self.3, &other.3) {
                (Some(a), Some(b)) => a.abs_diff_eq(b, epsilon),
                (None, None) => true,
                _ => false,
            }
    }
}

#[cfg(test)]
impl RelativeEq for CalculatedDirectDiffuseTotalIrradiance {
    fn default_max_relative() -> Self::Epsilon {
        1e-8
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        self.0.relative_eq(&other.0, epsilon, max_relative)
            && self.1.relative_eq(&other.1, epsilon, max_relative)
            && self.2.relative_eq(&other.2, epsilon, max_relative)
            && match (&self.3, &other.3) {
                (Some(a), Some(b)) => a.relative_eq(b, epsilon, max_relative),
                (None, None) => true,
                _ => false,
            }
    }
}

#[derive(PartialEq, Debug)]
struct DiffuseIrradiance(f64, f64, f64, f64);

// implement traits for comparing floats within diffuse radiances in tests
#[cfg(test)]
impl AbsDiffEq<Self> for DiffuseIrradiance {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-8
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.0.abs_diff_eq(&other.0, epsilon)
            && self.1.abs_diff_eq(&other.1, epsilon)
            && self.2.abs_diff_eq(&other.2, epsilon)
            && self.3.abs_diff_eq(&other.3, epsilon)
    }
}

#[cfg(test)]
impl RelativeEq for DiffuseIrradiance {
    fn default_max_relative() -> Self::Epsilon {
        1e-8
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        self.0.relative_eq(&other.0, epsilon, max_relative)
            && self.1.relative_eq(&other.1, epsilon, max_relative)
            && self.2.relative_eq(&other.2, epsilon, max_relative)
            && self.3.relative_eq(&other.3, epsilon, max_relative)
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
pub(crate) struct DiffuseBreakdown {
    sky: f64,
    circumsolar: f64,
    horiz: f64,
    ground_refl: f64,
}

// implement traits for comparing floats within diffuse breakdowns in tests
#[cfg(test)]
impl AbsDiffEq<Self> for DiffuseBreakdown {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-8
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.sky.abs_diff_eq(&other.sky, epsilon)
            && self.circumsolar.abs_diff_eq(&other.circumsolar, epsilon)
            && self.horiz.abs_diff_eq(&other.horiz, epsilon)
            && self.ground_refl.abs_diff_eq(&other.ground_refl, epsilon)
    }
}

#[cfg(test)]
impl RelativeEq for DiffuseBreakdown {
    fn default_max_relative() -> Self::Epsilon {
        1e-8
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        self.sky.relative_eq(&other.sky, epsilon, max_relative)
            && self
                .circumsolar
                .relative_eq(&other.circumsolar, epsilon, max_relative)
            && self.horiz.relative_eq(&other.horiz, epsilon, max_relative)
            && self
                .ground_refl
                .relative_eq(&other.ground_refl, epsilon, max_relative)
    }
}

pub fn create_external_conditions(
    input: ExternalConditionsInput,
    simulation_time: &SimulationTimeIterator,
) -> anyhow::Result<ExternalConditions> {
    // TODO (from Python) Some inputs are not currently used, so set to None here rather
    //       than requiring them in input file.
    // TODO (from Python) Read timezone from input file. For now, set timezone to 0 (GMT)

    // Let direct beam conversion input be optional, this will be set if comes from weather file.
    let dir_beam_conversion = input.direct_beam_conversion_needed.unwrap_or(false);

    Ok(ExternalConditions::new(
        simulation_time,
        input.air_temperatures.ok_or_else(|| anyhow!("Air temperatures for external conditions were not available when expected"))?,
        input.wind_speeds.ok_or_else(|| anyhow!("Wind speeds for external conditions were not available when expected"))?,
        input.wind_directions.ok_or_else(|| anyhow!("Wind directions for external conditions were not available when expected"))?,
        input
            .diffuse_horizontal_radiation
            .ok_or_else(|| anyhow!("Diffuse horizontal radiation values for external conditions were not available when expected"))?,
        input.direct_beam_radiation.ok_or_else(|| anyhow!("Direct beam radiation values for external conditions were not available when expected"))?,
        input
            .solar_reflectivity_of_ground
            .ok_or_else(|| anyhow!("Solar reflectivity of ground values for external conditions were not available when expected"))?,
        input.latitude.ok_or_else(|| anyhow!("Latitude for external conditions were not available when expected"))?,
        input.longitude.ok_or_else(|| anyhow!("Longitude for external conditions were not available when expected"))?,
        0,
        0,
        Some(365),
        1.,
        None,
        None,
        false,
        dir_beam_conversion,
        input.shading_segments,
    ))
}

fn init_direct_beam_radiation(
    direct_beam_conversion_needed: bool,
    raw_value: f64,
    solar_altitude: f64,
) -> f64 {
    // # if the climate data to only provide direct horizontal (rather than normal:
    // # If only direct (beam) solar irradiance at horizontal plane is available in the climatic data set,
    // # it shall be converted to normal incidence by dividing the value by the sine of the solar altitude.
    // """ ISO 52010 section 6.4.2
    // TODO (from Python) investigate the impact of these notes further. Applicable for weather from CIBSE file.
    // NOTE 1 If the solar altitude angle is low, this conversion is very sensative for tiny
    // errors in the calculation of the solar altitude. Such tiny errors are feasible given the
    // sensitivity for the parameters needed to calculate the solar angle and given the atmospheric
    // refraction of solar radiation near the ground. there fore the value at normal incidence is
    // preferred.
    // NOTE 2 method 1 proved to be most effective in mid-latitude climates
    // other models might be more suitable for tropical climates.
    // NOTE 3 if the solar altitude angle is low, the conversion from direct horizontal to direct
    // normal beam irradiance is very sensitive for tiny errors in the calculation of the
    // solar altitude."""
    if direct_beam_conversion_needed {
        let sin_asol = solar_altitude.to_radians().sin();
        if sin_asol > 0.0 {
            raw_value / sin_asol
        } else {
            raw_value // # TODO (from Python) should this be zero?
        }
    } else {
        raw_value
    }
}

fn init_earth_orbit_deviation(current_day: u32) -> f64 {
    let current_day = current_day + 1; //use 1-indexed day for this

    (360.0 / 365.0) * current_day as f64
}

fn init_solar_declination(earth_orbit_deviation: f64) -> f64 {
    //earth_orbit_deviation passed as degrees; shadow internally though as radians for trig functions
    let earth_orbit_deviation = earth_orbit_deviation.to_radians();

    0.33281
        - 22.984 * earth_orbit_deviation.cos()
        - 0.3499 * (2.0 * earth_orbit_deviation).cos()
        - 0.1398 * (3.0 * earth_orbit_deviation).cos()
        + 3.7872 * earth_orbit_deviation.sin()
        + 0.03205 * (2.0 * earth_orbit_deviation).sin()
        + 0.07187 * (3.0 * earth_orbit_deviation).sin()
}

fn init_equation_of_time(current_day: u32) -> f64 {
    // """ Calculate the equation of time """
    //
    // """
    // teq is the equation of time, in minutes;
    // nday is the day of the year, from 1 to 365 or 366 (leap year)
    // """
    let current_day = (current_day + 1) as i32; //use current_day as 1-indexed day of the year, and make signed

    // # note we convert the values inside the cos() to radians for the python function
    // # even though the 180 / pi is converting from radians into degrees
    // # this way the formula remains consistent with as written in the ISO document
    match current_day {
        nday if current_day < 21 => 2.6 + 0.44 * nday as f64,
        nday if current_day < 136 => 5.2 + 9.0 * ((nday - 43) as f64 * 0.0357).cos(),
        nday if current_day < 241 => 1.4 - 5.0 * ((nday - 135) as f64 * 0.0449).cos(),
        nday if current_day < 336 => -6.3 - 10.0 * ((nday - 306) as f64 * 0.036).cos(),
        nday if current_day <= 366 => 0.45 * (nday - 359) as f64,
        _ => panic!(),
    }
}

fn init_time_shift(timezone: i32, longitude: f64) -> f64 {
    // """ Calculate the time shift, in hours, resulting from the fact that the
    // longitude and the path of the sun are not equal
    //
    // NOTE Daylight saving time is disregarded in tshift which is time independent
    // """
    timezone as f64 - longitude / 15.0
}

fn init_solar_time(
    hour_of_day: u32,
    equation_of_time: f64,
    time_shift: f64,
    _current_hour: u32,
) -> f64 {
    // """ Calculate the solar time, tsol, as a function of the equation of time,
    // the time shift and the hour of the day """

    // #note we +1 here because the simulation hour of day starts at 0
    // #while the sun path standard hour of day starts at 1 (hour 0 to 1)

    let hour_of_day = hour_of_day + 1;

    hour_of_day as f64 - (equation_of_time / 60.0) - time_shift
}

fn init_solar_hour_angle(solar_time: f64) -> f64 {
    // """ Calculate the solar hour angle, in the middle of the
    // current hour as a function of the solar time """
    //
    // TODO (from Python) How is this to be adjusted for timesteps that are not hourly?
    //      would allowing solar_time to be a decimal be all that is needed?
    //
    // """
    // w is the solar hour angle, in degrees
    //
    // Notes from ISO 52020 6.4.1.5
    // NOTE 1 The limitation of angles ranging between -180 and +180 degrees is
    // needed to determine which shading objects are in the direction of the sun;
    // see also the calculation of the azimuth angle of the sun in 6.4.1.7.
    // NOTE 2 Explanation of "12.5": The hour numbers are actually hour sections:
    // the first hour section of a day runs from 0h to 1h. So, the average position
    // of the sun for the solar radiation measured during (solar) hour section N is
    // at (solar) time = (N -0,5) h of the (solar) day.
    // """
    let mut solar_angle = (180 / 12) as f64 * (12.5 - solar_time);

    if solar_angle > 180.0 {
        solar_angle -= 360.0;
    } else if solar_angle < -180.0 {
        solar_angle += 360.0;
    }

    solar_angle
}

fn init_solar_altitude(latitude: f64, solar_declination: f64, solar_hour_angle: f64) -> f64 {
    // """  the angle between the solar beam and the horizontal surface, determined
    //  in the middle of the current hour as a function of the solar hour angle,
    //  the solar declination and the latitude """
    //
    // TODO (from Python) How is this to be adjusted for timesteps that are not hourly?
    //      would allowing solar_time to be a decimal be all that is needed?
    //
    // """
    // asol is the solar altitude angle, the angle between the solar beam
    // and the horizontal surface, in degrees;
    // """

    //all three params provided as degrees, but we need to shadow each as radians for trig calcs
    let latitude = latitude.to_radians();
    let solar_declination = solar_declination.to_radians();
    let solar_hour_angle = solar_hour_angle.to_radians();

    let asol = (solar_declination.sin() * latitude.sin()
        + solar_declination.cos() * latitude.cos() * solar_hour_angle.cos())
    .asin()
    .to_degrees();

    if asol < 0.0001 {
        return 0.;
    }

    asol
}

fn init_solar_zenith_angle(solar_altitude: f64) -> f64 {
    90.0 - solar_altitude
}

fn init_solar_azimuth_angle(
    latitude: f64,
    solar_declination: f64,
    solar_hour_angle: f64,
    solar_altitude: f64,
) -> f64 {
    // """  calculates the solar azimuth angle,
    // angle from South, eastwards positive, westwards negative, in degrees """
    //
    // """
    // NOTE The azimuth angles range between −180 and +180 degrees; this is needed to determine which shading
    // objects are in the direction of the sun
    // """

    //shadow all degrees as radians for trig functions - solar_hour_angle needs to be subtracted from 180 before conversion
    let latitude = latitude.to_radians();
    let solar_declination = solar_declination.to_radians();
    let solar_hour_angle = (180.0 - solar_hour_angle).to_radians();
    let solar_altitude = solar_altitude.to_radians();

    //now do calculation!
    let sin_aux1_numerator = solar_declination.cos() * solar_hour_angle.sin();
    let cos_aux1_numerator = latitude.cos() * solar_declination.sin()
        + latitude.sin() * solar_declination.cos() * solar_hour_angle.cos();

    let denominator = solar_altitude.sin().asin().cos();

    let sin_aux1 = sin_aux1_numerator / denominator;
    let cos_aux1 = cos_aux1_numerator / denominator;
    let aux2 = (sin_aux1_numerator.asin() / denominator).to_degrees();

    // # BS EN ISO 52010-1:2017. Formula 16
    if sin_aux1 >= 0.0 && cos_aux1 > 0.0 {
        if aux2 > 180.0 {
            aux2 - 180.0
        } else {
            180.0 - aux2
        }
    } else if cos_aux1 < 0.0 {
        aux2
    } else {
        -(180.0 + aux2)
    }
}

fn init_air_mass(solar_altitude: f64) -> f64 {
    // """  calculates the air mass, m, the distance the solar beam travels through the earth atmosphere.
    // The air mass is determined as a function of the sine of the solar altitude angle """

    if solar_altitude >= 10.0 {
        1.0 / solar_altitude.to_radians().sin()
    } else {
        1.0 / (solar_altitude.to_radians().sin() + 0.15 * (solar_altitude + 3.885).powf(-1.253))
    }
}

fn init_extra_terrestrial_radiation(earth_orbit_deviation: f64) -> f64 {
    // #NOTE the ISO 52010 has an error in this formula.
    // #it lists Gsol,c as the solar angle of incidence on the inclined surface
    // #when it should be the solar constant, given elsewhere as 1367
    // #we use the correct version of the formula here
    1367.0 * (1.0 + 0.033 * earth_orbit_deviation.to_radians().cos())
}

enum BrightnessCoefficientName {
    F11,
    F12,
    F13,
    F21,
    F22,
    F23,
}

struct BrightnessCoefficientsRow {
    f11: f64,
    f12: f64,
    f13: f64,
    f21: f64,
    f22: f64,
    f23: f64,
}

// version of Table 8 in ISO 52010
static BRIGHTNESS_COEFFICIENTS: [BrightnessCoefficientsRow; 8] = [
    BrightnessCoefficientsRow {
        f11: -0.008,
        f12: 0.588,
        f13: -0.062,
        f21: -0.06,
        f22: 0.072,
        f23: -0.022,
    },
    BrightnessCoefficientsRow {
        f11: 0.13,
        f12: 0.683,
        f13: -0.151,
        f21: -0.019,
        f22: 0.066,
        f23: -0.029,
    },
    BrightnessCoefficientsRow {
        f11: 0.33,
        f12: 0.487,
        f13: -0.221,
        f21: 0.055,
        f22: -0.064,
        f23: -0.026,
    },
    BrightnessCoefficientsRow {
        f11: 0.568,
        f12: 0.187,
        f13: -0.295,
        f21: 0.109,
        f22: -0.152,
        f23: -0.014,
    },
    BrightnessCoefficientsRow {
        f11: 0.873,
        f12: -0.392,
        f13: -0.362,
        f21: 0.226,
        f22: -0.462,
        f23: 0.001,
    },
    BrightnessCoefficientsRow {
        f11: 1.132,
        f12: -1.237,
        f13: -0.412,
        f21: 0.288,
        f22: -0.823,
        f23: 0.056,
    },
    BrightnessCoefficientsRow {
        f11: 1.06,
        f12: -1.6,
        f13: -0.359,
        f21: 0.264,
        f22: -1.127,
        f23: 0.131,
    },
    BrightnessCoefficientsRow {
        f11: 0.678,
        f12: -0.327,
        f13: -0.25,
        f21: 0.156,
        f22: -1.377,
        f23: 0.251,
    },
];

fn brightness_coefficient(e: f64, fij: BrightnessCoefficientName) -> f64 {
    // """ returns brightness coefficient as a look up from Table 8 in ISO 52010
    //
    // Arguments:
    // e    -- dimensionless clearness parameter
    // fij  -- the coefficient to be returned. e.g. f12 or f23
    // """
    //
    // TODO (from Python) I've not had a need for the clearness index parameters contained in this table yet,
    //      if they are needed as input or output later then this function can be reworked
    let row = &BRIGHTNESS_COEFFICIENTS[if e < 1.065 {
        0usize
    } else if e < 1.23 {
        1usize
    } else if e < 1.5 {
        2usize
    } else if e < 1.95 {
        3usize
    } else if e < 2.8 {
        4usize
    } else if e < 4.5 {
        5usize
    } else if e < 6.2 {
        6usize
    } else {
        7usize
    }];
    match fij {
        BrightnessCoefficientName::F11 => row.f11,
        BrightnessCoefficientName::F12 => row.f12,
        BrightnessCoefficientName::F13 => row.f13,
        BrightnessCoefficientName::F21 => row.f21,
        BrightnessCoefficientName::F22 => row.f22,
        BrightnessCoefficientName::F23 => row.f23,
    }
}

fn init_f1_circumsolar_brightness_coefficient(e: f64, delta: f64, solar_zenith_angle: f64) -> f64 {
    // """ returns the circumsolar brightness coefficient, F1
    //
    // Arguments:
    // E -- dimensionless clearness parameter for the current timestep
    // delta -- dimensionless sky brightness parameter for the current timestep
    // solar_zenith_angle -- solar zenith angle for the current hour
    // """
    let f1: f64 = brightness_coefficient(e, BrightnessCoefficientName::F11)
        + brightness_coefficient(e, BrightnessCoefficientName::F12) * delta
        + brightness_coefficient(e, BrightnessCoefficientName::F13)
            * (std::f64::consts::PI * solar_zenith_angle / 180.0);
    if f1 < 0.0 {
        0.0
    } else {
        f1
    }
}

fn init_f2_horizontal_brightness_coefficient(e: f64, delta: f64, solar_zenith_angle: f64) -> f64 {
    // """ returns the horizontal brightness coefficient, F2
    //
    // Arguments:
    // E -- dimensionless clearness parameter
    // delta -- dimensionless sky brightness parameter
    // solar_zenith_angle -- solar zenith angle for the current hour
    // """
    brightness_coefficient(e, BrightnessCoefficientName::F21)
        + brightness_coefficient(e, BrightnessCoefficientName::F22) * delta
        + brightness_coefficient(e, BrightnessCoefficientName::F23)
            * (std::f64::consts::PI * solar_zenith_angle / 180.0)
}

const CLEARNESS_FORMULA_K: f64 = 1.014;

fn init_dimensionless_clearness_parameter(
    diffuse_horizontal_radiation: f64,
    direct_beam_radiation: f64,
    solar_altitude: f64,
) -> f64 {
    // returns the dimensionless clearness parameter, E, anisotropic sky conditions (Perez model)
    if diffuse_horizontal_radiation == 0.0 {
        999.0
    } else {
        (((diffuse_horizontal_radiation + direct_beam_radiation) / diffuse_horizontal_radiation)
            + CLEARNESS_FORMULA_K * (std::f64::consts::PI / 180.0 * solar_altitude).powi(3))
            / (1.0 + CLEARNESS_FORMULA_K * (std::f64::consts::PI / 180.0 * solar_altitude).powi(3))
    }
}

fn init_dimensionless_sky_brightness_parameter(
    air_mass: f64,
    diffuse_horizontal_radiation: f64,
    extra_terrestrial_radiation: f64,
) -> f64 {
    // """  calculates the dimensionless sky brightness parameter, delta
    //
    // Arguments:
    // air_mass -- air mass for the current hour
    // diffuse_horizontal_radiation -- diffuse horizontal radiation for the current timestep
    // extra_terrestrial_radiation -- extra-terrestrial radiation for the current day
    // """
    air_mass * diffuse_horizontal_radiation / extra_terrestrial_radiation
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::units::DAYS_IN_MONTH;
    use crate::external_conditions::DaylightSavingsConfig::NotApplicable;
    use crate::simulation_time::{SimulationTime, HOURS_IN_DAY};
    use approx::assert_relative_eq;
    use pretty_assertions::assert_eq;
    use rstest::*;

    const BASE_AIR_TEMPS: [f64; 24] = [
        0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 7.5, 10.0, 12.5, 15.0, 19.5, 17.0,
        15.0, 12.0, 10.0, 7.0, 5.0, 3.0, 1.0,
    ];

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(7.0, 15.0, 1.0)
    }

    #[fixture]
    fn air_temp_day_jan() -> [f64; 24] {
        BASE_AIR_TEMPS
    }

    #[fixture]
    fn air_temp_day_feb() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 1.0)
    }

    #[fixture]
    fn air_temp_day_mar() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 2.0)
    }

    #[fixture]
    fn air_temp_day_apr() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 3.0)
    }

    #[fixture]
    fn air_temp_day_may() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 4.0)
    }

    #[fixture]
    fn air_temp_day_jun() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 5.0)
    }

    #[fixture]
    fn air_temp_day_jul() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 6.0)
    }

    #[fixture]
    fn air_temp_day_aug() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 6.0)
    }

    #[fixture]
    fn air_temp_day_sep() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 5.0)
    }

    #[fixture]
    fn air_temp_day_oct() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 4.0)
    }

    #[fixture]
    fn air_temp_day_nov() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 3.0)
    }

    #[fixture]
    fn air_temp_day_dec() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 2.0)
    }

    #[fixture]
    fn air_temps() -> Vec<f64> {
        let mut temps: Vec<f64> = vec![];
        let months = [
            (air_temp_day_jan as fn() -> [f64; 24], DAYS_IN_MONTH[0]),
            (air_temp_day_feb as fn() -> [f64; 24], DAYS_IN_MONTH[1]),
            (air_temp_day_mar as fn() -> [f64; 24], DAYS_IN_MONTH[2]),
            (air_temp_day_apr as fn() -> [f64; 24], DAYS_IN_MONTH[3]),
            (air_temp_day_may as fn() -> [f64; 24], DAYS_IN_MONTH[4]),
            (air_temp_day_jun as fn() -> [f64; 24], DAYS_IN_MONTH[5]),
            (air_temp_day_jul as fn() -> [f64; 24], DAYS_IN_MONTH[6]),
            (air_temp_day_aug as fn() -> [f64; 24], DAYS_IN_MONTH[7]),
            (air_temp_day_sep as fn() -> [f64; 24], DAYS_IN_MONTH[8]),
            (air_temp_day_oct as fn() -> [f64; 24], DAYS_IN_MONTH[9]),
            (air_temp_day_nov as fn() -> [f64; 24], DAYS_IN_MONTH[10]),
            (air_temp_day_dec as fn() -> [f64; 24], DAYS_IN_MONTH[11]),
        ];
        for (temps_fn, month_days_count) in months {
            temps.extend_from_slice(
                temps_fn()
                    .to_vec()
                    .iter()
                    .cloned()
                    .cycle()
                    .take((month_days_count * HOURS_IN_DAY) as usize)
                    .collect::<Vec<f64>>()
                    .as_slice(),
            );
        }

        temps
    }

    const BASE_WIND_SPEEDS: [f64; 24] = [
        4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.7, 5.4, 5.6, 5.3, 5.1, 4.8, 4.7, 4.6, 4.5, 4.2,
        4.9, 4.3, 4.4, 4.5, 4.3, 4.6,
    ];

    #[fixture]
    fn wind_speed_day_jan() -> [f64; 24] {
        BASE_WIND_SPEEDS
    }

    #[fixture]
    fn wind_speed_day_feb() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 0.1)
    }

    #[fixture]
    fn wind_speed_day_mar() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 0.2)
    }

    #[fixture]
    fn wind_speed_day_apr() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 0.6)
    }

    #[fixture]
    fn wind_speed_day_may() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 0.8)
    }

    #[fixture]
    fn wind_speed_day_jun() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 1.1)
    }

    #[fixture]
    fn wind_speed_day_jul() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 1.2)
    }

    #[fixture]
    fn wind_speed_day_aug() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 1.2)
    }

    #[fixture]
    fn wind_speed_day_sep() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 1.1)
    }

    #[fixture]
    fn wind_speed_day_oct() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 0.7)
    }

    #[fixture]
    fn wind_speed_day_nov() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 0.5)
    }

    #[fixture]
    fn wind_speed_day_dec() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 0.3)
    }

    #[fixture]
    fn wind_speeds() -> Vec<f64> {
        let mut speeds: Vec<f64> = vec![];
        let months = [
            (wind_speed_day_jan as fn() -> [f64; 24], DAYS_IN_MONTH[0]),
            (wind_speed_day_feb as fn() -> [f64; 24], DAYS_IN_MONTH[1]),
            (wind_speed_day_mar as fn() -> [f64; 24], DAYS_IN_MONTH[2]),
            (wind_speed_day_apr as fn() -> [f64; 24], DAYS_IN_MONTH[3]),
            (wind_speed_day_may as fn() -> [f64; 24], DAYS_IN_MONTH[4]),
            (wind_speed_day_jun as fn() -> [f64; 24], DAYS_IN_MONTH[5]),
            (wind_speed_day_jul as fn() -> [f64; 24], DAYS_IN_MONTH[6]),
            (wind_speed_day_aug as fn() -> [f64; 24], DAYS_IN_MONTH[7]),
            (wind_speed_day_sep as fn() -> [f64; 24], DAYS_IN_MONTH[8]),
            (wind_speed_day_oct as fn() -> [f64; 24], DAYS_IN_MONTH[9]),
            (wind_speed_day_nov as fn() -> [f64; 24], DAYS_IN_MONTH[10]),
            (wind_speed_day_dec as fn() -> [f64; 24], DAYS_IN_MONTH[11]),
        ];
        for (speeds_fn, month_days_count) in months {
            speeds.extend_from_slice(
                speeds_fn()
                    .to_vec()
                    .iter()
                    .cloned()
                    .cycle()
                    .take((month_days_count * HOURS_IN_DAY) as usize)
                    .collect::<Vec<f64>>()
                    .as_slice(),
            );
        }

        speeds
    }

    const BASE_WIND_DIRECTIONS: [f64; 24] = [
        300., 250., 220., 180., 150., 120., 100., 80., 60., 40., 20., 10., 50., 100., 140., 190.,
        200., 320., 330., 340., 350., 355., 315., 5.,
    ];

    #[fixture]
    fn wind_directions() -> Vec<f64> {
        let wind_direction_day_jan = BASE_WIND_DIRECTIONS;
        let wind_direction_day_feb = BASE_WIND_DIRECTIONS.map(|d| d - 1.);
        let wind_direction_day_mar = BASE_WIND_DIRECTIONS.map(|d| d - 2.);
        let wind_direction_day_apr = BASE_WIND_DIRECTIONS.map(|d| d - 3.);
        let wind_direction_day_may = BASE_WIND_DIRECTIONS.map(|d| d - 4.);
        let wind_direction_day_jun = BASE_WIND_DIRECTIONS.map(|d| d + 1.);
        let wind_direction_day_jul = BASE_WIND_DIRECTIONS.map(|d| d + 2.);
        let wind_direction_day_aug = BASE_WIND_DIRECTIONS.map(|d| d + 3.);
        let wind_direction_day_sep = BASE_WIND_DIRECTIONS.map(|d| d + 4.);
        let wind_direction_day_oct = BASE_WIND_DIRECTIONS.map(|d| d - 5.);
        let wind_direction_day_nov = BASE_WIND_DIRECTIONS.map(|d| d + 5.);
        let wind_direction_day_dec = BASE_WIND_DIRECTIONS.map(|d| d - 0.);

        let mut wind_directions = Vec::with_capacity(8760);
        for (directions, days_in_month) in [
            (wind_direction_day_jan, DAYS_IN_MONTH[0]),
            (wind_direction_day_feb, DAYS_IN_MONTH[1]),
            (wind_direction_day_mar, DAYS_IN_MONTH[2]),
            (wind_direction_day_apr, DAYS_IN_MONTH[3]),
            (wind_direction_day_may, DAYS_IN_MONTH[4]),
            (wind_direction_day_jun, DAYS_IN_MONTH[5]),
            (wind_direction_day_jul, DAYS_IN_MONTH[6]),
            (wind_direction_day_aug, DAYS_IN_MONTH[7]),
            (wind_direction_day_sep, DAYS_IN_MONTH[8]),
            (wind_direction_day_oct, DAYS_IN_MONTH[9]),
            (wind_direction_day_nov, DAYS_IN_MONTH[10]),
            (wind_direction_day_dec, DAYS_IN_MONTH[11]),
        ] {
            wind_directions.extend_from_slice(
                directions
                    .iter()
                    .cloned()
                    .cycle()
                    .take((days_in_month * HOURS_IN_DAY) as usize)
                    .collect::<Vec<f64>>()
                    .as_slice(),
            );
        }

        wind_directions
    }

    #[fixture]
    fn diffuse_horizontal_radiation() -> [f64; 24] {
        [
            0., 0., 0., 0., 0., 0., 0., 0., 136., 308., 365., 300., 128., 90., 30., 0., 0., 0., 0.,
            0., 0., 0., 0., 0.,
        ]
    }

    #[fixture]
    fn direct_beam_radiation() -> [f64; 24] {
        [
            0., 0., 0., 0., 0., 0., 0., 0., 54., 113., 148., 149., 98., 50., 10., 0., 0., 0., 0.,
            0., 0., 0., 0., 0.,
        ]
    }

    #[fixture]
    fn solar_reflectivity_of_ground() -> [f64; 8760] {
        [0.2; 8760]
    }

    #[fixture]
    fn latitude() -> f64 {
        51.42
    }

    #[fixture]
    fn longitude() -> f64 {
        -0.75
    }

    #[fixture]
    fn timezone() -> i32 {
        0
    }

    #[fixture]
    fn start_day() -> u32 {
        0
    }

    #[fixture]
    fn end_day() -> Option<u32> {
        Some(0)
    }

    #[fixture]
    fn time_series_step() -> f64 {
        1.0
    }

    #[fixture]
    fn january_first() -> Option<u32> {
        Some(1)
    }

    #[fixture]
    fn daylight_savings() -> Option<DaylightSavingsConfig> {
        Some(NotApplicable)
    }

    #[fixture]
    fn leap_day_included() -> bool {
        false
    }

    #[fixture]
    fn direct_beam_conversion_needed() -> bool {
        false
    }

    #[fixture]
    fn shading_segments() -> Option<Vec<ShadingSegment>> {
        vec![
            ShadingSegment {
                start: 180.,
                end: 135.,
                ..Default::default()
            },
            ShadingSegment {
                start: 135.,
                end: 90.,
                ..Default::default()
            },
            ShadingSegment {
                start: 90.,
                end: 45.,
                ..Default::default()
            },
            ShadingSegment {
                start: 45.,
                end: 0.,
                shading_objects: Some(vec![ShadingObject {
                    object_type: ShadingObjectType::Obstacle,
                    height: 10.5,
                    distance: 12.,
                }]),
            },
            ShadingSegment {
                start: 0.,
                end: -45.,
                ..Default::default()
            },
            ShadingSegment {
                start: -45.,
                end: -90.,
                ..Default::default()
            },
            ShadingSegment {
                start: -90.,
                end: -135.,
                ..Default::default()
            },
            ShadingSegment {
                start: -135.,
                end: -180.,
                ..Default::default()
            },
        ]
        .into()
    }

    #[fixture]
    fn external_conditions() -> ExternalConditions {
        ExternalConditions::new(
            &simulation_time().iter(),
            air_temps(),
            wind_speeds(),
            wind_directions(),
            diffuse_horizontal_radiation().to_vec(),
            direct_beam_radiation().to_vec(),
            solar_reflectivity_of_ground().to_vec(),
            latitude(),
            longitude(),
            timezone(),
            start_day(),
            end_day(),
            time_series_step(),
            january_first(),
            daylight_savings(),
            leap_day_included(),
            direct_beam_conversion_needed(),
            shading_segments(),
        )
    }

    #[rstest]
    fn test_air_temp(external_conditions: ExternalConditions, simulation_time: SimulationTime) {
        for (i, simtime_step) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions.air_temp(&simtime_step),
                [3.5, 4.0, 4.5, 5.0, 7.5, 10.0, 12.5, 15.0][i],
                "failed on iteration index {} with step {:?}",
                i,
                simtime_step
            );
        }
    }

    #[rstest]
    fn test_air_temp_annual(external_conditions: ExternalConditions) {
        let precision = 1e-6;
        assert_relative_eq!(
            external_conditions.air_temp_annual().unwrap(),
            10.1801369863014,
            max_relative = precision
        );
    }

    #[rstest]
    fn test_air_temp_monthly(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        let expected_monthly_air_temps: [f64; 12] = [
            6.75, 7.75, 8.75, 9.75, 10.75, 11.75, 12.75, 12.75, 11.75, 10.75, 9.75, 8.75,
        ];
        let external_conditions = external_conditions.clone();
        for simtime_step in simulation_time.iter() {
            let month_idx = simtime_step.current_month().unwrap() as usize;
            assert_eq!(
                external_conditions.air_temp_monthly(simtime_step.current_month_start_end_hours()),
                expected_monthly_air_temps[month_idx]
            );
        }
    }

    #[rstest]
    fn test_wind_speed(external_conditions: ExternalConditions, simulation_time: SimulationTime) {
        let external_conditions = external_conditions.clone();
        for (i, simtime_step) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions.wind_speed(&simtime_step),
                [5.4, 5.7, 5.4, 5.6, 5.3, 5.1, 4.8, 4.7][i]
            );
        }
    }

    #[rstest]
    fn test_wind_speed_annual(external_conditions: ExternalConditions) {
        assert_relative_eq!(
            external_conditions.wind_speed_annual().unwrap(),
            4.23,
            max_relative = 0.01
        );
    }

    #[rstest]
    fn test_wind_direction(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        for t_it in simulation_time.iter() {
            assert_eq!(
                external_conditions.wind_direction(t_it),
                [80., 60., 40., 20., 10., 50., 100., 140.][t_it.index]
            );
        }
    }

    #[rstest]
    fn test_diffuse_horizontal_radiation(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        let external_conditions = external_conditions.clone();
        for (i, _simtime_step) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions.diffuse_horizontal_radiation(i),
                [0., 136., 308., 365., 300., 128., 90., 30.][i]
            );
        }
    }

    #[rstest]
    fn test_direct_beam_radiation(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        let external_conditions = external_conditions.clone();
        for (i, _simtime_step) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions.direct_beam_radiation(i),
                [0., 54., 113., 148., 149., 98., 50., 10.][i]
            );
        }
    }

    #[rstest]
    fn should_have_correct_solar_reflectivity_of_ground(
        external_conditions: ExternalConditions,
        solar_reflectivity_of_ground: [f64; 8760],
        simulation_time: SimulationTime,
    ) {
        let external_conditions = external_conditions.clone();
        for (i, simtime_step) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions.solar_reflectivity_of_ground(&simtime_step),
                solar_reflectivity_of_ground[i]
            );
        }
    }

    /// Test that using a reveal produces the same results as an equivalent overhang and side fins
    #[rstest]
    fn test_window_shading(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        let reveal_depth = 0.1;
        let reveal_distance = 0.2;
        let base_height = 0.;
        let height = 2.;
        let width = 2.;
        let tilt = 90.;
        let orientation = 180.;

        // Create shading objects with reveal
        let shading_with_reveal = vec![WindowShadingObject::Reveal {
            depth: reveal_depth,
            distance: reveal_distance,
        }];

        // Create shading objects with overhang and fins
        let shading_with_overhang_fin = vec![
            WindowShadingObject::Overhang {
                depth: reveal_depth,
                distance: reveal_distance,
            },
            WindowShadingObject::SideFinRight {
                depth: reveal_depth,
                distance: reveal_distance,
            },
            WindowShadingObject::SideFinLeft {
                depth: reveal_depth,
                distance: reveal_distance,
            },
        ];

        for t_it in simulation_time.iter() {
            // Calculate shading factors with reveal
            let shading_factor_reveal = external_conditions
                .shading_reduction_factor_direct_diffuse(
                    base_height,
                    height,
                    width,
                    tilt,
                    orientation,
                    &shading_with_reveal,
                    t_it,
                )
                .unwrap();

            // Calculate shading factors with overhang and fins
            let shading_factor_overhang_fin = external_conditions
                .shading_reduction_factor_direct_diffuse(
                    base_height,
                    height,
                    width,
                    tilt,
                    orientation,
                    &shading_with_overhang_fin,
                    t_it,
                )
                .unwrap();

            assert_relative_eq!(
                shading_factor_reveal.0,
                shading_factor_overhang_fin.0,
                max_relative = 1e-5
            );
            assert_relative_eq!(
                shading_factor_reveal.1,
                shading_factor_overhang_fin.1,
                max_relative = 1e-5
            );
        }
    }

    #[rstest]
    fn test_init_direct_beam_radiation(direct_beam_conversion_needed: bool) {
        // Check with direct_beam_conversion_needed = false
        assert_eq!(
            init_direct_beam_radiation(direct_beam_conversion_needed, 100., 5.0),
            100.0
        );

        // Check with direct_beam_conversion_needed = true
        assert_relative_eq!(
            init_direct_beam_radiation(true, 100., 5.0),
            1147.3713245669855,
            max_relative = 1e-8
        );

        // Check with solar altitude = 0
        assert_eq!(init_direct_beam_radiation(true, 10., 0.0), 10.);
    }

    #[rstest]
    fn test_init_earth_orbit_deviation() {
        // check for non-leap year
        assert_eq!(init_earth_orbit_deviation(364), 360.0);

        // check for leap year
        assert_relative_eq!(
            init_earth_orbit_deviation(365),
            360.986301369863,
            max_relative = 1e-10
        );
    }

    #[rstest]
    fn test_init_solar_declination() {
        assert_relative_eq!(
            init_solar_declination(100.),
            8.239299094353976,
            max_relative = 1e-8
        );
    }

    #[rstest]
    fn test_init_equation_of_time() {
        assert_eq!(init_equation_of_time(1), 3.48);
        assert_relative_eq!(
            init_equation_of_time(22),
            12.001736328751521,
            max_relative = 1e-8
        );
        assert_relative_eq!(
            init_equation_of_time(137),
            -3.5547083185334247,
            max_relative = 1e-8
        );
        assert_relative_eq!(
            init_equation_of_time(242),
            0.12076422612546622,
            max_relative = 1e-8
        );
        assert_eq!(init_equation_of_time(365), 3.15);
    }

    #[rstest]
    #[should_panic]
    fn test_init_equation_of_time_panics_with_out_of_bounds_param() {
        init_equation_of_time(366);
    }

    #[rstest]
    fn test_init_time_shift(external_conditions: ExternalConditions) {
        assert_eq!(
            init_time_shift(external_conditions.timezone, external_conditions.longitude),
            0.05
        );

        assert_relative_eq!(
            init_time_shift(-5, -73.),
            -0.13333333333333375,
            max_relative = 1e-8
        );
    }

    #[rstest]
    fn test_init_solar_time() {
        // Mid day without timeshift
        assert_relative_eq!(
            init_solar_time(12, 5., 0., 0),
            12.916666666666666,
            max_relative = 1e-8
        );
        // Mid day with timeshift
        assert_relative_eq!(
            init_solar_time(12, 5., 1., 0),
            11.916666666666666,
            max_relative = 1e-8
        );
        // End of day with -ve equation and +ve timeshift
        assert_relative_eq!(
            init_solar_time(23, -2.5, 1., 0),
            23.041666666666668,
            max_relative = 1e-8
        );
    }

    #[rstest]
    fn test_init_solar_hour_angle() {
        assert_eq!(init_solar_hour_angle(23.), -157.5);
        assert_eq!(init_solar_hour_angle(1.), 172.5);
    }

    #[rstest]
    fn test_init_solar_altitude(external_conditions: ExternalConditions) {
        assert_relative_eq!(
            init_solar_altitude(external_conditions.latitude, 10., 50.),
            32.03953794285834,
            max_relative = 1e-8
        );
        assert_eq!(
            init_solar_altitude(external_conditions.latitude, -23., 60.),
            0.0
        );
    }

    #[rstest]
    fn test_init_solar_azimuth_angle(external_conditions: ExternalConditions) {
        let lat = external_conditions.latitude;

        assert_relative_eq!(
            init_solar_azimuth_angle(lat, 23., 0., 45.),
            9.134285141104091e-15,
            max_relative = 1e-24
        );

        // check for East
        assert_relative_eq!(
            init_solar_azimuth_angle(lat, 23., -15., 45.),
            -19.49201220785029,
            max_relative = 1e-8
        );

        // check for West
        assert_relative_eq!(
            init_solar_azimuth_angle(lat, 23., 15., 45.),
            19.49201220785031,
            max_relative = 1e-8
        );

        // Negative declination
        assert_relative_eq!(
            init_solar_azimuth_angle(lat, -23., 15., 45.),
            19.49201220785031,
            max_relative = 1e-8
        );
    }

    #[rstest]
    fn test_init_air_mass() {
        assert_relative_eq!(init_air_mass(5.), 10.323080326274896, max_relative = 1e-8);
        assert_relative_eq!(init_air_mass(15.), 3.8637033051562737, max_relative = 1e-8);
    }

    #[rstest]
    fn test_solar_angle_of_incidence(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.solar_angle_of_incidence(10., 10., &t_it),
                [
                    89.28367858027447,
                    80.39193381264141,
                    73.0084830472421,
                    67.6707156423694,
                    64.91086719421193,
                    65.0693571287282,
                    68.1252145585689,
                    73.70764473968616
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.solar_angle_of_incidence(0., 10., &t_it),
                [
                    95.78366411883604,
                    88.23114711240953,
                    81.97931494689718,
                    77.41946173012103,
                    74.90762116550648,
                    74.67327369297139,
                    76.73922946774607,
                    80.91274624480684
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.solar_angle_of_incidence(90., -180., &t_it),
                [
                    120.13031074122472,
                    131.83510459862302,
                    143.4374125410406,
                    154.33449512110948,
                    162.6874363991092,
                    163.66695597946457,
                    156.3258508149546,
                    145.71870170993543
                ][t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    fn test_sun_surface_azimuth(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.sun_surface_azimuth(180., t_it),
                [
                    -110.99000000000001,
                    -125.99,
                    -140.99,
                    -155.99,
                    -170.98999999999998,
                    174.01000000000002,
                    159.01,
                    144.01
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.sun_surface_azimuth(0., t_it),
                [
                    69.00999999999999,
                    54.010000000000005,
                    39.010000000000005,
                    24.010000000000005,
                    9.010000000000007,
                    -5.989999999999993,
                    -20.989999999999995,
                    -35.989999999999995
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.sun_surface_azimuth(-180., t_it),
                [
                    -110.99000000000001,
                    -125.99000000000001,
                    -140.99,
                    -155.99,
                    -170.98999999999998,
                    174.01000000000002,
                    159.01,
                    144.01
                ][t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    fn test_sun_surface_tilt(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.sun_surface_tilt(0., t_it),
                [
                    -90.,
                    -88.23114711240953,
                    -81.97931494689718,
                    -77.41946173012103,
                    -74.90762116550647,
                    -74.67327369297139,
                    -76.73922946774607,
                    -80.91274624480684
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.sun_surface_tilt(90., t_it),
                [
                    0.,
                    1.7688528875904694,
                    8.020685053102824,
                    12.580538269878971,
                    15.09237883449353,
                    15.326726307028608,
                    13.260770532253929,
                    9.08725375519316
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.sun_surface_tilt(180., t_it),
                [
                    90.,
                    91.76885288759047,
                    98.02068505310282,
                    102.58053826987897,
                    105.09237883449353,
                    105.32672630702861,
                    103.26077053225393,
                    99.08725375519316
                ][t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    fn test_direct_irradiance(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.direct_irradiance(0., 180., &t_it),
                [
                    0.0,
                    1.6668397643248698,
                    15.766957881489741,
                    32.23613721421553,
                    38.79603659734242,
                    25.90364916630168,
                    11.469168196073362,
                    1.579383993776423
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions.direct_irradiance(65., 180., &t_it),
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0][t_idx],
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions.direct_irradiance(65., -180., &t_it),
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0][t_idx],
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.direct_irradiance(65., 0., &t_it),
                [
                    0.,
                    33.34729692480562,
                    88.92202737221507,
                    134.5232406821787,
                    145.31786982869673,
                    96.18110717866088,
                    46.34890065056773,
                    8.15613633989895
                ][t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    fn test_init_extra_terrestrial_radiation() {
        assert_relative_eq!(
            init_extra_terrestrial_radiation(0.),
            1412.1109999999999,
            max_relative = 1e-8
        );

        assert_eq!(init_extra_terrestrial_radiation(90.), 1367.0);

        assert_eq!(init_extra_terrestrial_radiation(180.), 1321.889);
    }

    #[rstest]
    fn test_brightness_coefficient() {
        assert_eq!(
            brightness_coefficient(1., BrightnessCoefficientName::F11),
            -0.008
        );
        assert_eq!(
            brightness_coefficient(1.23, BrightnessCoefficientName::F12),
            0.487
        );
        assert_eq!(
            brightness_coefficient(1.6, BrightnessCoefficientName::F13),
            -0.295
        );
        assert_eq!(
            brightness_coefficient(2., BrightnessCoefficientName::F21),
            0.226
        );
        assert_eq!(
            brightness_coefficient(2.9, BrightnessCoefficientName::F12),
            -1.237
        );
        assert_eq!(
            brightness_coefficient(4.5, BrightnessCoefficientName::F22),
            -1.127
        );
        assert_eq!(
            brightness_coefficient(7., BrightnessCoefficientName::F23),
            0.251
        );
    }

    #[rstest]
    fn test_init_f1() {
        assert_relative_eq!(
            init_f1_circumsolar_brightness_coefficient(1.23, 1., 30.),
            0.7012846705927759,
            max_relative = 1e-8
        );
        assert_eq!(init_f1_circumsolar_brightness_coefficient(3., 1., 30.), 0.0);
    }

    #[rstest]
    fn test_init_f2() {
        assert_relative_eq!(
            init_f2_horizontal_brightness_coefficient(1.23, 1., 180.),
            -0.09068140899333463,
            max_relative = 1e-8
        );
        assert_relative_eq!(
            init_f2_horizontal_brightness_coefficient(3., 0., 30.),
            0.3173215314335047,
            max_relative = 1e-8
        );
    }

    #[rstest]
    fn test_init_dimensionless_clearness_parameter() {
        assert_eq!(init_dimensionless_clearness_parameter(0., 5., 15.), 999.);
        assert_relative_eq!(
            init_dimensionless_clearness_parameter(1., 5., 15.),
            5.910652372233723,
            max_relative = 1e-8
        );
    }

    #[rstest]
    fn test_init_dimensionless_sky_brightness_parameter() {
        assert_relative_eq!(
            init_dimensionless_sky_brightness_parameter(36.0, 50., 1412.),
            1.274787535410765,
            max_relative = 1e-8
        );
    }

    #[rstest]
    fn test_a_over_b(external_conditions: ExternalConditions, simulation_time: SimulationTime) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.a_over_b(0., 0., &t_it),
                [
                    0.0,
                    0.354163731154509,
                    1.0,
                    1.0,
                    0.9999999999999986,
                    1.0,
                    1.0,
                    1.0
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions.a_over_b(90., 180., &t_it),
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,][t_idx],
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.a_over_b(180., -180., &t_it),
                [1.156236348588363, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0][t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    fn test_diffuse_irradiance(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.diffuse_irradiance(0., 0., &t_it),
                [
                    DiffuseIrradiance(0.0, 0.0, 0.0, 0.0),
                    DiffuseIrradiance(
                        51.050674242609674,
                        4.466157979564503,
                        46.584516263045174,
                        -0.0
                    ),
                    DiffuseIrradiance(308.0, 80.07131055825931, 227.9286894417407, -0.0),
                    DiffuseIrradiance(365.0, 142.60279131112304, 222.39720868887696, -0.0),
                    DiffuseIrradiance(
                        299.9999999999998,
                        168.47209561804453,
                        131.52790438195527,
                        -0.0
                    ),
                    DiffuseIrradiance(128.0, 96.29997383163165, 31.700026168368353, 0.0),
                    DiffuseIrradiance(90.0, 69.76354921067887, 20.236450789321122, 0.0),
                    DiffuseIrradiance(30.0, 27.570058662939754, 2.4299413370602463, 0.0)
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.diffuse_irradiance(90., 180., &t_it),
                [
                    DiffuseIrradiance(0.0, 0.0, 0.0, 0.0),
                    DiffuseIrradiance(
                        -13.202357427309334,
                        2.2330789897822516,
                        0.0,
                        -15.435436417091585
                    ),
                    DiffuseIrradiance(
                        16.122288185513806,
                        40.03565527912966,
                        0.0,
                        -23.91336709361585
                    ),
                    DiffuseIrradiance(
                        50.83171656478211,
                        71.30139565556152,
                        0.0,
                        -20.469679090779405
                    ),
                    DiffuseIrradiance(
                        74.87257439807136,
                        84.23604780902227,
                        0.0,
                        -9.363473410950908
                    ),
                    DiffuseIrradiance(53.09439354327976, 48.14998691581582, 0.0, 4.944406627463934),
                    DiffuseIrradiance(39.20317296034313, 34.88177460533944, 0.0, 4.321398355003694),
                    DiffuseIrradiance(
                        14.084774138824505,
                        13.785029331469877,
                        0.0,
                        0.29974480735462805
                    )
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.diffuse_irradiance(180., -180., &t_it),
                [
                    DiffuseIrradiance(0.0, 0.0, 0.0, 0.0),
                    DiffuseIrradiance(-1.89029578016337e-15, 0.0, 0.0, -1.89029578016337e-15),
                    DiffuseIrradiance(-2.9285428468032295e-15, 0.0, 0.0, -2.9285428468032295e-15),
                    DiffuseIrradiance(-2.50681269780965e-15, 0.0, 0.0, -2.50681269780965e-15),
                    DiffuseIrradiance(-1.1466947741622378e-15, 0.0, 0.0, -1.1466947741622378e-15),
                    DiffuseIrradiance(6.055151750006666e-16, 0.0, 0.0, 6.055151750006666e-16),
                    DiffuseIrradiance(5.292186663295911e-16, 0.0, 0.0, 5.292186663295911e-16),
                    DiffuseIrradiance(3.670815188878853e-17, 0.0, 0.0, 3.670815188878853e-17)
                ][t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    fn test_ground_reflection_irradiance(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions.ground_reflection_irradiance(0., &t_it),
                [0., 0., 0., 0., 0., 0., 0., 0.][t_idx]
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.ground_reflection_irradiance(90., &t_it),
                [
                    0.0,
                    13.766683976432486,
                    32.37669578814897,
                    39.72361372142155,
                    33.87960365973424,
                    15.390364916630167,
                    10.146916819607336,
                    3.157938399377642
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.ground_reflection_irradiance(180., &t_it),
                [
                    0.0,
                    27.533367952864975,
                    64.75339157629796,
                    79.44722744284311,
                    67.75920731946849,
                    30.780729833260338,
                    20.293833639214675,
                    6.315876798755285
                ][t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    fn test_circumsolar_irradiance(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.circumsolar_irradiance(0., 0., &t_it),
                [
                    0.0,
                    46.584516263045174,
                    227.9286894417407,
                    222.39720868887696,
                    131.52790438195527,
                    31.700026168368353,
                    20.236450789321122,
                    2.4299413370602463
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions.circumsolar_irradiance(90., 180., &t_it),
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,][t_idx],
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions.circumsolar_irradiance(180., -180., &t_it),
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,][t_idx],
            );
        }
    }

    #[rstest]
    fn test_calculated_direct_irradiance(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.calculated_direct_irradiance(0., 0., &t_it),
                [
                    0.0,
                    48.25135602737004,
                    243.69564732323042,
                    254.6333459030925,
                    170.32394097929767,
                    57.60367533467003,
                    31.705618985394484,
                    4.009325330836669
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions.calculated_direct_irradiance(90., 180., &t_it),
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,][t_idx],
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions.calculated_direct_irradiance(180., -180., &t_it),
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,][t_idx],
            );
        }
    }

    #[rstest]
    fn test_calculated_diffuse_irradiance(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.calculated_diffuse_irradiance(0., 0., &t_it),
                [
                    0.0,
                    4.4661579795645,
                    80.07131055825931,
                    142.60279131112304,
                    168.4720956180445,
                    96.29997383163165,
                    69.76354921067887,
                    27.570058662939754
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.calculated_diffuse_irradiance(90., 180., &t_it),
                [
                    0.0,
                    0.5643265491231517,
                    48.498983973662774,
                    90.55533028620366,
                    108.7521780578056,
                    68.48475845990993,
                    49.350089779950466,
                    17.24271253820215
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.calculated_diffuse_irradiance(180., -180., &t_it),
                [
                    0.0,
                    27.53336795286497,
                    64.75339157629796,
                    79.44722744284311,
                    67.75920731946849,
                    30.780729833260338,
                    20.293833639214675,
                    6.315876798755285
                ][t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    fn test_calculated_total_solar_irradiance(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.calculated_total_solar_irradiance(0., 0., &t_it),
                [
                    0.0,
                    52.71751400693454,
                    323.7669578814897,
                    397.23613721421555,
                    338.7960365973422,
                    153.90364916630168,
                    101.46916819607335,
                    31.579383993776425
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.calculated_total_solar_irradiance(90., 180., &t_it),
                [
                    0.0,
                    0.5643265491231517,
                    48.498983973662774,
                    90.55533028620366,
                    108.7521780578056,
                    68.48475845990993,
                    49.350089779950466,
                    17.24271253820215
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.calculated_total_solar_irradiance(180., -180., &t_it),
                [
                    0.0,
                    27.53336795286497,
                    64.75339157629796,
                    79.44722744284311,
                    67.75920731946849,
                    30.780729833260338,
                    20.293833639214675,
                    6.315876798755285
                ][t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    fn test_calculated_direct_diffuse_total_irradiance(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.calculated_direct_diffuse_total_irradiance(0., 0., true, &t_it),
                [
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        0.0,
                        0.0,
                        Some(DiffuseBreakdown {
                            sky: 0.0,
                            circumsolar: 0.0,
                            horiz: 0.0,
                            ground_refl: 0.0
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        48.25135602737004,
                        4.4661579795645,
                        52.71751400693454,
                        Some(DiffuseBreakdown {
                            sky: 4.466157979564503,
                            circumsolar: 46.584516263045174,
                            horiz: -0.0,
                            ground_refl: 0.0
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        243.69564732323042,
                        80.07131055825931,
                        323.7669578814897,
                        Some(DiffuseBreakdown {
                            sky: 80.07131055825931,
                            circumsolar: 227.9286894417407,
                            horiz: -0.0,
                            ground_refl: 0.0
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        254.6333459030925,
                        142.60279131112304,
                        397.23613721421555,
                        Some(DiffuseBreakdown {
                            sky: 142.60279131112304,
                            circumsolar: 222.39720868887696,
                            horiz: -0.0,
                            ground_refl: 0.0
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        170.32394097929767,
                        168.4720956180445,
                        338.7960365973422,
                        Some(DiffuseBreakdown {
                            sky: 168.47209561804453,
                            circumsolar: 131.52790438195527,
                            horiz: -0.0,
                            ground_refl: 0.0
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        57.60367533467003,
                        96.29997383163165,
                        153.90364916630168,
                        Some(DiffuseBreakdown {
                            sky: 96.29997383163165,
                            circumsolar: 31.700026168368353,
                            horiz: 0.0,
                            ground_refl: 0.0
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        31.705618985394484,
                        69.76354921067887,
                        101.46916819607335,
                        Some(DiffuseBreakdown {
                            sky: 69.76354921067887,
                            circumsolar: 20.236450789321122,
                            horiz: 0.0,
                            ground_refl: 0.0
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        4.009325330836669,
                        27.570058662939754,
                        31.579383993776425,
                        Some(DiffuseBreakdown {
                            sky: 27.570058662939754,
                            circumsolar: 2.4299413370602463,
                            horiz: 0.0,
                            ground_refl: 0.0
                        })
                    ),
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions
                    .calculated_direct_diffuse_total_irradiance(90., 180., true, &t_it),
                [
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        0.0,
                        0.0,
                        Some(DiffuseBreakdown {
                            sky: 0.0,
                            circumsolar: 0.0,
                            horiz: 0.0,
                            ground_refl: 0.0
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        0.5643265491231517,
                        0.5643265491231517,
                        Some(DiffuseBreakdown {
                            sky: 2.2330789897822516,
                            circumsolar: 0.0,
                            horiz: -15.435436417091585,
                            ground_refl: 13.766683976432486
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        48.498983973662774,
                        48.498983973662774,
                        Some(DiffuseBreakdown {
                            sky: 40.03565527912966,
                            circumsolar: 0.0,
                            horiz: -23.91336709361585,
                            ground_refl: 32.37669578814897
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        90.55533028620366,
                        90.55533028620366,
                        Some(DiffuseBreakdown {
                            sky: 71.30139565556152,
                            circumsolar: 0.0,
                            horiz: -20.469679090779405,
                            ground_refl: 39.72361372142155
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        108.7521780578056,
                        108.7521780578056,
                        Some(DiffuseBreakdown {
                            sky: 84.23604780902227,
                            circumsolar: 0.0,
                            horiz: -9.363473410950908,
                            ground_refl: 33.87960365973424
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        68.48475845990993,
                        68.48475845990993,
                        Some(DiffuseBreakdown {
                            sky: 48.14998691581582,
                            circumsolar: 0.0,
                            horiz: 4.944406627463934,
                            ground_refl: 15.390364916630167
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        49.350089779950466,
                        49.350089779950466,
                        Some(DiffuseBreakdown {
                            sky: 34.88177460533944,
                            circumsolar: 0.0,
                            horiz: 4.321398355003694,
                            ground_refl: 10.146916819607336
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        17.24271253820215,
                        17.24271253820215,
                        Some(DiffuseBreakdown {
                            sky: 13.785029331469877,
                            circumsolar: 0.0,
                            horiz: 0.29974480735462805,
                            ground_refl: 3.157938399377642
                        })
                    ),
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions
                    .calculated_direct_diffuse_total_irradiance(180., -180., true, &t_it),
                [
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        0.0,
                        0.0,
                        Some(DiffuseBreakdown {
                            sky: 0.0,
                            circumsolar: 0.0,
                            horiz: 0.0,
                            ground_refl: 0.0
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        27.53336795286497,
                        27.53336795286497,
                        Some(DiffuseBreakdown {
                            sky: 0.0,
                            circumsolar: 0.0,
                            horiz: -1.89029578016337e-15,
                            ground_refl: 27.533367952864975
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        64.75339157629796,
                        64.75339157629796,
                        Some(DiffuseBreakdown {
                            sky: 0.0,
                            circumsolar: 0.0,
                            horiz: -2.9285428468032295e-15,
                            ground_refl: 64.75339157629796
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        79.44722744284311,
                        79.44722744284311,
                        Some(DiffuseBreakdown {
                            sky: 0.0,
                            circumsolar: 0.0,
                            horiz: -2.50681269780965e-15,
                            ground_refl: 79.44722744284311
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        67.75920731946849,
                        67.75920731946849,
                        Some(DiffuseBreakdown {
                            sky: 0.0,
                            circumsolar: 0.0,
                            horiz: -1.1466947741622378e-15,
                            ground_refl: 67.75920731946849
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        30.780729833260338,
                        30.780729833260338,
                        Some(DiffuseBreakdown {
                            sky: 0.0,
                            circumsolar: 0.0,
                            horiz: 6.055151750006666e-16,
                            ground_refl: 30.780729833260338
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        20.293833639214675,
                        20.293833639214675,
                        Some(DiffuseBreakdown {
                            sky: 0.0,
                            circumsolar: 0.0,
                            horiz: 5.292186663295911e-16,
                            ground_refl: 20.293833639214675
                        })
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        6.315876798755285,
                        6.315876798755285,
                        Some(DiffuseBreakdown {
                            sky: 0.0,
                            circumsolar: 0.0,
                            horiz: 3.670815188878853e-17,
                            ground_refl: 6.315876798755285
                        })
                    ),
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions
                    .calculated_direct_diffuse_total_irradiance(0., 0., false, &t_it),
                [
                    CalculatedDirectDiffuseTotalIrradiance(0.0, 0.0, 0.0, None),
                    CalculatedDirectDiffuseTotalIrradiance(
                        48.25135602737004,
                        4.4661579795645,
                        52.71751400693454,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        243.69564732323042,
                        80.07131055825931,
                        323.7669578814897,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        254.6333459030925,
                        142.60279131112304,
                        397.23613721421555,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        170.32394097929767,
                        168.4720956180445,
                        338.7960365973422,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        57.60367533467003,
                        96.29997383163165,
                        153.90364916630168,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        31.705618985394484,
                        69.76354921067887,
                        101.46916819607335,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        4.009325330836669,
                        27.570058662939754,
                        31.579383993776425,
                        None
                    ),
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions
                    .calculated_direct_diffuse_total_irradiance(90., 180., false, &t_it),
                [
                    CalculatedDirectDiffuseTotalIrradiance(0.0, 0.0, 0.0, None),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        0.5643265491231517,
                        0.5643265491231517,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        48.498983973662774,
                        48.498983973662774,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        90.55533028620366,
                        90.55533028620366,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        108.7521780578056,
                        108.7521780578056,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        68.48475845990993,
                        68.48475845990993,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        49.350089779950466,
                        49.350089779950466,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        17.24271253820215,
                        17.24271253820215,
                        None
                    ),
                ][t_idx],
                max_relative = 1e-8
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions
                    .calculated_direct_diffuse_total_irradiance(180., -180., false, &t_it),
                [
                    CalculatedDirectDiffuseTotalIrradiance(0.0, 0.0, 0.0, None),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        27.53336795286497,
                        27.53336795286497,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        64.75339157629796,
                        64.75339157629796,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        79.44722744284311,
                        79.44722744284311,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        67.75920731946849,
                        67.75920731946849,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        30.780729833260338,
                        30.780729833260338,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        20.293833639214675,
                        20.293833639214675,
                        None
                    ),
                    CalculatedDirectDiffuseTotalIrradiance(
                        0.0,
                        6.315876798755285,
                        6.315876798755285,
                        None
                    ),
                ][t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    fn test_outside_solar_beam(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        for t_it in simulation_time.iter() {
            assert!(!external_conditions.outside_solar_beam(0., 0., &t_it));
        }

        for t_it in simulation_time.iter() {
            assert!(external_conditions.outside_solar_beam(90., 180., &t_it));
        }

        for t_it in simulation_time.iter() {
            assert!(external_conditions.outside_solar_beam(180., -180., &t_it));
        }
    }

    #[rstest]
    fn test_get_segment(
        mut external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions.get_segment(&t_it).unwrap(),
                [
                    ShadingSegment {
                        start: 90.,
                        end: 45.,
                        ..Default::default()
                    },
                    ShadingSegment {
                        start: 90.,
                        end: 45.,
                        ..Default::default()
                    },
                    ShadingSegment {
                        start: 45.,
                        end: 0.,
                        shading_objects: Some(vec![ShadingObject {
                            object_type: ShadingObjectType::Obstacle,
                            height: 10.5,
                            distance: 12.
                        }]),
                    },
                    ShadingSegment {
                        start: 45.,
                        end: 0.,
                        shading_objects: Some(vec![ShadingObject {
                            object_type: ShadingObjectType::Obstacle,
                            height: 10.5,
                            distance: 12.
                        }]),
                    },
                    ShadingSegment {
                        start: 45.,
                        end: 0.,
                        shading_objects: Some(vec![ShadingObject {
                            object_type: ShadingObjectType::Obstacle,
                            height: 10.5,
                            distance: 12.
                        }]),
                    },
                    ShadingSegment {
                        start: 0.,
                        end: -45.,
                        ..Default::default()
                    },
                    ShadingSegment {
                        start: 0.,
                        end: -45.,
                        ..Default::default()
                    },
                    ShadingSegment {
                        start: 0.,
                        end: -45.,
                        ..Default::default()
                    },
                ][t_idx]
            );
        }

        // For the gap in second shading segment
        external_conditions.shading_segments = vec![
            ShadingSegment {
                start: 180.,
                end: 135.,
                ..Default::default()
            },
            ShadingSegment {
                start: 50.,
                end: 90.,
                ..Default::default()
            },
            ShadingSegment {
                start: 90.,
                end: 45.,
                ..Default::default()
            },
        ]
        .into();

        assert!(external_conditions
            .get_segment(&simulation_time.iter().next().unwrap())
            .is_err());

        // For the value of end > start in second shading segment
        external_conditions.shading_segments = vec![
            ShadingSegment {
                start: 180.,
                end: 135.,
                ..Default::default()
            },
            ShadingSegment {
                start: 135.,
                end: 140.,
                ..Default::default()
            },
            ShadingSegment {
                start: 90.,
                end: 45.,
                ..Default::default()
            },
        ]
        .into();

        assert!(external_conditions
            .get_segment(&simulation_time.iter().next().unwrap())
            .is_err());
    }

    #[rstest]
    fn test_obstacle_shading_height(
        mut external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        // Test for peak hours in Jan
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.obstacle_shading_height(3., 5., 1., &t_it),
                [
                    2.0,
                    1.9691178812624324,
                    1.8590909935052966,
                    1.7768301326744398,
                    1.7303219853707263,
                    1.7259295208086824,
                    1.7643328437619492,
                    1.8400541145006308
                ][t_idx],
                max_relative = 1e-8
            );
        }

        // Test for peak hours in July
        let simulation_time_july = SimulationTime::new(4474., 4482., 1.);
        let base_zero_vec = vec![0.; 4474];
        let mut direct_beam_radiations = base_zero_vec.clone();
        direct_beam_radiations.extend_from_slice(&[70., 71., 39., 51., 44., 26., 108., 141.]);
        external_conditions.direct_beam_radiations = direct_beam_radiations;
        let mut diffuse_horizontal_radiations = base_zero_vec.clone();
        diffuse_horizontal_radiations
            .extend_from_slice(&[232., 310., 342., 393., 421., 426., 466., 424., 397.]);
        external_conditions.diffuse_horizontal_radiations = diffuse_horizontal_radiations;
        let mut solar_reflectivity_of_ground = base_zero_vec.clone();
        solar_reflectivity_of_ground.extend_from_slice(&[0.6; 8]);
        external_conditions.solar_reflectivity_of_ground = solar_reflectivity_of_ground;

        for (t_idx, t_it) in simulation_time_july.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.obstacle_shading_height(3., 5., 1., &t_it),
                [
                    0.5367256982895989,
                    0.240221480766996,
                    0.1963213106984325,
                    0.4473885431039979,
                    0.79264418501541,
                    1.1031372689550807,
                    1.3579235500301707,
                    1.5678670755252178
                ][t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    fn test_overhang_shading_height(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions.overhang_shading_height(5., 1., 3., 2., t_it),
                [
                    3.0,
                    3.0617642374751353,
                    3.281818012989407,
                    3.4463397346511204,
                    3.5393560292585473,
                    3.548140958382635,
                    3.4713343124761016,
                    3.3198917709987383
                ][t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    fn test_direct_shading_reduction_factor(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        // obstacle shading defined in segment and empty window shading
        let base_height = 1.;
        let height = 1.25;
        let width = 4.;
        let orientation = 90.;
        let window_shading = [];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions
                    .direct_shading_reduction_factor(
                        base_height,
                        height,
                        width,
                        orientation,
                        Some(&window_shading),
                        t_it
                    )
                    .unwrap(),
                [1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0][t_idx]
            );
        }

        // with window shading
        let window_shading_val = vec![
            WindowShadingObject::Overhang {
                depth: 0.5,
                distance: 0.5,
            },
            WindowShadingObject::SideFinLeft {
                depth: 0.25,
                distance: 0.1,
            },
            WindowShadingObject::SideFinRight {
                depth: 0.25,
                distance: 0.1,
            },
        ];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions
                    .direct_shading_reduction_factor(
                        base_height,
                        height,
                        width,
                        orientation,
                        Some(&window_shading_val),
                        t_it
                    )
                    .unwrap(),
                [
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    0.4002331427001323,
                    0.8511094450438235,
                    0.9292880457188379
                ][t_idx],
                max_relative = 1e-8
            );
        }

        // Test with zero orientation and with window shading
        let orientation = 0.;
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions
                    .direct_shading_reduction_factor(
                        base_height,
                        height,
                        width,
                        orientation,
                        Some(&window_shading_val),
                        t_it
                    )
                    .unwrap(),
                [
                    0.9201388647583034,
                    0.9552620660549097,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0
                ][t_idx],
                max_relative = 1e-8
            );
        }

        // Test with negative orientation and with window shading
        let orientation = -180.;
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions
                    .direct_shading_reduction_factor(
                        base_height,
                        height,
                        width,
                        orientation,
                        Some(&window_shading_val),
                        t_it
                    )
                    .unwrap(),
                [
                    0.9201388647583035,
                    0.9552620660549097,
                    0.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0
                ][t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    fn test_diffuse_shading_reduction_factor(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        // With window shading
        let diffuse_breakdown = DiffuseBreakdown {
            sky: 12.0,
            circumsolar: 0.0,
            horiz: -2.0164780331874668,
            ground_refl: 2.4,
        };
        let tilt = 90.;
        let height = 1.25;
        let base_height = 1.;
        let width = 4.;
        let orientation = 90.;
        let window_shading = vec![
            WindowShadingObject::Overhang {
                depth: 0.5,
                distance: 0.5,
            },
            WindowShadingObject::SideFinLeft {
                depth: 0.25,
                distance: 0.1,
            },
            WindowShadingObject::SideFinRight {
                depth: 0.25,
                distance: 0.1,
            },
        ];
        let f_sky = 0.5;

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions
                    .diffuse_shading_reduction_factor(
                        diffuse_breakdown,
                        tilt,
                        height,
                        base_height,
                        width,
                        orientation,
                        Some(&window_shading),
                        f_sky,
                        t_it
                    )
                    .unwrap(),
                [
                    0.5905238865770632,
                    0.5905238865770632,
                    0.5905238865770632,
                    0.5905238865770632,
                    0.5905238865770632,
                    0.5905238865770632,
                    0.5905238865770632,
                    0.5905238865770632
                ][t_idx],
                max_relative = 1e-8
            );
        }

        // With no window shading
        let window_shading = vec![];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions
                    .diffuse_shading_reduction_factor(
                        diffuse_breakdown,
                        tilt,
                        height,
                        base_height,
                        width,
                        orientation,
                        Some(&window_shading),
                        f_sky,
                        t_it
                    )
                    .unwrap(),
                [
                    0.9699932956653645,
                    0.9699932956653645,
                    0.9699932956653645,
                    0.9699932956653645,
                    0.9699932956653645,
                    0.9699932956653645,
                    0.9699932956653645,
                    0.9699932956653645
                ][t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    fn test_solar_irradiance(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        let base_height = 1.;
        let height = 1.25;
        let width = 4.;
        let tilt = 0.;
        let orientation = 90.;
        let window_shading = vec![
            WindowShadingObject::Overhang {
                depth: 0.5,
                distance: 0.5,
            },
            WindowShadingObject::SideFinLeft {
                depth: 0.25,
                distance: 0.1,
            },
            WindowShadingObject::SideFinRight {
                depth: 0.25,
                distance: 0.1,
            },
        ];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                external_conditions
                    .surface_irradiance(
                        base_height,
                        height,
                        width,
                        tilt,
                        orientation,
                        &window_shading,
                        t_it
                    )
                    .unwrap(),
                [
                    0.0,
                    50.88872899552963,
                    47.28402151418233,
                    84.21035456178225,
                    99.48679668415026,
                    114.47111015899463,
                    72.90266120669479,
                    20.290103525633484
                ][t_idx],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    fn test_shading_reduction_factor_direct_diffuse(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        // Test without window shading
        let base_height = 1.;
        let height = 1.25;
        let width = 4.;
        let tilt = 0.;
        let orientation = 0.;
        let window_shading = vec![];
        let _f_sky = 1.;

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions
                    .shading_reduction_factor_direct_diffuse(
                        base_height,
                        height,
                        width,
                        tilt,
                        orientation,
                        &window_shading,
                        t_it
                    )
                    .unwrap(),
                [
                    (0.0, 0.0),
                    (1.0, 1.0),
                    (0.0, 0.9948359023012616),
                    (0.0, 0.988044845054899),
                    (0.0, 0.983862781044833),
                    (1.0, 0.9816340106900042),
                    (1.0, 0.9808582136351668),
                    (1.0, 0.9791897009662986)
                ][t_idx],
            );
        }

        // Test window shading
        let window_shading_val = vec![
            WindowShadingObject::Overhang {
                depth: 0.5,
                distance: 0.5,
            },
            WindowShadingObject::SideFinLeft {
                depth: 0.25,
                distance: 0.1,
            },
            WindowShadingObject::SideFinRight {
                depth: 0.25,
                distance: 0.1,
            },
        ];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions
                    .shading_reduction_factor_direct_diffuse(
                        base_height,
                        height,
                        width,
                        tilt,
                        orientation,
                        &window_shading_val,
                        t_it
                    )
                    .unwrap(),
                [
                    (0.0, 0.0),
                    (0.9552620660549097, 0.5905238865770632),
                    (0.0, 0.5905238865770632),
                    (0.0, 0.5905238865770632),
                    (0.0, 0.5905238865770632),
                    (1.0, 0.5905238865770632),
                    (1.0, 0.5905238865770632),
                    (1.0, 0.5905238865770631)
                ][t_idx],
            );
        }

        // Test with different combination of tilt , orientation with and without shading
        let base_height = 1.;
        let height = 1.25;
        let width = 4.;
        let tilt = 90.;
        let orientation = 180.;
        let window_shading = vec![];
        let _f_sky = 0.5;

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions
                    .shading_reduction_factor_direct_diffuse(
                        base_height,
                        height,
                        width,
                        tilt,
                        orientation,
                        &window_shading,
                        t_it
                    )
                    .unwrap(),
                [
                    (0.0, 0.0),
                    (1.0, 1.0),
                    (1.0, 1.0),
                    (1.0, 1.0),
                    (1.0, 1.0),
                    (1.0, 1.0),
                    (1.0, 1.0),
                    (1.0, 1.0)
                ][t_idx],
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions
                    .shading_reduction_factor_direct_diffuse(
                        base_height,
                        height,
                        width,
                        tilt,
                        orientation,
                        &window_shading_val,
                        t_it
                    )
                    .unwrap(),
                [
                    (0.0, 0.0),
                    (1.0, 0.5905238865770648),
                    (1.0, 0.5905238865770633),
                    (1.0, 0.5905238865770633),
                    (1.0, 0.5905238865770633),
                    (1.0, 0.5905238865770632),
                    (1.0, 0.5905238865770633),
                    (1.0, 0.5905238865770632)
                ][t_idx],
            );
        }

        let base_height = 1.;
        let height = 1.25;
        let width = 4.;
        let tilt = 180.;
        let orientation = -180.;
        let window_shading = vec![];
        let _f_sky = 0.;

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions
                    .shading_reduction_factor_direct_diffuse(
                        base_height,
                        height,
                        width,
                        tilt,
                        orientation,
                        &window_shading,
                        t_it
                    )
                    .unwrap(),
                [
                    (0.0, 0.0),
                    (1.0, 1.0),
                    (1.0, 1.0),
                    (1.0, 1.0),
                    (1.0, 1.0),
                    (1.0, 1.0),
                    (1.0, 1.0),
                    (1.0, 1.0)
                ][t_idx],
            );
        }

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                external_conditions
                    .shading_reduction_factor_direct_diffuse(
                        base_height,
                        height,
                        width,
                        tilt,
                        orientation,
                        &window_shading_val,
                        t_it
                    )
                    .unwrap(),
                [
                    (0.0, 0.0),
                    (1.0, 0.5905238865770632),
                    (1.0, 0.5905238865770632),
                    (1.0, 0.5905238865770632),
                    (1.0, 0.5905238865770632),
                    (1.0, 0.5905238865770632),
                    (1.0, 0.5905238865770632),
                    (1.0, 0.5905238865770632)
                ][t_idx],
            );
        }
    }

    // test_no_shading_segments in Python is unnecessary to set up in Rust as already guaranteed by type system
}
