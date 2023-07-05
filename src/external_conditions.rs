use crate::simulation_time::{
    SimulationTime, SimulationTimeIteration, SimulationTimeIterator, HOURS_IN_DAY,
};
use itertools::Itertools;
use serde::Deserialize;

#[derive(Clone, Debug, Deserialize)]
pub enum DaylightSavingsConfig {
    #[serde(rename(deserialize = "applicable and taken into account"))]
    ApplicableAndTakenIntoAccount,
    #[serde(rename(deserialize = "applicable but not taken into account"))]
    ApplicableButNotTakenIntoAccount,
    #[serde(rename(deserialize = "not applicable"))]
    NotApplicable,
}

#[derive(Clone, Debug, Deserialize)]
pub struct ShadingSegment {
    pub number: usize,
    pub start: i32,
    pub end: i32,
    pub objects: Option<Vec<ShadingObject>>,
}

#[derive(Clone, Debug, Deserialize)]
pub struct ShadingObject {
    #[serde(rename(deserialize = "type"))]
    object_type: ShadingObjectType,
    height: f64,
    distance: f64,
}

#[derive(Clone, Debug, Deserialize)]
pub struct WindowShadingObject {
    #[serde(rename(deserialize = "type"))]
    object_type: WindowShadingObjectType,
    depth: f64,
    distance: f64,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ShadingObjectType {
    Obstacle,
    Overhang,
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum WindowShadingObjectType {
    Obstacle,
    Overhang,
    SideFinRight,
    SideFinLeft,
}

#[derive(Clone)]
pub struct ExternalConditions {
    simulation_time: SimulationTimeIterator,
    air_temps: Vec<f64>,
    wind_speeds: Vec<f64>,
    diffuse_horizontal_radiations: Vec<f64>,
    direct_beam_radiations: Vec<f64>,
    solar_reflectivity_of_ground: Vec<f64>,
    pub latitude: f64,
    pub longitude: f64,
    pub timezone: u32,
    pub start_day: u32,
    pub end_day: u32,
    time_series_step: f64,
    pub january_first: u32,
    pub daylight_savings: DaylightSavingsConfig,
    pub leap_day_included: bool,
    direct_beam_conversion_needed: bool,
    shading_segments: Vec<ShadingSegment>,
    extra_terrestrial_radiation: Vec<f64>,
    solar_declinations: Vec<f64>,
    solar_times: Vec<f64>,
    solar_hour_angles: Vec<f64>,
    solar_altitudes: Vec<f64>,
    solar_zenith_angles: Vec<f64>,
    solar_azimuth_angles: Vec<f64>,
    air_masses: Vec<f64>,
    f1_circumsolar_brightness_coefficients: Vec<f64>,
    f2_horizontal_brightness_coefficients: Vec<f64>,
}

impl ExternalConditions {
    pub fn new(
        simulation_time: SimulationTimeIterator,
        air_temps: Vec<f64>,
        wind_speeds: Vec<f64>,
        diffuse_horizontal_radiations: Vec<f64>,
        direct_beam_radiations: Vec<f64>,
        solar_reflectivity_of_ground: Vec<f64>,
        latitude: f64,
        longitude: f64,
        timezone: u32,
        start_day: u32,
        end_day: u32,
        time_series_step: f64,
        january_first: u32,
        daylight_savings: DaylightSavingsConfig,
        leap_day_included: bool,
        direct_beam_conversion_needed: bool,
        shading_segments: Vec<ShadingSegment>,
    ) -> Self {
        let days_in_year = if leap_day_included { 366 } else { 365 };
        let hours_in_year = days_in_year * HOURS_IN_DAY;
        let time_shift = init_time_shift(timezone, longitude);

        // # Calculate earth orbit deviation for each day of year
        let earth_orbit_deviations = (0..days_in_year)
            .map(|day| init_earth_orbit_deviation(day))
            .collect::<Vec<f64>>();

        let extra_terrestrial_radiation = (0..days_in_year)
            .map(|day| init_extra_terrestrial_radiation(earth_orbit_deviations[day as usize]))
            .collect::<Vec<f64>>();

        let solar_declinations = (0..days_in_year)
            .map(|day| init_solar_declination(earth_orbit_deviations[day as usize]))
            .collect::<Vec<f64>>();

        let equations_of_time = (0..days_in_year)
            .map(|day| init_equation_of_time(day))
            .collect::<Vec<f64>>();

        let solar_times = (0..hours_in_year)
            .map(|hour| {
                init_solar_time(
                    hour % 24,
                    equations_of_time[hour.div_euclid(24) as usize],
                    time_shift,
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
            simulation_time,
            air_temps,
            wind_speeds,
            diffuse_horizontal_radiations,
            direct_beam_radiations,
            solar_reflectivity_of_ground,
            latitude,
            longitude,
            timezone,
            start_day,
            end_day,
            time_series_step,
            january_first,
            daylight_savings,
            leap_day_included,
            direct_beam_conversion_needed,
            shading_segments,
            extra_terrestrial_radiation,
            solar_declinations,
            solar_times,
            solar_hour_angles,
            solar_altitudes,
            solar_zenith_angles,
            solar_azimuth_angles,
            air_masses,
            f1_circumsolar_brightness_coefficients,
            f2_horizontal_brightness_coefficients,
        }
    }

    pub fn next(&mut self) -> () {
        self.simulation_time.next();
    }

    pub fn air_temp(&self) -> f64 {
        self.air_temps[self
            .simulation_time
            .time_series_idx(self.start_day, self.time_series_step)]
    }

    pub fn air_temp_annual(&self) -> Option<f64> {
        if self.air_temps.len() != 8760 {
            return None;
        }
        let sum: f64 = self.air_temps.iter().sum();
        Some(sum / self.air_temps.len() as f64)
    }

    pub fn air_temp_monthly(&self) -> f64 {
        let (idx_start, idx_end) = self.simulation_time.current_month_start_end_hours();
        let (idx_start, idx_end) = (idx_start as usize, idx_end as usize);
        let air_temps_month = &self.air_temps[idx_start..idx_end];
        let sum: f64 = air_temps_month.iter().sum();
        sum / air_temps_month.len() as f64
    }

    pub fn wind_speed(&self) -> f64 {
        self.wind_speeds[self
            .simulation_time
            .time_series_idx(self.start_day, self.time_series_step)]
    }

    pub fn wind_speed_annual(&self) -> Option<f64> {
        if self.wind_speeds.len() != (8760.0 / self.time_series_step) as usize {
            return None;
        }
        let sum: f64 = self.wind_speeds.iter().sum();
        Some(sum / self.wind_speeds.len() as f64)
    }

    pub fn diffuse_horizontal_radiation(&self) -> f64 {
        self.diffuse_horizontal_radiations[self.simulation_time.current_index()]
    }

    pub fn direct_beam_radiation(&self) -> f64 {
        self.direct_beam_radiations[self.simulation_time.current_index()]
    }

    pub fn solar_reflectivity_of_ground(&self) -> f64 {
        self.solar_reflectivity_of_ground[self
            .simulation_time
            .time_series_idx(self.start_day, self.time_series_step)]
    }

    fn solar_angle_of_incidence(&self, tilt: f64, orientation: f64) -> f64 {
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
        // """

        //set up/ shadow some vars as radians for trig stuff
        let tilt = tilt.to_radians();
        let orientation = orientation.to_radians();
        let latitude = self.latitude.to_radians();
        let solar_declination: f64 =
            self.solar_declinations[self.simulation_time.current_day() as usize].to_radians();
        let solar_hour_angle: f64 =
            self.solar_hour_angles[self.simulation_time.current_hour() as usize].to_radians();

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

    // following references broken method in python
    // fn sun_surface_azimuth(&self, orientation: f64) -> f64 {
    //     // """  calculates the azimuth angle between sun and the inclined surface,
    //     // needed as input for the calculation of the irradiance in case of solar shading by objects
    //     //
    //     // Arguments:
    //     //
    //     // orientation    -- is the orientation angle of the inclined surface, expressed as the
    //     //                   geographical azimuth angle of the horizontal projection of the inclined
    //     //                   surface normal, -180 to 180, in degrees;
    //     //
    //     // """
    //     match self.solar_hour_angles()
    // }
    //
    // fn sun_surface_tilt(&self, tilt: f64) -> f64 {}

    fn direct_irradiance(&self, tilt: f64, orientation: f64) -> f64 {
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
        let direct_irradiance = self.direct_beam_radiation()
            * self
                .solar_angle_of_incidence(tilt, orientation)
                .to_radians()
                .cos();
        if direct_irradiance < 0.0 {
            0.0
        } else {
            direct_irradiance
        }
    }

    fn a_over_b(&self, tilt: f64, orientation: f64) -> f64 {
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
        let a = self
            .solar_angle_of_incidence(tilt, orientation)
            .to_radians()
            .cos();

        if a <= 0.0 {
            return 0.0;
        }

        let cosine_of_85_degrees = 85.0f64.to_radians().cos();
        let mut b = self.solar_zenith_angles[self.simulation_time.current_hour() as usize];
        if b < cosine_of_85_degrees {
            b = cosine_of_85_degrees;
        }

        a / b
    }

    fn diffuse_irradiance(&self, tilt: f64, orientation: f64) -> (f64, f64, f64, f64) {
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
        let gsol_d = self.diffuse_horizontal_radiation();
        let f1 = self.f1_circumsolar_brightness_coefficients[self.simulation_time.current_index()];
        let f2 = self.f2_horizontal_brightness_coefficients[self.simulation_time.current_index()];

        // # Calculate components of diffuse radiation
        let diffuse_irr_sky = gsol_d * (1.0 - f1) * ((1.0 + tilt.to_radians().cos()) / 2.0);
        let diffuse_irr_circumsolar = self.circumsolar_irradiance(tilt, orientation);
        let diffuse_irr_horiz = gsol_d * f2 * tilt.to_radians().sin();

        (
            diffuse_irr_sky + diffuse_irr_circumsolar + diffuse_irr_horiz,
            diffuse_irr_sky,
            diffuse_irr_circumsolar,
            diffuse_irr_horiz,
        )
    }

    fn ground_reflection_irradiance(&self, tilt: f64) -> f64 {
        // """  calculates the contribution of the ground reflection to the irradiance on the inclined surface,
        // determined as function of global horizontal irradiance, which in this case is calculated from the solar
        // altitude, diffuse and beam solar irradiance and the solar reflectivity of the ground
        //
        // Arguments:
        // tilt           -- is the tilt angle of the inclined surface from horizontal, measured
        //                   upwards facing, 0 to 180, in degrees;
        // """

        self.diffuse_horizontal_radiation()
            + self.direct_beam_radiation()
                * self.solar_altitudes[self.simulation_time.current_hour() as usize]
                    .to_radians()
                    .sin()
                * self.solar_reflectivity_of_ground()
                * ((1.0 - tilt.to_radians().cos()) / 2.0)
    }

    fn circumsolar_irradiance(&self, tilt: f64, orientation: f64) -> f64 {
        // """  calculates the circumsolar_irradiance
        //
        // Arguments:
        // tilt           -- is the tilt angle of the inclined surface from horizontal, measured
        //                   upwards facing, 0 to 180, in degrees;
        // orientation    -- is the orientation angle of the inclined surface, expressed as the
        //                   geographical azimuth angle of the horizontal projection of the inclined
        //                   surface normal, -180 to 180, in degrees;
        // """

        self.diffuse_horizontal_radiation()
            * self.f1_circumsolar_brightness_coefficients[self.simulation_time.current_index()]
            * self.a_over_b(tilt, orientation)
    }

    fn calculated_direct_irradiance(&self, tilt: f64, orientation: f64) -> f64 {
        // """  calculates the total direct irradiance on an inclined surface including circumsolar
        //
        // Arguments:
        // tilt           -- is the tilt angle of the inclined surface from horizontal, measured
        //                   upwards facing, 0 to 180, in degrees;
        // orientation    -- is the orientation angle of the inclined surface, expressed as the
        //                   geographical azimuth angle of the horizontal projection of the inclined
        //                   surface normal, -180 to 180, in degrees;
        // """

        self.direct_irradiance(tilt, orientation) + self.circumsolar_irradiance(tilt, orientation)
    }

    fn calculated_diffuse_irradiance(&self, tilt: f64, orientation: f64) -> f64 {
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

        let (diffuse_irr_total, _, diffuse_irr_circumsolar, _) =
            self.diffuse_irradiance(tilt, orientation);

        diffuse_irr_total - diffuse_irr_circumsolar + self.ground_reflection_irradiance(tilt)
    }

    pub fn calculated_total_solar_irradiance(&self, tilt: f64, orientation: f64) -> f64 {
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

        self.calculated_direct_irradiance(tilt, orientation)
            + self.calculated_diffuse_irradiance(tilt, orientation)
    }

    pub fn calculated_direct_diffuse_total_irradiance(
        &self,
        tilt: f64,
        orientation: f64,
        diffuse_breakdown: bool,
    ) -> CalculatedDirectDiffuseTotalIrradiance {
        // NB. the original Python implementation uses an internal cache for each timestep here
        // it only retains one set of results at a time so it may be of limited use, and i've skipped reimplementing it
        // in a first pass of conversion into Rust

        let (diffuse_irr_total, diffuse_irr_sky, diffuse_irr_circumsolar, diffuse_irr_horiz) =
            self.diffuse_irradiance(tilt, orientation);
        let ground_reflection_irradiance = self.ground_reflection_irradiance(tilt);

        let calculated_direct = self.direct_irradiance(tilt, orientation) + diffuse_irr_circumsolar;
        let calculated_diffuse =
            diffuse_irr_total - diffuse_irr_circumsolar + ground_reflection_irradiance;
        let total_irradiance = calculated_direct + calculated_diffuse;

        if diffuse_breakdown {
            CalculatedDirectDiffuseTotalIrradianceWithBreakdown((
                calculated_direct,
                calculated_diffuse,
                total_irradiance,
                DiffuseBreakdown {
                    sky: diffuse_irr_sky,
                    circumsolar: diffuse_irr_circumsolar,
                    horiz: diffuse_irr_horiz,
                    ground_refl: ground_reflection_irradiance,
                },
            ))
        } else {
            CalculatedDirectDiffuseTotalIrradianceWithoutBreakdown((
                calculated_direct,
                calculated_diffuse,
                total_irradiance,
            ))
        }
    }

    fn outside_solar_beam(&self, tilt: f64, orientation: f64) -> bool {
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

        let current_hour_idx = self.simulation_time.current_hour() as usize;

        let test1 = orientation - self.solar_azimuth_angles[current_hour_idx];
        let test2 = tilt - self.solar_altitudes[current_hour_idx];

        (-90.0 > test1 || test1 > 90.0) || (-90.0 > test2 || test2 > 90.0)
    }

    fn get_segment(&self) -> Result<ShadingSegment, &'static str> {
        // """ for complex (environment) shading objects, we need to know which
        // segment the azimuth of the sun occupies at each timestep
        //
        // """

        let current_hour_idx = self.simulation_time.current_hour() as usize;
        let azimuth = self.solar_azimuth_angles[current_hour_idx];

        match self
            .shading_segments
            .iter()
            .find(|&segment| azimuth < segment.start as f64 && azimuth > segment.end as f64)
        {
            Some(segment) => Ok(segment.clone()),
            None => Err("Solar segment was not found - this is an unexpected error"),
        }
    }

    fn obstacle_shading_height(
        &self,
        base_height_of_k: f64,
        height_of_obstacle: f64,
        horiz_distance_from_surface_to_obstacle: f64,
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
                * self.solar_altitudes[self.simulation_time.current_hour() as usize]
                    .to_radians()
                    .tan();
        if hshade < 0.0 {
            0.0
        } else {
            hshade
        }
    }

    // following references nonexistent solar_altitude method in the original python, so leaving incomplete
    // fn overhang_shading_height(
    //     &self,
    //     shaded_surface_height: f64,
    //     base_shaded_surface_height: f64,
    //     lowest_height_of_overhang: f64,
    //     horiz_distance_from_surface_to_overhang: f64,
    // ) -> f64 {
    //     // """ calculates the height of the shading on the shaded surface (k),
    //     // from the shading overhang in segment i at time t. Note that "overhang"
    //     // has a specific meaning in ISO 52016 Annex F
    //     //
    //     // Arguments:
    //     // shaded_surface_height - Hk            -- is the height of the shaded surface, k, in m
    //     // base_shaded_surface_height - Hkbase        -- is the base height of the shaded surface k, in m
    //     // lowest_height_of_overhang - Hovh          -- is the lowest height of the overhang q, in segment i, in m
    //     // horiz_distance_from_surface_to_overhang Lkovh         -- is the horizontal distance between the shaded surface k
    //     //                  and the shading overhang, q, in segment i, in m
    //     // """
    //     let hshade = shaded_surface_height + base_shaded_surface_height - lowest_height_of_overhang + horiz_distance_from_surface_to_overhang * self.solar_altitude().to_radians().tan();
    //     if hshade < 0.0 {
    //         0.0
    //     } else {
    //         hshade
    //     }
    // }

    // following method is incompletely transcoded from python as original method references a broken method
    // fn direct_shading_reduction_factor(
    //     &self,
    //     base_height: f64,
    //     height: f64,
    //     width: f64,
    //     orientation: f64,
    //     window_shading: f64,
    // ) -> Result<f64, &'static str> {
    //     // """ calculates the shading factor of direct radiation due to external
    //     // shading objects
    //     //
    //     // Arguments:
    //     // height         -- is the height of the shaded surface (if surface is tilted then
    //     //                   this must be the vertical projection of the height), in m
    //     // base_height    -- is the base height of the shaded surface k, in m
    //     // width          -- is the width of the shaded surface, in m
    //     // orientation    -- is the orientation angle of the inclined surface, expressed as the
    //     //                   geographical azimuth angle of the horizontal projection of the
    //     //                   inclined surface normal, -180 to 180, in degrees;
    //     // window_shading -- data on overhangs and side fins associated to this building element
    //     //                   includes the shading object type, depth, anf distance from element
    //     // """
    //
    //     // # start with default assumption of no shading
    //     let mut hshade_obst = 0.0;
    //     let mut hshade_ovh = 0.0;
    //     let mut wfinr = 0.0;
    //     let mut wfinl = 0.0;
    //
    //     // #first process the distant (environment) shading for this building element
    //     let segment = self.get_segment().unwrap();
    //
    //     if let Some(shading_objects) = segment.objects {
    //         for &shading_object in shading_objects {
    //             match shading_object.object_type {
    //                 ShadingObjectType::Obstacle => {
    //                     let new_shade_height = self.obstacle_shading_height(
    //                         base_height,
    //                         shading_object.height,
    //                         shading_object.distance,
    //                     );
    //                     if new_shade_height > hshade_obst {
    //                         hshade_obst = new_shade_height;
    //                     }
    //                 }
    //                 ShadingObjectType::Overhang => {
    //                     let new_shade_height = self.overhang_shading_height() //etc etc stopping as this method is incomplete
    //                 }
    //             }
    //         }
    //     }
    // }

    fn diffuse_shading_reduction_factor(
        &self,
        diffuse_breakdown: DiffuseBreakdown,
        tilt: f64,
        height: f64,
        width: f64,
        window_shading: Option<&Vec<WindowShadingObject>>,
    ) -> f64 {
        // """ calculates the shading factor of diffuse radiation due to external
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
        // window_shading -- data on overhangs and side fins associated to this building element
        //                   includes the shading object type, depth, and distance from element
        // """
        // # Note: Shading factor for circumsolar radiation is same as for direct.
        // #       As circumsolar radiation will be subtracted from diffuse and
        // #       added to direct later on, we don't need to do anything for
        // #       circumsolar radiation here and it is excluded.
        let diffuse_irr_sky = diffuse_breakdown.sky;
        let diffuse_irr_hor = diffuse_breakdown.horiz;
        let diffuse_irr_ref = diffuse_breakdown.ground_refl;
        let diffuse_irr_total = diffuse_irr_sky + diffuse_irr_hor + diffuse_irr_ref;

        // # TODO Calculate shading from distant objects (PD CEN ISO/TR 52016-2:2017 Section F.6.2)
        // # TODO Loop over segments
        // # TODO     Calculate sky-diffuse shading factor for the segment
        // # TODO     Calculate horiz-diffuse shading factor for the segment
        // # TODO Sum product of shading factor and irradiance for each element
        //
        // # Calculate shading from fins and overhangs (PD CEN ISO/TR 52016-2:2017 Section F.6.3)
        // # TODO Alpha is not defined in this standard but is possibly the angular
        // #      height of the horizon. This would make some sense as raising the
        // #      horizon angle would have essentially the same effect as increasing
        // #      the tilt of the building element so it would make sense to add it
        // #      to beta when calculating the sky view factor. The overall effect
        // #      of changing alpha when there are no fins or overhangs seems to be
        // #      to decrease the shading factor (meaning more shading) for diffuse
        // #      radiation from sky and horizon, and to increase the shading factor
        // #      for radiation reflected from the ground (sometimes to above 1),
        // #      which would also seem to make sense as in this case there is less
        // #      sky and more ground in view than in the basic assumption of
        // #      perfectly flat surroundings.
        // # TODO Should angular height of horizon be user input or derived from
        // #      calculation of distant shading objects above? Set to zero for now
        let angular_height_of_horizon = 0.0f64;
        let alpha = angular_height_of_horizon.to_radians();

        // # TODO Beta is not defined in this standard but is used for tilt of the
        // #      building element in BS EN ISO 52016-1:2017, so assuming the same.
        // #      This seems to give sensible numbers for the sky view factor
        // #      F_w_sky when alpha = 0
        let beta = tilt.to_radians();

        // #create lists of diffuse shading factors to keep the largest one
        // #in case there are multiple shading objects
        let mut fdiff_list: Vec<f64> = vec![];

        // # Unpack window shading details
        let mut ovh_D_L_ls = vec![(0.0, 1.0)]; // # [D,L] - L cannot be zero as this leads to divide-by-zero later on
        let mut finR_D_L_ls = vec![(0.0, 1.0)]; //# [D,L] - L cannot be zero as this leads to divide-by-zero later on
        let mut finL_D_L_ls = vec![(0.0, 1.0)]; // # [D,L] - L cannot be zero as this leads to divide-by-zero later on

        if let Some(shading_objects) = window_shading {
            for shading_object in shading_objects.iter() {
                //todo (maybe make some window shading types?)
                match shading_object.object_type {
                    WindowShadingObjectType::Obstacle => {
                        ovh_D_L_ls.push((shading_object.depth, shading_object.distance));
                    }
                    WindowShadingObjectType::SideFinLeft => {
                        finR_D_L_ls.push((shading_object.depth, shading_object.distance));
                    }
                    WindowShadingObjectType::SideFinRight => {
                        finL_D_L_ls.push((shading_object.depth, shading_object.distance));
                    }
                    WindowShadingObjectType::Overhang => {
                        // do nothing
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
            finR_D_L_ls.remove(0);
        }

        // #perform the diff shading calculation for each combination of overhangs and fins
        for iteration_vec in [ovh_D_L_ls, finR_D_L_ls, finL_D_L_ls]
            .iter()
            .multi_cartesian_product()
        {
            let (ovh_D_L, finR_D_L, finL_D_L) =
                (iteration_vec[0], iteration_vec[1], iteration_vec[2]);
            let d_ovh = ovh_D_L.0;
            let l_ovh = ovh_D_L.1;
            let d_finL = finL_D_L.0;
            let l_finL = finL_D_L.1;
            let d_finR = finR_D_L.0;
            let l_finR = finR_D_L.1;
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
            let f_w_s = (0.6514 * (1.0 - (p2_finL / (p1_finL.powi(2) + p2_finL.powi(2)).sqrt()))
                + 0.6514 * (1.0 - (p2_finR / (p1_finR.powi(2) + p2_finR.powi(2)).sqrt())))
                / 2.0;
            let f_w_o = 0.3282 * (1.0 - (p2_ovh / (p1_ovh.powi(2) + p2_ovh.powi(2)).sqrt()));
            let f_w_sky = (1.0 - (alpha + beta + 90.0f64.to_radians()).sin()) / 2.0;

            // # Calculate denominators of eqns F.9 to F.14
            let view_factor_sky_no_obstacles = (1.0 + beta.cos()) / 2.0;
            let view_factor_ground_no_obstacles = (1.0 - beta.cos()) / 2.0;

            // # Setback and remote obstacles (eqns F.9 and F.10): Top half of each eqn
            // # is view factor to sky (F.9) or ground (F.10) with setback and distant
            // # obstacles
            // # TODO Uncomment these lines when definitions of P1 and P2 in formula
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

            // # Overhangs and remote obstacles (eqns F.13 and F.14)
            // # Top half of eqn F.13 is view factor to sky with overhangs
            let f_sh_dif_overhangs = if view_factor_sky_no_obstacles == 0.0 {
                1.0
            } else {
                (f_w_sky - f_w_o) / view_factor_sky_no_obstacles
            };

            // # Top half of eqn F.14 is view factor to ground with distant obstacles,
            // # but does not account for overhangs blocking any part of the view of
            // # the ground, presumably because this will not happen in the vast
            // # majority of cases
            let f_sh_ref_overhangs = if view_factor_ground_no_obstacles == 0.0 {
                1.0
            } else {
                (1.0 - f_w_sky) / view_factor_ground_no_obstacles
            };

            // # Keep the smallest of the three shading reduction factors as the
            // # diffuse or reflected shading factor. Also enforce that these cannot be
            // # negative (which may happen with some extreme tilt values)
            // # TODO Add setback shading factors to the arguments to min function when
            // #      definitions of P1 and P2 in formula for F_w_r have been confirmed.
            // # F_sh_dif = max(0.0, min(F_sh_dif_setback, F_sh_dif_fins, F_sh_dif_overhangs))
            // # F_sh_ref = max(0.0, min(F_sh_ref_setback, F_sh_ref_fins, F_sh_ref_overhangs))
            let min_f_sh = if f_sh_dif_fins < f_sh_dif_overhangs {
                f_sh_dif_fins
            } else {
                f_sh_dif_overhangs
            };
            let f_sh_dif = if min_f_sh > 0.0 { min_f_sh } else { 0.0 };
            let min_f_sh = if f_sh_ref_fins < f_sh_ref_overhangs {
                f_sh_ref_fins
            } else {
                f_sh_ref_overhangs
            };
            let f_sh_ref = if min_f_sh > 0.0 { min_f_sh } else { 0.0 };

            let fdiff = (f_sh_dif * (diffuse_irr_sky + diffuse_irr_hor)
                + f_sh_ref * diffuse_irr_ref)
                / diffuse_irr_total;
            fdiff_list.push(fdiff);
        }

        // following is finding the max value of fdiff_list
        fdiff_list
            .iter()
            .max_by(|a, b| a.total_cmp(b))
            .unwrap()
            .clone()
    }

    // commenting out for now as uses broken method direct_shading_reduction_factor
    // pub fn shading_reduction_factor_direct_diffuse(
    //     &self,
    //     base_height: f64,
    //     height: f64,
    //     width: f64,
    //     tilt: f64,
    //     orientation: f64,
    //     window_shading: Vec<WindowShadingObject>,
    // ) -> (f64, f64) {
    //     // """ calculates the direct and diffuse shading factors due to external
    //     // shading objects
    //     //
    //     // Arguments:
    //     // height         -- is the height of the shaded surface (if surface is tilted then
    //     //                   this must be the vertical projection of the height), in m
    //     // base_height    -- is the base height of the shaded surface k, in m
    //     // width          -- is the width of the shaded surface, in m
    //     // orientation    -- is the orientation angle of the inclined surface, expressed as the
    //     //                   geographical azimuth angle of the horizontal projection of the
    //     //                   inclined surface normal, -180 to 180, in degrees;
    //     // tilt           -- is the tilt angle of the inclined surface from horizontal, measured
    //     //                   upwards facing, 0 to 180, in degrees;
    //     // window_shading -- data on overhangs and side fins associated to this building element
    //     //                   includes the shading object type, depth, anf distance from element
    //     // """
    //
    //     // # first chceck if there is any radiation. This is needed to prevent a potential
    //     // # divide by zero error in the final step, but also, if there is no radiation
    //     // # then shading is irrelevant and we can skip the whole calculation
    //     let CalculatedDirectDiffuseTotalIrradianceWithBreakdown((
    //         direct,
    //         diffuse,
    //         _,
    //         diffuse_breakdown,
    //     )) = self.calculated_direct_diffuse_total_irradiance(tilt, orientation, true);
    //     if direct + diffuse == 0.0 {
    //         return (0.0, 0.0);
    //     }
    //
    //     // # first check if the surface is outside the solar beam
    //     // # if so then direct shading is complete and we don't need to
    //     // # calculate shading from objects
    //     let fdir = if self.outside_solar_beam(tilt, orientation) {
    //         0.0
    //     } else {
    //         self.direct_shading_reduction_factor(
    //             base_height,
    //             height,
    //             width,
    //             orientation,
    //             window_shading,
    //         )
    //     };
    //
    //     let fdiff = self.diffuse_shading_reduction_factor(
    //         diffuse_breakdown,
    //         tilt,
    //         height,
    //         width,
    //         Some(&window_shading),
    //     );
    //
    //     (fdir, fdiff)
    // }
}

pub enum CalculatedDirectDiffuseTotalIrradiance {
    CalculatedDirectDiffuseTotalIrradianceWithBreakdown((f64, f64, f64, DiffuseBreakdown)),
    CalculatedDirectDiffuseTotalIrradianceWithoutBreakdown((f64, f64, f64)),
}

pub struct DiffuseBreakdown {
    sky: f64,
    circumsolar: f64,
    horiz: f64,
    ground_refl: f64,
}

pub fn init_direct_beam_radiation(
    direct_beam_conversion_needed: bool,
    raw_value: f64,
    solar_altitude: f64,
) -> f64 {
    // # if the climate data to only provide direct horizontal (rather than normal:
    // # If only direct (beam) solar irradiance at horizontal plane is available in the climatic data set,
    // # it shall be converted to normal incidence by dividing the value by the sine of the solar altitude.
    // """ ISO 52010 section 6.4.2
    // TODO investigate the impact of these notes further. Applicable for weather from CIBSE file.
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
            raw_value // # TODO should this be zero?
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

fn init_time_shift(timezone: u32, longitude: f64) -> f64 {
    // """ Calculate the time shift, in hours, resulting from the fact that the
    // longitude and the path of the sun are not equal
    //
    // NOTE Daylight saving time is disregarded in tshift which is time independent
    // """
    timezone as f64 - longitude / 15.0
}

fn init_solar_time(hour_of_day: u32, equation_of_time: f64, time_shift: f64) -> f64 {
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
    // # TODO How is this to be adjusted for timesteps that are not hourly?
    // # would allowing solar_time to be a decimal be all that is needed?
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
    // # TODO How is this to be adjusted for timesteps that are not hourly?
    // # would allowing solar_time to be a decimal be all that is needed?
    //
    // """
    // asol is the solar altitude angle, the angle between the solar beam
    // and the horizontal surface, in degrees;
    // """

    //all three params provided as degrees, but we need to shadow each as radians for trig calcs
    let latitude = latitude.to_radians();
    let solar_declination = solar_declination.to_radians();
    let solar_hour_angle = solar_hour_angle.to_radians();

    (solar_declination.sin() * latitude.sin()
        + solar_declination.cos() * latitude.cos() * solar_hour_angle.cos())
    .to_degrees()
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
    let sin_aux1_numerator = solar_declination.cos() * solar_hour_angle.cos();
    let cos_aux1_numerator =
        solar_declination.sin() + latitude.sin() * solar_declination.cos() * solar_hour_angle.cos();

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

use crate::external_conditions::BrightnessCoefficientName::{F11, F12, F13, F21, F22, F23};
use crate::external_conditions::CalculatedDirectDiffuseTotalIrradiance::{
    CalculatedDirectDiffuseTotalIrradianceWithBreakdown,
    CalculatedDirectDiffuseTotalIrradianceWithoutBreakdown,
};
use variants_struct::VariantsStruct;

#[derive(VariantsStruct)]
#[struct_name = "BrightnessCoefficientsRow"]
enum BrightnessCoefficientName {
    F11,
    F12,
    F13,
    F21,
    F22,
    F23,
}

type Bcn = BrightnessCoefficientName;

// version of Table 8 in ISO 52010
static BRIGHTNESS_COEFFICIENTS: [BrightnessCoefficientsRow<f64>; 8] = [
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

fn brightness_coefficient(e: f64, fij: BrightnessCoefficientName) -> &'static f64 {
    // """ returns brightness coefficient as a look up from Table 8 in ISO 52010
    //
    // Arguments:
    // e    -- dimensionless clearness parameter
    // fij  -- the coefficient to be returned. e.g. f12 or f23
    // """
    //
    // #TODO I've not had a need for the clearness index parameters contained in this table yet,
    // #if they are needed as input or output later then this function can be reworked
    BRIGHTNESS_COEFFICIENTS[if e < 1.065 {
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
    }]
    .get(&fij)
    .unwrap()
}

fn init_f1_circumsolar_brightness_coefficient(e: f64, delta: f64, solar_zenith_angle: f64) -> f64 {
    // """ returns the circumsolar brightness coefficient, F1
    //
    // Arguments:
    // E -- dimensionless clearness parameter for the current timestep
    // delta -- dimensionless sky brightness parameter for the current timestep
    // solar_zenith_angle -- solar zenith angle for the current hour
    // """
    let f1: f64 = brightness_coefficient(e, F11)
        + brightness_coefficient(e, F12) * delta
        + brightness_coefficient(e, F13) * (std::f64::consts::PI * solar_zenith_angle / 180.0);
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
    brightness_coefficient(e, F21)
        + brightness_coefficient(e, F22) * delta
        + brightness_coefficient(e, F23) * (std::f64::consts::PI * solar_zenith_angle / 180.0)
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
mod test {
    use super::*;
    use crate::external_conditions::DaylightSavingsConfig::NotApplicable;
    use crate::simulation_time::{SimulationTime, HOURS_IN_DAY};
    use rstest::*;

    const BASE_AIR_TEMPS: [f64; 24] = [
        0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 7.5, 10.0, 12.5, 15.0, 19.5, 17.0,
        15.0, 12.0, 10.0, 7.0, 5.0, 3.0, 1.0,
    ];

    #[fixture]
    pub fn simulation_time_iterator() -> SimulationTimeIterator {
        SimulationTime::new(0.0, 8.0, 1.0).iter()
    }

    #[fixture]
    pub fn air_temp_day_jan() -> [f64; 24] {
        BASE_AIR_TEMPS
    }

    #[fixture]
    pub fn air_temp_day_feb() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 1.0)
    }

    #[fixture]
    pub fn air_temp_day_mar() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 2.0)
    }

    #[fixture]
    pub fn air_temp_day_apr() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 3.0)
    }

    #[fixture]
    pub fn air_temp_day_may() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 4.0)
    }

    #[fixture]
    pub fn air_temp_day_jun() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 5.0)
    }

    #[fixture]
    pub fn air_temp_day_jul() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 6.0)
    }

    #[fixture]
    pub fn air_temp_day_aug() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 6.0)
    }

    #[fixture]
    pub fn air_temp_day_sep() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 5.0)
    }

    #[fixture]
    pub fn air_temp_day_oct() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 4.0)
    }

    #[fixture]
    pub fn air_temp_day_nov() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 3.0)
    }

    #[fixture]
    pub fn air_temp_day_dec() -> [f64; 24] {
        BASE_AIR_TEMPS.map(|temp| temp + 2.0)
    }

    #[fixture]
    pub fn air_temps() -> Vec<f64> {
        let mut temps: Vec<f64> = vec![];
        let months = [
            (air_temp_day_jan as fn() -> [f64; 24], 31),
            (air_temp_day_feb as fn() -> [f64; 24], 28),
            (air_temp_day_mar as fn() -> [f64; 24], 31),
            (air_temp_day_apr as fn() -> [f64; 24], 30),
            (air_temp_day_may as fn() -> [f64; 24], 31),
            (air_temp_day_jun as fn() -> [f64; 24], 30),
            (air_temp_day_jul as fn() -> [f64; 24], 31),
            (air_temp_day_aug as fn() -> [f64; 24], 31),
            (air_temp_day_sep as fn() -> [f64; 24], 30),
            (air_temp_day_oct as fn() -> [f64; 24], 31),
            (air_temp_day_nov as fn() -> [f64; 24], 30),
            (air_temp_day_dec as fn() -> [f64; 24], 31),
        ];
        for (temps_fn, month_days_count) in months {
            temps.extend_from_slice(
                temps_fn()
                    .to_vec()
                    .iter()
                    .cloned()
                    .cycle()
                    .take(month_days_count * HOURS_IN_DAY as usize)
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
    pub fn wind_speed_day_jan() -> [f64; 24] {
        BASE_WIND_SPEEDS
    }

    #[fixture]
    pub fn wind_speed_day_feb() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 0.1)
    }

    #[fixture]
    pub fn wind_speed_day_mar() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 0.2)
    }

    #[fixture]
    pub fn wind_speed_day_apr() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 0.6)
    }

    #[fixture]
    pub fn wind_speed_day_may() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 0.8)
    }

    #[fixture]
    pub fn wind_speed_day_jun() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 1.1)
    }

    #[fixture]
    pub fn wind_speed_day_jul() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 1.2)
    }

    #[fixture]
    pub fn wind_speed_day_aug() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 1.2)
    }

    #[fixture]
    pub fn wind_speed_day_sep() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 1.1)
    }

    #[fixture]
    pub fn wind_speed_day_oct() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 0.7)
    }

    #[fixture]
    pub fn wind_speed_day_nov() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 0.5)
    }

    #[fixture]
    pub fn wind_speed_day_dec() -> [f64; 24] {
        BASE_WIND_SPEEDS.map(|speed| speed - 0.3)
    }

    #[fixture]
    pub fn wind_speeds() -> Vec<f64> {
        let mut speeds: Vec<f64> = vec![];
        let months = [
            (wind_speed_day_jan as fn() -> [f64; 24], 31),
            (wind_speed_day_feb as fn() -> [f64; 24], 28),
            (wind_speed_day_mar as fn() -> [f64; 24], 31),
            (wind_speed_day_apr as fn() -> [f64; 24], 30),
            (wind_speed_day_may as fn() -> [f64; 24], 31),
            (wind_speed_day_jun as fn() -> [f64; 24], 30),
            (wind_speed_day_jul as fn() -> [f64; 24], 31),
            (wind_speed_day_aug as fn() -> [f64; 24], 31),
            (wind_speed_day_sep as fn() -> [f64; 24], 30),
            (wind_speed_day_oct as fn() -> [f64; 24], 31),
            (wind_speed_day_nov as fn() -> [f64; 24], 30),
            (wind_speed_day_dec as fn() -> [f64; 24], 31),
        ];
        for (speeds_fn, month_days_count) in months {
            speeds.extend_from_slice(
                speeds_fn()
                    .to_vec()
                    .iter()
                    .cloned()
                    .cycle()
                    .take(month_days_count * HOURS_IN_DAY as usize)
                    .collect::<Vec<f64>>()
                    .as_slice(),
            );
        }

        speeds
    }

    #[fixture]
    pub fn diffuse_horizontal_radiation() -> [f64; 8] {
        [333.0, 610.0, 572.0, 420.0, 0.0, 10.0, 90.0, 275.0]
    }

    #[fixture]
    pub fn direct_beam_radiation() -> [f64; 8] {
        [420.0, 750.0, 425.0, 500.0, 0.0, 40.0, 0.0, 388.0]
    }

    #[fixture]
    pub fn solar_reflectivity_of_ground() -> [f64; 8760] {
        [0.2; 8760]
    }

    #[fixture]
    pub fn latitude() -> f64 {
        51.42
    }

    #[fixture]
    pub fn longitude() -> f64 {
        -0.75
    }

    #[fixture]
    pub fn timezone() -> u32 {
        0
    }

    #[fixture]
    pub fn start_day() -> u32 {
        0
    }

    #[fixture]
    pub fn end_day() -> u32 {
        0
    }

    #[fixture]
    pub fn time_series_step() -> f64 {
        1.0
    }

    #[fixture]
    pub fn january_first() -> u32 {
        1
    }

    #[fixture]
    pub fn daylight_savings() -> DaylightSavingsConfig {
        NotApplicable
    }

    #[fixture]
    pub fn leap_day_included() -> bool {
        false
    }

    #[fixture]
    pub fn direct_beam_conversion_needed() -> bool {
        false
    }

    #[fixture]
    pub fn shading_segments() -> Vec<ShadingSegment> {
        vec![
            ShadingSegment {
                number: 1,
                start: 180,
                end: 135,
                objects: None,
            },
            ShadingSegment {
                number: 2,
                start: 135,
                end: 90,
                objects: None,
            },
            ShadingSegment {
                number: 3,
                start: 90,
                end: 45,
                objects: None,
            },
            ShadingSegment {
                number: 4,
                start: 45,
                end: 0,
                objects: None,
            },
            ShadingSegment {
                number: 5,
                start: 0,
                end: -45,
                objects: None,
            },
            ShadingSegment {
                number: 6,
                start: -45,
                end: -90,
                objects: None,
            },
            ShadingSegment {
                number: 7,
                start: -90,
                end: -135,
                objects: None,
            },
            ShadingSegment {
                number: 8,
                start: -135,
                end: -180,
                objects: None,
            },
        ]
    }

    #[fixture]
    pub fn external_conditions() -> ExternalConditions {
        ExternalConditions::new(
            simulation_time_iterator(),
            air_temps(),
            wind_speeds(),
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
    fn should_have_correct_air_temps(
        external_conditions: ExternalConditions,
        air_temps: Vec<f64>,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let mut external_conditions = external_conditions.clone();
        for (i, simtime_step) in simulation_time_iterator.enumerate() {
            external_conditions.next();
            assert_eq!(
                external_conditions.air_temp(),
                air_temps[i],
                "failed on iteration index {} with step {:?}",
                i,
                simtime_step
            );
        }
    }

    fn round_by_precision(src: f64, precision: f64) -> f64 {
        (precision * src).round() / precision
    }

    #[rstest]
    fn should_have_correct_air_temp_annual(external_conditions: ExternalConditions) {
        let precision = 1e6;
        assert_eq!(
            round_by_precision(external_conditions.air_temp_annual().unwrap(), 1e6),
            round_by_precision(10.1801369863014, 1e6)
        );
    }

    #[rstest]
    fn should_have_correct_air_temp_monthly(
        external_conditions: ExternalConditions,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let expected_monthly_air_temps: [f64; 12] = [
            6.75, 7.75, 8.75, 9.75, 10.75, 11.75, 12.75, 12.75, 11.75, 10.75, 9.75, 8.75,
        ];
        let mut external_conditions = external_conditions.clone();
        for (i, simtime_step) in simulation_time_iterator.enumerate() {
            external_conditions.next();
            let month_idx = simtime_step.current_month().unwrap() as usize;
            assert_eq!(
                external_conditions.air_temp_monthly(),
                expected_monthly_air_temps[month_idx]
            );
        }
    }

    #[rstest]
    fn should_have_correct_wind_speeds(
        external_conditions: ExternalConditions,
        wind_speeds: Vec<f64>,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let mut external_conditions = external_conditions.clone();
        for (i, simtime_step) in simulation_time_iterator.enumerate() {
            external_conditions.next();
            assert_eq!(external_conditions.wind_speed(), wind_speeds[i]);
        }
    }

    #[rstest]
    fn should_have_correct_wind_speed_annual(external_conditions: ExternalConditions) {
        assert_eq!(
            round_by_precision(external_conditions.wind_speed_annual().unwrap(), 0.01),
            round_by_precision(4.23, 0.01)
        );
    }

    #[rstest]
    fn should_have_correct_diffuse_horizontal_radiation(
        external_conditions: ExternalConditions,
        diffuse_horizontal_radiation: [f64; 8],
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let mut external_conditions = external_conditions.clone();
        for (i, simtime_step) in simulation_time_iterator.enumerate() {
            external_conditions.next();
            assert_eq!(
                external_conditions.diffuse_horizontal_radiation(),
                diffuse_horizontal_radiation[i]
            );
        }
    }

    #[rstest]
    fn should_have_correct_direct_beam_radiation(
        external_conditions: ExternalConditions,
        direct_beam_radiation: [f64; 8],
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let mut external_conditions = external_conditions.clone();
        for (i, simtime_step) in simulation_time_iterator.enumerate() {
            external_conditions.next();
            assert_eq!(
                external_conditions.direct_beam_radiation(),
                direct_beam_radiation[i]
            );
        }
    }

    #[rstest]
    fn should_have_correct_solar_reflectivity_of_ground(
        external_conditions: ExternalConditions,
        solar_reflectivity_of_ground: [f64; 8760],
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let mut external_conditions = external_conditions.clone();
        for (i, simtime_step) in simulation_time_iterator.enumerate() {
            external_conditions.next();
            assert_eq!(
                external_conditions.solar_reflectivity_of_ground(),
                solar_reflectivity_of_ground[i]
            );
        }
    }
}
