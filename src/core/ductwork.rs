use crate::input::MVHRLocation;
use std::error::Error;
use std::f64::consts::PI;

// Set default value for the heat transfer coefficient inside the duct, in W / m^2 K
const INTERNAL_HTC: f64 = 15.5; // CIBSE Guide C, Table 3.25, air flow rate approx 3 m/s

// Set default values for the heat transfer coefficient at the outer surface, in W / m^2 K
const EXTERNAL_REFLECTIVE_HTC: f64 = 5.7; // low emissivity reflective surface, CIBSE Guide C, Table 3.25
const EXTERNAL_NONREFLECTIVE_HTC: f64 = 10.0; // high emissivity non-reflective surface, CIBSE Guide C, Table 3.25

/// A struct to represent ductwork for mechanical ventilation with heat recovery
/// (MVHR), assuming steady state heat transfer in a hollow cyclinder (duct)
/// with radial heat flow. ISO 12241:2022
pub struct Ductwork {
    length_in_in_m: f64,
    length_out_in_m: f64,
    mvhr_location: MVHRLocation,
    diameter_including_insulation_in_m: f64,
    internal_surface_resistance: f64, // in K m / W
    insulation_resistance: f64,       // in K m / W
    external_surface_resistance: f64, // K m / W
}

impl Ductwork {
    /// Arguments:
    /// * `internal_diameter` - internal diameter of the duct, in m
    /// * `external_diameter` - external diameter of the duct, in m
    /// * `length_in` - length of intake duct, in m
    /// * `length_out` - length of exhaust duct, in m
    /// * `k_insulation` - thermal conductivity of the insulation, in W / m K
    /// * `thickness_insulation` - thickness of the duct insulation, in m
    /// * `reflective` - whether the outer surface of the duct is reflective (true) or not (false) (boolean input)
    /// * `mvhr_location` - location of the MVHR unit (inside or outside the thermal envelope)
    pub fn new(
        internal_diameter_in_m: f64,
        external_diameter_in_m: f64,
        length_in_in_m: f64,
        length_out_in_m: f64,
        k_insulation: f64,
        thickness_insulation_in_m: f64,
        reflective: bool,
        mvhr_location: MVHRLocation,
    ) -> Self {
        let external_htc = if reflective {
            EXTERNAL_REFLECTIVE_HTC
        } else {
            EXTERNAL_NONREFLECTIVE_HTC
        };

        // Calculate the diameter of the duct including the insulation (D_ins, here diameter_including_insulation_in_m), in m
        let diameter_including_insulation_in_m =
            external_diameter_in_m + (2. * thickness_insulation_in_m);

        // Calculate the interior linear surface resistance, in K m / W
        let internal_surface_resistance = 1. / (INTERNAL_HTC * PI * internal_diameter_in_m);

        // Calculate the insulation linear thermal resistance, in K m / W
        let insulation_resistance = (diameter_including_insulation_in_m / internal_diameter_in_m)
            .ln()
            / (2. * PI * k_insulation);

        // Calculate the exterior linear surface resistance, in K m / W
        let external_surface_resistance =
            1. / (external_htc * PI * diameter_including_insulation_in_m);

        Self {
            length_in_in_m,
            length_out_in_m,
            mvhr_location,
            diameter_including_insulation_in_m,
            internal_surface_resistance,
            insulation_resistance,
            external_surface_resistance,
        }
    }

    pub fn get_mvhr_location(&self) -> &MVHRLocation {
        &self.mvhr_location
    }

    /// Return the heat loss for air inside the duct for the current timestep
    /// Arguments:
    /// * `inside_temp_in_c` - temperature of air inside the duct, in degrees C
    /// * `outside_temp_in_c` - temperature outside the duct, in degrees C
    /// * `length`
    pub fn duct_heat_loss(
        &self,
        inside_temp_in_c: f64,
        outside_temp_in_c: f64,
        length: f64,
    ) -> f64 {
        // Calculate total thermal resistance
        let total_resistance = self.internal_surface_resistance
            + self.insulation_resistance
            + self.external_surface_resistance;

        // Calculate heat loss, in W
        (inside_temp_in_c - outside_temp_in_c) / total_resistance * length
    }

    /// Return the heat loss for air inside the duct for the current timestep
    /// Arguments:
    /// * `supply_duct_temp_in_c` - temperature of air inside the supply duct, in degrees C
    /// * `extract_duct_temp_in_c` - temperature of air inside the extract duct, in degrees C
    /// * `intake_duct_temp_in_c` - temperature of air inside the intake duct, in degrees C
    /// * `exhaust_duct_temp_in_c` - temperature of air inside the exhaust duct, in degrees C
    /// * `efficiency` - heat recovery efficiency of MVHR
    pub fn total_duct_heat_loss(
        &self,
        supply_duct_temp_in_c: Option<f64>,
        extract_duct_temp_in_c: Option<f64>,
        intake_duct_temp_in_c: Option<f64>,
        exhaust_duct_temp_in_c: Option<f64>,
        efficiency: f64,
    ) -> Result<f64, &str> {
        match self.mvhr_location {
            // Outside location
            // Air inside the duct loses heat, external environment gains heat
            // Loses energy to outside in extract duct - losses must be X by the efficiency of heat recovery
            // Loses energy to outside in supply duct - lose all because after MVHR unit
            MVHRLocation::Outside => {
                if supply_duct_temp_in_c.is_none()
                    || extract_duct_temp_in_c.is_none()
                    || intake_duct_temp_in_c.is_none()
                {
                    return Err("Duct temperatures not provided for outside MVHR.");
                }
                let outside_temp = intake_duct_temp_in_c.unwrap();
                let supply_heat_loss = self.duct_heat_loss(
                    supply_duct_temp_in_c.unwrap(),
                    outside_temp,
                    self.length_in_in_m,
                );
                let extract_heat_loss = self.duct_heat_loss(
                    extract_duct_temp_in_c.unwrap(),
                    outside_temp,
                    self.length_out_in_m,
                );
                Ok(-(supply_heat_loss + (extract_heat_loss * efficiency)))
            }
            // Inside location
            // This will be a negative heat loss i.e. air inside the duct gains heat, dwelling loses heat
            // Gains energy from zone in intake duct - benefit of gain must be X by the efficiency of heat recovery
            // Gains energy from zone in exhaust duct
            MVHRLocation::Inside => {
                if intake_duct_temp_in_c.is_none()
                    || exhaust_duct_temp_in_c.is_none()
                    || extract_duct_temp_in_c.is_none()
                {
                    return Err("Duct temperatures not provided for inside MVHR.");
                }
                let outside_temp = extract_duct_temp_in_c.unwrap();
                let intake_heat_loss = self.duct_heat_loss(
                    intake_duct_temp_in_c.unwrap(),
                    outside_temp,
                    self.length_in_in_m,
                );
                let exhaust_heat_loss = self.duct_heat_loss(
                    exhaust_duct_temp_in_c.unwrap(),
                    outside_temp,
                    self.length_out_in_m,
                );
                Ok((intake_heat_loss * efficiency) + exhaust_heat_loss)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::MVHRLocation::Inside;
    use crate::simulation_time::SimulationTime;
    use rstest::*;

    #[fixture]
    pub fn ductwork() -> Ductwork {
        Ductwork::new(0.025, 0.027, 0.4, 0.4, 0.02, 0.022, false, Inside)
    }

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 8., 1.)
    }

    #[rstest]
    pub fn should_have_correct_diameter(ductwork: Ductwork) {
        assert_eq!(
            round_by_precision(ductwork.diameter_including_insulation_in_m, 1e3),
            0.071,
            "incorrect diameter returned"
        );
    }

    #[rstest]
    pub fn should_have_correct_internal_surface_resistance(ductwork: Ductwork) {
        assert_eq!(
            round_by_precision(ductwork.internal_surface_resistance, 1e5),
            0.82144,
            "incorrect internal surface resistance returned"
        );
    }

    #[rstest]
    pub fn should_have_correct_insulation_resistance(ductwork: Ductwork) {
        assert_eq!(
            round_by_precision(ductwork.insulation_resistance, 1e5),
            8.30633,
            "incorrect insulation resistance returned"
        );
    }

    #[rstest]
    pub fn should_have_correct_external_surface_resistance(ductwork: Ductwork) {
        assert_eq!(
            round_by_precision(ductwork.external_surface_resistance, 1e5),
            0.44832,
            "incorrect external surface resistance returned"
        );
    }

    #[rstest]
    pub fn should_calc_correct_duct_heat_loss(ductwork: Ductwork, simulation_time: SimulationTime) {
        let outside_temp = [20.0, 19.5, 19.0, 18.5, 19.0, 19.5, 20.0, 20.5];
        let inside_temp = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0];

        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_eq!(
                round_by_precision(
                    ductwork.duct_heat_loss(inside_temp[t_idx], outside_temp[t_idx], 0.4),
                    1e5
                ),
                [-0.62656, -0.56390, -0.50125, -0.43859, -0.41771, -0.39682, -0.37594, -0.35505]
                    [t_idx],
                "incorrect heat loss returned"
            );
        }
    }

    #[rstest]
    pub fn should_calc_correct_total_duct_heat_loss(
        ductwork: Ductwork,
        simulation_time: SimulationTime,
    ) {
        // let outside_temp = [20.0, 19.5, 19.0, 18.5, 19.0, 19.5, 20.0, 20.5]; // unused param from Python code
        let intake_temp = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0];
        let exhaust_temp = [20.0, 19.5, 19.0, 18.5, 19.0, 19.5, 20.0, 20.5];
        let supply_temp = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0];
        let extract_temp = [20.0, 19.5, 19.0, 18.5, 19.0, 19.5, 20.0, 20.5];

        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_eq!(
                round_by_precision(
                    ductwork
                        .total_duct_heat_loss(
                            Some(supply_temp[t_idx]),
                            Some(extract_temp[t_idx]),
                            Some(intake_temp[t_idx]),
                            Some(exhaust_temp[t_idx]),
                            0.7
                        )
                        .unwrap(),
                    1e5
                ),
                [-0.43859, -0.39473, -0.35087, -0.30701, -0.29239, -0.27777, -0.26316, -0.24854]
                    [t_idx],
                "incorrect total heat loss returned"
            );
        }
    }

    fn round_by_precision(src: f64, precision: f64) -> f64 {
        (precision * src).round() / precision
    }
}
