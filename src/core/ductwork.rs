use crate::input::{DuctShape, DuctType};
use anyhow::anyhow;
use std::f64::consts::PI;

// Set default value for the heat transfer coefficient inside the duct, in W / m^2 K
const INTERNAL_HTC: f64 = 15.5; // CIBSE Guide C, Table 3.25, air flow rate approx 3 m/s

// Set default values for the heat transfer coefficient at the outer surface, in W / m^2 K
const EXTERNAL_REFLECTIVE_HTC: f64 = 5.7;
// low emissivity reflective surface, CIBSE Guide C, Table 3.25
const EXTERNAL_NONREFLECTIVE_HTC: f64 = 10.0; // high emissivity non-reflective surface, CIBSE Guide C, Table 3.25

/// A struct to represent ductwork for mechanical ventilation with heat recovery (MVHR)
///
/// Assumes steady state heat transfer in, 1. a hollow cylinder (duct)
/// with radial heat flow and 2. a rectangular cross-section ISO 12241:2022
#[derive(Clone, Copy, Debug)]
pub struct Ductwork {
    length_in_m: f64,
    internal_surface_resistance: f64,
    duct_type: DuctType,
    // in K m / W
    insulation_resistance: f64,
    // in K m / W
    external_surface_resistance: f64, // K m / W
    #[cfg(test)]
    pub diameter_including_insulation_in_m: Option<f64>,
}

impl Ductwork {
    /// Arguments:
    /// * `cross_section_shape` - whether cross-section of duct is circular or rectangular (square)
    /// * `duct_perimeter` - if ductwork is rectangular(square) enter perimeter, in m
    /// * `internal_diameter` - internal diameter of the duct, in m
    /// * `external_diameter` - external diameter of the duct, in m
    /// * `length` - length of duct, in m
    /// * `k_insulation` - thermal conductivity of the insulation, in W / m K
    /// * `thickness_insulation` - thickness of the duct insulation, in m
    /// * `reflective` - whether the outer surface of the duct is reflective
    /// * `duct_type` - intake, supply, extract or exhaust
    pub fn new(
        cross_section_shape: DuctShape,
        duct_perimeter_in_m: Option<f64>,
        internal_diameter_in_m: Option<f64>,
        external_diameter_in_m: Option<f64>,
        length_in_m: f64,
        k_insulation: f64,
        thickness_insulation_in_m: f64,
        reflective: bool,
        duct_type: DuctType,
    ) -> anyhow::Result<Self> {
        let external_htc = if reflective {
            EXTERNAL_REFLECTIVE_HTC
        } else {
            EXTERNAL_NONREFLECTIVE_HTC
        };

        let (
            _diameter_including_insulation_in_m,
            internal_surface_resistance,
            insulation_resistance,
            external_surface_resistance,
        ) = match cross_section_shape {
            DuctShape::Circular => {
                let external_diameter_in_m = external_diameter_in_m.ok_or_else(|| {
                    anyhow!("An external diameter is needed for circular ductwork.")
                })?;
                let internal_diameter_in_m = internal_diameter_in_m.ok_or_else(|| {
                    anyhow!("An internal diameter is needed for circular ductwork.")
                })?;

                // Calculate the diameter of the duct including the insulation (D_ins, here diameter_including_insulation_in_m), in m
                let diameter_including_insulation_in_m =
                    external_diameter_in_m + (2. * thickness_insulation_in_m);

                // Calculate the interior linear surface resistance, in K m / W
                let internal_surface_resistance = 1. / (INTERNAL_HTC * PI * internal_diameter_in_m);

                // Calculate the insulation linear thermal resistance, in K m / W
                let insulation_resistance =
                    (diameter_including_insulation_in_m / external_diameter_in_m).ln()
                        / (2. * PI * k_insulation);

                // Calculate the exterior linear surface resistance, in K m / W
                let external_surface_resistance =
                    1. / (external_htc * PI * diameter_including_insulation_in_m);

                (
                    Some(diameter_including_insulation_in_m),
                    internal_surface_resistance,
                    insulation_resistance,
                    external_surface_resistance,
                )
            }
            DuctShape::Rectangular => {
                let duct_perimeter_in_m = duct_perimeter_in_m.ok_or_else(|| {
                    anyhow!("Duct perimeter was expected to be provided for rectangular ductwork.")
                })?;

                // Calculate the perimeter of the duct including the insulation, in m
                // the value 8 is specified in the standard ISO 12241:2022 and not assigned a description
                let duct_perimeter_external =
                    duct_perimeter_in_m + (8. * thickness_insulation_in_m);

                // Calculate the interior linear surface resistance, in K m / W
                let internal_surface_resistance = 1.0 / (INTERNAL_HTC * duct_perimeter_in_m);

                // Calculate the insulation linear thermal resistance, in K m / W
                let insulation_resistance = (2.0 * thickness_insulation_in_m)
                    / (k_insulation * (duct_perimeter_in_m + duct_perimeter_external));

                // Calculate the exterior linear surface resistance, in K m / W
                let external_surface_resistance = 1.0 / (external_htc * duct_perimeter_external);

                (
                    None,
                    internal_surface_resistance,
                    insulation_resistance,
                    external_surface_resistance,
                )
            }
        };

        Ok(Self {
            length_in_m,
            internal_surface_resistance,
            insulation_resistance,
            external_surface_resistance,
            duct_type,
            #[cfg(test)]
            diameter_including_insulation_in_m: _diameter_including_insulation_in_m,
        })
    }

    pub fn duct_type(&self) -> DuctType {
        self.duct_type
    }

    /// Return the heat loss for air inside the duct for the current timestep
    /// Arguments:
    /// * `inside_temp_in_c` - temperature of air inside the duct, in degrees C
    /// * `outside_temp_in_c` - temperature outside the duct, in degrees C
    pub fn duct_heat_loss(&self, inside_temp_in_c: f64, outside_temp_in_c: f64) -> f64 {
        // Calculate total thermal resistance
        let total_resistance = self.internal_surface_resistance
            + self.insulation_resistance
            + self.external_surface_resistance;

        // Calculate heat loss, in W
        (inside_temp_in_c - outside_temp_in_c) / total_resistance * self.length_in_m
    }

    #[cfg(test)]
    pub fn internal_surface_resistance(&self) -> f64 {
        self.internal_surface_resistance
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation_time::SimulationTime;
    use approx::assert_relative_eq;
    use rstest::*;

    #[fixture]
    fn circular_ductwork() -> [Ductwork; 8] {
        [
            Ductwork::new(
                DuctShape::Circular,
                None,
                Some(0.025),
                Some(0.027),
                0.4,
                0.02,
                0.022,
                false,
                DuctType::Exhaust,
            )
            .unwrap(),
            Ductwork::new(
                DuctShape::Circular,
                None,
                Some(0.025),
                Some(0.027),
                0.4,
                0.02,
                0.022,
                false,
                DuctType::Intake,
            )
            .unwrap(),
            Ductwork::new(
                DuctShape::Circular,
                None,
                Some(0.025),
                Some(0.027),
                0.4,
                0.02,
                0.022,
                false,
                DuctType::Supply,
            )
            .unwrap(),
            Ductwork::new(
                DuctShape::Circular,
                None,
                Some(0.025),
                Some(0.027),
                0.4,
                0.02,
                0.022,
                false,
                DuctType::Extract,
            )
            .unwrap(),
            Ductwork::new(
                DuctShape::Circular,
                None,
                Some(0.025),
                Some(0.027),
                0.4,
                0.02,
                0.022,
                false,
                DuctType::Exhaust,
            )
            .unwrap(),
            Ductwork::new(
                DuctShape::Circular,
                None,
                Some(0.025),
                Some(0.027),
                0.4,
                0.02,
                0.022,
                false,
                DuctType::Intake,
            )
            .unwrap(),
            Ductwork::new(
                DuctShape::Circular,
                None,
                Some(0.025),
                Some(0.027),
                0.4,
                0.02,
                0.022,
                false,
                DuctType::Supply,
            )
            .unwrap(),
            Ductwork::new(
                DuctShape::Circular,
                None,
                Some(0.025),
                Some(0.027),
                0.4,
                0.02,
                0.022,
                false,
                DuctType::Extract,
            )
            .unwrap(),
        ]
    }

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 8., 1.)
    }

    #[rstest]
    fn test_get_duct_type(circular_ductwork: [Ductwork; 8]) {
        assert_eq!(circular_ductwork[0].duct_type(), DuctType::Exhaust);
        assert_eq!(circular_ductwork[1].duct_type(), DuctType::Intake);
        assert_eq!(circular_ductwork[2].duct_type(), DuctType::Supply);
        assert_eq!(circular_ductwork[3].duct_type(), DuctType::Extract);
    }

    #[rstest]
    fn test_duct_heat_loss(circular_ductwork: [Ductwork; 8], simulation_time: SimulationTime) {
        let outside_temp = [20.0, 19.5, 19.0, 18.5, 19.0, 19.5, 20.0, 20.5];
        let inside_temp = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0];
        let ductwork = circular_ductwork[0];

        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                ductwork.duct_heat_loss(inside_temp[t_idx], outside_temp[t_idx]),
                [
                    -0.669369, -0.602432, -0.535495, -0.468558, -0.446246, -0.423933, -0.401621,
                    -0.379309
                ][t_idx],
                max_relative = 1e-5
            );
        }
    }

    #[fixture]
    fn rectangular_ductwork() -> Ductwork {
        Ductwork::new(
            DuctShape::Rectangular,
            Some(0.1),
            None,
            None,
            0.4,
            0.02,
            0.022,
            false,
            DuctType::Exhaust,
        )
        .unwrap()
    }

    #[rstest]
    fn test_get_duct_type_rectangular(rectangular_ductwork: Ductwork) {
        assert_eq!(rectangular_ductwork.duct_type(), DuctType::Exhaust);
    }

    #[rstest]
    fn test_duct_heat_loss_rectangular(
        rectangular_ductwork: Ductwork,
        simulation_time: SimulationTime,
    ) {
        let outside_temp = [20.0, 19.5, 19.0, 18.5, 19.0, 19.5, 20.0, 20.5];
        let inside_temp = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0];
        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                rectangular_ductwork.duct_heat_loss(inside_temp[t_idx], outside_temp[t_idx]),
                [-0.87482, -0.787339, -0.69986, -0.61237, -0.58321, -0.55405, -0.52489, -0.49573]
                    [t_idx],
                max_relative = 1e-5
            );
        }
    }
}
