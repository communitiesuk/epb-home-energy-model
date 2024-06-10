use crate::core::material_properties::WATER;
use crate::core::units::{LITRES_PER_CUBIC_METRE, MILLIMETRES_IN_METRE};
use crate::input::{WaterPipeContentsType, WaterPipework};
use std::f64::consts::PI;

// Set default values for the heat transfer coefficients inside the pipe, in W / m^2 K
const INTERNAL_HTC_AIR: f64 = 15.5; // CIBSE Guide C, Table 3.25, air flow rate approx 3 m/s
const INTERNAL_HTC_WATER: f64 = 1500.0; // CIBSE Guide C, Table 3.32

// Set default values for the heat transfer coefficient at the outer surface, in W / m^2 K
const EXTERNAL_REFLECTIVE_HTC: f64 = 5.7; // low emissivity reflective surface, CIBSE Guide C, Table 3.25
const EXTERNAL_NONREFLECTIVE_HTC: f64 = 10.0; // high emissivity non-reflective surface, CIBSE Guide C, Table 3.25

#[derive(Debug)]
pub struct Pipework {
    length_in_m: f64,
    volume_in_litres: f64,
    interior_surface_resistance: f64, // in K m / W
    insulation_resistance: f64,       // in K m / W
    external_surface_resistance: f64, // in K m / W
}

impl From<WaterPipework> for Pipework {
    fn from(input: WaterPipework) -> Self {
        Self::new(
            input.internal_diameter_mm / MILLIMETRES_IN_METRE as f64,
            input.external_diameter_mm / MILLIMETRES_IN_METRE as f64,
            input.length,
            input.insulation_thermal_conductivity,
            input.insulation_thickness_mm / MILLIMETRES_IN_METRE as f64,
            input.surface_reflectivity,
            input.pipe_contents,
        )
    }
}

impl Pipework {
    /// Arguments:
    /// * `internal_diameter` - internal diameter of the pipe, in m
    /// * `external_diameter` - external diameter of the pipe, in m
    /// * `length` - length of pipe, in m
    /// * `k_insulation` - thermal conductivity of the insulation, in W / m K
    /// * `thickness_insulation` - thickness of the pipe insulation, in m
    /// * `reflective` - whether the surface is reflective or not (boolean input)
    /// * `contents` - whether the pipe is carrying air or water
    pub fn new(
        internal_diameter_in_m: f64,
        external_diameter_in_m: f64,
        length_in_m: f64,
        k_insulation: f64,
        thickness_insulation: f64,
        reflective: bool,
        contents: WaterPipeContentsType,
    ) -> Self {
        let volume_in_litres = PI
            * (internal_diameter_in_m / 2f64)
            * (internal_diameter_in_m / 2f64)
            * length_in_m
            * LITRES_PER_CUBIC_METRE as f64;

        // Set the heat transfer coefficient inside the pipe, in W / m^2 K
        let internal_htc = match contents {
            WaterPipeContentsType::Air => INTERNAL_HTC_AIR,
            WaterPipeContentsType::Water => INTERNAL_HTC_WATER,
        };

        // Set the heat transfer coefficient at the outer surface, in W / m^2 K
        let external_htc = if reflective {
            EXTERNAL_REFLECTIVE_HTC
        } else {
            EXTERNAL_NONREFLECTIVE_HTC
        };

        // Calculate the diameter of the pipe including the insulation (D_insulation), in m
        let d_insulation_in_m = external_diameter_in_m + (2f64 * thickness_insulation);

        // Calculate the interior surface resistance, in K m / W
        let interior_surface_resistance = 1f64 / (internal_htc * PI * internal_diameter_in_m);

        // Calculate the insulation resistance, in K m / W
        let insulation_resistance =
            (d_insulation_in_m / internal_diameter_in_m).ln() / (2f64 * PI * k_insulation);

        // Calculate the external surface resistance, in K m / W
        let external_surface_resistance = 1f64 / (external_htc * PI * d_insulation_in_m);

        Self {
            length_in_m,
            volume_in_litres,
            interior_surface_resistance,
            insulation_resistance,
            external_surface_resistance,
        }
    }

    pub fn volume_in_litres(&self) -> f64 {
        self.volume_in_litres
    }

    /// Return the heat loss from the pipe for the current timestep (in W)
    ///
    /// Arguments:
    /// * `inside_temp` - temperature of water (or air) inside the pipe, in degrees C
    /// * `outside_temp` - temperature outside the pipe, in degrees C
    pub fn heat_loss(&self, inside_temp: f64, outside_temp: f64) -> f64 {
        // Calculate total thermal resistance
        let total_resistance = self.interior_surface_resistance
            + self.insulation_resistance
            + self.external_surface_resistance;

        // Calculate the heat loss for the current timestep, in W
        (inside_temp - outside_temp) / total_resistance * self.length_in_m
    }

    #[cfg(test)]
    /// Calculates by how much the temperature of water in a full pipe will fall
    /// over the timestep.
    ///
    /// Arguments:
    /// * `inside_temp` - temperature of water (or air) inside the pipe, in degrees C
    /// * `outside_temp` - temperature outside the pipe, in degrees C
    pub fn temperature_drop(&self, inside_temp: f64, outside_temp: f64) -> f64 {
        use crate::core::units::{JOULES_PER_KILOWATT_HOUR, SECONDS_PER_HOUR, WATTS_PER_KILOWATT};
        let heat_loss_kwh = (SECONDS_PER_HOUR as f64 * self.heat_loss(inside_temp, outside_temp))
            / WATTS_PER_KILOWATT as f64; // heat loss for the one hour timestep in kWh

        *[
            (heat_loss_kwh * JOULES_PER_KILOWATT_HOUR as f64)
                / (WATER.volumetric_heat_capacity() * self.volume_in_litres),
            inside_temp - outside_temp,
        ]
        .iter()
        .min_by(|a, b| a.total_cmp(b))
        .unwrap()
    }

    /// Calculates the total heat loss from a full pipe from demand temp to ambient
    /// temp in kWh
    ///
    /// Arguments:
    /// * `inside_temp` - temperature of water (or air) inside the pipe, in degrees C
    /// * `outside_temp` - temperature outside the pipe, in degrees C
    pub fn cool_down_loss(&self, inside_temp: f64, outside_temp: f64) -> f64 {
        WATER.volumetric_energy_content_kwh_per_litre(inside_temp, outside_temp)
            * self.volume_in_litres
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
    use approx::assert_relative_eq;
    use rstest::*;

    #[fixture]
    pub fn simulation_time() -> SimulationTimeIterator {
        SimulationTime::new(0.0, 8.0, 1.0).iter()
    }

    #[fixture]
    pub fn pipework() -> Pipework {
        Pipework::new(
            0.025,
            0.027,
            1.0,
            0.035,
            0.038,
            false,
            WaterPipeContentsType::Water,
        )
    }

    #[rstest]
    pub fn should_have_correct_interior_surface_resistance(pipework: Pipework) {
        assert_relative_eq!(
            pipework.interior_surface_resistance,
            0.00849,
            max_relative = 0.0005
        );
    }

    #[rstest]
    pub fn should_have_correct_insulation_resistance(pipework: Pipework) {
        assert_relative_eq!(pipework.insulation_resistance, 6.43829, max_relative = 1e-4);
    }

    #[rstest]
    pub fn should_have_correct_external_surface_resistance(pipework: Pipework) {
        assert_relative_eq!(
            pipework.external_surface_resistance,
            0.30904,
            max_relative = 1e-4
        );
    }

    #[rstest]
    pub fn should_have_correct_heat_loss(
        pipework: Pipework,
        simulation_time: SimulationTimeIterator,
    ) {
        let temps_inside = [50.0, 51.0, 52.0, 52.0, 51.0, 50.0, 51.0, 52.0];
        let temps_outside = [15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 21.0];
        let expected_losses = [
            5.18072, 5.18072, 5.18072, 5.03270, 4.73666, 4.44062, 4.44062, 4.58864,
        ];
        for (idx, _) in simulation_time.enumerate() {
            assert_relative_eq!(
                pipework.heat_loss(temps_inside[idx], temps_outside[idx]),
                expected_losses[idx],
                max_relative = 1e-4
            );
        }
    }

    #[rstest]
    pub fn should_have_correct_temp_drop(
        pipework: Pipework,
        simulation_time: SimulationTimeIterator,
    ) {
        let temps_inside = [50.0, 51.0, 52.0, 52.0, 51.0, 50.0, 51.0, 52.0];
        let temps_outside = [15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 21.0];
        let expected_drops = [35.0, 35.0, 35.0, 34.0, 32.0, 30.0, 30.0, 31.0];
        for (idx, _) in simulation_time.enumerate() {
            assert_relative_eq!(
                pipework.temperature_drop(temps_inside[idx], temps_outside[idx]),
                expected_drops[idx],
                max_relative = 0.0005
            );
        }
    }

    #[rstest]
    pub fn should_have_correct_cool_down_loss(
        pipework: Pipework,
        simulation_time: SimulationTimeIterator,
    ) {
        let temps_inside = [50.0, 51.0, 52.0, 52.0, 51.0, 50.0, 51.0, 52.0];
        let temps_outside = [15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 21.0];
        let expected_losses = [
            0.01997, 0.01997, 0.01997, 0.01940, 0.01826, 0.01712, 0.01712, 0.01769,
        ];
        for (idx, _) in simulation_time.enumerate() {
            assert_relative_eq!(
                pipework.cool_down_loss(temps_inside[idx], temps_outside[idx]),
                expected_losses[idx],
                max_relative = 0.0005
            );
        }
    }
}
