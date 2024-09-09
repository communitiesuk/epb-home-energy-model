use crate::core::material_properties::{MaterialProperties, WATER};
use crate::core::units::{LITRES_PER_CUBIC_METRE, MILLIMETRES_IN_METRE};
use crate::input::{WaterPipeContentsType, WaterPipework, WaterPipeworkLocation};
use anyhow::bail;
use std::f64::consts::PI;

// Set default values for the heat transfer coefficients inside the pipe, in W / m^2 K
const INTERNAL_HTC_AIR: f64 = 15.5; // CIBSE Guide C, Table 3.25, air flow rate approx 3 m/s
const INTERNAL_HTC_WATER: f64 = 1500.0; // CIBSE Guide C, Table 3.32 #Note, consider changing to 1478.4

// Set default values for the heat transfer coefficient at the outer surface, in W / m^2 K
const EXTERNAL_REFLECTIVE_HTC: f64 = 5.7; // low emissivity reflective surface, CIBSE Guide C, Table 3.25
const EXTERNAL_NONREFLECTIVE_HTC: f64 = 10.0; // high emissivity non-reflective surface, CIBSE Guide C, Table 3.25

#[derive(Clone, Copy, Debug)]
pub enum PipeworkLocation {
    External,
    Internal,
}

impl From<WaterPipeworkLocation> for PipeworkLocation {
    fn from(value: WaterPipeworkLocation) -> Self {
        match value {
            WaterPipeworkLocation::Internal => PipeworkLocation::Internal,
            WaterPipeworkLocation::External => PipeworkLocation::External,
        }
    }
}

/// A trait to mark that which can be considered "pipeworkesque", analogous to the PipeworkSimple
/// type in the upstream Python, which both the concrete PipeworkSimple and Pipework can be considered
/// to belong to.
pub trait Pipeworkesque {
    fn location(&self) -> PipeworkLocation;

    fn volume_litres(&self) -> f64;

    fn cool_down_loss(&self, inside_temp: f64, outside_temp: f64) -> f64;
}

/// An object to represent heat loss from pipework after flow has stopped
#[derive(Debug)]
pub struct PipeworkSimple {
    location: PipeworkLocation,
    length_in_m: f64,
    volume_litres: f64,
    contents_properties: MaterialProperties,
}

impl PipeworkSimple {
    pub fn new(
        location: PipeworkLocation,
        internal_diameter: f64,
        length: f64,
        contents: WaterPipeContentsType,
    ) -> anyhow::Result<Self> {
        let volume_litres = PI
            * (internal_diameter / 2.)
            * (internal_diameter / 2.)
            * length
            * LITRES_PER_CUBIC_METRE as f64;
        let contents_properties = match contents {
            WaterPipeContentsType::Water => *WATER,
            _ => bail!("No properties available for specified pipe content"),
        };
        Ok(Self {
            location,
            length_in_m: length,
            volume_litres,
            contents_properties,
        })
    }
}

impl Pipeworkesque for PipeworkSimple {
    fn location(&self) -> PipeworkLocation {
        self.location
    }

    fn volume_litres(&self) -> f64 {
        self.volume_litres
    }

    /// Calculates the total heat loss from a full pipe from demand temp to ambient
    /// temp in kWh
    /// Arguments:
    /// * `inside_temp` - temperature of water (or air) inside the pipe, in degrees C
    /// * `outside_temp` - temperature outside the pipe, in degrees C
    fn cool_down_loss(&self, inside_temp: f64, outside_temp: f64) -> f64 {
        self.contents_properties
            .volumetric_energy_content_kwh_per_litre(inside_temp, outside_temp)
            * self.volume_litres
    }
}

#[derive(Debug)]
pub struct Pipework {
    simple_pipework: PipeworkSimple,
    interior_surface_resistance: f64, // in K m / W
    insulation_resistance: f64,       // in K m / W
    external_surface_resistance: f64, // in K m / W
}

impl TryFrom<WaterPipework> for Pipework {
    type Error = anyhow::Error;

    fn try_from(input: WaterPipework) -> anyhow::Result<Self> {
        Self::new(
            input.location.into(),
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
        location: PipeworkLocation,
        internal_diameter_in_m: f64,
        external_diameter_in_m: f64,
        length_in_m: f64,
        k_insulation: f64,
        thickness_insulation: f64,
        reflective: bool,
        contents: WaterPipeContentsType,
    ) -> anyhow::Result<Self> {
        if external_diameter_in_m <= internal_diameter_in_m {
            bail!("Pipework: external diameter must be greater than internal diameter");
        }

        let simple_pipework =
            PipeworkSimple::new(location, internal_diameter_in_m, length_in_m, contents)?;

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

        Ok(Self {
            simple_pipework,
            interior_surface_resistance,
            insulation_resistance,
            external_surface_resistance,
        })
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
        (inside_temp - outside_temp) / total_resistance * self.simple_pipework.length_in_m
    }

    #[cfg(test)]
    /// Calculates by how much the temperature of water in a full pipe will fall
    /// over the timestep.
    ///
    /// Arguments:
    /// * `inside_temp` - temperature of water (or air) inside the pipe, in degrees C
    /// * `outside_temp` - temperature outside the pipe, in degrees C
    pub fn temperature_drop(&self, inside_temp: f64, outside_temp: f64) -> f64 {
        use crate::compare_floats::min_of_2;
        use crate::core::units::{JOULES_PER_KILOWATT_HOUR, SECONDS_PER_HOUR, WATTS_PER_KILOWATT};
        let heat_loss_kwh = (SECONDS_PER_HOUR as f64 * self.heat_loss(inside_temp, outside_temp))
            / WATTS_PER_KILOWATT as f64; // heat loss for the one hour timestep in kWh

        min_of_2(
            (heat_loss_kwh * JOULES_PER_KILOWATT_HOUR as f64)
                / (self
                    .simple_pipework
                    .contents_properties
                    .volumetric_heat_capacity()
                    * self.volume_litres()),
            inside_temp - outside_temp,
        )
    }
}

impl Pipeworkesque for Pipework {
    fn location(&self) -> PipeworkLocation {
        self.simple_pipework.location()
    }

    fn volume_litres(&self) -> f64 {
        self.simple_pipework.volume_litres()
    }

    fn cool_down_loss(&self, inside_temp: f64, outside_temp: f64) -> f64 {
        self.simple_pipework
            .cool_down_loss(inside_temp, outside_temp)
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
            PipeworkLocation::Internal,
            0.025,
            0.027,
            1.0,
            0.035,
            0.038,
            false,
            WaterPipeContentsType::Water,
        )
        .unwrap()
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

    #[fixture]
    fn internal_simple_pipe() -> PipeworkSimple {
        PipeworkSimple::new(
            PipeworkLocation::Internal,
            0.025,
            2.5,
            WaterPipeContentsType::Water,
        )
        .unwrap()
    }

    #[fixture]
    fn external_simple_pipe() -> PipeworkSimple {
        PipeworkSimple::new(
            PipeworkLocation::External,
            0.05,
            7.,
            WaterPipeContentsType::Water,
        )
        .unwrap()
    }

    #[rstest]
    fn test_invalid_contents() {
        assert!(PipeworkSimple::new(
            PipeworkLocation::External,
            0.05,
            7.,
            WaterPipeContentsType::Air
        )
        .is_err());
    }

    #[rstest]
    fn test_get_location(
        internal_simple_pipe: PipeworkSimple,
        external_simple_pipe: PipeworkSimple,
    ) {
        assert!(matches!(
            internal_simple_pipe.location(),
            PipeworkLocation::Internal
        ));
        assert!(matches!(
            external_simple_pipe.location(),
            PipeworkLocation::External
        ));
    }

    #[rstest]
    fn test_get_volume_litres(
        internal_simple_pipe: PipeworkSimple,
        external_simple_pipe: PipeworkSimple,
    ) {
        assert_relative_eq!(
            internal_simple_pipe.volume_litres(),
            1.227184630308513,
            max_relative = 1e-7
        );
        assert_relative_eq!(
            external_simple_pipe.volume_litres(),
            13.744467859455346,
            max_relative = 1e-7
        );
    }

    #[rstest]
    fn test_cool_down_loss(internal_simple_pipe: PipeworkSimple) {
        assert_relative_eq!(
            internal_simple_pipe.cool_down_loss(20.0, 10.0),
            0.014262612481141163,
            max_relative = 1e-7
        );
        assert_relative_eq!(
            internal_simple_pipe.cool_down_loss(25.0, 30.0),
            -0.007131306240570581,
            max_relative = 1e-7
        );
    }
}
