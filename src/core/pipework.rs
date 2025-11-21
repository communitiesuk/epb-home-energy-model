use crate::core::material_properties::{MaterialProperties, GLYCOL25, WATER};
use crate::core::units::{LITRES_PER_CUBIC_METRE, MILLIMETRES_IN_METRE};
use crate::input::{PipeworkContents, WaterPipework, WaterPipeworkLocation};
use std::f64::consts::PI;
use thiserror::Error;

// Set default values for heat transfer coefficients inside pipework, in W / m^2 K
const INTERNAL_HTC_WATER: f64 = 1500.0; // CIBSE Guide C, Table 3.32 #Note, consider changing to 1478.4
const INTERNAL_HTC_GLYCOL25: f64 = INTERNAL_HTC_WATER;
// TODO (from Python) In the absence of a specific figure, use same value for water/glycol mix as for water.
//      Given the figure is relatively high (meaning little resistance to heat flow
//      between the fluid and the inside surface of the pipe) this is unlikely to
//      make a significant difference.

// Set default values for the heat transfer coefficient at the outer surface, in W / m^2 K
const EXTERNAL_REFLECTIVE_HTC: f64 = 5.7; // low emissivity reflective surface, CIBSE Guide C, Table 3.25
const EXTERNAL_NONREFLECTIVE_HTC: f64 = 10.0; // high emissivity non-reflective surface, CIBSE Guide C, Table 3.25

const C_N_CALCULATION_FACTOR: f64 = 0.02716;  // CIBSE Guide C empirical factor
const C_N_EXPONENT: f64 = 0.8254;  // CIBSE Guide C empirical exponent
const N_CALCULATION_FACTOR: f64 = -0.001793;  // CIBSE Guide C empirical factor
const N_OFFSET: f64 = 1.245135;  // CIBSE Guide C empirical offset

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
pub(crate) trait Pipeworkesque {
    fn location(&self) -> PipeworkLocation;

    fn volume(&self) -> f64;

    fn calculate_cool_down_loss(&self, inside_temp: f64, outside_temp: f64) -> f64;
}

/// An object to represent heat loss from pipework after flow has stopped
#[derive(Debug)]
pub(crate) struct PipeworkSimple {
    location: PipeworkLocation,
    length_in_m: f64,
    volume_litres: f64,
    contents_properties: MaterialProperties,
}

impl PipeworkSimple {
    pub(crate) fn new(
        location: PipeworkLocation,
        internal_diameter: f64,
        length: f64,
        contents: PipeworkContents,
    ) -> Result<Self, InvalidPipeworkInput> {
        let volume = Self::calculate_volume(internal_diameter, length);
        let contents_properties = match contents {
            PipeworkContents::Water => *WATER,
            PipeworkContents::Glycol25 => *GLYCOL25,
        };
        Ok(Self {
            location,
            length_in_m: length,
            volume_litres: volume,
            contents_properties,
        })
    }
    
    fn calculate_volume(internal_diameter: f64, length: f64) -> f64 {
        let radius = internal_diameter / 2.;
        return PI
            * radius * radius
            * length
            * LITRES_PER_CUBIC_METRE as f64;
    }
}

impl Pipeworkesque for PipeworkSimple {
    fn location(&self) -> PipeworkLocation {
        self.location
    }

    fn volume(&self) -> f64 {
        self.volume_litres
    }

    /// Calculates the total heat loss from a full pipe from demand temp to ambient
    /// temp in kWh
    /// Arguments:
    /// * `inside_temp` - temperature of water (or air) inside the pipe, in degrees C
    /// * `outside_temp` - temperature outside the pipe, in degrees C
    fn calculate_cool_down_loss(&self, inside_temp: f64, outside_temp: f64) -> f64 {
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
    type Error = InvalidPipeworkInput;

    fn try_from(input: WaterPipework) -> Result<Self, InvalidPipeworkInput> {
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

#[derive(Debug, Error)]
pub enum InvalidPipeworkInput {
    #[error("Pipework: external diameter {external_diameter_in_m}m was not greater than internal diameter {internal_diameter_in_m}m")]
    ExternalNotMoreThanInternal {
        external_diameter_in_m: f64,
        internal_diameter_in_m: f64,
    },
}

impl Pipework {
    /// Arguments:
    /// * `internal_diameter` - internal diameter of the pipe, in m
    /// * `external_diameter` - external diameter of the pipe, in m
    /// * `length` - length of pipe, in m
    /// * `k_insulation` - thermal conductivity of the insulation, in W / m K
    /// * `thickness_insulation` - thickness of the pipe insulation, in m
    /// * `reflective` - whether the surface is reflective or not (boolean input)
    /// * `contents` - whether the pipe is carrying air, water or glycol(25%)/water(75%)
    pub fn new(
        location: PipeworkLocation,
        internal_diameter_in_m: f64,
        external_diameter_in_m: f64,
        length_in_m: f64,
        insulation_thermal_conductivity: f64,
        insulation_thickness: f64,
        reflective: bool,
        contents: PipeworkContents,
    ) -> Result<Self, InvalidPipeworkInput> {
        if external_diameter_in_m <= internal_diameter_in_m {
            return Err(InvalidPipeworkInput::ExternalNotMoreThanInternal {
                external_diameter_in_m,
                internal_diameter_in_m,
            });
        }

        let simple_pipework =
            PipeworkSimple::new(location, internal_diameter_in_m, length_in_m, contents)?;

        let internal_htc = Self::get_internal_heat_transfer_coefficient(contents);
        let external_htc = Self::get_external_heat_transfer_coefficient(reflective);
        let d_insulation_in_m = Self::calculate_insulated_diameter(external_diameter_in_m, insulation_thickness);
        let interior_surface_resistance = Self::calculate_interior_surface_resistance(internal_diameter_in_m, internal_htc);
        let insulation_resistance = Self::calculate_insulation_resistance(insulation_thermal_conductivity, external_diameter_in_m, d_insulation_in_m);
        let external_surface_resistance = Self::calculate_external_surface_resistance(external_htc, d_insulation_in_m);

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

    /// Get the internal heat transfer coefficient based on pipe contents.
    fn get_internal_heat_transfer_coefficient(contents: PipeworkContents) -> f64 {
        match contents {
            PipeworkContents::Water => INTERNAL_HTC_WATER,
            PipeworkContents::Glycol25 => INTERNAL_HTC_GLYCOL25,
        }
    }

    /// Get the external heat transfer coefficient based on pipe reflective surface.
    fn get_external_heat_transfer_coefficient(reflective: bool) -> f64 {
        if reflective {
            EXTERNAL_REFLECTIVE_HTC
        } else {
            EXTERNAL_NONREFLECTIVE_HTC
        }
    }

    /// Calculate the total outer diameter including insulation (D_insulation).
    ///       Args:
    ///           external_diameter_in_m: Outer diameter of the bare pipe, in metres
    ///           thickness_insulation: Thickness of insulation layer, in metres
    ///       Returns:
    ///           Total outer diameter including insulation, in metres
    fn calculate_insulated_diameter(external_diameter_in_m: f64, thickness_insulation: f64) -> f64 {
        external_diameter_in_m + (2.0 * thickness_insulation)
    }

    /// Calculate the interior surface resistance.
    ///
    ///        Args:
    ///            internal_htc: Heat transfer coefficient inside the pipe, in W/m^2K
    ///            internal_diameter: Internal diameter of the pipe, in metres
    ///
    ///        Returns:
    ///            Interior surface resistance, in Km/W
    fn calculate_interior_surface_resistance(internal_diameter_in_m: f64, internal_htc: f64) -> f64 {
        1.0 / (internal_htc * PI * internal_diameter_in_m)
    }

    /// Calculate the insulation resistance.
    ///
    ///        Args:
    ///            insulation_thermal_conductivity: Thermal conductivity of insulation, in W/mK
    ///            external_diameter: External diameter of the pipe, in metres
    ///
    ///        Returns:
    ///            Insulation resistance, in Km/W
    fn calculate_insulation_resistance(insulation_thermal_conductivity: f64, external_diameter: f64, d_insulation_in_m: f64) -> f64 {
        // before (d_insulation_in_m / internal_diameter_in_m).ln() / (2f64 * PI * k_insulation);
        (d_insulation_in_m / external_diameter).ln() / (2. * PI * insulation_thermal_conductivity)
    }

    /// Calculate the external surface resistance.
    ///
    ///        Args:
    ///            external_htc: Heat transfer coefficient outside the pipe, in W/m^2K
    ///            D_insulation: Diameter of the pipe including insulation, in metres
    ///
    ///        Returns:
    ///            External surface resistance, in Km/W
    fn calculate_external_surface_resistance(external_htc: f64, d_insulation_in_m: f64) -> f64 {
        1f64 / (external_htc * PI * d_insulation_in_m)
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
                    * self.volume()),
            inside_temp - outside_temp,
        )
    }
}

impl Pipeworkesque for Pipework {
    fn location(&self) -> PipeworkLocation {
        self.simple_pipework.location()
    }

    fn volume(&self) -> f64 {
        self.simple_pipework.volume()
    }

    fn calculate_cool_down_loss(&self, inside_temp: f64, outside_temp: f64) -> f64 {
        self.simple_pipework
            .calculate_cool_down_loss(inside_temp, outside_temp)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
    use approx::assert_relative_eq;
    use rstest::*;

    #[fixture]
    fn simulation_time() -> SimulationTimeIterator {
        SimulationTime::new(0.0, 8.0, 1.0).iter()
    }

    #[fixture]
    fn pipework() -> Pipework {
        Pipework::new(
            PipeworkLocation::Internal,
            0.025,
            0.027,
            1.0,
            0.035,
            0.038,
            false,
            PipeworkContents::Water,
        )
        .unwrap()
    }

    #[rstest]
    fn should_have_correct_interior_surface_resistance(pipework: Pipework) {
        assert_relative_eq!(
            pipework.interior_surface_resistance,
            0.00849,
            max_relative = 0.0005
        );
    }

    #[rstest]
    fn should_have_correct_insulation_resistance(pipework: Pipework) {
        assert_relative_eq!(pipework.insulation_resistance, 6.08832, max_relative = 1e-4);
    }

    #[rstest]
    fn should_have_correct_external_surface_resistance(pipework: Pipework) {
        assert_relative_eq!(
            pipework.external_surface_resistance,
            0.30904,
            max_relative = 1e-4
        );
    }

    #[ignore = "work in progress - pipework migration"]
    #[rstest]
    fn should_have_correct_heat_loss(pipework: Pipework, simulation_time: SimulationTimeIterator) {
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
    fn should_have_correct_temp_drop(pipework: Pipework, simulation_time: SimulationTimeIterator) {
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
    fn should_have_correct_cool_down_loss(
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
                pipework.calculate_cool_down_loss(temps_inside[idx], temps_outside[idx]),
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
            PipeworkContents::Water,
        )
        .unwrap()
    }

    #[fixture]
    fn external_simple_pipe() -> PipeworkSimple {
        PipeworkSimple::new(
            PipeworkLocation::External,
            0.05,
            7.,
            PipeworkContents::Water,
        )
        .unwrap()
    }

    #[test]
    fn test_glycol() {
        let pipe = PipeworkSimple::new(
            PipeworkLocation::Internal,
            0.05,
            7.,
            PipeworkContents::Glycol25,
        )
        .unwrap();

        assert_eq!(pipe.contents_properties, *GLYCOL25);
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
            internal_simple_pipe.volume(),
            1.227184630308513,
            max_relative = 1e-7
        );
        assert_relative_eq!(
            external_simple_pipe.volume(),
            13.744467859455346,
            max_relative = 1e-7
        );
    }

    #[rstest]
    fn test_cool_down_loss(internal_simple_pipe: PipeworkSimple) {
        assert_relative_eq!(
            internal_simple_pipe.calculate_cool_down_loss(20.0, 10.0),
            0.014262612481141163,
            max_relative = 1e-7
        );
        assert_relative_eq!(
            internal_simple_pipe.calculate_cool_down_loss(25.0, 30.0),
            -0.007131306240570581,
            max_relative = 1e-7
        );
    }

    #[test]
    fn test_invalid_diameter() {
        assert!(Pipework::new(
            PipeworkLocation::Internal,
            2.,
            1.,
            1.0,
            0.035,
            0.038,
            false,
            PipeworkContents::Water,
        )
        .is_err());
    }

    #[test]
    fn test_interior_surface_resistance_from_contents() {
        let pipework = Pipework::new(
            PipeworkLocation::Internal,
            1.,
            2.,
            1.,
            1.,
            1.,
            false,
            PipeworkContents::Water,
        )
        .unwrap();

        assert_eq!(
            pipework.interior_surface_resistance,
            1. / (INTERNAL_HTC_WATER * PI)
        );

        let pipework = Pipework::new(
            PipeworkLocation::Internal,
            1.,
            2.,
            1.,
            1.,
            1.,
            false,
            PipeworkContents::Glycol25,
        )
        .unwrap();

        assert_eq!(
            pipework.interior_surface_resistance,
            1. / (INTERNAL_HTC_GLYCOL25 * PI)
        );
    }

    #[test]
    fn test_external_surface_resistance_from_reflective() {
        let pipework = Pipework::new(
            PipeworkLocation::Internal,
            1.,
            2.,
            1.,
            1.,
            1.,
            false,
            PipeworkContents::Water,
        )
        .unwrap();

        assert_eq!(
            pipework.external_surface_resistance,
            1. / (EXTERNAL_NONREFLECTIVE_HTC * PI * 4.)
        );

        let pipework = Pipework::new(
            PipeworkLocation::Internal,
            1.,
            2.,
            1.,
            1.,
            1.,
            true,
            PipeworkContents::Water,
        )
        .unwrap();

        assert_eq!(
            pipework.external_surface_resistance,
            1. / (EXTERNAL_REFLECTIVE_HTC * PI * 4.)
        );
    }
}
