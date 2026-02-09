use crate::core::material_properties::{MaterialProperties, GLYCOL25, WATER};
use crate::core::units::{JOULES_PER_KILOWATT_HOUR, LITRES_PER_CUBIC_METRE, MILLIMETRES_IN_METRE, WATTS_PER_KILOWATT};
use crate::input::{PipeworkContents, WaterPipework, WaterPipeworkLocation};
use ordered_float::Pow;
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

const C_N_CALCULATION_FACTOR: f64 = 0.02761; // CIBSE Guide C empirical factor
const C_N_EXPONENT: f64 = 0.8254; // CIBSE Guide C empirical exponent
const N_CALCULATION_FACTOR: f64 = -0.001793; // CIBSE Guide C empirical factor
const N_OFFSET: f64 = 1.245135; // CIBSE Guide C empirical offset

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

    fn length(&self) -> f64;

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
        PI * radius * radius * length * LITRES_PER_CUBIC_METRE as f64
    }
}

impl Pipeworkesque for PipeworkSimple {
    fn location(&self) -> PipeworkLocation {
        self.location
    }

    fn volume(&self) -> f64 {
        self.volume_litres
    }

    fn length(&self) -> f64 {
        self.length_in_m
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
    external_diameter_in_m: f64,
    total_resistance: f64, // in K m / W
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
        let d_insulation_in_m =
            Self::calculate_insulated_diameter(external_diameter_in_m, insulation_thickness);
        let interior_surface_resistance =
            Self::calculate_interior_surface_resistance(internal_diameter_in_m, internal_htc);
        let insulation_resistance = Self::calculate_insulation_resistance(
            insulation_thermal_conductivity,
            external_diameter_in_m,
            d_insulation_in_m,
        );
        let external_surface_resistance =
            Self::calculate_external_surface_resistance(external_htc, d_insulation_in_m);
        let total_resistance =
            interior_surface_resistance + insulation_resistance + external_surface_resistance;

        // Note Python calculates __linear_thermal_transmittance here but it is never used

        Ok(Self {
            simple_pipework,
            external_diameter_in_m,
            total_resistance,
        })
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
    fn calculate_interior_surface_resistance(
        internal_diameter_in_m: f64,
        internal_htc: f64,
    ) -> f64 {
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
    fn calculate_insulation_resistance(
        insulation_thermal_conductivity: f64,
        external_diameter: f64,
        d_insulation_in_m: f64,
    ) -> f64 {
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

    /// Return the c and n value equivalent and thermal mass for the pipework as if working like a emitter
    pub(crate) fn c_n_equivalence(&self) -> (f64, f64, f64) {
        // Calculate c and n values needed to give the same heat loss predicted by CIBSE Guide C
        // using equations 3.101, 3.107 and 3.108. The following equations give a good fit through
        // the data over the realistic working range of temperatures and pipework thicknesses used
        // in dwellings derived from the more detailed CIBSE equations,
        // outputting in the c and n format needed here.
        // The derived equation expects pipework diameter in mm (so unit conversion is required)
        // and outputs c and n coefficients which provide output in W (so conversion to kW required)

        let dop = self.external_diameter_in_m * MILLIMETRES_IN_METRE as f64;
        let c_watts = self.length() * C_N_CALCULATION_FACTOR * dop.pow(C_N_EXPONENT);
        let c = c_watts / WATTS_PER_KILOWATT as f64;
        let n = N_CALCULATION_FACTOR * dop.ln() + N_OFFSET;

        let thermal_mass =
            self.volume() * WATER.specific_heat_capacity() / JOULES_PER_KILOWATT_HOUR as f64;

        // # TODO (from Python): Add thermal mass of pipe itself

        (c, n, thermal_mass)
    }

    /// Return the heat loss from the pipe for the current timestep
    ///
    /// Args:
    ///     inside_temp    -- temperature of water inside the pipe, in degrees C
    ///     outside_temp   -- temperature outside the pipe, in degrees C
    ///
    /// Returns:
    ///     Heat loss in W
    pub(crate) fn calculate_steady_state_heat_loss(
        &self,
        inside_temp: f64,
        outside_temp: f64,
    ) -> f64 {
        (inside_temp - outside_temp) / (self.total_resistance) * self.length()
    }
}

impl Pipeworkesque for Pipework {
    fn location(&self) -> PipeworkLocation {
        self.simple_pipework.location()
    }

    fn volume(&self) -> f64 {
        self.simple_pipework.volume()
    }

    fn length(&self) -> f64 {
        self.simple_pipework.length()
    }

    fn calculate_cool_down_loss(&self, inside_temp: f64, outside_temp: f64) -> f64 {
        self.simple_pipework
            .calculate_cool_down_loss(inside_temp, outside_temp)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::units::WATTS_PER_KILOWATT;
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
    #[case(50.0, 15.0, 5.463755)]
    #[case(51.0, 16.0, 5.463755)]
    #[case(52.0, 17.0, 5.463755)]
    #[case(52.0, 18.0, 5.307648)]
    #[case(51.0, 19.0, 4.995433)]
    #[case(50.0, 20.0, 4.683219)]
    #[case(51.0, 21.0, 4.683219)]
    #[case(52.0, 21.0, 4.839326)]
    fn test_heat_loss(
        pipework: Pipework,
        #[case] t_i: f64,
        #[case] t_o: f64,
        #[case] expected: f64,
    ) {
        let result = pipework.calculate_steady_state_heat_loss(t_i, t_o);
        assert_relative_eq!(result, expected, epsilon = 1e-5);
    }

    #[rstest]
    #[case(50.0, 15.0, 0.01997)]
    #[case(51.0, 16.0, 0.01997)]
    #[case(52.0, 17.0, 0.01997)]
    #[case(52.0, 18.0, 0.01940)]
    #[case(51.0, 19.0, 0.01826)]
    #[case(50.0, 20.0, 0.01712)]
    #[case(51.0, 21.0, 0.01712)]
    #[case(52.0, 21.0, 0.01769)]
    fn test_cool_down_loss(
        pipework: Pipework,
        #[case] t_i: f64,
        #[case] t_o: f64,
        #[case] expected: f64,
    ) {
        let result = pipework.calculate_cool_down_loss(t_i, t_o);
        assert_relative_eq!(result, expected, epsilon = 1e-5);
    }

    #[rstest]
    #[case(2.0, 1.0)]
    #[case(1.0, 1.0)]
    #[case(0.5, 0.4)]
    fn test_invalid_diameter_combinations(
        #[case] internal_diameter: f64,
        #[case] external_diameter: f64,
    ) {
        let result = Pipework::new(
            PipeworkLocation::Internal,
            internal_diameter,
            external_diameter,
            1.0,
            0.035,
            0.038,
            false,
            PipeworkContents::Water,
        );

        assert!(result.is_err());
    }

    // Not porting test_invalid_contents - impossible state in Rust implementation

    // Not porting test_interior_surface_resistance_from_contents or test_external_surface_resistance_from_reflective because
    // they test internal properties which are not accessible in this implementation

    #[rstest]
    fn test_c_n_equivalence(pipework: Pipework) {
        let (c, n, thermal_mass) = pipework.c_n_equivalence();

        let expected_c = 1.0 * 0.02761 * 27.0.pow(0.8254) / WATTS_PER_KILOWATT as f64;
        let expected_n = -0.001793 * 27.0f64.ln() + 1.245135;

        let volume_litres = PI * (0.025 / 2.) * (0.025 / 2.) * 1.0 * LITRES_PER_CUBIC_METRE as f64;
        let expected_thermal_mass =
            volume_litres * WATER.specific_heat_capacity() / JOULES_PER_KILOWATT_HOUR as f64;

        assert_relative_eq!(c, expected_c, epsilon = 1e-5);
        assert_relative_eq!(n, expected_n, epsilon = 1e-5);
        assert_relative_eq!(thermal_mass, expected_thermal_mass, epsilon = 1e-5);

        // Test with a different pipework configuration to ensure formula works correctly
        let pipework2 = Pipework::new(
            PipeworkLocation::External,
            0.050,
            0.055,
            2.0,
            0.040,
            0.025,
            true,
            PipeworkContents::Water,
        )
        .unwrap();

        let (c2, n2, thermal_mass2) = pipework2.c_n_equivalence();
        let expected_c2 = 2.0 * 0.02761 * 55.0.pow(0.8254) / WATTS_PER_KILOWATT as f64;
        let expected_n2 = -0.001793 * 55.0f64.ln() + 1.245135;
        let volume_litres2 = PI * (0.050 / 2.) * (0.050 / 2.) * 2.0 * LITRES_PER_CUBIC_METRE as f64;
        let expected_thermal_mass2 =
            volume_litres2 * WATER.specific_heat_capacity() / JOULES_PER_KILOWATT_HOUR as f64;

        assert_relative_eq!(c2, expected_c2, epsilon = 1e-5);
        assert_relative_eq!(n2, expected_n2, epsilon = 1e-5);
        assert_relative_eq!(thermal_mass2, expected_thermal_mass2, epsilon = 1e-5);
    }

    #[rstest]
    fn test_heat_loss_increases_with_temperature_difference(pipework: Pipework) {
        let loss_small = pipework.calculate_steady_state_heat_loss(50.0, 45.0); // ΔT = 5°C
        let loss_large = pipework.calculate_steady_state_heat_loss(50.0, 20.0); // ΔT = 30°C
        assert!(loss_large > loss_small);
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
    fn test_volume_litres(
        internal_simple_pipe: PipeworkSimple,
        external_simple_pipe: PipeworkSimple,
    ) {
        let volume_ratio = external_simple_pipe.volume() / internal_simple_pipe.volume();
        let expected_ratio = (0.05 / 0.025).pow(2) * (7.0 / 2.5); // (2² × 2.8) = 11.2

        assert_relative_eq!(volume_ratio, expected_ratio, epsilon = 1e-5);
    }

    #[rstest]
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

    #[rstest]
    fn test_cool_down_loss_for_simple_pipework(internal_simple_pipe: PipeworkSimple) {
        // Test case 1: Positive temperature difference (heat loss)
        let loss_positive = internal_simple_pipe.calculate_cool_down_loss(20.0, 10.0);
        assert!(
            loss_positive > 0.,
            "Heat loss should be positive when inside is hotter"
        );

        // Test case 2: Negative temperature difference (heat gain)
        let loss_negative = internal_simple_pipe.calculate_cool_down_loss(25.0, 30.0);
        assert!(
            loss_negative < 0.,
            "Heat gain should be negative when outside is hotter"
        );

        // Test case 3: Zero temperature difference (no heat transfer)
        let loss_zero = internal_simple_pipe.calculate_cool_down_loss(25.0, 25.0);
        assert_relative_eq!(loss_zero, 0., epsilon = 1e-10); // "No heat transfer when temperatures are equal"

        // # Test that energy loss scales with temperature difference
        // Doubling the temperature difference should roughly double the energy loss

        let loss_small = internal_simple_pipe.calculate_cool_down_loss(30.0, 20.0); // 10°C diff
        let loss_large = internal_simple_pipe.calculate_cool_down_loss(40.0, 20.0); // 20°C diff

        assert_relative_eq!(loss_large, 2.0 * loss_small, epsilon = 1e-2) // "Energy loss should scale with temperature difference"
    }
}
