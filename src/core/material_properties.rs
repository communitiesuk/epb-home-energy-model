use crate::core::units::{JOULES_PER_KILOWATT_HOUR, LITRES_PER_CUBIC_METRE};
use std::sync::LazyLock;

/// This module contains data on the properties of materials, and classes to
/// organise this data.

#[derive(Clone, Copy, Debug)]
pub struct MaterialProperties {
    density: f64,                  // kg/litre
    specific_heat_capacity: f64,   // J/(kg.K)
    volumetric_heat_capacity: f64, // J/(litre.K)
}

impl MaterialProperties {
    pub fn new(density: f64, specific_heat_capacity: f64) -> Self {
        Self {
            density,
            specific_heat_capacity,
            volumetric_heat_capacity: specific_heat_capacity * density,
        }
    }

    pub fn density(&self) -> f64 {
        self.density
    }

    pub fn density_kg_per_m3(&self) -> f64 {
        self.density * LITRES_PER_CUBIC_METRE as f64
    }

    pub fn specific_heat_capacity(&self) -> f64 {
        self.specific_heat_capacity
    }

    pub fn specific_heat_capacity_kwh(&self) -> f64 {
        self.specific_heat_capacity / JOULES_PER_KILOWATT_HOUR as f64
    }

    pub fn volumetric_heat_capacity(&self) -> f64 {
        self.volumetric_heat_capacity
    }

    /// Return energy content of material, in J / litre
    ///
    /// Arguments:
    /// * `temp_high` - temperature for which energy content should be calculated, in deg C or K
    /// * `temp_base` - temperature which defines "zero energy", in same units as temp_high
    pub fn volumetric_energy_content_joules_per_litre(
        &self,
        temp_high: f64,
        temp_base: f64,
    ) -> f64 {
        (temp_high - temp_base) * self.volumetric_heat_capacity
    }

    /// Return energy content of material, in kWh / litre
    ///
    /// Arguments:
    /// * `temp_high` - temperature for which energy content should be calculated, in deg C or K
    /// * `temp_base` - temperature which defines "zero energy", in same units as temp_high
    pub fn volumetric_energy_content_kwh_per_litre(&self, temp_high: f64, temp_base: f64) -> f64 {
        self.volumetric_energy_content_joules_per_litre(temp_high, temp_base)
            / JOULES_PER_KILOWATT_HOUR as f64
    }
}

pub static WATER: LazyLock<MaterialProperties> =
    LazyLock::new(|| MaterialProperties::new(1.0, 4184.0));
pub static AIR: LazyLock<MaterialProperties> =
    LazyLock::new(|| MaterialProperties::new(0.001204, 1006.0));

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    pub fn matprop() -> MaterialProperties {
        MaterialProperties::new(1.5, 4184.0)
    }

    #[rstest]
    pub fn should_have_correct_density(matprop: MaterialProperties) {
        assert_eq!(matprop.density(), 1.5, "incorrect density returned");
    }

    #[rstest]
    pub fn should_have_correct_specific_heat_capacity(matprop: MaterialProperties) {
        assert_eq!(
            matprop.specific_heat_capacity(),
            4184.0,
            "incorrect specific heat capacity returned"
        );
    }

    #[rstest]
    pub fn should_have_correct_volumetric_energy_content(matprop: MaterialProperties) {
        assert_eq!(
            matprop.volumetric_heat_capacity(),
            6276.0,
            "incorrect volumetric heat capacity"
        );
    }

    #[rstest]
    pub fn should_provide_correct_volumetric_energy_content(matprop: MaterialProperties) {
        let temp_high = 30.0;
        let temp_low = 20.0;
        assert_eq!(
            matprop.volumetric_energy_content_joules_per_litre(temp_high, temp_low),
            62_760.0,
            "incorrect volumetric energy content (J per litre)"
        );
        assert_relative_eq!(
            matprop.volumetric_energy_content_kwh_per_litre(temp_high, temp_low),
            0.01743333333,
            max_relative = 1e-9
        );
    }
}
