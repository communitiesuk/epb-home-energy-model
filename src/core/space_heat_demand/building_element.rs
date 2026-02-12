use crate::bail;
use crate::core::controls::time_control::{Control, ControlBehaviour};
use crate::core::units::{
    average_monthly_to_annual, calculate_thermal_resistance_of_virtual_layer, Orientation360,
    JOULES_PER_KILOJOULE,
};
use crate::corpus::Controls;
use crate::external_conditions::{
    CalculatedDirectDiffuseTotalIrradiance, ExternalConditions, WindowShadingObject,
};
use crate::input::{
    EdgeInsulation, FloorData, MassDistributionClass, WindShieldLocation,
    WindowTreatment as WindowTreatmentInput, WindowTreatmentControl as WindowTreatmentControlInput,
    WindowTreatmentType, PITCH_LIMIT_HORIZ_CEILING, PITCH_LIMIT_HORIZ_FLOOR,
};
use crate::simulation_time::SimulationTimeIteration;
use anyhow::anyhow;
use atomic_float::AtomicF64;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

// Difference between external air temperature and sky temperature
// (default value for intermediate climatic region from BS EN ISO 52016-1:2017, Table B.19)
const TEMP_DIFF_SKY: f64 = 11.0; // Kelvin

/// calc the vertically projected height of a surface from
/// the actual height and tilt of the surface
pub(crate) fn projected_height(tilt: f64, height: f64) -> f64 {
    let mut ph = height * tilt.to_radians().sin();
    // BS EN ISO 52010-1 Table 7 geometric input data; shading. Footnote d
    // validity interval H1;ic > 0
    // if horizontal (height = 0): choose small value e.g. H1 = 0.01 m"""
    if ph < 0.01 {
        ph = 0.01;
    }

    ph
}

#[derive(PartialEq, Clone, Debug)]
enum PartyWallCavityType {
    Solid,             // Solid wall or structurally insulated panel
    UnfilledUnsealed,  // Unfilled cavity with no effective edge sealing
    UnfilledSealed,    // Unfilled cavity with effective sealing
    FilledSealed,      // Fully filled cavity with effective sealing
    FilledUnsealed,    // Fully filled cavity with no effective edge sealing
    DefinedResistance, // User-defined thermal resistance
}

#[derive(PartialEq, Clone, Debug)]
enum PartyWallLiningType {
    WetPlaster,
    DryLined,
}

fn calculate_area(height: f64, width: f64) -> f64 {
    height * width
}

/// Calculate the effective thermal resistance of the party wall cavity
///
/// For defined_resistance type, uses the provided value directly.
/// For other types, derives the cavity resistance from equivalent U-values
/// which are based on research into heat loss via air movement in party wall cavities.
///
/// The derivation assumes typical party wall constructions and calculates
/// what cavity resistance (R_cavity) is needed to achieve the equivalent
/// U-values (initially based on SAP guidance):
/// - Unsealed cavity, dry lined: U = 0.6 W/m².K → R_cavity ≈ 1.2 m².K/W
/// - Unsealed cavity, wet plaster: U = 0.2 W/m².K → R_cavity ≈ 4.5 m².K/W
/// - Edge-sealed cavity: U = 0.2 W/m².K → R_cavity ≈ 4.5 m².K/W
/// - Unsealed filled cavity: U = 0.2 W/m².K → R_cavity ≈ 4.5 m².K/W
/// - Solid or filled and edge-sealed: U = 0.0 W/m².K → R_cavity = 999999 m².K/W (effectively adiabatic)
///
/// Arguments:
/// * `party_wall_cavity_type` - type of party wall cavity construction
/// * `party_wall_lining_type` - type of party wall lining material
/// * `thermal_resistance_cavity` - user-defined cavity resistance (only for 'defined_resistance' type)
///
/// Returns:
/// * `r_cavity` - effective thermal resistance of the cavity, in m2.K / W
fn calculate_cavity_resistance(
    party_wall_cavity_type: &PartyWallCavityType,
    party_wall_lining_type: &Option<PartyWallLiningType>,
    thermal_resistance_cavity: Option<f64>,
) -> anyhow::Result<f64> {
    if *party_wall_cavity_type == PartyWallCavityType::DefinedResistance {
        return thermal_resistance_cavity.ok_or_else(|| {
            anyhow!(
                "thermal_resistance_cavity is validated by schema to be \
                not None for DEFINED_RESISTANCE type"
            )
        });
    }

    // Thermal resistance equivalents for deriving cavity resistance
    let thermal_resistance = match (party_wall_lining_type, party_wall_cavity_type) {
        (None, PartyWallCavityType::Solid) => 999999.,
        (None, PartyWallCavityType::FilledSealed) => 999999.,
        (None, _) => {
            bail!("invalid combination of party wall cavity type and party wall lining type")
        }
        (Some(PartyWallLiningType::DryLined), PartyWallCavityType::UnfilledUnsealed) => 1.2,
        (Some(PartyWallLiningType::DryLined), PartyWallCavityType::FilledUnsealed) => 4.5,
        (Some(PartyWallLiningType::DryLined), PartyWallCavityType::UnfilledSealed) => 4.5,
        (Some(PartyWallLiningType::DryLined), _) => {
            bail!("invalid combination of party wall cavity type and party wall lining type")
        }
        (Some(PartyWallLiningType::WetPlaster), PartyWallCavityType::UnfilledUnsealed) => 4.5,
        (Some(PartyWallLiningType::WetPlaster), PartyWallCavityType::FilledUnsealed) => 4.5,
        (Some(PartyWallLiningType::WetPlaster), PartyWallCavityType::UnfilledSealed) => 4.5,
        (Some(PartyWallLiningType::WetPlaster), _) => {
            bail!("invalid combination of party wall cavity type and party wall lining type")
        }
    };
    Ok(thermal_resistance)
}

#[derive(Debug, PartialEq)]
pub enum HeatFlowDirection {
    Horizontal,
    Upwards,
    Downwards,
}

// The different variations for each of these are:
// - Heat transfer with internal environment:
//    - Common
// - Heat transfer through and storage within the element:
//    - 2 nodes
//    - 5 nodes
//    - 3+2 nodes
// - Heat transfer with environment on other side of element:
//    - Ground
//    - Conditioned space
//    - Unconditioned space
//    - Outside
// - Interaction with solar radiation:
//    - Absorbed
//    - Transmitted
//    - Not exposed
//
// The different building element types are composed of the following combinations:
// - BuildingElementOpaque: Common, 5 nodes, Outside, Absorbed
// - BuildingElementTransparent: Common, 2 nodes, Outside, Transmitted
// - BuildingElementAdjacentConditionedSpace: Common, 5 nodes, Conditioned space, Not exposed
// - BuildingElementAdjacentUnconditionedSpaceSimple: Common, 5 nodes, Unconditioned space, Not exposed
// - BuildingElementPartyWall: Common, 5 nodes, Unconditioned space, Not exposed
// - BuildingElementGround: Common, 3+2 nodes, Ground, Not exposed
//
// BuildingElementTransparent also has the functions projected_height, mid_height and orientation, which I think are now
// unused. If this is confirmed to be the case, then these can be deleted.
//
// BuildingElementGround has several alternative sets of inputs, depending on the floor_type input. The current class
// requires all of these inputs and ignores the ones that are not relevant. This class could arguably be split into 5
// different classes (one for each floor_type option) to avoid this. The alternative sets of inputs are for calculating
// self.__h_pi and self.__h_pe, which are ultimately used in the temp_ext function. Therefore, one way to handle this would
// be to have 5 different classes, calculate h_pi and h_pe, then feed these into the constructor of the "Ground" object
// under the category "Heat transfer with environment on other side of element".
//
// The relevant functions for each set of characteristics are:
// - Heat transfer with internal environment:
//    - heat_flow_direction
//    - r_si
//    - h_ci
//    - h_ri
//    - pitch_class
// - Heat transfer through and storage within the element:
//    - no_of_nodes
//    - no_of_inside_nodes
//    - init_h_pli
//    - init_k_pli
// - Heat transfer with environment on other side of element:
//    - r_se
//    - h_ce
//    - h_re
//    - temp_ext
// - Interaction with solar radiation:
//    - i_sol_dir_dif
//    - solar_gains
//    - shading_factors_direct_diffuse

//Values from BS EN ISO 13789:2017, Table 8: Conventional surface heat
//transfer coefficients
const H_CI_UPWARDS: f64 = 5.0;
const H_CI_HORIZONTAL: f64 = 2.5;
const H_CI_DOWNWARDS: f64 = 0.7;
pub(crate) const H_CE: f64 = 20.0;
const H_RI: f64 = 5.13;
pub(crate) const H_RE: f64 = 4.14;

// Surface resistances of building elements, in m2 K / W
pub(crate) const R_SI_HORIZONTAL: f64 = 1.0 / (H_RI + H_CI_HORIZONTAL);
pub(crate) const R_SI_UPWARDS: f64 = 1.0 / (H_RI + H_CI_UPWARDS);
pub(crate) const R_SI_DOWNWARDS: f64 = 1.0 / (H_RI + H_CI_DOWNWARDS);
const R_SE: f64 = 1.0 / (H_CE + H_RE);

#[derive(Debug)]
pub(crate) enum BuildingElement {
    Opaque(BuildingElementOpaque),
    AdjacentConditionedSpace(BuildingElementAdjacentConditionedSpace),
    AdjacentUnconditionedSpaceSimple(BuildingElementAdjacentUnconditionedSpaceSimple),
    Ground(BuildingElementGround),
    Transparent(BuildingElementTransparent),
    PartyWall(BuildingElementPartyWall),
}

impl BuildingElement {
    pub(crate) fn area(&self) -> f64 {
        match self {
            BuildingElement::Opaque(el) => el.area(),
            BuildingElement::AdjacentConditionedSpace(el) => el.area(),
            BuildingElement::AdjacentUnconditionedSpaceSimple(el) => el.area(),
            BuildingElement::Ground(el) => el.area(),
            BuildingElement::Transparent(el) => el.area(),
            BuildingElement::PartyWall(el) => el.area(),
        }
    }

    fn as_heat_transfer_through(&self) -> &dyn HeatTransferThrough {
        match self {
            BuildingElement::Opaque(el) => el,
            BuildingElement::AdjacentConditionedSpace(el) => el,
            BuildingElement::AdjacentUnconditionedSpaceSimple(el) => el,
            BuildingElement::Ground(el) => el,
            BuildingElement::Transparent(el) => el,
            BuildingElement::PartyWall(el) => el,
        }
    }

    /// Return number of nodes including external and internal layers
    pub(crate) fn number_of_nodes(&self) -> usize {
        self.as_heat_transfer_through().number_of_nodes()
    }

    pub(crate) fn k_pli(&self) -> &[f64] {
        self.as_heat_transfer_through().k_pli()
    }

    pub(crate) fn solar_gains(&self, simtime: SimulationTimeIteration) -> anyhow::Result<f64> {
        match self {
            BuildingElement::Opaque(el) => Ok(el.solar_gains()),
            BuildingElement::AdjacentConditionedSpace(el) => Ok(el.solar_gains()),
            BuildingElement::AdjacentUnconditionedSpaceSimple(el) => Ok(el.solar_gains()),
            BuildingElement::Ground(el) => Ok(el.solar_gains()),
            BuildingElement::Transparent(el) => el.solar_gains(simtime),
            BuildingElement::PartyWall(el) => Ok(el.solar_gains()),
        }
    }

    pub(crate) fn h_ce(&self) -> f64 {
        match self {
            BuildingElement::Opaque(el) => el.h_ce(),
            BuildingElement::AdjacentConditionedSpace(el) => el.h_ce(),
            BuildingElement::AdjacentUnconditionedSpaceSimple(el) => el.h_ce(),
            BuildingElement::Ground(el) => el.h_ce(),
            BuildingElement::Transparent(el) => el.h_ce(),
            BuildingElement::PartyWall(el) => el.h_ce(),
        }
    }

    pub(crate) fn h_re(&self) -> f64 {
        match self {
            BuildingElement::Opaque(el) => el.h_re(),
            BuildingElement::AdjacentConditionedSpace(el) => el.h_re(),
            BuildingElement::AdjacentUnconditionedSpaceSimple(el) => el.h_re(),
            BuildingElement::Ground(el) => el.h_re(),
            BuildingElement::Transparent(el) => el.h_re(),
            BuildingElement::PartyWall(el) => el.h_re(),
        }
    }

    pub(crate) fn h_ri(&self) -> f64 {
        match self {
            BuildingElement::Opaque(el) => el.h_ri(),
            BuildingElement::AdjacentConditionedSpace(el) => el.h_ri(),
            BuildingElement::AdjacentUnconditionedSpaceSimple(el) => el.h_ri(),
            BuildingElement::Ground(el) => el.h_ri(),
            BuildingElement::Transparent(el) => el.h_ri(),
            BuildingElement::PartyWall(el) => el.h_ri(),
        }
    }

    pub(crate) fn solar_absorption_coeff(&self) -> f64 {
        match self {
            BuildingElement::Opaque(el) => el.solar_absorption_coeff(),
            BuildingElement::AdjacentConditionedSpace(el) => el.solar_absorption_coeff(),
            BuildingElement::AdjacentUnconditionedSpaceSimple(el) => el.solar_absorption_coeff(),
            BuildingElement::Ground(el) => el.solar_absorption_coeff(),
            BuildingElement::Transparent(el) => el.solar_absorption_coeff(),
            BuildingElement::PartyWall(el) => el.solar_absorption_coeff(),
        }
    }

    pub(crate) fn therm_rad_to_sky(&self) -> f64 {
        match self {
            BuildingElement::Opaque(el) => el.therm_rad_to_sky(),
            BuildingElement::AdjacentConditionedSpace(el) => el.therm_rad_to_sky(),
            BuildingElement::AdjacentUnconditionedSpaceSimple(el) => el.therm_rad_to_sky(),
            BuildingElement::Ground(el) => el.therm_rad_to_sky(),
            BuildingElement::Transparent(el) => el.therm_rad_to_sky(),
            BuildingElement::PartyWall(el) => el.therm_rad_to_sky(),
        }
    }

    pub(crate) fn h_pli_by_index_unchecked(
        &self,
        idx: usize,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        match self {
            BuildingElement::Opaque(el) => Ok(el.h_pli()[idx]),
            BuildingElement::AdjacentConditionedSpace(el) => Ok(el.h_pli()[idx]),
            BuildingElement::AdjacentUnconditionedSpaceSimple(el) => Ok(el.h_pli()[idx]),
            BuildingElement::Ground(el) => Ok(el.h_pli()[idx]),
            BuildingElement::Transparent(el) => el.h_pli_by_index(idx, simulation_time_iteration),
            BuildingElement::PartyWall(el) => Ok(el.h_pli()[idx]),
        }
    }

    pub(crate) fn i_sol_dir_dif(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        Ok(match self {
            BuildingElement::Opaque(el) => el.i_sol_dir_dif(simtime)?,
            BuildingElement::AdjacentConditionedSpace(el) => el.i_sol_dir_dif(simtime),
            BuildingElement::AdjacentUnconditionedSpaceSimple(el) => el.i_sol_dir_dif(simtime),
            BuildingElement::Ground(el) => el.i_sol_dir_dif(simtime),
            BuildingElement::Transparent(el) => el.i_sol_dir_dif(simtime),
            BuildingElement::PartyWall(el) => el.i_sol_dir_dif(simtime),
        })
    }

    pub(crate) fn shading_factors_direct_diffuse(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        match self {
            BuildingElement::Opaque(el) => el.shading_factors_direct_diffuse(simtime),
            BuildingElement::AdjacentConditionedSpace(el) => {
                Ok(el.shading_factors_direct_diffuse(simtime))
            }
            BuildingElement::AdjacentUnconditionedSpaceSimple(el) => {
                Ok(el.shading_factors_direct_diffuse(simtime))
            }
            BuildingElement::Ground(el) => Ok(el.shading_factors_direct_diffuse(simtime)),
            BuildingElement::Transparent(el) => el.shading_factors_direct_diffuse(simtime),
            BuildingElement::PartyWall(el) => Ok(el.shading_factors_direct_diffuse(simtime)),
        }
    }

    pub(crate) fn temp_ext(&self, simtime: SimulationTimeIteration) -> f64 {
        match self {
            BuildingElement::Opaque(el) => el.temp_ext(simtime),
            BuildingElement::AdjacentConditionedSpace(el) => el.temp_ext(simtime),
            BuildingElement::AdjacentUnconditionedSpaceSimple(el) => el.temp_ext(simtime),
            BuildingElement::Ground(el) => el.temp_ext(simtime),
            BuildingElement::Transparent(el) => el.temp_ext(simtime),
            BuildingElement::PartyWall(el) => el.temp_ext(simtime),
        }
    }

    /// Return number of nodes excluding external and internal layers
    pub(crate) fn number_of_inside_nodes(&self) -> usize {
        self.number_of_nodes() - 2
    }

    pub(crate) fn h_ci(&self, temp_int_air: f64, temp_int_surface: f64) -> f64 {
        match self {
            BuildingElement::Opaque(el) => el.h_ci(temp_int_air, temp_int_surface),
            BuildingElement::AdjacentConditionedSpace(el) => {
                el.h_ci(temp_int_air, temp_int_surface)
            }
            BuildingElement::AdjacentUnconditionedSpaceSimple(el) => {
                el.h_ci(temp_int_air, temp_int_surface)
            }
            BuildingElement::Ground(el) => el.h_ci(temp_int_air, temp_int_surface),
            BuildingElement::Transparent(el) => el.h_ci(temp_int_air, temp_int_surface),
            BuildingElement::PartyWall(el) => el.h_ci(temp_int_air, temp_int_surface),
        }
    }

    pub(crate) fn fabric_heat_loss(&self) -> f64 {
        match self {
            BuildingElement::Opaque(el) => el.fabric_heat_loss(),
            BuildingElement::AdjacentConditionedSpace(el) => el.fabric_heat_loss(),
            BuildingElement::AdjacentUnconditionedSpaceSimple(el) => el.fabric_heat_loss(),
            BuildingElement::Ground(el) => el.fabric_heat_loss(),
            BuildingElement::Transparent(el) => el.fabric_heat_loss(),
            BuildingElement::PartyWall(el) => el.fabric_heat_loss(),
        }
    }

    pub(crate) fn heat_capacity(&self) -> f64 {
        match self {
            BuildingElement::Opaque(el) => el.heat_capacity(),
            BuildingElement::AdjacentConditionedSpace(el) => el.heat_capacity(),
            BuildingElement::AdjacentUnconditionedSpaceSimple(el) => el.heat_capacity(),
            BuildingElement::Ground(el) => el.heat_capacity(),
            BuildingElement::Transparent(el) => el.heat_capacity(),
            BuildingElement::PartyWall(el) => el.heat_capacity(),
        }
    }
}

pub(crate) trait HeatTransferInternal {
    /// Determine direction of heat flow for a surface
    fn heat_flow_direction(&self, temp_int_air: f64, temp_int_surface: f64) -> HeatFlowDirection {
        let pitch = self.pitch();
        if (PITCH_LIMIT_HORIZ_CEILING..=PITCH_LIMIT_HORIZ_FLOOR).contains(&pitch) {
            HeatFlowDirection::Horizontal
        } else {
            let inwards_heat_flow = temp_int_air < temp_int_surface;
            let is_floor = pitch > PITCH_LIMIT_HORIZ_FLOOR;
            let is_ceiling = pitch < PITCH_LIMIT_HORIZ_CEILING;
            let upwards_heat_flow =
                (is_floor && inwards_heat_flow) || (is_ceiling && !inwards_heat_flow);
            if upwards_heat_flow {
                HeatFlowDirection::Upwards
            } else {
                HeatFlowDirection::Downwards
            }
        }
    }

    fn pitch(&self) -> f64;

    /// Return internal surface resistance, in m2 K / W
    fn r_si(&self) -> f64 {
        self.r_si_with_pitch(self.pitch())
    }

    fn r_si_with_pitch(&self, pitch: f64) -> f64 {
        match pitch {
            PITCH_LIMIT_HORIZ_CEILING..=PITCH_LIMIT_HORIZ_FLOOR => R_SI_HORIZONTAL,
            ..PITCH_LIMIT_HORIZ_CEILING => R_SI_UPWARDS,
            PITCH_LIMIT_HORIZ_FLOOR.. => R_SI_DOWNWARDS,
            _ => unreachable!("Rust cannot tell that above is exhaustive"),
        }
    }

    /// Return internal convective heat transfer coefficient, in W / (m2.K)
    fn h_ci(&self, temp_int_air: f64, temp_int_surface: f64) -> f64 {
        match self.heat_flow_direction(temp_int_air, temp_int_surface) {
            HeatFlowDirection::Horizontal => H_CI_HORIZONTAL,
            HeatFlowDirection::Upwards => H_CI_UPWARDS,
            HeatFlowDirection::Downwards => H_CI_DOWNWARDS,
        }
    }

    /// Return internal radiative heat transfer coefficient, in W / (m2.K)
    fn h_ri(&self) -> f64 {
        H_RI
    }

    /// Convert U-value from input data to thermal resistance of construction only
    /// (not incl. surface resistances)
    #[cfg(test)]
    fn convert_uvalue_to_resistance(&self, u_value: f64, pitch: f64) -> f64 {
        (1.0 / u_value) - self.r_si_with_pitch(pitch) - R_SE
    }

    #[cfg(test)]
    fn pitch_class(&self, pitch: f64) -> HeatFlowDirection {
        match pitch {
            PITCH_LIMIT_HORIZ_CEILING..=PITCH_LIMIT_HORIZ_FLOOR => HeatFlowDirection::Horizontal,
            ..PITCH_LIMIT_HORIZ_CEILING => HeatFlowDirection::Upwards,
            PITCH_LIMIT_HORIZ_FLOOR.. => HeatFlowDirection::Downwards,
            _ => unreachable!("Rust cannot tell that above is exhaustive"),
        }
    }
}

#[allow(dead_code)]
pub(crate) trait HeatTransferInternalCommon: HeatTransferInternal {}

pub(crate) trait HeatTransferThrough {
    fn init_heat_transfer_through(
        &mut self,
        thermal_resistance_construction: f64,
        areal_heat_capacity: f64,
    ) {
        self.set_r_c(thermal_resistance_construction);
        self.set_k_m(areal_heat_capacity);
    }

    fn set_r_c(&mut self, thermal_resistance_construction: f64);
    fn r_c(&self) -> f64;
    fn k_m(&self) -> f64;
    fn set_k_m(&mut self, areal_heat_capacity: f64);
    fn k_pli(&self) -> &[f64];

    fn number_of_nodes(&self) -> usize {
        self.k_pli().len()
    }

    #[cfg(test)]
    fn number_of_inside_nodes(&self) -> usize {
        self.number_of_nodes() - 2
    }

    fn r_se(&self) -> f64;
    fn r_si(&self) -> f64;

    fn area(&self) -> f64;

    /// Return the fabric heat loss for the building element
    fn fabric_heat_loss(&self) -> f64 {
        let u_value = 1.0 / (self.r_c() + self.r_se() + self.r_si());
        self.area() * u_value
    }

    /// Return the fabric heat capacity for the building element
    fn heat_capacity(&self) -> f64 {
        self.area() * (self.k_m() / JOULES_PER_KILOJOULE as f64)
    }

    fn h_pli(&self) -> &[f64];
}

pub(crate) trait HeatTransferThrough2Nodes: HeatTransferThrough {
    fn init_heat_transfer_through_2_nodes(&mut self, thermal_resistance_construction: f64) {
        self.init_heat_transfer_through(thermal_resistance_construction, 0.);
        self.set_h_pli(self.init_h_pli(thermal_resistance_construction));
        self.set_k_pli(self.init_k_pli());
    }

    fn set_h_pli(&mut self, h_pli: [f64; 1]);
    fn set_k_pli(&mut self, k_pli: [f64; 2]);

    // Calculate node conductances (h_pli) and node heat capacities (k_pli)
    // according to BS EN ISO 52016-1:2017, section 6.5.7.4

    fn init_h_pli(&self, thermal_resistance_construction: f64) -> [f64; 1] {
        [1.0 / thermal_resistance_construction]
    }

    fn init_k_pli(&self) -> [f64; 2] {
        [0.0, 0.0]
    }
}

pub(crate) trait HeatTransferThrough5Nodes: HeatTransferThrough {
    fn init_heat_transfer_through_5_nodes(
        &mut self,
        thermal_resistance_construction: f64,
        mass_distribution: MassDistributionClass,
        areal_heat_capacity: f64,
    ) {
        self.set_r_c(thermal_resistance_construction);
        self.set_k_m(areal_heat_capacity);
        self.set_h_pli(self.init_h_pli(thermal_resistance_construction));
        self.set_k_pli(self.init_k_pli(mass_distribution, areal_heat_capacity));
    }

    fn set_h_pli(&mut self, h_pli: [f64; 4]);
    fn set_k_pli(&mut self, k_pli: [f64; 5]);

    fn init_h_pli(&self, thermal_resistance_construction: f64) -> [f64; 4] {
        let h_outer = 6.0 / thermal_resistance_construction;
        let h_inner = 3.0 / thermal_resistance_construction;
        [h_outer, h_inner, h_inner, h_outer]
    }

    fn init_k_pli(
        &self,
        mass_distribution: MassDistributionClass,
        areal_heat_capacity: f64,
    ) -> [f64; 5] {
        match mass_distribution {
            MassDistributionClass::I => [0.0, 0.0, 0.0, 0.0, areal_heat_capacity],
            MassDistributionClass::E => [areal_heat_capacity, 0.0, 0.0, 0.0, 0.0],
            MassDistributionClass::IE => {
                let k_ie = areal_heat_capacity / 2.0;
                [k_ie, 0.0, 0.0, 0.0, k_ie]
            }
            MassDistributionClass::D => {
                let k_inner = areal_heat_capacity / 4.0;
                let k_outer = areal_heat_capacity / 8.0;
                [k_outer, k_inner, k_inner, k_inner, k_outer]
            }
            MassDistributionClass::M => [0.0, 0.0, areal_heat_capacity, 0.0, 0.0],
        }
    }
}

pub(crate) trait HeatTransferThrough3Plus2Nodes: HeatTransferThrough {
    fn init_heat_transfer_through_3_plus_2_nodes(
        &mut self,
        thermal_resistance_floor_construction: f64,
        r_gr: f64,
        mass_distribution: MassDistributionClass,
        k_gr: f64,
        areal_heat_capacity: f64,
    ) {
        self.set_k_m(areal_heat_capacity);
        // Calculate node conductances (h_pli) and node heat capacities (k_pli)
        // according to BS EN ISO 52016-1:2017, section 6.5.7.4
        self.set_h_pli(self.init_h_pli(thermal_resistance_floor_construction, r_gr));
        self.set_k_pli(self.init_k_pli(mass_distribution, k_gr, areal_heat_capacity));
    }

    fn set_h_pli(&mut self, h_pli: [f64; 4]);
    fn set_k_pli(&mut self, k_pli: [f64; 5]);

    fn init_h_pli(&self, thermal_resistance_floor_construction: f64, r_gr: f64) -> [f64; 4] {
        // BS EN ISO 52016:2017 states that the thermal_resistance_construction (resistance including the
        //         effect of the ground) should be used in the equations below. However,
        //         this leads to double-counting of r_si, r_gr and r_vi as these are already
        //         accounted for separately, so we have used thermal_resistance_floor_construction (resistance of the floor
        //         construction only) here instead
        let h_4 = 4.0 / thermal_resistance_floor_construction;
        let h_3 = 2.0 / thermal_resistance_floor_construction;
        let h_2 = 1.0 / (thermal_resistance_floor_construction / 4. + r_gr / 2.);
        let h_1 = 2.0 / r_gr;
        [h_1, h_2, h_3, h_4]
    }

    fn init_k_pli(
        &self,
        mass_distribution: MassDistributionClass,
        k_gr: f64,
        areal_heat_capacity: f64,
    ) -> [f64; 5] {
        match mass_distribution {
            MassDistributionClass::I => [0.0, k_gr, 0.0, 0.0, areal_heat_capacity],
            MassDistributionClass::E => [0.0, k_gr, areal_heat_capacity, 0.0, 0.0],
            MassDistributionClass::IE => {
                let k_ie = areal_heat_capacity / 2.0;
                [0.0, k_gr, k_ie, 0.0, k_ie]
            }
            MassDistributionClass::D => {
                let k_inner = areal_heat_capacity / 2.0;
                let k_outer = areal_heat_capacity / 4.0;
                [0.0, k_gr, k_outer, k_inner, k_outer]
            }
            MassDistributionClass::M => [0.0, k_gr, 0.0, areal_heat_capacity, 0.0],
        }
    }
}

pub(crate) trait HeatTransferOtherSide {
    fn init_heat_transfer_other_side(&mut self, f_sky: Option<f64>) {
        self.init_super(f_sky)
    }

    // The "_super" suffix naming here is to enable similar functionality to Python calling e.g. super.__init()
    // It's impossible to do this in Rust as implementations in subtraits simply override but we can work round
    // by using this kind of naming, even if it pollutes the interface of any implementing structs.
    // We could use "sealed traits" for this, but likely not worth it.
    fn init_super(&mut self, f_sky: Option<f64>) {
        let f_sky = f_sky.unwrap_or(0.0);
        self.set_f_sky(f_sky);
        self.set_therm_rad_to_sky(f_sky * self.h_re() * TEMP_DIFF_SKY);
    }

    fn set_f_sky(&mut self, f_sky: f64);
    fn set_therm_rad_to_sky(&mut self, therm_rad_to_sky: f64);
    fn therm_rad_to_sky(&self) -> f64;
    #[allow(dead_code)]
    fn f_sky(&self) -> f64;

    /// Return external surface resistance, in m2 K / W
    fn r_se(&self) -> f64 {
        R_SE
    }

    /// Return external convective heat transfer coefficient, in W / (m2.K)
    fn h_ce(&self) -> f64 {
        self.h_ce_super()
    }

    // See note against init_super above re purpose of this "_super" suffix here
    fn h_ce_super(&self) -> f64 {
        H_CE
    }

    fn h_re(&self) -> f64 {
        self.h_re_super()
    }

    // See note against init_super above re purpose of this "_super" suffix here
    fn h_re_super(&self) -> f64 {
        H_RE
    }

    fn temp_ext(&self, simtime: SimulationTimeIteration) -> f64 {
        self.external_conditions().air_temp(&simtime)
    }

    fn external_conditions(&self) -> &ExternalConditions;

    // there is a method here in the Python for fabric_heat_loss() that intentionally blows up
    // just eliding implementing this here as we should enforce at compile time that it is never called
}

// Assume values for temp_int_annual and temp_int_monthly
// These are based on SAP 10 notional building runs for 5 archetypes used
// for inter-model comparison/validation. The average of the monthly mean
// internal temperatures from each run was taken.
const HEAT_TRANSFER_OTHER_SIDE_GROUND_TEMP_INT_MONTHLY: [f64; 12] = [
    19.46399546,
    19.66940204,
    19.90785898,
    20.19719837,
    20.37461865,
    20.45679018,
    20.46767703,
    20.46860812,
    20.43505593,
    20.22266322,
    19.82726777,
    19.45430847,
];

pub(crate) trait HeatTransferOtherSideGround: HeatTransferOtherSide {
    fn init_heat_transfer_other_side_ground(
        &mut self,
        r_vi: f64,
        d_we: f64,
        thermal_conductivity: f64,
        r_si: f64,
        thermal_resistance_floor_construction: f64,
        floor_data: &FloorData,
        periodic_penetration_depth: f64,
        total_area: f64,
        perimeter: f64,
        u_value: f64,
        psi_wall_floor_junc: f64,
    ) -> anyhow::Result<()> {
        self.init_super(None);
        self.set_temp_int_annual(average_monthly_to_annual(
            HEAT_TRANSFER_OTHER_SIDE_GROUND_TEMP_INT_MONTHLY,
        ));

        // defining following closures out of order from the Python
        // as we are using closures which cannot be hoisted

        let total_equiv_thickness = || -> f64 {
            d_we + thermal_conductivity
                * (r_si + thermal_resistance_floor_construction + self.r_se())
        };

        let d_eq = total_equiv_thickness();

        // use sub-scope for calculating h_pi and h_pe so immutable references to self aren't
        // kept around for setters below.
        let (h_pi, h_pe) = {
            let init_unheated_basement = |u_f_s: f64, u_w: f64, z_b: f64, h_w: f64| -> (f64, f64) {
                let thermal_capacity_air = 0.33; // Wh/(m3·K)
                let air_vol_base = total_area * (h_w + z_b);

                // air changes per hour
                // From BS EN ISO 13370:2017 section 7.4
                let vent_rate_base = 0.3;

                // H.8.1. Internal temperature variation
                let h_pi: f64 = (1. / (total_area * u_f_s)
                    + 1. / ((total_area + z_b * perimeter) * thermal_conductivity
                        / periodic_penetration_depth
                        + h_w * perimeter * u_w
                        + thermal_capacity_air * vent_rate_base * air_vol_base))
                    .powi(-1);

                // H.8.2. External temperature variation
                // 0.37 is constant in the standard but not labelled
                let h_pe = total_area
                    * u_f_s
                    * (0.37
                        * perimeter
                        * thermal_conductivity
                        * (2. - (-z_b / periodic_penetration_depth).exp())
                        * (periodic_penetration_depth / d_eq + 1.).ln()
                        + h_w * perimeter * u_w
                        + thermal_capacity_air * vent_rate_base * air_vol_base)
                    / ((total_area + z_b * perimeter) * thermal_conductivity
                        / periodic_penetration_depth
                        + h_w * perimeter * u_w
                        + thermal_capacity_air * vent_rate_base * air_vol_base
                        + total_area * u_f_s);
                (h_pi, h_pe)
            };

            // Equivalent thickness for the basement walls
            let equiv_thick_base_wall = |r_w_b| -> f64 {
                // r_w_b is the thermal resistance of the walls
                thermal_conductivity * (r_si + r_w_b + self.r_se())
            };

            // Heated basement periodic coefficients
            let init_heated_basement = |z_b: f64, r_w_b: f64| -> (f64, f64) {
                // total equivalent thickness
                let d_w_b = equiv_thick_base_wall(r_w_b);

                // H.7.1. Internal temperature variation
                let h_pi = total_area
                    * ((thermal_conductivity / d_eq)
                        * (2. / (1. + periodic_penetration_depth / d_eq).powi(2) + 1.).powf(0.5))
                    + z_b
                        * perimeter
                        * (thermal_conductivity / d_w_b)
                        * (2. / ((1. + periodic_penetration_depth / d_w_b).powi(2) + 1.)).powf(0.5);

                // H.7.2. External temperature variation
                // 0.37 is constant in the standard but not labelled
                let h_pe = 0.37
                    * perimeter
                    * thermal_conductivity
                    * ((-z_b / periodic_penetration_depth).exp()
                        * (periodic_penetration_depth / d_eq + 1.).ln()
                        + 2. * (1. - (-z_b / periodic_penetration_depth).exp())
                            * (periodic_penetration_depth / d_w_b + 1.).ln());

                (h_pi, h_pe)
            };

            let charac_dimen_floor = || -> f64 {
                total_area / (0.5 * perimeter) // in m
            };

            // wind shielding factor
            let wind_shield_fact = |shield_fact_location| -> f64 {
                match shield_fact_location {
                    WindShieldLocation::Sheltered => 0.02,
                    WindShieldLocation::Average => 0.05,
                    WindShieldLocation::Exposed => 0.10,
                }
            };

            // equivalent thermal transmittance between the underfloor space and the outside
            let equiv_therma_trans = |h_upper,
                                      u_w,
                                      shield_fact_location: WindShieldLocation,
                                      area_per_perimeter_vent|
             -> anyhow::Result<f64> {
                // Characteristic dimension of floor
                let char_dimen = charac_dimen_floor();

                // 1450 is constant in the standard but not labelled
                Ok(2. * (h_upper * u_w / char_dimen)
                    + 1450.
                        * (area_per_perimeter_vent
                            * self.wind_speed()?
                            * wind_shield_fact(shield_fact_location))
                        / char_dimen)
            };

            let total_equiv_thickness_sus =
                |r_f_ins| -> f64 { d_we + thermal_conductivity * (r_si + r_f_ins + self.r_se()) };

            // thermal transmittance of suspended part of floor
            let thermal_transmittance_sus_floor =
                || -> f64 { 1. / (thermal_resistance_floor_construction + 2. * r_si) };

            // Suspended floor periodic coefficients
            let init_suspended_floor = |h_upper,
                                        u_w,
                                        area_vent,
                                        shield_fact_location,
                                        r_f_ins|
             -> anyhow::Result<(f64, f64)> {
                // H.6.1.
                // thermal transmittance of suspended part of floor, in W/(m2·K)
                let u_f = thermal_transmittance_sus_floor();

                // equivalent thermal transmittance, in W/(m2·K)
                let u_x = equiv_therma_trans(h_upper, u_w, shield_fact_location, area_vent)?;

                // equivalent thickness, in m
                let d_g = total_equiv_thickness_sus(r_f_ins);

                // H.6.2. Internal temperature variation
                let h_pi = total_area
                    * (1. / u_f + 1. / (thermal_conductivity / periodic_penetration_depth + u_x));

                // H.6.3. External temperature variation
                // 0.37 is constant in the standard but not labelled
                let h_pe = u_f
                    * ((0.37
                        * perimeter
                        * thermal_conductivity
                        * (periodic_penetration_depth / d_g + 1.).ln()
                        + u_x * total_area)
                        / (thermal_conductivity / (periodic_penetration_depth + u_x + u_f)));
                Ok((h_pi, h_pe))
            };

            // Additional equivalent thickness
            let add_eq_thickness = |d_n: f64, r_n: f64| -> f64 {
                // m2·K/W, thermal resistance
                let r_add_eq = r_n - d_n / thermal_conductivity;

                // m, thickness_edge-insulation or foundation
                r_add_eq * thermal_conductivity
            };

            // horizontal edge insulation
            let h_pe_h = |d_h: f64, r_n: f64| -> f64 {
                // 0.37 is constant in the standard but not labelled
                let eq_thick_additional = add_eq_thickness(d_h, r_n);

                0.37 * perimeter
                    * thermal_conductivity
                    * ((1. - (-d_h / periodic_penetration_depth).exp())
                        * (periodic_penetration_depth / (d_eq + eq_thick_additional) + 1.).ln()
                        + (-d_h / periodic_penetration_depth).exp()
                            * (periodic_penetration_depth / d_eq + 1.).ln())
            };

            // vertical edge insulation
            let h_pe_v = |d_v: f64, r_n: f64| -> f64 {
                // 0.37 is constant in the standard but not labelled
                let eq_thick_additional = add_eq_thickness(d_v, r_n);
                0.37 * perimeter
                    * thermal_conductivity
                    * ((1. - (-2. * d_v / periodic_penetration_depth).exp())
                        * (periodic_penetration_depth / (d_eq + eq_thick_additional) + 1.).ln()
                        + (-2. * d_v / periodic_penetration_depth).exp()
                            * (periodic_penetration_depth / d_eq + 1.).ln())
            };

            let edge_type = |edge_insulation: &[EdgeInsulation]| -> f64 {
                edge_insulation
                    .iter()
                    .map(|edge| match edge {
                        EdgeInsulation::Horizontal {
                            width,
                            edge_thermal_resistance,
                        } => h_pe_h(*width, *edge_thermal_resistance),
                        EdgeInsulation::Vertical {
                            depth,
                            edge_thermal_resistance,
                        } => h_pe_v(*depth, *edge_thermal_resistance),
                    })
                    .min_by(|a, b| a.total_cmp(b))
                    .expect("Edge insulation should not be empty")
            };

            let internal_temp_variation = || -> f64 {
                // H.4.1., H.5.1. Internal temperature variation
                total_area
                    * (thermal_conductivity / d_eq)
                    * (2. / ((1. + periodic_penetration_depth / d_eq).powi(2) + 1.0)).powf(0.5)
            };

            // Slab-on-ground floor uninsulated or with all-over insulated
            let init_slab_on_ground_floor_uninsulated_or_all_insulation = || -> (f64, f64) {
                // H.4.1. Internal temperature variation
                let h_pi = internal_temp_variation();

                // H.4.2. External temperature variation
                // 0.37 is constant in the standard but not labelled
                let h_pe = 0.37
                    * perimeter
                    * thermal_conductivity
                    * (periodic_penetration_depth / d_eq + 1.).ln();

                (h_pi, h_pe)
            };

            // Slab-on-ground-with-edge-insulation
            let init_slab_on_ground_floor_edge_insulated = |edge_insulation| -> (f64, f64) {
                // H.5.1. Internal temperature variation
                let h_pi = internal_temp_variation();

                // edge insulation (vertically or horizontally)
                let h_pe = edge_type(edge_insulation);

                (h_pi, h_pe)
            };

            // Return the periodic heat transfer coefficient for the building element
            //              h_pi     -- Internal periodic heat transfer coefficient, in W / K
            //                         BS EN ISO 13370:2017 Annex H
            //              h_pe     -- external periodic heat transfer coefficient, in W / K
            //                          BS EN ISO 13370:2017 Annex H
            let init_periodic_heat_transfer = || -> anyhow::Result<(f64, f64)> {
                Ok(match floor_data {
                    FloorData::SlabNoEdgeInsulation => {
                        init_slab_on_ground_floor_uninsulated_or_all_insulation()
                    }
                    FloorData::SlabEdgeInsulation { edge_insulation } => {
                        init_slab_on_ground_floor_edge_insulated(edge_insulation)
                    }
                    FloorData::SuspendedFloor {
                        thermal_transmission_walls,
                        area_per_perimeter_vent,
                        thermal_resistance_of_insulation,
                        height_upper_surface,
                        shield_fact_location,
                        ..
                    } => init_suspended_floor(
                        *height_upper_surface,
                        *thermal_transmission_walls,
                        *area_per_perimeter_vent,
                        *shield_fact_location,
                        *thermal_resistance_of_insulation,
                    )?,
                    FloorData::HeatedBasement {
                        depth_basement_floor,
                        thermal_resistance_of_basement_walls,
                        ..
                    } => init_heated_basement(
                        *depth_basement_floor,
                        *thermal_resistance_of_basement_walls,
                    ),
                    FloorData::UnheatedBasement {
                        thermal_transmittance_of_floor_above_basement,
                        thermal_transmission_walls,
                        depth_basement_floor,
                        height_basement_walls,
                        ..
                    } => init_unheated_basement(
                        *thermal_transmittance_of_floor_above_basement,
                        *thermal_transmission_walls,
                        *depth_basement_floor,
                        *height_basement_walls,
                    ),
                })
            };

            init_periodic_heat_transfer()?
        };

        self.set_total_area(total_area);
        self.set_perimeter(perimeter);
        self.set_u_value(u_value);
        self.set_psi_wall_floor_junc(psi_wall_floor_junc);
        self.set_floor_data(floor_data.clone());

        // Set external surface heat transfer coeffs as per BS EN ISO 52016-1:2017 eqn 49
        // Must be set before initialisation of base class, as these are referenced there
        // BS EN ISO 52016-1:2017 Table 14 validity interval h_ce 0 to 50
        self.set_h_ce(1.0 / r_vi); // in W/(m2.K)
        self.set_d_eq(d_eq);
        self.set_h_pi(h_pi);
        self.set_h_pe(h_pe);

        Ok(())
    }

    // setters to be implemented (when there is access to struct fields et al)
    fn set_temp_int_annual(&mut self, temp_int_annual: f64);
    fn temp_int_annual(&self) -> f64;

    fn wind_speed(&self) -> anyhow::Result<f64> {
        self.external_conditions().wind_speed_annual()
    }

    fn set_total_area(&mut self, total_area: f64);
    fn total_area(&self) -> f64;
    fn set_perimeter(&mut self, perimeter: f64);
    fn perimeter(&self) -> f64;
    fn set_u_value(&mut self, u_value: f64);
    fn set_psi_wall_floor_junc(&mut self, psi_wall_floor_junc: f64);
    fn psi_wall_floor_junc(&self) -> f64;
    fn set_floor_data(&mut self, floor_data: FloorData);
    fn set_h_ce(&mut self, h_ce: f64);
    fn set_d_eq(&mut self, d_eq: f64);
    fn set_h_pi(&mut self, h_pi: f64);
    fn h_pi(&self) -> f64;
    fn set_h_pe(&mut self, h_pe: f64);
    fn h_pe(&self) -> f64;

    fn u_value(&self) -> f64;

    /// Return external convective heat transfer coefficient, in W / (m2.K)
    fn h_ce(&self) -> f64;

    /// Return external radiative heat transfer coefficient, in W / (m2.K)
    fn h_re(&self) -> f64 {
        0.
    }

    /// Return the temperature on the other side of the building element
    fn temp_ext(&self, simtime: SimulationTimeIteration) -> f64 {
        let temp_ext_annual = self.external_conditions().air_temp_annual().expect("The external conditions provided do not contain data for an entire year, therefore an annual temperature cannot be calculated.");
        let temp_ext_month = self
            .external_conditions()
            .air_temp_monthly(simtime.current_month_start_end_hours());

        let current_month = simtime
            .current_month()
            .expect("A timestep should be able to derive its current month here.")
            as usize;
        let temp_int_month = HEAT_TRANSFER_OTHER_SIDE_GROUND_TEMP_INT_MONTHLY[current_month];

        // BS EN ISO 13370:2017 Eqn C.4
        let heat_flow_month =
            self.u_value() * self.total_area() * (self.temp_int_annual() - temp_ext_annual)
                + self.perimeter() * self.psi_wall_floor_junc() * (temp_int_month - temp_ext_month)
                - self.h_pi() * (self.temp_int_annual() - temp_int_month)
                + self.h_pe() * (temp_ext_annual - temp_ext_month);

        // BS EN ISO 13370:2017 Eqn F.2
        temp_int_month
            - (heat_flow_month
                - (self.perimeter()
                    * self.psi_wall_floor_junc()
                    * (self.temp_int_annual() - temp_ext_annual)))
                / (self.total_area() * self.u_value())
    }
}

pub(crate) trait HeatTransferOtherSideConditionedSpace: HeatTransferOtherSide {
    /// Return external convective heat transfer coefficient, in W / (m2.K)
    fn h_ce(&self) -> f64 {
        // Element is adjacent to another building / thermally conditioned zone
        // therefore according to BS EN ISO 52016-1:2017, section 6.5.6.3.6,
        // external heat transfer coefficients are zero
        0.0
    }

    /// Return external radiative heat transfer coefficient, in W / (m2.K)
    fn h_re(&self) -> f64 {
        // Element is adjacent to another building / thermally conditioned zone
        // therefore according to BS EN ISO 52016-1:2017, section 6.5.6.3.6,
        // external heat transfer coefficients are zero
        0.0
    }
}

pub(crate) trait HeatTransferOtherSideUnconditionedSpace: HeatTransferOtherSide {
    fn init_heat_transfer_other_side_unconditioned_space(
        &mut self,
        thermal_resistance_unconditioned_space: f64,
    ) {
        self.set_r_u(thermal_resistance_unconditioned_space);
        self.init_super(None);
    }

    fn set_r_u(&mut self, thermal_resistance_unconditioned_space: f64);
    fn r_u(&self) -> f64;

    /// Return external convective heat transfer coefficient, in W / (m2.K)
    fn h_ce(&self) -> f64 {
        // Add an additional thermal resistance to the outside of the wall and
        // incorporate this in the values for the external surface heat transfer
        // coefficient.
        // As this is an adjusted figure in this class, and the split between
        // h_ce and h_re does not affect the calculation results, assign entire
        // effective surface heat transfer to h_ce and set h_re to zero.
        let h_ce = self.h_ce_super();
        let h_re = self.h_re_super();
        let h_se = h_ce + h_re;
        let r_se = 1.0 / h_se;
        let r_se_effective = r_se + self.r_u();
        1.0 / r_se_effective
    }

    /// Return external radiative heat transfer coefficient, in W / (m2.K)
    fn h_re(&self) -> f64 {
        // As this is an adjusted figure in this class, and the split between
        // h_ce and h_re does not affect the calculation results, assign entire
        // effective surface heat transfer to h_ce and set h_re to zero.
        0.0
    }
}

pub(crate) trait HeatTransferOtherSideOutside: HeatTransferOtherSide {
    fn init_heat_transfer_other_side_outside(&mut self, pitch: f64) {
        let f_sky = ExternalConditions::sky_view_factor(&pitch);

        self.init_super(Some(f_sky));
    }
}

pub(crate) trait SolarRadiationInteraction {
    /// solar_absorption_coeff -- solar absorption coefficient at the external surface (dimensionless)
    /// according to BS EN ISO 52016-1:2017, section 6.5.7.3
    fn init_solar_radiation_interaction(
        &mut self,
        pitch: f64,
        orientation: Option<Orientation360>,
        shading: Option<Vec<WindowShadingObject>>,
        base_height: f64,
        projected_height: f64,
        width: f64,
        solar_absorption_coeff: f64,
    ) {
        self.set_external_pitch(pitch);
        if let Some(orientation) = orientation {
            self.set_orientation(orientation);
        }
        self.set_shading(shading);
        self.set_base_height(base_height);
        self.set_projected_height(projected_height);
        self.set_width(width);
        self.set_solar_absorption_coeff(solar_absorption_coeff);
    }

    fn set_external_pitch(&mut self, pitch: f64);
    fn external_pitch(&self) -> f64;
    fn set_orientation(&mut self, orientation: Orientation360);
    fn orientation(&self) -> Option<Orientation360>;
    fn set_shading(&mut self, shading: Option<Vec<WindowShadingObject>>);
    fn shading(&self) -> &[WindowShadingObject];
    fn set_base_height(&mut self, base_height: f64);
    fn base_height(&self) -> f64;
    fn set_projected_height(&mut self, projected_height: f64);
    fn projected_height(&self) -> f64;
    fn set_width(&mut self, width: f64);
    fn width(&self) -> f64;
    fn set_solar_absorption_coeff(&mut self, solar_absorption_coeff: f64);
    fn solar_absorption_coeff(&self) -> f64;

    fn i_sol_dir_dif(&self, _simtime: SimulationTimeIteration) -> (f64, f64) {
        // Return default of zero for i_sol_dir and i_sol_dif
        (0.0, 0.0)
    }

    fn solar_gains(&self) -> f64 {
        // Return default of zero for solar gains
        0.
    }

    fn shading_factors_direct_diffuse(&self, _simtime: SimulationTimeIteration) -> (f64, f64) {
        // Return default of one for shading factor (no shading)
        (1.0, 1.0)
    }
}

pub(crate) trait SolarRadiationInteractionAbsorbed: SolarRadiationInteraction {
    fn init_solar_radiation_interaction_absorbed(
        &mut self,
        pitch: f64,
        orientation: Option<Orientation360>,
        shading: Option<Vec<WindowShadingObject>>,
        base_height: f64,
        projected_height: f64,
        width: f64,
        solar_absorption_coeff: f64,
    ) {
        self.init_solar_radiation_interaction(
            pitch,
            orientation,
            shading,
            base_height,
            projected_height,
            width,
            solar_absorption_coeff,
        );
    }

    fn external_conditions(&self) -> &ExternalConditions;

    /// Return calculated i_sol_dir and i_sol_dif using pitch and orientation of element
    fn i_sol_dir_dif(&self, simtime: SimulationTimeIteration) -> anyhow::Result<(f64, f64)> {
        let CalculatedDirectDiffuseTotalIrradiance(i_sol_dir, i_sol_dif, _, _) = self
            .external_conditions()
            .calculated_direct_diffuse_total_irradiance(
                self.external_pitch(),
                self.orientation().ok_or_else(|| {
                    anyhow!(
                        "Cannot compute direct and diffuse irradiance because orientation is None."
                    )
                })?,
                false,
                &simtime,
            );
        Ok((i_sol_dir, i_sol_dif))
    }

    fn shading_factors_direct_diffuse(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        self.external_conditions()
            .shading_reduction_factor_direct_diffuse(
                self.base_height(),
                self.projected_height(),
                self.width(),
                self.external_pitch(),
                self.orientation().ok_or_else(|| {
                    anyhow!(
                        "Cannot compute direct and diffuse shading factors because orientation is None."
                    )
                })?,
                &[],
                simtime,
            )
    }
}

pub(crate) trait SolarRadiationInteractionTransmitted: SolarRadiationInteraction {
    // there is an init method in the Python for this that adds a simtime into the object there,
    // but as we do not inject simulation time references generally in the Rust, we can just inherit
    // the init method in SolarRadiationInteraction

    fn unconverted_g_value(&self) -> f64;

    /// return g_value corrected for angle of solar radiation
    fn convert_g_value(&self) -> f64 {
        // TODO (from Python) for windows with scattering glazing or solar shading provisions
        //      there is a different, more complex method for conversion that depends on
        //      timestep (via solar altitude).
        //      suggest this is implemented at the same time as window shading (devices
        //      rather than fixed features) as will also need to link to shading schedule.
        //      see ISO 52016 App E. Page 177
        //      How do we know whether a window has "scattering glazing"?
        //
        //      # g_value = agl * g_alt + (1 - agl) * g_dif

        let fw = 0.90;
        // default from ISO 52016 App B Table B.22
        fw * self.unconverted_g_value()
    }

    /// Return calculated solar gains using pitch and orientation of element
    fn solar_gains(&self, simtime: SimulationTimeIteration) -> anyhow::Result<f64> {
        let g_value = self.convert_g_value();
        let surf_irrad = self.external_conditions().surface_irradiance(
            self.base_height(),
            self.projected_height(),
            self.width(),
            self.pitch(),
            self.orientation().ok_or_else(|| {
                anyhow!("Cannot compute surface irradiance because orientation is None.")
            })?,
            self.shading(),
            simtime,
        )?;

        Ok(g_value * surf_irrad * self.area() * (1. - self.frame_area_fraction()))
    }

    fn external_conditions(&self) -> &ExternalConditions;
    fn pitch(&self) -> f64;
    fn area(&self) -> f64;
    fn frame_area_fraction(&self) -> f64;

    fn shading_factors_direct_diffuse(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        self.external_conditions()
            .shading_reduction_factor_direct_diffuse(
                self.base_height(),
                self.projected_height(),
                self.width(),
                self.pitch(),
                self.orientation().ok_or_else(|| {
                    anyhow!("Cannot compute shading factor because orientation is None.")
                })?,
                self.shading(),
                simtime,
            )
    }
}

#[allow(dead_code)]
pub(crate) trait SolarRadiationInteractionNotExposed: SolarRadiationInteraction {}

/// A type to represent opaque building elements (walls, roofs, etc.)
/// TODO make this into canonical BuildingElementOpaque
#[derive(Clone, Debug)]
pub(crate) struct BuildingElementOpaque {
    area: f64,
    external_conditions: Arc<ExternalConditions>,
    shading: Option<Vec<WindowShadingObject>>,
    r_c: f64,
    _internal_pitch: f64,
    external_pitch: f64,
    solar_absorption_coeff: f64,
    base_height: f64,
    projected_height: f64,
    width: f64,
    orientation: Option<Orientation360>,
    k_m: f64,
    k_pli: [f64; 5],
    h_pli: [f64; 4],
    f_sky: f64,
    therm_rad_to_sky: f64,
}

impl BuildingElementOpaque {
    /// Arguments (names based on those in BS EN ISO 52016-1:2017):
    /// * `area` - net area of the opaque building element (i.e. minus any windows / doors / etc.)
    /// * `is_unheated_pitched_roof`
    /// * `pitch` - tilt angle of the surface from horizontal, in degrees between 0 and 180,
    ///          where 0 means the external surface is facing up, 90 means the external
    ///          surface is vertical and 180 means the external surface is facing down
    /// * `solar_absorption_coeff`    - solar absorption coefficient at the external surface (dimensionless)
    /// * `thermal_resistance_construction`      - thermal resistance, in m2.K / W
    /// * `areal_heat_capacity`      - areal heat capacity, in J / (m2.K)
    /// * `mass_distribution_class`
    ///          - distribution of mass in building element, one of:
    ///             - 'I':  mass concentrated on internal side
    ///             - 'E':  mass concentrated on external side
    ///             - 'IE': mass divided over internal and external side
    ///             - 'D':  mass equally distributed
    ///             - 'M':  mass concentrated inside
    /// * `orientation` -- is the orientation angle of the inclined surface, expressed as the
    ///                geographical azimuth angle of the horizontal projection of the inclined
    ///                surface normal, 0 to 360, in degrees
    /// * `base_height` - is the distance between the ground and the lowest edge of the element, in m
    /// * `height`      - is the height of the building element, in m
    /// * `width`       - is the width of the building element, in m
    /// * `external_conditions` -- reference to ExternalConditions object
    pub(crate) fn new(
        area: f64,
        is_unheated_pitched_roof: bool,
        pitch: f64,
        solar_absorption_coeff: f64,
        thermal_resistance_construction: f64,
        areal_heat_capacity: f64,
        mass_distribution_class: MassDistributionClass,
        orientation: Option<Orientation360>,
        base_height: f64,
        height: f64,
        width: f64,
        external_conditions: Arc<ExternalConditions>,
    ) -> Self {
        let mut new_opaque: Self = Self {
            area,
            external_conditions,
            shading: None,
            r_c: Default::default(),
            _internal_pitch: if is_unheated_pitched_roof { 0. } else { pitch },
            external_pitch: pitch,
            solar_absorption_coeff,
            base_height,
            projected_height: Default::default(),
            width,
            orientation,
            k_m: Default::default(),
            k_pli: Default::default(),
            h_pli: Default::default(),
            f_sky: Default::default(),
            therm_rad_to_sky: Default::default(),
        };

        new_opaque.init_heat_transfer_through_5_nodes(
            thermal_resistance_construction,
            mass_distribution_class,
            areal_heat_capacity,
        );
        new_opaque.init_heat_transfer_other_side_outside(pitch);
        // shading is None because the model ignores nearby shading on opaque elements
        new_opaque.init_solar_radiation_interaction_absorbed(
            pitch,
            orientation,
            None,
            base_height,
            projected_height(pitch, height),
            width,
            solar_absorption_coeff,
        );

        new_opaque
    }

    pub(crate) fn solar_absorption_coeff(&self) -> f64 {
        self.solar_absorption_coeff
    }

    pub(crate) fn i_sol_dir_dif(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        SolarRadiationInteractionAbsorbed::i_sol_dir_dif(self, simtime)
    }

    pub(crate) fn shading_factors_direct_diffuse(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        SolarRadiationInteractionAbsorbed::shading_factors_direct_diffuse(self, simtime)
    }

    #[cfg(test)]
    pub(crate) fn r_si(&self) -> f64 {
        HeatTransferInternal::r_si(self)
    }
}

impl HeatTransferInternal for BuildingElementOpaque {
    fn pitch(&self) -> f64 {
        self.external_pitch()
    }
}
impl HeatTransferInternalCommon for BuildingElementOpaque {}

impl HeatTransferThrough for BuildingElementOpaque {
    fn set_r_c(&mut self, thermal_resistance_construction: f64) {
        self.r_c = thermal_resistance_construction;
    }

    fn r_c(&self) -> f64 {
        self.r_c
    }

    fn k_m(&self) -> f64 {
        self.k_m
    }

    fn set_k_m(&mut self, areal_heat_capacity: f64) {
        self.k_m = areal_heat_capacity;
    }

    fn k_pli(&self) -> &[f64] {
        &self.k_pli
    }

    fn r_se(&self) -> f64 {
        // upstream Python uses duck typing/ lookups to find which method this is, but Rust needs to be explicit
        <dyn HeatTransferOtherSide>::r_se(self)
    }

    fn r_si(&self) -> f64 {
        <dyn HeatTransferInternal>::r_si(self)
    }

    fn area(&self) -> f64 {
        self.area
    }

    fn h_pli(&self) -> &[f64] {
        &self.h_pli
    }
}

impl HeatTransferThrough5Nodes for BuildingElementOpaque {
    fn set_h_pli(&mut self, h_pli: [f64; 4]) {
        self.h_pli = h_pli;
    }

    fn set_k_pli(&mut self, k_pli: [f64; 5]) {
        self.k_pli = k_pli;
    }
}

impl HeatTransferOtherSide for BuildingElementOpaque {
    fn set_f_sky(&mut self, f_sky: f64) {
        self.f_sky = f_sky;
    }

    fn f_sky(&self) -> f64 {
        self.f_sky
    }

    fn set_therm_rad_to_sky(&mut self, therm_rad_to_sky: f64) {
        self.therm_rad_to_sky = therm_rad_to_sky;
    }

    fn therm_rad_to_sky(&self) -> f64 {
        self.therm_rad_to_sky
    }

    fn external_conditions(&self) -> &ExternalConditions {
        self.external_conditions.as_ref()
    }
}

impl HeatTransferOtherSideOutside for BuildingElementOpaque {}

impl SolarRadiationInteraction for BuildingElementOpaque {
    fn set_external_pitch(&mut self, pitch: f64) {
        self.external_pitch = pitch;
    }

    fn external_pitch(&self) -> f64 {
        self.external_pitch
    }

    fn set_orientation(&mut self, orientation: Orientation360) {
        self.orientation = orientation.into();
    }

    fn orientation(&self) -> Option<Orientation360> {
        self.orientation
    }

    fn set_shading(&mut self, shading: Option<Vec<WindowShadingObject>>) {
        self.shading = shading;
    }

    fn shading(&self) -> &[WindowShadingObject] {
        match self.shading {
            Some(ref shading) => &shading[..],
            None => &[],
        }
    }

    fn set_base_height(&mut self, base_height: f64) {
        self.base_height = base_height;
    }

    fn base_height(&self) -> f64 {
        self.base_height
    }

    fn set_projected_height(&mut self, projected_height: f64) {
        self.projected_height = projected_height;
    }

    fn projected_height(&self) -> f64 {
        self.projected_height
    }

    fn set_width(&mut self, width: f64) {
        self.width = width;
    }

    fn width(&self) -> f64 {
        self.width
    }

    fn set_solar_absorption_coeff(&mut self, solar_absorption_coeff: f64) {
        self.solar_absorption_coeff = solar_absorption_coeff;
    }

    fn solar_absorption_coeff(&self) -> f64 {
        self.solar_absorption_coeff
    }
}

impl SolarRadiationInteractionAbsorbed for BuildingElementOpaque {
    fn external_conditions(&self) -> &ExternalConditions {
        self.external_conditions.as_ref()
    }
}

/// A class to represent all party wall building elements
///
/// This class handles all party wall types, with specific provisions for cavity walls
/// where heat loss via air movement occurs. Research has shown this heat loss to be
/// significant in cavity party walls, contrary to historical assumptions of negligible
/// heat transfer.
///
/// For cavity party walls, heat loss is modelled by treating the cavity as an
/// unconditioned space with an effective thermal resistance that depends on the
/// cavity type. The thermal resistance values (R_cavity) are derived from equivalent
/// U-values which are based on research into air movement patterns within party wall cavities:
/// - Solid walls and filled-sealed cavities: R_cavity = 0 (no heat loss via cavity)
/// - Unfilled unsealed cavities: R_cavity ≈ 1.5 m².K/W
/// - Unfilled sealed cavities: R_cavity ≈ 4.5 m².K/W
///
/// Note: This approach does not explicitly model losses via air movement in party
/// ceilings/floors, though these may be implicitly included in the cavity resistance
/// values derived from the original research.
#[derive(Clone, Debug)]
pub(crate) struct BuildingElementPartyWall {
    area: f64,
    pitch: f64,
    thermal_resistance_construction: f64,
    party_wall_cavity_type: PartyWallCavityType,
    party_wall_lining_type: Option<PartyWallLiningType>,
    thermal_resistance_cavity: Option<f64>,
    areal_heat_capacity: f64,
    mass_distribution_class: MassDistributionClass,
    external_conditions: Arc<ExternalConditions>,
    h_pli: [f64; 4],
    k_pli: [f64; 5],
    r_c: f64,
    k_m: f64,
    f_sky: f64,
    therm_rad_to_sky: f64,
    r_u: f64,
    external_pitch: f64,
    orientation: Option<Orientation360>,
    shading: Option<Vec<WindowShadingObject>>,
    base_height: f64,
    projected_height: f64,
    width: f64,
    solar_absorption_coeff: f64,
}

impl BuildingElementPartyWall {
    /// Construct a BuildingElementPartyWall object
    /// Arguments (names based on those in BS EN ISO 52016-1:2017):
    /// * `area` - area (in m2) of this building element
    /// * `pitch` - tilt angle of the surface from horizontal, in degrees between 0 and 180,
    ///             where` 0 means the external surface is facing up, 90 means the external
    ///             surface` is vertical and 180 means the external surface is facing down
    /// * `thermal_resistance_construction` - thermal resistance of wall layers before the cavity, in m2.K / W
    /// * `party_wall_cavity_type` - type of party wall cavity (solid, unfilled_unsealed,
    ///                              unfilled_sealed, filled_sealed, or defined_resistance)
    /// * `party_wall_lining_type` - type of party wall lining (wet_plaster, or dry_lined)
    /// * `thermal_resistance_cavity` - effective thermal resistance of the cavity, in m2.K / W;
    ///                                 required only for 'defined_resistance' type, otherwise calculated automatically
    /// * `areal_heat_capacity` - areal heat capacity, in J / (m2.K)
    /// * `mass_distribution_class` - distribution of mass in building element, one of:
    ///     - 'I':  mass concentrated on internal side
    ///     - 'E':  mass concentrated on external side
    ///     - 'IE': mass divided over internal and external side
    ///     - 'D':  mass equally distributed
    ///     - 'M':  mass concentrated inside
    /// * `ext_cond` -- reference to ExternalConditions object
    ///
    /// Other variables:
    /// * `f_sky` - view factor to the sky (see BS EN ISO 52016-1:2017, section 6.5.13.3)
    /// * `h_ce` - external convective heat transfer coefficient, in W / (m2.K)
    /// * `h_re` - external radiative heat transfer coefficient, in W / (m2.K)
    /// * `solar_absorption_coeff` - solar absorption coefficient at the external surface (dimensionless)
    fn new(
        area: f64,
        pitch: f64,
        thermal_resistance_construction: f64,
        party_wall_cavity_type: PartyWallCavityType,
        party_wall_lining_type: Option<PartyWallLiningType>,
        thermal_resistance_cavity: Option<f64>,
        areal_heat_capacity: f64,
        mass_distribution_class: MassDistributionClass,
        external_conditions: Arc<ExternalConditions>,
    ) -> anyhow::Result<Self> {
        // Calculate the effective thermal resistance of the unconditioned space (cavity)
        // based on the party wall cavity type and party wall lining type

        let r_unconditioned = calculate_cavity_resistance(
            &party_wall_cavity_type,
            &party_wall_lining_type,
            thermal_resistance_cavity,
        )?;

        let mut party_wall = Self {
            area,
            pitch,
            thermal_resistance_construction,
            party_wall_cavity_type,
            party_wall_lining_type,
            thermal_resistance_cavity,
            areal_heat_capacity,
            mass_distribution_class,
            external_conditions,
            h_pli: Default::default(),
            k_pli: Default::default(),
            r_c: Default::default(),
            k_m: Default::default(),
            f_sky: Default::default(),
            therm_rad_to_sky: Default::default(),
            r_u: r_unconditioned,
            external_pitch: Default::default(),
            orientation: Default::default(),
            shading: None,
            base_height: Default::default(),
            projected_height: Default::default(),
            width: Default::default(),
            solar_absorption_coeff: Default::default(),
        };

        party_wall.init_heat_transfer_through_5_nodes(
            thermal_resistance_construction,
            mass_distribution_class,
            areal_heat_capacity,
        );
        party_wall.init_heat_transfer_other_side_unconditioned_space(r_unconditioned);

        // Solar absorption coefficient at the external surface is zero
        // (party walls are not exposed to solar radiation)
        party_wall.init_solar_radiation_interaction(pitch, None, None, 0., 0., 0., 0.);

        Ok(party_wall)
    }

    /// Return external convective heat transfer coefficient, in W / (m2.K)
    pub(crate) fn h_ce(&self) -> f64 {
        // For party walls with solid or filled_sealed cavity types, return zero
        // to represent adiabatic boundary condition (no heat loss)

        match self.party_wall_cavity_type {
            PartyWallCavityType::Solid | PartyWallCavityType::FilledSealed => 0.,
            _ => HeatTransferOtherSideUnconditionedSpace::h_ce(self),
        }
    }

    /// Return the fabric heat loss for the building element
    ///
    /// For party walls with solid or filled_sealed cavity types, return zero
    /// as these represent adiabatic boundaries with no heat loss to account for
    /// in the building's heat loss calculation.
    ///
    /// For cavity party walls, the fabric heat loss calculation must include
    /// the effective thermal resistance of the cavity (R_cavity) in addition to
    /// the construction resistance.
    fn fabric_heat_loss(&self) -> f64 {
        // For solid and filled_sealed types, there is no heat loss
        match self.party_wall_cavity_type {
            PartyWallCavityType::Solid | PartyWallCavityType::FilledSealed => 0.,
            _ => {
                // For cavity party walls, calculate U-value including cavity resistance
                // The parent class h_ce() already incorporates the cavity resistance,
                // so we need to calculate the effective R_se from h_ce and h_re
                let h_ce = self.h_ce();
                let h_re = self.h_re();
                let h_se_effective = h_ce + h_re;

                let r_se_effective = 1.0 / h_se_effective;

                // Calculate U-value with effective external surface resistance
                // skipping check on r_c as always set from thermal_resistance_construction as part of init_heat_transfer_through_5_nodes
                let u_value = 1.0 / (self.r_c + r_se_effective + HeatTransferInternal::r_si(self));
                self.area * u_value
            }
        }
    }

    pub(crate) fn h_re(&self) -> f64 {
        HeatTransferOtherSideUnconditionedSpace::h_re(self)
    }

    #[cfg(test)]
    pub(crate) fn r_si(&self) -> f64 {
        HeatTransferInternal::r_si(self)
    }
}

impl HeatTransferThrough for BuildingElementPartyWall {
    fn set_r_c(&mut self, thermal_resistance_construction: f64) {
        self.r_c = thermal_resistance_construction;
    }

    fn r_c(&self) -> f64 {
        self.r_c
    }

    fn k_m(&self) -> f64 {
        self.k_m
    }

    fn set_k_m(&mut self, areal_heat_capacity: f64) {
        self.k_m = areal_heat_capacity;
    }

    fn k_pli(&self) -> &[f64] {
        &self.k_pli
    }

    fn r_se(&self) -> f64 {
        // upstream Python uses duck typing/ lookups to find which method this is, but Rust needs to be explicit
        <dyn HeatTransferOtherSide>::r_se(self)
    }

    fn r_si(&self) -> f64 {
        <dyn HeatTransferInternal>::r_si(self)
    }

    fn area(&self) -> f64 {
        self.area
    }

    fn h_pli(&self) -> &[f64] {
        &self.h_pli
    }
}

impl HeatTransferOtherSide for BuildingElementPartyWall {
    fn set_f_sky(&mut self, f_sky: f64) {
        self.f_sky = f_sky;
    }

    fn set_therm_rad_to_sky(&mut self, therm_rad_to_sky: f64) {
        self.therm_rad_to_sky = therm_rad_to_sky;
    }

    fn therm_rad_to_sky(&self) -> f64 {
        self.therm_rad_to_sky
    }

    fn f_sky(&self) -> f64 {
        self.f_sky
    }

    fn external_conditions(&self) -> &ExternalConditions {
        self.external_conditions.as_ref()
    }
}

impl HeatTransferOtherSideUnconditionedSpace for BuildingElementPartyWall {
    fn set_r_u(&mut self, thermal_resistance_unconditioned_space: f64) {
        self.r_u = thermal_resistance_unconditioned_space;
    }

    fn r_u(&self) -> f64 {
        self.r_u
    }
}

impl HeatTransferThrough5Nodes for BuildingElementPartyWall {
    fn set_h_pli(&mut self, h_pli: [f64; 4]) {
        self.h_pli = h_pli;
    }

    fn set_k_pli(&mut self, k_pli: [f64; 5]) {
        self.k_pli = k_pli;
    }
}

impl HeatTransferInternal for BuildingElementPartyWall {
    fn pitch(&self) -> f64 {
        self.external_pitch()
    }
}
impl HeatTransferInternalCommon for BuildingElementPartyWall {}

impl SolarRadiationInteraction for BuildingElementPartyWall {
    fn set_external_pitch(&mut self, pitch: f64) {
        self.external_pitch = pitch;
    }

    fn external_pitch(&self) -> f64 {
        self.external_pitch
    }

    fn set_orientation(&mut self, orientation: Orientation360) {
        self.orientation = Some(orientation);
    }

    fn orientation(&self) -> Option<Orientation360> {
        self.orientation
    }

    fn set_shading(&mut self, shading: Option<Vec<WindowShadingObject>>) {
        self.shading = shading;
    }

    fn shading(&self) -> &[WindowShadingObject] {
        match self.shading {
            Some(ref shading) => &shading[..],
            None => &[],
        }
    }

    fn set_base_height(&mut self, base_height: f64) {
        self.base_height = base_height;
    }

    fn base_height(&self) -> f64 {
        self.base_height
    }

    fn set_projected_height(&mut self, projected_height: f64) {
        self.projected_height = projected_height;
    }

    fn projected_height(&self) -> f64 {
        self.projected_height
    }

    fn set_width(&mut self, width: f64) {
        self.width = width;
    }

    fn width(&self) -> f64 {
        self.width
    }

    fn set_solar_absorption_coeff(&mut self, solar_absorption_coeff: f64) {
        self.solar_absorption_coeff = solar_absorption_coeff;
    }

    fn solar_absorption_coeff(&self) -> f64 {
        self.solar_absorption_coeff
    }
}

impl SolarRadiationInteractionNotExposed for BuildingElementPartyWall {}

#[derive(Debug)]
pub(crate) struct BuildingElementAdjacentConditionedSpace {
    area: f64,
    pitch: f64,
    external_conditions: Arc<ExternalConditions>,
    r_c: f64,
    h_pli: [f64; 4],
    k_pli: [f64; 5],
    k_m: f64,
    f_sky: f64,
    therm_rad_to_sky: f64,
    external_pitch: f64,
    orientation: Orientation360,
}

/// A type to represent building elements adjacent to a thermally conditioned zone (ZTC)
impl BuildingElementAdjacentConditionedSpace {
    /// Arguments (names based on those in BS EN ISO 52016-1:2017):
    /// * `area`     - area (in m2) of this building element
    /// * `pitch` - tilt angle of the surface from horizontal, in degrees between 0 and 180,
    ///          where 0 means the external surface is facing up, 90 means the external
    ///          surface is vertical and 180 means the external surface is facing down
    /// * `thermal_resistance_construction` - thermal resistance, in m2.K / W
    /// * `areal_heat_capacity` - areal heat capacity, in J / (m2.K)
    /// * `mass_distribution_class` - distribution of mass in building element, one of:
    ///    - 'I':  mass concentrated on internal side
    ///    - 'E':  mass concentrated on external side
    ///    - 'IE': mass divided over internal and external side
    ///    - 'D':  mass equally distributed
    ///    - 'M':  mass concentrated inside
    /// * `external_conditions` - reference to ExternalConditions object
    ///
    /// Other variables:
    /// * `f_sky` - view factor to the sky (see BS EN ISO 52016-1:2017, section 6.5.13.3)
    /// * `h_ce` - external convective heat transfer coefficient, in W / (m2.K)
    /// * `h_re` - external radiative heat transfer coefficient, in W / (m2.K)
    /// * `solar_absorption_coeff` - solar absorption coefficient at the external surface (dimensionless)
    pub(crate) fn new(
        area: f64,
        pitch: f64,
        thermal_resistance_construction: f64,
        areal_heat_capacity: f64,
        mass_distribution_class: MassDistributionClass,
        external_conditions: Arc<ExternalConditions>,
    ) -> Self {
        let mut building_element = Self {
            area,
            pitch,
            external_conditions,
            r_c: Default::default(),
            k_m: Default::default(),
            h_pli: Default::default(),
            k_pli: Default::default(),
            f_sky: Default::default(),
            therm_rad_to_sky: Default::default(),
            external_pitch: Default::default(),
            orientation: Default::default(),
        };

        building_element.init_heat_transfer_through_5_nodes(
            thermal_resistance_construction,
            mass_distribution_class,
            areal_heat_capacity,
        );
        building_element.init_heat_transfer_other_side(None);
        building_element.init_solar_radiation_interaction(pitch, None, None, 0.0, 0.0, 0.0, 0.0);

        building_element
    }

    pub(crate) fn fabric_heat_loss(&self) -> f64 {
        0.0 // no heat loss to thermally conditioned zones
    }

    pub(crate) fn h_ce(&self) -> f64 {
        HeatTransferOtherSideConditionedSpace::h_ce(self)
    }

    pub(crate) fn h_re(&self) -> f64 {
        HeatTransferOtherSideConditionedSpace::h_re(self)
    }

    #[cfg(test)]
    pub(crate) fn r_si(&self) -> f64 {
        HeatTransferInternal::r_si(self)
    }
}

impl HeatTransferThrough for BuildingElementAdjacentConditionedSpace {
    fn set_r_c(&mut self, thermal_resistance_construction: f64) {
        self.r_c = thermal_resistance_construction;
    }

    fn r_c(&self) -> f64 {
        self.r_c
    }

    fn k_m(&self) -> f64 {
        self.k_m
    }

    fn set_k_m(&mut self, areal_heat_capacity: f64) {
        self.k_m = areal_heat_capacity;
    }

    fn k_pli(&self) -> &[f64] {
        &self.k_pli
    }

    fn r_se(&self) -> f64 {
        // upstream Python uses duck typing/ lookups to find which method this is, but Rust needs to be explicit
        <dyn HeatTransferOtherSide>::r_se(self)
    }

    fn r_si(&self) -> f64 {
        <dyn HeatTransferInternal>::r_si(self)
    }

    fn area(&self) -> f64 {
        self.area
    }

    fn h_pli(&self) -> &[f64] {
        &self.h_pli
    }
}

impl HeatTransferThrough5Nodes for BuildingElementAdjacentConditionedSpace {
    fn set_h_pli(&mut self, h_pli: [f64; 4]) {
        self.h_pli = h_pli;
    }

    fn set_k_pli(&mut self, k_pli: [f64; 5]) {
        self.k_pli = k_pli;
    }
}

impl HeatTransferInternal for BuildingElementAdjacentConditionedSpace {
    fn pitch(&self) -> f64 {
        self.pitch
    }
}

impl HeatTransferOtherSide for BuildingElementAdjacentConditionedSpace {
    fn set_f_sky(&mut self, f_sky: f64) {
        self.f_sky = f_sky;
    }

    fn set_therm_rad_to_sky(&mut self, therm_rad_to_sky: f64) {
        self.therm_rad_to_sky = therm_rad_to_sky;
    }

    fn therm_rad_to_sky(&self) -> f64 {
        self.therm_rad_to_sky
    }

    fn f_sky(&self) -> f64 {
        self.f_sky
    }

    fn external_conditions(&self) -> &ExternalConditions {
        self.external_conditions.as_ref()
    }
}

impl HeatTransferOtherSideConditionedSpace for BuildingElementAdjacentConditionedSpace {}

impl HeatTransferInternalCommon for BuildingElementAdjacentConditionedSpace {}

impl SolarRadiationInteraction for BuildingElementAdjacentConditionedSpace {
    fn set_external_pitch(&mut self, pitch: f64) {
        self.external_pitch = pitch;
    }

    fn external_pitch(&self) -> f64 {
        self.external_pitch
    }

    fn set_orientation(&mut self, orientation: Orientation360) {
        self.orientation = orientation;
    }

    fn orientation(&self) -> Option<Orientation360> {
        self.orientation.into()
    }

    fn set_shading(&mut self, _shading: Option<Vec<WindowShadingObject>>) {
        // do nothing
    }

    fn shading(&self) -> &[WindowShadingObject] {
        &[]
    }

    fn set_base_height(&mut self, _base_height: f64) {
        // do nothing
    }

    fn base_height(&self) -> f64 {
        0.0
    }

    fn set_projected_height(&mut self, _projected_height: f64) {
        // do nothing
    }

    fn projected_height(&self) -> f64 {
        0.0
    }

    fn set_width(&mut self, _width: f64) {
        // do nothing
    }

    fn width(&self) -> f64 {
        0.0
    }

    fn set_solar_absorption_coeff(&mut self, _solar_absorption_coeff: f64) {
        // do nothing
    }

    fn solar_absorption_coeff(&self) -> f64 {
        0.0
    }
}

impl SolarRadiationInteractionNotExposed for BuildingElementAdjacentConditionedSpace {}

/// A type to represent building elements adjacent to a thermally unconditioned zone (ZTU)
///
/// This class uses a simple calculation by adding an additional thermal
/// resistance to the outside of the wall and incorporating this in the values
/// for the external surface heat transfer coefficients. This differs from both
/// of the approaches (internal and external) in BS EN ISO 52016-1:2017 which
/// require detailed inputs for the unconditioned zone.
#[derive(Debug)]
pub(crate) struct BuildingElementAdjacentUnconditionedSpaceSimple {
    area: f64,
    pitch: f64,
    external_pitch: f64,
    external_conditions: Arc<ExternalConditions>,
    r_u: f64,
    f_sky: f64,
    therm_rad_to_sky: f64,
    k_m: f64,
    h_pli: [f64; 4],
    k_pli: [f64; 5],
    r_c: f64,
}

impl BuildingElementAdjacentUnconditionedSpaceSimple {
    /// Arguments (names based on those in BS EN ISO 52016-1:2017):
    /// * `area`     - area (in m2) of this building element
    /// * `pitch` - tilt angle of the surface from horizontal, in degrees between 0 and 180,
    ///          where 0 means the external surface is facing up, 90 means the external
    ///          surface is vertical and 180 means the external surface is facing down
    /// * `thermal_resistance_construction`      - thermal resistance, in m2.K / W
    /// * `thermal_resistance_unconditioned_space`      - effective thermal resistance of unheated space, in m2.K / W;
    ///             see SAP 10.2 section 3.3 for suggested values
    /// * `areal_heat_capacity`      - areal heat capacity, in J / (m2.K)
    /// * `mass_distribution_class`
    ///          - distribution of mass in building element, one of:
    ///             - 'I':  mass concentrated on internal side
    ///             - 'E':  mass concentrated on external side
    ///             - 'IE': mass divided over internal and external side
    ///             - 'D':  mass equally distributed
    ///             - 'M':  mass concentrated inside
    /// * `external_conditions` - reference to ExternalConditions object
    ///      Other variables:
    /// * `f_sky` -- view factor to the sky (see BS EN ISO 52016-1:2017, section 6.5.13.3)
    /// * `h_ce`     -- external convective heat transfer coefficient, in W / (m2.K)
    /// * `h_re`     -- external radiative heat transfer coefficient, in W / (m2.K)
    /// * `solar_absorption_coeff`    -- solar absorption coefficient at the external surface (dimensionless)
    pub(crate) fn new(
        area: f64,
        pitch: f64,
        thermal_resistance_construction: f64,
        thermal_resistance_unconditioned_space: f64,
        areal_heat_capacity: f64,
        mass_distribution_class: MassDistributionClass,
        external_conditions: Arc<ExternalConditions>,
    ) -> Self {
        let mut building_element = Self {
            area,
            pitch,
            external_pitch: 0.0,
            external_conditions,
            r_u: Default::default(),
            f_sky: Default::default(),
            therm_rad_to_sky: Default::default(),
            k_m: Default::default(),
            h_pli: Default::default(),
            k_pli: Default::default(),
            r_c: thermal_resistance_construction,
        };
        building_element.init_heat_transfer_through_5_nodes(
            thermal_resistance_construction,
            mass_distribution_class,
            areal_heat_capacity,
        );
        building_element.init_heat_transfer_other_side_unconditioned_space(
            thermal_resistance_unconditioned_space,
        );
        //  Solar absorption coefficient at the external surface is zero
        building_element.init_solar_radiation_interaction(pitch, None, None, 0.0, 0.0, 0.0, 0.0);

        building_element
    }

    pub(crate) fn h_ce(&self) -> f64 {
        HeatTransferOtherSideUnconditionedSpace::h_ce(self)
    }

    pub(crate) fn h_re(&self) -> f64 {
        HeatTransferOtherSideUnconditionedSpace::h_re(self)
    }
}

impl HeatTransferInternal for BuildingElementAdjacentUnconditionedSpaceSimple {
    fn pitch(&self) -> f64 {
        self.pitch
    }
}

impl HeatTransferInternalCommon for BuildingElementAdjacentUnconditionedSpaceSimple {}

impl HeatTransferThrough for BuildingElementAdjacentUnconditionedSpaceSimple {
    fn set_r_c(&mut self, thermal_resistance_construction: f64) {
        self.r_c = thermal_resistance_construction;
    }

    fn r_c(&self) -> f64 {
        self.r_c
    }

    fn k_m(&self) -> f64 {
        self.k_m
    }

    fn set_k_m(&mut self, areal_heat_capacity: f64) {
        self.k_m = areal_heat_capacity;
    }

    fn k_pli(&self) -> &[f64] {
        &self.k_pli
    }

    fn r_se(&self) -> f64 {
        // upstream Python uses duck typing/ lookups to find which method this is, but Rust needs to be explicit
        <dyn HeatTransferOtherSide>::r_se(self)
    }

    fn r_si(&self) -> f64 {
        <dyn HeatTransferInternal>::r_si(self)
    }

    fn area(&self) -> f64 {
        self.area
    }

    fn h_pli(&self) -> &[f64] {
        &self.h_pli
    }
}

impl HeatTransferThrough5Nodes for BuildingElementAdjacentUnconditionedSpaceSimple {
    fn set_h_pli(&mut self, h_pli: [f64; 4]) {
        self.h_pli = h_pli;
    }

    fn set_k_pli(&mut self, k_pli: [f64; 5]) {
        self.k_pli = k_pli;
    }
}

impl HeatTransferOtherSide for BuildingElementAdjacentUnconditionedSpaceSimple {
    fn set_f_sky(&mut self, f_sky: f64) {
        self.f_sky = f_sky;
    }

    fn set_therm_rad_to_sky(&mut self, therm_rad_to_sky: f64) {
        self.therm_rad_to_sky = therm_rad_to_sky;
    }

    fn therm_rad_to_sky(&self) -> f64 {
        self.therm_rad_to_sky
    }

    fn f_sky(&self) -> f64 {
        self.f_sky
    }

    fn external_conditions(&self) -> &ExternalConditions {
        self.external_conditions.as_ref()
    }
}

impl HeatTransferOtherSideUnconditionedSpace for BuildingElementAdjacentUnconditionedSpaceSimple {
    fn set_r_u(&mut self, thermal_resistance_unconditioned_space: f64) {
        self.r_u = thermal_resistance_unconditioned_space;
    }

    fn r_u(&self) -> f64 {
        self.r_u
    }
}

impl SolarRadiationInteraction for BuildingElementAdjacentUnconditionedSpaceSimple {
    fn set_external_pitch(&mut self, pitch: f64) {
        self.external_pitch = pitch;
    }

    fn external_pitch(&self) -> f64 {
        self.external_pitch
    }

    fn set_orientation(&mut self, _orientation: Orientation360) {
        // do nothing
    }

    fn orientation(&self) -> Option<Orientation360> {
        Some(0.0.into())
    }

    fn set_shading(&mut self, _shading: Option<Vec<WindowShadingObject>>) {
        // do nothing
    }

    fn shading(&self) -> &[WindowShadingObject] {
        &[]
    }

    fn set_base_height(&mut self, _base_height: f64) {
        // do nothing
    }

    fn base_height(&self) -> f64 {
        0.0
    }

    fn set_projected_height(&mut self, _projected_height: f64) {
        // do nothing
    }

    fn projected_height(&self) -> f64 {
        0.0
    }

    fn set_width(&mut self, _width: f64) {
        // do nothing
    }

    fn width(&self) -> f64 {
        0.0
    }

    fn set_solar_absorption_coeff(&mut self, _solar_absorption_coeff: f64) {
        // do nothing
    }

    fn solar_absorption_coeff(&self) -> f64 {
        0.0
    }
}

impl SolarRadiationInteractionNotExposed for BuildingElementAdjacentUnconditionedSpaceSimple {}

/// A type to represent ground building elements
///
/// There are various variable names used within this type with quite obscure names - the following is provided as a guide:
///
/// They have generally been included into the FloorData input type, whose purpose is to encapsulate the different
/// necessary fields dependent on the floor type.
///
/// * `floor_type`
///        - Slab_no_edge_insulation
///        - Slab_edge_insulation
///        - Suspended_floor
///        - Heated_basement
///        - Unheated_basement
/// * `edge_insulation`
///          - horizontal edge insulation
///          - vertical or external edge insulation
/// * `h_upper` - height of the floor upper surface, in m
///             average value is used if h varies
/// * `u_w` - thermal transmittance of walls above ground, in W/(m2·K)
///         in accordance with ISO 6946
/// * `u_f_s` - thermal transmittance of floor above basement), in W/(m2·K)
///           in accordance with ISO 6946
/// * `area_per_perimeter_vent` -  area of ventilation openings per perimeter, in m2/m
/// * `shield_fact_location` - wind shielding factor
///         - Sheltered
///         - Average
///         - Exposed
/// * `d_we` - thickness of the walls, in m
/// * `r_f_ins` - thermal resistance of insulation on base of underfloor space, in m2·K/W
/// * `z_b` - depth of basement floor below ground level, in m
/// * `r_w_b` - thermal resistance of walls of the basement, in m2·K/W
/// * `h_w` - height of the basement walls above ground level, in m
#[derive(Debug)]
pub(crate) struct BuildingElementGround {
    area: f64,
    total_area: f64,
    pitch: f64,
    external_pitch: f64,
    u_value: f64,
    external_conditions: Arc<ExternalConditions>,
    perimeter: f64,
    floor_data: FloorData,
    r_c: f64,
    psi_wall_floor_junc: f64,
    temp_int_annual: f64,
    f_sky: f64,
    therm_rad_to_sky: f64,
    h_pli: [f64; 4],
    k_pli: [f64; 5],
    k_m: f64,
    h_pi: f64,
    h_pe: f64,
    h_ce: f64,
    d_eq: f64,
}

impl BuildingElementGround {
    /// Arguments (names based on those in BS EN ISO 52016-1:2017):
    /// * `total_area` -- total area (in m2) of the building element across entire dwelling.
    ///     If the Floor is devided among several zones,
    ///     this is the total area across all zones.
    /// * `area`     -- area (in m2) of this building element within the zone.
    /// * `pitch` -- tilt angle of the surface from horizontal, in degrees between 0 and 180,
    ///     where 0 means the external surface is facing up, 90 means the external
    ///     surface is vertical and 180 means the external surface is facing down
    /// * `u_value`  -- steady-state thermal transmittance of floor, including the effect of the ground, in W / (m2.K)
    ///     Calculated for the entire ground floor, even if it is distributed among several zones.
    /// * `thermal_resistance_floor_construction`      -- total thermal resistance of all layers in the floor construction, in (m2.K) / W
    /// * `areal_heat_capacity`      -- areal heat capacity of the ground floor element, in J / (m2.K)
    /// * `perimeter` -- perimeter of the floor, in metres. Calculated for the entire ground floor,
    ///     even if it is distributed among several zones.
    /// * `psi_wall_floor_junc` -- linear thermal transmittance of the junction between the floor and the walls, in W / (m.K)
    /// * `ext_cond` -- reference to ExternalConditions object
    /// * `simulation_time` -- reference to SimulationTime object
    /// * `mass_distribution_class` -- distribution of mass in building element, one of:
    ///         - 'I':  mass concentrated on internal side
    ///         - 'E':  mass concentrated on external side
    ///         - 'IE': mass divided over internal and external side
    ///         - 'D':  mass equally distributed
    ///         - 'M':  mass concentrated inside
    /// * `floor_type`
    ///         - Slab_no_edge_insulation
    ///         - Slab_edge_insulation
    ///         - Suspended_floor
    ///         - Heated_basement
    ///         - Unheated_basement
    /// * `edge_insulation`
    ///         - horizontal edge insulation
    ///         - vertical or external edge insulation
    /// * `height_upper_surface` - height of the floor upper surface, in m
    ///     average value is used if h varies
    /// * `u_w` - thermal transmittance of walls above ground, in W/(m2·K)
    ///     in accordance with ISO 6946
    /// * `u_f_s` - thermal transmittance of floor above basement), in W/(m2·K)
    ///     in accordance with ISO 6946
    /// * `area_per_perimeter_vent` -  area of ventilation openings per perimeter, in m2/m
    /// * `shield_fact_location` - wind shielding factor
    ///         - Sheltered
    ///         - Average
    ///         - Exposed
    /// * `d_we` - thickness of the walls, in m
    /// * `r_f_ins` - thermal resistance of insulation on base of underfloor space, in m2·K/W
    /// * `z_b` - depth of basement floor below ground level, in m
    /// * `r_w_b` - thermal resistance of walls of the basement, in m2·K/W
    /// * `h_w` - height of the basement walls above ground level, in m
    ///
    /// Other variables:
    /// * `f_sky` -- view factor to the sky (see BS EN ISO 52016-1:2017, section 6.5.6.3.6)
    /// * `h_ce`     -- external convective heat transfer coefficient, in W / (m2.K)
    /// * `h_re`     -- external radiative heat transfer coefficient, in W / (m2.K)
    /// * `thermal_resistance_construction`      -- thermal resistance of the ground floor element including the
    ///     effect of the ground, in m2.K / W
    /// * `r_gr`     -- thermal resistance of the fixed ground layer, in m2.K / W
    /// * `k_gr`     -- areal heat capacity of the fixed ground layer, in J / (m2.K)
    pub(crate) fn new(
        total_area: f64,
        area: f64,
        pitch: f64,
        u_value: f64,
        thermal_resistance_floor_construction: f64,
        areal_heat_capacity: f64,
        mass_distribution_class: MassDistributionClass,
        floor_data: &FloorData,
        d_we: f64,
        perimeter: f64,
        psi_wall_floor_junc: f64,
        external_conditions: Arc<ExternalConditions>,
    ) -> anyhow::Result<Self> {
        let mut new_ground = Self {
            total_area,
            area,
            u_value,
            pitch,
            external_conditions,
            external_pitch: Default::default(),
            perimeter,
            floor_data: floor_data.clone(),
            r_c: Default::default(),
            psi_wall_floor_junc,
            temp_int_annual: Default::default(),
            f_sky: Default::default(),
            therm_rad_to_sky: Default::default(),
            h_pli: Default::default(),
            k_pli: Default::default(),
            k_m: Default::default(),
            h_pi: Default::default(),
            h_pe: Default::default(),
            h_ce: Default::default(),
            d_eq: Default::default(),
        };

        // # Thermal properties of ground from BS EN ISO 13370:2017 Table 7
        // Use values for clay or silt (same as BR 443 and SAP 10)
        let thermal_conductivity = 1.5; // in W/(m.K)
        let heat_capacity_per_vol = 3000000.; // in J/(m3.K)

        // Periodic penetration depth of ground from BS EN ISO 13370:2017 Table H.1
        // Use values for clay or silt (same as BR 443 and SAP 10)
        let periodic_penetration_depth = 2.2; // in m

        // Calculate thermal resistance and heat capacity of fixed ground layer
        // using BS EN ISO 13370:2017
        let thickness_ground_layer = 0.5; // in m. Specified in BS EN ISO 52016-1:2017 section 6.5.8.2

        // thermal resistance in (m2.K)/W
        let r_gr = thickness_ground_layer / thermal_conductivity;
        // areal heat capacity in J/(m2.K)
        let k_gr = thickness_ground_layer * heat_capacity_per_vol;

        // Calculate thermal resistance of virtual layer using BS EN ISO 13370:2017 Equation (F1)
        let r_si = 0.17; // ISO 6946 - internal surface resistance
        let r_vi = calculate_thermal_resistance_of_virtual_layer(
            u_value,
            thermal_resistance_floor_construction,
        )?;

        new_ground.init_heat_transfer_through_3_plus_2_nodes(
            thermal_resistance_floor_construction,
            r_gr,
            mass_distribution_class,
            k_gr,
            areal_heat_capacity,
        );
        new_ground.init_heat_transfer_other_side_ground(
            r_vi,
            d_we,
            thermal_conductivity,
            r_si,
            thermal_resistance_floor_construction,
            floor_data,
            periodic_penetration_depth,
            total_area,
            perimeter,
            u_value,
            psi_wall_floor_junc,
        )?;
        new_ground.init_solar_radiation_interaction(pitch, None, None, 0.0, 0.0, 0.0, 0.0);

        Ok(new_ground)
    }

    pub(crate) fn fabric_heat_loss(&self) -> f64 {
        self.area * self.u_value
    }

    pub(crate) fn h_ce(&self) -> f64 {
        HeatTransferOtherSideGround::h_ce(self)
    }

    pub(crate) fn h_re(&self) -> f64 {
        HeatTransferOtherSideGround::h_re(self)
    }

    #[cfg(test)]
    pub(crate) fn r_si(&self) -> f64 {
        HeatTransferInternal::r_si(self)
    }

    fn temp_ext(&self, simtime: SimulationTimeIteration) -> f64 {
        HeatTransferOtherSideGround::temp_ext(self, simtime)
    }
}

impl HeatTransferInternal for BuildingElementGround {
    fn pitch(&self) -> f64 {
        self.pitch
    }
}

impl HeatTransferInternalCommon for BuildingElementGround {}

impl HeatTransferThrough for BuildingElementGround {
    fn set_r_c(&mut self, thermal_resistance_construction: f64) {
        self.r_c = thermal_resistance_construction;
    }

    fn r_c(&self) -> f64 {
        self.r_c
    }

    fn k_m(&self) -> f64 {
        self.k_m
    }

    fn set_k_m(&mut self, areal_heat_capacity: f64) {
        self.k_m = areal_heat_capacity;
    }

    fn k_pli(&self) -> &[f64] {
        &self.k_pli
    }

    fn r_se(&self) -> f64 {
        // upstream Python uses duck typing/ lookups to find which method this is, but Rust needs to be explicit
        <dyn HeatTransferOtherSide>::r_se(self)
    }

    fn r_si(&self) -> f64 {
        <dyn HeatTransferInternal>::r_si(self)
    }

    fn area(&self) -> f64 {
        self.area
    }

    fn h_pli(&self) -> &[f64] {
        &self.h_pli
    }
}

impl HeatTransferThrough3Plus2Nodes for BuildingElementGround {
    fn set_h_pli(&mut self, h_pli: [f64; 4]) {
        self.h_pli = h_pli
    }

    fn set_k_pli(&mut self, k_pli: [f64; 5]) {
        self.k_pli = k_pli
    }
}

impl HeatTransferOtherSide for BuildingElementGround {
    fn set_f_sky(&mut self, f_sky: f64) {
        self.f_sky = f_sky;
    }

    fn set_therm_rad_to_sky(&mut self, therm_rad_to_sky: f64) {
        self.therm_rad_to_sky = therm_rad_to_sky;
    }

    fn therm_rad_to_sky(&self) -> f64 {
        self.therm_rad_to_sky
    }

    fn f_sky(&self) -> f64 {
        self.f_sky
    }

    fn external_conditions(&self) -> &ExternalConditions {
        self.external_conditions.as_ref()
    }
}

impl HeatTransferOtherSideGround for BuildingElementGround {
    fn set_temp_int_annual(&mut self, temp_int_annual: f64) {
        self.temp_int_annual = temp_int_annual;
    }

    fn temp_int_annual(&self) -> f64 {
        self.temp_int_annual
    }

    fn set_total_area(&mut self, total_area: f64) {
        self.total_area = total_area;
    }

    fn total_area(&self) -> f64 {
        self.total_area
    }

    fn set_perimeter(&mut self, perimeter: f64) {
        self.perimeter = perimeter;
    }

    fn perimeter(&self) -> f64 {
        self.perimeter
    }

    fn set_u_value(&mut self, u_value: f64) {
        self.u_value = u_value;
    }

    fn set_psi_wall_floor_junc(&mut self, psi_wall_floor_junc: f64) {
        self.psi_wall_floor_junc = psi_wall_floor_junc;
    }

    fn psi_wall_floor_junc(&self) -> f64 {
        self.psi_wall_floor_junc
    }

    fn set_floor_data(&mut self, floor_data: FloorData) {
        self.floor_data = floor_data;
    }

    fn set_h_ce(&mut self, h_ce: f64) {
        self.h_ce = h_ce;
    }

    fn set_d_eq(&mut self, d_eq: f64) {
        self.d_eq = d_eq;
    }

    fn set_h_pi(&mut self, h_pi: f64) {
        self.h_pi = h_pi;
    }

    fn h_pi(&self) -> f64 {
        self.h_pi
    }

    fn set_h_pe(&mut self, h_pe: f64) {
        self.h_pe = h_pe;
    }

    fn h_pe(&self) -> f64 {
        self.h_pe
    }

    fn u_value(&self) -> f64 {
        self.u_value
    }

    fn h_ce(&self) -> f64 {
        self.h_ce
    }
}

impl SolarRadiationInteraction for BuildingElementGround {
    fn set_external_pitch(&mut self, pitch: f64) {
        self.external_pitch = pitch;
    }

    fn external_pitch(&self) -> f64 {
        self.external_pitch
    }

    fn set_orientation(&mut self, _orientation: Orientation360) {
        // do nothing
    }

    fn orientation(&self) -> Option<Orientation360> {
        Default::default()
    }

    fn set_shading(&mut self, _shading: Option<Vec<WindowShadingObject>>) {
        // do nothing
    }

    fn shading(&self) -> &[WindowShadingObject] {
        &[]
    }

    fn set_base_height(&mut self, _base_height: f64) {
        // do nothing
    }

    fn base_height(&self) -> f64 {
        0.0
    }

    fn set_projected_height(&mut self, _projected_height: f64) {
        // do nothing
    }

    fn projected_height(&self) -> f64 {
        0.0
    }

    fn set_width(&mut self, _width: f64) {
        // do nothing
    }

    fn width(&self) -> f64 {
        0.0
    }

    fn set_solar_absorption_coeff(&mut self, _solar_absorption_coeff: f64) {
        // do nothing
    }

    fn solar_absorption_coeff(&self) -> f64 {
        0.0
    }
}

impl SolarRadiationInteractionNotExposed for BuildingElementGround {}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum WindowTreatmentControl {
    Manual,
    ManualMotorised,
    AutoMotorised,
    CombinedLightBlindHvac,
}

impl WindowTreatmentControl {
    pub(crate) fn is_manual(&self) -> bool {
        [Self::Manual, Self::ManualMotorised].contains(self)
    }

    pub(crate) fn is_automatic(&self) -> bool {
        !self.is_manual()
    }
}

impl From<WindowTreatmentControlInput> for WindowTreatmentControl {
    fn from(input: WindowTreatmentControlInput) -> Self {
        match input {
            WindowTreatmentControlInput::Manual => WindowTreatmentControl::Manual,
            WindowTreatmentControlInput::ManualMotorised => WindowTreatmentControl::ManualMotorised,
            WindowTreatmentControlInput::AutoMotorised => WindowTreatmentControl::AutoMotorised,
            WindowTreatmentControlInput::CombinedLightBlindHvac => {
                WindowTreatmentControl::CombinedLightBlindHvac
            }
        }
    }
}

#[derive(Debug)]
pub(crate) struct WindowTreatment {
    _treatment_type: WindowTreatmentType,
    controls: WindowTreatmentControl,
    #[allow(dead_code)]
    delta_r: f64,
    trans_red: f64,
    closing_irradiance_control: Option<Arc<Control>>,
    opening_irradiance_control: Option<Arc<Control>>,
    open_control: Option<Arc<Control>>,
    is_open: AtomicBool,
    opening_delay_hrs: f64,
    time_last_adjusted: AtomicF64,
}

impl WindowTreatment {
    pub(crate) fn from_input(
        input: &WindowTreatmentInput,
        controls: &Controls,
        current_hour: u32,
    ) -> Self {
        Self {
            _treatment_type: input.treatment_type,
            controls: input.controls.into(),
            delta_r: input.delta_r,
            trans_red: input.trans_red,
            closing_irradiance_control: input
                .control_closing_irrad
                .as_ref()
                .and_then(|ctrl| controls.get_with_string(ctrl)),
            opening_irradiance_control: input
                .control_opening_irrad
                .as_ref()
                .and_then(|ctrl| controls.get_with_string(ctrl)),
            open_control: input
                .control_open
                .as_ref()
                .and_then(|ctrl| controls.get_with_string(ctrl)),
            is_open: input.is_open.unwrap_or_default().into(),
            opening_delay_hrs: input.opening_delay_hrs,
            time_last_adjusted: (current_hour as f64).into(),
        }
    }
}

/// A type to represent transparent building elements (windows etc.)
#[derive(Debug)]
pub(crate) struct BuildingElementTransparent {
    area: f64,
    _mid_height: f64,
    g_value: f64,
    external_conditions: Arc<ExternalConditions>,
    frame_area_fraction: f64,
    pitch: f64,
    external_pitch: f64,
    orientation: Option<Orientation360>,
    base_height: f64,
    projected_height: f64,
    width: f64,
    r_c: f64,
    k_m: f64,
    h_pli: [f64; 1],
    k_pli: [f64; 2],
    f_sky: f64,
    therm_rad_to_sky: f64,
    shading: Vec<WindowShadingObject>,
    treatment: Vec<WindowTreatment>,
}

impl BuildingElementTransparent {
    /// Arguments (names based on those in BS EN ISO 52016-1:2017):
    /// * `pitch` - tilt angle of the surface from horizontal, in degrees between 0 and 180,
    ///          where 0 means the external surface is facing up, 90 means the external
    ///          surface is vertical and 180 means the external surface is facing down
    /// * `thermal_resistance_construction` - thermal resistance, in m2.K / W
    /// * `orientation` - is the orientation angle of the inclined surface, expressed
    ///                as the geographical azimuth angle of the horizontal projection
    ///                of the inclined surface normal, -180 to 180, in degrees
    /// * `base_height` - is the distance between the ground and the lowest edge of the element, in m
    /// * `height`      - is the height of the building element, in m
    /// * `width`       - is the width of the building element, in m
    /// * `g_value` - total solar energy transmittance of the transparent part of the window
    /// * `frame_area_fraction` - is the frame area fraction of window wi, ratio of the
    ///                        projected frame area to the overall projected area of
    ///                        the glazed element of the window
    /// * `treatment` - is additional window elements such as curtain or blinds
    /// * `external_conditions` -- reference to ExternalConditions object
    pub(crate) fn new(
        pitch: f64,
        thermal_resistance_construction: f64,
        orientation: Option<Orientation360>,
        g_value: f64,
        frame_area_fraction: f64,
        base_height: f64,
        height: f64,
        width: f64,
        shading: Option<Vec<WindowShadingObject>>,
        treatment: Vec<WindowTreatment>,
        external_conditions: Arc<ExternalConditions>,
    ) -> Self {
        let mut new_trans = Self {
            _mid_height: base_height + height / 2.,
            g_value,
            external_conditions,
            frame_area_fraction,
            area: calculate_area(height, width),
            pitch,
            external_pitch: Default::default(),
            orientation,
            base_height,
            projected_height: Default::default(),
            width,
            r_c: thermal_resistance_construction,
            k_m: Default::default(),
            h_pli: Default::default(),
            k_pli: Default::default(),
            f_sky: Default::default(),
            therm_rad_to_sky: Default::default(),
            shading: Default::default(),
            treatment,
        };

        new_trans.init_heat_transfer_through_2_nodes(thermal_resistance_construction);
        new_trans.init_heat_transfer_other_side_outside(pitch);
        new_trans.init_solar_radiation_interaction(
            pitch,
            orientation,
            shading,
            base_height,
            projected_height(pitch, height),
            width,
            0.0,
        );

        new_trans
    }

    fn adjust_treatment(&self, simtime: SimulationTimeIteration) -> anyhow::Result<()> {
        fn open_treatment(treatment: &WindowTreatment, current_time: f64) {
            if !treatment.is_open.load(Ordering::SeqCst) {
                treatment.is_open.store(true, Ordering::SeqCst);
                treatment
                    .time_last_adjusted
                    .store(current_time, Ordering::SeqCst);
            }
        }

        fn close_treatment(treatment: &WindowTreatment, current_time: f64) {
            if treatment.is_open.load(Ordering::SeqCst) {
                treatment.is_open.store(false, Ordering::SeqCst);
                treatment
                    .time_last_adjusted
                    .store(current_time, Ordering::SeqCst);
            }
        }

        // Operation and control logic for window treatments (windows, blinds, etc.)
        // as per Annex G and Tables B.23 and B.24 in BS EN ISO 52016-1:2017.
        // 'Manual' modes for curtains/shutters also requires occupancy driver, however
        // only time and solar based controls for now.
        // 'trans_red' specific to selected treatment as per BS EN 13125:2001.
        if !self.treatment.is_empty() {
            let surf_irrad = self.external_conditions.surface_irradiance(
                self.base_height,
                self.projected_height,
                self.width,
                self.pitch,
                self.orientation.ok_or_else(|| {
                    anyhow!("Cannot compute surface irradiance because orientation is None.")
                })?,
                self.shading(),
                simtime,
            )?;
            let time_current = simtime.time;
            for treatment in self.treatment.iter() {
                let ctrl_open: Option<bool> = treatment
                    .open_control
                    .as_ref()
                    .map(|ctrl| ctrl.is_on(&simtime));
                let closing_irrad_threshold: Option<f64> = treatment
                    .closing_irradiance_control
                    .as_ref()
                    .and_then(|ctrl| ctrl.setpnt(&simtime));
                let opening_irrad_threshold: Option<f64> = treatment
                    .opening_irradiance_control
                    .as_ref()
                    .and_then(|ctrl| ctrl.setpnt(&simtime));

                match ctrl_open {
                    Some(true) => open_treatment(treatment, time_current),
                    Some(false) => close_treatment(treatment, time_current),
                    None => {
                        if closing_irrad_threshold.is_some()
                            && closing_irrad_threshold.is_some_and(|closing_irrad_threshold| {
                                surf_irrad > closing_irrad_threshold
                            })
                        {
                            close_treatment(treatment, time_current);
                        } else if opening_irrad_threshold.is_some()
                            && opening_irrad_threshold.is_some_and(|opening_irrad_threshold| {
                                surf_irrad < opening_irrad_threshold
                            })
                            && (treatment.controls.is_manual()
                                || (treatment.controls.is_automatic()
                                    && time_current
                                        - treatment.time_last_adjusted.load(Ordering::SeqCst)
                                        >= treatment.opening_delay_hrs))
                        {
                            open_treatment(treatment, time_current);
                        }
                    }
                }
            }
        }

        Ok(())
    }

    pub(crate) fn solar_gains(&self, simtime: SimulationTimeIteration) -> anyhow::Result<f64> {
        let mut solar_gains = SolarRadiationInteractionTransmitted::solar_gains(self, simtime)?;

        if !self.treatment.is_empty() {
            self.adjust_treatment(simtime)?;
            for treatment in self.treatment.iter() {
                if !treatment.is_open.load(Ordering::SeqCst) {
                    solar_gains -= solar_gains * treatment.trans_red;
                }
            }
        }

        Ok(solar_gains)
    }

    pub(crate) fn fabric_heat_loss(&self) -> f64 {
        // Effective window U-value includes assumed use of curtains/blinds, see
        // SAP10.2 spec, paragraph 3.2
        // TODO (from Python) Confirm this is still the desired approach for SAP 11 (sic from Python)
        let r_curtains_blinds = 0.04;
        // Add standard surface resistances to resistance of construction when calculating U-value
        let u_value = 1.0 / ((self.r_c + self.r_se() + self.r_si()) + r_curtains_blinds);
        self.area * u_value
    }

    pub(crate) fn h_pli_by_index(
        &self,
        idx: usize,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // Account for resistance of window treatment in heat transfer coefficient
        // TODO (from Python) Check that idx refers to inside surface?
        let mut thermal_resistance_construction = 1.0
            / self.h_pli.get(idx).ok_or_else(|| {
                anyhow!("Could not get h_pli value in transparent building element for index {idx}")
            })?;
        for treatment in self.treatment.iter() {
            self.adjust_treatment(simtime)?;
            if !treatment.is_open.load(Ordering::SeqCst) {
                thermal_resistance_construction += treatment.delta_r;
            }
        }

        Ok(1.0 / thermal_resistance_construction)
    }

    pub(crate) fn r_se(&self) -> f64 {
        HeatTransferOtherSide::r_se(self)
    }

    pub(crate) fn r_si(&self) -> f64 {
        HeatTransferInternal::r_si(self)
    }

    pub(crate) fn area(&self) -> f64 {
        self.area
    }

    pub(crate) fn shading_factors_direct_diffuse(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        SolarRadiationInteractionTransmitted::shading_factors_direct_diffuse(self, simtime)
    }
}

impl HeatTransferInternal for BuildingElementTransparent {
    fn pitch(&self) -> f64 {
        self.pitch
    }
}

impl HeatTransferInternalCommon for BuildingElementTransparent {}

impl HeatTransferThrough for BuildingElementTransparent {
    fn set_r_c(&mut self, thermal_resistance_construction: f64) {
        self.r_c = thermal_resistance_construction;
    }

    fn r_c(&self) -> f64 {
        self.r_c
    }

    fn k_m(&self) -> f64 {
        self.k_m
    }

    fn set_k_m(&mut self, areal_heat_capacity: f64) {
        self.k_m = areal_heat_capacity;
    }

    fn k_pli(&self) -> &[f64] {
        &self.k_pli
    }

    fn r_se(&self) -> f64 {
        // upstream Python uses duck typing/ lookups to find which method this is, but Rust needs to be explicit
        <dyn HeatTransferOtherSide>::r_se(self)
    }

    fn r_si(&self) -> f64 {
        <dyn HeatTransferInternal>::r_si(self)
    }

    fn area(&self) -> f64 {
        self.area
    }

    fn h_pli(&self) -> &[f64] {
        &self.h_pli
    }
}

impl HeatTransferThrough2Nodes for BuildingElementTransparent {
    fn set_h_pli(&mut self, h_pli: [f64; 1]) {
        self.h_pli = h_pli;
    }

    fn set_k_pli(&mut self, k_pli: [f64; 2]) {
        self.k_pli = k_pli;
    }
}

impl HeatTransferOtherSide for BuildingElementTransparent {
    fn set_f_sky(&mut self, f_sky: f64) {
        self.f_sky = f_sky;
    }

    fn f_sky(&self) -> f64 {
        self.f_sky
    }

    fn set_therm_rad_to_sky(&mut self, therm_rad_to_sky: f64) {
        self.therm_rad_to_sky = therm_rad_to_sky;
    }

    fn therm_rad_to_sky(&self) -> f64 {
        self.therm_rad_to_sky
    }

    fn external_conditions(&self) -> &ExternalConditions {
        self.external_conditions.as_ref()
    }
}

impl HeatTransferOtherSideOutside for BuildingElementTransparent {}

impl SolarRadiationInteraction for BuildingElementTransparent {
    fn set_external_pitch(&mut self, pitch: f64) {
        self.external_pitch = pitch;
    }

    fn external_pitch(&self) -> f64 {
        self.external_pitch
    }

    fn set_orientation(&mut self, orientation: Orientation360) {
        self.orientation.replace(orientation);
    }

    fn orientation(&self) -> Option<Orientation360> {
        self.orientation
    }

    fn set_shading(&mut self, shading: Option<Vec<WindowShadingObject>>) {
        self.shading = shading.unwrap_or_default();
    }

    fn shading(&self) -> &[WindowShadingObject] {
        self.shading.as_slice()
    }

    fn set_base_height(&mut self, base_height: f64) {
        self.base_height = base_height;
    }

    fn base_height(&self) -> f64 {
        self.base_height
    }

    fn set_projected_height(&mut self, projected_height: f64) {
        self.projected_height = projected_height;
    }

    fn projected_height(&self) -> f64 {
        self.projected_height
    }

    fn set_width(&mut self, width: f64) {
        self.width = width;
    }

    fn width(&self) -> f64 {
        self.width
    }

    fn set_solar_absorption_coeff(&mut self, _solar_absorption_coeff: f64) {
        // do nothing
    }

    fn solar_absorption_coeff(&self) -> f64 {
        0.0
    }
}

impl SolarRadiationInteractionTransmitted for BuildingElementTransparent {
    fn unconverted_g_value(&self) -> f64 {
        self.g_value
    }

    fn external_conditions(&self) -> &ExternalConditions {
        self.external_conditions.as_ref()
    }

    fn pitch(&self) -> f64 {
        self.pitch
    }

    fn area(&self) -> f64 {
        self.area
    }

    fn frame_area_fraction(&self) -> f64 {
        self.frame_area_fraction
    }
}

pub fn pitch_class(pitch: f64) -> HeatFlowDirection {
    match pitch {
        PITCH_LIMIT_HORIZ_CEILING..=PITCH_LIMIT_HORIZ_FLOOR => HeatFlowDirection::Horizontal,
        ..PITCH_LIMIT_HORIZ_CEILING => HeatFlowDirection::Upwards,
        PITCH_LIMIT_HORIZ_FLOOR.. => HeatFlowDirection::Downwards,
        _ => unreachable!("Rust cannot tell that above is exhaustive"),
    }
}

pub(crate) fn convert_uvalue_to_resistance(u_value: f64, pitch: f64) -> f64 {
    (1.0 / u_value) - r_si_for_pitch(pitch) - R_SE
}

fn r_si_for_pitch(pitch: f64) -> f64 {
    match pitch {
        _ if (PITCH_LIMIT_HORIZ_CEILING..=PITCH_LIMIT_HORIZ_FLOOR).contains(&pitch) => {
            R_SI_HORIZONTAL
        }
        _ if pitch < PITCH_LIMIT_HORIZ_CEILING => R_SI_UPWARDS,
        _ if pitch > PITCH_LIMIT_HORIZ_FLOOR => R_SI_DOWNWARDS,
        _ => panic!("problem with pitch value"), // this case should never happen as above cases are exhaustive but rust can't tell
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::controls::time_control::{OnOffTimeControl, SetpointTimeControl};
    use crate::external_conditions::{DaylightSavingsConfig, ShadingSegment};
    use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
    use approx::assert_relative_eq;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[rstest]
    fn test_convert_uvalue_to_resistance() {
        assert_relative_eq!(convert_uvalue_to_resistance(2., 40.), 0.35985829616804244);
    }

    struct MockHeatTransferInternal(f64);

    impl HeatTransferInternal for MockHeatTransferInternal {
        fn pitch(&self) -> f64 {
            self.0
        }
    }

    #[fixture]
    fn heat_transfer_internal_a() -> MockHeatTransferInternal {
        MockHeatTransferInternal(0.)
    }

    #[fixture]
    fn heat_transfer_internal_b() -> MockHeatTransferInternal {
        MockHeatTransferInternal(45.)
    }

    #[fixture]
    fn heat_transfer_internal_c() -> MockHeatTransferInternal {
        MockHeatTransferInternal(90.)
    }

    #[fixture]
    fn heat_transfer_internal_d() -> MockHeatTransferInternal {
        MockHeatTransferInternal(180.)
    }

    #[rstest]
    fn test_heat_flow_direction_for_heat_transfer_internal(
        heat_transfer_internal_b: impl HeatTransferInternal,
        heat_transfer_internal_c: impl HeatTransferInternal,
        heat_transfer_internal_d: impl HeatTransferInternal,
    ) {
        assert_eq!(
            heat_transfer_internal_b.heat_flow_direction(20., 25.),
            HeatFlowDirection::Downwards
        );
        assert_eq!(
            heat_transfer_internal_c.heat_flow_direction(20., 25.),
            HeatFlowDirection::Horizontal
        );
        assert_eq!(
            heat_transfer_internal_d.heat_flow_direction(20., 25.),
            HeatFlowDirection::Upwards
        );
    }

    #[rstest]
    fn test_convert_uvalue_to_resistance_for_heat_transfer_internal_a(
        heat_transfer_internal_a: impl HeatTransferInternal,
    ) {
        assert_relative_eq!(
            heat_transfer_internal_a.convert_uvalue_to_resistance(1., 180.),
            0.7870483926665635,
            max_relative = 1e-8
        );
    }

    #[rstest]
    fn test_r_si(heat_transfer_internal_a: impl HeatTransferInternal) {
        assert_relative_eq!(
            heat_transfer_internal_a.r_si(),
            0.0987166831194472,
            max_relative = 1e-8
        );
    }

    #[rstest]
    fn test_r_si_pitch(heat_transfer_internal_a: impl HeatTransferInternal) {
        // r_si_horizontal
        assert_relative_eq!(
            heat_transfer_internal_a.r_si_with_pitch(90.0),
            0.1310615989515072,
            max_relative = 1e-8
        );

        // r_si_upwards
        assert_relative_eq!(
            heat_transfer_internal_a.r_si_with_pitch(30.0),
            0.0987166831194472,
            max_relative = 1e-8
        );

        // r_si_downwards
        assert_relative_eq!(
            heat_transfer_internal_a.r_si_with_pitch(150.0),
            0.17152658662092624,
            max_relative = 1e-8
        );
    }

    #[rstest]
    #[should_panic(expected = "internal error: entered unreachable code")]
    fn test_r_si_invalid_pitch(heat_transfer_internal_a: impl HeatTransferInternal) {
        heat_transfer_internal_a.r_si_with_pitch(f64::NAN);
    }

    #[rstest]
    fn test_pitch_class(heat_transfer_internal_a: impl HeatTransferInternal) {
        assert_eq!(
            heat_transfer_internal_a.pitch_class(90.0),
            HeatFlowDirection::Horizontal
        );
        assert_eq!(
            heat_transfer_internal_a.pitch_class(30.0),
            HeatFlowDirection::Upwards
        );
        assert_eq!(
            heat_transfer_internal_a.pitch_class(150.0),
            HeatFlowDirection::Downwards
        );
    }

    #[rstest]
    #[should_panic(expected = "internal error: entered unreachable code")]
    fn test_pitch_class_invalid(heat_transfer_internal_a: impl HeatTransferInternal) {
        heat_transfer_internal_a.pitch_class(f64::NAN);
    }

    #[rstest]
    fn test_h_ri(heat_transfer_internal_a: impl HeatTransferInternal) {
        assert_eq!(heat_transfer_internal_a.h_ri(), 5.13);
    }

    // skip two tests (`test_init_invalid_mass_distribution_class`) as Rust would not allow invalid mass_distribution_class values

    struct MockHeatTransferOtherSide(f64);

    impl HeatTransferOtherSide for MockHeatTransferOtherSide {
        fn set_f_sky(&mut self, f_sky: f64) {
            self.0 = f_sky;
        }

        fn set_therm_rad_to_sky(&mut self, _therm_rad_to_sky: f64) {
            unreachable!()
        }

        fn therm_rad_to_sky(&self) -> f64 {
            unreachable!()
        }

        fn f_sky(&self) -> f64 {
            self.0
        }

        fn external_conditions(&self) -> &ExternalConditions {
            unreachable!()
        }
    }

    #[fixture]
    fn heat_transfer_other_side_a() -> MockHeatTransferOtherSide {
        MockHeatTransferOtherSide(0.)
    }

    #[rstest]
    fn test_r_se(heat_transfer_other_side_a: impl HeatTransferOtherSide) {
        assert_relative_eq!(
            heat_transfer_other_side_a.r_se(),
            0.041425020712510356,
            max_relative = 1e-8
        );
    }

    #[rstest]
    fn test_h_ce(heat_transfer_other_side_a: impl HeatTransferOtherSide) {
        assert_eq!(heat_transfer_other_side_a.h_ce(), 20.0);
    }

    #[rstest]
    fn test_h_re_for_heat_transfer_other_side(
        heat_transfer_other_side_a: impl HeatTransferOtherSide,
    ) {
        assert_eq!(heat_transfer_other_side_a.h_re(), 4.14);
    }

    // skip test_fabric_heat_loss as Rust would not allow invalid method call

    #[fixture]
    fn simulation_time() -> SimulationTimeIterator {
        SimulationTime::new(0.0, 4.0, 1.0).iter()
    }

    #[fixture]
    fn external_conditions(simulation_time: SimulationTimeIterator) -> Arc<ExternalConditions> {
        Arc::new(ExternalConditions::new(
            &simulation_time,
            vec![0.0, 5.0, 10.0, 15.0],
            vec![],
            vec![],
            vec![0.0; 4],
            vec![0.0; 4],
            vec![],
            55.0,
            0.0,
            0,
            0,
            None,
            1.0,
            None,
            None,
            false,
            false,
            None,
        ))
    }

    #[fixture]
    fn be_i(external_conditions: Arc<ExternalConditions>) -> BuildingElementOpaque {
        BuildingElementOpaque::new(
            20.,
            false,
            180.,
            0.60,
            0.25,
            19000.0,
            MassDistributionClass::I,
            Orientation360::create_from_180(0.).unwrap().into(),
            0.,
            2.,
            10.,
            external_conditions,
        )
    }

    #[fixture]
    fn be_e(external_conditions: Arc<ExternalConditions>) -> BuildingElementOpaque {
        BuildingElementOpaque::new(
            22.5,
            false,
            135.,
            0.61,
            0.50,
            18000.0,
            MassDistributionClass::E,
            Orientation360::create_from_180(180.0).unwrap().into(),
            0.,
            2.25,
            10.,
            external_conditions,
        )
    }

    #[fixture]
    fn be_ie(external_conditions: Arc<ExternalConditions>) -> BuildingElementOpaque {
        BuildingElementOpaque::new(
            25.,
            false,
            90.,
            0.62,
            0.75,
            17000.0,
            MassDistributionClass::IE,
            Orientation360::create_from_180(90.).unwrap().into(),
            0.,
            2.5,
            10.,
            external_conditions,
        )
    }

    #[fixture]
    fn be_d(external_conditions: Arc<ExternalConditions>) -> BuildingElementOpaque {
        BuildingElementOpaque::new(
            27.5,
            true,
            45.,
            0.63,
            0.80,
            16000.0,
            MassDistributionClass::D,
            Orientation360::create_from_180(-90.0).unwrap().into(),
            0.,
            2.75,
            10.,
            external_conditions,
        )
    }

    #[fixture]
    fn be_m(external_conditions: Arc<ExternalConditions>) -> BuildingElementOpaque {
        BuildingElementOpaque::new(
            30.,
            false,
            0.,
            0.64,
            0.40,
            15000.0,
            MassDistributionClass::M,
            Orientation360::create_from_180(0.).unwrap().into(),
            0.,
            3.,
            10.,
            external_conditions,
        )
    }

    #[fixture]
    fn opaque_building_elements(
        be_i: BuildingElementOpaque,
        be_e: BuildingElementOpaque,
        be_ie: BuildingElementOpaque,
        be_d: BuildingElementOpaque,
        be_m: BuildingElementOpaque,
    ) -> [BuildingElementOpaque; 5] {
        [be_i, be_e, be_ie, be_d, be_m]
    }

    #[rstest]
    fn test_no_of_nodes_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        for be in opaque_building_elements.iter() {
            assert_eq!(be.number_of_nodes(), 5, "incorrect number of nodes");
            assert_eq!(
                be.number_of_inside_nodes(),
                3,
                "incorrect number of inside nodes"
            );
        }
    }

    #[rstest]
    fn test_area_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        // Define increment between test cases
        let area_inc = 2.5;
        for (i, be) in opaque_building_elements.iter().enumerate() {
            assert_eq!(
                be.area(),
                20.0 + i as f64 * area_inc,
                "incorrect area returned"
            );
        }
    }

    #[rstest]
    fn test_heat_flow_direction_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        let temp_int_air = 20.0;
        let temp_int_surface = [19.0, 21.0, 22.0, 21.0, 19.0];
        let results = [
            HeatFlowDirection::Downwards,
            HeatFlowDirection::Upwards,
            HeatFlowDirection::Horizontal,
            HeatFlowDirection::Downwards,
            HeatFlowDirection::Upwards,
        ];
        for (i, be) in opaque_building_elements.iter().enumerate() {
            assert_eq!(
                be.heat_flow_direction(temp_int_air, temp_int_surface[i]),
                results[i],
                "incorrect heat flow direction returned"
            );
        }
    }

    #[rstest]
    fn test_r_si_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        let results = [0.17, 0.17, 0.13, 0.10, 0.10];

        for (i, be) in opaque_building_elements.iter().enumerate() {
            assert_relative_eq!(be.r_si(), results[i], max_relative = 0.05);
        }
    }

    #[rstest]
    fn test_h_ci_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        let temp_int_air = 20.0;
        let temp_int_surface = [19.0, 21.0, 22.0, 21.0, 19.0];
        let results = [0.7, 5.0, 2.5, 0.7, 5.0];

        for (i, be) in opaque_building_elements.iter().enumerate() {
            assert_relative_eq!(
                be.h_ci(temp_int_air, temp_int_surface[i]),
                results[i],
                max_relative = 1e-2
            );
        }
    }

    #[rstest]
    fn test_h_ri_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        for be in opaque_building_elements.iter() {
            assert_relative_eq!(be.h_ri(), 5.13,);
        }
    }

    #[rstest]
    fn test_h_ce_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        for be in opaque_building_elements.iter() {
            assert_relative_eq!(be.h_ce(), 20.0,);
        }
    }

    #[rstest]
    fn test_h_re(opaque_building_elements: [BuildingElementOpaque; 5]) {
        for be in opaque_building_elements.iter() {
            assert_relative_eq!(be.h_re(), 4.14,);
        }
    }

    #[rstest]
    fn test_solar_absorption_coeff_for_opaque(
        opaque_building_elements: [BuildingElementOpaque; 5],
    ) {
        // Define increment between test cases
        let solar_absorption_coeff_inc = 0.01;

        for (i, be) in opaque_building_elements.iter().enumerate() {
            assert_relative_eq!(
                be.solar_absorption_coeff(),
                0.6 + i as f64 * solar_absorption_coeff_inc,
            );
        }
    }

    #[rstest]
    fn test_therm_rad_to_sky_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        let results = [0.0, 6.6691785923823135, 22.77, 38.87082140761768, 45.54];

        for (i, be) in opaque_building_elements.iter().enumerate() {
            assert_relative_eq!(be.therm_rad_to_sky(), results[i],);
        }
    }

    #[rstest]
    fn test_h_pli_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        let results = [
            [24.0, 12.0, 12.0, 24.0],
            [12.0, 6.0, 6.0, 12.0],
            [8.0, 4.0, 4.0, 8.0],
            [7.5, 3.75, 3.75, 7.5],
            [15.0, 7.5, 7.5, 15.0],
        ];

        for (i, be) in opaque_building_elements.iter().enumerate() {
            assert_eq!(be.h_pli(), &results[i], "incorrect h_pli list returned");
        }
    }

    #[rstest]
    pub fn test_k_pli_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        let results = [
            [0.0, 0.0, 0.0, 0.0, 19000.0],
            [18000.0, 0.0, 0.0, 0.0, 0.0],
            [8500.0, 0.0, 0.0, 0.0, 8500.0],
            [2000.0, 4000.0, 4000.0, 4000.0, 2000.0],
            [0.0, 0.0, 15000.0, 0.0, 0.0],
        ];

        for (i, be) in opaque_building_elements.iter().enumerate() {
            assert_eq!(be.k_pli(), &results[i], "incorrect k_pli returned");
        }
    }

    #[rstest]
    fn test_temp_ext_for_opaque(
        opaque_building_elements: [BuildingElementOpaque; 5],
        simulation_time: SimulationTimeIterator,
    ) {
        for be in opaque_building_elements.iter() {
            for (t_idx, t_it) in simulation_time.clone().enumerate() {
                assert_eq!(
                    be.temp_ext(t_it),
                    t_idx as f64 * 5.,
                    "incorrect ext temp returned"
                );
            }
        }
    }

    #[rstest]
    fn test_fabric_heat_loss_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        let results = [43.20, 31.56, 27.10, 29.25, 55.54];
        for (i, be) in opaque_building_elements.iter().enumerate() {
            assert_relative_eq!(be.fabric_heat_loss(), results[i], max_relative = 1e-2);
        }
    }

    #[rstest]
    fn test_heat_capacity_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        let results = [380., 405., 425., 440., 450.];
        for (i, be) in opaque_building_elements.iter().enumerate() {
            assert_eq!(
                be.heat_capacity(),
                results[i],
                "incorrect heat capacity returned"
            );
        }
    }

    #[fixture]
    fn conditioned_space_adjacent_building_elements(
        external_conditions: Arc<ExternalConditions>,
    ) -> [BuildingElementAdjacentConditionedSpace; 5] {
        let be_i = BuildingElementAdjacentConditionedSpace::new(
            20.0,
            180.,
            0.25,
            19000.0,
            MassDistributionClass::I,
            external_conditions.clone(),
        );
        let be_e = BuildingElementAdjacentConditionedSpace::new(
            22.5,
            135.,
            0.50,
            18000.0,
            MassDistributionClass::E,
            external_conditions.clone(),
        );
        let be_ie = BuildingElementAdjacentConditionedSpace::new(
            25.0,
            90.,
            0.75,
            17000.0,
            MassDistributionClass::IE,
            external_conditions.clone(),
        );
        let be_d = BuildingElementAdjacentConditionedSpace::new(
            27.5,
            45.,
            0.80,
            16000.0,
            MassDistributionClass::D,
            external_conditions.clone(),
        );
        let be_m = BuildingElementAdjacentConditionedSpace::new(
            30.0,
            0.,
            0.40,
            15000.0,
            MassDistributionClass::M,
            external_conditions,
        );
        [be_i, be_e, be_ie, be_d, be_m]
    }

    #[rstest]
    fn test_no_of_nodes_for_conditioned_space_adjacent_elements(
        conditioned_space_adjacent_building_elements: [BuildingElementAdjacentConditionedSpace; 5],
    ) {
        for be in conditioned_space_adjacent_building_elements {
            assert_eq!(be.number_of_nodes(), 5, "incorrect number of nodes");
            assert_eq!(
                be.number_of_inside_nodes(),
                3,
                "incorrect number of inside nodes"
            );
        }
    }

    #[rstest]
    fn test_area_for_conditioned_space_adjacent_elements(
        conditioned_space_adjacent_building_elements: [BuildingElementAdjacentConditionedSpace; 5],
    ) {
        // Define increment between test cases
        let area_inc = 2.5;

        for (i, be) in conditioned_space_adjacent_building_elements
            .iter()
            .enumerate()
        {
            assert_relative_eq!(be.area(), 20.0 + i as f64 * area_inc,);
        }
    }

    #[rstest]
    fn test_heat_flow_direction_for_conditioned_space_adjacent_elements(
        conditioned_space_adjacent_building_elements: [BuildingElementAdjacentConditionedSpace; 5],
    ) {
        let temp_int_air = 20.0;
        let temp_int_surface = [19.0, 21.0, 22.0, 21.0, 19.0];
        let results = [
            HeatFlowDirection::Downwards,
            HeatFlowDirection::Upwards,
            HeatFlowDirection::Horizontal,
            HeatFlowDirection::Downwards,
            HeatFlowDirection::Upwards,
        ];

        for (i, be) in conditioned_space_adjacent_building_elements
            .iter()
            .enumerate()
        {
            assert_eq!(
                be.heat_flow_direction(temp_int_air, temp_int_surface[i]),
                results[i],
                "incorrect heat flow direction returned"
            );
        }
    }

    #[rstest]
    fn test_r_si_for_conditioned_space_adjacent_elements(
        conditioned_space_adjacent_building_elements: [BuildingElementAdjacentConditionedSpace; 5],
    ) {
        let results = [0.17, 0.17, 0.13, 0.10, 0.10];

        for (i, be) in conditioned_space_adjacent_building_elements
            .iter()
            .enumerate()
        {
            assert_relative_eq!(be.r_si(), results[i], max_relative = 0.05);
        }
    }

    #[rstest]
    fn test_h_ci_for_conditioned_space_adjacent_elements(
        conditioned_space_adjacent_building_elements: [BuildingElementAdjacentConditionedSpace; 5],
    ) {
        let temp_int_air = 20.0;
        let temp_int_surface = [19.0, 21.0, 22.0, 21.0, 19.0];
        let results = [0.7, 5.0, 2.5, 0.7, 5.0];

        for (i, be) in conditioned_space_adjacent_building_elements
            .iter()
            .enumerate()
        {
            assert_relative_eq!(be.h_ci(temp_int_air, temp_int_surface[i]), results[i],);
        }
    }

    #[rstest]
    fn test_h_ri_for_conditioned_space_adjacent_elements(
        conditioned_space_adjacent_building_elements: [BuildingElementAdjacentConditionedSpace; 5],
    ) {
        for be in conditioned_space_adjacent_building_elements {
            assert_relative_eq!(be.h_ri(), 5.13,);
        }
    }

    #[rstest]
    fn test_h_ce_for_conditioned_space_adjacent_elements(
        conditioned_space_adjacent_building_elements: [BuildingElementAdjacentConditionedSpace; 5],
    ) {
        for be in conditioned_space_adjacent_building_elements {
            assert_relative_eq!(be.h_ce(), 0.0,);
        }
    }

    #[rstest]
    fn test_h_re_for_conditioned_space_adjacent_elements(
        conditioned_space_adjacent_building_elements: [BuildingElementAdjacentConditionedSpace; 5],
    ) {
        for be in conditioned_space_adjacent_building_elements {
            assert_eq!(be.h_re(), 0.0, "incorrect h_re returned");
        }
    }

    #[rstest]
    fn test_solar_absorption_coeff_for_conditioned_space_adjacent_elements(
        conditioned_space_adjacent_building_elements: [BuildingElementAdjacentConditionedSpace; 5],
    ) {
        for be in conditioned_space_adjacent_building_elements {
            assert_eq!(
                be.solar_absorption_coeff(),
                0.0,
                "incorrect solar_absorption_coeff returned"
            );
        }
    }

    #[rstest]
    fn test_therm_rad_to_sky_for_conditioned_space_adjacent_elements(
        conditioned_space_adjacent_building_elements: [BuildingElementAdjacentConditionedSpace; 5],
    ) {
        for be in conditioned_space_adjacent_building_elements {
            assert_eq!(
                be.therm_rad_to_sky(),
                0.0,
                "incorrect therm_rad_to_sky returned"
            );
        }
    }

    #[rstest]
    fn test_h_pli_for_conditioned_space_adjacent_elements(
        conditioned_space_adjacent_building_elements: [BuildingElementAdjacentConditionedSpace; 5],
    ) {
        let results = [
            [24.0, 12.0, 12.0, 24.0],
            [12.0, 6.0, 6.0, 12.0],
            [8.0, 4.0, 4.0, 8.0],
            [7.5, 3.75, 3.75, 7.5],
            [15.0, 7.5, 7.5, 15.0],
        ];
        for (i, be) in conditioned_space_adjacent_building_elements
            .iter()
            .enumerate()
        {
            assert_eq!(be.h_pli(), &results[i], "incorrect h_pli list returned");
        }
    }

    #[rstest]
    fn test_k_pli_for_conditioned_space_adjacent_elements(
        conditioned_space_adjacent_building_elements: [BuildingElementAdjacentConditionedSpace; 5],
    ) {
        let results = [
            [0.0, 0.0, 0.0, 0.0, 19000.0],
            [18000.0, 0.0, 0.0, 0.0, 0.0],
            [8500.0, 0.0, 0.0, 0.0, 8500.0],
            [2000.0, 4000.0, 4000.0, 4000.0, 2000.0],
            [0.0, 0.0, 15000.0, 0.0, 0.0],
        ];
        for (i, be) in conditioned_space_adjacent_building_elements
            .iter()
            .enumerate()
        {
            assert_eq!(be.k_pli(), &results[i], "incorrect k_pli list returned");
        }
    }

    #[rstest]
    fn test_fabric_heat_loss_for_conditioned_space_adjacent_elements(
        conditioned_space_adjacent_building_elements: [BuildingElementAdjacentConditionedSpace; 5],
    ) {
        for be in conditioned_space_adjacent_building_elements {
            assert_eq!(
                be.fabric_heat_loss(),
                0.0,
                "incorrect fabric heat loss returned"
            );
        }
    }

    #[rstest]
    fn test_heat_capacity_for_conditioned_space_adjacent_elements(
        conditioned_space_adjacent_building_elements: [BuildingElementAdjacentConditionedSpace; 5],
    ) {
        let results = [380., 405., 425., 440., 450.];

        for (i, be) in conditioned_space_adjacent_building_elements
            .iter()
            .enumerate()
        {
            assert_eq!(
                be.heat_capacity(),
                results[i],
                "incorrect heat capacity returned"
            );
        }
    }

    #[fixture]
    fn ground_building_elements(
        external_conditions_for_ground: Arc<ExternalConditions>,
    ) -> [BuildingElementGround; 5] {
        let be_i_floor_data = FloorData::SuspendedFloor {
            height_upper_surface: 0.5,
            thermal_transmission_walls: 0.5,
            area_per_perimeter_vent: 0.01,
            shield_fact_location: WindShieldLocation::Sheltered,
            thermal_resistance_of_insulation: 7.,
            edge_insulation: Default::default(),
        };
        let be_i = BuildingElementGround::new(
            20.0,
            20.0,
            180.,
            1.5,
            0.1,
            19000.0,
            MassDistributionClass::I,
            &be_i_floor_data,
            0.3,
            18.0,
            0.5,
            external_conditions_for_ground.clone(),
        )
        .unwrap();
        let be_e_floor_data = FloorData::SlabNoEdgeInsulation;
        let be_e = BuildingElementGround::new(
            22.5,
            22.5,
            135.,
            1.4,
            0.2,
            18000.0,
            MassDistributionClass::E,
            &be_e_floor_data,
            0.3,
            19.0,
            0.6,
            external_conditions_for_ground.clone(),
        )
        .unwrap();
        let edge_insulation_ie = vec![
            EdgeInsulation::Horizontal {
                width: 3.0,
                edge_thermal_resistance: 2.0,
            },
            EdgeInsulation::Vertical {
                depth: 1.0,
                edge_thermal_resistance: 2.0,
            },
        ];
        let be_ie_floor_data = FloorData::SlabEdgeInsulation {
            edge_insulation: edge_insulation_ie,
        };
        let be_ie = BuildingElementGround::new(
            25.0,
            25.0,
            90.,
            1.33,
            0.2,
            17000.0,
            MassDistributionClass::IE,
            &be_ie_floor_data,
            0.3,
            20.0,
            0.7,
            external_conditions_for_ground.clone(),
        )
        .unwrap();
        let be_d_floor_data = FloorData::HeatedBasement {
            depth_basement_floor: 2.3,
            thermal_resistance_of_basement_walls: 6.,
            edge_insulation: Default::default(),
        };
        let be_d = BuildingElementGround::new(
            27.5,
            27.5,
            45.,
            1.25,
            0.2,
            16000.0,
            MassDistributionClass::D,
            &be_d_floor_data,
            0.3,
            21.0,
            0.8,
            external_conditions_for_ground.clone(),
        )
        .unwrap();
        let be_m_floor_data = FloorData::UnheatedBasement {
            thermal_transmittance_of_floor_above_basement: 1.2,
            thermal_transmission_walls: 0.5,
            depth_basement_floor: 2.3,
            height_basement_walls: 2.3,
            thermal_resistance_of_basement_walls: 0.5,
            edge_insulation: Default::default(),
        };
        let be_m = BuildingElementGround::new(
            30.0,
            30.0,
            0.,
            1.0,
            0.3,
            15000.0,
            MassDistributionClass::M,
            &be_m_floor_data,
            0.3,
            22.0,
            0.9,
            external_conditions_for_ground,
        )
        .unwrap();
        [be_i, be_e, be_ie, be_d, be_m]
    }

    #[fixture]
    fn simulation_time_for_ground() -> SimulationTime {
        SimulationTime::new(742., 746., 1.)
    }

    #[fixture]
    fn external_conditions_for_ground(
        simulation_time_for_ground: SimulationTime,
    ) -> Arc<ExternalConditions> {
        let air_temp_day_jan = vec![
            0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 7.5, 10.0, 12.5, 15.0, 19.5,
            17.0, 15.0, 12.0, 10.0, 7.0, 5.0, 3.0, 1.0,
        ];
        let air_temp_day_feb: Vec<f64> = air_temp_day_jan.iter().map(|x| x + 1.0).collect();
        let air_temp_day_mar: Vec<f64> = air_temp_day_jan.iter().map(|x| x + 2.0).collect();
        let air_temp_day_apr: Vec<f64> = air_temp_day_jan.iter().map(|x| x + 3.0).collect();
        let air_temp_day_may: Vec<f64> = air_temp_day_jan.iter().map(|x| x + 4.0).collect();
        let air_temp_day_jun: Vec<f64> = air_temp_day_jan.iter().map(|x| x + 5.0).collect();
        let air_temp_day_jul: Vec<f64> = air_temp_day_jan.iter().map(|x| x + 6.0).collect();
        let air_temp_day_aug: Vec<f64> = air_temp_day_jan.iter().map(|x| x + 6.0).collect();
        let air_temp_day_sep: Vec<f64> = air_temp_day_jan.iter().map(|x| x + 5.0).collect();
        let air_temp_day_oct: Vec<f64> = air_temp_day_jan.iter().map(|x| x + 4.0).collect();
        let air_temp_day_nov: Vec<f64> = air_temp_day_jan.iter().map(|x| x + 3.0).collect();
        let air_temp_day_dec: Vec<f64> = air_temp_day_jan.iter().map(|x| x + 2.0).collect();

        let mut airtemp = vec![];
        for _ in 0..31 {
            airtemp.extend(&air_temp_day_jan);
        }
        for _ in 0..28 {
            airtemp.extend(&air_temp_day_feb);
        }
        for _ in 0..31 {
            airtemp.extend(&air_temp_day_mar);
        }
        for _ in 0..30 {
            airtemp.extend(&air_temp_day_apr);
        }
        for _ in 0..31 {
            airtemp.extend(&air_temp_day_may);
        }
        for _ in 0..30 {
            airtemp.extend(&air_temp_day_jun);
        }
        for _ in 0..31 {
            airtemp.extend(&air_temp_day_jul);
        }
        for _ in 0..31 {
            airtemp.extend(&air_temp_day_aug);
        }
        for _ in 0..30 {
            airtemp.extend(&air_temp_day_sep);
        }
        for _ in 0..31 {
            airtemp.extend(&air_temp_day_oct);
        }
        for _ in 0..30 {
            airtemp.extend(&air_temp_day_nov);
        }
        for _ in 0..31 {
            airtemp.extend(&air_temp_day_dec);
        }

        Arc::new(ExternalConditions::new(
            &simulation_time_for_ground.iter(),
            airtemp,
            vec![2.; 8760],
            vec![180.; 8760].into_iter().map(Into::into).collect(),
            vec![0.0; 8760],
            vec![0.0; 8760],
            vec![0.0; 8760],
            55.0,
            0.0,
            0,
            0,
            None,
            1.,
            None,
            Some(DaylightSavingsConfig::NotApplicable),
            false,
            false,
            None,
        ))
    }

    // skip test_init_invalid_floor_type and test_init_invalid_edge_type because Rust would not allow invalid values

    fn create_suspended_floor(wind_shield_location: WindShieldLocation) -> FloorData {
        FloorData::SuspendedFloor {
            height_upper_surface: 1.,
            thermal_transmission_walls: 0.5,
            area_per_perimeter_vent: 1.,
            shield_fact_location: wind_shield_location,
            thermal_resistance_of_insulation: 1.,
            edge_insulation: Default::default(),
        }
    }

    fn create_building_element_ground(
        simulation_time: SimulationTime,
        floor_data: FloorData,
    ) -> BuildingElementGround {
        let external_conditions = Arc::new(ExternalConditions::new(
            &simulation_time.iter(),
            vec![1.; 8760],
            vec![1.; 8760],
            vec![],
            vec![0.0; 8760],
            vec![0.0; 8760],
            vec![],
            55.0,
            0.0,
            0,
            0,
            None,
            1.,
            None,
            Some(DaylightSavingsConfig::NotApplicable),
            false,
            false,
            None,
        ));

        BuildingElementGround::new(
            30.0,
            30.0,
            0.,
            1.0,
            0.3,
            15000.0,
            MassDistributionClass::M,
            &floor_data,
            0.3,
            22.0,
            0.9,
            external_conditions,
        )
        .unwrap()
    }

    #[rstest]
    /// Test that the shield_fact_location affects the wind_s_factor from h_pi and h_pe values
    fn test_wind_shield_fact(simulation_time_for_ground: SimulationTime) {
        let suspended_floor_sheltered = create_suspended_floor(WindShieldLocation::Sheltered);
        let suspended_floor_average = create_suspended_floor(WindShieldLocation::Average);
        let suspended_floor_exposed = create_suspended_floor(WindShieldLocation::Exposed);

        let be_sheltered =
            create_building_element_ground(simulation_time_for_ground, suspended_floor_sheltered);
        assert_relative_eq!(be_sheltered.h_pi, 21.76809338521401);
        assert_relative_eq!(be_sheltered.h_pe, 5208.393950400433);

        let be_average =
            create_building_element_ground(simulation_time_for_ground, suspended_floor_average);
        assert_relative_eq!(be_average.h_pi, 20.285704885671986);
        assert_relative_eq!(be_average.h_pe, 26144.01441657715);

        let be_exposed =
            create_building_element_ground(simulation_time_for_ground, suspended_floor_exposed);
        assert_relative_eq!(be_exposed.h_pi, 19.75335084679448);
        assert_relative_eq!(be_exposed.h_pe, 96370.39574909392);
    }

    // skipping test_wind_shield_fact_invalid_location as Rust would not allow invalid shield factor

    #[rstest]
    fn test_no_of_nodes_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        for be in ground_building_elements {
            assert_eq!(be.number_of_nodes(), 5, "incorrect number of nodes");
            assert_eq!(
                be.number_of_inside_nodes(),
                3,
                "incorrect number of inside nodes"
            );
        }
    }

    #[rstest]
    fn test_area_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        // Define increment between test cases
        let area_inc = 2.5;

        for (i, be) in ground_building_elements.iter().enumerate() {
            assert_eq!(be.area(), 20.0 + i as f64 * area_inc,);
        }
    }

    #[rstest]
    fn test_heat_flow_direction_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        let temp_int_air = 20.0;
        let temp_int_surface = [19.0, 21.0, 22.0, 21.0, 19.0];
        let results = [
            HeatFlowDirection::Downwards,
            HeatFlowDirection::Upwards,
            HeatFlowDirection::Horizontal,
            HeatFlowDirection::Downwards,
            HeatFlowDirection::Upwards,
        ];

        for (i, be) in ground_building_elements.iter().enumerate() {
            assert_eq!(
                be.heat_flow_direction(temp_int_air, temp_int_surface[i]),
                results[i],
                "incorrect heat flow direction returned"
            );
        }
    }

    #[rstest]
    fn test_r_si_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        let results = [0.17, 0.17, 0.13, 0.10, 0.10];

        for (i, be) in ground_building_elements.iter().enumerate() {
            assert_relative_eq!(be.r_si(), results[i], max_relative = 0.05);
        }
    }

    #[rstest]
    fn test_h_ci_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        let temp_int_air = 20.0;
        let temp_int_surface = [19.0, 21.0, 22.0, 21.0, 19.0];
        let results = [0.7, 5.0, 2.5, 0.7, 5.0];

        for (i, be) in ground_building_elements.iter().enumerate() {
            assert_eq!(be.h_ci(temp_int_air, temp_int_surface[i]), results[i],);
        }
    }

    #[rstest]
    fn test_h_ri_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        for be in ground_building_elements.iter() {
            assert_eq!(be.h_ri(), 5.13, "incorrect h_ri returned");
        }
    }

    #[rstest]
    fn test_h_ce_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        let results = [
            15.78947368,
            91.30434783,
            20.59886422,
            10.34482759,
            5.084745763,
        ];

        for (i, be) in ground_building_elements.iter().enumerate() {
            assert_relative_eq!(be.h_ce(), results[i], max_relative = 1e-6);
        }
    }

    #[rstest]
    fn test_h_re_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        for be in ground_building_elements.iter() {
            assert_eq!(be.h_re(), 0.0, "incorrect h_re returned");
        }
    }

    #[rstest]
    fn test_solar_absorption_coeff_for_ground(
        ground_building_elements: [BuildingElementGround; 5],
    ) {
        for be in ground_building_elements.iter() {
            assert_eq!(
                be.solar_absorption_coeff(),
                0.0,
                "incorrect solar_absorption_coeff returned"
            );
        }
    }

    #[rstest]
    fn test_therm_rad_to_sky_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        for be in ground_building_elements.iter() {
            assert_eq!(
                be.therm_rad_to_sky(),
                0.0,
                "incorrect therm_rad_to_sky returned"
            );
        }
    }

    #[rstest]
    fn test_h_pli_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        let results = [
            [6.0, 5.217391304347826, 20.0, 40.0],
            [6.0, 4.615384615384615, 10.0, 20.0],
            [6.0, 4.615384615384615, 10.0, 20.0],
            [6.0, 4.615384615384615, 10.0, 20.0],
            [
                6.0,
                4.137931034482759,
                6.666666666666667,
                13.333333333333334,
            ],
        ];
        for (i, be) in ground_building_elements.iter().enumerate() {
            assert_eq!(be.h_pli(), &results[i], "incorrect h_pli list returned");
        }
    }

    #[rstest]
    fn test_k_pli_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        let results = [
            [0.0, 1500000.0, 0.0, 0.0, 19000.0],
            [0.0, 1500000.0, 18000.0, 0.0, 0.0],
            [0.0, 1500000.0, 8500.0, 0.0, 8500.0],
            [0.0, 1500000.0, 4000.0, 8000.0, 4000.0],
            [0.0, 1500000.0, 0.0, 15000.0, 0.0],
        ];
        for (i, be) in ground_building_elements.iter().enumerate() {
            assert_eq!(be.k_pli(), &results[i], "incorrect k_pli list returned");
        }
    }

    #[rstest]
    fn test_temp_ext_for_ground(
        ground_building_elements: [BuildingElementGround; 5],
        simulation_time_for_ground: SimulationTime,
    ) {
        let results = [
            [
                -0.6471789638993641,
                -0.6471789638993641,
                2.505131315506123,
                2.505131315506123,
            ],
            [
                7.428039862980361,
                7.428039862980361,
                8.234778286483786,
                8.234778286483786,
            ],
            [
                7.732888552541917,
                7.732888552541917,
                8.448604706949336,
                8.448604706949336,
            ],
            [
                8.366361777378224,
                8.366361777378224,
                8.86671506616569,
                8.86671506616569,
            ],
            [
                6.293446005722892,
                6.293446005722892,
                7.413004622032444,
                7.413004622032444,
            ],
        ];

        for (t_idx, t_it) in simulation_time_for_ground.iter().enumerate() {
            for (i, be) in ground_building_elements.iter().enumerate() {
                assert_eq!(be.temp_ext(t_it), results[i][t_idx],);
            }
        }
    }

    #[rstest]
    fn test_fabric_heat_loss_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        let expected = [30.0, 31.5, 33.25, 34.375, 30.0];
        for (i, be) in ground_building_elements.iter().enumerate() {
            assert_relative_eq!(be.fabric_heat_loss(), expected[i], max_relative = 1e-2);
        }
    }

    #[rstest]
    fn test_heat_capacity_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        let results = [380., 405., 425., 440., 450.];
        for (i, be) in ground_building_elements.iter().enumerate() {
            assert_eq!(
                be.heat_capacity(),
                results[i],
                "incorrect heat capacity returned"
            );
        }
    }

    // below window treatment tests are now in test_enums.py in the Python
    #[rstest]
    #[case(WindowTreatmentControl::Manual)]
    #[case(WindowTreatmentControl::ManualMotorised)]
    #[case(WindowTreatmentControl::AutoMotorised)]
    #[case(WindowTreatmentControl::CombinedLightBlindHvac)]
    fn test_is_either_automatic_or_manual(#[case] control: WindowTreatmentControl) {
        assert!(control.is_automatic() || control.is_manual());
    }

    #[rstest]
    #[case(WindowTreatmentControl::AutoMotorised)]
    #[case(WindowTreatmentControl::CombinedLightBlindHvac)]
    fn test_is_not_manual(#[case] control: WindowTreatmentControl) {
        assert!(!control.is_manual());
    }

    #[rstest]
    #[case(WindowTreatmentControl::Manual)]
    #[case(WindowTreatmentControl::ManualMotorised)]
    fn test_is_not_automatic(#[case] control: WindowTreatmentControl) {
        assert!(!control.is_automatic());
    }

    struct MockSolarRadiationInteraction(Option<f64>);
    impl SolarRadiationInteraction for MockSolarRadiationInteraction {
        fn base_height(&self) -> f64 {
            unreachable!()
        }
        fn external_pitch(&self) -> f64 {
            unreachable!()
        }
        fn orientation(&self) -> Option<Orientation360> {
            unreachable!()
        }
        fn projected_height(&self) -> f64 {
            unreachable!()
        }
        fn set_base_height(&mut self, _: f64) {
            unreachable!()
        }
        fn set_external_pitch(&mut self, _: f64) {
            unreachable!()
        }
        fn set_orientation(&mut self, _: Orientation360) {
            unreachable!()
        }
        fn set_projected_height(&mut self, _: f64) {
            unreachable!()
        }
        fn set_shading(&mut self, _shading: Option<Vec<WindowShadingObject>>) {
            unreachable!()
        }
        fn set_solar_absorption_coeff(&mut self, _: f64) {
            unreachable!()
        }
        fn set_width(&mut self, _: f64) {
            unreachable!()
        }
        fn shading(&self) -> &[WindowShadingObject] {
            unreachable!()
        }
        fn solar_absorption_coeff(&self) -> f64 {
            unreachable!()
        }
        fn width(&self) -> f64 {
            unreachable!()
        }
    }

    #[rstest]
    fn test_i_sol_dir_dif(mut simulation_time: SimulationTimeIterator) {
        assert_eq!(
            MockSolarRadiationInteraction(None).i_sol_dir_dif(simulation_time.next().unwrap()),
            (0.0, 0.0)
        );
    }

    #[rstest]
    fn test_solar_gains() {
        assert_eq!(
            SolarRadiationInteraction::solar_gains(&MockSolarRadiationInteraction(None)),
            0.0
        );
    }

    #[rstest]
    fn test_shading_factors_direct_diffuse(mut simulation_time: SimulationTimeIterator) {
        assert_eq!(
            SolarRadiationInteraction::shading_factors_direct_diffuse(
                &MockSolarRadiationInteraction(None),
                simulation_time.next().unwrap()
            ),
            (1.0, 1.0)
        );
    }

    // skip test_orientation_none as orientation is not optional

    impl SolarRadiationInteractionTransmitted for MockSolarRadiationInteraction {
        fn unconverted_g_value(&self) -> f64 {
            self.0.unwrap()
        }

        fn external_conditions(&self) -> &ExternalConditions {
            unreachable!()
        }

        fn pitch(&self) -> f64 {
            unreachable!()
        }

        fn area(&self) -> f64 {
            unreachable!()
        }

        fn frame_area_fraction(&self) -> f64 {
            unreachable!()
        }
    }

    #[rstest]
    fn test_convert_g_value() {
        assert_eq!(
            MockSolarRadiationInteraction(Some(0.5)).convert_g_value(),
            0.45
        );
    }

    #[rstest]
    fn test_solar_gains_for_transmitted(
        simulation_time: SimulationTimeIterator,
        mut transparent_building_element: BuildingElementTransparent,
    ) {
        let external_conditions = external_conditions_surface_irradiance(
            simulation_time.clone(),
            vec![0.5331371, 0., 0., 0.], // surface irradiance 0.5000000601251663
        );
        transparent_building_element.g_value = 0.5;
        transparent_building_element.pitch = 20.;
        transparent_building_element.area = 5.;
        transparent_building_element.frame_area_fraction = 0.2;
        transparent_building_element.external_conditions = external_conditions;

        let solar_gains = SolarRadiationInteractionTransmitted::solar_gains(
            &transparent_building_element,
            simulation_time.current_iteration(),
        )
        .unwrap();

        assert_relative_eq!(solar_gains, 0.9, max_relative = 1e-6);
    }

    // skip test_solar_gains_error as orientation and shading cannot be None at the point of the methods calls in this test

    #[fixture]
    fn transparent_building_element(
        simulation_time: SimulationTimeIterator,
    ) -> BuildingElementTransparent {
        let external_conditions = Arc::new(ExternalConditions::new(
            &simulation_time,
            vec![0.0, 5.0, 10.0, 15.0],
            vec![],
            vec![],
            vec![0.0; 4],
            vec![0.0; 4],
            vec![0.0; 4],
            55.0,
            0.0,
            0,
            0,
            None,
            1.0,
            None,
            None,
            false,
            false,
            None,
        ));

        BuildingElementTransparent::new(
            90.,
            0.4,
            Orientation360::create_from_180(180.).unwrap().into(),
            0.75,
            0.25,
            1.,
            1.25,
            4.,
            Default::default(),
            Default::default(),
            external_conditions,
        )
    }

    #[rstest]
    fn test_no_of_nodes_for_transparent(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.number_of_nodes(),
            2,
            "incorrect number of nodes"
        );
        assert_eq!(
            transparent_building_element.number_of_inside_nodes(),
            0,
            "incorrect number of inside nodes"
        );
    }

    #[rstest]
    fn test_area_for_transparent(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.area(),
            5.0,
            "incorrect area returned"
        );
    }

    #[rstest]
    fn test_heat_flow_direction(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.heat_flow_direction(10., 10.),
            HeatFlowDirection::Horizontal,
            "incorrect heat flow direction returned"
        );
    }

    #[rstest]
    fn test_r_si_for_transparent(transparent_building_element: BuildingElementTransparent) {
        assert_relative_eq!(
            transparent_building_element.r_si(),
            0.13,
            max_relative = 1e-2
        );
    }

    #[rstest]
    fn test_h_ci_for_transparent(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.h_ci(10., 10.),
            2.5,
            "incorrect h_ci returned"
        );
    }

    #[rstest]
    fn test_h_ri_for_transparent(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.h_ri(),
            5.13,
            "incorrect h_ri returned"
        );
    }

    #[rstest]
    fn test_h_ce_for_transparent(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.h_ce(),
            20.0,
            "incorrect h_ce returned"
        );
    }

    #[rstest]
    fn test_h_re_for_transparent(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.h_re(),
            4.14,
            "incorrect h_re returned"
        );
    }

    #[rstest]
    fn test_solar_absorption_coeff_for_transparent(
        transparent_building_element: BuildingElementTransparent,
    ) {
        assert_eq!(
            transparent_building_element.solar_absorption_coeff(),
            0.0,
            "incorrect solar_absorption_coeff returned"
        );
    }

    #[rstest]
    fn test_therm_rad_to_sky_for_transparent(
        transparent_building_element: BuildingElementTransparent,
    ) {
        assert_eq!(
            transparent_building_element.therm_rad_to_sky(),
            22.77,
            "incorrect therm_rad_to_sky returned"
        );
    }

    #[rstest]
    fn test_h_pli_for_transparent(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.h_pli(),
            &[2.5],
            "incorrect h_pli list returned"
        );
    }

    #[rstest]
    fn test_h_pli_with_treatment_for_transparent(
        simulation_time: SimulationTimeIterator,
        mut transparent_building_element: BuildingElementTransparent,
    ) {
        let window_treatment = WindowTreatment {
            _treatment_type: WindowTreatmentType::Curtains,
            controls: WindowTreatmentControl::Manual,
            delta_r: 0.2,
            trans_red: 0.0,
            closing_irradiance_control: None,
            opening_irradiance_control: None,
            open_control: None,
            is_open: Default::default(),
            opening_delay_hrs: 0.0,
            time_last_adjusted: Default::default(),
        };
        transparent_building_element.treatment = vec![window_treatment];

        assert_eq!(
            transparent_building_element
                .h_pli_by_index(0, simulation_time.current_iteration())
                .unwrap(),
            1.6666666666666665,
            "incorrect h_pli returned"
        );
    }

    #[rstest]
    fn test_k_pli_for_transparent(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.k_pli(),
            &[0.0, 0.0],
            "non-zero k_pli list returned"
        );
    }

    #[rstest]
    fn test_temp_ext_for_transparent(
        transparent_building_element: BuildingElementTransparent,
        simulation_time: SimulationTimeIterator,
    ) {
        for (t_idx, t_it) in simulation_time.enumerate() {
            assert_eq!(
                transparent_building_element.temp_ext(t_it),
                t_idx as f64 * 5.0,
                "incorrect temp ext returned"
            );
        }
    }

    // Can be used to return a specific surface_irradiance
    fn external_conditions_surface_irradiance(
        simulation_time: SimulationTimeIterator,
        diffuse_horizontal_radiations: Vec<f64>,
    ) -> Arc<ExternalConditions> {
        let shading_segments = vec![
            ShadingSegment {
                start360: Orientation360::create_from_180(180.).unwrap(),
                end360: Orientation360::create_from_180(135.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(135.).unwrap(),
                end360: Orientation360::create_from_180(90.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(90.).unwrap(),
                end360: Orientation360::create_from_180(45.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(45.).unwrap(),
                end360: Orientation360::create_from_180(0.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(0.).unwrap(),
                end360: Orientation360::create_from_180(-45.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(-45.).unwrap(),
                end360: Orientation360::create_from_180(-90.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(-90.).unwrap(),
                end360: Orientation360::create_from_180(-135.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(-135.).unwrap(),
                end360: Orientation360::create_from_180(-180.).unwrap(),
                ..Default::default()
            },
        ]
        .into();

        Arc::new(ExternalConditions::new(
            &simulation_time,
            vec![0.0, 5.0, 10.0, 15.0],
            vec![0.0; 4],
            vec![0.0; 4].into_iter().map(Into::into).collect(),
            diffuse_horizontal_radiations,
            vec![0.0, 0., 0., 0.],
            vec![0.0; 4],
            55.0,
            0.0,
            0,
            0,
            None,
            1.0,
            None,
            None,
            false,
            false,
            shading_segments,
        ))
    }
    #[rstest]
    fn test_solar_gains_for_transparent(
        simulation_time: SimulationTimeIterator,
        mut transparent_building_element: BuildingElementTransparent,
    ) {
        transparent_building_element.external_conditions = external_conditions_surface_irradiance(
            simulation_time.clone(),
            vec![11.87873, 0., 0., 0.], // surface irradiance 9.999995798989683
        );

        assert_relative_eq!(
            transparent_building_element
                .solar_gains(simulation_time.current_iteration())
                .unwrap(),
            25.3125,
            max_relative = 1e-6,
        );
    }

    #[rstest]
    fn test_solar_gains_with_treatment_for_transparent(
        simulation_time: SimulationTimeIterator,
        mut transparent_building_element: BuildingElementTransparent,
    ) {
        let window_treatment = WindowTreatment {
            _treatment_type: WindowTreatmentType::Curtains,
            controls: WindowTreatmentControl::Manual,
            delta_r: 0.2,
            trans_red: 0.3,
            closing_irradiance_control: None,
            opening_irradiance_control: None,
            open_control: None,
            is_open: Default::default(),
            opening_delay_hrs: 0.0,
            time_last_adjusted: Default::default(),
        };
        transparent_building_element.treatment = vec![window_treatment];
        transparent_building_element.external_conditions = external_conditions_surface_irradiance(
            simulation_time.clone(),
            vec![11.87873, 0., 0., 0.], // surface irradiance 9.999995798989683
        );

        assert_relative_eq!(
            transparent_building_element
                .solar_gains(simulation_time.current_iteration())
                .unwrap(),
            17.71875,
            max_relative = 1e-6,
        );

        transparent_building_element.treatment[0].is_open = AtomicBool::from(true);

        assert_relative_eq!(
            transparent_building_element
                .solar_gains(simulation_time.current_iteration())
                .unwrap(),
            25.3125,
            max_relative = 1e-6,
        );
    }

    fn create_setpoint_time_control(setpnt: f64) -> Arc<Control> {
        Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(setpnt)], // causes control.setpnt() to return specified value
            0,
            1.0,
            Default::default(),
            Default::default(),
            1.0,
        )))
    }

    #[rstest]
    fn test_adjust_treatment_open(
        simulation_time: SimulationTimeIterator,
        mut transparent_building_element: BuildingElementTransparent,
    ) {
        // Test that adjust_treatment opens when control is on
        let control = Arc::new(Control::OnOffTime(OnOffTimeControl::new(
            vec![Some(true)], // control is on
            0,
            1.,
        )));
        let setpoint_time_control = create_setpoint_time_control(20.);
        let window_treatment = WindowTreatment {
            _treatment_type: WindowTreatmentType::Curtains,
            controls: WindowTreatmentControl::Manual,
            delta_r: 0.2,
            trans_red: 0.3,
            closing_irradiance_control: Some(setpoint_time_control.clone()),
            opening_irradiance_control: Some(setpoint_time_control),
            open_control: Some(control),
            is_open: Default::default(),
            opening_delay_hrs: 0.0,
            time_last_adjusted: Default::default(),
        };
        transparent_building_element.treatment = vec![window_treatment];
        transparent_building_element.external_conditions = external_conditions_surface_irradiance(
            simulation_time.clone(),
            vec![11.87873, 0., 0., 0.], // surface irradiance 9.999995798989683
        );

        transparent_building_element
            .adjust_treatment(simulation_time.current_iteration())
            .unwrap();

        assert!(transparent_building_element.treatment[0]
            .is_open
            .load(Ordering::SeqCst))
    }

    #[rstest]
    fn test_adjust_treatment_not_opened(
        simulation_time: SimulationTimeIterator,
        mut transparent_building_element: BuildingElementTransparent,
    ) {
        // Test that adjust_treatment doesn't open when control is off
        let control = Arc::new(Control::OnOffTime(OnOffTimeControl::new(
            vec![Some(false)], // control is off
            0,
            1.,
        )));
        let setpoint_time_control = create_setpoint_time_control(20.);
        let window_treatment = WindowTreatment {
            _treatment_type: WindowTreatmentType::Curtains,
            controls: WindowTreatmentControl::Manual,
            delta_r: 0.2,
            trans_red: 0.3,
            closing_irradiance_control: Some(setpoint_time_control.clone()),
            opening_irradiance_control: Some(setpoint_time_control),
            open_control: Some(control),
            is_open: Default::default(),
            opening_delay_hrs: 0.0,
            time_last_adjusted: Default::default(),
        };
        transparent_building_element.treatment = vec![window_treatment];
        transparent_building_element.external_conditions = external_conditions_surface_irradiance(
            simulation_time.clone(),
            vec![11.87873, 0., 0., 0.], // surface irradiance 9.999995798989683
        );

        transparent_building_element
            .adjust_treatment(simulation_time.current_iteration())
            .unwrap();

        assert!(!transparent_building_element.treatment[0]
            .is_open
            .load(Ordering::SeqCst))
    }

    #[rstest]
    fn test_adjust_treatment_close(
        simulation_time: SimulationTimeIterator,
        mut transparent_building_element: BuildingElementTransparent,
    ) {
        // Test that adjust_treatment closes when control is off
        let control = Arc::new(Control::OnOffTime(OnOffTimeControl::new(
            vec![Some(false)], // control is off
            0,
            1.,
        )));
        let setpoint_time_control = create_setpoint_time_control(20.);
        let window_treatment = WindowTreatment {
            _treatment_type: WindowTreatmentType::Curtains,
            controls: WindowTreatmentControl::Manual,
            delta_r: 0.2,
            trans_red: 0.3,
            closing_irradiance_control: Some(setpoint_time_control.clone()),
            opening_irradiance_control: Some(setpoint_time_control),
            open_control: Some(control),
            is_open: AtomicBool::from(true),
            opening_delay_hrs: 0.0,
            time_last_adjusted: Default::default(),
        };
        transparent_building_element.treatment = vec![window_treatment];
        transparent_building_element.external_conditions = external_conditions_surface_irradiance(
            simulation_time.clone(),
            vec![11.87873, 0., 0., 0.], // surface irradiance 9.999995798989683
        );

        transparent_building_element
            .adjust_treatment(simulation_time.current_iteration())
            .unwrap();

        assert!(!transparent_building_element.treatment[0]
            .is_open
            .load(Ordering::SeqCst))
    }

    #[rstest]
    fn test_adjust_treatment_not_closed(
        simulation_time: SimulationTimeIterator,
        mut transparent_building_element: BuildingElementTransparent,
    ) {
        let external_conditions = external_conditions_surface_irradiance(
            simulation_time.clone(),
            vec![19.77, 0., 0., 0.], // surface irradiance 29.99986997757254
        );
        // Test that adjust_treatment doesn't close when control is on
        let control = Arc::new(Control::OnOffTime(OnOffTimeControl::new(
            vec![Some(true)], // control is on
            0,
            1.,
        )));
        let setpoint_time_control = create_setpoint_time_control(20.);
        let window_treatment = WindowTreatment {
            _treatment_type: WindowTreatmentType::Curtains,
            controls: WindowTreatmentControl::Manual,
            delta_r: 0.2,
            trans_red: 0.3,
            closing_irradiance_control: Some(setpoint_time_control.clone()),
            opening_irradiance_control: Some(setpoint_time_control),
            open_control: Some(control),
            is_open: AtomicBool::from(true),
            opening_delay_hrs: 0.0,
            time_last_adjusted: Default::default(),
        };
        transparent_building_element.treatment = vec![window_treatment];
        transparent_building_element.external_conditions = external_conditions;

        transparent_building_element
            .adjust_treatment(simulation_time.current_iteration())
            .unwrap();

        assert!(transparent_building_element.treatment[0]
            .is_open
            .load(Ordering::SeqCst))
    }

    #[rstest]
    fn test_adjust_treatment_open_irrad(
        simulation_time: SimulationTimeIterator,
        mut transparent_building_element: BuildingElementTransparent,
    ) {
        // Test that adjust_treatment opens from irradiance values
        let setpoint_time_control = create_setpoint_time_control(40.);
        let window_treatment = WindowTreatment {
            _treatment_type: WindowTreatmentType::Curtains,
            controls: WindowTreatmentControl::Manual,
            delta_r: 0.2,
            trans_red: 0.3,
            closing_irradiance_control: Some(setpoint_time_control.clone()),
            opening_irradiance_control: Some(setpoint_time_control),
            open_control: None,
            is_open: Default::default(),
            opening_delay_hrs: 0.0,
            time_last_adjusted: Default::default(),
        };

        transparent_building_element.treatment = vec![window_treatment];
        // The Python test sets the surface irradiance return value to 30, our building element gets
        // 29.99986997757254 for surface irradiance. The same conditions are triggered in both tests
        // because surface_irradiance < control.setpnt()
        transparent_building_element.external_conditions = external_conditions_surface_irradiance(
            simulation_time.clone(),
            vec![19.77, 0., 0., 0.], // surface irradiance 29.99986997757254
        );

        transparent_building_element
            .adjust_treatment(simulation_time.current_iteration())
            .unwrap();

        assert!(transparent_building_element.treatment[0]
            .is_open
            .load(Ordering::SeqCst));
    }

    #[rstest]
    fn test_adjust_treatment_close_irrad(
        simulation_time: SimulationTimeIterator,
        mut transparent_building_element: BuildingElementTransparent,
    ) {
        // Test that adjust_treatment closes from irradiance values
        let setpoint_time_control = create_setpoint_time_control(20.);
        let window_treatment = WindowTreatment {
            _treatment_type: WindowTreatmentType::Curtains,
            controls: WindowTreatmentControl::Manual,
            delta_r: 0.2,
            trans_red: 0.3,
            closing_irradiance_control: Some(setpoint_time_control.clone()),
            opening_irradiance_control: Some(setpoint_time_control),
            open_control: None,
            is_open: AtomicBool::from(true),
            opening_delay_hrs: 0.0,
            time_last_adjusted: Default::default(),
        };
        transparent_building_element.treatment = vec![window_treatment];
        // The Python test sets the surface irradiance return value to 30, our building element gets
        // 29.99986997757254 for surface irradiance. The same conditions are triggered in both tests
        // because surface_irradiance > control.setpnt()
        transparent_building_element.external_conditions = external_conditions_surface_irradiance(
            simulation_time.clone(),
            vec![19.77, 0., 0., 0.], // surface irradiance 29.99986997757254
        );

        transparent_building_element
            .adjust_treatment(simulation_time.current_iteration())
            .unwrap();

        assert!(!transparent_building_element.treatment[0]
            .is_open
            .load(Ordering::SeqCst))
    }

    // skip test_adjust_treatment_orientation_shading_error as orientation and shading cannot be None at the point of the methods calls in this test

    #[rstest]
    fn test_fabric_heat_loss_for_transparent(
        transparent_building_element: BuildingElementTransparent,
    ) {
        assert_relative_eq!(
            transparent_building_element.fabric_heat_loss(),
            8.16,
            max_relative = 1e-2
        );
    }

    #[rstest]
    fn test_heat_capacity(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.heat_capacity(),
            0.,
            "incorrect heat capacity returned"
        );
    }

    #[rstest]
    fn test_projected_height(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.projected_height(),
            1.25,
            "incorrect projected height returned"
        );
    }

    #[rstest]
    fn test_orientation(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element
                .orientation()
                .unwrap()
                .transform_to_180(),
            180.,
            "incorrect orientation returned"
        );
    }

    #[fixture]
    fn unconditioned_space_adjacent_building_elements(
        external_conditions: Arc<ExternalConditions>,
    ) -> [BuildingElementAdjacentUnconditionedSpaceSimple; 5] {
        let be_i = BuildingElementAdjacentUnconditionedSpaceSimple::new(
            20.0,
            180.,
            0.25,
            0.5,
            19000.0,
            MassDistributionClass::I,
            external_conditions.clone(),
        );
        let be_e = BuildingElementAdjacentUnconditionedSpaceSimple::new(
            22.5,
            135.,
            0.50,
            1.,
            18000.0,
            MassDistributionClass::E,
            external_conditions.clone(),
        );
        let be_ie = BuildingElementAdjacentUnconditionedSpaceSimple::new(
            25.0,
            90.,
            0.75,
            1.5,
            17000.0,
            MassDistributionClass::IE,
            external_conditions.clone(),
        );
        let be_d = BuildingElementAdjacentUnconditionedSpaceSimple::new(
            27.5,
            45.,
            0.80,
            2.,
            16000.0,
            MassDistributionClass::D,
            external_conditions.clone(),
        );
        let be_m = BuildingElementAdjacentUnconditionedSpaceSimple::new(
            30.0,
            0.,
            0.40,
            2.5,
            15000.0,
            MassDistributionClass::M,
            external_conditions.clone(),
        );
        [be_i, be_e, be_ie, be_d, be_m]
    }

    #[rstest]
    fn test_h_ce_for_unconditioned_space_adjacent_elements(
        unconditioned_space_adjacent_building_elements: [BuildingElementAdjacentUnconditionedSpaceSimple;
            5],
    ) {
        for (i, be) in unconditioned_space_adjacent_building_elements
            .iter()
            .enumerate()
        {
            assert_relative_eq!(
                be.h_ce(),
                [
                    1.8469778117827087,
                    0.960222752585521,
                    0.6487503359312012,
                    0.4898538961038961,
                    0.39348003259983705,
                ][i],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    fn test_h_re_for_for_unconditioned_space_adjacent_elements(
        unconditioned_space_adjacent_building_elements: [BuildingElementAdjacentUnconditionedSpaceSimple;
            5],
    ) {
        for (i, be) in unconditioned_space_adjacent_building_elements
            .iter()
            .enumerate()
        {
            assert_relative_eq!(
                be.h_re(),
                [0.0, 0.0, 0.0, 0.0, 0.0,][i],
                max_relative = 1e-8
            );
        }
    }

    #[rstest]
    fn test_temp_ext_for_for_unconditioned_space_adjacent_elements(
        unconditioned_space_adjacent_building_elements: [BuildingElementAdjacentUnconditionedSpaceSimple;
            5],
        simulation_time: SimulationTimeIterator,
    ) {
        let results = [
            [0.0, 5.0, 10.0, 15.0],
            [0.0, 5.0, 10.0, 15.0],
            [0.0, 5.0, 10.0, 15.0],
            [0.0, 5.0, 10.0, 15.0],
            [0.0, 5.0, 10.0, 15.0],
        ];
        for (t_idx, t_it) in simulation_time.enumerate() {
            for (i, be) in unconditioned_space_adjacent_building_elements
                .iter()
                .enumerate()
            {
                assert_eq!(be.temp_ext(t_it), results[i][t_idx]);
            }
        }
    }

    #[rstest]
    fn test_fabric_heat_loss_for_for_unconditioned_space_adjacent_elements(
        unconditioned_space_adjacent_building_elements: [BuildingElementAdjacentUnconditionedSpaceSimple;
            5],
    ) {
        for (i, be) in unconditioned_space_adjacent_building_elements
            .iter()
            .enumerate()
        {
            assert_relative_eq!(
                be.fabric_heat_loss(),
                [43.20, 31.56, 27.10, 29.25, 55.54][i],
                max_relative = 1e-4
            );
        }
    }

    #[rstest]
    fn test_heat_capacity_for_for_unconditioned_space_adjacent_elements(
        unconditioned_space_adjacent_building_elements: [BuildingElementAdjacentUnconditionedSpaceSimple;
            5],
    ) {
        let results = [380., 405., 425., 440., 450.];
        for (i, be) in unconditioned_space_adjacent_building_elements
            .iter()
            .enumerate()
        {
            assert_eq!(be.heat_capacity(), results[i]);
        }
    }
    mod test_building_element_party_wall {
        use super::*;
        use std::assert_eq;

        #[fixture]
        fn simtime() -> SimulationTime {
            SimulationTime::new(0., 4., 1.)
        }

        #[fixture]
        fn area() -> f64 {
            10.
        }

        #[fixture]
        fn pitch() -> f64 {
            90.
        }

        #[fixture]
        fn thermal_resistance_construction() -> f64 {
            0.5
        }

        #[fixture]
        fn areal_heat_capacity() -> f64 {
            10_000.
        }

        #[fixture]
        fn mass_distribution_class() -> MassDistributionClass {
            MassDistributionClass::D
        }

        #[rstest]
        fn test_calculate_cavity_resistance_defined_resistance_with_value(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let thermal_resistance_cavity = 3.0;
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::DefinedResistance,
                None,
                Some(thermal_resistance_cavity),
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            // Check that h_ce reflects the custom resistance value
            // h_ce = 1 / (R_se + R_cavity) where R_se = 1/(H_CE + H_RE) = 1/24.14 ≈ 0.0414
            // h_ce = 1 / (0.0414 + 3.0) = 1 / 3.0414 ≈ 0.3287932443475892
            assert_relative_eq!(party_wall.h_ce(), 0.3287932443475892);
        }

        // skipped test_calculate_cavity_resistance_invalid_type_raises_error as Rust typing takes of this
        // skipped test_calculate_cavity_resistance_invalid_lining_type_raises_error as Rust typing takes of this

        #[rstest]
        fn test_calculate_cavity_resistance_incompatible_lining_type_raises_error(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::Solid,
                Some(PartyWallLiningType::WetPlaster),
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            );
            assert!(party_wall.is_err());
            assert_eq!(
                party_wall.unwrap_err().to_string(),
                "invalid combination of party wall cavity type and party wall lining type"
            )
        }

        #[rstest]
        fn test_calculate_cavity_resistance_solid_type(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::Solid,
                None,
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            // For solid type, cavity resistance should be 999999, making h_ce effectively zero
            assert_eq!(party_wall.h_ce(), 0.);
            // Fabric heat loss should be zero (no heat loss through party wall)
            assert_relative_eq!(party_wall.fabric_heat_loss(), 0.);
        }

        #[rstest]
        fn test_calculate_cavity_resistance_unfilled_unsealed_dry_lined(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::UnfilledUnsealed,
                Some(PartyWallLiningType::DryLined),
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            // Check h_ce value for unsealed cavity (R_cavity = 1.2)
            // h_ce = 1 / (R_se + R_cavity) = 1 / (0.0414 + 1.2) ≈ 0.805542
            assert_relative_eq!(party_wall.h_ce(), 0.8055258942872398);

            // Check fabric heat loss calculation
            let expected_heat_loss = 5.340492100175494;
            assert_relative_eq!(party_wall.fabric_heat_loss(), expected_heat_loss);
        }

        #[rstest]
        fn test_calculate_cavity_resistance_unfilled_sealed_dry_lined(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::UnfilledSealed,
                Some(PartyWallLiningType::DryLined),
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            // Check h_ce value for sealed cavity (R_cavity = 4.5)
            assert_relative_eq!(party_wall.h_ce(), 0.2201952, max_relative = 1e-7);

            // Check fabric heat loss calculation
            let expected_heat_loss = 1.933306112766621;
            assert_relative_eq!(party_wall.fabric_heat_loss(), expected_heat_loss);
        }

        #[rstest]
        fn test_calculate_cavity_resistance_filled_unsealed_dry_lined(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::FilledUnsealed,
                Some(PartyWallLiningType::DryLined),
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            // Check h_ce value for sealed cavity (R_cavity = 4.5)
            assert_relative_eq!(party_wall.h_ce(), 0.2201952, max_relative = 1e-7);

            // Check fabric heat loss calculation
            let expected_heat_loss = 1.933306112766621;
            assert_relative_eq!(party_wall.fabric_heat_loss(), expected_heat_loss);
        }

        #[rstest]
        fn test_calculate_cavity_resistance_unfilled_unsealed_wet_plaster(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::UnfilledUnsealed,
                Some(PartyWallLiningType::WetPlaster),
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            // Check h_ce value for unsealed cavity (R_cavity = 4.5)
            assert_relative_eq!(party_wall.h_ce(), 0.2201952, max_relative = 1e-7);

            // Check fabric heat loss calculation
            let expected_heat_loss = 1.933306112766621;
            assert_relative_eq!(party_wall.fabric_heat_loss(), expected_heat_loss);
        }

        #[rstest]
        fn test_calculate_cavity_resistance_unfilled_sealed_wet_plaster(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::UnfilledSealed,
                Some(PartyWallLiningType::WetPlaster),
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            // Check h_ce value for sealed cavity (R_cavity = 4.5)
            assert_relative_eq!(party_wall.h_ce(), 0.2201952, max_relative = 1e-7);

            // Check fabric heat loss calculation
            let expected_heat_loss = 1.933306112766621;
            assert_relative_eq!(party_wall.fabric_heat_loss(), expected_heat_loss);
        }

        #[rstest]
        fn test_calculate_cavity_resistance_filled_unsealed_wet_plaster(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::FilledUnsealed,
                Some(PartyWallLiningType::WetPlaster),
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            // Check h_ce value for sealed cavity (R_cavity = 4.5)
            assert_relative_eq!(party_wall.h_ce(), 0.2201952, max_relative = 1e-7);

            // Check fabric heat loss calculation
            let expected_heat_loss = 1.933306112766621;
            assert_relative_eq!(party_wall.fabric_heat_loss(), expected_heat_loss);
        }

        #[rstest]
        fn test_calculate_cavity_resistance_filled_sealed(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::FilledSealed,
                None,
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            // For filled_sealed type, cavity resistance should be 999999, making h_ce zero
            assert_eq!(party_wall.h_ce(), 0.);

            // Fabric heat loss should be zero (no heat loss through party wall)
            assert_relative_eq!(party_wall.fabric_heat_loss(), 0.);
        }

        #[rstest]
        fn test_all_party_wall_cavity_types_valid(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            const VALID_TYPES: [PartyWallCavityType; 6] = [
                PartyWallCavityType::Solid,
                PartyWallCavityType::UnfilledUnsealed,
                PartyWallCavityType::UnfilledSealed,
                PartyWallCavityType::FilledSealed,
                PartyWallCavityType::FilledUnsealed,
                PartyWallCavityType::DefinedResistance,
            ];

            for cavity_type in VALID_TYPES {
                let thermal_resistance = match cavity_type {
                    PartyWallCavityType::DefinedResistance => Some(2.5),
                    _ => None,
                };
                let lining_type = match cavity_type {
                    PartyWallCavityType::UnfilledUnsealed
                    | PartyWallCavityType::FilledUnsealed
                    | PartyWallCavityType::UnfilledSealed => Some(PartyWallLiningType::DryLined),
                    _ => None,
                };

                let party_wall = BuildingElementPartyWall::new(
                    area,
                    pitch,
                    thermal_resistance_construction,
                    cavity_type,
                    lining_type,
                    thermal_resistance,
                    areal_heat_capacity,
                    mass_distribution_class,
                    external_conditions.clone(),
                )
                .unwrap();
                assert_relative_eq!(party_wall.area, area);
                assert_relative_eq!(party_wall.pitch, pitch);
            }
        }

        #[rstest]
        fn test_h_ce_returns_zero_for_solid_type(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::Solid,
                None,
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            assert_eq!(party_wall.h_ce(), 0.);
            // Also verify h_re is zero (no radiative transfer)
            assert_relative_eq!(party_wall.h_re(), 0.);
        }

        #[rstest]
        fn test_h_ce_returns_zero_for_filled_sealed_type(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::FilledSealed,
                None,
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            assert_eq!(party_wall.h_ce(), 0.);
            // Also verify h_re is zero (no radiative transfer)
            assert_relative_eq!(party_wall.h_re(), 0.);
        }

        #[rstest]
        fn test_h_ce_returns_parent_value_for_unfilled_unsealed_type(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::UnfilledUnsealed,
                Some(PartyWallLiningType::DryLined),
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            assert_eq!(party_wall.h_ce(), 0.8055258942872398);
            // skipping second assertion as duplicated
        }

        #[rstest]
        fn test_h_ce_returns_parent_value_for_unfilled_sealed_type(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::UnfilledSealed,
                Some(PartyWallLiningType::DryLined),
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            assert_eq!(party_wall.h_ce(), 0.22019520204323634);
            // skipping second assertion as duplicated
        }

        #[rstest]
        fn test_h_ce_returns_parent_value_for_defined_resistance_type(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let thermal_resistance_cavity = 2.;

            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::DefinedResistance,
                None,
                Some(thermal_resistance_cavity),
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            assert_eq!(party_wall.h_ce(), 0.4898538961038961);
            // skipping second assertion as duplicated
        }

        #[rstest]
        fn test_fabric_heat_loss_comparison_across_types(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall_unsealed = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::UnfilledUnsealed,
                Some(PartyWallLiningType::DryLined),
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions.clone(),
            )
            .unwrap();

            let party_wall_sealed = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::UnfilledSealed,
                Some(PartyWallLiningType::DryLined),
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions.clone(),
            )
            .unwrap();

            let party_wall_solid = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::Solid,
                None,
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            let heat_loss_unsealed = party_wall_unsealed.fabric_heat_loss();
            let heat_loss_sealed = party_wall_sealed.fabric_heat_loss();
            let heat_loss_solid = party_wall_solid.fabric_heat_loss();

            // Verify heat loss decreases as cavity resistance increases
            assert!(heat_loss_unsealed > heat_loss_sealed);
            assert!(heat_loss_sealed > heat_loss_solid);
            assert_relative_eq!(heat_loss_solid, 0.);
        }

        #[rstest]
        fn test_heat_capacity_calculation(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::UnfilledUnsealed,
                Some(PartyWallLiningType::DryLined),
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            // Heat capacity = area * (areal_heat_capacity / 1000)
            // Heat capacity = 10.0 * (10000 / 1000) = 100.0 kJ/K
            let expected_heat_capacity = area * (areal_heat_capacity / 1000.);
            assert_eq!(party_wall.heat_capacity(), expected_heat_capacity);
        }

        #[rstest]
        fn test_r_si_calculation(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::UnfilledUnsealed,
                Some(PartyWallLiningType::DryLined),
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            // For vertical walls (60 < pitch < 120), R_si = R_SI_HORIZONTAL
            // R_SI_HORIZONTAL = 1 / (H_RI + H_CI_HORIZONTAL) = 1 / (5.13 + 2.5) ≈ 0.131
            let expected_r_si = 1.0 / (5.13 + 2.5);
            assert_eq!(party_wall.r_si(), expected_r_si);
        }

        #[rstest]
        fn test_no_of_nodes(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::UnfilledUnsealed,
                Some(PartyWallLiningType::DryLined),
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            // 5-node model has 5 nodes total, 3 inside nodes
            assert_eq!(party_wall.number_of_nodes(), 5);
            assert_eq!(party_wall.number_of_inside_nodes(), 3);
        }

        #[rstest]
        fn test_node_conductances_h_pli(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            let party_wall = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::UnfilledUnsealed,
                Some(PartyWallLiningType::DryLined),
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            // For 5-node model:
            // h_pli[0] = 6 / R_c (outer)
            // h_pli[1] = 3 / R_c (inner)
            // h_pli[2] = 3 / R_c (inner)
            // h_pli[3] = 6 / R_c (outer)
            let h_outer = 6.0 / thermal_resistance_construction;
            let h_inner = 3.0 / thermal_resistance_construction;

            assert_relative_eq!(party_wall.h_pli()[0], h_outer);
            assert_relative_eq!(party_wall.h_pli()[1], h_inner);
            assert_relative_eq!(party_wall.h_pli()[2], h_inner);
            assert_relative_eq!(party_wall.h_pli()[3], h_outer);
        }

        #[rstest]
        fn test_equivalent_u_values_match_sap_guidance(
            area: f64,
            pitch: f64,
            thermal_resistance_construction: f64,
            areal_heat_capacity: f64,
            mass_distribution_class: MassDistributionClass,
            external_conditions: Arc<ExternalConditions>,
        ) {
            // Create party walls with standard construction assumptions
            // Assume R_c represents typical party wall construction

            // Test unsealed cavity: should give equivalent U ≈ 0.5 W/m²K
            let party_wall_unsealed = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::UnfilledUnsealed,
                Some(PartyWallLiningType::DryLined),
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions.clone(),
            )
            .unwrap();

            let u_unsealed = 0.5340492100175493;
            let heat_loss_unsealed = party_wall_unsealed.fabric_heat_loss();
            let calculated_u_unsealed = heat_loss_unsealed / area;

            // Check U-value is in reasonable range for unsealed cavity
            assert_relative_eq!(calculated_u_unsealed, u_unsealed);

            let party_wall_sealed = BuildingElementPartyWall::new(
                area,
                pitch,
                thermal_resistance_construction,
                PartyWallCavityType::UnfilledSealed,
                Some(PartyWallLiningType::DryLined),
                None,
                areal_heat_capacity,
                mass_distribution_class,
                external_conditions,
            )
            .unwrap();

            let u_sealed = 0.1933306112766621;
            let heat_loss_sealed = party_wall_sealed.fabric_heat_loss();
            let calculated_u_sealed = heat_loss_sealed / area;

            // Check U-value is in reasonable range for sealed cavity
            assert_relative_eq!(calculated_u_sealed, u_sealed);
        }
    }
}
