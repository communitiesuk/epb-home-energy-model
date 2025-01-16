use crate::core::units::{average_monthly_to_annual, JOULES_PER_KILOJOULE};
use crate::external_conditions::{ExternalConditions, WindowShadingObject};
use crate::input::{
    BuildingElement as BuildingElementInput, EdgeInsulation, FloorType, MassDistributionClass,
    WindShieldLocation, WindowTreatment,
};
use crate::simulation_time::{SimulationTimeIteration, SimulationTimeIterator};
use anyhow::{anyhow, bail};
use std::f64::consts::PI;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

// Difference between external air temperature and sky temperature
// (default value for intermediate climatic region from BS EN ISO 52016-1:2017, Table B.19)
const TEMP_DIFF_SKY: f64 = 11.0; // Kelvin

// Values from BS EN ISO 13789:2017, Table 8: Conventional surface heat
// transfer coefficients
const H_CI_UPWARDS: f64 = 5.0;
const H_CI_HORIZONTAL: f64 = 2.5;
const H_CI_DOWNWARDS: f64 = 0.7;
const H_CE: f64 = 20.0;
const H_RI: f64 = 5.13;
const H_RE: f64 = 4.14;

// Surface resistances of building elements, in m2 K / W
const R_SI_HORIZONTAL: f64 = 1.0 / (H_RI + H_CI_HORIZONTAL);
const R_SI_UPWARDS: f64 = 1.0 / (H_RI + H_CI_UPWARDS);
const R_SI_DOWNWARDS: f64 = 1.0 / (H_RI + H_CI_DOWNWARDS);
const R_SE: f64 = 1.0 / (H_CE + H_RE);

// From BR 443: The values under "horizontal" apply to heat flow
// directions +/- 30 degrees from horizontal plane.
pub const PITCH_LIMIT_HORIZ_CEILING: f64 = 60.0;
const PITCH_LIMIT_HORIZ_FLOOR: f64 = 120.0;

#[derive(Clone, Debug)]
pub enum BuildingElement {
    Opaque(BuildingElementOpaque),
    AdjacentZTC(BuildingElementAdjacentZTC),
    AdjacentZTUSimple(BuildingElementAdjacentZTUSimple),
    Ground(BuildingElementGround),
    Transparent(BuildingElementTransparent),
}

// macro so accessing individual building elements through the enum isn't so repetitive
macro_rules! per_element {
    ($val:expr, $pattern:pat => { $res:expr }) => {
        match $val {
            BuildingElement::Opaque($pattern) => $res,
            BuildingElement::AdjacentZTC($pattern) => $res,
            BuildingElement::AdjacentZTUSimple($pattern) => $res,
            BuildingElement::Ground($pattern) => $res,
            BuildingElement::Transparent($pattern) => $res,
        }
    };
}

pub trait BuildingElementBehaviour {
    fn area(&self) -> f64;

    fn a_sol(&self) -> f64;

    fn therm_rad_to_sky(&self) -> f64;

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
        let pitch = self.pitch();
        match pitch {
            PITCH_LIMIT_HORIZ_CEILING..=PITCH_LIMIT_HORIZ_FLOOR => R_SI_HORIZONTAL,
            ..PITCH_LIMIT_HORIZ_CEILING => R_SI_UPWARDS,
            PITCH_LIMIT_HORIZ_FLOOR.. => R_SI_DOWNWARDS,
            _ => unreachable!("Rust cannot tell that above is exhaustive"),
        }
    }

    /// Return external surface resistance, in m2 K / W
    fn r_se(&self) -> f64 {
        R_SE
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

    /// Return external convective heat transfer coefficient, in W / (m2.K)
    fn h_ce(&self) -> f64 {
        H_CE
    }

    /// Return external radiative heat transfer coefficient, in W / (m2.K)
    fn h_re(&self) -> f64 {
        H_RE
    }

    /// Return default of zero for i_sol_dir and i_sol_dif
    fn i_sol_dir_dif(&self, _simtime: SimulationTimeIteration) -> (f64, f64) {
        Default::default()
    }

    /// Return default of zero for solar gains
    fn solar_gains(&self, _simtime: SimulationTimeIteration) -> anyhow::Result<f64> {
        Ok(Default::default())
    }

    /// Return default of one for shading factor (no shading)
    fn shading_factors_direct_diffuse(
        &self,
        _simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        Ok((1.0, 1.0))
    }

    fn k_pli(&self) -> &[f64];

    fn h_pli(&self) -> &[f64];

    fn h_pli_by_index_unchecked(&self, idx: usize) -> f64 {
        self.h_pli()[idx]
    }

    /// Return number of nodes including external and internal layers
    fn number_of_nodes(&self) -> usize {
        self.k_pli().len()
    }

    /// Return number of nodes excluding external and internal layers
    fn number_of_inside_nodes(&self) -> usize {
        self.number_of_nodes() - 2
    }

    /// Return the temperature of the air on the other side of the building element
    fn temp_ext(&self, simtime: SimulationTimeIteration) -> f64;

    /// Return the fabric heat loss for the building element
    fn fabric_heat_loss(&self) -> f64;

    /// Return the fabric heat capacity for the building element
    fn heat_capacity(&self) -> f64;
}

pub(crate) use per_element;

/// Implement common interface on BuildingElement wrapper to delegate down to specialised building element struct types
impl BuildingElementBehaviour for BuildingElement {
    fn area(&self) -> f64 {
        per_element!(self, el => { el.area() })
    }

    fn a_sol(&self) -> f64 {
        per_element!(self, el => { el.a_sol() })
    }

    fn therm_rad_to_sky(&self) -> f64 {
        per_element!(self, el => { el.therm_rad_to_sky() })
    }

    fn heat_flow_direction(&self, temp_int_air: f64, temp_int_surface: f64) -> HeatFlowDirection {
        per_element!(self, el => { el.heat_flow_direction(temp_int_air, temp_int_surface) })
    }

    fn pitch(&self) -> f64 {
        per_element!(self, el => { el.pitch() })
    }

    fn r_si(&self) -> f64 {
        per_element!(self, el => { el.r_si() })
    }

    fn r_se(&self) -> f64 {
        per_element!(self, el => { el.r_se() })
    }

    fn h_ci(&self, temp_int_air: f64, temp_int_surface: f64) -> f64 {
        per_element!(self, el => { el.h_ci(temp_int_air, temp_int_surface) })
    }

    fn h_ri(&self) -> f64 {
        per_element!(self, el => { el.h_ri() })
    }

    fn h_ce(&self) -> f64 {
        per_element!(self, el => { el.h_ce() })
    }

    fn h_re(&self) -> f64 {
        per_element!(self, el => { el.h_re() })
    }

    fn i_sol_dir_dif(&self, simtime: SimulationTimeIteration) -> (f64, f64) {
        per_element!(self, el => { el.i_sol_dir_dif(simtime) })
    }

    fn solar_gains(&self, simtime: SimulationTimeIteration) -> anyhow::Result<f64> {
        per_element!(self, el => { el.solar_gains(simtime) })
    }

    fn shading_factors_direct_diffuse(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        per_element!(self, el => { el.shading_factors_direct_diffuse(simtime) })
    }

    fn k_pli(&self) -> &[f64] {
        per_element!(self, el => { el.k_pli() })
    }

    fn h_pli(&self) -> &[f64] {
        per_element!(self, el => { el.h_pli() })
    }

    fn temp_ext(&self, simtime: SimulationTimeIteration) -> f64 {
        per_element!(self, el => { el.temp_ext(simtime) })
    }

    fn fabric_heat_loss(&self) -> f64 {
        per_element!(self, el => { el.fabric_heat_loss() })
    }

    fn heat_capacity(&self) -> f64 {
        per_element!(self, el => { el.heat_capacity() })
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

fn init_therm_rad_to_sky(f_sky: f64) -> f64 {
    f_sky * H_RE * TEMP_DIFF_SKY
}

/// A struct to represent opaque building elements (walls, roofs, etc.)
#[derive(Clone, Debug)]
pub struct BuildingElementOpaque {
    base_height: f64,
    width: f64,
    projected_height: f64,
    orientation: f64,
    external_conditions: Arc<ExternalConditions>,
    area: f64,
    r_c: f64,
    k_m: f64,
    pub a_sol: f64,
    external_pitch: f64,
    pitch: f64,
    pub therm_rad_to_sky: f64,
    h_pli: [f64; 4],
    k_pli: [f64; 5],
}

/// Arguments (names based on those in BS EN ISO 52016-1:2017):
/// * `area` - net area of the opaque building element (i.e. minus any windows / doors / etc.)
/// * `is_unheated_pitched_roof`
/// * `pitch` - tilt angle of the surface from horizontal, in degrees between 0 and 180,
///          where 0 means the external surface is facing up, 90 means the external
///          surface is vertical and 180 means the external surface is facing down
/// * `a_sol`    - solar absorption coefficient at the external surface (dimensionless)
/// * `r_c`      - thermal resistance, in m2.K / W
/// * `k_m`      - areal heat capacity, in J / (m2.K)
/// * `mass_distribution_class`
///          - distribution of mass in building element, one of:
///             - 'I':  mass concentrated on internal side
///             - 'E':  mass concentrated on external side
///             - 'IE': mass divided over internal and external side
///             - 'D':  mass equally distributed
///             - 'M':  mass concentrated inside
/// * `orientation` -- is the orientation angle of the inclined surface, expressed as the
///                geographical azimuth angle of the horizontal projection of the inclined
///                surface normal, -180 to 180, in degrees
/// * `base_height` - is the distance between the ground and the lowest edge of the element, in m
/// * `height`      - is the height of the building element, in m
/// * `width`       - is the width of the building element, in m
/// * `external_conditions` -- reference to ExternalConditions object
impl BuildingElementOpaque {
    pub fn new(
        area: f64,
        is_unheated_pitched_roof: bool,
        pitch: f64,
        a_sol: f64,
        r_c: f64,
        k_m: f64,
        mass_distribution_class: MassDistributionClass,
        orientation: f64,
        base_height: f64,
        height: f64,
        width: f64,
        external_conditions: Arc<ExternalConditions>,
    ) -> Self {
        // To determine if element is an unheated pitched roof
        let (external_pitch, internal_pitch) = if is_unheated_pitched_roof {
            (pitch, 0.)
        } else {
            (pitch, pitch)
        };

        Self {
            base_height,
            width,
            projected_height: projected_height(pitch, height),
            orientation,
            external_conditions,
            area,
            r_c,
            k_m,
            a_sol,
            pitch: internal_pitch,
            external_pitch,
            therm_rad_to_sky: init_therm_rad_to_sky(sky_view_factor(&external_pitch)),
            // Calculate node conductances (h_pli) and node heat capacities (k_pli)
            // according to BS EN ISO 52016-1:2017, section 6.5.7.2
            h_pli: {
                let h_outer = 6.0 / r_c;
                let h_inner = 3.0 / r_c;

                [h_outer, h_inner, h_inner, h_outer]
            },
            k_pli: match mass_distribution_class {
                MassDistributionClass::I => [0.0, 0.0, 0.0, 0.0, k_m],
                MassDistributionClass::E => [k_m, 0.0, 0.0, 0.0, 0.0],
                MassDistributionClass::IE => {
                    let k_ie = k_m / 2.0;
                    [k_ie, 0.0, 0.0, 0.0, k_ie]
                }
                MassDistributionClass::D => {
                    let k_inner = k_m / 4.0;
                    let k_outer = k_m / 8.0;
                    [k_outer, k_inner, k_inner, k_inner, k_outer]
                }
                MassDistributionClass::M => [0.0, 0.0, k_m, 0.0, 0.0],
            },
        }
    }
}

impl BuildingElementBehaviour for BuildingElementOpaque {
    fn area(&self) -> f64 {
        self.area
    }

    fn a_sol(&self) -> f64 {
        self.a_sol
    }

    fn therm_rad_to_sky(&self) -> f64 {
        self.therm_rad_to_sky
    }

    fn pitch(&self) -> f64 {
        self.pitch
    }

    fn i_sol_dir_dif(&self, simtime: SimulationTimeIteration) -> (f64, f64) {
        let (i_sol_dir, i_sol_dif, _, _) = self
            .external_conditions
            .calculated_direct_diffuse_total_irradiance(
                self.external_pitch,
                self.orientation,
                false,
                &simtime,
            );

        (i_sol_dir, i_sol_dif)
    }

    fn shading_factors_direct_diffuse(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        self.external_conditions
            .shading_reduction_factor_direct_diffuse(
                self.base_height,
                self.projected_height,
                self.width,
                self.external_pitch,
                self.orientation,
                &Default::default(),
                simtime,
            )
    }

    fn k_pli(&self) -> &[f64] {
        &self.k_pli
    }

    fn h_pli(&self) -> &[f64] {
        &self.h_pli
    }

    fn temp_ext(&self, simtime: SimulationTimeIteration) -> f64 {
        self.external_conditions.air_temp(&simtime)
    }

    fn fabric_heat_loss(&self) -> f64 {
        let u_value = 1. / (self.r_c + self.r_se() + self.r_si());
        self.area * u_value
    }

    fn heat_capacity(&self) -> f64 {
        self.area * (self.k_m / JOULES_PER_KILOJOULE as f64)
    }
}

/// A struct to represent building elements adjacent to a thermally conditioned zone (ZTC)
#[derive(Clone, Debug)]
pub struct BuildingElementAdjacentZTC {
    area: f64,
    pitch: f64,
    k_m: f64,
    external_conditions: Arc<ExternalConditions>,
    pub a_sol: f64,
    pub therm_rad_to_sky: f64,
    h_pli: [f64; 4],
    k_pli: [f64; 5],
}

impl BuildingElementAdjacentZTC {
    /// Arguments (names based on those in BS EN ISO 52016-1:2017):
    /// * `area`     - area (in m2) of this building element
    /// * `pitch` - tilt angle of the surface from horizontal, in degrees between 0 and 180,
    ///          where 0 means the external surface is facing up, 90 means the external
    ///          surface is vertical and 180 means the external surface is facing down
    /// * `r_c`      - thermal resistance, in m2.K / W
    /// * `k_m`      - areal heat capacity, in J / (m2.K)
    /// * `mass_distribution_class`
    ///          - distribution of mass in building element, one of:
    ///             - 'I':  mass concentrated on internal side
    ///             - 'E':  mass concentrated on external side
    ///             - 'IE': mass divided over internal and external side
    ///             - 'D':  mass equally distributed
    ///             - 'M':  mass concentrated inside
    /// * `external_conditions` - reference to ExternalConditions object
    pub fn new(
        area: f64,
        pitch: f64,
        r_c: f64,
        k_m: f64,
        mass_distribution_class: MassDistributionClass,
        external_conditions: Arc<ExternalConditions>,
    ) -> Self {
        // Element is adjacent to another building / thermally conditioned zone therefore
        // according to BS EN ISO 52016-1:2017, section 6.5.6.3.6:
        // View factor to the sky is zero
        let f_sky = 0.;
        // Solar absorption coefficient at the external surface is zero
        let a_sol = 0.;

        Self {
            area,
            pitch,
            k_m,
            external_conditions,
            a_sol,
            therm_rad_to_sky: init_therm_rad_to_sky(f_sky),
            // Calculate node conductances (h_pli) and node heat capacities (k_pli)
            // according to BS EN ISO 52016-1:2017, section 6.5.7.2
            h_pli: {
                let h_outer = 6.0 / r_c;
                let h_inner = 3.0 / r_c;

                [h_outer, h_inner, h_inner, h_outer]
            },
            k_pli: match mass_distribution_class {
                MassDistributionClass::I => [0.0, 0.0, 0.0, 0.0, k_m],
                MassDistributionClass::E => [k_m, 0.0, 0.0, 0.0, 0.0],
                MassDistributionClass::IE => {
                    let k_ie = k_m / 2.0;
                    [k_ie, 0.0, 0.0, 0.0, k_ie]
                }
                MassDistributionClass::D => {
                    let k_inner = k_m / 4.0;
                    let k_outer = k_m / 8.0;
                    [k_outer, k_inner, k_inner, k_inner, k_outer]
                }
                MassDistributionClass::M => [0.0, 0.0, k_m, 0.0, 0.0],
            },
        }
    }
}

impl BuildingElementBehaviour for BuildingElementAdjacentZTC {
    fn area(&self) -> f64 {
        self.area
    }

    fn a_sol(&self) -> f64 {
        self.a_sol
    }

    fn therm_rad_to_sky(&self) -> f64 {
        self.therm_rad_to_sky
    }

    fn pitch(&self) -> f64 {
        self.pitch
    }

    fn h_ce(&self) -> f64 {
        // Element is adjacent to another building / thermally conditioned zone
        // therefore according to BS EN ISO 52016-1:2017, section 6.5.6.3.6,
        // external heat transfer coefficients are zero
        0.0
    }

    fn h_re(&self) -> f64 {
        // Element is adjacent to another building / thermally conditioned zone
        // therefore according to BS EN ISO 52016-1:2017, section 6.5.6.3.6,
        // external heat transfer coefficients are zero
        0.0
    }

    fn k_pli(&self) -> &[f64] {
        &self.k_pli
    }

    fn h_pli(&self) -> &[f64] {
        &self.h_pli
    }

    fn temp_ext(&self, simtime: SimulationTimeIteration) -> f64 {
        self.external_conditions.air_temp(&simtime)
    }

    fn fabric_heat_loss(&self) -> f64 {
        // no heat loss to thermally conditioned zones
        0.0
    }

    fn heat_capacity(&self) -> f64 {
        self.area * (self.k_m / JOULES_PER_KILOJOULE as f64)
    }
}

/// A struct to represent building elements adjacent to a thermally unconditioned zone (ZTU)
#[derive(Clone, Debug)]
pub struct BuildingElementAdjacentZTUSimple {
    area: f64,
    pitch: f64,
    r_c: f64,
    r_u: f64,
    k_m: f64,
    external_conditions: Arc<ExternalConditions>,
    h_pli: [f64; 4],
    k_pli: [f64; 5],
    pub a_sol: f64,
    pub therm_rad_to_sky: f64,
}

impl BuildingElementAdjacentZTUSimple {
    /// Arguments (names based on those in BS EN ISO 52016-1:2017):
    /// * `area`     - area (in m2) of this building element
    /// * `pitch` - tilt angle of the surface from horizontal, in degrees between 0 and 180,
    ///          where 0 means the external surface is facing up, 90 means the external
    ///          surface is vertical and 180 means the external surface is facing down
    /// * `r_c`      - thermal resistance, in m2.K / W
    /// * `r_u`      - effective thermal resistance of unheated space, in m2.K / W;
    ///             see SAP 10.2 section 3.3 for suggested values
    /// * `k_m`      - areal heat capacity, in J / (m2.K)
    /// * `mass_distribution_class`
    ///          - distribution of mass in building element, one of:
    ///             - 'I':  mass concentrated on internal side
    ///             - 'E':  mass concentrated on external side
    ///             - 'IE': mass divided over internal and external side
    ///             - 'D':  mass equally distributed
    ///             - 'M':  mass concentrated inside
    /// * `external_conditions` - reference to ExternalConditions object
    pub fn new(
        area: f64,
        pitch: f64,
        r_c: f64,
        r_u: f64,
        k_m: f64,
        mass_distribution_class: MassDistributionClass,
        external_conditions: Arc<ExternalConditions>,
    ) -> Self {
        // Element is adjacent to another building / thermally conditioned zone therefore
        // according to BS EN ISO 52016-1:2017, section 6.5.6.3.6:
        // View factor to the sky is zero
        let f_sky = 0.;
        // Solar absorption coefficient at the external surface is zero
        let a_sol = 0.;

        Self {
            area,
            pitch,
            r_c,
            r_u,
            k_m,
            external_conditions,
            // Calculate node conductances (h_pli) and node heat capacities (k_pli)
            // according to BS EN ISO 52016-1:2017, section 6.5.7.2
            h_pli: {
                let h_outer = 6.0 / r_c;
                let h_inner = 3.0 / r_c;

                [h_outer, h_inner, h_inner, h_outer]
            },
            k_pli: match mass_distribution_class {
                MassDistributionClass::I => [0.0, 0.0, 0.0, 0.0, k_m],
                MassDistributionClass::E => [k_m, 0.0, 0.0, 0.0, 0.0],
                MassDistributionClass::IE => {
                    let k_ie = k_m / 2.0;
                    [k_ie, 0.0, 0.0, 0.0, k_ie]
                }
                MassDistributionClass::D => {
                    let k_inner = k_m / 4.0;
                    let k_outer = k_m / 8.0;
                    [k_outer, k_inner, k_inner, k_inner, k_outer]
                }
                MassDistributionClass::M => [0.0, 0.0, k_m, 0.0, 0.0],
            },
            a_sol,
            therm_rad_to_sky: init_therm_rad_to_sky(f_sky),
        }
    }
}

impl BuildingElementBehaviour for BuildingElementAdjacentZTUSimple {
    fn area(&self) -> f64 {
        self.area
    }

    fn a_sol(&self) -> f64 {
        self.a_sol
    }

    fn therm_rad_to_sky(&self) -> f64 {
        self.therm_rad_to_sky
    }

    fn pitch(&self) -> f64 {
        self.pitch
    }

    fn h_ce(&self) -> f64 {
        // Add an additional thermal resistance to the outside of the wall and
        // incorporate this in the values for the external surface heat transfer
        // coefficient.
        // As this is an adjusted figure in this class, and the split between
        // h_ce and h_re does not affect the calculation results, assign entire
        // effective surface heat transfer to h_ce and set h_re to zero.
        1.0 / ((1.0 / (H_CE + H_RE)) + self.r_u)
    }

    fn h_re(&self) -> f64 {
        // As this is an adjusted figure in this class, and the split between
        // h_ce and h_re does not affect the calculation results, assign entire
        // effective surface heat transfer to h_ce and set h_re to zero.
        0.0
    }

    fn k_pli(&self) -> &[f64] {
        &self.k_pli
    }

    fn h_pli(&self) -> &[f64] {
        &self.h_pli
    }

    fn temp_ext(&self, simtime: SimulationTimeIteration) -> f64 {
        self.external_conditions.air_temp(&simtime)
    }

    fn fabric_heat_loss(&self) -> f64 {
        let u_value = 1. / (self.r_c + self.r_se() + self.r_si());
        self.area * u_value
    }

    fn heat_capacity(&self) -> f64 {
        self.area * (self.k_m / JOULES_PER_KILOJOULE as f64)
    }
}

/// A struct to represent ground building elements
#[derive(Clone, Debug)]
pub struct BuildingElementGround {
    area: f64,
    pitch: f64,
    u_value: f64,
    k_m: f64,
    h_pi: f64,
    h_pe: f64,
    perimeter: f64,
    psi_wall_floor_junc: f64,
    external_conditions: Arc<ExternalConditions>,
    temp_int_annual: f64,
    h_ce: f64,
    h_re: f64,
    pub a_sol: f64,
    pub therm_rad_to_sky: f64,
    h_pli: [f64; 4],
    k_pli: [f64; 5],
}

impl BuildingElementGround {
    /// Arguments (names based on those in BS EN ISO 52016-1:2017):
    /// * `total_area` - total area (in m2) of the building element across entire dwelling.
    ///                  If the Floor is divided among several zones,
    ///                  this is the total area across all zones.
    /// * `area`     - area (in m2) of this building element within the zone
    /// * `pitch` - tilt angle of the surface from horizontal, in degrees between 0 and 180,
    ///          where 0 means the external surface is facing up, 90 means the external
    ///          surface is vertical and 180 means the external surface is facing down
    /// * `u_value`  - steady-state thermal transmittance of floor, including the
    ///             effect of the ground, in W / (m2.K)
    ///             Calculated for the entire ground floor,
    ///             even if it is distributed among several zones.
    /// * `r_f`      - total thermal resistance of all layers in the floor construction, in (m2.K) / W
    /// * `k_m`      - areal heat capacity of the ground floor element, in J / (m2.K)
    /// * `mass_distribution_class`
    ///          - distribution of mass in building element, one of:
    ///             - 'I':  mass concentrated on internal side
    ///             - 'E':  mass concentrated on external side
    ///             - 'IE': mass divided over internal and external side
    ///             - 'D':  mass equally distributed
    ///             - 'M':  mass concentrated inside
    /// * `perimeter` - perimeter of the floor, in metres. Calculated for the entire ground floor,
    ///                 even if it is distributed among several zones.
    /// * `psi_wall_floor_junc` - linear thermal transmittance of the junction
    ///                        between the floor and the walls, in W / (m.K)
    /// * `external_conditions` - reference to ExternalConditions object
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
    pub fn new(
        total_area: f64,
        area: f64,
        pitch: f64,
        u_value: f64,
        r_f: f64,
        k_m: f64,
        mass_distribution_class: MassDistributionClass,
        floor_type: FloorType,
        edge_insulation: Option<&[EdgeInsulation]>,
        h_upper: Option<f64>,
        u_f_s: Option<f64>,
        u_w: Option<f64>,
        area_per_perimeter_vent: Option<f64>,
        shield_fact_location: Option<WindShieldLocation>,
        d_we: f64,
        r_f_ins: Option<f64>,
        z_b: Option<f64>,
        r_w_b: Option<f64>,
        h_w: Option<f64>,
        perimeter: f64,
        psi_wall_floor_junc: f64,
        external_conditions: Arc<ExternalConditions>,
    ) -> anyhow::Result<Self> {
        // Solar absorption coefficient at the external surface of the ground element is zero
        // according to BS EN ISO 52016-1:2017, section 6.5.7.3
        let a_sol = 0.0;

        // View factor to the sky is zero because element is in contact with the ground
        let f_sky = 0.0;

        let d_eq = Self::total_equiv_thickness(d_we, r_f);
        let (h_pi, h_pe) = Self::init_periodic_heat_transfer(
            floor_type,
            total_area,
            d_eq,
            perimeter,
            r_f,
            d_we,
            r_f_ins,
            h_upper,
            u_w,
            h_w,
            z_b,
            u_f_s,
            r_w_b,
            area_per_perimeter_vent,
            shield_fact_location,
            external_conditions.as_ref(),
            edge_insulation,
        )?;

        Ok(Self {
            area,
            pitch,
            u_value,
            k_m,
            h_pi,
            h_pe,
            perimeter,
            psi_wall_floor_junc,
            external_conditions,
            temp_int_annual: average_monthly_to_annual(TEMP_INT_MONTHLY_FOR_GROUND),
            h_ce: {
                // in m2.K/W
                let r_vi = (1.0 / u_value) - R_SI_FOR_GROUND - r_f - R_GR_FOR_GROUND;

                // BS EN ISO 13370:2017 Table 2 validity interval r_vi > 0
                if r_vi <= 0.0 {
                    bail!("r_vi should be greater than zero. check u-value and r_f inputs for floors - this should be checked at input validation");
                }

                1.0 / r_vi
            },
            h_re: 0.0,
            a_sol,
            therm_rad_to_sky: init_therm_rad_to_sky(f_sky),
            // Calculate node conductances (h_pli) and node heat capacities (k_pli)
            // according to BS EN ISO 52016-1:2017, section 6.5.7.2
            //
            // BS EN ISO 52016:2017 states that the r_c (resistance including the
            // effect of the ground) should be used in the equations below. However,
            // this leads to double-counting of r_si, r_gr and r_vi as these are already
            // accounted for separately, so we have used r_f (resistance of the floor
            // construction only) here instead
            h_pli: {
                let r_gr = R_GR_FOR_GROUND;

                [
                    2.0 / r_gr,
                    1.0 / (r_f / 4.0 + r_gr / 2.0),
                    2.0 / r_f,
                    4.0 / r_f,
                ]
            },
            k_pli: {
                let k_gr = K_GR_FOR_GROUND;
                match mass_distribution_class {
                    MassDistributionClass::I => [0.0, k_gr, 0.0, 0.0, k_m],
                    MassDistributionClass::E => [0.0, k_gr, k_m, 0.0, 0.0],
                    MassDistributionClass::IE => {
                        let k_ie = k_m / 2.0;
                        [0.0, k_gr, k_ie, 0.0, k_ie]
                    }
                    MassDistributionClass::D => {
                        let k_inner = k_m / 2.0;
                        let k_outer = k_m / 4.0;
                        [0.0, k_gr, k_outer, k_inner, k_outer]
                    }
                    MassDistributionClass::M => [0.0, k_gr, 0.0, k_m, 0.0],
                }
            },
        })
    }

    fn total_equiv_thickness(d_we: f64, r_f: f64) -> f64 {
        d_we + THERMAL_CONDUCTIVITY_OF_GROUND * (R_SI_FOR_GROUND + r_f + R_SE)
    }

    /// Return the periodic heat transfer coefficient for the building element
    ///             h_pi     -- Internal periodic heat transfer coefficient, in W / K
    ///                         BS EN ISO 13370:2017 Annex H
    ///             h_pe     -- external periodic heat transfer coefficient, in W / K
    ///                         BS EN ISO 13370:2017 Annex H
    fn init_periodic_heat_transfer(
        floor_type: FloorType,
        total_area: f64,
        d_eq: f64,
        perimeter: f64,
        r_f: f64,
        d_we: f64,
        r_f_ins: Option<f64>,
        h_upper: Option<f64>,
        u_w: Option<f64>,
        h_w: Option<f64>,
        z_b: Option<f64>,
        u_f_s: Option<f64>,
        r_w_b: Option<f64>,
        area_vent: Option<f64>,
        shield_fact_location: Option<WindShieldLocation>,
        external_conditions: &ExternalConditions,
        edge_insulation: Option<&[EdgeInsulation]>,
    ) -> anyhow::Result<(f64, f64)> {
        Ok(match floor_type {
            FloorType::SlabNoEdgeInsulation => {
                Self::init_slab_on_ground_floor_uninsulated_or_all_insulation(
                    total_area, d_eq, perimeter,
                )
            }
            FloorType::SlabEdgeInsulation => Self::init_slab_on_ground_floor_edge_insulated(
                total_area,
                d_eq,
                perimeter,
                edge_insulation,
            )?,
            FloorType::SuspendedFloor => Self::init_suspended_floor(
                r_f,
                d_we,
                r_f_ins,
                total_area,
                perimeter,
                h_upper,
                u_w,
                area_vent,
                shield_fact_location,
                external_conditions,
            )?,
            FloorType::HeatedBasement => {
                Self::init_heated_basement(total_area, z_b, perimeter, d_eq, r_w_b)?
            }
            FloorType::UnheatedBasement => {
                Self::init_unheated_basement(total_area, h_w, z_b, u_f_s, u_w, perimeter, d_eq)?
            }
        })
    }

    fn internal_temp_variation(total_area: f64, d_eq: f64) -> f64 {
        // H.4.1., H.5.1. Internal temperature variation
        total_area
            * (THERMAL_CONDUCTIVITY_OF_GROUND / d_eq)
            * (2. / ((1. + PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES / d_eq).powi(2) + 1.))
                .powf(0.5)
    }

    /// Slab-on-ground floor uninsulated or with all-over insulated
    fn init_slab_on_ground_floor_uninsulated_or_all_insulation(
        total_area: f64,
        d_eq: f64,
        perimeter: f64,
    ) -> (f64, f64) {
        // H.4.1. Internal temperature variation
        let h_pi = Self::internal_temp_variation(total_area, d_eq);

        // H.4.2. External temperature variation
        // 0.37 is constant in the standard but not labelled
        let h_pe = 0.37
            * perimeter
            * THERMAL_CONDUCTIVITY_OF_GROUND
            * (PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES / d_eq + 1.).ln();

        (h_pi, h_pe)
    }

    /// Slab-on-ground-with-edge-insulation
    fn init_slab_on_ground_floor_edge_insulated(
        total_area: f64,
        d_eq: f64,
        perimeter: f64,
        edge_insulation: Option<&[EdgeInsulation]>,
    ) -> anyhow::Result<(f64, f64)> {
        // H.5.1. Internal temperature variation
        let h_pi = Self::internal_temp_variation(total_area, d_eq);

        // edge insulation (vertically or horizontally)
        let h_pe = Self::edge_type(
            edge_insulation.ok_or_else(|| {
                anyhow!(
                    "Edge insulation was expected to be provided for this ground building element."
                )
            })?,
            d_eq,
            perimeter,
        )?;

        Ok((h_pi, h_pe))
    }

    fn h_pe_h(d_h: f64, r_n: f64, d_eq: f64, perimeter: f64) -> f64 {
        // horizontal edge insulation
        // 0.37 is constant in the standard but not labelled
        let eq_thick_additional = Self::add_eq_thickness(d_h, r_n);

        0.37 * perimeter
            * THERMAL_CONDUCTIVITY_OF_GROUND
            * ((1. - (-d_h / PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES).exp())
                * (PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES / (d_eq + eq_thick_additional)
                    + 1.)
                    .ln()
                + (-d_h / PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES).exp()
                    * (PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES / d_eq + 1.).ln())
    }

    fn h_pe_v(d_v: f64, r_n: f64, d_eq: f64, perimeter: f64) -> f64 {
        // vertical edge insulation
        // 0.37 is constant in the standard but not labelled
        let eq_thick_additional = Self::add_eq_thickness(d_v, r_n);

        0.37 * perimeter
            * THERMAL_CONDUCTIVITY_OF_GROUND
            * ((1. - (-2. * d_v / PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES).exp())
                * (PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES / (d_eq + eq_thick_additional)
                    + 1.)
                    .ln()
                + (-2. * d_v / PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES).exp()
                    * (PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES / d_eq + 1.).ln())
    }

    /// Additional equivalent thickness
    fn add_eq_thickness(d_n: f64, r_n: f64) -> f64 {
        // m2·K/W, thermal resistance
        let r_add_eq = r_n - d_n / THERMAL_CONDUCTIVITY_OF_GROUND;
        // m, thickness_edge-insulation or foundation
        r_add_eq * THERMAL_CONDUCTIVITY_OF_GROUND
    }

    /// edge insulation vertically or horizontally
    fn edge_type(
        edge_insulation: &[EdgeInsulation],
        d_eq: f64,
        perimeter: f64,
    ) -> anyhow::Result<f64> {
        // Initialise edge width and depth
        let mut h_pe_list = vec![];
        for edge in edge_insulation {
            match edge {
                EdgeInsulation::Horizontal {
                    width,
                    edge_thermal_resistance,
                } => {
                    h_pe_list.push(Self::h_pe_h(
                        *width,
                        *edge_thermal_resistance,
                        d_eq,
                        perimeter,
                    ));
                }
                EdgeInsulation::Vertical {
                    depth,
                    edge_thermal_resistance,
                } => {
                    h_pe_list.push(Self::h_pe_v(
                        *depth,
                        *edge_thermal_resistance,
                        d_eq,
                        perimeter,
                    ));
                }
            }
        }

        h_pe_list
            .iter()
            .min_by(|a, b| a.total_cmp(b))
            .ok_or_else(|| {
                anyhow!("There was no edge insulation provided for a ground building element.")
            })
            .copied()
    }

    /// Suspended floor periodic coefficients
    fn init_suspended_floor(
        r_f: f64,
        d_we: f64,
        r_f_ins: Option<f64>,
        total_area: f64,
        perimeter: f64,
        h_upper: Option<f64>,
        u_w: Option<f64>,
        area_vent: Option<f64>,
        shield_fact_location: Option<WindShieldLocation>,
        external_conditions: &ExternalConditions,
    ) -> anyhow::Result<(f64, f64)> {
        // H.6.1.
        // thermal transmittance of suspended part of floor, in W/(m2·K)
        let u_f = Self::thermal_transmittance_sus_floor(r_f, R_SI_FOR_GROUND);

        // equivalent thermal transmittance, in W/(m2·K)
        let u_x = Self::equiv_therma_trans(
            total_area,
            perimeter,
            h_upper,
            u_w,
            area_vent,
            shield_fact_location,
            external_conditions,
        )?;

        // equivalent thickness, in m
        let d_g = Self::total_equiv_thickness_sus(d_we, r_f_ins)?;

        // H.6.2. Internal temperature variation
        let h_pi = total_area
            * (1. / u_f
                + 1. / (THERMAL_CONDUCTIVITY_OF_GROUND
                    / PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES
                    + u_x));

        // H.6.3. External temperature variation
        // 0.37 is constant in the standard but not labelled
        let h_pe = u_f
            * ((0.37
                * perimeter
                * THERMAL_CONDUCTIVITY_OF_GROUND
                * (PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES / d_g + 1.).ln()
                + u_x * total_area)
                / (THERMAL_CONDUCTIVITY_OF_GROUND
                    / (PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES + u_x + u_f)));

        Ok((h_pi, h_pe))
    }

    /// thermal transmittance of suspended part of floor
    fn thermal_transmittance_sus_floor(r_f: f64, r_si: f64) -> f64 {
        1. / (r_f + 2. * r_si)
    }

    /// equivalent thermal transmittance between the underfloor space and the outside
    fn equiv_therma_trans(
        total_area: f64,
        perimeter: f64,
        h_upper: Option<f64>,
        u_w: Option<f64>,
        area_vent: Option<f64>,
        shield_fact_location: Option<WindShieldLocation>,
        external_conditions: &ExternalConditions,
    ) -> anyhow::Result<f64> {
        let h_upper = h_upper.ok_or_else(|| anyhow!("A value for the height of the floor upper surface was needed for a ground building element."))?;
        let u_w = u_w.ok_or_else(|| anyhow!("A value for the thermal transmittance of walls above ground was needed for a ground building element."))?;
        let area_vent = area_vent.ok_or_else(|| {
            anyhow!("An area per perimeter vent value was needed for a ground building element.")
        })?;
        let shield_fact_location = shield_fact_location.ok_or_else(|| {
            anyhow!("A wind shielding factor indicator was needed for a ground building element.")
        })?;

        // Characteristic dimension of floor
        let char_dimen = Self::charac_dimen_floor(total_area, perimeter);

        // 1450 is constant in the standard but not labelled
        Ok(2. * (h_upper * u_w / char_dimen)
            + 1450.
                * (area_vent
                    * Self::wind_speed(external_conditions)?
                    * Self::wind_shield_fact(shield_fact_location))
                / char_dimen)
    }

    /// Equivalent thickness for the ground
    fn total_equiv_thickness_sus(d_we: f64, r_f_ins: Option<f64>) -> anyhow::Result<f64> {
        Ok(d_we + THERMAL_CONDUCTIVITY_OF_GROUND * (R_SI_FOR_GROUND + r_f_ins.ok_or_else(|| anyhow!("A value for thermal resistance of insulation was expected to be given for a ground building element."))? + R_SE))
    }

    /// Characteristic dimension of floor, in metres
    fn charac_dimen_floor(total_area: f64, perimeter: f64) -> f64 {
        total_area / (0.5 * perimeter)
    }

    fn wind_speed(external_conditions: &ExternalConditions) -> anyhow::Result<f64> {
        external_conditions.wind_speed_annual().ok_or_else(|| anyhow!("Annual wind speed was expected to be available when instantiating a ground building element for a calculation."))
    }

    /// wind shielding factor
    fn wind_shield_fact(shield_fact_location: WindShieldLocation) -> f64 {
        // Values from BS EN ISO 13370:2017 Table 8
        match shield_fact_location {
            WindShieldLocation::Sheltered => 0.02,
            WindShieldLocation::Average => 0.05,
            WindShieldLocation::Exposed => 0.10,
        }
    }

    /// Heated basement periodic coefficients
    fn init_heated_basement(
        total_area: f64,
        z_b: Option<f64>,
        perimeter: f64,
        d_eq: f64,
        r_w_b: Option<f64>,
    ) -> anyhow::Result<(f64, f64)> {
        let z_b = z_b.ok_or_else(|| anyhow!("A value for the depth of the basement floor below ground level was needed for this ground building element."))?;

        // total equivalent thickness
        let d_w_b = Self::equiv_thick_base_wall(r_w_b)?;

        // H.7.1. Internal temperature variation
        let h_pi = total_area
            * ((THERMAL_CONDUCTIVITY_OF_GROUND / d_eq)
                * (2. / (1. + PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES / d_eq).powi(2)
                    + 1.)
                    .powf(0.5))
            + z_b
                * perimeter
                * (THERMAL_CONDUCTIVITY_OF_GROUND / d_w_b)
                * (2.
                    / ((1. + PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES / d_w_b).powi(2)
                        + 1.))
                    .powf(0.5);

        // H.7.2. External temperature variation
        // 0.37 is constant in the standard but not labelled
        let h_pe = 0.37
            * perimeter
            * THERMAL_CONDUCTIVITY_OF_GROUND
            * ((-z_b / PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES).exp()
                * (PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES / d_eq + 1.).ln()
                + 2. * (1. - (-z_b / PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES).exp())
                    * (PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES / d_w_b + 1.).ln());

        Ok((h_pi, h_pe))
    }

    /// Equivalent thickness for the basement walls
    fn equiv_thick_base_wall(r_w_b: Option<f64>) -> anyhow::Result<f64> {
        // r_w_b is the thermal resistance of the walls
        Ok(THERMAL_CONDUCTIVITY_OF_GROUND * (R_SI_FOR_GROUND + r_w_b.ok_or_else(|| anyhow!("A value for thermal resistance of walls for the basement is needed for this ground building element."))? + R_SE))
    }

    /// Unheated basement
    fn init_unheated_basement(
        total_area: f64,
        h_w: Option<f64>,
        z_b: Option<f64>,
        u_f_s: Option<f64>,
        u_w: Option<f64>,
        perimeter: f64,
        d_eq: f64,
    ) -> anyhow::Result<(f64, f64)> {
        let h_w = h_w.ok_or_else(|| anyhow!("A value for the height of the basement walls above ground level was needed for this ground building element."))?;
        let z_b = z_b.ok_or_else(|| anyhow!("A value for the depth of the basement floor below ground level was needed for this ground building element."))?;
        let u_f_s = u_f_s.ok_or_else(|| anyhow!("A value for thermal transmittance of floor above basement was needed for this ground building element."))?;
        let u_w = u_w.ok_or_else(|| anyhow!("A value for the thermal transmittance of walls above ground was needed for this ground building element."))?;

        // Wh/(m3·K)
        let thermal_capacity_air = 0.33;
        let air_vol_base = total_area * (h_w + z_b);

        // air changes per hour
        // From BS EN ISO 13370:2017 section 7.4
        let vent_rate_base = 0.3;

        // H.8.1. Internal temperature variation
        let h_pi = (1. / (total_area * u_f_s)
            + 1. / ((total_area + z_b * perimeter) * THERMAL_CONDUCTIVITY_OF_GROUND
                / PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES
                + h_w * perimeter * u_w
                + thermal_capacity_air * vent_rate_base * air_vol_base))
            .powi(-1);

        // H.8.2. External temperature variation
        // 0.37 is constant in the standard but not labelled
        let h_pe = total_area
            * u_f_s
            * (0.37
                * perimeter
                * THERMAL_CONDUCTIVITY_OF_GROUND
                * (2. - (-z_b / PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES).exp())
                * (PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES / d_eq + 1.).ln()
                + h_w * perimeter * u_w
                + thermal_capacity_air * vent_rate_base * air_vol_base)
            / ((total_area + z_b * perimeter) * THERMAL_CONDUCTIVITY_OF_GROUND
                / PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES
                + h_w * perimeter * u_w
                + thermal_capacity_air * vent_rate_base * air_vol_base
                + total_area * u_f_s);

        Ok((h_pi, h_pe))
    }
}

impl BuildingElementBehaviour for BuildingElementGround {
    fn area(&self) -> f64 {
        self.area
    }

    fn a_sol(&self) -> f64 {
        self.a_sol
    }

    fn therm_rad_to_sky(&self) -> f64 {
        self.therm_rad_to_sky
    }

    fn pitch(&self) -> f64 {
        self.pitch
    }

    fn h_ce(&self) -> f64 {
        self.h_ce
    }

    fn h_re(&self) -> f64 {
        self.h_re
    }

    fn k_pli(&self) -> &[f64] {
        &self.k_pli
    }

    fn h_pli(&self) -> &[f64] {
        &self.h_pli
    }

    fn temp_ext(&self, simtime: SimulationTimeIteration) -> f64 {
        let temp_ext_annual = self
            .external_conditions
            .air_temp_annual()
            .expect("no annual air temp available");
        let temp_ext_month = self
            .external_conditions
            .air_temp_monthly(simtime.current_month_start_end_hours());

        let current_month = simtime.current_month().unwrap_or(0);
        let temp_int_month = TEMP_INT_MONTHLY_FOR_GROUND[current_month as usize];

        // BS EN ISO 13370:2017 Eqn C.4
        let heat_flow_month = self.u_value * self.area * (self.temp_int_annual - temp_ext_annual)
            + self.perimeter * self.psi_wall_floor_junc * (temp_int_month - temp_ext_month)
            - self.h_pi * (self.temp_int_annual - temp_int_month)
            + self.h_pe * (temp_ext_annual - temp_ext_month);

        // BS EN ISO 13370:2017 Eqn F.2
        temp_int_month
            - (heat_flow_month
                - (self.perimeter
                    * self.psi_wall_floor_junc
                    * (self.temp_int_annual - temp_ext_annual)))
                / (self.area * self.u_value)
    }

    fn fabric_heat_loss(&self) -> f64 {
        self.u_value * self.area
    }

    fn heat_capacity(&self) -> f64 {
        self.area * (self.k_m / JOULES_PER_KILOJOULE as f64)
    }
}

/// A class to represent transparent building elements (windows etc.)
#[derive(Clone, Debug)]
pub struct BuildingElementTransparent {
    area: f64,
    r_c: f64,
    pitch: f64,
    orientation: f64,
    g_value: f64,
    frame_area_fraction: f64,
    base_height: f64,
    projected_height: f64,
    mid_height: f64,
    width: f64,
    shading: Vec<WindowShadingObject>,
    h_pli: [f64; 1],
    k_pli: [f64; 2],
    pub a_sol: f64,
    pub therm_rad_to_sky: f64,
    external_conditions: Arc<ExternalConditions>,
}

impl BuildingElementTransparent {
    pub fn new(
        pitch: f64,
        r_c: f64,
        orientation: f64,
        g_value: f64,
        frame_area_fraction: f64,
        base_height: f64,
        height: f64,
        width: f64,
        shading: Vec<WindowShadingObject>,
        treatment: Option<Vec<WindowTreatment>>,
        external_conditions: Arc<ExternalConditions>,
        simulation_time: &SimulationTimeIterator,
    ) -> Self {
        // Solar absorption coefficient is zero because element is transparent
        let a_sol = 0.0;

        // This is the f_sky value for an unshaded surface
        let f_sky = sky_view_factor(&pitch);

        Self {
            area: height * width,
            r_c,
            pitch,
            orientation,
            g_value,
            frame_area_fraction,
            base_height,
            projected_height: projected_height(pitch, height),
            mid_height: base_height + height / 2.0,
            width,
            shading,
            h_pli: [1.0 / r_c],
            k_pli: [0.0, 0.0],
            a_sol,
            therm_rad_to_sky: init_therm_rad_to_sky(f_sky),
            external_conditions,
        }
    }

    /// return g_value corrected for angle of solar radiation
    fn convert_g_value(&self) -> f64 {
        // TODO (from Python) for windows with scattering glazing or solar shading provisions
        // there is a different, more complex method for conversion that depends on
        // timestep (via solar altitude).
        // suggest this is implemented at the same time as window shading (devices
        // rather than fixed features) as will also need to link to shading schedule.
        // see ISO 52016 App E. Page 177
        // How do we know whether a window has "scattering glazing"?
        //
        // g_value = agl * g_alt + (1 - agl) * g_dif

        let fw = 0.90;
        // default from ISO 52016 App B Table B.22
        fw * self.g_value
    }

    pub fn projected_height(&self) -> f64 {
        self.projected_height
    }

    pub fn mid_height(&self) -> f64 {
        self.mid_height
    }

    pub fn orientation(&self) -> f64 {
        self.orientation
    }
}

impl BuildingElementBehaviour for BuildingElementTransparent {
    fn area(&self) -> f64 {
        self.area
    }

    fn a_sol(&self) -> f64 {
        self.a_sol
    }

    fn therm_rad_to_sky(&self) -> f64 {
        self.therm_rad_to_sky
    }

    fn pitch(&self) -> f64 {
        self.pitch
    }

    fn solar_gains(&self, simtime: SimulationTimeIteration) -> anyhow::Result<f64> {
        let (i_sol_dir, i_sol_dif, _, _) = self
            .external_conditions
            .calculated_direct_diffuse_total_irradiance(
                self.pitch,
                self.orientation,
                false,
                &simtime,
            );
        let g_value = self.convert_g_value();

        let (f_sh_dir, f_sh_dif) = self.shading_factors_direct_diffuse(simtime)?;
        Ok(g_value
            * (i_sol_dif * f_sh_dif + i_sol_dir * f_sh_dir)
            * self.area
            * (1. - self.frame_area_fraction))
    }

    fn shading_factors_direct_diffuse(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        self.external_conditions
            .shading_reduction_factor_direct_diffuse(
                self.base_height,
                self.projected_height,
                self.width,
                self.pitch,
                self.orientation,
                &self.shading,
                simtime,
            )
    }

    fn k_pli(&self) -> &[f64] {
        &self.k_pli
    }

    fn h_pli(&self) -> &[f64] {
        &self.h_pli
    }

    fn h_pli_by_index_unchecked(&self, idx: usize) -> f64 {
        // Account for resistance of window treatment in heat transfer coefficient
        // TODO (from Python) Check that idx refers to inside surface?
        let r_c = 1.0 / self.h_pli[idx];
        // TODO: port the below as part of migration to 0.32
        //        if self._treatment is not None:
        //             for t in self._treatment:
        //                 self._adjust_treatment()
        //                 if t['is_open'] == False:
        //                     # delta_r for window treatment (curtains, blinds, etc.) as per
        //                     # BS EN 13125:2001
        //                     r_c += t['delta_r']
        1.0 / r_c
    }

    fn temp_ext(&self, simtime: SimulationTimeIteration) -> f64 {
        self.external_conditions.air_temp(&simtime)
    }

    fn fabric_heat_loss(&self) -> f64 {
        // Effective window U-value includes assumed use of curtains/blinds, see
        // SAP10.2 spec, paragraph 3.2
        // TODO (from Python) Confirm this is still the desired approach for SAP 11
        let r_curtains_blinds = 0.04;
        let u_value = 1. / ((self.r_c + self.r_si() + self.r_se()) + r_curtains_blinds);
        self.area * u_value
    }

    fn heat_capacity(&self) -> f64 {
        // Set to zero as not included in heat loss calculations
        0.0
    }
}

pub struct NamedBuildingElementTransparent {
    pub name: String,
    pub window: BuildingElementTransparent,
}

// equality and hashing based on name for identity

impl PartialEq for NamedBuildingElementTransparent {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
    }
}

impl Eq for NamedBuildingElementTransparent {}

impl Hash for NamedBuildingElementTransparent {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.name.hash(state);
    }
}

pub fn area_for_building_element_input(element: &BuildingElementInput) -> f64 {
    match *element {
        BuildingElementInput::Opaque { area: a, .. } => a,
        BuildingElementInput::Transparent { height, width, .. } => height * width,
        BuildingElementInput::Ground { area: a, .. } => a,
        BuildingElementInput::AdjacentZTC { area: a, .. } => a,
        BuildingElementInput::AdjacentZTUSimple { area: a, .. } => a,
    }
}

/// calc the vertically projected height of a surface from
/// the actual height and tilt of the surface
pub fn projected_height(tilt: f64, height: f64) -> f64 {
    let mut ph = height * tilt.to_radians().sin();
    // BS EN ISO 52010-1 Table 7 geometric input data; shading. Footnote d
    // validity interval H1;ic > 0
    // if horizontal (height = 0): choose small value e.g. H1 = 0.01 m"""
    if ph < 0.01 {
        ph = 0.01;
    }

    ph
}

pub fn convert_uvalue_to_resistance(u_value: f64, pitch: f64) -> f64 {
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

/// Calculate longwave sky view factor from pitch in degrees
pub(crate) fn sky_view_factor(pitch: &f64) -> f64 {
    // TODO (from Python) account for shading
    // TODO (from Python) check longwave is correct
    let pitch_rads = pitch * PI / 180.0;

    0.5 * (1.0 + pitch_rads.cos())
}

#[derive(Debug, PartialEq)]
pub enum HeatFlowDirection {
    Horizontal,
    Upwards,
    Downwards,
}

// Thermal properties of ground from BS EN ISO 13370:2017 Table 7
// Use values for clay or silt (same as BR 443 and SAP 10)
const THERMAL_CONDUCTIVITY_OF_GROUND: f64 = 1.5;
// in W/(m.K)
const HEAT_CAPACITY_PER_VOLUME_OF_GROUND: f64 = 3_000_000.;

// Periodic penetration depth of ground from BS EN ISO 13370:2017 Table H.1
// Use values for clay or silt (same as BR 443 and SAP 10)
const PERIODIC_PENETRATION_DEPTH_FOR_GROUND_IN_METRES: f64 = 2.2;

// in J/(m3.K)
const THICKNESS_GROUND_LAYER: f64 = 0.5; // in m. Specified in BS EN ISO 52016-1:2017 section 6.5.8.2

// thermal resistance in (m2.K)/W
const R_GR_FOR_GROUND: f64 = THICKNESS_GROUND_LAYER / THERMAL_CONDUCTIVITY_OF_GROUND;
// areal heat capacity in J/(m2.K)
const K_GR_FOR_GROUND: f64 = THICKNESS_GROUND_LAYER * HEAT_CAPACITY_PER_VOLUME_OF_GROUND;

const R_SI_FOR_GROUND: f64 = 0.17; // ISO 6946 - internal surface resistance

// Assume values for temp_int_annual and temp_int_monthly
// These are based on SAP 10 notional building runs for 5 archetypes used
// for inter-model comparison/validation. The average of the monthly mean
// internal temperatures from each run was taken.
const TEMP_INT_MONTHLY_FOR_GROUND: [f64; 12] = [
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external_conditions::DaylightSavingsConfig;
    use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
    use approx::assert_relative_eq;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    pub fn simulation_time() -> SimulationTimeIterator {
        SimulationTime::new(0.0, 4.0, 1.0).iter()
    }

    #[fixture]
    pub fn external_conditions(simulation_time: SimulationTimeIterator) -> Arc<ExternalConditions> {
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
            vec![],
        ))
    }

    #[fixture]
    pub fn be_i(external_conditions: Arc<ExternalConditions>) -> BuildingElementOpaque {
        BuildingElementOpaque::new(
            20.,
            false,
            180.,
            0.60,
            0.25,
            19000.0,
            MassDistributionClass::I,
            0.,
            0.,
            2.,
            10.,
            external_conditions,
        )
    }

    #[fixture]
    pub fn be_e(external_conditions: Arc<ExternalConditions>) -> BuildingElementOpaque {
        BuildingElementOpaque::new(
            22.5,
            false,
            135.,
            0.61,
            0.50,
            18000.0,
            MassDistributionClass::E,
            180.,
            0.,
            2.25,
            10.,
            external_conditions,
        )
    }

    #[fixture]
    pub fn be_ie(external_conditions: Arc<ExternalConditions>) -> BuildingElementOpaque {
        BuildingElementOpaque::new(
            25.,
            false,
            90.,
            0.62,
            0.75,
            17000.0,
            MassDistributionClass::IE,
            90.,
            0.,
            2.5,
            10.,
            external_conditions,
        )
    }

    #[fixture]
    pub fn be_d(external_conditions: Arc<ExternalConditions>) -> BuildingElementOpaque {
        BuildingElementOpaque::new(
            27.5,
            true,
            45.,
            0.63,
            0.80,
            16000.0,
            MassDistributionClass::D,
            -90.,
            0.,
            2.75,
            10.,
            external_conditions,
        )
    }

    #[fixture]
    pub fn be_m(external_conditions: Arc<ExternalConditions>) -> BuildingElementOpaque {
        BuildingElementOpaque::new(
            30.,
            false,
            0.,
            0.64,
            0.40,
            15000.0,
            MassDistributionClass::M,
            0.,
            0.,
            3.,
            10.,
            external_conditions,
        )
    }

    #[fixture]
    pub fn opaque_building_elements(
        be_i: BuildingElementOpaque,
        be_e: BuildingElementOpaque,
        be_ie: BuildingElementOpaque,
        be_d: BuildingElementOpaque,
        be_m: BuildingElementOpaque,
    ) -> [BuildingElementOpaque; 5] {
        [be_i, be_e, be_ie, be_d, be_m]
    }

    #[rstest]
    pub fn test_no_of_nodes_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
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
    pub fn test_area_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
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
    pub fn test_heat_flow_direction_for_opaque(
        opaque_building_elements: [BuildingElementOpaque; 5],
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
        for (i, be) in opaque_building_elements.iter().enumerate() {
            assert_eq!(
                be.heat_flow_direction(temp_int_air, temp_int_surface[i]),
                results[i],
                "incorrect heat flow direction returned"
            );
        }
    }

    #[rstest]
    pub fn test_r_si_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        let results = [0.17, 0.17, 0.13, 0.10, 0.10];

        for (i, be) in opaque_building_elements.iter().enumerate() {
            assert_relative_eq!(be.r_si(), results[i], max_relative = 0.05);
        }
    }

    #[rstest]
    pub fn test_h_ci_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
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
    pub fn test_h_ri_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        for be in opaque_building_elements.iter() {
            assert_relative_eq!(be.h_ri(), 5.13,);
        }
    }

    #[rstest]
    pub fn test_h_ce_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        for be in opaque_building_elements.iter() {
            assert_relative_eq!(be.h_ce(), 20.0,);
        }
    }

    #[rstest]
    pub fn test_h_re(opaque_building_elements: [BuildingElementOpaque; 5]) {
        for be in opaque_building_elements.iter() {
            assert_relative_eq!(be.h_re(), 4.14,);
        }
    }

    #[rstest]
    pub fn test_a_sol_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        // Define increment between test cases
        let a_sol_inc = 0.01;

        for (i, be) in opaque_building_elements.iter().enumerate() {
            assert_relative_eq!(be.a_sol(), 0.6 + i as f64 * a_sol_inc,);
        }
    }

    #[rstest]
    pub fn test_therm_rad_to_sky_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        let results = [0.0, 6.6691785923823135, 22.77, 38.87082140761768, 45.54];

        for (i, be) in opaque_building_elements.iter().enumerate() {
            assert_relative_eq!(be.therm_rad_to_sky(), results[i],);
        }
    }

    #[rstest]
    pub fn test_h_pli_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
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
    pub fn test_temp_ext_for_opaque(
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

    #[ignore = "the assertion values here cause failures - upstream fix has been committed to"]
    #[rstest]
    pub fn test_fabric_heat_loss_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
        let results = [43.20, 35.15, 27.10, 27.15, 55.54];

        for (i, be) in opaque_building_elements.iter().enumerate() {
            assert_relative_eq!(be.fabric_heat_loss(), results[i], max_relative = 1e-2);
        }
    }

    #[rstest]
    pub fn test_heat_capacity_for_opaque(opaque_building_elements: [BuildingElementOpaque; 5]) {
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
    pub fn adjacent_ztc_building_elements(
        external_conditions: Arc<ExternalConditions>,
    ) -> [BuildingElementAdjacentZTC; 5] {
        let be_i = BuildingElementAdjacentZTC::new(
            20.0,
            180.,
            0.25,
            19000.0,
            MassDistributionClass::I,
            external_conditions.clone(),
        );
        let be_e = BuildingElementAdjacentZTC::new(
            22.5,
            135.,
            0.50,
            18000.0,
            MassDistributionClass::E,
            external_conditions.clone(),
        );
        let be_ie = BuildingElementAdjacentZTC::new(
            25.0,
            90.,
            0.75,
            17000.0,
            MassDistributionClass::IE,
            external_conditions.clone(),
        );
        let be_d = BuildingElementAdjacentZTC::new(
            27.5,
            45.,
            0.80,
            16000.0,
            MassDistributionClass::D,
            external_conditions.clone(),
        );
        let be_m = BuildingElementAdjacentZTC::new(
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
    pub fn test_no_of_nodes_for_adjacent_ztc(
        adjacent_ztc_building_elements: [BuildingElementAdjacentZTC; 5],
    ) {
        for be in adjacent_ztc_building_elements {
            assert_eq!(be.number_of_nodes(), 5, "incorrect number of nodes");
            assert_eq!(
                be.number_of_inside_nodes(),
                3,
                "incorrect number of inside nodes"
            );
        }
    }

    #[rstest]
    pub fn test_area_for_adjacent_ztc(
        adjacent_ztc_building_elements: [BuildingElementAdjacentZTC; 5],
    ) {
        // Define increment between test cases
        let area_inc = 2.5;

        for (i, be) in adjacent_ztc_building_elements.iter().enumerate() {
            assert_relative_eq!(be.area(), 20.0 + i as f64 * area_inc,);
        }
    }

    #[rstest]
    pub fn test_heat_flow_direction_for_adjacent_ztc(
        adjacent_ztc_building_elements: [BuildingElementAdjacentZTC; 5],
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

        for (i, be) in adjacent_ztc_building_elements.iter().enumerate() {
            assert_eq!(
                be.heat_flow_direction(temp_int_air, temp_int_surface[i]),
                results[i],
                "incorrect heat flow direction returned"
            );
        }
    }

    #[rstest]
    pub fn test_r_si_for_adjacent_ztc(
        adjacent_ztc_building_elements: [BuildingElementAdjacentZTC; 5],
    ) {
        let results = [0.17, 0.17, 0.13, 0.10, 0.10];

        for (i, be) in adjacent_ztc_building_elements.iter().enumerate() {
            assert_relative_eq!(be.r_si(), results[i], max_relative = 0.05);
        }
    }

    #[rstest]
    pub fn test_h_ci_for_adjacent_ztc(
        adjacent_ztc_building_elements: [BuildingElementAdjacentZTC; 5],
    ) {
        let temp_int_air = 20.0;
        let temp_int_surface = [19.0, 21.0, 22.0, 21.0, 19.0];
        let results = [0.7, 5.0, 2.5, 0.7, 5.0];

        for (i, be) in adjacent_ztc_building_elements.iter().enumerate() {
            assert_relative_eq!(be.h_ci(temp_int_air, temp_int_surface[i]), results[i],);
        }
    }

    #[rstest]
    pub fn test_h_ri_for_adjacent_ztc(
        adjacent_ztc_building_elements: [BuildingElementAdjacentZTC; 5],
    ) {
        for be in adjacent_ztc_building_elements {
            assert_relative_eq!(be.h_ri(), 5.13,);
        }
    }

    #[rstest]
    pub fn test_h_ce_for_adjacent_ztc(
        adjacent_ztc_building_elements: [BuildingElementAdjacentZTC; 5],
    ) {
        for be in adjacent_ztc_building_elements {
            assert_relative_eq!(be.h_ce(), 0.0,);
        }
    }

    #[rstest]
    pub fn test_h_re_for_adjacent_ztc(
        adjacent_ztc_building_elements: [BuildingElementAdjacentZTC; 5],
    ) {
        for be in adjacent_ztc_building_elements {
            assert_eq!(be.h_re(), 0.0, "incorrect h_re returned");
        }
    }

    #[rstest]
    pub fn test_a_sol_for_adjacent_ztc(
        adjacent_ztc_building_elements: [BuildingElementAdjacentZTC; 5],
    ) {
        for be in adjacent_ztc_building_elements {
            assert_eq!(be.a_sol(), 0.0, "incorrect a_sol returned");
        }
    }

    #[rstest]
    pub fn test_therm_rad_to_sky_for_adjacent_ztc(
        adjacent_ztc_building_elements: [BuildingElementAdjacentZTC; 5],
    ) {
        for be in adjacent_ztc_building_elements {
            assert_eq!(be.therm_rad_to_sky(), 0.0, "incorrect a_sol returned");
        }
    }

    #[rstest]
    fn test_h_pli_for_adjacent_ztc(
        adjacent_ztc_building_elements: [BuildingElementAdjacentZTC; 5],
    ) {
        let results = [
            [24.0, 12.0, 12.0, 24.0],
            [12.0, 6.0, 6.0, 12.0],
            [8.0, 4.0, 4.0, 8.0],
            [7.5, 3.75, 3.75, 7.5],
            [15.0, 7.5, 7.5, 15.0],
        ];
        for (i, be) in adjacent_ztc_building_elements.iter().enumerate() {
            assert_eq!(be.h_pli(), &results[i], "incorrect h_pli list returned");
        }
    }

    #[rstest]
    pub fn test_k_pli_for_adjacent_ztc(
        adjacent_ztc_building_elements: [BuildingElementAdjacentZTC; 5],
    ) {
        let results = [
            [0.0, 0.0, 0.0, 0.0, 19000.0],
            [18000.0, 0.0, 0.0, 0.0, 0.0],
            [8500.0, 0.0, 0.0, 0.0, 8500.0],
            [2000.0, 4000.0, 4000.0, 4000.0, 2000.0],
            [0.0, 0.0, 15000.0, 0.0, 0.0],
        ];
        for (i, be) in adjacent_ztc_building_elements.iter().enumerate() {
            assert_eq!(be.k_pli(), &results[i], "incorrect k_pli list returned");
        }
    }

    #[rstest]
    pub fn test_fabric_heat_loss_for_adjacent_ztc(
        adjacent_ztc_building_elements: [BuildingElementAdjacentZTC; 5],
    ) {
        for be in adjacent_ztc_building_elements {
            assert_eq!(
                be.fabric_heat_loss(),
                0.0,
                "incorrect fabric heat loss returned"
            );
        }
    }

    #[rstest]
    pub fn test_heat_capacity_for_adjacent_ztc(
        adjacent_ztc_building_elements: [BuildingElementAdjacentZTC; 5],
    ) {
        let results = [380., 405., 425., 440., 450.];

        for (i, be) in adjacent_ztc_building_elements.iter().enumerate() {
            assert_eq!(
                be.heat_capacity(),
                results[i],
                "incorrect heat capacity returned"
            );
        }
    }

    #[fixture]
    pub fn ground_building_elements(
        external_conditions_for_ground: Arc<ExternalConditions>,
    ) -> [BuildingElementGround; 5] {
        let be_i = BuildingElementGround::new(
            20.0,
            20.0,
            180.,
            1.5,
            0.1,
            19000.0,
            MassDistributionClass::I,
            FloorType::SuspendedFloor,
            None,
            Some(0.5),
            None,
            Some(0.5),
            Some(0.01),
            Some(WindShieldLocation::Sheltered),
            0.3,
            Some(7.),
            None,
            None,
            None,
            18.0,
            0.5,
            external_conditions_for_ground.clone(),
        )
        .unwrap();
        let be_e = BuildingElementGround::new(
            22.5,
            22.5,
            135.,
            1.4,
            0.2,
            18000.0,
            MassDistributionClass::E,
            FloorType::SlabNoEdgeInsulation,
            None,
            None,
            None,
            None,
            None,
            None,
            0.3,
            None,
            None,
            None,
            None,
            19.0,
            0.6,
            external_conditions_for_ground.clone(),
        )
        .unwrap();
        let edge_insulation_ie = [
            EdgeInsulation::Horizontal {
                width: 3.0,
                edge_thermal_resistance: 2.0,
            },
            EdgeInsulation::Vertical {
                depth: 1.0,
                edge_thermal_resistance: 2.0,
            },
        ];
        let be_ie = BuildingElementGround::new(
            25.0,
            25.0,
            90.,
            1.33,
            0.2,
            17000.0,
            MassDistributionClass::IE,
            FloorType::SlabEdgeInsulation,
            Some(&edge_insulation_ie),
            None,
            None,
            None,
            None,
            None,
            0.3,
            None,
            None,
            None,
            None,
            20.0,
            0.7,
            external_conditions_for_ground.clone(),
        )
        .unwrap();
        let be_d = BuildingElementGround::new(
            27.5,
            27.5,
            45.,
            1.25,
            0.2,
            16000.0,
            MassDistributionClass::D,
            FloorType::HeatedBasement,
            None,
            None,
            None,
            None,
            None,
            None,
            0.3,
            None,
            Some(2.3),
            Some(6.),
            None,
            21.0,
            0.8,
            external_conditions_for_ground.clone(),
        )
        .unwrap();
        let be_m = BuildingElementGround::new(
            30.0,
            30.0,
            0.,
            1.0,
            0.3,
            15000.0,
            MassDistributionClass::M,
            FloorType::UnheatedBasement,
            None,
            None,
            Some(1.2),
            Some(0.5),
            None,
            None,
            0.3,
            None,
            Some(2.3),
            Some(0.15),
            Some(2.3),
            22.0,
            0.9,
            external_conditions_for_ground,
        )
        .unwrap();
        [be_i, be_e, be_ie, be_d, be_m]
    }

    #[fixture]
    pub fn simulation_time_for_ground() -> SimulationTime {
        SimulationTime::new(742., 746., 1.)
    }

    #[fixture]
    pub fn external_conditions_for_ground(
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
            vec![180.; 8760],
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
            vec![],
        ))
    }

    #[rstest]
    pub fn test_no_of_nodes_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
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
    pub fn test_area_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        // Define increment between test cases
        let area_inc = 2.5;

        for (i, be) in ground_building_elements.iter().enumerate() {
            assert_eq!(be.area(), 20.0 + i as f64 * area_inc,);
        }
    }

    #[rstest]
    pub fn test_heat_flow_direction_for_ground(
        ground_building_elements: [BuildingElementGround; 5],
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

        for (i, be) in ground_building_elements.iter().enumerate() {
            assert_eq!(
                be.heat_flow_direction(temp_int_air, temp_int_surface[i]),
                results[i],
                "incorrect heat flow direction returned"
            );
        }
    }

    #[rstest]
    pub fn test_r_si_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        let results = [0.17, 0.17, 0.13, 0.10, 0.10];

        for (i, be) in ground_building_elements.iter().enumerate() {
            assert_relative_eq!(be.r_si(), results[i], max_relative = 0.05);
        }
    }

    #[rstest]
    pub fn test_h_ci_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        let temp_int_air = 20.0;
        let temp_int_surface = [19.0, 21.0, 22.0, 21.0, 19.0];
        let results = [0.7, 5.0, 2.5, 0.7, 5.0];

        for (i, be) in ground_building_elements.iter().enumerate() {
            assert_eq!(be.h_ci(temp_int_air, temp_int_surface[i]), results[i],);
        }
    }

    #[rstest]
    pub fn test_h_ri_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        for be in ground_building_elements.iter() {
            assert_eq!(be.h_ri(), 5.13, "incorrect h_ri returned");
        }
    }

    #[rstest]
    pub fn test_h_ce_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
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
    pub fn test_h_re_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        for be in ground_building_elements.iter() {
            assert_eq!(be.h_re(), 0.0, "incorrect h_re returned");
        }
    }

    #[rstest]
    pub fn test_a_sol_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        for be in ground_building_elements.iter() {
            assert_eq!(be.a_sol(), 0.0, "incorrect a_sol returned");
        }
    }

    #[rstest]
    pub fn test_therm_rad_to_sky_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        for be in ground_building_elements.iter() {
            assert_eq!(
                be.therm_rad_to_sky(),
                0.0,
                "incorrect therm_rad_to_sky returned"
            );
        }
    }

    #[rstest]
    pub fn test_h_pli_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
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
    pub fn test_k_pli_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
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

    // following test seems incomplete - Python test only seems to exercise the first building element in the list, and other elements come out with different values

    #[ignore = "this faulty test has been superseded upstream by a more complex (and working!) test"]
    #[rstest]
    pub fn test_temp_ext_for_ground(
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
    pub fn test_fabric_heat_loss_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        let expected = [30.0, 31.5, 33.25, 34.375, 30.0];
        for (i, be) in ground_building_elements.iter().enumerate() {
            assert_relative_eq!(be.fabric_heat_loss(), expected[i], max_relative = 1e-2);
        }
    }

    #[rstest]
    pub fn test_heat_capacity_for_ground(ground_building_elements: [BuildingElementGround; 5]) {
        let results = [380., 405., 425., 440., 450.];
        for (i, be) in ground_building_elements.iter().enumerate() {
            assert_eq!(
                be.heat_capacity(),
                results[i],
                "incorrect heat capacity returned"
            );
        }
    }

    #[fixture]
    pub fn transparent_building_element(
        external_conditions: Arc<ExternalConditions>,
        simulation_time: SimulationTimeIterator,
    ) -> BuildingElementTransparent {
        BuildingElementTransparent::new(
            90.,
            0.4,
            180.,
            0.75,
            0.25,
            1.,
            1.25,
            4.,
            vec![],
            None, // TODO: check if this needs updating as part of migration to 0.32
            external_conditions,
            &simulation_time,
        )
    }

    #[rstest]
    pub fn test_no_of_nodes_for_transparent(
        transparent_building_element: BuildingElementTransparent,
    ) {
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
    pub fn test_area_for_transparent(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.area(),
            5.0,
            "incorrect area returned"
        );
    }

    #[rstest]
    pub fn test_heat_flow_direction(transparent_building_element: BuildingElementTransparent) {
        // Python test uses None here, but reasonable to expecta temperature to be passed in because
        // that example only works because of the pitch set on this building element, and would fail
        // for a different value
        let stock_temperature = 20.;
        assert_eq!(
            transparent_building_element.heat_flow_direction(stock_temperature, stock_temperature),
            HeatFlowDirection::Horizontal,
            "incorrect heat flow direction returned"
        );
    }

    #[rstest]
    pub fn test_r_si_for_transparent(transparent_building_element: BuildingElementTransparent) {
        assert_relative_eq!(
            transparent_building_element.r_si(),
            0.13,
            max_relative = 1e-2
        );
    }

    #[rstest]
    pub fn test_h_ci_for_transparent(transparent_building_element: BuildingElementTransparent) {
        // Python test uses None here, but reasonable to expecta temperature to be passed in because
        // that example only works because of the pitch set on this building element, and would fail
        // for a different value
        let stock_temperature = 20.;
        assert_eq!(
            transparent_building_element.h_ci(stock_temperature, stock_temperature),
            2.5,
            "incorrect h_ci returned"
        );
    }

    #[rstest]
    pub fn test_h_ri_for_transparent(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.h_ri(),
            5.13,
            "incorrect h_ri returned"
        );
    }

    #[rstest]
    pub fn test_h_ce_for_transparent(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.h_ce(),
            20.0,
            "incorrect h_ce returned"
        );
    }

    #[rstest]
    pub fn test_h_re_for_transparent(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.h_re(),
            4.14,
            "incorrect h_re returned"
        );
    }

    #[rstest]
    pub fn test_a_sol_for_transparent(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.a_sol(),
            0.0,
            "incorrect a_sol returned"
        );
    }

    #[rstest]
    pub fn test_therm_rad_to_sky_for_transparent(
        transparent_building_element: BuildingElementTransparent,
    ) {
        assert_eq!(
            transparent_building_element.therm_rad_to_sky(),
            22.77,
            "incorrect therm_rad_to_sky returned"
        );
    }

    #[rstest]
    pub fn test_h_pli_for_transparent(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.h_pli(),
            &[2.5],
            "incorrect h_pli list returned"
        );
    }

    #[rstest]
    pub fn test_k_pli_for_transparent(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.k_pli(),
            &[0.0, 0.0],
            "non-zero k_pli list returned"
        );
    }

    #[rstest]
    pub fn test_temp_ext_for_transparent(
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

    #[rstest]
    pub fn test_fabric_heat_loss_for_transparent(
        transparent_building_element: BuildingElementTransparent,
    ) {
        assert_relative_eq!(
            transparent_building_element.fabric_heat_loss(),
            8.16,
            max_relative = 1e-2
        );
    }

    #[rstest]
    pub fn test_heat_capacity(transparent_building_element: BuildingElementTransparent) {
        assert_eq!(
            transparent_building_element.heat_capacity(),
            0.,
            "incorrect heat capacity returned"
        );
    }
}
