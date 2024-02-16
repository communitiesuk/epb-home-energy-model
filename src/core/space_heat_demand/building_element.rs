use crate::core::space_heat_demand::zone::NamedBuildingElement;
use crate::core::units::{average_monthly_to_annual, JOULES_PER_KILOJOULE};
use crate::external_conditions::ExternalConditions;
use crate::input::{BuildingElement, MassDistributionClass};
use crate::simulation_time::{SimulationTimeIteration, SimulationTimeIterator};
use std::f64::consts::PI;

impl BuildingElement {
    pub fn fabric_heat_loss(&self) -> f64 {
        match *self {
            BuildingElement::Opaque {
                u_value,
                r_c,
                area,
                pitch,
                ..
            } => {
                let u_value =
                    u_value.unwrap_or_else(|| 1. / (r_c.unwrap() + R_SE + r_si_for_pitch(pitch)));
                area * u_value
            }
            BuildingElement::AdjacentZTC { .. } => 0.0, // no heat loss to thermally conditioned zones
            BuildingElement::AdjacentZTUSimple {
                u_value,
                r_c,
                area,
                pitch,
                ..
            } => {
                let u_value =
                    u_value.unwrap_or_else(|| 1. / (r_c.unwrap() + R_SE + r_si_for_pitch(pitch)));
                area * u_value
            }
            BuildingElement::Ground { u_value, area, .. } => u_value * area,
            BuildingElement::Transparent {
                u_value,
                r_c,
                pitch,
                ..
            } => {
                // Effective window U-value includes assumed use of curtains/blinds, see
                // SAP10.2 spec, paragraph 3.2
                // TODO Confirm this is still the desired approach for SAP 11
                let r_curtains_blinds = 0.04;
                let u_value = u_value.unwrap_or_else(|| {
                    1. / ((r_c.unwrap() + r_si_for_pitch(pitch) + R_SE) + r_curtains_blinds)
                });
                area_for_building_element_input(self) * u_value
            }
        }
    }

    pub fn heat_capacity(&self) -> f64 {
        match *self {
            BuildingElement::Opaque { area, k_m, .. } => area * (k_m / JOULES_PER_KILOJOULE as f64),
            BuildingElement::AdjacentZTC { area, k_m, .. } => {
                area * (k_m / JOULES_PER_KILOJOULE as f64)
            }
            BuildingElement::AdjacentZTUSimple { area, k_m, .. } => {
                area * (k_m / JOULES_PER_KILOJOULE as f64)
            }
            BuildingElement::Ground { area, k_m, .. } => area * (k_m / JOULES_PER_KILOJOULE as f64),
            BuildingElement::Transparent { .. } => 0.0, // Set to zero as not included in heat loss calculations
        }
    }

    /// Return calculated solar gains using pitch and orientation of element
    pub fn solar_gains(
        &self,
        external_conditions: &ExternalConditions,
        simulation_time: SimulationTimeIteration,
    ) -> f64 {
        match *self {
            BuildingElement::Transparent {
                area,
                pitch,
                orientation,
                g_value,
                frame_area_fraction,
                ..
            } => {
                let (i_sol_dir, i_sol_dif, _, _) = external_conditions
                    .calculated_direct_diffuse_total_irradiance(
                        pitch,
                        orientation,
                        false,
                        &simulation_time,
                    );
                let g_value = Self::convert_g_value(g_value);

                let (f_sh_dir, f_sh_dif) =
                    shading_factors_direct_diffuse_for(self, external_conditions);
                g_value
                    * (i_sol_dif * f_sh_dif + i_sol_dir * f_sh_dir)
                    * area.expect("area expected to be available for transparent building element")
                    * (1. - frame_area_fraction)
            }
            _ => 0.,
        }
    }

    /// return g_value corrected for angle of solar radiation
    fn convert_g_value(g_value: f64) -> f64 {
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
        fw * g_value
    }
}

pub fn area_for_building_element_input(element: &BuildingElement) -> f64 {
    match *element {
        BuildingElement::Opaque { area: a, .. } => a,
        BuildingElement::Transparent {
            area: a,
            height,
            width,
            ..
        } => match a {
            Some(a) => a,
            None => height * width, // just to give some nominal value
        },
        BuildingElement::Ground { area: a, .. } => a,
        BuildingElement::AdjacentZTC { area: a, .. } => a,
        BuildingElement::AdjacentZTUSimple { area: a, .. } => a,
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

/// Calculate the number of nodes for a building element
/// (this is a shortcut to properly implementing building elements!)
pub fn number_of_building_element_nodes(element: &BuildingElement) -> u32 {
    k_pli_for(element).len() as u32
}

/// Return number of nodes excluding external and internal layers
pub fn number_of_inside_nodes(element: &BuildingElement) -> u32 {
    number_of_building_element_nodes(element) - 2
}

/// Calculate node heat capacities (k_pli)
/// according to BS EN ISO 52016-1:2017, section 6.5.7.3
pub fn k_pli_for(element: &BuildingElement) -> Vec<f64> {
    match element {
        BuildingElement::Opaque {
            mass_distribution_class: dist_class,
            k_m,
            ..
        } => match dist_class {
            MassDistributionClass::I => vec![0.0, 0.0, 0.0, 0.0, *k_m],
            MassDistributionClass::E => vec![*k_m, 0.0, 0.0, 0.0, 0.0],
            MassDistributionClass::IE => {
                let k_ie = k_m / 2.0;
                vec![k_ie, 0.0, 0.0, 0.0, k_ie]
            }
            MassDistributionClass::D => {
                let k_inner = *k_m / 4.0;
                let k_outer = *k_m / 8.0;
                vec![k_outer, k_inner, k_inner, k_inner, k_outer]
            }
            MassDistributionClass::M => vec![0.0, 0.0, *k_m, 0.0, 0.0],
        },
        BuildingElement::AdjacentZTC {
            mass_distribution_class: dist_class,
            k_m,
            ..
        } => match dist_class {
            MassDistributionClass::I => vec![0.0, 0.0, 0.0, 0.0, *k_m],
            MassDistributionClass::E => vec![*k_m, 0.0, 0.0, 0.0, 0.0],
            MassDistributionClass::IE => {
                let k_ie = *k_m / 2.0;
                vec![k_ie, 0.0, 0.0, 0.0, k_ie]
            }
            MassDistributionClass::D => {
                let k_inner = *k_m / 4.0;
                let k_outer = *k_m / 8.0;
                vec![k_outer, k_inner, k_inner, k_inner, k_outer]
            }
            MassDistributionClass::M => vec![0.0, 0.0, *k_m, 0.0, 0.0],
        },
        BuildingElement::AdjacentZTUSimple {
            mass_distribution_class: dist_class,
            k_m,
            ..
        } => match dist_class {
            MassDistributionClass::I => vec![0.0, 0.0, 0.0, 0.0, *k_m],
            MassDistributionClass::E => vec![*k_m, 0.0, 0.0, 0.0, 0.0],
            MassDistributionClass::IE => {
                let k_ie = k_m / 2.0;
                vec![k_ie, 0.0, 0.0, 0.0, k_ie]
            }
            MassDistributionClass::D => {
                let k_inner = *k_m / 4.0;
                let k_outer = *k_m / 8.0;
                vec![k_outer, k_inner, k_inner, k_inner, k_outer]
            }
            MassDistributionClass::M => vec![0.0, 0.0, *k_m, 0.0, 0.0],
        },
        BuildingElement::Ground {
            mass_distribution_class: dist_class,
            k_m,
            ..
        } => {
            let k_gr = K_GR_FOR_GROUND;
            match dist_class {
                MassDistributionClass::I => vec![0.0, k_gr, 0.0, 0.0, *k_m],
                MassDistributionClass::E => vec![0.0, k_gr, *k_m, 0.0, 0.0],
                MassDistributionClass::IE => {
                    let k_ie = *k_m / 2.0;
                    vec![0.0, k_gr, k_ie, 0.0, k_ie]
                }
                MassDistributionClass::D => {
                    let k_inner = *k_m / 2.0;
                    let k_outer = *k_m / 4.0;
                    vec![0.0, k_gr, k_outer, k_inner, k_outer]
                }
                MassDistributionClass::M => vec![0.0, k_gr, 0.0, *k_m, 0.0],
            }
        }
        &BuildingElement::Transparent { .. } => vec![0.0, 0.0],
    }
}

/// Calculate node conductances (h_pli)
/// according to BS EN ISO 52016-1:2017, section 6.5.7.3
pub fn h_pli_for(element: &BuildingElement) -> Vec<f64> {
    match element {
        BuildingElement::Opaque { .. } => {
            let r_c = r_c_for(element);
            let h_outer = 6.0 / r_c;
            let h_inner = 3.0 / r_c;

            vec![h_outer, h_inner, h_inner, h_outer]
        }
        BuildingElement::AdjacentZTC { .. } => {
            let r_c = r_c_for(element);
            let h_outer = 6.0 / r_c;
            let h_inner = 3.0 / r_c;

            vec![h_outer, h_inner, h_inner, h_outer]
        }
        BuildingElement::AdjacentZTUSimple { .. } => {
            let r_c = r_c_for(element);
            let h_outer = 6.0 / r_c;
            let h_inner = 3.0 / r_c;

            vec![h_outer, h_inner, h_inner, h_outer]
        }
        BuildingElement::Ground { u_value, .. } => {
            let r_c = 1.0 / u_value;
            let r_gr = R_GR_FOR_GROUND;

            vec![
                2.0 / r_gr,
                1.0 / (r_c / 4.0 + r_gr / 2.0),
                2.0 / r_c,
                4.0 / r_c,
            ]
        }
        BuildingElement::Transparent { .. } => vec![1.0 / r_c_for(element)],
    }
}

// this is based on __init_resistance_or_uvalue in project.py
fn r_c_for(element: &BuildingElement) -> f64 {
    match element {
        BuildingElement::Opaque {
            r_c,
            u_value,
            pitch,
            ..
        } => match (r_c, u_value) {
            (Some(r_c), _) => *r_c,
            (None, Some(u_value)) => convert_uvalue_to_resistance(*u_value, *pitch),
            _ => panic!("either r_c or uvalue needed to be provided for this element and this should have been caught at input validation")
        },
        BuildingElement::Transparent {
            r_c,
            u_value,
            pitch,
            ..
        } => match (r_c, u_value) {
            (Some(r_c), _) => *r_c,
            (None, Some(u_value)) => convert_uvalue_to_resistance(*u_value, *pitch),
            _ => panic!("either r_c or uvalue needed to be provided for this element and this should have been caught at input validation")
        },
        BuildingElement::Ground {
            u_value,
            pitch,
            ..
        } => convert_uvalue_to_resistance(*u_value, *pitch),
        BuildingElement::AdjacentZTC {
            r_c,
            u_value,
            pitch,
            ..
        } => match (r_c, u_value) {
            (Some(r_c), _) => *r_c,
            (None, Some(u_value)) => convert_uvalue_to_resistance(*u_value, *pitch),
            _ => panic!("either r_c or uvalue needed to be provided for this element and this should have been caught at input validation")
        },
        BuildingElement::AdjacentZTUSimple {
            r_c,
            u_value,
            pitch,
            ..
        } => match (r_c, u_value) {
            (Some(r_c), _) => *r_c,
            (None, Some(u_value)) => convert_uvalue_to_resistance(*u_value, *pitch),
            _ => panic!("either r_c or uvalue needed to be provided for this element and this should have been caught at input validation")
        },
    }
}

fn convert_uvalue_to_resistance(u_value: f64, pitch: f64) -> f64 {
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

/// Return external convective heat transfer coefficient, in W / (m2.K)
pub fn h_ce_for(element: &BuildingElement) -> f64 {
    match element {
        BuildingElement::Opaque { .. } => H_CE,
        BuildingElement::AdjacentZTC { .. } => {
            // Element is adjacent to another building / thermally conditioned zone
            // therefore according to BS EN ISO 52016-1:2017, section 6.5.6.3.6,
            // external heat transfer coefficients are zero
            0.0
        }
        BuildingElement::AdjacentZTUSimple { r_u, .. } => {
            // Add an additional thermal resistance to the outside of the wall and
            // incorporate this in the values for the external surface heat transfer
            // coefficient.
            // As this is an adjusted figure in this class, and the split between
            // h_ce and h_re does not affect the calculation results, assign entire
            // effective surface heat transfer to h_ce and set h_re to zero.
            1.0 / ((1.0 / (H_CE + H_RE)) + r_u)
        }
        BuildingElement::Ground { u_value, r_f, .. } => {
            // Calculate thermal resistance of virtual layer using BS EN ISO 13370:2017 Equation (F1)
            let r_vi = (1.0 / u_value) - R_SI_FOR_GROUND - r_f - R_GR_FOR_GROUND; // in m2.K/W
            assert!(r_vi > 0.0, "r_vi should be greater than zero. check u-value and r_f inputs for floors - this should be checked at input validation");

            1.0 / r_vi
        }
        BuildingElement::Transparent { .. } => H_CE,
    }
}

/// Return external radiative heat transfer coefficient, in W / (m2.K)
pub fn h_re_for(element: &BuildingElement) -> f64 {
    match element {
        BuildingElement::Opaque { .. } => H_RE,
        BuildingElement::AdjacentZTC { .. } => {
            // Element is adjacent to another building / thermally conditioned zone
            // therefore according to BS EN ISO 52016-1:2017, section 6.5.6.3.6,
            // external heat transfer coefficients are zero
            0.0
        }
        BuildingElement::AdjacentZTUSimple { .. } => {
            // As this is an adjusted figure in this class (=variant of building element),
            // and the split between h_ce and h_re does not affect the calculation results,
            // assign entire effective surface heat transfer to h_ce and set h_re to zero.
            0.0
        }
        BuildingElement::Ground { .. } => 0.0,
        BuildingElement::Transparent { .. } => H_RE,
    }
}

pub fn i_sol_dir_dif_for(
    element: &BuildingElement,
    external_conditions: &ExternalConditions,
    simulation_time: &SimulationTimeIteration,
) -> (f64, f64) {
    match element {
        BuildingElement::Opaque {
            pitch, orientation, ..
        } => {
            let (i_sol_dir, i_sol_dif, _, _) = external_conditions
                .calculated_direct_diffuse_total_irradiance(
                    *pitch,
                    *orientation,
                    false,
                    simulation_time,
                );

            (i_sol_dir, i_sol_dif)
        }
        _ => (0.0, 0.0),
    }
}

pub fn shading_factors_direct_diffuse_for(
    element: &BuildingElement,
    _external_conditions: &ExternalConditions,
) -> (f64, f64) {
    match element {
        BuildingElement::Opaque { .. } => {
            // external_conditions --- argh uses currently broken method
            // TODO: fix this!!
            (1.0, 1.0) // fixme
        }
        BuildingElement::Transparent { .. } => {
            (1.0, 1.0) // fixme
        }
        _ => (1.0, 1.0),
    }
}

/// Return the temperature on the other side of the building element
pub fn temp_ext_for(
    element: &BuildingElement,
    external_conditions: &ExternalConditions,
    simulation_time: &SimulationTimeIteration,
) -> f64 {
    match element {
        BuildingElement::Ground {
            perimeter,
            psi_wall_floor_junc,
            u_value,
            area,
            h_pi,
            h_pe,
            ..
        } => {
            let temp_ext_annual = external_conditions
                .air_temp_annual()
                .expect("no annual air temp available");
            let temp_ext_month = external_conditions
                .air_temp_monthly(simulation_time.current_month_start_end_hours());

            let current_month = simulation_time.current_month().unwrap_or(0);
            let temp_int_month = TEMP_INT_MONTHLY_FOR_GROUND[current_month as usize];

            let temp_int_annual = average_monthly_to_annual(TEMP_INT_MONTHLY_FOR_GROUND);

            // BS EN ISO 13370:2017 Eqn C.4
            let heat_flow_month = u_value * area * (temp_int_annual - temp_ext_annual)
                + perimeter * psi_wall_floor_junc * (temp_int_month - temp_ext_month)
                - h_pi * (temp_int_annual - temp_int_month)
                + h_pe * (temp_ext_annual - temp_ext_month);

            // BS EN ISO 13370:2017 Eqn F.2
            temp_int_month
                - (heat_flow_month
                    - (perimeter * psi_wall_floor_junc * (temp_int_annual - temp_ext_annual)))
                    / (area * u_value)
        }
        _ => external_conditions.air_temp_for_timestep_idx(simulation_time.index),
    }
}

pub fn a_sol_for(element: &BuildingElement) -> f64 {
    match element {
        BuildingElement::Opaque { a_sol, .. } => *a_sol,
        _ => 0.0,
    }
}

fn f_sky_for(element: &BuildingElement) -> f64 {
    match element {
        BuildingElement::Opaque { pitch, .. } => sky_view_factor(pitch),
        BuildingElement::Transparent { pitch, .. } => sky_view_factor(pitch),
        _ => 0.0,
    }
}

/// Calculate longwave sky view factor from pitch in degrees
fn sky_view_factor(pitch: &f64) -> f64 {
    // # TODO account for shading
    // # TODO check longwave is correct
    let pitch_rads = pitch * PI / 180.0;

    0.5 * (1.0 + pitch_rads.cos())
}

pub fn therm_rad_to_sky_for(element: &BuildingElement) -> f64 {
    f_sky_for(element) * h_re_for(element) * TEMP_DIFF_SKY
}

/// Return internal convective heat transfer coefficient, in W / (m2.K)
pub fn h_ci_for(element: &BuildingElement, temp_int_air: f64, temp_int_surface: f64) -> f64 {
    match heat_flow_direction_for(element, temp_int_air, temp_int_surface) {
        HeatFlowDirection::Horizontal => H_CI_HORIZONTAL,
        HeatFlowDirection::Upwards => H_CI_UPWARDS,
        HeatFlowDirection::Downwards => H_CI_DOWNWARDS,
    }
}

/// Determine direction of heat flow for a surface
fn heat_flow_direction_for(
    element: &BuildingElement,
    temp_int_air: f64,
    temp_int_surface: f64,
) -> HeatFlowDirection {
    let pitch = pitch_for(element);
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

fn pitch_for(element: &BuildingElement) -> f64 {
    match element {
        BuildingElement::Opaque { pitch, .. } => *pitch,
        BuildingElement::AdjacentZTC { pitch, .. } => *pitch,
        BuildingElement::AdjacentZTUSimple { pitch, .. } => *pitch,
        BuildingElement::Ground { pitch, .. } => *pitch,
        BuildingElement::Transparent { pitch, .. } => *pitch,
    }
}

/// Return internal radiative heat transfer coefficient, in W / (m2.K)
pub fn h_ri_for(_: &BuildingElement) -> f64 {
    H_RI
}

pub fn mid_height_for(element: &BuildingElement) -> Option<f64> {
    match element {
        BuildingElement::Transparent {
            base_height,
            height,
            ..
        } => Some(base_height + height / 2.0),
        _ => None,
    }
}

pub fn orientation_for(element: &BuildingElement) -> Option<f64> {
    match element {
        BuildingElement::Transparent { orientation, .. } => Some(*orientation),
        _ => None,
    }
}

pub fn projected_height_for_transparent_element(element: &BuildingElement) -> Option<f64> {
    match element {
        BuildingElement::Transparent { pitch, height, .. } => {
            Some(projected_height(*pitch, *height))
        }
        _ => None,
    }
}

#[derive(Debug, PartialEq)]
enum HeatFlowDirection {
    Horizontal,
    Upwards,
    Downwards,
}

pub fn element_from_named(named_element: &NamedBuildingElement) -> &BuildingElement {
    let NamedBuildingElement { element, .. } = named_element;
    element
}

pub fn cloned_element_from_named(named_element: &NamedBuildingElement) -> BuildingElement {
    let NamedBuildingElement { element, .. } = named_element;
    (*element).clone()
}

// Thermal properties of ground from BS EN ISO 13370:2017 Table 7
// Use values for clay or silt (same as BR 443 and SAP 10)
const THERMAL_CONDUCTIVITY_OF_GROUND: f64 = 1.5;
// in W/(m.K)
const HEAT_CAPACITY_PER_VOLUME_OF_GROUND: f64 = 3000000.0;
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
const PITCH_LIMIT_HORIZ_CEILING: f64 = 60.0;
const PITCH_LIMIT_HORIZ_FLOOR: f64 = 120.0;

#[cfg(test)]
mod test {
    use super::*;
    use crate::external_conditions::DaylightSavingsConfig;
    use crate::simulation_time::SimulationTime;
    use rstest::*;

    #[fixture]
    pub fn simulation_time() -> SimulationTimeIterator {
        SimulationTime::new(0.0, 4.0, 1.0).iter()
    }

    #[fixture]
    pub fn external_conditions(simulation_time: SimulationTimeIterator) -> ExternalConditions {
        ExternalConditions::new(
            &simulation_time,
            vec![0.0, 5.0, 10.0, 15.0],
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
            DaylightSavingsConfig::NotApplicable,
            false,
            false,
            vec![],
        )
    }

    #[fixture]
    pub fn be_i() -> BuildingElement {
        BuildingElement::Opaque {
            area: 20.0,
            pitch: 180.0,
            a_sol: 0.6,
            u_value: None,
            r_c: Some(0.25),
            k_m: 19000.0,
            mass_distribution_class: MassDistributionClass::I,
            is_external_door: None,
            orientation: 0.0,
            base_height: 0.0,
            height: 2.0,
            width: 10.0,
            h_re: None,
            h_ci: None,
            h_ri: None,
            h_ce: None,
        }
    }

    #[fixture]
    pub fn be_e() -> BuildingElement {
        BuildingElement::Opaque {
            area: 22.5,
            pitch: 135.0,
            a_sol: 0.61,
            u_value: None,
            r_c: Some(0.5),
            k_m: 18000.0,
            mass_distribution_class: MassDistributionClass::E,
            is_external_door: None,
            orientation: 180.0,
            base_height: 0.0,
            height: 2.25,
            width: 10.0,
            h_re: None,
            h_ci: None,
            h_ri: None,
            h_ce: None,
        }
    }

    #[fixture]
    pub fn be_ie() -> BuildingElement {
        BuildingElement::Opaque {
            area: 25.0,
            pitch: 90.0,
            a_sol: 0.62,
            u_value: None,
            r_c: Some(0.75),
            k_m: 17000.0,
            mass_distribution_class: MassDistributionClass::IE,
            is_external_door: None,
            orientation: 90.0,
            base_height: 0.0,
            height: 2.5,
            width: 10.0,
            h_re: None,
            h_ci: None,
            h_ri: None,
            h_ce: None,
        }
    }

    #[fixture]
    pub fn be_d() -> BuildingElement {
        BuildingElement::Opaque {
            area: 27.5,
            pitch: 45.0,
            a_sol: 0.63,
            u_value: None,
            r_c: Some(0.8),
            k_m: 16000.0,
            mass_distribution_class: MassDistributionClass::D,
            is_external_door: None,
            orientation: -90.0,
            base_height: 0.0,
            height: 2.75,
            width: 10.0,
            h_re: None,
            h_ci: None,
            h_ri: None,
            h_ce: None,
        }
    }

    #[fixture]
    pub fn be_m() -> BuildingElement {
        BuildingElement::Opaque {
            area: 30.0,
            pitch: 0.0,
            a_sol: 0.64,
            u_value: None,
            r_c: Some(0.4),
            k_m: 15000.0,
            mass_distribution_class: MassDistributionClass::M,
            is_external_door: None,
            orientation: 0.0,
            base_height: 0.0,
            height: 3.0,
            width: 10.0,
            h_re: None,
            h_ci: None,
            h_ri: None,
            h_ce: None,
        }
    }

    #[fixture]
    pub fn building_elements(
        be_i: BuildingElement,
        be_e: BuildingElement,
        be_ie: BuildingElement,
        be_d: BuildingElement,
        be_m: BuildingElement,
    ) -> [BuildingElement; 5] {
        [be_i, be_e, be_ie, be_d, be_m]
    }

    fn round_by_precision(src: f64, precision: f64) -> f64 {
        (precision * src).round() / precision
    }

    #[rstest]
    pub fn should_have_correct_no_of_nodes(building_elements: [BuildingElement; 5]) {
        for be in building_elements.iter() {
            assert_eq!(
                number_of_building_element_nodes(be),
                5,
                "incorrect number of nodes"
            );
            assert_eq!(
                number_of_inside_nodes(be),
                3,
                "incorrect number of inside nodes"
            );
        }
    }

    #[rstest]
    pub fn should_have_correct_areas(building_elements: [BuildingElement; 5]) {
        // Define increment between test cases
        let area_inc = 2.5;
        for (i, be) in building_elements.iter().enumerate() {
            assert_eq!(
                area_for_building_element_input(be),
                20.0 + i as f64 * area_inc,
                "incorrect area returned"
            );
        }
    }

    #[rstest]
    pub fn should_have_correct_heat_flow_direction(building_elements: [BuildingElement; 5]) {
        let temp_int_air = 20.0;
        let temp_int_surface = [19.0, 21.0, 22.0, 21.0, 19.0];
        let results = [
            HeatFlowDirection::Downwards,
            HeatFlowDirection::Upwards,
            HeatFlowDirection::Horizontal,
            HeatFlowDirection::Downwards,
            HeatFlowDirection::Upwards,
        ];
        for (i, be) in building_elements.iter().enumerate() {
            assert_eq!(
                heat_flow_direction_for(be, temp_int_air, temp_int_surface[i]),
                results[i],
                "incorrect heat flow direction returned"
            )
        }
    }

    #[rstest]
    pub fn should_have_correct_r_si(building_elements: [BuildingElement; 5]) {
        let results = [0.17, 0.17, 0.13, 0.10, 0.10];

        for (i, be) in building_elements.iter().enumerate() {
            assert_eq!(
                round_by_precision(r_si_for_pitch(pitch_for(be)), 1e2),
                round_by_precision(results[i], 1e2),
                "incorrect r_si returned"
            )
        }
    }

    #[rstest]
    pub fn should_have_correct_h_ci(building_elements: [BuildingElement; 5]) {
        let temp_int_air = 20.0;
        let temp_int_surface = [19.0, 21.0, 22.0, 21.0, 19.0];
        let results = [0.7, 5.0, 2.5, 0.7, 5.0];

        for (i, be) in building_elements.iter().enumerate() {
            assert_eq!(
                round_by_precision(h_ci_for(be, temp_int_air, temp_int_surface[i]), 1e1),
                round_by_precision(results[i], 1e1),
                "incorrect h_ci returned"
            );
        }
    }

    // incomplete - lots more tests to migrate here!
}
