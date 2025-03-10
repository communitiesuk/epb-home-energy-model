// This module provides objects to represent Infiltration and Ventilation.
// The calculations are based on Method 1 of BS EN 16798-7.

use crate::compare_floats::max_of_2;
use crate::core::controls::time_control::{Control, ControlBehaviour};
use crate::core::ductwork::Ductwork;
use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::core::material_properties::AIR;
use crate::core::solvers::fsolve;
use crate::core::space_heat_demand::building_element::{pitch_class, HeatFlowDirection};
use crate::core::units::{
    celsius_to_kelvin, LITRES_PER_CUBIC_METRE, MILLIMETRES_IN_METRE, SECONDS_PER_HOUR,
    WATTS_PER_KILOWATT,
};
use crate::corpus::{CompletedVentilationLeaks, Controls, ReportingFlag};
use crate::errors::NotImplementedError;
use crate::input::{
    init_orientation, BuildingElement, CombustionAirSupplySituation, CombustionApplianceType,
    CombustionFuelType, DuctShape, FlueGasExhaustSituation,
    InfiltrationVentilation as InfiltrationVentilationInput, SupplyAirFlowRateControlType,
    SupplyAirTemperatureControlType, TerrainClass, VentType, VentilationShieldClass,
    WindowPart as WindowPartInput, ZoneDictionary,
};
use crate::simulation_time::SimulationTimeIteration;
use anyhow::{anyhow, bail, Error};
use argmin::{
    core::{CostFunction, Executor},
    solver::brent::BrentRoot,
};
use indexmap::IndexMap;
use itertools::Itertools;
use parking_lot::RwLock;
use std::sync::Arc;
use thiserror::Error;

fn p_a_ref() -> f64 {
    AIR.density_kg_per_m3()
}

//referenced for unimplementable LOAD arm further down in module - remove following directive when this is implemented
#[allow(dead_code)]
fn c_a() -> f64 {
    AIR.specific_heat_capacity_kwh()
}

// (Default values from BS EN 16798-7, Table 11)
// Coefficient to take into account stack effect in airing calculation in (m/s)/(m*K)
const _C_STACK: f64 = 0.0035;
// Coefficient to take into account wind speed in airing calculation in 1/(m/s)
const _C_WND: f64 = 0.001;
// Gravitational constant in m/s2
const G: f64 = 9.81;
//Room temperature in degrees K
const T_E_REF: f64 = 293.15;
//Absolute zero in degrees K
//referenced for unimplementable LOAD arm further down in module - remove following directive when this is implemented
#[allow(dead_code)]
const T_0_ABS: f64 = 273.15;

// In Python this is defined in InfiltrationVentilation.calculate_internal_reference_pressure
const INTERVAL_EXPANSION_LIST: [f64; 9] = [1., 5., 10., 15., 20., 40., 50., 100., 200.];

/// Calculate pressure difference between the exterior and the interior of the dwelling
/// for a flow path (at it's elevavation above the vent zone floor)
/// Arguments:
/// * `h_path` - height of air flow path (m)
/// * `c_p_path` - wind pressure coefficient
/// * `u_site` - wind velocity at zone level (m/s)
/// * `t_e` - external air temperature (K)
/// * `t_z` - thermal zone air temperature (K)
/// * `p_z_ref` - internal reference pressure (Pa)
fn calculate_pressure_difference_at_an_airflow_path(
    h_path: f64,
    c_p_path: f64,
    u_site: f64,
    t_e: f64,
    t_z: f64,
    p_z_ref: f64,
) -> f64 {
    let p_e_path = p_a_ref() * T_E_REF / t_e * (0.5 * c_p_path * u_site.powi(2) - h_path * G); //(5)
    let p_z_path = p_z_ref - p_a_ref() * h_path * G * T_E_REF / t_z; //(6)
    p_e_path - p_z_path //(4)
}

/// Convert infiltration rate from ach to m^3/s
fn _air_change_rate_to_flow_rate(air_change_rate: f64, zone_volume: f64) -> f64 {
    air_change_rate * zone_volume / SECONDS_PER_HOUR as f64
}

/// Table B.3 fuel flow factors
/// Arguments:
/// fuel_type -- options are 'wood', 'gas', 'oil' or 'coal'.
/// appliance_type -- options are 'open_fireplace', 'closed_with_fan',
/// open_gas_flue_balancer', 'open_gas_kitchen_stove', 'open_gas_fire' or 'closed_fire'
fn get_fuel_flow_factor(
    fuel_type: CombustionFuelType,
    appliance_type: CombustionApplianceType,
) -> f64 {
    match (fuel_type, appliance_type) {
        (CombustionFuelType::Wood, CombustionApplianceType::OpenFireplace) => 2.8,
        (CombustionFuelType::Gas, CombustionApplianceType::ClosedWithFan) => 0.38,
        (CombustionFuelType::Gas, CombustionApplianceType::OpenGasFlueBalancer) => 0.78,
        (CombustionFuelType::Gas, CombustionApplianceType::OpenGasKitchenStove) => 3.35,
        (CombustionFuelType::Gas, CombustionApplianceType::OpenGasFire) => 3.35,
        (CombustionFuelType::Oil, CombustionApplianceType::ClosedFire) => 0.32,
        (CombustionFuelType::Coal, CombustionApplianceType::ClosedFire) => 0.52,
        (_, _) => panic!("Invalid combination of fuel and appliance types ({fuel_type:?} and {appliance_type:?}."),
    }
}

/// Interpreted from Table B.2 from BS EN 16798-7, get the appliance system factor
/// for a combustion appliance.
/// Arguments:
/// supply_situation -- Combustion air supply situation: 'room_air' or 'outside'
/// exhaust_situation -- flue gas exhaust situation: 'into_room', 'into_separate_duct' or 'into_mech_vent'
fn get_appliance_system_factor(
    supply_situation: CombustionAirSupplySituation,
    exhaust_situation: FlueGasExhaustSituation,
) -> f64 {
    match (supply_situation, exhaust_situation) {
        (CombustionAirSupplySituation::Outside, _) => 0.,
        (CombustionAirSupplySituation::RoomAir, FlueGasExhaustSituation::IntoRoom) => 0.,
        (CombustionAirSupplySituation::RoomAir, FlueGasExhaustSituation::IntoSeparateDuct) => 1.,
        (CombustionAirSupplySituation::RoomAir, FlueGasExhaustSituation::IntoMechVent) => {
            panic!("Invalid combination of supply situation ({supply_situation:?}) and exhaust situation ({exhaust_situation:?})")
        }
    }
}

/// Adjust air density for altitude above sea level.
/// Arguments:
/// h_alt -- altitude above sea level (m)
fn adjust_air_density_for_altitude(h_alt: f64) -> f64 {
    p_a_ref() * (1. - ((0.00651 * h_alt) / 293.)).powf(4.255)
}

/// Recalculate air density based on the current temperature
/// Arguments:
/// temperature -- temperature to adjust (K)
/// air_density_adjusted_for_alt - The air density after adjusting for altitude (Kg/m3)
fn air_density_at_temp(temperature: f64, air_density_adjusted_for_alt: f64) -> f64 {
    T_E_REF / temperature * air_density_adjusted_for_alt
}

/// Converts volume air flow rate (qv) to mass air flow rate (qm).
/// (Equations 65 & 66 from BS EN 16798-7)
/// Arguments:
/// qv_in -- volume flow rate of air entering the dwelling
/// qv_out -- volume flow rate of air leaving the dwelling
/// T_e -- External air temperature (K)
/// T_e -- Thermal zone air temperature (K)
/// p_a_alt -- The air density after adjusting for altitude (Kg/m3)
fn convert_to_mass_air_flow_rate(
    qv_in: f64,
    qv_out: f64,
    t_e: f64,
    t_z: f64,
    p_a_alt: f64,
) -> (f64, f64) {
    let qm_in = convert_volume_flow_rate_to_mass_flow_rate(qv_in, t_e, p_a_alt);
    let qm_out = convert_volume_flow_rate_to_mass_flow_rate(qv_out, t_z, p_a_alt);
    (qm_in, qm_out)
}

/// Convert volume flow rate in m3/hr to mass flow rate in kg/hr, at temperature in Kelvin
/// Arguments:
/// qv -- volume flow rate (m3/h)
/// temperature -- air temperature (K)
/// p_a_alt -- The air density after adjusting for altitude (Kg/m3)
fn convert_volume_flow_rate_to_mass_flow_rate(qv: f64, temperature: f64, p_a_alt: f64) -> f64 {
    qv * air_density_at_temp(temperature, p_a_alt)
}

/// Convert mass flow rate in kg/hr to volume flow rate in m3/hr, at temperature in Kelvin
/// Arguments:
/// qm -- mass flow rate (Kg/h)
/// temperature -- air temperature (K)
/// p_a_alt -- The air density after adjusting for altitude (Kg/m3)
fn convert_mass_flow_rate_to_volume_flow_rate(qm: f64, temperature: f64, p_a_alt: f64) -> f64 {
    qm / air_density_at_temp(temperature, p_a_alt)
}

/// Retrieves the roughness parameters and calculates the roughness coefficient (CR)
///     based on the terrain type and height of airflow path.
///
///     Args:
///     * `terrain_class` - The terrain type ('OpenWater', 'OpenField', 'Suburban', 'Urban').
///     * `relative_airflow_path_height` - Height of airflow path relative to the ground (m).
///
///    Returns:
///        float: Calculated roughness coefficient CR.
fn ter_class_to_roughness_coeff(terrain: &TerrainClass, relative_airflow_path_height: f64) -> f64 {
    let (kr, z0, zmin) = match terrain {
        TerrainClass::OpenWater => (0.17, 0.01, 2.),
        TerrainClass::OpenField => (0.19, 0.05, 4.),
        TerrainClass::Suburban => (0.22, 0.3, 8.),
        TerrainClass::Urban => (0.24, 1.0, 16.),
    };

    // ensure z is at least zmin
    let z = relative_airflow_path_height.max(zmin);

    // calculate the roughness coefficient
    kr * (z / z0).ln()
}

/// Meteorological wind speed at 10 m corrected to reference wind speed at zone level of the dwelling
/// Arguments:
/// C_rgh_site -- roughness coefficient at building site
/// u_10 -- wind velocity at 10m (m/s)
/// C_top_site -- topography coefficient at building site
/// C_rgh_met -- roughness coefficient at 10m depending on meteorological station
/// C_top_met -- topography coefficient at building height depending on meteorological station
fn wind_speed_at_zone_level(
    c_rgh_site: f64,
    u_10: f64,
    c_top_site: Option<f64>,
    c_rgh_met: Option<f64>,
    c_top_met: Option<f64>,
) -> f64 {
    let (c_top_site, c_rgh_met, c_top_met) = (
        c_top_site.unwrap_or(1.),
        c_rgh_met.unwrap_or(1.),
        c_top_met.unwrap_or(1.),
    );

    ((c_rgh_site * c_top_site) / (c_rgh_met * c_top_met)) * u_10
}

/// Determine difference between two bearings, taking shortest route around circle
fn orientation_difference(orientation1: f64, orientation2: f64) -> f64 {
    if !(0. ..=360.).contains(&orientation1) || !(0. ..=360.).contains(&orientation2) {
        panic!("Orientation values must be between 0 and 360 degrees"); // panicking here, but we should be able to previously enforce this constraint
    }
    let op_rel_orientation = (orientation1 - orientation2).abs();

    if op_rel_orientation > 180. {
        360. - op_rel_orientation
    } else {
        op_rel_orientation
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum FacadeDirection {
    Roof,
    Roof10,
    Roof10_30,
    Roof30,
    Windward,
    Leeward,
    Neither,
}

/// Gets direction of the facade from pitch and orientation
/// Arguments:
/// f_cross -- boolean, dependant on if cross ventilation is possible or not
/// orientation -- orientation of the facade (degrees)
/// pitch -- pitch of the facade (degrees)
/// wind_direction -- direction the wind is blowing (degrees)
fn get_facade_direction(
    f_cross: bool,
    orientation: f64,
    pitch: f64,
    wind_direction: f64,
) -> FacadeDirection {
    if f_cross {
        if pitch < 10. {
            FacadeDirection::Roof10
        } else if pitch <= 30. {
            FacadeDirection::Roof10_30
        } else if pitch < 60. {
            FacadeDirection::Roof30
        } else {
            let orientation_diff = orientation_difference(orientation, wind_direction);
            if orientation_diff <= 60. {
                FacadeDirection::Windward
            } else if orientation_diff < 120. {
                FacadeDirection::Neither
            } else {
                FacadeDirection::Leeward
            }
        }
    } else if pitch < 60. {
        FacadeDirection::Roof
    } else {
        let orientation_diff = orientation_difference(orientation, wind_direction);
        if orientation_diff <= 60. {
            FacadeDirection::Windward
        } else if orientation_diff < 120. {
            FacadeDirection::Neither
        } else {
            FacadeDirection::Leeward
        }
    }
}

// we split the python get_c_p_path method into two methods below:
fn get_c_p_path_from_pitch_and_orientation(
    f_cross: bool,
    shield_class: VentilationShieldClass,
    relative_airflow_path_height: f64,
    wind_direction: f64,
    orientation: f64,
    pitch: f64,
) -> f64 {
    let facade_direction = get_facade_direction(f_cross, orientation, pitch, wind_direction);
    get_c_p_path(
        f_cross,
        shield_class,
        relative_airflow_path_height,
        facade_direction,
    )
}

/// Interpreted from Table B.7 for determining dimensionless wind pressure coefficients
/// Arguments:
/// * `f_cross` - boolean, dependant on if cross ventilation is possible or not
/// * `shield_class` - indicates exposure to wind
/// * `relative_airflow_path_height` - height of air flow path relative to ground (m)
/// * `h_path` - height of flow path (m)
/// * `wind_direction` - direction the wind is blowing (degrees)
/// * `orientation` - orientation of the facade (degrees)
/// * `pitch` - pitch of the facade (degrees)
/// * `facade_direction` - direction of the facade (from get_facade_direction or manual entry)
fn get_c_p_path(
    f_cross: bool,
    shield_class: VentilationShieldClass,
    relative_airflow_path_height: f64,
    facade_direction: FacadeDirection,
) -> f64 {
    if f_cross {
        if relative_airflow_path_height < 15. {
            match shield_class {
                VentilationShieldClass::Open => match facade_direction {
                    FacadeDirection::Windward => 0.50,
                    FacadeDirection::Leeward => -0.70,
                    FacadeDirection::Neither => 0.0,
                    FacadeDirection::Roof10 => -0.70,
                    FacadeDirection::Roof10_30 => -0.60,
                    FacadeDirection::Roof30 => -0.20,
                    _ => panic!("Invalid combination of shield_class and facade_direction"),
                },
                VentilationShieldClass::Normal => match facade_direction {
                    FacadeDirection::Windward => 0.25,
                    FacadeDirection::Leeward => -0.50,
                    FacadeDirection::Neither => 0.0,
                    FacadeDirection::Roof10 => -0.60,
                    FacadeDirection::Roof10_30 => -0.50,
                    FacadeDirection::Roof30 => -0.20,
                    _ => panic!("Invalid combination of shield_class and facade_direction"),
                },
                VentilationShieldClass::Shielded => match facade_direction {
                    FacadeDirection::Windward => 0.05,
                    FacadeDirection::Leeward => -0.30,
                    FacadeDirection::Neither => 0.0,
                    FacadeDirection::Roof10 => -0.50,
                    FacadeDirection::Roof10_30 => -0.40,
                    FacadeDirection::Roof30 => -0.20,
                    _ => panic!("Invalid combination of shield_class and facade_direction"),
                },
            }
        } else if (15. ..50.).contains(&relative_airflow_path_height) {
            match shield_class {
                VentilationShieldClass::Open => match facade_direction {
                    FacadeDirection::Windward => 0.65,
                    FacadeDirection::Leeward => -0.70,
                    FacadeDirection::Neither => 0.0,
                    FacadeDirection::Roof10 => -0.70,
                    FacadeDirection::Roof10_30 => -0.60,
                    FacadeDirection::Roof30 => -0.20,
                    _ => panic!("Invalid combination of shield_class and facade_direction"),
                },
                VentilationShieldClass::Normal => match facade_direction {
                    FacadeDirection::Windward => 0.45,
                    FacadeDirection::Leeward => -0.50,
                    FacadeDirection::Neither => 0.0,
                    FacadeDirection::Roof10 => -0.60,
                    FacadeDirection::Roof10_30 => -0.50,
                    FacadeDirection::Roof30 => -0.20,
                    _ => panic!("Invalid combination of shield_class and facade_direction"),
                },
                VentilationShieldClass::Shielded => match facade_direction {
                    FacadeDirection::Windward => 0.25,
                    FacadeDirection::Leeward => -0.30,
                    FacadeDirection::Neither => 0.0,
                    FacadeDirection::Roof10 => -0.50,
                    FacadeDirection::Roof10_30 => -0.40,
                    FacadeDirection::Roof30 => -0.20,
                    _ => panic!("Invalid combination of shield_class and facade_direction"),
                },
            }
        } else {
            // In python this is an elif h_path >= 50.
            match shield_class {
                VentilationShieldClass::Open => match facade_direction {
                    FacadeDirection::Windward => 0.80,
                    FacadeDirection::Leeward => -0.70,
                    FacadeDirection::Neither => 0.0,
                    FacadeDirection::Roof10 => -0.70,
                    FacadeDirection::Roof10_30 => -0.60,
                    FacadeDirection::Roof30 => -0.20,
                    _ => panic!("Invalid combination of shield_class and facade_direction"),
                },
                _ => panic!("Invalid combination of shield_class and facade_direction"),
            }
        }
    } else {
        match facade_direction {
            FacadeDirection::Windward => 0.05,
            FacadeDirection::Leeward => -0.05,
            FacadeDirection::Neither => 0.0,
            FacadeDirection::Roof => 0.,
            _ => panic!("Invalid combination of shield_class and facade_direction"),
        }
    }
}

/// this is our implementation of the numpy sign function
/// we might get precision issues with numbers really close to 0
fn sign(value: f64) -> i8 {
    if value < 0. {
        -1
    } else if value == 0. {
        0
    } else {
        1
    }
}

#[derive(Debug)]
pub(crate) struct Window {
    a_w_max: f64,
    c_d_w: f64,
    n_w: f64,
    orientation: f64,
    pitch: f64,
    on_off_ctrl_obj: Option<Arc<Control>>,
    _altitude: f64,
    p_a_alt: f64,
    z: f64,
    window_parts: Vec<WindowPart>,
}

impl Window {
    pub(crate) fn new(
        h_w_fa: f64,
        h_w_path: f64,
        a_w_max: f64,
        window_part_list: Vec<WindowPartInput>,
        orientation: f64,
        pitch: f64,
        altitude: f64,
        on_off_ctrl_obj: Option<Arc<Control>>,
        ventilation_zone_base_height: f64,
    ) -> Self {
        let n_w_div = max_of_2(window_part_list.len() - 1, 0usize) as f64;
        Self {
            a_w_max,
            c_d_w: 0.67,
            n_w: 0.5,
            orientation,
            pitch,
            on_off_ctrl_obj,
            _altitude: altitude,
            p_a_alt: adjust_air_density_for_altitude(altitude),
            z: h_w_path + ventilation_zone_base_height,
            window_parts: window_part_list
                .iter()
                .enumerate()
                .map(|(window_part_number, window_part_input)| {
                    WindowPart::new(
                        window_part_input.mid_height_air_flow_path,
                        h_w_fa,
                        n_w_div,
                        window_part_number + 1,
                        ventilation_zone_base_height,
                    )
                })
                .collect(),
        }
    }

    /// The window opening free area A_w for a window
    /// Equation 40 in BS EN 16798-7.
    /// Arguments:
    ///     R_w_arg -- ratio of window opening (0-1)
    fn calculate_window_opening_free_area(
        &self,
        r_w_arg: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        // Assume windows are shut if the control object is empty
        match &self.on_off_ctrl_obj {
            None => 0.,
            Some(ctrl) => {
                if ctrl.is_on(simulation_time_iteration) {
                    r_w_arg * self.a_w_max
                } else {
                    0.
                }
            }
        }
    }

    /// The C_w_path flow coefficient for a window
    /// Equation 54 from BS EN 16798-7
    /// Arguments:
    ///     R_w_arg -- ratio of window opening (0-1)
    fn calculate_flow_coeff_for_window(
        &self,
        r_w_arg: f64,
        simulation_time: SimulationTimeIteration,
    ) -> f64 {
        // Assume windows are shut if the control object is empty
        match &self.on_off_ctrl_obj {
            None => 0.,
            Some(ctrl) => {
                if ctrl.is_on(simulation_time) {
                    let a_w = self.calculate_window_opening_free_area(r_w_arg, simulation_time);
                    3600. * self.c_d_w * a_w * (2. / p_a_ref()).powf(self.n_w)
                } else {
                    0.
                }
            }
        }
    }

    /// Calculate the airflow through window opening based on how open the window is and internal pressure
    /// Arguments:
    /// * `wind_direction` - direction wind is blowing from, in clockwise degrees from North
    /// * `u_site` - wind velocity at zone level (m/s)
    /// * `t_e` - external air temperature (K)
    /// * `t_z` - thermal zone air temperature (K)
    /// * `p_z_ref` - internal reference pressure (Pa)
    /// * `f_cross` - boolean, dependent on if cross ventilation is possible or not
    /// * `shield_class` - indicates exposure to wind
    /// * `r_w_arg` - ratio of window opening (0-1)
    /// * `simulation_time`
    fn calculate_flow_from_internal_p(
        &self,
        wind_direction: f64,
        u_site: f64,
        t_e: f64,
        t_z: f64,
        p_z_ref: f64,
        f_cross: bool,
        shield_class: VentilationShieldClass,
        r_w_arg: Option<f64>,
        simulation_time: SimulationTimeIteration,
    ) -> (f64, f64) {
        // Assume windows are shut if the control object is empty
        let r_w_arg = match &self.on_off_ctrl_obj {
            None => 0.,
            Some(control) => {
                if !control.is_on(simulation_time) {
                    0.
                } else {
                    r_w_arg.expect("r_w_arg was None")
                }
            }
        };
        // Wind pressure coefficient for the window
        let c_p_path = get_c_p_path_from_pitch_and_orientation(
            f_cross,
            shield_class,
            self.z,
            wind_direction,
            self.orientation,
            self.pitch,
        );

        // Airflow coefficient of the window
        let c_w_path = self.calculate_flow_coeff_for_window(r_w_arg, simulation_time);

        //  Sum airflow through each window part entering and leaving - based on Equation 56 and 57
        let mut qv_in_through_window_opening = 0.;
        let mut qv_out_through_window_opening = 0.;
        for window_part in &self.window_parts {
            let air_flow = window_part.calculate_ventilation_through_windows_using_internal_p(
                u_site, t_e, t_z, c_w_path, p_z_ref, c_p_path,
            );
            if air_flow >= 0. {
                qv_in_through_window_opening += air_flow;
            } else {
                qv_out_through_window_opening += air_flow;
            }
        }

        //  Convert volume air flow rate to mass air flow rate
        convert_to_mass_air_flow_rate(
            qv_in_through_window_opening,
            qv_out_through_window_opening,
            t_e,
            t_z,
            self.p_a_alt,
        )
    }
}

#[derive(Clone, Copy, Debug)]
struct WindowPart {
    n_w_div: f64,
    h_w_div_path: f64,
    n_w: f64,
    _z: f64,
}

impl WindowPart {
    fn new(
        h_w_path: f64,
        h_w_fa: f64,
        n_w_div: f64,
        window_part_number: usize,
        ventilation_zone_base_height: f64,
    ) -> Self {
        Self {
            n_w_div,
            h_w_div_path: Self::calculate_height_for_delta_p_w_div_path(
                h_w_path,
                h_w_fa,
                n_w_div,
                window_part_number,
            ),
            n_w: 0.5,
            _z: h_w_path + ventilation_zone_base_height,
        }
    }

    /// Calculate the airflow through window parts from internal pressure
    /// Arguments:
    /// u_site -- wind velocity at zone level (m/s)
    /// T_e -- external air temperature (K)
    /// T_z -- thermal zone air temperature (K)
    /// C_w_path -- wind pressure coefficient at height of the window
    /// p_z_ref -- internal reference pressure (Pa)
    /// C_p_path -- wind pressure coefficient at the height of the window part
    fn calculate_ventilation_through_windows_using_internal_p(
        &self,
        u_site: f64,
        t_e: f64,
        t_z: f64,
        c_w_path: f64,
        p_z_ref: f64,
        c_p_path: f64,
    ) -> f64 {
        let delta_p_path = calculate_pressure_difference_at_an_airflow_path(
            self.h_w_div_path,
            c_p_path,
            u_site,
            t_e,
            t_z,
            p_z_ref,
        );

        // Based on Equation 53
        c_w_path / (self.n_w_div + 1.)
            * f64::from(sign(delta_p_path))
            * delta_p_path.abs().powf(self.n_w)
    }

    /// The height to be considered for delta_p_w_div_path
    /// Equation 55 from BS EN 16798-7
    fn calculate_height_for_delta_p_w_div_path(
        h_w_path: f64,
        h_w_fa: f64,
        n_w_div: f64,
        window_part_number: usize,
    ) -> f64 {
        h_w_path - h_w_fa / 2.
            + h_w_fa / (2. * (n_w_div + 1.))
            + (h_w_fa / (n_w_div + 1.)) * (window_part_number - 1) as f64
    }
}

#[derive(Debug)]
pub(crate) struct Vent {
    h_path: f64,
    a_vent: f64,
    delta_p_vent_ref: f64,
    orientation: f64,
    pitch: f64,
    _altitude: f64,
    n_vent: f64,
    c_d_vent: f64,
    p_a_alt: f64,
    z: f64,
    // NOTE - in Python we have C_vent_path as an instance variable but here we calculate it when needed instead
}

impl Vent {
    /// Construct a Vent object
    ///
    /// Arguments:
    ///    h_path -- mid-height of air flow path relative to ventilation zone (m)
    ///    A_vent - Equivalent area of a vent (m2)
    ///    delta_p_vent_ref -- reference pressure difference for vent (Pa)
    ///    orientation -- The orientation of the vent (degrees)
    ///    pitch -- The pitch of the vent (degrees)
    ///    altitude -- altitude of dwelling above sea level (m)
    ///
    /// Method:
    ///    - Based on Section 6.4.3.6 Airflow through vents from BS EN 16798-7
    pub(crate) fn new(
        h_path: f64,
        a_vent: f64,
        delta_p_vent_ref: f64,
        orientation: f64,
        pitch: f64,
        altitude: f64,
        ventilation_zone_base_height: f64, // TODO: added as part of the 0.32 migration, still WIP
    ) -> Self {
        Self {
            h_path,
            a_vent,
            delta_p_vent_ref,
            orientation,
            pitch,
            _altitude: altitude,
            n_vent: 0.5, // Flow exponent for vents based on Section B.3.2.2 from BS EN 16798-7
            c_d_vent: 0.6, // Discharge coefficient of vents based on B.3.2.1 from BS EN 16798-7
            p_a_alt: adjust_air_density_for_altitude(altitude),
            z: h_path + ventilation_zone_base_height,
        }
    }

    /// The vent opening free area A_vent for a vent
    /// Arguments:
    /// * `r_v_arg` - ratio of vent opening (0-1)
    fn calculate_vent_opening_free_area(&self, r_v_arg: f64) -> f64 {
        r_v_arg * self.a_vent
    }

    /// The airflow coefficient of the vent calculated from equivalent area A_vent_i
    /// according to EN 13141-1 and EN 13141-2.
    /// Based on Equation 59 from BS EN 16798-7.
    fn calculate_flow_coeff_for_vent(&self, r_v_arg: f64) -> f64 {
        // NOTE: The standard does not define what the below 3600 and 10000 are.

        let a_vent = self.calculate_vent_opening_free_area(r_v_arg);
        (3600. / 10000.)
            * self.c_d_vent
            * a_vent
            * (2. / p_a_ref()).powf(0.5)
            * (1. / self.delta_p_vent_ref).powf(self.n_vent - 0.5)
    }

    /// Calculate the airflow through vents from internal pressure
    /// Arguments:
    /// * `u_site` - wind velocity at zone level (m/s)
    /// * `t_e` - external air temperature (K)
    /// * `t_z` -- thermal zone air temperature (K)
    /// * `c_vent_path` - wind pressure coefficient at height of the vent
    /// * `c_p_path` - wind pressure coefficient at the height of the window part
    /// * `p_z_ref` - internal reference pressure (Pa)
    fn calculate_ventilation_through_vents_using_internal_p(
        &self,
        u_site: f64,
        t_e: f64,
        t_z: f64,
        c_vent_path: f64,
        c_p_path: f64,
        p_z_ref: f64,
    ) -> f64 {
        // Pressure_difference at the vent level
        let delta_p_path = calculate_pressure_difference_at_an_airflow_path(
            self.h_path,
            c_p_path,
            u_site,
            t_e,
            t_z,
            p_z_ref,
        );

        // Air flow rate for each couple of height and wind pressure coeficient associated with vents.
        // Based on Equation 58
        c_vent_path * sign(delta_p_path) as f64 * delta_p_path.abs().powf(self.n_vent)
    }

    /// Calculate the airflow through vents from internal pressure
    ///
    /// Arguments:
    /// * `wind_direction` - direction wind is blowing from, in clockwise degrees from North
    /// * `u_site` - wind velocity at zone level (m/s)
    /// * `t_e` - external air temperature (K)
    /// * `t_z` - thermal zone air temperature (K)
    /// * `p_z_ref` - internal reference pressure (Pa)
    /// * `f_cross` - boolean, dependent on if cross ventilation is possible or not
    /// * `shield_class` - indicates exposure to wind
    /// * `r_v_arg`
    fn calculate_flow_from_internal_p(
        &self,
        wind_direction: f64,
        u_site: f64,
        t_e: f64,
        t_z: f64,
        p_z_ref: f64,
        f_cross: bool,
        shield_class: VentilationShieldClass,
        r_v_arg: f64,
    ) -> (f64, f64) {
        // Wind pressure coefficient for the air flow path
        let c_p_path = get_c_p_path_from_pitch_and_orientation(
            f_cross,
            shield_class,
            self.z,
            wind_direction,
            self.orientation,
            self.pitch,
        );

        let c_vent_path = self.calculate_flow_coeff_for_vent(r_v_arg);
        // Calculate airflow through each vent
        let air_flow = self.calculate_ventilation_through_vents_using_internal_p(
            u_site,
            t_e,
            t_z,
            c_vent_path,
            c_p_path,
            p_z_ref,
        );

        // Sum airflows entering and leaving - based on Equation 60 and 61
        let mut qv_in_through_vent = 0.;
        let mut qv_out_through_vent = 0.;

        if air_flow >= 0. {
            qv_in_through_vent += air_flow
        } else {
            qv_out_through_vent += air_flow
        }

        // Convert volume air flow rate to mass air flow rate
        convert_to_mass_air_flow_rate(
            qv_in_through_vent,
            qv_out_through_vent,
            t_e,
            t_z,
            self.p_a_alt,
        )
    }
}

// NOTE - In the python implementation this is a property of Leaks
// low exponent through leaks based on value in B.3.3.14
const N_LEAK: f64 = 0.667;

/// An object to represent Leaks
#[derive(Debug)]
struct Leaks {
    h_path: f64,
    delta_p_leak_ref: f64,
    a_roof: f64,
    a_facades: f64,
    a_leak: f64,
    qv_delta_p_leak_ref: f64,
    facade_direction: FacadeDirection,
    _altitude: f64,
    p_a_alt: f64,
    z: f64,
    // In Python there are extra properties:
    // n_leak - this is now N_LEAK as it is constant
    // c_leak_path - this is now calculated when needed with calculate_flow_coeff_for_leak
}

impl Leaks {
    /// Arguments:
    /// * `h_path` - mid-height of the air flow path relative to ventilation zone floor level
    /// * `delta_p_leak_ref` - Reference pressure difference (From pressure test e.g. blower door = 50Pa)
    /// * `qv_delta_p_leak_ref` - flow rate through
    /// * `facade_direction` - The direction of the facade the leak is on.
    /// * `a_roof` - Surface area of the roof of the ventilation zone (m2)
    /// * `a_facades` - Surface area of facades (m2)
    /// * `a_leak` - Reference area of the envelope airtightness index qv_delta_p_leak_ref (depends on national context)
    /// * `altitude` - altitude of dwelling above sea level (m)
    /// * `ventilation_zone_base_height` - Base height of the ventilation zone relative to ground (m)
    ///
    /// Based on Section 6.4.3.6 Airflow through leaks from BS EN 16798-7.
    //
    fn new(
        h_path: f64,
        delta_p_leak_ref: f64,
        qv_delta_p_leak_ref: f64,
        facade_direction: FacadeDirection,
        a_roof: f64,
        a_facades: f64,
        a_leak: f64,
        altitude: f64,
        ventilation_zone_base_height: f64,
    ) -> Self {
        Self {
            h_path,
            delta_p_leak_ref,
            a_roof,
            a_facades,
            a_leak,
            qv_delta_p_leak_ref,
            facade_direction,
            _altitude: altitude,
            p_a_alt: adjust_air_density_for_altitude(altitude),
            z: h_path + ventilation_zone_base_height,
        }
    }

    fn calculate_flow_coeff_for_leak(&self) -> f64 {
        //  c_leak - Leakage coefficient of ventilation zone

        let c_leak = self.qv_delta_p_leak_ref * self.a_leak / (self.delta_p_leak_ref).powf(N_LEAK);

        // Leakage coefficient of roof, estimated to be proportional to ratio
        // of surface area of the facades to that of the facades plus the roof.
        if self.facade_direction != FacadeDirection::Windward
            && self.facade_direction != FacadeDirection::Leeward
        {
            // leak in roof

            let c_leak_roof = c_leak * self.a_roof / (self.a_facades + self.a_roof);
            return c_leak_roof;
        }

        // Leakage coefficient of facades, estimated to be proportional to ratio
        // of surface area of the roof to that of the facades plus the roof.
        let c_leak_facades = c_leak * self.a_facades / (self.a_facades + self.a_roof);
        0.25 * c_leak_facades // Table B.12
    }

    /// Calculate the airflow through leaks from internal pressure
    /// Arguments:
    ///      u_site -- wind velocity at zone level (m/s)
    ///      t_e -- external air temperature (K)
    ///      t_z -- thermal zone air temperature (K)
    ///      c_p_path -- wind pressure coefficient at the height of the window part
    ///      p_z_ref -- internal reference pressure (Pa)
    fn calculate_ventilation_through_leaks_using_internal_p(
        &self,
        u_site: f64,
        t_e: f64,
        t_z: f64,
        c_p_path: f64,
        p_z_ref: f64,
    ) -> f64 {
        // For each couple of height and wind pressure coefficient associated with vents,
        // the air flow rate.
        let delta_p_path = calculate_pressure_difference_at_an_airflow_path(
            self.h_path,
            c_p_path,
            u_site,
            t_e,
            t_z,
            p_z_ref,
        );

        let c_leak_path = Self::calculate_flow_coeff_for_leak(self);

        // Airflow through leaks based on Equation 62
        c_leak_path * f64::from(sign(delta_p_path)) * delta_p_path.abs().powf(N_LEAK)
    }

    fn calculate_flow_from_internal_p(
        &self,
        u_site: f64,
        t_e: f64,
        t_z: f64,
        p_z_ref: f64,
        f_cross: bool,
        shield_class: VentilationShieldClass,
    ) -> (f64, f64) {
        // Wind pressure coefficient for the air flow path
        let c_p_path = get_c_p_path(f_cross, shield_class, self.z, self.facade_direction); // #TABLE from annex B

        // Calculate airflow through each leak
        let mut qv_in_through_leak = 0.;
        let mut qv_out_through_leak = 0.;
        let air_flow = self.calculate_ventilation_through_leaks_using_internal_p(
            u_site, t_e, t_z, c_p_path, p_z_ref,
        );

        // Add airflow entering and leaving through leak
        if air_flow >= 0. {
            qv_in_through_leak += air_flow;
        } else {
            qv_out_through_leak += air_flow;
        }

        // Convert volume air flow rate to mass air flow rate
        let (qm_in_through_leak, qm_out_through_leak) = convert_to_mass_air_flow_rate(
            qv_in_through_leak,
            qv_out_through_leak,
            t_e,
            t_z,
            self.p_a_alt,
        );

        (qm_in_through_leak, qm_out_through_leak)
    }
}

/// An object to represent AirTerminalDevices
#[derive(Clone, Copy, Debug)]
pub(crate) struct AirTerminalDevices {
    c_d_atd: f64,
    n_atd: f64,
    a_atd: f64,
    delta_p_atd_ref: f64,
    // NOTE - in Python we have c_atd_path as an instance variable but here we calculate it when needed instead
}

// remove following directive when things start referencing AirTerminalDevices (upstream is not fully implemented as of 0.30)
#[allow(dead_code)]
impl AirTerminalDevices {
    /// Construct a AirTerminalDevices object
    /// Arguments:
    /// a_atd -- equivalent area of the air terminal device (m2)
    /// delta_p_atd_ref -- Reference pressure difference for an air terminal device (Pa)
    /// Method: Based on Section 6.4.3.2.2 from BS EN 16798-7
    fn new(a_atd: f64, delta_p_atd_ref: f64) -> Self {
        Self {
            c_d_atd: 0.6, // Discharge coefficient for air terminal devices based on B.3.2.1
            n_atd: 0.5,   // Flow exponent of air terminal devices based on B.3.2.2
            a_atd,
            delta_p_atd_ref,
        }
    }

    /// The airflow coefficient of the ATD is calculated
    /// from the equivalent area A_vent value, according to
    /// EN 13141-1 and EN 13141-2.
    /// Equation 26 from BS EN 16798-7.
    fn calculate_flow_coeff_for_atd(&self) -> f64 {
        // NOTE: The standard does not define what the below 3600 and 10000 are.
        (3600. / 10000.)
            * self.c_d_atd
            * self.a_atd
            * (2. / p_a_ref()).powf(0.5)
            * (1. / self.delta_p_atd_ref).powf(self.n_atd - 0.5)
    }

    /// The pressure loss at internal air terminal devices is calculated from
    /// the total air flow rate passing through the device.
    /// Equation 25 from BS EN 16798-7.
    /// Solving for qv_pdu.
    /// Arguments:
    /// qv_pdu - volume flow rate through passive and hybrid ducts.
    fn calculate_pressure_difference_atd(&self, qv_pdu: f64) -> f64 {
        let c_atd_path = self.calculate_flow_coeff_for_atd();

        -(f64::from(sign(qv_pdu))) * (qv_pdu.abs() / c_atd_path).powf(1. / self.n_atd)
    }
}

/// An object to represent Cowls
// remove following directive when something references Cowls type
#[allow(dead_code)]
struct Cowls {
    c_p_cowl_roof: f64,
    height: f64,
    // NOTE - in Python we have delta_cowl_height as an instance variable but here we calculate it when needed from the height
}

// remove following directive when something references Cowls type
#[allow(dead_code)]
impl Cowls {
    /// Construct a Cowls object
    /// Arguments:
    /// height - height Between the top of the roof and the roof outlet in m (m)
    fn new(height: f64) -> Self {
        Self {
            c_p_cowl_roof: 0., // Default B.3.3.5
            height,
        }
    }

    /// Interpreted Table B.9 from BS EN 16798-7
    /// Get values for delta_C_cowl_height.
    /// Arguments:
    /// height - height Between the top of the roof and the roof outlet in m (m)
    fn get_delta_cowl_height(&self) -> f64 {
        if self.height < 0.5 {
            -0.0
        } else if (0.5..=1.).contains(&self.height) {
            -0.1
        } else {
            0.2
        }
    }
}

/// An object to represent CombustionAppliances
#[derive(Clone, Copy, Debug)]
pub(crate) struct CombustionAppliances {
    f_as: f64,
    f_ff: f64,
}

impl CombustionAppliances {
    /// Construct a CombustionAppliances object
    /// Arguments:
    /// f_op_comp -- Operation requirement signal (combustion appliance) (0 =  OFF ; 1 = ON)
    /// supply_situation - Combustion air supply situation: 'room_air' or 'outside'
    /// exhaust_situation - flue gas exhaust situation: 'into_room', 'into_separate_duct' or 'into_mech_vent'
    /// f_ff -- combustion air flow factor
    /// p_h_fi - Combustion appliance heating fuel input power (kW)
    pub(crate) fn new(
        supply_situation: CombustionAirSupplySituation,
        exhaust_situation: FlueGasExhaustSituation,
        fuel_type: CombustionFuelType,
        appliance_type: CombustionApplianceType,
    ) -> Self {
        Self {
            f_as: get_appliance_system_factor(supply_situation, exhaust_situation), // Combustion appliance system factor (0 or 1)
            f_ff: get_fuel_flow_factor(fuel_type, appliance_type), // Fuel flow factor (0-5)
        }
    }

    /// Calculate additional air flow rate required for the operation of
    /// combustion appliance q_v_comb.
    /// Arguments:
    /// f_op_comp --  Operation requirement signal (combustion appliance) (0 =  OFF ; 1 = ON)
    /// p_h_fi -- Combustion appliance heating fuel input power (kW)
    fn calculate_air_flow_req_for_comb_appliance(&self, f_op_comp: f64, p_h_fi: f64) -> (f64, f64) {
        // TODO (from the Python) flue is considered as vertical passive duct, standard formulas in CIBSE guide.
        let q_v_in_through_comb = 0.; // (37)
        let mut q_v_out_through_comb = 0.; // (38)

        if f_op_comp == 1. {
            let q_v_comb = 3.6 * f_op_comp * self.f_as * self.f_ff * p_h_fi; // (35)
            q_v_out_through_comb = -q_v_comb; // temp associated with q_v_out_through_comb is ventilation zone temperature tz
        }

        (q_v_in_through_comb, q_v_out_through_comb)
    }
}

/// An object to represent Mechanical Ventilation
#[derive(Debug)]
pub(crate) struct MechanicalVentilation {
    // theta_z_t will be referenced once upstream reference is implementable
    _theta_z_t: f64,
    sup_air_flw_ctrl: SupplyAirFlowRateControlType,
    _sup_air_temp_ctrl: SupplyAirTemperatureControlType,
    _q_h_des: f64,
    _q_c_des: f64,
    _theta_ctrl_sys: Option<f64>,
    vent_type: VentType,
    _total_volume: f64,
    ctrl_intermittent_mev: Option<Arc<Control>>,
    sfp: f64,
    energy_supply_conn: EnergySupplyConnection,
    _altitude: f64,
    pub(crate) design_outdoor_air_flow_rate_m3_h: f64,
    mvhr_eff: f64,
    qv_oda_req_design: f64,
    p_a_alt: f64,
}

impl MechanicalVentilation {
    /// Construct a Mechanical Ventilation object
    /// Arguments:
    /// sup_air_flw_ctrl -- supply air flow rate control
    /// sup_air_temp_ctrl --supply air temperature control
    /// q_h_des -- design zone heating need to be covered by the mechanical ventilation system
    /// q_c_des -- design zone cooling need to be covered by the mechanical ventilation system
    /// vent_type -- ventilation system type
    /// specific_fan_power -- in W / (litre / second), inclusive of any in use factors
    /// design_outdoor_air_flow_rate -- design outdoor air flow rate in m3/h
    /// simulation_time -- reference to Simulation time object
    /// energy_supply_conn -- Energy supply connection
    /// total_volume  -- Total zone volume (m3)
    /// altitude -- altitude of dwelling above sea level (m)
    /// ctrl_intermittent_MEV -- reference to Control object with boolean schedule
    /// defining when the MechVent should be on.
    /// mvhr_eff -- MVHR efficiency
    /// theta_ctrl_sys -- Temperature variation based on control system (K)
    pub(crate) fn new(
        _sup_air_flw_ctrl: SupplyAirFlowRateControlType,
        _sup_air_temp_ctrl: SupplyAirTemperatureControlType,
        q_h_des: f64,
        q_c_des: f64,
        vent_type: VentType,
        specific_fan_power: f64,
        design_outdoor_air_flow_rate: f64,
        energy_supply_conn: EnergySupplyConnection,
        total_volume: f64,
        altitude: f64,
        ctrl_intermittent_mev: Option<Arc<Control>>,
        mvhr_eff: Option<f64>,
        theta_ctrl_sys: Option<f64>, // Only required if sup_air_temp_ctrl = LOAD_COM
    ) -> Self {
        let f_ctrl = 1.; // From table B.4, for residential buildings, default f_ctrl = 1
        let f_sys = 1.1; // From table B.5, f_sys = 1.1
        let e_v = 1.; // Section B.3.3.7 defaults E_v = 1 (this is the assumption for perfect mixing)

        Self {
            _theta_z_t: 0., // TODO (from Python) get Thermal zone temperature - used for LOAD
            sup_air_flw_ctrl: SupplyAirFlowRateControlType::ODA, // TODO (from Python) currently hard coded until load comp implemented
            _sup_air_temp_ctrl: SupplyAirTemperatureControlType::NoControl, // TODO (from Python) currently hard coded until load comp implemented
            _q_h_des: q_h_des,
            _q_c_des: q_c_des,
            _theta_ctrl_sys: theta_ctrl_sys,
            vent_type,
            _total_volume: total_volume,
            ctrl_intermittent_mev,
            sfp: specific_fan_power,
            energy_supply_conn,
            _altitude: altitude,
            design_outdoor_air_flow_rate_m3_h: design_outdoor_air_flow_rate, // in m3/h
            mvhr_eff: mvhr_eff.unwrap_or(0.0),
            // Calculated variables
            qv_oda_req_design: Self::calculate_required_outdoor_air_flow_rate(
                f_ctrl,
                f_sys,
                e_v,
                design_outdoor_air_flow_rate,
            ),
            p_a_alt: adjust_air_density_for_altitude(altitude),
        }
    }

    /// Calculate required outdoor ventilation air flow rates.
    /// Equation 9 from BS EN 16798-7.
    fn calculate_required_outdoor_air_flow_rate(
        f_ctrl: f64,
        f_sys: f64,
        e_v: f64,
        design_outdoor_air_flow_rate_m3_h: f64,
    ) -> f64 {
        // Required outdoor air flow rate in m3/h
        ((f_ctrl * f_sys) / e_v) * design_outdoor_air_flow_rate_m3_h
    }

    /// Calculate required outdoor air flow rates at the air terminal devices
    /// Equations 10-17 from BS EN 16798-7
    /// Adjusted to be based on ventilation type instead of vent_sys_op.
    fn calc_req_oda_flow_rates_at_atds(&self) -> (f64, f64) {
        match &self.vent_type {
            VentType::Mvhr => (self.qv_oda_req_design, -self.qv_oda_req_design),
            VentType::IntermittentMev
            | VentType::CentralisedContinuousMev
            | VentType::DecentralisedContinuousMev => {
                // NOTE: Calculation of effective flow rate of external air (in func
                // calc_mech_vent_air_flw_rates_req_to_supply_vent_zone) assumes that
                // supply and extract are perfectly balanced (as defined above), so
                // any future change to this assumption will need to be considered
                // with that in mind
                (0., -self.qv_oda_req_design)
            }
            VentType::Piv => (self.qv_oda_req_design, 0.),
        }
    }

    /// Returns the fraction of the timestep for which the ventilation is running
    fn f_op_v(&self, simulation_time: &SimulationTimeIteration) -> f64 {
        match &self.vent_type {
            VentType::IntermittentMev => {
                let f_op_v = self
                    .ctrl_intermittent_mev
                    .as_ref()
                    .expect("ctrl_intermittent_mev was expected to be set")
                    .setpnt(simulation_time)
                    .expect("A setpoint was expected to be derivable for a control.");

                if !(0. ..=1.).contains(&f_op_v) {
                    panic!("Error f_op_v is not between 0 and 1")
                }

                f_op_v
            }
            VentType::DecentralisedContinuousMev
            | VentType::CentralisedContinuousMev
            | VentType::Mvhr => {
                // Assumed to operate continuously
                1.
            }
            // NOTE - this will happen for VentType::Piv
            // same behaviour as Python
            _ => panic!("Unknown mechanical ventilation system type"),
        }
    }

    /// Calculate the air flow rates to and from the ventilation zone required from mechanical ventilation.
    /// T_z -- thermal zone temperature (K)
    fn calc_mech_vent_air_flw_rates_req_to_supply_vent_zone(
        &self,
        t_z: f64,
        t_e: f64,
        simulation_time: &SimulationTimeIteration,
    ) -> Result<(f64, f64, f64), NotImplementedError> {
        // Required air flow at air terminal devices
        let (qv_sup_req, qv_eta_req) = self.calc_req_oda_flow_rates_at_atds();

        // Amount of air flow depends on controls
        let (qv_sup_dis_req, qv_eta_dis_req) = match self.sup_air_flw_ctrl {
            SupplyAirFlowRateControlType::ODA => {
                let f_op_v = self.f_op_v(simulation_time);
                let qv_sup_dis_req = f_op_v * qv_sup_req;
                let qv_eta_dis_req = f_op_v * qv_eta_req;

                (qv_sup_dis_req, qv_eta_dis_req)
            }
            SupplyAirFlowRateControlType::Load => {
                // NOTE - this is not currently implemented in the Python code
                // reported up to BRE https://dev.azure.com/BreGroup/SAP%2011/_workitems/edit/45523
                // and not yet fixed as of 0.33
                return Err(NotImplementedError::new("calc_mech_vent_air_flw_rates_req_to_supply_vent_zone is not implemented for SupplyAirFlowRateControlType::Load as outstanding bug in upstream Python"));
            }
        };

        // Calculate effective flow rate of external air
        // NOTE: Technically, the MVHR system supplies air at a higher
        // temperature than the outside air. However, it is simpler to
        // account for the heat recovery effect using an "equivalent" or
        // "effective" flow rate of external air
        let qv_effective_heat_recovery_saving = qv_sup_dis_req * self.mvhr_eff;

        // Convert volume air flow rate to mass air flow rate
        let (qm_sup_dis_req, qm_eta_dis_req) =
            convert_to_mass_air_flow_rate(qv_sup_dis_req, qv_eta_dis_req, t_e, t_z, self.p_a_alt);

        let qm_in_effective_heat_recovery_saving = convert_volume_flow_rate_to_mass_flow_rate(
            qv_effective_heat_recovery_saving,
            t_e,
            self.p_a_alt,
        );

        Ok((
            qm_sup_dis_req,
            qm_eta_dis_req,
            qm_in_effective_heat_recovery_saving,
        ))
    }

    /// Calculate gains and energy use due to fans
    /// zone_volume -- volume of the zone (m3)
    /// total_volume -- volume of the dwelling (m3)
    /// vent_type -- one of "Intermittent MEV", "Centralised continuous MEV",
    /// "Decentralised continuous MEV", "MVHR" or "PIV".
    pub(crate) fn fans(
        &self,
        zone_volume: f64,
        total_volume: f64,
        throughput_factor: Option<f64>,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        let _throughput_factor = throughput_factor.unwrap_or(1.0);
        // Calculate energy use by fans
        let fan_power_w = (self.sfp
            * (self.qv_oda_req_design / f64::from(SECONDS_PER_HOUR))
            * f64::from(LITRES_PER_CUBIC_METRE))
            * (zone_volume / total_volume);
        let fan_energy_use_kwh = (fan_power_w / f64::from(WATTS_PER_KILOWATT))
            * simulation_time_iteration.timestep
            * self.f_op_v(simulation_time_iteration);

        let (supply_fan_energy_use_kwh, extract_fan_energy_use_in_kwh) = match self.vent_type {
            VentType::IntermittentMev
            | VentType::CentralisedContinuousMev
            | VentType::DecentralisedContinuousMev => {
                // Fan energy use = 0
                (0.0, fan_energy_use_kwh)
            }
            VentType::Mvhr => {
                // Balanced, therefore split power between extract and supply fans
                (fan_energy_use_kwh / 2., fan_energy_use_kwh / 2.)
            }
            VentType::Piv => {
                // Positive input, supply fans only
                (fan_energy_use_kwh, 0.)
            }
        };
        self.energy_supply_conn
            .demand_energy(supply_fan_energy_use_kwh, simulation_time_iteration.index)
            .unwrap();
        self.energy_supply_conn
            .demand_energy(
                extract_fan_energy_use_in_kwh,
                simulation_time_iteration.index,
            )
            .unwrap();

        supply_fan_energy_use_kwh
            / (f64::from(WATTS_PER_KILOWATT) * simulation_time_iteration.timestep)
    }

    pub fn vent_type(&self) -> VentType {
        self.vent_type
    }
}

/// A class to represent Infiltration and Ventilation object
#[derive(Debug)]
pub(crate) struct InfiltrationVentilation {
    f_cross: bool,
    shield_class: VentilationShieldClass,
    c_rgh_site: f64,
    ventilation_zone_height: f64,
    windows: Vec<Window>,
    vents: Vec<Vent>,
    leaks: Vec<Leaks>,
    combustion_appliances: Vec<CombustionAppliances>,
    air_terminal_devices: Vec<AirTerminalDevices>,
    mech_vents: Vec<Arc<MechanicalVentilation>>,
    space_heating_ductworks: Arc<IndexMap<String, Vec<Ductwork>>>,
    detailed_output_heating_cooling: bool,
    p_a_alt: f64,
    total_volume: f64,
    detailed_results: Arc<RwLock<Vec<VentilationDetailedResult>>>,
    #[cfg(test)] // optional behaviour override for tests, akin to mocking
    calc_air_changes_fn: Option<CalcAirChangesFn>,
}

type CalcAirChangesFn = fn(
    &InfiltrationVentilation,
    f64,
    f64,
    f64,
    f64,
    f64,
    Option<f64>,
    f64,
    Option<ReportingFlag>,
    SimulationTimeIteration,
) -> anyhow::Result<f64>;

/// Arguments:
/// * `f_cross` - cross-ventilation factor
/// * `shield_class` - indicates the exposure to wind of an air flow path on a facade
///                    (can can be open, normal or shielded)
/// * `terrain class`
/// * `average_roof_pitch`
/// * `windows` - list of windows
/// * `vents` - list of vents
/// * `leaks` - required inputs for leaks
/// * `combustion_appliances`
/// * `air_terminal_devices` - list of air terminal devices
/// * `mech_vents` - list of mech vents
/// * `space_heating_ductworks`
/// * `detailed_output_heating_cooling` - whether to output detailed heating/cooling data
/// * `altitude` - altitude of dwelling above sea level (m)
/// * `total_volume` - total zone volume
/// * `ventilation_zone_base_height` -- base height of the ventilation zone (m)
impl InfiltrationVentilation {
    pub(crate) fn new(
        f_cross: bool,
        shield_class: VentilationShieldClass,
        terrain_class: &TerrainClass,
        average_roof_pitch: f64,
        windows: Vec<Window>,
        vents: Vec<Vent>,
        leaks: CompletedVentilationLeaks,
        combustion_appliances: Vec<CombustionAppliances>,
        air_terminal_devices: Vec<AirTerminalDevices>,
        mech_vents: Vec<Arc<MechanicalVentilation>>,
        space_heating_ductworks: Arc<IndexMap<String, Vec<Ductwork>>>,
        detailed_output_heating_cooling: bool,
        altitude: f64,
        total_volume: f64,
        ventilation_zone_base_height: f64,
    ) -> Self {
        let ventilation_zone_height = leaks.ventilation_zone_height;
        Self {
            f_cross,
            shield_class,
            c_rgh_site: ter_class_to_roughness_coeff(
                terrain_class,
                ventilation_zone_base_height + ventilation_zone_height / 2.,
            ),
            ventilation_zone_height,
            windows,
            vents,
            leaks: Self::make_leak_objects(
                leaks,
                average_roof_pitch,
                ventilation_zone_base_height,
                f_cross,
            ),
            combustion_appliances,
            air_terminal_devices,
            mech_vents,
            space_heating_ductworks,
            detailed_output_heating_cooling,
            p_a_alt: adjust_air_density_for_altitude(altitude),
            total_volume,
            detailed_results: Default::default(),
            #[cfg(test)]
            calc_air_changes_fn: None,
        }
    }

    #[cfg(test)] // set override of ach function under test scenario
    fn set_calc_air_changes_fn(&mut self, calc_air_changes_fn: CalcAirChangesFn) {
        self.calc_air_changes_fn.replace(calc_air_changes_fn);
    }

    pub(crate) fn mech_vents(&self) -> &Vec<Arc<MechanicalVentilation>> {
        &self.mech_vents
    }

    pub(crate) fn space_heating_ductworks(&self) -> &Arc<IndexMap<String, Vec<Ductwork>>> {
        &self.space_heating_ductworks
    }

    /// Calculate total volume air flow rate entering ventilation zone
    /// Equation 68 from BS EN 16798-7
    #[cfg(test)]
    fn calculate_total_volume_air_flow_rate_in(qm_in: f64, external_air_density: f64) -> f64 {
        qm_in / external_air_density // from weather file?
    }

    /// Calculate total volume air flow rate leaving ventilation zone
    /// Equation 69 from BS EN 16798-7
    #[cfg(test)]
    fn calculate_total_volume_air_flow_rate_out(qm_out: f64, zone_air_density: f64) -> f64 {
        qm_out / zone_air_density
    }

    /// Distribute leaks around the dwelling according to Table B.12 from BS EN 16798-7.
    /// Create 5 leak objects:
    ///     At 0.25*Height of the Ventilation Zone in the Windward facade
    ///     At 0.25*Height of the Ventilation Zone in the Leeward facade
    ///     At 0.75*Height of the Ventilation Zone in the Windward facade
    ///     At 0.75*Height of the Ventilation Zone in the Leeward facade
    ///     At the Height of the Ventilation Zone in the roof
    /// Arguments:
    /// leak - dict of leaks input data from JSON file
    /// average_roof_pitch - calculated in project.py, average pitch of all roof elements weighted by area (degrees)
    fn make_leak_objects(
        leaks: CompletedVentilationLeaks,
        average_roof_pitch: f64,
        ventilation_zone_base_height: f64,
        f_cross: bool,
    ) -> Vec<Leaks> {
        let h_path1_2 = 0.25 * leaks.ventilation_zone_height;
        let h_path3_4 = 0.75 * leaks.ventilation_zone_height;
        let h_path5 = leaks.ventilation_zone_height;
        let h_path_list = [h_path1_2, h_path1_2, h_path3_4, h_path3_4, h_path5];

        let roof_pitch = if f_cross {
            match average_roof_pitch {
                ..10.0 => FacadeDirection::Roof10,
                10.0..=30.0 => FacadeDirection::Roof10_30,
                30.0..60.0 => FacadeDirection::Roof30,
                _ => panic!("Average roof pitch was not expected to be greater than 60 degrees."),
            }
        } else {
            FacadeDirection::Roof
        };

        let facade_direction = [
            FacadeDirection::Windward,
            FacadeDirection::Leeward,
            FacadeDirection::Windward,
            FacadeDirection::Leeward,
            roof_pitch,
        ];

        (0..5)
            .map(|i| {
                Leaks::new(
                    h_path_list[i],
                    leaks.test_pressure,
                    leaks.test_result,
                    facade_direction[i],
                    leaks.area_roof,
                    leaks.area_facades,
                    leaks.env_area,
                    leaks.altitude,
                    ventilation_zone_base_height,
                )
            })
            .collect()
    }

    /// Implicit solver for qv_pdu
    fn calculate_qv_pdu(
        &self,
        qv_pdu: f64,
        p_z_ref: f64,
        t_z: f64,
        t_e: f64,
        h_z: f64,
    ) -> anyhow::Result<f64> {
        let func = |qv_pdu, [p_z_ref, t_z, h_z]: [f64; 3]| {
            Ok(self.implicit_formula_for_qv_pdu(qv_pdu, p_z_ref, t_z, t_e, h_z))
        };

        fsolve(func, qv_pdu, [p_z_ref, t_z, h_z]) // returns qv_pdu
    }

    /// Implicit formula solving for qv_pdu as unknown.
    /// Equation 30 from BS EN 16798-7
    /// Arguments:
    /// qv_pdu -- volume flow rate from passive and hybrid ducts (m3/h)
    /// p_z_ref -- internal reference pressure (Pa)
    /// T_z -- thermal zone temperature (K)
    /// h_z -- height of ventilation zone (m)
    fn implicit_formula_for_qv_pdu(
        &self,
        qv_pdu: f64,
        p_z_ref: f64,
        t_z: f64,
        t_e: f64,
        h_z: f64,
    ) -> f64 {
        let external_air_density = air_density_at_temp(t_e, self.p_a_alt);
        let zone_air_density = air_density_at_temp(t_z, self.p_a_alt);

        // TODO (from Python) Standard isn't clear if delta_p_ATD can be totalled or not.
        let delta_p_atd_list: Vec<f64> = self
            .air_terminal_devices
            .iter()
            .map(|atd| atd.calculate_pressure_difference_atd(qv_pdu))
            .collect();

        let delta_p_atd: f64 = delta_p_atd_list.iter().sum();

        // Stack effect in passive and hybrid duct. As there is no air transfer
        // between levels of the ventilation zone Equation B.1 is used.
        let h_pdu_stack = h_z + 2.;

        // TODO (from Python) include delta_p_dpu and delta_p_cowl in the return.
        delta_p_atd - p_z_ref - h_pdu_stack * G * (external_air_density - zone_air_density)
    }

    /// The root scalar function will iterate until it finds a value of p_z_ref
    /// that satisfies the mass balance equation.
    /// The root scalar solver allows a range of intervals to be entered.
    /// The loop begins with a small interval to start with and if no solution is
    /// found or the boundary is too small for to cause a sign change then a wider
    /// interval is used until a solution is found.
    pub(crate) fn calculate_internal_reference_pressure(
        &self,
        initial_p_z_ref_guess: f64,
        wind_speed: f64,
        wind_direction: f64,
        temp_int_air: f64,
        temp_ext_air: f64,
        r_v_arg: f64,
        r_w_arg: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> Result<f64, InternalReferencePressureCalculationError> {
        for interval_expansion in INTERVAL_EXPANSION_LIST {
            let bracket = (
                initial_p_z_ref_guess - interval_expansion,
                initial_p_z_ref_guess + interval_expansion,
            );

            let result = root_scalar_for_implicit_mass_balance(
                self,
                wind_speed,
                wind_direction,
                temp_int_air,
                temp_ext_air,
                r_v_arg,
                r_w_arg,
                simtime,
                bracket,
            );

            if let Ok(root) = result {
                let p_z_ref = root;
                return Ok(p_z_ref);
            }
        }

        Err(InternalReferencePressureCalculationError {
            initial_p_z_ref_guess,
            temp_int_air,
            r_w_arg,
        })
    }

    /// Used in calculate_internal_reference_pressure function for p_z_ref solve
    pub(crate) fn implicit_mass_balance_for_internal_reference_pressure(
        &self,
        p_z_ref: f64,
        wind_speed: f64,
        wind_direction: f64,
        temp_int_air: f64,
        temp_ext_air: f64,
        r_v_arg: f64,
        r_w_arg_min_max: Option<f64>,
        flag: Option<ReportingFlag>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let (qm_in, qm_out, _) = self
            .implicit_mass_balance_for_internal_reference_pressure_components(
                p_z_ref,
                wind_speed,
                wind_direction,
                temp_int_air,
                temp_ext_air,
                r_v_arg,
                r_w_arg_min_max,
                flag,
                simtime,
            )?;
        Ok(qm_in + qm_out)
    }

    /// Calculate incoming air flow, in m3/hr, at specified conditions
    pub(crate) fn incoming_air_flow(
        &self,
        p_z_ref: f64,
        wind_speed: f64,
        wind_direction: f64,
        temp_int_air: f64,
        temp_ext_air: f64,
        r_v_arg: f64,
        r_w_arg_min_max: Option<f64>,
        reporting_flag: Option<ReportingFlag>,
        report_effective_flow_rate: Option<bool>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let report_effective_flow_rate = report_effective_flow_rate.unwrap_or(false);
        let (mut qm_in, _, qm_effective_flow_rate) = self
            .implicit_mass_balance_for_internal_reference_pressure_components(
                p_z_ref,
                wind_speed,
                wind_direction,
                temp_int_air,
                temp_ext_air,
                r_v_arg,
                r_w_arg_min_max,
                reporting_flag,
                simtime,
            )?;

        if report_effective_flow_rate {
            qm_in -= qm_effective_flow_rate
        }
        Ok(convert_mass_flow_rate_to_volume_flow_rate(
            qm_in,
            celsius_to_kelvin(temp_ext_air)?,
            self.p_a_alt,
        ))
    }

    /// Implicit mass balance for calculation of the internal reference pressure
    /// Equation 67 from BS EN 16798-7.
    ///
    /// Arguments:
    /// * `p_z_ref` - internal reference pressure (Pa)
    /// * `wind_speed` - wind speed, in m/s
    /// * `wind_direction` - direction wind is blowing from, in clockwise degrees from North
    /// * `temp_int_air` - temperature of air in the zone (C)
    /// * `temp_ext_air` - temperature of external air (C)
    /// * `reporting_flag` - flag used to give more detailed ventilation outputs (None = no additional reporting)
    ///
    /// Key Variables:
    /// qm_SUP_to_vent_zone - Supply air mass flow rate going to ventilation zone
    /// qm_ETA_from_vent_zone - Extract air mass flow rate from a ventilation zone
    /// qm_in_through_comb - Air mass flow rate entering through combustion appliances
    /// qm_out_through_comb - Air mass flow rate leaving through combustion appliances
    /// qm_in_through_passive_hybrid_ducts - Air mass flow rate entering through passive or hybrid duct
    /// qm_out_through_passive_hybrid_ducts - Air mass flow rate leaving through passive or hybrid duct
    /// qm_in_through_window_opening - Air mass flow rate entering through window opening
    /// qm_out_through_window_opening - Air mass flow rate leaving through window opening
    /// qm_in_through_vents - Air mass flow rate entering through vents (openings in the external envelope)
    /// qm_out_through_vents - Air mass flow rate leaving through vents (openings in the external envelope)
    /// qm_in_through_leaks - Air mass flow rate entering through envelope leakage
    /// qm_out_through_leaks - Air mass flow rate leaving through envelope leakage
    fn implicit_mass_balance_for_internal_reference_pressure_components(
        &self,
        p_z_ref: f64,
        wind_speed: f64,
        wind_direction: f64,
        temp_int_air: f64,
        temp_ext_air: f64,
        r_v_arg: f64,
        r_w_arg_min_max: Option<f64>,
        reporting_flag: Option<ReportingFlag>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64, f64)> {
        let u_site = wind_speed_at_zone_level(self.c_rgh_site, wind_speed, None, None, None);
        let t_e = celsius_to_kelvin(temp_ext_air)?;
        let t_z = celsius_to_kelvin(temp_int_air)?;
        let mut qm_in_through_window_opening = 0.;
        let mut qm_out_through_window_opening = 0.;
        let mut qm_in_through_vents = 0.;
        let mut qm_out_through_vents = 0.;
        let mut qm_in_through_leaks = 0.;
        let mut qm_out_through_leaks = 0.;
        let mut qm_in_through_comb = 0.;
        let mut qm_out_through_comb = 0.;
        let mut qm_in_through_passive_hybrid_ducts = 0.;
        let mut qm_out_through_passive_hybrid_ducts = 0.;
        let mut qm_sup_to_vent_zone = 0.;
        let mut qm_eta_from_vent_zone = 0.;
        let mut qm_in_effective_heat_recovery_saving_total = 0.0;

        for window in &self.windows {
            let (qm_in, qm_out) = window.calculate_flow_from_internal_p(
                wind_direction,
                u_site,
                t_e,
                t_z,
                p_z_ref,
                self.f_cross,
                self.shield_class,
                r_w_arg_min_max,
                simtime,
            );
            qm_in_through_window_opening += qm_in;
            qm_out_through_window_opening += qm_out;
        }

        for vent in &self.vents {
            let (qm_in, qm_out) = vent.calculate_flow_from_internal_p(
                wind_direction,
                u_site,
                t_e,
                t_z,
                p_z_ref,
                self.f_cross,
                self.shield_class,
                r_v_arg,
            );
            qm_in_through_vents += qm_in;
            qm_out_through_vents += qm_out;
        }

        for leak in &self.leaks {
            let (qm_in, qm_out) = leak.calculate_flow_from_internal_p(
                u_site,
                t_e,
                t_z,
                p_z_ref,
                self.f_cross,
                self.shield_class,
            );
            qm_in_through_leaks += qm_in;
            qm_out_through_leaks += qm_out;
        }

        for _atd in &self.air_terminal_devices {
            let qv_pdu_initial = 0.; // TODO (from Python) get from prev timestep
            let h_z = self.ventilation_zone_height;
            let qv_pdu = self.calculate_qv_pdu(qv_pdu_initial, p_z_ref, t_z, t_e, h_z)?;

            let (qv_pdu_in, qv_pdu_out) = if qv_pdu >= 0. {
                (qv_pdu, 0.)
            } else {
                (0., qv_pdu)
            };

            let (qm_in_through_phds, qm_out_through_phds) =
                convert_to_mass_air_flow_rate(qv_pdu_in, qv_pdu_out, t_e, t_z, self.p_a_alt);

            qm_in_through_passive_hybrid_ducts += qm_in_through_phds;
            qm_out_through_passive_hybrid_ducts += qm_out_through_phds;
        }

        for combustion_appliance in &self.combustion_appliances {
            let p_h_fi = 0.; // TODO (from Python) to work out from previous zone temperature? - Combustion appliance heating fuel input power
            let f_op_comb = 1.; // TODO (from Python) work out what turns the appliance on or off. Schedule or Logic?
            let (qv_in, qv_out) =
                combustion_appliance.calculate_air_flow_req_for_comb_appliance(f_op_comb, p_h_fi);
            let (qm_in_comb, qm_out_comb) =
                convert_to_mass_air_flow_rate(qv_in, qv_out, t_e, t_z, self.p_a_alt);
            qm_in_through_comb += qm_in_comb;
            qm_out_through_comb += qm_out_comb;
        }

        for mech_vent in &self.mech_vents {
            let (qm_sup, qm_eta, qm_in_effective_heat_recovery_saving) = mech_vent
                .calc_mech_vent_air_flw_rates_req_to_supply_vent_zone(t_z, t_e, &simtime)?;
            qm_sup_to_vent_zone += qm_sup;
            qm_eta_from_vent_zone += qm_eta;
            qm_in_effective_heat_recovery_saving_total += qm_in_effective_heat_recovery_saving;
        }

        let qm_in = qm_in_through_window_opening
            + qm_in_through_vents
            + qm_in_through_leaks
            + qm_in_through_comb
            + qm_in_through_passive_hybrid_ducts
            + qm_sup_to_vent_zone;

        let qm_out = qm_out_through_window_opening
            + qm_out_through_vents
            + qm_out_through_leaks
            + qm_out_through_comb
            + qm_out_through_passive_hybrid_ducts
            + qm_eta_from_vent_zone;

        // Output detailed ventilation file
        if self.detailed_output_heating_cooling {
            if let Some(reporting_flag) = reporting_flag {
                let incoming_air_flow =
                    convert_mass_flow_rate_to_volume_flow_rate(qm_in, t_e, self.p_a_alt);
                let air_changes_per_hour = incoming_air_flow / self.total_volume;

                self.detailed_results
                    .write()
                    .push(VentilationDetailedResult {
                        timestep_index: simtime.index,
                        reporting_flag,
                        r_v_arg,
                        incoming_air_flow,
                        total_volume: self.total_volume,
                        air_changes_per_hour,
                        temp_int_air,
                        p_z_ref,
                        qm_in_through_window_opening,
                        qm_out_through_window_opening,
                        qm_in_through_vents,
                        qm_out_through_vents,
                        qm_in_through_leaks,
                        qm_out_through_leaks,
                        qm_in_through_comb,
                        qm_out_through_comb,
                        qm_in_through_passive_hybrid_ducts,
                        qm_out_through_passive_hybrid_ducts,
                        qm_sup_to_vent_zone,
                        qm_eta_from_vent_zone,
                        qm_in_effective_heat_recovery_saving_total,
                        qm_in,
                        qm_out,
                    });
            }
        }

        Ok((qm_in, qm_out, qm_in_effective_heat_recovery_saving_total))
    }

    pub(crate) fn output_vent_results(&self) -> Arc<RwLock<Vec<VentilationDetailedResult>>> {
        Arc::clone(&self.detailed_results)
    }

    pub(crate) fn calc_air_changes_per_hour(
        &self,
        wind_speed: f64,
        wind_direction: f64,
        temp_int_air: f64,
        temp_ext_air: f64,
        r_v_arg: f64,
        r_w_arg: Option<f64>,
        initial_p_z_ref_guess: f64,
        reporting_flag: Option<ReportingFlag>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        #[cfg(test)]
        {
            // if there is an override implementation for calculating air changes per hour, use that
            if let Some(calc_ach_fn) = self.calc_air_changes_fn {
                return calc_ach_fn(
                    self,
                    wind_speed,
                    wind_direction,
                    temp_int_air,
                    temp_ext_air,
                    r_v_arg,
                    r_w_arg,
                    initial_p_z_ref_guess,
                    reporting_flag,
                    simtime,
                );
            }
        }

        let internal_reference_pressure = self.calculate_internal_reference_pressure(
            initial_p_z_ref_guess,
            wind_speed,
            wind_direction,
            temp_int_air,
            temp_ext_air,
            r_v_arg,
            r_w_arg,
            simtime,
        )?;

        let incoming_air_flow = self.incoming_air_flow(
            internal_reference_pressure,
            wind_speed,
            wind_direction,
            temp_int_air,
            temp_ext_air,
            r_v_arg,
            r_w_arg,
            reporting_flag,
            Some(true),
            simtime,
        )?;

        Ok(incoming_air_flow / self.total_volume)
    }

    /// Calculates the difference between the target air changes per hour (ACH) and the current ACH.
    ///
    /// Arguments:
    /// * `r_v_arg` - Current vent position, where 0 means vents are fully closed and 1 means vents are fully open.
    /// * `wind_speed` - Speed of the wind.
    /// * `wind_direction` - Direction of the wind.
    /// * `temp_int_air` - Interior air temperature.
    /// * `temp_ext_air` - Exterior air temperature.
    /// * `ach_target` - The desired target ACH value that needs to be achieved.
    /// * `r_w_arg` - Parameter related to the wind or building ventilation.
    /// * `initial_p_z_ref_guess` -Initial guess for reference pressure.
    /// * `reporting_flag` - Flag indicating whether to report detailed output
    ///
    /// Returns:
    ///     The adjusted absolute difference between the calculated ACH and the target ACH.
    ///         The difference is rounded to the 10th decimal place and a small gradient adjustment is applied
    ///         to help avoid numerical issues and local minima.
    fn calc_diff_ach_target(
        &self,
        r_v_arg: f64,
        wind_speed: f64,
        wind_direction: f64,
        temp_int_air: f64,
        temp_ext_air: f64,
        ach_target: f64,
        r_w_arg: Option<f64>,
        initial_p_z_ref_guess: f64,
        reporting_flag: Option<ReportingFlag>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let ach = self.calc_air_changes_per_hour(
            wind_speed,
            wind_direction,
            temp_int_air,
            temp_ext_air,
            r_v_arg,
            r_w_arg,
            initial_p_z_ref_guess,
            reporting_flag,
            simtime,
        )?;

        // To avoid the solver finding local minimums & numerically stagnating:
        // 1. The residuals of this function (ach- ach_target) are rounded to the 10th decimal place
        // 2. A very small gradient (1e-10 * R_v_arg) is added to the residuals to slightly 'tilt'
        //    the surface of the function towards higher R_v_arg values. This is because it can be flat
        //    at low R_v_arg values when flow is dominated by other components.
        Ok(((ach - ach_target) * 1e10).round() / 1e10 - 1e-10 * r_v_arg)
    }

    /// Determines the optimal vent position (R_v_arg) to achieve a desired air
    /// changes per hour (ACH) within specified bounds.
    ///
    /// Arguments:
    /// * `ach_min` - Minimum ACH limit.
    /// * `ach_max` - Maximum ACH limit.
    /// * `initial_r_v_arg` - Initial vent position, 0 = vents closed and 1 = vents fully open.
    /// * `wind_speed` - Speed of the wind.
    /// * `wind_direction` - Direction of the wind.
    /// * `temp_int_air` - Interior air temperature.
    /// * `temp_ext_air` - Exterior air temperature.
    /// * `r_w_arg` - Parameter related to the wind or building ventilation.
    /// * `initial_p_z_ref_guess` - Initial guess for reference pressure.
    /// * `reporting_flag` - Flag indicating whether to report detailed output.
    /// * `simtime`
    ///
    /// Returns:
    ///     The optimal vent position (R_v_arg) that brings the ACH within the specified bounds.
    pub(crate) fn find_r_v_arg_within_bounds(
        &self,
        ach_min: Option<f64>,
        ach_max: Option<f64>,
        initial_r_v_arg: f64,
        wind_speed: f64,
        wind_direction: f64,
        temp_int_air: f64,
        temp_ext_air: f64,
        r_w_arg: Option<f64>,
        initial_p_z_ref_guess: f64,
        reporting_flag: Option<ReportingFlag>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let initial_ach = self.calc_air_changes_per_hour(
            wind_speed,
            wind_direction,
            temp_int_air,
            temp_ext_air,
            initial_r_v_arg,
            r_w_arg,
            initial_p_z_ref_guess,
            reporting_flag,
            simtime,
        )?;

        // Determine if initial_ach is within the bounds
        if let (Some(ach_min), Some(ach_max)) = (ach_min, ach_max) {
            if ach_min > ach_max {
                bail!("ach_min must be less than ach_max");
            }
            if ach_min <= initial_ach && initial_ach <= ach_max {
                return Ok(initial_r_v_arg);
            }
        }

        let mut ach_target: Option<f64> = None;
        // Check extremes with fully open or closed vents
        match ach_min {
            Some(ach_min) if initial_ach < ach_min => {
                // If initial ACH is less than ach_min, check ach with vents fully open
                let ach_vent_open = self.calc_air_changes_per_hour(
                    wind_speed,
                    wind_direction,
                    temp_int_air,
                    temp_ext_air,
                    1., // vents fully open
                    r_w_arg,
                    initial_p_z_ref_guess,
                    reporting_flag,
                    simtime,
                )?;
                if ach_vent_open < ach_min {
                    // If the maximum achievable ACH with vents fully open is less than ach_min
                    return Ok(1.0);
                }

                // If current ACH is too low but ACH with vents fully open is higher than the threshold, set ach_target to ach_min
                ach_target.replace(ach_min);
            }
            _ => {}
        }

        match ach_max {
            Some(ach_max) if initial_ach > ach_max => {
                let ach_vent_closed = self.calc_air_changes_per_hour(
                    wind_speed,
                    wind_direction,
                    temp_int_air,
                    temp_ext_air,
                    0., // vents fully closed
                    r_w_arg,
                    initial_p_z_ref_guess,
                    reporting_flag,
                    simtime,
                )?;
                if ach_vent_closed > ach_max {
                    // If the minimum achievable ACH with vents fully closed is less than ach_max
                    return Ok(0.0);
                }

                // If current ACH is too high but ACH with vents fully closed is lower than the threshold, set ach_target to ach_max
                ach_target.replace(ach_max);
            }
            _ => {}
        }

        // if ach_target is still None, no need for further adjustment
        let ach_target = match ach_target {
            Some(ach_target) => ach_target,
            None => return Ok(initial_r_v_arg),
        };

        // With ach_target set to either ach_min or ach_max, run a Brent search (as equivalent of the minimize_scalar solver in Python)
        let cost = FindRVArgProblem {
            infiltration_ventilation: self,
            wind_speed,
            wind_direction,
            temp_int_air,
            temp_ext_air,
            ach_target,
            r_w_arg,
            initial_p_z_ref_guess,
            reporting_flag,
            simtime,
        };
        let solver = BrentRoot::new(0., 1., 1e-10);

        let optimization = Executor::new(cost, solver).run()?;

        optimization
            .state()
            .best_param
            .ok_or_else(|| anyhow!("No best param available in solver result"))
    }

    /// Equivalent of create_infiltration_ventilation in upstream
    pub(crate) fn create(
        input: &InfiltrationVentilationInput,
        zones: &ZoneDictionary,
        detailed_output_heating_cooling: bool,
        energy_supplies: &IndexMap<String, Arc<RwLock<EnergySupply>>>,
        controls: &Controls,
    ) -> anyhow::Result<Self> {
        let ventilation_zone_base_height = input.ventilation_zone_base_height;

        let windows = zones
            .values()
            .flat_map(|zone| zone.building_elements.values())
            .map(|building_element| {
                anyhow::Ok(if let BuildingElement::Transparent {
                    window_openable_control,
                    free_area_height,
                    mid_height,
                    max_window_open_area,
                    window_part_list,
                    orientation,
                    pitch,
                    ..
                } = building_element
                {
                    Some({
                        let on_off_ctrl =
                            window_openable_control
                                .as_ref()
                                .and_then(|window_openable_control| {
                                    controls.get_with_string(window_openable_control)
                                });
                        anyhow::Ok(Window::new(free_area_height.ok_or_else(|| anyhow!("A free_area_height value was expected for a transparent building element."))?, mid_height.ok_or_else(|| anyhow!("A mid_height value was expected for a transparent building element."))?,
                                       max_window_open_area.ok_or_else(|| anyhow!("A max_window_open_area value was expected for a transparent building element."))?,
                                       window_part_list.as_ref().unwrap_or(&vec![]).clone(), init_orientation(*orientation),
                                       *pitch,
                                       input.altitude,
                                       on_off_ctrl, ventilation_zone_base_height))
                    })
                } else {
                    None
                })
            })
            .filter_map(|x| x.ok())
            .flatten()
            .try_collect()?;

        let (pitches, areas): (Vec<f64>, Vec<f64>) = zones
            .values()
            .flat_map(|zone| zone.building_elements.values())
            .flat_map(|building_element| match building_element {
                BuildingElement::Opaque { pitch, area, .. }
                    if pitch_class(*pitch) == HeatFlowDirection::Upwards =>
                {
                    Some((pitch, area))
                }
                _ => None,
            })
            .unzip();
        // Work out the average pitch, weighted by area.
        let area_total = areas.iter().sum::<f64>();
        let average_pitch = if !pitches.is_empty() {
            areas
                .iter()
                .map(|x| x / area_total)
                .zip(pitches.iter())
                .map(|(x, &y)| x * y)
                .sum::<f64>()
        } else {
            0.
        };

        let (surface_area_facades_list, surface_area_roof_list) = zones
            .values()
            .flat_map(|zone| zone.building_elements.values())
            .fold((vec![], vec![]), |(mut facades, mut roofs), item| {
                match pitch_class(item.pitch()) {
                    HeatFlowDirection::Horizontal => match item {
                        BuildingElement::Opaque { area, .. } => {
                            facades.push(*area);
                        }
                        BuildingElement::Transparent { height, width, .. } => {
                            facades.push(*height * *width);
                        }
                        _ => {}
                    },
                    HeatFlowDirection::Upwards => match item {
                        BuildingElement::Opaque { area, .. } => {
                            roofs.push(*area);
                        }
                        BuildingElement::Transparent { height, width, .. } => {
                            roofs.push(*height * *width)
                        }
                        _ => {}
                    },
                    _ => {}
                }

                (facades, roofs)
            });

        let surface_area_facades = surface_area_facades_list.iter().sum::<f64>();
        let surface_area_roof = surface_area_roof_list.iter().sum::<f64>();

        let total_volume = zones.values().map(|z| z.volume).sum::<f64>();

        let vents = input
            .vents
            .values()
            .map(|vent| {
                Vent::new(
                    vent.mid_height_air_flow_path,
                    vent.area_cm2,
                    vent.pressure_difference_ref,
                    // Python uses "orientation360" value here
                    init_orientation(vent.orientation),
                    vent.pitch,
                    input.altitude,
                    ventilation_zone_base_height,
                )
            })
            .collect();

        let leaks = CompletedVentilationLeaks::complete_input(
            input,
            surface_area_facades,
            surface_area_roof,
        );

        let combustion_appliances = input
            .combustion_appliances
            .values()
            .map(|combustion_appliances_data| {
                CombustionAppliances::new(
                    combustion_appliances_data.supply_situation,
                    combustion_appliances_data.exhaust_situation,
                    combustion_appliances_data.fuel_type,
                    combustion_appliances_data.appliance_type,
                )
            })
            .collect();

        // (unfinished in upstream Python)
        let atds = Default::default();

        let mut mechanical_ventilations: Vec<Arc<MechanicalVentilation>> = Default::default();
        let mut space_heating_ductwork: IndexMap<String, Vec<Ductwork>> = Default::default();

        for (mech_vents_name, mech_vents_data) in input.mechanical_ventilation.iter() {
            let ctrl_intermittent_mev = mech_vents_data
                .control
                .as_ref()
                .and_then(|ctrl_name| controls.get_with_string(ctrl_name));

            let energy_supply = energy_supplies
                .get(&mech_vents_data.energy_supply)
                .ok_or_else(|| {
                    anyhow!(
                    "The energy supply '{}' indicated for mechanical ventilation was not declared.",
                    mech_vents_data.energy_supply
                )
                })?;
            let energy_supply_connection =
                EnergySupply::connection(energy_supply.clone(), mech_vents_name)?;

            mechanical_ventilations.push(
                Arc::new(MechanicalVentilation::new(mech_vents_data.supply_air_flow_rate_control, mech_vents_data.supply_air_temperature_control_type, 0., 0., mech_vents_data.vent_type, mech_vents_data.sfp.ok_or_else(|| anyhow!("A specific fan power value is expected for a mechanical ventilation unit."))?, mech_vents_data.design_outdoor_air_flow_rate, energy_supply_connection, total_volume, input.altitude, ctrl_intermittent_mev, match mech_vents_data.vent_type {
                    VentType::Mvhr => mech_vents_data.mvhr_efficiency,
                    VentType::IntermittentMev
                    | VentType::CentralisedContinuousMev
                    | VentType::DecentralisedContinuousMev => {
                        None
                    }
                    VentType::Piv => bail!("PIV vent type is not currently recognised when building up mechanical ventilation values for calculation"),
                }, None)),
            );

            // TODO (from Python) not all dwellings have mech vents - update to make mech vents optional
            if mech_vents_data.vent_type == VentType::Mvhr {
                space_heating_ductwork.insert(
                    mech_vents_name.to_owned(),
                    mech_vents_data
                        .ductwork
                        .as_ref()
                        .iter()
                        .flat_map(|ductworks| {
                            ductworks.iter().map(|ductwork| -> anyhow::Result<Ductwork> {
                                let (duct_perimeter, internal_diameter, external_diameter) =
                                    match ductwork.cross_section_shape {
                                        DuctShape::Circular => (None, Some(ductwork.internal_diameter_mm.ok_or_else(|| anyhow!("Expected an internal diameter value for ductwork with a circular cross-section."))? / MILLIMETRES_IN_METRE as f64), Some(ductwork.external_diameter_mm.ok_or_else(|| anyhow!("Expected an internal diameter value for ductwork with a circular cross-section."))? / MILLIMETRES_IN_METRE as f64)),
                                        DuctShape::Rectangular => (Some(ductwork.duct_perimeter_mm.ok_or_else(|| anyhow!("Expected a duct perimeter value for ductwork with a rectangular cross-section."))?), None, None),
                                    };

                                Ductwork::new(ductwork.cross_section_shape, duct_perimeter, internal_diameter, external_diameter, ductwork.length, ductwork.insulation_thermal_conductivity, ductwork.insulation_thickness_mm, ductwork.reflective, ductwork.duct_type, mech_vents_data.mvhr_location.ok_or_else(|| anyhow!("An MVHR location was expected for mechanical ventilation with an MVHR vent type."))?, mech_vents_data.mvhr_efficiency.ok_or_else(|| anyhow!("An MVHR efficiency value was expected for mechanical ventilation with an MVHR vent type."))?)
                            })
                        })
                        .collect::<anyhow::Result<Vec<Ductwork>>>()?,
                );
            }
        }

        Ok(InfiltrationVentilation::new(
            input.cross_vent_factor,
            input.shield_class,
            &input.terrain_class,
            average_pitch,
            windows,
            vents,
            leaks,
            combustion_appliances,
            atds,
            mechanical_ventilations,
            Arc::new(space_heating_ductwork),
            detailed_output_heating_cooling,
            input.altitude,
            total_volume,
            ventilation_zone_base_height,
        ))
    }
}

struct FindRVArgProblem<'a> {
    infiltration_ventilation: &'a InfiltrationVentilation,
    wind_speed: f64,
    wind_direction: f64,
    temp_int_air: f64,
    temp_ext_air: f64,
    ach_target: f64,
    r_w_arg: Option<f64>,
    initial_p_z_ref_guess: f64,
    reporting_flag: Option<ReportingFlag>,
    simtime: SimulationTimeIteration,
}

impl CostFunction for FindRVArgProblem<'_> {
    type Param = f64;
    type Output = f64;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, Error> {
        InfiltrationVentilation::calc_diff_ach_target(
            self.infiltration_ventilation,
            *param,
            self.wind_speed,
            self.wind_direction,
            self.temp_int_air,
            self.temp_ext_air,
            self.ach_target,
            self.r_w_arg,
            self.initial_p_z_ref_guess,
            self.reporting_flag,
            self.simtime,
        )
    }
}

struct ImplicitMassBalanceProblem<'a> {
    wind_speed: f64,
    wind_direction: f64,
    temp_int_air: f64,
    temp_ext_air: f64,
    r_v_arg: f64,
    r_w_arg: Option<f64>,
    simtime: SimulationTimeIteration,
    infiltration_ventilation: &'a InfiltrationVentilation,
}

impl CostFunction for ImplicitMassBalanceProblem<'_> {
    type Param = f64;
    type Output = f64;

    fn cost(&self, p_z_ref: &Self::Param) -> Result<Self::Output, Error> {
        let cost = self
            .infiltration_ventilation
            .implicit_mass_balance_for_internal_reference_pressure(
                *p_z_ref,
                self.wind_speed,
                self.wind_direction,
                self.temp_int_air,
                self.temp_ext_air,
                self.r_v_arg,
                self.r_w_arg,
                None,
                self.simtime,
            )?;
        Ok(cost)
    }
}

fn root_scalar_for_implicit_mass_balance(
    infiltration_ventilation: &InfiltrationVentilation,
    wind_speed: f64,
    wind_direction: f64,
    temp_int_air: f64,
    temp_ext_air: f64,
    r_v_arg: f64,
    r_w_arg: Option<f64>,
    simtime: SimulationTimeIteration,
    bracket: (f64, f64),
) -> Result<f64, &'static str> {
    let problem = ImplicitMassBalanceProblem {
        wind_speed,
        wind_direction,
        temp_int_air,
        temp_ext_air,
        r_v_arg,
        r_w_arg,
        simtime,
        infiltration_ventilation,
    };

    let tol = 0.;

    let (min, max) = bracket;
    let solver = BrentRoot::new(min, max, tol);

    let executor = Executor::new(problem, solver);
    let res = executor.run();

    let best_p_z_ref = match res {
        Ok(res) => res.state().best_param,
        Err(_) => return Err("Error calculating root for implicit mass balance"),
    };

    match best_p_z_ref {
        Some(p_z_ref) => Ok(p_z_ref),
        None => Err("No best_param in result"),
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct VentilationDetailedResult {
    timestep_index: usize,
    reporting_flag: ReportingFlag,
    r_v_arg: f64,
    incoming_air_flow: f64,
    total_volume: f64,
    air_changes_per_hour: f64,
    temp_int_air: f64,
    p_z_ref: f64,
    qm_in_through_window_opening: f64,
    qm_out_through_window_opening: f64,
    qm_in_through_vents: f64,
    qm_out_through_vents: f64,
    qm_in_through_leaks: f64,
    qm_out_through_leaks: f64,
    qm_in_through_comb: f64,
    qm_out_through_comb: f64,
    qm_in_through_passive_hybrid_ducts: f64,
    qm_out_through_passive_hybrid_ducts: f64,
    qm_sup_to_vent_zone: f64,
    qm_eta_from_vent_zone: f64,
    qm_in_effective_heat_recovery_saving_total: f64,
    qm_in: f64,
    qm_out: f64,
}

impl VentilationDetailedResult {
    pub(crate) fn as_string_values(&self) -> Vec<String> {
        vec![
            self.timestep_index.to_string(),
            self.reporting_flag.to_string(),
            self.r_v_arg.to_string(),
            self.incoming_air_flow.to_string(),
            self.total_volume.to_string(),
            self.air_changes_per_hour.to_string(),
            self.temp_int_air.to_string(),
            self.p_z_ref.to_string(),
            self.qm_in_through_window_opening.to_string(),
            self.qm_out_through_window_opening.to_string(),
            self.qm_in_through_vents.to_string(),
            self.qm_out_through_vents.to_string(),
            self.qm_in_through_leaks.to_string(),
            self.qm_out_through_leaks.to_string(),
            self.qm_in_through_comb.to_string(),
            self.qm_out_through_comb.to_string(),
            self.qm_in_through_passive_hybrid_ducts.to_string(),
            self.qm_out_through_passive_hybrid_ducts.to_string(),
            self.qm_sup_to_vent_zone.to_string(),
            self.qm_eta_from_vent_zone.to_string(),
            self.qm_in_effective_heat_recovery_saving_total.to_string(),
            self.qm_in.to_string(),
            self.qm_out.to_string(),
        ]
    }
}

#[derive(Debug, Error)]
#[error("Could not resolve an internal reference pressure for infiltration ventilation. Initial p_z_ref_guess: {initial_p_z_ref_guess}, temp_int_air: {temp_int_air}, r_w_arg: {r_w_arg:?}"
)]
pub struct InternalReferencePressureCalculationError {
    initial_p_z_ref_guess: f64,
    temp_int_air: f64,
    r_w_arg: Option<f64>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::controls::time_control::OnOffTimeControl;
    use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyBuilder};

    use crate::external_conditions::{DaylightSavingsConfig, ExternalConditions, ShadingSegment};
    use crate::input::FuelType;
    use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
    use approx::assert_relative_eq;
    use parking_lot::lock_api::RwLock;
    use rstest::{fixture, rstest};

    const EIGHT_DECIMAL_PLACES: f64 = 1e-7;

    #[test]
    fn test_calculate_pressure_difference_at_an_airflow_path() {
        let h_path: f64 = 0.4;
        let c_p_path: f64 = 0.45;
        let u_site: f64 = 1.;
        let t_e: f64 = 294.95;
        let t_z: f64 = 299.15;
        let p_z_ref: f64 = 2.5;
        let result = calculate_pressure_difference_at_an_airflow_path(
            h_path, c_p_path, u_site, t_e, t_z, p_z_ref,
        );
        assert_relative_eq!(result, -2.2966793114, max_relative = EIGHT_DECIMAL_PLACES);
        // Use spreadsheet to find answer.
    }

    #[test]
    fn test_wind_speed_at_zone_level() {
        let c_rgh_site = 0.8;
        let u_10 = 10.;
        let result = wind_speed_at_zone_level(c_rgh_site, u_10, None, None, None);
        assert_eq!(result, 8.)
    }

    #[rstest]
    #[case(CombustionFuelType::Wood, CombustionApplianceType::OpenFireplace, 2.8)]
    #[case(CombustionFuelType::Gas, CombustionApplianceType::ClosedWithFan, 0.38)]
    #[case(
        CombustionFuelType::Gas,
        CombustionApplianceType::OpenGasFlueBalancer,
        0.78
    )]
    #[case(
        CombustionFuelType::Gas,
        CombustionApplianceType::OpenGasKitchenStove,
        3.35
    )]
    #[case(CombustionFuelType::Gas, CombustionApplianceType::OpenGasFire, 3.35)]
    #[case(CombustionFuelType::Oil, CombustionApplianceType::ClosedFire, 0.32)]
    #[case(CombustionFuelType::Coal, CombustionApplianceType::ClosedFire, 0.52)]
    fn test_get_fuel_flow_factor(
        #[case] fuel_type: CombustionFuelType,
        #[case] appliance_type: CombustionApplianceType,
        #[case] expected: f64,
    ) {
        assert_eq!(get_fuel_flow_factor(fuel_type, appliance_type), expected);
    }

    #[rstest]
    #[case(
        CombustionAirSupplySituation::Outside,
        FlueGasExhaustSituation::IntoRoom,
        0.
    )]
    #[case(
        CombustionAirSupplySituation::RoomAir,
        FlueGasExhaustSituation::IntoRoom,
        0.
    )]
    #[case(
        CombustionAirSupplySituation::RoomAir,
        FlueGasExhaustSituation::IntoSeparateDuct,
        1.
    )]
    fn test_get_appliance_system_factor(
        #[case] supply_situation: CombustionAirSupplySituation,
        #[case] exhaust_situation: FlueGasExhaustSituation,
        #[case] expected: f64,
    ) {
        assert_eq!(
            get_appliance_system_factor(supply_situation, exhaust_situation),
            expected
        );
    }

    #[test]
    #[should_panic]
    fn test_get_appliance_system_factor_with_invalid_combination() {
        get_appliance_system_factor(
            CombustionAirSupplySituation::RoomAir,
            FlueGasExhaustSituation::IntoMechVent,
        );
    }

    #[test]
    fn test_adjust_air_density_for_altitude() {
        let h_alt = 10.; // meters
        let expected = 1.2028621569154314; // Pa
        let result = adjust_air_density_for_altitude(h_alt);
        assert_relative_eq!(result, expected); // Use spreadsheet to find answer.
    }

    #[test]
    fn test_air_density_at_temp() {
        let temperature = 300.; // K
        let air_density_adjusted_for_alt = 1.2; // kg/m^3
        let expected = 1.1725999999999999; // kg/m^3
        let result = air_density_at_temp(temperature, air_density_adjusted_for_alt);
        assert_relative_eq!(result, expected);
    }

    #[test]
    fn test_convert_volume_flow_rate_to_mass_flow_rate() {
        let qv = 1000.; // m ^ 3 / h
        let temperature = 300.; // K
        let p_a_alt = p_a_ref();
        let expected = 1176.5086666666666; // kg / h
        let result = convert_volume_flow_rate_to_mass_flow_rate(qv, temperature, p_a_alt);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_convert_mass_flow_rate_to_volume_flow_rate() {
        let qm = 1200.; // kg / h
        let temperature = 300.; // K
        let p_a_alt = p_a_ref();
        let expected = 1019.9669870685186; // m ^ 3 / h
        let result = convert_mass_flow_rate_to_volume_flow_rate(qm, temperature, p_a_alt);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_convert_to_mass_air_flow_rate() {
        let qv_in = 30.; // m ^ 3 / h
        let qv_out = 40.; // m ^ 3 / h
        let t_e = 300.; // K
        let t_z = 295.; // K
        let p_a_alt = p_a_ref();
        let expected_qm_in = 35.29526; // kg / h
        let expected_qm_out = 47.85797966101694; // kg / h
        let (qm_in, qm_out) = convert_to_mass_air_flow_rate(qv_in, qv_out, t_e, t_z, p_a_alt);
        assert_relative_eq!(qm_in, expected_qm_in);
        assert_relative_eq!(qm_out, expected_qm_out);
    }

    #[test]
    fn test_ter_class_to_roughness_coeff() {
        let z = 2.5;
        assert_eq!(
            ter_class_to_roughness_coeff(&TerrainClass::OpenWater, z),
            0.9386483560365819
        );
        assert_eq!(
            ter_class_to_roughness_coeff(&TerrainClass::OpenField, z),
            0.8325850605880374
        );
        assert_eq!(
            ter_class_to_roughness_coeff(&TerrainClass::Suburban, z),
            0.7223511561212699
        );
        assert_eq!(
            ter_class_to_roughness_coeff(&TerrainClass::Urban, z),
            0.6654212933375474
        );
    }

    #[test]
    fn test_orientation_difference() {
        // test simple cases
        assert_eq!(orientation_difference(0., 90.), 90.);
        assert_eq!(orientation_difference(100., 90.), 10.);
        // test handling of out of range input
        // (see test_orientation_difference_with_out_of_range_input below)
        // test cases where shortest angle crosses North
        assert_eq!(orientation_difference(0., 310.), 50.);
        assert_eq!(orientation_difference(300., 10.), 70.);
    }

    #[rstest]
    #[case(0., 450.)]
    #[case(540., 180.)]
    #[case(90., -290.)]
    #[case(-90., 90.)]
    #[should_panic]
    fn test_orientation_difference_with_out_of_range_input(
        #[case] orientation_1: f64,
        #[case] orientation_2: f64,
    ) {
        orientation_difference(orientation_1, orientation_2);
    }

    #[test]
    fn test_get_facade_direction() {
        assert_eq!(
            get_facade_direction(true, 0., 5., 0.),
            FacadeDirection::Roof10
        );
        assert_eq!(
            get_facade_direction(true, 0., 20., 0.),
            FacadeDirection::Roof10_30
        );
        assert_eq!(
            get_facade_direction(true, 0., 45., 0.),
            FacadeDirection::Roof30
        );
        assert_eq!(
            get_facade_direction(true, 0., 70., 0.),
            FacadeDirection::Windward
        );
        assert_eq!(
            get_facade_direction(true, 180., 70., 0.),
            FacadeDirection::Leeward
        );
        assert_eq!(
            get_facade_direction(false, 0., 45., 0.),
            FacadeDirection::Roof
        );
        assert_eq!(
            get_facade_direction(false, 0., 70., 0.),
            FacadeDirection::Windward
        );
        assert_eq!(
            get_facade_direction(false, 180., 70., 0.),
            FacadeDirection::Leeward
        );
    }

    #[test]
    fn test_get_c_p_path() {
        assert_relative_eq!(
            get_c_p_path_from_pitch_and_orientation(
                true,
                VentilationShieldClass::Open,
                10.,
                0.,
                0.,
                70.
            ),
            0.50
        );
        assert_relative_eq!(
            get_c_p_path_from_pitch_and_orientation(
                true,
                VentilationShieldClass::Normal,
                10.,
                0.,
                0.,
                70.
            ),
            0.25
        );
        assert_relative_eq!(
            get_c_p_path_from_pitch_and_orientation(
                true,
                VentilationShieldClass::Shielded,
                10.,
                0.,
                0.,
                70.
            ),
            0.05
        );
        assert_relative_eq!(
            get_c_p_path_from_pitch_and_orientation(
                true,
                VentilationShieldClass::Open,
                30.,
                0.,
                0.,
                70.
            ),
            0.65
        );
        assert_relative_eq!(
            get_c_p_path_from_pitch_and_orientation(
                true,
                VentilationShieldClass::Normal,
                30.,
                0.,
                0.,
                70.
            ),
            0.45
        );
        assert_relative_eq!(
            get_c_p_path_from_pitch_and_orientation(
                true,
                VentilationShieldClass::Shielded,
                30.,
                0.,
                0.,
                70.
            ),
            0.25
        );
        assert_relative_eq!(
            get_c_p_path_from_pitch_and_orientation(
                true,
                VentilationShieldClass::Open,
                60.,
                0.,
                0.,
                70.
            ),
            0.80
        );
        assert_relative_eq!(
            get_c_p_path_from_pitch_and_orientation(
                false,
                VentilationShieldClass::Normal,
                10.,
                0.,
                0.,
                70.
            ),
            0.05
        );
        assert_relative_eq!(
            get_c_p_path_from_pitch_and_orientation(
                false,
                VentilationShieldClass::Normal,
                10.,
                0.,
                0.,
                45.
            ),
            0.00
        );
    }

    #[fixture]
    fn simulation_time_iterator() -> SimulationTimeIterator {
        SimulationTime::new(0.0, 2.0, 1.0).iter()
    }

    #[fixture]
    fn wind_speeds() -> Vec<f64> {
        vec![3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4]
    }

    #[fixture]
    fn wind_directions() -> Vec<f64> {
        vec![200., 220., 230., 240., 250., 260., 260., 270.]
    }

    #[fixture]
    fn air_temps() -> Vec<f64> {
        vec![0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0]
    }

    #[fixture]
    fn external_conditions(
        simulation_time_iterator: SimulationTimeIterator,
    ) -> Arc<ExternalConditions> {
        let wind_speeds = vec![3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4];
        let wind_directions = vec![200., 220., 230., 240., 250., 260., 260., 270.];
        let air_temps = vec![0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0];
        let diffuse_horizontal_radiations = vec![333., 610., 572., 420., 0., 10., 90., 275.];
        let direct_beam_radiations = vec![420., 750., 425., 500., 0., 40., 0., 388.];
        let shading_segments = vec![
            ShadingSegment {
                number: 1,
                start: 180.,
                end: 135.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 2,
                start: 135.,
                end: 90.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 3,
                start: 90.,
                end: 45.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 4,
                start: 45.,
                end: 0.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 5,
                start: 0.,
                end: -45.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 6,
                start: -45.,
                end: -90.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 7,
                start: -90.,
                end: -135.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 8,
                start: -135.,
                end: -180.,
                shading_objects: None,
                ..Default::default()
            },
        ];
        Arc::new(ExternalConditions::new(
            &simulation_time_iterator,
            air_temps,
            wind_speeds,
            wind_directions,
            diffuse_horizontal_radiations,
            direct_beam_radiations,
            vec![0.2; 8760],
            51.42,
            -0.75,
            0,
            0,
            None,
            1.0,
            Some(1),
            Some(DaylightSavingsConfig::NotApplicable),
            false,
            false,
            shading_segments,
        ))
    }

    fn create_window(ctrl: Control, altitude: f64) -> Window {
        Window::new(
            1.6,
            1.5,
            3.,
            vec![WindowPartInput {
                mid_height_air_flow_path: 1.5,
            }],
            0.,
            90.,
            altitude,
            Some(Arc::new(ctrl)),
            0.,
        )
    }

    pub fn ctrl_that_is_on(simulation_time_iterator: SimulationTimeIterator) -> Control {
        Control::OnOffTime(OnOffTimeControl::new(
            vec![Some(true)],
            simulation_time_iterator.current_day(),
            1.,
        ))
    }

    pub fn ctrl_that_is_off(simulation_time_iterator: SimulationTimeIterator) -> Control {
        Control::OnOffTime(OnOffTimeControl::new(
            vec![Some(false)],
            simulation_time_iterator.current_day(),
            1.,
        ))
    }

    #[rstest]
    fn test_calculate_window_opening_free_area_ctrl_off(
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let ctrl = ctrl_that_is_off(simulation_time_iterator.clone());
        let window = create_window(ctrl, 0.);
        assert_eq!(
            window.calculate_window_opening_free_area(
                0.5,
                simulation_time_iterator.current_iteration()
            ),
            0.
        )
    }

    #[rstest]
    fn test_calculate_window_opening_free_area_ctrl_on(
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let ctrl = ctrl_that_is_on(simulation_time_iterator.clone());
        let window = create_window(ctrl, 0.);
        assert_eq!(
            window.calculate_window_opening_free_area(
                0.5,
                simulation_time_iterator.current_iteration()
            ),
            1.5
        )
    }

    #[rstest]
    fn test_calculate_flow_coeff_for_window_ctrl_off(
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let ctrl = ctrl_that_is_off(simulation_time_iterator.clone());
        let window = create_window(ctrl, 0.);
        assert_relative_eq!(
            window
                .calculate_flow_coeff_for_window(0.5, simulation_time_iterator.current_iteration()),
            0.
        )
    }

    #[rstest]
    fn test_calculate_flow_coeff_for_window_ctrl_on(
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let ctrl = ctrl_that_is_on(simulation_time_iterator.clone());
        let window = create_window(ctrl, 0.);
        let expected_a_w = 1.5;
        let expected_flow_coeff =
            3600. * window.c_d_w * expected_a_w * (2. / p_a_ref()).powf(window.n_w);
        assert_relative_eq!(
            window
                .calculate_flow_coeff_for_window(0.5, simulation_time_iterator.current_iteration()),
            expected_flow_coeff
        )
    }

    #[rstest]
    fn test_calculate_flow_from_internal_p(
        air_temps: Vec<f64>,
        wind_directions: Vec<f64>,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let u_site = 5.0;
        let _p_a_alt = p_a_ref();
        let _t_e = 290.15;
        let t_z = 293.15;
        let p_z_ref = 1.;
        let f_cross = true;
        let shield_class = VentilationShieldClass::Open;
        let r_w_arg = 0.5;
        let ctrl = ctrl_that_is_on(simulation_time_iterator.clone());
        let window = create_window(ctrl, 0.);

        let (qm_in, qm_out) = window.calculate_flow_from_internal_p(
            wind_directions[0],
            u_site,
            celsius_to_kelvin(air_temps[0]).unwrap(),
            t_z,
            p_z_ref,
            f_cross,
            shield_class,
            Some(r_w_arg),
            simulation_time_iterator.current_iteration(),
        );

        // qm_in returns 0.0 and qm_out returns -20707.309683335046
        assert_relative_eq!(qm_in, 0.);
        assert_relative_eq!(
            qm_out,
            -20707.309683335,
            max_relative = EIGHT_DECIMAL_PLACES
        );
    }

    #[fixture]
    fn window_part() -> WindowPart {
        WindowPart::new(1., 1.6, 0., 1, 0.)
    }

    #[rstest]
    fn test_calculate_ventilation_through_windows_using_internal_p(window_part: WindowPart) {
        let u_site = 3.7;
        let t_e = 273.15;
        let t_z = 293.15;
        let c_w_path = 4663.05;
        let c_p_path = -0.7;
        let p_z_ref = 1.;
        let expected_output = -13235.33116157;

        assert_relative_eq!(
            window_part.calculate_ventilation_through_windows_using_internal_p(
                u_site, t_e, t_z, c_w_path, p_z_ref, c_p_path
            ),
            expected_output,
            max_relative = EIGHT_DECIMAL_PLACES
        );
    }

    #[test]
    fn test_calculate_height_for_delta_p_w_div_path() {
        let expected_output = 1.;
        assert_relative_eq!(
            WindowPart::calculate_height_for_delta_p_w_div_path(1., 1.6, 0., 1usize),
            expected_output
        );
    }

    #[fixture]
    fn vent() -> Vent {
        Vent::new(1., 100., 20., 0., 90., 0., 0.)
    }

    #[rstest]
    fn test_calculate_vent_opening_free_area(vent: Vent) {
        let r_v_arg = 0.5;
        let expected_output = 50.;
        assert_eq!(
            vent.calculate_vent_opening_free_area(r_v_arg),
            expected_output,
        );
    }

    #[rstest]
    fn test_calculate_flow_coeff_for_vent(vent: Vent) {
        let r_v_arg = 1.;
        let expected_output = 27.8391201602292;
        assert_relative_eq!(
            vent.calculate_flow_coeff_for_vent(r_v_arg),
            expected_output,
            max_relative = EIGHT_DECIMAL_PLACES
        );
    }

    #[rstest]
    fn test_calculate_ventilation_through_vents_using_internal_p(vent: Vent) {
        let u_site = 3.7;
        let t_e = 273.15;
        let t_z = 293.15;
        let c_vent_path = 27.8391201602292;
        let c_p_path = -0.7;
        let p_z_ref = 1.;
        let expected_output = -79.01694696980;

        assert_relative_eq!(
            vent.calculate_ventilation_through_vents_using_internal_p(
                u_site,
                t_e,
                t_z,
                c_vent_path,
                c_p_path,
                p_z_ref
            ),
            expected_output,
            max_relative = EIGHT_DECIMAL_PLACES
        );
    }

    #[rstest]
    // in Python this is test_calculate_flow_from_internal_p
    fn test_calculate_flow_from_internal_p_for_vents(
        vent: Vent,
        wind_directions: Vec<f64>,
        air_temps: Vec<f64>,
    ) {
        let u_site = 3.7;
        let t_z = 293.15;
        let p_z_ref = 1.;
        let f_cross = true;
        let shield_class = VentilationShieldClass::Open;
        let r_v_arg = 1.;

        let (qm_in_through_vent, qm_out_through_vent) = vent.calculate_flow_from_internal_p(
            wind_directions[0],
            u_site,
            celsius_to_kelvin(air_temps[0]).unwrap(),
            t_z,
            p_z_ref,
            f_cross,
            shield_class,
            r_v_arg,
        );

        assert_relative_eq!(qm_in_through_vent, 0.);
        assert_relative_eq!(
            qm_out_through_vent,
            -95.136404151646,
            max_relative = EIGHT_DECIMAL_PLACES
        );
    }

    #[fixture]
    fn leaks() -> Leaks {
        Leaks::new(
            1.,
            50.,
            1.2,
            FacadeDirection::Leeward,
            100.,
            120.,
            220.,
            0.,
            0.,
        )
    }

    #[rstest]
    fn test_calculate_flow_coeff_for_leak(leaks: Leaks) {
        let expected_result = 2.6490460494125543;
        assert_relative_eq!(leaks.calculate_flow_coeff_for_leak(), expected_result)
    }

    #[rstest]
    fn test_calculate_ventilation_through_leaks_using_internal_p(leaks: Leaks) {
        let u_site = 3.7;
        let t_e = 273.15;
        let t_z = 293.15;
        let c_p_path = -0.7;
        let p_z_ref = 1.;
        let expected_output = -10.6531458050959;

        assert_relative_eq!(
            leaks.calculate_ventilation_through_leaks_using_internal_p(
                u_site, t_e, t_z, c_p_path, p_z_ref
            ),
            expected_output,
            max_relative = EIGHT_DECIMAL_PLACES
        );
    }

    #[rstest]
    // in Python this test is named test_calculate_flow_from_internal_p
    fn test_calculate_flow_from_internal_p_for_leaks(leaks: Leaks, air_temps: Vec<f64>) {
        let u_site = 3.7;
        let t_z = 293.15;
        let p_z_ref = 1.;
        let f_cross = true;
        let shield_class = VentilationShieldClass::Open;

        let (qm_in_through_leaks, qm_out_through_leaks) = leaks.calculate_flow_from_internal_p(
            u_site,
            celsius_to_kelvin(air_temps[0]).unwrap(),
            t_z,
            p_z_ref,
            f_cross,
            shield_class,
        );

        assert_relative_eq!(qm_in_through_leaks, 0.);
        assert_relative_eq!(qm_out_through_leaks, -12.826387549335472);
    }

    #[fixture]
    fn combustion_appliances() -> CombustionAppliances {
        CombustionAppliances::new(
            CombustionAirSupplySituation::RoomAir,
            FlueGasExhaustSituation::IntoSeparateDuct,
            CombustionFuelType::Wood,
            CombustionApplianceType::OpenFireplace,
        )
    }

    #[rstest]
    fn test_calculate_air_flow_req_for_comb_appliance(combustion_appliances: CombustionAppliances) {
        let f_op_comp = 1.;
        let p_h_fi = 1.;
        let (q_in_comb, q_out_comb) =
            combustion_appliances.calculate_air_flow_req_for_comb_appliance(f_op_comp, p_h_fi);

        assert_relative_eq!(q_in_comb, 0.);
        assert_relative_eq!(q_out_comb, -10.08);
    }

    #[fixture]
    pub fn energy_supply(simulation_time_iterator: SimulationTimeIterator) -> EnergySupply {
        EnergySupplyBuilder::new(
            FuelType::Electricity,
            simulation_time_iterator.total_steps(),
        )
        .build()
    }
    #[fixture]
    fn mechanical_ventilation(energy_supply: EnergySupply) -> MechanicalVentilation {
        let energy_supply = Arc::new(RwLock::new(energy_supply));
        let energy_supply_connection =
            EnergySupply::connection(energy_supply.clone(), "mech_vent_fans").unwrap();

        MechanicalVentilation::new(
            SupplyAirFlowRateControlType::ODA,
            SupplyAirTemperatureControlType::Constant,
            1.,
            3.4,
            VentType::Mvhr,
            1.5,
            0.5,
            energy_supply_connection,
            250.,
            0.,
            None,
            Some(0.),
            None,
        )
    }

    #[rstest]
    // In Python this tests calls 'calculate_required_outdoor_air_flow_rate' in the assertion,
    // we've implemented the 'new' function on MechanicalVentilation so that it sets
    // qv_oda_req_design by calling 'calculate_required_outdoor_air_flow_rate'
    fn test_calculate_required_outdoor_air_flow_rate(
        mechanical_ventilation: MechanicalVentilation,
    ) {
        let expected_result = 0.55;
        assert_relative_eq!(mechanical_ventilation.qv_oda_req_design, expected_result)
    }

    #[rstest]
    fn test_calc_req_oda_flow_rates_at_atds(mechanical_ventilation: MechanicalVentilation) {
        let (qv_sup_req, qv_eta_req) = mechanical_ventilation.calc_req_oda_flow_rates_at_atds();
        assert_relative_eq!(qv_sup_req, 0.55);
        assert_relative_eq!(qv_eta_req, -0.55);
    }

    #[rstest]
    fn test_calc_mech_vent_air_flw_rates_req_to_supply_vent_zone(
        mechanical_ventilation: MechanicalVentilation,
        air_temps: Vec<f64>,
        mut simulation_time_iterator: SimulationTimeIterator,
    ) {
        let (qm_sup_dis_req, qm_eta_dis_req, qm_in_effective_heat_recovery_saving) =
            mechanical_ventilation
                .calc_mech_vent_air_flw_rates_req_to_supply_vent_zone(
                    293.15,
                    celsius_to_kelvin(air_temps[0]).unwrap(),
                    &simulation_time_iterator.next().unwrap(),
                )
                .unwrap();
        assert_relative_eq!(qm_sup_dis_req, 0.7106861797547136);
        assert_relative_eq!(qm_eta_dis_req, -0.6622);
        assert_relative_eq!(qm_in_effective_heat_recovery_saving, 0.);
    }

    // NOTE - Python has a commented out test here

    #[fixture]
    fn infiltration_ventilation(
        simulation_time_iterator: SimulationTimeIterator,
        combustion_appliances: CombustionAppliances,
        mechanical_ventilation: MechanicalVentilation,
    ) -> InfiltrationVentilation {
        let ctrl = ctrl_that_is_on(simulation_time_iterator.clone());
        let windows = vec![create_window(ctrl, 30.)];
        let vents = vec![Vent::new(1.5, 100., 20., 0., 90., 30., 2.5)];
        let leaks = CompletedVentilationLeaks {
            ventilation_zone_height: 6.,
            test_pressure: 50.,
            test_result: 1.2,
            area_roof: 25.,
            area_facades: 85.,
            env_area: 220.,
            altitude: 30.,
        };
        let combustion_appliances_list = vec![combustion_appliances];
        let air_terminal_devices = Vec::<AirTerminalDevices>::new();
        let mechanical_ventilations = vec![Arc::new(mechanical_ventilation)];

        InfiltrationVentilation::new(
            true,
            VentilationShieldClass::Open,
            &TerrainClass::OpenField,
            20.0,
            windows,
            vents,
            leaks,
            combustion_appliances_list,
            air_terminal_devices,
            mechanical_ventilations,
            Default::default(),
            false,
            0.,
            250.,
            2.5,
        )
    }

    #[test]
    fn test_calculate_total_volume_air_flow_rate_in() {
        let qm_in = 0.5;
        let external_air_density = 1.;
        assert_relative_eq!(
            InfiltrationVentilation::calculate_total_volume_air_flow_rate_in(
                qm_in,
                external_air_density
            ),
            0.5
        );
    }

    #[test]
    fn test_calculate_total_volume_air_flow_rate_out() {
        let qm_out = 0.5;
        let zone_air_density = 1.;
        assert_relative_eq!(
            InfiltrationVentilation::calculate_total_volume_air_flow_rate_out(
                qm_out,
                zone_air_density
            ),
            0.5
        )
    }

    // Python has a make_leaks_object test here which isn't required for Rust

    // NOTE - Python has a commented out test here for test_calculate_qv_pdu
    // NOTE - Python has a commented out test here for test_implicit_formula_for_qv_pdu

    #[rstest]
    #[case(20.,  0.5, -6.235527862635629)]
    fn test_calculate_internal_reference_pressure(
        #[case] temp_int_air: f64,
        #[case] r_w_arg: f64,
        #[case] expected: f64,
        infiltration_ventilation: InfiltrationVentilation,
        wind_speeds: Vec<f64>,
        wind_directions: Vec<f64>,
        air_temps: Vec<f64>,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let initial_p_z_ref_guess = 0.;
        let r_v_arg = 1.;
        assert_relative_eq!(
            infiltration_ventilation
                .calculate_internal_reference_pressure(
                    initial_p_z_ref_guess,
                    wind_speeds[0],
                    wind_directions[0],
                    temp_int_air,
                    air_temps[0],
                    r_v_arg,
                    Some(r_w_arg),
                    simulation_time_iterator.current_iteration()
                )
                .unwrap(),
            expected,
            max_relative = EIGHT_DECIMAL_PLACES
        )
    }

    // TODO more cases here / different states
    // copy expected values across from Python
    #[rstest]
    fn test_implicit_mass_balance_for_internal_reference_pressure(
        infiltration_ventilation: InfiltrationVentilation,
        wind_speeds: Vec<f64>,
        wind_directions: Vec<f64>,
        air_temps: Vec<f64>,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let p_z_ref = 1.;
        let temp_int_air = 20.;
        let r_v_arg = 1.;
        let r_w_arg_min_max = 1.;
        assert_relative_eq!(
            infiltration_ventilation
                .implicit_mass_balance_for_internal_reference_pressure(
                    p_z_ref,
                    wind_speeds[0],
                    wind_directions[0],
                    temp_int_air,
                    air_temps[0],
                    r_v_arg,
                    Some(r_w_arg_min_max),
                    None,
                    simulation_time_iterator.current_iteration()
                )
                .unwrap(),
            -30270.984047975235
        )
    }

    #[rstest]
    fn test_incoming_air_flow(
        infiltration_ventilation: InfiltrationVentilation,
        wind_speeds: Vec<f64>,
        wind_directions: Vec<f64>,
        air_temps: Vec<f64>,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let p_z_ref = 1.;
        let temp_int_air = 20.;
        let r_v_arg = 1.;
        let r_w_arg_min_max = 1.;

        assert_relative_eq!(
            infiltration_ventilation
                .incoming_air_flow(
                    p_z_ref,
                    wind_speeds[0],
                    wind_directions[0],
                    temp_int_air,
                    air_temps[0],
                    r_v_arg,
                    Some(r_w_arg_min_max),
                    None,
                    None,
                    simulation_time_iterator.current_iteration()
                )
                .unwrap(),
            4.846594835429536
        )
    }

    #[rstest]
    fn test_find_r_v_arg_within_bounds(
        infiltration_ventilation: InfiltrationVentilation,
        air_temps: Vec<f64>,
        wind_directions: Vec<f64>,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        // Checking for ach_target = ach_max
        let ach_min = 0.3;
        let ach_max = 1.;
        let temp_int_air = 20.;
        let initial_r_v_arg = 1.;
        let expected_output = 0.5359731535118643;
        let actual_output = infiltration_ventilation
            .find_r_v_arg_within_bounds(
                Some(ach_min),
                Some(ach_max),
                initial_r_v_arg,
                20.,
                wind_directions[0],
                temp_int_air,
                air_temps[0],
                Some(0.),
                0.,
                None,
                simulation_time_iterator.current_iteration(),
            )
            .unwrap();
        assert_relative_eq!(
            actual_output,
            expected_output,
            max_relative = EIGHT_DECIMAL_PLACES
        );

        let ach_min = 1.0;
        let ach_max = 1.4;
        let temp_int_air = 20.;
        let initial_r_v_arg = 0.4;
        let expected_output = 0.5359731535118643;
        let actual_output = infiltration_ventilation
            .find_r_v_arg_within_bounds(
                Some(ach_min),
                Some(ach_max),
                initial_r_v_arg,
                20.,
                wind_directions[0],
                temp_int_air,
                air_temps[0],
                Some(0.),
                0.,
                None,
                simulation_time_iterator.current_iteration(),
            )
            .unwrap();
        assert_relative_eq!(
            actual_output,
            expected_output,
            max_relative = EIGHT_DECIMAL_PLACES
        );
    }

    #[fixture]
    fn infiltration_ventilation_with_patched_ach_fn(
        mut infiltration_ventilation: InfiltrationVentilation,
    ) -> InfiltrationVentilation {
        infiltration_ventilation.set_calc_air_changes_fn(|_, _, _, _, _, _, _, _, _, _| Ok(2.0));
        infiltration_ventilation
    }

    #[rstest]
    fn test_ach_within_bounds(
        infiltration_ventilation_with_patched_ach_fn: InfiltrationVentilation,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let result = infiltration_ventilation_with_patched_ach_fn
            .find_r_v_arg_within_bounds(
                Some(1.5),
                Some(2.5),
                0.5,
                5.0,
                90.0,
                20.0,
                10.0,
                Some(1.0),
                0.5,
                None,
                simulation_time_iterator.current_iteration(),
            )
            .unwrap();
        assert_eq!(result, 0.5);
    }

    #[rstest]
    fn test_no_ach_target(
        infiltration_ventilation_with_patched_ach_fn: InfiltrationVentilation,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        let result = infiltration_ventilation_with_patched_ach_fn
            .find_r_v_arg_within_bounds(
                None,
                None,
                0.5,
                5.0,
                90.0,
                20.0,
                10.0,
                Some(1.0),
                0.5,
                None,
                simulation_time_iterator.current_iteration(),
            )
            .unwrap();
        assert_eq!(result, 0.5);
    }
}
