// This module provides objects to represent Infiltration and Ventilation.
// The calculations are based on Method 1 of BS EN 16798-7.

use crate::compare_floats::max_of_2;
use crate::core::controls::time_control::{Control, ControlBehaviour};
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::material_properties::AIR;
use crate::core::units::{
    celsius_to_kelvin, LITRES_PER_CUBIC_METRE, SECONDS_PER_HOUR, WATTS_PER_KILOWATT,
};
use crate::external_conditions::ExternalConditions;
use crate::input::{
    CombustionAirSupplySituation, CombustionApplianceType, CombustionFuelType,
    FlueGasExhaustSituation, SupplyAirFlowRateControlType, SupplyAirTemperatureControlType,
    TerrainClass, VentType, VentilationLeaks, VentilationShieldClass,
    WindowPart as WindowPartInput,
};
use crate::simulation_time::SimulationTimeIteration;
use anyhow::Error;
use rand_distr::num_traits::abs;
use std::sync::Arc;

fn p_a_ref() -> f64 {
    AIR.density_kg_per_m3()
}
fn c_a() -> f64 {
    AIR.specific_heat_capacity_kwh()
}

// (Default values from BS EN 16798-7, Table 11)
// Coefficient to take into account stack effect in airing calculation in (m/s)/(m*K)
const C_STACK: f64 = 0.0035;
// Coefficient to take into account wind speed in airing calculation in 1/(m/s)
const C_WND: f64 = 0.001;
// Gravitational constant in m/s2
const G: f64 = 9.81;
//Room temperature in degrees K
const T_E_REF: f64 = 293.15;
//Absolute zero in degrees K
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
fn air_change_rate_to_flow_rate(air_change_rate: f64, zone_volume: f64) -> f64 {
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

/// Interpreted from Table B.13 in BS EN 16798-7.
/// Terrain Class input to roughness coefficient at building site at 10m
/// Arguments:
/// TER_CLASS -- Terrain class, one of 'Open terrain', 'Country' or 'Urban'
fn ter_class_to_roughness_coeff(terrain: TerrainClass) -> f64 {
    match terrain {
        TerrainClass::OpenTerrain => 1.0,
        TerrainClass::Country => 0.9,
        TerrainClass::Urban => 0.8,
    }
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

/// Determine orientation of other windows relative to largest
fn orientation_difference(orientation1: f64, orientation2: f64) -> f64 {
    let op_rel_orientation = abs(orientation1 - orientation2);

    if op_rel_orientation > 360. {
        return op_rel_orientation - 360.;
    }
    op_rel_orientation
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum FacadeDirection {
    Roof,
    Roof10,
    Roof10_30,
    Roof30,
    Windward,
    Leeward,
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
            if orientation_diff < 90. {
                FacadeDirection::Windward
            } else {
                FacadeDirection::Leeward
            }
        }
    } else if pitch < 60. {
        FacadeDirection::Roof
    } else {
        let orientation_diff = orientation_difference(orientation, wind_direction);
        if orientation_diff < 90. {
            FacadeDirection::Windward
        } else {
            FacadeDirection::Leeward
        }
    }
}

// we split the python get_c_p_path method into two methods below:
fn get_c_p_path_from_pitch_and_orientation(
    f_cross: bool,
    shield_class: VentilationShieldClass,
    h_path: f64,
    wind_direction: f64,
    orientation: f64,
    pitch: f64,
) -> f64 {
    let facade_direction = get_facade_direction(f_cross, orientation, pitch, wind_direction);
    get_c_p_path(f_cross, shield_class, h_path, facade_direction)
}

/// Interpreted from Table B.7 for determining dimensionless wind pressure coefficients
/// Arguments:
/// f_cross -- boolean, dependant on if cross ventilation is possible or not
/// shield_class -- indicates exposure to wind
/// h_path - height of flow path (m)
/// wind_direction -- direction the wind is blowing (degrees)
/// orientation -- orientation of the facade (degrees)
/// pitch -- pitch of the facade (degrees)
/// facade_direction -- direction of the facade (from get_facade_direction or manual entry)
fn get_c_p_path(
    f_cross: bool,
    shield_class: VentilationShieldClass,
    h_path: f64,
    facade_direction: FacadeDirection,
) -> f64 {
    if f_cross {
        if h_path < 15. {
            match shield_class {
                VentilationShieldClass::Open => match facade_direction {
                    FacadeDirection::Windward => 0.50,
                    FacadeDirection::Leeward => -0.70,
                    FacadeDirection::Roof10 => -0.70,
                    FacadeDirection::Roof10_30 => -0.60,
                    FacadeDirection::Roof30 => -0.20,
                    _ => panic!("Invalid combination of shield_class and facade_direction"),
                },
                VentilationShieldClass::Normal => match facade_direction {
                    FacadeDirection::Windward => 0.25,
                    FacadeDirection::Leeward => -0.50,
                    FacadeDirection::Roof10 => -0.60,
                    FacadeDirection::Roof10_30 => -0.50,
                    FacadeDirection::Roof30 => -0.20,
                    _ => panic!("Invalid combination of shield_class and facade_direction"),
                },
                VentilationShieldClass::Shielded => match facade_direction {
                    FacadeDirection::Windward => 0.05,
                    FacadeDirection::Leeward => -0.30,
                    FacadeDirection::Roof10 => -0.50,
                    FacadeDirection::Roof10_30 => -0.40,
                    FacadeDirection::Roof30 => -0.20,
                    _ => panic!("Invalid combination of shield_class and facade_direction"),
                },
            }
        } else if (15. ..50.).contains(&h_path) {
            match shield_class {
                VentilationShieldClass::Open => match facade_direction {
                    FacadeDirection::Windward => 0.65,
                    FacadeDirection::Leeward => -0.70,
                    FacadeDirection::Roof10 => -0.70,
                    FacadeDirection::Roof10_30 => -0.60,
                    FacadeDirection::Roof30 => -0.20,
                    _ => panic!("Invalid combination of shield_class and facade_direction"),
                },
                VentilationShieldClass::Normal => match facade_direction {
                    FacadeDirection::Windward => 0.45,
                    FacadeDirection::Leeward => -0.50,
                    FacadeDirection::Roof10 => -0.60,
                    FacadeDirection::Roof10_30 => -0.50,
                    FacadeDirection::Roof30 => -0.20,
                    _ => panic!("Invalid combination of shield_class and facade_direction"),
                },
                VentilationShieldClass::Shielded => match facade_direction {
                    FacadeDirection::Windward => 0.25,
                    FacadeDirection::Leeward => -0.30,
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

struct Window {
    h_w_fa: f64,
    h_w_path: f64,
    a_w_max: f64,
    c_d_w: f64,
    n_w: f64,
    orientation: f64,
    pitch: f64,
    external_conditions: Arc<ExternalConditions>,
    n_w_div: f64,
    on_off_ctrl_obj: Option<Arc<Control>>,
    _altitude: f64,
    p_a_alt: f64,
    window_parts: Vec<WindowPart>,
}

impl Window {
    fn new(
        external_conditions: Arc<ExternalConditions>,
        h_w_fa: f64,
        h_w_path: f64,
        a_w_max: f64,
        window_part_list: Vec<WindowPartInput>,
        orientation: f64,
        pitch: f64,
        altitude: f64,
        on_off_ctrl_obj: Option<Arc<Control>>,
    ) -> Self {
        let n_w_div = max_of_2(window_part_list.len() - 1, 0usize) as f64;
        Self {
            h_w_fa,
            h_w_path,
            a_w_max,
            c_d_w: 0.67,
            n_w: 0.5,
            orientation,
            pitch,
            external_conditions,
            n_w_div,
            on_off_ctrl_obj,
            _altitude: altitude,
            p_a_alt: adjust_air_density_for_altitude(altitude),
            window_parts: window_part_list
                .iter()
                .enumerate()
                .map(|(window_part_number, window_part_input)| {
                    WindowPart::new(
                        window_part_input.mid_height_air_flow_path,
                        h_w_fa,
                        n_w_div,
                        window_part_number + 1,
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

    /// Calculate the airflow through window opening based on the how open the window is and internal pressure
    /// Arguments:
    /// u_site -- wind velocity at zone level (m/s)
    /// T_z -- thermal zone air temperature (K)
    /// p_z_ref -- internal reference pressure (Pa)
    /// f_cross -- boolean, dependant on if cross ventilation is possible or not
    /// shield_class -- indicates exposure to wind
    /// R_w_arg -- ratio of window opening (0-1)
    fn calculate_flow_from_internal_p(
        &self,
        u_site: f64,
        t_z: f64,
        p_z_ref: f64,
        f_cross: bool,
        shield_class: VentilationShieldClass,
        r_w_arg: Option<f64>,
        simulation_time: SimulationTimeIteration,
    ) -> (f64, f64) {
        let wind_direction = self.external_conditions.wind_direction(simulation_time);
        let t_e = celsius_to_kelvin(self.external_conditions.air_temp(&simulation_time));
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
            self.h_w_path,
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

struct WindowPart {
    h_w_path: f64,
    h_w_fa: f64,
    n_w_div: f64,
    h_w_div_path: f64,
    n_w: f64,
}

impl WindowPart {
    fn new(h_w_path: f64, h_w_fa: f64, n_w_div: f64, window_part_number: usize) -> Self {
        Self {
            h_w_path,
            h_w_fa,
            n_w_div,
            h_w_div_path: Self::calculate_height_for_delta_p_w_div_path(
                h_w_path,
                h_w_fa,
                n_w_div,
                window_part_number,
            ),
            n_w: 0.5,
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
            * abs(delta_p_path).powf(self.n_w)
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

struct Vent {
    external_conditions: ExternalConditions,
    h_path: f64,
    a_vent: f64,
    delta_p_vent_ref: f64,
    orientation: f64,
    pitch: f64,
    altitude: f64,
    n_vent: f64,
    c_d_vent: f64,
    p_a_alt: f64,
    // NOTE - in Python we have C_vent_path as an instance variable but here we calculate it when needed instead
}

impl Vent {
    /// Construct a Vent object
    ///
    /// Arguments:
    ///    external_conditions -- reference to ExternalConditions object
    ///    h_path -- mid height of air flow path relative to ventilation zone (m)
    ///    A_vent - Equivalent area of a vent (m2)
    ///    delta_p_vent_ref -- reference pressure difference for vent (Pa)
    ///    orientation -- The orientation of the vent (degrees)
    ///    pitch -- The pitch of the vent (degrees)
    ///    altitude -- altitude of dwelling above sea level (m)
    ///
    /// Method:
    ///    - Based on Section 6.4.3.6 Airflow through vents from BS EN 16798-7
    fn new(
        external_conditions: ExternalConditions,
        h_path: f64,
        a_vent: f64,
        delta_p_vent_ref: f64,
        orientation: f64,
        pitch: f64,
        altitude: f64,
    ) -> Self {
        Self {
            h_path,
            a_vent,
            delta_p_vent_ref,
            orientation,
            pitch,
            altitude,
            external_conditions,
            n_vent: 0.5, // Flow exponent for vents based on Section B.3.2.2 from BS EN 16798-7
            c_d_vent: 0.6, // Discharge coefficient of vents based on B.3.2.1 from BS EN 16798-7
            p_a_alt: adjust_air_density_for_altitude(altitude),
        }
    }

    /// The airflow coefficient of the vent calculated from equivalent area A_vent_i
    /// according to EN 13141-1 and EN 13141-2.
    /// Based on Equation 59 from BS EN 16798-7.
    fn calculate_flow_coeff_for_vent(&self) -> f64 {
        // NOTE: The standard does not define what the below 3600 and 10000 are.

        (3600. / 10000.)
            * self.c_d_vent
            * self.a_vent
            * (2. / p_a_ref()).powf(0.5)
            * (1. / self.delta_p_vent_ref).powf(self.n_vent - 0.5)
    }

    /// Calculate the airflow through vents from internal pressure
    /// Arguments:
    /// u_site -- wind velocity at zone level (m/s)
    /// T_e -- external air temperature (K)
    /// T_z -- thermal zone air temperature (K)
    /// C_vent_path -- wind pressure coefficient at height of the vent
    /// C_p_path -- wind pressure coefficient at the height of the window part
    /// p_z_ref -- internal reference pressure (Pa)
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
        c_vent_path * sign(delta_p_path) as f64 * abs(delta_p_path).powf(self.n_vent)
    }

    /// Calculate the airflow through vents from internal pressure
    ///
    ///     Arguments:
    ///     u_site -- wind velocity at zone level (m/s)
    ///     T_z -- thermal zone air temperature (K)
    ///     p_z_ref -- internal reference pressure (Pa)
    ///     f_cross -- boolean, dependant on if cross ventilation is possible or not
    ///     shield_class -- indicates exposure to wind
    fn calculate_flow_from_internal_p(
        &self,
        u_site: f64,
        t_z: f64,
        p_z_ref: f64,
        f_cross: bool,
        shield_class: VentilationShieldClass,
        simulation_time: SimulationTimeIteration,
    ) -> (f64, f64) {
        let wind_direction = self.external_conditions.wind_direction(simulation_time);
        let t_e = celsius_to_kelvin(self.external_conditions.air_temp(&simulation_time));

        // Wind pressure coefficient for the air flow path
        let c_p_path = get_c_p_path_from_pitch_and_orientation(
            f_cross,
            shield_class,
            self.h_path,
            wind_direction,
            self.orientation,
            self.pitch,
        );

        // Calculate airflow through each vent

        let c_vent_path = self.calculate_flow_coeff_for_vent();
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
struct Leaks {
    h_path: f64,
    delta_p_leak_ref: f64,
    a_roof: f64,
    a_facades: f64,
    a_leak: f64,
    qv_delta_p_leak_ref: f64,
    facade_direction: FacadeDirection,
    altitude: f64,
    external_conditions: Arc<ExternalConditions>,
    p_a_alt: f64,
    // In Python there are extra properties:
    // n_leak - this is now N_LEAK as it is constant
    // c_leak_path - this is now calculated when needed with calculate_flow_coeff_for_leak
}

impl Leaks {
    /// Arguments:
    ///      extcond -- reference to ExternalConditions object
    ///      h_path -- mid height of the air flow path relative to ventlation zone floor level
    ///      delta_p_leak_ref -- Reference pressure difference (From pressure test e.g blower door = 50Pa)
    ///      qv_delta_p_leak_ref -- flow rate through
    ///      facade_direction -- The direction of the facade the leak is on.
    ///      A_roof -- Surface area of the roof of the ventilation zone (m2)
    ///      A_facades -- Surface area of facades (m2)
    ///      A_leak - Reference area of the envelope airtightness index qv_delta_p_leak_ref (depends on national context)
    ///      altitude -- altitude of dwelling above sea level (m)
    fn new(
        extcond: Arc<ExternalConditions>,
        h_path: f64,
        delta_p_leak_ref: f64,
        qv_delta_p_leak_ref: f64,
        facade_direction: FacadeDirection,
        a_roof: f64,
        a_facades: f64,
        a_leak: f64,
        altitude: f64,
    ) -> Self {
        Self {
            h_path,
            delta_p_leak_ref,
            a_roof,
            a_facades,
            a_leak,
            qv_delta_p_leak_ref,
            facade_direction,
            altitude,
            external_conditions: extcond,
            p_a_alt: adjust_air_density_for_altitude(altitude),
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

        (c_leak_path * f64::from(sign(delta_p_path)) * abs(delta_p_path)).powf(N_LEAK)
    }

    fn calculate_flow_from_internal_p(
        &self,
        u_site: f64,
        t_z: f64,
        p_z_ref: f64,
        f_cross: bool,
        shield_class: VentilationShieldClass,
        simtime: SimulationTimeIteration,
    ) -> (f64, f64) {
        let t_e = celsius_to_kelvin(self.external_conditions.air_temp(&simtime));

        // Wind pressure coefficient for the air flow path
        let c_p_path = get_c_p_path(f_cross, shield_class, self.h_path, self.facade_direction); // #TABLE from annex B

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
struct AirTerminalDevices {
    c_d_atd: f64,
    n_atd: f64,
    a_atd: f64,
    delta_p_atd_ref: f64,
    // NOTE - in Python we have c_atd_path as an instance variable but here we calculate it when needed instead
}

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

        -(f64::from(sign(qv_pdu))) * (abs(qv_pdu) / c_atd_path).powf(1. / self.n_atd)
    }
}

/// An object to represent Cowls
struct Cowls {
    c_p_cowl_roof: f64,
    height: f64,
    // NOTE - in Python we have delta_cowl_height as an instance variable but here we calculate it when needed from the height
}

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
struct CombustionAppliances {
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
    fn new(
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
struct MechanicalVentilation {
    f_ctrl: f64,
    f_sys: f64,
    e_v: f64,
    theta_z_t: f64,
    sup_air_flw_ctrl: SupplyAirFlowRateControlType,
    sup_air_temp_ctrl: SupplyAirTemperatureControlType,
    external_conditions: ExternalConditions,
    q_h_des: f64,
    q_c_des: f64,
    theta_ctrl_sys: Option<f64>,
    vent_type: VentType,
    total_volume: f64,
    ctrl_intermittent_mev: Option<Arc<Control>>,
    sfp: f64,
    simulation_timestep: f64,
    energy_supply_conn: EnergySupplyConnection,
    altitude: f64,
    design_outdoor_air_flow_rate_m3_h: f64,
    mvhr_eff: f64,
    qv_oda_req_design: f64,
    p_a_alt: f64,
}

impl MechanicalVentilation {
    /// Construct a Mechanical Ventilation object
    /// Arguments:
    /// external_conditions -- reference to ExternalConditions object
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
    fn new(
        extcond: ExternalConditions,
        _sup_air_flw_ctrl: SupplyAirFlowRateControlType,
        _sup_air_temp_ctrl: SupplyAirTemperatureControlType,
        q_h_des: f64,
        q_c_des: f64,
        vent_type: VentType,
        specific_fan_power: f64,
        design_outdoor_air_flow_rate: f64,
        // NOTE - in Python this is simtime
        simulation_timestep: f64,
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
            // Hard coded variables
            f_ctrl,
            f_sys,
            e_v,
            theta_z_t: 0., // (From Python) TODO get Thermal zone temperature - used for LOAD
            sup_air_flw_ctrl: SupplyAirFlowRateControlType::ODA, // (From Python) TODO currently hard coded until load comp implemented
            sup_air_temp_ctrl: SupplyAirTemperatureControlType::NoControl, // (From Python) TODO currently hard coded until load comp implemented
            // Arguments
            external_conditions: extcond,
            q_h_des,
            q_c_des,
            theta_ctrl_sys,
            vent_type,
            total_volume,
            ctrl_intermittent_mev,
            sfp: specific_fan_power,
            simulation_timestep,
            energy_supply_conn,
            altitude,
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
                    .setpnt(&simulation_time)
                    .expect("A setpoint was expected to be derivable for a control.");

                if f_op_v < 0. || f_op_v > 1. {
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
        simulation_time: &SimulationTimeIteration,
    ) -> (f64, f64, f64) {
        let t_e = celsius_to_kelvin(self.external_conditions.air_temp(simulation_time));

        // Required air flow at air terminal devices
        let (qv_sup_req, qv_eta_req) = self.calc_req_oda_flow_rates_at_atds();

        // Amount of air flow depends on controls
        let (f_op_v, qv_sup_dis_req, qv_eta_dis_req) = match self.sup_air_flw_ctrl {
            SupplyAirFlowRateControlType::ODA => {
                let f_op_v = self.f_op_v(simulation_time);
                let qv_sup_dis_req = f_op_v * qv_sup_req;
                let qv_eta_dis_req = f_op_v * qv_eta_req;

                (f_op_v, qv_sup_dis_req, qv_eta_dis_req)
            }
            SupplyAirFlowRateControlType::Load => {
                // NOTE - this is not currently implemented in the Python code
                unimplemented!("calc_mech_vent_air_flw_rates_req_to_supply_vent_zone is not implemented for SupplyAirFlowRateControlType::Load")
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

        (
            qm_sup_dis_req,
            qm_eta_dis_req,
            qm_in_effective_heat_recovery_saving,
        )
    }

    /// Calculate gains and energy use due to fans
    /// zone_volume -- volume of the zone (m3)
    /// total_volume -- volume of the dwelling (m3)
    /// vent_type -- one of "Intermittent MEV", "Centralised continuous MEV",
    /// "Decentralised continuous MEV", "MVHR" or "PIV".
    fn fans(
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
            * self.simulation_timestep
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

        supply_fan_energy_use_kwh / (f64::from(WATTS_PER_KILOWATT) * self.simulation_timestep)
    }

    pub fn vent_type(&self) -> VentType {
        self.vent_type
    }
}

/// A class to represent Infiltration and Ventilation object
struct InfiltrationVentilation {
    external_conditions: Arc<ExternalConditions>,
    f_cross: bool,
    shield_class: VentilationShieldClass,
    terrain_class: TerrainClass,
    c_rgh_site: f64,
    ventilation_zone_height: f64,
    windows: Vec<Window>,
    vents: Vec<Vent>,
    leaks: Vec<Leaks>,
    combustion_appliances: Vec<CombustionAppliances>,
    air_terminal_devices: Vec<AirTerminalDevices>,
    mech_vents: Vec<MechanicalVentilation>,
    detailed_output_heating_cooling: bool,
    p_a_alt: f64,
    total_volume: f64,
}

///Constructs a InfiltrationVentilation object
///
///         Arguments:
///             external_conditions -- reference to ExternalConditions object
///             simulation_time -- reference to SimulationTime object
///             f_cross -- cross-ventilation factor
///             shield_class -- indicates the exposure to wind of an air flow path on a facade
///                 (can can be open, normal and shielded)
///             ventilation_zone_height -- height of ventilation zone (m)
///             windows -- list of windows
///             vents -- list of vents
///             leaks -- required inputs for leaks
///             ATDs -- list of air terminal devices
///             mech_vents -- list of mech vents
///             altitude -- altitude of dwelling above sea level (m)
///             total_volume -- total zone volume
impl InfiltrationVentilation {
    fn new(
        external_conditions: Arc<ExternalConditions>,
        f_cross: bool,
        shield_class: VentilationShieldClass,
        terrain_class: TerrainClass,
        average_roof_pitch: f64,
        windows: Vec<Window>,
        vents: Vec<Vent>,
        leaks: VentilationLeaks,
        combustion_appliances: Vec<CombustionAppliances>,
        air_terminal_devices: Vec<AirTerminalDevices>,
        mech_vents: Vec<MechanicalVentilation>,
        detailed_output_heating_cooling: bool,
        altitude: f64,
        total_volume: f64,
    ) -> Self {
        Self {
            external_conditions: external_conditions.clone(),
            f_cross,
            shield_class,
            terrain_class,
            c_rgh_site: ter_class_to_roughness_coeff(terrain_class),
            ventilation_zone_height: leaks.ventilation_zone_height,
            windows,
            vents,
            leaks: Self::make_leak_objects(leaks, average_roof_pitch, external_conditions),
            combustion_appliances,
            air_terminal_devices,
            mech_vents,
            detailed_output_heating_cooling,
            p_a_alt: adjust_air_density_for_altitude(altitude),
            total_volume,
        }
    }

    /// Calculate supply temperature of the air flow element
    pub fn temp_supply(&self, simtime: SimulationTimeIteration) -> f64 {
        // NOTE: Technically, the MVHR system supplies air at a higher temperature
        // than the outside air, i.e.:
        //     temp_supply = self.__efficiency * temp_int_air \
        //                 + (1 - self.__efficiency) * self.__external_conditions.air_temp()
        // However, calculating this requires the internal air temperature, which
        // has not been calculated yet. Calculating this properly would require
        // the equation above to be added to the heat balance solver. Therefore,
        // it is simpler to adjust the heat transfer coefficient h_ve to account
        // for the heat recovery effect using an "equivalent" flow rate of
        // external air, which is done elsewhere
        self.external_conditions.air_temp(&simtime)
    }

    /// Calculate total volume air flow rate entering ventilation zone
    /// Equation 68 from BS EN 16798-7
    fn calculate_total_volume_air_flow_rate_in(qm_in: f64, external_air_density: f64) -> f64 {
        qm_in / external_air_density // from weather file?
    }

    /// Calculate total volume air flow rate leaving ventilation zone
    /// Equation 69 from BS EN 16798-7
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
        leaks: VentilationLeaks,
        average_roof_pitch: f64,
        external_conditions: Arc<ExternalConditions>,
    ) -> Vec<Leaks> {
        let h_path1_2 = 0.25 * leaks.ventilation_zone_height;
        let h_path3_4 = 0.75 * leaks.ventilation_zone_height;
        let h_path5 = leaks.ventilation_zone_height;
        let h_path_list = [h_path1_2, h_path1_2, h_path3_4, h_path3_4, h_path5];

        let roof_pitch = match average_roof_pitch {
            ..10.0 => FacadeDirection::Roof10,
            10.0..=30.0 => FacadeDirection::Roof10_30,
            30.0..60.0 => FacadeDirection::Roof30,
            _ => panic!("Average roof pitch was not expected to be greater than 60 degrees."),
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
                    external_conditions.clone(),
                    h_path_list[i],
                    leaks.test_pressure,
                    leaks.test_result,
                    facade_direction[i],
                    leaks
                        .area_roof
                        .expect("VentilationLeaks did not have an area_roof defined"),
                    leaks
                        .area_facades
                        .expect("VentilationLeaks did not have an area_facades defined"),
                    leaks.env_area,
                    leaks
                        .altitude
                        .expect("VentilationLeaks did not have an altitude defined"),
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
        h_z: f64,
        simtime: &SimulationTimeIteration,
    ) -> f64 {
        let func = |qv_pdu, p_z_ref, t_z, h_z| {
            self.implicit_formula_for_qv_pdu(qv_pdu, p_z_ref, t_z, h_z, simtime)
        };

        fsolve(func, qv_pdu, (p_z_ref, t_z, h_z)) // returns qv_pdu
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
        h_z: f64,
        simtime: &SimulationTimeIteration,
    ) -> f64 {
        let t_e = celsius_to_kelvin(self.external_conditions.air_temp(simtime));
        let external_air_density = air_density_at_temp(t_e, self.p_a_alt);
        let zone_air_density = air_density_at_temp(t_z, self.p_a_alt);

        // (From Python) TODO Standard isn't clear if delta_p_ATD can be totalled or not.
        let delta_p_atd_list: Vec<f64> = self
            .air_terminal_devices
            .iter()
            .map(|atd| atd.calculate_pressure_difference_atd(qv_pdu))
            .collect();

        let delta_p_atd: f64 = delta_p_atd_list.iter().sum();

        // Stack effect in passive and hybrid duct. As there is no air transfer
        // between levels of the ventilation zone Equation B.1 is used.
        let h_pdu_stack = h_z + 2.;

        // (From Python) TODO include delta_p_dpu and delta_p_cowl in the return.
        delta_p_atd - p_z_ref - h_pdu_stack * G * (external_air_density - zone_air_density)
    }

    /// The root scalar function will iterate until it finds a value of p_z_ref
    /// that satisfies the mass balance equation.
    /// The root scalar solver allows a range of intervals to be entered.
    /// The loop begins with a small interval to start with and if no solution is
    /// found or the boundary is too small for to cause a sign change then a wider
    /// interval is used until a solution is found.
    fn calculate_internal_reference_pressure(
        &self,
        intial_p_z_ref_guess: f64,
        temp_int_air: f64,
        r_w_arg: Option<f64>,
        simtime: &SimulationTimeIteration,
    ) -> f64 {
        let func = |qv_pdu, p_z_ref: f64, t_z: f64, h_z: f64| {
            self.implicit_formula_for_qv_pdu(qv_pdu, p_z_ref, t_z, h_z, simtime)
        };

        for interval_expansion in INTERVAL_EXPANSION_LIST {
            let result = root_scalar(
                func,
                [
                    intial_p_z_ref_guess - interval_expansion,
                    intial_p_z_ref_guess + interval_expansion,
                ],
                (temp_int_air, r_w_arg),
                "brentq",
            );

            match result {
                Ok(sol) => {
                    let p_z_ref = sol.root;
                    return p_z_ref;
                }
                Err(_) => {}
            }
        }

        panic!("Solver failed");
    }

    /// Used in calculate_internal_reference_pressure function for p_z_ref solve
    fn implicit_mass_balance_for_internal_reference_pressure(
        self,
        p_z_ref: f64,
        temp_int_air: f64,
        r_w_arg_min_max: f64,
        // flag = None,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        let (qm_in, qm_out, _) = self
            .implicit_mass_balance_for_internal_reference_pressure_components(
                p_z_ref,
                temp_int_air,
                r_w_arg_min_max,
                // flag,
                simtime,
            );
        qm_in + qm_out
    }

    /// Calculate incoming air flow, in m3/hr, at specified conditions
    fn incoming_air_flow(
        self,
        p_z_ref: f64,
        temp_int_air: f64,
        r_w_arg_min_max: f64,
        // reporting_flag = None,
        report_effective_flow_rate: bool,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        let (mut qm_in, _, qm_effective_flow_rate) = self
            .implicit_mass_balance_for_internal_reference_pressure_components(
                p_z_ref,
                temp_int_air,
                r_w_arg_min_max,
                //reporting_flag,
                simtime,
            );

        if report_effective_flow_rate {
            qm_in -= qm_effective_flow_rate
        }
        convert_mass_flow_rate_to_volume_flow_rate(
            qm_in,
            celsius_to_kelvin(self.external_conditions.air_temp(&simtime)),
            self.p_a_alt,
        )
    }

    /// Implicit mass balance for calculation of the internal reference pressure
    /// Equation 67 from BS EN 16798-7.
    ///
    /// Arguments:
    /// p_z_ref -- internal reference pressure (Pa)
    /// temp_int_air -- temperature of the intake air (K)
    /// reporting_flag -- flag used to give more detailed ventilation outputs (None = no additional reporting)
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
        temp_int_air: f64,
        r_w_arg_min_max: f64,
        // reporting_flag: bool,
        simtime: SimulationTimeIteration,
    ) -> (f64, f64, f64) {
        let wind_speed = self.external_conditions.wind_speed(&simtime);
        let u_site = wind_speed_at_zone_level(self.c_rgh_site, wind_speed, None, None, None);
        let t_e = celsius_to_kelvin(self.external_conditions.air_temp(&simtime));
        let t_z = celsius_to_kelvin(temp_int_air);
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
                u_site,
                t_z,
                p_z_ref,
                self.f_cross,
                self.shield_class,
                Some(r_w_arg_min_max),
                simtime,
            );
            qm_in_through_window_opening += qm_in;
            qm_out_through_window_opening += qm_out;
        }

        for vent in &self.vents {
            let (qm_in, qm_out) = vent.calculate_flow_from_internal_p(
                u_site,
                t_z,
                p_z_ref,
                self.f_cross,
                self.shield_class,
                simtime,
            );
            qm_in_through_vents += qm_in;
            qm_out_through_vents += qm_out;
        }

        for leak in &self.leaks {
            let (qm_in, qm_out) = leak.calculate_flow_from_internal_p(
                u_site,
                t_z,
                p_z_ref,
                self.f_cross,
                self.shield_class,
                simtime,
            );
            qm_in_through_leaks += qm_in;
            qm_out_through_leaks += qm_out;
        }

        for atd in &self.air_terminal_devices {
            let qv_pdu_initial = 0.; // (From Python) TODO get from prev timestep
            let h_z = self.ventilation_zone_height;
            let qv_pdu = self.calculate_qv_pdu(qv_pdu_initial, p_z_ref, t_z, h_z, &simtime);

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
            let p_h_fi = 0.; // (From Python) TODO to work out from previous zone temperature? - Combustion appliance heating fuel input power
            let f_op_comb = 1.; // (From Python) TODO work out what turns the appliance on or off. Schedule or Logic?
            let (qv_in, qv_out) =
                combustion_appliance.calculate_air_flow_req_for_comb_appliance(f_op_comb, p_h_fi);
            let (qm_in_comb, qm_out_comb) =
                convert_to_mass_air_flow_rate(qv_in, qv_out, t_e, t_z, self.p_a_alt);
            qm_in_through_comb += qm_in_comb;
            qm_out_through_comb += qm_out_comb;
        }

        for mech_vent in &self.mech_vents {
            let (qm_sup, qm_eta, qm_in_effective_heat_recovery_saving) =
                mech_vent.calc_mech_vent_air_flw_rates_req_to_supply_vent_zone(t_z, &simtime);
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
            let incoming_air_flow = convert_mass_flow_rate_to_volume_flow_rate(
                qm_in,
                celsius_to_kelvin(self.external_conditions.air_temp(&simtime)),
                self.p_a_alt,
            );
            let air_changes_per_hour = incoming_air_flow / self.total_volume;

            // TODO implement detailed reporting if required
            // Python has optional detailed reporting
            // which is currently not implemented here

            // if reporting_flag {
            //      ...
            // }
        }

        (qm_in, qm_out, qm_in_effective_heat_recovery_saving_total)
    }
}

// TODO this is from scipy.
// Find equivalent function in a Rust library or implement
fn root_scalar(
    func: impl FnOnce(f64, f64, f64, f64) -> f64,
    bracket: [f64; 2],
    args: (f64, Option<f64>),
    method: &str,
) -> Result<RootScalarResult, Error> {
    todo!()
}

struct RootScalarResult {
    pub root: f64,
}

// TODO this is from scipy
// Find equivalent function in a Rust library or implement
fn fsolve(func: impl FnOnce(f64, f64, f64, f64) -> f64, x0: f64, args: (f64, f64, f64)) -> f64 {
    // Stub implementation for the timebeing
    x0
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use rstest::rstest;
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
        assert_relative_eq!(result, -2.2966793114, max_relative = 1e-7); // Use spreadsheet to find answer.
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
        assert_eq!(ter_class_to_roughness_coeff(TerrainClass::OpenTerrain), 1.0);
        assert_eq!(ter_class_to_roughness_coeff(TerrainClass::Country), 0.9);
        assert_eq!(ter_class_to_roughness_coeff(TerrainClass::Urban), 0.8);
    }

    #[test]
    fn test_orientation_difference() {
        assert_eq!(orientation_difference(0., 90.), 90.);
        assert_eq!(orientation_difference(0., 450.), 90.); // 450 - 360 = 90
        assert_eq!(orientation_difference(180., 540.), 360.);
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
}
