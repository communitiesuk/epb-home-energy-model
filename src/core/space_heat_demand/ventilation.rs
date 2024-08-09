// This module provides objects to represent Infiltration and Ventilation.
// The calculations are based on Method 1 of BS EN 16798-7.

use crate::core::material_properties::AIR;
use crate::core::units::{celsius_to_kelvin, SECONDS_PER_HOUR};
use crate::external_conditions::ExternalConditions;
use crate::input::{
    CombustionAirSupplySituation, CombustionApplianceType, CombustionFuelType,
    FlueGasExhaustSituation, TerrainClass, VentilationShieldClass,
};
use crate::simulation_time::{self, SimulationTimeIteration};
use rand_distr::num_traits::abs;

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
    p_a_ref() * ((1. - ((0.00651 * h_alt) / 293.)) as f64).powf(4.255)
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

#[derive(Clone, Copy, PartialEq)]
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
    } else {
        if pitch > 60. {
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
}

// we split the python get_c_p_path method into two methods below:
fn get_c_p_path_from_pitch_and_orientation(
    f_cross: bool,
    shield_class: VentilationShieldClass,
    h_path: f64,
    wind_direction: f64,
    pitch: f64,
    orientation: f64,
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
        } else if 15. <= h_path && h_path < 50. {
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
    external_conditions: ExternalConditions,
    p_a_alt: f64, // In Python there are extra properties:
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
        extcond: ExternalConditions,
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

        let c_leak =
            (self.qv_delta_p_leak_ref * self.a_leak / (self.delta_p_leak_ref)).powf(N_LEAK);

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
    ///      T_e -- external air temperature (K)
    ///      T_z -- thermal zone air temperature (K)
    ///      C_p_path -- wind pressure coefficient at the height of the window part
    ///      p_z_ref -- internal reference pressure (Pa)
    fn calculate_ventilation_through_leaks_using_internal_p(
        &self,
        u_site: f64,
        t_e: f64,
        t_z: f64,
        c_p_path: f64,
        p_z_ref: f64,
    ) -> f64 {
        // For each couple of height and wind pressure coeficient associated with vents,
        // the air flow rate.
        let delta_p_path = calculate_pressure_difference_at_an_airflow_path(
            self.h_path,
            c_p_path,
            u_site,
            t_e,
            t_z,
            p_z_ref,
        );

        let c_leak_path = Self::calculate_flow_coeff_for_leak(&self);

        // Airflow through leaks based on Equation 62
        let qv_leak_path =
            (c_leak_path * f64::from(sign(delta_p_path)) * abs(delta_p_path)).powf(N_LEAK);
        qv_leak_path
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
            qv_in_through_leak += air_flow
        } else {
            qv_out_through_leak += air_flow
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

// TODO:
// WIP: a bunch of top level functions (called from other parts of ventilation.py)
// a Window class
// a WindowPart class
// a Vent class
// a Leaks class
// a AirTerminalDevices class
// a Cowls class
// a CombustionAppliances class
// a MechanicalVentilation class
// InfiltrationVentilation class
