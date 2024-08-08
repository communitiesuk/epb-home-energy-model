// This module provides objects to represent Infiltration and Ventilation.
// The calculations are based on Method 1 of BS EN 16798-7.

use crate::core::material_properties::AIR;
use crate::core::units::SECONDS_PER_HOUR;
use crate::input::{
    CombustionAirSupplySituation, CombustionApplianceType, CombustionFuelType,
    DiverterHeatSourceType, FlueGasExhaustSituation, StorageTankType, TerrainClass,
};
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
