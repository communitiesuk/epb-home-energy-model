use crate::core::space_heat_demand::building_element::{
    a_sol_for, area_for_building_element_input, h_ce_for, h_ci_for, h_pli_for, h_re_for, h_ri_for,
    i_sol_dir_dif_for, k_pli_for, number_of_building_element_nodes, number_of_inside_nodes,
    shading_factors_direct_diffuse_for, temp_ext_for, therm_rad_to_sky_for,
};
use crate::core::space_heat_demand::thermal_bridge::{
    heat_transfer_coefficient_for_thermal_bridge, ThermalBridging,
};
use crate::core::space_heat_demand::ventilation_element::{
    temp_supply_for_window_opening, MechanicalVentilationHeatRecovery, NaturalVentilation,
    VentilationElement, VentilationElementInfiltration, WholeHouseExtractVentilation,
    WindowOpeningForCooling,
};
use crate::core::units::{kelvin_to_celsius, SECONDS_PER_HOUR, WATTS_PER_KILOWATT};
use crate::external_conditions::ExternalConditions;
use crate::input::BuildingElement;
use crate::simulation_time::SimulationTimeIterator;
use nalgebra::{DMatrix, DVector};
use serde::Deserialize;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

// Convective fractions
// (default values from BS EN ISO 52016-1:2017, Table B.11)
const F_INT_C: f64 = 0.4; // Can be different for each source of internal gains
const F_SOL_C: f64 = 0.1;

// Areal thermal capacity of air and furniture
// (default value from BS EN ISO 52016-1:2017, Table B.17)
const K_M_INT: f64 = 10000.0; // J / (m2.K)

pub struct Zone<'a> {
    useful_area: f64,
    volume: f64,
    building_elements: Vec<NamedBuildingElement>,
    thermal_bridging: ThermalBridging,
    vent_elements: Vec<Box<dyn VentilationElement>>,
    vent_cool_extra: Option<WindowOpeningForCooling<'a>>,
    tb_heat_trans_coeff: f64,
    /// total area of all building elements associated with this zone, in m2
    area_el_total: f64,
    /// internal thermal capacity of the zone, in J / K
    c_int: f64,
    /// dictionary where key is building element (name) and
    ///                      values are 2-element tuples storing matrix row and
    ///                      column numbers (both same) where the first element
    ///                      of the tuple gives the position of the heat
    ///                      balance eqn (row) and node temperature (column)
    ///                      for the external surface and the second element
    ///                      gives the position for the internal surface.
    ///                      Positions in between will be for the heat balance
    ///                      and temperature of the inside nodes of the
    ///                      building element
    element_positions: Vec<(usize, usize)>,
    /// matrix row and column number (both same)
    ///                      corresponding to heat balance eqn for zone (row)
    ///                      and temperature of internal air (column)
    zone_idx: usize,
    /// number of unknown temperatures (each node in each
    ///                      building element + 1 for internal air) to be
    ///                      solved for
    no_of_temps: u32,
    /// list of temperatures (nodes and internal air) from
    ///                      previous timestep. Positions in list defined in
    ///                      element_positions and zone_idx
    temp_prev: Vec<f64>,
    external_conditions: ExternalConditions,
}

impl<'a> Zone<'a> {
    /// Construct a Zone object
    ///
    /// ## Arguments
    ///
    /// * `area` - useful floor area of the zone, in m2
    /// * `volume` - total volume of the zone, m3
    /// * `building_elements` - list of BuildingElement objects (walls, floors, windows etc.)
    /// * `thermal_bridging` - Either:
    ///                      - overall heat transfer coefficient for thermal
    ///                        bridges in the zone, in W / K
    ///                      - list of ThermalBridge objects for this zone
    /// * `vent_elements` - list of ventilation elements (infiltration, mech vent etc.)
    /// * `vent_cool_extra` - element providing additional ventilation in response to high
    ///                      internal temperature
    /// * `temp_ext_air_init` - external air temperature to use during initialisation, in Celsius
    /// * `temp_setpnt_init` - setpoint temperature to use during initialisation, in Celsius
    pub fn new(
        area: f64,
        volume: f64,
        building_elements: HashMap<String, BuildingElement>,
        thermal_bridging: ThermalBridging,
        vent_elements: Vec<Box<dyn VentilationElement>>,
        vent_cool_extra: Option<WindowOpeningForCooling>,
        temp_ext_air_init: f64,
        temp_setpnt_init: f64,
        external_conditions: ExternalConditions,
        simulation_time: SimulationTimeIterator,
    ) -> Zone {
        let tb_heat_trans_coeff = match thermal_bridging {
            ThermalBridging::Number(heat_coeff) => heat_coeff,
            ThermalBridging::Bridges(ref bridges) => bridges
                .values()
                .map(|bridge| heat_transfer_coefficient_for_thermal_bridge(bridge))
                .sum::<f64>(),
        };

        let area_el_total = building_elements
            .values()
            .map(|el| area_for_building_element_input(el))
            .sum::<f64>();
        let c_int = K_M_INT * area;

        // Calculate:
        // - size of required matrix/vectors (total number of nodes across all
        //                                      #   building elements + 1 for internal air)
        // - positions of heat balance eqns and temperatures in matrix for each node
        let mut element_positions = vec![];
        let mut n = 0;
        let mut named_building_elements = vec![];
        for (name, building_element) in building_elements.iter() {
            let start_idx = n;
            n = n + number_of_building_element_nodes(&building_element);
            let end_idx = n - 1;
            element_positions.push((start_idx as usize, end_idx as usize));
            named_building_elements.push(NamedBuildingElement {
                name: name.clone(),
                element: building_element.clone(),
            })
        }
        let zone_idx = n as usize;
        let no_of_temps = n + 1;

        let temp_prev = init_node_temps(
            temp_ext_air_init,
            temp_setpnt_init,
            no_of_temps,
            area_el_total,
            area,
            volume,
            &simulation_time,
            &vent_cool_extra,
            &external_conditions,
            &named_building_elements,
            &element_positions,
            &vent_elements,
            zone_idx,
            c_int,
            tb_heat_trans_coeff,
        );

        Zone {
            useful_area: area,
            volume,
            building_elements: named_building_elements,
            thermal_bridging,
            vent_elements,
            vent_cool_extra,
            tb_heat_trans_coeff,
            area_el_total,
            c_int,
            element_positions,
            zone_idx,
            no_of_temps,
            temp_prev,
            external_conditions,
        }
    }

    pub fn volume(&self) -> f64 {
        self.volume
    }
}

/// Initialise temperatures of heat balance nodes
//
/// ## Arguments
/// * `temp_ext_air_init` - external air temperature to use during initialisation, in Celsius
/// * `temp_setpnt_init` - setpoint temperature to use during initialisation, in Celsius
/// * `no_of_temps` - number of unknown temperatures (each node in each
///                             building element + 1 for internal air) to be
///                             solved for
pub fn init_node_temps(
    temp_ext_air_init: f64,
    temp_setpnt_init: f64,
    no_of_temps: u32,
    area_el_total: f64,
    area: f64,
    volume: f64,
    simulation_time: &SimulationTimeIterator,
    vent_cool_extra: &Option<WindowOpeningForCooling>,
    external_conditions: &ExternalConditions,
    building_elements: &Vec<NamedBuildingElement>,
    element_positions: &Vec<(usize, usize)>,
    vent_elements: &Vec<Box<dyn VentilationElement>>,
    passed_zone_idx: usize,
    c_int: f64,
    tb_heat_trans_coeff: f64,
) -> Vec<f64> {
    let no_of_temps = no_of_temps as usize;
    // Set starting point for all node temperatures (elements of
    //                                               # self.__temp_prev) as average of external air temp and setpoint. This
    // is somewhat arbitrary, but of all options for a uniform initial
    // temperature, this should lead to relatively fast stabilisation of
    // fabric temperatures, which are expected to be close to the external
    // air temperature towards the external surface nodes and close to the
    // setpoint temperature towards the internal surface nodes, and therefore
    // node temperatures on average should be close to the average of the
    // external air and setpoint temperatures.
    let temp_start = (temp_ext_air_init + temp_setpnt_init) / 2.0;
    let mut temp_prev = vec![temp_start; no_of_temps];

    // Iterate over space heating calculation and meet all space heating
    // demand until temperatures stabilise, under steady-state conditions
    // using specified constant setpoint and external air temperatures.
    loop {
        let (space_heat_demand, space_cool_demand, _) = space_heat_cool_demand(
            DELTA_T_H as f64,
            temp_ext_air_init,
            0.0,
            0.0,
            FRAC_CONVECTIVE,
            FRAC_CONVECTIVE,
            temp_setpnt_init,
            temp_setpnt_init,
            None,
            building_elements,
            element_positions,
            vent_elements,
            vent_cool_extra,
            simulation_time,
            external_conditions,
            &temp_prev,
            no_of_temps,
            area_el_total,
            area,
            volume,
            c_int,
            tb_heat_trans_coeff,
            passed_zone_idx,
        );

        // Note: space_cool_demand returned by function above is negative,
        // and only one of space_heat_demand and space_cool_demand will be
        // non-zero.
        let gains_heat_cool =
            (space_heat_demand + space_cool_demand) * WATTS_PER_KILOWATT as f64 / DELTA_T_H as f64;

        let temps_updated = calc_temperatures(
            DELTA_T as f64,
            &temp_prev,
            temp_ext_air_init,
            0.0,
            0.0,
            gains_heat_cool,
            FRAC_CONVECTIVE,
            None,
            None,
            no_of_temps,
            building_elements,
            element_positions,
            external_conditions,
            simulation_time,
            passed_zone_idx,
            area_el_total,
            volume,
            c_int,
            tb_heat_trans_coeff,
            vent_elements,
            vent_cool_extra,
        );

        if !isclose(&temps_updated, &temp_prev, Some(1e-08), None) {
            temp_prev = temps_updated;
        } else {
            break;
        }
    }

    temp_prev
}

const DELTA_T_H: u32 = 8760; // hours in a non leap year
const DELTA_T: u32 = DELTA_T_H * SECONDS_PER_HOUR;

// # Assume default convective fraction for heating/cooling suggested in
// # BS EN ISO 52016-1:2017 Table B.11
const FRAC_CONVECTIVE: f64 = 0.4;

/// Calculate heating and cooling demand in the zone for the current timestep
///
/// According to the procedure in BS EN ISO 52016-1:2017, section 6.5.5.2, steps 1 to 4.
///
/// # Arguments
/// * `delta_t_h` - calculation timestep, in hours
/// * `temp_ext_air` - temperature of the external air for the current timestep, in deg C
/// * `gains_internal` - internal gains for the current timestep, in W
/// * `gains_solar` - directly transmitted solar gains, in W
/// * `frac_convective_heat` - convective fraction for heating
/// * `frac_convective_cool` - convective fraction for cooling
/// * `temp_setpnt_heat` - temperature setpoint for heating, in deg C
/// * `temp_setpnt_cool` - temperature setpoint for cooling, in deg C
/// * `throughput_factor` - proportional increase in ventilation rate due to
///                         over-ventilation requirement
/// * `vent_cool_extra` - (optional) window cooling
/// * `a simulation_time`
pub fn space_heat_cool_demand(
    delta_t_h: f64,
    temp_ext_air: f64,
    gains_internal: f64,
    gains_solar: f64,
    frac_convective_heat: f64,
    frac_convective_cool: f64,
    temp_setpnt_heat: f64,
    temp_setpnt_cool: f64,
    throughput_factor: Option<f64>,
    building_elements: &Vec<NamedBuildingElement>,
    element_positions: &Vec<(usize, usize)>,
    vent_elements: &Vec<Box<dyn VentilationElement>>,
    vent_cool_extra: &Option<WindowOpeningForCooling>,
    simulation_time: &SimulationTimeIterator,
    external_conditions: &ExternalConditions,
    temp_prev: &Vec<f64>,
    no_of_temps: usize,
    area_el_total: f64,
    area: f64,
    volume: f64,
    c_int: f64,
    tb_heat_trans_coeff: f64,
    passed_zone_idx: usize,
) -> (f64, f64, f64) {
    let throughput_factor = throughput_factor.unwrap_or(1.0);
    assert!(
        temp_setpnt_cool >= temp_setpnt_heat,
        "ERROR: Cooling setpoint is below heating setpoint."
    );

    let mut temp_setpnt_cool_vent: Option<f64> = None;
    if let Some(window_cooling) = vent_cool_extra {
        // temp_setpnt_cool_vent =
        // Set cooling setpoint to Planck temperature to ensure no cooling demand
        let temp_setpnt_cool_vent = Some(kelvin_to_celsius(1.4e32));
        assert!(
            temp_setpnt_cool_vent.is_none() || (temp_setpnt_cool_vent.unwrap() >= temp_setpnt_heat),
            "ERROR: Setpoint for additional ventilation is below heating setpoint."
        );
    }

    // Calculate timestep in seconds
    let delta_t = delta_t_h * SECONDS_PER_HOUR as f64;

    // For calculation of demand, set heating/cooling gains to zero
    let gains_heat_cool = 0.0;

    // Calculate node and internal air temperatures with heating/cooling gains of zero
    let temp_vector_no_heat_cool = calc_temperatures(
        delta_t,
        temp_prev,
        temp_ext_air,
        gains_internal,
        gains_solar,
        gains_heat_cool,
        1.0, // Value does not matter as gains_heat_cool = 0.0
        None,
        Some(throughput_factor),
        no_of_temps,
        building_elements,
        element_positions,
        external_conditions,
        simulation_time,
        passed_zone_idx,
        area_el_total,
        volume,
        c_int,
        tb_heat_trans_coeff,
        vent_elements,
        vent_cool_extra,
    );

    // Calculate internal operative temperature at free-floating conditions
    // i.e. with no heating/cooling
    let mut temp_operative_free = temp_operative(
        &temp_vector_no_heat_cool,
        building_elements,
        element_positions,
        passed_zone_idx,
    );
    let temp_int_air_free = temp_vector_no_heat_cool[passed_zone_idx];

    // Check setpoint for additional ventilation. If above setpoint:
    // First calculate temps with max. additional ventilation, then check
    // setpoints again. If still above cooling setpoint, do not use additional
    // ventilation - just use cooling instead. Otherwise, cooling demand is zero
    // and need to use interpolation to work out additional ventilation required
    // (just like calc for heat_cool_load_unrestricted below)
    let mut h_ve_cool_extra = 0.0;
    if vent_cool_extra.is_some() && temp_operative_free > temp_setpnt_cool_vent.unwrap() {
        // Calculate node and internal air temperatures with maximum additional ventilation
        let h_ve_cool_max = vent_cool_extra.as_ref().unwrap().h_ve_max(
            volume,
            temp_operative_free,
            simulation_time.current_index(),
        );
        let temp_vector_vent_max = calc_temperatures(
            delta_t,
            temp_prev,
            temp_ext_air,
            gains_internal,
            gains_solar,
            gains_heat_cool,
            1.,
            Some(h_ve_cool_max),
            Some(throughput_factor),
            no_of_temps,
            building_elements,
            element_positions,
            external_conditions,
            simulation_time,
            passed_zone_idx,
            area_el_total,
            volume,
            c_int,
            tb_heat_trans_coeff,
            vent_elements,
            vent_cool_extra,
        );

        // Calculate internal operative temperature with maximum ventilation
        let temp_operative_vent_max = temp_operative(
            &temp_vector_vent_max,
            building_elements,
            element_positions,
            passed_zone_idx,
        );
        let temp_int_air_vent_max = temp_vector_vent_max[passed_zone_idx];

        let vent_cool_extra_temp_supply = temp_supply_for_window_opening(
            vent_cool_extra.as_ref().unwrap(),
            simulation_time.current_index(),
        );

        // If there is cooling potential from additional ventilation
        if temp_operative_vent_max < temp_operative_free
            && temp_int_air_free > vent_cool_extra_temp_supply
        {
            // Calculate ventilation required to reach cooling setpoint for ventilation
            let h_ve_cool_req = h_ve_cool_max
                * (temp_setpnt_cool_vent.unwrap() - temp_operative_free)
                / (temp_operative_vent_max - temp_operative_free)
                * ((temp_int_air_vent_max - vent_cool_extra_temp_supply)
                    / (temp_int_air_free - vent_cool_extra_temp_supply));

            // Calculate additional ventilation rate achieved
            let mut h_ve_cool_extra = match h_ve_cool_req < h_ve_cool_max {
                true => h_ve_cool_req,
                false => h_ve_cool_max,
            };

            // Calculate node and internal air temperatures with heating/cooling gains of zero
            let temp_vector_no_heat_cool_vent_extra = calc_temperatures(
                delta_t,
                temp_prev,
                temp_ext_air,
                gains_internal,
                gains_solar,
                gains_heat_cool,
                1.0,
                Some(h_ve_cool_extra),
                Some(throughput_factor),
                no_of_temps,
                building_elements,
                element_positions,
                external_conditions,
                simulation_time,
                passed_zone_idx,
                area_el_total,
                volume,
                c_int,
                tb_heat_trans_coeff,
                vent_elements,
                vent_cool_extra,
            );

            // Calculate internal operative temperature at free-floating conditions
            // i.e. with no heating/cooling
            let temp_operative_free_vent_extra = temp_operative(
                &temp_vector_no_heat_cool_vent_extra,
                building_elements,
                element_positions,
                passed_zone_idx,
            );

            // If temperature achieved by additional ventilation is above setpoint
            // for active cooling, assume cooling system will be used instead of
            // additional ventilation. Otherwise, use resultant operative temperature
            // in calculation of space heating/cooling demand.
            if temp_operative_free_vent_extra > temp_setpnt_cool {
                h_ve_cool_extra = 0.0;
            } else {
                temp_operative_free = temp_operative_free_vent_extra;
            }
        }
    }

    // Determine relevant setpoint (if neither, then return space heating/cooling demand of zero)
    // Determine maximum heating/cooling
    let mut temp_setpnt = 0.0;
    let mut heat_cool_load_upper = 0.0;
    let mut frac_convective = 0.0;
    if temp_operative_free > temp_setpnt_cool {
        // Cooling
        // TODO Implement eqn 26 "if max power available" case rather than just "otherwise" case?
        //      Could max. power be available at this point for all heating/cooling systems?
        temp_setpnt = temp_setpnt_cool;
        heat_cool_load_upper = -10. * area;
        frac_convective = frac_convective_cool;
    } else if temp_operative_free < temp_setpnt_heat {
        // Heating
        // TODO Implement eqn 26 "if max power available" case rather than just "otherwise" case?
        //      Could max. power be available at this point for all heating/cooling systems?
        temp_setpnt = temp_setpnt_heat;
        heat_cool_load_upper = 10. * area;
        frac_convective = frac_convective_heat;
    } else {
        return (0.0, 0.0, h_ve_cool_extra);
    }

    // Calculate node and internal air temperatures with maximum heating/cooling
    let temp_vector_upper_heat_cool = calc_temperatures(
        delta_t,
        temp_prev,
        temp_ext_air,
        gains_internal,
        gains_solar,
        heat_cool_load_upper,
        frac_convective,
        Some(h_ve_cool_extra),
        Some(throughput_factor),
        no_of_temps,
        building_elements,
        element_positions,
        external_conditions,
        simulation_time,
        passed_zone_idx,
        area_el_total,
        volume,
        c_int,
        tb_heat_trans_coeff,
        vent_elements,
        vent_cool_extra,
    );

    // Calculate internal operative temperature with maximum heating/cooling
    let temp_operative_upper = temp_operative(
        &temp_vector_upper_heat_cool,
        building_elements,
        element_positions,
        passed_zone_idx,
    );

    // Calculate heating (positive) or cooling (negative) required to reach setpoint
    let heat_cool_load_unrestricted = heat_cool_load_upper * (temp_setpnt - temp_operative_free)
        / (temp_operative_upper - temp_operative_free);

    // Convert from W to kWh
    let heat_cool_demand = heat_cool_load_unrestricted / WATTS_PER_KILOWATT as f64 * delta_t_h;

    let mut space_heat_demand = 0.0;
    let mut space_cool_demand = 0.0;
    if heat_cool_demand < 0.0 {
        space_cool_demand = heat_cool_demand;
    } else if heat_cool_demand > 0.0 {
        space_heat_demand = heat_cool_demand;
    }

    (space_heat_demand, space_cool_demand, h_ve_cool_extra)
}

/// Calculate temperatures according to procedure in BS EN ISO 52016-1:2017, section 6.5.6
//
/// ## Arguments
/// * `delta_t` - calculation timestep, in seconds
/// * `temp_prev` - temperature vector X (see below) from previous timestep
/// * `temp_ext_air` - temperature of external air, in deg C
/// * `gains_internal` - total internal heat gains, in W
/// * `gains_solar` - directly transmitted solar gains, in W
/// * `gains_heat_cool` - gains from heating (positive) or cooling (negative), in W
/// * `f_hc_c` - convective fraction for heating/cooling
/// * `vent_extra_h_ve` - additional ventilation heat transfer coeff in response
///                       to high internal temperature
/// * `throughput_factor` - proportional increase in ventilation rate due to
///                         over-ventilation requirement
/// * `no_of_temps` - number of unknown temperatures (each node in each
///                             building element + 1 for internal air) to be
///                             solved for
/// * `building_elements` - the building elements of the zone in question
/// * `element_positions` - dictionary where key is building element (name) and
///                      values are 2-element tuples storing matrix row and
///                      column numbers (both same) where the first element
///                      of the tuple gives the position of the heat
///                      balance eqn (row) and node temperature (column)
///                      for the external surface and the second element
///                      gives the position for the internal surface.
///                      Positions in between will be for the heat balance
///                      and temperature of the inside nodes of the
///                      building element
///
/// Temperatures are calculated by solving (for X) a matrix equation A.X = B, where:
/// A is a matrix of known coefficients
/// X is a vector of unknown temperatures
/// B is a vector of known quantities
///
/// Each row in vector X is a temperature variable - one for each node in each
/// building element plus the internal air temperature in the zone.
///
/// Each row of matrix A contains the coefficients from the heat balance equations
/// for each of the nodes in each building element, plus one row for the heat
/// balance equation of the zone.
///
/// Each column of matrix A contains the coefficients for a particular temperature
/// variable (in same order that they appear in vector X). Where the particular
/// temperature does not appear in the equation this coefficient will be zero.
///
/// Note that for this implementation, the columns and rows will be in corresponding
/// order, so the heat balance equation for node i will be in row i and the
/// coefficients in each row for the temperature at node i will be in column i.
///
/// Each row of vector B contains the other quantities (i.e. those that are not
/// coefficients of the temperature variables) from the heat balance equations
/// for each of the nodes in each building element, plus one row for the heat
/// balance equation of the zone, in the same order that the rows appear in matrix
/// A.
fn calc_temperatures(
    delta_t: f64,
    temp_prev: &Vec<f64>,
    temp_ext_air: f64,
    gains_internal: f64,
    gains_solar: f64,
    gains_heat_cool: f64,
    f_hc_c: f64,
    vent_extra_h_ve: Option<f64>,
    throughput_factor: Option<f64>,
    no_of_temps: usize,
    building_elements: &Vec<NamedBuildingElement>,
    element_positions: &Vec<(usize, usize)>,
    external_conditions: &ExternalConditions,
    simulation_time: &SimulationTimeIterator,
    passed_zone_idx: usize,
    area_el_total: f64,
    volume: f64,
    c_int: f64,
    tb_heat_trans_coeff: f64,
    vent_elements: &Vec<Box<dyn VentilationElement>>,
    vent_cool_extra: &Option<WindowOpeningForCooling>,
) -> Vec<f64> {
    let throughput_factor = throughput_factor.unwrap_or(1.0);
    let vent_extra_h_ve = vent_extra_h_ve.unwrap_or(0.0);

    // Init matrix with zeroes
    // Number of rows in matrix = number of columns
    // = total number of nodes + 1 for overall zone heat balance (and internal air temp)
    let mut matrix_a: DMatrix<f64> = DMatrix::zeros(no_of_temps, no_of_temps);

    // Init vector_b with zeroes (length = number of nodes + 1 for overall zone heat balance)
    let mut vector_b: DVector<f64> = DVector::zeros(no_of_temps);

    // One term in eqn 39 is sum from k = 1 to n of (A_elk / A_tot). Given
    // that A_tot is defined as the sum of A_elk from k = 1 to n, this term
    // will always evaluate to 1.
    // TODO Check this is correct. It seems a bit pointless if it is but we
    //      should probably retain it as an explicit term anyway to match
    //      the standard.
    let sum_area_frac = 1.0;

    // Node heat balances - loop through building elements and their nodes:
    // - Construct row of matrix_a for each node energy balance eqn
    // - Calculate RHS of node energy balance eqn and add to vector_b
    for (eli_idx, NamedBuildingElement { element: eli, .. }) in building_elements.iter().enumerate()
    {
        // External surface node (eqn 41)
        // Get position (row == column) in matrix previously calculated for the first (external) node
        let mut idx = element_positions[eli_idx].0;
        // Position of first (external) node within element is zero
        let mut i = 0usize;

        // load in k_pli, h_pli, h_ce and h_re for this element
        let (k_pli, h_pli, h_ce, h_re, h_ri, a_sol, therm_rad_to_sky) = (
            k_pli_for(&eli),
            h_pli_for(&eli),
            h_ce_for(&eli),
            h_re_for(&eli),
            h_ri_for(&eli),
            a_sol_for(&eli),
            therm_rad_to_sky_for(&eli),
        );

        // Coeff for temperature of this node
        matrix_a[(idx, idx)] = (k_pli[i] / delta_t) + h_ce + h_re + h_pli[i];
        // Coeff for temperature of next node
        matrix_a[(idx, idx + 1)] = -h_pli[i];
        // RHS of heat balance eqn for this node
        let (i_sol_dir, i_sol_dif) = i_sol_dir_dif_for(eli, external_conditions);
        let (f_sh_dir, f_sh_dif) = shading_factors_direct_diffuse_for(&eli, external_conditions);
        vector_b[idx] = (k_pli[i] / delta_t) * temp_prev[idx]
            + (h_ce + h_re)
                * temp_ext_for(eli, temp_ext_air, external_conditions, &simulation_time)
            + a_sol * (i_sol_dif * f_sh_dif + i_sol_dir * f_sh_dir)
            - therm_rad_to_sky;

        // Inside node(s), if any (eqn 40)
        for _ in 1..(number_of_inside_nodes(&eli) + 1) {
            i += 1;
            idx += 1;
            // Coeff for temperature of prev node
            matrix_a[(idx, idx - 1)] = -h_pli[i - 1];
            // Coeff for temperature of this node
            matrix_a[(idx, idx)] = (k_pli[i] / delta_t) + h_pli[i] + h_pli[i - 1];
            // Coeff for temperature of next node
            matrix_a[(idx, idx + 1)] = -h_pli[i];
            // RHS of heat balance eqn for this node
            vector_b[idx] = (k_pli[i] / delta_t) * temp_prev[idx];
        }

        // Internal surface node (eqn 39)
        idx += 1;
        assert_eq!(idx, element_positions[eli_idx].1);
        i += 1;
        assert_eq!(i as u32, number_of_building_element_nodes(&eli) - 1);
        // Get internal convective surface heat transfer coefficient, which
        // depends on direction of heat flow, which depends in temperature of
        // zone and internal surface
        let h_ci = h_ci_for(&eli, temp_prev[passed_zone_idx], temp_prev[idx]);
        // Coeff for temperature of prev node
        matrix_a[(idx, idx - 1)] = -h_pli[i - 1];
        // Coeff for temperature of this node
        matrix_a[(idx, idx)] = (k_pli[i] / delta_t) + h_ci + h_ri * sum_area_frac + h_pli[i - 1];
        // Add final sum term for LHS of eqn 39 in loop below.
        // These are coeffs for temperatures of internal surface nodes of
        // all building elements in the zone
        for (elk_idx, NamedBuildingElement { element: elk, .. }) in
            building_elements.iter().enumerate()
        {
            let col = element_positions[elk_idx].1;
            // The line below must be an adjustment to the existing value
            // to handle the case where col = idx (i.e. where we have
            // already partially set the value of the matrix element above
            // (before this loop) and do not want to overwrite it)
            matrix_a[(idx, col)] -= (area_for_building_element_input(elk) / area_el_total) * h_ri;
        }
        // Coeff for temperature of thermal zone
        matrix_a[(idx, passed_zone_idx)] = -h_ci;
        // RHS of heat balance eqn for this node
        vector_b[idx] = (k_pli[i] / delta_t)
            * temp_prev[idx]
            * ((1.0 - F_INT_C) * gains_internal
                + (1.0 - F_SOL_C) * gains_solar
                + (1.0 - f_hc_c) * gains_heat_cool)
            / area_el_total;
    }
    // Zone heat balance:
    // - Construct row of matrix A for zone heat balance eqn
    // - Calculate RHS of zone heat balance eqn and add to vector_b
    //
    // Coeff for temperature of thermal zone
    // TODO Throughput factor only applies to MVHR and WHEV, therefore only
    //      these systems accept throughput_factor as an argument to the h_ve
    //      function, hence the branch on the type in the loop below. This
    //      means that the MVHR and WHEV classes no longer have the same
    //      interface as other ventilation element classes, which could make
    //      future development more difficult. Ideally, we would find a
    //      cleaner way to implement this difference.
    let mut sum_vent_elements_h_ve = vent_extra_h_ve;
    for vei in vent_elements.iter() {
        sum_vent_elements_h_ve += vei.h_ve_heat_transfer_coefficient(
            volume,
            Some(throughput_factor),
            Some(simulation_time.current_index()),
        );
    }
    matrix_a[(passed_zone_idx, passed_zone_idx)] = (c_int / delta_t)
        + building_elements
            .iter()
            .enumerate()
            .map(|(eli_idx, nel)| {
                let NamedBuildingElement { element: eli, .. } = nel;
                h_ci_for(
                    eli,
                    temp_prev[passed_zone_idx],
                    temp_prev[element_positions[eli_idx].1],
                )
            })
            .sum::<f64>()
        + sum_vent_elements_h_ve
        + tb_heat_trans_coeff;
    // Add final sum term for LHS of eqn 38 in loop below.
    // These are coeffs for temperatures of internal surface nodes of
    // all building elements in the zone
    for (eli_idx, NamedBuildingElement { element: eli, .. }) in building_elements.iter().enumerate()
    {
        let col = element_positions[eli_idx].1; // Column for internal surface node temperature
        matrix_a[(passed_zone_idx, col)] = -area_for_building_element_input(eli)
            * h_ci_for(
                eli,
                temp_prev[passed_zone_idx],
                temp_prev[element_positions[eli_idx].1],
            );
    }
    // RHS of heat balance eqn for zone
    // TODO Throughput factor only applies to MVHR and WHEV, therefore only
    //      these systems accept throughput_factor as an argument to the h_ve
    //      function, hence the branch on the type in the loop below. This
    //      means that the MVHR and WHEV classes no longer have the same
    //      interface as other ventilation element classes, which could make
    //      future development more difficult. Ideally, we would find a
    //      cleaner way to implement this difference.
    let mut sum_vent_elements_h_ve_times_temp_supply = 0.0;
    if vent_extra_h_ve != 0.0 {
        sum_vent_elements_h_ve_times_temp_supply += vent_extra_h_ve
            * temp_supply_for_window_opening(
                vent_cool_extra.as_ref().expect(
                    "TODO: correct this - we are assuming there is always a window opening",
                ),
                simulation_time.current_index(),
            );
    }
    for vei in vent_elements.iter() {
        sum_vent_elements_h_ve_times_temp_supply += vei.h_ve_heat_transfer_coefficient(
            volume,
            Some(throughput_factor),
            Some(simulation_time.current_index()),
        ) * vei
            .temp_supply(simulation_time.current_index());
    }
    vector_b[passed_zone_idx] = (c_int * delta_t) * temp_prev[passed_zone_idx]
        + sum_vent_elements_h_ve_times_temp_supply
        + tb_heat_trans_coeff * temp_ext_air
        + F_INT_C * gains_internal
        + F_SOL_C * gains_solar
        + f_hc_c * gains_heat_cool;

    fast_solver(
        matrix_a,
        vector_b,
        no_of_temps,
        building_elements,
        element_positions,
        passed_zone_idx,
    )
}

/// Optimised heat balance solver
//
/// ## Arguments
/// * `coeffs` - full matrix of coefficients for the heat balance eqns
/// * `rhs` - full vector of values that are not temperatures or coefficients
///        (i.e. terms on right hand side of heat balance eqns)
/// * `no_of_temps` - number of unknown temperatures (each node in each
///                             building element + 1 for internal air) to be
///                             solved for
/// * `building_elements` - the building elements of the zone in question
/// * `element_positions` - dictionary where key is building element (name) and
///                      values are 2-element tuples storing matrix row and
///                      column numbers (both same) where the first element
///                      of the tuple gives the position of the heat
///                      balance eqn (row) and node temperature (column)
///                      for the external surface and the second element
///                      gives the position for the internal surface.
///                      Positions in between will be for the heat balance
///                      and temperature of the inside nodes of the
///                      building element
/// * `passed_zone_idx` - the index (to be) set on the zone
///
/// The heat balance equations from BS EN ISO 52016-1:2017 are expressed as a matrix equation and
/// solved simultaneously. While this provides a generic calculation procedure that works for an
/// arbitrary number of nodes (N), it also has a runtime proportional to N^3 which means that more
/// complex buildings can take a long time to simulate. However, many of the nodes are known not to
/// interact (e.g. the node in the middle of one wall has no heat transfer with the node in another
/// wall) and therefore we do not require the full flexibility of the matrix approach to solve for
/// every node temperature. The only part of the heat balance calculation where this flexibility is
/// needed is in the interaction between internal air and internal surfaces, so the calculation of the
/// other node temperatures can be removed from the matrix equation using algebraic substitution.
///
/// Consider generic heat balance eqns for a 4-node element:
///
/// A1a + B1b                               = Z1    # Heat balance at node a (external surface node)
/// A2a + B2b + C2c                         = Z2    # Heat balance at node b (inside node)
///       B3b + C3c + D3d                   = Z3    # Heat balance at node c (inside node)
///             C4c + D4d + J4j + K4k + Y4y = Z4    # Heat balance at node d (internal surface node)
/// where:
/// - a, b, c and d are the node temperatures in the building element to be solved for
/// - j and k are the node temperatures for the internal surfaces of other elements
/// - y is the internal air temperature
/// - A1 is the coefficient for temperature a in equation 1, A2 is the coefficient for temperature a in
///   equation 2, etc.
/// - Z1, Z2, etc. are the terms in equation 1, 2, etc. that are not the temperatures to be solved for
///   or their coefficients (i.e. Z1, Z2 etc. are the terms on the RHS of the equations)
///
/// The heat balance equation for node a (external surface node) can be rearranged to solve for a:
///
/// A1a + B1b = Z1
/// A1a = Z1 - B1b
/// a = (Z1 - B1b) / A1
///
/// Using the rearranged heat balance equation for node a, we can substitute a in the heat balance
/// equation for node b (next inside node) to eliminate a as a variable:
///
/// A2a + B2b + C2c = Z2
/// A2 * (Z1 - B1b) / A1 + B2b + C2c = Z2
///
/// Rearranging to consolidate the occurances of b gives a new heat balance equation for b:
///
/// A2 * Z1 / A1 + A2 * (- B1b) / A1 + B2b + C2c = Z2
/// b * (B2 - A2 * B1 / A1) + A2 * Z1 / A1 + C2c = Z2
/// (B2 - A2 * B1 / A1) * b + C2c = Z2 - A2 * Z1 / A1
///
/// This new heat balance equation can then be expressed in terms of modified versions of B2 and Z2:
///
/// B2'b + C2c = Z2'
/// where:
/// B2' = B2 - B1 * A2 / A1
/// Z2' = Z2 - Z1 * A2 / A1
///
/// The process can then be repeated, rearranging this new heat balance equation to solve for b and
/// then substituting into the heat balance equation for c:
///
/// b = (Z2' - C2c) / B2'
///
/// C3'c + D3d = Z3'
/// where:
/// C3' = C3 - C2 * B3 / B2'
/// Z3' = Z3 - Z2' * B3 / B2'
///
/// And repeated again to generate a new equation for node d:
///
/// c = (Z3' - D3d) / C3'
///
/// D4'd + J4j + K4k + Y4y = Z4'
/// where:
/// D4' = D4 - D3 * C4 / C3'
/// Z4' = Z4 - Z3' * C4 / C3'
///
/// At this point, we have reached the internal surface and we need the flexibility of the matrix
/// approach, but we have reduced the number of nodes to be solved for by 3. The process can be
/// repeated for each building element, and once the matrix solver has solved for the internal surface
/// temperatures, we can then go back through the other nodes, using the rearranged heat balance
/// equations that solve for (in this case) c, b and a.
///
/// In order to deal with different building elements having different numbers of nodes, we can express
/// the above relationships generically:
///
/// temperature[i] = (Z_adjusted[i] - coeff[i][i+1] * temperature[i+1]) / coeff_adjusted[i][i]
///
/// coeff_adjusted[i][i] = coeff[i][i] - coeff[i-1][i] * coeff[i][i-1] / coeff_adjusted[i-1][i-1]
/// Z_adjusted[i] = Z[i] - Z_adjusted[i-1] * coeff[i][i-1] / coeff_adjusted[i-1][i-1]
///
/// where i is the number of the node and its heat balance equation (counting from the external surface
/// to the internal surface), e.g. coeff[i-1][i] would be the coeffient for the temperature of the
/// current node in the heat balance equation for the previous node.
///
/// The optimised calculation procedure is therefore:
/// - Loop over nodes, from external surface to internal surface, and calculate adjusted coeffs and RHS
///   for each heat balance eqn
/// - Construct matrix eqn for inside and air nodes only
/// - Solve heat balance eqns for inside and air nodes using normal matrix solver
/// - Loop over nodes, from internal inside node (i.e. inside node nearest to the internal surface) to
///   external surface, and calculate temperatures in sequence
fn fast_solver(
    coeffs: DMatrix<f64>,
    rhs: DVector<f64>,
    no_of_temps: usize,
    building_elements: &Vec<NamedBuildingElement>,
    element_positions: &Vec<(usize, usize)>,
    passed_zone_idx: usize,
) -> Vec<f64> {
    // Init matrix with zeroes
    // Number of rows in matrix = number of columns
    // = total number of nodes + 1 for overall zone heat balance (and internal air temp)
    let mut coeffs_adj: DMatrix<f64> = DMatrix::zeros(no_of_temps, no_of_temps);

    // Init rhs_adj with zeroes (length = number of nodes + 1 for overall zone heat balance)
    let mut rhs_adj: DVector<f64> = DVector::zeros(no_of_temps);

    // Init matrix with zeroes
    // Number of rows in matrix = number of columns
    // = total number of internal surface nodes + 1 for internal air node
    let num_rows_cols_optimised = building_elements.len() + 1;
    let zone_idx = num_rows_cols_optimised - 1;
    let mut matrix_a: DMatrix<f64> =
        DMatrix::zeros(num_rows_cols_optimised, num_rows_cols_optimised);

    // Init vector_b with zeroes (length = number of internal surfaces + 1 for air node)
    let mut vector_b: DVector<f64> = DVector::zeros(num_rows_cols_optimised);

    for (el_idx, eli) in building_elements.iter().enumerate() {
        let (idx_ext_surface, idx_int_surface) = element_positions[el_idx];

        // No adjusted coeffs and RHS for external surface heat balance eqn
        coeffs_adj[(idx_ext_surface, idx_int_surface)] = coeffs[(idx_ext_surface, idx_int_surface)];
        rhs_adj[idx_ext_surface] = rhs[idx_ext_surface];

        // Loop over nodes, from inside node adjacent to external surface, to internal surface
        for idx in (idx_ext_surface + 1)..(idx_int_surface + 1) {
            // Calculate adjusted coeffs and RHS for each heat balance eqn
            coeffs_adj[(idx, idx)] = coeffs[(idx, idx)]
                - coeffs[(idx - 1, idx)] * coeffs[(idx, idx - 1)] / coeffs_adj[(idx - 1, idx - 1)];
            rhs_adj[idx] = rhs[idx]
                - rhs_adj[idx - 1] * coeffs[(idx, idx - 1)] / coeffs_adj[(idx - 1, idx - 1)];
        }

        // Construct matrix eqn for internal surface nodes only (and air node, after this loop)
        matrix_a[(el_idx, el_idx)] = coeffs_adj[(idx_int_surface, idx_int_surface)];
        vector_b[el_idx] = rhs_adj[idx_int_surface];

        for (el_idx_other, _) in building_elements.iter().enumerate() {
            if el_idx == el_idx_other {
                continue;
            }

            let idx_other_int_surface = element_positions[el_idx_other].1;
            matrix_a[(el_idx, el_idx_other)] = coeffs[(idx_int_surface, idx_other_int_surface)];
        }

        // Add coeff for air temperature to this element's internal surface heat balance eqn
        matrix_a[(el_idx, zone_idx)] = coeffs[(idx_int_surface, passed_zone_idx)];
        // Add coeff for this element's internal surface temp to the air node heat balance eqn
        matrix_a[(zone_idx, el_idx)] = coeffs[(passed_zone_idx, idx_int_surface)];
    }

    // Add rest of air node heat balance eqn to matrix
    // Coeffs for temperatures other than the air temp are added in the loop above
    matrix_a[(zone_idx, zone_idx)] = coeffs[(passed_zone_idx, passed_zone_idx)];
    vector_b[zone_idx] = rhs[passed_zone_idx];

    println!("{:?}", matrix_a);
    println!("{:?}", vector_b);

    // Solve heat balance eqns for inside and air nodes using normal matrix solver
    // Solve matrix eqn A.X = B to calculate vector_x (temperatures)
    // use LU solver with partial pivoting
    // NB. .full_piv_lu() would give full pivoting which is slower but theoretically more
    // numerically stable - may be able to speed this up with static matrices
    let vector_x = matrix_a.lu().solve(&vector_b).unwrap();

    // Init temperature with zeroes (length = number of nodes + 1 for overall zone heat balance)
    let mut temperatures = vec![0.0; no_of_temps as usize];
    temperatures[passed_zone_idx] = vector_x[zone_idx];

    // Populate node temperature results for each building element
    for (el_idx, _) in building_elements.iter().enumerate() {
        let (idx_ext_surface, idx_int_surface) = element_positions[el_idx];

        // Populate internal surface temperature result
        temperatures[idx_int_surface] = vector_x[el_idx];

        // Loop over nodes, from internal inside node (i.e. inside node nearest to the
        // internal surface) to external surface, and calculate temperatures in sequence
        for idx in (idx_ext_surface..idx_int_surface).rev() {
            temperatures[idx] = (rhs_adj[idx] - coeffs[(idx, idx + 1)] * temperatures[idx + 1])
                / coeffs_adj[(idx, idx)];
        }
    }

    temperatures
}

/// Calculate the operative temperature, in deg C
///
/// According to the procedure in BS EN ISO 52016-1:2017, section 6.5.5.3.
///
/// ## Arguments
/// * `temp_vector` - vector (list) of temperatures calculated from the heat balance equations
fn temp_operative(
    temp_vector: &Vec<f64>,
    building_elements: &Vec<NamedBuildingElement>,
    element_positions: &Vec<(usize, usize)>,
    passed_zone_id: usize,
) -> f64 {
    let temp_int_air = temp_vector[passed_zone_id];

    // Mean radiant temperature is weighted average of internal surface temperatures
    let temp_mean_radiant = building_elements
        .iter()
        .enumerate()
        .map(|(eli_idx, nel)| {
            let NamedBuildingElement { element: eli, .. } = nel;
            area_for_building_element_input(eli) * temp_vector[element_positions[eli_idx].1]
        })
        .sum::<f64>();

    (temp_int_air + temp_mean_radiant) / 2.0
}

/// Close-enough port of numpy's isclose function for use here
///
/// Returns a single boolean where two arrays are element-wise all equal within a tolerance.
///
/// See [docs for numpy.isclose](https://docs.rs/is_close/latest/is_close/)
pub fn isclose(a: &Vec<f64>, b: &Vec<f64>, rtol: Option<f64>, atol: Option<f64>) -> bool {
    let rtol = rtol.unwrap_or(1e-5);
    let atol = atol.unwrap_or(1e-8);
    all_close!(
        (*a).clone(),
        (*b).clone(),
        rel_tol = rtol,
        abs_tol = atol,
        method = is_close::ASYMMETRIC
    )
}

#[derive(Clone)]
pub struct NamedBuildingElement {
    pub name: String,
    pub element: BuildingElement,
}

// equality and hashing based on name for identity

impl PartialEq for NamedBuildingElement {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
    }
}

impl Eq for NamedBuildingElement {}

impl Hash for NamedBuildingElement {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.name.hash(state);
    }
}

struct HeatBalanceAirNode {
    solar_gains: f64,
    internal_gains: f64,
    heating_or_cooling_system_gains: f64,
    energy_to_change_internal_temperature: f64,
    heat_loss_through_thermal_bridges: f64,
    heat_loss_through_infiltration: f64,
    heat_loss_through_ventilation: f64,
    fabric_heat_loss: f64,
}

struct HeatBalanceInternalBoundary {
    fabric_int_air_convective: f64,
    fabric_int_sol: f64,
    fabric_int_int_gains: f64,
    fabric_int_heat_cool: f64,
}

struct HeatBalanceExternalBoundary {
    solar_gains: f64,
    internal_gains: f64,
    heating_or_cooling_system_gains: f64,
    thermal_bridges: f64,
    ventilation: f64,
    infiltration: f64,
    fabric_ext_air_convective: f64,
    fabric_ext_air_radiative: f64,
    fabric_ext_sol: f64,
    fabric_ext_sky: f64,
    opaque_fabric_ext: f64,
    transparent_fabric_ext: f64,
    ground_fabric_ext: f64,
    ztc_fabric_ext: f64,
    ztu_fabric_ext: f64,
}

struct HeatBalance {
    air_node: HeatBalanceAirNode,
    internal_boundary: HeatBalanceInternalBoundary,
    external_boundary: HeatBalanceExternalBoundary,
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::core::space_heat_demand::thermal_bridge::ThermalBridge;
    use crate::core::units::DAYS_IN_MONTH;
    use crate::external_conditions::DaylightSavingsConfig;
    use crate::input::{
        InfiltrationBuildType, InfiltrationShelterType, InfiltrationTestType,
        MassDistributionClass, ZoneInput,
    };
    use crate::simulation_time::{SimulationTime, HOURS_IN_DAY};
    use rstest::*;

    const BASE_AIR_TEMPS: [f64; 24] = [
        0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 7.5, 10.0, 12.5, 15.0, 19.5, 17.0,
        15.0, 12.0, 10.0, 7.0, 5.0, 3.0, 1.0,
    ];
    const BASE_WIND_SPEEDS: [f64; 24] = [
        4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.7, 5.4, 5.6, 5.3, 5.1, 4.8, 4.7, 4.6, 4.5, 4.2,
        4.9, 4.3, 4.4, 4.5, 4.3, 4.6,
    ];

    #[fixture]
    pub fn zone<'a>() -> Zone<'a> {
        let simulation_time = SimulationTime::new(0.0, 4.0, 1.0);

        let air_temp_day_jan = BASE_AIR_TEMPS;
        let air_temp_day_feb = BASE_AIR_TEMPS.map(|t| t + 1.0);
        let air_temp_day_mar = BASE_AIR_TEMPS.map(|t| t + 2.0);
        let air_temp_day_apr = BASE_AIR_TEMPS.map(|t| t + 3.0);
        let air_temp_day_may = BASE_AIR_TEMPS.map(|t| t + 4.0);
        let air_temp_day_jun = BASE_AIR_TEMPS.map(|t| t + 5.0);
        let air_temp_day_jul = BASE_AIR_TEMPS.map(|t| t + 6.0);
        let air_temp_day_aug = BASE_AIR_TEMPS.map(|t| t + 6.0);
        let air_temp_day_sep = BASE_AIR_TEMPS.map(|t| t + 5.0);
        let air_temp_day_oct = BASE_AIR_TEMPS.map(|t| t + 4.0);
        let air_temp_day_nov = BASE_AIR_TEMPS.map(|t| t + 3.0);
        let air_temp_day_dec = BASE_AIR_TEMPS.map(|t| t + 2.0);

        let mut air_temps = vec![];
        for (temps, days_in_month) in [
            (air_temp_day_jan, DAYS_IN_MONTH[0]),
            (air_temp_day_feb, DAYS_IN_MONTH[1]),
            (air_temp_day_mar, DAYS_IN_MONTH[2]),
            (air_temp_day_apr, DAYS_IN_MONTH[3]),
            (air_temp_day_may, DAYS_IN_MONTH[4]),
            (air_temp_day_jun, DAYS_IN_MONTH[5]),
            (air_temp_day_jul, DAYS_IN_MONTH[6]),
            (air_temp_day_aug, DAYS_IN_MONTH[7]),
            (air_temp_day_sep, DAYS_IN_MONTH[8]),
            (air_temp_day_oct, DAYS_IN_MONTH[9]),
            (air_temp_day_nov, DAYS_IN_MONTH[10]),
            (air_temp_day_dec, DAYS_IN_MONTH[11]),
        ] {
            air_temps.extend_from_slice(
                temps
                    .iter()
                    .cloned()
                    .cycle()
                    .take((days_in_month * HOURS_IN_DAY) as usize)
                    .collect::<Vec<f64>>()
                    .as_slice(),
            );
        }

        let wind_speed_day_jan = BASE_WIND_SPEEDS;
        let wind_speed_day_feb = BASE_WIND_SPEEDS.map(|t| t - 0.1);
        let wind_speed_day_mar = BASE_WIND_SPEEDS.map(|t| t - 0.2);
        let wind_speed_day_apr = BASE_WIND_SPEEDS.map(|t| t - 0.6);
        let wind_speed_day_may = BASE_WIND_SPEEDS.map(|t| t - 0.8);
        let wind_speed_day_jun = BASE_WIND_SPEEDS.map(|t| t - 1.1);
        let wind_speed_day_jul = BASE_WIND_SPEEDS.map(|t| t - 1.2);
        let wind_speed_day_aug = BASE_WIND_SPEEDS.map(|t| t - 1.2);
        let wind_speed_day_sep = BASE_WIND_SPEEDS.map(|t| t - 1.1);
        let wind_speed_day_oct = BASE_WIND_SPEEDS.map(|t| t - 0.7);
        let wind_speed_day_nov = BASE_WIND_SPEEDS.map(|t| t - 0.5);
        let wind_speed_day_dec = BASE_WIND_SPEEDS.map(|t| t - 0.3);

        let mut wind_speeds = vec![];
        for (temps, days_in_month) in [
            (wind_speed_day_jan, DAYS_IN_MONTH[0]),
            (wind_speed_day_feb, DAYS_IN_MONTH[1]),
            (wind_speed_day_mar, DAYS_IN_MONTH[2]),
            (wind_speed_day_apr, DAYS_IN_MONTH[3]),
            (wind_speed_day_may, DAYS_IN_MONTH[4]),
            (wind_speed_day_jun, DAYS_IN_MONTH[5]),
            (wind_speed_day_jul, DAYS_IN_MONTH[6]),
            (wind_speed_day_aug, DAYS_IN_MONTH[7]),
            (wind_speed_day_sep, DAYS_IN_MONTH[8]),
            (wind_speed_day_oct, DAYS_IN_MONTH[9]),
            (wind_speed_day_nov, DAYS_IN_MONTH[10]),
            (wind_speed_day_dec, DAYS_IN_MONTH[11]),
        ] {
            wind_speeds.extend_from_slice(
                temps
                    .iter()
                    .cloned()
                    .cycle()
                    .take((days_in_month * HOURS_IN_DAY) as usize)
                    .collect::<Vec<f64>>()
                    .as_slice(),
            );
        }

        let external_conditions = ExternalConditions::new(
            simulation_time.iter(),
            air_temps,
            wind_speeds,
            vec![0.0; 4],
            vec![0.0; 4],
            vec![0.2; 4],
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
        );

        // Create objects for the different building elements in the zone
        let be_opaque_i = BuildingElement::Opaque {
            area: 20.0,
            pitch: 180.,
            a_sol: 0.6,
            r_c: Some(0.25),
            k_m: 19000.,
            mass_distribution_class: MassDistributionClass::I,
            orientation: 0.,
            base_height: 0.,
            height: 2.,
            width: 10.,
            u_value: None,
            h_ci: None,
            h_ri: None,
            h_ce: None,
            h_re: None,
        };
        let be_opaque_d = BuildingElement::Opaque {
            area: 26.0,
            pitch: 180.,
            a_sol: 0.55,
            r_c: Some(0.33),
            k_m: 16000.,
            mass_distribution_class: MassDistributionClass::D,
            orientation: 0.,
            base_height: 0.,
            height: 2.,
            width: 10.,
            u_value: None,
            h_ci: None,
            h_ri: None,
            h_ce: None,
            h_re: None,
        };
        let be_ztc = BuildingElement::AdjacentZTC {
            area: 22.5,
            pitch: 135.,
            r_c: Some(0.5),
            k_m: 18000.,
            mass_distribution_class: MassDistributionClass::E,
            u_value: None,
        };
        let be_ground = BuildingElement::Ground {
            area: 25.0,
            pitch: 90.,
            u_value: 1.33,
            r_f: 0.2,
            k_m: 17000.,
            mass_distribution_class: MassDistributionClass::IE,
            h_pi: 2.2,
            h_pe: 2.7,
            perimeter: 20.0,
            psi_wall_floor_junc: 0.7,
        };
        let be_transparent = BuildingElement::Transparent {
            pitch: 90.,
            r_c: Some(0.4),
            orientation: 180.,
            g_value: 0.75,
            frame_area_fraction: 0.25,
            base_height: 1.0,
            height: 1.25,
            width: 4.0,
            shading: vec![],
            u_value: None,
            area: None,
        };
        let be_ztu = BuildingElement::AdjacentZTUSimple {
            area: 30.0,
            pitch: 130.,
            r_c: Some(0.5),
            r_u: 0.6,
            k_m: 18000.,
            mass_distribution_class: MassDistributionClass::E,
            u_value: None,
        };

        // Put building element objects in a list that can be iterated over
        let be_objs = HashMap::from([
            ("be_opaque_i".to_string(), be_opaque_i),
            ("be_opaque_d".to_string(), be_opaque_d),
            ("be_ztc".to_string(), be_ztc),
            ("be_ground".to_string(), be_ground),
            ("be_transparent".to_string(), be_transparent),
            ("be_ztu".to_string(), be_ztu),
        ]);

        // Create objects for thermal bridges
        let tb_linear_1 = ThermalBridge::Linear {
            linear_thermal_transmittance: 0.28,
            length: 5.,
        };
        let tb_linear_2 = ThermalBridge::Linear {
            linear_thermal_transmittance: 0.25,
            length: 6.,
        };
        let tb_point = ThermalBridge::Point {
            heat_transfer_coefficient: 1.4,
        };

        // Put thermal bridge objects in a list that can be iterated over
        let thermal_bridging = ThermalBridging::Bridges(HashMap::from([
            ("tb_linear_1".to_string(), tb_linear_1),
            ("tb_linear_2".to_string(), tb_linear_2),
            ("tb_point".to_string(), tb_point),
        ]));

        // Create ventilation objects
        let ve = VentilationElementInfiltration::new(
            1,
            InfiltrationShelterType::Sheltered,
            InfiltrationBuildType::House,
            4.5,
            InfiltrationTestType::FiftyPascals,
            40.,
            75.,
            2,
            2,
            2,
            1,
            0,
            0,
            0,
            3,
            6,
            0,
            external_conditions.clone(),
        );

        // Put thermal ventilation objects in a list that can be iterated over
        let ve_objs: Vec<Box<dyn VentilationElement>> = vec![Box::new(ve)];

        let temp_ext_air_init = 17.;
        let temp_setpnt_init = 21.;

        let external_conditions = external_conditions.clone();

        Zone::new(
            50.,
            125.,
            be_objs,
            thermal_bridging,
            ve_objs,
            None,
            temp_ext_air_init,
            temp_setpnt_init,
            external_conditions,
            simulation_time.iter(),
        )
    }

    // commented out for now as zone instantiation does not resolve at the moment
    // #[rstest]
    // pub fn should_have_correct_volume(zone: Zone) {
    //     assert_eq!(zone.volume(), 125.);
    // }

    #[rstest]
    pub fn should_replicate_numpy_isclose() {
        // test cases for python doctests
        assert!(!isclose(
            &vec![1e10, 1e-7],
            &vec![1.00001e10, 1e-8],
            None,
            None,
        ));
        assert!(isclose(
            &vec![1e10, 1e-8],
            &vec![1.00001e10, 1e-9],
            None,
            None
        ));
        assert!(!isclose(
            &vec![1e10, 1e-8],
            &vec![1.0001e10, 1e-9],
            None,
            None
        ));
        assert!(!isclose(&vec![1e-8, 1e-7], &vec![0.0, 0.0], None, None),);
        assert!(!isclose(
            &vec![1e-100, 1e-7],
            &vec![0.0, 0.0],
            None,
            Some(0.0)
        ));
        assert!(isclose(&vec![1e-10, 1e-10], &vec![1e-20, 0.0], None, None),);
        assert!(!isclose(
            &vec![1e-10, 1e-10],
            &vec![1e-20, 0.999999e-10],
            None,
            Some(0.0),
        ));
    }
}
