// This module provides objects to represent the thermal zones in the building,
// and to calculate the temperatures in the zone and associated building elements.

use crate::compare_floats::is_close;
use crate::core::controls::time_control::ControlBehaviour;
use crate::core::material_properties::AIR;
use crate::core::space_heat_demand::building_element::BuildingElement;
use crate::core::space_heat_demand::thermal_bridge::{
    heat_transfer_coefficient_for_thermal_bridge, ThermalBridging,
};
use crate::core::space_heat_demand::ventilation::InfiltrationVentilation;
use crate::core::units::{kelvin_to_celsius, SECONDS_PER_HOUR, WATTS_PER_KILOWATT};
use crate::corpus::TempInternalAirFn;
use crate::input::ZoneTemperatureControlBasis;
use crate::simulation_time::{SimulationTimeIteration, SimulationTimeIterator};
use anyhow::bail;
use field_types::FieldName;
use indexmap::IndexMap;
use nalgebra::{DMatrix, DVector};
use parking_lot::RwLock;
use serde_enum_str::Serialize_enum_str;
use smartstring::alias::String;
use std::hash::{Hash, Hasher};
use std::mem;
use std::sync::Arc;
use thiserror::Error;

// Convective fractions
// (default values from BS EN ISO 52016-1:2017, Table B.11)
const F_INT_C: f64 = 0.4;
// Can be different for each source of internal gains
const F_SOL_C: f64 = 0.1;

// Areal thermal capacity of air and furniture
// (default value from BS EN ISO 52016-1:2017, Table B.17)
const K_M_INT: f64 = 10000.0; // J / (m2.K)

/// Calculate ventilation heat transfer co-efficient from air changes per hour
pub(crate) fn calc_vent_heat_transfer_coeff(volume: f64, air_changes_per_hour: f64) -> f64 {
    let q_ve = air_changes_per_hour * volume / SECONDS_PER_HOUR as f64;
    AIR.density_kg_per_m3() * AIR.specific_heat_capacity() * q_ve
}

#[derive(Debug)]
pub struct Zone {
    useful_area: f64,
    volume: f64,
    building_elements: Vec<NamedBuildingElement>,
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
    no_of_temps: usize,
    temp_prev: Arc<RwLock<Vec<f64>>>,
    print_heat_balance: bool,
    // Python has a use_fast_solver field that we don't need because we always use the equivalent fast solver in Rust
    _ventilation: Arc<InfiltrationVentilation>,
    control: Option<Arc<dyn ControlBehaviour>>,
    temp_setpnt_basis: ZoneTemperatureControlBasis,
    /// list of temperatures (nodes and internal air) from
    ///                      previous timestep. Positions in list defined in
    ///                      element_positions and zone_idx
    temp_setpnt_init: f64,
}

impl Zone {
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
    /// * `ventilation` - reference to ventilation object
    /// * `temp_ext_air_init` - external air temperature to use during initialisation, in Celsius
    /// * `temp_setpnt_init` - setpoint temperature to use during initialisation, in Celsius
    /// * `control` - reference to a control (generally setpoint time control)
    /// * `simulation_time`
    pub(crate) fn new(
        area: f64,
        volume: f64,
        building_elements: IndexMap<String, Arc<BuildingElement>>,
        thermal_bridging: ThermalBridging,
        ventilation: Arc<InfiltrationVentilation>,
        temp_ext_air_init: f64,
        temp_setpnt_init: f64,
        temp_setpnt_basis: ZoneTemperatureControlBasis,
        control: Option<Arc<dyn ControlBehaviour>>,
        print_heat_balance: bool,
        simulation_time: &SimulationTimeIterator,
    ) -> anyhow::Result<Self> {
        let tb_heat_trans_coeff = match thermal_bridging {
            ThermalBridging::Number(heat_coeff) => heat_coeff,
            ThermalBridging::Bridges(ref bridges) => bridges
                .values()
                .map(heat_transfer_coefficient_for_thermal_bridge)
                .sum::<f64>(),
        };

        let area_el_total = building_elements.values().map(|el| el.area()).sum::<f64>();
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
            n += building_element.number_of_nodes();
            let end_idx = n - 1;
            element_positions.push((start_idx, end_idx));
            named_building_elements.push(NamedBuildingElement {
                name: name.clone(),
                element: building_element.clone(),
            })
        }
        let zone_idx = n;
        let no_of_temps = n + 1;

        let zone = Zone {
            useful_area: area,
            volume,
            building_elements: named_building_elements,
            tb_heat_trans_coeff,
            area_el_total,
            c_int,
            element_positions,
            zone_idx,
            no_of_temps,
            temp_prev: Arc::new(RwLock::new(Vec::new())),
            print_heat_balance,
            _ventilation: ventilation,
            control,
            temp_setpnt_basis,
            temp_setpnt_init,
        };

        zone.init_node_temps(
            temp_ext_air_init,
            temp_setpnt_init,
            *simulation_time
                .clone()
                .peekable()
                .peek()
                .expect("Non-zero simulation period expected."),
        )?;

        Ok(zone)
    }

    /// Return temp_setpnt_init
    pub(crate) fn setpnt_init(&self) -> f64 {
        self.temp_setpnt_init
    }

    fn init_node_temps(
        &self,
        temp_ext_air_init: f64,
        temp_setpnt_init: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<()> {
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
        *self.temp_prev.write() = vec![temp_start; self.no_of_temps];

        // Iterate over space heating calculation and meet all space heating
        // demand until temperatures stabilise, under steady-state conditions
        // using specified constant setpoint and external air temperatures.
        loop {
            let (space_heat_demand, space_cool_demand, _, _) = self.space_heat_cool_demand(
                DELTA_T_H as f64,
                temp_ext_air_init,
                0.0,
                0.0,
                FRAC_CONVECTIVE,
                FRAC_CONVECTIVE,
                temp_setpnt_init,
                temp_setpnt_init,
                temp_ext_air_init,
                None,
                None,
                AirChangesPerHourArgument::from_ach_cooling(0.0),
                simtime,
            )?;

            // Note: space_cool_demand returned by function above is negative,
            // and only one of space_heat_demand and space_cool_demand will be
            // non-zero.
            let gains_heat_cool = (space_heat_demand + space_cool_demand)
                * WATTS_PER_KILOWATT as f64
                / DELTA_T_H as f64;

            let (temps_updated, _) = self.calc_temperatures(
                DELTA_T as f64,
                self.temp_prev.read().as_ref(),
                temp_ext_air_init,
                0.0,
                0.0,
                gains_heat_cool,
                FRAC_CONVECTIVE,
                0.0,
                temp_ext_air_init,
                simtime,
                self.print_heat_balance,
            )?;

            if !isclose(
                &temps_updated,
                self.temp_prev.read().as_ref(),
                Some(1e-08),
                None,
            ) {
                *self.temp_prev.write() = temps_updated;
            } else {
                break;
            }
        }

        Ok(())
    }

    pub(crate) fn area(&self) -> f64 {
        self.useful_area
    }

    pub(crate) fn volume(&self) -> f64 {
        self.volume
    }

    /// sum solar gains for all elements in the zone
    /// only transparent elements will have solar gains > 0
    pub(crate) fn gains_solar(&self, simulation_time: SimulationTimeIteration) -> f64 {
        self.building_elements
            .iter()
            .map(|el| el.element.solar_gains(simulation_time).unwrap())
            .sum::<f64>()
    }

    /// Calculate ventilation heat transfer coefficient from air changes per hour
    fn calc_vent_heat_transfer_coeff(&self, air_changes_per_hour: f64) -> f64 {
        calc_vent_heat_transfer_coeff(self.volume, air_changes_per_hour)
    }

    /// Calculate temperatures according to procedure in BS EN ISO 52016-1:2017, section 6.5.6
    ///
    ///     ## Arguments:
    ///     * `delta_t`         -- calculation timestep, in seconds
    ///     * `temp_prev`       -- temperature vector X (see below) from previous timestep
    ///     * `temp_ext_air`    -- temperature of external air, in deg C
    ///     * `gains_internal`  -- total internal heat gains, in W
    ///     * `gains_solar`     -- directly transmitted solar gains, in W
    ///     * `gains_heat_cool` -- gains from heating (positive) or cooling (negative), in W
    ///     * `f_hc_c`          -- convective fraction for heating/cooling
    ///     * `ach`             -- air changes per hour
    ///     * `print_heat_balance` -- flag to record whether to return the heat balance outputs
    ///     * `avg_supply_temp` -- average supply temperature
    ///
    ///     Temperatures are calculated by solving (for X) a matrix equation A.X = B, where:
    ///     A is a matrix of known coefficients
    ///     X is a vector of unknown temperatures
    ///     B is a vector of known quantities
    ///
    ///     Each row in vector X is a temperature variable - one for each node in each
    ///     building element plus the internal air temperature in the zone.
    ///
    ///     Each row of matrix A contains the coefficients from the heat balance equations
    ///     for each of the nodes in each building element, plus one row for the heat
    ///     balance equation of the zone.
    ///
    ///     Each column of matrix A contains the coefficients for a particular temperature
    ///     variable (in same order that they appear in vector X). Where the particular
    ///     temperature does not appear in the equation this coefficient will be zero.
    ///
    ///     Note that for this implementation, the columns and rows will be in corresponding
    ///     order, so the heat balance equation for node i will be in row i and the
    ///     coefficients in each row for the temperature at node i will be in column i.
    ///
    ///     Each row of vector B contains the other quantities (i.e. those that are not
    ///     coefficients of the temperature variables) from the heat balance equations
    ///     for each of the nodes in each building element, plus one row for the heat
    ///     balance equation of the zone, in the same order that the rows appear in matrix
    ///     A.
    fn calc_temperatures(
        &self,
        delta_t: f64,
        temp_prev: &[f64],
        temp_ext_air: f64,
        gains_internal: f64,
        gains_solar: f64,
        gains_heat_cool: f64,
        f_hc_c: f64,
        ach: f64,
        avg_supply_temp: f64,
        simtime: SimulationTimeIteration,
        print_heat_balance: bool,
    ) -> anyhow::Result<(Vec<f64>, Option<HeatBalance>)> {
        let h_ve = self.calc_vent_heat_transfer_coeff(ach);

        // Init matrix with zeroes
        // Number of rows in matrix = number of columns
        // = total number of nodes + 1 for overall zone heat balance (and internal air temp)
        let mut matrix_a: DMatrix<f64> = DMatrix::zeros(self.no_of_temps, self.no_of_temps);

        // Init vector_b with zeroes (length = number of nodes + 1 for overall zone heat balance)
        let mut vector_b: DVector<f64> = DVector::zeros(self.no_of_temps);

        // One term in eqn 39 is sum from k = 1 to n of (A_elk / A_tot). Given
        // that A_tot is defined as the sum of A_elk from k = 1 to n, this term
        // will always evaluate to 1.
        // TODO (from Python) Check this is correct. It seems a bit pointless if it is but we
        //      should probably retain it as an explicit term anyway to match
        //      the standard.
        let sum_area_frac = 1.0;

        // Node heat balances - loop through building elements and their nodes:
        // - Construct row of matrix_a for each node energy balance eqn
        // - Calculate RHS of node energy balance eqn and add to vector_b
        for (eli_idx, NamedBuildingElement { element: eli, .. }) in
            self.building_elements.iter().enumerate()
        {
            // External surface node (eqn 41)
            // Get position (row == column) in matrix previously calculated for the first (external) node
            let mut idx = self.element_positions[eli_idx].0;
            // Position of first (external) node within element is zero
            let mut i = 0usize;

            // load in k_pli, h_ce and h_re for this element
            let (k_pli, h_ce, h_re, h_ri, solar_absorption_coeff, therm_rad_to_sky) = (
                eli.k_pli(),
                eli.h_ce(),
                eli.h_re(),
                eli.h_ri(),
                eli.solar_absorption_coeff(),
                eli.therm_rad_to_sky(),
            );

            // Coeff for temperature of this node
            matrix_a[(idx, idx)] =
                (k_pli[i] / delta_t) + h_ce + h_re + eli.h_pli_by_index_unchecked(i, simtime)?;

            // Coeff for temperature of next node
            matrix_a[(idx, idx + 1)] = -eli.h_pli_by_index_unchecked(i, simtime)?;
            // RHS of heat balance eqn for this node
            let (i_sol_dir, i_sol_dif) = eli.i_sol_dir_dif(simtime);
            let (f_sh_dir, f_sh_dif) = eli.shading_factors_direct_diffuse(simtime).unwrap();
            vector_b[idx] = (k_pli[i] / delta_t) * temp_prev[idx]
                + (h_ce + h_re) * eli.temp_ext(simtime)
                + solar_absorption_coeff * (i_sol_dif * f_sh_dif + i_sol_dir * f_sh_dir)
                - therm_rad_to_sky;

            // Inside node(s), if any (eqn 40)
            for _ in 1..(eli.number_of_inside_nodes() + 1) {
                i += 1;
                idx += 1;
                // Coeff for temperature of prev node
                matrix_a[(idx, idx - 1)] = -eli.h_pli_by_index_unchecked(i - 1, simtime)?;
                // Coeff for temperature of this node
                matrix_a[(idx, idx)] = (k_pli[i] / delta_t)
                    + eli.h_pli_by_index_unchecked(i, simtime)?
                    + eli.h_pli_by_index_unchecked(i - 1, simtime)?;
                // Coeff for temperature of next node
                matrix_a[(idx, idx + 1)] = -eli.h_pli_by_index_unchecked(i, simtime)?;
                // RHS of heat balance eqn for this node
                vector_b[idx] = (k_pli[i] / delta_t) * temp_prev[idx];
            }

            // Internal surface node (eqn 39)
            idx += 1;
            debug_assert_eq!(idx, self.element_positions[eli_idx].1);
            i += 1;
            debug_assert_eq!(i, eli.number_of_nodes() - 1);
            // Get internal convective surface heat transfer coefficient, which
            // depends on direction of heat flow, which depends on temperature of
            // zone and internal surface
            let h_ci = eli.h_ci(temp_prev[self.zone_idx], temp_prev[idx]);
            // Coeff for temperature of prev node
            matrix_a[(idx, idx - 1)] = -eli.h_pli_by_index_unchecked(i - 1, simtime)?;
            // Coeff for temperature of this node
            matrix_a[(idx, idx)] = (k_pli[i] / delta_t)
                + h_ci
                + h_ri * sum_area_frac
                + eli.h_pli_by_index_unchecked(i - 1, simtime)?;
            // Add final sum term for LHS of eqn 39 in loop below.
            // These are coeffs for temperatures of internal surface nodes of
            // all building elements in the zone
            for (elk_idx, NamedBuildingElement { element: elk, .. }) in
                self.building_elements.iter().enumerate()
            {
                let col = self.element_positions[elk_idx].1;
                // The line below must be an adjustment to the existing value
                // to handle the case where col = idx (i.e. where we have
                // already partially set the value of the matrix element above
                // (before this loop) and do not want to overwrite it)
                matrix_a[(idx, col)] -= (elk.area() / self.area_el_total) * h_ri;
            }
            // Coeff for temperature of thermal zone
            matrix_a[(idx, self.zone_idx)] = -h_ci;
            // RHS of heat balance eqn for this node
            vector_b[idx] = (k_pli[i] / delta_t) * temp_prev[idx]
                + ((1.0 - F_INT_C) * gains_internal
                    + (1.0 - F_SOL_C) * gains_solar
                    + (1.0 - f_hc_c) * gains_heat_cool)
                    / self.area_el_total;
        }
        // Zone heat balance:
        // - Construct row of matrix A for zone heat balance eqn
        // - Calculate RHS of zone heat balance eqn and add to vector_b
        //
        // Coeff for temperature of thermal zone
        // TODO (from Python) Throughput factor only applies to MVHR and WHEV, therefore only
        //      these systems accept throughput_factor as an argument to the h_ve
        //      function, hence the branch on the type in the loop below. This
        //      means that the MVHR and WHEV classes no longer have the same
        //      interface as other ventilation element classes, which could make
        //      future development more difficult. Ideally, we would find a
        //      cleaner way to implement this difference.
        matrix_a[(self.zone_idx, self.zone_idx)] = (self.c_int / delta_t)
            + self
                .building_elements
                .iter()
                .enumerate()
                .map(|(eli_idx, nel)| {
                    let NamedBuildingElement { element: eli, .. } = nel;
                    eli.area()
                        * eli.h_ci(
                            temp_prev[self.zone_idx],
                            temp_prev[self.element_positions[eli_idx].1],
                        )
                })
                .sum::<f64>()
            + h_ve
            + self.tb_heat_trans_coeff;

        // Add final sum term for LHS of eqn 38 in loop below.
        // These are coeffs for temperatures of internal surface nodes of
        // all building elements in the zone
        for (eli_idx, NamedBuildingElement { element: eli, .. }) in
            self.building_elements.iter().enumerate()
        {
            let col = self.element_positions[eli_idx].1; // Column for internal surface node temperature
            matrix_a[(self.zone_idx, col)] = -eli.area()
                * eli.h_ci(
                    temp_prev[self.zone_idx],
                    temp_prev[self.element_positions[eli_idx].1],
                );
        }
        // RHS of heat balance eqn for zone
        vector_b[self.zone_idx] = (self.c_int / delta_t) * temp_prev[self.zone_idx]
            + h_ve * avg_supply_temp
            + self.tb_heat_trans_coeff * temp_ext_air
            + F_INT_C * gains_internal
            + F_SOL_C * gains_solar
            + f_hc_c * gains_heat_cool;

        let vector_x = self.fast_solver(matrix_a, vector_b)?;

        let heat_balance = print_heat_balance.then(|| {
            // Collect outputs, in W, for heat balance at air node
            let temp_internal = vector_x[self.zone_idx];
            let hb_gains_solar = F_SOL_C * gains_solar;
            let hb_gains_internal = F_INT_C * gains_internal;
            let hb_gains_heat_cool = f_hc_c * gains_heat_cool;
            let hb_energy_to_change_temp =
                -(self.c_int / delta_t) * (temp_internal - temp_prev[self.zone_idx]);
            let hb_loss_thermal_bridges = self.tb_heat_trans_coeff * (temp_internal - temp_ext_air);
            let hb_loss_infiltration_ventilation = h_ve * (temp_internal - avg_supply_temp);
            let hb_loss_fabric = (hb_gains_solar
                + hb_gains_internal
                + hb_gains_heat_cool
                + hb_energy_to_change_temp)
                - (hb_loss_thermal_bridges + hb_loss_infiltration_ventilation);
            let air_node = HeatBalanceAirNode {
                solar_gains: hb_gains_solar,
                internal_gains: hb_gains_internal,
                heating_or_cooling_system_gains: hb_gains_heat_cool,
                energy_to_change_internal_temperature: hb_energy_to_change_temp,
                thermal_bridges: -hb_loss_thermal_bridges,
                infiltration_ventilation: -hb_loss_infiltration_ventilation,
                fabric_heat_loss: -hb_loss_fabric,
            };

            // Collect outputs, in W, for heat balance at internal fabric boundary
            let fabric_int_sol = (1.0 - F_SOL_C) * gains_solar;
            let fabric_int_int_gains = (1.0 - F_INT_C) * gains_internal;
            let fabric_int_heat_cool = (1.0 - f_hc_c) * gains_heat_cool;

            let fabric_int_air_convective = self
                .building_elements
                .iter()
                .enumerate()
                .map(|(eli_idx, eli)| {
                    let idx = self.element_positions[eli_idx].1;
                    let temp_int_surface = vector_x[idx];
                    let air_node_temp = vector_x[self.zone_idx];

                    eli.element.area()
                        * ((eli.element.h_ci(
                            temp_prev[self.zone_idx],
                            temp_prev[self.element_positions[eli_idx].1],
                        )) * (air_node_temp - temp_int_surface))
                })
                .sum::<f64>();

            let internal_boundary = HeatBalanceInternalBoundary {
                fabric_int_air_convective,
                fabric_int_sol,
                fabric_int_int_gains,
                fabric_int_heat_cool,
            };

            // Collect outputs, in W, for heat balance at external boundary
            let mut hb_fabric_ext_air_convective = 0.0;
            let mut hb_fabric_ext_air_radiative = 0.0;
            let mut hb_fabric_ext_sol = 0.0;
            let mut hb_fabric_ext_sky = 0.0;
            let mut hb_fabric_ext_opaque = 0.0;
            let mut hb_fabric_ext_transparent = 0.0;
            let mut hb_fabric_ext_ground = 0.0;
            let mut hb_fabric_ext_ztc = 0.0;
            let mut hb_fabric_ext_ztu = 0.0;

            for (eli_idx, NamedBuildingElement { element: eli, .. }) in
                self.building_elements.iter().enumerate()
            {
                // Get position in vector for the first (external) node of the building element
                let idx = self.element_positions[eli_idx].0;
                let temp_ext_surface = vector_x[idx];
                let (i_sol_dir, i_sol_dif) = eli.i_sol_dir_dif(simtime);
                let (f_sh_dir, f_sh_dif) = eli
                    .shading_factors_direct_diffuse(simtime)
                    .expect("Expected shading factors direct diffuse to be calculable.");
                hb_fabric_ext_air_convective +=
                    eli.area() * eli.h_ce() * (eli.temp_ext(simtime) - temp_ext_surface);
                hb_fabric_ext_air_radiative +=
                    eli.area() * eli.h_re() * (eli.temp_ext(simtime) - temp_ext_surface);
                hb_fabric_ext_sol += eli.area()
                    * eli.solar_absorption_coeff()
                    * (i_sol_dif * f_sh_dif + i_sol_dir * f_sh_dir);
                hb_fabric_ext_sky += eli.area() * (-eli.therm_rad_to_sky());
                // fabric heat loss per building element type
                let hb_fabric_ext = eli.area()
                    * ((eli.h_ce()) * (eli.temp_ext(simtime) - temp_ext_surface))
                    + eli.area() * (eli.h_re()) * (eli.temp_ext(simtime) - temp_ext_surface)
                    + eli.area()
                        * eli.solar_absorption_coeff()
                        * (i_sol_dif * f_sh_dif + i_sol_dir * f_sh_dir)
                    + eli.area() * (-eli.therm_rad_to_sky());
                match eli.as_ref() {
                    BuildingElement::Opaque(_) => {
                        hb_fabric_ext_opaque += hb_fabric_ext;
                    }
                    BuildingElement::Transparent(_) => {
                        hb_fabric_ext_transparent += hb_fabric_ext;
                    }
                    BuildingElement::Ground(_) => {
                        hb_fabric_ext_ground += hb_fabric_ext;
                    }
                    BuildingElement::AdjacentConditionedSpace(_) => {
                        hb_fabric_ext_ztc += hb_fabric_ext;
                    }
                    BuildingElement::AdjacentUnconditionedSpaceSimple(_) => {
                        hb_fabric_ext_ztu += hb_fabric_ext;
                    }
                };
            }
            let external_boundary = HeatBalanceExternalBoundary {
                solar_gains: gains_solar,
                internal_gains: gains_internal,
                heating_or_cooling_system_gains: gains_heat_cool,
                thermal_bridges: -hb_loss_thermal_bridges,
                infiltration_ventilation: -hb_loss_infiltration_ventilation,
                fabric_ext_air_convective: hb_fabric_ext_air_convective,
                fabric_ext_air_radiative: hb_fabric_ext_air_radiative,
                fabric_ext_sol: hb_fabric_ext_sol,
                fabric_ext_sky: hb_fabric_ext_sky,
                opaque_fabric_ext: hb_fabric_ext_opaque,
                transparent_fabric_ext: hb_fabric_ext_transparent,
                ground_fabric_ext: hb_fabric_ext_ground,
                ztc_fabric_ext: hb_fabric_ext_ztc,
                ztu_fabric_ext: hb_fabric_ext_ztu,
            };

            HeatBalance {
                air_node,
                internal_boundary,
                external_boundary,
            }
        });

        Ok((vector_x, heat_balance)) // pass empty heat balance map for now
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
        &self,
        coeffs: DMatrix<f64>,
        rhs: DVector<f64>,
    ) -> Result<Vec<f64>, SolverError> {
        // Init matrix with zeroes
        // Number of rows in matrix = number of columns
        // = total number of nodes + 1 for overall zone heat balance (and internal air temp)
        let mut coeffs_adj: DMatrix<f64> = DMatrix::zeros(self.no_of_temps, self.no_of_temps);

        // Init rhs_adj with zeroes (length = number of nodes + 1 for overall zone heat balance)
        let mut rhs_adj: DVector<f64> = DVector::zeros(self.no_of_temps);

        // Init matrix with zeroes
        // Number of rows in matrix = number of columns
        // = total number of internal surface nodes + 1 for internal air node
        let num_rows_cols_optimised = self.building_elements.len() + 1;
        let zone_idx = num_rows_cols_optimised - 1;
        let mut matrix_a: DMatrix<f64> =
            DMatrix::zeros(num_rows_cols_optimised, num_rows_cols_optimised);

        // Init vector_b with zeroes (length = number of internal surfaces + 1 for air node)
        let mut vector_b: DVector<f64> = DVector::zeros(num_rows_cols_optimised);

        for (el_idx, _eli) in self.building_elements.iter().enumerate() {
            let (idx_ext_surface, idx_int_surface) = self.element_positions[el_idx];

            // No adjusted coeffs and RHS for external surface heat balance eqn
            coeffs_adj[(idx_ext_surface, idx_ext_surface)] =
                coeffs[(idx_ext_surface, idx_ext_surface)];
            rhs_adj[idx_ext_surface] = rhs[idx_ext_surface];

            // Loop over nodes, from inside node adjacent to external surface, to internal surface
            for idx in (idx_ext_surface + 1)..(idx_int_surface + 1) {
                // Calculate adjusted coeffs and RHS for each heat balance eqn
                coeffs_adj[(idx, idx)] = coeffs[(idx, idx)]
                    - coeffs[(idx - 1, idx)] * coeffs[(idx, idx - 1)]
                        / coeffs_adj[(idx - 1, idx - 1)];
                rhs_adj[idx] = rhs[idx]
                    - rhs_adj[idx - 1] * coeffs[(idx, idx - 1)] / coeffs_adj[(idx - 1, idx - 1)];
            }

            // Construct matrix eqn for internal surface nodes only (and air node, after this loop)
            matrix_a[(el_idx, el_idx)] = coeffs_adj[(idx_int_surface, idx_int_surface)];
            vector_b[el_idx] = rhs_adj[idx_int_surface];

            for (el_idx_other, _) in self.building_elements.iter().enumerate() {
                if el_idx == el_idx_other {
                    continue;
                }

                let idx_other_int_surface = self.element_positions[el_idx_other].1;
                matrix_a[(el_idx, el_idx_other)] = coeffs[(idx_int_surface, idx_other_int_surface)];
            }

            // Add coeff for air temperature to this element's internal surface heat balance eqn
            matrix_a[(el_idx, zone_idx)] = coeffs[(idx_int_surface, self.zone_idx)];
            // Add coeff for this element's internal surface temp to the air node heat balance eqn
            matrix_a[(zone_idx, el_idx)] = coeffs[(self.zone_idx, idx_int_surface)];
        }

        // Add rest of air node heat balance eqn to matrix
        // Coeffs for temperatures other than the air temp are added in the loop above
        matrix_a[(zone_idx, zone_idx)] = coeffs[(self.zone_idx, self.zone_idx)];
        vector_b[zone_idx] = rhs[self.zone_idx];

        // Solve heat balance eqns for inside and air nodes using normal matrix solver
        // Solve matrix eqn A.X = B to calculate vector_x (temperatures)
        // use LU solver with partial pivoting
        // NB. .full_piv_lu() would give full pivoting which is slower but theoretically more
        // numerically stable - may be able to speed this up with static matrices
        let vector_x = matrix_a.lu().solve(&vector_b).unwrap();

        if vector_x.iter().any(|&val| val.is_nan()) {
            return Err(SolverError::SolverWillNotResolve(
                "A NaN value was detected in the output of a linear algebra calculation.",
            ));
        }

        // Init temperature with zeroes (length = number of nodes + 1 for overall zone heat balance)
        let mut temperatures = vec![0.0; self.no_of_temps];
        temperatures[self.zone_idx] = vector_x[zone_idx];

        // Populate node temperature results for each building element
        for (el_idx, _) in self.building_elements.iter().enumerate() {
            let (idx_ext_surface, idx_int_surface) = self.element_positions[el_idx];

            // Populate internal surface temperature result
            temperatures[idx_int_surface] = vector_x[el_idx];

            // Loop over nodes, from internal inside node (i.e. inside node nearest to the
            // internal surface) to external surface, and calculate temperatures in sequence
            for idx in (idx_ext_surface..idx_int_surface).rev() {
                temperatures[idx] = (rhs_adj[idx] - coeffs[(idx, idx + 1)] * temperatures[idx + 1])
                    / coeffs_adj[(idx, idx)];
            }
        }

        Ok(temperatures)
    }

    // Calculate the operative temperature, in deg C
    ///
    /// According to the procedure in BS EN ISO 52016-1:2017, section 6.5.5.3.
    ///
    /// ## Arguments
    /// * `temp_vector` - vector (list) of temperatures calculated from the heat balance equations
    fn temp_operative_by_temp_vector(&self, temp_vector: &[f64]) -> f64 {
        let temp_int_air = temp_vector[self.zone_idx];

        // Mean radiant temperature is weighted average of internal surface temperatures
        let temp_mean_radiant = self
            .building_elements
            .iter()
            .enumerate()
            .map(|(eli_idx, nel)| {
                let NamedBuildingElement { element: eli, .. } = nel;
                eli.area() * temp_vector[self.element_positions[eli_idx].1]
            })
            .sum::<f64>()
            / self.area_el_total;
        (temp_int_air + temp_mean_radiant) / 2.0
    }

    pub(crate) fn temp_operative(&self) -> f64 {
        self.temp_operative_by_temp_vector(&self.temp_prev.read())
    }

    /// Return internal air temperature, in deg C
    pub(crate) fn temp_internal_air(&self) -> f64 {
        self.temp_prev.read()[self.zone_idx]
    }

    fn ach_req_to_reach_temperature(
        temp_target: f64,
        ach_min: f64,
        ach_max: f64,
        temp_ach_min: f64,
        temp_ach_max: f64,
        temp_int_air_min: f64,
        temp_int_air_max: f64,
        temp_supply: f64,
    ) -> f64 {
        if temp_ach_max >= temp_ach_min
            || temp_int_air_max >= temp_int_air_min
            || temp_int_air_min <= temp_supply
        {
            return ach_min;
        }

        let frac_interp = (temp_target - temp_ach_min) / (temp_ach_max - temp_ach_min);
        let temp_int_air_req =
            temp_int_air_min + frac_interp * (temp_int_air_max - temp_int_air_min);

        if temp_int_air_req <= temp_supply {
            return ach_max;
        }

        let ach_req = (ach_max
            * frac_interp
            * ((temp_int_air_max - temp_supply) / (temp_int_air_min - temp_supply))
            + ach_min * (1.0 - frac_interp))
            / ((temp_int_air_req - temp_supply) / (temp_int_air_min - temp_supply));

        ach_req.max(ach_min).min(ach_max)
    }

    fn calc_cooling_potential_from_ventilation(
        &self,
        delta_t: f64,
        temp_ext_air: f64,
        gains_internal: f64,
        gains_solar: f64,
        gains_heat_cool: f64,
        frac_conv_gains_heat_cool: f64,
        temp_setpnt_heat: f64,
        temp_setpnt_cool: f64,
        temp_setpnt_cool_vent: f64,
        temp_free: f64,
        temp_int_air_free: f64,
        ach_args: AirChangesPerHourArgument,
        avg_supply_temp: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64, Option<f64>)> {
        let mut temp_free = temp_free;

        // If ach required for cooling has not been provided, check if the
        // maximum air changes when window shut and open are different
        let (ach_to_trigger_heating, ach_cooling) = match ach_args {
            AirChangesPerHourArgument::TargetAndWindowsOpen {
                ach_target,
                ach_windows_open,
            } if ach_windows_open != ach_target => {
                // Calculate node and internal air temperatures with maximum additional ventilation
                let (temp_vector_vent_max, _) = self.calc_temperatures(
                    delta_t,
                    self.temp_prev.read().as_ref(),
                    temp_ext_air,
                    gains_internal,
                    gains_solar,
                    gains_heat_cool,
                    frac_conv_gains_heat_cool,
                    ach_windows_open,
                    avg_supply_temp,
                    simtime,
                    self.print_heat_balance,
                )?;
                let temp_int_air_vent_max = temp_vector_vent_max[self.zone_idx];

                let temp_vent_max = match self.temp_setpnt_basis {
                    ZoneTemperatureControlBasis::Operative => {
                        self.temp_operative_by_temp_vector(&temp_vector_vent_max)
                    }
                    ZoneTemperatureControlBasis::Air => temp_int_air_vent_max,
                };

                let ach_to_trigger_heating = Self::ach_req_to_reach_temperature(
                    temp_setpnt_heat,
                    ach_target,
                    ach_windows_open,
                    temp_free,
                    temp_vent_max,
                    temp_int_air_free,
                    temp_int_air_vent_max,
                    avg_supply_temp,
                );

                // If there is cooling potential from additional ventilation, and
                // free-floating temperature exceeds setpoint for additional ventilation
                let ach_cooling = if temp_vent_max < temp_free
                    && temp_free > temp_setpnt_cool_vent
                    && temp_int_air_free > avg_supply_temp
                {
                    // Calculate ventilation required to reach cooling setpoint for ventilation
                    let mut ach_cooling = Self::ach_req_to_reach_temperature(
                        temp_setpnt_cool_vent,
                        ach_target,
                        ach_windows_open,
                        temp_free,
                        temp_vent_max,
                        temp_int_air_free,
                        temp_int_air_vent_max,
                        avg_supply_temp,
                    );

                    // Calculate node and internal air temperatures with heating/cooling gains of zero
                    let (temp_vector_no_heat_cool_vent_extra, _) = self.calc_temperatures(
                        delta_t,
                        self.temp_prev.read().as_ref(),
                        temp_ext_air,
                        gains_internal,
                        gains_solar,
                        gains_heat_cool,
                        frac_conv_gains_heat_cool,
                        ach_cooling,
                        avg_supply_temp,
                        simtime,
                        self.print_heat_balance,
                    )?;

                    // Calculate internal operative temperature at free-floating conditions
                    // i.e. with no heating/cooling
                    let temp_free_vent_extra = match self.temp_setpnt_basis {
                        ZoneTemperatureControlBasis::Operative => {
                            self.temp_operative_by_temp_vector(&temp_vector_no_heat_cool_vent_extra)
                        }
                        ZoneTemperatureControlBasis::Air => {
                            temp_vector_no_heat_cool_vent_extra[self.zone_idx]
                        }
                    };

                    // If temperature achieved by additional ventilation is above setpoint
                    // for active cooling, assume cooling system will be used instead of
                    // additional ventilation. Otherwise, use resultant operative temperature
                    // in calculation of space heating/cooling demand.
                    if temp_free_vent_extra > temp_setpnt_cool {
                        ach_cooling = ach_target;
                    } else {
                        temp_free = temp_free_vent_extra;
                    }

                    Some(ach_cooling)
                } else {
                    None
                };

                (
                    Some(ach_to_trigger_heating),
                    ach_cooling.unwrap_or(ach_target),
                )
            }
            AirChangesPerHourArgument::TargetAndWindowsOpen { ach_target, .. } => {
                (None, ach_target)
            }
            AirChangesPerHourArgument::Cooling { ach_cooling } => (None, ach_cooling),
        };

        Ok((temp_free, ach_cooling, ach_to_trigger_heating))
    }

    fn interp_heat_cool_demand(
        &self,
        delta_t_h: f64,
        temp_setpnt: f64,
        heat_cool_load_upper: f64,
        temp_free: f64,
        temp_upper: f64,
    ) -> anyhow::Result<f64> {
        if temp_upper - temp_free == 0.0 {
            bail!(
                "Divide-by-zero in calculation of heating/cooling demand.
            This may be caused by the specification of very low overall
            areal heat capacity of BuildingElements and/or very high thermal
            mass of WetDistribution."
            )
        }

        let heat_cool_load_unrestricted =
            heat_cool_load_upper * (temp_setpnt - temp_free) / (temp_upper - temp_free);

        // Convert from W to kWh
        let heat_cool_demand = heat_cool_load_unrestricted / WATTS_PER_KILOWATT as f64 * delta_t_h;

        Ok(heat_cool_demand)
    }

    /// Calculate heating and cooling demand in the zone for the current timestep
    ///
    /// According to the procedure in BS EN ISO 52016-1:2017, section 6.5.5.2, steps 1 to 4.
    ///
    /// Arguments:
    /// delta_t_h -- calculation timestep, in hours
    /// temp_ext_air -- temperature of the external air for the current timestep, in deg C
    /// gains_internal -- internal gains for the current timestep, in W
    /// gains_solar -- directly transmitted solar gains, in W
    /// frac_convective_heat -- convective fraction for heating
    /// frac_convective_cool -- convective fraction for cooling
    /// temp_setpnt_heat -- temperature setpoint for heating, in deg C
    /// temp_setpnt_cool -- temperature setpoint for cooling, in deg C
    /// ach_windows_open -- air changes per hour when all windows open
    /// ach_target -- air changes per hour required for ventilation requirement/target
    /// ach_cooling -- air changes per hour required to meet cooling requirement/target
    pub(crate) fn space_heat_cool_demand(
        &self,
        delta_t_h: f64,
        temp_ext_air: f64,
        gains_internal: f64,
        gains_solar: f64,
        frac_convective_heat: f64,
        frac_convective_cool: f64,
        temp_setpnt_heat: f64,
        temp_setpnt_cool: f64,
        avg_air_supply_temp: f64,
        gains_heat_cool_convective: Option<f64>,
        gains_heat_cool_radiative: Option<f64>,
        ach_args: AirChangesPerHourArgument,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64, f64, Option<f64>)> {
        let gains_heat_cool_convective = gains_heat_cool_convective.unwrap_or(0.0);
        let gains_heat_cool_radiative = gains_heat_cool_radiative.unwrap_or(0.0);

        if temp_setpnt_cool < temp_setpnt_heat {
            bail!("Cooling setpoint is below heating setpoint.");
        }

        let temp_setpnt_cool_vent = if let Some(control) = self.control.clone() {
            let temp_setpnt_cool_vent_response = control
                .setpnt(&simulation_time_iteration)
                .unwrap_or_else(|| kelvin_to_celsius(1.4e32).expect("Not below absolute zero"));
            if temp_setpnt_cool_vent_response < temp_setpnt_heat {
                bail!("Setpoint for additional ventilation is below heating setpoint.");
            }
            temp_setpnt_cool_vent_response
        } else {
            kelvin_to_celsius(1.4e32).expect("Not below absolute zero")
        };

        if temp_setpnt_cool_vent < temp_setpnt_heat {
            bail!("Cooling ventilation setpoint is below heating setpoint.");
        }

        // Calculate timestep in seconds
        let delta_t = delta_t_h * SECONDS_PER_HOUR as f64;

        // For calculation of demand, set heating/cooling gains to zero
        let gains_heat_cool = gains_heat_cool_convective + gains_heat_cool_radiative;
        let frac_conv_gains_heat_cool = if gains_heat_cool == 0.0 {
            0.0
        } else {
            gains_heat_cool_convective / gains_heat_cool
        };

        // Calculate node and internal air temperatures with heating/cooling gains of zero
        let (temp_vector_no_heat_cool, _) = self.calc_temperatures(
            delta_t,
            &self.temp_prev.read(),
            temp_ext_air,
            gains_internal,
            gains_solar,
            gains_heat_cool,
            frac_conv_gains_heat_cool,
            match ach_args {
                AirChangesPerHourArgument::Cooling { ach_cooling } => ach_cooling,
                AirChangesPerHourArgument::TargetAndWindowsOpen { ach_target, .. } => ach_target,
            },
            avg_air_supply_temp,
            simulation_time_iteration,
            self.print_heat_balance,
        )?;

        // Calculate internal operative temperature at free-floating conditions
        // i.e. with no heating/cooling
        let temp_int_air_free = temp_vector_no_heat_cool[self.zone_idx];
        let temp_free = match self.temp_setpnt_basis {
            ZoneTemperatureControlBasis::Operative => {
                self.temp_operative_by_temp_vector(&temp_vector_no_heat_cool)
            }
            ZoneTemperatureControlBasis::Air => temp_int_air_free,
        };

        let (temp_free, ach_cooling, ach_to_trigger_heating) = self
            .calc_cooling_potential_from_ventilation(
                delta_t,
                temp_ext_air,
                gains_internal,
                gains_solar,
                gains_heat_cool,
                frac_conv_gains_heat_cool,
                temp_setpnt_heat,
                temp_setpnt_cool,
                temp_setpnt_cool_vent,
                temp_free,
                temp_int_air_free,
                ach_args,
                avg_air_supply_temp,
                simulation_time_iteration,
            )?;

        // Determine relevant setpoint (if neither, then return space heating/cooling demand of zero)
        // Determine maximum heating/cooling
        let (temp_setpnt, heat_cool_load_upper, frac_convective) = if temp_free > temp_setpnt_cool
            && !is_close(temp_free, temp_setpnt_cool, 1e-10)
        {
            // Cooling
            // TODO (from Python) Implement eqn 26 "if max power available" case rather than just "otherwise" case?
            //      Could max. power be available at this point for all heating/cooling systems?
            (temp_setpnt_cool, -10. * self.area(), frac_convective_cool)
        } else if temp_free < temp_setpnt_heat && !is_close(temp_free, temp_setpnt_cool, 1e-10) {
            // Heating
            // TODO (from Python) Implement eqn 26 "if max power available" case rather than just "otherwise" case?
            //      Could max. power be available at this point for all heating/cooling systems?
            (temp_setpnt_heat, 10. * self.area(), frac_convective_heat)
        } else {
            return Ok((0.0, 0.0, ach_cooling, ach_to_trigger_heating));
        };

        // Calculate node and internal air temperatures with maximum heating/cooling
        let gains_heat_cool_upper = gains_heat_cool + heat_cool_load_upper;
        let frac_convective_heat_cool_upper = (gains_heat_cool * frac_conv_gains_heat_cool
            + heat_cool_load_upper * frac_convective)
            / gains_heat_cool_upper;
        let (temp_vector_upper_heat_cool, _) = self.calc_temperatures(
            delta_t,
            &self.temp_prev.read(),
            temp_ext_air,
            gains_internal,
            gains_solar,
            gains_heat_cool_upper,
            frac_convective_heat_cool_upper,
            ach_cooling,
            avg_air_supply_temp,
            simulation_time_iteration,
            self.print_heat_balance,
        )?;

        // Calculate internal operative temperature with maximum heating/cooling
        let temp_upper = match self.temp_setpnt_basis {
            ZoneTemperatureControlBasis::Operative => {
                self.temp_operative_by_temp_vector(&temp_vector_upper_heat_cool)
            }
            ZoneTemperatureControlBasis::Air => temp_vector_upper_heat_cool[self.zone_idx],
        };

        // Calculate heating (positive) or cooling (negative) required to reach setpoint
        let heat_cool_demand = self.interp_heat_cool_demand(
            delta_t_h,
            temp_setpnt,
            heat_cool_load_upper,
            temp_free,
            temp_upper,
        )?;

        let mut space_heat_demand = 0.0;
        let mut space_cool_demand = 0.0;
        if heat_cool_demand < 0.0 {
            space_cool_demand = heat_cool_demand;
        } else if heat_cool_demand > 0.0 {
            space_heat_demand = heat_cool_demand;
        }

        Ok((
            space_heat_demand,
            space_cool_demand,
            ach_cooling,
            ach_to_trigger_heating,
        ))
    }

    /// Update node and internal air temperatures for calculation of next timestep
    ///
    /// Arguments:
    /// * `delta_t`         - calculation timestep, in seconds
    /// * `temp_ext_air`    - temperature of external air, in deg C
    /// * `gains_internal`  - total internal heat gains, in W
    /// * `gains_solar`     - directly transmitted solar gains, in W
    /// * `gains_heat_cool` - gains from heating (positive) or cooling (negative), in W
    /// * `frac_convective` - convective fraction for heating/cooling (as appropriate)
    /// * `ach`             - air changes per hour
    pub fn update_temperatures(
        &self,
        delta_t: f64,
        temp_ext_air: f64,
        gains_internal: f64,
        gains_solar: f64,
        gains_heat_cool: f64,
        frac_convective: f64,
        ach: f64,
        avg_supply_temp: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<Option<HeatBalance>> {
        let (temp_prev, heat_balance_map) = self.calc_temperatures(
            delta_t,
            &self.temp_prev.read(),
            temp_ext_air,
            gains_internal,
            gains_solar,
            gains_heat_cool,
            frac_convective,
            ach,
            avg_supply_temp,
            simulation_time_iteration,
            self.print_heat_balance,
        )?;

        let _ = mem::replace(&mut *self.temp_prev.write(), temp_prev);

        Ok(heat_balance_map)
    }
    pub fn total_fabric_heat_loss(&self) -> f64 {
        self.building_elements
            .iter()
            .map(|nel| {
                let NamedBuildingElement { element, .. } = nel;
                element.fabric_heat_loss()
            })
            .sum::<f64>()
    }

    /// Return the total heat loss area, in m2
    pub fn total_heat_loss_area(&self) -> f64 {
        self.building_elements
            .iter()
            .filter_map(|el| {
                if matches!(
                    el.element.as_ref(),
                    BuildingElement::Opaque { .. }
                        | BuildingElement::Transparent { .. }
                        | BuildingElement::Ground { .. }
                        | BuildingElement::AdjacentUnconditionedSpaceSimple { .. }
                ) {
                    Some(el.element.area())
                } else {
                    None
                }
            })
            .sum::<f64>()
    }

    pub fn total_heat_capacity(&self) -> f64 {
        self.building_elements
            .iter()
            .map(|nel| {
                let NamedBuildingElement { element, .. } = nel;
                element.heat_capacity()
            })
            .sum::<f64>()
    }

    /// Return the total heat transfer coefficient for all
    /// thermal bridges in a zone, in W / K
    pub fn total_thermal_bridges(&self) -> f64 {
        self.tb_heat_trans_coeff
    }
}

/// Trait to give access to a constrained set of informational data from the zone.
pub(crate) trait SimpleZone: Send + Sync {
    fn temp_internal_air(&self) -> f64;
    fn area(&self) -> f64;
}

impl SimpleZone for Zone {
    fn temp_internal_air(&self) -> f64 {
        self.temp_internal_air()
    }

    fn area(&self) -> f64 {
        self.area()
    }
}

#[derive(Clone, Debug)]
pub(crate) struct ZoneTempInternalAir(pub(crate) Arc<Zone>);

// utility struct for when a zone needs to be injected into another module only to access temp_internal_air
impl ZoneTempInternalAir {
    pub(crate) fn as_fn(&self) -> TempInternalAirFn {
        let zone = self.0.clone();
        Arc::from(move || zone.temp_internal_air())
    }
}

const DELTA_T_H: u32 = 8760;
// hours in a non leap year
const DELTA_T: u32 = DELTA_T_H * SECONDS_PER_HOUR;

// # Assume default convective fraction for heating/cooling suggested in
// # BS EN ISO 52016-1:2017 Table B.11
const FRAC_CONVECTIVE: f64 = 0.4;

/// Close-enough port of numpy's isclose function for use here
///
/// Returns a single boolean where two arrays are element-wise all equal within a tolerance.
///
/// See [docs for numpy.isclose](https://docs.rs/is_close/latest/is_close/)
fn isclose(a: &Vec<f64>, b: &Vec<f64>, rtol: Option<f64>, atol: Option<f64>) -> bool {
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

#[derive(Debug, Error)]
pub(crate) enum SolverError {
    #[error("Detected that temperature calculation for a zone will not resolve: {}", .0)]
    SolverWillNotResolve(&'static str),
}

#[derive(Clone, Debug)]
struct NamedBuildingElement {
    pub name: String,
    pub element: Arc<BuildingElement>,
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

// There are cases where functions in the upstream Python code have some optional arguments to do with
// air changes per hour (ach), but EITHER ach_cooling OR both ach_target and ach_windows_open need to
// have values given (i.e. not None). The relation between these arguments is not explicit in the Python,
// but de facto the above both always holds and must hold in order for there not to be errors present
// such as arithmetic attempted with None values. This enum makes this relationship between arguments
// explicit in the Rust, so that we aren't forced to make assertions about the "someness" of a value
// when we need to use the underlying argument values.
#[derive(Clone, Copy, Debug)]
pub enum AirChangesPerHourArgument {
    Cooling {
        ach_cooling: f64,
    },
    TargetAndWindowsOpen {
        ach_target: f64,
        ach_windows_open: f64,
    },
}

impl AirChangesPerHourArgument {
    pub fn from_ach_cooling(ach_cooling: f64) -> Self {
        AirChangesPerHourArgument::Cooling { ach_cooling }
    }

    pub fn from_ach_target_windows_open(ach_target: f64, ach_windows_open: f64) -> Self {
        AirChangesPerHourArgument::TargetAndWindowsOpen {
            ach_target,
            ach_windows_open,
        }
    }
}

pub(crate) trait GainsLossesAsIndexMap {
    fn as_index_map(&self) -> IndexMap<String, f64>;
}

#[derive(Debug, FieldName, PartialEq)]
#[field_name_derive(Debug, Eq, Hash, PartialEq, Serialize_enum_str)]
pub struct HeatBalanceAirNode {
    pub solar_gains: f64,
    pub internal_gains: f64,
    pub heating_or_cooling_system_gains: f64,
    pub energy_to_change_internal_temperature: f64,
    pub thermal_bridges: f64,
    pub infiltration_ventilation: f64,
    pub fabric_heat_loss: f64,
}

impl From<HeatBalanceAirNodeFieldName> for String {
    fn from(value: HeatBalanceAirNodeFieldName) -> Self {
        serde_json::to_value(&value)
            .unwrap()
            .as_str()
            .unwrap()
            .into()
    }
}

impl GainsLossesAsIndexMap for HeatBalanceAirNode {
    fn as_index_map(&self) -> IndexMap<String, f64> {
        let Self {
            solar_gains,
            internal_gains,
            heating_or_cooling_system_gains,
            energy_to_change_internal_temperature,
            thermal_bridges,
            infiltration_ventilation,
            fabric_heat_loss,
        } = self;
        IndexMap::from([
            (HeatBalanceAirNodeFieldName::SolarGains.into(), *solar_gains),
            (
                HeatBalanceAirNodeFieldName::InternalGains.into(),
                *internal_gains,
            ),
            (
                HeatBalanceAirNodeFieldName::HeatingOrCoolingSystemGains.into(),
                *heating_or_cooling_system_gains,
            ),
            (
                HeatBalanceAirNodeFieldName::EnergyToChangeInternalTemperature.into(),
                *energy_to_change_internal_temperature,
            ),
            (
                HeatBalanceAirNodeFieldName::ThermalBridges.into(),
                *thermal_bridges,
            ),
            (
                HeatBalanceAirNodeFieldName::InfiltrationVentilation.into(),
                *infiltration_ventilation,
            ),
            (
                HeatBalanceAirNodeFieldName::FabricHeatLoss.into(),
                *fabric_heat_loss,
            ),
        ])
    }
}

#[derive(Debug, FieldName, PartialEq)]
#[field_name_derive(Debug, Eq, Hash, PartialEq, Serialize_enum_str)]
pub struct HeatBalanceInternalBoundary {
    pub fabric_int_air_convective: f64,
    pub fabric_int_sol: f64,
    pub fabric_int_int_gains: f64,
    pub fabric_int_heat_cool: f64,
}

impl From<HeatBalanceInternalBoundaryFieldName> for String {
    fn from(value: HeatBalanceInternalBoundaryFieldName) -> Self {
        serde_json::to_value(&value)
            .unwrap()
            .as_str()
            .unwrap()
            .into()
    }
}

impl GainsLossesAsIndexMap for HeatBalanceInternalBoundary {
    fn as_index_map(&self) -> IndexMap<String, f64> {
        let Self {
            fabric_int_air_convective,
            fabric_int_sol,
            fabric_int_int_gains,
            fabric_int_heat_cool,
        } = self;
        IndexMap::from([
            (
                HeatBalanceInternalBoundaryFieldName::FabricIntAirConvective.into(),
                *fabric_int_air_convective,
            ),
            (
                HeatBalanceInternalBoundaryFieldName::FabricIntSol.into(),
                *fabric_int_sol,
            ),
            (
                HeatBalanceInternalBoundaryFieldName::FabricIntIntGains.into(),
                *fabric_int_int_gains,
            ),
            (
                HeatBalanceInternalBoundaryFieldName::FabricIntHeatCool.into(),
                *fabric_int_heat_cool,
            ),
        ])
    }
}

#[derive(Debug, FieldName, PartialEq)]
#[field_name_derive(Debug, Eq, Hash, PartialEq, Serialize_enum_str)]
pub struct HeatBalanceExternalBoundary {
    pub solar_gains: f64,
    pub internal_gains: f64,
    pub heating_or_cooling_system_gains: f64,
    pub thermal_bridges: f64,
    pub infiltration_ventilation: f64,
    pub fabric_ext_air_convective: f64,
    pub fabric_ext_air_radiative: f64,
    pub fabric_ext_sol: f64,
    pub fabric_ext_sky: f64,
    pub opaque_fabric_ext: f64,
    pub transparent_fabric_ext: f64,
    pub ground_fabric_ext: f64,
    pub ztc_fabric_ext: f64,
    pub ztu_fabric_ext: f64,
}

impl From<HeatBalanceExternalBoundaryFieldName> for String {
    fn from(value: HeatBalanceExternalBoundaryFieldName) -> Self {
        serde_json::to_value(&value)
            .unwrap()
            .as_str()
            .unwrap()
            .into()
    }
}

impl GainsLossesAsIndexMap for HeatBalanceExternalBoundary {
    fn as_index_map(&self) -> IndexMap<String, f64> {
        let Self {
            solar_gains,
            internal_gains,
            heating_or_cooling_system_gains,
            thermal_bridges,
            infiltration_ventilation,
            fabric_ext_air_convective,
            fabric_ext_air_radiative,
            fabric_ext_sol,
            fabric_ext_sky,
            opaque_fabric_ext,
            transparent_fabric_ext,
            ground_fabric_ext,
            ztc_fabric_ext,
            ztu_fabric_ext,
        } = self;
        IndexMap::from([
            (
                HeatBalanceExternalBoundaryFieldName::SolarGains.into(),
                *solar_gains,
            ),
            (
                HeatBalanceExternalBoundaryFieldName::InternalGains.into(),
                *internal_gains,
            ),
            (
                HeatBalanceExternalBoundaryFieldName::HeatingOrCoolingSystemGains.into(),
                *heating_or_cooling_system_gains,
            ),
            (
                HeatBalanceExternalBoundaryFieldName::ThermalBridges.into(),
                *thermal_bridges,
            ),
            (
                HeatBalanceExternalBoundaryFieldName::InfiltrationVentilation.into(),
                *infiltration_ventilation,
            ),
            (
                HeatBalanceExternalBoundaryFieldName::FabricExtAirConvective.into(),
                *fabric_ext_air_convective,
            ),
            (
                HeatBalanceExternalBoundaryFieldName::FabricExtAirRadiative.into(),
                *fabric_ext_air_radiative,
            ),
            (
                HeatBalanceExternalBoundaryFieldName::FabricExtSol.into(),
                *fabric_ext_sol,
            ),
            (
                HeatBalanceExternalBoundaryFieldName::FabricExtSky.into(),
                *fabric_ext_sky,
            ),
            (
                HeatBalanceExternalBoundaryFieldName::OpaqueFabricExt.into(),
                *opaque_fabric_ext,
            ),
            (
                HeatBalanceExternalBoundaryFieldName::TransparentFabricExt.into(),
                *transparent_fabric_ext,
            ),
            (
                HeatBalanceExternalBoundaryFieldName::GroundFabricExt.into(),
                *ground_fabric_ext,
            ),
            (
                HeatBalanceExternalBoundaryFieldName::ZtcFabricExt.into(),
                *ztc_fabric_ext,
            ),
            (
                HeatBalanceExternalBoundaryFieldName::ZtuFabricExt.into(),
                *ztu_fabric_ext,
            ),
        ])
    }
}

#[derive(Debug, FieldName, PartialEq)]
#[field_name_derive(Debug, Eq, Hash, PartialEq, Serialize_enum_str)]
pub struct HeatBalance {
    pub air_node: HeatBalanceAirNode,
    pub internal_boundary: HeatBalanceInternalBoundary,
    pub external_boundary: HeatBalanceExternalBoundary,
}

impl HeatBalance {
    pub(crate) fn as_index_map(&self) -> IndexMap<HeatBalanceFieldName, IndexMap<String, f64>> {
        let Self {
            air_node,
            internal_boundary,
            external_boundary,
        } = self;
        IndexMap::from([
            (HeatBalanceFieldName::AirNode, (*air_node).as_index_map()),
            (
                HeatBalanceFieldName::InternalBoundary,
                (*internal_boundary).as_index_map(),
            ),
            (
                HeatBalanceFieldName::ExternalBoundary,
                (*external_boundary).as_index_map(),
            ),
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::space_heat_demand::building_element::{
        BuildingElementAdjacentConditionedSpace, BuildingElementAdjacentUnconditionedSpaceSimple,
        BuildingElementGround, BuildingElementOpaque, BuildingElementTransparent,
    };
    use crate::core::space_heat_demand::thermal_bridge::ThermalBridge;
    use crate::core::space_heat_demand::ventilation::{Vent, Window};
    use crate::core::units::DAYS_IN_MONTH;
    use crate::corpus::CompletedVentilationLeaks;
    use crate::external_conditions::{DaylightSavingsConfig, ExternalConditions};
    use crate::hem_core::external_conditions::ShadingSegment;
    use crate::input::{
        FloorData, MassDistributionClass, TerrainClass, VentilationShieldClass, WindShieldLocation,
        WindowPart,
    };
    use crate::simulation_time::{SimulationTime, HOURS_IN_DAY};
    use approx::assert_relative_eq;
    use indexmap::IndexMap;
    use pretty_assertions::assert_eq;
    use rstest::{fixture, rstest};
    use std::collections::HashMap;

    const BASE_AIR_TEMPS: [f64; 24] = [
        0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 7.5, 10.0, 12.5, 15.0, 19.5, 17.0,
        15.0, 12.0, 10.0, 7.0, 5.0, 3.0, 1.0,
    ];
    const BASE_WIND_SPEEDS: [f64; 24] = [
        4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.7, 5.4, 5.6, 5.3, 5.1, 4.8, 4.7, 4.6, 4.5, 4.2,
        4.9, 4.3, 4.4, 4.5, 4.3, 4.6,
    ];
    const BASE_WIND_DIRECTIONS: [f64; 24] = [
        300., 250., 220., 180., 150., 120., 100., 80., 60., 40., 20., 10., 50., 100., 140., 190.,
        200., 320., 330., 340., 350., 355., 315., 5.,
    ];

    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 4., 1.)
    }

    fn external_conditions(simulation_time: SimulationTime) -> Arc<ExternalConditions> {
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

        let mut wind_speeds = Vec::with_capacity(8760);
        for (speeds, days_in_month) in [
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
                speeds
                    .iter()
                    .cloned()
                    .cycle()
                    .take((days_in_month * HOURS_IN_DAY) as usize)
                    .collect::<Vec<f64>>()
                    .as_slice(),
            );
        }

        let wind_direction_day_jan = BASE_WIND_DIRECTIONS;
        let wind_direction_day_feb = BASE_WIND_DIRECTIONS.map(|d| d - 1.);
        let wind_direction_day_mar = BASE_WIND_DIRECTIONS.map(|d| d - 2.);
        let wind_direction_day_apr = BASE_WIND_DIRECTIONS.map(|d| d - 3.);
        let wind_direction_day_may = BASE_WIND_DIRECTIONS.map(|d| d - 4.);
        let wind_direction_day_jun = BASE_WIND_DIRECTIONS.map(|d| d + 1.);
        let wind_direction_day_jul = BASE_WIND_DIRECTIONS.map(|d| d + 2.);
        let wind_direction_day_aug = BASE_WIND_DIRECTIONS.map(|d| d + 3.);
        let wind_direction_day_sep = BASE_WIND_DIRECTIONS.map(|d| d + 4.);
        let wind_direction_day_oct = BASE_WIND_DIRECTIONS.map(|d| d - 5.);
        let wind_direction_day_nov = BASE_WIND_DIRECTIONS.map(|d| d + 5.);
        let wind_direction_day_dec = BASE_WIND_DIRECTIONS.map(|d| d - 0.);

        let mut wind_directions = Vec::with_capacity(8760);
        for (directions, days_in_month) in [
            (wind_direction_day_jan, DAYS_IN_MONTH[0]),
            (wind_direction_day_feb, DAYS_IN_MONTH[1]),
            (wind_direction_day_mar, DAYS_IN_MONTH[2]),
            (wind_direction_day_apr, DAYS_IN_MONTH[3]),
            (wind_direction_day_may, DAYS_IN_MONTH[4]),
            (wind_direction_day_jun, DAYS_IN_MONTH[5]),
            (wind_direction_day_jul, DAYS_IN_MONTH[6]),
            (wind_direction_day_aug, DAYS_IN_MONTH[7]),
            (wind_direction_day_sep, DAYS_IN_MONTH[8]),
            (wind_direction_day_oct, DAYS_IN_MONTH[9]),
            (wind_direction_day_nov, DAYS_IN_MONTH[10]),
            (wind_direction_day_dec, DAYS_IN_MONTH[11]),
        ] {
            wind_directions.extend_from_slice(
                directions
                    .iter()
                    .cloned()
                    .cycle()
                    .take((days_in_month * HOURS_IN_DAY) as usize)
                    .collect::<Vec<f64>>()
                    .as_slice(),
            );
        }

        let shading_segments = vec![
            ShadingSegment {
                start: 180.,
                end: 135.,
                ..Default::default()
            },
            ShadingSegment {
                start: 135.,
                end: 90.,
                ..Default::default()
            },
            ShadingSegment {
                start: 90.,
                end: 45.,
                ..Default::default()
            },
            ShadingSegment {
                start: 45.,
                end: 0.,
                ..Default::default()
            },
            ShadingSegment {
                start: 0.,
                end: -45.,
                ..Default::default()
            },
            ShadingSegment {
                start: -45.,
                end: -90.,
                ..Default::default()
            },
            ShadingSegment {
                start: -90.,
                end: -135.,
                ..Default::default()
            },
            ShadingSegment {
                start: -135.,
                end: -180.,
                ..Default::default()
            },
        ];

        Arc::new(ExternalConditions::new(
            &simulation_time.iter(),
            air_temps,
            wind_speeds,
            wind_directions,
            vec![333.0, 610.0, 572.0, 420.0, 0.0, 10.0, 90.0, 275.0],
            vec![420.0, 750.0, 425.0, 500.0, 0.0, 40.0, 0.0, 388.0],
            vec![0.2; 8760],
            51.42,
            -0.75,
            0,
            0,
            None,
            1.0,
            None,
            Some(DaylightSavingsConfig::NotApplicable),
            false,
            false,
            shading_segments.into(),
        ))
    }

    fn zone(
        thermal_bridging: ThermalBridging,
        control: Option<Arc<dyn ControlBehaviour>>,
    ) -> anyhow::Result<Zone> {
        let simulation_time = simulation_time();
        let external_conditions = external_conditions(simulation_time);
        // Create objects for the different building elements in the zone
        let be_opaque_i = BuildingElement::Opaque(BuildingElementOpaque::new(
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
            external_conditions.clone(),
        ));
        let be_opaque_d = BuildingElement::Opaque(BuildingElementOpaque::new(
            26.,
            true,
            45.,
            0.55,
            0.33,
            16000.0,
            MassDistributionClass::D,
            0.,
            0.,
            2.,
            10.,
            external_conditions.clone(),
        ));
        let be_ztc = BuildingElement::AdjacentConditionedSpace(
            BuildingElementAdjacentConditionedSpace::new(
                22.5,
                135.,
                0.50,
                18000.0,
                MassDistributionClass::E,
                external_conditions.clone(),
            ),
        );
        let be_ground_floor_data = FloorData::SuspendedFloor {
            height_upper_surface: 0.5,
            thermal_transmission_walls: 0.5,
            area_per_perimeter_vent: 0.01,
            shield_fact_location: WindShieldLocation::Sheltered,
            thermal_resistance_of_insulation: 7.,
        };
        let be_ground = BuildingElement::Ground(
            BuildingElementGround::new(
                25.0,
                25.0,
                90.,
                1.33,
                0.2,
                17000.0,
                MassDistributionClass::IE,
                &be_ground_floor_data,
                0.3,
                20.0,
                0.7,
                external_conditions.clone(),
            )
            .unwrap(),
        );
        let be_transparent = BuildingElement::Transparent(BuildingElementTransparent::new(
            90.,
            0.4,
            180.,
            0.75,
            0.25,
            1.,
            1.25,
            4.,
            Default::default(),
            Default::default(),
            external_conditions.clone(),
        ));
        let be_ztu = BuildingElement::AdjacentUnconditionedSpaceSimple(
            BuildingElementAdjacentUnconditionedSpaceSimple::new(
                30.,
                130.,
                0.50,
                0.6,
                18000.0,
                MassDistributionClass::E,
                external_conditions.clone(),
            ),
        );

        // Put building element objects in a list that can be iterated over
        let be_objs = IndexMap::from([
            ("be_opaque_i".into(), be_opaque_i.into()),
            ("be_opaque_d".into(), be_opaque_d.into()),
            ("be_ztc".into(), be_ztc.into()),
            ("be_ground".into(), be_ground.into()),
            ("be_transparent".into(), be_transparent.into()),
            ("be_ztu".into(), be_ztu.into()),
        ]);

        let temp_ext_air_init = 2.2;
        let temp_setpnt_init = 21.;
        let temp_setpnt_basis = ZoneTemperatureControlBasis::Air;

        Zone::new(
            80.,
            250.,
            be_objs,
            thermal_bridging,
            Arc::new(infiltration_ventilation()),
            temp_ext_air_init,
            temp_setpnt_init,
            temp_setpnt_basis,
            control,
            true,
            &simulation_time.iter(),
        )
    }

    #[fixture]
    fn thermal_bridging_objects() -> ThermalBridging {
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
        ThermalBridging::Bridges(IndexMap::from([
            ("tb_linear_1".into(), tb_linear_1),
            ("tb_linear_2".into(), tb_linear_2),
            ("tb_point".into(), tb_point),
        ]))
    }

    fn infiltration_ventilation() -> InfiltrationVentilation {
        let window_part_list = vec![WindowPart {
            mid_height_air_flow_path: 1.5,
        }];
        let window = Window::new(1.6, 1., 3., window_part_list, Some(0.), 0., 30., None, 2.5);
        let windows = HashMap::from([("window 0".to_string(), window)]);

        let vent = Vent::new(1.5, 100.0, 20.0, 180.0, 60.0, 30.0, 2.5);
        let vents = HashMap::from([("vent 1".to_string(), vent)]);

        let leaks = CompletedVentilationLeaks {
            ventilation_zone_height: 6.,
            test_pressure: 50.,
            test_result: 1.2,
            area_roof: 25.0,
            area_facades: 85.0,
            env_area: 220.,
            altitude: 30.,
        };

        InfiltrationVentilation::new(
            true,
            VentilationShieldClass::Normal,
            &TerrainClass::OpenField,
            20.0,
            windows.into_values().collect(),
            vents.into_values().collect(),
            leaks,
            vec![],
            vec![],
            vec![],
            Default::default(),
            false,
            30.0,
            250.0,
            2.5,
        )
    }

    #[test]
    fn test_init_single_thermal_bridging_value() {
        let thermal_bridging = ThermalBridging::Number(4.);
        let zone = zone(thermal_bridging, None).unwrap();

        assert_eq!(zone.tb_heat_trans_coeff, 4.)
    }

    #[rstest]
    fn test_setpnt_init(thermal_bridging_objects: ThermalBridging) {
        assert_eq!(
            zone(thermal_bridging_objects, None).unwrap().setpnt_init(),
            21.
        );
    }

    #[rstest]
    fn test_area(thermal_bridging_objects: ThermalBridging) {
        assert_eq!(zone(thermal_bridging_objects, None).unwrap().area(), 80.);
    }

    #[rstest]
    fn test_gains_solar(thermal_bridging_objects: ThermalBridging) {
        let simulation_time_iteration = simulation_time().iter().next().unwrap();
        assert_relative_eq!(
            zone(thermal_bridging_objects, None)
                .unwrap()
                .gains_solar(simulation_time_iteration),
            -2154.583062153444
        );
    }

    #[rstest]
    fn test_volume(thermal_bridging_objects: ThermalBridging) {
        assert_eq!(zone(thermal_bridging_objects, None).unwrap().volume(), 250.);
    }

    #[rstest]
    fn test_total_fabric_heat_loss(thermal_bridging_objects: ThermalBridging) {
        assert_relative_eq!(
            zone(thermal_bridging_objects, None)
                .unwrap()
                .total_fabric_heat_loss(),
            181.99557093947166,
            max_relative = 1e-2
        );
    }

    #[rstest]
    fn test_total_heat_capacity(thermal_bridging_objects: ThermalBridging) {
        assert_eq!(
            zone(thermal_bridging_objects, None)
                .unwrap()
                .total_heat_capacity(),
            2166.
        );
    }

    #[rstest]
    fn test_total_heat_loss_area(thermal_bridging_objects: ThermalBridging) {
        assert_eq!(
            zone(thermal_bridging_objects, None)
                .unwrap()
                .total_heat_loss_area(),
            106.
        );
    }

    #[rstest]
    fn test_total_thermal_bridges(thermal_bridging_objects: ThermalBridging) {
        assert_relative_eq!(
            zone(thermal_bridging_objects, None)
                .unwrap()
                .total_thermal_bridges(),
            4.3,
            max_relative = 1e-2
        );
    }

    #[rstest]
    fn test_temp_operative(thermal_bridging_objects: ThermalBridging) {
        assert_relative_eq!(
            zone(thermal_bridging_objects, None)
                .unwrap()
                .temp_operative(),
            18.92809674634258,
        );
    }

    #[rstest]
    fn test_temp_internal_air(thermal_bridging_objects: ThermalBridging) {
        assert_relative_eq!(
            zone(thermal_bridging_objects, None)
                .unwrap()
                .temp_internal_air(),
            20.999999999999996,
        );
    }

    #[rstest]
    #[ignore = "TODO: Python uses np.linalg.solve, what would be valuable in Rust?"]
    fn test_fast_solver() {}

    fn maps_approx_equal(
        actual: &IndexMap<String, f64>,
        expected: &IndexMap<String, f64>,
        tol: f64,
    ) -> bool {
        if actual.len() != expected.len() {
            return false;
        }

        for (key, &expected_value) in expected.iter() {
            match actual.get(key) {
                Some(&actual_value) => {
                    if (expected_value - actual_value).abs() > tol {
                        eprintln!(
                            "Field '{}' differs. Expected {}, got {}",
                            key, expected_value, actual_value
                        );
                        return false;
                    }
                }
                None => {
                    eprintln!("Could not find expected key '{}'", key);
                    return false;
                }
            }
        }

        true
    }

    #[rstest]
    fn test_update_temperatures(thermal_bridging_objects: ThermalBridging) {
        let expected_heat_balance = HeatBalance {
            air_node: HeatBalanceAirNode {
                solar_gains: 22.,
                internal_gains: 80.,
                heating_or_cooling_system_gains: 0.,
                energy_to_change_internal_temperature: 1617.0201164499708,
                thermal_bridges: -31.655330373346523,
                infiltration_ventilation: -247.6853738767846,
                fabric_heat_loss: -1439.6794121998396,
            },
            internal_boundary: HeatBalanceInternalBoundary {
                fabric_int_air_convective: 1439.679412199842,
                fabric_int_sol: 198.,
                fabric_int_int_gains: 120.,
                fabric_int_heat_cool: 0.,
            },
            external_boundary: HeatBalanceExternalBoundary {
                solar_gains: 220.,
                internal_gains: 200.,
                heating_or_cooling_system_gains: 0.,
                thermal_bridges: -31.655330373346523,
                infiltration_ventilation: -247.6853738767846,
                fabric_ext_air_convective: 465.59408303760006,
                fabric_ext_air_radiative: 355.39880327150297,
                fabric_ext_sol: -3204.6068977063023,
                fabric_ext_sky: -1124.4913565980596,
                opaque_fabric_ext: -2118.880733022755,
                transparent_fabric_ext: -137.91628674680447,
                ground_fabric_ext: -816.8733439646372,
                ztc_fabric_ext: 0.,
                ztu_fabric_ext: -434.43500426106175,
            },
        }
        .as_index_map();
        let simulatio_time_iteration = simulation_time().iter().next().unwrap();
        let actual_heat_balance = zone(thermal_bridging_objects, None)
            .unwrap()
            .update_temperatures(
                1800.,
                10.,
                200.,
                220.,
                0.,
                1.,
                0.4,
                10.,
                simulatio_time_iteration,
            )
            .unwrap()
            .unwrap()
            .as_index_map();

        for (key, actual_value) in actual_heat_balance {
            let expected_value = expected_heat_balance.get(&key).unwrap();
            assert!(maps_approx_equal(&actual_value, expected_value, 1e-8));
        }
    }

    #[test]
    fn test_ach_req_to_reach_temperature() {
        let ach = Zone::ach_req_to_reach_temperature(21., 0.08, 48., 9.9, 8.8, 9.8, 8.2, 8.1);
        assert_relative_eq!(ach, 0.08);

        let ach = Zone::ach_req_to_reach_temperature(21., 0.08, 48., 9.9, 12., 9.8, 8.2, 8.1);
        assert_relative_eq!(ach, 0.08);

        let ach = Zone::ach_req_to_reach_temperature(21., 0.08, 48., 9.9, 8.8, 12., 8.2, 10.);
        assert_relative_eq!(ach, 21.65371789094187);

        let ach = Zone::ach_req_to_reach_temperature(5., 0.08, 48., 9.9, 8.8, 12., 11., 10.);
        assert_relative_eq!(ach, 48.);
    }

    #[rstest]
    fn test_calc_cooling_potential_from_ventilation_1(thermal_bridging_objects: ThermalBridging) {
        let simulation_time_iteration = simulation_time().iter().next().unwrap();
        let (temp_free, ach_cooling, ach_to_trigger_heating) = zone(thermal_bridging_objects, None)
            .unwrap()
            .calc_cooling_potential_from_ventilation(
                1800.0,
                17.8,
                6.6,
                0.0,
                0.0,
                0.0,
                21.0,
                24.0,
                22.0,
                20.0,
                20.0,
                AirChangesPerHourArgument::from_ach_target_windows_open(0.13, 0.16),
                17.8,
                simulation_time_iteration,
            )
            .unwrap();

        assert_relative_eq!(temp_free, 20.0);
        assert_relative_eq!(ach_cooling, 0.13);
        assert_relative_eq!(ach_to_trigger_heating.unwrap(), 0.13105245458346534);
    }

    #[rstest]
    fn test_calc_cooling_potential_from_ventilation_2(thermal_bridging_objects: ThermalBridging) {
        let simulation_time_iteration = simulation_time().iter().next().unwrap();
        let (temp_free, ach_cooling, ach_to_trigger_heating) = zone(thermal_bridging_objects, None)
            .unwrap()
            .calc_cooling_potential_from_ventilation(
                1800.0,
                17.8,
                6.6,
                0.0,
                0.0,
                0.0,
                21.0,
                24.0,
                18.0,
                20.0,
                20.0,
                AirChangesPerHourArgument::from_ach_target_windows_open(0.13, 0.16),
                17.8,
                simulation_time_iteration,
            )
            .unwrap();

        assert_relative_eq!(temp_free, 17.52067994302452, max_relative = 1e-8);
        assert_relative_eq!(ach_cooling, 0.13);
        assert_relative_eq!(ach_to_trigger_heating.unwrap(), 0.13105245458346534);
    }

    #[rstest]
    fn test_calc_cooling_potential_from_ventilation_3(thermal_bridging_objects: ThermalBridging) {
        let simulation_time_iteration = simulation_time().iter().next().unwrap();
        let (temp_free, ach_cooling, ach_to_trigger_heating) = zone(thermal_bridging_objects, None)
            .unwrap()
            .calc_cooling_potential_from_ventilation(
                1800.0,
                17.8,
                6.6,
                0.0,
                0.0,
                0.0,
                21.0,
                15.0,
                18.0,
                20.0,
                20.0,
                AirChangesPerHourArgument::from_ach_target_windows_open(0.13, 0.16),
                17.8,
                simulation_time_iteration,
            )
            .unwrap();

        assert_relative_eq!(temp_free, 20.);
        assert_relative_eq!(ach_cooling, 0.13);
        assert_relative_eq!(ach_to_trigger_heating.unwrap(), 0.13105245458346534);
    }

    #[rstest]
    fn test_calc_cooling_potential_from_ventilation_4(thermal_bridging_objects: ThermalBridging) {
        let simulation_time_iteration = simulation_time().iter().next().unwrap();
        let mut zone = zone(thermal_bridging_objects, None).unwrap();
        zone.temp_setpnt_basis = ZoneTemperatureControlBasis::Operative;
        let (temp_free, ach_cooling, ach_to_trigger_heating) = zone
            .calc_cooling_potential_from_ventilation(
                1800.0,
                17.8,
                6.6,
                0.0,
                0.0,
                0.0,
                21.0,
                24.0,
                22.0,
                19.9,
                20.0,
                AirChangesPerHourArgument::from_ach_target_windows_open(0.13, 0.16),
                17.8,
                simulation_time_iteration,
            )
            .unwrap();

        assert_relative_eq!(temp_free, 19.9);
        assert_relative_eq!(ach_cooling, 0.13);
        assert_relative_eq!(ach_to_trigger_heating.unwrap(), 0.1306962441148573);
    }

    #[rstest]
    fn test_calc_cooling_potential_from_ventilation_5(thermal_bridging_objects: ThermalBridging) {
        let simulation_time_iteration = simulation_time().iter().next().unwrap();
        let mut zone = zone(thermal_bridging_objects, None).unwrap();
        zone.temp_setpnt_basis = ZoneTemperatureControlBasis::Operative;
        let (temp_free, ach_cooling, ach_to_trigger_heating) = zone
            .calc_cooling_potential_from_ventilation(
                1800.0,
                17.8,
                6.6,
                0.0,
                0.0,
                0.0,
                21.0,
                24.0,
                18.0,
                19.9,
                20.0,
                AirChangesPerHourArgument::from_ach_target_windows_open(0.13, 0.16),
                17.8,
                simulation_time_iteration,
            )
            .unwrap();

        assert_relative_eq!(temp_free, 15.144632024928118, max_relative = 1e-8);
        assert_relative_eq!(ach_cooling, 0.13);
        assert_relative_eq!(ach_to_trigger_heating.unwrap(), 0.1306962441148573);
    }

    #[rstest]
    fn test_calc_cooling_potential_from_ventilation_6(thermal_bridging_objects: ThermalBridging) {
        let simulation_time_iteration = simulation_time().iter().next().unwrap();
        let mut zone = zone(thermal_bridging_objects, None).unwrap();
        zone.temp_setpnt_basis = ZoneTemperatureControlBasis::Operative;
        let (temp_free, ach_cooling, ach_to_trigger_heating) = zone
            .calc_cooling_potential_from_ventilation(
                1800.0,
                17.8,
                6.6,
                0.0,
                0.0,
                0.0,
                21.0,
                15.0,
                18.0,
                19.9,
                20.0,
                AirChangesPerHourArgument::from_ach_target_windows_open(0.13, 0.16),
                17.8,
                simulation_time_iteration,
            )
            .unwrap();

        assert_relative_eq!(temp_free, 19.9);
        assert_relative_eq!(ach_cooling, 0.13);
        assert_relative_eq!(ach_to_trigger_heating.unwrap(), 0.1306962441148573);
    }

    #[rstest]
    fn test_interp_heat_cool_demand(thermal_bridging_objects: ThermalBridging) {
        let mut zone = zone(thermal_bridging_objects, None).unwrap();

        zone.temp_setpnt_basis = ZoneTemperatureControlBasis::Air;
        let heat_cool_demand = zone.interp_heat_cool_demand(0.5, 20., 4000., 18., 21.2);

        assert_relative_eq!(heat_cool_demand.unwrap(), 1.2500000000000002);

        zone.temp_setpnt_basis = ZoneTemperatureControlBasis::Operative;
        let heat_cool_demand = zone.interp_heat_cool_demand(0.5, 20., 4000., 20., 19.);

        assert_eq!(heat_cool_demand.unwrap(), 0.);
    }

    #[rstest]
    /// Cases where temp upper and temp free are the same
    fn test_interp_heat_cool_demand_invalid(thermal_bridging_objects: ThermalBridging) {
        let mut zone = zone(thermal_bridging_objects, None).unwrap();

        zone.temp_setpnt_basis = ZoneTemperatureControlBasis::Operative;
        let heat_cool_demand = zone.interp_heat_cool_demand(0.5, 20., 4000., 19., 19.);

        assert!(heat_cool_demand.is_err());

        zone.temp_setpnt_basis = ZoneTemperatureControlBasis::Air;
        let heat_cool_demand = zone.interp_heat_cool_demand(0.5, 20., 4000., 18., 18.);

        assert!(heat_cool_demand.is_err());
    }

    #[rstest]
    /// errors when cooling setpoint is below heating setpoint
    fn test_space_heat_cool_demand_1(thermal_bridging_objects: ThermalBridging) {
        let simulation_time_iteration = simulation_time().iter().next().unwrap();
        let zone = zone(thermal_bridging_objects, None).unwrap();
        let space_heat_cool_demand = zone.space_heat_cool_demand(
            0.5,
            2.8,
            13.5,
            9.1,
            0.4,
            0.95,
            24.0,
            21.0,
            2.8,
            Some(0.0),
            Some(0.0),
            AirChangesPerHourArgument::from_ach_target_windows_open(0.14, 0.17),
            simulation_time_iteration,
        );
        assert!(space_heat_cool_demand.is_err());
    }

    #[rstest]
    fn test_space_heat_cool_demand_2(thermal_bridging_objects: ThermalBridging) {
        let simulation_time_iteration = simulation_time().iter().next().unwrap();
        let zone = zone(thermal_bridging_objects, None).unwrap();
        let (space_heat_demand, space_cool_demand, ach_cooling, ach_to_trigger_heating) = zone
            .space_heat_cool_demand(
                0.5,
                2.8,
                13.5,
                9.1,
                0.4,
                0.95,
                21.0,
                24.0,
                2.8,
                Some(0.0),
                Some(0.0),
                AirChangesPerHourArgument::from_ach_target_windows_open(0.14, 0.17),
                simulation_time_iteration,
            )
            .unwrap();

        assert_relative_eq!(space_heat_demand, 2.1541345392835387);
        assert_eq!(space_cool_demand, 0.);
        assert_relative_eq!(ach_cooling, 0.14);
        assert_relative_eq!(ach_to_trigger_heating.unwrap(), 0.14);
    }

    struct FakeControl(f64);

    impl ControlBehaviour for FakeControl {
        fn setpnt(&self, _simulation_time_iteration: &SimulationTimeIteration) -> Option<f64> {
            Some(self.0)
        }
    }

    #[rstest]
    fn test_space_heat_cool_demand_3(thermal_bridging_objects: ThermalBridging) {
        let simulation_time_iteration = simulation_time().iter().next().unwrap();
        let fake_control = Arc::new(FakeControl(25.));
        let zone = zone(thermal_bridging_objects, Some(fake_control)).unwrap();
        let (space_heat_demand, space_cool_demand, ach_cooling, ach_to_trigger_heating) = zone
            .space_heat_cool_demand(
                0.5,
                2.8,
                13.5,
                9.1,
                0.4,
                0.95,
                21.0,
                24.0,
                2.8,
                Some(0.0),
                Some(0.0),
                AirChangesPerHourArgument::from_ach_target_windows_open(0.14, 0.17),
                simulation_time_iteration,
            )
            .unwrap();

        assert_relative_eq!(space_heat_demand, 2.1541345392835387);
        assert_eq!(space_cool_demand, 0.);
        assert_relative_eq!(ach_cooling, 0.14);
        assert_relative_eq!(ach_to_trigger_heating.unwrap(), 0.14);
    }

    #[rstest]
    fn test_space_heat_cool_demand_4(thermal_bridging_objects: ThermalBridging) {
        let fake_control = Arc::new(FakeControl(20.));
        let zone = zone(thermal_bridging_objects, Some(fake_control));
        // In the Python test the error is raised when space_cool_heat_demand is called.
        // In Rust we get the error earlier, when creating the Zone with the fake Control.
        // This is because in the Python test set up, the Control object is None when Zone
        // is created and is only added just before space_cool_heat_demand is called.
        // Because of the above, this test deviates a bit from the Python and asserts that
        // Zone itself returns an Error.

        assert!(zone.is_err());
    }

    #[test]
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
            None,
        ));
        assert!(!isclose(
            &vec![1e10, 1e-8],
            &vec![1.0001e10, 1e-9],
            None,
            None,
        ));
        assert!(!isclose(&vec![1e-8, 1e-7], &vec![0.0, 0.0], None, None),);
        assert!(!isclose(
            &vec![1e-100, 1e-7],
            &vec![0.0, 0.0],
            None,
            Some(0.0),
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
