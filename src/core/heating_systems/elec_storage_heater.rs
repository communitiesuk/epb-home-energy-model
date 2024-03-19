use crate::core::controls::time_control::{
    ControlBehaviour, SetpointTimeControl, ToUChargeControl,
};
use crate::core::space_heat_demand::zone::InternalAirTempAccess;
use crate::core::units::WATTS_PER_KILOWATT;
use crate::input::ElectricStorageHeaterAirFlowType;
use crate::simulation_time::SimulationTimeIteration;
use interp::interp;
use mathru::algebra::linear::vector::Vector;
use mathru::analysis::differential_equation::ordinary::problem::Euler;
use mathru::analysis::differential_equation::ordinary::solver::implicit::BDF;
use mathru::analysis::differential_equation::ordinary::ImplicitInitialValueProblemBuilder;
use mathru::vector;
use nalgebra::SVector;
use std::sync::Arc;

/// This module provides object(s) to model the behaviour of electric storage heaters.

const C_P: f64 = 1.0054; // J/kg/K air specific heat

// The value of R (resistance of air gap) depends on several factors such as the thickness of the air gap, the
// temperature difference across the gap, and the air flow characteristics.
// Here are some typical values for R per inch of air gap thickness:
// Still air: 0.17 m²·K/W
// Air movement (average): 0.07 m²·K/W
// Air movement (high): 0.04 m²·K/W
// Note: These values are for a temperature difference of 24°F (14°C) and a pressure difference of 1 inch of water
// column. The values will change with changes in temperature and pressure differences.

// 0.17 Thermal resistance of the air layer between the insulation and the wall of the device.
// Assuming no change of resistance with temp diff. (Rough approx)
const R_AIR_OFF: f64 = 0.17;
const R_AIR_ON: f64 = 0.07;

// labs_test for electric storage heater reaching 300 degC
// This represents the temperature difference between the core and the room on the first column
// and the fraction of air flow relating to the nominal as defined above on the second column
const LABS_TESTS_400: [[f64; 2]; 16] = [
    [85.07, 1.6],
    [91.65, 1.6],
    [98.73, 1.6],
    [106.36, 1.6],
    [114.57, 2.6],
    [123.42, 2.62],
    [133.25, 2.68],
    [144.29, 2.75],
    [156.79, 2.83],
    [171.03, 2.92],
    [187.41, 3.02],
    [206.41, 3.13],
    [228.63, 3.25],
    [254.93, 3.4],
    [286.52, 3.58],
    [324.91, 3.77],
];
const LABS_TESTS_400_FAN: [[f64; 2]; 11] = [
    [0.0, 0.0],
    [8.31, 3.1],
    [21.1, 3.5],
    [34.84, 3.8],
    [49.84, 3.5],
    [106.53, 4.53],
    [235.47, 4.53],
    [347.15, 3.53],
    [463.39, 2.53],
    [584.73, 2.53],
    [713.24, 0.03],
];
const LABS_TESTS_400_MOD: [[f64; 2]; 20] = [
    [2.62, 4.02],
    [3.49, 4.03],
    [4.66, 4.02],
    [6.22, 4.03],
    [8.31, 4.03],
    [11.1, 4.03],
    [14.84, 4.03],
    [19.84, 4.03],
    [26.53, 4.03],
    [35.47, 4.03],
    [47.42, 4.03],
    [63.39, 4.03],
    [84.73, 4.03],
    [113.24, 4.03],
    [151.33, 4.03],
    [202.2, 4.03],
    [270.15, 6.03],
    [470.15, 6.03],
    [570.15, 6.03],
    [670.15, 6.03],
];

const TEMP_CHARGE_CUT_DELTA: [f64; 12] = [
    -1.2, -0.6, 0.0, 0.6, 1.2, 1.2, 1.2, 1.2, 0.6, 0.0, -0.6, -1.2,
];

/// Type to represent electric storage heaters
pub struct ElecStorageHeater<T>
where
    T: InternalAirTempAccess,
{
    rated_power_in_w: f64,
    rated_power_instant_in_w: f64,
    air_flow_type: ElectricStorageHeaterAirFlowType,
    temp_dis_safe: f64,
    thermal_mass: f64,
    frac_convective: f64,
    u_ins: f64,
    temp_charge_cut: f64,
    n_units: usize,
    zone: Arc<T>,
    control: Arc<SetpointTimeControl>,
    charge_control: Arc<ToUChargeControl>,
    mass: f64,
    c_pcore: f64,
    t_core_target: f64,
    a: f64,
    c: f64,
    n: f64,
    thermal_mass_wall: f64,
    fan_power_in_w: f64,
    temp_air: f64,
    labs_tests: Vec<[f64; 2]>,
    t_core: f64,
    t_wall: f64,
    damper_fraction: f64,
    energy_in: f64,
}

impl<T> ElecStorageHeater<T>
where
    T: InternalAirTempAccess,
{
    pub fn new(
        rated_power_in_kw: f64,
        rated_power_instant_in_kw: f64,
        air_flow_type: ElectricStorageHeaterAirFlowType,
        temp_dis_safe: f64,
        thermal_mass: f64,
        frac_convective: f64,
        u_ins: f64,
        temp_charge_cut: f64,
        mass_core: f64,
        c_pcore: f64,
        temp_core_target: f64,
        a_core: f64,
        c_wall: f64,
        n_wall: f64,
        thermal_mass_wall: f64,
        fan_power_in_w: f64,
        n_units: usize,
        zone: Arc<T>,
        control: Arc<SetpointTimeControl>,
        charge_control: Arc<ToUChargeControl>,
    ) -> Self
    where
        T: InternalAirTempAccess,
    {
        let temp_air = zone.temp_internal_air();
        Self {
            rated_power_in_w: rated_power_in_kw * WATTS_PER_KILOWATT as f64,
            rated_power_instant_in_w: rated_power_instant_in_kw * WATTS_PER_KILOWATT as f64,
            air_flow_type,
            temp_dis_safe,
            thermal_mass,
            frac_convective,
            u_ins,
            temp_charge_cut,
            n_units,
            zone,
            control,
            charge_control,
            mass: mass_core, // kg of core
            c_pcore,
            t_core_target: temp_core_target,
            a: a_core, // m^2 transfer area between core and case or wall
            c: c_wall,
            n: n_wall,
            thermal_mass_wall,
            fan_power_in_w,
            temp_air,
            labs_tests: match air_flow_type {
                ElectricStorageHeaterAirFlowType::FanAssisted => LABS_TESTS_400_FAN.to_vec(),
                ElectricStorageHeaterAirFlowType::DamperOnly => LABS_TESTS_400.to_vec(),
            },
            t_core: 200.0,
            t_wall: 50.0,
            damper_fraction: 1.0,
            energy_in: 0.0,
        }
    }

    pub fn temp_setpnt(&self, simulation_time_iteration: &SimulationTimeIteration) -> Option<f64> {
        self.control.setpnt(simulation_time_iteration)
    }

    pub fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool> {
        self.control.in_required_period(simulation_time_iteration)
    }

    pub fn frac_convective(&self) -> f64 {
        self.frac_convective
    }

    /// Converts power value supplied to the correct energy unit
    /// Arguments:
    /// * `power` - Power value in watts
    /// * `timestep` - length of the timestep
    ///
    /// returns -- Energy in kWh
    fn convert_to_kwh(&self, power: f64, timestep: f64) -> f64 {
        return (power / WATTS_PER_KILOWATT as f64) * timestep * self.n_units as f64;
    }

    /// Correct nominal/json temp_charge_cut with monthly table
    fn temp_charge_cut_corr(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<f64> {
        simulation_time_iteration
            .current_month()
            .and_then(|month| Some(self.temp_charge_cut + TEMP_CHARGE_CUT_DELTA[month as usize]))
    }

    /// Calculates power required for unit
    /// Arguments
    /// * `time` - current time period that we are looking at
    /// * `t_core` - current temperature of the core
    ///
    /// returns -- Power required in watts
    fn electric_charge(
        &self,
        _time: f64,
        t_core: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        let temp_charge_cut_corr = self.temp_charge_cut_corr(simulation_time_iteration);

        if temp_charge_cut_corr.is_some() && self.temp_air >= temp_charge_cut_corr.unwrap() {
            return 0.0;
        }

        let target_charge = self.charge_control.target_charge(simulation_time_iteration);

        if self.charge_control.is_on(simulation_time_iteration)
            && t_core <= self.t_core_target * target_charge
        {
            let power_required =
                (self.t_core_target * target_charge - t_core) * self.mass * self.c_pcore
                    / simulation_time_iteration.timestep;
            if power_required > self.rated_power_in_w {
                self.rated_power_in_w
            } else {
                power_required
            }
        } else {
            0.0
        }
    }

    fn lab_test_ha(&self, t_core_rm_diff: f64) -> f64 {
        let x: Vec<f64> = self.labs_tests.iter().map(|row| row[0]).collect();
        let y: Vec<f64> = self.labs_tests.iter().map(|row| row[1]).collect();

        interp(&x, &y, t_core_rm_diff)
    }

    fn calculate_q_dis(&self, time: f64, t_core: f64, q_out_wall: f64, q_dis_modo: f64) -> f64 {
        if q_dis_modo == f64::INFINITY {
            self.lab_test_ha(t_core - self.temp_air) * (t_core - self.temp_air)
        } else if q_dis_modo == 0. {
            q_dis_modo
        } else {
            q_dis_modo - q_out_wall
        }
    }

    fn return_q_released(
        &mut self,
        new_temp_core_and_wall: [f64; 2],
        q_released: f64,
        q_dis: f64,
        q_in: f64,
        timestep: f64,
        q_instant: Option<f64>,
    ) -> f64 {
        let q_instant = q_instant.unwrap_or(0.0);

        // setting core and wall temperatures to new values, for next iteration
        self.t_core = new_temp_core_and_wall[0];
        self.t_wall = new_temp_core_and_wall[1];

        // the purpose of this calculation is to calculate fan energy required by the device
        let (energy_for_fan_kwh, power_for_fan) = if matches!(
            self.air_flow_type,
            ElectricStorageHeaterAirFlowType::FanAssisted
        ) && q_dis > 0.
        {
            (
                self.fan_power_in_w,
                self.convert_to_kwh(self.fan_power_in_w, timestep),
            )
        } else {
            (0., 0.)
        };

        // convert values to correct kwh unit
        let _q_in_kwh = self.convert_to_kwh(q_in, timestep);
        let _q_instant_kwh = self.convert_to_kwh(q_instant, timestep);

        // TODO save demand energy to energy supply

        self.convert_to_kwh(q_released + q_instant, timestep)
    }

    fn heat_balance(
        &mut self,
        temp_core_and_wall: [f64; 2],
        time: f64,
        q_dis_modo: Option<f64>,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> ([f64; 2], f64, f64, f64) {
        let q_dis_modo = q_dis_modo.unwrap_or(0.);

        let t_core = temp_core_and_wall[0];
        let t_wall = temp_core_and_wall[1];
        let q_in = self.electric_charge(time, t_core, simulation_time_iteration);
        self.energy_in = q_in;

        let q_out_wall = if t_wall >= self.temp_air {
            self.c * (t_wall - self.temp_air).powf(self.n)
        } else {
            self.c * (t_wall - self.temp_air)
        };

        // equation for calculating q_dis
        let q_dis = self.calculate_q_dis(time, t_core, q_out_wall, q_dis_modo);

        // Calculation of the U value between core and wall/case as
        // U value for the insulation and resistance of the air layer between the insulation and the wall/case
        let insulation = if q_dis > 0. {
            1. / (1. / self.u_ins + R_AIR_ON)
        } else {
            1. / (1. / self.u_ins + R_AIR_OFF)
        };

        // Equation for the heat transfer between the core and the wall/case of the heater
        let q_out_ins = insulation * self.a * (t_core - t_wall);

        // Variation of Core temperature as per heat balance inside the heater
        let d_t_core = (1. / (self.mass * self.c_pcore)) * (q_in - q_out_ins - q_dis);

        // Variation of Wall/case temperature as per heat balance in surface of the heater
        let d_t_wall = (1. / self.thermal_mass_wall) * (q_out_ins - q_out_wall);

        let q_released = q_dis + q_out_wall;

        ([d_t_core, d_t_wall], q_released, q_dis, q_in)
    }

    fn func_core_temperature_change_rate<'a, 'b>(
        &'a mut self,
        q_dis_modo: f64,
        simulation_time_iteration: &'a SimulationTimeIteration,
    ) -> impl FnMut(f64, &'a Vector<f64>) -> Result<Vector<f64>, anyhow::Error> {
        let simulation_time_iteration = (*simulation_time_iteration).clone();
        move |time: f64, t_core_and_wall: &'a Vector<f64>, _p: &mut ()| {
            let [d_t_core, d_t_wall] = self
                .heat_balance(
                    t_core_and_wall.try_into().unwrap(),
                    time,
                    Some(q_dis_modo),
                    &simulation_time_iteration,
                )
                .0;
            Ok(vector![d_t_core, d_t_wall])
        }
    }

    fn calculate_sol_and_q_released(
        &mut self,
        time_range: [f64; 2],
        temp_core_and_wall: [f64; 2],
        q_dis_modo: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Result<([f64; 2], f64, f64, f64), String> {
        let ode = Euler::default();
        let problem = ImplicitInitialValueProblemBuilder::new(
            &ode,
            time_range[0],
            vector![temp_core_and_wall[0], temp_core_and_wall[1]],
        )
        .t_end(time_range[1])
        .callback(self.func_core_temperature_change_rate(q_dis_modo, simulation_time_iteration))
        .build();
        let step_size = time_range[1] - time_range[0];
        let solver: BDF<f64> = BDF::new(6, step_size);
        let (_x, y): (Vec<f64>, Vec<Vector<f64>>) = solver.solve(&problem).unwrap();

        let new_temp_core_and_wall: [f64; 2] =
            [*y[0].iter().last().unwrap(), *y[1].iter().last().unwrap()];

        // let bdf = BDF6::new()
        //     .with_start(time_range[0])?
        //     .with_end(time_range[1])?
        //     .with_initial_conditions(&temp_core_and_wall)?
        //     .build();
        // let path = bdf.solve_ivp(
        //     self.func_core_temperature_change_rate(q_dis_modo, simulation_time_iteration),
        //     &mut (),
        // )?;
        // println!("BDF path result: {path:?}");
        let (_, q_released, q_dis, q_in) = self.heat_balance(
            new_temp_core_and_wall,
            time_range[1],
            Some(q_dis_modo),
            simulation_time_iteration,
        );

        Ok((new_temp_core_and_wall, q_released, q_dis, q_in))
    }

    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> f64 {
        self.temp_air = self.zone.temp_internal_air();
        let timestep = simulation_time_iteration.timestep;
        let time_unit = 3_600. * timestep;
        let current_hour = simulation_time_iteration.current_hour();
        let time_range = [
            current_hour as f64 * time_unit,
            (current_hour + 1) as f64 * time_unit,
        ];
        let temp_core_and_wall = [self.t_core, self.t_wall];

        // Converting energy_demand from kWh to Wh and distributing it through all units
        let energy_demand = energy_demand * WATTS_PER_KILOWATT as f64 / self.n_units as f64;

        // #################################################
        // # Step 1                                        #
        // #################################################
        // first calculate how much the system is leaking without active discharging
        let (new_temp_core_and_wall, q_released, q_dis, q_in) = self
            .calculate_sol_and_q_released(
                time_range,
                temp_core_and_wall,
                0.,
                simulation_time_iteration,
            )
            .expect("expected ivp solver to succeed at Step 1");

        // if Q_released is more than what the zone wants, that's it
        if q_released >= energy_demand / timestep {
            return self.return_q_released(
                new_temp_core_and_wall,
                q_released,
                q_dis,
                q_in,
                timestep,
                None,
            );
        }

        // #################################################
        // # Step 2                                        #
        // #################################################
        // Zone needs more than leaked, let's calculate with max discharging capability
        let (new_temp_core_and_wall, q_released, q_dis, q_in) = self
            .calculate_sol_and_q_released(
                time_range,
                temp_core_and_wall,
                f64::INFINITY,
                simulation_time_iteration,
            )
            .expect("expected ivh solver to succeed at Step 2");

        // If Q_released is not sufficient for zone demand, that's it
        // unless there is instant backup that can top up the energy provided
        if q_released < energy_demand / timestep {
            let power_supplied_instant = if self.rated_power_instant_in_w > 0. {
                *[
                    energy_demand / timestep - q_released,
                    self.rated_power_instant_in_w,
                ]
                .iter()
                .max_by(|a, b| a.total_cmp(b).reverse())
                .unwrap()
            } else {
                0.0
            };

            // The system can only discharge the maximum amount, zone doesn't get everything it needs
            return self.return_q_released(
                new_temp_core_and_wall,
                q_released,
                q_dis,
                q_in,
                timestep,
                Some(power_supplied_instant),
            );
        }

        // #################################################
        // # Step 3                                        #
        // #################################################
        // Zone actually needs an amount of energy that can be released by the system:
        // Let's call the heat balance forcing that amount (assuming perfect damper or
        // fan assisted control of the unit)
        let q_dis = energy_demand / timestep;
        let (new_temp_core_and_wall, q_released, q_dis, q_in) = self
            .calculate_sol_and_q_released(
                time_range,
                temp_core_and_wall,
                q_dis,
                simulation_time_iteration,
            )
            .expect("expected ivh solver to work at Step 3");

        self.return_q_released(
            new_temp_core_and_wall,
            q_released,
            q_dis,
            q_in,
            timestep,
            None,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::controls::time_control::ToUChargeControl;
    use crate::simulation_time::SimulationTime;
    use rstest::*;

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 24., 1.)
    }

    struct MockZone();

    impl InternalAirTempAccess for MockZone {
        fn temp_internal_air(&self) -> f64 {
            20.0
        }
    }

    #[fixture]
    pub fn elec_storage_heater(simulation_time: SimulationTime) -> ElecStorageHeater<MockZone> {
        let zone = MockZone();
        let control = SetpointTimeControl::new(
            vec![
                Some(15.0),
                Some(15.0),
                Some(15.0),
                Some(15.0),
                Some(15.0),
                Some(15.0),
                Some(15.0),
                Some(21.0),
                Some(21.0),
                Some(21.0),
                Some(21.0),
                Some(21.0),
                Some(21.0),
                Some(21.0),
                Some(21.0),
                Some(21.0),
                Some(21.0),
                Some(21.0),
                Some(21.0),
                Some(21.0),
                Some(15.0),
                Some(15.0),
                Some(15.0),
                Some(15.0),
            ],
            0,
            1.,
            None,
            None,
            None,
            None,
            1.,
        )
        .unwrap();
        let charge_control = ToUChargeControl {
            schedule: vec![
                true, true, true, true, true, true, true, true, false, false, false, false, false,
                false, false, false, true, true, true, true, false, false, false, false,
            ],
            start_day: 0,
            time_series_step: 1.0,
            charge_level: vec![1.0, 0.8],
        };

        ElecStorageHeater::new(
            4.0,
            0.75,
            ElectricStorageHeaterAirFlowType::FanAssisted,
            60.0,
            0.01278,
            0.7,
            0.3,
            15.5,
            130.0,
            920.0,
            450.0,
            4.0,
            8.0,
            0.9,
            23.0,
            11.0,
            2,
            Arc::new(zone),
            Arc::new(control),
            Arc::new(charge_control),
        )
    }

    #[rstest]
    pub fn test_demand_energy(
        mut elec_storage_heater: ElecStorageHeater<MockZone>,
        simulation_time: SimulationTime,
    ) {
        let inputs = [
            4.69, 3.59, 4.26, 2.82, 0.31, 3.72, 2.11, 6.55, 7.59, 7.55, 4.52, 2.92, 3.42, 5.83,
            4.26, 3.63, 4.38, 5.34, 4.65, 3.85, 0.0, 1.86, 2.27, 2.62,
        ];
        let expected = [
            3.18, 2.93, 2.71, 2.49, 0.31, 2.26, 2.11, 2.03, 1.95, 1.9, 1.86, 1.81, 1.77, 1.73, 1.7,
            1.67, 1.65, 1.63, 1.61, 1.6, 0.02, 1.58, 1.57, 1.56,
        ];
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                round_by_precision(elec_storage_heater.demand_energy(inputs[t_idx], &t_it), 1e2),
                expected[t_idx]
            );
        }
    }

    fn round_by_precision(src: f64, precision: f64) -> f64 {
        (precision * src).round() / precision
    }
}
