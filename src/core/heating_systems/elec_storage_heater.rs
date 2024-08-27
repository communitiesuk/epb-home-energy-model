use crate::core::controls::time_control::{ControlBehaviour, ToUChargeControl};
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::space_heat_demand::zone::Zone;
use crate::core::units::WATTS_PER_KILOWATT;
use crate::input::AirFlowType;
use crate::simulation_time::SimulationTimeIteration;
use anyhow::anyhow;
use interp::interp;
use std::sync::Arc;

/// This module provides object(s) to model the behaviour of electric storage heaters.

// J/kg/K air specific heat
const C_P: f64 = 1.0054;

// The value of R (resistance of air gap) depends on several factors such as the thickness of the air gap, the
// temperature difference across the gap, and the air flow characteristics.
// Here are some typical values for R per inch of air gap thickness:
// Still air: 0.17 m²·K/W
// Air movement (average): 0.07 m²·K/W
// Air movement (high): 0.04 m²·K/W
// Note: These values are for a temperature difference of 24°F (14°C) and a pressure difference of 1 inch of water
// column. The values will change with changes in temperature and pressure differences.

// 0.17  Thermal resistance of the air layer between the insulation and the wall of the device.
//       Assuming no change of resistance with temp diff. (Rough approx)
const R_AIR_OFF: f64 = 0.17;
// Same as above when the damper is on
const R_AIR_ON: f64 = 0.17;

// labs_test for electric storage heater reaching 300 degC
// This represents the temperature difference between the core and the room on the first column
// and the fraction of air flow relating to the nominal as defined above on the second column
const LABS_TESTS_400: [(f64, f64); LABS_TESTS_400_LEN] = [
    (85.07, 1.6),
    (91.65, 1.6),
    (98.73, 1.6),
    (106.36, 1.6),
    (114.57, 2.6),
    (123.42, 2.62),
    (133.25, 2.68),
    (144.29, 2.75),
    (156.79, 2.83),
    (171.03, 2.92),
    (187.41, 3.02),
    (206.41, 3.13),
    (228.63, 3.25),
    (254.93, 3.4),
    (286.52, 3.58),
    (324.91, 3.77),
];
const LABS_TESTS_400_LEN: usize = 16;

const LABS_TESTS_400_X: [f64; LABS_TESTS_400_LEN] = {
    let mut arr = [0.0; 16];
    let mut i = 0usize;
    while i < LABS_TESTS_400_LEN {
        arr[i] = LABS_TESTS_400[i].0;
        i += 1;
    }
    arr
};

const LABS_TESTS_400_Y: [f64; LABS_TESTS_400_LEN] = {
    let mut arr = [0.0; LABS_TESTS_400_LEN];
    let mut i = 0usize;
    while i < LABS_TESTS_400_LEN {
        arr[i] = LABS_TESTS_400[i].1;
        i += 1;
    }
    arr
};

const LABS_TESTS_400_FAN: [(f64, f64); 11] = [
    (0.0, 0.0),
    (8.31, 3.1),
    (21.1, 3.5),
    (34.84, 3.8),
    (49.84, 3.5),
    (106.53, 4.53),
    (235.47, 4.53),
    (347.15, 3.53),
    (463.39, 2.53),
    (584.73, 2.53),
    (713.24, 0.03),
];
const LABS_TESTS_400_FAN_LEN: usize = 11;

const LABS_TESTS_400_FAN_X: [f64; LABS_TESTS_400_FAN_LEN] = {
    let mut arr = [0.0; LABS_TESTS_400_FAN_LEN];
    let mut i = 0usize;
    while i < LABS_TESTS_400_FAN_LEN {
        arr[i] = LABS_TESTS_400_FAN[i].0;
        i += 1;
    }
    arr
};

const LABS_TESTS_400_FAN_Y: [f64; LABS_TESTS_400_FAN_LEN] = {
    let mut arr = [0.0; LABS_TESTS_400_FAN_LEN];
    let mut i = 0usize;
    while i < LABS_TESTS_400_FAN_LEN {
        arr[i] = LABS_TESTS_400_FAN[i].1;
        i += 1;
    }
    arr
};

const LABS_TESTS_400_MOD: [(f64, f64); 20] = [
    (2.62, 4.02),
    (3.49, 4.03),
    (4.66, 4.02),
    (6.22, 4.03),
    (8.31, 4.03),
    (11.1, 4.03),
    (14.84, 4.03),
    (19.84, 4.03),
    (26.53, 4.03),
    (35.47, 4.03),
    (47.42, 4.03),
    (63.39, 4.03),
    (84.73, 4.03),
    (113.24, 4.03),
    (151.33, 4.03),
    (202.2, 4.03),
    (270.15, 6.03),
    (470.15, 6.03),
    (570.15, 6.03),
    (670.15, 6.03),
];

const TEMP_CHARGE_CUT_DELTA: [f64; 12] = [
    -1.2, -0.6, 0.0, 0.6, 1.2, 1.2, 1.2, 1.2, 0.6, 0.0, -0.6, -1.2,
];

/// Struct to represent electric storage heaters
pub(crate) struct ElecStorageHeater {
    pwr_in: f64,
    pwr_instant: f64,
    air_flow_type: AirFlowType,
    t_dis_safe: f64,
    thermal_mass: f64,
    frac_convective: f64,
    u_ins: f64,
    temp_charge_cut: f64,
    n_units: usize,
    zone: Arc<Zone>,
    energy_supply_connection: EnergySupplyConnection,
    control: Arc<Box<dyn ControlBehaviour>>,
    charge_control: Arc<ToUChargeControl>,
    mass: f64,
    c_pcore: f64,
    t_core_target: f64,
    a: f64,
    c: f64,
    n: f64,
    thermal_mass_wall: f64,
    fan_pwr: f64,
    temp_air: f64,
    damper_fraction: f64,
    t_core: f64,
    t_wall: f64,
    energy_in: f64,
}

impl ElecStorageHeater {
    /// Arguments
    ///
    /// * `rated_power`          - in kW (Charging)
    /// * `rated_power_instant`  - in kW (Instant backup)
    /// * `air_flow_type`        - air flow type of Electric Storage Heater:
    ///                              - fan-assisted
    ///                              - damper-only
    /// * `temp_dis_safe`        - safe temperature to discharge hot air from device (60 degC)
    /// * `thermal_mass`         - thermal mass of emitters, in kWh / K
    /// * `frac_convective`      - convective fraction for heating (TODO: Check if necessary)
    /// * `u_ins`                - U-value insulation between core and case [W/m^2/K]
    /// * `temp_charge_cut`      - Room temperature at which, if sensed during a charging hour, the heater won't charge
    /// * `mass_core`            - kg mass core material [kg]
    /// * `c_pcore`              - thermal capacity of core material [J/kg/K]
    /// * `temp_core_target`     - target temperature for the core material on charging mode
    ///                              - this might include weather compensation with future more
    ///                              - advances controls
    /// * `a_core`               - Transfer area between the core and air [m2]
    /// * `c_wall`               - constant from characteristic equation of emitters (e.g. derived from BS EN 442 tests)
    /// * `n_wall`               - exponent from characteristic equation of emitters (e.g. derived from BS EN 442 tests)
    /// * `thermal_mass_wall`    - thermal mass of the case
    /// * `fan_pwr`              - Fan power [W]
    /// * `n_units`              - number of units install in zone
    /// * `zone`                 - zone where the unit(s) is/are installed
    /// * `energy_supply_connection`   - reference to EnergySupplyConnection object
    /// * `control`              - reference to a control object which must implement is_on() and setpnt() funcs
    /// * `charge_control`
    pub(crate) fn new(
        rated_power: f64,
        rated_power_instant: f64,
        air_flow_type: AirFlowType,
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
        fan_pwr: f64,
        n_units: usize,
        zone: Arc<Zone>,
        energy_supply_connection: EnergySupplyConnection,
        control: Arc<Box<dyn ControlBehaviour>>,
        charge_control: Arc<ToUChargeControl>,
    ) -> Self {
        let temp_air = zone.temp_internal_air();
        Self {
            pwr_in: rated_power * WATTS_PER_KILOWATT as f64,
            pwr_instant: rated_power_instant * WATTS_PER_KILOWATT as f64,
            air_flow_type,
            t_dis_safe: temp_dis_safe,
            thermal_mass,
            frac_convective,
            u_ins,
            temp_charge_cut,
            n_units,
            zone,
            energy_supply_connection,
            control,
            charge_control,
            mass: mass_core,
            c_pcore,
            t_core_target: temp_core_target,
            a: a_core,
            c: c_wall,
            n: n_wall,
            thermal_mass_wall,
            fan_pwr,
            temp_air,
            damper_fraction: 1.0,
            // initial conditions
            t_core: 200.0,
            t_wall: 50.0,
            energy_in: Default::default(),
        }
    }

    fn temp_setpnt(&self, simtime: SimulationTimeIteration) -> Option<f64> {
        self.control.setpnt(&simtime)
    }

    fn in_required_period(&self, simtime: SimulationTimeIteration) -> Option<bool> {
        self.control.in_required_period(&simtime)
    }

    fn frac_convective(&self) -> f64 {
        self.frac_convective
    }

    /// Converts power value supplied to the correct energy unit
    ///
    /// Arguments
    ///
    /// * `power` - Power value in watts
    /// * `timestep` - length of the timestep
    ///
    /// returns -- Energy in kWh
    fn convert_to_kwh(&self, power: f64, timestep: f64) -> f64 {
        power / WATTS_PER_KILOWATT as f64 * timestep * self.n_units as f64
    }

    /// Correct nominal/json temp_charge_cut with monthly table
    ///
    /// returns -- temp_charge_cut (corrected)
    fn temp_charge_cut_corr(&self, simtime: SimulationTimeIteration) -> anyhow::Result<f64> {
        let current_month = simtime.current_month().ok_or_else(|| anyhow!("Could not get the current month in temp_charge_cut_corr within elec_storage_heater module."))?;
        Ok(self.temp_charge_cut + TEMP_CHARGE_CUT_DELTA[current_month as usize])
    }

    /// Calculates power required for unit
    ///
    /// Arguments
    ///
    /// * `time` - current time period that we are looking at
    /// * `t_core` - current temperature of the core
    ///
    /// returns -- Power required in watts
    fn electric_charge(
        &self,
        _time: f64,
        t_core: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        if self.temp_air >= self.temp_charge_cut_corr(simtime)? {
            return Ok(0.0);
        }

        let target_charge = self.charge_control.target_charge(&simtime);
        Ok(
            if self.charge_control.is_on(&simtime) && t_core <= self.t_core_target * target_charge {
                let pwr_required =
                    (self.t_core_target * target_charge - t_core) * self.mass * self.c_pcore
                        / simtime.timestep;
                if pwr_required > self.pwr_in {
                    self.pwr_in
                } else {
                    pwr_required
                }
            } else {
                0.0
            },
        )
    }

    /// accessor method to the relevant labs_tests values, rather than storing these against self
    fn labs_tests(&self) -> (&[f64], &[f64]) {
        match self.air_flow_type {
            AirFlowType::FanAssisted => (&LABS_TESTS_400_FAN_X, &LABS_TESTS_400_FAN_Y),
            AirFlowType::DamperOnly => (&LABS_TESTS_400_X, &LABS_TESTS_400_Y),
        }
    }

    fn lab_test_ha(&self, t_core_rm_diff: f64) -> f64 {
        let (x, y) = self.labs_tests();
        interp(x, y, t_core_rm_diff)
    }

    fn calculate_q_dis(
        &self,
        _time: f64,
        t_core: f64,
        q_out_wall: f64,
        q_dis_modo: QDisModoArg,
    ) -> f64 {
        match q_dis_modo {
            QDisModoArg::Max => self.lab_test_ha(t_core - self.temp_air) * (t_core - self.temp_air),
            QDisModoArg::Float(0.0) => 0.0,
            QDisModoArg::Float(q_dis_modo) => q_dis_modo - q_out_wall,
        }
    }

    // TODO: complete implementation and tests
}

#[derive(Clone, Copy, Debug)]
enum QDisModoArg {
    Max,
    Float(f64),
}
