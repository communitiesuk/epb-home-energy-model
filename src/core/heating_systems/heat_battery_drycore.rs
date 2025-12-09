//! This module provides the base behaviour for dry core heat storage systems.
//! This includes the common functionality for electrical storage and discharge
//! that is shared between Electric Storage Heaters and Dry Core Heat Batteries.

use crate::core::controls::time_control::Control;
use crate::core::units::HOURS_PER_DAY;
use crate::hem_core::simulation_time::SimulationTimeIteration;
use crate::input::ControlLogicType;
use crate::statistics::{np_interp, np_interp_with_extrapolate};
use anyhow::bail;
use atomic_float::AtomicF64;
use itertools::Itertools;
use nalgebra::{Vector1, Vector3};
use ode_solvers::dop_shared::OutputType;
use ode_solvers::{Dopri5, System};
use std::sync::atomic::Ordering;
use std::sync::{Arc, Weak};

type State = Vector1<f64>;
type Time = f64;

type EnergyOutputState = Vector3<f64>;

pub(crate) enum OutputMode {
    Min,
    Max,
}

pub(crate) struct HeatStorageDryCore {
    pwr_in: f64,
    storage_capacity: f64,
    n_units: u32,
    charge_control: Arc<Control>,
    dry_core_min_output: Vec<[f64; 2]>,
    dry_core_max_output: Vec<[f64; 2]>,
    state_of_charge: AtomicF64,
    energy_in: AtomicF64,
    demand_met: AtomicF64,
    demand_unmet: AtomicF64,
    soc_max_array: Vec<f64>,
    power_max_array: Vec<f64>,
    soc_min_array: Vec<f64>,
    power_min_array: Vec<f64>,
    heat_retention_ratio: f64,
    // weak reference back to value that is composing this struct, as we need two-way references
    owner: Weak<dyn HeatBatteryDryCoreCommonBehaviour>,
}

impl HeatStorageDryCore {
    /// * `pwr_in`               -- in kW (Charging)
    /// * `storage_capacity`     -- in kWh
    /// * `n_units`              -- number of units installed
    /// * `charge_control`       -- reference to a ChargeControl object
    /// * `dry_core_min_output`       -- Data from test showing the output from the storage heater when not actively
    ///                                outputting heat, i.e. case losses only (with units kW)
    /// * `dry_core_max_output`       -- Data from test showing the output from the storage heater when it is actively
    ///                                outputting heat, e.g. damper open / fan running (with units kW)
    /// * `state_of_charge_init` -- state of charge at initialisation of dry core heat storage
    pub(crate) fn new(
        pwr_in: f64,
        storage_capacity: f64,
        n_units: u32,
        charge_control: Arc<Control>,
        dry_core_min_output: Vec<[f64; 2]>,
        dry_core_max_output: Vec<[f64; 2]>,
        state_of_charge_init: f64,
        owner: Weak<dyn HeatBatteryDryCoreCommonBehaviour>,
    ) -> anyhow::Result<Self> {
        // Convert dry_core_max_output to vecs without sorting
        let soc_max_array = dry_core_max_output.iter().map(|f| f[0]).collect_vec();
        let power_max_array = dry_core_max_output.iter().map(|f| f[1]).collect_vec();

        // Convert dry_core_min_output to vecs without sorting
        let soc_min_array = dry_core_min_output.iter().map(|f| f[0]).collect_vec();
        let power_min_array = dry_core_min_output.iter().map(|f| f[1]).collect_vec();

        // Validate that both SOC arrays are in strictly increasing order
        if !soc_max_array
            .iter()
            .zip(soc_max_array.iter().skip(1))
            .all(|(a, b)| a <= b)
        {
            bail!("dry_core_max_output SOC values must be in increasing order (from 0.0 to 1.0).");
        }

        if !soc_min_array
            .iter()
            .zip(soc_min_array.iter().skip(1))
            .all(|(a, b)| a <= b)
        {
            bail!("dry_core_min_output SOC values must be in increasing order (from 0.0 to 1.0).");
        }

        // Validate that both SOC arrays start at 0.0 and end at 1.0
        if !is_close!(*soc_max_array.first().unwrap(), 0.) {
            bail!("The first SOC value in dry_core_max_output must be 0.0 (fully discharged).");
        }

        if !is_close!(*soc_max_array.last().unwrap(), 1.) {
            bail!("The last SOC value in dry_core_max_output must be 1.0 (fully charged).");
        }

        if !is_close!(*soc_min_array.first().unwrap(), 0.) {
            bail!("The first SOC value in dry_core_min_output must be 0.0 (fully discharged).");
        }

        if !is_close!(*soc_min_array.last().unwrap(), 1.) {
            bail!("The last SOC value in dry_core_min_output must be 1.0 (fully charged).");
        }

        // Validate that for any SOC, power_max >= power_min
        // Sample a fine grid of SOCs and ensure power_max >= power_min
        let fine_soc: Vec<f64> =
            crate::core::heating_systems::elec_storage_heater::linspace(0., 1., 100);

        let power_max_fine: Vec<f64> = fine_soc
            .iter()
            .map(|s| np_interp(*s, &soc_max_array, &power_max_array))
            .collect();

        // TODO in Python a fill_value is used to make this return 0 when out of bounds
        let power_min_fine: Vec<f64> = fine_soc
            .iter()
            .map(|s| np_interp(*s, &soc_min_array, &power_min_array))
            .collect();

        for i in 0..fine_soc.len() {
            if power_max_fine[i] < power_min_fine[i] {
                bail!("At all SOC levels, dry_core_max_output must be >= dry_core_min_output.")
            }
        }

        let heat_retention_ratio =
            Self::heat_retention_output(&soc_min_array, &power_min_array, storage_capacity);

        Ok(Self {
            pwr_in,
            storage_capacity,
            n_units,
            charge_control,
            dry_core_min_output: dry_core_min_output.clone(),
            dry_core_max_output: dry_core_max_output.clone(),
            state_of_charge: AtomicF64::new(state_of_charge_init),
            energy_in: AtomicF64::new(0.0),
            demand_met: AtomicF64::new(0.0),
            demand_unmet: AtomicF64::new(0.0),
            soc_max_array,
            power_max_array,
            soc_min_array,
            power_min_array,
            heat_retention_ratio,
            owner,
        })
    }

    pub(super) fn heat_retention_output(
        soc_array: &[f64],
        power_array: &[f64],
        storage_capacity: f64,
    ) -> f64 {
        // Simulates the heat retention over 16 hours in OutputMode.MIN.

        // Starts with a SOC of 1.0 and calculates the SOC after 16 hours.
        // This is a self-contained function, and the SOC is not stored in self.__state_of_charge.

        // :return: Final SOC after 16 hours.

        // Set initial state of charge to 1.0 (fully charged)
        let initial_soc = 1.0;

        // Total time for the simulation (16 hours)
        let total_time = 16.0; // This is the value from BS EN 60531 for determining heat retention ability

        // Select the SOC and power arrays for OutputMode.MIN
        let soc_ode = SocOdeFunction {
            soc_array,
            power_array,
            storage_capacity,
        };

        let f = soc_ode; // f - Structure implementing the System trait
        let x = 0.; // x - Initial value of the independent variable (usually time)
        let x_end = total_time; // x_end - Final value of the independent variable
        let dx = 0.; // dx - Increment in the dense output. This argument has no effect if the output type is Sparse
        let y0: State = State::new(initial_soc); // y - Initial value of the dependent variable(s)

        // scipy implementation for reference:
        // https://github.com/scipy/scipy/blob/6b657ede0c3c4cffef3156229afddf02a2b1d99a/scipy/integrate/_ivp/rk.py#L293
        let rtol = 1e-3; // rtol - set from scipy docs - Relative tolerance used in the computation of the adaptive step size
        let atol = 1e-6; // atol - set from scipy docs - Absolute tolerance used in the computation of the adaptive step size
        let h = 0.; // initial step size - 0
        let safety_factor = 0.9; // matches scipy implementation
        let beta = 0.; // setting this to 0 gives us an alpha of 0.2 and matches scipy's adaptive step size logic (default was 0.04)
        let fac_min = 0.2; // matches scipy implementation
        let fac_max = 10.; // matches scipy implementation
        let h_max = x_end - x;
        let n_max = 100000;
        let n_stiff = 1000;
        let mut stepper = Dopri5::from_param(
            f,
            x,
            x_end,
            dx,
            y0,
            rtol,
            atol,
            safety_factor,
            beta,
            fac_min,
            fac_max,
            h_max,
            h,
            n_max,
            n_stiff,
            OutputType::Sparse,
        );

        // Solve the ODE for SOC and cumulative energy delivered
        let _ = stepper.integrate();

        // sol = solve_ivp(soc_ode, [0, total_time], [initial_soc], method='RK45', rtol=1e-1, atol=1e-3)

        // Final state of charge after 16 hours
        let final_soc = stepper.y_out().last().unwrap()[0];

        // Clip the final SOC to ensure it's between 0 and 1

        // Return the final state of charge after 16 hours
        clip(final_soc, 0., 1.)
    }

    pub(crate) fn energy_output(
        &self,
        mode: OutputMode,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64, f64, f64)> {
        let (soc_array, power_array) = match mode {
            OutputMode::Min => (&self.soc_min_array, &self.power_min_array),
            OutputMode::Max => (&self.soc_max_array, &self.power_max_array),
        };

        let target_charge = self.target_electric_charge(*simulation_time_iteration)?;

        let (charge_rate, soc_max) = if target_charge > 0. {
            (self.pwr_in, target_charge)
        } else {
            (0., 1.)
        };

        // TODO Stop function!!
        let energy_output_soc_ode = EnergyOutputSocOdeFunction {
            soc_array,
            power_array,
            storage_capacity: self.storage_capacity,
            soc_max,
            charge_rate,
            pwr_in: self.pwr_in,
            target_charge,
        };

        // Set initial conditions
        let current_soc = self.state_of_charge.load(Ordering::SeqCst);
        let initial_energy_charged = 0.; // No energy charged initially
        let initial_energy_delivered = 0.; // No energy delivered initially
        let time_remaining = simulation_time_iteration.timestep; // in hours

        let f = energy_output_soc_ode; // f - Structure implementing the System trait
        let x = 0.; // x - Initial value of the independent variable (usually time)
        let x_end = time_remaining; // x_end - Final value of the independent variable
        let dx = 0.; // dx - Increment in the dense output. This argument has no effect if the output type is Sparse
        let y0: EnergyOutputState = EnergyOutputState::new(
            current_soc,
            initial_energy_charged,
            initial_energy_delivered,
        ); // y - Initial value of the dependent variable(s)

        // scipy implementation for reference:
        // https://github.com/scipy/scipy/blob/6b657ede0c3c4cffef3156229afddf02a2b1d99a/scipy/integrate/_ivp/rk.py#L293
        let rtol = 1e-4; // rtol - set to match Python - Relative tolerance used in the computation of the adaptive step size
        let atol = 1e-6; // atol - set from scipy docs - Absolute tolerance used in the computation of the adaptive step size
        let h = 0.; // initial step size - 0
        let safety_factor = 0.9; // matches scipy implementation
        let beta = 0.; // setting this to 0 gives us an alpha of 0.2 and matches scipy's adaptive step size logic (default was 0.04)
        let fac_min = 0.2; // matches scipy implementation
        let fac_max = 10.; // matches scipy implementation
        let h_max = x_end - x;
        let n_max = 100000;
        let n_stiff = 1000;
        let mut stepper = Dopri5::from_param(
            f,
            x,
            x_end,
            dx,
            y0,
            rtol,
            atol,
            safety_factor,
            beta,
            fac_min,
            fac_max,
            h_max,
            h,
            n_max,
            n_stiff,
            OutputType::Sparse,
        );

        // Solve the ODE for SOC and cumulative energy delivered
        let _ = stepper.integrate();

        // Final state of charge after 16 hours
        let final_soc = stepper.y_out().last().unwrap()[0];

        // Total energy charged during the timestep
        let total_energy_charged = stepper.y_out().last().unwrap()[1];

        // Total energy delivered during the timestep
        let total_energy_delivered = stepper.y_out().last().unwrap()[2];
        // Total time used in delivering energy
        let time_used = stepper.x_out().last().unwrap(); // TODO implement with root finder

        // Return the total energy delivered, time used, and total energy charged
        Ok((
            total_energy_delivered,
            *time_used,
            total_energy_charged,
            final_soc,
        ))
    }

    fn energy_output_with_losses(
        &self,
        mode: OutputMode,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64, f64, f64, f64)> {
        // TODO: complete porting this function!

        let (soc_array, power_array) = match mode {
            OutputMode::Min => (&self.soc_min_array, &self.power_min_array),
            OutputMode::Max => (&self.soc_max_array, &self.power_max_array),
        };

        let target_charge = self.target_electric_charge(*simulation_time_iteration)?;

        let (charge_rate, soc_max) = if target_charge > 0. {
            (self.pwr_in, target_charge)
        } else {
            (0., 1.)
        };

        // TODO Stop function!!
        let energy_output_soc_ode = EnergyOutputSocOdeFunction {
            soc_array,
            power_array,
            storage_capacity: self.storage_capacity,
            soc_max,
            charge_rate,
            pwr_in: self.pwr_in,
            target_charge,
        };

        // Set initial conditions
        let current_soc = self.state_of_charge.load(Ordering::SeqCst);
        let initial_energy_charged = 0.; // No energy charged initially
        let initial_energy_delivered = 0.; // No energy delivered initially
        let time_remaining = simulation_time_iteration.timestep; // in hours

        let f = energy_output_soc_ode; // f - Structure implementing the System trait
        let x = 0.; // x - Initial value of the independent variable (usually time)
        let x_end = time_remaining; // x_end - Final value of the independent variable
        let dx = 0.; // dx - Increment in the dense output. This argument has no effect if the output type is Sparse
        let y0: EnergyOutputState = EnergyOutputState::new(
            current_soc,
            initial_energy_charged,
            initial_energy_delivered,
        ); // y - Initial value of the dependent variable(s)

        // scipy implementation for reference:
        // https://github.com/scipy/scipy/blob/6b657ede0c3c4cffef3156229afddf02a2b1d99a/scipy/integrate/_ivp/rk.py#L293
        let rtol = 1e-4; // rtol - set to match Python - Relative tolerance used in the computation of the adaptive step size
        let atol = 1e-6; // atol - set from scipy docs - Absolute tolerance used in the computation of the adaptive step size
        let h = 0.; // initial step size - 0
        let safety_factor = 0.9; // matches scipy implementation
        let beta = 0.; // setting this to 0 gives us an alpha of 0.2 and matches scipy's adaptive step size logic (default was 0.04)
        let fac_min = 0.2; // matches scipy implementation
        let fac_max = 10.; // matches scipy implementation
        let h_max = x_end - x;
        let n_max = 100000;
        let n_stiff = 1000;
        let mut stepper = Dopri5::from_param(
            f,
            x,
            x_end,
            dx,
            y0,
            rtol,
            atol,
            safety_factor,
            beta,
            fac_min,
            fac_max,
            h_max,
            h,
            n_max,
            n_stiff,
            OutputType::Sparse,
        );

        // Solve the ODE for SOC and cumulative energy delivered
        let _ = stepper.integrate();

        // Final state of charge after 16 hours
        let final_soc = stepper.y_out().last().unwrap()[0];

        // Total energy charged during the timestep
        let total_energy_charged = stepper.y_out().last().unwrap()[1];

        // Total energy delivered during the timestep
        let total_energy_delivered = stepper.y_out().last().unwrap()[2];
        // Total time used in delivering energy
        let time_used = stepper.x_out().last().unwrap(); // TODO implement with root finder

        // Return the total energy delivered, time used, and total energy charged
        Ok((
            total_energy_delivered,
            *time_used,
            total_energy_charged,
            final_soc,
            Default::default(),
        ))
    }

    pub(crate) fn target_electric_charge(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // Calculates target charge from potential to charge system

        let charge_control = match self.charge_control.as_ref() {
            Control::Charge(charge_control) => charge_control,
            _ => unreachable!("charge_control must be of type ChargeControl"),
        };

        let logic_type = charge_control.logic_type();

        let temp_air = self.owner().get_temp_for_charge_control();

        match logic_type {
            ControlLogicType::Manual => {
                // Implements the "Manual" control logic for ESH
                charge_control.target_charge(simulation_time_iteration, None)
            }
            ControlLogicType::Automatic => {
                // Implements the "Automatic" control logic for ESH
                // Automatic charge control can be achieved using internal thermostat(s) to
                // control the extent of charging of the heaters.
                // Availability of electricity to the heaters may be controlled by the electricity
                // supplier on the basis of daily weather predictions (see 24-hour tariff, 12.4.3);
                // this should be treated as automatic charge control.
                // This is currently included by the schedule parameter in charge control object

                // TODO (from Python): Check and implement if external temperature sensors are also used for Automatic controls.
                charge_control.target_charge(simulation_time_iteration, temp_air)
            }
            ControlLogicType::Celect => {
                // Implements the "CELECT" control logic for ESH
                // A CELECT-type controller has electronic sensors throughout the dwelling linked
                // to a central control device. It monitors the individual room sensors and optimises
                // the charging of all the storage heaters individually (and may select direct acting
                // heaters in preference to storage heaters).
                charge_control.target_charge(simulation_time_iteration, temp_air)
            }
            ControlLogicType::Hhrsh => {
                // Implements the "HHRSH" control logic for ESH
                // A ‘high heat retention storage heater’ is one with heat retention not less
                // than 45% measured according to BS EN 60531. It incorporates a timer, electronic
                // room thermostat and fan to control the heat output. It is also able to estimate
                // the next day’s heating demand based on external temperature, room temperature
                // settings and heat demand periods.

                let mut energy_to_store = charge_control.energy_to_store(
                    self.demand_met.load(Ordering::SeqCst)
                        + self.demand_unmet.load(Ordering::SeqCst),
                    self.owner().get_zone_setpoint(),
                    simulation_time_iteration,
                );

                // None means not enough past data to do the calculation (Initial 24h of the calculation)
                // We go for a full load of the hhrsh
                if energy_to_store.is_nan() {
                    energy_to_store = self.pwr_in * HOURS_PER_DAY as f64;
                }

                let target_charge_hhrsh = if energy_to_store > 0. {
                    let current_state_of_charge = self.state_of_charge.load(Ordering::SeqCst);
                    let energy_stored = current_state_of_charge * self.storage_capacity; // kWh

                    let energy_to_add = if self.heat_retention_ratio <= 0. {
                        self.storage_capacity - energy_stored // kWh
                    } else {
                        (1.0 / self.heat_retention_ratio) * (energy_to_store - energy_stored)
                        // kWh
                    };
                    let target_charge_hhrsh =
                        current_state_of_charge + energy_to_add / self.storage_capacity;
                    clip(target_charge_hhrsh, 0., 1.)
                } else {
                    0.
                };

                // target_charge (from input file, or zero when control is off) applied here
                // is treated as an upper limit for target charge
                charge_control
                    .target_charge(simulation_time_iteration, None)
                    .map(|tc| tc.min(target_charge_hhrsh))
            }
            ControlLogicType::HeatBattery => {
                unimplemented!("HeatBattery control logic not implemented for ESH")
            }
        }
    }

    /// Resolve a reference to the owner
    fn owner(&self) -> Arc<dyn HeatBatteryDryCoreCommonBehaviour> {
        // we don't expect this reference to be broken as this and owner should have roughly same lifetimes
        self.owner.upgrade().unwrap()
    }
}

// replicate numpy's clip function
pub(crate) fn clip(n: f64, min: f64, max: f64) -> f64 {
    if n < min {
        min
    } else if n > max {
        max
    } else {
        n
    }
}

struct SocOdeFunction<'a> {
    soc_array: &'a [f64],
    power_array: &'a [f64],
    storage_capacity: f64,
}

impl System<Time, State> for SocOdeFunction<'_> {
    fn system(&self, _x: Time, y: &State, dy: &mut State) {
        // Define the ODE for SOC and energy delivered (no charging, only discharging)

        // Ensure SOC stays within bounds
        let soc = clip(y[0], 0., 1.);

        // Discharging: calculate power used based on SOC
        let discharge_rate = -np_interp(soc, self.soc_array, self.power_array);

        // Track the total energy delivered (discharged energy)
        let ddelivered_dt = -discharge_rate; // Energy delivered (positive value)

        // SOC rate of change (discharging), divided by storage capacity
        let dsoc_dt = -ddelivered_dt / self.storage_capacity;

        dy[0] = dsoc_dt;
    }
}

pub(crate) struct EnergyOutputSocOdeFunction<'a> {
    pub(crate) soc_array: &'a Vec<f64>,
    pub(crate) power_array: &'a Vec<f64>,
    pub(crate) storage_capacity: f64,
    pub(crate) soc_max: f64,
    pub(crate) charge_rate: f64,
    pub(crate) pwr_in: f64,
    pub(crate) target_charge: f64,
}

impl System<Time, EnergyOutputState> for EnergyOutputSocOdeFunction<'_> {
    fn system(&self, _x: Time, y: &EnergyOutputState, dy: &mut EnergyOutputState) {
        // Ensure SOC stays within bounds
        let soc = clip(y[0], 0., 1.);

        // Discharging: calculate power used based on SOC

        // // Python does two interpolations, which we've copied here
        // TODO confirm this is intended
        let power_max_func_result_arr = self
            .soc_array
            .iter()
            .map(|s| np_interp(*s, self.soc_array, self.power_array))
            .collect_vec();
        let discharge_rate =
            -np_interp_with_extrapolate(soc, self.soc_array, &power_max_func_result_arr);

        // Single interpolation version:
        // let discharge_rate =
        //     -np_interp_with_extrapolate(soc, self.soc_array, &self.power_array);

        // Track the total energy delivered (discharged energy)
        let ddelivered_dt = -discharge_rate; // Energy delivered (positive value)

        let dcharged_dt = if soc < self.soc_max {
            self.charge_rate
        } else if self.target_charge > 0. {
            ddelivered_dt.min(self.pwr_in)
        } else {
            0.0
        };

        // Net SOC rate of change (discharge + charge), divided by storage capacity
        let dsoc_dt = (-ddelivered_dt + dcharged_dt) / self.storage_capacity;

        dy[0] = dsoc_dt;
        dy[1] = dcharged_dt;
        dy[2] = ddelivered_dt;
    }

    // TODO implement
    // fn solout(&mut self, _x: Time, y: &EnergyOutputState, _dy: &EnergyOutputState) -> bool {
    //     // TODO we want to check if we've passed this value, not that we are equal to it
    //     // see Emitters for example of this
    //     todo!()
    // }
}

/// A trait that will need to be implemented by any type that uses HeatBatteryDryCore in
/// every case. This is written as the equivalent of the mechanism in the upstream Python
/// whereby there are abstract methods that are called within HeatBatteryDryCore.
pub(crate) trait HeatBatteryDryCoreCommonBehaviour {
    /// Get temperature for charge control calculations
    fn get_temp_for_charge_control(&self) -> Option<f64>;

    /// Get zone setpoint for HHRSH calculations.
    fn get_zone_setpoint(&self) -> f64;

    /// Process energy demand.
    fn demand_energy(&self, energy_demand: f64) -> anyhow::Result<f64>;
}
