//! This module provides the base behaviour for dry core heat storage systems.
//! This includes the common functionality for electrical storage and discharge
//! that is shared between Electric Storage Heaters and Dry Core Heat Batteries.

use crate::core::common::{WaterSupply, WaterSupplyBehaviour};
use crate::core::controls::time_control::{Control, ControlBehaviour};
use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::core::heating_systems::common::HeatingServiceType;
use crate::core::material_properties::WATER;
use crate::core::units::{
    HOURS_PER_DAY, KILOJOULES_PER_KILOWATT_HOUR, SECONDS_PER_HOUR, WATTS_PER_KILOWATT,
};
use crate::core::water_heat_demand::misc::{
    calculate_volume_weighted_average_temperature, water_demand_to_kwh, WaterEventResult,
};
use crate::corpus::{ResultParamValue, ResultsAnnual, ResultsPerTimestep};
use crate::hem_core::simulation_time::SimulationTimeIteration;
use crate::input::{ControlLogicType, HeatBattery};
use crate::statistics::{np_interp, np_interp_with_extrapolate};
use anyhow::bail;
use atomic_float::AtomicF64;
use indexmap::IndexMap;
use itertools::Itertools;
use nalgebra::{Const, OVector, SVector, Vector1};
use ode_solvers::dop_shared::OutputType;
use ode_solvers::{Dopri5, System};
use parking_lot::RwLock;
use smartstring::alias::String;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Weak};

type State = Vector1<f64>;
type Time = f64;

type EnergyOutputState<const D: usize> = SVector<f64, D>;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum OutputMode {
    Min,
    Max,
}

#[derive(Debug)]
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
    owner: Option<Weak<dyn HeatBatteryDryCoreCommonBehaviour>>,
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
            owner: None,
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
        time_remaining: Option<f64>,
        _target_energy: Option<f64>,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64, f64, f64)> {
        let time_remaining = time_remaining.unwrap_or(simulation_time_iteration.timestep);

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

        let mut stepper = create_stepper(
            self.state_of_charge.load(Ordering::SeqCst),
            time_remaining,
            energy_output_soc_ode,
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
        time_remaining: Option<f64>,
        _target_energy: Option<f64>,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64, f64, f64, f64)> {
        let time_remaining = time_remaining.unwrap_or(simulation_time_iteration.timestep);

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
        let energy_output_soc_ode = EnergyOutputWithLossesSocOdeFunction {
            soc_array,
            power_array,
            soc_min_array: &self.soc_min_array,
            output_mode: mode,
            storage_capacity: self.storage_capacity,
            soc_max,
            charge_rate,
            pwr_in: self.pwr_in,
            target_charge,
        };

        let mut stepper = create_stepper(
            self.state_of_charge.load(Ordering::SeqCst),
            time_remaining,
            energy_output_soc_ode,
        );

        // Solve the ODE for SOC and cumulative energy delivered
        let _ = stepper.integrate();

        // Final state of charge after 16 hours
        let final_soc = stepper.y_out().last().unwrap()[0];

        // Total energy charged during the timestep
        let total_energy_charged = stepper.y_out().last().unwrap()[1];

        // Total energy delivered during the timestep
        let total_energy_delivered = stepper.y_out().last().unwrap()[2];

        // Total energy lost during the timestep
        let total_energy_lost = stepper.y_out().last().unwrap()[3];

        // Total time used in delivering energy
        let time_used = stepper.x_out().last().unwrap(); // TODO implement with root finder

        // Return the total energy delivered, time used, and total energy charged
        Ok((
            total_energy_delivered,
            *time_used,
            total_energy_charged,
            final_soc,
            total_energy_lost,
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
                // Implements the "HEAT_BATTERY" control logic
                // A 'high heat retention storage battery' is one with heat retention not less
                // than 45% measured according to BS EN 60531. It incorporates a timer
                // and fan to control the heat output. It is also able to estimate
                // the next day's heating demand based on external temperature, room temperature
                // settings and heat demand periods.

                let energy_to_store = charge_control.energy_to_store(
                    self.demand_met() + self.demand_unmet(),
                    self.get_zone_setpoint(),
                    simulation_time_iteration,
                );

                // Python here allows for energy_to_store to be optional - here it's not, so skipping this logic

                let target_charge_hb = if energy_to_store > 0. {
                    let heat_retention_ratio = self.heat_retention_ratio;

                    let energy_stored = self.state_of_charge() * self.storage_capacity;

                    let energy_to_add = if self.heat_retention_ratio <= 0. {
                        self.storage_capacity - energy_stored
                    } else {
                        (1.0 / heat_retention_ratio) * (energy_to_store - energy_stored)
                    };

                    let target_charge_hb =
                        self.state_of_charge() + energy_to_add / self.storage_capacity;
                    clip(target_charge_hb, 0., 1.)
                } else {
                    0.
                };

                // target_charge (from input, or zero when control is off) applied here
                // is treated as an upper limit for target charge
                charge_control
                    .target_charge(simulation_time_iteration, None)
                    .map(|target_charge| target_charge.min(target_charge_hb))
            }
        }
    }

    /// Resolve a reference to the owner
    fn owner(&self) -> Arc<dyn HeatBatteryDryCoreCommonBehaviour> {
        // we don't expect this reference to be broken as this and owner should have roughly same lifetimes
        self.owner
            .as_ref()
            .expect("storage struct should not be used without having set an owner on it")
            .upgrade()
            .unwrap()
    }

    pub(super) fn set_owner(&mut self, owner: Arc<dyn HeatBatteryDryCoreCommonBehaviour>) {
        self.owner.replace(Arc::downgrade(&owner));
    }

    fn state_of_charge(&self) -> f64 {
        self.state_of_charge.load(Ordering::SeqCst)
    }

    fn set_state_of_charge(&self, soc: f64) {
        self.state_of_charge
            .store(clip(soc, 0.0, 1.0), Ordering::SeqCst);
    }

    fn demand_met(&self) -> f64 {
        self.demand_met.load(Ordering::SeqCst)
    }

    fn set_demand_met(&self, demand_met: f64) {
        self.demand_met.store(demand_met, Ordering::SeqCst);
    }

    fn demand_unmet(&self) -> f64 {
        self.demand_unmet.load(Ordering::SeqCst)
    }

    fn set_demand_unmet(&self, demand_unmet: f64) {
        self.demand_unmet.store(demand_unmet, Ordering::SeqCst);
    }

    fn power_max_func(&self, soc: f64) -> f64 {
        // TODO: confirm this logic - it's quite different from the Python
        np_interp(soc, &self.soc_max_array, &self.power_max_array)
    }

    fn power_min_func(&self, soc: f64) -> f64 {
        // TODO: confirm this logic - it's quite different from the Python
        np_interp(soc, &self.soc_min_array, &self.power_min_array)
    }

    fn get_zone_setpoint(&self) -> f64 {
        self.owner().get_zone_setpoint()
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

fn create_stepper<T, const D: usize>(
    initial_state_of_change: f64,
    time_remaining: f64,
    energy_output_ode_func: T,
) -> Dopri5<Time, OVector<Time, Const<D>>, T>
where
    T: System<Time, EnergyOutputState<D>>,
{
    let current_soc = initial_state_of_change;

    let f = energy_output_ode_func; // f - Structure implementing the System trait
    let x = 0.; // x - Initial value of the independent variable (usually time)
    let x_end = time_remaining; // x_end - Final value of the independent variable
    let dx = 0.; // dx - Increment in the dense output. This argument has no effect if the output type is Sparse
    let mut y0 = EnergyOutputState::<D>::zeros();
    y0[(0, 0)] = current_soc;

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
    Dopri5::from_param(
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
    )
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

type EnergyOutputStateWithoutLosses = EnergyOutputState<3>;

impl System<Time, EnergyOutputStateWithoutLosses> for EnergyOutputSocOdeFunction<'_> {
    fn system(
        &self,
        _x: Time,
        y: &EnergyOutputStateWithoutLosses,
        dy: &mut EnergyOutputStateWithoutLosses,
    ) {
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

pub(crate) struct EnergyOutputWithLossesSocOdeFunction<'a> {
    pub(crate) soc_array: &'a Vec<f64>,
    pub(crate) power_array: &'a Vec<f64>,
    pub(crate) soc_min_array: &'a Vec<f64>,
    pub(crate) output_mode: OutputMode,
    pub(crate) storage_capacity: f64,
    pub(crate) soc_max: f64,
    pub(crate) charge_rate: f64,
    pub(crate) pwr_in: f64,
    pub(crate) target_charge: f64,
}

type EnergyOutputStateWithLosses = EnergyOutputState<4>;

impl System<Time, EnergyOutputStateWithLosses> for EnergyOutputWithLossesSocOdeFunction<'_> {
    fn system(
        &self,
        _x: Time,
        y: &EnergyOutputStateWithLosses,
        dy: &mut EnergyOutputStateWithLosses,
    ) {
        // Ensure SOC stays within bounds
        let soc = clip(y[0], 0., 1.);

        // Calculate the instantaneous loss rate (always based on MIN output)
        let power_min_func_result_arr = self
            .soc_min_array
            .iter()
            .map(|s| np_interp(*s, self.soc_array, self.power_array))
            .collect_vec();
        let loss_rate =
            np_interp_with_extrapolate(soc, self.soc_min_array, &power_min_func_result_arr);
        let mut dlost_dt = loss_rate;

        // Calculate the active discharge rate
        let ddelivered_dt = if self.output_mode == OutputMode::Max {
            let power_max_func_result_arr = self
                .soc_array
                .iter()
                .map(|s| np_interp(*s, self.soc_array, self.power_array))
                .collect_vec();
            let total_discharge_rate =
                np_interp_with_extrapolate(soc, self.soc_array, &power_max_func_result_arr);
            total_discharge_rate - loss_rate
        } else {
            dlost_dt = {
                let power_max_func_result_arr = self
                    .soc_array
                    .iter()
                    .map(|s| np_interp(*s, self.soc_array, self.power_array))
                    .collect_vec();
                np_interp_with_extrapolate(soc, self.soc_array, &power_max_func_result_arr)
            };
            0.0
        };

        let dcharged_dt = if soc < self.soc_max {
            self.charge_rate
        } else if self.target_charge > 0. {
            (ddelivered_dt + dlost_dt).min(self.charge_rate)
        } else {
            0.0
        };

        // Net SOC rate of change
        let dsoc_dt = (-ddelivered_dt - dlost_dt + dcharged_dt) / self.storage_capacity;

        dy[0] = dsoc_dt;
        dy[1] = dcharged_dt;
        dy[2] = ddelivered_dt;
        dy[3] = dlost_dt;
    }
}

/// A trait that will need to be implemented by any type that uses HeatBatteryDryCore in
/// every case. This is written as the equivalent of the mechanism in the upstream Python
/// whereby there are abstract methods that are called within HeatBatteryDryCore.
pub(crate) trait HeatBatteryDryCoreCommonBehaviour: Send + Sync {
    /// Get temperature for charge control calculations
    fn get_temp_for_charge_control(&self) -> Option<f64>;

    /// Get zone setpoint for HHRSH calculations.
    fn get_zone_setpoint(&self) -> f64;

    /// Process energy demand.
    fn demand_energy(
        &self,
        service_name: &str,
        service_type: HeatingServiceType,
        energy_output_required: f64,
        temp_return_feed: f64,
        temp_output: Option<f64>,
        service_on: bool,
        time_start: Option<f64>,
        update_heat_source_state: Option<bool>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64>;

    fn energy_output_max(
        &self,
        _temp_output: f64,
        time_start: Option<f64>,
        simtime: &SimulationTimeIteration,
    ) -> anyhow::Result<f64>;

    fn get_temp_hot_water(&self, inlet_temp: f64, volume: f64, setpoint_temp: f64) -> f64;
}

pub(crate) struct HeatBatteryDryCoreService {
    control: Option<Arc<Control>>,
}

impl HeatBatteryDryCoreService {
    pub(crate) fn new(control: Option<Arc<Control>>) -> Self {
        Self { control }
    }

    fn is_on(&self, simtime: SimulationTimeIteration) -> bool {
        if let Some(control) = &self.control {
            control.is_on(&simtime)
        } else {
            true
        }
    }
}

trait HeatBatteryDryCoreServiceBehaviour {
    fn is_on(&self, simtime: SimulationTimeIteration) -> bool;
}

pub(crate) struct HeatBatteryDryCoreServiceWaterRegular {
    core_service: HeatBatteryDryCoreService,
    heat_battery: Arc<dyn HeatBatteryDryCoreCommonBehaviour>,
    service_name: String,
    cold_feed: WaterSupply,
    control_min: Arc<Control>,
    control_max: Arc<Control>,
}

impl HeatBatteryDryCoreServiceWaterRegular {
    pub(crate) fn new(
        heat_battery: Arc<dyn HeatBatteryDryCoreCommonBehaviour>,
        service_name: String,
        cold_feed: WaterSupply,
        control_min: Arc<Control>,
        control_max: Arc<Control>,
    ) -> Self {
        Self {
            core_service: HeatBatteryDryCoreService::new(control_min.clone().into()),
            heat_battery,
            service_name,
            cold_feed,
            control_min,
            control_max,
        }
    }

    pub(crate) fn setpnt(&self, simtime: SimulationTimeIteration) -> (Option<f64>, Option<f64>) {
        (
            self.control_min.setpnt(&simtime),
            self.control_max.setpnt(&simtime),
        )
    }

    pub(crate) fn demand_energy(
        &self,
        energy_demand: f64,
        temp_flow: f64,
        temp_return: f64,
        update_heat_source_state: Option<bool>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);

        let service_on = self.core_service.is_on(simtime);
        let energy_demand = if service_on { energy_demand } else { 0.0 };

        self.heat_battery.demand_energy(
            self.service_name.as_str(),
            HeatingServiceType::DomesticHotWaterRegular,
            energy_demand,
            temp_return,
            temp_flow.into(),
            service_on,
            None,
            update_heat_source_state.into(),
            simtime,
        )
    }

    /// Calculate the maximum energy output of the heat_battery
    pub(crate) fn energy_output_max(
        &self,
        temp_flow: f64,
        _temp_return: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let service_on = self.core_service.is_on(simtime);
        if !service_on {
            return Ok(0.0);
        }

        self.heat_battery
            .energy_output_max(temp_flow, None, &simtime)
    }
}

/// A struct to represent a direct water heating service provided by a dry core heat battery.
///
/// This is similar to a combi boiler or HIU providing hot water on demand.
pub(crate) struct HeatBatteryDryCoreServiceWaterDirect {
    core_service: HeatBatteryDryCoreService,
    heat_battery: Arc<dyn HeatBatteryDryCoreCommonBehaviour>,
    service_name: String,
    setpoint_temp: f64,
    cold_feed: WaterSupply,
}

impl HeatBatteryDryCoreServiceWaterDirect {
    fn new(
        heat_battery: Arc<dyn HeatBatteryDryCoreCommonBehaviour>,
        service_name: &str,
        setpoint_temp: f64,
        cold_feed: WaterSupply,
    ) -> Self {
        Self {
            core_service: HeatBatteryDryCoreService::new(None),
            heat_battery,
            service_name: service_name.into(),
            setpoint_temp,
            cold_feed,
        }
    }

    pub(crate) fn get_cold_water_source(&self) -> &WaterSupply {
        &self.cold_feed
    }

    /// Return temperature of hot water at outlet
    fn get_temp_hot_water(
        &self,
        volume_req: f64,
        volume_req_already: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<[(Option<f64>, f64); 1]> {
        let volume_req_already = volume_req_already.unwrap_or(0.0);

        if is_close!(volume_req, 0.0, abs_tol = 1e-10) {
            return Ok([(None, volume_req)]);
        }

        let temp_hot_water = |volume| -> anyhow::Result<f64> {
            let list_temp_vol = self.cold_feed.get_temp_cold_water(volume, simtime)?;
            let inlet_temp =
                calculate_volume_weighted_average_temperature(list_temp_vol, volume.into(), None)?;

            Ok(self
                .heat_battery
                .get_temp_hot_water(inlet_temp, volume, self.setpoint_temp))
        };

        let volume_req_cumulative = volume_req + volume_req_already;
        let temp_hot_water_cumulative = temp_hot_water(volume_req_cumulative)?;

        // Base temperature on the part of the draw-off for volume_req, and
        // ignore any volume previously considered
        let temp_hot_water_req = if is_close!(volume_req_already, 0.0, abs_tol = 1e-10) {
            temp_hot_water_cumulative
        } else {
            let temp_hot_water_req_already = temp_hot_water(volume_req_already)?;
            (temp_hot_water_cumulative * volume_req_cumulative
                - temp_hot_water_req_already * volume_req_already)
                / volume_req
        };

        Ok([(temp_hot_water_req.into(), volume_req)])
    }

    /// Process hot water demand directly from dry core heat battery
    pub(crate) fn demand_hot_water(
        &self,
        usage_events: Option<Vec<WaterEventResult>>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let mut energy_demand = 0.0;
        let mut total_volume = 0.0;
        let mut weighted_cold_temp_sum = 0.0;

        if let Some(usage_events) = usage_events {
            for event in usage_events {
                let hot_temp = self.get_temp_hot_water(event.volume_hot, None, simtime)?[0].0;
                let hot_temp = match hot_temp {
                    Some(hot_temp) => hot_temp,
                    None => continue,
                };

                if is_close!(event.volume_hot, 0.0, abs_tol = 1e-10) {
                    continue;
                }

                let list_temp_vol = self.cold_feed.draw_off_water(event.volume_hot, simtime)?;
                let cold_temp = calculate_volume_weighted_average_temperature(
                    list_temp_vol,
                    event.volume_hot.into(),
                    None,
                )?; // this validates the volume

                // Calculate energy needed to heat water
                energy_demand += water_demand_to_kwh(event.volume_hot, hot_temp, cold_temp);

                // Accumulate for weighted average cold water temperature
                total_volume += event.volume_hot;
                weighted_cold_temp_sum += cold_temp * event.volume_hot;
            }
        }

        // Calculate weighted average cold water temperature
        let cold_water_temp = if total_volume > 0.0 {
            weighted_cold_temp_sum / total_volume
        } else {
            // Fallback to sampling method if no events processed
            let cold_water_temp_vol = self.cold_feed.get_temp_cold_water(1.0, simtime)?;
            calculate_volume_weighted_average_temperature(cold_water_temp_vol, 1.0.into(), None)?
        };

        // Demand energy from heat battery
        self.heat_battery.demand_energy(
            self.service_name.as_str(),
            HeatingServiceType::DomesticHotWaterDirect,
            energy_demand,
            cold_water_temp,
            None,
            true,
            None,
            true.into(),
            simtime,
        )
    }
}

/// Wrapper for space heating service from dry core heat battery.
pub(crate) struct HeatBatteryDryCoreServiceSpace {
    core_service: HeatBatteryDryCoreService,
    heat_battery: Arc<dyn HeatBatteryDryCoreCommonBehaviour>,
    service_name: String,
    control: Option<Arc<Control>>,
}

impl HeatBatteryDryCoreServiceSpace {
    fn new(
        heat_battery: Arc<dyn HeatBatteryDryCoreCommonBehaviour>,
        service_name: &str,
        control: Option<Arc<Control>>,
    ) -> Self {
        Self {
            core_service: HeatBatteryDryCoreService::new(control.clone()),
            heat_battery,
            service_name: String::from(service_name),
            control,
        }
    }

    pub(crate) fn temp_setpnt(&self, simtime: SimulationTimeIteration) -> Option<f64> {
        self.control.as_ref().and_then(|c| c.setpnt(&simtime))
    }

    pub(crate) fn in_required_period(&self, simtime: SimulationTimeIteration) -> Option<bool> {
        self.control
            .as_ref()
            .and_then(|c| c.in_required_period(&simtime))
    }

    /// Process space heating demand through dry core heat battery.
    pub(crate) fn demand_energy(
        &self,
        energy_demand: f64,
        temp_flow: f64,
        temp_return: f64,
        _time_start: Option<f64>,
        update_heat_source_state: Option<bool>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);

        let service_on = self.core_service.is_on(simtime);
        if !service_on {
            return Ok(0.0);
        }

        self.heat_battery.demand_energy(
            self.service_name.as_str(),
            HeatingServiceType::Space,
            energy_demand,
            temp_return,
            temp_flow.into(),
            service_on,
            None,
            update_heat_source_state.into(),
            simtime,
        )
    }

    /// Calculate maximum energy output for space heating.
    pub(crate) fn energy_output_max(
        &self,
        temp_output: f64,
        _temp_return_feed: f64,
        time_start: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let time_start = time_start.unwrap_or(0.0);
        if !self.core_service.is_on(simtime) {
            return Ok(0.0);
        }

        self.heat_battery
            .energy_output_max(temp_output, time_start.into(), &simtime)
    }
}

/// Struct to represent dry core heat batteries.
///
/// These batteries use electrical storage similar to ESH but provide
/// heating through water services (space heating via (e.g.) radiators and DHW).
#[derive(Debug)]
struct HeatBatteryDryCore {
    storage: Arc<RwLock<HeatStorageDryCore>>,
    energy_supply: Arc<RwLock<EnergySupply>>,
    energy_supply_connection: Arc<EnergySupplyConnection>,
    energy_supply_connections: Arc<RwLock<IndexMap<String, EnergySupplyConnection>>>,
    fan_power: f64,
    power_instant: f64,
    power_circ_pump: f64,
    power_standby: f64,
    detailed_results: Option<Arc<RwLock<Vec<DetailedResult>>>>,
    service_results: Arc<RwLock<Vec<ServiceResult>>>,
    total_time_running_current_timestep: AtomicF64,
    flag_first_call: AtomicBool,
    battery_losses: AtomicF64,
    simulation_timestep: f64,
}

const ZONE_TEMP_INIT: f64 = 21.0;

impl HeatBatteryDryCore {
    pub(crate) fn new(
        heat_battery_input: HeatBattery,
        charge_control: Arc<Control>,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_connection: Arc<EnergySupplyConnection>,
        n_units: Option<u32>,
        simulation_timestep: f64,
        output_detailed_results: Option<bool>,
    ) -> anyhow::Result<Arc<Self>> {
        let n_units = n_units.unwrap_or(1);
        let output_detailed_results = output_detailed_results.unwrap_or(false);
        let (
            pwr_in,
            storage_capacity,
            dry_core_min_output,
            dry_core_max_output,
            fan_power,
            rated_power_instant,
            state_of_charge_init,
            power_circ_pump,
            power_standby,
        ) = match heat_battery_input {
            HeatBattery::DryCore {
                pwr_in,
                heat_storage_capacity,
                dry_core_min_output,
                dry_core_max_output,
                fan_power,
                rated_power_instant,
                state_of_charge_init,
                electricity_circ_pump,
                electricity_standby,
                ..
            } => (
                pwr_in,
                heat_storage_capacity,
                dry_core_min_output,
                dry_core_max_output,
                fan_power,
                rated_power_instant,
                state_of_charge_init,
                electricity_circ_pump,
                electricity_standby,
            ),
            HeatBattery::Pcm { .. } => {
                unreachable!("HeatBatteryDryCore should only be used for dry core heat batteries")
            }
        };

        let storage = Arc::new(RwLock::new(HeatStorageDryCore::new(
            pwr_in,
            storage_capacity,
            n_units,
            charge_control,
            dry_core_min_output,
            dry_core_max_output,
            state_of_charge_init,
        )?));

        let battery = Arc::new(Self {
            storage: storage.clone(),
            energy_supply,
            energy_supply_connection,
            energy_supply_connections: Arc::new(RwLock::new(IndexMap::new())),
            fan_power,
            power_instant: rated_power_instant,
            power_circ_pump,
            power_standby,
            detailed_results: output_detailed_results.then_some(Arc::new(RwLock::new(vec![]))),
            service_results: Arc::new(RwLock::new(Vec::new())),
            total_time_running_current_timestep: Default::default(),
            flag_first_call: AtomicBool::new(true),
            battery_losses: Default::default(),
            simulation_timestep,
        });

        storage.write().set_owner(battery.clone());

        Ok(battery)
    }

    fn create_service_connection(&self, service_name: &str) -> anyhow::Result<()> {
        if self
            .energy_supply_connections
            .read()
            .contains_key(service_name)
        {
            bail!("Service name already used: {service_name}");
        }

        self.energy_supply_connections.write().insert(
            String::from(service_name),
            EnergySupply::connection(self.energy_supply.clone(), service_name)?,
        );

        Ok(())
    }

    /// Return a HeatBatteryDryCoreServiceWaterRegular object for DHW.
    pub(crate) fn create_service_hot_water_regular(
        battery: Arc<Self>,
        service_name: &str,
        cold_feed: WaterSupply,
        control_min: Arc<Control>,
        control_max: Arc<Control>,
    ) -> anyhow::Result<HeatBatteryDryCoreServiceWaterRegular> {
        battery.create_service_connection(service_name)?;

        Ok(HeatBatteryDryCoreServiceWaterRegular::new(
            battery.clone(),
            service_name.into(),
            cold_feed,
            control_min,
            control_max,
        ))
    }

    /// Return a HeatBatteryDryCoreServiceWaterDirect object and create an EnergySupplyConnection for it
    pub(crate) fn create_service_hot_water_direct(
        battery: Arc<Self>,
        service_name: &str,
        setpoint_temp: f64,
        cold_feed: WaterSupply,
    ) -> anyhow::Result<HeatBatteryDryCoreServiceWaterDirect> {
        battery.create_service_connection(service_name)?;

        Ok(HeatBatteryDryCoreServiceWaterDirect::new(
            battery.clone(),
            service_name,
            setpoint_temp,
            cold_feed,
        ))
    }

    /// Return a HeatBatteryDryCoreServiceSpace object for space heating.
    pub(crate) fn create_service_space_heating(
        battery: Arc<Self>,
        service_name: &str,
        control: Option<Arc<Control>>,
    ) -> anyhow::Result<HeatBatteryDryCoreServiceSpace> {
        battery.create_service_connection(service_name)?;

        Ok(HeatBatteryDryCoreServiceSpace::new(
            battery.clone(),
            service_name,
            control,
        ))
    }

    fn get_battery_losses(&self) -> f64 {
        let battery_losses = self.battery_losses.load(Ordering::SeqCst) * self.n_units() as f64;
        self.battery_losses.store(0.0, Ordering::SeqCst);
        battery_losses
    }

    fn n_units(&self) -> u32 {
        self.storage.read().n_units
    }

    fn time_available(&self, time_start: f64, timestep: f64) -> f64 {
        // Assumes that time spent on other services is evenly spread throughout
        // the timestep so the adjustment for start time below is a proportional
        // reduction of the overall time available, not simply a subtraction
        (timestep
            - self
                .total_time_running_current_timestep
                .load(Ordering::SeqCst))
            * (1.0 - time_start / timestep)
    }

    /// Calculate maximum temperature that can be delivered based on SOC and inlet conditions.
    ///
    /// This method calculates the maximum outlet temperature achievable based on:
    //
    //  Args:
    //    * `inlet_temp`      -- Inlet water temperature (°C)
    //    * `volume`          -- Volume of DHW required (l)
    //    * `setpoint_temp`   -- temperature of hot water to be provided (°C)
    //  Returns:
    //    Maximum outlet temperature achievable (°C)
    fn calculate_max_deliverable_temp(
        &self,
        inlet_temp: f64,
        volume: f64,
        setpoint_temp: f64,
    ) -> f64 {
        // Get maximum power output at current SOC using the accessor method
        let max_power_kw = {
            let storage = self.storage.read();
            storage.power_max_func(storage.state_of_charge())
        };

        let flow_rate_kg_per_s = volume / self.simulation_timestep * WATER.density();

        // Calculate maximum temperature rise using heat transfer equation
        // Q = ṁ × c_p × ΔT
        // Rearranged: ΔT = Q / (ṁ × c_p)
        // Note: specific_heat_capacity_kWh() gives kWh/(kg·K), multiply by 3600 to get kJ/(kg·K)
        let max_outlet_temp = if flow_rate_kg_per_s > 0. {
            let specific_heat_kj_per_kg_k =
                WATER.specific_heat_capacity_kwh() * KILOJOULES_PER_KILOWATT_HOUR as f64;
            let max_temp_rise = max_power_kw
                / (flow_rate_kg_per_s * specific_heat_kj_per_kg_k / WATTS_PER_KILOWATT as f64);
            inlet_temp + max_temp_rise
        } else {
            // No flow means no heat transfer possible
            inlet_temp
        };

        // Cap at maximum design temperature to prevent unrealistic values
        // and ensure system safety limits are respected
        max_outlet_temp.min(setpoint_temp)
    }

    fn energy_output_with_losses(
        &self,
        mode: OutputMode,
        time_remaining: Option<f64>,
        target_energy: Option<f64>,
        simtime: &SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64, f64, f64, f64)> {
        self.storage
            .read()
            .energy_output_with_losses(mode, time_remaining, target_energy, simtime)
    }

    fn energy_output(
        &self,
        mode: OutputMode,
        time_remaining: Option<f64>,
        target_energy: Option<f64>,
        simtime: &SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64, f64, f64)> {
        self.storage
            .read()
            .energy_output(mode, time_remaining, target_energy, simtime)
    }

    /// Calculations to be done at the end of each timestep.
    pub(crate) fn timestep_end(&self, simtime: SimulationTimeIteration) -> anyhow::Result<()> {
        let timestep = simtime.timestep;
        let total_time_running_current_timestep = self
            .total_time_running_current_timestep
            .load(Ordering::SeqCst);
        let time_remaining = timestep - total_time_running_current_timestep;

        // Calculate auxiliary energy
        let mut energy_aux = total_time_running_current_timestep * self.power_circ_pump;
        energy_aux += self.power_standby * time_remaining;
        self.energy_supply_connection
            .demand_energy(energy_aux, simtime.index)?;

        let (energy_charged, final_losses) = if time_remaining > 0. {
            let (_, _, energy_charged, _, final_losses) =
                self.storage.read().energy_output_with_losses(
                    OutputMode::Min,
                    time_remaining.into(),
                    None,
                    &simtime,
                )?;
            (energy_charged, final_losses)
        } else {
            (0., 0.)
        };

        self.battery_losses
            .fetch_add(final_losses, Ordering::SeqCst);

        // save detailed results if required
        if let Some(detailed_results) = self.detailed_results.as_ref() {
            let n_units: f64 = self.n_units() as f64;
            detailed_results.write().push(DetailedResult {
                timestep: simtime.index,
                services: self.service_results.read().clone(),
                soc: self.storage.read().state_of_charge(),
                energy_aux: energy_aux * n_units,
                non_service_energy_lost: final_losses * n_units,
                non_service_charge: energy_charged * n_units,
            });
        }

        self.total_time_running_current_timestep
            .store(0.0, Ordering::SeqCst);
        self.service_results.write().clear();
        self.flag_first_call.store(true, Ordering::SeqCst);

        Ok(())
    }

    /// Output detailed results of heat battery calculation.
    pub(crate) fn output_detailed_results(&self) -> (ResultsPerTimestep, ResultsAnnual) {
        let mut results_per_timestep: ResultsPerTimestep =
            [("auxiliary".into(), Default::default())].into();
        // Report auxiliary parameters (not specific to a service)
        for (parameter, param_unit, _) in AUX_PARAMETERS {
            let auxiliary_results = results_per_timestep["auxiliary"]
                .entry((parameter.into(), param_unit.map(Into::into)))
                .or_default();
            if let Some(detailed_results) =
                self.detailed_results.as_ref().map(|results| results.read())
            {
                for timestep_data in detailed_results.iter() {
                    // For dry core, auxiliary data is stored directly in the timestep map
                    // NB. unlike assumption written out in Python, energy_aux is expected to be available
                    let result = timestep_data.param(parameter);
                    auxiliary_results.push(result);
                }
            }
        }

        // For each service, report required output parameters
        let service_names: Vec<String> = self
            .energy_supply_connections
            .read()
            .keys()
            .cloned()
            .collect();
        for service_name in service_names.iter() {
            results_per_timestep.insert(service_name.clone(), Default::default());
            // Look up each required parameter
            for (parameter, param_unit, _) in OUTPUT_PARAMETERS {
                let param_results_per_timestep = results_per_timestep[service_name]
                    .entry((parameter.into(), param_unit.map(Into::into)))
                    .or_default();
                // Look up value of required parameter in each timestep
                if let Some(detailed_results) =
                    self.detailed_results.as_ref().map(|results| results.read())
                {
                    for timestep_data in detailed_results.iter() {
                        let services_list = &timestep_data.services;
                        let service_data = services_list
                            .iter()
                            .find(|service| service.service_name == *service_name);
                        let result = if let Some(service_data) = service_data {
                            service_data.param(parameter)
                        } else {
                            // Default value if parameter not found
                            match param_unit {
                                Some(_) => ResultParamValue::from(0.0),
                                None => ResultParamValue::Empty,
                            }
                        };
                        param_results_per_timestep.push(result);
                    }
                }
            }
        }

        let mut results_annual: ResultsAnnual = [
            (
                "Overall".into(),
                OUTPUT_PARAMETERS
                    .into_iter()
                    .filter_map(|(parameter, param_unit, incl_in_annual)| {
                        incl_in_annual.then_some((
                            (String::from(parameter), param_unit.map(Into::into)),
                            0.0f64.into(),
                        ))
                    })
                    .collect(),
            ),
            ("auxiliary".into(), Default::default()),
        ]
        .into();

        // Report auxiliary parameters (not specific to a service)
        if let Some(auxiliary_results_annual) = results_annual.get_mut("auxiliary") {
            let auxiliary_results_per_timestep = &results_per_timestep["auxiliary"];
            for (parameter, param_unit, incl_in_annual) in AUX_PARAMETERS {
                if incl_in_annual {
                    auxiliary_results_annual.insert(
                        (String::from(parameter), param_unit.map(Into::into)),
                        auxiliary_results_per_timestep
                            [&(String::from(parameter), param_unit.map(Into::into))]
                            .iter()
                            .cloned()
                            .sum::<ResultParamValue>(),
                    );
                }
            }
        }

        // For each service, report required output parameters
        for service_name in service_names {
            results_annual.insert(service_name.clone(), Default::default());
            for (parameter, param_unit, incl_in_annual) in OUTPUT_PARAMETERS {
                if incl_in_annual {
                    let parameter_annual_total = results_per_timestep[&service_name]
                        [&(parameter.into(), param_unit.map(Into::into))]
                        .iter()
                        .cloned()
                        .sum::<ResultParamValue>();
                    results_annual[&service_name].insert(
                        (parameter.into(), param_unit.map(Into::into)),
                        parameter_annual_total.clone(),
                    );
                    *results_annual["Overall"]
                        .entry((parameter.into(), param_unit.map(Into::into)))
                        .or_insert(ResultParamValue::Number(0.)) += parameter_annual_total;
                }
            }
        }

        (results_per_timestep, results_annual)
    }

    pub(crate) fn set_state_of_charge(&self, state_of_charge: f64) {
        self.storage.read().set_state_of_charge(state_of_charge);
    }

    pub(crate) fn state_of_charge(&self) -> f64 {
        self.storage.read().state_of_charge()
    }

    #[cfg(test)]
    pub(crate) fn get_storage_capacity(&self) -> f64 {
        self.storage.read().storage_capacity
    }

    #[cfg(test)]
    pub(crate) fn get_pwr_in(&self) -> f64 {
        self.storage.read().pwr_in
    }

    #[cfg(test)]
    pub(crate) fn get_demand_met(&self) -> f64 {
        self.storage.read().demand_met()
    }

    #[cfg(test)]
    pub(crate) fn get_demand_unmet(&self) -> f64 {
        self.storage.read().demand_unmet()
    }

    #[cfg(test)]
    pub(crate) fn set_detailed_results(&self, new_results: Vec<DetailedResult>) {
        if let Some(detailed_results) = self.detailed_results.as_ref() {
            *detailed_results.write() = new_results;
        }
    }
}

impl HeatBatteryDryCoreCommonBehaviour for HeatBatteryDryCore {
    fn get_temp_for_charge_control(&self) -> Option<f64> {
        // For heat batteries, return None as they don't have direct zone temperature sensing
        None
    }

    fn get_zone_setpoint(&self) -> f64 {
        ZONE_TEMP_INIT
    }

    fn demand_energy(
        &self,
        service_name: &str,
        service_type: HeatingServiceType,
        energy_output_required: f64,
        temp_return_feed: f64,
        temp_output: Option<f64>,
        service_on: bool,
        time_start: Option<f64>,
        update_heat_source_state: Option<bool>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let time_start = time_start.unwrap_or(0.0);
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);

        let timestep = simtime.timestep;
        let time_remaining = self.time_available(time_start, timestep);

        let energy_output_required = energy_output_required / self.n_units() as f64;
        let mut energy_instant = 0.;
        let mut energy_for_fan = 0.;
        let mut time_running_current_service = 0.;

        // Process energy demand from a specific service
        if !service_on
            || energy_output_required < 0.
            || is_close!(energy_output_required, 0.0, abs_tol = 1e-10)
        {
            if update_heat_source_state {
                self.service_results.write().push(ServiceResult {
                    service_name: service_name.into(),
                    service_type,
                    service_on,
                    energy_output_required: energy_output_required * self.n_units() as f64,
                    temp_output,
                    temp_inlet: temp_return_feed,
                    time_running: 0.,
                    demand_unmet: self.storage.read().demand_unmet() * self.n_units() as f64,
                    energy_delivered_hb: 0.0,
                    energy_delivered_backup: 0.0,
                    energy_delivered_total: 0.0,
                    energy_charged_during_service: 0.0,
                    energy_for_fans: ResultParamValue::Empty,
                    dry_core_soc: self.storage.read().state_of_charge(),
                    current_hb_power: ResultParamValue::Empty,
                    energy_lost: ResultParamValue::Empty,
                });
            }
            return Ok(0.0);
        }

        let (energy_delivered_hb, energy_lost, energy_charged) = if time_remaining < 0.
            || is_close!(time_remaining, 0.0, abs_tol = 1e-10)
        {
            // No time left to run this service
            let energy_delivered_hb = 0.0;
            let energy_lost = 0.0;

            // Update demand tracking
            {
                let storage = self.storage.read();
                storage.set_demand_met(0.0);
                storage.set_demand_unmet(energy_output_required)
            }

            (energy_delivered_hb, energy_lost, 0.0)
        } else {
            // Use the enhanced energy output method with loss tracking

            // First check maximum available energy
            let (q_released_max, time_used_max, energy_charged_max, final_soc, losses_max) =
                self.storage.read().energy_output_with_losses(
                    OutputMode::Max,
                    time_remaining.into(),
                    None,
                    &simtime,
                )?;

            // For DHW direct, no charging during same timestep

            // Determine how to deliver the energy
            let (
                energy_delivered_hb,
                time_running_current_service,
                energy_charged,
                final_soc,
                energy_lost,
            ) = if q_released_max > energy_output_required
                || is_close!(q_released_max, energy_output_required, abs_tol = 1e-10)
            {
                let (
                    energy_delivered_hb,
                    time_running_current_service,
                    energy_charged,
                    final_soc,
                    energy_lost,
                ) = self.storage.read().energy_output_with_losses(
                    OutputMode::Max,
                    time_remaining.into(),
                    energy_output_required.into(),
                    &simtime,
                )?;
                let energy_delivered_hb = energy_delivered_hb.min(energy_output_required);
                (
                    energy_delivered_hb,
                    time_running_current_service,
                    energy_charged,
                    final_soc,
                    energy_lost,
                )
            } else {
                // Not enough energy in storage - deliver what we can
                let energy_delivered_hb = q_released_max;
                time_running_current_service = time_used_max;
                let energy_charged = energy_charged_max;
                let energy_lost = losses_max;

                // Top up with instant heater if available
                if self.power_instant != 0.0 {
                    energy_instant = (energy_output_required - energy_delivered_hb)
                        .min(self.power_instant * time_remaining);
                    let time_instant = energy_instant / self.power_instant;
                    time_running_current_service += time_instant;
                    time_running_current_service = time_running_current_service.min(time_remaining);
                }

                (
                    energy_delivered_hb,
                    time_running_current_service,
                    energy_charged,
                    final_soc,
                    energy_lost,
                )
            };

            // The losses are now accurately integrated during the service delivery
            self.battery_losses.fetch_add(energy_lost, Ordering::SeqCst);

            {
                let storage = self.storage.read();

                // Update state of charge (the ODE has already integrated everything accurately)
                storage.set_state_of_charge(final_soc);

                // Update demand tracking
                storage.set_demand_met(energy_delivered_hb - energy_instant);
                storage.set_demand_unmet(
                    0.0f64.max(energy_output_required - energy_delivered_hb - energy_instant),
                );
            }

            // Calculate fan energy
            energy_for_fan = convert_to_kwh(self.fan_power, time_running_current_service);

            // (from Python) Add energy for fan to internal gains or core or service... TBD

            if update_heat_source_state {
                self.total_time_running_current_timestep
                    .fetch_add(time_running_current_service, Ordering::SeqCst);
            }

            (energy_delivered_hb, energy_lost, energy_charged)
        };

        if update_heat_source_state {
            // Log the energy charged, fan energy, and total energy delivered
            self.energy_supply_connection.demand_energy(
                self.n_units() as f64 * (energy_charged + energy_instant + energy_for_fan),
                simtime.index,
            )?;

            let current_hb_power = if time_running_current_service > 0. {
                energy_delivered_hb * SECONDS_PER_HOUR as f64 / time_running_current_service
            } else {
                0.
            };

            // Record service results with accurate loss tracking
            let n_units: f64 = self.n_units() as f64;
            self.service_results.write().push(ServiceResult {
                service_name: service_name.into(),
                service_type,
                service_on,
                energy_output_required: energy_output_required * n_units,
                temp_output,
                temp_inlet: temp_return_feed,
                time_running: time_running_current_service,
                demand_unmet: self.storage.read().demand_unmet() * n_units,
                energy_delivered_hb: energy_delivered_hb * n_units,
                energy_delivered_backup: energy_instant * n_units,
                energy_delivered_total: (energy_delivered_hb + energy_instant) * n_units,
                energy_charged_during_service: energy_charged * n_units,
                energy_for_fans: (energy_for_fan * n_units).into(),
                dry_core_soc: self.storage.read().state_of_charge(),
                current_hb_power: (current_hb_power * n_units).into(),
                energy_lost: (energy_lost * n_units).into(),
            });
        }

        Ok((energy_delivered_hb + energy_instant) * self.n_units() as f64)
    }

    /// Calculate maximum energy output for current SOC and temperature requirements.
    ///
    /// Args:
    ///   * `temp_output`: Required output temperature (°C)
    ///   * `time_start`: Start time within timestep (unused currently)
    ///
    /// Returns:
    ///   Maximum energy output (kWh)
    fn energy_output_max(
        &self,
        _temp_output: f64,
        time_start: Option<f64>,
        simtime: &SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let time_start = time_start.unwrap_or(0.0);

        let timestep = self.simulation_timestep;
        let time_remaining = self.time_available(time_start, timestep);

        // First, get the base storage calculation which includes charging logic
        let (q_released_max, _, _, _) = self.storage.read().energy_output(
            OutputMode::Max,
            time_remaining.into(),
            None,
            simtime,
        )?;

        let energy_instant = if self.power_instant != 0.0 {
            self.power_instant + time_remaining
        } else {
            0.0
        };

        // Can only provide energy if we can achieve the required output temperature
        Ok((q_released_max + energy_instant) * self.n_units() as f64)
    }

    fn get_temp_hot_water(&self, inlet_temp: f64, volume: f64, setpoint_temp: f64) -> f64 {
        self.calculate_max_deliverable_temp(inlet_temp, volume, setpoint_temp)
    }
}

pub(super) fn convert_to_kwh(power_in_watts: f64, time_in_hours: f64) -> f64 {
    power_in_watts / WATTS_PER_KILOWATT as f64 * time_in_hours
}

#[derive(Debug, Clone)]
struct ServiceResult {
    service_name: String,
    service_type: HeatingServiceType,
    service_on: bool,
    energy_output_required: f64,
    temp_output: Option<f64>,
    temp_inlet: f64,
    time_running: f64,
    demand_unmet: f64,
    energy_delivered_hb: f64,
    energy_delivered_backup: f64,
    energy_delivered_total: f64,
    energy_charged_during_service: f64,
    energy_for_fans: ResultParamValue,
    dry_core_soc: f64,
    current_hb_power: ResultParamValue,
    energy_lost: ResultParamValue,
}

impl ServiceResult {
    fn param(&self, param: &str) -> ResultParamValue {
        match param {
            "service_name" => ResultParamValue::from(self.service_name.clone()),
            "service_type" => ResultParamValue::from(String::from(self.service_type.to_string())),
            "service_on" => self.service_on.into(),
            "energy_output_required" => self.energy_output_required.into(),
            "temp_output" => self.temp_output.into(),
            "temp_inlet" => self.temp_inlet.into(),
            "time_running" => self.time_running.into(),
            "unmet_demand" => self.demand_unmet.into(),
            "energy_delivered_HB" => self.energy_delivered_hb.into(),
            "energy_delivered_backup" => self.energy_delivered_backup.into(),
            "energy_delivered_total" => self.energy_delivered_total.into(),
            "energy_charged_during_service" => self.energy_charged_during_service.into(),
            "energy_for_fans" => self.energy_for_fans.clone(),
            "dry_core_soc" => self.dry_core_soc.into(),
            "current_hb_power" => self.current_hb_power.clone(),
            "energy_lost" => self.energy_lost.clone(),
            _ => panic!("Parameter {param} not recognised"),
        }
    }
}

const OUTPUT_PARAMETERS: [(&str, Option<&str>, bool); 16] = [
    ("service_name", None, false),
    ("service_type", None, false),
    ("service_on", None, false),
    ("energy_output_required", Some("kWh"), true),
    ("temp_output", Some("degC"), false),
    ("temp_inlet", Some("degC"), false),
    ("time_running", Some("secs"), true),
    ("unmet_demand", Some("kWh"), false),
    ("energy_delivered_HB", Some("kWh"), true),
    ("energy_delivered_backup", Some("kWh"), true),
    ("energy_delivered_total", Some("kWh"), true),
    ("energy_charged_during_service", Some("kWh"), true),
    ("energy_for_fans", Some("kWh"), true),
    ("dry_core_soc", Some("ratio"), true),
    ("current_hb_power", Some("kW"), true),
    ("energy_lost", Some("kWh"), true),
];

const AUX_PARAMETERS: [(&str, Option<&str>, bool); 4] = [
    ("energy_aux", Some("kWh"), true),
    ("soc", Some("ratio"), false),
    ("non_service_energy_lost", Some("kWh"), true),
    ("non_service_charge", Some("kWh"), true),
];

#[derive(Debug)]
struct DetailedResult {
    timestep: usize,
    services: Vec<ServiceResult>,
    soc: f64,
    energy_aux: f64,
    non_service_energy_lost: f64,
    non_service_charge: f64,
}

impl DetailedResult {
    fn param(&self, param: &str) -> ResultParamValue {
        match param {
            "energy_aux" => self.energy_aux.into(),
            "soc" => self.soc.into(),
            "non_service_energy_lost" => self.non_service_energy_lost.into(),
            "non_service_charge" => self.non_service_charge.into(),
            _ => panic!("Parameter {param} not recognised"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::common::MockWaterSupply;
    use crate::core::controls::time_control::{ChargeControl, MockControl, SetpointTimeControl};
    use crate::core::water_heat_demand::misc::WaterEventResultType;
    use crate::hem_core::external_conditions::{DaylightSavingsConfig, ExternalConditions};
    use crate::hem_core::simulation_time::SimulationTime;
    use crate::input::{ExternalSensor, ExternalSensorCorrelation, FuelType};
    use approx::assert_relative_eq;
    use rstest::*;

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 5., 1.)
    }

    #[fixture]
    fn schedule() -> Vec<bool> {
        vec![true; 24]
    }

    #[fixture]
    fn external_conditions(simulation_time: SimulationTime) -> Arc<ExternalConditions> {
        ExternalConditions::new(
            &simulation_time.iter(),
            vec![15.0; 24],
            vec![4.0; 24],
            vec![180.; 24],
            vec![100.; 24],
            vec![200.; 24],
            vec![0.2; 24],
            51.5,
            -0.1,
            0,
            0,
            0.into(),
            1.,
            1.into(),
            DaylightSavingsConfig::NotApplicable.into(),
            false,
            false,
            vec![].into(),
        )
        .into()
    }

    #[fixture]
    fn external_sensor() -> ExternalSensor {
        ExternalSensor {
            correlation: vec![
                ExternalSensorCorrelation {
                    temperature: 0.0,
                    max_charge: 1.0,
                },
                ExternalSensorCorrelation {
                    temperature: 10.0,
                    max_charge: 0.9,
                },
                ExternalSensorCorrelation {
                    temperature: 18.0,
                    max_charge: 0.5,
                },
            ],
        }
    }

    #[fixture]
    fn charge_control(
        schedule: Vec<bool>,
        simulation_time: SimulationTime,
        external_conditions: Arc<ExternalConditions>,
        external_sensor: ExternalSensor,
    ) -> Arc<Control> {
        Arc::new(Control::Charge(
            ChargeControl::new(
                ControlLogicType::HeatBattery,
                schedule,
                &simulation_time.iter().current_iteration(),
                0,
                1.,
                vec![Some(1.0), Some(1.8)],
                Some(22.),
                None,
                external_conditions.into(),
                external_sensor.into(),
                None,
            )
            .unwrap(),
        ))
    }

    #[fixture]
    fn charge_control_target_0(
        schedule: Vec<bool>,
        simulation_time: SimulationTime,
        external_conditions: Arc<ExternalConditions>,
        external_sensor: ExternalSensor,
    ) -> Arc<Control> {
        Arc::new(Control::Charge(
            ChargeControl::new(
                ControlLogicType::HeatBattery,
                schedule,
                &simulation_time.iter().current_iteration(),
                0,
                1.,
                vec![Some(0.0), Some(0.0)],
                Some(22.),
                None,
                external_conditions.into(),
                external_sensor.into(),
                None,
            )
            .unwrap(),
        ))
    }

    #[fixture]
    fn energy_supply(simulation_time: SimulationTime) -> Arc<RwLock<EnergySupply>> {
        Arc::new(RwLock::new(
            EnergySupply::new(
                FuelType::Electricity,
                simulation_time.total_steps(),
                None,
                None,
                None,
                None,
            )
            .unwrap(),
        ))
    }

    #[fixture]
    fn energy_supply_connection(
        energy_supply: Arc<RwLock<EnergySupply>>,
    ) -> Arc<EnergySupplyConnection> {
        EnergySupply::connection(energy_supply, "heat_battery")
            .unwrap()
            .into()
    }

    #[fixture]
    fn original_setpoint_temp_water() -> f64 {
        70.
    }

    #[fixture]
    fn heat_battery_input() -> HeatBattery {
        HeatBattery::DryCore {
            control_charge: Default::default(),
            energy_supply: "heat_battery".into(),
            electricity_circ_pump: 0.05,
            electricity_standby: 0.02,
            pwr_in: 3.0,
            state_of_charge_init: 0.0,
            rated_power_instant: 2.0,
            heat_storage_capacity: 12.0,
            number_of_units: 1,
            dry_core_min_output: vec![[0.0, 0.0], [0.5, 0.03], [1.0, 0.06]],
            dry_core_max_output: vec![[0.0, 0.0], [0.5, 2.0], [1.0, 4.0]],
            fan_power: 0.02,
        }
    }

    #[fixture]
    fn heat_battery_input1() -> HeatBattery {
        HeatBattery::DryCore {
            control_charge: Default::default(),
            energy_supply: "heat_battery".into(),
            electricity_circ_pump: 0.05,
            electricity_standby: 0.02,
            pwr_in: 3.0,
            state_of_charge_init: 0.0,
            rated_power_instant: 0.0,
            heat_storage_capacity: 12.0,
            number_of_units: 1,
            dry_core_min_output: vec![[0.0, 0.0], [0.5, 0.03], [1.0, 0.06]],
            dry_core_max_output: vec![[0.0, 0.0], [0.5, 2.0], [1.0, 4.0]],
            fan_power: 0.02,
        }
    }

    #[fixture]
    fn heat_battery(
        heat_battery_input: HeatBattery,
        charge_control: Arc<Control>,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_connection: Arc<EnergySupplyConnection>,
        simulation_time: SimulationTime,
    ) -> Arc<HeatBatteryDryCore> {
        let n_units: u32 = match heat_battery_input {
            HeatBattery::DryCore {
                number_of_units, ..
            } => number_of_units as u32,
            _ => unreachable!(),
        };
        let battery = HeatBatteryDryCore::new(
            heat_battery_input,
            charge_control,
            energy_supply,
            energy_supply_connection,
            n_units.into(),
            simulation_time.step,
            true.into(),
        )
        .unwrap();

        battery.set_state_of_charge(0.7);

        battery
    }

    #[fixture]
    fn heat_battery1(
        heat_battery_input: HeatBattery,
        charge_control: Arc<Control>,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_connection: Arc<EnergySupplyConnection>,
        simulation_time: SimulationTime,
    ) -> Arc<HeatBatteryDryCore> {
        let n_units: u32 = match heat_battery_input {
            HeatBattery::DryCore {
                number_of_units, ..
            } => number_of_units as u32,
            _ => unreachable!(),
        };

        HeatBatteryDryCore::new(
            heat_battery_input,
            charge_control,
            energy_supply,
            energy_supply_connection,
            n_units.into(),
            simulation_time.step,
            true.into(),
        )
        .unwrap()
    }

    fn mock_cold_feed(temperature: Option<f64>) -> WaterSupply {
        WaterSupply::Mock(MockWaterSupply::new(temperature.unwrap_or(10.)))
    }

    #[fixture]
    fn mock_control_dhw() -> Arc<Control> {
        Arc::new(Control::Mock(MockControl::with_is_on(true)))
    }

    #[fixture]
    fn mock_control_dhw_off() -> Arc<Control> {
        Arc::new(Control::Mock(MockControl::with_is_on(false)))
    }

    #[fixture]
    fn mock_control_space() -> Arc<Control> {
        Arc::new(Control::Mock(MockControl::new(
            Some(21.0),
            Some(true),
            Some(true),
        )))
    }

    #[fixture]
    fn default_control_max(simulation_time: SimulationTime) -> Arc<Control> {
        Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(65.), Some(66.)],
            0,
            0.0,
            None,
            None,
            simulation_time.step,
        )))
    }

    // redundant to port Python tests for abstract methods

    #[rstest]
    fn heat_battery_initialization(heat_battery: Arc<HeatBatteryDryCore>) {
        assert_eq!(heat_battery.state_of_charge(), 0.7);
        assert_eq!(heat_battery.get_storage_capacity(), 12.0);
        assert_eq!(heat_battery.get_pwr_in(), 3.0);
        assert!(heat_battery.detailed_results.is_some());
    }

    #[rstest]
    #[ignore = "until ode solving code is corrected in heat battery drycore module"]
    fn test_create_service_hot_water_regular(
        heat_battery: Arc<HeatBatteryDryCore>,
        mock_control_dhw_off: Arc<Control>,
        simulation_time: SimulationTime,
    ) {
        let control_min = Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(45.), Some(46.)],
            0,
            1.,
            None,
            None,
            simulation_time.step,
        )));
        let control_max = Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(65.), Some(66.)],
            0,
            1.,
            None,
            None,
            simulation_time.step,
        )));
        let mock_cold_feed = mock_cold_feed(None); // we can just set up a mock cold water source here - it isn't used
        let service = HeatBatteryDryCore::create_service_hot_water_regular(
            heat_battery.clone(),
            "dhw_service",
            mock_cold_feed.clone(),
            control_min,
            control_max.clone(),
        )
        .unwrap();

        let simtime = simulation_time.iter().current_iteration();

        assert_eq!(service.service_name, "dhw_service");
        let (setpntmin, setpntmax) = service.setpnt(simtime);
        assert_eq!(setpntmin.unwrap(), 45.);
        assert_eq!(setpntmax.unwrap(), 65.);
        // Demand energy
        assert_relative_eq!(
            service.energy_output_max(55., 34., simtime).unwrap(),
            4.829918824420231
        );
        assert_relative_eq!(
            service
                .demand_energy(4.829918824420231, 55., 34., None, simtime)
                .unwrap(),
            4.787470041454905
        );

        let service1 = HeatBatteryDryCore::create_service_hot_water_regular(
            heat_battery,
            "dhw_service1",
            mock_cold_feed,
            mock_control_dhw_off,
            control_max,
        )
        .unwrap();
        // Demand energy
        assert_eq!(service1.energy_output_max(55., 34., simtime).unwrap(), 0.0);
        assert_eq!(
            service1
                .demand_energy(100., 55., 34., None, simtime)
                .unwrap(),
            0.0
        );
    }

    #[rstest]
    #[ignore = "will not resolve until ode solving code is corrected in heat battery drycore module"]
    fn test_heat_battery_dhw_temperature_edge_case(
        charge_control: Arc<Control>,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_connection: Arc<EnergySupplyConnection>,
        simulation_time: SimulationTime,
    ) {
        // Create a heat battery with very low power output capabilities
        let input = HeatBattery::DryCore {
            pwr_in: 0.5,
            heat_storage_capacity: 2.0,
            state_of_charge_init: 0.0,
            dry_core_min_output: vec![[0.0, 0.0], [0.5, 0.001], [1.0, 0.002]], // Very low power
            dry_core_max_output: vec![[0.0, 0.0], [0.5, 2.0], [1.0, 4.0]],     // Very low max power
            fan_power: 0.01,
            rated_power_instant: 0.0, // No instant backup
            energy_supply: "heat_battery".into(),
            number_of_units: 1,
            electricity_circ_pump: 0.01,
            electricity_standby: 0.01,
            control_charge: Default::default(),
        };

        let heat_battery = HeatBatteryDryCore::new(
            input,
            charge_control,
            energy_supply,
            energy_supply_connection,
            Some(1),
            simulation_time.step,
            Some(true),
        )
        .unwrap();

        heat_battery.set_state_of_charge(1.);

        // Create control that requires high temperature
        let _control_min = Arc::new(Control::Mock(MockControl::new(
            Some(40.0),
            Some(true),
            Some(true),
        )));
        let _control_max = Arc::new(Control::Mock(MockControl::new(
            Some(85.0), // High temperature requirement
            None,
            None,
        )));

        // cold feed temperature not relevant

        // Directly call the private __demand_energy method to ensure we hit the right code path
        // This bypasses the service wrapper and gives us more control
        let energy_delivered = heat_battery
            .demand_energy(
                "dhw_high_temp_test",
                HeatingServiceType::DomesticHotWaterRegular,
                2.0,       // Request 2 kWh
                10.0,      // Cold inlet temp
                Some(45.), // Required output temp (high)
                true,
                Some(0.),
                Some(true),
                simulation_time.iter().current_iteration(),
            )
            .unwrap();

        // Energy should be 0 because temperature requirement cannot be met
        // With such low max power (0.02 kW at SOC=0.3), the heat battery
        // cannot raise water temperature from 10°C to 85°C
        assert_eq!(energy_delivered, 2.);

        // Verify demand tracking was updated correctly
        assert_eq!(heat_battery.get_demand_met(), 2.0);
        assert_eq!(heat_battery.get_demand_unmet(), 0.0);
    }

    #[rstest]
    fn test_create_service_space_heating(
        heat_battery: Arc<HeatBatteryDryCore>,
        mock_control_space: Arc<Control>,
        simulation_time: SimulationTime,
    ) {
        let service = HeatBatteryDryCore::create_service_space_heating(
            heat_battery,
            "space_service",
            Some(mock_control_space),
        );

        assert!(service.is_ok());

        let service = service.unwrap();
        let simtime = simulation_time.iter().current_iteration();

        assert_eq!(service.service_name, "space_service");
        assert_eq!(service.temp_setpnt(simtime), Some(21.));
        assert!(service.in_required_period(simtime).unwrap());
    }

    #[rstest]
    #[should_panic = "Service name already used"]
    fn test_duplicate_service_name_error(
        heat_battery: Arc<HeatBatteryDryCore>,
        mock_control_dhw: Arc<Control>,
        default_control_max: Arc<Control>,
        mock_control_space: Arc<Control>,
        simulation_time: SimulationTime,
    ) {
        let mock_cold_feed = mock_cold_feed(None); // we can just set up a normal cold water source here - it isn't used

        // Create first service
        HeatBatteryDryCore::create_service_hot_water_regular(
            heat_battery.clone(),
            "test_service",
            mock_cold_feed,
            mock_control_dhw,
            default_control_max,
        )
        .unwrap();

        // Try to create another service with same name
        HeatBatteryDryCore::create_service_space_heating(
            heat_battery,
            "test_service",
            Some(mock_control_space),
        )
        .unwrap();
    }

    #[rstest]
    fn test_dhw_service_demand_hot_water(
        heat_battery: Arc<HeatBatteryDryCore>,
        simulation_time: SimulationTime,
    ) {
        let mock_cold_feed = mock_cold_feed(None);

        let service = HeatBatteryDryCore::create_service_hot_water_direct(
            heat_battery.clone(),
            "dhw_complex",
            65.,
            mock_cold_feed.clone(),
        )
        .unwrap();

        // skipping check that same cold feed is used as not useful

        // Test with usage events
        let usage_events = vec![
            WaterEventResult {
                event_result_type: WaterEventResultType::Other,
                temperature_warm: 40.0,
                volume_warm: 50.0,
                volume_hot: 8.0,
            },
            WaterEventResult {
                event_result_type: WaterEventResultType::Other,
                temperature_warm: 35.0,
                volume_warm: 0.0,
                volume_hot: 0.0,
            },
        ];

        // Set high SOC to ensure temperature can be met
        heat_battery.set_state_of_charge(0.9);

        let energy = service
            .demand_hot_water(
                usage_events.into(),
                simulation_time.iter().current_iteration(),
            )
            .unwrap();

        assert_relative_eq!(energy, 0.511377777777777, epsilon = 1e-7);
    }

    #[rstest]
    fn test_dhw_service_demand_hot_water_fallback_path(
        heat_battery: Arc<HeatBatteryDryCore>,
        simulation_time: SimulationTime,
    ) {
        let mock_cold_feed = mock_cold_feed(Some(15.));

        let service = HeatBatteryDryCore::create_service_hot_water_direct(
            heat_battery.clone(),
            "dhw_fallback",
            65.,
            mock_cold_feed,
        )
        .unwrap();

        // Set reasonable SOC
        heat_battery.set_state_of_charge(0.5);

        // test with None usage_events (triggers fallback path)
        let energy = service
            .demand_hot_water(None, simulation_time.iter().current_iteration())
            .unwrap();

        // can't verify fallback method was called as mock is not spy

        // Should return 0 energy since no events processed
        assert_eq!(energy, 0.0);
    }

    #[rstest]
    fn test_space_service_demand_energy(
        heat_battery: Arc<HeatBatteryDryCore>,
        mock_control_space: Arc<Control>,
        mock_control_dhw_off: Arc<Control>,
        simulation_time: SimulationTime,
    ) {
        let service = HeatBatteryDryCore::create_service_space_heating(
            heat_battery.clone(),
            "space_service",
            Some(mock_control_space),
        )
        .unwrap();

        // Test normal demand
        let energy = service
            .demand_energy(
                2.,
                60.,
                40.,
                None,
                None,
                simulation_time.iter().current_iteration(),
            )
            .unwrap();

        assert!(energy > 0.);

        // Test with service off
        let service = HeatBatteryDryCore::create_service_space_heating(
            heat_battery,
            "space_service_2",
            Some(mock_control_dhw_off),
        )
        .unwrap();

        let energy_off = service
            .demand_energy(
                2.,
                60.,
                40.,
                None,
                None,
                simulation_time.iter().current_iteration(),
            )
            .unwrap();

        assert_eq!(energy_off, 0.);
    }

    #[rstest]
    fn test_space_service_energy_output_max(
        heat_battery: Arc<HeatBatteryDryCore>,
        mock_control_space: Arc<Control>,
        mock_control_dhw_off: Arc<Control>,
        simulation_time: SimulationTime,
    ) {
        let service = HeatBatteryDryCore::create_service_space_heating(
            heat_battery.clone(),
            "space_service",
            Some(mock_control_space),
        )
        .unwrap();

        // Test normal operation
        let max_energy = service
            .energy_output_max(60., 40., None, simulation_time.iter().current_iteration())
            .unwrap();

        assert!(max_energy > 0.);

        // Test with service off
        let service = HeatBatteryDryCore::create_service_space_heating(
            heat_battery,
            "space_service_2",
            Some(mock_control_dhw_off),
        )
        .unwrap();

        let max_energy_off = service
            .energy_output_max(60., 40., None, simulation_time.iter().current_iteration())
            .unwrap();

        assert_eq!(max_energy_off, 0.);
    }

    // Skipping Python's test_dhw_service_get_temp_hot_water due to mocking

    #[rstest]
    fn test_dhw_service_get_temp_hot_water(
        heat_battery: Arc<HeatBatteryDryCore>,
        heat_battery_input: HeatBattery,
        simulation_time: SimulationTime,
        original_setpoint_temp_water: f64,
    ) {
        let mock_cold_feed_temperature = 20.;
        let mock_cold_feed = mock_cold_feed(mock_cold_feed_temperature.into());

        let service = HeatBatteryDryCore::create_service_hot_water_direct(
            heat_battery.clone(),
            "dhw_complex",
            65.,
            mock_cold_feed,
        )
        .unwrap();

        let hot_water_temp = service
            .get_temp_hot_water(20.0, None, simulation_time.iter().current_iteration())
            .unwrap()[0]
            .0
            .unwrap();

        assert!(hot_water_temp > mock_cold_feed_temperature);
        assert!(hot_water_temp < original_setpoint_temp_water);

        let hot_water_temp = service
            .get_temp_hot_water(20.0, Some(10.), simulation_time.iter().current_iteration())
            .unwrap()[0]
            .0
            .unwrap();
        assert!(hot_water_temp > mock_cold_feed_temperature);
        assert!(hot_water_temp < original_setpoint_temp_water);

        let results = service
            .get_temp_hot_water(0.0, None, simulation_time.iter().current_iteration())
            .unwrap()[0];
        assert!(results.0.is_none());
        assert_eq!(results.1, 0.0);
    }

    // Skipping Python's test_heat_battery_direct_demand_energy_error as not relevant in the Rust

    #[rstest]
    fn test_timestep_end(
        heat_battery: Arc<HeatBatteryDryCore>,
        mock_control_space: Arc<Control>,
        simulation_time: SimulationTime,
    ) {
        let service = HeatBatteryDryCore::create_service_space_heating(
            heat_battery.clone(),
            "space_service",
            Some(mock_control_space),
        )
        .unwrap();

        service
            .demand_energy(
                1.,
                60.,
                0.,
                None,
                None,
                simulation_time.iter().current_iteration(),
            )
            .unwrap();

        heat_battery
            .timestep_end(simulation_time.iter().current_iteration())
            .unwrap();

        assert_eq!(
            heat_battery
                .total_time_running_current_timestep
                .load(Ordering::SeqCst),
            0.,
        );
        assert_eq!(heat_battery.service_results.read().len(), 0);
    }

    // TODO Skipping test_heat_battery_temperature_control as looks unfinished, check in later migrations (unchanged in 1.0.0a2) whether it can be migrated

    #[rstest]
    fn test_electric_charge_heat_battery(
        heat_battery: Arc<HeatBatteryDryCore>,
        simulation_time: SimulationTime,
    ) {
        for t_it in simulation_time.iter() {
            let target_charge = heat_battery
                .storage
                .read()
                .target_electric_charge(t_it)
                .unwrap();

            // Should be a valid charge value between 0 and 1
            assert!(target_charge >= 0.);
            assert!(target_charge <= 1.);
        }
    }

    #[rstest]
    fn test_heat_battery_with_instant_power(
        heat_battery: Arc<HeatBatteryDryCore>,
        mock_control_space: Arc<Control>,
        simulation_time: SimulationTime,
    ) {
        // The heat battery already has instant power configured
        let service = HeatBatteryDryCore::create_service_space_heating(
            heat_battery.clone(),
            "space_service",
            Some(mock_control_space),
        )
        .unwrap();

        // Set low SOC and demand high energy
        heat_battery.set_state_of_charge(0.1);
        let energy = service
            .demand_energy(
                10.0,
                60.0,
                40.0,
                None,
                None,
                simulation_time.iter().current_iteration(),
            )
            .unwrap(); // High demand

        // Should get some energy even with low SOC due to instant backup
        assert!(energy > 0.);
    }

    #[rstest]
    fn test_heat_battery_without_instant_power(
        heat_battery: Arc<HeatBatteryDryCore>,
        mock_control_space: Arc<Control>,
        simulation_time: SimulationTime,
    ) {
        // The heat battery already has instant power configured
        let service = HeatBatteryDryCore::create_service_space_heating(
            heat_battery.clone(),
            "space_service",
            Some(mock_control_space),
        )
        .unwrap();

        // Set low SOC and demand high energy
        heat_battery.set_state_of_charge(0.1);
        let energy = service
            .energy_output_max(50., 10., None, simulation_time.iter().current_iteration())
            .unwrap(); // High demand

        // Should get some energy even with low SOC due to instant backup
        assert!(energy > 0.);
    }

    #[rstest]
    fn test_output_detailed_results(
        heat_battery: Arc<HeatBatteryDryCore>,
        mock_control_space: Arc<Control>,
        simulation_time: SimulationTime,
    ) {
        let mock_cold_feed = mock_cold_feed(Some(10.));

        let dhw_service = HeatBatteryDryCore::create_service_hot_water_direct(
            heat_battery.clone(),
            "dhw_complex",
            65.0,
            mock_cold_feed,
        )
        .unwrap();

        let space_service = HeatBatteryDryCore::create_service_space_heating(
            heat_battery.clone(),
            "space_service",
            mock_control_space.into(),
        )
        .unwrap();

        // Simulate a few timesteps

        // the variable hot_water_energy in the upstream Python is unnecessary as that parameter does not get used

        let mut timestep_idx: usize = 0;
        for t_it in simulation_time.iter() {
            if timestep_idx >= 3 {
                break;
            }

            // DHW demand
            let usage_events = vec![WaterEventResult {
                event_result_type: WaterEventResultType::Other,
                temperature_warm: 40.0,
                volume_warm: 30.0,
                volume_hot: 20.0,
            }];
            let _dhw_energy = dhw_service
                .demand_hot_water(usage_events.into(), t_it)
                .unwrap();

            // Space heating demand
            space_service
                .demand_energy(1.5, 55.0, 45.0, None, None, t_it)
                .unwrap();

            // End timestep
            heat_battery.timestep_end(t_it).unwrap();
            timestep_idx += 1;
        }

        let (results_per_timestep, results_annual) = heat_battery.output_detailed_results();

        // Check structure of results
        assert!(results_per_timestep.contains_key("auxiliary"));
        assert!(results_per_timestep.contains_key("dhw_complex"));
        assert!(results_per_timestep.contains_key("space_service"));

        // Check annual results
        assert!(results_annual.contains_key("Overall"));
        assert!(results_annual.contains_key("auxiliary"));
        assert!(results_annual.contains_key("dhw_complex"));
        assert!(results_annual.contains_key("space_service"));

        // unnecessary to check other test case from Python as parameters to output_detailed_results are not used
    }

    // test_output_detailed_results_without_all_aux_parameter_keys is unneeded as case it is testing is impossible in the Rust

    #[rstest]
    fn test_heat_battery_without_detailed_results(
        heat_battery_input: HeatBattery,
        charge_control: Arc<Control>,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_connection: Arc<EnergySupplyConnection>,
    ) {
        let heat_battery = HeatBatteryDryCore::new(
            heat_battery_input,
            charge_control,
            energy_supply,
            energy_supply_connection,
            Some(1),
            1.0,
            Some(false),
        )
        .unwrap();

        // Check that detailed results are None
        assert!(heat_battery.detailed_results.is_none());
    }

    #[rstest]
    fn test_multiple_units(
        heat_battery_input: HeatBattery,
        charge_control: Arc<Control>,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_connection: Arc<EnergySupplyConnection>,
        mock_control_space: Arc<Control>,
        simulation_time: SimulationTime,
    ) {
        let heat_battery = HeatBatteryDryCore::new(
            heat_battery_input,
            charge_control,
            energy_supply,
            energy_supply_connection,
            Some(3),
            1.0,
            Some(false),
        )
        .unwrap();

        heat_battery.calculate_max_deliverable_temp(5., 0.0, 65.0);

        let service = HeatBatteryDryCore::create_service_space_heating(
            heat_battery.clone(),
            "space_service",
            mock_control_space.into(),
        )
        .unwrap();

        // Energy output should scale with number of units
        let energy = service
            .demand_energy(
                2.0,
                60.0,
                40.0,
                None,
                None,
                simulation_time.iter().current_iteration(),
            )
            .unwrap();
        assert!(energy > 0.);
    }

    #[rstest]
    fn test_base_service_class_without_control(simulation_time: SimulationTime) {
        let base_service = HeatBatteryDryCoreService::new(None);

        // Should return true when no control
        assert!(base_service.is_on(simulation_time.iter().current_iteration()));
    }

    #[rstest]
    fn test_get_temp_for_charge_control(heat_battery: Arc<HeatBatteryDryCore>) {
        let temperature = heat_battery.get_temp_for_charge_control();
        assert!(temperature.is_none());
    }

    #[rstest]
    fn test_get_zone_setpoint(heat_battery: Arc<HeatBatteryDryCore>) {
        let setpoint = heat_battery.get_zone_setpoint();
        assert_eq!(setpoint, 21.);
    }

    // test_heat_battery_charge_with_nonpositive_energy_to_store,
    // test_heat_battery_charge_zero_heat_retention, and
    // test_heat_battery_charge_no_heat_retention_ratio
    // are tricksier to port
    // (as they call on ChargeControl's API (energy_to_store)) than is beneficial

    #[rstest]
    fn test_heat_battery_dry_core_service_off_conditions(
        heat_battery_input: HeatBattery,
        charge_control: Arc<Control>,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_connection: Arc<EnergySupplyConnection>,
        simulation_time: SimulationTime,
    ) {
        let heat_battery = HeatBatteryDryCore::new(
            heat_battery_input,
            charge_control,
            energy_supply,
            energy_supply_connection,
            Some(1),
            1.0,
            Some(false),
        )
        .unwrap();

        let ctrl_off = Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            vec![None; 24],
            0,
            1.0,
            None,
            None,
            1.0,
        )));

        let service = HeatBatteryDryCore::create_service_space_heating(
            heat_battery.clone(),
            "test_space",
            ctrl_off.into(),
        )
        .unwrap();

        // Test energy_output_max when service is off
        let result = service
            .energy_output_max(
                50.0,
                40.0,
                Some(0.0),
                simulation_time.iter().current_iteration(),
            )
            .unwrap();
        assert_eq!(result, 0.0);
    }

    #[rstest]
    fn test_heat_battery_dry_core_no_detailed_results(
        heat_battery_input: HeatBattery,
        charge_control: Arc<Control>,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_connection: Arc<EnergySupplyConnection>,
        simulation_time: SimulationTime,
    ) {
        let heat_battery = HeatBatteryDryCore::new(
            heat_battery_input,
            charge_control,
            energy_supply,
            energy_supply_connection,
            Some(1),
            1.0,
            Some(false),
        )
        .unwrap();

        // Run through timestep - no parameters
        heat_battery
            .timestep_end(simulation_time.iter().current_iteration())
            .unwrap();

        // Test that detailed results are None (internal attribute)
        assert!(heat_battery.detailed_results.is_none());
    }

    #[rstest]
    fn test_heat_battery_dry_core_water_service_edge_cases(
        heat_battery_input: HeatBattery,
        charge_control: Arc<Control>,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_connection: Arc<EnergySupplyConnection>,
        simulation_time: SimulationTime,
    ) {
        let mock_cold_feed = mock_cold_feed(Some(10.));

        let heat_battery = HeatBatteryDryCore::new(
            heat_battery_input,
            charge_control,
            energy_supply,
            energy_supply_connection,
            Some(1),
            1.0,
            Some(false),
        )
        .unwrap();

        // Create DHW service
        let service = HeatBatteryDryCore::create_service_hot_water_direct(
            heat_battery.clone(),
            "dhw_complex",
            65.0, // High setpoint
            mock_cold_feed,
        )
        .unwrap();

        // Test edge cases

        // Test with zero volume request
        let temp_vol_list = service
            .get_temp_hot_water(0.0, None, simulation_time.iter().current_iteration())
            .unwrap();
        assert_eq!(temp_vol_list[0].0, None);
        assert_eq!(temp_vol_list[0].1, 0.0);

        // Test with very small volume (close to zero but not exactly zero)
        let temperature_volume = service
            .get_temp_hot_water(0.0000001, None, simulation_time.iter().current_iteration())
            .unwrap();
        assert_eq!(temperature_volume.len(), 1);

        // Test with normal volume request - reset mock to return correct volume
        let temp_vol_list = service
            .get_temp_hot_water(10.0, None, simulation_time.iter().current_iteration())
            .unwrap();
        let (_temperature, volume) = temp_vol_list[0];
        assert_eq!(volume, 10.0);
    }

    // test_heat_battery_dry_core_extreme_temperatures is incomplete/ tests nothing, so not ported

    /// Test HeatBatteryDryCore with multiple services demanding energy
    #[rstest]
    #[ignore = "will not pass until ode solver code is correct"]
    fn test_heat_battery_dry_core_multiple_services_interaction(
        heat_battery_input: HeatBattery,
        charge_control: Arc<Control>,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_connection: Arc<EnergySupplyConnection>,
        mock_control_space: Arc<Control>,
        mock_control_dhw: Arc<Control>,
        default_control_max: Arc<Control>,
        simulation_time: SimulationTime,
    ) {
        let simtime = simulation_time.iter().current_iteration();

        let heat_battery = HeatBatteryDryCore::new(
            heat_battery_input,
            charge_control,
            energy_supply,
            energy_supply_connection,
            Some(2),
            1.0,
            Some(false),
        )
        .unwrap();

        let service1 = HeatBatteryDryCore::create_service_space_heating(
            heat_battery.clone(),
            "space1",
            mock_control_space.into(),
        );

        let service2 = HeatBatteryDryCore::create_service_hot_water_regular(
            heat_battery.clone(),
            "dhw1",
            mock_cold_feed(None),
            mock_control_dhw,
            default_control_max,
        );

        // Test services work correctly (i.e. we can access them
        let service1 = if let Ok(service) = service1 {
            service
        } else {
            panic!("Failed to create service1")
        };
        if service2.is_err() {
            panic!("Failed to create service2")
        };

        // energy_output_max from space service
        let result1 = service1.energy_output_max(35., 10., None, simtime).unwrap();
        assert_relative_eq!(result1, 4.897563801409152, epsilon = 1e-7);

        // demand energy from services
        let result1 = service1
            .demand_energy(1.0, 50.0, 40.0, None, None, simtime)
            .unwrap();
        assert!(result1 >= 0.0);
    }

    // test_heat_battery_dry_core_base_service_is_on_without_control duplicates test_base_service_class_without_control

    #[rstest]
    fn test_heat_battery_dry_core_get_temp_hot_water_edge_cases(
        heat_battery_input: HeatBattery,
        charge_control: Arc<Control>,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_connection: Arc<EnergySupplyConnection>,
        simulation_time: SimulationTime,
    ) {
        let simtime = simulation_time.iter().current_iteration();

        let heat_battery = HeatBatteryDryCore::new(
            heat_battery_input,
            charge_control,
            energy_supply,
            energy_supply_connection,
            Some(1),
            1.0,
            Some(false),
        )
        .unwrap();

        let mock_cold_feed_10 = mock_cold_feed(Some(10.));

        // Test regular DHW service
        let service_regular = HeatBatteryDryCore::create_service_hot_water_direct(
            heat_battery.clone(),
            "dhw_complex",
            65.0,
            mock_cold_feed_10,
        )
        .unwrap();

        // Test with very small volume (close to zero but not exactly zero)
        let temperature_volume = service_regular
            .get_temp_hot_water(0.0000001, None, simtime)
            .unwrap();
        assert_eq!(temperature_volume.len(), 1);

        // Test with multiple temperature/volume pairs (though only exercising the first simulation timestep)
        let mock_cold_feed_8 = mock_cold_feed(Some(8.0));
        let temperature_volume = service_regular
            .get_temp_hot_water(10.0, None, simtime)
            .unwrap();
        assert_eq!(temperature_volume.len(), 1);
        assert_eq!(temperature_volume[0].1, 10.0);
    }

    // skip test_heat_battery_dry_core_water_regular_demand_edge_cases as uses spies

    #[rstest]
    fn test_heat_battery_dry_core_output_detailed_results_edge_cases(
        heat_battery_input: HeatBattery,
        charge_control: Arc<Control>,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_connection: Arc<EnergySupplyConnection>,
        mock_control_space: Arc<Control>,
        simulation_time: SimulationTime,
    ) {
        let heat_battery = HeatBatteryDryCore::new(
            heat_battery_input,
            charge_control,
            energy_supply,
            energy_supply_connection,
            Some(2),
            1.0,
            Some(false),
        )
        .unwrap();

        // Create and use a service to generate some results
        let service = HeatBatteryDryCore::create_service_space_heating(
            heat_battery.clone(),
            "space_service",
            mock_control_space.into(),
        )
        .unwrap();

        let simtime = simulation_time.iter().current_iteration();

        // Demand energy multiple times
        for i in 0..3 {
            service
                .demand_energy(0.5 + i as f64 * 0.1, 50.0, 40.0, None, None, simtime)
                .unwrap();
        }

        // Call timestep_end to save results
        heat_battery.timestep_end(simtime).unwrap();

        // Use the correct way to get number of timesteps
        let num_timesteps = simulation_time.total_steps();

        // skipping setting up different lengths of inputs for output_detailed_results as these are ignored/unused
        let (results1, _) = heat_battery.output_detailed_results();

        assert!(results1.contains_key("auxiliary"));
    }

    #[rstest]
    fn test_heat_battery_dry_core_timestep_end_with_charging(
        heat_battery_input: HeatBattery,
        charge_control: Arc<Control>,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_connection: Arc<EnergySupplyConnection>,
        simulation_time: SimulationTime,
    ) {
        let simtime = simulation_time.iter().current_iteration();

        let mut heat_battery = HeatBatteryDryCore::new(
            heat_battery_input,
            charge_control,
            energy_supply,
            energy_supply_connection,
            Some(1),
            1.0,
            Some(true),
        )
        .unwrap();

        // Set initial state
        heat_battery.set_state_of_charge(0.3);

        // Don't run any services, so time_remaining > 0
        // This should trigger the charging logic in timestep_end

        heat_battery.timestep_end(simtime).unwrap();
        let final_soc = heat_battery.state_of_charge();

        // SOC might increase due to charging (depends on charge control logic)
        assert!(final_soc >= 0.0);
        assert!(final_soc <= 1.0);

        // Time available
        assert_eq!(heat_battery.time_available(0., 1.), 1.);

        heat_battery
            .total_time_running_current_timestep
            .store(1., Ordering::SeqCst);
        heat_battery.set_detailed_results(vec![]);
        heat_battery.timestep_end(simtime).unwrap();

        // Get detailed results
        let (_, _) = heat_battery.output_detailed_results();

        assert_eq!(
            heat_battery.detailed_results.as_ref().unwrap().read()[0].non_service_energy_lost,
            0.0
        );
    }

    // skipped test_heat_battery_dry_core_demand_energy_not_implemented as this is enforced statically by the Rust compiler
}
