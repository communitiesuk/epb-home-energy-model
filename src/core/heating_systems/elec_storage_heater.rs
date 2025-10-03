use crate::core::controls::time_control::{per_control, ControlBehaviour};
use crate::{
    core::{
        controls::time_control::Control, energy_supply::energy_supply::EnergySupplyConnection,
        units::WATTS_PER_KILOWATT,
    },
    external_conditions::ExternalConditions,
    input::{ControlLogicType, ElectricStorageHeaterAirFlowType},
    simulation_time::{SimulationTimeIteration, SimulationTimeIterator},
    statistics::{np_interp, np_interp_with_extrapolate},
};
use anyhow::bail;
use atomic_float::AtomicF64;
use derivative::Derivative;
use itertools::Itertools;
use nalgebra::{Vector1, Vector3};
use ode_solvers::{dop_shared::OutputType, Dopri5, System};
use parking_lot::RwLock;
use std::sync::atomic::Ordering;
use std::sync::Arc;

type State = Vector1<f64>;
type Time = f64;

type EnergyOutputState = Vector3<f64>;

// replicates numpys's linspace function
fn linspace(start: f64, end: f64, num: i32) -> Vec<f64> {
    let step = (end - start) / f64::from(num - 1);
    (0..num).map(|n| start + (f64::from(n) * step)).collect()
}

// replicate numpy's clip function
fn clip(n: f64, min: f64, max: f64) -> f64 {
    if n < min {
        min
    } else if n > max {
        max
    } else {
        n
    }
}

#[derive(Derivative)]
#[derivative(Debug)]
pub(crate) struct ElecStorageHeater {
    pwr_in: f64,
    pwr_instant: f64,
    storage_capacity: f64,
    air_flow_type: ElectricStorageHeaterAirFlowType,
    frac_convective: f64,
    n_units: i32,
    energy_supply_conn: EnergySupplyConnection,
    control: Arc<Control>,
    charge_control: Arc<Control>,
    fan_pwr: f64,
    external_conditions: Arc<ExternalConditions>,
    temp_air: f64,
    state_of_charge: AtomicF64,
    demand_met: AtomicF64,
    demand_unmet: AtomicF64,
    zone_setpoint_init: f64,
    #[derivative(Debug = "ignore")]
    zone_internal_air_func: Arc<dyn Fn() -> f64 + Send + Sync>,
    soc_max_array: Vec<f64>,
    power_max_array: Vec<f64>,
    soc_min_array: Vec<f64>,
    power_min_array: Vec<f64>,
    heat_retention_ratio: f64,
    current_energy_profile: RwLock<CurrentEnergyProfile>,
    esh_detailed_results: Option<Arc<RwLock<Vec<StorageHeaterDetailedResult>>>>,
}

#[derive(Clone, Copy, Debug, Default)]
/// A struct to encapsulate energy values for a current step.
struct CurrentEnergyProfile {
    energy_for_fan: f64,
    energy_instant: f64,
    energy_charged: f64,
    energy_delivered: f64,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct StorageHeaterDetailedResult {
    timestep_idx: usize,
    n_units: i32,
    energy_demand: f64,
    energy_delivered: f64,
    energy_instant: f64,
    energy_charged: f64,
    energy_for_fan: f64,
    state_of_charge: f64,
    final_soc: f64,
    time_used_max: f64,
}

impl StorageHeaterDetailedResult {
    pub(crate) fn as_string_values(&self) -> Vec<String> {
        vec![
            self.timestep_idx.to_string(),
            self.n_units.to_string(),
            self.energy_demand.to_string(),
            self.energy_delivered.to_string(),
            self.energy_instant.to_string(),
            self.energy_charged.to_string(),
            self.energy_for_fan.to_string(),
            self.state_of_charge.to_string(),
            self.final_soc.to_string(),
            self.time_used_max.to_string(),
        ]
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

struct EnergyOutputSocOdeFunction<'a> {
    soc_array: &'a Vec<f64>,
    power_array: &'a Vec<f64>,
    storage_capacity: f64,
    soc_max: f64,
    charge_rate: f64,
    pwr_in: f64,
    target_charge: f64,
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

pub(crate) enum OutputMode {
    Min,
    Max,
}

impl ElecStorageHeater {
    /// Arguments:
    /// * `pwr_in` - in kW (Charging)
    /// * `rated_power_instant` - in kW (Instant backup)
    /// * `storage_capacity` - in kWh
    /// * `air_flow_type` - string specifying type of Electric Storage Heater:
    ///                   - fan-assisted
    ///                   - damper-only
    /// * `frac_convective`      - convective fraction for heating (TODO: Check if necessary)
    /// * `fan_pwr`              - Fan power [W]
    /// * `n_units`              - number of units install in zone
    /// * `zone_internal_air_func`  - function that provides access to the internal air temp of the zone
    /// * `energy_supply_conn`   - reference to EnergySupplyConnection object
    /// * `simulation_time`      - reference to SimulationTime object
    /// * `control`              - reference to a control object which must implement is_on() and setpnt() funcs
    /// * `charge_control`       - reference to a ChargeControl object which must implement different logic types
    ///                         for charging the Electric Storage Heaters.
    /// * `esh_min_output`       - Data from test showing the output from the storage heater when not actively
    ///                         outputting heat, i.e. case losses only (with units kW)
    /// * `esh_max_output`       - Data from test showing the output from the storage heater when it is actively
    ///                         outputting heat, e.g. damper open / fan running (with units kW)
    /// * `external_conditions`  - reference to ExternalConditions object
    pub(crate) fn new(
        pwr_in: f64,
        rated_power_instant: f64,
        storage_capacity: f64,
        air_flow_type: ElectricStorageHeaterAirFlowType,
        frac_convective: f64,
        fan_pwr: f64,
        n_units: i32,
        zone_setpoint_init: f64,
        zone_internal_air_func: Arc<dyn Fn() -> f64 + Send + Sync>,
        energy_supply_conn: EnergySupplyConnection,
        simulation_time: &SimulationTimeIterator,
        control: Arc<Control>,
        charge_control: Arc<Control>,
        esh_min_output: Vec<(f64, f64)>,
        esh_max_output: Vec<(f64, f64)>,
        external_conditions: Arc<ExternalConditions>,
        output_detailed_results: Option<bool>,
    ) -> anyhow::Result<Self> {
        let output_detailed_results = output_detailed_results.unwrap_or(false);

        if !matches!(charge_control.as_ref(), Control::Charge(_)) {
            bail!("charge_control must be a ChargeControl");
        }

        let temp_air = zone_internal_air_func();

        // Convert ESH_max_output to NumPy arrays without sorting
        let soc_max_array = esh_max_output.iter().map(|f| f.0).collect_vec();
        let power_max_array = esh_max_output.iter().map(|f| f.1).collect_vec();

        // Convert ESH_min_output to NumPy arrays without sorting
        let soc_min_array = esh_min_output.iter().map(|f| f.0).collect_vec();
        let power_min_array = esh_min_output.iter().map(|f| f.1).collect_vec();

        // Validate that both SOC arrays are in strictly increasing order
        if !soc_max_array
            .iter()
            .zip(soc_max_array.iter().skip(1))
            .all(|(a, b)| a <= b)
        {
            bail!("esh_max_output SOC values must be in increasing order (from 0.0 to 1.0).");
        }

        if !soc_min_array
            .iter()
            .zip(soc_min_array.iter().skip(1))
            .all(|(a, b)| a <= b)
        {
            bail!("esh_min_output SOC values must be in increasing order (from 0.0 to 1.0).");
        }

        // Validate that both SOC arrays start at 0.0 and end at 1.0
        if !is_close!(*soc_max_array.first().unwrap(), 0.) {
            bail!("The first SOC value in esh_max_output must be 0.0 (fully discharged).");
        }

        if !is_close!(*soc_max_array.last().unwrap(), 1.) {
            bail!("The last SOC value in esh_max_output must be 1.0 (fully charged).");
        }

        if !is_close!(*soc_min_array.first().unwrap(), 0.) {
            bail!("The first SOC value in esh_min_output must be 0.0 (fully discharged).");
        }

        if !is_close!(*soc_min_array.last().unwrap(), 1.) {
            bail!("The last SOC value in esh_min_output must be 1.0 (fully charged).");
        }

        // Validate that for any SOC, power_max >= power_min
        // Sample a fine grid of SOCs and ensure power_max >= power_min
        let fine_soc: Vec<f64> = linspace(0., 1., 100);

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
                bail!("At all SOC levels, ESH_max_output must be >= ESH_min_output.")
            }
        }

        // TODO can we pass these vecs/arrays by reference instead
        let heat_retention_ratio =
            Self::heat_retention_output(&soc_min_array, &power_min_array, storage_capacity);

        Ok(Self {
            pwr_in,
            pwr_instant: rated_power_instant,
            storage_capacity,
            air_flow_type,
            frac_convective,
            n_units,
            energy_supply_conn,
            control,
            charge_control,
            fan_pwr,
            external_conditions,
            temp_air,
            state_of_charge: Default::default(),
            demand_met: Default::default(),
            demand_unmet: Default::default(),
            zone_setpoint_init,
            zone_internal_air_func,
            soc_max_array,
            power_max_array,
            soc_min_array,
            power_min_array,
            heat_retention_ratio,
            current_energy_profile: Default::default(),
            esh_detailed_results: output_detailed_results.then(|| {
                Arc::new(RwLock::new(Vec::with_capacity(
                    simulation_time.total_steps(),
                )))
            }),
        })
    }

    pub(crate) fn temp_setpnt(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<f64> {
        per_control!(self.control.as_ref(), ctrl => { ctrl.setpnt(simulation_time_iteration) })
    }

    pub(crate) fn in_required_period(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> Option<bool> {
        per_control!(self.control.as_ref(), ctrl => { ctrl.in_required_period(simulation_time_iteration) })
    }

    pub(crate) fn frac_convective(&self) -> f64 {
        self.frac_convective
    }

    fn convert_to_kwh(power: f64, time: f64) -> f64 {
        // Converts power value supplied to the correct energy unit
        // Arguments
        // power -- Power value in watts
        // time -- length of the time active
        // returns -- Energy in kWh
        power / f64::from(WATTS_PER_KILOWATT) * time
    }

    pub fn heat_retention_output(
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

    pub fn energy_output_min(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // Calculates the minimum energy that must be delivered based on ESH_min_output.
        // :return: Tuple containing (minimum energy deliverable in kWh, time used in hours).
        let (soc, _, _, _) = self.energy_output(OutputMode::Min, simulation_time_iteration)?;
        Ok(soc * f64::from(self.n_units))
    }

    #[cfg(test)]
    pub(crate) fn energy_output_max(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64, f64, f64)> {
        // Calculates the maximum energy that can be delivered based on ESH_max_output.
        // :return: Tuple containing (maximum energy deliverable in kWh, time used in hours).
        self.energy_output(OutputMode::Max, simulation_time_iteration)
    }

    pub(crate) fn demand_energy(
        &self,
        energy_demand: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // Determines the amount of energy to release based on energy demand, while also handling the
        // energy charging and logging fan energy.
        // :param energy_demand: Energy demand in kWh.
        // :return: Total net energy delivered (including instant heating and fan energy).
        let mut current_profile = self.current_energy_profile.write();

        let timestep = simulation_time_iteration.timestep;
        let energy_demand = energy_demand / f64::from(self.n_units);
        current_profile.energy_instant = 0.;

        // Initialize time_used_max and energy_charged_max to default values
        let mut time_used_max = 0.;
        let _energy_charged_max = 0.;

        // Calculate minimum energy that can be delivered
        let (q_released_min, _, energy_charged, mut final_soc) =
            self.energy_output(OutputMode::Min, simulation_time_iteration)?;
        current_profile.energy_charged = energy_charged;

        let mut q_released_max: Option<f64> = None;

        if q_released_min > energy_demand {
            // Deliver at least the minimum energy
            current_profile.energy_delivered = q_released_min;
            self.demand_met.store(q_released_min, Ordering::SeqCst);
            self.demand_unmet.store(0., Ordering::SeqCst);
        } else {
            // Calculate maximum energy that can be delivered
            let (q_released_max_value, time_used_max_tmp, energy_charged, final_soc_override) =
                self.energy_output(OutputMode::Max, simulation_time_iteration)?;
            final_soc = final_soc_override;

            q_released_max = Some(q_released_max_value);
            time_used_max = time_used_max_tmp;
            current_profile.energy_charged = energy_charged;

            if q_released_max_value < energy_demand {
                // Deliver as much as possible up to the maximum energy
                current_profile.energy_delivered = q_released_max_value;
                self.demand_met
                    .store(q_released_max_value, Ordering::SeqCst);
                self.demand_unmet
                    .store(energy_demand - q_released_max_value, Ordering::SeqCst);

                // For now, we assume demand not met from storage is topped-up by
                // the direct top-up heater (if applicable). If still some unmet,
                // this is reported as unmet demand.
                if self.pwr_instant != 0. {
                    current_profile.energy_instant = self
                        .demand_unmet
                        .load(Ordering::SeqCst)
                        .min(self.pwr_instant * timestep); // kWh
                    let time_instant = current_profile.energy_instant / self.pwr_instant;
                    time_used_max += time_instant;
                    time_used_max = time_used_max.min(timestep);
                }
            } else {
                // Deliver the demanded energy
                current_profile.energy_delivered = energy_demand;

                if q_released_max_value > 0. {
                    time_used_max *= energy_demand / q_released_max_value;
                }

                self.demand_met.store(energy_demand, Ordering::SeqCst);
                self.demand_unmet.store(0., Ordering::SeqCst);
            }
        }

        // Ensure energy_delivered does not exceed q_released_max
        let max = q_released_max.unwrap_or(q_released_min);

        current_profile.energy_delivered = current_profile.energy_delivered.min(max);

        let new_state_of_charge = self.state_of_charge.load(Ordering::SeqCst)
            + (current_profile.energy_charged - current_profile.energy_delivered)
                / self.storage_capacity;
        let new_state_of_charge = clip(new_state_of_charge, 0., 1.);
        self.state_of_charge
            .store(new_state_of_charge, Ordering::SeqCst);

        // Calculate fan energy
        current_profile.energy_for_fan = 0.;

        if self.air_flow_type == ElectricStorageHeaterAirFlowType::FanAssisted
            && q_released_max.is_some()
        {
            let power_for_fan = self.fan_pwr;
            current_profile.energy_for_fan = Self::convert_to_kwh(power_for_fan, time_used_max);
        }

        // Log the energy charged, fan energy, and total energy delivered
        let amount_demanded = f64::from(self.n_units)
            * (current_profile.energy_charged
                + current_profile.energy_instant
                + current_profile.energy_for_fan);
        self.energy_supply_conn
            .demand_energy(amount_demanded, simulation_time_iteration.index)?;

        // If detailed results flag is set populate with values
        if let Some(esh_detailed_results) = &self.esh_detailed_results {
            let CurrentEnergyProfile {
                energy_for_fan,
                energy_instant,
                energy_charged,
                energy_delivered,
            } = *self.current_energy_profile.read();
            let result = StorageHeaterDetailedResult {
                timestep_idx: simulation_time_iteration.index,
                n_units: self.n_units,
                energy_delivered,
                energy_demand,
                energy_instant,
                energy_charged,
                energy_for_fan,
                state_of_charge: self.state_of_charge.load(Ordering::SeqCst),
                final_soc,
                time_used_max,
            };
            esh_detailed_results
                .write()
                .insert(simulation_time_iteration.index, result);
        }

        // Return total net energy delivered (discharged + instant heat + fan energy)
        Ok(f64::from(self.n_units)
            * (current_profile.energy_delivered + current_profile.energy_instant))
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

        let temp_air: f64 = (self.zone_internal_air_func)();

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
                charge_control.target_charge(simulation_time_iteration, Some(temp_air))
            }
            ControlLogicType::Celect => {
                // Implements the "CELECT" control logic for ESH
                // A CELECT-type controller has electronic sensors throughout the dwelling linked
                // to a central control device. It monitors the individual room sensors and optimises
                // the charging of all the storage heaters individually (and may select direct acting
                // heaters in preference to storage heaters).
                charge_control.target_charge(simulation_time_iteration, Some(temp_air))
            }
            ControlLogicType::Hhrsh => {
                // Implements the "HHRSH" control logic for ESH
                // A ‘high heat retention storage heater’ is one with heat retention not less
                // than 45% measured according to BS EN 60531. It incorporates a timer, electronic
                // room thermostat and fan to control the heat output. It is also able to estimate
                // the next day’s heating demand based on external temperature, room temperature
                // settings and heat demand periods.

                let energy_to_store = charge_control.energy_to_store(
                    self.demand_met.load(Ordering::SeqCst)
                        + self.demand_unmet.load(Ordering::SeqCst),
                    self.zone_setpoint_init,
                    simulation_time_iteration,
                );

                // None means not enough past data to do the calculation (Initial 24h of the calculation)
                // We go for a full load of the hhrsh

                // TODO Python handles a None case from energy_to_store here, which is currently not possible in the Rust
                // if energy_to_store is None:
                // energy_to_store = self.__pwr_in * units.hours_per_day

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

    pub(crate) fn output_esh_results(&self) -> Option<Vec<StorageHeaterDetailedResult>> {
        self.esh_detailed_results
            .as_ref()
            .map(|results| (*results.read()).clone())
    }

    #[cfg(test)]
    fn energy_for_fan(&self) -> f64 {
        self.current_energy_profile.read().energy_for_fan
    }

    #[cfg(test)]
    fn energy_instant(&self) -> f64 {
        self.current_energy_profile.read().energy_instant
    }

    #[cfg(test)]
    fn energy_charged(&self) -> f64 {
        self.current_energy_profile.read().energy_charged
    }

    #[cfg(test)]
    fn energy_delivered(&self) -> f64 {
        self.current_energy_profile.read().energy_delivered
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::excessive_precision)]
    use super::*;
    use crate::{
        core::{
            controls::time_control::{ChargeControl, SetpointTimeControl},
            energy_supply::energy_supply::{EnergySupply, EnergySupplyBuilder},
        },
        external_conditions::{DaylightSavingsConfig, ExternalConditions},
        input::{ControlLogicType, ExternalSensor, FuelType},
        simulation_time::{SimulationTime, SimulationTimeIteration, SimulationTimeIterator},
    };
    use approx::assert_relative_eq;
    use parking_lot::RwLock;
    use rstest::{fixture, rstest};
    use serde_json::json;

    const EIGHT_DECIMAL_PLACES: f64 = 1e-7;

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 24., 1.)
    }

    #[fixture]
    pub fn simulation_time_iterator(simulation_time: SimulationTime) -> SimulationTimeIterator {
        simulation_time.iter()
    }

    #[fixture]
    pub fn simulation_time_iteration(
        simulation_time_iterator: SimulationTimeIterator,
    ) -> SimulationTimeIteration {
        simulation_time_iterator.current_iteration()
    }

    #[fixture]
    fn external_conditions(simulation_time: SimulationTime) -> Arc<ExternalConditions> {
        Arc::new(ExternalConditions::new(
            &simulation_time.iter(),
            vec![
                19.0, 0.0, 1.0, 2.0, 5.0, 7.0, 6.0, 12.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0,
                19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0,
            ],
            vec![
                3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1, 3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1,
                3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1,
            ],
            vec![
                300., 250., 220., 180., 150., 120., 100., 80., 60., 40., 20., 10., 50., 100., 140.,
                190., 200., 320., 330., 340., 350., 355., 315., 5.,
            ],
            vec![
                0., 0., 0., 0., 35., 73., 139., 244., 320., 361., 369., 348., 318., 249., 225.,
                198., 121., 68., 19., 0., 0., 0., 0., 0.,
            ],
            vec![
                0., 0., 0., 0., 0., 0., 7., 53., 63., 164., 339., 242., 315., 577., 385., 285.,
                332., 126., 7., 0., 0., 0., 0., 0.,
            ],
            vec![
                0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
            ],
            51.383,
            -0.783,
            0,
            0,
            Some(0),
            1.,
            Some(1),
            Some(DaylightSavingsConfig::NotApplicable),
            false,
            false,
            serde_json::from_value(json!(
                [
                    {"number": 1, "start360": 180, "end360": 135},
                    {"number": 2, "start360": 135, "end360": 90},
                    {"number": 3, "start360": 90, "end360": 45},
                    {"number": 4, "start360": 45, "end360": 0,
                        "shading": [
                            {"type": "obstacle", "height": 10.5, "distance": 12}
                        ]
                    },
                    {"number": 5, "start360": 0, "end360": -45},
                    {"number": 6, "start360": -45, "end360": -90},
                    {"number": 7, "start360": -90, "end360": -135},
                    {"number": 8, "start360": -135, "end360": -180}
                ]
            ))
            .unwrap(),
        ))
    }

    #[fixture]
    fn external_sensor() -> ExternalSensor {
        serde_json::from_value(json!({
            "correlation": [
                {"temperature": 0.0, "max_charge": 1.0},
                {"temperature": 10.0, "max_charge": 0.9},
                {"temperature": 18.0, "max_charge": 0.5}
            ]
        }))
        .unwrap()
    }

    #[fixture]
    fn charge_control(
        external_conditions: Arc<ExternalConditions>,
        external_sensor: ExternalSensor,
    ) -> Arc<Control> {
        Arc::new(Control::Charge(
            ChargeControl::new(
                ControlLogicType::Automatic,
                vec![
                    true, true, true, true, true, true, true, true, false, false, false, false,
                    false, false, false, false, true, true, true, true, false, false, false, false,
                ],
                1.,
                0,
                1.,
                [1.0, 0.8].into_iter().map(Into::into).collect(),
                Some(22.),
                None,
                None,
                None,
                external_conditions,
                Some(external_sensor),
            )
            .unwrap(),
        ))
    }

    #[fixture]
    fn control() -> Arc<Control> {
        Arc::new(Control::SetpointTime(
            SetpointTimeControl::new(
                vec![Some(21.), Some(21.), None, Some(21.)],
                0,
                1.,
                None,
                None,
                None,
                Default::default(),
                1.,
            )
            .unwrap(),
        ))
    }

    fn create_elec_storage_heater(
        simulation_time: SimulationTime,
        charge_control: Arc<Control>,
        control: Arc<Control>,
        external_conditions: Arc<ExternalConditions>,
        esh_min_output: Vec<(f64, f64)>,
        esh_max_output: Vec<(f64, f64)>,
        output_detailed_results: Option<bool>,
    ) -> ElecStorageHeater {
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
        let energy_supply_conn =
            EnergySupply::connection(energy_supply.clone(), "storage_heater").unwrap();

        let elec_storage_heater = ElecStorageHeater::new(
            3.5,
            2.5,
            10.0,
            ElectricStorageHeaterAirFlowType::FanAssisted,
            0.7,
            11.,
            1,
            21.,
            Arc::new(|| 20.),
            energy_supply_conn,
            &simulation_time.iter(),
            control,
            charge_control,
            esh_min_output,
            esh_max_output,
            external_conditions, // NOTE this is None in Python
            output_detailed_results,
        )
        .unwrap();

        elec_storage_heater
            .state_of_charge
            .store(0.5, Ordering::SeqCst);

        elec_storage_heater
    }

    #[fixture]
    fn elec_storage_heater(
        simulation_time: SimulationTime,
        charge_control: Arc<Control>,
        control: Arc<Control>,
        external_conditions: Arc<ExternalConditions>,
    ) -> ElecStorageHeater {
        let esh_min_output = vec![(0.0, 0.0), (0.5, 0.02), (1.0, 0.05)];
        let esh_max_output = vec![(0.0, 0.0), (0.5, 1.5), (1.0, 3.0)];

        create_elec_storage_heater(
            simulation_time,
            charge_control,
            control,
            external_conditions,
            esh_min_output,
            esh_max_output,
            None,
        )
    }

    #[rstest]
    pub fn test_initialisation(elec_storage_heater: ElecStorageHeater) {
        assert_eq!(elec_storage_heater.pwr_in, 3.5);
        assert_eq!(elec_storage_heater.storage_capacity, 10.0);
        assert_eq!(
            elec_storage_heater.air_flow_type,
            ElectricStorageHeaterAirFlowType::FanAssisted
        );
        assert_eq!(
            elec_storage_heater.state_of_charge.load(Ordering::SeqCst),
            0.5
        );

        assert_relative_eq!(
            elec_storage_heater.heat_retention_ratio,
            0.92372001,
            max_relative = EIGHT_DECIMAL_PLACES
        );
    }

    #[rstest]
    fn test_energy_output_min_single(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: ElecStorageHeater,
    ) {
        let min_energy_output = elec_storage_heater
            .energy_output_min(&simulation_time_iterator.current_iteration())
            .unwrap();
        let _ =
            elec_storage_heater.demand_energy(5.0, &simulation_time_iterator.current_iteration());

        assert_relative_eq!(
            min_energy_output,
            0.019999999999999997,
            max_relative = EIGHT_DECIMAL_PLACES
        );
    }

    #[rstest]
    #[ignore = "known issue"]
    fn test_energy_output_min(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: ElecStorageHeater,
    ) {
        // Test minimum energy output calculation across all timesteps.

        let expected_min_energy_output = [
            0.019999999999999997,
            0.030419151282454364,
            0.0406826518009459,
            0.046014105108884346,
            0.04674323706081289,
            0.045800000000000014,
            0.046400000000000004,
            0.038,
            0.046886897763894056,
            0.03215061095009046,
            0.021233700713726503,
            0.01542628227371725,
            0.011428072222071864,
            0.008466125322130943,
            0.006271860545638567,
            0.004646308882570838,
            0.010432746398675731,
            0.01897277720547578,
            0.019999999999999997,
            0.019999999999999997,
            0.020056174317410178,
            0.014844117518306485,
            0.010996794239857175,
            0.008146626650951819,
        ]; // Actual minimum energy output for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let min_energy_output = elec_storage_heater.energy_output_min(&t_it).unwrap();
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);

            assert_relative_eq!(
                min_energy_output,
                expected_min_energy_output[t_idx],
                max_relative = 0.1
            );
        }
    }

    #[rstest]
    #[ignore = "known issue"]
    fn test_energy_output_max(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: ElecStorageHeater,
    ) {
        // Test maximum energy output calculation across all timesteps.
        let expected_max_energy_output = [
            1.5,
            1.772121660521405,
            2.2199562136927717,
            2.5517202117781994,
            2.7913851590672585,
            2.7899999999999996,
            2.8200000000000003,
            2.4000000000000004,
            2.463423313846487,
            1.8249489529640162,
            1.3519554011630448,
            1.0015529506734968,
            0.7419686857505708,
            0.5496640579327374,
            0.40720123344887615,
            0.30166213975932143,
            0.6996897293886958,
            1.3814284569589004,
            1.5,
            1.5,
            1.3009346098448467,
            0.9637557636923015,
            0.713967931076402,
            0.5289205810615784,
        ]; // Expected max energy output for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let max_energy_output = elec_storage_heater.energy_output_max(&t_it).unwrap();
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);
            let (energy, _, _, _) = max_energy_output;

            assert_relative_eq!(
                energy,
                expected_max_energy_output[t_idx],
                max_relative = 1e-1
            );
        }
    }

    #[rstest]
    fn test_energy_output_max_with_zero_event_single(
        simulation_time: SimulationTime,
        simulation_time_iterator: SimulationTimeIterator,
        charge_control: Arc<Control>,
        control: Arc<Control>,
        external_conditions: Arc<ExternalConditions>,
    ) {
        let esh_min_output = vec![(0.0, 0.0), (0.5, 0.02), (1.0, 0.05)];
        let esh_max_output = vec![(0.0, 0.0), (0.5, 30.0), (1.0, 50.0)];
        let elec_storage_heater = create_elec_storage_heater(
            simulation_time,
            charge_control,
            control,
            external_conditions,
            esh_min_output,
            esh_max_output,
            None,
        );

        let expected_max_energy_output = 7.905696339321716;

        let expected_time_used = 1.0;

        let max_energy_output = elec_storage_heater
            .energy_output_max(&simulation_time_iterator.current_iteration())
            .unwrap();
        let _ =
            elec_storage_heater.demand_energy(5.0, &simulation_time_iterator.current_iteration());
        let (energy, time_used, _, _) = max_energy_output;

        assert_relative_eq!(energy, expected_max_energy_output, max_relative = 1e-2);

        assert_relative_eq!(
            time_used,
            expected_time_used,
            max_relative = EIGHT_DECIMAL_PLACES
        );
    }

    #[rstest]
    #[ignore = "known issue"]
    fn test_energy_output_max_with_zero_event(
        simulation_time: SimulationTime,
        simulation_time_iterator: SimulationTimeIterator,
        charge_control: Arc<Control>,
        control: Arc<Control>,
        external_conditions: Arc<ExternalConditions>,
    ) {
        let esh_min_output = vec![(0.0, 0.0), (0.5, 0.02), (1.0, 0.05)];
        let esh_max_output = vec![(0.0, 0.0), (0.5, 30.0), (1.0, 50.0)];
        let elec_storage_heater = create_elec_storage_heater(
            simulation_time,
            charge_control,
            control,
            external_conditions,
            esh_min_output,
            esh_max_output,
            None,
        );

        // Test maximum energy output calculation across all timesteps.
        let expected_max_energy_output = [
            7.905696339321716,
            6.409423918606012,
            4.913144554873739,
            3.503509067065079,
            3.4999662250517263,
            3.5000272321002064,
            3.4999664471092706,
            3.5000359941220256,
            0.5819017195075272,
            0.0014406990152162622,
            7.97047184775446e-06,
            9.068336550842638e-08,
            0.0,
            0.0,
            0.0,
            0.0,
            2.9181276338320536,
            3.4985510038277714,
            3.5000309931442937,
            3.4999785113707134,
            0.5818631932808774,
            0.0014406043531798149,
            7.969521833559643e-06,
            9.066927928048405e-08,
        ];

        let expected_time_used = [
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            0.37318890798358917,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            0.37318890798725507,
        ];

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let max_energy_output = elec_storage_heater.energy_output_max(&t_it).unwrap();
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);
            let (energy, time_used, _, _) = max_energy_output;

            assert_relative_eq!(
                energy,
                expected_max_energy_output[t_idx],
                max_relative = 1e-2
            );

            assert_relative_eq!(
                time_used,
                expected_time_used[t_idx],
                max_relative = EIGHT_DECIMAL_PLACES
            );
        }
    }

    #[rstest]
    fn test_electric_charge(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: ElecStorageHeater,
    ) {
        // Test electric charge calculation across all timesteps.
        let expected_target_elec_charge = [
            0.5,
            1.0,
            0.99,
            0.98,
            0.95,
            0.93,
            0.9400000000000001,
            0.8,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.5,
            0.5,
            0.5,
            0.5,
            0.0,
            0.0,
            0.0,
            0.0,
        ]; // Expected target charge for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let target_elec_charge = elec_storage_heater.target_electric_charge(t_it).unwrap();
            assert_relative_eq!(target_elec_charge, expected_target_elec_charge[t_idx]);
        }
    }

    #[rstest]
    pub fn test_demand_energy(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: ElecStorageHeater,
    ) {
        let expected_energy = [
            4.0,
            4.272121660521405,
            4.719956213692772,
            5.0,
            5.0,
            5.0,
            5.0,
            4.9,
            4.963423313846487,
            4.324948952964016,
            3.851955401163045,
            3.5015529506734966,
            3.241968685750571,
            3.0496640579327376,
            2.907201233448876,
            2.8016621397593213,
            3.199689729388696,
            3.8814284569589006,
            4.0,
            4.0,
            3.8009346098448464,
            3.4637557636923013,
            3.213967931076402,
            3.0289205810615782,
        ]; // Expected energy for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let energy_out = elec_storage_heater.demand_energy(5.0, &t_it).unwrap();
            assert_relative_eq!(energy_out, expected_energy[t_idx], max_relative = 1e-1);
        }
    }

    #[rstest]
    #[ignore = "known issue"]
    pub fn test_energy_for_fan(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: ElecStorageHeater,
    ) {
        let expected_energy_for_fan = [
            0.003666666666666666,
            0.002707621094790285,
            0.0019580976410457644,
            0.0018707482993197276,
            0.0019298245614035089,
            0.0019713261648745518,
            0.001950354609929078,
            0.0022916666666666662,
            0.0021220628632204496,
            0.002233750362976654,
            0.0023578475867206904,
            0.0024965444862427343,
            0.0026525785014320812,
            0.0028294170564121318,
            0.003031518268419362,
            0.0032647119841350417,
            0.003666666666666666,
            0.003666666666666666,
            0.003666666666666666,
            0.003666666666666666,
            0.0021220628632204505,
            0.0022337503629766544,
            0.0023578475867206913,
            0.002496544486242736,
        ]; // Expected energy for fan for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let _ = elec_storage_heater.demand_energy(0.5, &t_it);
            let energy_for_fan = elec_storage_heater.energy_for_fan();
            assert_relative_eq!(
                energy_for_fan,
                expected_energy_for_fan[t_idx],
                max_relative = EIGHT_DECIMAL_PLACES
            );
        }
    }

    #[rstest]
    #[ignore = "known issue"]
    pub fn test_energy_instant(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: ElecStorageHeater,
    ) {
        let expected_energy_instant = [
            2.5,
            2.5,
            2.5,
            2.4482797882218006,
            2.208614840932741,
            2.2100000000000004,
            2.18,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
            2.5,
        ]; // Expected backup energy instant for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);
            let energy_instant = elec_storage_heater.energy_instant();
            assert_relative_eq!(
                energy_instant,
                expected_energy_instant[t_idx],
                max_relative = 0.1
            );
        }
    }

    #[rstest]
    #[ignore = "known issue"]
    pub fn test_energy_charged(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: ElecStorageHeater,
    ) {
        let expected_energy_charged = [
            1.5,
            3.500000000000001,
            3.5,
            3.5,
            3.3397997646920756,
            2.7899999999999996,
            2.82,
            2.4000000000000004,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            3.500000000000001,
            2.738269398472182,
            1.5,
            1.5,
            0.0,
            0.0,
            0.0,
            0.0,
        ]; // Expected energy charged for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);
            let energy_charged = elec_storage_heater.energy_charged();
            assert_relative_eq!(
                energy_charged,
                expected_energy_charged[t_idx],
                max_relative = 0.1
            );
        }
    }

    #[rstest]
    #[ignore = "known issue"]
    pub fn test_energy_stored_delivered(
        simulation_time_iterator: SimulationTimeIterator,
        elec_storage_heater: ElecStorageHeater,
    ) {
        let expected_energy_delivered = [
            1.5,
            1.772121660521405,
            2.219956213692772,
            2.5517202117781994,
            2.791385159067259,
            2.7899999999999996,
            2.82,
            2.4000000000000004,
            2.4634233138464876,
            1.8249489529640166,
            1.3519554011630448,
            1.0015529506734968,
            0.7419686857505706,
            0.5496640579327375,
            0.4072012334488761,
            0.30166213975932155,
            0.6996897293886956,
            1.3814284569589008,
            1.5,
            1.5,
            1.300934609844847,
            0.9637557636923016,
            0.713967931076402,
            0.5289205810615786,
        ]; // Expected energy stored delivered for each timestep

        for (t_idx, t_it) in simulation_time_iterator.enumerate() {
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);
            let energy_delivered = elec_storage_heater.energy_delivered();
            assert_relative_eq!(
                energy_delivered,
                expected_energy_delivered[t_idx],
                max_relative = 1e-2
            );
        }
    }

    #[test]
    pub fn test_heat_retention_output() {
        let soc_array = vec![0., 0.5, 1.];
        let power_array = vec![0., 0.02, 0.05];
        let storage_capacity = 10.;
        let actual =
            ElecStorageHeater::heat_retention_output(&soc_array, &power_array, storage_capacity);

        assert_relative_eq!(actual, 0.92372001, max_relative = EIGHT_DECIMAL_PLACES);
    }

    #[rstest]
    fn test_initialisation_invalid_soc_arrays(
        simulation_time: SimulationTime,
        charge_control: Arc<Control>,
        control: Arc<Control>,
        external_conditions: Arc<ExternalConditions>,
    ) {
        let test_cases = vec![
            (
                vec![(0.0, 0.0), (0.5, 0.02), (0.3, 0.02), (1.0, 0.05)],
                vec![(0.0, 0.0), (0.5, 1.5), (0.7, 0.02), (1.0, 3.)],
                "shouldn't allow esh_min_output values in non-increasing order",
            ),
            (
                vec![(0.0, 0.0), (0.5, 0.02), (0.7, 0.02), (1.0, 0.05)],
                vec![(0.0, 0.0), (0.5, 1.5), (0.3, 0.02), (1.0, 3.)],
                "shouldn't allow esh_max_output values in non-increasing order",
            ),
            (
                vec![(0.0, 0.0), (0.5, 0.02), (0.7, 0.02), (0.9, 0.05)],
                vec![(0.0, 0.0), (0.5, 1.5), (0.7, 0.02), (1.0, 3.)],
                "shouldn't allow esh_min_output values not ending in 1.0",
            ),
            (
                vec![(0.0, 0.0), (0.5, 0.02), (0.7, 0.02), (1.0, 0.05)],
                vec![(0.0, 0.0), (0.5, 1.5), (0.7, 0.02), (0.9, 3.)],
                "shouldn't allow esh_max_output values not ending in 1.0",
            ),
            (
                vec![(0.2, 0.0), (0.5, 0.02), (0.7, 0.02), (1.0, 0.05)],
                vec![(0.0, 0.0), (0.5, 1.5), (0.7, 0.02), (1.0, 3.0)],
                "shouldn't allow esh_min_output values not starting at 1.0",
            ),
            (
                vec![(0.0, 0.0), (0.5, 0.02), (0.7, 0.02), (1.0, 0.05)],
                vec![(0.2, 0.0), (0.5, 1.5), (0.7, 0.02), (1.0, 3.0)],
                "shouldn't allow esh_max_output values not starting at 1.0",
            ),
            (
                vec![(0.0, 0.0), (1.0, 1.0)],
                vec![(0.0, 0.0), (1.0, 0.5)],
                "shouldn't allow any power_max values below power_min",
            ),
            (
                vec![(0.0, 0.0), (1.0, 1.0)],
                vec![(0.0, 0.0), (0.5, 0.4), (1.0, 1.0)],
                "shouldn't allow any power_max values below power_min",
            ),
        ];

        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
        let energy_supply_conn =
            EnergySupply::connection(energy_supply.clone(), "storage_heater").unwrap();

        for (esh_min_output, esh_max_output, debug_msg) in test_cases.iter() {
            assert!(
                ElecStorageHeater::new(
                    1.,
                    1.,
                    10.,
                    ElectricStorageHeaterAirFlowType::FanAssisted,
                    1.,
                    10.,
                    1,
                    21.,
                    Arc::new(|| 20.),
                    energy_supply_conn.clone(),
                    &simulation_time.iter(),
                    control.clone(),
                    charge_control.clone(),
                    esh_min_output.clone(),
                    esh_max_output.clone(),
                    external_conditions.clone(), // NOTE this is None in Python
                    None,
                )
                .is_err(),
                "{}",
                debug_msg
            );
        }
    }

    fn test_initialisation_detailed_results(
        simulation_time: SimulationTime,
        simulation_time_iterator: SimulationTimeIterator,
        charge_control: Arc<Control>,
        control: Arc<Control>,
        external_conditions: Arc<ExternalConditions>,
    ) {
        let esh_min_output = vec![(0.0, 0.0), (0.5, 0.02), (1.0, 0.05)];
        let esh_max_output = vec![(0.0, 0.0), (0.5, 30.0), (1.0, 50.0)];
        let heater_with_detailed_results = create_elec_storage_heater(
            simulation_time,
            charge_control.clone(),
            control.clone(),
            external_conditions.clone(),
            esh_min_output.clone(),
            esh_max_output.clone(),
            Some(true),
        );
        let heater_without_detailed_results = create_elec_storage_heater(
            simulation_time,
            charge_control,
            control,
            external_conditions,
            esh_min_output,
            esh_max_output,
            None,
        );

        // TODO add assertions
    }
}
