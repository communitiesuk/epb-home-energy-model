use std::sync::Arc;

use anyhow::bail;
use derivative::Derivative;
use itertools::Itertools;
use nalgebra::{Vector1, Vector3};
use ode_solvers::{dop_shared::OutputType, Dopri5, System};

use crate::{
    core::{
        controls::time_control::{ChargeControl, Control, SetpointTimeControl},
        energy_supply::energy_supply::EnergySupplyConnection,
    },
    external_conditions::ExternalConditions,
    input::ElectricStorageHeaterAirFlowType,
    simulation_time::{SimulationTime, SimulationTimeIteration},
    statistics::np_interp,
};

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
pub struct ElecStorageHeater {
    pwr_in: f64,
    pwr_instant: f64,
    storage_capacity: f64,
    air_flow_type: ElectricStorageHeaterAirFlowType,
    frac_convective: f64,
    n_units: i32,
    energy_supply_conn: EnergySupplyConnection,
    simulation_time: SimulationTime,
    control: Arc<Control>,
    charge_control: Arc<Control>,
    fan_pwr: f64,
    external_conditions: ExternalConditions,
    temp_air: f64,
    state_of_charge: f64,
    esh_min_output: Vec<(f64, f64)>,
    esh_max_output: Vec<(f64, f64)>,
    demand_met: f64,
    demand_unmet: f64,
    zone_setpoint_init: f64,
    #[derivative(Debug = "ignore")]
    zone_internal_air_func: Arc<dyn Fn() -> f64>,
    soc_max_array: Vec<f64>,
    power_max_array: Vec<f64>,
    soc_min_array: Vec<f64>,
    power_min_array: Vec<f64>,
    heat_retention_ratio: f64,
    // TODO review - do we need to keep these as public properties?
    pub energy_for_fan: f64,
    pub energy_instant: f64,
    pub energy_charged: f64,
    pub energy_delivered: f64,
}

struct SocOdeFunction {
    soc_array: Vec<f64>,
    power_array: Vec<f64>,
    storage_capacity: f64,
}

impl System<Time, State> for SocOdeFunction {
    fn system(&self, _x: Time, y: &State, dy: &mut State) {
        // Define the ODE for SOC and energy delivered (no charging, only discharging)

        // Ensure SOC stays within bounds
        let soc = clip(y[0], 0., 1.);

        // Discharging: calculate power used based on SOC
        let discharge_rate = -np_interp(soc, &self.soc_array, &self.power_array);

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
        let discharge_rate = -np_interp(soc, &self.soc_array, &self.power_array);

        // Track the total energy delivered (discharged energy)
        let ddelivered_dt = -discharge_rate; // Energy delivered (positive value)

        let dcharged_dt = if soc > self.soc_max {
            self.charge_rate
        } else {
            if self.target_charge > 0. {
                ddelivered_dt.min(self.pwr_in)
            } else {
                0.0
            }
        };

        // Net SOC rate of change (discharge + charge), divided by storage capacity
        let dsoc_dt = (-ddelivered_dt + dcharged_dt) / self.storage_capacity;

        dy[0] = dsoc_dt;
        dy[1] = dcharged_dt;
        dy[2] = ddelivered_dt;
    }

    fn solout(&mut self, _x: Time, y: &EnergyOutputState, _dy: &EnergyOutputState) -> bool {
        // TODO we want to check if we've passed this value, not that we are equal to it
        // see Emitters for example of this
        y[0] == self.soc_max
    }
}

pub(crate) enum OutputMode {
    Min,
    Max,
}

impl ElecStorageHeater {
    pub fn new(
        pwr_in: f64,
        rated_power_instant: f64,
        storage_capacity: f64,
        air_flow_type: ElectricStorageHeaterAirFlowType,
        frac_convective: f64,
        fan_pwr: f64,
        n_units: i32,
        zone_setpoint_init: f64,
        zone_internal_air_func: Arc<dyn Fn() -> f64>,
        energy_supply_conn: EnergySupplyConnection,
        simulation_time: SimulationTime,
        control: Arc<Control>,
        charge_control: Arc<Control>,
        esh_min_output: Vec<(f64, f64)>,
        esh_max_output: Vec<(f64, f64)>,
        ext_cond: ExternalConditions,
    ) -> anyhow::Result<Self> {
        // Arguments:
        // pwr_in               -- in kW (Charging)
        // rated_power_instant  -- in kW (Instant backup)
        // storage_capacity     -- in kWh
        // air_flow_type        -- string specifying type of Electric Storage Heater:
        //                      -- fan-assisted
        //                      -- damper-only
        // frac_convective      -- convective fraction for heating (TODO: Check if necessary)
        // fan_pwr              -- Fan power [W]
        // n_units              -- number of units install in zone
        // zone                 -- zone where the unit(s) is/are installed
        // energy_supply_conn   -- reference to EnergySupplyConnection object
        // simulation_time      -- reference to SimulationTime object
        // control              -- reference to a control object which must implement is_on() and setpnt() funcs
        // charge_control       -- reference to a ChargeControl object which must implement different logic types
        //                         for charging the Electric Storage Heaters.
        // ESH_min_output       -- Data from test showing the output from the storage heater when not actively
        //                         outputting heat, i.e. case losses only (with units kW)
        // ESH_max_output       -- Data from test showing the output from the storage heater when it is actively
        //                         outputting heat, e.g. damper open / fan running (with units kW)
        // extcond              -- reference to ExternalConditions object

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

        // Validate that both SOC arrays start at 0.0 and end at 1.0
        // TODO Result or more specific panic
        if !is_close!(*soc_max_array.first().unwrap(), 0.) {
            panic!("The first SOC value in esh_max_output must be 0.0 (fully discharged).");
        }

        if !is_close!(*soc_max_array.last().unwrap(), 1.) {
            panic!("The last SOC value in esh_max_output must be 1.0 (fully charged).");
        }

        if !is_close!(*soc_min_array.first().unwrap(), 0.) {
            panic!("The first SOC value in esh_min_output must be 0.0 (fully discharged).");
        }

        if !is_close!(*soc_min_array.last().unwrap(), 1.) {
            panic!("The last SOC value in esh_min_output must be 1.0 (fully charged).");
        }

        // Validate that for any SOC, power_max >= power_min
        // Sample a fine grid of SOCs and ensure power_max >= power_min
        let fine_soc: Vec<f64> = linspace(0., 1., 100);

        let power_max_fine: Vec<f64> = fine_soc
            .iter()
            .map(|s| np_interp(*s, &soc_max_array, &power_max_array))
            .collect();
        let power_min_fine: Vec<f64> = fine_soc
            .iter()
            .map(|s| np_interp(*s, &soc_min_array, &power_min_array))
            .collect();

        for i in 0..fine_soc.len() {
            if power_max_fine[i] < power_min_fine[i] {
                panic!("At all SOC levels, ESH_max_output must be >= ESH_min_output.")
            }
        }

        // TODO can we pass these vecs/arrays by reference instead
        let heat_retention_ratio = Self::heat_retention_output(
            soc_min_array.clone(),
            power_min_array.clone(),
            storage_capacity,
        );

        Ok(Self {
            pwr_in,
            pwr_instant: rated_power_instant,
            storage_capacity,
            air_flow_type,
            frac_convective,
            n_units,
            energy_supply_conn,
            simulation_time,
            control,
            charge_control,
            fan_pwr,
            external_conditions: ext_cond,
            temp_air,
            state_of_charge: 0.,
            esh_min_output,
            esh_max_output,
            demand_met: 0.,
            demand_unmet: 0.,
            zone_setpoint_init,
            zone_internal_air_func,
            soc_max_array,
            power_max_array,
            soc_min_array,
            power_min_array,
            heat_retention_ratio,
            // TODO ...
            energy_for_fan: 0.,
            energy_instant: 0.,
            energy_charged: 0.,
            energy_delivered: 0.,
        })
    }

    pub fn heat_retention_output(
        soc_array: Vec<f64>,
        power_array: Vec<f64>,
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
            soc_array: soc_array,
            power_array: power_array,
            storage_capacity: storage_capacity,
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
        let final_soc = clip(final_soc, 0., 1.);

        // Return the final state of charge after 16 hours
        return final_soc;
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
        let current_soc = self.state_of_charge;
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

        // TODO handle events

        // scipy implementation for reference:
        // https://github.com/scipy/scipy/blob/6b657ede0c3c4cffef3156229afddf02a2b1d99a/scipy/integrate/_ivp/rk.py#L293
        let rtol = 1e-4; // rtol - set to match PythonPython - Relative tolerance used in the computation of the adaptive step size
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
        let time_used = 0.; // TODO implement with root solver

        // Return the total energy delivered, time used, and total energy charged
        Ok((
            total_energy_delivered,
            time_used,
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

    pub fn energy_output_max(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64, f64, f64)> {
        // Calculates the maximum energy that can be delivered based on ESH_max_output.
        // :return: Tuple containing (maximum energy deliverable in kWh, time used in hours).
        self.energy_output(OutputMode::Max, simulation_time_iteration)
    }

    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        
        // Determines the amount of energy to release based on energy demand, while also handling the
        // energy charging and logging fan energy.
        // :param energy_demand: Energy demand in kWh.
        // :return: Total net energy delivered (including instant heating and fan energy).
        
        let timestep = simulation_time_iteration.timestep;
        let energy_demand = energy_demand / f64::from(self.n_units);
        self.energy_instant = 0.;
    
        // Initialize time_used_max and energy_charged_max to default values
        let time_used_max = 0.;
        let energy_charged_max = 0.;
    
        // Calculate minimum energy that can be delivered
        let (q_released_min, _, energy_charged, final_soc) = self.energy_output(OutputMode::Min, simulation_time_iteration)?;
        self.energy_charged = energy_charged;

        // TODO double check nesting matches Python

        let mut time_instant: f64;
        let mut time_used_max: f64;

        if q_released_min > energy_demand {
            // Deliver at least the minimum energy
            self.energy_delivered = q_released_min;
            self.demand_met = q_released_min;
            self.demand_unmet = 0.;
        }
        else {
            // Calculate maximum energy that can be delivered
            let (q_released_max, time_used_max_tmp, energy_charged, final_soc) = self.energy_output(OutputMode::Max, simulation_time_iteration)?;
            time_used_max = time_used_max_tmp;
            self.energy_charged = energy_charged;
            
            if q_released_max < energy_demand {
                // Deliver as much as possible up to the maximum energy
                self.energy_delivered = q_released_max;
                self.demand_met = q_released_max;
                self.demand_unmet = energy_demand - q_released_max;

                // For now, we assume demand not met from storage is topped-up by
                // the direct top-up heater (if applicable). If still some unmet, 
                // this is reported as unmet demand.
                // if self.pwr_instant {

                self.energy_instant = self.demand_unmet.min(self.pwr_instant * f64::from(timestep)); // kWh
                time_instant = self.energy_instant / self.pwr_instant;
                time_used_max += time_instant;
                time_used_max = time_used_max.min(timestep);
                
                //}
            }
                    
            else {
                // Deliver the demanded energy
                self.energy_delivered = energy_demand;

                if q_released_max > 0. {
                    time_used_max *= energy_demand / q_released_max;
                }

                self.demand_met = energy_demand;
                self.demand_unmet = 0.;
            }
        }
    
        todo!();
        // Ensure energy_delivered does not exceed q_released_max
        // self.energy_delivered = min(self.energy_delivered, q_released_max if 'q_released_max' in locals() else q_released_min)

        // ...
    }

    pub fn target_electric_charge(
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

        let target_charge = match logic_type {
            crate::input::ControlLogicType::Manual => {
                // Implements the "Manual" control logic for ESH
                charge_control.target_charge(simulation_time_iteration, None)
            }
            crate::input::ControlLogicType::Automatic => {
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
            crate::input::ControlLogicType::Celect => {
                // Implements the "CELECT" control logic for ESH
                // A CELECT-type controller has electronic sensors throughout the dwelling linked
                // to a central control device. It monitors the individual room sensors and optimises
                // the charging of all the storage heaters individually (and may select direct acting
                // heaters in preference to storage heaters).
                charge_control.target_charge(simulation_time_iteration, Some(temp_air))
            }
            crate::input::ControlLogicType::Hhrsh => {
                // Implements the "HHRSH" control logic for ESH
                // A ‘high heat retention storage heater’ is one with heat retention not less
                // than 45% measured according to BS EN 60531. It incorporates a timer, electronic
                // room thermostat and fan to control the heat output. It is also able to estimate
                // the next day’s heating demand based on external temperature, room temperature
                // settings and heat demand periods.

                let energy_to_store = charge_control.energy_to_store(
                    self.demand_met + self.demand_unmet,
                    self.zone_setpoint_init,
                    simulation_time_iteration,
                );

                // None means not enough past data to do the calculation (Initial 24h of the calculation)
                // We go for a full load of the hhrsh

                // TODO Python handles a None case from energy_to_store here, which is currently not possible in the Rust
                // if energy_to_store is None:
                // energy_to_store = self.__pwr_in * units.hours_per_day

                let target_charge_hhrsh = if energy_to_store > 0. {
                    let energy_stored = self.state_of_charge * self.storage_capacity; // kWh

                    let energy_to_add = if self.heat_retention_ratio <= 0. {
                        self.storage_capacity - energy_stored // kWh
                    } else {
                        (1.0 / self.heat_retention_ratio) * (energy_to_store - energy_stored)
                        // kWh
                    };
                    let target_charge_hhrsh =
                        self.state_of_charge + energy_to_add / self.storage_capacity;
                    clip(target_charge_hhrsh, 0., 1.)
                } else {
                    0.
                };

                // target_charge (from input file, or zero when control is off) applied here
                // is treated as an upper limit for target charge
                charge_control
                    .target_charge(simulation_time_iteration, None)
                    .and_then(|tc| Ok(tc.min(target_charge_hhrsh)))
            }
        };

        target_charge
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        core::{
            controls::time_control::SetpointTimeControl,
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
        SimulationTime::new(0., 2., 1.)
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
    fn external_conditions(simulation_time: SimulationTime) -> ExternalConditions {
        ExternalConditions::new(
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
            // TODO check these are correct
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
        )
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
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
    ) -> Control {
        Control::Charge(
            ChargeControl::new(
                ControlLogicType::Automatic,
                vec![
                    true, true, true, true, true, true, true, true, false, false, false, false,
                    false, false, false, false, true, true, true, true, false, false, false, false,
                ],
                1.,
                0,
                1.,
                vec![1.0, 0.8],
                Some(22.),
                None,
                None,
                None,
                external_conditions.into(),
                Some(external_sensor),
            )
            .unwrap(),
        )
    }

    #[fixture]
    fn control() -> Control {
        Control::SetpointTime(
            SetpointTimeControl::new(
                vec![Some(21.), Some(21.), None, Some(21.)],
                0,
                1.,
                None,
                None,
                None,
                None,
                1.,
            )
            .unwrap(),
        )
    }

    #[fixture]
    fn elec_storage_heater(
        simulation_time: SimulationTime,
        charge_control: Control,
        control: Control,
        external_conditions: ExternalConditions,
    ) -> ElecStorageHeater {
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
        let energy_supply_conn =
            EnergySupply::connection(energy_supply.clone(), "storage_heater").unwrap();

        let esh_min_output = vec![(0.0, 0.0), (0.5, 0.02), (1.0, 0.05)];
        let esh_max_output = vec![(0.0, 0.0), (0.5, 1.5), (1.0, 3.0)];

        ElecStorageHeater::new(
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
            simulation_time,
            Arc::new(control),
            Arc::new(charge_control),
            esh_min_output,
            esh_max_output,
            external_conditions, // TODO this is None in Python
        )
        .unwrap()
    }

    #[rstest]
    pub fn test_initialisation(elec_storage_heater: ElecStorageHeater) {
        assert_eq!(elec_storage_heater.pwr_in, 3.5);
        assert_eq!(elec_storage_heater.storage_capacity, 10.0);
        assert_eq!(
            elec_storage_heater.air_flow_type,
            ElectricStorageHeaterAirFlowType::FanAssisted
        );
        assert_relative_eq!(
            elec_storage_heater.heat_retention_ratio,
            0.92372001,
            max_relative = EIGHT_DECIMAL_PLACES
        );
    }

    #[rstest]
    #[ignore = "not yet implemented"]
    fn test_energy_output_min(
        simulation_time_iterator: SimulationTimeIterator,
        mut elec_storage_heater: ElecStorageHeater,
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

            // TODO is this line needed?
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);

            assert_relative_eq!(min_energy_output, expected_min_energy_output[t_idx]);
        }
        assert!(false);
    }

    #[rstest]
    #[ignore = "not yet implemented"]
    fn test_energy_output_max(
        simulation_time_iterator: SimulationTimeIterator,
        mut elec_storage_heater: ElecStorageHeater,
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

            // TODO is this line needed
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);
            let (energy, _, _, _) = max_energy_output;

            assert_relative_eq!(energy, expected_max_energy_output[t_idx]);
        }
    }

    #[rstest]
    #[ignore = "not yet implemented"]
    fn test_electric_charge(
        simulation_time_iterator: SimulationTimeIterator,
        mut elec_storage_heater: ElecStorageHeater,
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
    #[ignore = "not yet implemented"]
    pub fn test_demand_energy(
        simulation_time_iterator: SimulationTimeIterator,
        mut elec_storage_heater: ElecStorageHeater,
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
            assert_relative_eq!(energy_out, expected_energy[t_idx]);
        }
    }

    #[rstest]
    #[ignore = "not yet implemented"]
    pub fn test_energy_for_fan(
        simulation_time_iterator: SimulationTimeIterator,
        mut elec_storage_heater: ElecStorageHeater,
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
            let _ = elec_storage_heater.demand_energy(5.0, &t_it);
            let energy_for_fan = elec_storage_heater.energy_for_fan;
            assert_relative_eq!(energy_for_fan, expected_energy_for_fan[t_idx]);
        }
    }

    #[rstest]
    #[ignore = "not yet implemented"]
    pub fn test_energy_instant(
        simulation_time_iterator: SimulationTimeIterator,
        mut elec_storage_heater: ElecStorageHeater,
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
            let energy_instant = elec_storage_heater.energy_instant;
            assert_relative_eq!(energy_instant, expected_energy_instant[t_idx]);
        }
    }

    #[rstest]
    #[ignore = "not yet implemented"]
    pub fn test_energy_charged(
        simulation_time_iterator: SimulationTimeIterator,
        mut elec_storage_heater: ElecStorageHeater,
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
            let energy_charged = elec_storage_heater.energy_charged;
            assert_relative_eq!(energy_charged, expected_energy_charged[t_idx]);
        }
    }

    #[rstest]
    #[ignore = "not yet implemented"]
    pub fn test_energy_stored_delivered(
        simulation_time_iterator: SimulationTimeIterator,
        mut elec_storage_heater: ElecStorageHeater,
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
            let energy_delivered = elec_storage_heater.energy_delivered;
            assert_relative_eq!(energy_delivered, expected_energy_delivered[t_idx]);
        }
    }

    #[test]
    pub fn test_heat_retention_output() {
        let soc_array = vec![0., 0.5, 1.];
        let power_array = vec![0., 0.02, 0.05];
        let storage_capacity = 10.;
        let actual =
            ElecStorageHeater::heat_retention_output(soc_array, power_array, storage_capacity);

        assert_relative_eq!(actual, 0.92372001, max_relative = EIGHT_DECIMAL_PLACES);
    }
}
