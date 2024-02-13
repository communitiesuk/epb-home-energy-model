use crate::compare_floats::min_of_2;
use crate::core::controls::time_control::{per_control, Control, ControlBehaviour};
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::external_conditions::ExternalConditions;
use crate::input::{HeatSourceLocation, HeatSourceWetDetails};
use crate::simulation_time::{SimulationTimeIteration, SimulationTimeIterator};
use interp::interp;
use std::sync::Arc;

/// This module provides object(s) to model the behaviour of heat batteries.

enum ServiceType {
    WaterRegular,
    Space,
}

/// An object to represent a water heating service provided by a regular heat battery.
///
/// This object contains the parts of the heat battery calculation that are
/// specific to providing hot water.
pub struct HeatBatteryServiceWaterRegular {
    heat_battery: HeatBattery,
    service_name: String,
    control: Option<Arc<Control>>,
    temp_hot_water: f64,
    cold_feed: Arc<ColdWaterSource>,
    temp_return: f64,
}

impl HeatBatteryServiceWaterRegular {
    pub fn new(
        heat_battery: HeatBattery,
        service_name: String,
        temp_hot_water: f64,
        cold_feed: Arc<ColdWaterSource>,
        temp_return: f64,
        control: Option<Arc<Control>>,
    ) -> Self {
        Self {
            heat_battery,
            service_name,
            control,
            temp_hot_water,
            cold_feed,
            temp_return,
        }
    }

    /// Demand energy (in kWh) from the heat_battery
    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        let service_on = self.is_on(simulation_time_iteration);
        let energy_demand = if !service_on { 0.0 } else { energy_demand };

        self.heat_battery.demand_energy(
            &self.service_name,
            ServiceType::WaterRegular,
            energy_demand,
            self.temp_return,
        )
    }

    pub fn energy_output_max(&mut self, simulation_time_iteration: SimulationTimeIteration) -> f64 {
        let service_on = self.is_on(simulation_time_iteration);
        if !service_on {
            return 0.0;
        }

        self.heat_battery.energy_output_max(self.temp_hot_water)
    }

    fn is_on(&self, simulation_time_iteration: SimulationTimeIteration) -> bool {
        if self.control.is_some() {
            per_control!(self.control.as_ref().unwrap().as_ref(), ctrl => { ctrl.is_on(&simulation_time_iteration) })
        } else {
            true
        }
    }
}

pub struct HeatBatteryServiceSpace {
    heat_battery: HeatBattery,
    service_name: String,
    control: Arc<Control>,
}

/// An object to represent a space heating service provided by a heat_battery to e.g. radiators.
///
/// This object contains the parts of the heat battery calculation that are
/// specific to providing space heating.
impl HeatBatteryServiceSpace {
    pub fn new(heat_battery: HeatBattery, service_name: String, control: Arc<Control>) -> Self {
        Self {
            heat_battery,
            service_name,
            control,
        }
    }

    pub fn temp_setpnt(&self, simulation_time_iteration: SimulationTimeIteration) -> Option<f64> {
        per_control!(self.control.as_ref(), ctrl => { ctrl.setpnt(&simulation_time_iteration) })
    }

    pub fn in_required_period(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> Option<bool> {
        per_control!(self.control.as_ref(), ctrl => { ctrl.in_required_period(&simulation_time_iteration) })
    }

    /// Demand energy (in kWh) from the heat battery
    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        _temp_flow: f64,
        temp_return: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        if !self.is_on(simulation_time_iteration) {
            return 0.0;
        }

        self.heat_battery.demand_energy(
            &self.service_name,
            ServiceType::Space,
            energy_demand,
            temp_return,
        )
    }

    fn is_on(&self, simulation_time_iteration: SimulationTimeIteration) -> bool {
        per_control!(self.control.as_ref(), ctrl => { ctrl.is_on(&simulation_time_iteration) })
    }
}

const LABS_TESTS_RATED_OUTPUT: [[f64; 2]; 21] = [
    [0.0, 0.0],
    [0.08, 0.00],
    [0.16, 0.03],
    [0.17, 0.05],
    [0.19, 0.10],
    [0.21, 0.15],
    [0.23, 0.21],
    [0.25, 0.23],
    [0.28, 0.26],
    [0.31, 0.29],
    [0.34, 0.32],
    [0.38, 0.36],
    [0.42, 0.41],
    [0.47, 0.45],
    [0.52, 0.51],
    [0.58, 0.57],
    [0.64, 0.64],
    [0.72, 0.71],
    [0.8, 0.8],
    [0.89, 0.89],
    [1.0, 1.0],
];

const LABS_TESTS_RATED_OUTPUT_ENHANCED: [[f64; 2]; 16] = [
    [0.0, 0.0],
    [0.101, 0.0],
    [0.12, 0.18],
    [0.144, 0.235],
    [0.175, 0.313],
    [0.215, 0.391],
    [0.266, 0.486],
    [0.328, 0.607],
    [0.406, 0.728],
    [0.494, 0.795],
    [0.587, 0.825],
    [0.683, 0.875],
    [0.781, 0.906],
    [0.891, 0.953],
    [0.981, 0.992],
    [1.0, 1.0],
];

const LABS_TESTS_LOSSES: [[f64; 2]; 20] = [
    [0.0, 0.0],
    [0.16, 0.13],
    [0.17, 0.15],
    [0.19, 0.17],
    [0.21, 0.18],
    [0.23, 0.21],
    [0.25, 0.23],
    [0.28, 0.26],
    [0.31, 0.29],
    [0.34, 0.32],
    [0.38, 0.36],
    [0.42, 0.41],
    [0.47, 0.45],
    [0.52, 0.51],
    [0.58, 0.57],
    [0.64, 0.64],
    [0.72, 0.71],
    [0.8, 0.8],
    [0.89, 0.89],
    [1.0, 1.0],
];

const HEAT_BATTERY_TIME_UNIT: u32 = 3_600;

struct HeatBatteryResult {
    service_name: String,
    time_running: f64,
    current_hb_power: f64,
}

pub struct HeatBattery {
    simulation_time: Arc<SimulationTimeIterator>,
    external_conditions: Arc<ExternalConditions>,
    heat_battery_location: HeatSourceLocation,
    pwr_in: f64,
    heat_storage_capacity: f64,
    max_rated_heat_output: f64,
    max_rated_losses: f64,
    power_circ_pump: f64,
    power_standby: f64,
    n_units: usize,
    charge_control: Arc<Control>, // ToUChargeControl variant expected
    service_results: Vec<HeatBatteryResult>,
    total_time_running_current_timestamp: f64,
    flag_first_call: bool,
    charge_level: f64,
    q_in_ts: Option<f64>,
    q_out_ts: Option<f64>,
    q_loss_ts: Option<f64>,
}

impl HeatBattery {
    pub fn new(
        heat_battery_details: &HeatSourceWetDetails,
        charge_control: Arc<Control>,
        simulation_time: Arc<SimulationTimeIterator>,
        external_conditions: Arc<ExternalConditions>,
    ) -> Self {
        let (
            heat_battery_location,
            pwr_in,
            heat_storage_capacity,
            max_rated_heat_output,
            max_rated_losses,
            power_circ_pump,
            power_standby,
            n_units,
        ) = if let HeatSourceWetDetails::HeatBattery {
            heat_battery_location,
            rated_charge_power: pwr_in,
            heat_storage_capacity,
            max_rated_heat_output,
            max_rated_losses,
            electricity_circ_pump: power_circ_pump,
            electricity_standby: power_standby,
            number_of_units: n_units,
            ..
        } = heat_battery_details
        {
            (
                *heat_battery_location,
                *pwr_in,
                *heat_storage_capacity,
                *max_rated_heat_output,
                *max_rated_losses,
                *power_circ_pump,
                *power_standby,
                *n_units,
            )
        } else {
            unreachable!()
        };

        Self {
            simulation_time,
            external_conditions,
            heat_battery_location,
            pwr_in,
            heat_storage_capacity,
            max_rated_heat_output,
            max_rated_losses,
            power_circ_pump,
            power_standby,
            n_units,
            charge_control,
            service_results: Default::default(),
            total_time_running_current_timestamp: Default::default(),
            flag_first_call: true,
            charge_level: Default::default(),
            q_in_ts: Default::default(),
            q_out_ts: Default::default(),
            q_loss_ts: Default::default(),
        }
    }

    /// Converts power value supplied to the correct units
    ///        
    /// Arguments:
    /// * `power` - Power in watts
    /// * `timestep` - length of the timestep
    ///
    /// returns  -- Energy in kWH
    fn convert_to_energy(&self, power: f64, timestep: f64) -> f64 {
        power * timestep * self.n_units as f64
    }

    /// Calculates power required for unit
    ///        
    /// Arguments
    /// * `time` - current time period that we are looking at
    /// * `simulation_time_iteration` - an iteration of the contextual simulation time
    ///
    /// returns -- Power required in watts
    fn electric_charge(&self, _time: f64) -> f64 {
        if per_control!(self.charge_control.as_ref(), ctrl => { ctrl.is_on(&self.simulation_time.current_iteration()) })
        {
            self.pwr_in
        } else {
            0.0
        }
    }

    fn lab_test_rated_output(&self, charge_level: f64) -> f64 {
        // labs_test for heat battery
        let x = LABS_TESTS_RATED_OUTPUT_ENHANCED
            .iter()
            .map(|row| row[0])
            .collect::<Vec<_>>();
        let y = LABS_TESTS_RATED_OUTPUT_ENHANCED
            .iter()
            .map(|row| row[1])
            .collect::<Vec<_>>();
        interp(&x, &y, charge_level) * self.max_rated_heat_output
    }

    fn lab_test_losses(&self, charge_level: f64) -> f64 {
        let x = LABS_TESTS_LOSSES
            .iter()
            .map(|row| row[0])
            .collect::<Vec<_>>();
        let y = LABS_TESTS_LOSSES
            .iter()
            .map(|row| row[1])
            .collect::<Vec<_>>();
        interp(&x, &y, charge_level) * self.max_rated_losses
    }

    fn first_call(&mut self) {
        let timestep = self.simulation_time.step_in_hours();
        let current_hour = self.simulation_time.current_hour();
        let time_range = current_hour * HEAT_BATTERY_TIME_UNIT;
        let charge_level = self.charge_level;
        let mut charge_level_qin = charge_level;
        let target_charge = self.target_charge();

        self.q_in_ts = Some(self.electric_charge(time_range as f64));

        // Calculate max charge level possible in next timestep
        if charge_level_qin < target_charge {
            let delta_charge_level = self.q_in_ts.unwrap() * timestep / self.heat_storage_capacity;
            charge_level_qin += delta_charge_level;
            if charge_level_qin > target_charge {
                charge_level_qin = target_charge;
            }
        }

        // Estimating output rate at average of capacity in timestep
        let max_output = self.lab_test_rated_output(charge_level_qin);
        let delta_charge_level = max_output * timestep / self.heat_storage_capacity;
        self.q_out_ts =
            Some(self.lab_test_rated_output(charge_level_qin - delta_charge_level / 2.));
        self.q_loss_ts = Some(self.lab_test_losses(charge_level_qin - delta_charge_level / 2.));
        self.flag_first_call = false;
    }

    pub fn demand_energy(
        &mut self,
        service_name: &str,
        _service_type: ServiceType,
        energy_output_required: f64,
        _temp_return_feed: f64,
    ) -> f64 {
        let timestep = self.simulation_time.step_in_hours();
        let mut charge_level = self.charge_level;

        // Picking target charge level from control
        let target_charge = self.target_charge();
        // __demand_energy is called for each service in each timestep
        // Some calculations are only required once per timestep
        // For example the amount of charge added to the system
        // Perform these calculations here
        if self.flag_first_call {
            self.first_call();
        }

        // Distributing energy demand through all units
        let energy_demand = energy_output_required / self.n_units as f64;

        let q_out_ts = self.q_out_ts.unwrap();

        // Create power variables and assign the values just calculated at average of timestep
        let e_out = min_of_2(
            energy_demand,
            q_out_ts * (timestep - self.total_time_running_current_timestamp),
        );
        let time_running_current_service = if q_out_ts > 0. { e_out / q_out_ts } else { 0. };

        let q_loss_ts = self.q_loss_ts.unwrap();
        let q_in_ts = self.q_in_ts.unwrap();

        let e_loss = q_loss_ts * time_running_current_service;
        let mut e_in = q_in_ts * time_running_current_service;

        let mut delta_charge_level = (e_in - e_out - e_loss) / self.heat_storage_capacity;

        // Calculate new charge level after accounting for energy in and out and cap at target_charge
        charge_level += delta_charge_level;
        if charge_level > target_charge {
            e_in -= (charge_level - target_charge) * self.heat_storage_capacity;
            if e_in < 0.0 {
                e_in = 0.;
                charge_level -= delta_charge_level;
                delta_charge_level = -(e_out + e_loss) / self.heat_storage_capacity;
                charge_level += delta_charge_level;
            } else {
                charge_level = target_charge;
            }
        }

        let energy_output_provided = e_out;

        self.charge_level = charge_level;

        // TODO record energy supply demand

        self.total_time_running_current_timestamp += time_running_current_service;

        // Save results that are needed later (in the timestep_end function)
        self.service_results.push(HeatBatteryResult {
            service_name: service_name.to_string(),
            time_running: time_running_current_service,
            current_hb_power: q_out_ts,
        });

        energy_output_provided
    }

    /// Calculation of heat battery auxilary energy consumption
    fn calc_auxiliary_energy(&self, _timestep: f64, time_remaining_current_timestep: f64) {
        // Energy used by circulation pump
        let mut energy_aux = self.total_time_running_current_timestamp * self.power_circ_pump;

        // Energy used in standby mode
        energy_aux += self.power_standby * time_remaining_current_timestep;

        // TODO report energy supply demand
    }

    /// Calculations to be done at the end of each timestep
    pub fn timestep_end(&mut self) {
        let timestep = self.simulation_time.step_in_hours();
        let time_remaining_current_timestep = timestep - self.total_time_running_current_timestamp;

        if self.flag_first_call {
            self.first_call();
        }

        // Calculating auxiliary energy to provide services during timestep
        self.calc_auxiliary_energy(timestep, time_remaining_current_timestep);

        // Completing any charging left in the timestep and removing all losses from the charge level
        // Calculating heat battery losses in timestep to correct charge level
        // Currently assumed all losses are to the exterior independently of the
        // heat battery location
        let q_loss_ts = self.q_loss_ts.unwrap();
        let q_in_ts = self.q_in_ts.unwrap();

        let e_loss = q_loss_ts * time_remaining_current_timestep;
        let mut e_in = q_in_ts * time_remaining_current_timestep;

        let mut charge_level = self.charge_level;
        let target_charge = self.target_charge();
        let mut delta_charge_level = (e_in - e_loss) / self.heat_storage_capacity;

        // Calculate new charge level after accounting for energy in and out and cap at target_charge
        charge_level += delta_charge_level;
        if charge_level > target_charge {
            e_in -= (charge_level - target_charge) * self.heat_storage_capacity;
            if e_in < 0.0 {
                e_in = 0.;
                charge_level -= delta_charge_level;
                delta_charge_level = -e_loss / self.heat_storage_capacity;
                charge_level += delta_charge_level;
            } else {
                charge_level = target_charge;
            }
        }

        self.charge_level = charge_level;

        // TODO report energy supply demand

        let current_hour = self.simulation_time.current_hour();

        // Preparing Heat battery for next time step
        // Variables below need to be reset at the end of each timestep
        // Picking target charge level from control
        let time_range = (current_hour + 1) * HEAT_BATTERY_TIME_UNIT;

        let target_charge = self.target_charge();
        let mut charge_level_qin = self.charge_level;
        self.q_in_ts = Some(self.electric_charge(time_range as f64));
        let q_in_ts = self.q_in_ts.unwrap();
        // Calculate max charge level possible in next timestep
        if charge_level_qin < target_charge {
            delta_charge_level = q_in_ts * timestep / self.heat_storage_capacity;
            charge_level_qin += delta_charge_level;
            if charge_level_qin > target_charge {
                charge_level_qin = target_charge;
            }
        }

        // Estimating output rate at average of capacity in timestep
        let max_output = self.lab_test_rated_output(charge_level_qin);
        let delta_charge_level = max_output * timestep / self.heat_storage_capacity;
        self.q_out_ts =
            Some(self.lab_test_rated_output(charge_level_qin - delta_charge_level / 2.));
        self.q_loss_ts = Some(self.lab_test_losses(charge_level_qin - delta_charge_level / 2.));

        self.total_time_running_current_timestamp = Default::default();
        self.service_results = Default::default();
    }

    pub fn energy_output_max(&mut self, _temp_output: f64) -> f64 {
        let timestep = self.simulation_time.step_in_hours();
        let current_hour = self.simulation_time.current_hour();
        let time_range = current_hour * HEAT_BATTERY_TIME_UNIT;

        // Picking target charge level from control
        let target_charge = self.target_charge();
        let mut charge_level_qin = self.charge_level;
        self.q_in_ts = Some(self.electric_charge(time_range as f64));
        // Calculate max charge level possible in next timestep
        if charge_level_qin < target_charge {
            let delta_charge_level = self.q_in_ts.unwrap() * timestep / self.heat_storage_capacity;
            charge_level_qin += delta_charge_level;
            if charge_level_qin > target_charge {
                charge_level_qin = target_charge;
            }
        }

        // Estimating output rate at average of capacity in timestep
        let max_output = self.lab_test_rated_output(charge_level_qin);
        let delta_charge_level = max_output * timestep / self.heat_storage_capacity;

        self.lab_test_rated_output(charge_level_qin - delta_charge_level / 2.) * timestep
    }

    fn target_charge(&self) -> f64 {
        match self.charge_control.as_ref() {
            Control::ToUChargeControl(ctrl) => {
                ctrl.target_charge(&self.simulation_time.current_iteration())
            }
            _ => unreachable!(),
        }
    }
}
