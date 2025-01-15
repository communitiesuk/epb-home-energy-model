/// This module provides object(s) to model the behaviour of heat batteries.
use crate::compare_floats::min_of_2;
use crate::core::controls::time_control::{per_control, Control, ControlBehaviour};
use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::input::HeatSourceWetDetails;
use crate::simulation_time::{SimulationTimeIteration, SimulationTimeIterator};
use crate::statistics::np_interp;
use anyhow::bail;
use parking_lot::{Mutex, RwLock};
use std::collections::HashMap;
use std::sync::Arc;

pub enum ServiceType {
    WaterRegular,
    Space,
}

/// An object to represent a water heating service provided by a regular heat battery.
///
/// This object contains the parts of the heat battery calculation that are
/// specific to providing hot water.
#[derive(Clone, Debug)]
pub struct HeatBatteryServiceWaterRegular {
    heat_battery: Arc<Mutex<HeatBattery>>,
    service_name: String,
    control: Option<Arc<Control>>,
    control_min: Option<Arc<Control>>,
    control_max: Option<Arc<Control>>,
}

impl HeatBatteryServiceWaterRegular {
    /// Arguments:
    /// * `heat_battery` - reference to the Heat Battery object providing the service
    /// * `service_name` - name of the service demanding energy
    /// * `control_min` - optional reference to a control object which must select current the minimum timestep temperature
    /// * `control_max` - optional reference to a control object which must select current the maximum timestep temperature
    pub(crate) fn new(
        heat_battery: Arc<Mutex<HeatBattery>>,
        service_name: String,
        control_min: Option<Arc<Control>>,
        control_max: Option<Arc<Control>>,
    ) -> Self {
        let control = control_min.clone();

        Self {
            heat_battery,
            service_name,
            control,
            control_min,
            control_max,
        }
    }

    pub(crate) fn temp_setpnt(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> (Option<f64>, Option<f64>) {
        let control_min_setpnt = self
            .control_min
            .as_ref()
            .and_then(|control| control.setpnt(simulation_time_iteration));
        let control_max_setpnt = self
            .control_max
            .as_ref()
            .and_then(|control| control.setpnt(simulation_time_iteration));

        (control_min_setpnt, control_max_setpnt)
    }

    /// Demand energy (in kWh) from the heat_battery
    pub fn demand_energy(
        &self,
        energy_demand: f64,
        temp_return: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let service_on = self.is_on(simulation_time_iteration);
        let energy_demand = if !service_on { 0.0 } else { energy_demand };

        self.heat_battery.lock().demand_energy(
            &self.service_name,
            ServiceType::WaterRegular,
            energy_demand,
            temp_return,
            None,
            None,
            simulation_time_iteration.index,
        )
    }

    pub fn energy_output_max(
        &self,
        _temp_return: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let service_on = self.is_on(simulation_time_iteration);
        if !service_on {
            return Ok(0.0);
        }

        let control_max_setpnt = self
            .control_max
            .as_ref()
            .and_then(|control| control.setpnt(&simulation_time_iteration));

        self.heat_battery
            .lock()
            .energy_output_max(control_max_setpnt, None)
    }

    fn is_on(&self, simulation_time_iteration: SimulationTimeIteration) -> bool {
        if self.control.is_some() {
            per_control!(self.control.as_ref().unwrap().as_ref(), ctrl => { ctrl.is_on(&simulation_time_iteration) })
        } else {
            true
        }
    }
}

#[derive(Clone, Debug)]
pub struct HeatBatteryServiceSpace {
    heat_battery: Arc<Mutex<HeatBattery>>,
    service_name: String,
    control: Arc<Control>,
}

/// An object to represent a space heating service provided by a heat_battery to e.g. radiators.
///
/// This object contains the parts of the heat battery calculation that are
/// specific to providing space heating.
impl HeatBatteryServiceSpace {
    pub(crate) fn new(
        heat_battery: Arc<Mutex<HeatBattery>>,
        service_name: String,
        control: Arc<Control>,
    ) -> Self {
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
        &self,
        energy_demand: f64,
        _temp_flow: f64,
        temp_return: f64,
        time_start: Option<f64>,
        update_heat_source_state: Option<bool>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let _time_start = time_start.unwrap_or(0.);
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);

        if !self.is_on(simulation_time_iteration) {
            return Ok(0.0);
        }

        self.heat_battery.lock().demand_energy(
            &self.service_name,
            ServiceType::Space,
            energy_demand,
            temp_return,
            None,
            Some(update_heat_source_state),
            simulation_time_iteration.index,
        )
    }

    fn is_on(&self, simulation_time_iteration: SimulationTimeIteration) -> bool {
        per_control!(self.control.as_ref(), ctrl => { ctrl.is_on(&simulation_time_iteration) })
    }

    pub(crate) fn energy_output_max(
        &self,
        temp_output: f64,
        _temp_return_feed: f64,
        time_start: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let time_start = time_start.unwrap_or(0.);

        if !self.is_on(simtime) {
            return Ok(0.0);
        }

        self.heat_battery
            .lock()
            .energy_output_max(Some(temp_output), Some(time_start))
    }
}

const HEAT_BATTERY_TIME_UNIT: u32 = 3_600;

// nothing seems to read this - check upstream whether service_results field is necessary
#[derive(Clone, Debug)]
#[allow(dead_code)]
struct HeatBatteryResult {
    service_name: String,
    time_running: f64,
    current_hb_power: f64,
}

#[derive(Clone, Debug)]
pub struct HeatBattery {
    simulation_time: Arc<SimulationTimeIterator>,
    energy_supply: Arc<RwLock<EnergySupply>>,
    energy_supply_connection: EnergySupplyConnection,
    energy_supply_connections: HashMap<String, EnergySupplyConnection>,
    pwr_in: f64,
    heat_storage_capacity: f64,
    max_rated_heat_output: f64,
    max_rated_losses: f64,
    power_circ_pump: f64,
    power_standby: f64,
    n_units: usize,
    charge_control: Arc<Control>, // ChargeControl variant expected
    // nothing external seems to read this - check upstream whether service_results field is necessary
    service_results: Vec<HeatBatteryResult>,
    total_time_running_current_timestep: f64,
    flag_first_call: bool,
    charge_level: f64,
    q_in_ts: Option<f64>,
    q_out_ts: Option<f64>,
    q_loss_ts: Option<f64>,
    // TODO - check upstream whether it's an error that these numbers are not used
    labs_tests_rated_output: Vec<(f64, f64)>,
    labs_tests_rated_output_enhanced: Vec<(f64, f64)>,
    labs_tests_losses: Vec<(f64, f64)>,
}

impl HeatBattery {
    pub(crate) fn new(
        heat_battery_details: &HeatSourceWetDetails,
        charge_control: Arc<Control>,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_connection: EnergySupplyConnection,
        simulation_time: Arc<SimulationTimeIterator>,
    ) -> Self {
        let (
            pwr_in,
            heat_storage_capacity,
            max_rated_heat_output,
            max_rated_losses,
            power_circ_pump,
            power_standby,
            n_units,
            labs_tests_rated_output,
            labs_tests_rated_output_enhanced,
            labs_tests_losses,
            ..,
        ) = if let HeatSourceWetDetails::HeatBattery {
            rated_charge_power: pwr_in,
            heat_storage_capacity,
            max_rated_heat_output,
            max_rated_losses,
            electricity_circ_pump: power_circ_pump,
            electricity_standby: power_standby,
            number_of_units: n_units,
            labs_tests_rated_output,
            labs_tests_rated_output_enhanced,
            labs_tests_losses,
            ..
        } = heat_battery_details
        {
            (
                *pwr_in,
                *heat_storage_capacity,
                *max_rated_heat_output,
                *max_rated_losses,
                *power_circ_pump,
                *power_standby,
                *n_units,
                labs_tests_rated_output.clone(),
                labs_tests_rated_output_enhanced.clone(),
                labs_tests_losses.clone(),
            )
        } else {
            unreachable!()
        };

        Self {
            simulation_time,
            energy_supply,
            energy_supply_connection,
            energy_supply_connections: Default::default(),
            pwr_in,
            heat_storage_capacity,
            max_rated_heat_output,
            max_rated_losses,
            power_circ_pump,
            power_standby,
            n_units,
            charge_control,
            service_results: Default::default(),
            total_time_running_current_timestep: Default::default(),
            flag_first_call: true,
            charge_level: Default::default(),
            q_in_ts: Default::default(),
            q_out_ts: Default::default(),
            q_loss_ts: Default::default(),
            labs_tests_rated_output,
            labs_tests_rated_output_enhanced,
            labs_tests_losses,
        }
    }

    pub fn create_service_connection(
        heat_battery: Arc<Mutex<Self>>,
        service_name: &str,
    ) -> anyhow::Result<()> {
        if heat_battery
            .lock()
            .energy_supply_connections
            .contains_key(service_name)
        {
            bail!("Error: Service name already used: {service_name}");
        }
        let energy_supply = heat_battery.lock().energy_supply.clone();

        // Set up EnergySupplyConnection for this service
        heat_battery.lock().energy_supply_connections.insert(
            service_name.to_string(),
            EnergySupply::connection(energy_supply, service_name).unwrap(),
        );

        Ok(())
    }

    /// Return a HeatBatteryServiceWaterRegular object and create an EnergySupplyConnection for it
    ///
    /// Arguments:
    /// * `heat_battery` - reference to heat battery
    /// * `service_name` - name of the service demanding energy from the heat battery
    /// * `temp_hot_water` - temperature of the hot water to be provided, in deg C
    /// * `temp_limit_upper` - upper operating limit for temperature, in deg C
    /// * `cold_feed` - reference to ColdWaterSource object
    /// * `control` - reference to a control object
    pub(crate) fn create_service_hot_water_regular(
        heat_battery: Arc<Mutex<Self>>,
        service_name: &str,
        control_min: Option<Arc<Control>>,
        control_max: Option<Arc<Control>>,
    ) -> HeatBatteryServiceWaterRegular {
        Self::create_service_connection(heat_battery.clone(), service_name).unwrap();
        HeatBatteryServiceWaterRegular::new(
            heat_battery,
            service_name.to_string(),
            control_min,
            control_max,
        )
    }

    pub(crate) fn create_service_space_heating(
        heat_battery: Arc<Mutex<Self>>,
        service_name: &str,
        control: Arc<Control>,
    ) -> HeatBatteryServiceSpace {
        Self::create_service_connection(heat_battery.clone(), service_name).unwrap();
        HeatBatteryServiceSpace::new(heat_battery, service_name.to_string(), control)
    }

    /// Converts power value supplied to the correct units
    ///
    /// Arguments:
    /// * `power` - Power in watts
    /// * `timestep` - length of the timestep
    ///
    /// returns  -- Energy in kWH
    pub fn convert_to_energy(&self, power: f64, timestep: f64) -> f64 {
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
        let x = self
            .labs_tests_rated_output_enhanced
            .iter()
            .map(|row| row.0)
            .collect::<Vec<_>>();
        let y = self
            .labs_tests_rated_output_enhanced
            .iter()
            .map(|row| row.1)
            .collect::<Vec<_>>();

        np_interp(charge_level, &x, &y) * self.max_rated_heat_output
    }

    fn lab_test_losses(&self, charge_level: f64) -> f64 {
        let x = self
            .labs_tests_losses
            .iter()
            .map(|row| row.0)
            .collect::<Vec<_>>();
        let y = self
            .labs_tests_losses
            .iter()
            .map(|row| row.1)
            .collect::<Vec<_>>();
        np_interp(charge_level, &x, &y) * self.max_rated_losses
    }

    fn first_call(&mut self) -> anyhow::Result<()> {
        let timestep = self.simulation_time.step_in_hours();
        let current_hour = self.simulation_time.current_hour();
        let time_range = current_hour * HEAT_BATTERY_TIME_UNIT;
        let charge_level = self.charge_level;
        let mut charge_level_qin = charge_level;
        let target_charge = self.target_charge()?;

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

        Ok(())
    }

    /// Calculate time available for the current service
    fn time_available(&self, time_start: f64, timestep: f64) -> f64 {
        // Assumes that time spent on other services is evenly spread throughout
        // the timestep so the adjustment for start time below is a proportional
        // reduction of the overall time available, not simply a subtraction
        (timestep - self.total_time_running_current_timestep) * (1. - time_start / timestep)
    }

    pub fn demand_energy(
        &mut self,
        service_name: &str,
        _service_type: ServiceType,
        energy_output_required: f64,
        _temp_return_feed: f64,
        time_start: Option<f64>,
        update_heat_source_state: Option<bool>,
        timestep_idx: usize,
    ) -> anyhow::Result<f64> {
        let time_start = time_start.unwrap_or(0.);
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);
        let timestep = self.simulation_time.step_in_hours();
        let time_available = self.time_available(time_start, timestep);
        let mut charge_level = self.charge_level;

        // Picking target charge level from control
        let target_charge = self.target_charge()?;
        // __demand_energy is called for each service in each timestep
        // Some calculations are only required once per timestep
        // For example the amount of charge added to the system
        // Perform these calculations here
        if self.flag_first_call {
            self.first_call()?;
        }

        // Distributing energy demand through all units
        let energy_demand = energy_output_required / self.n_units as f64;

        let q_out_ts = self.q_out_ts.unwrap();

        // Create power variables and assign the values just calculated at average of timestep
        let e_out = min_of_2(energy_demand, q_out_ts * time_available);
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

        if update_heat_source_state {
            self.charge_level = charge_level;

            self.energy_supply_connection
                .demand_energy(e_in * self.n_units as f64, timestep_idx)
                .unwrap();
            self.energy_supply_connections[service_name]
                .energy_out(e_out * self.n_units as f64, timestep_idx)
                .unwrap();

            self.total_time_running_current_timestep += time_running_current_service;

            // Save results that are needed later (in the timestep_end function)
            self.service_results.push(HeatBatteryResult {
                service_name: service_name.to_string(),
                time_running: time_running_current_service,
                current_hb_power: q_out_ts,
            });
        }

        Ok(energy_output_provided)
    }

    /// Calculation of heat battery auxilary energy consumption
    fn calc_auxiliary_energy(
        &self,
        _timestep: f64,
        time_remaining_current_timestep: f64,
        timestep_idx: usize,
    ) {
        // Energy used by circulation pump
        let mut energy_aux = self.total_time_running_current_timestep * self.power_circ_pump;

        // Energy used in standby mode
        energy_aux += self.power_standby * time_remaining_current_timestep;

        self.energy_supply_connection
            .demand_energy(energy_aux, timestep_idx)
            .unwrap();
    }

    /// Calculations to be done at the end of each timestep
    pub fn timestep_end(&mut self, timestep_idx: usize) -> anyhow::Result<()> {
        let timestep = self.simulation_time.step_in_hours();
        let time_remaining_current_timestep = timestep - self.total_time_running_current_timestep;

        if self.flag_first_call {
            self.first_call()?;
        }

        // Calculating auxiliary energy to provide services during timestep
        self.calc_auxiliary_energy(timestep, time_remaining_current_timestep, timestep_idx);

        // Completing any charging left in the timestep and removing all losses from the charge level
        // Calculating heat battery losses in timestep to correct charge level
        // Currently assumed all losses are to the exterior independently of the
        // heat battery location
        let q_loss_ts = self.q_loss_ts.unwrap();
        let q_in_ts = self.q_in_ts.unwrap();

        let e_loss = q_loss_ts * time_remaining_current_timestep;
        let mut e_in = q_in_ts * time_remaining_current_timestep;

        let mut charge_level = self.charge_level;
        let target_charge = self.target_charge()?;
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

        self.energy_supply_connection
            .demand_energy(e_in * self.n_units as f64, timestep_idx)
            .unwrap();

        let current_hour = self.simulation_time.current_hour();

        // Preparing Heat battery for next time step
        // Variables below need to be reset at the end of each timestep
        // Picking target charge level from control
        let time_range = (current_hour + 1) * HEAT_BATTERY_TIME_UNIT;

        let target_charge = self.target_charge()?;
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

        self.total_time_running_current_timestep = Default::default();
        self.service_results = Default::default();

        Ok(())
    }

    /// Calculate the maximum energy output of the heat battery, accounting
    /// for time spent on higher-priority services.
    pub fn energy_output_max(
        &mut self,
        _temp_output: Option<f64>,
        time_start: Option<f64>,
    ) -> anyhow::Result<f64> {
        let time_start = time_start.unwrap_or(0.);
        let timestep = self.simulation_time.step_in_hours();
        let time_available = self.time_available(time_start, timestep);
        let current_hour = self.simulation_time.current_hour();
        let time_range = current_hour * HEAT_BATTERY_TIME_UNIT;

        // Picking target charge level from control
        let target_charge = self.target_charge()?;
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

        // Estimating output rate at average of capacity in time_available
        let max_output = self.lab_test_rated_output(charge_level_qin);
        let delta_charge_level = max_output * time_available / self.heat_storage_capacity;

        Ok(self.lab_test_rated_output(charge_level_qin - delta_charge_level / 2.) * time_available)
    }

    fn target_charge(&self) -> anyhow::Result<f64> {
        match self.charge_control.as_ref() {
            Control::Charge(ctrl) => {
                ctrl.target_charge(self.simulation_time.current_iteration(), None)
            }
            _ => unreachable!(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::core::controls::time_control::SetpointTimeControl;
    use crate::core::controls::time_control::{ChargeControl, Control};
    use crate::core::energy_supply::energy_supply::{
        EnergySupply, EnergySupplyBuilder, EnergySupplyConnection,
    };
    use crate::core::heating_systems::heat_battery::HeatBattery;
    use crate::core::heating_systems::heat_battery::HeatBatteryServiceSpace;
    use crate::core::heating_systems::heat_battery::HeatBatteryServiceWaterRegular;
    use crate::core::heating_systems::heat_battery::ServiceType;
    use crate::external_conditions::{DaylightSavingsConfig, ExternalConditions};
    use crate::input::HeatSourceWetDetails;
    use crate::input::{ControlLogicType, HeatSourceLocation};
    use crate::input::{ExternalSensor, FuelType};
    use crate::simulation_time::SimulationTimeIteration;
    use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
    use approx::assert_relative_eq;
    use itertools::Itertools;
    use parking_lot::{Mutex, RwLock};
    use rstest::fixture;
    use rstest::rstest;
    use serde_json::json;
    use std::sync::Arc;

    const SERVICE_NAME: &str = "TestService";
    const TEMP_HOT_WATER: f64 = 55.;

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 2., 1.)
    }

    #[fixture]
    pub fn simulation_time_iterator(
        simulation_time: SimulationTime,
    ) -> Arc<SimulationTimeIterator> {
        simulation_time.iter().into()
    }

    #[fixture]
    pub fn simulation_time_iteration(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) -> SimulationTimeIteration {
        simulation_time_iterator.current_iteration()
    }

    #[fixture]
    fn external_sensor() -> ExternalSensor {
        serde_json::from_value(json!({
            "correlation": [
                {"temperature": 0.0, "max_charge": 1.0},
                {"temperature": 10.0, "max_charge": 0.9},
                {"temperature": 18.0, "max_charge": 0.0}
            ]
        }))
        .unwrap()
    }

    #[fixture]
    fn external_conditions(simulation_time: SimulationTime) -> ExternalConditions {
        ExternalConditions::new(
            &simulation_time.iter(),
            vec![0.0, 2.5],
            vec![3.7, 3.8],
            vec![200., 220.],
            vec![333., 610.],
            vec![420., 750.],
            vec![0.2; 8760],
            51.42,
            -0.75,
            0,
            0,
            Some(0),
            1.,
            Some(1),
            Some(DaylightSavingsConfig::NotApplicable),
            false,
            false,
            // following shading segments are corrected from upstream Python, which uses angles measured from wrong origin
            serde_json::from_value(json!(
                [
                    {"number": 1, "start360": 0, "end360": 45},
                    {"number": 2, "start360": 45, "end360": 90},
                ]
            ))
            .unwrap(),
        )
    }

    #[fixture]
    fn battery_control_off(
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
    ) -> Control {
        create_control_with_value(false, external_conditions, external_sensor)
    }

    #[fixture]
    fn battery_control_on(
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
    ) -> Control {
        create_control_with_value(true, external_conditions, external_sensor)
    }

    fn create_control_with_value(
        boolean: bool,
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
    ) -> Control {
        Control::Charge(
            ChargeControl::new(
                ControlLogicType::Manual,
                vec![boolean],
                1.,
                0,
                1.,
                vec![0.2],
                None,
                None,
                None,
                None,
                external_conditions.into(),
                Some(external_sensor),
            )
            .unwrap(),
        )
    }

    fn create_heat_battery(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        control: Control,
    ) -> Arc<Mutex<HeatBattery>> {
        let labs_tests_rated_output = vec![
            (0.0, 0.0),
            (0.08, 0.00),
            (0.17, 0.05),
            (0.19, 0.10),
            (0.21, 0.15),
            (0.23, 0.21),
            (0.25, 0.23),
            (0.28, 0.26),
            (0.31, 0.29),
            (0.34, 0.32),
            (0.38, 0.36),
            (0.42, 0.41),
            (0.47, 0.45),
            (0.52, 0.51),
            (0.58, 0.57),
            (0.64, 0.64),
            (0.72, 0.71),
            (0.8, 0.8),
            (0.89, 0.89),
            (1.0, 1.0),
        ];
        let labs_tests_rated_output_enhanced = vec![
            (0.0, 0.0),
            (0.101, 0.0),
            (0.12, 0.18),
            (0.144, 0.235),
            (0.175, 0.313),
            (0.215, 0.391),
            (0.266, 0.486),
            (0.328, 0.607),
            (0.406, 0.728),
            (0.494, 0.795),
            (0.587, 0.825),
            (0.683, 0.875),
            (0.781, 0.906),
            (0.891, 0.953),
            (0.981, 0.992),
            (1.0, 1.0),
        ];
        let labs_tests_losses = vec![
            (0.0, 0.),
            (0.16, 0.13),
            (0.17, 0.15),
            (0.19, 0.17),
            (0.21, 0.18),
            (0.23, 0.21),
            (0.25, 0.23),
            (0.28, 0.26),
            (0.31, 0.29),
            (0.34, 0.32),
            (0.38, 0.36),
            (0.42, 0.41),
            (0.47, 0.45),
            (0.52, 0.51),
            (0.58, 0.57),
            (0.64, 0.64),
            (0.72, 0.71),
            (0.8, 0.8),
            (0.89, 0.89),
            (1.0, 1.0),
        ];
        let heat_battery_details: &HeatSourceWetDetails = &HeatSourceWetDetails::HeatBattery {
            energy_supply: "mains elec".to_string(),
            heat_battery_location: HeatSourceLocation::Internal,
            electricity_circ_pump: 0.06,
            electricity_standby: 0.0244,
            rated_charge_power: 20.0,
            heat_storage_capacity: 80.0,
            max_rated_heat_output: 15.0,
            max_rated_losses: 0.22,
            number_of_units: 1,
            control_charge: "hb_charge_control".into(),
            labs_tests_rated_output,
            labs_tests_rated_output_enhanced,
            labs_tests_losses,
        };

        let energy_supply: Arc<RwLock<EnergySupply>> = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time_iterator.total_steps())
                .build(),
        ));

        let energy_supply_connection: EnergySupplyConnection =
            EnergySupply::connection(energy_supply.clone(), "WaterHeating").unwrap();

        let heat_battery = Arc::new(Mutex::new(HeatBattery::new(
            heat_battery_details,
            control.into(),
            energy_supply,
            energy_supply_connection,
            simulation_time_iterator,
        )));

        HeatBattery::create_service_connection(heat_battery.clone(), SERVICE_NAME).unwrap();

        heat_battery
    }

    fn create_setpoint_time_control(schedule: Vec<Option<f64>>) -> Control {
        Control::SetpointTime(
            SetpointTimeControl::new(schedule, 0, 1., None, None, None, None, 1.).unwrap(),
        )
    }

    fn get_service_names_from_results(heat_battery: Arc<Mutex<HeatBattery>>) -> Vec<String> {
        heat_battery
            .lock()
            .service_results
            .iter()
            .map(|result| result.service_name.clone())
            .collect_vec()
    }

    // in Python this test is called test_service_is_on_with_control
    #[rstest]
    fn test_service_is_on_when_service_control_is_on(
        simulation_time_iteration: SimulationTimeIteration,
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        // Test when controlvent is provided and returns True
        let service_control_on: Control = create_setpoint_time_control(vec![Some(21.0)]);

        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off);

        let heat_battery_service = HeatBatteryServiceSpace::new(
            heat_battery.clone(),
            SERVICE_NAME.into(),
            service_control_on.into(),
        );

        assert!(heat_battery_service.is_on(simulation_time_iteration));

        let service_control_off: Control = create_setpoint_time_control(vec![None]);

        let heat_battery_service: HeatBatteryServiceSpace = HeatBatteryServiceSpace::new(
            heat_battery,
            SERVICE_NAME.into(),
            service_control_off.into(),
        );

        assert!(!heat_battery_service.is_on(simulation_time_iteration));
    }

    // test_service_is_on_without_control
    #[rstest]
    fn test_service_with_no_service_control_is_always_on_for_water_regular(
        simulation_time_iteration: SimulationTimeIteration,
        battery_control_off: Control,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off);
        let control_min = create_setpoint_time_control(vec![
            Some(52.),
            None,
            None,
            None,
            Some(52.),
            Some(52.),
            Some(52.),
            Some(52.),
        ]);
        let control_max = create_setpoint_time_control(vec![
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
        ]);

        let heat_battery_service: HeatBatteryServiceWaterRegular =
            HeatBatteryServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                Some(Arc::new(control_min)),
                Some(Arc::new(control_max)),
            );

        assert!(heat_battery_service.is_on(simulation_time_iteration));
    }

    #[rstest]
    fn test_demand_energy_when_service_control_on_for_water_regular(
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
        let energy_demand = 10.;
        let temp_return = 40.;

        let control_min = create_setpoint_time_control(vec![
            Some(52.),
            None,
            None,
            None,
            Some(52.),
            Some(52.),
            Some(52.),
            Some(52.),
        ]);
        let control_max = create_setpoint_time_control(vec![
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
        ]);

        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        let heat_battery_service: HeatBatteryServiceWaterRegular =
            HeatBatteryServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                Some(Arc::new(control_min)),
                Some(Arc::new(control_max)),
            );

        let result = heat_battery_service
            .demand_energy(energy_demand, temp_return, simulation_time_iteration)
            .unwrap();

        assert_relative_eq!(result, 4.358566028225806);
    }

    // In Python this is test_demand_energy_service_off
    #[rstest]
    fn test_demand_energy_returns_zero_when_service_control_is_off_for_water_regular(
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
        let energy_demand = 10.;
        let temp_return = 40.;

        let service_control_off = create_setpoint_time_control(vec![None]);

        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        let heat_battery_service: HeatBatteryServiceWaterRegular =
            HeatBatteryServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                Some(service_control_off.into()),
                None,
            );

        let result = heat_battery_service
            .demand_energy(energy_demand, temp_return, simulation_time_iteration)
            .unwrap();

        assert_eq!(result, 0.);
    }

    // In Python this is test_energy_output_max_service_on
    #[rstest]
    fn test_energy_output_max_when_service_control_on_for_water_regular(
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let temp_return = 40.;
        let battery_control_on: Control = Control::Charge(
            ChargeControl::new(
                ControlLogicType::Manual,
                vec![true, true, true],
                1.,
                0,
                1.,
                vec![1.5, 1.6], // these values change the result
                None,
                None,
                None,
                None,
                external_conditions.into(),
                Some(external_sensor),
            )
            .unwrap(),
        );
        let control_min = create_setpoint_time_control(vec![
            Some(52.),
            None,
            None,
            None,
            Some(52.),
            Some(52.),
            Some(52.),
            Some(52.),
        ]);
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        let heat_battery_service: HeatBatteryServiceWaterRegular =
            HeatBatteryServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                None,
                Some(Arc::new(control_min)),
            );

        let result = heat_battery_service
            .energy_output_max(temp_return, simulation_time_iteration)
            .unwrap();

        assert_relative_eq!(result, 5.637774816176471);
    }

    #[rstest]
    fn test_energy_output_max_service_off_for_water_regular(
        // In Python this is test_energy_output_max_service_off
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_off: Control,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off);
        let control_min = create_setpoint_time_control(vec![
            Some(52.),
            None,
            None,
            None,
            Some(52.),
            Some(52.),
            Some(52.),
            Some(52.),
        ]);
        let control_max = create_setpoint_time_control(vec![
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
        ]);

        let heat_battery_service: HeatBatteryServiceWaterRegular =
            HeatBatteryServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                Some(Arc::new(control_min)),
                Some(Arc::new(control_max)),
            );

        let temp_return = 40.;
        let result = heat_battery_service
            .energy_output_max(temp_return, simulation_time_iteration)
            .unwrap();

        assert_eq!(result, 0.);
    }

    #[rstest]
    fn test_temp_setpnt_for_space(
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_off: Control,
    ) {
        let first_scheduled_temp = Some(21.);
        let ctrl: Control = create_setpoint_time_control(vec![first_scheduled_temp]);
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off);
        let heat_battery_space =
            HeatBatteryServiceSpace::new(heat_battery, SERVICE_NAME.into(), ctrl.into());

        assert_eq!(
            heat_battery_space.temp_setpnt(simulation_time_iteration),
            first_scheduled_temp
        );
    }

    #[rstest]
    fn test_in_required_period_for_space(
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_off: Control,
    ) {
        let ctrl: Control = create_setpoint_time_control(vec![Some(21.)]);
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off);
        let heat_battery_space =
            HeatBatteryServiceSpace::new(heat_battery, SERVICE_NAME.into(), ctrl.into());

        assert_eq!(
            heat_battery_space.in_required_period(simulation_time_iteration),
            Some(true)
        );
    }

    // TODO test_demand_energy_for_space

    #[rstest]
    fn test_demand_energy_service_off_for_space(
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_off: Control,
    ) {
        let energy_demand = 10.;
        let temp_return = 40.;
        let temp_flow = 1.;
        let ctrl: Control = create_setpoint_time_control(vec![None]);
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off);
        let heat_battery_space =
            HeatBatteryServiceSpace::new(heat_battery, SERVICE_NAME.into(), ctrl.into());
        let result = heat_battery_space
            .demand_energy(
                energy_demand,
                temp_return,
                temp_flow,
                None,
                None,
                simulation_time_iteration,
            )
            .unwrap();
        assert_eq!(result, 0.);
    }

    // in Python this test is called test_energy_output_max_service_on
    #[rstest]
    fn test_energy_output_max_service_on_for_space(
        battery_control_on: Control,
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let temp_output = 70.;
        let temp_return = 40.;
        let time_start = 0.1;
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        let service_control_on: Control =
            create_setpoint_time_control(vec![Some(21.0), Some(21.0)]);

        let heat_battery_service: HeatBatteryServiceSpace = HeatBatteryServiceSpace::new(
            heat_battery,
            SERVICE_NAME.into(),
            service_control_on.into(),
        );

        let result = heat_battery_service
            .energy_output_max(
                temp_output,
                temp_return,
                Some(time_start),
                simulation_time_iteration,
            )
            .unwrap();

        assert_relative_eq!(result, 4.037907837701614);
    }

    // in Python this test is called test_energy_output_max_service_off
    #[rstest]
    fn test_energy_output_max_service_off_for_space(
        battery_control_on: Control,
        simulation_time_iteration: SimulationTimeIteration,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        let temp_output = 70.;
        let temp_return = 40.;
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        let service_control_off: Control = create_setpoint_time_control(vec![None]);

        let heat_battery_service: HeatBatteryServiceSpace = HeatBatteryServiceSpace::new(
            heat_battery,
            SERVICE_NAME.into(),
            service_control_off.into(),
        );

        let result = heat_battery_service
            .energy_output_max(temp_output, temp_return, None, simulation_time_iteration)
            .unwrap();

        assert_relative_eq!(result, 0.);
    }

    #[rstest]
    fn test_create_service_connection(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        let create_connection_result =
            HeatBattery::create_service_connection(heat_battery.clone(), "new service");
        assert!(create_connection_result.is_ok());
        let create_connection_result =
            HeatBattery::create_service_connection(heat_battery, "new service");
        assert!(create_connection_result.is_err()) // second attempt to create a service connection with same name should error
    }

    #[rstest]
    fn test_convert_to_energy(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
        let power = 10.;
        let timestep = 0.25;

        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        let result = heat_battery.lock().convert_to_energy(power, timestep);
        assert_relative_eq!(result, 2.5);
    }

    #[rstest]
    fn test_electric_charge(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_off: Control,
        battery_control_on: Control,
    ) {
        // electric charge should be 0 when battery control is off
        let heat_battery =
            create_heat_battery(simulation_time_iterator.clone(), battery_control_off);
        assert_relative_eq!(heat_battery.lock().electric_charge(3600.), 0.0);

        // electric charge should be calculated when battery control is on
        let heat_battery =
            create_heat_battery(simulation_time_iterator.clone(), battery_control_on);
        assert_relative_eq!(heat_battery.lock().electric_charge(10.), 20.0);
    }

    #[rstest]
    fn test_lab_test_rated_output(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_off: Control,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_off);

        assert_relative_eq!(heat_battery.lock().lab_test_rated_output(5.), 15.);
    }

    #[rstest]
    fn test_first_call(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
        simulation_time: SimulationTime,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        for (t_idx, _) in simulation_time.iter().enumerate() {
            heat_battery.lock().first_call().unwrap();

            assert!(!heat_battery.lock().flag_first_call);
            assert_relative_eq!(heat_battery.lock().q_in_ts.unwrap(), 20.);
            assert_relative_eq!(heat_battery.lock().q_out_ts.unwrap(), 4.358566028225806);
            assert_relative_eq!(heat_battery.lock().q_loss_ts.unwrap(), 0.031277812499999995);

            heat_battery.lock().timestep_end(t_idx).unwrap();
        }
    }

    #[rstest]
    fn test_demand_energy(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        simulation_time: SimulationTime,
        battery_control_on: Control,
    ) {
        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);

        let service_name = "new_service";
        HeatBattery::create_service_connection(heat_battery.clone(), service_name).unwrap();

        for (t_idx, _) in simulation_time.iter().enumerate() {
            let demand_energy_actual = heat_battery
                .clone()
                .lock()
                .demand_energy(
                    service_name,
                    ServiceType::WaterRegular,
                    5.,
                    40.,
                    None,
                    None,
                    t_idx,
                )
                .unwrap();

            assert_relative_eq!(demand_energy_actual, 4.358566028225806);

            let service_names_in_results = get_service_names_from_results(heat_battery.clone());

            assert!(service_names_in_results.contains(&service_name.into()));

            assert_relative_eq!(
                heat_battery.lock().charge_level,
                [0.19512695199092742, 0.2][t_idx]
            );

            assert_relative_eq!(heat_battery.lock().total_time_running_current_timestep, 1.);

            heat_battery.lock().timestep_end(t_idx).unwrap();
        }
    }

    #[rstest]
    /// Check heat battery auxilary energy consumption
    fn test_calc_auxiliary_energy(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        battery_control_on: Control,
    ) {
        let heat_battery =
            create_heat_battery(simulation_time_iterator.clone(), battery_control_on);

        heat_battery.lock().calc_auxiliary_energy(
            1.0,
            0.5,
            simulation_time_iterator.current_index(),
        );

        let results_by_end_user = heat_battery
            .lock()
            .energy_supply
            .read()
            .results_by_end_user();

        let end_user_name = heat_battery
            .lock()
            .energy_supply_connection
            .end_user_name
            .clone();

        let results_by_end_user = results_by_end_user.get(&end_user_name).unwrap();

        assert_eq!(*results_by_end_user, vec![0.0122, 0.]);
    }

    #[rstest]
    fn test_timestep_end(
        external_sensor: ExternalSensor,
        external_conditions: ExternalConditions,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) {
        // not using the fixture here
        // because we need to set different charge_levels
        let battery_control_on: Control = Control::Charge(
            ChargeControl::new(
                ControlLogicType::Manual,
                vec![true, true, true],
                1.,
                0,
                1.,
                vec![1.0, 1.5],
                None,
                None,
                None,
                None,
                external_conditions.into(),
                Some(external_sensor),
            )
            .unwrap(),
        );

        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);
        let service_name = "new_timestep_end_service";
        HeatBattery::create_service_connection(heat_battery.clone(), service_name).unwrap();

        let t_idx = 0;
        heat_battery
            .lock()
            .demand_energy(
                service_name,
                ServiceType::WaterRegular,
                5.0,
                40.,
                None,
                None,
                t_idx,
            )
            .unwrap();

        assert_relative_eq!(heat_battery.lock().q_in_ts.unwrap(), 20.);
        assert_relative_eq!(heat_battery.lock().q_out_ts.unwrap(), 5.637774816176471);
        assert_relative_eq!(heat_battery.lock().q_loss_ts.unwrap(), 0.03929547794117647);
        assert_relative_eq!(
            heat_battery.lock().total_time_running_current_timestep,
            0.8868747268254661
        );

        let service_names_in_results = get_service_names_from_results(heat_battery.clone());

        assert!(service_names_in_results.contains(&service_name.into()));

        heat_battery.lock().timestep_end(t_idx).unwrap();

        assert!(!heat_battery.lock().flag_first_call);
        assert_relative_eq!(heat_battery.lock().q_in_ts.unwrap(), 20.);
        assert_relative_eq!(heat_battery.lock().q_out_ts.unwrap(), 10.001923317091928);
        assert_relative_eq!(heat_battery.lock().q_loss_ts.unwrap(), 0.07624000227068732);
        assert_relative_eq!(heat_battery.lock().total_time_running_current_timestep, 0.0);
        assert_eq!(heat_battery.lock().service_results.len(), 0);
    }

    #[rstest]
    fn test_energy_output_max(
        external_conditions: ExternalConditions,
        external_sensor: ExternalSensor,
        simulation_time_iterator: Arc<SimulationTimeIterator>,
        simulation_time: SimulationTime,
    ) {
        // not using the fixture here
        // because we need to set different charge_levels
        let battery_control_on: Control = Control::Charge(
            ChargeControl::new(
                ControlLogicType::Manual,
                vec![true, true, true],
                1.,
                0,
                1.,
                vec![1.5, 1.6], // these values change the result
                None,
                None,
                None,
                None,
                external_conditions.into(),
                Some(external_sensor),
            )
            .unwrap(),
        );

        let heat_battery = create_heat_battery(simulation_time_iterator, battery_control_on);

        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                heat_battery
                    .lock()
                    .energy_output_max(Some(0.), None)
                    .unwrap(),
                [5.637774816176471, 11.13482970854502][t_idx]
            );

            heat_battery.lock().timestep_end(t_idx).unwrap();
        }
    }
}
