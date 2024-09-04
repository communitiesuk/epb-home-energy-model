use crate::compare_floats::min_of_2;
use crate::core::controls::time_control::{per_control, Control, ControlBehaviour};
use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::input::HeatSourceWetDetails;
use crate::simulation_time::{SimulationTimeIteration, SimulationTimeIterator};
use anyhow::bail;
use interp::interp;
use parking_lot::{Mutex, RwLock};
use std::collections::HashMap;
use std::sync::Arc;

/// This module provides object(s) to model the behaviour of heat batteries.

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
    temp_hot_water: f64,
}

impl HeatBatteryServiceWaterRegular {
    pub fn new(
        heat_battery: Arc<Mutex<HeatBattery>>,
        service_name: String,
        temp_hot_water: f64,
        control: Option<Arc<Control>>,
    ) -> Self {
        Self {
            heat_battery,
            service_name,
            control,
            temp_hot_water,
        }
    }

    /// Demand energy (in kWh) from the heat_battery
    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        temp_return: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        let service_on = self.is_on(simulation_time_iteration);
        let energy_demand = if !service_on { 0.0 } else { energy_demand };

        self.heat_battery.lock().demand_energy(
            &self.service_name,
            ServiceType::WaterRegular,
            energy_demand,
            temp_return,
            simulation_time_iteration.index,
        )
    }

    pub fn energy_output_max(
        &self,
        _temp_return: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        let service_on = self.is_on(simulation_time_iteration);
        if !service_on {
            return 0.0;
        }

        self.heat_battery
            .lock()
            .energy_output_max(self.temp_hot_water)
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
    pub fn new(
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
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        if !self.is_on(simulation_time_iteration) {
            return 0.0;
        }

        self.heat_battery.lock().demand_energy(
            &self.service_name,
            ServiceType::Space,
            energy_demand,
            temp_return,
            simulation_time_iteration.index,
        )
    }

    fn is_on(&self, simulation_time_iteration: SimulationTimeIteration) -> bool {
        per_control!(self.control.as_ref(), ctrl => { ctrl.is_on(&simulation_time_iteration) })
    }
}

// TODO - check upstream whether it's an error that these numbers are not used
const _LABS_TESTS_RATED_OUTPUT: [[f64; 2]; 21] = [
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
    charge_control: Arc<Control>, // ToUChargeControl variant expected
    // nothing external seems to read this - check upstream whether service_results field is necessary
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
            ..,
        ) = if let HeatSourceWetDetails::HeatBattery {
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
            total_time_running_current_timestamp: Default::default(),
            flag_first_call: true,
            charge_level: Default::default(),
            q_in_ts: Default::default(),
            q_out_ts: Default::default(),
            q_loss_ts: Default::default(),
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
    pub fn create_service_hot_water_regular(
        heat_battery: Arc<Mutex<Self>>,
        service_name: &str,
        temp_hot_water: f64,
        control: Option<Arc<Control>>,
    ) -> HeatBatteryServiceWaterRegular {
        Self::create_service_connection(heat_battery.clone(), service_name).unwrap();
        HeatBatteryServiceWaterRegular::new(
            heat_battery,
            service_name.to_string(),
            temp_hot_water,
            control,
        )
    }

    pub fn create_service_space_heating(
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
    fn _convert_to_energy(&self, power: f64, timestep: f64) -> f64 {
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
        timestep_idx: usize,
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

        self.energy_supply_connection
            .demand_energy(e_in * self.n_units as f64, timestep_idx)
            .unwrap();
        self.energy_supply_connections[service_name]
            .energy_out(e_out * self.n_units as f64, timestep_idx)
            .unwrap();

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
    fn calc_auxiliary_energy(
        &self,
        _timestep: f64,
        time_remaining_current_timestep: f64,
        timestep_idx: usize,
    ) {
        // Energy used by circulation pump
        let mut energy_aux = self.total_time_running_current_timestamp * self.power_circ_pump;

        // Energy used in standby mode
        energy_aux += self.power_standby * time_remaining_current_timestep;

        self.energy_supply_connection
            .demand_energy(energy_aux, timestep_idx)
            .unwrap();
    }

    /// Calculations to be done at the end of each timestep
    pub fn timestep_end(&mut self, timestep_idx: usize) {
        let timestep = self.simulation_time.step_in_hours();
        let time_remaining_current_timestep = timestep - self.total_time_running_current_timestamp;

        if self.flag_first_call {
            self.first_call();
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

        self.energy_supply_connection
            .demand_energy(e_in * self.n_units as f64, timestep_idx)
            .unwrap();

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

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use crate::core::controls::time_control::Control;
    use crate::core::controls::time_control::SetpointTimeControl;
    use crate::core::controls::time_control::ToUChargeControl;
    use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
    use crate::core::heating_systems::heat_battery::HeatBattery;
    use crate::core::heating_systems::heat_battery::HeatBatteryServiceSpace;
    use crate::core::heating_systems::heat_battery::HeatBatteryServiceWaterRegular;
    use crate::input::EnergySupplyType;
    use crate::input::FuelType;
    use crate::input::HeatSourceLocation;
    use crate::input::HeatSourceWetDetails;
    use crate::simulation_time::SimulationTimeIteration;
    use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
    use parking_lot::{Mutex, RwLock};
    use rstest::fixture;
    use rstest::rstest;

    const SERVICE_NAME: &str = "TestService";
    const TEMP_HOT_WATER: f64 = 55.;

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 8., 1.)
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
    pub fn heat_battery(
        simulation_time_iterator: Arc<SimulationTimeIterator>,
    ) -> Arc<Mutex<HeatBattery>> {
        let heat_battery_details: &HeatSourceWetDetails = &HeatSourceWetDetails::HeatBattery {
            energy_supply: EnergySupplyType::Electricity,
            heat_battery_location: HeatSourceLocation::Internal,
            electricity_circ_pump: 0.06,
            electricity_standby: 0.0244,
            rated_charge_power: 20.0,
            heat_storage_capacity: 80.0,
            max_rated_heat_output: 15.0,
            max_rated_losses: 0.22,
            number_of_units: 1,
            control_charge: "hb_charge_control".into(),
        };

        let charge_control: Control = Control::ToUChargeControl(ToUChargeControl {
            schedule: vec![false],
            start_day: 0,
            time_series_step: 1.,
            charge_level: vec![0.2],
        });

        let energy_supply: Arc<RwLock<EnergySupply>> = Arc::new(RwLock::new(EnergySupply::new(
            FuelType::MainsGas,
            1,
            None,
            None,
            None,
        )));
        let energy_supply_connection: EnergySupplyConnection =
            EnergySupplyConnection::new(energy_supply.clone(), "WaterHeating".into());

        let heat_battery = HeatBattery::new(
            heat_battery_details,
            charge_control.into(),
            energy_supply,
            energy_supply_connection,
            simulation_time_iterator,
        );
        Mutex::new(heat_battery.into()).into()
    }

    fn create_setpoint_time_control(schedule: Vec<Option<f64>>) -> Control {
        Control::SetpointTimeControl(
            SetpointTimeControl::new(schedule, 0, 1., None, None, None, None, 1.)
                .unwrap(),
        )
    }

    #[rstest]
    fn test_service_is_on_with_control(
        simulation_time_iteration: SimulationTimeIteration,
        heat_battery: Arc<Mutex<HeatBattery>>,
    ) {
        // Test when controlvent is provided and returns True
        let ctrl_with_scheduled_temps: Control = create_setpoint_time_control(vec![Some(21.0)]);

        let heat_battery_service = HeatBatteryServiceSpace::new(
            heat_battery.clone(),
            SERVICE_NAME.into(),
            ctrl_with_scheduled_temps.into(),
        );

        assert_eq!(
            heat_battery_service.is_on(simulation_time_iteration),
            true
        );

        let ctrl_with_no_scheduled_temps: Control = create_setpoint_time_control(vec![None]);

        let heat_battery_service: HeatBatteryServiceSpace = HeatBatteryServiceSpace::new(
            heat_battery,
            SERVICE_NAME.into(),
            ctrl_with_no_scheduled_temps.into(),
        );

        assert_eq!(
            heat_battery_service.is_on(simulation_time_iteration),
            false
        );
    }

    #[rstest]
    fn test_service_is_on_without_control(
        simulation_time_iteration: SimulationTimeIteration,
        heat_battery: Arc<Mutex<HeatBattery>>,
    ) {
        let heat_battery_service: HeatBatteryServiceWaterRegular =
            HeatBatteryServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                TEMP_HOT_WATER,
                None,
            );

        assert_eq!(
            heat_battery_service.is_on(simulation_time_iteration),
            true
        );
    }

    // TODO test_demand_energy_for_water_regular
    // TODO test_demand_energy_service_off_for_water_regular
    // TODO test_energy_output_max_service_on_for_water_regular

    #[rstest]
    fn test_energy_output_max_service_off_for_water_regular(
        simulation_time_iteration: SimulationTimeIteration,
        heat_battery: Arc<Mutex<HeatBattery>>,
    ) {
        let heat_battery_service: HeatBatteryServiceWaterRegular =
            HeatBatteryServiceWaterRegular::new(
                heat_battery,
                SERVICE_NAME.into(),
                TEMP_HOT_WATER,
                None,
            );

        let temp_return = 40.;
        let result = heat_battery_service
            .energy_output_max(temp_return, simulation_time_iteration);

        assert_eq!(result, 0.);
    }

    #[rstest]
    fn test_temp_setpnt_for_space(
        simulation_time_iteration: SimulationTimeIteration,
        heat_battery: Arc<Mutex<HeatBattery>>,
    ) {
        let first_scheduled_temp = Some(21.);
        let ctrl: Control = create_setpoint_time_control(vec![first_scheduled_temp]);
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
        heat_battery: Arc<Mutex<HeatBattery>>,
    ) {
        let ctrl: Control = create_setpoint_time_control(vec![Some(21.)]);
        let heat_battery_space =
            HeatBatteryServiceSpace::new(heat_battery, SERVICE_NAME.into(), ctrl.into());

        assert_eq!(
            heat_battery_space.in_required_period(simulation_time_iteration),
            Some(true)
        );
    }

    // TODO test_demand_energy_for_space

    #[rstest]
    fn test_demand_energy_service_off_for_space(simulation_time_iteration: SimulationTimeIteration, heat_battery: Arc<Mutex<HeatBattery>>) {
        let energy_demand = 10.;
        let temp_return = 40.;
        let temp_flow = 1.;
        let ctrl: Control = create_setpoint_time_control(vec![None]);
        let heat_battery_space =
            HeatBatteryServiceSpace::new(heat_battery, SERVICE_NAME.into(), ctrl.into());
        let result = heat_battery_space.demand_energy(energy_demand, temp_return, temp_flow, simulation_time_iteration);
        assert_eq!(result, 0.);
    }
    
    // TODO in Rust we don't have energy_output_max on HeatBatteryServiceSpace
    // confirm this is expected. Following tests affected: 
    // test_energy_output_max_service_on_for_space
    // test_energy_output_max_service_off_for_space

    #[rstest]
    fn test_create_service_connection(heat_battery: Arc<Mutex<HeatBattery>>) {
        let service_name = "new_service";
        let result = HeatBattery::create_service_connection(heat_battery, service_name);

        // TODO Python tests that the service_name is added to the energy_supply_connections
        // but this is private in Rust
    }

}
