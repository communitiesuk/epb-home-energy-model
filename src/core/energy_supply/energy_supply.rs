use crate::input::{ElectricBattery, EnergyDiverter, EnergySupplyType};
use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
use std::collections::HashMap;
use std::error::Error;

pub struct EnergySupply {
    fuel_type: EnergySupplyType,
    simulation_time: SimulationTimeIterator,
    electric_battery: Option<ElectricBattery>,
    diverter: Option<EnergyDiverter>,
    demand_total: Vec<f64>,
    demand_by_end_user: HashMap<String, Vec<f64>>,
    energy_out_by_end_user: HashMap<String, Vec<f64>>,
    beta_factor: Vec<f64>,
    supply_surplus: Vec<f64>,
    demand_not_met: Vec<f64>,
    energy_generated_consumed: Vec<f64>,
}

impl EnergySupply {
    /// Arguments:
    /// * `fuel_type` - string denoting type of fuel
    /// * `simulation_time` - reference to SimulationTimeIterator object
    /// * `electric_battery` - reference to an ElectricBattery object
    pub fn new(
        fuel_type: EnergySupplyType,
        simulation_time: SimulationTimeIterator,
        electric_battery: Option<ElectricBattery>,
    ) -> Self {
        let total_steps = &simulation_time.total_steps();
        Self {
            fuel_type,
            simulation_time,
            electric_battery,
            diverter: None,
            demand_total: init_demand_list(total_steps),
            demand_by_end_user: Default::default(),
            energy_out_by_end_user: Default::default(),
            beta_factor: init_demand_list(total_steps),
            supply_surplus: init_demand_list(total_steps),
            demand_not_met: init_demand_list(total_steps),
            energy_generated_consumed: init_demand_list(total_steps),
        }
    }

    pub fn fuel_type(&self) -> EnergySupplyType {
        self.fuel_type
    }

    pub fn energy_out(
        &mut self,
        end_user_name: String,
        amount_demanded: f64,
        timestep_index: usize,
    ) -> Result<(), &'static str> {
        if !self.demand_by_end_user.contains_key(&end_user_name) {
            return Err(
                "Error: End user name not already registered by calling connection function.",
            );
        }
        self.energy_out_by_end_user.get_mut(&end_user_name).unwrap()[timestep_index] +=
            amount_demanded;

        Ok(())
    }

    pub fn connect_diverter(&mut self, diverter: EnergyDiverter) -> Result<(), &'static str> {
        if self.diverter.is_some() {
            return Err("diverter was already connected");
        }

        self.diverter = Some(diverter);

        Ok(())
    }

    /// This method is used in place of calling .connection() in the Python codebase in order to register an end user name
    pub fn register_end_user_name(&mut self, end_user_name: String) -> () {
        self.demand_by_end_user.insert(
            end_user_name.clone(),
            init_demand_list(&self.simulation_time.total_steps()),
        );
        self.energy_out_by_end_user.insert(
            end_user_name,
            init_demand_list(&self.simulation_time.total_steps()),
        );
    }

    pub fn demand_energy(
        &mut self,
        end_user_name: String,
        amount_demanded: f64,
        timestep_index: usize,
    ) -> Result<(), &'static str> {
        if !self.demand_by_end_user.contains_key(&end_user_name) {
            return Err(
                "Error: End user name not already registered by calling connection function.",
            );
        }
        self.demand_total[timestep_index] += amount_demanded;
        self.demand_by_end_user.get_mut(&end_user_name).unwrap()[timestep_index] += amount_demanded;

        Ok(())
    }

    /// Record energy produced (in kWh) for the end user specified.
    ///
    /// Note: this is energy generated so it is subtracted from demand.
    /// Treat as negative
    pub fn supply_energy(
        &mut self,
        end_user_name: String,
        amount_produced: f64,
        timestep_index: usize,
    ) -> Result<(), &'static str> {
        self.demand_energy(end_user_name, amount_produced * -1.0, timestep_index)
    }

    /// Return list of the total demand on this energy source for each timestep
    pub fn results_total(&self) -> &Vec<f64> {
        &self.demand_total
    }

    // Return the demand from each end user on this energy source for each timestep.
    //
    // Returns dictionary of lists, where dictionary keys are names of end users.
    // pub fn results_by_end_user(&self) -> HashMap<String, Vec<f64>> {
    //     if self.demand_by_end_user.keys().collect::<Vec<&str>>()
    //         == self.energy_out_by_end_user.keys().collect::<Vec<&str>>()
    //     {
    //         return self.demand_by_end_user;
    //     }
    //
    //     HashMap::from([])
    // }
}

fn init_demand_list(timestep_count: &i32) -> Vec<f64> {
    vec![0.0; *timestep_count as usize]
}

// EnergySupplyConnection is a delegating object in the Python codebase - working round implementing it
// in Rust as it seems non-essential and involves some unnecessary fighting with the borrow checker
//
// pub struct EnergySupplyConnection<'a> {
//     energy_supply: &'a EnergySupply,
//     end_user_name: String,
// }
//
// impl<'a> EnergySupplyConnection<'a> {
//     pub fn energy_out(&self, amount_demanded: f64) -> () {}
//
//     pub fn demand_energy(&self, amount_demanded: f64) -> () {}
//
//     pub fn supply_energy(&self, amount_demanded: f64) -> () {}
// }

#[cfg(test)]
mod test {
    use super::*;
    use crate::input::EnergySupplyType::MainsGas;
    use crate::simulation_time::SimulationTime;
    use rstest::*;

    #[fixture]
    pub fn simulation_time() -> SimulationTimeIterator {
        SimulationTime::new(0.0, 8.0, 1.0).iter()
    }

    #[fixture]
    pub fn energy_supply(simulation_time: SimulationTimeIterator) -> EnergySupply {
        let mut energy_supply = EnergySupply::new(MainsGas, simulation_time, None);
        energy_supply.register_end_user_name("shower".to_string());
        energy_supply.register_end_user_name("bath".to_string());

        energy_supply
    }

    // #[fixture]
    // pub fn energy_supply_connection_1<'a>(
    //     energy_supply: &'a EnergySupply,
    // ) -> EnergySupplyConnection<'a> {
    //     EnergySupplyConnection {
    //         energy_supply,
    //         end_user_name: "shower".to_string(),
    //     }
    // }
    //
    // #[fixture]
    // pub fn energy_supply_connection_2<'a>(
    //     energy_supply: &'a EnergySupply,
    // ) -> EnergySupplyConnection<'a> {
    //     EnergySupplyConnection {
    //         energy_supply,
    //         end_user_name: "bath".to_string(),
    //     }
    // }

    const EXPECTED_TOTAL_DEMANDS: [f64; 8] =
        [50.0, 120.0, 190.0, 260.0, 330.0, 400.0, 470.0, 540.0];

    #[rstest]
    pub fn should_have_correct_results_total(
        mut energy_supply: EnergySupply,
        simulation_time: SimulationTimeIterator,
    ) {
        for simtime in simulation_time {
            energy_supply.demand_energy(
                "shower".to_string(),
                (simtime.index as f64 + 1.0) * 50.0,
                simtime.index,
            );
            energy_supply.demand_energy(
                "bath".to_string(),
                simtime.index as f64 * 20.0,
                simtime.index,
            );
            assert_eq!(
                energy_supply.results_total()[simtime.index],
                EXPECTED_TOTAL_DEMANDS[simtime.index],
                "incorrect total demand energy returned on iteration {} (1-indexed)",
                simtime.index + 1
            )
        }
    }

    const EXPECTED_TOTAL_DEMANDS_BY_END_USER: [[f64; 8]; 2] = [
        [50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0],
        [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0],
    ];

    // #[rstest]
    // pub fn should_have_correct_results_by_end_user(
    //     mut energy_supply: EnergySupply,
    //     simulation_time: SimulationTime,
    // ) {
    //     for simtime in simulation_time.iter() {
    //         energy_supply.demand_energy("shower".to_string(), (simtime.index as f64 + 1.0)*50.0, simtime.index);
    //         energy_supply.demand_energy("bath".to_string(), simtime.index as f64 * 20.0, simtime.index);
    //         assert_eq!(
    //             energy_supply.results_by_end_user
    //         );
    //     }
    // }
}
