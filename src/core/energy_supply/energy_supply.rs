use crate::compare_floats::min_of_2;
use crate::core::energy_supply::elec_battery::ElectricBattery;
use crate::core::heating_systems::storage_tank::PVDiverter;
use crate::input::{EnergySupplyDetails, EnergySupplyInput, EnergySupplyType, HeatNetwork};
use crate::simulation_time::SimulationTimeIteration;
use anyhow::bail;
use indexmap::{indexmap, IndexMap};
use parking_lot::Mutex;
use std::sync::Arc;

#[derive(Debug)]
pub struct EnergySupplies {
    pub mains_electricity: Option<Arc<Mutex<EnergySupply>>>,
    pub mains_gas: Option<Arc<Mutex<EnergySupply>>>,
    pub bulk_lpg: Option<Arc<Mutex<EnergySupply>>>,
    pub bottled_lpg: Option<Arc<Mutex<EnergySupply>>>,
    pub condition_11f_lpg: Option<Arc<Mutex<EnergySupply>>>,
    pub custom: Option<Arc<Mutex<EnergySupply>>>,
    pub heat_network: Option<Arc<Mutex<EnergySupply>>>,
    pub unmet_demand: Arc<Mutex<EnergySupply>>,
}

impl EnergySupplies {
    pub fn calc_energy_import_export_betafactor(&mut self, simtime: SimulationTimeIteration) {
        if let Some(ref mut supply) = self.mains_electricity {
            supply.lock().calc_energy_import_export_betafactor(simtime);
        }
        if let Some(ref mut supply) = self.mains_gas {
            supply.lock().calc_energy_import_export_betafactor(simtime);
        }
        if let Some(ref mut supply) = self.bulk_lpg {
            supply.lock().calc_energy_import_export_betafactor(simtime);
        }
        self.unmet_demand
            .lock()
            .calc_energy_import_export_betafactor(simtime);
    }

    pub fn ensured_get_for_type(
        &mut self,
        energy_supply_type: EnergySupplyType,
        timesteps: usize,
    ) -> Arc<Mutex<EnergySupply>> {
        let energy_supply = match energy_supply_type {
            EnergySupplyType::Electricity => &mut self.mains_electricity,
            EnergySupplyType::MainsGas => &mut self.mains_gas,
            EnergySupplyType::UnmetDemand => return self.unmet_demand.clone(),
            EnergySupplyType::Custom => &mut self.custom,
            EnergySupplyType::LpgBulk => &mut self.bulk_lpg,
            EnergySupplyType::LpgBottled => &mut self.bottled_lpg,
            EnergySupplyType::LpgCondition11F => &mut self.condition_11f_lpg,
            EnergySupplyType::HeatNetwork => {
                unimplemented!("Undetermined what to do when heat network is requested.")
            }
        };
        match energy_supply {
            Some(supply) => supply.clone(),
            None => Arc::new(Mutex::new(EnergySupply::new(
                energy_supply_type,
                timesteps,
                None,
            ))),
        }
    }

    pub fn supplies_by_name(&self) -> IndexMap<&str, Arc<Mutex<EnergySupply>>> {
        let mut supplies: IndexMap<&str, Arc<Mutex<EnergySupply>>> = Default::default();
        supplies.insert("unmet_demand", self.unmet_demand.clone());
        if let Some(elec) = &self.mains_electricity {
            supplies.insert("electricity", elec.clone());
        }
        if let Some(gas) = &self.mains_gas {
            supplies.insert("mains gas", gas.clone());
        }
        if let Some(lpg) = &self.bulk_lpg {
            supplies.insert("LPG_bulk", lpg.clone());
        }
        if let Some(lpg) = &self.bottled_lpg {
            supplies.insert("LPG_bottled", lpg.clone());
        }
        if let Some(condition_11f) = &self.condition_11f_lpg {
            supplies.insert("LPG_condition_11F", condition_11f.clone());
        }
        if let Some(custom) = &self.custom {
            supplies.insert("custom", custom.clone());
        }
        if let Some(heat_network) = &self.heat_network {
            supplies.insert("heat network", heat_network.clone());
        }

        supplies
    }
}

/// An object to represent the connection of a system that consumes energy to the energy supply
///
/// This object encapsulates the name of the connection, meaning that the
/// system consuming the energy does not have to specify these on every call,
/// and helping to enforce that each connection to a single supply has a unique
/// name.
#[derive(Clone, Debug)]
pub struct EnergySupplyConnection {
    energy_supply: Arc<Mutex<EnergySupply>>,
    end_user_name: String,
}

impl EnergySupplyConnection {
    pub fn new(energy_supply: Arc<Mutex<EnergySupply>>, end_user_name: String) -> Self {
        Self {
            energy_supply,
            end_user_name,
        }
    }

    /// Forwards the amount of energy out (in kWh) to the relevant EnergySupply object
    pub fn energy_out(
        &self,
        amount_demanded: f64,
        timestep_idx: usize,
    ) -> Result<(), anyhow::Error> {
        self.energy_supply.lock().energy_out(
            self.end_user_name.as_str(),
            amount_demanded,
            timestep_idx,
        )
    }

    /// Forwards the amount of energy demanded (in kWh) to the relevant EnergySupply object
    pub fn demand_energy(
        &self,
        amount_demanded: f64,
        timestep_idx: usize,
    ) -> Result<(), anyhow::Error> {
        self.energy_supply.lock().demand_energy(
            self.end_user_name.as_str(),
            amount_demanded,
            timestep_idx,
        )
    }

    pub fn supply_energy(
        &self,
        amount_produced: f64,
        timestep_idx: usize,
    ) -> Result<(), anyhow::Error> {
        self.energy_supply.lock().supply_energy(
            self.end_user_name.as_str(),
            amount_produced,
            timestep_idx,
        )
    }

    pub fn fuel_type(&self) -> EnergySupplyType {
        self.energy_supply.lock().fuel_type()
    }
}

#[derive(Clone, Debug)]
pub struct EnergySupply {
    fuel_type: EnergySupplyType,
    simulation_timesteps: usize,
    electric_battery: Option<ElectricBattery>,
    diverter: Option<Arc<Mutex<PVDiverter>>>,
    demand_total: Vec<f64>,
    demand_by_end_user: IndexMap<String, Vec<f64>>,
    energy_out_by_end_user: IndexMap<String, Vec<f64>>,
    beta_factor: Vec<f64>,
    supply_surplus: Vec<f64>,
    demand_not_met: Vec<f64>,
    energy_into_battery: Vec<f64>,
    energy_out_of_battery: Vec<f64>,
    energy_diverted: Vec<f64>,
    energy_generated_consumed: Vec<f64>,
}

impl EnergySupply {
    /// Arguments:
    /// * `fuel_type` - string denoting type of fuel
    /// * `simulation_timesteps` - the number of steps in the simulation time being used
    /// * `electric_battery` - reference to an ElectricBattery object
    pub fn new(
        fuel_type: EnergySupplyType,
        simulation_timesteps: usize,
        electric_battery: Option<ElectricBattery>,
    ) -> Self {
        Self {
            fuel_type,
            simulation_timesteps,
            electric_battery,
            diverter: None,
            demand_total: init_demand_list(simulation_timesteps),
            demand_by_end_user: Default::default(),
            energy_out_by_end_user: Default::default(),
            beta_factor: init_demand_list(simulation_timesteps),
            supply_surplus: init_demand_list(simulation_timesteps),
            demand_not_met: init_demand_list(simulation_timesteps),
            energy_into_battery: init_demand_list(simulation_timesteps),
            energy_out_of_battery: init_demand_list(simulation_timesteps),
            energy_diverted: init_demand_list(simulation_timesteps),
            energy_generated_consumed: init_demand_list(simulation_timesteps),
        }
    }

    pub fn fuel_type(&self) -> EnergySupplyType {
        self.fuel_type
    }

    pub fn connection(
        energy_supply: Arc<Mutex<EnergySupply>>,
        end_user_name: &str,
    ) -> Result<EnergySupplyConnection, anyhow::Error> {
        let mut supply = energy_supply.lock();
        if supply.demand_by_end_user.contains_key(end_user_name) {
            bail!("The end user name '{end_user_name}' was already used.");
        }
        let timesteps = supply.simulation_timesteps;
        supply
            .demand_by_end_user
            .entry(end_user_name.into())
            .or_insert(init_demand_list(timesteps));
        supply
            .energy_out_by_end_user
            .entry(end_user_name.into())
            .or_insert(init_demand_list(timesteps));

        Ok(EnergySupplyConnection {
            energy_supply: energy_supply.clone(),
            end_user_name: end_user_name.to_string(),
        })
    }

    pub fn energy_out(
        &mut self,
        end_user_name: &str,
        amount_demanded: f64,
        timestep_index: usize,
    ) -> Result<(), anyhow::Error> {
        if !self.demand_by_end_user.contains_key(end_user_name) {
            bail!("Error: End user name not already registered by calling connection function.",);
        }
        self.energy_out_by_end_user.get_mut(end_user_name).unwrap()[timestep_index] +=
            amount_demanded;

        Ok(())
    }

    pub fn connect_diverter(
        &mut self,
        diverter: Arc<Mutex<PVDiverter>>,
    ) -> Result<(), &'static str> {
        if self.diverter.is_some() {
            return Err("diverter was already connected");
        }

        self.diverter = Some(diverter);

        Ok(())
    }

    /// This method is used in place of calling .connection() in the Python codebase in order to register an end user name
    pub fn register_end_user_name(&mut self, end_user_name: String) {
        self.demand_by_end_user.insert(
            end_user_name.clone(),
            init_demand_list(self.simulation_timesteps),
        );
        self.energy_out_by_end_user
            .insert(end_user_name, init_demand_list(self.simulation_timesteps));
    }

    pub fn demand_energy(
        &mut self,
        end_user_name: &str,
        amount_demanded: f64,
        timestep_index: usize,
    ) -> Result<(), anyhow::Error> {
        if !self.demand_by_end_user.contains_key(end_user_name) {
            bail!("Error: End user name not already registered by calling connection function.",);
        }
        self.demand_total[timestep_index] += amount_demanded;
        self.demand_by_end_user.get_mut(end_user_name).unwrap()[timestep_index] += amount_demanded;

        Ok(())
    }

    /// Record energy produced (in kWh) for the end user specified.
    ///
    /// Note: this is energy generated so it is subtracted from demand.
    /// Treat as negative
    pub fn supply_energy(
        &mut self,
        end_user_name: &str,
        amount_produced: f64,
        timestep_index: usize,
    ) -> Result<(), anyhow::Error> {
        self.demand_energy(end_user_name, amount_produced * -1.0, timestep_index)
    }

    /// Return list of the total demand on this energy source for each timestep
    pub fn results_total(&self) -> &Vec<f64> {
        &self.demand_total
    }

    /// Return the demand from each end user on this energy source for each timestep.
    ///
    /// Returns dictionary of lists, where dictionary keys are names of end users.
    pub fn results_by_end_user(&self) -> IndexMap<String, Vec<f64>> {
        if self
            .demand_by_end_user
            .keys()
            .cloned()
            .collect::<Vec<String>>()
            == self
                .energy_out_by_end_user
                .keys()
                .cloned()
                .collect::<Vec<String>>()
        {
            return self.demand_by_end_user.clone();
        }

        let mut all_results_by_end_user = indexmap! {};
        for (demand, energy_out) in self
            .demand_by_end_user
            .iter()
            .zip(self.energy_out_by_end_user.iter())
        {
            if demand.0 == energy_out.0 {
                let user_name = demand.0.clone(); // can use demand.0 or energy_out.0 to get end user name
                all_results_by_end_user.insert(
                    user_name,
                    demand
                        .1
                        .iter()
                        .enumerate()
                        .map(|(i, demand_val)| demand_val + energy_out.1[i])
                        .collect(),
                );
            }
        }

        all_results_by_end_user
    }

    pub fn get_energy_import(&self) -> &[f64] {
        &self.demand_not_met
    }

    pub fn get_energy_export(&self) -> &[f64] {
        &self.supply_surplus
    }

    /// Return the amount of generated energy consumed in the building for all timesteps
    pub fn get_energy_generated_consumed(&self) -> &[f64] {
        &self.energy_generated_consumed
    }

    /// Return the amount of generated energy sent to battery and drawn from battery
    pub fn get_energy_to_from_battery(&self) -> (&[f64], &[f64]) {
        (&self.energy_into_battery, &self.energy_out_of_battery)
    }

    /// Return the amount of generated energy diverted to minimise export
    pub fn get_energy_diverted(&self) -> &[f64] {
        &self.energy_diverted
    }

    pub fn get_beta_factor(&self) -> &[f64] {
        &self.beta_factor
    }

    /// Calculate how much of that supply can be offset against demand.
    /// And then calculate what demand and supply is left after offsetting, which are the amount exported imported
    pub fn calc_energy_import_export_betafactor(&mut self, simtime: SimulationTimeIteration) {
        let mut supplies = vec![];
        let mut demands = vec![];
        let timestep_idx = simtime.index;
        for user in self.demand_by_end_user.keys() {
            let demand = self.demand_by_end_user[user][timestep_idx];
            // if energy is negative that means it's actually a supply, we
            // need to separate the two for beta factor calc. If we had
            // multiple different supplies they would have to be separated
            // here
            if demand < 0. {
                supplies.push(demand);
            } else {
                demands.push(demand);
            }
        }

        *self.beta_factor.get_mut(timestep_idx).unwrap() = self.beta_factor_function(
            -supplies.iter().sum::<f64>(),
            demands.iter().sum::<f64>(),
            BetaFactorFunction::Pv,
        );

        let current_beta_factor = self.beta_factor[timestep_idx];
        let supplies_sum = supplies.iter().sum::<f64>();
        // PV elec consumed within dwelling in absence of battery storage or diverter (kWh)
        // if there were multiple sources they would each have their own beta factors
        let supply_consumed = supplies_sum * current_beta_factor;
        // Surplus PV elec generation (kWh) - ie amount to be exported to the grid or batteries
        let mut supply_surplus = supplies_sum * (1. - current_beta_factor);
        // Elec demand not met by PV (kWh) - ie amount to be imported from the grid or batteries
        let mut demand_not_met = demands.iter().sum::<f64>() + supply_consumed;
        // See if there is a net supply/demand for the timestep
        if let Some(ref mut battery) = &mut self.electric_battery {
            // See if the battery can deal with excess supply/demand for this timestep
            // supply_surplus is -ve by convention and demand_not_met is +ve
            let energy_out_of_battery = battery.charge_discharge_battery(supply_surplus);
            supply_surplus -= energy_out_of_battery;
            *self.energy_into_battery.get_mut(timestep_idx).unwrap() = -energy_out_of_battery;
            let energy_out_of_battery = battery.charge_discharge_battery(demand_not_met);
            demand_not_met -= energy_out_of_battery;
            *self.energy_out_of_battery.get_mut(timestep_idx).unwrap() = -energy_out_of_battery;
        }

        if let Some(ref mut diverter) = &mut self.diverter {
            *self.energy_diverted.get_mut(timestep_idx).unwrap() =
                diverter.lock().divert_surplus(supply_surplus, simtime);
            supply_surplus += self.energy_diverted[timestep_idx];
        }

        *self.supply_surplus.get_mut(timestep_idx).unwrap() += supply_surplus;
        *self.demand_not_met.get_mut(timestep_idx).unwrap() += demand_not_met;
        // Report energy generated and consumed as positive number, so subtract negative number
        *self
            .energy_generated_consumed
            .get_mut(timestep_idx)
            .unwrap() -= supply_consumed;
    }

    /// wrapper that applies relevant function to obtain
    /// beta factor from energy supply+demand at a given timestep
    fn beta_factor_function(
        &self,
        supply: f64,
        demand: f64,
        beta_factor_function: BetaFactorFunction,
    ) -> f64 {
        if supply == 0. {
            return 1.;
        }
        if demand == 0. {
            return 0.;
        }

        let demand_ratio = supply / demand;
        let beta_factor = match beta_factor_function {
            BetaFactorFunction::Pv => min_of_2(0.6748 * demand_ratio.powf(-0.703), 1.),
            BetaFactorFunction::Wind => {
                unimplemented!("Wind beta factor function is not implemented")
            } // wind is mentioned in Python but currently commented out
        };

        min_of_2(beta_factor, 1. / demand_ratio)
    }
}

enum BetaFactorFunction {
    Pv,
    Wind,
}

fn init_demand_list(timestep_count: usize) -> Vec<f64> {
    vec![0.0; timestep_count]
}

pub fn from_input(input: EnergySupplyInput, simulation_timesteps: usize) -> EnergySupplies {
    EnergySupplies {
        mains_electricity: input
            .mains_electricity
            .map(|s| supply_from_details(s, simulation_timesteps)),
        mains_gas: input
            .mains_gas
            .map(|s| supply_from_details(s, simulation_timesteps)),
        bulk_lpg: input
            .bulk_lpg
            .map(|s| supply_from_details(s, simulation_timesteps)),
        heat_network: input
            .heat_network
            .map(|hn| supply_from_heat_network_details(hn, simulation_timesteps)),
        unmet_demand: Arc::new(Mutex::new(EnergySupply::new(
            EnergySupplyType::UnmetDemand,
            simulation_timesteps,
            Default::default(),
        ))),
        custom: None,
        condition_11f_lpg: None,
        bottled_lpg: None,
    }
}

fn supply_from_details(
    energy_supply_details: EnergySupplyDetails,
    simulation_timesteps: usize,
) -> Arc<Mutex<EnergySupply>> {
    Arc::new(Mutex::new(EnergySupply::new(
        energy_supply_details.fuel,
        simulation_timesteps,
        energy_supply_details
            .electric_battery
            .map(ElectricBattery::from_input),
    )))
}

fn supply_from_heat_network_details(
    heat_network: HeatNetwork,
    simulation_timesteps: usize,
) -> Arc<Mutex<EnergySupply>> {
    Arc::new(Mutex::new(EnergySupply::new(
        heat_network.fuel,
        simulation_timesteps,
        None,
    )))
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::input::EnergySupplyType::MainsGas;
    use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
    use rstest::*;

    #[fixture]
    pub fn simulation_time() -> SimulationTimeIterator {
        SimulationTime::new(0.0, 8.0, 1.0).iter()
    }

    #[fixture]
    pub fn energy_supply<'a>(simulation_time: SimulationTimeIterator) -> EnergySupply {
        let mut energy_supply = EnergySupply::new(MainsGas, simulation_time.total_steps(), None);
        energy_supply.register_end_user_name("shower".to_string());
        energy_supply.register_end_user_name("bath".to_string());

        energy_supply
    }

    #[fixture]
    pub fn energy_supply_connections(
        energy_supply: EnergySupply,
    ) -> (
        EnergySupplyConnection,
        EnergySupplyConnection,
        Arc<Mutex<EnergySupply>>,
    ) {
        let shared_supply = Arc::new(Mutex::new(energy_supply));
        let energy_connection_1 = EnergySupplyConnection {
            energy_supply: shared_supply.clone(),
            end_user_name: "shower".to_string(),
        };
        let energy_connection_2 = EnergySupplyConnection {
            energy_supply: shared_supply.clone(),
            end_user_name: "bath".to_string(),
        };
        (energy_connection_1, energy_connection_2, shared_supply)
    }

    #[fixture]
    pub fn energy_supply_connection_1<'a>(energy_supply: EnergySupply) -> EnergySupplyConnection {
        EnergySupplyConnection {
            energy_supply: Arc::new(Mutex::new(energy_supply)),
            end_user_name: "shower".to_string(),
        }
    }

    #[fixture]
    pub fn energy_supply_connection_2<'a>(energy_supply: EnergySupply) -> EnergySupplyConnection {
        EnergySupplyConnection {
            energy_supply: Arc::new(Mutex::new(energy_supply)),
            end_user_name: "bath".to_string(),
        }
    }

    const EXPECTED_TOTAL_DEMANDS: [f64; 8] =
        [50.0, 120.0, 190.0, 260.0, 330.0, 400.0, 470.0, 540.0];

    #[rstest]
    pub fn test_results_total(
        mut energy_supply: EnergySupply,
        simulation_time: SimulationTimeIterator,
    ) {
        for simtime in simulation_time {
            let _ = energy_supply.demand_energy(
                "shower",
                (simtime.index as f64 + 1.0) * 50.0,
                simtime.index,
            );
            let _ = energy_supply.demand_energy("bath", simtime.index as f64 * 20.0, simtime.index);
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

    #[rstest]
    pub fn test_results_by_end_user(
        energy_supply_connections: (
            EnergySupplyConnection,
            EnergySupplyConnection,
            Arc<Mutex<EnergySupply>>,
        ),
        simulation_time: SimulationTimeIterator,
    ) {
        let (energy_connection_1, energy_connection_2, energy_supply) = energy_supply_connections;
        for simtime in simulation_time {
            let _ = energy_connection_1
                .demand_energy((simtime.index as f64 + 1.0) * 50.0, simtime.index);
            let _ = energy_connection_2.demand_energy(simtime.index as f64 * 20.0, simtime.index);
            assert_eq!(
                energy_supply.lock().results_by_end_user()["shower"][simtime.index],
                EXPECTED_TOTAL_DEMANDS_BY_END_USER[0][simtime.index]
            );
            assert_eq!(
                energy_supply.lock().results_by_end_user()["bath"][simtime.index],
                EXPECTED_TOTAL_DEMANDS_BY_END_USER[1][simtime.index]
            );
        }
    }

    const EXPECTED_BETA_FACTORS: [f64; 8] = [
        1.0,
        0.8973610789278808,
        0.4677549807236648,
        0.3297589507351858,
        0.2578125,
        0.2,
        0.16319444444444445,
        0.1377551020408163,
    ];
    const EXPECTED_SURPLUSES: [f64; 8] = [
        0.0,
        -8.21111368576954,
        -170.3184061684273,
        -482.57355547066624,
        -950.0,
        -1600.0,
        -2410.0,
        -3380.0,
    ];
    const EXPECTED_DEMANDS_NOT_MET: [f64; 8] = [
        50.0,
        48.21111368576953,
        40.31840616842726,
        22.573555470666236,
        0.0,
        0.0,
        0.0,
        0.0,
    ];

    #[rstest]
    pub fn test_beta_factor(
        energy_supply_connections: (
            EnergySupplyConnection,
            EnergySupplyConnection,
            Arc<Mutex<EnergySupply>>,
        ),
        simulation_time: SimulationTimeIterator,
    ) {
        let (energy_connection_1, energy_connection_2, energy_supply) = energy_supply_connections;
        let energy_connection_3 = EnergySupply::connection(energy_supply.clone(), "PV").unwrap();
        for (t_idx, t_it) in simulation_time.enumerate() {
            energy_connection_1
                .demand_energy((t_idx as f64 + 1.) * 50., t_idx)
                .unwrap();
            energy_connection_2
                .demand_energy(t_idx as f64 * 20., t_idx)
                .unwrap();
            energy_connection_3
                .supply_energy(t_idx as f64 * t_idx as f64 * 80., t_idx)
                .unwrap();

            let mut energy_supply = energy_supply.lock();
            energy_supply.calc_energy_import_export_betafactor(t_it);

            assert_eq!(
                energy_supply.get_beta_factor()[t_idx],
                EXPECTED_BETA_FACTORS[t_idx],
                "incorrect beta factor returned"
            );
            assert_eq!(
                energy_supply.get_energy_export()[t_idx],
                EXPECTED_SURPLUSES[t_idx],
                "incorrect energy export returned"
            );
            assert_eq!(
                energy_supply.get_energy_import()[t_idx],
                EXPECTED_DEMANDS_NOT_MET[t_idx],
                "incorrect energy import returned"
            );
        }
    }
}
