use crate::compare_floats::min_of_2;
use crate::core::energy_supply::elec_battery::ElectricBattery;
use crate::core::heating_systems::storage_tank::PVDiverter;
use crate::errors::NotImplementedError;
use crate::input::{EnergySupplyType, FuelType, SecondarySupplyType};
use crate::simulation_time::SimulationTimeIteration;
use anyhow::bail;
use atomic_float::AtomicF64;
use indexmap::{indexmap, IndexMap};
use parking_lot::RwLock;
use std::sync::atomic::Ordering;
use std::sync::Arc;

pub(crate) const UNMET_DEMAND_SUPPLY_NAME: &str = "_unmet_demand";

#[derive(Debug)]
pub struct EnergySupplies {
    pub mains_electricity: Option<Arc<RwLock<EnergySupply>>>,
    pub mains_gas: Option<Arc<RwLock<EnergySupply>>>,
    pub bulk_lpg: Option<Arc<RwLock<EnergySupply>>>,
    pub bottled_lpg: Option<Arc<RwLock<EnergySupply>>>,
    pub condition_11f_lpg: Option<Arc<RwLock<EnergySupply>>>,
    pub custom: Option<Arc<RwLock<EnergySupply>>>,
    pub heat_network: Option<Arc<RwLock<EnergySupply>>>,
    pub unmet_demand: Arc<RwLock<EnergySupply>>,
}

impl EnergySupplies {
    pub fn calc_energy_import_export_betafactor(
        &self,
        simtime: SimulationTimeIteration,
    ) -> Result<(), NotImplementedError> {
        if let Some(ref supply) = self.mains_electricity {
            supply
                .read()
                .calc_energy_import_export_betafactor(simtime)?;
        }
        if let Some(ref supply) = self.mains_gas {
            supply
                .read()
                .calc_energy_import_export_betafactor(simtime)?;
        }
        if let Some(ref supply) = self.bulk_lpg {
            supply
                .read()
                .calc_energy_import_export_betafactor(simtime)?;
        }
        self.unmet_demand
            .read()
            .calc_energy_import_export_betafactor(simtime)?;

        Ok(())
    }

    pub fn ensured_get_for_type(
        &mut self,
        energy_supply_type: EnergySupplyType,
        timesteps: usize,
    ) -> anyhow::Result<Arc<RwLock<EnergySupply>>> {
        let energy_supply = match energy_supply_type {
            EnergySupplyType::Electricity => &mut self.mains_electricity,
            EnergySupplyType::MainsGas => &mut self.mains_gas,
            EnergySupplyType::UnmetDemand => return Ok(self.unmet_demand.clone()),
            EnergySupplyType::Custom => &mut self.custom,
            EnergySupplyType::LpgBulk => &mut self.bulk_lpg,
            EnergySupplyType::LpgBottled => &mut self.bottled_lpg,
            EnergySupplyType::LpgCondition11F => &mut self.condition_11f_lpg,
            EnergySupplyType::HeatNetwork => &mut self.heat_network,
            #[cfg(feature = "fhs")]
            EnergySupplyType::NotionalHeatNetwork => &mut None, // nothing seems to request this, so match with nothing
        };
        match energy_supply {
            Some(supply) => Ok(supply.clone()),
            None => Ok(Arc::new(RwLock::new(EnergySupply::new(
                energy_supply_type.try_into()?,
                timesteps,
                None,
                None,
                None,
            )))),
        }
    }

    pub fn supplies_by_name(&self) -> IndexMap<&str, Arc<RwLock<EnergySupply>>> {
        let mut supplies: IndexMap<&str, Arc<RwLock<EnergySupply>>> = Default::default();
        supplies.insert("_unmet_demand", self.unmet_demand.clone());
        if let Some(elec) = &self.mains_electricity {
            supplies.insert("mains elec", elec.clone());
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

impl Default for EnergySupplies {
    fn default() -> Self {
        Self {
            mains_electricity: None,
            mains_gas: None,
            bulk_lpg: None,
            bottled_lpg: None,
            condition_11f_lpg: None,
            custom: None,
            heat_network: None,
            unmet_demand: Arc::new(RwLock::new(EnergySupply::new(
                FuelType::UnmetDemand,
                0,
                None,
                None,
                None,
            ))),
        }
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
    energy_supply: Arc<RwLock<EnergySupply>>,
    pub(crate) end_user_name: String,
}

impl EnergySupplyConnection {
    pub fn new(energy_supply: Arc<RwLock<EnergySupply>>, end_user_name: String) -> Self {
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
        self.energy_supply.read().energy_out(
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
        self.energy_supply.read().demand_energy(
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
        self.energy_supply.read().supply_energy(
            self.end_user_name.as_str(),
            amount_produced,
            timestep_idx,
        )
    }

    pub fn fuel_type(&self) -> FuelType {
        self.energy_supply.read().fuel_type()
    }
}

#[derive(Debug)]
pub struct EnergySupply {
    fuel_type: FuelType,
    simulation_timesteps: usize,
    electric_battery: Option<ElectricBattery>,
    diverter: Option<Arc<RwLock<PVDiverter>>>,
    priority: Option<Vec<SecondarySupplyType>>,
    is_export_capable: bool,
    demand_total: Vec<AtomicF64>,
    demand_by_end_user: IndexMap<String, Vec<AtomicF64>>,
    energy_out_by_end_user: IndexMap<String, Vec<AtomicF64>>,
    beta_factor: Vec<AtomicF64>,
    supply_surplus: Vec<AtomicF64>,
    demand_not_met: Vec<AtomicF64>,
    energy_into_battery: Vec<AtomicF64>,
    energy_out_of_battery: Vec<AtomicF64>,
    energy_diverted: Vec<AtomicF64>,
    energy_generated_consumed: Vec<AtomicF64>,
}

impl EnergySupply {
    /// Arguments:
    /// * `fuel_type` - string denoting type of fuel
    /// * `simulation_timesteps` - the number of steps in the simulation time being used
    /// * `electric_battery` - reference to an ElectricBattery object
    /// * `priority`
    /// * `is_export_capable` - denotes that this Energy Supply can export its surplus supply
    pub fn new(
        fuel_type: FuelType,
        simulation_timesteps: usize,
        electric_battery: Option<ElectricBattery>,
        priority: Option<Vec<SecondarySupplyType>>,
        is_export_capable: Option<bool>,
    ) -> Self {
        Self {
            fuel_type,
            simulation_timesteps,
            electric_battery,
            diverter: None,
            priority,
            is_export_capable: is_export_capable.unwrap_or(true),
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

    pub fn fuel_type(&self) -> FuelType {
        self.fuel_type
    }

    pub fn connection(
        energy_supply: Arc<RwLock<EnergySupply>>,
        end_user_name: &str,
    ) -> Result<EnergySupplyConnection, anyhow::Error> {
        let mut supply = energy_supply.write();
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
        &self,
        end_user_name: &str,
        amount_demanded: f64,
        timestep_index: usize,
    ) -> Result<(), anyhow::Error> {
        if !self.demand_by_end_user.contains_key(end_user_name) {
            bail!("Error: End user name not already registered by calling connection function.",);
        }
        self.energy_out_by_end_user.get(end_user_name).unwrap()[timestep_index]
            .fetch_add(amount_demanded, Ordering::SeqCst);

        Ok(())
    }

    pub fn connect_diverter(
        &mut self,
        diverter: Arc<RwLock<PVDiverter>>,
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
        &self,
        end_user_name: &str,
        amount_demanded: f64,
        timestep_index: usize,
    ) -> Result<(), anyhow::Error> {
        if !self.demand_by_end_user.contains_key(end_user_name) {
            bail!("Error: End user name not already registered by calling connection function.",);
        }
        self.demand_total[timestep_index].fetch_add(amount_demanded, Ordering::SeqCst);
        self.demand_by_end_user.get(end_user_name).unwrap()[timestep_index]
            .fetch_add(amount_demanded, Ordering::SeqCst);

        Ok(())
    }

    /// Record energy produced (in kWh) for the end user specified.
    ///
    /// Note: this is energy generated so it is subtracted from demand.
    /// Treat as negative
    pub fn supply_energy(
        &self,
        end_user_name: &str,
        amount_produced: f64,
        timestep_index: usize,
    ) -> Result<(), anyhow::Error> {
        self.demand_energy(end_user_name, amount_produced * -1.0, timestep_index)
    }

    /// Return list of the total demand on this energy source for each timestep
    pub fn results_total(&self) -> Vec<f64> {
        self.demand_total
            .iter()
            .map(|d| d.load(Ordering::SeqCst))
            .collect()
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
            return self
                .demand_by_end_user
                .iter()
                .map(|(end_user, demand)| {
                    (
                        end_user.clone(),
                        demand.iter().map(|d| d.load(Ordering::SeqCst)).collect(),
                    )
                })
                .collect();
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
                        .map(|(i, demand_val)| {
                            demand_val.load(Ordering::SeqCst)
                                + energy_out.1[i].load(Ordering::SeqCst)
                        })
                        .collect(),
                );
            }
        }

        all_results_by_end_user
    }

    pub fn get_energy_import(&self) -> Vec<f64> {
        Self::vec_of_floats_from_atomics(&self.demand_not_met)
    }

    pub fn get_energy_export(&self) -> Vec<f64> {
        Self::vec_of_floats_from_atomics(&self.supply_surplus)
    }

    /// Return the amount of generated energy consumed in the building for all timesteps
    pub fn get_energy_generated_consumed(&self) -> Vec<f64> {
        Self::vec_of_floats_from_atomics(&self.energy_generated_consumed)
    }

    /// Return the amount of generated energy sent to battery and drawn from battery
    pub fn get_energy_to_from_battery(&self) -> (Vec<f64>, Vec<f64>) {
        (
            Self::vec_of_floats_from_atomics(&self.energy_into_battery),
            Self::vec_of_floats_from_atomics(&self.energy_out_of_battery),
        )
    }

    /// Return the amount of generated energy diverted to minimise export
    pub fn get_energy_diverted(&self) -> Vec<f64> {
        Self::vec_of_floats_from_atomics(&self.energy_diverted)
    }

    pub fn get_beta_factor(&self) -> Vec<f64> {
        Self::vec_of_floats_from_atomics(&self.beta_factor)
    }

    fn vec_of_floats_from_atomics(atomics: &[AtomicF64]) -> Vec<f64> {
        atomics
            .iter()
            .map(|v| v.load(Ordering::SeqCst))
            .collect::<Vec<_>>()
    }

    /// Calculate how much of that supply can be offset against demand.
    /// And then calculate what demand and supply is left after offsetting, which are the amount exported imported
    pub fn calc_energy_import_export_betafactor(
        &self,
        simtime: SimulationTimeIteration,
    ) -> Result<(), NotImplementedError> {
        let end_user_count = self.demand_by_end_user.len();
        let mut supplies = Vec::with_capacity(end_user_count);
        let mut demands = Vec::with_capacity(end_user_count);
        let timestep_idx = simtime.index;
        for user in self.demand_by_end_user.keys() {
            let demand = self.demand_by_end_user[user].get(timestep_idx).unwrap();
            // if energy is negative that means it's actually a supply, we
            // need to separate the two for beta factor calc. If we had
            // multiple different supplies they would have to be separated
            // here
            if demand.load(Ordering::SeqCst) < 0. {
                supplies.push(demand);
            } else {
                demands.push(demand);
            }
        }

        let supplies_sum = supplies
            .iter()
            .map(|d| d.load(Ordering::SeqCst))
            .sum::<f64>();
        let demands_sum = demands
            .iter()
            .map(|d| d.load(Ordering::SeqCst))
            .sum::<f64>();

        self.beta_factor.get(timestep_idx).unwrap().store(
            self.beta_factor_function(-supplies_sum, demands_sum, BetaFactorFunction::Pv)?,
            Ordering::SeqCst,
        );

        let current_beta_factor = self.beta_factor[timestep_idx].load(Ordering::SeqCst);

        // PV elec consumed within dwelling in absence of battery storage or diverter (kWh)
        // if there were multiple sources they would each have their own beta factors
        let supply_consumed = supplies_sum * current_beta_factor;
        // Surplus PV elec generation (kWh) - ie amount to be exported to the grid or batteries
        let mut supply_surplus = supplies_sum * (1. - current_beta_factor);
        // Elec demand not met by PV (kWh) - ie amount to be imported from the grid or batteries
        let mut demand_not_met = demands_sum + supply_consumed;

        match &self.priority {
            None => {
                if let Some(ref battery) = &self.electric_battery {
                    // See if the battery can deal with excess supply/demand for this timestep
                    // supply_surplus is -ve by convention and demand_not_met is +ve
                    let energy_out_of_battery =
                        battery.charge_discharge_battery(supply_surplus, simtime);
                    supply_surplus -= energy_out_of_battery;
                    self.energy_into_battery
                        .get(timestep_idx)
                        .unwrap()
                        .store(-energy_out_of_battery, Ordering::SeqCst);
                    let energy_out_of_battery =
                        battery.charge_discharge_battery(demand_not_met, simtime);
                    demand_not_met -= energy_out_of_battery;
                    self.energy_out_of_battery
                        .get(timestep_idx)
                        .unwrap()
                        .store(-energy_out_of_battery, Ordering::SeqCst);
                }

                if let Some(ref diverter) = &self.diverter {
                    self.energy_diverted.get(timestep_idx).unwrap().store(
                        diverter.read().divert_surplus(supply_surplus, simtime),
                        Ordering::SeqCst,
                    );
                    supply_surplus += self.energy_diverted[timestep_idx].load(Ordering::SeqCst);
                }
            }
            Some(priority) => {
                for item in priority {
                    if matches!(item, SecondarySupplyType::ElectricBattery)
                        && self.electric_battery.is_some()
                    {
                        let electric_battery = self.electric_battery.as_ref().unwrap();
                        let energy_out_of_battery =
                            electric_battery.charge_discharge_battery(supply_surplus, simtime);
                        supply_surplus -= energy_out_of_battery;
                        self.energy_into_battery[simtime.index]
                            .store(-energy_out_of_battery, Ordering::SeqCst);
                        let energy_out_of_battery =
                            electric_battery.charge_discharge_battery(demand_not_met, simtime);
                        demand_not_met -= energy_out_of_battery;
                        self.energy_out_of_battery[simtime.index]
                            .store(-energy_out_of_battery, Ordering::SeqCst);
                    } else if matches!(item, SecondarySupplyType::Diverter)
                        && self.diverter.is_some()
                    {
                        let diverter = self.diverter.as_ref().unwrap();
                        self.energy_diverted[simtime.index].store(
                            diverter.read().divert_surplus(supply_surplus, simtime),
                            Ordering::SeqCst,
                        );
                        supply_surplus +=
                            self.energy_diverted[simtime.index].load(Ordering::SeqCst);
                    }
                }
            }
        }

        if self.is_export_capable {
            self.supply_surplus
                .get(timestep_idx)
                .unwrap()
                .fetch_add(supply_surplus, Ordering::SeqCst);
        }

        self.demand_not_met
            .get(timestep_idx)
            .unwrap()
            .fetch_add(demand_not_met, Ordering::SeqCst);
        // Report energy generated and consumed as positive number, so subtract negative number
        self.energy_generated_consumed
            .get(timestep_idx)
            .unwrap()
            .fetch_sub(supply_consumed, Ordering::SeqCst);

        Ok(())
    }

    /// wrapper that applies relevant function to obtain
    /// beta factor from energy supply+demand at a given timestep
    fn beta_factor_function(
        &self,
        supply: f64,
        demand: f64,
        beta_factor_function: BetaFactorFunction,
    ) -> Result<f64, NotImplementedError> {
        if supply == 0. {
            return Ok(1.);
        }
        if demand == 0. {
            return Ok(0.);
        }

        let demand_ratio = supply / demand;
        let beta_factor = match beta_factor_function {
            BetaFactorFunction::Pv => min_of_2(0.6748 * demand_ratio.powf(-0.703), 1.),
            BetaFactorFunction::Wind => {
                return Err(NotImplementedError::new(
                    "Wind beta factor function is not implemented upstream",
                ));
            } // wind is mentioned in Python but currently commented out
        };

        Ok(min_of_2(beta_factor, 1. / demand_ratio))
    }
}

enum BetaFactorFunction {
    Pv,
    // variant currently commented out in upstream
    #[allow(dead_code)]
    Wind,
}

fn init_demand_list(timestep_count: usize) -> Vec<AtomicF64> {
    (0..timestep_count)
        .map(|_| Default::default())
        .collect::<Vec<_>>()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    pub fn simulation_time() -> SimulationTimeIterator {
        SimulationTime::new(0.0, 8.0, 1.0).iter()
    }

    #[fixture]
    pub fn energy_supply<'a>(simulation_time: SimulationTimeIterator) -> EnergySupply {
        let mut energy_supply = EnergySupply::new(
            FuelType::MainsGas,
            simulation_time.total_steps(),
            None,
            None,
            None,
        );
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
        Arc<RwLock<EnergySupply>>,
    ) {
        let shared_supply = Arc::new(RwLock::new(energy_supply));
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
            energy_supply: Arc::new(RwLock::new(energy_supply)),
            end_user_name: "shower".to_string(),
        }
    }

    #[fixture]
    pub fn energy_supply_connection_2<'a>(energy_supply: EnergySupply) -> EnergySupplyConnection {
        EnergySupplyConnection {
            energy_supply: Arc::new(RwLock::new(energy_supply)),
            end_user_name: "bath".to_string(),
        }
    }

    const EXPECTED_TOTAL_DEMANDS: [f64; 8] =
        [50.0, 120.0, 190.0, 260.0, 330.0, 400.0, 470.0, 540.0];

    #[rstest]
    pub fn test_results_total(
        energy_supply: EnergySupply,
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
            Arc<RwLock<EnergySupply>>,
        ),
        simulation_time: SimulationTimeIterator,
    ) {
        let (energy_connection_1, energy_connection_2, energy_supply) = energy_supply_connections;
        for simtime in simulation_time {
            let _ = energy_connection_1
                .demand_energy((simtime.index as f64 + 1.0) * 50.0, simtime.index);
            let _ = energy_connection_2.demand_energy(simtime.index as f64 * 20.0, simtime.index);
            assert_eq!(
                energy_supply.read().results_by_end_user()["shower"][simtime.index],
                EXPECTED_TOTAL_DEMANDS_BY_END_USER[0][simtime.index]
            );
            assert_eq!(
                energy_supply.read().results_by_end_user()["bath"][simtime.index],
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
            Arc<RwLock<EnergySupply>>,
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

            let energy_supply = energy_supply.read();
            energy_supply
                .calc_energy_import_export_betafactor(t_it)
                .unwrap();

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

    #[rstest]
    pub fn test_energy_supply_without_export(simulation_time: SimulationTimeIterator) {
        let energy_supply = EnergySupply::new(
            FuelType::MainsGas,
            simulation_time.total_steps(),
            None,
            None,
            Some(false),
        );
        let shared_supply = Arc::new(RwLock::new(energy_supply));
        let energy_connection_1 =
            EnergySupply::connection(shared_supply.clone(), "shower").unwrap();
        let energy_connection_2 = EnergySupply::connection(shared_supply.clone(), "bath").unwrap();
        for t_it in simulation_time {
            let t_idx = t_it.index;
            energy_connection_1
                .demand_energy(((t_idx + 1) * 50) as f64, t_idx)
                .unwrap();
            energy_connection_2
                .demand_energy((t_idx * 20) as f64, t_idx)
                .unwrap();
            assert_eq!(
                shared_supply.read().get_energy_export()[t_idx],
                0.,
                "incorrect energy export returned"
            );
        }
    }
}
