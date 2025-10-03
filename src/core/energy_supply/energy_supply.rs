use crate::compare_floats::min_of_2;
use crate::core::energy_supply::elec_battery::ElectricBattery;
use crate::core::energy_supply::tariff_data::TariffData;
use crate::core::heating_systems::storage_tank::SurplusDiverting;
use crate::errors::NotImplementedError;
use crate::input::{EnergySupplyPriorityEntry, EnergySupplyTariff, FuelType};
use crate::simulation_time::SimulationTimeIteration;
use anyhow::{anyhow, bail};
use atomic_float::AtomicF64;
use derivative::Derivative;
use indexmap::{indexmap, IndexMap};
use parking_lot::RwLock;
use smartstring::alias::String;
use std::io::Read;
use std::sync::atomic::Ordering;
use std::sync::Arc;

pub(crate) const UNMET_DEMAND_SUPPLY_NAME: &str = "_unmet_demand";
pub(crate) const ENERGY_FROM_ENVIRONMENT_SUPPLY_NAME: &str = "_energy_from_environment";
/// An object to represent the connection of a system that consumes energy to the energy supply
///
/// This object encapsulates the name of the connection, meaning that the
/// system consuming the energy does not have to specify these on every call,
/// and helping to enforce that each connection to a single supply has a unique
/// name.
#[derive(Clone, Debug)]
pub(crate) struct EnergySupplyConnection {
    energy_supply: Arc<RwLock<EnergySupply>>,
    pub(crate) end_user_name: String,
}

impl EnergySupplyConnection {
    pub(crate) fn new(energy_supply: Arc<RwLock<EnergySupply>>, end_user_name: String) -> Self {
        Self {
            energy_supply,
            end_user_name,
        }
    }

    /// Forwards the amount of energy out (in kWh) to the relevant EnergySupply object
    pub(crate) fn energy_out(
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
    pub(crate) fn demand_energy(
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

    pub(crate) fn supply_energy(
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
}

pub(crate) struct EnergySupplyTariffInput {
    tariff: EnergySupplyTariff,
    tariff_data: Box<dyn Read>,
    threshold_charges: Vec<f64>,
    threshold_prices: Vec<f64>,
}

impl EnergySupplyTariffInput {
    pub(crate) fn new(
        tariff: EnergySupplyTariff,
        tariff_data: Box<dyn Read>,
        threshold_charges: Vec<f64>,
        threshold_prices: Vec<f64>,
    ) -> Self {
        Self {
            tariff,
            tariff_data,
            threshold_charges,
            threshold_prices,
        }
    }
}

#[derive(Debug)]
struct EnergySupplyTariffInfo {
    tariff: EnergySupplyTariff,
    tariff_data: TariffData,
    threshold_charges: Vec<f64>,
    threshold_prices: Vec<f64>,
}

impl TryFrom<EnergySupplyTariffInput> for EnergySupplyTariffInfo {
    type Error = anyhow::Error;

    fn try_from(input: EnergySupplyTariffInput) -> Result<Self, Self::Error> {
        Ok(Self {
            tariff: input.tariff,
            tariff_data: TariffData::new(input.tariff_data)?,
            threshold_charges: input.threshold_charges,
            threshold_prices: input.threshold_prices,
        })
    }
}

#[derive(Derivative)]
#[derivative(Debug)]
pub(crate) struct EnergySupply {
    fuel_type: FuelType,
    tariff_info: Option<EnergySupplyTariffInfo>,
    simulation_timesteps: usize,
    electric_battery: Option<ElectricBattery>,
    #[derivative(Debug = "ignore")]
    diverter: Option<Arc<RwLock<dyn SurplusDiverting>>>,
    priority: Option<Vec<EnergySupplyPriorityEntry>>,
    is_export_capable: bool,
    demand_total: Vec<AtomicF64>,
    demand_by_end_user: IndexMap<String, Vec<AtomicF64>>,
    energy_out_by_end_user: IndexMap<String, Vec<AtomicF64>>,
    beta_factor: Vec<AtomicF64>,
    supply_surplus: Vec<AtomicF64>,
    demand_not_met: Vec<AtomicF64>,
    energy_into_battery_from_generation: Vec<AtomicF64>,
    energy_out_of_battery: Vec<AtomicF64>,
    energy_into_battery_from_grid: Vec<AtomicF64>,
    battery_state_of_charge: Vec<AtomicF64>,
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
    pub(crate) fn new(
        fuel_type: FuelType,
        simulation_timesteps: usize,
        tariff_input: Option<EnergySupplyTariffInput>,
        electric_battery: Option<ElectricBattery>,
        priority: Option<Vec<EnergySupplyPriorityEntry>>,
        is_export_capable: Option<bool>,
    ) -> anyhow::Result<Self> {
        let tariff_info = if electric_battery
            .as_ref()
            .is_some_and(|battery| battery.is_grid_charging_possible())
        {
            if let Some(tariff_input) = tariff_input {
                Some(tariff_input.try_into()?)
            } else {
                bail!("A battery that can be charged from the grid is present but no tariff data source was provided.");
            }
        } else {
            None
        };

        Ok(Self {
            fuel_type,
            simulation_timesteps,
            tariff_info,
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
            energy_into_battery_from_generation: init_demand_list(simulation_timesteps),
            energy_out_of_battery: init_demand_list(simulation_timesteps),
            energy_into_battery_from_grid: init_demand_list(simulation_timesteps),
            battery_state_of_charge: init_demand_list(simulation_timesteps),
            energy_diverted: init_demand_list(simulation_timesteps),
            energy_generated_consumed: init_demand_list(simulation_timesteps),
        })
    }

    pub(crate) fn fuel_type(&self) -> FuelType {
        self.fuel_type
    }

    pub(crate) fn has_battery(&self) -> bool {
        self.electric_battery.is_some()
    }

    pub(crate) fn get_battery_max_capacity(&self) -> Option<f64> {
        self.electric_battery
            .as_ref()
            .map(|battery| battery.get_max_capacity())
    }

    pub(crate) fn get_battery_charge_efficiency(
        &self,
        simtime: SimulationTimeIteration,
    ) -> Option<f64> {
        self.electric_battery
            .as_ref()
            .map(|battery| battery.get_charge_efficiency(simtime))
    }

    pub(crate) fn get_battery_discharge_efficiency(
        &self,
        simtime: SimulationTimeIteration,
    ) -> Option<f64> {
        self.electric_battery
            .as_ref()
            .map(|battery| battery.get_discharge_efficiency(simtime))
    }

    pub(crate) fn get_battery_max_discharge(&self, charge: f64) -> Option<f64> {
        self.electric_battery
            .as_ref()
            .map(|battery| battery.calculate_max_discharge(charge))
    }

    pub(crate) fn get_battery_available_charge(&self) -> Option<f64> {
        self.electric_battery
            .as_ref()
            .map(|battery| battery.get_state_of_charge() * battery.get_max_capacity())
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

        Ok(EnergySupplyConnection::new(
            energy_supply.clone(),
            end_user_name.into(),
        ))
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
        diverter: Arc<RwLock<dyn SurplusDiverting>>,
    ) -> Result<(), &'static str> {
        if self.diverter.is_some() {
            return Err("diverter was already connected");
        }

        self.diverter = Some(diverter);

        Ok(())
    }

    /// This method is used in place of calling .connection() in the Python codebase in order to register an end user name
    #[cfg(test)]
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
        self.demand_energy(end_user_name, -amount_produced, timestep_index)
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

    /// Return the demand from each end user on this energy source for this timestep.
    /// Returns dictionary of floats, where dictionary keys are names of end users.
    pub(crate) fn results_by_end_user_single_step(&self, t_idx: usize) -> IndexMap<String, f64> {
        self.demand_by_end_user
            .keys()
            .map(|user_name| {
                (
                    user_name.clone(),
                    if self.energy_out_by_end_user.contains_key(user_name) {
                        self.demand_by_end_user[user_name][t_idx].load(Ordering::SeqCst)
                            + self.energy_out_by_end_user[user_name][t_idx].load(Ordering::SeqCst)
                    } else {
                        self.demand_by_end_user[user_name][t_idx].load(Ordering::SeqCst)
                    },
                )
            })
            .collect()
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
    pub fn get_energy_to_from_battery(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        (
            Self::vec_of_floats_from_atomics(&self.energy_into_battery_from_generation),
            Self::vec_of_floats_from_atomics(&self.energy_out_of_battery),
            Self::vec_of_floats_from_atomics(&self.energy_into_battery_from_grid),
            Self::vec_of_floats_from_atomics(&self.battery_state_of_charge),
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

    /// Check whether the Electric Battery is in a state where we allow charging from the grid
    ///       This function is called at two different stages in the calculation:
    ///       1. When considering discharging from the battery (electric demand from house)
    ///       2. When considering charging from the grid
    ///
    /// return parameters are:
    ///       charging_condition      -- charging condition combining the price and charge thresholds criteria
    ///       threshold_charge        -- threshold charge for current timestep
    ///       can_charge_if_not_full  -- just the price threshold criteria for charging
    pub(crate) fn is_charging_from_grid(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(bool, f64, bool)> {
        // TODO (from Python): Additional logic for grid charging decision
        //      Negative prices - Priority over PV? That would mean calling the function twice
        //                        Once before PV and again after but flagging if charging was
        //                        done in the first call.
        //      PV generation   - Currently set as priority for battery charging
        //      Seasonal threshold - Improve approach for charge threshold to cut grid charging when more PV available
        //
        let t_idx = simtime.index;
        let month = simtime.current_month().ok_or_else(|| {
            anyhow!("Month could not be resolved for current simulation timestep.")
        })? as usize;
        let EnergySupplyTariffInfo {
            tariff,
            tariff_data,
            threshold_charges,
            threshold_prices,
        } = self
            .tariff_info
            .as_ref()
            .ok_or_else(|| anyhow!("Tariff info not set when expected."))?;
        let threshold_charge = *threshold_charges
            .get(month)
            .ok_or_else(|| anyhow!("Threshold charge not set for month {month}."))?;
        let threshold_price = *threshold_prices
            .get(month)
            .ok_or_else(|| anyhow!("Threshold price not set for month {month}."))?;
        // For tariff selected look up price etc and decide whether to charge
        let elec_price = tariff_data.price(tariff, t_idx)?;
        let (current_charge, charge_discharge_efficiency) = {
            let battery = self
                .electric_battery
                .as_ref()
                .expect("Electric battery expected to be set if tariff data is.");
            (
                battery.get_state_of_charge(),
                battery.get_charge_discharge_efficiency(),
            )
        };

        Ok(
            match (
                elec_price / charge_discharge_efficiency < threshold_price,
                current_charge < threshold_charge,
            ) {
                (false, _) => (false, threshold_charge, false),
                (true, charging_condition) => (charging_condition, threshold_charge, true),
            },
        )
    }

    pub(crate) fn calc_energy_import_from_grid_to_battery(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Option<f64>> {
        if let Some(electric_battery) = &self.electric_battery {
            let t_idx = simtime.index;
            if electric_battery.is_grid_charging_possible() {
                // Current conditions of the battery
                let current_charge = electric_battery.get_state_of_charge();
                let max_capacity = electric_battery.get_max_capacity();

                let (charging_condition, threshold_charge, _) = self.is_charging_from_grid(simtime).expect("Expected to be able to determine whether charging from grid if grid charging is possible on battery.");
                let energy_accepted = if charging_condition {
                    // Create max elec_demand from grid to complete battery charging if battery conditions allow
                    let elec_demand = -max_capacity * (threshold_charge - current_charge)
                        / electric_battery.get_charge_efficiency(simtime);
                    // Attempt charging battery and retrieving energy_accepted
                    let energy_accepted =
                        -electric_battery.charge_discharge_battery(elec_demand, false, simtime);
                    self.energy_into_battery_from_grid[t_idx]
                        .store(energy_accepted, Ordering::SeqCst);

                    // Informing EnergyImport of imported electricity
                    self.demand_not_met[t_idx].fetch_add(elec_demand, Ordering::SeqCst);
                    energy_accepted
                } else {
                    0.0
                };

                // this function is called at the end of the timestep to reset time charging etc:
                electric_battery.timestep_end();
                self.battery_state_of_charge[t_idx]
                    .store(electric_battery.get_state_of_charge(), Ordering::SeqCst);
                Ok(Some(energy_accepted))
            } else {
                electric_battery.timestep_end();
                self.battery_state_of_charge[t_idx]
                    .store(electric_battery.get_state_of_charge(), Ordering::SeqCst);
                Ok(None)
            }
        } else {
            Ok(None)
        }
    }

    /// Calculate how much of that supply can be offset against demand.
    /// And then calculate what demand and supply is left after offsetting, which are the amount exported imported
    pub fn calc_energy_import_export_betafactor(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<()> {
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

        // See if there is a net supply/demand for the timestep
        match &self.priority {
            None => {
                if let Some(ref battery) = &self.electric_battery {
                    // See if the battery can deal with excess supply/demand for this timestep
                    // supply_surplus is -ve by convention and demand_not_met is +ve
                    let (charging_condition, _, can_charge_if_not_full) = if battery
                        .is_grid_charging_possible()
                    {
                        self.is_charging_from_grid(simtime).expect("Expected to be able to determine whether charging from grid if grid charging is possible on battery.")
                    } else {
                        (false, Default::default(), false)
                    };
                    if supply_surplus < 0. {
                        let energy_out_of_battery = battery.charge_discharge_battery(
                            supply_surplus,
                            charging_condition,
                            simtime,
                        );
                        supply_surplus -= energy_out_of_battery;
                        self.energy_into_battery_from_generation[timestep_idx]
                            .store(-energy_out_of_battery, Ordering::SeqCst);
                    }
                    if demand_not_met > 0. {
                        // Calling is_charging_from_grid threshold level to avoid
                        // discharging from the electric battery and
                        // triggering lots of small grid recharge events
                        // when the level of charge is close to the threshold
                        if !can_charge_if_not_full {
                            let energy_out_of_battery =
                                battery.charge_discharge_battery(demand_not_met, false, simtime);
                            demand_not_met -= energy_out_of_battery;
                            self.energy_out_of_battery[timestep_idx]
                                .store(-energy_out_of_battery, Ordering::SeqCst);
                        }
                    }
                }

                if let Some(ref diverter) = &self.diverter {
                    self.energy_diverted.get(timestep_idx).unwrap().store(
                        diverter.read().divert_surplus(supply_surplus, simtime)?,
                        Ordering::SeqCst,
                    );
                    supply_surplus += self.energy_diverted[timestep_idx].load(Ordering::SeqCst);
                }
            }
            Some(priority) => {
                for item in priority {
                    if matches!(item, EnergySupplyPriorityEntry::ElectricBattery)
                        && self.electric_battery.is_some()
                    {
                        let electric_battery = self.electric_battery.as_ref().unwrap();
                        let (charging_condition, _, can_charge_if_not_full) = if electric_battery
                            .is_grid_charging_possible()
                        {
                            self.is_charging_from_grid(simtime).expect("Expected to be able to determine whether charging from grid if grid charging is possible on battery.")
                        } else {
                            (false, Default::default(), false)
                        };
                        let energy_out_of_battery = electric_battery.charge_discharge_battery(
                            supply_surplus,
                            charging_condition,
                            simtime,
                        );
                        supply_surplus -= energy_out_of_battery;
                        self.energy_into_battery_from_generation[simtime.index]
                            .store(-energy_out_of_battery, Ordering::SeqCst);
                        let energy_out_of_battery = electric_battery.charge_discharge_battery(
                            demand_not_met,
                            can_charge_if_not_full,
                            simtime,
                        );
                        demand_not_met -= energy_out_of_battery;
                        self.energy_out_of_battery[simtime.index]
                            .store(-energy_out_of_battery, Ordering::SeqCst);
                    } else if matches!(item, EnergySupplyPriorityEntry::Diverter)
                        && self.diverter.is_some()
                    {
                        let diverter = self.diverter.as_ref().unwrap();
                        self.energy_diverted[simtime.index].store(
                            diverter.read().divert_surplus(supply_surplus, simtime)?,
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

pub(crate) struct EnergySupplyBuilder {
    energy_supply: EnergySupply,
}

impl EnergySupplyBuilder {
    pub(crate) fn new(fuel_type: FuelType, simulation_timesteps: usize) -> Self {
        Self {
            energy_supply: EnergySupply::new(
                fuel_type,
                simulation_timesteps,
                None,
                None,
                None,
                None,
            )
            .unwrap(),
        }
    }

    pub(crate) fn with_export_capable(mut self, is_export_capable: bool) -> Self {
        self.energy_supply.is_export_capable = is_export_capable;
        self
    }

    pub(crate) fn with_tariff_input(
        mut self,
        tariff_input: EnergySupplyTariffInput,
    ) -> anyhow::Result<Self> {
        self.energy_supply.tariff_info = Some(tariff_input.try_into()?);
        Ok(self)
    }

    pub(crate) fn with_electric_battery(mut self, electric_battery: ElectricBattery) -> Self {
        self.energy_supply.electric_battery = Some(electric_battery);
        self
    }

    pub(crate) fn with_priority(mut self, priority: Vec<EnergySupplyPriorityEntry>) -> Self {
        self.energy_supply.priority = Some(priority);
        self
    }

    pub(crate) fn build(self) -> EnergySupply {
        self.energy_supply
    }

    // write other builder methods
}

fn init_demand_list(timestep_count: usize) -> Vec<AtomicF64> {
    (0..timestep_count)
        .map(|_| Default::default())
        .collect::<Vec<_>>()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::external_conditions::{DaylightSavingsConfig, ExternalConditions};
    use crate::input::BatteryLocation;
    use crate::simulation_time::SimulationTime;
    use itertools::Itertools;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use serde_json::json;

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0.0, 8.0, 1.0)
    }

    #[fixture]
    pub fn energy_supply<'a>(simulation_time: SimulationTime) -> EnergySupply {
        let mut energy_supply =
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time.total_steps()).build();
        energy_supply.register_end_user_name("shower".into());
        energy_supply.register_end_user_name("bath".into());

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
            end_user_name: "shower".into(),
        };
        let energy_connection_2 = EnergySupplyConnection {
            energy_supply: shared_supply.clone(),
            end_user_name: "bath".into(),
        };
        (energy_connection_1, energy_connection_2, shared_supply)
    }

    #[fixture]
    pub fn energy_supply_connection_1<'a>(energy_supply: EnergySupply) -> EnergySupplyConnection {
        EnergySupplyConnection {
            energy_supply: Arc::new(RwLock::new(energy_supply)),
            end_user_name: "shower".into(),
        }
    }

    #[fixture]
    pub fn energy_supply_connection_2<'a>(energy_supply: EnergySupply) -> EnergySupplyConnection {
        EnergySupplyConnection {
            energy_supply: Arc::new(RwLock::new(energy_supply)),
            end_user_name: "bath".into(),
        }
    }

    #[rstest]
    pub fn test_init_demand_list(simulation_time: SimulationTime) {
        assert_eq!(
            init_demand_list(simulation_time.total_steps()),
            [0.; 8].into_iter().map(AtomicF64::new).collect::<Vec<_>>()
        );
    }

    #[rstest]
    fn test_energy_out(energy_supply: EnergySupply, simulation_time: SimulationTime) {
        // Check with existing end user name
        let amount_demand = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0];

        for (t_idx, _) in simulation_time.iter().enumerate() {
            energy_supply
                .energy_out("shower", amount_demand[t_idx], t_idx)
                .unwrap();
            assert_eq!(
                energy_supply.energy_out_by_end_user["shower"][t_idx].load(Ordering::SeqCst),
                amount_demand[t_idx]
            );
        }
        // Check an error is raised with new end user name
        assert!(energy_supply.energy_out("electricshower", 10., 0).is_err());
    }

    #[fixture]
    fn pv_diverter() -> Arc<RwLock<dyn SurplusDiverting>> {
        struct NullDiverter;

        impl SurplusDiverting for NullDiverter {
            fn divert_surplus(
                &self,
                _surplus: f64,
                _simtime: SimulationTimeIteration,
            ) -> anyhow::Result<f64> {
                Ok(0.)
            }
        }

        Arc::new(RwLock::new(NullDiverter))
    }

    #[rstest]
    fn test_connect_diverter(
        mut energy_supply: EnergySupply,
        pv_diverter: Arc<RwLock<dyn SurplusDiverting>>,
    ) {
        assert!(energy_supply.diverter.is_none());
        energy_supply.connect_diverter(pv_diverter.clone()).unwrap();
        assert!(energy_supply.diverter.is_some());
        assert!(energy_supply.connect_diverter(pv_diverter.clone()).is_err());
    }

    #[rstest]
    fn test_demand_energy(energy_supply: EnergySupply, simulation_time: SimulationTime) {
        let amount_demanded = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0];
        for (t_idx, _) in simulation_time.iter().enumerate() {
            energy_supply
                .demand_energy("shower", amount_demanded[t_idx], t_idx)
                .unwrap();
            assert_eq!(
                energy_supply.demand_total[t_idx].load(Ordering::SeqCst),
                amount_demanded[t_idx]
            );
            assert_eq!(
                energy_supply.demand_by_end_user["shower"][t_idx].load(Ordering::SeqCst),
                amount_demanded[t_idx]
            );
        }
    }

    #[rstest]
    fn test_supply_energy(energy_supply: EnergySupply, simulation_time: SimulationTime) {
        let amount_produced = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0];
        for (t_idx, _) in simulation_time.iter().enumerate() {
            energy_supply
                .supply_energy("shower", amount_produced[t_idx], t_idx)
                .unwrap();
            assert_eq!(
                energy_supply.demand_total[t_idx].load(Ordering::SeqCst),
                [-10.0, -20.0, -30.0, -40.0, -50.0, -60.0, -70.0, -80.0][t_idx]
            );
            assert_eq!(
                energy_supply.demand_by_end_user["shower"][t_idx].load(Ordering::SeqCst),
                [-10.0, -20.0, -30.0, -40.0, -50.0, -60.0, -70.0, -80.0][t_idx]
            );
        }
    }

    const EXPECTED_TOTAL_DEMANDS: [f64; 8] =
        [50.0, 120.0, 190.0, 260.0, 330.0, 400.0, 470.0, 540.0];

    #[rstest]
    pub fn test_results_total(energy_supply: EnergySupply, simulation_time: SimulationTime) {
        for simtime in simulation_time.iter() {
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
        simulation_time: SimulationTime,
    ) {
        let (energy_connection_1, energy_connection_2, energy_supply) = energy_supply_connections;
        for simtime in simulation_time.iter() {
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
        simulation_time: SimulationTime,
    ) {
        let (energy_connection_1, energy_connection_2, energy_supply) = energy_supply_connections;
        let energy_connection_3 = EnergySupply::connection(energy_supply.clone(), "PV").unwrap();
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
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

    #[fixture]
    fn external_conditions(simulation_time: SimulationTime) -> ExternalConditions {
        ExternalConditions::new(
            &simulation_time.iter(),
            vec![0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0],
            vec![3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1],
            vec![0., 20., 40., 60., 0., 20., 40., 60.],
            vec![11., 25., 42., 52., 60., 44., 28., 15.],
            vec![11., 25., 42., 52., 60., 44., 28., 15.],
            vec![0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2],
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
            serde_json::from_value(json!([
                // upstream Python gives old 'start' fields, but we need 'start360' here
                {"number": 1, "start360": 0, "end360": 45},
                {"number": 2, "start360": 45, "end360": 90,
                 "shading": [
                     {"type": "overhang", "height": 2.2, "distance": 6}
                     ]
                 },
                {"number": 3, "start360": 90, "end360": 135},
                {"number": 4, "start360": 135, "end360": 180,
                 "shading": [
                     {"type": "obstacle", "height": 40, "distance": 4},
                     {"type": "overhang", "height": 3, "distance": 7}
                     ]
                 },
                {"number": 5, "start360": 180, "end360": 225,
                 "shading": [
                     {"type": "obstacle", "height": 3, "distance": 8},
                     ]
                 },
                {"number": 6, "start360": 225, "end360": 270},
                {"number": 7, "start360": 270, "end360": 315},
                {"number": 8, "start360": 315, "end360": 360}
            ]))
            .unwrap(),
        )
    }

    #[rstest]
    fn test_calc_energy_import_export_betafactor(
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        let amount_demanded = [50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0];
        let amount_produced = [50.0, 90.0, 130.0, 210.0, 2300.0, 290.0, 300.0, 350.0];

        let elec_battery = ElectricBattery::new(
            2.,
            0.8,
            3.,
            0.001,
            1.5,
            1.5,
            BatteryLocation::Outside,
            false,
            simulation_time.step,
            Arc::new(external_conditions.clone()),
        );

        let builder =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps());
        let energy_supply = builder.with_electric_battery(elec_battery).build();

        let energy_supply = Arc::new(RwLock::new(energy_supply));

        // test with elec battery
        let _shower_connection = EnergySupply::connection(energy_supply.clone(), "shower").unwrap();
        let _bath_connection = EnergySupply::connection(energy_supply.clone(), "bath").unwrap();

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let energy_supply = energy_supply.read();
            energy_supply
                .demand_energy("shower", amount_demanded[t_idx], t_idx)
                .unwrap();
            energy_supply
                .supply_energy("bath", amount_produced[t_idx], t_idx)
                .unwrap();
            energy_supply
                .calc_energy_import_export_betafactor(t_it)
                .unwrap();
        }

        {
            let energy_supply = energy_supply.read();

            assert_eq!(
                energy_supply
                    .demand_total
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![0., 10., 20., -10., -2050., 10., 50., 50.]
            );

            assert_eq!(
                energy_supply
                    .demand_not_met
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    15.138528000000004,
                    34.37872423953049,
                    52.991809292243516,
                    63.07009986960006,
                    -2.842170943040401e-14,
                    99.5880926222444,
                    124.3891811736907,
                    140.57522007095577
                ]
            );

            assert_eq!(
                energy_supply
                    .supply_surplus
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    -14.582949016875158,
                    -24.59889302603037,
                    -32.991809292243516,
                    -73.07009986960004,
                    -2050.,
                    -89.5880926222444,
                    -74.38918117369072,
                    -90.5752200709558
                ]
            );

            assert_eq!(
                energy_supply
                    .energy_generated_consumed
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    33.739999999999995,
                    65.40110697396963,
                    97.00819070775648,
                    136.92990013039994,
                    250.00000000000003,
                    200.4119073777556,
                    225.6108188263093,
                    259.42477992904423
                ]
            );

            assert_eq!(
                energy_supply
                    .energy_into_battery_from_generation
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![1.6770509831248424, -0., -0., -0., -0., -0., -0., -0.]
            );

            assert_eq!(
                energy_supply
                    .energy_out_of_battery
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![-1.121472, -0.2201687864998738, -0., -0., 0., -0., -0., -0.]
            );

            assert_eq!(
                energy_supply
                    .energy_diverted
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![0., 0., 0., 0., 0., 0., 0., 0.]
            );

            assert_eq!(
                energy_supply
                    .beta_factor
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    0.6748,
                    0.7266789663774403,
                    0.7462168515981268,
                    0.652047143478095,
                    0.10869565217391305,
                    0.6910755426819158,
                    0.7520360627543643,
                    0.7412136569401263
                ]
            );
        }

        // Test with PV diverter
        struct MockDiverter;

        impl SurplusDiverting for MockDiverter {
            fn divert_surplus(
                &self,
                _supply_surplus: f64,
                _simulation_time_iteration: SimulationTimeIteration,
            ) -> anyhow::Result<f64> {
                Ok(10.)
            }
        }

        let diverter = Arc::new(RwLock::new(MockDiverter));
        energy_supply.write().connect_diverter(diverter).unwrap();

        for (t_idx, simtime) in simulation_time.iter().enumerate() {
            let energy_supply = energy_supply.read();

            energy_supply
                .demand_energy("shower", amount_demanded[t_idx], t_idx)
                .unwrap();
            energy_supply
                .supply_energy("bath", amount_produced[t_idx], t_idx)
                .unwrap();
            energy_supply
                .calc_energy_import_export_betafactor(simtime)
                .unwrap();
        }

        {
            let energy_supply = energy_supply.read();

            assert_eq!(
                energy_supply
                    .demand_total
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![0., 20., 40., -20., -4100., 20., 100., 100.]
            );

            assert_eq!(
                energy_supply
                    .demand_not_met
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    47.65852800000002,
                    103.57651029159123,
                    158.97542787673055,
                    189.21029960880017,
                    -8.526512829121202e-14,
                    298.7642778667332,
                    373.1675435210721,
                    421.7256602128673
                ]
            );

            assert_eq!(
                energy_supply
                    .supply_surplus
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    -37.10294901687516,
                    -63.79667907809112,
                    -88.97542787673055,
                    -209.21029960880014,
                    -6140.0,
                    -258.7642778667332,
                    -213.16754352107216,
                    -261.7256602128674
                ]
            );

            assert_eq!(
                energy_supply
                    .energy_generated_consumed
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    101.21999999999998,
                    196.2033209219089,
                    291.0245721232694,
                    410.78970039119986,
                    750.0000000000001,
                    601.2357221332668,
                    676.832456478928,
                    778.2743397871327
                ]
            );

            assert_eq!(
                energy_supply
                    .energy_into_battery_from_generation
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![-0., -0., -0., -0., -0., -0., -0., -0.]
            );

            assert_eq!(
                energy_supply
                    .energy_out_of_battery
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![-0., -0., -0., -0., 0., -0., -0., -0.]
            );

            assert_eq!(
                energy_supply
                    .energy_diverted
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![10., 10., 10., 10., 10., 10., 10., 10.]
            );

            assert_eq!(
                energy_supply
                    .beta_factor
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    0.6748,
                    0.7266789663774403,
                    0.7462168515981268,
                    0.652047143478095,
                    0.10869565217391305,
                    0.6910755426819158,
                    0.7520360627543643,
                    0.7412136569401263
                ]
            );
        }

        // important so the energy supply Arc only has one strong reference to the energy supply
        // below where Arc::into_inner is called, otherwise that call would fail
        drop(_shower_connection);
        drop(_bath_connection);

        // LOOK AWAY 👀👀👀👀👀
        // (the upstream Python shared the same electric battery across energy supplies in this test,
        // so its internal state is not isolated - therefore we need to cannibalise the previous energy
        // supply here for scraps (the electric battery) for use in the next set of assertions)
        let elec_battery = Arc::into_inner(energy_supply)
            .unwrap()
            .into_inner()
            .electric_battery
            .unwrap();

        // Set priority
        let priority = vec![
            EnergySupplyPriorityEntry::Diverter,
            EnergySupplyPriorityEntry::ElectricBattery,
        ];

        let mut builder =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps());
        builder = builder
            .with_electric_battery(elec_battery)
            .with_priority(priority);

        let energy_supply = Arc::new(RwLock::new(builder.build()));

        let _shower_connection = EnergySupply::connection(energy_supply.clone(), "shower").unwrap();
        let _bath_connection = EnergySupply::connection(energy_supply.clone(), "bath").unwrap();

        let diverter = Arc::new(RwLock::new(MockDiverter));
        energy_supply.write().connect_diverter(diverter).unwrap();

        for (t_idx, simtime) in simulation_time.iter().enumerate() {
            let energy_supply = energy_supply.read();

            energy_supply
                .demand_energy("shower", amount_demanded[t_idx], t_idx)
                .unwrap();
            energy_supply
                .supply_energy("bath", amount_produced[t_idx], t_idx)
                .unwrap();
            energy_supply
                .calc_energy_import_export_betafactor(simtime)
                .unwrap();
        }

        {
            let energy_supply = energy_supply.read();

            assert_eq!(
                energy_supply
                    .demand_total
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![0., 10., 20., -10., -2050., 10., 50., 50.]
            );

            assert_eq!(
                energy_supply
                    .demand_not_met
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    16.260000000000005,
                    34.59889302603037,
                    52.991809292243516,
                    63.07009986960006,
                    -2.842170943040401e-14,
                    99.5880926222444,
                    124.3891811736907,
                    140.57522007095577
                ]
            );

            assert_eq!(
                energy_supply
                    .supply_surplus
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    -6.260000000000002,
                    -14.598893026030371,
                    -22.991809292243516,
                    -63.07009986960004,
                    -2040.0,
                    -79.5880926222444,
                    -64.38918117369072,
                    -80.5752200709558
                ]
            );

            assert_eq!(
                energy_supply
                    .energy_generated_consumed
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    33.739999999999995,
                    65.40110697396963,
                    97.00819070775648,
                    136.92990013039994,
                    250.00000000000003,
                    200.4119073777556,
                    225.6108188263093,
                    259.42477992904423
                ]
            );

            assert_eq!(
                energy_supply
                    .energy_into_battery_from_generation
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![-0., -0., -0., -0., -0., -0., -0., -0.]
            );

            assert_eq!(
                energy_supply
                    .energy_out_of_battery
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![-0., -0., -0., -0., -0., -0., -0., -0.]
            );

            assert_eq!(
                energy_supply
                    .energy_diverted
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![10., 10., 10., 10., 10., 10., 10., 10.]
            );

            assert_eq!(
                energy_supply
                    .beta_factor
                    .iter()
                    .map(|x| x.load(Ordering::SeqCst))
                    .collect_vec(),
                vec![
                    0.6748,
                    0.7266789663774403,
                    0.7462168515981268,
                    0.652047143478095,
                    0.10869565217391305,
                    0.6910755426819158,
                    0.7520360627543643,
                    0.7412136569401263
                ]
            );
        }
    }

    #[rstest]
    pub fn test_energy_supply_without_export(simulation_time: SimulationTime) {
        let mut builder =
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time.total_steps());
        builder = builder.with_export_capable(false);
        let energy_supply = builder.build();
        let shared_supply = Arc::new(RwLock::new(energy_supply));
        let energy_connection_1 =
            EnergySupply::connection(shared_supply.clone(), "shower").unwrap();
        let energy_connection_2 = EnergySupply::connection(shared_supply.clone(), "bath").unwrap();
        for t_it in simulation_time.iter() {
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
