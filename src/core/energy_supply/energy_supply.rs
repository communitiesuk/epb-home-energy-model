use crate::compare_floats::min_of_2;
use crate::core::energy_supply::elec_battery::ElectricBattery;
use crate::core::energy_supply::tariff_data::TariffData;
use crate::core::heating_systems::storage_tank::SurplusDiverting;
use crate::errors::NotImplementedError;
use crate::input::{EnergySupplyTariff, FuelType, SecondarySupplyType};
use crate::simulation_time::SimulationTimeIteration;
use anyhow::{anyhow, bail};
use atomic_float::AtomicF64;
use derivative::Derivative;
use indexmap::{indexmap, IndexMap};
use parking_lot::RwLock;
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
    pub(crate) energy_supply: Arc<RwLock<EnergySupply>>,
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
    pub(crate) electric_battery: Option<ElectricBattery>,
    #[derivative(Debug = "ignore")]
    pub(crate) diverter: Option<Arc<RwLock<dyn SurplusDiverting>>>,
    priority: Option<Vec<SecondarySupplyType>>,
    is_export_capable: bool,
    pub(crate) demand_total: Vec<AtomicF64>,
    pub(crate) demand_by_end_user: IndexMap<String, Vec<AtomicF64>>,
    pub(crate) energy_out_by_end_user: IndexMap<String, Vec<AtomicF64>>,
    pub(crate) beta_factor: Vec<AtomicF64>,
    pub(crate) supply_surplus: Vec<AtomicF64>,
    pub(crate) demand_not_met: Vec<AtomicF64>,
    pub(crate) energy_into_battery_from_generation: Vec<AtomicF64>,
    pub(crate) energy_out_of_battery: Vec<AtomicF64>,
    energy_into_battery_from_grid: Vec<AtomicF64>,
    battery_state_of_charge: Vec<AtomicF64>,
    pub(crate) energy_diverted: Vec<AtomicF64>,
    pub(crate) energy_generated_consumed: Vec<AtomicF64>,
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
        priority: Option<Vec<SecondarySupplyType>>,
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
                    if matches!(item, SecondarySupplyType::ElectricBattery)
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
                    } else if matches!(item, SecondarySupplyType::Diverter)
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

    pub(crate) fn with_priority(mut self, priority: Vec<SecondarySupplyType>) -> Self {
        self.energy_supply.priority = Some(priority);
        self
    }

    pub(crate) fn build(self) -> EnergySupply {
        self.energy_supply
    }

    // write other builder methods
}

pub(crate) fn init_demand_list(timestep_count: usize) -> Vec<AtomicF64> {
    (0..timestep_count)
        .map(|_| Default::default())
        .collect::<Vec<_>>()
}
