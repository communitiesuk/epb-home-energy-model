use crate::compare_floats::{max_of_2, min_of_2};
use crate::core::energy_supply::elec_battery::ElectricBattery;
use crate::core::energy_supply::tariff_data::TariffData;
use crate::core::heating_systems::storage_tank::SurplusDiverting;
use crate::errors::NotImplementedError;
use crate::input::{EnergySupplyTariff, FuelType};
use crate::simulation_time::SimulationTimeIteration;
use anyhow::{anyhow, bail};
use approx::relative_eq;
use atomic_float::AtomicF64;
use educe::Educe;
use fsum::FSum;
use indexmap::{indexmap, IndexMap};
use itertools::Itertools;
use parking_lot::RwLock;
use smartstring::alias::String;
use std::collections::HashSet;
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
    pub(crate) fn _energy_out(
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

// TODO 1.0.0a9 migration - add new EnergySupply threshold fields to this struct?
#[derive(Debug)]
pub struct EnergySupplyTariffInfo {
    pub(crate) tariff: EnergySupplyTariff,
    pub(crate) threshold_charges: Option<Vec<f64>>,
    pub(crate) threshold_prices: Option<Vec<f64>>,
}

#[derive(Educe)]
#[educe(Debug)]
pub struct EnergySupply {
    fuel_type: FuelType,
    tariff_info: Option<EnergySupplyTariffInfo>,
    tariff_data: Option<TariffData>,
    simulation_timesteps: usize,
    electric_batteries: IndexMap<String, Arc<ElectricBattery>>,
    #[educe(Debug(ignore))]
    diverters: IndexMap<String, Arc<RwLock<dyn SurplusDiverting>>>,
    priority: Option<Vec<String>>,
    is_export_capable: bool,
    power_limit_export: Option<f64>,
    tariff_export: Option<EnergySupplyTariff>,
    threshold_charges_export: Option<[f64; 12]>,
    threshold_prices_export: Option<[f64; 12]>,
    demand_total: Vec<AtomicF64>,
    demand_by_end_user: IndexMap<String, Vec<AtomicF64>>,
    energy_out_by_end_user: IndexMap<String, Vec<AtomicF64>>,
    beta_factor: Vec<AtomicF64>,
    supply_surplus: Vec<AtomicF64>,
    demand_not_met: Vec<AtomicF64>,
    grid_to_consumption: Vec<AtomicF64>,
    energy_into_battery_from_generation: Vec<AtomicF64>,
    energy_battery_to_consumption: Vec<AtomicF64>,
    energy_into_battery_from_grid: Vec<AtomicF64>,
    energy_into_grid_from_battery: Vec<AtomicF64>,
    battery_state_of_charge: Vec<AtomicF64>,
    energy_diverted: Vec<AtomicF64>,
    energy_generated_consumed: Vec<AtomicF64>,
    generation_curtailed: Vec<AtomicF64>,
    power_limit_battery_import: Option<f64>,
}

impl EnergySupply {
    /// Arguments:
    /// * `fuel_type` - string denoting type of fuel
    /// * `simulation_timesteps` - the number of steps in the simulation time being used
    /// * `electric_battery` - reference to a map from name to an ElectricBattery object
    /// * `is_export_capable` - denotes that this Energy Supply can export its surplus supply
    /// * `power_limit_export` - maximum AC power (kW) exportable to the grid, e.g. a
    //                           Distribution Network Operator (DNO) export limit; applies to
    //                           the whole supply (generation surplus and battery discharge
    //                           combined); None means no limit
    /// * `tariff_export` - energy tariff for export
    /// * `threshold_charges_export` - level of battery charge below which battery prohibited from exporting to grid (0 - 1)
    /// * `threshold_prices_export` - grid price above which battery is permitted to export to grid (p/kWh)
    /// * `tariff_data` - tariff data containing electricity prices
    /// * `power_limit_battery_import` - the maximum power limit for charging batteries from the
    //                                   energy supply connection (not limited to grid — also applies
    //                                   to on-site generation), shared across all batteries (kW)
    pub(crate) fn new(
        fuel_type: FuelType,
        simulation_timesteps: usize,
        tariff_info: Option<EnergySupplyTariffInfo>,
        tariff_data: Option<TariffData>,
        electric_batteries: IndexMap<String, ElectricBattery>,
        priority: Option<Vec<String>>,
        is_export_capable: Option<bool>,
        power_limit_export: Option<f64>,
        tariff_export: Option<EnergySupplyTariff>,
        threshold_charges_export: Option<[f64; 12]>,
        threshold_prices_export: Option<[f64; 12]>,
        power_limit_battery_import: Option<f64>,
    ) -> anyhow::Result<Self> {
        if electric_batteries
            .iter()
            .any(|(_, battery)| battery.is_grid_charging_possible())
            && (tariff_data.is_none() || tariff_info.is_none())
        {
            bail!(
                "A battery that can be charged from the grid is present but tariff data is missing"
            )
        };

        Ok(Self {
            fuel_type,
            simulation_timesteps,
            tariff_info,
            tariff_data,
            electric_batteries: electric_batteries
                .into_iter()
                .map(|(k, v)| (k, v.into()))
                .collect(),
            diverters: Default::default(),
            priority,
            is_export_capable: is_export_capable.unwrap_or(true),
            power_limit_export,
            tariff_export,
            threshold_charges_export,
            threshold_prices_export,
            demand_total: init_demand_list(simulation_timesteps),
            demand_by_end_user: Default::default(),
            energy_out_by_end_user: Default::default(),
            beta_factor: init_demand_list(simulation_timesteps),
            supply_surplus: init_demand_list(simulation_timesteps),
            demand_not_met: init_demand_list(simulation_timesteps),
            grid_to_consumption: init_demand_list(simulation_timesteps),
            energy_into_battery_from_generation: init_demand_list(simulation_timesteps),
            energy_battery_to_consumption: init_demand_list(simulation_timesteps),
            energy_into_battery_from_grid: init_demand_list(simulation_timesteps),
            energy_into_grid_from_battery: init_demand_list(simulation_timesteps),
            battery_state_of_charge: init_demand_list(simulation_timesteps),
            energy_diverted: init_demand_list(simulation_timesteps),
            energy_generated_consumed: init_demand_list(simulation_timesteps),
            generation_curtailed: init_demand_list(simulation_timesteps),
            power_limit_battery_import,
        })
    }

    pub(crate) fn timestep_end(&self) -> anyhow::Result<()> {
        for battery in self.get_batteries()? {
            battery.timestep_end()
        }

        Ok(())
    }

    pub(crate) fn fuel_type(&self) -> FuelType {
        self.fuel_type
    }

    pub(crate) fn get_diverters(&self) -> anyhow::Result<Vec<Arc<RwLock<dyn SurplusDiverting>>>> {
        self.sort_by_priority(&self.diverters)
    }

    pub(crate) fn get_batteries(&self) -> anyhow::Result<Vec<Arc<ElectricBattery>>> {
        self.sort_by_priority(&self.electric_batteries)
    }

    /// Returns the values from items sorted in the order that the keys appear in the priority list
    pub(crate) fn sort_by_priority<T: Clone>(
        &self,
        items: &IndexMap<String, T>,
    ) -> anyhow::Result<Vec<T>> {
        if let Some(priority) = self.priority.as_ref() {
            let priority_set: HashSet<&String> = HashSet::from_iter(priority);
            let items_set: HashSet<&String> = HashSet::from_iter(items.keys());
            if !priority_set.is_superset(&items_set) {
                bail!("Energy supply items missing from priority list")
            };

            Ok(priority
                .iter()
                .filter_map(|k| items.get(k).cloned())
                .collect())
        } else {
            Ok(items.values().cloned().collect())
        }
    }

    pub(crate) fn has_battery(&self) -> anyhow::Result<bool> {
        Ok(!self.get_batteries()?.is_empty())
    }

    #[cfg(test)]
    pub(crate) fn get_battery_max_capacity(&self) -> anyhow::Result<Option<f64>> {
        let batteries = self.get_batteries()?;
        if batteries.is_empty() {
            return Ok(None);
        };

        Ok(Some(
            batteries
                .iter()
                .map(|battery| battery.get_max_capacity())
                .sum(),
        ))
    }

    #[cfg(test)] // TODO 1.0.0a9 migration - this is only used in tests now, are these tests useful?
    pub(crate) fn get_battery_charge_efficiency(
        &self,
        simtime: SimulationTimeIteration,
        battery: Option<ElectricBattery>,
    ) -> anyhow::Result<Option<f64>> {
        match battery {
            None => {
                let batteries = self.get_batteries()?;

                if batteries.is_empty() {
                    Ok(None)
                } else if batteries.len() == 1 {
                    Ok(Some(batteries[0].get_charge_efficiency(simtime)))
                } else {
                    bail!("Battery not specified for function 'get_battery_charge_efficiency'")
                }
            }
            Some(battery) => Ok(Some(battery.get_charge_efficiency(simtime))),
        }
    }

    #[cfg(test)] // TODO 1.0.0a9 migration - this is only used in tests now, are these tests useful?
    pub(crate) fn get_battery_discharge_efficiency(
        &self,
        simtime: SimulationTimeIteration,
        battery: Option<ElectricBattery>,
    ) -> anyhow::Result<Option<f64>> {
        match battery {
            None => {
                let batteries = self.get_batteries()?;

                if batteries.is_empty() {
                    Ok(None)
                } else if batteries.len() == 1 {
                    Ok(Some(batteries[0].get_discharge_efficiency(simtime)))
                } else {
                    bail!("Battery not specified for function 'get_battery_discharge_efficiency'")
                }
            }
            Some(battery) => Ok(Some(battery.get_discharge_efficiency(simtime))),
        }
    }

    #[cfg(test)] // TODO 1.0.0a9 migration - this is only used in tests now, are these tests useful?
    pub(crate) fn get_battery_max_discharge(
        &self,
        charge: f64,
        battery: Option<ElectricBattery>,
    ) -> anyhow::Result<Option<f64>> {
        match battery {
            None => {
                let batteries = self.get_batteries()?;

                if batteries.is_empty() {
                    Ok(None)
                } else if batteries.len() == 1 {
                    Ok(Some(batteries[0].calculate_max_discharge(charge)))
                } else {
                    bail!("Battery not specified for function 'get_battery_max_discharge'")
                }
            }
            Some(battery) => Ok(Some(battery.calculate_max_discharge(charge))),
        }
    }

    pub(crate) fn get_battery_available_charge(&self) -> anyhow::Result<Option<f64>> {
        let batteries = self.get_batteries()?;
        if batteries.is_empty() {
            return Ok(None);
        };

        Ok(Some(
            batteries
                .iter()
                .map(|battery| battery.get_state_of_charge() * battery.get_max_capacity())
                .sum(),
        ))
    }

    pub(crate) fn connection(
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

    #[allow(dead_code)]
    fn energy_out(
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
        name: Option<String>,
    ) -> anyhow::Result<()> {
        let name = name.unwrap_or("diverter".into());

        if self.diverters.keys().contains(&name) {
            bail!("diverter was already connected");
        }

        self.diverters.insert(name, diverter);

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
    pub fn results_by_end_user(&self) -> IndexMap<Arc<str>, Vec<f64>> {
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
                        end_user.to_string().into(),
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
                    user_name.to_string().into(),
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
        let supply_surplus = Self::vec_of_floats_from_atomics(&self.supply_surplus);
        let energy_into_grid_from_battery =
            Self::vec_of_floats_from_atomics(&self.energy_into_grid_from_battery);

        (0..self.simulation_timesteps)
            .map(|i| supply_surplus[i] + energy_into_grid_from_battery[i])
            .collect()
    }

    pub fn get_energy_export_from_generation(&self) -> Vec<f64> {
        Self::vec_of_floats_from_atomics(&self.supply_surplus)
    }

    /// Return the amount of generated energy consumed in the building for all timesteps
    pub fn get_energy_generated_consumed(&self) -> Vec<f64> {
        Self::vec_of_floats_from_atomics(&self.energy_generated_consumed)
    }

    pub(crate) fn get_grid_to_consumption(&self) -> Vec<f64> {
        Self::vec_of_floats_from_atomics(&self.grid_to_consumption)
    }

    #[allow(clippy::type_complexity)]
    /// Return the amount of generated energy sent to battery and drawn from battery
    pub fn get_battery_energy_flows(&self) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        (
            Self::vec_of_floats_from_atomics(&self.energy_into_battery_from_generation),
            Self::vec_of_floats_from_atomics(&self.energy_battery_to_consumption),
            Self::vec_of_floats_from_atomics(&self.energy_into_battery_from_grid),
            Self::vec_of_floats_from_atomics(&self.energy_into_grid_from_battery),
            Self::vec_of_floats_from_atomics(&self.battery_state_of_charge),
        )
    }

    /// Return the amount of generated energy diverted to minimise export
    pub fn get_energy_diverted(&self) -> Vec<f64> {
        Self::vec_of_floats_from_atomics(&self.energy_diverted)
    }

    /// Return the generated energy curtailed by the export power limit for all timesteps.
    //  Curtailed energy is generation that could be neither consumed, stored, diverted
    //  nor exported, because export was capped at the DNO export power limit. It is
    //  reported so the generation balance closes; it is zero when no limit is set.
    pub fn get_energy_generation_curtailed(&self) -> Vec<f64> {
        Self::vec_of_floats_from_atomics(&self.generation_curtailed)
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

    // Check whether the Electric Battery is in a state where we allow exporting to grid
    /// return parameters are:
    ///       exporting_condition     -- exporting condition combining the price and charge thresholds criteria
    ///       threshold_charge        -- threshold charge for current timestep
    ///       can_export_if_not_empty -- just the price threshold criteria for charging
    pub(crate) fn is_exporting_to_grid(
        &self,
        battery: &ElectricBattery,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(bool, Option<f64>, bool)> {
        let month = simtime.current_month().ok_or_else(|| {
            anyhow!("Month could not be resolved for current simulation timestep.")
        })? as usize;
        let threshold_charge_export = self
            .threshold_charges_export
            .and_then(|charges| charges.get(month).copied());
        let threshold_price_export = self
            .threshold_prices_export
            .and_then(|charges| charges.get(month).copied());

        // For tariff selected look up price, etc, and decide whether to charge
        let elec_price = if let Some(tariff_export) = self.tariff_export {
            Some(
                self.tariff_data
                    .as_ref()
                    .ok_or_else(|| anyhow!("Tariff data expected to be set on energy supply"))?
                    .price(&tariff_export, simtime)?,
            )
        } else {
            None
        };

        let current_charge = battery.get_state_of_charge();
        let charge_discharge_efficiency = battery.get_charge_discharge_efficiency();

        Ok(match (elec_price, threshold_price_export) {
            (Some(elec_price), Some(threshold_price_export))
                if elec_price * charge_discharge_efficiency > threshold_price_export =>
            {
                let is_over_threshold =
                    threshold_charge_export.is_some_and(|threshold| current_charge > threshold);

                (is_over_threshold, threshold_charge_export, true)
            }
            _ => (false, threshold_charge_export, false),
        })
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
        battery: &ElectricBattery,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(bool, Option<f64>, bool)> {
        // TODO (from Python): Additional logic for grid charging decision
        //      Negative prices - Priority over PV? That would mean calling the function twice
        //                        Once before PV and again after but flagging if charging was
        //                        done in the first call.
        //      PV generation   - Currently set as priority for battery charging
        //      Seasonal threshold - Improve approach for charge threshold to cut grid charging when more PV available
        //
        let month = simtime.current_month().ok_or_else(|| {
            anyhow!("Month could not be resolved for current simulation timestep.")
        })? as usize;
        let EnergySupplyTariffInfo {
            tariff,
            threshold_charges,
            threshold_prices,
        } = self
            .tariff_info
            .as_ref()
            .ok_or_else(|| anyhow!("Tariff info not set when expected."))?;
        let threshold_charge = threshold_charges
            .as_ref()
            .and_then(|threshold_charges| threshold_charges.get(month).copied());
        let threshold_price = threshold_prices
            .as_ref()
            .and_then(|threshold_charges| threshold_charges.get(month).copied());
        // For tariff selected look up price etc and decide whether to charge
        let elec_price = self
            .tariff_data
            .as_ref()
            .ok_or_else(|| anyhow!("Tariff data expected to be set on energy supply"))?
            .price(tariff, simtime)?;

        let current_charge = battery.get_state_of_charge();
        let charge_discharge_efficiency = battery.get_charge_discharge_efficiency();

        Ok(match threshold_price {
            Some(threshold_price) if elec_price / charge_discharge_efficiency < threshold_price => {
                match threshold_charge {
                    Some(threshold_charge) if current_charge < threshold_charge => {
                        (true, Some(threshold_charge), true)
                    }
                    _ => (false, threshold_charge, true),
                }
            }
            _ => (false, threshold_charge, false),
        })
    }

    pub(crate) fn calc_energy_import_from_grid_to_battery(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<()> {
        let mut total_charge = 0.;
        let mut total_max_capacity = 0.;
        let t_idx = simtime.index;

        // Connection-level budget: total import power is shared across all batteries
        let mut remaining_import_energy = self
            .power_limit_battery_import
            .map(|power_limit_battery_import| power_limit_battery_import * simtime.timestep);

        for battery in &self.get_batteries()? {
            // Current conditions of the battery
            let current_charge = battery.get_state_of_charge();
            let max_capacity = battery.get_max_capacity();

            if battery.is_grid_charging_possible() {
                let (charging_condition, threshold_charge, _) =
                    self.is_charging_from_grid(battery, simtime)?;

                if let Some(threshold_charge) = threshold_charge {
                    if charging_condition {
                        // Create max elec_demand from grid to complete battery charging if battery conditions allow
                        let mut elec_demand = -max_capacity * (threshold_charge - current_charge)
                            / battery.get_charge_efficiency(simtime);
                        if let Some(remaining_import_energy) = remaining_import_energy {
                            elec_demand = max_of_2(elec_demand, -remaining_import_energy);
                        };

                        // Attempt charging battery and retrieving energy_accepted
                        let energy_accepted =
                            -battery.charge_discharge_battery(elec_demand, false, simtime);

                        if let Some(remaining_import_energy) = &mut remaining_import_energy {
                            *remaining_import_energy -= energy_accepted;
                        };

                        self.energy_into_battery_from_grid[t_idx]
                            .fetch_add(energy_accepted, Ordering::SeqCst);

                        // Informing EnergyImport of imported electricity
                        self.demand_not_met[t_idx].fetch_add(energy_accepted, Ordering::SeqCst);
                    }
                };
            }

            total_charge += battery.get_state_of_charge() * max_capacity;
            total_max_capacity += max_capacity;
        }

        if relative_eq!(total_max_capacity, 0., epsilon = 1e-10, max_relative = 1e-9) {
            self.battery_state_of_charge[t_idx].store(0., Ordering::SeqCst);
        } else {
            self.battery_state_of_charge[t_idx]
                .store(total_charge / total_max_capacity, Ordering::SeqCst);
        }

        Ok(())
    }

    pub(crate) fn calc_energy_export_from_battery_to_grid(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<()> {
        if !self.is_export_capable {
            return Ok(());
        }

        let mut total_charge = 0.;
        let mut total_max_capacity = 0.;
        let t_idx = simtime.index;

        // The whole-house export power limit is shared between generation surplus and
        // battery discharge. Generation surplus has already been capped and recorded in
        // calc_energy_import_export_betafactor for this timestep, so the budget left for
        // battery discharge is the limit less the surplus already exported. supply_surplus
        // is stored negative by convention (export is negative demand), so it is added, and
        // the headroom is floored at zero for the case where surplus already met the limit.
        // Capping the discharge request (not the accepted energy) means any energy the
        // battery cannot export is retained as charge rather than discarded.
        let mut remaining_export_energy = self.power_limit_export.map(|power_limit_export| {
            max_of_2(
                0.,
                power_limit_export * simtime.timestep
                    + self.supply_surplus[t_idx].load(Ordering::SeqCst),
            )
        });

        for battery in &self.get_batteries()? {
            // Current conditions of the battery
            let current_charge = battery.get_state_of_charge();
            let max_capacity = battery.get_max_capacity();

            if battery.is_grid_exporting_possible() {
                let (exporting_condition, threshold_charge, _) =
                    self.is_exporting_to_grid(battery, simtime)?;

                if let Some(threshold_charge) = threshold_charge {
                    if exporting_condition {
                        for battery in self.get_batteries()? {
                            if self.is_charging_from_grid(&battery, simtime)?.0 {
                                bail!("Battery export conditions met while importing from grid")
                            }
                        }

                        // Create max elec_supply for battery discharge to grid if conditions allow
                        let mut elec_supply = max_capacity * (current_charge - threshold_charge)
                            / battery.get_discharge_efficiency(simtime);
                        if let Some(remaining_export_energy) = remaining_export_energy {
                            elec_supply = min_of_2(elec_supply, remaining_export_energy);
                        };

                        // Attempt charging battery and retrieving energy_accepted
                        let energy_accepted =
                            -battery.charge_discharge_battery(elec_supply, false, simtime);

                        if let Some(remaining_export_energy) = &mut remaining_export_energy {
                            // energy_accepted is negative for discharge, so this reduces the budget
                            *remaining_export_energy += energy_accepted;
                        };

                        self.energy_into_grid_from_battery[t_idx]
                            .fetch_add(energy_accepted, Ordering::SeqCst);
                    }
                };
            }

            total_charge += battery.get_state_of_charge() * max_capacity;
            total_max_capacity += max_capacity;
        }

        if relative_eq!(total_max_capacity, 0., epsilon = 1e-10, max_relative = 1e-9) {
            self.battery_state_of_charge[t_idx].store(0., Ordering::SeqCst);
        } else {
            self.battery_state_of_charge[t_idx]
                .store(total_charge / total_max_capacity, Ordering::SeqCst);
        }

        Ok(())
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

        self.energy_into_battery_from_generation
            .get(timestep_idx)
            .unwrap()
            .store(0., Ordering::SeqCst);
        self.energy_battery_to_consumption
            .get(timestep_idx)
            .unwrap()
            .store(0., Ordering::SeqCst);

        let supplies_sum =
            FSum::with_all(supplies.iter().map(|d| d.load(Ordering::SeqCst))).value();
        let demands_sum = FSum::with_all(demands.iter().map(|d| d.load(Ordering::SeqCst))).value();

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

        // If the priority order of energy surplus is specified, calculate the same according to the order
        let items_by_priority: Vec<BatteryOrDiverter> = match &self.priority {
            None => self
                .get_batteries()?
                .into_iter()
                .map(BatteryOrDiverter::from)
                .chain(
                    self.get_diverters()?
                        .into_iter()
                        .map(BatteryOrDiverter::from),
                )
                .collect(),
            Some(_) => {
                let batteries_and_diverters: IndexMap<String, BatteryOrDiverter> = self
                    .electric_batteries
                    .iter()
                    .map(|(k, v)| (k.clone(), BatteryOrDiverter::from(v.clone())))
                    .chain(
                        self.diverters
                            .iter()
                            .map(|(k, v)| (k.clone(), BatteryOrDiverter::from(v.clone()))),
                    )
                    .collect();

                self.sort_by_priority(&batteries_and_diverters)?
            }
        };

        for item in items_by_priority {
            match item {
                BatteryOrDiverter::Battery(battery) => {
                    (supply_surplus, demand_not_met) = self.charge_discharge_battery(
                        battery,
                        supply_surplus,
                        demand_not_met,
                        simtime,
                    )?
                }
                BatteryOrDiverter::Diverter(diverter) => {
                    supply_surplus = self.divert_surplus_to_pv(diverter, supply_surplus, simtime)?
                }
            }
        }

        if self.is_export_capable {
            if let Some(power_limit_export) = self.power_limit_export {
                // Cap grid export at the whole-house export power limit. Convert the power
                // limit (kW) to the maximum energy exportable in this timestep (kWh).
                // supply_surplus is negative by convention, so limiting the export magnitude
                // means taking the less negative of the surplus and the negative limit (the
                // max of the two). Any surplus above the limit is curtailed: it is neither
                // exported nor stored, having already been offered to self-consumption,
                // battery and diverter above. Battery discharge later in the timestep shares
                // this same limit (see calc_energy_export_from_battery_to_grid).
                let max_energy_export = power_limit_export * simtime.timestep;
                let surplus_before_cap = supply_surplus;
                let supply_surplus = max_of_2(supply_surplus, -max_energy_export);
                // Curtailed generation: the share that could be neither used nor exported.
                // supply_surplus is negative by convention and the cap makes it less negative,
                // so (supply_surplus - surplus_before_cap) is the reduction in export
                // magnitude, i.e. the curtailed energy (kWh, >= 0).
                self.generation_curtailed
                    .get(timestep_idx)
                    .unwrap()
                    .fetch_add(supply_surplus - surplus_before_cap, Ordering::SeqCst);
            }
            self.supply_surplus
                .get(timestep_idx)
                .unwrap()
                .fetch_add(supply_surplus, Ordering::SeqCst);
        }

        self.demand_not_met
            .get(timestep_idx)
            .unwrap()
            .fetch_add(demand_not_met, Ordering::SeqCst);
        self.grid_to_consumption
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

    /// Diverts surplus energy to a diverter and returns the new surplus
    fn divert_surplus_to_pv(
        &self,
        diverter: Arc<RwLock<dyn SurplusDiverting>>,
        supply_surplus: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let diverted = diverter.read().divert_surplus(supply_surplus, simtime)?;
        self.energy_diverted
            .get(simtime.index)
            .unwrap()
            .fetch_add(diverted, Ordering::SeqCst);

        Ok(supply_surplus + diverted)
    }

    /// Adjusts supply_surplus and demand_not_met by charging or discharging the battery. Returns the new supply_surplus, and demand_not_met
    fn charge_discharge_battery(
        &self,
        battery: Arc<ElectricBattery>,
        supply_surplus: f64,
        demand_not_met: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        // See if the battery can deal with excess supply/demand for this timestep
        // supply_surplus is -ve by convention and demand_not_met is +ve
        // TODO (from Python): assumption made here that supply is done before demand, could
        // revise in future if more evidence becomes available.
        let (charging_condition, _, can_charge_if_not_full) = if battery.is_grid_charging_possible()
        {
            self.is_charging_from_grid(&battery, simtime)?
        } else {
            (false, Default::default(), false)
        };

        let mut supply_surplus = supply_surplus;
        if supply_surplus < 0. {
            let energy_battery_to_consumption =
                battery.charge_discharge_battery(supply_surplus, charging_condition, simtime);
            supply_surplus -= energy_battery_to_consumption;
            self.energy_into_battery_from_generation[simtime.index]
                .fetch_add(-energy_battery_to_consumption, Ordering::SeqCst);
        }

        let mut demand_not_met = demand_not_met;
        if demand_not_met > 0. {
            // Calling is_charging_from_grid threshold level to avoid
            // discharging from the electric battery and
            // triggering lots of small grid recharge events
            // when the level of charge is close to the threshold
            if !can_charge_if_not_full {
                let energy_battery_to_consumption =
                    battery.charge_discharge_battery(demand_not_met, false, simtime);
                demand_not_met -= energy_battery_to_consumption;
                self.energy_battery_to_consumption[simtime.index]
                    .fetch_add(-energy_battery_to_consumption, Ordering::SeqCst);
            }
        }

        Ok((supply_surplus, demand_not_met))
    }

    /// wrapper that applies relevant function to obtain
    /// beta factor from energy supply+demand at a given timestep
    fn beta_factor_function(
        &self,
        supply: f64,
        demand: f64,
        beta_factor_function: BetaFactorFunction,
    ) -> Result<f64, NotImplementedError> {
        if relative_eq!(supply, 0., epsilon = 1e-10, max_relative = 1e-9) {
            return Ok(1.);
        }
        if relative_eq!(demand, 0., epsilon = 1e-10, max_relative = 1e-9) {
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

    #[cfg(test)]
    pub(crate) fn set_fuel_type(&mut self, fuel_type: FuelType) {
        self.fuel_type = fuel_type;
    }
}

enum BetaFactorFunction {
    Pv,
    // variant currently commented out in upstream
    #[allow(dead_code)]
    Wind,
}

pub struct EnergySupplyBuilder {
    energy_supply: EnergySupply,
}

impl EnergySupplyBuilder {
    pub fn new(fuel_type: FuelType, simulation_timesteps: usize) -> Self {
        Self {
            energy_supply: EnergySupply::new(
                fuel_type,
                simulation_timesteps,
                None,
                None,
                Default::default(),
                None,
                None,
                None,
                None,
                None,
                None,
                None,
            )
            .unwrap(),
        }
    }

    pub fn with_export_capable(mut self, is_export_capable: bool) -> Self {
        self.energy_supply.is_export_capable = is_export_capable;
        self
    }

    pub fn with_tariff_info(mut self, tariff_info: EnergySupplyTariffInfo) -> anyhow::Result<Self> {
        self.energy_supply.tariff_info = Some(tariff_info);
        Ok(self)
    }

    pub fn with_tariff_export(
        mut self,
        threshold_charges_export: [f64; 12],
        threshold_prices_export: [f64; 12],
    ) -> Self {
        self.energy_supply.tariff_export = Some(EnergySupplyTariff::ExportTariff);
        self.energy_supply.threshold_charges_export = Some(threshold_charges_export);
        self.energy_supply.threshold_prices_export = Some(threshold_prices_export);
        self
    }

    pub fn with_tariff_data(mut self, tariff_data: TariffData) -> Self {
        self.energy_supply.tariff_data = Some(tariff_data);
        self
    }

    pub fn with_electric_battery(
        mut self,
        electric_batteries: IndexMap<String, ElectricBattery>,
    ) -> Self {
        self.energy_supply.electric_batteries = electric_batteries
            .into_iter()
            .map(|(k, v)| (k, v.into()))
            .collect();
        self
    }

    pub fn with_priority(mut self, priority: Vec<impl Into<String>>) -> Self {
        self.energy_supply.priority = Some(priority.into_iter().map(|s| s.into()).collect());
        self
    }

    pub fn build(self) -> EnergySupply {
        self.energy_supply
    }

    // write other builder methods
}

#[derive(Clone)]
enum BatteryOrDiverter {
    Battery(Arc<ElectricBattery>),
    Diverter(Arc<RwLock<dyn SurplusDiverting>>),
}

impl From<Arc<RwLock<dyn SurplusDiverting>>> for BatteryOrDiverter {
    fn from(value: Arc<RwLock<dyn SurplusDiverting>>) -> Self {
        Self::Diverter(value)
    }
}

impl From<Arc<ElectricBattery>> for BatteryOrDiverter {
    fn from(value: Arc<ElectricBattery>) -> Self {
        Self::Battery(value)
    }
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
    use approx::assert_relative_eq;
    use itertools::Itertools;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use serde_json::json;
    use std::io::{BufReader, Cursor};

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0.0, 8.0, 1.0)
    }

    #[fixture]
    pub fn energy_supply<'a>(simulation_time: SimulationTime) -> EnergySupply {
        let mut energy_supply =
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time.iter().total_steps())
                .build();
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

    #[fixture]
    pub fn tariff_data(simulation_time: SimulationTime) -> TariffData {
        let prices = TariffData::load_data_from_file(BufReader::new(Cursor::new(include_str!(
            "../../../examples/tariff_data/tariff_data_25-06-2024.csv"
        ))))
        .unwrap();

        TariffData::new(
            &simulation_time.iter(),
            Some(0),
            1.,
            TariffData::expand_prices_schedule(prices).unwrap(),
        )
        .unwrap()
    }

    #[fixture]
    fn tariff_info() -> EnergySupplyTariffInfo {
        let threshold_charges = vec![0.8, 0.7, 0.7, 0.8, 0.6, 0.8, 0.7, 0.7, 0.8, 0.7, 0.8, 0.8];
        let threshold_prices = vec![16., 16., 16., 20., 20., 20., 20., 20., 20., 20., 20., 20.];

        EnergySupplyTariffInfo {
            tariff: EnergySupplyTariff::VariableTimeOfDay,
            threshold_charges: Some(threshold_charges),
            threshold_prices: Some(threshold_prices),
        }
    }

    fn create_elec_battery(
        grid_charging_possible: bool,
        grid_exporting_possible: bool,
        battery_location: BatteryLocation,
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) -> ElectricBattery {
        ElectricBattery::new(
            2.,
            0.8,
            3.,
            0.001,
            1.5,
            1.5,
            battery_location,
            grid_charging_possible,
            grid_exporting_possible,
            simulation_time.step,
            Arc::new(external_conditions),
        )
    }

    #[rstest]
    /// Tests for EnergySupply with no battery or invalid battery charging due to no tariff data
    fn test_no_battery_or_invalid_charging(
        simulation_time: SimulationTime,
        external_conditions: ExternalConditions,
        energy_supply: EnergySupply,
    ) {
        let elec_battery = create_elec_battery(
            true,
            true,
            BatteryLocation::Inside,
            external_conditions,
            simulation_time,
        );

        assert!(EnergySupply::new(
            FuelType::Electricity,
            simulation_time.iter().total_steps(),
            None,
            None,
            indexmap! {"Electric_battery".into() => elec_battery},
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .is_err());

        assert!(energy_supply.get_batteries().unwrap().is_empty());
    }

    #[rstest]
    /// Test that the state of charge is 0 where there are no batteries
    fn test_calc_energy_import_from_grid_to_battery_no_battery(
        simulation_time: SimulationTime,
        tariff_data: TariffData,
        tariff_info: EnergySupplyTariffInfo,
    ) {
        let energy_supply =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps())
                .with_tariff_data(tariff_data)
                .with_tariff_info(tariff_info)
                .unwrap()
                .build();

        energy_supply
            .calc_energy_import_from_grid_to_battery(simulation_time.iter().current_iteration())
            .unwrap();

        let (_, _, _, _, state_of_charge) = energy_supply.get_battery_energy_flows();

        assert_eq!(state_of_charge, [0.; 8]);
        assert_eq!(energy_supply.get_batteries().unwrap().len(), 0);
    }

    #[rstest]
    /// Test that the state of charge is 0 where there are no batteries
    fn test_calc_energy_export_from_battery_to_grid_no_battery(
        simulation_time: SimulationTime,
        tariff_data: TariffData,
    ) {
        let energy_supply =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps())
                .with_tariff_data(tariff_data)
                .with_tariff_export(
                    [0.8, 0.7, 0.7, 0.8, 0.6, 0.8, 0.7, 0.7, 0.8, 0.7, 0.8, 0.8],
                    [16., 16., 16., 20., 20., 20., 20., 20., 20., 20., 20., 20.],
                )
                .build();

        energy_supply
            .calc_energy_export_from_battery_to_grid(simulation_time.iter().current_iteration())
            .unwrap();

        let (_, _, _, _, state_of_charge) = energy_supply.get_battery_energy_flows();

        assert_eq!(state_of_charge, [0.; 8]);
        assert_eq!(energy_supply.get_batteries().unwrap().len(), 0);
    }

    #[rstest]
    fn test_is_exporting_to_grid(
        simulation_time: SimulationTime,
        external_conditions: ExternalConditions,
        tariff_data: TariffData,
    ) {
        let elec_battery = create_elec_battery(
            true,
            true,
            BatteryLocation::Inside,
            external_conditions.clone(),
            simulation_time,
        );

        let builder =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.iter().total_steps());
        let energy_supply = builder
            .with_electric_battery(indexmap! {"battery".into() => elec_battery})
            .with_tariff_data(tariff_data.clone())
            .with_tariff_export([0.8; 12], [5.; 12])
            .build();

        let battery = &energy_supply.electric_batteries[0];
        battery.charge_discharge_battery(-100., false, simulation_time.iter().current_iteration());

        // Meets threshold charge
        let (exporting_condition, threshold_charge, can_export_if_not_empty) = energy_supply
            .is_exporting_to_grid(battery, simulation_time.iter().current_iteration())
            .unwrap();

        assert!(exporting_condition);
        assert_eq!(threshold_charge, Some(0.8));
        assert!(can_export_if_not_empty);

        let elec_battery = create_elec_battery(
            true,
            true,
            BatteryLocation::Inside,
            external_conditions.clone(),
            simulation_time,
        );

        let builder =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.iter().total_steps());
        let energy_supply = builder
            .with_electric_battery(indexmap! {"battery".into() => elec_battery})
            .with_tariff_data(tariff_data.clone())
            .with_tariff_export([0.9; 12], [5.; 12])
            .build();

        let battery = &energy_supply.electric_batteries[0];

        // Meets threshold price
        let (exporting_condition, threshold_charge, can_export_if_not_empty) = energy_supply
            .is_exporting_to_grid(battery, simulation_time.iter().current_iteration())
            .unwrap();

        assert!(!exporting_condition);
        assert_eq!(threshold_charge, Some(0.9));
        assert!(can_export_if_not_empty);

        let elec_battery = create_elec_battery(
            true,
            true,
            BatteryLocation::Inside,
            external_conditions,
            simulation_time,
        );

        let builder =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.iter().total_steps());
        let energy_supply = builder
            .with_electric_battery(indexmap! {"battery".into() => elec_battery})
            .with_tariff_data(tariff_data)
            .with_tariff_export([0.8; 12], [20.; 12])
            .build();

        let battery = &energy_supply.electric_batteries[0];

        // can't export
        let (exporting_condition, threshold_charge, can_export_if_not_empty) = energy_supply
            .is_exporting_to_grid(battery, simulation_time.iter().current_iteration())
            .unwrap();

        assert!(!exporting_condition);
        assert_eq!(threshold_charge, Some(0.8));
        assert!(!can_export_if_not_empty);
    }

    #[rstest]
    /// Test that calc_energy_export_from_battery_to_grid doesn't export if the energy supply is not export capable
    fn test_calc_energy_export_from_battery_to_grid_not_export_capable(
        tariff_data: TariffData,
        simulation_time: SimulationTime,
        external_conditions: ExternalConditions,
    ) {
        let elec_battery = create_elec_battery(
            true,
            true,
            BatteryLocation::Inside,
            external_conditions,
            simulation_time,
        );

        let builder =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.iter().total_steps());
        let energy_supply = builder
            .with_electric_battery(indexmap! {"battery".into() => elec_battery})
            .with_tariff_data(tariff_data)
            .with_tariff_export([0.8; 12], [5.; 12])
            .with_export_capable(false)
            .build();

        energy_supply
            .calc_energy_export_from_battery_to_grid(simulation_time.iter().current_iteration())
            .unwrap();

        let (_, _, _, energy_into_grid_from_battery, _) = energy_supply.get_battery_energy_flows();

        assert_eq!(energy_into_grid_from_battery[0], 0.);
    }

    #[rstest]
    fn test_calc_energy_export_from_battery_to_grid(
        tariff_data: TariffData,
        simulation_time: SimulationTime,
        external_conditions: ExternalConditions,
    ) {
        let elec_battery = create_elec_battery(
            true,
            true,
            BatteryLocation::Inside,
            external_conditions,
            simulation_time,
        );
        elec_battery.charge_discharge_battery(
            -10.,
            false,
            simulation_time.iter().current_iteration(),
        );

        let builder =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.iter().total_steps());
        let energy_supply = builder
            .with_electric_battery(indexmap! {"battery".into() => elec_battery})
            .with_tariff_data(tariff_data)
            .with_tariff_export([0.8; 12], [5.; 12])
            .build();

        for (_, t_it) in simulation_time.iter().enumerate() {
            energy_supply
                .calc_energy_export_from_battery_to_grid(t_it)
                .unwrap();
            energy_supply.timestep_end().unwrap();
        }
    }

    #[rstest]
    fn test_existing_user_name(energy_supply: EnergySupply) {
        assert!(EnergySupply::connection(Arc::new(energy_supply.into()), "shower").is_err());
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

            assert_eq!(energy_supply.demand_total[t_idx].load(Ordering::SeqCst), 0.);
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
        assert!(energy_supply.diverters.is_empty());
        energy_supply
            .connect_diverter(pv_diverter.clone(), None)
            .unwrap();
        assert!(!energy_supply.diverters.is_empty());
        assert!(energy_supply
            .connect_diverter(pv_diverter.clone(), None)
            .is_err());
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
            assert!(energy_supply
                .demand_energy("others", amount_demanded[t_idx], t_idx)
                .is_err());
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
    pub fn test_results_by_end_user_and_step(
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
            assert_eq!(
                energy_supply
                    .read()
                    .results_by_end_user_single_step(simtime.index),
                IndexMap::from([
                    (
                        energy_connection_1.clone().end_user_name,
                        EXPECTED_TOTAL_DEMANDS_BY_END_USER[0][simtime.index]
                    ),
                    (
                        energy_connection_2.clone().end_user_name,
                        EXPECTED_TOTAL_DEMANDS_BY_END_USER[1][simtime.index]
                    )
                ])
            );

            // Case where end_user is not in energy_out_by_end_user,
            // (demand_by_end_user.keys() != energy_out_by_end_user.keys())
            energy_supply.write().energy_out_by_end_user.insert(
                "others".into(),
                vec![EXPECTED_TOTAL_DEMANDS_BY_END_USER[0][simtime.index].into()],
            );
            assert_eq!(
                energy_supply.read().results_by_end_user()["shower"][simtime.index],
                EXPECTED_TOTAL_DEMANDS_BY_END_USER[0][simtime.index]
            );

            // Testing the edge case at the last timestep,
            // to check when an end_user exists only in demand_by_end_user, not in energy_out_by_end_user
            if simtime.index == 7 {
                energy_supply.write().demand_by_end_user.insert(
                    "others1".into(),
                    [0.0; 8].into_iter().map(AtomicF64::new).collect::<Vec<_>>(),
                );
                energy_supply.write().demand_by_end_user["others1"][simtime.index] =
                    ((simtime.index as f64 + 1.0) * 50.0).into();

                assert_eq!(
                    energy_supply
                        .read()
                        .results_by_end_user_single_step(simtime.index),
                    IndexMap::from([
                        (
                            energy_connection_1.clone().end_user_name,
                            EXPECTED_TOTAL_DEMANDS_BY_END_USER[0][simtime.index]
                        ),
                        (
                            energy_connection_2.clone().end_user_name,
                            EXPECTED_TOTAL_DEMANDS_BY_END_USER[1][simtime.index]
                        ),
                        ("others1".into(), ((simtime.index as f64 + 1.0) * 50.0))
                    ])
                );
            }
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

        // When beta_factor_function is not PV (not captured when calling get_beta_factor())
        assert!(energy_supply
            .read()
            .beta_factor_function(1., 1., BetaFactorFunction::Wind)
            .is_err());
        // When there is no demand, beta_factor_function returns 0
        assert_eq!(
            energy_supply
                .read()
                .beta_factor_function(5., 0., BetaFactorFunction::Pv)
                .unwrap(),
            0.
        );
    }

    #[rstest]
    fn test_battery_with_grid_charging_and_priority(
        simulation_time: SimulationTime,
        external_conditions: ExternalConditions,
        tariff_info: EnergySupplyTariffInfo,
        tariff_data: TariffData,
    ) {
        // Valid battery where there is grid charging and tariff_path is set
        let battery_age = 3.;
        let elec_battery = create_elec_battery(
            true,
            true,
            BatteryLocation::Inside,
            external_conditions,
            simulation_time,
        );
        let builder =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.iter().total_steps());
        let energy_supply = builder
            .with_electric_battery(indexmap! {"ElectricBattery".into() => elec_battery})
            .with_tariff_info(tariff_info)
            .unwrap()
            .with_tariff_data(tariff_data)
            .with_priority(vec!["ElectricBattery", "diverter"])
            .build();

        assert!(energy_supply.tariff_data.is_some());
        assert_eq!(energy_supply.get_batteries().unwrap().len(), 1);
        assert!(energy_supply.has_battery().unwrap());

        let battery_state_of_health = -0.04 * battery_age + 1.;

        assert_eq!(
            energy_supply.get_battery_max_capacity().unwrap().unwrap(),
            2. * battery_state_of_health
        ); // max capacity * state of health
        assert_eq!(
            energy_supply
                .get_battery_charge_efficiency(simulation_time.iter().current_iteration(), None)
                .unwrap()
                .unwrap(),
            0.8_f64.powf(0.5) * battery_state_of_health * 1.
        ); // one way efficiency * state of health * air_temp_capacity_factor
        assert_eq!(
            energy_supply
                .get_battery_discharge_efficiency(simulation_time.iter().current_iteration(), None)
                .unwrap()
                .unwrap(),
            0.8_f64.powf(0.5) * battery_state_of_health * 1.
        ); // one way efficiency * state of health * air_temp_capacity_factor
        assert_eq!(
            energy_supply
                .get_battery_max_discharge(0.7, None)
                .unwrap()
                .unwrap(),
            -(1.5 * 1.)
        ); // max discharge rate * discharge factor * timestep * -1
        assert_eq!(
            energy_supply
                .get_battery_available_charge()
                .unwrap()
                .unwrap(),
            0.
        );

        let expected_charging_state = [
            (true, Some(0.8), true), // elec_price/efficiency=13.5877158875 < 16; current charge=0
            (false, Some(0.8), true), // elec_price/efficiency=12.16550081375 < 16;
            (false, Some(0.8), false), // elec_price/efficiency=18.1486140875 !< 16
            (false, Some(0.8), true), // elec_price/efficiency=11.6733265175 < 16
            (false, Some(0.8), false), // elec_price/efficiency=25.5426676375 !< 16
            (false, Some(0.8), false), // elec_price/efficiency=24.914735325 !< 16
            (false, Some(0.8), false), // elec_price/efficiency=24.1577282 !< 16
            (false, Some(0.8), false), // elec_price/efficiency=17.394483375 !< 16
        ];
        let expected_energy_import_from_grid = vec![1.788854381999832, 0., 0., 0., 0., 0., 0., 0.];
        let expected_energy_export_to_grid = vec![0.; 8];
        let expected_diverted_energy = [0.; 8];
        let expected_generated_energy_into_battery = vec![0.; 8];
        let expected_energy_out_of_battery = vec![0.; 8];
        let expected_battery_state_of_charge = vec![0.8; 8];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let battery = energy_supply.get_batteries().unwrap()[0].clone();

            assert_eq!(
                energy_supply.is_charging_from_grid(&battery, t_it).unwrap(),
                expected_charging_state[t_idx]
            );
            energy_supply
                .calc_energy_import_from_grid_to_battery(t_it)
                .unwrap();
            energy_supply.timestep_end().unwrap();

            energy_supply
                .calc_energy_import_export_betafactor(t_it)
                .unwrap();
        }

        assert_eq!(
            energy_supply.get_energy_diverted(),
            expected_diverted_energy
        );
        assert_eq!(
            energy_supply.get_battery_energy_flows(),
            (
                expected_generated_energy_into_battery,
                expected_energy_out_of_battery,
                expected_energy_import_from_grid,
                expected_energy_export_to_grid,
                expected_battery_state_of_charge,
            )
        );
    }

    #[rstest]
    fn test_battery_with_grid_charging_no_priority(
        simulation_time: SimulationTime,
        external_conditions: ExternalConditions,
        tariff_info: EnergySupplyTariffInfo,
        tariff_data: TariffData,
    ) {
        let elec_battery = create_elec_battery(
            true,
            true,
            BatteryLocation::Inside,
            external_conditions,
            simulation_time,
        );
        let builder =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.iter().total_steps());
        let energy_supply = builder
            .with_electric_battery(indexmap! {"Electric_battery".into() => elec_battery})
            .with_tariff_info(tariff_info)
            .unwrap()
            .with_tariff_data(tariff_data)
            .build();

        for t_idx in simulation_time.iter() {
            energy_supply
                .calc_energy_import_export_betafactor(t_idx)
                .unwrap();
        }

        assert_eq!(
            energy_supply.get_battery_energy_flows(),
            (
                vec![0.; 8],
                vec![0.; 8],
                vec![0.; 8],
                vec![0.; 8],
                vec![0.; 8],
            )
        )
    }

    #[rstest]
    fn test_battery_without_grid_charging(
        simulation_time: SimulationTime,
        external_conditions: ExternalConditions,
    ) {
        let elec_battery = create_elec_battery(
            false,
            true,
            BatteryLocation::Inside,
            external_conditions,
            simulation_time,
        );
        let builder =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.iter().total_steps());
        let energy_supply = builder
            .with_electric_battery(indexmap! {"Electric_battery".into() => elec_battery})
            .build();

        for t_idx in simulation_time.iter() {
            energy_supply
                .calc_energy_import_from_grid_to_battery(t_idx)
                .unwrap();
            energy_supply.timestep_end().unwrap();
        }

        assert_eq!(
            energy_supply.get_battery_energy_flows(),
            (
                vec![0.; 8],
                vec![0.; 8],
                vec![0.; 8],
                vec![0.; 8],
                vec![0.; 8],
            )
        )
    }

    #[fixture]
    fn external_conditions(simulation_time: SimulationTime) -> ExternalConditions {
        ExternalConditions::new(
            &simulation_time.iter(),
            vec![0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0],
            vec![3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1],
            vec![0., 20., 40., 60., 0., 20., 40., 60.]
                .into_iter()
                .map(Into::into)
                .collect(),
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
                {"start360": 0, "end360": 45},
                {"start360": 45, "end360": 90,
                 "shading": [
                     {"type": "overhang", "height": 2.2, "distance": 6}
                     ]
                 },
                {"start360": 90, "end360": 135},
                {"start360": 135, "end360": 180,
                 "shading": [
                     {"type": "obstacle", "height": 40, "distance": 4},
                     {"type": "overhang", "height": 3, "distance": 7}
                     ]
                 },
                {"start360": 180, "end360": 225,
                 "shading": [
                     {"type": "obstacle", "height": 3, "distance": 8},
                     ]
                 },
                {"start360": 225, "end360": 270},
                {"start360": 270, "end360": 315},
                {"start360": 315, "end360": 360}
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

        let elec_battery = create_elec_battery(
            false,
            true,
            BatteryLocation::Outside,
            external_conditions,
            simulation_time,
        );

        let builder =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.iter().total_steps());
        let energy_supply = builder
            .with_electric_battery(indexmap! {"ElectricBattery".into() => elec_battery})
            .build();

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
                    15.256924949254948,
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
                    -14.016897653541719,
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
                energy_supply.get_energy_generated_consumed(),
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
                energy_supply.get_grid_to_consumption(),
                vec![
                    15.256924949254948,
                    34.59889302603037,
                    52.991809292243516,
                    63.07009986960006,
                    -2.842170943040401e-14,
                    99.5880926222444,
                    124.3891811736907,
                    140.57522007095577,
                ]
            );

            assert_eq!(
                energy_supply.get_battery_energy_flows(),
                (
                    vec![2.2431023464582824, 0., 0., 0., 0., 0., 0., 0.,],
                    vec![-1.0030750507450579, -0., 0., 0., 0., 0., 0., 0.,],
                    vec![0.; 8],
                    vec![0.; 8],
                    vec![0.; 8]
                )
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
        energy_supply
            .write()
            .connect_diverter(diverter, Some("diverter".into()))
            .unwrap();

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
                    47.77692494925496,
                    103.7966790780911,
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
                    -36.53689765354172,
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
                    .energy_battery_to_consumption
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
        let elec_batteries = Arc::into_inner(energy_supply)
            .unwrap()
            .into_inner()
            .electric_batteries
            .into_iter()
            .map(|(k, v)| (k, Arc::into_inner(v).unwrap()))
            .collect();

        // Set priority
        let priority = vec!["diverter", "ElectricBattery"];

        let mut builder =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.iter().total_steps());
        builder = builder
            .with_electric_battery(elec_batteries)
            .with_priority(priority);

        let energy_supply = Arc::new(RwLock::new(builder.build()));

        let _shower_connection = EnergySupply::connection(energy_supply.clone(), "shower").unwrap();
        let _bath_connection = EnergySupply::connection(energy_supply.clone(), "bath").unwrap();

        let diverter = Arc::new(RwLock::new(MockDiverter));
        energy_supply
            .write()
            .connect_diverter(diverter, None)
            .unwrap();

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
                    .energy_battery_to_consumption
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

    // skipping python's test_is_charging_from_grid_exception as in rust is_charging_from_grid requires a battery to be passed in

    #[rstest]
    /// Test that the battery state of charge is correct after calling calc_energy_import_from_grid_to_battery
    fn test_calc_energy_import_from_grid_to_battery_multiple_items(
        external_conditions: ExternalConditions,
    ) {
        let simtime = SimulationTime::new(0.0, 4.0, 1.0);
        let battery_a = create_elec_battery(
            false,
            true,
            BatteryLocation::Outside,
            external_conditions.clone(),
            simtime,
        );
        let battery_b = create_elec_battery(
            false,
            true,
            BatteryLocation::Outside,
            external_conditions,
            simtime,
        );

        battery_a.charge_discharge_battery(-3., false, simtime.iter().current_iteration());
        battery_b.charge_discharge_battery(-2., false, simtime.iter().current_iteration());

        let mut builder =
            EnergySupplyBuilder::new(FuelType::MainsGas, simtime.iter().total_steps());
        builder = builder
            .with_electric_battery(indexmap! {"A".into() => battery_a, "B".into() => battery_b})
            .with_priority(vec!["B", "A"]);
        let energy_supply = Arc::new(RwLock::new(builder.build()));

        let shower_connection = EnergySupply::connection(energy_supply.clone(), "shower").unwrap();
        shower_connection
            .demand_energy(100., simtime.iter().current_index())
            .unwrap();
        shower_connection
            .supply_energy(120., simtime.iter().current_index())
            .unwrap();

        energy_supply
            .read()
            .calc_energy_import_from_grid_to_battery(simtime.iter().current_iteration())
            .unwrap();

        assert_relative_eq!(
            energy_supply.read().get_battery_energy_flows().4[0],
            0.806089034373128
        );
        assert_relative_eq!(
            energy_supply
                .read()
                .electric_batteries
                .get("A")
                .unwrap()
                .get_state_of_charge(),
            0.8522727272727273
        );
        assert_relative_eq!(
            energy_supply
                .read()
                .electric_batteries
                .get("B")
                .unwrap()
                .get_state_of_charge(),
            0.7599053414735286
        );
    }

    // skipping python's test_calc_energy_import_export_betafactor_multiple_items as mocking/assertions difficult to replicate in rust

    #[rstest]
    fn test_sort_by_priority(energy_supply: EnergySupply, simulation_time: SimulationTime) {
        assert_eq!(
            energy_supply
                .sort_by_priority(&indexmap! {"B".into() => 2, "A".into() => 1, "C".into() => 3})
                .unwrap(),
            vec![2, 1, 3]
        );

        let mut builder =
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time.total_steps());
        builder = builder.with_priority(vec!["A", "C", "B"]);
        let energy_supply = builder.build();

        assert_eq!(
            energy_supply
                .sort_by_priority(&indexmap! {"B".into() => 2, "A".into() => 1, "C".into() => 3})
                .unwrap(),
            vec![1, 3, 2]
        );
        assert_eq!(
            energy_supply
                .sort_by_priority(&indexmap! {"A".into() => 1, "C".into() => 2})
                .unwrap(),
            vec![1, 2]
        );
        assert!(energy_supply
            .sort_by_priority(&indexmap! {"D".into() => 2})
            .is_err());
        assert!(energy_supply
            .sort_by_priority(
                &indexmap! {"A".into() => 1, "B".into() => 2, "C".into() => 3, "D".into() => 4}
            )
            .is_err());
    }

    #[rstest]
    // in python this test checks that timestep_end is called on the mock battery - difficult to replicate
    // in rust so instead checking the result of the function being called (total_time_charging_current_timestep reset to 0)
    fn test_timestep_end(simulation_time: SimulationTime, external_conditions: ExternalConditions) {
        let battery = create_elec_battery(
            false,
            false,
            BatteryLocation::Outside,
            external_conditions,
            simulation_time,
        );
        battery.set_total_time_charging_current_timestep(10.);
        let builder =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps())
                .with_electric_battery(indexmap! {"ElectricBattery".into() => battery})
                .with_priority(vec!["diverter", "ElectricBattery"]);

        let energy_supply = builder.build();

        assert_eq!(
            energy_supply
                .electric_batteries
                .get("ElectricBattery")
                .unwrap()
                .get_total_time_charging_current_timestep(),
            10.
        );

        energy_supply.timestep_end().unwrap();

        assert_eq!(
            energy_supply
                .electric_batteries
                .get("ElectricBattery")
                .unwrap()
                .get_total_time_charging_current_timestep(),
            0.
        );
    }

    #[rstest]
    fn test_no_battery(simulation_time: SimulationTime) {
        let energy_supply =
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build();

        assert!(!energy_supply.has_battery().unwrap());
        assert!(energy_supply.get_battery_max_capacity().unwrap().is_none());
        assert!(energy_supply
            .get_battery_available_charge()
            .unwrap()
            .is_none());
    }

    #[rstest]
    pub fn test_energy_supply_without_export(simulation_time: SimulationTime) {
        let mut builder =
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time.iter().total_steps());
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
