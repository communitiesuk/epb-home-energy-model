use crate::compare_floats::{max_of_2, min_of_2};
use crate::core::common::WaterSourceWithTemperature;
use crate::core::controls::time_control::{Control, ControlBehaviour};
use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::core::material_properties::WATER;
use crate::core::units::{DAYS_PER_YEAR, HOURS_PER_DAY, WATTS_PER_KILOWATT};
use crate::core::water_heat_demand::dhw_demand::{DemandVolTargetKey, VolumeReference};
use crate::external_conditions::ExternalConditions;
use crate::input::{BoilerHotWaterTest, HotWaterSourceDetails};
use crate::input::{EnergySupplyType, HeatSourceLocation, HeatSourceWetDetails};
use crate::simulation_time::SimulationTimeIteration;
use crate::statistics::np_interp;
use anyhow::bail;
use arrayvec::ArrayString;
use atomic_float::AtomicF64;
use indexmap::IndexMap;
use parking_lot::RwLock;
use std::borrow::Cow;
use std::collections::HashMap;
use std::fmt;
use std::sync::atomic::Ordering;
use std::sync::Arc;

#[derive(Clone, Copy, PartialEq)]
pub enum ServiceType {
    WaterCombi,
    WaterRegular,
    Space,
}

#[derive(Debug)]
pub struct BoilerServiceWaterCombi {
    boiler: Arc<RwLock<Boiler>>,
    service_name: String,
    temperature_hot_water_in_c: f64,
    cold_feed: WaterSourceWithTemperature,
    separate_dhw_tests: BoilerHotWaterTest,
    rejected_energy_1: Option<f64>,
    storage_loss_factor_2: Option<f64>,
    rejected_factor_3: Option<f64>,
    daily_hot_water_usage: f64,
    simulation_timestep: f64,
    combi_loss: AtomicF64,
}

#[derive(Debug)]
pub struct IncorrectBoilerDataType;

impl fmt::Display for IncorrectBoilerDataType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Incorrect boiler data type provided (expected details for a combi boiler)"
        )
    }
}

impl std::error::Error for IncorrectBoilerDataType {}

impl BoilerServiceWaterCombi {
    pub fn new(
        boiler: Arc<RwLock<Boiler>>,
        boiler_data: HotWaterSourceDetails,
        service_name: String,
        temperature_hot_water_in_c: f64,
        cold_feed: WaterSourceWithTemperature,
        simulation_timestep: f64,
    ) -> Result<Self, IncorrectBoilerDataType> {
        match boiler_data {
            HotWaterSourceDetails::CombiBoiler {
                separate_dhw_tests,
                rejected_energy_1,
                storage_loss_factor_2,
                rejected_factor_3,
                daily_hot_water_usage,
                ..
            } => {
                let (rejected_energy_1, storage_loss_factor_2, rejected_factor_3) =
                    match separate_dhw_tests {
                        BoilerHotWaterTest::ML | BoilerHotWaterTest::MS => (
                            Some(rejected_energy_1),
                            Some(storage_loss_factor_2),
                            Some(rejected_factor_3),
                        ),
                        BoilerHotWaterTest::MOnly => {
                            (Some(rejected_energy_1), Some(storage_loss_factor_2), None)
                        }
                        BoilerHotWaterTest::NoAdditionalTests => (None, None, None),
                    };
                Ok(Self {
                    boiler,
                    service_name,
                    temperature_hot_water_in_c,
                    separate_dhw_tests,
                    rejected_energy_1,
                    storage_loss_factor_2,
                    rejected_factor_3,
                    daily_hot_water_usage,
                    cold_feed,
                    simulation_timestep,
                    combi_loss: Default::default(),
                })
            }
            _ => Err(IncorrectBoilerDataType),
        }
    }

    pub fn get_cold_water_source(&self) -> &WaterSourceWithTemperature {
        &self.cold_feed
    }

    pub fn temperature_hot_water_in_c(&self) -> f64 {
        self.temperature_hot_water_in_c
    }

    pub fn demand_hot_water(
        &self,
        volume_demanded_target: IndexMap<DemandVolTargetKey, VolumeReference>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let timestep = self.simulation_timestep;
        let return_temperature = 60.;

        let energy_content_kwh_per_litre = WATER.volumetric_energy_content_kwh_per_litre(
            self.temperature_hot_water_in_c,
            self.cold_feed.temperature(simtime),
        );

        let volume_demanded = volume_demanded_target
            .get(&DemandVolTargetKey::TempHotWater)
            .map(|volume_reference| volume_reference.warm_vol)
            .unwrap_or(0.0);

        let mut energy_demand = volume_demanded * energy_content_kwh_per_litre;

        let combi_loss = self.boiler_combi_loss(energy_demand, timestep);
        energy_demand += combi_loss;

        self.boiler
            .write()
            .demand_energy(
                &self.service_name,
                ServiceType::WaterCombi,
                energy_demand,
                return_temperature,
                None,
                None,
                simtime,
            )
            .map(|res| res.0)
    }

    fn boiler_combi_loss(&self, energy_demand: f64, timestep: f64) -> f64 {
        // daily hot water usage factor
        let threshold_volume = 100.;
        let fu = if self.daily_hot_water_usage < threshold_volume {
            self.daily_hot_water_usage / threshold_volume
        } else {
            1.
        };

        // Equivalent hot water litres at 60C for HW load profiles
        let hw_litres_s_profile = 36.0;
        let hw_litres_m_profile = 100.2;
        let hw_litres_l_profile = 199.8;

        let daily_vol_factor = if matches!(self.separate_dhw_tests, BoilerHotWaterTest::MS)
            && self.daily_hot_water_usage < hw_litres_s_profile
        {
            64.2
        } else if (matches!(self.separate_dhw_tests, BoilerHotWaterTest::ML)
            && self.daily_hot_water_usage < hw_litres_m_profile)
            || (matches!(self.separate_dhw_tests, BoilerHotWaterTest::MS)
                && self.daily_hot_water_usage > hw_litres_m_profile)
        {
            0.
        } else if matches!(self.separate_dhw_tests, BoilerHotWaterTest::ML)
            && self.daily_hot_water_usage > hw_litres_l_profile
        {
            -99.6
        } else {
            hw_litres_m_profile - self.daily_hot_water_usage
        };

        let combi_loss = match self.separate_dhw_tests {
            BoilerHotWaterTest::ML | BoilerHotWaterTest::MS => {
                (energy_demand
                    * (self.rejected_energy_1.unwrap()
                        + daily_vol_factor * self.rejected_factor_3.unwrap()))
                    * fu
                    + self.storage_loss_factor_2.unwrap() * (timestep / HOURS_PER_DAY as f64)
            }
            BoilerHotWaterTest::MOnly => {
                (energy_demand * self.rejected_energy_1.unwrap()) * fu
                    + self.storage_loss_factor_2.unwrap() * (timestep / HOURS_PER_DAY as f64)
            }
            BoilerHotWaterTest::NoAdditionalTests => {
                let default_combi_loss = 600;
                (default_combi_loss / DAYS_PER_YEAR) as f64 * (timestep / HOURS_PER_DAY as f64)
            }
        };

        self.combi_loss.store(combi_loss, Ordering::SeqCst);

        combi_loss
    }

    pub fn internal_gains(&self) -> f64 {
        // TODO (from the Python) Fraction of hot water energy resulting in internal gains should
        // ideally be defined in one place, but it is duplicated here and in
        // main hot water demand calculation for now.
        let frac_dhw_energy_internal_gains = 0.25;
        let gain_internal = frac_dhw_energy_internal_gains
            * self.combi_loss.load(Ordering::SeqCst)
            * WATTS_PER_KILOWATT as f64
            / self.simulation_timestep;

        self.combi_loss.store(Default::default(), Ordering::SeqCst);

        gain_internal
    }

    pub fn energy_output_max(&self) -> f64 {
        self.boiler
            .read()
            .energy_output_max(self.temperature_hot_water_in_c, None)
    }

    pub fn is_on(&self) -> bool {
        // always true as there is no associated control
        true
    }
}

#[derive(Debug)]
pub struct BoilerServiceWaterRegular {
    boiler: Arc<RwLock<Boiler>>,
    service_name: String,
    temperature_hot_water_in_c: f64,
    control: Option<Arc<Control>>,
}

impl BoilerServiceWaterRegular {
    pub fn new(
        boiler: Arc<RwLock<Boiler>>,
        service_name: String,
        temperature_hot_water_in_c: f64,
        control: Option<Arc<Control>>,
    ) -> Self {
        Self {
            boiler,
            service_name,
            temperature_hot_water_in_c,
            control,
        }
    }

    pub fn demand_energy(
        &mut self,
        mut energy_demand: f64,
        temp_return: f64,
        hybrid_service_bool: Option<bool>,
        time_elapsed_hp: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, Option<f64>)> {
        let hybrid_service_bool = hybrid_service_bool.unwrap_or(false);

        if !self.is_on(simtime) {
            energy_demand = 0.;
        }

        self.boiler.write().demand_energy(
            &self.service_name,
            ServiceType::WaterRegular,
            energy_demand,
            temp_return,
            Some(hybrid_service_bool),
            time_elapsed_hp,
            simtime,
        )
    }

    pub fn energy_output_max(
        &self,
        _temp_return: f64,
        time_elapsed_hp: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        if !self.is_on(simtime) {
            return 0.;
        }

        self.boiler
            .read()
            .energy_output_max(self.temperature_hot_water_in_c, time_elapsed_hp)
    }

    fn is_on(&self, simtime: SimulationTimeIteration) -> bool {
        match &self.control {
            Some(c) => c.is_on(simtime),
            None => true,
        }
    }
}

/// A struct representing a space heating service provided by a boiler to e.g. a cylinder.
#[derive(Debug)]
pub struct BoilerServiceSpace {
    boiler: Arc<RwLock<Boiler>>,
    service_name: String,
    control: Arc<Control>,
}

impl BoilerServiceSpace {
    pub fn new(boiler: Arc<RwLock<Boiler>>, service_name: String, control: Arc<Control>) -> Self {
        Self {
            boiler,
            service_name,
            control,
        }
    }

    pub fn temp_setpnt(&self, simtime: SimulationTimeIteration) -> Option<f64> {
        self.control.setpnt(&simtime)
    }

    pub fn in_required_period(&self, simtime: SimulationTimeIteration) -> Option<bool> {
        self.control.in_required_period(&simtime)
    }

    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        _temp_flow: f64,
        temp_return: f64,
        hybrid_service_bool: Option<bool>,
        time_elapsed_hp: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, Option<f64>)> {
        let hybrid_service_bool = hybrid_service_bool.unwrap_or(false);

        if !self.is_on(simtime) {
            return Ok(if !hybrid_service_bool {
                (0.0, None)
            } else {
                (0.0, Some(0.0))
            });
        }

        self.boiler.write().demand_energy(
            &self.service_name,
            ServiceType::Space,
            energy_demand,
            temp_return,
            Some(hybrid_service_bool),
            time_elapsed_hp,
            simtime,
        )
    }

    pub fn energy_output_max(
        &self,
        temp_output: f64,
        _temp_return_feed: f64,
        time_elapsed_hp: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        if !self.is_on(simtime) {
            0.0
        } else {
            self.boiler
                .read()
                .energy_output_max(temp_output, time_elapsed_hp)
        }
    }

    fn is_on(&self, simtime: SimulationTimeIteration) -> bool {
        self.control.is_on(simtime)
    }
}

#[derive(Debug)]
pub struct Boiler {
    energy_supply: Arc<RwLock<EnergySupply>>,
    simulation_timestep: f64,
    external_conditions: Arc<ExternalConditions>,
    energy_supply_connections: HashMap<String, EnergySupplyConnection>,
    energy_supply_connection_aux: EnergySupplyConnection,
    energy_supply_type: EnergySupplyType,
    // service_results: (),
    boiler_location: HeatSourceLocation,
    min_modulation_load: f64,
    boiler_power: f64,
    // fuel_code: (), // derived from energy supply
    power_circ_pump: f64,
    power_part_load: f64,
    power_full_load: f64,
    power_standby: f64,
    // some kind of memoisation? review this
    total_time_running_current_timestep: AtomicF64,
    corrected_full_load_gross: f64,
    room_temperature: f64,
    temp_rise_standby_loss: f64,
    standby_loss_index: f64,
    ebv_curve_offset: f64,
    service_results: RwLock<Vec<ServiceResult>>,
}

impl Boiler {
    /// Arguments:
    /// * `boiler_data` - boiler characteristics
    /// * `external_conditions` - reference to an ExternalConditions value
    pub fn new(
        boiler_data: HeatSourceWetDetails,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_conn_aux: EnergySupplyConnection,
        external_conditions: Arc<ExternalConditions>,
        simulation_timestep: f64,
    ) -> anyhow::Result<Self> {
        match boiler_data {
            HeatSourceWetDetails::Boiler {
                energy_supply: energy_supply_type,
                energy_supply_auxiliary: _,
                rated_power: boiler_power,
                efficiency_full_load: full_load_gross,
                efficiency_part_load: part_load_gross,
                boiler_location,
                // NB. there is a validation check in the Python here that modulation_load <= 1 - in this project the value has already been validated on the way in
                modulation_load: min_modulation_load,
                electricity_circ_pump: power_circ_pump,
                electricity_part_load: power_part_load,
                electricity_full_load: power_full_load,
                electricity_standby: power_standby,
            } => {
                let total_time_running_current_timestep = 0.;

                let net_to_gross = Self::net_to_gross(energy_supply_type)?;
                let full_load_net = full_load_gross / net_to_gross;
                let part_load_net = part_load_gross / net_to_gross;
                let corrected_full_load_net = Self::high_value_correction_full_load(full_load_net);
                let corrected_part_load_net =
                    Self::high_value_correction_part_load(energy_supply_type, part_load_net)?;
                let corrected_full_load_gross = corrected_full_load_net * net_to_gross;
                let corrected_part_load_gross = corrected_part_load_net * net_to_gross;

                // SAP model properties
                let room_temperature = 19.5; // TODO in the Python there is a todo to make this less hard-coded

                // 30 is the nominal temperature difference between boiler and test room
                // during standby loss test (EN15502-1 or EN15034)
                let temp_rise_standby_loss = 30.;
                // boiler standby heat loss power law index
                let standby_loss_index = 1.25;

                // Calculate offset for EBV curves
                let average_measured_eff =
                    (corrected_part_load_gross + corrected_full_load_gross) / 2.;
                // test conducted at return temperature 30C
                let temp_part_load_test = 30.;
                // test conducted at return temperature 60C
                let temp_full_load_test = 60.;
                let offset_for_theoretical_eff = 0.;
                let theoretical_eff_part_load = Self::efficiency_over_return_temperatures(
                    energy_supply_type,
                    temp_part_load_test,
                    offset_for_theoretical_eff,
                )?;
                let theoretical_eff_full_load = Self::efficiency_over_return_temperatures(
                    energy_supply_type,
                    temp_full_load_test,
                    offset_for_theoretical_eff,
                )?;
                let average_theoretical_eff =
                    (theoretical_eff_part_load + theoretical_eff_full_load) / 2.;
                let ebv_curve_offset = average_theoretical_eff - average_measured_eff;

                Ok(Self {
                    external_conditions,
                    energy_supply,
                    energy_supply_connection_aux: energy_supply_conn_aux,
                    energy_supply_connections: Default::default(),
                    simulation_timestep,
                    energy_supply_type,
                    boiler_location,
                    min_modulation_load,
                    boiler_power,
                    power_circ_pump,
                    power_part_load,
                    power_full_load,
                    power_standby,
                    total_time_running_current_timestep: total_time_running_current_timestep.into(),
                    corrected_full_load_gross,
                    room_temperature,
                    temp_rise_standby_loss,
                    standby_loss_index,
                    ebv_curve_offset,
                    service_results: Default::default(),
                })
            }
            _ => unreachable!("Expected boiler data"),
        }
    }

    /// Return boiler efficiency at different return temperatures
    /// In Python this is effvsreturntemp
    fn efficiency_over_return_temperatures(
        energy_supply_type: EnergySupplyType,
        return_temp: f64,
        offset: f64,
    ) -> anyhow::Result<f64> {
        let mains_gas_dewpoint = 52.2;
        let lpg_dewpoint = 48.3;
        let theoretical_eff = match energy_supply_type {
            EnergySupplyType::MainsGas => {
                if return_temp < mains_gas_dewpoint {
                    -0.00007 * return_temp.powi(2) + 0.0017 * return_temp + 0.979
                } else {
                    -0.0006 * return_temp + 0.9129
                }
            }
            EnergySupplyType::LpgBulk
            | EnergySupplyType::LpgBottled
            | EnergySupplyType::LpgCondition11F => {
                if return_temp < lpg_dewpoint {
                    -0.00006 * return_temp.powi(2) + 0.0013 * return_temp + 0.9859
                } else {
                    -0.0006 * return_temp + 0.933
                }
            }
            _ => bail!("Unexpected energy supply type {energy_supply_type:?} encountered"),
        };

        Ok(theoretical_eff - offset)
    }

    pub fn boiler_efficiency_over_return_temperatures(
        &self,
        return_temp: f64,
        offset: f64,
    ) -> anyhow::Result<f64> {
        Self::efficiency_over_return_temperatures(self.energy_supply_type, return_temp, offset)
    }

    pub fn high_value_correction_part_load(
        energy_supply_type: EnergySupplyType,
        net_efficiency_part_load: f64,
    ) -> anyhow::Result<f64> {
        let maximum_part_load_eff = match energy_supply_type {
            EnergySupplyType::MainsGas => 1.08,
            EnergySupplyType::LpgBulk
            | EnergySupplyType::LpgBottled
            | EnergySupplyType::LpgCondition11F => 1.06,
            _ => bail!("could not calculate maximum_part_load_eff for energy supply type {energy_supply_type:?}"),
        };

        Ok(min_of_2(
            net_efficiency_part_load - 0.213 * (net_efficiency_part_load - 0.966),
            maximum_part_load_eff,
        ))
    }

    pub fn high_value_correction_full_load(net_efficiency_full_load: f64) -> f64 {
        min_of_2(
            net_efficiency_full_load - 0.673 * (net_efficiency_full_load - 0.955),
            0.98,
        )
    }

    fn net_to_gross(energy_supply_type: EnergySupplyType) -> anyhow::Result<f64> {
        match energy_supply_type {
            EnergySupplyType::MainsGas => Ok(0.901),
            EnergySupplyType::LpgBulk
            | EnergySupplyType::LpgBottled
            | EnergySupplyType::LpgCondition11F => Ok(0.921),
            _ => bail!(
                "could not convert net to gross for energy supply type '{energy_supply_type:?}'"
            ),
        }
    }

    /// Create an EnergySupplyConnection for the service name given
    pub fn create_service_connection(
        &mut self,
        service_name: Cow<'static, str>,
    ) -> anyhow::Result<()> {
        if self
            .energy_supply_connections
            .contains_key(service_name.as_ref())
        {
            bail!("Error: Service name already used: {service_name}");
        }

        self.energy_supply_connections.insert(
            service_name.to_string(),
            EnergySupply::connection(self.energy_supply.clone(), service_name.as_ref()).unwrap(),
        );

        Ok(())
    }

    pub fn create_service_hot_water_combi(
        boiler: Arc<RwLock<Self>>,
        boiler_data: HotWaterSourceDetails,
        service_name: String,
        temperature_hot_water_in_c: f64,
        cold_feed: WaterSourceWithTemperature,
    ) -> Result<BoilerServiceWaterCombi, IncorrectBoilerDataType> {
        boiler
            .write()
            .create_service_connection(service_name.clone().into())
            .unwrap();
        BoilerServiceWaterCombi::new(
            boiler.clone(),
            boiler_data,
            service_name,
            temperature_hot_water_in_c,
            cold_feed,
            boiler.read().simulation_timestep,
        )
    }

    pub fn create_service_hot_water_regular(
        boiler: Arc<RwLock<Self>>,
        service_name: String,
        temperature_hot_water_in_c: f64,
        control: Option<Arc<Control>>,
    ) -> BoilerServiceWaterRegular {
        boiler
            .write()
            .create_service_connection(service_name.clone().into())
            .unwrap();
        BoilerServiceWaterRegular::new(
            boiler.clone(),
            service_name,
            temperature_hot_water_in_c,
            control,
        )
    }

    pub fn create_service_space_heating(
        boiler: Arc<RwLock<Self>>,
        service_name: String,
        control: Arc<Control>,
    ) -> BoilerServiceSpace {
        boiler
            .write()
            .create_service_connection(service_name.clone().into())
            .unwrap();
        BoilerServiceSpace::new(boiler.clone(), service_name, control)
    }

    fn cycling_adjustment(
        &self,
        temperature_return_feed: f64,
        standing_loss: f64,
        prop_of_timestep_at_min_rate: f64,
        temperature_boiler_loc: f64,
    ) -> f64 {
        let ton_toff = (1. - prop_of_timestep_at_min_rate) / prop_of_timestep_at_min_rate;

        standing_loss
            * ton_toff
            * ((temperature_return_feed - temperature_boiler_loc) / (self.temp_rise_standby_loss))
                .powf(self.standby_loss_index)
    }

    fn location_adjustment(
        &self,
        temperature_return_feed: f64,
        standing_loss: f64,
        temperature_boiler_loc: f64,
    ) -> f64 {
        max_of_2(
            standing_loss
                * (temperature_return_feed - self.room_temperature).powf(self.standby_loss_index)
                - (temperature_return_feed - temperature_boiler_loc).powf(self.standby_loss_index),
            0.,
        )
    }

    fn calc_current_boiler_power(&self, energy_output_provided: f64, time_available: f64) -> f64 {
        if time_available <= 0. {
            return 0.0;
        }

        let min_power = self.boiler_power * self.min_modulation_load;

        max_of_2(energy_output_provided / time_available, min_power)
    }

    pub fn calc_boiler_eff(
        &self,
        service_type_is_water_combi: bool,
        temp_return_feed: f64,
        energy_output_required: f64,
        time_elapsed_hp: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let time_available = self.time_available(time_elapsed_hp);
        let energy_output_provided =
            self.calc_energy_output_provided(energy_output_required, time_available);

        let current_boiler_power =
            self.calc_current_boiler_power(energy_output_provided, time_available);

        // The efficiency of the boiler depends on whether it cycles on/off.
        // If this occurs, an adjustment is calculated for the calculation
        // timestep as follows (when the boiler is firing continuously no
        // adjustment is necessary so cycling_adjustment=0).
        let prop_of_timestep_at_min_rate = if time_available <= 0. {
            0.0
        } else {
            min_of_2(
                energy_output_required
                    / (self.boiler_power * self.min_modulation_load * time_available),
                1.0,
            )
        };

        // Default value for the stand-by heat losses as a function of the current boiler power
        // Equation 5 in EN15316-4-1
        // fgen = (c5*(Pn)^c6)/100
        // where c5 = 4.0, c6 = -0.4 and Pn is the current boiler power
        let standing_loss = if current_boiler_power == 0.0 {
            0.0
        } else {
            (4.0 * current_boiler_power.powf(-0.4)) / 100.0
        };

        // use weather temperature at timestep
        let outside_temp = self.external_conditions.air_temp(&simtime);

        // A boiler’s efficiency reduces when installed outside due to an increase in case heat loss.
        // The following adjustment is made when the boiler is located outside
        // (when installed inside no adjustment is necessary so location_adjustment=0)
        let temp_boiler_loc = match self.boiler_location {
            HeatSourceLocation::External => outside_temp,
            HeatSourceLocation::Internal => self.room_temperature,
        };

        // Calculate location adjustment
        let location_adjustment = if let HeatSourceLocation::External = self.boiler_location {
            self.location_adjustment(temp_return_feed, standing_loss, temp_boiler_loc)
        } else {
            0.0
        };

        // Calculate cycling adjustment
        let cycling_adjustment = if 0.0 < prop_of_timestep_at_min_rate
            && prop_of_timestep_at_min_rate < 1.0
            && !service_type_is_water_combi
        {
            self.cycling_adjustment(
                temp_return_feed,
                standing_loss,
                prop_of_timestep_at_min_rate,
                temp_boiler_loc,
            )
        } else {
            0.0
        };

        // Calculate combined cyclic and location adjustment
        let cyclic_location_adjustment = cycling_adjustment + location_adjustment;

        // Calculate boiler efficiency based on the return temperature and offset
        // If boiler starts cycling use the corrected full load efficiency
        // as the boiler eff before cycling adjustment is applied.
        let boiler_eff = if cycling_adjustment > 0.0 {
            self.corrected_full_load_gross
        } else {
            self.boiler_efficiency_over_return_temperatures(
                temp_return_feed,
                self.ebv_curve_offset,
            )?
        };

        let blr_eff_final = 1.0 / ((1.0 / boiler_eff) + cyclic_location_adjustment);
        Ok(blr_eff_final)
    }

    fn calc_energy_output_provided(&self, energy_output_required: f64, time_available: f64) -> f64 {
        let energy_output_max_power = self.boiler_power * time_available;

        min_of_2(energy_output_required, energy_output_max_power)
    }

    fn time_available(&self, time_elapsed_hp: Option<f64>) -> f64 {
        let timestep = self.simulation_timestep;

        match time_elapsed_hp {
            Some(time_elapsed) => timestep - time_elapsed,
            None => {
                timestep
                    - self
                        .total_time_running_current_timestep
                        .load(Ordering::SeqCst)
            }
        }
    }

    /// Calculate energy required by boiler to satisfy demand for the service indicated.
    pub fn demand_energy(
        &self,
        service_name: &str,
        service_type: ServiceType,
        energy_output_required: f64,
        temperature_return_feed: f64,
        hybrid_service_bool: Option<bool>,
        time_elapsed_hp: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, Option<f64>)> {
        let hybrid_service_bool = hybrid_service_bool.unwrap_or(false);

        // Account for time control where present. If no control present, assume
        // system is always active (except for basic thermostatic control, which
        // is implicit in demand calculation).
        let time_available = self.time_available(time_elapsed_hp);

        // If there is no demand on the boiler or no remaining time then no energy should be provided
        if energy_output_required <= 0.0 || time_available <= 0.0 {
            self.energy_supply_connections
                .get(service_name)
                .unwrap()
                .demand_energy(0.0, simtime.index)?;
        }

        let blr_eff_final = self.calc_boiler_eff(
            service_type == ServiceType::WaterCombi,
            temperature_return_feed,
            energy_output_required,
            time_elapsed_hp,
            simtime,
        )?;

        let energy_output_provided =
            self.calc_energy_output_provided(energy_output_required, time_available);
        let current_boiler_power =
            self.calc_current_boiler_power(energy_output_provided, time_available);

        let fuel_demand = energy_output_provided / blr_eff_final;

        self.energy_supply_connections
            .get(service_name)
            .unwrap()
            .demand_energy(fuel_demand, simtime.index)?;

        // Calculate running time of boiler
        let time_running_current_service = if current_boiler_power <= 0.0 {
            0.0
        } else {
            min_of_2(
                energy_output_provided / current_boiler_power,
                time_available,
            )
        };
        self.total_time_running_current_timestep
            .fetch_add(time_running_current_service, Ordering::SeqCst);

        // Save results that are needed later (in the timestep_end function)
        let mut result_service_name = ArrayString::<64>::new(); // ArrayStrings are a pain, therefore this small song and dance
        result_service_name.push_str(service_name);
        self.service_results.write().push(ServiceResult {
            _service_name: result_service_name,
            time_running: time_running_current_service,
            current_boiler_power,
        });

        Ok(if hybrid_service_bool {
            (energy_output_provided, Some(time_running_current_service))
        } else {
            (energy_output_provided, None)
        })
    }

    /// Calculation of boiler electrical consumption
    fn calc_auxiliary_energy(&mut self, time_remaining_current_timestep: f64, timestep_idx: usize) {
        // Energy used by circulation pump
        let mut energy_aux = self
            .total_time_running_current_timestep
            .load(Ordering::SeqCst)
            * self.power_circ_pump;

        // Energy used in standby mode
        energy_aux += self.power_standby * time_remaining_current_timestep;

        // Energy used by flue fan electricity for on-off boilers
        let _elec_energy_flue_fan = self
            .total_time_running_current_timestep
            .load(Ordering::SeqCst)
            * self.power_full_load;

        // Overwrite (sic from Python - no overwrite actually happens with current logic) flue fan if boiler modulates
        for service_data in self.service_results.read().iter() {
            let modulation_ratio =
                min_of_2(service_data.current_boiler_power / self.boiler_power, 1.);
            if self.min_modulation_load < 1. {
                let x_axis = [0.3, 1.0];
                let y_axis = [self.power_part_load, self.power_full_load];

                let flue_fan_el = np_interp(modulation_ratio, &x_axis, &y_axis);
                let elec_energy_flue_fan = service_data.time_running * flue_fan_el;
                energy_aux += elec_energy_flue_fan;
            }
        }

        self.energy_supply_connection_aux
            .demand_energy(energy_aux, timestep_idx)
            .unwrap();
    }

    /// Calculations to be done at the end of each timestep
    pub fn timestep_end(&mut self, simtime: SimulationTimeIteration) {
        let timestep = simtime.timestep;
        let time_remaining_current_timestep = timestep
            - self
                .total_time_running_current_timestep
                .load(Ordering::SeqCst);

        self.calc_auxiliary_energy(time_remaining_current_timestep, simtime.index);

        self.total_time_running_current_timestep = Default::default();
        self.service_results = Default::default();
    }

    pub fn energy_output_max(&self, _temp_output: f64, time_elapsed_hp: Option<f64>) -> f64 {
        let time_available = self.time_available(time_elapsed_hp);

        self.boiler_power * time_available
    }
}

#[derive(Copy, Clone, Debug)]
struct ServiceResult {
    _service_name: ArrayString<64>,
    time_running: f64,
    current_boiler_power: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::controls::time_control::SetpointTimeControl;
    use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
    use crate::external_conditions::{DaylightSavingsConfig, ShadingSegment};
    use crate::input::{ColdWaterSourceType, FuelType, HeatSourceControlType, HeatSourceWetType};
    use crate::simulation_time::SimulationTime;
    use approx::{assert_relative_eq, assert_ulps_eq};
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    pub fn boiler_data() -> HeatSourceWetDetails {
        HeatSourceWetDetails::Boiler {
            energy_supply: EnergySupplyType::MainsGas,
            energy_supply_auxiliary: EnergySupplyType::Electricity,
            rated_power: 24.0,
            efficiency_full_load: 0.88,
            efficiency_part_load: 0.986,
            boiler_location: HeatSourceLocation::Internal,
            modulation_load: 0.2,
            electricity_circ_pump: 0.0600,
            electricity_part_load: 0.0131,
            electricity_full_load: 0.0388,
            electricity_standby: 0.0244,
        }
    }

    #[fixture]
    pub fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 2., 1.)
    }

    #[fixture]
    pub fn boiler_energy_output_required() -> [f64; 2] {
        [2.0, 10.0]
    }

    #[fixture]
    pub fn temp_return_feed() -> [f64; 2] {
        [51.05, 60.00]
    }

    #[fixture]
    pub fn external_conditions(simulation_time: SimulationTime) -> ExternalConditions {
        ExternalConditions::new(
            &simulation_time.iter(),
            vec![0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0],
            vec![3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4],
            vec![200., 220., 230., 240., 250., 260., 260., 270.],
            vec![333., 610., 572., 420., 0., 10., 90., 275.],
            vec![420., 750., 425., 500., 0., 40., 0., 388.],
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
            vec![
                ShadingSegment {
                    number: 1,
                    start: 180.,
                    end: 135.,
                    shading_objects: None,
                },
                ShadingSegment {
                    number: 2,
                    start: 135.,
                    end: 90.,
                    shading_objects: None,
                },
                ShadingSegment {
                    number: 3,
                    start: 90.,
                    end: 90.,
                    shading_objects: None,
                },
                ShadingSegment {
                    number: 4,
                    start: 45.,
                    end: 0.,
                    shading_objects: None,
                },
                ShadingSegment {
                    number: 5,
                    start: 0.,
                    end: -45.,
                    shading_objects: None,
                },
                ShadingSegment {
                    number: 6,
                    start: -45.,
                    end: -90.,
                    shading_objects: None,
                },
                ShadingSegment {
                    number: 7,
                    start: -90.,
                    end: -135.,
                    shading_objects: None,
                },
                ShadingSegment {
                    number: 8,
                    start: -135.,
                    end: -180.,
                    shading_objects: None,
                },
            ],
        )
    }

    #[fixture]
    pub fn boiler(
        boiler_data: HeatSourceWetDetails,
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) -> (Boiler, Arc<RwLock<EnergySupply>>) {
        let energy_supply = Arc::new(RwLock::new(EnergySupply::new(
            FuelType::MainsGas,
            simulation_time.total_steps(),
            None,
            None,
            None,
        )));
        let energy_supply_aux = Arc::new(RwLock::new(EnergySupply::new(
            FuelType::Electricity,
            simulation_time.total_steps(),
            None,
            None,
            None,
        )));
        let energy_supply_conn_aux =
            EnergySupply::connection(energy_supply_aux, "Boiler_auxiliary").unwrap();

        let mut boiler = Boiler::new(
            boiler_data,
            energy_supply.clone(),
            energy_supply_conn_aux,
            Arc::new(external_conditions),
            simulation_time.step,
        )
        .unwrap();
        boiler
            .create_service_connection("boiler_test".into())
            .unwrap();

        (boiler, energy_supply)
    }

    #[rstest]
    pub fn should_provide_correct_energy_output(
        boiler: (Boiler, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
        boiler_energy_output_required: [f64; 2],
        temp_return_feed: [f64; 2],
    ) {
        let (boiler, energy_supply) = boiler;
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_ulps_eq!(
                boiler
                    .demand_energy(
                        "boiler_test",
                        ServiceType::WaterCombi,
                        boiler_energy_output_required[t_idx],
                        temp_return_feed[t_idx],
                        None,
                        None,
                        t_it,
                    )
                    .unwrap()
                    .0,
                [2.0, 10.0][t_idx],
            );
            assert_relative_eq!(
                energy_supply.read().results_by_end_user()["boiler_test"][t_idx],
                [2.2843673926764496, 11.5067107][t_idx],
                max_relative = 1e-7
            );
        }
    }

    #[rstest]
    pub fn should_provide_correct_efficiency_over_return_temp(
        boiler: (Boiler, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
    ) {
        let (boiler, _) = boiler;
        let return_temp = [30., 60.];
        for (idx, _) in simulation_time.iter().enumerate() {
            assert_eq!(
                boiler
                    .boiler_efficiency_over_return_temperatures(return_temp[idx], 0.)
                    .unwrap(),
                [0.967, 0.8769][idx],
                "incorrect theoretical boiler efficiency returned"
            );
        }
    }

    #[rstest]
    pub fn should_have_correct_high_value_correction(boiler: (Boiler, Arc<RwLock<EnergySupply>>)) {
        let (boiler, _) = boiler;
        assert_eq!(
            Boiler::high_value_correction_full_load(0.980),
            0.963175,
            "incorrect high value correction for full load"
        );
        assert_eq!(
            Boiler::high_value_correction_part_load(boiler.energy_supply_type, 1.081).unwrap(),
            1.056505,
            "incorrect high value correction for part load"
        );
    }

    #[rstest]
    pub fn should_calc_correct_net_to_gross(boiler: (Boiler, Arc<RwLock<EnergySupply>>)) {
        let (boiler, _) = boiler;
        assert_eq!(
            Boiler::net_to_gross(boiler.energy_supply_type).unwrap(),
            0.901,
            "incorrect net to gross"
        );
    }

    #[fixture]
    pub fn boiler_data_for_combi() -> HeatSourceWetDetails {
        HeatSourceWetDetails::Boiler {
            rated_power: 16.85,
            energy_supply: EnergySupplyType::MainsGas,
            energy_supply_auxiliary: EnergySupplyType::Electricity,
            efficiency_full_load: 0.868,
            efficiency_part_load: 0.952,
            boiler_location: HeatSourceLocation::Internal,
            modulation_load: 1.,
            electricity_circ_pump: 0.0600,
            electricity_part_load: 0.0131,
            electricity_full_load: 0.0388,
            electricity_standby: 0.0244,
        }
    }

    #[fixture]
    pub fn boiler_for_combi(
        boiler_data_for_combi: HeatSourceWetDetails,
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) -> Boiler {
        let energy_supply = Arc::new(RwLock::new(EnergySupply::new(
            FuelType::MainsGas,
            simulation_time.total_steps(),
            None,
            None,
            None,
        )));
        let energy_supply_aux = Arc::new(RwLock::new(EnergySupply::new(
            FuelType::Electricity,
            simulation_time.total_steps(),
            None,
            None,
            None,
        )));
        let energy_supply_conn_aux =
            EnergySupply::connection(energy_supply_aux, "Boiler_auxiliary").unwrap();

        let mut boiler = Boiler::new(
            boiler_data_for_combi,
            energy_supply,
            energy_supply_conn_aux,
            Arc::new(external_conditions),
            simulation_time.step,
        )
        .unwrap();
        boiler
            .create_service_connection("boiler_test".into())
            .unwrap();

        boiler
    }

    #[fixture]
    pub fn combi_boiler_data() -> HotWaterSourceDetails {
        HotWaterSourceDetails::CombiBoiler {
            separate_dhw_tests: BoilerHotWaterTest::ML,
            // fuel_energy_1: 7.099, // we don't have this field currently - unsure whether this is a mistake in the test fixture
            rejected_energy_1: 0.0004,
            // storage_loss_factor_1: 0.98328, // we don't have this field currently - unsure whether this is a mistake in the test fixture
            fuel_energy_2: 13.078,
            rejected_energy_2: 0.0004,
            storage_loss_factor_2: 0.91574,
            rejected_factor_3: 0.,
            setpoint_temp: None,
            daily_hot_water_usage: 132.5802,
            cold_water_source: ColdWaterSourceType::MainsWater,
            heat_source_wet: HeatSourceWetType::Boiler,
            control: HeatSourceControlType::HotWaterTimer,
        }
    }

    #[fixture]
    pub fn cold_water_source(simulation_time: SimulationTime) -> ColdWaterSource {
        ColdWaterSource::new(vec![1.0, 1.2], 0, simulation_time.step)
    }

    #[fixture]
    pub fn combi_boiler(
        boiler_for_combi: Boiler,
        combi_boiler_data: HotWaterSourceDetails,
        cold_water_source: ColdWaterSource,
        simulation_time: SimulationTime,
    ) -> BoilerServiceWaterCombi {
        BoilerServiceWaterCombi::new(
            Arc::new(RwLock::new(boiler_for_combi)),
            combi_boiler_data,
            "boiler_test".to_string(),
            60.,
            WaterSourceWithTemperature::ColdWaterSource(Arc::new(cold_water_source)),
            simulation_time.step,
        )
        .unwrap()
    }

    #[fixture]
    pub fn volume_demanded() -> [IndexMap<DemandVolTargetKey, VolumeReference>; 2] {
        [
            IndexMap::from([
                (
                    41.0.into(),
                    VolumeReference {
                        warm_temp: 41.0,
                        warm_vol: 48.0,
                    },
                ),
                (
                    43.0.into(),
                    VolumeReference {
                        warm_temp: 43.0,
                        warm_vol: 100.0,
                    },
                ),
                (
                    40.0.into(),
                    VolumeReference {
                        warm_temp: 40.0,
                        warm_vol: 0.0,
                    },
                ),
                (
                    DemandVolTargetKey::TempHotWater,
                    VolumeReference {
                        warm_temp: 55.0,
                        warm_vol: 110.59194954841298,
                    },
                ),
            ]),
            IndexMap::from([
                (
                    41.0.into(),
                    VolumeReference {
                        warm_temp: 41.0,
                        warm_vol: 48.0,
                    },
                ),
                (
                    DemandVolTargetKey::TempHotWater,
                    VolumeReference {
                        warm_temp: 55.0,
                        warm_vol: 32.60190808710678,
                    },
                ),
            ]),
        ]
    }

    #[rstest]
    pub fn combi_boiler_should_provide_demand_hot_water(
        combi_boiler: BoilerServiceWaterCombi,
        simulation_time: SimulationTime,
        volume_demanded: [IndexMap<DemandVolTargetKey, VolumeReference>; 2],
    ) {
        for (idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                combi_boiler
                    .demand_hot_water(volume_demanded[idx].clone(), t_it)
                    .unwrap(),
                [7.624602058956146, 2.267017951167212][idx],
                max_relative = 1e-6
            );
        }
    }

    #[fixture]
    pub fn boiler_data_for_regular() -> HeatSourceWetDetails {
        HeatSourceWetDetails::Boiler {
            rated_power: 24.0,
            energy_supply: EnergySupplyType::MainsGas,
            energy_supply_auxiliary: EnergySupplyType::Electricity,
            efficiency_full_load: 0.891,
            efficiency_part_load: 0.991,
            boiler_location: HeatSourceLocation::Internal,
            modulation_load: 0.3,
            electricity_circ_pump: 0.0600,
            electricity_part_load: 0.0131,
            electricity_full_load: 0.0388,
            electricity_standby: 0.0244,
        }
    }

    #[fixture]
    pub fn boiler_for_regular(
        boiler_data_for_regular: HeatSourceWetDetails,
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) -> Boiler {
        let energy_supply = Arc::new(RwLock::new(EnergySupply::new(
            FuelType::MainsGas,
            simulation_time.total_steps(),
            None,
            None,
            None,
        )));
        let energy_supply_aux = Arc::new(RwLock::new(EnergySupply::new(
            FuelType::Electricity,
            simulation_time.total_steps(),
            None,
            None,
            None,
        )));
        let energy_supply_conn_aux =
            EnergySupply::connection(energy_supply_aux, "Boiler_auxiliary").unwrap();

        let mut boiler = Boiler::new(
            boiler_data_for_regular,
            energy_supply,
            energy_supply_conn_aux,
            Arc::new(external_conditions),
            simulation_time.step,
        )
        .unwrap();
        boiler
            .create_service_connection("boiler_test".into())
            .unwrap();

        boiler
    }

    #[fixture]
    pub fn regular_boiler<'a>(boiler_for_regular: Boiler) -> BoilerServiceWaterRegular {
        BoilerServiceWaterRegular::new(
            Arc::new(RwLock::new(boiler_for_regular)),
            "boiler_test".to_string(),
            60.,
            None,
        )
    }

    #[rstest]
    pub fn regular_boiler_should_provide_demand_hot_water(
        mut regular_boiler: BoilerServiceWaterRegular,
        simulation_time: SimulationTime,
    ) {
        let temp_return_feed = [51.05, 60.00];
        for (idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                regular_boiler
                    .demand_energy(
                        [0.7241412, 0.1748878][idx],
                        temp_return_feed[idx],
                        None,
                        None,
                        t_it
                    )
                    .unwrap()
                    .0,
                [0.7241412, 0.1748878][idx],
                max_relative = 1e-7
            );
        }
    }

    #[fixture]
    pub fn boiler_data_for_service_space() -> HeatSourceWetDetails {
        HeatSourceWetDetails::Boiler {
            rated_power: 16.85,
            energy_supply: EnergySupplyType::MainsGas,
            energy_supply_auxiliary: EnergySupplyType::Electricity,
            efficiency_full_load: 0.868,
            efficiency_part_load: 0.952,
            boiler_location: HeatSourceLocation::Internal,
            modulation_load: 1.0,
            electricity_circ_pump: 0.0600,
            electricity_part_load: 0.0131,
            electricity_full_load: 0.0388,
            electricity_standby: 0.0244,
        }
    }

    #[fixture]
    pub fn simulation_time_for_service_space() -> SimulationTime {
        SimulationTime::new(0., 3., 1.)
    }

    #[fixture]
    pub fn boiler_for_service_space(
        boiler_data_for_service_space: HeatSourceWetDetails,
        external_conditions: ExternalConditions,
        simulation_time_for_service_space: SimulationTime,
    ) -> Boiler {
        let energy_supply = Arc::new(RwLock::new(EnergySupply::new(
            FuelType::MainsGas,
            simulation_time_for_service_space.total_steps(),
            None,
            None,
            None,
        )));
        let energy_supply_aux = Arc::new(RwLock::new(EnergySupply::new(
            FuelType::Electricity,
            simulation_time_for_service_space.total_steps(),
            None,
            None,
            None,
        )));
        let energy_supply_conn_aux =
            EnergySupply::connection(energy_supply_aux, "Boiler_auxiliary").unwrap();

        let mut boiler = Boiler::new(
            boiler_data_for_service_space,
            energy_supply,
            energy_supply_conn_aux,
            Arc::new(external_conditions),
            simulation_time_for_service_space.step,
        )
        .unwrap();
        boiler
            .create_service_connection("boiler_test".into())
            .unwrap();

        boiler
    }

    #[fixture]
    pub fn control_for_service_space() -> Control {
        Control::SetpointTimeControl(
            SetpointTimeControl::new(
                vec![Some(21.0), Some(21.0), None],
                0,
                1.0,
                None,
                None,
                None,
                None,
                1.0,
            )
            .unwrap(),
        )
    }

    #[fixture]
    pub fn service_space_boiler(
        boiler_for_service_space: Boiler,
        control_for_service_space: Control,
    ) -> BoilerServiceSpace {
        BoilerServiceSpace::new(
            Arc::new(RwLock::new(boiler_for_service_space)),
            "boiler_test".to_string(),
            Arc::new(control_for_service_space),
        )
    }

    #[rstest]
    pub fn service_space_boiler_should_provide_demand_hot_water(
        mut service_space_boiler: BoilerServiceSpace,
        simulation_time_for_service_space: SimulationTime,
    ) {
        let energy_demanded = [10.0, 2.0, 2.0];
        let temp_flow = [55.0, 65.0, 65.0];
        let temp_return_feed = [50.0, 60.0, 60.0];
        for (idx, t_it) in simulation_time_for_service_space.iter().enumerate() {
            assert_ulps_eq!(
                service_space_boiler
                    .demand_energy(
                        energy_demanded[idx],
                        temp_flow[idx],
                        temp_return_feed[idx],
                        None,
                        None,
                        t_it,
                    )
                    .unwrap()
                    .0,
                [10.0, 2.0, 0.0][idx]
            );
        }
    }
}
