use crate::compare_floats::{max_of_2, min_of_2};
use crate::core::common::WaterSourceWithTemperature;
use crate::core::controls::time_control::{Control, ControlBehaviour};
use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::core::material_properties::WATER;
use crate::core::units::{DAYS_PER_YEAR, HOURS_PER_DAY, WATTS_PER_KILOWATT};
use crate::core::water_heat_demand::dhw_demand::{DemandVolTargetKey, VolumeReference};
use crate::external_conditions::ExternalConditions;
use crate::input::{BoilerHotWaterTest, FuelType, HotWaterSourceDetails};
use crate::input::{HeatSourceLocation, HeatSourceWetDetails};
use crate::simulation_time::SimulationTimeIteration;
use crate::statistics::np_interp;
use anyhow::bail;
use atomic_float::AtomicF64;
use indexmap::IndexMap;
use parking_lot::RwLock;
use smartstring::alias::String;
use std::collections::HashMap;
use std::fmt;
use std::sync::atomic::Ordering;
use std::sync::Arc;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ServiceType {
    WaterCombi,
    WaterRegular,
    Space,
}

#[derive(Debug)]
pub(crate) struct BoilerServiceWaterCombi {
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
                daily_hw_usage: daily_hot_water_usage,
                ..
            } => {
                let (rejected_energy_1, storage_loss_factor_2, rejected_factor_3) =
                    match separate_dhw_tests {
                        BoilerHotWaterTest::ML | BoilerHotWaterTest::MS => {
                            (rejected_energy_1, storage_loss_factor_2, rejected_factor_3)
                        }
                        BoilerHotWaterTest::MOnly => {
                            (rejected_energy_1, storage_loss_factor_2, None)
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
            self.cold_feed.temperature(simtime, None),
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
                Some(return_temperature),
                None,
                None,
                None,
                None,
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
        self.boiler.read().energy_output_max(None, None)
    }
}

#[derive(Debug)]
pub struct BoilerServiceWaterRegular {
    boiler: Arc<RwLock<Boiler>>,
    service_name: String,
    control_min: Arc<Control>,
    _control_max: Arc<Control>,
}

impl BoilerServiceWaterRegular {
    pub(crate) fn new(
        boiler: Arc<RwLock<Boiler>>,
        service_name: String,
        control_min: Arc<Control>,
        control_max: Arc<Control>,
    ) -> anyhow::Result<Self> {
        Ok(Self {
            boiler,
            service_name,
            control_min,
            _control_max: control_max.clone(),
        })
    }

    /// Return setpoint (not necessarily temperature)
    pub(crate) fn setpnt(&self, simtime: SimulationTimeIteration) -> (Option<f64>, Option<f64>) {
        (
            self.control_min.setpnt(&simtime),
            self._control_max.setpnt(&simtime),
        )
    }

    pub fn demand_energy(
        &self,
        mut energy_demand: f64,
        _temp_flow: f64,
        temp_return: Option<f64>,
        hybrid_service_bool: Option<bool>,
        time_elapsed_hp: Option<f64>,
        update_heat_source_state: Option<bool>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, Option<f64>)> {
        let hybrid_service_bool = hybrid_service_bool.unwrap_or(false);

        if !self.is_on(simtime) {
            energy_demand = 0.;
        }

        if temp_return.is_none() && energy_demand != 0. {
            bail!("temp_return is None and energy_demand is not 0.0");
        }

        self.boiler.write().demand_energy(
            &self.service_name,
            ServiceType::WaterRegular,
            energy_demand,
            temp_return,
            None,
            Some(hybrid_service_bool),
            time_elapsed_hp,
            Some(update_heat_source_state.unwrap_or(true)),
        )
    }

    pub fn energy_output_max(
        &self,
        _temp_flow: f64,
        _temp_return: f64,
        time_elapsed_hp: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        if !self.is_on(simtime) {
            return 0.;
        }

        self.boiler.read().energy_output_max(None, time_elapsed_hp)
    }

    fn is_on(&self, simtime: SimulationTimeIteration) -> bool {
        self.control_min.is_on(simtime)
    }
}

/// A struct representing a space heating service provided by a boiler to e.g. a cylinder.
#[derive(Clone, Debug)]
pub struct BoilerServiceSpace {
    boiler: Arc<RwLock<Boiler>>,
    service_name: String,
    control: Arc<Control>,
}

impl BoilerServiceSpace {
    pub(crate) fn new(
        boiler: Arc<RwLock<Boiler>>,
        service_name: String,
        control: Arc<Control>,
    ) -> Self {
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
        &self,
        energy_demand: f64,
        _temp_flow: f64,
        temp_return: Option<f64>,
        time_start: Option<f64>,
        hybrid_service_bool: Option<bool>,
        time_elapsed_hp: Option<f64>,
        update_heat_source_state: Option<bool>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, Option<f64>)> {
        let time_start = time_start.unwrap_or(0.);
        let hybrid_service_bool = hybrid_service_bool.unwrap_or(false);
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);

        let energy_demand = if !self.is_on(simtime) {
            0.0
        } else {
            energy_demand
        };

        self.boiler.write().demand_energy(
            &self.service_name,
            ServiceType::Space,
            energy_demand,
            temp_return,
            Some(time_start),
            Some(hybrid_service_bool),
            time_elapsed_hp,
            Some(update_heat_source_state),
        )
    }

    pub fn energy_output_max(
        &self,
        _temp_output: f64,
        _temp_return_feed: f64,
        time_start: Option<f64>,
        time_elapsed_hp: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> f64 {
        if !self.is_on(simtime) {
            0.0
        } else {
            self.boiler
                .read()
                .energy_output_max(time_start, time_elapsed_hp)
        }
    }

    fn is_on(&self, simtime: SimulationTimeIteration) -> bool {
        self.control.is_on(simtime)
    }
}

#[derive(Debug)]
pub(crate) struct Boiler {
    energy_supply: Arc<RwLock<EnergySupply>>,
    simulation_timestep: f64,
    external_conditions: Arc<ExternalConditions>,
    energy_supply_connections: HashMap<String, EnergySupplyConnection>,
    energy_supply_connection_aux: EnergySupplyConnection,
    _energy_supply_type: String,
    // service_results: (),
    boiler_location: HeatSourceLocation,
    min_modulation_load: f64,
    boiler_power: f64,
    fuel_code: FuelType,
    power_circ_pump: f64,
    power_part_load: f64,
    power_full_load: f64,
    power_standby: f64,
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
    pub(crate) fn new(
        boiler_data: HeatSourceWetDetails,
        energy_supply: Arc<RwLock<EnergySupply>>,
        energy_supply_conn_aux: EnergySupplyConnection,
        external_conditions: Arc<ExternalConditions>,
        simulation_timestep: f64,
    ) -> anyhow::Result<Self> {
        match boiler_data {
            HeatSourceWetDetails::Boiler {
                energy_supply: energy_supply_type,
                energy_supply_aux: _,
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

                let fuel_code = energy_supply.read().fuel_type();

                let net_to_gross = Self::net_to_gross(&fuel_code)?;
                let full_load_net = full_load_gross / net_to_gross;
                let part_load_net = part_load_gross / net_to_gross;
                let corrected_full_load_net = Self::high_value_correction_full_load(full_load_net);
                let corrected_part_load_net =
                    Self::high_value_correction_part_load(&fuel_code, part_load_net)?;
                let corrected_full_load_gross = corrected_full_load_net * net_to_gross;
                let corrected_part_load_gross = corrected_part_load_net * net_to_gross;

                // SAP model properties
                let room_temperature = 19.5; // TODO (from Python) use actual room temp instead of hard coding

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
                    &fuel_code,
                    temp_part_load_test,
                    offset_for_theoretical_eff,
                )?;
                let theoretical_eff_full_load = Self::efficiency_over_return_temperatures(
                    &fuel_code,
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
                    _energy_supply_type: energy_supply_type,
                    boiler_location,
                    min_modulation_load,
                    boiler_power,
                    fuel_code,
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
        fuel_code: &FuelType,
        return_temp: f64,
        offset: f64,
    ) -> anyhow::Result<f64> {
        let mains_gas_dewpoint = 52.2;
        let lpg_dewpoint = 48.3;
        let theoretical_eff = match fuel_code {
            FuelType::MainsGas => {
                if return_temp < mains_gas_dewpoint {
                    -0.00007 * return_temp.powi(2) + 0.0017 * return_temp + 0.979
                } else {
                    -0.0006 * return_temp + 0.9129
                }
            }
            FuelType::LpgBulk | FuelType::LpgBottled | FuelType::LpgCondition11F => {
                if return_temp < lpg_dewpoint {
                    -0.00006 * return_temp.powi(2) + 0.0013 * return_temp + 0.9859
                } else {
                    -0.0006 * return_temp + 0.933
                }
            }
            _ => bail!("Unexpected fuel code {fuel_code:?} encountered"),
        };

        Ok(theoretical_eff - offset)
    }

    pub fn boiler_efficiency_over_return_temperatures(
        &self,
        return_temp: f64,
        offset: f64,
    ) -> anyhow::Result<f64> {
        Self::efficiency_over_return_temperatures(&self.fuel_code, return_temp, offset)
    }

    fn high_value_correction_part_load(
        fuel_code: &FuelType,
        net_efficiency_part_load: f64,
    ) -> anyhow::Result<f64> {
        let maximum_part_load_eff = match fuel_code {
            FuelType::MainsGas => 1.08,
            FuelType::LpgBulk | FuelType::LpgBottled | FuelType::LpgCondition11F => 1.06,
            _ => bail!("could not calculate maximum_part_load_eff for fuel_code {fuel_code:?}"),
        };

        Ok(min_of_2(
            net_efficiency_part_load - 0.213 * (net_efficiency_part_load - 0.966),
            maximum_part_load_eff,
        ))
    }

    fn high_value_correction_full_load(net_efficiency_full_load: f64) -> f64 {
        min_of_2(
            net_efficiency_full_load - 0.673 * (net_efficiency_full_load - 0.955),
            0.98,
        )
    }

    fn net_to_gross(fuel_code: &FuelType) -> anyhow::Result<f64> {
        match fuel_code {
            FuelType::MainsGas => Ok(0.901),
            FuelType::LpgBulk | FuelType::LpgBottled | FuelType::LpgCondition11F => Ok(0.921),
            _ => bail!("could not convert net to gross for fuel code '{fuel_code:?}'"),
        }
    }

    /// Create an EnergySupplyConnection for the service name given
    pub fn create_service_connection(&mut self, service_name: &str) -> anyhow::Result<()> {
        if self.energy_supply_connections.contains_key(service_name) {
            bail!("Error: Service name already used: {service_name}");
        }

        self.energy_supply_connections.insert(
            service_name.into(),
            EnergySupply::connection(self.energy_supply.clone(), service_name).unwrap(),
        );

        Ok(())
    }

    pub(crate) fn create_service_hot_water_combi(
        boiler: Arc<RwLock<Self>>,
        boiler_data: HotWaterSourceDetails,
        service_name: &str,
        temperature_hot_water_in_c: f64,
        cold_feed: WaterSourceWithTemperature,
    ) -> Result<BoilerServiceWaterCombi, IncorrectBoilerDataType> {
        boiler
            .write()
            .create_service_connection(service_name)
            .unwrap();
        BoilerServiceWaterCombi::new(
            boiler.clone(),
            boiler_data,
            service_name.into(),
            temperature_hot_water_in_c,
            cold_feed,
            boiler.read().simulation_timestep,
        )
    }

    pub(crate) fn create_service_hot_water_regular(
        boiler: Arc<RwLock<Self>>,
        service_name: &str,
        control_min: Arc<Control>,
        control_max: Arc<Control>,
    ) -> anyhow::Result<BoilerServiceWaterRegular> {
        boiler.write().create_service_connection(service_name)?;
        BoilerServiceWaterRegular::new(
            boiler.clone(),
            service_name.into(),
            control_min,
            control_max,
        )
    }

    pub(crate) fn create_service_space_heating(
        boiler: Arc<RwLock<Self>>,
        service_name: &str,
        control: Arc<Control>,
    ) -> BoilerServiceSpace {
        boiler
            .write()
            .create_service_connection(service_name)
            .unwrap();
        BoilerServiceSpace::new(boiler.clone(), service_name.into(), control)
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

    pub(crate) fn calc_boiler_eff(
        &self,
        service_type_is_water_combi: bool,
        temp_return_feed: f64,
        energy_output_required: f64,
        time_start: Option<f64>,
        time_elapsed_hp: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let time_start = time_start.unwrap_or(0.0);
        let time_available = self.time_available(time_start, time_elapsed_hp);

        self.calc_boiler_eff_internal(
            service_type_is_water_combi,
            temp_return_feed,
            energy_output_required,
            time_available,
            simtime,
        )
    }

    fn calc_boiler_eff_internal(
        &self,
        service_type_is_water_combi: bool,
        temp_return_feed: f64,
        energy_output_required: f64,
        time_available: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
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

    /// Calculate time available for the current service
    fn time_available(&self, time_start: f64, time_elapsed_hp: Option<f64>) -> f64 {
        // Assumes that time spent on other services is evenly spread throughout
        // the timestep so the adjustment for start time below is a proportional
        // reduction of the overall time available, not simply a subtraction
        //         timestep = self.__simulation_time.timestep()
        //         total_time_running_current_timestep \
        //             = time_elapsed_hp if time_elapsed_hp is not None else self.__total_time_running_current_timestep
        //         time_available \
        //             = (timestep - total_time_running_current_timestep) * (1.0 - time_start / timestep)
        //         return time_available
        let timestep = self.simulation_timestep;
        let total_time_running_current_timestep = if let Some(time_elapsed_hp) = time_elapsed_hp {
            time_elapsed_hp
        } else {
            self.total_time_running_current_timestep
                .load(Ordering::SeqCst)
        };
        (timestep - total_time_running_current_timestep) * (1.0 - time_start / timestep)
    }

    /// Calculate running time of Boiler
    fn time_running(&self, energy_output_provided: f64, time_available: f64) -> f64 {
        let current_boiler_power =
            self.calc_current_boiler_power(energy_output_provided, time_available);
        if current_boiler_power <= 0.0 {
            0.0
        } else {
            time_available.min(energy_output_provided / current_boiler_power)
        }
    }

    /// Calculate energy required by boiler to satisfy demand for the service indicated.
    pub(crate) fn demand_energy(
        &self,
        service_name: &str,
        service_type: ServiceType,
        energy_output_required: f64,
        temp_return_feed: Option<f64>,
        time_start: Option<f64>,
        hybrid_service_bool: Option<bool>,
        time_elapsed_hp: Option<f64>,
        update_heat_source_state: Option<bool>,
    ) -> anyhow::Result<(f64, Option<f64>)> {
        let time_start = time_start.unwrap_or(0.0);
        let hybrid_service_bool = hybrid_service_bool.unwrap_or(false);
        let update_heat_source_state = update_heat_source_state.unwrap_or(true);

        // Account for time control where present. If no control present, assume
        // system is always active (except for basic thermostatic control, which
        // is implicit in demand calculation).
        let time_available = self.time_available(time_start, time_elapsed_hp);

        let energy_output_provided =
            self.calc_energy_output_provided(energy_output_required, time_available);

        // TODO (from Python) Ideally, the boiler power used for the running time calculation
        //      would account for space heating demand for all zones, but the
        //      calculation flow does not allow for this without circularity.
        //      Therefore, the value for time running returned from this function
        //      (used in the hybrid HP calculation) will be slightly inaccurate.
        let time_running_current_service =
            self.time_running(energy_output_provided, time_available);

        if update_heat_source_state {
            self.total_time_running_current_timestep
                .fetch_add(time_running_current_service, Ordering::SeqCst);

            // Save results that are needed later (in the timestep_end function)
            self.service_results.write().push(ServiceResult {
                service_name: service_name.into(),
                service_type,
                temp_return_feed,
                energy_output_required,
                energy_output_provided,
                time_available,
                _time_start: time_start,
                _time_elapsed_hp: time_elapsed_hp,
            });
        }

        Ok(if hybrid_service_bool {
            (energy_output_provided, Some(time_running_current_service))
        } else {
            (energy_output_provided, None)
        })
    }

    fn sum_space_heating_service_results_energy_output_required(&self) -> f64 {
        self.service_results
            .read()
            .iter()
            .map(|r| r.energy_output_required)
            .sum()
    }

    fn sum_space_heating_service_results_energy_output_provided(&self) -> f64 {
        self.service_results
            .read()
            .iter()
            .map(|r| r.energy_output_provided)
            .sum()
    }

    /// Calculate boiler fuel demand for all services (excl. auxiliary),
    /// and request this from relevant EnergySupplyConnection
    fn fuel_demand(&self, simtime: SimulationTimeIteration) -> anyhow::Result<()> {
        // pre-calc max time available outside the loop
        // (Get max time available from all space heating services to use
        // as overall time available for all space heating services. Note
        // that for this assumption to be valid, the space heating
        // services must be called consecutively, with no services of
        // another type called in between.)
        let max_time_available = self
            .service_results
            .read()
            .iter()
            .filter_map(|x| (x.service_type == ServiceType::Space).then_some(x.time_available))
            .sum::<f64>();

        for service_data in self.service_results.read().iter() {
            let service_name = service_data.service_name.as_str();
            let service_type = service_data.service_type;
            let temp_return_feed = service_data.temp_return_feed;
            let energy_output_provided = service_data.energy_output_provided;

            // Aggregate space heating services
            // TODO (from Python) This is only necessary because the model cannot handle an
            //                    emitter circuit that serves more than one zone. If/when this
            //                    capability is added, there will no longer be separate space
            //                    heating services for each zone and this aggregation can be
            //                    removed as it will not be necessary. At that point, the other
            //                    contents of this function could also be moved back to their
            //                    original locations
            let (combined_energy_output_required, time_available) =
                if service_type == ServiceType::Space {
                    (
                        self.sum_space_heating_service_results_energy_output_required(),
                        max_time_available,
                    )
                } else {
                    (
                        service_data.energy_output_required,
                        service_data.time_available,
                    )
                };

            let fuel_demand = if let Some(temp_return_feed) = temp_return_feed {
                let blr_eff_final = self.calc_boiler_eff_internal(
                    service_type == ServiceType::WaterCombi,
                    temp_return_feed,
                    combined_energy_output_required,
                    time_available,
                    simtime,
                )?;
                energy_output_provided / blr_eff_final
            } else {
                0.0
            };
            self.energy_supply_connections[service_name]
                .demand_energy(fuel_demand, simtime.index)?;
        }

        Ok(())
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

        // Overwrite flue fan electricity if boiler modulates
        let mut space_heat_services_processed = false;

        // pre-calc max time available outside the loop
        // (Get max time available from all space heating services to use
        // as overall time available for all space heating services. Note
        // that for this assumption to be valid, the space heating
        // services must be called consecutively, with no services of
        // another type called in between.)
        let max_time_available = self
            .service_results
            .read()
            .iter()
            .filter_map(|x| (x.service_type == ServiceType::Space).then_some(x.time_available))
            .sum::<f64>();

        for service_data in self.service_results.read().iter() {
            // Aggregate space heating services
            // TODO (from Python) This is only necessary because the model cannot handle an
            //                    emitter circuit that serves more than one zone. If/when this
            //                    capability is added, there will no longer be separate space
            //                    heating services for each zone and this aggregation can be
            //                    removed as it will not be necessary. At that point, the other
            //                    contents of this function could also be moved back to their
            //                    original locations
            let (combined_energy_output_required, time_available) =
                if service_data.service_type == ServiceType::Space {
                    if space_heat_services_processed {
                        continue;
                    }

                    space_heat_services_processed = true;

                    (
                        self.sum_space_heating_service_results_energy_output_provided(),
                        max_time_available,
                    )
                } else {
                    (
                        service_data.energy_output_provided,
                        service_data.time_available,
                    )
                };

            let current_boiler_power =
                self.calc_current_boiler_power(combined_energy_output_required, time_available);
            let modulation_ratio = (current_boiler_power / self.boiler_power).min(1.0);
            if self.min_modulation_load < 1. {
                let x_axis = [self.min_modulation_load, 1.];
                let y_axis = [self.power_part_load, self.power_full_load];

                let flue_fan_el = np_interp(modulation_ratio, &x_axis, &y_axis);
                let time_running =
                    self.time_running(combined_energy_output_required, time_available);
                let elec_energy_flue_fan = time_running * flue_fan_el;
                energy_aux += elec_energy_flue_fan;
            }
        }

        self.energy_supply_connection_aux
            .demand_energy(energy_aux, timestep_idx)
            .unwrap();
    }

    /// Calculations to be done at the end of each timestep
    pub(crate) fn timestep_end(&mut self, simtime: SimulationTimeIteration) -> anyhow::Result<()> {
        self.fuel_demand(simtime)?;

        let timestep = simtime.timestep;
        let time_remaining_current_timestep = timestep
            - self
                .total_time_running_current_timestep
                .load(Ordering::SeqCst);

        self.calc_auxiliary_energy(time_remaining_current_timestep, simtime.index);

        self.total_time_running_current_timestep = Default::default();
        self.service_results = Default::default();

        Ok(())
    }

    pub fn energy_output_max(&self, time_start: Option<f64>, time_elapsed_hp: Option<f64>) -> f64 {
        let time_start = time_start.unwrap_or(0.0);
        let time_available = self.time_available(time_start, time_elapsed_hp);

        self.boiler_power * time_available
    }
}

#[derive(Clone, Debug)]
struct ServiceResult {
    service_name: String,
    service_type: ServiceType,
    temp_return_feed: Option<f64>,
    energy_output_required: f64,
    energy_output_provided: f64,
    time_available: f64,
    _time_start: f64,
    _time_elapsed_hp: Option<f64>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::controls::time_control::SetpointTimeControl;
    use crate::core::energy_supply::energy_supply::EnergySupplyBuilder;
    use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
    use crate::external_conditions::{DaylightSavingsConfig, ShadingSegment};
    use crate::input::FuelType;
    use crate::simulation_time::SimulationTime;
    use approx::{assert_relative_eq, assert_ulps_eq};
    use itertools::Itertools;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use serde_json::json;
    use std::any::type_name;

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 2., 1.)
    }

    #[fixture]
    fn external_conditions(simulation_time: SimulationTime) -> ExternalConditions {
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
                    start: 180.,
                    end: 135.,
                    ..Default::default()
                },
                ShadingSegment {
                    start: 135.,
                    end: 90.,
                    ..Default::default()
                },
                ShadingSegment {
                    start: 90.,
                    end: 90.,
                    ..Default::default()
                },
                ShadingSegment {
                    start: 45.,
                    end: 0.,
                    ..Default::default()
                },
                ShadingSegment {
                    start: 0.,
                    end: -45.,
                    ..Default::default()
                },
                ShadingSegment {
                    start: -45.,
                    end: -90.,
                    ..Default::default()
                },
                ShadingSegment {
                    start: -90.,
                    end: -135.,
                    ..Default::default()
                },
                ShadingSegment {
                    start: -135.,
                    end: -180.,
                    ..Default::default()
                },
            ]
            .into(),
        )
    }

    mod test_boiler_service_water_combi {
        use crate::core::common::WaterSourceWithTemperature;
        use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyBuilder};
        use crate::core::heating_systems::boiler::tests::{external_conditions, simulation_time};
        use crate::core::heating_systems::boiler::{Boiler, BoilerServiceWaterCombi};
        use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
        use crate::core::water_heat_demand::dhw_demand::{DemandVolTargetKey, VolumeReference};
        use crate::hem_core::external_conditions::ExternalConditions;
        use crate::hem_core::simulation_time::SimulationTime;
        use crate::input::{
            BoilerHotWaterTest, FuelType, HeatSourceLocation, HeatSourceWetDetails,
            HotWaterSourceDetails,
        };
        use approx::assert_relative_eq;
        use indexmap::IndexMap;
        use parking_lot::RwLock;
        use rstest::{fixture, rstest};
        use std::sync::Arc;

        #[fixture]
        pub fn boiler_data() -> HeatSourceWetDetails {
            HeatSourceWetDetails::Boiler {
                rated_power: 16.85,
                energy_supply: "mains gas".into(),
                energy_supply_aux: "mains elec".into(),
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
        fn boiler_service_water_combi_data() -> HotWaterSourceDetails {
            HotWaterSourceDetails::CombiBoiler {
                separate_dhw_tests: BoilerHotWaterTest::ML,
                // fuel_energy_1: 7.099, // we don't have this field currently - unsure whether this is a mistake in the test fixture
                rejected_energy_1: Some(0.0004),
                // storage_loss_factor_1: 0.98328, // we don't have this field currently - unsure whether this is a mistake in the test fixture
                storage_loss_factor_2: Some(0.91574),
                rejected_factor_3: Some(0.),
                setpoint_temp: None,
                daily_hw_usage: 132.5802,
                cold_water_source: "mains water".into(),
                heat_source_wet: "boiler".into(),
            }
        }

        #[fixture]
        fn volume_demanded() -> [IndexMap<DemandVolTargetKey, VolumeReference>; 2] {
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

        #[fixture]
        pub fn boiler(
            boiler_data: HeatSourceWetDetails,
            external_conditions: ExternalConditions,
            simulation_time: SimulationTime,
        ) -> Boiler {
            let energy_supply = Arc::new(RwLock::new(
                EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time.total_steps()).build(),
            ));
            let energy_supply_aux = Arc::new(RwLock::new(
                EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps())
                    .build(),
            ));
            let energy_supply_conn_aux =
                EnergySupply::connection(energy_supply_aux, "Boiler_auxiliary").unwrap();

            let mut boiler = Boiler::new(
                boiler_data,
                energy_supply,
                energy_supply_conn_aux,
                Arc::new(external_conditions),
                simulation_time.step,
            )
            .unwrap();
            boiler.create_service_connection("boiler_test").unwrap();

            boiler
        }

        #[fixture]
        fn boiler_service_water(
            boiler: Boiler,
            boiler_service_water_combi_data: HotWaterSourceDetails,
            simulation_time: SimulationTime,
        ) -> BoilerServiceWaterCombi {
            let cold_water_source = ColdWaterSource::new(vec![1.0, 1.2], 0, simulation_time.step);

            BoilerServiceWaterCombi::new(
                Arc::new(RwLock::new(boiler)),
                boiler_service_water_combi_data,
                "boiler_test".into(),
                60.,
                WaterSourceWithTemperature::ColdWaterSource(Arc::new(cold_water_source)),
                simulation_time.step,
            )
            .unwrap()
        }

        #[rstest]
        fn test_boiler_service_water(
            boiler_service_water: BoilerServiceWaterCombi,
            simulation_time: SimulationTime,
            volume_demanded: [IndexMap<DemandVolTargetKey, VolumeReference>; 2],
        ) {
            for (idx, t_it) in simulation_time.iter().enumerate() {
                assert_relative_eq!(
                    boiler_service_water
                        .demand_hot_water(volume_demanded[idx].clone(), t_it)
                        .unwrap(),
                    [7.624602058956146, 2.267017951167212][idx],
                    max_relative = 1e-6
                );
            }
        }
    }

    #[fixture]
    pub fn boiler_data() -> HeatSourceWetDetails {
        HeatSourceWetDetails::Boiler {
            energy_supply: "mains gas".into(),
            energy_supply_aux: "mains elec".into(),
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
    pub fn boiler_energy_output_required() -> [f64; 2] {
        [2.0, 10.0]
    }

    #[fixture]
    pub fn temp_return_feed() -> [f64; 2] {
        [51.05, 60.00]
    }

    #[fixture]
    fn boiler(
        boiler_data: HeatSourceWetDetails,
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) -> (Boiler, Arc<RwLock<EnergySupply>>) {
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time.total_steps()).build(),
        ));
        let energy_supply_aux = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
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
        boiler.create_service_connection("boiler_test").unwrap();

        (boiler, energy_supply)
    }

    #[fixture]
    fn boiler_data_for_regular() -> HeatSourceWetDetails {
        HeatSourceWetDetails::Boiler {
            rated_power: 24.0,
            energy_supply: "mains gas".into(),
            energy_supply_aux: "mains elec".into(),
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
    fn boiler_for_regular(
        boiler_data_for_regular: HeatSourceWetDetails,
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) -> Boiler {
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time.total_steps()).build(),
        ));
        let energy_supply_aux = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
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
        boiler.create_service_connection("boiler_test").unwrap();

        boiler
    }

    #[fixture]
    fn control_min() -> Arc<Control> {
        Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(52.), Some(52.)],
            0,
            1.,
            Default::default(),
            Default::default(),
            1.,
        )))
    }

    #[fixture]
    fn control_max() -> Arc<Control> {
        Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(60.), Some(60.)],
            0,
            1.,
            Default::default(),
            Default::default(),
            1.,
        )))
    }

    #[fixture]
    fn regular_boiler<'a>(
        boiler_for_regular: Boiler,
        control_min: Arc<Control>,
        control_max: Arc<Control>,
    ) -> BoilerServiceWaterRegular {
        BoilerServiceWaterRegular::new(
            Arc::new(RwLock::new(boiler_for_regular)),
            "boiler_test".into(),
            control_min,
            control_max,
        )
        .unwrap()
    }

    #[rstest]
    fn regular_boiler_should_provide_demand_hot_water(
        regular_boiler: BoilerServiceWaterRegular,
        simulation_time: SimulationTime,
    ) {
        let temp_return_feed = [51.05, 60.00];
        for (idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                regular_boiler
                    .demand_energy(
                        [0.7241412, 0.1748878][idx],
                        Default::default(),
                        Some(temp_return_feed[idx]),
                        None,
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

    #[rstest]
    fn test_temp_setpnt(
        regular_boiler: BoilerServiceWaterRegular,
        simulation_time: SimulationTime,
    ) {
        for t_it in simulation_time.iter() {
            assert_eq!(regular_boiler.setpnt(t_it), (Some(52.), Some(60.)));
        }
    }

    #[fixture]
    fn boiler_data_for_service_space() -> HeatSourceWetDetails {
        HeatSourceWetDetails::Boiler {
            rated_power: 16.85,
            energy_supply: "mains gas".into(),
            energy_supply_aux: "mains elec".into(),
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
    fn simulation_time_for_service_space() -> SimulationTime {
        SimulationTime::new(0., 3., 1.)
    }

    #[fixture]
    fn boiler_for_service_space(
        boiler_data_for_service_space: HeatSourceWetDetails,
        external_conditions: ExternalConditions,
        simulation_time_for_service_space: SimulationTime,
    ) -> Boiler {
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(
                FuelType::MainsGas,
                simulation_time_for_service_space.total_steps(),
            )
            .build(),
        ));
        let energy_supply_aux = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(
                FuelType::Electricity,
                simulation_time_for_service_space.total_steps(),
            )
            .build(),
        ));
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
        boiler.create_service_connection("boiler_test").unwrap();

        boiler
    }

    #[fixture]
    fn control_for_service_space() -> Control {
        Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(21.0), Some(21.0), None],
            0,
            1.0,
            Default::default(),
            Default::default(),
            1.0,
        ))
    }

    #[fixture]
    fn service_space_boiler(
        boiler_for_service_space: Boiler,
        control_for_service_space: Control,
    ) -> BoilerServiceSpace {
        BoilerServiceSpace::new(
            Arc::new(RwLock::new(boiler_for_service_space)),
            "boiler_test".into(),
            Arc::new(control_for_service_space),
        )
    }

    #[rstest]
    fn service_space_boiler_should_provide_demand_hot_water(
        service_space_boiler: BoilerServiceSpace,
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
                        Some(temp_return_feed[idx]),
                        None,
                        None,
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

    #[rstest]
    fn test_create_service_connection(
        #[from(boiler)] (mut boiler, _energy_supply): (Boiler, Arc<RwLock<EnergySupply>>),
    ) {
        let service_name = "new_service";
        // Ensure the service name does not already exist in energy supply connections
        assert!(!boiler.energy_supply_connections.contains_key(service_name));
        // Call the method under test
        boiler.create_service_connection(service_name).unwrap();
        // Check that the service name was added to enercy supply connections
        assert!(boiler.energy_supply_connections.contains_key(service_name));
        // Check there is an error when connection is attempted with existing service name
        assert!(boiler.create_service_connection(service_name).is_err());
    }

    #[rstest]
    fn test_create_service_hot_water_combi(
        #[from(boiler)] (boiler, _): (Boiler, Arc<RwLock<EnergySupply>>),
    ) {
        let service_name = "service_hot_water_combi";
        let cold_feed = ColdWaterSource::new(vec![1.0, 1.2], 0, 1.);
        let temp_hot_water = 50.;
        let boiler_data: HotWaterSourceDetails = serde_json::from_value(json!({
            "type": "CombiBoiler",
            "ColdWaterSource": "mains water",
            "HeatSourceWet": "hp",
            "separate_DHW_tests": "M&L",
            "rejected_energy_1": 0.0004,
            "storage_loss_factor_2": 0.91574,
            "rejected_factor_3": 0,
            "daily_HW_usage": 120,
            "setpoint_temp": 60.0
        }))
        .unwrap();

        let boiler = Arc::new(RwLock::new(boiler));

        let boiler_service_result = Boiler::create_service_hot_water_combi(
            boiler,
            boiler_data,
            service_name,
            temp_hot_water,
            WaterSourceWithTemperature::ColdWaterSource(Arc::new(cold_feed)),
        );
        assert!(boiler_service_result.is_ok());
    }

    #[rstest]
    fn test_create_service_hot_water_regular(
        #[from(boiler)] (boiler, _): (Boiler, Arc<RwLock<EnergySupply>>),
    ) {
        let service_name = "service_hot_water_regular";
        let control_min = Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            vec![None, None],
            0,
            1.0,
            Default::default(),
            Default::default(),
            1.0,
        )));
        let control_max = Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            vec![None, None],
            0,
            1.0,
            Default::default(),
            Default::default(),
            1.0,
        )));

        let boiler = Arc::new(RwLock::new(boiler));

        let boiler_hotwater_regular_result = Boiler::create_service_hot_water_regular(
            boiler,
            service_name,
            control_min,
            control_max,
        );
        assert!(boiler_hotwater_regular_result.is_ok());
    }

    #[rstest]
    fn test_create_service_space_heating(
        #[from(boiler)] (boiler, _): (Boiler, Arc<RwLock<EnergySupply>>),
    ) {
        let boiler = Arc::new(RwLock::new(boiler));

        let boiler_service_space_heating = Boiler::create_service_space_heating(
            boiler,
            "BoilerServiceSpace",
            Arc::new(Control::SetpointTime(SetpointTimeControl::new(
                vec![None, None],
                0,
                1.0,
                Default::default(),
                Default::default(),
                1.0,
            ))),
        );
        assert_eq!(
            type_of(boiler_service_space_heating),
            type_name::<BoilerServiceSpace>()
        );
    }

    // auxiliary method to check type - this is a little against the spirit of rust, but given for parity with the Python
    fn type_of<T>(_: T) -> &'static str {
        type_name::<T>()
    }

    #[rstest]
    fn test_cycling_adjustment(#[from(boiler)] (boiler, _): (Boiler, Arc<RwLock<EnergySupply>>)) {
        assert_relative_eq!(
            boiler.cycling_adjustment(40.0, 0.05, 0.5, 20.),
            0.030120066786994828,
            max_relative = 1e-7
        );
    }

    #[rstest]
    fn test_location_adjustment(
        #[from(boiler)] (boiler, energy_supply): (Boiler, Arc<RwLock<EnergySupply>>),
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) {
        // Internal boiler settings
        assert_relative_eq!(
            boiler.location_adjustment(30., 5., 20.),
            76.72260667117162,
            max_relative = 1e-7
        );

        let boiler_external_data: HeatSourceWetDetails = serde_json::from_value(json!({
            "type": "Boiler",
            "rated_power": 24.0,
            "EnergySupply": "mains_gas",
            "EnergySupply_aux": "Boiler_auxiliary", // added into test data here (compared to Python) as required for input
            "efficiency_full_load": 0.88,
            "efficiency_part_load": 0.986,
            "boiler_location": "external",
            "modulation_load" : 0.2,
            "electricity_circ_pump": 0.0600,
            "electricity_part_load" : 0.0131,
            "electricity_full_load" : 0.0388,
            "electricity_standby" : 0.0244
        }))
        .unwrap();
        let energy_supply_conn_aux =
            EnergySupply::connection(energy_supply.clone(), "Boiler_auxiliary").unwrap();

        let boiler_external = Boiler::new(
            boiler_external_data,
            energy_supply,
            energy_supply_conn_aux,
            Arc::new(external_conditions),
            simulation_time.step,
        )
        .unwrap();

        assert_relative_eq!(
            boiler_external.location_adjustment(30., 5., 2.5),
            31.530735570835994,
            max_relative = 1e-7
        );
    }

    #[rstest]
    fn test_calc_current_boiler_power(
        #[from(boiler)] (boiler, _): (Boiler, Arc<RwLock<EnergySupply>>),
    ) {
        assert_eq!(boiler.calc_current_boiler_power(10., 0.), 0.0);
        assert_relative_eq!(
            boiler.calc_current_boiler_power(10., 3.),
            4.800000000000001,
            max_relative = 1e-7
        );
    }

    #[rstest]
    fn test_calc_boiler_eff(
        #[from(boiler)] (boiler, energy_supply): (Boiler, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
        external_conditions: ExternalConditions,
    ) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                boiler
                    .calc_boiler_eff(false, 37., 3., None, Some(0.), t_it)
                    .unwrap(),
                [0.8619648380446757, 0.8619648380446757][t_idx],
                max_relative = 1e-7
            );
        }

        let boiler_external_data: HeatSourceWetDetails = serde_json::from_value(json!({
            "type": "Boiler",
            "rated_power": 24.0,
            "EnergySupply": "mains_gas",
            "EnergySupply_aux": "Boiler_auxiliary", // added into test data here (compared to Python) as required for input
            "efficiency_full_load": 0.88,
            "efficiency_part_load": 0.986,
            "boiler_location": "external",
            "modulation_load" : 0.2,
            "electricity_circ_pump": 0.0600,
            "electricity_part_load" : 0.0131,
            "electricity_full_load" : 0.0388,
            "electricity_standby" : 0.0244
        }))
        .unwrap();
        let energy_supply_conn_aux =
            EnergySupply::connection(energy_supply.clone(), "Boiler_auxiliary").unwrap();

        let boiler_external = Boiler::new(
            boiler_external_data,
            energy_supply,
            energy_supply_conn_aux,
            Arc::new(external_conditions),
            simulation_time.step,
        )
        .unwrap();

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_relative_eq!(
                boiler_external
                    .calc_boiler_eff(false, 37., 3., None, Some(0.), t_it)
                    .unwrap(),
                [0.8545088068138385, 0.8555283757048464][t_idx],
                max_relative = 1e-7
            );
        }
    }

    #[rstest]
    fn test_calc_energy_output_provided(
        #[from(boiler)] (boiler, _): (Boiler, Arc<RwLock<EnergySupply>>),
    ) {
        assert_eq!(boiler.calc_energy_output_provided(5.0, 1.), 5.0);
        assert_eq!(boiler.calc_energy_output_provided(25.0, 1.), 24.0);
    }

    #[rstest]
    fn test_time_available(
        #[from(boiler)] (boiler, _): (Boiler, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
    ) {
        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_eq!(boiler.time_available(0.0, None), [1.0, 1.0][t_idx]);
        }
        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_eq!(boiler.time_available(0.2, Some(0.5)), [0.4, 0.4][t_idx])
        }
    }

    #[rstest]
    fn test_demand_energy(
        #[from(boiler)] (mut boiler, _): (Boiler, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
    ) {
        boiler
            .create_service_connection("boiler_demand_energy")
            .unwrap();

        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_eq!(
                boiler
                    .demand_energy(
                        "boiler_demand_energy",
                        ServiceType::WaterCombi,
                        10.,
                        Some(37.),
                        None,
                        Some(false),
                        None,
                        None,
                    )
                    .unwrap()
                    .0,
                [10.0, 0.0][t_idx]
            );
        }

        boiler
            .create_service_connection("boiler_demand_energy_with_hybrid")
            .unwrap();

        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_eq!(
                boiler
                    .demand_energy(
                        "boiler_demand_energy_with_hybrid",
                        ServiceType::Space,
                        100.,
                        Some(37.),
                        None,
                        Some(true),
                        Some(0.),
                        None,
                    )
                    .unwrap(),
                [(24.0, Some(1.0)), (24.0, Some(1.0))][t_idx]
            );
        }

        // Test with time_elapsed_hp
        boiler
            .create_service_connection("boiler_demand_energy_hybrid_time_elapsed")
            .unwrap();

        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_eq!(
                boiler
                    .demand_energy(
                        "boiler_demand_energy_hybrid_time_elapsed",
                        ServiceType::Space,
                        100.,
                        Some(37.),
                        None,
                        Some(true),
                        Some(0.5),
                        None,
                    )
                    .unwrap(),
                [(12.0, Some(0.5)), (12.0, Some(0.5))][t_idx]
            );
        }

        let required_services = [
            "boiler_demand_energy",
            "boiler_demand_energy_with_hybrid",
            "boiler_demand_energy_hybrid_time_elapsed",
        ];
        let service_names_in_list = boiler
            .service_results
            .read()
            .iter()
            .map(|result| result.service_name.clone())
            .collect_vec();
        assert!(required_services
            .iter()
            .all(|&service| service_names_in_list.iter().any(|x| x.as_str() == service)));
    }

    // Python contains some further assertions using mocked boiler methods, which is difficult to do in Rust
    // without littering the implementation with test-specific overrides - deciding that this isn't worth
    // the trade-off here, at least for now
    // the Python method is called test_fuel_demand

    #[rstest]
    fn test_calc_auxiliary_energy(
        simulation_time: SimulationTime,
        boiler_data: HeatSourceWetDetails,
        #[from(boiler)] (_, energy_supply): (Boiler, Arc<RwLock<EnergySupply>>),
        external_conditions: ExternalConditions,
    ) {
        let external_conditions = Arc::new(external_conditions);

        let energy_supply_aux = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
        let energy_supply_conn_auxiliary =
            EnergySupply::connection(energy_supply_aux.clone(), "boiler_auxiliary").unwrap();

        let mut boiler = Boiler::new(
            boiler_data,
            energy_supply.clone(),
            energy_supply_conn_auxiliary,
            external_conditions.clone(),
            simulation_time.step,
        )
        .unwrap();

        // Check the function runs without panicking
        boiler.calc_auxiliary_energy(1., 0);

        // in Python there is some use of mocking here, which does not seem worth porting due to
        // the disproportionate difficulty in doing this vs the benefit of the assertion provided
    }

    #[rstest]
    fn test_timestep_end(
        #[from(boiler)] (mut boiler, _): (Boiler, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
    ) {
        boiler
            .create_service_connection("boiler_demand_energy")
            .unwrap();

        boiler
            .demand_energy(
                "boiler_demand_energy",
                ServiceType::WaterCombi,
                10.,
                Some(60.),
                None,
                Some(false),
                None,
                None,
            )
            .unwrap();

        assert_eq!(
            boiler
                .total_time_running_current_timestep
                .load(Ordering::SeqCst),
            1.0
        );
        assert_eq!(
            boiler.service_results.read()[0].service_name.as_str(),
            "boiler_demand_energy"
        );

        // Call the method under test
        boiler
            .timestep_end(simulation_time.iter().next().unwrap())
            .unwrap();

        // Assertions to check if the internal state was updated correctly
        assert_eq!(
            boiler
                .total_time_running_current_timestep
                .load(Ordering::SeqCst),
            0.0
        );
        assert!(boiler.service_results.read().is_empty());
    }

    #[rstest]
    fn test_energy_output_max(#[from(boiler)] (boiler, _): (Boiler, Arc<RwLock<EnergySupply>>)) {
        assert_eq!(boiler.energy_output_max(None, Some(0.)), 24.0);
        assert_eq!(boiler.energy_output_max(None, Some(0.5)), 12.0);
    }

    #[rstest]
    fn test_effvsreturntemp(
        #[from(boiler)] (boiler, _): (Boiler, Arc<RwLock<EnergySupply>>),
        simulation_time: SimulationTime,
        boiler_data: HeatSourceWetDetails,
        external_conditions: ExternalConditions,
    ) {
        let external_conditions = Arc::new(external_conditions);

        let efficiency = boiler
            .boiler_efficiency_over_return_temperatures(50., 0.5)
            .unwrap();
        let expected_efficiency = (-0.00007 * 50.0f64.powi(2) + 0.0017 * 50. + 0.979) - 0.5;
        assert_relative_eq!(efficiency, expected_efficiency, max_relative = 1e-7);

        let efficiency = boiler
            .boiler_efficiency_over_return_temperatures(60., 0.5)
            .unwrap();
        let expected_efficiency = (-0.0006 * 60.0 + 0.9129) - 0.5;
        assert_relative_eq!(efficiency, expected_efficiency, max_relative = 1e-7);

        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::LpgBulk, simulation_time.total_steps()).build(),
        ));
        let energy_supply_connection_aux =
            EnergySupply::connection(energy_supply.clone(), "boiler_lpg_bulk").unwrap();
        let boiler_lpg = Boiler::new(
            boiler_data.clone(),
            energy_supply.clone(),
            energy_supply_connection_aux.clone(),
            external_conditions.clone(),
            simulation_time.step,
        )
        .unwrap();

        let efficiency = boiler_lpg
            .boiler_efficiency_over_return_temperatures(45., 0.5)
            .unwrap();
        let expected_efficiency = (-0.00006 * 45.0f64.powi(2) + 0.0013 * 45. + 0.9859) - 0.5;
        assert_relative_eq!(efficiency, expected_efficiency, max_relative = 1e-7);

        let efficiency = boiler_lpg
            .boiler_efficiency_over_return_temperatures(50., 0.5)
            .unwrap();
        let expected_efficiency = (-0.0006 * 50.0 + 0.933) - 0.5;
        assert_relative_eq!(efficiency, expected_efficiency, max_relative = 1e-7);

        // Python here has a check for handling bad fuel codes, which are inexpressible in Rust due to use of enum (good thing!)

        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::MainsGas, simulation_time.total_steps()).build(),
        ));
        let boiler = Boiler::new(
            boiler_data,
            energy_supply.clone(),
            energy_supply_connection_aux,
            external_conditions.clone(),
            simulation_time.step,
        )
        .unwrap();
        assert_relative_eq!(
            boiler
                .boiler_efficiency_over_return_temperatures(50., 0.)
                .unwrap(),
            0.889,
            max_relative = 1e-3
        );
    }

    #[rstest]
    fn test_high_value_correction_part_load(
        #[from(boiler)] (boiler, _): (Boiler, Arc<RwLock<EnergySupply>>),
    ) {
        assert_eq!(
            Boiler::high_value_correction_part_load(&boiler.fuel_code, 5.).unwrap(),
            1.08
        );
        assert_relative_eq!(
            Boiler::high_value_correction_part_load(&boiler.fuel_code, 1.).unwrap(),
            0.99275,
            max_relative = 1e-5
        );

        assert_eq!(
            Boiler::high_value_correction_part_load(&FuelType::LpgBulk, 5.).unwrap(),
            1.06
        );
        assert_relative_eq!(
            Boiler::high_value_correction_part_load(&FuelType::LpgBulk, 1.).unwrap(),
            0.992758,
            max_relative = 1e-5
        );

        // Python here has a check for handling bad fuel codes, which are inexpressible in Rust due to use of enum (good thing!)
    }

    #[rstest]
    fn test_high_value_correction_full_load() {
        assert_relative_eq!(
            Boiler::high_value_correction_full_load(1.0),
            0.969715,
            max_relative = 1e-5
        );
    }

    #[rstest]
    fn test_net_to_gross(#[from(boiler)] (boiler, _): (Boiler, Arc<RwLock<EnergySupply>>)) {
        assert_eq!(Boiler::net_to_gross(&boiler.fuel_code).unwrap(), 0.901);

        assert_eq!(Boiler::net_to_gross(&FuelType::LpgBulk).unwrap(), 0.921);

        // Python here has a check for handling bad fuel codes, which are inexpressible in Rust due to use of enum (good thing!)
    }
}
