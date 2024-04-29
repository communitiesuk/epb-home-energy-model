use crate::compare_floats::{max_of_2, min_of_2};
use crate::core::common::WaterSourceWithTemperature;
use crate::core::controls::time_control::Control;
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::material_properties::MaterialProperties;
use crate::core::pipework::Pipework;
use crate::core::units::WATTS_PER_KILOWATT;
use crate::corpus::{HeatSource, PositionedHeatSource};
use crate::external_conditions::ExternalConditions;
use crate::input::{SolarCellLocation, WaterPipework};
use crate::simulation_time::SimulationTimeIteration;
use indexmap::IndexMap;
use parking_lot::Mutex;
use std::collections::HashMap;
use std::iter::zip;
use std::sync::Arc;

// BS EN 15316-5:2017 Appendix B default input data
// Model Information
// number of volumes the storage is modelled with
// see App.C (C.1.2 selection of the number of volumes to model the storage unit)
// for more details if this wants to be changed.
const STORAGE_TANK_NB_VOL: usize = 4;

// Product Description Data
// factors for energy recovery Table B.3
// part of the auxiliary energy transmitted to the medium
const STORAGE_TANK_F_RVD_AUX: f64 = 0.25;

// part of the thermal losses transmitted to the room
const STORAGE_TANK_F_STO_M: f64 = 0.75;

// ambient temperature - degrees
const STORAGE_TANK_TEMP_AMB: f64 = 16.;

#[derive(Clone, Debug)]
pub enum HeatSourceWithStorageTank {
    Immersion(Arc<Mutex<ImmersionHeater>>),
    Solar(Arc<Mutex<SolarThermalSystem>>),
}

/// An object to represent a hot water storage tank/cylinder
///
/// Models the case where hot water is drawn off and replaced by fresh cold
/// water which is then heated in the tank by a heat source. Assumes the water
/// is stratified by temperature.
///
/// Implements function demand_hot_water(volume_demanded) which all hot water
/// source objects must implement.
#[derive(Debug)]
pub struct StorageTank {
    q_std_ls_ref: f64, // measured standby losses due to cylinder insulation at standardised conditions, in kWh/24h
    temp_out_w_min: f64, // minimum temperature required for DHW (domestic hot water)
    temp_set_on: f64,  // set point temperature
    cold_feed: WaterSourceWithTemperature,
    simulation_timestep: f64,
    energy_supply_connection_unmet_demand: Option<EnergySupplyConnection>,
    control_hold_at_setpoint: Option<Arc<Control>>,
    volume_total_in_litres: f64,
    vol_n: [f64; STORAGE_TANK_NB_VOL],
    cp: f64,  // contents (usually water) specific heat in kWh/kg.K
    rho: f64, // volumic mass in kg/litre
    temp_n: [f64; STORAGE_TANK_NB_VOL],
    input_energy_adj_prev_timestep: f64,
    primary_pipework: Option<Pipework>,
    heat_source_data: IndexMap<String, PositionedHeatSource>, // heat sources, sorted by heater position
    heating_active: HashMap<String, bool>,
    q_ls_n_prev_heat_source: [f64; STORAGE_TANK_NB_VOL],
    q_sto_h_ls_rbl: Option<f64>, // total recoverable heat losses for heating in kWh, memoised between steps
    primary_gains: f64,          // primary pipework gains for a timestep (mutates over lifetime)
    #[cfg(test)]
    energy_demand_test: f64,
}

impl StorageTank {
    /// Arguments:
    /// * `volume` - total volume of the tank, in litres
    /// * `losses` - measured standby losses due to cylinder insulation
    ///                                at standardised conditions, in kWh/24h
    /// * `min_temp` - minimum temperature required for DHW
    /// * `setpoint_temp` - set point temperature
    /// * `cold_feed` - reference to ColdWaterSource object
    /// * `simulation_timestep` - the timestep for the simulation time being used in the calculation
    /// * `heat_sources`     -- hashmap of names and heat source objects
    /// * `primary_pipework` - optional reference to pipework
    /// * `energy_supply_connection_unmet_demand` - an energy supply connection representing unmet demand
    /// * `control_hold_at_setpnt` - reference to Control object with Boolean schedule
    ///                               defining when the StorageTank should be held at
    ///                                the setpoint temperature and not allowed to fall
    ///                               to the minimum before recharging
    /// * `contents` - MaterialProperties object
    pub fn new(
        volume: f64,
        losses: f64,
        min_temp: f64,
        setpoint_temp: f64,
        cold_feed: WaterSourceWithTemperature,
        simulation_timestep: f64,
        heat_sources: IndexMap<String, PositionedHeatSource>,
        primary_pipework: Option<WaterPipework>,
        energy_supply_connection_unmet_demand: Option<EnergySupplyConnection>,
        control_hold_at_setpoint: Option<Arc<Control>>,
        contents: MaterialProperties,
    ) -> Self {
        let q_std_ls_ref = losses;
        let temp_out_w_min = min_temp;
        let temp_set_on = setpoint_temp;

        let volume_total_in_litres = volume;
        // list of volume of layers in litres
        let vol_n = [volume_total_in_litres / STORAGE_TANK_NB_VOL as f64; STORAGE_TANK_NB_VOL];
        // water specific heat in kWh/kg.K
        let cp = contents.specific_heat_capacity_kwh();
        let rho = contents.density();

        // 6.4.3.2 STEP 0 Initialization
        //   for initial conditions all temperatures in the thermal storage unit(s)
        //   are equal to the set point temperature in degrees.
        //   We are expecting to run a "warm-up" period for the main calculation so this doesn't matter.
        let temp_n = [temp_set_on; STORAGE_TANK_NB_VOL];

        #[cfg(test)]
        let energy_demand_test = 0.;

        let input_energy_adj_prev_timestep = 0.;

        let primary_pipework: Option<Pipework> = primary_pipework.map(|pipework| pipework.into());

        let heating_active: HashMap<String, bool> = heat_sources
            .iter()
            .map(|(name, _heat_source)| ((*name).clone(), false))
            .collect();

        Self {
            q_std_ls_ref,
            temp_out_w_min,
            temp_set_on,
            cold_feed,
            simulation_timestep,
            energy_supply_connection_unmet_demand,
            control_hold_at_setpoint,
            volume_total_in_litres,
            vol_n,
            cp,
            rho,
            temp_n,
            input_energy_adj_prev_timestep,
            primary_pipework,
            heat_source_data: heat_sources,
            heating_active,
            q_ls_n_prev_heat_source: Default::default(),
            q_sto_h_ls_rbl: Default::default(),
            primary_gains: Default::default(),
            #[cfg(test)]
            energy_demand_test,
        }
    }

    pub fn get_cold_water_source(&self) -> &WaterSourceWithTemperature {
        &self.cold_feed
    }

    fn get_setpoint_min(&self, simtime: SimulationTimeIteration) -> f64 {
        match &self.control_hold_at_setpoint {
            Some(control) if control.is_on(simtime) => self.temp_set_on,
            _ => self.temp_out_w_min,
        }
    }

    /// Appendix B B.2.8 Stand-by losses are usually determined in terms of energy losses during
    /// a 24h period. Formula (B.2) allows the calculation of _sto_stbl_ls_tot based on a reference
    /// value of the daily thermal energy losses.
    ///
    /// h_sto_ls is the stand-by losses, in kW/K
    pub fn stand_by_losses_coefficient(&self) -> f64 {
        // BS EN 12897:2016 appendix B B.2.2
        // temperature of the water in the storage for the standardized conditions - degrees
        // these are reference (ref) temperatures from the standard test conditions for cylinder loss.
        let temp_set_ref = 65.;
        let temp_amb_ref = 20.;

        (1000. * self.q_std_ls_ref) / (24. * (temp_set_ref - temp_amb_ref))
    }

    /// Calculate the energy stored for each layer in the storage volume - kWh
    ///
    /// The energy stored is calculated, for information, accordingly to the limit value of
    /// temperature for domestic hot water.
    pub fn energy_stored(&self, timestep_idx: usize) -> [f64; STORAGE_TANK_NB_VOL] {
        self
            .temp_n
            .iter()
            .enumerate()
            .map(|(i, temp_i)| {
                let cold_feed_temperature = self.cold_feed.temperature(timestep_idx);
                if *temp_i > cold_feed_temperature {
                    self.rho * self.cp * self.vol_n[i] * (*temp_i - cold_feed_temperature)
                } else {
                    0.
                }
            })
            .collect::<Vec<f64>>().try_into().expect("Unexpected difficulty in calculating energy stored for a storage tank encountered.")
    }

    /// Convert the volume (in litres) demanded into an energy required in kWh
    pub fn energy_required(&self, volume_demanded: f64, timestep_idx: usize) -> f64 {
        self.rho
            * self.cp
            * volume_demanded
            * (self.temp_out_w_min - self.cold_feed.temperature(timestep_idx))
    }

    /// The calculation of the volume to be withdrawn is made accordingly with the energy to be
    /// delivered to the distribution system with a threshold value for the minimum available
    /// temperature according to the scenarios for domestic hot water.
    ///
    /// In this model only domestic hot water is considered for now.
    ///
    /// The volume of water withdrawn is based on contribution of the homogenous volumes of the
    /// storage unit, from the volume connected to the water output to the volume connected
    /// to the water input.
    pub fn energy_withdrawn(
        &self,
        q_out_w_dis_req: f64,
        q_out_w_n: [f64; STORAGE_TANK_NB_VOL],
        timestep_idx: usize,
    ) -> ([f64; STORAGE_TANK_NB_VOL], f64, [f64; STORAGE_TANK_NB_VOL]) {
        // initialise list of volume(s) to be withdrawn in litres
        let mut vol_use_w_n = [0.; STORAGE_TANK_NB_VOL];
        // initialise list of energy used for DHW in kWh
        let mut q_use_w_n = [0.; STORAGE_TANK_NB_VOL];
        // initialise tracker for energy required remainder / unmet energy
        let mut q_out_w_dis_req_rem = q_out_w_dis_req;

        //iterate in reverse order - from top of tank where draw off happens
        for (i, vol_i) in self.vol_n.iter().enumerate().rev() {
            if i == STORAGE_TANK_NB_VOL - 1 {
                // threshold minimum temperature
                if self.temp_n[i] >= self.temp_out_w_min {
                    // total energy to be delivered can be met by top layer
                    if q_out_w_dis_req <= q_out_w_n[i] {
                        vol_use_w_n[i] = q_out_w_dis_req
                            / (self.rho
                                * self.cp
                                * (self.temp_n[i] - self.cold_feed.temperature(timestep_idx)));
                        q_use_w_n[i] = q_out_w_dis_req;
                        q_out_w_dis_req_rem = 0.;
                        // no need to carry on as energy required has been met
                        break;
                    } else {
                        // top layer cannot meet all energy required
                        // so all of top layer volume will be withdrawn
                        vol_use_w_n[i] = *vol_i;
                        q_use_w_n[i] = q_out_w_n[i];
                        // update remaining energy still required from lower layers
                        q_out_w_dis_req_rem -= q_out_w_n[i];
                    }
                } else {
                    // temperature not met by top layer so no volume will be withdrawn
                    break;
                }
            } else if self.temp_n[i] >= self.temp_out_w_min {
                // now iterate over lower layers in turn
                if q_out_w_dis_req_rem <= q_out_w_n[i] {
                    // this layer can meet and/or exceed remainder energy required for distribution
                    vol_use_w_n[i] = q_out_w_dis_req_rem
                        / (self.rho
                            * self.cp
                            * (self.temp_n[i] - self.cold_feed.temperature(timestep_idx)));
                    q_use_w_n[i] = q_out_w_dis_req_rem;
                    q_out_w_dis_req_rem = 0.;
                    // no need to carry on as energy required has been met
                    break;
                } else if q_out_w_n[i] > 0. {
                    // this layer cannot meet remainder energy required
                    // so all of this layer volume will be withdrawn
                    vol_use_w_n[i] = *vol_i;
                    q_use_w_n[i] = q_out_w_n[i];
                    // update remaining energy still required from lower layers
                    q_out_w_dis_req_rem -= q_out_w_n[i];
                }
            }
        }

        (q_use_w_n, q_out_w_dis_req_rem, vol_use_w_n)
    }

    /// Principles: the volume withdrawn is replaced with the identical quantity of water
    /// provided to the input of the storage heater (bottom). The water of the upper volume is
    /// melted with the quantity of withdrawn water at the temperature of the lower level.
    pub fn volume_withdrawn_replaced(
        &self,
        vol_use_w_n: [f64; STORAGE_TANK_NB_VOL],
        timestep_idx: usize,
    ) -> [f64; STORAGE_TANK_NB_VOL] {
        // initialise list of temperature of layers AFTER volume withdrawn in degrees
        let mut temp_s3_n = self.temp_n;
        // initialise volume in each layer remaining after draw-off
        let mut v_sto_rem_n: [f64; STORAGE_TANK_NB_VOL] = zip(self.vol_n, vol_use_w_n)
            .map(|(x, y)| x - y)
            .collect::<Vec<f64>>()
            .try_into()
            .expect("Could not form a fixed-sized list of values while calculating volume_withdrawn_replaced");

        // temperature change only applicable if there is any volume withdrawn
        if vol_use_w_n.iter().sum::<f64>() > 0. {
            // determine how much water is displaced
            // IMPORTANT to iterate in reverse order - from top of tank
            for (i, _vol_i) in self.vol_n.iter().enumerate().rev() {
                let mut vol_to_replace = self.vol_n[i];
                // set list of flags for which layers need mixing for this layer
                let mut vol_mix_n = [0.; STORAGE_TANK_NB_VOL];
                // loop through layers i and below
                // (cleaner not to follow the Python here and iterate over v_sto_rem_n slice
                // as values from enumeration aren't used yet v_sto_rem_n is mutated in loop,
                // so we can just loop over a reversed range of the correct size instead)
                for j in (0..(i + 1)).rev() {
                    if v_sto_rem_n[j] == 0. {
                        // do nothing
                    } else if v_sto_rem_n[j] >= vol_to_replace {
                        vol_mix_n[j] = vol_to_replace;
                        v_sto_rem_n[j] -= vol_to_replace;
                        vol_to_replace = 0.;
                        break;
                    } else if v_sto_rem_n[j] < vol_to_replace {
                        vol_mix_n[j] = v_sto_rem_n[j];
                        v_sto_rem_n[j] = 0.;
                        vol_to_replace -= vol_mix_n[j];
                    }
                }
                // any volume remainder after looping through all layers will be at cold water temp

                // calculate new temperature of layer
                // note 6.4.3.5 equation 9 has an error as adding temps to volumes (sic from Python comment)
                temp_s3_n[i] = ((self.cold_feed.temperature(timestep_idx) * vol_to_replace)
                    + ((0..vol_mix_n.len())
                        .map(|k| self.temp_n[k] * vol_mix_n[k])
                        .sum::<f64>()))
                    / self.vol_n[i];
            }
        }

        temp_s3_n
    }

    pub fn potential_energy_input(
        &mut self,
        temp_s3_n: [f64; STORAGE_TANK_NB_VOL],
        heat_source: Arc<Mutex<HeatSource>>,
        heat_source_name: &str,
        heater_layer: usize,
        thermostat_layer: usize,
        simulation_time: SimulationTimeIteration,
    ) -> [f64; STORAGE_TANK_NB_VOL] {
        // initialise list of potential energy input for each layer
        let mut q_x_in_n = [0.; STORAGE_TANK_NB_VOL];

        let heat_source = &mut *(heat_source.lock());

        let energy_potential = match heat_source {
            HeatSource::Storage(HeatSourceWithStorageTank::Solar(ref solar_heat_source)) => {
                // we are passing the storage tank object to the SolarThermal as this needs to call back the storage tank (sic from Python)
                solar_heat_source
                    .lock()
                    .energy_output_max(self, temp_s3_n, &simulation_time)
            }
            HeatSource::Storage(HeatSourceWithStorageTank::Immersion(immersion_heater)) => {
                // no demand from heat source if the temperature of the tank at the thermostat position is below the set point

                // trigger heating to start when temperature falls below the minimum
                if temp_s3_n[thermostat_layer] <= self.get_setpoint_min(simulation_time) {
                    self.heating_active
                        .entry(heat_source_name.to_string())
                        .and_modify(|e| {
                            *e = true;
                        });
                }

                if self.heating_active[heat_source_name] {
                    let mut energy_potential = immersion_heater
                        .lock()
                        .energy_output_max(simulation_time, false);

                    if !matches!(
                        heat_source,
                        HeatSource::Storage(HeatSourceWithStorageTank::Immersion(_))
                    ) {
                        let (primary_pipework_losses_kwh, _) =
                            self.primary_pipework_losses(energy_potential);
                        energy_potential -= primary_pipework_losses_kwh;
                    }

                    energy_potential
                } else {
                    0.
                }
            }
            // the source Python code expects to run the above arm for any heat source
            // that is not solar that is passed here - we may have to implement some kind of common trait/
            // generic handling if the heat source is ever not immersion or solar here (like boilers)
            _ => todo!(),
        };

        q_x_in_n[heater_layer] += energy_potential;

        q_x_in_n
    }

    /// Function added into Storage tank to be called by the Solar Thermal object.
    /// Calculates the impact on storage tank temperature due to the proposed energy input
    pub fn storage_tank_potential_effect(
        &self,
        energy_proposed: f64,
        temp_s3_n: [f64; STORAGE_TANK_NB_VOL],
    ) -> (f64, f64) {
        // assuming initially no water draw-off

        // initialise list of potential energy input for each layer
        let mut q_x_in_n: [f64; STORAGE_TANK_NB_VOL] = Default::default();

        // TODO (from Python) - ensure we are feeding in the correct volume
        q_x_in_n[0] = energy_proposed;

        let (_q_s6, temp_s6_n) = self.energy_input(temp_s3_n, q_x_in_n);

        // 6.4.3.9 STEP 7 Re-arrange the temperatures in the storage after energy input
        let (_q_h_sto_s7, temp_s7_n) = self.rearrange_temperatures(temp_s6_n);

        // TODO (from Python) Check [0] is bottom layer temp and that solar thermal inlet is top layer NB_VOL-1
        (temp_s7_n[0], temp_s7_n[STORAGE_TANK_NB_VOL - 1])
    }

    /// The input of energy(s) is (are) allocated to the specific location(s)
    /// of the input of energy.
    /// Note: for energy withdrawn from a heat exchanger, the energy is accounted negatively.
    ///
    /// For step 6, the addition of the temperature of volume 'i' and theoretical variation of
    /// temperature calculated according to formula (10) can exceed the set temperature defined
    /// by the control system of the storage unit.
    fn energy_input(
        &self,
        temp_s3_n: [f64; STORAGE_TANK_NB_VOL],
        q_x_in_n: [f64; STORAGE_TANK_NB_VOL],
    ) -> (f64, [f64; STORAGE_TANK_NB_VOL]) {
        // initialise list of theoretical variation of temperature of layers in degrees
        let mut delta_temp_n: [f64; STORAGE_TANK_NB_VOL] = Default::default();
        // initialise list of theoretical temperature of layers after input in degrees
        let mut temp_s6_n: [f64; STORAGE_TANK_NB_VOL] = Default::default();
        // output energy delivered by the storage in kWh - timestep dependent
        let q_sto_h_out_n: [f64; STORAGE_TANK_NB_VOL] = Default::default();

        for i in 0..self.vol_n.len() {
            delta_temp_n[i] =
                (q_x_in_n[i] + q_sto_h_out_n[i]) / (self.rho * self.cp * self.vol_n[i]);
            temp_s6_n[i] = temp_s3_n[i] + delta_temp_n[i];
        }

        let q_s6 = self.rho
            * self.cp
            * (0..self.vol_n.len())
                .map(|i| self.vol_n[i] * temp_s6_n[i])
                .sum::<f64>();

        (q_s6, temp_s6_n)
    }

    /// When the temperature of the volume i is higher than the one of the upper volume,
    /// then the 2 volumes are melted. This iterative process is maintained until the temperature
    /// of the volume i is lower or equal to the temperature of the volume i+1.
    fn rearrange_temperatures(
        &self,
        temp_s6_n: [f64; STORAGE_TANK_NB_VOL],
    ) -> ([f64; STORAGE_TANK_NB_VOL], [f64; STORAGE_TANK_NB_VOL]) {
        // set list of flags for which layers need mixing
        let mut mix_layer_n: [u8; STORAGE_TANK_NB_VOL] = Default::default();
        let mut temp_s7_n = temp_s6_n;

        // loop through layers from bottom to top, without including top layer;
        // this is because the top layer has no upper layer to compare to
        for i in 0..self.vol_n.len() - 1 {
            if temp_s7_n[i] > temp_s7_n[i + 1] {
                // set layers to mix
                mix_layer_n[i] = 1;
                mix_layer_n[i + 1] = 1;
                // mix temperatures of all applicable layers
                // note error in formula 12 in standard as adding temperature to volume
                // this is what I think they intended from the description (comment sic from Python code)
                let temp_mix = (0..self.vol_n.len())
                    .map(|k| self.vol_n[k] * temp_s7_n[k] * mix_layer_n[k] as f64)
                    .sum::<f64>()
                    / (0..self.vol_n.len())
                        .map(|l| self.vol_n[l] * mix_layer_n[l] as f64)
                        .sum::<f64>();
                // set same temperature for all applicable layers
                for j in 0..i + 2 {
                    if mix_layer_n[j] == 1 {
                        temp_s7_n[j] = temp_mix;
                    }
                }
            } else {
                // reset mixing as lower levels now stabilised
                mix_layer_n = Default::default();
            }
        }

        let q_h_sto_end: [f64; STORAGE_TANK_NB_VOL] = (0..self.vol_n.len())
            .map(|i| self.rho * self.cp * self.vol_n[i] * temp_s7_n[i])
            .collect::<Vec<f64>>()
            .try_into()
            .unwrap();

        let _q_h_sto_end_no = self.rho
            * self.cp
            * (0..self.vol_n.len())
                .map(|i| self.vol_n[i] * temp_s7_n[i])
                .sum::<f64>();

        (q_h_sto_end, temp_s7_n)
    }

    /// Thermal losses are calculated with respect to the impact of the temperature set point
    pub fn thermal_losses(
        &self,
        temp_s7_n: [f64; STORAGE_TANK_NB_VOL],
        q_x_in_n: [f64; STORAGE_TANK_NB_VOL],
        q_h_sto_s7: [f64; STORAGE_TANK_NB_VOL],
        heater_layer: usize,
        q_ls_n_prev_heat_source: [f64; STORAGE_TANK_NB_VOL],
    ) -> (
        f64,
        f64,
        [f64; STORAGE_TANK_NB_VOL],
        [f64; STORAGE_TANK_NB_VOL],
    ) {
        // standby losses coefficient - kW/K
        let h_sto_ls = self.stand_by_losses_coefficient();

        // standby losses correction factor - dimensionless
        // note from Python code: "do not think these are applicable so used: f_sto_dis_ls = 1, f_sto_bac_acc = 1"

        // initialise list of thermal losses in kWh
        let mut q_ls_n: [f64; STORAGE_TANK_NB_VOL] = Default::default();
        // initialise list of final temperature of layers after thermal losses in degrees
        let mut temp_s8_n: [f64; STORAGE_TANK_NB_VOL] = Default::default();

        // Thermal losses
        // Note: Eqn 13 from BS EN 15316-5:2017 does not explicitly multiply by
        // timestep (it seems to assume a 1 hour timestep implicitly), but it is
        // necessary to convert the rate of heat loss to a total heat loss over
        // the time period
        for i in 0..self.vol_n.len() {
            q_ls_n[i] = (h_sto_ls * self.rho * self.cp)
                * (self.vol_n[i] / self.volume_total_in_litres)
                * (min_of_2(temp_s7_n[i], self.temp_set_on) - STORAGE_TANK_TEMP_AMB)
                * self.simulation_timestep;
            q_ls_n[i] = max_of_2(0., q_ls_n[i] - q_ls_n_prev_heat_source[i]);
        }

        // total thermal losses kWh
        let q_ls = q_ls_n.iter().sum();

        // the final value of the temperature is reduced due to the effect of the thermal losses.
        // check temperature compared to set point
        // the temperature for each volume are limited to the set point for any volume controlled
        for i in 0..self.vol_n.len() {
            temp_s8_n[i] = if temp_s7_n[i] > self.temp_set_on {
                // Case 2 - Temperature exceeding the set point
                self.temp_set_on
            } else {
                // Case 1 - Temperature below the set point

                // the final value of the temperature
                // is reduced due to the effect of the thermal losses
                // Formula (14) in the standard appears to have error as addition not multiply
                // and P instead of rho
                temp_s7_n[i] - (q_ls_n[i] / (self.rho * self.cp * self.vol_n[i]))
            };
        }

        // excess energy/ energy surplus
        // excess energy is calculated as the difference from the energy stored, Qsto,step7, and
        // energy stored once the set temperature is obtained, Qsto,step8, with addition of the
        // thermal losses.
        // Note: The surplus must be calculated only for those layers that the
        //       heat source currently being considered is capable of heating,
        //       i.e. excluding those below the heater position.
        let energy_surplus = if temp_s7_n[heater_layer] > self.temp_set_on {
            (heater_layer..STORAGE_TANK_NB_VOL).fold(0., |acc, i| {
                acc + q_h_sto_s7[i]
                    - q_ls_n[i]
                    - (self.rho * self.cp * self.vol_n[i] * self.temp_set_on)
            })
        } else {
            0.
        };

        // the thermal energy provided to the system (from heat sources) shall be limited
        // adjustment of the energy delivered to the storage according with the set temperature
        // potential input from generation
        let q_x_in_adj = q_x_in_n.iter().sum::<f64>();
        // TODO (from Python code) - find in standard - availability of back-up - where is this from?

        // also referred to as electrical power on
        let sto_bu_on = 1.;
        let q_in_h_w = min_of_2(q_x_in_adj - energy_surplus, q_x_in_adj * sto_bu_on);

        (q_in_h_w, q_ls, temp_s8_n, q_ls_n)
    }

    // NB. there is a testoutput() function here in the Python to output to a test file - will not reimplement unless seen as necessary

    fn run_heat_sources(
        &mut self,
        temp_s3_n: [f64; STORAGE_TANK_NB_VOL],
        heat_source: Arc<Mutex<HeatSource>>,
        heat_source_name: &str,
        heater_layer: usize,
        thermostat_layer: usize,
        q_ls_prev_heat_source: [f64; STORAGE_TANK_NB_VOL],
        simulation_time: SimulationTimeIteration,
    ) -> TemperatureCalculation {
        // 6.4.3.8 STEP 6 Energy input into the storage
        // input energy delivered to the storage in kWh - timestep dependent
        let q_x_in_n = self.potential_energy_input(
            temp_s3_n,
            heat_source.clone(),
            heat_source_name,
            heater_layer,
            thermostat_layer,
            simulation_time,
        );

        self.calculate_temperatures(
            temp_s3_n,
            heat_source,
            q_x_in_n,
            heater_layer,
            q_ls_prev_heat_source,
            simulation_time,
        )
    }

    fn calculate_temperatures(
        &mut self,
        temp_s3_n: [f64; STORAGE_TANK_NB_VOL],
        heat_source: Arc<Mutex<HeatSource>>,
        q_x_in_n: [f64; STORAGE_TANK_NB_VOL],
        heater_layer: usize,
        q_ls_n_prev_heat_source: [f64; STORAGE_TANK_NB_VOL],
        simulation_time_iteration: SimulationTimeIteration,
    ) -> TemperatureCalculation {
        let (q_s6, temp_s6_n) = self.energy_input(temp_s3_n, q_x_in_n);

        // 6.4.3.9 STEP 7 Re-arrange the temperatures in the storage after energy input
        let (q_h_sto_s7, temp_s7_n) = self.rearrange_temperatures(temp_s6_n);

        // STEP 8 Thermal losses and final temperature
        let (q_in_h_w, q_ls, temp_s8_n, q_ls_n) = self.thermal_losses(
            temp_s7_n,
            q_x_in_n,
            q_h_sto_s7,
            heater_layer,
            q_ls_n_prev_heat_source,
        );

        // TODO 6.4.3.11 Heat exchanger

        // demand adjusted energy from heat source (before was just using potential without taking it)
        let input_energy_adj = q_in_h_w;

        #[cfg(test)]
        {
            self.energy_demand_test = input_energy_adj;
        }

        let _heat_source_output =
            self.heat_source_output(heat_source, input_energy_adj, simulation_time_iteration);
        // variable is updated in upstream but then never read
        // input_energy_adj -= _heat_source_output;

        (
            temp_s8_n, q_x_in_n, q_s6, temp_s6_n, temp_s7_n, q_in_h_w, q_ls, q_ls_n,
        )
    }

    /// Draw off hot water from the tank
    /// Energy calculation as per BS EN 15316-5:2017 Method A sections 6.4.3, 6.4.6, 6.4.7
    ///
    /// Arguments:
    /// * `volume_demanded` - volume of hot water required, in litres
    pub fn demand_hot_water(
        &mut self,
        volume_demanded: f64,
        simulation_time: SimulationTimeIteration,
    ) -> f64 {
        // 6.4.3.3 STEP 1 Calculate energy stored
        // energy stored for domestic hot water - kWh
        let q_out_w_n = self.energy_stored(simulation_time.index);

        // 6.4.3.4 STEP 2 Volume (and energy) to be withdrawn from the storage (for DHW)
        // energy required for domestic hot water in kWh
        let q_out_w_dis_req = self.energy_required(volume_demanded, simulation_time.index);
        // energy withdrawn, unmet energy required, volume withdrawn
        let (q_use_w_n, q_out_w_dis_req_rem, vol_use_w_n) =
            self.energy_withdrawn(q_out_w_dis_req, q_out_w_n, simulation_time.index);

        // if tank cannot provide enough hot water report unmet demand
        if let Some(energy_supply) = &self.energy_supply_connection_unmet_demand {
            energy_supply
                .demand_energy(q_out_w_dis_req_rem, simulation_time.index)
                .unwrap();
        }

        // 6.4.3.5 STEP 3 Temperature of the storage after volume withdrawn (for DHW)
        let temp_s3_n = self.volume_withdrawn_replaced(vol_use_w_n, simulation_time.index);

        // Run over multiple heat sources
        let mut temp_after_prev_heat_source = temp_s3_n;
        let mut q_ls: f64 = Default::default();
        self.q_ls_n_prev_heat_source = Default::default();
        // we need the last value of temp_s8_n after the loop for later in the function
        let mut temp_s8_n: [f64; STORAGE_TANK_NB_VOL] = Default::default();
        let mut heat_source_data = self.heat_source_data.clone(); // TODO see if we can avoid this allocation
        for (heat_source_name, positioned_heat_source) in heat_source_data.iter_mut() {
            let heater_layer =
                (positioned_heat_source.heater_position * STORAGE_TANK_NB_VOL as f64) as usize;
            let thermostat_layer =
                (positioned_heat_source.thermostat_position * STORAGE_TANK_NB_VOL as f64) as usize;

            let (
                temp_s8_n_step,
                _q_x_in_n,
                _q_s6,
                _temp_s6_n,
                _temp_s7_n,
                _q_in_h_w,
                q_ls_this_heat_source,
                q_ls_n_this_heat_source,
            ) = self.run_heat_sources(
                temp_after_prev_heat_source,
                positioned_heat_source.heat_source.clone(),
                heat_source_name,
                heater_layer,
                thermostat_layer,
                self.q_ls_n_prev_heat_source,
                simulation_time,
            );

            temp_after_prev_heat_source = temp_s8_n_step;
            q_ls += q_ls_this_heat_source;
            for (i, q_ls_n) in q_ls_n_this_heat_source.iter().enumerate() {
                self.q_ls_n_prev_heat_source[i] += *q_ls_n;
            }

            // Trigger heating to stop when setpoint is reached
            if temp_s8_n_step[thermostat_layer] >= self.temp_set_on {
                self.heating_active
                    .entry(heat_source_name.to_string())
                    .and_modify(|e| {
                        *e = false;
                    });
            }

            temp_s8_n = temp_s8_n_step;
        }

        // Additional calculations
        // 6.4.6 Calculation of the auxiliary energy
        // accounted for elsewhere so not included here
        let w_sto_aux = 0.;

        // 6.4.7 Recoverable, recovered thermal losses
        // recovered auxiliary energy to the heating medium - kWh
        let _q_sto_h_aux_rvd = w_sto_aux * STORAGE_TANK_F_RVD_AUX;
        // recoverable auxiliary energy transmitted to the heated space - kWh
        let q_sto_h_rbl_aux = w_sto_aux * STORAGE_TANK_F_STO_M * (1. - STORAGE_TANK_F_RVD_AUX);
        // recoverable heat losses (storage) - kWh
        let q_sto_h_rbl_env = q_ls * STORAGE_TANK_F_STO_M;
        // total recoverable heat losses for heating - kWh
        self.q_sto_h_ls_rbl = Some(q_sto_h_rbl_env + q_sto_h_rbl_aux);

        // set temperatures calculated to be initial temperatures of volumes for the next timestep
        self.temp_n = temp_s8_n;

        q_use_w_n.iter().sum()
    }

    fn additional_energy_input(
        &mut self,
        heat_source: Arc<Mutex<HeatSource>>,
        heat_source_name: &str,
        energy_input: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        if energy_input == 0. {
            return 0.;
        }

        let heat_source_data = &self.heat_source_data[heat_source_name];

        let heater_layer = (heat_source_data.heater_position * STORAGE_TANK_NB_VOL as f64) as usize;
        let _thermostat_layer =
            (heat_source_data.thermostat_position * STORAGE_TANK_NB_VOL as f64) as usize;

        let mut q_x_in_n: [f64; STORAGE_TANK_NB_VOL] = Default::default();
        q_x_in_n[heater_layer] = energy_input;
        let (temp_s8_n, _, _, _, _, q_in_h_w, _, q_ls_n_this_heat_source) = self
            .calculate_temperatures(
                self.temp_n,
                heat_source,
                q_x_in_n,
                heater_layer,
                self.q_ls_n_prev_heat_source,
                simulation_time_iteration,
            );

        for (i, q_ls_n) in q_ls_n_this_heat_source.iter().enumerate() {
            self.q_ls_n_prev_heat_source[i] += *q_ls_n;
        }

        self.temp_n = temp_s8_n;

        q_in_h_w
    }

    #[cfg(test)]
    fn test_energy_demand(&self) -> f64 {
        self.energy_demand_test
    }

    /// Return the DHW recoverable heat losses as internal gain for the current timestep in W
    pub fn internal_gains(&mut self) -> f64 {
        let primary_gains_timestep = self.primary_gains;
        self.primary_gains = Default::default();

        (self.q_sto_h_ls_rbl.expect(
            "storage tank logic expects q_sto_h_ls_rbl to have been set internally at this point",
        ) * WATTS_PER_KILOWATT as f64
            / self.simulation_timestep)
            + primary_gains_timestep
    }

    fn primary_pipework_losses(&self, input_energy_adj: f64) -> (f64, f64) {
        let mut primary_pipework_losses_kwh = Default::default();
        let mut primary_gains_w = Default::default();

        // start of heating event
        if input_energy_adj > 0. && self.input_energy_adj_prev_timestep == 0. {
            primary_pipework_losses_kwh += self.primary_pipework.as_ref().expect("primary pipework is expected to have been set on the storage tank at this point").cool_down_loss(self.temp_set_on, STORAGE_TANK_TEMP_AMB);
        }

        // during heating event
        if input_energy_adj > 0. {
            // primary losses for the timestep calculated from temperature difference
            let primary_pipework_losses_w = self
                .primary_pipework
                .as_ref()
                .expect("pipework expected to have been set")
                .heat_loss(self.temp_set_on, STORAGE_TANK_TEMP_AMB);
            primary_gains_w += primary_pipework_losses_w;
            primary_pipework_losses_kwh +=
                primary_pipework_losses_w * self.simulation_timestep / WATTS_PER_KILOWATT as f64;
        }

        // end of heating event
        if input_energy_adj == 0. && self.input_energy_adj_prev_timestep > 0. {
            primary_gains_w += self
                .primary_pipework
                .as_ref()
                .expect("pipework expected")
                .cool_down_loss(self.temp_set_on, STORAGE_TANK_TEMP_AMB);
        }

        (primary_pipework_losses_kwh, primary_gains_w)
    }

    fn heat_source_output(
        &mut self,
        heat_source: Arc<Mutex<HeatSource>>,
        input_energy_adj: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        match &mut *(heat_source.clone().lock()) {
            HeatSource::Storage(HeatSourceWithStorageTank::Immersion(immersion)) => immersion
                .lock()
                .demand_energy(input_energy_adj, simulation_time_iteration),
            HeatSource::Storage(HeatSourceWithStorageTank::Solar(solar)) => solar
                .lock()
                .demand_energy(input_energy_adj, simulation_time_iteration.index),
            _ => {
                // TODO need to be able to call demand_energy on the other heat sources
                let (primary_pipework_losses_kwh, primary_gains) =
                    self.primary_pipework_losses(input_energy_adj);
                let input_energy_adj = input_energy_adj + primary_pipework_losses_kwh;
                let heat_source_output = heat_source
                    .lock()
                    .demand_energy(input_energy_adj, simulation_time_iteration)
                    - primary_pipework_losses_kwh;
                self.input_energy_adj_prev_timestep = input_energy_adj;
                self.primary_gains = primary_gains;

                heat_source_output
            }
        }
    }
}

type TemperatureCalculation = (
    [f64; STORAGE_TANK_NB_VOL],
    [f64; STORAGE_TANK_NB_VOL],
    f64,
    [f64; STORAGE_TANK_NB_VOL],
    [f64; STORAGE_TANK_NB_VOL],
    f64,
    f64,
    [f64; STORAGE_TANK_NB_VOL],
);

#[derive(Clone, Debug)]
pub struct ImmersionHeater {
    pwr: f64, // rated power
    energy_supply_connection: EnergySupplyConnection,
    simulation_timestep: f64,
    control: Option<Arc<Control>>,
    diverter: Option<Arc<Mutex<PVDiverter>>>,
}

impl ImmersionHeater {
    pub fn new(
        rated_power: f64,
        energy_supply_connection: EnergySupplyConnection,
        simulation_timestep: f64,
        control: Option<Arc<Control>>,
    ) -> Self {
        Self {
            pwr: rated_power,
            energy_supply_connection,
            simulation_timestep,
            control,
            diverter: Default::default(),
        }
    }

    pub fn connect_diverter(&mut self, diverter: Arc<Mutex<PVDiverter>>) {
        if self.diverter.is_some() {
            panic!("diverter was already connected");
        }

        self.diverter = Some(diverter);
    }

    /// Demand energy (in kWh) from the heater
    pub fn demand_energy(&mut self, energy_demand: f64, simtime: SimulationTimeIteration) -> f64 {
        // Account for time control where present. If no control present, assume
        // system is always active (except for basic thermostatic control, which
        // is implicit in demand calculation).
        let energy_supplied =
            if self.control.is_none() || self.control.as_ref().unwrap().is_on(simtime) {
                min_of_2(energy_demand, self.pwr * self.simulation_timestep)
            } else {
                0.
            };

        // If there is a diverter to this immersion heater, then any heating
        // capacity already in use is not available to the diverter.
        if let Some(ref mut diverter) = &mut self.diverter {
            diverter.lock().capacity_already_in_use(energy_supplied);
        }

        self.energy_supply_connection
            .demand_energy(energy_supplied, simtime.index)
            .unwrap();

        energy_supplied
    }

    /// Calculate the maximum energy output (in kWh) from the heater
    pub fn energy_output_max(
        &self,
        simtime: SimulationTimeIteration,
        ignore_standard_control: bool,
    ) -> f64 {
        // Account for time control where present. If no control present, assume
        // system is always active (except for basic thermostatic control, which
        // is implicit in demand calculation).
        if self.control.is_none()
            || self.control.as_ref().unwrap().is_on(simtime)
            || ignore_standard_control
        {
            self.pwr * self.simulation_timestep
        } else {
            0.
        }
    }
}

#[derive(Clone, Debug)]
pub struct PVDiverter {
    storage_tank: Arc<Mutex<StorageTank>>,
    immersion_heater: Arc<Mutex<ImmersionHeater>>,
    heat_source_name: String,
    capacity_already_in_use: f64,
}

impl PVDiverter {
    pub fn new(
        storage_tank: Arc<Mutex<StorageTank>>,
        heat_source: Arc<Mutex<ImmersionHeater>>,
        heat_source_name: String,
    ) -> Arc<Mutex<Self>> {
        let diverter = Arc::new(Mutex::new(Self {
            storage_tank,
            heat_source_name,
            immersion_heater: heat_source.clone(),
            capacity_already_in_use: Default::default(),
        }));

        heat_source.lock().connect_diverter(diverter.clone());

        diverter
    }

    /// Record heater output that would happen anyway, to avoid double-counting
    pub fn capacity_already_in_use(&mut self, energy_supplied: f64) {
        self.capacity_already_in_use += energy_supplied;
    }

    /// Divert as much surplus as possible to the heater
    ///
    /// Arguments:
    /// * `supply_surplus` - surplus energy, in kWh, available to be diverted (negative by convention)
    pub fn divert_surplus(
        &mut self,
        supply_surplus: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        // check how much spare capacity the immersion heater has
        let imm_heater_max_capacity_spare = self
            .immersion_heater
            .lock()
            .energy_output_max(simulation_time_iteration, true)
            - self.capacity_already_in_use;

        // Calculate the maximum energy that could be diverted
        // Note: supply_surplus argument is negative by convention, so negate it here
        let energy_diverted_max = min_of_2(imm_heater_max_capacity_spare, -supply_surplus);

        // Add additional energy to storage tank and calculate how much energy was accepted
        let energy_diverted = self.storage_tank.lock().additional_energy_input(
            Arc::new(Mutex::new(HeatSource::Storage(
                HeatSourceWithStorageTank::Immersion(self.immersion_heater.clone()),
            ))),
            &self.heat_source_name,
            energy_diverted_max,
            simulation_time_iteration,
        );

        energy_diverted
    }

    pub fn timestep_end(&mut self) {
        self.capacity_already_in_use = Default::default();
    }
}

/// The following code contains objects that represent solar thermal systems.
/// Method 3 in BS EN 15316-4-3:2017.
#[derive(Clone, Debug)]
pub struct SolarThermalSystem {
    sol_loc: SolarCellLocation,
    area: f64,
    peak_collector_efficiency: f64,
    incidence_angle_modifier: f64,
    first_order_hlc: f64,
    second_order_hlc: f64,
    collector_mass_flow_rate: f64,
    power_pump: f64,
    power_pump_control: f64,
    energy_supply_connection: EnergySupplyConnection,
    tilt: f64,
    orientation: f64,
    solar_loop_piping_hlc: f64,
    external_conditions: Arc<ExternalConditions>,
    simulation_timestep: f64,
    heat_output_collector_loop: f64,
    energy_supplied: f64,
    cp: f64,
    air_temp_coll_loop: f64, //mutating internally
    inlet_temp: f64,
}

// BS EN 15316-4-3:2017 Appendix B default input data
// Model Information
// Air temperature in a heated space in the building
// Default taken from Table B20 of standard
const AIR_TEMP_HEATED_ROOM: f64 = 20.;

impl SolarThermalSystem {
    pub fn new(
        sol_loc: SolarCellLocation,
        area_module: f64,
        modules: usize,
        peak_collector_efficiency: f64,
        incidence_angle_modifier: f64,
        first_order_hlc: f64,
        second_order_hlc: f64,
        collector_mass_flow_rate: f64,
        power_pump: f64,
        power_pump_control: f64,
        energy_supply_connection: EnergySupplyConnection,
        tilt: f64,
        orientation: f64,
        solar_loop_piping_hlc: f64,
        external_conditions: Arc<ExternalConditions>,
        simulation_timestep: f64,
        contents: MaterialProperties,
    ) -> Self {
        Self {
            sol_loc,
            area: area_module * modules as f64,
            peak_collector_efficiency,
            incidence_angle_modifier,
            first_order_hlc,
            second_order_hlc,
            collector_mass_flow_rate,
            power_pump,
            power_pump_control,
            energy_supply_connection,
            tilt,
            orientation,
            solar_loop_piping_hlc,
            external_conditions,
            simulation_timestep,
            heat_output_collector_loop: 0.0,
            energy_supplied: 0.0,
            // Water specific heat in J/kg.K
            // (defined under eqn 51 on page 40 of BS EN ISO 15316-4-3:2017)
            cp: contents.specific_heat_capacity(),
            air_temp_coll_loop: Default::default(),
            inlet_temp: Default::default(),
        }
    }

    /// Calculate collector loop heat output
    /// eq 49 to 58 of STANDARD
    pub fn energy_output_max(
        &mut self,
        storage_tank: &StorageTank,
        temp_storage_tank_s3_n: [f64; STORAGE_TANK_NB_VOL],
        simulation_time: &SimulationTimeIteration,
    ) -> f64 {
        self.air_temp_coll_loop = match self.sol_loc {
            SolarCellLocation::Out => self.external_conditions.air_temp(simulation_time),
            SolarCellLocation::Hs => AIR_TEMP_HEATED_ROOM,
            SolarCellLocation::Nhs => {
                (AIR_TEMP_HEATED_ROOM + self.external_conditions.air_temp(simulation_time)) / 2.
            }
        };

        // First estimation of average collector water temperature. Eq 51
        // initialise temperature
        // if first time step, pick bottom of the tank temperature as inlet_temp_s1
        let inlet_temp_s1 = if simulation_time.index == 0 {
            let inlet_temp_s1 = temp_storage_tank_s3_n[0];
            self.inlet_temp = inlet_temp_s1;

            inlet_temp_s1
        } else {
            self.inlet_temp
        };

        // solar irradiance in W/m2
        let solar_irradiance = self.external_conditions.calculated_total_solar_irradiance(
            self.tilt,
            self.orientation,
            simulation_time,
        );

        if solar_irradiance == 0. {
            self.heat_output_collector_loop = 0.;
            return 0.;
        }

        let mut avg_collector_water_temp = inlet_temp_s1
            + (0.4 * solar_irradiance * self.area) / (self.collector_mass_flow_rate * self.cp * 2.);

        let mut inlet_temp2: f64 = Default::default(); // need a running slot in the loop for this to be overridden each time

        // calculation of collector efficiency
        for _ in 0..4 {
            // Eq 53
            let th = (avg_collector_water_temp
                - self.external_conditions.air_temp(simulation_time))
                / solar_irradiance;

            // Eq 52
            let collector_efficiency = self.peak_collector_efficiency
                * self.incidence_angle_modifier
                - self.first_order_hlc * th
                - self.second_order_hlc * th.powi(2) * solar_irradiance;

            // Eq 54
            let _collector_absorber_heat_input = self.peak_collector_efficiency
                * solar_irradiance
                * self.area
                * simulation_time.timestep
                / WATTS_PER_KILOWATT as f64;

            // Eq 55
            let collector_output_heat =
                collector_efficiency * solar_irradiance * self.area * simulation_time.timestep
                    / WATTS_PER_KILOWATT as f64;

            // Eq 56
            let heat_loss_collector_loop_piping = self.solar_loop_piping_hlc
                * (avg_collector_water_temp - self.air_temp_coll_loop)
                * simulation_time.timestep
                / WATTS_PER_KILOWATT as f64;

            // Eq 57
            self.heat_output_collector_loop =
                collector_output_heat - heat_loss_collector_loop_piping;
            if self.heat_output_collector_loop
                < self.power_pump * simulation_time.timestep * 3. / WATTS_PER_KILOWATT as f64
            {
                self.heat_output_collector_loop = 0.;
            }

            // Call to the storage tank
            let (_temp_layer_0, inlet_temp2_temp) = storage_tank.storage_tank_potential_effect(
                self.heat_output_collector_loop,
                temp_storage_tank_s3_n,
            );
            inlet_temp2 = inlet_temp2_temp;

            // Eq 58
            avg_collector_water_temp = (self.inlet_temp + inlet_temp2) / 2.
                + self.heat_output_collector_loop / (self.collector_mass_flow_rate * self.cp * 2.);
        }

        self.inlet_temp = inlet_temp2;

        self.heat_output_collector_loop
    }

    pub fn demand_energy(&mut self, energy_demand: f64, timestep_idx: usize) -> f64 {
        self.energy_supplied = min_of_2(energy_demand, self.heat_output_collector_loop);

        // Eq 59 and 60 to calculate auxiliary energy - note that the if condition
        // is the wrong way round in BS EN 15316-4-3:2017
        let auxiliary_energy_consumption = if self.energy_supplied == 0. {
            self.power_pump_control * self.simulation_timestep
        } else {
            (self.power_pump_control + self.power_pump) * self.simulation_timestep
        };

        self.energy_supply_connection
            .demand_energy(auxiliary_energy_consumption, timestep_idx)
            .unwrap();

        self.energy_supplied
    }

    #[cfg(test)]
    fn test_energy_potential(&self) -> f64 {
        self.heat_output_collector_loop
    }

    #[cfg(test)]
    fn test_energy_supplied(&self) -> f64 {
        self.energy_supplied
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::controls::time_control::OnOffTimeControl;
    use crate::core::energy_supply::energy_supply::EnergySupply;
    use crate::core::material_properties::WATER;
    use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
    use crate::corpus::HeatSource;
    use crate::external_conditions::{
        DaylightSavingsConfig, ShadingObject, ShadingObjectType, ShadingSegment,
    };
    use crate::input::FuelType;
    use crate::simulation_time::SimulationTime;
    use rstest::*;

    fn round_by_precision(src: f64, precision: f64) -> f64 {
        (precision * src).round() / precision
    }

    fn round_each_by_precision(
        src: Vec<[f64; STORAGE_TANK_NB_VOL]>,
        precision: f64,
    ) -> Vec<[f64; STORAGE_TANK_NB_VOL]> {
        src.iter()
            .map(|timestep_values| {
                timestep_values
                    .iter()
                    .map(|num| round_by_precision(*num, precision))
                    .collect::<Vec<f64>>()
                    .try_into()
                    .unwrap()
            })
            .collect()
    }

    #[fixture]
    pub fn simulation_time_for_storage_tank() -> SimulationTime {
        SimulationTime::new(0., 8., 1.)
    }

    #[fixture]
    pub fn control_for_storage_tank() -> Arc<Control> {
        Arc::new(Control::OnOffTimeControl(OnOffTimeControl::new(
            vec![true, false, false, false, true, true, true, true],
            0,
            1.,
        )))
    }

    #[fixture]
    pub fn cold_water_source(
        simulation_time_for_storage_tank: SimulationTime,
    ) -> Arc<ColdWaterSource> {
        Arc::new(ColdWaterSource::new(
            vec![10.0, 10.1, 10.2, 10.5, 10.6, 11.0, 11.5, 12.1],
            &simulation_time_for_storage_tank,
            1.,
        ))
    }

    #[fixture]
    pub fn storage_tank(
        simulation_time_for_storage_tank: SimulationTime,
        cold_water_source: Arc<ColdWaterSource>,
        control_for_storage_tank: Arc<Control>,
    ) -> ((StorageTank, StorageTank), Arc<Mutex<EnergySupply>>) {
        let energy_supply = Arc::new(Mutex::new(EnergySupply::new(
            FuelType::Electricity,
            simulation_time_for_storage_tank.total_steps(),
            None,
        )));
        let energy_supply_conns = (
            EnergySupply::connection(energy_supply.clone(), "immersion").unwrap(),
            EnergySupply::connection(energy_supply.clone(), "immersion2").unwrap(),
        );

        let storage_tanks = (
            StorageTank::new(
                150.0,
                1.68,
                52.0,
                55.0,
                WaterSourceWithTemperature::ColdWaterSource(cold_water_source.clone()),
                simulation_time_for_storage_tank.step,
                IndexMap::from([(
                    "imheater".to_string(),
                    PositionedHeatSource {
                        heat_source: Arc::new(Mutex::new(HeatSource::Storage(
                            HeatSourceWithStorageTank::Immersion(Arc::new(Mutex::new(
                                ImmersionHeater::new(
                                    50.,
                                    energy_supply_conns.0.clone(),
                                    simulation_time_for_storage_tank.step,
                                    Some(control_for_storage_tank.clone()),
                                ),
                            ))),
                        ))),
                        heater_position: 0.1,
                        thermostat_position: 0.33,
                    },
                )]),
                None,
                Some(energy_supply_conns.0),
                None,
                WATER.clone(),
            ),
            // Also test case where heater does not heat all layers, to ensure this is handled correctly
            StorageTank::new(
                210.0,
                1.61,
                52.0,
                60.0,
                WaterSourceWithTemperature::ColdWaterSource(cold_water_source),
                simulation_time_for_storage_tank.step,
                IndexMap::from([(
                    "imheater2".to_string(),
                    PositionedHeatSource {
                        heat_source: Arc::new(Mutex::new(HeatSource::Storage(
                            HeatSourceWithStorageTank::Immersion(Arc::new(Mutex::new(
                                ImmersionHeater::new(
                                    5.,
                                    energy_supply_conns.1.clone(),
                                    simulation_time_for_storage_tank.step,
                                    Some(control_for_storage_tank),
                                ),
                            ))),
                        ))),
                        heater_position: 0.6,
                        thermostat_position: 0.6,
                    },
                )]),
                None,
                Some(energy_supply_conns.1),
                None,
                WATER.clone(),
            ),
        );

        (storage_tanks, energy_supply)
    }

    #[rstest]
    pub fn should_calc_demand_hot_water_for_storage_tank(
        storage_tank: ((StorageTank, StorageTank), Arc<Mutex<EnergySupply>>),
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        let ((mut storage_tank1, mut storage_tank2), energy_supply) = storage_tank;
        //collect all values for steps and compare at the same time
        let mut temp_n_values_1: Vec<[f64; STORAGE_TANK_NB_VOL]> = vec![];
        let mut temp_n_values_2: Vec<[f64; STORAGE_TANK_NB_VOL]> = vec![];
        let expected_temp_n_values_1 = vec![
            [
                43.5117037037037,
                54.595555555555556,
                54.595555555555556,
                54.595555555555556,
            ],
            [
                34.923351362284535,
                51.44088940589104,
                54.19530534979424,
                54.19530534979424,
            ],
            [
                25.428671888696492,
                44.86111831060492,
                52.763271736704276,
                53.79920588690749,
            ],
            [
                17.778914378539547,
                34.731511258769736,
                48.38455458241966,
                52.883165319588585,
            ],
            [55., 55., 55., 55.],
            [
                32.955654320987655,
                54.595555555555556,
                54.595555555555556,
                54.595555555555556,
            ],
            [55., 55., 55., 55.],
            [
                33.53623703703703,
                54.595555555555556,
                54.595555555555556,
                54.595555555555556,
            ],
        ];
        let expected_energy_supply_results_1 = [
            0.0,
            0.0,
            0.0,
            0.0,
            3.9189973050595626,
            0.0,
            2.0255553251028573,
            0.0,
        ];
        let expected_temp_n_values_2 = vec![
            [
                51.74444444444445,
                59.687654320987654,
                59.687654320987654,
                59.687654320987654,
            ],
            [
                44.83576096913369,
                58.10817048730805,
                59.37752591068435,
                59.37752591068435,
            ],
            [
                36.279411505184825,
                54.60890513377094,
                58.76352191705448,
                59.06959902921961,
            ],
            [
                27.803758539213316,
                48.41088769491589,
                57.11721566595131,
                58.66493643832885,
            ],
            [
                22.115012458237494,
                41.46704433740872,
                53.98882801141131,
                57.857823384416925,
            ],
            [18.392953648519935, 34.88146733500239, 60.0, 60.0],
            [
                16.198781370486113,
                29.539425498912564,
                51.75379869179794,
                59.687654320987654,
            ],
            [14.889587258686573, 25.21241834280409, 60.0, 60.0],
        ];
        let expected_energy_supply_results_2 = [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.8689721305845337,
            0.0,
            1.1479005355748102,
        ];
        for (t_idx, iteration) in simulation_time_for_storage_tank.iter().enumerate() {
            storage_tank1.demand_hot_water(
                [10.0, 10.0, 15.0, 20.0, 20.0, 20.0, 20.0, 20.0][t_idx],
                iteration,
            );
            temp_n_values_1.push(storage_tank1.temp_n);
            assert_eq!(
                round_by_precision(
                    energy_supply.lock().results_by_end_user()["immersion"][t_idx],
                    1e7
                ),
                round_by_precision(expected_energy_supply_results_1[t_idx], 1e7),
                "incorrect energy supplied returned"
            );

            storage_tank2.demand_hot_water(
                [10.0, 10.0, 15.0, 20.0, 20.0, 20.0, 20.0, 20.0][t_idx],
                iteration,
            );
            temp_n_values_2.push(storage_tank2.temp_n);
            assert_eq!(
                round_by_precision(
                    energy_supply.lock().results_by_end_user()["immersion2"][t_idx],
                    1e7
                ),
                round_by_precision(expected_energy_supply_results_2[t_idx], 1e7),
                "incorrect energy supplied returned in case where heater does not heat all layers, iteration {t_idx}"
            );
        }
        assert_eq!(
            round_each_by_precision(temp_n_values_1, 1e7),
            round_each_by_precision(expected_temp_n_values_1, 1e7),
            "incorrect temperatures returned"
        );
        assert_eq!(
            round_each_by_precision(temp_n_values_2, 1e7),
            round_each_by_precision(expected_temp_n_values_2, 1e7),
            "incorrect temperatures returned in case where heater does not heat all layers"
        );
    }

    #[fixture]
    pub fn simulation_time_for_immersion_heater() -> SimulationTime {
        SimulationTime::new(0., 4., 1.)
    }

    #[fixture]
    pub fn immersion_heater(
        simulation_time_for_immersion_heater: SimulationTime,
    ) -> ImmersionHeater {
        ImmersionHeater::new(
            50.,
            EnergySupply::connection(
                Arc::new(Mutex::new(EnergySupply::new(
                    FuelType::MainsGas,
                    simulation_time_for_immersion_heater.total_steps(),
                    None,
                ))),
                "shower",
            )
            .unwrap(),
            simulation_time_for_immersion_heater.step,
            Some(Arc::new(Control::OnOffTimeControl(OnOffTimeControl::new(
                vec![true, true, false, true],
                0,
                1.,
            )))),
        )
    }

    #[rstest]
    pub fn should_calc_demand_energy_for_immersion_heater(
        mut immersion_heater: ImmersionHeater,
        simulation_time_for_immersion_heater: SimulationTime,
    ) {
        let energy_inputs = [40., 100., 30., 20.];
        let expected_energy = [40., 50., 0., 20.];
        for (t_idx, t_it) in simulation_time_for_immersion_heater.iter().enumerate() {
            assert_eq!(
                immersion_heater.demand_energy(energy_inputs[t_idx], t_it),
                expected_energy[t_idx],
                "incorrect energy demand calculated"
            );
        }
    }

    // following tests are from a separate test file in the Python test_storage_tank_with_solar_thermal.py

    #[fixture]
    pub fn storage_tank_with_solar_thermal() -> (
        StorageTank,
        Arc<Mutex<SolarThermalSystem>>,
        SimulationTime,
        Arc<Mutex<EnergySupply>>,
    ) {
        let cold_water_temps = [
            17.0, 17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7, 17.0, 17.1, 17.2, 17.3, 17.4, 17.5,
            17.6, 17.7, 17.0, 17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7,
        ];
        let simulation_time = SimulationTime::new(5088., 5112., 1.);
        let cold_feed = WaterSourceWithTemperature::ColdWaterSource(Arc::new(
            ColdWaterSource::new(cold_water_temps.to_vec(), &simulation_time, 1.),
        ));
        let energy_supply = Arc::new(Mutex::new(EnergySupply::new(
            FuelType::Electricity,
            simulation_time.total_steps(),
            None,
        )));
        let energy_supply_conn =
            EnergySupply::connection(energy_supply.clone(), "solarthermal").unwrap();

        let external_conditions = Arc::new(ExternalConditions::new(
            &simulation_time.clone().iter(),
            vec![
                19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0,
                19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0,
            ],
            vec![
                3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1, 3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1,
                3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1,
            ],
            vec![
                0, 0, 0, 0, 35, 73, 139, 244, 320, 361, 369, 348, 318, 249, 225, 198, 121, 68, 19,
                0, 0, 0, 0, 0,
            ]
            .iter()
            .map(|x| *x as f64)
            .collect(),
            vec![
                0, 0, 0, 0, 0, 0, 7, 53, 63, 164, 339, 242, 315, 577, 385, 285, 332, 126, 7, 0, 0,
                0, 0, 0,
            ]
            .iter()
            .map(|x| *x as f64)
            .collect(),
            vec![
                0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
                0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
            ],
            51.383,
            -0.783,
            0,
            212,
            Some(212),
            1.,
            Some(1),
            DaylightSavingsConfig::NotApplicable,
            false,
            false,
            vec![
                ShadingSegment {
                    number: 1,
                    start: 180,
                    end: 135,
                    objects: None,
                },
                ShadingSegment {
                    number: 2,
                    start: 135,
                    end: 90,
                    objects: None,
                },
                ShadingSegment {
                    number: 3,
                    start: 90,
                    end: 45,
                    objects: None,
                },
                ShadingSegment {
                    number: 4,
                    start: 45,
                    end: 0,
                    objects: Some(vec![ShadingObject {
                        object_type: ShadingObjectType::Obstacle,
                        height: 10.5,
                        distance: 120.,
                    }]),
                },
                ShadingSegment {
                    number: 5,
                    start: 0,
                    end: -45,
                    objects: None,
                },
                ShadingSegment {
                    number: 6,
                    start: -45,
                    end: -90,
                    objects: None,
                },
                ShadingSegment {
                    number: 7,
                    start: -90,
                    end: -135,
                    objects: None,
                },
                ShadingSegment {
                    number: 8,
                    start: -135,
                    end: -180,
                    objects: None,
                },
            ],
        ));
        let solar_thermal = Arc::new(Mutex::new(SolarThermalSystem::new(
            SolarCellLocation::Out,
            3.,
            1,
            0.8,
            0.9,
            3.5,
            0.,
            1.,
            100.,
            10.,
            energy_supply_conn,
            30.,
            0.,
            0.5,
            external_conditions,
            simulation_time.step,
            WATER.clone(),
        )));

        let storage_tank = StorageTank::new(
            150.0,
            1.68,
            52.0,
            55.0,
            cold_feed,
            simulation_time.step,
            IndexMap::from([(
                "solthermal".to_string(),
                PositionedHeatSource {
                    heat_source: Arc::new(Mutex::new(HeatSource::Storage(
                        HeatSourceWithStorageTank::Solar(solar_thermal.clone()),
                    ))),
                    heater_position: 0.1,
                    thermostat_position: 0.33,
                },
            )]),
            None,
            None,
            None,
            WATER.clone(),
        );

        (storage_tank, solar_thermal, simulation_time, energy_supply)
    }

    #[rstest]
    pub fn test_demand_hot_water(
        storage_tank_with_solar_thermal: (
            StorageTank,
            Arc<Mutex<SolarThermalSystem>>,
            SimulationTime,
            Arc<Mutex<EnergySupply>>,
        ),
    ) {
        let (mut storage_tank, solar_thermal, simulation_time, energy_supply) =
            storage_tank_with_solar_thermal;
        let demands = [
            100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let expected_energy_demands = [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.3944013153635651,
            0.7205382866008986,
            1.1792815529120688,
            0.9563670953583516,
            1.066201484260018,
            0.2842009512268733,
            0.07050814814814632,
            0.07050814814814682,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ];
        let expected_energy_potentials = [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.3944013153635651,
            0.7205382866008986,
            1.1792815529120688,
            0.9563670953583516,
            1.066201484260018,
            1.3754941274949404,
            0.788682346923819,
            0.4490991945005249,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ];
        let expected_energy_supplied = [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.3944013153635651,
            0.7205382866008986,
            1.1792815529120688,
            0.9563670953583516,
            1.066201484260018,
            0.2842009512268733,
            0.07050814814814632,
            0.07050814814814682,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ];
        let expected_energy_supply_results = [
            10, 10, 10, 10, 10, 10, 10, 10, 110, 110, 110, 110, 110, 110, 110, 110, 10, 10, 10, 10,
            10, 10, 10, 10,
        ]
        .map(|x| x as f64);

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            storage_tank.demand_hot_water(demands[t_idx], t_it);
            assert_eq!(
                round_by_precision(storage_tank.test_energy_demand(), 1e7),
                round_by_precision(expected_energy_demands[t_idx], 1e7),
                "incorrect energy demand from tank"
            );
            assert_eq!(
                round_by_precision(solar_thermal.lock().test_energy_potential(), 1e7),
                round_by_precision(expected_energy_potentials[t_idx], 1e7),
                "incorrect energy potential by solar thermal returned"
            );
            assert_eq!(
                round_by_precision(solar_thermal.lock().test_energy_supplied(), 1e7),
                round_by_precision(expected_energy_supplied[t_idx], 1e7),
                "incorrect energy supplied by solar thermal returned"
            );
            assert_eq!(
                round_by_precision(
                    energy_supply.lock().results_by_end_user()["solarthermal"][t_idx],
                    1e7
                ),
                round_by_precision(expected_energy_supply_results[t_idx], 1e7),
                "incorrect electric energy consumed returned"
            );
        }
    }
}
