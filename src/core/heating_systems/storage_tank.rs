use crate::compare_floats::{max_of_2, min_of_2};
use crate::core::common::WaterSourceWithTemperature;
use crate::core::controls::time_control::Control;
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::material_properties::MaterialProperties;
use crate::core::pipework::{Pipework, PipeworkLocation, Pipeworkesque};
use crate::core::schedule::TypedScheduleEvent;
use crate::core::units::WATTS_PER_KILOWATT;
use crate::core::water_heat_demand::misc::frac_hot_water;
use crate::corpus::{HeatSource, PositionedHeatSource, TempInternalAirAccessor};
use crate::external_conditions::ExternalConditions;
use crate::input::{SolarCellLocation, WaterPipework};
use crate::simulation_time::SimulationTimeIteration;
use atomic_float::AtomicF64;
use indexmap::IndexMap;
use itertools::Itertools;
use parking_lot::{Mutex, RwLock};
use std::collections::HashMap;
use std::iter;
use std::sync::atomic::Ordering;
use std::sync::Arc;

// BS EN 15316-5:2017 Appendix B default input data
// Model Information
// Product Description Data
// factors for energy recovery Table B.3
// part of the auxiliary energy transmitted to the medium
const STORAGE_TANK_F_RVD_AUX: f64 = 0.25;

// part of the thermal losses transmitted to the room. Note same approach in BufferTank if this is
// modified in future
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
    energy_supply_conn_unmet_demand: Option<EnergySupplyConnection>,
    control_hold_at_setpoint: Option<Arc<Control>>,
    nb_vol: usize,
    temp_internal_air_accessor: TempInternalAirAccessor,
    external_conditions: Arc<ExternalConditions>,
    volume_total_in_litres: f64,
    vol_n: Vec<f64>,
    cp: f64,  // contents (usually water) specific heat in kWh/kg.K
    rho: f64, // volumic mass in kg/litre
    temp_n: Vec<f64>,
    input_energy_adj_prev_timestep: f64,
    primary_pipework_lst: Option<Vec<Pipework>>,
    primary_pipework_losses_kwh: f64,
    storage_losses_kwh: f64,
    heat_source_data: IndexMap<String, PositionedHeatSource>, // heat sources, sorted by heater position
    heating_active: HashMap<String, bool>,
    q_ls_n_prev_heat_source: Vec<f64>,
    q_sto_h_ls_rbl: Option<f64>, // total recoverable heat losses for heating in kWh, memoised between steps
    primary_gains: f64,          // primary pipework gains for a timestep (mutates over lifetime)
    #[cfg(test)]
    energy_demand_test: f64,
    temp_final_drawoff: Option<f64>, // In Python this is created from inside allocate_hot_water()
    temp_average_drawoff: Option<f64>, // In Python this is created from inside allocate_hot_water()
    temp_average_drawoff_volweighted: Option<f64>, // In Python this is created from inside allocate_hot_water()
    total_volume_drawoff: Option<f64>, // In Python this is created from inside allocate_hot_water()
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
    /// *  `nb_vol` -number of volumes the storage is modelled with
    ///              see App.C (C.1.2 selection of the number of volumes to model the storage unit)
    ///              for more details if this wants to be changed.
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
        // In Python this is "project" but only temp_internal_air is accessed from it
        temp_internal_air_accessor: TempInternalAirAccessor,
        external_conditions: Arc<ExternalConditions>,
        nb_vol: Option<usize>,
        primary_pipework_lst: Option<&Vec<WaterPipework>>,
        energy_supply_conn_unmet_demand: Option<EnergySupplyConnection>,
        control_hold_at_setpoint: Option<Arc<Control>>,
        contents: MaterialProperties,
    ) -> Self {
        let q_std_ls_ref = losses;
        let temp_out_w_min = min_temp;
        let temp_set_on = setpoint_temp;

        let volume_total_in_litres = volume;
        let nb_vol = nb_vol.unwrap_or(4);
        // list of volume of layers in litres
        let vol_n = iter::repeat(volume_total_in_litres / nb_vol as f64)
            .take(nb_vol)
            .collect_vec();
        // water specific heat in kWh/kg.K
        let cp = contents.specific_heat_capacity_kwh();
        let rho = contents.density();

        // 6.4.3.2 STEP 0 Initialization
        //   for initial conditions all temperatures in the thermal storage unit(s)
        //   are equal to the set point temperature in degrees.
        //   We are expecting to run a "warm-up" period for the main calculation so this doesn't matter.
        let temp_n = iter::repeat(temp_set_on).take(nb_vol).collect_vec();

        #[cfg(test)]
        let energy_demand_test = 0.;

        // primary_pipework_losses_kwh added for reporting
        let primary_pipework_losses_kwh = 0.;
        let storage_losses_kwh = 0.;

        let input_energy_adj_prev_timestep = 0.;

        let primary_pipework_lst: Option<Vec<Pipework>> = if primary_pipework_lst.is_some() {
            primary_pipework_lst
                .unwrap()
                .into_iter()
                .map(|pipework| pipework.to_owned().try_into().unwrap())
                .collect::<Vec<Pipework>>()
                .into()
        } else {
            None
        };

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
            energy_supply_conn_unmet_demand,
            control_hold_at_setpoint,
            nb_vol,
            temp_internal_air_accessor,
            external_conditions,
            volume_total_in_litres,
            vol_n,
            cp,
            rho,
            temp_n,
            input_energy_adj_prev_timestep,
            primary_pipework_lst,
            primary_pipework_losses_kwh,
            storage_losses_kwh,
            heat_source_data: heat_sources,
            heating_active,
            q_ls_n_prev_heat_source: Default::default(),
            q_sto_h_ls_rbl: Default::default(),
            primary_gains: Default::default(),
            #[cfg(test)]
            energy_demand_test,
            temp_final_drawoff: None,
            temp_average_drawoff: None,
            temp_average_drawoff_volweighted: None,
            total_volume_drawoff: None,
        }
    }

    fn temp_surrounding_primary_pipework(
        &self,
        pipework_data: &Pipework,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        match pipework_data.location() {
            PipeworkLocation::External => self
                .external_conditions
                .air_temp(&simulation_time_iteration),
            PipeworkLocation::Internal => self.temp_internal_air_accessor.call(),
        }
    }

    pub fn get_cold_water_source(&self) -> &WaterSourceWithTemperature {
        &self.cold_feed
    }

    /// This is only used to calculate the equivalent volume of water for IES showers
    /// in order to get the energy content for the internal gains.
    /// Therefore the actual value used is not critical.
    /// It has been suggested/considered the use of the top layer of the storage tank
    /// but this could be similar to the cold feed temperature after big draw-offs
    /// To avoid any issues in those situations we use the setpoing temperature of the
    /// tank.
    fn get_temp_hot_water(&self) -> f64 {
        self.temp_set_on
    }

    /// Return temp_out_w_min unless tank is being held at setpnt, in which case return that
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
    /// h_sto_ls is the stand-by losses, in W/K
    ///
    /// TODO (from Python) there are alternative methods listed in App B (B.2.8) which are not included here.
    pub fn stand_by_losses_coefficient(&self) -> f64 {
        // BS EN 12897:2016 appendix B B.2.2
        // temperature of the water in the storage for the standardized conditions - degrees
        // these are reference (ref) temperatures from the standard test conditions for cylinder loss.
        let temp_set_ref = 65.;
        let temp_amb_ref = 20.;

        (1000. * self.q_std_ls_ref) / (24. * (temp_set_ref - temp_amb_ref))
    }

    /// Energy input for the storage from the generation system
    /// (expressed per energy carrier X)
    /// Heat Source = energy carrier
    pub fn potential_energy_input(
        // Heat source. Addition of temp_s3_n as an argument
        &mut self,
        temp_s3_n: &[f64],
        heat_source: &Arc<Mutex<HeatSource>>,
        heat_source_name: &str,
        heater_layer: usize,
        thermostat_layer: usize,
        simulation_time: SimulationTimeIteration,
    ) -> Vec<f64> {
        // initialise list of potential energy input for each layer
        let mut q_x_in_n = iter::repeat(0.).take(self.nb_vol).collect_vec();

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
                    let mut energy_potential = immersion_heater.lock().energy_output_max(
                        simulation_time,
                        false,
                    );
                    // TODO (from Python) Consolidate checks for systems with/without primary pipework

                    if !matches!(
                        heat_source,
                        HeatSource::Storage(HeatSourceWithStorageTank::Immersion(_))
                    ) {
                        let (primary_pipework_losses_kwh, _) =
                            self.primary_pipework_losses(energy_potential, simulation_time);
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
        temp_s3_n: &[f64],
    ) -> (f64, f64) {
        // assuming initially no water draw-off

        // initialise list of potential energy input for each layer
        let mut q_x_in_n = vec![0.; self.nb_vol];

        // TODO (from Python) - ensure we are feeding in the correct volume
        q_x_in_n[0] = energy_proposed;

        let (_q_s6, temp_s6_n) = self.energy_input(&temp_s3_n, &q_x_in_n);

        // 6.4.3.9 STEP 7 Re-arrange the temperatures in the storage after energy input
        let (_q_h_sto_s7, temp_s7_n) = self.rearrange_temperatures(&temp_s6_n);

        // TODO (from Python) Check [0] is bottom layer temp and that solar thermal inlet is top layer NB_VOL-1
        (temp_s7_n[0], temp_s7_n[self.nb_vol - 1])
    }

    /// The input of energy(s) is (are) allocated to the specific location(s)
    /// of the input of energy.
    /// Note: for energy withdrawn from a heat exchanger, the energy is accounted negatively.
    ///
    /// For step 6, the addition of the temperature of volume 'i' and theoretical variation of
    /// temperature calculated according to formula (10) can exceed the set temperature defined
    /// by the control system of the storage unit.
    fn energy_input(&self, temp_s3_n: &[f64], q_x_in_n: &[f64]) -> (f64, Vec<f64>) {
        // initialise list of theoretical variation of temperature of layers in degrees
        let mut delta_temp_n = vec![0.; self.nb_vol];
        // initialise list of theoretical temperature of layers after input in degrees
        let mut temp_s6_n = vec![0.; self.nb_vol];
        // output energy delivered by the storage in kWh - timestep dependent
        let q_sto_h_out_n: Vec<f64> = vec![0.; self.nb_vol];

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
    /// then the 2 volumes are melded. This iterative process is maintained until the temperature
    /// of the volume i is lower or equal to the temperature of the volume i+1.
    fn rearrange_temperatures(&self, temp_s6_n: &[f64]) -> (Vec<f64>, Vec<f64>) {
        // set list of flags for which layers need mixing
        let mut mix_layer_n: Vec<u8> = vec![0; self.nb_vol];
        let mut temp_s7_n = temp_s6_n.to_vec();

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
                mix_layer_n = vec![0; self.nb_vol];
            }
        }

        let q_h_sto_end = (0..self.vol_n.len())
            .map(|i| self.rho * self.cp * self.vol_n[i] * temp_s7_n[i])
            .collect::<Vec<f64>>()
            .try_into()
            .unwrap();

        let _q_h_sto_end_no = self.rho
            * self.cp
            * (0..self.vol_n.len())
                .map(|i| self.vol_n[i] * temp_s7_n[i])
                .sum::<f64>();

        (q_h_sto_end, temp_s7_n.to_owned())
    }

    /// Thermal losses are calculated with respect to the impact of the temperature set point
    pub fn thermal_losses(
        &mut self,
        temp_s7_n: &[f64],
        q_x_in_n: &[f64],
        q_h_sto_s7: Vec<f64>,
        heater_layer: usize,
        q_ls_n_prev_heat_source: &[f64],
    ) -> (f64, f64, Vec<f64>, Vec<f64>) {
        // standby losses coefficient - W/K
        let h_sto_ls = self.stand_by_losses_coefficient();

        // standby losses correction factor - dimensionless
        // note from Python code: "do not think these are applicable so used: f_sto_dis_ls = 1, f_sto_bac_acc = 1"

        // initialise list of thermal losses in kWh
        let mut q_ls_n: Vec<f64> = vec![];
        // initialise list of final temperature of layers after thermal losses in degrees
        let mut temp_s8_n: Vec<f64> = vec![];

        // Thermal losses
        // Note: Eqn 13 from BS EN 15316-5:2017 does not explicitly multiply by
        // timestep (it seems to assume a 1 hour timestep implicitly), but it is
        // necessary to convert the rate of heat loss to a total heat loss over
        // the time period
        for i in 0..self.vol_n.len() {
            let q_ls_n_step = (h_sto_ls * self.rho * self.cp)
                * (self.vol_n[i] / self.volume_total_in_litres)
                * (min_of_2(temp_s7_n[i], self.temp_set_on) - STORAGE_TANK_TEMP_AMB)
                * self.simulation_timestep;

            let q_ls_n_step = max_of_2(0., q_ls_n_step - q_ls_n_prev_heat_source[i]);

            q_ls_n.push(q_ls_n_step);
        }

        // total thermal losses kWh
        let q_ls = q_ls_n.iter().sum();

        self.storage_losses_kwh = q_ls;

        // the final value of the temperature is reduced due to the effect of the thermal losses.
        // check temperature compared to set point
        // the temperature for each volume are limited to the set point for any volume controlled
        for i in 0..self.vol_n.len() {
            let temp_s8_n_step = if temp_s7_n[i] > self.temp_set_on {
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
            temp_s8_n.push(temp_s8_n_step);
        }

        // excess energy/ energy surplus
        // excess energy is calculated as the difference from the energy stored, Qsto,step7, and
        // energy stored once the set temperature is obtained, Qsto,step8, with addition of the
        // thermal losses.
        // Note: The surplus must be calculated only for those layers that the
        //       heat source currently being considered is capable of heating,
        //       i.e. excluding those below the heater position.
        let energy_surplus = if temp_s7_n[heater_layer] > self.temp_set_on {
            (heater_layer..self.nb_vol).fold(0., |acc, i| {
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
        temp_s3_n: Vec<f64>,
        heat_source: Arc<Mutex<HeatSource>>,
        heat_source_name: &str,
        heater_layer: usize,
        thermostat_layer: usize,
        q_ls_prev_heat_source: &[f64],
        simulation_time: SimulationTimeIteration,
    ) -> TemperatureCalculation {
        // 6.4.3.8 STEP 6 Energy input into the storage
        // input energy delivered to the storage in kWh - timestep dependent
        let q_x_in_n = self.potential_energy_input(
            &temp_s3_n,
            &heat_source,
            heat_source_name,
            heater_layer,
            thermostat_layer,
            simulation_time,
        );

        self.calculate_temperatures(
            &temp_s3_n,
            heat_source,
            q_x_in_n,
            heater_layer,
            q_ls_prev_heat_source,
            simulation_time,
        )
    }

    fn calculate_temperatures(
        &mut self,
        temp_s3_n: &[f64],
        heat_source: Arc<Mutex<HeatSource>>,
        q_x_in_n: Vec<f64>,
        heater_layer: usize,
        q_ls_n_prev_heat_source: &[f64],
        simulation_time_iteration: SimulationTimeIteration,
    ) -> TemperatureCalculation {
        let (q_s6, temp_s6_n) = self.energy_input(&temp_s3_n, &q_x_in_n);

        // 6.4.3.9 STEP 7 Re-arrange the temperatures in the storage after energy input
        let (q_h_sto_s7, temp_s7_n) = self.rearrange_temperatures(&temp_s6_n);

        // STEP 8 Thermal losses and final temperature
        let (q_in_h_w, q_ls, temp_s8_n, q_ls_n) = self.thermal_losses(
            &temp_s7_n,
            &q_x_in_n,
            q_h_sto_s7,
            heater_layer,
            q_ls_n_prev_heat_source,
        );

        // TODO (from Python) 6.4.3.11 Heat exchanger

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

    /// Allocate hot water layers to meet a single temperature demand.
    ///
    /// Arguments:
    /// * `event` -- Dictionary containing information about the draw-off event
    /// (e.g. {'start': 18, 'duration': 1, 'temperature': 41.0, 'type': 'Other', 'name': 'other', 'warm_volume': 8.0})
    fn allocate_hot_water(
        &mut self,
        event: TypedScheduleEvent,
        simulation_time: SimulationTimeIteration,
    ) -> (f64, Vec<f64>, f64, f64) {
        // Make a copy of the volume list to keep track of remaining volumes
        // Remaining volume of water in storage tank layers
        let mut remaining_vols = self.vol_n.clone();

        // Extract the temperature and required warm volume from the event
        let mut warm_temp = event.temperature;
        let warm_volume = event.warm_volume;
        // Remaining volume of warm water to be satisfied for current event
        let mut remaining_demanded_warm_volume = warm_volume.unwrap();

        // Initialize the unmet and met energies
        let mut energy_unmet = 0.0;
        let mut energy_withdrawn = 0.0;
        let mut pipework_temp = self.cold_feed.temperature(simulation_time.index); // This value set to initialise, but is never used - overwritten later.

        let mut pipework_considered = if event.pipework_volume.unwrap() <= 0.0 {
            true
        } else {
            false
        };

        let mut temp_average_drawoff_volweighted: f64 = Default::default();
        let mut total_volume_drawoff: f64 = Default::default();
        let mut last_layer_index: usize = Default::default();
        //  Loop through storage layers (starting from the top)
        for layer_index in (0..self.temp_n.len()).rev() {
            last_layer_index = layer_index;
            let layer_temp = self.temp_n[layer_index];
            let layer_vol = remaining_vols[layer_index];

            if remaining_demanded_warm_volume <= 0. {
                if pipework_considered {
                    // Event inclusive of pipework is completed at this layer temp
                    self.temp_final_drawoff = Some(pipework_temp);
                    break;
                } else {
                    remaining_demanded_warm_volume = event.pipework_volume.unwrap();
                    warm_temp = layer_temp;
                    pipework_considered = true;
                }
            }
            // If event is finished and we are serving the pipework, this is the temperature
            // of the water stranded
            pipework_temp = layer_temp;

            // Skip this layer if its remaining volume is already zero
            if remaining_vols[layer_index] <= 0.0 {
                continue;
            }

            // Skip this layer if its temperature is lower than the target temperature
            if layer_temp < warm_temp {
                break;
            }

            // Calculate the fraction of hot water required
            let fraction = frac_hot_water(
                warm_temp,
                layer_temp,
                self.cold_feed.temperature(simulation_time.index),
            );

            let _warm_vol_removed: f64;
            let required_vol: f64;
            // Volume of hot water required at this layer
            if layer_vol <= remaining_demanded_warm_volume * fraction {
                // This is the case where layer cannot meet all remaining demand for this event
                required_vol = layer_vol;
                _warm_vol_removed = layer_vol / fraction;
                // Deduct the required volume from the remaining demand and update the layer's volume
                remaining_vols[layer_index] -= layer_vol;
                remaining_demanded_warm_volume -= _warm_vol_removed;
            } else {
                //This is the case where layer can meet all remaining demand for this event
                required_vol = remaining_demanded_warm_volume * fraction;
                _warm_vol_removed = remaining_demanded_warm_volume;
                // Deduct the required volume from the remaining demand and update the layer's volume
                remaining_vols[layer_index] -= required_vol;
                remaining_demanded_warm_volume = 0.0;
            }

            temp_average_drawoff_volweighted += required_vol * layer_temp;
            total_volume_drawoff += required_vol;

            // Record the met volume demand for the current temperature target
            // warm_vol_removed is the volume of warm water that has been satisfied from hot water in this layer
            energy_withdrawn +=
                // Calculation with event water parameters
                // self.__rho * self.__Cp * warm_vol_removed * (warm_temp - self.__cold_feed.temperature())
                // Calculation with layer water parameters
                self.rho
                    * self.cp
                    * required_vol
                    * (layer_temp - self.cold_feed.temperature(simulation_time.index))
        }

        self.temp_average_drawoff_volweighted = Some(temp_average_drawoff_volweighted);
        self.total_volume_drawoff = Some(total_volume_drawoff);

        // When the event has not been fully met or has been exactly met with the last of the hot water
        // in the tank, there's only cold water from the feed left to fill the pipework after the event.
        if !pipework_considered {
            self.temp_final_drawoff = Some(self.temp_n[last_layer_index]);
        }

        // Record the unmet energy for the current event
        energy_unmet += self.rho
            * self.cp
            * remaining_demanded_warm_volume
            * (warm_temp - self.cold_feed.temperature(simulation_time.index));

        //  Calculate the remaining total volume
        let remaining_total_volume: f64 = remaining_vols.iter().sum();

        //  Calculate the total volume used
        let volume_used = self.volume_total_in_litres - remaining_total_volume;

        //  Determine the new temperature distribution after displacement
        let new_temp_distribution =
            self.calculate_new_temperatures(remaining_vols, simulation_time);

        //  Return the remaining storage volumes, volume used, new temperature distribution, and the met/unmet targets
        (
            volume_used,
            new_temp_distribution,
            energy_withdrawn,
            energy_unmet,
        )
    }

    /// Calculate the new temperature distribution after displacement.
    /// Arguments:
    /// * `remaining_vols` -- List of remaining volumes for each storage layer after draw-off
    /// * `temp_cold` -- Temperature of the cold water being added
    fn calculate_new_temperatures(
        &self,
        mut remaining_vols: Vec<f64>,
        simulation_time: SimulationTimeIteration,
    ) -> Vec<f64> {
        let mut new_temps = self.temp_n.clone();

        // Iterate from the top layer downwards
        for i in (0..self.vol_n.len()).rev() {
            // Determine how much volume needs to be added to this layer
            let mut needed_volume = self.vol_n[i] - remaining_vols[i];
            // If this layer is already full, continue to the next
            if needed_volume <= 0. {
                break;
            }

            // Initialize the variables for mixing temperatures
            let mut total_volume = remaining_vols[i];
            let mut volume_weighted_temperature = remaining_vols[i] * self.temp_n[i];

            // Add water from the layers below to this layer
            for j in (0..i).rev() {
                let available_volume = remaining_vols[j];
                if available_volume > 0. {
                    // Determine the volume to move up from this layer
                    let move_volume = f64::min(needed_volume, available_volume);
                    remaining_vols[j] -= move_volume;

                    // Adjust the temperature by mixing in the moved volume
                    total_volume += move_volume;
                    volume_weighted_temperature += move_volume * self.temp_n[j];

                    // Decrease the amount of volume needed for the current layer
                    needed_volume -= move_volume;
                    if needed_volume <= 0. {
                        break;
                    }
                }
            }

            // If not enough water is available from the lower layers, use the cold supply
            if needed_volume > 0. {
                total_volume += needed_volume;
                volume_weighted_temperature +=
                    needed_volume * self.cold_feed.temperature(simulation_time.index);
            }

            // Calculate the new temperature for the current layer
            // Round to 2 decimals to match instrumentation limits and significant figures,
            // ensuring practical accuracy, computational efficiency, and avoiding minute e-18 differences.
            new_temps[i] = (100. * (volume_weighted_temperature / total_volume)).round() / 100.;
            remaining_vols[i] = total_volume;
        }
        new_temps
    }

    /// Draw off hot water from the tank
    /// Energy calculation as per BS EN 15316-5:2017 Method A sections 6.4.3, 6.4.6, 6.4.7
    /// Modification of calculation based on volumes and actual temperatures for each layer of water in the tank
    /// instead of the energy stored in the layer and a generic temperature (self.temp_out_w_min) = min_temp
    /// to decide if the tank can satisfy the demand (this was producing unnecesary unmet demand for strict high
    /// temp_out_w_min values
    /// Arguments:
    /// * `usage_events` -- All draw off events for the timestep
    pub fn demand_hot_water(
        &mut self,
        usage_events: Vec<TypedScheduleEvent>,
        simulation_time: SimulationTimeIteration,
    ) -> (f64, f64, f64, f64, f64) {
        let mut q_use_w = 0.;
        let mut q_unmet_w = 0.;
        let mut _volume_demanded = 0.;

        let mut temp_s3_n = self.temp_n.clone();

        self.temp_final_drawoff = Some(self.get_temp_hot_water());
        self.temp_average_drawoff_volweighted = Some(0.);
        self.total_volume_drawoff = Some(0.);
        self.temp_average_drawoff = Some(self.get_temp_hot_water());

        // Filtering out IES events that don't get added a 'warm_volume' when processing
        // the dhw_demand calculation
        let filtered_events = usage_events
            .into_iter()
            .filter(|e| e.warm_volume.is_some())
            .collect_vec();

        for mut event in filtered_events {
            // Check if 'pipework_volume' key exists in the event dictionary
            if event.pipework_volume.is_none() {
                // If 'pipework_volume' is not found, add it with a default value of 0.0
                event.pipework_volume = Some(0.0);
            }

            // Decision no to include yet the overlapping of events for pipework losses
            // even if applying pipework losses to all events might be overstimating
            // the following overlapping processing could be understimating for multiple
            // branches of the pipework system
            // TODO (from Python) Improve approach for avoiding double counting of genuine overlapping
            // events
            // Avoid double counting pipework loses when events overlap
            // time_start_current_event = event['start']
            // if self.__time_end_previous_event >= time_start_current_event:
            // event['pipework_volume'] = 0.0
            // 0.0 can be modified for additional minutes when pipework could be considered still warm/hot
            // self.__time_end_previous_event = deepcopy(time_start_current_event + (event['duration'] + 0.0) / 60.0)

            let (volume_used, temp_s3_n_step, energy_withdrawn, energy_unmet) =
                self.allocate_hot_water(event.clone(), simulation_time);

            temp_s3_n = temp_s3_n_step;
            self.temp_n = temp_s3_n.clone();

            _volume_demanded += volume_used;
            q_unmet_w += energy_unmet;
            q_use_w += energy_withdrawn;
        }

        // if tank cannot provide enough hot water report unmet demand
        if self.energy_supply_conn_unmet_demand.is_some() {
            self.energy_supply_conn_unmet_demand
                .as_ref()
                .unwrap()
                .demand_energy(q_unmet_w, simulation_time.index)
                .expect("expected to be able to demand energy");
        }

        self.temp_average_drawoff = match self.total_volume_drawoff {
            Some(value) if value != 0. => {
                let temp_average_drawoff_volweighted = self
                    .temp_average_drawoff_volweighted
                    .expect("temp_average_drawoff_volweighted was not set");
                Some(temp_average_drawoff_volweighted / value)
            },
            _ => {
                temp_s3_n.last().copied()
            }
        };
        
        // TODO (from Python) 6.4.3.6 STEP 4 Volume to be withdrawn from the storage (for Heating)
        // TODO (from Python) - 6.4.3.7 STEP 5 Temperature of the storage after volume withdrawn (for Heating)

        // Run over multiple heat sources
        let mut temp_after_prev_heat_source = temp_s3_n;
        let mut q_ls = 0.0;
        let mut temp_s8_n = vec![0.; self.nb_vol];

        self.q_ls_n_prev_heat_source = vec![0.0; self.nb_vol];
        for (heat_source_name, positioned_heat_source) in self.heat_source_data.clone() {
            let heater_layer =
                (positioned_heat_source.heater_position * self.nb_vol as f64) as usize;
            let thermostat_layer =
                (positioned_heat_source.thermostat_position * self.nb_vol as f64) as usize;

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
                temp_after_prev_heat_source.clone(),
                positioned_heat_source.heat_source.clone(),
                &heat_source_name,
                heater_layer,
                thermostat_layer,
                &self.q_ls_n_prev_heat_source.clone(),
                simulation_time,
            );

            temp_after_prev_heat_source = temp_s8_n_step.clone();
            q_ls += q_ls_this_heat_source;

            for (i, q_ls_n) in q_ls_n_this_heat_source.iter().enumerate() {
                self.q_ls_n_prev_heat_source[i] += q_ls_n;
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

        // TODO (from Python) recoverable heat losses for heating should impact heating

        // Return total energy of hot water supplied and unmet
        (
            q_use_w,
            q_unmet_w,
            self.temp_final_drawoff.unwrap(),
            self.temp_average_drawoff.unwrap(),
            self.total_volume_drawoff.unwrap(),
        )
        // Sending temp_final_drawoff, temp_average_drawoff
        // for pipework loss and internal gains calculations
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

        let heater_layer = (heat_source_data.heater_position * self.nb_vol as f64) as usize;

        let mut q_x_in_n = vec![0.; self.nb_vol];
        q_x_in_n[heater_layer] = energy_input;
        let (temp_s8_n, _, _, _, _, q_in_h_w, _, q_ls_n_this_heat_source) = self
            .calculate_temperatures(
                &self.temp_n.clone(),
                heat_source,
                q_x_in_n,
                heater_layer,
                &self.q_ls_n_prev_heat_source.clone(),
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

    fn primary_pipework_losses(
        &mut self,
        input_energy_adj: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (f64, f64) {
        let mut primary_pipework_losses_kwh = Default::default();
        let mut primary_gains_w = Default::default();
        let primary_pipework_lst = self.primary_pipework_lst.as_ref().expect(
            "primary pipeworks are expected to have been set on the storage tank at this point",
        );

        // start of heating event
        if input_energy_adj > 0. && self.input_energy_adj_prev_timestep == 0. {
            for pipework_data in primary_pipework_lst {
                primary_pipework_losses_kwh += pipework_data.cool_down_loss(
                    self.temp_set_on,
                    self.temp_surrounding_primary_pipework(
                        pipework_data,
                        simulation_time_iteration,
                    ),
                )
            }
        }

        // during heating event
        if input_energy_adj > 0. {
            for pipework_data in primary_pipework_lst {
                // Primary losses for the timestep calculated from temperature difference
                let primary_pipework_losses_w = pipework_data.heat_loss(
                    self.temp_set_on,
                    self.temp_surrounding_primary_pipework(
                        pipework_data,
                        simulation_time_iteration,
                    ),
                );

                // Check if pipework location is internal
                let location = pipework_data.location();
                match location {
                    PipeworkLocation::Internal => primary_gains_w += primary_pipework_losses_w,
                    PipeworkLocation::External => {
                        primary_pipework_losses_kwh += primary_pipework_losses_w
                            * self.simulation_timestep
                            / WATTS_PER_KILOWATT as f64
                    }
                }
            }
        }

        // end of heating event
        if input_energy_adj == 0. && self.input_energy_adj_prev_timestep > 0. {
            for pipework_data in primary_pipework_lst {
                let location = pipework_data.location();
                match location {
                    PipeworkLocation::External => {}
                    PipeworkLocation::Internal => {
                        primary_gains_w +=
                            pipework_data.cool_down_loss(self.temp_set_on, STORAGE_TANK_TEMP_AMB)
                    }
                }
            }
        }

        // keeping primary_pipework_losses_kWh for reporting as part of investigation of issue #31225: FDEV A082
        self.primary_pipework_losses_kwh = primary_pipework_losses_kwh;

        (primary_pipework_losses_kwh, primary_gains_w)
    }

    fn heat_source_output(
        &mut self,
        heat_source: Arc<Mutex<HeatSource>>,
        input_energy_adj: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        match &mut *(heat_source.clone().lock()) {
            HeatSource::Storage(HeatSourceWithStorageTank::Immersion(immersion)) => Ok(immersion
                .lock()
                .demand_energy(input_energy_adj, simulation_time_iteration)),
            HeatSource::Storage(HeatSourceWithStorageTank::Solar(solar)) => Ok(solar
                .lock()
                .demand_energy(input_energy_adj, simulation_time_iteration.index)),
            _ => {
                // TODO need to be able to call demand_energy on the other heat sources
                let (primary_pipework_losses_kwh, primary_gains) =
                    self.primary_pipework_losses(input_energy_adj, simulation_time_iteration);
                let input_energy_adj = input_energy_adj + primary_pipework_losses_kwh;

                let heat_source_output = heat_source.lock().demand_energy(
                    input_energy_adj,
                    self.temp_set_on,
                    simulation_time_iteration,
                )? - primary_pipework_losses_kwh;
                self.input_energy_adj_prev_timestep = input_energy_adj;
                self.primary_gains = primary_gains;

                // TODO (from Python) - how are these gains reflected in the calculations? allocation by zone?
                Ok(heat_source_output)
            }
        }
    }
}

type TemperatureCalculation = (
    Vec<f64>,
    Vec<f64>,
    f64,
    Vec<f64>,
    Vec<f64>,
    f64,
    f64,
    Vec<f64>,
);

#[derive(Clone, Debug)]
pub struct ImmersionHeater {
    pwr: f64, // rated power
    energy_supply_connection: EnergySupplyConnection,
    simulation_timestep: f64,
    control: Option<Arc<Control>>,
    diverter: Option<Arc<RwLock<PVDiverter>>>,
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

    pub fn connect_diverter(&mut self, diverter: Arc<RwLock<PVDiverter>>) {
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
            diverter.read().capacity_already_in_use(energy_supplied);
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

#[derive(Debug)]
pub struct PVDiverter {
    storage_tank: Arc<Mutex<StorageTank>>,
    immersion_heater: Arc<Mutex<ImmersionHeater>>,
    heat_source_name: String,
    capacity_already_in_use: AtomicF64,
}

impl PVDiverter {
    pub fn new(
        storage_tank: Arc<Mutex<StorageTank>>,
        heat_source: Arc<Mutex<ImmersionHeater>>,
        heat_source_name: String,
    ) -> Arc<RwLock<Self>> {
        let diverter = Arc::new(RwLock::new(Self {
            storage_tank,
            heat_source_name,
            immersion_heater: heat_source.clone(),
            capacity_already_in_use: Default::default(),
        }));

        heat_source.lock().connect_diverter(diverter.clone());

        diverter
    }

    /// Record heater output that would happen anyway, to avoid double-counting
    pub fn capacity_already_in_use(&self, energy_supplied: f64) {
        self.capacity_already_in_use
            .fetch_add(energy_supplied, Ordering::SeqCst);
    }

    /// Divert as much surplus as possible to the heater
    ///
    /// Arguments:
    /// * `supply_surplus` - surplus energy, in kWh, available to be diverted (negative by convention)
    pub fn divert_surplus(
        &self,
        supply_surplus: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        // check how much spare capacity the immersion heater has
        let imm_heater_max_capacity_spare =
            self.immersion_heater
                .lock()
                .energy_output_max(simulation_time_iteration, true)
                - self.capacity_already_in_use.load(Ordering::SeqCst);

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
    temp_internal_air_accessor: TempInternalAirAccessor,
    heat_output_collector_loop: f64,
    energy_supplied: f64,
    cp: f64,
    air_temp_coll_loop: f64, //mutating internally
    inlet_temp: f64,
}

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
        temp_internal_air_accessor: TempInternalAirAccessor,
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
            temp_internal_air_accessor,
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
        temp_storage_tank_s3_n: &[f64],
        simulation_time: &SimulationTimeIteration,
    ) -> f64 {
        // Air temperature in a heated space in the building
        let air_temp_heated_room = self.temp_internal_air_accessor.call();

        self.air_temp_coll_loop = match self.sol_loc {
            SolarCellLocation::Hs => air_temp_heated_room,
            SolarCellLocation::Nhs => {
                (air_temp_heated_room + self.external_conditions.air_temp(simulation_time)) / 2.
            }
            SolarCellLocation::Out => self.external_conditions.air_temp(simulation_time),
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
    use crate::core::schedule::WaterScheduleEventType;
    use crate::core::space_heat_demand::thermal_bridge::ThermalBridging;
    use crate::core::space_heat_demand::zone::Zone;
    use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
    use crate::corpus::HeatSource;
    use crate::external_conditions::{
        DaylightSavingsConfig, ShadingObject, ShadingObjectType, ShadingSegment,
    };
    use crate::input::FuelType;
    use crate::simulation_time::SimulationTime;
    use approx::assert_relative_eq;
    use pretty_assertions::assert_eq;
    use rstest::*;

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
    pub fn external_conditions(
        simulation_time_for_storage_tank: SimulationTime,
    ) -> Arc<ExternalConditions> {
        let air_temps = vec![0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0];
        let wind_speeds = vec![3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4];
        let wind_directions = vec![0.0; 8];
        let diffuse_horizontal_radiations = vec![333., 610., 572., 420., 0., 10., 90., 275.];
        let direct_beam_radiations = vec![420., 750., 425., 500., 0., 40., 0., 388.];
        let solar_reflectivity_of_ground = vec![0.2; 8760];
        let shading_segments = vec![
            ShadingSegment {
                number: 1,
                start: 180.,
                end: 135.,
                objects: None,
            },
            ShadingSegment {
                number: 2,
                start: 135.,
                end: 90.,
                objects: None,
            },
            ShadingSegment {
                number: 3,
                start: 90.,
                end: 45.,
                objects: None,
            },
            ShadingSegment {
                number: 4,
                start: 45.,
                end: 0.,
                objects: Some(vec![ShadingObject {
                    object_type: ShadingObjectType::Obstacle,
                    height: 10.5,
                    distance: 120.,
                }]),
            },
            ShadingSegment {
                number: 5,
                start: 0.,
                end: -45.,
                objects: None,
            },
            ShadingSegment {
                number: 6,
                start: -45.,
                end: -90.,
                objects: None,
            },
            ShadingSegment {
                number: 7,
                start: -90.,
                end: -135.,
                objects: None,
            },
            ShadingSegment {
                number: 8,
                start: -135.,
                end: -180.,
                objects: None,
            },
        ];

        Arc::new(ExternalConditions::new(
            &simulation_time_for_storage_tank.iter(),
            air_temps,
            wind_speeds,
            wind_directions, // change to 8-length of 0.0 values when migrating to 0.30
            diffuse_horizontal_radiations,
            direct_beam_radiations,
            solar_reflectivity_of_ground,
            51.383,
            -0.783,
            0,
            212,
            Some(212),
            1.0,
            Some(1),
            DaylightSavingsConfig::NotApplicable,
            false,
            false,
            shading_segments,
        ))
    }

    #[fixture]
    pub fn temp_internal_air_accessor(
        external_conditions: Arc<ExternalConditions>,
        simulation_time_for_storage_tank: SimulationTime,
    ) -> TempInternalAirAccessor {
        let be_objs = IndexMap::from([]);
        let ve_objs = vec![];
        let thermal_bridging = ThermalBridging::Bridges(IndexMap::from([]));

        let zone = Zone::new(
            500., // any number above 0
            0.,   // any number
            be_objs,
            thermal_bridging,
            ve_objs,
            None,
            0., // any number
            0., // any number
            external_conditions,
            &simulation_time_for_storage_tank.iter(),
        );

        let zones = IndexMap::from([("zone one".to_string(), zone)]).into();
        TempInternalAirAccessor {
            zones,
            total_volume: 0., // any number
        }
    }

    #[fixture]
    pub fn storage_tank(
        simulation_time_for_storage_tank: SimulationTime,
        cold_water_source: Arc<ColdWaterSource>,
        control_for_storage_tank: Arc<Control>,
        external_conditions: Arc<ExternalConditions>,
        temp_internal_air_accessor: TempInternalAirAccessor,
    ) -> ((StorageTank, StorageTank), Arc<RwLock<EnergySupply>>) {
        let energy_supply = Arc::new(RwLock::new(EnergySupply::new(
            FuelType::Electricity,
            simulation_time_for_storage_tank.total_steps(),
            None,
            None,
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
                temp_internal_air_accessor.clone(),
                external_conditions.clone(),
                None,
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
                                    50.,
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
                temp_internal_air_accessor.clone(), // this may need to be a different instance
                external_conditions,
                None,
                None,
                Some(energy_supply_conns.1),
                None,
                WATER.clone(),
            ),
        );

        (storage_tanks, energy_supply)
    }

    #[rstest]
    pub fn test_demand_hot_water(
        simulation_time_for_storage_tank: SimulationTime,
        storage_tank: ((StorageTank, StorageTank), Arc<RwLock<EnergySupply>>),
    ) {
        let ((mut storage_tank1, mut storage_tank2), energy_supply) = storage_tank;
        let usage_events = vec![
            vec![
                TypedScheduleEvent {
                    start: 6.,
                    duration: Some(6.),
                    temperature: 41.0,
                    name: "IES".to_string(),
                    event_type: WaterScheduleEventType::Shower,
                    warm_volume: None,
                    pipework_volume: None,
                },
                TypedScheduleEvent {
                    start: 6.,
                    duration: Some(6.),
                    temperature: 41.0,
                    name: "mixer".to_string(),
                    event_type: WaterScheduleEventType::Shower,
                    warm_volume: Some(48.),
                    pipework_volume: None,
                },
                TypedScheduleEvent {
                    start: 6.,
                    duration: Some(20.),
                    temperature: 43.0,
                    name: "medium".to_string(),
                    event_type: WaterScheduleEventType::Bath,
                    warm_volume: Some(100.),
                    pipework_volume: None,
                },
                TypedScheduleEvent {
                    start: 6.,
                    duration: Some(1.),
                    temperature: 40.0,
                    name: "other".to_string(),
                    event_type: WaterScheduleEventType::Other,
                    warm_volume: Some(8.),
                    pipework_volume: None,
                },
            ],
            vec![TypedScheduleEvent {
                start: 7.,
                duration: Some(6.),
                temperature: 41.0,
                name: "mixer".to_string(),
                event_type: WaterScheduleEventType::Shower,
                warm_volume: Some(48.),
                pipework_volume: None,
            }],
            vec![],
            vec![TypedScheduleEvent {
                start: 9.,
                duration: Some(6.),
                temperature: 45.0,
                name: "mixer".to_string(),
                event_type: WaterScheduleEventType::Shower,
                warm_volume: Some(48.),
                pipework_volume: None,
            }],
            vec![],
            vec![TypedScheduleEvent {
                start: 11.,
                duration: Some(6.5),
                temperature: 41.0,
                name: "mixer".to_string(),
                event_type: WaterScheduleEventType::Shower,
                warm_volume: Some(52.),
                pipework_volume: None,
            }],
            vec![],
            vec![],
        ];

        //  Expected results for the unit test
        let expected_temperatures_1 = [
            [55.0, 55.0, 55.0, 55.0],
            [
                15.45,
                54.595555555555556,
                54.595555555555556,
                54.595555555555556,
            ],
            [
                15.45,
                54.19530534979424,
                54.19530534979424,
                54.19530534979424,
            ],
            [10.5, 15.4, 53.38820740740741, 53.80385185185185],
            [55.0, 55.0, 55.0, 55.0],
            [
                13.4,
                54.595555555555556,
                54.595555555555556,
                54.595555555555556,
            ],
            [
                13.4,
                54.19530534979424,
                54.19530534979424,
                54.19530534979424,
            ],
            [
                13.4,
                53.79920588690749,
                53.79920588690749,
                53.79920588690749,
            ],
        ];

        let expected_temperatures_2 = [
            [10.0, 24.55880864197531, 60.0, 60.0],
            [
                10.06,
                16.317728395061728,
                39.76012654320988,
                59.687654320987654,
            ],
            [
                10.06,
                16.315472915714068,
                39.59145897824265,
                59.37752591068435,
            ],
            [10.34, 12.28, 24.509163580246913, 46.392706790123455],
            [10.34, 12.28, 60.0, 60.0],
            [10.74, 11.1, 60.0, 60.0],
            [10.74, 11.1, 59.687654320987654, 59.687654320987654],
            [10.74, 11.1, 59.37752591068435, 59.37752591068435],
        ];

        let expected_energy_supplied_1 = [
            5.913725648148166,
            0.0,
            0.0,
            0.0,
            3.858245898765432,
            0.0,
            0.0,
            0.0,
        ];

        let expected_energy_supplied_2 = [
            0.6720797510288063,
            0.0,
            0.0,
            0.0,
            3.033920793930041,
            1.8039389176954712,
            0.0,
            0.0,
        ];

        // Loop through the timesteps and the associated data pairs using `subTest`
        for (t_idx, t_it) in simulation_time_for_storage_tank.iter().enumerate() {
            let usage_events_for_iteration = usage_events[t_idx].clone();
            storage_tank1.demand_hot_water(usage_events_for_iteration.clone(), t_it);

            // Verify the temperatures against expected results
            assert_eq!(
                storage_tank1.temp_n, expected_temperatures_1[t_idx],
                "incorrect temperatures returned"
            );

            assert_relative_eq!(
                energy_supply.read().results_by_end_user()["immersion"][t_idx],
                expected_energy_supplied_1[t_idx],
                max_relative = 1e-6
            );

            storage_tank2.demand_hot_water(usage_events_for_iteration, t_it);

            assert_eq!(
                storage_tank2.temp_n, expected_temperatures_2[t_idx],
                "incorrect temperatures returned"
            );

            assert_relative_eq!(
                energy_supply.read().results_by_end_user()["immersion2"][t_idx],
                expected_energy_supplied_2[t_idx],
                max_relative = 1e-6
            );
        }
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
                Arc::new(RwLock::new(EnergySupply::new(
                    FuelType::MainsGas,
                    simulation_time_for_immersion_heater.total_steps(),
                    None,
                    None,
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
    pub fn test_demand_energy_for_immersion_heater(
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
    fn external_conditions_for_solar_thermal() -> Arc<ExternalConditions> {
        let simulation_time = SimulationTime::new(5088., 5112., 1.);

        let air_temps = vec![
            19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0,
            19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0,
        ];
        let wind_speeds = vec![
            3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1, 3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1, 3.9,
            3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1,
        ];
        let wind_directions = vec![
            30.0, 250., 220., 180., 150., 120., 100., 80., 60., 40., 20., 10., 50., 100., 140.,
            190., 200., 320., 330., 340., 350., 355., 315., 5.,
        ];
        let diffuse_horizontal_radiations = vec![
            0., 0., 0., 0., 35., 73., 139., 244., 320., 361., 369., 348., 318., 249., 225., 198.,
            121., 68., 19., 0., 0., 0., 0., 0.,
        ];
        let direct_beam_radiations = vec![
            0., 0., 0., 0., 0., 0., 7., 53., 63., 164., 339., 242., 315., 577., 385., 285., 332.,
            126., 7., 0., 0., 0., 0., 0.,
        ];
        let solar_reflectivity_of_ground = vec![
            0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
            0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
        ];
        let shading_segments = vec![
            ShadingSegment {
                number: 1,
                start: 180.,
                end: 135.,
                objects: None,
            },
            ShadingSegment {
                number: 2,
                start: 135.,
                end: 90.,
                objects: None,
            },
            ShadingSegment {
                number: 3,
                start: 90.,
                end: 45.,
                objects: None,
            },
            ShadingSegment {
                number: 4,
                start: 45.,
                end: 0.,
                objects: Some(vec![ShadingObject {
                    object_type: ShadingObjectType::Obstacle,
                    height: 10.5,
                    distance: 12.,
                }]),
            },
            ShadingSegment {
                number: 5,
                start: 0.,
                end: -45.,
                objects: None,
            },
            ShadingSegment {
                number: 6,
                start: -45.,
                end: -90.,
                objects: None,
            },
            ShadingSegment {
                number: 7,
                start: -90.,
                end: -135.,
                objects: None,
            },
            ShadingSegment {
                number: 8,
                start: -135.,
                end: -180.,
                objects: None,
            },
        ];

        Arc::new(ExternalConditions::new(
            &simulation_time.iter(),
            air_temps,
            wind_speeds,
            wind_directions,
            diffuse_horizontal_radiations,
            direct_beam_radiations,
            solar_reflectivity_of_ground,
            51.383,
            -0.783,
            0,
            212,
            Some(212),
            1.0,
            Some(1),
            DaylightSavingsConfig::NotApplicable,
            false,
            false,
            shading_segments,
        ))
    }

    #[fixture]
    pub fn storage_tank_with_solar_thermal(
        external_conditions_for_solar_thermal: Arc<ExternalConditions>,
        temp_internal_air_accessor: TempInternalAirAccessor,
    ) -> (
        StorageTank,
        Arc<Mutex<SolarThermalSystem>>,
        SimulationTime,
        Arc<RwLock<EnergySupply>>,
    ) {
        let cold_water_temps = [
            17.0, 17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7, 17.0, 17.1, 17.2, 17.3, 17.4, 17.5,
            17.6, 17.7, 17.0, 17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7,
        ];
        let simulation_time = SimulationTime::new(5088., 5112., 1.);
        let cold_feed = WaterSourceWithTemperature::ColdWaterSource(Arc::new(
            ColdWaterSource::new(cold_water_temps.to_vec(), &simulation_time, 1.),
        ));
        let energy_supply = Arc::new(RwLock::new(EnergySupply::new(
            FuelType::Electricity,
            simulation_time.total_steps(),
            None,
            None,
            None,
        )));
        let energy_supply_conn =
            EnergySupply::connection(energy_supply.clone(), "solarthermal").unwrap();
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
            external_conditions_for_solar_thermal.clone(),
            temp_internal_air_accessor.clone(),
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
            temp_internal_air_accessor.clone(),
            external_conditions_for_solar_thermal,
            None,
            None,
            None,
            None,
            WATER.clone(),
        );

        (storage_tank, solar_thermal, simulation_time, energy_supply)
    }

    #[rstest]
    // in Python this test is called test_demand_hot_water and is from test_storage_tank_with_solar_thermal.py
    pub fn test_demand_hot_water_for_storage_tank_with_solar_thermal(
        storage_tank_with_solar_thermal: (
            StorageTank,
            Arc<Mutex<SolarThermalSystem>>,
            SimulationTime,
            Arc<RwLock<EnergySupply>>,
        ),
    ) {
        let (mut storage_tank, solar_thermal, simulation_time, energy_supply) =
            storage_tank_with_solar_thermal;

        let expected_energy_demands = [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.3943936789277888,
            0.8431700006423382,
            1.3873931365189958,
            1.0919832923582113,
            1.1503273689665232,
            1.484482745066628,
            0.9003339693974624,
            0.49807749362100157,
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
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.3943936789277888,
            0.8431700006423382,
            1.3873931365189958,
            1.0919832923582113,
            1.1503273689665232,
            1.484482745066628,
            0.9003339693974624,
            0.49807749362100157,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
        ];
        let expected_energy_supplied = [
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.3943936789277888,
            0.8431700006423382,
            1.3873931365189958,
            1.0919832923582113,
            1.1503273689665232,
            1.484482745066628,
            0.9003339693974624,
            0.49807749362100157,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
        ];
        let expected_energy_supply_results = [
            10., 10., 10., 10., 10., 10., 10., 10., 110., 110., 110., 110., 110., 110., 110., 110.,
            10., 10., 10., 10., 10., 10., 10., 10.,
        ];

        let usage_events: [Vec<TypedScheduleEvent>; 24] = [
            vec![],
            vec![TypedScheduleEvent {
                start: 7.,
                duration: Some(6.),
                temperature: 41.0,
                event_type: WaterScheduleEventType::Shower,
                name: "mixer".to_owned(),
                warm_volume: Some(48.0),
                pipework_volume: None,
            }],
            vec![],
            vec![TypedScheduleEvent {
                start: 9.,
                duration: Some(6.),
                temperature: 45.0,
                event_type: WaterScheduleEventType::Shower,
                name: "mixer".to_owned(),
                warm_volume: Some(48.0),
                pipework_volume: None,
            }],
            vec![],
            vec![TypedScheduleEvent {
                start: 11.,
                duration: Some(6.5),
                temperature: 41.0,
                event_type: WaterScheduleEventType::Shower,
                name: "mixer".to_owned(),
                warm_volume: Some(52.0),
                pipework_volume: None,
            }],
            vec![],
            vec![],
            vec![],
            vec![TypedScheduleEvent {
                start: 7.,
                duration: Some(6.),
                temperature: 41.0,
                event_type: WaterScheduleEventType::Shower,
                name: "mixer".to_owned(),
                warm_volume: Some(48.0),
                pipework_volume: None,
            }],
            vec![],
            vec![TypedScheduleEvent {
                start: 9.,
                duration: Some(6.),
                temperature: 45.0,
                event_type: WaterScheduleEventType::Shower,
                name: "mixer".to_owned(),
                warm_volume: Some(48.0),
                pipework_volume: None,
            }],
            vec![],
            vec![TypedScheduleEvent {
                start: 11.,
                duration: Some(6.5),
                temperature: 41.0,
                event_type: WaterScheduleEventType::Shower,
                name: "mixer".to_owned(),
                warm_volume: Some(52.0),
                pipework_volume: None,
            }],
            vec![],
            vec![],
            vec![],
            vec![TypedScheduleEvent {
                start: 7.,
                duration: Some(6.),
                temperature: 41.0,
                event_type: WaterScheduleEventType::Shower,
                name: "mixer".to_owned(),
                warm_volume: Some(48.0),
                pipework_volume: None,
            }],
            vec![],
            vec![TypedScheduleEvent {
                start: 9.,
                duration: Some(6.),
                temperature: 45.0,
                event_type: WaterScheduleEventType::Shower,
                name: "mixer".to_owned(),
                warm_volume: Some(48.0),
                pipework_volume: None,
            }],
            vec![],
            vec![TypedScheduleEvent {
                start: 11.,
                duration: Some(6.5),
                temperature: 41.0,
                event_type: WaterScheduleEventType::Shower,
                name: "mixer".to_owned(),
                warm_volume: Some(52.0),
                pipework_volume: None,
            }],
            vec![],
            vec![],
        ];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            storage_tank.demand_hot_water(usage_events.get(t_idx).unwrap().clone(), t_it);
            assert_relative_eq!(
                storage_tank.test_energy_demand(),
                expected_energy_demands[t_idx],
                max_relative = 1e-7
            );
            assert_relative_eq!(
                solar_thermal.lock().test_energy_potential(),
                expected_energy_potentials[t_idx],
                max_relative = 1e-7
            );
            assert_relative_eq!(
                solar_thermal.lock().test_energy_supplied(),
                expected_energy_supplied[t_idx],
                max_relative = 1e-7
            );
            assert_relative_eq!(
                energy_supply.read().results_by_end_user()["solarthermal"][t_idx],
                expected_energy_supply_results[t_idx],
                max_relative = 1e-7
            );
        }
    }
}
