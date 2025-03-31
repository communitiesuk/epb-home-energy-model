use crate::compare_floats::{max_of_2, min_of_2};
use crate::core::common::WaterSourceWithTemperature;
use crate::core::controls::time_control::{Control, ControlBehaviour};
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::material_properties::MaterialProperties;
use crate::core::pipework::{Pipework, PipeworkLocation, Pipeworkesque};
use crate::core::schedule::TypedScheduleEvent;
use crate::core::units::WATTS_PER_KILOWATT;
use crate::core::water_heat_demand::misc::frac_hot_water;
use crate::corpus::{HeatSource, TempInternalAirFn};
use crate::external_conditions::ExternalConditions;
use crate::input::{SolarCellLocation, WaterPipework};
use crate::simulation_time::SimulationTimeIteration;
use anyhow::{anyhow, bail};
use atomic_float::AtomicF64;
use derivative::Derivative;
use indexmap::IndexMap;
use itertools::Itertools;
use ordered_float::OrderedFloat;
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
pub(crate) struct PositionedHeatSource {
    pub heat_source: Arc<Mutex<HeatSource>>,
    pub heater_position: f64,
    pub thermostat_position: Option<f64>,
}

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
#[derive(Derivative)]
#[derivative(Debug)]
pub struct StorageTank {
    init_temp: f64,
    q_std_ls_ref: f64, // measured standby losses due to cylinder insulation at standardised conditions, in kWh/24h
    cold_feed: WaterSourceWithTemperature,
    simulation_timestep: f64,
    energy_supply_conn_unmet_demand: Option<EnergySupplyConnection>,
    nb_vol: usize,
    #[derivative(Debug = "ignore")]
    temp_internal_air_fn: TempInternalAirFn,
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
    primary_gains: AtomicF64,    // primary pipework gains for a timestep (mutates over lifetime)
    #[cfg(test)]
    energy_demand_test: f64,
    temp_final_drawoff: Option<f64>, // In Python this is created from inside allocate_hot_water()
    temp_average_drawoff: Option<f64>, // In Python this is created from inside allocate_hot_water()
    temp_average_drawoff_volweighted: Option<f64>, // In Python this is created from inside allocate_hot_water()
    total_volume_drawoff: Option<f64>, // In Python this is created from inside allocate_hot_water()
}

pub(crate) struct StorageTankDetailedResult {} // TODO implement detailed results for StorageTank

impl StorageTank {
    /// Arguments:
    /// * `volume` - total volume of the tank, in litres
    /// * `losses` - measured standby losses due to cylinder insulation
    ///                                at standardised conditions, in kWh/24h
    /// * `init_temp` - initial temperature required for DHW
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
    pub(crate) fn new(
        volume: f64,
        losses: f64,
        init_temp: f64,
        cold_feed: WaterSourceWithTemperature,
        simulation_timestep: f64,
        heat_sources: IndexMap<String, PositionedHeatSource>,
        // In Python this is "project" but only temp_internal_air is accessed from it
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions: Arc<ExternalConditions>,
        nb_vol: Option<usize>,
        primary_pipework_lst: Option<&Vec<WaterPipework>>,
        energy_supply_conn_unmet_demand: Option<EnergySupplyConnection>,
        contents: MaterialProperties,
        _detailed_output_heating_cooling: bool, // TODO implement logic for this to match Python 0.32
    ) -> Self {
        let q_std_ls_ref = losses;

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
        let temp_n = vec![init_temp; nb_vol];

        #[cfg(test)]
        let energy_demand_test = 0.;

        // primary_pipework_losses_kwh added for reporting
        let primary_pipework_losses_kwh = 0.;
        let storage_losses_kwh = 0.;

        let input_energy_adj_prev_timestep = 0.;

        let primary_pipework_lst: Option<Vec<Pipework>> =
            if let Some(primary_pipework_lst) = primary_pipework_lst {
                primary_pipework_lst
                    .iter()
                    .map(|pipework| pipework.to_owned().try_into().unwrap())
                    .collect::<Vec<Pipework>>()
                    .into()
            } else {
                None
            };

        // With pre-heatd storage tanks, there could be the situation of tanks without heat sources
        // They could just get warmed up with WWHRS water.
        let mut heat_source_data = heat_sources.clone();

        if !heat_sources.is_empty() {
            // sort heat source data in order from the bottom of the tank based on heater position
            heat_source_data = heat_source_data
                .iter()
                .sorted_by(|a, b| {
                    OrderedFloat(a.1.heater_position).cmp(&OrderedFloat(b.1.heater_position))
                })
                .rev()
                .map(|x| (x.0.to_owned(), x.1.to_owned()))
                .collect();
        }

        let heating_active: HashMap<String, bool> = heat_sources
            .iter()
            .map(|(name, _heat_source)| ((*name).clone(), false))
            .collect();

        Self {
            init_temp,
            q_std_ls_ref,
            cold_feed,
            simulation_timestep,
            energy_supply_conn_unmet_demand,
            nb_vol,
            temp_internal_air_fn,
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
            heat_source_data,
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

    /// Draw off hot water from the tank
    /// Energy calculation as per BS EN 15316-5:2017 Method A sections 6.4.3, 6.4.6, 6.4.7
    /// Modification of calculation based on volumes and actual temperatures for each layer of water in the tank
    /// instead of the energy stored in the layer and a generic temperature (self.temp_out_w_min) = min_temp
    /// to decide if the tank can satisfy the demand (this was producing unnecesary unmet demand for strict high
    /// temp_out_w_min values
    /// Arguments:
    /// * `usage_events` -- All draw off events for the timestep
    pub(crate) fn demand_hot_water(
        &mut self,
        mut usage_events: Option<Vec<TypedScheduleEvent>>,
        simulation_time: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64, f64, f64, f64)> {
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
            .iter_mut()
            .flatten()
            .filter(|e| e.warm_volume.is_some())
            .collect_vec();

        for event in filtered_events {
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

            let (volume_used, temp_s3_n_step, energy_withdrawn, energy_unmet, rearrange) =
                self.allocate_hot_water(event.clone(), simulation_time);

            // Re-arrange the temperatures in the storage after energy input from pre-heated tank
            temp_s3_n = temp_s3_n_step;

            if rearrange {
                temp_s3_n = self.rearrange_temperatures(&temp_s3_n).1
            }

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
            }
            _ => temp_s3_n.last().copied(),
        };

        // TODO (from Python) 6.4.3.6 STEP 4 Volume to be withdrawn from the storage (for Heating)
        // TODO (from Python) - 6.4.3.7 STEP 5 Temperature of the storage after volume withdrawn (for Heating)

        // Run over multiple heat sources
        let mut temp_after_prev_heat_source = temp_s3_n.clone();
        let mut q_ls = 0.0;
        self.q_ls_n_prev_heat_source = vec![0.0; self.nb_vol];
        // In Python extra variables initialized and assigned here
        // for the purpose of passing them to the testoutput method
        // which we have decided to port (for now)
        let mut temp_s8_n = vec![0.; self.nb_vol];

        for (heat_source_name, positioned_heat_source) in self.heat_source_data.clone() {
            let (_, _setpntmax) = positioned_heat_source
                .heat_source
                .lock()
                .temp_setpnt(simulation_time)?;
            let heater_layer =
                (positioned_heat_source.heater_position * self.nb_vol as f64) as usize;
            let thermostat_layer =
                (positioned_heat_source.thermostat_position.ok_or_else(|| {
                    anyhow!("expected thermostat position on storage tank heat source")
                })? * self.nb_vol as f64) as usize;
            let TemperatureCalculation {
                temp_s8_n: temp_s8_n_step,
                q_ls: q_ls_this_heat_source,
                q_ls_n: q_ls_n_this_heat_source,
                ..
            } = self.run_heat_sources(
                temp_after_prev_heat_source.clone(),
                positioned_heat_source.heat_source.clone(),
                &heat_source_name,
                heater_layer,
                thermostat_layer,
                &self.q_ls_n_prev_heat_source.clone(),
                simulation_time,
            )?;

            temp_after_prev_heat_source = temp_s8_n_step.clone();
            q_ls += q_ls_this_heat_source;

            for (i, q_ls_n) in q_ls_n_this_heat_source.iter().enumerate() {
                self.q_ls_n_prev_heat_source[i] += q_ls_n;
            }

            temp_s8_n = temp_s8_n_step;

            // Trigger heating to stop
            self.determine_heat_source_switch_off(
                temp_s8_n.clone(),
                heat_source_name,
                positioned_heat_source,
                heater_layer,
                thermostat_layer,
                simulation_time,
            )?;
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
        Ok((
            q_use_w,
            q_unmet_w,
            self.temp_final_drawoff.unwrap(),
            self.temp_average_drawoff.unwrap(),
            self.total_volume_drawoff.unwrap(),
        ))
        // Sending temp_final_drawoff, temp_average_drawoff
        // for pipework loss and internal gains calculations
    }

    /// Allocate hot water layers to meet a single temperature demand.
    ///
    /// Arguments:
    /// * `event` -- Dictionary containing information about the draw-off event
    ///              (e.g. {'start': 18, 'duration': 1, 'temperature': 41.0, 'type': 'Other', 'name': 'other', 'warm_volume': 8.0})
    fn allocate_hot_water(
        &mut self,
        event: TypedScheduleEvent,
        simulation_time: SimulationTimeIteration,
    ) -> (f64, Vec<f64>, f64, f64, bool) {
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
        let mut pipework_temp = self.cold_feed.temperature(simulation_time, None); // This value set to initialise, but is never used - overwritten later.

        let mut pipework_considered = event.pipework_volume.unwrap() <= 0.0;

        let mut temp_average_drawoff_volweighted: f64 =
            self.temp_average_drawoff_volweighted.unwrap_or(0.);
        let mut total_volume_drawoff: f64 = self.total_volume_drawoff.unwrap_or(0.);
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
                self.cold_feed.temperature(simulation_time, None),
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
                    * (layer_temp - self.cold_feed.temperature(simulation_time, None))
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
            * (warm_temp - self.cold_feed.temperature(simulation_time, None));

        //  Calculate the remaining total volume
        let remaining_total_volume: f64 = remaining_vols.iter().sum();

        //  Calculate the total volume used
        let volume_used = self.volume_total_in_litres - remaining_total_volume;

        // Determine the new temperature distribution after displacement
        // Now that pre-heated sources can be the 'cold' feed, rearrangement of temperaturs, that used to
        // only happen before after the input from heat sources, could be required after the displacement
        // of water bringing new water from the 'cold' feed that could be warmer than the existing one.
        // flag is calculated for that purpose.
        let (new_temp_distribution, flag_rearrange_layers) =
            self.calculate_new_temperatures(remaining_vols, simulation_time);

        // Return the remaining storage volumes, volume used, new temperature distribution, the met/unmet targets, and flag to rearrange layers
        (
            volume_used,
            new_temp_distribution,
            energy_withdrawn,
            energy_unmet,
            flag_rearrange_layers,
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
    ) -> (Vec<f64>, bool) {
        let mut new_temps = self.temp_n.clone();

        // If the 'cold' feed water is hotter than the existing water in the tank, rearrange will be needed.
        // as if it was a heat source coming from the cold feed.
        let mut flag_rearrange_layers = false;

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

            // Initialisation of min temperature of tank layers to compare eventually against
            // the 'cold' feed temperature to check if rearrangement is needed.
            let mut temp_layer_min = self.temp_n[i];

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

                    // Update min temperature of the tank so far.
                    if self.temp_n[j] < temp_layer_min {
                        temp_layer_min = self.temp_n[j]
                    }

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
                // This is when the tank gets refilled with 'cold' water from the 'cold' feed.
                // Amount/Volume wasn't important before as it was assumed an infinite amount at
                // the cold feed was available.
                // The pre-heated tank is limited in the amount of water that can be provided at
                // a given temperature, eventually resourting to its own cold feed. So cold feed
                // temperature for the tank depends on the volume required.
                let temp_cold_feed = self
                    .cold_feed
                    .temperature(simulation_time, Some(needed_volume)); // In Python the temperature methods have changed to take in a needed_volume but then it is never used.
                volume_weighted_temperature += needed_volume * temp_cold_feed;
                flag_rearrange_layers = temp_cold_feed > temp_layer_min;
            }

            // Calculate the new temperature for the current layer
            // Round to 2 decimals to match instrumentation limits and significant figures,
            // ensuring practical accuracy, computational efficiency, and avoiding minute e-18 differences.
            new_temps[i] = (100. * (volume_weighted_temperature / total_volume)).round() / 100.;
            remaining_vols[i] = total_volume;
        }
        (new_temps, flag_rearrange_layers)
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
            .collect::<Vec<f64>>();

        (q_h_sto_end, temp_s7_n.to_owned())
    }

    fn run_heat_sources(
        &mut self,
        temp_s3_n: Vec<f64>,
        heat_source: Arc<Mutex<HeatSource>>,
        heat_source_name: &str,
        heater_layer: usize,
        thermostat_layer: usize,
        q_ls_prev_heat_source: &[f64],
        simulation_time: SimulationTimeIteration,
    ) -> anyhow::Result<TemperatureCalculation> {
        // 6.4.3.8 STEP 6 Energy input into the storage
        // input energy delivered to the storage in kWh - timestep dependent
        let q_x_in_n = self.potential_energy_input(
            &temp_s3_n,
            heat_source.clone(),
            heat_source_name,
            heater_layer,
            thermostat_layer,
            simulation_time,
        )?;

        self.calculate_temperatures(
            &temp_s3_n,
            heat_source,
            q_x_in_n,
            heater_layer,
            q_ls_prev_heat_source,
            simulation_time,
        )
    }

    /// Energy input for the storage from the generation system
    /// (expressed per energy carrier X)
    /// Heat Source = energy carrier
    pub(crate) fn potential_energy_input(
        // Heat source. Addition of temp_s3_n as an argument
        &mut self,
        temp_s3_n: &[f64],
        heat_source: Arc<Mutex<HeatSource>>,
        heat_source_name: &str,
        heater_layer: usize,
        thermostat_layer: usize,
        simulation_time: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<f64>> {
        // initialise list of potential energy input for each layer
        let mut q_x_in_n = iter::repeat(0.).take(self.nb_vol).collect_vec();
        let expect_message = format!(
            "Expected temp flow to be set for storage tank with heat source: {heat_source_name}"
        );
        let temp_flow = self
            .temp_flow(heat_source.clone(), simulation_time)?
            .expect(&expect_message);

        let heat_source = &mut *(heat_source.lock()); // TODO: can we refactor this? use a RwLock?

        let energy_potential =
            if let HeatSource::Storage(HeatSourceWithStorageTank::Solar(ref solar_heat_source)) =
                heat_source
            {
                // we are passing the storage tank object to the SolarThermal as this needs to call back the storage tank (sic from Python)
                solar_heat_source
                    .lock()
                    .energy_output_max(self, temp_s3_n, &simulation_time)
            } else {
                self.determine_heat_source_switch_on(
                    temp_s3_n,
                    heat_source_name.to_string(),
                    heat_source,
                    heater_layer,
                    thermostat_layer,
                    simulation_time,
                )?;

                if self.heating_active[heat_source_name] {
                    // upstream Python uses duck-typing/ polymorphism here, but we need to be more explicit
                    let mut energy_potential = match heat_source {
                        HeatSource::Storage(HeatSourceWithStorageTank::Immersion(
                            immersion_heater,
                        )) => immersion_heater
                            .lock()
                            .energy_output_max(simulation_time, false),
                        HeatSource::Storage(HeatSourceWithStorageTank::Solar(_)) => unreachable!(), // this case was already covered in the first arm of this if let clause, so can't repeat here
                        HeatSource::Wet(heat_source_wet) => {
                            heat_source_wet.energy_output_max(temp_flow, simulation_time)?
                        }
                    };

                    // TODO (from Python) Consolidate checks for systems with/without primary pipework

                    if !matches!(
                        heat_source,
                        HeatSource::Storage(HeatSourceWithStorageTank::Immersion(_))
                    ) {
                        let (primary_pipework_losses_kwh, _) = self.primary_pipework_losses(
                            energy_potential,
                            temp_flow,
                            simulation_time,
                        );
                        energy_potential -= primary_pipework_losses_kwh;
                    }

                    energy_potential
                } else {
                    0.
                }
            };

        q_x_in_n[heater_layer] += energy_potential;

        Ok(q_x_in_n)
    }

    fn calculate_temperatures(
        &mut self,
        temp_s3_n: &[f64],
        heat_source: Arc<Mutex<HeatSource>>,
        q_x_in_n: Vec<f64>,
        heater_layer: usize,
        q_ls_n_prev_heat_source: &[f64],
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<TemperatureCalculation> {
        let temp_flow = self
            .temp_flow(heat_source.clone(), simulation_time_iteration)?
            .expect("Expected temp flow to be set for storage tank with heat source");

        let (q_s6, temp_s6_n) = self.energy_input(temp_s3_n, &q_x_in_n);

        // 6.4.3.9 STEP 7 Re-arrange the temperatures in the storage after energy input
        let (q_h_sto_s7, temp_s7_n) = self.rearrange_temperatures(&temp_s6_n);

        // STEP 8 Thermal losses and final temperature
        let (q_in_h_w, q_ls, temp_s8_n, q_ls_n) = self.thermal_losses(
            &temp_s7_n,
            &q_x_in_n,
            q_h_sto_s7,
            heater_layer,
            q_ls_n_prev_heat_source,
            temp_flow,
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

        Ok(TemperatureCalculation {
            temp_s8_n,
            q_x_in_n,
            q_s6,
            temp_s6_n,
            temp_s7_n,
            q_in_h_w,
            q_ls,
            q_ls_n,
        })
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

    /// Thermal losses are calculated with respect to the impact of the temperature set point
    pub fn thermal_losses(
        &mut self,
        temp_s7_n: &[f64],
        q_x_in_n: &[f64],
        q_h_sto_s7: Vec<f64>,
        heater_layer: usize,
        q_ls_n_prev_heat_source: &[f64],
        setpntmax: f64,
    ) -> (f64, f64, Vec<f64>, Vec<f64>) {
        let q_x_in_adj: f64 = q_x_in_n.iter().sum();

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
                * (min_of_2(temp_s7_n[i], setpntmax) - STORAGE_TANK_TEMP_AMB)
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
            // Need to check for heating input below, otherwise temperature will be wrongly
            // capped at setpntmax in cases where tank temperature already exceeded setpntmax
            // without any contribution from the heat source. This could happen if the setpntmax
            // setting is lower in the current timestep than the previous timestep, or if
            // another heat source has a higher setpntmax (and therefore heated the tank to a
            // higher temperature) than the heat source currently being considered
            let temp_s8_n_step = if q_x_in_adj > 0.0 && temp_s7_n[i] > setpntmax {
                // Case 2 - Temperature exceeding the set point
                setpntmax
            } else {
                // Case 1 - Temperature below the set point
                // TODO (from Python) - spreadsheet accounts for total thermal losses not just layer

                // the final value of the temperature
                // is reduced due to the effect of the thermal losses
                // Formula (14) in the standard appears to have error as addition not multiply
                // and P instead of rho
                temp_s7_n[i] - (q_ls_n[i] / (self.rho * self.cp * self.vol_n[i]))
            };
            temp_s8_n.push(temp_s8_n_step);
        }

        let q_in_h_w = if q_x_in_adj > 0.0 {
            // excess energy/ energy surplus
            // excess energy is calculated as the difference from the energy stored, Qsto,step7, and
            // energy stored once the set temperature is obtained, Qsto,step8, with addition of the
            // thermal losses.
            // Note: The surplus must be calculated only for those layers that the
            //       heat source currently being considered is capable of heating,
            //       i.e. excluding those below the heater position.
            let mut energy_surplus = 0.0;
            if temp_s7_n[heater_layer] > setpntmax {
                for i in heater_layer..self.nb_vol {
                    energy_surplus += q_h_sto_s7[i]
                        - q_ls_n[i]
                        - (self.rho * self.cp * self.vol_n[i] * setpntmax);
                }
            }
            // the thermal energy provided to the system (from heat sources) shall be limited
            // adjustment of the energy delivered to the storage according with the set temperature
            // potential input from generation
            // TODO (from Python code) - find in standard - availability of back-up - where is this from?
            // also referred to as electrical power on
            let sto_bu_on = 1.;
            min_of_2(q_x_in_adj - energy_surplus, q_x_in_adj * sto_bu_on)
        } else {
            0.
        };

        (q_in_h_w, q_ls, temp_s8_n, q_ls_n)
    }

    fn heat_source_output(
        &mut self,
        heat_source: Arc<Mutex<HeatSource>>,
        input_energy_adj: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let expect_message =
            "Expected set point max to be set for storage tank with this heat source".to_string();
        let temp_flow = self
            .temp_flow(heat_source.clone(), simulation_time_iteration)?
            .expect(&expect_message);

        match &mut *(heat_source.clone().lock()) {
            HeatSource::Storage(HeatSourceWithStorageTank::Immersion(immersion)) => immersion
                .lock()
                .demand_energy(input_energy_adj, simulation_time_iteration),
            HeatSource::Storage(HeatSourceWithStorageTank::Solar(solar)) => Ok(solar
                .lock()
                .demand_energy(input_energy_adj, simulation_time_iteration.index)),
            HeatSource::Wet(ref mut wet_heat_source) => {
                let (primary_pipework_losses_kwh, primary_gains) = self.primary_pipework_losses(
                    input_energy_adj,
                    temp_flow,
                    simulation_time_iteration,
                );
                let input_energy_adj = input_energy_adj + primary_pipework_losses_kwh;

                let heat_source_output = wet_heat_source.demand_energy(
                    input_energy_adj,
                    temp_flow,
                    simulation_time_iteration,
                )? - primary_pipework_losses_kwh;
                self.input_energy_adj_prev_timestep = input_energy_adj;
                self.primary_gains.store(primary_gains, Ordering::SeqCst);

                // TODO (from Python) - how are these gains reflected in the calculations? allocation by zone?
                Ok(heat_source_output)
            }
        }
    }

    /// No demand from heat source if the temperature of the tank at the
    /// thermostat position is below the set point
    /// Trigger heating to start when temperature falls below the minimum
    fn retrieve_setpnt(
        &self,
        heat_source: &HeatSource,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<(Option<f64>, Option<f64>)> {
        heat_source.temp_setpnt(simulation_time_iteration)
    }

    fn determine_heat_source_switch_on(
        &mut self,
        temp_s3_n: &[f64],
        heat_source_name: String,
        heat_source: &mut HeatSource,
        _heater_layer: usize,
        thermostat_layer: usize,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<()> {
        // TODO (from Python) - update for smarthot water tank to use SOC instead. - take out all areas where temp_setpnt functions is used and deals with systems switches on and off needs to be taken out into separate function - write as soc.
        //  look at heating ssystems that doesn't use if conditions implement like this.
        //  l205-217 separate function and call new function. in the new smart control object [same inputs ]
        //  TODO (from Python) remove unused parameters
        let (setpntmin, _) = self.retrieve_setpnt(heat_source, simulation_time_iteration)?;
        if setpntmin.is_some() && temp_s3_n[thermostat_layer] <= setpntmin.unwrap() {
            self.heating_active.insert(heat_source_name, true);
        };
        Ok(())
    }

    fn determine_heat_source_switch_off(
        &mut self,
        temp_s8_n: Vec<f64>,
        heat_source_name: String,
        heat_source: PositionedHeatSource,
        _heater_layer: usize,
        thermostat_layer: usize,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<()> {
        let heat_source = heat_source.heat_source;
        let (_, setpntmax) =
            self.retrieve_setpnt(&(heat_source.lock()), simulation_time_iteration)?;
        let expect_message = format!(
            "Expected set point max to be set for storage tank with heat source: {heat_source_name}"
        );
        if temp_s8_n[thermostat_layer] >= setpntmax.ok_or_else(|| anyhow!(expect_message))? {
            self.heating_active.insert(heat_source_name, false);
        };

        Ok(())
    }

    fn temp_flow(
        &self,
        heat_source: Arc<Mutex<HeatSource>>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<Option<f64>> {
        let (_, setpntmax) =
            self.retrieve_setpnt(&(heat_source.lock()), simulation_time_iteration)?;
        Ok(setpntmax)
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
            PipeworkLocation::Internal => (self.temp_internal_air_fn)(),
        }
    }

    pub(crate) fn get_cold_water_source(&self) -> &WaterSourceWithTemperature {
        &self.cold_feed
    }

    /// This is only used to calculate the equivalent volume of water for IES showers
    /// in order to get the energy content for the internal gains.
    /// Therefore the actual value used is not critical.
    /// It has been suggested/considered the use of the top layer of the storage tank
    /// but this could be similar to the cold feed temperature after big draw-offs
    /// To avoid any issues in those situations we use the setpoing temperature of the
    /// tank.
    pub(crate) fn get_temp_hot_water(&self) -> f64 {
        self.init_temp // Use intial temperature of tank as reference value.
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

        let (_q_s6, temp_s6_n) = self.energy_input(temp_s3_n, &q_x_in_n);

        // 6.4.3.9 STEP 7 Re-arrange the temperatures in the storage after energy input
        let (_q_h_sto_s7, temp_s7_n) = self.rearrange_temperatures(&temp_s6_n);

        // TODO (from Python) Check [0] is bottom layer temp and that solar thermal inlet is top layer NB_VOL-1
        (temp_s7_n[0], temp_s7_n[self.nb_vol - 1])
    }

    pub(crate) fn to_report(&self) -> (f64, f64) {
        (self.primary_pipework_losses_kwh, self.storage_losses_kwh)
    }

    // NB. there is a testoutput() function here in the Python to output to a test file - will not reimplement unless seen as necessary

    /// draw off hot water layers until required volume is provided.
    ///
    /// Arguments:
    /// * volume    -- volume of water required
    fn draw_off_hot_water(
        &mut self,
        volume: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        // Remaining volume of water in storage tank layers
        let mut remaining_vols = self.vol_n.clone();

        let mut remaining_demanded_volume = volume;

        // Initialize the unmet and met energies
        let mut _energy_withdrawn = 0.0;
        self.temp_average_drawoff_volweighted = Some(0.0);
        self.total_volume_drawoff = Some(0.0);
        self.temp_average_drawoff =
            Some(self.cold_feed.temperature(simulation_time_iteration, None));
        self.temp_final_drawoff = Some(self.cold_feed.temperature(simulation_time_iteration, None));
        let _temp_ini_n = self.temp_n.clone();
        let temp_s3_n = self.temp_n.clone();

        // Loop through storage layers (starting from the top)
        for layer_index in (0..self.temp_n.len()).rev() {
            let layer_temp = self.temp_n[layer_index];
            let layer_vol = remaining_vols[layer_index];

            //  This cannot happen in the preheated tank. Check!
            //  Skip this layer if its remaining volume is already zero
            // if remaining_vols[layer_index] <= 0.0:
            //     continue

            // Volume of water required at this layer

            let required_vol;
            if layer_vol <= remaining_demanded_volume {
                // This is the case where layer cannot meet all remaining demanded volume
                required_vol = layer_vol;
                remaining_vols[layer_index] -= layer_vol;
                remaining_demanded_volume -= layer_vol;
            } else {
                // This is the case where layer can meet all remaining demanded volume
                required_vol = remaining_demanded_volume;
                // Deduct the required volume from the remaining demand and update the layer's volume
                remaining_vols[layer_index] -= required_vol;
                remaining_demanded_volume = 0.0;
            }

            self.temp_average_drawoff_volweighted =
                Some(self.temp_average_drawoff_volweighted.unwrap() + required_vol * layer_temp);
            self.total_volume_drawoff = Some(self.total_volume_drawoff.unwrap() + required_vol);

            //  Record the met volume demand for the current temperature target
            //  warm_vol_removed is the volume of warm water that has been satisfied from hot water in this layer
            _energy_withdrawn +=
                //  Calculation with event water parameters
                // self.__rho * self.__Cp * warm_vol_removed * (warm_temp - self.__cold_feed.temperature())
                //  Calculation with layer water parameters
                self.rho
                    * self.cp
                    * required_vol
                    * (layer_temp - self.cold_feed.temperature(simulation_time_iteration, None));

            if remaining_demanded_volume <= 0.0 {
                break;
            }
        }

        if remaining_demanded_volume > 0.0 {
            self.temp_average_drawoff_volweighted = Some(
                self.temp_average_drawoff_volweighted.unwrap()
                    + remaining_demanded_volume
                        * self.cold_feed.temperature(
                            simulation_time_iteration,
                            remaining_demanded_volume.into(),
                        ),
            );
            self.total_volume_drawoff =
                Some(self.total_volume_drawoff.unwrap() + remaining_demanded_volume);
        }

        if self.total_volume_drawoff != Some(0.0) {
            self.temp_average_drawoff = Some(
                self.temp_average_drawoff_volweighted.unwrap() / self.total_volume_drawoff.unwrap(),
            );
        } else {
            self.temp_average_drawoff = temp_s3_n.last().cloned();
        }

        // Determine the new temperature distribution after displacement
        let (mut new_temp_distribution, flag_rearrange_layers) =
            self.calculate_new_temperatures(remaining_vols, simulation_time_iteration);

        if flag_rearrange_layers {
            // Re-arrange the temperatures in the storage after energy input from pre-heated tank
            (_, new_temp_distribution) = self.rearrange_temperatures(&new_temp_distribution);
        }

        self.temp_n = new_temp_distribution.clone();

        // Return the remaining storage volumes, volume used, new temperature distribution, and the met/unmet targets
        self.temp_average_drawoff.unwrap()
    }
    fn additional_energy_input(
        &mut self,
        heat_source: Arc<Mutex<HeatSource>>,
        heat_source_name: &str,
        energy_input: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        if energy_input == 0. {
            return Ok(0.);
        }

        let heat_source_data = &self.heat_source_data[heat_source_name];

        let heater_layer = (heat_source_data.heater_position * self.nb_vol as f64) as usize;

        let mut q_x_in_n = vec![0.; self.nb_vol];
        q_x_in_n[heater_layer] = energy_input;
        let TemperatureCalculation {
            temp_s8_n,
            q_in_h_w,
            q_ls_n: q_ls_n_this_heat_source,
            ..
        } = self.calculate_temperatures(
            &self.temp_n.clone(),
            heat_source,
            q_x_in_n,
            heater_layer,
            &self.q_ls_n_prev_heat_source.clone(),
            simulation_time_iteration,
        )?;

        for (i, q_ls_n) in q_ls_n_this_heat_source.iter().enumerate() {
            self.q_ls_n_prev_heat_source[i] += *q_ls_n;
        }

        self.temp_n = temp_s8_n;

        Ok(q_in_h_w)
    }

    #[cfg(test)]
    fn test_energy_demand(&self) -> f64 {
        self.energy_demand_test
    }

    /// Return the DHW recoverable heat losses as internal gain for the current timestep in W
    pub fn internal_gains(&self) -> f64 {
        let primary_gains_timestep = self.primary_gains.load(Ordering::SeqCst);
        self.primary_gains
            .store(Default::default(), Ordering::SeqCst);

        (self.q_sto_h_ls_rbl.expect(
            "storage tank logic expects q_sto_h_ls_rbl to have been set internally at this point",
        ) * WATTS_PER_KILOWATT as f64
            / self.simulation_timestep)
            + primary_gains_timestep
    }

    fn primary_pipework_losses(
        &mut self,
        input_energy_adj: f64,
        setpnt_max: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (f64, f64) {
        let mut primary_pipework_losses_kwh = Default::default();
        let mut primary_gains_w = Default::default();
        if let Some(primary_pipework_lst) = self.primary_pipework_lst.as_ref() {
            // start of heating event
            if input_energy_adj > 0. && self.input_energy_adj_prev_timestep == 0. {
                for pipework_data in primary_pipework_lst {
                    primary_pipework_losses_kwh += pipework_data.cool_down_loss(
                        setpnt_max,
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
                        setpnt_max,
                        self.temp_surrounding_primary_pipework(
                            pipework_data,
                            simulation_time_iteration,
                        ),
                    );

                    // Check if pipework location is internal
                    let location = pipework_data.location();

                    if matches!(location, PipeworkLocation::Internal) {
                        primary_gains_w += primary_pipework_losses_w
                    }

                    primary_pipework_losses_kwh += primary_pipework_losses_w
                        * self.simulation_timestep
                        / WATTS_PER_KILOWATT as f64
                }
            }

            // end of heating event
            if input_energy_adj == 0. && self.input_energy_adj_prev_timestep > 0. {
                for pipework_data in primary_pipework_lst {
                    let location = pipework_data.location();
                    match location {
                        PipeworkLocation::External => {}
                        PipeworkLocation::Internal => {
                            primary_gains_w += pipework_data.cool_down_loss(
                                setpnt_max,
                                self.temp_surrounding_primary_pipework(
                                    pipework_data,
                                    simulation_time_iteration,
                                ),
                            ) * WATTS_PER_KILOWATT as f64
                                / self.simulation_timestep
                        }
                    }
                }
            }
        }

        // keeping primary_pipework_losses_kWh for reporting as part of investigation of issue #31225: FDEV A082
        self.primary_pipework_losses_kwh = primary_pipework_losses_kwh;

        (primary_pipework_losses_kwh, primary_gains_w)
    }

    /// Return the pre-heated water temperature for the current timestep
    pub(crate) fn temperature(
        &mut self,
        volume_needed: Option<f64>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        // This is only relevant when the storage tank is working as a pre-heated source.
        // If the volume required is 0.0 or not provided, the calculation assumes the requirement is
        // for the actual cold feed temperature of the pre-heated tank (eventually the real cold-feed)
        // Otherwise, it calculates the average water the tank can provided for the required volume.
        let volume_needed = volume_needed.unwrap_or(0.);
        if volume_needed == 0. {
            self.cold_feed.temperature(simulation_time_iteration, None)
        } else {
            self.draw_off_hot_water(volume_needed, simulation_time_iteration)
        }
    }

    pub(crate) fn output_results(&self) -> Option<Vec<StorageTankDetailedResult>> {
        todo!("Will be completed after 0.32 migration")
    }
}

#[derive(Debug, PartialEq)]
struct TemperatureCalculation {
    temp_s8_n: Vec<f64>,
    q_x_in_n: Vec<f64>,
    q_s6: f64,
    temp_s6_n: Vec<f64>,
    temp_s7_n: Vec<f64>,
    q_in_h_w: f64,
    q_ls: f64,
    q_ls_n: Vec<f64>,
}

#[derive(Debug)]
pub(crate) struct SmartHotWaterTank {}

impl SmartHotWaterTank {
    pub(crate) fn new() -> Self {
        todo!()
    }
}

#[derive(Clone, Debug)]
pub struct ImmersionHeater {
    pwr: f64, // rated power
    energy_supply_connection: EnergySupplyConnection,
    simulation_timestep: f64,
    control_min: Arc<Control>,
    control_max: Arc<Control>,
    diverter: Option<Arc<RwLock<PVDiverter>>>,
}

/// An object to represent an immersion heater
impl ImmersionHeater {
    /// Construct an ImmersionHeater object
    /// Arguments:
    /// * rated_power        -- in kW
    /// * energy_supply_conn -- reference to EnergySupplyConnection object
    /// * simulation_time    -- reference to SimulationTime object
    /// * controlmin            -- reference to a control object which must select current
    ///                         the minimum timestep temperature
    /// * controlmax            -- reference to a control object which must select current
    ///                         the maximum timestep temperature
    /// * diverter           -- reference to a PV diverter object
    pub(crate) fn new(
        rated_power: f64,
        energy_supply_connection: EnergySupplyConnection,
        simulation_timestep: f64,
        control_min: Arc<Control>,
        control_max: Arc<Control>,
    ) -> Self {
        Self {
            pwr: rated_power,
            energy_supply_connection,
            simulation_timestep,
            control_min,
            control_max,
            diverter: None,
        }
    }
    pub(crate) fn temp_setpnt(
        &self,
        simtime: SimulationTimeIteration,
    ) -> (Option<f64>, Option<f64>) {
        (
            self.control_min.setpnt(&simtime),
            self.control_max.setpnt(&simtime),
        )
    }

    pub fn connect_diverter(&mut self, diverter: Arc<RwLock<PVDiverter>>) {
        if self.diverter.is_some() {
            panic!("diverter was already connected");
        }

        self.diverter = Some(diverter);
    }

    /// Demand energy (in kWh) from the heater
    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        if energy_demand < 0.0 {
            bail!("Negative energy demand on ImmersionHeater");
        };

        // Account for time control. In the Python they also check here whether control_min is None
        // but this is not possible as it's a required field for an ImmersionHeater.
        let energy_supplied = if self.control_min.as_ref().is_on(simtime) {
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

        Ok(energy_supplied)
    }

    /// Calculate the maximum energy output (in kWh) from the heater
    pub fn energy_output_max(
        &self,
        simtime: SimulationTimeIteration,
        ignore_standard_control: bool,
    ) -> f64 {
        // Account for time control. In the Python they also check here whether control_min is None
        // but this is not possible as it's a required field for an ImmersionHeater.
        if self.control_min.as_ref().is_on(simtime) || ignore_standard_control {
            self.pwr * self.simulation_timestep
        } else {
            0.
        }
    }
}

/// Trait to represent a thing that can divert a surplus, like a PV diverter.
pub(crate) trait SurplusDiverting: Send + Sync {
    fn divert_surplus(
        &self,
        supply_surplus: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64>;
}

#[derive(Debug, Clone)]
pub(crate) enum HotWaterStorageTank {
    StorageTank(Arc<RwLock<StorageTank>>),
    SmartHotWaterTank(Arc<RwLock<SmartHotWaterTank>>),
}
#[derive(Debug)]
pub struct PVDiverter {
    pre_heated_water_source: HotWaterStorageTank,
    immersion_heater: Arc<Mutex<ImmersionHeater>>,
    heat_source_name: String,
    control_max: Option<Arc<Control>>,
    capacity_already_in_use: AtomicF64,
}

impl PVDiverter {
    pub(crate) fn new(
        storage_tank: &HotWaterStorageTank,
        heat_source: Arc<Mutex<ImmersionHeater>>,
        heat_source_name: String,
        control_max: Option<Arc<Control>>,
    ) -> Arc<RwLock<Self>> {
        let diverter = Arc::new(RwLock::new(Self {
            pre_heated_water_source: storage_tank.clone(),
            heat_source_name,
            control_max,
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

    pub fn timestep_end(&mut self) {
        self.capacity_already_in_use = Default::default();
    }
}

impl SurplusDiverting for PVDiverter {
    /// Divert as much surplus as possible to the heater
    ///
    /// Arguments:
    /// * `supply_surplus` - surplus energy, in kWh, available to be diverted (negative by convention)
    fn divert_surplus(
        &self,
        supply_surplus: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // check how much spare capacity the immersion heater has
        let imm_heater_max_capacity_spare = self
            .immersion_heater
            .lock()
            .energy_output_max(simulation_time_iteration, true)
            - self.capacity_already_in_use.load(Ordering::SeqCst);

        // Calculate the maximum energy that could be diverted
        // Note: supply_surplus argument is negative by convention, so negate it here
        let energy_diverted_max = min_of_2(imm_heater_max_capacity_spare, -supply_surplus);

        // Add additional energy to storage tank and calculate how much energy was accepted

        let energy_diverted = match &self.pre_heated_water_source {
            HotWaterStorageTank::StorageTank(storage_tank) => {
                storage_tank.write().additional_energy_input(
                    Arc::new(Mutex::new(HeatSource::Storage(
                        HeatSourceWithStorageTank::Immersion(self.immersion_heater.clone()),
                    ))),
                    &self.heat_source_name,
                    energy_diverted_max,
                    simulation_time_iteration,
                )?
            }
            HotWaterStorageTank::SmartHotWaterTank(_) => todo!("as part of migration 0.34"),
        };
        Ok(energy_diverted)
    }
}

/// The following code contains objects that represent solar thermal systems.
/// Method 3 in BS EN 15316-4-3:2017.
/// An object to represent a solar thermal system
#[derive(Derivative)]
#[derivative(Debug)]
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
    #[derivative(Debug = "ignore")]
    temp_internal_air_fn: TempInternalAirFn,
    heat_output_collector_loop: f64,
    energy_supplied: f64,
    control_max: Arc<Control>,
    cp: f64,
    air_temp_coll_loop: f64, //mutating internally
    inlet_temp: f64,
}

impl SolarThermalSystem {
    /// Construct a SolarThermalSystem object
    /// Arguments:
    /// * sol_loc         -- Location of the main part of the collector loop piping
    /// * area_module     -- Collector module reference area
    /// * modules         -- Number of collector modules installed
    /// * peak_collector_efficiency -- Peak collector efficiency
    /// * incidence_angle_modifier -- Hemispherical incidence angle modifier
    /// * first_order_hlc -- First order heat loss coefficient
    /// * second_order_hlc -- Second order heat loss coefficient
    /// * collector_mass_flow_rate -- Mass flow rate solar loop
    /// * power_pump      -- Power of collector pump
    /// * power_pump_control -- Power of collector pump controller
    /// * energy_supply_conn -- reference to EnergySupplyConnection object
    /// * tilt            -- is the tilt angle (inclination) of the PV panel from horizontal,
    ///                   measured upwards facing, 0 to 90, in degrees.
    ///                   0=horizontal surface, 90=vertical surface.
    ///                   Needed to calculate solar irradiation at the panel surface.
    /// * orientation     -- is the orientation angle of the inclined surface, expressed as the
    ///                   geographical azimuth angle of the horizontal projection of the inclined
    ///                   surface normal, -180 to 180, in degrees;
    ///                   Assumed N 180 or -180, E 90, S 0, W -90
    ///                   TODO - PV standard refers to angle as between 0 to 360?
    ///                   Needed to calculate solar irradiation at the panel surface.
    /// * solar_loop_piping_hlc -- Heat loss coefficient of the collector loop piping
    /// * ext_cond        -- reference to ExternalConditions object
    /// * simulation_time -- reference to SimulationTime object
    /// * contents        -- reference to MaterialProperties object
    /// * overshading     -- TODO could add at a later date. Feed into solar module
    /// * controlmax      -- reference to a control object which must select current the maximum timestep temperature
    pub(crate) fn new(
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
        temp_internal_air_fn: TempInternalAirFn,
        simulation_timestep: f64,
        control_max: Arc<Control>,
        contents: MaterialProperties,
        energy_supply_from_environment_conn: Option<EnergySupplyConnection>,
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
            temp_internal_air_fn,
            heat_output_collector_loop: 0.0,
            energy_supplied: 0.0,
            control_max,
            // Water specific heat in J/kg.K
            // (defined under eqn 51 on page 40 of BS EN ISO 15316-4-3:2017)
            cp: contents.specific_heat_capacity(),
            air_temp_coll_loop: Default::default(),
            inlet_temp: Default::default(),
        }
    }

    pub(crate) fn temp_setpnt(
        &self,
        simtime: SimulationTimeIteration,
    ) -> (Option<f64>, Option<f64>) {
        let temp_setpnt = self.control_max.setpnt(&simtime);
        (temp_setpnt, temp_setpnt)
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
        let air_temp_heated_room = (self.temp_internal_air_fn)();

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
    use crate::core::controls::time_control::SetpointTimeControl;
    use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyBuilder};
    use crate::core::material_properties::WATER;
    use crate::core::schedule::WaterScheduleEventType;
    use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
    use crate::corpus::HeatSource;
    use crate::external_conditions::{
        DaylightSavingsConfig, ShadingObject, ShadingObjectType, ShadingSegment,
    };
    use crate::input::{FuelType, WaterPipeContentsType};
    use crate::simulation_time::SimulationTime;
    use approx::assert_relative_eq;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[fixture]
    fn simulation_time_for_storage_tank() -> SimulationTime {
        SimulationTime::new(0., 8., 1.)
    }

    #[fixture]
    fn cold_water_source() -> Arc<ColdWaterSource> {
        Arc::new(ColdWaterSource::new(
            vec![10.0, 10.1, 10.2, 10.5, 10.6, 11.0, 11.5, 12.1],
            0,
            1.,
        ))
    }

    #[fixture]
    fn external_conditions(
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
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 2,
                start: 135.,
                end: 90.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 3,
                start: 90.,
                end: 45.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 4,
                start: 45.,
                end: 0.,
                shading_objects: Some(vec![ShadingObject {
                    object_type: ShadingObjectType::Obstacle,
                    height: 10.5,
                    distance: 120.,
                }]),
                ..Default::default()
            },
            ShadingSegment {
                number: 5,
                start: 0.,
                end: -45.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 6,
                start: -45.,
                end: -90.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 7,
                start: -90.,
                end: -135.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 8,
                start: -135.,
                end: -180.,
                shading_objects: None,
                ..Default::default()
            },
        ]
        .into();

        Arc::new(ExternalConditions::new(
            &simulation_time_for_storage_tank.iter(),
            air_temps,
            wind_speeds,
            wind_directions,
            diffuse_horizontal_radiations,
            direct_beam_radiations,
            solar_reflectivity_of_ground,
            51.383,
            -0.783,
            0,
            0,
            Some(0),
            1.0,
            Some(1),
            Some(DaylightSavingsConfig::NotApplicable),
            false,
            false,
            shading_segments,
        ))
    }

    #[fixture]
    fn temp_internal_air_fn() -> TempInternalAirFn {
        Arc::new(|| 20.)
    }

    #[fixture]
    fn energy_supply(
        simulation_time_for_storage_tank: SimulationTime,
    ) -> Arc<RwLock<EnergySupply>> {
        let energy_supply = EnergySupplyBuilder::new(
            FuelType::Electricity,
            simulation_time_for_storage_tank.total_steps(),
        )
        .build();

        Arc::new(RwLock::new(energy_supply))
    }

    fn heat_source(
        simulation_time_for_storage_tank: SimulationTime,
        energy_supply_connection: EnergySupplyConnection,
        rated_power: f64,
        heater_position: f64,
        thermostat_position: f64,
        control_min_schedule: Vec<Option<f64>>,
        control_max_schedule: Vec<Option<f64>>,
    ) -> PositionedHeatSource {
        let simulation_timestep = simulation_time_for_storage_tank.step;
        let control_min = SetpointTimeControl::new(
            control_min_schedule,
            0,
            1.,
            None,
            None,
            None,
            Default::default(),
            simulation_timestep,
        )
        .unwrap();

        let control_max = SetpointTimeControl::new(
            control_max_schedule,
            0,
            1.,
            None,
            None,
            None,
            Default::default(),
            simulation_timestep,
        )
        .unwrap();
        let immersion_heater = ImmersionHeater::new(
            rated_power,
            energy_supply_connection.clone(),
            simulation_timestep,
            Arc::new(Control::SetpointTime(control_min)),
            Arc::new(Control::SetpointTime(control_max)),
        );
        PositionedHeatSource {
            heat_source: Arc::new(Mutex::new(HeatSource::Storage(
                HeatSourceWithStorageTank::Immersion(Arc::new(Mutex::new(immersion_heater))),
            ))),
            heater_position,
            thermostat_position: Some(thermostat_position),
        }
    }

    #[fixture]
    fn storage_tank1(
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time_for_storage_tank: SimulationTime,
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions: Arc<ExternalConditions>,
    ) -> (StorageTank, Arc<RwLock<EnergySupply>>) {
        let control_min_schedule = vec![
            Some(52.),
            None,
            None,
            None,
            Some(52.),
            Some(52.),
            Some(52.),
            Some(52.),
        ];
        let control_max_schedule = vec![
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
            Some(55.),
        ];
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(
                FuelType::Electricity,
                simulation_time_for_storage_tank.total_steps(),
            )
            .build(),
        ));
        let energy_supply_connection =
            EnergySupply::connection(energy_supply.clone(), "immersion").unwrap();
        let heat_source = heat_source(
            simulation_time_for_storage_tank,
            energy_supply_connection.clone(),
            50.0,
            0.1,
            0.33,
            control_min_schedule,
            control_max_schedule,
        );

        let cold_feed = WaterSourceWithTemperature::ColdWaterSource(cold_water_source.clone());
        let simulation_timestep = simulation_time_for_storage_tank.step;

        let heat_sources = IndexMap::from([("imheater".to_string(), heat_source)]);
        let storage_tank = StorageTank::new(
            150.0,
            1.68,
            55.0,
            cold_feed,
            simulation_timestep,
            heat_sources,
            temp_internal_air_fn.clone(),
            external_conditions.clone(),
            None,
            None,
            Some(energy_supply_connection),
            *WATER,
            false,
        );

        (storage_tank, energy_supply)
    }

    #[fixture]
    fn storage_tank2(
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time_for_storage_tank: SimulationTime,
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions: Arc<ExternalConditions>,
    ) -> (StorageTank, Arc<RwLock<EnergySupply>>) {
        let control_min_schedule = vec![
            Some(52.),
            None,
            None,
            None,
            Some(52.),
            Some(52.),
            Some(52.),
            Some(52.),
        ];
        let control_max_schedule = vec![
            Some(60.),
            Some(60.),
            Some(60.),
            Some(60.),
            Some(60.),
            Some(60.),
            Some(60.),
            Some(60.),
        ];
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(
                FuelType::Electricity,
                simulation_time_for_storage_tank.total_steps(),
            )
            .build(),
        ));
        let energy_supply_connection =
            EnergySupply::connection(energy_supply.clone(), "immersion2").unwrap();
        let heat_source = heat_source(
            simulation_time_for_storage_tank,
            energy_supply_connection.clone(),
            5.0,
            0.6,
            0.6,
            control_min_schedule,
            control_max_schedule,
        );

        let cold_feed = WaterSourceWithTemperature::ColdWaterSource(cold_water_source.clone());
        let simulation_timestep = simulation_time_for_storage_tank.step;

        let heat_sources = IndexMap::from([("imheater2".to_string(), heat_source)]);
        let storage_tank = StorageTank::new(
            210.0,
            1.61,
            60.0,
            cold_feed,
            simulation_timestep,
            heat_sources,
            temp_internal_air_fn.clone(),
            external_conditions.clone(),
            None,
            None,
            Some(energy_supply_connection),
            *WATER,
            false,
        );

        (storage_tank, energy_supply)
    }

    #[fixture]
    fn external_conditions_for_pv_diverter(
        simulation_time_for_storage_tank: SimulationTime,
    ) -> Arc<ExternalConditions> {
        let air_temps = vec![0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0];
        let wind_speeds = vec![3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4];
        let wind_directions = vec![0.0; 8];
        let diffuse_horizontal_radiations = vec![333., 610., 572., 420., 0., 10., 90., 275.];
        let direct_beam_radiations = vec![420., 750., 425., 500., 0., 40., 0., 388.];
        let solar_reflectivity_of_ground = vec![0.2; 8760];
        let latitude = 51.42;
        let longitude = -0.75;
        let timezone = 0;
        let start_day = 0;
        let end_day = 0;
        let time_series_step = 1.;
        let january_first = 1;
        let leap_day_included = false;
        let direct_beam_conversion_needed = false;

        let shading_segments = vec![
            ShadingSegment {
                number: 1,
                start: 180.,
                end: 135.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 2,
                start: 135.,
                end: 90.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 3,
                start: 90.,
                end: 45.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 4,
                start: 45.,
                end: 0.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 5,
                start: 0.,
                end: -45.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 6,
                start: -45.,
                end: -90.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 7,
                start: -90.,
                end: -135.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 8,
                start: -135.,
                end: -180.,
                shading_objects: None,
                ..Default::default()
            },
        ]
        .into();

        Arc::new(ExternalConditions::new(
            &simulation_time_for_storage_tank.iter(),
            air_temps,
            wind_speeds,
            wind_directions,
            diffuse_horizontal_radiations,
            direct_beam_radiations,
            solar_reflectivity_of_ground,
            latitude,
            longitude,
            timezone,
            start_day,
            Some(end_day),
            time_series_step,
            Some(january_first),
            Some(DaylightSavingsConfig::NotApplicable),
            leap_day_included,
            direct_beam_conversion_needed,
            shading_segments,
        ))
    }

    #[fixture]
    fn storage_tank_for_pv_diverter(
        simulation_time_for_storage_tank: SimulationTime,
        immersion_heater: ImmersionHeater,
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions_for_pv_diverter: Arc<ExternalConditions>,
    ) -> StorageTank {
        let heater_position = 0.1;
        let thermostat_position = 0.33;
        let heat_source = PositionedHeatSource {
            heat_source: Arc::new(Mutex::new(HeatSource::Storage(
                HeatSourceWithStorageTank::Immersion(Arc::new(Mutex::new(immersion_heater))),
            ))),
            heater_position,
            thermostat_position: Some(thermostat_position),
        };
        let start_day = 0;
        let time_series_step = 1.;
        let cold_water_temps = vec![10.6, 11.0, 11.5, 12.1];
        let cold_feed = WaterSourceWithTemperature::ColdWaterSource(Arc::new(
            ColdWaterSource::new(cold_water_temps, start_day, time_series_step),
        ));
        let simulation_timestep = simulation_time_for_storage_tank.step;

        let heat_sources = IndexMap::from([("imheater".to_string(), heat_source)]);

        StorageTank::new(
            150.0,
            1.68,
            55.0,
            cold_feed,
            simulation_timestep,
            heat_sources,
            temp_internal_air_fn.clone(),
            external_conditions_for_pv_diverter.clone(),
            None,
            None,
            None,
            *WATER,
            false,
        )
    }

    #[rstest]
    fn test_divert_surplus(
        mut storage_tank_for_pv_diverter: StorageTank,
        immersion_heater: ImmersionHeater,
    ) {
        // _StorageTank__Q_ls_n_prev_heat_source is needed for the functions to
        // run the test but have no bearing in the results

        storage_tank_for_pv_diverter.q_ls_n_prev_heat_source = vec![0.0, 0.1, 0.2, 0.3];
        let pvdiverter = PVDiverter::new(
            &HotWaterStorageTank::StorageTank(Arc::new(RwLock::new(storage_tank_for_pv_diverter))),
            Arc::new(Mutex::new(immersion_heater)),
            "imheater".to_string(),
            None, // TODO (migration 0.34)
        );
        let sim_time = SimulationTime::new(0., 4., 1.);

        let supply_surplus = -1.0;
        assert_relative_eq!(
            pvdiverter
                .read()
                .divert_surplus(supply_surplus, sim_time.iter().current_iteration())
                .unwrap(),
            0.891553580246915
        );

        let supply_surplus = 0.0;
        assert_relative_eq!(
            pvdiverter
                .read()
                .divert_surplus(supply_surplus, sim_time.iter().current_iteration())
                .unwrap(),
            0.0
        );

        let supply_surplus = 1.0;
        assert_relative_eq!(
            pvdiverter
                .read()
                .divert_surplus(supply_surplus, sim_time.iter().current_iteration())
                .unwrap(),
            0.0
        );
    }

    #[fixture]
    fn simulation_time_for_solar_thermal() -> SimulationTime {
        SimulationTime::new(5088., 5112., 1.)
    }

    #[fixture]
    fn external_conditions_for_solar_thermal(
        simulation_time_for_solar_thermal: SimulationTime,
    ) -> Arc<ExternalConditions> {
        let air_temps = vec![
            19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0,
            19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0, 19.0,
        ];
        let wind_speeds = vec![
            3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1, 3.9, 3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1, 3.9,
            3.8, 3.9, 4.1, 3.8, 4.2, 4.3, 4.1,
        ];
        let wind_directions = vec![
            300.0, 250., 220., 180., 150., 120., 100., 80., 60., 40., 20., 10., 50., 100., 140.,
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
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 2,
                start: 135.,
                end: 90.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 3,
                start: 90.,
                end: 45.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 4,
                start: 45.,
                end: 0.,
                shading_objects: Some(vec![ShadingObject {
                    object_type: ShadingObjectType::Obstacle,
                    height: 10.5,
                    distance: 12.,
                }]),
                ..Default::default()
            },
            ShadingSegment {
                number: 5,
                start: 0.,
                end: -45.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 6,
                start: -45.,
                end: -90.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 7,
                start: -90.,
                end: -135.,
                shading_objects: None,
                ..Default::default()
            },
            ShadingSegment {
                number: 8,
                start: -135.,
                end: -180.,
                shading_objects: None,
                ..Default::default()
            },
        ]
        .into();

        Arc::new(ExternalConditions::new(
            &simulation_time_for_solar_thermal.iter(),
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
            Some(DaylightSavingsConfig::NotApplicable),
            false,
            false,
            shading_segments,
        ))
    }

    #[fixture]
    fn storage_tank_with_solar_thermal(
        external_conditions_for_solar_thermal: Arc<ExternalConditions>,
        temp_internal_air_fn: TempInternalAirFn,
        simulation_time_for_solar_thermal: SimulationTime,
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
        let cold_feed = WaterSourceWithTemperature::ColdWaterSource(Arc::new(
            ColdWaterSource::new(cold_water_temps.to_vec(), 212, 1.),
        ));
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(
                FuelType::Electricity,
                simulation_time_for_solar_thermal.total_steps(),
            )
            .build(),
        ));
        let energy_supply_conn =
            EnergySupply::connection(energy_supply.clone(), "solarthermal").unwrap();
        let control_max = SetpointTimeControl::new(
            vec![
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
                Some(55.),
            ],
            212,
            1.,
            None,
            None,
            None,
            Default::default(),
            simulation_time_for_solar_thermal.step,
        );

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
            temp_internal_air_fn.clone(),
            simulation_time_for_solar_thermal.step,
            Arc::new(Control::SetpointTime(control_max.unwrap())),
            *WATER,
            None, // TODO (migration 0.34) check if this is correct
        )));

        let storage_tank = StorageTank::new(
            150.0,
            1.68,
            55.0,
            cold_feed,
            simulation_time_for_solar_thermal.step,
            IndexMap::from([(
                "solthermal".to_string(),
                PositionedHeatSource {
                    heat_source: Arc::new(Mutex::new(HeatSource::Storage(
                        HeatSourceWithStorageTank::Solar(solar_thermal.clone()),
                    ))),
                    heater_position: 0.1,
                    thermostat_position: Some(0.33),
                },
            )]),
            temp_internal_air_fn,
            external_conditions_for_solar_thermal,
            None,
            None,
            None,
            *WATER,
            false,
        );

        (
            storage_tank,
            solar_thermal,
            simulation_time_for_solar_thermal,
            energy_supply,
        )
    }

    #[rstest]
    fn test_demand_hot_water(
        simulation_time_for_storage_tank: SimulationTime,
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        storage_tank2: (StorageTank, Arc<RwLock<EnergySupply>>),
    ) {
        let (mut storage_tank1, energy_supply1) = storage_tank1;
        let (mut storage_tank2, energy_supply2) = storage_tank2;
        let usage_events = [
            vec![
                TypedScheduleEvent {
                    start: 6.,
                    duration: Some(6.),
                    temperature: 41.0,
                    name: "IES".to_string(),
                    event_type: WaterScheduleEventType::Shower,
                    volume: None,
                    warm_volume: None,
                    pipework_volume: None,
                },
                TypedScheduleEvent {
                    start: 6.,
                    duration: Some(6.),
                    temperature: 41.0,
                    name: "mixer".to_string(),
                    event_type: WaterScheduleEventType::Shower,
                    volume: None,
                    warm_volume: Some(48.),
                    pipework_volume: None,
                },
                TypedScheduleEvent {
                    start: 6.,
                    duration: Some(20.),
                    temperature: 43.0,
                    name: "medium".to_string(),
                    event_type: WaterScheduleEventType::Bath,
                    volume: None,
                    warm_volume: Some(100.),
                    pipework_volume: None,
                },
                TypedScheduleEvent {
                    start: 6.,
                    duration: Some(1.),
                    temperature: 40.0,
                    name: "other".to_string(),
                    event_type: WaterScheduleEventType::Other,
                    volume: None,
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
                volume: None,
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
                volume: None,
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
                volume: None,
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

        // Also test case where heater does not heat all layers, to ensure this is handled correctly

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
            let usage_events_for_iteration = Some(usage_events[t_idx].clone());
            storage_tank1
                .demand_hot_water(usage_events_for_iteration.clone(), t_it)
                .unwrap();

            // Verify the temperatures against expected results
            assert_eq!(
                storage_tank1.temp_n, expected_temperatures_1[t_idx],
                "incorrect temperatures returned"
            );

            assert_relative_eq!(
                energy_supply1.read().results_by_end_user()["immersion"][t_idx],
                expected_energy_supplied_1[t_idx],
                max_relative = 1e-6
            );

            storage_tank2
                .demand_hot_water(usage_events_for_iteration.clone(), t_it)
                .unwrap();

            assert_eq!(
                storage_tank2.temp_n, expected_temperatures_2[t_idx],
                "incorrect temperatures returned"
            );

            assert_relative_eq!(
                energy_supply2.read().results_by_end_user()["immersion2"][t_idx],
                expected_energy_supplied_2[t_idx],
                max_relative = 1e-6
            );
        }
    }

    #[rstest]
    fn test_temp_surrounding_primary_pipework(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        let (storage_tank1, _) = storage_tank1;
        // External Pipe
        let pipework = Pipework::new(
            PipeworkLocation::External,
            0.025,
            0.027,
            1.0,
            0.035,
            0.038,
            false,
            WaterPipeContentsType::Water,
        )
        .unwrap();
        for (t_idx, t_it) in simulation_time_for_storage_tank.iter().enumerate() {
            assert_eq!(
                storage_tank1.temp_surrounding_primary_pipework(&pipework, t_it),
                [0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0][t_idx]
            )
        }
        // Internal Pipe
        let pipework = Pipework::new(
            PipeworkLocation::Internal,
            0.025,
            0.027,
            1.0,
            0.035,
            0.038,
            false,
            WaterPipeContentsType::Water,
        )
        .unwrap();
        for (t_idx, t_it) in simulation_time_for_storage_tank.iter().enumerate() {
            assert_eq!(
                storage_tank1.temp_surrounding_primary_pipework(&pipework, t_it),
                [20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0][t_idx]
            )
        }
    }

    #[rstest]
    fn test_get_cold_water_source(storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>)) {
        let (storage_tank1, _) = storage_tank1;
        let result = storage_tank1.get_cold_water_source();

        assert!(matches!(
            result,
            WaterSourceWithTemperature::ColdWaterSource(_)
        ));
    }

    #[rstest]
    fn test_get_temp_hot_water(storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>)) {
        let (storage_tank1, _) = storage_tank1;

        assert_eq!(storage_tank1.get_temp_hot_water(), 55.0);
    }

    #[rstest]
    fn test_stand_by_losses_coefficient(storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>)) {
        let (storage_tank1, _) = storage_tank1;

        assert_relative_eq!(
            storage_tank1.stand_by_losses_coefficient(),
            1.5555555555555556
        );
    }

    #[rstest]
    fn test_potential_energy_input(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        storage_tank_with_solar_thermal: (
            StorageTank,
            Arc<Mutex<SolarThermalSystem>>,
            SimulationTime,
            Arc<RwLock<EnergySupply>>,
        ),
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        // ImmersionHeater as heat source
        let (mut storage_tank1, _) = storage_tank1;
        let temp_s3_n = [55.0, 55.0, 55.0, 55.0, 55.0, 55.0, 55.0, 55.0];
        let heat_source = storage_tank1.heat_source_data["imheater"]
            .clone()
            .heat_source;

        assert_eq!(
            storage_tank1
                .potential_energy_input(
                    &temp_s3_n,
                    heat_source,
                    "imheater",
                    0,
                    7,
                    simulation_time_for_storage_tank.iter().current_iteration()
                )
                .unwrap(),
            [0.0, 0., 0., 0.]
        );

        // SolarThermal as heat source
        let (mut storage_tank_solar_thermal, _, simtime, _) = storage_tank_with_solar_thermal;
        let temp_s3_n = [
            25.0, 15.0, 35.0, 45.0, 55.0, 50.0, 30.0, 20.0, 25.0, 15.0, 35.0, 45.0, 55.0, 50.0,
            30.0, 20.0, 25.0, 15.0, 35.0, 45.0, 55.0, 50.0, 30.0, 20.0, 25.0, 15.0, 35.0, 45.0,
            55.0, 50.0, 30.0, 20.0,
        ];
        let heat_source = storage_tank_solar_thermal.heat_source_data["solthermal"]
            .clone()
            .heat_source;

        for (t_idx, t_it) in simtime.iter().enumerate() {
            let actual_result = storage_tank_solar_thermal
                .potential_energy_input(&temp_s3_n, heat_source.clone(), "solthermal", 0, 7, t_it)
                .unwrap();
            let expected_result = [
                [0., 0., 0., 0.],
                [0., 0., 0., 0.],
                [0., 0., 0., 0.],
                [0., 0., 0., 0.],
                [0., 0., 0., 0.],
                [0., 0., 0., 0.],
                [0., 0., 0., 0.],
                [0., 0., 0., 0.],
                [0.47214338269526945, 0., 0., 0.],
                [0.794165996101526, 0., 0., 0.],
                [1.2488375719961642, 0., 0., 0.],
                [1.0218936489635675, 0., 0., 0.],
                [1.1483985152150102, 0., 0., 0.],
                [1.5175839864027383, 0., 0., 0.],
                [0.9602170463493307, 0., 0., 0.],
                [0.5981490998786696, 0., 0., 0.],
                [0.3454397002046902, 0., 0., 0.],
                [0., 0., 0., 0.],
                [0., 0., 0., 0.],
                [0., 0., 0., 0.],
                [0., 0., 0., 0.],
                [0., 0., 0., 0.],
                [0., 0., 0., 0.],
                [0., 0., 0., 0.],
            ][t_idx];

            // Compare each element using assert_relative_eq
            for (expected_value, actual_value) in expected_result.iter().zip(actual_result) {
                assert_relative_eq!(expected_value, &actual_value, max_relative = 1e-13);
            }
        }
    }

    #[rstest]
    fn test_storage_tank_potential_effect(storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>)) {
        let (storage_tank1, _) = storage_tank1;
        let energy_proposed = 0.;
        let temp_s3_n = [25.0, 15.0, 35.0, 45.0, 55.0, 50.0, 30.0, 20.0];
        assert_eq!(
            storage_tank1.storage_tank_potential_effect(energy_proposed, &temp_s3_n),
            (20.0, 45.0)
        );
    }

    #[rstest]
    fn test_energy_input(storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>)) {
        let (storage_tank1, _) = storage_tank1;
        let temp_s3_n = [25.0, 15.0, 35.0, 45.0, 55.0, 50.0, 30.0, 20.0];
        let q_x_in_n = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7];

        let result = storage_tank1.energy_input(&temp_s3_n, &q_x_in_n);

        assert_eq!(
            result,
            (
                5.83,
                vec![
                    25.0,
                    17.294455066921607,
                    39.588910133843214,
                    51.883365200764814
                ]
            )
        );
    }

    #[rstest]
    fn test_rearrange_temperatures(storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>)) {
        let (storage_tank1, _) = storage_tank1;
        let temp_s6_n = [2.5, 3.7, 10.36, 17.43, 32.95, 35.91, 35.91, 42.2];
        assert_eq!(
            storage_tank1.rearrange_temperatures(&temp_s6_n),
            (
                vec![
                    0.10895833333333334,
                    0.16125833333333334,
                    0.45152333333333333,
                    0.7596575
                ],
                vec![2.5, 3.7, 10.36, 17.43, 32.95, 35.91, 35.91, 42.2]
            )
        )
    }

    #[rstest]
    fn test_thermal_losses(storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>)) {
        let (mut storage_tank1, _) = storage_tank1;
        let temp_s7_n = [12.0, 18.0, 25.0, 32.0, 37.0, 45.0, 49.0, 58.0];
        let q_x_in_n = [0., 1., 2., 3., 4., 5., 6., 7., 8.];
        let q_h_sto_s7 = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8];
        let heater_layer = 2;
        let q_ls_n_prev_heat_source = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let setpntmax = 55.0;

        assert_eq!(
            storage_tank1.thermal_losses(
                &temp_s7_n,
                &q_x_in_n,
                q_h_sto_s7,
                heater_layer,
                &q_ls_n_prev_heat_source,
                setpntmax
            ),
            (
                36.0,
                0.012203333333333333,
                vec![
                    12.0,
                    17.97925925925926,
                    24.906666666666666,
                    31.834074074074074
                ],
                vec![
                    0.0,
                    0.0009039506172839507,
                    0.004067777777777778,
                    0.007231604938271605
                ]
            )
        );
    }

    #[rstest]
    fn test_run_heat_sources(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        let (mut storage_tank1, _) = storage_tank1;
        let temp_s3_n = vec![5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0];
        let heat_source = storage_tank1.heat_source_data["imheater"]
            .clone()
            .heat_source;
        let heater_layer = 2;
        let thermostat_layer = 7;
        let q_ls_prev_heat_source = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert_eq!(
            storage_tank1
                .run_heat_sources(
                    temp_s3_n,
                    heat_source,
                    "imheater",
                    heater_layer,
                    thermostat_layer,
                    &q_ls_prev_heat_source,
                    simulation_time_for_storage_tank.iter().current_iteration()
                )
                .unwrap(),
            TemperatureCalculation {
                temp_s8_n: vec![5.0, 10.0, 55.0, 55.0],
                q_x_in_n: vec![0., 0., 50.0, 0.],
                q_s6: 52.17916666666666,
                temp_s6_n: vec![5.0, 10.0, 1162.227533460803, 20.0],
                temp_s7_n: vec![5.0, 10.0, 591.1137667304015, 591.1137667304015],
                q_in_h_w: 3.3040040740740793,
                q_ls: 0.03525407407407408,
                q_ls_n: vec![0.0, 0.0, 0.01762703703703704, 0.01762703703703704]
            }
        );
    }

    #[rstest]
    fn test_calculate_temperatures(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        let (mut storage_tank1, _) = storage_tank1;
        let temp_s3_n = vec![10., 15., 20., 25., 25., 30., 35., 50.];
        let heat_source = storage_tank1.heat_source_data["imheater"]
            .clone()
            .heat_source;
        let q_x_in_n = vec![0., 1., 2., 3., 4., 5., 6., 7., 8.];
        let heater_layer = 2;
        let q_ls_n_prev_heat_source = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

        assert_eq!(
            storage_tank1
                .calculate_temperatures(
                    &temp_s3_n,
                    heat_source,
                    q_x_in_n,
                    heater_layer,
                    &q_ls_n_prev_heat_source,
                    simulation_time_for_storage_tank.iter().current_iteration()
                )
                .unwrap(),
            TemperatureCalculation {
                temp_s8_n: vec![10.0, 37.71697755116492, 55.0, 55.0],
                q_x_in_n: vec![0., 1., 2., 3., 4., 5., 6., 7., 8.],
                q_s6: 9.050833333333333,
                temp_s6_n: vec![
                    10.0,
                    37.944550669216056,
                    65.88910133843211,
                    93.83365200764818
                ],
                temp_s7_n: vec![
                    10.0,
                    37.944550669216056,
                    65.88910133843211,
                    93.83365200764818
                ],
                q_in_h_w: 33.86817074074074,
                q_ls: 0.04517246913580247,
                q_ls_n: vec![
                    0.0,
                    0.009918395061728393,
                    0.01762703703703704,
                    0.01762703703703704
                ]
            }
        )
    }

    #[rstest]
    fn test_allocate_hot_water(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        let (mut storage_tank1, _) = storage_tank1;
        storage_tank1.temp_average_drawoff_volweighted = Some(0.0);
        storage_tank1.total_volume_drawoff = Some(0.0);
        let event = TypedScheduleEvent {
            start: 0.,
            duration: Some(1.),
            temperature: 41.0,
            event_type: WaterScheduleEventType::Other,
            name: "other".to_string(),
            volume: None,
            warm_volume: Some(8.0),
            pipework_volume: Some(5.),
        };
        assert_eq!(
            storage_tank1.allocate_hot_water(
                event,
                simulation_time_for_storage_tank.iter().current_iteration()
            ),
            (
                10.51111111111112,
                vec![42.39, 55.0, 55.0, 55.0],
                0.5497311111111112,
                0.0,
                false
            )
        );

        // With no pipework vol
        let event = TypedScheduleEvent {
            start: 0.,
            duration: Some(1.),
            temperature: 41.0,
            event_type: WaterScheduleEventType::Other,
            name: "other".to_string(),
            volume: None,
            warm_volume: Some(8.0),
            pipework_volume: Some(0.),
        };
        assert_eq!(
            storage_tank1.allocate_hot_water(
                event,
                simulation_time_for_storage_tank.iter().current_iteration()
            ),
            (
                5.51111111111112,
                vec![48.39, 55.0, 55.0, 55.0],
                0.2882311111111111,
                0.0,
                false
            )
        );
    }

    #[rstest]
    fn test_calculate_new_temperatures(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        let (storage_tank1, _) = storage_tank1;
        let remaining_vol = vec![0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4.];

        assert_eq!(
            storage_tank1.calculate_new_temperatures(
                remaining_vol,
                simulation_time_for_storage_tank.iter().current_iteration()
            ),
            (vec![10.0, 10.0, 10.0, 16.0], false)
        );
    }

    #[rstest]
    fn test_additional_energy_input(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        storage_tank2: (StorageTank, Arc<RwLock<EnergySupply>>),
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        let (mut storage_tank1, _) = storage_tank1;
        let heat_source = storage_tank1.heat_source_data["imheater"]
            .clone()
            .heat_source;
        let energy_input = 5.0;
        storage_tank1.q_ls_n_prev_heat_source = vec![0.0, 0.1, 0.2, 0.3];
        assert_eq!(
            storage_tank1
                .additional_energy_input(
                    heat_source,
                    "imheater",
                    energy_input,
                    simulation_time_for_storage_tank.iter().current_iteration()
                )
                .unwrap(),
            0.01762703703703572
        );

        // Test with no energy input
        let (mut storage_tank2, _) = storage_tank2;
        let heat_source = storage_tank2.heat_source_data["imheater2"]
            .clone()
            .heat_source;
        let energy_input = 0.;
        storage_tank2.q_ls_n_prev_heat_source = vec![0.0, 0.1, 0.2, 0.3];
        assert_eq!(
            storage_tank2
                .additional_energy_input(
                    heat_source,
                    "imheater2",
                    energy_input,
                    simulation_time_for_storage_tank.iter().current_iteration()
                )
                .unwrap(),
            0.0
        );
    }

    #[rstest]
    fn test_internal_gains(storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>)) {
        let (mut storage_tank1, _) = storage_tank1;
        storage_tank1.q_sto_h_ls_rbl = Some(0.05);

        assert_eq!(storage_tank1.internal_gains(), 50.);
    }

    #[rstest]
    fn test_primary_pipework_losses(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        let (mut storage_tank1, _) = storage_tank1;
        let input_energy_adj = 0.0;
        let setpnt_max = 55.0;
        let primary_pipework_lst = vec![
            Pipework::new(
                PipeworkLocation::Internal,
                0.024,
                0.027,
                2.0,
                0.035,
                0.04,
                false,
                WaterPipeContentsType::Water,
            )
            .unwrap(),
            Pipework::new(
                PipeworkLocation::External,
                0.025,
                0.027,
                0.0,
                0.035,
                0.038,
                false,
                WaterPipeContentsType::Water,
            )
            .unwrap(),
        ];

        storage_tank1.primary_pipework_lst = Some(primary_pipework_lst);

        for (t_idx, t_it) in simulation_time_for_storage_tank.iter().enumerate() {
            assert_eq!(
                storage_tank1.primary_pipework_losses(input_energy_adj, setpnt_max, t_it),
                [
                    (0.0, 0.0),
                    (0.0, 0.0),
                    (0.0, 0.0),
                    (0.0, 0.0),
                    (0.0, 0.0),
                    (0.0, 0.0),
                    (0.0, 0.0),
                    (0.0, 0.0)
                ][t_idx]
            )
        }

        // With value for input_energy_adj
        let input_energy_adj = 3.;

        for (t_idx, t_it) in simulation_time_for_storage_tank.iter().enumerate() {
            assert_eq!(
                storage_tank1.primary_pipework_losses(input_energy_adj, setpnt_max, t_it),
                [
                    (0.04665869119015863, 9.854304934823482),
                    (0.04665869119015863, 9.854304934823482),
                    (0.04665869119015863, 9.854304934823482),
                    (0.04665869119015863, 9.854304934823482),
                    (0.04665869119015863, 9.854304934823482),
                    (0.04665869119015863, 9.854304934823482),
                    (0.04665869119015863, 9.854304934823482),
                    (0.04665869119015863, 9.854304934823482)
                ][t_idx]
            )
        }
    }

    #[fixture]
    fn simulation_time_for_immersion_heater() -> SimulationTime {
        SimulationTime::new(0., 4., 1.)
    }

    #[fixture]
    fn immersion_heater(simulation_time_for_immersion_heater: SimulationTime) -> ImmersionHeater {
        let rated_power = 50.;
        let energy_supply = EnergySupplyBuilder::new(
            FuelType::MainsGas,
            simulation_time_for_immersion_heater.total_steps(),
        )
        .build();
        let energy_supply_connection =
            EnergySupply::connection(Arc::new(RwLock::new(energy_supply)), "shower").unwrap();
        let timestep = simulation_time_for_immersion_heater.step;

        let control_min = Arc::new(Control::SetpointTime(
            SetpointTimeControl::new(
                vec![Some(52.), Some(52.), None, Some(52.)],
                0,
                1.,
                None,
                None,
                None,
                Default::default(),
                timestep,
            )
            .unwrap(),
        ));

        let control_max = Arc::new(Control::SetpointTime(
            SetpointTimeControl::new(
                vec![Some(60.), Some(60.), Some(60.), Some(60.)],
                0,
                1.,
                None,
                None,
                None,
                Default::default(),
                timestep,
            )
            .unwrap(),
        ));

        ImmersionHeater::new(
            rated_power,
            energy_supply_connection,
            timestep,
            control_min,
            control_max,
        )
    }

    #[rstest]
    fn test_demand_energy_for_immersion_heater(
        mut immersion_heater: ImmersionHeater,
        simulation_time_for_immersion_heater: SimulationTime,
    ) {
        let energy_inputs = [40., 100., 30., 20.];
        let expected_energy = [40., 50., 0., 20.];
        for (t_idx, t_it) in simulation_time_for_immersion_heater.iter().enumerate() {
            assert_eq!(
                immersion_heater
                    .demand_energy(energy_inputs[t_idx], t_it)
                    .unwrap(),
                expected_energy[t_idx],
                "incorrect energy demand calculated"
            );
        }
        assert!(immersion_heater
            .demand_energy(
                -1.,
                simulation_time_for_immersion_heater
                    .iter()
                    .current_iteration()
            )
            .is_err());
    }

    #[rstest]
    fn test_energy_output_max_for_immersion_heater(
        immersion_heater: ImmersionHeater,
        simulation_time_for_immersion_heater: SimulationTime,
    ) {
        for t_it in simulation_time_for_immersion_heater.iter() {
            assert_eq!(
                immersion_heater.energy_output_max(t_it, true), // In Python another parameter (return_temp = 55.0) is passed in to energy_output_max but never used so we have skipped this in Rust
                50.,
                "incorrect energy output max calculated"
            );
        }

        for (t_idx, t_it) in simulation_time_for_immersion_heater.iter().enumerate() {
            assert_eq!(
                immersion_heater.energy_output_max(t_it, false), // In Python another parameter (return_temp = 40.0) is passed in to energy_output_max but never used so we have skipped this in Rust
                [50.0, 50.0, 0.0, 50.0][t_idx],
                "incorrect energy output max calculated"
            );
        }
    }

    #[rstest]
    fn test_energy_output_max_with_solar_thermal(
        storage_tank_with_solar_thermal: (
            StorageTank,
            Arc<Mutex<SolarThermalSystem>>,
            SimulationTime,
            Arc<RwLock<EnergySupply>>,
        ),
    ) {
        let temp_storage_tank_s3_n = [
            17.2, 17.2, 17.2, 17.2, 17.43, 32.95, 35.91, 35.91, 35.91, 42.25, 43.46, 43.46, 43.46,
            43.46, 43.46, 43.46, 43.46, 43.46, 43.46, 43.46, 43.46, 43.46, 43.46, 43.46,
        ];

        let (storage_tank_solar_thermal, solar_thermal, simulation_time, _) =
            storage_tank_with_solar_thermal;

        let expected = [
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.5441009409757523,
            0.7375096332994749,
            1.043768276267769,
            1.4751675743361337,
            1.2419712751344847,
            1.3717388178387737,
            1.7256641787433769,
            1.1745201463530213,
            0.8403815010507395,
            0.6056200608960752,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
        ];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let actual = solar_thermal.lock().energy_output_max(
                &storage_tank_solar_thermal,
                &temp_storage_tank_s3_n,
                &t_it,
            );
            assert_relative_eq!(actual, expected[t_idx], max_relative = 1e-7);
        }

        solar_thermal.lock().sol_loc = SolarCellLocation::Nhs;
        let expected = [
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.5443432944582923,
            0.7377445749305638,
            1.044003444574131,
            1.4754027357101285,
            1.2422064367204904,
            1.371973979418296,
            1.7258993403230973,
            1.1747553079327355,
            0.8406166626304539,
            0.6058552224757895,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
        ];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let actual = solar_thermal.lock().energy_output_max(
                &storage_tank_solar_thermal,
                &temp_storage_tank_s3_n,
                &t_it,
            );
            assert_relative_eq!(actual, expected[t_idx], max_relative = 1e-7);
        }

        solar_thermal.lock().sol_loc = SolarCellLocation::Hs;
        let expected = [
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.5445856479408322,
            0.7379795165616525,
            1.044238612880493,
            1.475637897084123,
            1.2424415983064963,
            1.372209140997818,
            1.7261345019028178,
            1.1749904695124498,
            0.8408518242101682,
            0.6060903840555039,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
        ];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let actual = solar_thermal.lock().energy_output_max(
                &storage_tank_solar_thermal,
                &temp_storage_tank_s3_n,
                &t_it,
            );
            assert_relative_eq!(actual, expected[t_idx], max_relative = 1e-7);
        }
    }

    #[rstest]
    fn test_demand_energy_with_solar_thermal(
        #[from(storage_tank_with_solar_thermal)]
        (storage_tank_solar_thermal, solar_thermal, simulation_time, _): (
            StorageTank,
            Arc<Mutex<SolarThermalSystem>>,
            SimulationTime,
            Arc<RwLock<EnergySupply>>,
        ),
    ) {
        let temp_storage_tank_s3_n = [
            17.2, 17.2, 17.2, 17.2, 17.43, 32.95, 35.91, 35.91, 35.91, 42.25, 43.46, 43.46, 43.46,
            43.46, 43.46, 43.46, 43.46, 43.46, 43.46, 43.46, 43.46, 43.46, 43.46, 43.46,
        ];

        solar_thermal.lock().energy_output_max(
            &storage_tank_solar_thermal,
            &temp_storage_tank_s3_n,
            &simulation_time.iter().current_iteration(),
        );

        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_eq!(solar_thermal.lock().demand_energy(100., t_idx), 0.);
        }
    }

    #[rstest]
    // in Python this test is called test_demand_hot_water and is from test_storage_tank_with_solar_thermal.py
    fn test_demand_hot_water_for_storage_tank_with_solar_thermal(
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
                volume: None,
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
                volume: None,
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
                volume: None,
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
                volume: None,
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
                volume: None,
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
                volume: None,
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
                volume: None,
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
                volume: None,
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
                volume: None,
                warm_volume: Some(52.0),
                pipework_volume: None,
            }],
            vec![],
            vec![],
        ];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            storage_tank
                .demand_hot_water(Some(usage_events.get(t_idx).unwrap().clone()), t_it)
                .unwrap();
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
