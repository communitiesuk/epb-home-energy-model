use crate::compare_floats::{max_of_2, min_of_2};
use crate::core::common::WaterSourceWithTemperature;
use crate::core::controls::time_control::{Control, ControlBehaviour};
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::material_properties::{MaterialProperties, WATER};
use crate::core::pipework::{Pipework, PipeworkLocation, Pipeworkesque};
use crate::core::schedule::TypedScheduleEvent;
use crate::core::units::{MINUTES_PER_HOUR, WATTS_PER_KILOWATT};
use crate::core::water_heat_demand::misc::frac_hot_water;
use crate::corpus::{HeatSource, TempInternalAirFn};
use crate::external_conditions::ExternalConditions;
use crate::input::{SolarCollectorLoopLocation, WaterPipework};
use crate::simulation_time::SimulationTimeIteration;
use crate::statistics::np_interp;
use anyhow::{anyhow, bail};
use arc_swap::ArcSwapOption;
use atomic_float::AtomicF64;
use derivative::Derivative;
use indexmap::IndexMap;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use parking_lot::{Mutex, RwLock};
use smartstring::alias::String;
use std::collections::HashMap;
use std::iter;
use std::ops::Deref;
use std::sync::atomic::{AtomicBool, Ordering};
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

// TODO (from Python) - link to zone temp at timestep possibly and location of tank (in or out of heated space)
const DEFAULT_AMBIENT_TEMPERATURE: f64 = 16.;

// Primary pipework gains for the timestep
const DEFAULT_PIPEWORK_PRIMARY_GAINS_FOR_TIMESTEP: f64 = 0.;

// Time of finalisation of the previous hot water event
const DEFAULT_PREVIOUS_EVENT_TIME_END: f64 = 0.;

// Auxiliary energy recovery factor
const THERMAL_CONSTANTS_F_RVD_AUX: f64 = 0.25;

// Thermal loss recovery factor
const THERMAL_CONSTANTS_F_STO_M: f64 = 0.75;

// Standby losses adaptation factor
const THERMAL_CONSTANTS_F_STO_BAC_ACC: f64 = 1.;

// utility method to check if an array is sorted
fn is_sorted(vec: &[f64]) -> bool {
    vec.windows(2).all(|w| w[0] <= w[1])
}

// utility method for rounding
fn round_by_precision(src: f64, precision: f64) -> f64 {
    (precision * src).round() / precision
}

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
    initial_temperature: f64,
    q_std_ls_ref: f64, // measured standby losses due to cylinder insulation at standardised conditions, in kWh/24h
    cold_feed: WaterSourceWithTemperature,
    simulation_timestep: f64,
    number_of_volumes: usize,
    temp_flow_prev: AtomicF64,
    #[derivative(Debug = "ignore")]
    temp_internal_air_fn: TempInternalAirFn,
    external_conditions: Arc<ExternalConditions>,
    volume_total_in_litres: f64,
    vol_n: Vec<f64>,
    cp: f64,  // contents (usually water) specific heat in kWh/kg.K
    rho: f64, // volumic mass in kg/litre
    temp_n: Arc<RwLock<Vec<f64>>>,
    input_energy_adj_prev_timestep: AtomicF64,
    primary_pipework: Option<Vec<Pipework>>,
    primary_pipework_losses_kwh: AtomicF64,
    storage_losses_kwh: AtomicF64,
    heat_source_data: IndexMap<String, PositionedHeatSource>, // heat sources, sorted by heater position
    heating_active: HashMap<String, AtomicBool>,
    q_ls_n_prev_heat_source: Arc<RwLock<Vec<f64>>>,
    q_sto_h_ls_rbl: AtomicF64, // total recoverable heat losses for heating in kWh, memoised between steps
    pipework_primary_gains_for_timestep: AtomicF64, // primary pipework gains for a timestep (mutates over lifetime)
    #[cfg(test)]
    energy_demand_test: AtomicF64,
    temp_final_drawoff: AtomicF64, // In Python this is created from inside extract_hot_water()
    temp_average_drawoff: AtomicF64, // In Python this is created from inside extract_hot_water()
    temp_average_drawoff_volweighted: AtomicF64, // In Python this is created from inside extract_hot_water()
    total_volume_drawoff: AtomicF64, // In Python this is created from inside extract_hot_water()
    ambient_temperature: f64,        // TODO should these be AtomicF64
    previous_event_time_end: f64,
}

#[derive(Debug)]
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
    /// *  `number_of_volumes` -number of volumes the storage is modelled with
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
        initial_temperature: f64,
        cold_feed: WaterSourceWithTemperature,
        simulation_timestep: f64,
        heat_sources: IndexMap<String, PositionedHeatSource>,
        // In Python this is "project" but only temp_internal_air is accessed from it
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions: Arc<ExternalConditions>,
        number_of_volumes: Option<usize>,
        primary_pipework_lst: Option<&Vec<WaterPipework>>,
        contents: MaterialProperties,
        ambient_temperature: Option<f64>,
        pipework_primary_gains_for_timestep: Option<f64>, // TODO check we need this
        previous_event_time_end: Option<f64>,
        _detailed_output_heating_cooling: bool, // TODO implement logic for this to match Python 0.32
    ) -> anyhow::Result<Self> {
        let q_std_ls_ref = losses;
        let ambient_temperature = ambient_temperature.unwrap_or(DEFAULT_AMBIENT_TEMPERATURE);
        let pipework_primary_gains_for_timestep = pipework_primary_gains_for_timestep
            .unwrap_or(DEFAULT_PIPEWORK_PRIMARY_GAINS_FOR_TIMESTEP);
        let previous_event_time_end =
            previous_event_time_end.unwrap_or(DEFAULT_PREVIOUS_EVENT_TIME_END);

        let volume_total_in_litres = volume;
        let number_of_volumes = number_of_volumes.unwrap_or(4);
        // list of volume of layers in litres
        let vol_n = iter::repeat_n(
            volume_total_in_litres / number_of_volumes as f64,
            number_of_volumes,
        )
        .collect_vec();
        // water specific heat in kWh/kg.K
        let cp = contents.specific_heat_capacity_kwh();
        let rho = contents.density();

        // 6.4.3.2 STEP 0 Initialization
        let temp_n = Arc::new(RwLock::new(vec![initial_temperature; number_of_volumes]));

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
                    .map(|pipework| pipework.to_owned().try_into().map_err(anyhow::Error::msg))
                    .collect::<anyhow::Result<Vec<Pipework>>>()?
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

        let heating_active = heat_sources
            .iter()
            .map(|(name, _heat_source)| ((*name).clone(), false.into()))
            .collect();

        Ok(Self {
            initial_temperature,
            q_std_ls_ref,
            cold_feed,
            simulation_timestep,
            number_of_volumes,
            temp_flow_prev: Default::default(),
            temp_internal_air_fn,
            external_conditions,
            volume_total_in_litres,
            vol_n,
            cp,
            rho,
            temp_n,
            input_energy_adj_prev_timestep: input_energy_adj_prev_timestep.into(),
            primary_pipework: primary_pipework_lst,
            primary_pipework_losses_kwh: primary_pipework_losses_kwh.into(),
            storage_losses_kwh: storage_losses_kwh.into(),
            heat_source_data,
            heating_active,
            q_ls_n_prev_heat_source: Default::default(),
            q_sto_h_ls_rbl: Default::default(),
            pipework_primary_gains_for_timestep: pipework_primary_gains_for_timestep.into(),
            #[cfg(test)]
            energy_demand_test: energy_demand_test.into(),
            temp_final_drawoff: Default::default(),
            temp_average_drawoff: Default::default(),
            temp_average_drawoff_volweighted: Default::default(),
            total_volume_drawoff: Default::default(),
            ambient_temperature,
            previous_event_time_end,
        })
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
        &self,
        usage_events: Option<Vec<TypedScheduleEvent>>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let mut q_use_w = 0.;
        let mut _volume_demanded = 0.;

        let mut temp_s3_n = self.temp_n.read().clone();

        self.temp_average_drawoff_volweighted
            .store(0., Ordering::SeqCst);
        self.temp_final_drawoff.store(0., Ordering::SeqCst);
        self.total_volume_drawoff.store(0., Ordering::SeqCst);
        self.temp_average_drawoff
            .store(self.initial_temperature, Ordering::SeqCst);

        for event in usage_events.iter().flatten() {
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
            // self.previous_event_time_end = deepcopy(time_start_current_event + (event['duration'] + 0.0) / 60.0)

            let (volume_used, energy_withdrawn, remaining_vols) =
                self.extract_hot_water(event.clone(), simtime)?;

            // Determine the new temperature distribution after displacement
            // Now that pre-heated sources can be the 'cold' feed, rearrangement of temperaturs, that used to
            // only happen before after the input from heat sources, could be required after the displacement
            // of water bringing new water from the 'cold' feed that could be warmer than the existing one.
            // flag is calculated for that purpose.
            let (temp_s3_n_new, rearrange) =
                self.calc_temps_after_extraction(remaining_vols, simtime);
            temp_s3_n = temp_s3_n_new;

            if rearrange {
                // Re-arrange the temperatures in the storage after energy input from pre-heated tank
                temp_s3_n = self.rearrange_temperatures(&temp_s3_n).1
            }

            *self.temp_n.write() = temp_s3_n.clone();

            _volume_demanded += volume_used;
            q_use_w += energy_withdrawn;
        }

        self.temp_average_drawoff.store(
            match self.total_volume_drawoff.load(Ordering::SeqCst) {
                value if value != 0. => {
                    let temp_average_drawoff_volweighted =
                        self.temp_average_drawoff_volweighted.load(Ordering::SeqCst);
                    temp_average_drawoff_volweighted / value
                }
                _ => temp_s3_n
                    .last()
                    .copied()
                    .ok_or_else(|| anyhow!("temp_s3_n was unexpectedly empty"))?,
            },
            Ordering::SeqCst,
        );

        // TODO (from Python) 6.4.3.6 STEP 4 Volume to be withdrawn from the storage (for Heating)
        // TODO (from Python) - 6.4.3.7 STEP 5 Temperature of the storage after volume withdrawn (for Heating)

        // Run over multiple heat sources
        let mut temp_after_prev_heat_source = temp_s3_n.clone();
        let mut q_ls = 0.0;
        *self.q_ls_n_prev_heat_source.write() = vec![0.0; self.number_of_volumes];
        // In Python extra variables initialized and assigned here
        // for the purpose of passing them to the testoutput method
        // which we have decided to port (for now)
        let mut temp_s8_n = temp_s3_n;

        for (heat_source_name, positioned_heat_source) in self.heat_source_data.clone() {
            let (_, _setpntmax) = positioned_heat_source.heat_source.lock().setpnt(simtime)?;
            let heater_layer =
                (positioned_heat_source.heater_position * self.number_of_volumes as f64) as usize;

            // In cases where there is no thermostat or tank is one layer, set the thermostat layer to the heater layer
            let thermostat_layer = match positioned_heat_source.thermostat_position {
                Some(thermostat_position) => {
                    (thermostat_position * self.number_of_volumes as f64) as usize
                }
                None => heater_layer,
            };

            let TemperatureCalculation {
                temp_s8_n: temp_s8_n_step,
                q_ls: q_ls_this_heat_source,
                q_ls_n: q_ls_n_this_heat_source,
                ..
            } = self.run_heat_sources(
                temp_after_prev_heat_source.clone(),
                &positioned_heat_source.heat_source.lock(),
                &heat_source_name,
                heater_layer,
                thermostat_layer,
                &self.q_ls_n_prev_heat_source.read().clone(),
                simtime,
            )?;

            temp_after_prev_heat_source = temp_s8_n_step.clone();
            q_ls += q_ls_this_heat_source;

            for (i, q_ls_n) in q_ls_n_this_heat_source.iter().enumerate() {
                let mut q_ls_n_prev = self.q_ls_n_prev_heat_source.write();
                q_ls_n_prev[i] += q_ls_n;
            }

            temp_s8_n = temp_s8_n_step;

            // Trigger heating to stop
            self.determine_heat_source_switch_off(
                &temp_s8_n,
                &heat_source_name,
                positioned_heat_source,
                heater_layer,
                thermostat_layer,
                simtime,
            )?;
        }

        // Additional calculations
        // 6.4.6 Calculation of the auxiliary energy
        // accounted for elsewhere so not included here
        let w_sto_aux = 0.;

        // 6.4.7 Recoverable, recovered thermal losses
        // recoverable auxiliary energy transmitted to the heated space - kWh
        let q_sto_h_rbl_aux =
            w_sto_aux * THERMAL_CONSTANTS_F_STO_M * (1. - THERMAL_CONSTANTS_F_RVD_AUX);
        // recoverable heat losses (storage) - kWh
        let q_sto_h_rbl_env = q_ls * THERMAL_CONSTANTS_F_STO_M;
        // total recoverable heat losses for heating - kWh
        self.q_sto_h_ls_rbl
            .store(q_sto_h_rbl_env + q_sto_h_rbl_aux, Ordering::SeqCst);

        // set temperatures calculated to be initial temperatures of volumes for the next timestep
        *self.temp_n.write() = temp_s8_n;

        // TODO (from Python) recoverable heat losses for heating should impact heating

        // Return total energy of hot water supplied and unmet
        Ok(q_use_w)
    }

    /// Allocate hot water layers to meet a single temperature demand.
    ///
    /// Arguments:
    /// * `event` -- Dictionary containing information about the draw-off event
    ///              (e.g. {'start': 18, 'duration': 1, 'temperature': 41.0, 'type': 'Other', 'name': 'other', 'warm_volume': 8.0})
    fn extract_hot_water(
        &self,
        event: TypedScheduleEvent,
        simulation_time: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64, Vec<f64>)> {
        // Make a copy of the volume list to keep track of remaining volumes
        // Remaining volume of water in storage tank layers
        let mut remaining_vols = self.vol_n.clone();

        // Extract the temperature and required hot volume from the event
        let hot_volume = event.volume_hot.ok_or_else(|| anyhow!("Storage tank hot water events must have a volume_hot."));
        let hot_volume = event.volume_hot.unwrap();

        // # Remaining volume of hot water to be satisfied for current event
        let mut remaining_demanded_volume = hot_volume;
        let mut energy_withdrawn = 0.;

        let mut temp_average_drawoff_volweighted: f64 =
            self.temp_average_drawoff_volweighted.load(Ordering::SeqCst);
        let mut total_volume_drawoff: f64 = self.total_volume_drawoff.load(Ordering::SeqCst);
        let mut last_layer_index: usize = Default::default();

        //  Loop through storage layers (starting from the top)
        for (layer_index, &layer_temp) in self.temp_n.read().iter().enumerate().rev() {
            last_layer_index = layer_index;
            let layer_vol = remaining_vols[layer_index];

            if remaining_demanded_volume <= 0. {
                break;
            }

            // Skip this layer if its remaining volume is already zero
            if remaining_vols[layer_index] <= 0. {
                continue;
            }

            let required_vol: f64;
            // Volume of hot water required at this layer
            if layer_vol <= remaining_demanded_volume {
                // This is the case where layer cannot meet all remaining demand for this event
                required_vol = layer_vol;
                // Deduct the required volume from the remaining demand and update the layer's volume
                remaining_vols[layer_index] -= layer_vol;
                remaining_demanded_volume -= layer_vol;
            } else {
                //This is the case where layer can meet all remaining demand for this event
                required_vol = remaining_demanded_volume;
                // Deduct the required volume from the remaining demand and update the layer's volume
                remaining_vols[layer_index] -= required_vol;
                remaining_demanded_volume = 0.0;
            }

            temp_average_drawoff_volweighted += required_vol * layer_temp;
            total_volume_drawoff += required_vol;

            // Record the met volume demand for the current temperature target
            // vol_removed is the volume of warm water that has been satisfied from hot water in this layer

            let list_temp_vol = self
                .cold_feed
                .get_temp_cold_water(hot_volume, simulation_time);
            let sum_t_by_v: f64 = list_temp_vol.iter().map(|(t, v)| t * v).sum();
            let sum_v: f64 = list_temp_vol.iter().map(|(_t, v)| v).sum();
            let temp_cold = sum_t_by_v / sum_v;

            energy_withdrawn +=
                // Calculation with event water parameters
                // self.__rho * self.__Cp * warm_vol_removed * (warm_temp - self.__cold_feed.temperature())
                // Calculation with layer water parameters
                self.rho
                    * self.cp
                    * required_vol
                    * (layer_temp - temp_cold)
        }

        self.temp_average_drawoff_volweighted
            .store(temp_average_drawoff_volweighted, Ordering::SeqCst);
        self.total_volume_drawoff
            .store(total_volume_drawoff, Ordering::SeqCst);

        //  Calculate the remaining total volume
        let remaining_total_volume: f64 = remaining_vols.iter().sum();

        //  Calculate the total volume used
        let volume_used = self.volume_total_in_litres - remaining_total_volume;

        Ok((volume_used, energy_withdrawn, remaining_vols))
    }

    /// Calculate the new temperature distribution after displacement.
    /// Arguments:
    /// * `remaining_vols` -- List of remaining volumes for each storage layer after draw-off
    /// * `temp_cold` -- Temperature of the cold water being added
    fn calc_temps_after_extraction(
        &self,
        mut remaining_vols: Vec<f64>,
        simulation_time: SimulationTimeIteration,
    ) -> (Vec<f64>, bool) {
        let mut new_temps = self.temp_n.read().clone();

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
            let mut volume_weighted_temperature = remaining_vols[i] * self.temp_n.read()[i];

            // Initialisation of min temperature of tank layers to compare eventually against
            // the 'cold' feed temperature to check if rearrangement is needed.
            let mut temp_layer_min = self.temp_n.read()[i];

            // Add water from the layers below to this layer
            for j in (0..i).rev() {
                let available_volume = remaining_vols[j];
                if available_volume > 0. {
                    // Determine the volume to move up from this layer
                    let move_volume = f64::min(needed_volume, available_volume);
                    remaining_vols[j] -= move_volume;

                    // Adjust the temperature by mixing in the moved volume
                    total_volume += move_volume;
                    volume_weighted_temperature += move_volume * self.temp_n.read()[j];

                    // Update min temperature of the tank so far.
                    {
                        let current_temp = self.temp_n.read()[j];
                        if current_temp < temp_layer_min {
                            temp_layer_min = current_temp;
                        }
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

                let list_temp_vol = self
                    .cold_feed
                    .draw_off_water(needed_volume, simulation_time);
                let sum_t_by_v: f64 = list_temp_vol.iter().map(|(t, v)| t * v).sum();
                let sum_v: f64 = list_temp_vol.iter().map(|(_t, v)| v).sum();

                let temp_cold_feed = sum_t_by_v / sum_v;
                volume_weighted_temperature += needed_volume * temp_cold_feed;
                flag_rearrange_layers = temp_cold_feed > temp_layer_min;
            }

            new_temps[i] = volume_weighted_temperature / total_volume;
            remaining_vols[i] = total_volume;
        }
        (new_temps, flag_rearrange_layers)
    }

    /// When the temperature of the volume i is higher than the one of the upper volume,
    /// then the 2 volumes are melded. This iterative process is maintained until the temperature
    /// of the volume i is lower or equal to the temperature of the volume i+1.
    fn rearrange_temperatures(&self, temp_s6_n: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let mut temp_s7_n = temp_s6_n.to_vec();

        loop {
            // Flag for which layers need mixing
            let mut mix_layer_n: Vec<u8> = vec![0; self.number_of_volumes];

            // #for loop :-1 is important here!
            // #loop through layers from bottom to top, without including top layer.
            // #this is because the top layer has no upper layer to compare too
            // loop through layers from bottom to top, without including top layer;
            for i in 0..self.vol_n.len() - 1 {
                if temp_s7_n[i] >= temp_s7_n[i + 1] {
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
                    mix_layer_n = vec![0; self.number_of_volumes];
                }
            }

            if is_sorted(&temp_s7_n) {
                break;
            }
        }

        let q_h_sto_end = (0..self.vol_n.len())
            .map(|i| self.rho * self.cp * self.vol_n[i] * temp_s7_n[i])
            .collect::<Vec<f64>>();

        (q_h_sto_end, temp_s7_n.to_owned())
    }

    fn run_heat_sources(
        &self,
        temp_s3_n: Vec<f64>,
        heat_source: &HeatSource,
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
            heat_source,
            heat_source_name,
            heater_layer,
            thermostat_layer,
            simulation_time,
        )?;

        self.calc_final_temps(
            &temp_s3_n,
            heat_source,
            q_x_in_n,
            heater_layer,
            q_ls_prev_heat_source,
            simulation_time,
            None,
        )
    }

    /// Energy input for the storage from the generation system
    /// (expressed per energy carrier X)
    /// Heat Source = energy carrier
    fn potential_energy_input(
        // Heat source. Addition of temp_s3_n as an argument
        &self,
        temp_s3_n: &[f64],
        heat_source: &HeatSource,
        heat_source_name: &str,
        heater_layer: usize,
        thermostat_layer: usize,
        simulation_time: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<f64>> {
        // initialise list of potential energy input for each layer
        let mut q_x_in_n = vec![0.; self.number_of_volumes];

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
                    heat_source_name,
                    heat_source,
                    heater_layer,
                    thermostat_layer,
                    simulation_time,
                )?;

                let default_temp_flow = self.temp_n.read()[heater_layer];
                let temp_flow = self
                    .temp_flow(heat_source, simulation_time)
                    .unwrap_or(default_temp_flow);
                if self.heating_active[heat_source_name].load(Ordering::SeqCst) {
                    // upstream Python uses duck-typing/ polymorphism here, but we need to be more explicit
                    let mut energy_potential = match heat_source {
                        HeatSource::Storage(HeatSourceWithStorageTank::Immersion(
                            immersion_heater,
                        )) => immersion_heater
                            .lock()
                            .energy_output_max(simulation_time, false),
                        HeatSource::Storage(HeatSourceWithStorageTank::Solar(_)) => unreachable!(), // this case was already covered in the first arm of this if let clause, so can't repeat here
                        HeatSource::Wet(heat_source_wet) => {
                            // TODO Use different temperatures for flow and return in the call to
                            // heat_source.energy_output_max below
                            // Fallback to current tank temperature at heater layer when heat source has no setpoint
                            heat_source_wet.energy_output_max(
                                Some(temp_flow),
                                temp_flow,
                                simulation_time,
                            )?
                        }
                    };

                    // TODO (from Python) Consolidate checks for systems with/without primary pipework
                    if !matches!(
                        heat_source,
                        HeatSource::Storage(HeatSourceWithStorageTank::Immersion(_))
                    ) {
                        let (primary_pipework_losses_kwh, _) = self
                            .calculate_primary_pipework_losses(
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

    fn calc_final_temps(
        &self,
        temp_s3_n: &[f64],
        heat_source: &HeatSource,
        q_x_in_n: Vec<f64>,
        heater_layer: usize,
        q_ls_n_prev_heat_source: &[f64],
        simtime: SimulationTimeIteration,
        control_max_diverter: Option<&Control>,
    ) -> anyhow::Result<TemperatureCalculation> {
        let setpntmax = if let Some(control_max_diverter) = control_max_diverter {
            control_max_diverter.setpnt(&simtime)
        } else {
            let (_, setpntmax) = self.retrieve_setpnt(heat_source, simtime)?;
            setpntmax
        };

        let (q_s6, temp_s6_n) = self.calc_temps_with_energy_input(temp_s3_n, &q_x_in_n);

        // 6.4.3.9 STEP 7 Re-arrange the temperatures in the storage after energy input
        let (q_h_sto_s7, temp_s7_n) = self.rearrange_temperatures(&temp_s6_n);

        // STEP 8 Thermal losses and final temperature
        let (q_in_h_w, q_ls, temp_s8_n, q_ls_n) = self.calc_temps_after_thermal_losses(
            &temp_s7_n,
            &q_x_in_n,
            q_h_sto_s7,
            heater_layer,
            q_ls_n_prev_heat_source,
            setpntmax,
        );

        // TODO (from Python) 6.4.3.11 Heat exchanger

        // demand adjusted energy from heat source (before was just using potential without taking it)
        let input_energy_adj = q_in_h_w;

        #[cfg(test)]
        {
            self.energy_demand_test
                .store(input_energy_adj, Ordering::SeqCst);
        }

        let _heat_source_output =
            self.heat_source_output(heat_source, input_energy_adj, heater_layer, simtime, None);
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
    fn calc_temps_with_energy_input(&self, temp_s3_n: &[f64], q_x_in_n: &[f64]) -> (f64, Vec<f64>) {
        // initialise list of theoretical variation of temperature of layers in degrees
        let mut delta_temp_n = vec![0.; self.number_of_volumes];
        // initialise list of theoretical temperature of layers after input in degrees
        let mut temp_s6_n = vec![0.; self.number_of_volumes];
        // output energy delivered by the storage in kWh - timestep dependent
        let q_sto_h_out_n: Vec<f64> = vec![0.; self.number_of_volumes];

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
    fn calc_temps_after_thermal_losses(
        &self,
        temp_s7_n: &[f64],
        q_x_in_n: &[f64],
        q_h_sto_s7: Vec<f64>,
        heater_layer: usize,
        q_ls_n_prev_heat_source: &[f64],
        temp_setpntmax: Option<f64>,
    ) -> (f64, f64, Vec<f64>, Vec<f64>) {
        let q_x_in_adj: f64 = q_x_in_n.iter().sum();

        // standby losses coefficient - W/K
        let h_sto_ls = self.stand_by_losses_coefficient();

        // standby losses correction factor - dimensionless
        // note from Python code: "do not think these are applicable so used: f_sto_dis_ls = 1, f_sto_bac_acc = 1"

        // initialise list of thermal losses in kWh
        let mut q_ls_n: Vec<f64> = Vec::with_capacity(self.number_of_volumes);
        // initialise list of final temperature of layers after thermal losses in degrees
        let mut temp_s8_n: Vec<f64> = Vec::with_capacity(self.number_of_volumes);

        // Thermal losses
        // Note: Eqn 13 from BS EN 15316-5:2017 does not explicitly multiply by
        // timestep (it seems to assume a 1 hour timestep implicitly), but it is
        // necessary to convert the rate of heat loss to a total heat loss over
        // the time period
        for i in 0..self.vol_n.len() {
            let temp_before_losses = if let Some(temp_setpntmax) = temp_setpntmax {
                min_of_2(temp_s7_n[i], temp_setpntmax)
            } else {
                temp_s7_n[i]
            };

            let q_ls_n_step = (h_sto_ls * self.rho * self.cp)
                * (self.vol_n[i] / self.volume_total_in_litres)
                * (temp_before_losses - self.ambient_temperature)
                * self.simulation_timestep;

            let q_ls_n_step = max_of_2(0., q_ls_n_step - q_ls_n_prev_heat_source[i]);

            q_ls_n.push(q_ls_n_step);
        }

        // total thermal losses kWh
        let q_ls = q_ls_n.iter().sum();

        self.storage_losses_kwh.store(q_ls, Ordering::SeqCst);

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
            let temp_s8_n_step = if let (true, Some(temp_setpntmax), true) = (
                q_x_in_adj > 0.0,
                temp_setpntmax,
                temp_setpntmax.is_some_and(|t| temp_s7_n[i] > t),
            ) {
                // Case 2 - Temperature exceeding the set point
                temp_setpntmax
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
            if let (Some(temp_setpntmax), true) = (
                temp_setpntmax,
                temp_setpntmax.is_some_and(|t| temp_s7_n[heater_layer] > t),
            ) {
                for i in heater_layer..self.number_of_volumes {
                    energy_surplus += q_h_sto_s7[i]
                        - q_ls_n[i]
                        - (self.rho * self.cp * self.vol_n[i] * temp_setpntmax);
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

    /// Calculates pipework loss before sending on the demand energy
    fn heat_source_output(
        &self,
        heat_source: &HeatSource,
        input_energy_adj: f64,
        _heater_layer: usize,
        simulation_time_iteration: SimulationTimeIteration,
        smart_hot_water_tank: Option<&SmartHotWaterTank>, // the temp_flow method might need to be called as a smart hot water tank if this is a storage tank composed by a smart hot water tank
    ) -> anyhow::Result<f64> {
        // if immersion heater, no pipework losses
        // TODO (from Python):  Critical - temp_flow cannot be None for downstream method calculate_primary_pipework_losses
        // but providing a fallback value will change the e2e test results
        let temp_flow = match smart_hot_water_tank {
            None => self.temp_flow(heat_source, simulation_time_iteration),
            Some(smart_hot_water_tank) => smart_hot_water_tank.temp_flow(simulation_time_iteration),
        }?;

        // Input energy rounded so that almost zero negative numbers (caused by
        // floating point error) do not cause errors in subsequent code

        let input_energy_adj = round_by_precision(input_energy_adj, 1e10);

        match heat_source {
            HeatSource::Storage(HeatSourceWithStorageTank::Immersion(immersion)) => immersion
                .lock()
                .demand_energy(input_energy_adj, simulation_time_iteration),
            HeatSource::Storage(HeatSourceWithStorageTank::Solar(solar)) => Ok(solar
                .lock()
                .demand_energy(input_energy_adj, simulation_time_iteration.index)),
            HeatSource::Wet(ref wet_heat_source) => {
                let (primary_pipework_losses_kwh, primary_gains) = self
                    .calculate_primary_pipework_losses(
                        input_energy_adj,
                        temp_flow,
                        simulation_time_iteration,
                    );
                let input_energy_adj = input_energy_adj + primary_pipework_losses_kwh;

                // TODO Use different temperatures for flow and return in the call to
                // heat_source.demand_energy below
                let heat_source_output = wet_heat_source.demand_energy(
                    input_energy_adj,
                    Some(temp_flow),
                    temp_flow,
                    simulation_time_iteration,
                )? - primary_pipework_losses_kwh;
                self.input_energy_adj_prev_timestep
                    .store(input_energy_adj, Ordering::SeqCst);
                self.pipework_primary_gains_for_timestep
                    .store(primary_gains, Ordering::SeqCst);

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
        let (setpntmin, setpntmax) = heat_source.setpnt(simulation_time_iteration)?;

        match (setpntmax, setpntmin) {
            (None, Some(_)) => bail!("setpntmin must be None if setpntmax is None"),
            (Some(setpointmax), Some(setpointmin)) => {
                if setpointmin > setpointmax {
                    bail!("setpntmin: {setpointmin} must not be greater than setpntmax: {setpointmax}");
                }
            }
            _ => {}
        }

        Ok((setpntmin, setpntmax))
    }

    fn determine_heat_source_switch_on(
        &self,
        temp_s3_n: &[f64],
        heat_source_name: &str,
        heat_source: &HeatSource,
        _heater_layer: usize,
        thermostat_layer: usize,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<()> {
        let (setpntmin, _) = self.retrieve_setpnt(heat_source, simulation_time_iteration)?;
        if setpntmin.is_some() && temp_s3_n[thermostat_layer] <= setpntmin.unwrap() {
            self.heating_active[heat_source_name].store(true, Ordering::SeqCst);
        };
        Ok(())
    }

    fn determine_heat_source_switch_off(
        &self,
        temp_s8_n: &[f64],
        heat_source_name: &str,
        heat_source: PositionedHeatSource,
        _heater_layer: usize,
        thermostat_layer: usize,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<()> {
        let heat_source = heat_source.heat_source;
        let (_, setpntmax) =
            self.retrieve_setpnt(&(heat_source.lock()), simulation_time_iteration)?;

        if setpntmax.is_none() || temp_s8_n[thermostat_layer] >= setpntmax.unwrap() {
            self.heating_active[heat_source_name].store(false, Ordering::SeqCst);
        };

        Ok(())
    }

    fn temp_flow(
        &self,
        heat_source: &HeatSource,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let (_, setpntmax) = self.retrieve_setpnt(heat_source, simulation_time_iteration)?;

        let setpntmax = if let Some(setpntmax) = setpntmax {
            setpntmax
        } else {
            self.temp_flow_prev.load(Ordering::SeqCst)
        };

        self.temp_flow_prev.store(setpntmax, Ordering::SeqCst);

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
    pub(crate) fn get_temp_hot_water(
        &self,
        volume_req: f64,
        volume_req_already: Option<f64>,
    ) -> Vec<(f64, f64)> {
        let mut volume_req = volume_req;
        let volume_req_already = volume_req_already.unwrap_or(0.);
        let mut volume_req_cumulative = volume_req + volume_req_already;

        let mut list_temp_vol: Vec<(f64, f64)> = vec![];
        // Loop through storage layers (starting from the top)
        // TODO (from Python) Handle case where we reach bottom of tank
        for (layer_index, &layer_temp) in self.temp_n.read().iter().enumerate().rev() {
            let layer_vol = self.vol_n[layer_index];
            let volume_from_current_layer = volume_req_cumulative.min(layer_vol);

            list_temp_vol.push((layer_temp, volume_from_current_layer));
            volume_req_cumulative -= volume_from_current_layer;

            if volume_req_cumulative <= 0. {
                break;
            }
        }

        // Base temperature on the part of the draw-off for volume_req, and
        // ignore any volume previously considered
        let mut list_temp_vol_req: Vec<(f64, f64)> = vec![];
        for (layer_temp, layer_vol) in list_temp_vol.iter().rev() {
            let volume_from_current_layer = volume_req.min(*layer_vol);
            list_temp_vol_req.push((*layer_temp, volume_from_current_layer));
            volume_req -= volume_from_current_layer;

            if volume_req < 0. {
                break;
            }
        }

        list_temp_vol_req.into_iter().rev().collect_vec()
    }

    /// Appendix B B.2.8 Stand-by losses are usually determined in terms of energy losses during
    /// a 24h period. Formula (B.2) allows the calculation of _sto_stbl_ls_tot based on a reference
    /// value of the daily thermal energy losses.
    ///
    /// h_sto_ls is the stand-by losses, in W/K
    ///
    /// TODO (from Python) there are alternative methods listed in App B (B.2.8) which are not included here.
    fn stand_by_losses_coefficient(&self) -> f64 {
        // BS EN 12897:2016 appendix B B.2.2
        // temperature of the water in the storage for the standardized conditions - degrees
        // these are reference (ref) temperatures from the standard test conditions for cylinder loss.
        let temp_set_ref = 65.;
        let temp_amb_ref = 20.;

        (1000. * self.q_std_ls_ref) / (24. * (temp_set_ref - temp_amb_ref))
    }

    /// Function added into Storage tank to be called by the Solar Thermal object.
    /// Calculates the impact on storage tank temperature due to the proposed energy input
    fn storage_tank_potential_effect(&self, energy_proposed: f64, temp_s3_n: &[f64]) -> (f64, f64) {
        // assuming initially no water draw-off

        // initialise list of potential energy input for each layer
        let mut q_x_in_n = vec![0.; self.number_of_volumes];

        // TODO (from Python) - ensure we are feeding in the correct volume
        q_x_in_n[0] = energy_proposed;

        let (_q_s6, temp_s6_n) = self.calc_temps_with_energy_input(temp_s3_n, &q_x_in_n);

        // 6.4.3.9 STEP 7 Re-arrange the temperatures in the storage after energy input
        let (_q_h_sto_s7, temp_s7_n) = self.rearrange_temperatures(&temp_s6_n);

        // TODO (from Python) Check [0] is bottom layer temp and that solar thermal inlet is top layer NB_VOL-1
        (temp_s7_n[0], temp_s7_n[self.number_of_volumes - 1])
    }

    /// Send more intermediate output parameters to report
    pub(crate) fn get_losses_from_primary_pipework_and_storage(&self) -> (f64, f64) {
        (
            self.primary_pipework_losses_kwh.load(Ordering::SeqCst),
            self.storage_losses_kwh.load(Ordering::SeqCst),
        )
    }

    // NB. there is a testoutput() function here in the Python to output to a test file - will not reimplement unless seen as necessary

    /// draw off hot water layers until required volume is provided.
    ///
    /// Arguments:
    /// * volume    -- volume of water required
    pub(crate) fn draw_off_hot_water(
        &self,
        volume: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (Option<f64>, f64) {
        if volume.abs() <= 1e-10 {
            return (None, volume);
        }

        // Remaining volume of water in storage tank layers
        let mut remaining_vols = self.vol_n.clone();

        let mut remaining_demanded_volume = volume;

        // Initialize the unmet and met energies
        let mut _energy_withdrawn = 0.0;
        self.temp_average_drawoff_volweighted
            .store(0.0, Ordering::SeqCst);
        self.total_volume_drawoff.store(0.0, Ordering::SeqCst);

        let list_temp_vol = self
            .cold_feed
            .get_temp_cold_water(volume, simulation_time_iteration);
        let sum_t_by_v: f64 = list_temp_vol.iter().map(|(t, v)| t * v).sum();
        let sum_v: f64 = list_temp_vol.iter().map(|(_t, v)| v).sum();

        self.temp_average_drawoff
            .store(sum_t_by_v / sum_v, Ordering::SeqCst);

        let _temp_ini_n = self.temp_n.clone();
        let temp_s3_n = self.temp_n.clone();

        // Loop through storage layers (starting from the top)
        for (layer_index, &layer_temp) in self.temp_n.read().iter().enumerate().rev() {
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

            self.temp_average_drawoff_volweighted
                .fetch_add(required_vol * layer_temp, Ordering::SeqCst);
            self.total_volume_drawoff
                .fetch_add(required_vol, Ordering::SeqCst);

            let list_temp_vol = self
                .cold_feed
                .get_temp_cold_water(required_vol, simulation_time_iteration);
            let sum_t_by_v: f64 = list_temp_vol.iter().map(|(t, v)| t * v).sum();
            let sum_v: f64 = list_temp_vol.iter().map(|(_t, v)| v).sum();
            let temp_cold_water = sum_t_by_v / sum_v;
            //  Record the met volume demand for the current temperature target
            //  warm_vol_removed is the volume of warm water that has been satisfied from hot water in this layer
            _energy_withdrawn +=
                //  Calculation with event water parameters
                // self.__rho * self.__Cp * warm_vol_removed * (warm_temp - self.__cold_feed.temperature())
                //  Calculation with layer water parameters
                self.rho
                    * self.cp
                    * required_vol
                    * (layer_temp - temp_cold_water);

            if remaining_demanded_volume <= 0.0 {
                break;
            }
        }

        if remaining_demanded_volume > 0.0 {
            let list_temp_vol = self
                .cold_feed
                .get_temp_cold_water(remaining_demanded_volume, simulation_time_iteration);
            let sum_t_by_v: f64 = list_temp_vol.iter().map(|(t, v)| t * v).sum();
            let sum_v: f64 = list_temp_vol.iter().map(|(_t, v)| v).sum();
            let temp_cold_water = sum_t_by_v / sum_v;

            self.temp_average_drawoff_volweighted.fetch_add(
                remaining_demanded_volume * temp_cold_water,
                Ordering::SeqCst,
            );
            self.total_volume_drawoff
                .fetch_add(remaining_demanded_volume, Ordering::SeqCst);
        }

        self.temp_average_drawoff.store(
            self.temp_average_drawoff_volweighted.load(Ordering::SeqCst)
                / self.total_volume_drawoff.load(Ordering::SeqCst),
            Ordering::SeqCst,
        );
        // Determine the new temperature distribution after displacement
        let (mut new_temp_distribution, flag_rearrange_layers) =
            self.calc_temps_after_extraction(remaining_vols, simulation_time_iteration);

        if flag_rearrange_layers {
            // Re-arrange the temperatures in the storage after energy input from pre-heated tank
            (_, new_temp_distribution) = self.rearrange_temperatures(&new_temp_distribution);
        }

        *self.temp_n.write() = new_temp_distribution;

        // Return the average temperature and volume drawn
        (
            Some(self.temp_average_drawoff.load(Ordering::SeqCst)),
            self.total_volume_drawoff.load(Ordering::SeqCst),
        )
    }
    fn additional_energy_input(
        &self,
        heat_source: &HeatSource,
        heat_source_name: &str,
        energy_input: f64,
        control_max_diverter: Option<&Control>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        if energy_input == 0. {
            return Ok(0.);
        }

        let heat_source_data = &self.heat_source_data[heat_source_name];

        let heater_layer =
            (heat_source_data.heater_position * self.number_of_volumes as f64) as usize;

        let mut q_x_in_n = vec![0.; self.number_of_volumes];
        q_x_in_n[heater_layer] = energy_input;
        let TemperatureCalculation {
            temp_s8_n,
            q_in_h_w,
            q_ls_n: q_ls_n_this_heat_source,
            ..
        } = self.calc_final_temps(
            &self.temp_n.read(),
            heat_source,
            q_x_in_n,
            heater_layer,
            &self.q_ls_n_prev_heat_source.read(),
            simulation_time_iteration,
            control_max_diverter,
        )?;

        for (i, q_ls_n) in q_ls_n_this_heat_source.iter().enumerate() {
            let mut q_ls_n_prev = self.q_ls_n_prev_heat_source.write();
            q_ls_n_prev[i] += *q_ls_n;
        }

        *self.temp_n.write() = temp_s8_n;

        Ok(q_in_h_w)
    }

    #[cfg(test)]
    fn test_energy_demand(&self) -> f64 {
        self.energy_demand_test.load(Ordering::SeqCst)
    }

    /// Return the DHW recoverable heat losses as internal gain for the current timestep in W
    pub(crate) fn internal_gains(&self) -> f64 {
        let primary_gains_timestep = self
            .pipework_primary_gains_for_timestep
            .load(Ordering::SeqCst);
        self.pipework_primary_gains_for_timestep
            .store(0., Ordering::SeqCst);

        (self.q_sto_h_ls_rbl.load(Ordering::SeqCst) * WATTS_PER_KILOWATT as f64
            / self.simulation_timestep)
            + primary_gains_timestep
    }

    fn calculate_primary_pipework_losses(
        &self,
        input_energy_adj: f64,
        temp_flow: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (f64, f64) {
        let mut primary_pipework_losses_kwh: f64 = Default::default();
        let mut primary_gains_w = Default::default();
        if let Some(primary_pipework) = self.primary_pipework.as_ref() {
            // start of heating event
            if input_energy_adj > 0.
                && self.input_energy_adj_prev_timestep.load(Ordering::SeqCst) == 0.
            {
                for pipework_data in primary_pipework {
                    let outside_temperature = self.temp_surrounding_primary_pipework(
                        &pipework_data,
                        simulation_time_iteration,
                    );
                    let cool_down_loss =
                        pipework_data.calculate_cool_down_loss(temp_flow, outside_temperature);

                    primary_pipework_losses_kwh += cool_down_loss;
                }
            }

            // during heating event
            if input_energy_adj > 0. {
                for pipework_data in primary_pipework {
                    // Primary losses for the timestep calculated from temperature difference

                    let outside_temperature = self.temp_surrounding_primary_pipework(
                        &pipework_data,
                        simulation_time_iteration,
                    );
                    let primary_pipework_losses_w = pipework_data
                        .calculate_steady_state_heat_loss(temp_flow, outside_temperature);

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
            if input_energy_adj == 0.
                && self.input_energy_adj_prev_timestep.load(Ordering::SeqCst) > 0.
            {
                for pipework_data in primary_pipework {
                    let location = pipework_data.location();
                    match location {
                        PipeworkLocation::External => {}
                        PipeworkLocation::Internal => {
                            primary_gains_w += pipework_data.calculate_cool_down_loss(
                                temp_flow,
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
        self.primary_pipework_losses_kwh
            .store(primary_pipework_losses_kwh, Ordering::SeqCst);

        (primary_pipework_losses_kwh, primary_gains_w)
    }

    // TODO Python has get_temp_cold_water and draw_off_water defined here
    // which are called but currently unsure where from

    /// Return the pre-heated water temperature for the current timestep and the volume drawn
    pub(crate) fn get_temp_cold_water(&self, volume_needed: f64) -> Vec<(f64, f64)> {
        // TODO this matches Python - is it correct?
        self.get_temp_hot_water(volume_needed, None)
    }

    /// Return the pre-heated water temperature for the current timestep and the volume drawn
    pub(crate) fn draw_off_water(
        &self,
        volume_needed: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> Vec<(f64, f64)> {
        let list_temp_vol = self.get_temp_cold_water(volume_needed);
        self.draw_off_hot_water(volume_needed, simulation_time_iteration);
        list_temp_vol
    }

    pub(crate) fn output_results(&self) -> Option<Vec<StorageTankDetailedResult>> {
        None // TODO implement this
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

/// A struct to represent a smart hot water storage tank/cylinder
#[derive(Debug)]
pub(crate) struct SmartHotWaterTank {
    storage_tank: StorageTank,
    power_pump_kw: f64,
    max_flow_rate_pump_l_per_min: f64,
    temp_usable: f64,
    temp_setpnt_max: Arc<Control>,
    energy_supply_connection_pump: EnergySupplyConnection,
}

impl SmartHotWaterTank {
    /// Construct a SmartHotWaterTank object
    ///
    /// Arguments:
    /// * `volume` - total volume of the tank, in litres
    /// * `losses` - measured standby losses due to cylinder insulation
    ///                                at standardised conditions, in kWh/24h
    /// * `init_temp` - initial temperature required for DHW
    /// * `power_pump_kw` - power of pump used to pump water from the bottom
    ///                                 to the top of the tank in kW
    /// * `max_flow_rate_pump_l_per_min` - maximum flow rate that pump can provide in l/min
    /// * `temp_usable` - lowest water temperature that the water can be useable
    /// * `temp_setpnt_max` - maximum set point temperature
    /// * `cold_feed` - reference to ColdWaterSource object
    /// * `heat_sources` - dict where keys are heat source objects and
    ///                                values are tuples of heater and thermostat
    ///                                position
    /// * `number_of_volumes` -
    ///                               number of volumes the storage is modelled with
    ///                               see App.C (C.1.2 selection of the number of volumes to model the storage unit)
    ///                               for more details if this wants to be changed.
    /// * `energy_supply_conn_pump`
    /// * `contents` - reference to MaterialProperties object
    pub(crate) fn new(
        volume: f64,
        losses: f64,
        init_temp: f64,
        power_pump_kw: f64,
        max_flow_rate_pump_l_per_min: f64,
        temp_usable: f64,
        temp_setpnt_max: Arc<Control>,
        cold_feed: WaterSourceWithTemperature,
        simulation_timestep: f64,
        heat_sources: IndexMap<String, PositionedHeatSource>,
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions: Arc<ExternalConditions>,
        detailed_output: Option<bool>,
        number_of_volumes: Option<usize>,
        primary_pipework_lst: Option<&Vec<WaterPipework>>,
        energy_supply_conn_pump: EnergySupplyConnection,
        contents: Option<MaterialProperties>,
    ) -> anyhow::Result<Self> {
        let detailed_output = detailed_output.unwrap_or(false);
        let number_of_volumes = number_of_volumes.unwrap_or(100);
        let contents = contents.unwrap_or(*WATER);

        let storage_tank = StorageTank::new(
            volume,
            losses,
            init_temp,
            cold_feed,
            simulation_timestep,
            heat_sources,
            temp_internal_air_fn,
            external_conditions,
            number_of_volumes.into(),
            primary_pipework_lst,
            contents,
            None,
            None,
            None,
            detailed_output,
        )?;

        Ok(Self {
            storage_tank,
            power_pump_kw,
            max_flow_rate_pump_l_per_min,
            temp_usable,
            temp_setpnt_max,
            energy_supply_connection_pump: energy_supply_conn_pump,
        })
    }

    fn retrieve_setpnt(
        &self,
        heat_source: &HeatSource,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<(Option<f64>, Option<f64>)> {
        // N.B. implementation from StorageTank:
        let (setpntmin, setpntmax) = self
            .storage_tank
            .retrieve_setpnt(heat_source, simulation_time_iteration)?;

        // N.B. extra checks specific to SmartHotWaterTank
        if let Some(setpntmin) = setpntmin {
            if !(0. ..=1.).contains(&setpntmin) {
                bail!(">= 0. and <= 1. required for setpoints");
            }
        }

        if let Some(setpntmax) = setpntmax {
            if !(0. ..=1.).contains(&setpntmax) {
                bail!(">= 0. and <= 1. required for setpoints");
            }
        }

        Ok((setpntmin, setpntmax))
    }

    // Inherited methods from StorageTank (NB. in Python, SmartHotWaterTank inherits from StorageTank)

    /// Draw off hot water from the tank
    /// Energy calculation as per BS EN 15316-5:2017 Method A sections 6.4.3, 6.4.6, 6.4.7
    /// Modification of calculation based on volumes and actual temperatures for each layer of water in the tank
    /// instead of the energy stored in the layer and a generic temperature (self.temp_out_w_min) = min_temp
    /// to decide if the tank can satisfy the demand (this was producing unnecesary unmet demand for strict high
    /// temp_out_w_min values
    /// Arguments:
    /// * `usage_events` -- All draw off events for the timestep
    pub(crate) fn demand_hot_water(
        &self,
        usage_events: Option<Vec<TypedScheduleEvent>>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // N.B. implementation from StorageTank but calling SmartHotWaterTank specific methods further down
        let mut q_use_w = 0.;
        let mut _volume_demanded = 0.;

        let mut temp_s3_n = self.storage_tank.temp_n.read().clone();

        self.storage_tank
            .temp_average_drawoff_volweighted
            .store(0., Ordering::SeqCst);
        self.storage_tank
            .temp_final_drawoff
            .store(0., Ordering::SeqCst);
        self.storage_tank
            .total_volume_drawoff
            .store(0., Ordering::SeqCst);
        self.storage_tank
            .temp_average_drawoff
            .store(self.storage_tank.initial_temperature, Ordering::SeqCst);

        for event in usage_events.iter().flatten() {
            let (volume_used, energy_withdrawn, remaining_vols) =
                self.storage_tank.extract_hot_water(event.clone(), simtime)?;

            let (temp_s3_n_new, rearrange) = self
                .storage_tank
                .calc_temps_after_extraction(remaining_vols, simtime);
            temp_s3_n = temp_s3_n_new;

            if rearrange {
                // Re-arrange the temperatures in the storage after energy input from pre-heated tank
                temp_s3_n = self.storage_tank.rearrange_temperatures(&temp_s3_n).1
            }

            *self.storage_tank.temp_n.write() = temp_s3_n.clone();

            _volume_demanded += volume_used;
            q_use_w += energy_withdrawn;
        }

        self.storage_tank.temp_average_drawoff.store(
            match self
                .storage_tank
                .total_volume_drawoff
                .load(Ordering::SeqCst)
            {
                value if value != 0. => {
                    let temp_average_drawoff_volweighted = self
                        .storage_tank
                        .temp_average_drawoff_volweighted
                        .load(Ordering::SeqCst);
                    temp_average_drawoff_volweighted / value
                }
                _ => temp_s3_n
                    .last()
                    .copied()
                    .ok_or_else(|| anyhow!("temp_s3_n was unexpectedly empty"))?,
            },
            Ordering::SeqCst,
        );

        // Run over multiple heat sources
        let mut temp_after_prev_heat_source = temp_s3_n.clone();
        let mut q_ls = 0.0;
        *self.storage_tank.q_ls_n_prev_heat_source.write() =
            vec![0.0; self.storage_tank.number_of_volumes];

        let mut temp_s8_n = vec![0.; self.storage_tank.number_of_volumes];

        for (heat_source_name, positioned_heat_source) in self.storage_tank.heat_source_data.clone()
        {
            let (_, _setpntmax) = positioned_heat_source.heat_source.lock().setpnt(simtime)?;
            let heater_layer = (positioned_heat_source.heater_position
                * self.storage_tank.number_of_volumes as f64)
                as usize;

            // In cases where there is no thermostat or tank is one layer, set the thermostat layer to the heater layer
            let thermostat_layer = match positioned_heat_source.thermostat_position {
                Some(thermostat_position) => {
                    (thermostat_position * self.storage_tank.number_of_volumes as f64) as usize
                }
                None => heater_layer,
            };

            // N.B run_heat_sources is the SmartStorageTank specific call
            let TemperatureCalculation {
                temp_s8_n: temp_s8_n_step,
                q_ls: q_ls_this_heat_source,
                q_ls_n: q_ls_n_this_heat_source,
                ..
            } = self.run_heat_sources(
                temp_after_prev_heat_source.clone(),
                &positioned_heat_source.heat_source.lock(),
                &heat_source_name,
                heater_layer,
                thermostat_layer,
                &self.storage_tank.q_ls_n_prev_heat_source.read().clone(),
                simtime,
            )?;

            temp_after_prev_heat_source = temp_s8_n_step.clone();
            q_ls += q_ls_this_heat_source;

            for (i, q_ls_n) in q_ls_n_this_heat_source.iter().enumerate() {
                let mut q_ls_n_prev = self.storage_tank.q_ls_n_prev_heat_source.write();
                q_ls_n_prev[i] += q_ls_n;
            }

            temp_s8_n = temp_s8_n_step;

            // Trigger heating to stop
            self.storage_tank.determine_heat_source_switch_off(
                &temp_s8_n,
                &heat_source_name,
                positioned_heat_source,
                heater_layer,
                thermostat_layer,
                simtime,
            )?;
        }

        // Additional calculations
        // 6.4.6 Calculation of the auxiliary energy
        // accounted for elsewhere so not included here
        let w_sto_aux = 0.;

        // 6.4.7 Recoverable, recovered thermal losses
        // recoverable auxiliary energy transmitted to the heated space - kWh
        let q_sto_h_rbl_aux =
            w_sto_aux * THERMAL_CONSTANTS_F_STO_M * (1. - THERMAL_CONSTANTS_F_RVD_AUX);
        // recoverable heat losses (storage) - kWh
        let q_sto_h_rbl_env = q_ls * THERMAL_CONSTANTS_F_STO_M;
        // total recoverable heat losses for heating - kWh
        self.storage_tank
            .q_sto_h_ls_rbl
            .store(q_sto_h_rbl_env + q_sto_h_rbl_aux, Ordering::SeqCst);

        // set temperatures calculated to be initial temperatures of volumes for the next timestep
        *self.storage_tank.temp_n.write() = temp_s8_n;

        // TODO (from Python) recoverable heat losses for heating should impact heating

        // Return total energy of hot water supplied and unmet
        Ok(q_use_w)
    }

    fn run_heat_sources(
        &self,
        temp_s3_n: Vec<f64>,
        heat_source: &HeatSource,
        heat_source_name: &str,
        heater_layer: usize,
        thermostat_layer: usize,
        q_ls_prev_heat_source: &[f64],
        simulation_time: SimulationTimeIteration,
    ) -> anyhow::Result<TemperatureCalculation> {
        // N.B.: implementation from StorageTank but without thermostat_layer:

        // 6.4.3.8 STEP 6 Energy input into the storage
        // input energy delivered to the storage in kWh - timestep dependent

        // N.B. we're calling the SmartHotWaterTank specific method here
        let q_x_in_n = self.potential_energy_input(
            &temp_s3_n,
            heat_source,
            heat_source_name,
            heater_layer,
            thermostat_layer,
            simulation_time,
        )?;

        // N.B. we're calling the SmartHotWaterTank specific method here
        self.calc_final_temps(
            temp_s3_n,
            heat_source,
            q_x_in_n,
            heater_layer,
            q_ls_prev_heat_source,
            None,
            simulation_time,
        )
    }

    /// Energy input for the storage from the generation system
    /// (expressed per energy carrier X)
    /// Heat Source = energy carrier
    fn potential_energy_input(
        // Heat source. Addition of temp_s3_n as an argument
        &self,
        temp_s3_n: &[f64],
        heat_source: &HeatSource,
        heat_source_name: &str,
        heater_layer: usize,
        thermostat_layer: usize,
        simulation_time: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<f64>> {
        // N.B. implementation from StorageTank but without thermostat_layer & with calling a SmartHotWaterTank specific method
        // initialise list of potential energy input for each layer
        // initialise list of potential energy input for each layer
        let mut q_x_in_n = vec![0.; self.storage_tank.number_of_volumes];

        let energy_potential =
            if let HeatSource::Storage(HeatSourceWithStorageTank::Solar(ref solar_heat_source)) =
                heat_source
            {
                // we are passing the storage tank object to the SolarThermal as this needs to call back the storage tank (sic from Python)
                solar_heat_source.lock().energy_output_max(
                    &self.storage_tank,
                    temp_s3_n,
                    &simulation_time,
                )
            } else {
                // N.B calling the SmartStorageTank specific method here
                self.determine_heat_source_switch_on(
                    temp_s3_n,
                    heat_source_name,
                    heat_source,
                    heater_layer,
                    thermostat_layer,
                    simulation_time,
                )?;

                let default_temp_flow = self.storage_tank.temp_n.read()[heater_layer];
                let temp_flow = self.temp_flow(simulation_time).unwrap_or(default_temp_flow);
                if self.storage_tank.heating_active[heat_source_name].load(Ordering::SeqCst) {
                    // upstream Python uses duck-typing/ polymorphism here, but we need to be more explicit
                    let mut energy_potential = match heat_source {
                        HeatSource::Storage(HeatSourceWithStorageTank::Immersion(
                            immersion_heater,
                        )) => immersion_heater
                            .lock()
                            .energy_output_max(simulation_time, false),
                        HeatSource::Storage(HeatSourceWithStorageTank::Solar(_)) => unreachable!(), // this case was already covered in the first arm of this if let clause, so can't repeat here
                        HeatSource::Wet(heat_source_wet) => {
                            // TODO Use different temperatures for flow and return in the call to
                            // heat_source.energy_output_max below
                            // Fallback to current tank temperature at heater layer when heat source has no setpoint
                            heat_source_wet.energy_output_max(
                                Some(temp_flow),
                                temp_flow,
                                simulation_time,
                            )?
                        }
                    };

                    // TODO (from Python) Consolidate checks for systems with/without primary pipework
                    if !matches!(
                        heat_source,
                        HeatSource::Storage(HeatSourceWithStorageTank::Immersion(_))
                    ) {
                        let (primary_pipework_losses_kwh, _) =
                            self.storage_tank.calculate_primary_pipework_losses(
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

    /// Return the DHW recoverable heat losses as internal gain for the current timestep in W
    pub(crate) fn internal_gains(&self) -> f64 {
        self.storage_tank.internal_gains()
    }

    pub(crate) fn output_results(&self) -> Option<Vec<StorageTankDetailedResult>> {
        self.storage_tank.output_results()
    }

    pub(crate) fn get_cold_water_source(&self) -> &WaterSourceWithTemperature {
        self.storage_tank.get_cold_water_source()
    }

    fn determine_heat_source_switch_on(
        &self,
        temp_s3_n: &[f64],
        heat_source_name: &str,
        heat_source: &HeatSource,
        _heater_layer: usize,
        _thermostat_layer: usize,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<()> {
        let (setpntmin, _) = self.retrieve_setpnt(heat_source, simtime)?;

        // Calculates state of charge
        let state_of_charge = self.calc_state_of_charge(temp_s3_n, simtime)?;

        // Turn heater on if state of charge is less than minimum state of charge
        if setpntmin.is_some_and(|setpntmin| state_of_charge <= setpntmin) {
            self.storage_tank.heating_active[heat_source_name].store(true, Ordering::SeqCst);
        }

        Ok(())
    }

    fn determine_heat_source_switch_off(
        &self,
        temp_s8_n: &[f64],
        heat_source_name: &str,
        _heater_layer: usize,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<()> {
        let heat_source = self.storage_tank.heat_source_data[heat_source_name]
            .heat_source
            .clone();
        let (_, setpntmax) = self.retrieve_setpnt(heat_source.lock().deref(), simtime)?;

        // Calculates state of charge
        let state_of_charge = self.calc_state_of_charge(temp_s8_n, simtime)?;

        // Turn heater off if max temp is None or state of charge has reached maximum state of charge
        if setpntmax.is_some_and(|setpntmax| state_of_charge >= setpntmax) {
            self.storage_tank.heating_active[heat_source_name].store(false, Ordering::SeqCst);
        }

        Ok(())
    }

    // Making this method return a Result as the corresponding method on StorageTank does, and in the original Python
    // SmartHotWaterTank subclasses StorageTank. We're making the assumption here that .setpnt() will always return a
    // `Some` value in normal functioning. If that isn't the case, we would need to address this differently.
    fn temp_flow(&self, simtime: SimulationTimeIteration) -> anyhow::Result<f64> {
        self.temp_setpnt_max.setpnt(&simtime).ok_or_else(|| {
            anyhow!("Expected to be able to access a setpoint value in SmartHotWaterTank")
        })
    }

    fn calc_state_of_charge(
        &self,
        t_h: &[f64],
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // Thermocline sensors calculate temperatures at all layers in the tank
        let number_of_layers = t_h.len();
        let height_of_layer = 1.0 / number_of_layers as f64;

        // Usable temperature
        let t_u = self.temp_usable;

        // Cold inlet temperature
        let list_temp_vol: Vec<(f64, f64)> = self
            .storage_tank
            .cold_feed
            .get_temp_cold_water(self.storage_tank.volume_total_in_litres, simtime);
        let sum_t_by_v: f64 = list_temp_vol.iter().map(|(t, v)| t * v).sum();
        let sum_v: f64 = list_temp_vol.iter().map(|(_t, v)| v).sum();
        let t_c = sum_t_by_v / sum_v;
        // TODO (from Python) Maybe use underlying cold feed?

        // Max set point temperature
        let t_sp = self
            .temp_setpnt_max
            .setpnt(&simtime)
            .unwrap_or(self.temp_usable);

        // Calculate state of charge
        let mut soc_numerator_total = 0.0;
        for &t_h_i in t_h {
            if t_h_i >= t_u {
                soc_numerator_total += (1. + (t_h_i - t_u) / (t_u - t_c)) * height_of_layer;
            }
        }
        let soc_denominator = 1. + (t_sp - t_u) / (t_u - t_c);

        // Rounding to avoid floating point errors
        let soc = round_by_precision(soc_numerator_total / soc_denominator, 1e5);

        // Raise error if below 0
        // TODO (from Python) add an error message if state of charge above 1 when function called.
        // The error should be raised when appropriate as there are instances when
        // the SOC can be above 1 which may not be invalid such as when temp_setpnt_max
        // is decreased from one timestep to another. To determine whether when it's
        // appropriate to call an error the soc from the pre timestep is needed
        // which is currently not recorded.
        if soc < 0.0 {
            bail!("State of charge should not be below 0, instead SOC is {soc}");
        }

        Ok(soc)
    }

    /// Calculate energy required to hit target temperature
    fn energy_req_target_temp(
        &self,
        initial_temps: &[f64],
        energy_available: &[f64],
        heater_layer: usize,
        q_ls_n_prev_heat_source: &[f64],
        temp_setpntmax: Option<f64>,
    ) -> (f64, Vec<f64>, Vec<f64>) {
        // Calculate temps with energy_available
        let (_, temp_after_input) = self
            .storage_tank
            .calc_temps_with_energy_input(initial_temps, energy_available);

        // Rearrange tank temperatures
        let (stored_heat, rearranged_temps) =
            self.storage_tank.rearrange_temperatures(&temp_after_input);

        let (energy_input, _q_ls, final_temps, q_ls_n) =
            self.storage_tank.calc_temps_after_thermal_losses(
                &rearranged_temps,
                energy_available,
                stored_heat,
                heater_layer,
                q_ls_n_prev_heat_source,
                temp_setpntmax,
            );

        (energy_input, final_temps, q_ls_n)
    }

    fn calculate_energy_for_state_of_charge(
        &self,
        heat_source: &HeatSource,
        initial_temps: &[f64],
        q_x_in_n: &[f64],
        heater_layer: usize,
        q_ls_n_prev_heat_source: &[f64],
        control_max_diverter: Option<&Control>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, Vec<f64>)> {
        let soc_max = if let Some(control_max_diverter) = control_max_diverter {
            control_max_diverter.setpnt(&simtime)
        } else {
            let (_, soc_max) = self.retrieve_setpnt(heat_source, simtime)?;
            soc_max
        };

        let mut temp_layers = initial_temps.to_vec();
        let mut energy_available = q_x_in_n.to_vec();
        let mut q_in_h_w = vec![0.; self.storage_tank.vol_n.len()];
        let mut q_ls_n_already_considered = q_ls_n_prev_heat_source.to_vec();

        for i in 0..self.storage_tank.vol_n.len() {
            if energy_available.iter().sum::<f64>() <= 0. {
                break;
            }

            // Calculate energy required for usable and max temperatures
            let (energy_req_usable, temp_simulation_usable, _q_ls_n_usable) = self
                .energy_req_target_temp(
                    &temp_layers,
                    &energy_available,
                    heater_layer,
                    &q_ls_n_already_considered,
                    self.temp_usable.into(),
                );
            let (energy_req_max, temp_simulation_max, q_ls_n_max) = self.energy_req_target_temp(
                &temp_layers,
                &energy_available,
                heater_layer,
                &q_ls_n_already_considered,
                self.temp_setpnt_max.setpnt(&simtime),
            );

            // Ensure energy required for usable and max temperatures are not negative
            let energy_req_usable = energy_req_usable.max(0.);
            let energy_req_max = energy_req_max.max(0.);

            // Calculate state of charge for usable and max temperatures
            let soc_temp_usable = self.calc_state_of_charge(&temp_simulation_usable, simtime)?;
            let soc_temp_max = self.calc_state_of_charge(&temp_simulation_max, simtime)?;
            if soc_max.is_some() {
                let soc_max = soc_max.unwrap();
                if soc_temp_usable >= soc_max {
                    q_in_h_w[i] += energy_req_usable;
                    break;
                } else if soc_temp_max > soc_max {
                    let energy_required_for_soc = np_interp(
                        soc_max,
                        &[soc_temp_usable, soc_temp_max],
                        &[energy_req_usable, energy_req_max],
                    );
                    q_in_h_w[heater_layer] += energy_required_for_soc;
                    break;
                }
            }

            q_ls_n_already_considered = q_ls_n_max
                .iter()
                .zip(q_ls_n_already_considered.iter())
                .map(|(x1, x2)| x1 + x2)
                .collect();
            q_in_h_w[heater_layer] += energy_req_max;
            energy_available[heater_layer] -= energy_req_max;
            temp_layers = temp_simulation_max;

            let energy_req_bottom_layer_to_setpnt = self.storage_tank.rho
                * self.storage_tank.cp
                * self.storage_tank.vol_n[0]
                * temp_layers[0];
            if energy_available.iter().sum::<f64>() <= energy_req_bottom_layer_to_setpnt {
                // Pump partial layer to the top
                let fraction_to_pump =
                    energy_available.iter().sum::<f64>() / energy_req_bottom_layer_to_setpnt;
                let volume_to_pump = fraction_to_pump * self.storage_tank.vol_n[0];
                let mut remaining_vols = self.storage_tank.vol_n.clone();
                temp_layers =
                    self.temps_after_pumping(volume_to_pump, &mut remaining_vols, &temp_layers);
            } else {
                // Pump one layer to the top
                let temp_pumped_layer = temp_layers.remove(0);
                temp_layers.push(temp_pumped_layer);
                let q_ls_layer = q_ls_n_already_considered.remove(0);
                q_ls_n_already_considered.push(q_ls_layer);
            }
        }

        // Calculate total energy required to meet max state of charge
        let energy_req_for_soc = q_in_h_w.iter().sum::<f64>();

        Ok((energy_req_for_soc, q_in_h_w))
    }

    fn calc_final_temps(
        &self,
        temp_s3_n: Vec<f64>,
        heat_source: &HeatSource,
        q_x_in_n: Vec<f64>,
        heater_layer: usize,
        q_ls_n_prev_heat_source: &[f64],
        control_max_diverter: Option<&Control>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<TemperatureCalculation> {
        let temp_setpntmax = self.temp_setpnt_max.setpnt(&simtime);

        // Tank with energy required for state of charge
        let (energy_req_for_soc, q_in_h_w_n) = self.calculate_energy_for_state_of_charge(
            heat_source,
            temp_s3_n.as_slice(),
            q_x_in_n.as_slice(),
            heater_layer,
            q_ls_n_prev_heat_source,
            control_max_diverter,
            simtime,
        )?;

        // Calculate temperatures after energy required to hit state of charge input
        let (q_s6, temp_s6_n) = self
            .storage_tank
            .calc_temps_with_energy_input(temp_s3_n.as_slice(), &q_in_h_w_n);

        // Rearrange tank
        let (_q_h_sto_s7, temp_s7_n) = self.storage_tank.rearrange_temperatures(&temp_s6_n);

        // Calculate new temperatures after operation of top up pump
        let temp_s7_n = self.calc_temps_after_top_up_pump(
            &temp_s7_n,
            energy_req_for_soc,
            heater_layer,
            simtime,
        )?;

        // Rearrange tank
        let (q_h_sto_s7, temp_s7_n) = self.storage_tank.rearrange_temperatures(&temp_s7_n);

        // STEP 8 Thermal losses and final temperature
        let (q_in_h_w, q_ls, temp_s8_n, q_ls_n) =
            self.storage_tank.calc_temps_after_thermal_losses(
                &temp_s7_n,
                &q_in_h_w_n,
                q_h_sto_s7,
                heater_layer,
                q_ls_n_prev_heat_source,
                temp_setpntmax,
            );

        // Adjust energy input based on actual usage
        let input_energy_adj = q_in_h_w;

        #[cfg(test)]
        self.storage_tank
            .energy_demand_test
            .store(input_energy_adj, Ordering::SeqCst);

        // Actual heat source output
        let heat_source_output = self.storage_tank.heat_source_output(
            heat_source,
            input_energy_adj,
            heater_layer,
            simtime,
            Some(self),
        )?;

        // calculate volume pumped using actual heat source output
        let volumes = self.storage_tank.vol_n.clone();
        let volume_pumped = self.bottom_to_top_pump_volume(
            temp_s3_n.as_slice(),
            heat_source_output,
            heater_layer,
            &volumes,
            simtime,
        )?;

        // Calculate pump energy consumption
        let energy_per_litre =
            self.power_pump_kw / (self.max_flow_rate_pump_l_per_min * MINUTES_PER_HOUR as f64);
        let pump_energy_kwh = energy_per_litre * volume_pumped;

        // Record pump energy consumption
        self.energy_supply_connection_pump
            .demand_energy(pump_energy_kwh, simtime.index)?;

        Ok(TemperatureCalculation {
            temp_s8_n,
            q_x_in_n: q_x_in_n.to_vec(),
            q_s6,
            temp_s6_n,
            temp_s7_n,
            q_in_h_w,
            q_ls,
            q_ls_n,
        })
    }

    /// Calculate new temperatures after top up pump of Smart hot water tank
    fn calc_temps_after_top_up_pump(
        &self,
        temp_s7_n: &[f64],
        q_x_in_n: f64,
        heater_layer: usize,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<f64>> {
        // Init for remaining volume of water in storage layers
        let mut remaining_vols = self.storage_tank.vol_n.clone();
        // Temperature of water in storage tank layers
        let tank_layer_temperatures = temp_s7_n.to_vec();

        // Volume pumped using top up pump
        let volume_pumped = self.bottom_to_top_pump_volume(
            &tank_layer_temperatures,
            q_x_in_n,
            heater_layer,
            &remaining_vols,
            simtime,
        )?;

        Ok(self.temps_after_pumping(volume_pumped, &mut remaining_vols, &tank_layer_temperatures))
    }

    /// Calculate the temperatures of the tank after volume is pumped
    fn temps_after_pumping(
        &self,
        volume_pumped: f64,
        remaining_vols: &mut [f64],
        tank_layer_temperatures: &[f64],
    ) -> Vec<f64> {
        let mut tank_layer_temperatures = tank_layer_temperatures.to_vec();
        let mut remaining_vols = remaining_vols.to_vec();
        if volume_pumped > 0. {
            // Calculate water removed
            // ---------------
            // If there is water to be pumped, remove water from bottom layers
            // starting from bottom layer. This will keep removing until there
            // is no more water to be removed.
            let mut volume_pumped_remaining = volume_pumped;
            for remaining_vol in remaining_vols.iter_mut() {
                if volume_pumped_remaining <= 0. {
                    break;
                }
                let volume_removed = volume_pumped_remaining.min(*remaining_vol);
                *remaining_vol -= volume_removed;
                volume_pumped_remaining -= volume_removed;
            }

            // Carry out water redistribution
            // ---------------
            // Iterate from the bottom layer upwards. Calculate the amount of
            // water needed to refill each layer
            for i in 0..self.storage_tank.vol_n.len() {
                // Determine how much volume needs to be added to this layer
                let mut needed_volume = self.storage_tank.vol_n[i] - remaining_vols[i];

                // If this layer is already full, continue to the next
                if needed_volume <= 0. {
                    continue;
                }

                // Initialise the variables for mixing temperatures
                let mut volume_weighted_temperature =
                    remaining_vols[i] * tank_layer_temperatures[i];

                // Filling layer
                // ---------------
                // For each layer that needs water, it looks at layer above
                // it to find available water. As it finds available water, it
                // moves it to the current layer and mixes the temperature.
                // Code allows circular movement of water where if it reaches the
                // top layer and still requires more water, it will circle back to
                // the bottom layer and check again.
                for mut j in (i + 1)..(i + self.storage_tank.vol_n.len()) {
                    if j >= self.storage_tank.vol_n.len() {
                        j -= self.storage_tank.vol_n.len();
                    }
                    // remaining_vols list is the volume of water available to replenish layer i
                    if remaining_vols[j] > 0. {
                        // Determine the volume to move down from this layer
                        let move_volume = needed_volume.min(remaining_vols[j]);
                        remaining_vols[j] -= move_volume;
                        remaining_vols[i] += move_volume;

                        // Adjust the temperature by mixing in the moved volume
                        volume_weighted_temperature += move_volume * tank_layer_temperatures[j];

                        // Decrease the amount of volume needed for the current layer
                        needed_volume -= move_volume;
                        if needed_volume <= 0. {
                            break;
                        }
                    }
                }

                debug_assert!(
                    remaining_vols[i] == self.storage_tank.vol_n[i],
                    "Volume mismatch in layer {i}"
                );

                // Temperatures after moving
                // ----------------
                // After moving water to a layer, calculate the new temperature
                // for the current layer based on vol and temperature.
                tank_layer_temperatures[i] = volume_weighted_temperature / remaining_vols[i];
            }
        }

        tank_layer_temperatures
    }

    /// Calculate the volume of water pumped from bottom to top of the tank
    fn bottom_to_top_pump_volume(
        &self,
        temp_s7_n: &[f64],
        qin: f64,
        heater_layer: usize,
        volumes: &[f64],
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // Initialise list of thermal losses in kWh
        let mut q_ls_n = vec![0.; self.storage_tank.vol_n.len()];

        // Standby losses coefficient - W/K
        let h_sto_ls = self.storage_tank.stand_by_losses_coefficient();

        let setpnt = self
            .temp_setpnt_max
            .setpnt(&simtime)
            .unwrap_or(self.temp_usable);

        let setpnt = self
            .temp_setpnt_max
            .setpnt(&simtime)
            .unwrap_or(self.temp_usable);

        // Calculate heat losses difference for all layers
        for (i, &vol_i) in self.storage_tank.vol_n.iter().enumerate() {
            q_ls_n[i] = (h_sto_ls * self.storage_tank.rho * self.storage_tank.cp)
                * (vol_i / self.storage_tank.volume_total_in_litres)
                * (setpnt - self.storage_tank.ambient_temperature)
                * simtime.timestep;
        }

        // The heat losses list is used to calculate the temperature difference
        // required for the top layer.
        let temp_diff_losses = *q_ls_n
            .iter()
            .last()
            .expect("q_ls_n was not expected to be empty")
            / (self.storage_tank.rho
                * self.storage_tank.cp
                * self
                    .storage_tank
                    .vol_n
                    .iter()
                    .last()
                    .expect("vol_n was not expected to be empty"));

        // Top layer temperature
        let top_layer_temp = *temp_s7_n
            .iter()
            .last()
            .expect("temp_s7_n was not expected to be empty");

        // Target temperature is increased to account for thermal losses.
        let temp_target = setpnt + temp_diff_losses;

        if top_layer_temp <= temp_target || qin <= 0. {
            // No pumping needed if top layer is below setpoint or no energy available
            return Ok(0.);
        }

        // Split volumes into below the heater layer
        let bottom_volumes = volumes[..heater_layer].to_vec();

        // Initialize fractions
        // 0 for layers below heater layer, 1 for heater layer and above
        let mut temp_factors = Vec::with_capacity(volumes.len());
        temp_factors.extend(iter::repeat_n(0., heater_layer));
        temp_factors.extend(iter::repeat_n(1., volumes.len() - heater_layer));

        // Iterates through the tank up to heater layer to determine how much each layer needs to be pumped
        for current_layer in 0..heater_layer {
            // Calculate the fraction of the current layer that needs to be pumped to maintain the overall temperature
            // Note: strictly speaking, the sums in the formula below should
            // exclude the current layer, but as the initial value of the temperature
            // factor is zero, this makes no difference in practice

            let numerator = temp_s7_n
                .iter()
                .zip(temp_factors.iter())
                .map(|(&t, &f)| t * f)
                .sum::<f64>()
                - temp_target * temp_factors.iter().copied().sum::<f64>();
            let denominator = temp_target - temp_s7_n[current_layer];
            temp_factors[current_layer] = if denominator <= 0. {
                // If the current layer is at or above target temperature, pump all of it
                1.0
            } else {
                // Calculate the fraction of the current layer to be pumped
                numerator / denominator
            };

            if temp_factors[current_layer] < 1. {
                // If we don't need to pump the entire layer, stop iteration
                break;
            } else {
                // If entire layer needs to be pumped, set factor to 1
                // and continue to next layer
                temp_factors[current_layer] = 1.0;
            }
        }

        // Calculate volume to be pumped (only from layers below heater_layer)
        let volume_pumped = bottom_volumes
            .iter()
            .zip(temp_factors[..heater_layer].iter())
            .map(|(&v, &f)| v * f)
            .sum::<f64>();

        // Check that the volume pumped doesn't exceed the volume of water up to the heater layer
        debug_assert!(volume_pumped <= bottom_volumes.iter().sum::<f64>());

        // Cap volume pumped based on pump max flow rate in timestep
        let max_volume_pumped = self.max_flow_rate_pump_l_per_min
            * self.storage_tank.simulation_timestep
            * MINUTES_PER_HOUR as f64;
        let volume_pumped = volume_pumped.min(max_volume_pumped);

        Ok(volume_pumped)
    }

    pub(crate) fn get_temp_hot_water(
        &self,
        volume_req: f64,
        volume_req_already: Option<f64>,
    ) -> Vec<(f64, f64)> {
        self.storage_tank
            .get_temp_hot_water(volume_req, volume_req_already)
    }

    pub(crate) fn draw_off_hot_water(
        &self,
        volume: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (Option<f64>, f64) {
        self.storage_tank
            .draw_off_hot_water(volume, simulation_time_iteration)
    }

    pub(crate) fn draw_off_water(&self,
        volume_needed: f64,
        simulation_time_iteration: SimulationTimeIteration) -> Vec<(f64, f64)> {
            self.storage_tank.draw_off_water(volume_needed, simulation_time_iteration)
        }

     pub(crate) fn get_temp_cold_water(&self, volume_needed: f64) -> Vec<(f64, f64)> {
        // TODO this matches Python - is it correct?
        self.storage_tank.get_temp_hot_water(volume_needed, None)
    }
}

#[derive(Debug)]
pub struct ImmersionHeater {
    pwr: f64, // rated power
    energy_supply_connection: EnergySupplyConnection,
    simulation_timestep: f64,
    control_min: Option<Arc<Control>>,
    control_max: Option<Arc<Control>>,
    diverter: ArcSwapOption<RwLock<PVDiverter>>,
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
        control_min: Option<Arc<Control>>,
        control_max: Option<Arc<Control>>,
    ) -> Self {
        Self {
            pwr: rated_power,
            energy_supply_connection,
            simulation_timestep,
            control_min,
            control_max,
            diverter: Default::default(),
        }
    }
    pub(crate) fn setpnt(&self, simtime: SimulationTimeIteration) -> (Option<f64>, Option<f64>) {
        (
            if self.control_min.is_some() {
                self.control_min.as_ref().unwrap().setpnt(&simtime)
            } else {
                None
            },
            if self.control_max.is_some() {
                self.control_max.as_ref().unwrap().setpnt(&simtime)
            } else {
                None
            },
        )
    }

    pub(crate) fn connect_diverter(&self, diverter: Arc<RwLock<PVDiverter>>) {
        if self.diverter.load_full().is_some() {
            panic!("diverter was already connected");
        }

        self.diverter.swap(Some(diverter));
    }

    /// Demand energy (in kWh) from the heater
    pub fn demand_energy(
        &self,
        energy_demand: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        if energy_demand < 0.0 {
            bail!("Negative energy demand on ImmersionHeater");
        };

        let energy_supplied =
            if self.control_min.is_some() && self.control_min.as_ref().unwrap().is_on(simtime) {
                min_of_2(energy_demand, self.pwr * self.simulation_timestep)
            } else {
                0.
            };

        // If there is a diverter to this immersion heater, then any heating
        // capacity already in use is not available to the diverter.
        if let Some(ref diverter) = &self.diverter.load_full() {
            diverter.read().increment_capacity_used(energy_supplied);
        }

        self.energy_supply_connection
            .demand_energy(energy_supplied, simtime.index)?;

        Ok(energy_supplied)
    }

    /// Calculate the maximum energy output (in kWh) from the heater
    pub fn energy_output_max(
        &self,
        simtime: SimulationTimeIteration,
        ignore_standard_control: bool,
    ) -> f64 {
        if self.control_min.is_some() && self.control_min.as_ref().unwrap().is_on(simtime)
            || ignore_standard_control
        {
            self.pwr * self.simulation_timestep
        } else {
            0.
        }
    }
}

/// Trait to represent a thing that can divert a surplus, like a PV diverter.
pub trait SurplusDiverting: Send + Sync {
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
    capacity_used: AtomicF64,
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
            capacity_used: Default::default(),
        }));

        heat_source.lock().connect_diverter(diverter.clone());

        diverter
    }

    /// Record heater output that would happen anyway, to avoid double-counting
    pub fn increment_capacity_used(&self, energy_supplied: f64) {
        self.capacity_used
            .fetch_add(energy_supplied, Ordering::SeqCst);
    }

    pub fn timestep_end(&self) {
        self.capacity_used
            .store(Default::default(), Ordering::SeqCst);
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
            - self.capacity_used.load(Ordering::SeqCst);

        // Calculate the maximum energy that could be diverted
        // Note: supply_surplus argument is negative by convention, so negate it here
        let energy_diverted_max = min_of_2(imm_heater_max_capacity_spare, -supply_surplus);

        // Add additional energy to storage tank and calculate how much energy was accepted

        let energy_diverted = match &self.pre_heated_water_source {
            HotWaterStorageTank::StorageTank(storage_tank) => {
                storage_tank.read().additional_energy_input(
                    &HeatSource::Storage(HeatSourceWithStorageTank::Immersion(
                        self.immersion_heater.clone(),
                    )),
                    &self.heat_source_name,
                    energy_diverted_max,
                    self.control_max.as_ref().map(|control| control.as_ref()),
                    simulation_time_iteration,
                )?
            }
            HotWaterStorageTank::SmartHotWaterTank(smart_hot_water_tank) => smart_hot_water_tank
                .read()
                .storage_tank
                .additional_energy_input(
                    &HeatSource::Storage(HeatSourceWithStorageTank::Immersion(
                        self.immersion_heater.clone(),
                    )),
                    &self.heat_source_name,
                    energy_diverted_max,
                    self.control_max.as_ref().map(|control| control.as_ref()),
                    simulation_time_iteration,
                )?,
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
    sol_loc: SolarCollectorLoopLocation,
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
    heat_output_collector_loop: AtomicF64,
    energy_supplied: AtomicF64,
    control_max: Arc<Control>,
    cp: f64,
    air_temp_coll_loop: AtomicF64,
    inlet_temp: AtomicF64,
    energy_supply_from_environment_conn: Option<EnergySupplyConnection>,
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
        sol_loc: SolarCollectorLoopLocation,
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
            heat_output_collector_loop: Default::default(),
            energy_supplied: Default::default(),
            control_max,
            // Water specific heat in J/kg.K
            // (defined under eqn 51 on page 40 of BS EN ISO 15316-4-3:2017)
            cp: contents.specific_heat_capacity(),
            air_temp_coll_loop: Default::default(),
            inlet_temp: Default::default(),
            energy_supply_from_environment_conn,
        }
    }

    pub(crate) fn energy_potential(&self) -> f64 {
        self.heat_output_collector_loop.load(Ordering::SeqCst)
    }

    pub(crate) fn setpnt(
        &self,
        simulation_time_iteration: &SimulationTimeIteration,
    ) -> (Option<f64>, Option<f64>) {
        let control_max_setpnt = self.control_max.setpnt(simulation_time_iteration);
        (control_max_setpnt, control_max_setpnt)
    }

    /// Calculate collector loop heat output
    /// eq 49 to 58 of STANDARD
    pub fn energy_output_max(
        &self,
        storage_tank: &StorageTank,
        temp_storage_tank_s3_n: &[f64],
        simulation_time: &SimulationTimeIteration,
    ) -> f64 {
        // Air temperature in a heated space in the building
        let air_temp_heated_room = (self.temp_internal_air_fn)();

        self.air_temp_coll_loop.store(
            match self.sol_loc {
                SolarCollectorLoopLocation::Hs => air_temp_heated_room,
                SolarCollectorLoopLocation::Nhs => {
                    (air_temp_heated_room + self.external_conditions.air_temp(simulation_time)) / 2.
                }
                SolarCollectorLoopLocation::Out => {
                    self.external_conditions.air_temp(simulation_time)
                }
            },
            Ordering::SeqCst,
        );

        // First estimation of average collector water temperature. Eq 51
        // initialise temperature
        // if first time step, pick bottom of the tank temperature as inlet_temp_s1
        let inlet_temp_s1 = if simulation_time.index == 0 {
            let inlet_temp_s1 = temp_storage_tank_s3_n[0];
            self.inlet_temp.store(inlet_temp_s1, Ordering::SeqCst);

            inlet_temp_s1
        } else {
            self.inlet_temp.load(Ordering::SeqCst)
        };

        // solar irradiance in W/m2
        let solar_irradiance = self.external_conditions.calculated_total_solar_irradiance(
            self.tilt,
            self.orientation,
            simulation_time,
        );

        if solar_irradiance == 0. {
            self.heat_output_collector_loop
                .store(Default::default(), Ordering::SeqCst);
            return 0.;
        }

        let mut avg_collector_water_temp = inlet_temp_s1
            + (0.4 * solar_irradiance * self.area) / (self.collector_mass_flow_rate * self.cp * 2.);

        // calculation of collector efficiency
        // Initialize inlet_temp2 before the loop using the initial inlet temperature
        let mut inlet_temp2: f64 = Default::default();
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

            // Eq 55
            let collector_output_heat =
                collector_efficiency * solar_irradiance * self.area * simulation_time.timestep
                    / WATTS_PER_KILOWATT as f64;

            // Eq 56
            let heat_loss_collector_loop_piping = self.solar_loop_piping_hlc
                * (avg_collector_water_temp - self.air_temp_coll_loop.load(Ordering::SeqCst))
                * simulation_time.timestep
                / WATTS_PER_KILOWATT as f64;

            // Eq 57
            self.heat_output_collector_loop.store(
                collector_output_heat - heat_loss_collector_loop_piping,
                Ordering::SeqCst,
            );
            if self.heat_output_collector_loop.load(Ordering::SeqCst)
                < self.power_pump * simulation_time.timestep * 3. / WATTS_PER_KILOWATT as f64
            {
                self.heat_output_collector_loop.store(0., Ordering::SeqCst);
            }

            // Call to the storage tank
            let (_temp_layer_0, inlet_temp2_temp) = storage_tank.storage_tank_potential_effect(
                self.heat_output_collector_loop.load(Ordering::SeqCst),
                temp_storage_tank_s3_n,
            );
            inlet_temp2 = inlet_temp2_temp;

            // Eq 58
            avg_collector_water_temp = (self.inlet_temp.load(Ordering::SeqCst) + inlet_temp2) / 2.
                + self.heat_output_collector_loop.load(Ordering::SeqCst)
                    / (self.collector_mass_flow_rate * self.cp * 2.);
        }

        self.inlet_temp.store(inlet_temp2, Ordering::SeqCst);

        self.heat_output_collector_loop.load(Ordering::SeqCst)
    }

    pub fn demand_energy(&self, energy_demand: f64, timestep_idx: usize) -> f64 {
        self.energy_supplied.store(
            min_of_2(
                energy_demand,
                self.heat_output_collector_loop.load(Ordering::SeqCst),
            ),
            Ordering::SeqCst,
        );

        // Eq 59 and 60 to calculate auxiliary energy - note that the if condition
        // is the wrong way round in BS EN 15316-4-3:2017
        let auxiliary_energy_consumption = if self.energy_supplied.load(Ordering::SeqCst) == 0. {
            self.power_pump_control * self.simulation_timestep
        } else {
            (self.power_pump_control + self.power_pump) * self.simulation_timestep
        };

        self.energy_supply_connection
            .demand_energy(auxiliary_energy_consumption, timestep_idx)
            .unwrap();

        if let Some(supply) = &self.energy_supply_from_environment_conn {
            supply
                .demand_energy(self.energy_supplied.load(Ordering::SeqCst), timestep_idx)
                .unwrap();
        }

        self.energy_supplied.load(Ordering::SeqCst)
    }

    #[cfg(test)]
    fn test_energy_potential(&self) -> f64 {
        self.heat_output_collector_loop.load(Ordering::SeqCst)
    }

    #[cfg(test)]
    fn test_energy_supplied(&self) -> f64 {
        self.energy_supplied.load(Ordering::SeqCst)
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
    use crate::input::{FuelType, PipeworkContents};
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
                end: 45.,
                ..Default::default()
            },
            ShadingSegment {
                start: 45.,
                end: 0.,
                shading_objects: vec![ShadingObject {
                    object_type: ShadingObjectType::Obstacle,
                    height: 10.5,
                    distance: 120.,
                }],
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
            Default::default(),
            Default::default(),
            simulation_timestep,
        );

        let control_max = SetpointTimeControl::new(
            control_max_schedule,
            0,
            1.,
            Default::default(),
            Default::default(),
            simulation_timestep,
        );
        let immersion_heater = ImmersionHeater::new(
            rated_power,
            energy_supply_connection.clone(),
            simulation_timestep,
            Some(Arc::new(Control::SetpointTime(control_min))),
            Some(Arc::new(Control::SetpointTime(control_max))),
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

        let heat_sources = IndexMap::from([("imheater".into(), heat_source)]);
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
            *WATER,
            None,
            None,
            None,
            false
        )
        .unwrap();

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

        let heat_sources = IndexMap::from([(String::from("imheater2"), heat_source)]);
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
            *WATER,
            None,
            None,
            None,
            false
        )
        .unwrap();

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
                end: 45.,
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

        let heat_sources = IndexMap::from([(String::from("imheater"), heat_source)]);

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
            *WATER,
            None,
            None,
            None,
            false,
        )
        .unwrap()
    }

    #[fixture]
    fn diverter_control() -> Arc<Control> {
        Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(60.), Some(60.), Some(60.), Some(60.)],
            0,
            1.,
            Default::default(),
            Default::default(),
            1.,
        ))
        .into()
    }

    #[rstest]
    fn test_divert_surplus(
        mut storage_tank_for_pv_diverter: StorageTank,
        immersion_heater: ImmersionHeater,
        diverter_control: Arc<Control>,
    ) {
        // _StorageTank__Q_ls_n_prev_heat_source is needed for the functions to
        // run the test but have no bearing in the results

        storage_tank_for_pv_diverter.q_ls_n_prev_heat_source =
            Arc::new(RwLock::new(vec![0.0, 0.1, 0.2, 0.3]));
        let pvdiverter = PVDiverter::new(
            &HotWaterStorageTank::StorageTank(Arc::new(RwLock::new(storage_tank_for_pv_diverter))),
            Arc::new(Mutex::new(immersion_heater)),
            "imheater".into(),
            diverter_control.into(),
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
                end: 45.,
                ..Default::default()
            },
            ShadingSegment {
                start: 45.,
                end: 0.,
                shading_objects: vec![ShadingObject {
                    object_type: ShadingObjectType::Obstacle,
                    height: 10.5,
                    distance: 12.,
                }],
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
            Default::default(),
            Default::default(),
            simulation_time_for_solar_thermal.step,
        );

        let solar_thermal = Arc::new(Mutex::new(SolarThermalSystem::new(
            SolarCollectorLoopLocation::Out,
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
            Arc::new(Control::SetpointTime(control_max)),
            *WATER,
            None,
        )));

        let storage_tank = StorageTank::new(
            150.0,
            1.68,
            55.0,
            cold_feed,
            simulation_time_for_solar_thermal.step,
            IndexMap::from([(
                String::from("solthermal"),
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
            *WATER,
            None,
            None,
            None,
            false
        )
        .unwrap();

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
        let (storage_tank1, energy_supply1) = storage_tank1;
        let (storage_tank2, energy_supply2) = storage_tank2;
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
                    volume_hot: None, // TODO
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
                    volume_hot: None, // TODO
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
                    volume_hot: None, // TODO
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
                    volume_hot: None, // TODO
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
                volume_hot: None, // TODO
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
                volume_hot: None, // TODO
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
                volume_hot: None, // TODO
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
                storage_tank1.temp_n.read().clone(),
                expected_temperatures_1[t_idx],
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
                storage_tank2.temp_n.read().clone(),
                expected_temperatures_2[t_idx],
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
            PipeworkContents::Water,
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
            PipeworkContents::Water,
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
        let expected = vec![(55.0, 37.5), (55.0, 37.5), (55.0, 25.0)];
        assert_eq!(storage_tank1.get_temp_hot_water(100.0, None), expected);
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
        let (storage_tank1, _) = storage_tank1;
        let temp_s3_n = [55.0, 55.0, 55.0, 55.0, 55.0, 55.0, 55.0, 55.0];
        let heat_source = storage_tank1.heat_source_data["imheater"]
            .clone()
            .heat_source;

        assert_eq!(
            storage_tank1
                .potential_energy_input(
                    &temp_s3_n,
                    &heat_source.lock(),
                    "imheater",
                    0,
                    7,
                    simulation_time_for_storage_tank.iter().current_iteration()
                )
                .unwrap(),
            [0.0, 0., 0., 0.]
        );

        // SolarThermal as heat source
        let (storage_tank_solar_thermal, _, simtime, _) = storage_tank_with_solar_thermal;
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
                .potential_energy_input(&temp_s3_n, &heat_source.lock(), "solthermal", 0, 7, t_it)
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

        let result = storage_tank1.calc_temps_with_energy_input(&temp_s3_n, &q_x_in_n);

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
        let (storage_tank1, _) = storage_tank1;
        let temp_s7_n = [12.0, 18.0, 25.0, 32.0, 37.0, 45.0, 49.0, 58.0];
        let q_x_in_n = [0., 1., 2., 3., 4., 5., 6., 7., 8.];
        let q_h_sto_s7 = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8];
        let heater_layer = 2;
        let q_ls_n_prev_heat_source = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let setpntmax = 55.0;

        assert_eq!(
            storage_tank1.calc_temps_after_thermal_losses(
                &temp_s7_n,
                &q_x_in_n,
                q_h_sto_s7,
                heater_layer,
                &q_ls_n_prev_heat_source,
                setpntmax.into()
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
        let (storage_tank1, _) = storage_tank1;
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
                    &heat_source.lock(),
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
        let (storage_tank1, _) = storage_tank1;
        let temp_s3_n = vec![10., 15., 20., 25., 25., 30., 35., 50.];
        let heat_source = storage_tank1.heat_source_data["imheater"]
            .clone()
            .heat_source;
        let q_x_in_n = vec![0., 1., 2., 3., 4., 5., 6., 7., 8.];
        let heater_layer = 2;
        let q_ls_n_prev_heat_source = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

        assert_eq!(
            storage_tank1
                .calc_final_temps(
                    &temp_s3_n,
                    &heat_source.lock(),
                    q_x_in_n,
                    heater_layer,
                    &q_ls_n_prev_heat_source,
                    simulation_time_for_storage_tank.iter().current_iteration(),
                    None
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
    fn test_extract_hot_water(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        let (storage_tank1, _) = storage_tank1;
        storage_tank1
            .temp_average_drawoff_volweighted
            .store(0.0, Ordering::SeqCst);
        storage_tank1
            .total_volume_drawoff
            .store(0.0, Ordering::SeqCst);
        let event = TypedScheduleEvent {
            start: 0.,
            duration: Some(1.),
            temperature: 41.0,
            event_type: WaterScheduleEventType::Other,
            name: "other".to_string(),
            volume: None,
            warm_volume: Some(8.0),
            pipework_volume: Some(5.),
            volume_hot: Some(5.511111111111113)
        };

        assert_eq!(
            storage_tank1.extract_hot_water(
                event,
                simulation_time_for_storage_tank.iter().current_iteration()
            ).unwrap(),
            (
                5.51111111111112,
                0.2882311111111112,
                vec![37.5, 37.5, 32.5, 31.988888888888887]
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
            storage_tank1.calc_temps_after_extraction(
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
        storage_tank1.q_ls_n_prev_heat_source = Arc::new(RwLock::new(vec![0.0, 0.1, 0.2, 0.3]));
        assert_eq!(
            storage_tank1
                .additional_energy_input(
                    &heat_source.lock(),
                    "imheater",
                    energy_input,
                    None,
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
        storage_tank2.q_ls_n_prev_heat_source = Arc::new(RwLock::new(vec![0.0, 0.1, 0.2, 0.3]));
        assert_eq!(
            storage_tank2
                .additional_energy_input(
                    &heat_source.lock(),
                    "imheater2",
                    energy_input,
                    None,
                    simulation_time_for_storage_tank.iter().current_iteration()
                )
                .unwrap(),
            0.0
        );
    }

    #[rstest]
    fn test_internal_gains(storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>)) {
        let (storage_tank1, _) = storage_tank1;
        storage_tank1.q_sto_h_ls_rbl.store(0.05, Ordering::SeqCst);

        assert_eq!(storage_tank1.internal_gains(), 50.);
    }

    #[ignore = "work in progress - pipework migration"]
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
                PipeworkContents::Water,
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
                PipeworkContents::Water,
            )
            .unwrap(),
        ];

        storage_tank1.primary_pipework = Some(primary_pipework_lst);

        for (t_idx, t_it) in simulation_time_for_storage_tank.iter().enumerate() {
            assert_eq!(
                storage_tank1.calculate_primary_pipework_losses(input_energy_adj, setpnt_max, t_it),
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
                storage_tank1.calculate_primary_pipework_losses(input_energy_adj, setpnt_max, t_it),
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

        let control_min = Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(52.), Some(52.), None, Some(52.)],
            0,
            1.,
            Default::default(),
            Default::default(),
            timestep,
        )));

        let control_max = Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            vec![Some(60.), Some(60.), Some(60.), Some(60.)],
            0,
            1.,
            Default::default(),
            Default::default(),
            timestep,
        )));

        ImmersionHeater::new(
            rated_power,
            energy_supply_connection,
            timestep,
            Some(control_min),
            Some(control_max),
        )
    }

    #[rstest]
    fn test_demand_energy_for_immersion_heater(
        immersion_heater: ImmersionHeater,
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

        solar_thermal.lock().sol_loc = SolarCollectorLoopLocation::Nhs;
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

        solar_thermal.lock().sol_loc = SolarCollectorLoopLocation::Hs;
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

    #[fixture]
    fn simulation_time_for_smart_hot_water_tank() -> SimulationTime {
        SimulationTime::new(0., 8., 1.)
    }

    #[fixture]
    fn simulation_time_iteration_for_smart_hot_water_tank(
        simulation_time_for_smart_hot_water_tank: SimulationTime,
    ) -> SimulationTimeIteration {
        simulation_time_for_smart_hot_water_tank
            .iter()
            .current_iteration()
    }

    #[fixture]
    fn energy_supply_for_smart_hot_water_tank_immersion(
        simulation_time_for_smart_hot_water_tank: SimulationTime,
    ) -> Arc<RwLock<EnergySupply>> {
        Arc::from(RwLock::from(
            EnergySupplyBuilder::new(
                FuelType::Electricity,
                simulation_time_for_smart_hot_water_tank.total_steps(),
            )
            .build(),
        ))
    }

    #[fixture]
    fn energy_supply_for_smart_hot_water_tank_pump(
        simulation_time_for_smart_hot_water_tank: SimulationTime,
    ) -> Arc<RwLock<EnergySupply>> {
        Arc::from(RwLock::from(
            EnergySupplyBuilder::new(
                FuelType::Electricity,
                simulation_time_for_smart_hot_water_tank.total_steps(),
            )
            .build(),
        ))
    }

    #[fixture]
    fn external_conditions_for_smart_hot_water_tank(
        external_conditions_for_pv_diverter: Arc<ExternalConditions>,
    ) -> Arc<ExternalConditions> {
        external_conditions_for_pv_diverter // external_conditions_for_pv_diverter has the same data & set up as what we need for smart hot water tank
    }

    fn create_smart_hot_water_tank(
        simulation_time_for_smart_hot_water_tank: SimulationTime,
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions_for_smart_hot_water_tank: Arc<ExternalConditions>,
        energy_supply_for_smart_hot_water_tank_immersion: Arc<RwLock<EnergySupply>>,
        energy_supply_for_smart_hot_water_tank_pump: Arc<RwLock<EnergySupply>>,
        heat_source_name: &str,
        volume: f64,
        losses: f64,
        init_temp: f64,
    ) -> SmartHotWaterTank {
        let power_pump_kw = 5.;
        let max_flow_rate_pump_l_per_min = 1000.;
        let temp_usable = 40.;
        let temp_setpnt_max = Arc::new(Control::SetpointTime(SetpointTimeControl::new(
            vec![
                Some(50.0),
                Some(40.0),
                Some(30.0),
                Some(20.0),
                Some(50.0),
                Some(50.0),
                Some(50.0),
                Some(50.0),
            ],
            0,
            1.,
            Default::default(),
            0.,
            1.,
        )));
        let cold_water_temps = vec![10.0, 10.1, 10.2, 10.5, 10.6, 11.0, 11.5, 12.1];
        let cold_feed = WaterSourceWithTemperature::ColdWaterSource(Arc::new(
            ColdWaterSource::new(cold_water_temps, 0, 1.),
        ));
        let energy_supply_connection = EnergySupply::connection(
            energy_supply_for_smart_hot_water_tank_immersion.clone(),
            heat_source_name,
        )
        .unwrap();
        let control_min = Control::SetpointTime(SetpointTimeControl::new(
            vec![
                Some(0.5),
                None,
                None,
                None,
                Some(0.5),
                Some(0.5),
                Some(0.5),
                Some(0.5),
            ],
            0,
            1.,
            Default::default(),
            0.,
            1.,
        ));
        let control_max = Control::SetpointTime(SetpointTimeControl::new(
            vec![
                Some(1.0),
                Some(1.0),
                Some(0.9),
                Some(0.8),
                Some(0.7),
                Some(1.0),
                Some(0.9),
                Some(0.8),
            ],
            0,
            1.,
            Default::default(),
            0.,
            1.,
        ));
        let immersion_heater = ImmersionHeater::new(
            5.,
            energy_supply_connection,
            1.,
            Some(control_min.into()),
            Some(control_max.into()),
        );
        let heat_source = HeatSource::Storage(HeatSourceWithStorageTank::Immersion(Arc::new(
            Mutex::new(immersion_heater),
        )));
        let heat_sources = IndexMap::from([(
            String::from(heat_source_name),
            PositionedHeatSource {
                heat_source: Arc::new(Mutex::new(heat_source)),
                heater_position: 0.6,
                thermostat_position: None,
            },
        )]);

        let energy_supply_conn_pump =
            EnergySupply::connection(energy_supply_for_smart_hot_water_tank_pump.clone(), "pump")
                .unwrap(); // N.B. this is a MagicMock in Python

        SmartHotWaterTank::new(
            volume,
            losses,
            init_temp,
            power_pump_kw,
            max_flow_rate_pump_l_per_min,
            temp_usable,
            temp_setpnt_max,
            cold_feed,
            simulation_time_for_smart_hot_water_tank
                .iter()
                .step_in_hours(),
            heat_sources,
            temp_internal_air_fn,
            external_conditions_for_smart_hot_water_tank,
            None,
            Some(4),
            None,
            energy_supply_conn_pump,
            None,
        )
        .unwrap()
    }

    #[fixture]
    fn smart_hot_water_tank(
        simulation_time_for_smart_hot_water_tank: SimulationTime,
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions_for_smart_hot_water_tank: Arc<ExternalConditions>,
        energy_supply_for_smart_hot_water_tank_immersion: Arc<RwLock<EnergySupply>>,
        energy_supply_for_smart_hot_water_tank_pump: Arc<RwLock<EnergySupply>>,
    ) -> SmartHotWaterTank {
        create_smart_hot_water_tank(
            simulation_time_for_smart_hot_water_tank,
            temp_internal_air_fn,
            external_conditions_for_smart_hot_water_tank,
            energy_supply_for_smart_hot_water_tank_immersion,
            energy_supply_for_smart_hot_water_tank_pump,
            "imheater",
            300.0,
            1.68,
            50.,
        )
    }

    // TODO use a new WaterEventResult type instead of TypedScheduleEvent
    fn get_event_data_immersion() -> Vec<Option<Vec<TypedScheduleEvent>>> {
        &[
            Some(vec![
                TypedScheduleEvent {
                    start: 6.,
                    duration: Some(6.),
                    temperature: 41.,
                    name: "IES".into(),
                    event_type: WaterScheduleEventType::Shower,
                    volume: None,
                    warm_volume: Some(48.0),
                    volume_hot: Some(33.0666666666667),
                    pipework_volume: None,
                },
                TypedScheduleEvent {
                    start: 6.,
                    duration: Some(6.),
                    temperature: 41.,
                    name: "mixer".into(),
                    event_type: WaterScheduleEventType::Shower,
                    volume: None,
                    warm_volume: Some(48.),
                    pipework_volume: None,
                },
                TypedScheduleEvent {
                    start: 6.,
                    duration: Some(20.),
                    temperature: 43.,
                    name: "medium".into(),
                    event_type: WaterScheduleEventType::Bath,
                    volume: None,
                    warm_volume: Some(100.),
                    pipework_volume: None,
                },
                TypedScheduleEvent {
                    start: 6.,
                    duration: Some(1.),
                    temperature: 40.,
                    name: "other".into(),
                    event_type: WaterScheduleEventType::Other,
                    volume: None,
                    warm_volume: Some(8.),
                    pipework_volume: None,
                },
            ]),
            Some(vec![TypedScheduleEvent {
                start: 7.,
                duration: Some(6.),
                temperature: 41.,
                name: "mixer".into(),
                event_type: WaterScheduleEventType::Shower,
                volume: None,
                warm_volume: Some(48.),
                pipework_volume: None,
            }]),
            None,
            Some(vec![TypedScheduleEvent {
                start: 9.,
                duration: Some(6.),
                temperature: 45.,
                name: "mixer".into(),
                event_type: WaterScheduleEventType::Shower,
                volume: None,
                warm_volume: Some(48.),
                pipework_volume: None,
            }]),
            None,
            Some(vec![TypedScheduleEvent {
                start: 11.,
                duration: Some(6.5),
                temperature: 41.,
                name: "mixer".into(),
                event_type: WaterScheduleEventType::Shower,
                volume: None,
                warm_volume: Some(52.),
                pipework_volume: None,
            }]),
            None,
            None,
        ];
    }

    const TWO_DECIMAL_PLACES: f64 = 1e-3;
    const FIVE_DECIMAL_PLACES: f64 = 1e-6;

    #[rstest]
    fn test_calc_state_of_charge_for_smart_hot_water_tank(
        smart_hot_water_tank: SmartHotWaterTank,
        simulation_time_iteration_for_smart_hot_water_tank: SimulationTimeIteration,
    ) {
        let t_h = [
            43.984858220267675,
            43.984858220267675,
            43.984858220267675,
            43.984858220267725,
            43.984858220267725,
            43.98485822026773,
            43.984858220267775,
            43.984858220267775,
        ];
        let soc = smart_hot_water_tank
            .calc_state_of_charge(&t_h, simulation_time_iteration_for_smart_hot_water_tank);
        assert_relative_eq!(soc.unwrap(), 0.850, max_relative = TWO_DECIMAL_PLACES);
    }

    #[rstest]
    fn test_calc_state_of_charge_low_high_for_smart_hot_water_tank(
        smart_hot_water_tank: SmartHotWaterTank,
        simulation_time_iteration_for_smart_hot_water_tank: SimulationTimeIteration,
    ) {
        let t_h_low = [10., 20., 25., 30., 35., 35., 35., 35.];
        let t_h_high = [50., 50., 50., 50., 50., 50., 50., 50.];
        let soc_low = smart_hot_water_tank
            .calc_state_of_charge(&t_h_low, simulation_time_iteration_for_smart_hot_water_tank);
        let soc_high = smart_hot_water_tank.calc_state_of_charge(
            &t_h_high,
            simulation_time_iteration_for_smart_hot_water_tank,
        );

        assert_eq!(soc_low.unwrap(), 0.0);
        assert_eq!(soc_high.unwrap(), 1.0);
    }

    #[rstest]
    fn test_calc_state_of_charge_varied_for_smart_hot_water_tank(
        smart_hot_water_tank: SmartHotWaterTank,
        simulation_time_iteration_for_smart_hot_water_tank: SimulationTimeIteration,
    ) {
        let t_h_high = [50., 40., 30., 20., 50., 50., 50., 50.];
        let soc_high = smart_hot_water_tank.calc_state_of_charge(
            &t_h_high,
            simulation_time_iteration_for_smart_hot_water_tank,
        );
        assert_eq!(soc_high.unwrap(), 0.71875);
    }

    #[rstest]
    fn test_bottom_to_top_pump_volume_no_pumping_for_smart_hot_water_tank(
        smart_hot_water_tank: SmartHotWaterTank,
        simulation_time_iteration_for_smart_hot_water_tank: SimulationTimeIteration,
    ) {
        let temp_s6_n = &[10., 20., 25., 30., 35., 35., 35., 35.];
        let qin = 0.;
        let heater_layer = 7;
        let volumes = &[10., 10., 10., 10., 10., 10., 10., 10.];
        let volume_pumped = smart_hot_water_tank.bottom_to_top_pump_volume(
            temp_s6_n,
            qin,
            heater_layer,
            volumes,
            simulation_time_iteration_for_smart_hot_water_tank,
        );
        assert_eq!(volume_pumped.unwrap(), 0.);
    }

    #[rstest]
    fn test_soc_over_time_for_smart_hot_water_tank(
        smart_hot_water_tank: SmartHotWaterTank,
        simulation_time_for_smart_hot_water_tank: SimulationTime,
    ) {
        let expected_soc = &[
            0.5625, 0.75042, 1.13005, 2.33553, 0.56155, 0.5609, 0.56006, 0.55904,
        ];
        let t_h_high = &[50., 40., 30., 20., 30., 40., 50., 50.];

        for (t_idx, t_it) in simulation_time_for_smart_hot_water_tank.iter().enumerate() {
            let soc = smart_hot_water_tank.calc_state_of_charge(t_h_high, t_it);
            assert_relative_eq!(
                soc.unwrap(),
                expected_soc[t_idx],
                max_relative = FIVE_DECIMAL_PLACES
            );
        }
    }

    #[rstest]
    fn test_demand_hot_water_for_smart_hot_water_tank(
        simulation_time_for_smart_hot_water_tank: SimulationTime,
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions_for_smart_hot_water_tank: Arc<ExternalConditions>,
    ) {
        let energy_supply_for_smart_hot_water_tank_immersion_1 = Arc::from(RwLock::from(
            EnergySupplyBuilder::new(
                FuelType::Electricity,
                simulation_time_for_smart_hot_water_tank.total_steps(),
            )
            .build(),
        ));
        let energy_supply_for_smart_hot_water_tank_pump_1 = Arc::from(RwLock::from(
            EnergySupplyBuilder::new(
                FuelType::Electricity,
                simulation_time_for_smart_hot_water_tank.total_steps(),
            )
            .build(),
        ));

        let energy_supply_for_smart_hot_water_tank_immersion_2 = Arc::from(RwLock::from(
            EnergySupplyBuilder::new(
                FuelType::Electricity,
                simulation_time_for_smart_hot_water_tank.total_steps(),
            )
            .build(),
        ));
        let energy_supply_for_smart_hot_water_tank_pump_2 = Arc::from(RwLock::from(
            EnergySupplyBuilder::new(
                FuelType::Electricity,
                simulation_time_for_smart_hot_water_tank.total_steps(),
            )
            .build(),
        ));

        let smart_hot_water_tank = create_smart_hot_water_tank(
            simulation_time_for_smart_hot_water_tank,
            temp_internal_air_fn.clone(),
            external_conditions_for_smart_hot_water_tank.clone(),
            energy_supply_for_smart_hot_water_tank_immersion_1.clone(),
            energy_supply_for_smart_hot_water_tank_pump_1.clone(),
            "imheater",
            300.,
            1.68,
            50.,
        );

        let smart_hot_water_tank_2 = create_smart_hot_water_tank(
            simulation_time_for_smart_hot_water_tank,
            temp_internal_air_fn,
            external_conditions_for_smart_hot_water_tank,
            energy_supply_for_smart_hot_water_tank_immersion_2.clone(),
            energy_supply_for_smart_hot_water_tank_pump_2.clone(),
            "immersion2",
            210.,
            1.61,
            60.,
        );

        let usage_events = get_event_data_immersion();

        let expected_temperatures_1 = &[
            vec![42.06224020071341, 50.0, 50.0, 50.0],
            vec![
                26.167007407407407,
                45.94555555555556,
                49.87555555555556,
                49.87555555555556,
            ],
            vec![
                26.11428959122085,
                45.872962962962966,
                49.802962962962965,
                49.802962962962965,
            ],
            vec![
                17.333051851851852,
                34.74925925925926,
                47.57925925925926,
                49.779259259259256,
            ],
            vec![31.873616537191992, 50.0, 50.0, 50.0],
            vec![
                20.71542222222222,
                40.20384444444444,
                49.82370370370371,
                49.82370370370371,
            ],
            vec![
                20.69097188477366,
                40.07834302880658,
                49.64832153635117,
                49.64832153635117,
            ],
            vec![
                20.666648326852613,
                39.9534923612498,
                49.47384875801453,
                49.47384875801453,
            ],
        ];

        let expected_temperatures_2 = [
            vec![
                10.0,
                24.55880864197531,
                50.03864197530864,
                59.08864197530864,
            ],
            vec![
                10.06,
                16.158864197530864,
                35.20270987654321,
                53.69962962962963,
            ],
            vec![
                10.06,
                16.157736457857034,
                35.103327160493826,
                53.60024691358024,
            ],
            vec![10.38, 11.7, 21.211604938271602, 40.031604938271606],
            vec![
                11.525423857571475,
                48.64021695897895,
                48.64021695897895,
                48.64021695897895,
            ],
            vec![50.0, 50.0, 50.0, 50.0],
            vec![
                49.75864197530864,
                49.75864197530864,
                49.75864197530864,
                49.75864197530864,
            ],
            vec![
                49.518997294619716,
                49.518997294619716,
                49.518997294619716,
                49.518997294619716,
            ],
        ];

        let expected_results_by_end_user_1 =
            [2.2100280434, 0.0, 0.0, 0.0, 2.0949861665, 0.0, 0.0, 0.0];

        let expected_results_by_end_user_2 =
            [0.0, 0.0, 0.0, 0.0, 4.5048003443, 0.1954190576, 0.0, 0.0];

        for (t_idx, t_it) in simulation_time_for_smart_hot_water_tank.iter().enumerate() {
            let _ = smart_hot_water_tank
                .demand_hot_water(usage_events[t_idx].clone(), t_it)
                .unwrap();

            let temp_n = smart_hot_water_tank.storage_tank.temp_n.read();
            for (i, expected_temp) in expected_temperatures_1[t_idx].iter().enumerate() {
                // note - high max_relative here. We expect differences due to Emitters
                assert_relative_eq!(temp_n[i], *expected_temp, max_relative = 0.03);
                // 3% rel tolerance
            }

            let results_by_end_user = energy_supply_for_smart_hot_water_tank_immersion_1
                .read()
                .results_by_end_user();
            let actual_results_by_end_user_1 = results_by_end_user.get("imheater").unwrap();
            assert_relative_eq!(
                actual_results_by_end_user_1[t_idx],
                expected_results_by_end_user_1[t_idx],
                max_relative = 0.05
            ); // 5% rel tolerance

            // smart_hot_water_tank_2 tests for case where heater does not heat all layers
            let _ = smart_hot_water_tank_2
                .demand_hot_water(usage_events[t_idx].clone(), t_it)
                .unwrap();

            let temp_n = smart_hot_water_tank_2.storage_tank.temp_n.read();
            for (i, expected_temp) in expected_temperatures_2[t_idx].iter().enumerate() {
                assert_relative_eq!(temp_n[i], *expected_temp, max_relative = 0.0019);
                // 0.19% rel tolerance
            }

            let results_by_end_user_2 = energy_supply_for_smart_hot_water_tank_immersion_2
                .read()
                .results_by_end_user();
            let actual_results_by_end_user_2 = results_by_end_user_2.get("immersion2").unwrap();
            assert_relative_eq!(
                actual_results_by_end_user_2[t_idx],
                expected_results_by_end_user_2[t_idx],
                max_relative = 0.007
            ); // 0.7% rel tolerance
        }
    }

    #[rstest]
    fn test_calc_temps_after_extraction(
        smart_hot_water_tank: SmartHotWaterTank,
        simulation_time_iteration_for_smart_hot_water_tank: SimulationTimeIteration,
    ) {
        let remaining_vol = vec![0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4.];
        let actual_remaining_vol = smart_hot_water_tank
            .storage_tank
            .calc_temps_after_extraction(
                remaining_vol,
                simulation_time_iteration_for_smart_hot_water_tank,
            );
        assert_eq!(actual_remaining_vol, (vec![10.0, 10.0, 10.0, 12.67], false));
    }
}
