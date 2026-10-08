use crate::compare_floats::{max_of_2, min_of_2};
#[cfg(test)]
use crate::core::common::MockWaterSupply;
use crate::core::common::{WaterSupply, WaterSupplyBehaviour};
use crate::core::controls::time_control::{Control, ControlBehaviour, RangeTimeControl};
use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::heating_systems::primary_pipework_losses_mixin::PrimaryPipeworkLossesMixin;
use crate::core::material_properties::{MaterialProperties, WATER};
use crate::core::pipework::Pipework;
use crate::core::units::{Orientation360, MINUTES_PER_HOUR, WATTS_PER_KILOWATT};
#[cfg(test)]
use crate::core::water_heat_demand::dhw_demand::tests::HotWaterSourceMockKind;
use crate::core::water_heat_demand::misc::{summarise_events, WaterEventResult};
use crate::corpus::{HeatSource, HotWaterSourceBehaviour, TempInternalAirFn};
use crate::external_conditions::ExternalConditions;
use crate::input::{SolarCollectorLoopLocation, WaterPipework};
use crate::simulation_time::SimulationTimeIteration;
use crate::StringOrNumber;
use anyhow::{anyhow, bail};
use approx::relative_eq;
use arc_swap::ArcSwapOption;
use arcstr::ArcStr;
use atomic_float::AtomicF64;
use educe::Educe;
use fsum::FSum;
use indexmap::IndexMap;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use parking_lot::{Mutex, RwLock};
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
pub(crate) const DEFAULT_AMBIENT_TEMPERATURE: f64 = 16.;

// Primary pipework gains for the timestep
const DEFAULT_PIPEWORK_PRIMARY_GAINS_FOR_TIMESTEP: f64 = 0.;

// Time of finalisation of the previous hot water event
const DEFAULT_PREVIOUS_EVENT_TIME_END: f64 = 0.;

// Auxiliary energy recovery factor
const THERMAL_CONSTANTS_F_RVD_AUX: f64 = 0.25;

// Thermal loss recovery factor
pub const THERMAL_CONSTANTS_F_STO_M: f64 = 0.75;

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
#[derive(Educe)]
#[educe(Debug)]
pub struct StorageTank {
    initial_temperature: f64,
    q_std_ls_ref: f64, // measured standby losses due to cylinder insulation at standardised conditions, in kWh/24h
    cold_feed: WaterSupply,
    simulation_timestep: f64,
    number_of_volumes: usize,
    temp_flow_prev: Arc<RwLock<Option<f64>>>,
    #[educe(Debug(ignore))]
    temp_internal_air_fn: TempInternalAirFn,
    volume_total_in_litres: f64,
    vol_n: Vec<f64>,
    cp: f64,  // contents (usually water) specific heat in kWh/kg.K
    rho: f64, // volumic mass in kg/litre
    temp_n: Arc<RwLock<Vec<f64>>>,
    primary_pipework_losses_kwh: AtomicF64,
    storage_losses_kwh: AtomicF64,
    heat_source_data: IndexMap<ArcStr, PositionedHeatSource>, // heat sources, sorted by heater position
    heating_active: IndexMap<ArcStr, AtomicBool>,
    q_ls_n_prev_heat_source: Arc<RwLock<Vec<f64>>>,
    q_sto_h_ls_rbl: AtomicF64, // total recoverable heat losses for heating in kWh, memoised between steps
    pipework_primary_gains_for_timestep: AtomicF64, // primary pipework gains for a timestep (mutates over lifetime)
    #[cfg(test)]
    energy_demand_test: AtomicF64,
    temp_final_drawoff: AtomicF64, // In Python this is created from inside extract_hot_water()
    temp_average_drawoff: AtomicF64, // In Python this is created from inside extract_hot_water()
    temp_average_drawoff_volweighted: AtomicF64, // In Python this is created from inside extract_hot_water()
    total_volume_drawoff: AtomicF64, // In Python this is created from inside extract_hot_water()
    ambient_temperature: f64,
    detailed_results: Option<Arc<RwLock<Vec<Vec<StringOrNumber>>>>>,
    pipework: PrimaryPipeworkLossesMixin,
}

impl StorageTank {
    /// Arguments:
    /// * `volume` - total volume of the tank, in litres
    /// * `losses` - measured standby losses due to cylinder insulation
    ///                                at standardised conditions, in kWh/24h
    /// * `init_temp` - initial temperature required for DHW
    /// * `cold_feed` - reference to ColdWaterSource object
    /// * `simulation_time_iteration` - reference to SimulationTime iteration
    /// * `heat_sources`     -- hashmap of names and heat source objects
    /// * `number_of_volumes` -number of volumes the storage is modelled with
    ///              see App.C (C.1.2 selection of the number of volumes to model the storage unit)
    ///              for more details if this wants to be changed.
    /// * `primary_pipework` - optional reference to pipework
    /// * `contents` - MaterialProperties object
    pub(crate) fn new(
        volume: f64,
        losses: f64,
        initial_temperature: f64,
        cold_feed: WaterSupply,
        simulation_time_iteration: &SimulationTimeIteration,
        heat_sources: IndexMap<ArcStr, PositionedHeatSource>,
        // In Python this is "project" but only temp_internal_air is accessed from it
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions: Arc<ExternalConditions>,
        detailed_output: bool,
        number_of_volumes: Option<usize>,
        primary_pipework_lst: Option<&Vec<WaterPipework>>,
        contents: MaterialProperties,
        ambient_temperature: Option<f64>,
        pipework_primary_gains_for_timestep: Option<f64>,
        previous_event_time_end: Option<f64>,
    ) -> anyhow::Result<Self> {
        let q_std_ls_ref = losses;
        let ambient_temperature = ambient_temperature.unwrap_or(DEFAULT_AMBIENT_TEMPERATURE);
        let pipework_primary_gains_for_timestep = pipework_primary_gains_for_timestep
            .unwrap_or(DEFAULT_PIPEWORK_PRIMARY_GAINS_FOR_TIMESTEP);
        let _previous_event_time_end =
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

        if !heat_sources.is_empty() {
            // Disallow multiple heat sources until per-heat-source pipework is modelled.
            let wet_heat_source_count = heat_sources
                .values()
                .filter(|source| matches!(*source.heat_source.lock(), HeatSource::Wet(_)))
                .count();
            if wet_heat_source_count > 1 {
                bail!("Only one wet heat source is allowed on a storage tank")
            }
        }

        // Build Pipework objects from raw input data, then initialise the
        // mixin state (event tracking, surrounding temperatures)
        let mut pipework_lst: Vec<Pipework> = Vec::new();

        if let Some(primary_pipework_lst) = primary_pipework_lst {
            for pipework_data in primary_pipework_lst {
                let new_pipework: Pipework = pipework_data
                    .to_owned()
                    .try_into()
                    .map_err(anyhow::Error::msg)?;

                pipework_lst.push(new_pipework);
            }
        };

        let pipework = PrimaryPipeworkLossesMixin::new(
            pipework_lst,
            // TODO review 1.0.0a9
            Arc::new(move |simtime| external_conditions.air_temp(simtime)),
            temp_internal_air_fn.clone(),
            simulation_time_iteration,
        );

        // With pre-heated storage tanks, there could be the situation of tanks without heat sources
        // They could just get warmed up with WWHRS water.
        let mut heat_source_data = heat_sources.clone();

        if !heat_sources.is_empty() {
            // sort heat source data in order from the bottom of the tank based on heater position
            heat_source_data = heat_source_data
                .iter()
                .sorted_by(|a, b| {
                    OrderedFloat(a.1.heater_position).cmp(&OrderedFloat(b.1.heater_position))
                })
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
            simulation_timestep: simulation_time_iteration.timestep,
            number_of_volumes,
            temp_flow_prev: Default::default(),
            temp_internal_air_fn,
            volume_total_in_litres,
            vol_n,
            cp,
            rho,
            temp_n,
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
            detailed_results: detailed_output.then_some(Default::default()),
            pipework,
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
        usage_events: Option<Vec<WaterEventResult>>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let mut q_use_w = 0.;
        let q_unmet_w = 0.;
        let mut volume_demanded = 0.;

        let mut temp_s3_n = self.temp_n.read().clone();
        let temp_ini_n = temp_s3_n.clone();

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
                self.extract_hot_water(*event, simtime)?;

            // Determine the new temperature distribution after displacement
            // Now that pre-heated sources can be the 'cold' feed, rearrangement of temperaturs, that used to
            // only happen before after the input from heat sources, could be required after the displacement
            // of water bringing new water from the 'cold' feed that could be warmer than the existing one.
            // flag is calculated for that purpose.
            let (temp_s3_n_new, rearrange) =
                self.calc_temps_after_extraction(remaining_vols, simtime)?;
            temp_s3_n = temp_s3_n_new;

            if rearrange {
                // Re-arrange the temperatures in the storage after energy input from pre-heated tank
                temp_s3_n = self.rearrange_temperatures(&temp_s3_n).1
            }

            *self.temp_n.write() = temp_s3_n.clone();

            volume_demanded += volume_used;
            q_use_w += energy_withdrawn;
        }

        self.temp_average_drawoff.store(
            match self.total_volume_drawoff.load(Ordering::SeqCst) {
                value if !relative_eq!(value, 0., epsilon = 1e-10, max_relative = 1e-9) => {
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

        // With the possibility of not having heat sources now, some parameters might not be defined now
        // in the for loop before and wouldn't be available for the testoutput unless initialised here.
        let mut q_x_in_n = vec![0.; self.number_of_volumes];
        let mut q_s6 = 0.;
        let mut q_in_h_w = 0.;
        let mut temp_s6_n = temp_s3_n.clone();
        let mut temp_s7_n = temp_s3_n.clone();
        let mut temp_s8_n = temp_s3_n.clone();
        let mut q_ls_this_heat_source = 0.;

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

            let calc = self.run_heat_sources(
                temp_after_prev_heat_source.clone(),
                &positioned_heat_source.heat_source.lock(),
                &heat_source_name,
                heater_layer,
                thermostat_layer,
                &self.q_ls_n_prev_heat_source.read().clone(),
                simtime,
            )?;
            let _ = std::mem::replace(&mut temp_s8_n, calc.temp_s8_n);
            let _ = std::mem::replace(&mut q_x_in_n, calc.q_x_in_n);
            q_s6 = calc.q_s6;
            let _ = std::mem::replace(&mut temp_s6_n, calc.temp_s6_n);
            let _ = std::mem::replace(&mut temp_s7_n, calc.temp_s7_n);
            q_in_h_w = calc.q_in_h_w;
            q_ls_this_heat_source = calc.q_ls;
            let q_ls_n_this_heat_source = calc.q_ls_n;

            temp_after_prev_heat_source = temp_s8_n.clone();
            q_ls += q_ls_this_heat_source;

            {
                let mut q_ls_n_prev = self.q_ls_n_prev_heat_source.write();
                for (i, q_ls_n) in q_ls_n_this_heat_source.iter().enumerate() {
                    q_ls_n_prev[i] += q_ls_n;
                }
            }

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

        self.testoutput(
            usage_events.as_ref().unwrap_or(&vec![]),
            volume_demanded,
            q_use_w,
            q_unmet_w,
            &temp_ini_n,
            &temp_s3_n,
            &q_x_in_n,
            q_s6,
            &temp_s6_n,
            &temp_s7_n,
            q_in_h_w,
            q_ls_this_heat_source,
            &temp_s8_n,
            self.temp_average_drawoff.load(Ordering::SeqCst),
            simtime,
        )?;

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
        event: WaterEventResult,
        simulation_time: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64, Vec<f64>)> {
        // Make a copy of the volume list to keep track of remaining volumes
        // Remaining volume of water in storage tank layers
        let mut remaining_vols = self.vol_n.clone();

        // Extract the temperature and required hot volume from the event
        let hot_volume = event.volume_hot;

        // # Remaining volume of hot water to be satisfied for current event
        let mut remaining_demanded_volume = hot_volume;
        let mut energy_withdrawn = 0.;

        let cold_water_source = &self.cold_feed.ultimate_cold_water_source();

        let mut temp_average_drawoff_volweighted: f64 =
            self.temp_average_drawoff_volweighted.load(Ordering::SeqCst);
        let mut total_volume_drawoff: f64 = self.total_volume_drawoff.load(Ordering::SeqCst);

        //  Loop through storage layers (starting from the top)
        for (layer_index, &layer_temp) in self.temp_n.read().iter().enumerate().rev() {
            let layer_vol = remaining_vols[layer_index];

            if remaining_demanded_volume < 0.
                || relative_eq!(
                    remaining_demanded_volume,
                    0.,
                    max_relative = 1e-09,
                    epsilon = 1e-10
                )
            {
                break;
            }

            // Skip this layer if its remaining volume is already zero
            if remaining_vols[layer_index] < 0.
                || relative_eq!(
                    remaining_vols[layer_index],
                    0.,
                    max_relative = 1e-09,
                    epsilon = 1e-10
                )
            {
                continue;
            }

            let required_vol: f64;
            // Volume of hot water required at this layer
            if layer_vol < remaining_demanded_volume
                || relative_eq!(layer_vol, remaining_demanded_volume, max_relative = 1e-09)
            {
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
            // Use ultimate cold water source temperature for consistent energy accounting
            let list_temp_vol =
                cold_water_source.get_temp_cold_water(hot_volume, simulation_time)?;
            let sum_t_by_v = FSum::with_all(list_temp_vol.iter().map(|(t, v)| t * v)).value();
            let sum_v = FSum::with_all(list_temp_vol.iter().map(|(_t, v)| v)).value();
            let temp_cold = sum_t_by_v / sum_v;

            energy_withdrawn +=
                // Calculation with event water parameters
                // self.__rho * self.__Cp * warm_vol_removed * (warm_temp - self.__cold_feed.temperature())
                // Calculation with layer water parameters
                self.rho
                    * self.cp
                    * required_vol
                    * (layer_temp - temp_cold);
        }

        self.temp_average_drawoff_volweighted
            .store(temp_average_drawoff_volweighted, Ordering::SeqCst);
        self.total_volume_drawoff
            .store(total_volume_drawoff, Ordering::SeqCst);

        // Handle case where demand exceeds tank capacity
        // Draw remaining volume from cold feed (which may be a pre-heat tank)
        if remaining_demanded_volume > 0. {
            //  Get water from cold feed for the remaining demand
            //  This triggers draw-off from pre-heat tank if cold_feed is a StorageTank
            let list_temp_vol_drawn = self
                .cold_feed
                .draw_off_water(remaining_demanded_volume, simulation_time)?;

            // Get ultimate cold water temperature for energy calculation
            let list_temp_vol = cold_water_source
                .get_temp_cold_water(remaining_demanded_volume, simulation_time)?;
            let sum_t_by_v = FSum::with_all(list_temp_vol.iter().map(|(t, v)| t * v)).value();
            let sum_v = FSum::with_all(list_temp_vol.iter().map(|(_t, v)| v)).value();
            let temp_cold = sum_t_by_v / sum_v;

            // Calculate volume-weighted temperature contribution from cold feed
            for (temp_feed, vol_feed) in list_temp_vol_drawn {
                self.temp_average_drawoff_volweighted
                    .fetch_add(vol_feed * temp_feed, Ordering::SeqCst);
                self.total_volume_drawoff
                    .fetch_add(vol_feed, Ordering::SeqCst);

                // Energy content is the difference between the drawn water temperature
                // and the ultimate cold water temperature reference
                energy_withdrawn += self.rho * self.cp * vol_feed * (temp_feed - temp_cold);
            }
        }

        //  Calculate the remaining total volume
        let remaining_total_volume: f64 = FSum::with_all(&remaining_vols).value();

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
    ) -> anyhow::Result<(Vec<f64>, bool)> {
        let mut new_temps = self.temp_n.read().clone();

        // If the 'cold' feed water is hotter than the existing water in the tank, rearrange will be needed.
        // as if it was a heat source coming from the cold feed.
        let mut flag_rearrange_layers = false;

        // Iterate from the top layer downwards
        for i in (0..self.vol_n.len()).rev() {
            // Determine how much volume needs to be added to this layer
            let mut needed_volume = self.vol_n[i] - remaining_vols[i];
            // If this layer is already full, continue to the next
            if needed_volume < 0.
                || relative_eq!(needed_volume, 0., max_relative = 1e-09, epsilon = 1e-10)
            {
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
                    if needed_volume < 0.
                        || relative_eq!(needed_volume, 0., max_relative = 1e-09, epsilon = 1e-10)
                    {
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
                    .draw_off_water(needed_volume, simulation_time)?;
                let sum_t_by_v = FSum::with_all(list_temp_vol.iter().map(|(t, v)| t * v)).value();
                let sum_v = FSum::with_all(list_temp_vol.iter().map(|(_t, v)| v)).value();

                let temp_cold_feed = sum_t_by_v / sum_v;
                volume_weighted_temperature += needed_volume * temp_cold_feed;
                if temp_cold_feed > temp_layer_min {
                    flag_rearrange_layers = true;
                }
            }

            new_temps[i] = volume_weighted_temperature / total_volume;
            remaining_vols[i] = total_volume;
        }
        Ok((new_temps, flag_rearrange_layers))
    }

    /// When the temperature of the volume i is higher than the one of the upper volume,
    /// then the 2 volumes are melded. This iterative process is maintained until the temperature
    /// of the volume i is lower or equal to the temperature of the volume i+1.
    fn rearrange_temperatures(&self, temp_s6_n: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let mut temp_s7_n = temp_s6_n.to_vec();

        loop {
            // Flag for which layers need mixing
            let mut mix_layer_n: Vec<u8> = vec![0; self.number_of_volumes];

            // for loop :-1 is important here!
            // loop through layers from bottom to top, without including top layer.
            // this is because the top layer has no upper layer to compare to
            for i in 0..self.vol_n.len() - 1 {
                if temp_s7_n[i] > temp_s7_n[i + 1]
                    || relative_eq!(temp_s7_n[i], temp_s7_n[i + 1], max_relative = 1e-09)
                {
                    // set layers to mix
                    mix_layer_n[i] = 1;
                    mix_layer_n[i + 1] = 1;
                    // mix temperatures of all applicable layers
                    // note error in formula 12 in standard as adding temperature to volume
                    // this is what I think they intended from the description (comment sic from Python code)
                    let temp_mix = FSum::with_all(
                        (0..self.vol_n.len())
                            .map(|k| self.vol_n[k] * temp_s7_n[k] * mix_layer_n[k] as f64),
                    )
                    .value()
                        / FSum::with_all(
                            (0..self.vol_n.len()).map(|l| self.vol_n[l] * mix_layer_n[l] as f64),
                        )
                        .value();
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
            heat_source_name.into(),
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
                    .temp_flow(heat_source, simulation_time)?
                    .unwrap_or(default_temp_flow);
                if self.heating_active[heat_source_name].load(Ordering::SeqCst) {
                    // upstream Python uses duck-typing/ polymorphism here, but we need to be more explicit

                    // The charging deadband (_heating_active) is the authority on
                    // whether to charge, so bypass the heat source's own minimum-setpoint
                    // schedule gate: once charging is active it must continue to the
                    // maximum setpoint even after the minimum-setpoint schedule
                    // deactivates mid-cycle.
                    let mut energy_potential = match heat_source {
                        HeatSource::Storage(HeatSourceWithStorageTank::Immersion(
                            immersion_heater,
                        )) => immersion_heater
                            .lock()
                            .energy_output_max(simulation_time, false),
                        HeatSource::Storage(HeatSourceWithStorageTank::Solar(_)) => unreachable!(), // this case was already covered in the first arm of this if let clause, so can't repeat here
                        HeatSource::Wet(heat_source_wet) => {
                            // TODO (from Python) Use different temperatures for flow and return in the call to
                            // heat_source.energy_output_max below
                            // Fallback to current tank temperature at heater layer when heat source has no setpoint
                            heat_source_wet.energy_output_max(
                                temp_flow,
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
                            self.pipework.calculate_primary_pipework_losses(
                                energy_potential,
                                temp_flow.into(),
                                Some(false),
                                &simulation_time,
                            )?;
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
        heat_source_name: ArcStr,
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
            temp_s3_n,
            &temp_s7_n,
            &q_x_in_n,
            &q_h_sto_s7,
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

        let _heat_source_output = self.heat_source_output(
            heat_source,
            heat_source_name,
            input_energy_adj,
            heater_layer,
            simtime,
            None,
            Some(control_max_diverter.is_some()),
        )?;
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
            * FSum::with_all((0..self.vol_n.len()).map(|i| self.vol_n[i] * temp_s6_n[i])).value();

        (q_s6, temp_s6_n)
    }

    /// Apply standby thermal losses to each layer and clamp to the setpoint.
    ///
    /// Both the loss-driving temperature and the final-temperature clamp are limited to
    /// the setpoint only for layers this heat source actually heated. A layer left above
    /// the setpoint by an earlier source - one with a higher setpoint, or the same source
    /// at a higher setpoint in the previous timestep - is not held there by this source,
    /// so it loses heat, and settles, at its actual temperature.
    ///
    /// Args:
    /// temp_s3_n: Layer temperatures before this source's energy input (°C), used to
    /// detect which layers this source heated this timestep.
    /// temp_s7_n: Layer temperatures after energy input and rearrangement (°C).
    /// q_x_in_n: Energy input to each layer from this source (kWh).
    /// q_h_sto_s7: Stored energy per layer after rearrangement (kWh).
    /// heater_layer: Index of the layer the heat source feeds.
    /// q_ls_n_prev_heat_source: Losses already attributed to earlier sources this
    /// timestep (kWh), subtracted to avoid double-counting.
    /// temp_setpntmax: Maximum setpoint of this source (°C), or None when uncontrolled.
    fn calc_temps_after_thermal_losses(
        &self,
        temp_s3_n: &[f64],
        temp_s7_n: &[f64],
        q_x_in_n: &[f64],
        q_h_sto_s7: &[f64],
        heater_layer: usize,
        q_ls_n_prev_heat_source: &[f64],
        temp_setpntmax: Option<f64>,
    ) -> (f64, f64, Vec<f64>, Vec<f64>) {
        let q_x_in_adj: f64 = FSum::with_all(q_x_in_n).value();

        // A layer is held at the setpoint by this source only if the source raised its
        // temperature this timestep. Energy enters at the heater layer and spreads upward by
        // buoyancy, so the per-layer temperature rise - not the heat-source input, which is
        // non-zero only at the heater layer - identifies which layers this source heated.
        let mut layer_warmed_n = Vec::with_capacity(self.number_of_volumes);

        for i in 0..self.number_of_volumes {
            layer_warmed_n.push(
                temp_s7_n[i] > temp_s3_n[i]
                    && !relative_eq!(
                        temp_s7_n[i],
                        temp_s3_n[i],
                        epsilon = 1e-10,
                        max_relative = 1e-9
                    ),
            );
        }

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
            // Cap the loss-driving temperature at the setpoint only for layers this source
            // warmed this timestep, and therefore holds at the setpoint. An unwarmed layer
            // above the setpoint was heated earlier and loses heat at its actual temperature.
            let temp_before_losses = match temp_setpntmax {
                Some(temp_setpntmax) if layer_warmed_n[i] => min_of_2(temp_s7_n[i], temp_setpntmax),
                _ => temp_s7_n[i],
            };

            let q_ls_n_step = (h_sto_ls * self.rho * self.cp)
                * (self.vol_n[i] / self.volume_total_in_litres)
                * (temp_before_losses - self.ambient_temperature)
                * self.simulation_timestep;

            let q_ls_n_step = max_of_2(0., q_ls_n_step - q_ls_n_prev_heat_source[i]);

            q_ls_n.push(q_ls_n_step);
        }

        // total thermal losses kWh
        let q_ls = FSum::with_all(&q_ls_n).value();

        self.storage_losses_kwh.store(q_ls, Ordering::SeqCst);

        // the final value of the temperature is reduced due to the effect of the thermal losses.
        // check temperature compared to set point
        // the temperature for each volume are limited to the set point for any volume controlled
        for i in 0..self.vol_n.len() {
            // Clamp to the setpoint only for layers this source warmed this timestep, and
            // therefore holds at the setpoint. Clamping an unwarmed layer would wrongly pull
            // it down to temp_setpntmax when it already exceeded the setpoint without any
            // contribution from this source - because the setpoint is lower this timestep than
            // the previous one, or because another source with a higher setpoint heated it.
            let temp_s8_n_step = match temp_setpntmax {
                Some(temp_setpntmax)
                    if layer_warmed_n[i]
                        && temp_s7_n[i] > temp_setpntmax
                        && !relative_eq!(
                            temp_s7_n[i],
                            temp_setpntmax,
                            epsilon = 1e-10,
                            max_relative = 1e-9
                        ) =>
                {
                    // Case 2 - Temperature exceeding the set point. Compared with a tolerance
                    // (as for layer_warmed_n) because the source drives the layer temperature
                    // onto the setpoint, so a last-bit difference would otherwise flip the clamp
                    // between platforms.
                    temp_setpntmax
                }
                _ =>
                // Case 1 - Temperature below the set point
                // TODO (from Python) - spreadsheet accounts for total thermal losses not just layer

                // the final value of the temperature
                // is reduced due to the effect of the thermal losses
                // Formula (14) in the standard appears to have error as addition not multiply
                // and P instead of rho
                {
                    temp_s7_n[i] - (q_ls_n[i] / (self.rho * self.cp * self.vol_n[i]))
                }
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
            if let Some(temp_setpntmax) = temp_setpntmax {
                if temp_s7_n[heater_layer] > temp_setpntmax
                    && !relative_eq!(
                        temp_s7_n[heater_layer],
                        temp_setpntmax,
                        epsilon = 1e-10,
                        max_relative = 1e-9
                    )
                {
                    for i in heater_layer..self.number_of_volumes {
                        energy_surplus += q_h_sto_s7[i]
                            - q_ls_n[i]
                            - (self.rho * self.cp * self.vol_n[i] * temp_setpntmax);
                    }
                }
            }
            // the thermal energy provided to the system (from heat sources) shall be limited
            // adjustment of the energy delivered to the storage according with the set temperature
            // potential input from generation
            // TODO (from Python code) - find in standard - availability of back-up - where is this from?
            // also referred to as electrical power on
            let sto_bu_on = 1.;
            // BS EN 15316-5 Formula (16) can yield a negative result in multi-source
            // tanks when a prior heat source has already driven the stored energy above
            // the setpoint. The standard does not directly address the possible negative
            // result, but describes this equation as limiting "thermal energy
            // provided to the system", which cannot be negative, so clamp to zero here.
            max_of_2(
                0.,
                min_of_2(q_x_in_adj - energy_surplus, q_x_in_adj * sto_bu_on),
            )
        } else {
            0.
        };

        (q_in_h_w, q_ls, temp_s8_n, q_ls_n)
    }

    /// Calculates pipework loss before sending on the demand energy
    ///
    /// Args:
    /// heat_source: The heat source to demand energy from.
    /// input_energy_adj: Adjusted energy input in kWh.
    /// heater_layer: Index of the heater layer in the tank.
    /// ignore_standard_ctrl: If True, bypass the standard time control check
    /// when demanding energy from an ImmersionHeater.
    /// Set when PV diverter is active.
    fn heat_source_output(
        &self,
        heat_source: &HeatSource,
        heat_source_name: ArcStr,
        input_energy_adj: f64,
        _heater_layer: usize,
        simulation_time_iteration: SimulationTimeIteration,
        smart_hot_water_tank: Option<&SmartHotWaterTank>, // the temp_flow method might need to be called as a smart hot water tank if this is a storage tank composed by a smart hot water tank
        ignore_standard_ctrl: Option<bool>,
    ) -> anyhow::Result<f64> {
        let ignore_standard_ctrl = ignore_standard_ctrl.unwrap_or(false);
        // if immersion heater, no pipework losses
        // TODO (from Python):  Critical - temp_flow cannot be None for downstream method calculate_primary_pipework_losses
        // but providing a fallback value will change the e2e test results
        let temp_flow = match smart_hot_water_tank {
            None => self.temp_flow(heat_source, simulation_time_iteration)?,
            Some(smart_hot_water_tank) => {
                smart_hot_water_tank.temp_flow(simulation_time_iteration)?
            }
        };

        // Input energy clamped to zero if within 1e-10 of zero
        // so that almost zero negative numbers (caused by
        // floating point error) do not cause errors in subsequent code
        let input_energy_adj =
            if relative_eq!(input_energy_adj, 0.0, epsilon = 1e-10, max_relative = 1e-9) {
                0.
            } else {
                input_energy_adj
            };

        // The charging deadband (_heating_active) - or an active PV diverter - is
        // the authority on whether to charge, so bypass the heat source's own
        // minimum-setpoint schedule gate: once charging is active it must continue
        // to the maximum setpoint even after the minimum-setpoint schedule
        // deactivates mid-cycle. A heat source absent from _heating_active is not
        // deadband-active, so it gets no bypass.
        let bypass_min_ctrl = ignore_standard_ctrl
            || self
                .heating_active
                .get(&heat_source_name)
                .map(|bool| bool.load(Ordering::SeqCst))
                .unwrap_or_default();

        match heat_source {
            HeatSource::Storage(HeatSourceWithStorageTank::Immersion(immersion)) => {
                immersion.lock().demand_energy(
                    input_energy_adj,
                    Some(bypass_min_ctrl),
                    simulation_time_iteration,
                )
            }
            HeatSource::Storage(HeatSourceWithStorageTank::Solar(solar)) => Ok(solar
                .lock()
                .demand_energy(input_energy_adj, simulation_time_iteration.index)),
            HeatSource::Wet(ref wet_heat_source) => {
                let (primary_pipework_losses_kwh, primary_gains) =
                    self.pipework.calculate_primary_pipework_losses(
                        input_energy_adj,
                        temp_flow,
                        None,
                        &simulation_time_iteration,
                    )?;
                // Save for reporting
                self.primary_pipework_losses_kwh
                    .store(primary_pipework_losses_kwh, Ordering::SeqCst);
                let input_energy_adj = input_energy_adj + primary_pipework_losses_kwh;

                // TODO Use different temperatures for flow and return in the call to
                // heat_source.demand_energy below
                let heat_source_output = wet_heat_source.demand_energy(
                    input_energy_adj,
                    temp_flow,
                    temp_flow,
                    Some(bypass_min_ctrl),
                    simulation_time_iteration,
                )? - primary_pipework_losses_kwh;
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
            (Some(setpointmax), Some(setpointmin)) if setpointmin > setpointmax => {
                bail!("setpntmin: {setpointmin} must not be greater than setpntmax: {setpointmax}");
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
        let (setpntmin, setpntmax) =
            self.retrieve_setpnt(heat_source, simulation_time_iteration)?;

        // In an off period (no setpoint) deactivate a source left active by a
        // previous timestep, before the charging block, so _heating_active stays
        // consistent with the schedule and a stale active flag cannot drive
        // charging across the on-to-off transition. The temperature cut-out is
        // in _determine_heat_source_switch_off, evaluated on the post-charge
        // temperatures, to preserve the min/max charging hysteresis.
        if setpntmax.is_none() {
            self.heating_active[heat_source_name].store(false, Ordering::SeqCst);
        } else if setpntmin.is_some_and(|setpntmin| {
            temp_s3_n[thermostat_layer] < setpntmin
                || relative_eq!(temp_s3_n[thermostat_layer], setpntmin, max_relative = 1e-09)
        }) {
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

        if setpntmax.is_none_or(|setpntmax| {
            temp_s8_n[thermostat_layer] > setpntmax
                || relative_eq!(temp_s8_n[thermostat_layer], setpntmax, max_relative = 1e-09)
        }) {
            self.heating_active[heat_source_name].store(false, Ordering::SeqCst);
        };

        Ok(())
    }

    fn temp_flow(
        &self,
        heat_source: &HeatSource,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<Option<f64>> {
        let (_, setpntmax) = self.retrieve_setpnt(heat_source, simulation_time_iteration)?;

        let setpntmax = setpntmax
            .inspect(|&setpntmax| {
                *self.temp_flow_prev.write() = Some(setpntmax);
            })
            .or_else(|| *self.temp_flow_prev.read());

        Ok(setpntmax)
    }

    /// Return the ambient temperature surrounding a primary pipework segment.
    ///
    /// Delegates to PrimaryPipeworkLossesMixin temp_surrounding_pipework.
    ///
    /// Args:
    ///     pipework_data: Pipework object to query location from.
    /// Returns:
    ///     Surrounding temperature in °C.
    #[cfg(test)]
    fn temperature_surrounding_primary_pipework(
        &self,
        pipework_data: &Pipework,
        simtime: &SimulationTimeIteration,
    ) -> f64 {
        PrimaryPipeworkLossesMixin::temp_surrounding_pipework(
            pipework_data,
            self.pipework.temp_external_air_fn().clone(),
            self.pipework.temp_internal_air_fn().clone(),
            simtime,
        )
    }

    pub(crate) fn get_cold_water_source(&self) -> &WaterSupply {
        &self.cold_feed
    }

    pub(crate) fn get_temp_hot_water(
        &self,
        volume_req: f64,
        volume_req_already: Option<f64>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        let volume_req_already = volume_req_already.unwrap_or(0.);
        let mut volume_req_cumulative = volume_req + volume_req_already;

        let mut list_temp_vol: Vec<(f64, f64)> = vec![];
        // Loop through storage layers (starting from the top)
        for (layer_index, &layer_temp) in self.temp_n.read().iter().enumerate().rev() {
            let layer_vol = self.vol_n[layer_index];
            let volume_from_current_layer = volume_req_cumulative.min(layer_vol);

            list_temp_vol.push((layer_temp, volume_from_current_layer));
            volume_req_cumulative -= volume_from_current_layer;

            if volume_req_cumulative < 0.
                || relative_eq!(
                    volume_req_cumulative,
                    0.,
                    max_relative = 1e-09,
                    epsilon = 1e-10
                )
            {
                break;
            }
        }

        // If requested volume exceeds tank capacity, get remaining from cold feed
        if volume_req_cumulative > 0. {
            let cold_feed_temp_vol = self
                .cold_feed
                .get_temp_cold_water(volume_req_cumulative, simulation_time_iteration)?;
            list_temp_vol.extend(cold_feed_temp_vol)
        }

        // Base temperature on the part of the draw-off for volume_req, and
        // ignore any volume previously considered
        let mut list_temp_vol_req: Vec<(f64, f64)> = vec![];
        let mut volume_still_to_satisfy = volume_req;
        for (layer_temp, layer_vol) in list_temp_vol.iter().rev() {
            let volume_from_current_layer = volume_still_to_satisfy.min(*layer_vol);
            list_temp_vol_req.push((*layer_temp, volume_from_current_layer));
            volume_still_to_satisfy -= volume_from_current_layer;

            if volume_still_to_satisfy < 0.
                || relative_eq!(
                    volume_still_to_satisfy,
                    0.,
                    max_relative = 1e-09,
                    epsilon = 1e-10
                )
            {
                break;
            }
        }

        Ok(list_temp_vol_req.into_iter().rev().collect_vec())
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

    fn testoutput(
        &self,
        usage_events: &[WaterEventResult],
        volume_extracted: f64,
        q_use_w: f64,
        q_unmet_w: f64,
        temp_ini_n: &[f64],
        temp_s3_n: &[f64],
        q_x_in_n: &[f64],
        q_s6: f64,
        temp_s6_n: &[f64],
        temp_s7_n: &[f64],
        q_in_h_w: f64,
        q_ls: f64,
        temp_s8_n: &[f64],
        temp_average: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<()> {
        let mut detailed_output = match self.detailed_results.as_ref() {
            None => return Ok(()),
            Some(detailed_output) => detailed_output.write(),
        };

        let demand = summarise_events(usage_events);
        if simtime.index == 0 {
            fn header_dup(header: &str, n: usize) -> Vec<StringOrNumber> {
                (0..n)
                    .map(|i| format!("{header} {}", i + 1).into())
                    .collect()
            }

            let mut header_row: Vec<StringOrNumber> = [
                "time",
                "volume total",
                "specific heat",
                "density",
                "cold water",
                "events",
                "volume extracted",
            ]
            .into_iter()
            .map(Into::into)
            .collect();
            header_row.extend(header_dup("initial temp.", temp_ini_n.len()));
            header_row.extend(["energy withdrawn", "energy unmet"].map(Into::into));
            header_row.extend(header_dup("temp. after volume withdrawn", temp_s3_n.len()));
            header_row.extend(header_dup("potential energy input", q_x_in_n.len()));
            header_row.push("theoretical energy stored after energy input".into());
            header_row.extend(header_dup(
                "theoretical temp. after energy input",
                temp_s6_n.len(),
            ));
            header_row.extend(header_dup("temp. after volume mixing", temp_s7_n.len()));
            header_row.extend(
                ["energy input (adjusted)", "thermal losses"]
                    .into_iter()
                    .map(Into::into),
            );
            header_row.extend(header_dup("temp. after thermal losses", temp_s8_n.len()));
            header_row.push("temp_average_drawoff".into());

            detailed_output.push(header_row);

            let mut units_row: Vec<StringOrNumber> = [
                "h",
                "litres",
                "kWh/kgK",
                "kg/l",
                "oC",
                "Type: litres hot (litres @ oC)",
                "litres",
            ]
            .into_iter()
            .map(Into::into)
            .collect();
            units_row.extend(iter::repeat_n("oC".into(), temp_ini_n.len()));
            units_row.extend(["kWh", "kWh"].into_iter().map(Into::into));
            units_row.extend(iter::repeat_n("oC".into(), temp_s3_n.len()));
            units_row.extend(iter::repeat_n("kWh".into(), q_x_in_n.len()));
            units_row.push("kWh".into());
            units_row.extend(iter::repeat_n("oC".into(), temp_s6_n.len()));
            units_row.extend(iter::repeat_n("oC".into(), temp_s7_n.len()));
            units_row.extend(["kWh", "kWh"].into_iter().map(Into::into));
            units_row.extend(iter::repeat_n("oC".into(), temp_s8_n.len()));
            units_row.push("oC".into());

            detailed_output.push(units_row);
        }

        let temp_cold_water: StringOrNumber =
            if relative_eq!(volume_extracted, 0.0, epsilon = 1e-10, max_relative = 1e-9) {
                "".into() // using empty string to represent None
            } else {
                let list_temp_vol = self
                    .cold_feed
                    .get_temp_cold_water(volume_extracted, simtime)?;
                (FSum::with_all(list_temp_vol.iter().map(|(t, v)| t * v)).value()
                    / FSum::with_all(list_temp_vol.iter().map(|(_, v)| v)).value())
                .into()
            };

        let mut values_row: Vec<StringOrNumber> = vec![
            simtime.hour_of_day().into(),
            self.volume_total_in_litres.into(),
            self.cp.into(),
            self.rho.into(),
            temp_cold_water,
            demand.into(),
            volume_extracted.into(),
        ];
        values_row.extend(temp_ini_n.iter().map(|t| t.into()));
        values_row.extend([q_use_w.into(), q_unmet_w.into()]);
        values_row.extend(temp_s3_n.iter().map(|t| t.into()));
        values_row.extend(q_x_in_n.iter().map(|q| q.into()));
        values_row.push(q_s6.into());
        values_row.extend(temp_s6_n.iter().map(|t| t.into()));
        values_row.extend(temp_s7_n.iter().map(|t| t.into()));
        values_row.extend([q_in_h_w.into(), q_ls.into()]);
        values_row.extend(temp_s8_n.iter().map(|t| t.into()));
        values_row.push(temp_average.into());

        detailed_output.push(values_row);

        Ok(())
    }

    /// draw off hot water layers until required volume is provided.
    ///
    /// Arguments:
    /// * volume    -- volume of water required
    pub(crate) fn draw_off_hot_water(
        &self,
        volume: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<(Option<f64>, f64)> {
        if relative_eq!(volume, 0., epsilon = 1e-10, max_relative = 1e-9) {
            return Ok((None, volume));
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
            .get_temp_cold_water(volume, simulation_time_iteration)?;
        let sum_t_by_v = FSum::with_all(list_temp_vol.iter().map(|(t, v)| t * v)).value();
        let sum_v = FSum::with_all(list_temp_vol.iter().map(|(_t, v)| v)).value();

        self.temp_average_drawoff
            .store(sum_t_by_v / sum_v, Ordering::SeqCst);

        let _temp_ini_n = self.temp_n.clone();
        let _temp_s3_n = self.temp_n.clone();

        // Loop through storage layers (starting from the top)
        for (layer_index, &layer_temp) in self.temp_n.read().iter().enumerate().rev() {
            let layer_vol = remaining_vols[layer_index];

            // This cannot happen in the preheated tank. Check!
            // Skip this layer if its remaining volume is already zero
            // if remaining_vols[layer_index] <= 0.0:
            //     continue

            // Volume of water required at this layer
            let required_vol;
            if layer_vol < remaining_demanded_volume
                || relative_eq!(
                    layer_vol,
                    remaining_demanded_volume,
                    max_relative = 1e-09,
                    epsilon = 1e-10
                )
            {
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
                .get_temp_cold_water(required_vol, simulation_time_iteration)?;
            let sum_t_by_v = FSum::with_all(list_temp_vol.iter().map(|(t, v)| t * v)).value();
            let sum_v = FSum::with_all(list_temp_vol.iter().map(|(_t, v)| v)).value();
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

            if remaining_demanded_volume < 0.
                || relative_eq!(
                    remaining_demanded_volume,
                    0.,
                    max_relative = 1e-09,
                    epsilon = 1e-10
                )
            {
                break;
            }
        }

        if remaining_demanded_volume > 0.0 {
            let list_temp_vol = self
                .cold_feed
                .draw_off_water(remaining_demanded_volume, simulation_time_iteration)?;
            let sum_t_by_v = FSum::with_all(list_temp_vol.iter().map(|(t, v)| t * v)).value();
            let sum_v = FSum::with_all(list_temp_vol.iter().map(|(_t, v)| v)).value();
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
            self.calc_temps_after_extraction(remaining_vols, simulation_time_iteration)?;

        if flag_rearrange_layers {
            // Re-arrange the temperatures in the storage after energy input from pre-heated tank
            (_, new_temp_distribution) = self.rearrange_temperatures(&new_temp_distribution);
        }

        *self.temp_n.write() = new_temp_distribution;

        // Return the average temperature and volume drawn
        Ok((
            Some(self.temp_average_drawoff.load(Ordering::SeqCst)),
            self.total_volume_drawoff.load(Ordering::SeqCst),
        ))
    }

    fn additional_energy_input(
        &self,
        heat_source: &HeatSource,
        heat_source_name: &str,
        energy_input: f64,
        control_max_diverter: Option<&Control>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        if relative_eq!(energy_input, 0., epsilon = 1e-10, max_relative = 1e-9) {
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
            heat_source_name.into(),
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

    // TODO Python has get_temp_cold_water and draw_off_water defined here
    // which are called but currently unsure where from

    /// Return the pre-heated water temperature for the current timestep and the volume drawn
    pub(crate) fn get_temp_cold_water(
        &self,
        volume_needed: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        // TODO this matches Python - is it correct?
        self.get_temp_hot_water(volume_needed, None, simulation_time_iteration)
    }

    /// Return the pre-heated water temperature for the current timestep and the volume drawn
    pub(crate) fn draw_off_water(
        &self,
        volume_needed: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        let list_temp_vol = self.get_temp_cold_water(volume_needed, simulation_time_iteration);
        self.draw_off_hot_water(volume_needed, simulation_time_iteration)?;
        list_temp_vol
    }

    pub(crate) fn output_results(&self) -> Option<Vec<Vec<StringOrNumber>>> {
        self.detailed_results
            .as_ref()
            .map(|results| results.read().clone())
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

#[derive(Debug, PartialEq)]
struct PartialTemperatureCalculation {
    temp_s8_n: Vec<f64>,
    q_s6: f64,
    temp_s6_n: Vec<f64>,
    temp_s7_n: Vec<f64>,
    q_in_h_w: f64,
    q_ls: f64,
    q_ls_n: Vec<f64>,
}

/// A struct to represent a smart hot water storage tank/cylinder
#[derive(Debug)]
pub struct SmartHotWaterTank {
    storage_tank: StorageTank,
    power_pump_kw: f64,
    max_flow_rate_pump_l_per_min: f64,
    temp_usable: f64,
    temp_setpnt_max: Control,
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
        temp_setpnt_max: Control,
        cold_feed: WaterSupply,
        simulation_time_iteration: SimulationTimeIteration,
        heat_sources: IndexMap<ArcStr, PositionedHeatSource>,
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
            &simulation_time_iteration,
            heat_sources,
            temp_internal_air_fn,
            external_conditions,
            detailed_output,
            number_of_volumes.into(),
            primary_pipework_lst,
            contents,
            None,
            None,
            None,
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
        usage_events: Option<Vec<WaterEventResult>>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // N.B. implementation from StorageTank but calling SmartHotWaterTank specific methods further down
        let mut q_use_w = 0.;
        let q_unmet_w = 0.;
        let mut volume_demanded = 0.;

        let mut temp_s3_n = self.storage_tank.temp_n.read().clone();
        let temp_ini_n = temp_s3_n.clone();

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
                self.storage_tank.extract_hot_water(*event, simtime)?;

            let (temp_s3_n_new, rearrange) = self
                .storage_tank
                .calc_temps_after_extraction(remaining_vols, simtime)?;
            temp_s3_n = temp_s3_n_new;

            if rearrange {
                // Re-arrange the temperatures in the storage after energy input from pre-heated tank
                temp_s3_n = self.storage_tank.rearrange_temperatures(&temp_s3_n).1
            }

            *self.storage_tank.temp_n.write() = temp_s3_n.clone();

            volume_demanded += volume_used;
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

        // With the possibility of not having heat sources now, some parameters might not be defined now
        // in the for loop before and wouldn't be available for the testoutput unless initialised here.
        let mut q_x_in_n = vec![0.; self.storage_tank.number_of_volumes];
        let mut q_s6 = 0.;
        let mut q_in_h_w = 0.;
        let mut temp_s6_n = temp_s3_n.clone();
        let mut temp_s7_n = temp_s3_n.clone();
        let mut temp_s8_n = temp_s3_n.clone();
        let mut q_ls_this_heat_source = 0.;

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

            let calc = self.run_heat_sources(
                temp_after_prev_heat_source.clone(),
                &positioned_heat_source.heat_source.lock(),
                &heat_source_name,
                heater_layer,
                thermostat_layer,
                &self.storage_tank.q_ls_n_prev_heat_source.read().clone(),
                simtime,
            )?;
            let _ = std::mem::replace(&mut temp_s8_n, calc.temp_s8_n);
            let _ = std::mem::replace(&mut q_x_in_n, calc.q_x_in_n);
            q_s6 = calc.q_s6;
            let _ = std::mem::replace(&mut temp_s6_n, calc.temp_s6_n);
            let _ = std::mem::replace(&mut temp_s7_n, calc.temp_s7_n);
            q_in_h_w = calc.q_in_h_w;
            q_ls_this_heat_source = calc.q_ls;
            let q_ls_n_this_heat_source = calc.q_ls_n;

            temp_after_prev_heat_source = temp_s8_n.clone();
            q_ls += q_ls_this_heat_source;

            {
                let mut q_ls_n_prev = self.storage_tank.q_ls_n_prev_heat_source.write();
                for (i, q_ls_n) in q_ls_n_this_heat_source.iter().enumerate() {
                    q_ls_n_prev[i] += q_ls_n;
                }
            }

            // Trigger heating to stop
            self.determine_heat_source_switch_off(
                &temp_s8_n,
                &heat_source_name,
                heater_layer,
                simtime,
            )?;
        }

        self.testoutput(
            usage_events.as_ref().unwrap_or(&vec![]),
            volume_demanded,
            q_use_w,
            q_unmet_w,
            &temp_ini_n,
            &temp_s3_n,
            &q_x_in_n,
            q_s6,
            &temp_s6_n,
            &temp_s7_n,
            q_in_h_w,
            q_ls_this_heat_source,
            &temp_s8_n,
            self.storage_tank
                .temp_average_drawoff
                .load(Ordering::SeqCst),
            simtime,
        )?;

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
            &temp_s3_n,
            heat_source,
            heat_source_name.into(),
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
                let temp_flow = self
                    .temp_flow(simulation_time)?
                    .unwrap_or(default_temp_flow);
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
                                temp_flow,
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
                            .storage_tank
                            .pipework
                            .calculate_primary_pipework_losses(
                                energy_potential,
                                Some(temp_flow),
                                Some(false),
                                &simulation_time,
                            )?;
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

    fn additional_energy_input(
        &self,
        heat_source: &HeatSource,
        heat_source_name: &str,
        energy_input: f64,
        control_max_diverter: Option<&Control>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // N.B. implementation from StorageTank but calling SmartHotWaterTank specific methods further down

        if relative_eq!(energy_input, 0., epsilon = 1e-10, max_relative = 1e-9) {
            return Ok(0.);
        }

        let heat_source_data = &self.storage_tank.heat_source_data[heat_source_name];

        let heater_layer = (heat_source_data.heater_position
            * self.storage_tank.number_of_volumes as f64) as usize;

        let mut q_x_in_n = vec![0.; self.storage_tank.number_of_volumes];
        q_x_in_n[heater_layer] = energy_input;

        // N.B. we're calling the SmartHotWaterTank specific method here
        let TemperatureCalculation {
            temp_s8_n,
            q_in_h_w,
            q_ls_n: q_ls_n_this_heat_source,
            ..
        } = self.calc_final_temps(
            &self.storage_tank.temp_n.read(),
            heat_source,
            heat_source_name.into(),
            q_x_in_n,
            heater_layer,
            &self.storage_tank.q_ls_n_prev_heat_source.read(),
            control_max_diverter,
            simulation_time_iteration,
        )?;

        for (i, q_ls_n) in q_ls_n_this_heat_source.iter().enumerate() {
            let mut q_ls_n_prev = self.storage_tank.q_ls_n_prev_heat_source.write();
            q_ls_n_prev[i] += *q_ls_n;
        }

        *self.storage_tank.temp_n.write() = temp_s8_n;

        Ok(q_in_h_w)
    }

    /// Return the DHW recoverable heat losses as internal gain for the current timestep in W
    pub(crate) fn internal_gains(&self) -> f64 {
        self.storage_tank.internal_gains()
    }

    pub(crate) fn output_results(&self) -> Option<Vec<Vec<StringOrNumber>>> {
        self.storage_tank.output_results()
    }

    pub(crate) fn get_cold_water_source(&self) -> &WaterSupply {
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
        if setpntmax.is_none_or(|setpntmax| {
            state_of_charge > setpntmax
                || relative_eq!(state_of_charge, setpntmax, max_relative = 1e-09)
        }) {
            self.storage_tank.heating_active[heat_source_name].store(false, Ordering::SeqCst);
        }

        Ok(())
    }

    // Making this method return a Result as the corresponding method on StorageTank does, and in the original Python
    // SmartHotWaterTank subclasses StorageTank. We're making the assumption here that .setpnt() will always return a
    // `Some` value in normal functioning. If that isn't the case, we would need to address this differently.
    fn temp_flow(&self, simtime: SimulationTimeIteration) -> anyhow::Result<Option<f64>> {
        Ok(self.temp_setpnt_max.setpnt(&simtime))
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
            .get_temp_cold_water(self.storage_tank.volume_total_in_litres, simtime)?;
        let sum_t_by_v = FSum::with_all(list_temp_vol.iter().map(|(t, v)| t * v)).value();
        let sum_v = FSum::with_all(list_temp_vol.iter().map(|(_t, v)| v)).value();
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
            if t_h_i > t_u || relative_eq!(t_h_i, t_u, max_relative = 1e-09) {
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

    /// Charge the storage and return the temperatures and energy for the timestep.
    ///
    /// Sizes the charge so the realised state of charge - after the top-up pump and
    /// thermal losses - meets the target, by solving for it directly, applies it
    /// through the full charging path, and records the pump energy consumed.
    ///
    /// Args:
    ///     temp_s3_n: storage layer temperatures after volume withdrawal.
    ///     heat_source: heat source charging the storage this timestep.
    ///     Q_x_in_n: potential energy input per layer, in kWh.
    ///     heater_layer: index of the layer the heat source charges.
    ///     Q_ls_n_prev_heat_source: thermal losses already attributed to earlier
    ///     heat sources this timestep, per layer.
    ///     controlmax_diverter: diverter control selecting the maximum state of
    ///     charge, or None to use the heat source's own maximum.
    ///
    /// Returns:
    ///     Final temperatures, potential energy input, theoretical stored energy
    ///     after input, temperatures after input, temperatures after pumping and
    ///     rearrangement, adjusted energy input, total thermal losses, and
    ///     per-layer thermal losses.
    fn calc_final_temps(
        &self,
        temp_s3_n: &[f64],
        heat_source: &HeatSource,
        heat_source_name: ArcStr,
        q_x_in_n: Vec<f64>,
        heater_layer: usize,
        q_ls_n_prev_heat_source: &[f64],
        control_max_diverter: Option<&Control>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<TemperatureCalculation> {
        let temp_setpntmax = self.temp_setpnt_max.setpnt(&simtime);
        let energy_available = FSum::with_all(q_x_in_n.iter().copied()).value();

        // Target state of charge for the heat source being considered
        let soc_max = if let Some(control_max_diverter) = control_max_diverter {
            control_max_diverter.setpnt(&simtime)
        } else {
            let (_, soc_max) = self.retrieve_setpnt(heat_source, simtime)?;

            soc_max
        };

        // Size the charge so the realised state of charge - after the top-up pump and
        // thermal losses - meets the target, by solving for it directly. With no
        // maximum set the source has no target to charge towards, so it delivers
        // nothing and the storage only coasts through pumping and losses.
        let energy_charge = if let Some(soc_max) = soc_max {
            self.solve_charge_energy_for_state_of_charge(
                temp_s3_n,
                soc_max,
                energy_available,
                heater_layer,
                q_ls_n_prev_heat_source,
                temp_setpntmax,
                simtime,
            )?
        } else {
            0.
        };

        let PartialTemperatureCalculation {
            temp_s8_n,
            q_s6,
            temp_s6_n,
            temp_s7_n,
            q_in_h_w,
            q_ls,
            q_ls_n,
        } = self.charge_to_temps(
            temp_s3_n,
            energy_charge,
            heater_layer,
            q_ls_n_prev_heat_source,
            temp_setpntmax,
            simtime,
        )?;

        // Adjust energy input based on actual usage
        let input_energy_adj = q_in_h_w;

        #[cfg(test)]
        self.storage_tank
            .energy_demand_test
            .store(input_energy_adj, Ordering::SeqCst);

        // Actual heat source output
        let heat_source_output = self.storage_tank.heat_source_output(
            heat_source,
            heat_source_name,
            input_energy_adj,
            heater_layer,
            simtime,
            Some(self),
            Some(control_max_diverter.is_some()),
        )?;

        // calculate volume pumped using actual heat source output
        let volumes = self.storage_tank.vol_n.clone();
        let volume_pumped = self.bottom_to_top_pump_volume(
            temp_s3_n,
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

    /// Charge the storage with a heater-layer energy through the full charging path.
    ///
    /// Runs energy input, buoyancy rearrangement, the top-up pump, a further
    /// rearrangement and the thermal-loss step - the sequence that fixes the
    /// realised storage temperatures for a given charge - so the same path can be
    /// evaluated both for the chosen charge and when solving for it.
    ///
    /// Args:
    ///     temp_s3_n: storage layer temperatures after volume withdrawal.
    ///     energy: charge delivered to the heater layer this timestep, in kWh.
    ///     heater_layer: index of the layer the heat source charges.
    ///     Q_ls_n_prev_heat_source: thermal losses already attributed to earlier
    ///     heat sources this timestep, per layer, to avoid double-counting.
    ///     temp_setpntmax: maximum storage temperature, or None when uncontrolled.
    ///
    /// Returns:
    ///     Final temperatures after thermal losses, theoretical stored energy after
    ///     input, temperatures after input, temperatures after pumping and
    ///     rearrangement, adjusted energy input, total thermal losses, and
    ///     per-layer thermal losses.
    fn charge_to_temps(
        &self,
        temp_s3_n: &[f64],
        energy: f64,
        heater_layer: usize,
        q_ls_n_prev_heat_source: &[f64],
        temp_setpntmax: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<PartialTemperatureCalculation> {
        let mut q_x_in_n = vec![0.; self.storage_tank.number_of_volumes];
        q_x_in_n[heater_layer] = energy;

        let (q_s6, temp_s6_n) = self
            .storage_tank
            .calc_temps_with_energy_input(temp_s3_n, &q_x_in_n);
        let (_, temp_s7_n) = self.storage_tank.rearrange_temperatures(&temp_s6_n);
        let temp_s7_n =
            self.calc_temps_after_top_up_pump(&temp_s7_n, energy, heater_layer, simtime)?;
        let (q_h_sto_s7, temp_s7_n) = self.storage_tank.rearrange_temperatures(&temp_s7_n);
        let (q_in_h_w, q_ls, temp_s8_n, q_ls_n) =
            self.storage_tank.calc_temps_after_thermal_losses(
                temp_s3_n,
                &temp_s7_n,
                &q_x_in_n,
                &q_h_sto_s7,
                heater_layer,
                q_ls_n_prev_heat_source,
                temp_setpntmax,
            );

        Ok(PartialTemperatureCalculation {
            temp_s8_n,
            q_s6,
            temp_s6_n,
            temp_s7_n,
            q_in_h_w,
            q_ls,
            q_ls_n,
        })
    }

    fn realised_state_of_charge(
        &self,
        energy: f64,
        temp_s3_n: &[f64],
        heater_layer: usize,
        q_ls_n_prev_heat_source: &[f64],
        temp_setpntmax: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let temp_s8_n = self
            .charge_to_temps(
                temp_s3_n,
                energy,
                heater_layer,
                q_ls_n_prev_heat_source,
                temp_setpntmax,
                simtime,
            )?
            .temp_s8_n;

        self.calc_state_of_charge(&temp_s8_n, simtime)
    }

    /// Solve for the smallest heater charge whose realised state of charge meets the target.
    ///
    /// The realised state of charge rises with the charge delivered, so the energy
    /// is found by bisection between no charge and all the energy available this
    /// timestep. Returns no charge when the target is already met without any
    /// input, or all the available energy when the target cannot be reached.
    ///
    /// Args:
    ///     temp_s3_n: storage layer temperatures after volume withdrawal.
    ///     soc_target: state of charge the charge should achieve.
    ///     energy_available: energy the heat source can deliver this timestep, in kWh.
    ///     heater_layer: index of the layer the heat source charges.
    ///     Q_ls_n_prev_heat_source: thermal losses already attributed to earlier
    ///     heat sources this timestep, per layer.
    ///     temp_setpntmax: maximum storage temperature, or None when uncontrolled.
    ///
    /// Returns:
    ///     Charge to deliver to the heater layer this timestep, in kWh.
    fn solve_charge_energy_for_state_of_charge(
        &self,
        temp_s3_n: &[f64],
        soc_target: f64,
        energy_available: f64,
        heater_layer: usize,
        q_ls_n_prev_heat_source: &[f64],
        temp_setpntmax: Option<f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let mut energy_lower = 0.;
        let mut energy_upper = max_of_2(energy_available, 0.);

        // No interior solution when the target is already met without input or
        // cannot be reached with all the available energy: return the bound.
        let soc_lower = self.realised_state_of_charge(
            energy_lower,
            temp_s3_n,
            heater_layer,
            q_ls_n_prev_heat_source,
            temp_setpntmax,
            simtime,
        )?;
        if soc_lower > soc_target
            || relative_eq!(soc_lower, soc_target, epsilon = 1e-10, max_relative = 1e-9)
        {
            return Ok(energy_lower);
        }

        let soc_upper = self.realised_state_of_charge(
            energy_upper,
            temp_s3_n,
            heater_layer,
            q_ls_n_prev_heat_source,
            temp_setpntmax,
            simtime,
        )?;
        if soc_upper < soc_target
            || relative_eq!(soc_upper, soc_target, epsilon = 1e-10, max_relative = 1e-9)
        {
            return Ok(energy_upper);
        }

        // Bisection on the charge: the realised state of charge rises as the charge
        // increases. energy_lower is the largest charge whose realised state of
        // charge is still below the target; energy_upper the smallest that reaches
        // it. 50 halvings settle the charge far below any precision a tank size
        // needs. Return energy_upper (the smallest charge that reaches the target):
        // the state of charge steps rather than varying smoothly - it jumps as each
        // layer crosses the usable temperature - so where the target falls inside a
        // step it is reached by charging through to the next step rather than
        // stopping short of it.
        for _ in 0..50 {
            let energy_mid = 0.5 * (energy_lower + energy_upper);
            if self.realised_state_of_charge(
                energy_mid,
                temp_s3_n,
                heater_layer,
                q_ls_n_prev_heat_source,
                temp_setpntmax,
                simtime,
            )? < soc_target
            {
                energy_lower = energy_mid;
            } else {
                energy_upper = energy_mid;
            }
        }

        Ok(energy_upper)
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

        self.temps_after_pumping(volume_pumped, &mut remaining_vols, &tank_layer_temperatures)
    }

    /// Calculate the temperatures of the tank after volume is pumped
    fn temps_after_pumping(
        &self,
        volume_pumped: f64,
        remaining_vols: &mut [f64],
        tank_layer_temperatures: &[f64],
    ) -> anyhow::Result<Vec<f64>> {
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
                if volume_pumped_remaining < 0.
                    || relative_eq!(
                        volume_pumped_remaining,
                        0.,
                        epsilon = 1e-10,
                        max_relative = 1e-9
                    )
                {
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
                if needed_volume < 0.
                    || relative_eq!(needed_volume, 0., epsilon = 1e-10, max_relative = 1e-9)
                {
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
                        if needed_volume < 0.
                            || relative_eq!(
                                needed_volume,
                                0.,
                                max_relative = 1e-09,
                                epsilon = 1e-10
                            )
                        {
                            break;
                        }
                    }
                }

                if remaining_vols[i] != self.storage_tank.vol_n[i] {
                    bail!("Volume mismatch in layer {i}");
                }

                // Temperatures after moving
                // ----------------
                // After moving water to a layer, calculate the new temperature
                // for the current layer based on vol and temperature.
                tank_layer_temperatures[i] = volume_weighted_temperature / remaining_vols[i];
            }
        }

        Ok(tank_layer_temperatures)
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

        let _setpnt = self
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

        if top_layer_temp < temp_target
            || relative_eq!(top_layer_temp, temp_target, max_relative = 1e-09)
            || qin < 0.
            || relative_eq!(qin, 0., max_relative = 1e-09, epsilon = 1e-10)
        {
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

            let numerator = FSum::with_all(
                temp_s7_n
                    .iter()
                    .zip(temp_factors.iter())
                    .map(|(&t, &f)| t * f),
            )
            .value()
                - temp_target * FSum::with_all(temp_factors.iter().copied()).value();
            let denominator = temp_target - temp_s7_n[current_layer];
            temp_factors[current_layer] = if denominator < 0.
                || relative_eq!(denominator, 0., max_relative = 1e-09, epsilon = 1e-10)
            {
                // If the current layer is at or above target temperature, pump all of it
                1.0
            } else {
                // Calculate the fraction of the current layer to be pumped
                numerator / denominator
            };

            if temp_factors[current_layer] < 1. {
                // If we don't need to pump the entire layer, stop iteration
                break;
            }
            // If entire layer needs to be pumped, set factor to 1
            // and continue to next layer
            temp_factors[current_layer] = 1.0;
        }

        // Calculate volume to be pumped (only from layers below heater_layer)
        let volume_pumped = FSum::with_all(
            bottom_volumes
                .iter()
                .zip(temp_factors[..heater_layer].iter())
                .map(|(&v, &f)| v * f),
        )
        .value();

        // Check that the volume pumped doesn't exceed the volume of water up to the heater layer
        if volume_pumped > FSum::with_all(bottom_volumes).value() {
            bail!("Volume pumped is higher than total bottom volumes");
        }

        // Cap volume pumped based on pump max flow rate in timestep
        let max_volume_pumped = self.max_flow_rate_pump_l_per_min
            * self.storage_tank.simulation_timestep
            * MINUTES_PER_HOUR as f64;

        let volume_pumped = volume_pumped.min(max_volume_pumped);

        Ok(volume_pumped)
    }

    fn get_losses_from_primary_pipework_and_storage(&self) -> (f64, f64) {
        self.storage_tank
            .get_losses_from_primary_pipework_and_storage()
    }

    fn testoutput(
        &self,
        usage_events: &[WaterEventResult],
        volume_extracted: f64,
        q_use_w: f64,
        q_unmet_w: f64,
        temp_ini_n: &[f64],
        temp_s3_n: &[f64],
        q_x_in_n: &[f64],
        q_s6: f64,
        temp_s6_n: &[f64],
        temp_s7_n: &[f64],
        q_in_h_w: f64,
        q_ls: f64,
        temp_s8_n: &[f64],
        temp_average: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<()> {
        let mut detailed_output = match self.storage_tank.detailed_results.as_ref() {
            None => return Ok(()),
            Some(detailed_output) => detailed_output.write(),
        };

        let demand = summarise_events(usage_events);

        // Calculates state of charge
        let state_of_charge_draw_off =
            (self.calc_state_of_charge(temp_s3_n, simtime)? * 1e6).round() / 1e6;
        let state_of_charge_final =
            (self.calc_state_of_charge(temp_s8_n, simtime)? * 1e6).round() / 1e6;

        if simtime.index == 0 {
            #[derive(Clone, Copy, Debug, Default)]
            enum WithIndex {
                #[default]
                Yes,
                No,
            }

            fn header_dup(header: &str, n: usize, with_index: WithIndex) -> Vec<StringOrNumber> {
                (0..n)
                    .map(|i| match with_index {
                        WithIndex::Yes => format!("{header} {}", i + 1).into(),
                        WithIndex::No => header.into(),
                    })
                    .collect()
            }

            let mut header_row: Vec<StringOrNumber> = [
                "time",
                "volume total",
                "specific heat",
                "density",
                "cold water",
                "events",
                "volume extracted",
            ]
            .into_iter()
            .map(Into::into)
            .collect();
            header_row.extend(header_dup(
                "initial temp.",
                temp_ini_n.len(),
                WithIndex::Yes,
            ));
            header_row.extend(["energy withdrawn", "energy unmet"].map(Into::into));
            header_row.extend(header_dup(
                "temp. after volume withdrawn",
                temp_s3_n.len(),
                WithIndex::Yes,
            ));
            header_row.push("state of charge after volume withdrawn".into());
            header_row.extend(header_dup(
                "potential energy input",
                q_x_in_n.len(),
                WithIndex::Yes,
            ));
            header_row.push("theoretical energy stored after energy input".into());
            header_row.extend(header_dup(
                "theoretical temp. after energy input",
                temp_s6_n.len(),
                WithIndex::Yes,
            ));
            header_row.extend(header_dup(
                "temp. after volume mixing",
                temp_s7_n.len(),
                WithIndex::Yes,
            ));
            header_row.extend(
                ["energy input (adjusted)", "thermal losses"]
                    .into_iter()
                    .map(Into::into),
            );
            header_row.extend(header_dup(
                "temp. after thermal losses",
                temp_s8_n.len(),
                WithIndex::Yes,
            ));
            header_row.extend(
                ["temp_average_drawoff", "state of charge final"]
                    .into_iter()
                    .map(Into::into),
            );

            detailed_output.push(header_row);

            let mut units_row: Vec<StringOrNumber> = [
                "h",
                "litres",
                "kWh/kgK",
                "kg/l",
                "oC",
                "Type: litres hot (litres @ oC)",
                "litres",
            ]
            .into_iter()
            .map(Into::into)
            .collect();
            units_row.extend(header_dup("oC", temp_ini_n.len(), WithIndex::No));
            units_row.extend(["kWh", "kWh"].into_iter().map(Into::into));
            units_row.extend(header_dup("oC", temp_s3_n.len(), WithIndex::No));
            units_row.extend(["fraction"].into_iter().map(Into::into));
            units_row.extend(header_dup("kWh", q_x_in_n.len(), WithIndex::No));
            units_row.extend(["kWh"].into_iter().map(Into::into));
            units_row.extend(header_dup("oC", temp_s6_n.len(), WithIndex::No));
            units_row.extend(header_dup("oC", temp_s7_n.len(), WithIndex::No));
            units_row.extend(["kWh", "kWh"].into_iter().map(Into::into));
            units_row.extend(header_dup("oC", temp_s8_n.len(), WithIndex::No));
            units_row.extend(["oC", "fraction"].into_iter().map(Into::into));

            detailed_output.push(units_row);
        }

        let temp_cold_water: StringOrNumber =
            if relative_eq!(volume_extracted, 0.0, epsilon = 1e-10, max_relative = 1e-9) {
                "".into()
            } else {
                let list_temp_vol = self
                    .storage_tank
                    .cold_feed
                    .get_temp_cold_water(volume_extracted, simtime)?;
                (FSum::with_all(list_temp_vol.iter().map(|(t, v)| t * v)).value()
                    / FSum::with_all(list_temp_vol.iter().map(|(_, v)| v)).value())
                .into()
            };

        let mut values_row: Vec<StringOrNumber> = vec![
            simtime.hour_of_day().into(),
            self.storage_tank.volume_total_in_litres.into(),
            self.storage_tank.cp.into(),
            self.storage_tank.rho.into(),
            temp_cold_water,
            demand.into(),
            volume_extracted.into(),
        ];
        values_row.extend(temp_ini_n.iter().map(|t| t.into()));
        values_row.extend([q_use_w.into(), q_unmet_w.into()]);
        values_row.extend(temp_s3_n.iter().map(|t| t.into()));
        values_row.push(state_of_charge_draw_off.into());
        values_row.extend(q_x_in_n.iter().map(|q| q.into()));
        values_row.push(q_s6.into());
        values_row.extend(temp_s6_n.iter().map(|t| t.into()));
        values_row.extend(temp_s7_n.iter().map(|t| t.into()));
        values_row.extend([q_in_h_w.into(), q_ls.into()]);
        values_row.extend(temp_s8_n.iter().map(|t| t.into()));
        values_row.push(temp_average.into());
        values_row.push(state_of_charge_final.into());

        detailed_output.push(values_row);

        Ok(())
    }

    pub(crate) fn get_temp_hot_water(
        &self,
        volume_req: f64,
        volume_req_already: Option<f64>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        self.storage_tank.get_temp_hot_water(
            volume_req,
            volume_req_already,
            simulation_time_iteration,
        )
    }

    pub(crate) fn draw_off_hot_water(
        &self,
        volume: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<(Option<f64>, f64)> {
        self.storage_tank
            .draw_off_hot_water(volume, simulation_time_iteration)
    }

    pub(crate) fn draw_off_water(
        &self,
        volume_needed: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        self.storage_tank
            .draw_off_water(volume_needed, simulation_time_iteration)
    }

    pub(crate) fn get_temp_cold_water(
        &self,
        volume_needed: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        // TODO this matches Python - is it correct?
        self.storage_tank
            .get_temp_hot_water(volume_needed, None, simulation_time_iteration)
    }
}

#[derive(Debug)]
pub struct ImmersionHeater {
    pwr: f64, // rated power
    energy_supply_connection: EnergySupplyConnection,
    simulation_timestep: f64,
    control: Option<Arc<RangeTimeControl>>,
    diverter: ArcSwapOption<RwLock<PVDiverter>>,
}

/// An object to represent an immersion heater
impl ImmersionHeater {
    /// Construct an ImmersionHeater object
    /// Arguments:
    /// * `rated_power` - in kW
    /// * `energy_supply_conn`- reference to EnergySupplyConnection object
    /// * `simulation_time` - reference to SimulationTime object
    /// * `controlmin` - reference to a control object which must select current
    ///                  the minimum timestep temperature
    /// * `controlmax` - reference to a control object which must select current
    ///                  the maximum timestep temperature
    /// * `control` - Reference to a RangeTimeControl object, combining controlmax and controlmin.
    ///               Takes precedence if set.
    /// * `diverter` - reference to a PV diverter object
    pub(crate) fn new(
        rated_power: f64,
        energy_supply_connection: EnergySupplyConnection,
        simulation_timestep: f64,
        control: Option<Arc<RangeTimeControl>>,
    ) -> Self {
        Self {
            pwr: rated_power,
            energy_supply_connection,
            simulation_timestep,
            control,
            diverter: Default::default(),
        }
    }
    pub(crate) fn setpnt(&self, simtime: SimulationTimeIteration) -> (Option<f64>, Option<f64>) {
        if let Some(control) = self.control.as_ref() {
            control.setpnt_range_time_control(&simtime)
        } else {
            (None, None)
        }
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
        ignore_standard_ctrl: Option<bool>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let _ignore_standard_ctrl = ignore_standard_ctrl.unwrap_or(false); // TODO 1.0.0a9 migration
        if energy_demand < 0.0 {
            bail!("Negative energy demand on ImmersionHeater");
        };

        let energy_supplied =
            if self.control.is_none() || self.control.as_ref().unwrap().is_on(&simtime) {
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
        if self.control.is_some() && self.control.as_ref().unwrap().is_on(&simtime)
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
pub enum HotWaterStorageTank {
    StorageTank(Arc<RwLock<StorageTank>>),
    SmartHotWaterTank(Arc<RwLock<SmartHotWaterTank>>),
    #[cfg(test)]
    Mock(Box<HotWaterSourceMockKind>),
}

impl HotWaterSourceBehaviour for HotWaterStorageTank {
    fn get_cold_water_source(&self) -> WaterSupply {
        match self {
            HotWaterStorageTank::StorageTank(storage_tank) => {
                storage_tank.read().get_cold_water_source().clone()
            }
            HotWaterStorageTank::SmartHotWaterTank(smart_storage_tank) => {
                smart_storage_tank.read().get_cold_water_source().clone()
            }
            #[cfg(test)]
            HotWaterStorageTank::Mock(source) => source.get_cold_water_source(),
        }
    }

    fn demand_hot_water(
        &self,
        usage_events: Vec<WaterEventResult>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        match &self {
            HotWaterStorageTank::StorageTank(rw_lock) => rw_lock
                .read()
                .demand_hot_water(usage_events.into(), simtime),
            HotWaterStorageTank::SmartHotWaterTank(rw_lock) => rw_lock
                .read()
                .demand_hot_water(usage_events.into(), simtime),
            #[cfg(test)]
            HotWaterStorageTank::Mock(source) => source.demand_hot_water(usage_events, simtime),
        }
    }

    fn get_temp_hot_water(
        &self,
        volume_required: f64,
        volume_required_already: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<Vec<(f64, f64)>> {
        match self {
            HotWaterStorageTank::StorageTank(rw_lock) => rw_lock.read().get_temp_hot_water(
                volume_required,
                Some(volume_required_already),
                simtime,
            ),
            HotWaterStorageTank::SmartHotWaterTank(rw_lock) => rw_lock.read().get_temp_hot_water(
                volume_required,
                Some(volume_required_already),
                simtime,
            ),
            #[cfg(test)]
            HotWaterStorageTank::Mock(source) => {
                source.get_temp_hot_water(volume_required, volume_required_already, simtime)
            }
        }
    }

    fn internal_gains(&self) -> Option<f64> {
        match self {
            HotWaterStorageTank::StorageTank(rw_lock) => Some(rw_lock.read().internal_gains()),
            HotWaterStorageTank::SmartHotWaterTank(rw_lock) => {
                Some(rw_lock.read().internal_gains())
            }
            #[cfg(test)]
            HotWaterStorageTank::Mock(source) => source.internal_gains(),
        }
    }

    fn get_losses_from_primary_pipework_and_storage(&self) -> (f64, f64) {
        match self {
            HotWaterStorageTank::StorageTank(rw_lock) => rw_lock
                .read()
                .get_losses_from_primary_pipework_and_storage(),
            HotWaterStorageTank::SmartHotWaterTank(rw_lock) => rw_lock
                .read()
                .get_losses_from_primary_pipework_and_storage(),
            #[cfg(test)]
            HotWaterStorageTank::Mock(source) => {
                source.get_losses_from_primary_pipework_and_storage()
            }
        }
    }

    fn is_point_of_use(&self) -> bool {
        false
    }
}

impl HotWaterStorageTank {
    pub(crate) fn ultimate_cold_water_source(&self) -> WaterSupply {
        match self {
            HotWaterStorageTank::StorageTank(rw_lock) => {
                let cold_feed = rw_lock.read().cold_feed.clone();
                if let WaterSupply::Preheated(tank) = cold_feed {
                    tank.ultimate_cold_water_source()
                } else {
                    cold_feed
                }
            }
            HotWaterStorageTank::SmartHotWaterTank(rw_lock) => {
                let cold_feed = rw_lock.read().storage_tank.cold_feed.clone().clone();
                if let WaterSupply::Preheated(tank) = cold_feed {
                    tank.ultimate_cold_water_source()
                } else {
                    cold_feed
                }
            }
            #[cfg(test)]
            HotWaterStorageTank::Mock(_) => WaterSupply::Mock(MockWaterSupply::default()),
        }
    }

    fn additional_energy_input(
        &self,
        heat_source: &HeatSource,
        heat_source_name: &str,
        energy_input: f64,
        control_max_diverter: Option<&Control>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        match self {
            HotWaterStorageTank::StorageTank(storage_tank) => {
                storage_tank.read().additional_energy_input(
                    heat_source,
                    heat_source_name,
                    energy_input,
                    control_max_diverter,
                    simulation_time_iteration,
                )
            }
            HotWaterStorageTank::SmartHotWaterTank(smart_hot_water_tank) => {
                smart_hot_water_tank.read().additional_energy_input(
                    heat_source,
                    heat_source_name,
                    energy_input,
                    control_max_diverter,
                    simulation_time_iteration,
                )
            }
            #[cfg(test)]
            HotWaterStorageTank::Mock(_source) => Ok(0.),
        }
    }
}

#[derive(Debug)]
pub struct PVDiverter {
    pre_heated_water_source: HotWaterStorageTank,
    immersion_heater: Arc<Mutex<ImmersionHeater>>,
    heat_source_name: ArcStr,
    control_max: Option<Control>,
    capacity_used: AtomicF64,
}

impl PVDiverter {
    pub(crate) fn new(
        storage_tank: &HotWaterStorageTank,
        heat_source: Arc<Mutex<ImmersionHeater>>,
        heat_source_name: ArcStr,
        control_max: Option<Control>,
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
        let energy_diverted = self.pre_heated_water_source.additional_energy_input(
            &HeatSource::Storage(HeatSourceWithStorageTank::Immersion(
                self.immersion_heater.clone(),
            )),
            &self.heat_source_name,
            energy_diverted_max,
            self.control_max.as_ref(),
            simulation_time_iteration,
        )?;

        Ok(energy_diverted)
    }
}

/// The following code contains objects that represent solar thermal systems.
/// Method 3 in BS EN 15316-4-3:2017.
/// An object to represent a solar thermal system
#[derive(Educe)]
#[educe(Debug)]
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
    orientation: Orientation360,
    solar_loop_piping_hlc: f64,
    external_conditions: Arc<ExternalConditions>,
    simulation_timestep: f64,
    #[educe(Debug(ignore))]
    temp_internal_air_fn: TempInternalAirFn,
    heat_output_collector_loop: AtomicF64,
    energy_supplied: AtomicF64,
    control_max: Control,
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
    /// * orientation     -- The orientation angle of the inclined surface, expressed as the geographical azimuth angle of the
    ///                 horizontal projection of the inclined surface normal, 0 to 360, in degrees;
    ///                 Needed to calculate solar irradiation at the panel surface.
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
        orientation: Orientation360,
        solar_loop_piping_hlc: f64,
        external_conditions: Arc<ExternalConditions>,
        temp_internal_air_fn: TempInternalAirFn,
        simulation_timestep: f64,
        control_max: Control,
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
    use crate::core::controls::time_control::{
        MockControl, ScheduleOrControl, SetpointTimeControl,
    };
    use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyBuilder};
    use crate::core::material_properties::WATER;
    use crate::core::pipework::PipeworkLocation;
    use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
    use crate::core::water_heat_demand::misc::WaterEventResultType;
    use crate::corpus::HeatSource;
    use crate::external_conditions::{
        DaylightSavingsConfig, ShadingObject, ShadingObjectType, ShadingSegment,
    };
    use crate::input::{FuelType, PipeworkContents, WaterPipeworkLocation};
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
        let wind_directions = vec![0.0; 8].into_iter().map(Into::into).collect();
        let diffuse_horizontal_radiations = vec![333., 610., 572., 420., 0., 10., 90., 275.];
        let direct_beam_radiations = vec![420., 750., 425., 500., 0., 40., 0., 388.];
        let solar_reflectivity_of_ground = vec![0.2; 8760];
        let shading_segments = vec![
            ShadingSegment {
                start360: Orientation360::create_from_180(180.).unwrap(),
                end360: Orientation360::create_from_180(135.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(135.).unwrap(),
                end360: Orientation360::create_from_180(90.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(90.).unwrap(),
                end360: Orientation360::create_from_180(45.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(45.).unwrap(),
                end360: Orientation360::create_from_180(0.).unwrap(),
                shading_objects: vec![ShadingObject {
                    object_type: ShadingObjectType::Obstacle,
                    height: 10.5,
                    distance: 120.,
                }],
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(0.).unwrap(),
                end360: Orientation360::create_from_180(-45.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(-45.).unwrap(),
                end360: Orientation360::create_from_180(-90.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(-90.).unwrap(),
                end360: Orientation360::create_from_180(-135.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(-135.).unwrap(),
                end360: Orientation360::create_from_180(-180.).unwrap(),
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
            simulation_time_for_storage_tank.iter().total_steps(),
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
        let control_min = ScheduleOrControl::Schedule(control_min_schedule);

        let control_max = ScheduleOrControl::Schedule(control_max_schedule);

        let immersion_heater = ImmersionHeater::new(
            rated_power,
            energy_supply_connection.clone(),
            simulation_timestep,
            Some(Arc::from(
                RangeTimeControl::new(
                    control_min,
                    control_max,
                    simulation_time_for_storage_tank.iter(),
                    0,
                    1.,
                    None,
                )
                .unwrap(),
            )),
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
        energy_supply: Arc<RwLock<EnergySupply>>,
    ) -> (StorageTank, Arc<RwLock<EnergySupply>>) {
        let cold_feed = WaterSupply::ColdWaterSource(cold_water_source.clone());
        let simtime = simulation_time_for_storage_tank.iter().current_iteration();

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

        let energy_supply_connection =
            EnergySupply::connection(energy_supply.clone(), "immersion").unwrap();

        let heat_source_imheater = heat_source(
            simulation_time_for_storage_tank,
            energy_supply_connection.clone(),
            50.0,
            0.1,
            0.33,
            control_min_schedule,
            control_max_schedule,
        );

        let heat_sources = IndexMap::from([("imheater".into(), heat_source_imheater)]);
        let storage_tank = StorageTank::new(
            150.0,
            1.68,
            55.0,
            cold_feed,
            &simtime,
            heat_sources,
            temp_internal_air_fn.clone(),
            external_conditions.clone(),
            false,
            None,
            None,
            *WATER,
            None,
            None,
            None,
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
                simulation_time_for_storage_tank.iter().total_steps(),
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

        let cold_feed = WaterSupply::ColdWaterSource(cold_water_source.clone());
        let simtime = simulation_time_for_storage_tank.iter().current_iteration();

        let heat_sources = IndexMap::from([("imheater2".into(), heat_source)]);
        let storage_tank = StorageTank::new(
            210.0,
            1.61,
            60.0,
            cold_feed,
            &simtime,
            heat_sources,
            temp_internal_air_fn.clone(),
            external_conditions.clone(),
            false,
            None,
            None,
            *WATER,
            None,
            None,
            None,
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
        let wind_directions = vec![0.0; 8].into_iter().map(Into::into).collect();
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
                start360: Orientation360::create_from_180(180.).unwrap(),
                end360: Orientation360::create_from_180(135.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(135.).unwrap(),
                end360: Orientation360::create_from_180(90.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(90.).unwrap(),
                end360: Orientation360::create_from_180(45.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(45.).unwrap(),
                end360: Orientation360::create_from_180(0.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(0.).unwrap(),
                end360: Orientation360::create_from_180(-45.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(-45.).unwrap(),
                end360: Orientation360::create_from_180(-90.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(-90.).unwrap(),
                end360: Orientation360::create_from_180(-135.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(-135.).unwrap(),
                end360: Orientation360::create_from_180(-180.).unwrap(),
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
        let cold_feed = WaterSupply::ColdWaterSource(Arc::new(ColdWaterSource::new(
            cold_water_temps,
            start_day,
            time_series_step,
        )));
        let simtime = simulation_time_for_storage_tank.iter().current_iteration();

        let heat_sources = IndexMap::from([("imheater".into(), heat_source)]);

        StorageTank::new(
            150.0,
            1.68,
            55.0,
            cold_feed,
            &simtime,
            heat_sources,
            temp_internal_air_fn.clone(),
            external_conditions_for_pv_diverter.clone(),
            false,
            None,
            None,
            *WATER,
            None,
            None,
            None,
        )
        .unwrap()
    }

    #[fixture]
    fn diverter_control() -> Control {
        Control::SetpointTime(
            SetpointTimeControl::new(
                vec![Some(60.), Some(60.), Some(60.), Some(60.)],
                0,
                1.,
                None,
                None,
                1.,
            )
            .into(),
        )
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
        ]
        .into_iter()
        .map(Into::into)
        .collect();
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
                start360: Orientation360::create_from_180(180.).unwrap(),
                end360: Orientation360::create_from_180(135.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(135.).unwrap(),
                end360: Orientation360::create_from_180(90.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(90.).unwrap(),
                end360: Orientation360::create_from_180(45.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(45.).unwrap(),
                end360: Orientation360::create_from_180(0.).unwrap(),
                shading_objects: vec![ShadingObject {
                    object_type: ShadingObjectType::Obstacle,
                    height: 10.5,
                    distance: 12.,
                }],
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(0.).unwrap(),
                end360: Orientation360::create_from_180(-45.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(-45.).unwrap(),
                end360: Orientation360::create_from_180(-90.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(-90.).unwrap(),
                end360: Orientation360::create_from_180(-135.).unwrap(),
                ..Default::default()
            },
            ShadingSegment {
                start360: Orientation360::create_from_180(-135.).unwrap(),
                end360: Orientation360::create_from_180(-180.).unwrap(),
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
        let cold_feed = WaterSupply::ColdWaterSource(Arc::new(ColdWaterSource::new(
            cold_water_temps.to_vec(),
            212,
            1.,
        )));
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
            Orientation360::create_from_180(0.).unwrap(),
            0.5,
            external_conditions_for_solar_thermal.clone(),
            temp_internal_air_fn.clone(),
            simulation_time_for_solar_thermal.step,
            Control::SetpointTime(control_max.into()),
            *WATER,
            None,
        )));

        let storage_tank = StorageTank::new(
            150.0,
            1.68,
            55.0,
            cold_feed,
            &simulation_time_for_solar_thermal.iter().current_iteration(),
            IndexMap::from([(
                "solthermal".into(),
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
            false,
            None,
            None,
            *WATER,
            None,
            None,
            None,
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
        let event_data = get_event_data_immersion();

        //  Expected results for the unit test
        let expected_temperatures_1 = [
            [55.0, 55.0, 55.0, 55.0],
            [
                15.448000000000006,
                54.595555555555556,
                54.595555555555556,
                54.595555555555556,
            ],
            [
                15.448000000000006,
                54.19530534979424,
                54.19530534979424,
                54.19530534979424,
            ],
            [
                10.5,
                15.39537857738237,
                53.39140601130916,
                53.79920588690748,
            ],
            [55.0, 55.0, 55.0, 55.0],
            [
                13.400000000000002,
                54.595555555555556,
                54.595555555555556,
                54.595555555555556,
            ],
            [
                13.400000000000002,
                54.19530534979424,
                54.19530534979424,
                54.19530534979424,
            ],
            [
                13.400000000000002,
                53.79920588690749,
                53.79920588690749,
                53.79920588690749,
            ],
        ];

        // Also test case where heater does not heat all layers, to ensure this is handled correctly

        let expected_temperatures_2 = [
            [10.0, 24.55607367670878, 60.0, 60.0],
            [
                10.056616089321501,
                16.312757933665832,
                39.76314001211786,
                59.687654320987654,
            ],
            [
                10.056616089321501,
                16.31053773845771,
                39.59445105524171,
                59.37752591068435,
            ],
            [
                10.342751590114434,
                12.274601927758184,
                24.50747434931655,
                46.393323819752034,
            ],
            [10.342751590114434, 12.274601927758184, 60.0, 60.0],
            [10.741316223514428, 11.103100848370719, 60.0, 60.0],
            [
                10.741316223514428,
                11.103100848370719,
                59.687654320987654,
                59.687654320987654,
            ],
            [
                10.741316223514428,
                11.103100848370719,
                59.37752591068435,
                59.37752591068435,
            ],
        ];

        let expected_energy_supplied_1 = [5.9141614815, 0.0, 0.0, 0.0, 3.8585103966, 0.0, 0.0, 0.0];

        let expected_energy_supplied_2 = [
            0.6719988461,
            0.0,
            0.0,
            0.0,
            3.0339862161,
            1.8040212455,
            0.0,
            0.0,
        ];

        // Loop through the timesteps and the associated data pairs using `subTest`
        for (t_idx, t_it) in simulation_time_for_storage_tank.iter().enumerate() {
            let usage_events = event_data[t_idx].clone();

            storage_tank1
                .demand_hot_water(usage_events.clone(), t_it)
                .unwrap();

            // Verify the temperatures against expected results
            let actual_temperatues_1 = storage_tank1.temp_n.read().clone();
            for i in 0..actual_temperatues_1.len() {
                // TODO decrease max_relative here
                assert_relative_eq!(
                    actual_temperatues_1[i],
                    expected_temperatures_1[t_idx][i],
                    max_relative = 1e-2
                );
            }

            assert_relative_eq!(
                energy_supply1.read().results_by_end_user()["immersion"][t_idx],
                expected_energy_supplied_1[t_idx],
                max_relative = 1e-3
            );

            let temp_hot = if t_idx == 0 {
                60.
            } else {
                expected_temperatures_2[t_idx - 1][3]
            };

            // Convert usage events based on HW temp of 55 to equivalent 60:
            let mut usage_events2 = vec![];

            if usage_events.is_some() {
                for event in usage_events.unwrap() {
                    let volume_hot = event.volume_warm
                        * (event.temperature_warm - COLD_WATER_TEMPS[t_idx])
                        / (temp_hot - COLD_WATER_TEMPS[t_idx]);
                    usage_events2.push(WaterEventResult {
                        event_result_type: event.event_result_type,
                        temperature_warm: event.temperature_warm,
                        volume_warm: event.volume_warm,
                        volume_hot,
                        event_duration: 0.,
                    });
                }
            }

            storage_tank2
                .demand_hot_water(Some(usage_events2), t_it)
                .unwrap();

            let actual_temperatues_2 = storage_tank2.temp_n.read().clone();
            for i in 0..actual_temperatues_2.len() {
                assert_relative_eq!(
                    actual_temperatues_2[i],
                    expected_temperatures_2[t_idx][i],
                    max_relative = 1e-7
                );
            }

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
                storage_tank1.temperature_surrounding_primary_pipework(&pipework, &t_it),
                [0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0][t_idx]
            );
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
                storage_tank1.temperature_surrounding_primary_pipework(&pipework, &t_it),
                [20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0][t_idx]
            );
        }
    }

    #[rstest]
    fn test_get_cold_water_source(storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>)) {
        let (storage_tank1, _) = storage_tank1;
        let result = storage_tank1.get_cold_water_source();

        assert!(matches!(result, WaterSupply::ColdWaterSource(_)));
    }

    #[rstest]
    fn test_get_temp_hot_water(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        let simtime = simulation_time_for_storage_tank.iter().current_iteration();
        let (storage_tank1, _) = storage_tank1;
        let expected = vec![(55.0, 37.5), (55.0, 37.5), (55.0, 25.0)];
        assert_eq!(
            storage_tank1
                .get_temp_hot_water(100.0, None, simtime)
                .unwrap(),
            expected
        );
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
        );
    }

    #[rstest]
    fn test_thermal_losses(storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>)) {
        // The tank has four layers (default nb_vol), so the fixtures supply one value per layer.
        let (storage_tank1, _) = storage_tank1;
        let temp_s7_n = [12.0, 18.0, 25.0, 32.0];
        let q_x_in_n = [0., 1., 2., 3.];
        let q_h_sto_s7 = [0.1, 0.2, 0.3, 0.4];
        let heater_layer = 2;
        let q_ls_n_prev_heat_source = [0.0, 0.0, 0.0, 0.0];
        let setpntmax = 55.0;

        // Every layer sits below the setpoint, so the over-setpoint cap never applies and the
        // loss is independent of any temperature rise; pass temp_s3_n == temp_s7_n.
        assert_eq!(
            storage_tank1.calc_temps_after_thermal_losses(
                &temp_s7_n,
                &temp_s7_n,
                &q_x_in_n,
                &q_h_sto_s7,
                heater_layer,
                &q_ls_n_prev_heat_source,
                setpntmax.into()
            ),
            (
                6.0,
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

    /// When a layer is above the setpoint but this heat source did not warm it, the loss is
    /// computed at the layer's actual temperature, not the setpoint.
    ///
    /// Each layer ends at the same 60°C it started at (temp_s7_n == temp_s3_n), so this source
    /// did not warm it and is not holding it at the setpoint - the excess pre-dates it. Each
    /// layer therefore loses heat against (60 - ambient) and cools below 60°C. Were the loss
    /// capped at the 55°C setpoint it would be underestimated and the layer left too warm. (The
    /// tank's loss-reference ambient is 16°C, so the loss is proportional to 60 - 16 = 44 K
    /// rather than 55 - 16 = 39 K.)
    #[rstest]
    fn test_thermal_losses_above_setpoint_without_heat_input(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
    ) {
        let (storage_tank1, _) = storage_tank1;

        assert_eq!(
            storage_tank1.calc_temps_after_thermal_losses(
                &[60.0, 60.0, 60.0, 60.0],
                &[60.0, 60.0, 60.0, 60.0],
                &[0.0, 0.0, 0.0, 0.0],
                &[0.0, 0.0, 0.0, 0.0],
                2,
                &[0.0, 0.0, 0.0, 0.0],
                Some(55.)
            ),
            (
                0.0,
                0.07954765432098766,
                vec![
                    59.543703703703706,
                    59.543703703703706,
                    59.543703703703706,
                    59.543703703703706
                ],
                vec![
                    0.019886913580246916,
                    0.019886913580246916,
                    0.019886913580246916,
                    0.019886913580246916
                ]
            )
        );
    }

    /// When this heat source warmed a layer above the setpoint, it is held at the setpoint
    /// (Case 2) and its loss is computed at the setpoint - unchanged by the fix.
    ///
    /// Each layer rose from 10°C to 60°C this timestep (temp_s7_n > temp_s3_n), so the source
    /// is holding it at the setpoint: each 60°C layer is clamped to the 55°C setpoint and loses
    /// heat against (55 - ambient), confirming the over-setpoint cap still applies in this case.
    #[rstest]
    fn test_thermal_losses_above_setpoint_with_heat_input(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
    ) {
        let (storage_tank1, _) = storage_tank1;

        assert_eq!(
            storage_tank1.calc_temps_after_thermal_losses(
                &[10.0, 10.0, 10.0, 10.0],
                &[60.0, 60.0, 60.0, 60.0],
                &[0.0, 0.0, 1.0, 0.0],
                &[0.1, 0.2, 0.3, 0.4],
                2,
                &[0.0, 0.0, 0.0, 0.0],
                Some(55.)
            ),
            (
                1.0,
                0.07050814814814815,
                vec![55.0, 55.0, 55.0, 55.0],
                vec![
                    0.01762703703703704,
                    0.01762703703703704,
                    0.01762703703703704,
                    0.01762703703703704
                ]
            )
        );
    }

    /// A warmed layer within tolerance of the setpoint is not clamped to it.
    ///
    /// Each layer rose from 10°C this timestep (so it is warmed) and ends a fraction
    /// above the 55°C setpoint, within abs_tol=1e-10. It is therefore treated as at the
    /// setpoint and loses heat against its actual temperature (Case 1), settling below
    /// 55°C - the same as a layer warmed to exactly the setpoint. Without the tolerance
    /// the layer would be clamped to exactly 55.0 (Case 2), so the result would depend on
    /// the last bit of the computed layer temperature.
    #[rstest]
    fn test_thermal_losses_warmed_within_tolerance_of_setpoint_not_clamped(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
    ) {
        let (storage_tank1, _) = storage_tank1;

        assert_eq!(
            storage_tank1.calc_temps_after_thermal_losses(
                &[10.0, 10.0, 10.0, 10.0],
                &[55.0 + 5e-11, 55.0 + 5e-11, 55.0 + 5e-11, 55.0 + 5e-11],
                &[0.0, 0.0, 1.0, 0.0],
                &[0.1, 0.2, 0.3, 0.4],
                2,
                &[0.0, 0.0, 0.0, 0.0],
                Some(55.)
            ),
            (
                1.0,
                0.07050814814814815,
                vec![
                    54.59555555560556,
                    54.59555555560556,
                    54.59555555560556,
                    54.59555555560556
                ],
                vec![
                    0.01762703703703704,
                    0.01762703703703704,
                    0.01762703703703704,
                    0.01762703703703704
                ]
            )
        );
    }

    /// The uncapped loss is still reduced by losses already attributed to an earlier heat
    /// source, avoiding double-counting across heat sources.
    ///
    /// Each layer ends at the 60°C it started at (temp_s7_n == temp_s3_n), so the loss is taken
    /// at the actual temperature. The gross loss at 60°C (0.019886... kWh) has the previous heat
    /// source's per-layer loss (0.01 kWh) subtracted, leaving 0.009886... kWh per layer.
    #[rstest]
    fn test_thermal_losses_above_setpoint_prev_heat_source(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
    ) {
        let (storage_tank1, _) = storage_tank1;

        assert_eq!(
            storage_tank1.calc_temps_after_thermal_losses(
                &[60.0, 60.0, 60.0, 60.0],
                &[60.0, 60.0, 60.0, 60.0],
                &[0.0, 0.0, 0.0, 0.0],
                &[0.0, 0.0, 0.0, 0.0],
                2,
                &[0.01, 0.01, 0.01, 0.01],
                Some(55.)
            ),
            (
                0.0,
                0.03954765432098766,
                vec![
                    59.773149210395864,
                    59.773149210395864,
                    59.773149210395864,
                    59.773149210395864
                ],
                vec![
                    0.009886913580246915,
                    0.009886913580246915,
                    0.009886913580246915,
                    0.009886913580246915,
                ]
            )
        );
    }

    /// A single source heating only some layers must not cap the losses of layers it did
    /// not warm.
    ///
    /// The lower two layers rise from 10°C to 55°C this timestep, so this source holds them at
    /// the 55°C setpoint and they lose heat against (55 - 16) = 39 K. The upper two layers were
    /// heated to 60°C by an earlier, higher-setpoint source and do not rise here, so this source
    /// is not holding them at its setpoint: they lose heat against their actual (60 - 16) = 44 K
    /// and coast below 60°C rather than being pulled down to 55°C. The cap is keyed on each
    /// layer's own temperature rise, so it applies to the warmed lower layers but not the
    /// unwarmed upper layers, even though the source adds energy elsewhere in the tank.
    #[rstest]
    fn test_thermal_losses_mixed_warmed_and_unwarmed_layers(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
    ) {
        let (storage_tank1, _) = storage_tank1;

        assert_eq!(
            storage_tank1.calc_temps_after_thermal_losses(
                &[10.0, 10.0, 60.0, 60.0],
                &[55.0, 55.0, 60.0, 60.0],
                &[1.0, 0.0, 0.0, 0.0],
                &[0.1, 0.2, 0.3, 0.4],
                0,
                &[0.0, 0.0, 0.0, 0.0],
                Some(55.)
            ),
            (
                // Q_in_H_W: 1.0 kWh input, no surplus (heater layer 0 is at the setpoint)
                1.0,
                // Total loss: 2 layers at 39 K (0.01762703703703704) + 2 at 44 K (0.019886913580246916)
                0.0750279012345679,
                // Warmed layers settle at 55 - loss/(rho*Cp*Vol) = 54.595...; unwarmed layers
                // coast from 60 at the 44 K loss to 59.543... (not pulled down to the 55 setpoint)
                vec![
                    54.595555555555556,
                    54.595555555555556,
                    59.543703703703706,
                    59.543703703703706
                ],
                vec![
                    0.01762703703703704,
                    0.01762703703703704,
                    0.019886913580246916,
                    0.019886913580246916,
                ]
            )
        );
    }

    /// Q_in_H_W is clamped to zero when a prior source drove the tank above this
    /// source's setpoint, making BS EN 15316-5 Formula (16) negative.
    ///
    /// Layers 2 and 3 are at 65 °C, well above the 55 °C solar-thermal setpoint,
    /// because an earlier heat source with a higher setpoint heated them. Solar thermal
    /// contributes only 0.1 kWh. The cumulative stored energy in layers 2–3 produces
    /// an energy_surplus (~0.96 kWh) far exceeding that 0.1 kWh, so Formula (16) yields
    /// a value of -6978472/8100000. The test reproduces the unclamped
    /// Formula (16) result from the returned Q_ls_n values and asserts it equals that
    /// expected negative value, confirming the scenario genuinely requires the clamp
    /// before asserting Q_in_H_W is 0.0.
    #[rstest]
    fn test_thermal_losses_multi_source_surplus_clamped_to_zero(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
    ) {
        let (storage_tank1, _) = storage_tank1;

        let q_x_in_n = [0.0, 0.0, 0.1, 0.0];
        let q_h_sto_s7 = vec![2.4, 2.4, 2.9, 2.9];
        let heater_layer = 2;
        let temp_setpntmax = 55.0;

        let (q_in_h_w, _, _, q_ls_n) = storage_tank1.calc_temps_after_thermal_losses(
            &[55.0, 55.0, 65.0, 65.0],
            &[55.0, 55.0, 65.0, 65.0],
            &q_x_in_n,
            &q_h_sto_s7,
            heater_layer,
            &[0.; 4],
            Some(temp_setpntmax),
        );

        // Reproduce Formula (16) from BS EN 15316-5 using the returned Q_ls_n values.
        // The unclamped result is exactly -6978472/8100000 ≈ -0.8615397530864198.
        let q_x_in_adj = q_x_in_n.iter().sum::<f64>();
        // Note - python uses fsum here but the test passes with sum (cheaper) so using that for now
        let energy_surplus = (heater_layer..q_h_sto_s7.len())
            .map(|i| {
                q_h_sto_s7[i]
                    - q_ls_n[i]
                    - storage_tank1.rho * storage_tank1.cp * storage_tank1.vol_n[i] * temp_setpntmax
            })
            .sum::<f64>();

        assert_relative_eq!(q_x_in_adj - energy_surplus, -0.8615397530864198);
        assert_eq!(q_in_h_w, 0.);
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
        let heat_source_name = "imheater";
        let temp_s3_n = vec![10., 15., 20., 25., 25., 30., 35., 50.];
        let heat_source = storage_tank1.heat_source_data[heat_source_name]
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
                    heat_source_name.into(),
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
        );
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
        let event = WaterEventResult {
            event_result_type: WaterEventResultType::Other,
            temperature_warm: 41.0,
            volume_warm: 8.0,
            volume_hot: 5.511111111111113,
            event_duration: 0.,
        };

        assert_eq!(
            storage_tank1
                .extract_hot_water(
                    event,
                    simulation_time_for_storage_tank.iter().current_iteration()
                )
                .unwrap(),
            (
                5.51111111111112,
                0.2882311111111112,
                vec![37.5, 37.5, 37.5, 31.988888888888887]
            )
        );
    }

    // Python test test_extract_hot_water_skips_empty_layer skipped as it modifies private state of storage tank instance
    // which would be difficult to replicate here

    #[rstest]
    fn test_calc_temps_ater_extraction(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        let (storage_tank1, _) = storage_tank1;

        let remaining_vols = vec![0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];

        let expected_new_temps = vec![10.0, 10.0, 10.0, 16.0];
        let (actual_new_temps, flag) = storage_tank1
            .calc_temps_after_extraction(
                remaining_vols,
                simulation_time_for_storage_tank.iter().current_iteration(),
            )
            .unwrap();

        assert_eq!(actual_new_temps, expected_new_temps);
        assert_eq!(flag, false);

        let remaining_vols = vec![40.0, 40.0, 40.0, 40.0];

        let expected_new_temps = vec![55.0, 55.0, 55.0, 55.0];
        let (actual_new_temps, flag) = storage_tank1
            .calc_temps_after_extraction(
                remaining_vols,
                simulation_time_for_storage_tank.iter().current_iteration(),
            )
            .unwrap();

        assert_eq!(actual_new_temps, expected_new_temps);
        assert_eq!(flag, false);
    }

    #[rstest]
    fn test_calculate_new_temperatures(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        let (storage_tank1, _) = storage_tank1;
        let remaining_vol = vec![0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4.];

        assert_eq!(
            storage_tank1
                .calc_temps_after_extraction(
                    remaining_vol,
                    simulation_time_for_storage_tank.iter().current_iteration()
                )
                .unwrap(),
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
        let setpnt_diverter = Control::SetpointTime(
            SetpointTimeControl::new(vec![Some(60.); 8], 0, 1.0, None, None, 1.0).into(),
        );
        storage_tank1.q_ls_n_prev_heat_source = Arc::new(RwLock::new(vec![0.0, 0.1, 0.2, 0.3]));
        assert_eq!(
            storage_tank1
                .additional_energy_input(
                    &heat_source.lock(),
                    "imheater",
                    energy_input,
                    Some(&setpnt_diverter),
                    simulation_time_for_storage_tank.iter().current_iteration()
                )
                .unwrap(),
            0.8915535802469137
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
                    Some(&setpnt_diverter),
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

    #[fixture]
    fn storage_tank_with_pipework(
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time_for_storage_tank: SimulationTime,
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions: Arc<ExternalConditions>,
        energy_supply: Arc<RwLock<EnergySupply>>,
    ) -> StorageTank {
        let cold_feed = WaterSupply::ColdWaterSource(cold_water_source.clone());
        let simtime = simulation_time_for_storage_tank.iter().current_iteration();

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

        let energy_supply_connection =
            EnergySupply::connection(energy_supply.clone(), "immersion").unwrap();

        let heat_source_imheater = heat_source(
            simulation_time_for_storage_tank,
            energy_supply_connection.clone(),
            50.0,
            0.1,
            0.33,
            control_min_schedule,
            control_max_schedule,
        );

        let primary_pipework_lst = vec![
            WaterPipework {
                location: WaterPipeworkLocation::Internal,
                internal_diameter_mm: 24.,
                external_diameter_mm: 27.,
                length: 2.,
                insulation_thermal_conductivity: 0.035,
                insulation_thickness_mm: 40.,
                surface_reflectivity: false,
                pipe_contents: PipeworkContents::Water,
            },
            WaterPipework {
                location: WaterPipeworkLocation::External,
                internal_diameter_mm: 25.,
                external_diameter_mm: 27.,
                length: 0.,
                insulation_thermal_conductivity: 0.035,
                insulation_thickness_mm: 38.,
                surface_reflectivity: false,
                pipe_contents: PipeworkContents::Water,
            },
        ];

        let heat_sources = IndexMap::from([("imheater".into(), heat_source_imheater)]);

        StorageTank::new(
            150.0,
            1.68,
            55.0,
            cold_feed,
            &simtime,
            heat_sources,
            temp_internal_air_fn.clone(),
            external_conditions.clone(),
            false,
            Some(4),
            Some(&primary_pipework_lst),
            *WATER,
            None,
            None,
            None,
        )
        .unwrap()
    }

    #[rstest]
    fn test_primary_pipework_losses(
        storage_tank_with_pipework: StorageTank,
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        let input_energy_adj = 0.0;
        let setpnt_max = 55.0;

        for (t_idx, t_it) in simulation_time_for_storage_tank.iter().enumerate() {
            assert_eq!(
                storage_tank_with_pipework
                    .pipework
                    .calculate_primary_pipework_losses(
                        input_energy_adj,
                        setpnt_max.into(),
                        None,
                        &t_it
                    )
                    .unwrap(),
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
            );
        }

        // With value for input_energy_adj
        let input_energy_adj = 3.;

        // First timestep triggers Phase 1 (cool-down) + Phase 2 (steady-state);
        // subsequent timesteps are Phase 2 only because update_tracking records
        // the non-zero energy_input within the mixin.
        let expected_first = (0.04746228058715814, 10.657894331822993);
        let expected_steady = (0.010657894331822992, 10.657894331822993);

        for (t_idx, t_it) in simulation_time_for_storage_tank.iter().enumerate() {
            assert_eq!(
                storage_tank_with_pipework
                    .pipework
                    .calculate_primary_pipework_losses(
                        input_energy_adj,
                        setpnt_max.into(),
                        None,
                        &t_it
                    )
                    .unwrap(),
                if t_idx == 0 {
                    expected_first
                } else {
                    expected_steady
                },
            );
        }
    }

    // Python test test_primary_pipework_losses_end_of_heating skipped as it would required mocking Pipeworkesque::calculate_cool_down_loss

    #[rstest]
    fn test_get_losses_from_primary_pipework_and_storage(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        let (storage_tank1, _) = storage_tank1;
        let usage_event = get_usage_event();
        let _ = storage_tank1.demand_hot_water(
            Some(usage_event),
            simulation_time_for_storage_tank.iter().current_iteration(),
        );

        let expected = (0., 0.07050814814814815);
        let actual = storage_tank1.get_losses_from_primary_pipework_and_storage();

        assert_eq!(actual, expected);
    }

    #[rstest]
    fn test_energy_demand(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        let (storage_tank1, _) = storage_tank1;
        let usage_event = get_usage_event();
        let _ = storage_tank1.demand_hot_water(
            Some(usage_event),
            simulation_time_for_storage_tank.iter().current_iteration(),
        );

        assert_eq!(storage_tank1.test_energy_demand(), 5.9141614814815);
    }

    #[rstest]
    fn test_temperature_and_draw_off_hot_water(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        simulation_time_for_storage_tank: SimulationTime,
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions: Arc<ExternalConditions>,
    ) {
        let cold_water_temps = [60.0, 10.1, 10.2, 10.5, 10.6, 11.0, 11.5, 12.1];
        let cold_feed = WaterSupply::ColdWaterSource(Arc::new(ColdWaterSource::new(
            cold_water_temps.to_vec(),
            0,
            1.,
        )));

        // create storage tank (same as storage_tank1) but with our cold feed
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(
                FuelType::Electricity,
                simulation_time_for_storage_tank.iter().total_steps(),
            )
            .build(),
        ));
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
        let simtime = simulation_time_for_storage_tank.iter().current_iteration();
        let heat_sources = IndexMap::from([("imheater".into(), heat_source)]);
        let storage_tank = StorageTank::new(
            150.0,
            1.68,
            55.0,
            cold_feed,
            &simtime,
            heat_sources,
            temp_internal_air_fn.clone(),
            external_conditions.clone(),
            false,
            None,
            None,
            *WATER,
            None,
            None,
            None,
        )
        .unwrap();

        let expected = (Some(56.875), 240.0);
        let iteration = simulation_time_for_storage_tank.iter().current_iteration();
        let actual = storage_tank.draw_off_hot_water(240., iteration).unwrap();

        assert_eq!(actual, expected);

        // use the default storage tank for the next tests
        // Python re-runs setUp to achieve the same
        let (storage_tank1, _) = storage_tank1;

        let iteration = simulation_time_for_storage_tank.iter().current_iteration();
        assert_eq!(
            (None, 0.0),
            storage_tank1.draw_off_hot_water(0., iteration).unwrap()
        );

        // Note the below are draw_off_water, not draw_off_hot_water
        assert_eq!(
            vec![(55.0, 0.0)],
            storage_tank1.draw_off_water(0., iteration).unwrap()
        );

        assert_eq!(
            vec![(55.0, 22.3)],
            storage_tank1.draw_off_water(22.3, iteration).unwrap()
        );

        assert_eq!(
            vec![
                (55.0, 37.5),
                (55.0, 37.5),
                (55.0, 37.5),
                (28.24, 37.5),
                (10.0, 15.0)
            ],
            storage_tank1.draw_off_water(165., iteration).unwrap()
        );
    }

    #[rstest]
    fn test_heat_source_output(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        storage_tank_with_solar_thermal: (
            StorageTank,
            Arc<Mutex<SolarThermalSystem>>,
            SimulationTime,
            Arc<RwLock<EnergySupply>>,
        ),
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        let (storage_tank1, _) = storage_tank1;
        let heat_source_name = "imheater";
        let positioned_heat_source = storage_tank1.heat_source_data[heat_source_name].clone();
        let heat_source = &*positioned_heat_source.heat_source.lock();

        let iteration = simulation_time_for_storage_tank.iter().current_iteration();
        assert_eq!(
            storage_tank1
                .heat_source_output(
                    heat_source,
                    heat_source_name.into(),
                    43.2,
                    0,
                    iteration,
                    None,
                    None
                )
                .unwrap(),
            43.2
        );

        let (storage_tank_solar_thermal, _, _, _) = storage_tank_with_solar_thermal;
        let heat_source_name = "solthermal";
        let heat_source = storage_tank_solar_thermal.heat_source_data[heat_source_name]
            .clone()
            .heat_source;

        assert_eq!(
            storage_tank1
                .heat_source_output(
                    &heat_source.lock(),
                    "solthermal".into(),
                    43.2,
                    0,
                    iteration,
                    None,
                    None
                )
                .unwrap(),
            0.
        );

        // TODO 1.0.0a9 migration can assertion with heat pump heat source now be replicated?
    }

    /// Test that when hot water demand exceeds tank capacity, remaining volume is drawn from cold feed (e.g. pre-heat tank).
    #[rstest]
    fn test_extract_hot_water_demand_exceeds_tank_capacity(
        storage_tank1: (StorageTank, Arc<RwLock<EnergySupply>>),
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        let (storage_tank1, _) = storage_tank1;
        // Reset draw-off tracking variables
        storage_tank1
            .temp_average_drawoff_volweighted
            .store(0., Ordering::SeqCst);
        storage_tank1
            .total_volume_drawoff
            .store(0., Ordering::SeqCst);

        // Create an event that demands more hot water than the tank can provide
        // Tank volume is 150 litres (4 layers of 37.5 litres each at 55°C)
        // Request 200 litres of hot water - this exceeds tank capacity
        let event = WaterEventResult {
            event_result_type: WaterEventResultType::Bath,
            temperature_warm: 41.,
            volume_warm: 200.,
            volume_hot: 200., // Demand exceeds tank capacity of 150 litres
            event_duration: 0.,
        };

        let (volume_used, energy_withdrawn, remaining_vols) = storage_tank1
            .extract_hot_water(
                event,
                simulation_time_for_storage_tank.iter().current_iteration(),
            )
            .unwrap();

        // Volume used from tank should be entire tank capacity
        assert_relative_eq!(volume_used, 150.);

        // All tank layers should be depleted
        for vol in remaining_vols.iter() {
            assert_relative_eq!(*vol, 0.);
        }

        // Total draw-off should include both tank water and cold feed water
        // 150 litres from tank + 50 litres from cold feed = 200 litres total
        assert_relative_eq!(
            storage_tank1.total_volume_drawoff.load(Ordering::SeqCst),
            200.
        );

        // Energy withdrawn should be positive (hot water from tank + any pre-heated water)
        assert_relative_eq!(energy_withdrawn, 7.845);
    }

    #[fixture]
    fn storage_tank3(
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time_for_storage_tank: SimulationTime,
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions: Arc<ExternalConditions>,
        energy_supply: Arc<RwLock<EnergySupply>>,
    ) -> StorageTank {
        let cold_feed = WaterSupply::ColdWaterSource(cold_water_source.clone());
        let simtime = simulation_time_for_storage_tank.iter().current_iteration();

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

        let energy_supply_connection =
            EnergySupply::connection(energy_supply.clone(), "immersion3").unwrap();

        let imheater3 = heat_source(
            simulation_time_for_storage_tank,
            energy_supply_connection.clone(),
            50.0,
            0.1,
            0.33,
            control_min_schedule,
            control_max_schedule,
        );

        let heat_sources = IndexMap::from([("imheater3".into(), imheater3)]);
        StorageTank::new(
            210.0,
            1.61,
            52.0,
            cold_feed,
            &simtime,
            heat_sources,
            temp_internal_air_fn.clone(),
            external_conditions.clone(),
            false,
            None,
            None,
            *WATER,
            None,
            None,
            None,
        )
        .unwrap()
    }

    /// Test that when demand exceeds tank capacity, water is properly drawn from a pre-heat tank configured as the cold feed.
    #[rstest]
    fn test_extract_hot_water_demand_exceeds_capacity_with_preheat_tank(
        storage_tank3: StorageTank,
        simulation_time_for_storage_tank: SimulationTime,
    ) {
        // Use storagetank3 which has a pre-heated storage tank as its cold feed
        // storagetank3: 210 litres, init_temp=52°C
        // preheatfeed: 80 litres, init_temp=30°C
        storage_tank3
            .temp_average_drawoff_volweighted
            .store(0., Ordering::SeqCst);
        storage_tank3
            .total_volume_drawoff
            .store(0., Ordering::SeqCst);

        // Create an event that demands more than storagetank3 capacity (210 litres)
        // but less than storagetank3 + preheatfeed combined (290 litres)
        let event = WaterEventResult {
            event_result_type: WaterEventResultType::Bath,
            temperature_warm: 41.,
            volume_warm: 250.,
            volume_hot: 250., // Exceeds 210 litre tank, needs 40 litres from pre-heat
            event_duration: 0.,
        };

        let (volume_used, energy_withdrawn, _) = storage_tank3
            .extract_hot_water(
                event,
                simulation_time_for_storage_tank.iter().current_iteration(),
            )
            .unwrap();

        // Volume used from main tank should be entire tank capacity
        assert_relative_eq!(volume_used, 210.);

        // Total draw-off should include water from pre-heat tank
        // 210 litres from main tank + 40 litres from pre-heat tank = 250 litres
        assert_relative_eq!(
            storage_tank3.total_volume_drawoff.load(Ordering::SeqCst),
            250.
        );

        // Verify energy calculation accounts for pre-heated water
        // Energy should be greater than zero since we're drawing hot/warm water
        assert!(energy_withdrawn > 0.);
    }

    // Skipping Python's test_primary_pipework_losses_between_events due to mocking of calculate_cool_down_loss return value

    /// A source left active by a previous timestep is switched off on entering
    /// an off period, before any charging.
    ///
    /// When the schedule has no setpoint (off period) _determine_heat_source_switch_on
    /// deactivates a source that was active in the previous timestep, before the
    /// charging block runs. Without this the source would charge for one extra
    /// timestep across the on-to-off transition. _heating_active is the method's
    /// only output, so it is asserted directly.
    #[rstest]
    fn test_determine_heat_source_switch_on_off_period_deactivates(
        cold_water_source: Arc<ColdWaterSource>,
        simulation_time_for_storage_tank: SimulationTime,
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions: Arc<ExternalConditions>,
        energy_supply: Arc<RwLock<EnergySupply>>,
    ) {
        // Controls with no setpoint at every timestep represent an off period
        let cold_feed = WaterSupply::ColdWaterSource(cold_water_source.clone());
        let simtime = simulation_time_for_storage_tank.iter().current_iteration();

        let control_off = vec![None; 8];
        let heat_source_name = "immersion_off";
        let energy_supply_connection =
            EnergySupply::connection(energy_supply.clone(), heat_source_name).unwrap();

        let heat_source_imheater = heat_source(
            simulation_time_for_storage_tank,
            energy_supply_connection.clone(),
            50.0,
            0.1,
            0.33,
            control_off.clone(),
            control_off,
        );

        let heat_sources =
            IndexMap::from([(heat_source_name.into(), heat_source_imheater.clone())]);

        let tank = StorageTank::new(
            150.0,
            1.68,
            55.0,
            cold_feed,
            &simtime,
            heat_sources,
            temp_internal_air_fn.clone(),
            external_conditions.clone(),
            false,
            Some(4),
            None,
            *WATER,
            None,
            None,
            None,
        )
        .unwrap();

        // Source was left active by the previous (on) timestep
        tank.heating_active[heat_source_name].store(true, Ordering::SeqCst);

        tank.determine_heat_source_switch_on(
            &vec![55.; tank.number_of_volumes],
            heat_source_name,
            &heat_source_imheater.heat_source.lock(),
            (0.1 * tank.number_of_volumes as f64) as usize,
            (0.33 * tank.number_of_volumes as f64) as usize,
            simulation_time_for_storage_tank.iter().current_iteration(),
        )
        .unwrap();

        assert!(!tank.heating_active[heat_source_name].load(Ordering::SeqCst));
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
            simulation_time_for_immersion_heater.iter().total_steps(),
        )
        .build();
        let energy_supply_connection =
            EnergySupply::connection(Arc::new(RwLock::new(energy_supply)), "shower").unwrap();
        let timestep = simulation_time_for_immersion_heater.step;

        let control_min = ScheduleOrControl::Schedule(vec![Some(52.), Some(52.), None, Some(52.)]);
        let control_max =
            ScheduleOrControl::Schedule(vec![Some(60.), Some(60.), Some(60.), Some(60.)]);

        let control = Some(Arc::from(
            RangeTimeControl::new(
                control_min,
                control_max,
                simulation_time_for_immersion_heater.iter(),
                0,
                1.,
                None,
            )
            .unwrap(),
        ));

        ImmersionHeater::new(rated_power, energy_supply_connection, timestep, control)
    }

    #[ignore = "Update as part of migration 1.0.0a9"]
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
                    .demand_energy(energy_inputs[t_idx], None, t_it)
                    .unwrap(),
                expected_energy[t_idx],
                "incorrect energy demand calculated"
            );
        }
        assert!(immersion_heater
            .demand_energy(
                -1.,
                None,
                simulation_time_for_immersion_heater
                    .iter()
                    .current_iteration()
            )
            .is_err());
    }

    #[ignore = "Update as part of migration 1.0.0a9"]
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
                simulation_time_for_smart_hot_water_tank
                    .iter()
                    .total_steps(),
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
                simulation_time_for_smart_hot_water_tank
                    .iter()
                    .total_steps(),
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

    static COLD_WATER_TEMPS: [f64; 8] = [10.0, 10.1, 10.2, 10.5, 10.6, 11.0, 11.5, 12.1];

    fn create_smart_hot_water_tank_with_defaults(
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
        let temp_setpnt_max = Control::SetpointTime(
            SetpointTimeControl::new(
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
                None,
                None,
                1.,
            )
            .into(),
        );

        create_smart_hot_water_tank(
            simulation_time_for_smart_hot_water_tank,
            temp_internal_air_fn,
            external_conditions_for_smart_hot_water_tank,
            energy_supply_for_smart_hot_water_tank_immersion,
            energy_supply_for_smart_hot_water_tank_pump,
            heat_source_name,
            volume,
            losses,
            init_temp,
            power_pump_kw,
            max_flow_rate_pump_l_per_min,
            temp_usable,
            temp_setpnt_max,
            4,
        )
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
        power_pump_kw: f64,
        max_flow_rate_pump_l_per_min: f64,
        temp_usable: f64,
        temp_setpnt_max: Control,
        nb_vol: usize,
    ) -> SmartHotWaterTank {
        let cold_feed = WaterSupply::ColdWaterSource(Arc::new(ColdWaterSource::new(
            COLD_WATER_TEMPS.to_vec(),
            0,
            1.,
        )));
        let energy_supply_connection = EnergySupply::connection(
            energy_supply_for_smart_hot_water_tank_immersion.clone(),
            heat_source_name,
        )
        .unwrap();

        let control_min = ScheduleOrControl::Schedule(vec![
            Some(0.5),
            None,
            None,
            None,
            Some(0.5),
            Some(0.5),
            Some(0.5),
            Some(0.5),
        ]);

        let control_max = ScheduleOrControl::Schedule(vec![
            Some(1.0),
            Some(1.0),
            Some(0.9),
            Some(0.8),
            Some(0.7),
            Some(1.0),
            Some(0.9),
            Some(0.8),
        ]);

        let control = Some(Arc::from(
            RangeTimeControl::new(
                control_min,
                control_max,
                simulation_time_for_smart_hot_water_tank.iter(),
                0,
                1.,
                None,
            )
            .unwrap(),
        ));
        let immersion_heater = ImmersionHeater::new(5., energy_supply_connection, 1., control);
        let heat_source = HeatSource::Storage(HeatSourceWithStorageTank::Immersion(Arc::new(
            Mutex::new(immersion_heater),
        )));
        let heat_sources = IndexMap::from([(
            heat_source_name.into(),
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
                .current_iteration(),
            heat_sources,
            temp_internal_air_fn,
            external_conditions_for_smart_hot_water_tank,
            None,
            Some(nb_vol),
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
        create_smart_hot_water_tank_with_defaults(
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

    fn get_event_data_immersion() -> Vec<Option<Vec<WaterEventResult>>> {
        vec![
            Some(vec![
                WaterEventResult {
                    event_result_type: WaterEventResultType::Shower,
                    temperature_warm: 41.0,
                    volume_warm: 48.0,
                    volume_hot: 33.0666666666667,
                    event_duration: 0.,
                },
                WaterEventResult {
                    event_result_type: WaterEventResultType::Bath,
                    temperature_warm: 43.0,
                    volume_warm: 100.0,
                    volume_hot: 73.3333333333333,
                    event_duration: 0.,
                },
                WaterEventResult {
                    event_result_type: WaterEventResultType::Other,
                    temperature_warm: 40.0,
                    volume_warm: 8.0,
                    volume_hot: 5.3333333333333,
                    event_duration: 0.,
                },
            ]),
            Some(vec![WaterEventResult {
                event_result_type: WaterEventResultType::Shower,
                temperature_warm: 41.0,
                volume_warm: 48.0,
                volume_hot: 33.0334075723831,
                event_duration: 0.,
            }]),
            None,
            Some(vec![WaterEventResult {
                event_result_type: WaterEventResultType::Shower,
                temperature_warm: 45.0,
                volume_warm: 48.0,
                volume_hot: 37.8988082756996,
                event_duration: 0.,
            }]),
            None,
            Some(vec![WaterEventResult {
                event_result_type: WaterEventResultType::Shower,
                temperature_warm: 41.0,
                volume_warm: 52.0,
                volume_hot: 35.4545454545455,
                event_duration: 0.,
            }]),
            None,
            None,
        ]
    }

    fn get_usage_event() -> Vec<WaterEventResult> {
        get_event_data_immersion().first().unwrap().clone().unwrap()
    }

    fn get_event_data_solthermal() -> Vec<Option<Vec<WaterEventResult>>> {
        vec![
            None,
            Some(vec![WaterEventResult {
                event_result_type: WaterEventResultType::Shower,
                temperature_warm: 41.0,
                volume_warm: 48.0,
                volume_hot: 30.5956261482843,
                event_duration: 0.,
            }]),
            None,
            Some(vec![WaterEventResult {
                event_result_type: WaterEventResultType::Shower,
                temperature_warm: 45.0,
                volume_warm: 48.0,
                volume_hot: 36.4281898110265,
                event_duration: 0.,
            }]),
            None,
            Some(vec![WaterEventResult {
                event_result_type: WaterEventResultType::Shower,
                temperature_warm: 41.0,
                volume_warm: 52.0,
                volume_hot: 34.4038433055010,
                event_duration: 0.,
            }]),
            None,
            None,
            None,
            Some(vec![WaterEventResult {
                event_result_type: WaterEventResultType::Shower,
                temperature_warm: 41.0,
                volume_warm: 48.0,
                volume_hot: 33.3416695316938,
                event_duration: 0.,
            }]),
            None,
            Some(vec![]),
            None,
            Some(vec![WaterEventResult {
                event_result_type: WaterEventResultType::Shower,
                temperature_warm: 41.0,
                volume_warm: 52.0,
                volume_hot: 40.521971319747124,
                event_duration: 0.,
            }]),
            None,
            None,
            None,
            Some(vec![WaterEventResult {
                event_result_type: WaterEventResultType::Shower,
                temperature_warm: 41.0,
                volume_warm: 48.0,
                volume_hot: 30.5956261482843,
                event_duration: 0.,
            }]),
            None,
            Some(vec![WaterEventResult {
                event_result_type: WaterEventResultType::Shower,
                temperature_warm: 45.0,
                volume_warm: 48.0,
                volume_hot: 36.42818981102645,
                event_duration: 0.,
            }]),
            None,
            Some(vec![WaterEventResult {
                event_result_type: WaterEventResultType::Shower,
                temperature_warm: 41.0,
                volume_warm: 52.0,
                volume_hot: 34.40384330550096,
                event_duration: 0.,
            }]),
            None,
            None,
        ]
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
        assert_relative_eq!(soc.unwrap(), 0.850, max_relative = 1e-3);
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

    #[ignore = "Update as part of migration 1.0.0a9"]
    #[rstest]
    fn test_demand_hot_water_for_smart_hot_water_tank(
        simulation_time_for_smart_hot_water_tank: SimulationTime,
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions_for_smart_hot_water_tank: Arc<ExternalConditions>,
    ) {
        // this is TestSmartHotWaterTank.test_demand_hot_water in Python
        let energy_supply_for_smart_hot_water_tank_immersion_1 = Arc::from(RwLock::from(
            EnergySupplyBuilder::new(
                FuelType::Electricity,
                simulation_time_for_smart_hot_water_tank
                    .iter()
                    .total_steps(),
            )
            .build(),
        ));
        let energy_supply_for_smart_hot_water_tank_pump_1 = Arc::from(RwLock::from(
            EnergySupplyBuilder::new(
                FuelType::Electricity,
                simulation_time_for_smart_hot_water_tank
                    .iter()
                    .total_steps(),
            )
            .build(),
        ));

        let energy_supply_for_smart_hot_water_tank_immersion_2 = Arc::from(RwLock::from(
            EnergySupplyBuilder::new(
                FuelType::Electricity,
                simulation_time_for_smart_hot_water_tank
                    .iter()
                    .total_steps(),
            )
            .build(),
        ));
        let energy_supply_for_smart_hot_water_tank_pump_2 = Arc::from(RwLock::from(
            EnergySupplyBuilder::new(
                FuelType::Electricity,
                simulation_time_for_smart_hot_water_tank
                    .iter()
                    .total_steps(),
            )
            .build(),
        ));

        let smart_hot_water_tank = create_smart_hot_water_tank_with_defaults(
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

        let smart_hot_water_tank_2 = create_smart_hot_water_tank_with_defaults(
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

        let event_data = get_event_data_immersion();

        let expected_temperatures_1 = &[
            vec![42.06412979639594, 50.0, 50.0, 50.0],
            vec![
                26.168457194735883,
                45.942228008024884,
                49.87555555555556,
                49.87555555555556,
            ],
            vec![
                26.115731861133547,
                45.86963541543229,
                49.802962962962965,
                49.802962962962965,
            ],
            vec![
                17.336010881790145,
                34.751355008536194,
                47.57251929648306,
                49.782222222222224,
            ],
            vec![31.875588805820787, 50.0, 50.0, 50.0],
            vec![
                20.717353598198578,
                40.20747289529573,
                49.82370370370371,
                49.82370370370371,
            ],
            vec![
                20.69289324620792,
                40.08195266546827,
                49.64832153635117,
                49.64832153635117,
            ],
            vec![
                20.668559725672026,
                39.95708328127696,
                49.47384875801453,
                49.47384875801453,
            ],
        ];

        let expected_temperatures_2 = [
            vec![
                10.0,
                24.55607367670878,
                50.03631427851564,
                59.092295619623506,
            ],
            vec![
                10.057665043481068,
                16.16115527929868,
                35.205810167360966,
                53.69978967126601,
            ],
            vec![
                10.057665043481068,
                16.160011275772792,
                35.10642745131158,
                53.60040695521663,
            ],
            vec![
                10.381386078348926,
                11.69403406025759,
                21.212174941529444,
                40.037268379332176,
            ],
            vec![
                11.520424219588156,
                48.63806146445693,
                48.63806146445693,
                48.63806146445693,
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
            [2.2101151057, 0.0, 0.0, 0.0, 2.0951108105, 0.0, 0.0, 0.0];

        let expected_results_by_end_user_2 =
            [0.0, 0.0, 0.0, 0.0, 4.5043354264, 0.1956556236, 0.0, 0.0];

        for (t_idx, t_it) in simulation_time_for_smart_hot_water_tank.iter().enumerate() {
            // # Convert usage events based on HW temp of 55 to equivalent 50:
            let usage_events = event_data[t_idx].clone();
            let temp_hot = if t_idx == 0 {
                50.
            } else {
                *expected_temperatures_1[t_idx - 1].last().unwrap()
            };

            let mut usage_events1 = vec![];
            if usage_events.is_some() {
                for event in usage_events.clone().unwrap() {
                    let volume_hot = event.volume_warm
                        * (event.temperature_warm - COLD_WATER_TEMPS[t_idx])
                        / (temp_hot - COLD_WATER_TEMPS[t_idx]);
                    usage_events1.push(WaterEventResult {
                        event_result_type: event.event_result_type,
                        temperature_warm: event.temperature_warm,
                        volume_warm: event.volume_warm,
                        volume_hot,
                        event_duration: 0.,
                    });
                }
            }

            let _ = smart_hot_water_tank
                .demand_hot_water(Some(usage_events1), t_it)
                .unwrap();

            let temp_n = smart_hot_water_tank.storage_tank.temp_n.read();

            for (i, expected_temp) in expected_temperatures_1[t_idx].iter().enumerate() {
                assert_relative_eq!(temp_n[i], *expected_temp, max_relative = 1e-9);
            }

            let results_by_end_user = energy_supply_for_smart_hot_water_tank_immersion_1
                .read()
                .results_by_end_user();
            let actual_results_by_end_user_1 = results_by_end_user.get("imheater").unwrap();
            assert_relative_eq!(
                actual_results_by_end_user_1[t_idx],
                expected_results_by_end_user_1[t_idx],
                max_relative = FIVE_DECIMAL_PLACES
            );

            // smart_hot_water_tank_2 tests for case where heater does not heat all layers

            // Convert usage events based on HW temp of 55 to equivalent 60:
            let temp_hot = if t_idx == 0 {
                60.
            } else {
                *expected_temperatures_2[t_idx - 1].last().unwrap()
            };
            let mut usage_events2 = vec![];

            if usage_events.is_some() {
                for event in usage_events.unwrap() {
                    let volume_hot = event.volume_warm
                        * (event.temperature_warm - COLD_WATER_TEMPS[t_idx])
                        / (temp_hot - COLD_WATER_TEMPS[t_idx]);
                    usage_events2.push(WaterEventResult {
                        event_result_type: event.event_result_type,
                        temperature_warm: event.temperature_warm,
                        volume_warm: event.volume_warm,
                        volume_hot,
                        event_duration: 0.,
                    });
                }
            }
            let _ = smart_hot_water_tank_2
                .demand_hot_water(Some(usage_events2.clone()), t_it)
                .unwrap();

            let temp_n = smart_hot_water_tank_2.storage_tank.temp_n.read();
            for (i, expected_temp) in expected_temperatures_2[t_idx].iter().enumerate() {
                assert_relative_eq!(
                    temp_n[i],
                    *expected_temp,
                    max_relative = FIVE_DECIMAL_PLACES
                );
            }

            let results_by_end_user_2 = energy_supply_for_smart_hot_water_tank_immersion_2
                .read()
                .results_by_end_user();
            let actual_results_by_end_user_2 = results_by_end_user_2.get("immersion2").unwrap();
            assert_relative_eq!(
                actual_results_by_end_user_2[t_idx],
                expected_results_by_end_user_2[t_idx],
                max_relative = FIVE_DECIMAL_PLACES
            );
        }
    }

    // Python test_demand_hot_water_edge_cases skipped

    #[rstest]
    fn test_calc_final_temps(
        mut smart_hot_water_tank: SmartHotWaterTank,
        simulation_time_iteration_for_smart_hot_water_tank: SimulationTimeIteration,
    ) {
        let temp_s3_n = vec![10.0, 15.0, 20.0, 25.0, 25.0, 30.0, 35.0, 50.0];
        let q_x_in_n = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let heater_layer = 2;
        let q_ls_n_prev_heat_source = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

        let heat_source_name = "imheater";
        let positioned_heat_source =
            smart_hot_water_tank.storage_tank.heat_source_data[heat_source_name].clone();
        let heat_source = &*positioned_heat_source.heat_source.lock();

        let control_max_diverter = Control::SetpointTime(
            SetpointTimeControl::new(
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
                None,
                None,
                0.,
            )
            .into(),
        );

        let expected = TemperatureCalculation {
            temp_s8_n: vec![50.0, 50.0, 50.0, 50.0],
            q_x_in_n: q_x_in_n.clone(),
            q_s6: 42.10166666666667,
            temp_s6_n: vec![10.0, 15.0, 433.00191204588907, 25.0],
            temp_s7_n: vec![
                229.00095602294454,
                229.00095602294454,
                229.00095602294454,
                229.00095602294454,
            ],
            q_in_h_w: 4.824900987654324,
            q_ls: 0.061468641975308644,
            q_ls_n: vec![
                0.015367160493827161,
                0.015367160493827161,
                0.015367160493827161,
                0.015367160493827161,
            ],
        };

        let actual = smart_hot_water_tank
            .calc_final_temps(
                &temp_s3_n,
                heat_source,
                heat_source_name.into(),
                q_x_in_n.clone(),
                heater_layer,
                &q_ls_n_prev_heat_source,
                Some(&control_max_diverter),
                simulation_time_iteration_for_smart_hot_water_tank,
            )
            .unwrap();

        assert_eq!(actual, expected);

        smart_hot_water_tank.temp_usable = 100.0;
        let control_max_diverter = Control::SetpointTime(
            SetpointTimeControl::new(
                vec![
                    Some(0.0),
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
                None,
                None,
                0.,
            )
            .into(),
        );

        // NOTE - these are the same expected values as above. Same behaviour in Python
        let expected = TemperatureCalculation {
            temp_s8_n: vec![10.0, 15.0, 19.97925925925926, 24.953333333333333],
            q_x_in_n: q_x_in_n.clone(),
            q_s6: 6.101666666666667,
            temp_s6_n: vec![10.0, 15.0, 20.0, 25.0],
            temp_s7_n: vec![10.0, 15.0, 20.0, 25.0],
            q_in_h_w: 0.,
            q_ls: 0.005875679012345679,
            q_ls_n: vec![0.0, 0.0, 0.0018079012345679013, 0.004067777777777778],
        };

        let actual = smart_hot_water_tank
            .calc_final_temps(
                &temp_s3_n,
                heat_source,
                heat_source_name.into(),
                q_x_in_n,
                heater_layer,
                &q_ls_n_prev_heat_source,
                Some(&control_max_diverter),
                simulation_time_iteration_for_smart_hot_water_tank,
            )
            .unwrap();
        assert_eq!(actual, expected);
    }

    #[rstest]
    fn test_temps_after_pumping(smart_hot_water_tank: SmartHotWaterTank) {
        let mut volumes = vec![120.0, 37.5, 37.5, 37.5];
        let guard = smart_hot_water_tank.storage_tank.temp_n.read();
        let tank_layer_temperatures = guard.as_slice();

        let expected = vec![50.0, 50.0, 50.0, 50.0];
        let actual = smart_hot_water_tank
            .temps_after_pumping(10., &mut volumes, tank_layer_temperatures)
            .unwrap();

        assert_eq!(actual, expected);
    }

    #[rstest]
    fn test_bottom_to_top_pump_volume_none_setpoint(
        simulation_time_for_smart_hot_water_tank: SimulationTime,
        temp_internal_air_fn: TempInternalAirFn,
        external_conditions_for_smart_hot_water_tank: Arc<ExternalConditions>,
        energy_supply_for_smart_hot_water_tank_immersion: Arc<RwLock<EnergySupply>>,
        energy_supply_for_smart_hot_water_tank_pump: Arc<RwLock<EnergySupply>>,
    ) {
        let temp_setpnt_max = Control::Mock(MockControl::new(None, None, None));

        let tank_with_none_setpoint = create_smart_hot_water_tank(
            simulation_time_for_smart_hot_water_tank,
            temp_internal_air_fn,
            external_conditions_for_smart_hot_water_tank,
            energy_supply_for_smart_hot_water_tank_immersion,
            energy_supply_for_smart_hot_water_tank_pump,
            "imheater",
            80.,
            1.0,
            10.,
            0.1,
            10.0,
            40.0,
            temp_setpnt_max,
            8,
        );

        let temp_s7_n = &[60.0, 20.0, 25.0, 30.0, 35.0, 35.0, 35.0, 55.0];
        let qin = 1.0;
        let heater_layer = 7;
        let volumes = &[10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0];

        // This should use temp_usable (40.0) instead of the None setpoint
        let volume_pumped = tank_with_none_setpoint
            .bottom_to_top_pump_volume(
                temp_s7_n,
                qin,
                heater_layer,
                volumes,
                simulation_time_for_smart_hot_water_tank
                    .iter()
                    .current_iteration(),
            )
            .unwrap();

        assert!(volume_pumped >= 0.0);

        // NOTE in Python they assert that the mocked setpoint is called
        // But we can't replicate that easily in Rust
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
            )
            .unwrap();
        assert_eq!(
            actual_remaining_vol,
            (vec![10.0, 10.0, 10.0, 12.666666666666666], false)
        );
    }

    #[rstest]
    fn test_demand_hot_water_solthermal(
        storage_tank_with_solar_thermal: (
            StorageTank,
            Arc<Mutex<SolarThermalSystem>>,
            SimulationTime,
            Arc<RwLock<EnergySupply>>,
        ),
    ) {
        let (storage_tank_solar_thermal, _, simulation_time, _) = storage_tank_with_solar_thermal;
        let event_data = get_event_data_solthermal();

        let expected_energy_demand = [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.39440131536356493,
            0.8431945125549533,
            1.3874298880308749,
            1.092014226211686,
            1.1503560996860809,
            1.484510483919223,
            0.9003607869563452,
            0.4981024012117776,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let usage_events = event_data[t_idx].clone();
            let _ = storage_tank_solar_thermal.demand_hot_water(usage_events, t_it);

            let actual = storage_tank_solar_thermal
                .energy_demand_test
                .load(Ordering::SeqCst);

            assert_relative_eq!(actual, expected_energy_demand[t_idx], max_relative = 1e-7);
        }
    }

    #[rstest]
    fn test_capacity_used(
        mut storage_tank_for_pv_diverter: StorageTank,
        immersion_heater: ImmersionHeater,
        diverter_control: Control,
    ) {
        storage_tank_for_pv_diverter.q_ls_n_prev_heat_source =
            Arc::new(RwLock::new(vec![0.0, 0.1, 0.2, 0.3]));
        let pvdiverter = PVDiverter::new(
            &HotWaterStorageTank::StorageTank(Arc::new(RwLock::new(storage_tank_for_pv_diverter))),
            Arc::new(Mutex::new(immersion_heater)),
            "imheater".into(),
            diverter_control.into(),
        );

        let capacity_used = 2.3;
        pvdiverter.read().increment_capacity_used(capacity_used);

        let actual: f64 = pvdiverter.read().capacity_used.load(Ordering::SeqCst);
        assert_eq!(actual, capacity_used);
    }

    #[rstest]
    fn test_timestep_end(
        mut storage_tank_for_pv_diverter: StorageTank,
        immersion_heater: ImmersionHeater,
        diverter_control: Control,
    ) {
        storage_tank_for_pv_diverter.q_ls_n_prev_heat_source =
            Arc::new(RwLock::new(vec![0.0, 0.1, 0.2, 0.3]));
        let pvdiverter = PVDiverter::new(
            &HotWaterStorageTank::StorageTank(Arc::new(RwLock::new(storage_tank_for_pv_diverter))),
            Arc::new(Mutex::new(immersion_heater)),
            "imheater".into(),
            diverter_control.into(),
        );

        let capacity_used = 2.3;
        pvdiverter.read().increment_capacity_used(capacity_used);

        pvdiverter.read().timestep_end();

        let actual: f64 = pvdiverter.read().capacity_used.load(Ordering::SeqCst);
        assert_eq!(actual, 0.);
    }

    #[rstest]
    fn test_divert_surplus(
        mut storage_tank_for_pv_diverter: StorageTank,
        immersion_heater: ImmersionHeater,
        diverter_control: Control,
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

    // Python test test_energy_potential and test_energy_supply skipped as they only test initialisation
}
