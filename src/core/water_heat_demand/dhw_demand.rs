use crate::core::common::WaterSupplyBehaviour;
use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::core::heating_systems::wwhrs::Wwhrs;
use crate::core::pipework::{PipeworkLocation, PipeworkSimple, Pipeworkesque};
use crate::core::schedule::{TypedScheduleEvent, WaterScheduleEventType};
use crate::core::units::{MILLIMETRES_IN_METRE, WATTS_PER_KILOWATT};
use crate::core::water_heat_demand::bath::Bath;
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::misc::{
    water_demand_to_kwh, WaterEventResult, WaterEventResultType, FRAC_DHW_ENERGY_INTERNAL_GAINS,
};
use crate::core::water_heat_demand::other_hot_water_uses::OtherHotWater;
use crate::core::water_heat_demand::shower::Shower;
use crate::core::water_heat_demand::shower::{InstantElectricShower, MixerShower};
use crate::corpus::{ColdWaterSources, EventSchedule, HotWaterSourceBehaviour};
use crate::input::{
    BathDetails, Baths as BathInput, OtherWaterUse, OtherWaterUses as OtherWaterUseInput,
    PipeworkContents, Shower as ShowerInput, Showers as ShowersInput,
    WaterDistribution as WaterDistributionInput, WaterHeatingEvent, WaterPipeworkSimple,
};
use crate::simulation_time::SimulationTimeIteration;
use anyhow::{anyhow, bail};
use indexmap::IndexMap;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use parking_lot::{Mutex, RwLock};
use smartstring::alias::String;
use std::collections::HashMap;
use std::sync::Arc;

const ELECTRIC_SHOWERS_HWS_NAME: &str = "_electric_showers";

#[derive(Debug)]
pub(crate) struct DomesticHotWaterDemand<T: HotWaterSourceBehaviour> {
    showers: HashMap<String, Shower>,
    baths: HashMap<String, Bath>,
    other: HashMap<String, OtherHotWater>,
    hot_water_sources: IndexMap<String, T>,
    energy_supply_conn_unmet_demand: IndexMap<String, EnergySupplyConnection>,
    source_supplying_outlet: HashMap<(OutletType, String), String>,
    hot_water_distribution_pipework: IndexMap<String, Vec<PipeworkSimple>>,
    event_schedules: EventSchedule,
    pre_heated_water_sources: IndexMap<String, T>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VolumeReference {
    pub warm_temp: f64,
    pub warm_vol: f64,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum DemandVolTargetKey {
    Number(OrderedFloat<f64>),
    TempHotWater,
}

impl From<f64> for DemandVolTargetKey {
    fn from(value: f64) -> Self {
        Self::Number(OrderedFloat(value))
    }
}

// Note this type is not in Python
pub(crate) enum TappingPoint<'a> {
    Shower(&'a Shower),
    Bath(&'a Bath),
    Other(&'a OtherHotWater),
}

impl TappingPoint<'_> {
    pub fn hot_water_demand<'a>(
        &'a self,
        event: WaterHeatingEvent,
        func_temp_hot_water: &'a Box<dyn Fn(f64) -> f64 + 'a>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(Option<f64>, f64)> {
        match self {
            TappingPoint::Shower(shower) => {
                shower.hot_water_demand(event, func_temp_hot_water, simtime)
            }
            TappingPoint::Bath(bath) => bath.hot_water_demand(event, func_temp_hot_water, simtime),
            TappingPoint::Other(other_hot_water) => {
                other_hot_water.hot_water_demand(&event, func_temp_hot_water, simtime)
            }
        }
    }

    pub fn get_cold_water_source(&self) -> &ColdWaterSource {
        match self {
            TappingPoint::Shower(shower) => shower.get_cold_water_source(),
            TappingPoint::Bath(bath) => bath.get_cold_water_source(),
            TappingPoint::Other(other_hot_water) => other_hot_water.get_cold_water_source(),
        }
    }
}

#[derive(Eq, Hash, PartialEq, Debug)]
pub enum OutletType {
    Shower,
    Bath,
    Other,
}

impl<T: HotWaterSourceBehaviour> DomesticHotWaterDemand<T> {
    // TODO (from Python) Enhance analysis for overlapping events
    // Part of draft code for future overlapping analysis of events
    // For pipework losses count only none overlapping events
    // Time of finalisation of the previous hot water event
    // __time_end_previous_event = 0.0

    pub(crate) fn new(
        showers_input: ShowersInput,
        bath_input: BathInput,
        other_hot_water_input: OtherWaterUseInput,
        hw_pipework_inputs: WaterDistributionInput,
        cold_water_sources: &ColdWaterSources,
        wwhrs: &IndexMap<String, Arc<Mutex<Wwhrs>>>,
        energy_supplies: &IndexMap<String, Arc<RwLock<EnergySupply>>>,
        event_schedules: EventSchedule,
        hot_water_sources: IndexMap<String, T>,
        pre_heated_water_sources: IndexMap<String, T>,
    ) -> anyhow::Result<Self> {
        let showers: HashMap<String, Shower> = showers_input
            .0
            .iter()
            .map(|(name, shower)| {
                Ok((
                    name.into(),
                    shower_from_input(name, shower, cold_water_sources, energy_supplies, wwhrs)?,
                ))
            })
            .collect::<anyhow::Result<HashMap<_, _>>>()?;
        let baths: HashMap<String, Bath> = bath_input
            .0
            .iter()
            .map(|(name, bath)| (name.into(), input_to_bath(bath, cold_water_sources)))
            .collect();
        let other: HashMap<String, OtherHotWater> = other_hot_water_input
            .0
            .iter()
            .map(|(name, other)| {
                (
                    name.into(),
                    input_to_other_water_events(other, cold_water_sources),
                )
            })
            .collect();
        let mixer_shower_count = showers
            .iter()
            .filter(|(_, shower)| !matches!(shower, Shower::InstantElectricShower(_)))
            .count();
        let total_number_tapping_points = mixer_shower_count + baths.len() + other.len();

        let mut hot_water_distribution_pipework: IndexMap<String, Vec<PipeworkSimple>> =
            hot_water_sources
                .keys()
                .map(|key| -> (String, Vec<PipeworkSimple>) { (key.clone(), vec![]) })
                .collect();

        // if we have a list (and only one heat source) convert it into map
        let mut hw_pipework_inputs: IndexMap<String, Vec<WaterPipeworkSimple>> =
            match hw_pipework_inputs {
                WaterDistributionInput::List(pipeworks) => {
                    if hot_water_sources.len() == 1 {
                        hot_water_sources
                            .keys()
                            .map(|key| -> (String, Vec<WaterPipeworkSimple>) {
                                (key.clone(), pipeworks.clone())
                            })
                            .collect()
                    } else {
                        bail!("If more than one HotWaterSource is defined, then distribution pipework must be defined for each one");
                    }
                }
                WaterDistributionInput::Map(index_map) => index_map
                    .iter()
                    .map(|(key, value)| -> (String, Vec<WaterPipeworkSimple>) {
                        (key.into(), value.clone())
                    })
                    .collect(),
            };

        // pipework without a valid hot water source
        let pws_without_hws: Vec<_> = hw_pipework_inputs
            .keys()
            .filter(|key| !hot_water_sources.keys().contains(key))
            .collect();
        if !pws_without_hws.is_empty() {
            // TODO include names in error message
            bail!("Distribution pipework defined for non-existent HotWaterSource(s)");
        }

        // hot water sources (not including point of use) without any pipework
        let hws_without_pws: Vec<_> = hot_water_sources
            .keys()
            .filter(|key| !hw_pipework_inputs.keys().contains(key))
            .collect();
        for hws_name in hws_without_pws {
            let hot_water_source = hot_water_sources.get(hws_name).unwrap();

            if hot_water_source.is_point_of_use() {
                // point of use doesn't need pipework - just add an empty vec
                hw_pipework_inputs.insert(hws_name.clone(), vec![]);
            } else {
                // TODO include name in error message
                bail!("Distribution pipework not specified for HotWaterSource");
            };
        }

        for (hws_name, hws) in &hot_water_sources {
            if hws.is_point_of_use() {
                continue;
            }

            for data in hw_pipework_inputs.get(hws_name).unwrap() {
                let pipework =
                    input_to_water_distribution_pipework(data, total_number_tapping_points)?;
                let entry = hot_water_distribution_pipework.get_mut(hws_name).unwrap();
                entry.push(pipework);
            }
        }

        // Set up unmet demand connection for each hot water source
        let energy_supply_conn_unmet_demand: IndexMap<String, EnergySupplyConnection> =
            hot_water_sources
                .keys()
                .map(|name| {
                    let energy_supply = energy_supplies
                        .get("_unmet_demand")
                        .expect("_unmmet_demand energy supply expected");
                    let energy_suppy_conn_unmet_demand =
                        EnergySupply::connection(energy_supply.clone(), name)
                            .expect("_unmmet_demand energy supply connection expected");
                    (name.clone(), energy_suppy_conn_unmet_demand)
                })
                .collect();

        let source_supplying_outlet = Self::init_outlet_to_source_mapping(
            showers_input,
            bath_input,
            other_hot_water_input,
            &hot_water_sources,
        );

        Ok(Self {
            showers,
            baths,
            other,
            hot_water_sources: hot_water_sources.clone(),
            energy_supply_conn_unmet_demand,
            source_supplying_outlet,
            hot_water_distribution_pipework,
            event_schedules,
            pre_heated_water_sources,
        })
    }

    fn init_outlet_to_source_mapping(
        showers_dict: ShowersInput,
        baths_dict: BathInput,
        other_hw_users_dict: OtherWaterUseInput,
        hot_water_sources: &IndexMap<String, impl HotWaterSourceBehaviour>,
    ) -> HashMap<(OutletType, String), String> {
        let mut mapping = HashMap::<(OutletType, String), String>::default();
        for (shower_name, shower) in showers_dict.0.iter() {
            match shower {
                ShowerInput::InstantElectricShower { .. } => {
                    mapping.insert(
                        (OutletType::Shower, shower_name.into()),
                        ELECTRIC_SHOWERS_HWS_NAME.into(),
                    );
                    continue;
                }
                ShowerInput::MixerShower {
                    hot_water_source, ..
                } => match hot_water_source {
                    Some(hot_water_source) => {
                        mapping.insert(
                            (OutletType::Shower, shower_name.into()),
                            hot_water_source.clone(),
                        );
                    }
                    None => {
                        if hot_water_sources.len() == 1 {
                            let default_hot_water_source = hot_water_sources.first().unwrap().0;
                            mapping.insert(
                                (OutletType::Shower, shower_name.into()),
                                default_hot_water_source.clone(),
                            );
                        } else {
                            panic!("HotWaterSource not specified for tapping point")
                        }
                    }
                },
            }
        }

        for (bath_name, bath) in baths_dict.0.iter() {
            match &bath.hot_water_source {
                Some(hot_water_source) => {
                    mapping.insert(
                        (OutletType::Bath, bath_name.into()),
                        hot_water_source.clone(),
                    );
                }
                None => {
                    if hot_water_sources.len() == 1 {
                        let default_hot_water_source = hot_water_sources.first().unwrap().0;
                        mapping.insert(
                            (OutletType::Bath, bath_name.into()),
                            default_hot_water_source.clone(),
                        );
                    } else {
                        panic!("HotWaterSource not specified for tapping point")
                    }
                }
            }
        }

        for (other_name, other) in other_hw_users_dict.0.iter() {
            match &other.hot_water_source {
                Some(hot_water_source) => {
                    mapping.insert(
                        (OutletType::Other, other_name.into()),
                        hot_water_source.clone(),
                    );
                }
                None => {
                    if hot_water_sources.len() == 1 {
                        let default_hot_water_source = hot_water_sources.first().unwrap().0;
                        mapping.insert(
                            (OutletType::Other, other_name.into()),
                            default_hot_water_source.clone(),
                        );
                    } else {
                        panic!("HotWaterSource not specified for tapping point")
                    }
                }
            }
        }

        mapping
    }

    pub(crate) fn temp_hot_water(
        &self,
        hot_water_source: T,
        volume_required_already: f64,
        volume_required: f64,
    ) -> f64 {
        let list_temperature_for_required_volume =
            hot_water_source.get_temp_hot_water(volume_required, volume_required_already);
        let sum_t_by_v: f64 = list_temperature_for_required_volume
            .iter()
            .map(|(t, v)| t * v)
            .sum();
        let sum_v: f64 = list_temperature_for_required_volume
            .iter()
            .map(|(_, v)| v)
            .sum();

        // Return average hot water temperature for the required volume
        sum_t_by_v / sum_v
    }

    pub(crate) fn get_tapping_point_for_event(
        &'_ self,
        event: TypedScheduleEvent,
    ) -> (TappingPoint<'_>, OutletType, String) {
        // TODO Results instead of panics
        match event.event_type {
            WaterScheduleEventType::Shower => {
                let shower = self
                    .showers
                    .get(&*event.name)
                    .unwrap_or_else(|| panic!("Tapping point not found for event"));
                (
                    TappingPoint::Shower(shower),
                    OutletType::Shower,
                    event.name.into(),
                )
            }
            WaterScheduleEventType::Bath => {
                let bath = self
                    .baths
                    .get(&*event.name)
                    .unwrap_or_else(|| panic!("Tapping point not found for event"));
                (
                    TappingPoint::Bath(bath),
                    OutletType::Bath,
                    event.name.into(),
                )
            }
            WaterScheduleEventType::Other => {
                let other = self
                    .other
                    .get(&*event.name)
                    .unwrap_or_else(|| panic!("Tapping point not found for event"));
                (
                    TappingPoint::Other(other),
                    OutletType::Other,
                    event.name.into(),
                )
            }
        }
    }

    pub(crate) fn hot_water_demand<'a>(
        &'a self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(
        IndexMap<String, f64>,
        IndexMap<String, f64>,
        IndexMap<String, u32>,
        IndexMap<String, f64>,
        IndexMap<String, Vec<WaterEventResult>>,
    )> {
        let hot_water_source_keys = self.hot_water_sources.keys();

        let mut hw_demand_volume: IndexMap<String, f64> = hot_water_source_keys
            .clone()
            .map(|key| (key.clone(), 0.))
            .collect();
        let mut hw_energy_demand: IndexMap<String, f64> = hot_water_source_keys
            .clone()
            .map(|key| (key.clone(), 0.))
            .collect();
        let mut hw_duration: IndexMap<String, f64> = hot_water_source_keys
            .clone()
            .map(|key| (key.clone(), 0.))
            .collect();
        let mut all_events: IndexMap<String, u32> = hot_water_source_keys
            .clone()
            .map(|key| (key.clone(), 0))
            .collect();

        hw_demand_volume.insert(ELECTRIC_SHOWERS_HWS_NAME.into(), 0.);
        hw_energy_demand.insert(ELECTRIC_SHOWERS_HWS_NAME.into(), 0.);
        hw_duration.insert(ELECTRIC_SHOWERS_HWS_NAME.into(), 0.);
        all_events.insert(ELECTRIC_SHOWERS_HWS_NAME.into(), 0);

        let mut volume_hot_water_left_in_pipework: IndexMap<String, f64> = hot_water_source_keys
            .clone()
            .map(|key| (key.clone(), 0.))
            .collect();

        for hws_name in self.hot_water_sources.keys() {
            for pipework in self.hot_water_distribution_pipework.get(hws_name).unwrap() {
                *volume_hot_water_left_in_pipework.get_mut(hws_name).unwrap() += pipework.volume()
            }
        }

        // NOTE - this is not in the Python code
        // Needed in Rust because we access by the ELECTRIC_SHOWERS_HWS_NAME key outside of a conditional
        volume_hot_water_left_in_pipework.insert(ELECTRIC_SHOWERS_HWS_NAME.into(), 0.);

        // Events have been organised now so that they are structured by simple step t_idx and
        // sorted for each time step from start to end.
        //
        // TODO (from Python) No overlapping is currently considered in terms of interaction with hot water
        // source (tank). The first event that starts is Served before the second event
        // is considered even if this starts before the previous event has finished.

        let mut usage_events: Option<Vec<TypedScheduleEvent>> =
            self.event_schedules[simtime.index].clone();

        let mut usage_events_with_flushes: IndexMap<String, Vec<WaterEventResult>> =
            hot_water_source_keys
                .map(|key| (key.clone(), vec![]))
                .collect();
        usage_events_with_flushes.insert(ELECTRIC_SHOWERS_HWS_NAME.into(), vec![]);

        if let Some(usage_events) = &mut usage_events {
            for event in usage_events.iter() {
                let (tapping_point, tapping_point_type, tapping_point_name) =
                    self.get_tapping_point_for_event(event.clone());

                let (
                    hot_water_source_name,
                    hot_water_source,
                    hw_demand_i,
                    hw_demand_target_i,
                    energy_supply_conn_unmet_demand,
                ): (
                    String,
                    Option<&T>,
                    Option<f64>,
                    f64,
                    Option<&EnergySupplyConnection>,
                ) =
                    match tapping_point {
                        TappingPoint::Shower(shower)
                            if matches!(shower, Shower::InstantElectricShower(_)) =>
                        {
                            let hot_water_source_name = ELECTRIC_SHOWERS_HWS_NAME;

                            if let Shower::InstantElectricShower(instant_electric_shower) = shower {
                                let (hw_demand_i, hw_demand_target_i) = instant_electric_shower
                                    .hot_water_demand(event.into(), simtime)?;
                                (
                                    hot_water_source_name.into(),
                                    None,
                                    Some(hw_demand_i),
                                    hw_demand_target_i,
                                    None,
                                )
                            } else {
                                unreachable!()
                            }
                        }
                        _ => {
                            let hot_water_source_name = self
                                .source_supplying_outlet
                                .get(&(tapping_point_type, tapping_point_name))
                                .unwrap();

                            let hot_water_source =
                                self.hot_water_sources.get(hot_water_source_name).unwrap();

                            let energy_supply_conn_unmet_demand = self
                                .energy_supply_conn_unmet_demand
                                .get(hot_water_source_name);

                            let volume_required_already = hw_demand_volume[hot_water_source_name];

                            let func = move |volume_required: f64| -> f64 {
                                self.temp_hot_water(
                                    hot_water_source.clone(),
                                    volume_required_already,
                                    volume_required,
                                )
                            };

                            let func_temp_hot_water: Box<dyn Fn(f64) -> f64 + 'a> = Box::new(func);
                            let (hw_demand_i, hw_demand_target_i) = tapping_point
                                .hot_water_demand(event.into(), &func_temp_hot_water, simtime)?;

                            (
                                hot_water_source_name.clone(),
                                Some(hot_water_source),
                                hw_demand_i,
                                hw_demand_target_i,
                                energy_supply_conn_unmet_demand,
                            )
                        }
                    };

                let cold_water_source = tapping_point.get_cold_water_source();

                let cold_water_temperature = if let Some(hw_demand_i) = hw_demand_i {
                    *hw_demand_volume.get_mut(&hot_water_source_name).unwrap() += hw_demand_i;
                    let list_temperature_volume = cold_water_source
                        .get_temp_cold_water(hw_demand_target_i - hw_demand_i, simtime)?;
                    let sum_t_by_v: f64 = list_temperature_volume.iter().map(|(t, v)| t * v).sum();
                    let sum_v: f64 = list_temperature_volume.iter().map(|(_, v)| v).sum();
                    sum_t_by_v / sum_v
                } else {
                    let list_temperature_volume =
                        cold_water_source.get_temp_cold_water(hw_demand_target_i, simtime)?;
                    let sum_t_by_v: f64 = list_temperature_volume.iter().map(|(t, v)| t * v).sum();
                    let sum_v: f64 = list_temperature_volume.iter().map(|(_, v)| v).sum();
                    sum_t_by_v / sum_v
                };

                let hw_energy_demand_i = water_demand_to_kwh(
                    hw_demand_target_i,
                    event.temperature,
                    cold_water_temperature,
                );

                *hw_duration.get_mut(&hot_water_source_name).unwrap() +=
                    Self::get_duration_for_tapping_point_event(&tapping_point, event);
                *all_events.get_mut(&hot_water_source_name).unwrap() += 1;

                if hw_demand_i.is_none() && energy_supply_conn_unmet_demand.is_some() {
                    // TODO check simtime.index is correct to pass here
                    energy_supply_conn_unmet_demand
                        .unwrap()
                        .demand_energy(hw_energy_demand_i, simtime.index)?;
                }

                // If event demand cannot be met, skip to the next one
                if hw_demand_i.is_none() {
                    continue;
                }

                let event_result = WaterEventResult {
                    event_result_type: event.event_type.into(),
                    temperature_warm: event.temperature,
                    volume_warm: hw_demand_target_i,
                    volume_hot: hw_demand_i.unwrap(),
                };

                // Add pipework flushes after every event (except for IES)
                let volume_hot_water_left = *volume_hot_water_left_in_pipework
                    .get(&hot_water_source_name)
                    .unwrap();

                usage_events_with_flushes
                    .get_mut(&hot_water_source_name)
                    .unwrap()
                    .push(event_result);

                if event_result.volume_hot.abs() > 1e-10
                    && volume_hot_water_left > 0.
                    && hot_water_source.is_some()
                {
                    let volume_required_already =
                        *hw_demand_volume.get(&hot_water_source_name).unwrap();
                    let volume_required = volume_hot_water_left;
                    let temperature_pipe_flush = self.temp_hot_water(
                        hot_water_source.unwrap().clone(),
                        volume_required_already,
                        volume_required,
                    );
                    usage_events_with_flushes
                        .get_mut(&hot_water_source_name)
                        .unwrap()
                        .push(WaterEventResult {
                            event_result_type: WaterEventResultType::PipeFlush,
                            temperature_warm: temperature_pipe_flush,
                            volume_warm: volume_required,
                            volume_hot: volume_required,
                        });

                    *hw_demand_volume.get_mut(&hot_water_source_name).unwrap() +=
                        volume_hot_water_left;
                }
                //  TODO (from Python)   Enhance analysis for overlapping events
                // Part of draft code for future overlapping analysis of events
                // For pipework losses count only none overlapping events
            }
        }

        // TODO (from Python) Refine pipework losses by considering overlapping of events
        // and shared pipework between serving tap points
        // none_overlapping_events calculated above is a lower bound(ish)
        // approximation for this

        // Return:
        // - litres hot water per timestep (demand on hw system)
        // - minutes demand per timestep,
        // - number of events in timestep
        // - hot water energy demand (kWh)
        // - usage_events updated to reflect pipework volumes and bath durations
        Ok((
            hw_demand_volume,
            hw_duration,
            all_events,
            hw_energy_demand,
            usage_events_with_flushes,
        ))
    }

    pub fn calc_water_heating(
        &self,
        simtime: SimulationTimeIteration,
        internal_air_temperature: f64,
        external_air_temperature: f64,
    ) -> anyhow::Result<(
        IndexMap<String, f64>,
        IndexMap<String, f64>,
        IndexMap<String, u32>,
        IndexMap<String, f64>,
        IndexMap<String, f64>,
        IndexMap<String, f64>,
        IndexMap<String, f64>,
        IndexMap<String, f64>,
        IndexMap<String, f64>,
        IndexMap<String, f64>,
    )> {
        let (
            hw_demand_vol,
            hw_duration,
            no_events,
            hw_energy_demand_at_tapping_points,
            usage_events,
        ) = self.hot_water_demand(simtime)?;

        // Running heat sources of pre-heated tanks and updating thermal losses, etc.
        for hot_water_source in self.pre_heated_water_sources.values() {
            // Note Python passes None - we are assuming an empty vec is the same behaviour
            hot_water_source.demand_hot_water(vec![], simtime)?;
        }

        let mut hw_energy_demand_at_hot_water_source: IndexMap<String, f64> = Default::default();
        let mut hw_energy_output: IndexMap<String, f64> = Default::default();
        let mut pw_losses_total: IndexMap<String, f64> = Default::default();
        let mut gains_internal_dhw: IndexMap<String, f64> = Default::default();
        let mut primary_pw_losses: IndexMap<String, f64> = Default::default();
        let mut storage_losses: IndexMap<String, f64> = Default::default();

        let mut all_keys: Vec<String> = self.hot_water_sources.keys().cloned().collect();
        all_keys.push(ELECTRIC_SHOWERS_HWS_NAME.into());

        for hws_name in all_keys {
            let (
                pw_losses_internal_for_hws,
                pw_losses_external_for_hws,
                gains_internal_dhw_use_for_hws,
            ) = self.pipework_losses_and_internal_gains_from_hot_water_events(
                hws_name.clone(),
                usage_events.get(&hws_name.clone()).unwrap(),
                internal_air_temperature,
                external_air_temperature,
            );
            pw_losses_total.insert(
                hws_name.clone(),
                pw_losses_internal_for_hws + pw_losses_external_for_hws,
            );

            // TODO check timestep is correct here
            let gains_internal_dhw_for_hws = (pw_losses_internal_for_hws
                + gains_internal_dhw_use_for_hws)
                * (WATTS_PER_KILOWATT as f64)
                / simtime.timestep;
            gains_internal_dhw.insert(hws_name, gains_internal_dhw_for_hws);
        }

        for (hws_name, hws) in &self.hot_water_sources {
            // Filtering out IES events that don't get added a 'hot_volume' when processing
            // the dhw_demand calculation
            let filtered_events: Vec<WaterEventResult> = usage_events
                .get(hws_name)
                .unwrap()
                .iter()
                .filter(|event| event.volume_hot.abs() > 1e-10)
                .copied()
                .collect();

            // TODO update demand_hot_water to accept usage_events
            hw_energy_output.insert(
                hws_name.clone(),
                hws.demand_hot_water(filtered_events.clone(), simtime)?,
            );

            // Convert from litres to kWh
            // Find underlying cold water source, ignoring pre-heat tanks

            // NOTE - Python has some logic here to find a cold water source - assumption is that we don't need that here
            let cold_water_source = hws.get_cold_water_source();
            hw_energy_demand_at_hot_water_source.insert(hws_name.clone(), 0.);
            for event in &filtered_events {
                let list_temperature_volume =
                    cold_water_source.get_temp_cold_water(event.volume_hot, simtime)?;
                let sum_t_by_v: f64 = list_temperature_volume.iter().map(|(t, v)| t * v).sum();
                let sum_v: f64 = list_temperature_volume.iter().map(|(_, v)| v).sum();
                let cold_water_temperature = sum_t_by_v / sum_v;

                *hw_energy_demand_at_hot_water_source
                    .get_mut(hws_name)
                    .unwrap() += water_demand_to_kwh(
                    event.volume_warm,
                    event.temperature_warm,
                    cold_water_temperature,
                );
            }

            let internal_gains = hws.internal_gains();
            if let Some(internal_gains) = internal_gains {
                *gains_internal_dhw.get_mut(hws_name).unwrap() += internal_gains;
            }

            let (losses, storage) = hws.get_losses_from_primary_pipework_and_storage();
            primary_pw_losses.insert(hws_name.clone(), losses);
            storage_losses.insert(hws_name.clone(), storage);
        }

        Ok((
            hw_demand_vol,
            hw_duration,
            no_events,
            hw_energy_demand_at_tapping_points,
            hw_energy_demand_at_hot_water_source,
            hw_energy_output,
            pw_losses_total,
            primary_pw_losses,
            storage_losses,
            gains_internal_dhw,
        ))
    }

    pub fn calc_pipework_losses(
        &self,
        hot_water_source_name: &String,
        no_of_hw_events: u32,
        demand_water_temperature: f64,
        internal_air_temperature: f64,
        external_air_temperature: f64,
    ) -> (f64, f64) {
        if self.hot_water_distribution_pipework.is_empty() {
            return (0., 0.);
        }

        let mut cool_down_loss_internal = 0.0;
        let mut cool_down_loss_external = 0.0;

        for pipework in self
            .hot_water_distribution_pipework
            .get(hot_water_source_name)
            .unwrap()
        {
            match pipework.location() {
                PipeworkLocation::Internal => {
                    cool_down_loss_internal += pipework.calculate_cool_down_loss(
                        demand_water_temperature,
                        internal_air_temperature,
                    );
                }
                PipeworkLocation::External => {
                    cool_down_loss_external += pipework.calculate_cool_down_loss(
                        demand_water_temperature,
                        external_air_temperature,
                    );
                }
            }
        }

        let pipework_heat_loss_internal = no_of_hw_events as f64 * cool_down_loss_internal;
        let pipework_heat_loss_external = no_of_hw_events as f64 * cool_down_loss_external;

        (pipework_heat_loss_internal, pipework_heat_loss_external)
    }

    fn pipework_losses_and_internal_gains_from_hot_water_events(
        &self,
        hot_water_source_name: String,
        usage_events: &Vec<WaterEventResult>,
        internal_air_temperature: f64,
        external_air_temperature: f64,
    ) -> (f64, f64, f64) {
        let mut pw_losses_internal = 0.;
        let mut pw_losses_external = 0.;
        let mut gains_internal_dhw_use = 0.;

        for event in usage_events {
            match event.event_result_type {
                WaterEventResultType::PipeFlush => {
                    let (pw_losses_internal_i, pw_losses_external_i) = self.calc_pipework_losses(
                        &hot_water_source_name,
                        1,
                        event.temperature_warm,
                        internal_air_temperature,
                        external_air_temperature,
                    );
                    pw_losses_internal += pw_losses_internal_i;
                    pw_losses_external += pw_losses_external_i;
                }
                _ => {
                    gains_internal_dhw_use += FRAC_DHW_ENERGY_INTERNAL_GAINS
                        * water_demand_to_kwh(
                            event.volume_warm,
                            event.temperature_warm,
                            internal_air_temperature,
                        );
                }
            }
        }

        (
            pw_losses_internal,
            pw_losses_external,
            gains_internal_dhw_use,
        )
    }

    fn get_duration_for_tapping_point_event(
        tapping_point: &TappingPoint,
        event: &TypedScheduleEvent,
    ) -> f64 {
        // In Python, for Baths, the hot_water_demand function mutates the event we pass in
        // to avoid this in the Rust we replicate the logic here
        match tapping_point {
            TappingPoint::Bath(bath) => {
                match event.duration {
                    Some(duration) => duration,
                    None => {
                        // if no duration is specified for a Bath then a volume is required
                        // to calculate the duration
                        event.volume.unwrap() / bath.get_flowrate()
                    }
                }
            }
            _ => event.duration.unwrap(),
        }
    }
}

fn shower_from_input(
    name: &str,
    input: &ShowerInput,
    cold_water_sources: &ColdWaterSources,
    energy_supplies: &IndexMap<String, Arc<RwLock<EnergySupply>>>,
    wwhrs: &IndexMap<String, Arc<Mutex<Wwhrs>>>,
) -> anyhow::Result<Shower> {
    Ok(match input {
        ShowerInput::MixerShower {
            cold_water_source,
            wwhrs_config,
            flowrate,
            hot_water_source,
        } => {
            let cold_water_source = cold_water_sources.get(cold_water_source).unwrap().clone();
            let _wwhrs_instance: Option<Arc<Mutex<Wwhrs>>> = wwhrs_config
                .as_ref()
                .map(|config| &config.waste_water_heat_recovery_system)
                .and_then(|w| wwhrs.get(w).cloned());

            Shower::MixerShower(MixerShower::new(
                *flowrate,
                cold_water_source,
                None, // TODO (migration 1.0.0a1)
                None, // TODO (migration 1.0.0a1)
            ))
        }
        ShowerInput::InstantElectricShower {
            cold_water_source,
            energy_supply,
            rated_power,
        } => {
            let cold_water_source = cold_water_sources.get(cold_water_source).unwrap().clone();

            let energy_supply = energy_supplies
                .get(energy_supply)
                .ok_or_else(|| {
                    anyhow!("The energy supply with name '{energy_supply}' is not known.")
                })?
                .clone();
            let energy_supply_conn = EnergySupply::connection(energy_supply, name)?;

            Shower::InstantElectricShower(InstantElectricShower::new(
                *rated_power,
                cold_water_source,
                energy_supply_conn,
            ))
        }
    })
}

fn input_to_bath(input: &BathDetails, cold_water_sources: &ColdWaterSources) -> Bath {
    let cold_water_source = cold_water_sources
        .get(&input.cold_water_source)
        .unwrap()
        .clone();

    Bath::new(input.size, cold_water_source, input.flowrate)
}

fn input_to_other_water_events(
    input: &OtherWaterUse,
    cold_water_sources: &ColdWaterSources,
) -> OtherHotWater {
    let cold_water_source = cold_water_sources
        .get(&input.cold_water_source)
        .unwrap()
        .clone();

    OtherHotWater::new(input.flowrate, cold_water_source)
}

fn input_to_water_distribution_pipework(
    input: &WaterPipeworkSimple,
    total_number_tapping_points: usize,
) -> anyhow::Result<PipeworkSimple> {
    // TODO check this is up to date
    // Calculate average length of pipework between HW system and tapping point
    let length_average = input.length / total_number_tapping_points as f64;

    PipeworkSimple::new(
        input.location.into(),
        input.internal_diameter_mm / MILLIMETRES_IN_METRE as f64,
        length_average,
        PipeworkContents::Water,
    )
    .map_err(anyhow::Error::msg)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compare_floats::max_of_2;
    use crate::core::common::WaterSupply;
    use crate::core::energy_supply::energy_supply::EnergySupplyBuilder;
    use crate::core::heating_systems::wwhrs::{WWHRSInstantaneousSystemB, Wwhrs};
    use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
    use crate::input::{
        Baths, FuelType, MixerShowerWwhrsConfiguration, OtherWaterUses, Showers,
        WaterPipeworkLocation,
    };
    use crate::simulation_time::SimulationTime;
    use parking_lot::RwLock;
    use pretty_assertions::assert_eq;
    use rstest::*;

    #[derive(Debug, Clone)]
    struct HotWaterSourceMock {
        cold_feed: WaterSupply,
    }

    impl HotWaterSourceBehaviour for HotWaterSourceMock {
        fn get_cold_water_source(&self) -> WaterSupply {
            self.cold_feed.clone()
        }

        fn temp_hot_water(&self) -> anyhow::Result<f64> {
            unimplemented!();
        }

        fn demand_hot_water(
            &self,
            usage_events: Vec<WaterEventResult>,
            _: SimulationTimeIteration,
        ) -> anyhow::Result<f64> {
            Ok(usage_events
                .iter()
                .map(|e| e.temperature_warm * e.volume_warm / 4200.0)
                .sum())
        }

        fn get_temp_hot_water(
            &self,
            volume_required: f64,
            volume_required_already: f64,
        ) -> Vec<(f64, f64)> {
            let volume_req_cumulative = volume_required + volume_required_already;
            let frac_volume_req_already = volume_required_already / volume_req_cumulative;
            let frac_layer_1 = max_of_2(0.0, 0.6 - frac_volume_req_already);
            let frac_layer_2 = 0.4 - max_of_2(0.0, frac_volume_req_already - 0.6);

            vec![
                (55.0, volume_req_cumulative * frac_layer_1),
                (45.0, volume_req_cumulative * frac_layer_2),
            ]
        }

        fn internal_gains(&self) -> Option<f64> {
            None
        }

        fn get_losses_from_primary_pipework_and_storage(&self) -> (f64, f64) {
            (0., 0.)
        }

        fn is_point_of_use(&self) -> bool {
            false
        }
    }

    #[derive(Debug, Clone)]
    struct HotWaterSourceMockWithUniqueHotWaterTemperature {}

    impl HotWaterSourceBehaviour for HotWaterSourceMockWithUniqueHotWaterTemperature {
        fn get_cold_water_source(&self) -> WaterSupply {
            unimplemented!()
        }

        fn temp_hot_water(&self) -> anyhow::Result<f64> {
            todo!()
        }

        fn demand_hot_water(
            &self,
            _usage_events: Vec<WaterEventResult>,
            _simtime: SimulationTimeIteration,
        ) -> anyhow::Result<f64> {
            Ok(0.) // Python doesn't mock this
        }

        fn get_temp_hot_water(
            &self,
            volume_required: f64,
            volume_required_already: f64,
        ) -> Vec<(f64, f64)> {
            let volume_req_cumulative = volume_required + volume_required_already;
            vec![(55.0, volume_req_cumulative)]
        }

        fn internal_gains(&self) -> Option<f64> {
            Some(0.) // Python doesn't mock this
        }

        fn get_losses_from_primary_pipework_and_storage(&self) -> (f64, f64) {
            (0., 0.) // Python doesn't mock this
        }

        fn is_point_of_use(&self) -> bool {
            false
        }
    }

    #[derive(Debug, Clone)]
    struct HotWaterSourceMockWithInternalGains {
        cold_feed: WaterSupply,
    }

    impl HotWaterSourceBehaviour for HotWaterSourceMockWithInternalGains {
        fn get_cold_water_source(&self) -> WaterSupply {
            self.cold_feed.clone()
        }

        fn temp_hot_water(&self) -> anyhow::Result<f64> {
            unimplemented!();
        }

        fn demand_hot_water(
            &self,
            usage_events: Vec<WaterEventResult>,
            _: SimulationTimeIteration,
        ) -> anyhow::Result<f64> {
            Ok(usage_events
                .iter()
                .map(|e| e.temperature_warm * e.volume_warm / 4200.0)
                .sum())
        }

        fn get_temp_hot_water(
            &self,
            volume_required: f64,
            volume_required_already: f64,
        ) -> Vec<(f64, f64)> {
            let volume_req_cumulative = volume_required + volume_required_already;
            let frac_volume_req_already = volume_required_already / volume_req_cumulative;
            let frac_layer_1 = max_of_2(0.0, 0.6 - frac_volume_req_already);
            let frac_layer_2 = 0.4 - max_of_2(0.0, frac_volume_req_already - 0.6);

            vec![
                (55.0, volume_req_cumulative * frac_layer_1),
                (45.0, volume_req_cumulative * frac_layer_2),
            ]
        }

        fn internal_gains(&self) -> Option<f64> {
            Some(20.)
        }

        fn get_losses_from_primary_pipework_and_storage(&self) -> (f64, f64) {
            (5., 15.)
        }

        fn is_point_of_use(&self) -> bool {
            false
        }
    }

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 24., 1.)
    }

    #[fixture]
    fn cold_water_source() -> Arc<ColdWaterSource> {
        let cold_water_temps = vec![
            2.0, 3.0, 4.0, 2.0, 3.0, 4.0, 2.0, 3.0, 4.0, 2.0, 3.0, 4.0, 2.0, 3.0, 4.0, 2.0, 3.0,
            4.0, 2.0, 3.0, 4.0, 2.0, 3.0, 4.0,
        ];
        Arc::from(ColdWaterSource::new(cold_water_temps, 0, 1.))
    }

    #[fixture]
    fn event_schedules() -> Vec<Option<Vec<TypedScheduleEvent>>> {
        vec![
            None,
            None,
            None,
            None,
            Some(vec![
                TypedScheduleEvent {
                    start: 4.1,
                    duration: Some(6.),
                    temperature: 41.0,
                    name: "IES".to_string(),
                    event_type: WaterScheduleEventType::Shower,
                    volume: None,
                    warm_volume: None,
                    pipework_volume: None,
                },
                TypedScheduleEvent {
                    start: 4.5,
                    duration: Some(6.),
                    temperature: 41.0,
                    name: "IES".to_string(),
                    event_type: WaterScheduleEventType::Shower,
                    volume: None,
                    warm_volume: None,
                    pipework_volume: None,
                },
            ]),
            None,
            Some(vec![
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
                    duration: None, // TODO match Python
                    // duration: Some(1.), // set duration to Some to avoid unwrap error
                    temperature: 41.0,
                    name: "medium".to_string(),
                    event_type: WaterScheduleEventType::Bath,
                    volume: Some(100.0),
                    warm_volume: None,
                    pipework_volume: None,
                },
            ]),
            Some(vec![
                TypedScheduleEvent {
                    start: 7.,
                    duration: Some(6.),
                    temperature: 41.0,
                    name: "mixer".to_string(),
                    event_type: WaterScheduleEventType::Shower,
                    volume: None,
                    warm_volume: None,
                    pipework_volume: None,
                },
                TypedScheduleEvent {
                    start: 7.,
                    duration: Some(1.),
                    temperature: 41.0,
                    name: "other".to_string(),
                    event_type: WaterScheduleEventType::Other,
                    volume: None,
                    warm_volume: None,
                    pipework_volume: None,
                },
            ]),
            Some(vec![
                TypedScheduleEvent {
                    start: 8.,
                    duration: Some(6.),
                    temperature: 41.0,
                    name: "mixer".to_string(),
                    event_type: WaterScheduleEventType::Shower,
                    volume: None,
                    warm_volume: None,
                    pipework_volume: None,
                },
                TypedScheduleEvent {
                    start: 8.,
                    duration: Some(3.),
                    temperature: 41.0,
                    name: "medium".to_string(),
                    event_type: WaterScheduleEventType::Bath,
                    volume: None,
                    warm_volume: None,
                    pipework_volume: None,
                },
            ]),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        ]
    }

    fn create_dhw_demand<T: HotWaterSourceBehaviour>(
        simulation_time: SimulationTime,
        hot_water_sources: IndexMap<String, T>,
        pre_heated_water_sources: IndexMap<String, T>,
        cold_water_source: Arc<ColdWaterSource>,
        event_schedules: Vec<Option<Vec<TypedScheduleEvent>>>,
    ) -> DomesticHotWaterDemand<T> {
        let flow_rates = vec![5., 7., 9., 11., 13.];
        let efficiencies = vec![44.8, 39.1, 34.8, 31.4, 28.6];
        let utilisation_factor = 0.7;
        let cold_water_sources =
            ColdWaterSources::from([("mains water".into(), cold_water_source.clone())]);
        let wwhrsb = Arc::new(Mutex::new(Wwhrs::WWHRSInstantaneousSystemB(
            WWHRSInstantaneousSystemB::new(
                cold_water_source.clone(),
                flow_rates,
                efficiencies,
                utilisation_factor,
            ),
        )));
        let wwhrs = IndexMap::from([(String::from("Example_Inst_WWHRS"), wwhrsb.clone())]);

        let electricity_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));

        let unmet_demand = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::UnmetDemand, simulation_time.total_steps()).build(),
        ));

        let energy_supplies = IndexMap::from([
            ("_unmet_demand".into(), unmet_demand.clone()),
            ("mains elec".into(), electricity_supply.clone()),
        ]);

        let showers_input = Showers(IndexMap::from([
            (
                "mixer".into(),
                ShowerInput::MixerShower {
                    flowrate: 8.0,
                    cold_water_source: "mains water".into(),
                    wwhrs_config: Some(MixerShowerWwhrsConfiguration {
                        waste_water_heat_recovery_system: "Example_Inst_WWHRS".into(),
                        ..Default::default()
                    }),
                    hot_water_source: None,
                },
            ),
            (
                "IES".into(),
                ShowerInput::InstantElectricShower {
                    rated_power: 9.0,
                    cold_water_source: "mains water".into(),
                    energy_supply: "mains elec".into(),
                },
            ),
        ]));

        let baths_input = Baths(IndexMap::from([(
            "medium".into(),
            BathDetails {
                size: 100.,
                cold_water_source: "mains water".into(),
                flowrate: 8.0,
                hot_water_source: None,
            },
        )]));

        let other_input = OtherWaterUses(IndexMap::from([(
            "other".into(),
            OtherWaterUse {
                flowrate: 8.0,
                cold_water_source: "mains water".into(),
                hot_water_source: None,
            },
        )]));

        let hw_pipework = WaterDistributionInput::List(vec![
            WaterPipeworkSimple {
                location: WaterPipeworkLocation::Internal,
                internal_diameter_mm: 30.,
                length: 10.0,
                external_diameter_mm: None,
                insulation_thermal_conductivity: None,
                insulation_thickness_mm: None,
                surface_reflectivity: None,
                pipe_contents: None,
            },
            WaterPipeworkSimple {
                location: WaterPipeworkLocation::Internal,
                internal_diameter_mm: 28.,
                length: 9.0,
                external_diameter_mm: None,
                insulation_thermal_conductivity: None,
                insulation_thickness_mm: None,
                surface_reflectivity: None,
                pipe_contents: None,
            },
            WaterPipeworkSimple {
                location: WaterPipeworkLocation::External,
                internal_diameter_mm: 32.,
                length: 5.0,
                external_diameter_mm: None,
                insulation_thermal_conductivity: None,
                insulation_thickness_mm: None,
                surface_reflectivity: None,
                pipe_contents: None,
            },
            WaterPipeworkSimple {
                location: WaterPipeworkLocation::External,
                internal_diameter_mm: 31.,
                length: 8.0,
                external_diameter_mm: None,
                insulation_thermal_conductivity: None,
                insulation_thickness_mm: None,
                surface_reflectivity: None,
                pipe_contents: None,
            },
        ]);

        DomesticHotWaterDemand::new(
            showers_input,
            baths_input,
            other_input,
            hw_pipework,
            &cold_water_sources,
            &wwhrs,
            &energy_supplies,
            event_schedules,
            hot_water_sources,
            pre_heated_water_sources,
        )
        .unwrap()
    }

    #[rstest]
    #[ignore = "work in progress for migration to 1.0.0a1"]
    fn test_hot_water_demand(
        simulation_time: SimulationTime,
        event_schedules: Vec<Option<Vec<TypedScheduleEvent>>>,
        cold_water_source: Arc<ColdWaterSource>,
    ) {
        let hot_water_sources: IndexMap<String, HotWaterSourceMockWithUniqueHotWaterTemperature> =
            IndexMap::from([(
                "hw cylinder".into(),
                HotWaterSourceMockWithUniqueHotWaterTemperature {},
            )]);
        let pre_heated_water_sources: IndexMap<
            String,
            HotWaterSourceMockWithUniqueHotWaterTemperature,
        > = Default::default();

        let dhw_demand = create_dhw_demand(
            simulation_time,
            hot_water_sources,
            pre_heated_water_sources,
            cold_water_source,
            event_schedules,
        );

        let expected_results: Vec<(
            IndexMap<String, f64>,
            IndexMap<String, f64>,
            IndexMap<String, u32>,
            IndexMap<String, f64>,
            IndexMap<String, Vec<WaterEventResult>>,
        )> = vec![
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), 12.),
                    ("hw cylinder".into(), 0.),
                ]),
                IndexMap::from([("_electric_showers".into(), 2), ("hw cylinder".into(), 0)]),
                IndexMap::from([
                    ("_electric_showers".into(), 1.7999999999999998),
                    ("hw cylinder".into(), 0.),
                ]),
                IndexMap::from([
                    (
                        "_electric_showers".into(),
                        vec![
                            WaterEventResult {
                                event_result_type: WaterEventResultType::Shower,
                                temperature_warm: 41.,
                                volume_warm: 20.378383818053738,
                                volume_hot: 0.,
                            },
                            WaterEventResult {
                                event_result_type: WaterEventResultType::Shower,
                                temperature_warm: 41.,
                                volume_warm: 20.378383818053738,
                                volume_hot: 0.,
                            },
                        ],
                    ),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([
                    ("_electric_showers".into(), 0.),
                    ("hw cylinder".into(), 81.141483189812),
                ]),
                IndexMap::from([
                    ("_electric_showers".into(), 0.),
                    ("hw cylinder".into(), 12.5),
                ]),
                IndexMap::from([("_electric_showers".into(), 1), ("hw cylinder".into(), 1)]),
                IndexMap::from([
                    ("_electric_showers".into(), 0.9),
                    ("hw cylinder".into(), 4.532666666666667),
                ]),
                IndexMap::from([
                    (
                        "_electric_showers".into(),
                        vec![WaterEventResult {
                            event_result_type: WaterEventResultType::Shower,
                            temperature_warm: 41.,
                            volume_warm: 19.85586115605236,
                            volume_hot: 0.,
                        }],
                    ),
                    (
                        "hw cylinder".into(),
                        vec![
                            WaterEventResult {
                                event_result_type: WaterEventResultType::Bath,
                                temperature_warm: 41.,
                                volume_warm: 100.,
                                volume_hot: 73.58490566037736,
                            },
                            WaterEventResult {
                                event_result_type: WaterEventResultType::PipeFlush,
                                temperature_warm: 55.,
                                volume_warm: 7.556577529434648,
                                volume_hot: 7.556577529434648,
                            },
                        ],
                    ),
                ]),
            ),
            (
                IndexMap::from([
                    ("_electric_showers".into(), 0.),
                    ("hw cylinder".into(), 53.58989404060333),
                ]),
                IndexMap::from([
                    ("_electric_showers".into(), 0.),
                    ("hw cylinder".into(), 7.0),
                ]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 2)]),
                IndexMap::from([
                    ("_electric_showers".into(), 0.),
                    ("hw cylinder".into(), 2.473208888888889),
                ]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    (
                        "hw cylinder".into(),
                        vec![
                            WaterEventResult {
                                event_result_type: WaterEventResultType::Shower,
                                temperature_warm: 41.,
                                volume_warm: 48.,
                                volume_hot: 32.63058513558019,
                            },
                            WaterEventResult {
                                event_result_type: WaterEventResultType::PipeFlush,
                                temperature_warm: 55.,
                                volume_warm: 7.556577529434648,
                                volume_hot: 7.556577529434648,
                            },
                            WaterEventResult {
                                event_result_type: WaterEventResultType::PipeFlush,
                                temperature_warm: 41.,
                                volume_warm: 8.,
                                volume_hot: 5.846153846153845,
                            },
                            WaterEventResult {
                                event_result_type: WaterEventResultType::PipeFlush,
                                temperature_warm: 55.,
                                volume_warm: 7.556577529434648,
                                volume_hot: 7.556577529434648,
                            },
                        ],
                    ),
                ]),
            ),
            (
                IndexMap::from([
                    ("_electric_showers".into(), 0.),
                    ("hw cylinder".into(), 64.89041357202146),
                ]),
                IndexMap::from([
                    ("_electric_showers".into(), 0.),
                    ("hw cylinder".into(), 2.0),
                ]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 2)]),
                IndexMap::from([
                    ("_electric_showers".into(), 0.),
                    ("hw cylinder".into(), 3.09616),
                ]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    (
                        "hw cylinder".into(),
                        vec![
                            WaterEventResult {
                                event_result_type: WaterEventResultType::PipeFlush,
                                temperature_warm: 41.,
                                volume_warm: 48.,
                                volume_hot: 32.365493807269814,
                            },
                            WaterEventResult {
                                event_result_type: WaterEventResultType::PipeFlush,
                                temperature_warm: 55.,
                                volume_warm: 7.556577529434648,
                                volume_hot: 7.556577529434648,
                            },
                            WaterEventResult {
                                event_result_type: WaterEventResultType::Bath,
                                temperature_warm: 55.,
                                volume_warm: 24.,
                                volume_hot: 17.41176470588235,
                            },
                            WaterEventResult {
                                event_result_type: WaterEventResultType::PipeFlush,
                                temperature_warm: 54.99999999999999,
                                volume_warm: 7.556577529434648,
                                volume_hot: 7.556577529434648,
                            },
                        ],
                    ),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
            (
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 0)]),
                IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
                IndexMap::from([
                    ("_electric_showers".into(), vec![]),
                    ("hw cylinder".into(), vec![]),
                ]),
            ),
        ];

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let (
                hw_demand_volume,
                hw_duration,
                all_events,
                hw_energy_demand,
                _usage_events_with_flushes,
            ) = dhw_demand.hot_water_demand(t_it).unwrap();

            let (
                expected_hw_demand_volume,
                expected_hw_duration,
                expected_all_events,
                expected_hw_energy_demand,
                _expected_usage_events_with_flushes,
            ) = &expected_results[t_idx];

            assert_eq!(hw_demand_volume, *expected_hw_demand_volume);
            assert_eq!(hw_duration, *expected_hw_duration);
            assert_eq!(all_events, *expected_all_events);
            assert_eq!(hw_energy_demand, *expected_hw_energy_demand);

            // TODO allow comparison of this
            // assert_eq!(_usage_events_with_flushes, *_expected_usage_events_with_flushes);
        }
    }

    // Skipping Python test test_unmet_demand_connection_failure

    #[rstest]
    #[ignore = "work in progress for migration to 1.0.0a1"]
    fn test_hot_water_demand_unmet(
        simulation_time: SimulationTime,
        event_schedules: Vec<Option<Vec<TypedScheduleEvent>>>,
        cold_water_source: Arc<ColdWaterSource>,
    ) {
        let mut event_schedules_modified = event_schedules.clone();
        event_schedules_modified[0] = Some(vec![
            TypedScheduleEvent {
                start: 0.,
                duration: None,
                temperature: 90.,
                name: "medium".into(),
                event_type: WaterScheduleEventType::Bath,
                volume: Some(100.),
                warm_volume: None,
                pipework_volume: None,
            },
            TypedScheduleEvent {
                start: 0.,
                duration: Some(6.),
                temperature: 70.,
                name: "mixer".into(),
                event_type: WaterScheduleEventType::Shower,
                volume: None,
                warm_volume: None,
                pipework_volume: None,
            },
            TypedScheduleEvent {
                start: 0.,
                duration: Some(1.),
                temperature: 80.,
                name: "other".into(),
                event_type: WaterScheduleEventType::Other,
                volume: None,
                warm_volume: None,
                pipework_volume: None,
            },
        ]);

        let hot_water_sources: IndexMap<String, HotWaterSourceMockWithUniqueHotWaterTemperature> =
            IndexMap::from([(
                "hw cylinder".into(),
                HotWaterSourceMockWithUniqueHotWaterTemperature {},
            )]);

        let pre_heated_water_sources: IndexMap<
            String,
            HotWaterSourceMockWithUniqueHotWaterTemperature,
        > = Default::default();

        let dhw_demand = create_dhw_demand(
            simulation_time,
            hot_water_sources,
            pre_heated_water_sources,
            cold_water_source,
            event_schedules_modified,
        );

        let expected: (
            IndexMap<String, f64>,
            IndexMap<String, f64>,
            IndexMap<String, u32>,
            IndexMap<String, f64>,
            IndexMap<String, Vec<WaterEventResult>>,
        ) = (
            IndexMap::from([("_electric_showers".into(), 0.), ("hw cylinder".into(), 0.)]),
            IndexMap::from([
                ("_electric_showers".into(), 0.),
                ("hw cylinder".into(), 19.5),
            ]),
            IndexMap::from([("_electric_showers".into(), 0), ("hw cylinder".into(), 3)]),
            IndexMap::from([
                ("_electric_showers".into(), 0.),
                ("hw cylinder".into(), 14.746275555555554),
            ]),
            IndexMap::from([
                ("_electric_showers".into(), vec![]),
                ("hw cylinder".into(), vec![]),
            ]),
        );

        let (
            expected_hw_demand_volume,
            expected_hw_duration,
            expected_all_events,
            expected_hw_energy_demand,
            _expected_usage_events_with_flushes,
        ) = expected;

        let actual = dhw_demand
            .hot_water_demand(simulation_time.iter().current_iteration())
            .unwrap();

        let (
            hw_demand_volume,
            hw_duration,
            all_events,
            hw_energy_demand,
            _usage_events_with_flushes,
        ) = actual;

        assert_eq!(hw_demand_volume, expected_hw_demand_volume);
        assert_eq!(hw_duration, expected_hw_duration);
        assert_eq!(all_events, expected_all_events);
        assert_eq!(hw_energy_demand, expected_hw_energy_demand);

        // TODO allow comparison of this
        // assert_eq!(_usage_events_with_flushes, *_expected_usage_events_with_flushes);
    }

    #[rstest]
    #[ignore = "work in progress for migration to 1.0.0a1"]
    fn test_calc_water_heating(
        simulation_time: SimulationTime,
        event_schedules: Vec<Option<Vec<TypedScheduleEvent>>>,
        cold_water_source: Arc<ColdWaterSource>,
    ) {
        let hot_water_sources: IndexMap<String, HotWaterSourceMock> = IndexMap::from([(
            "hw cylinder".into(),
            HotWaterSourceMock {
                cold_feed: WaterSupply::ColdWaterSource(cold_water_source.clone()),
            },
        )]);

        let pre_heated_water_sources: IndexMap<String, HotWaterSourceMock> = IndexMap::from([(
            "pre-heat tank".into(),
            HotWaterSourceMock {
                cold_feed: WaterSupply::ColdWaterSource(cold_water_source.clone()),
            },
        )]);

        let temp_int_air = 20.;
        let temp_ext_air = 5.;

        let dhw_demand = create_dhw_demand(
            simulation_time,
            hot_water_sources,
            pre_heated_water_sources,
            cold_water_source,
            event_schedules,
        );

        let hw_demand_vol_expected: IndexMap<String, Vec<f64>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    87.14841426412852,
                    58.267631655968856,
                    72.45826885901587,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            ),
        ]);
        let hw_duration_expected: IndexMap<String, Vec<f64>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.5, 7.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0.0, 0.0, 0.0, 0.0, 12.0, 0.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            ),
        ]);
        let no_events_expected: IndexMap<String, Vec<u32>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    0, 0, 0, 0, 0, 0, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                ],
            ),
        ]);
        let hw_energy_demand_expected: IndexMap<String, Vec<f64>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    4.532666666666667,
                    2.473208888888889,
                    3.09616,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.7999999999999998,
                    0.0,
                    0.9,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ),
        ]);
        let hw_energy_demand_incl_pipework_loss_expected: IndexMap<String, Vec<f64>> =
            IndexMap::from([
                (
                    "hw cylinder".into(),
                    vec![
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        4.91031082679879,
                        3.2109323644958288,
                        3.8163186309496315,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                    ],
                ), // NOTE - no _electric_showers entry expected
            ]);
        let hw_energy_output_expected: IndexMap<String, Vec<f64>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0571538068629902,
                    0.7085933280116948,
                    0.864783804202171,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            ),
        ]);
        let dist_pw_losses_expected: IndexMap<String, Vec<f64>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.2780167312270571,
                    0.5560334624541142,
                    0.5560334624541142,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            ),
        ]);
        let primary_pw_losses_expected: IndexMap<String, Vec<f64>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            ),
        ]);
        let storage_losses_expected: IndexMap<String, Vec<f64>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            ),
        ]);
        let gains_internal_dhw_expected: IndexMap<String, Vec<f64>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    732.3002698651746,
                    585.9605397303493,
                    683.5872063970161,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    248.68421052631578,
                    0.0,
                    121.15384615384616,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ),
        ]);

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let actual = dhw_demand
                .calc_water_heating(t_it, temp_int_air, temp_ext_air)
                .unwrap();

            let (
                hw_demand_vol,
                hw_duration,
                no_events,
                hw_energy_demand,
                hw_energy_demand_incl_pipework_loss,
                hw_energy_output,
                dist_pw_losses,
                primary_pw_losses,
                storage_losses,
                gains_internal_dhw,
            ) = actual;

            let keys: Vec<String> = vec!["hw cylinder".into(), "_electric_showers".into()];
            for key in keys {
                assert_eq!(
                    *hw_demand_vol.get(&key).unwrap(),
                    hw_demand_vol_expected.get(&key).unwrap()[t_idx]
                );
                assert_eq!(
                    *hw_duration.get(&key).unwrap(),
                    hw_duration_expected.get(&key).unwrap()[t_idx]
                );
                assert_eq!(
                    *no_events.get(&key).unwrap(),
                    no_events_expected.get(&key).unwrap()[t_idx]
                );
                assert_eq!(
                    *hw_energy_demand.get(&key).unwrap(),
                    hw_energy_demand_expected.get(&key).unwrap()[t_idx]
                );
                assert_eq!(
                    *dist_pw_losses.get(&key).unwrap(),
                    dist_pw_losses_expected.get(&key).unwrap()[t_idx]
                );
                assert_eq!(
                    *gains_internal_dhw.get(&key).unwrap(),
                    gains_internal_dhw_expected.get(&key).unwrap()[t_idx]
                );

                // don't check for _electric_showers entry for the below
                if &key == "_electric_showers" {
                    continue;
                };
                assert_eq!(
                    *hw_energy_demand_incl_pipework_loss.get(&key).unwrap(),
                    hw_energy_demand_incl_pipework_loss_expected
                        .get(&key)
                        .unwrap()[t_idx]
                );
                assert_eq!(
                    *hw_energy_output.get(&key).unwrap(),
                    hw_energy_output_expected.get(&key).unwrap()[t_idx]
                );
                assert_eq!(
                    *primary_pw_losses.get(&key).unwrap(),
                    primary_pw_losses_expected.get(&key).unwrap()[t_idx]
                );
                assert_eq!(
                    *storage_losses.get(&key).unwrap(),
                    storage_losses_expected.get(&key).unwrap()[t_idx]
                );
            }
        }
    }

    #[rstest]
    #[ignore = "work in progress for migration to 1.0.0a1"]
    fn test_calc_water_heating_with_internal_gains(
        simulation_time: SimulationTime,
        event_schedules: Vec<Option<Vec<TypedScheduleEvent>>>,
        cold_water_source: Arc<ColdWaterSource>,
    ) {
        let hot_water_sources: IndexMap<String, HotWaterSourceMockWithInternalGains> =
            IndexMap::from([(
                "hw cylinder".into(),
                HotWaterSourceMockWithInternalGains {
                    cold_feed: WaterSupply::ColdWaterSource(cold_water_source.clone()),
                },
            )]);

        let pre_heated_water_sources: IndexMap<String, HotWaterSourceMockWithInternalGains> =
            IndexMap::from([]);

        let temp_int_air = 20.;
        let temp_ext_air = 5.;

        let dhw_demand = create_dhw_demand(
            simulation_time,
            hot_water_sources,
            pre_heated_water_sources,
            cold_water_source,
            event_schedules,
        );

        let hw_demand_vol_expected: IndexMap<String, Vec<f64>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    87.14841426412852,
                    58.267631655968856,
                    72.45826885901587,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            ),
        ]);
        let hw_duration_expected: IndexMap<String, Vec<f64>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.5, 7.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0.0, 0.0, 0.0, 0.0, 12.0, 0.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            ),
        ]);
        let no_events_expected: IndexMap<String, Vec<u32>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    0, 0, 0, 0, 0, 0, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                ],
            ),
        ]);
        let hw_energy_demand_expected: IndexMap<String, Vec<f64>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    4.532666666666667,
                    2.473208888888889,
                    3.09616,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.7999999999999998,
                    0.0,
                    0.9,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ),
        ]);
        let hw_energy_demand_incl_pipework_loss_expected: IndexMap<String, Vec<f64>> =
            IndexMap::from([
                (
                    "hw cylinder".into(),
                    vec![
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        4.91031082679879,
                        3.2109323644958288,
                        3.8163186309496315,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                    ],
                ), // NOTE - no _electric_showers entry expected
            ]);
        let hw_energy_output_expected: IndexMap<String, Vec<f64>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0571538068629902,
                    0.7085933280116948,
                    0.864783804202171,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            ),
        ]);
        let dist_pw_losses_expected: IndexMap<String, Vec<f64>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.2780167312270571,
                    0.5560334624541142,
                    0.5560334624541142,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            ),
        ]);
        let primary_pw_losses_expected: IndexMap<String, Vec<f64>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
                    5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            ),
        ]);
        let storage_losses_expected: IndexMap<String, Vec<f64>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0,
                    15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                ],
            ),
        ]);
        let gains_internal_dhw_expected: IndexMap<String, Vec<f64>> = IndexMap::from([
            (
                "hw cylinder".into(),
                vec![
                    20.0,
                    20.0,
                    20.0,
                    20.0,
                    20.0,
                    20.0,
                    752.3002698651746,
                    605.9605397303493,
                    703.5872063970161,
                    20.0,
                    20.0,
                    20.0,
                    20.0,
                    20.0,
                    20.0,
                    20.0,
                    20.0,
                    20.0,
                    20.0,
                    20.0,
                    20.0,
                    20.0,
                    20.0,
                    20.0,
                ],
            ),
            (
                "_electric_showers".into(),
                vec![
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    248.68421052631578,
                    0.0,
                    121.15384615384616,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ),
        ]);

        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            let actual = dhw_demand
                .calc_water_heating(t_it, temp_int_air, temp_ext_air)
                .unwrap();

            let (
                hw_demand_vol,
                hw_duration,
                no_events,
                hw_energy_demand,
                hw_energy_demand_incl_pipework_loss,
                hw_energy_output,
                dist_pw_losses,
                primary_pw_losses,
                storage_losses,
                gains_internal_dhw,
            ) = actual;

            let keys: Vec<String> = vec!["hw cylinder".into(), "_electric_showers".into()];
            for key in keys {
                assert_eq!(
                    *hw_demand_vol.get(&key).unwrap(),
                    hw_demand_vol_expected.get(&key).unwrap()[t_idx]
                );
                assert_eq!(
                    *hw_duration.get(&key).unwrap(),
                    hw_duration_expected.get(&key).unwrap()[t_idx]
                );
                assert_eq!(
                    *no_events.get(&key).unwrap(),
                    no_events_expected.get(&key).unwrap()[t_idx]
                );
                assert_eq!(
                    *hw_energy_demand.get(&key).unwrap(),
                    hw_energy_demand_expected.get(&key).unwrap()[t_idx]
                );
                assert_eq!(
                    *dist_pw_losses.get(&key).unwrap(),
                    dist_pw_losses_expected.get(&key).unwrap()[t_idx]
                );
                assert_eq!(
                    *gains_internal_dhw.get(&key).unwrap(),
                    gains_internal_dhw_expected.get(&key).unwrap()[t_idx]
                );

                // don't check for _electric_showers entry for the below
                if &key == "_electric_showers" {
                    continue;
                };
                assert_eq!(
                    *hw_energy_demand_incl_pipework_loss.get(&key).unwrap(),
                    hw_energy_demand_incl_pipework_loss_expected
                        .get(&key)
                        .unwrap()[t_idx]
                );
                assert_eq!(
                    *hw_energy_output.get(&key).unwrap(),
                    hw_energy_output_expected.get(&key).unwrap()[t_idx]
                );
                assert_eq!(
                    *primary_pw_losses.get(&key).unwrap(),
                    primary_pw_losses_expected.get(&key).unwrap()[t_idx]
                );
                assert_eq!(
                    *storage_losses.get(&key).unwrap(),
                    storage_losses_expected.get(&key).unwrap()[t_idx]
                );
            }
        }
    }

    #[rstest]
    fn test_calc_pipework_losses(
        simulation_time: SimulationTime,
        event_schedules: Vec<Option<Vec<TypedScheduleEvent>>>,
        cold_water_source: Arc<ColdWaterSource>,
    ) {
        let hot_water_sources: IndexMap<String, HotWaterSourceMockWithInternalGains> =
            IndexMap::from([(
                "hw cylinder".into(),
                HotWaterSourceMockWithInternalGains {
                    cold_feed: WaterSupply::ColdWaterSource(cold_water_source.clone()),
                },
            )]);

        let pre_heated_water_sources: IndexMap<String, HotWaterSourceMockWithInternalGains> =
            IndexMap::from([]);

        let dhw_demand = create_dhw_demand(
            simulation_time,
            hot_water_sources,
            pre_heated_water_sources,
            cold_water_source,
            event_schedules,
        );

        let no_of_hw_events: [u32; 24] = [
            0, 0, 0, 2, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let demand_water_temperature: [f64; 24] = [
            55.0, 55.0, 55.0, 57.0, 55.0, 58.0, 60.0, 40.0, 55.0, 55.0, 55.0, 55.0, 55.0, 55.0,
            55.0, 55.0, 55.0, 55.0, 55.0, 55.0, 55.0, 55.0, 55.0, 55.0,
        ];
        let expected: [(f64, f64); 24] = [
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.3615154654675836, 0.4052961328742277),
            (0.0, 0.0),
            (0.3712861537234643, 0.41309028927565516),
            (0.3908275302352256, 0.42867860207851005),
            (0.1954137651176128, 0.27279547404996096),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 0.0),
        ];

        for (t_idx, _) in simulation_time.iter().enumerate() {
            let actual = dhw_demand.calc_pipework_losses(
                &"hw cylinder".into(),
                no_of_hw_events[t_idx],
                demand_water_temperature[t_idx],
                20.,
                5.,
            );
            assert_eq!(actual, expected[t_idx]);
        }
    }

    #[rstest]
    fn test_calc_pipework_losses_with_no_pipework(
        simulation_time: SimulationTime,
        event_schedules: Vec<Option<Vec<TypedScheduleEvent>>>,
        cold_water_source: Arc<ColdWaterSource>,
    ) {
        let hot_water_sources: IndexMap<String, HotWaterSourceMockWithInternalGains> =
            IndexMap::from([(
                "hw cylinder".into(),
                HotWaterSourceMockWithInternalGains {
                    cold_feed: WaterSupply::ColdWaterSource(cold_water_source.clone()),
                },
            )]);

        let pre_heated_water_sources: IndexMap<String, HotWaterSourceMockWithInternalGains> =
            IndexMap::from([]);

        let dhw_demand = create_dhw_demand(
            simulation_time,
            hot_water_sources,
            pre_heated_water_sources,
            cold_water_source,
            event_schedules,
        );

        let actual = dhw_demand.calc_pipework_losses(&"hw cylinder".into(), 0, 55., 20., 5.);
        assert_eq!(actual, (0., 0.));
    }

    #[rstest]
    fn test_hot_water_demand_with_other_event_new_temp(
        simulation_time: SimulationTime,
        cold_water_source: Arc<ColdWaterSource>,
    ) {
        let event_schedules = vec![
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            Some(vec![TypedScheduleEvent {
                start: 13.,
                duration: Some(2.),
                temperature: 38.,
                name: "other".into(),
                event_type: WaterScheduleEventType::Other,
                volume: None,
                warm_volume: None,
                pipework_volume: None,
            }]),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        ];

        let hot_water_sources: IndexMap<String, HotWaterSourceMockWithUniqueHotWaterTemperature> =
            IndexMap::from([(
                "hw cylinder".into(),
                HotWaterSourceMockWithUniqueHotWaterTemperature {},
            )]);

        let pre_heated_water_sources: IndexMap<
            String,
            HotWaterSourceMockWithUniqueHotWaterTemperature,
        > = IndexMap::from([]);

        let dhw_demand = create_dhw_demand(
            simulation_time,
            hot_water_sources,
            pre_heated_water_sources,
            cold_water_source,
            event_schedules,
        );

        // Test timestep 13 with Other event
        let actual = dhw_demand
            .hot_water_demand(SimulationTimeIteration {
                index: 13,
                time: 0.,
                timestep: 1.,
            })
            .unwrap();
        let (
            _hw_demand_volume,
            _hw_duration,
            _all_events,
            _hw_energy_demand,
            usage_events_with_flushes,
        ) = actual;

        let usage_events_with_flushes_for_hw_cylinder = usage_events_with_flushes
            .get::<String>(&"hw cylinder".into())
            .unwrap();

        // Verify that the Other event was processed
        assert_eq!(usage_events_with_flushes_for_hw_cylinder.len(), 2); // should have 1 event + 1 pipe flush
        assert_eq!(
            usage_events_with_flushes_for_hw_cylinder[0].event_result_type,
            WaterEventResultType::Other
        );
    }

    #[rstest]
    fn test_hot_water_demand_mixer_new_temp(
        simulation_time: SimulationTime,
        cold_water_source: Arc<ColdWaterSource>,
    ) {
        let event_schedules = vec![
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            Some(vec![TypedScheduleEvent {
                start: 14.,
                duration: Some(5.),
                temperature: 43.,
                name: "mixer".into(),
                event_type: WaterScheduleEventType::Shower,
                volume: None,
                warm_volume: None,
                pipework_volume: None,
            }]),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        ];

        let hot_water_sources: IndexMap<String, HotWaterSourceMockWithUniqueHotWaterTemperature> =
            IndexMap::from([(
                "hw cylinder".into(),
                HotWaterSourceMockWithUniqueHotWaterTemperature {},
            )]);

        let pre_heated_water_sources: IndexMap<
            String,
            HotWaterSourceMockWithUniqueHotWaterTemperature,
        > = IndexMap::from([]);

        let dhw_demand = create_dhw_demand(
            simulation_time,
            hot_water_sources,
            pre_heated_water_sources,
            cold_water_source,
            event_schedules,
        );

        // Test timestep 14 with mixer shower at new temperature
        let actual = dhw_demand
            .hot_water_demand(SimulationTimeIteration {
                index: 14,
                time: 0.,
                timestep: 1.,
            })
            .unwrap();
        let (
            _hw_demand_volume,
            _hw_duration,
            _all_events,
            _hw_energy_demand,
            usage_events_with_flushes,
        ) = actual;

        let usage_events_with_flushes_for_hw_cylinder = usage_events_with_flushes
            .get::<String>(&"hw cylinder".into())
            .unwrap();

        // Verify that the mixer shower event was processed
        assert_eq!(usage_events_with_flushes_for_hw_cylinder.len(), 2); // should have 1 event + 1 pipe flush
        assert_eq!(
            usage_events_with_flushes_for_hw_cylinder[0].event_result_type,
            WaterEventResultType::Shower
        );
    }

    #[rstest]
    fn test_hot_water_demand_multiple_mixer_same_temp(
        simulation_time: SimulationTime,
        cold_water_source: Arc<ColdWaterSource>,
    ) {
        let event_schedules = vec![
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            Some(vec![
                TypedScheduleEvent {
                    start: 13.,
                    duration: Some(3.),
                    temperature: 42.,
                    name: "mixer".into(),
                    event_type: WaterScheduleEventType::Shower,
                    volume: None,
                    warm_volume: None,
                    pipework_volume: None,
                },
                TypedScheduleEvent {
                    start: 13.2,
                    duration: Some(3.),
                    temperature: 42.,
                    name: "mixer".into(),
                    event_type: WaterScheduleEventType::Shower,
                    volume: None,
                    warm_volume: None,
                    pipework_volume: None,
                },
            ]),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        ];

        let hot_water_sources: IndexMap<String, HotWaterSourceMockWithUniqueHotWaterTemperature> =
            IndexMap::from([(
                "hw cylinder".into(),
                HotWaterSourceMockWithUniqueHotWaterTemperature {},
            )]);

        let pre_heated_water_sources: IndexMap<
            String,
            HotWaterSourceMockWithUniqueHotWaterTemperature,
        > = IndexMap::from([]);

        let dhw_demand = create_dhw_demand(
            simulation_time,
            hot_water_sources,
            pre_heated_water_sources,
            cold_water_source,
            event_schedules,
        );

        // Test timestep 13 with multiple mixer showers at same temperature
        let actual = dhw_demand
            .hot_water_demand(SimulationTimeIteration {
                index: 13,
                time: 0.,
                timestep: 1.,
            })
            .unwrap();
        let (
            _hw_demand_volume,
            _hw_duration,
            _all_events,
            _hw_energy_demand,
            usage_events_with_flushes,
        ) = actual;

        let usage_events_with_flushes_for_hw_cylinder = usage_events_with_flushes
            .get::<String>(&"hw cylinder".into())
            .unwrap();

        // Verify that both mixer shower events were processed
        assert_eq!(usage_events_with_flushes_for_hw_cylinder.len(), 4); // should have 2 event + 2 pipe flushes
        assert_eq!(
            usage_events_with_flushes_for_hw_cylinder[0].event_result_type,
            WaterEventResultType::Shower
        );
        assert_eq!(
            usage_events_with_flushes_for_hw_cylinder[1].event_result_type,
            WaterEventResultType::PipeFlush
        );
        assert_eq!(
            usage_events_with_flushes_for_hw_cylinder[2].event_result_type,
            WaterEventResultType::Shower
        );
        assert_eq!(
            usage_events_with_flushes_for_hw_cylinder[3].event_result_type,
            WaterEventResultType::PipeFlush
        );
    }

    // The following Python tests were skipped as they are not necessary to validate the Rust port:
    // test_hot_water_source_assignment
    // test_default_hot_water_source_assignment
    // test_missing_hot_water_source_assignment
    // test_hot_water_demand_invalid_events
    // test_init_distribution_pointofuse
    // test_init_distribution_missing_dict
    // test_init_distribution_part_missing_list
    // test_init_distribution_extra
    // test_calc_pipework_losses_with_invalid_location
    // test_init_invalid_shower_type
    // test_init_invalid_wwhrs_name_shower
    // test_init_invalid_energysupply_name_shower
    // test_init_missing_cold_water_source_shower
    // test_init_missing_cold_water_source_bath
    // test_init_missing_cold_water_source_other_water_use
    // test_create_water_distribution_pipework_no_tapping_points
}
