use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyConnection};
use crate::core::heating_systems::wwhrs::Wwhrs;
use crate::core::pipework::{PipeworkLocation, PipeworkSimple, Pipeworkesque};
use crate::core::schedule::{TypedScheduleEvent, WaterScheduleEventType};
use crate::core::units::{MILLIMETRES_IN_METRE, WATTS_PER_KILOWATT};
use crate::core::water_heat_demand::bath::Bath;
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::misc::{
    FRAC_DHW_ENERGY_INTERNAL_GAINS, WaterEventResult, WaterEventResultType, water_demand_to_kwh
};
use crate::core::water_heat_demand::other_hot_water_uses::OtherHotWater;
use crate::core::water_heat_demand::shower::Shower;
use crate::core::water_heat_demand::shower::{InstantElectricShower, MixerShower};
use crate::corpus::{ColdWaterSources, EventSchedule, HotWaterSource};
use crate::input::{
    BathDetails, Baths as BathInput, OtherWaterUse, OtherWaterUses as OtherWaterUseInput,
    PipeworkContents, Shower as ShowerInput, Showers as ShowersInput,
    WaterDistribution as WaterDistributionInput, WaterHeatingEvent,
    WaterPipeworkSimple,
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
pub struct DomesticHotWaterDemand {
    showers: HashMap<String, Shower>,
    baths: HashMap<String, Bath>,
    other: HashMap<String, OtherHotWater>,
    hot_water_sources: IndexMap<String, HotWaterSource>,
    energy_supply_conn_unmet_demand: IndexMap<String, EnergySupplyConnection>,
    source_supplying_outlet: HashMap<(OutletType, String), String>,
    hot_water_distribution_pipework:  IndexMap<String, Vec<PipeworkSimple>> ,
    event_schedules: EventSchedule,
    pre_heated_water_sources: IndexMap<String, HotWaterSource>,
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

impl DomesticHotWaterDemand {
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
        hot_water_sources: IndexMap<String, HotWaterSource>,
        pre_heated_water_sources: IndexMap<String, HotWaterSource>,
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
        let mut hw_pipework_inputs: IndexMap<String, Vec<WaterPipeworkSimple>> = match hw_pipework_inputs {
            WaterDistributionInput::List(pipeworks) => {
                if hot_water_sources.len() == 1 {
                    hot_water_sources.keys().map(|key| -> (String, Vec<WaterPipeworkSimple>) { (key.clone(), pipeworks.clone()) }).collect()
                } else {
                    bail!("If more than one HotWaterSource is defined, then distribution pipework must be defined for each one");
                }
            },
            WaterDistributionInput::Map(index_map) => {
                index_map.iter().map(|(key, value)| -> (String, Vec<WaterPipeworkSimple>) { (key.into(), value.clone()) }).collect()
            },
        };

        // pipework without a valid hot water source
        let pws_without_hws: Vec<_> = hw_pipework_inputs.keys().filter(|key| { !hot_water_sources.keys().contains(key) }).collect();
        if !pws_without_hws.is_empty() {
            // TODO include names in error message
            bail!("Distribution pipework defined for non-existent HotWaterSource(s)");
        }

        // hot water sources (not including point of use) without any pipework
        let hws_without_pws: Vec<_> = hot_water_sources.keys().filter(|key| { !hw_pipework_inputs.keys().contains(key) }).collect();
        for hws_name in hws_without_pws {
            let hot_water_source = hot_water_sources.get(hws_name).unwrap();
            match hot_water_source {
                HotWaterSource::PointOfUse(_) => {
                    // point of use doesn't need pipework - just add an empty vec
                    hw_pipework_inputs.insert(hws_name.clone(), vec![]);
                },
                _ => {
                    // TODO include name in error message
                    bail!("Distribution pipework not specified for HotWaterSource");
                },
            }
        }

        for (hws_name, hws) in &hot_water_sources {
            match hws {
                HotWaterSource::PointOfUse(_) => continue,
                _ => {
                    for data in hw_pipework_inputs.get(hws_name).unwrap() {
                        let pipework = input_to_water_distribution_pipework(data, total_number_tapping_points)?;
                        let entry = hot_water_distribution_pipework.get_mut(hws_name).unwrap();
                        entry.push(pipework);
                    }
                }
            }
        }

        // Set up unmet demand connection for each hot water source
        let energy_supply_conn_unmet_demand: IndexMap<String, EnergySupplyConnection> =
            hot_water_sources
                .keys()
                .map(|name| {
                    let energy_supply = energy_supplies.get("unmet_demand").unwrap();
                    let energy_suppy_conn_unmet_demand =
                        EnergySupply::connection(energy_supply.clone(), &name).unwrap(); // TODO avoid unwrap here
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
        hot_water_sources: &IndexMap<String, HotWaterSource>,
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
                        (OutletType::Shower, bath_name.into()),
                        hot_water_source.clone(),
                    );
                }
                None => {
                    if hot_water_sources.len() == 1 {
                        let default_hot_water_source = hot_water_sources.first().unwrap().0;
                        mapping.insert(
                            (OutletType::Shower, bath_name.into()),
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
                        (OutletType::Shower, other_name.into()),
                        hot_water_source.clone(),
                    );
                }
                None => {
                    if hot_water_sources.len() == 1 {
                        let default_hot_water_source = hot_water_sources.first().unwrap().0;
                        mapping.insert(
                            (OutletType::Shower, other_name.into()),
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
        hot_water_source: HotWaterSource,
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

    pub fn hot_water_demand<'a>(
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
                    Option<&HotWaterSource>,
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
                                    volume_required.clone(),
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

                let cold_water_temperature = if hw_demand_i.is_some() {
                    let hw_demand_i = hw_demand_i.unwrap();
                    *hw_demand_volume.get_mut(&hot_water_source_name).unwrap() += hw_demand_i;
                    let list_temperature_volume = cold_water_source
                        .get_temp_cold_water(hw_demand_target_i - hw_demand_i, simtime)?;
                    let sum_t_by_v: f64 = list_temperature_volume.iter().map(|(t, v)| t * v).sum();
                    let sum_v: f64 = list_temperature_volume.iter().map(|(t, v)| v).sum();
                    sum_t_by_v / sum_v
                } else {
                    let list_temperature_volume =
                        cold_water_source.get_temp_cold_water(hw_demand_target_i, simtime)?;
                    let sum_t_by_v: f64 = list_temperature_volume.iter().map(|(t, v)| t * v).sum();
                    let sum_v: f64 = list_temperature_volume.iter().map(|(t, v)| v).sum();
                    sum_t_by_v / sum_v
                };

                let hw_energy_demand_i = water_demand_to_kwh(
                    hw_demand_target_i,
                    event.temperature,
                    cold_water_temperature,
                );
                *hw_duration.get_mut(&hot_water_source_name).unwrap() += event.duration.unwrap();
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
                    .push(event.into());

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
                        })
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
            match hot_water_source {
                HotWaterSource::PreHeated(storage_tank) => {
                    storage_tank.demand_hot_water(None, simtime);
                }
                // TODO is this true? can we just always call demand_hot_water regardless?
                _ => {
                    bail!("Preheated hot water sources must be storage tanks");
                }
            }
        }

        let mut hw_energy_demand_at_hot_water_source: IndexMap<String, f64> = Default::default();
        let mut hw_energy_output: IndexMap<String, f64> = Default::default();
        let mut pw_losses_total: IndexMap<String, f64> = Default::default();
        let mut gains_internal_dhw: IndexMap<String, f64> = Default::default();
        let mut primary_pw_losses: IndexMap<String, f64> = Default::default();
        let mut storage_losses: IndexMap<String, f64> = Default::default();

        let mut all_keys: Vec<String> = self
            .hot_water_sources
            .keys()
            .into_iter()
            .map(|x| x.clone())
            .collect();
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
                .map(|x| *x)
                .collect();

            // TODO update demand_hot_water to accept usage_events
            hw_energy_output.insert(hws_name.clone(), hws.demand_hot_water(filtered_events.clone(), simtime)?);

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
            if internal_gains.is_some() {
                *gains_internal_dhw.get_mut(hws_name).unwrap() += internal_gains.unwrap();
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
        no_of_hw_events: usize,
        demand_water_temperature: f64,
        internal_air_temperature: f64,
        external_air_temperature: f64,
    ) -> (f64, f64) {
        if self.hot_water_distribution_pipework.is_empty() {
            return (0., 0.);
        }

        let mut cool_down_loss_internal = 0.0;
        let mut cool_down_loss_external = 0.0;

        for pipework in self.hot_water_distribution_pipework.get(hot_water_source_name).unwrap() {
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
                    let (pw_losses_internal_i, pw_losses_external_i) = self.calc_pipework_losses(&hot_water_source_name, 1, event.temperature_warm, internal_air_temperature, external_air_temperature);
                    pw_losses_internal += pw_losses_internal_i;
                    pw_losses_external += pw_losses_external_i;
                },
                _ => {
                    gains_internal_dhw_use += FRAC_DHW_ENERGY_INTERNAL_GAINS * water_demand_to_kwh(event.volume_warm, event.temperature_warm, internal_air_temperature);
                }
            }
        }

        (pw_losses_internal, pw_losses_external, gains_internal_dhw_use)
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
            let energy_supply_conn = EnergySupply::connection(energy_supply, name).unwrap();

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

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0., 24., 1.)
    }

    #[fixture]
    fn dhw_demand(simulation_time: SimulationTime) -> DomesticHotWaterDemand {
        let cold_water_temps = vec![
            2.0, 3.0, 4.0, 2.0, 3.0, 4.0, 2.0, 3.0, 4.0, 2.0, 3.0, 4.0, 2.0, 3.0, 4.0, 2.0, 3.0,
            4.0, 2.0, 3.0, 4.0, 2.0, 3.0, 4.0,
        ];
        let cold_water_source = Arc::from(ColdWaterSource::new(cold_water_temps, 0, 1.));
        let cold_water_sources =
            ColdWaterSources::from([("mains water".into(), cold_water_source.clone())]);
        let flow_rates = vec![5., 7., 9., 11., 13.];
        let efficiencies = vec![44.8, 39.1, 34.8, 31.4, 28.6];
        let utilisation_factor = 0.7;
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
        let energy_supplies = IndexMap::from([("mains elec".into(), electricity_supply.clone())]);

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

        let hw_pipework = Some(WaterDistribution::List(vec![
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
        ]));

        let event_schedules = vec![
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
                    duration: None,
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
        ];

        let hot_water_sources = todo!();
        let pre_heated_water_sources = todo!();

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
    #[ignore = "not yet implemented for 1_0_a1"]
    fn test_hot_water_demand(dhw_demand: DomesticHotWaterDemand, simulation_time: SimulationTime) {
        for (t_idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                dhw_demand.hot_water_demand(t_it, 55.0).unwrap(),
                [
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 29.783791734078537,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: Some(vec![
                            TypedScheduleEvent {
                                start: 4.1,
                                duration: Some(6.),
                                temperature: 41.0,
                                name: "IES".to_string(),
                                event_type: WaterScheduleEventType::Shower,
                                volume: None,
                                warm_volume: None,
                                pipework_volume: Some(7.556577529434648),
                            },
                            TypedScheduleEvent {
                                start: 4.5,
                                duration: Some(6.),
                                temperature: 41.0,
                                name: "IES".to_string(),
                                event_type: WaterScheduleEventType::Shower,
                                volume: None,
                                warm_volume: None,
                                pipework_volume: Some(7.556577529434648),
                            }
                        ]),
                        vol_hot_water_equiv_elec_shower: 29.783791734078537
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 81.141483189812,
                        hw_demand_vol_target: IndexMap::from([
                            (
                                41.0.into(),
                                VolumeReference {
                                    warm_temp: 41.0,
                                    warm_vol: 100.
                                }
                            ),
                            (
                                DemandVolTargetKey::TempHotWater,
                                VolumeReference {
                                    warm_temp: 55.0,
                                    warm_vol: 81.141483189812
                                }
                            )
                        ]),
                        hw_vol_at_tapping_points: 88.195822360114,
                        hw_duration: 12.5,
                        all_events: 1usize,
                        hw_energy_demand: 4.532666666666667,
                        usage_events: Some(vec![
                            TypedScheduleEvent {
                                start: 6.,
                                duration: Some(6.),
                                temperature: 41.0,
                                name: "IES".to_string(),
                                event_type: WaterScheduleEventType::Shower,
                                volume: None,
                                warm_volume: None,
                                pipework_volume: Some(7.556577529434648),
                            },
                            TypedScheduleEvent {
                                start: 6.,
                                duration: Some(12.5),
                                temperature: 41.0,
                                name: "medium".to_string(),
                                event_type: WaterScheduleEventType::Bath,
                                volume: Some(100.0),
                                warm_volume: Some(100.),
                                pipework_volume: Some(7.556577529434648),
                            }
                        ]),
                        vol_hot_water_equiv_elec_shower: 14.610916699736642
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 53.58989404060333,
                        hw_demand_vol_target: IndexMap::from([
                            (
                                41.0.into(),
                                VolumeReference {
                                    warm_temp: 41.0,
                                    warm_vol: 56.0
                                }
                            ),
                            (
                                DemandVolTargetKey::TempHotWater,
                                VolumeReference {
                                    warm_temp: 55.0,
                                    warm_vol: 53.58989404060333
                                }
                            )
                        ]),
                        hw_vol_at_tapping_points: 38.47673898173404,
                        hw_duration: 7.0,
                        all_events: 2,
                        hw_energy_demand: 2.325363096327197,
                        usage_events: Some(vec![
                            TypedScheduleEvent {
                                start: 7.,
                                duration: Some(6.),
                                temperature: 41.0,
                                name: "mixer".to_string(),
                                event_type: WaterScheduleEventType::Shower,
                                volume: None,
                                warm_volume: Some(48.0),
                                pipework_volume: Some(7.556577529434648),
                            },
                            TypedScheduleEvent {
                                start: 7.,
                                duration: Some(1.),
                                temperature: 41.0,
                                name: "other".to_string(),
                                event_type: WaterScheduleEventType::Other,
                                volume: None,
                                warm_volume: Some(8.0),
                                pipework_volume: Some(7.556577529434648),
                            }
                        ]),
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 64.89041357202146,
                        hw_demand_vol_target: IndexMap::from([
                            (
                                41.0.into(),
                                VolumeReference {
                                    warm_temp: 41.0,
                                    warm_vol: 72.0
                                }
                            ),
                            (
                                DemandVolTargetKey::TempHotWater,
                                VolumeReference {
                                    warm_temp: 55.0,
                                    warm_vol: 64.89041357202146
                                }
                            ),
                        ]),
                        hw_vol_at_tapping_points: 49.77725851315216,
                        hw_duration: 9.0,
                        all_events: 2,
                        hw_energy_demand: 2.9504640362695724,
                        usage_events: Some(vec![
                            TypedScheduleEvent {
                                start: 8.,
                                duration: Some(6.),
                                temperature: 41.0,
                                name: "mixer".to_string(),
                                event_type: WaterScheduleEventType::Shower,
                                volume: None,
                                warm_volume: Some(48.0),
                                pipework_volume: Some(7.556577529434648),
                            },
                            TypedScheduleEvent {
                                start: 8.,
                                duration: Some(3.),
                                temperature: 41.0,
                                name: "medium".to_string(),
                                event_type: WaterScheduleEventType::Bath,
                                volume: None,
                                warm_volume: Some(24.),
                                pipework_volume: Some(7.556577529434648),
                            }
                        ]),
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    },
                    DomesticHotWaterDemandData {
                        hw_demand_vol: 0.0,
                        hw_demand_vol_target: Default::default(),
                        hw_vol_at_tapping_points: 0.0,
                        hw_duration: 0.0,
                        all_events: Default::default(),
                        hw_energy_demand: 0.0,
                        usage_events: None,
                        vol_hot_water_equiv_elec_shower: 0.0
                    }
                ][t_idx]
            );
        }
    }

    #[rstest]
    #[ignore = "not yet implemented for 1_0_a1"]
    fn test_calc_pipework_losses(
        dhw_demand: DomesticHotWaterDemand,
        simulation_time: SimulationTime,
    ) {
        for (t_idx, _) in simulation_time.iter().enumerate() {
            assert_eq!(
                dhw_demand.calc_pipework_losses(
                    1.,
                    0.,
                    [0, 0, 0, 2, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0][t_idx],
                    [
                        55.0, 55.0, 55.0, 57.0, 55.0, 58.0, 60.0, 40.0, 55.0, 55.0, 55.0, 55.0,
                        55.0, 55.0, 55.0, 55.0, 55.0, 55.0, 55.0, 55.0, 55.0, 55.0, 55.0, 55.0
                    ][t_idx],
                    20.0,
                    5.0
                ),
                [
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
                    (0.0, 0.0)
                ][t_idx]
            );
        }
    }

    #[rstest]
    #[ignore = "not yet implemented for 1_0_a1"]
    fn test_calc_pipework_losses_with_no_pipework(
        mut dhw_demand: DomesticHotWaterDemand,
        simulation_time: SimulationTime,
    ) {
        dhw_demand.hot_water_distribution_pipework = vec![];

        for _ in simulation_time.iter() {
            assert_eq!(
                dhw_demand.calc_pipework_losses(1., 0., 0, 55., 20., 5.),
                (0., 0.)
            );
        }
    }
}
