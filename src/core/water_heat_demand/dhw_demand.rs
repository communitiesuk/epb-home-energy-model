use crate::core::energy_supply::energy_supply::EnergySupply;
use crate::core::heating_systems::wwhrs::Wwhrs;
use crate::core::pipework::{PipeworkLocation, PipeworkSimple, Pipeworkesque};
use crate::core::schedule::{TypedScheduleEvent, WaterScheduleEventType};
use crate::core::units::MILLIMETRES_IN_METRE;
use crate::core::water_heat_demand::bath::Bath;
use crate::core::water_heat_demand::misc::water_demand_to_kwh;
use crate::core::water_heat_demand::other_hot_water_uses::OtherHotWater;
use crate::core::water_heat_demand::shower::Shower;
use crate::core::water_heat_demand::shower::{InstantElectricShower, MixerShower};
use crate::corpus::{ColdWaterSources, EventSchedule};
use crate::input::{
    BathDetails, Baths as BathInput, OtherWaterUse, OtherWaterUses as OtherWaterUseInput,
    PipeworkContents, Shower as ShowerInput, Showers as ShowersInput,
    WaterDistribution as WaterDistributionInput, WaterDistribution, WaterPipeworkSimple,
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

#[derive(Debug)]
pub struct DomesticHotWaterDemand {
    showers: HashMap<String, Shower>,
    baths: HashMap<String, Bath>,
    other: HashMap<String, OtherHotWater>,
    hot_water_distribution_pipework: Vec<PipeworkSimple>,
    event_schedules: EventSchedule,
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

impl DomesticHotWaterDemand {
    // TODO (from Python) Enhance analysis for overlapping events
    // Part of draft code for future overlapping analysis of events
    // For pipework losses count only none overlapping events
    // Time of finalisation of the previous hot water event
    // __time_end_previous_event = 0.0

    pub fn new(
        showers_input: ShowersInput,
        bath_input: BathInput,
        other_hot_water_input: OtherWaterUseInput,
        water_distribution_input: Option<WaterDistributionInput>,
        cold_water_sources: &ColdWaterSources,
        wwhrs: &IndexMap<String, Arc<Mutex<Wwhrs>>>,
        energy_supplies: &IndexMap<String, Arc<RwLock<EnergySupply>>>,
        event_schedules: EventSchedule,
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

        let hot_water_distribution_pipework = water_distribution_input
            .iter()
            .flat_map(|input| {
                // TODO: revise this to reflect upstream logic - below is just to make compiler happy for now after changing input for 1.0.0a1
                let simple_pipework = match input {
                    WaterDistribution::List(simple_pipework) => simple_pipework.clone(),
                    WaterDistribution::Map(map_of_simple_pipework) => map_of_simple_pipework
                        .values()
                        .flatten()
                        .cloned()
                        .collect_vec(),
                };
                simple_pipework
                    .iter()
                    .map(|input| {
                        input_to_water_distribution_pipework(input, total_number_tapping_points)
                    })
                    .collect_vec()
            })
            .collect::<anyhow::Result<Vec<PipeworkSimple>>>()?;

        Ok(Self {
            showers,
            baths,
            other,
            hot_water_distribution_pipework,
            event_schedules,
        })
    }

    pub fn hot_water_demand(
        &self,
        simtime: SimulationTimeIteration,
        temp_hot_water: f64,
    ) -> anyhow::Result<DomesticHotWaterDemandData> {
        let mut hw_demand_vol = 0.;
        let mut hw_demand_vol_target: IndexMap<DemandVolTargetKey, VolumeReference> =
            Default::default();
        let mut hw_energy_demand = 0.;
        let mut hw_duration = 0.;
        let mut all_events = 0usize;
        let mut vol_hot_water_equiv_elec_shower = 0.;

        // Events have been organised now so that they are structured by simple step t_idx and
        // sorted for each time step from start to end.
        //
        // TODO (from Python) No overlapping is currently considered in terms of interaction with hot water
        // source (tank). The first event that starts is Served before the second event
        // is considered even if this starts before the previous event has finished.

        let mut usage_events = self.event_schedules[simtime.index].clone();

        if let Some(usage_events) = &mut usage_events {
            for event in usage_events.iter_mut() {
                match event.event_type {
                    WaterScheduleEventType::Shower => {
                        for (name, shower) in self.showers.iter() {
                            if name != &event.name {
                                continue;
                            }
                            // If shower is used in the current timestep, get details of use
                            // and calculate HW demand from shower
                            let the_cold_water_temp = shower.get_cold_water_source();
                            let cold_water_temperature = the_cold_water_temp.temperature(simtime);

                            let shower_temp = event.temperature;
                            let label_temp = shower_temp.into();
                            let shower_duration = event
                                .duration
                                .expect("A duration is expected to be defined for a shower event.");

                            let (hw_demand_i, hw_demand_target_i) = shower.hot_water_demand(
                                shower_temp,
                                temp_hot_water,
                                shower_duration,
                                simtime,
                            );

                            if let Shower::InstantElectricShower(_) = shower {
                                vol_hot_water_equiv_elec_shower += hw_demand_i;
                            } else {
                                event.warm_volume = Some(hw_demand_target_i);
                                // don't add hw demand and pipework loss from electric shower
                                hw_demand_vol += hw_demand_i;
                                hw_energy_demand += water_demand_to_kwh(
                                    hw_demand_i,
                                    temp_hot_water,
                                    cold_water_temperature,
                                );
                                hw_duration += shower_duration;
                                all_events += 1;

                                hw_demand_vol_target
                                    .entry(label_temp)
                                    .and_modify(|vol| vol.warm_vol += hw_demand_target_i)
                                    .or_insert(VolumeReference {
                                        warm_temp: shower_temp,
                                        warm_vol: hw_demand_target_i,
                                    });
                            }
                        }
                    }
                    WaterScheduleEventType::Other => {
                        for (name, other) in self.other.iter() {
                            if name != &event.name {
                                continue;
                            }
                            // If other is used in the current timestep, get details of use
                            // and calculate HW demand from other
                            let the_cold_water_temp = other.get_cold_water_source();
                            let cold_water_temperature = the_cold_water_temp.temperature(simtime);

                            let other_temp = event.temperature;
                            let label_temp = other_temp.into();
                            let other_duration = event
                                .duration
                                .expect("A duration is expected for an 'other' water use event.");
                            let (hw_demand_i, hw_demand_target_i) = other.hot_water_demand(
                                other_temp,
                                temp_hot_water,
                                other_duration,
                                simtime,
                            );
                            event.warm_volume = Some(hw_demand_target_i);
                            hw_demand_vol_target
                                .entry(label_temp)
                                .and_modify(|vol| vol.warm_vol += hw_demand_target_i)
                                .or_insert(VolumeReference {
                                    warm_temp: other_temp,
                                    warm_vol: hw_demand_target_i,
                                });

                            hw_demand_vol += hw_demand_i;

                            // Check if it makes sense to call again the hot_water_demand function instead of sending hw_demand_i previously calculated
                            hw_energy_demand += water_demand_to_kwh(
                                hw_demand_i,
                                temp_hot_water,
                                cold_water_temperature,
                            );
                            hw_duration += other_duration;
                            all_events += 1;
                        }
                    }
                    WaterScheduleEventType::Bath => {
                        for (name, bath) in self.baths.iter() {
                            if name != &event.name {
                                continue;
                            }
                            // If bath is used in the current timestep, get details of use
                            // and calculate HW demand from bath
                            let the_cold_water_temp = bath.get_cold_water_source();
                            let cold_water_temperature = the_cold_water_temp.temperature(simtime);

                            // Assume flow rate for bath event is the same as other hot water events
                            let peak_flowrate = bath.get_flowrate();
                            // litres bath  / litres per minute flowrate = minutes
                            let (bath_volume, bath_duration) = if let Some(volume) = event.volume {
                                let bath_volume = volume;
                                let bath_duration = volume / peak_flowrate;
                                event.duration.replace(bath_duration);
                                (bath_volume, bath_duration)
                            } else if let Some(duration) = event.duration {
                                let bath_duration = duration;
                                let bath_volume = bath_duration * peak_flowrate;
                                (bath_volume, bath_duration)
                            } else {
                                bail!("Water event '{name}' has no volume or duration defined.");
                            };
                            let bath_temp = event.temperature;
                            let label_temp = bath_temp.into();
                            let (hw_demand_i, hw_demand_target_i) = bath.hot_water_demand(
                                bath_temp,
                                temp_hot_water,
                                bath_volume,
                                simtime,
                            );
                            event.warm_volume = Some(hw_demand_target_i);
                            hw_demand_vol_target
                                .entry(label_temp)
                                .and_modify(|vol| vol.warm_vol += hw_demand_target_i)
                                .or_insert(VolumeReference {
                                    warm_temp: bath_temp,
                                    warm_vol: hw_demand_target_i,
                                });

                            hw_demand_vol += hw_demand_i;
                            // Check if it makes sense to call again the hot_water_demand function instead of sending hw_demand_i previously calculated
                            // bath.hot_water_demand(bath_temp, temp_hot_water)[0],
                            hw_energy_demand += water_demand_to_kwh(
                                hw_demand_i,
                                temp_hot_water,
                                cold_water_temperature,
                            );
                            hw_duration += bath_duration;
                            all_events += 1;
                        }
                    }
                }
            }
        }

        let hw_vol_at_tapping_points = hw_demand_vol + vol_hot_water_equiv_elec_shower;

        let vol_hot_water_left_in_pipework = self
            .hot_water_distribution_pipework
            .iter()
            .map(|pipework| pipework.volume())
            .sum::<f64>();
        hw_demand_vol += all_events as f64 * vol_hot_water_left_in_pipework;

        // TODO (from Python) Refine pipework losses by considering overlapping of events
        //                    and shared pipework between serving tap points
        //                    none_overlapping_events calculated above is a lower bound(ish)
        //                    approximation for this
        if let Some(usage_events) = &mut usage_events {
            for event in usage_events.iter_mut() {
                event.pipework_volume = Some(vol_hot_water_left_in_pipework);
            }
        }

        if hw_demand_vol > 0.0 {
            hw_demand_vol_target.insert(
                DemandVolTargetKey::TempHotWater,
                VolumeReference {
                    warm_temp: temp_hot_water,
                    warm_vol: hw_demand_vol,
                },
            );
        }

        // Return:
        // - litres hot water per timestep (demand on hw system)
        // - map containing litres of warm water required at different temperature levels
        // - litres hot water per timestep (output at tapping points)
        // - minutes demand per timestep,
        // - number of events in timestep
        // - hot water energy demand (kWh)
        // - usage_events updated to reflect pipework volumes and bath durations
        Ok(DomesticHotWaterDemandData {
            hw_demand_vol,
            hw_demand_vol_target,
            hw_vol_at_tapping_points,
            hw_duration,
            all_events,
            hw_energy_demand,
            usage_events,
            vol_hot_water_equiv_elec_shower,
        })
    }

    pub fn calc_pipework_losses(
        &self,
        _delta_t_h: f64,
        _hw_duration: f64,
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

        for pipework in &self.hot_water_distribution_pipework {
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
}

#[derive(Clone, Debug, PartialEq)]
pub struct DomesticHotWaterDemandData {
    pub hw_demand_vol: f64,
    pub(crate) hw_demand_vol_target: IndexMap<DemandVolTargetKey, VolumeReference>,
    pub(crate) hw_vol_at_tapping_points: f64,
    pub(crate) hw_duration: f64,
    pub(crate) all_events: usize,
    pub(crate) hw_energy_demand: f64,
    pub(crate) usage_events: Option<Vec<TypedScheduleEvent>>,
    pub(crate) vol_hot_water_equiv_elec_shower: f64,
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
            ..
        } => {
            let cold_water_source = cold_water_sources.get(cold_water_source).unwrap().clone();
            let wwhrs_instance: Option<Arc<Mutex<Wwhrs>>> = wwhrs_config
                .as_ref()
                .map(|config| &config.waste_water_heat_recovery_system)
                .and_then(|w| wwhrs.get(w).cloned());

            Shower::MixerShower(MixerShower::new(
                *flowrate,
                cold_water_source,
                wwhrs_instance,
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
                    volume_hot: None,
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
                    volume_hot: None,
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
                    volume_hot: None,
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
                    volume_hot: None,
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
                    volume_hot: None,
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
                    volume_hot: None,
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
                    volume_hot: None,
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
                    volume_hot: None,
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

        DomesticHotWaterDemand::new(
            showers_input,
            baths_input,
            other_input,
            hw_pipework,
            &cold_water_sources,
            &wwhrs,
            &energy_supplies,
            event_schedules,
        )
        .unwrap()
    }

    #[rstest]
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
                                volume_hot: None
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
                                volume_hot: None
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
                                volume_hot: None
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
                                volume_hot: None
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
                                volume_hot: None
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
                                volume_hot: None
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
                                volume_hot: None
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
                                volume_hot: None
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
