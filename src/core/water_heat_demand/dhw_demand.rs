use indexmap::IndexMap;

use crate::core::energy_supply::energy_supply::{EnergySupplies, EnergySupply};
use crate::core::heating_systems::wwhrs::Wwhrs;
use crate::core::pipework::Pipework;
use crate::core::schedule::ScheduleEvent;
use crate::core::units::{MILLIMETRES_IN_METRE, MINUTES_PER_HOUR};
use crate::core::water_heat_demand::bath::Bath;
use crate::core::water_heat_demand::misc::water_demand_to_kwh;
use crate::core::water_heat_demand::other_hot_water_uses::OtherHotWater;
use crate::core::water_heat_demand::shower::Shower;
use crate::core::water_heat_demand::shower::{InstantElectricShower, MixerShower};
use crate::corpus::{ColdWaterSources, HotWaterEventSchedules};
use crate::input::{
    Bath as BathInput, BathDetails, ColdWaterSourceType, EnergySupplyType,
    InstantElectricShower as InstantElectricShowerInput, MixerShower as MixerShowerInput,
    OtherWaterUse as OtherWaterUseInput, OtherWaterUseDetails, Shower as ShowerInput,
    WaterDistribution as WaterDistributionInput, WaterPipework,
};
use std::collections::HashMap;

pub struct DomesticHotWaterDemand {
    showers: HashMap<String, Shower>,
    baths: HashMap<String, Bath>,
    other: HashMap<String, OtherHotWater>,
    hot_water_distribution_pipework: HashMap<String, Pipework>,
    event_schedules: HotWaterEventSchedules,
}

impl DomesticHotWaterDemand {
    pub fn new(
        shower_input: Option<ShowerInput>,
        bath_input: Option<BathInput>,
        other_hot_water_input: Option<OtherWaterUseInput>,
        water_distribution_input: Option<WaterDistributionInput>,
        cold_water_sources: &ColdWaterSources,
        wwhrs: &IndexMap<String, Wwhrs>,
        energy_supplies: &EnergySupplies,
        event_schedules: HotWaterEventSchedules,
    ) -> Self {
        let showers: HashMap<String, Shower> = shower_input
            .iter()
            .flat_map(|input| {
                let mut showers = vec![(
                    "mixer".to_owned(),
                    mixer_shower_input_to_shower(&input.mixer, cold_water_sources, wwhrs),
                )];
                if let Some(ies) = &input.ies {
                    showers.push((
                        "IES".to_owned(),
                        instant_electric_shower_input_to_shower(
                            "IES".to_owned(),
                            ies,
                            cold_water_sources,
                            energy_supplies,
                        ),
                    ))
                }

                showers
            })
            .collect();
        let baths: HashMap<String, Bath> = bath_input
            .iter()
            .flat_map(|input| match &input.medium {
                Some(details) => vec![(
                    "medium".to_owned(),
                    input_to_bath(details, cold_water_sources),
                )],
                None => vec![],
            })
            .collect();
        let other: HashMap<String, OtherHotWater> = other_hot_water_input
            .iter()
            .flat_map(|input| match &input.other {
                Some(details) => vec![(
                    "other".to_owned(),
                    input_to_other_water_events(details, cold_water_sources),
                )],
                None => vec![],
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
                vec![
                    (
                        "internal".to_owned(),
                        input_to_water_distribution_pipework(
                            &input.internal,
                            total_number_tapping_points,
                        ),
                    ),
                    (
                        "external".to_owned(),
                        input_to_water_distribution_pipework(
                            &input.external,
                            total_number_tapping_points,
                        ),
                    ),
                ]
            })
            .collect();

        Self {
            showers,
            baths,
            other,
            hot_water_distribution_pipework,
            event_schedules,
        }
    }

    pub fn hot_water_demand(&mut self, timestep_idx: usize) -> (f64, f64, f64, usize, f64) {
        let mut hw_demand_vol = 0.;
        let mut hw_energy_demand = 0.;
        let mut hw_duration = 0.;
        let mut all_events = 0;
        let mut vol_hot_water_equiv_elec_shower = 0.;

        for (name, shower) in self.showers.iter_mut() {
            let usage_events: Option<Vec<ScheduleEvent>> = self
                .event_schedules
                .shower
                .get(name)
                .and_then(|schedule| schedule[timestep_idx].clone());

            let cold_water_source = shower.get_cold_water_source();
            let cold_water_temperature = cold_water_source.temperature(timestep_idx);

            for event in usage_events.iter().flatten() {
                let shower_temp = event
                    .temperature
                    .expect("This usage event is expected to have an associated temperature.");
                let shower_duration = event
                    .duration
                    .expect("This usage event is expected to have an associated duration.");
                let hw_demand_i =
                    shower.hot_water_demand(shower_temp, shower_duration, timestep_idx);
                match &shower {
                    Shower::InstantElectricShower(_) => {
                        vol_hot_water_equiv_elec_shower += hw_demand_i;
                    }
                    // NB. the Python code here assumes the below is any shower that is not instant electric, rather than specifically mixer
                    Shower::MixerShower(s) => {
                        // don't add hw demand and pipework loss from electric shower
                        hw_demand_vol += hw_demand_i;
                        hw_energy_demand += water_demand_to_kwh(
                            hw_demand_i,
                            s.get_temp_hot(),
                            cold_water_temperature,
                        );
                        hw_duration += shower_duration;
                        all_events += 1;
                    }
                }
            }
        }

        for (name, other) in &self.other {
            // Get all other use events for the current timestep
            let usage_events: Option<Vec<ScheduleEvent>> = self
                .event_schedules
                .other
                .get(name)
                .and_then(|schedule| schedule[timestep_idx].clone());
            let cold_water_source = other.get_cold_water_source();
            let cold_water_temperature = cold_water_source.temperature(timestep_idx);

            // If other is used in the current timestep, get details of use
            // and calculate HW demand from other
            for event in usage_events.iter().flatten() {
                let other_temp = event
                    .temperature
                    .expect("This usage event is expected to have an associated temperature.");
                let other_duration = event
                    .duration
                    .expect("This usage event is expected to have an associated duration.");
                hw_demand_vol += other.hot_water_demand(other_temp, other_duration, timestep_idx);
                hw_energy_demand += water_demand_to_kwh(
                    other.hot_water_demand(other_temp, other_duration, timestep_idx),
                    other.get_temp_hot(),
                    cold_water_temperature,
                );
                hw_duration += other_duration;
                all_events += 1;
            }
        }

        for (name, bath) in &self.baths {
            let usage_events: Option<Vec<ScheduleEvent>> = self
                .event_schedules
                .bath
                .get(name)
                .and_then(|schedule| schedule[timestep_idx].clone());
            let cold_water_source = bath.get_cold_water_source();
            let cold_water_temperature = cold_water_source.temperature(timestep_idx);

            // Assume flow rate for bath event is the same as other hot water events
            let peak_flowrate = bath.get_flowrate();

            // If bath is used in the current timestep, get details of use
            // and calculate HW demand from bath
            // Note that bath size is the total water used per bath, not the total capacity of the bath
            for event in usage_events.iter().flatten() {
                let bath_temp = event
                    .temperature
                    .expect("This usage event is expected to have an associated temperature.");
                hw_demand_vol += bath.hot_water_demand(bath_temp, timestep_idx);
                let bath_duration = bath.get_size() / peak_flowrate;
                hw_energy_demand += water_demand_to_kwh(
                    bath.hot_water_demand(bath_temp, timestep_idx),
                    bath.get_temp_hot(),
                    cold_water_temperature,
                );
                hw_duration += bath_duration;
                all_events += 1;
            }
        }

        let hw_vol_at_tapping_points = hw_demand_vol + vol_hot_water_equiv_elec_shower;

        if !self.hot_water_distribution_pipework.is_empty() {
            let vol_hot_water_left_in_pipework = self
                .hot_water_distribution_pipework
                .get("internal")
                .expect("Internal pipework is expected to be defined at this point.")
                .volume_in_litres()
                + self
                    .hot_water_distribution_pipework
                    .get("external")
                    .expect("External pipework is expected to be defined at this point.")
                    .volume_in_litres();
            hw_demand_vol += all_events as f64 * vol_hot_water_left_in_pipework;
        }

        // Return:
        // - litres hot water per timestep (demand on hw system)
        // - litres hot water per timestep (output at tapping points)
        // - minutes demand per timestep,
        // - number of events in timestep
        // - hot water energy demand (kWh)
        (
            hw_demand_vol,
            hw_vol_at_tapping_points,
            hw_duration,
            all_events,
            hw_energy_demand,
        )
    }

    pub fn calc_pipework_losses(
        &self,
        delta_t_h: f64,
        hw_duration: f64,
        no_of_hw_events: usize,
        demand_water_temperature: f64,
        internal_air_temperature: f64,
        external_air_temperature: f64,
    ) -> (f64, f64) {
        if self.hot_water_distribution_pipework.is_empty() {
            return (0., 0.);
        }

        let mut _hot_water_time_fraction = hw_duration / (delta_t_h * MINUTES_PER_HOUR as f64);
        if _hot_water_time_fraction > 1. {
            _hot_water_time_fraction = 1.;
        }

        // TODO For now, ignore heat loss from pipes while water is flowing, as
        //              this is is not currently added to hot water demand, but is added
        //              to the internal gains. This would mean that reducing insulation
        //              would reduce overall energy demand, which would not be correct.
        //         pipework_watts_heat_loss_internal = self.__hw_distribution_pipework["internal"].heat_loss(
        //             demand_water_temperature,
        //             internal_air_temperature,
        //             )
        //         pipework_watts_heat_loss_external = self.__hw_distribution_pipework["external"].heat_loss(
        //             demand_water_temperature,
        //             external_air_temperature,
        //             )
        //
        //         # only calculate loss for times when there is hot water in the pipes - multiply by time fraction to get to kWh
        //         pipework_heat_loss_internal \
        //             = pipework_watts_heat_loss_internal \
        //             * hot_water_time_fraction \
        //             * delta_t_h \
        //             / units.W_per_kW # convert to kWh
        //         pipework_heat_loss_external \
        //             = pipework_watts_heat_loss_external \
        //             * hot_water_time_fraction \
        //             * delta_t_h \
        //             / units.W_per_kW # convert to kWh
        let pipework_heat_loss_internal = no_of_hw_events as f64
            * self
                .hot_water_distribution_pipework
                .get("internal")
                .expect("expected internal pipework to be defined")
                .cool_down_loss(demand_water_temperature, internal_air_temperature);
        let pipework_heat_loss_external = no_of_hw_events as f64
            * self
                .hot_water_distribution_pipework
                .get("external")
                .expect("expected external pipework to be defined")
                .cool_down_loss(demand_water_temperature, external_air_temperature);

        (pipework_heat_loss_internal, pipework_heat_loss_external)
    }
}

fn mixer_shower_input_to_shower(
    input: &MixerShowerInput,
    cold_water_sources: &ColdWaterSources,
    wwhrs: &IndexMap<String, Wwhrs>,
) -> Shower {
    let cold_water_source = match input.cold_water_source {
        ColdWaterSourceType::MainsWater => cold_water_sources.ref_for_mains_water().unwrap(),
        ColdWaterSourceType::HeaderTank => cold_water_sources.ref_for_header_tank().unwrap(),
    };
    let wwhrs_instance: Option<Wwhrs> = input
        .waste_water_heat_recovery
        .as_ref()
        .and_then(|w| wwhrs.get(&w.to_string()).cloned());

    Shower::MixerShower(MixerShower::new(
        input.flowrate,
        cold_water_source,
        wwhrs_instance,
    ))
}

fn instant_electric_shower_input_to_shower(
    name: String,
    input: &InstantElectricShowerInput,
    cold_water_sources: &ColdWaterSources,
    energy_supplies: &EnergySupplies,
) -> Shower {
    let cold_water_source = match input.cold_water_source {
        ColdWaterSourceType::MainsWater => cold_water_sources.ref_for_mains_water().unwrap(),
        ColdWaterSourceType::HeaderTank => cold_water_sources.ref_for_header_tank().unwrap(),
    };

    let energy_supply = match &input.energy_supply {
        EnergySupplyType::Electricity => energy_supplies.mains_electricity.clone().unwrap(),
        EnergySupplyType::MainsGas => energy_supplies.mains_gas.clone().unwrap(),
        EnergySupplyType::UnmetDemand => energy_supplies.unmet_demand.clone(),
        EnergySupplyType::LpgBulk => energy_supplies.bulk_lpg.clone().unwrap(),
        EnergySupplyType::HeatNetwork => energy_supplies.heat_network.clone().unwrap(),
        _ => panic!("Unexpected energy supply type for a shower"),
    };
    let energy_supply_conn = EnergySupply::connection(energy_supply, name.as_str()).unwrap();

    Shower::InstantElectricShower(InstantElectricShower::new(
        input.rated_power,
        cold_water_source,
        energy_supply_conn,
    ))
}

fn input_to_bath(input: &BathDetails, cold_water_sources: &ColdWaterSources) -> Bath {
    let cold_water_source = match input.cold_water_source {
        ColdWaterSourceType::MainsWater => cold_water_sources.ref_for_mains_water().unwrap(),
        ColdWaterSourceType::HeaderTank => cold_water_sources.ref_for_header_tank().unwrap(),
    };

    Bath::new(input.size, cold_water_source, input.flowrate)
}

fn input_to_other_water_events(
    input: &OtherWaterUseDetails,
    cold_water_sources: &ColdWaterSources,
) -> OtherHotWater {
    let cold_water_source = match input.cold_water_source {
        ColdWaterSourceType::MainsWater => cold_water_sources.ref_for_mains_water().unwrap(),
        ColdWaterSourceType::HeaderTank => cold_water_sources.ref_for_header_tank().unwrap(),
    };

    OtherHotWater::new(input.flowrate, cold_water_source)
}

fn input_to_water_distribution_pipework(
    input: &WaterPipework,
    total_number_tapping_points: usize,
) -> Pipework {
    // Calculate average length of pipework between HW system and tapping point
    let length_average = input.length / total_number_tapping_points as f64;

    Pipework::new(
        input.internal_diameter_mm / MILLIMETRES_IN_METRE as f64,
        input.external_diameter_mm / MILLIMETRES_IN_METRE as f64,
        length_average,
        input.insulation_thermal_conductivity,
        input.insulation_thickness_mm / MILLIMETRES_IN_METRE as f64,
        input.surface_reflectivity,
        input.pipe_contents,
    )
}
