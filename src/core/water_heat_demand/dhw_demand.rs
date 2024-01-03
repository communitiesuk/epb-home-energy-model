use crate::core::energy_supply::energy_supply::{EnergySupplies, EnergySupplyLoose};
use crate::core::heating_systems::wwhrs::Wwhrs;
use crate::core::pipework::{Pipework, PipeworkContentsType};
use crate::core::schedule::ScheduleEvent;
use crate::core::units::MILLIMETRES_IN_METRE;
use crate::core::water_heat_demand::bath::Bath;
use crate::core::water_heat_demand::misc::water_demand_to_kWh;
use crate::core::water_heat_demand::other_hot_water_uses::OtherHotWater;
use crate::core::water_heat_demand::shower::Shower;
use crate::core::water_heat_demand::shower::{InstantElectricShower, MixerShower};
use crate::corpus::{ColdWaterSources, HotWaterEventSchedules};
use crate::input::{
    Bath as BathInput, BathDetails, ColdWaterSourceType, EnergySupplyType,
    InstantElectricShower as InstantElectricShowerInput, MixerShower as MixerShowerInput,
    OtherWaterUse as OtherWaterUseInput, OtherWaterUseDetails, Shower as ShowerInput,
    WaterDistribution as WaterDistributionInput, WaterPipeContentsType, WaterPipework,
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
        wwhrs: &HashMap<String, Wwhrs>,
        energy_supplies: &EnergySupplies,
        event_schedules: HotWaterEventSchedules,
    ) -> Self {
        let showers: HashMap<String, Shower> = shower_input
            .iter()
            .map(|input| {
                let mut showers = vec![(
                    "mixer".to_owned(),
                    mixer_shower_input_to_shower(&input.mixer, cold_water_sources, &wwhrs),
                )];
                if let Some(ies) = &input.ies {
                    showers.push((
                        "ies".to_owned(),
                        instant_electric_shower_input_to_shower(
                            "ies".to_owned(),
                            ies,
                            cold_water_sources,
                            energy_supplies,
                        ),
                    ))
                }

                showers
            })
            .flatten()
            .collect();
        let baths = bath_input
            .iter()
            .map(|input| match &input.medium {
                Some(details) => vec![(
                    "medium".to_owned(),
                    input_to_bath(details, cold_water_sources),
                )],
                None => vec![],
            })
            .flatten()
            .collect();
        let other = other_hot_water_input
            .iter()
            .map(|input| match &input.other {
                Some(details) => vec![(
                    "other".to_owned(),
                    input_to_other_water_events(details, cold_water_sources),
                )],
                None => vec![],
            })
            .flatten()
            .collect();
        let total_number_tapping_points = showers
            .iter()
            .filter(|(_, shower)| match shower {
                Shower::InstantElectricShower(_) => false,
                _ => true,
            })
            .count();
        let hot_water_distribution_pipework = water_distribution_input
            .iter()
            .map(|input| {
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
            .flatten()
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
                        hw_energy_demand += water_demand_to_kWh(
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
                hw_energy_demand += water_demand_to_kWh(
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
                hw_energy_demand += water_demand_to_kWh(
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
}

fn mixer_shower_input_to_shower<'a>(
    input: &MixerShowerInput,
    cold_water_sources: &ColdWaterSources,
    wwhrs: &HashMap<String, Wwhrs>,
) -> Shower {
    let cold_water_source = match input.cold_water_source {
        ColdWaterSourceType::MainsWater => cold_water_sources.ref_for_mains_water().unwrap(),
        ColdWaterSourceType::HeaderTank => cold_water_sources.ref_for_header_tank().unwrap(),
    };
    let wwhrs_instance: Option<Wwhrs> = input.waste_water_heat_recovery.as_ref().and_then(|w| {
        wwhrs
            .get(&w.to_string())
            .and_then(|system| Some(system.clone()))
    });

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
        EnergySupplyType::Electricity => {
            EnergySupplyLoose::EnergySupply(energy_supplies.mains_electricity.clone().unwrap())
        }
        EnergySupplyType::MainsGas => {
            EnergySupplyLoose::EnergySupply(energy_supplies.mains_gas.clone().unwrap())
        }
        EnergySupplyType::UnmetDemand => {
            EnergySupplyLoose::EnergySupply(energy_supplies.unmet_demand.clone())
        }
        EnergySupplyType::LpgBulk => {
            EnergySupplyLoose::EnergySupply(energy_supplies.bulk_lpg.clone().unwrap())
        }
        EnergySupplyType::HeatNetwork => {
            EnergySupplyLoose::HeatNetwork(energy_supplies.heat_network.clone().unwrap())
        }
        _ => panic!("Unexpected energy supply type for a shower"),
    };

    Shower::InstantElectricShower(InstantElectricShower::new(
        input.rated_power,
        cold_water_source,
        energy_supply,
        name,
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
        input.insulation_thickness_mm,
        input.surface_reflectivity,
        match input.pipe_contents {
            WaterPipeContentsType::Water => PipeworkContentsType::Water,
        },
    )
}
