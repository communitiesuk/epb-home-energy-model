use crate::core::common::WaterSourceWithTemperature;
use crate::core::controls::time_control::{
    Control, HeatSourceControl, OnOffMinimisingTimeControl, OnOffTimeControl, SetpointTimeControl,
    ToUChargeControl,
};
use crate::core::cooling_systems::air_conditioning::AirConditioning;
use crate::core::ductwork::Ductwork;
use crate::core::energy_supply::elec_battery::ElectricBattery;
use crate::core::energy_supply::energy_supply::{
    EnergySupplies, EnergySupply, EnergySupplyConnection, UNMET_DEMAND_SUPPLY_NAME,
};
use crate::core::energy_supply::pv::PhotovoltaicSystem;
use crate::core::heating_systems::boiler::{Boiler, BoilerServiceWaterCombi};
use crate::core::heating_systems::common::{HeatSourceWet, SpaceHeatSystem, SpaceHeatingService};
use crate::core::heating_systems::emitters::{Emitters, EmittersDetailedResult};
use crate::core::heating_systems::heat_battery::HeatBattery;
use crate::core::heating_systems::heat_network::{HeatNetwork, HeatNetworkServiceWaterDirect};
use crate::core::heating_systems::heat_pump::{
    HeatPump, HeatPumpHotWaterOnly, ResultsAnnual, ResultsPerTimestep,
};
use crate::core::heating_systems::instant_elec_heater::InstantElecHeater;
use crate::core::heating_systems::point_of_use::PointOfUse;
use crate::core::heating_systems::storage_tank::{
    HeatSourceWithStorageTank, ImmersionHeater, PVDiverter, SolarThermalSystem, StorageTank,
};
use crate::core::heating_systems::wwhrs::{
    WWHRSInstantaneousSystemA, WWHRSInstantaneousSystemB, WWHRSInstantaneousSystemC, Wwhrs,
};
use crate::core::material_properties::WATER;
use crate::core::schedule::{
    expand_boolean_schedule, expand_events, expand_numeric_schedule, reject_nulls, ScheduleEvent,
    TypedScheduleEvent, WaterScheduleEventType,
};
use crate::core::space_heat_demand::building_element::{
    convert_uvalue_to_resistance, pitch_class, BuildingElement, BuildingElementAdjacentZTC,
    BuildingElementAdjacentZTUSimple, BuildingElementGround, BuildingElementOpaque,
    BuildingElementTransparent, HeatFlowDirection, PITCH_LIMIT_HORIZ_CEILING,
};
use crate::core::space_heat_demand::internal_gains::{
    ApplianceGains, EventApplianceGains, Gains, InternalGains,
};
use crate::core::space_heat_demand::thermal_bridge::{ThermalBridge, ThermalBridging};
use crate::core::space_heat_demand::ventilation::{
    AirTerminalDevices, CombustionAppliances, InfiltrationVentilation, MechanicalVentilation, Vent,
    VentilationDetailedResult, Window,
};
use crate::core::space_heat_demand::zone::{
    AirChangesPerHourArgument, HeatBalance, HeatBalanceFieldName, Zone,
};
use crate::core::units::{
    kelvin_to_celsius, MILLIMETRES_IN_METRE, SECONDS_PER_HOUR, WATTS_PER_KILOWATT,
};
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::dhw_demand::{
    DemandVolTargetKey, DomesticHotWaterDemand, DomesticHotWaterDemandData, VolumeReference,
};
use crate::core::water_heat_demand::misc::water_demand_to_kwh;
use crate::errors::NotImplementedError;
use crate::external_conditions::ExternalConditions;
use crate::input::{
    init_orientation, ApplianceGains as ApplianceGainsInput, ApplianceGainsDetails,
    BuildingElement as BuildingElementInput, ChargeLevel, ColdWaterSourceDetails,
    ColdWaterSourceInput, ColdWaterSourceReference, ColdWaterSourceType, Control as ControlInput,
    ControlDetails, DuctShape, DuctType, EnergyDiverter, EnergySupplyDetails, EnergySupplyInput,
    EnergySupplyKey, EnergySupplyType, ExternalConditionsInput, FloorType, FuelType,
    HeatPumpSourceType, HeatSource as HeatSourceInput, HeatSourceControlType, HeatSourceWetDetails,
    HeatSourceWetType, HotWaterSourceDetails,
    InfiltrationVentilation as InfiltrationVentilationInput, Input,
    InternalGains as InternalGainsInput, InternalGainsDetails, OnSiteGeneration,
    OnSiteGenerationDetails, SpaceCoolSystem as SpaceCoolSystemInput, SpaceCoolSystemDetails,
    SpaceCoolSystemType, SpaceHeatSystem as SpaceHeatSystemInput, SpaceHeatSystemDetails,
    SystemReference, ThermalBridging as ThermalBridgingInput, ThermalBridgingDetails, VentType,
    VentilationLeaks, WasteWaterHeatRecovery, WasteWaterHeatRecoveryDetails, WaterHeatingEvent,
    WaterHeatingEventType, WaterHeatingEvents, WwhrsType, ZoneDictionary, ZoneInput,
};
use crate::simulation_time::{SimulationTimeIteration, SimulationTimeIterator};
use crate::ProjectFlags;
use anyhow::{anyhow, bail};
use arrayvec::ArrayString;
use indexmap::IndexMap;
#[cfg(feature = "indicatif")]
use indicatif::ProgressIterator;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use parking_lot::{Mutex, RwLock};
use serde_enum_str::{Deserialize_enum_str, Serialize_enum_str};
use std::borrow::Borrow;
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::fmt::{Display, Formatter};
use std::hash::Hash;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

/// As of adopting Rust 1.82 as an MSRV we'll be able to declare this using constants as it supports floating-point arithmetic at compile time
fn temp_setpnt_heat_none() -> f64 {
    kelvin_to_celsius(0.0).expect("Not below absolute zero")
}

fn temp_setpnt_cool_none() -> f64 {
    kelvin_to_celsius(1.4e32).expect("Not below absolute zero")
}

// used for calculations re internal gains from pipework
const FRAC_DHW_ENERGY_INTERNAL_GAINS: f64 = 0.25;

#[derive(Debug)]
pub struct Corpus {
    pub simulation_time: Arc<SimulationTimeIterator>,
    pub external_conditions: Arc<ExternalConditions>,
    pub cold_water_sources: ColdWaterSources,
    pub energy_supplies: IndexMap<String, Arc<RwLock<EnergySupply>>>,
    pub internal_gains: InternalGainsCollection,
    pub controls: Controls,
    pub wwhrs: IndexMap<String, Arc<Mutex<Wwhrs>>>,
    pub domestic_hot_water_demand: DomesticHotWaterDemand,
    pub ventilation: Arc<InfiltrationVentilation>,
    mechanical_ventilations: IndexMap<String, Arc<MechanicalVentilation>>,
    pub space_heating_ductwork: IndexMap<String, Vec<Ductwork>>,
    pub zones: Arc<IndexMap<String, Zone>>,
    pub energy_supply_conn_unmet_demand_zone: IndexMap<String, Arc<EnergySupplyConnection>>,
    pub heat_system_name_for_zone: IndexMap<String, String>,
    pub cool_system_name_for_zone: IndexMap<String, String>,
    pub total_floor_area: f64,
    pub total_volume: f64,
    pub wet_heat_sources: IndexMap<String, WetHeatSource>,
    pub hot_water_sources: IndexMap<String, HotWaterSource>,
    pub heat_system_names_requiring_overvent: Vec<String>,
    pub heat_sources_wet_with_buffer_tank: Vec<String>,
    pub space_heat_systems: IndexMap<String, Arc<Mutex<SpaceHeatSystem>>>,
    pub space_cool_systems: IndexMap<String, AirConditioning>,
    pub on_site_generation: IndexMap<String, PhotovoltaicSystem>,
    pub diverters: Vec<Arc<RwLock<PVDiverter>>>,
    required_vent_data: Option<RequiredVentData>,
    energy_supply_conn_names_for_hot_water_source: IndexMap<String, Vec<String>>,
    energy_supply_conn_names_for_heat_systems: IndexMap<String, String>,
    timestep_end_calcs: Arc<RwLock<Vec<WetHeatSource>>>,
    initial_loop: AtomicBool,
    detailed_output_heating_cooling: bool,
}

impl Corpus {
    pub fn from_inputs(
        input: &Input,
        external_conditions: Option<&ExternalConditions>,
        output_options: &OutputOptions,
    ) -> anyhow::Result<Self> {
        let simulation_time_iterator = Arc::new(input.simulation_time.iter());

        let external_conditions = Arc::new(match external_conditions {
            Some(external_conditions) => external_conditions.clone(),
            None => external_conditions_from_input(
                input.external_conditions.clone(),
                &simulation_time_iterator,
            ),
        });

        let diverter_types: DiverterTypes = input
            .energy_supply
            .iter()
            .flat_map(|(name, supply)| {
                diverter_from_energy_supply(supply).map(|diverter| (name.to_owned(), diverter))
            })
            .collect();
        let mut diverters: Vec<Arc<RwLock<PVDiverter>>> = Default::default();

        let cold_water_sources = cold_water_sources_from_input(&input.cold_water_source);
        let wwhrs = wwhrs_from_input(
            input.waste_water_heat_recovery.as_ref(),
            &cold_water_sources,
            simulation_time_iterator.current_iteration(),
        );

        let mut energy_supplies = energy_supplies_from_input(
            &input.energy_supply,
            simulation_time_iterator.clone().as_ref(),
            external_conditions.clone(),
        );

        let controls =
            control_from_input(&input.control, simulation_time_iterator.clone().as_ref())?;

        let event_schedules = event_schedules_from_input(
            &input.water_heating_events,
            simulation_time_iterator.as_ref(),
        )?;

        let domestic_hot_water_demand = DomesticHotWaterDemand::new(
            input.hot_water_demand.shower.clone().unwrap_or_default(),
            input.hot_water_demand.bath.clone().unwrap_or_default(),
            input
                .hot_water_demand
                .other_water_use
                .clone()
                .unwrap_or_default(),
            match &input.hot_water_source.hot_water_cylinder {
                HotWaterSourceDetails::PointOfUse { .. } => None,
                _ => input.hot_water_demand.water_distribution.clone(),
            },
            &cold_water_sources,
            &wwhrs,
            &energy_supplies,
            event_schedules,
        )?;

        let _opening_area_total_from_zones = opening_area_total_from_zones(&input.zone);

        let total_volume = input.zone.values().map(|zone| zone.volume).sum::<f64>();

        let mut heat_system_name_for_zone: IndexMap<String, String> = Default::default();
        let mut cool_system_name_for_zone: IndexMap<String, String> = Default::default();

        // infiltration ventilation
        let (
            infiltration_ventilation,
            window_adjust_control,
            mechanical_ventilations,
            space_heating_ductwork,
        ) = infiltration_ventilation_from_input(
            &input.zone,
            &input.infiltration_ventilation,
            &controls,
            &mut energy_supplies,
            simulation_time_iterator.as_ref(),
            total_volume,
            external_conditions.clone(),
            output_options.detailed_output_heating_cooling,
        )?;

        let infiltration_ventilation = Arc::from(infiltration_ventilation);

        let required_vent_data = required_vent_data_from_input(&input.control)?;

        let zones: IndexMap<String, Zone> = input
            .zone
            .iter()
            .map(|(i, zone)| -> anyhow::Result<(String, Zone)> {
                Ok(((*i).clone(), {
                    let (zone_for_corpus, heat_system_name, cool_system_name) = zone_from_input(
                        zone,
                        external_conditions.clone(),
                        infiltration_ventilation.clone(),
                        window_adjust_control.clone(),
                        output_options.print_heat_balance,
                        simulation_time_iterator.clone().as_ref(),
                    )?;
                    // TODO 0.32 correct logic here for possible multiple heat system names
                    heat_system_name_for_zone.insert(
                        (*i).clone(),
                        match heat_system_name {
                            SystemReference::None(_) => String::default(),
                            SystemReference::Single(name) => name.clone(),
                            SystemReference::Multiple(names) => {
                                names.first().cloned().unwrap_or_default()
                            }
                        },
                    );
                    cool_system_name_for_zone.insert(
                        (*i).clone(),
                        match cool_system_name {
                            SystemReference::None(_) => String::default(),
                            SystemReference::Single(name) => name.clone(),
                            SystemReference::Multiple(names) => {
                                names.first().cloned().unwrap_or_default()
                            }
                        },
                    );

                    zone_for_corpus
                }))
            })
            .collect::<anyhow::Result<_>>()?;
        let zones = Arc::new(zones);

        if !has_unique_non_default_values(&heat_system_name_for_zone)
            || !has_unique_non_default_values(&cool_system_name_for_zone)
        {
            bail!("the heat or cool systems do not have unique names in the inputs");
        }

        let energy_supply_conn_unmet_demand_zone = set_up_energy_supply_unmet_demand_zones(
            energy_supplies[UNMET_DEMAND_SUPPLY_NAME].clone(),
            &input.zone,
        );
        // TODO: there needs to be some equivalent here of the Python code that builds the dict __energy_supply_conn_unmet_demand_zone

        let total_floor_area = zones.values().fold(0., |acc, zone| zone.area() + acc);

        let mut internal_gains =
            internal_gains_from_input(&input.internal_gains, total_floor_area)?;

        apply_appliance_gains_from_input(
            &mut internal_gains,
            &input.appliance_gains,
            &mut energy_supplies,
            total_floor_area,
            simulation_time_iterator.as_ref(),
        )?;

        let timestep_end_calcs: Arc<RwLock<Vec<WetHeatSource>>> = Default::default();
        let mut heat_sources_wet_with_buffer_tank = vec![];

        let mut wet_heat_sources: IndexMap<String, WetHeatSource> = input
            .heat_source_wet
            .clone()
            .unwrap_or_default()
            .iter()
            .map(|(name, heat_source_wet_details)| {
                let heat_source = heat_source_wet_from_input(
                    name,
                    (*heat_source_wet_details).clone(),
                    external_conditions.clone(),
                    simulation_time_iterator.clone(),
                    &mechanical_ventilations,
                    zones.len(),
                    TempInternalAirAccessor {
                        zones: zones.clone(),
                        total_volume,
                    },
                    &controls,
                    &mut energy_supplies,
                    output_options.detailed_output_heating_cooling,
                )?;
                timestep_end_calcs.write().push(heat_source.clone());
                if let HeatSourceWetDetails::HeatPump {
                    buffer_tank: Some(_),
                    ..
                } = heat_source_wet_details
                {
                    heat_sources_wet_with_buffer_tank.push(name.clone());
                }
                anyhow::Ok(((*name).clone(), heat_source))
            })
            .collect::<anyhow::Result<IndexMap<String, WetHeatSource>>>()?;

        let mut hot_water_sources: IndexMap<String, HotWaterSource> = Default::default();
        let mut energy_supply_conn_names_for_hot_water_source: IndexMap<String, Vec<String>> =
            Default::default();
        let (hot_water_source, hw_cylinder_conn_names) = hot_water_source_from_input(
            "hw cylinder".to_string(),
            &input.hot_water_source.hot_water_cylinder,
            &cold_water_sources,
            &mut wet_heat_sources,
            &wwhrs,
            &controls,
            &mut energy_supplies,
            &diverter_types,
            &mut diverters,
            TempInternalAirAccessor {
                zones: zones.clone(),
                total_volume,
            },
            simulation_time_iterator.clone().as_ref(),
            external_conditions.clone(),
        )?;
        hot_water_sources.insert("hw cylinder".to_string(), hot_water_source);
        energy_supply_conn_names_for_hot_water_source
            .insert("hw cylinder".to_string(), hw_cylinder_conn_names);

        let mut heat_system_names_requiring_overvent: Vec<String> = Default::default();

        let (space_heat_systems, energy_supply_conn_names_for_heat_systems) = input
            .space_heat_system
            .as_ref()
            .map(|system| {
                anyhow::Ok(space_heat_systems_from_input(
                    system,
                    &controls,
                    &mut energy_supplies,
                    simulation_time_iterator.as_ref(),
                    &wet_heat_sources,
                    &mut heat_system_names_requiring_overvent,
                    &heat_system_name_for_zone,
                    &zones,
                    &heat_sources_wet_with_buffer_tank,
                    TempInternalAirAccessor {
                        zones: zones.clone(),
                        total_volume,
                    },
                    external_conditions.clone(),
                    output_options.detailed_output_heating_cooling,
                )?)
            })
            .transpose()?
            .unwrap_or_default();

        let space_cool_systems = input
            .space_cool_system
            .as_ref()
            .map(|system| {
                anyhow::Ok(space_cool_systems_from_input(
                    system,
                    cool_system_name_for_zone
                        .values()
                        .map(|s| s.as_str())
                        .collect::<Vec<_>>(),
                    &controls,
                    &mut energy_supplies,
                    &simulation_time_iterator,
                )?)
            })
            .transpose()?
            .unwrap_or_default();

        let on_site_generation = input
            .on_site_generation
            .as_ref()
            .map(|on_site_generation| {
                anyhow::Ok(on_site_generation_from_input(
                    on_site_generation,
                    &mut energy_supplies,
                    external_conditions.clone(),
                    &simulation_time_iterator,
                )?)
            })
            .transpose()?
            .unwrap_or_default();

        Ok(Self {
            simulation_time: simulation_time_iterator,
            external_conditions,
            cold_water_sources,
            energy_supplies,
            internal_gains,
            controls,
            wwhrs,
            domestic_hot_water_demand,
            ventilation: infiltration_ventilation,
            mechanical_ventilations,
            space_heating_ductwork,
            zones,
            energy_supply_conn_unmet_demand_zone,
            heat_system_name_for_zone,
            cool_system_name_for_zone,
            total_floor_area,
            total_volume,
            wet_heat_sources,
            hot_water_sources,
            heat_system_names_requiring_overvent,
            heat_sources_wet_with_buffer_tank,
            space_heat_systems,
            space_cool_systems,
            on_site_generation,
            diverters,
            required_vent_data,
            energy_supply_conn_names_for_hot_water_source,
            energy_supply_conn_names_for_heat_systems,
            timestep_end_calcs,
            initial_loop: AtomicBool::new(false),
            detailed_output_heating_cooling: output_options.detailed_output_heating_cooling,
        })
    }

    pub fn total_floor_area(&self) -> f64 {
        self.total_floor_area
    }

    /// Calculate heat transfer coefficient (HTC) and heat loss parameter (HLP)
    /// according to the SAP10.2 specification
    pub fn calc_htc_hlp(&self) -> (f64, f64, HashMap<String, f64>, HashMap<String, f64>) {
        let mut htc_map: HashMap<String, f64> = Default::default();
        let mut hlp_map: HashMap<String, f64> = Default::default();

        // Calculate the total fabric heat loss, total heat capacity, total ventilation heat
        // loss and total heat transfer coeffient for thermal bridges across all zones
        for (z_name, zone) in self.zones.iter() {
            let fabric_heat_loss = zone.total_fabric_heat_loss();
            let thermal_bridges = zone.total_thermal_bridges();
            let vent_heat_loss = zone.total_vent_heat_loss();

            // Calculate the heat transfer coefficent (HTC), in W / K
            // TODO (from Python) check ventilation losses are correct
            let htc = fabric_heat_loss + thermal_bridges + vent_heat_loss;

            // Calculate the HLP, in W/m2 K
            let hlp = htc / zone.area();

            htc_map.insert((*z_name).clone(), htc);
            hlp_map.insert((*z_name).clone(), hlp);
        }

        let total_htc = htc_map.values().sum();
        let total_hlp = total_htc / self.total_floor_area;

        (total_htc, total_hlp, htc_map, hlp_map)
    }

    /// Calculate the total heat capacity normalised for floor area
    pub fn calc_hcp(&self) -> f64 {
        // TODO (from Python) party walls and solid doors should be exluded according to SAP spec - if party walls are
        // assumed to be ZTU building elements this could be set to zero?
        let total_heat_capacity = self
            .zones
            .values()
            .map(|zone| zone.total_heat_capacity())
            .sum::<f64>();

        total_heat_capacity / self.total_floor_area
    }

    /// Calculate the heat loss form factor, defined as exposed area / floor area
    pub fn calc_hlff(&self) -> f64 {
        let total_heat_loss_area = self
            .zones
            .values()
            .map(|zone| zone.total_heat_loss_area())
            .sum::<f64>();

        total_heat_loss_area / self.total_floor_area
    }

    pub fn temp_internal_air(&self) -> f64 {
        temp_internal_air_for_zones(self.zones.clone(), self.total_volume)
    }

    fn pipework_losses_and_internal_gains_from_hw_storage_tank(
        &self,
        delta_t_h: f64,
        volume_water_remove_from_tank: f64,
        hw_duration: f64,
        no_of_hw_events: usize,
        temp_final_drawoff: f64,
        temp_average_drawoff: f64,
        temp_hot_water: f64,
        vol_hot_water_equiv_elec_shower: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (f64, f64, f64) {
        let (pw_losses_internal, pw_losses_external) = self.calc_pipework_losses(
            delta_t_h,
            hw_duration,
            no_of_hw_events,
            temp_final_drawoff,
            simulation_time_iteration,
        );
        let gains_internal_dhw_use_storagetank = FRAC_DHW_ENERGY_INTERNAL_GAINS
            * water_demand_to_kwh(
                volume_water_remove_from_tank,
                temp_average_drawoff,
                self.temp_internal_air(),
            );

        let gains_internal_dhw_use_ies = FRAC_DHW_ENERGY_INTERNAL_GAINS
            * water_demand_to_kwh(
                vol_hot_water_equiv_elec_shower,
                temp_hot_water,
                self.temp_internal_air(),
            );

        let gains_internal_dhw_use =
            gains_internal_dhw_use_storagetank + gains_internal_dhw_use_ies;

        // Return:
        // losses from internal distribution pipework (kWh)
        // losses from external distribution pipework (kWh)
        // internal gains due to hot water use (kWh)
        (
            pw_losses_internal,
            pw_losses_external,
            gains_internal_dhw_use,
        )
    }

    /// Return:
    ///  # - losses from internal distribution pipework (kWh)
    ///  # - losses from external distribution pipework (kWh)
    ///  # - internal gains due to hot water use (kWh)
    fn pipework_losses_and_internal_gains_from_hw(
        &self,
        delta_t_h: f64,
        vol_hot_water_at_tapping_point: f64,
        hw_duration: f64,
        no_of_hw_events: usize,
        temp_hot_water: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (f64, f64, f64) {
        let (pw_losses_internal, pw_losses_external) = self.calc_pipework_losses(
            delta_t_h,
            hw_duration,
            no_of_hw_events,
            temp_hot_water,
            simulation_time_iteration,
        );

        let gains_internal_dhw_use = FRAC_DHW_ENERGY_INTERNAL_GAINS
            * water_demand_to_kwh(
                vol_hot_water_at_tapping_point,
                temp_hot_water,
                self.temp_internal_air(),
            );

        (
            pw_losses_internal,
            pw_losses_external,
            gains_internal_dhw_use,
        )
    }

    fn calc_pipework_losses(
        &self,
        delta_t_h: f64,
        hw_duration: f64,
        no_of_hw_events: usize,
        temp_hot_water: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (f64, f64) {
        let demand_water_temperature = temp_hot_water;
        let internal_air_temperature = self.temp_internal_air();
        let external_air_temperature = self
            .external_conditions
            .air_temp(&simulation_time_iteration);

        self.domestic_hot_water_demand.calc_pipework_losses(
            delta_t_h,
            hw_duration,
            no_of_hw_events,
            demand_water_temperature,
            internal_air_temperature,
            external_air_temperature,
        )
    }

    /// Calculate the losses in the buffer tank
    fn calc_internal_gains_buffer_tank(&self) -> f64 {
        self.heat_sources_wet_with_buffer_tank
            .iter()
            .map(
                |heat_source_name| match self.wet_heat_sources.get(heat_source_name) {
                    Some(heat_source) => match heat_source {
                        WetHeatSource::HeatPump(heat_pump) => heat_pump.lock().buffer_int_gains(),
                        _ => unreachable!(),
                    },
                    None => 0.,
                },
            )
            .sum::<f64>()
    }

    /// Calculate the losses/gains in the MVHR ductwork
    fn calc_internal_gains_ductwork(&self, simulation_time: SimulationTimeIteration) -> f64 {
        let mut internal_gains_ductwork_watts = 0.0;
        for mvhr_ductwork in self.space_heating_ductwork.values() {
            // assume MVHR unit is running 100% of the time
            for duct in mvhr_ductwork {
                match duct.duct_type() {
                    DuctType::Intake | DuctType::Exhaust => {
                        // Heat loss from intake or exhaust ducts is to zone, so add
                        // to internal gains (may be negative gains)
                        internal_gains_ductwork_watts += duct.total_duct_heat_loss(
                            self.temp_internal_air(),
                            self.external_conditions.air_temp(&simulation_time),
                        );
                    }
                    DuctType::Supply | DuctType::Extract => {
                        // Heat loss from supply and extract ducts is to outside, so
                        // subtract from internal gains
                        internal_gains_ductwork_watts -= duct.total_duct_heat_loss(
                            self.temp_internal_air(),
                            self.external_conditions.air_temp(&simulation_time),
                        );
                    }
                }
            }
        }
        internal_gains_ductwork_watts
    }

    fn space_heat_internal_gains_for_zone(
        &self,
        zone: &Zone,
        gains_internal_dhw: f64,
        internal_gains_ductwork_per_m3: f64,
        gains_internal_buffer_tank: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        // Initialise to dhw internal gains split proportionally to zone floor area
        let mut gains_internal_zone =
            (gains_internal_buffer_tank + gains_internal_dhw) * zone.area() / self.total_floor_area;

        for internal_gains in self.internal_gains.values() {
            gains_internal_zone += internal_gains.total_internal_gain_in_w(zone.area(), simtime)?;
        }

        // Add gains from ventilation fans (also calculates elec demand from fans)
        // TODO (from Python) Remove the branch on the type of ventilation (find a better way)
        for mech_vent in self.mechanical_ventilations.values() {
            gains_internal_zone += mech_vent.fans(zone.volume(), self.total_volume, None, &simtime);
            gains_internal_zone += internal_gains_ductwork_per_m3 * zone.volume();
        }

        Ok(gains_internal_zone)
    }

    /// Look up relevant heating and cooling systems for the specified zone
    fn heat_cool_systems_for_zone(
        &self,
        z_name: &str,
        simtime: SimulationTimeIteration,
    ) -> (Vec<String>, Vec<String>, SetpointsAndConvectiveFractions) {
        // TODO (from Python) For now, the existing single system inputs are each added to a
        //      list. This will eventually need to handle a list being specified
        //      in the inputs.
        let h_name_list = self
            .heat_system_name_for_zone
            .get(z_name)
            .iter()
            .map(|x| (*x).clone())
            .collect::<Vec<String>>();
        let c_name_list = self
            .cool_system_name_for_zone
            .get(z_name)
            .iter()
            .map(|x| (*x).clone())
            .collect::<Vec<String>>();

        let SetpointsAndConvectiveFractions {
            temp_setpnt_heat: temp_setpnt_heat_system,
            temp_setpnt_cool: temp_setpnt_cool_system,
            frac_convective_heat: frac_convective_heat_system,
            frac_convective_cool: frac_convective_cool_system,
        } = self.setpoints_and_convective_fractions(
            h_name_list.as_slice(),
            c_name_list.as_slice(),
            simtime,
        );

        // Sort heating and cooling systems by setpoint (highest first for
        // heating, lowest first for cooling)
        // TODO (from Python) In the event of two systems having the same setpoint, make
        //      sure the one listed first by the user takes priority
        let h_name_list_sorted: Vec<String> = temp_setpnt_heat_system
            .iter()
            .sorted_by(|a, b| OrderedFloat(*a.1).cmp(&OrderedFloat(*b.1)))
            .rev()
            .map(|x| x.0.to_owned())
            .collect();
        let c_name_list_sorted: Vec<String> = temp_setpnt_cool_system
            .iter()
            .sorted_by(|a, b| OrderedFloat(*a.1).cmp(&OrderedFloat(*b.1)))
            .map(|x| x.0.to_owned())
            .collect();

        (
            h_name_list_sorted,
            c_name_list_sorted,
            SetpointsAndConvectiveFractions {
                temp_setpnt_heat: temp_setpnt_heat_system,
                temp_setpnt_cool: temp_setpnt_cool_system,
                frac_convective_heat: frac_convective_heat_system,
                frac_convective_cool: frac_convective_cool_system,
            },
        )
    }

    fn setpoints_and_convective_fractions(
        &self,
        h_name_list: &[String],
        c_name_list: &[String],
        simtime: SimulationTimeIteration,
    ) -> SetpointsAndConvectiveFractions {
        let mut frac_convective_heat: IndexMap<String, f64> = Default::default();
        let mut frac_convective_cool: IndexMap<String, f64> = Default::default();
        let mut temp_setpnt_heat: IndexMap<String, f64> = Default::default();
        let mut temp_setpnt_cool: IndexMap<String, f64> = Default::default();

        for h_name in h_name_list {
            match h_name.as_str() {
                h_name @ "" => {
                    frac_convective_heat.insert((*h_name).to_owned(), 1.0);
                    temp_setpnt_heat.insert((*h_name).to_owned(), temp_setpnt_heat_none());
                }
                h_name => {
                    let space_heat_system = self.space_heat_systems.get(h_name).unwrap().lock();
                    frac_convective_heat
                        .insert((*h_name).to_owned(), space_heat_system.frac_convective());
                    temp_setpnt_heat.insert(
                        (*h_name).to_owned(),
                        space_heat_system
                            .temp_setpnt(simtime)
                            .unwrap_or_else(temp_setpnt_heat_none),
                    );
                }
            }
        }

        for c_name in c_name_list {
            match c_name.as_str() {
                c_name @ "" => {
                    frac_convective_cool.insert((*c_name).to_owned(), 1.0);
                    temp_setpnt_cool.insert((*c_name).to_owned(), temp_setpnt_cool_none());
                }
                c_name => {
                    let space_cool_system = self.space_cool_systems.get(c_name).unwrap();
                    frac_convective_cool
                        .insert((*c_name).to_owned(), space_cool_system.frac_convective());
                    temp_setpnt_cool.insert(
                        (*c_name).to_owned(),
                        space_cool_system
                            .temp_setpnt(&simtime)
                            .unwrap_or_else(temp_setpnt_cool_none),
                    );
                }
            }
        }

        SetpointsAndConvectiveFractions {
            temp_setpnt_heat,
            temp_setpnt_cool,
            frac_convective_heat,
            frac_convective_cool,
        }
    }

    fn gains_heat_cool(
        &self,
        delta_t_h: f64,
        hc_output_convective: &IndexMap<&str, f64>,
        hc_output_radiative: &IndexMap<&str, f64>,
    ) -> (f64, f64) {
        let gains_heat_cool_convective =
            hc_output_convective.values().sum::<f64>() * WATTS_PER_KILOWATT as f64 / delta_t_h;
        let gains_heat_cool_radiative =
            hc_output_radiative.values().sum::<f64>() * WATTS_PER_KILOWATT as f64 / delta_t_h;

        (gains_heat_cool_convective, gains_heat_cool_radiative)
    }

    /// Calculate the incoming air changes per hour
    /// initial_p_z_ref_guess is used for calculation in first timestep.
    /// Later timesteps use the previous timesteps p_z_ref of max and min ACH ,respective to calc.
    fn calc_air_changes_per_hour(
        &self,
        temp_int_air: f64,
        r_w_arg: f64,
        initial_p_z_ref_guess: f64,
        reporting_flag: ReportingFlag,
        internal_pressure_window: &mut HashMap<ReportingFlag, f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let current_internal_pressure_window = if self.initial_loop.load(Ordering::SeqCst) {
            self.ventilation.calculate_internal_reference_pressure(
                initial_p_z_ref_guess,
                temp_int_air,
                Some(r_w_arg),
                simtime,
            )?
        } else {
            self.ventilation.calculate_internal_reference_pressure(
                internal_pressure_window[&reporting_flag],
                temp_int_air,
                Some(r_w_arg),
                simtime,
            )?
        };

        internal_pressure_window.insert(reporting_flag, current_internal_pressure_window);

        let incoming_air_flow = self.ventilation.incoming_air_flow(
            current_internal_pressure_window,
            temp_int_air,
            r_w_arg,
            Some(reporting_flag),
            Some(true),
            simtime,
        )?;

        Ok(incoming_air_flow / self.total_volume)
    }

    /// Calculate space heating demand, heating system output and temperatures
    ///
    /// Arguments:
    /// * `delta_t_h` - calculation timestep, in hours
    /// * `gains_internal_dhw` - internal gains from hot water system for this timestep, in W
    fn calc_space_heating(
        &self,
        delta_t_h: f64,
        gains_internal_dhw: f64,
        internal_pressure_window: &mut HashMap<ReportingFlag, f64>,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<SpaceHeatingCalculation> {
        let temp_ext_air = self.external_conditions.air_temp(&simtime);
        let temp_int_air = self.temp_internal_air();
        // Calculate timestep in seconds
        let delta_t = delta_t_h * SECONDS_PER_HOUR as f64;

        let internal_gains_ductwork = self.calc_internal_gains_ductwork(simtime);
        let internal_gains_ductwork_per_m3 = internal_gains_ductwork / self.total_volume;

        let internal_gains_buffer_tank = self.calc_internal_gains_buffer_tank();

        // Windows shut
        let ach_windows_shut = self.calc_air_changes_per_hour(
            temp_int_air,
            0.,
            0.,
            ReportingFlag::Min,
            internal_pressure_window,
            simtime,
        )?;

        // Windows fully open
        let ach_windows_open = self.calc_air_changes_per_hour(
            temp_int_air,
            1.,
            0.,
            ReportingFlag::Max,
            internal_pressure_window,
            simtime,
        )?;

        // To indicate the future loop should involve the p_Z_ref from previous calc
        self.initial_loop.store(false, Ordering::SeqCst);

        let ach_target = if let Some(required_vent_data) = self.required_vent_data.as_ref() {
            let ach_target = required_vent_data.schedule[simtime.time_series_idx(
                required_vent_data.start_day,
                required_vent_data.time_series_step,
            )];

            ach_windows_shut.max(ach_target.unwrap_or(ach_windows_open).min(ach_windows_open))
        } else {
            ach_windows_shut
        };

        let mut gains_internal_zone: HashMap<String, f64> = Default::default();
        let mut gains_solar_zone: HashMap<String, f64> = Default::default();
        let mut h_name_list_sorted_zone: HashMap<&str, Vec<String>> = Default::default();
        let mut c_name_list_sorted_zone: HashMap<&str, Vec<String>> = Default::default();
        let mut temp_setpnt_heat_zone_system: HashMap<&str, IndexMap<String, f64>> =
            Default::default();
        let mut temp_setpnt_cool_zone_system: HashMap<&str, IndexMap<String, f64>> =
            Default::default();
        let mut frac_convective_heat_zone_system: HashMap<&str, IndexMap<String, f64>> =
            Default::default();
        let mut frac_convective_cool_zone_system: HashMap<&str, IndexMap<String, f64>> =
            Default::default();
        let mut ach_cooling_zone: HashMap<&str, f64> = Default::default();
        let mut ach_to_trigger_heating_zone: HashMap<&str, Option<f64>> = Default::default();
        let mut internal_air_temp: HashMap<String, f64> = Default::default();
        let mut operative_temp: HashMap<String, f64> = Default::default();
        let mut space_heat_demand_zone: HashMap<String, f64> = Default::default();
        let mut space_cool_demand_zone: HashMap<String, f64> = Default::default();
        let mut space_heat_demand_system: HashMap<String, f64> = Default::default();
        let mut space_cool_demand_system: HashMap<String, f64> = Default::default();
        let mut space_heat_provided_system: HashMap<String, f64> = Default::default();
        let mut space_cool_provided_system: HashMap<String, f64> = Default::default();
        let mut heat_balance_map: HashMap<String, Option<HeatBalance>> = Default::default();

        let avg_air_supply_temp = self.ventilation.temp_supply(simtime);

        for (z_name, zone) in self.zones.iter() {
            let z_name = z_name.as_str();
            // Calculate internal and solar gains
            gains_internal_zone.insert(
                z_name.to_string(),
                self.space_heat_internal_gains_for_zone(
                    zone,
                    gains_internal_dhw,
                    internal_gains_ductwork_per_m3,
                    internal_gains_buffer_tank,
                    simtime,
                )?,
            );
            gains_solar_zone.insert(z_name.to_string(), zone.gains_solar(simtime));

            // Get heating and cooling characteristics for the current zone
            let (
                h_name_list_sorted_zone_current,
                c_name_list_sorted_zone_current,
                SetpointsAndConvectiveFractions {
                    temp_setpnt_heat: temp_setpnt_heat_zone_system_current,
                    temp_setpnt_cool: temp_setpnt_cool_zone_system_current,
                    frac_convective_heat: frac_convective_heat_zone_system_current,
                    frac_convective_cool: frac_convective_cool_zone_system_current,
                },
            ) = self.heat_cool_systems_for_zone(z_name, simtime);

            h_name_list_sorted_zone.insert(z_name, h_name_list_sorted_zone_current);
            c_name_list_sorted_zone.insert(z_name, c_name_list_sorted_zone_current);
            temp_setpnt_heat_zone_system.insert(z_name, temp_setpnt_heat_zone_system_current);
            temp_setpnt_cool_zone_system.insert(z_name, temp_setpnt_cool_zone_system_current);
            frac_convective_heat_zone_system
                .insert(z_name, frac_convective_heat_zone_system_current);
            frac_convective_cool_zone_system
                .insert(z_name, frac_convective_cool_zone_system_current);

            // Calculate space heating demand based on highest-priority systems,
            // assuming no output from any other systems

            let (
                space_heat_demand_zone_current,
                space_cool_demand_zone_current,
                ach_cooling_zone_current,
                ach_to_trigger_heating_zone_current,
            ) = zone.space_heat_cool_demand(
                delta_t_h,
                temp_ext_air,
                gains_internal_zone[z_name],
                gains_solar_zone[z_name],
                frac_convective_heat_zone_system[z_name][&h_name_list_sorted_zone[z_name][0]],
                frac_convective_cool_zone_system[z_name][&c_name_list_sorted_zone[z_name][0]],
                temp_setpnt_heat_zone_system[z_name][&h_name_list_sorted_zone[z_name][0]],
                temp_setpnt_cool_zone_system[z_name][&c_name_list_sorted_zone[z_name][0]],
                avg_air_supply_temp,
                None,
                None,
                AirChangesPerHourArgument::TargetAndWindowsOpen {
                    ach_target,
                    ach_windows_open,
                },
                simtime,
            )?;

            space_heat_demand_zone.insert(z_name.to_owned(), space_heat_demand_zone_current);
            space_cool_demand_zone.insert(z_name.to_owned(), space_cool_demand_zone_current);
            ach_cooling_zone.insert(z_name, ach_cooling_zone_current);
            ach_to_trigger_heating_zone.insert(z_name, ach_to_trigger_heating_zone_current);
        }

        // Ventilation required, including for cooling
        let is_heating_demand = space_heat_demand_zone.values().any(|&demand| demand > 0.0);
        let is_cooling_demand = space_cool_demand_zone.values().any(|&demand| demand < 0.0);

        let ach_cooling = if is_heating_demand {
            // Do not open windows any further than required for ventilation
            // requirement if there is any heating demand in any zone
            ach_target
        } else if is_cooling_demand {
            let ach_cooling = ach_target;

            // In this case, will need to recalculate space cooling demand for
            // each zone, this time assuming no window opening for all zones
            // TODO (from Python) There might be a way to make this more efficient and reduce
            //      the number of times the heat balance solver has to run, but
            //      this would require a wider refactoring of the zone module's
            //      space_heat_cool_demand function
            for (z_name, zone) in self.zones.iter() {
                let z_name = z_name.as_str();
                let (space_heat_demand_zone_current, space_cool_demand_zone_current, _, _) = zone
                    .space_heat_cool_demand(
                    delta_t_h,
                    temp_ext_air,
                    gains_internal_zone[z_name],
                    gains_solar_zone[z_name],
                    frac_convective_heat_zone_system[z_name][&h_name_list_sorted_zone[z_name][0]],
                    frac_convective_cool_zone_system[z_name][&c_name_list_sorted_zone[z_name][0]],
                    temp_setpnt_heat_zone_system[z_name][&h_name_list_sorted_zone[z_name][0]],
                    temp_setpnt_cool_zone_system[z_name][&c_name_list_sorted_zone[z_name][0]],
                    avg_air_supply_temp,
                    None,
                    None,
                    AirChangesPerHourArgument::Cooling { ach_cooling },
                    simtime,
                )?;
                space_heat_demand_zone.insert(z_name.to_owned(), space_heat_demand_zone_current);
                space_cool_demand_zone.insert(z_name.to_owned(), space_cool_demand_zone_current);
            }

            ach_cooling
        } else {
            // Subject to the above/below limits, take the maximum required window
            // opening from across all the zones
            let mut ach_cooling = *ach_cooling_zone
                .values()
                .max_by(|a, b| a.total_cmp(b))
                .unwrap();

            // Do not open windows to an extent where it would cause any zone
            // temperature to fall below the heating setpoint for that zone
            let ach_to_trigger_heating_list: Vec<f64> = ach_to_trigger_heating_zone
                .values()
                .filter_map(|x| *x)
                .collect();
            if !ach_to_trigger_heating_list.is_empty() {
                ach_cooling = ach_to_trigger_heating_list
                    .iter()
                    .max_by(|a, b| a.total_cmp(b).reverse())
                    .unwrap()
                    .min(ach_cooling);
            }

            // Do not reduce air change rate below ventilation requirement even
            // if it would help with temperature regulation

            ach_cooling.max(ach_target)
        };

        // Calculate heating/cooling system response and temperature achieved in each zone
        for (z_name, zone) in self.zones.iter() {
            let z_name = z_name.as_str();
            let c_name_list_hashset = c_name_list_sorted_zone[z_name]
                .iter()
                .collect::<HashSet<_>>();
            let h_name_list_hashset = h_name_list_sorted_zone[z_name]
                .iter()
                .collect::<HashSet<_>>();
            let intersection = h_name_list_hashset
                .intersection(&c_name_list_hashset)
                .collect::<Vec<_>>();
            if !intersection.is_empty() {
                bail!("All heating and cooling systems must have unique names")
            };
            // drop temporary values for working out intersection from scope
            drop(intersection);
            drop(c_name_list_hashset);
            drop(h_name_list_hashset);

            // TODO (from Python) Populate the lists below with minimum output of each heating and cooling system
            let mut hc_output_convective = h_name_list_sorted_zone[z_name]
                .iter()
                .chain(c_name_list_sorted_zone[z_name].iter())
                .map(|hc_name| (hc_name.as_str(), 0.0))
                .collect::<IndexMap<_, _>>();
            let mut hc_output_radiative = h_name_list_sorted_zone[z_name]
                .iter()
                .chain(c_name_list_sorted_zone[z_name].iter())
                .map(|hc_name| (hc_name.as_str(), 0.0))
                .collect::<IndexMap<_, _>>();
            let mut space_heat_demand_zone_system = h_name_list_sorted_zone[z_name]
                .iter()
                .map(|h_name| (h_name.as_str(), 0.0))
                .collect::<IndexMap<_, _>>();
            let mut space_cool_demand_zone_system = c_name_list_sorted_zone[z_name]
                .iter()
                .map(|c_name| (c_name.as_str(), 0.0))
                .collect::<IndexMap<_, _>>();
            let mut space_heat_provided_zone_system = h_name_list_sorted_zone[z_name]
                .iter()
                .map(|h_name| (h_name.as_str(), 0.0))
                .collect::<IndexMap<_, _>>();
            let mut space_cool_provided_zone_system = c_name_list_sorted_zone[z_name]
                .iter()
                .map(|c_name| (c_name.as_str(), 0.0))
                .collect::<IndexMap<_, _>>();

            let mut h_idx: usize = 0;
            let mut c_idx: usize = 0;
            let _space_heat_running_time_cumulative = 0.0;

            // these following variables need to be referenced outside the while loop with their last value retained

            let mut frac_convective_heat = 1.0;
            let mut frac_convective_cool = 1.0;

            while h_idx < h_name_list_sorted_zone[z_name].len()
                && c_idx < c_name_list_sorted_zone[z_name].len()
            {
                let h_name = &h_name_list_sorted_zone[z_name][h_idx].as_str();
                let c_name = &c_name_list_sorted_zone[z_name][c_idx].as_str();
                frac_convective_heat = frac_convective_heat_zone_system[z_name][h_name.to_owned()];
                frac_convective_cool = frac_convective_cool_zone_system[z_name][c_name.to_owned()];
                let temp_setpnt_heat = temp_setpnt_heat_zone_system[z_name][h_name.to_owned()];
                let temp_setpnt_cool = temp_setpnt_cool_zone_system[z_name][c_name.to_owned()];

                // Calculate space heating/cooling demand, accounting for any
                // output from systems (either output already calculated for
                // higher-priority systems, or min output of current or
                // lower-priority systems).
                // TODO (from Python) This doesn't yet handle minimum output from current or lower-priority systems
                let (gains_heat_cool_convective, gains_heat_cool_radiative) =
                    self.gains_heat_cool(delta_t_h, &hc_output_convective, &hc_output_radiative);
                if gains_heat_cool_convective == 0.0 && gains_heat_cool_radiative == 0.0 {
                    // If there is no output from any systems, then don't need to
                    // calculate demand again
                    space_heat_demand_zone_system.insert(h_name, space_heat_demand_zone[z_name]);
                    space_cool_demand_zone_system.insert(c_name, space_cool_demand_zone[z_name]);
                } else {
                    let (
                        space_heat_demand_zone_system_current,
                        space_cool_demand_zone_system_current,
                        ach_cooling_zone_current,
                        _,
                    ) = zone.space_heat_cool_demand(
                        delta_t_h,
                        temp_ext_air,
                        gains_internal_zone[z_name],
                        gains_solar_zone[z_name],
                        frac_convective_heat,
                        frac_convective_cool,
                        temp_setpnt_heat,
                        temp_setpnt_cool,
                        avg_air_supply_temp,
                        Some(gains_heat_cool_convective),
                        Some(gains_heat_cool_radiative),
                        AirChangesPerHourArgument::Cooling { ach_cooling },
                        simtime,
                    )?;
                    space_heat_demand_zone_system
                        .insert(h_name, space_heat_demand_zone_system_current);
                    space_cool_demand_zone_system
                        .insert(c_name, space_cool_demand_zone_system_current);
                    ach_cooling_zone.insert(z_name, ach_cooling_zone_current);
                }

                // If any heating systems potentially require overventilation,
                // calculate running time and throughput factor for current service
                // based on space heating demand assuming only overventilation
                // required for DHW

                // NB. there is a large chunk of Python in the upstream code here that is commented out there currently, so not brought over

                // Calculate heating/cooling provided
                if space_heat_demand_zone_system[h_name] > 0.0 {
                    space_heat_provided_zone_system.insert(
                        h_name,
                        self.space_heat_systems[h_name.to_owned()]
                            .lock()
                            .demand_energy(space_heat_demand_zone_system[h_name], simtime)?,
                    );
                    hc_output_convective.insert(
                        h_name,
                        space_heat_provided_zone_system[h_name] * frac_convective_heat,
                    );
                    hc_output_radiative.insert(
                        h_name,
                        space_heat_provided_zone_system[h_name] * (1.0 - frac_convective_heat),
                    );
                    // If heating has been provided, then next iteration of loop
                    // should use next-priority heating system
                    h_idx += 1;
                }
                if space_cool_demand_zone_system[c_name] < 0.0 {
                    space_cool_provided_zone_system.insert(
                        c_name,
                        self.space_cool_systems[c_name.to_owned()]
                            .demand_energy(space_cool_demand_zone_system[c_name], simtime),
                    );
                    hc_output_convective.insert(
                        c_name,
                        space_cool_provided_zone_system[c_name] * frac_convective_cool,
                    );
                    hc_output_radiative.insert(
                        c_name,
                        space_cool_provided_zone_system[c_name] * (1.0 - frac_convective_cool),
                    );
                    // If cooling has been provided, then next iteration of loop
                    // should use next-priority cooling system
                    c_idx += 1;
                }

                // Terminate loop if there is no more demand
                if space_heat_demand_zone_system[h_name] <= 0.0
                    && space_cool_demand_zone_system[c_name] >= 0.0
                {
                    break;
                }
            }

            // Call any remaining heating and cooling systems with zero demand
            for h_name in h_name_list_sorted_zone[z_name][h_idx..].iter() {
                let h_name = h_name.as_str();
                space_heat_provided_zone_system.insert(
                    h_name,
                    if !h_name.is_empty() {
                        self.space_heat_systems[h_name]
                            .lock()
                            .demand_energy(0.0, simtime)?
                    } else {
                        0.0
                    },
                );
                hc_output_convective.insert(
                    h_name,
                    space_heat_provided_zone_system[h_name] * frac_convective_heat,
                );
                hc_output_radiative.insert(
                    h_name,
                    space_heat_provided_zone_system[h_name] * (1.0 - frac_convective_heat),
                );
            }
            for c_name in c_name_list_sorted_zone[z_name][c_idx..].iter() {
                let c_name = c_name.as_str();
                space_cool_provided_zone_system.insert(
                    c_name,
                    if !c_name.is_empty() {
                        self.space_cool_systems[c_name].demand_energy(0.0, simtime)
                    } else {
                        0.0
                    },
                );
                hc_output_convective.insert(
                    c_name,
                    space_cool_provided_zone_system[c_name] * frac_convective_cool,
                );
                hc_output_radiative.insert(
                    c_name,
                    space_cool_provided_zone_system[c_name] * (1.0 - frac_convective_cool),
                );
            }
            // Calculate unmet demand
            self.unmet_demand(
                delta_t_h,
                temp_ext_air,
                z_name,
                zone,
                gains_internal_zone[z_name],
                gains_solar_zone[z_name],
                &temp_setpnt_heat_zone_system[z_name],
                &temp_setpnt_cool_zone_system[z_name],
                &frac_convective_heat_zone_system[z_name],
                &frac_convective_cool_zone_system[z_name],
                &h_name_list_sorted_zone[z_name],
                &c_name_list_sorted_zone[z_name],
                space_heat_demand_zone[z_name],
                space_cool_demand_zone[z_name],
                &hc_output_convective,
                &hc_output_radiative,
                ach_windows_open,
                ach_target,
                avg_air_supply_temp,
                simtime,
            )?;

            // Sum heating gains (+ve) and cooling gains (-ve) and convert from kWh to W
            let hc_output_convective_total = hc_output_convective.values().sum::<f64>();
            let hc_output_radiative_total = hc_output_radiative.values().sum::<f64>();

            let gains_heat_cool = (hc_output_convective_total + hc_output_radiative_total)
                * WATTS_PER_KILOWATT as f64
                / delta_t_h;
            let frac_convective = if gains_heat_cool != 0.0 {
                hc_output_convective_total
                    / (hc_output_convective_total + hc_output_radiative_total)
            } else {
                1.0
            };

            // Calculate final temperatures achieved
            heat_balance_map.insert(
                z_name.to_owned(),
                zone.update_temperatures(
                    delta_t,
                    temp_ext_air,
                    gains_internal_zone[z_name],
                    gains_solar_zone[z_name],
                    gains_heat_cool,
                    frac_convective,
                    ach_cooling,
                    avg_air_supply_temp,
                    simtime,
                ),
            );
            internal_air_temp.insert(z_name.to_owned(), zone.temp_internal_air());
            operative_temp.insert(z_name.to_owned(), zone.temp_operative());

            for h_name in h_name_list_sorted_zone[z_name].iter() {
                *space_heat_demand_system
                    .entry(h_name.to_owned())
                    .or_insert(0.0) += space_heat_demand_zone_system[h_name.as_str()];
                *space_heat_provided_system
                    .entry(h_name.to_owned())
                    .or_insert(0.0) += space_heat_provided_zone_system[h_name.as_str()];
            }
            for c_name in c_name_list_sorted_zone[z_name].iter() {
                *space_cool_demand_system
                    .entry(c_name.to_owned())
                    .or_insert(0.0) += space_cool_demand_zone_system[c_name.as_str()];
                *space_cool_provided_system
                    .entry(c_name.to_owned())
                    .or_insert(0.0) += space_cool_provided_zone_system[c_name.as_str()];
            }
        }

        Ok(SpaceHeatingCalculation {
            gains_internal_zone,
            gains_solar_zone,
            operative_temp,
            internal_air_temp,
            space_heat_demand_zone,
            space_cool_demand_zone,
            space_heat_demand_system,
            space_cool_demand_system,
            space_heat_provided_system,
            space_cool_provided_system,
            internal_gains_ductwork,
            heat_balance_map,
        })
    }

    /// Determine highest-priority system that is in its required heating
    /// or cooling period (and not just the setback period)
    ///
    /// Arguments:
    /// * `hc_name_list_sorted` - list of heating or cooling systems (not combined
    ///                                list), sorted in order of priority
    /// * `space_heat_cool_systems` - dict of space heating or cooling system objects
    ///                                   (not combined list)
    fn highest_priority_required_system(
        &self,
        hc_name_list_sorted: &[String],
        space_heat_cool_systems: SpaceHeatCoolSystems,
        simtime: SimulationTimeIteration,
    ) -> Option<String> {
        let mut hc_name_highest_req = Default::default();
        for hc_name in hc_name_list_sorted {
            if !hc_name.is_empty()
                && space_heat_cool_systems
                    .in_required_period_for_name(hc_name, simtime)
                    .unwrap_or(false)
            {
                hc_name_highest_req = Some(hc_name.to_owned());
                break;
            }
        }

        hc_name_highest_req
    }

    /// Calculate how much space heating / cooling demand is unmet
    fn unmet_demand(
        &self,
        delta_t_h: f64,
        temp_ext_air: f64,
        z_name: &str,
        zone: &Zone,
        gains_internal: f64,
        gains_solar: f64,
        temp_setpnt_heat_system: &IndexMap<String, f64>,
        temp_setpnt_cool_system: &IndexMap<String, f64>,
        frac_convective_heat_system: &IndexMap<String, f64>,
        frac_convective_cool_system: &IndexMap<String, f64>,
        h_name_list_sorted: &[String],
        c_name_list_sorted: &[String],
        space_heat_demand: f64,
        space_cool_demand: f64,
        hc_output_convective: &IndexMap<&str, f64>,
        hc_output_radiative: &IndexMap<&str, f64>,
        ach_max: f64,
        ach_target: f64,
        avg_air_supply_temp: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<()> {
        // Note: Use demand calculated based on highest-priority systems
        // Note: Demand is not considered unmet if it is outside the
        //       required heating/cooling period (which does not include
        //       times when the system is on due to setback or advanced
        //       start). If different systems have different required
        //       heating/cooling periods, unmet demand will be based on the
        //       system with the highest setpoint, ignoring any systems that
        //       are not in required periods (e.g. systems that are in
        //       setback or advanced start periods).
        // Note: Need to check that demand is non-zero, to avoid
        //       reporting unmet demand when heating system is absorbing
        //       energy from zone or cooling system is releasing energy
        //       to zone, which may be the case in some timesteps for
        //       systems with significant thermal mass.

        // Determine highest-priority system that is in its required heating
        // or cooling period (and not just the setback period)
        let h_name_highest_req = self.highest_priority_required_system(
            h_name_list_sorted,
            SpaceHeatCoolSystems::Heat(&self.space_heat_systems),
            simtime,
        );
        let c_name_highest_req = self.highest_priority_required_system(
            c_name_list_sorted,
            SpaceHeatCoolSystems::Cool(&self.space_cool_systems),
            simtime,
        );

        let gains_heat = h_name_list_sorted
            .iter()
            .map(|h_name| {
                hc_output_convective[h_name.as_str()] + hc_output_radiative[h_name.as_str()]
            })
            .sum::<f64>();
        let gains_cool = c_name_list_sorted
            .iter()
            .map(|c_name| {
                hc_output_convective[c_name.as_str()] + hc_output_radiative[c_name.as_str()]
            })
            .sum::<f64>();
        let energy_shortfall_heat = 0.0f64.max(space_heat_demand - gains_heat);
        let energy_shortfall_cool = 0.0f64.max(-(space_cool_demand - gains_cool));

        if (h_name_highest_req.is_some() && space_heat_demand > 0.0 && energy_shortfall_heat > 0.0)
            || (c_name_highest_req.is_some()
                && space_cool_demand < 0.0
                && energy_shortfall_cool > 0.0)
        {
            let (unmet_demand_heat, unmet_demand_cool) = if (energy_shortfall_heat > 0.0
                && h_name_highest_req
                    .as_ref()
                    .is_some_and(|h_name| h_name != &h_name_list_sorted[0]))
                || (energy_shortfall_cool > 0.0
                    && c_name_highest_req
                        .as_ref()
                        .is_some_and(|c_name| c_name != &c_name_list_sorted[0]))
            {
                // If the highest-priority system is not in required heating
                // period, but a lower-priority system is, calculate demand
                // based on the highest-priority system that is in required
                // heating period
                //
                // Handle case where no heating/cooling system is in required
                // period. In this case, there will be no heat output anyway so
                // the convective fraction doesn't matter
                let (frac_convective_heat, temp_setpnt_heat) =
                    if let Some(h_name) = h_name_highest_req.as_ref() {
                        (
                            frac_convective_heat_system[h_name],
                            temp_setpnt_heat_system[h_name],
                        )
                    } else {
                        (1.0, temp_setpnt_heat_none())
                    };
                let (frac_convective_cool, temp_setpnt_cool) =
                    if let Some(c_name) = c_name_highest_req.as_ref() {
                        (
                            frac_convective_cool_system[c_name],
                            temp_setpnt_cool_system[c_name],
                        )
                    } else {
                        (1.0, temp_setpnt_cool_none())
                    };

                let (space_heat_demand_req, space_cool_demand_req, _, _) = zone
                    .space_heat_cool_demand(
                        delta_t_h,
                        temp_ext_air,
                        gains_internal,
                        gains_solar,
                        frac_convective_heat,
                        frac_convective_cool,
                        temp_setpnt_heat,
                        temp_setpnt_cool,
                        avg_air_supply_temp,
                        None,
                        None,
                        AirChangesPerHourArgument::TargetAndWindowsOpen {
                            ach_target,
                            ach_windows_open: ach_max,
                        },
                        simtime,
                    )?;

                (
                    0.0f64.max(space_heat_demand_req - gains_heat),
                    0.0f64.max(-(space_cool_demand_req - gains_cool)),
                )
            } else {
                // If highest-priority system is in required heating period,
                // use the demand already calculated for the zone
                (energy_shortfall_heat, energy_shortfall_cool)
            };

            self.energy_supply_conn_unmet_demand_zone[z_name]
                .demand_energy(unmet_demand_heat + unmet_demand_cool, simtime.index)?;
        }

        Ok(())
    }

    pub fn run(&self) -> anyhow::Result<RunResults> {
        let simulation_time = self.simulation_time.as_ref().to_owned();
        let vec_capacity = || Vec::with_capacity(simulation_time.total_steps());

        let mut timestep_array = vec_capacity();
        let mut gains_internal_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut gains_solar_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut operative_temp_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut internal_air_temp_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut space_heat_demand_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut space_cool_demand_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut space_heat_demand_system_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut space_cool_demand_system_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut space_heat_provided_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut space_cool_provided_dict: IndexMap<String, Vec<f64>> = Default::default();
        let mut zone_list: Vec<KeyString> = Default::default();
        let mut hot_water_demand_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut hot_water_energy_demand_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut hot_water_energy_demand_dict_incl_pipework: IndexMap<KeyString, Vec<f64>> =
            Default::default();
        let mut hot_water_energy_output_dict: IndexMap<&str, Vec<f64>> = Default::default();
        let mut hot_water_duration_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut hot_water_no_events_dict: IndexMap<KeyString, Vec<usize>> = Default::default();
        let mut hot_water_pipework_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut ductwork_gains_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut hot_water_primary_pipework_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut hot_water_storage_losses_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut heat_balance_all_dict: HeatBalanceAllResults = IndexMap::from([
            (HeatBalanceFieldName::AirNode, Default::default()),
            (HeatBalanceFieldName::InternalBoundary, Default::default()),
            (HeatBalanceFieldName::ExternalBoundary, Default::default()),
        ]);
        let mut heat_source_wet_results_dict: IndexMap<String, ResultsPerTimestep> =
            Default::default();
        let mut heat_source_wet_results_annual_dict: IndexMap<String, ResultsAnnual> =
            Default::default();
        let mut emitters_output_dict: IndexMap<String, Vec<EmittersDetailedResult>> =
            Default::default();
        let mut vent_output_list: Vec<VentilationDetailedResult> = Default::default();

        for z_name in self.zones.keys() {
            let z_name = z_name.as_str().try_into().unwrap();
            gains_internal_dict.insert(z_name, vec_capacity());
            gains_solar_dict.insert(z_name, vec_capacity());
            operative_temp_dict.insert(z_name, vec_capacity());
            internal_air_temp_dict.insert(z_name, vec_capacity());
            space_heat_demand_dict.insert(z_name, vec_capacity());
            space_cool_demand_dict.insert(z_name, vec_capacity());
            zone_list.push(z_name);
            for heat_balance_value in heat_balance_all_dict.values_mut() {
                heat_balance_value.insert(z_name.to_string(), Default::default());
            }
        }

        for h_name in self.heat_system_name_for_zone.values() {
            space_heat_demand_system_dict.insert(h_name.to_owned(), vec_capacity());
            space_heat_provided_dict.insert(h_name.to_owned(), vec_capacity());
        }

        for c_name in self.cool_system_name_for_zone.values() {
            space_cool_demand_system_dict.insert(c_name.to_owned(), vec_capacity());
            space_cool_provided_dict.insert(c_name.to_owned(), vec_capacity());
        }

        hot_water_demand_dict.insert("demand".try_into().unwrap(), vec_capacity());
        hot_water_energy_demand_dict.insert("energy_demand".try_into().unwrap(), vec_capacity());
        hot_water_energy_demand_dict_incl_pipework.insert(
            "energy_demand_incl_pipework_loss".try_into().unwrap(),
            vec_capacity(),
        );
        hot_water_energy_output_dict.insert("energy_output", vec_capacity());
        hot_water_duration_dict.insert("duration".try_into().unwrap(), vec_capacity());
        hot_water_no_events_dict.insert(
            "no_events".try_into().unwrap(),
            Vec::with_capacity(simulation_time.total_steps()),
        );
        hot_water_pipework_dict.insert("pw_losses".try_into().unwrap(), vec_capacity());
        ductwork_gains_dict.insert("ductwork_gains".try_into().unwrap(), vec_capacity());
        hot_water_primary_pipework_dict
            .insert("primary_pw_losses".try_into().unwrap(), vec_capacity());
        hot_water_storage_losses_dict.insert("storage_losses".try_into().unwrap(), vec_capacity());
        self.initial_loop.store(true, Ordering::SeqCst);
        let mut internal_pressure_window: HashMap<ReportingFlag, f64> = Default::default();

        let delta_t_h = simulation_time.step_in_hours();

        #[cfg(feature = "indicatif")]
        let simulation_time_iter = simulation_time.progress();
        #[cfg(not(feature = "indicatif"))]
        let simulation_time_iter = simulation_time;

        for t_it in simulation_time_iter {
            timestep_array.push(t_it.time);
            let temp_hot_water = self.hot_water_sources["hw cylinder"].temp_hot_water();
            let _temp_final_drawoff = temp_hot_water;
            let _temp_average_drawoff = temp_hot_water;

            let DomesticHotWaterDemandData {
                hw_demand_vol,
                hw_demand_vol_target,
                hw_vol_at_tapping_points,
                hw_duration,
                all_events: no_events,
                hw_energy_demand,
                mut usage_events,
                vol_hot_water_equiv_elec_shower,
            } = self
                .domestic_hot_water_demand
                .hot_water_demand(t_it, temp_hot_water);

            let (hw_energy_output, pw_losses_internal, pw_losses_external, gains_internal_dhw_use) =
                if let HotWaterSource::StorageTank(storage_tank) =
                    &self.hot_water_sources["hw cylinder"]
                {
                    let (
                        hw_energy_output,
                        _unmet_energy,
                        temp_final_drawoff,
                        temp_average_drawoff,
                        volume_water_remove_from_tank,
                    ) = storage_tank
                        .lock()
                        .demand_hot_water(&mut usage_events, t_it)?;

                    let (pw_losses_internal, pw_losses_external, gains_internal_dhw_use) = self
                        .pipework_losses_and_internal_gains_from_hw_storage_tank(
                            delta_t_h,
                            volume_water_remove_from_tank,
                            hw_duration,
                            no_events,
                            temp_final_drawoff,
                            temp_average_drawoff,
                            temp_hot_water,
                            vol_hot_water_equiv_elec_shower,
                            t_it,
                        );

                    (
                        hw_energy_output,
                        pw_losses_internal,
                        pw_losses_external,
                        gains_internal_dhw_use,
                    )
                } else {
                    let hw_energy_output = self.hot_water_sources["hw cylinder"]
                        .demand_hot_water(hw_demand_vol_target, t_it);

                    let (pw_losses_internal, pw_losses_external, gains_internal_dhw_use) = self
                        .pipework_losses_and_internal_gains_from_hw(
                            delta_t_h,
                            hw_vol_at_tapping_points,
                            hw_duration,
                            no_events,
                            temp_hot_water,
                            t_it,
                        );

                    (
                        hw_energy_output,
                        pw_losses_internal,
                        pw_losses_external,
                        gains_internal_dhw_use,
                    )
                };

            // Convert from litres to kWh
            let cold_water_source = self.hot_water_sources["hw cylinder"].get_cold_water_source();
            let cold_water_temperature = cold_water_source.temperature(t_it);
            let hw_energy_demand_incl_pipework_loss =
                water_demand_to_kwh(hw_demand_vol, temp_hot_water, cold_water_temperature);
            let mut gains_internal_dhw = (pw_losses_internal + gains_internal_dhw_use)
                * WATTS_PER_KILOWATT as f64
                / t_it.timestep;
            match self.hot_water_sources.get("hw cylinder").unwrap() {
                HotWaterSource::StorageTank(ref source) => {
                    gains_internal_dhw += source.lock().internal_gains();
                }
                HotWaterSource::CombiBoiler(ref source) => {
                    gains_internal_dhw += source.internal_gains();
                }
                _ => {}
            }

            // loop through on-site energy generation
            for pv in self.on_site_generation.values() {
                // Get energy produced for the current timestep
                let (_energy_produced, energy_lost) = pv.produce_energy(t_it);
                // Add the energy lost figure to the internal gains if it is considered inside the building
                if pv.inverter_is_inside() {
                    gains_internal_dhw += energy_lost * WATTS_PER_KILOWATT as f64 / delta_t_h;
                }
            }

            // Addition of primary_pipework_losses_kWh for reporting as part of investigation of (upstream BRE) issue #31225: FDEV A082
            let (primary_pw_losses, storage_losses) =
                if let HotWaterSource::StorageTank(storage_tank) =
                    &self.hot_water_sources["hw cylinder"]
                {
                    storage_tank.lock().to_report()
                } else {
                    (0.0, 0.0)
                };

            let SpaceHeatingCalculation {
                gains_internal_zone,
                gains_solar_zone,
                operative_temp,
                internal_air_temp,
                space_heat_demand_zone,
                space_cool_demand_zone,
                space_heat_demand_system,
                space_cool_demand_system,
                space_heat_provided_system: space_heat_provided,
                space_cool_provided_system: space_cool_provided,
                internal_gains_ductwork: ductwork_gains,
                heat_balance_map: heat_balance_dict,
            } = self.calc_space_heating(
                t_it.timestep,
                gains_internal_dhw,
                &mut internal_pressure_window,
                t_it,
            )?;

            // Perform calculations that can only be done after all heating
            // services have been calculated
            for system in self.timestep_end_calcs.read().iter() {
                system.timestep_end(t_it);
            }

            for (z_name, gains_internal) in gains_internal_zone {
                gains_internal_dict
                    .get_mut(z_name.as_str())
                    .unwrap()
                    .push(gains_internal);
            }

            for (z_name, gains_solar) in gains_solar_zone {
                gains_solar_dict
                    .get_mut(z_name.as_str())
                    .unwrap()
                    .push(gains_solar);
            }

            for (z_name, temp) in operative_temp {
                operative_temp_dict
                    .get_mut(z_name.as_str())
                    .unwrap()
                    .push(temp);
            }

            for (z_name, temp) in internal_air_temp {
                internal_air_temp_dict
                    .get_mut(z_name.as_str())
                    .unwrap()
                    .push(temp);
            }

            for (z_name, demand) in space_heat_demand_zone {
                space_heat_demand_dict
                    .get_mut(z_name.as_str())
                    .unwrap()
                    .push(demand);
            }

            for (z_name, demand) in space_cool_demand_zone {
                space_cool_demand_dict
                    .get_mut(z_name.as_str())
                    .unwrap()
                    .push(demand);
            }

            for (h_name, demand) in space_heat_demand_system {
                space_heat_demand_system_dict
                    .get_mut(h_name.as_str())
                    .unwrap()
                    .push(demand);
            }

            for (c_name, demand) in space_cool_demand_system {
                space_cool_demand_system_dict
                    .get_mut(c_name.as_str())
                    .unwrap()
                    .push(demand);
            }

            for (h_name, output) in space_heat_provided {
                space_heat_provided_dict
                    .get_mut(&h_name)
                    .unwrap()
                    .push(output);
            }

            for (c_name, output) in space_cool_provided {
                space_cool_provided_dict
                    .get_mut(&c_name)
                    .unwrap()
                    .push(output);
            }

            for (z_name, hb_dict) in heat_balance_dict {
                if let Some(hb_dict) = hb_dict {
                    for (hb_name, gains_losses) in hb_dict.as_index_map() {
                        for (heat_gains_losses_name, heat_gains_losses_value) in gains_losses {
                            heat_balance_all_dict
                                .get_mut(&hb_name)
                                .unwrap()
                                .get_mut(&z_name)
                                .unwrap()
                                .entry(heat_gains_losses_name)
                                .or_default()
                                .push(heat_gains_losses_value);
                        }
                    }
                }
            }

            hot_water_demand_dict
                .get_mut("demand")
                .unwrap()
                .push(hw_demand_vol);
            hot_water_energy_demand_dict
                .get_mut("energy_demand")
                .unwrap()
                .push(hw_energy_demand);
            hot_water_energy_demand_dict_incl_pipework
                .get_mut("energy_demand_incl_pipework_loss")
                .unwrap()
                .push(hw_energy_demand_incl_pipework_loss);
            hot_water_energy_output_dict
                .get_mut("energy_output")
                .unwrap()
                .push(hw_energy_output);
            hot_water_duration_dict
                .get_mut("duration")
                .unwrap()
                .push(hw_duration);
            hot_water_no_events_dict
                .get_mut("no_events")
                .unwrap()
                .push(no_events);
            hot_water_pipework_dict
                .get_mut("pw_losses")
                .unwrap()
                .push(pw_losses_internal + pw_losses_external);
            ductwork_gains_dict
                .get_mut("ductwork_gains")
                .unwrap()
                .push(ductwork_gains);
            hot_water_primary_pipework_dict
                .get_mut("primary_pw_losses")
                .unwrap()
                .push(primary_pw_losses);
            hot_water_storage_losses_dict
                .get_mut("storage_losses")
                .unwrap()
                .push(storage_losses);

            self.energy_supplies
                .values()
                .map(|supply| Ok(supply.read().calc_energy_import_export_betafactor(t_it)?))
                .collect::<anyhow::Result<()>>()?;

            for diverter in &self.diverters {
                diverter.write().timestep_end();
            }
        }

        // Return results from all energy supplies
        let mut results_totals: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut results_end_user: IndexMap<KeyString, IndexMap<String, Vec<f64>>> =
            Default::default();
        let mut energy_import: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut energy_export: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut energy_generated_consumed: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut energy_to_storage: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut energy_from_storage: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut energy_diverted: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut betafactor: IndexMap<KeyString, Vec<f64>> = Default::default();
        for (name, supply) in self.energy_supplies.iter() {
            let supply = supply.read();
            let name: KeyString = KeyString::from(name)?;
            results_totals.insert(name, supply.results_total());
            results_end_user.insert(name, supply.results_by_end_user().to_owned());
            energy_import.insert(name, supply.get_energy_import().to_owned());
            energy_export.insert(name, supply.get_energy_export().to_owned());
            energy_generated_consumed
                .insert(name, supply.get_energy_generated_consumed().to_owned());
            let (energy_to, energy_from) = supply.get_energy_to_from_battery();
            energy_to_storage.insert(name, energy_to.to_owned());
            energy_from_storage.insert(name, energy_from.to_owned());
            energy_diverted.insert(name, supply.get_energy_diverted().to_owned());
            betafactor.insert(name, supply.get_beta_factor().to_owned());
        }

        let hot_water_energy_out: IndexMap<String, Vec<f64>> = IndexMap::from([(
            "hw cylinder".to_owned(),
            hot_water_energy_output_dict
                .get("energy_output")
                .unwrap()
                .to_owned(),
        )]);
        let dhw_cop_dict = self.heat_cool_cop(
            &hot_water_energy_out,
            &results_end_user,
            self.energy_supply_conn_names_for_hot_water_source.clone(),
        );
        let heat_cop_dict = self.heat_cool_cop(
            &space_heat_provided_dict,
            &results_end_user,
            self.energy_supply_conn_names_for_heat_systems
                .iter()
                .map(|(k, v)| (k.clone(), vec![v.clone()]))
                .collect(),
        );
        let cool_cop_dict = self.heat_cool_cop(
            &space_cool_provided_dict,
            &results_end_user,
            self.space_cool_systems
                .keys()
                .map(|system_name| (system_name.clone(), vec![system_name.clone()]))
                .collect(),
        );

        let zone_dict = IndexMap::from([
            (ZoneResultKey::InternalGains, gains_internal_dict),
            (ZoneResultKey::SolarGains, gains_solar_dict),
            (ZoneResultKey::OperativeTemp, operative_temp_dict),
            (ZoneResultKey::InternalAirTemp, internal_air_temp_dict),
            (ZoneResultKey::SpaceHeatDemand, space_heat_demand_dict),
            (ZoneResultKey::SpaceCoolDemand, space_cool_demand_dict),
        ]);

        let hc_system_dict = IndexMap::from([
            (
                HeatingCoolingSystemResultKey::HeatingSystem,
                space_heat_demand_system_dict,
            ),
            (
                HeatingCoolingSystemResultKey::HeatingSystemOutput,
                space_heat_provided_dict,
            ),
            (
                HeatingCoolingSystemResultKey::CoolingSystem,
                space_cool_demand_system_dict,
            ),
            (
                HeatingCoolingSystemResultKey::CoolingSystemOutput,
                space_cool_provided_dict,
            ),
        ]);

        let hot_water_dict = IndexMap::from([
            (
                HotWaterResultKey::HotWaterDemand,
                HotWaterResultMap::Float(hot_water_demand_dict),
            ),
            (
                HotWaterResultKey::HotWaterEnergyDemandIncludingPipeworkLoss,
                HotWaterResultMap::Float(hot_water_energy_demand_dict_incl_pipework),
            ),
            (
                HotWaterResultKey::HotWaterEnergyDemand,
                HotWaterResultMap::Float(hot_water_energy_demand_dict),
            ),
            (
                HotWaterResultKey::HotWaterDuration,
                HotWaterResultMap::Float(hot_water_duration_dict),
            ),
            (
                HotWaterResultKey::HotWaterEvents,
                HotWaterResultMap::Int(hot_water_no_events_dict),
            ),
            (
                HotWaterResultKey::PipeworkLosses,
                HotWaterResultMap::Float(hot_water_pipework_dict),
            ),
            (
                HotWaterResultKey::PrimaryPipeworkLosses,
                HotWaterResultMap::Float(hot_water_primary_pipework_dict),
            ),
            (
                HotWaterResultKey::StorageLosses,
                HotWaterResultMap::Float(hot_water_storage_losses_dict),
            ),
        ]);

        // Report detailed outputs from heat source wet objects, if requested and available
        // TODO (from Python) Note that the below assumes that there is only one water
        //      heating service and therefore that all hot water energy
        //      output is assigned to that service. If the model changes in
        //      future to allow more than one hot water system, this code may
        //      need to be revised to handle that scenario.
        if self.detailed_output_heating_cooling {
            for (name, heat_source_wet) in self.wet_heat_sources.iter() {
                if let Some((results, results_annual)) = heat_source_wet
                    .output_detailed_results(&hot_water_energy_output_dict["energy_output"])
                {
                    heat_source_wet_results_dict.insert(name.clone(), results);
                    heat_source_wet_results_annual_dict.insert(name.clone(), results_annual);
                }
            }

            // Emitter detailed output results are stored with respect to heat_system_name
            for (heat_system_name, heat_system) in self.space_heat_systems.iter() {
                if let Some(emitters_output) = heat_system.lock().output_emitter_results() {
                    emitters_output_dict.insert(heat_system_name.to_owned(), emitters_output);
                }
            }

            vent_output_list = self.ventilation.output_vent_results().read().clone();
        }

        Ok(RunResults {
            timestep_array,
            results_totals,
            results_end_user,
            energy_import,
            energy_export,
            energy_generated_consumed,
            energy_to_storage,
            energy_from_storage,
            energy_diverted,
            betafactor,
            zone_dict,
            zone_list,
            hc_system_dict,
            hot_water_dict,
            heat_cop_dict,
            cool_cop_dict,
            dhw_cop_dict,
            ductwork_gains: ductwork_gains_dict,
            heat_balance_dict: heat_balance_all_dict,
            heat_source_wet_results_dict,
            heat_source_wet_results_annual_dict,
            emitters_output_dict,
            vent_output_list,
        })
    }

    /// Calculate overall CoP over calculation period for each heating and cooling system
    fn heat_cool_cop(
        &self,
        energy_provided: &IndexMap<String, Vec<f64>>,
        results_end_user: &IndexMap<KeyString, IndexMap<String, Vec<f64>>>,
        energy_supply_conn_name_for_space_hc_system: IndexMap<String, Vec<String>>,
    ) -> IndexMap<String, NumberOrDivisionByZero> {
        let mut hc_output_overall: IndexMap<String, f64> = Default::default();
        let mut hc_input_overall: IndexMap<String, f64> = Default::default();
        let mut cop_dict: IndexMap<String, NumberOrDivisionByZero> = Default::default();
        for (hc_name, hc_output) in energy_provided {
            hc_output_overall.insert(hc_name.to_owned(), hc_output.iter().sum::<f64>().abs());
            hc_input_overall.insert(hc_name.to_owned(), 0.);
            let energy_supply_conn_names =
                match energy_supply_conn_name_for_space_hc_system.get(hc_name.as_str()) {
                    Some(hc_name) => hc_name.clone(),
                    None => vec![],
                };
            for (fuel_name, fuel_summary) in results_end_user {
                if fuel_name == "_unmet_demand" {
                    continue;
                }
                for (conn_name, energy_cons) in fuel_summary {
                    if energy_supply_conn_names.contains(conn_name) {
                        *hc_input_overall.get_mut(hc_name).unwrap() +=
                            energy_cons.iter().sum::<f64>();
                    }
                }
            }

            cop_dict.insert(
                hc_name.to_owned(),
                if hc_input_overall[hc_name] > 0. {
                    NumberOrDivisionByZero::Number(
                        hc_output_overall[hc_name] / hc_input_overall[hc_name],
                    )
                } else {
                    NumberOrDivisionByZero::DivisionByZero
                },
            );
        }

        cop_dict
    }
}

#[derive(Debug, Default)]
pub struct OutputOptions {
    print_heat_balance: bool,
    detailed_output_heating_cooling: bool,
}

impl From<&ProjectFlags> for OutputOptions {
    fn from(flags: &ProjectFlags) -> Self {
        Self {
            print_heat_balance: flags.contains(ProjectFlags::HEAT_BALANCE),
            detailed_output_heating_cooling: flags
                .contains(ProjectFlags::DETAILED_OUTPUT_HEATING_COOLING),
        }
    }
}

struct SetpointsAndConvectiveFractions {
    temp_setpnt_heat: IndexMap<String, f64>,
    temp_setpnt_cool: IndexMap<String, f64>,
    frac_convective_heat: IndexMap<String, f64>,
    frac_convective_cool: IndexMap<String, f64>,
}

#[derive(Clone, Debug)]
pub enum HotWaterResultMap {
    Float(IndexMap<KeyString, Vec<f64>>),
    Int(IndexMap<KeyString, Vec<usize>>),
}

#[derive(Clone, Copy, Debug)]
pub enum NumberOrDivisionByZero {
    Number(f64),
    DivisionByZero,
}

impl Display for NumberOrDivisionByZero {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                NumberOrDivisionByZero::Number(number) => Cow::Owned(format!("{}", number)),
                NumberOrDivisionByZero::DivisionByZero => Cow::Borrowed("DIV/0"),
            }
        )
    }
}

pub type ResultsEndUser = IndexMap<KeyString, IndexMap<String, Vec<f64>>>;

enum SpaceHeatCoolSystems<'a> {
    Heat(&'a IndexMap<String, Arc<Mutex<SpaceHeatSystem>>>),
    Cool(&'a IndexMap<String, AirConditioning>),
}

impl SpaceHeatCoolSystems<'_> {
    fn in_required_period_for_name(
        &self,
        system_name: &str,
        simtime: SimulationTimeIteration,
    ) -> Option<bool> {
        match self {
            SpaceHeatCoolSystems::Heat(heat) => {
                heat[system_name].lock().in_required_period(simtime)
            }
            SpaceHeatCoolSystems::Cool(cool) => cool[system_name].in_required_period(&simtime),
        }
    }
}

fn has_unique_non_default_values<K, V: Default + Eq + Hash>(map: &IndexMap<K, V>) -> bool {
    let values: Vec<&V> = map.values().filter(|x| x != &&V::default()).collect();
    let value_set: HashSet<&&V> = values.iter().collect();
    values.len() == value_set.len()
}

fn external_conditions_from_input(
    input: Arc<ExternalConditionsInput>,
    simulation_time: &SimulationTimeIterator,
) -> ExternalConditions {
    ExternalConditions::new(
        simulation_time,
        input.air_temperatures.clone().unwrap_or_default(),
        input.wind_speeds.clone().unwrap_or_default(),
        input.wind_directions.clone().unwrap_or_default(),
        input
            .diffuse_horizontal_radiation
            .clone()
            .unwrap_or_default(),
        input.direct_beam_radiation.clone().unwrap_or_default(),
        input
            .solar_reflectivity_of_ground
            .clone()
            .unwrap_or_default(),
        input.latitude.unwrap_or(51.5), // default to London (though we expect this to be set!)
        input.longitude.unwrap_or(-0.13), // default to London (though we expect this to be set!)
        input.timezone.unwrap_or(0),    // default to GMT/ UTC
        input.start_day.unwrap_or(0),
        input.end_day,
        input.time_series_step.unwrap_or(1.0),
        input.january_first,
        input.daylight_savings,
        input.leap_day_included.unwrap_or(false),
        input.direct_beam_conversion_needed.unwrap_or(false),
        input.shading_segments.clone(),
    )
}

pub(crate) type ColdWaterSources = IndexMap<ColdWaterSourceType, ColdWaterSource>;

fn cold_water_sources_from_input(input: &ColdWaterSourceInput) -> ColdWaterSources {
    input
        .iter()
        .map(|(source_type, source_details)| {
            (
                *source_type,
                cold_water_source_from_input_details(source_details),
            )
        })
        .collect()
}

fn cold_water_source_from_input_details(details: &ColdWaterSourceDetails) -> ColdWaterSource {
    ColdWaterSource::new(
        details.temperatures.clone(),
        details.start_day,
        details.time_series_step,
    )
}

fn energy_supplies_from_input(
    input: &EnergySupplyInput,
    simulation_time_iterator: &SimulationTimeIterator,
    external_conditions: Arc<ExternalConditions>,
) -> IndexMap<String, Arc<RwLock<EnergySupply>>> {
    let mut supplies: IndexMap<String, Arc<RwLock<EnergySupply>>> = input
        .iter()
        .map(|(name, supply)| {
            (
                name.to_owned(),
                energy_supply_from_input(
                    supply,
                    simulation_time_iterator,
                    external_conditions.clone(),
                ),
            )
        })
        .collect();
    // set up supply representing unmet demand
    supplies.insert(
        UNMET_DEMAND_SUPPLY_NAME.to_string(),
        Arc::new(RwLock::new(EnergySupply::new(
            FuelType::UnmetDemand,
            0,
            None,
            None,
            None,
        ))),
    );

    supplies
}

fn energy_supply_from_input(
    input: &EnergySupplyDetails,
    simulation_time_iterator: &SimulationTimeIterator,
    external_conditions: Arc<ExternalConditions>,
) -> Arc<RwLock<EnergySupply>> {
    Arc::new(RwLock::new(EnergySupply::new(
        input.fuel,
        simulation_time_iterator.total_steps(),
        input.electric_battery.as_ref().map(|battery_input| {
            ElectricBattery::from_input(
                battery_input,
                simulation_time_iterator.step_in_hours(),
                external_conditions,
            )
        }),
        input.priority.as_ref().cloned(),
        input.is_export_capable,
    )))
}

type DiverterTypes = IndexMap<String, EnergyDiverter>;

// struct DiverterTypes {
//     pub mains_electricity: Option<EnergyDiverter>,
//     pub mains_gas: Option<EnergyDiverter>,
//     pub bulk_lpg: Option<EnergyDiverter>,
// }

// impl From<&EnergySupplyInput> for DiverterTypes {
//     fn from(input: &EnergySupplyInput) -> Self {
//         Self {
//             mains_electricity: diverter_from_energy_supply(
//                 input.get(&EnergySupplyKey::MainsElectricity),
//             ),
//             mains_gas: diverter_from_energy_supply(input.get(&EnergySupplyKey::MainsGas)),
//             bulk_lpg: diverter_from_energy_supply(input.get(&EnergySupplyKey::BulkLpg)),
//         }
//     }
// }
//
// impl DiverterTypes {
//     pub fn get_for_supply_type(
//         &self,
//         energy_supply_type: EnergySupplyType,
//     ) -> Option<&EnergyDiverter> {
//         match energy_supply_type {
//             EnergySupplyType::Electricity => self.mains_electricity.as_ref(),
//             EnergySupplyType::MainsGas => self.mains_gas.as_ref(),
//             EnergySupplyType::LpgBulk => self.bulk_lpg.as_ref(),
//             _ => None,
//         }
//     }
// }

fn diverter_from_energy_supply(supply: &EnergySupplyDetails) -> Option<EnergyDiverter> {
    supply.diverter.as_ref().map(|diverter| diverter.clone())
}

pub type InternalGainsCollection = IndexMap<String, Gains>;

fn internal_gains_from_input(
    input: &InternalGainsInput,
    total_floor_area: f64,
) -> anyhow::Result<InternalGainsCollection> {
    let mut gains_collection = InternalGainsCollection::from([]);
    if let Some(internal_gains) = input.total_internal_gains.as_ref() {
        gains_collection.insert(
            "total_internal_gains".to_string(),
            Gains::Internal(internal_gains_from_details(
                internal_gains,
                total_floor_area,
            )?),
        );
    }
    if let Some(internal_gains) = input.metabolic_gains.as_ref() {
        gains_collection.insert(
            "metabolic_gains".to_string(),
            Gains::Internal(internal_gains_from_details(
                internal_gains,
                total_floor_area,
            )?),
        );
    }
    if let Some(internal_gains) = input.evaporative_losses.as_ref() {
        gains_collection.insert(
            "evaporative_losses".to_string(),
            Gains::Internal(internal_gains_from_details(
                internal_gains,
                total_floor_area,
            )?),
        );
    }
    if let Some(internal_gains) = input.other.as_ref() {
        gains_collection.insert(
            "other".to_string(),
            Gains::Internal(internal_gains_from_details(
                internal_gains,
                total_floor_area,
            )?),
        );
    }

    Ok(gains_collection)
}

fn internal_gains_from_details(
    details: &InternalGainsDetails,
    total_floor_area: f64,
) -> anyhow::Result<InternalGains> {
    Ok(InternalGains::new(
        convert_energy_to_wm2(details, total_floor_area)?,
        details.start_day,
        details.time_series_step,
    ))
}

fn convert_energy_to_wm2(
    internal_gains_details: &InternalGainsDetails,
    total_floor_area: f64,
) -> anyhow::Result<Vec<f64>> {
    let schedule = &internal_gains_details.schedule;
    Ok(reject_nulls(expand_numeric_schedule(schedule))?
        .iter()
        .map(|energy_data| energy_data / total_floor_area)
        .collect())
}

#[derive(Debug)]
pub struct Controls {
    core: Vec<HeatSourceControl>,
    extra: HashMap<String, Arc<Control>>,
}

impl Controls {
    pub fn new(core: Vec<HeatSourceControl>, extra: HashMap<String, Arc<Control>>) -> Self {
        Self { core, extra }
    }

    pub fn get(&self, control_type: &HeatSourceControlType) -> Option<Arc<Control>> {
        self.core
            .iter()
            .find(|heat_source_control| heat_source_control.has_type(*control_type))
            .map(|heat_source_control| heat_source_control.get())
    }

    // access a control using a string, possibly because it is one of the "extra" controls
    pub fn get_with_string(&self, control_name: &str) -> Option<Arc<Control>> {
        match control_name {
            // hard-code ways of resolving to core control types (for now)
            "hw timer" => self.get(&HeatSourceControlType::HotWaterTimer),
            "window opening" => self.get(&HeatSourceControlType::WindowOpening),
            other => self.extra.get(other).cloned(),
        }
    }
}

fn control_from_input(
    control_input: &ControlInput,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<Controls> {
    let mut core: Vec<HeatSourceControl> = Default::default();
    let mut extra: HashMap<String, Arc<Control>> = Default::default();

    // this is very ugly(!) but is just a reflection of the lack of clarity in the schema
    // and the way the variants-struct crate works;
    // we should be able to improve it in time
    if let Some(control) = &control_input.hot_water_timer.as_ref() {
        core.push(HeatSourceControl::HotWaterTimer(Arc::new(
            single_control_from_details(control, simulation_time_iterator)?,
        )));
    }
    if let Some(control) = &control_input.window_opening.as_ref() {
        core.push(HeatSourceControl::WindowOpening(Arc::new(
            single_control_from_details(control, simulation_time_iterator)?,
        )));
    }
    for (name, control) in &control_input.extra {
        extra.insert(
            name.to_string(),
            Arc::new(single_control_from_details(
                control,
                simulation_time_iterator,
            )?),
        );
    }

    Ok(Controls { core, extra })
}

fn single_control_from_details(
    details: &ControlDetails,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<Control> {
    Ok(match details {
        ControlDetails::OnOffTime {
            start_day,
            time_series_step,
            schedule,
            ..
        } => Control::OnOffTimeControl(OnOffTimeControl::new(
            reject_nulls(expand_boolean_schedule(schedule))?,
            *start_day,
            *time_series_step,
        )),
        ControlDetails::OnOffCostMinimisingTime {
            start_day,
            time_series_step,
            time_on_daily,
            schedule,
            ..
        } => Control::OnOffMinimisingTimeControl(OnOffMinimisingTimeControl::new(
            reject_nulls(expand_numeric_schedule(schedule))?,
            *start_day,
            *time_series_step,
            time_on_daily.unwrap_or_default(),
        )),
        ControlDetails::SetpointTime {
            start_day,
            time_series_step,
            advanced_start,
            setpoint_min,
            setpoint_max,
            default_to_max,
            schedule,
            ..
        } => Control::SetpointTimeControl(
            SetpointTimeControl::new(
                expand_numeric_schedule(schedule),
                *start_day,
                *time_series_step,
                *setpoint_min,
                *setpoint_max,
                *default_to_max,
                *advanced_start,
                simulation_time_iterator.step_in_hours(),
            )
            .unwrap(),
        ),
        ControlDetails::Charge {
            start_day,
            time_series_step,
            charge_level,
            schedule,
            ..
        } => {
            // Simulation manual charge control
            // Set charge level to 1.0 (max) for each day of simulation (plus 1)
            let vec_size = ((simulation_time_iterator.total_steps() as f64
                * simulation_time_iterator.step_in_hours()
                / 24.0)
                + 1.0)
                .ceil() as usize;
            let mut charge_level_vec: Vec<f64> = vec![1.0; vec_size];
            // if charge_level is present in input, overwrite initial vector
            // user can specify a vector with all days (plus 1), or as a single float value to be used for each day
            if let Some(charge) = charge_level {
                match charge {
                    ChargeLevel::List(charge_vec) => {
                        charge_level_vec = charge_vec.to_vec();
                    }
                    ChargeLevel::Single(charge) => {
                        charge_level_vec = vec![*charge; vec_size];
                    }
                    ChargeLevel::Schedule(_) => {
                        unimplemented!("TODO 0.32 add support for schedules in charge level")
                    }
                }
            }

            Control::ToUChargeControl(ToUChargeControl {
                schedule: reject_nulls(expand_boolean_schedule(schedule))?,
                start_day: *start_day,
                time_series_step: *time_series_step,
                charge_level: charge_level_vec,
            })
        }
        ControlDetails::CombinationTime { .. } => {
            unimplemented!("CombinationTime control to be implemented during 0.32 migration")
        }
        ControlDetails::SmartAppliance { .. } => {
            unimplemented!("SmartAppliance control to be implemented during 0.32 migration")
        }
    })
}

fn wwhrs_from_input(
    wwhrs: Option<&WasteWaterHeatRecovery>,
    cold_water_sources: &ColdWaterSources,
    initial_simtime: SimulationTimeIteration,
) -> IndexMap<String, Arc<Mutex<Wwhrs>>> {
    let mut wwhr_systems: IndexMap<String, Arc<Mutex<Wwhrs>>> = IndexMap::from([]);
    if let Some(systems) = wwhrs {
        for (name, system) in systems {
            wwhr_systems
                .entry(name.clone())
                .or_insert(Arc::new(Mutex::new(wwhr_system_from_details(
                    system.clone(),
                    cold_water_sources,
                    initial_simtime,
                ))));
        }
    }

    wwhr_systems
}

fn wwhr_system_from_details(
    system: WasteWaterHeatRecoveryDetails,
    cold_water_sources: &ColdWaterSources,
    initial_simtime: SimulationTimeIteration,
) -> Wwhrs {
    match system.system_type {
        WwhrsType::SystemA => Wwhrs::WWHRSInstantaneousSystemA(WWHRSInstantaneousSystemA::new(
            system.flow_rates,
            system.efficiencies,
            get_cold_water_source_ref_for_type(system.cold_water_source, cold_water_sources)
                .unwrap(),
            system.utilisation_factor,
            initial_simtime,
        )),
        WwhrsType::SystemB => Wwhrs::WWHRSInstantaneousSystemB(WWHRSInstantaneousSystemB::new(
            get_cold_water_source_ref_for_type(system.cold_water_source, cold_water_sources)
                .unwrap(),
            system.flow_rates,
            system.efficiencies,
            system.utilisation_factor,
        )),
        WwhrsType::SystemC => Wwhrs::WWHRSInstantaneousSystemC(WWHRSInstantaneousSystemC::new(
            system.flow_rates,
            system.efficiencies,
            get_cold_water_source_ref_for_type(system.cold_water_source, cold_water_sources)
                .unwrap(),
            system.utilisation_factor,
            initial_simtime,
        )),
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) enum ReportingFlag {
    Min,
    Max,
}

impl Display for ReportingFlag {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                ReportingFlag::Min => "min",
                ReportingFlag::Max => "max",
            }
        )
    }
}

fn temp_internal_air_for_zones(zones: Arc<IndexMap<String, Zone>>, total_volume: f64) -> f64 {
    let internal_air_temperature = zones
        .values()
        .map(|zone| zone.temp_internal_air() * zone.volume())
        .sum::<f64>();

    internal_air_temperature / total_volume
}

/// This is a struct that encapsulates a shared part of the corpus, just enough to be able to provide temperature of internal air
/// as an equivalent to Corpus::temp_internal_air
#[derive(Clone, Debug)]
pub struct TempInternalAirAccessor {
    pub zones: Arc<IndexMap<String, Zone>>,
    pub total_volume: f64,
}

impl TempInternalAirAccessor {
    pub fn call(&self) -> f64 {
        temp_internal_air_for_zones(self.zones.clone(), self.total_volume)
    }
}

pub(crate) type TempInternalAirFn = Arc<dyn Fn() -> f64 + Send + Sync>;

fn temp_internal_air_fn(accessor: impl Into<Arc<TempInternalAirAccessor>>) -> TempInternalAirFn {
    let accessor: Arc<TempInternalAirAccessor> = accessor.into();
    Arc::new(move || accessor.call())
}

pub type KeyString = ArrayString<64>;

/// A struct definition to encapsulate results from a corpus run.
pub struct RunResults {
    pub(crate) timestep_array: Vec<f64>,
    pub(crate) results_totals: IndexMap<KeyString, Vec<f64>>,
    pub(crate) results_end_user: ResultsEndUser,
    pub(crate) energy_import: IndexMap<KeyString, Vec<f64>>,
    pub(crate) energy_export: IndexMap<KeyString, Vec<f64>>,
    pub(crate) energy_generated_consumed: IndexMap<KeyString, Vec<f64>>,
    pub(crate) energy_to_storage: IndexMap<KeyString, Vec<f64>>,
    pub(crate) energy_from_storage: IndexMap<KeyString, Vec<f64>>,
    pub(crate) energy_diverted: IndexMap<KeyString, Vec<f64>>,
    pub(crate) betafactor: IndexMap<KeyString, Vec<f64>>,
    pub(crate) zone_dict: IndexMap<ZoneResultKey, IndexMap<KeyString, Vec<f64>>>,
    pub(crate) zone_list: Vec<KeyString>,
    pub(crate) hc_system_dict: IndexMap<HeatingCoolingSystemResultKey, IndexMap<String, Vec<f64>>>,
    pub(crate) hot_water_dict: IndexMap<HotWaterResultKey, HotWaterResultMap>,
    pub(crate) heat_cop_dict: IndexMap<String, NumberOrDivisionByZero>,
    pub(crate) cool_cop_dict: IndexMap<String, NumberOrDivisionByZero>,
    pub(crate) dhw_cop_dict: IndexMap<String, NumberOrDivisionByZero>,
    pub(crate) ductwork_gains: IndexMap<KeyString, Vec<f64>>,
    pub(crate) heat_balance_dict: HeatBalanceAllResults,
    pub(crate) heat_source_wet_results_dict: IndexMap<String, ResultsPerTimestep>,
    pub(crate) heat_source_wet_results_annual_dict: IndexMap<String, ResultsAnnual>,
    pub(crate) emitters_output_dict: IndexMap<String, Vec<EmittersDetailedResult>>,
    pub(crate) vent_output_list: Vec<VentilationDetailedResult>,
}

impl RunResults {
    pub(crate) fn space_heat_demand_total(&self) -> f64 {
        self.zone_dict[&ZoneResultKey::SpaceHeatDemand]
            .values()
            .map(|v| v.iter().sum::<f64>())
            .sum::<f64>()
    }

    pub(crate) fn space_cool_demand_total(&self) -> f64 {
        self.zone_dict[&ZoneResultKey::SpaceCoolDemand]
            .values()
            .map(|v| v.iter().sum::<f64>())
            .sum::<f64>()
    }
}

pub(crate) type HeatBalanceAllResults =
    IndexMap<HeatBalanceFieldName, IndexMap<String, IndexMap<String, Vec<f64>>>>;

struct SpaceHeatingCalculation {
    gains_internal_zone: HashMap<String, f64>,
    gains_solar_zone: HashMap<String, f64>,
    operative_temp: HashMap<String, f64>,
    internal_air_temp: HashMap<String, f64>,
    space_heat_demand_zone: HashMap<String, f64>,
    space_cool_demand_zone: HashMap<String, f64>,
    space_heat_demand_system: HashMap<String, f64>,
    space_cool_demand_system: HashMap<String, f64>,
    space_heat_provided_system: HashMap<String, f64>,
    space_cool_provided_system: HashMap<String, f64>,
    internal_gains_ductwork: f64,
    heat_balance_map: HashMap<String, Option<HeatBalance>>,
}

#[derive(Clone, Copy, Deserialize_enum_str, Debug, Eq, Hash, PartialEq, Serialize_enum_str)]
pub(crate) enum ZoneResultKey {
    #[serde(rename = "internal gains")]
    InternalGains,
    #[serde(rename = "solar gains")]
    SolarGains,
    #[serde(rename = "operative temp")]
    OperativeTemp,
    #[serde(rename = "internal air temp")]
    InternalAirTemp,
    #[serde(rename = "space heat demand")]
    SpaceHeatDemand,
    #[serde(rename = "space cool demand")]
    SpaceCoolDemand,
}

impl ZoneResultKey {
    pub(crate) fn as_str(&self) -> &'static str {
        // if there's a clean way of reusing the serde deserialization map above, this would be preferable
        // but for now, this minor duplication will do
        match self {
            ZoneResultKey::InternalGains => "internal gains",
            ZoneResultKey::SolarGains => "solar gains",
            ZoneResultKey::OperativeTemp => "operative temp",
            ZoneResultKey::InternalAirTemp => "internal air temp",
            ZoneResultKey::SpaceHeatDemand => "space heat demand",
            ZoneResultKey::SpaceCoolDemand => "space cool demand",
        }
    }
}

#[derive(Clone, Copy, Deserialize_enum_str, Debug, Eq, Hash, PartialEq, Serialize_enum_str)]
pub(crate) enum HeatingCoolingSystemResultKey {
    #[serde(rename = "Heating system")]
    HeatingSystem,
    #[serde(rename = "Heating system output")]
    HeatingSystemOutput,
    #[serde(rename = "Cooling system")]
    CoolingSystem,
    #[serde(rename = "Cooling system output")]
    CoolingSystemOutput,
}

#[derive(Clone, Copy, Deserialize_enum_str, Debug, Eq, Hash, PartialEq, Serialize_enum_str)]
pub(crate) enum HotWaterResultKey {
    #[serde(rename = "Hot water demand")]
    HotWaterDemand,
    #[serde(rename = "Hot water energy demand incl pipework_loss")]
    HotWaterEnergyDemandIncludingPipeworkLoss,
    #[serde(rename = "Hot water energy demand")]
    HotWaterEnergyDemand,
    #[serde(rename = "Hot water duration")]
    HotWaterDuration,
    #[serde(rename = "Hot Water Events")]
    HotWaterEvents,
    #[serde(rename = "Pipework losses")]
    PipeworkLosses,
    #[serde(rename = "Primary pipework losses")]
    PrimaryPipeworkLosses,
    #[serde(rename = "Storage losses")]
    StorageLosses,
}

fn get_cold_water_source_ref_for_type(
    source_type: ColdWaterSourceType,
    cold_water_sources: &ColdWaterSources,
) -> Option<ColdWaterSource> {
    cold_water_sources.get(&source_type).cloned()
}

pub type EventSchedule = Vec<Option<Vec<TypedScheduleEvent>>>;

fn event_schedules_from_input(
    events: &WaterHeatingEvents,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<EventSchedule> {
    let mut schedule: EventSchedule = vec![None; simulation_time_iterator.total_steps()];
    if let Some(&shower_events) = events.0.get(&WaterHeatingEventType::Shower).as_ref() {
        for (name, events) in shower_events {
            schedule = schedule_event_from_input(
                events.iter().collect(),
                name,
                WaterScheduleEventType::Shower,
                schedule,
                simulation_time_iterator,
            )?;
        }
    }

    if let Some(&bath_events) = events.0.get(&WaterHeatingEventType::Bath).as_ref() {
        for (name, events) in bath_events {
            schedule = schedule_event_from_input(
                events.iter().collect(),
                name,
                WaterScheduleEventType::Bath,
                schedule,
                simulation_time_iterator,
            )?;
        }
    }

    if let Some(&other_events) = events.0.get(&WaterHeatingEventType::Other).as_ref() {
        for (name, events) in other_events {
            schedule = schedule_event_from_input(
                events.iter().collect(),
                name,
                WaterScheduleEventType::Other,
                schedule,
                simulation_time_iterator,
            )?;
        }
    }

    Ok(schedule)
}

fn schedule_event_from_input(
    events_input: Vec<&WaterHeatingEvent>,
    name: &str,
    event_type: WaterScheduleEventType,
    existing_schedule: EventSchedule,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<EventSchedule> {
    let sim_timestep = simulation_time_iterator.step_in_hours();
    let total_timesteps = simulation_time_iterator.total_steps();

    expand_events(
        events_input
            .iter()
            .map(|event| ScheduleEvent::from(*event))
            .collect::<Vec<_>>(),
        sim_timestep,
        total_timesteps,
        name,
        event_type,
        existing_schedule,
    )
}

fn opening_area_total_from_zones(zones: &ZoneDictionary) -> f64 {
    zones
        .iter()
        .flat_map(|(_, zone)| {
            zone.building_elements
                .iter()
                .map(|(_, building_element)| match building_element {
                    BuildingElementInput::Transparent { height, width, .. } => height * width,
                    _ => 0.,
                })
        })
        .sum()
}

fn thermal_bridging_from_input(input: &ThermalBridgingInput) -> ThermalBridging {
    match input {
        ThermalBridgingInput::ThermalBridgingElements(input_bridges) => ThermalBridging::Bridges({
            let mut bridges = IndexMap::new();
            bridges.extend(input_bridges.iter().map(|(name, details)| {
                (
                    name.clone(),
                    match details {
                        ThermalBridgingDetails::Linear {
                            linear_thermal_transmittance,
                            length,
                            ..
                        } => ThermalBridge::Linear {
                            linear_thermal_transmittance: *linear_thermal_transmittance,
                            length: *length,
                        },
                        ThermalBridgingDetails::Point {
                            heat_transfer_coefficient,
                        } => ThermalBridge::Point {
                            heat_transfer_coefficient: *heat_transfer_coefficient,
                        },
                    },
                )
            }));
            bridges
        }),
        ThermalBridgingInput::ThermalBridgingNumber(num) => ThermalBridging::Number(*num),
    }
}

fn zone_from_input(
    input: &ZoneInput,
    external_conditions: Arc<ExternalConditions>,
    infiltration_ventilation: Arc<InfiltrationVentilation>,
    window_adjust_control: Option<Arc<Control>>,
    print_heat_balance: bool,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<(Zone, SystemReference, SystemReference)> {
    let heat_system_name = input.space_heat_system.clone();
    let cool_system_name = input.space_cool_system.clone();

    Ok((
        Zone::new(
            input.area,
            input.volume,
            input
                .building_elements
                .iter()
                .map(|(element_name, el)| {
                    Ok((
                        element_name.to_owned(),
                        building_element_from_input(el, external_conditions.clone())?,
                    ))
                })
                .collect::<anyhow::Result<IndexMap<String, BuildingElement>>>()?,
            thermal_bridging_from_input(&input.thermal_bridging),
            infiltration_ventilation,
            external_conditions.air_temp(&simulation_time_iterator.current_iteration()),
            input.temp_setpnt_init.unwrap(),
            window_adjust_control,
            print_heat_balance,
            simulation_time_iterator,
        )?,
        heat_system_name,
        cool_system_name,
    ))
}

#[allow(clippy::type_complexity)]
fn infiltration_ventilation_from_input(
    zones: &ZoneDictionary,
    input: &InfiltrationVentilationInput,
    controls: &Controls,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    simulation_time: &SimulationTimeIterator,
    total_volume: f64,
    external_conditions: Arc<ExternalConditions>,
    detailed_output_heating_cooling: bool,
) -> anyhow::Result<(
    InfiltrationVentilation,
    Option<Arc<Control>>,
    IndexMap<String, Arc<MechanicalVentilation>>,
    IndexMap<String, Vec<Ductwork>>,
)> {
    let windows: IndexMap<String, Window> = zones
        .values()
        .flat_map(|zone| {
            zone.building_elements
                .iter()
                .filter_map(|(building_element_name, building_element)| {
                    let window: Option<anyhow::Result<Window>> = if let BuildingElementInput::Transparent {
                        window_openable_control,
                        free_area_height,
                        mid_height,
                        max_window_open_area,
                        window_part_list,
                        orientation,
                        pitch,
                        ..
                    } = building_element
                    {
                        // Check control for openable window
                        let on_off_ctrl = window_openable_control.as_ref().and_then(|ctrl_name| {
                            controls.get_with_string(ctrl_name)
                        });

                        let window_result_fn = || {
                            Ok(Window::new(
                                external_conditions.clone(),
                                free_area_height.ok_or_else(|| anyhow!("A free_area_height value was expected for a transparent building element."))?,
                                mid_height.ok_or_else(|| anyhow!("A mid_height value was expected for a transparent building element."))?,
                                max_window_open_area.ok_or_else(|| anyhow!("A max_window_open_area value was expected for a transparent building element."))?,
                                window_part_list.as_ref().unwrap_or(&vec![]).clone(),
                                // running orientation through init_orientation again to convert it back to "orientation360" value, which is expected by the Python
                                // (orientation_difference function this is used within expects values in 0-360 degree range)
                                init_orientation(*orientation),
                                *pitch,
                                input.altitude,
                                on_off_ctrl
                            ))
                        };

                        Some(window_result_fn())
                    } else {
                        None
                    };
                    window.map(|window: anyhow::Result<Window>| anyhow::Ok((building_element_name.clone(), window?)))
                })
        })
        .collect::<anyhow::Result<IndexMap<String, Window>>>()?;

    let InfiltrationVentilationInput {
        cross_vent_factor: f_cross,
        shield_class,
        terrain_class,
        altitude,
        ..
    } = input;

    let window_adjust_control = input
        .window_adjust_control
        .as_ref()
        .and_then(|ctrl_name| controls.get_with_string(ctrl_name));

    let vents: IndexMap<String, Vent> = input
        .vents
        .iter()
        .map(|(vent_name, vent)| {
            (
                vent_name.clone(),
                Vent::new(
                    external_conditions.clone(),
                    vent.mid_height_air_flow_path,
                    vent.area_cm2,
                    vent.pressure_difference_ref,
                    // Python uses "orientation360" value here
                    init_orientation(vent.orientation),
                    vent.pitch,
                    *altitude,
                ),
            )
        })
        .collect();

    let (pitches, areas): (Vec<f64>, Vec<f64>) = zones
        .values()
        .flat_map(|zone| {
            zone.building_elements
                .values()
                .filter_map(|building_element| {
                    if let BuildingElementInput::Opaque { pitch, area, .. } = building_element {
                        (pitch_class(*pitch) == HeatFlowDirection::Upwards)
                            .then_some((*pitch, *area))
                    } else {
                        None
                    }
                })
        })
        .unzip();
    // Work out the average pitch, weighted by area
    let area_total = areas.iter().sum::<f64>();
    let average_pitch = if !pitches.is_empty() {
        areas
            .iter()
            .map(|x| x / area_total)
            .zip(pitches.iter())
            .map(|(x, &y)| x * y)
            .sum::<f64>()
    } else {
        // This case doesn't matter as if the area of roof = 0, the leakage coefficient = 0 anyway.
        0.
    };

    let (surface_area_facades_list, surface_area_roof_list) = zones
        .values()
        .flat_map(|zone| zone.building_elements.values())
        .fold((vec![], vec![]), |(mut facades, mut roofs), item| {
            match pitch_class(item.pitch()) {
                HeatFlowDirection::Horizontal => match item {
                    BuildingElementInput::Opaque { area, .. } => {
                        facades.push(*area);
                    }
                    BuildingElementInput::Transparent { height, width, .. } => {
                        facades.push(*height * *width);
                    }
                    _ => {}
                },
                HeatFlowDirection::Upwards => match item {
                    BuildingElementInput::Opaque { area, .. } => {
                        roofs.push(*area);
                    }
                    BuildingElementInput::Transparent { height, width, .. } => {
                        roofs.push(*height * *width)
                    }
                    _ => {}
                },
                _ => {}
            }

            (facades, roofs)
        });

    let surface_area_facades = surface_area_facades_list.iter().sum::<f64>();
    let surface_area_roof = surface_area_roof_list.iter().sum::<f64>();

    let leaks =
        CompletedVentilationLeaks::complete_input(input, surface_area_facades, surface_area_roof);

    // Empty map for air terminal devices until passive ducts work
    let atds: IndexMap<String, AirTerminalDevices> = Default::default();

    let mut mechanical_ventilations: IndexMap<String, Arc<MechanicalVentilation>> =
        Default::default();
    let mut space_heating_ductwork: IndexMap<String, Vec<Ductwork>> = Default::default();

    for (mech_vents_name, mech_vents_data) in input.mechanical_ventilation.iter() {
        let ctrl_intermittent_mev = mech_vents_data
            .control
            .as_ref()
            .and_then(|ctrl_name| controls.get_with_string(ctrl_name));

        let energy_supply = energy_supplies
            .get(&mech_vents_data.energy_supply)
            .ok_or_else(|| {
                anyhow!(
                    "The energy supply '{}' indicated for mechanical ventilation was not declared.",
                    mech_vents_data.energy_supply
                )
            })?;
        let energy_supply_connection =
            EnergySupply::connection(energy_supply.clone(), mech_vents_name)?;

        mechanical_ventilations.insert(
                mech_vents_name.clone(),
                Arc::new(MechanicalVentilation::new(external_conditions.clone(), mech_vents_data.supply_air_flow_rate_control, mech_vents_data.supply_air_temperature_control_type, 0., 0., mech_vents_data.vent_type, mech_vents_data.sfp.ok_or_else(|| anyhow!("A specific fan power value is expected for a mechanical ventilation unit."))?, mech_vents_data.design_outdoor_air_flow_rate, energy_supply_connection, total_volume, *altitude, ctrl_intermittent_mev, match mech_vents_data.vent_type {
                    VentType::Mvhr => mech_vents_data.mvhr_efficiency,
                    VentType::IntermittentMev
                    | VentType::CentralisedContinuousMev
                    | VentType::DecentralisedContinuousMev => {
                        None
                    }
                    VentType::Piv => bail!("PIV vent type is not currently recognised when building up mechanical ventilation values for calculation"),
                }, None)),
            );

        // TODO (from Python) not all dwellings have mech vents - update to make mech vents optional
        if mech_vents_data.vent_type == VentType::Mvhr {
            // the upstream Python incorrectly nixes the ductwork map here. it's fixed on more-recent-than-0.30 versions of upstream HEM,
            // so keep this in place for now and delete this comment when fixing 😉
            space_heating_ductwork = Default::default();
            space_heating_ductwork.insert(
                    mech_vents_name.to_owned(),
                    mech_vents_data
                        .ductwork
                        .as_ref()
                        .iter()
                        .flat_map(|ductworks| {
                            ductworks.iter().map(|ductwork| -> anyhow::Result<Ductwork> {
                                let (duct_perimeter, internal_diameter, external_diameter) =
                                    match ductwork.cross_section_shape {
                                        DuctShape::Circular => (None, Some(ductwork.internal_diameter_mm.ok_or_else(|| anyhow!("Expected an internal diameter value for ductwork with a circular cross-section."))? / MILLIMETRES_IN_METRE as f64), Some(ductwork.external_diameter_mm.ok_or_else(|| anyhow!("Expected an internal diameter value for ductwork with a circular cross-section."))? / MILLIMETRES_IN_METRE as f64)),
                                        DuctShape::Rectangular => (Some(ductwork.duct_perimeter_mm.ok_or_else(|| anyhow!("Expected a duct perimeter value for ductwork with a rectangular cross-section."))?), None, None),
                                    };

                                Ductwork::new(ductwork.cross_section_shape, duct_perimeter, internal_diameter, external_diameter, ductwork.length, ductwork.insulation_thermal_conductivity, ductwork.insulation_thickness_mm, ductwork.reflective, ductwork.duct_type, mech_vents_data.mvhr_location.ok_or_else(|| anyhow!("An MVHR location was expected for mechanical ventilation with an MVHR vent type."))?, mech_vents_data.mvhr_efficiency.ok_or_else(|| anyhow!("An MVHR efficiency value was expected for mechanical ventilation with an MVHR vent type."))?)
                            })
                        })
                        .collect::<anyhow::Result<Vec<Ductwork>>>()?,
                );
        }
    }

    let combustion_appliances: IndexMap<String, CombustionAppliances> = input
        .combustion_appliances
        .iter()
        .map(|(combustion_appliances_name, combustion_appliances_data)| {
            (
                combustion_appliances_name.to_owned(),
                CombustionAppliances::new(
                    combustion_appliances_data.supply_situation,
                    combustion_appliances_data.exhaust_situation,
                    combustion_appliances_data.fuel_type,
                    combustion_appliances_data.appliance_type,
                ),
            )
        })
        .collect();

    let ventilation = InfiltrationVentilation::new(
        external_conditions.clone(),
        *f_cross,
        *shield_class,
        *terrain_class,
        average_pitch,
        windows.into_values().collect(),
        vents.into_values().collect(),
        leaks,
        combustion_appliances.into_values().collect(),
        atds.into_values().collect(),
        mechanical_ventilations.values().cloned().collect(),
        detailed_output_heating_cooling,
        *altitude,
        zones.values().map(|zone| zone.area).sum::<f64>(),
    );

    Ok((
        ventilation,
        window_adjust_control,
        mechanical_ventilations,
        space_heating_ductwork,
    ))
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct CompletedVentilationLeaks {
    pub(crate) ventilation_zone_height: f64,
    pub(crate) test_pressure: f64,
    pub(crate) test_result: f64,
    pub(crate) area_roof: f64,
    pub(crate) area_facades: f64,
    pub(crate) env_area: f64,
    pub(crate) altitude: f64,
}

impl CompletedVentilationLeaks {
    fn complete_input(
        input: &InfiltrationVentilationInput,
        area_facades: f64,
        area_roof: f64,
    ) -> Self {
        let VentilationLeaks {
            ventilation_zone_height,
            test_pressure,
            test_result,
            env_area,
            ..
        } = input.leaks;
        Self {
            ventilation_zone_height,
            test_pressure,
            test_result,
            area_roof,
            area_facades,
            env_area,
            altitude: input.altitude,
        }
    }
}

fn building_element_from_input(
    input: &BuildingElementInput,
    external_conditions: Arc<ExternalConditions>,
) -> anyhow::Result<BuildingElement> {
    Ok(match input {
        BuildingElementInput::Opaque {
            is_unheated_pitched_roof,
            area,
            pitch,
            a_sol,
            r_c,
            k_m,
            mass_distribution_class,
            orientation,
            base_height,
            height,
            width,
            u_value,
            ..
        } => {
            let is_unheated_pitched_roof = if *pitch < PITCH_LIMIT_HORIZ_CEILING {
                is_unheated_pitched_roof
                    .ok_or_else(|| anyhow!("Pitch of opaque building element was {pitch} degrees, so it is necessary for this element to indicate whether this is an unheated pitched roof."))?
            } else {
                false
            };

            BuildingElement::Opaque(BuildingElementOpaque::new(
                *area,
                is_unheated_pitched_roof,
                *pitch,
                *a_sol,
                init_r_c_for_building_element(*r_c, *u_value, *pitch)?,
                *k_m,
                *mass_distribution_class,
                *orientation,
                *base_height,
                *height,
                *width,
                external_conditions,
            ))
        }
        BuildingElementInput::Transparent {
            u_value,
            r_c,
            pitch,
            orientation,
            g_value,
            frame_area_fraction,
            base_height,
            height,
            width,
            shading,
            ..
        } => BuildingElement::Transparent(BuildingElementTransparent::new(
            *pitch,
            init_r_c_for_building_element(*r_c, *u_value, *pitch)?,
            *orientation,
            *g_value,
            *frame_area_fraction,
            *base_height,
            *height,
            *width,
            shading.clone(),
            external_conditions,
        )),
        BuildingElementInput::Ground {
            area,
            total_area,
            pitch,
            u_value,
            r_f,
            k_m,
            mass_distribution_class,
            floor_type,
            height_upper_surface,
            thermal_transmission_walls,
            thermal_resistance_of_insulation,
            area_per_perimeter_vent,
            shield_fact_location,
            thickness_walls,
            depth_basement_floor,
            thermal_resistance_of_basement_walls,
            thermal_transmittance_of_floor_above_basement,
            height_basement_walls,
            perimeter,
            psi_wall_floor_junc,
            edge_insulation,
        } => {
            let (
                edge_insulation,
                height_upper_surface,
                thermal_transm_envi_base,
                thermal_transm_walls,
                area_per_perimeter_vent,
                shield_fact_location,
                thickness_walls,
                thermal_resist_insul,
                depth_basement_floor,
                thermal_resist_walls_base,
                height_basement_walls,
            ) = match floor_type {
                FloorType::SlabNoEdgeInsulation => (
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    *thickness_walls,
                    None,
                    None,
                    None,
                    None,
                ),
                FloorType::SlabEdgeInsulation => (
                    edge_insulation
                        .as_ref()
                        .map(|insulation| insulation.as_slice()),
                    None,
                    None,
                    None,
                    None,
                    None,
                    *thickness_walls,
                    None,
                    None,
                    None,
                    None,
                ),
                FloorType::SuspendedFloor => (
                    None,
                    *height_upper_surface,
                    None,
                    *thermal_transmission_walls,
                    *area_per_perimeter_vent,
                    *shield_fact_location,
                    *thickness_walls,
                    *thermal_resistance_of_insulation,
                    None,
                    None,
                    None,
                ),
                FloorType::HeatedBasement => (
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    *thickness_walls,
                    None,
                    *depth_basement_floor,
                    *thermal_resistance_of_basement_walls,
                    None,
                ),
                FloorType::UnheatedBasement => (
                    None,
                    None,
                    *thermal_transmittance_of_floor_above_basement,
                    *thermal_transmission_walls,
                    None,
                    None,
                    *thickness_walls,
                    None,
                    *depth_basement_floor,
                    *thermal_resistance_of_basement_walls,
                    *height_basement_walls,
                ),
            };

            BuildingElement::Ground(BuildingElementGround::new(
                *total_area,
                *area,
                *pitch,
                *u_value,
                *r_f,
                *k_m,
                *mass_distribution_class,
                *floor_type,
                edge_insulation,
                height_upper_surface,
                thermal_transm_envi_base,
                thermal_transm_walls,
                area_per_perimeter_vent,
                shield_fact_location,
                thickness_walls,
                thermal_resist_insul,
                depth_basement_floor,
                thermal_resist_walls_base,
                height_basement_walls,
                *perimeter,
                *psi_wall_floor_junc,
                external_conditions,
            )?)
        }
        BuildingElementInput::AdjacentZTC {
            area,
            pitch,
            u_value,
            r_c,
            k_m,
            mass_distribution_class,
        } => BuildingElement::AdjacentZTC(BuildingElementAdjacentZTC::new(
            *area,
            *pitch,
            init_r_c_for_building_element(*r_c, *u_value, *pitch)?,
            *k_m,
            *mass_distribution_class,
            external_conditions,
        )),
        BuildingElementInput::AdjacentZTUSimple {
            area,
            pitch,
            u_value,
            r_c,
            r_u,
            k_m,
            mass_distribution_class,
        } => BuildingElement::AdjacentZTUSimple(BuildingElementAdjacentZTUSimple::new(
            *area,
            *pitch,
            init_r_c_for_building_element(*r_c, *u_value, *pitch)?,
            *r_u,
            *k_m,
            *mass_distribution_class,
            external_conditions,
        )),
    })
}

fn init_r_c_for_building_element(
    r_c: Option<f64>,
    u_value: Option<f64>,
    pitch: f64,
) -> anyhow::Result<f64> {
    Ok(if let Some(r_c) = r_c {
        r_c
    } else {
        convert_uvalue_to_resistance(
            u_value.ok_or_else(|| {
                anyhow!(
                    "Neither r_c nor u_value were provided for one of the building element inputs."
                )
            })?,
            pitch,
        )
    })
}

fn set_up_energy_supply_unmet_demand_zones(
    unmet_demand_supply: Arc<RwLock<EnergySupply>>,
    zones: &ZoneDictionary,
) -> IndexMap<String, Arc<EnergySupplyConnection>> {
    let mut energy_supplies: IndexMap<String, Arc<EnergySupplyConnection>> = Default::default();

    for name in zones.keys() {
        energy_supplies.insert(
            (*name).clone(),
            Arc::new(EnergySupply::connection(unmet_demand_supply.clone(), name.as_str()).unwrap()),
        );
    }

    energy_supplies
}

fn apply_appliance_gains_from_input(
    internal_gains_collection: &mut InternalGainsCollection,
    input: &ApplianceGainsInput,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    total_floor_area: f64,
    simulation_time: &SimulationTimeIterator,
) -> anyhow::Result<()> {
    for (name, gains_details) in input {
        let energy_supply_conn = EnergySupply::connection(
            energy_supplies
                .get(&gains_details.energy_supply)
                .ok_or_else(|| {
                    anyhow!(
                        "Appliance gains reference an undeclared energy supply '{}'.",
                        gains_details.energy_supply
                    )
                })?
                .clone(),
            name.as_str(),
        )?;

        let gains = if gains_details.events.is_some() && gains_details.standby.is_some() {
            Gains::Event(EventApplianceGains::new(
                energy_supply_conn,
                simulation_time,
                gains_details,
                total_floor_area,
            ))
        } else {
            Gains::Appliance(appliance_gains_from_single_input(
                gains_details,
                energy_supply_conn,
                total_floor_area,
            )?)
        };

        internal_gains_collection.insert(name.clone(), gains);
    }

    Ok(())
}

fn appliance_gains_from_single_input(
    input: &ApplianceGainsDetails,
    energy_supply_connection: EnergySupplyConnection,
    total_floor_area: f64,
) -> anyhow::Result<ApplianceGains> {
    let total_energy_supply = reject_nulls(expand_numeric_schedule(
        input
            .schedule
            .as_ref()
            .ok_or_else(|| anyhow!("Appliance gains did not have schedule when expected."))?,
    ))?
    .iter()
    .map(|energy_data| energy_data / total_floor_area)
    .collect();

    Ok(ApplianceGains::new(
        total_energy_supply,
        input.gains_fraction,
        input.start_day,
        input.time_series_step,
        energy_supply_connection,
    ))
}

#[derive(Debug)]
pub enum HeatSource {
    Storage(HeatSourceWithStorageTank),
    Wet(Box<HeatSourceWet>),
}

impl HeatSource {
    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        temp_return: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        match self {
            HeatSource::Storage(ref mut storage) => match storage {
                HeatSourceWithStorageTank::Immersion(imm) => Ok(imm
                    .lock()
                    .demand_energy(energy_demand, simulation_time_iteration)),
                HeatSourceWithStorageTank::Solar(ref solar) => Ok(solar
                    .lock()
                    .demand_energy(energy_demand, simulation_time_iteration.index)),
            },
            HeatSource::Wet(ref mut wet) => {
                wet.demand_energy(energy_demand, temp_return, simulation_time_iteration)
            }
        }
    }

    pub fn as_immersion_heater(&self) -> Option<Arc<Mutex<ImmersionHeater>>> {
        match self {
            HeatSource::Storage(HeatSourceWithStorageTank::Immersion(immersion)) => {
                Some(immersion.clone())
            }
            _ => None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct PositionedHeatSource {
    pub heat_source: Arc<Mutex<HeatSource>>,
    pub heater_position: f64,
    pub thermostat_position: f64,
}

#[derive(Clone, Debug)]
pub enum WetHeatSource {
    HeatPump(Arc<Mutex<HeatPump>>),
    Boiler(Arc<RwLock<Boiler>>),
    Hiu(Arc<Mutex<HeatNetwork>>),
    HeatBattery(Arc<Mutex<HeatBattery>>),
}

impl WetHeatSource {
    pub fn timestep_end(&self, simtime: SimulationTimeIteration) {
        match self {
            WetHeatSource::HeatPump(heat_pump) => heat_pump.lock().timestep_end(simtime.index),
            WetHeatSource::Boiler(boiler) => boiler.write().timestep_end(simtime),
            WetHeatSource::Hiu(heat_network) => heat_network.lock().timestep_end(simtime.index),
            WetHeatSource::HeatBattery(heat_battery) => {
                heat_battery.lock().timestep_end(simtime.index)
            }
        }
    }

    pub(crate) fn create_service_hot_water_combi(
        &mut self,
        boiler_data: HotWaterSourceDetails,
        service_name: &str,
        temp_hot_water: f64,
        cold_feed: WaterSourceWithTemperature,
    ) -> anyhow::Result<BoilerServiceWaterCombi> {
        match self {
            WetHeatSource::HeatPump(heat_pump) => heat_pump.lock().create_service_hot_water_combi(
                boiler_data,
                service_name,
                temp_hot_water,
                cold_feed,
            ),
            WetHeatSource::Boiler(ref mut boiler) => Boiler::create_service_hot_water_combi(
                boiler.clone(),
                boiler_data,
                service_name.to_string(),
                temp_hot_water,
                cold_feed,
            )
            .map_err(|err| anyhow!(format!("{err}"))),
            _ => {
                bail!("Expect to only be able to create a hot water combi service for boilers and heat pumps.")
            }
        }
    }

    fn output_detailed_results(
        &self,
        hot_water_energy_output: &[f64],
    ) -> Option<(ResultsPerTimestep, ResultsAnnual)> {
        if let WetHeatSource::HeatPump(heat_pump) = self {
            Some(
                heat_pump.lock().clone().output_detailed_results(
                    hot_water_energy_output
                        .iter()
                        .map(|&x| x.into())
                        .collect_vec(),
                ),
            )
        } else {
            None
        }
    }
}

fn heat_source_wet_from_input(
    energy_supply_name: &str,
    input: HeatSourceWetDetails,
    external_conditions: Arc<ExternalConditions>,
    simulation_time: Arc<SimulationTimeIterator>,
    mechanical_ventilations: &IndexMap<String, Arc<MechanicalVentilation>>,
    number_of_zones: usize,
    temp_internal_air_accessor: TempInternalAirAccessor,
    controls: &Controls,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    detailed_output_heating_cooling: bool,
) -> anyhow::Result<WetHeatSource> {
    match &input {
        HeatSourceWetDetails::HeatPump {
            source_type,
            energy_supply_heat_network,
            energy_supply,
            boiler,
            ..
        } => {
            let throughput_exhaust_air = if source_type.is_exhaust_air() {
                // Check that ventilation system is compatible with exhaust air HP
                // the Python code here will assign to throughput_exhaust_air on the last iteration, though mech vents collection is not ordered - this may be erroneous
                // reported to BRE as possible bug https://dev.azure.com/BreGroup/SAP%2011/_workitems/edit/45508
                let mut throughput_exhaust_air: Option<f64> = Default::default();
                for mech_vent in mechanical_ventilations.values() {
                    match mech_vent.vent_type() {
                        VentType::IntermittentMev | VentType::DecentralisedContinuousMev => {
                            bail!("Exhaust air heat pump does not work with Intermittent MEV or Decentralised continuous MEV.")
                        }
                        _ => {
                            throughput_exhaust_air =
                                Some(mech_vent.as_ref().design_outdoor_air_flow_rate_m3_h);
                            // In Python code data is appended here to a project instance variable called __heat_source_wet_names_requiring_overvent
                            // but this is never read, so we are leaving out the implementation here
                        }
                    }
                }
                throughput_exhaust_air
            } else {
                None
            };

            let energy_supply_hn = if matches!(source_type, HeatPumpSourceType::HeatNetwork) {
                let energy_supply_heat_network = energy_supply_heat_network.as_ref().ok_or_else(|| anyhow!("A heat pump with a heat network source is expected to reference an energy supply for the heat network."))?;
                Some(energy_supplies.get(energy_supply_heat_network).ok_or_else(|| anyhow!("A heat network with a heat network source references an undeclared energy supply '{energy_supply_heat_network}'."))?.clone())
            } else {
                None
            };

            let (boiler, cost_schedule_hybrid_hp) = if let Some(boiler) = boiler {
                let energy_supply_boiler = energy_supplies
                    .get(&boiler.energy_supply)
                    .ok_or_else(|| {
                        anyhow!(
                            "A boiler references an undeclared energy supply '{}'.",
                            boiler.energy_supply
                        )
                    })?
                    .clone();
                let energy_supply_aux_boiler = energy_supplies
                    .get(&boiler.energy_supply_auxiliary)
                    .ok_or_else(|| {
                        anyhow!(
                            "A boiler references an undeclared energy supply '{}'.",
                            boiler.energy_supply_auxiliary
                        )
                    })?
                    .clone();
                let energy_supply_conn_aux_boiler = EnergySupply::connection(
                    energy_supply_aux_boiler,
                    format!("Boiler_auxiliary: {energy_supply_name}").as_str(),
                )?;

                let cost_schedule_hybrid_hp = boiler.cost_schedule_hybrid.clone();

                let boiler = Boiler::new(
                    boiler.as_ref().into(),
                    energy_supply_boiler,
                    energy_supply_conn_aux_boiler,
                    external_conditions.clone(),
                    simulation_time.step_in_hours(),
                )?;

                (Some(boiler), cost_schedule_hybrid_hp)
            } else {
                Default::default()
            };

            let energy_supply = energy_supplies
                .get(energy_supply)
                .ok_or_else(|| {
                    anyhow!("Heat pump references an undeclared energy supply '{energy_supply}'.")
                })?
                .clone();
            let energy_supply_conn_name_auxiliary =
                format!("HeatPump_auxiliary: {energy_supply_name}");

            Ok(WetHeatSource::HeatPump(Arc::new(Mutex::new(
                HeatPump::new(
                    &input,
                    energy_supply,
                    &energy_supply_conn_name_auxiliary,
                    simulation_time.step_in_hours(),
                    external_conditions.clone(),
                    number_of_zones,
                    throughput_exhaust_air,
                    energy_supply_hn,
                    detailed_output_heating_cooling,
                    boiler.map(|boiler: Boiler| Arc::new(RwLock::new(boiler))),
                    cost_schedule_hybrid_hp,
                    temp_internal_air_fn(temp_internal_air_accessor),
                )?,
            ))))
        }
        HeatSourceWetDetails::Boiler {
            energy_supply,
            energy_supply_auxiliary,
            ..
        } => {
            let energy_supply = energy_supplies
                .get(energy_supply)
                .ok_or_else(|| {
                    anyhow!("Boiler references undeclared energy supply '{energy_supply}'.")
                })?
                .clone();
            let energy_supply_aux = energy_supplies.get(energy_supply_auxiliary).ok_or_else(|| anyhow!("Boiler references undeclared auxiliary energy supply '{energy_supply_auxiliary}'."))?.clone();
            let aux_supply_name = format!("Boiler_auxiliary: {energy_supply_name}");
            let energy_supply_conn_aux =
                EnergySupply::connection(energy_supply_aux.clone(), aux_supply_name.as_str())?;

            Ok(WetHeatSource::Boiler(Arc::new(RwLock::new(
                Boiler::new(
                    input,
                    energy_supply,
                    energy_supply_conn_aux,
                    external_conditions.clone(),
                    simulation_time.step_in_hours(),
                )
                .expect("could not construct boiler value from provided data"),
            ))))
        }
        HeatSourceWetDetails::Hiu {
            power_max,
            hiu_daily_loss,
            building_level_distribution_losses,
            energy_supply,
            ..
        } => {
            let energy_supply = energy_supplies
                .get(energy_supply)
                .ok_or_else(|| {
                    anyhow!("HIU references an undeclared energy supply '{energy_supply}'.")
                })?
                .clone();
            let energy_supply_conn_name_auxiliary =
                format!("HeatNetwork_auxiliary: {energy_supply_name}");
            let energy_supply_conn_name_building_level_distribution_losses =
                format!("HeatNetwork_building_level_distribution_losses: {energy_supply_name}");

            Ok(WetHeatSource::Hiu(Arc::new(Mutex::new(HeatNetwork::new(
                *power_max,
                *hiu_daily_loss,
                *building_level_distribution_losses,
                energy_supply,
                energy_supply_conn_name_auxiliary,
                energy_supply_conn_name_building_level_distribution_losses,
                simulation_time.step_in_hours(),
            )))))
        }
        HeatSourceWetDetails::HeatBattery {
            control_charge,
            energy_supply,
            ..
        } => {
            let energy_supply = energy_supplies
                .get(energy_supply)
                .ok_or_else(|| {
                    anyhow!(
                        "Heat battery references an undeclared energy supply '{energy_supply}'."
                    )
                })?
                .clone();
            let energy_supply_conn =
                EnergySupply::connection(energy_supply.clone(), energy_supply_name).unwrap();

            let heat_source = WetHeatSource::HeatBattery(Arc::new(Mutex::new(HeatBattery::new(
                &input,
                controls
                    .get_with_string(control_charge)
                    .unwrap_or_else(|| {
                        panic!(
                            "expected a control to be registered with the name '{control_charge}'"
                        )
                    })
                    .clone(),
                energy_supply,
                energy_supply_conn,
                simulation_time,
            ))));
            Ok(heat_source)
        }
    }
}

fn heat_source_from_input(
    name: &str,
    input: &HeatSourceInput,
    cold_water_source: &WaterSourceWithTemperature,
    temp_setpoint: f64,
    volume: f64,
    daily_losses: f64,
    heat_exchanger_surface_area: Option<f64>,
    wet_heat_sources: &IndexMap<String, WetHeatSource>,
    simulation_time: &SimulationTimeIterator,
    controls: &Controls,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    temp_internal_air_accessor: TempInternalAirAccessor,
    external_conditions: Arc<ExternalConditions>,
) -> anyhow::Result<(HeatSource, String)> {
    match input {
        HeatSourceInput::ImmersionHeater {
            power,
            control,
            energy_supply,
            ..
        } => {
            let energy_supply = energy_supplies.get(energy_supply).ok_or_else(|| anyhow!("Immersion heater references an undeclared energy supply '{energy_supply}'."))?.clone();
            let energy_supply_conn = EnergySupply::connection(energy_supply.clone(), name)?;

            Ok((
                HeatSource::Storage(HeatSourceWithStorageTank::Immersion(Arc::new(Mutex::new(
                    ImmersionHeater::new(
                        *power,
                        energy_supply_conn,
                        simulation_time.step_in_hours(),
                        (*control).and_then(|ctrl| controls.get_with_string(&ctrl.to_string())),
                    ),
                )))),
                name.into(),
            ))
        }
        HeatSourceInput::SolarThermalSystem {
            solar_cell_location,
            area_module,
            modules,
            peak_collector_efficiency,
            incidence_angle_modifier,
            first_order_hlc,
            second_order_hlc,
            collector_mass_flow_rate,
            power_pump,
            power_pump_control,
            tilt,
            orientation,
            solar_loop_piping_hlc,
            energy_supply,
            ..
        } => {
            let energy_supply = energy_supplies.get(energy_supply).ok_or_else(|| anyhow!("Solar thermal system references an undeclared energy supply '{energy_supply}'."))?.clone();
            let energy_supply_conn = EnergySupply::connection(energy_supply.clone(), name)?;

            Ok((
                HeatSource::Storage(HeatSourceWithStorageTank::Solar(Arc::new(Mutex::new(
                    SolarThermalSystem::new(
                        *solar_cell_location,
                        *area_module,
                        *modules,
                        *peak_collector_efficiency,
                        *incidence_angle_modifier,
                        *first_order_hlc,
                        *second_order_hlc,
                        *collector_mass_flow_rate,
                        *power_pump,
                        *power_pump_control,
                        energy_supply_conn,
                        *tilt,
                        *orientation,
                        *solar_loop_piping_hlc,
                        external_conditions.clone(),
                        temp_internal_air_fn(temp_internal_air_accessor),
                        simulation_time.step_in_hours(),
                        *WATER,
                    ),
                )))),
                name.into(),
            ))
        }
        HeatSourceInput::Wet {
            name,
            control,
            temp_flow_limit_upper,
            ..
        } => {
            let energy_supply_conn_name = format!("{name}_water_heating");
            let heat_source_wet = wet_heat_sources
                .get(name)
                .unwrap_or_else(|| {
                    panic!("Expected a wet heat source registered with the name '{name}'.")
                })
                .clone();
            let source_control = (*control).and_then(|ctrl| controls.get(&ctrl));

            let mut heat_source_wet_clone = heat_source_wet.clone();

            Ok((
                match heat_source_wet_clone {
                    WetHeatSource::HeatPump(heat_pump) => HeatSource::Wet(Box::new(
                        HeatSourceWet::HeatPumpWater(HeatPump::create_service_hot_water(
                            heat_pump.clone(),
                            energy_supply_conn_name.clone(),
                            55.,
                            temp_flow_limit_upper
                                .expect("temp_flow_limit_upper field was expected to be set"),
                            cold_water_source.as_cold_water_source()?,
                            source_control,
                        )),
                    )),
                    WetHeatSource::Boiler(ref mut boiler) => HeatSource::Wet(Box::new(
                        HeatSourceWet::WaterRegular(Boiler::create_service_hot_water_regular(
                            boiler.clone(),
                            energy_supply_conn_name.clone(),
                            temp_setpoint,
                            source_control,
                        )),
                    )),
                    WetHeatSource::Hiu(heat_network) => {
                        HeatSource::Wet(Box::new(HeatSourceWet::HeatNetworkWaterStorage(
                            HeatNetwork::create_service_hot_water_storage(
                                heat_network,
                                energy_supply_conn_name.clone(),
                                temp_setpoint,
                                source_control,
                            ),
                        )))
                    }
                    WetHeatSource::HeatBattery(battery) => {
                        HeatSource::Wet(Box::new(HeatSourceWet::HeatBatteryHotWater(
                            HeatBattery::create_service_hot_water_regular(
                                battery,
                                &energy_supply_conn_name,
                                temp_setpoint,
                                source_control,
                            ),
                        )))
                    }
                },
                energy_supply_conn_name,
            ))
        }
        HeatSourceInput::HeatPumpHotWaterOnly {
            power_max,
            vol_hw_daily_average,
            tank_volume_declared,
            heat_exchanger_surface_area_declared,
            daily_losses_declared,
            ref test_data,
            energy_supply,
            control,
            in_use_factor_mismatch,
            ..
        } => {
            let energy_supply = energy_supplies.get(energy_supply).ok_or_else(|| anyhow!("Heat pump (hot water only) references an undeclared energy supply '{energy_supply}'."))?.clone();
            let energy_supply_conn_name = name;
            let energy_supply_connection =
                EnergySupply::connection(energy_supply.clone(), energy_supply_conn_name)?;

            Ok((
                HeatSource::Wet(Box::new(HeatSourceWet::HeatPumpWaterOnly(
                    HeatPumpHotWaterOnly::new(
                        *power_max,
                        energy_supply_connection,
                        test_data,
                        *vol_hw_daily_average,
                        volume,
                        daily_losses,
                        heat_exchanger_surface_area.ok_or(anyhow!("A heat exchanger surface_area is expected for a HeatPumpHotWaterOnly heat source"))?,
                        *in_use_factor_mismatch,
                        *tank_volume_declared,
                        *heat_exchanger_surface_area_declared,
                        *daily_losses_declared,
                        simulation_time.step_in_hours(),
                        // TODO fix control logic as part of 0.32 migration
                        controls.get(&control.expect("This may not be present during migration to 0.32")),
                    ),
                ))),
                energy_supply_conn_name.into(),
            ))
        }
    }
}

#[derive(Debug)]
pub enum HotWaterSource {
    StorageTank(Arc<Mutex<StorageTank>>),
    CombiBoiler(BoilerServiceWaterCombi),
    PointOfUse(PointOfUse),
    HeatNetwork(HeatNetworkServiceWaterDirect),
}

impl HotWaterSource {
    pub fn get_cold_water_source(&self) -> WaterSourceWithTemperature {
        match self {
            HotWaterSource::StorageTank(source) => source.lock().get_cold_water_source().clone(),
            HotWaterSource::CombiBoiler(source) => source.get_cold_water_source().clone(),
            HotWaterSource::PointOfUse(source) => source.get_cold_water_source().clone(),
            HotWaterSource::HeatNetwork(source) => source.get_cold_water_source().clone(),
        }
    }

    pub(crate) fn temp_hot_water(&self) -> f64 {
        match self {
            HotWaterSource::StorageTank(storage_tank) => storage_tank.lock().get_temp_hot_water(),
            HotWaterSource::CombiBoiler(combi) => combi.temperature_hot_water_in_c(),
            HotWaterSource::PointOfUse(point_of_use) => point_of_use.get_temp_hot_water(),
            HotWaterSource::HeatNetwork(heat_network) => heat_network.temp_hot_water(),
        }
    }

    pub fn demand_hot_water(
        &self,
        vol_demand_target: IndexMap<DemandVolTargetKey, VolumeReference>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        match self {
            HotWaterSource::StorageTank(_) => {
                // StorageTank does not match the same method signature or return type as all other Hot Water sources
                panic!("demand_hot_water for HotWaterSource::StorageTank should be called directly on the HotWaterSource::StorageTank");
            }
            HotWaterSource::CombiBoiler(ref source) => source
                .demand_hot_water(vol_demand_target, simulation_time_iteration)
                .expect("Combi boiler could not calc demand hot water."),
            HotWaterSource::PointOfUse(ref source) => {
                source.demand_hot_water(vol_demand_target, &simulation_time_iteration)
            }
            HotWaterSource::HeatNetwork(ref source) => {
                source.demand_hot_water(vol_demand_target, simulation_time_iteration)
            }
        }
    }
}

fn hot_water_source_from_input(
    source_name: String,
    input: &HotWaterSourceDetails,
    cold_water_sources: &ColdWaterSources,
    wet_heat_sources: &mut IndexMap<String, WetHeatSource>,
    wwhrs: &IndexMap<String, Arc<Mutex<Wwhrs>>>,
    controls: &Controls,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    diverter_types: &DiverterTypes,
    diverters: &mut Vec<Arc<RwLock<PVDiverter>>>,
    temp_internal_air_accessor: TempInternalAirAccessor,
    simulation_time: &SimulationTimeIterator,
    external_conditions: Arc<ExternalConditions>,
) -> anyhow::Result<(HotWaterSource, Vec<String>)> {
    let mut energy_supply_conn_names = vec![];
    let cloned_input = input.clone();
    let hot_water_source = match input {
        HotWaterSourceDetails::StorageTank {
            volume,
            daily_losses,
            heat_exchanger_surface_area,
            min_temp,
            setpoint_temp,
            control_hold_at_setpoint,
            cold_water_source: cold_water_source_type,
            primary_pipework,
            heat_source,
            ..
        } => {
            let mut cold_water_source: WaterSourceWithTemperature = cold_water_source_for_type(
                match cold_water_source_type {
                    ColdWaterSourceReference::Type(source_type) => source_type,
                    ColdWaterSourceReference::Reference(_) => {
                        unimplemented!("TODO 0.32 allow for free references to cold water sources")
                    }
                },
                cold_water_sources,
            )?;
            // TODO (from Python) Need to handle error if ColdWaterSource name is invalid.
            // TODO (from Python) assuming here there is only one WWHRS
            if !wwhrs.is_empty() {
                for heat_recovery_system in wwhrs.values() {
                    match *heat_recovery_system.lock() {
                        Wwhrs::WWHRSInstantaneousSystemC(_)
                        | Wwhrs::WWHRSInstantaneousSystemA(_) => {
                            cold_water_source =
                                WaterSourceWithTemperature::Wwhrs(heat_recovery_system.clone());
                        }
                        _ => {}
                    }
                }
            }
            // At this point in the Python, the internal_diameter and external_diameter fields on
            // primary_pipework are updated, this is done in Pipework.rs in the Rust
            let primary_pipework_lst = primary_pipework.as_ref();
            let mut heat_sources: IndexMap<String, PositionedHeatSource> = Default::default();

            for (name, hs) in heat_source {
                let heater_position = hs.heater_position();
                let thermostat_position = hs.thermostat_position();

                //heat exchanger area
                let heat_exchanger_surface_area =
                    if let HeatSourceInput::HeatPumpHotWaterOnly { .. } = hs {
                        *heat_exchanger_surface_area
                    } else {
                        None
                    };

                let (heat_source, energy_supply_conn_name) = heat_source_from_input(
                    name.as_str(),
                    hs,
                    &cold_water_source,
                    setpoint_temp.ok_or_else(|| anyhow!("A setpoint temp for a storage tank was expected to be available when building the corpus for the HEM calculation."))?,
                    *volume,
                    *daily_losses,
                    heat_exchanger_surface_area,
                    wet_heat_sources,
                    simulation_time,
                    controls,
                    energy_supplies,
                    temp_internal_air_accessor.clone(),
                    external_conditions.clone(),
                )?;
                let heat_source = Arc::new(Mutex::new(heat_source));
                heat_sources.insert(
                    name.clone(),
                    PositionedHeatSource {
                        heat_source: heat_source.clone(),
                        heater_position,
                        thermostat_position,
                    },
                );
                energy_supply_conn_names.push(energy_supply_conn_name);
            }
            let ctrl_hold_at_setpoint = control_hold_at_setpoint
                .as_ref()
                .and_then(|ctrl| controls.get_with_string(ctrl.as_str()));

            let storage_tank = Arc::new(Mutex::new(StorageTank::new(
                *volume,
                *daily_losses,
                min_temp.ok_or_else(|| anyhow!("A min temp for a storage tank was expected to be available when building the corpus for the HEM calculation."))?,
                setpoint_temp.ok_or_else(|| anyhow!("A setpoint temp for a storage tank was expected to be available when building the corpus for the HEM calculation."))?,
                cold_water_source,
                simulation_time.step_in_hours(),
                heat_sources.clone(),
                temp_internal_air_fn(temp_internal_air_accessor),
                external_conditions,
                Some(24),
                primary_pipework_lst,
                Some(EnergySupply::connection(
                    energy_supplies.get(UNMET_DEMAND_SUPPLY_NAME).expect("Energy supply representing unmet demand was expected to have been declared.").clone(),
                    &source_name,
                )?),
                ctrl_hold_at_setpoint,
                *WATER,
            )));
            for (heat_source_name, hs) in heat_source {
                let energy_supply_type = hs.energy_supply_type();
                if let Some(diverter) = diverter_types.get_for_supply_type(energy_supply_type) {
                    if diverter.storage_tank.matches(&source_name)
                        && diverter.heat_source.matches(heat_source_name)
                    {
                        let energy_supply = energy_supplies.ensured_get_for_type(
                            hs.energy_supply_type(),
                            simulation_time.total_steps(),
                        )?;

                        let positioned_heat_source =
                            &heat_sources.get(&heat_source_name.to_owned());
                        let immersion_heater = positioned_heat_source
                            .unwrap()
                            .heat_source
                            .lock()
                            .as_immersion_heater();

                        if let Some(im) = immersion_heater {
                            let pv_diverter =
                                PVDiverter::new(storage_tank.clone(), im, heat_source_name.clone());
                            energy_supply
                                .write()
                                .connect_diverter(pv_diverter.clone())
                                .unwrap();
                            diverters.push(pv_diverter);
                        }
                    }
                }
            }
            HotWaterSource::StorageTank(storage_tank)
        }
        HotWaterSourceDetails::CombiBoiler {
            cold_water_source: cold_water_source_type,
            heat_source_wet: heat_source_wet_type,
            setpoint_temp,
            ..
        } => {
            let cold_water_source =
                cold_water_source_for_type(cold_water_source_type, cold_water_sources)?;
            let energy_supply_conn_name = format!(
                "{}_water_heating",
                heat_source_wet_type.to_canonical_string()
            );
            energy_supply_conn_names.push(energy_supply_conn_name.clone());
            let heat_source_wet = match heat_source_wet_type {
                HeatSourceWetType::Boiler => wet_heat_sources
                    .get_mut(&heat_source_wet_type.to_canonical_string())
                    .ok_or_else(|| {
                        anyhow!("Expected a boiler to have been defined as a wet heat source")
                    })?,
                HeatSourceWetType::HeatPump => wet_heat_sources
                    .get_mut(&heat_source_wet_type.to_canonical_string())
                    .ok_or_else(|| {
                        anyhow!("Expected a heat pump to have been defined as a wet heat source")
                    })?,
                _ => {
                    panic!("Did not expect a heat source type that was not a boiler or a heat pump")
                }
            };
            HotWaterSource::CombiBoiler(
                heat_source_wet
                    .create_service_hot_water_combi(
                        cloned_input,
                        energy_supply_conn_name.as_str(),
                        setpoint_temp.ok_or_else(|| {
                            anyhow!("A setpoint temp was expected on a combi boiler input.")
                        })?,
                        cold_water_source,
                    )
                    .expect("expected to be able to instantiate a combi boiler object"),
            )
        }
        HotWaterSourceDetails::PointOfUse {
            efficiency,
            cold_water_source: cold_water_source_type,
            energy_supply,
            setpoint_temp,
            ..
        } => {
            let energy_supply = energy_supplies
                .ensured_get_for_type(*energy_supply, simulation_time.total_steps())?;
            let energy_supply_conn_name = source_name;
            energy_supply_conn_names.push(energy_supply_conn_name.clone());
            let energy_supply_conn =
                EnergySupply::connection(energy_supply.clone(), &energy_supply_conn_name)?;
            let cold_water_source =
                cold_water_source_for_type(cold_water_source_type, cold_water_sources)?;
            HotWaterSource::PointOfUse(PointOfUse::new(
                *efficiency,
                energy_supply_conn,
                cold_water_source,
                (*setpoint_temp).ok_or_else(|| {
                    anyhow!(
                        "A setpoint_temp value was expected on a point of use hot water source."
                    )
                })?,
            ))
        }
        HotWaterSourceDetails::Hiu {
            cold_water_source: cold_water_source_type,
            heat_source_wet: heat_source_wet_type,
            setpoint_temp,
            ..
        } => {
            let energy_supply_conn_name = format!(
                "{}_water_heating",
                heat_source_wet_type.to_canonical_string()
            );
            energy_supply_conn_names.push(energy_supply_conn_name.clone());
            let cold_water_source =
                cold_water_source_for_type(cold_water_source_type, cold_water_sources)?;
            let heat_source_wet = match heat_source_wet_type {
                HeatSourceWetType::HeatNetwork => {
                    match wet_heat_sources
                        .get("HeatNetwork")
                        .expect("expected a heat network in this context")
                    {
                        WetHeatSource::Hiu(heat_network) => heat_network.clone(),
                        _ => panic!("expected a heat network in this context"),
                    }
                }
                _ => panic!("expected a heat network in this context"),
            };
            HotWaterSource::HeatNetwork(HeatNetwork::create_service_hot_water_direct(
                heat_source_wet.clone(),
                energy_supply_conn_name,
                (*setpoint_temp).ok_or_else(|| {
                    anyhow!(
                        "A setpoint_temp value was expected on a point of use hot water source."
                    )
                })?,
                cold_water_source,
            ))
        }
    };

    Ok((hot_water_source, energy_supply_conn_names))
}

fn cold_water_source_for_type(
    cold_water_source_type: &ColdWaterSourceType,
    cold_water_sources: &ColdWaterSources,
) -> anyhow::Result<WaterSourceWithTemperature> {
    Ok(WaterSourceWithTemperature::ColdWaterSource(Arc::new(
        cold_water_sources
            .get(cold_water_source_type)
            .ok_or_else(|| anyhow!("referenced cold water source was expected to exist"))?
            .clone(),
    )))
}

fn space_heat_systems_from_input(
    input: &SpaceHeatSystemInput,
    controls: &Controls,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    simulation_time: &SimulationTimeIterator,
    heat_sources_wet: &IndexMap<String, WetHeatSource>,
    heat_system_names_requiring_overvent: &mut Vec<String>,
    heat_system_name_for_zone: &IndexMap<String, String>,
    zones: &Arc<IndexMap<String, Zone>>,
    heat_sources_wet_with_buffer_tank: &[String],
    temp_internal_air_accessor: TempInternalAirAccessor,
    external_conditions: Arc<ExternalConditions>,
    detailed_output_heating_cooling: bool,
) -> anyhow::Result<SpaceHeatSystemsWithEnergyConnections> {
    let mut energy_conn_names_for_systems: IndexMap<String, String> = Default::default();
    let space_heat_systems = input
        .iter()
        .filter(|(system_name, _)| heat_system_name_for_zone.values().any(|heat_system_name| heat_system_name == system_name.as_str()))
        .map(|(system_name, space_heat_system_details)| {
            Ok((
                (*system_name).clone(),
                Arc::new(Mutex::new(match space_heat_system_details {
                    SpaceHeatSystemDetails::InstantElectricHeater {
                        rated_power,
                        control,
                        frac_convective,
                        energy_supply,
                        ..
                    } => {
                        let energy_supply = energy_supplies.ensured_get_for_type(*energy_supply, simulation_time.total_steps())?;
                        let energy_supply_conn_name = system_name;
                        energy_conn_names_for_systems.insert(system_name.clone(), energy_supply_conn_name.clone());
                        let energy_supply_conn = EnergySupply::connection(energy_supply, energy_supply_conn_name.as_str()).unwrap();
                        SpaceHeatSystem::Instant(InstantElecHeater::new(
                            *rated_power,
                            *frac_convective,
                            energy_supply_conn,
                            simulation_time.step_in_hours(),
                            control
                                .as_ref()
                                .and_then(|ctrl| controls.get_with_string(ctrl)),
                        ))
                    },
                    SpaceHeatSystemDetails::ElectricStorageHeater { .. } => return Err(NotImplementedError::new("Electric storage heater module not yet implemented.").into()), // requires implementation of ElecStorageHeater, make sure to add energy supply conn name to energy_conn_names_for_systems collection
                    SpaceHeatSystemDetails::WetDistribution { heat_source, temp_diff_emit_dsgn, control, thermal_mass, ecodesign_controller, design_flow_temp, .. } => {
                        // TODO 0.32 following are placeholder variables that we expect to come from emitters during migration to 0.32
                        let (c, n, frac_convective) = (0., 0., 0.);

                        let heat_source_name = &heat_source.name;
                        let temp_flow_limit_upper = &heat_source.temp_flow_limit_upper;

                        let energy_supply_conn_name = format!("{heat_source_name}_space_heating: {system_name}");
                        energy_conn_names_for_systems.insert(system_name.clone(), energy_supply_conn_name.clone());

                        let heat_source = heat_sources_wet.get(&heat_source.name).unwrap_or_else(|| panic!("A heat source name provided under the name '{heat_source_name}' was expected when setting up space heat systems in the calculation corpus."));
                        let mut with_buffer_tank = false;

                        let control = control
                            .as_ref()
                            .and_then(|ctrl| controls.get_with_string(ctrl)).expect("A control object was expected for a heat pump system");

                        let heat_source_service: SpaceHeatingService =
                        match heat_source {
                            WetHeatSource::HeatPump(heat_pump) => {
                                // TODO (from Python) If EAHP, feed zone volume into function below

                                // For HPs, checking if there's a buffer tank to inform both the service space heating
                                // and the emitters of its presence.
                                if heat_sources_wet_with_buffer_tank.contains(heat_source_name) {
                                    with_buffer_tank = true;
                                }

                                let volume_heated = total_volume_heated_by_system(zones, heat_system_name_for_zone, system_name);

                                let heat_source_service = HeatPump::create_service_space_heating(
                                    heat_pump.clone(),
                                    &energy_supply_conn_name,
                                    temp_flow_limit_upper.expect("Expected a temp_flow_limit_upper to be present for a heat pump"),
                                    *temp_diff_emit_dsgn, control,
                                    volume_heated);

                                if heat_pump.lock().source_is_exhaust_air() {
                                    // Record heating system as potentially requiring overventilation
                                    heat_system_names_requiring_overvent.push((*system_name).clone());
                                }
                                SpaceHeatingService::HeatPump(heat_source_service)
                            }
                            WetHeatSource::Boiler(boiler) => {
                                let heat_source_service = Boiler::create_service_space_heating(boiler.clone(), energy_supply_conn_name, control);
                                SpaceHeatingService::Boiler(heat_source_service)
                            }
                            WetHeatSource::Hiu(heat_network) => {
                                let heat_source_service = HeatNetwork::create_service_space_heating(heat_network.clone(), energy_supply_conn_name, control);
                                SpaceHeatingService::HeatNetwork(heat_source_service)
                            }
                            WetHeatSource::HeatBattery(heat_battery) => {
                                let heat_source_service = HeatBattery::create_service_space_heating(heat_battery.clone(), &energy_supply_conn_name, control);
                                SpaceHeatingService::HeatBattery(heat_source_service)
                            }
                        };
                        let temp_internal_air_fn = temp_internal_air_fn(temp_internal_air_accessor.clone());
                        // TODO Fix thermal mass logic as part of 0.32 migration
                        let space_heater = Emitters::new(thermal_mass.expect("Thermal mass may not be present while migrating to 0.32"), c, n, *temp_diff_emit_dsgn, frac_convective, Arc::new(RwLock::new(heat_source_service)), temp_internal_air_fn, external_conditions.clone(), *ecodesign_controller, *design_flow_temp as f64, with_buffer_tank, detailed_output_heating_cooling);
                        SpaceHeatSystem::WetDistribution(space_heater)
                    }
                    SpaceHeatSystemDetails::WarmAir {
                        frac_convective,
                        heat_source,
                        control,
                        ..
                    } => {
                        let heat_source_name = &heat_source.name;
                        let energy_supply_conn_name = format!("{heat_source_name}_space_heating: {system_name}");
                        energy_conn_names_for_systems.insert(system_name.clone(), energy_supply_conn_name.clone());
                        let heat_source = heat_sources_wet.get(&heat_source.name).unwrap_or_else(|| panic!("A heat source name provided under the name '{heat_source_name}' was expected when setting up space heat systems in the calculation corpus."));
                        let control = control
                            .as_ref()
                            .and_then(|ctrl| controls.get_with_string(ctrl)).expect("A control object was expected for a heat pump warm air system");

                        match heat_source {
                            WetHeatSource::HeatPump(heat_pump) => {
                                if heat_pump.lock().source_is_exhaust_air() {
                                    heat_system_names_requiring_overvent.push((*system_name).clone());
                                }
                                let volume_heated = total_volume_heated_by_system(zones, heat_system_name_for_zone, system_name);
                                SpaceHeatSystem::WarmAir(HeatPump::create_service_space_heating_warm_air(heat_pump.clone(), energy_supply_conn_name, control, *frac_convective, volume_heated).unwrap())
                            }
                            _ => panic!("The heat source referenced by details about warm air space heating with the name '{heat_source_name}' was expected to be a heat pump."),
                        }
                    }
                })),
            ))
        })
        .collect::<anyhow::Result<IndexMap<_, _>>>()?;
    Ok((space_heat_systems, energy_conn_names_for_systems))
}

type SpaceHeatSystemsWithEnergyConnections = (
    IndexMap<String, Arc<Mutex<SpaceHeatSystem>>>,
    IndexMap<String, String>,
);

fn space_cool_systems_from_input(
    input: &SpaceCoolSystemInput,
    cool_system_names_for_zone: Vec<&str>,
    controls: &Controls,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<IndexMap<String, AirConditioning>> {
    input
        .iter()
        .filter(|(system_name, _)| cool_system_names_for_zone.contains(&system_name.as_str()))
        .map(|(system_name, space_cool_system_details)| {
            if !matches!(
                space_cool_system_details.system_type,
                SpaceCoolSystemType::AirConditioning
            ) {
                unreachable!(
                    "There are no known space cool system types other than air conditioning."
                )
            }
            let SpaceCoolSystemDetails {
                cooling_capacity,
                efficiency,
                frac_convective,
                control,
                energy_supply,
                ..
            } = space_cool_system_details;
            let energy_supply = energy_supplies
                .ensured_get_for_type(*energy_supply, simulation_time_iterator.total_steps())?;
            let energy_supply_conn_name = system_name;
            let energy_supply_conn =
                EnergySupply::connection(energy_supply, energy_supply_conn_name).unwrap();
            let control = control
                .as_ref()
                .and_then(|ctrl| controls.get_with_string(ctrl));

            Ok((
                (*energy_supply_conn_name).clone(),
                AirConditioning::new(
                    *cooling_capacity,
                    *efficiency,
                    *frac_convective,
                    energy_supply_conn,
                    simulation_time_iterator.step_in_hours(),
                    control,
                ),
            ))
        })
        .collect::<anyhow::Result<IndexMap<_, _>>>()
}

fn on_site_generation_from_input(
    input: &OnSiteGeneration,
    energy_supplies: &mut IndexMap<String, Arc<RwLock<EnergySupply>>>,
    external_conditions: Arc<ExternalConditions>,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<IndexMap<String, PhotovoltaicSystem>> {
    input
        .iter()
        .map(|(name, generation_details)| {
            Ok(((*name).clone(), {
                let OnSiteGenerationDetails::PhotovoltaicSystem {
                    peak_power,
                    ventilation_strategy,
                    pitch,
                    orientation,
                    base_height,
                    height,
                    width,
                    energy_supply,
                    shading,
                    inverter_peak_power,
                    inverter_is_inside,
                    ..
                } = generation_details;
                let energy_supply = energy_supplies
                    .ensured_get_for_type(*energy_supply, simulation_time_iterator.total_steps())?;
                let energy_supply_conn = EnergySupply::connection(energy_supply, name).unwrap();
                PhotovoltaicSystem::new(
                    *peak_power,
                    *ventilation_strategy,
                    *pitch,
                    *orientation,
                    *base_height,
                    *height,
                    *width,
                    external_conditions.clone(),
                    energy_supply_conn,
                    simulation_time_iterator.step_in_hours(),
                    shading.clone(),
                    inverter_peak_power
                        .expect("Became optional during migration to 0.32 - correct for this"),
                    *inverter_is_inside,
                )
            }))
        })
        .collect::<anyhow::Result<IndexMap<_, _>>>()
}

fn total_volume_heated_by_system(
    zones: &Arc<IndexMap<String, Zone>>,
    heat_system_name_for_zone: &IndexMap<String, String>,
    heat_system_name: &str,
) -> f64 {
    zones
        .iter()
        .filter_map(|(z_name, zone)| {
            if let Some(system_name) = heat_system_name_for_zone.get(z_name) {
                (system_name == heat_system_name).then(|| zone.volume())
            } else {
                None
            }
        })
        .sum::<f64>()
}

fn required_vent_data_from_input(input: &ControlInput) -> anyhow::Result<Option<RequiredVentData>> {
    input
        .extra
        .get("required_vent")
        .map(|ctrl| {
            anyhow::Ok(RequiredVentData {
                schedule: expand_numeric_schedule(ctrl.numeric_schedule()?),
                start_day: ctrl.start_day(),
                time_series_step: ctrl.time_series_step(),
            })
        })
        .transpose()
}

#[derive(Clone, Debug)]
struct RequiredVentData {
    schedule: Vec<Option<f64>>,
    start_day: u32,
    time_series_step: f64,
}
