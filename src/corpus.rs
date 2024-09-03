use crate::compare_floats::max_of_2;
use crate::core::common::WaterSourceWithTemperature;
use crate::core::controls::time_control::{
    Control, ControlBehaviour, HeatSourceControl, OnOffMinimisingTimeControl, OnOffTimeControl,
    SetpointTimeControl, ToUChargeControl,
};
use crate::core::cooling_systems::air_conditioning::AirConditioning;
use crate::core::ductwork::Ductwork;
use crate::core::energy_supply::elec_battery::ElectricBattery;
use crate::core::energy_supply::energy_supply::{
    EnergySupplies, EnergySupply, EnergySupplyConnection,
};
use crate::core::energy_supply::pv::PhotovoltaicSystem;
use crate::core::heating_systems::boiler::{Boiler, BoilerServiceWaterCombi};
use crate::core::heating_systems::common::{HeatSourceWet, SpaceHeatSystem};
use crate::core::heating_systems::heat_battery::HeatBattery;
use crate::core::heating_systems::heat_network::{HeatNetwork, HeatNetworkServiceWaterDirect};
use crate::core::heating_systems::heat_pump::{HeatPump, HeatPumpHotWaterOnly};
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
    expand_boolean_schedule, expand_events, expand_numeric_schedule, NumericSchedule,
    ScheduleEvent, TypedScheduleEvent, WaterScheduleEventType,
};
use crate::core::space_heat_demand::building_element::{
    area_for_building_element_input, convert_uvalue_to_resistance, pitch_class, BuildingElement,
    BuildingElementAdjacentZTC, BuildingElementAdjacentZTUSimple, BuildingElementGround,
    BuildingElementOpaque, BuildingElementTransparent, HeatFlowDirection,
    NamedBuildingElementTransparent, PITCH_LIMIT_HORIZ_CEILING,
};
use crate::core::space_heat_demand::internal_gains::{
    ApplianceGains, EventApplianceGains, Gains, InternalGains,
};
use crate::core::space_heat_demand::thermal_bridge::{ThermalBridge, ThermalBridging};
use crate::core::space_heat_demand::ventilation::{
    AirTerminalDevices, CombustionAppliances, InfiltrationVentilation, MechanicalVentilation, Vent,
    Window,
};
use crate::core::space_heat_demand::ventilation_element::{
    air_change_rate_to_flow_rate, NaturalVentilation, VentilationElement,
    VentilationElementInfiltration, WholeHouseExtractVentilation, WindowOpeningForCooling,
};
use crate::core::space_heat_demand::zone::{AirChangesPerHourArgument, HeatBalance, Zone};
use crate::core::units::{
    kelvin_to_celsius, LITRES_PER_CUBIC_METRE, MILLIMETRES_IN_METRE, SECONDS_PER_HOUR,
    WATTS_PER_KILOWATT,
};
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::dhw_demand::{
    DemandVolTargetKey, DomesticHotWaterDemand, VolumeReference,
};
use crate::core::water_heat_demand::misc::water_demand_to_kwh;
use crate::external_conditions::ExternalConditions;
use crate::input::{
    init_orientation, ApplianceGains as ApplianceGainsInput, ApplianceGainsDetails,
    BuildingElement as BuildingElementInput, ColdWaterSourceDetails, ColdWaterSourceInput,
    ColdWaterSourceType, Control as ControlInput, ControlDetails, DuctShape, EnergyDiverter,
    EnergySupplyDetails, EnergySupplyInput, EnergySupplyKey, EnergySupplyType,
    ExternalConditionsInput, FloorType, FuelType, HeatPumpSourceType,
    HeatSource as HeatSourceInput, HeatSourceControl as HeatSourceControlInput,
    HeatSourceControlType, HeatSourceWetDetails, HeatSourceWetType, HotWaterSourceDetails,
    Infiltration, InfiltrationVentilation as InfiltrationVentilationInput, Input,
    InternalGains as InternalGainsInput, InternalGainsDetails, OnSiteGeneration,
    OnSiteGenerationDetails, SpaceCoolSystem as SpaceCoolSystemInput, SpaceCoolSystemDetails,
    SpaceCoolSystemType, SpaceHeatSystem as SpaceHeatSystemInput, SpaceHeatSystemDetails,
    TerrainClass, ThermalBridging as ThermalBridgingInput, ThermalBridgingDetails, VentType,
    Ventilation, VentilationLeaks, VentilationShieldClass, WasteWaterHeatRecovery,
    WasteWaterHeatRecoveryDetails, WaterHeatingEvent, WaterHeatingEvents,
    WindowOpeningForCooling as WindowOpeningForCoolingInput, WwhrsType, ZoneDictionary, ZoneInput,
};
use crate::simulation_time::{SimulationTime, SimulationTimeIteration, SimulationTimeIterator};
use anyhow::{anyhow, bail};
use arrayvec::ArrayString;
use indexmap::IndexMap;
#[cfg(feature = "indicatif")]
use indicatif::ProgressIterator;
use parking_lot::{Mutex, RwLock};
use serde_json::Value;
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::fmt::{Display, Formatter};
use std::hash::Hash;
use std::sync::Arc;

// TODO make this a runtime parameter?
const DETAILED_OUTPUT_HEATING_COOLING: bool = true;

/// As of Rust 1.82 we'll be able to declare this using constants as it is due to support floating-point arithmetic at compile time
fn temp_setpnt_heat_none() -> f64 {
    kelvin_to_celsius(0.0)
}

fn temp_setpnt_cool_none() -> f64 {
    kelvin_to_celsius(1.4e32)
}

// used for calculations re internal gains from pipework
const FRAC_DHW_ENERGY_INTERNAL_GAINS: f64 = 0.25;

pub struct Corpus {
    pub simulation_time: Arc<SimulationTimeIterator>,
    pub external_conditions: Arc<ExternalConditions>,
    pub infiltration: VentilationElementInfiltration,
    pub cold_water_sources: ColdWaterSources,
    pub energy_supplies: EnergySupplies,
    pub internal_gains: InternalGainsCollection,
    pub controls: Controls,
    pub wwhrs: IndexMap<String, Arc<Mutex<Wwhrs>>>,
    pub event_schedules: HotWaterEventSchedules,
    pub domestic_hot_water_demand: DomesticHotWaterDemand,
    pub ventilation: Option<Arc<Mutex<VentilationElement>>>,
    pub space_heating_ductwork: Option<Ductwork>,
    pub zones: Arc<IndexMap<String, Zone>>,
    pub energy_supply_conn_unmet_demand_zone: IndexMap<String, Arc<EnergySupplyConnection>>,
    pub heat_system_name_for_zone: IndexMap<String, String>,
    pub cool_system_name_for_zone: IndexMap<String, String>,
    pub total_floor_area: f64,
    pub total_volume: f64,
    pub wet_heat_sources: IndexMap<String, Arc<Mutex<WetHeatSource>>>,
    pub hot_water_sources: IndexMap<String, HotWaterSource>,
    pub heat_system_names_requiring_overvent: Vec<String>,
    pub space_heat_systems: IndexMap<String, Arc<Mutex<SpaceHeatSystem>>>,
    pub space_cool_systems: IndexMap<String, AirConditioning>,
    pub on_site_generation: IndexMap<String, PhotovoltaicSystem>,
    pub diverters: Vec<Arc<RwLock<PVDiverter>>>,
    required_vent_data: Option<RequiredVentData>,
    energy_supply_conn_names_for_hot_water_source: IndexMap<String, Vec<String>>,
    energy_supply_conn_names_for_heat_systems: IndexMap<String, String>,
    timestep_end_calcs: Vec<Arc<Mutex<WetHeatSource>>>,
}

impl Corpus {
    pub fn from_inputs(
        input: &Input,
        external_conditions: Option<ExternalConditions>,
    ) -> anyhow::Result<Self> {
        let simulation_time_iterator = Arc::new(input.simulation_time.iter());

        let external_conditions = Arc::new(match external_conditions {
            Some(external_conditions) => external_conditions,
            None => external_conditions_from_input(
                input.external_conditions.clone(),
                &simulation_time_iterator,
            ),
        });

        let diverter_types: DiverterTypes = (&input.energy_supply).into();
        let mut diverters: Vec<Arc<RwLock<PVDiverter>>> = Default::default();

        let cold_water_sources =
            cold_water_sources_from_input(&input.cold_water_source, &input.simulation_time);
        let wwhrs = wwhrs_from_input(
            input.waste_water_heat_recovery.as_ref(),
            &cold_water_sources,
        );

        let mut energy_supplies = energy_supplies_from_input(
            &input.energy_supply,
            simulation_time_iterator.clone().as_ref(),
            external_conditions.clone(),
        );

        let controls =
            control_from_input(&input.control, simulation_time_iterator.clone().as_ref());

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
            vec![], // use empty while migrating to 0.30
        )?;

        let infiltration = infiltration_from_input(input.infiltration.as_ref().unwrap());

        let space_heating_ductwork = ductwork_from_ventilation_input(&input.ventilation);

        let ventilation = input
            .ventilation
            .as_ref()
            .map(|v| {
                anyhow::Ok(Arc::new(Mutex::new(ventilation_from_input(
                    "Ventilation system",
                    v,
                    &infiltration,
                    simulation_time_iterator.clone().as_ref(),
                    &mut energy_supplies,
                )?)))
            })
            .transpose()?;

        let opening_area_total_from_zones = opening_area_total_from_zones(&input.zone);

        let total_volume = input.zone.values().map(|zone| zone.volume).sum::<f64>();

        let mut heat_system_name_for_zone: IndexMap<String, String> = Default::default();
        let mut cool_system_name_for_zone: IndexMap<String, String> = Default::default();

        // infiltration ventilation
        let (
            infiltration_ventilation,
            window_adjust_control,
            mechanical_ventilations,
            space_heating_ductwork_new,
        ) = infiltration_ventilation_from_input(
            &input.zone,
            &input.infiltration_ventilation,
            &controls,
            &mut energy_supplies,
            simulation_time_iterator.as_ref(),
            total_volume,
            external_conditions.clone(),
        )?;

        let infiltration_ventilation = Arc::from(infiltration_ventilation);
        let mechanical_ventilations: IndexMap<String, Arc<MechanicalVentilation>> =
            mechanical_ventilations
                .into_iter()
                .map(|(name, mech_vent)| (name.to_owned(), Arc::from(mech_vent)))
                .collect();

        let required_vent_data = required_vent_data_from_input(&input.control);

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
                        simulation_time_iterator.clone().as_ref(),
                    )?;
                    if let Some(heat_system_name) = heat_system_name {
                        heat_system_name_for_zone.insert((*i).clone(), heat_system_name);
                    }
                    if let Some(cool_system_name) = cool_system_name {
                        cool_system_name_for_zone.insert((*i).clone(), cool_system_name);
                    }

                    zone_for_corpus
                }))
            })
            .collect::<anyhow::Result<_>>()?;
        let zones = Arc::new(zones);

        if !has_unique_values(&heat_system_name_for_zone)
            || !has_unique_values(&cool_system_name_for_zone)
        {
            bail!("the heat or cool systems do not have unique names in the inputs");
        }

        let energy_supply_conn_unmet_demand_zone = set_up_energy_supply_unmet_demand_zones(
            energy_supplies.unmet_demand.clone(),
            &input.zone,
        );
        // TODO: there needs to be some equivalent here of the Python code that builds the dict __energy_supply_conn_unmet_demand_zone

        let total_floor_area = zones.values().fold(0., |acc, zone| zone.area() + acc);

        let mut internal_gains = internal_gains_from_input(&input.internal_gains, total_floor_area);

        apply_appliance_gains_from_input(
            &mut internal_gains,
            &input.appliance_gains,
            &mut energy_supplies,
            total_floor_area,
            simulation_time_iterator.as_ref(),
        )?;

        let mut timestep_end_calcs = vec![];

        let wet_heat_sources: IndexMap<String, Arc<Mutex<WetHeatSource>>> = input
            .heat_source_wet
            .clone()
            .unwrap_or_default()
            .iter()
            .map(|(name, heat_source_wet_details)| {
                let ventilation = ventilation.as_ref().map(|ventilation| ventilation.lock());
                let heat_source = Arc::new(Mutex::new(heat_source_wet_from_input(
                    name,
                    (*heat_source_wet_details).clone(),
                    external_conditions.clone(),
                    simulation_time_iterator.clone(),
                    ventilation.map(|v| (*v).clone()),
                    input.ventilation.as_ref().map(|v| v.req_ach()),
                    &mechanical_ventilations,
                    zones.len(),
                    TempInternalAirAccessor {
                        zones: zones.clone(),
                        total_volume,
                    },
                    total_volume,
                    &controls,
                    &mut energy_supplies,
                )?));
                match *heat_source.lock() {
                    WetHeatSource::HeatPump(_)
                    | WetHeatSource::Boiler(_)
                    | WetHeatSource::HeatBattery(_) => {
                        timestep_end_calcs.push(heat_source.clone());
                    }
                    _ => {}
                }
                anyhow::Ok(((*name).clone(), heat_source))
            })
            .collect::<anyhow::Result<IndexMap<String, Arc<Mutex<WetHeatSource>>>>>()?;

        let mut hot_water_sources: IndexMap<String, HotWaterSource> = Default::default();
        let mut energy_supply_conn_names_for_hot_water_source: IndexMap<String, Vec<String>> =
            Default::default();
        let (hot_water_source, hw_cylinder_conn_names) = hot_water_source_from_input(
            "hw cylinder".to_string(),
            &input.hot_water_source.hot_water_cylinder,
            &cold_water_sources,
            &wet_heat_sources,
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
                    &Default::default(),
                    &mut heat_system_names_requiring_overvent,
                    &heat_system_name_for_zone,
                    &zones,
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
            infiltration,
            cold_water_sources,
            energy_supplies,
            internal_gains,
            controls,
            wwhrs,
            event_schedules,
            domestic_hot_water_demand,
            ventilation,
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
            space_heat_systems,
            space_cool_systems,
            on_site_generation,
            diverters,
            required_vent_data,
            energy_supply_conn_names_for_hot_water_source,
            energy_supply_conn_names_for_heat_systems,
            timestep_end_calcs,
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
        no_events: usize,
        temp_final_drawoff: f64,
        temp_average_drawoff: f64,
        temp_hot_water: f64,
        vol_hot_water_equiv_elec_shower: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (f64, f64, f64) {
        todo!()
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
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (f64, f64, f64) {
        let frac_dhw_energy_internal_gains = 0.25;

        let (pw_losses_internal, pw_losses_external) = self.calc_pipework_losses(
            delta_t_h,
            hw_duration,
            no_of_hw_events,
            simulation_time_iteration,
        );

        let gains_internal_dhw_use = frac_dhw_energy_internal_gains
            * water_demand_to_kwh(
                vol_hot_water_at_tapping_point,
                52.,
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
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (f64, f64) {
        let demand_water_temperature = 52.;
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

    fn calc_ductwork_losses(
        &self,
        _t_idx: usize,
        _delta_t_h: f64,
        efficiency: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        // assume 100% efficiency
        // i.e. temp inside the supply and extract ducts is room temp and temp inside exhaust and intake is external temp
        // assume MVHR unit is running 100% of the time
        let internal_air_temperature = self.temp_internal_air();

        // Calculate heat loss from ducts when unit is inside
        // Air temp inside ducts increases, heat lost from dwelling
        let ductwork = &self.space_heating_ductwork;
        match ductwork {
            None => 0.,
            Some(_ductwork) => {
                // MVHR duct temperatures:
                // extract_duct_temp - indoor air temperature
                // intake_duct_temp - outside air temperature

                let intake_duct_temp = self
                    .external_conditions
                    .air_temp(&simulation_time_iteration);

                let temp_diff = internal_air_temperature - intake_duct_temp;

                // Supply duct contains what the MVHR could recover
                let _supply_duct_temp = intake_duct_temp + (efficiency * temp_diff);

                // Exhaust duct contans the heat that couldn't be recovered
                let _exhaust_duct_temp = intake_duct_temp + ((1. - efficiency) * temp_diff);

                // comment out for now while migrating to 0.30
                // ductwork
                //     .total_duct_heat_loss(
                //         Some(internal_air_temperature),
                //         Some(supply_duct_temp),
                //         Some(internal_air_temperature),
                //         Some(intake_duct_temp),
                //         Some(exhaust_duct_temp),
                //         efficiency,
                //     )
                //     .unwrap()
                0.
            }
        }
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
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<SpaceHeatingCalculation> {
        let temp_ext_air = self
            .external_conditions
            .air_temp(&simulation_time_iteration);
        // Calculate timestep in seconds
        let delta_t = delta_t_h * SECONDS_PER_HOUR as f64;

        let (ductwork_losses, ductwork_losses_per_m3) = self
            .ventilation
            .as_ref()
            .map(|ventilation| match &*ventilation.lock() {
                VentilationElement::Mvhr(mvhr) => {
                    let ductwork_losses = self.calc_ductwork_losses(
                        0,
                        delta_t_h,
                        mvhr.efficiency(),
                        simulation_time_iteration,
                    );
                    let ductwork_losses_per_m3 = ductwork_losses / self.total_volume;

                    (ductwork_losses, ductwork_losses_per_m3)
                }
                _ => Default::default(),
            })
            .unwrap_or_default();

        // Calculate internal and and solar gains for each zone
        let mut gains_internal_zone: HashMap<&str, f64> = Default::default();
        let mut gains_solar_zone: HashMap<&str, f64> = Default::default();

        for (z_name, zone) in self.zones.iter() {
            // Initialise to dhw internal gains split proportionally to zone floor area
            let mut gains_internal_zone_inner =
                gains_internal_dhw * zone.area() / self.total_floor_area;
            for gains in self.internal_gains.values() {
                gains_internal_zone_inner +=
                    gains.total_internal_gain_in_w(zone.area(), simulation_time_iteration)?;
            }
            let gains_internal_zone_entry = gains_internal_zone
                .entry(z_name)
                .or_insert(gains_internal_zone_inner);
            // Add gains from ventilation fans (also calculates elec demand from fans)
            // TODO (from Python) Remove the branch on the type of ventilation (find a better way)
            match &self.ventilation {
                None => {}
                Some(ventilation) => {
                    let mut ventilation = ventilation.lock();
                    if !matches!(*ventilation, VentilationElement::Natural(_)) {
                        *gains_internal_zone_entry +=
                            ventilation.fans(zone.volume(), simulation_time_iteration.index, None)
                                + ductwork_losses_per_m3 * zone.volume();
                    }
                }
            }
            gains_solar_zone.insert(z_name, zone.gains_solar(simulation_time_iteration));
        }

        // Calculate space heating and cooling demand for each zone and sum
        // Keep track of how much is from each zone, so that energy provided
        // can be split between them in same proportion later
        let (
            mut space_heat_demand_system,
            mut space_cool_demand_system,
            mut space_heat_demand_zone,
            mut space_cool_demand_zone,
            mut _h_ve_cool_extra_zone,
        ) = self.space_heat_cool_demand_by_system_and_zone(
            delta_t_h,
            temp_ext_air,
            &gains_internal_zone,
            &gains_solar_zone,
            None,
            simulation_time_iteration,
        );

        // If any heating systems potentially require overventilation,
        // calculate running time and throughput factor for all services
        // combined based on space heating demand assuming no overventilation
        let mut space_heat_running_time_cumulative = 0.0;
        let mut throughput_factor = 1.0;
        for (heat_system_name, heat_system) in &self.space_heat_systems {
            if self
                .heat_system_names_requiring_overvent
                .contains(heat_system_name)
            {
                (space_heat_running_time_cumulative, throughput_factor) = heat_system
                    .lock()
                    .running_time_throughput_factor(
                        space_heat_demand_system[heat_system_name],
                        space_heat_running_time_cumulative,
                        simulation_time_iteration,
                    )
                    .unwrap();
            }
        }

        // If there is overventilation due to heating or hot water system (e.g.
        // exhaust air heat pump) then recalculate space heating/cooling demand
        // with additional ventilation calculated based on throughput factor
        // based on original space heating demand calculation. Note the
        // additional ventilation throughput is the result of the HP running
        // to satisfy both space and water heating demand but will affect
        // space heating demand only
        // TODO (from Python) The space heating demand is only recalculated once, rather
        //                    than feeding back in to the throughput factor calculation
        //                    above to get a further-refined space heating demand. This is
        //                    consistent with the approach in SAP 10.2 and keeps the
        //                    execution time of the calculation bounded. However, the
        //                    merits of iterating over this calculation until converging on
        //                    a solution should be considered in the future.
        if throughput_factor > 1.0 {
            for (z_name, zone) in self.zones.iter() {
                // Add additional gains from ventilation fans
                match &self.ventilation {
                    None => {}
                    Some(ventilation) => {
                        let mut ventilation = ventilation.lock();
                        if !matches!(*ventilation, VentilationElement::Natural(_)) {
                            *gains_internal_zone.get_mut(z_name.as_str()).unwrap() += ventilation
                                .fans(
                                    zone.volume(),
                                    simulation_time_iteration.index,
                                    Some(throughput_factor - 1.0),
                                );
                        }
                    }
                }
            }
            (
                space_heat_demand_system,
                space_cool_demand_system,
                space_heat_demand_zone,
                space_cool_demand_zone,
                _h_ve_cool_extra_zone,
            ) = self.space_heat_cool_demand_by_system_and_zone(
                delta_t_h,
                temp_ext_air,
                &gains_internal_zone,
                &gains_solar_zone,
                Some(throughput_factor),
                simulation_time_iteration,
            );
        }

        // pre-calc whether different space heat/cool systems are in required period, and frac_convective before self.space_heat_systems is mutably borrowed
        let space_heat_systems_in_required_period =
            self.space_heat_systems_in_required_period(simulation_time_iteration);
        let space_heat_systems_frac_convective = self.space_heat_systems_frac_convective();
        let space_cool_systems_in_required_period =
            self.space_cool_systems_in_required_period(simulation_time_iteration);
        let space_cool_systems_frac_convective = self.space_cool_systems_frac_convective();

        // Calculate how much heating the systems can provide
        let mut space_heat_provided: HashMap<&str, f64> = Default::default();
        for (heat_system_name, heat_system) in &self.space_heat_systems {
            space_heat_provided.insert(
                heat_system_name.as_str(),
                heat_system
                    .lock()
                    .demand_energy(
                        space_heat_demand_system[heat_system_name.as_str()],
                        simulation_time_iteration,
                    )
                    .unwrap(),
            );
        }

        // Calculate how much cooling the systems can provide
        let mut space_cool_provided: HashMap<&str, f64> = Default::default();
        for (cool_system_name, cool_system) in &self.space_cool_systems {
            space_cool_provided.insert(
                cool_system_name.as_str(),
                cool_system.demand_energy(
                    space_cool_demand_system[cool_system_name.as_str()],
                    simulation_time_iteration,
                ),
            );
        }

        // Apportion the provided heating/cooling between the zones in
        // proportion to the heating/cooling demand in each zone. Then
        // update resultant temperatures in zones.
        let mut internal_air_temp: HashMap<&str, f64> = Default::default();
        let mut operative_temp: HashMap<&str, f64> = Default::default();
        let mut heat_balance_map: HashMap<&str, Option<HeatBalance>> = Default::default(); // using unit type here as placeholder
        for (z_name, zone) in self.zones.iter() {
            // Look up names of relevant heating and cooling systems for this zone
            let h_name = self.heat_system_name_for_zone.get(z_name.as_str());
            let c_name = self.cool_system_name_for_zone.get(z_name.as_str());

            // If zone is unheated or there was no demand on heating system,
            // set heating gains for zone to zero, else calculate
            let gains_heat = match h_name {
                None => 0.0,
                Some(h_name) => space_heat_provided[h_name.as_str()],
            };

            // If zone is uncooled or there was no demand on cooling system,
            // set cooling gains for zone to zero, else calculate
            let gains_cool = match c_name {
                None => 0.0,
                Some(c_name) => space_cool_provided[c_name.as_str()],
            };

            // Sum heating gains (+ve) and cooling gains (-ve) and convert from kWh to W
            let gains_heat_cool = (gains_heat + gains_cool) * WATTS_PER_KILOWATT as f64 / delta_t_h;

            // Calculate how much space heating / cooling demand is unmet
            // Note: Demand is not considered unmet if it is outside the
            //        required heating/cooling period (which does not include
            //        times when the system is on due to setback or advanced
            //        start)
            // Note: Need to check that demand is non-zero, to avoid
            //        reporting unmet demand when heating system is absorbing
            //        energy from zone or cooling system is releasing energy
            //        to zone, which may be the case in some timesteps for
            //        systems with significant thermal mass.
            let in_req_heat_period = match h_name {
                None => false,
                Some(h_name) => {
                    space_heat_systems_in_required_period[h_name.as_str()].unwrap_or(false)
                }
            };
            let space_heat_demand_zone_current = space_heat_demand_zone[z_name.as_str()];
            if in_req_heat_period && space_heat_demand_zone_current > 0. {
                let _energy_shortfall_heat =
                    max_of_2(0., space_heat_demand_zone_current - gains_heat);
                // TODO report energy supply unmet demand
            }
            let in_req_cool_period = match c_name {
                None => false,
                Some(c_name) => {
                    space_cool_systems_in_required_period[c_name.as_str()].unwrap_or(false)
                }
            };
            let space_cool_demand_zone_current = space_cool_demand_zone[z_name.as_str()];
            if in_req_cool_period && space_cool_demand_zone_current > 0. {
                let _energy_shortfall_cool =
                    max_of_2(0., -(space_cool_demand_zone_current - gains_cool));
                // TODO report energy supply unmet demand
            }

            // Look up convective fraction for heating/cooling for this zone
            // Note: gains_heat could be negative (or gains_cool could be
            //       positive) if thermal mass of emitters causes e.g. the
            //       heat emitters to absorb energy from the zone.
            let frac_convective = if gains_heat != 0. {
                space_heat_systems_frac_convective
                    [h_name.expect("name for space heat system expected to be set")]
            } else if gains_cool != 0. {
                space_cool_systems_frac_convective
                    [c_name.expect("name for space cool system expected to be set")]
            } else {
                1.0
            };

            heat_balance_map.insert(
                z_name.as_str(),
                zone.update_temperatures(
                    delta_t,
                    temp_ext_air,
                    gains_internal_zone[z_name.as_str()],
                    gains_solar_zone[z_name.as_str()],
                    gains_heat_cool,
                    frac_convective,
                    Default::default(),
                    Default::default(), // temporary defaults during migrating to 0.30
                    simulation_time_iteration,
                ),
            );

            // In the Python code we handle an edge case here
            // and add a None key to these dictionaries
            // This is currently handlded in this codebase
            // when we generate the output file instead.
            // TODO move the edge case logic here to match the Python code
            // as closely as possible

            internal_air_temp.insert(z_name.as_str(), zone.temp_internal_air());
            operative_temp.insert(z_name.as_str(), zone.temp_operative());
        }

        Ok((
            gains_internal_zone,
            gains_solar_zone,
            operative_temp,
            internal_air_temp,
            space_heat_demand_zone,
            space_cool_demand_zone,
            space_heat_demand_system,
            space_cool_demand_system,
            space_heat_provided,
            space_cool_provided,
            ductwork_losses,
            heat_balance_map,
        ))
    }

    pub fn run(&mut self) -> RunResults {
        let simulation_time = self.simulation_time.as_ref().to_owned();
        let vec_capacity = || Vec::with_capacity(simulation_time.total_steps());

        let mut timestep_array = vec_capacity();
        let mut gains_internal_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut gains_solar_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut operative_temp_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut internal_air_temp_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut space_heat_demand_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut space_cool_demand_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut space_heat_demand_system_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut space_cool_demand_system_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut space_heat_provided_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
        let mut space_cool_provided_dict: IndexMap<KeyString, Vec<f64>> = Default::default();
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
        let mut heat_balance_all_dict: IndexMap<
            KeyString,
            IndexMap<KeyString, IndexMap<KeyString, f64>>,
        > = IndexMap::from([
            ("air_node".try_into().unwrap(), Default::default()),
            ("internal_boundary".try_into().unwrap(), Default::default()),
            ("external_boundary".try_into().unwrap(), Default::default()),
        ]);
        let heat_source_wet_results_dict = Default::default();
        let heat_source_wet_results_annual_dict = Default::default();

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
                heat_balance_value.insert(z_name, Default::default());
            }
        }

        for h_name in self.heat_system_name_for_zone.values() {
            let h_name = h_name.as_str().try_into().unwrap();
            space_heat_demand_system_dict.insert(h_name, vec_capacity());
            space_heat_provided_dict.insert(h_name, vec_capacity());
        }

        for c_name in self.cool_system_name_for_zone.values() {
            let c_name = c_name.as_str().try_into().unwrap();
            space_cool_demand_system_dict.insert(c_name, vec_capacity());
            space_cool_provided_dict.insert(c_name, vec_capacity());
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

        #[cfg(feature = "indicatif")]
        let simulation_time_iter = simulation_time.progress();
        #[cfg(not(feature = "indicatif"))]
        let simulation_time_iter = simulation_time;

        for t_it in simulation_time_iter {
            timestep_array.push(t_it.time);
            let temp_hot_water = 0.; // TODO as part of migration: implement hw_source.get_temp_hot_water()
            let mut temp_final_drawoff = temp_hot_water;
            let mut temp_average_drawoff = temp_hot_water;
            let mut pw_losses_internal;
            let mut pw_losses_external;
            let mut gains_internal_dhw_use;
            let mut hw_energy_output;

            let (
                hw_demand_vol,
                hw_demand_vol_target,
                hw_vol_at_tapping_points,
                hw_duration,
                no_events,
                hw_energy_demand,
                usage_events,
                vol_hot_water_equiv_elec_shower,
            ) = self
                .domestic_hot_water_demand
                .hot_water_demand(t_it.index, 52.0); // temporary value to be changed to hot water temp from hw cylinder source while migrating to 0.30

            let hw_source = self.hot_water_sources.get_mut("hw cylinder").unwrap();
            match hw_source {
                HotWaterSource::StorageTank(source) => {
                    let volume_water_remove_from_tank;

                    (
                        hw_energy_output,
                        _, // Python has an unused unmet_demand variable here
                        temp_final_drawoff,
                        temp_average_drawoff,
                        volume_water_remove_from_tank,
                    ) = source
                        .lock()
                        .demand_hot_water(usage_events.expect("usage_events was not set"), t_it);

                    (
                        pw_losses_internal,
                        pw_losses_external,
                        gains_internal_dhw_use,
                    ) = self.pipework_losses_and_internal_gains_from_hw_storage_tank(
                        t_it.timestep,
                        volume_water_remove_from_tank,
                        hw_duration,
                        no_events,
                        temp_final_drawoff,
                        temp_average_drawoff,
                        temp_hot_water,
                        vol_hot_water_equiv_elec_shower,
                        t_it,
                    );
                }
                _ => {
                    hw_energy_output = hw_source.demand_hot_water(hw_demand_vol_target, t_it);

                    (
                        pw_losses_internal,
                        pw_losses_external,
                        gains_internal_dhw_use,
                    ) = self.pipework_losses_and_internal_gains_from_hw(
                        t_it.timestep,
                        hw_vol_at_tapping_points,
                        hw_duration,
                        no_events,
                        t_it,
                    );
                }
            }

            // Convert from litres to kWh
            let cold_water_source = self.hot_water_sources["hw cylinder"]
                .get_cold_water_source()
                .expect("expected cold water source to be available on hot water cylinder");
            let cold_water_temperature = cold_water_source.temperature(t_it.index);
            let hw_energy_demand_incl_pipework_loss = water_demand_to_kwh(
                hw_demand_vol,
                // assumed cold water temperature
                52.0,
                cold_water_temperature,
            );

            let mut gains_internal_dhw = (pw_losses_internal + gains_internal_dhw_use)
                * WATTS_PER_KILOWATT as f64
                / t_it.timestep;
            match self.hot_water_sources.get_mut("hw cylinder").unwrap() {
                HotWaterSource::StorageTank(ref mut source) => {
                    gains_internal_dhw += source.lock().internal_gains();
                }
                HotWaterSource::CombiBoiler(ref mut source) => {
                    gains_internal_dhw += source.internal_gains();
                }
                _ => {}
            }

            let (
                gains_internal_zone,
                gains_solar_zone,
                operative_temp,
                internal_air_temp,
                space_heat_demand_zone,
                space_cool_demand_zone,
                space_heat_demand_system,
                space_cool_demand_system,
                space_heat_provided,
                space_cool_provided,
                ductwork_gains,
                heat_balance_dict,
            ) = self
                .calc_space_heating(t_it.timestep, gains_internal_dhw, t_it)
                .expect("Expected a space heating calculation to be possible.");

            // Perform calculations that can only be done after all heating
            // services have been calculated
            for system in &self.timestep_end_calcs {
                system.lock().timestep_end(t_it);
            }

            for (z_name, gains_internal) in gains_internal_zone {
                gains_internal_dict
                    .get_mut(z_name)
                    .unwrap()
                    .push(gains_internal);
            }

            for (z_name, gains_solar) in gains_solar_zone {
                gains_solar_dict.get_mut(z_name).unwrap().push(gains_solar);
            }

            for (z_name, temp) in operative_temp {
                operative_temp_dict.get_mut(z_name).unwrap().push(temp);
            }

            for (z_name, temp) in internal_air_temp {
                internal_air_temp_dict.get_mut(z_name).unwrap().push(temp);
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
                    .get_mut(h_name)
                    .unwrap()
                    .push(output);
            }

            for (c_name, output) in space_cool_provided {
                space_cool_provided_dict
                    .get_mut(c_name)
                    .unwrap()
                    .push(output);
            }

            for (_z_name, hb_dict) in heat_balance_dict {
                if hb_dict.is_some() {
                    // TODO complete implementation here
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

            // loop through on-site energy generation
            for gen in self.on_site_generation.values() {
                // Get energy produced for the current timestep
                gen.produce_energy(t_it);
            }

            self.energy_supplies
                .calc_energy_import_export_betafactor(t_it);

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
        for (name, supply) in self.energy_supplies.supplies_by_name() {
            let supply = supply.read();
            let name: KeyString = name.try_into().unwrap();
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

        let hot_water_energy_out: IndexMap<KeyString, Vec<f64>> = IndexMap::from([(
            "hw cylinder".try_into().unwrap(),
            hot_water_energy_output_dict
                .get("energy_output")
                .unwrap()
                .to_owned(),
        )]);
        // TODO replace in energy supply names when available
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
            ("Internal gains".try_into().unwrap(), gains_internal_dict),
            ("Solar gains".try_into().unwrap(), gains_solar_dict),
            ("Operative temp".try_into().unwrap(), operative_temp_dict),
            (
                "Internal air temp".try_into().unwrap(),
                internal_air_temp_dict,
            ),
            (
                "Space heat demand".try_into().unwrap(),
                space_heat_demand_dict,
            ),
            (
                "Space cool demand".try_into().unwrap(),
                space_cool_demand_dict,
            ),
        ]);
        let hc_system_dict = IndexMap::from([
            (
                "Heating system".try_into().unwrap(),
                space_heat_demand_system_dict,
            ),
            (
                "Cooling system".try_into().unwrap(),
                space_cool_demand_system_dict,
            ),
            (
                "Heating system output".try_into().unwrap(),
                space_heat_provided_dict,
            ),
            (
                "Cooling system output".try_into().unwrap(),
                space_cool_provided_dict,
            ),
        ]);
        let hot_water_dict = IndexMap::from([
            (
                "Hot water demand".try_into().unwrap(),
                HotWaterResultMap::Float(hot_water_demand_dict),
            ),
            (
                "Hot water energy demand".try_into().unwrap(),
                HotWaterResultMap::Float(hot_water_energy_demand_dict),
            ),
            (
                "Hot water energy demand incl pipework_loss"
                    .try_into()
                    .unwrap(),
                HotWaterResultMap::Float(hot_water_energy_demand_dict_incl_pipework),
            ),
            (
                "Hot water duration".try_into().unwrap(),
                HotWaterResultMap::Float(hot_water_duration_dict),
            ),
            (
                "Hot Water Events".try_into().unwrap(),
                HotWaterResultMap::Int(hot_water_no_events_dict),
            ),
            (
                "Pipework losses".try_into().unwrap(),
                HotWaterResultMap::Float(hot_water_pipework_dict),
            ),
        ]);

        // Report detailed outputs from heat source wet objects, if requested and available
        // TODO implement once detailed_output_heating_cooling instance var implemented

        (
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
            ductwork_gains_dict,
            heat_balance_all_dict,
            heat_source_wet_results_dict,
            heat_source_wet_results_annual_dict,
        )
    }

    fn space_heat_systems_in_required_period(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> HashMap<String, Option<bool>> {
        self.space_heat_systems
            .iter()
            .map(|(system_name, system)| {
                let system = system.lock();
                (
                    system_name.clone(),
                    system.in_required_period(simulation_time_iteration),
                )
            })
            .collect()
    }

    fn space_heat_systems_frac_convective(&self) -> HashMap<String, f64> {
        self.space_heat_systems
            .iter()
            .map(|(system_name, system)| {
                let system = system.lock();
                (system_name.clone(), system.frac_convective())
            })
            .collect()
    }

    fn space_cool_systems_in_required_period(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> HashMap<String, Option<bool>> {
        self.space_cool_systems
            .iter()
            .map(|(system_name, system)| {
                (
                    system_name.clone(),
                    system.in_required_period(&simulation_time_iteration),
                )
            })
            .collect()
    }

    fn space_cool_systems_frac_convective(&self) -> HashMap<String, f64> {
        self.space_cool_systems
            .iter()
            .map(|(system_name, system)| (system_name.clone(), system.frac_convective()))
            .collect()
    }

    /// Calculate overall CoP over calculation period for each heating and cooling system
    fn heat_cool_cop(
        &self,
        energy_provided: &IndexMap<KeyString, Vec<f64>>,
        results_end_user: &IndexMap<KeyString, IndexMap<String, Vec<f64>>>,
        energy_supply_conn_name_for_space_hc_system: IndexMap<String, Vec<String>>,
    ) -> IndexMap<KeyString, NumberOrDivisionByZero> {
        let mut hc_output_overall: IndexMap<KeyString, f64> = Default::default();
        let mut hc_input_overall: IndexMap<KeyString, f64> = Default::default();
        let mut cop_dict: IndexMap<KeyString, NumberOrDivisionByZero> = Default::default();
        for (hc_name, hc_output) in energy_provided {
            hc_output_overall.insert(*hc_name, hc_output.iter().sum::<f64>().abs());
            hc_input_overall.insert(*hc_name, 0.);
            let energy_supply_conn_names =
                energy_supply_conn_name_for_space_hc_system[hc_name.as_str()].clone();
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
                *hc_name,
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

    /// Calculate space heating and cooling demand for each zone and sum.
    ///
    /// Keep track of how much is from each zone, so that energy provided
    /// can be split between them in same proportion later
    fn space_heat_cool_demand_by_system_and_zone(
        &self,
        delta_t_h: f64,
        temp_ext_air: f64,
        gains_internal_zone: &HashMap<&str, f64>,
        gains_solar_zone: &HashMap<&str, f64>,
        _throughput_factor: Option<f64>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (NumberMap, NumberMap, NumberMap, NumberMap, NumberMap) {
        let mut space_heat_demand_system: HashMap<String, f64> = Default::default();
        for heat_system_name in self.space_heat_systems.keys() {
            space_heat_demand_system.insert((*heat_system_name).clone(), 0.0);
        }

        let mut space_cool_demand_system: HashMap<String, f64> = Default::default();
        for cool_system_name in self.space_cool_systems.keys() {
            space_cool_demand_system.insert((*cool_system_name).clone(), 0.0);
        }

        let mut space_heat_demand_zone: HashMap<String, f64> = Default::default();
        let mut space_cool_demand_zone: HashMap<String, f64> = Default::default();
        let mut h_ve_cool_extra_zone: HashMap<String, f64> = Default::default();
        for (z_name, zone) in self.zones.iter() {
            // Look up names of relevant heating and cooling systems for this zone
            let h_name = self.heat_system_name_for_zone.get(z_name);
            let c_name = self.cool_system_name_for_zone.get(z_name);

            // Look up convective fraction for heating/cooling for this zone
            let (frac_convective_heat, temp_setpnt_heat) = match h_name {
                Some(h_name) => {
                    let system = &self.space_heat_systems[h_name].lock();
                    (
                        system.frac_convective(),
                        system.temp_setpnt(simulation_time_iteration),
                    )
                }
                None => (1.0, None),
            };
            let (frac_convective_cool, temp_setpnt_cool) = match c_name {
                Some(c_name) => {
                    let system = &self.space_cool_systems[c_name];
                    (
                        system.frac_convective(),
                        system.temp_setpnt(&simulation_time_iteration),
                    )
                }
                None => (1.0, None),
            };

            // Use default setpoints when there is no heat/cool system or
            // there is no setpoint for the current timestep

            // set heating setpoint to absolute zero to ensure no heating demand
            let temp_setpnt_heat = temp_setpnt_heat.unwrap_or(kelvin_to_celsius(0.));
            // set cooling setpoint to Planck temperature to ensure no cooling demand
            let temp_setpnt_cool = temp_setpnt_cool.unwrap_or(kelvin_to_celsius(1.4e32));

            let (
                space_heat_demand_zone_current,
                space_cool_demand_zone_current,
                h_ve_cool_extra_zone_current,
                _,
            ) = zone
                .space_heat_cool_demand(
                    delta_t_h,
                    temp_ext_air,
                    gains_internal_zone[z_name.as_str()],
                    gains_solar_zone[z_name.as_str()],
                    frac_convective_heat,
                    frac_convective_cool,
                    temp_setpnt_heat,
                    temp_setpnt_cool,
                    0.0,
                    None, // temporary nones while migrating to 0.30
                    None,
                    AirChangesPerHourArgument::from_ach_target_windows_open(0.0, 0.0), // temporary arg while migrating to 0.30
                    simulation_time_iteration,
                )
                .unwrap();

            space_heat_demand_zone.insert((*z_name).clone(), space_heat_demand_zone_current);
            space_cool_demand_zone.insert((*z_name).clone(), space_cool_demand_zone_current);
            h_ve_cool_extra_zone.insert((*z_name).clone(), h_ve_cool_extra_zone_current);

            if let Some(h_name) = h_name {
                *space_heat_demand_system.get_mut(h_name.as_str()).unwrap() +=
                    space_heat_demand_zone[z_name];
            }
            if let Some(c_name) = c_name {
                *space_cool_demand_system.get_mut(c_name.as_str()).unwrap() +=
                    space_cool_demand_zone[z_name];
            }
        }

        (
            space_heat_demand_system,
            space_cool_demand_system,
            space_heat_demand_zone,
            space_cool_demand_zone,
            h_ve_cool_extra_zone,
        )
    }
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

type NumberMap = HashMap<String, f64>;

pub type ResultsEndUser = IndexMap<KeyString, IndexMap<String, Vec<f64>>>;

fn has_unique_values<K, V: Eq + Hash>(map: &IndexMap<K, V>) -> bool {
    let values: Vec<&V> = map.values().collect();
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
        input.latitude.unwrap(),
        input.longitude.unwrap(),
        input.timezone.unwrap(),
        input.start_day.unwrap_or(0),
        input.end_day,
        input.time_series_step.unwrap_or(1.0),
        input.january_first,
        input.daylight_savings.clone().unwrap(),
        input.leap_day_included.unwrap_or(false),
        input.direct_beam_conversion_needed.unwrap_or(false),
        input.shading_segments.clone(),
    )
}

fn infiltration_from_input(input: &Infiltration) -> VentilationElementInfiltration {
    VentilationElementInfiltration::new(
        input.storeys_in_building,
        input.shelter,
        input.build_type,
        input.test_result,
        input.test_type,
        input.env_area,
        input.volume,
        input.sheltered_sides,
        input.open_chimneys,
        input.open_flues,
        input.closed_fire,
        input.flues_d,
        input.flues_e,
        input.blocked_chimneys,
        input.extract_fans,
        input.passive_vents,
        input.gas_fires,
        input.storey_of_dwelling,
    )
}

pub struct ColdWaterSources {
    mains_water: Option<ColdWaterSource>,
    header_tank: Option<ColdWaterSource>,
}

impl ColdWaterSources {
    #[cfg(test)]
    pub fn new(mains_water: Option<ColdWaterSource>, header_tank: Option<ColdWaterSource>) -> Self {
        Self {
            mains_water,
            header_tank,
        }
    }

    pub fn ref_for_mains_water(&self) -> Option<ColdWaterSource> {
        self.mains_water.clone()
    }

    pub fn ref_for_header_tank(&self) -> Option<ColdWaterSource> {
        self.header_tank.clone()
    }

    pub fn ref_for_type(&self, source_type: ColdWaterSourceType) -> Option<ColdWaterSource> {
        match source_type {
            ColdWaterSourceType::MainsWater => self.ref_for_mains_water(),
            ColdWaterSourceType::HeaderTank => self.ref_for_header_tank(),
        }
    }
}

fn cold_water_sources_from_input(
    input: &ColdWaterSourceInput,
    simulation_time: &SimulationTime,
) -> ColdWaterSources {
    ColdWaterSources {
        mains_water: input
            .mains_water
            .as_ref()
            .map(|details| cold_water_source_from_input_details(details, simulation_time)),
        header_tank: input
            .header_tank
            .as_ref()
            .map(|details| cold_water_source_from_input_details(details, simulation_time)),
    }
}

fn cold_water_source_from_input_details(
    details: &ColdWaterSourceDetails,
    simulation_time: &SimulationTime,
) -> ColdWaterSource {
    ColdWaterSource::new(
        details.temperatures.clone(),
        simulation_time,
        details.time_series_step,
    )
}

fn energy_supplies_from_input(
    input: &EnergySupplyInput,
    simulation_time_iterator: &SimulationTimeIterator,
    external_conditions: Arc<ExternalConditions>,
) -> EnergySupplies {
    EnergySupplies {
        mains_electricity: energy_supply_from_input(
            input.get(&EnergySupplyKey::MainsElectricity),
            simulation_time_iterator,
            external_conditions.clone(),
        ),
        mains_gas: energy_supply_from_input(
            input.get(&EnergySupplyKey::MainsGas),
            simulation_time_iterator,
            external_conditions.clone(),
        ),
        bulk_lpg: energy_supply_from_input(
            input.get(&EnergySupplyKey::BulkLpg),
            simulation_time_iterator,
            external_conditions.clone(),
        ),
        heat_network: energy_supply_from_input(
            input.get(&EnergySupplyKey::HeatNetwork),
            simulation_time_iterator,
            external_conditions,
        ),
        unmet_demand: Arc::new(RwLock::new(EnergySupply::new(
            FuelType::UnmetDemand,
            simulation_time_iterator.total_steps(),
            Default::default(),
            Default::default(),
            Default::default(),
        ))),
        bottled_lpg: None,
        condition_11f_lpg: None,
        custom: None,
    }
}

fn energy_supply_from_input(
    input: Option<&EnergySupplyDetails>,
    simulation_time_iterator: &SimulationTimeIterator,
    external_conditions: Arc<ExternalConditions>,
) -> Option<Arc<RwLock<EnergySupply>>> {
    input.map(|details| {
        Arc::new(RwLock::new(EnergySupply::new(
            details.fuel,
            simulation_time_iterator.total_steps(),
            details.electric_battery.as_ref().map(|battery_input| {
                ElectricBattery::from_input(
                    battery_input,
                    simulation_time_iterator.step_in_hours(),
                    external_conditions,
                )
            }),
            details.priority.as_ref().cloned(),
            details.is_export_capable,
        )))
    })
}

struct DiverterTypes {
    pub mains_electricity: Option<EnergyDiverter>,
    pub mains_gas: Option<EnergyDiverter>,
    pub bulk_lpg: Option<EnergyDiverter>,
}

impl From<&EnergySupplyInput> for DiverterTypes {
    fn from(input: &EnergySupplyInput) -> Self {
        Self {
            mains_electricity: diverter_from_energy_supply(
                input.get(&EnergySupplyKey::MainsElectricity),
            ),
            mains_gas: diverter_from_energy_supply(input.get(&EnergySupplyKey::MainsGas)),
            bulk_lpg: diverter_from_energy_supply(input.get(&EnergySupplyKey::BulkLpg)),
        }
    }
}

impl DiverterTypes {
    pub fn get_for_supply_type(
        &self,
        energy_supply_type: EnergySupplyType,
    ) -> Option<&EnergyDiverter> {
        match energy_supply_type {
            EnergySupplyType::Electricity => self.mains_electricity.as_ref(),
            EnergySupplyType::MainsGas => self.mains_gas.as_ref(),
            EnergySupplyType::LpgBulk => self.bulk_lpg.as_ref(),
            _ => unimplemented!(),
        }
    }
}

fn diverter_from_energy_supply(supply: Option<&EnergySupplyDetails>) -> Option<EnergyDiverter> {
    supply.and_then(|supply| supply.diverter.clone())
}

// #[derive(Default)]
// pub struct InternalGainsCollection {
//     total_internal_gains: Option<InternalGains>,
//     metabolic_gains: Option<InternalGains>,
//     _evaporative_losses: Option<InternalGains>,
//     lighting: Option<ApplianceGains>,
//     cooking: Option<ApplianceGains>,
//     cooking1: Option<ApplianceGains>,
//     cooking2: Option<ApplianceGains>,
//     other: Option<InternalGains>,
// }

pub type InternalGainsCollection = IndexMap<String, Gains>;

fn internal_gains_from_input(
    input: &InternalGainsInput,
    total_floor_area: f64,
) -> InternalGainsCollection {
    let mut gains_collection = InternalGainsCollection::from([]);
    input
        .total_internal_gains
        .as_ref()
        .and_then(|internal_gains| {
            gains_collection.insert(
                "total_internal_gains".to_string(),
                Gains::Internal(internal_gains_from_details(
                    internal_gains,
                    total_floor_area,
                )),
            )
        });
    input.metabolic_gains.as_ref().and_then(|internal_gains| {
        gains_collection.insert(
            "metabolic_gains".to_string(),
            Gains::Internal(internal_gains_from_details(
                internal_gains,
                total_floor_area,
            )),
        )
    });
    input
        .evaporative_losses
        .as_ref()
        .and_then(|internal_gains| {
            gains_collection.insert(
                "evaporative_losses".to_string(),
                Gains::Internal(internal_gains_from_details(
                    internal_gains,
                    total_floor_area,
                )),
            )
        });
    input.other.as_ref().and_then(|internal_gains| {
        gains_collection.insert(
            "other".to_string(),
            Gains::Internal(internal_gains_from_details(
                internal_gains,
                total_floor_area,
            )),
        )
    });

    gains_collection
}

fn internal_gains_from_details(
    details: &InternalGainsDetails,
    total_floor_area: f64,
) -> InternalGains {
    InternalGains::new(
        convert_energy_to_wm2(details, total_floor_area),
        details.start_day,
        details.time_series_step,
    )
}

fn convert_energy_to_wm2(
    internal_gains_details: &InternalGainsDetails,
    total_floor_area: f64,
) -> Vec<f64> {
    let schedule: IndexMap<String, Value> = (&internal_gains_details.schedule).into();
    expand_numeric_schedule(&schedule, false)
        .iter()
        .map(|energy_data| energy_data.unwrap_or_default() / total_floor_area)
        .collect()
}

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
) -> Controls {
    let mut core: Vec<HeatSourceControl> = Default::default();
    let mut extra: HashMap<String, Arc<Control>> = Default::default();

    // this is very ugly(!) but is just a reflection of the lack of clarity in the schema
    // and the way the variants-struct crate works;
    // we should be able to improve it in time
    for control in &control_input.core {
        match control {
            HeatSourceControlInput::HotWaterTimer(control) => {
                core.push(HeatSourceControl::HotWaterTimer(Arc::new(
                    single_control_from_details(control, simulation_time_iterator),
                )));
            }
            HeatSourceControlInput::WindowOpening(control) => {
                core.push(HeatSourceControl::WindowOpening(Arc::new(
                    single_control_from_details(control, simulation_time_iterator),
                )));
            }
            unknown => panic!(
                "incorrectly formed HeatSourceControl struct encountered: {:?}",
                unknown
            ),
        }
    }
    for (name, control) in &control_input.extra {
        extra.insert(
            name.to_string(),
            Arc::new(single_control_from_details(
                control,
                simulation_time_iterator,
            )),
        );
    }

    Controls { core, extra }
}

fn single_control_from_details(
    details: &ControlDetails,
    simulation_time_iterator: &SimulationTimeIterator,
) -> Control {
    match details {
        ControlDetails::OnOffTime {
            start_day,
            time_series_step,
            schedule,
            ..
        } => Control::OnOffTimeControl(OnOffTimeControl::new(
            expand_boolean_schedule(schedule, false),
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
            expand_numeric_schedule(schedule, false)
                .into_iter()
                .flatten()
                .collect(),
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
                expand_numeric_schedule(schedule, true).to_vec(),
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
        ControlDetails::ToUCharge {
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
                    Value::Array(charge_vec) => {
                        charge_level_vec = charge_vec.iter().map(|v| v.as_f64().unwrap()).collect();
                    }
                    Value::Number(charge) => {
                        charge_level_vec = vec![charge.as_f64().unwrap(); vec_size];
                    }
                    _ => {
                        panic!("Control charge value must be either a number of a list of numbers")
                    }
                }
            }

            Control::ToUChargeControl(ToUChargeControl {
                schedule: expand_boolean_schedule(schedule, false),
                start_day: *start_day,
                time_series_step: *time_series_step,
                charge_level: charge_level_vec,
            })
        }
    }
}

fn wwhrs_from_input(
    wwhrs: Option<&WasteWaterHeatRecovery>,
    cold_water_sources: &ColdWaterSources,
) -> IndexMap<String, Arc<Mutex<Wwhrs>>> {
    let mut wwhr_systems: IndexMap<String, Arc<Mutex<Wwhrs>>> = IndexMap::from([]);
    if let Some(systems) = wwhrs {
        for (name, system) in systems {
            wwhr_systems
                .entry(name.clone())
                .or_insert(Arc::new(Mutex::new(wwhr_system_from_details(
                    system.clone(),
                    cold_water_sources,
                ))));
        }
    }

    wwhr_systems
}

fn wwhr_system_from_details(
    system: WasteWaterHeatRecoveryDetails,
    cold_water_sources: &ColdWaterSources,
) -> Wwhrs {
    match system.system_type {
        WwhrsType::SystemA => Wwhrs::WWHRSInstantaneousSystemA(WWHRSInstantaneousSystemA::new(
            system.flow_rates,
            system.efficiencies,
            get_cold_water_source_ref_for_type(system.cold_water_source, cold_water_sources)
                .unwrap(),
            system.utilisation_factor,
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
        )),
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

pub type KeyString = ArrayString<64>;

pub type RunResults = (
    Vec<f64>,
    IndexMap<KeyString, Vec<f64>>,
    ResultsEndUser,
    IndexMap<KeyString, Vec<f64>>,
    IndexMap<KeyString, Vec<f64>>,
    IndexMap<KeyString, Vec<f64>>,
    IndexMap<KeyString, Vec<f64>>,
    IndexMap<KeyString, Vec<f64>>,
    IndexMap<KeyString, Vec<f64>>,
    IndexMap<KeyString, Vec<f64>>,
    IndexMap<KeyString, IndexMap<KeyString, Vec<f64>>>,
    Vec<KeyString>,
    IndexMap<KeyString, IndexMap<KeyString, Vec<f64>>>,
    IndexMap<KeyString, HotWaterResultMap>,
    IndexMap<KeyString, NumberOrDivisionByZero>,
    IndexMap<KeyString, NumberOrDivisionByZero>,
    IndexMap<KeyString, NumberOrDivisionByZero>,
    IndexMap<KeyString, Vec<f64>>,
    IndexMap<KeyString, IndexMap<KeyString, IndexMap<KeyString, f64>>>,
    IndexMap<KeyString, f64>,
    IndexMap<KeyString, f64>,
);

type SpaceHeatingCalculation<'a> = (
    HashMap<&'a str, f64>,
    HashMap<&'a str, f64>,
    HashMap<&'a str, f64>,
    HashMap<&'a str, f64>,
    HashMap<String, f64>,
    HashMap<String, f64>,
    HashMap<String, f64>,
    HashMap<String, f64>,
    HashMap<&'a str, f64>,
    HashMap<&'a str, f64>,
    f64,
    HashMap<&'a str, Option<HeatBalance>>,
);

fn get_cold_water_source_ref_for_type(
    source_type: ColdWaterSourceType,
    cold_water_sources: &ColdWaterSources,
) -> Option<ColdWaterSource> {
    match source_type {
        ColdWaterSourceType::MainsWater => cold_water_sources.ref_for_mains_water(),
        ColdWaterSourceType::HeaderTank => cold_water_sources.ref_for_header_tank(),
    }
}

pub type EventSchedule = Vec<Option<Vec<TypedScheduleEvent>>>;

#[derive(Clone)]
pub struct HotWaterEventSchedules {
    pub shower: HashMap<String, EventSchedule>,
    pub bath: HashMap<String, EventSchedule>,
    pub other: HashMap<String, EventSchedule>,
}

fn event_schedules_from_input(
    events: &WaterHeatingEvents,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<HotWaterEventSchedules> {
    let mut shower_schedules: HashMap<String, EventSchedule> = Default::default();
    let shower_events = &events.shower;
    for (name, events) in shower_events {
        shower_schedules.insert(
            name.to_owned(),
            schedule_event_from_input(
                events.iter().collect(),
                name,
                WaterScheduleEventType::Shower,
                None,
                simulation_time_iterator,
            )?,
        );
    }

    let mut bath_schedules: HashMap<String, EventSchedule> = Default::default();
    let bath_events = &events.bath;
    for (name, events) in bath_events {
        bath_schedules.insert(
            name.to_owned(),
            schedule_event_from_input(
                events.iter().collect(),
                name,
                WaterScheduleEventType::Bath,
                None,
                simulation_time_iterator,
            )?,
        );
    }

    let mut other_schedules: HashMap<String, EventSchedule> = Default::default();
    let other_events = &events.other;
    for (name, events) in other_events {
        other_schedules.insert(
            name.to_owned(),
            schedule_event_from_input(
                events.iter().collect(),
                name,
                WaterScheduleEventType::Other,
                None,
                simulation_time_iterator,
            )?,
        );
    }

    Ok(HotWaterEventSchedules {
        shower: shower_schedules,
        bath: bath_schedules,
        other: other_schedules,
    })
}

fn schedule_event_from_input(
    events_input: Vec<&WaterHeatingEvent>,
    name: &str,
    event_type: WaterScheduleEventType,
    existing_schedule: Option<Vec<Option<Vec<TypedScheduleEvent>>>>,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<EventSchedule> {
    let sim_timestep = simulation_time_iterator.step_in_hours();
    let total_timesteps = simulation_time_iterator.total_steps();

    let schedule = if let Some(existing_schedule) = existing_schedule {
        existing_schedule
    } else {
        vec![None; total_timesteps]
    };

    expand_events(
        events_input
            .iter()
            .map(|event| ScheduleEvent::from(*event))
            .collect::<Vec<_>>(),
        sim_timestep,
        total_timesteps,
        name,
        event_type,
        schedule,
    )
}

fn ductwork_from_ventilation_input(_ventilation: &Option<Ventilation>) -> Option<Ductwork> {
    None
}

fn ventilation_from_input<'a>(
    energy_supply_name: &str,
    ventilation: &'a Ventilation,
    infiltration: &VentilationElementInfiltration,
    simulation_time: &'a SimulationTimeIterator,
    energy_supplies: &mut EnergySupplies,
) -> anyhow::Result<VentilationElement> {
    Ok(match ventilation {
        Ventilation::Whev {
            req_ach,
            sfp,
            energy_supply,
        } => {
            let energy_supply = energy_supplies
                .ensured_get_for_type(*energy_supply, simulation_time.total_steps())?;
            let energy_supply_conn =
                EnergySupply::connection(energy_supply, energy_supply_name).unwrap();

            VentilationElement::Whev(WholeHouseExtractVentilation::new(
                *req_ach,
                *sfp,
                infiltration.infiltration_rate(),
                energy_supply_conn,
                // energy_supply_from_type_for_ventilation(energy_supply, energy_supplies),
                // "Ventilation system".to_string(),
                simulation_time.step_in_hours(),
            ))
        }
        Ventilation::Natural { req_ach } => VentilationElement::Natural(NaturalVentilation::new(
            *req_ach,
            infiltration.infiltration_rate(),
        )),
    })
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

fn zone_from_input<'a>(
    input: &ZoneInput,
    external_conditions: Arc<ExternalConditions>,
    infiltration_ventilation: Arc<InfiltrationVentilation>,
    window_adjust_control: Option<Arc<Control>>,
    simulation_time_iterator: &SimulationTimeIterator,
) -> anyhow::Result<(Zone, Option<String>, Option<String>)> {
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
            simulation_time_iterator,
        )?,
        heat_system_name,
        cool_system_name,
    ))
}

fn infiltration_ventilation_from_input(
    zones: &ZoneDictionary,
    input: &InfiltrationVentilationInput,
    controls: &Controls,
    energy_supplies: &mut EnergySupplies,
    simulation_time: &SimulationTimeIterator,
    total_volume: f64,
    external_conditions: Arc<ExternalConditions>,
) -> anyhow::Result<(
    InfiltrationVentilation,
    Option<Arc<Control>>,
    IndexMap<String, MechanicalVentilation>,
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
                                // running orientation through init_orientation again to convert it back to "orientation360" value, which is expected by the Python (possibly erroneously?)
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
    let average_pitch = if pitches.len() > 0 {
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

    let mut mechanical_ventilations: IndexMap<String, MechanicalVentilation> = Default::default();
    let mut space_heating_ductwork: IndexMap<String, Vec<Ductwork>> = Default::default();

    if let Some(mech_vent_input) = input.mechanical_ventilation.as_ref() {
        for (mech_vents_name, mech_vents_data) in mech_vent_input {
            let ctrl_intermittent_mev = mech_vents_data
                .control
                .as_ref()
                .map(|ctrl_name| controls.get_with_string(ctrl_name));

            let energy_supply = energy_supplies.ensured_get_for_type(
                mech_vents_data.energy_supply,
                simulation_time.total_steps(),
            )?;
            let energy_supply_connection =
                EnergySupply::connection(energy_supply.clone(), mech_vents_name)?;

            mechanical_ventilations.insert(
                mech_vents_name.clone(),
                MechanicalVentilation::new(external_conditions.clone(), mech_vents_data.supply_air_flow_rate_control, mech_vents_data.supply_air_temperature_control_type, 0., 0., mech_vents_data.vent_type, mech_vents_data.sfp.ok_or_else(|| anyhow!("A specific fan power value is expected for a mechanical ventilation unit."))?, mech_vents_data.design_outdoor_air_flow_rate, energy_supply_connection, total_volume, *altitude, None, match mech_vents_data.vent_type {
                    VentType::Mvhr => mech_vents_data.mvhr_efficiency,
                    VentType::IntermittentMev
                    | VentType::CentralisedContinuousMev
                    | VentType::DecentralisedContinuousMev => {
                        None
                    }
                    VentType::Piv => bail!("PIV vent type is not currently recognised when building up mechanical ventilation values for calculation"),
                }, None),
            );

            // TODO (from Python) not all dwellings have mech vents - update to make mech vents optional
            if mech_vents_data.vent_type == VentType::Mvhr {
                // the Python nixes the ductwork map here, perhaps erroneously?
                space_heating_ductwork = Default::default();
                space_heating_ductwork.insert(
                    mech_vents_name.to_owned(),
                    mech_vents_data
                        .ductwork
                        .as_ref()
                        .iter()
                        .map(|ductworks| {
                            ductworks.iter().map(|ductwork| -> anyhow::Result<Ductwork> {
                                let (duct_perimeter, internal_diameter, external_diameter) =
                                    match ductwork.cross_section_shape {
                                        DuctShape::Circular => (None, Some(ductwork.internal_diameter_mm.ok_or_else(|| anyhow!("Expected an internal diameter value for ductwork with a circular cross-section."))? / MILLIMETRES_IN_METRE as f64), Some(ductwork.external_diameter_mm.ok_or_else(|| anyhow!("Expected an internal diameter value for ductwork with a circular cross-section."))? / MILLIMETRES_IN_METRE as f64)),
                                        DuctShape::Rectangular => (Some(ductwork.duct_perimeter_mm.ok_or_else(|| anyhow!("Expected a duct perimeter value for ductwork with a rectangular cross-section."))?), None, None),
                                    };

                                Ductwork::new(ductwork.cross_section_shape, duct_perimeter, internal_diameter, external_diameter, ductwork.length, ductwork.insulation_thermal_conductivity, ductwork.insulation_thickness_mm, ductwork.reflective, ductwork.duct_type, mech_vents_data.mvhr_location.ok_or_else(|| anyhow!("An MVHR location was expected for mechanical ventilation with an MVHR vent type."))?, mech_vents_data.mvhr_efficiency.ok_or_else(|| anyhow!("An MVHR efficiency value was expected for mechanical ventilation with an MVHR vent type."))?)
                            })
                        })
                        .flatten()
                        .collect::<anyhow::Result<Vec<Ductwork>>>()?,
                );
            }
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
        vec![],
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
    energy_supplies: &mut EnergySupplies,
    total_floor_area: f64,
    simulation_time: &SimulationTimeIterator,
) -> anyhow::Result<()> {
    for (name, gains_details) in input {
        let energy_supply_conn = EnergySupply::connection(
            energy_supplies
                .ensured_get_for_type(gains_details.energy_supply, simulation_time.total_steps())?,
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
            ))
        };

        internal_gains_collection.insert(name.clone(), gains);
    }

    Ok(())
}

fn appliance_gains_from_single_input(
    input: &ApplianceGainsDetails,
    energy_supply_connection: EnergySupplyConnection,
    total_floor_area: f64,
) -> ApplianceGains {
    let total_energy_supply = expand_numeric_schedule(&input.schedule, false)
        .iter()
        .map(|energy_data| energy_data.unwrap() / total_floor_area)
        .collect();

    ApplianceGains::new(
        total_energy_supply,
        input.gains_fraction,
        input.start_day,
        input.time_series_step,
        energy_supply_connection,
    )
}

#[derive(Clone, Debug)]
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
            HeatSource::Wet(ref mut wet) => match wet.as_mut() {
                HeatSourceWet::WaterCombi(_) => {
                    unimplemented!("not expected? this value does not have a demand_energy method")
                    // the Python uses duck-typing here but there is no method for this type
                }
                HeatSourceWet::WaterRegular(ref mut r) => Ok(r
                    .demand_energy(
                        energy_demand,
                        temp_return,
                        None,
                        None,
                        simulation_time_iteration,
                    )
                    .expect("Regular water boiler could not register energy demand")
                    .0),
                HeatSourceWet::Space(_) => {
                    unimplemented!("not expected? this value does not have a demand_energy method")
                    // the Python uses duck-typing here but there is no method for this type
                }
                HeatSourceWet::HeatNetworkWaterStorage(ref mut h) => {
                    Ok(h.demand_energy(energy_demand, temp_return, &simulation_time_iteration))
                }
                HeatSourceWet::HeatBatteryHotWater(ref mut h) => {
                    Ok(h.demand_energy(energy_demand, temp_return, simulation_time_iteration))
                }
                HeatSourceWet::HeatPumpWater(ref mut h) => {
                    h.demand_energy(energy_demand, temp_return, simulation_time_iteration)
                }
                HeatSourceWet::HeatPumpWaterOnly(h) => {
                    Ok(h.demand_energy(energy_demand, temp_return, simulation_time_iteration))
                }
            },
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

#[derive(Clone)]
pub enum WetHeatSource {
    HeatPump(Arc<Mutex<HeatPump>>),
    Boiler(Boiler),
    Hiu(HeatNetwork),
    HeatBattery(HeatBattery),
}

impl WetHeatSource {
    pub fn timestep_end(&mut self, simtime: SimulationTimeIteration) {
        match self {
            WetHeatSource::HeatPump(heat_pump) => heat_pump.lock().timestep_end(simtime.index),
            WetHeatSource::Boiler(boiler) => boiler.timestep_end(simtime),
            WetHeatSource::Hiu(heat_network) => heat_network.timestep_end(simtime.index),
            WetHeatSource::HeatBattery(heat_battery) => heat_battery.timestep_end(simtime.index),
        }
    }
}

fn heat_source_wet_from_input(
    energy_supply_name: &str,
    input: HeatSourceWetDetails,
    external_conditions: Arc<ExternalConditions>,
    simulation_time: Arc<SimulationTimeIterator>,
    ventilation: Option<VentilationElement>,
    ventilation_req_ach: Option<f64>,
    mechanical_ventilations: &IndexMap<String, Arc<MechanicalVentilation>>,
    number_of_zones: usize,
    temp_internal_air_accessor: TempInternalAirAccessor,
    total_volume: f64,
    controls: &Controls,
    energy_supplies: &mut EnergySupplies,
) -> anyhow::Result<WetHeatSource> {
    match &input {
        HeatSourceWetDetails::HeatPump {
            source_type,
            energy_supply_heat_network,
            energy_supply,
            ..
        } => {
            let throughput_exhaust_air = if source_type.is_exhaust_air() {
                // Check that ventilation system is compatible with exhaust air HP
                // the Python code here will assign to throughput_exhaust_air on the last iteration, though mech vents collection is not ordered - this may be erroneous
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
                debug_assert_eq!(
                    energy_supply_heat_network,
                    &Some("heat network".to_string()),
                    "value for EnergySupply_heat_network is always expected to be 'heat network'"
                );
                energy_supplies.heat_network.clone()
            } else {
                None
            };

            let energy_supply = energy_supplies
                .ensured_get_for_type(*energy_supply, simulation_time.total_steps())?;
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
                    DETAILED_OUTPUT_HEATING_COOLING,
                    None,
                    None, // temporary during migrating to 0.30
                    temp_internal_air_accessor,
                )?,
            ))))
        }
        HeatSourceWetDetails::Boiler {
            energy_supply,
            energy_supply_auxiliary,
            ..
        } => {
            let energy_supply = energy_supplies
                .ensured_get_for_type(*energy_supply, simulation_time.total_steps())?;
            let energy_supply_aux = energy_supplies
                .ensured_get_for_type(*energy_supply_auxiliary, simulation_time.total_steps())?;
            let aux_supply_name = format!("Boiler_auxiliary: {energy_supply_name}");
            let energy_supply_conn_aux =
                EnergySupply::connection(energy_supply_aux.clone(), aux_supply_name.as_str())?;

            Ok(WetHeatSource::Boiler(
                Boiler::new(
                    input,
                    energy_supply,
                    energy_supply_conn_aux,
                    external_conditions.clone(),
                    simulation_time.step_in_hours(),
                )
                .expect("could not construct boiler value from provided data"),
            ))
        }
        HeatSourceWetDetails::Hiu {
            power_max,
            hiu_daily_loss,
            building_level_distribution_losses,
            energy_supply,
            ..
        } => {
            let energy_supply = energy_supplies
                .ensured_get_for_type(*energy_supply, simulation_time.total_steps())?;
            let energy_supply_conn_name_auxiliary =
                format!("HeatNetwork_auxiliary: {energy_supply_name}");
            let energy_supply_conn_name_building_level_distribution_losses =
                format!("HeatNetwork_building_level_distribution_losses: {energy_supply_name}");

            Ok(WetHeatSource::Hiu(HeatNetwork::new(
                *power_max,
                *hiu_daily_loss,
                *building_level_distribution_losses,
                energy_supply,
                energy_supply_conn_name_auxiliary,
                energy_supply_conn_name_building_level_distribution_losses,
                simulation_time.step_in_hours(),
            )))
            // TODO add heat network to timestep_end_calcs
        }
        HeatSourceWetDetails::HeatBattery {
            control_charge,
            energy_supply,
            ..
        } => {
            let energy_supply = energy_supplies
                .ensured_get_for_type(*energy_supply, simulation_time.total_steps())?;
            let energy_supply_conn =
                EnergySupply::connection(energy_supply.clone(), energy_supply_name).unwrap();

            let heat_source = WetHeatSource::HeatBattery(HeatBattery::new(
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
            ));
            Ok(heat_source)
        }
    }
}

fn heat_source_from_input(
    name: &str,
    input: &HeatSourceInput,
    temp_setpoint: f64,
    volume: f64,
    daily_losses: f64,
    heat_exchanger_surface_area: Option<f64>,
    wet_heat_sources: &IndexMap<String, Arc<Mutex<WetHeatSource>>>,
    simulation_time: &SimulationTimeIterator,
    controls: &Controls,
    energy_supplies: &mut EnergySupplies,
    cold_water_sources: &ColdWaterSources,
    temp_internal_air_accessor: TempInternalAirAccessor,
    external_conditions: Arc<ExternalConditions>,
) -> anyhow::Result<(HeatSource, String)> {
    // TODO add in all the stuff to do with energy supply

    match input {
        HeatSourceInput::ImmersionHeater {
            power,
            control,
            energy_supply,
            ..
        } => {
            let energy_supply = energy_supplies
                .ensured_get_for_type(*energy_supply, simulation_time.total_steps())?;
            let energy_supply_conn = EnergySupply::connection(energy_supply.clone(), name).unwrap();

            Ok((
                HeatSource::Storage(HeatSourceWithStorageTank::Immersion(Arc::new(Mutex::new(
                    ImmersionHeater::new(
                        *power,
                        energy_supply_conn,
                        simulation_time.step_in_hours(),
                        (*control).and_then(|ctrl| controls.get(&ctrl)),
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
            let energy_supply = energy_supplies
                .ensured_get_for_type(*energy_supply, simulation_time.total_steps())?;
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
                        temp_internal_air_accessor,
                        simulation_time.step_in_hours(),
                        *WATER,
                    ),
                )))),
                name.into(),
            ))
        }
        HeatSourceInput::Wet {
            name,
            cold_water_source: cold_water_source_type,
            control,
            temp_flow_limit_upper,
            ..
        } => {
            let cold_water_source = cold_water_sources
                .ref_for_type(
                    cold_water_source_type
                        .expect("Expect a cold water source to be defined on a wet heat source"),
                )
                .expect("Expected a cold water source to be available to a boiler heat source.");
            let energy_supply_conn_name = format!("{name}_water_heating");
            let heat_source_wet = wet_heat_sources
                .get(name)
                .unwrap_or_else(|| {
                    panic!("Expected a wet heat source registered with the name '{name}'.")
                })
                .clone();
            let source_control = (*control).and_then(|ctrl| controls.get(&ctrl));

            let lock = heat_source_wet.lock();
            let mut heat_source_wet_clone = (*lock).clone();

            Ok((
                match heat_source_wet_clone {
                    WetHeatSource::HeatPump(heat_pump) => HeatSource::Wet(Box::new(
                        HeatSourceWet::HeatPumpWater(HeatPump::create_service_hot_water(
                            heat_pump.clone(),
                            energy_supply_conn_name.clone(),
                            55.,
                            temp_flow_limit_upper
                                .expect("temp_flow_limit_upper field was expected to be set"),
                            Arc::new(cold_water_source),
                            source_control,
                        )),
                    )),
                    WetHeatSource::Boiler(ref mut boiler) => HeatSource::Wet(Box::new(
                        HeatSourceWet::WaterRegular(boiler.create_service_hot_water_regular(
                            energy_supply_conn_name.clone(),
                            temp_setpoint,
                            source_control,
                        )),
                    )),
                    WetHeatSource::Hiu(heat_network) => {
                        HeatSource::Wet(Box::new(HeatSourceWet::HeatNetworkWaterStorage(
                            HeatNetwork::create_service_hot_water_storage(
                                Arc::new(Mutex::new(heat_network)),
                                energy_supply_conn_name.clone(),
                                temp_setpoint,
                                source_control,
                            ),
                        )))
                    }
                    WetHeatSource::HeatBattery(battery) => {
                        HeatSource::Wet(Box::new(HeatSourceWet::HeatBatteryHotWater(
                            HeatBattery::create_service_hot_water_regular(
                                Arc::new(Mutex::new(battery)),
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
            let energy_supply = energy_supplies
                .ensured_get_for_type(*energy_supply, simulation_time.total_steps())?;
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
                        controls.get(control),
                    ),
                ))),
                energy_supply_conn_name.into(),
            ))
        }
    }
}

pub enum HotWaterSource {
    StorageTank(Arc<Mutex<StorageTank>>),
    CombiBoiler(BoilerServiceWaterCombi),
    PointOfUse(PointOfUse),
    HeatNetwork(HeatNetworkServiceWaterDirect),
    HeatBattery(()),
}

impl HotWaterSource {
    pub fn get_cold_water_source(&self) -> Option<WaterSourceWithTemperature> {
        match self {
            HotWaterSource::StorageTank(source) => {
                Some(source.lock().get_cold_water_source().clone())
            }
            HotWaterSource::CombiBoiler(source) => Some(source.get_cold_water_source().clone()),
            HotWaterSource::PointOfUse(source) => Some(source.get_cold_water_source().clone()),
            HotWaterSource::HeatNetwork(source) => Some(source.get_cold_water_source().clone()),
            HotWaterSource::HeatBattery(_) => None,
        }
    }

    pub fn demand_hot_water(
        &mut self,
        vol_demand_target: IndexMap<DemandVolTargetKey, VolumeReference>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        match self {
            HotWaterSource::StorageTank(_) => {
                // StorageTank does not match the same method signature or return type as as all other Hot Water sources
                panic!("demand_hot_water for HotWaterSource::StorageTank should be called directly on the HotWaterSource::StorageTank");
            }
            HotWaterSource::CombiBoiler(ref mut source) => source
                .demand_hot_water(vol_demand_target, simulation_time_iteration)
                .expect("Combi boiler could not calc demand hot water."),
            HotWaterSource::PointOfUse(ref mut source) => {
                source.demand_hot_water(vol_demand_target, &simulation_time_iteration)
            }
            HotWaterSource::HeatNetwork(ref mut source) => {
                source.demand_hot_water(vol_demand_target, simulation_time_iteration.index)
            }
            HotWaterSource::HeatBattery(_) => Default::default(),
        }
    }
}

fn hot_water_source_from_input(
    source_name: String,
    input: &HotWaterSourceDetails,
    cold_water_sources: &ColdWaterSources,
    wet_heat_sources: &IndexMap<String, Arc<Mutex<WetHeatSource>>>,
    wwhrs: &IndexMap<String, Arc<Mutex<Wwhrs>>>,
    controls: &Controls,
    energy_supplies: &mut EnergySupplies,
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
        } => {
            let mut cold_water_source: WaterSourceWithTemperature =
                cold_water_source_for_type(cold_water_source_type, cold_water_sources);
            if !wwhrs.is_empty() {
                for heat_recovery_system in wwhrs.values() {
                    cold_water_source =
                        WaterSourceWithTemperature::Wwhrs(heat_recovery_system.clone());
                }
            }
            let pipework = primary_pipework.as_ref().and_then(|p| p.into());
            let mut heat_sources: IndexMap<String, PositionedHeatSource> = Default::default();
            let mut heat_source_for_diverter: Option<Arc<Mutex<HeatSource>>> = Default::default();
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
                    *setpoint_temp,
                    *volume,
                    *daily_losses,
                    heat_exchanger_surface_area,
                    wet_heat_sources,
                    simulation_time,
                    controls,
                    energy_supplies,
                    cold_water_sources,
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
                heat_source_for_diverter = Some(heat_source);
                energy_supply_conn_names.push(energy_supply_conn_name);
            }
            let ctrl_hold_at_setpoint = control_hold_at_setpoint
                .as_ref()
                .and_then(|ctrl| controls.get_with_string(ctrl.as_str()));
            let storage_tank = Arc::new(Mutex::new(StorageTank::new(
                *volume,
                *daily_losses,
                *min_temp,
                *setpoint_temp,
                cold_water_source,
                simulation_time.step_in_hours(),
                heat_sources,
                temp_internal_air_accessor,
                external_conditions,
                Some(24),
                pipework, // TODO as part of migration 0.28 to 0.30: double check
                Some(
                    EnergySupply::connection(energy_supplies.unmet_demand.clone(), &source_name)
                        .unwrap(),
                ),
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
                        let immersion_heater = heat_source_for_diverter
                            .clone()
                            .expect("More than one heat source was expected to be present")
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
            ..
        } => {
            let cold_water_source =
                cold_water_source_for_type(cold_water_source_type, cold_water_sources);
            let energy_supply_conn_name = format!(
                "{}_water_heating",
                heat_source_wet_type.to_canonical_string()
            );
            energy_supply_conn_names.push(energy_supply_conn_name.clone());
            let mut heat_source_wet = match heat_source_wet_type {
                HeatSourceWetType::Boiler => {
                    match &*wet_heat_sources
                        .get("boiler")
                        .expect("Expected a boiler as wet heat source")
                        .lock()
                    {
                        WetHeatSource::Boiler(boiler) => boiler.clone(),
                        _ => panic!("Expected a boiler here"),
                    }
                }
                _ => panic!("Did not expect a heat source type that was not a boiler"),
            };
            HotWaterSource::CombiBoiler(
                heat_source_wet
                    .create_service_hot_water_combi(
                        cloned_input,
                        energy_supply_conn_name,
                        60.,
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
        } => {
            let energy_supply = energy_supplies
                .ensured_get_for_type(*energy_supply, simulation_time.total_steps())?;
            let energy_supply_conn_name = source_name;
            energy_supply_conn_names.push(energy_supply_conn_name.clone());
            let energy_supply_conn =
                EnergySupply::connection(energy_supply.clone(), &energy_supply_conn_name)?;
            let cold_water_source =
                cold_water_source_for_type(cold_water_source_type, cold_water_sources);
            HotWaterSource::PointOfUse(PointOfUse::new(
                *efficiency,
                energy_supply_conn,
                cold_water_source,
                *setpoint_temp,
            ))
        }
        HotWaterSourceDetails::Hiu {
            cold_water_source: cold_water_source_type,
            heat_source_wet: heat_source_wet_type,
            ..
        } => {
            let energy_supply_conn_name = format!(
                "{}_water_heating",
                heat_source_wet_type.to_canonical_string()
            );
            energy_supply_conn_names.push(energy_supply_conn_name.clone());
            let cold_water_source =
                cold_water_source_for_type(cold_water_source_type, cold_water_sources);
            let heat_source_wet = match heat_source_wet_type {
                HeatSourceWetType::HeatNetwork => {
                    match &*wet_heat_sources
                        .get("HeatNetwork")
                        .expect("expected a heat network in this context")
                        .lock()
                    {
                        WetHeatSource::Hiu(heat_network) => heat_network.clone(),
                        _ => panic!("expected a heat network in this context"),
                    }
                }
                _ => panic!("expected a heat network in this context"),
            };
            HotWaterSource::HeatNetwork(HeatNetwork::create_service_hot_water_direct(
                Arc::new(Mutex::new(heat_source_wet.clone())),
                energy_supply_conn_name,
                60.,
                cold_water_source,
            ))
        }
        HotWaterSourceDetails::HeatBattery { .. } => todo!(), // TODO is from Python
    };

    Ok((hot_water_source, energy_supply_conn_names))
}

fn cold_water_source_for_type(
    cold_water_source_type: &ColdWaterSourceType,
    cold_water_sources: &ColdWaterSources,
) -> WaterSourceWithTemperature {
    WaterSourceWithTemperature::ColdWaterSource(Arc::new(match cold_water_source_type {
        ColdWaterSourceType::MainsWater => cold_water_sources
            .mains_water
            .as_ref()
            .expect("referenced cold water source was expected to exist")
            .clone(),
        ColdWaterSourceType::HeaderTank => cold_water_sources
            .header_tank
            .as_ref()
            .expect("referenced cold water source was expected to exist")
            .clone(),
    }))
}

fn space_heat_systems_from_input(
    input: &SpaceHeatSystemInput,
    controls: &Controls,
    energy_supplies: &mut EnergySupplies,
    simulation_time: &SimulationTimeIterator,
    heat_sources_wet: &HashMap<String, Arc<WetHeatSource>>,
    heat_system_names_requiring_overvent: &mut Vec<String>,
    heat_system_name_for_zone: &IndexMap<String, String>,
    zones: &Arc<IndexMap<String, Zone>>,
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
                    SpaceHeatSystemDetails::ElectricStorageHeater { .. } => unimplemented!(), // requires implementation of ElecStorageHeater, make sure to add energy supply conn name to energy_conn_names_for_systems collection
                    SpaceHeatSystemDetails::WetDistribution { .. } => unimplemented!(), // requires implementation of Emitters, make sure to add energy supply conn name to energy_conn_names_for_systems collection
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
                        match heat_source.as_ref() {
                            WetHeatSource::HeatPump(heat_pump) => {
                                if heat_pump.lock().source_is_exhaust_air() {
                                    heat_system_names_requiring_overvent.push((*system_name).clone());
                                }
                                let volume_heated = total_volume_heated_by_system(zones, heat_system_name_for_zone, system_name);
                                SpaceHeatSystem::WarmAir(HeatPump::create_service_space_heating_warm_air(heat_pump.clone(), energy_supply_conn_name, control
                                    .as_ref()
                                    .and_then(|ctrl| controls.get_with_string(ctrl)).expect("A control object was expected for a heat pump warm air system"), *frac_convective, volume_heated).unwrap())
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
    energy_supplies: &mut EnergySupplies,
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
    energy_supplies: &mut EnergySupplies,
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
                    *inverter_peak_power,
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
                (system_name == heat_system_name).then_some(zone.volume())
            } else {
                None
            }
        })
        .sum::<f64>()
}

fn required_vent_data_from_input(input: &ControlInput) -> Option<RequiredVentData> {
    input
        .extra
        .get("required_vent")
        .map(|ctrl| RequiredVentData {
            schedule: expand_numeric_schedule(ctrl.schedule(), true),
            start_day: ctrl.start_day(),
            time_series_step: ctrl.time_series_step(),
        })
}

#[derive(Clone, Debug)]
struct RequiredVentData {
    schedule: NumericSchedule,
    start_day: u32,
    time_series_step: f64,
}
