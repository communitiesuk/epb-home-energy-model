use crate::compare_floats::max_of_2;
use crate::core::common::WaterSourceWithTemperature;
use crate::core::controls::time_control::{
    Control, OnOffMinimisingTimeControl, OnOffTimeControl, SetpointTimeControl, ToUChargeControl,
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
    expand_boolean_schedule, expand_numeric_schedule, expand_water_heating_events, ScheduleEvent,
};
use crate::core::space_heat_demand::building_element::area_for_building_element_input;
use crate::core::space_heat_demand::internal_gains::{ApplianceGains, Gains, InternalGains};
use crate::core::space_heat_demand::thermal_bridge::{ThermalBridge, ThermalBridging};
use crate::core::space_heat_demand::ventilation_element::{
    air_change_rate_to_flow_rate, MechanicalVentilationHeatRecovery, NaturalVentilation,
    VentilationElement, VentilationElementInfiltration, WholeHouseExtractVentilation,
    WindowOpeningForCooling,
};
use crate::core::space_heat_demand::zone::{NamedBuildingElement, Zone};
use crate::core::units::{
    kelvin_to_celsius, LITRES_PER_CUBIC_METRE, MILLIMETRES_IN_METRE, SECONDS_PER_HOUR,
    WATTS_PER_KILOWATT,
};
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::dhw_demand::DomesticHotWaterDemand;
use crate::core::water_heat_demand::misc::water_demand_to_kwh;
use crate::external_conditions::ExternalConditions;
use crate::input::{
    ApplianceGains as ApplianceGainsInput, ApplianceGainsDetails, BuildingElement,
    ColdWaterSourceDetails, ColdWaterSourceInput, ColdWaterSourceType, Control as ControlInput,
    ControlDetails, EnergyDiverter, EnergySupplyDetails, EnergySupplyInput, EnergySupplyType,
    ExternalConditionsInput, HeatNetwork as HeatNetworkInput, HeatPumpSourceType,
    HeatSource as HeatSourceInput, HeatSourceControl, HeatSourceControlType, HeatSourceWetDetails,
    HeatSourceWetType, HotWaterSourceDetails, Infiltration, Input,
    InternalGains as InternalGainsInput, InternalGainsDetails, OnSiteGeneration,
    OnSiteGenerationDetails, SpaceCoolSystem as SpaceCoolSystemInput, SpaceCoolSystemDetails,
    SpaceCoolSystemType, SpaceHeatSystem as SpaceHeatSystemInput, SpaceHeatSystemDetails,
    ThermalBridging as ThermalBridgingInput, ThermalBridgingDetails, Ventilation,
    WasteWaterHeatRecovery, WasteWaterHeatRecoveryDetails, WaterHeatingEvent, WaterHeatingEvents,
    WindowOpeningForCooling as WindowOpeningForCoolingInput, WwhrsType, ZoneDictionary, ZoneInput,
};
use crate::simulation_time::{SimulationTime, SimulationTimeIteration, SimulationTimeIterator};
use anyhow::bail;
use arrayvec::ArrayString;
use indexmap::IndexMap;
use indicatif::ProgressBar;
use parking_lot::Mutex;
use serde_json::Value;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::sync::Arc;

// TODO make this a runtime parameter?
const DETAILED_OUTPUT_HEATING_COOLING: bool = true;

pub struct Corpus {
    pub simulation_time: Arc<SimulationTimeIterator>,
    pub external_conditions: Arc<ExternalConditions>,
    pub infiltration: VentilationElementInfiltration,
    pub cold_water_sources: ColdWaterSources,
    pub energy_supplies: EnergySupplies,
    pub internal_gains: InternalGainsCollection,
    pub controls: Controls,
    pub wwhrs: HashMap<String, Wwhrs>,
    pub event_schedules: HotWaterEventSchedules,
    pub domestic_hot_water_demand: DomesticHotWaterDemand,
    pub ventilation: Option<Arc<Mutex<VentilationElement>>>,
    pub space_heating_ductwork: Option<Ductwork>,
    pub zones: HashMap<String, Zone>,
    pub energy_supply_conn_unmet_demand_zone: HashMap<String, Arc<EnergySupplyConnection>>,
    pub heat_system_name_for_zone: HashMap<String, Option<String>>,
    pub cool_system_name_for_zone: HashMap<String, Option<String>>,
    pub total_floor_area: f64,
    pub total_volume: f64,
    pub wet_heat_sources: HashMap<String, Arc<Mutex<WetHeatSource>>>,
    pub hot_water_sources: HashMap<String, HotWaterSource>,
    pub heat_system_names_requiring_overvent: Vec<String>,
    pub space_heat_systems: HashMap<String, Arc<Mutex<SpaceHeatSystem>>>,
    pub space_cool_systems: HashMap<String, AirConditioning>,
    pub on_site_generation: HashMap<String, PhotovoltaicSystem>,
    pub diverters: Vec<Arc<Mutex<PVDiverter>>>,
    energy_supply_conn_names_for_hot_water_source: HashMap<String, Vec<String>>,
    energy_supply_conn_names_for_heat_systems: HashMap<String, String>,
    timestep_end_calcs: Vec<Arc<Mutex<WetHeatSource>>>,
}

impl Corpus {
    pub fn from_inputs(
        input: Input,
        external_conditions: ExternalConditions,
    ) -> Result<Self, anyhow::Error> {
        let simulation_time_iterator = Arc::new(input.simulation_time.iter());

        let external_conditions = Arc::new(external_conditions);

        let diverter_types: DiverterTypes = (&input.energy_supply).into();
        let mut diverters: Vec<Arc<Mutex<PVDiverter>>> = Default::default();

        let cold_water_sources =
            cold_water_sources_from_input(input.cold_water_source, &input.simulation_time);
        let wwhrs = wwhrs_from_input(input.waste_water_heat_recovery, &cold_water_sources);

        let mut energy_supplies = energy_supplies_from_input(
            input.energy_supply,
            simulation_time_iterator.clone().as_ref(),
        );

        let controls = control_from_input(input.control, simulation_time_iterator.clone().as_ref());

        let event_schedules = event_schedules_from_input(
            input.water_heating_events,
            simulation_time_iterator.as_ref(),
        );

        let domestic_hot_water_demand = DomesticHotWaterDemand::new(
            input.shower,
            input.bath,
            input.other_water_use,
            match &input.hot_water_source.hot_water_cylinder {
                HotWaterSourceDetails::PointOfUse { .. } => None,
                _ => input.water_distribution,
            },
            &cold_water_sources,
            &wwhrs,
            &energy_supplies,
            event_schedules.clone(),
        );

        let infiltration = infiltration_from_input(input.infiltration);

        let space_heating_ductwork = ductwork_from_ventilation_input(&input.ventilation);

        let ventilation = input.ventilation.as_ref().map(|v| {
            Arc::new(Mutex::new(ventilation_from_input(
                "Ventilation system",
                v,
                &infiltration,
                simulation_time_iterator.clone().as_ref(),
                &mut energy_supplies,
            )))
        });

        let opening_area_total_from_zones = opening_area_total_from_zones(&input.zone);

        let mut heat_system_name_for_zone: HashMap<String, Option<String>> = Default::default();
        let mut cool_system_name_for_zone: HashMap<String, Option<String>> = Default::default();

        let zones: HashMap<String, Zone> = input
            .zone
            .iter()
            .map(|(i, zone)| {
                let ventilation = ventilation.as_ref().map(|ventilation| ventilation.lock());
                ((*i).clone(), {
                    let (zone_for_corpus, heat_system_name, cool_system_name) = zone_from_input(
                        zone,
                        opening_area_total_from_zones,
                        &input.window_opening_for_cooling,
                        &controls,
                        ventilation.map(|v| (*v).clone()),
                        external_conditions.clone(),
                        &infiltration,
                        simulation_time_iterator.clone().as_ref(),
                    );
                    heat_system_name_for_zone.insert((*i).clone(), heat_system_name);
                    cool_system_name_for_zone.insert((*i).clone(), cool_system_name);

                    zone_for_corpus
                })
            })
            .collect();

        if !has_unique_some_values(&heat_system_name_for_zone)
            || !has_unique_some_values(&cool_system_name_for_zone)
        {
            bail!("the heat or cool systems do not have unique names in the inputs");
        }

        let energy_supply_conn_unmet_demand_zone = set_up_energy_supply_unmet_demand_zones(
            energy_supplies.unmet_demand.clone(),
            &input.zone,
        );
        // TODO: there needs to be some equivalent here of the Python code that builds the dict __energy_supply_conn_unmet_demand_zone

        let (total_floor_area, total_volume) = zones.values().fold((0., 0.), |acc, zone| {
            (zone.area() + acc.0, zone.volume() + acc.1)
        });

        let mut internal_gains = internal_gains_from_input(input.internal_gains);

        apply_appliance_gains_from_input(
            &mut internal_gains,
            input.appliance_gains,
            &mut energy_supplies,
            total_floor_area,
            simulation_time_iterator.total_steps(),
        );

        let mut timestep_end_calcs = vec![];

        let wet_heat_sources: HashMap<String, Arc<Mutex<WetHeatSource>>> = input
            .heat_source_wet
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
                    total_volume,
                    &controls,
                    &mut energy_supplies,
                )));
                match *heat_source.lock() {
                    WetHeatSource::HeatPump(_)
                    | WetHeatSource::Boiler(_)
                    | WetHeatSource::HeatBattery(_) => {
                        timestep_end_calcs.push(heat_source.clone());
                    }
                    _ => {}
                }
                ((*name).clone(), heat_source)
            })
            .collect();

        let mut hot_water_sources: HashMap<String, HotWaterSource> = Default::default();
        let mut energy_supply_conn_names_for_hot_water_source: HashMap<String, Vec<String>> =
            Default::default();
        let (hot_water_source, hw_cylinder_conn_names) = hot_water_source_from_input(
            "hw cylinder".to_string(),
            input.hot_water_source.hot_water_cylinder,
            &cold_water_sources,
            &wet_heat_sources,
            &wwhrs,
            &controls,
            &mut energy_supplies,
            &diverter_types,
            &mut diverters,
            simulation_time_iterator.clone().as_ref(),
            external_conditions.clone(),
        );
        hot_water_sources.insert("hw cylinder".to_string(), hot_water_source);
        energy_supply_conn_names_for_hot_water_source
            .insert("hw cylinder".to_string(), hw_cylinder_conn_names);

        let mut heat_system_names_requiring_overvent: Vec<String> = Default::default();

        let (space_heat_systems, energy_supply_conn_names_for_heat_systems) = input
            .space_heat_system
            .as_ref()
            .map(|system| {
                space_heat_systems_from_input(
                    system,
                    &controls,
                    &mut energy_supplies,
                    simulation_time_iterator.as_ref(),
                    &Default::default(),
                    &mut heat_system_names_requiring_overvent,
                    heat_system_name_for_zone
                        .values()
                        .flatten()
                        .map(|s| s.as_str())
                        .collect::<Vec<_>>(),
                )
            })
            .unwrap_or_default();

        let space_cool_systems = input
            .space_cool_system
            .as_ref()
            .map(|system| {
                space_cool_systems_from_input(
                    system,
                    cool_system_name_for_zone
                        .values()
                        .flatten()
                        .map(|s| s.as_str())
                        .collect::<Vec<_>>(),
                    &controls,
                    &mut energy_supplies,
                    &simulation_time_iterator,
                )
            })
            .unwrap_or_default();

        let on_site_generation = input
            .on_site_generation
            .map(|on_site_generation| {
                on_site_generation_from_input(
                    &on_site_generation,
                    &mut energy_supplies,
                    external_conditions.clone(),
                    &simulation_time_iterator,
                )
            })
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
        for (z_name, zone) in &self.zones {
            let fabric_heat_loss = zone.total_fabric_heat_loss();
            let thermal_bridges = zone.total_thermal_bridges();
            let vent_heat_loss = zone.total_vent_heat_loss(self.external_conditions.as_ref());

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
        let internal_air_temperature = self
            .zones
            .values()
            .map(|zone| zone.temp_internal_air() * zone.volume())
            .sum::<f64>();

        internal_air_temperature / self.total_volume
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
            Some(ductwork) => {
                // MVHR duct temperatures:
                // extract_duct_temp - indoor air temperature
                // intake_duct_temp - outside air temperature

                let intake_duct_temp = self
                    .external_conditions
                    .air_temp(&simulation_time_iteration);

                let temp_diff = internal_air_temperature - intake_duct_temp;

                // Supply duct contains what the MVHR could recover
                let supply_duct_temp = intake_duct_temp + (efficiency * temp_diff);

                // Exhaust duct contans the heat that couldn't be recovered
                let exhaust_duct_temp = intake_duct_temp + ((1. - efficiency) * temp_diff);

                ductwork
                    .total_duct_heat_loss(
                        Some(internal_air_temperature),
                        Some(supply_duct_temp),
                        Some(internal_air_temperature),
                        Some(intake_duct_temp),
                        Some(exhaust_duct_temp),
                        efficiency,
                    )
                    .unwrap()
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
    ) -> SpaceHeatingCalculation {
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

        for (z_name, zone) in &self.zones {
            // Initialise to dhw internal gains split proportionally to zone floor area
            let mut gains_internal_zone_inner =
                gains_internal_dhw * zone.area() / self.total_floor_area;
            for (_, gains) in [
                (
                    "total_internal_gains",
                    &self
                        .internal_gains
                        .total_internal_gains
                        .as_ref()
                        .map(Gains::Internal),
                ),
                (
                    "metabolic_gains",
                    &self
                        .internal_gains
                        .metabolic_gains
                        .as_ref()
                        .map(Gains::Internal),
                ),
                (
                    "lighting",
                    &self.internal_gains.lighting.as_ref().map(Gains::Appliance),
                ),
                (
                    "cooking",
                    &self.internal_gains.cooking.as_ref().map(Gains::Appliance),
                ),
                (
                    "cooking1",
                    &self.internal_gains.cooking1.as_ref().map(Gains::Appliance),
                ),
                (
                    "cooking2",
                    &self.internal_gains.cooking2.as_ref().map(Gains::Appliance),
                ),
                (
                    "other",
                    &self.internal_gains.other.as_ref().map(Gains::Internal),
                ),
            ]
            .iter()
            .filter_map(|(name, option)| option.as_ref().map(|gains| (name, gains)))
            {
                gains_internal_zone_inner +=
                    gains.total_internal_gain_in_w(zone.area(), simulation_time_iteration);
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
            gains_solar_zone.insert(
                z_name,
                zone.gains_solar(self.external_conditions.as_ref(), simulation_time_iteration),
            );
        }

        // Calculate space heating and cooling demand for each zone and sum
        // Keep track of how much is from each zone, so that energy provided
        // can be split between them in same proportion later
        let (
            mut space_heat_demand_system,
            mut space_cool_demand_system,
            mut space_heat_demand_zone,
            mut space_cool_demand_zone,
            mut h_ve_cool_extra_zone,
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
                (space_heat_running_time_cumulative, throughput_factor) =
                    heat_system.lock().running_time_throughput_factor(
                        space_heat_demand_system[heat_system_name],
                        space_heat_running_time_cumulative,
                        simulation_time_iteration,
                    );
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
            for (z_name, zone) in &self.zones {
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
                h_ve_cool_extra_zone,
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
                heat_system.lock().demand_energy(
                    space_heat_demand_system[heat_system_name.as_str()],
                    simulation_time_iteration,
                ),
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
        let mut heat_balance_map: HashMap<&str, Option<()>> = Default::default(); // using unit type here as placeholder
        for (z_name, zone) in &self.zones {
            // Look up names of relevant heating and cooling systems for this zone
            let h_name = self.heat_system_name_for_zone[z_name.as_str()].as_ref();
            let c_name = self.cool_system_name_for_zone[z_name.as_str()].as_ref();

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
                    h_ve_cool_extra_zone.get(z_name.as_str()).copied(),
                    Some(throughput_factor),
                    simulation_time_iteration,
                    self.external_conditions.as_ref(),
                ),
            );

            if let Some(h_name) = h_name {
                *space_heat_demand_system.get_mut(h_name.as_str()).unwrap() = 0.0;
                *space_heat_provided.get_mut(h_name.as_str()).unwrap() = 0.0;
            }
            if let Some(c_name) = c_name {
                *space_cool_demand_system.get_mut(c_name.as_str()).unwrap() = 0.0;
                *space_cool_provided.get_mut(c_name.as_str()).unwrap() = 0.0;
            }

            internal_air_temp.insert(z_name.as_str(), zone.temp_internal_air());
            operative_temp.insert(z_name.as_str(), zone.temp_operative());
        }

        (
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
        )
    }

    pub fn run(&mut self) -> RunResults {
        let simulation_time = self.simulation_time.as_ref().to_owned();
        let vec_capacity = || Vec::with_capacity(simulation_time.total_steps());

        let mut timestep_array = vec_capacity();
        let mut gains_internal_dict: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut gains_solar_dict: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut operative_temp_dict: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut internal_air_temp_dict: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut space_heat_demand_dict: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut space_cool_demand_dict: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut space_heat_demand_system_dict: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut space_cool_demand_system_dict: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut space_heat_provided_dict: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut space_cool_provided_dict: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut zone_list: Vec<KeyString> = Default::default();
        let mut hot_water_demand_dict: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut hot_water_energy_demand_dict: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut hot_water_energy_demand_dict_incl_pipework: HashMap<KeyString, Vec<f64>> =
            Default::default();
        let mut hot_water_energy_output_dict: HashMap<&str, Vec<f64>> = Default::default();
        let mut hot_water_duration_dict: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut hot_water_no_events_dict: HashMap<KeyString, Vec<usize>> = Default::default();
        let mut hot_water_pipework_dict: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut ductwork_gains_dict: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut heat_balance_all_dict: HashMap<
            KeyString,
            HashMap<KeyString, HashMap<KeyString, f64>>,
        > = HashMap::from([
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

        for h_name in self.heat_system_name_for_zone.values().flatten() {
            let h_name = h_name.as_str().try_into().unwrap();
            space_heat_demand_system_dict.insert(h_name, vec_capacity());
            space_heat_provided_dict.insert(h_name, vec_capacity());
        }

        for c_name in self.cool_system_name_for_zone.values().flatten() {
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

        let progress_bar = ProgressBar::new(simulation_time.total_steps() as u64);

        for t_it in simulation_time {
            timestep_array.push(t_it.time);
            let (hw_demand_vol, hw_vol_at_tapping_points, hw_duration, no_events, hw_energy_demand) =
                self.domestic_hot_water_demand.hot_water_demand(t_it.index);

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

            let hw_energy_output = self
                .hot_water_sources
                .get_mut("hw cylinder")
                .unwrap()
                .demand_hot_water(hw_demand_vol, t_it);

            let (pw_losses_internal, pw_losses_external, gains_internal_dhw_use) = self
                .pipework_losses_and_internal_gains_from_hw(
                    t_it.timestep,
                    hw_vol_at_tapping_points,
                    hw_duration,
                    no_events,
                    t_it,
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
            ) = self.calc_space_heating(t_it.timestep, gains_internal_dhw, t_it);

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
                diverter.lock().timestep_end();
            }

            progress_bar.inc(1);
        }

        progress_bar.finish();

        // Return results from all energy supplies
        let mut results_totals: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut results_end_user: HashMap<KeyString, IndexMap<String, Vec<f64>>> =
            Default::default();
        let mut energy_import: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut energy_export: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut energy_generated_consumed: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut energy_to_storage: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut energy_from_storage: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut energy_diverted: HashMap<KeyString, Vec<f64>> = Default::default();
        let mut betafactor: HashMap<KeyString, Vec<f64>> = Default::default();
        for (name, supply) in self.energy_supplies.supplies_by_name() {
            let supply = supply.lock();
            let name: KeyString = name.try_into().unwrap();
            results_totals.insert(name, supply.results_total().to_owned());
            results_end_user.insert(name, supply.results_by_end_user());
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

        let hot_water_energy_out: HashMap<KeyString, Vec<f64>> = HashMap::from([(
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
            &space_cool_provided_dict,
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

        let zone_dict = HashMap::from([
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
        let hc_system_dict = HashMap::from([
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
        let hot_water_dict = HashMap::from([
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
        energy_provided: &HashMap<KeyString, Vec<f64>>,
        results_end_user: &HashMap<KeyString, IndexMap<String, Vec<f64>>>,
        energy_supply_conn_name_for_space_hc_system: HashMap<String, Vec<String>>,
    ) -> HashMap<KeyString, NumberOrDivisionByZero> {
        let mut hc_output_overall: HashMap<KeyString, f64> = Default::default();
        let mut hc_input_overall: HashMap<KeyString, f64> = Default::default();
        let mut cop_dict: HashMap<KeyString, NumberOrDivisionByZero> = Default::default();
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
        throughput_factor: Option<f64>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> (
        HashMap<String, f64>,
        HashMap<String, f64>,
        HashMap<String, f64>,
        HashMap<String, f64>,
        HashMap<String, f64>,
    ) {
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
        for (z_name, zone) in &self.zones {
            // Look up names of relevant heating and cooling systems for this zone
            let h_name = self
                .heat_system_name_for_zone
                .get(z_name)
                .and_then(|h| h.as_ref());
            let c_name = self
                .cool_system_name_for_zone
                .get(z_name)
                .and_then(|c| c.as_ref());

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
            ) = zone.space_heat_cool_demand(
                delta_t_h,
                temp_ext_air,
                gains_internal_zone[z_name.as_str()],
                gains_solar_zone[z_name.as_str()],
                frac_convective_heat,
                frac_convective_cool,
                temp_setpnt_heat,
                temp_setpnt_cool,
                throughput_factor,
                simulation_time_iteration,
                self.external_conditions.as_ref(),
            );

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
    Float(HashMap<KeyString, Vec<f64>>),
    Int(HashMap<KeyString, Vec<usize>>),
}

#[derive(Clone, Debug)]
pub enum NumberOrDivisionByZero {
    Number(f64),
    DivisionByZero,
}

fn has_unique_some_values<K, V: Eq + Hash>(map: &HashMap<K, Option<V>>) -> bool {
    let some_values: Vec<&V> = map.values().flat_map(|v| v.iter()).collect();
    let value_set: HashSet<&&V> = some_values.iter().collect();
    some_values.len() == value_set.len()
}

fn external_conditions_from_input(
    input: Arc<ExternalConditionsInput>,
    simulation_time: &SimulationTimeIterator,
) -> ExternalConditions {
    ExternalConditions::new(
        simulation_time,
        input.air_temperatures.clone().unwrap_or_default(),
        input.wind_speeds.clone().unwrap_or_default(),
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

fn infiltration_from_input(input: Infiltration) -> VentilationElementInfiltration {
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
    input: ColdWaterSourceInput,
    simulation_time: &SimulationTime,
) -> ColdWaterSources {
    ColdWaterSources {
        mains_water: input
            .mains_water
            .map(|details| cold_water_source_from_input_details(details, simulation_time)),
        header_tank: input
            .header_tank
            .map(|details| cold_water_source_from_input_details(details, simulation_time)),
    }
}

fn cold_water_source_from_input_details(
    details: ColdWaterSourceDetails,
    simulation_time: &SimulationTime,
) -> ColdWaterSource {
    ColdWaterSource::new(
        details.temperatures,
        simulation_time,
        details.time_series_step,
    )
}

fn energy_supplies_from_input(
    input: EnergySupplyInput,
    simulation_time_iterator: &SimulationTimeIterator,
) -> EnergySupplies {
    EnergySupplies {
        mains_electricity: energy_supply_from_input(
            input.mains_electricity,
            simulation_time_iterator,
        ),
        mains_gas: energy_supply_from_input(input.mains_gas, simulation_time_iterator),
        bulk_lpg: energy_supply_from_input(input.bulk_lpg, simulation_time_iterator),
        heat_network: energy_supply_from_heat_network_input(
            input.heat_network,
            simulation_time_iterator,
        ),
        unmet_demand: Arc::new(Mutex::new(EnergySupply::new(
            EnergySupplyType::UnmetDemand,
            simulation_time_iterator.total_steps(),
            Default::default(),
        ))),
        bottled_lpg: None,
        condition_11f_lpg: None,
        custom: None,
    }
}

fn energy_supply_from_input(
    input: Option<EnergySupplyDetails>,
    simulation_time_iterator: &SimulationTimeIterator,
) -> Option<Arc<Mutex<EnergySupply>>> {
    input.map(|details| {
        Arc::new(Mutex::new(EnergySupply::new(
            details.fuel,
            simulation_time_iterator.total_steps(),
            details.electric_battery.map(ElectricBattery::from_input),
        )))
    })
}

fn energy_supply_from_heat_network_input(
    input: Option<HeatNetworkInput>,
    simulation_time_iterator: &SimulationTimeIterator,
) -> Option<Arc<Mutex<EnergySupply>>> {
    input.map(|heat_network| {
        Arc::new(Mutex::new(EnergySupply::new(
            heat_network.fuel,
            simulation_time_iterator.total_steps(),
            None,
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
            mains_electricity: diverter_from_energy_supply(&input.mains_electricity),
            mains_gas: diverter_from_energy_supply(&input.mains_gas),
            bulk_lpg: diverter_from_energy_supply(&input.bulk_lpg),
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

fn diverter_from_energy_supply(supply: &Option<EnergySupplyDetails>) -> Option<EnergyDiverter> {
    supply
        .as_ref()
        .map(|supply| supply.diverter.clone().unwrap_or_default())
}

pub struct InternalGainsCollection {
    total_internal_gains: Option<InternalGains>,
    metabolic_gains: Option<InternalGains>,
    lighting: Option<ApplianceGains>,
    cooking: Option<ApplianceGains>,
    cooking1: Option<ApplianceGains>,
    cooking2: Option<ApplianceGains>,
    other: Option<InternalGains>,
}

fn internal_gains_from_input(input: InternalGainsInput) -> InternalGainsCollection {
    InternalGainsCollection {
        total_internal_gains: input.total_internal_gains.map(internal_gains_from_details),
        metabolic_gains: input.metabolic_gains.map(internal_gains_from_details),
        lighting: None,
        cooking: None,
        cooking1: None,
        cooking2: None,
        other: input.other.map(internal_gains_from_details),
    }
}

fn internal_gains_from_details(details: InternalGainsDetails) -> InternalGains {
    let mut schedule = HashMap::from([("main".to_string(), details.schedule.main)]);
    if let Some(day) = details.schedule.day {
        schedule.insert("day".to_string(), day);
    }
    InternalGains::new(
        expand_numeric_schedule(schedule, false),
        details.start_day,
        details.time_series_step,
    )
}

pub struct Controls {
    core: Vec<HeatSourceControl<Option<Arc<Control>>>>,
    extra: HashMap<String, Arc<Control>>,
}

impl Controls {
    pub fn new(
        core: Vec<HeatSourceControl<Option<Arc<Control>>>>,
        extra: HashMap<String, Arc<Control>>,
    ) -> Self {
        Self { core, extra }
    }

    pub fn get(&self, control_type: &HeatSourceControlType) -> Option<&Arc<Control>> {
        self.core
            .iter()
            .find(|heat_source_control| heat_source_control.get(control_type).is_some())
            .and_then(|heat_source_control| heat_source_control.get(control_type).unwrap().as_ref())
    }

    // access a control using a string, possibly because it is one of the "extra" controls
    pub fn get_with_string(&self, control_name: &str) -> Option<&Arc<Control>> {
        match control_name {
            // hard-code ways of resolving to core control types (for now)
            "hw timer" => self.get(&HeatSourceControlType::HotWaterTimer),
            "window opening" => self.get(&HeatSourceControlType::WindowOpening),
            other => self.extra.get(other),
        }
    }
}

fn control_from_input(
    control_input: ControlInput,
    simulation_time_iterator: &SimulationTimeIterator,
) -> Controls {
    let mut core: Vec<HeatSourceControl<Option<Arc<Control>>>> = Default::default();
    let mut extra: HashMap<String, Arc<Control>> = Default::default();

    // this is very ugly(!) but is just a reflection of the lack of clarity in the schema
    // and the way the variants-struct crate works;
    // we should be able to improve it in time
    for control in control_input.core {
        match control {
            HeatSourceControl {
                hot_water_timer: Some(control),
                ..
            } => {
                core.push(HeatSourceControl::new(
                    Some(Arc::new(single_control_from_details(
                        control,
                        simulation_time_iterator,
                    ))),
                    None,
                ));
            }
            HeatSourceControl {
                window_opening: Some(control),
                ..
            } => {
                core.push(HeatSourceControl::new(
                    None,
                    Some(Arc::new(single_control_from_details(
                        control,
                        simulation_time_iterator,
                    ))),
                ));
            }
            unknown => panic!(
                "incorrectly formed HeatSourceControl struct encountered: {:?}",
                unknown
            ),
        }
    }
    for (name, control) in control_input.extra {
        extra.insert(
            name,
            Arc::new(single_control_from_details(
                control,
                simulation_time_iterator,
            )),
        );
    }

    Controls { core, extra }
}

fn single_control_from_details(
    details: ControlDetails,
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
            start_day,
            time_series_step,
        )),
        ControlDetails::OnOffCostMinimisingTime {
            start_day,
            time_series_step,
            time_on_daily,
            schedule,
            ..
        } => Control::OnOffMinimisingTimeControl(OnOffMinimisingTimeControl::new(
            expand_numeric_schedule(schedule, false),
            start_day,
            time_series_step,
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
                expand_numeric_schedule(schedule, true)
                    .iter()
                    .map(|s| Some(*s))
                    .collect(),
                start_day,
                time_series_step,
                setpoint_min,
                setpoint_max,
                default_to_max,
                advanced_start,
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
                start_day,
                time_series_step,
                charge_level: charge_level_vec,
            })
        }
    }
}

fn wwhrs_from_input(
    wwhrs: Option<WasteWaterHeatRecovery>,
    cold_water_sources: &ColdWaterSources,
) -> HashMap<String, Wwhrs> {
    let mut wwhr_systems: HashMap<String, Wwhrs> = HashMap::from([]);
    if let Some(systems) = wwhrs {
        for (name, system) in systems {
            wwhr_systems
                .entry(name)
                .or_insert(wwhr_system_from_details(system, cold_water_sources));
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

pub type KeyString = ArrayString<48>;

pub type RunResults = (
    Vec<f64>,
    HashMap<KeyString, Vec<f64>>,
    HashMap<KeyString, IndexMap<String, Vec<f64>>>,
    HashMap<KeyString, Vec<f64>>,
    HashMap<KeyString, Vec<f64>>,
    HashMap<KeyString, Vec<f64>>,
    HashMap<KeyString, Vec<f64>>,
    HashMap<KeyString, Vec<f64>>,
    HashMap<KeyString, Vec<f64>>,
    HashMap<KeyString, Vec<f64>>,
    HashMap<KeyString, HashMap<KeyString, Vec<f64>>>,
    Vec<KeyString>,
    HashMap<KeyString, HashMap<KeyString, Vec<f64>>>,
    HashMap<KeyString, HotWaterResultMap>,
    HashMap<KeyString, NumberOrDivisionByZero>,
    HashMap<KeyString, NumberOrDivisionByZero>,
    HashMap<KeyString, NumberOrDivisionByZero>,
    HashMap<KeyString, Vec<f64>>,
    HashMap<KeyString, HashMap<KeyString, HashMap<KeyString, f64>>>,
    HashMap<KeyString, f64>,
    HashMap<KeyString, f64>,
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
    HashMap<&'a str, Option<()>>,
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

pub type EventSchedule = Vec<Option<Vec<ScheduleEvent>>>;

#[derive(Clone)]
pub struct HotWaterEventSchedules {
    pub shower: HashMap<String, EventSchedule>,
    pub bath: HashMap<String, EventSchedule>,
    pub other: HashMap<String, EventSchedule>,
}

fn event_schedules_from_input(
    events: WaterHeatingEvents,
    simulation_time_iterator: &SimulationTimeIterator,
) -> HotWaterEventSchedules {
    let mut shower_schedules: HashMap<String, EventSchedule> = Default::default();
    if let Some(shower_events) = events.shower {
        shower_schedules.insert(
            "ies".to_string(),
            schedule_event_from_input(shower_events.ies.iter().collect(), simulation_time_iterator),
        );
        shower_schedules.insert(
            "mixer".to_string(),
            schedule_event_from_input(
                shower_events.mixer.iter().collect(),
                simulation_time_iterator,
            ),
        );
    }

    let mut bath_schedules: HashMap<String, EventSchedule> = Default::default();
    if let Some(bath_events) = events.bath {
        bath_schedules.insert(
            "medium".to_string(),
            schedule_event_from_input(
                bath_events.medium.iter().collect(),
                simulation_time_iterator,
            ),
        );
    }

    let mut other_schedules: HashMap<String, EventSchedule> = Default::default();
    if let Some(other_events) = events.other {
        other_schedules.insert(
            "other".to_string(),
            schedule_event_from_input(
                other_events.other.iter().collect(),
                simulation_time_iterator,
            ),
        );
    }

    HotWaterEventSchedules {
        shower: shower_schedules,
        bath: bath_schedules,
        other: other_schedules,
    }
}

fn schedule_event_from_input(
    events_input: Vec<&WaterHeatingEvent>,
    simulation_time_iterator: &SimulationTimeIterator,
) -> EventSchedule {
    let sim_timestep = simulation_time_iterator.step_in_hours();
    let total_timesteps = simulation_time_iterator.total_steps();
    expand_water_heating_events(events_input, sim_timestep, total_timesteps)
}

fn ductwork_from_ventilation_input(ventilation: &Option<Ventilation>) -> Option<Ductwork> {
    ventilation.as_ref().and_then(|v| match v {
        Ventilation::Mvhr { ductwork, .. } => Some(Ductwork::new(
            ductwork.internal_diameter_mm / MILLIMETRES_IN_METRE as f64,
            ductwork.external_diameter_mm / MILLIMETRES_IN_METRE as f64,
            ductwork.length_in,
            ductwork.length_out,
            ductwork.insulation_thermal_conductivity,
            ductwork.insulation_thickness_mm / MILLIMETRES_IN_METRE as f64,
            ductwork.reflective,
            ductwork.mvhr_location,
        )),
        _ => None,
    })
}

fn ventilation_from_input<'a>(
    energy_supply_name: &str,
    ventilation: &'a Ventilation,
    infiltration: &VentilationElementInfiltration,
    simulation_time: &'a SimulationTimeIterator,
    energy_supplies: &mut EnergySupplies,
) -> VentilationElement {
    match ventilation {
        Ventilation::Whev {
            req_ach,
            sfp,
            energy_supply,
        } => {
            let energy_supply =
                energy_supplies.ensured_get_for_type(*energy_supply, simulation_time.total_steps());
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
        Ventilation::Mvhr {
            req_ach,
            sfp,
            efficiency,
            energy_supply,
            ..
        } => {
            let energy_supply =
                energy_supplies.ensured_get_for_type(*energy_supply, simulation_time.total_steps());
            let energy_supply_conn =
                EnergySupply::connection(energy_supply, energy_supply_name).unwrap();

            VentilationElement::Mvhr(MechanicalVentilationHeatRecovery::new(
                *req_ach,
                *sfp,
                *efficiency,
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
    }
}

fn energy_supply_from_type_for_ventilation(
    energy_supply_type: &EnergySupplyType,
    energy_supplies: EnergySupplies,
) -> Arc<Mutex<EnergySupply>> {
    match energy_supply_type {
        EnergySupplyType::Electricity => energy_supplies.mains_electricity.unwrap().clone(),
        EnergySupplyType::MainsGas => energy_supplies.mains_gas.unwrap().clone(),
        EnergySupplyType::UnmetDemand => energy_supplies.unmet_demand.clone(),
        EnergySupplyType::LpgBulk => energy_supplies.bulk_lpg.unwrap().clone(),
        // commenting out for now as this is a different type and might not be used in this context
        // EnergySupplyType::HeatNetwork => energy_supplies.heat_network.unwrap(),
        _ => panic!("Unexpected energy supply type listed for ventilation."),
    }
}

fn opening_area_total_from_zones(zones: &ZoneDictionary) -> f64 {
    zones
        .iter()
        .flat_map(|(_, zone)| {
            zone.building_elements
                .iter()
                .map(|(_, building_element)| match building_element {
                    BuildingElement::Transparent { height, width, .. } => height * width,
                    _ => 0.,
                })
        })
        .sum()
}

fn check_space_heat_systems_unique_to_zones(zones: &ZoneDictionary) -> Result<(), &'static str> {
    let res = {
        let mut name_set = HashSet::new();
        zones.iter().all(|(_, zone)| match &zone.space_heat_system {
            Some(system_name) => name_set.insert(system_name),
            None => true,
        })
    };
    if res {
        Ok(())
    } else {
        Err("A space heat system was declared as used in more than one zone.")
    }
}

fn check_space_cool_systems_unique_to_zones(zones: &ZoneDictionary) -> Result<(), &'static str> {
    let res = {
        let mut name_set = HashSet::new();
        zones.iter().all(|(_, zone)| match &zone.space_cool_system {
            Some(system_name) => name_set.insert(system_name),
            None => true,
        })
    };
    if res {
        Ok(())
    } else {
        Err("A space cool system was declared as used in more than one zone.")
    }
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
    opening_area_total: f64,
    window_opening_for_cooling: &Option<WindowOpeningForCoolingInput>,
    controls: &'a Controls,
    ventilation: Option<VentilationElement>,
    external_conditions: Arc<ExternalConditions>,
    infiltration: &'a VentilationElementInfiltration,
    simulation_time_iterator: &'a SimulationTimeIterator,
) -> (Zone, Option<String>, Option<String>) {
    let ventilation = ventilation.as_ref();
    let heat_system_name = input.space_heat_system.clone();
    let cool_system_name = input.space_cool_system.clone();

    let vent_cool_extra = window_opening_for_cooling.as_ref().map(|opening| {
        let openings: HashMap<&String, &BuildingElement> = input
            .building_elements
            .iter()
            .filter(|(_, el)| matches!(el, BuildingElement::Transparent { .. }))
            .collect();
        let opening_area_zone: f64 = openings
            .values()
            .map(|op| area_for_building_element_input(op))
            .sum();
        let opening_area_equivalent =
            opening.equivalent_area * opening_area_zone / opening_area_total;
        let control = input
            .control_window_opening
            .as_ref()
            .and_then(|opening_control| {
                controls
                    .get(opening_control)
                    .and_then(|c| match c.as_ref() {
                        Control::SetpointTimeControl(ctrl) => Some((*ctrl).clone()),
                        _ => None,
                    })
            });
        let natvent = ventilation.and_then(|v| match v {
            VentilationElement::Natural(natural_ventilation) => {
                Some((*natural_ventilation).clone())
            }
            _ => None,
        });
        let named_openings: Vec<NamedBuildingElement> = openings
            .iter()
            .map(|(name, element)| NamedBuildingElement {
                name: name.to_string(),
                element: (*element).clone(),
            })
            .collect::<Vec<_>>();
        WindowOpeningForCooling::new(
            opening_area_equivalent,
            external_conditions.clone(),
            named_openings,
            control,
            natvent,
        )
    });

    let infiltration_ventilation = VentilationElement::Infiltration((*infiltration).clone());
    let mut vent_elements: Vec<VentilationElement> = vec![infiltration_ventilation];

    if let Some(v) = ventilation {
        vent_elements.push((*v).clone());
    }

    (
        Zone::new(
            input.area,
            input.volume,
            input.building_elements.clone(),
            thermal_bridging_from_input(&input.thermal_bridging),
            vent_elements,
            vent_cool_extra,
            external_conditions.air_temp(&simulation_time_iterator.current_iteration()),
            input.temp_setpnt_init.unwrap(),
            external_conditions.clone(),
            simulation_time_iterator,
        ),
        heat_system_name,
        cool_system_name,
    )
}

fn set_up_energy_supply_unmet_demand_zones(
    unmet_demand_supply: Arc<Mutex<EnergySupply>>,
    zones: &ZoneDictionary,
) -> HashMap<String, Arc<EnergySupplyConnection>> {
    let mut energy_supplies: HashMap<String, Arc<EnergySupplyConnection>> = Default::default();

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
    input: ApplianceGainsInput,
    energy_supplies: &mut EnergySupplies,
    total_floor_area: f64,
    simulation_timesteps: usize,
) {
    if let Some(details) = input.lighting {
        let energy_supply_conn = EnergySupply::connection(
            energy_supplies.ensured_get_for_type(details.energy_supply, simulation_timesteps),
            "lighting",
        )
        .unwrap();

        internal_gains_collection.lighting = Some(appliance_gains_from_single_input(
            details,
            energy_supply_conn,
            total_floor_area,
        ));
    }
    if let Some(details) = input.cooking {
        let energy_supply_conn = EnergySupply::connection(
            energy_supplies.ensured_get_for_type(details.energy_supply, simulation_timesteps),
            "cooking",
        )
        .unwrap();

        internal_gains_collection.cooking = Some(appliance_gains_from_single_input(
            details,
            energy_supply_conn,
            total_floor_area,
        ));
    }
    if let Some(details) = input.cooking1 {
        let energy_supply_conn = EnergySupply::connection(
            energy_supplies.ensured_get_for_type(details.energy_supply, simulation_timesteps),
            "cooking1",
        )
        .unwrap();

        internal_gains_collection.cooking1 = Some(appliance_gains_from_single_input(
            details,
            energy_supply_conn,
            total_floor_area,
        ));
    }
    if let Some(details) = input.cooking2 {
        let energy_supply_conn = EnergySupply::connection(
            energy_supplies.ensured_get_for_type(details.energy_supply, simulation_timesteps),
            "cooking2",
        )
        .unwrap();

        internal_gains_collection.cooking2 = Some(appliance_gains_from_single_input(
            details,
            energy_supply_conn,
            total_floor_area,
        ));
    }
}

fn appliance_gains_from_single_input(
    input: ApplianceGainsDetails,
    energy_supply_connection: EnergySupplyConnection,
    total_floor_area: f64,
) -> ApplianceGains {
    let total_energy_supply = expand_numeric_schedule(input.schedule, false)
        .iter()
        .map(|energy_data| energy_data / total_floor_area)
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
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        match self {
            HeatSource::Storage(ref mut storage) => match storage {
                HeatSourceWithStorageTank::Immersion(imm) => imm
                    .lock()
                    .demand_energy(energy_demand, simulation_time_iteration),
                HeatSourceWithStorageTank::Solar(ref mut solar) => {
                    solar.demand_energy(energy_demand, simulation_time_iteration.index)
                }
            },
            HeatSource::Wet(ref mut wet) => match wet.as_mut() {
                HeatSourceWet::WaterCombi(_) => {
                    unimplemented!("not expected? this value does not have a demand_energy method")
                    // the Python uses duck-typing here but there is no method for this type
                }
                HeatSourceWet::WaterRegular(ref mut r) => {
                    r.demand_energy(energy_demand, simulation_time_iteration)
                }
                HeatSourceWet::Space(_) => {
                    unimplemented!("not expected? this value does not have a demand_energy method")
                    // the Python uses duck-typing here but there is no method for this type
                }
                HeatSourceWet::HeatNetworkWaterStorage(ref mut h) => {
                    h.demand_energy(energy_demand, &simulation_time_iteration)
                }
                HeatSourceWet::HeatBatteryHotWater(ref mut h) => {
                    h.demand_energy(energy_demand, simulation_time_iteration)
                }
                HeatSourceWet::HeatPumpWater(ref mut h) => {
                    h.demand_energy(energy_demand, simulation_time_iteration)
                }
                HeatSourceWet::HeatPumpWaterOnly(h) => {
                    h.demand_energy(energy_demand, simulation_time_iteration)
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
    total_volume: f64,
    controls: &Controls,
    energy_supplies: &mut EnergySupplies,
) -> WetHeatSource {
    match &input {
        HeatSourceWetDetails::HeatPump {
            source_type,
            energy_supply_heat_network,
            energy_supply,
            ..
        } => {
            let throughput_exhaust_air = if source_type.is_exhaust_air() {
                // Check that ventilation system is compatible with exhaust air HP
                if ventilation.is_none()
                    || !matches!(
                        ventilation.unwrap(),
                        VentilationElement::Mvhr(_) | VentilationElement::Whev(_)
                    )
                {
                    panic!("Exhaust air heat pump requires ventilation to be MVHR or WHEV.")
                }
                Some(
                    air_change_rate_to_flow_rate(ventilation_req_ach.unwrap(), total_volume)
                        * LITRES_PER_CUBIC_METRE as f64,
                )
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

            let energy_supply =
                energy_supplies.ensured_get_for_type(*energy_supply, simulation_time.total_steps());
            let energy_supply_conn_name_auxiliary =
                format!("HeatPump_auxiliary: {energy_supply_name}");

            WetHeatSource::HeatPump(Arc::new(Mutex::new(
                HeatPump::new(
                    &input,
                    energy_supply,
                    &energy_supply_conn_name_auxiliary,
                    simulation_time.step_in_hours(),
                    external_conditions.clone(),
                    throughput_exhaust_air,
                    energy_supply_hn,
                    DETAILED_OUTPUT_HEATING_COOLING,
                )
                .unwrap(),
            )))
        }
        HeatSourceWetDetails::Boiler {
            energy_supply,
            energy_supply_auxiliary,
            ..
        } => {
            let energy_supply =
                energy_supplies.ensured_get_for_type(*energy_supply, simulation_time.total_steps());
            let energy_supply_aux = energy_supplies
                .ensured_get_for_type(*energy_supply_auxiliary, simulation_time.total_steps());
            let aux_supply_name = format!("Boiler_auxiliary: {energy_supply_name}");
            let energy_supply_conn_aux =
                EnergySupply::connection(energy_supply_aux.clone(), aux_supply_name.as_str())
                    .unwrap();

            WetHeatSource::Boiler(
                Boiler::new(
                    input,
                    energy_supply,
                    energy_supply_conn_aux,
                    external_conditions.clone(),
                    simulation_time.step_in_hours(),
                )
                .expect("could not construct boiler value from provided data"),
            )
        }
        HeatSourceWetDetails::Hiu {
            power_max,
            hiu_daily_loss,
            building_level_distribution_losses,
            energy_supply,
            ..
        } => {
            let energy_supply =
                energy_supplies.ensured_get_for_type(*energy_supply, simulation_time.total_steps());
            let energy_supply_conn_name_auxiliary =
                format!("HeatNetwork_auxiliary: {energy_supply_name}");
            let energy_supply_conn_name_building_level_distribution_losses =
                format!("HeatNetwork_building_level_distribution_losses: {energy_supply_name}");

            WetHeatSource::Hiu(HeatNetwork::new(
                *power_max,
                *hiu_daily_loss,
                *building_level_distribution_losses,
                energy_supply,
                energy_supply_conn_name_auxiliary,
                energy_supply_conn_name_building_level_distribution_losses,
                simulation_time.step_in_hours(),
            ))
            // TODO add heat network to timestep_end_calcs
        }
        HeatSourceWetDetails::HeatBattery {
            control_charge,
            energy_supply,
            ..
        } => {
            let energy_supply =
                energy_supplies.ensured_get_for_type(*energy_supply, simulation_time.total_steps());
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
                external_conditions.clone(),
            ));
            heat_source
        }
    }
}

fn heat_source_from_input(
    name: &str,
    input: &HeatSourceInput,
    temp_setpoint: f64,
    wet_heat_sources: &HashMap<String, Arc<Mutex<WetHeatSource>>>,
    simulation_time: &SimulationTimeIterator,
    controls: &Controls,
    energy_supplies: &mut EnergySupplies,
    cold_water_sources: &ColdWaterSources,
    external_conditions: Arc<ExternalConditions>,
) -> (HeatSource, String) {
    // TODO add in all the stuff to do with energy supply

    match input {
        HeatSourceInput::ImmersionHeater {
            power,
            control,
            energy_supply,
            ..
        } => {
            let energy_supply =
                energy_supplies.ensured_get_for_type(*energy_supply, simulation_time.total_steps());
            let energy_supply_conn = EnergySupply::connection(energy_supply.clone(), name).unwrap();

            (
                HeatSource::Storage(HeatSourceWithStorageTank::Immersion(Arc::new(Mutex::new(
                    ImmersionHeater::new(
                        *power,
                        energy_supply_conn,
                        simulation_time.step_in_hours(),
                        control
                            .clone()
                            .and_then(|ctrl| controls.get(&ctrl).map(|c| (*c).clone())),
                    ),
                )))),
                name.into(),
            )
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
            let energy_supply =
                energy_supplies.ensured_get_for_type(*energy_supply, simulation_time.total_steps());
            let energy_supply_conn = EnergySupply::connection(energy_supply.clone(), name).unwrap();

            (
                HeatSource::Storage(HeatSourceWithStorageTank::Solar(SolarThermalSystem::new(
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
                    simulation_time.step_in_hours(),
                    WATER.clone(),
                ))),
                name.into(),
            )
        }
        HeatSourceInput::Wet {
            name,
            cold_water_source: cold_water_source_type,
            control,
            temp_flow_limit_upper,
            ..
        } => {
            let cold_water_source = cold_water_sources
                .ref_for_type(*cold_water_source_type)
                .expect("Expected a cold water source to be available to a boiler heat source.");
            let energy_supply_conn_name = format!("{name}_water_heating");
            let heat_source_wet = wet_heat_sources
                .get(name)
                .unwrap_or_else(|| {
                    panic!("Expected a wet heat source registered with the name '{name}'.")
                })
                .clone();
            let source_control = control
                .clone()
                .and_then(|ctrl| controls.get(&ctrl).map(|c| (*c).clone()));

            let lock = heat_source_wet.lock();
            let mut heat_source_wet_clone = (*lock).clone();

            (
                match heat_source_wet_clone {
                    WetHeatSource::HeatPump(heat_pump) => HeatSource::Wet(Box::new(
                        HeatSourceWet::HeatPumpWater(HeatPump::create_service_hot_water(
                            heat_pump.clone(),
                            energy_supply_conn_name.clone(),
                            temp_setpoint,
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
                            cold_water_source,
                            55.,
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
                                Arc::new(cold_water_source),
                                55.,
                                source_control,
                            ),
                        )))
                    }
                },
                energy_supply_conn_name,
            )
        }
        HeatSourceInput::HeatPumpHotWaterOnly {
            power_max,
            vol_hw_daily_average,
            ref test_data,
            energy_supply,
            control,
            ..
        } => {
            let energy_supply =
                energy_supplies.ensured_get_for_type(*energy_supply, simulation_time.total_steps());
            let energy_supply_conn_name = name;
            let energy_supply_connection =
                EnergySupply::connection(energy_supply.clone(), energy_supply_conn_name).unwrap();

            (
                HeatSource::Wet(Box::new(HeatSourceWet::HeatPumpWaterOnly(
                    HeatPumpHotWaterOnly::new(
                        *power_max,
                        energy_supply_connection,
                        test_data,
                        *vol_hw_daily_average,
                        simulation_time.step_in_hours(),
                        controls.get(control).map(|c| (*c).clone()),
                    ),
                ))),
                energy_supply_conn_name.into(),
            )
        }
    }
}

enum HotWaterSource {
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
        vol_demanded: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> f64 {
        match self {
            HotWaterSource::StorageTank(ref mut source) => source
                .lock()
                .demand_hot_water(vol_demanded, simulation_time_iteration),
            HotWaterSource::CombiBoiler(ref mut source) => {
                source.demand_hot_water(vol_demanded, simulation_time_iteration)
            }
            HotWaterSource::PointOfUse(ref mut source) => {
                source.demand_hot_water(vol_demanded, &simulation_time_iteration)
            }
            HotWaterSource::HeatNetwork(ref mut source) => {
                source.demand_hot_water(vol_demanded, simulation_time_iteration.index)
            }
            HotWaterSource::HeatBattery(_) => Default::default(),
        }
    }
}

fn hot_water_source_from_input(
    source_name: String,
    input: HotWaterSourceDetails,
    cold_water_sources: &ColdWaterSources,
    wet_heat_sources: &HashMap<String, Arc<Mutex<WetHeatSource>>>,
    wwhrs: &HashMap<String, Wwhrs>,
    controls: &Controls,
    energy_supplies: &mut EnergySupplies,
    diverter_types: &DiverterTypes,
    diverters: &mut Vec<Arc<Mutex<PVDiverter>>>,
    simulation_time: &SimulationTimeIterator,
    external_conditions: Arc<ExternalConditions>,
) -> (HotWaterSource, Vec<String>) {
    let mut energy_supply_conn_names = vec![];
    let cloned_input = input.clone();
    let hot_water_source = match input {
        HotWaterSourceDetails::StorageTank {
            volume,
            daily_losses,
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
                    match heat_recovery_system {
                        Wwhrs::WWHRSInstantaneousSystemC(c) => {
                            cold_water_source =
                                WaterSourceWithTemperature::WwhrsC(Arc::new((*c).clone()));
                        }
                        Wwhrs::WWHRSInstantaneousSystemA(a) => {
                            cold_water_source =
                                WaterSourceWithTemperature::WwhrsA(Arc::new((*a).clone()));
                        }
                        _ => {}
                    }
                }
            }
            let pipework = primary_pipework.and_then(|p| p.into());
            let mut heat_sources: IndexMap<String, PositionedHeatSource> = Default::default();
            let mut heat_source_for_diverter: Option<Arc<Mutex<HeatSource>>> = Default::default();
            for (name, hs) in &heat_source {
                let heater_position = hs.heater_position();
                let thermostat_position = hs.thermostat_position();
                let (heat_source, energy_supply_conn_name) = heat_source_from_input(
                    name.as_str(),
                    hs,
                    setpoint_temp,
                    wet_heat_sources,
                    simulation_time,
                    controls,
                    energy_supplies,
                    cold_water_sources,
                    external_conditions.clone(),
                );
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
                .and_then(|ctrl| controls.get_with_string(ctrl.as_str()).cloned());
            let storage_tank = Arc::new(Mutex::new(StorageTank::new(
                volume,
                daily_losses,
                min_temp,
                setpoint_temp,
                cold_water_source,
                simulation_time.step_in_hours(),
                heat_sources,
                pipework,
                Some(
                    EnergySupply::connection(energy_supplies.unmet_demand.clone(), &source_name)
                        .unwrap(),
                ),
                ctrl_hold_at_setpoint,
                WATER.clone(),
            )));
            for (heat_source_name, hs) in heat_source {
                let energy_supply_type = hs.energy_supply_type();
                if let Some(diverter) = diverter_types.get_for_supply_type(energy_supply_type) {
                    if diverter.storage_tank.matches(&source_name)
                        && diverter.heat_source.matches(&heat_source_name)
                    {
                        let energy_supply = energy_supplies.ensured_get_for_type(
                            hs.energy_supply_type(),
                            simulation_time.total_steps(),
                        );
                        let immersion_heater = heat_source_for_diverter
                            .clone()
                            .expect("More than one heat source was expected to be present")
                            .lock()
                            .as_immersion_heater();
                        if let Some(im) = immersion_heater {
                            let pv_diverter =
                                PVDiverter::new(storage_tank.clone(), im, heat_source_name);
                            energy_supply
                                .lock()
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
            power,
            efficiency,
            cold_water_source: cold_water_source_type,
            energy_supply,
        } => {
            let energy_supply =
                energy_supplies.ensured_get_for_type(energy_supply, simulation_time.total_steps());
            let energy_supply_conn_name = source_name;
            energy_supply_conn_names.push(energy_supply_conn_name.clone());
            let energy_supply_conn =
                EnergySupply::connection(energy_supply.clone(), &energy_supply_conn_name).unwrap();
            let cold_water_source =
                cold_water_source_for_type(cold_water_source_type, cold_water_sources);
            HotWaterSource::PointOfUse(PointOfUse::new(
                power,
                efficiency,
                energy_supply_conn,
                cold_water_source,
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
    (hot_water_source, energy_supply_conn_names)
}

fn cold_water_source_for_type(
    cold_water_source_type: ColdWaterSourceType,
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
    heat_system_names_for_zone: Vec<&str>,
) -> (
    HashMap<String, Arc<Mutex<SpaceHeatSystem>>>,
    HashMap<String, String>,
) {
    let mut energy_conn_names_for_systems: HashMap<String, String> = Default::default();
    let space_heat_systems = input
        .iter()
        .filter(|(system_name, _)| heat_system_names_for_zone.contains(&system_name.as_str()))
        .map(|(system_name, space_heat_system_details)| {
            (
                (*system_name).clone(),
                Arc::new(Mutex::new(match space_heat_system_details {
                    SpaceHeatSystemDetails::InstantElectricHeater {
                        rated_power,
                        control,
                        frac_convective,
                        energy_supply,
                        ..
                    } => {
                        let energy_supply = energy_supplies.ensured_get_for_type(*energy_supply, simulation_time.total_steps());
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
                                .and_then(|ctrl| controls.get_with_string(ctrl).map(|c| (*c).clone())),
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
                                SpaceHeatSystem::WarmAir(HeatPump::create_service_space_heating_warm_air(heat_pump.clone(), energy_supply_conn_name, control
                                    .as_ref()
                                    .and_then(|ctrl| controls.get_with_string(ctrl).map(|c| (*c).clone())).expect("A control object was expected for a heat pump warm air system"), *frac_convective).unwrap())
                            }
                            _ => panic!("The heat source referenced by details about warm air space heating with the name '{heat_source_name}' was expected to be a heat pump."),
                        }
                    }
                })),
            )
        })
        .collect::<HashMap<_, _>>();
    (space_heat_systems, energy_conn_names_for_systems)
}

fn space_cool_systems_from_input(
    input: &SpaceCoolSystemInput,
    cool_system_names_for_zone: Vec<&str>,
    controls: &Controls,
    energy_supplies: &mut EnergySupplies,
    simulation_time_iterator: &SimulationTimeIterator,
) -> HashMap<String, AirConditioning> {
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
                .ensured_get_for_type(*energy_supply, simulation_time_iterator.total_steps());
            let energy_supply_conn_name = system_name;
            let energy_supply_conn =
                EnergySupply::connection(energy_supply, energy_supply_conn_name).unwrap();
            let control = control
                .as_ref()
                .and_then(|ctrl| controls.get_with_string(ctrl).map(|c| (*c).clone()));

            (
                (*energy_supply_conn_name).clone(),
                AirConditioning::new(
                    *cooling_capacity,
                    *efficiency,
                    *frac_convective,
                    energy_supply_conn,
                    simulation_time_iterator.step_in_hours(),
                    control,
                ),
            )
        })
        .collect::<HashMap<_, _>>()
}

fn on_site_generation_from_input(
    input: &OnSiteGeneration,
    energy_supplies: &mut EnergySupplies,
    external_conditions: Arc<ExternalConditions>,
    simulation_time_iterator: &SimulationTimeIterator,
) -> HashMap<String, PhotovoltaicSystem> {
    input
        .iter()
        .map(|(name, generation_details)| {
            ((*name).clone(), {
                let OnSiteGenerationDetails {
                    peak_power,
                    ventilation_strategy,
                    pitch,
                    orientation,
                    base_height,
                    height,
                    width,
                    energy_supply,
                    ..
                } = generation_details;
                let energy_supply = energy_supplies
                    .ensured_get_for_type(*energy_supply, simulation_time_iterator.total_steps());
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
                )
            })
        })
        .collect::<HashMap<_, _>>()
}
