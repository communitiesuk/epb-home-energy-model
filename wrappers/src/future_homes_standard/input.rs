use anyhow::{anyhow, bail};
use hem::core::schedule::NumericSchedule;
use hem::input::{
    ApplianceGainsEvent, BuildingElement, ColdWaterSourceInput, ExternalConditionsInput,
    HeatSourceWetDetails, HeatingControlType, Input, JsonAccessResult, ReducedInputForCalcHtcHlp,
    SmartApplianceBattery, SpaceHeatSystemHeatSource, WasteWaterHeatRecovery, WaterDistribution,
    WaterHeatingEvent, WaterPipework, json_error,
};
use hem::simulation_time::SimulationTime;
use indexmap::IndexMap;
use itertools::Itertools;
use jsonschema::{BasicOutput, Validator};
use serde_json::{Map, Value as JsonValue, json};
use serde_valid::json::ToJsonString;
use std::collections::HashSet;
use std::io::{BufReader, Read};
use std::sync::LazyLock;

static FHS_SCHEMA_VALIDATOR: LazyLock<Validator> = LazyLock::new(|| {
    let schema =
        serde_json::from_str(include_str!("../../../schemas/input_fhs.schema.json")).unwrap();
    jsonschema::validator_for(&schema).unwrap()
});

pub fn ingest_for_processing(json: impl Read) -> Result<InputForProcessing, anyhow::Error> {
    InputForProcessing::init_with_json(json)
}

#[derive(Clone, Debug)]
pub struct InputForProcessing {
    pub(crate) input: JsonValue,
}

/// This type makes methods available for restricted access by wrappers,
/// in order to work towards a reasonable API for wrappers to interact with inputs rather than
/// the more brittle approach of allowing full access to the input data structure.
/// If the full access is encapsulated within methods here, it becomes possible to update the
/// underlying structure without breaking wrappers.
impl InputForProcessing {
    pub fn init_with_json(json: impl Read) -> Result<Self, anyhow::Error> {
        let input_for_processing = Self::init_with_json_skip_validation(json)?;

        let validator = &FHS_SCHEMA_VALIDATOR;

        if let BasicOutput::Invalid(errors) = validator.apply(&input_for_processing.input).basic() {
            bail!(
                "Invalid JSON against the FHS schema: {}",
                serde_json::to_value(errors)?.to_json_string_pretty()?
            ); // TODO build this handling logic out
        }

        Ok(input_for_processing)
    }

    pub(crate) fn init_with_json_skip_validation(json: impl Read) -> Result<Self, anyhow::Error> {
        let reader = BufReader::new(json);

        let input: JsonValue = serde_json::from_reader(reader)?;

        Ok(Self { input })
    }

    pub fn as_input(&self) -> anyhow::Result<Input> {
        serde_json::from_value(self.input.to_owned()).map_err(|err| anyhow!(err))
    }

    pub(crate) fn as_input_for_calc_htc_hlp(&self) -> anyhow::Result<ReducedInputForCalcHtcHlp> {
        serde_json::from_value(self.input.to_owned()).map_err(|err| anyhow!(err))
    }

    pub fn finalize(self) -> anyhow::Result<Input> {
        // NB. this _might_ in time be a good point to perform a validation against the core schema - or it might not
        // if let BasicOutput::Invalid(errors) =
        //     CORE_INCLUDING_FHS_VALIDATOR.apply(&self.input).basic()
        // {
        //     bail!(
        //         "Wrapper formed invalid JSON for the core schema: {}",
        //         serde_json::to_value(errors)?.to_json_string_pretty()?
        //     );
        // }
        serde_json::from_value(self.input).map_err(|err| anyhow!(err))
    }

    fn root(&self) -> JsonAccessResult<&Map<std::string::String, JsonValue>> {
        self.input
            .as_object()
            .ok_or(json_error("Root document is not an object"))
    }

    fn root_mut(&mut self) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.input
            .as_object_mut()
            .ok_or(json_error("Root document is not an object"))
    }

    fn set_on_root_key(&mut self, root_key: &str, value: JsonValue) -> JsonAccessResult<&mut Self> {
        self.root_mut()?.insert(root_key.into(), value);

        Ok(self)
    }

    fn remove_root_key(&mut self, root_key: &str) -> JsonAccessResult<&mut Self> {
        self.root_mut()?.shift_remove(root_key);

        Ok(self)
    }

    fn root_object(
        &self,
        root_key: &str,
    ) -> JsonAccessResult<&Map<std::string::String, JsonValue>> {
        self.root()?
            .get(root_key)
            .ok_or(json_error(format!("No {root_key} node found")))?
            .as_object()
            .ok_or(json_error(format!("{root_key} node was not an object")))
    }

    fn root_object_mut(
        &mut self,
        root_key: &str,
    ) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.root_mut()?
            .get_mut(root_key)
            .ok_or(json_error(format!("No {root_key} node found")))?
            .as_object_mut()
            .ok_or(json_error(format!("{root_key} node was not an object")))
    }

    /// Uses entry API to ensure that root key is created if it does not already exist.
    fn root_object_entry_mut(
        &mut self,
        root_key: &str,
    ) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.root_mut()?
            .entry(root_key)
            .or_insert(json!({}))
            .as_object_mut()
            .ok_or(json_error(format!("{root_key} node was not an object")))
    }

    fn optional_root_object(
        &self,
        root_key: &str,
    ) -> JsonAccessResult<Option<&Map<std::string::String, JsonValue>>> {
        Ok(self.root()?.get(root_key).and_then(|v| v.as_object()))
    }

    pub fn set_simulation_time(
        &mut self,
        simulation_time: SimulationTime,
    ) -> anyhow::Result<&mut Self> {
        self.set_on_root_key("SimulationTime", serde_json::to_value(simulation_time)?)
            .map_err(Into::into)
    }

    pub fn set_temp_internal_air_static_calcs(
        &mut self,
        temp_internal_air_static_calcs: Option<f64>,
    ) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key(
            "temp_internal_air_static_calcs",
            temp_internal_air_static_calcs.into(),
        )
    }

    pub(crate) fn merge_external_conditions_data(
        &mut self,
        external_conditions_data: Option<ExternalConditionsInput>,
    ) -> anyhow::Result<()> {
        if let Some(external_conditions) = external_conditions_data {
            let shading_segments = self
                .root_object("ExternalConditions")?
                .get("shading_segments")
                .cloned()
                .unwrap_or(json!([]));
            let mut new_external_conditions = serde_json::to_value(external_conditions)?;
            let new_external_conditions_map = new_external_conditions.as_object_mut().ok_or(json_error("External conditions was not a JSON object when it was expected to be provided as one"))?;
            new_external_conditions_map.insert("shading_segments".into(), shading_segments);
            self.set_on_root_key("ExternalConditions", new_external_conditions)?;
        }

        Ok(())
    }

    pub fn reset_internal_gains(&mut self) -> JsonAccessResult<&Self> {
        self.root_mut()?.insert("InternalGains".into(), json!({}));

        Ok(self)
    }

    fn zone_node(&self) -> JsonAccessResult<&serde_json::Map<std::string::String, JsonValue>> {
        self.root_object("Zone")
    }

    fn zone_node_mut(
        &mut self,
    ) -> JsonAccessResult<&mut serde_json::Map<std::string::String, JsonValue>> {
        self.root_object_mut("Zone")
    }

    fn specific_zone(
        &self,
        zone_key: &str,
    ) -> JsonAccessResult<&serde_json::Map<std::string::String, JsonValue>> {
        self.zone_node()?
            .get(zone_key)
            .ok_or(json_error(format!("Zone key {zone_key} did not exist")))?
            .as_object()
            .ok_or(json_error("Zone node was not an object"))
    }

    fn specific_zone_mut(
        &mut self,
        zone_key: &str,
    ) -> JsonAccessResult<&mut serde_json::Map<std::string::String, JsonValue>> {
        self.zone_node_mut()?
            .get_mut(zone_key)
            .ok_or(json_error(format!("Zone key {zone_key} did not exist")))?
            .as_object_mut()
            .ok_or(json_error("Zone node was not an object"))
    }

    pub fn total_zone_area(&self) -> JsonAccessResult<f64> {
        self.zone_node()?
            .values()
            .map(|z| {
                z.get("area")
                    .ok_or(json_error("Area field not found on zone"))?
                    .as_f64()
                    .ok_or(json_error("Area field not a number"))
            })
            .sum::<JsonAccessResult<f64>>()
    }

    pub fn total_zone_volume(&self) -> JsonAccessResult<f64> {
        Ok(self
            .zone_node()?
            .values()
            .map(|z| {
                z.get("volume")
                    .ok_or(json_error("Volume field not found on zone"))?
                    .as_number()
                    .ok_or(json_error("Volume field not a number"))?
                    .as_f64()
                    .ok_or(json_error("Volume field not a number"))
            })
            .collect::<JsonAccessResult<Vec<_>>>()?
            .into_iter()
            .sum::<f64>())
    }

    pub fn area_for_zone(&self, zone: &str) -> anyhow::Result<f64> {
        Ok(self
            .zone_node()?
            .get(zone)
            .ok_or(anyhow!("Used zone key for a zone that does not exist"))?
            .get("area")
            .ok_or(json_error("Area not found on zone"))?
            .as_number()
            .ok_or(json_error("Area on zone was not a number"))?
            .as_f64()
            .ok_or(json_error("Area number could not be read as a number"))?)
    }

    #[cfg(test)]
    pub(crate) fn all_thermal_bridgings(&self) -> JsonAccessResult<Vec<&JsonValue>> {
        Ok(self
            .zone_node()?
            .values()
            .flat_map(|z| z.get("ThermalBridging"))
            .collect::<Vec<_>>())
    }

    pub(crate) fn all_thermal_bridging_elements(
        &mut self,
    ) -> JsonAccessResult<Vec<&mut Map<std::string::String, JsonValue>>> {
        let zones = self.zone_node_mut()?;
        let mut result = Vec::new();
        for zone in zones.values_mut() {
            if let Some(thermal_bridging) = zone.get_mut("ThermalBridging") {
                if let Some(obj) = thermal_bridging.as_object_mut() {
                    result.push(obj);
                }
            }
        }

        Ok(result)
    }

    pub fn number_of_bedrooms(&self) -> JsonAccessResult<Option<usize>> {
        match self.input.get("NumberOfBedrooms") {
            None => Ok(None),
            Some(JsonValue::Number(n)) => Ok(Some(
                n.as_u64()
                    .ok_or(json_error("NumberOfBedrooms not a positive integer"))?
                    as usize,
            )),
            Some(_) => Err(json_error("NumberOfBedrooms not a number")),
        }
    }

    pub(crate) fn number_of_wet_rooms(&self) -> JsonAccessResult<Option<usize>> {
        match self.input.get("NumberOfWetRooms") {
            None => Ok(None),
            Some(JsonValue::Number(n)) => Ok(Some(
                n.as_u64()
                    .ok_or(json_error("NumberOfWetRooms not a positive integer"))?
                    as usize,
            )),
            Some(_) => Err(json_error("NumberOfWetRooms not a number")),
        }
    }

    fn internal_gains_mut(&mut self) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.root_object_entry_mut("InternalGains")
    }

    pub fn set_metabolic_gains(
        &mut self,
        start_day: u32,
        time_series_step: f64,
        schedule_json: JsonValue,
    ) -> anyhow::Result<&Self> {
        self.internal_gains_mut()?.insert(
            "metabolic gains".into(),
            json!({
                "start_day": start_day,
                "time_series_step": time_series_step,
                "schedule": schedule_json,
            }),
        );

        Ok(self)
    }

    pub fn set_evaporative_losses(
        &mut self,
        start_day: u32,
        time_series_step: f64,
        schedule_json: JsonValue,
    ) -> anyhow::Result<&Self> {
        self.internal_gains_mut()?.insert(
            "EvaporativeLosses".into(),
            json!({
                "start_day": start_day,
                "time_series_step": time_series_step,
                "schedule": schedule_json,
            }),
        );

        Ok(self)
    }

    pub fn set_cold_water_losses(
        &mut self,
        start_day: u32,
        time_series_step: f64,
        schedule_json: JsonValue,
    ) -> anyhow::Result<&Self> {
        self.internal_gains_mut()?.insert(
            "ColdWaterLosses".into(),
            json!({
                "start_day": start_day,
                "time_series_step": time_series_step,
                "schedule": schedule_json,
            }),
        );

        Ok(self)
    }

    pub fn heating_control_type(&self) -> JsonAccessResult<Option<HeatingControlType>> {
        self.root()?
            .get("HeatingControlType")
            .map(
                |node| match serde_json::from_value::<HeatingControlType>(node.to_owned()) {
                    Ok(t) => Ok(t),
                    Err(_) => Err(json_error(
                        "Could not parse HeatingControlType into a known value",
                    )),
                },
            )
            .transpose()
    }

    pub fn set_heating_control_type(
        &mut self,
        heating_control_type_value: JsonValue,
    ) -> anyhow::Result<&mut Self> {
        self.set_on_root_key("HeatingControlType", heating_control_type_value)
            .map_err(Into::into)
    }

    pub fn add_control(
        &mut self,
        control_key: &str,
        control_json: JsonValue,
    ) -> JsonAccessResult<&Self> {
        self.root_object_entry_mut("Control")?
            .insert(control_key.into(), control_json);

        Ok(self)
    }

    pub fn remove_all_smart_appliance_controls(&mut self) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key("SmartApplianceControls", json!({}))
    }

    fn smart_appliance_controls_mut(
        &mut self,
    ) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.root_object_entry_mut("SmartApplianceControls")
    }

    pub fn add_smart_appliance_control(
        &mut self,
        smart_control_name: &str,
        control: JsonValue,
    ) -> JsonAccessResult<&Self> {
        self.smart_appliance_controls_mut()?
            .insert(smart_control_name.into(), control);

        Ok(self)
    }

    pub fn set_non_appliance_demand_24hr_on_smart_appliance_control(
        &mut self,
        smart_control_name: &str,
        non_appliance_demand_24hr_input: IndexMap<smartstring::alias::String, Vec<f64>>,
    ) -> JsonAccessResult<&Self> {
        if let Some(ref mut control) = self
            .smart_appliance_controls_mut()?
            .get_mut(smart_control_name)
            .and_then(|v| v.as_object_mut())
        {
            control.insert(
                "non_appliance_demand_24hr".into(),
                json!(non_appliance_demand_24hr_input),
            );
        }

        Ok(self)
    }

    pub fn set_battery24hr_on_smart_appliance_control(
        &mut self,
        smart_control_name: &str,
        battery24hr_input: SmartApplianceBattery,
    ) -> JsonAccessResult<&Self> {
        if let Some(ref mut control) = self
            .smart_appliance_controls_mut()?
            .get_mut(smart_control_name)
            .and_then(|v| v.as_object_mut())
        {
            control.insert("battery24hr".into(), json!(battery24hr_input));
        }

        Ok(self)
    }

    pub fn zone_keys(&self) -> JsonAccessResult<Vec<smartstring::alias::String>> {
        Ok(self
            .zone_node()?
            .keys()
            .map(smartstring::alias::String::from)
            .collect())
    }

    #[cfg(test)]
    pub(crate) fn all_init_temp_setpoints(&self) -> JsonAccessResult<Vec<Option<f64>>> {
        Ok(self
            .zone_node()?
            .values()
            .map(|zone| zone.get("temp_setpnt_init").and_then(|t| t.as_f64()))
            .collect())
    }

    pub fn set_init_temp_setpoint_for_zone(
        &mut self,
        zone: &str,
        temperature: f64,
    ) -> JsonAccessResult<&Self> {
        self.specific_zone_mut(zone)?
            .insert("temp_setpnt_init".into(), json!(temperature));
        Ok(self)
    }

    pub fn space_heat_control_for_zone(
        &self,
        zone: &str,
    ) -> anyhow::Result<Option<smartstring::alias::String>> {
        Ok(self
            .specific_zone(zone)?
            .get("SpaceHeatControl")
            .and_then(|field| field.as_str())
            .map(smartstring::alias::String::from))
    }

    pub fn space_heat_system_for_zone(
        &self,
        zone: &str,
    ) -> JsonAccessResult<Vec<smartstring::alias::String>> {
        Ok(match self.specific_zone(zone)?.get("SpaceHeatSystem") {
            Some(JsonValue::String(system)) => vec![smartstring::alias::String::from(system)],
            Some(JsonValue::Array(systems)) => systems
                .iter()
                .map(|system| {
                    Ok(smartstring::alias::String::from(system.as_str().ok_or(
                        json_error("Space heat system list contained a non-string"),
                    )?))
                })
                .collect::<Result<Vec<_>, _>>()?,
            _ => vec![],
        })
    }

    pub fn set_space_heat_system_for_zone(
        &mut self,
        zone: &str,
        system_name: &str,
    ) -> anyhow::Result<&Self> {
        let zone = self.specific_zone_mut(zone)?;
        zone.insert("SpaceHeatSystem".into(), system_name.into());

        Ok(self)
    }

    pub fn space_cool_system_for_zone(
        &self,
        zone: &str,
    ) -> JsonAccessResult<Vec<smartstring::alias::String>> {
        Ok(match self.specific_zone(zone)?.get("SpaceCoolSystem") {
            Some(JsonValue::String(system)) => vec![smartstring::alias::String::from(system)],
            Some(JsonValue::Array(systems)) => systems
                .iter()
                .map(|s| {
                    Ok(smartstring::alias::String::from(s.as_str().ok_or(
                        json_error("SpaceCoolSystem list contained non-strings"),
                    )?))
                })
                .collect::<Result<_, _>>()?,
            _ => vec![],
        })
    }

    pub fn set_space_cool_system_for_zone(
        &mut self,
        zone: &str,
        system_name: &str,
    ) -> anyhow::Result<&Self> {
        let zone = self.specific_zone_mut(zone)?;
        zone.insert("SpaceCoolSystem".into(), system_name.into());

        Ok(self)
    }

    #[cfg(test)]
    pub(crate) fn lighting_efficacy_for_zone(&self, zone: &str) -> JsonAccessResult<Option<f64>> {
        Ok(self
            .specific_zone(zone)?
            .get("Lighting")
            .and_then(|v| v.as_object())
            .and_then(|lighting| lighting.get("efficacy"))
            .and_then(|efficacy| efficacy.as_f64()))
    }

    pub fn set_lighting_efficacy_for_all_zones(
        &mut self,
        efficacy: f64,
    ) -> JsonAccessResult<&Self> {
        for lighting in self
            .zone_node_mut()?
            .values_mut()
            .filter_map(|zone| zone.get_mut("Lighting").and_then(|l| l.as_object_mut()))
        {
            lighting.insert("efficacy".into(), efficacy.into());
        }

        Ok(self)
    }

    pub fn all_zones_have_bulbs(&self) -> JsonAccessResult<bool> {
        Ok(self.zone_node()?.values().all(|zone| {
            zone.get("Lighting")
                .and_then(|l| l.as_object())
                .and_then(|l| l.get("bulbs"))
                .is_some_and(|bulbs| bulbs.is_object())
        }))
    }

    pub fn light_bulbs_for_each_zone(
        &self,
    ) -> JsonAccessResult<IndexMap<smartstring::alias::String, Map<std::string::String, JsonValue>>>
    {
        Ok(self
            .zone_node()?
            .iter()
            .map(|(zone_name, zone)| {
                let bulbs = zone
                    .get("Lighting")
                    .and_then(|lighting| lighting.get("bulbs"))
                    .and_then(|bulbs| bulbs.as_object());
                (
                    smartstring::alias::String::from(zone_name),
                    bulbs.map(ToOwned::to_owned).unwrap_or_default(),
                )
            })
            .collect())
    }

    pub fn set_control_window_opening_for_zone(
        &mut self,
        zone: &str,
        opening_type: Option<&str>,
    ) -> anyhow::Result<&Self> {
        self.specific_zone_mut(zone)?
            .insert("Control_WindowOpening".into(), json!(opening_type));

        Ok(self)
    }

    pub fn set_control_string_for_space_heat_system(
        &mut self,
        space_heat_system: &str,
        control_string: &str,
    ) -> anyhow::Result<&Self> {
        self.root_object_mut("SpaceHeatSystem")?
            .get_mut(space_heat_system)
            .ok_or(anyhow!(
                "There is no provided space heat system with the name '{space_heat_system}'"
            ))?
            .as_object_mut()
            .ok_or(json_error("Space heat system was not an object"))?
            .insert("Control".into(), json!(control_string));

        Ok(self)
    }

    pub fn set_control_string_for_space_cool_system(
        &mut self,
        space_cool_system: &str,
        control_string: &str,
    ) -> anyhow::Result<&Self> {
        self.root_object_mut("SpaceCoolSystem")?
            .get_mut(space_cool_system)
            .ok_or(anyhow!(
                "There is no provided space cool system with the name '{space_cool_system}'"
            ))?
            .as_object_mut()
            .ok_or(json_error("Space cool system was not an object"))?
            .insert("Control".into(), json!(control_string));

        Ok(self)
    }

    pub fn has_named_smart_appliance_control(&self, name: &str) -> JsonAccessResult<bool> {
        Ok(self
            .optional_root_object("SmartApplianceControls")?
            .is_some_and(|controls| controls.contains_key(name)))
    }

    pub(crate) fn set_efficiency_for_all_space_cool_systems(
        &mut self,
        efficiency: f64,
    ) -> JsonAccessResult<()> {
        let systems = self.root_object_entry_mut("SpaceCoolSystem")?;
        for system in systems.values_mut().flat_map(|s| s.as_object_mut()) {
            system.insert("efficiency".into(), json!(efficiency));
        }

        Ok(())
    }

    pub(crate) fn set_frac_convective_for_all_space_cool_systems(
        &mut self,
        frac_convective: f64,
    ) -> JsonAccessResult<()> {
        let systems = self.root_object_entry_mut("SpaceCoolSystem")?;
        for system in systems.values_mut().flat_map(|s| s.as_object_mut()) {
            system.insert("frac_convective".into(), json!(frac_convective));
        }

        Ok(())
    }

    pub(crate) fn set_energy_supply_for_all_space_cool_systems(
        &mut self,
        energy_supply_name: &str,
    ) -> JsonAccessResult<()> {
        let systems = self.root_object_entry_mut("SpaceCoolSystem")?;
        for system in systems.values_mut().flat_map(|s| s.as_object_mut()) {
            system.insert("EnergySupply".into(), json!(energy_supply_name));
        }

        Ok(())
    }

    #[cfg(test)]
    pub(crate) fn space_cool_system(
        &self,
    ) -> JsonAccessResult<Option<&Map<std::string::String, JsonValue>>> {
        self.optional_root_object("SpaceCoolSystem")
    }

    pub(crate) fn space_heat_system_keys(
        &self,
    ) -> JsonAccessResult<Vec<smartstring::alias::String>> {
        Ok(match self.optional_root_object("SpaceHeatSystem")? {
            Some(space_heat_system) => space_heat_system
                .keys()
                .map(smartstring::alias::String::from)
                .collect(),
            None => vec![],
        })
    }

    pub fn temperature_setback_for_space_heat_system(
        &self,
        space_heat_system: &str,
    ) -> JsonAccessResult<Option<f64>> {
        let space_heat_systems = self.optional_root_object("SpaceHeatSystem")?;
        let space_heat_systems = match space_heat_systems {
            Some(ref space_heat_systems) => space_heat_systems,
            None => return Ok(None),
        };
        let space_heat_system = space_heat_systems.get(space_heat_system);

        Ok(space_heat_system.and_then(|space_heat_system| {
            space_heat_system.as_object().and_then(|space_heat_system| {
                space_heat_system
                    .get("temp_setback")
                    .and_then(|temp_setback| temp_setback.as_f64())
            })
        }))
    }

    pub fn temperature_setback_for_space_cool_system(
        &self,
        space_cool_system: &str,
    ) -> JsonAccessResult<Option<f64>> {
        let space_cool_systems = self.optional_root_object("SpaceCoolSystem")?;
        let space_cool_systems = match space_cool_systems {
            Some(ref space_heat_systems) => space_heat_systems,
            None => return Ok(None),
        };
        let space_cool_system = space_cool_systems.get(space_cool_system);

        Ok(space_cool_system.and_then(|space_cool_system| {
            space_cool_system.as_object().and_then(|space_cool_system| {
                space_cool_system
                    .get("temp_setback")
                    .and_then(|temp_setback| temp_setback.as_f64())
            })
        }))
    }

    pub fn advanced_start_for_space_cool_system(
        &self,
        space_cool_system: &str,
    ) -> JsonAccessResult<Option<f64>> {
        let space_cool_systems = self.optional_root_object("SpaceCoolSystem")?;
        let space_cool_systems = match space_cool_systems {
            Some(ref space_heat_systems) => space_heat_systems,
            None => return Ok(None),
        };
        let space_cool_system = space_cool_systems.get(space_cool_system);

        Ok(space_cool_system.and_then(|space_cool_system| {
            space_cool_system.as_object().and_then(|space_cool_system| {
                space_cool_system
                    .get("advanced_start")
                    .and_then(|temp_setback| temp_setback.as_f64())
            })
        }))
    }

    pub fn advanced_start_for_space_heat_system(
        &self,
        space_heat_system: &str,
    ) -> JsonAccessResult<Option<f64>> {
        let space_heat_systems = self.optional_root_object("SpaceHeatSystem")?;
        let space_heat_systems = match space_heat_systems {
            Some(ref space_heat_systems) => space_heat_systems,
            None => return Ok(None),
        };
        let space_heat_system = space_heat_systems.get(space_heat_system);

        Ok(space_heat_system.and_then(|space_heat_system| {
            space_heat_system.as_object().and_then(|space_heat_system| {
                space_heat_system
                    .get("advanced_start")
                    .and_then(|temp_setback| temp_setback.as_f64())
            })
        }))
    }

    pub(crate) fn set_advance_start_for_space_heat_system(
        &mut self,
        space_heat_system: &str,
        new_advanced_start: f64,
    ) -> anyhow::Result<&Self> {
        self.root_object_entry_mut("SpaceHeatSystem")?
            .get_mut(space_heat_system)
            .ok_or(anyhow!(
                "There is no provided space heat system with the name '{space_heat_system}'"
            ))?
            .as_object_mut()
            .ok_or(json_error(
                "The indicated space heat system was not an object",
            ))?
            .insert("advanced_start".into(), new_advanced_start.into());
        Ok(self)
    }

    pub(crate) fn set_temperature_setback_for_space_heat_systems(
        &mut self,
        new_temperature_setback: Option<f64>,
    ) -> anyhow::Result<()> {
        self.root_object_entry_mut("SpaceHeatSystem")?
            .values_mut()
            .flat_map(|system| system.as_object_mut())
            .for_each(|system_details| {
                system_details.insert("temp_setback".into(), new_temperature_setback.into());
            });
        Ok(())
    }

    #[cfg(test)]
    pub(crate) fn heat_source_for_space_heat_system(
        &self,
        space_heat_system: &str,
    ) -> JsonAccessResult<Option<&JsonValue>> {
        let space_heat_systems = self.optional_root_object("SpaceHeatSystem")?;
        let space_heat_systems = match space_heat_systems {
            Some(ref space_heat_systems) => space_heat_systems,
            None => return Ok(None),
        };
        let space_heat_system = space_heat_systems.get(space_heat_system);

        Ok(space_heat_system.and_then(|space_heat_system| {
            space_heat_system
                .as_object()
                .and_then(|space_heat_system| space_heat_system.get("HeatSource"))
        }))
    }

    pub(crate) fn set_heat_source_for_all_space_heat_systems(
        &mut self,
        heat_source: SpaceHeatSystemHeatSource,
    ) -> anyhow::Result<()> {
        self.root_object_entry_mut("SpaceHeatSystem")?
            .values_mut()
            .flat_map(|system| system.as_object_mut())
            .for_each(|system_details| {
                system_details.insert("HeatSource".into(), json!(heat_source));
            });
        Ok(())
    }

    pub(crate) fn set_hot_water_source(
        &mut self,
        hot_water_source: JsonValue,
    ) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key("HotWaterSource", hot_water_source)
    }

    pub fn hot_water_source(&self) -> JsonAccessResult<&Map<std::string::String, JsonValue>> {
        self.root_object("HotWaterSource")
    }

    pub fn hot_water_source_mut(
        &mut self,
    ) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.root_object_mut("HotWaterSource")
    }

    pub fn names_of_energy_supplies_with_diverters(
        &self,
    ) -> JsonAccessResult<Vec<smartstring::alias::String>> {
        Ok(self
            .root_object("EnergySupply")?
            .iter()
            .filter_map(|(energy_supply_name, energy_supply)| {
                energy_supply.as_object().and_then(|energy_supply| {
                    energy_supply
                        .get("diverter")
                        .map(|_| smartstring::alias::String::from(energy_supply_name))
                })
            })
            .collect_vec())
    }

    pub fn set_control_max_name_for_energy_supply_diverter(
        &mut self,
        energy_supply_name: &str,
        control_max_name: &str,
    ) -> JsonAccessResult<&Self> {
        self.root_object_mut("EnergySupply")?
            .get_mut(energy_supply_name)
            .ok_or(json_error(format!(
                "There is no provided energy supply with the name '{energy_supply_name}'"
            )))?
            .as_object_mut()
            .ok_or(json_error("The indicated energy supply was not an object"))?
            .get_mut(energy_supply_name)
            .map(|energy_supply| {
                energy_supply.get("diverter").as_mut().map(|diverter| {
                    diverter.get("Controlmax").replace(&json!(control_max_name));
                })
            });
        Ok(self)
    }

    pub fn set_lighting_gains(&mut self, gains_details: JsonValue) -> JsonAccessResult<&Self> {
        self.set_gains_for_field("lighting", gains_details)
    }

    pub fn set_topup_gains(&mut self, gains_details: JsonValue) -> JsonAccessResult<&Self> {
        self.set_gains_for_field("topup", gains_details)
    }

    pub fn set_gains_for_field(
        &mut self,
        field: impl Into<std::string::String>,
        gains_details: JsonValue,
    ) -> JsonAccessResult<&Self> {
        self.root_object_entry_mut("ApplianceGains")?
            .insert(field.into(), gains_details);

        Ok(self)
    }

    pub fn energy_supply_type_for_appliance_gains_field(
        &self,
        field: &str,
    ) -> Option<smartstring::alias::String> {
        self.root_object("ApplianceGains")
            .ok()
            .and_then(|appliance_gains| appliance_gains.get(field))
            .and_then(|details| {
                details
                    .get("EnergySupply")
                    .and_then(|energy_supply| energy_supply.as_str())
                    .map(smartstring::alias::String::from)
            })
    }

    pub fn clear_appliance_gains(&mut self) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key("ApplianceGains", json!({}))
    }

    pub fn set_priority_for_gains_appliance(
        &mut self,
        priority: isize,
        appliance: &str,
    ) -> anyhow::Result<()> {
        self.root_object_entry_mut("ApplianceGains")?
            .get_mut(appliance)
            .ok_or_else(|| anyhow!("Encountered bad appliance gains reference {appliance:?}"))?
            .as_object_mut()
            .ok_or_else(|| anyhow!("Appliance gains reference was not a JSON object"))?
            .insert("priority".into(), json!(priority));

        Ok(())
    }

    pub fn fuel_type_for_energy_supply_reference(
        &self,
        reference: &str,
    ) -> anyhow::Result<smartstring::alias::String> {
        Ok(self
            .root_object("EnergySupply")?
            .get(reference)
            .ok_or(anyhow!(
                "Energy supply with reference '{reference}' could not be found"
            ))?
            .get("fuel")
            .ok_or(json_error("Energy supply object did not have a fuel field"))?
            .as_str()
            .ok_or(json_error(
                "Energy supply fuel field expected to have a string value",
            ))?
            .into())
    }

    pub fn shower_flowrates(&self) -> JsonAccessResult<IndexMap<smartstring::alias::String, f64>> {
        let showers = match self
            .hot_water_demand()?
            .get("Shower")
            .and_then(|s| s.as_object())
        {
            None => return Ok(Default::default()),
            Some(showers) => showers,
        };

        Ok(showers
            .iter()
            .filter_map(|(name, shower)| {
                shower
                    .get("flowrate")
                    .and_then(|s| s.as_f64())
                    .map(|flow_rate| (smartstring::alias::String::from(name), flow_rate))
            })
            .collect())
    }

    pub fn reset_water_heating_events(&mut self) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key("Events", json!({}))
    }

    pub(crate) fn showers(&self) -> JsonAccessResult<Option<&Map<std::string::String, JsonValue>>> {
        self.hot_water_demand()?
            .get("Shower")
            .map(|showers| {
                showers
                    .as_object()
                    .ok_or(json_error("Shower was not an object"))
            })
            .transpose()
    }

    pub fn shower_keys(&self) -> JsonAccessResult<Vec<smartstring::alias::String>> {
        Ok(self
            .hot_water_demand()?
            .get("Shower")
            .map_or(vec![], |showers| match showers.as_object() {
                Some(showers) => showers
                    .keys()
                    .map(smartstring::alias::String::from)
                    .collect(),
                None => vec![],
            }))
    }

    pub fn shower_name_refers_to_instant_electric(&self, name: &str) -> bool {
        self.hot_water_demand()
            .ok()
            .and_then(|demand| demand.get("Shower"))
            .and_then(|showers| showers.get(name))
            .and_then(|shower| shower.get("type"))
            .and_then(|shower_type| shower_type.as_str())
            .is_some_and(|shower_type| shower_type == "InstantElecShower")
    }

    pub(crate) fn register_wwhrs_name_on_mixer_shower(
        &mut self,
        wwhrs: &str,
    ) -> anyhow::Result<()> {
        let mixer_shower = self
            .hot_water_demand_mut()?
            .get_mut("Shower")
            .ok_or(json_error("Shower node not set on HotWaterDemand"))?
            .get_mut("mixer")
            .ok_or(json_error(
                "A mixer shower with the reference 'mixer' was expected to have been set",
            ))?
            .as_object_mut()
            .ok_or(json_error("Mixer shower was not a JSON object"))?;
        mixer_shower.insert("WWHRS".into(), json!(wwhrs));

        Ok(())
    }

    pub(crate) fn baths(&self) -> JsonAccessResult<Option<&Map<std::string::String, JsonValue>>> {
        Ok(self
            .hot_water_demand()?
            .get("Bath")
            .and_then(|baths| baths.as_object()))
    }

    pub fn bath_keys(&self) -> JsonAccessResult<Vec<smartstring::alias::String>> {
        Ok(self
            .hot_water_demand()?
            .get("Bath")
            .map_or(vec![], |baths| match baths.as_object() {
                Some(baths) => baths.keys().map(smartstring::alias::String::from).collect(),
                None => vec![],
            }))
    }

    pub fn size_for_bath_field(&self, field: &str) -> JsonAccessResult<Option<f64>> {
        Ok(self
            .hot_water_demand()?
            .get("Bath")
            .and_then(|baths| baths.as_object())
            .and_then(|bath| bath.get(field))
            .and_then(|bath| bath.get("size"))
            .and_then(|size| size.as_f64()))
    }

    pub fn flowrate_for_bath_field(&self, field: &str) -> JsonAccessResult<Option<f64>> {
        Ok(self
            .hot_water_demand()?
            .get("Bath")
            .and_then(|baths| baths.as_object())
            .and_then(|bath| bath.get(field))
            .and_then(|bath| bath.get("flowrate"))
            .and_then(|flowrate| flowrate.as_f64()))
    }

    pub(crate) fn other_water_uses(
        &self,
    ) -> JsonAccessResult<Option<&Map<std::string::String, JsonValue>>> {
        Ok(self
            .hot_water_demand()?
            .get("Other")
            .and_then(|other| other.as_object()))
    }

    pub fn other_water_use_keys(&self) -> JsonAccessResult<Vec<smartstring::alias::String>> {
        Ok(self
            .other_water_uses()?
            .map(|other| other.keys().map(smartstring::alias::String::from).collect())
            .unwrap_or(vec![]))
    }

    pub fn flow_rate_for_other_water_use_field(
        &self,
        field: &str,
    ) -> JsonAccessResult<Option<f64>> {
        Ok(self
            .hot_water_demand()?
            .get("Other")
            .and_then(|others| others.as_object())
            .and_then(|other| other.get(field))
            .and_then(|other| other.get("flowrate"))
            .and_then(|flowrate| flowrate.as_f64()))
    }

    pub fn set_other_water_use_details(
        &mut self,
        cold_water_source_type: &str,
        flowrate: f64,
    ) -> JsonAccessResult<()> {
        let other_details = json!({
            "flowrate": flowrate,
            "ColdWaterSource": cold_water_source_type,
        });

        let other_water_uses = self
            .hot_water_demand_mut()?
            .entry("Other")
            .or_insert(json!({}))
            .as_object_mut()
            .ok_or(json_error("Other water uses not provided as an object"))?;
        other_water_uses.insert("other".into(), other_details);

        Ok(())
    }

    pub(crate) fn water_distribution(&self) -> anyhow::Result<Option<WaterDistribution>> {
        Ok(self
            .root_object("HotWaterDemand")?
            .get("Distribution")
            .map(|node| serde_json::from_value::<WaterDistribution>(node.to_owned()))
            .transpose()?)
    }

    /// Override all the vol_hw_daily_average values on the heat pump hot water only heat sources.
    pub fn override_vol_hw_daily_average_on_heat_pumps(&mut self, vol_hw_daily_average: f64) {
        let heat_sources = match self
            .root_object_mut("HotWaterSource")
            .ok()
            .and_then(|hot_water_source| hot_water_source.get_mut("hw cylinder"))
            .and_then(|cylinder| cylinder.as_object_mut())
            .and_then(|cylinder| {
                cylinder
                    .get("type")
                    .and_then(|source_type| source_type.as_str())
                    .is_some_and(|source_type| source_type == "StorageTank")
                    .then_some(cylinder)
            }) {
            None => return,
            Some(heat_sources) => heat_sources,
        };

        for heat_source in heat_sources.values_mut().flat_map(|hs| hs.as_object_mut()) {
            if heat_source
                .get("type")
                .and_then(|heat_source_type| heat_source_type.as_str())
                .is_some_and(|heat_source_type| heat_source_type == "HeatPump_HWOnly")
            {
                heat_source.insert("vol_hw_daily_average".into(), json!(vol_hw_daily_average));
            }
        }
    }

    pub fn part_g_compliance(&self) -> JsonAccessResult<Option<bool>> {
        self.root()?
            .get("PartGcompliance")
            .map(|node| {
                node.as_bool()
                    .ok_or(json_error("Part G compliance was not passed as a boolean"))
            })
            .transpose()
    }

    pub fn set_part_g_compliance(&mut self, is_compliant: bool) -> JsonAccessResult<&Self> {
        self.root_mut()?
            .insert("PartGcompliance".into(), is_compliant.into());

        Ok(self)
    }

    pub fn add_water_heating_event(
        &mut self,
        event_type: &str,
        subtype_name: &str,
        event: JsonValue,
    ) -> JsonAccessResult<&Self> {
        let node_for_type = self
            .root_object_entry_mut("Events")?
            .entry(event_type)
            .or_insert(json!({}))
            .as_object_mut()
            .ok_or(json_error("Events node was not an object"))?;
        let node_for_subtype = node_for_type.entry(subtype_name).or_insert(json!([]));
        if let Some(events) = node_for_subtype.as_array_mut() {
            events.push(event);
        } else {
            return Err(json_error(format!(
                "Events node at '{event_type}' -> '{subtype_name}' was not an array"
            )));
        }

        Ok(self)
    }

    pub fn water_heating_events_of_types(
        &self,
        event_types: &[&str],
    ) -> JsonAccessResult<Vec<JsonValue>> {
        Ok(self
            .root_object("Events")?
            .iter()
            .filter(|(event_type, _)| event_types.contains(&&***event_type))
            .flat_map(|(_, events)| {
                events
                    .as_object()
                    .map(|events| events.values().filter_map(JsonValue::as_array))
            })
            .flatten()
            .flatten()
            .cloned()
            .collect_vec())
    }

    pub fn cold_water_source_has_header_tank(&self) -> JsonAccessResult<bool> {
        Ok(self
            .root_object("ColdWaterSource")?
            .contains_key("header tank"))
    }

    pub fn set_cold_water_source_by_key(
        &mut self,
        key: &str,
        source_details: JsonValue,
    ) -> JsonAccessResult<&Self> {
        self.root_object_mut("ColdWaterSource")?
            .insert(key.into(), source_details);

        Ok(self)
    }

    pub fn set_hot_water_cylinder(&mut self, source_value: JsonValue) -> JsonAccessResult<&Self> {
        let hot_water_source = self.root_object_mut("HotWaterSource")?;
        hot_water_source.insert("hw cylinder".into(), source_value);

        Ok(self)
    }

    fn hot_water_demand(&self) -> JsonAccessResult<&Map<std::string::String, JsonValue>> {
        self.root_object("HotWaterDemand")
    }

    fn hot_water_demand_mut(
        &mut self,
    ) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.root_object_mut("HotWaterDemand")
    }

    pub fn set_water_distribution(
        &mut self,
        distribution_value: JsonValue,
    ) -> JsonAccessResult<&Self> {
        self.hot_water_demand_mut()?
            .insert("Distribution".into(), distribution_value);

        Ok(self)
    }

    pub fn set_shower(&mut self, shower_value: JsonValue) -> anyhow::Result<&Self> {
        self.hot_water_demand_mut()?
            .insert("Shower".into(), shower_value);

        Ok(self)
    }

    pub fn set_bath(&mut self, bath_value: JsonValue) -> anyhow::Result<&Self> {
        self.hot_water_demand_mut()?
            .insert("Bath".into(), bath_value);

        Ok(self)
    }

    pub fn set_other_water_use(
        &mut self,
        other_water_use_value: JsonValue,
    ) -> anyhow::Result<&Self> {
        self.hot_water_demand_mut()?
            .insert("Other".into(), other_water_use_value);

        Ok(self)
    }

    pub fn remove_wwhrs(&mut self) -> JsonAccessResult<&mut Self> {
        self.remove_root_key("WWHRS")
    }

    pub(crate) fn wwhrs(&self) -> anyhow::Result<Option<WasteWaterHeatRecovery>> {
        Ok(self
            .root()?
            .get("WWHRS")
            .map(|wwhrs| match serde_json::from_value(wwhrs.to_owned()) {
                Ok(wwhrs) => Ok(wwhrs),
                Err(err) => Err(err),
            })
            .transpose()?)
    }

    pub(crate) fn set_wwhrs(&mut self, wwhrs: JsonValue) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key("WWHRS", wwhrs)
    }

    pub fn remove_space_heat_systems(&mut self) -> JsonAccessResult<&mut Self> {
        self.remove_root_key("SpaceHeatSystem")
    }

    #[cfg(test)]
    pub(crate) fn space_heat_system_for_key(
        &self,
        key: &str,
    ) -> JsonAccessResult<Option<&JsonValue>> {
        Ok(self.root_object("SpaceHeatSystem")?.get(key))
    }

    pub fn set_space_heat_system_for_key(
        &mut self,
        key: &str,
        space_heat_system_value: JsonValue,
    ) -> JsonAccessResult<&Self> {
        self.root_object_entry_mut("SpaceHeatSystem")?
            .insert(key.into(), space_heat_system_value);
        Ok(self)
    }

    pub fn remove_space_cool_systems(&mut self) -> JsonAccessResult<&mut Self> {
        self.remove_root_key("SpaceCoolSystem")
    }

    pub fn set_space_cool_system_for_key(
        &mut self,
        key: &str,
        space_cool_system_value: JsonValue,
    ) -> JsonAccessResult<&Self> {
        self.root_object_entry_mut("SpaceCoolSystem")?
            .insert(key.into(), space_cool_system_value);
        Ok(self)
    }

    pub(crate) fn set_on_site_generation(
        &mut self,
        on_site_generation: JsonValue,
    ) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key("OnSiteGeneration", on_site_generation)
    }

    pub fn remove_on_site_generation(&mut self) -> JsonAccessResult<&mut Self> {
        self.remove_root_key("OnSiteGeneration")
    }

    #[cfg(test)]
    pub(crate) fn on_site_generation(
        &self,
    ) -> JsonAccessResult<Option<&Map<std::string::String, JsonValue>>> {
        self.optional_root_object("OnSiteGeneration")
    }

    pub fn remove_all_diverters_from_energy_supplies(&mut self) -> JsonAccessResult<&mut Self> {
        self.root_object_entry_mut("EnergySupply")?
            .values_mut()
            .filter_map(|value| value.as_object_mut())
            .for_each(|energy_supply| {
                energy_supply.shift_remove("diverter");
            });
        Ok(self)
    }

    pub(crate) fn add_energy_supply_for_key(
        &mut self,
        energy_supply_key: &str,
        energy_supply_details: JsonValue,
    ) -> JsonAccessResult<()> {
        self.root_object_entry_mut("EnergySupply")?
            .insert(energy_supply_key.into(), energy_supply_details);

        Ok(())
    }

    #[cfg(test)]
    pub(crate) fn energy_supply_by_key(
        &self,
        energy_supply_key: &str,
    ) -> JsonAccessResult<Option<&Map<std::string::String, JsonValue>>> {
        Ok(self
            .root_object("EnergySupply")?
            .get(energy_supply_key)
            .and_then(|energy_supply| energy_supply.as_object()))
    }

    #[cfg(test)]
    pub(crate) fn add_diverter_to_energy_supply(
        &mut self,
        energy_supply_key: &str,
        diverter: JsonValue,
    ) -> JsonAccessResult<()> {
        if let Some(energy_supply) = self
            .root_object_mut("EnergySupply")?
            .get_mut(energy_supply_key)
            .and_then(|energy_supply| energy_supply.as_object_mut())
        {
            energy_supply.insert("diverter".into(), diverter);
        }
        Ok(())
    }

    #[cfg(test)]
    pub(crate) fn add_electric_battery_to_energy_supply(
        &mut self,
        energy_supply_key: &str,
        electric_battery: JsonValue,
    ) -> JsonAccessResult<()> {
        if let Some(energy_supply) = self
            .root_object_mut("EnergySupply")?
            .get_mut(energy_supply_key)
            .and_then(|energy_supply| energy_supply.as_object_mut())
        {
            energy_supply.insert("ElectricBattery".into(), electric_battery);
        }

        Ok(())
    }

    pub fn remove_all_batteries_from_energy_supplies(&mut self) -> JsonAccessResult<&mut Self> {
        for energy_supply in self
            .root_object_mut("EnergySupply")?
            .values_mut()
            .flat_map(|v| v.as_object_mut())
        {
            energy_supply.shift_remove("ElectricBattery");
        }
        Ok(self)
    }

    pub fn external_conditions(&self) -> anyhow::Result<ExternalConditionsInput> {
        serde_json::from_value(
            self.root()?
                .get("ExternalConditions")
                .ok_or(json_error("ExternalConditions not found"))?
                .to_owned(),
        )
        .map_err(Into::into)
    }

    fn all_building_elements_mut_of_types(
        &mut self,
        types: &[&str],
    ) -> JsonAccessResult<Vec<&mut Map<std::string::String, JsonValue>>> {
        Ok(self
            .zone_node_mut()?
            .values_mut()
            .filter_map(|zone| {
                zone.get_mut("BuildingElement")
                    .and_then(|building_element_node| building_element_node.as_object_mut())
            })
            .flat_map(|building_elements| building_elements.values_mut())
            .filter(|building_element| {
                building_element
                    .get("type")
                    .and_then(|building_element_type| building_element_type.as_str())
                    .is_some_and(|building_element_type| types.contains(&building_element_type))
            })
            .filter_map(|element| element.as_object_mut())
            .collect())
    }

    pub fn all_transparent_building_elements_mut(
        &mut self,
    ) -> JsonAccessResult<Vec<&mut Map<std::string::String, JsonValue>>> {
        self.all_building_elements_mut_of_types(&["BuildingElementTransparent"])
    }

    pub(crate) fn all_ground_building_elements_mut(
        &mut self,
    ) -> JsonAccessResult<Vec<&mut Map<std::string::String, JsonValue>>> {
        self.all_building_elements_mut_of_types(&["BuildingElementGround"])
    }

    pub(crate) fn all_opaque_and_adjztu_building_elements_mut_u_values(
        &mut self,
    ) -> JsonAccessResult<Vec<&mut Map<std::string::String, JsonValue>>> {
        self.all_building_elements_mut_of_types(&[
            "BuildingElementOpaque",
            "BuildingElementAdjacentUnconditionedSpace_Simple",
        ])
    }

    pub(crate) fn max_base_height_from_building_elements(&self) -> JsonAccessResult<Option<f64>> {
        Ok(self
            .zone_node()?
            .values()
            .filter_map(|zone| zone.get("BuildingElement"))
            .filter_map(|building_element| building_element.as_object())
            .flat_map(|building_element_node| building_element_node.values())
            .filter_map(|building_element| {
                building_element.get("base_height").and_then(|h| h.as_f64())
            })
            .max_by(|a, b| a.total_cmp(b)))
    }

    pub(crate) fn set_numeric_field_for_building_element(
        &mut self,
        building_element_reference: &str,
        field: &str,
        value: f64,
    ) -> anyhow::Result<()> {
        *self.zone_node_mut()?
            .values_mut()
            .filter_map(|zone| zone.get_mut("BuildingElement").and_then(|el| el.as_object_mut()))
            .flatten()
            .find(|(name, _value)| *name == building_element_reference)
            .ok_or(anyhow!("Could not find building element with reference '{building_element_reference}'"))?
            .1
            .get_mut(field)
            .ok_or(anyhow!("Could not find field '{field}' on building element with reference '{building_element_reference}'"))? = json!(value);

        Ok(())
    }

    pub fn all_building_elements(
        &self,
    ) -> anyhow::Result<IndexMap<smartstring::alias::String, BuildingElement>> {
        self.zone_node()?
            .values()
            .filter_map(|zone| zone.get("BuildingElement").and_then(|el| el.as_object()))
            .flatten()
            .map(|(name, el)| {
                Ok((
                    smartstring::alias::String::from(name),
                    serde_json::from_value(el.to_owned())?,
                ))
            })
            .collect()
    }

    #[cfg(test)]
    pub(crate) fn building_element_by_key(
        &self,
        zone_key: &str,
        key: &str,
    ) -> JsonAccessResult<&Map<std::string::String, JsonValue>> {
        self.specific_zone(zone_key)?
            .get("BuildingElement")
            .ok_or(json_error("BuildingElement node not present"))?
            .as_object()
            .ok_or(json_error("BuildingElement node was not an object"))?
            .get(key)
            .ok_or(json_error(format!(
                "BuildingElement with name {key} was not present"
            )))?
            .as_object()
            .ok_or(json_error(
                "Building element with name {key} not provided as an object",
            ))
    }

    pub fn all_energy_supply_fuel_types(
        &self,
    ) -> JsonAccessResult<HashSet<smartstring::alias::String>> {
        let mut fuel_types = HashSet::new();
        for fuel in self
            .root_object("EnergySupply")?
            .values()
            .flat_map(|supply| supply.get("fuel").map(|fuel| fuel.as_str()))
            .flatten()
        {
            fuel_types.insert(smartstring::alias::String::from(fuel));
        }

        Ok(fuel_types)
    }

    pub fn has_appliances(&self) -> JsonAccessResult<bool> {
        Ok(self.root()?.contains_key("Appliances"))
    }

    pub fn merge_in_appliances(
        &mut self,
        appliances: &IndexMap<&str, JsonValue>,
    ) -> anyhow::Result<()> {
        let mut appliances_value = serde_json::to_value(appliances.to_owned())?;
        let appliances = appliances_value
            .as_object_mut()
            .ok_or(anyhow!("Appliances were not an object when expected to be"))?;
        let existing_appliances = self.root_object_entry_mut("Appliances")?;
        existing_appliances.append(appliances);

        Ok(())
    }

    pub fn remove_appliance(&mut self, appliance_key: &str) -> JsonAccessResult<&Self> {
        // we use .shift_remove instead of remove here
        // to preserve the relative order of the appliances
        self.root_object_entry_mut("Appliances")?
            .shift_remove(appliance_key);

        Ok(self)
    }

    pub fn appliances_contain_key(&self, name: &str) -> bool {
        self.root_object("Appliances")
            .ok()
            .is_some_and(|appliances| appliances.contains_key(name))
    }

    pub fn appliance_key_has_reference(
        &self,
        key: &str,
        reference: &str,
    ) -> JsonAccessResult<bool> {
        let empty_map = Map::new();
        Ok(self
            .root_object("Appliances")
            .unwrap_or(&empty_map)
            .get(key)
            .and_then(|value| value.as_str())
            .is_some_and(|appliance_reference| appliance_reference == reference))
    }

    pub fn appliance_keys(&self) -> JsonAccessResult<Vec<smartstring::alias::String>> {
        let empty_map = Map::new();
        Ok(self
            .root_object("Appliances")
            .unwrap_or(&empty_map)
            .keys()
            .map(smartstring::alias::String::from)
            .collect())
    }

    pub fn appliance_with_key(&self, key: &str) -> JsonAccessResult<Option<&JsonValue>> {
        Ok(match self.root_object("Appliances") {
            Err(_) => return Ok(None),
            Ok(appliances) => appliances.get(key),
        })
    }

    pub(crate) fn appliance_with_key_mut(
        &mut self,
        key: &str,
    ) -> JsonAccessResult<Option<&mut JsonValue>> {
        Ok(self.root_object_entry_mut("Appliances")?.get_mut(key))
    }

    pub fn clone_appliances(&self) -> Map<std::string::String, JsonValue> {
        self.root_object("Appliances")
            .cloned()
            .unwrap_or(Map::new())
    }

    pub fn tariff_schedule(&self) -> anyhow::Result<Option<NumericSchedule>> {
        self.root_object("Tariff")
            .ok()
            .cloned()
            .and_then(|tariff| tariff.get("schedule").cloned())
            .map(|schedule| serde_json::from_value(schedule.clone()))
            .transpose()
            .map_err(|err| anyhow!(err))
    }

    pub fn energy_supply_for_appliance(&self, key: &str) -> anyhow::Result<&str> {
        let appliances = self.root_object("Appliances")?;

        appliances
            .get(key)
            .and_then(|appliance| appliance.as_object())
            .ok_or_else(|| anyhow!("No {key} object in appliances input"))?
            .get("Energysupply")
            .and_then(|supply| supply.as_str())
            .ok_or_else(|| anyhow!("No energy supply for appliance '{key}'"))
    }

    pub fn loadshifting_for_appliance(
        &self,
        appliance_key: &str,
    ) -> JsonAccessResult<Option<Map<std::string::String, JsonValue>>> {
        let appliance = self.appliance_with_key(appliance_key)?;

        Ok(appliance
            .and_then(|appliance| appliance.get("loadshifting"))
            .and_then(|load_shifting| load_shifting.as_object())
            .cloned())
    }

    pub fn set_loadshifting_for_appliance(
        &mut self,
        appliance_key: &str,
        new_load_shifting: JsonValue,
    ) -> JsonAccessResult<()> {
        let mut appliance = self.appliance_with_key_mut(appliance_key)?;
        if let Some(appliance) = appliance
            .as_mut()
            .and_then(|appliance| appliance.as_object_mut())
        {
            appliance.insert("loadshifting".into(), new_load_shifting);
        }

        Ok(())
    }

    fn infiltration_ventilation_node_mut(
        &mut self,
    ) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.root_object_mut("InfiltrationVentilation")
    }

    pub fn mechanical_ventilations_for_processing(
        &mut self,
    ) -> JsonAccessResult<Vec<&mut Map<std::string::String, JsonValue>>> {
        let mech_vents = match self
            .infiltration_ventilation_node_mut()?
            .get_mut("MechanicalVentilation")
            .and_then(|v| v.as_object_mut())
        {
            None => return Ok(Vec::new()),
            Some(mech_vents) => mech_vents,
        };
        Ok(mech_vents
            .values_mut()
            .filter_map(|v| v.as_object_mut())
            .collect())
    }

    pub fn keyed_mechanical_ventilations_for_processing(
        &mut self,
    ) -> JsonAccessResult<
        IndexMap<smartstring::alias::String, &mut Map<std::string::String, JsonValue>>,
    > {
        let mech_vents = match self
            .infiltration_ventilation_node_mut()?
            .get_mut("MechanicalVentilation")
            .and_then(|v| v.as_object_mut())
        {
            None => return Ok(Default::default()),
            Some(mech_vents) => mech_vents,
        };
        Ok(mech_vents
            .iter_mut()
            .filter_map(|(name, v)| {
                let mech_vent = match v.as_object_mut() {
                    None => return None,
                    Some(mech_vent) => mech_vent,
                };
                Some((smartstring::alias::String::from(name), mech_vent))
            })
            .collect())
    }

    pub fn has_mechanical_ventilation(&self) -> bool {
        self.root_object("InfiltrationVentilation")
            .ok()
            .is_some_and(|node| node.contains_key("MechanicalVentilation"))
    }

    pub fn reset_mechanical_ventilation(&mut self) -> JsonAccessResult<&Self> {
        self.root_object_entry_mut("InfiltrationVentilation")?
            .shift_remove("MechanicalVentilation");
        Ok(self)
    }

    pub fn add_mechanical_ventilation(
        &mut self,
        vent_name: &str,
        mech_vent: JsonValue,
    ) -> anyhow::Result<()> {
        let infiltration_ventilation_node = self
            .input
            .get_mut("InfiltrationVentilation")
            .ok_or(json_error("InfiltrationVentilation node not found"))?
            .as_object_mut()
            .ok_or(json_error("InfiltrationVentilation node is not an object"))?;
        let mech_vent_map = infiltration_ventilation_node
            .entry("MechanicalVentilation")
            .or_insert(json!({}))
            .as_object_mut()
            .ok_or(json_error("MechanicalVentilation node is not an object"))?;
        mech_vent_map.insert(vent_name.into(), mech_vent);

        Ok(())
    }

    pub fn appliance_gains_events(
        &self,
    ) -> anyhow::Result<IndexMap<smartstring::alias::String, Vec<ApplianceGainsEvent>>> {
        let appliance_gains = match self.root_object("ApplianceGains") {
            Ok(appliance_gains) => appliance_gains,
            Err(_) => return Ok(IndexMap::new()),
        };
        appliance_gains
            .iter()
            .map(
                |(name, gain)| -> Result<(smartstring::alias::String, Vec<ApplianceGainsEvent>), _> {
                    Ok((
                        smartstring::alias::String::from(name),
                        serde_json::from_value(
                            gain.get("Events")
                                .and_then(|events| events.is_array().then_some(events))
                                .cloned()
                                .unwrap_or(json!([])),
                        )?,
                    ))
                },
            )
            .collect::<anyhow::Result<_>>()
    }

    pub fn set_window_adjust_control_for_infiltration_ventilation(
        &mut self,
        control: &str,
    ) -> JsonAccessResult<&Self> {
        self.infiltration_ventilation_node_mut()?
            .insert("Control_WindowAdjust".into(), control.into());
        Ok(self)
    }

    pub fn set_vent_adjust_min_control_for_infiltration_ventilation(
        &mut self,
        control: &str,
    ) -> JsonAccessResult<&Self> {
        self.infiltration_ventilation_node_mut()?
            .insert("Control_VentAdjustMin".into(), control.into());
        Ok(self)
    }

    pub fn set_vent_adjust_max_control_for_infiltration_ventilation(
        &mut self,
        control: &str,
    ) -> JsonAccessResult<&Self> {
        self.infiltration_ventilation_node_mut()?
            .insert("Control_VentAdjustMax".into(), control.into());
        Ok(self)
    }

    pub fn infiltration_ventilation_is_noise_nuisance(&self) -> bool {
        self.root_object("InfiltrationVentilation")
            .ok()
            .and_then(|infiltration| infiltration.get("noise_nuisance"))
            .and_then(|nuisance| nuisance.as_bool())
            .unwrap_or(false)
    }

    pub(crate) fn infiltration_ventilation_mut(
        &mut self,
    ) -> JsonAccessResult<&mut Map<std::string::String, JsonValue>> {
        self.root_object_entry_mut("InfiltrationVentilation")
    }

    pub(crate) fn set_heat_source_wet(
        &mut self,
        heat_source_wet: JsonValue,
    ) -> JsonAccessResult<()> {
        self.root_mut()?
            .insert("HeatSourceWet".into(), heat_source_wet);
        Ok(())
    }

    pub(crate) fn heat_source_wet(
        &self,
    ) -> anyhow::Result<IndexMap<smartstring::alias::String, HeatSourceWetDetails>> {
        self.root()?
            .get("HeatSourceWet")
            .and_then(|value| value.as_object())
            .into_iter()
            .flatten()
            .map(|(name, source)| {
                Ok((
                    smartstring::alias::String::from(name),
                    serde_json::from_value(source.clone())?,
                ))
            })
            .collect::<anyhow::Result<_, _>>()
    }

    pub(crate) fn cold_water_source(&self) -> anyhow::Result<ColdWaterSourceInput> {
        Ok(serde_json::from_value(
            self.root()?
                .get("ColdWaterSource")
                .cloned()
                .ok_or(json_error("ColdWaterSource was not present"))?,
        )?)
    }

    #[cfg(test)]
    pub(crate) fn set_storeys_in_building(&mut self, storeys: usize) -> JsonAccessResult<&Self> {
        self.root_object_mut("General")?
            .insert("storeys_in_building".into(), json!(storeys));

        Ok(self)
    }

    pub(crate) fn storeys_in_building(&self) -> JsonAccessResult<usize> {
        Ok(self
            .input
            .get("General")
            .ok_or(json_error("General node not found"))?
            .get("storeys_in_building")
            .ok_or(json_error("storeys_in_building field not found"))?
            .as_u64()
            .ok_or(json_error(
                "storeys_in_building field is not a positive integer",
            ))? as usize)
    }

    pub(crate) fn build_type(&self) -> JsonAccessResult<smartstring::alias::String> {
        Ok(self
            .root_object("General")?
            .get("build_type")
            .ok_or(json_error(
                "There was no build_type field on the General input object",
            ))?
            .as_str()
            .ok_or(json_error("The build_type field was not a string"))?
            .into())
    }

    pub(crate) fn hot_water_cylinder_volume(&self) -> JsonAccessResult<Option<f64>> {
        Ok(self
            .root_object("HotWaterSource")?
            .get("hw cylinder")
            .and_then(|cylinder| cylinder.get("volume"))
            .and_then(|v| v.as_f64()))
    }

    pub(crate) fn ground_floor_area(&self) -> JsonAccessResult<Option<f64>> {
        Ok(self
            .root()?
            .get("GroundFloorArea")
            .and_then(|area| area.as_f64()))
    }

    pub(crate) fn primary_pipework_clone(&self) -> anyhow::Result<Option<Vec<WaterPipework>>> {
        Ok(self
            .hot_water_source()?
            .get("hw cylinder")
            .and_then(|cylinder| cylinder.get("primary_pipework"))
            .and_then(|primary_pipework| primary_pipework.is_array().then_some(primary_pipework))
            .map(|primary_pipework| serde_json::from_value(primary_pipework.to_owned()))
            .transpose()?)
    }

    pub(crate) fn water_heating_event_by_type_and_name(
        &self,
        event_type: &str,
        event_name: &str,
    ) -> anyhow::Result<Option<Vec<WaterHeatingEvent>>> {
        let result = Ok(self
            .root_object("Events")?
            .get(event_type)
            .and_then(|event_group| event_group.get(event_name))
            .map(|events| serde_json::from_value(events.to_owned()))
            .transpose()?);

        if let Ok(None) = result {
            println!(
                "No events found for event type {} and name {}",
                event_type, event_name
            );
            println!("Events: {:?}", self.root_object("Events")?.keys());
        }

        result
    }

    pub(crate) fn part_o_active_cooling_required(&self) -> JsonAccessResult<Option<bool>> {
        Ok(match self.input.get("PartO_active_cooling_required") {
            None => None,
            Some(JsonValue::Bool(whether)) => Some(*whether),
            Some(_) => {
                return Err(json_error(
                    "PartO_active_cooling_required field not a boolean",
                ));
            }
        })
    }

    #[cfg(test)]
    pub fn set_part_o_active_cooling_required(
        &mut self,
        required: bool,
    ) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key("PartO_active_cooling_required", json!(required))
    }

    #[cfg(test)]
    pub(crate) fn set_zone(&mut self, zone: JsonValue) -> JsonAccessResult<&mut Self> {
        self.set_on_root_key("Zone", zone)
    }
}

#[cfg(test)]
mod tests {
    use hem::input::ingest_for_processing;
    use itertools::Itertools;
    use rstest::{fixture, rstest};
    use std::fs::File;
    use walkdir::{DirEntry, WalkDir};

    fn files_with_root(root: &str) -> Vec<DirEntry> {
        WalkDir::new(root)
            .into_iter()
            .filter_map(Result::ok)
            .filter(|e| {
                !e.file_type().is_dir()
                    && e.file_name().to_str().unwrap().ends_with("json")
                    && !e
                        .path()
                        .parent()
                        .unwrap()
                        .to_str()
                        .unwrap()
                        .ends_with("results") // don't test against files in results output directories
            })
            .collect_vec()
    }

    #[fixture]
    fn fhs_files() -> Vec<DirEntry> {
        files_with_root("./examples/input/wrappers/future_homes_standard")
    }

    #[rstest]
    fn should_successfully_parse_all_fhs_demo_files(fhs_files: Vec<DirEntry>) {
        for entry in fhs_files {
            let parsed = ingest_for_processing(File::open(entry.path()).unwrap());
            assert!(
                parsed.is_ok(),
                "error was {:?} when parsing file {}",
                parsed.err().unwrap(),
                entry.file_name().to_str().unwrap()
            );
        }
    }
}

#[cfg(test)]
mod accessors_tests {
    use super::*;
    use rstest::*;

    #[fixture]
    fn events_input() -> InputForProcessing {
        let events_input_json = json!({
            "Events": {
                "Shower": {
                  "IES": [
                    {
                      "start": 4.1,
                      "duration": 6,
                      "temperature": 41.0
                    },
                    {
                      "start": 4.5,
                      "duration": 6,
                      "temperature": 41.0
                    },
                    {
                      "start": 6,
                      "duration": 6,
                      "temperature": 41.0
                    }
                  ],
                  "mixer": [
                    {
                      "start": 7,
                      "duration": 6,
                      "temperature": 41.0
                    }
                  ]
                }
          }
        });

        InputForProcessing {
            input: events_input_json,
        }
    }

    #[rstest]
    fn test_water_heating_event_by_type_and_name_when_exists(events_input: InputForProcessing) {
        assert_eq!(
            events_input
                .water_heating_event_by_type_and_name("Shower", "mixer")
                .unwrap(),
            Some(vec![WaterHeatingEvent {
                start: 7.,
                duration: Some(6.),
                volume: None,
                temperature: 41.0
            }])
        );
    }

    #[fixture]
    fn hot_water_cylinder_input() -> InputForProcessing {
        let hot_water_source_json = json!({
            "HotWaterSource": {
                "hw cylinder": {
                  "type": "StorageTank",
                  "volume": 80.0,
                  "daily_losses": 1.68,
                  "min_temp": 52.0,
                  "setpoint_temp": 55.0,
                  "ColdWaterSource": "mains water",
                  "HeatSource": {
                    "hp": {
                      "type": "HeatSourceWet",
                      "name": "hp",
                      "temp_flow_limit_upper": 65,
                      "ColdWaterSource": "mains water",
                      "EnergySupply": "mains elec",
                      "Control": "hw timer",
                      "heater_position": 0.1,
                      "thermostat_position": 0.33
                    }
                  },
                  "primary_pipework": [
                    {
                      "location": "external",
                      "internal_diameter_mm": 26.,
                      "external_diameter_mm": 28.,
                      "length": 3.0,
                      "insulation_thermal_conductivity": 0.037,
                      "insulation_thickness_mm": 25.,
                      "surface_reflectivity": false,
                      "pipe_contents": "water"
                    }
                  ]
                }
              }
        });

        InputForProcessing {
            input: hot_water_source_json,
        }
    }

    #[rstest]
    fn test_hot_water_cylinder_volume(hot_water_cylinder_input: InputForProcessing) {
        assert_eq!(
            hot_water_cylinder_input
                .hot_water_cylinder_volume()
                .unwrap(),
            Some(80.0)
        )
    }

    #[rstest]
    fn test_set_gains_for_field() {
        let mut input = InputForProcessing {
            input: json!({
                "ApplianceGains": {
                    "Clothes_washing": 2,
                }
            }),
        };

        let expected_appliance_gains = json!({
            "Clothes_washing": 2,
            "Clothes_drying": 42,
        });

        input
            .set_gains_for_field("Clothes_drying", json!(42))
            .unwrap();

        assert_eq!(
            json!(input.root_object("ApplianceGains").unwrap()),
            expected_appliance_gains
        );
    }

    #[rstest]
    fn test_water_heating_events_of_types(events_input: InputForProcessing) {
        let actual = events_input
            .water_heating_events_of_types(&["Shower"])
            .unwrap();
        let expected = vec![
            json!({
              "start": 4.1,
              "duration": 6,
              "temperature": 41.0
            }),
            json!({
              "start": 4.5,
              "duration": 6,
              "temperature": 41.0
            }),
            json!({
              "start": 6,
              "duration": 6,
              "temperature": 41.0
            }),
            json!({
              "start": 7,
              "duration": 6,
              "temperature": 41.0
            }),
        ];
        assert_eq!(actual, expected);
    }

    #[rstest]
    fn test_reset_internal_gains() {
        let base_input = json!({
            "InternalGains": {
                "metabolic gains": {
                    "start_day": 0,
                    "time_series_step": 1,
                    "schedule": {
                        "main": [1305.6, 1876.8, 2978.4, 2121.6, 3631.2, 2284.8, 4161.6, 3304.8]
                    }
                }
            }
        });
        let mut input = InputForProcessing { input: base_input };
        input.reset_internal_gains().unwrap();
        assert_eq!(input.input, json!({"InternalGains": {}}));
    }
}
