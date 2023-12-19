use crate::core::controls::time_control::SetpointTimeControl;
use crate::core::energy_supply::energy_supply::EnergySupply;
use crate::core::space_heat_demand::building_element::{
    area_for_building_element_input, element_from_named, mid_height_for, orientation_for,
    projected_height_for_transparent_element,
};
use crate::core::space_heat_demand::zone::NamedBuildingElement;
use crate::core::units::{
    celsius_to_kelvin, LITRES_PER_CUBIC_METRE, SECONDS_PER_HOUR, WATTS_PER_KILOWATT,
};
use crate::external_conditions::ExternalConditions;
use crate::input::{
    BuildingElement, InfiltrationBuildType, InfiltrationShelterType, InfiltrationTestType,
};
use crate::simulation_time::SimulationTimeIterator;
use std::collections::HashSet;

fn air_change_rate_to_flow_rate(air_change_rate: f64, zone_volume: f64) -> f64 {
    air_change_rate * zone_volume / SECONDS_PER_HOUR as f64
}

pub trait VentilationElement {
    fn h_ve_heat_transfer_coefficient(
        &self,
        zone_volume: f64,
        throughput_factor: Option<f64>,
        timestep_idx: Option<usize>,
    ) -> f64;

    fn temp_supply(&self, timestep_index: usize) -> f64;

    /// Calculate the heat transfer coefficient (h_ve), in W/K,
    ///         according to ISO 52016-1:2017, Section 6.5.10.1, for a constant average windspeed
    ///
    /// # Arguments
    ///
    ///  * `zone_volume` - volume of zone, in m3
    fn h_ve_average_heat_transfer_coefficient(&self, zone_volume: f64) -> f64;
}

pub struct VentilationElementInfiltration {
    external_conditions: ExternalConditions,
    volume: f64,
    infiltration_rate_from_openings: f64,
    q50_divisor: f64,
    shelter_factor: f64,
    infiltration_rate: f64,
}

impl VentilationElementInfiltration {
    /// # Arguments
    ///
    /// * `storeys_in_building` - total number of storeys in building
    /// * `shelter` - exposure level of the building i.e. very sheltered, sheltered, normal, or exposed
    /// * `build_type` - type of building e.g. house, flat, etc.
    /// * `pressure_test_result_ach` - result of pressure test, in ach
    /// * `pressure_test_type` - measurement used for pressure test i.e. based on air change rate value at 50 Pa (50Pa) or 4 Pa (4Pa)
    /// * `envelope_area` - total envelope area of the building including party walls and floors, in m^2
    /// * `volume` - total volume of dwelling, m^3
    /// * `sheltered_sides` - number of sides of the building which are sheltered (max 4)
    /// * `open_chimneys` - number of chimneys / flues attached to closed fire
    /// * `closed_fire` - number of chimneys / flues attached to closed fire
    /// * `flues_d` - number of flues attached to other heater
    /// * `flues_e` - number of flues attached to other heater
    /// * `blocked_chimneys` - number of blocked chimneys
    /// * `extract_fans` - number of intermittent extract fans
    /// * `passive_vents` - number of passive vents
    /// * `flueless_gas_fires` - number of flueless gas fires
    /// * `external_conditions`
    /// * `storey_of_dwelling` - the storey in the building of the dwelling, if this is a flat
    pub fn new(
        storeys_in_building: u32,
        shelter: InfiltrationShelterType,
        build_type: InfiltrationBuildType,
        pressure_test_result_ach: f64,
        pressure_test_type: InfiltrationTestType,
        envelope_area: f64,
        volume: f64,
        sheltered_sides: u32,
        open_chimneys: u32,
        open_flues: u32,
        closed_fire: u32,
        flues_d: u32,
        flues_e: u32,
        blocked_chimneys: u32,
        extract_fans: u32,
        passive_vents: u32,
        flueless_gas_fires: u32,
        external_conditions: ExternalConditions,
        storey_of_dwelling: Option<u32>,
    ) -> VentilationElementInfiltration {
        let infiltration_rate_from_openings = ((open_chimneys * INFILTRATION_RATE_CHIMNEY_OPEN)
            + (open_flues * INFILTRATION_RATE_FLUE_OPEN)
            + (closed_fire * INFILTRATION_RATE_FIRE_CLOSED)
            + (flues_d * INFILTRATION_RATE_FLUE_SOLID_FUEL_BOILER)
            + (flues_e * INFILTRATION_RATE_FLUE_OTHER_HEATER)
            + (blocked_chimneys * INFILTRATION_RATE_CHIMNEY_BLOCKED)
            + (extract_fans * INFILTRATION_RATE_EXTRACT_FAN)
            + (passive_vents * INFILTRATION_RATE_PASSIVE_STACK_VENT)
            + (flueless_gas_fires * INFILTRATION_RATE_FIRE_GAS))
            as f64
            / volume;

        let storey = match build_type {
            InfiltrationBuildType::House => storeys_in_building,
            InfiltrationBuildType::Flat => storey_of_dwelling.expect(
                "The storey of the dwelling was expected to have been provided as this is a flat.",
            ),
        };

        let q50_divisor = init_divisor(build_type, storey, shelter);

        let shelter_factor = init_shelter_factor(sheltered_sides);

        let infiltration_rate = init_infiltration(
            pressure_test_type,
            pressure_test_result_ach,
            envelope_area,
            volume,
            q50_divisor,
            infiltration_rate_from_openings,
            shelter_factor,
        );

        VentilationElementInfiltration {
            external_conditions,
            volume,
            infiltration_rate_from_openings,
            q50_divisor,
            shelter_factor,
            infiltration_rate,
        }
    }

    pub fn infiltration_rate_from_openings(&self) -> f64 {
        self.infiltration_rate_from_openings
    }

    pub fn q50_divisor(&self) -> f64 {
        self.q50_divisor
    }

    pub fn shelter_factor(&self) -> f64 {
        self.shelter_factor
    }

    pub fn infiltration_rate(&self) -> f64 {
        self.infiltration_rate
    }
}

impl VentilationElement for VentilationElementInfiltration {
    /// Calculate the heat transfer coefficient (h_ve), in W/K,
    ///         according to ISO 52016-1:2017, Section 6.5.10.1
    ///
    /// # Arguments
    ///
    /// * `zone_volume` - volume of zone, in m3
    /// * `timestep_idx` - the index of the current timestamp
    fn h_ve_heat_transfer_coefficient(
        &self,
        zone_volume: f64,
        _throughput_factor: Option<f64>,
        timestep_idx: Option<usize>,
    ) -> f64 {
        P_A * C_A
            * (self.infiltration_rate()
                * self
                    .external_conditions
                    .wind_speed_for_timestep_idx(timestep_idx.unwrap())
                / 4.0)
            * (zone_volume / SECONDS_PER_HOUR as f64)
    }

    fn temp_supply(&self, timestep_idx: usize) -> f64 {
        // Calculate the supply temperature of the air flow element
        // according to ISO 52016-1:2017, Section 6.5.10.2
        self.external_conditions
            .air_temp_for_timestep_idx(timestep_idx)
    }

    fn h_ve_average_heat_transfer_coefficient(&self, zone_volume: f64) -> f64 {
        P_A * C_A
            * (self.infiltration_rate() * self.external_conditions.wind_speed_annual().unwrap()
                / 4.0
                * zone_volume
                / SECONDS_PER_HOUR as f64)
    }
}

// # Infiltration rates for openings (m3 per hour)
const INFILTRATION_RATE_CHIMNEY_OPEN: u32 = 80;
const INFILTRATION_RATE_CHIMNEY_BLOCKED: u32 = 20;
const INFILTRATION_RATE_FLUE_OPEN: u32 = 20;
const INFILTRATION_RATE_FLUE_SOLID_FUEL_BOILER: u32 = 20;
const INFILTRATION_RATE_FLUE_OTHER_HEATER: u32 = 35;
const INFILTRATION_RATE_FIRE_CLOSED: u32 = 10;
const INFILTRATION_RATE_FIRE_GAS: u32 = 40;
const INFILTRATION_RATE_EXTRACT_FAN: u32 = 10;
const INFILTRATION_RATE_PASSIVE_STACK_VENT: u32 = 10;

fn init_divisor(
    build_type: InfiltrationBuildType,
    storey: u32,
    shelter: InfiltrationShelterType,
) -> f64 {
    divisor_for(dwelling_type(build_type, storey), shelter)
}

/// Divisors to convert air change rate at 50 Pa to infiltration
/// Values for "Normal" House 1-2 storey and Flat storeys 1-10 are from CIBSE Guide A
/// Values for "Exposed" based on CIBSE Guida A: "on severely exposed sites, a
/// 50% increase to the tabulated values should be allowed.
/// Values for "Sheltered" based on CIBSE Guide A: "on sheltered sites, the
/// infiltration rate may be reduced by 33%".
/// Values for "Very sheltered" assume a reduction of 50%.
/// Values for Flat storeys 11+ are extrapolated based on profiles of wind
/// speed vs. height, assuming storey height of 3.5 metres.
fn divisor_for(dwelling_type: DwellingType, shelter_type: InfiltrationShelterType) -> f64 {
    match (dwelling_type, shelter_type) {
        (DwellingType::House1Storey, InfiltrationShelterType::VerySheltered) => 41.2,
        (DwellingType::House1Storey, InfiltrationShelterType::Sheltered) => 30.7,
        (DwellingType::House1Storey, InfiltrationShelterType::Normal) => 20.6,
        (DwellingType::House1Storey, InfiltrationShelterType::Exposed) => 13.7,
        (DwellingType::House2Storey, InfiltrationShelterType::VerySheltered) => 34.0,
        (DwellingType::House2Storey, InfiltrationShelterType::Sheltered) => 25.4,
        (DwellingType::House2Storey, InfiltrationShelterType::Normal) => 17.0,
        (DwellingType::House2Storey, InfiltrationShelterType::Exposed) => 11.3,
        (DwellingType::FlatStorey1To5, InfiltrationShelterType::VerySheltered) => 34.6,
        (DwellingType::FlatStorey1To5, InfiltrationShelterType::Sheltered) => 25.8,
        (DwellingType::FlatStorey1To5, InfiltrationShelterType::Normal) => 17.3,
        (DwellingType::FlatStorey1To5, InfiltrationShelterType::Exposed) => 11.5,
        (DwellingType::FlatStorey6To10, InfiltrationShelterType::VerySheltered) => 30.2,
        (DwellingType::FlatStorey6To10, InfiltrationShelterType::Sheltered) => 22.5,
        (DwellingType::FlatStorey6To10, InfiltrationShelterType::Normal) => 15.1,
        (DwellingType::FlatStorey6To10, InfiltrationShelterType::Exposed) => 10.1,
        (DwellingType::FlatStorey11Plus, InfiltrationShelterType::VerySheltered) => 29.3,
        (DwellingType::FlatStorey11Plus, InfiltrationShelterType::Sheltered) => 19.9,
        (DwellingType::FlatStorey11Plus, InfiltrationShelterType::Normal) => 13.7,
        (DwellingType::FlatStorey11Plus, InfiltrationShelterType::Exposed) => 9.3,
    }
}

fn dwelling_type(build_type: InfiltrationBuildType, storey: u32) -> DwellingType {
    match (build_type, storey) {
        (InfiltrationBuildType::House, s) if s == 1 => DwellingType::House1Storey,
        (InfiltrationBuildType::House, _) => DwellingType::House2Storey,
        (InfiltrationBuildType::Flat, s) if s <= 5 => DwellingType::FlatStorey1To5,
        (InfiltrationBuildType::Flat, s) if s <= 10 => DwellingType::FlatStorey6To10,
        (InfiltrationBuildType::Flat, _) => DwellingType::FlatStorey11Plus,
    }
}

enum DwellingType {
    House1Storey,
    House2Storey,
    FlatStorey1To5,
    FlatStorey6To10,
    FlatStorey11Plus,
}

fn init_shelter_factor(sheltered_sides: u32) -> f64 {
    if sheltered_sides > 4 {
        // assume in time we'll be able to apply validation further up to make this impossible
        panic!("There were not expected to be more than 4 sheltered sides");
    }

    1.0 - (0.075 * sheltered_sides as f64)
}

fn init_infiltration(
    pressure_test_type: InfiltrationTestType,
    pressure_test_result_ach: f64,
    envelope_area: f64,
    volume: f64,
    divisor: f64,
    infiltration_rate_from_openings: f64,
    shelter_factor: f64,
) -> f64 {
    ((match pressure_test_type {
        /// If test results are at 4 Pa, convert to equivalent 50 Pa result
        /// before applying divisor.
        /// SAP 10 Technical Paper S10TP-19 "Use of low pressure pulse
        /// test data in SAP" gives the relationship between air
        /// permeability measured at 50 Pa and 4 Pa. The equation below is
        /// based on this but has been converted to work with test results
        /// expressed in ach rather than m3/m2/h.
        InfiltrationTestType::FourPascals => {
            5.254
                * (pressure_test_result_ach.powf(0.9241)
                    * (envelope_area / volume).powf(1.0 - 0.9241))
        }
        InfiltrationTestType::FiftyPascals => pressure_test_result_ach,
    }) / divisor)
        + (infiltration_rate_from_openings * shelter_factor)
}

const P_A: f64 = 1.204; // Air density at 20 degrees C, in kg/m^3 , BS EN ISO 52016-1:2017, Section 6.3.6
const C_A: f64 = 1006.0; // Specific heat of air at constant pressure, in J/(kg K), BS EN ISO 52016-1:2017, Section 6.3.6

pub struct MechanicalVentilationHeatRecovery {
    air_change_rate: f64,
    specific_fan_power: f64,
    efficiency_hr: f64,
    energy_supply: EnergySupply,
    energy_supply_end_user_name: String, // rather than using an EnergySupplyConnection object that encapsulates this
    external_conditions: ExternalConditions,
    simulation_time: SimulationTimeIterator,
}

impl MechanicalVentilationHeatRecovery {
    pub fn new(
        required_air_change_rate: f64,
        specific_fan_power: f64,
        efficiency_hr: f64,
        mut energy_supply: EnergySupply,
        energy_supply_end_user_name: String,
        external_conditions: ExternalConditions,
        simulation_time: SimulationTimeIterator,
    ) -> Self {
        energy_supply.register_end_user_name(energy_supply_end_user_name.clone());
        Self {
            air_change_rate: required_air_change_rate,
            specific_fan_power,
            efficiency_hr,
            energy_supply,
            energy_supply_end_user_name,
            external_conditions,
            simulation_time,
        }
    }

    pub fn fans(
        &mut self,
        zone_volume: f64,
        timestep_index: usize,
        throughput_factor: Option<f64>,
    ) -> f64 {
        let throughput_factor = throughput_factor.unwrap_or(1.0);
        // Calculate energy use by fans (only fans on intake/supply side
        // contribute to internal gains - assume that this is half of the fan power)
        let q_v =
            air_change_rate_to_flow_rate(self.air_change_rate, zone_volume) * throughput_factor;
        let fan_power_w = self.specific_fan_power * (q_v * LITRES_PER_CUBIC_METRE as f64);
        let fan_energy_use_kwh =
            (fan_power_w / WATTS_PER_KILOWATT as f64) * self.simulation_time.step_in_hours();

        let _ = self.energy_supply.demand_energy(
            self.energy_supply_end_user_name.clone(),
            fan_energy_use_kwh,
            timestep_index,
        );

        fan_energy_use_kwh / 2.0
    }

    pub fn efficiency(&self) -> f64 {
        self.efficiency_hr
    }
}

impl VentilationElement for MechanicalVentilationHeatRecovery {
    /// Calculate the heat transfer coefficient (h_ve), in W/K,
    /// according to ISO 52016-1:2017, Section 6.5.10.1
    ///
    /// # Arguments
    /// * `zone_volume` - volume of zone, in m3
    /// * `throughput_factor` - proportional increase in ventilation rate due to over-ventilation requirement
    fn h_ve_heat_transfer_coefficient(
        &self,
        zone_volume: f64,
        throughput_factor: Option<f64>,
        _timestep_idx: Option<usize>,
    ) -> f64 {
        let throughput_factor = throughput_factor.unwrap_or(1.0);

        let q_v =
            air_change_rate_to_flow_rate(self.air_change_rate, zone_volume) * throughput_factor;

        // # Calculate effective flow rate of external air
        // # NOTE: Technically, the MVHR system supplies air at a higher temperature
        // # than the outside air. However, it is simpler to adjust the heat
        // # transfer coefficient h_ve to account for the heat recovery effect
        // # using an "equivalent" or "effective" flow rate of external air.
        let q_v_effective = q_v * (1.0 - self.efficiency_hr);

        // Calculate h_ve according to BS EN ISO 52016-1:2017 section 6.5.10 equation 61
        P_A * C_A * q_v_effective
    }

    fn temp_supply(&self, timestep_index: usize) -> f64 {
        self.external_conditions
            .air_temp_for_timestep_idx(timestep_index)
    }

    fn h_ve_average_heat_transfer_coefficient(&self, zone_volume: f64) -> f64 {
        self.h_ve_heat_transfer_coefficient(zone_volume, None, None)
    }
}

pub struct WholeHouseExtractVentilation {
    air_change_rate: f64,
    specific_fan_power: f64,
    infiltration_rate: f64,
    energy_supply: EnergySupply,
    energy_supply_end_user_name: String,
    external_conditions: ExternalConditions,
    simulation_time: SimulationTimeIterator,
}

impl WholeHouseExtractVentilation {
    pub fn new(
        required_air_change_rate: f64,
        specific_fan_power: f64,
        infiltration_rate: f64,
        energy_supply: EnergySupply,
        energy_supply_end_user_name: String,
        external_conditions: ExternalConditions,
        simulation_time: SimulationTimeIterator,
    ) -> Self {
        Self {
            air_change_rate: required_air_change_rate,
            infiltration_rate,
            specific_fan_power,
            energy_supply,
            energy_supply_end_user_name,
            external_conditions,
            simulation_time,
        }
    }

    /// Calculate air change rate for the system
    //
    /// # Arguments
    /// * `infiltration_rate` - in ach, adjusted for wind speed
    fn air_change_rate(&self, infiltration_rate: f64) -> f64 {
        // # The calculation below is based on SAP 10.2 equations, but with the
        // # sharp "elbow" at infiltration_rate == 0.5 * air_change_rate_req
        // # replaced by a smooth curve between infiltration_rate == 0 and
        // # infiltration_rate = air_change_rate_req.
        // # As we are already accounting for infiltration separately, it is
        // # subtracted from the totals here, compared to the equations in SAP 10.2
        if infiltration_rate < self.air_change_rate {
            self.air_change_rate - infiltration_rate
                + (infiltration_rate.powi(2) * 0.5 / self.air_change_rate)
        } else {
            0.5 * self.air_change_rate
        }
    }

    /// Calculate gains and energy use due to fans
    pub fn fans(
        &mut self,
        zone_volume: f64,
        throughput_factor: Option<f64>,
        timestep_index: usize,
    ) -> f64 {
        let throughput_factor = throughput_factor.unwrap_or(1.0);

        // # Calculate energy use by fans (does not contribute to internal gains as
        // # this is extract-only ventilation)
        let q_v =
            air_change_rate_to_flow_rate(self.air_change_rate, zone_volume) * throughput_factor;
        let fan_power_w = self.specific_fan_power * (q_v * LITRES_PER_CUBIC_METRE as f64);
        let fan_power_use_kwh =
            (fan_power_w / WATTS_PER_KILOWATT as f64) * self.simulation_time.step_in_hours();

        self.energy_supply.demand_energy(
            self.energy_supply_end_user_name.clone(),
            fan_power_use_kwh,
            timestep_index,
        );

        0.0
    }
}

impl VentilationElement for WholeHouseExtractVentilation {
    /// Calculate the heat transfer coefficient (h_ve), in W/K,
    /// according to ISO 52016-1:2017, Section 6.5.10.1
    //
    /// # Arguments
    /// * `zone_volume` - volume of zone, in m3
    /// * `throughput_factor` - proportional increase in ventilation rate due to over-ventilation requirement
    /// * `timestep_index` - index of timestep being calculated
    fn h_ve_heat_transfer_coefficient(
        &self,
        zone_volume: f64,
        throughput_factor: Option<f64>,
        timestep_index: Option<usize>,
    ) -> f64 {
        let throughput_factor = throughput_factor.unwrap_or(1.0);

        let infiltration_rate_adj = self.infiltration_rate
            * self
                .external_conditions
                .wind_speed_for_timestep_idx(timestep_index.unwrap())
            / 4.0;
        let ach = self.air_change_rate(infiltration_rate_adj);
        let q_v = air_change_rate_to_flow_rate(ach, zone_volume) * throughput_factor;

        // Calculate h_ve according to BS EN ISO 52016-1:2017 section 6.5.10 equation 61
        P_A * C_A * q_v
    }

    fn temp_supply(&self, timestep_index: usize) -> f64 {
        self.external_conditions
            .air_temp_for_timestep_idx(timestep_index)
    }

    fn h_ve_average_heat_transfer_coefficient(&self, _zone_volume: f64) -> f64 {
        todo!()
    }
}

//tbc
pub struct NaturalVentilation {
    required_air_change_rate: f64,
    infiltration_rate: f64,
    external_conditions: ExternalConditions,
}

impl VentilationElement for NaturalVentilation {
    fn h_ve_heat_transfer_coefficient(
        &self,
        zone_volume: f64,
        _throughput_factor: Option<f64>,
        timestep_idx: Option<usize>,
    ) -> f64 {
        let infiltration_rate_adj = self.infiltration_rate
            * self
                .external_conditions
                .wind_speed_for_timestep_idx(timestep_idx.unwrap())
            / 4.0;
        let ach = self.air_change_rate(infiltration_rate_adj);
        let q_v = air_change_rate_to_flow_rate(ach, zone_volume);

        // Calculate h_ve according to BS EN ISO 52016-1:2017 section 6.5.10 equation 61
        P_A * C_A * q_v
    }

    fn temp_supply(&self, timestep_index: usize) -> f64 {
        self.external_conditions
            .air_temp_for_timestep_idx(timestep_index)
    }

    fn h_ve_average_heat_transfer_coefficient(&self, _zone_volume: f64) -> f64 {
        todo!()
    }
}

impl NaturalVentilation {
    fn air_change_rate(&self, infiltration_rate: f64) -> f64 {
        // The calculation below is based on SAP 10.2 equations, but with the curve between
        // infiltration_rate == 0 and infiltration_rate == 2 * air_change_rate_req
        // adjusted to handle values of air_change_rate_req other than 0.5.
        // As we are already accounting for infiltration separately, it is
        // subtracted from the totals here, compared to the equations in SAP 10.2
        if infiltration_rate < 2.0 * self.required_air_change_rate {
            self.required_air_change_rate - infiltration_rate
                + (infiltration_rate.powi(2) * 0.25 / self.required_air_change_rate)
        } else {
            0.0
        }
    }
}

const G: f64 = 9.81; // m/s
const C_D: f64 = 0.62; // discharge coefficient
const DC_P: f64 = 0.2 - (-0.25); // Difference in wind pressure coeff from CIBSE Guide A Table 4.12 for urban environment

pub struct WindowOpeningForCooling<'a> {
    window_area_equivalent: f64,
    external_conditions: ExternalConditions,
    openings: Option<Vec<&'a BuildingElement>>, // actually only meaningfully contains BuildingElement::Transparent
    control: &'a SetpointTimeControl,
    natural_ventilation: Option<NaturalVentilation>,
    a_b: Option<f64>,
    a_w: Option<f64>,
    opening_height_diff: f64,
    opening_area_ratio: Option<f64>,
    cross_vent: bool,
    stack_vent: bool,
}

impl<'a> WindowOpeningForCooling<'a> {
    /// Arguments:
    /// `window_area_equivalent` - maximum equivalent area of all openings in the relevant zone
    /// `external_conditions` - reference to ExternalConditions object
    /// `named_openings` - list of openings to be considered
    /// `control` - reference to control object (must implement setpnt function)
    /// `natural_ventilation` - reference to NaturalVentilation object, if building is naturally ventilated
    pub fn new(
        window_area_equivalent: f64,
        external_conditions: ExternalConditions,
        named_openings: &'a Vec<NamedBuildingElement>,
        control: &'a SetpointTimeControl,
        natural_ventilation: Option<NaturalVentilation>,
    ) -> Self {
        let openings = named_openings
            .iter()
            .map(element_from_named)
            .collect::<Vec<&BuildingElement>>();
        // Assign equivalent areas to each window/group in proportion to actual area
        let opening_area_total = openings
            .iter()
            .map(|element| area_for_building_element_input(element))
            .sum::<f64>();
        let opening_area_equiv_total_ratio = window_area_equivalent / opening_area_total;

        // Find orientation of largest window
        // Find height of highest and lowest windows
        let mut largest_op = &openings[0];
        let mut highest_op = &openings[0];
        let mut lowest_op = &openings[0];
        for op in openings[1..].iter() {
            if area_for_building_element_input(op) > area_for_building_element_input(largest_op) {
                largest_op = op;
            }
            if mid_height_for(op) > mid_height_for(highest_op) {
                highest_op = op;
            }
            if mid_height_for(op) < mid_height_for(lowest_op) {
                lowest_op = op;
            }
        }
        let largest_op_orientation = orientation_for(largest_op).unwrap();
        let op_height_threshold =
            (mid_height_for(highest_op).unwrap() - mid_height_for(lowest_op).unwrap()) / 2.0;

        let mut openings_same_side = vec![];
        let mut openings_opp_side = vec![];
        let mut openings_low = vec![];
        let mut openings_high = vec![];
        for op in named_openings.iter() {
            // Determine orientation of other windows relative to largest
            let mut op_rel_orientation =
                (orientation_for(element_from_named(op)).unwrap() - largest_op_orientation).abs();
            if op_rel_orientation > 360.0 {
                op_rel_orientation -= 360.0;
            }
            // Group windows into same, opposite and adjacent sides
            if op_rel_orientation <= 45.0 {
                openings_same_side.push(op);
            } else if op_rel_orientation >= 135.0 {
                openings_opp_side.push(op);
            }
            // Else opening is on adjacent side, so ignore

            // Assign windows to high and low groups based on which they are closest to
            if mid_height_for(element_from_named(op)).unwrap() < op_height_threshold {
                openings_low.push(op);
            } else {
                openings_high.push(op);
            }
        }

        let cross_vent = !openings_opp_side.is_empty();

        let mut a_b = None;
        let mut a_w = None;
        let mut stack_vent = false;
        let mut opening_height_diff = 0.0;
        let mut opening_area_ratio = None;
        let mut openings_for_struct = None;

        if cross_vent {
            let unique_openings_same_side = openings_same_side.iter().collect::<HashSet<_>>();
            let unique_openings_opp_side = openings_opp_side.iter().collect::<HashSet<_>>();
            let unique_openings_high = openings_high.iter().collect::<HashSet<_>>();
            let unique_openings_low = openings_low.iter().collect::<HashSet<_>>();
            let openings_same_side_high = unique_openings_same_side
                .intersection(&unique_openings_high)
                .collect::<Vec<_>>();
            let openings_same_side_low = unique_openings_same_side
                .intersection(&unique_openings_low)
                .collect::<Vec<_>>();
            let openings_opp_side_high = unique_openings_opp_side
                .intersection(&unique_openings_high)
                .collect::<Vec<_>>();
            let openings_opp_side_low = unique_openings_same_side
                .intersection(&unique_openings_low)
                .collect::<Vec<_>>();

            // Calculate high and low opening areas on same and opposite sides of building
            let a1 = opening_area_equiv_total_ratio
                * openings_same_side_high
                    .iter()
                    .map(|nel| area_for_building_element_input(element_from_named(nel)))
                    .sum::<f64>();
            let a2 = opening_area_equiv_total_ratio
                * openings_same_side_low
                    .iter()
                    .map(|nel| area_for_building_element_input(element_from_named(nel)))
                    .sum::<f64>();
            let a3 = opening_area_equiv_total_ratio
                * openings_opp_side_high
                    .iter()
                    .map(|nel| area_for_building_element_input(element_from_named(nel)))
                    .sum::<f64>();
            let a4 = opening_area_equiv_total_ratio
                * openings_opp_side_low
                    .iter()
                    .map(|nel| area_for_building_element_input(element_from_named(nel)))
                    .sum::<f64>();
            a_w = Some((1.0 / ((1.0 / ((a1 + a3).powi(2))) + 1.0 / ((a2 + a4).powi(2)))).sqrt());
            if a2 + a4 == 0.0 {
                a_b = Some(0.0);
                opening_height_diff = 0.0;
            } else {
                a_b =
                    Some((1.0 / ((1.0 / ((a1 + a3).powi(2))) + 1.0 / ((a2 + a4).powi(2)))).sqrt());

                // Calculate area-weighted average height of windows in high and low groups
                let opening_mid_height_ave_upper = (openings_same_side_high
                    .iter()
                    .map(|nel| {
                        let el = element_from_named(nel);
                        mid_height_for(el).unwrap() * area_for_building_element_input(el)
                    })
                    .sum::<f64>()
                    + openings_opp_side_high
                        .iter()
                        .map(|nel| {
                            let el = element_from_named(nel);
                            mid_height_for(el).unwrap() * area_for_building_element_input(el)
                        })
                        .sum::<f64>())
                    / (openings_same_side_high
                        .iter()
                        .map(|nel| area_for_building_element_input(element_from_named(nel)))
                        .sum::<f64>()
                        + openings_opp_side_high
                            .iter()
                            .map(|nel| area_for_building_element_input(element_from_named(nel)))
                            .sum::<f64>());
                let opening_mid_height_ave_lower = (openings_same_side_low
                    .iter()
                    .map(|nel| {
                        let el = element_from_named(nel);
                        mid_height_for(el).unwrap() * area_for_building_element_input(el)
                    })
                    .sum::<f64>()
                    + openings_opp_side_low
                        .iter()
                        .map(|nel| {
                            let el = element_from_named(nel);
                            mid_height_for(el).unwrap() * area_for_building_element_input(el)
                        })
                        .sum::<f64>())
                    / (openings_same_side_low
                        .iter()
                        .map(|nel| area_for_building_element_input(element_from_named(nel)))
                        .sum::<f64>()
                        + openings_opp_side_low
                            .iter()
                            .map(|nel| area_for_building_element_input(element_from_named(nel)))
                            .sum::<f64>());
                opening_height_diff = opening_mid_height_ave_upper - opening_mid_height_ave_lower;
            }
        } else if openings_high.len() > 1 && openings_low.len() > 1 {
            stack_vent = true;
            let opening_area_upper = openings_high
                .iter()
                .map(|nel| area_for_building_element_input(element_from_named(nel)))
                .sum::<f64>();
            let opening_area_lower = openings_low
                .iter()
                .map(|nel| area_for_building_element_input(element_from_named(nel)))
                .sum::<f64>();

            // Calculate opening area ratio
            opening_area_ratio = Some(opening_area_upper / opening_area_lower);

            let opening_mid_height_ave_upper = openings_high
                .iter()
                .map(|nel| {
                    let el = element_from_named(nel);
                    mid_height_for(el).unwrap() * area_for_building_element_input(el)
                })
                .sum::<f64>()
                / opening_area_upper;
            let opening_mid_height_ave_lower = openings_low
                .iter()
                .map(|nel| {
                    let el = element_from_named(nel);
                    mid_height_for(el).unwrap() * area_for_building_element_input(el)
                })
                .sum::<f64>()
                / opening_area_lower;
            // Calculate opening height difference
            opening_height_diff = opening_mid_height_ave_upper - opening_mid_height_ave_lower;
        } else {
            stack_vent = false;
            openings_for_struct = Some(openings);
        }

        Self {
            window_area_equivalent,
            external_conditions,
            openings: openings_for_struct,
            control,
            natural_ventilation,
            a_b,
            a_w,
            opening_height_diff,
            opening_area_ratio,
            cross_vent,
            stack_vent,
        }
    }

    pub fn h_ve_max(&self, zone_volume: f64, temp_int: f64, timestep_idx: usize) -> f64 {
        let wind_speed = self
            .external_conditions
            .wind_speed_for_timestep_idx(timestep_idx);
        let temp_ext = self
            .external_conditions
            .air_temp_for_timestep_idx(timestep_idx);
        let temp_diff = (temp_int - temp_ext).abs();
        let temp_average_c = (temp_int + temp_ext) / 2.0;
        let temp_average_k = celsius_to_kelvin(temp_average_c);

        let mut q_v = 0.0;

        if self.cross_vent {
            let q_v_wind = C_D * self.a_w.unwrap() * wind_speed * DC_P.powf(0.5);
            let q_v_stack = C_D
                * self.a_b.unwrap()
                * ((2.0 * temp_diff * self.opening_height_diff * G) / temp_average_k).powf(0.5);
            if wind_speed / temp_diff.sqrt()
                < 0.26
                    * (self.a_b.unwrap() / self.a_w.unwrap())
                    * (self.opening_height_diff * DC_P).powf(0.5)
            {
                q_v = q_v_stack;
            } else {
                q_v = q_v_wind;
            }
        } else {
            let q_v_wind = 0.025 * self.window_area_equivalent * wind_speed;
            let mut q_v_stack = 0.0;
            if self.stack_vent {
                q_v_stack = C_D
                    * self.window_area_equivalent
                    * (self.opening_area_ratio.unwrap() * 2.0f64.sqrt()
                        / (1.0 + self.opening_area_ratio.unwrap())
                        * (1.0
                            + self
                                .opening_area_ratio
                                .expect("opening area ratio should be set for natural ventilation")
                                .powi(2))
                        .powf(0.5))
                    * ((temp_diff * self.opening_height_diff * G) / temp_average_k).powf(0.5);
            } else {
                // Stack effect is only between top and bottom of each opening
                if let Some(openings) = &self.openings {
                    for op in openings.iter() {
                        q_v_stack += C_D * self.window_area_equivalent / 3.0
                            * ((temp_diff
                                * projected_height_for_transparent_element(op).unwrap()
                                * G)
                                / temp_average_k)
                                .powf(0.5);
                    }
                }
            }
            q_v = match q_v_wind > q_v_stack {
                true => q_v_wind,
                false => q_v_stack,
            };
        }

        // Calculate h_ve according to BS EN ISO 52016-1:2017 section 6.5.10 equation 61
        let h_ve = P_A * C_A * q_v;

        // Calculate max h_ve achievable for window opening
        let h_ve_nat_vent = match &self.natural_ventilation {
            Some(nv) => nv.h_ve_heat_transfer_coefficient(zone_volume, None, Some(timestep_idx)),
            None => 0.0,
        };

        match 0.0 > h_ve - h_ve_nat_vent {
            true => 0.0,
            false => h_ve - h_ve_nat_vent,
        }
    }
}

/// Calculate the supply temperature of the air flow element
///         according to ISO 52016-1:2017, Section 6.5.10.2
pub fn temp_supply_for_window_opening(
    window_opening: &WindowOpeningForCooling,
    timestep_idx: usize,
) -> f64 {
    window_opening
        .external_conditions
        .air_temp_for_timestep_idx(timestep_idx)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::core::energy_supply::energy_supply::EnergySupply;
    use crate::external_conditions::{DaylightSavingsConfig, ExternalConditions};
    use crate::input::EnergySupplyType;
    use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
    use rstest::*;

    #[fixture]
    pub fn simulation_time_iterator() -> SimulationTimeIterator {
        SimulationTime::new(0.0, 8.0, 1.0).iter()
    }

    #[fixture]
    pub fn external_conditions(
        simulation_time_iterator: SimulationTimeIterator,
    ) -> ExternalConditions {
        let air_temps = vec![0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5];
        let wind_speeds = vec![3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4];

        ExternalConditions::new(
            simulation_time_iterator,
            air_temps,
            wind_speeds,
            vec![0.0; 8],
            vec![0.0; 8],
            vec![],
            55.0,
            0.0,
            0,
            0,
            None,
            1.0,
            None,
            DaylightSavingsConfig::NotApplicable,
            false,
            false,
            vec![],
        )
    }

    #[fixture]
    pub fn infiltration_element(
        external_conditions: ExternalConditions,
    ) -> VentilationElementInfiltration {
        VentilationElementInfiltration::new(
            1,
            InfiltrationShelterType::Sheltered,
            InfiltrationBuildType::House,
            4.5,
            InfiltrationTestType::FiftyPascals,
            40.0,
            75.0,
            2,
            2,
            2,
            1,
            0,
            0,
            0,
            3,
            6,
            0,
            external_conditions,
            None,
        )
    }

    #[rstest]
    pub fn should_have_correct_infiltration_rate_from_openings(
        infiltration_element: VentilationElementInfiltration,
    ) {
        assert_eq!(
            infiltration_element.infiltration_rate_from_openings, 4.0,
            "incorrect infiltration rate for openings returned"
        )
    }

    #[rstest]
    pub fn should_have_correct_q50_divisor(infiltration_element: VentilationElementInfiltration) {
        assert_eq!(
            infiltration_element.q50_divisor(),
            30.7,
            "incorrect Q50 divisor returned"
        )
    }

    #[rstest]
    pub fn should_have_correct_shelter_factor(
        infiltration_element: VentilationElementInfiltration,
    ) {
        assert_eq!(
            infiltration_element.shelter_factor(),
            0.85,
            "incorrect shelter factor returned"
        )
    }

    #[rstest]
    pub fn should_have_correct_infiltration_rate(
        infiltration_element: VentilationElementInfiltration,
    ) {
        assert_eq!(
            round_by_precision(infiltration_element.infiltration_rate(), 1e6),
            round_by_precision(3.5465798045602606, 1e6),
            "incorrect infiltration rate returned"
        )
    }

    fn round_by_precision(src: f64, precision: f64) -> f64 {
        (precision * src).round() / precision
    }

    const EXPECTED_HEAT_TRANSFER_COEFFICIENTS: [f64; 8] = [
        82.78176841476655,
        85.01911350705754,
        87.25645859934853,
        89.49380369163951,
        91.7311487839305,
        93.9684938762215,
        96.20583896851247,
        98.4431840608035,
    ];

    #[rstest]
    pub fn should_have_correct_h_ve_heat_transfer_coefficient(
        infiltration_element: VentilationElementInfiltration,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        for simtime_step in simulation_time_iterator {
            // external_conditions.next();
            assert_eq!(
                round_by_precision(
                    infiltration_element.h_ve_heat_transfer_coefficient(
                        75.0,
                        None,
                        Some(simtime_step.index)
                    ),
                    1e6
                ),
                round_by_precision(EXPECTED_HEAT_TRANSFER_COEFFICIENTS[simtime_step.index], 1e6),
                "incorrect heat transfer coeffient (h_ve) returned"
            )
        }
    }

    #[rstest]
    pub fn should_have_correct_temp_supply(
        infiltration_element: VentilationElementInfiltration,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        for simtime_step in simulation_time_iterator {
            assert_eq!(
                infiltration_element.temp_supply(simtime_step.index),
                2.5 * simtime_step.index as f64,
                "incorrect external temperature returned on iteration {} (1-indexed)",
                simtime_step.index + 1
            )
        }
    }

    #[fixture]
    pub fn energy_supply(simulation_time_iterator: SimulationTimeIterator) -> EnergySupply {
        EnergySupply::new(
            EnergySupplyType::Electricity,
            simulation_time_iterator,
            None,
        )
    }

    #[fixture]
    pub fn mvhr(
        external_conditions: ExternalConditions,
        mut energy_supply: EnergySupply,
        simulation_time_iterator: SimulationTimeIterator,
    ) -> MechanicalVentilationHeatRecovery {
        MechanicalVentilationHeatRecovery::new(
            0.5,
            2.0,
            0.66,
            energy_supply,
            "MVHR".to_string(),
            external_conditions,
            simulation_time_iterator,
        )
    }

    #[rstest]
    pub fn should_have_correct_h_ve_for_mechanical(mvhr: MechanicalVentilationHeatRecovery) {
        // for (i, _) in simulation_time_iterator.enumerate() {
        assert_eq!(
            round_by_precision(mvhr.h_ve_heat_transfer_coefficient(75.0, None, None), 1e6),
            round_by_precision(4.28975166666666, 1e6),
        );
        assert_eq!(
            round_by_precision(
                mvhr.h_ve_heat_transfer_coefficient(75.0, Some(1.2), None),
                1e6
            ),
            round_by_precision(5.147701999999999, 1e6),
        );
        // }
    }

    #[rstest]
    pub fn should_have_correct_fan_gains(
        mut mvhr: MechanicalVentilationHeatRecovery,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        for (i, _) in simulation_time_iterator.enumerate() {
            assert_eq!(
                round_by_precision(mvhr.fans(75.0, i, None), 1e6),
                round_by_precision(0.010416666666666666, 1e6),
                "incorrect fan gains for MVHR on iteration {} (1-indexed)",
                i + 1
            );
            // assert_eq!(mvhr.energy_supply().results_by_end_user().get("MVHR".to_string())) // result_by_end_user not yet implemented
        }
    }

    #[rstest]
    pub fn should_have_correct_temp_supply_for_mechanical(
        mvhr: MechanicalVentilationHeatRecovery,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        for (i, _) in simulation_time_iterator.enumerate() {
            assert_eq!(
                mvhr.temp_supply(i),
                i as f64 * 2.5,
                "incorrect supply temp returned"
            );
        }
    }

    #[fixture]
    pub fn whole_house_extract_ventilation(
        external_conditions: ExternalConditions,
        mut energy_supply: EnergySupply,
        simulation_time_iterator: SimulationTimeIterator,
    ) -> WholeHouseExtractVentilation {
        WholeHouseExtractVentilation::new(
            0.5,
            2.0,
            0.25,
            energy_supply,
            "WHEV".to_string(),
            external_conditions,
            simulation_time_iterator,
        )
    }

    const WHOLE_HOUSE_H_VE_RESULTS: [f64; 8] = [
        8.131011373697916,
        8.047227161458334,
        7.965414342447915,
        7.885572916666667,
        7.807702884114583,
        7.731804244791666,
        7.657876998697917,
        7.585921145833332,
    ];
    const WHOLE_HOUSE_H_VE_RESULTS_WITH_THROUGHPUT_FACTOR: [f64; 8] = [
        9.757213648437498,
        9.656672593749999,
        9.5584972109375,
        9.4626875,
        9.3692434609375,
        9.27816509375,
        9.1894523984375,
        9.103105374999997,
    ];

    #[rstest]
    pub fn should_have_correct_h_ve_for_whole_house(
        whole_house_extract_ventilation: WholeHouseExtractVentilation,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        for (i, _) in simulation_time_iterator.enumerate() {
            assert_eq!(
                round_by_precision(
                    whole_house_extract_ventilation.h_ve_heat_transfer_coefficient(
                        75.0,
                        None,
                        Some(i)
                    ),
                    1e6
                ),
                round_by_precision(WHOLE_HOUSE_H_VE_RESULTS[i], 1e6),
                "incorrect heat transfer coefficient (h_ve) returned for iteration {} (1-indexed)",
                i + 1
            );
            assert_eq!(
                round_by_precision(
                    whole_house_extract_ventilation.h_ve_heat_transfer_coefficient(
                        75.0,
                        Some(1.2),
                        Some(i)
                    ),
                    1e6
                ),
                round_by_precision(WHOLE_HOUSE_H_VE_RESULTS_WITH_THROUGHPUT_FACTOR[i], 1e6),
                "incorrect heat transfer coefficient with throughput factor (h_ve) returned for iteration {} (1-indexed)",
                i + 1
            );
        }
    }

    #[rstest]
    pub fn should_have_correct_fan_gains_for_whole_house(
        mut whole_house_extract_ventilation: WholeHouseExtractVentilation,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        for (i, _) in simulation_time_iterator.enumerate() {
            assert_eq!(
                round_by_precision(whole_house_extract_ventilation.fans(75.0, None, i), 1e6),
                round_by_precision(0.0, 1e6),
            );
            // skip test of results_by_end_user from the energy supply as this is not yet implemented
        }
    }

    #[rstest]
    pub fn should_have_correct_temp_supply_for_whole_house(
        whole_house_extract_ventilation: WholeHouseExtractVentilation,
        simulation_time_iterator: SimulationTimeIterator,
    ) {
        for (i, _) in simulation_time_iterator.enumerate() {
            assert_eq!(
                whole_house_extract_ventilation.temp_supply(i),
                i as f64 * 2.5,
                "incorrect supply temp returned"
            );
        }
    }
}
