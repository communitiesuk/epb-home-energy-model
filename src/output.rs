use crate::{EnergySupplyStatKey, StringOrNumber};
use indexmap::IndexMap;
use itertools::Itertools;
use serde::Serialize;
use serde_enum_str::{Deserialize_enum_str, Serialize_enum_str};
use smartstring::alias::String;
use std::sync::Arc;

#[derive(Debug, Serialize)]
pub struct OutputStatic {
    /// Heat transfer coefficient (unit: W/K)
    pub(crate) heat_transfer_coefficient: f64,
    /// Heat loss parameter (unit: W/m².K)
    pub(crate) heat_loss_param: f64,
    /// Heat capacity parameter (unit: kJ/m².K)
    pub(crate) heat_capacity_param: f64,
    /// Heat loss form factor
    pub(crate) heat_loss_form_factor: f64,
    /// Internal air temperature (unit: ˚C)
    pub(crate) temperature_air_internal: f64,
    /// External air temperature (unit: ˚C)
    pub(crate) temperature_air_external: f64,
}

/// Zone output data. Each field is keyed by zone name.
/// These fields names are used in the core output file (substituting " " for "_")
#[derive(Debug, Serialize)]
pub struct OutputZoneData {
    // NB: Field order is important as it affects the output files.
    /// Internal gains (unit: W)
    pub(crate) internal_gains: IndexMap<Arc<str>, Vec<f64>>,
    /// Solar gains (unit: W)
    pub(crate) solar_gains: IndexMap<Arc<str>, Vec<f64>>,
    /// Operative temperature (unit: ˚C)
    pub(crate) operative_temp: IndexMap<Arc<str>, Vec<f64>>,
    /// Internal air temperature (unit: ˚C)
    pub(crate) internal_air_temp: IndexMap<Arc<str>, Vec<f64>>,
    /// Space heat demand (unit: kWh)
    pub(crate) space_heat_demand: IndexMap<Arc<str>, Vec<f64>>,
    /// Space cool demand (unit: kWh)
    pub(crate) space_cool_demand: IndexMap<Arc<str>, Vec<f64>>,
}

impl OutputZoneData {
    pub(crate) fn ordered_values_for_zone(&self, zone: &str) -> [&[f64]; 6] {
        [
            &self.internal_gains[zone],
            &self.solar_gains[zone],
            &self.operative_temp[zone],
            &self.internal_air_temp[zone],
            &self.space_heat_demand[zone],
            &self.space_cool_demand[zone],
        ]
    }
}

pub(crate) const OUTPUT_ZONE_DATA_FIELD_HEADINGS: &[&str] = &[
    "internal gains",
    "solar gains",
    "operative temp",
    "internal air temp",
    "space heat demand",
    "space cool demand",
];

/// Contains the output timeseries for the heating and cooling systems, respectively.
#[derive(Debug, Serialize)]
pub struct OutputHeatingCoolingSystem {
    /// Heating system output, for each heating system (unit: kWh)
    pub(crate) heating_system_output: IndexMap<Option<Arc<str>>, Vec<f64>>,
    /// Cooling system output, keyed by cooling system name (unit: kWh)
    pub(crate) cooling_system_output: IndexMap<Option<Arc<str>>, Vec<f64>>,
}

/// Hot water systems data for every time step.
#[derive(Debug, Serialize)]
pub struct OutputHotWaterSystems {
    // NB: Field order is important as it affects the output files.
    // Each field's alias is the heading to use in the core output CSV.
    /// Hot water volume required from hot water source, for each hot water source (unit: litres)
    pub(crate) demand: IndexMap<Arc<str>, Vec<f64>>,
    /// Hot water energy demand at hot water source, for each hot water source (unit: kWh)
    pub(crate) energy_demand_at_hot_water_source: IndexMap<Arc<str>, Vec<f64>>,
    /// Hot water energy demand at connected tapping points, for each hot water source (unit: kWh)
    pub(crate) energy_demand_at_tapping_points: IndexMap<Arc<str>, Vec<f64>>,
    /// Total hot water event duration, for each hot water source (unit: minutes)
    pub(crate) duration: IndexMap<Arc<str>, Vec<f64>>,
    /// Number of hot water events, for each hot water source (unit: count)
    pub(crate) events_count: IndexMap<Arc<str>, Vec<f64>>,
    /// Pipework losses, for each hot water source (unit: kWh)
    pub(crate) losses_pipework: IndexMap<Arc<str>, Vec<f64>>,
    /// Primary pipework losses, for each hot water source (unit: kWh)
    pub(crate) losses_primary_pipework: IndexMap<Arc<str>, Vec<f64>>,
    /// Storage losses, for each hot water source (unit: kWh)
    pub(crate) losses_storage: IndexMap<Arc<str>, Vec<f64>>,
}

#[derive(Clone, Copy, Deserialize_enum_str, Debug, Eq, Hash, PartialEq, Serialize_enum_str)]
pub(crate) enum OutputHotWaterSystemsAlias {
    #[serde(rename = "hot water volume required from hot water source")]
    VolumeRequired,
    #[serde(rename = "hot water energy demand at hot water source")]
    EnergyDemandSource,
    #[serde(rename = "hot water energy demand at connected tapping points")]
    EnergyDemandTappingPoints,
    #[serde(rename = "total event duration")]
    TotalEventDuration,
    #[serde(rename = "number of events")]
    NumberOfEvents,
    #[serde(rename = "distribution pipework losses")]
    DistributionPipeworkLosses,
    #[serde(rename = "primary pipework losses")]
    PrimaryPipeworkLosses,
    #[serde(rename = "storage losses")]
    StorageLosses,
}

impl OutputHotWaterSystems {
    pub(crate) fn fields(&self) -> [(OutputHotWaterSystemsAlias, Vec<Arc<str>>); 8] {
        [
            (
                OutputHotWaterSystemsAlias::VolumeRequired,
                self.demand.keys().cloned().collect_vec(),
            ),
            (
                OutputHotWaterSystemsAlias::EnergyDemandSource,
                self.energy_demand_at_hot_water_source
                    .keys()
                    .cloned()
                    .collect_vec(),
            ),
            (
                OutputHotWaterSystemsAlias::EnergyDemandTappingPoints,
                self.energy_demand_at_tapping_points
                    .keys()
                    .cloned()
                    .collect_vec(),
            ),
            (
                OutputHotWaterSystemsAlias::TotalEventDuration,
                self.duration.keys().cloned().collect_vec(),
            ),
            (
                OutputHotWaterSystemsAlias::NumberOfEvents,
                self.events_count.keys().cloned().collect_vec(),
            ),
            (
                OutputHotWaterSystemsAlias::DistributionPipeworkLosses,
                self.losses_pipework.keys().cloned().collect_vec(),
            ),
            (
                OutputHotWaterSystemsAlias::PrimaryPipeworkLosses,
                self.losses_primary_pipework.keys().cloned().collect_vec(),
            ),
            (
                OutputHotWaterSystemsAlias::StorageLosses,
                self.losses_storage.keys().cloned().collect_vec(),
            ),
        ]
    }

    pub(crate) fn ordered_values(&self) -> [&IndexMap<Arc<str>, Vec<f64>>; 8] {
        [
            &self.demand,
            &self.energy_demand_at_hot_water_source,
            &self.energy_demand_at_tapping_points,
            &self.duration,
            &self.events_count,
            &self.losses_pipework,
            &self.losses_primary_pipework,
            &self.losses_storage,
        ]
    }
}

/// Coefficients of performance for space heating & cooling, and hot water systems.
/// Null-values are used when a divide-by-zero occurs in the results. These are output as "DIV/0" in the CSVs.
#[derive(Debug, Serialize)]
pub struct OutputCop {
    /// Overall coefficient of performance for each heating system (unitless)
    pub(crate) space_heating_system: IndexMap<Arc<str>, Option<f64>>,
    /// Overall coefficient of performance for each heating system (unitless)
    pub(crate) space_cooling_system: IndexMap<Arc<str>, Option<f64>>,
    /// Overall coefficient of performance for each heating system (unitless)
    pub(crate) hot_water_system: IndexMap<Arc<str>, Option<f64>>,
}

/// Emitters data for every time step.
#[derive(Debug, Serialize)]
pub struct OutputEmitters {
    /// Current time step
    pub(crate) simulation_time_idx: usize,
    /// Energy demand (unit: kWh)
    pub(crate) energy_demand: f64,
    /// Required emitter temperature to satisfy zone setpoint (unit: Celsius)
    pub(crate) temp_emitter_required: f64,
    /// Time at which the emitter begins operating
    pub(crate) time_heating_start: f64,
    /// Energy provided by heat source (unit: kWh)
    pub(crate) energy_provided_by_heat_source: f64,
    /// Temperature of the emitter (unit: Celsius)
    pub(crate) temp_emitter: f64,
    /// Maximum temperature the emitter can safely reach (unit: Celsius)
    pub(crate) temp_emitter_max: f64,
    /// Energy released from emitters (unit: kWh)
    pub(crate) energy_released_from_emitters: f64,
    /// Desired supply temperature to the emitter (unit: Celsius)
    pub(crate) temp_flow_target: f64,
    /// Desired return temperature from the emitter (unit: Celsius)
    pub(crate) temp_return_target: f64,
    /// Whether the maximum emitter temperature is used as the final supply temperature
    pub(crate) temp_emitter_max_is_final_temp: bool,
    /// Energy required from the heat source (unit: kWh)
    pub(crate) energy_required_from_heat_source: f64,
    /// Fan electricity consumption (unit: kWh)
    #[serde(rename = "fan_energy_kWh")]
    pub(crate) fan_energy_kwh: f64,
}

/// The raw HEM Core outputs, from Corpus::run()
///
/// For the _unmet_demand energy supply object, end-uses are:
///        - space heating or cooling per zone, where the naming convention is _unmet_demand: [zone name from input]
///       - water heating unmet demand, where the naming convention is _unmet_demand: [hot water source names from input]
///
/// For all other energy supplies, end-use names refer to specific systems or appliances specified in the input file.
#[derive(Debug, Serialize)]
pub struct OutputCore {
    /// The list of timesteps, serves as an index for the other outputs (unit: Hours)
    pub(crate) timestep_array: Vec<f64>,
    /// Total energy (unit: kWh)
    pub(crate) results_totals: IndexMap<Arc<str>, Vec<f64>>,
    /// Energy per supply, per end use (unit: kWh)
    pub(crate) results_end_user: IndexMap<Arc<str>, IndexMap<Arc<str>, Vec<f64>>>,
    /// Energy import (unit: kWh)
    pub(crate) energy_import: IndexMap<Arc<str>, Vec<f64>>,
    /// Energy export (unit: kWh)
    pub(crate) energy_export: IndexMap<Arc<str>, Vec<f64>>,
    /// Energy from grid to consumption (unit: kWh)
    pub(crate) grid_to_consumption: IndexMap<Arc<str>, Vec<f64>>,
    /// Energy export from generation (unit: kWh)
    pub(crate) generation_to_grid: IndexMap<Arc<str>, Vec<f64>>,
    /// Energy generated and consumed (unit: kWh)
    pub(crate) energy_generated_consumed: IndexMap<Arc<str>, Vec<f64>>,
    /// Energy to storage (unit: kWh)
    pub(crate) energy_to_storage: IndexMap<Arc<str>, Vec<f64>>,
    /// Energy from storage (unit: kWh)
    pub(crate) energy_from_storage: IndexMap<Arc<str>, Vec<f64>>,
    /// Imported energy to storage, (unit: kWh)
    pub(crate) storage_from_grid: IndexMap<Arc<str>, Vec<f64>>,
    /// Battery charge level (unit: ratio 0 to 1 as a percentage)
    pub(crate) battery_state_of_charge: IndexMap<Arc<str>, Vec<f64>>,
    /// Energy diverted (unit: kWh)
    pub(crate) energy_diverted: IndexMap<Arc<str>, Vec<f64>>,
    /// Energy supply beta factor (unit: ratio 0 to 1)
    pub(crate) beta_factor: IndexMap<Arc<str>, Vec<f64>>,
    /// List of the unique zone names in the zone data
    pub(crate) zone_list: Vec<Arc<str>>,
    pub(crate) zone_data: OutputZoneData,
    pub(crate) heating_cooling_system: OutputHeatingCoolingSystem,
    pub(crate) cop: OutputCop,
    /// Ventilation ductwork gains (unit: kWh)
    pub(crate) ductwork_gains: Vec<f64>,
    /// Heat balance data for each zone.
    pub(crate) heat_balance_all:
        IndexMap<Arc<str>, IndexMap<Arc<str>, IndexMap<Arc<str>, Vec<f64>>>>,
    /// Heat source wet detailed results.
    pub(crate) heat_source_wet_results: IndexMap<Arc<str>, IndexMap<Arc<str>, Vec<f64>>>,
    /// Annual heat source wet detailed results.
    pub(crate) heat_source_wet_results_annual: IndexMap<Arc<str>, IndexMap<Arc<str>, f64>>,
    /// Hot water source results summary.
    /// Currently unstructured, see CSV for column-order.
    pub(crate) hot_water_source_results_summary:
        IndexMap<Arc<str>, Vec<Vec<Option<StringOrNumber>>>>,
    /// Heating system emitters detailed outputs.
    /// Currently unstructured, see CSV for column-order.
    pub(crate) emitters: IndexMap<Arc<str>, IndexMap<usize, OutputEmitters>>,
    /// Electric storage heaters detailed outputs.
    /// Currently unstructured, see CSV for column-order.
    pub(crate) electric_storage_heaters: IndexMap<Arc<str>, IndexMap<usize, Vec<f64>>>,
    /// Ventilation detailed outputs.
    /// Currently unstructured, see CSV for column-order.
    pub(crate) ventilation: Vec<Vec<StringOrNumber>>,
    pub(crate) hot_water_systems: OutputHotWaterSystems,
}

/// Details of the time-step with the peak electricity consumption within the simulation.
#[derive(Debug, Serialize)]
pub struct OutputSummaryPeakElectricityConsumption {
    /// The maximum electricity consumption (unit: kWh)
    pub(crate) peak: f64,
    /// The time step the peak occurred at
    pub(crate) index: usize,
    /// The month the peak occurred in
    pub(crate) month: u8,
    /// The day of the month the peak occurred on
    pub(crate) day: u8,
    /// The hour of the day the peak occurred on. Decimal-time is used if the simulation step is less than 1 hour.
    pub(crate) hour: u8,
}

#[derive(Debug, Serialize)]
pub struct OutputSummaryEnergySupply {
    /// All energy generated on site (unit: kWh)
    pub(crate) electricity_generated: f64,
    /// All energy consumed for space and water heating, appliances and cooking (unit: kWh)
    pub(crate) electricity_consumed: f64,
    /// Total energy generated and consumed immediately, excluding PV diverter (unit: kWh)
    pub(crate) generation_to_consumption: f64,
    /// Total energy exported to grid immediately from generation (unit: kWh)
    pub(crate) generation_to_grid: f64,
    /// Total energy generated and sent to PV diverter (unit: kWh)
    pub(crate) generation_to_diverter: f64,
    /// Total energy imported from the grid (unit: kWh)
    pub(crate) grid_to_consumption: f64,
    /// Total energy imported and put into storage (unit: kWh)
    pub(crate) grid_to_storage: f64,
    /// Total energy generated and put into storage (unit: kWh)
    pub(crate) generation_to_storage: f64,
    /// Total energy consumed from storage (unit: kWh)
    pub(crate) storage_to_consumption: f64,
    /// Storage round-trip efficiency.
    /// Total energy consumed from storage divided by total energy put into storage, from both grid and on site generation.
    /// (unit: kWh)
    pub(crate) storage_efficiency: f64,
    /// Net import grid to consumption minus generation to grid (unit: kWh)
    pub(crate) net_import: f64,
    /// Gross import from grid (unit: kWh)
    pub(crate) total_gross_import: f64,
    /// Gross export to grid (unit: kWh)
    pub(crate) total_gross_export: f64,
}

impl OutputSummaryEnergySupply {
    pub(crate) fn field(&self, key: &EnergySupplyStatKey) -> f64 {
        match key {
            EnergySupplyStatKey::ElectricityGenerated => self.electricity_generated,
            EnergySupplyStatKey::ElectricityConsumed => self.electricity_consumed,
            EnergySupplyStatKey::GenerationToConsumption => self.generation_to_consumption,
            EnergySupplyStatKey::GridToConsumption => self.grid_to_consumption,
            EnergySupplyStatKey::GenerationToGrid => self.generation_to_grid,
            EnergySupplyStatKey::NetImport => self.net_import,
            EnergySupplyStatKey::GenerationToStorage => self.generation_to_storage,
            EnergySupplyStatKey::StorageToConsumption => self.storage_to_consumption,
            EnergySupplyStatKey::StorageFromGrid => self.grid_to_storage,
            EnergySupplyStatKey::GenerationToDiverter => self.generation_to_diverter,
            EnergySupplyStatKey::StorageEfficiency => self.storage_efficiency,
        }
    }
}

/// Summary outputs from HEM Core.
#[derive(Debug, Serialize)]
pub struct OutputSummary {
    /// Total floor-area of all zones (unit: m²)
    pub(crate) total_floor_area: f64,
    /// Space heating demand total (unit: kWh)
    pub(crate) space_heat_demand_total: f64,
    /// Space cooling demand total (unit: kWh)
    pub(crate) space_cool_demand_total: f64,
    pub(crate) electricity_peak_consumption: OutputSummaryPeakElectricityConsumption,
    pub(crate) energy_supply: IndexMap<Arc<str>, OutputSummaryEnergySupply>,
    /// Delivered energy summary, total energy per fuel and end-use (unit: kWh)
    pub(crate) delivered_energy: IndexMap<Arc<str>, IndexMap<Arc<str>, f64>>,
    /// 75th percentile of hot water demand summed over each 24 hour segment of the simulation.
    pub(crate) hot_water_demand_daily_75th_percentile: IndexMap<Arc<str>, f64>,
}

impl OutputSummary {
    /// Space heating demand total per floor-area (unit: kWh/m²)
    pub(crate) fn space_heat_demand_by_floor_area(&self) -> f64 {
        self.space_heat_demand_total / self.total_floor_area
    }

    /// Space cooling demand total per floor-area (unit: kWh/m²)
    pub(crate) fn space_cool_demand_by_floor_area(&self) -> f64 {
        self.space_cool_demand_total / self.total_floor_area
    }

    pub(crate) fn delivered_energy_by_floor_area(
        &self,
    ) -> IndexMap<Arc<str>, IndexMap<Arc<str>, f64>> {
        self.delivered_energy
            .iter()
            .map(|(fuel, end_use_map)| {
                (
                    fuel.clone(),
                    end_use_map
                        .iter()
                        .map(|(end_use, energy)| (end_use.clone(), energy / self.total_floor_area))
                        .collect(),
                )
            })
            .collect()
    }
}

#[derive(Debug, Serialize)]
pub struct OutputMetadata {
    pub(crate) hem_core_version: String,
}

#[derive(Debug, Serialize)]
pub struct Output {
    pub(crate) static_: OutputStatic,
    pub(crate) core: OutputCore,
    pub(crate) summary: OutputSummary,
    pub(crate) metadata: OutputMetadata,
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use rstest::*;

    #[rstest]
    fn test_output_summary_computed_fields() {
        let mut base_summary = OutputSummary {
            total_floor_area: 0.0,
            space_heat_demand_total: 0.0,
            space_cool_demand_total: 0.0,
            electricity_peak_consumption: OutputSummaryPeakElectricityConsumption {
                peak: 0.0,
                index: 0,
                month: 0,
                day: 0,
                hour: 0,
            },
            energy_supply: Default::default(),
            delivered_energy: Default::default(),
            hot_water_demand_daily_75th_percentile: Default::default(),
        };

        let space_heat_demand_data = [
            (100.0, 12_000.0, 120.0),
            (50.0, 12_000.0, 240.0),
            (85.0, 12_000.0, 141.1764705882),
        ];

        for (total_floor_area, space_heat_demand_total, expected) in space_heat_demand_data {
            base_summary.total_floor_area = total_floor_area;
            base_summary.space_heat_demand_total = space_heat_demand_total;
            assert_relative_eq!(
                base_summary.space_heat_demand_by_floor_area(),
                expected,
                epsilon = 1e-7
            );
        }

        let space_cool_demand_data = [
            (100.0, 6_000.0, 60.0),
            (50.0, 6_000.0, 120.0),
            (85.0, 2_000.0, 23.529411764705),
        ];

        for (total_floor_area, space_cool_demand_total, expected) in space_cool_demand_data {
            base_summary.total_floor_area = total_floor_area;
            base_summary.space_cool_demand_total = space_cool_demand_total;
            assert_relative_eq!(
                base_summary.space_cool_demand_by_floor_area(),
                expected,
                epsilon = 1e-7
            );
        }

        let delivered_energy_data: [(
            f64,
            IndexMap<Arc<str>, IndexMap<Arc<str>, f64>>,
            IndexMap<Arc<str>, IndexMap<Arc<str>, f64>>,
        ); 1] = [(
            100.0,
            IndexMap::from([
                (
                    "total".into(),
                    IndexMap::from([
                        ("total".into(), 20_000.0),
                        ("_unmet_demand".into(), 0.0),
                        ("cooking".into(), 3_000.0),
                        ("lighting".into(), 500.0),
                    ]),
                ),
                (
                    "mains_elec".into(),
                    IndexMap::from([
                        ("total".into(), 20_000.0),
                        ("_unmet_demand".into(), 0.0),
                        ("cooking".into(), 1_000.0),
                        ("lighting".into(), 500.0),
                    ]),
                ),
            ]),
            IndexMap::from([
                (
                    "total".into(),
                    IndexMap::from([
                        ("total".into(), 200.0),
                        ("_unmet_demand".into(), 0.0),
                        ("cooking".into(), 30.0),
                        ("lighting".into(), 5.0),
                    ]),
                ),
                (
                    "mains_elec".into(),
                    IndexMap::from([
                        ("total".into(), 200.0),
                        ("_unmet_demand".into(), 0.0),
                        ("cooking".into(), 10.0),
                        ("lighting".into(), 5.0),
                    ]),
                ),
            ]),
        )];

        for (total_floor_area, delivered_energy, expected) in delivered_energy_data {
            base_summary.total_floor_area = total_floor_area;
            base_summary.delivered_energy = delivered_energy;
            assert_eq!(base_summary.delivered_energy_by_floor_area(), expected);
        }
    }
}
