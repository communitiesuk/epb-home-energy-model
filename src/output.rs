use crate::StringOrNumber;
use indexmap::IndexMap;
use serde::Serialize;
use smartstring::alias::String;

#[derive(Debug, Serialize)]
struct OutputStatic {
    /// Heat transfer coefficient (unit: W/K)
    heat_transfer_coefficient: f64,
    /// Heat loss parameter (unit: W/m².K)
    heat_loss_param: f64,
    /// Heat capacity parameter (unit: kJ/m².K)
    heat_capacity_param: f64,
    /// Heat loss form factor
    heat_loss_form_factor: f64,
    /// Internal air temperature (unit: ˚C)
    temperature_air_internal: f64,
    /// External air temperature (unit: ˚C)
    temperature_air_external: f64,
}

/// Zone output data. Each field is keyed by zone name.
/// These fields names are used in the core output file (substituting " " for "_")
#[derive(Debug, Serialize)]
struct OutputZoneData {
    // NB: Field order is important as it affects the output files.
    /// Internal gains (unit: W)
    internal_gains: IndexMap<String, Vec<f64>>,
    /// Solar gains (unit: W)
    solar_gains: IndexMap<String, Vec<f64>>,
    /// Operative temperature (unit: ˚C)
    operative_temp: IndexMap<String, Vec<f64>>,
    /// Internal air temperature (unit: ˚C)
    internal_air_temp: IndexMap<String, Vec<f64>>,
    /// Space heat demand (unit: kWh)
    space_heat_demand: IndexMap<String, Vec<f64>>,
    /// Space cool demand (unit: kWh)
    space_cool_demand: IndexMap<String, Vec<f64>>,
}

/// Contains the output timeseries for the heating and cooling systems, respectively.
#[derive(Debug, Serialize)]
struct OutputHeatingCoolingSystem {
    /// Heating system output, for each heating system (unit: kWh)
    heating_system_output: IndexMap<String, Vec<f64>>,
    /// Cooling system output, keyed by cooling system name (unit: kWh)
    cooling_system_output: IndexMap<String, Vec<f64>>,
}

/// Hot water systems data for every time step.
#[derive(Debug, Serialize)]
struct OutputHotWaterSystems {
    // NB: Field order is important as it affects the output files.
    // Each field's alias is the heading to use in the core output CSV.
    /// Hot water volume required from hot water source, for each hot water source (unit: litres)
    #[serde(alias = "hot water volume required from hot water source")]
    demand: IndexMap<String, Vec<f64>>,
    /// Hot water energy demand at hot water source, for each hot water source (unit: kWh)
    #[serde(alias = "hot water energy demand at hot water source")]
    energy_demand_at_hot_water_source: IndexMap<String, Vec<f64>>,
    /// Hot water energy demand at connected tapping points, for each hot water source (unit: kWh)
    #[serde(alias = "hot water energy demand at connected tapping points")]
    energy_demand_at_tapping_points: IndexMap<String, Vec<f64>>,
    /// Total hot water event duration, for each hot water source (unit: minutes)
    #[serde(alias = "total event duration")]
    duration: IndexMap<String, Vec<f64>>,
    /// Number of hot water events, for each hot water source (unit: count)
    #[serde(alias = "number of events")]
    events_count: IndexMap<String, Vec<i32>>,
    /// Pipework losses, for each hot water source (unit: kWh)
    #[serde(alias = "distribution pipework losses")]
    losses_pipework: IndexMap<String, Vec<f64>>,
    /// Primary pipework losses, for each hot water source (unit: kWh)
    #[serde(alias = "primary pipework losses")]
    losses_primary_pipework: IndexMap<String, Vec<f64>>,
    /// Storage losses, for each hot water source (unit: kWh)
    #[serde(alias = "storage losses")]
    losses_storage: IndexMap<String, Vec<f64>>,
}

/// Coefficients of performance for space heating & cooling, and hot water systems.
/// Null-values are used when a divide-by-zero occurs in the results. These are output as "DIV/0" in the CSVs.
#[derive(Debug, Serialize)]
struct OutputCop {
    /// Overall coefficient of performance for each heating system (unitless)
    space_heating_system: IndexMap<String, Option<f64>>,
    /// Overall coefficient of performance for each heating system (unitless)
    space_cooling_system: IndexMap<String, Option<f64>>,
    /// Overall coefficient of performance for each heating system (unitless)
    hot_water_system: IndexMap<String, Option<f64>>,
}

/// Emitters data for every time step.
#[derive(Debug, Serialize)]
struct OutputEmitters {
    /// Current time step
    simulation_time_idx: usize,
    /// Energy demand (unit: kWh)
    energy_demand: f64,
    /// Required emitter temperature to satisfy zone setpoint (unit: Celsius)
    temp_emitter_required: f64,
    /// Time at which the emitter begins operating
    time_heating_start: f64,
    /// Energy provided by heat source (unit: kWh)
    energy_provided_by_heat_source: f64,
    /// Temperature of the emitter (unit: Celsius)
    temp_emitter: f64,
    /// Maximum temperature the emitter can safely reach (unit: Celsius)
    temp_emitter_max: f64,
    /// Energy released from emitters (unit: kWh)
    energy_released_from_emitters: f64,
    /// Desired supply temperature to the emitter (unit: Celsius)
    temp_flow_target: f64,
    /// Desired return temperature from the emitter (unit: Celsius)
    temp_return_target: f64,
    /// Whether the maximum emitter temperature is used as the final supply temperature
    temp_emitter_max_is_final_temp: bool,
    /// Energy required from the heat source (unit: kWh)
    energy_required_from_heat_source: f64,
    /// Fan electricity consumption (unit: kWh)
    #[serde(rename = "fan_energy_kWh")]
    fan_energy_kwh: f64,
}

/// The raw HEM Core outputs, from Corpus::run()
///
/// For the _unmet_demand energy supply object, end-uses are:
///        - space heating or cooling per zone, where the naming convention is _unmet_demand: [zone name from input]
///       - water heating unmet demand, where the naming convention is _unmet_demand: [hot water source names from input]
///
/// For all other energy supplies, end-use names refer to specific systems or appliances specified in the input file.
#[derive(Debug, Serialize)]
struct OutputCore {
    /// The list of timesteps, serves as an index for the other outputs (unit: Hours)
    timestep_array: Vec<f64>,
    /// Total energy (unit: kWh)
    results_totals: IndexMap<String, Vec<f64>>,
    /// Energy per supply, per end use (unit: kWh)
    results_end_user: IndexMap<String, IndexMap<String, Vec<f64>>>,
    /// Energy import (unit: kWh)
    energy_import: IndexMap<String, Vec<f64>>,
    /// Energy export (unit: kWh)
    energy_export: IndexMap<String, Vec<f64>>,
    /// Energy from grid to consumption (unit: kWh)
    grid_to_consumption: IndexMap<String, Vec<f64>>,
    /// Energy export from generation (unit: kWh)
    generation_to_grid: IndexMap<String, Vec<f64>>,
    /// Energy generated and consumed (unit: kWh)
    energy_generated_consumed: IndexMap<String, Vec<f64>>,
    /// Energy to storage (unit: kWh)
    energy_to_storage: IndexMap<String, Vec<f64>>,
    /// Energy from storage (unit: kWh)
    energy_from_storage: IndexMap<String, Vec<f64>>,
    /// Imported energy to storage, (unit: kWh)
    storage_from_grid: IndexMap<String, Vec<f64>>,
    /// Battery charge level (unit: ratio 0 to 1 as a percentage)
    battery_state_of_charge: IndexMap<String, Vec<f64>>,
    /// Energy diverted (unit: kWh)
    energy_diverted: IndexMap<String, Vec<f64>>,
    /// Energy supply beta factor (unit: ratio 0 to 1)
    beta_factor: IndexMap<String, Vec<f64>>,
    /// List of the unique zone names in the zone data
    zone_list: Vec<String>,
    output_zone_data: OutputZoneData,
    heating_cooling_system: OutputHeatingCoolingSystem,
    cop: OutputCop,
    /// Ventilation ductwork gains (unit: kWh)
    ductwork_gains: Vec<f64>,
    /// Heat balance data for each zone.
    heat_balance_all: IndexMap<String, IndexMap<String, Vec<f64>>>,
    /// Heat source wet detailed results.
    heat_source_wet_results: IndexMap<String, IndexMap<String, Vec<f64>>>,
    /// Annual heat source wet detailed results.
    heat_source_wet_results_annual: IndexMap<String, IndexMap<String, f64>>,
    /// Hot water source results summary.
    /// Currently unstructured, see CSV for column-order.
    hot_water_source_results_summary: IndexMap<String, Vec<Vec<Option<StringOrNumber>>>>,
    /// Heating system emitters detailed outputs.
    /// Currently unstructured, see CSV for column-order.
    emitters: IndexMap<String, IndexMap<usize, OutputEmitters>>,
    /// Electric storage heaters detailed outputs.
    /// Currently unstructured, see CSV for column-order.
    electric_storage_heaters: IndexMap<String, IndexMap<usize, Vec<f64>>>,
    /// Ventilation detailed outputs.
    /// Currently unstructured, see CSV for column-order.
    ventilation: Vec<Vec<StringOrNumber>>,
    hot_water_systems: OutputHotWaterSystems,
}

/// Details of the time-step with the peak electricity consumption within the simulation.
#[derive(Debug, Serialize)]
struct OutputSummaryPeakElectricityConsumption {
    /// The maximum electricity consumption (unit: kWh)
    peak: f64,
    /// The time step the peak occurred at
    index: usize,
    /// The month the peak occurred in
    month: u8,
    /// The day of the month the peak occurred on
    day: u8,
    /// The hour of the day the peak occurred on. Decimal-time is used if the simulation step is less than 1 hour.
    hour: u8,
}

#[derive(Debug, Serialize)]
struct OutputSummaryEnergySupply {
    /// All energy generated on site (unit: kWh)
    electricity_generated: f64,
    /// All energy consumed for space and water heating, appliances and cooking (unit: kWh)
    electricity_consumed: f64,
    /// Total energy generated and consumed immediately, excluding PV diverter (unit: kWh)
    generation_to_consumption: f64,
    /// Total energy exported to grid immediately from generation (unit: kWh)
    generation_to_grid: f64,
    /// Total energy generated and sent to PV diverter (unit: kWh)
    generation_to_diverter: f64,
    /// Total energy imported from the grid (unit: kWh)
    grid_to_consumption: f64,
    /// Total energy imported and put into storage (unit: kWh)
    grid_to_storage: f64,
    /// Total energy generated and put into storage (unit: kWh)
    generation_to_storage: f64,
    /// Total energy consumed from storage (unit: kWh)
    storage_to_consumption: f64,
    /// Storage round-trip efficiency.
    /// Total energy consumed from storage divided by total energy put into storage, from both grid and on site generation.
    /// (unit: kWh)
    storage_efficiency: f64,
    /// Net import grid to consumption minus generation to grid (unit: kWh)
    net_import: f64,
    /// Gross import from grid (unit: kWh)
    total_gross_import: f64,
    /// Gross export to grid (unit: kWh)
    total_gross_export: f64,
}

/// Summary outputs from HEM Core.
#[derive(Debug, Serialize)]
struct OutputSummary {
    /// Total floor-area of all zones (unit: m²)
    total_floor_area: f64,
    /// Space heating demand total (unit: kWh)
    space_heat_demand_total: f64,
    /// Space cooling demand total (unit: kWh)
    space_cool_demand_total: f64,
    electricity_peak_consumption: OutputSummaryPeakElectricityConsumption,
    energy_supply: IndexMap<String, OutputSummaryEnergySupply>,
    /// Delivered energy summary, total energy per fuel and end-use (unit: kWh)
    delivered_energy: IndexMap<String, IndexMap<String, f64>>,
    /// 75th percentile of hot water demand summed over each 24 hour segment of the simulation.
    hot_water_demand_daily_75th_percentile: IndexMap<String, f64>,
}

impl OutputSummary {
    /// Space heating demand total per floor-area (unit: kWh/m²)
    fn space_heat_demand_by_floor_area(&self) -> f64 {
        self.space_heat_demand_total / self.total_floor_area
    }

    /// Space cooling demand total per floor-area (unit: kWh/m²)
    fn space_cool_demand_by_floor_area(&self) -> f64 {
        self.space_cool_demand_total / self.total_floor_area
    }

    fn delivered_energy_by_floor_area(&self) -> IndexMap<String, IndexMap<String, f64>> {
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
struct OutputMetadata {
    hem_core_version: String,
}

#[derive(Debug, Serialize)]
struct Output {
    static_: OutputStatic,
    core: OutputCore,
    summary: OutputSummary,
    metadata: OutputMetadata,
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
            IndexMap<String, IndexMap<String, f64>>,
            IndexMap<String, IndexMap<String, f64>>,
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
