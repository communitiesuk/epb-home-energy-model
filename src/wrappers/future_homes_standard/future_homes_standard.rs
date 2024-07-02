use crate::corpus::{KeyString, ResultsEndUser};
use crate::input::{
    EnergySupplyDetails, FuelType, HeatSourceControlType, HeatingControlType,
    HotWaterSourceDetailsForProcessing, Input, InputForProcessing, SpaceHeatControlType,
    WaterHeatingEvent, WaterHeatingEventType,
};
use crate::output::Output;
use crate::simulation_time::SimulationTime;
use crate::wrappers::future_homes_standard::fhs_hw_events::{
    reset_events_and_provide_drawoff_generator, HotWaterEventGenerator,
};
use anyhow::{anyhow, bail};
use arrayvec::{ArrayString, ArrayVec};
use csv::{Reader, WriterBuilder};
use indexmap::IndexMap;
use lazy_static::lazy_static;
use log::warn;
use serde::Deserialize;
use serde_json::{json, Number, Value};
use std::collections::{HashMap, HashSet};
use std::io::{BufReader, Cursor};
use std::iter::{repeat, zip};

const _EMIS_FACTOR_NAME: &str = "Emissions Factor kgCO2e/kWh";
const _EMIS_OOS_FACTOR_NAME: &str = "Emissions Factor kgCO2e/kWh including out-of-scope emissions";
const _PE_FACTOR_NAME: &str = "Primary Energy Factor kWh/kWh delivered";

pub const ENERGY_SUPPLY_NAME_ELECTRICITY: &str = "mains elec";
const APPL_OBJ_NAME: &str = "appliances";
const ELEC_COOK_OBJ_NAME: &str = "Eleccooking";
const GAS_COOK_OBJ_NAME: &str = "Gascooking";
const HW_TIMER_MAIN_NAME: &str = "hw timer";
const HW_TIMER_HOLD_AT_SETPNT_NAME: &str = "hw timer eco7";

const LIVING_ROOM_SETPOINT_FHS: f64 = 21.0;
const REST_OF_DWELLING_SETPOINT_FHS: f64 = 20.0;

const SIMTIME_START: f64 = 0.;
const SIMTIME_END: f64 = 8760.;
const SIMTIME_STEP: f64 = 0.5;
fn simtime() -> SimulationTime {
    SimulationTime::new(SIMTIME_START, SIMTIME_END, SIMTIME_STEP)
}

pub fn apply_fhs_preprocessing(input: &mut InputForProcessing) -> anyhow::Result<()> {
    input.set_simulation_time(simtime());

    input.reset_internal_gains();

    let tfa = calc_tfa(input);

    let nbeds = calc_nbeds(input)?;

    let n_occupants = calc_n_occupants(tfa, nbeds)?;

    // construct schedules
    let (schedule_occupancy_weekday, schedule_occupancy_weekend) = create_occupancy(n_occupants);

    create_metabolic_gains(
        input,
        tfa,
        schedule_occupancy_weekday,
        schedule_occupancy_weekend,
    )?;
    create_water_heating_pattern(input)?;
    create_heating_pattern(input)?;
    create_evaporative_losses(input, tfa, n_occupants)?;
    create_lighting_gains(input, tfa, n_occupants)?;
    create_cooking_gains(input, tfa, n_occupants)?;
    create_appliance_gains(input, tfa, n_occupants)?;

    for source_key in input.hot_water_source_keys() {
        let source = input.hot_water_source_details_for_key(&source_key);
        if source.is_storage_tank() {
            source.set_min_temp_and_setpoint_temp_if_storage_tank(52.0, 60.0);
        }
    }

    let cold_water_feed_temps = create_cold_water_feed_temps(input)?;
    create_hot_water_use_pattern(input, n_occupants, &cold_water_feed_temps)?;
    create_cooling(input)?;
    create_window_opening_schedule(input)?;

    Ok(())
}

lazy_static! {
    static ref EMIS_PE_FACTORS: HashMap<String, FactorData> = {
        let mut factors: HashMap<String, FactorData> = Default::default();

        let mut factors_reader = Reader::from_reader(BufReader::new(Cursor::new(include_str!(
            "./FHS_emisPEfactors_07-06-2023.csv"
        ))));
        for factor_data in factors_reader.deserialize() {
            let factor_data: FactorData = factor_data.expect("Reading the PE factors file failed.");
            if let Some(fuel_code) = &factor_data.fuel_code {
                factors.insert(fuel_code.clone(), factor_data);
            }
        }

        factors
    };
}

#[derive(Clone, Debug, Deserialize)]
struct FactorData {
    #[serde(rename = "Fuel Code")]
    fuel_code: Option<String>,
    #[serde(rename = "Fuel")]
    _fuel: String,
    #[serde(rename = "Emissions Factor kgCO2e/kWh")]
    emissions_factor: f64,
    #[serde(rename = "Emissions Factor kgCO2e/kWh including out-of-scope emissions")]
    emissions_factor_including_out_of_scope_emissions: f64,
    #[serde(rename = "Primary Energy Factor kWh/kWh delivered")]
    primary_energy_factor: f64,
}

pub fn apply_fhs_postprocessing(
    input: &Input,
    output: &impl Output,
    energy_import: &IndexMap<ArrayString<64>, Vec<f64>>,
    energy_export: &IndexMap<ArrayString<64>, Vec<f64>>,
    results_end_user: &ResultsEndUser,
    timestep_array: &[f64],
    notional: bool,
) -> anyhow::Result<()> {
    let no_of_timesteps = timestep_array.len();

    // Add unmet demand to list of EnergySupply objects

    // For each EnergySupply object:
    // look up relevant factors for import/export from csv or custom factors
    // from input file
    // - look up relevant factors for generation from csv
    // - apply relevant factors for import, export and generation
    // Applying factors in this way rather than applying a net export factor to
    // exported energy accounts for energy generated and used on site and also
    // accounts for battery storage losses
    let mut emis_results: IndexMap<String, FhsCalculationResult> = Default::default();
    let mut emis_oos_results: IndexMap<String, FhsCalculationResult> = Default::default();
    let mut pe_results: IndexMap<String, FhsCalculationResult> = Default::default();

    for (energy_supply_key, energy_supply_details) in input
        .energy_supply
        .iter()
        .map(|(key, value)| (String::from(key), value))
        .chain(
            [(
                "_unmet_demand".to_string(),
                &EnergySupplyDetails {
                    fuel: FuelType::UnmetDemand,
                    diverter: None,
                    electric_battery: None,
                    factor: None,
                },
            )]
            .into_iter(),
        )
    {
        let supply_emis_result = emis_results.entry(energy_supply_key.clone()).or_default();
        let supply_emis_oos_result = emis_oos_results
            .entry(energy_supply_key.clone())
            .or_default();
        let supply_pe_result = pe_results.entry(energy_supply_key.clone()).or_default();

        let fuel_code = energy_supply_details.fuel;

        // Get emissions/PE factors for import/export
        let (emis_factor_import_export, emis_oos_factor_import_export, pe_factor_import_export) =
            if fuel_code == FuelType::Custom {
                let factor = energy_supply_details.factor.expect("Expected custom fuel type to have associated factor values as part of energy supply input.");
                (
                    factor.emissions,
                    factor.emissions_including_out_of_scope,
                    factor.primary_energy_factor,
                )
            } else {
                let factor = EMIS_PE_FACTORS
                    .get(&fuel_code.to_string())
                    .unwrap_or_else(|| {
                        panic!("Expected factor values in the table for the fuel code {fuel_code} were not present.");
                    });
                (
                    factor.emissions_factor,
                    factor.emissions_factor_including_out_of_scope_emissions,
                    factor.primary_energy_factor,
                )
            };

        // Calculate energy imported and associated emissions/PE
        supply_emis_result.import = energy_import[&KeyString::from(&energy_supply_key).unwrap()]
            .iter()
            .map(|x| x * emis_factor_import_export)
            .collect::<Vec<_>>();
        supply_emis_oos_result.import = energy_import
            [&KeyString::from(&energy_supply_key).unwrap()]
            .iter()
            .map(|x| x * emis_oos_factor_import_export)
            .collect::<Vec<_>>();
        supply_pe_result.import = energy_import[&KeyString::from(&energy_supply_key).unwrap()]
            .iter()
            .map(|x| x * pe_factor_import_export)
            .collect::<Vec<_>>();

        // If there is any export, Calculate energy exported and associated emissions/PE
        // Note that by convention, exported energy is negative
        (
            supply_emis_result.export,
            supply_emis_oos_result.export,
            supply_pe_result.export,
        ) = if energy_export[&KeyString::from(&energy_supply_key).unwrap()]
            .iter()
            .sum::<f64>()
            < 0.
        {
            (
                energy_export[&KeyString::from(&energy_supply_key).unwrap()]
                    .iter()
                    .map(|x| x * emis_factor_import_export)
                    .collect::<Vec<_>>(),
                energy_export[&KeyString::from(&energy_supply_key).unwrap()]
                    .iter()
                    .map(|x| x * emis_oos_factor_import_export)
                    .collect::<Vec<_>>(),
                energy_export[&KeyString::from(&energy_supply_key).unwrap()]
                    .iter()
                    .map(|x| x * pe_factor_import_export)
                    .collect::<Vec<_>>(),
            )
        } else {
            (
                vec![0.; no_of_timesteps],
                vec![0.; no_of_timesteps],
                vec![0.; no_of_timesteps],
            )
        };

        // Calculate energy generated and associated emissions/PE
        let mut energy_generated = vec![0.; no_of_timesteps];
        for end_user_energy in
            results_end_user[&KeyString::from(&energy_supply_key).unwrap()].values()
        {
            if end_user_energy.iter().sum::<f64>() < 0. {
                for t_idx in 0..no_of_timesteps {
                    // Subtract here because generation is represented as negative demand
                    *energy_generated.get_mut(t_idx).unwrap() -= end_user_energy[t_idx];
                }
            }
        }

        (
            supply_emis_result.generated,
            supply_emis_oos_result.generated,
            supply_pe_result.generated,
        ) = if energy_generated.iter().sum::<f64>() > 0. {
            // TODO (from Python) Allow custom (user-defined) factors for generated energy?
            let fuel_code_generated = format!("{}_generated", fuel_code);
            let generated_factor = EMIS_PE_FACTORS.get(&fuel_code_generated).unwrap_or_else(|| panic!("Fuel code '{fuel_code}' does not have a generated row in the EMIS factors file."));
            let FactorData {
                emissions_factor: emis_factor_generated,
                emissions_factor_including_out_of_scope_emissions: emis_oos_factor_generated,
                primary_energy_factor: pe_factor_generated,
                ..
            } = generated_factor;

            (
                energy_generated
                    .iter()
                    .map(|x| x * emis_factor_generated)
                    .collect::<Vec<_>>(),
                energy_generated
                    .iter()
                    .map(|x| x * emis_oos_factor_generated)
                    .collect::<Vec<_>>(),
                energy_generated
                    .iter()
                    .map(|x| x * pe_factor_generated)
                    .collect::<Vec<_>>(),
            )
        } else {
            (
                vec![0.; no_of_timesteps],
                vec![0.; no_of_timesteps],
                vec![0.; no_of_timesteps],
            )
        };

        let mut energy_unregulated = vec![0.; no_of_timesteps];
        for (end_user_name, end_user_energy) in
            results_end_user[&KeyString::from(&energy_supply_key).unwrap()].iter()
        {
            if [APPL_OBJ_NAME, ELEC_COOK_OBJ_NAME, GAS_COOK_OBJ_NAME]
                .contains(&end_user_name.as_str())
            {
                for t_idx in 0..no_of_timesteps {
                    *energy_unregulated.get_mut(t_idx).unwrap() += end_user_energy[t_idx];
                }
            }
        }

        supply_emis_result.unregulated = energy_unregulated
            .iter()
            .map(|x| x * emis_factor_import_export)
            .collect::<Vec<_>>();
        supply_emis_oos_result.unregulated = energy_unregulated
            .iter()
            .map(|x| x * emis_oos_factor_import_export)
            .collect::<Vec<_>>();
        supply_pe_result.unregulated = energy_unregulated
            .iter()
            .map(|x| x * pe_factor_import_export)
            .collect::<Vec<_>>();

        // Calculate total CO2/PE for each EnergySupply based on import and export,
        // subtracting unregulated
        supply_emis_result.total = Vec::with_capacity(no_of_timesteps);
        supply_emis_oos_result.total = Vec::with_capacity(no_of_timesteps);
        supply_pe_result.total = Vec::with_capacity(no_of_timesteps);
        for t_idx in 0..no_of_timesteps {
            supply_emis_result.total.push(
                supply_emis_result.import[t_idx]
                    + supply_emis_result.export[t_idx]
                    + supply_emis_result.generated[t_idx]
                    - supply_emis_result.unregulated[t_idx],
            );
            supply_emis_oos_result.total.push(
                supply_emis_oos_result.import[t_idx]
                    + supply_emis_oos_result.export[t_idx]
                    + supply_emis_oos_result.generated[t_idx]
                    - supply_emis_oos_result.unregulated[t_idx],
            );
            supply_pe_result.total.push(
                supply_pe_result.import[t_idx]
                    + supply_pe_result.export[t_idx]
                    + supply_pe_result.generated[t_idx]
                    - supply_pe_result.unregulated[t_idx],
            );
        }
    }

    let tfa = calc_tfa_from_finalised_input(input);
    let total_emissions_rate = emis_results
        .values()
        .map(|emis| emis.total.iter().sum::<f64>())
        .sum::<f64>()
        / tfa;
    let total_pe_rate = pe_results
        .values()
        .map(|pe| pe.total.iter().sum::<f64>())
        .sum::<f64>()
        / tfa;

    // Write results to output files
    write_postproc_file(output, "emissions", emis_results, no_of_timesteps)?;
    write_postproc_file(
        output,
        "emissions_incl_out_of_scope",
        emis_oos_results,
        no_of_timesteps,
    )?;
    write_postproc_file(output, "primary_energy", pe_results, no_of_timesteps)?;
    write_postproc_summary_file(output, total_emissions_rate, total_pe_rate, notional)?;

    Ok(())
}

#[derive(Default)]
struct FhsCalculationResult {
    import: Vec<f64>,
    export: Vec<f64>,
    generated: Vec<f64>,
    unregulated: Vec<f64>,
    total: Vec<f64>,
}

impl FhsCalculationResult {
    fn labels(&self) -> [&'static str; 5] {
        ["import", "export", "generated", "unregulated", "total"]
    }

    fn printable_values_for_index(&self, index: usize) -> [String; 5] {
        [
            self.import[index].to_string(),
            self.export[index].to_string(),
            self.generated[index].to_string(),
            self.unregulated[index].to_string(),
            self.total[index].to_string(),
        ]
    }
}

fn write_postproc_file(
    output: &impl Output,
    file_location: &str,
    results: IndexMap<String, FhsCalculationResult>,
    no_of_timesteps: usize,
) -> anyhow::Result<()> {
    let file_location = format!("postproc_{file_location}");

    let mut row_headers: Vec<String> = Default::default();
    let mut rows_results: Vec<Vec<String>> = Default::default();

    // Loop over each EnergySupply object and add headers and results to rows
    for (energy_supply, energy_supply_results) in &results {
        for result_name in energy_supply_results.labels() {
            // Create header row
            row_headers.push(format!("{energy_supply} {result_name}"));
        }
    }

    // Create results rows
    for t_idx in 0..no_of_timesteps {
        let mut row = vec![];
        for energy_supply_results in results.values() {
            row.push(energy_supply_results.printable_values_for_index(t_idx));
        }
        rows_results.push(row.iter().flatten().cloned().collect());
    }

    let writer = output.writer_for_location_key(&file_location)?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    writer.write_record(row_headers)?;
    for record in rows_results {
        writer.write_record(record)?;
    }

    writer.flush()?;

    Ok(())
}

fn write_postproc_summary_file(
    output: &impl Output,
    total_emissions_rate: f64,
    total_pe_rate: f64,
    notional: bool,
) -> anyhow::Result<()> {
    let (emissions_rate_name, pe_rate_name) = if notional {
        ("TER", "TPER")
    } else {
        ("DER", "DPER")
    };

    let writer = output.writer_for_location_key("postproc_summary")?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    writer.write_record(["", "", "Total"])?;
    writer.write_record([
        emissions_rate_name,
        "kgCO2/m2",
        total_emissions_rate.to_string().as_str(),
    ])?;
    writer.write_record([pe_rate_name, "kWh/m2", total_pe_rate.to_string().as_str()])?;

    writer.flush()?;

    Ok(())
}

pub fn calc_tfa(input: &InputForProcessing) -> f64 {
    input.total_zone_area()
}

fn calc_tfa_from_finalised_input(input: &Input) -> f64 {
    input.zone.values().map(|z| z.area).sum::<f64>()
}

fn calc_nbeds(input: &InputForProcessing) -> anyhow::Result<usize> {
    match input.number_of_bedrooms() {
        Some(bedrooms) => Ok(bedrooms),
        None => bail!("missing NumberOfBedrooms - required for FHS calculation"),
    }
}

fn calc_n_occupants(total_floor_area: f64, number_of_bedrooms: usize) -> anyhow::Result<f64> {
    if total_floor_area <= 0. {
        bail!("Invalid floor area: {total_floor_area}");
    }

    // sigmoid curve is only used for one bedroom occupancy.
    // Therefore, sigmoid parameters only used if there is one bedroom
    Ok(match number_of_bedrooms {
        1 => {
            1. + ONE_BED_SIGMOID_PARAMS.j
                * (1. - (ONE_BED_SIGMOID_PARAMS.k * total_floor_area.powi(2)).exp())
        }
        2 => TWO_BED_OCCUPANCY,
        3 => THREE_BED_OCCUPANCY,
        4 => FOUR_BED_OCCUPANCY,
        n if n >= 5 => FIVE_BED_OCCUPANCY,
        _ => bail!("Invalid number of bedrooms: {number_of_bedrooms}"),
    })
}

struct SigmoidParams {
    j: f64,
    k: f64,
}

const ONE_BED_SIGMOID_PARAMS: SigmoidParams = SigmoidParams {
    j: 0.4373,
    k: -0.001902,
};
const TWO_BED_OCCUPANCY: f64 = 2.2472;
const THREE_BED_OCCUPANCY: f64 = 2.9796;
const FOUR_BED_OCCUPANCY: f64 = 3.3715;
const FIVE_BED_OCCUPANCY: f64 = 3.8997;

fn create_occupancy(n_occupants: f64) -> ([f64; 24], [f64; 24]) {
    let schedule_occupancy_weekday = OCCUPANCY_WEEKDAY_FHS.map(|factor| factor * n_occupants);
    let schedule_occupancy_weekend = OCCUPANCY_WEEKEND_FHS.map(|factor| factor * n_occupants);

    (schedule_occupancy_weekday, schedule_occupancy_weekend)
}

const OCCUPANCY_WEEKDAY_FHS: [f64; 24] = [
    1., 1., 1., 1., 1., 1., 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.5, 0.5, 0.5, 0.8,
    0.8, 1., 1., 1.,
];
const OCCUPANCY_WEEKEND_FHS: [f64; 24] = [
    1., 1., 1., 1., 1., 1., 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
    0.8, 1., 1., 1.,
];

fn create_metabolic_gains(
    input: &mut InputForProcessing,
    total_floor_area: f64,
    schedule_occupancy_weekday: [f64; 24],
    schedule_occupancy_weekend: [f64; 24],
) -> anyhow::Result<(Vec<f64>, Vec<f64>)> {
    // Profile below is in Watts/m^2 body surface area, average adult has 1.8m^2 surface area
    // Nighttime metabolic rate based on figure for sleeping from CIBSE Guide A
    // Daytime metabolic rate based on figures for "seated quiet" from CIBSE Guide A
    let body_area_average = 1.8;
    let night = 41.0;
    let daytime = 58.0;
    // 7 hours using night figure, 17 hours using daytime
    let mut metabolic_gains_fhs = ArrayVec::<f64, 24>::new();
    metabolic_gains_fhs.extend(repeat(night).take(17));
    metabolic_gains_fhs.extend(repeat(daytime).take(7));
    let metabolic_gains_fhs: [f64; 24] = metabolic_gains_fhs.into_inner().unwrap();
    let schedule_metabolic_gains_weekday: Vec<f64> =
        zip(schedule_occupancy_weekday, metabolic_gains_fhs)
            .map(|(occupancy, gains)| occupancy * body_area_average * gains / total_floor_area)
            .collect();
    let schedule_metabolic_gains_weekend: Vec<f64> =
        zip(schedule_occupancy_weekend, metabolic_gains_fhs)
            .map(|(occupancy, gains)| occupancy * body_area_average * gains / total_floor_area)
            .collect();

    input.set_metabolic_gains(
        0,
        1.,
        json!(
            {
                "main": [{"repeat": 53, "value": "week"}],
                "week": [{"repeat": 5, "value": "weekday"}, {"repeat": 2, "value": "weekend"}],
                "weekday": schedule_metabolic_gains_weekday,
                "weekend": schedule_metabolic_gains_weekend,
            }
        ),
    )?;

    Ok((
        schedule_metabolic_gains_weekday,
        schedule_metabolic_gains_weekend,
    ))
}

fn create_heating_pattern(input: &mut InputForProcessing) -> anyhow::Result<()> {
    // 07:00-09:30 and then 16:30-22:00
    let mut heating_fhs_weekday = Vec::with_capacity(48);
    heating_fhs_weekday.extend(repeat(false).take(14));
    heating_fhs_weekday.extend(repeat(true).take(5));
    heating_fhs_weekday.extend(repeat(false).take(14));
    heating_fhs_weekday.extend(repeat(true).take(11));
    heating_fhs_weekday.extend(repeat(false).take(4));
    let heating_fhs_weekday: [bool; 48] = heating_fhs_weekday.try_into().unwrap();

    // Start all-day HW schedule 1 hour before space heating
    let mut _sched_allday_weekday = Vec::with_capacity(48);
    _sched_allday_weekday.extend(repeat(false).take(13));
    _sched_allday_weekday.extend(repeat(true).take(31));
    _sched_allday_weekday.extend(repeat(false).take(4));
    let _sched_allday_weekday: [bool; 48] = _sched_allday_weekday.try_into().unwrap();

    // 07:00-09:30 and then 18:30-22:00
    let mut heating_nonlivingarea_fhs_weekday = Vec::with_capacity(48);
    heating_nonlivingarea_fhs_weekday.extend(repeat(false).take(14));
    heating_nonlivingarea_fhs_weekday.extend(repeat(true).take(5));
    heating_nonlivingarea_fhs_weekday.extend(repeat(false).take(18));
    heating_nonlivingarea_fhs_weekday.extend(repeat(true).take(7));
    heating_nonlivingarea_fhs_weekday.extend(repeat(false).take(4));
    let heating_nonlivingarea_fhs_weekday: [bool; 48] =
        heating_nonlivingarea_fhs_weekday.try_into().unwrap();

    // 08:30 - 22:00
    let mut heating_fhs_weekend = Vec::with_capacity(48);
    heating_fhs_weekend.extend(repeat(false).take(17));
    heating_fhs_weekend.extend(repeat(true).take(27));
    heating_fhs_weekend.extend(repeat(false).take(4));
    let heating_fhs_weekend: [bool; 48] = heating_fhs_weekend.try_into().unwrap();

    // Start all-day HW schedule 1 hour before space heating
    let mut _hw_sched_allday_weekend = Vec::with_capacity(48);
    _hw_sched_allday_weekend.extend(repeat(false).take(15));
    _hw_sched_allday_weekend.extend(repeat(true).take(29));
    _hw_sched_allday_weekend.extend(repeat(false).take(4));
    let _hw_sched_allday_weekend: [bool; 48] = _hw_sched_allday_weekend.try_into().unwrap();

    // if there is no separate time control of the non-living rooms
    // (i.e. control type 3 in SAP 10 terminology),
    // the heating times are necessarily the same as the living room,
    // so the evening heating period would also start at 16:30 on weekdays.
    let control_type = match input.heating_control_type() {
        Some(HeatingControlType::SeparateTimeAndTemperatureControl) => ControlType::Type3,
        Some(HeatingControlType::SeparateTemperatureControl) => ControlType::Type2,
        None => {
            bail!("missing HeatingControlType (SeparateTempControl or SeparateTimeAndTempControl)")
        }
    };

    let living_room_space_heat_system_name = "HeatingPattern_LivingRoom";
    let rest_of_dwelling_space_heat_system_name = "HeatingPattern_RestOfDwelling";

    for zone in input.zone_keys() {
        match (
            input.space_heat_control_for_zone(zone.as_str())?,
            control_type,
        ) {
            (Some(SpaceHeatControlType::LivingRoom), _) => {
                input.set_init_temp_setpoint_for_zone(zone.as_str(), LIVING_ROOM_SETPOINT_FHS)?;
                let mut living_room_control = json!(
                    {
                        "type": "SetpointTimeControl",
                        "start_day": 0,
                        "time_series_step": 0.5,
                        "schedule": {
                            "main": [{"repeat": 53, "value": "week"}],
                            "week": [{"repeat": 5, "value": "weekday"},
                                    {"repeat": 2, "value": "weekend"}],
                            "weekday": heating_fhs_weekday.iter().map(|on| on.then_some(LIVING_ROOM_SETPOINT_FHS)).collect::<Vec<Option<f64>>>(),
                            "weekend": heating_fhs_weekend.iter().map(|on| on.then_some(LIVING_ROOM_SETPOINT_FHS)).collect::<Vec<Option<f64>>>(),
                        }
                    }
                );
                let space_heat_system = input.space_heat_system_for_zone(zone.as_str())?;
                if let Some(space_heat_system) = space_heat_system {
                    input.set_control_string_for_space_heat_system(space_heat_system.as_str(), living_room_space_heat_system_name)?;
                    let control_schedule = living_room_control.as_object_mut().unwrap().get_mut("schedule").unwrap().as_object_mut().unwrap();
                    if let Some(temp_setback) = input.temperature_setback_for_space_heat_system(space_heat_system.as_str())? {
                        control_schedule.insert("setpoint_min".to_string(), temp_setback.into());
                    }
                    if let Some(advanced_start) = input.advanced_start_for_space_heat_system(space_heat_system.as_str())? {
                        control_schedule.insert("advanced_start".to_string(), advanced_start.into());
                    }
                }
                input.add_control(living_room_space_heat_system_name, living_room_control)?;
            }
            (Some(SpaceHeatControlType::RestOfDwelling), control_type) => {
                input.set_init_temp_setpoint_for_zone(zone.as_str(), REST_OF_DWELLING_SETPOINT_FHS)?;
                let mut rest_of_dwelling_control = json!(
                    {
                        "type": "SetpointTimeControl",
                        "start_day": 0,
                        "time_series_step": 0.5,
                        "schedule": {
                            "main": [{"repeat": 53, "value": "week"}],
                            "week": [{"repeat": 5, "value": "weekday"},
                                    {"repeat": 2, "value": "weekend"}],
                            "weekday": match control_type {
                                ControlType::Type2 => heating_fhs_weekend,
                                ControlType::Type3 => heating_nonlivingarea_fhs_weekday,
                            }.iter().map(|on| on.then_some(REST_OF_DWELLING_SETPOINT_FHS)).collect::<Vec<Option<f64>>>(),
                            "weekend": heating_fhs_weekend.iter().map(|on| on.then_some(REST_OF_DWELLING_SETPOINT_FHS)).collect::<Vec<Option<f64>>>(),
                        }
                    }
                );
                let space_heat_system = input.space_heat_system_for_zone(zone.as_str())?;
                if let Some(space_heat_system) = space_heat_system {
                    input.set_control_string_for_space_heat_system(space_heat_system.as_str(), rest_of_dwelling_space_heat_system_name)?;
                    let control_schedule = rest_of_dwelling_control.as_object_mut().unwrap().get_mut("schedule").unwrap().as_object_mut().unwrap();
                    if let Some(temp_setback) = input.temperature_setback_for_space_heat_system(space_heat_system.as_str())? {
                        control_schedule.insert("setpoint_min".to_string(), temp_setback.into());
                    }
                    if let Some(advanced_start) = input.advanced_start_for_space_heat_system(space_heat_system.as_str())? {
                        control_schedule.insert("advanced_start".to_string(), advanced_start.into());
                    }
                }
                input.add_control(rest_of_dwelling_space_heat_system_name, rest_of_dwelling_control)?;
            }
            (None, _) => todo!("condition to deal with zone that doesnt have specified living room/rest of dwelling"),
        }
    }

    Ok(())
}

#[derive(Clone, Copy, Debug)]
enum ControlType {
    Type2,
    Type3,
}

fn create_water_heating_pattern(input: &mut InputForProcessing) -> anyhow::Result<()> {
    // water heating pattern - if system is not instantaneous, hold at setpoint
    // 00:00-07:00 and then reheat as necessary 24/7
    input.add_control(
        HW_TIMER_MAIN_NAME,
        json!({
            "type": "OnOffTimeControl",
            "start_day": 0,
            "time_series_step": 0.5,
            "schedule": {
                "main": [{"value": "day", "repeat": 365}],
                "day": [{"value": true, "repeat": 48}]
            }
        }),
    )?;
    input.add_control(
        HW_TIMER_HOLD_AT_SETPNT_NAME,
        json!({
            "type": "OnOffTimeControl",
            "start_day": 0,
            "time_series_step": 0.5,
            "schedule": {
                "main": [{"value": "day", "repeat": 365}],
                "day": [
                    {"value": true, "repeat": 14},
                    {"value": false, "repeat": 34}
                ]
            }
        }),
    )?;

    for hwsource in input.hot_water_source_keys() {
        let source = input.hot_water_source_details_for_key(hwsource.as_str());
        if source.is_storage_tank() {
            source.set_control_hold_at_setpoint(HW_TIMER_HOLD_AT_SETPNT_NAME);
            source.set_control_name_for_heat_sources(HW_TIMER_MAIN_NAME)?;
        } else if source.is_combi_boiler() || source.is_point_of_use() || source.is_hiu() {
            // do nothing
        } else {
            bail!("Standard water heating schedule not defined for HotWaterSource type")
        }
    }

    Ok(())
}

fn create_evaporative_losses(
    input: &mut InputForProcessing,
    total_floor_area: f64,
    number_of_occupants: f64,
) -> anyhow::Result<()> {
    let evaporative_losses_fhs = -40. * number_of_occupants / total_floor_area;

    input.set_evaporative_losses(
        0,
        1.,
        json!({
            "main": [{"value": evaporative_losses_fhs, "repeat": 8760}]
        }),
    )?;

    Ok(())
}

/// Calculate the annual energy requirement in kWh using the procedure described in SAP 10.2 up to and including step 9.
/// Divide this by 365 to get the average daily energy use.
/// Multiply the daily energy consumption figure by the following profiles to
/// create a daily profile for each month of the year (to be applied to all days in that month).
fn create_lighting_gains(
    input: &mut InputForProcessing,
    total_floor_area: f64,
    number_of_occupants: f64,
) -> anyhow::Result<()> {
    // here we calculate an overall lighting efficacy as
    // the average of zone lighting efficacies weighted by zone
    // floor area.
    let mut lighting_efficacy = 0.;
    for zone_key in input.zone_keys() {
        let zone_lighting_efficacy =
            input
                .lighting_efficacy_for_zone(zone_key.as_str())?
                .ok_or(anyhow!(
                    "Lighting efficacy for zone {zone_key} not provided"
                ))?;
        lighting_efficacy += zone_lighting_efficacy * input.area_for_zone(zone_key.as_str())?;
    }
    if lighting_efficacy == 0. {
        bail!("invalid/missing lighting efficacy for all zones");
    }

    // fron analysis of EFUS 2017 data
    let lumens = 1_418. * (total_floor_area * number_of_occupants).powf(0.41);

    // dropped 1/3 - 2/3 split based on SAP2012 assumptions about portable lighting
    let kwh_per_year = lumens / lighting_efficacy;
    let kwh_per_day = kwh_per_year / 365.;

    // To obtain the lighting gains,
    // the above should be converted to Watts by multiplying the individual half-hourly figure by (2 x 1000).
    // Since some lighting energy will be used in external light
    // (e.g. outdoor security lights or lights in unheated spaces like garages and sheds)
    // a factor of 0.85 is also applied to get the internal gains from lighting.
    let [lighting_gains_w_jan, lighting_gains_w_feb, lighting_gains_w_mar, lighting_gains_w_apr, lighting_gains_w_may, lighting_gains_w_jun, lighting_gains_w_jul, lighting_gains_w_aug, lighting_gains_w_sep, lighting_gains_w_oct, lighting_gains_w_nov, lighting_gains_w_dec] =
        AVERAGE_MONTHLY_LIGHTING_HALF_HOUR_PROFILES
            .map(|monthly_profile| monthly_profile.map(|frac| (frac * kwh_per_day) * 2. * 1_000.));

    input.set_lighting_gains(json!({
        "type": "lighting",
        "start_day": 0,
        "time_series_step": 0.5,
        "gains_fraction": 0.85,
        "EnergySupply": ENERGY_SUPPLY_NAME_ELECTRICITY,
        "schedule": {
            "main": [
                {"value": "jan", "repeat": 31},
                {"value": "feb", "repeat": 28},
                {"value": "mar", "repeat": 31},
                {"value": "apr", "repeat": 30},
                {"value": "may", "repeat": 31},
                {"value": "jun", "repeat": 30},
                {"value": "jul", "repeat": 31},
                {"value": "aug", "repeat": 31},
                {"value": "sep", "repeat": 30},
                {"value": "oct", "repeat": 31},
                {"value": "nov", "repeat": 30},
                {"value": "dec", "repeat": 31},
            ],
            "jan": lighting_gains_w_jan.to_vec(),
            "feb": lighting_gains_w_feb.to_vec(),
            "mar": lighting_gains_w_mar.to_vec(),
            "apr": lighting_gains_w_apr.to_vec(),
            "may": lighting_gains_w_may.to_vec(),
            "jun": lighting_gains_w_jun.to_vec(),
            "jul": lighting_gains_w_jul.to_vec(),
            "aug": lighting_gains_w_aug.to_vec(),
            "sep": lighting_gains_w_sep.to_vec(),
            "oct": lighting_gains_w_oct.to_vec(),
            "nov": lighting_gains_w_nov.to_vec(),
            "dec": lighting_gains_w_dec.to_vec(),
        }
    }))?;

    Ok(())
}

fn create_cooking_gains(
    input: &mut InputForProcessing,
    _total_floor_area: f64,
    number_of_occupants: f64,
) -> anyhow::Result<()> {
    // check for gas and/or electric cooking. Remove any existing objects
    // so that we can add our own (just one for gas and one for elec)
    let cooking_fields = input.appliance_gains_fields_for_cooking();
    let mut cooking_energy_supplies = HashSet::new();
    for field in cooking_fields.iter() {
        if let Some(energy_supply_type) = input.energy_supply_type_for_appliance_gains_field(field)
        {
            cooking_energy_supplies.insert(energy_supply_type);
        }
        input.reset_appliance_gains_field(field)?;
    }

    // from the cooking energy supplies, need to find the associated fuel they use
    let mut cooking_fuels = HashSet::new();
    for supply in cooking_energy_supplies.iter() {
        if let Ok(fuel_type) = input.fuel_type_for_energy_supply_field(supply) {
            cooking_fuels.insert(fuel_type);
        }
    }

    let (ec1_elec, ec2_elec, ec1_gas, ec2_gas) = match (
        cooking_fuels.contains("electricity"),
        cooking_fuels.contains("mains_gas"),
    ) {
        (true, true) => (86, 49, 150, 86),
        (_, true) => (0, 0, 299, 171),
        (true, _) => (171, 98, 0, 0),
        _ => (0, 0, 0, 0),
    };

    let annual_cooking_elec_kwh = ec1_elec as f64 + ec2_elec as f64 * number_of_occupants;
    let annual_cooking_gas_kwh = ec1_gas as f64 + ec2_gas as f64 * number_of_occupants;

    // energy consumption, W_m2, gains factor not applied
    let cooking_elec_profile_w = COOKING_PROFILE_FHS
        .map(|half_hour| (1_000 * 2) as f64 * annual_cooking_elec_kwh / 365. * half_hour);
    let cooking_gas_profile_w = COOKING_PROFILE_FHS
        .map(|half_hour| (1_000 * 2) as f64 * annual_cooking_gas_kwh / 365. * half_hour);

    // add back gas and electric cooking gains if they are present
    if cooking_fuels.contains("mains gas") {
        input.set_gains_for_field(
            GAS_COOK_OBJ_NAME,
            json!({
                "type": "cooking",
                "EnergySupply": "mains gas",
                "start_day": 0,
                "time_series_step": 0.5,
                "gains_fraction": 0.5,
                "schedule": {
                    "main": [{"repeat": 365, "value": "day"}],
                    "day": cooking_gas_profile_w.to_vec()
                }
            }),
        )?;
    }
    if cooking_fuels.contains("electricity") {
        input.set_gains_for_field(
            ELEC_COOK_OBJ_NAME,
            json!({
                "type": "cooking",
                "EnergySupply": "mains gas",
                "start_day": 0,
                "time_series_step": 0.5,
                "gains_fraction": 0.5,
                "schedule": {
                    "main": [{"repeat": 365, "value": "day"}],
                    "day": cooking_elec_profile_w.to_vec()
                }
            }),
        )?;
    }

    Ok(())
}

fn create_appliance_gains(
    input: &mut InputForProcessing,
    total_floor_area: f64,
    number_of_occupants: f64,
) -> anyhow::Result<()> {
    // old relation based on sap2012, efus 1998 data verified in 2013
    // EA_annual_kWh = 207.8 * (TFA * N_occupants) ** 0.4714

    // new relation based on analysis of EFUS 2017 monitoring data
    let ea_annual_kwh = 145. * (total_floor_area * number_of_occupants).powf(0.49);

    let appliance_gains_w = AVERAGE_MONTHLY_APPLIANCES_HALF_HOUR_PROFILES
        .map(|month| month.map(|frac| 1_000. * ea_annual_kwh * frac / 365.));

    input.set_gains_for_field(
        APPL_OBJ_NAME,
        json!({
            "type": "appliances",
            "EnergySupply": ENERGY_SUPPLY_NAME_ELECTRICITY,
            "start_day": 0,
            "time_series_step": 1,
            // Internal gains are reduced from washer/dryers and dishwasher waste heat losses.
            // Assume 70% of their heat is lost as waste heat in waste water or vented hot air,
            // or 30% of total appliance energy, leaving 70% appliance gains fraction
            "gains_fraction": 0.7,
            "schedule": {
                // watts
                "main": [
                    {"value": "jan", "repeat": 31},
                    {"value": "feb", "repeat": 28},
                    {"value": "mar", "repeat": 31},
                    {"value": "apr", "repeat": 30},
                    {"value": "may", "repeat": 31},
                    {"value": "jun", "repeat": 30},
                    {"value": "jul", "repeat": 31},
                    {"value": "aug", "repeat": 31},
                    {"value": "sep", "repeat": 30},
                    {"value": "oct", "repeat": 31},
                    {"value": "nov", "repeat": 30},
                    {"value": "dec", "repeat": 31},
                ],
                "jan": appliance_gains_w[0].to_vec(),
                "feb": appliance_gains_w[1].to_vec(),
                "mar": appliance_gains_w[2].to_vec(),
                "apr": appliance_gains_w[3].to_vec(),
                "may": appliance_gains_w[4].to_vec(),
                "jun": appliance_gains_w[5].to_vec(),
                "jul": appliance_gains_w[6].to_vec(),
                "aug": appliance_gains_w[7].to_vec(),
                "sep": appliance_gains_w[8].to_vec(),
                "oct": appliance_gains_w[9].to_vec(),
                "nov": appliance_gains_w[10].to_vec(),
                "dec": appliance_gains_w[11].to_vec()
            }
        }),
    )?;

    Ok(())
}

/// Check (almost an assert) whether the shower flow rate is not less than the minimum allowed.
fn check_shower_flowrate(input: &InputForProcessing) -> anyhow::Result<()> {
    let min_flowrate = 8.0;

    if let Some(flowrate) = input.shower_flowrate() {
        if flowrate < min_flowrate {
            // only currently known shower name that can have a flowrate is "mixer"
            bail!("Invalid flow rate: {flowrate} l/s in shower with name 'mixer'");
        }
    }

    Ok(())
}

fn create_hot_water_use_pattern(
    input: &mut InputForProcessing,
    number_of_occupants: f64,
    cold_water_feed_temps: &[f64],
) -> anyhow::Result<()> {
    check_shower_flowrate(input)?;

    // temperature of mixed hot water for event
    let event_temperature = 41.0;
    let hw_temperature = 52.0;
    let mean_feedtemp =
        cold_water_feed_temps.iter().sum::<f64>() / cold_water_feed_temps.len() as f64;
    let _mean_delta_t = hw_temperature - mean_feedtemp;

    let _annual_hw_events: Vec<()> = vec![];
    let _annual_hw_events_energy: Vec<()> = vec![];
    let startmod = 0;

    // SAP 2012 relation
    // vol_daily_average = (25 * N_occupants) + 36

    // new relation based on Boiler Manufacturer data and EST surveys
    // reduced by 15% to account for pipework losses present in the source data
    let vol_hw_daily_average = 0.85 * 60.3 * number_of_occupants.powf(0.71);

    let mut hw_event_gen = HotWaterEventGenerator::new(vol_hw_daily_average, None, None)?;
    let ref_event_list = hw_event_gen.build_annual_hw_events(startmod)?;
    let mut ref_hw_vol = 0.;

    for event in &ref_event_list {
        // NB while calibration is done by event volumes we use the event durations from the HW csv data for showers
        // so the actual hw use predicted by sap depends on shower flowrates in dwelling, but this value does not
        ref_hw_vol += event.volume;
    }

    // Add daily average hot water use to hot water only heat pump (HWOHP) object, if present
    // TODO (from Python) This is probably only valid if HWOHP is the only heat source for the
    // storage tank. Make this more robust/flexible in future.
    input.override_vol_hw_daily_average_on_heat_pumps(vol_hw_daily_average);

    let fhw = (365. * vol_hw_daily_average) / ref_hw_vol;

    // if part G has been complied with, apply 5% reduction to duration of Other events
    let part_g_bonus = if let Some(part_g_compliance) = input.part_g_compliance() {
        if part_g_compliance {
            0.95
        } else {
            1.0
        }
    } else {
        bail!("Part G compliance missing from input file");
    };

    let mut hw_event_aa = reset_events_and_provide_drawoff_generator(
        input,
        fhw,
        event_temperature,
        hw_temperature,
        cold_water_feed_temps,
        part_g_bonus,
    )?;

    // now create lists of events
    // Shower events should be evenly spread across all showers in dwelling
    // and so on for baths etc
    let mut hourly_events: Vec<Vec<HourlyHotWaterEvent>> =
        std::iter::repeat_with(Vec::new).take(8760).collect();
    for event in &ref_event_list {
        let drawoff = if event.event_type.is_shower_type() {
            hw_event_aa.get_shower()
        } else if event.event_type.is_bath_type() {
            hw_event_aa.get_bath()
        } else {
            hw_event_aa.get_other()
        };

        let event_start = event.time;
        let duration = drawoff.call_duration_fn(*event);
        if !input.shower_name_refers_to_instant_electric(&drawoff.name) {
            // IES can overlap with anything so ignore them entirely
            // TODO (from Python) - implies 2 uses of the same IES may overlap, could check them separately
            hw_event_gen.overlap_check(
                &mut hourly_events,
                &[WaterHeatingEventType::Bath, WaterHeatingEventType::Shower],
                event_start,
                duration,
            );
            hourly_events
                .get_mut(event_start.floor() as usize)
                .unwrap()
                .push(HourlyHotWaterEvent {
                    event_type: WaterHeatingEventType::Shower,
                    start: event_start,
                    end: event_start + duration / 60.,
                });
        }

        input.add_water_heating_event(
            drawoff.event_type,
            &drawoff.name,
            WaterHeatingEvent {
                start: event_start,
                duration: Some(duration),
                temperature: event_temperature,
            },
        );
    }

    Ok(())
}

fn create_window_opening_schedule(input: &mut InputForProcessing) -> anyhow::Result<()> {
    if !input.defines_window_opening_for_cooling() {
        warn!("Warning: No window opening for cooling has been specified. The calculation will assume that there are no openable windows.");
        return Ok(());
    }

    let window_opening_setpoint = 22.0;

    // 09:00 - 22:00
    input.add_control(
        "WindowOpening_LivingRoom",
        json!({
            "type": "SetpointTimeControl",
            "start_day": 0,
            "time_series_step": 0.5,
            "schedule": {
                "main": [{"repeat": 365, "value": "day"}],
                "day": [
                    {"repeat": 18, "value": Value::Null},
                    {"repeat": 26, "value": window_opening_setpoint},
                    {"repeat": 4, "value": Value::Null},
                ]
            }
        }),
    )?;

    // 08:00 - 23:00
    input.add_control(
        "WindowOpening_RestOfDwelling",
        json!({
            "type": "SetpointTimeControl",
            "start_day": 0,
            "time_series_step": 0.5,
            "schedule": {
                "main": [{"repeat": 365, "value": "day"}],
                "day": [
                    {"repeat": 16, "value": Value::Null},
                    {"repeat": 30, "value": window_opening_setpoint},
                    {"repeat": 2, "value": Value::Null},
                ]
            }
        }),
    )?;

    let zone_keys = input.zone_keys();
    for zone_key in zone_keys {
        match input.space_heat_control_for_zone(&zone_key)? {
            Some(SpaceHeatControlType::LivingRoom) => {
                input.set_control_window_opening_for_zone(
                    &zone_key,
                    Some(HeatSourceControlType::WindowOpeningLivingRoom),
                )?;
            }
            Some(SpaceHeatControlType::RestOfDwelling) => {
                input.set_control_window_opening_for_zone(
                    &zone_key,
                    Some(HeatSourceControlType::WindowOpeningRestOfDwelling),
                )?;
            }
            None => {
                bail!("Space heat control for zone '{zone_key}' was not of expected type.");
            }
        }
    }

    Ok(())
}

fn create_cooling(input: &mut InputForProcessing) -> anyhow::Result<()> {
    let zone_keys = input.zone_keys();
    for zone_key in &zone_keys {
        if let Some(space_heat_control) = input.space_heat_control_for_zone(zone_key)? {
            match space_heat_control {
                SpaceHeatControlType::LivingRoom => {
                    if let Some(space_cool_system) = input.space_cool_system_for_zone(zone_key)? {
                        let mut living_room_control = json!({
                            "type": "SetpointTimeControl",
                            "start_day": 0,
                            "time_series_step": 0.5,
                            "schedule": {
                                "main": [{"repeat": 53, "value": "week"}],
                                "week": [{"repeat": 5, "value": "weekday"},
                                        {"repeat": 2, "value": "weekend"}],
                                "weekday": COOLING_SUBSCHEDULE_LIVINGROOM_WEEKDAY.to_vec(),
                                "weekend": COOLING_SUBSCHEDULE_LIVINGROOM_WEEKEND.to_vec(),
                            }
                        });
                        input.set_control_string_for_space_cool_system(
                            &space_cool_system,
                            "Cooling_LivingRoom",
                        )?;
                        if let Some(temp_setback) =
                            input.temperature_setback_for_space_cool_system(&space_cool_system)?
                        {
                            match living_room_control {
                                Value::Object(ref mut control_map) => {
                                    control_map.insert(
                                        "setpoint_max".into(),
                                        Value::Number(Number::from_f64(temp_setback).unwrap()),
                                    );
                                }
                                _ => unreachable!(),
                            }
                        }
                        input.add_control("Cooling_LivingRoom", living_room_control)?;
                    }
                }
                SpaceHeatControlType::RestOfDwelling => {
                    if let Some(space_cool_system) = input.space_cool_system_for_zone(zone_key)? {
                        let mut rest_of_dwelling_control = json!({
                            "type": "SetpointTimeControl",
                            "start_day": 0,
                            "time_series_step": 0.5,
                            "schedule": {
                                "main": [{"repeat": 365, "value": "day"}],
                                "day": COOLING_SUBSCHEDULE_RESTOFDWELLING.to_vec(),
                            }
                        });
                        input.set_control_string_for_space_cool_system(
                            &space_cool_system,
                            "Cooling_RestOfDwelling",
                        )?;
                        if let Some(temp_setback) =
                            input.temperature_setback_for_space_cool_system(&space_cool_system)?
                        {
                            match rest_of_dwelling_control {
                                Value::Object(ref mut control_map) => {
                                    control_map.insert(
                                        "setpoint_max".into(),
                                        Value::Number(Number::from_f64(temp_setback).unwrap()),
                                    );
                                }
                                _ => unreachable!(),
                            }
                        }
                        input.add_control("Cooling_RestOfDwelling", rest_of_dwelling_control)?;
                    }
                }
            }
        }
    }

    Ok(())
}

const COOLING_SETPOINT: f64 = 24.0;

// 07:00-09:30 and then 18:30-22:00
const COOLING_SUBSCHEDULE_LIVINGROOM_WEEKDAY: [Option<f64>; 48] = [
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
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
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
    None,
    None,
    None,
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    None,
    None,
    None,
    None,
];

// 08:30-22:30
const COOLING_SUBSCHEDULE_LIVINGROOM_WEEKEND: [Option<f64>; 48] = [
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
    None,
    None,
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    None,
    None,
    None,
];

// 22:00-07:00 - i.e. nighttime only
const COOLING_SUBSCHEDULE_RESTOFDWELLING: [Option<f64>; 48] = [
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
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
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
    Some(COOLING_SETPOINT),
];

fn create_cold_water_feed_temps(input: &mut InputForProcessing) -> anyhow::Result<Vec<f64>> {
    // 24-hour average feed temperature (degrees Celsius) per month m. SAP 10.2 Table J1
    let t24m_header_tank = [
        11.1, 11.3, 12.3, 14.5, 16.2, 18.8, 21.3, 19.3, 18.7, 16.2, 13.2, 11.2,
    ];
    let t24m_mains = [
        8.0, 8.2, 9.3, 12.7, 14.6, 16.7, 18.4, 17.6, 16.6, 14.3, 11.1, 8.5,
    ];
    // typical fall in feed temp from midnight to 6am
    let delta = 1.5;

    let (t24m, feed_type) = if input.cold_water_source_has_header_tank() {
        (t24m_header_tank, "header tank")
    } else {
        (t24m_mains, "mains water")
    };

    let mut cold_feed_schedule_m: Vec<Vec<f64>> = Vec::with_capacity(12 * 24);

    for t in t24m {
        // typical cold feed temp between 3pm and midnight
        let t_evening_m = t + (delta * 15. / 48.);

        // variation throughout the day
        cold_feed_schedule_m.push(
            (0..6)
                .map(|t| t_evening_m - delta * t as f64 / 6.)
                .chain((6..15).map(|t| t_evening_m - (15 - t) as f64 * delta / 9.))
                .chain((15..24).map(|_| t_evening_m))
                .collect(),
        );
    }

    let output_feed_temp = repeat(&cold_feed_schedule_m[0])
        .take(31)
        .flatten()
        .chain(repeat(&cold_feed_schedule_m[1]).take(28).flatten())
        .chain(repeat(&cold_feed_schedule_m[2]).take(31).flatten())
        .chain(repeat(&cold_feed_schedule_m[3]).take(30).flatten())
        .chain(repeat(&cold_feed_schedule_m[4]).take(31).flatten())
        .chain(repeat(&cold_feed_schedule_m[5]).take(30).flatten())
        .chain(repeat(&cold_feed_schedule_m[6]).take(31).flatten())
        .chain(repeat(&cold_feed_schedule_m[7]).take(31).flatten())
        .chain(repeat(&cold_feed_schedule_m[8]).take(30).flatten())
        .chain(repeat(&cold_feed_schedule_m[9]).take(31).flatten())
        .chain(repeat(&cold_feed_schedule_m[10]).take(30).flatten())
        .chain(repeat(&cold_feed_schedule_m[11]).take(31).flatten())
        .cloned()
        .collect::<Vec<_>>();

    input.set_cold_water_source_by_key(
        feed_type,
        json!({
            "start_day": 0,
            "time_series_step": 1.,
            "temperatures": output_feed_temp.to_vec(),
        }),
    )?;

    Ok(output_feed_temp)
}

#[derive(Clone, Copy)]
pub struct HourlyHotWaterEvent {
    pub event_type: WaterHeatingEventType,
    pub start: f64,
    pub end: f64,
}

const AVERAGE_MONTHLY_LIGHTING_HALF_HOUR_PROFILES: [[f64; 48]; 12] = [
    [
        0.029235831,
        0.02170637,
        0.016683155,
        0.013732757,
        0.011874713,
        0.010023118,
        0.008837131,
        0.007993816,
        0.007544302,
        0.007057335,
        0.007305208,
        0.007595198,
        0.009170401,
        0.013592425,
        0.024221707,
        0.034538234,
        0.035759809,
        0.02561524,
        0.019538678,
        0.017856399,
        0.016146846,
        0.014341097,
        0.013408345,
        0.013240894,
        0.013252628,
        0.013314013,
        0.013417126,
        0.01429735,
        0.014254224,
        0.014902582,
        0.017289786,
        0.023494947,
        0.035462982,
        0.050550653,
        0.065124006,
        0.072629223,
        0.073631053,
        0.074451912,
        0.074003097,
        0.073190397,
        0.071169797,
        0.069983033,
        0.06890179,
        0.066130187,
        0.062654436,
        0.056634675,
        0.047539646,
        0.037801233,
    ],
    [
        0.026270349,
        0.01864863,
        0.014605535,
        0.01133541,
        0.009557625,
        0.008620514,
        0.007385915,
        0.00674999,
        0.006144089,
        0.005812534,
        0.005834644,
        0.006389013,
        0.007680219,
        0.013106226,
        0.021999709,
        0.027144574,
        0.02507541,
        0.0179487,
        0.014855879,
        0.012930469,
        0.011690622,
        0.010230198,
        0.00994897,
        0.009668602,
        0.00969183,
        0.010174279,
        0.011264866,
        0.011500069,
        0.011588248,
        0.011285427,
        0.012248949,
        0.014420402,
        0.01932017,
        0.027098032,
        0.044955369,
        0.062118024,
        0.072183735,
        0.075100799,
        0.075170654,
        0.072433133,
        0.070588417,
        0.069756433,
        0.068356831,
        0.06656098,
        0.06324827,
        0.055573729,
        0.045490296,
        0.035742204,
    ],
    [
        0.02538112,
        0.018177936,
        0.012838313,
        0.00961673,
        0.007914015,
        0.006844738,
        0.00611386,
        0.005458354,
        0.00508359,
        0.004864933,
        0.004817922,
        0.005375289,
        0.006804643,
        0.009702514,
        0.013148583,
        0.013569968,
        0.01293754,
        0.009183378,
        0.007893734,
        0.00666975,
        0.006673791,
        0.006235776,
        0.006096299,
        0.006250229,
        0.006018285,
        0.00670324,
        0.006705105,
        0.006701531,
        0.006893458,
        0.006440525,
        0.006447363,
        0.007359989,
        0.009510975,
        0.011406472,
        0.017428875,
        0.026635564,
        0.042951415,
        0.057993474,
        0.066065305,
        0.067668248,
        0.067593187,
        0.067506237,
        0.065543759,
        0.063020652,
        0.06004127,
        0.052838397,
        0.043077683,
        0.033689246,
    ],
    [
        0.029044978,
        0.020558675,
        0.014440871,
        0.010798435,
        0.008612364,
        0.007330799,
        0.006848797,
        0.006406058,
        0.00602619,
        0.005718987,
        0.005804901,
        0.006746423,
        0.007160898,
        0.008643678,
        0.010489867,
        0.011675722,
        0.011633729,
        0.008939881,
        0.007346857,
        0.007177037,
        0.007113926,
        0.007536109,
        0.007443049,
        0.006922747,
        0.00685514,
        0.006721853,
        0.006695838,
        0.005746367,
        0.005945173,
        0.005250153,
        0.005665752,
        0.006481695,
        0.006585193,
        0.00751989,
        0.009038481,
        0.009984259,
        0.011695555,
        0.014495872,
        0.018177089,
        0.027110627,
        0.042244993,
        0.056861545,
        0.064008071,
        0.062680016,
        0.060886258,
        0.055751568,
        0.048310205,
        0.038721632,
    ],
    [
        0.023835444,
        0.016876637,
        0.012178456,
        0.009349274,
        0.007659691,
        0.006332517,
        0.005611274,
        0.005650048,
        0.005502101,
        0.005168442,
        0.005128425,
        0.005395259,
        0.004998272,
        0.005229362,
        0.006775116,
        0.007912694,
        0.008514274,
        0.006961449,
        0.00630672,
        0.00620858,
        0.005797218,
        0.005397357,
        0.006006318,
        0.005593869,
        0.005241095,
        0.005212189,
        0.00515531,
        0.004906504,
        0.004757624,
        0.004722969,
        0.004975738,
        0.005211879,
        0.005684004,
        0.006331507,
        0.007031149,
        0.008034144,
        0.008731998,
        0.010738922,
        0.013170262,
        0.016638631,
        0.021708313,
        0.0303703,
        0.043713685,
        0.051876584,
        0.054591464,
        0.05074126,
        0.043109775,
        0.033925231,
    ],
    [
        0.023960632,
        0.016910619,
        0.012253193,
        0.009539031,
        0.007685214,
        0.006311553,
        0.00556675,
        0.005140391,
        0.004604673,
        0.004352551,
        0.004156956,
        0.004098101,
        0.00388452,
        0.00433039,
        0.005658606,
        0.006828804,
        0.007253075,
        0.005872749,
        0.004923197,
        0.004521087,
        0.004454765,
        0.004304616,
        0.004466648,
        0.004178716,
        0.004186183,
        0.003934784,
        0.004014114,
        0.003773073,
        0.003469885,
        0.003708517,
        0.003801095,
        0.004367245,
        0.004558263,
        0.005596378,
        0.005862632,
        0.006068665,
        0.006445161,
        0.007402661,
        0.007880006,
        0.009723385,
        0.012243076,
        0.016280074,
        0.023909324,
        0.03586776,
        0.046595858,
        0.047521241,
        0.041417407,
        0.03322265,
    ],
    [
        0.024387138,
        0.017950032,
        0.01339296,
        0.010486231,
        0.008634325,
        0.00752814,
        0.006562675,
        0.006180296,
        0.00566116,
        0.005092682,
        0.004741384,
        0.004680853,
        0.00479228,
        0.004921812,
        0.005950605,
        0.007010479,
        0.007057257,
        0.005651136,
        0.004813649,
        0.00454666,
        0.004121156,
        0.003793481,
        0.004122788,
        0.004107635,
        0.004363668,
        0.004310674,
        0.004122943,
        0.004014391,
        0.004009496,
        0.003805058,
        0.004133355,
        0.004188447,
        0.005268291,
        0.005964825,
        0.005774607,
        0.006292344,
        0.006813734,
        0.007634982,
        0.008723529,
        0.009855823,
        0.012318322,
        0.017097237,
        0.026780014,
        0.037823534,
        0.046797578,
        0.045940354,
        0.039472789,
        0.033058217,
    ],
    [
        0.023920296,
        0.01690733,
        0.012917415,
        0.010191735,
        0.008787867,
        0.007681138,
        0.006600128,
        0.006043227,
        0.005963814,
        0.005885256,
        0.006164212,
        0.005876554,
        0.005432168,
        0.00580157,
        0.00641092,
        0.007280576,
        0.00811752,
        0.007006283,
        0.006505718,
        0.005917892,
        0.005420978,
        0.005527121,
        0.005317478,
        0.004793601,
        0.004577663,
        0.004958332,
        0.005159584,
        0.004925386,
        0.005192686,
        0.0054453,
        0.005400465,
        0.005331386,
        0.005994507,
        0.006370203,
        0.006800758,
        0.007947816,
        0.009005592,
        0.010608225,
        0.012905449,
        0.015976909,
        0.024610768,
        0.036414926,
        0.04680022,
        0.050678553,
        0.051188831,
        0.046725936,
        0.03998602,
        0.032496965,
    ],
    [
        0.022221313,
        0.016428778,
        0.01266253,
        0.010569518,
        0.008926713,
        0.007929788,
        0.007134802,
        0.006773883,
        0.006485147,
        0.006766094,
        0.007202971,
        0.007480145,
        0.008460127,
        0.011414527,
        0.014342431,
        0.01448993,
        0.012040415,
        0.008520428,
        0.0077578,
        0.006421555,
        0.005889369,
        0.005915144,
        0.006229011,
        0.005425193,
        0.005094464,
        0.005674584,
        0.005898523,
        0.006504338,
        0.005893063,
        0.005967896,
        0.0061056,
        0.006017598,
        0.007500459,
        0.008041236,
        0.0099079,
        0.012297435,
        0.01592606,
        0.021574549,
        0.032780393,
        0.04502082,
        0.054970312,
        0.05930568,
        0.060189471,
        0.057269758,
        0.05486585,
        0.047401041,
        0.038520417,
        0.029925316,
    ],
    [
        0.023567522,
        0.016304584,
        0.012443113,
        0.009961033,
        0.008395854,
        0.007242191,
        0.006314956,
        0.005722235,
        0.005385313,
        0.005197814,
        0.005444756,
        0.0064894,
        0.008409762,
        0.015347201,
        0.025458901,
        0.028619409,
        0.023359044,
        0.014869014,
        0.011900433,
        0.010931316,
        0.010085903,
        0.009253621,
        0.008044246,
        0.007866149,
        0.007665985,
        0.007218414,
        0.00797338,
        0.008005782,
        0.007407311,
        0.008118996,
        0.008648934,
        0.010378068,
        0.013347814,
        0.018541666,
        0.026917161,
        0.035860046,
        0.049702909,
        0.063560224,
        0.069741764,
        0.070609245,
        0.069689625,
        0.069439031,
        0.068785313,
        0.065634051,
        0.062207874,
        0.053986076,
        0.043508937,
        0.033498873,
    ],
    [
        0.025283869,
        0.018061868,
        0.013832406,
        0.01099122,
        0.009057752,
        0.007415348,
        0.006415533,
        0.006118688,
        0.005617255,
        0.005084989,
        0.005552217,
        0.006364787,
        0.00792208,
        0.014440148,
        0.02451,
        0.02993728,
        0.024790064,
        0.016859553,
        0.013140437,
        0.012181571,
        0.010857371,
        0.010621789,
        0.010389982,
        0.010087677,
        0.00981219,
        0.0097001,
        0.01014589,
        0.01052881,
        0.01044948,
        0.011167223,
        0.013610154,
        0.02047533,
        0.035335895,
        0.05409712,
        0.067805633,
        0.074003571,
        0.077948793,
        0.078981046,
        0.077543712,
        0.074620225,
        0.072631194,
        0.070886175,
        0.06972224,
        0.068354439,
        0.063806373,
        0.055709895,
        0.045866391,
        0.035248054,
    ],
    [
        0.030992394,
        0.022532047,
        0.016965296,
        0.013268634,
        0.010662773,
        0.008986943,
        0.007580978,
        0.006707669,
        0.00646337,
        0.006180296,
        0.006229094,
        0.006626391,
        0.00780049,
        0.013149437,
        0.022621172,
        0.033064744,
        0.035953213,
        0.029010413,
        0.023490829,
        0.020477646,
        0.018671663,
        0.017186751,
        0.016526661,
        0.015415424,
        0.014552683,
        0.014347935,
        0.014115058,
        0.013739051,
        0.014944386,
        0.017543021,
        0.021605977,
        0.032100988,
        0.049851633,
        0.063453382,
        0.072579104,
        0.076921792,
        0.079601317,
        0.079548711,
        0.078653413,
        0.076225647,
        0.073936893,
        0.073585752,
        0.071911165,
        0.069220452,
        0.065925982,
        0.059952377,
        0.0510938,
        0.041481111,
    ],
];

const COOKING_PROFILE_FHS: [f64; 48] = [
    0.001192419,
    0.000825857,
    0.000737298,
    0.000569196,
    0.000574409,
    0.000573778,
    0.000578369,
    0.000574619,
    0.000678235,
    0.000540799,
    0.000718043,
    0.002631192,
    0.002439288,
    0.003263445,
    0.003600656,
    0.005743044,
    0.011250675,
    0.015107564,
    0.014475307,
    0.016807917,
    0.018698336,
    0.018887283,
    0.021856976,
    0.047785397,
    0.08045051,
    0.099929701,
    0.042473353,
    0.02361216,
    0.015650513,
    0.014345379,
    0.015951211,
    0.01692045,
    0.037738026,
    0.066195428,
    0.062153502,
    0.073415686,
    0.077486476,
    0.069093846,
    0.046706527,
    0.024924648,
    0.014783978,
    0.009192004,
    0.005617715,
    0.0049381,
    0.003529689,
    0.002365773,
    0.001275927,
    0.001139293,
];

const AVERAGE_MONTHLY_APPLIANCES_HALF_HOUR_PROFILES: [[f64; 24]; 12] = [
    [
        0.025995114,
        0.023395603,
        0.022095847,
        0.020796091,
        0.019496336,
        0.022095847,
        0.02729487,
        0.040292427,
        0.048090962,
        0.049390717,
        0.050690473,
        0.049390717,
        0.053289984,
        0.049390717,
        0.050690473,
        0.053289984,
        0.074086076,
        0.087083633,
        0.08188461,
        0.070186809,
        0.064987786,
        0.057189252,
        0.046791206,
        0.033793649,
    ],
    [
        0.025995114,
        0.023395603,
        0.022095847,
        0.020796091,
        0.019496336,
        0.022095847,
        0.02729487,
        0.032493893,
        0.046791206,
        0.051990229,
        0.049390717,
        0.046791206,
        0.048090962,
        0.046791206,
        0.04549145,
        0.049390717,
        0.062388274,
        0.074086076,
        0.080584854,
        0.067587297,
        0.059788763,
        0.050690473,
        0.044191694,
        0.032493893,
    ],
    [
        0.024695359,
        0.020796091,
        0.020796091,
        0.019496336,
        0.020796091,
        0.022095847,
        0.029894381,
        0.041592183,
        0.04549145,
        0.048090962,
        0.04549145,
        0.04549145,
        0.049390717,
        0.048090962,
        0.048090962,
        0.049390717,
        0.057189252,
        0.070186809,
        0.07278632,
        0.067587297,
        0.061088519,
        0.051990229,
        0.041592183,
        0.029894381,
    ],
    [
        0.022095847,
        0.022095847,
        0.022095847,
        0.022095847,
        0.023395603,
        0.029894381,
        0.038992672,
        0.046791206,
        0.046791206,
        0.044191694,
        0.046791206,
        0.048090962,
        0.044191694,
        0.042891939,
        0.044191694,
        0.051990229,
        0.062388274,
        0.061088519,
        0.058489007,
        0.057189252,
        0.050690473,
        0.041592183,
        0.033793649,
        0.024695359,
    ],
    [
        0.024695359,
        0.022095847,
        0.020796091,
        0.020796091,
        0.023395603,
        0.031194137,
        0.038992672,
        0.044191694,
        0.048090962,
        0.046791206,
        0.044191694,
        0.04549145,
        0.041592183,
        0.037692916,
        0.038992672,
        0.049390717,
        0.05458974,
        0.058489007,
        0.051990229,
        0.055889496,
        0.050690473,
        0.041592183,
        0.031194137,
        0.024695359,
    ],
    [
        0.022095847,
        0.020796091,
        0.020796091,
        0.019496336,
        0.020796091,
        0.024695359,
        0.032493893,
        0.042891939,
        0.044191694,
        0.041592183,
        0.040292427,
        0.042891939,
        0.040292427,
        0.038992672,
        0.040292427,
        0.044191694,
        0.053289984,
        0.057189252,
        0.048090962,
        0.048090962,
        0.04549145,
        0.041592183,
        0.031194137,
        0.024695359,
    ],
    [
        0.022095847,
        0.020796091,
        0.020796091,
        0.019496336,
        0.020796091,
        0.024695359,
        0.032493893,
        0.041592183,
        0.042891939,
        0.042891939,
        0.041592183,
        0.041592183,
        0.040292427,
        0.037692916,
        0.037692916,
        0.044191694,
        0.051990229,
        0.05458974,
        0.046791206,
        0.046791206,
        0.04549145,
        0.042891939,
        0.031194137,
        0.024695359,
    ],
    [
        0.022095847,
        0.020796091,
        0.020796091,
        0.019496336,
        0.020796091,
        0.024695359,
        0.032493893,
        0.044191694,
        0.044191694,
        0.044191694,
        0.044191694,
        0.044191694,
        0.042891939,
        0.040292427,
        0.041592183,
        0.044191694,
        0.051990229,
        0.055889496,
        0.050690473,
        0.051990229,
        0.049390717,
        0.042891939,
        0.031194137,
        0.024695359,
    ],
    [
        0.022095847,
        0.020796091,
        0.020796091,
        0.019496336,
        0.023395603,
        0.029894381,
        0.040292427,
        0.041592183,
        0.044191694,
        0.044191694,
        0.04549145,
        0.044191694,
        0.042891939,
        0.042891939,
        0.042891939,
        0.051990229,
        0.059788763,
        0.064987786,
        0.061088519,
        0.058489007,
        0.051990229,
        0.038992672,
        0.031194137,
        0.023395603,
    ],
    [
        0.022095847,
        0.020796091,
        0.019496336,
        0.022095847,
        0.023395603,
        0.029894381,
        0.040292427,
        0.046791206,
        0.049390717,
        0.04549145,
        0.046791206,
        0.049390717,
        0.04549145,
        0.044191694,
        0.04549145,
        0.053289984,
        0.067587297,
        0.07278632,
        0.066287542,
        0.059788763,
        0.053289984,
        0.042891939,
        0.031194137,
        0.023395603,
    ],
    [
        0.024695359,
        0.022095847,
        0.020796091,
        0.020796091,
        0.020796091,
        0.024695359,
        0.029894381,
        0.042891939,
        0.048090962,
        0.049390717,
        0.04549145,
        0.04549145,
        0.046791206,
        0.046791206,
        0.044191694,
        0.051990229,
        0.064987786,
        0.08188461,
        0.076685587,
        0.067587297,
        0.061088519,
        0.05458974,
        0.04549145,
        0.032493893,
    ],
    [
        0.025995114,
        0.023395603,
        0.022095847,
        0.020796091,
        0.019496336,
        0.022095847,
        0.02729487,
        0.032493893,
        0.048090962,
        0.053289984,
        0.051990229,
        0.05458974,
        0.057189252,
        0.051990229,
        0.055889496,
        0.058489007,
        0.075385832,
        0.083184366,
        0.08188461,
        0.068887053,
        0.062388274,
        0.055889496,
        0.046791206,
        0.033793649,
    ],
];

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use rstest::*;

    #[ignore = "useless test reported up to BRE"]
    #[rstest]
    fn test_check_invalid_shower_flowrate() {
        assert!(!false);
    }

    #[ignore = "useless test reported up to BRE"]
    #[rstest]
    fn test_check_valid_shower_flowrate() {
        assert!(true);
    }

    #[ignore = "useless test reported up to BRE"]
    #[rstest]
    fn test_check_minimum_shower_flowrate() {
        assert!(true);
    }

    #[rstest]
    fn test_calc_1_occupant() {
        // test with on occupant and a range of floor areas
        assert_relative_eq!(
            1.075,
            calc_n_occupants(10., 1).unwrap(),
            max_relative = 1e-2
        );
        assert_relative_eq!(
            1.232,
            calc_n_occupants(20., 1).unwrap(),
            max_relative = 1e-2
        );
        assert_relative_eq!(
            1.433,
            calc_n_occupants(50., 1).unwrap(),
            max_relative = 1e-2
        );
        assert_relative_eq!(
            1.437,
            calc_n_occupants(100., 1).unwrap(),
            max_relative = 1e-2
        );
    }

    #[rstest]
    fn test_calc_n_occupants() {
        assert_eq!(2.2472, calc_n_occupants(100., 2).unwrap(),);
        assert_eq!(2.9796, calc_n_occupants(100., 3).unwrap(),);
        assert_eq!(3.3715, calc_n_occupants(100., 4).unwrap(),);
        assert_eq!(3.8997, calc_n_occupants(100., 5).unwrap(),);
        assert_eq!(3.8997, calc_n_occupants(100., 6).unwrap(),);
    }

    #[rstest]
    fn test_calc_n_occupants_invalid_bedrooms() {
        assert!(calc_n_occupants(100., 0).is_err());
    }

    #[rstest]
    fn test_calc_n_occupants_invalid_floor_area() {
        assert!(calc_n_occupants(0., 1).is_err());
        assert!(calc_n_occupants(-1., 1).is_err());
    }
}
