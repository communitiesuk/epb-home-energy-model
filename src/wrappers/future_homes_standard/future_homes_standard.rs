use crate::core::schedule::{expand_numeric_schedule, reject_nulls};
use crate::core::units::{
    DAYS_IN_MONTH, DAYS_PER_YEAR, HOURS_PER_DAY, LITRES_PER_CUBIC_METRE, MINUTES_PER_HOUR,
    SECONDS_PER_HOUR, WATTS_PER_KILOWATT,
};
use crate::corpus::{Corpus, KeyString, OutputOptions, ResultsEndUser};
use crate::external_conditions::{
    create_external_conditions, ExternalConditions, WindowShadingObject,
};
use crate::input::{
    Appliance, ApplianceEntry, ApplianceKey, ApplianceReference, ColdWaterSourceType,
    ControlDetails, EnergySupplyDetails, EnergySupplyType, FuelType, HeatingControlType,
    HotWaterSourceDetailsForProcessing, Input, InputForProcessing,
    MechanicalVentilationForProcessing, SmartApplianceBattery, SpaceHeatControlType,
    SystemReference, TransparentBuildingElement, VentType, WaterHeatingEvent,
    WaterHeatingEventType, WindowTreatmentType, ZoneLightingBulbs,
};
use crate::output::Output;
use crate::simulation_time::SimulationTime;
use crate::wrappers::future_homes_standard::fhs_appliance::FhsAppliance;
use crate::wrappers::future_homes_standard::fhs_hw_events::{
    reset_events_and_provide_drawoff_generator, HotWaterEventGenerator,
};
use crate::HOURS_TO_END_DEC;
use anyhow::{anyhow, bail};
use arrayvec::ArrayString;
use csv::{Reader, WriterBuilder};
use indexmap::IndexMap;
use itertools::Itertools;
use serde::Deserialize;
use serde_json::{json, Number, Value};
use std::collections::HashMap;
use std::io::{BufReader, Cursor, Read};
use std::iter::repeat;
use std::marker::PhantomData;
use std::sync::LazyLock;

const _EMIS_FACTOR_NAME: &str = "Emissions Factor kgCO2e/kWh";
const _EMIS_OOS_FACTOR_NAME: &str = "Emissions Factor kgCO2e/kWh including out-of-scope emissions";
const _PE_FACTOR_NAME: &str = "Primary Energy Factor kWh/kWh delivered";

pub(crate) const ENERGY_SUPPLY_NAME_GAS: &str = "mains gas";
pub(crate) const ENERGY_SUPPLY_NAME_ELECTRICITY: &str = "mains elec";
const APPL_OBJ_NAME: &str = "appliances";
const ELEC_COOK_OBJ_NAME: &str = "Eleccooking";
const GAS_COOK_OBJ_NAME: &str = "Gascooking";

pub(super) const LIVING_ROOM_SETPOINT_FHS: f64 = 21.0;
pub(super) const REST_OF_DWELLING_SETPOINT_FHS: f64 = 20.0;

pub(crate) const SIMTIME_START: f64 = 0.;
pub(crate) const SIMTIME_END: f64 = 8760.;
pub(crate) const SIMTIME_STEP: f64 = 0.5;
fn simtime() -> SimulationTime {
    SimulationTime::new(SIMTIME_START, SIMTIME_END, SIMTIME_STEP)
}

// Central point for hot water temperature (temp_hot_water) across the code
pub(super) const HW_TEMPERATURE: f64 = 52.0;
const HW_SETPOINT_MAX: f64 = 60.0;

// Occupant sleep+wake hours as per Part O
const OCCUPANT_WAKING_HR: usize = 7;
const OCCUPANT_SLEEPING_HR: usize = 23;
pub(crate) struct SimSettings {
    heat_balance: bool,
    detailed_output_heating_cooling: bool,
    _use_fast_solver: bool,
    tariff_data_filename: Option<String>,
}

pub fn apply_fhs_preprocessing(
    input: &mut InputForProcessing,
    is_fee: Option<bool>,
    sim_settings: Option<SimSettings>,
) -> anyhow::Result<()> {
    let is_fee = is_fee.unwrap_or(false);
    let default_sim_settings = SimSettings {
        heat_balance: false,
        detailed_output_heating_cooling: false,
        _use_fast_solver: false,
        tariff_data_filename: None,
    };

    let sim_settings = sim_settings.unwrap_or(default_sim_settings);

    static APPLIANCE_PROPENSITIES: LazyLock<AppliancePropensities<Normalised>> =
        LazyLock::new(|| {
            load_appliance_propensities(Cursor::new(include_str!("./appliance_propensities.csv")))
                .expect("Could not read and parse appliance_propensities.csv")
        });

    static EVAP_PROFILE_DATA: LazyLock<HalfHourWeeklyProfileData> = LazyLock::new(|| {
        load_evaporative_profile(Cursor::new(include_str!("./evap_loss_profile.csv")))
            .expect("Could not read evap_loss_profile.csv.")
    });

    static COLD_WATER_LOSS_PROFILE_DATA: LazyLock<HalfHourWeeklyProfileData> =
        LazyLock::new(|| {
            load_evaporative_profile(Cursor::new(include_str!("./cold_water_loss_profile.csv")))
                .expect("Could not read cold_water_loss_profile.csv")
        });

    input.set_simulation_time(simtime());

    input.reset_internal_gains();

    let tfa = calc_tfa(input);

    let nbeds = calc_nbeds(input)?;

    let n_occupants = calc_n_occupants(tfa, nbeds)?;

    // construct schedules
    let (_schedule_occupancy_weekday, _schedule_occupancy_weekend) =
        create_occupancy(n_occupants, APPLIANCE_PROPENSITIES.occupied);

    create_metabolic_gains(n_occupants, input)?;
    create_water_heating_pattern(input)?;
    create_heating_pattern(input)?;
    create_evaporative_losses(input, tfa, n_occupants, &EVAP_PROFILE_DATA)?;
    create_cold_water_losses(input, tfa, n_occupants, &COLD_WATER_LOSS_PROFILE_DATA)?;
    create_lighting_gains(input, tfa, n_occupants)?;
    create_appliance_gains(input, tfa, n_occupants, &APPLIANCE_PROPENSITIES)?;

    for source_key in input.hot_water_source_keys() {
        let source = input.hot_water_source_details_for_key(&source_key);
        if source.is_storage_tank() {
            source.set_init_temp_if_storage_tank(HW_SETPOINT_MAX);
        } else {
            source.set_setpoint_temp(HW_TEMPERATURE);
        }
    }

    let cold_water_feed_temps = create_cold_water_feed_temps(input)?;
    create_hot_water_use_pattern(input, n_occupants, &cold_water_feed_temps)?;
    create_cooling(input)?;
    create_window_opening_schedule(input)?;
    create_vent_opening_schedule(input)?;
    window_treatment(input)?;
    if !is_fee {
        calc_sfp_mech_vent(input)?;
    }
    if input.has_mechanical_ventilation() {
        create_mev_pattern(input)?;
    }

    set_temp_internal_static_calcs(input);

    if input.clone().has_control_for_loadshifting() {
        // run project for 24 hours to obtain initial estimate for daily heating demand
        sim_24h(input, sim_settings)?;
    }

    Ok(())
}

pub(super) fn set_temp_internal_static_calcs(input: &mut InputForProcessing) {
    input.set_temp_internal_air_static_calcs(Some(LIVING_ROOM_SETPOINT_FHS));
}

static EMIS_PE_FACTORS: LazyLock<HashMap<String, FactorData>> = LazyLock::new(|| {
    let mut factors: HashMap<String, FactorData> = Default::default();

    let mut factors_reader = Reader::from_reader(BufReader::new(Cursor::new(include_str!(
        "./FHS_emisPEfactors_05-08-2024.csv"
    ))));
    for factor_data in factors_reader.deserialize() {
        let factor_data: FactorData = factor_data.expect("Reading the PE factors file failed.");
        if let Some(fuel_code) = &factor_data.fuel_code {
            factors.insert(fuel_code.clone(), factor_data);
        }
    }

    factors
});

static EMIS_PE_FACTORS_ELEC: LazyLock<HashMap<usize, ElectricityFactorData>> =
    LazyLock::new(|| {
        // Load emissions factors and primary energy factors from data file for electricity
        let mut emis_pe_factors_elec: HashMap<usize, ElectricityFactorData> = Default::default();

        let mut factors_reader = Reader::from_reader(BufReader::new(Cursor::new(include_str!(
            "./DEMO_variable_grid_model.csv"
        ))));

        for factor_data in factors_reader.deserialize() {
            let factor_data: ElectricityFactorData =
                factor_data.expect("Reading the PE factors elec file failed.");
            let timestep = &factor_data.timestep;
            emis_pe_factors_elec.insert(*timestep, factor_data);
        }

        emis_pe_factors_elec
    });

static METABOLIC_GAINS: LazyLock<MetabolicGains> = LazyLock::new(|| {
    let (weekday, weekend) = load_metabolic_gains_profile(Cursor::new(include_str!(
        "./dry_metabolic_gains_profile_Wperm2.csv"
    )))
    .expect("Could not load in metabolic gains file.");
    MetabolicGains { weekday, weekend }
});

struct MetabolicGains {
    weekday: [f64; 48],
    weekend: [f64; 48],
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
#[derive(Clone, Debug, Deserialize)]
struct ElectricityFactorData {
    #[serde(rename = "Timestep")]
    timestep: usize,
    #[serde(rename = "Primary Energy Factor kWh/kWh delivered")]
    primary_energy_factor: f64,
    #[serde(rename = "Emissions Factor kgCO2e/kWh")]
    emissions_factor: f64,
    #[serde(rename = "Emissions Factor kgCO2e/kWh including out-of-scope emissions")]
    emissions_factor_including_out_of_scope_emissions: f64,
}

fn apply_energy_factor_series(energy_data: &[f64], factors: &Vec<f64>) -> anyhow::Result<Vec<f64>> {
    if energy_data.len() != factors.len() {
        bail!("Both energy_data and factors list must be of the same length.");
    }
    Ok(energy_data
        .iter()
        .zip(factors)
        .map(|(energy, factor)| energy * factor)
        .collect_vec())
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

    let FinalRates {
        emission_rate: total_emissions_rate,
        primary_energy_rate: total_pe_rate,
        emissions_results: emis_results,
        emissions_out_of_scope_results: emis_oos_results,
        primary_energy_results: pe_results,
    } = calc_final_rates(
        input,
        energy_import,
        energy_export,
        results_end_user,
        no_of_timesteps,
    )?;

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

pub(super) fn calc_final_rates(
    input: &Input,
    energy_import: &IndexMap<ArrayString<64>, Vec<f64>>,
    energy_export: &IndexMap<ArrayString<64>, Vec<f64>>,
    results_end_user: &ResultsEndUser,
    number_of_timesteps: usize,
) -> anyhow::Result<FinalRates> {
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
                &EnergySupplyDetails::with_fuel(FuelType::UnmetDemand),
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
            match fuel_code {
                FuelType::Custom => {
                    let factor = energy_supply_details.factor.expect("Expected custom fuel type to have associated factor values as part of energy supply input.");
                    (
                        vec![factor.emissions],
                        vec![factor.emissions_including_out_of_scope],
                        vec![factor.primary_energy_factor],
                    )
                }
                FuelType::Electricity => {
                    let emis_factor_import_export = EMIS_PE_FACTORS_ELEC
                        .values()
                        .map(|factor| factor.emissions_factor)
                        .collect_vec();
                    let emis_oos_factor_import_export = EMIS_PE_FACTORS_ELEC
                        .values()
                        .map(|factor| factor.emissions_factor_including_out_of_scope_emissions)
                        .collect_vec();
                    let pe_factor_import_export = EMIS_PE_FACTORS_ELEC
                        .values()
                        .map(|factor| factor.primary_energy_factor)
                        .collect_vec();
                    (
                        emis_factor_import_export,
                        emis_oos_factor_import_export,
                        pe_factor_import_export,
                    )
                }
                _ => {
                    let factor = EMIS_PE_FACTORS
                        .get(&fuel_code.to_string())
                        .unwrap_or_else(|| {
                            panic!("Expected factor values in the table for the fuel code {fuel_code} were not present.");
                        });
                    (
                        vec![factor.emissions_factor],
                        vec![factor.emissions_factor_including_out_of_scope_emissions],
                        vec![factor.primary_energy_factor],
                    )
                }
            };

        let energy_supply_key = &KeyString::from(&energy_supply_key).unwrap();

        // Calculate energy imported and associated emissions/PE
        if fuel_code == FuelType::Electricity {
            supply_emis_result.import = apply_energy_factor_series(
                &energy_import[energy_supply_key],
                &emis_factor_import_export,
            )?;
            supply_emis_oos_result.import = apply_energy_factor_series(
                &energy_import[energy_supply_key],
                &emis_oos_factor_import_export,
            )?;
            supply_pe_result.import = apply_energy_factor_series(
                &energy_import[energy_supply_key],
                &pe_factor_import_export,
            )?;
        } else {
            supply_emis_result.import = energy_import[energy_supply_key]
                .iter()
                .map(|x| x * emis_factor_import_export[0])
                .collect::<Vec<_>>();
            supply_emis_oos_result.import = energy_import[energy_supply_key]
                .iter()
                .map(|x| x * emis_oos_factor_import_export[0])
                .collect::<Vec<_>>();
            supply_pe_result.import = energy_import[energy_supply_key]
                .iter()
                .map(|x| x * pe_factor_import_export[0])
                .collect::<Vec<_>>();
        }

        // If there is any export, Calculate energy exported and associated emissions/PE
        // Note that by convention, exported energy is negative
        (
            supply_emis_result.export,
            supply_emis_oos_result.export,
            supply_pe_result.export,
        ) = if energy_export[energy_supply_key].iter().sum::<f64>() < 0. {
            match fuel_code {
                FuelType::Electricity => (
                    apply_energy_factor_series(
                        &energy_export[energy_supply_key],
                        &emis_factor_import_export,
                    )?,
                    apply_energy_factor_series(
                        &energy_export[energy_supply_key],
                        &emis_oos_factor_import_export,
                    )?,
                    apply_energy_factor_series(
                        &energy_export[energy_supply_key],
                        &pe_factor_import_export,
                    )?,
                ),
                _ => (
                    energy_export[energy_supply_key]
                        .iter()
                        .map(|x| x * emis_factor_import_export[0])
                        .collect::<Vec<_>>(),
                    energy_export[energy_supply_key]
                        .iter()
                        .map(|x| x * emis_oos_factor_import_export[0])
                        .collect::<Vec<_>>(),
                    energy_export[energy_supply_key]
                        .iter()
                        .map(|x| x * pe_factor_import_export[0])
                        .collect::<Vec<_>>(),
                ),
            }
        } else {
            (
                vec![0.; number_of_timesteps],
                vec![0.; number_of_timesteps],
                vec![0.; number_of_timesteps],
            )
        };

        // Calculate energy generated and associated emissions/PE
        let mut energy_generated = vec![0.; number_of_timesteps];
        for end_user_energy in results_end_user[energy_supply_key].values() {
            if end_user_energy.iter().sum::<f64>() < 0. {
                for (t_idx, energy_generated_value) in energy_generated.iter_mut().enumerate() {
                    *energy_generated_value -= end_user_energy[t_idx];
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
                vec![0.; number_of_timesteps],
                vec![0.; number_of_timesteps],
                vec![0.; number_of_timesteps],
            )
        };

        // Calculate unregulated energy demand and associated emissions/PE
        let mut energy_unregulated = vec![0.; number_of_timesteps];
        for (end_user_name, end_user_energy) in results_end_user[energy_supply_key].iter() {
            if [APPL_OBJ_NAME, ELEC_COOK_OBJ_NAME, GAS_COOK_OBJ_NAME]
                .contains(&end_user_name.as_str())
            {
                for (t_idx, energy_unregulated_value) in energy_unregulated.iter_mut().enumerate() {
                    *energy_unregulated_value += end_user_energy[t_idx];
                }
            }
        }
        if fuel_code == FuelType::Electricity {
            supply_emis_result.unregulated =
                apply_energy_factor_series(&energy_unregulated, &emis_factor_import_export)?;
            supply_emis_oos_result.unregulated =
                apply_energy_factor_series(&energy_unregulated, &emis_oos_factor_import_export)?;
            supply_pe_result.unregulated =
                apply_energy_factor_series(&energy_unregulated, &pe_factor_import_export)?;
        } else {
            supply_emis_result.unregulated = energy_unregulated
                .iter()
                .map(|x| x * emis_factor_import_export[0])
                .collect::<Vec<_>>();
            supply_emis_oos_result.unregulated = energy_unregulated
                .iter()
                .map(|x| x * emis_oos_factor_import_export[0])
                .collect::<Vec<_>>();
            supply_pe_result.unregulated = energy_unregulated
                .iter()
                .map(|x| x * pe_factor_import_export[0])
                .collect::<Vec<_>>();
        }

        // Calculate total CO2/PE for each EnergySupply based on import and export,
        // subtracting unregulated
        supply_emis_result.total = Vec::with_capacity(number_of_timesteps);
        supply_emis_oos_result.total = Vec::with_capacity(number_of_timesteps);
        supply_pe_result.total = Vec::with_capacity(number_of_timesteps);
        for t_idx in 0..number_of_timesteps {
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

    Ok(FinalRates {
        emission_rate: total_emissions_rate,
        primary_energy_rate: total_pe_rate,
        emissions_results: emis_results,
        emissions_out_of_scope_results: emis_oos_results,
        primary_energy_results: pe_results,
    })
}

pub(super) struct FinalRates {
    pub(super) emission_rate: f64,
    pub(super) primary_energy_rate: f64,
    pub(super) emissions_results: IndexMap<String, FhsCalculationResult>,
    pub(super) emissions_out_of_scope_results: IndexMap<String, FhsCalculationResult>,
    pub(super) primary_energy_results: IndexMap<String, FhsCalculationResult>,
}

#[derive(Default)]
pub(super) struct FhsCalculationResult {
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

    let writer = output.writer_for_location_key(&file_location, "csv")?;
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

    let writer = output.writer_for_location_key("postproc_summary", "csv")?;
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

pub(super) fn calc_nbeds(input: &InputForProcessing) -> anyhow::Result<usize> {
    match input.number_of_bedrooms() {
        Some(bedrooms) => Ok(bedrooms),
        None => bail!("missing NumberOfBedrooms - required for FHS calculation"),
    }
}

pub(super) fn calc_n_occupants(
    total_floor_area: f64,
    number_of_bedrooms: usize,
) -> anyhow::Result<f64> {
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

fn create_occupancy(n_occupants: f64, occupancy_fhs: [f64; 24]) -> ([f64; 24], [f64; 24]) {
    let schedule_occupancy_weekday = occupancy_fhs.map(|factor| factor * n_occupants);
    let schedule_occupancy_weekend = occupancy_fhs.map(|factor| factor * n_occupants);

    (schedule_occupancy_weekday, schedule_occupancy_weekend)
}

fn create_metabolic_gains(
    number_of_occupants: f64,
    input: &mut InputForProcessing,
) -> anyhow::Result<()> {
    // Calculate total body surface area of occupants
    let a = 2.0001;
    let b = 0.8492;
    let total_body_surface_area_occupants = a * number_of_occupants.powf(b);

    let metabolic_gains_weekday_absolute = METABOLIC_GAINS
        .weekday
        .map(|gains| gains * total_body_surface_area_occupants)
        .to_vec();
    let metabolic_gains_weekend_absolute = METABOLIC_GAINS
        .weekend
        .map(|gains| gains * total_body_surface_area_occupants)
        .to_vec();

    input.set_metabolic_gains(
        0,
        0.5,
        json!(
            {
                "main": [{"repeat": 53, "value": "week"}],
                "week": [{"repeat": 5, "value": "weekday"}, {"repeat": 2, "value": "weekend"}],
                "weekday": metabolic_gains_weekday_absolute,
                "weekend": metabolic_gains_weekend_absolute,
            }
        ),
    )?;

    Ok(())
}

fn load_metabolic_gains_profile(file: impl Read) -> anyhow::Result<([f64; 48], [f64; 48])> {
    let mut metabolic_gains_reader = Reader::from_reader(BufReader::new(file));
    let rows: Vec<DryMetabolicGainsRow> = metabolic_gains_reader
        .deserialize()
        .collect::<Result<Vec<DryMetabolicGainsRow>, _>>()?;
    Ok(rows
        .iter()
        .enumerate()
        .fold(([0.; 48], [0.; 48]), |mut acc, (i, item)| {
            acc.0[i] = item.weekday;
            acc.1[i] = item.weekend;
            acc
        }))
}

#[derive(Deserialize)]
#[serde(rename = "lowercase")]
struct DryMetabolicGainsRow {
    #[serde(rename = "half_hour")]
    _half_hour: usize,
    #[serde(alias = "Weekday")]
    weekday: f64,
    #[serde(alias = "Weekend")]
    weekend: f64,
}

/// Space heating.
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
                match space_heat_system {
                    SystemReference::Single(space_heat_system) => {
                        input.set_control_string_for_space_heat_system(space_heat_system.as_str(), living_room_space_heat_system_name)?;
                        let control_schedule = living_room_control.as_object_mut().unwrap().get_mut("schedule").unwrap().as_object_mut().unwrap();
                        if let Some(temp_setback) = input.temperature_setback_for_space_heat_system(space_heat_system.as_str())? {
                            control_schedule.insert("setpoint_min".to_string(), temp_setback.into());
                        }
                        if let Some(advanced_start) = input.advanced_start_for_space_heat_system(space_heat_system.as_str())? {
                            control_schedule.insert("advanced_start".to_string(), advanced_start.into());
                        }
                    }
                    SystemReference::Multiple(_) => bail!("Multiple space heat system references under zone not currently supported for FHS inputs"),
                    SystemReference::None(_) => {}
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
                match space_heat_system {
                    SystemReference::Single(space_heat_system) => {
                        input.set_control_string_for_space_heat_system(space_heat_system.as_str(), rest_of_dwelling_space_heat_system_name)?;
                        let control_schedule = rest_of_dwelling_control.as_object_mut().unwrap().get_mut("schedule").unwrap().as_object_mut().unwrap();
                        if let Some(temp_setback) = input.temperature_setback_for_space_heat_system(space_heat_system.as_str())? {
                            control_schedule.insert("setpoint_min".to_string(), temp_setback.into());
                        }
                        if let Some(advanced_start) = input.advanced_start_for_space_heat_system(space_heat_system.as_str())? {
                            control_schedule.insert("advanced_start".to_string(), advanced_start.into());
                        }
                    }
                    SystemReference::Multiple(_) => bail!("Multiple space heat systems defined on a zone are not currently supported for FHS inputs"),
                    SystemReference::None(_) => {}
                }
                input.add_control(rest_of_dwelling_space_heat_system_name, rest_of_dwelling_control)?;
            }
            (None, _) => bail!("FHS does not yet have a condition to deal with zone that doesn't have specified living room/rest of dwelling"),
        }
    }

    Ok(())
}

#[derive(Clone, Copy, Debug)]
enum ControlType {
    Type2,
    Type3,
}

/// water heating pattern - if system is not instantaneous, hold at setpoint
/// 00:00-02:00 and then reheat as necessary 24/7
/// Note: Holding at setpoint for two hours has been chosen because
/// typical setting is for sterilisation cycle to last one hour, but the
/// model can only set a maximum and minimum setpoint temperaure, not
/// guarantee that the temperature is actually reached. Therefore, setting
/// the minimum to the maximum for two hours allows time for the tank
/// to heat up to the required temperature before being held there
fn create_water_heating_pattern(input: &mut InputForProcessing) -> anyhow::Result<()> {
    let hw_min_temp = "_HW_min_temp";
    let hw_max_temp = "_HW_max_temp";

    input.add_control(
        hw_min_temp,
        json!({
            "type": "SetpointTimeControl",
            "start_day": 0,
            "time_series_step": 0.5,
            "schedule": {
                "main": [{"value": "day", "repeat": 365}],
                "day": [{"value": HW_SETPOINT_MAX, "repeat": 4},{"value": HW_TEMPERATURE, "repeat": 44}]
            }
        }),
    )?;
    input.add_control(
        hw_max_temp,
        json!({
            "type": "SetpointTimeControl",
            "start_day": 0,
            "time_series_step": 0.5,
            "schedule": {
                "main": [{"value": "day", "repeat": 365}],
                "day": [{"value": HW_SETPOINT_MAX, "repeat": 48}]
            }
        }),
    )?;

    for hwsource in input.hot_water_source_keys() {
        let source = input.hot_water_source_details_for_key(hwsource.as_str());
        if source.is_storage_tank() {
            source.set_control_min_name_for_storage_tank_heat_sources(hw_min_temp)?;
            source.set_control_max_name_for_storage_tank_heat_sources(hw_max_temp)?;
        } else if source.is_combi_boiler() || source.is_point_of_use() || source.is_hiu() {
            // do nothing
        } else {
            bail!("Standard water heating schedule not defined for HotWaterSource type")
        }
    }

    Ok(())
}

/// Load the daily evaporative profile from a CSV file.
///
/// This function reads a CSV file containing time-of-day factors for evaporative losses
/// for each day of the week. It constructs a dictionary mapping days of the week to
/// lists of evaporative loss factors.
///
/// Arguments:
///
///  * `file` - The name of the CSV file containing the evaporative profile data.
///
///  Returns:
///     dict: A dictionary with days of the week as keys and lists of float factors as values.
fn load_evaporative_profile(file: impl Read) -> anyhow::Result<HalfHourWeeklyProfileData> {
    let mut profile_reader = Reader::from_reader(BufReader::new(file));

    let rows = profile_reader
        .deserialize()
        .collect::<Result<Vec<HalfHourWeeklyProfile>, _>>()
        .map_err(|_| anyhow!("Could not read evaporative profile file."))?;

    let (monday, tuesday, wednesday, thursday, friday, saturday, sunday) =
        rows.iter().enumerate().fold(
            (
                [0.; 48], [0.; 48], [0.; 48], [0.; 48], [0.; 48], [0.; 48], [0.; 48],
            ),
            |mut acc, (i, item)| {
                acc.0[i] = item.monday;
                acc.1[i] = item.tuesday;
                acc.2[i] = item.wednesday;
                acc.3[i] = item.thursday;
                acc.4[i] = item.friday;
                acc.5[i] = item.saturday;
                acc.6[i] = item.sunday;
                acc
            },
        );

    Ok(HalfHourWeeklyProfileData {
        monday,
        tuesday,
        wednesday,
        thursday,
        friday,
        saturday,
        sunday,
    })
}

#[derive(Debug, Deserialize)]
struct HalfHourWeeklyProfile {
    #[serde(rename = "Half_hour")]
    _half_hour: usize,
    #[serde(rename = "Mon")]
    monday: f64,
    #[serde(rename = "Tue")]
    tuesday: f64,
    #[serde(rename = "Wed")]
    wednesday: f64,
    #[serde(rename = "Thu")]
    thursday: f64,
    #[serde(rename = "Fri")]
    friday: f64,
    #[serde(rename = "Sat")]
    saturday: f64,
    #[serde(rename = "Sun")]
    sunday: f64,
}

struct HalfHourWeeklyProfileData {
    monday: [f64; 48],
    tuesday: [f64; 48],
    wednesday: [f64; 48],
    thursday: [f64; 48],
    friday: [f64; 48],
    saturday: [f64; 48],
    sunday: [f64; 48],
}

/// Apply the evaporative loss profile to modify the base evaporative loss across a full year.
///
/// This function takes the base evaporative loss and modifies it according to the provided
/// daily profile for each day of the week. It extends this profile throughout the year,
/// adjusting for any discrepancies in the week cycle (e.g., leap years).
///
/// Arguments:
///     * `input` - The main project dictionary where results are stored.
///     * `total_floor_area` - Total floor area used in the base loss calculation.
///     * `number_of_occupants` - Number of occupants used in the base loss calculation.
///     * `evaporative_profile_data` - Daily evaporative loss profiles loaded from a CSV file.
///
/// Effects:
///     Modifies the input in-place by setting a detailed schedule for evaporative losses.
fn create_evaporative_losses(
    input: &mut InputForProcessing,
    _total_floor_area: f64,
    number_of_occupants: f64,
    evaporative_profile_data: &HalfHourWeeklyProfileData,
) -> anyhow::Result<()> {
    // Base evaporative loss calculation
    let evaporative_losses_fhs = -25. * number_of_occupants;

    // Prepare to populate a full-year schedule of gains adjusted by the profile
    let mut evaporative_losses_schedule: Vec<f64> = Vec::with_capacity(18000);

    // Repeat for each week in a standard year
    evaporative_losses_schedule.extend(
        evaporative_profile_data
            .monday
            .iter()
            .chain(evaporative_profile_data.tuesday.iter())
            .chain(evaporative_profile_data.wednesday.iter())
            .chain(evaporative_profile_data.thursday.iter())
            .chain(evaporative_profile_data.friday.iter())
            .chain(evaporative_profile_data.saturday.iter())
            .chain(evaporative_profile_data.sunday.iter())
            .map(|factor| evaporative_losses_fhs * factor)
            .cycle()
            .take(48 * 7 * 52), // number of half-hour periods in 52 weeks
    );

    // Handle the extra days in the year not covered by the full weeks
    // Adjust based on the year (e.g., extra Monday for leap years)
    evaporative_losses_schedule.extend(
        evaporative_profile_data
            .monday
            .iter()
            .map(|factor| evaporative_losses_fhs * factor),
    );

    input.set_evaporative_losses(
        0,
        0.5,
        json!({
            "main": evaporative_losses_schedule,
        }),
    )?;

    Ok(())
}

/// Apply the cold water loss profile to modify the base cold water loss across a full year.
///
/// This function takes the base cold water loss and modifies it according to the provided
/// daily profile for each day of the week. It extends this profile throughout the year,
/// adjusting for any discrepancies in the weekly cycle (e.g., leap years).
///
/// Arguments:
///     * `input` - The main project dictionary where results are stored.
///     * `total_floor_area` - Total floor area used in the base loss calculation.
///     * `number_of_occupants` - Number of occupants used in the base loss calculation.
///     * `cold_water_loss_profile_data` - Daily cold water loss profiles loaded from a CSV file.
///
/// Effects:
///    Modifies the project_dict in-place by setting a detailed schedule for cold water losses.
fn create_cold_water_losses(
    input: &mut InputForProcessing,
    _total_floor_area: f64,
    number_of_occupants: f64,
    cold_water_loss_profile_data: &HalfHourWeeklyProfileData,
) -> anyhow::Result<()> {
    // Base cold water loss calculation
    let cold_water_losses_fhs = -20. * number_of_occupants;

    // Prepare to populate a full-year schedule of gains adjusted by the profile
    let mut cold_water_losses_schedule: Vec<f64> = Vec::with_capacity(18000);

    // Repeat for each week in a standard year
    cold_water_losses_schedule.extend(
        cold_water_loss_profile_data
            .monday
            .iter()
            .chain(cold_water_loss_profile_data.tuesday.iter())
            .chain(cold_water_loss_profile_data.wednesday.iter())
            .chain(cold_water_loss_profile_data.thursday.iter())
            .chain(cold_water_loss_profile_data.friday.iter())
            .chain(cold_water_loss_profile_data.saturday.iter())
            .chain(cold_water_loss_profile_data.sunday.iter())
            .map(|factor| cold_water_losses_fhs * factor)
            .cycle()
            .take(48 * 7 * 52), // number of half-hour periods in 52 weeks
    );

    // Handle the extra days in the year not covered by the full weeks
    // Adjust based on the year (e.g., extra Monday for leap years)
    cold_water_losses_schedule.extend(
        cold_water_loss_profile_data
            .monday
            .iter()
            .map(|factor| cold_water_losses_fhs * factor),
    );

    input.set_cold_water_losses(
        0,
        0.5,
        json!({
            "main": cold_water_losses_schedule
        }),
    )?;

    Ok(())
}

fn load_appliance_propensities(
    file: impl Read,
) -> anyhow::Result<AppliancePropensities<Normalised>> {
    let mut propensities_reader = Reader::from_reader(BufReader::new(file));
    let appliance_propensities_rows: Vec<AppliancePropensityRow> = propensities_reader
        .deserialize()
        .collect::<Result<Vec<AppliancePropensityRow>, _>>()
        .expect("Could not parse out appliance propensities CSV file correctly.");

    let (
        hour,
        occupied,
        cleaning_washing_machine,
        cleaning_tumble_dryer,
        cleaning_dishwasher,
        cooking_electric_oven,
        cooking_microwave,
        cooking_kettle,
        cooking_gas_cooker,
        consumer_electronics,
    ): AppliancePropensitiesUnderConstruction = appliance_propensities_rows
        .iter()
        .enumerate()
        .fold(Default::default(), |acc, (i, item)| {
            let (
                mut hour,
                mut occupied,
                mut cleaning_washing_machine,
                mut cleaning_tumble_dryer,
                mut cleaning_dishwasher,
                mut cooking_electric_oven,
                mut cooking_microwave,
                mut cooking_kettle,
                mut cooking_gas_cooker,
                mut consumer_electronics,
            ) = acc;
            hour[i] = item.hour;
            occupied[i] = item.occupied;
            cleaning_washing_machine[i] = item.cleaning_washing_machine;
            cleaning_tumble_dryer[i] = item.cleaning_tumble_dryer;
            cleaning_dishwasher[i] = item.cleaning_dishwasher;
            cooking_electric_oven[i] = item.cooking_electric_oven;
            cooking_microwave[i] = item.cooking_microwave;
            cooking_kettle[i] = item.cooking_kettle;
            cooking_gas_cooker[i] = item.cooking_gas_cooker;
            consumer_electronics[i] = item.consumer_electronics;
            (
                hour,
                occupied,
                cleaning_washing_machine,
                cleaning_tumble_dryer,
                cleaning_dishwasher,
                cooking_electric_oven,
                cooking_microwave,
                cooking_kettle,
                cooking_gas_cooker,
                consumer_electronics,
            )
        });
    Ok(AppliancePropensities {
        hour,
        occupied,
        cleaning_washing_machine,
        cleaning_tumble_dryer,
        cleaning_dishwasher,
        cooking_electric_oven,
        cooking_microwave,
        cooking_kettle,
        cooking_gas_cooker,
        consumer_electronics,
        state: Default::default(),
    }
    .normalise())
}

type AppliancePropensitiesUnderConstruction = (
    [usize; 24],
    [f64; 24],
    [f64; 24],
    [f64; 24],
    [f64; 24],
    [f64; 24],
    [f64; 24],
    [f64; 24],
    [f64; 24],
    [f64; 24],
);

#[derive(Copy, Clone)]
struct AppliancePropensities<T> {
    hour: [usize; 24],
    occupied: [f64; 24],
    cleaning_washing_machine: [f64; 24],
    cleaning_tumble_dryer: [f64; 24],
    cleaning_dishwasher: [f64; 24],
    cooking_electric_oven: [f64; 24],
    cooking_microwave: [f64; 24],
    cooking_kettle: [f64; 24],
    cooking_gas_cooker: [f64; 24],
    consumer_electronics: [f64; 24],
    state: PhantomData<T>,
}

impl AppliancePropensities<AsDataFile> {
    fn normalise(self) -> AppliancePropensities<Normalised> {
        let AppliancePropensities {
            cleaning_washing_machine,
            cleaning_tumble_dryer,
            cleaning_dishwasher,
            cooking_electric_oven,
            cooking_microwave,
            cooking_kettle,
            cooking_gas_cooker,
            consumer_electronics,
            ..
        } = self;

        let [cleaning_washing_machine, cleaning_tumble_dryer, cleaning_dishwasher, cooking_electric_oven, cooking_microwave, cooking_kettle, cooking_gas_cooker, consumer_electronics] =
            [
                cleaning_washing_machine,
                cleaning_tumble_dryer,
                cleaning_dishwasher,
                cooking_electric_oven,
                cooking_microwave,
                cooking_kettle,
                cooking_gas_cooker,
                consumer_electronics,
            ]
            .into_iter()
            .map(|probabilities| -> [f64; 24] {
                let sumcol = probabilities.iter().sum::<f64>();
                probabilities.map(|x| x / sumcol)
            })
            .collect::<Vec<_>>()
            .try_into()
            .expect("Problem normalising appliance propensities.");

        AppliancePropensities {
            hour: self.hour,
            occupied: self.occupied,
            cleaning_washing_machine,
            cleaning_tumble_dryer,
            cleaning_dishwasher,
            cooking_electric_oven,
            cooking_microwave,
            cooking_kettle,
            cooking_gas_cooker,
            consumer_electronics,
            state: Default::default(),
        }
    }
}

struct AsDataFile;
struct Normalised;

#[derive(Deserialize)]
struct AppliancePropensityRow {
    #[serde(rename = "Hour")]
    hour: usize,
    #[serde(rename = "Occupied prop ( Chance the house is occupied)")]
    occupied: f64,
    #[serde(rename = "Cleaning Washing machine Prop")]
    cleaning_washing_machine: f64,
    #[serde(rename = "Cleaning Tumble dryer")]
    cleaning_tumble_dryer: f64,
    #[serde(rename = "Cleaning Dishwasher")]
    cleaning_dishwasher: f64,
    #[serde(rename = "Cooking Electric Oven")]
    cooking_electric_oven: f64,
    #[serde(rename = "Cooking Microwave")]
    cooking_microwave: f64,
    #[serde(rename = "Cooking Kettle")]
    cooking_kettle: f64,
    #[serde(rename = "Cooking Gas Cooker")]
    cooking_gas_cooker: f64,
    #[serde(rename = "Consumer Electronics")]
    consumer_electronics: f64,
}

struct FlatAnnualPropensities {
    cleaning_washing_machine: Vec<f64>,
    cleaning_tumble_dryer: Vec<f64>,
    cleaning_dishwasher: Vec<f64>,
    cooking_electric_oven: Vec<f64>,
    cooking_microwave: Vec<f64>,
    cooking_kettle: Vec<f64>,
    cooking_gas_cooker: Vec<f64>,
    consumer_electronics: Vec<f64>,
}

impl From<&AppliancePropensities<Normalised>> for FlatAnnualPropensities {
    fn from(value: &AppliancePropensities<Normalised>) -> Self {
        let hours_in_year = HOURS_TO_END_DEC as usize;
        Self {
            cleaning_washing_machine: value
                .cleaning_washing_machine
                .into_iter()
                .cycle()
                .take(hours_in_year)
                .collect::<Vec<_>>(),
            cleaning_tumble_dryer: value
                .cleaning_tumble_dryer
                .into_iter()
                .cycle()
                .take(hours_in_year)
                .collect::<Vec<_>>(),
            cleaning_dishwasher: value
                .cleaning_dishwasher
                .into_iter()
                .cycle()
                .take(hours_in_year)
                .collect::<Vec<_>>(),
            cooking_electric_oven: value
                .cooking_electric_oven
                .into_iter()
                .cycle()
                .take(hours_in_year)
                .collect::<Vec<_>>(),
            cooking_microwave: value
                .cooking_microwave
                .into_iter()
                .cycle()
                .take(hours_in_year)
                .collect::<Vec<_>>(),
            cooking_kettle: value
                .cooking_kettle
                .into_iter()
                .cycle()
                .take(hours_in_year)
                .collect::<Vec<_>>(),
            cooking_gas_cooker: value
                .cooking_gas_cooker
                .into_iter()
                .cycle()
                .take(hours_in_year)
                .collect::<Vec<_>>(),
            consumer_electronics: value
                .consumer_electronics
                .into_iter()
                .cycle()
                .take(hours_in_year)
                .collect::<Vec<_>>(),
        }
    }
}

/// Calculate the annual energy requirement in kWh using the procedure described in SAP 10.2 up to and including step 9.
/// Divide this by 365 to get the average daily energy use.
/// Multiply the daily energy consumption figure by the following profiles to
/// create a daily profile for each month of the year (to be applied to all days in that month).
/// Multiply by the daylighting at each half hourly timestep to correct for incidence of daylight.
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
        lighting_efficacy +=
            zone_lighting_efficacy * input.area_for_zone(zone_key.as_str())? / total_floor_area;
    }
    if lighting_efficacy == 0. {
        bail!("invalid/missing lighting efficacy for all zones");
    }

    // from analysis of EFUS 2017 data (updated to derive from harmonic mean)
    let lumens = 1_139. * (total_floor_area * number_of_occupants).powf(0.39);
    let mut topup = top_up_lighting(input, lumens)?;
    topup /= 21.3; // assumed efficacy of top up lighting
    let topup_per_day = topup / 365_f64;

    // dropped 1/3 - 2/3 split based on SAP2012 assumptions about portable lighting
    let kwh_per_year = lumens / lighting_efficacy;
    let kwh_per_day = kwh_per_year / 365.;
    let factor = daylight_factor(input, total_floor_area)?;

    // Need to expand the monthly profiles to get an annual profile
    let annual_half_hour_profile: Vec<f64> = DAYS_IN_MONTH
        .iter()
        .enumerate()
        .flat_map(|(month, days)| (0..*days).map(move |_| month))
        .flat_map(|month| AVERAGE_MONTHLY_LIGHTING_HALF_HOUR_PROFILES[month])
        .collect();

    // for each half hour time step in annual_halfhr_profiles:
    // To obtain the lighting gains,
    // the above should be converted to Watts by multiplying the individual half-hourly figure by (2 x 1000).
    // Since some lighting energy will be used in external light
    // (e.g. outdoor security lights or lights in unheated spaces like garages and sheds)
    // a factor of 0.85 is also applied to get the internal gains from lighting.
    let (lighting_gains_w, topup_gains_w): (Vec<f64>, Vec<f64>) = annual_half_hour_profile
        .into_iter()
        .enumerate()
        .map(|(i, profile)| {
            (
                (profile * kwh_per_day * factor[i]) * 2. * 1_000.,
                (profile * topup_per_day * factor[i]) * 2. * 1_000.,
            )
        })
        .collect();

    input.clear_appliance_gains();
    input.set_lighting_gains(json!({
        "type": "lighting",
        "start_day": 0,
        "time_series_step": 0.5,
        "gains_fraction": 0.85,
        "EnergySupply": ENERGY_SUPPLY_NAME_ELECTRICITY,
        "schedule": {
            "main": lighting_gains_w
        },
        "priority": -1
    }))?;

    input.set_topup_gains(json!({
        "type": "lighting",
        "start_day": 0,
        "time_series_step": 0.5,
        "gains_fraction": 0.85,
        "EnergySupply": ENERGY_SUPPLY_NAME_ELECTRICITY,
        "schedule": {
            "main": topup_gains_w
        }
    }))?;

    Ok(())
}

fn create_appliance_gains(
    input: &mut InputForProcessing,
    total_floor_area: f64,
    number_of_occupants: f64,
    appliance_propensities: &AppliancePropensities<Normalised>,
) -> anyhow::Result<()> {
    // create flattened year-long versions of energy use profiles
    let flat_efus_profile: Vec<f64> = DAYS_IN_MONTH
        .iter()
        .enumerate()
        .flat_map(|(month, days)| (0..*days).map(move |_| month))
        .flat_map(|month| AVERAGE_MONTHLY_APPLIANCES_HOURLY_PROFILES[month])
        .collect();

    let flat_annual_propensities: FlatAnnualPropensities = appliance_propensities.into();

    // add any missing required appliances to the assessment,
    // get default demand figures for any unknown appliances
    appliance_cooking_defaults(input, number_of_occupants, total_floor_area);
    let cookparams = cooking_demand(input, number_of_occupants)?;

    // TODO (from Python) change to enum
    // TODO (from Python) check appliances are named correctly and what to do if not?

    let appliance_map: IndexMap<ApplianceKey, ApplianceUseProfile> = IndexMap::from([
        (
            ApplianceKey::Fridge,
            ApplianceUseProfile::simple(1., 0., 1.0, flat_efus_profile.clone()),
        ),
        (
            ApplianceKey::Freezer,
            ApplianceUseProfile::simple(1., 0., 1.0, flat_efus_profile.clone()),
        ),
        (
            ApplianceKey::FridgeFreezer,
            ApplianceUseProfile::simple(1., 0., 1.0, flat_efus_profile.clone()),
        ),
        (
            ApplianceKey::OtherDevices,
            ApplianceUseProfile::simple(
                1.,
                0.,
                1.0,
                flat_annual_propensities.consumer_electronics.clone(),
            ),
        ),
        (
            ApplianceKey::Dishwasher,
            ApplianceUseProfile::complex(
                number_of_occupants,
                132,       // HES 2012 final report table 22
                Some(280), // EU standard
                0.75,
                0.3,
                flat_annual_propensities.cleaning_dishwasher.clone(),
                1.5,
                0.,
            ),
        ),
        (
            ApplianceKey::ClothesWashing,
            ApplianceUseProfile::clothes(
                number_of_occupants,
                174, // HES 2012 final report table 22
                220, // EU standard
                7.,
                0.75,
                0.3,
                flat_annual_propensities.cleaning_washing_machine.clone(),
                2.5,
                0.,
            ),
        ),
        (
            ApplianceKey::ClothesDrying,
            ApplianceUseProfile::clothes(
                number_of_occupants,
                145, // HES 2012 final report table 22
                160, // EU standard
                7.,
                0.50,
                0.7,
                flat_annual_propensities.cleaning_tumble_dryer.clone(),
                0.75,
                0.,
            ),
        ),
        (
            ApplianceKey::Oven,
            ApplianceUseProfile::complex(
                1.,
                cookparams.get(&ApplianceKey::Oven).unwrap().event_count, // analysis of HES - see folder
                None,
                0.50,
                0.5,
                flat_annual_propensities.cooking_electric_oven.clone(),
                0.5,
                0.7,
            ),
        ),
        (
            ApplianceKey::Hobs,
            ApplianceUseProfile::complex(
                1.,
                cookparams.get(&ApplianceKey::Hobs).unwrap().event_count, // analysis of HES - see folder
                None,
                0.50,
                0.5,
                flat_annual_propensities.cooking_gas_cooker.clone(),
                0.1,
                0.7,
            ),
        ),
        (
            ApplianceKey::Microwave,
            ApplianceUseProfile::complex(
                1.,
                cookparams
                    .get(&ApplianceKey::Microwave)
                    .unwrap()
                    .event_count, // analysis of HES - see folder
                None,
                0.50,
                1.,
                flat_annual_propensities.cooking_microwave.clone(),
                0.05,
                0.3,
            ),
        ),
        (
            ApplianceKey::Kettle,
            ApplianceUseProfile::complex(
                1.,
                cookparams.get(&ApplianceKey::Kettle).unwrap().event_count, // analysis of HES - see folder
                None,
                0.50,
                1.,
                flat_annual_propensities.cooking_kettle.clone(),
                0.05,
                0.3,
            ),
        ),
    ]);

    // add any missing required appliances to the assessment,
    // get default demand figures for any unknown appliances
    let mut priority: IndexMap<ApplianceKey, (Option<isize>, f64)> = Default::default();
    let mut power_scheds: IndexMap<ApplianceKey, Vec<f64>> = Default::default();
    let mut weight_scheds: IndexMap<ApplianceKey, Vec<f64>> = Default::default();
    // loop through appliances in the assessment.
    let input_appliances = input.clone_appliances();

    for (appliance_key, appliance) in input_appliances {
        // if it needs to be modelled per use
        let map_appliance = appliance_map
            .get(&appliance_key)
            .expect("Appliance key was not in appliance map");

        if let Some(use_data) = map_appliance.use_data {
            // value on energy label is defined differently between appliance types
            // TODO (from Python) - translation of efficiencies should be its own function
            let (kwhcycle, loadingfactor) =
                appliance_kwh_cycle_loading_factor(input, &appliance_key, &appliance_map)?;

            let app = FhsAppliance::new(
                map_appliance.util_unit,
                use_data.use_metric as f64 * loadingfactor,
                kwhcycle,
                use_data.duration,
                map_appliance.standby,
                map_appliance.gains_frac,
                &map_appliance.prof,
                None,
                Some(use_data.duration_deviation),
            )?;

            let appliance_energy_supply = if let ApplianceEntry::Object(appliance_obj) = &appliance
            {
                appliance_obj.energy_supply
            } else {
                None
            };

            // if the appliance specifies load shifting, add it to the appliance gains details
            let load_shifting = if let ApplianceEntry::Object(Appliance {
                load_shifting: Some(load_shifting),
                ..
            }) = &appliance
            {
                if load_shifting.max_shift_hrs >= 24. {
                    // could instead change length of buffers/initial simulation match this, but unclear what benefit this would have
                    bail!("{} max_shift_hrs too high, FHS wrapper cannot handle max shift >= 24 hours", appliance_key);
                }

                // establish priority between appliances based on user defined priority,
                // and failing that, demand per cycle
                priority.insert(appliance_key, (load_shifting.priority, kwhcycle));

                let mut load_shifting = load_shifting.clone();
                // create year long cost profile
                // loadshifting is also intended to respond to CO2, primary energy factors instead of cost, for example
                // so the weight timeseries is generic.

                // TODO (Python) - create weight timeseries as combination of PE, CO2, cost factors.
                // could also multiply by propensity factor
                let weight_timeseries = reject_nulls(expand_numeric_schedule(
                    input.tariff_schedule().ok_or_else(|| {
                        anyhow!(
                            "A tariff schedule was expected to have been provided in the input."
                        )
                    })?,
                ))?;
                load_shifting.weight_timeseries = Some(weight_timeseries.clone());
                weight_scheds.insert(appliance_key, weight_timeseries);

                Some(load_shifting)
            } else {
                // only add demand from appliances that DO NOT have loadshifting to the demands
                power_scheds.insert(appliance_key, app.flat_schedule.clone());
                priority.insert(appliance_key, (None, kwhcycle));
                None
            };

            input.set_gains_for_field(String::from(appliance_key), json!({
                "type": appliance_key,
                "EnergySupply": if [ApplianceKey::Hobs, ApplianceKey::Oven].contains(&appliance_key) {
                    appliance_energy_supply.ok_or_else(|| anyhow!("Could not get energy supply type for appliance with key {appliance_key}"))?.to_string()
                } else {
                    ENERGY_SUPPLY_NAME_ELECTRICITY.to_owned()
                },
                "start_day": 0,
                // TODO (from Python) - variable timestep
                "time_series_step": 1,
                "gains_fraction": app.gains_frac,
                "Events": app.event_list,
                "Standby": app.standby_w,
                "loadshifting": load_shifting
            }))?;
        } else {
            // model as yearlong time series schedule of demand in W
            let annual_kwh = if let ApplianceEntry::Object(Appliance {
                kwh_per_annum: Some(ref kwh_per_annum),
                ..
            }) = &appliance
            {
                kwh_per_annum * map_appliance.util_unit
            } else {
                continue;
            };
            // TODO (from Python) - check normalisation of flat profile here
            let flat_schedule: Vec<f64> = flat_efus_profile
                .iter()
                .map(|&frac| WATTS_PER_KILOWATT as f64 / DAYS_PER_YEAR as f64 * frac * annual_kwh)
                .collect();
            power_scheds.insert(appliance_key, flat_schedule.clone());

            priority.insert(appliance_key, (None, 1.)); // Python passes in kwhcycle here instead of 1. but kwhcycle hasn't been assigned at this point, may be erroneous.

            let appliance_uses_gas: bool = false; // upstream Python checks appliance key contains substring 'gas', may be erroneous

            input.set_gains_for_field(String::from(appliance_key), json!({
                "type": appliance_key,
                "EnergySupply": if appliance_uses_gas { ENERGY_SUPPLY_NAME_GAS } else { ENERGY_SUPPLY_NAME_ELECTRICITY },
                "start_day": 0,
                "time_series_step": 1,
                "gains_fraction": map_appliance.gains_frac,
                "schedule": {
                   // watts
                   "main": flat_schedule
                }
            }))?;
        }
    }
    // sum schedules for use with loadshifting
    // will this work with variable timestep?
    let sched_len = power_scheds
        .values()
        .next()
        .ok_or_else(|| anyhow!("Demand schedules are empty"))?
        .len();

    let sched_zeros: Vec<f64> = vec![0.; sched_len];

    let mut main_power_sched: IndexMap<String, Vec<f64>> = IndexMap::from([
        (ENERGY_SUPPLY_NAME_GAS.to_string(), sched_zeros.clone()),
        (
            ENERGY_SUPPLY_NAME_ELECTRICITY.to_string(),
            sched_zeros.clone(),
        ),
    ]);

    let mut main_weight_sched: IndexMap<String, Vec<f64>> = IndexMap::from([
        (ENERGY_SUPPLY_NAME_GAS.to_string(), sched_zeros.clone()),
        (
            ENERGY_SUPPLY_NAME_ELECTRICITY.to_string(),
            sched_zeros.clone(),
        ),
    ]);

    for appliance_key in power_scheds.keys() {
        let energy_supply_name = input
            .energy_supply_type_for_appliance_gains_field(&appliance_key.to_string())
            .ok_or_else(|| {
                anyhow!(
                    "No energy supply type for appliance gains for {}",
                    appliance_key
                )
            })?;

        let main_power_schedule_for_energy_supply: &Vec<f64> =
            main_power_sched.get(&energy_supply_name).ok_or_else(|| {
                anyhow!(
                    "There was no main power schedule for energy supply {}",
                    energy_supply_name
                )
            })?;

        main_power_sched.insert(
            energy_supply_name,
            main_power_schedule_for_energy_supply
                .iter()
                .enumerate()
                .map(|(i, main_power)| main_power + power_scheds.get(appliance_key).unwrap()[i])
                .collect(),
        );
    }

    for appliance_key in weight_scheds.keys() {
        let energy_supply_name = input
            .energy_supply_type_for_appliance_gains_field(&appliance_key.to_string())
            .ok_or_else(|| {
                anyhow!(
                    "No energy supply type for appliance gains for {}",
                    appliance_key
                )
            })?;

        let main_weight_sched_for_energy_supply: &Vec<f64> =
            main_weight_sched.get(&energy_supply_name).ok_or_else(|| {
                anyhow!(
                    "There was no main power schedule for energy supply {}",
                    energy_supply_name
                )
            })?;

        main_weight_sched.insert(
            energy_supply_name,
            main_weight_sched_for_energy_supply
                .iter()
                .enumerate()
                .map(|(i, main_weight)| main_weight + weight_scheds.get(appliance_key).unwrap()[i])
                .collect(),
        );
    }

    let smart_control: ControlDetails = ControlDetails::SmartAppliance {
        battery_24hr: None,
        non_appliance_demand_24hr: None,
        power_timeseries: main_power_sched,
        time_series_step: 1.,
        weight_timeseries: main_weight_sched,
    };

    input.set_loadshifting_control(smart_control);

    // work out order in which to process loadshifting appliances
    let defined_priority = priority
        .iter()
        .filter_map(|(appliance_name, priorities)| priorities.0.map(|_| appliance_name))
        .collect::<Vec<&ApplianceKey>>();

    let mut first_priority_ranks: Vec<isize> = defined_priority
        .iter()
        .filter_map(|appliance_name| priority.get(appliance_name.to_owned()))
        .filter_map(|p| p.0)
        .collect_vec();

    first_priority_ranks.append(&mut vec![0]);

    let lowest_priority = first_priority_ranks.iter().max().unwrap();

    let priority_kwhcycle: Vec<ApplianceKey> = priority
        .clone()
        .sorted_by(|_, (_, kwhcycle1), _, (_, kwhcycle2)| kwhcycle1.total_cmp(kwhcycle2))
        .filter(|(_, (priority, _))| priority.is_none())
        .rev()
        .map(|x| x.0)
        .collect();

    for appliance in priority.keys() {
        let new_priority = if defined_priority.contains(&appliance) {
            defined_priority
                .iter()
                .position(|&a| a == appliance)
                .unwrap() as isize
        } else {
            priority_kwhcycle
                .iter()
                .position(|a| a == appliance)
                .unwrap() as isize
                + *lowest_priority
        };
        input.set_priority_for_gains_appliance(new_priority, appliance)?;
    }

    Ok(())
}

#[derive(Clone, Debug)]
struct ApplianceUseProfile {
    util_unit: f64,
    use_data: Option<ApplianceUseData>,
    standby: f64,
    gains_frac: f64,
    prof: Vec<f64>,
}

impl ApplianceUseProfile {
    fn simple(util_unit: f64, standby: f64, gains_frac: f64, prof: Vec<f64>) -> Self {
        Self {
            util_unit,
            use_data: None,
            standby,
            gains_frac,
            prof,
        }
    }

    fn complex(
        util_unit: f64,
        use_metric: usize,
        standard_use: Option<usize>,
        standby: f64,
        gains_frac: f64,
        prof: Vec<f64>,
        duration: f64,
        duration_deviation: f64,
    ) -> Self {
        Self {
            util_unit,
            use_data: Some(ApplianceUseData {
                use_metric,
                clothes_use_data: None,
                _standard_use: standard_use,
                duration,
                duration_deviation,
            }),
            standby,
            gains_frac,
            prof,
        }
    }

    fn clothes(
        util_unit: f64,
        use_metric: usize,
        standard_use: usize,
        standard_load_kg: f64,
        standby: f64,
        gains_frac: f64,
        prof: Vec<f64>,
        duration: f64,
        duration_deviation: f64,
    ) -> Self {
        Self {
            util_unit,
            use_data: Some(ApplianceUseData {
                use_metric,
                clothes_use_data: Some(ClothesUseData { standard_load_kg }),
                _standard_use: Some(standard_use),
                duration,
                duration_deviation,
            }),
            standby,
            gains_frac,
            prof,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct ApplianceUseData {
    // maps to "use" field in upstream, though 'use' is a keywork in Rust so calling this "use_metric"
    use_metric: usize,
    clothes_use_data: Option<ClothesUseData>,
    _standard_use: Option<usize>,
    duration: f64,
    duration_deviation: f64,
}

#[derive(Clone, Copy, Debug)]
struct ClothesUseData {
    standard_load_kg: f64,
}
struct ApplianceCookingDemand {
    mean_annual_demand: f64,
    _mean_annual_events: f64,
    mean_event_demand: f64,
    fuel: Option<String>,
    event_count: usize,
}

fn cooking_demand(
    input: &mut InputForProcessing,
    number_of_occupants: f64,
) -> anyhow::Result<IndexMap<ApplianceKey, ApplianceCookingDemand>> {
    let oven_energy_supply = input.energy_supply_for_appliance(&ApplianceKey::Oven);
    let oven_fuel = match oven_energy_supply {
        Ok(energy_supply) => {
            Some(input.fuel_type_for_energy_supply_field(&energy_supply.to_string())?)
        }
        Err(_) => None,
    };
    let oven = ApplianceCookingDemand {
        mean_annual_demand: 285.14,
        _mean_annual_events: 441.11,
        mean_event_demand: 0.762,
        fuel: oven_fuel,
        event_count: Default::default(),
    };

    let hobs_energy_supply = input.energy_supply_for_appliance(&ApplianceKey::Hobs);
    let hobs_fuel = match hobs_energy_supply {
        Ok(energy_supply) => {
            Some(input.fuel_type_for_energy_supply_field(&energy_supply.to_string())?)
        }
        Err(_) => None,
    };
    let hobs = ApplianceCookingDemand {
        mean_annual_demand: 352.53,
        _mean_annual_events: 520.86,
        mean_event_demand: 0.810,
        fuel: hobs_fuel,
        event_count: Default::default(),
    };

    let microwave_fuel = match input.appliances_contain_key(&ApplianceKey::Microwave) {
        true => Some(EnergySupplyType::Electricity.to_string()),
        false => None,
    };
    let microwave = ApplianceCookingDemand {
        mean_annual_demand: 44.11,
        _mean_annual_events: 710.65,
        mean_event_demand: 0.0772,
        fuel: microwave_fuel,
        event_count: Default::default(),
    };

    let kettle_fuel = match input.appliances_contain_key(&ApplianceKey::Kettle) {
        true => Some(EnergySupplyType::Electricity.to_string()),
        false => None,
    };
    let kettle = ApplianceCookingDemand {
        mean_annual_demand: 173.03,
        _mean_annual_events: 1782.5,
        mean_event_demand: 0.0985,
        fuel: kettle_fuel,
        event_count: Default::default(),
    };
    let mut cook_params = IndexMap::from([
        (ApplianceKey::Oven, oven),
        (ApplianceKey::Hobs, hobs),
        (ApplianceKey::Microwave, microwave),
        (ApplianceKey::Kettle, kettle),
    ]);

    let gas_total: f64 = cook_params
        .values()
        .filter(|appliance_details| {
            appliance_details
                .fuel
                .as_ref()
                .is_some_and(|fuel| *fuel == FuelType::MainsGas.to_string())
        })
        .map(|appliance_details| appliance_details.mean_annual_demand)
        .sum();

    let elec_total: f64 = cook_params
        .values()
        .filter(|appliance_details| {
            appliance_details
                .fuel
                .as_ref()
                .is_some_and(|fuel| *fuel == FuelType::Electricity.to_string())
        })
        .map(|appliance_details| appliance_details.mean_annual_demand)
        .sum();

    // top down cooking demand estimate based on analysis of EFUS 2017 electricity monitoring data
    // and HES 2012
    let annual_cooking_elec_kwh = 448. * 0.8 + (171. + 98. * number_of_occupants) * 0.2;

    for cooking_demand in cook_params.values_mut() {
        // for each appliance, work out number of usage events based on
        // average HES annual demand and demand per cycle
        // do not consider gas and electricity separately for this purpose
        let demand_prop = cooking_demand.mean_annual_demand / (elec_total + gas_total);
        let annual_kwh = demand_prop * annual_cooking_elec_kwh;
        let events = annual_kwh / cooking_demand.mean_event_demand;
        cooking_demand.event_count = events as usize;
    }

    Ok(cook_params)
}

fn appliance_cooking_defaults(
    input: &mut InputForProcessing,
    number_of_occupants: f64,
    total_floor_area: f64,
) -> (
    IndexMap<ApplianceKey, Appliance>,
    IndexMap<ApplianceKey, Appliance>,
) {
    let cooking_fuels = input.all_energy_supply_fuel_types();

    // (from Python) also check gas/elec cooker/oven  together - better to have energysupply as a dict entry?
    let mut cooking_defaults: IndexMap<ApplianceKey, Appliance> = match (
        cooking_fuels.contains(&FuelType::Electricity),
        cooking_fuels.contains(&FuelType::MainsGas),
    ) {
        (true, true) => IndexMap::from([
            (
                ApplianceKey::Oven,
                Appliance::with_energy_supply(EnergySupplyType::Electricity, 0.59),
            ),
            (
                ApplianceKey::Hobs,
                Appliance::with_energy_supply(EnergySupplyType::MainsGas, 0.72),
            ),
        ]),
        (_, true) => IndexMap::from([
            (
                ApplianceKey::Oven,
                Appliance::with_energy_supply(EnergySupplyType::MainsGas, 1.57),
            ),
            (
                ApplianceKey::Hobs,
                Appliance::with_energy_supply(EnergySupplyType::MainsGas, 0.72),
            ),
        ]),
        (true, _) => IndexMap::from([
            (
                ApplianceKey::Oven,
                Appliance::with_energy_supply(EnergySupplyType::Electricity, 0.59),
            ),
            (
                ApplianceKey::Hobs,
                Appliance::with_energy_supply(EnergySupplyType::Electricity, 0.72),
            ),
        ]),
        _ => IndexMap::from([
            (
                ApplianceKey::Oven,
                Appliance::with_energy_supply(EnergySupplyType::Electricity, 0.59),
            ),
            (
                ApplianceKey::Hobs,
                Appliance::with_energy_supply(EnergySupplyType::Electricity, 0.72),
            ),
        ]),
    };

    let mut additional_cooking_defaults = IndexMap::from([
        (ApplianceKey::Kettle, Appliance::with_kwh_per_cycle(0.1)),
        (ApplianceKey::Microwave, Appliance::with_kwh_per_cycle(0.08)),
    ]);

    let appliance_defaults = IndexMap::from([
        (
            ApplianceKey::OtherDevices,
            Appliance::with_kwh_per_annum(
                30.0 * (number_of_occupants * total_floor_area).powf(0.49),
            ),
        ),
        (
            ApplianceKey::Dishwasher,
            Appliance::with_kwh_per_100_cycle(53.0, None),
        ),
        (
            ApplianceKey::ClothesWashing,
            Appliance::with_kwh_per_100_cycle(53.0, Some(7.0)),
        ),
        (
            ApplianceKey::ClothesDrying,
            Appliance::with_kwh_per_100_cycle(98.0, Some(7.0)),
        ),
        (ApplianceKey::Fridge, Appliance::with_kwh_per_annum(76.7)),
        (ApplianceKey::Freezer, Appliance::with_kwh_per_annum(128.2)),
        (
            ApplianceKey::FridgeFreezer,
            Appliance::with_kwh_per_annum(137.4),
        ),
    ]);

    if !input.has_appliances() {
        input.merge_in_appliances(&appliance_defaults);
        input.merge_in_appliances(&cooking_defaults);
        input.merge_in_appliances(&additional_cooking_defaults);
    } else {
        for appliance_name in appliance_defaults.keys() {
            if !input.appliances_contain_key(appliance_name)
                || input.appliance_key_has_reference(appliance_name, &ApplianceReference::Default)
            {
                input.merge_in_appliances(&IndexMap::from([(
                    appliance_name.to_owned(),
                    appliance_defaults[appliance_name].clone(),
                )]));
            } else if input
                .appliance_key_has_reference(appliance_name, &ApplianceReference::NotInstalled)
            {
                input.remove_appliance(appliance_name);
            } else {
                // user has specified appliance efficiency, overwrite efficiency with default
                let original_load_shifting_value = input.loadshifting_for_appliance(appliance_name);

                input.merge_in_appliances(&IndexMap::from([(
                    appliance_name.to_owned(),
                    appliance_defaults[appliance_name].clone(),
                )]));
                if let Some(load_shifting) = original_load_shifting_value {
                    input.set_loadshifting_for_appliance(appliance_name, load_shifting);
                }
            }
        }
        if !cooking_defaults
            .keys()
            .any(|cooking_appliance_name| input.appliances_contain_key(cooking_appliance_name))
        {
            // neither cooker nor oven specified, add cooker as minimum requirement
            input.merge_in_appliances(&IndexMap::from([(
                ApplianceKey::Hobs,
                cooking_defaults[&ApplianceKey::Hobs].clone(),
            )]));
        }
        cooking_defaults.append(&mut additional_cooking_defaults);
        for (cooking_name, cooking_appliance) in cooking_defaults.iter() {
            if !input.appliances_contain_key(cooking_name)
                || input.appliance_key_has_reference(cooking_name, &ApplianceReference::Default)
            {
                input.merge_in_appliances(&IndexMap::from([(
                    cooking_name.to_owned(),
                    cooking_appliance.clone(),
                )]));
            } else if input
                .appliance_key_has_reference(cooking_name, &ApplianceReference::NotInstalled)
            {
                input.remove_appliance(cooking_name);
            } else {
                // NB: there is a possible issue in the Python here where the wrong key is used
                input.merge_in_appliances(&IndexMap::from([(
                    cooking_name.to_owned(),
                    cooking_appliance.clone(),
                )]));
            }
        }
    }

    (appliance_defaults, cooking_defaults)
}

fn appliance_kwh_cycle_loading_factor(
    input: &InputForProcessing,
    appliance_key: &ApplianceKey,
    appliance_map: &IndexMap<ApplianceKey, ApplianceUseProfile>,
) -> anyhow::Result<(f64, f64)> {
    // value on energy label is defined differently between appliance types,
    // convert any different input types to simple kWh per cycle

    let (kwh_cycle, appliance) = if let Some(ApplianceEntry::Object(appliance)) =
        input.appliance_with_key(appliance_key)
    {
        (
            (if let Some(kwh_per_cycle) = appliance.kwh_per_cycle {
                kwh_per_cycle
            } else if let Some(kwh_per_100_cycle) = appliance.kwh_per_100_cycle {
                kwh_per_100_cycle / 100.
            } else if let (Some(kwh_per_annum), Some(standard_use)) =
                (appliance.kwh_per_annum, appliance.standard_use)
            {
                kwh_per_annum / standard_use
            } else {
                bail!("{} demand must be specified as one of 'kWh_per_cycle', 'kWh_per_100cycle' or 'kWh_per_annum'", appliance_key)
            }),
            appliance,
        )
    } else {
        bail!("Appliance with name '{appliance_key}' must exist.")
    };

    let map_appliance = appliance_map.get(appliance_key).ok_or_else(|| anyhow!("The appliance name '{appliance_key}' was expected to be found within the appliance map: {appliance_map:?}."))?;

    let mut loading_factor = 1.;
    if appliance_key.is_clothes_appliance() {
        // additionally, laundry appliances have variable load size,
        // which affects the required number of uses to do all the occupants' laundry for the year
        loading_factor = {
            map_appliance
                .use_data
                .as_ref()
                .ok_or_else(|| anyhow!("Appliance is expected to have clothes use data"))?
                .clothes_use_data
                .as_ref()
                .ok_or_else(|| anyhow!("Appliance is expected to have clothes use data"))?
                .standard_load_kg
                / appliance.kg_load.as_ref().ok_or_else(|| {
                    anyhow!("Passed in appliance is expected to have a kg_load value.")
                })?
        };

        // There is some unreachable code in the Python here around spin dry efficiency class
        // TODO - implement in future if needed
    }

    Ok((kwh_cycle, loading_factor))
}

fn sim_24h(input: &mut InputForProcessing, sim_settings: SimSettings) -> anyhow::Result<()> {
    let mut input_24h = input.clone();
    let range = (HOURS_PER_DAY as f64 / SIMTIME_STEP).ceil() as usize;
    let zeros_24h_by_supply = IndexMap::from([
        (ENERGY_SUPPLY_NAME_ELECTRICITY.to_string(), vec![0.; range]),
        (ENERGY_SUPPLY_NAME_GAS.to_string(), vec![0.; range]),
    ]);

    input_24h.set_non_appliance_demand_24hr(zeros_24h_by_supply.clone())?;

    input_24h.set_battery24hr(SmartApplianceBattery {
        energy_into_battery_from_generation: zeros_24h_by_supply.clone(),
        energy_out_of_battery: zeros_24h_by_supply.clone(),
        energy_into_battery_from_grid: zeros_24h_by_supply.clone(),
        battery_state_of_charge: zeros_24h_by_supply.clone(),
    })?;

    input_24h.set_simulation_time(SimulationTime::new(
        SIMTIME_START,
        SIMTIME_START + HOURS_PER_DAY as f64,
        SIMTIME_STEP,
    ));

    // create a corpus instance
    let output_options = OutputOptions {
        print_heat_balance: sim_settings.heat_balance,
        detailed_output_heating_cooling: sim_settings.detailed_output_heating_cooling,
    };

    let corpus = Corpus::from_inputs(
        input_24h.as_input(),
        None,
        sim_settings.tariff_data_filename.as_deref(),
        &output_options,
    )?;

    // Run main simulation sim
    let results = corpus.run()?;

    // sum results for electricity demand other than appliances to get 24h demand buffer for loadshifting
    let electricity_users = results
        .results_end_user
        .get(ENERGY_SUPPLY_NAME_ELECTRICITY)
        .ok_or_else(|| anyhow!("Expected one or more users of mains elec energy supply"))?;

    let min_demand_length = electricity_users
        .values()
        .map(|demand| demand.len())
        .min()
        .unwrap_or(0);

    let mut non_appliance_electricity_demand = vec![];
    for i in 0..min_demand_length {
        for (name, user) in electricity_users {
            let name = ApplianceKey::try_from(name.as_str());
            let do_increment = match name {
                Ok(name) => !input.appliances_contain_key(&name),
                Err(_) => true,
            };
            if do_increment {
                non_appliance_electricity_demand.insert(
                    i,
                    non_appliance_electricity_demand
                        .get(i)
                        .unwrap_or(&0.)
                        .clone()
                        + user[i],
                );
            }
        }
    }

    let non_appliance_demand_24hr = IndexMap::from([
        (
            ENERGY_SUPPLY_NAME_ELECTRICITY.to_string(),
            non_appliance_electricity_demand,
        ),
        (ENERGY_SUPPLY_NAME_GAS.to_string(), vec![0.; range]),
    ]);
    input.set_non_appliance_demand_24hr(non_appliance_demand_24hr)?;

    let energy_into_battery_from_generation = results
        .energy_to_storage
        .iter()
        .map(|(key, value)| (key.to_string(), value.clone()))
        .collect::<IndexMap<String, Vec<f64>>>();
    let energy_out_of_battery = results
        .energy_from_storage
        .iter()
        .map(|(key, value)| (key.to_string(), value.clone()))
        .collect::<IndexMap<String, Vec<f64>>>();
    let energy_into_battery_from_grid = results
        .storage_from_grid
        .iter()
        .map(|(key, value)| (key.to_string(), value.clone()))
        .collect::<IndexMap<String, Vec<f64>>>();
    let battery_state_of_charge = results
        .battery_state_of_charge
        .iter()
        .map(|(key, value)| (key.to_string(), value.clone()))
        .collect::<IndexMap<String, Vec<f64>>>();

    input.set_battery24hr(SmartApplianceBattery {
        energy_into_battery_from_generation,
        energy_out_of_battery,
        energy_into_battery_from_grid,
        battery_state_of_charge,
    })?;

    Ok(())
}

/// Check (almost an assert) whether the shower flow rate is not less than the minimum allowed.
fn check_shower_flowrate(input: &InputForProcessing) -> anyhow::Result<()> {
    let min_flowrate = 8.0;

    for (name, flowrate) in input.shower_flowrates() {
        if flowrate < min_flowrate {
            // only currently known shower name that can have a flowrate is "mixer"
            bail!("Invalid flow rate: {flowrate} l/s in shower with name '{name}'");
        }
    }

    Ok(())
}

pub(super) fn create_hot_water_use_pattern(
    input: &mut InputForProcessing,
    number_of_occupants: f64,
    cold_water_feed_temps: &[f64],
) -> anyhow::Result<()> {
    check_shower_flowrate(input)?;

    // temperature of mixed hot water for event
    let event_temperature_showers = 41.0;
    let event_temperature_bath = 41.0;
    let event_temperature_others = 41.0;

    let mean_feedtemp =
        cold_water_feed_temps.iter().sum::<f64>() / cold_water_feed_temps.len() as f64;
    let _mean_delta_t = HW_TEMPERATURE - mean_feedtemp;

    let _annual_hw_events: Vec<()> = vec![];
    let _annual_hw_events_energy: Vec<()> = vec![];
    let startmod = 0;

    // SAP 2012 relation
    // vol_daily_average = (25 * N_occupants) + 36

    // new relation based on Boiler Manufacturer data and EST surveys
    // reduced by 30% to account for pipework losses present in the source data
    let mut vol_hw_daily_average = 0.70 * 60.3 * number_of_occupants.powf(0.71);

    // The hot water data set only included hot water use via the central hot water system
    // Electric showers are common in the UK, sometimes in addition to a central shower.
    // It is therefore very likely more showers were taken than are recorded in our main dataset.
    // To attempt to correct for this additional shower events (and their equivalent volume)
    // need to be added for use in generating the correct list of water use events.
    // It was assumed that 30% of the homes had an additional electric shower and these were
    // used half as often as showers from the central water heating system (due to lower flow).
    // This would mean that about 15% of showers taken were missing from the data.
    // The proportion of total hot water volume due to with showers in the original sample
    // was 60.685%. Increasing this by 15%, then re-adding it to the non-shower total gives
    // 109.10%. So we need to multiply the hot water use by 1.0910 to correct for the missing showers.
    // (Note that this is only being used to generate the correct events list so does not assume
    // the dwelling being modelled actually has an electric shower, or a central shower. Allocation
    // of events to the actual showers types present in the home is done later.)
    let prop_with_elec_shower = 0.3; // 30% of homes had an additional electric shower
    let elec_shower_use_prop_of_main = 0.5; // they are used half as often as the main shower
    let correction_for_missing_elec_showers =
        1. + prop_with_elec_shower * elec_shower_use_prop_of_main; // 1.15
    let original_prop_hot_water_showers = 0.60685; // from original data set
    let uplifted_prop_hot_water_showers =
        original_prop_hot_water_showers * correction_for_missing_elec_showers;
    let elec_shower_correction_factor =
        1. - original_prop_hot_water_showers + uplifted_prop_hot_water_showers;
    vol_hw_daily_average *= elec_shower_correction_factor;

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
        number_of_occupants,
        input,
        fhw,
        event_temperature_others,
        HW_TEMPERATURE,
        cold_water_feed_temps,
        part_g_bonus,
    )?;

    // now create lists of events
    // Shower events should be evenly spread across all showers in dwelling
    // and so on for baths etc
    let mut hourly_events: Vec<Vec<HourlyHotWaterEvent>> =
        std::iter::repeat_with(Vec::new).take(8760).collect();
    for event in &ref_event_list {
        // assign HW usage events to end users and work out their durations
        // note that if there are no baths in the dwelling "bath" events are
        // assigned to showers, and vice versa
        let drawoff = if event.event_type.is_shower_type() {
            hw_event_aa.get_shower()
        } else if event.event_type.is_bath_type() {
            hw_event_aa.get_bath()
        } else {
            hw_event_aa.get_other()
        };
        let duration = drawoff.call_duration_fn(*event);

        let event_start = event.time;
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
                volume: if event.event_type.is_bath_type() {
                    // if the end user the event is being assigned to has a defined flowrate
                    // we are able to supply a volume
                    input
                        .flowrate_for_bath_field(&drawoff.name)
                        .map(|flowrate| duration * flowrate)
                } else {
                    None
                },
                temperature: if event.event_type.is_shower_type() {
                    event_temperature_showers
                } else if event.event_type.is_bath_type() {
                    event_temperature_bath
                } else {
                    event_temperature_others
                },
            },
        );
    }

    Ok(())
}

fn window_treatment(input: &mut InputForProcessing) -> anyhow::Result<()> {
    let simtime = simtime();
    let extcond = create_external_conditions(
        input.external_conditions().as_ref().clone(),
        &simtime.iter(),
    )?;
    let mut curtain_opening_sched_manual: Vec<Option<bool>> = Default::default();
    let mut curtain_opening_sched_auto: Vec<bool> = Default::default();
    let mut blinds_closing_irrad_manual: Vec<Option<f64>> = Default::default();
    let mut blinds_opening_irrad_manual: Vec<Option<f64>> = Default::default();

    for t_it in simtime.iter() {
        let hour_of_day = t_it.hour_of_day() as usize;
        // TODO (from Python) Are these waking hours correct? Check consistency with other parts of calculation
        let waking_hours = (OCCUPANT_WAKING_HR..OCCUPANT_SLEEPING_HR).contains(&hour_of_day);
        let sun_above_horizon = extcond.sun_above_horizon(t_it);

        curtain_opening_sched_manual.push(if waking_hours && sun_above_horizon {
            Some(true) // Open during waking hours after sunrise
        } else if waking_hours && !sun_above_horizon {
            Some(false) // Close during waking hours after sunset
        } else {
            None // Do not adjust outside waking hours
        });
        curtain_opening_sched_auto.push(sun_above_horizon);
        blinds_closing_irrad_manual.push(if waking_hours { Some(300.) } else { None });
        blinds_opening_irrad_manual.push(if waking_hours { Some(200.) } else { None });
    }

    input.add_control(
        "_curtains_open_manual",
        json!({
            "type": "OnOffTimeControl",
            "allow_null": true,
            "start_day": 0,
            "time_series_step": SIMTIME_STEP,
            "schedule": {
                "main": curtain_opening_sched_manual,
            }
        }),
    )?;

    input.add_control(
        "_curtains_open_auto",
        json!({
            "type": "OnOffTimeControl",
            "start_day": 0,
            "time_series_step": SIMTIME_STEP,
            "schedule": {
                "main": curtain_opening_sched_auto,
            }
        }),
    )?;

    input.add_control(
        "_blinds_closing_irrad_manual",
        json!({
            "type": "SetpointTimeControl",
            "start_day": 0,
            "time_series_step": SIMTIME_STEP,
            "schedule": {
                "main": blinds_closing_irrad_manual,
            }
        }),
    )?;

    input.add_control(
        "_blinds_closing_irrad_auto",
        json!({
            "type": "SetpointTimeControl",
            "start_day": 0,
            "time_series_step": 1.,
            "schedule": {
                "main": [{"repeat": SIMTIME_END as usize, "value": 200.}],
            }
        }),
    )?;

    input.add_control(
        "_blinds_opening_irrad_manual",
        json!({
            "type": "SetpointTimeControl",
            "start_day": 0,
            "time_series_step": SIMTIME_STEP,
            "schedule": {
                "main": blinds_opening_irrad_manual,
            }
        }),
    )?;

    input.add_control(
        "_blinds_opening_irrad_auto",
        json!({
            "type": "SetpointTimeControl",
            "start_day": 0,
            "time_series_step": 1.,
            "schedule": {
                "main": [{"repeat": SIMTIME_END as usize, "value": 200.}],
            }
        }),
    )?;

    let transparent_building_elements = input.all_transparent_building_elements_mut();

    for building_element in transparent_building_elements.iter() {
        if let Some(window_treatments) = building_element.treatment() {
            for mut treatment in window_treatments {
                treatment.set_is_open(false);

                match treatment.treatment_type {
                    WindowTreatmentType::Curtains => {
                        if treatment.controls.is_manual() {
                            treatment.set_open_control("_curtains_open_manual");
                        } else {
                            treatment.set_open_control("_curtains_open_auto");
                        }
                    }
                    // blinds are opened and closed in response to solar irradiance incident upon them
                    WindowTreatmentType::Blinds => {
                        if treatment.controls.is_manual() {
                            // manual control - Table B.24 in BS EN ISO 52016-1:2017.
                            treatment
                                .set_closing_irradiance_control("_blinds_closing_irrad_manual");
                            treatment
                                .set_opening_irradiance_control("_blinds_opening_irrad_manual");
                        } else {
                            // automatic control - Table B.24 in BS EN ISO 52016-1:2017.
                            treatment.set_closing_irradiance_control("_blinds_closing_irrad_auto");
                            treatment.set_opening_irradiance_control("_blinds_opening_irrad_auto");
                            treatment.set_opening_delay_hrs(2.);
                        }
                    }
                }
            }
        }
    }

    Ok(())
}

pub(super) fn create_window_opening_schedule(input: &mut InputForProcessing) -> anyhow::Result<()> {
    let window_opening_setpoint = 22.0;

    input.add_control(
        "_window_opening_adjust",
        json!({
            "type": "SetpointTimeControl",
            "start_day": 0,
            "time_series_step": 1.0,
            "schedule": {
                "main": [{"repeat": SIMTIME_END as usize, "value": window_opening_setpoint}],
            }
        }),
    )?;
    input.set_window_adjust_control_for_infiltration_ventilation("_window_opening_adjust");

    input.add_control(
        "_window_opening_openablealways",
        json!({
            "type": "OnOffTimeControl",
            "start_day": 0,
            "time_series_step": 1.0,
            "schedule": {
                "main": [{"repeat": SIMTIME_END as usize, "value": true}]
            }
        }),
    )?;

    input.add_control(
        "_window_opening_closedsleeping",
        json!({
            "type": "OnOffTimeControl",
            "start_day": 0,
            "time_series_step": 1.0,
            "schedule": {
                "main": [{"repeat": 365, "value": "day"}],
                "day": [
                    {"repeat": OCCUPANT_WAKING_HR, "value": false},
                    {"repeat": OCCUPANT_SLEEPING_HR - OCCUPANT_WAKING_HR, "value": true},
                    {"repeat": 24 - OCCUPANT_SLEEPING_HR, "value": false},
                ]
            }
        }),
    )?;

    let noise_nuisance = input.infiltration_ventilation_is_noise_nuisance();

    for transparent_building_element in input.all_transparent_building_elements_mut() {
        let element_is_security_risk = transparent_building_element.is_security_risk();
        transparent_building_element.set_window_openable_control(
            if noise_nuisance || element_is_security_risk {
                "_window_opening_closedsleeping"
            } else {
                "_window_opening_openablealways"
            },
        );
    }

    Ok(())
}

/// Calculate effective air change rate accoring to according to Part F 1.24 a
pub(crate) fn minimum_air_change_rate(
    _input: &InputForProcessing,
    total_floor_area: f64,
    total_volume: f64,
    bedroom_number: usize,
) -> f64 {
    // minimum ventilation rates method B
    let min_ventilation_rates_b = [19, 25, 31, 37, 43];

    // Calculate minimum whole dwelling ventilation rate l/s method A
    let min_ventilation_rate_a = total_floor_area * 0.3;

    // Calculate minimum whole dwelling ventilation rate l/s method B
    let min_ventilation_rate_b = if bedroom_number <= 5 {
        min_ventilation_rates_b[bedroom_number - 1]
    } else {
        min_ventilation_rates_b.last().unwrap() + (bedroom_number - 5) * 6
    };

    // Calculate air change rate ACH
    let highest_min_ventilation_rate =
        f64::max(min_ventilation_rate_a, min_ventilation_rate_b as f64);

    highest_min_ventilation_rate / total_volume * SECONDS_PER_HOUR as f64
        / LITRES_PER_CUBIC_METRE as f64
}

/// Set min and max vent opening thresholds
fn create_vent_opening_schedule(input: &mut InputForProcessing) -> anyhow::Result<()> {
    let tfa = calc_tfa(input);
    let number_of_bedrooms = input
        .number_of_bedrooms()
        .ok_or_else(|| anyhow!("Expected number of bedrooms to be indicated."))?;
    let total_volume = input.total_zone_volume();

    let vent_adjust_min_ach = minimum_air_change_rate(input, tfa, total_volume, number_of_bedrooms);
    let vent_adjust_max_ach = 2.;

    input.add_control(
        "_vent_adjust_min_ach",
        json!({
            "type": "SetpointTimeControl",
            "start_day": 0,
            "time_series_step": 1.0,
            "schedule": {
                "main": [{"repeat": (SIMTIME_END - SIMTIME_START) as usize, "value": vent_adjust_min_ach}],
            }
        }),
    )?;
    input.set_vent_adjust_min_control_for_infiltration_ventilation("_vent_adjust_min_ach");

    input.add_control(
        "_vent_adjust_max_ach",
        json!({
            "type": "SetpointTimeControl",
            "start_day": 0,
            "time_series_step": 1.0,
            "schedule": {
                "main": [{"repeat": (SIMTIME_END - SIMTIME_START) as usize, "value": vent_adjust_max_ach}],
            }
        }),
    )?;
    input.set_vent_adjust_max_control_for_infiltration_ventilation("_vent_adjust_max_ach");

    Ok(())
}

fn create_mev_pattern(input: &mut InputForProcessing) -> anyhow::Result<()> {
    // intermittent extract fans are assumed to turn on whenever cooking, bath or shower events occur

    let shower_and_bath_events = input.water_heating_events_of_types(&[
        WaterHeatingEventType::Shower,
        WaterHeatingEventType::Bath,
    ]);
    let appliance_gains_events = input.appliance_gains_events();

    let mech_vents = input.keyed_mechanical_ventilations_for_processing();
    let mut intermittent_mev: IndexMap<String, Vec<f64>> = mech_vents
        .iter()
        .filter(|(_, vent)| vent.vent_type() == VentType::IntermittentMev)
        .fold(IndexMap::from([]), |mut acc, (vent, _)| {
            acc.insert(
                vent.to_owned(),
                vec![0.; ((SIMTIME_END - SIMTIME_START) / SIMTIME_STEP).ceil() as usize],
            );
            acc
        });

    let mev_names = intermittent_mev.keys().cloned().collect::<Vec<_>>();
    if mev_names.is_empty() {
        return Ok(());
    }

    let mut cycle_mev = CycleMev::new(mev_names.iter().map(String::as_str).collect());

    for event in shower_and_bath_events {
        let mev_name = cycle_mev.mev();
        let idx = (event.start / SIMTIME_STEP).floor() as usize;
        let tsfrac = event.duration.ok_or_else(|| {
            anyhow!("Water heating event was expected to have a defined duration in FHS transform.")
        })? / (MINUTES_PER_HOUR as f64 * SIMTIME_STEP);
        // add fraction of the timestep for which appliance is turned on
        // to the fraction of the timestep for which the fan is turned on,
        // and cap that fraction at 1.
        let mut integralx: f64 = Default::default();
        let start_offset = event.start / SIMTIME_STEP - idx as f64;
        while integralx < tsfrac {
            let segment = (start_offset.ceil() - start_offset).min(tsfrac - integralx);
            let step_idx = (idx + (start_offset + integralx).floor() as usize)
                % intermittent_mev[mev_name].len();
            intermittent_mev[mev_name][step_idx] =
                (intermittent_mev[mev_name][step_idx] + segment).min(1.);
            integralx += segment;
        }
    }

    // these names are the same as those already defined in create_appliance_gains
    // NB. have reported possible bug here https://dev.azure.com/BreGroup/SAP%2011/_workitems/edit/45690 as these names don't seem to match with known
    // TODO (from Python) - define them at top level of wrapper
    // kettles and microwaves are assumed not to activate the extract fan
    for cook_enduse in ["Oven", "Hobs"] {
        if let Some(events) = appliance_gains_events.get(cook_enduse) {
            for event in events {
                let mev_name = cycle_mev.mev();
                let idx = (event.start / SIMTIME_STEP).floor() as usize;
                let tsfrac = event.duration / (MINUTES_PER_HOUR as f64 * SIMTIME_STEP);
                // add fraction of the timestep for which appliance is turned on
                // to the fraction of the timestep for which the fan is turned on,
                // and cap that fraction at 1.
                let mut integralx: f64 = Default::default();
                let start_offset = event.start / SIMTIME_STEP - idx as f64;
                while integralx < tsfrac {
                    let segment = (start_offset.ceil() - start_offset).min(tsfrac - integralx);
                    let step_idx = (idx + (start_offset + integralx).floor() as usize)
                        % intermittent_mev[mev_name].len();
                    intermittent_mev[mev_name][step_idx] =
                        (intermittent_mev[mev_name][step_idx] + segment).min(1.);
                    integralx += segment;
                }
            }
        }
    }

    let control_names: HashMap<String, String> = intermittent_mev
        .keys()
        .map(|name| (name.clone(), format!("_intermittent_MEV_control: {}", name)))
        .collect();

    for vent in intermittent_mev.keys() {
        let control_name = &control_names[vent];
        mech_vents.get_mut(vent).unwrap().set_control(control_name);
    }

    // loop through again as can't write to two different mutable refs based on input in one loop
    for vent in intermittent_mev.keys() {
        let control_name = &control_names[vent];
        input.add_control(
            control_name,
            json!({
                "type": "SetpointTimeControl",
                "start_day": 0,
                "time_series_step": SIMTIME_STEP,
                "schedule": {
                    "main": intermittent_mev[vent]
                }
            }),
        )?;
    }

    Ok(())
}

// if there are multiple extract fans they are cycled sequentially
// in order that they all be used an approximately equal amount,
// so a different extract fan could be activated by the same shower,
// and likewise the same extract fan could be activated by cooking as by a shower
struct CycleMev<'a> {
    names: Vec<&'a str>,
    cycle_count: usize,
}

impl<'a> CycleMev<'a> {
    fn new(names: Vec<&'a str>) -> Self {
        Self {
            names,
            cycle_count: Default::default(),
        }
    }

    fn mev(&mut self) -> &str {
        let res = self.names[self.cycle_count];
        self.cycle_count = (self.cycle_count + 1) % self.names.len();
        res
    }
}

fn calc_sfp_mech_vent(input: &mut InputForProcessing) -> anyhow::Result<()> {
    for mech_vents_data in input.mechanical_ventilations_for_processing().iter_mut() {
        match mech_vents_data.vent_type() {
            crate::input::VentType::CentralisedContinuousMev | crate::input::VentType::Mvhr => {
                let measured_fan_power = mech_vents_data.measured_fan_power().ok_or_else(|| anyhow!("Measured fan power was not given for a mechanical ventilation that expected one to be present."))?;
                let measured_air_flow_rate = mech_vents_data.measured_air_flow_rate().ok_or_else(|| anyhow!("Measured air flow rate was not given for a mechanical ventilation that expected one to be present."))?;
                // Specific fan power is total measured electrical power in Watts divided by air flow rate
                let measured_sfp = measured_fan_power / measured_air_flow_rate; // in W/l/s
                mech_vents_data.set_sfp(measured_sfp);
            }
            crate::input::VentType::IntermittentMev
            | crate::input::VentType::DecentralisedContinuousMev => {
                continue;
            }
            crate::input::VentType::Piv => {
                bail!("Mechanical ventilation type of PIV not recognised")
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
                    match input.space_cool_system_for_zone(zone_key)? {
                        SystemReference::Single(space_cool_system) => {
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
                            if let Some(advanced_start) =
                                input.advanced_start_for_space_cool_system(&space_cool_system)?
                            {
                                match living_room_control {
                                    Value::Object(ref mut control_map) => {
                                        control_map.insert(
                                            "advanced_start".into(),
                                            Value::Number(Number::from_f64(advanced_start).unwrap()),
                                        );
                                    }
                                    _ => unreachable!(),
                                }
                            }
                            input.add_control("Cooling_LivingRoom", living_room_control)?;
                        }
                        SystemReference::Multiple(_) => bail!("Multiple space heat systems references in zones not currently supported for FHS inputs"),
                        SystemReference::None(_) => {}
                    }
                }
                SpaceHeatControlType::RestOfDwelling => {
                    match input.space_cool_system_for_zone(zone_key)? {
                        SystemReference::Single(space_cool_system) => {
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

                            if let Some(advanced_start) =
                                input.advanced_start_for_space_cool_system(&space_cool_system)?
                            {
                                match rest_of_dwelling_control {
                                    Value::Object(ref mut control_map) => {
                                        control_map.insert(
                                            "advanced_start".into(),
                                            Value::Number(Number::from_f64(advanced_start).unwrap()),
                                        );
                                    }
                                    _ => unreachable!(),
                                }
                            }
                            input.add_control("Cooling_RestOfDwelling", rest_of_dwelling_control)?;
                        }
                        SystemReference::Multiple(_) => bail!("Multiple space cool systems references in zones not currently supported for FHS inputs"),
                        SystemReference::None(_) => {}
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

pub(super) fn create_cold_water_feed_temps(
    input: &mut InputForProcessing,
) -> anyhow::Result<Vec<f64>> {
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
        (t24m_header_tank, ColdWaterSourceType::HeaderTank)
    } else {
        (t24m_mains, ColdWaterSourceType::MainsWater)
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

fn daylight_factor(input: &InputForProcessing, total_floor_area: f64) -> anyhow::Result<Vec<f64>> {
    let mut total_area = vec![0.; simtime().total_steps()];

    let data: Vec<Vec<f64>> = input
        .all_building_elements()
        .iter()
        .filter_map(|el| match el {
            crate::input::BuildingElement::Transparent {
                orientation,
                g_value,
                frame_area_fraction,
                base_height,
                height,
                width,
                shading,
                ..
            } => {
                let ff = *frame_area_fraction;
                let g_val = *g_value;
                let width = *width;
                let height = *height;
                let base_height = *base_height;
                let orientation = *orientation;
                let w_area = width * height;
                // retrieve half-hourly shading factor
                let direct_result =
                    shading_factor(input, base_height, height, width, orientation, shading);

                let area = 0.9 * w_area * (1. - ff) * g_val;

                match direct_result {
                    Ok(direct) => Some(Ok(direct.iter().map(|factor| factor * area).collect())),
                    Err(err) => Some(Err(err)),
                }
            }
            _ => None,
        })
        .collect::<anyhow::Result<Vec<_>>>()?;

    for idx in data {
        for (t, gl) in idx.into_iter().enumerate() {
            total_area[t] += gl;
        }
    }

    // calculate Gl for each half hourly timestep
    Ok((0..(total_area.len()))
        .map(|i| {
            let gl = total_area[i] / total_floor_area;

            if gl > 0.095 {
                0.96
            } else {
                52.2 * gl.powi(2) - 9.94 * gl + 1.433
            }
        })
        .collect())
}

fn shading_factor(
    input: &InputForProcessing,
    base_height: f64,
    height: f64,
    width: f64,
    orientation: f64,
    shading: &[WindowShadingObject],
) -> anyhow::Result<Vec<f64>> {
    // there is code in the upstream Python to convert orientations from -180 to +180 (anticlockwise) to 0-360 (clockwise)
    // but the Rust input code has already implicitly performed this conversion on the way in, so we don't need to do it here

    let time = simtime();

    let input_external_conditions = input.external_conditions();

    let dir_beam_conversion = input_external_conditions
        .direct_beam_conversion_needed
        .is_some_and(|x| x);

    let conditions = ExternalConditions::new(
        &time.iter(),
        input_external_conditions
            .air_temperatures
            .as_ref()
            .ok_or_else(|| anyhow!("Air temps were expected in input and not provided."))?
            .to_vec(),
        input_external_conditions
            .wind_speeds
            .as_ref()
            .ok_or_else(|| anyhow!("Wind speeds were expected in input and not provided."))?
            .to_vec(),
        input_external_conditions
            .wind_directions
            .as_ref()
            .ok_or_else(|| anyhow!("Wind directions were expected in input and not provided."))?
            .to_vec(),
        input_external_conditions
            .diffuse_horizontal_radiation
            .as_ref()
            .ok_or_else(|| {
                anyhow!("Diffuse horizontal radiations were expected in input and not provided.")
            })?
            .to_vec(),
        input_external_conditions
            .direct_beam_radiation
            .as_ref()
            .ok_or_else(|| {
                anyhow!("Direct beam radiations were expected in input and not provided.")
            })?
            .to_vec(),
        input_external_conditions
            .solar_reflectivity_of_ground
            .as_ref()
            .ok_or_else(|| {
                anyhow!("Solar reflectivity of ground was expected in input and not provided.")
            })?
            .to_vec(),
        input_external_conditions
            .latitude
            .ok_or_else(|| anyhow!("Latitude was expected in input and not provided."))?,
        input_external_conditions
            .longitude
            .ok_or_else(|| anyhow!("Longitude was expected in input and not provided."))?,
        0,
        0,
        Some(365),
        1.,
        None,
        None,
        false,
        dir_beam_conversion,
        input_external_conditions.shading_segments.to_vec(),
    );

    time.iter()
        .map(|t_it| {
            conditions.direct_shading_reduction_factor(
                base_height,
                height,
                width,
                orientation,
                Some(shading),
                t_it,
            )
        })
        .collect()
}

fn top_up_lighting(input: &InputForProcessing, l_req: f64) -> anyhow::Result<f64> {
    if !input.all_zones_have_bulbs() {
        bail!("At least one zone has lighting that does not have bulbs defined.");
    }

    let capacity_tot = input
        .light_bulbs_for_all_zones()
        .iter()
        .map(ZoneLightingBulbs::capacity)
        .sum::<f64>();

    let tfa = calc_tfa(input);
    let capacity_ref = 330. * tfa;

    let l_prov = l_req * (capacity_tot / capacity_ref);

    let l_topup = if l_prov < (l_req / 3.) {
        (l_req / 3.) - l_prov
    } else {
        0.
    };

    Ok(l_topup)
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

const AVERAGE_MONTHLY_APPLIANCES_HOURLY_PROFILES: [[f64; 24]; 12] = [
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
    // remove following lint escape once test is made good
    #[allow(clippy::assertions_on_constants)]
    #[allow(clippy::nonminimal_bool)]
    #[rstest]
    fn test_check_invalid_shower_flowrate() {
        assert!(!false);
    }

    #[ignore = "useless test reported up to BRE"]
    // remove following lint escape once test is made good
    #[allow(clippy::assertions_on_constants)]
    #[rstest]
    fn test_check_valid_shower_flowrate() {
        assert!(true);
    }

    #[ignore = "useless test reported up to BRE"]
    // remove following lint escape once test is made good
    #[allow(clippy::assertions_on_constants)]
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
