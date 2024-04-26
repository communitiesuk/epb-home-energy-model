extern crate hem;

use clap::{Args, Parser};
use csv::WriterBuilder;
use hem::corpus::{HotWaterResultMap, KeyString};
use hem::read_weather_file::{weather_data_to_vec, ExternalConditions};
use hem::{run_project, UNITS_MAP};
use indexmap::IndexMap;
use std::borrow::Cow;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

#[derive(Parser, Default, Debug)]
#[clap(author, version, about, long_about = None)]
struct SapArgs {
    input_file: String,
    #[command(flatten)]
    weather_file: WeatherFileType,
    #[arg(long, short, default_value_t = false)]
    preprocess_only: bool,
    #[command(flatten)]
    wrapper_choice: WrapperChoice,
    #[clap(long, default_value_t = false)]
    heat_balance: bool,
}

#[derive(Args, Clone, Default, Debug)]
#[group(required = false, multiple = false)]
struct WrapperChoice {
    #[arg(long)]
    future_homes_standard: bool,
    #[arg(long)]
    future_homes_standard_fee: bool,
}

#[derive(Args, Clone, Default, Debug)]
#[group(required = false, multiple = false)]
struct WeatherFileType {
    #[arg(long, short)]
    epw_file: Option<String>,
    #[arg(long, short)]
    cibse_weather_file: Option<String>,
}

fn main() -> anyhow::Result<()> {
    let args = SapArgs::parse();

    let input_file = args.input_file.as_str();
    let input_file_ext = Path::new(input_file).extension().and_then(OsStr::to_str);
    let input_file_stem = match input_file_ext {
        Some(ext) => &input_file[..(input_file.len() - ext.len() - 1)],
        None => input_file,
    };
    let output_file_detailed = format!("{input_file_stem}_results.csv");
    let _output_file_static = format!("{input_file_stem}_results_static.csv");
    let _output_file_summary = format!("{input_file_stem}_results_summary.csv");

    let external_conditions: Option<ExternalConditions> = match args.weather_file {
        WeatherFileType {
            epw_file: Some(file),
            cibse_weather_file: None,
        } => {
            let external_conditions_data = weather_data_to_vec(file.as_str());
            match external_conditions_data {
                Ok(data) => Some(data),
                Err(_) => panic!("Could not parse the weather file!"),
            }
        }
        WeatherFileType {
            epw_file: None,
            cibse_weather_file: Some(_),
        } => None,
        _ => None,
    };

    let (
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
        _heat_cop_dict,
        _cool_cop_dict,
        _dhw_cop_dict,
        ductwork_gains,
        _heat_balance_dict,
        _heat_source_wet_results_dict,
        _heat_source_wet_results_annual_dict,
    ) = run_project(
        BufReader::new(File::open(Path::new(input_file))?),
        external_conditions,
        args.preprocess_only,
        false,
        false,
        false,
        false,
        args.heat_balance,
    )?;

    write_core_output_file(OutputFileArgs {
        output_file: output_file_detailed,
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
        ductwork_gains,
    })
}

struct OutputFileArgs {
    output_file: String,
    timestep_array: Vec<f64>,
    results_totals: IndexMap<KeyString, Vec<f64>>,
    results_end_user: IndexMap<KeyString, IndexMap<String, Vec<f64>>>,
    energy_import: IndexMap<KeyString, Vec<f64>>,
    energy_export: IndexMap<KeyString, Vec<f64>>,
    energy_generated_consumed: IndexMap<KeyString, Vec<f64>>,
    energy_to_storage: IndexMap<KeyString, Vec<f64>>,
    energy_from_storage: IndexMap<KeyString, Vec<f64>>,
    energy_diverted: IndexMap<KeyString, Vec<f64>>,
    betafactor: IndexMap<KeyString, Vec<f64>>,
    zone_dict: IndexMap<KeyString, IndexMap<KeyString, Vec<f64>>>,
    zone_list: Vec<KeyString>,
    hc_system_dict: IndexMap<KeyString, IndexMap<KeyString, Vec<f64>>>,
    hot_water_dict: IndexMap<KeyString, HotWaterResultMap>,
    ductwork_gains: IndexMap<KeyString, Vec<f64>>,
}

fn write_core_output_file(args: OutputFileArgs) -> Result<(), anyhow::Error> {
    let OutputFileArgs {
        output_file,
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
        ductwork_gains,
    } = args;
    println!("writing out to {output_file}");
    let file = File::create(output_file)?;
    let writer = BufWriter::new(file);
    let mut writer = WriterBuilder::new().flexible(true).from_writer(writer);

    let mut headings: Vec<Cow<'static, str>> = vec!["Timestep".into()];
    let mut units_row = vec!["[count]"];
    for totals_key in results_totals.keys() {
        let totals_header = format!("{totals_key} total");
        headings.push(totals_header.into());
        units_row.push("[kWh]");
        for end_user_key in results_end_user[totals_key].keys() {
            headings.push((*end_user_key).clone().into());
            units_row.push("[kWh]");
        }
        headings.push(format!("{totals_key} import").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key} export").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key} generated and consumed").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key} beta factor").into());
        units_row.push("[ratio]");
        headings.push(format!("{totals_key} to storage").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key} from storage").into());
        units_row.push("[kWh]");
        headings.push(format!("{totals_key} diverted").into());
        units_row.push("[kWh]");

        for zone in &zone_list {
            for zone_outputs in zone_dict.keys() {
                let zone_headings = format!("{zone_outputs} {zone}");
                headings.push(zone_headings.into());
                if UNITS_MAP.contains_key(zone_outputs) {
                    units_row.push(UNITS_MAP[zone_outputs]);
                } else {
                    units_row.push("Unit not defined");
                }
            }
        }

        for system in hc_system_dict.keys() {
            for hc_name in hc_system_dict[system].keys() {
                let hc_system_headings = format!("{system} {hc_name}");
                headings.push(hc_system_headings.into());
                units_row.push("[kWh]");
            }
        }

        headings.push("Ductwork gains".into());
        units_row.push("[kWh]");

        // Write headings and units to output file
        writer.write_record(headings.iter().map(|heading| heading.as_ref()))?;
        writer.write_record(&units_row)?;

        for (t_idx, _timestep) in timestep_array.iter().enumerate() {
            let mut energy_use_row = vec![];
            let mut zone_row = vec![];
            let mut hc_system_row = vec![];
            let mut hw_system_row = vec![];
            let mut hw_system_row_energy = vec![];
            let mut hw_system_row_duration = vec![];
            let mut hw_system_row_events = vec![];
            let mut pw_losses_row = vec![];
            let mut ductwork_row = vec![];
            let energy_shortfall: Vec<f64> = vec![];
            for totals_key in results_totals.keys() {
                energy_use_row.push(results_totals[totals_key][t_idx]);
                for (end_user_key, _) in results_end_user[totals_key].iter().enumerate() {
                    energy_use_row.push(results_end_user[totals_key][end_user_key][t_idx]);
                }
                energy_use_row.push(energy_import[totals_key][t_idx]);
                energy_use_row.push(energy_export[totals_key][t_idx]);
                energy_use_row.push(energy_generated_consumed[totals_key][t_idx]);
                energy_use_row.push(betafactor[totals_key][t_idx]);
                energy_use_row.push(energy_to_storage[totals_key][t_idx]);
                energy_use_row.push(energy_from_storage[totals_key][t_idx]);
                energy_use_row.push(energy_diverted[totals_key][t_idx]);
            }

            // Loop over results separated by zone
            for zone in &zone_list {
                for zone_outputs in zone_dict.keys() {
                    zone_row.push(zone_dict[zone_outputs][zone][t_idx]);
                }
            }

            // Loop over heading and cooling system demand
            for system in hc_system_dict.keys() {
                for hc_name in hc_system_dict[system].keys() {
                    hc_system_row.push(hc_system_dict[system][hc_name][t_idx]);
                }
            }

            // Loop over hot water demand
            if let HotWaterResultMap::Float(map) = &hot_water_dict["Hot water demand"] {
                hw_system_row.push(map["demand"][t_idx]);
            }
            if let HotWaterResultMap::Float(map) = &hot_water_dict["Hot water energy demand"] {
                hw_system_row_energy.push(map["energy_demand"][t_idx]);
            }
            if let HotWaterResultMap::Float(map) = &hot_water_dict["Hot water duration"] {
                hw_system_row_duration.push(map["duration"][t_idx]);
            }
            if let HotWaterResultMap::Float(map) = &hot_water_dict["Pipework losses"] {
                pw_losses_row.push(map["pw_losses"][t_idx]);
            }
            if let HotWaterResultMap::Int(map) = &hot_water_dict["Hot Water Events"] {
                hw_system_row_events.push(map["no_events"][t_idx]);
            }
            ductwork_row.push(ductwork_gains["ductwork_gains"][t_idx]);

            // create row of outputs and write to output file
            let mut row: Vec<String> = vec![];
            row.append(&mut vec![t_idx.to_string()]);
            row.append(
                &mut energy_use_row
                    .into_iter()
                    .map(|val| val.to_string())
                    .collect(),
            );
            row.append(&mut zone_row.into_iter().map(|val| val.to_string()).collect());
            row.append(
                &mut hc_system_row
                    .into_iter()
                    .map(|val| val.to_string())
                    .collect(),
            );
            row.append(
                &mut hw_system_row
                    .into_iter()
                    .map(|val| val.to_string())
                    .collect(),
            );
            row.append(
                &mut hw_system_row_energy
                    .into_iter()
                    .map(|val| val.to_string())
                    .collect(),
            );
            row.append(
                &mut hw_system_row_duration
                    .into_iter()
                    .map(|val| val.to_string())
                    .collect(),
            );
            row.append(
                &mut hw_system_row_events
                    .into_iter()
                    .map(|val| val.to_string())
                    .collect(),
            );
            row.append(
                &mut pw_losses_row
                    .into_iter()
                    .map(|val| val.to_string())
                    .collect(),
            );
            row.append(
                &mut ductwork_row
                    .into_iter()
                    .map(|val| val.to_string())
                    .collect(),
            );
            row.append(
                &mut energy_shortfall
                    .into_iter()
                    .map(|val| val.to_string())
                    .collect(),
            );

            writer.write_record(&row)?;
        }
    }

    println!("flushing out CSV");
    writer.flush()?;

    Ok(())
}
