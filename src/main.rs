extern crate hem;

use clap::{Args, Parser};
use hem::output::FileOutput;
use hem::read_weather_file::{weather_data_to_vec, ExternalConditions};
use hem::run_project;
use std::ffi::OsStr;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

#[derive(Parser, Default, Debug)]
#[clap(author, version, about, long_about = None)]
struct SapArgs {
    input_file: String,
    #[command(flatten)]
    weather_file: WeatherFileType,
    #[arg(
        long,
        short,
        default_value_t = false,
        help = "Run preprocessing step only"
    )]
    preprocess_only: bool,
    #[command(flatten)]
    wrapper_choice: WrapperChoice,
    #[clap(
        long,
        default_value_t = false,
        help = "Output heat balance for each zone"
    )]
    heat_balance: bool,
}

#[derive(Args, Clone, Default, Debug)]
#[group(required = false, multiple = false)]
struct WrapperChoice {
    #[arg(long, help = "Use Future Homes Standard calculation assumptions")]
    future_homes_standard: bool,
    #[arg(
        long = "future-homes-standard-FEE",
        help = "Use Future Homes Standard Fabric Energy Efficiency assumptions"
    )]
    future_homes_standard_fee: bool,
    #[arg(
        long = "future-homes-standard-notA",
        help = "Use Future Homes Standard calculation assumptions for notional option A"
    )]
    future_homes_standard_not_a: bool,
    #[arg(
        long = "future-homes-standard-notB",
        help = "Use Future Homes Standard calculation assumptions for notional option B"
    )]
    future_homes_standard_not_b: bool,
    #[arg(
        long = "future-homes-standard-FEE-notA",
        help = "Use Future Homes Standard Fabric Energy Efficiency assumptions for notional option A"
    )]
    future_homes_standard_fee_not_a: bool,
    #[arg(
        long = "future-homes-standard-FEE-notB",
        help = "Use Future Homes Standard Fabric Energy Efficiency assumptions for notional option B"
    )]
    future_homes_standard_fee_not_b: bool,
}

#[derive(Args, Clone, Default, Debug)]
#[group(required = false, multiple = false)]
struct WeatherFileType {
    #[arg(long, short, help = "Path to weather file in .epw format")]
    epw_file: Option<String>,
    #[arg(
        long = "CIBSE-weather-file",
        short,
        help = "Path to CIBSE weather file in .csv format"
    )]
    cibse_weather_file: Option<String>,
}

fn main() -> anyhow::Result<()> {
    // set up basic tracing
    let tracing_subscriber = tracing_subscriber::fmt::fmt()
        .with_max_level(tracing::Level::TRACE)
        .finish();
    tracing::subscriber::set_global_default(tracing_subscriber)
        .expect("setting tracing subscriber failed");

    let args = SapArgs::parse();

    let input_file = args.input_file.as_str();
    let input_file_ext = Path::new(input_file).extension().and_then(OsStr::to_str);
    let input_file_stem = match input_file_ext {
        Some(ext) => &input_file[..(input_file.len() - ext.len() - 1)],
        None => input_file,
    };
    let input_file_stem = PathBuf::from(input_file_stem);

    // let output_file_detailed = format!("{input_file_stem}_results.csv");
    // let _output_file_static = format!("{input_file_stem}_results_static.csv");
    // let _output_file_summary = format!("{input_file_stem}_results_summary.csv");

    let mut output_path = PathBuf::new();
    output_path.push(format!("{}__results", input_file_stem.to_str().unwrap()));
    fs::create_dir_all(&output_path)?;
    let input_file_name = input_file_stem.file_name().unwrap().to_str().unwrap();
    // following is rough initial mapping given existing fhs options
    let output_type = if args.wrapper_choice.future_homes_standard {
        "FHS"
    } else if args.wrapper_choice.future_homes_standard_fee {
        "FHS_FEE"
    } else {
        "core"
    };
    let file_output = FileOutput::new(
        output_path,
        format!("{input_file_name}__{output_type}__{{}}.csv"),
    );

    let external_conditions: Option<ExternalConditions> = match args.weather_file {
        WeatherFileType {
            epw_file: Some(file),
            cibse_weather_file: None,
        } => {
            let external_conditions_data = weather_data_to_vec(File::open(file)?);
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

    run_project(
        BufReader::new(File::open(Path::new(input_file))?),
        file_output,
        external_conditions,
        args.preprocess_only,
        args.wrapper_choice.future_homes_standard,
        args.wrapper_choice.future_homes_standard_fee,
        args.wrapper_choice.future_homes_standard_not_a,
        args.wrapper_choice.future_homes_standard_not_b,
        args.wrapper_choice.future_homes_standard_fee_not_a,
        args.wrapper_choice.future_homes_standard_fee_not_b,
        args.heat_balance,
        false, // TODO implement CLI arg
    )
}
