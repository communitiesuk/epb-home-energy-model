use clap::{Args, Parser};

use home_energy_model::output::FileOutput;
use home_energy_model::read_weather_file::{
    epw_weather_data_to_external_conditions, ExternalConditions,
};
use home_energy_model::{run_project, RunInput};
use std::ffi::OsStr;
use std::fs;
use std::fs::File;
use std::path::{Path, PathBuf};
use tracing_subscriber::fmt::format::FmtSpan;

#[derive(Parser, Default, Debug)]
#[clap(author, version, about, long_about = None)]
struct SapArgs {
    input_file: String,
    #[command(flatten)]
    weather_file: WeatherFileType,
    #[arg(long, short, help = "Path to tariff data file in .csv format")]
    tariff_file: Option<String>,
    #[clap(
        long,
        default_value_t = false,
        help = "Output heat balance for each zone"
    )]
    heat_balance: bool,
    #[clap(long, default_value_t = false, help = "Whether to log out spans")]
    log_spans: bool,
    #[clap(
        long,
        default_value_t = false,
        help = "Whether to output detailed information about heating and cooling"
    )]
    detailed_output_heating_cooling: bool,
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
    let args = SapArgs::parse();

    // set up basic tracing
    let tracing_subscriber = {
        let mut builder = tracing_subscriber::fmt::fmt().with_max_level(tracing::Level::TRACE);

        if args.log_spans {
            builder = builder.with_span_events(FmtSpan::CLOSE);
        }

        builder.finish()
    };
    tracing::subscriber::set_global_default(tracing_subscriber)
        .expect("setting tracing subscriber failed");

    let input_file = args.input_file.as_str();
    let input_file_ext = Path::new(input_file).extension().and_then(OsStr::to_str);
    let input_file_stem = match input_file_ext {
        Some(ext) => &input_file[..(input_file.len() - ext.len() - 1)],
        None => input_file,
    };
    let input_file_stem = PathBuf::from(input_file_stem);

    let mut output_path = PathBuf::new();
    output_path.push(format!("{}__results", input_file_stem.to_str().unwrap()));
    fs::create_dir_all(&output_path)?;
    let input_file_name = input_file_stem.file_name().unwrap().to_str().unwrap();
    let output_type = "core";
    let file_output = FileOutput::new(
        output_path,
        format!("{input_file_name}__{output_type}__{{}}.{{}}"),
    );

    let external_conditions: Option<ExternalConditions> = match args.weather_file {
        WeatherFileType {
            epw_file: Some(ref file),
            cibse_weather_file: None,
        } => {
            let external_conditions_data =
                epw_weather_data_to_external_conditions(File::open(file)?);
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
        RunInput::Read(Box::new(File::open(Path::new(input_file))?)),
        &file_output,
        external_conditions,
        args.tariff_file.as_ref().map(|f| f.as_str()),
        false,
        false,
    )?;

    Ok(())
}
