use clap::{Args, Parser};
use hem::ProjectFlags;
use hem::output::FileOutput;
use hem::read_weather_file::{ExternalConditions, weather_data_to_vec};
use std::ffi::OsStr;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use tracing::debug;
use tracing_subscriber::fmt::format::FmtSpan;
use wrappers::run_wrappers;

#[derive(Parser, Default, Debug)]
#[clap(author, version, about, long_about = None)]
struct WrappersArgs {
    input_file: String,
    #[command(flatten)]
    weather_file: WeatherFileType,
    #[arg(long, short, help = "Path to tariff data file in .csv format")]
    tariff_file: Option<String>,
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
    #[clap(long, default_value_t = false, help = "Whether to log out spans")]
    log_spans: bool,
    #[clap(
        long,
        default_value_t = false,
        help = "Whether to output detailed information about heating and cooling"
    )]
    detailed_output_heating_cooling: bool,
}

#[derive(Args, Clone, Copy, Default, Debug)]
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
    #[arg(
        long = "fhs-compliance",
        help = "Run an FHS compliance calculation. This overrides all other FHS related flags"
    )]
    fhs_compliance: bool,
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
    let args = WrappersArgs::parse();

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
    // following is rough initial mapping given existing fhs options
    let output_type = output_type_from_wrapper_choice(&args.wrapper_choice);
    let file_output = FileOutput::new(
        output_path,
        format!("{input_file_name}__{output_type}__{{}}.{{}}"),
    );

    let external_conditions: Option<ExternalConditions> = match args.weather_file {
        WeatherFileType {
            epw_file: Some(ref file),
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

    let project_flags = (&args).into();

    let response = run_wrappers(
        BufReader::new(File::open(Path::new(input_file))?),
        &file_output,
        external_conditions,
        args.tariff_file.as_ref().map(|f| f.as_str()),
        &project_flags,
        args.preprocess_only,
    )?;

    if let Some(response) = response {
        debug!(
            "JSON response: {}",
            serde_json::to_string_pretty(&response)?
        );
    }

    Ok(())
}

fn output_type_from_wrapper_choice(wrapper_choice: &WrapperChoice) -> &str {
    if wrapper_choice.future_homes_standard {
        "FHS"
    } else if wrapper_choice.future_homes_standard_fee {
        "FHS_FEE"
    } else if wrapper_choice.future_homes_standard_not_a {
        "FHS_notA"
    } else if wrapper_choice.future_homes_standard_not_b {
        "FHS_notB"
    } else if wrapper_choice.future_homes_standard_fee_not_a {
        "FHS_FEE_notA"
    } else if wrapper_choice.future_homes_standard_fee_not_b {
        "FHS_FEE_notB"
    } else {
        "compliance"
    }
}

impl From<&WrappersArgs> for ProjectFlags {
    fn from(args: &WrappersArgs) -> Self {
        let mut flags = ProjectFlags::empty();
        if args.heat_balance {
            flags.insert(ProjectFlags::HEAT_BALANCE);
        }
        if args.detailed_output_heating_cooling {
            flags.insert(ProjectFlags::DETAILED_OUTPUT_HEATING_COOLING);
        }
        {
            let fhs = args.wrapper_choice;
            if fhs.future_homes_standard {
                flags.insert(ProjectFlags::FHS_ASSUMPTIONS);
            }
            if fhs.future_homes_standard_fee {
                flags.insert(ProjectFlags::FHS_FEE_ASSUMPTIONS);
            }
            if fhs.future_homes_standard_not_a {
                flags.insert(ProjectFlags::FHS_NOT_A_ASSUMPTIONS);
            }
            if fhs.future_homes_standard_not_b {
                flags.insert(ProjectFlags::FHS_NOT_B_ASSUMPTIONS);
            }
            if fhs.future_homes_standard_fee_not_a {
                flags.insert(ProjectFlags::FHS_FEE_NOT_A_ASSUMPTIONS)
            }
            if fhs.future_homes_standard_fee_not_b {
                flags.insert(ProjectFlags::FHS_FEE_NOT_B_ASSUMPTIONS)
            }
            if fhs.fhs_compliance {
                flags.insert(ProjectFlags::FHS_COMPLIANCE);
            }
        }

        flags
    }
}
