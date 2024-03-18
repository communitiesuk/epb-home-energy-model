extern crate hem;

mod read_weather_file;

use clap::{Args, Parser};
use hem::read_weather_file::{weather_data_to_vec, ExternalConditions};
use hem::run_project;

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

fn main() {
    let args = SapArgs::parse();

    println!("{:?}", args);

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

    run_project(
        args.input_file.as_str(),
        external_conditions,
        args.preprocess_only,
        false,
        false,
        false,
        false,
        args.heat_balance,
    )
    .unwrap_or_else(|err| println!("{:?}", err))
}
