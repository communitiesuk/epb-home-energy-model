mod external_conditions;
mod read_weather_file;
mod simulation_time;

use crate::external_conditions::ExternalConditions;
use clap::{Args, Parser};

#[derive(Parser, Default, Debug)]
#[clap(author, version, about, long_about = None)]
struct SapArgs {
    input_file: Vec<String>,
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
        } => None,
        WeatherFileType {
            epw_file: None,
            cibse_weather_file: Some(file),
        } => None,
        _ => None,
    };

    println!("about to loop over the input files and run the calculation on each one!");
}
