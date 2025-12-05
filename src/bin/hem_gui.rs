use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::net::SocketAddr;
use std::path::{Path, PathBuf};

use axum::{
    response::Html,
    routing::{get, post},
    Json, Router,
};
use hem::output::FileOutput;
use hem::read_weather_file::{weather_data_to_vec, ExternalConditions};
use hem::{run_project, ProjectFlags};
use serde::{Deserialize, Serialize};
use tokio::net::TcpListener;
use tracing::{error, info};
use tracing_subscriber::FmtSubscriber;

#[derive(Deserialize)]
struct RunRequest {
    /// input JSON
    input_path: String,
    /// weather inputs
    weather_path: String,
    /// nice optional flags for heat balance and heating and cooling
    heat_balance: bool,
    detailed_output: bool,
}

#[derive(Serialize)]
struct RunResponse {
    ok: bool,
    message: String,
    // demo__results_gui
    output_folder: Option<String>,
    // absolute path example C:/Users/.../demo__results_gui
    output_folder_abs: Option<String>,
}

#[derive(Serialize)]
struct FileListResponse {
    files: Vec<String>,
}

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    // errror logging
    let subscriber = FmtSubscriber::builder()
        .with_max_level(tracing::Level::INFO)
        .finish();
    tracing::subscriber::set_global_default(subscriber)
        .expect("setting default subscriber failed");

    // flow
    //  GUI HTML
    //  /examples - list example JSONs
    //  /weather-files - list EPW files
    //  run it
    let app = Router::new()
        .route("/", get(index_handler))
        .route("/api/run", post(run_handler))
        .route("/api/examples", get(list_examples))
        .route("/api/weather-files", get(list_weather_files));

    let addr = SocketAddr::from(([127, 0, 0, 1], 3000));
    info!("Starting HEM GUI on http://{addr}");

    let listener = TcpListener::bind(addr).await?;
    axum::serve(listener, app.into_make_service()).await?;

    Ok(())
}

async fn index_handler() -> Html<&'static str> {
    // serve the gui
    Html(include_str!("../../gui/index.html"))
}

async fn run_handler(Json(payload): Json<RunRequest>) -> Json<RunResponse> {
    match run_hem_from_paths(&payload) {
        Ok((output_folder, output_folder_abs)) => Json(RunResponse {
            ok: true,
            message: "Simulation completed successfully".to_string(),
            output_folder: Some(output_folder),
            output_folder_abs: Some(output_folder_abs),
        }),
        Err(e) => {
            error!("Error running HEM: {e:?}");
            Json(RunResponse {
                ok: false,
                message: format!("Error: {e}"),
                output_folder: None,
                output_folder_abs: None,
            })
        }
    }
}

/// Calls the HEM engine properly
fn run_hem_from_paths(req: &RunRequest) -> anyhow::Result<(String, String)> {
    // path resolver
    let input_path = PathBuf::from(&req.input_path);
    if !input_path.exists() {
        anyhow::bail!("Input file not found: {}", input_path.display());
    }

    // same output logic as main.rs but for the gui
    let input_stem = input_path
        .file_stem()
        .ok_or_else(|| anyhow::anyhow!("Could not determine input file name"))?
        .to_string_lossy()
        .to_string();

    let output_folder = format!("{input_stem}__results_gui");
    let output_path = PathBuf::from(&output_folder);
    fs::create_dir_all(&output_path)?;

    // same logic as main.rs again
    let file_output = FileOutput::new(
        output_path.clone(),
        format!("{input_stem}__gui__{{}}.{{}}"),
    );

    // read weather
    let weather_path = PathBuf::from(&req.weather_path);
    if !weather_path.exists() {
        anyhow::bail!("Weather file not found: {}", weather_path.display());
    }

    let external_conditions_data: Option<ExternalConditions> =
        match weather_data_to_vec(File::open(&weather_path)?) {
            Ok(data) => Some(data),
            Err(_) => anyhow::bail!("Could not parse the weather file!"),
        };

    // nice boolean options
    let mut flags = ProjectFlags::empty();
    if req.heat_balance {
        flags.insert(ProjectFlags::HEAT_BALANCE);
    }
    if req.detailed_output {
        flags.insert(ProjectFlags::DETAILED_OUTPUT_HEATING_COOLING);
    }

    // run it!
    let response = run_project(
        BufReader::new(File::open(&input_path)?),
        file_output,
        external_conditions_data,
        None,
        &flags,
    )?;

    if response.is_some() {
        info!("Received JSON response from HEM engine");
    }

    // file path for convenience
    let abs = fs::canonicalize(&output_path)?;
    let mut abs_str = abs.to_string_lossy().to_string();

    abs_str = abs_str.replace("\\", "/");


    Ok((output_folder, abs_str))
}

/// We are going to list JSON input files under examples/input/core
async fn list_examples() -> Json<FileListResponse> {
    let mut files = Vec::new();
    let dir = Path::new("examples").join("input").join("core");

    if let Ok(entries) = fs::read_dir(&dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path
                .extension()
                .and_then(|s| s.to_str())
                .map(|ext| ext.eq_ignore_ascii_case("json"))
                .unwrap_or(false)
            {
                if let Some(name) = path.file_name().and_then(|s| s.to_str()) {
                    // Return the path in the format users can pass straight back into the API
                    files.push(format!("examples/input/core/{}", name));
                }
            }
        }
    }

    Json(FileListResponse { files })
}

/// list the weather files EPW ones
async fn list_weather_files() -> Json<FileListResponse> {
    let mut files = Vec::new();
    let dir = Path::new("weather");

    if let Ok(entries) = fs::read_dir(&dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path
                .extension()
                .and_then(|s| s.to_str())
                .map(|ext| ext.eq_ignore_ascii_case("epw"))
                .unwrap_or(false)
            {
                if let Some(name) = path.file_name().and_then(|s| s.to_str()) {
                    files.push(format!("weather/{}", name));
                }
            }
        }
    }

    Json(FileListResponse { files })
}
