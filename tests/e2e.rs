use home_energy_model::output_writer::OutputWriter;
use home_energy_model::read_weather_file::cibse_weather_data_to_external_conditions;
use home_energy_model::{run_project_from_input_file, OutputFormat};
use indexmap::IndexMap;
use itertools::Itertools;
use parking_lot::{Mutex, RwLock};
use rayon::prelude::*;
use rstest::*;
use std::fmt::Display;
use std::fs::File;
use std::io::{BufReader, Cursor, Write};
use std::path::Path;
use std::str::from_utf8;
use std::sync::Arc;
use walkdir::{DirEntry, WalkDir};

#[fixture]
fn files() -> Vec<DirEntry> {
    WalkDir::new("./examples/input/core")
        .into_iter()
        .filter_map(Result::ok)
        .filter(|e| {
            !e.file_type().is_dir()
                && e.file_name().to_str().unwrap().ends_with("json")
                && !e
                    .path()
                    .parent()
                    .unwrap()
                    .to_str()
                    .unwrap()
                    .ends_with("results") // don't test against files in results output directories
        })
        .collect::<Vec<_>>()
}

// there is one short file that is treated differently to the rest in the e2e tests, i.e. it is given different parameters for a calculation run
const EXCEPTIONAL_DEMO_FILE: &str = "demo_hp_with_setback_separate_ieh_plus_cooling";

#[rstest]
fn test_run_all_files(files: Vec<DirEntry>) {
    let london_weather = cibse_weather_data_to_external_conditions(BufReader::new(Cursor::new(
        include_str!("../examples/weather_data/London_weather_CIBSE_format.csv"),
    )))
    .unwrap();

    files.par_iter().for_each(move |file| {
        let output_writer = InMemoryDirectoryOutputWriter::new(file.file_name().to_str().unwrap());
        let use_additional_options = use_additional_options(file);
        let result = run_project_from_input_file(
            BufReader::new(File::open(file.path()).unwrap()).into(),
            &output_writer,
            Some(london_weather.clone()),
            Some(&vec![OutputFormat::Csv]),
            tariff_file(file),
            use_additional_options,
            use_additional_options,
        );
        assert!(
            result.is_ok(),
            "Error running project for file: {}",
            file.path().display()
        );
        let output_files = output_writer.files();
        println!(
            "Successfully processed file: {}\n{} captured output files compared to expected {}\nEmitted files from run were: {}\n\n",
            file.file_name().display(),
            output_writer.files().len(),
            expected_directory(file)
                .into_iter()
                .filter_map(Result::ok)
                .filter(|f| !f.file_type().is_dir())
                .count(),
            output_files.keys().sorted().join(", ")
        );
    });

    fn use_additional_options(file: &DirEntry) -> bool {
        !file
            .file_name()
            .to_str()
            .unwrap()
            .contains(EXCEPTIONAL_DEMO_FILE)
            && !file.path().to_str().unwrap().contains("long/")
    }

    enum FileKind {
        Short,
        Long,
    }

    impl Display for FileKind {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "{}", self.as_str())
        }
    }

    impl FileKind {
        fn as_str(&self) -> &'static str {
            match self {
                FileKind::Short => "short",
                FileKind::Long => "long",
            }
        }
    }

    impl From<&DirEntry> for FileKind {
        fn from(file: &DirEntry) -> Self {
            if file.path().to_str().unwrap().contains("long/") {
                FileKind::Long
            } else {
                FileKind::Short
            }
        }
    }

    fn tariff_file(file: &DirEntry) -> Option<&str> {
        Path::new(
            if file.file_name().to_str().unwrap().contains("demo_FHS")
                || file
                    .file_name()
                    .to_str()
                    .unwrap()
                    .contains(EXCEPTIONAL_DEMO_FILE)
            {
                "./examples/tariff_data/tariff_data_25-06-2024.csv"
            } else {
                "./examples/tariff_data/tariff_data_demo_files_24timesteps.csv"
            },
        )
        .to_str()
    }

    fn expected_directory(file: &DirEntry) -> WalkDir {
        WalkDir::new(format!(
            "./tests/e2e/expected_results/{}/{}__results",
            FileKind::from(file),
            file.file_name()
                .to_str()
                .unwrap()
                .split('.')
                .next()
                .unwrap()
        ))
    }
}

#[derive(Clone, Debug)]
struct InMemoryDirectoryOutputWriter {
    input_filename: String,
    files: Arc<Mutex<IndexMap<String, FileWriter>>>,
}

impl InMemoryDirectoryOutputWriter {
    fn new(input_filename: &str) -> Self {
        Self {
            input_filename: input_filename.split('.').next().unwrap().to_string(),
            files: Arc::new(Mutex::new(IndexMap::new())),
        }
    }

    fn output_file_index(&self, location_key: &str, file_extension: &str) -> String {
        format!(
            "{}__{}.{}",
            self.input_filename, location_key, file_extension
        )
    }

    pub fn files(&self) -> IndexMap<String, String> {
        self.files
            .lock()
            .iter()
            .map(|(k, v)| (k.clone(), v.0.read().clone()))
            .collect()
    }
}

impl OutputWriter for InMemoryDirectoryOutputWriter {
    fn writer_for_location_key(
        &self,
        location_key: &str,
        file_extension: &str,
    ) -> anyhow::Result<impl Write> {
        Ok(self
            .files
            .lock()
            .entry(self.output_file_index(location_key, file_extension))
            .or_insert_with(FileWriter::new)
            .clone())
    }
}

#[derive(Clone, Debug)]
struct FileWriter(Arc<RwLock<String>>);

impl FileWriter {
    fn new() -> Self {
        Self(Arc::new(RwLock::new(String::with_capacity(2usize.pow(14)))))
    }
}

impl Write for FileWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let utf8 =
            from_utf8(buf).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        self.0.write().push_str(utf8);

        Ok(utf8.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
