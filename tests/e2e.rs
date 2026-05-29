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
        let actual_file_count = output_files.len();
        let expected_files: IndexMap<String, DirEntry> = expected_directory(file)
            .into_iter()
            .filter_map(Result::ok)
            .filter(|f| !f.file_type().is_dir())
            .map(|f| (f.file_name().to_str().unwrap().to_owned(), f))
            .collect();
        let expected_file_count = expected_files.len();
        for (actual_file_name, actual_file) in output_files.iter() {
            let expected_file = expected_files.get(actual_file_name);
            let expected_file = if let Some(expected_file) = expected_file {
                expected_file
            } else {
                println!("🐙 Unexpected file emitted: {}", actual_file_name);
                continue;
            };
            let mut rust_file_read = BufReader::new(Cursor::new(actual_file));
            let rust_headers = csv_access::csv_reader(&mut rust_file_read).headers().unwrap().clone();
            let mut python_file_read = BufReader::new(File::open(expected_file.path()).unwrap());
            let python_headers = csv_access::csv_reader(&mut python_file_read).headers().unwrap().clone();
            let differences = compare::compare(python_headers, rust_headers);
            if let Err(differences) = differences {
                println!("❌ Headers differ for file: {}", actual_file_name);
                println!("Differences: {}", differences.iter().join("\n"));
            }
        }
        println!(
            "Successfully processed file: {}\n{} {} captured output files compared to expected {}{}\n\n",
            file.file_name().display(),
            if actual_file_count < expected_file_count {
                "❌"
            } else {
                "✅"
            },
            actual_file_count,
            expected_file_count,
            if actual_file_count < expected_file_count {
                format!("\nEmitted files from run were: {}", output_files.keys().sorted().join(", "))
            } else {
                "".into()
            },
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

mod csv_access {
    use std::io::Read;

    pub fn csv_reader<T: Read>(read: &mut T) -> csv::Reader<&mut T> {
        csv::ReaderBuilder::new().flexible(true).from_reader(read)
    }
}

mod compare {
    use csv::StringRecord;
    use std::fmt;
    use std::mem::discriminant;

    pub fn compare(left: StringRecord, right: StringRecord) -> ComparisonResult {
        OutputRecord::from(left).equiv(&OutputRecord::from(right))
    }

    const FLOAT_THRESHOLD: f64 = 1e-6; // 0.000001

    #[derive(Debug, Clone)]
    pub enum Difference {
        StringDifference {
            left: String,
            right: String,
            field_index: usize,
        },
        NumberDifference {
            left: f64,
            right: f64,
            numerical_difference: f64,
            field_index: usize,
        },
        RecordDifference {
            message: String,
        },
    }

    impl fmt::Display for Difference {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            // assumption here that left is Python and right is Rust
            match self {
                Difference::StringDifference {
                    left,
                    right,
                    field_index,
                } => {
                    write!(f, "column {field_index}, 🐍: \"{left}\", 🦀: \"{right}\"")
                }
                Difference::NumberDifference {
                    left,
                    right,
                    numerical_difference,
                    field_index,
                } => {
                    write!(
                        f,
                        "column {field_index}, 🐍: {left}, 🦀: {right}, Diff: {numerical_difference}"
                    )
                }
                Difference::RecordDifference { message } => {
                    write!(f, "{}", message)
                }
            }
        }
    }

    pub struct OutputRecord {
        record: StringRecord,
    }

    impl From<StringRecord> for OutputRecord {
        fn from(value: StringRecord) -> Self {
            Self { record: value }
        }
    }

    impl PartialEq for OutputRecord {
        fn eq(&self, other: &Self) -> bool {
            self.equiv(other).is_ok()
        }
    }

    type ComparisonResult = Result<(), Vec<Difference>>;

    impl OutputRecord {
        pub fn equiv(&self, other: &OutputRecord) -> ComparisonResult {
            if self.len() != other.len() {
                let message = format!(
                    "Record has unequal number of fields, {} to {}.",
                    self.len(),
                    other.len()
                );
                return Err(vec![Difference::RecordDifference { message }]);
            }
            let differences = self
                .record
                .iter()
                .zip(&other.record)
                .map(|(left_str, right_str)| {
                    (
                        OutputCellValue::from(left_str),
                        OutputCellValue::from(right_str),
                    )
                })
                .enumerate()
                // Result<(), Difference>
                .filter_map(|(index, (left, right))| {
                    let comparison = left.equiv(&right, index);
                    match comparison {
                        Ok(_) => None,
                        Err(difference) => Some(difference),
                    }
                })
                .collect::<Vec<_>>();

            if differences.is_empty() {
                Ok(())
            } else {
                Err(differences)
            }
        }

        fn len(&self) -> usize {
            self.record.len()
        }
    }

    #[derive(Debug)]
    enum OutputCellValue {
        Number(f64),
        String(String),
    }

    impl OutputCellValue {
        fn equiv(&self, other: &OutputCellValue, field_index: usize) -> Result<(), Difference> {
            if discriminant(self) != discriminant(other) {
                return Err(Difference::StringDifference {
                    left: format!("{:?}", self),
                    right: format!("{:?}", other),
                    field_index,
                });
            }
            match self {
                OutputCellValue::Number(float) => {
                    let other_float = match other {
                        OutputCellValue::Number(other_float) => *other_float,
                        _ => unreachable!(),
                    };
                    let numerical_difference = (*float - other_float).abs();
                    if numerical_difference < FLOAT_THRESHOLD {
                        Ok(())
                    } else {
                        Err(Difference::NumberDifference {
                            left: *float,
                            right: other_float,
                            numerical_difference,
                            field_index,
                        })
                    }
                }
                OutputCellValue::String(string) => {
                    let other_string = match other {
                        OutputCellValue::String(string) => string,
                        _ => unreachable!(),
                    };

                    if string == other_string {
                        Ok(())
                    } else {
                        Err(Difference::StringDifference {
                            left: string.clone(),
                            right: other_string.clone(),
                            field_index,
                        })
                    }
                }
            }
        }
    }

    impl From<&str> for OutputCellValue {
        fn from(value: &str) -> Self {
            if let Ok(float) = value.parse::<f64>() {
                Self::Number(float)
            } else {
                Self::String(value.to_string())
            }
        }
    }
}
