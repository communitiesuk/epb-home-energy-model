use crate::compare::DifferenceKind;
use home_energy_model::output_writer::OutputWriter;
use home_energy_model::read_weather_file::cibse_weather_data_to_external_conditions;
use home_energy_model::{run_project_from_input_file, OutputFormat};
use indexmap::IndexMap;
use itertools::Itertools;
use parking_lot::{Mutex, RwLock};
use rayon::prelude::*;
use rstest::*;
use std::collections::HashSet;
use std::fmt::Display;
use std::fs::File;
use std::io::{BufReader, Cursor, Write};
use std::path::Path;
use std::str::from_utf8;
use std::sync::Arc;
use walkdir::{DirEntry, WalkDir};

const PASSING_FILES: &[&str] = &[
    "demo_24hrs_January.json",
    "demo_24hrs_January_ieh.json",
    "demo_24hrs_January_MVHR_external.json",
    "demo.json",
    "demo_elec_battery_8.json",
    "demo_curtains_and_blinds.json",
    "demo_hwohp.json",
    "demo_suspended_floor.json",
    "demo_new_infiltration_vent_model_MVHR.json",
    "SAP11_deck_smart_hot_water_tank.json",
    "SAP11_deck_flat_nat_vent_with_window_opening_for_cooling.json",
    "demo_new_infiltration_vent_model_MVHR_inside.json",
    "demo_new_infiltration_vent_model.json",
    "demo_new_infiltration_vent_model_intermit_mev.json",
    "demo_hp_warm_air.json",
    "demo_FHS_with_setback.json",
    "demo_FHS_storage_tank_held_at_setpnt.json",
    "demo_FHS_smart_hot_water_tank.json",
    "demo_24hrs_January_instant_elec_heater_comp.json",
    "demo_24hrs_August_WWHRS_A.json",
    "demo_24hrs_August_WWHRS_B.json",
    "demo_24hrs_August_WWHRS_C.json",
    "demo_24hrs_August_nearby_shading_obstacle.json",
    "demo_elec_battery_inverter_inside.json",
    "demo_24hrs_August_pvdiverter.json",
    "demo_168hrs_loadshifting_appliances.json",
    "demo_168hrs_loadshifting_appliances_flat_weights+limit_0.json",
    "demo_FHS_default_hw_sched_allday.json",
    "demo_FHS_default_hw_sched_heatinghours.json",
    "demo_int_air_setpoint.json",
    "demo_24hrs_August.json",
    "demo_24hrs_August_remote_shading_objects.json",
    "demo_24hrs_August_cost_minimising_time_ctrl.json",
    "demo_24hrs_January_ieh_6_1kW.json",
    "demo_two_hot_water_sources_associated_outlets.json",
    "demo_elec_battery.json",
    "demo_24hrs_January_MVHR.json",
    "demo_FHS_U_values.json",
    "demo_FHS.json",
    "demo_FHS_battery_no_scope_for_charging.json",
    "demo_FHS_cooling_demand_beyond_free_ventilation.json",
    "demo_FHS_battery_capacity_reduction_due_to_temp.json",
    "demo_unconditioned.json",
    "demo_24hrs_January_point_of_use.json",
    "demo_24hrs_August_SolarThermal.json",
    "demo_24hrs_August_SolarThermal_to_preheat.json",
];

const PASSING_FILES_IN_USE_PYTHON_ONLY: &[&str] = &[
    "demo_hp_smart_hot_water_tank.json",
    "demo_hp_default_to_max.json",
    "demo_combiBoiler.json",
    "demo_hp_buffer_tank.json",
    "demo_hp_ufh.json",
    "demo_hp_with_setback_separate_ieh_same_setpoint.json",
    "demo_24hrs_August_WWHRS.json",
    "demo_hp_with_setback.json",
    "demo_hp_with_setback_separate_ieh_diff_setpoint.json",
    "demo_emitter_pipework.json",
    "demo_hp_surfacewater.json",
    "demo_hp_bypass.json",
    "demo_eahp_mixed.json",
    "demo_hp_with_advancedstart.json",
    "demo_heat_network_5G.json",
    "demo_eahp_single_zone.json",
    "demo_eahp.json",
    "demo_24hrs_January_esh_automatic.json",
    "demo_24hrs_January_esh_celect.json",
    "demo_heat_battery_space_heat.json",
    "demo_24hrs_January_esh_hhrsh.json",
    "demo_24hrs_January_esh_manual.json",
    "demo_hp_buffer_tank_fancoils.json",
    "demo_24hrs_August_pvdiverter_and_hp.json",
    "demo_heat_battery_water_only.json",
    "demo_heat_battery_all.json",
    "demo_heat_battery_charge_level.json",
    "demo_FHS_heat_battery.json",
    "demo_FHS_heating_system_priority.json",
    "demo_FHS.json",
    "demo_hp_primary_pipework.json",
    "demo_FHS_emitters_outside_temp_over_maximum.json",
    "demo_24hrs_August_pvdiverter_and_hp_smart_hot_water.json",
    "demo_combiBoilerLPG.json",
];

#[fixture]
fn files() -> Vec<DirEntry> {
    WalkDir::new("./examples/input/core")
        .into_iter()
        .filter_map(Result::ok)
        .filter(|e| {
            !e.file_type().is_dir()
                && e.file_name().to_str().unwrap().ends_with("json")
                && !PASSING_FILES.contains(&e.file_name().to_str().unwrap())
                && !PASSING_FILES_IN_USE_PYTHON_ONLY.contains(&e.file_name().to_str().unwrap())
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

const NON_TABULAR_FILE_SUFFIXES: &[&str] = &[
    "results_summary.csv",
    "results_static.csv",
    "postproc.csv",
    "postproc_summary.csv",
];

const NON_TABULAR_FILE_FRAGMENTS: &[&str] = &["_summary__"];

#[rstest]
fn test_run_all_files(files: Vec<DirEntry>) {
    let london_weather = cibse_weather_data_to_external_conditions(BufReader::new(Cursor::new(
        include_str!("../examples/weather_data/London_weather_CIBSE_format.csv"),
    )))
    .unwrap();

    let tracked_files: Arc<RwLock<HashSet<String>>> = Default::default();

    let difference_count: usize = files.par_iter().map(move |file| {
        println!("\n🎬 starting to run HEM calculation on file {}\n", file.file_name().display());
        let mut difference_count = 0usize;
        let output_writer = InMemoryDirectoryOutputWriter::new(file.file_name().to_str().unwrap());
        let use_additional_options = use_additional_options(file);
        let file_name_string = file.file_name().to_str().unwrap().to_string();
        {
            let mut tracked_files = tracked_files.write();
            tracked_files.insert(file_name_string.clone());
            println!("🏃‍♂️ starting new calc - current count of calculations under way: {}", tracked_files.len());
        }
        let result = run_project_from_input_file(
            BufReader::new(File::open(file.path()).unwrap()).into(),
            &output_writer,
            Some(london_weather.clone()),
            Some(&vec![OutputFormat::Csv]),
            tariff_file(file),
            use_additional_options,
            use_additional_options,
        );
        {
            let mut tracked_files = tracked_files.write();
            println!("🏁 finished calculation for {}", &file_name_string);
            tracked_files.remove(&file_name_string);
            println!("📋 {} files left being calculated: {}", tracked_files.len(), tracked_files.iter().join(", "));
        }
        if let Err(e) = result {
            println!("💥 Error running project for file (100,000 difference penalty!) {}: {}", file.path().display(), e);
            difference_count += 100000;
            return difference_count;
        }
        let output_files = output_writer.files();
        let actual_file_count = output_files.len();
        let expected_files: IndexMap<String, DirEntry> = expected_directory(file)
            .into_iter()
            .filter_map(Result::ok)
            .filter(|f| !f.file_type().is_dir())
            .map(|f| (f.file_name().to_str().unwrap().to_owned(), f))
            .collect();
        let expected_file_count = expected_files.len();
        for (file_name, actual_file) in output_files.iter() {
            let expected_file = expected_files.get(file_name);
            let expected_file = if let Some(expected_file) = expected_file {
                expected_file
            } else {
                println!("🐙 Unexpected file emitted: {}", file_name);
                continue;
            };
            let mut rust_file_read = BufReader::new(Cursor::new(actual_file));
            let mut rust_reader = csv_access::csv_reader(&mut rust_file_read);
            let mut python_file_read = BufReader::new(File::open(expected_file.path()).unwrap());
            let mut python_reader = csv_access::csv_reader(&mut python_file_read);
            let header_differences = compare::compare_headers(&mut python_reader, &mut rust_reader);
            if let Err(differences) = header_differences {
                println!("❌ Headers differ for file: {}", file_name);
                println!("Differences: {}", differences.iter().join("\n"));
                difference_count += differences.len();
            }
            let difference_kind = DifferenceKind::CountOnly;
            if is_tabular(file_name) {
                let file_differences = compare::compare_tabular_records_within_threshold(&mut python_reader, &mut rust_reader, difference_kind);
                if let Err(comparison_error) = file_differences {
                    let file_difference_count = comparison_error.differences.len();
                    println!("❌ Tabular records differ for file: {} - difference count is {}", file_name, file_difference_count);
                    difference_count += file_difference_count;
                } else {
                    println!("✅ Tabular records match for file: {}", file_name);
                }
            } else {
                let file_differences = compare::compare_non_tabular_files(&mut python_reader, &mut rust_reader, difference_kind);
                if let Err(comparison_error) = file_differences {
                    let file_difference_count = comparison_error.differences.len();
                    println!("❌ Non-tabular records differ for file: {} - difference count is {}", file_name, file_difference_count);
                    difference_count += file_difference_count;
                } else {
                    println!("✅ Non-tabular records match for file: {}", file_name);
                }
            }
        }
        println!(
            "Successfully processed file: {}\n{} {} captured output files compared to expected {}{}\n\n",
            file.file_name().display(),
            if actual_file_count < expected_file_count {
                "🚮"
            } else {
                "🎉"
            },
            actual_file_count,
            expected_file_count,
            if actual_file_count < expected_file_count {
                format!("\nEmitted files from run were: {}", output_files.keys().sorted().join(", "))
            } else {
                "".into()
            },
        );

        difference_count
    }).sum();

    assert_eq!(
        difference_count, 0,
        "Total difference count: {}",
        difference_count
    );

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

    fn is_tabular(file_name: &str) -> bool {
        !NON_TABULAR_FILE_SUFFIXES
            .iter()
            .any(|suffix| file_name.ends_with(suffix))
            && !NON_TABULAR_FILE_FRAGMENTS
                .iter()
                .any(|fragment| file_name.contains(fragment))
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
    use csv::{Reader, StringRecord};
    use std::fmt;
    use std::io::Read;
    use std::mem::discriminant;
    use thiserror::Error;

    pub fn compare(left: StringRecord, right: StringRecord) -> ComparisonResult {
        OutputRecord::from(left).equiv(&OutputRecord::from(right))
    }

    const FLOAT_THRESHOLD: f64 = 1e-6; // 0.000001

    #[derive(Debug, Clone)]
    pub enum Difference {
        String {
            left: String,
            right: String,
            field_index: usize,
        },
        Number {
            left: f64,
            right: f64,
            numerical_difference: f64,
            field_index: usize,
        },
        Record {
            message: String,
        },
    }

    impl Difference {
        fn field_index(&self) -> Option<usize> {
            match self {
                Difference::String { field_index, .. } => Some(*field_index),
                Difference::Number { field_index, .. } => Some(*field_index),
                Difference::Record { .. } => None,
            }
        }
    }

    impl fmt::Display for Difference {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            // assumption here that left is Python and right is Rust
            match self {
                Difference::String {
                    left,
                    right,
                    field_index,
                } => {
                    write!(f, "column {field_index}, 🐍: \"{left}\", 🦀: \"{right}\"")
                }
                Difference::Number {
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
                Difference::Record { message } => {
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

    #[derive(Debug)]
    #[allow(dead_code)]
    pub struct DifferenceInFile {
        pub line_number: usize,
        pub column_field: Option<String>,
        pub difference: Difference,
    }

    impl DifferenceInFile {
        fn new(difference: Difference, line_number: usize, column_field: Option<String>) -> Self {
            Self {
                difference,
                line_number,
                column_field,
            }
        }
    }

    type FileComparisonResult = Result<(), FileComparisonError>;

    impl OutputRecord {
        pub fn equiv(&self, other: &OutputRecord) -> ComparisonResult {
            if self.len() != other.len() {
                let message = format!(
                    "Record has unequal number of fields, {} to {}.",
                    self.len(),
                    other.len()
                );
                return Err(vec![Difference::Record { message }]);
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
                .filter_map(|(index, (left, right))| left.equiv(&right, index).err())
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
                return Err(Difference::String {
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
                        Err(Difference::Number {
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
                        // deal with boolean "True"/"False" values rendered as strings, as these should be compared case-insensitively
                        let lowercase_bools = ["true", "false"];
                        if lowercase_bools.contains(&string.to_lowercase().as_str())
                            && lowercase_bools.contains(&other_string.to_lowercase().as_str())
                            && string.to_lowercase() == other_string.to_lowercase()
                        {
                            Ok(())
                        } else {
                            Err(Difference::String {
                                left: string.clone(),
                                right: other_string.clone(),
                                field_index,
                            })
                        }
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

    pub fn compare_headers<T, U>(
        py_reader: &mut Reader<T>,
        rust_reader: &mut Reader<U>,
    ) -> ComparisonResult
    where
        T: Read,
        U: Read,
    {
        let rust_headers = rust_reader
            .headers()
            .expect("Failed to read Rust CSV headers");
        let py_headers = py_reader
            .headers()
            .expect("Failed to read Python CSV headers");

        compare(py_headers.clone(), rust_headers.clone())
    }

    pub fn compare_tabular_records_within_threshold<T, U>(
        py_reader: &mut Reader<T>,
        rust_reader: &mut Reader<U>,
        difference_kind: DifferenceKind,
    ) -> FileComparisonResult
    where
        T: Read,
        U: Read,
    {
        let headers;
        {
            // read headers in own scope so we can re-use py_reader later
            let python_headers = py_reader.headers().expect("Failed to read Python headers");
            headers = python_headers.clone();
        }

        let mut rust_records = rust_reader.records().enumerate();
        let python_records = py_reader.records().enumerate();

        let mut file_differences: Vec<DifferenceInFile> = vec![];
        let mut difference_count = usize::default();

        let mut warnings: Vec<String> = vec![];

        for (record_index, python_record) in python_records {
            let rust_record = if let Some((_, rust_record)) = rust_records.next() {
                rust_record
            } else {
                warnings.push(format!(
                    "Rust records are missing from index {}",
                    record_index
                ));
                break;
            };

            if let Err(ref rust_record) = &rust_record {
                warnings.push(format!(
                    "Failed to parse Rust record at index {}: {}",
                    record_index, rust_record
                ));
            }
            if let Err(ref python_record) = &python_record {
                warnings.push(format!(
                    "Failed to parse Python record at index {}: {}",
                    record_index, python_record
                ));
            }

            let (mut rust_record, python_record) =
                if let (Ok(rust_record), Ok(python_record)) = (rust_record, python_record) {
                    (rust_record, python_record)
                } else {
                    continue;
                };

            // hack to handle python blank lines (no data at all) being ignored
            // but rust blank lines (double quotes) being included
            // skip over rust records with just "" on them
            let mut blank_lines = usize::default();
            while rust_record.as_slice() == "" && python_record.as_slice() != "" {
                blank_lines += 1;
                rust_record = rust_records.next().map(|(_, res)| res.unwrap()).unwrap();
            }

            if let Err(differences) = compare(python_record, rust_record) {
                match difference_kind {
                    DifferenceKind::Full => {
                        file_differences.extend(differences.into_iter().map(
                            |difference: Difference| {
                                let column_field = difference
                                    .field_index()
                                    .and_then(|index| headers.get(index))
                                    .map(ToOwned::to_owned);

                                DifferenceInFile::new(
                                    difference,
                                    record_index + blank_lines + 2,
                                    column_field,
                                )
                            },
                        ));
                    }
                    DifferenceKind::CountOnly => {
                        difference_count += differences.len();
                    }
                }
            }
        }

        if match difference_kind {
            DifferenceKind::Full => file_differences.is_empty(),
            DifferenceKind::CountOnly => difference_count == 0,
        } {
            Ok(())
        } else {
            Err(FileComparisonError {
                differences: if difference_kind == DifferenceKind::CountOnly {
                    difference_count.into()
                } else {
                    file_differences.into()
                },
                warnings,
            })
        }
    }

    #[derive(Debug, Error)]
    #[error("File comparison failed with {} differences and {} warnings", self.differences.len(), self.warnings.len())]
    pub struct FileComparisonError {
        pub differences: FileDifferences,
        pub warnings: Vec<String>,
    }

    #[derive(Debug)]
    pub enum FileDifferences {
        List(Vec<DifferenceInFile>),
        Count(usize),
    }

    impl FileDifferences {
        pub fn len(&self) -> usize {
            match self {
                FileDifferences::List(differences) => differences.len(),
                FileDifferences::Count(count) => *count,
            }
        }
    }

    impl From<Vec<DifferenceInFile>> for FileDifferences {
        fn from(differences: Vec<DifferenceInFile>) -> Self {
            FileDifferences::List(differences)
        }
    }

    impl From<usize> for FileDifferences {
        fn from(count: usize) -> Self {
            FileDifferences::Count(count)
        }
    }

    #[allow(dead_code)]
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub enum DifferenceKind {
        Full,
        CountOnly,
    }

    pub fn compare_non_tabular_files<T, U>(
        py_reader: &mut Reader<T>,
        rust_reader: &mut Reader<U>,
        difference_kind: DifferenceKind,
    ) -> FileComparisonResult
    where
        T: Read,
        U: Read,
    {
        let mut rust_records = rust_reader.records();
        let python_records = py_reader.records();

        let mut blank_lines = usize::default();

        let mut file_differences: Vec<DifferenceInFile> = vec![];
        let mut difference_count = usize::default();

        for (record_index, python_record) in python_records.enumerate() {
            let python_record = python_record.unwrap();

            let mut rust_record = if let Some(next_record) = rust_records.next() {
                next_record.unwrap()
            } else {
                match difference_kind {
                    DifferenceKind::Full => {
                        file_differences.push(DifferenceInFile::new(
                            Difference::Record {
                                message: format!("Rust file only had {} lines", record_index + 1),
                            },
                            0,
                            None,
                        ));
                    }
                    DifferenceKind::CountOnly => {
                        difference_count += 1;
                    }
                }

                break;
            };

            // hack to handle python blank lines (no data at all) being ignored
            // but rust blank lines (double quotes) being included
            //
            // skip over rust records with just "" on them
            while rust_record.as_slice().trim().is_empty() {
                blank_lines += 1;
                rust_record = rust_records.next().unwrap().unwrap();
            }

            if let Err(differences) = compare(python_record, rust_record) {
                for difference in differences {
                    match difference_kind {
                        DifferenceKind::Full => {
                            file_differences.push(DifferenceInFile::new(
                                difference,
                                record_index + blank_lines + 1,
                                None,
                            ));
                        }
                        DifferenceKind::CountOnly => {
                            difference_count += 1;
                        }
                    }
                }
            }
        }

        if match difference_kind {
            DifferenceKind::Full => file_differences.is_empty(),
            DifferenceKind::CountOnly => difference_count == 0,
        } {
            Ok(())
        } else {
            Err(FileComparisonError {
                differences: if difference_kind == DifferenceKind::CountOnly {
                    difference_count.into()
                } else {
                    file_differences.into()
                },
                warnings: vec![],
            })
        }
    }
}
