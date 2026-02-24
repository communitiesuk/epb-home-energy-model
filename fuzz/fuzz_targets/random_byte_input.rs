#![no_main]

use home_energy_model::output_writer::OutputWriter;
use home_energy_model::run_project_from_input_file;
use libfuzzer_sys::fuzz_target;
use std::io;
use std::io::{BufReader, Cursor, Write};

fuzz_target!(|data: &[u8]| {
    let _run = run_project_from_input_file(
        BufReader::new(Cursor::new(data)).into(),
        &SinkOutput::default(),
        None,
        None,
        None,
        false,
        false,
    );
});

/// An output that goes to nowhere/ a "sink"/ /dev/null.
#[derive(Debug, Default)]
pub struct SinkOutput;

impl OutputWriter for SinkOutput {
    fn writer_for_location_key(
        &self,
        _location_key: &str,
        _file_extension: &str,
    ) -> anyhow::Result<impl Write> {
        Ok(io::sink())
    }

    fn is_noop(&self) -> bool {
        // make the output pretend it's not a no-op so fuzzing exercises code that calls it
        false
    }
}
