#![no_main]

use hem::output::Output;
use hem::run_project;
use hem::ProjectFlags;
use libfuzzer_sys::fuzz_target;
use std::io;
use std::io::{BufReader, Cursor, Write};

fuzz_target!(|data: &[u8]| {
    let _run = run_project(
        BufReader::new(Cursor::new(data)),
        SinkOutput::default(),
        None,
        &ProjectFlags::empty(),
    );
});

/// An output that goes to nowhere/ a "sink"/ /dev/null.
#[derive(Debug, Default)]
pub struct SinkOutput;

impl Output for SinkOutput {
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
