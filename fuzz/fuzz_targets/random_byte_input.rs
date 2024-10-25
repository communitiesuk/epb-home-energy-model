#![no_main]

use hem::output::Output;
use hem::run_project;
use libfuzzer_sys::fuzz_target;
use std::io;
use std::io::{BufReader, Cursor, Write};

fuzz_target!(|data: &[u8]| {
    let _run = run_project(
        BufReader::new(Cursor::new(data)),
        SinkOutput::default(),
        None,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
    );
});

/// An output that goes to nowhere/ a "sink"/ /dev/null.
#[derive(Debug, Default)]
pub struct SinkOutput;

impl Output for SinkOutput {
    fn writer_for_location_key(&self, _location_key: &str) -> anyhow::Result<impl Write> {
        Ok(io::sink())
    }

    fn is_noop(&self) -> bool {
        // make the output pretend it's a no-op so fuzzing exercises code that calls it
        false
    }
}
