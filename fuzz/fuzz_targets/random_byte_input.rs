#![no_main]

use hem::output::SinkOutput;
use hem::run_project;
use libfuzzer_sys::fuzz_target;
use std::io::{BufReader, Cursor};

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
