#![no_main]

use home_energy_model::corpus::Corpus;
use home_energy_model::input::Input;
use libfuzzer_sys::fuzz_target;
use std::sync::Arc;

fuzz_target!(|input: Input| {
    let _ = Corpus::from_inputs(Arc::new(input), None, None, &Default::default());
});
