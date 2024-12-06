#![no_main]

use hem::corpus::Corpus;
use hem::input::Input;
use libfuzzer_sys::fuzz_target;

fuzz_target!(|input: Input| {
    let _ = Corpus::from_inputs(&input, None, &Default::default());
});
