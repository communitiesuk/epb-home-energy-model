# Home Energy Model (HEM) engine written in Rust

Requires the `rustup` toolchain to use. ([Instructions](https://rustup.rs) for installation. For MacOS, don't use
Homebrew to install Rust.)

To run tests:

```bash
cargo test
```

Fully running the engine requires a weather file in EPW format and an input JSON file (format not yet documented/
stable):

```
cargo run --release --features="clap indicatif" -- path/to/input.json -e path/to/weather_file.epw
```

The `clap` feature above is mandatory. The `indicatif` feature switches on the output of a progress bar.
