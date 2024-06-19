# To use these recipes you need to have just installed.
# You can do this using homebrew with `brew install just`.

validate:
    #!/usr/bin/env bash
    cd {{justfile_directory()}}/../epb-home-energy-model-output-validator && cargo run -- --rust-only && cargo test && cd {{justfile_directory()}}
