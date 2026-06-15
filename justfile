# To use these recipes you need to have just installed.
# You can do this using homebrew with `brew install just`.

validate:
    #!/usr/bin/env bash
    cd {{justfile_directory()}}/../epb-home-energy-model-output-validator && cargo run -- --rust-only && cargo test && cd {{justfile_directory()}}

unit:
    cargo test --lib

e2e:
    cargo test --test e2e -- --nocapture

IMAGE := "dev"

INPUT_FILE := "examples/input/core/short/demo_hp_buffer_tank_fancoils.json"

# Build the image (cached)
build:
    docker build -t {{IMAGE}} .


run: build
    docker run --rm -it \
        -v $(pwd):/app \
        -w /app \
        {{IMAGE}} \
        cargo run --release --features="clap indicatif" -- {{INPUT_FILE}} --detailed-output-heating-cooling --heat-balance -e examples/weather_data/London_weather_EnergyPlus_format.epw


docker-unit: build
    docker run --rm -it \
        --memory=16g \
        --memory-swap=16g \
        -v $(pwd):/app \
        -w /app \
        {{IMAGE}} \
        cargo test --lib

docker-e2e: build
    docker run --rm -it \
        --memory=16g \
        --memory-swap=16g \
        -v $(pwd):/app \
        -w /app \
        {{IMAGE}} \
        cargo test --test e2e -- --nocapture

docker-test: build
    docker run --rm -it \
        -v $(pwd):/app \
        -w /app \
        {{IMAGE}} \
        cargo test -- --nocapture

# Remove local Rust build artifacts + Docker leftovers
clean:
    cargo clean
    docker system prune -f

prune:
    docker system prune -f
