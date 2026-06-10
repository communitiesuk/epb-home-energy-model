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

# Build the image (cached)
build:
    docker build -t {{IMAGE}} .

docker-unit: build
    docker run --rm -it \
        -v $(pwd):/app \
        -w /app \
        {{IMAGE}} \
        cargo test --lib

docker-e2e: build
    docker run --rm -it \
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
    rm -rf target
    docker system prune -f
