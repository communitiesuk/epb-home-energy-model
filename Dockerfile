FROM rust:slim

# Install dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-dev python3-pip python3-venv \
    python3-numpy python3-scipy \
    pkg-config build-essential gfortran \
    libopenblas-dev liblapack-dev \
    lldb just \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV PYO3_PYTHON=/usr/bin/python3
ENV PYTHON_SYS_EXECUTABLE=/usr/bin/python3

WORKDIR /app

# Copy only Cargo.toml first to cache dependencies
COPY Cargo.toml Cargo.lock ./

# Create an empty src to allow cargo to resolve deps
RUN mkdir src && echo "fn main() {}" > src/main.rs

# Pre-fetch and compile dependencies (cached)
RUN cargo fetch && cargo build --release || true

# Now copy the actual source
COPY . .

CMD ["/bin/bash"]
