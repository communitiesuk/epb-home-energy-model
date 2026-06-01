FROM mcr.microsoft.com/devcontainers/rust:1-trixie

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    python3-numpy \
    python3-scipy \
    pkg-config \
    build-essential \
    gfortran \
    libopenblas-dev \
    liblapack-dev \
    lldb \
    just \
    && rm -rf /var/lib/apt/lists/*

ENV PYO3_PYTHON=/usr/bin/python3
ENV PYTHON_SYS_EXECUTABLE=/usr/bin/python3

USER vscode
