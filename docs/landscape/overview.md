# Overview

## Deployment

Deployment is performed by a front-end application maintained in a [separate repository](https://github.com/communitiesuk/epb-ecaas-frontend).

## As-Is Architectural Map

This map outlines the primary domains within the existing codebase and their responsibilities.

```mermaid
graph TD
    subgraph "Application Entrypoints"
        A1["src/main.rs"]
        A2["hem-lambda/src/main.rs"]
    end

    subgraph "Core Simulation Engine"
        B1["src/lib.rs (run_project)"]
        B2["src/corpus.rs (God Object)"]
        B3["src/core/*"]
        B4["src/hem_core/*"]
    end

    subgraph "Input Processing & Validation"
        C1["schemas/core-input.schema.json"]
        C2["src/input.rs"]
        C3["examples/input/core/*"]
    end

    subgraph "External Data Integration"
        D1["src/read_weather_file.rs"]
        D2["src/core/energy_supply/tariff_data.rs"]
        D3["examples/weather_data/*"]
        D4["examples/tariff_data/*"]
    end

    subgraph "Output & Data Handling"
        E1["src/output.rs"]
        E2["src/output_writer.rs"]
        E3["src/lib.rs (write_core_output_files)"]
    end

    subgraph "Testing & Quality Assurance"
        F1["fuzz/*"]
        F2[".github/workflows/test.yml"]
        F3["src/compare_floats.rs"]
    end

    A1 --> B1
    A2 --> B1
    B1 -- uses --> B2
    B2 -- uses --> B3
    B2 -- uses --> B4
    B1 -- uses --> C2
    C1 -- "Defines structure for" --> C2
    C2 -- "Reads from" --> C3
    B2 -- uses --> D1
    B2 -- uses --> D2
    D1 -- "Reads from" --> D3
    D2 -- "Reads from" --> D4
    B1 -- uses --> E1
    B1 -- uses --> E2
    E2 -- "Writes output using format from" --> E3
   ```

### 1. Core Simulation Engine
- **Description:** The heart of the application, responsible for executing the energy model simulation based on the provided inputs. The current implementation is highly monolithic.
- **Key Files:**
    - `src/corpus.rs`: A "God Object" that holds the entire simulation state and orchestrates all steps.
    - `src/lib.rs`: The main library entry point that kicks off the simulation (`run_project`).
    - `src/core/*`: Modules representing the physical components and concepts of the energy model (e.g., `heating_systems`, `space_heat_demand`).
    - `src/hem_core/*`: Modules for fundamental simulation concepts like time (`simulation_time.rs`) and external conditions (`external_conditions.rs`).

### 2. Input Processing & Validation
- **Description:** Responsible for deserializing and validating the complex simulation inputs from JSON files.
- **Key Files:**
    - `src/input.rs`: Defines the large, deeply-nested Rust structs that mirror the input JSON.
    - `schemas/core-input.schema.json`: The formal JSON schema that defines the structure and constraints of the input files.
    - `examples/input/core/*`: A collection of example JSON input files that demonstrate various simulation scenarios.

### 3. Output & Data Handling
- **Description:** Responsible for collecting, formatting, and writing the simulation results to CSV and other formats.
- **Key Files:**
    - `src/output.rs`: Defines the data structures for the simulation output.
    - `src/output_writer.rs`: Contains the logic for writing data to files.
    - `src/lib.rs` (specifically `write_core_output_files`): Contains complex, imperative logic for transforming internal state into the desired output format.

### 4. External Data Integration
- **Description:** Handles the loading and integration of external data required for the simulation, such as weather and electricity tariffs.
- **Key Files:**
    - `src/read_weather_file.rs`: Logic for parsing weather data files.
    - `src/core/energy_supply/tariff_data.rs`: Logic for processing electricity tariff information.
    - `examples/weather_data/*`: Example weather data files (both `.csv` and `.epw` formats).
    - `examples/tariff_data/*`: Example tariff data files.

### 5. Application Entrypoints
- **Description:** The executable "runners" that consume the core simulation library.
- **Key Files:**
    - `src/main.rs`: The primary command-line interface for running simulations.
    - `hem-lambda/src/main.rs`: An adapter to run the simulation within an AWS Lambda environment.

### 6. Testing & Quality Assurance
- **Description:** The collection of tools and code dedicated to ensuring the correctness and stability of the model.
- **Key Files:**
    - `fuzz/*`: Fuzz testing harnesses to discover panics and bugs with random inputs.
    - `.github/workflows/test.yml`: The CI pipeline definition for running tests automatically.
    - `src/compare_floats.rs`: A utility module to handle the complexities of comparing floating-point numbers in tests.