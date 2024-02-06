mod compare_floats;
pub mod core;
mod corpus;
mod external_conditions;
mod input;
pub mod read_weather_file;
mod simulation_time;

#[macro_use]
extern crate is_close;
extern crate lazy_static;

use crate::input::{parse_input_file, ExternalConditionsInput, Input, Ventilation};
use crate::read_weather_file::ExternalConditions as ExternalConditionsFromFile;
use std::collections::HashMap;

use crate::core::units::{SECONDS_PER_HOUR, WATTS_PER_KILOWATT};
use crate::external_conditions::{DaylightSavingsConfig, ExternalConditions};
use crate::simulation_time::{SimulationTime, SimulationTimeIterator};
use std::error::Error;
use std::ffi::OsStr;
use std::ops::Deref;
use std::path::Path;
use std::sync::Arc;

pub fn run_project(
    input_file: &str,
    external_conditions_data: Option<ExternalConditionsFromFile>,
    _preprocess_only: bool,
    _fhs_assumptions: bool,
    _fhs_fee_assumptions: bool,
    _fhs_notA_assumptions: bool,
    _fhs_notB_assumptions: bool,
    _heat_balance: bool,
) -> Result<(), Box<dyn Error>> {
    let input_file_ext = Path::new(input_file).extension().and_then(OsStr::to_str);
    let input_file_stem = match input_file_ext {
        Some(ext) => &input_file[..(input_file.len() - ext.len() - 1)],
        None => input_file,
    };
    let _output_file = format_args!("{input_file_stem}_results.csv");
    let _output_file_static = format_args!("{input_file_stem}_results_static.csv");
    let _output_file_summary = format_args!("{input_file_stem}_results_summary.csv");

    println!("about to try and open {}", input_file);

    let project_data = parse_input_file(Path::new(input_file));

    println!("{:?}", project_data);

    let input = project_data.unwrap();

    let external_conditions = external_conditions_from_input(
        input.external_conditions.clone(),
        external_conditions_data,
        input.simulation_time.clone(),
    );

    let simulation_time = *input.simulation_time.clone().deref();

    Calculation::new(input, external_conditions, simulation_time).run();

    Ok(())
}

fn external_conditions_from_input<'a>(
    input: Arc<ExternalConditionsInput>,
    external_conditions_data: Option<ExternalConditionsFromFile>,
    simulation_time: Arc<SimulationTime>,
) -> ExternalConditions {
    match external_conditions_data {
        Some(ec) => ExternalConditions::new(
            &simulation_time.iter(),
            ec.air_temperatures,
            ec.wind_speeds,
            ec.diffuse_horizontal_radiation,
            ec.direct_beam_radiation,
            ec.solar_reflectivity_of_ground,
            ec.latitude,
            ec.longitude,
            0,
            0,
            Some(365),
            1.0,
            None,
            DaylightSavingsConfig::NotApplicable,
            false,
            ec.direct_beam_conversion_needed,
            input.shading_segments.clone(),
        ),
        None => ExternalConditions::new(
            &simulation_time.iter(),
            input.air_temperatures.clone().unwrap_or(vec![]),
            input.wind_speeds.clone().unwrap_or(vec![]),
            input.diffuse_horizontal_radiation.clone().unwrap_or(vec![]),
            input.direct_beam_radiation.clone().unwrap_or(vec![]),
            input.solar_reflectivity_of_ground.clone().unwrap_or(vec![]),
            input.latitude.unwrap_or(55.0),
            input.longitude.unwrap_or(0.0),
            0,
            0,
            Some(365),
            1.0,
            None,
            DaylightSavingsConfig::NotApplicable,
            false,
            input.direct_beam_conversion_needed.unwrap_or(false),
            input.shading_segments.clone(), //imperfect but this should be quite small...
        ),
    }
}

struct Calculation {
    corpus: Input,
    external_conditions: ExternalConditions,
    simulation_time: SimulationTime,
}

impl Calculation {
    pub fn new(
        corpus: Input,
        external_conditions: ExternalConditions,
        simulation_time: SimulationTime,
    ) -> Self {
        Self {
            corpus,
            external_conditions,
            simulation_time,
        }
    }

    pub fn run(&self) -> () {
        let simulation_time = self.simulation_time;

        for iteration in simulation_time.iter() {
            let (
                hw_demand,
                hw_duration,
                no_events,
                pw_losses_internal,
                gains_internal_dhw_use,
                hw_energy_demand,
            ) = hot_water_demand(iteration.index);

            let gains_internal_dhw = (pw_losses_internal + gains_internal_dhw_use)
                * WATTS_PER_KILOWATT as f64
                / iteration.timestep;
            // TODO: python has clauses to add to this value for particular types of hot water source

            let SpaceHeatingCalculation {
                gains_internal_zone,
                gains_solar_zone,
                operative_temp,
                internal_air_temp,
                space_heat_demand_zone,
                space_cool_demand_zone,
                space_heat_demand_system,
                space_cool_demand_system,
                ductwork_gains,
                heat_balance,
            } = self.calc_space_heating(iteration.timestep, gains_internal_dhw, iteration.index);
        }

        fn hot_water_demand(timestep_idx: usize) -> (f64, f64, u32, f64, f64, f64) {
            (10., 10., 10, 10., 10., 10.)
        }
    }

    /// Calculate the losses/gains in the MVHR ductwork
    ///
    /// ## Arguments
    /// * `timestep_idx` -- timestep index/count
    /// * `delta_t_h` -- calculation timestep, in hours
    /// * `efficiency` -- MVHR heat recovery efficiency
    fn calc_ductwork_losses(&self, timestep_idx: usize, delta_t_h: f64, efficiency: f64) -> f64 {
        // assume 100% efficiency
        // i.e. temp inside the supply and extract ducts is room temp and temp inside exhaust and intake is external temp
        // assume MVHR unit is running 100% of the time
        // let internal_air_temperature = self.temp_internal_air();

        // Calculate heat loss from ducts when unit is inside
        // Air temp inside ducts increases, heat lost from dwelling
        // TODO: complete
        0.0
    }

    /// Calculate space heating demand, heating system output and temperatures
    ///
    /// ## Arguments
    /// * `delta_t_h` - calculation timestep, in hours
    /// * `gains_internal_dhw` - internal gains from hot water system for this timestep, in W
    fn calc_space_heating(
        &self,
        delta_t_h: f64,
        gains_internal_dhw: f64,
        timestep_idx: usize,
    ) -> SpaceHeatingCalculation {
        let temp_ext_air = self
            .external_conditions
            .air_temp_for_timestep_idx(timestep_idx);

        // calculate timestep in seconds
        let delta_t = delta_t_h * SECONDS_PER_HOUR as f64;

        let (mut ductwork_losses, mut ductwork_losses_per_m3) = (0.0, 0.0);

        // ductwork gains/losses only for MVHR
        if let Some(Ventilation::MVHR { .. }) = self.corpus.ventilation {
            // ductwork_losses =
            //     calc_ductwork_losses(0, delta_t_h, self.corpus.ventilation.efficiency());
            // ductwork_losses_per_m3 = ductwork_losses /
        }

        // TODO: complete

        Default::default()
    }
}

#[derive(Default)]
struct SpaceHeatingCalculation {
    gains_internal_zone: f64,
    gains_solar_zone: f64,
    operative_temp: f64,
    internal_air_temp: f64,
    space_heat_demand_zone: f64,
    space_cool_demand_zone: f64,
    space_heat_demand_system: f64,
    space_cool_demand_system: f64,
    ductwork_gains: f64,
    heat_balance: HashMap<String, f64>,
}
