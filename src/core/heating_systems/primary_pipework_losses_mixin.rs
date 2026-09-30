/// Mixin for three-step primary pipework heat loss calculations.
///
/// Provides shared logic for calculating heat losses through primary pipework
/// between a heat source and a storage vessel or heat battery. Both StorageTank
/// and HeatBatteryPCM use this three-step model:
/// 1. Start of heating event: cool-down loss to fill cold pipes + between-event losses
/// 2. During heating event: steady-state conduction losses
/// 3. End of heating event: record surrounding temps for next event's between-event calc
///
/// The mixin resolves surrounding temperatures via two callbacks supplied at
/// initialisation — one for external pipework (outdoor air temperature) and one
/// for internal pipework (zone air temperature). This keeps the mixin
/// self-contained while allowing each subclass to provide its own temperature
/// sources.
use crate::core::pipework::{Pipework, PipeworkLocation, Pipeworkesque};
use crate::corpus::TempInternalAirFn;
use atomic_float::AtomicF64;
use std::sync::Arc;

/// Shared three-step primary pipework loss calculation.
///
/// Temperature lookup is handled internally via two callbacks passed to
/// init_pipework_state: one for external pipework (outdoor air) and one
/// for internal pipework (zone air).
///
/// State variables:
/// pipework_energy_input_prev_timestep: Energy input from the previous
/// timestep, used to detect start/end of heating events.
/// temp_surrounding_prev_heating_event: Surrounding temperature at each
/// pipe segment when the previous heating event ended. Used for
/// between-event cool-down loss calculation.
/// flag_first_pipework_heating_event: True until the first heating event
/// completes. Between-event losses are not calculated before the
/// first event ends.
struct PrimaryPipeworkLossesMixin {
    primary_pipework: Vec<Pipework>,
    pipework_energy_input_prev_timestep: AtomicF64,
    temp_surrounding_prev_heating_event: Vec<f64>,
    flag_first_pipework_heating_event: bool,
    temp_external_air_fn: Arc<dyn Fn() -> f64 + Send + Sync>, // TODO review type
    temp_internal_air_fn: TempInternalAirFn,
}

impl PrimaryPipeworkLossesMixin {
    /// Arguments
    /// * `pipework_list` - List of Pipework objects for primary circuit
    /// * `temp_external_air_fn` - Returns the current outdoor air temperature (°C) for external pipework segments
    /// * `temp_internal_air_fn` - Returns the current internal air temperature (°C) for internal pipework segments
    pub(crate) fn new(
        pipework_list: Vec<Pipework>,
        temp_external_air_fn: Arc<dyn Fn() -> f64 + Send + Sync>,
        temp_internal_air_fn: TempInternalAirFn,
    ) -> Self {
        let temp_surrounding_prev_heating_event: Vec<f64> = pipework_list
            .iter()
            .map(|pw| {
                Self::get_temp_surrounding_pipework(
                    pw,
                    temp_external_air_fn.clone(),
                    temp_internal_air_fn.clone(),
                )
            })
            .collect();

        Self {
            primary_pipework: pipework_list,
            temp_external_air_fn,
            temp_internal_air_fn,
            pipework_energy_input_prev_timestep: AtomicF64::new(0.),
            temp_surrounding_prev_heating_event,
            flag_first_pipework_heating_event: Default::default(),
        }
    }

    pub(crate) fn get_temp_surrounding_pipework(
        pipework: &Pipework,
        temp_external_air_fn: Arc<dyn Fn() -> f64 + Send + Sync>,
        temp_internal_air_fn: TempInternalAirFn,
    ) -> f64 {
        match pipework.location() {
            PipeworkLocation::External => temp_external_air_fn(),
            PipeworkLocation::Internal => temp_internal_air_fn(),
        }
    }
}
