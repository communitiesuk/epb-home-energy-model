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
use crate::core::units::WATTS_PER_KILOWATT;
use crate::corpus::TempInternalAirFn;
use anyhow::anyhow;
use approx::relative_eq;
use atomic_float::AtomicF64;
use std::sync::atomic::Ordering;
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
                Self::temp_surrounding_pipework(
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

    pub(crate) fn temp_surrounding_pipework(
        pipework: &Pipework,
        temp_external_air_fn: Arc<dyn Fn() -> f64 + Send + Sync>,
        temp_internal_air_fn: TempInternalAirFn,
    ) -> f64 {
        match pipework.location() {
            PipeworkLocation::External => temp_external_air_fn(), // TODO WaterPipeworkLocation?
            PipeworkLocation::Internal => temp_internal_air_fn(),
        }
    }

    /// Calculate heat losses through primary pipework.
    ///
    /// Follows the three-step model:
    /// 1. Start of heating event: cool-down loss + between-event losses
    /// 2. During heating event: steady-state conduction losses
    /// 3. End of heating event: record surrounding temps
    ///
    /// Internal pipework losses are returned separately as dwelling heat gains.
    ///
    /// Args:
    ///     energy_input: Energy being delivered this timestep (kWh). Zero
    ///         means no active heating (used to detect event boundaries).
    ///     temp_flow: Flow temperature of water in the pipework (°C).
    ///     update_tracking: Whether to update the internal energy tracking
    ///         state after this call. Set to False for exploratory calls
    ///         where the result may be discarded (e.g. StorageTank calls
    ///         this twice per timestep — exploratory then definitive).
    ///
    /// Returns:
    ///     Tuple of (pipework_losses_kWh, primary_gains_W).
    pub(crate) fn calculate_primary_pipework_losses(
        &self,
        energy_input: f64,
        temp_flow: f64,
        update_tracking: Option<bool>,
        timestep: f64,
    ) -> anyhow::Result<(f64, f64)> {
        let update_tracking = update_tracking.unwrap_or(true);

        if self.primary_pipework.is_empty() {
            return Ok((0., 0.));
        }

        let mut pipework_losses_kwh = 0.;
        let mut primary_gains_w = 0.;
        let energy_input_prev = self
            .pipework_energy_input_prev_timestep
            .load(Ordering::SeqCst);

        // Phase 1: Start of heating event — pipes are cold and need filling
        if energy_input > 0.
            && relative_eq!(energy_input_prev, 0., epsilon = 1e-10, max_relative = 1e-9)
        {
            for (pipe_idx, pipework) in self.primary_pipework.iter().enumerate() {
                let temp_surrounding = Self::temp_surrounding_pipework(
                    pipework,
                    self.temp_external_air_fn.clone(),
                    self.temp_internal_air_fn.clone(),
                );
                let cool_down_loss = pipework.calculate_cool_down_loss(temp_flow, temp_surrounding);
                pipework_losses_kwh += cool_down_loss;

                // Between-event losses: pipe cooled from previous event's
                // surrounding temp to current surrounding temp
                if !self.flag_first_pipework_heating_event {
                    let inside_temp = self.temp_surrounding_prev_heating_event.get(pipe_idx).ok_or_else(|| anyhow!("Index ({pipe_idx}) out of bounds for temp_surrounding_prev_heating_event"))?;
                    let between_events_loss =
                        pipework.calculate_cool_down_loss(*inside_temp, temp_surrounding);
                    pipework_losses_kwh += between_events_loss;
                    if matches!(pipework.location(), PipeworkLocation::Internal) {
                        primary_gains_w +=
                            between_events_loss * WATTS_PER_KILOWATT as f64 / timestep;
                    }
                }
            }
        }

        // Phase 2: During heating event — steady-state conduction losses
        if energy_input > 0. {
            for pipework in self.primary_pipework.iter() {
                let temp_surrounding = Self::temp_surrounding_pipework(
                    pipework,
                    self.temp_external_air_fn.clone(),
                    self.temp_internal_air_fn.clone(),
                );
                let steady_state_loss_w =
                    pipework.calculate_steady_state_heat_loss(temp_flow, temp_surrounding);
                if matches!(pipework.location(), PipeworkLocation::Internal) {
                    primary_gains_w += steady_state_loss_w;
                }
                pipework_losses_kwh += steady_state_loss_w * timestep / WATTS_PER_KILOWATT as f64;
            }
        }

        // TODO complete function

        if update_tracking {
            self.pipework_energy_input_prev_timestep
                .store(energy_input, Ordering::SeqCst);
        }

        Ok((pipework_losses_kwh, primary_gains_w))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::pipework::Pipework;
    use crate::hem_core::simulation_time::SimulationTime;
    use crate::input::PipeworkContents;
    use approx::assert_relative_eq;
    use rstest::{fixture, rstest};

    fn simtime(n_timesteps: f64) -> SimulationTime {
        SimulationTime::new(0., n_timesteps, 1.)
    }

    fn concrete_pipework_user(
        pipework_list: Vec<Pipework>,
        surrounding_temp: Option<f64>,
    ) -> PrimaryPipeworkLossesMixin {
        let surrounding_temp = surrounding_temp.unwrap_or(20.);

        PrimaryPipeworkLossesMixin::new(
            pipework_list,
            Arc::new(move || surrounding_temp),
            Arc::new(move || surrounding_temp),
        )
    }

    #[fixture]
    fn internal_pipework() -> Pipework {
        Pipework::new(
            PipeworkLocation::Internal,
            0.024,
            0.027,
            2.,
            0.035,
            0.04,
            false,
            PipeworkContents::Water,
        )
        .unwrap()
    }

    #[rstest]
    /// No pipework → zero losses and zero gains regardless of energy input.
    fn test_empty_pipework_list_returns_zero() {
        let simtime = simtime(1.);
        let mixin = concrete_pipework_user(vec![], None);
        let mut losses = 0.;
        let mut gains = 0.;

        for _ in simtime.iter() {
            (losses, gains) = mixin
                .calculate_primary_pipework_losses(5., 55., Some(true), simtime.step)
                .unwrap();
        }

        assert_eq!(losses, 0.);
        assert_eq!(gains, 0.);
    }

    #[rstest]
    /// No energy input → zero losses and zero gains (no active heating).
    fn test_zero_energy_returns_zero(internal_pipework: Pipework) {
        let simtime = simtime(1.);
        let mixin = concrete_pipework_user(vec![internal_pipework], None);
        let mut losses = 0.;
        let mut gains = 0.;

        for _ in simtime.iter() {
            (losses, gains) = mixin
                .calculate_primary_pipework_losses(0., 55., Some(true), simtime.step)
                .unwrap();
        }

        assert_eq!(losses, 0.);
        assert_eq!(gains, 0.);
    }

    #[rstest]
    /// First timestep with energy > 0 after zero triggers Phase 1 warm-up loss.
    /// Phase 1 (warm-up) fires once when transitioning from energy_input_prev=0
    /// to energy_input>0, adding warm-up loss on top of steady-state. Subsequent
    /// timesteps are Phase 2 only (steady-state).
    fn test_phase1_start_of_heating_adds_warm_up_loss(internal_pipework: Pipework) {
        let simtime = simtime(3.);
        let mixin = concrete_pipework_user(vec![internal_pipework], Some(20.));

        let mut results = Vec::new();
        for _ in simtime.iter() {
            let result = mixin
                .calculate_primary_pipework_losses(3., 55., Some(true), simtime.step)
                .unwrap();
            results.push(result);
        }

        // First timestep: Phase 1 (warm-up) + Phase 2 (steady-state)
        // warm_up(20→55) + steady_state_kWh(55→20)
        assert_relative_eq!(results[0].0, 0.04746228058715814);
        // Subsequent timesteps: Phase 2 only — steady_state_kWh(55→20)
        assert_relative_eq!(results[1].0, 0.010657894331822992);
        assert_relative_eq!(results[2].0, 0.010657894331822992);
    }
}
