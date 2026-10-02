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
use std::sync::atomic::{AtomicBool, Ordering};
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
    temp_surrounding_prev_heating_event: Vec<AtomicF64>,
    flag_first_pipework_heating_event: AtomicBool,
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
        let temp_surrounding_prev_heating_event: Vec<AtomicF64> = pipework_list
            .iter()
            .map(|pw| {
                AtomicF64::new(Self::temp_surrounding_pipework(
                    pw,
                    temp_external_air_fn.clone(),
                    temp_internal_air_fn.clone(),
                ))
            })
            .collect();

        Self {
            primary_pipework: pipework_list,
            temp_external_air_fn,
            temp_internal_air_fn,
            pipework_energy_input_prev_timestep: AtomicF64::new(0.),
            temp_surrounding_prev_heating_event,
            flag_first_pipework_heating_event: true.into(),
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
                if !self
                    .flag_first_pipework_heating_event
                    .load(Ordering::SeqCst)
                {
                    let inside_temp = self.temp_surrounding_prev_heating_event.get(pipe_idx).ok_or_else(|| anyhow!("Index ({pipe_idx}) out of bounds for temp_surrounding_prev_heating_event"))?.load(Ordering::SeqCst);
                    let between_events_loss =
                        pipework.calculate_cool_down_loss(inside_temp, temp_surrounding);
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

        // Phase 3: End of heating event — record surrounding temps
        if relative_eq!(energy_input, 0., epsilon = 1e-10, max_relative = 1e-9)
            && energy_input_prev > 0.
        {
            for (pipe_idx, pipework) in self.primary_pipework.iter().enumerate() {
                let temp_surrounding = Self::temp_surrounding_pipework(
                    pipework,
                    self.temp_external_air_fn.clone(),
                    self.temp_internal_air_fn.clone(),
                );
                self.temp_surrounding_prev_heating_event
                    .get(pipe_idx).ok_or_else(|| anyhow!("Index ({pipe_idx}) out of bounds for temp_surrounding_prev_heating_event"))?.store(temp_surrounding, Ordering::SeqCst);
                if matches!(pipework.location(), PipeworkLocation::Internal) {
                    primary_gains_w += pipework
                        .calculate_cool_down_loss(temp_flow, temp_surrounding)
                        * WATTS_PER_KILOWATT as f64
                        / timestep;
                }
            }
            if self
                .flag_first_pipework_heating_event
                .load(Ordering::SeqCst)
            {
                self.flag_first_pipework_heating_event
                    .store(false, Ordering::SeqCst);
            }
        }

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

    #[fixture]
    fn external_pipework() -> Pipework {
        Pipework::new(
            PipeworkLocation::External,
            0.025,
            0.027,
            1.5,
            0.035,
            0.038,
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
            results.push(
                mixin
                    .calculate_primary_pipework_losses(3., 55., Some(true), simtime.step)
                    .unwrap(),
            );
        }

        // First timestep: Phase 1 (warm-up) + Phase 2 (steady-state)
        // warm_up(20→55) + steady_state_kWh(55→20)
        assert_relative_eq!(results[0].0, 0.04746228058715814);
        // Subsequent timesteps: Phase 2 only — steady_state_kWh(55→20)
        assert_relative_eq!(results[1].0, 0.010657894331822992);
        assert_relative_eq!(results[2].0, 0.010657894331822992);
    }

    #[rstest]
    /// Internal pipework steady-state losses are also returned as dwelling gains.
    fn test_phase2_steady_state_internal_contributes_gains(internal_pipework: Pipework) {
        let simtime = simtime(2.);
        let mixin = concrete_pipework_user(vec![internal_pipework], Some(20.));

        // Skip first timestep (Phase 1 fires), check second (Phase 2 only)
        let mut results = Vec::new();
        for _ in simtime.iter() {
            results.push(
                mixin
                    .calculate_primary_pipework_losses(3., 55., Some(true), simtime.step)
                    .unwrap(),
            );
        }
        let (_, gains_steady) = results[1];

        // steady_state_heat_loss(55→20) in W
        assert_relative_eq!(gains_steady, 10.657894331822993);
    }

    #[rstest]
    /// External pipework losses do NOT contribute to dwelling heat gains.
    fn test_phase2_external_pipework_no_gains(external_pipework: Pipework) {
        let simtime = simtime(2.);
        let mixin = concrete_pipework_user(vec![external_pipework], Some(5.));

        // Skip first timestep (Phase 1 fires), check second (Phase 2 only)
        let mut results = Vec::new();
        for _ in simtime.iter() {
            results.push(
                mixin
                    .calculate_primary_pipework_losses(3., 55., Some(true), simtime.step)
                    .unwrap(),
            );
        }
        let (losses, gains) = results[1];

        // steady_state_kWh(55→5) for external pipe, no dwelling gains
        assert_relative_eq!(losses, 0.011708048420277326);
        assert_eq!(gains, 0.);
    }

    #[rstest]
    /// End of heating (energy→0 after prev>0) records surrounding temps.
    /// Internal pipework should contribute cool-down gains at the end of the event.
    fn test_phase3_end_of_heating_records_temps_and_reports_internal_gains(
        internal_pipework: Pipework,
    ) {
        let simtime = simtime(3.);
        let mixin = concrete_pipework_user(vec![internal_pipework], Some(20.));

        // Timestep 0: heating active
        let mut results = vec![mixin
            .calculate_primary_pipework_losses(3., 55., Some(true), simtime.step)
            .unwrap()];

        // Timestep 1: heating ends → Phase 3
        results.push(
            mixin
                .calculate_primary_pipework_losses(0., 55., Some(true), simtime.step)
                .unwrap(),
        );

        // Timestep 2: still off
        results.push(
            mixin
                .calculate_primary_pipework_losses(0., 55., Some(true), simtime.step)
                .unwrap(),
        );

        let (losses_end, gains_end) = results[1];

        // Phase 3 cool-down gains: cool_down(55→20) * W_per_kW / timestep
        assert_relative_eq!(gains_end, 36.804386255335146);
        // No losses when energy is zero
        assert_eq!(losses_end, 0.);
        // First heating event flag should now be cleared
        assert!(!mixin
            .flag_first_pipework_heating_event
            .load(Ordering::SeqCst));
    }

    #[rstest]
    /// Between-event losses are only calculated after the first heating event ends.
    /// Sequence: heating on → off (Phase 3, records surrounding temp at 20°C) →
    /// change surrounding temp to 15°C → on again (Phase 1 includes between-event
    /// cool-down from 20→15 on top of warm-up and steady-state).
    fn test_between_events_loss_only_after_first_event_ends(internal_pipework: Pipework) {
        let simtime = simtime(4.);
        let mut mixin = concrete_pipework_user(vec![internal_pipework], Some(20.));

        // Timestep 0: first heating event starts
        let (losses_first_start, _) = mixin
            .calculate_primary_pipework_losses(3., 55., Some(true), simtime.step)
            .unwrap();

        // Timestep 1: first heating event ends (Phase 3 records surrounding=20)
        mixin
            .calculate_primary_pipework_losses(0., 55., Some(true), simtime.step)
            .unwrap();

        // Change surrounding temp so between-event cool-down is non-zero
        mixin.temp_external_air_fn = Arc::new(move || 15.);
        mixin.temp_internal_air_fn = Arc::new(move || 15.);

        // Timestep 2: second heating event starts → should include between-event loss
        let (losses_second_start, _) = mixin
            .calculate_primary_pipework_losses(3., 55., Some(true), simtime.step)
            .unwrap();

        // First start: warm_up(20→55) + steady_state_kWh(55→20), no between-event
        assert_relative_eq!(losses_first_start, 0.04746228058715814);
        // Second start: warm_up(15→55) + between_event(20→15) + steady_state_kWh(55→15)
        assert_eq!(losses_second_start, 0.05950037585037147);
    }

    #[rstest]
    /// Between-event losses on internal pipework contribute to dwelling gains.
    /// With surrounding temp rising from 18°C to 22°C between events, the
    /// between-event cool-down is negative (pipe absorbs heat from the warmer
    /// surroundings). This reduces gains compared to a scenario without
    /// between-event losses, verifying the between-event term is applied.
    fn test_between_events_internal_gain(internal_pipework: Pipework) {
        let simtime = simtime(4.);
        let mut mixin = concrete_pipework_user(vec![internal_pipework], Some(18.));

        // Heat on
        mixin
            .calculate_primary_pipework_losses(3., 55., Some(true), simtime.step)
            .unwrap();

        // Heat off (end of first event, records surrounding=18)
        mixin
            .calculate_primary_pipework_losses(0., 55., Some(true), simtime.step)
            .unwrap();

        // Surrounding rises to 22°C between events
        mixin.temp_external_air_fn = Arc::new(|| 22.);
        mixin.temp_internal_air_fn = Arc::new(|| 22.);

        // Heat on again — Phase 1 between-event + Phase 2 steady-state
        let (_, gains) = mixin
            .calculate_primary_pipework_losses(3., 55., Some(true), simtime.step)
            .unwrap();

        // Total gains = between_event(18→22) * W/kW / timestep + ss(55→22)
        // between_event cool_down(18→22) is negative (pipe absorbs heat from
        // warmer surroundings), so total gains are less than ss(55→22) alone
        assert_relative_eq!(gains, 5.842656226537662);
        // Verify the between-event term specifically reduces gains below
        // what steady-state alone would give (10.049 W)
        let steady_state_only =
            mixin.primary_pipework[0].calculate_steady_state_heat_loss(55., 22.);

        assert!(gains < steady_state_only);
    }

    #[rstest]
    /// Energy within abs_tol=1e-10 of zero triggers Phase 3 (end of heating).
    /// This exercises the math.isclose boundary: a tiny energy_input after a
    /// real heating timestep should be treated as end-of-event.
    fn test_phase3_fires_with_tiny_energy_close_to_zero(internal_pipework: Pipework) {
        let simtime = simtime(3.);
        let mixin = concrete_pipework_user(vec![internal_pipework], Some(20.));

        // Timestep 0: heating active
        mixin
            .calculate_primary_pipework_losses(3., 55., Some(true), simtime.step)
            .unwrap();

        // Timestep 1: tiny energy (effectively zero) → Phase 3 should fire
        let (_, gains) = mixin
            .calculate_primary_pipework_losses(1e-11, 55., Some(true), simtime.step)
            .unwrap();

        // Phase 3 cool-down gains + Phase 2 steady-state gains:
        // warm_up(20→55) * W_per_kW / timestep + ss(55→20)
        assert_relative_eq!(gains, 47.46228058715814);

        // First-event flag should be cleared
        assert!(!mixin
            .flag_first_pipework_heating_event
            .load(Ordering::SeqCst));
    }

    #[rstest]
    /// Both internal and external pipes contribute to losses; only internal to gains.
    fn test_mixed_internal_and_external_pipework(
        internal_pipework: Pipework,
        external_pipework: Pipework,
    ) {
        let simtime = simtime(2.);
        let mixin = concrete_pipework_user(vec![internal_pipework, external_pipework], Some(20.));

        // Skip first timestep (Phase 1), check second (Phase 2 steady-state only)
        let mut results = Vec::new();
        for _ in simtime.iter() {
            results.push(
                mixin
                    .calculate_primary_pipework_losses(3., 55., None, simtime.step)
                    .unwrap(),
            );
        }
        let (losses_both, gains_both) = results[1];

        // Both pipes' steady-state losses combined
        // ss_kWh(int, 55→20) + ss_kWh(ext, 55→20)
        assert_relative_eq!(losses_both, 0.01885352822601712);
        // Only internal pipe's steady-state contributes to gains
        assert_relative_eq!(gains_both, 10.657894331822993);
    }

    #[rstest]
    /// initialising PrimaryPipeworkLossesMixin calls temp_surrounding_pipework for each pipe.
    fn test_init_records_initial_surrounding_temps(internal_pipework: Pipework) {
        let mixin = concrete_pipework_user(vec![internal_pipework], Some(18.));

        assert_eq!(
            mixin.temp_surrounding_prev_heating_event,
            [AtomicF64::new(18.)]
        );
        assert!(mixin
            .flag_first_pipework_heating_event
            .load(Ordering::SeqCst));
        assert_eq!(
            mixin
                .pipework_energy_input_prev_timestep
                .load(Ordering::SeqCst),
            0.
        );
    }
}
