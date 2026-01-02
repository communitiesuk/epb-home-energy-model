use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::heating_systems::wwhrs::{WwhrsConfiguration, WwhrsInstantaneous};
use crate::core::material_properties::WATER;
use crate::core::units::MINUTES_PER_HOUR;
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::misc::volume_hot_water_required;
use crate::core::water_heat_demand::misc::CallableGetHotWaterTemperature;
use crate::input::WaterHeatingEvent;
use crate::simulation_time::SimulationTimeIteration;
use anyhow::anyhow;
use parking_lot::Mutex;
use std::sync::Arc;

#[derive(Debug)]
pub(crate) enum Shower {
    MixerShower(MixerShower),
    InstantElectricShower(InstantElectricShower),
}

impl Shower {
    pub(crate) fn get_cold_water_source(&self) -> &ColdWaterSource {
        match self {
            Shower::MixerShower(s) => s.get_cold_water_source(),
            Shower::InstantElectricShower(s) => s.get_cold_water_source(),
        }
    }

    pub(crate) fn hot_water_demand(
        &self,
        temp_target: f64,
        temp_hot_water: f64,
        total_shower_duration: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        match self {
            Shower::MixerShower(_s) => {
                todo!("as part of migration to 1.0.0a1");
                // s.hot_water_demand(temp_target, temp_hot_water, total_shower_duration, simtime)
            }
            Shower::InstantElectricShower(s) => {
                todo!("as part of migration to 1.0.0a1");
                // s.hot_water_demand(temp_target, temp_hot_water, total_shower_duration, simtime)
            }
        }
    }
}

#[derive(Debug)]
/// An object to model mixer showers i.e. those that mix hot and cold water
pub struct MixerShower {
    flowrate: f64,
    cold_water_source: Arc<ColdWaterSource>,
    wwhrs: Option<Arc<Mutex<WwhrsInstantaneous>>>,
    wwhrs_configuration: WwhrsConfiguration,
}

impl MixerShower {
    /// Construct a MixerShower object
    ///
    /// Arguments:
    /// * `flowrate` - shower's flow rate, in litres/minute
    /// * `cold_water_source` - reference to ColdWaterSource object representing the
    ///                     cold water feed attached to the shower
    /// * `wwhrs` - reference to WWHRS object (optional)
    /// * `wwhrs_configuration` - WWHRS system configuration ('A', 'B', or 'C') for this shower
    pub(crate) fn new(
        flowrate: f64,
        cold_water_source: Arc<ColdWaterSource>,
        wwhrs: Option<Arc<Mutex<WwhrsInstantaneous>>>,
        wwhrs_configuration: Option<WwhrsConfiguration>,
    ) -> Self {
        Self {
            flowrate,
            cold_water_source,
            wwhrs,
            wwhrs_configuration: wwhrs_configuration
                .unwrap_or(WwhrsConfiguration::AShowerAndWaterHeatingSystem),
        }
    }

    pub(crate) fn get_cold_water_source(&self) -> &ColdWaterSource {
        &self.cold_water_source
    }

    /// Calculate volume of hot water required
    /// (and volume of warm water draining to WWHRS, if applicable)
    ///
    ///
    /// Arguments:
    /// * `event` - dict containing "temperature" and "duration" keys
    ///           - temperature: temperature of warm water delivered at shower head, in Celsius
    ///           - duration: cumulative running time of this shower during the current timestep, in minutes
    /// * `func_temp_hot_water` - callable
    pub(crate) fn hot_water_demand(
        &self,
        event: WaterHeatingEvent,
        func_temp_hot_water: &CallableGetHotWaterTemperature,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(Option<f64>, f64)> {
        let total_shower_duration = event
            .duration
            .ok_or(anyhow!("Expected an event duration for Mixer Shower event"))?;
        let temperature_target = event.temperature;

        // TODO (Python) Account for behavioural variation factor fbeh
        let volume_warm_water = self.flowrate * total_shower_duration;
        // ^^^ litres = litres/minute * minutes
        let volume_hot_water = volume_hot_water_required(
            volume_warm_water,
            temperature_target,
            func_temp_hot_water,
            |x, simtime| self.cold_water_source.get_temp_cold_water(x, simtime),
            simtime,
        )?;
        // first calculate the volume of hot water needed if heating from cold water source

        let mut volume_hot_water = match volume_hot_water {
            None => return Ok((None, volume_warm_water)),
            Some(v) => v,
        };

        let volume_cold_water = volume_warm_water - volume_hot_water;

        if let Some(wwhrs) = &self.wwhrs {
            let temperature_hot_water = func_temp_hot_water(volume_hot_water); // temperature of hot water supply, in Celsius

            //  Use the unified WWHRS interface
            let wwhrs_result = wwhrs.lock().calculate_performance(
                self.wwhrs_configuration,
                temperature_target,
                self.flowrate,
                volume_cold_water,
                temperature_hot_water,
                simtime,
            )?;
            let wwhrs_cyl_feed_temperature = wwhrs_result.t_cyl_feed;
            volume_hot_water = wwhrs_result.flowrate_hot.ok_or_else(|| {
                anyhow!("Flowrate hot water not available from WWHRS calculation")
            })? * total_shower_duration;

            // Set the actual return temperature given the temperature and flowrate of the waste water.
            if !matches!(self.wwhrs_configuration, WwhrsConfiguration::BShower) {
                if let Some(wwhrs) = self.wwhrs.as_ref() {
                    wwhrs
                        .lock()
                        .set_temperature_for_return(wwhrs_cyl_feed_temperature);
                }
            }
        }

        let volume_cold_water = volume_warm_water - volume_hot_water;

        match (self.wwhrs.as_ref(), self.wwhrs_configuration) {
            (
                Some(wwhrs),
                WwhrsConfiguration::AShowerAndWaterHeatingSystem
                | WwhrsConfiguration::CWaterHeatingSystem,
            ) => {
                wwhrs.lock().draw_off_water(volume_cold_water);
            }
            _ => {
                self.cold_water_source
                    .draw_off_water(volume_cold_water, simtime)?;
            }
        }

        Ok((volume_hot_water.into(), volume_warm_water))
    }
}

#[derive(Debug)]
pub struct InstantElectricShower {
    power_in_kilowatts: f64,
    cold_water_source: Arc<ColdWaterSource>,
    energy_supply_connection: EnergySupplyConnection,
}

impl InstantElectricShower {
    pub(crate) fn new(
        power_in_kilowatts: f64,
        cold_water_source: Arc<ColdWaterSource>,
        energy_supply_connection: EnergySupplyConnection,
    ) -> Self {
        Self {
            power_in_kilowatts,
            cold_water_source,
            energy_supply_connection,
        }
    }

    pub(crate) fn get_cold_water_source(&self) -> &ColdWaterSource {
        &self.cold_water_source
    }

    /// Calculate electrical energy required
    /// (and volume of warm water draining to WWHRS, if applicable)
    pub(crate) fn hot_water_demand(
        &self,
        event: WaterHeatingEvent,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        let total_shower_duration = event
            .duration
            .ok_or(anyhow!("Expected an event duration for shower"))?;
        let temperature_target = event.temperature;

        let electricity_demand =
            self.power_in_kilowatts * (total_shower_duration / MINUTES_PER_HOUR as f64);
        let mut volume_warm_water = total_shower_duration * 6.0;
        let mut volume_warm_water_prev = 0.0;

        while !is_close!(volume_warm_water, volume_warm_water_prev, abs_tol = 1e-10) {
            let list_temperature_volume = self
                .cold_water_source
                .get_temp_cold_water(volume_warm_water, simtime)?;
            let temperature_cold = list_temperature_volume
                .iter()
                .map(|(t, v)| t * v)
                .sum::<f64>()
                / list_temperature_volume.iter().map(|(_, v)| v).sum::<f64>();
            volume_warm_water_prev = volume_warm_water;
            volume_warm_water = electricity_demand
                / WATER
                    .volumetric_energy_content_kwh_per_litre(temperature_target, temperature_cold);
        }

        self.energy_supply_connection
            .demand_energy(electricity_demand, simtime.index)?;

        // Instantaneous electric shower heats its own water, so no demand on
        // the water heating system.
        Ok((0.0, volume_warm_water))
        // TODO (from Python) Should this return hot water demand or send message to HW system?
        //      The latter would allow for different showers to be connected to
        //      different HW systems, but complicates the implementation of the
        //      HW system as it will have to deal with calls from several
        //      different objects and work out when to amalgamate the figures to
        //      do its own calculation, rather than being given a single overall
        //      figure for each timestep.
        // TODO (from Python) Also send volume_warm_water to connected WWHRS object? Account for
        //      heat loss between shower head and drain?
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    mod mixer_shower {
        use super::*;
        use crate::simulation_time::SimulationTime;
        use approx::assert_relative_eq;
        use itertools::Itertools;
        use parking_lot::Mutex;
        use pretty_assertions::assert_eq;
        use rstest::*;

        #[fixture]
        fn simulation_time() -> SimulationTime {
            SimulationTime::new(0f64, 3f64, 1f64)
        }

        #[fixture]
        fn cold_water_source() -> ColdWaterSource {
            ColdWaterSource::new(vec![2.0, 3.0, 4.0], 0, 1.0)
        }

        #[fixture]
        fn mixer_shower(cold_water_source: ColdWaterSource) -> MixerShower {
            MixerShower::new(6.5, cold_water_source.into(), None, None)
        }

        fn func_temp_hot_water_fixed(_f: f64) -> f64 {
            52.0
        }

        #[fixture]
        fn func_temp_hot_water() -> CallableGetHotWaterTemperature {
            Box::new(func_temp_hot_water_fixed)
        }

        #[rstest]
        fn test_hot_water_demand(
            simulation_time: SimulationTime,
            mixer_shower: MixerShower,
            func_temp_hot_water: CallableGetHotWaterTemperature,
        ) {
            let expected_demands = [24.7, 24.54081632653061, 24.375];
            for (idx, t_it) in simulation_time.iter().enumerate() {
                let event = WaterHeatingEvent {
                    start: 0.,
                    duration: Some(5.0),
                    temperature: 40.0,
                    volume: None,
                };
                let result = mixer_shower
                    .hot_water_demand(event, &func_temp_hot_water, t_it)
                    .unwrap()
                    .0
                    .unwrap();
                assert_eq!(
                    result, expected_demands[idx],
                    "incorrect volume of hot water returned"
                );
            }
        }

        // upstream Python test_get_cold_water_source test here is redundant as it only checks for type

        #[rstest]
        fn test_vol_warm_water_for_mixer(
            simulation_time: SimulationTime,
            mixer_shower: MixerShower,
            func_temp_hot_water: CallableGetHotWaterTemperature,
        ) {
            let expected_volumes = [32.5, 32.5, 32.5];
            let event = WaterHeatingEvent {
                start: 0.,
                duration: Some(5.),
                temperature: 40.0,
                volume: None,
            };
            for (idx, t_it) in simulation_time.iter().enumerate() {
                assert_eq!(
                    mixer_shower
                        .hot_water_demand(event, &func_temp_hot_water, t_it)
                        .unwrap()
                        .1,
                    expected_volumes[idx],
                    "incorrect volume of warm water returned"
                );
            }
        }

        #[rstest]
        fn test_wwhrs_instantaneous_system_b(
            simulation_time: SimulationTime,
            cold_water_source: ColdWaterSource,
            func_temp_hot_water: CallableGetHotWaterTemperature,
        ) {
            let cold_water_source = Arc::new(cold_water_source);
            let flow_rates = vec![5., 7., 9., 11., 13.];
            let system_a_efficiencies = vec![44.8, 39.1, 34.8, 31.4, 28.6];
            let wwhrs = Arc::new(Mutex::new(
                WwhrsInstantaneous::new(
                    flow_rates,
                    system_a_efficiencies,
                    cold_water_source.clone(),
                    Some(0.7),
                    None,
                    Some(0.7),
                    None,
                    None,
                    Some(0.81),
                    None,
                    simulation_time.iter().current_iteration(),
                )
                .unwrap(),
            ));
            let mixer_shower = MixerShower::new(
                6.5,
                cold_water_source,
                wwhrs.into(),
                WwhrsConfiguration::BShower.into(),
            );

            // With Technical Recommendations calculations for System B
            // temp_shower = 40°C, temp_hot = 52°C, flowrate = 6.5 L/min
            // For each timestep with different cold water temps
            let cold_water_temps = [2.0, 3.0, 4.0];

            let expected_hot_water = (0..3)
                .map(|t_idx| {
                    let temp_cold = cold_water_temps[t_idx];

                    // Interpolate efficiency at 6.5 L/min ≈ 40.525%
                    // With System B factor (0.81) and utilisation factor (0.7): 40.525 * 0.81 * 0.7 = 22.977675%
                    // T_drain = 40 - 6 = 34°C
                    // temp = 0.22977675 * (34 - temp_cold) / (52 - 40)
                    // T_pre = (temp_cold + 52 * temp) / (1 + temp)

                    let temp_drain = 34.0;
                    let eta_uf = 0.40525 * 0.81 * 0.7; // ≈ 0.22977675
                    let temp = eta_uf * (temp_drain - temp_cold) / (52.0 - 40.0);
                    let temp_pre = (temp_cold + 52.0 * temp) / (1. + temp);

                    // Volume of hot water = warm_water * (temp_shower - temp_pre) / (temp_hot - temp_pre)
                    let vol_warm = 32.5;

                    vol_warm * (40.0 - temp_pre) / (52.0 - temp_pre)
                })
                .collect_vec();

            for (t_idx, t_it) in simulation_time.iter().enumerate() {
                assert_relative_eq!(
                    mixer_shower
                        .hot_water_demand(
                            WaterHeatingEvent {
                                temperature: 40.0,
                                duration: Some(5.0),
                                start: 0.,
                                volume: None,
                            },
                            &func_temp_hot_water,
                            t_it
                        )
                        .unwrap()
                        .0
                        .unwrap(),
                    expected_hot_water[t_idx],
                    epsilon = 1e-3
                );
            }
        }

        #[rstest]
        fn test_wwhrs_instantaneous_system_c(
            simulation_time: SimulationTime,
            mut mixer_shower: MixerShower,
            func_temp_hot_water: CallableGetHotWaterTemperature,
        ) {
            let flow_rates = vec![5., 7., 9., 11., 13.];
            let system_a_efficiencies = vec![44.8, 39.1, 34.8, 31.4, 28.6];
            let wwhrs = Arc::new(Mutex::new(
                WwhrsInstantaneous::new(
                    flow_rates,
                    system_a_efficiencies,
                    mixer_shower.cold_water_source.clone(),
                    Some(0.7),
                    None,
                    None,
                    None,
                    Some(0.7),
                    None,
                    Some(0.88),
                    simulation_time.iter().current_iteration(),
                )
                .unwrap(),
            ));
            mixer_shower.wwhrs = Some(wwhrs);
            mixer_shower.wwhrs_configuration = WwhrsConfiguration::CWaterHeatingSystem;

            let expected_volumes = [24.7, 24.54081632653061, 24.375];

            for (idx, t_it) in simulation_time.iter().enumerate() {
                assert_relative_eq!(
                    mixer_shower
                        .hot_water_demand(
                            WaterHeatingEvent {
                                temperature: 40.,
                                duration: Some(5.),
                                start: 0.,
                                volume: None
                            },
                            &func_temp_hot_water,
                            t_it
                        )
                        .unwrap()
                        .0
                        .unwrap(),
                    expected_volumes[idx],
                    epsilon = 1e-6
                );
            }
        }

        #[rstest]
        fn test_wwhrs_instantaneous_system_a(
            simulation_time: SimulationTime,
            mut mixer_shower: MixerShower,
            func_temp_hot_water: CallableGetHotWaterTemperature,
        ) {
            let flow_rates = vec![5., 7., 9., 11., 13.];
            let system_a_efficiencies = vec![44.8, 39.1, 34.8, 31.4, 28.6];
            let wwhrs = Arc::new(Mutex::new(
                WwhrsInstantaneous::new(
                    flow_rates,
                    system_a_efficiencies,
                    mixer_shower.cold_water_source.clone(),
                    Some(0.7),
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    simulation_time.iter().current_iteration(),
                )
                .unwrap(),
            ));
            mixer_shower.wwhrs = Some(wwhrs);

            let cold_water_temps = [2.0, 3.0, 4.0];

            let expected_hot_water = (0..3)
                .map(|t_idx| {
                    let temp_cold = cold_water_temps[t_idx];

                    // Interpolate efficiency at 6.5 L/min ≈ 40.525%
                    // With utilisation factor (0.7): 40.525 * 0.7 = 28.3675%
                    // T_drain = 40 - 6 = 34°C
                    // T_pre = temp_cold + 0.287 * (34 - temp_cold)

                    let temp_drain = 34.0;
                    let eta_uf = 0.40525 * 0.7; //≈ 0.283675
                    let temp_pre = temp_cold + eta_uf * (temp_drain - temp_cold);

                    // Volume of hot water = warm_water * (temp_shower - temp_pre) / (temp_hot - temp_pre)
                    let vol_warm = 32.5;
                    vol_warm * (40.0 - temp_pre) / (52.0 - temp_pre)
                })
                .collect_vec();

            for (t_idx, t_it) in simulation_time.iter().enumerate() {
                assert_relative_eq!(
                    mixer_shower
                        .hot_water_demand(
                            WaterHeatingEvent {
                                temperature: 40.0,
                                duration: Some(5.0),
                                start: 0.,
                                volume: None,
                            },
                            &func_temp_hot_water,
                            t_it
                        )
                        .unwrap()
                        .0
                        .unwrap(),
                    expected_hot_water[t_idx],
                    epsilon = 1e-3
                );
            }
        }
    }

    mod instant_elec_shower {
        use super::*;
        use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyBuilder};
        use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
        use crate::hem_core::simulation_time::SimulationTime;
        use crate::input::{FuelType, WaterHeatingEvent};
        use parking_lot::RwLock;
        use rstest::rstest;
        use std::sync::Arc;

        #[rstest]
        fn test_hot_water_demand() {
            let simulation_time = SimulationTime::new(0f64, 3f64, 1f64);
            let cold_water_temps = [2.0, 3.0, 4.0];
            let cold_water_source = ColdWaterSource::new(cold_water_temps.into(), 0, 1.0);
            let energy_supply = Arc::new(RwLock::new(
                EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps())
                    .build(),
            ));
            let energy_supply_conn =
                EnergySupply::connection(energy_supply.clone(), "shower").unwrap();
            let instant_shower =
                InstantElectricShower::new(50.0, cold_water_source.into(), energy_supply_conn);
            let expected_results_by_end_user = [5.0, 10.0, 15.0];
            let expected_demands = [0.0, 0.0, 0.0];
            for (idx, t_it) in simulation_time.iter().enumerate() {
                let event = WaterHeatingEvent {
                    start: 0.,
                    duration: Some(((idx + 1) * 6) as f64),
                    temperature: 40.0,
                    volume: None,
                };
                let _ = instant_shower.hot_water_demand(event, t_it);
                pretty_assertions::assert_eq!(
                    energy_supply.read().results_by_end_user()["shower"][idx],
                    expected_results_by_end_user[idx],
                    "correct electricity demand not returned"
                );
                pretty_assertions::assert_eq!(
                    instant_shower.hot_water_demand(event, t_it).unwrap().0,
                    expected_demands[idx]
                );
            }
        }
    }
}
