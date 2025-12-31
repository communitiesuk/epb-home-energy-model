use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::heating_systems::wwhrs::{WwhrsConfiguration, WwhrsInstantaneous};
use crate::core::material_properties::WATER;
use crate::core::schedule::TypedScheduleEvent;
use crate::core::units::MINUTES_PER_HOUR;
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::misc::volume_hot_water_required;
use crate::core::water_heat_demand::misc::{
    calc_fraction_hot_water, CallableGetHotWaterTemperature,
};
use crate::simulation_time::SimulationTimeIteration;
use anyhow::anyhow;
use parking_lot::Mutex;
use std::sync::Arc;

#[derive(Debug)]
pub enum Shower {
    MixerShower(MixerShower),
    InstantElectricShower(InstantElectricShower),
}

impl Shower {
    pub fn get_cold_water_source(&self) -> &ColdWaterSource {
        match self {
            Shower::MixerShower(s) => s.get_cold_water_source(),
            Shower::InstantElectricShower(s) => s.get_cold_water_source(),
        }
    }

    pub fn hot_water_demand(
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
                s.hot_water_demand(temp_target, temp_hot_water, total_shower_duration, simtime)
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
    wwhrs_configuration: Option<WwhrsConfiguration>,
}

impl MixerShower {
    /// Construct a MixerShower object
    ///
    /// Arguments:
    /// * `flowrate` - shower's flow rate, in litres/minute
    /// * `cold_water_source` - reference to ColdWaterSource object representing the
    ///                     cold water feed attached to the shower
    /// * `wwhrs` - reference to WWHRS object (optional)
    /// *`wwhrs_configuration` - WWHRS system configuration ('A', 'B', or 'C') for this shower
    pub(crate) fn new(
        flowrate: f64,
        cold_water_source: Arc<ColdWaterSource>,
        wwhrs: Option<Arc<Mutex<WwhrsInstantaneous>>>,
        wwhrs_configuration: Option<WwhrsConfiguration>, // TODO: migration to 1.0.0a1 can only be A_SHOWER_AND_WATER_HEATING_SYSTEM variant, add check somewhere?
    ) -> Self {
        Self {
            flowrate,
            cold_water_source,
            wwhrs,
            wwhrs_configuration,
        }
    }

    pub fn get_cold_water_source(&self) -> &ColdWaterSource {
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
    pub fn hot_water_demand(
        &self,
        event: TypedScheduleEvent,
        func_temp_hot_water: CallableGetHotWaterTemperature,
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

        if volume_hot_water.is_none() {
            // Unmet demand
            return Ok((None, volume_warm_water));
        }

        let _volume_cold_water = volume_warm_water - volume_hot_water.unwrap();

        if let Some(_wwhrs) = &self.wwhrs {
            // let temperature_hot_water = func_temp_hot_water(
            //     volume_hot_water
            // );  // temperature of hot water supply, in Celsius
            //
            // //  Use the unified WWHRS interface
            // let wwhrs_configuration = self.wwhrs_configuration.ok_or(|| anyhow!("WWHRS configuration expeected when WWHRS is provided"))?;
            // let wwhrs_result = wwhrs.lock().calculate_performance(wwhrs_configuration, temperature_target, self.flowrate, volume_cold_water, temperature_hot_water, simtime.iter().current_iteration());
            todo!("migratin to 1.0.0a1");
        }
        Ok((volume_hot_water, volume_warm_water))
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

    pub fn get_cold_water_source(&self) -> &ColdWaterSource {
        &self.cold_water_source
    }

    /// Calculate electrical energy required
    /// (and volume of warm water draining to WWHRS, if applicable)
    ///
    /// Arguments:
    /// * `temp_target` - temperature of warm water delivered at shower head, in Celcius
    /// * `total_shower_duration` - cumulative running time of this shower during
    ///                             the current timestep, in minutes
    pub fn hot_water_demand(
        &self,
        temp_target: f64,
        temp_hot_water: f64,
        total_shower_duration: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        let temp_cold = self.cold_water_source.temperature(simtime);

        let elec_demand =
            self.power_in_kilowatts * (total_shower_duration / MINUTES_PER_HOUR as f64);
        let vol_warm_water =
            elec_demand / WATER.volumetric_energy_content_kwh_per_litre(temp_target, temp_cold);

        let vol_hot_water_equiv =
            vol_warm_water * calc_fraction_hot_water(temp_target, temp_hot_water, temp_cold)?;

        self.energy_supply_connection
            .demand_energy(elec_demand, simtime.index)
            .unwrap();

        Ok((vol_hot_water_equiv, vol_warm_water))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // use crate::core::schedule::WaterScheduleEventType;
    use crate::simulation_time::SimulationTime;
    // use pretty_assertions::assert_eq;
    use rstest::*;
    // use std::sync::Arc;

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0f64, 3f64, 1f64)
    }

    #[fixture]
    fn mixer_shower() -> MixerShower {
        let cold_water_temps = [2.0, 3.0, 4.0];
        let cold_water_source = ColdWaterSource::new(cold_water_temps.into(), 0, 1.0);
        MixerShower::new(6.5, cold_water_source.into(), None, None)
    }

    // #[rstest]
    // fn test_hot_water_demand(simulation_time: SimulationTime, mixer_shower: MixerShower) {
    //     let expected_demands = [24.7, 24.54081632653061, 24.375];
    //     for (idx, t_it) in simulation_time.iter().enumerate() {
    //         let event = TypedScheduleEvent {
    //             start: 0.,
    //             duration: Some(((idx + 1) * 6) as f64),
    //             temperature: 40.0,
    //             name: "instantelec".to_string(),
    //             event_type: WaterScheduleEventType::Shower,
    //             volume: None,
    //             warm_volume: None,
    //             pipework_volume: None,
    //         };
    //         let result = mixer_shower.hot_water_demand(event, None, t_it).unwrap();
    //         assert_eq!(
    //             mixer_shower.hot_water_demand(event, None, t_it).unwrap().0.unwrap(),
    //             expected_demands[idx],
    //             "incorrect volume of hot water returned"
    //         );
    //     }
    // }

    // #[rstest]
    // fn test_vol_warm_water_for_mixer(simulation_time: SimulationTime, mixer_shower: MixerShower) {
    //     let expected_volumes = [32.5, 32.5, 32.5];
    //     let event = TypedScheduleEvent {
    //         start: 0.,
    //         duration: Some(5.),
    //         temperature: 40.0,
    //         name: "mixer".to_string(),
    //         event_type: WaterScheduleEventType::Shower,
    //         volume: None,
    //         warm_volume: None,
    //         pipework_volume: None,
    //     };
    //     for (idx, t_it) in simulation_time.iter().enumerate() {
    //         assert_eq!(
    //             mixer_shower.hot_water_demand(event, None, t_it).unwrap().1,
    //             expected_volumes[idx],
    //             "incorrect volume of warm water returned"
    //         );
    //     }
    // }

    // #[rstest]
    // fn test_wwhrs_instantaneous_system_b_for_mixer(
    //     simulation_time: SimulationTime,
    //     mut mixer_shower: MixerShower,
    // ) {
    //     let flow_rates = vec![5., 7., 9., 11., 13.];
    //     let efficiencies = vec![44.8, 39.1, 34.8, 31.4, 28.6];
    //     let wwhrs = Arc::new(Mutex::new(WwhrsInstantaneous::new(
    //             mixer_shower.cold_water_source.clone(),
    //             flow_rates,
    //             efficiencies,
    //             0.7,
    //     )));
    //     mixer_shower.wwhrs = Some(wwhrs);
    //     let expected_volumes = [22.903242227702766, 22.731048233573134, 22.552562007290966];
    //
    //     for (idx, t_it) in simulation_time.iter().enumerate() {
    //         assert_eq!(
    //             mixer_shower.hot_water_demand(40., 52., 5., t_it).0,
    //             expected_volumes[idx],
    //             "incorrect volume of hot water returned"
    //         );
    //     }
    // }

    // #[rstest]
    // fn test_wwhrs_instantaneous_system_c_for_mixer(
    //     simulation_time: SimulationTime,
    //     mut mixer_shower: MixerShower,
    // ) {
    //     let flow_rates = vec![5., 7., 9., 11., 13.];
    //     let efficiencies = vec![44.8, 39.1, 34.8, 31.4, 28.6];
    //     let wwhrs = Arc::new(Mutex::new(
    //         WwhrsInstantaneous::new(
    //             flow_rates,
    //             efficiencies,
    //             mixer_shower.cold_water_source.clone(),
    //             0.7,
    //             simulation_time.iter().current_iteration(),
    //         )));
    //     mixer_shower.wwhrs = Some(wwhrs);
    //     let expected_volumes = [24.7, 24.54081632653061, 24.375];
    //
    //     for (idx, t_it) in simulation_time.iter().enumerate() {
    //         assert_eq!(
    //             mixer_shower.hot_water_demand(40., 52., 5., t_it).0,
    //             expected_volumes[idx],
    //             "incorrect volume of hot water returned"
    //         );
    //     }
    // }

    // #[rstest]
    // fn test_wwhrs_instantaneous_system_a_for_mixer(
    //     simulation_time: SimulationTime,
    //     mut mixer_shower: MixerShower,
    // ) {
    //     let flow_rates = vec![5., 7., 9., 11., 13.];
    //     let efficiencies = vec![44.8, 39.1, 34.8, 31.4, 28.6];
    //     let wwhrs = Arc::new(Mutex::new(
    //         WwhrsInstantaneous::new(
    //             flow_rates,
    //             efficiencies,
    //             mixer_shower.cold_water_source.clone(),
    //             0.7,
    //             simulation_time.iter().current_iteration(),
    //         )));
    //     mixer_shower.wwhrs = Some(wwhrs);
    //     let expected_volumes = [22.903242227702766, 22.731048233573134, 22.552562007290966];
    //     let event = TypedScheduleEvent {
    //         start: 0.,
    //         duration: Some(5.),
    //         temperature: 40.,
    //         name: "mixer".to_string(),
    //         event_type: WaterScheduleEventType::Shower,
    //         volume: None,
    //         warm_volume: None,
    //         pipework_volume: None,
    //     };
    //     for (idx, t_it) in simulation_time.iter().enumerate() {
    //         assert_eq!(
    //             mixer_shower.hot_water_demand(event, None, t_it).0,
    //             expected_volumes[idx],
    //             "incorrect volume of hot water returned"
    //         );
    //     }
    // }

    mod test_instant_elec_shower {
        use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyBuilder};
        use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
        use crate::core::water_heat_demand::shower::InstantElectricShower;
        use crate::hem_core::simulation_time::SimulationTime;
        use crate::input::FuelType;
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
            let expected_demands = [86.04206500956023, 175.59605103991885, 268.8814531548757];
            for (idx, t_it) in simulation_time.iter().enumerate() {
                let _ = instant_shower.hot_water_demand(40.0, 52.0, ((idx + 1) * 6) as f64, t_it);
                pretty_assertions::assert_eq!(
                    energy_supply.read().results_by_end_user()["shower"][idx],
                    expected_results_by_end_user[idx],
                    "correct electricity demand not returned"
                );
                pretty_assertions::assert_eq!(
                    instant_shower
                        .hot_water_demand(40.0, 52.0, ((idx + 1) * 6) as f64, t_it)
                        .unwrap()
                        .0,
                    expected_demands[idx]
                );
            }
        }
    }
}
