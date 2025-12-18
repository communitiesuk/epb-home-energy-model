use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::heating_systems::wwhrs::Wwhrs;
use crate::core::material_properties::WATER;
use crate::core::units::MINUTES_PER_HOUR;
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::misc::frac_hot_water;
use crate::simulation_time::SimulationTimeIteration;
use parking_lot::Mutex;
use std::ops::DerefMut;
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
    ) -> (f64, f64) {
        match self {
            Shower::MixerShower(s) => {
                s.hot_water_demand(temp_target, temp_hot_water, total_shower_duration, simtime)
            }
            Shower::InstantElectricShower(s) => {
                s.hot_water_demand(temp_target, temp_hot_water, total_shower_duration, simtime)
            }
        }
    }
}

#[derive(Debug)]
pub struct MixerShower {
    flowrate: f64,
    cold_water_source: Arc<ColdWaterSource>,
    wwhrs: Option<Arc<Mutex<Wwhrs>>>,
}

impl MixerShower {
    pub fn new(
        flowrate: f64,
        cold_water_source: Arc<ColdWaterSource>,
        wwhrs: Option<Arc<Mutex<Wwhrs>>>,
    ) -> Self {
        Self {
            flowrate,
            cold_water_source,
            wwhrs,
        }
    }

    pub fn get_cold_water_source(&self) -> &ColdWaterSource {
        &self.cold_water_source
    }

    /// Calculate volume of hot water required
    /// (and volume of warm water draining to WWHRS, if applicable)
    ///
    /// Arguments:
    /// * `temp_target` - temperature of warm water delivered at shower head, in Celcius
    /// * `temp_hot_water`
    /// * `total_shower_duration` - cumulative running time of this shower during the current
    ///                             timestep, in minutes
    /// * `timestep_idx` - the index of the timestep for which we are querying the hot water demand
    pub fn hot_water_demand(
        &self,
        temp_target: f64,
        temp_hot_water: f64,
        total_shower_duration: f64,
        simtime: SimulationTimeIteration,
    ) -> (f64, f64) {
        let temp_cold = self.cold_water_source.temperature(simtime);

        let vol_warm_water = self.flowrate * total_shower_duration;
        let mut vol_hot_water =
            vol_warm_water * frac_hot_water(temp_target, temp_hot_water, temp_cold);

        if let Some(wwhrs) = &self.wwhrs {
            let mut wwhrs = wwhrs.lock();

            // Assumed temperature entering WWHRS
            let temp_drain = 35.0;

            let wwhrs_return_temperature =
                wwhrs.return_temperature(temp_drain, self.flowrate, simtime);
            match wwhrs.deref_mut() {
                Wwhrs::WWHRSInstantaneousSystemB(_) => {
                    vol_hot_water = vol_warm_water
                        * frac_hot_water(temp_target, temp_hot_water, wwhrs_return_temperature);
                }
                Wwhrs::WWHRSInstantaneousSystemC(ref mut system_c) => {
                    system_c.set_temperature_for_return(wwhrs_return_temperature)
                }
                Wwhrs::WWHRSInstantaneousSystemA(ref mut system_a) => {
                    system_a.set_temperature_for_return(wwhrs_return_temperature);

                    vol_hot_water = vol_warm_water
                        * frac_hot_water(temp_target, temp_hot_water, wwhrs_return_temperature);
                }
            }
        }

        (vol_hot_water, vol_warm_water)
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
    ) -> (f64, f64) {
        let temp_cold = self.cold_water_source.temperature(simtime);

        let elec_demand =
            self.power_in_kilowatts * (total_shower_duration / MINUTES_PER_HOUR as f64);
        let vol_warm_water =
            elec_demand / WATER.volumetric_energy_content_kwh_per_litre(temp_target, temp_cold);

        let vol_hot_water_equiv =
            vol_warm_water * frac_hot_water(temp_target, temp_hot_water, temp_cold);

        self.energy_supply_connection
            .demand_energy(elec_demand, simtime.index)
            .unwrap();

        (vol_hot_water_equiv, vol_warm_water)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::heating_systems::wwhrs::WWHRSInstantaneousSystemA;
    use crate::core::heating_systems::wwhrs::WWHRSInstantaneousSystemB;
    use crate::core::heating_systems::wwhrs::WWHRSInstantaneousSystemC;
    use crate::simulation_time::SimulationTime;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use std::sync::Arc;

    #[fixture]
    fn simulation_time() -> SimulationTime {
        SimulationTime::new(0f64, 3f64, 1f64)
    }

    #[fixture]
    fn mixer_shower() -> MixerShower {
        let cold_water_temps = [2.0, 3.0, 4.0];
        let cold_water_source = ColdWaterSource::new(cold_water_temps.into(), 0, 1.0);
        MixerShower::new(6.5, cold_water_source.into(), None)
    }

    #[rstest]
    fn test_hot_water_demand(simulation_time: SimulationTime, mixer_shower: MixerShower) {
        let expected_demands = [24.7, 24.54081632653061, 24.375];
        for (idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                mixer_shower.hot_water_demand(40.0, 52.0, 5.0, t_it).0,
                expected_demands[idx],
                "incorrect volume of hot water returned"
            );
        }
    }

    #[rstest]
    fn test_vol_warm_water_for_mixer(simulation_time: SimulationTime, mixer_shower: MixerShower) {
        let expected_volumes = [32.5, 32.5, 32.5];

        for (idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                mixer_shower.hot_water_demand(40., 52., 5., t_it).1,
                expected_volumes[idx],
                "incorrect volume of warm water returned"
            );
        }
    }

    #[rstest]
    fn test_wwhrs_instantaneous_system_b_for_mixer(
        simulation_time: SimulationTime,
        mut mixer_shower: MixerShower,
    ) {
        let flow_rates = vec![5., 7., 9., 11., 13.];
        let efficiencies = vec![44.8, 39.1, 34.8, 31.4, 28.6];
        let wwhrs = Arc::new(Mutex::new(Wwhrs::WWHRSInstantaneousSystemB(
            WWHRSInstantaneousSystemB::new(
                mixer_shower.cold_water_source.clone(),
                flow_rates,
                efficiencies,
                0.7,
            ),
        )));
        mixer_shower.wwhrs = Some(wwhrs);
        let expected_volumes = [22.903242227702766, 22.731048233573134, 22.552562007290966];

        for (idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                mixer_shower.hot_water_demand(40., 52., 5., t_it).0,
                expected_volumes[idx],
                "incorrect volume of hot water returned"
            );
        }
    }

    #[rstest]
    fn test_wwhrs_instantaneous_system_c_for_mixer(
        simulation_time: SimulationTime,
        mut mixer_shower: MixerShower,
    ) {
        let flow_rates = vec![5., 7., 9., 11., 13.];
        let efficiencies = vec![44.8, 39.1, 34.8, 31.4, 28.6];
        let wwhrs = Arc::new(Mutex::new(Wwhrs::WWHRSInstantaneousSystemC(
            WWHRSInstantaneousSystemC::new(
                flow_rates,
                efficiencies,
                mixer_shower.cold_water_source.clone(),
                0.7,
                simulation_time.iter().current_iteration(),
            ),
        )));
        mixer_shower.wwhrs = Some(wwhrs);
        let expected_volumes = [24.7, 24.54081632653061, 24.375];

        for (idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                mixer_shower.hot_water_demand(40., 52., 5., t_it).0,
                expected_volumes[idx],
                "incorrect volume of hot water returned"
            );
        }
    }

    #[rstest]
    fn test_wwhrs_instantaneous_system_a_for_mixer(
        simulation_time: SimulationTime,
        mut mixer_shower: MixerShower,
    ) {
        let flow_rates = vec![5., 7., 9., 11., 13.];
        let efficiencies = vec![44.8, 39.1, 34.8, 31.4, 28.6];
        let wwhrs = Arc::new(Mutex::new(Wwhrs::WWHRSInstantaneousSystemA(
            WWHRSInstantaneousSystemA::new(
                flow_rates,
                efficiencies,
                mixer_shower.cold_water_source.clone(),
                0.7,
                simulation_time.iter().current_iteration(),
            ),
        )));
        mixer_shower.wwhrs = Some(wwhrs);
        let expected_volumes = [22.903242227702766, 22.731048233573134, 22.552562007290966];

        for (idx, t_it) in simulation_time.iter().enumerate() {
            assert_eq!(
                mixer_shower.hot_water_demand(40., 52., 5., t_it).0,
                expected_volumes[idx],
                "incorrect volume of hot water returned"
            );
        }
    }

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
                instant_shower.hot_water_demand(40.0, 52.0, ((idx + 1) * 6) as f64, t_it);
                pretty_assertions::assert_eq!(
                    energy_supply.read().results_by_end_user()["shower"][idx],
                    expected_results_by_end_user[idx],
                    "correct electricity demand not returned"
                );
                pretty_assertions::assert_eq!(
                    instant_shower
                        .hot_water_demand(40.0, 52.0, ((idx + 1) * 6) as f64, t_it)
                        .0,
                    expected_demands[idx]
                );
            }
        }
    }
}
