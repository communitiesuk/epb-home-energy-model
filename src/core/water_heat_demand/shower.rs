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
    use crate::core::energy_supply::energy_supply::{EnergySupply, EnergySupplyBuilder};
    use crate::input::FuelType;
    use crate::simulation_time::SimulationTime;
    use parking_lot::RwLock;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use std::sync::Arc;

    #[rstest]
    pub fn should_calculate_correct_hot_water_demand_for_mixer() {
        let simulation_time = SimulationTime::new(0f64, 3f64, 1f64);
        let cold_water_temps = [2.0, 3.0, 4.0];
        let cold_water_source = ColdWaterSource::new(cold_water_temps.into(), 0, 1.0);
        let mixer_shower = MixerShower::new(6.5, cold_water_source.into(), None);
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
    pub fn should_calculate_correct_hot_water_demand_for_instant() {
        let simulation_time = SimulationTime::new(0f64, 3f64, 1f64);
        let cold_water_temps = [2.0, 3.0, 4.0];
        let cold_water_source = ColdWaterSource::new(cold_water_temps.into(), 0, 1.0);
        let energy_supply = Arc::new(RwLock::new(
            EnergySupplyBuilder::new(FuelType::Electricity, simulation_time.total_steps()).build(),
        ));
        let energy_supply_conn = EnergySupply::connection(energy_supply.clone(), "shower").unwrap();
        let instant_shower =
            InstantElectricShower::new(50.0, cold_water_source.into(), energy_supply_conn);
        let expected_results_by_end_user = [5.0, 10.0, 15.0];
        let expected_demands = [86.04206500956023, 175.59605103991885, 268.8814531548757];
        for (idx, t_it) in simulation_time.iter().enumerate() {
            instant_shower.hot_water_demand(40.0, 52.0, ((idx + 1) * 6) as f64, t_it);
            assert_eq!(
                energy_supply.read().results_by_end_user()["shower"][idx],
                expected_results_by_end_user[idx],
                "correct electricity demand not returned"
            );
            assert_eq!(
                instant_shower
                    .hot_water_demand(40.0, 52.0, ((idx + 1) * 6) as f64, t_it)
                    .0,
                expected_demands[idx]
            );
        }
    }
}
