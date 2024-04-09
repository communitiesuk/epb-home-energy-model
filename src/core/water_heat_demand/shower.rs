use crate::core::energy_supply::energy_supply::EnergySupplyConnection;
use crate::core::heating_systems::wwhrs::Wwhrs;
use crate::core::material_properties::WATER;
use crate::core::units::MINUTES_PER_HOUR;
use crate::core::water_heat_demand::cold_water_source::ColdWaterSource;
use crate::core::water_heat_demand::misc::frac_hot_water;

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
        &mut self,
        temp_target: f64,
        total_shower_duration: f64,
        timestep_idx: usize,
    ) -> f64 {
        match self {
            Shower::MixerShower(s) => {
                s.hot_water_demand(temp_target, total_shower_duration, timestep_idx)
            }
            Shower::InstantElectricShower(s) => {
                s.hot_water_demand(temp_target, total_shower_duration, timestep_idx)
            }
        }
    }
}

pub struct MixerShower {
    flowrate: f64,
    cold_water_source: ColdWaterSource,
    wwhrs: Option<Wwhrs>,
    temp_hot: f64,
}

impl MixerShower {
    pub fn new(flowrate: f64, cold_water_source: ColdWaterSource, wwhrs: Option<Wwhrs>) -> Self {
        Self {
            flowrate,
            cold_water_source,
            wwhrs,
            temp_hot: 52.0, // TODO Python (intent is to not hard-code this)
        }
    }

    pub fn get_cold_water_source(&self) -> &ColdWaterSource {
        &self.cold_water_source
    }

    pub fn get_temp_hot(&self) -> f64 {
        self.temp_hot
    }

    /// Calculate volume of hot water required
    /// (and volume of warm water draining to WWHRS, if applicable)
    ///
    /// Arguments:
    /// * `temp_target` - temperature of warm water delivered at shower head, in Celcius
    /// * `total_shower_duration` - cumulative running time of this shower during the current
    ///                             timestep, in minutes
    /// * `timestep_idx` - the index of the timestep for which we are querying the hot water demand
    pub fn hot_water_demand(
        &mut self,
        temp_target: f64,
        total_shower_duration: f64,
        timestep_idx: usize,
    ) -> f64 {
        let temp_cold = self.cold_water_source.temperature(timestep_idx);

        let vol_warm_water = self.flowrate * total_shower_duration;
        let mut vol_hot_water =
            vol_warm_water * frac_hot_water(temp_target, self.temp_hot, temp_cold);

        if let Some(wwhrs) = &mut self.wwhrs {
            let wwhrs_return_temperature =
                wwhrs.return_temperature(temp_target, self.flowrate, timestep_idx);
            match wwhrs {
                Wwhrs::WWHRSInstantaneousSystemB(_) => {
                    vol_hot_water = vol_warm_water
                        * frac_hot_water(temp_target, self.temp_hot, wwhrs_return_temperature);
                }
                Wwhrs::WWHRSInstantaneousSystemC(system_c) => {
                    system_c.set_temperature_for_return(wwhrs_return_temperature)
                }
                Wwhrs::WWHRSInstantaneousSystemA(system_a) => {
                    system_a.set_temperature_for_return(wwhrs_return_temperature);

                    vol_hot_water = vol_warm_water
                        * frac_hot_water(temp_target, self.temp_hot, wwhrs_return_temperature);
                }
            }
        }

        vol_hot_water
    }
}

pub struct InstantElectricShower {
    power_in_kilowatts: f64,
    cold_water_source: ColdWaterSource,
    energy_supply_connection: EnergySupplyConnection,
}

impl InstantElectricShower {
    pub fn new(
        power_in_kilowatts: f64,
        cold_water_source: ColdWaterSource,
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
        &mut self,
        temp_target: f64,
        total_shower_duration: f64,
        timestep_idx: usize,
    ) -> f64 {
        let temp_cold = self.cold_water_source.temperature(timestep_idx);

        let elec_demand =
            self.power_in_kilowatts * (total_shower_duration / MINUTES_PER_HOUR as f64);
        let vol_warm_water =
            elec_demand / WATER.volumetric_energy_content_kwh_per_litre(temp_target, temp_cold);
        let temp_hot = 52f64; // TODO note in Python to change source of this instead of hard-coding

        let vol_hot_water_equiv = vol_warm_water * frac_hot_water(temp_target, temp_hot, temp_cold);

        self.energy_supply_connection
            .demand_energy(elec_demand, timestep_idx)
            .unwrap();

        vol_hot_water_equiv
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::energy_supply::energy_supply::EnergySupply;
    use crate::input::EnergySupplyType;
    use crate::simulation_time::SimulationTime;
    use parking_lot::Mutex;
    use rstest::*;
    use std::sync::Arc;

    #[rstest]
    pub fn should_calculate_correct_hot_water_demand_for_mixer() {
        let simulation_time = SimulationTime::new(0f64, 3f64, 1f64);
        let cold_water_temps = [2.0, 3.0, 4.0];
        let cold_water_source =
            ColdWaterSource::new(cold_water_temps.into(), &simulation_time, 1.0);
        let mut mixer_shower = MixerShower::new(6.5, cold_water_source, None);
        let expected_demands = [24.7, 24.54081632653061, 24.375];
        for (idx, _) in simulation_time.iter().enumerate() {
            assert_eq!(
                mixer_shower.hot_water_demand(40.0, 5.0, idx),
                expected_demands[idx],
                "incorrect volume of hot water returned"
            );
        }
    }

    pub fn should_calculate_correct_hot_water_demand_for_instant() {
        let simulation_time = SimulationTime::new(0f64, 3f64, 1f64);
        let cold_water_temps = [2.0, 3.0, 4.0];
        let cold_water_source =
            ColdWaterSource::new(cold_water_temps.into(), &simulation_time, 1.0);
        let energy_supply = Arc::new(Mutex::new(EnergySupply::new(
            EnergySupplyType::Electricity,
            simulation_time.total_steps(),
            None,
        )));
        let energy_supply_conn = EnergySupply::connection(energy_supply.clone(), "shower").unwrap();
        let mut instant_shower =
            InstantElectricShower::new(50.0, cold_water_source, energy_supply_conn);
        let expected_results_by_end_user = [5.0, 10.0, 15.0];
        let expected_demands = [86.04206500956023, 175.59605103991885, 268.8814531548757];
        for (idx, _) in simulation_time.iter().enumerate() {
            assert_eq!(
                energy_supply.lock().results_by_end_user()["shower"][idx],
                expected_results_by_end_user[idx],
                "correct electricity demand not returned"
            );
            assert_eq!(
                instant_shower.hot_water_demand(40.0, ((idx + 1) * 6) as f64, idx),
                expected_demands[idx]
            );
        }
    }
}
