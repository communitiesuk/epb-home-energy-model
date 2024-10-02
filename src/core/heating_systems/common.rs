use crate::core::heating_systems::boiler::{
    BoilerServiceSpace, BoilerServiceWaterCombi, BoilerServiceWaterRegular,
};
use crate::core::heating_systems::heat_battery::HeatBatteryServiceWaterRegular;
use crate::core::heating_systems::heat_network::{
    HeatNetworkServiceSpace, HeatNetworkServiceWaterStorage,
};
use crate::core::heating_systems::heat_pump::{
    HeatPumpHotWaterOnly, HeatPumpServiceSpace, HeatPumpServiceSpaceWarmAir, HeatPumpServiceWater,
};
use crate::core::heating_systems::instant_elec_heater::InstantElecHeater;
use crate::simulation_time::SimulationTimeIteration;

#[derive(Clone, Debug)]
pub enum HeatSourceWet {
    WaterCombi(BoilerServiceWaterCombi),
    WaterRegular(BoilerServiceWaterRegular),
    Space(BoilerServiceSpace),
    HeatNetworkWaterStorage(HeatNetworkServiceWaterStorage),
    HeatBatteryHotWater(HeatBatteryServiceWaterRegular),
    HeatPumpWater(HeatPumpServiceWater),
    HeatPumpWaterOnly(HeatPumpHotWaterOnly),
}

impl HeatSourceWet {
    /// Common way of calling energy_output_max() on heat sources, implementing equivalent of duck-typing happening in upstream Python.
    pub(crate) fn energy_output_max(
        &self,
        temperature: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        match self {
            HeatSourceWet::WaterCombi(combi) => Ok(combi.energy_output_max()),
            HeatSourceWet::WaterRegular(regular) => {
                Ok(regular.energy_output_max(temperature, None, simtime))
            }
            HeatSourceWet::Space(space) => {
                Ok(space.energy_output_max(temperature, Default::default(), None, simtime))
            }
            HeatSourceWet::HeatNetworkWaterStorage(storage) => {
                Ok(storage.energy_output_max(temperature, &simtime))
            }
            HeatSourceWet::HeatBatteryHotWater(battery) => {
                Ok(battery.energy_output_max(temperature, simtime))
            }
            HeatSourceWet::HeatPumpWater(hp_water) => hp_water
                .energy_output_max(temperature, simtime)
                .map(|x| x.0),
            HeatSourceWet::HeatPumpWaterOnly(hp_water_only) => {
                Ok(hp_water_only.energy_output_max(temperature, simtime))
            }
        }
    }

    pub(crate) fn demand_energy(
        &mut self,
        energy_demand: f64,
        temperature: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        match self {
            HeatSourceWet::WaterCombi(_combi) => {
                unimplemented!("BoilerServiceWaterCombi does not implement demand_energy")
            }
            HeatSourceWet::WaterRegular(regular) => regular
                .demand_energy(energy_demand, temperature, None, None, simtime)
                .map(|x| x.0),
            HeatSourceWet::Space(space) => space
                .demand_energy(
                    energy_demand,
                    Default::default(),
                    temperature,
                    None,
                    None,
                    simtime,
                )
                .map(|x| x.0),
            HeatSourceWet::HeatNetworkWaterStorage(water_storage) => {
                Ok(water_storage.demand_energy(energy_demand, temperature, &simtime))
            }
            HeatSourceWet::HeatBatteryHotWater(battery) => {
                Ok(battery.demand_energy(energy_demand, temperature, simtime))
            }
            HeatSourceWet::HeatPumpWater(water) => {
                water.demand_energy(energy_demand, temperature, simtime)
            }
            HeatSourceWet::HeatPumpWaterOnly(water_only) => {
                Ok(water_only.demand_energy(energy_demand, temperature, simtime))
            }
        }
    }
}

#[derive(Clone)]
pub enum SpaceHeatSystem {
    Instant(InstantElecHeater),
    WarmAir(HeatPumpServiceSpaceWarmAir),
}

impl SpaceHeatSystem {
    pub fn temp_setpnt(&self, simulation_time_iteration: SimulationTimeIteration) -> Option<f64> {
        match self {
            SpaceHeatSystem::Instant(instant) => instant.temp_setpnt(&simulation_time_iteration),
            SpaceHeatSystem::WarmAir(warm_air) => warm_air.temp_setpnt(&simulation_time_iteration),
        }
    }

    pub fn frac_convective(&self) -> f64 {
        match self {
            SpaceHeatSystem::Instant(instant) => instant.frac_convective(),
            SpaceHeatSystem::WarmAir(warm_air) => warm_air.frac_convective(),
        }
    }

    pub fn running_time_throughput_factor(
        &self,
        energy_demand: f64,
        space_heat_running_time_cumulative: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        match self {
            SpaceHeatSystem::Instant(_instant) => unreachable!(), // it isn't expected that this will be called on instant heaters
            SpaceHeatSystem::WarmAir(warm_air) => warm_air.running_time_throughput_factor(
                energy_demand,
                space_heat_running_time_cumulative,
                simulation_time_iteration,
            ),
        }
    }

    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        Ok(match self {
            SpaceHeatSystem::Instant(ref mut instant) => {
                instant.demand_energy(energy_demand, simulation_time_iteration)
            }
            SpaceHeatSystem::WarmAir(ref mut warm_air) => {
                warm_air.demand_energy(energy_demand, simulation_time_iteration)?
            }
        })
    }

    pub fn in_required_period(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> Option<bool> {
        match self {
            SpaceHeatSystem::Instant(instant) => {
                instant.in_required_period(&simulation_time_iteration)
            }
            SpaceHeatSystem::WarmAir(warm_air) => {
                warm_air.in_required_period(&simulation_time_iteration)
            }
        }
    }
}

pub enum SpaceHeatingService {
    HeatPump(HeatPumpServiceSpace),
    Boiler(BoilerServiceSpace),
    HeatNetwork(HeatNetworkServiceSpace),
    HeatBattery(()),
}

impl SpaceHeatingService {
    pub(crate) fn energy_output_max(
        &self,
        temp_output: f64,
        temp_return_feed: f64,
        emitters_data_for_buffer_tank: Option<BufferTankEmittersData>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> Result<(f64, Option<BufferTankEmittersDataWithResult>), Error> {
        match self {
            SpaceHeatingService::HeatPump(heat_pump_service_space) => heat_pump_service_space
                .energy_output_max(
                    temp_output,
                    temp_return_feed,
                    emitters_data_for_buffer_tank,
                    simulation_time_iteration,
                ),
            SpaceHeatingService::Boiler(boiler_service_space) => {
                let time_elapsed_hp: Option<f64> = None; // TODO is this assumption correct?
                Ok((
                    boiler_service_space.energy_output_max(
                        temp_output,
                        temp_return_feed,
                        time_elapsed_hp,
                        simulation_time_iteration,
                    ),
                    None,
                ))
            }
            SpaceHeatingService::HeatNetwork(heat_network_service_space) => Ok((
                heat_network_service_space.energy_output_max(
                    temp_output,
                    temp_return_feed,
                    &simulation_time_iteration,
                ),
                None,
            )),
            SpaceHeatingService::HeatBattery(_) => unimplemented!(),
        }
    }
}

// macro so accessing individual controls through the enum isn't so repetitive
#[macro_use]
macro_rules! per_space_heating {
    ($val:expr, $pattern:pat => { $res:expr }) => {
        match $val {
            SpaceHeatingService::HeatPump($pattern) => $res,
            SpaceHeatingService::Boiler($pattern) => $res,
            SpaceHeatingService::HeatNetwork($pattern) => $res,
            SpaceHeatingService::HeatBattery($pattern) => unreachable!(),
        }
    };
}

use anyhow::Error;

use super::heat_pump::{BufferTankEmittersData, BufferTankEmittersDataWithResult};
