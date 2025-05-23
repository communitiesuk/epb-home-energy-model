use super::emitters::{Emitters, EmittersDetailedResult};
use super::heat_pump::{BufferTankEmittersData, BufferTankEmittersDataWithResult};
use crate::core::heating_systems::boiler::{
    BoilerServiceSpace, BoilerServiceWaterCombi, BoilerServiceWaterRegular,
};
use crate::core::heating_systems::elec_storage_heater::{
    ElecStorageHeater, StorageHeaterDetailedResult,
};
use crate::core::heating_systems::heat_battery::{
    HeatBatteryServiceSpace, HeatBatteryServiceWaterRegular,
};
use crate::core::heating_systems::heat_network::{
    HeatNetworkServiceSpace, HeatNetworkServiceWaterStorage,
};
use crate::core::heating_systems::heat_pump::{
    HeatPumpHotWaterOnly, HeatPumpServiceSpace, HeatPumpServiceSpaceWarmAir, HeatPumpServiceWater,
};
use crate::core::heating_systems::instant_elec_heater::InstantElecHeater;
use crate::simulation_time::SimulationTimeIteration;
use anyhow::{bail, Error};

#[derive(Debug)]
pub(crate) enum HeatSourceWet {
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
        temp_flow: Option<f64>,
        temperature: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        match self {
            HeatSourceWet::WaterCombi(combi) => Ok(combi.energy_output_max()),
            HeatSourceWet::WaterRegular(regular) => {
                // passing default for _temp_flow and _temp_return as they are unused
                Ok(
                    regular.energy_output_max(
                        Default::default(),
                        Default::default(),
                        None,
                        simtime,
                    ),
                )
            }
            HeatSourceWet::Space(space) => {
                // passing default for _temp_output and _temp_return_feed as they are unused
                Ok(space.energy_output_max(
                    Default::default(),
                    Default::default(),
                    None,
                    None,
                    simtime,
                ))
            }
            HeatSourceWet::HeatNetworkWaterStorage(storage) => {
                // passing default for _temp_flow and _temp_return as they are unused
                Ok(storage.energy_output_max(Default::default(), Default::default(), &simtime))
            }
            HeatSourceWet::HeatBatteryHotWater(battery) => {
                // passing default for _temp_return as it is unused
                battery.energy_output_max(
                    temp_flow.expect(
                        "HeatBatteryHotWater requires a temp_flow when calling energy_output_max",
                    ),
                    Default::default(),
                    simtime,
                )
            }
            HeatSourceWet::HeatPumpWater(hp_water) => hp_water
                .energy_output_max(
                    temp_flow.expect(
                        "HeatPumpWater requires a temp_flow when calling energy_output_max",
                    ),
                    temperature,
                    simtime,
                )
                .map(|x| x.0),
            HeatSourceWet::HeatPumpWaterOnly(hp_water_only) => {
                Ok(hp_water_only.energy_output_max(temperature, simtime))
            }
        }
    }

    pub(crate) fn demand_energy(
        &self,
        energy_demand: f64,
        temp_flow: Option<f64>,
        temperature: f64,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        match self {
            HeatSourceWet::WaterCombi(_combi) => {
                bail!("BoilerServiceWaterCombi does not implement demand_energy")
            }
            HeatSourceWet::WaterRegular(regular) => regular
                .demand_energy(
                    energy_demand,
                    Default::default(),
                    Some(temperature),
                    None,
                    None,
                    None,
                    simtime,
                )
                .map(|x| x.0),
            HeatSourceWet::Space(space) => space
                .demand_energy(
                    energy_demand,
                    Default::default(),
                    Some(temperature),
                    None,
                    None,
                    None,
                    None,
                    simtime,
                )
                .map(|x| x.0),
            HeatSourceWet::HeatNetworkWaterStorage(water_storage) => Ok(water_storage
                .demand_energy(energy_demand, Default::default(), temperature, &simtime)),
            HeatSourceWet::HeatBatteryHotWater(battery) => {
                battery.demand_energy(energy_demand, temperature, None, simtime)
            }
            HeatSourceWet::HeatPumpWater(water) => {
                water.demand_energy(energy_demand, temp_flow, temperature, simtime)
            }
            HeatSourceWet::HeatPumpWaterOnly(water_only) => {
                Ok(water_only.demand_energy(energy_demand, temperature, simtime))
            }
        }
    }

    pub(crate) fn setpnt(
        &self,
        simtime: SimulationTimeIteration,
    ) -> anyhow::Result<(Option<f64>, Option<f64>)> {
        Ok(match self {
            HeatSourceWet::WaterCombi(_) => {
                bail!("BoilerServiceWaterCombi does not implement setpnt")
            }
            HeatSourceWet::WaterRegular(regular) => regular.setpnt(simtime),
            HeatSourceWet::Space(space) => {
                // NOTE - this is temp_setpnt in Python, rather than setpnt
                // investigation needed to see if this is an issue
                let minmax = space.temp_setpnt(simtime);
                (minmax, minmax)
            }
            HeatSourceWet::HeatNetworkWaterStorage(storage) => storage.setpnt(simtime),
            HeatSourceWet::HeatBatteryHotWater(battery) => battery.setpnt(simtime)?,
            HeatSourceWet::HeatPumpWater(heat_pump_water) => heat_pump_water.setpnt(simtime),
            HeatSourceWet::HeatPumpWaterOnly(heat_pump_water_only) => {
                heat_pump_water_only.setpnt(simtime)
            }
        })
    }
}

#[derive(Debug)]
pub(crate) enum SpaceHeatSystem {
    ElecStorage(ElecStorageHeater),
    Instant(InstantElecHeater),
    WarmAir(HeatPumpServiceSpaceWarmAir),
    WetDistribution(Emitters),
}

impl SpaceHeatSystem {
    pub fn temp_setpnt(&self, simulation_time_iteration: SimulationTimeIteration) -> Option<f64> {
        match self {
            SpaceHeatSystem::ElecStorage(elec_storage) => {
                elec_storage.temp_setpnt(&simulation_time_iteration)
            }
            SpaceHeatSystem::Instant(instant) => instant.temp_setpnt(&simulation_time_iteration),
            SpaceHeatSystem::WarmAir(warm_air) => warm_air.temp_setpnt(&simulation_time_iteration),
            SpaceHeatSystem::WetDistribution(wet_distribution) => {
                wet_distribution.temp_setpnt(&simulation_time_iteration)
            }
        }
    }

    pub fn frac_convective(&self, simtime: SimulationTimeIteration) -> f64 {
        match self {
            SpaceHeatSystem::ElecStorage(elec_storage) => elec_storage.frac_convective(),
            SpaceHeatSystem::Instant(instant) => instant.frac_convective(),
            SpaceHeatSystem::WarmAir(warm_air) => warm_air.frac_convective(),
            SpaceHeatSystem::WetDistribution(wet_distribution) => {
                wet_distribution.frac_convective(simtime)
            }
        }
    }

    pub fn _running_time_throughput_factor(
        &self,
        energy_demand: f64,
        space_heat_running_time_cumulative: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        match self {
            SpaceHeatSystem::ElecStorage(..) => unreachable!(), // it isn't expected that this will be called on electric storage heaters
            SpaceHeatSystem::Instant(_instant) => unreachable!(), // it isn't expected that this will be called on instant heaters
            SpaceHeatSystem::WarmAir(warm_air) => warm_air.running_time_throughput_factor(
                energy_demand,
                space_heat_running_time_cumulative,
                simulation_time_iteration,
            ),
            SpaceHeatSystem::WetDistribution(_wet_distribution) => unreachable!(),
        }
    }

    pub fn demand_energy(
        &mut self,
        energy_demand: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        Ok(match self {
            SpaceHeatSystem::ElecStorage(elec_storage) => {
                elec_storage.demand_energy(energy_demand, &simulation_time_iteration)?
            }
            SpaceHeatSystem::Instant(ref mut instant) => {
                instant.demand_energy(energy_demand, simulation_time_iteration)
            }
            SpaceHeatSystem::WarmAir(ref mut warm_air) => {
                warm_air.demand_energy(energy_demand, simulation_time_iteration)?
            }
            SpaceHeatSystem::WetDistribution(ref mut wet_distribution) => {
                wet_distribution.demand_energy(energy_demand, simulation_time_iteration)?
            }
        })
    }

    pub fn in_required_period(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> Option<bool> {
        match self {
            SpaceHeatSystem::ElecStorage(elec_storage) => {
                elec_storage.in_required_period(&simulation_time_iteration)
            }
            SpaceHeatSystem::Instant(instant) => {
                instant.in_required_period(&simulation_time_iteration)
            }
            SpaceHeatSystem::WarmAir(warm_air) => {
                warm_air.in_required_period(&simulation_time_iteration)
            }
            SpaceHeatSystem::WetDistribution(emitters) => {
                emitters.in_required_period(&simulation_time_iteration)
            }
        }
    }

    pub(crate) fn output_emitter_results(&self) -> Option<Vec<EmittersDetailedResult>> {
        if let SpaceHeatSystem::WetDistribution(emitters) = self {
            emitters.output_emitter_results()
        } else {
            None
        }
    }

    pub(crate) fn output_esh_results(&self) -> Option<Vec<StorageHeaterDetailedResult>> {
        if let SpaceHeatSystem::ElecStorage(elec_storage) = self {
            elec_storage.output_esh_results()
        } else {
            None
        }
    }

    pub fn energy_output_min(
        &self,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        match self {
            SpaceHeatSystem::ElecStorage(elec_storage) => {
                elec_storage.energy_output_min(&simulation_time_iteration)
            }
            SpaceHeatSystem::Instant(instant) => Ok(instant.energy_output_min()),
            SpaceHeatSystem::WarmAir(warm_air) => Ok(warm_air.energy_output_min()),
            SpaceHeatSystem::WetDistribution(wet_distribution) => {
                Ok(wet_distribution.energy_output_min(simulation_time_iteration))
            }
        }
    }
}

#[derive(Debug)]
pub enum SpaceHeatingService {
    HeatPump(HeatPumpServiceSpace),
    Boiler(BoilerServiceSpace),
    HeatNetwork(HeatNetworkServiceSpace),
    HeatBattery(HeatBatteryServiceSpace),
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
                    None,
                    emitters_data_for_buffer_tank,
                    simulation_time_iteration,
                ),
            SpaceHeatingService::Boiler(boiler_service_space) => Ok((
                boiler_service_space.energy_output_max(
                    temp_output,
                    temp_return_feed,
                    None,
                    None,
                    simulation_time_iteration,
                ),
                None,
            )),
            SpaceHeatingService::HeatNetwork(heat_network_service_space) => Ok((
                heat_network_service_space.energy_output_max(
                    temp_output,
                    temp_return_feed,
                    None,
                    &simulation_time_iteration,
                ),
                None,
            )),
            SpaceHeatingService::HeatBattery(heat_battery_service_space) => Ok((
                heat_battery_service_space.energy_output_max(
                    temp_output,
                    temp_return_feed,
                    None,
                    simulation_time_iteration,
                )?,
                None,
            )),
        }
    }

    pub(crate) fn demand_energy(
        &mut self,
        energy_demand: f64,
        temp_flow: f64,
        temp_return: f64,
        time_start: Option<f64>,
        emitters_data_for_buffer_tank: Option<BufferTankEmittersDataWithResult>,
        update_heat_source_state: Option<bool>,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> Result<(f64, Option<f64>), Error> {
        match self {
            SpaceHeatingService::HeatPump(heat_pump_service_space) => Ok((
                heat_pump_service_space.demand_energy(
                    energy_demand,
                    temp_flow,
                    temp_return,
                    time_start,
                    emitters_data_for_buffer_tank,
                    update_heat_source_state,
                    simulation_time_iteration,
                )?,
                None,
            )),
            SpaceHeatingService::Boiler(ref mut boiler_service_space) => boiler_service_space
                .demand_energy(
                    energy_demand,
                    temp_flow,
                    Some(temp_return),
                    time_start,
                    None,
                    None,
                    update_heat_source_state,
                    simulation_time_iteration,
                ),
            SpaceHeatingService::HeatNetwork(ref mut heat_network_service_space) => Ok((
                heat_network_service_space.demand_energy(
                    energy_demand,
                    temp_flow,
                    temp_return,
                    time_start,
                    update_heat_source_state,
                    &simulation_time_iteration,
                ),
                None,
            )),
            SpaceHeatingService::HeatBattery(ref heat_battery_service_space) => Ok((
                heat_battery_service_space.demand_energy(
                    energy_demand,
                    temp_flow,
                    temp_return,
                    time_start,
                    update_heat_source_state,
                    simulation_time_iteration,
                )?,
                None,
            )),
        }
    }

    pub(crate) fn _running_time_throughput_factor(
        &self,
        space_heat_running_time_cumulative: f64,
        energy_demand: f64,
        temp_flow: f64,
        temp_return: f64,
        simulation_time_iteration: SimulationTimeIteration,
    ) -> anyhow::Result<(f64, f64)> {
        match self {
            SpaceHeatingService::HeatPump(heat_pump_service_space) => heat_pump_service_space
                .running_time_throughput_factor(
                    space_heat_running_time_cumulative,
                    energy_demand,
                    temp_flow,
                    temp_return,
                    None,
                    simulation_time_iteration,
                ),
            _ => unreachable!(),
        }
    }
}
