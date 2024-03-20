use crate::core::heating_systems::boiler::{
    BoilerServiceSpace, BoilerServiceWaterCombi, BoilerServiceWaterRegular,
};
use crate::core::heating_systems::heat_battery::HeatBatteryServiceWaterRegular;
use crate::core::heating_systems::heat_network::HeatNetworkServiceWaterStorage;
use crate::core::heating_systems::heat_pump::{
    HeatPumpHotWaterOnly, HeatPumpServiceSpaceWarmAir, HeatPumpServiceWater,
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
    ) -> (f64, f64) {
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
    ) -> f64 {
        match self {
            SpaceHeatSystem::Instant(ref mut instant) => {
                instant.demand_energy(energy_demand, simulation_time_iteration)
            }
            SpaceHeatSystem::WarmAir(ref mut warm_air) => {
                warm_air.demand_energy(energy_demand, simulation_time_iteration)
            }
        }
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
