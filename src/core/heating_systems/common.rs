use crate::core::heating_systems::boiler::{
    BoilerServiceSpace, BoilerServiceWaterCombi, BoilerServiceWaterRegular,
};
use crate::core::heating_systems::heat_battery::HeatBatteryServiceWaterRegular;
use crate::core::heating_systems::heat_network::HeatNetworkServiceWaterStorage;
use crate::core::heating_systems::heat_pump::{HeatPumpHotWaterOnly, HeatPumpServiceSpaceWarmAir, HeatPumpServiceWater};
use crate::core::heating_systems::instant_elec_heater::InstantElecHeater;

#[derive(Clone)]
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
