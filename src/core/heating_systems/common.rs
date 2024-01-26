use crate::core::heating_systems::boiler::{
    BoilerServiceSpace, BoilerServiceWaterCombi, BoilerServiceWaterRegular,
};
use crate::core::heating_systems::heat_network::HeatNetworkServiceWaterStorage;

#[derive(Clone)]
pub enum HeatSourceWet {
    WaterCombi(BoilerServiceWaterCombi),
    WaterRegular(BoilerServiceWaterRegular),
    Space(BoilerServiceSpace),
    HeatNetworkWaterStorage(HeatNetworkServiceWaterStorage),
}
