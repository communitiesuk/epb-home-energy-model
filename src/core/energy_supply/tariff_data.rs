use crate::input::EnergySupplyTariff;
use anyhow::anyhow;
use serde::Deserialize;
use std::io::Read;

/// This module contains data on the energy tariffs.

#[derive(Clone, Debug)]
pub(crate) struct TariffData {
    elec_prices: Vec<TariffRow>,
}

#[derive(Clone, Debug, Deserialize)]
struct TariffRow {
    #[serde(rename = "Standard Tariff")]
    standard: f64,
    #[serde(rename = "7-Hour Off Peak Tariff")]
    seven_hour_off_peak: f64,
    #[serde(rename = "10-Hour Off Peak Tariff")]
    ten_hour_off_peak: f64,
    #[serde(rename = "Variable Time of Day Tariff")]
    variable_time_of_day: f64,
}

impl TariffRow {
    fn get(&self, tariff: &EnergySupplyTariff) -> f64 {
        match tariff {
            EnergySupplyTariff::Standard => self.standard,
            EnergySupplyTariff::SevenHourOffPeak => self.seven_hour_off_peak,
            EnergySupplyTariff::TenHourOffPeak => self.ten_hour_off_peak,
            EnergySupplyTariff::VariableTimeOfDay => self.variable_time_of_day,
        }
    }
}

impl TariffData {
    pub(crate) fn new(csv: impl Read) -> anyhow::Result<Self> {
        // we're making the assumption here that the timestep column will always start with zero and
        // increment by one on each row, so we can just infer the indexes for the row records
        Ok(Self {
            elec_prices: csv::Reader::from_reader(csv)
                .deserialize::<TariffRow>()
                .collect::<Result<_, _>>()?,
        })
    }

    pub(crate) fn price(
        &self,
        tariff: &EnergySupplyTariff,
        timestep_id: usize,
    ) -> anyhow::Result<f64> {
        // TODO (from Python) Current solution is based on tariff data file that has at least as many timesteps (whatever length)
        //                    as the simtime object in the json file. This will need to be revisited when tariffs move to
        //                    PCDB and its processing is integrated alongside other database objects.
        Ok(self
            .elec_prices
            .get(timestep_id)
            .ok_or_else(|| {
                anyhow!("There were no tariff electricity prices for the timestep ID {timestep_id}")
            })?
            .get(tariff))
    }
}
