use anyhow::anyhow;
use serde::Deserialize;
use std::io::Read;

/// This module contains data on the energy tariffs.

#[derive(Clone, Copy, Debug)]
pub(super) enum Tariff {
    Standard,
    SevenHourOffPeak,
    TenHourOffPeak,
    VariableTimeOfDay,
}

#[derive(Clone, Debug)]
pub(super) struct TariffData {
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
    fn get(&self, tariff: &Tariff) -> f64 {
        match tariff {
            Tariff::Standard => self.standard,
            Tariff::SevenHourOffPeak => self.seven_hour_off_peak,
            Tariff::TenHourOffPeak => self.ten_hour_off_peak,
            Tariff::VariableTimeOfDay => self.variable_time_of_day,
        }
    }
}

impl TariffData {
    pub(super) fn new(csv: impl Read) -> anyhow::Result<Self> {
        // we're making the assumption here that the timestep column will always start with zero and
        // increment by one on each row, so we can just infer the indexes for the row records
        Ok(Self {
            elec_prices: csv::Reader::from_reader(csv)
                .deserialize::<TariffRow>()
                .collect::<Result<_, _>>()?,
        })
    }

    pub(super) fn price(&self, tariff: &Tariff, timestep_id: usize) -> anyhow::Result<f64> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;
    use std::io::{BufReader, Cursor};

    #[rstest]
    fn test_parse_fixture_file() {
        let data = BufReader::new(Cursor::new(include_str!(
            "../../../examples/tariff_data/tariff_data_25-06-2024.csv"
        )));
        let tariff_data = TariffData::new(data).unwrap();
        assert_eq!(
            tariff_data
                .price(&Tariff::VariableTimeOfDay, 10644)
                .unwrap(),
            22.92342657
        );
    }
}
