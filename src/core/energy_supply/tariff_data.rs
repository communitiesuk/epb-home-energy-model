use crate::core::schedule::{
    expand_numeric_schedule, reject_nulls, validate_schedule_length, NumericSchedule,
};
use crate::hem_core::simulation_time::{SimulationTimeIteration, SimulationTimeIterator};
use crate::input::EnergySupplyTariff;
use anyhow::anyhow;
use indexmap::IndexMap;
use std::io::Read;

/// This module contains data on the energy tariffs.

#[derive(Clone, Debug)]
pub(super) struct TariffData {
    start_day: Option<u32>,
    time_series_step: f64,
    electricity_prices: IndexMap<EnergySupplyTariff, Vec<f64>>,
}

impl TariffData {
    pub(super) fn new(
        simulation_time: SimulationTimeIterator,
        start_day: Option<u32>,
        time_series_step: f64,
        electricity_prices: IndexMap<EnergySupplyTariff, Vec<f64>>,
    ) -> anyhow::Result<Self> {
        let expected_length = match start_day {
            Some(start_day) => {
                simulation_time.total_steps_based_on_step(start_day, Some(time_series_step))?
            }
            None => simulation_time.total_steps(),
        };

        for prices in electricity_prices.values() {
            validate_schedule_length(prices, expected_length)?;
        }

        Ok(TariffData {
            start_day,
            time_series_step,
            electricity_prices,
        })
    }

    pub(super) fn expand_prices_schedule(
        prices: IndexMap<EnergySupplyTariff, NumericSchedule>,
    ) -> anyhow::Result<IndexMap<EnergySupplyTariff, Vec<f64>>> {
        prices
            .iter()
            .map(|(name, schedule)| {
                let prices = reject_nulls(expand_numeric_schedule(schedule))?;

                Ok((*name, prices))
            })
            .collect()
    }

    pub(super) fn load_data_from_file(
        csv: impl Read,
    ) -> anyhow::Result<IndexMap<EnergySupplyTariff, NumericSchedule>> {
        let mut reader = csv::Reader::from_reader(csv);
        let headers = reader.headers()?.clone();
        if !headers.iter().any(|h| h == "timestep") {
            anyhow::bail!("Tariff CSV must contain a 'timestep' column");
        }
        let mut electricity_prices: IndexMap<EnergySupplyTariff, Vec<f64>> = IndexMap::new();

        for row in reader.records().flatten() {
            for (tariff, price) in headers.iter().zip(row.iter()) {
                if tariff == "timestep" {
                    continue;
                }

                if let Ok(tariff) = tariff.parse::<EnergySupplyTariff>() {
                    let price: f64 = price.parse()?;

                    electricity_prices.entry(tariff).or_default().push(price);
                }
            }
        }

        let electricity_price_schedules = electricity_prices
            .into_iter()
            .map(|(tariff, prices)| (tariff, prices.into()))
            .collect();

        Ok(electricity_price_schedules)
    }

    pub(super) fn price(
        &self,
        tariff: &EnergySupplyTariff,
        simulation_time: SimulationTimeIteration,
    ) -> anyhow::Result<f64> {
        let prices = self
            .electricity_prices
            .get(tariff)
            .ok_or_else(|| anyhow!("Tariff ({tariff}) not found"))?;

        let idx = match self.start_day {
            // no start_day means the tariff starts alongside the simulation
            None => simulation_time.index,
            Some(start_day) => simulation_time.time_series_idx(start_day, self.time_series_step),
        };

        prices
            .get(idx)
            .copied()
            .ok_or_else(|| anyhow!("Index ({idx}) out of bounds for {tariff}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hem_core::simulation_time::SimulationTime;
    use rstest::*;
    use std::io::{BufReader, Cursor};

    #[fixture]
    pub fn simulation_time() -> SimulationTimeIterator {
        SimulationTime::new(0.0, 24.0, 1.0).iter()
    }

    #[fixture]
    fn loaded_prices() -> IndexMap<EnergySupplyTariff, NumericSchedule> {
        let data = BufReader::new(Cursor::new(include_str!(
            "../../../examples/tariff_data/tariff_data_demo_files_24timesteps.csv"
        )));

        TariffData::load_data_from_file(data).unwrap()
    }

    #[fixture]
    fn tariff_data(
        loaded_prices: IndexMap<EnergySupplyTariff, NumericSchedule>,
        simulation_time: SimulationTimeIterator,
    ) -> TariffData {
        let prices = TariffData::expand_prices_schedule(loaded_prices).unwrap();

        TariffData::new(simulation_time, Some(0), 1., prices).unwrap()
    }

    #[fixture]
    fn tariff_data_no_start_day(
        loaded_prices: IndexMap<EnergySupplyTariff, NumericSchedule>,
        simulation_time: SimulationTimeIterator,
    ) -> TariffData {
        let prices = TariffData::expand_prices_schedule(loaded_prices).unwrap();

        TariffData::new(simulation_time, None, 1., prices).unwrap()
    }

    #[rstest]
    fn test_price(tariff_data: TariffData, simulation_time: SimulationTimeIterator) {
        for (t_idx, t_it) in simulation_time.enumerate() {
            if t_idx == 0 {
                assert_eq!(
                    tariff_data
                        .price(&EnergySupplyTariff::Standard, t_it)
                        .unwrap(),
                    25.16
                );
                assert_eq!(
                    tariff_data
                        .price(&EnergySupplyTariff::SevenHourOffPeak, t_it)
                        .unwrap(),
                    14.6
                );
                assert_eq!(
                    tariff_data
                        .price(&EnergySupplyTariff::TenHourOffPeak, t_it)
                        .unwrap(),
                    16.04
                );
                assert_eq!(
                    tariff_data
                        .price(&EnergySupplyTariff::VariableTimeOfDay, t_it)
                        .unwrap(),
                    10.87017271
                );
            }
            if t_idx == 23 {
                assert_eq!(
                    tariff_data
                        .price(&EnergySupplyTariff::Standard, t_it)
                        .unwrap(),
                    25.16
                );
                assert_eq!(
                    tariff_data
                        .price(&EnergySupplyTariff::SevenHourOffPeak, t_it)
                        .unwrap(),
                    29.8
                );
                assert_eq!(
                    tariff_data
                        .price(&EnergySupplyTariff::TenHourOffPeak, t_it)
                        .unwrap(),
                    35.01
                );
                assert_eq!(
                    tariff_data
                        .price(&EnergySupplyTariff::VariableTimeOfDay, t_it)
                        .unwrap(),
                    23.91911304
                );
            }
        }
    }

    #[rstest]
    fn test_price_no_start_day(
        tariff_data_no_start_day: TariffData,
        simulation_time: SimulationTimeIterator,
    ) {
        for (t_idx, t_it) in simulation_time.enumerate() {
            if t_idx == 0 {
                assert_eq!(
                    tariff_data_no_start_day
                        .price(&EnergySupplyTariff::Standard, t_it)
                        .unwrap(),
                    25.16
                );
                assert_eq!(
                    tariff_data_no_start_day
                        .price(&EnergySupplyTariff::SevenHourOffPeak, t_it)
                        .unwrap(),
                    14.6
                );
                assert_eq!(
                    tariff_data_no_start_day
                        .price(&EnergySupplyTariff::TenHourOffPeak, t_it)
                        .unwrap(),
                    16.04
                );
                assert_eq!(
                    tariff_data_no_start_day
                        .price(&EnergySupplyTariff::VariableTimeOfDay, t_it)
                        .unwrap(),
                    10.87017271
                );
            }
            if t_idx == 23 {
                assert_eq!(
                    tariff_data_no_start_day
                        .price(&EnergySupplyTariff::Standard, t_it)
                        .unwrap(),
                    25.16
                );
                assert_eq!(
                    tariff_data_no_start_day
                        .price(&EnergySupplyTariff::SevenHourOffPeak, t_it)
                        .unwrap(),
                    29.8
                );
                assert_eq!(
                    tariff_data_no_start_day
                        .price(&EnergySupplyTariff::TenHourOffPeak, t_it)
                        .unwrap(),
                    35.01
                );
                assert_eq!(
                    tariff_data_no_start_day
                        .price(&EnergySupplyTariff::VariableTimeOfDay, t_it)
                        .unwrap(),
                    23.91911304
                );
            }
        }
    }

    // skipping python's test_get_price_out_of_range as not possible to pass invalid tariff in rust
}
