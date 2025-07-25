#![allow(dead_code)]

use crate::wrappers::future_homes_standard::future_homes_standard::{calc_final_rates, FinalRates};
use crate::wrappers::future_homes_standard::future_homes_standard_fee::calc_fabric_energy_efficiency;
use crate::{build_summary_data, CalculationKey, CalculationResultsWithContext, SummaryData};
use anyhow::anyhow;
use indexmap::IndexMap;
use itertools::Itertools;
use serde::Serialize;
use smartstring::alias::String;
use std::collections::HashMap;

#[derive(Serialize)]
pub struct FhsComplianceResponse {
    dwelling_emission_rate: f64,
    target_emission_rate: f64,
    emission_rate_compliant: bool,
    dwelling_primary_energy_rate: f64,
    target_primary_energy_rate: f64,
    primary_energy_rate_compliant: bool,
    dwelling_fabric_energy_efficiency: f64,
    target_fabric_energy_efficiency: f64,
    fabric_energy_efficiency_compliant: bool,
    energy_demand: EnergyDemand,
    delivered_energy_use: DeliveredEnergyUse,
    energy_use_by_fuel: IndexMap<String, PerformanceValue>,
}

impl FhsComplianceResponse {
    pub(super) fn build_from<T>(result: &T) -> anyhow::Result<Self>
    where
        T: FhsComplianceCalculationResult,
    {
        let dwelling_emission_rate = result.dwelling_emission_rate();
        let target_emission_rate = result.target_emission_rate();
        // compliance if dwelling rate is less than or equal to target rate
        let emission_rate_compliant = dwelling_emission_rate <= target_emission_rate;
        let dwelling_primary_energy_rate = result.dwelling_primary_energy_rate();
        let target_primary_energy_rate = result.target_primary_energy_rate();
        // compliance if dwelling rate is less than or equal to target rate
        let primary_energy_rate_compliant =
            dwelling_primary_energy_rate <= target_primary_energy_rate;
        let dwelling_fabric_energy_efficiency = result.dwelling_fabric_energy_efficiency();
        let target_fabric_energy_efficiency = result.target_fabric_energy_efficiency();
        let fabric_energy_efficiency_compliant =
            dwelling_fabric_energy_efficiency <= target_fabric_energy_efficiency;

        Ok(Self {
            dwelling_emission_rate,
            target_emission_rate,
            emission_rate_compliant,
            dwelling_primary_energy_rate,
            target_primary_energy_rate,
            primary_energy_rate_compliant,
            dwelling_fabric_energy_efficiency,
            target_fabric_energy_efficiency,
            fabric_energy_efficiency_compliant,
            energy_demand: result.energy_demand().clone(),
            delivered_energy_use: result.delivered_energy_use().clone(),
            energy_use_by_fuel: result.energy_use_by_fuel().clone(),
        })
    }
}

/// A trait to define data on the underlying calculation result that an FHS compliance result requires.
pub(super) trait FhsComplianceCalculationResult {
    fn dwelling_emission_rate(&self) -> f64;
    fn target_emission_rate(&self) -> f64;
    fn dwelling_primary_energy_rate(&self) -> f64;
    fn target_primary_energy_rate(&self) -> f64;
    fn dwelling_fabric_energy_efficiency(&self) -> f64;
    fn target_fabric_energy_efficiency(&self) -> f64;
    fn energy_demand(&self) -> &EnergyDemand;
    fn delivered_energy_use(&self) -> &DeliveredEnergyUse;
    fn energy_use_by_fuel(&self) -> &IndexMap<String, PerformanceValue>;
}

#[derive(Clone, Serialize)]
pub(super) struct EnergyDemand {
    space_heating: PerformanceValue,
    space_cooling: PerformanceValue,
}

impl<'a>
    From<(
        &'a CalculationResultsWithContext<'a>,
        &'a CalculationResultsWithContext<'a>,
        f64,
    )> for EnergyDemand
{
    fn from(
        (dwelling_results, target_results, total_floor_area): (
            &CalculationResultsWithContext,
            &CalculationResultsWithContext,
            f64,
        ),
    ) -> Self {
        EnergyDemand {
            space_heating: PerformanceValue {
                actual: dwelling_results.results.space_heat_demand_total() / total_floor_area,
                notional: target_results.results.space_heat_demand_total() / total_floor_area,
            },
            space_cooling: PerformanceValue {
                actual: dwelling_results.results.space_cool_demand_total() / total_floor_area,
                notional: target_results.results.space_cool_demand_total() / total_floor_area,
            },
        }
    }
}

#[derive(Clone, Serialize)]
pub(super) struct DeliveredEnergyUse {
    total: PerformanceValue,
    by_system: IndexMap<String, PerformanceValue>,
}

impl
    From<(
        &IndexMap<String, IndexMap<String, f64>>,
        &IndexMap<String, IndexMap<String, f64>>,
        f64,
    )> for DeliveredEnergyUse
{
    fn from(
        (dwelling_energy_use, target_energy_use, total_floor_area): (
            &IndexMap<String, IndexMap<String, f64>>,
            &IndexMap<String, IndexMap<String, f64>>,
            f64,
        ),
    ) -> Self {
        let total = PerformanceValue {
            actual: dwelling_energy_use
                .values()
                .map(|fuel_energy_use| fuel_energy_use["total"])
                .sum::<f64>()
                / total_floor_area,
            notional: target_energy_use
                .values()
                .map(|fuel_energy_use| fuel_energy_use["total"])
                .sum::<f64>()
                / total_floor_area,
        };
        let by_system = dwelling_energy_use
            .iter()
            .flat_map(|(_, energy_use)| energy_use.keys())
            .unique()
            .map(|key| {
                (key.clone(), {
                    let dwelling_use = dwelling_energy_use
                        .values()
                        .map(|fuel_energy_use| fuel_energy_use.get(key).unwrap_or(&0.))
                        .sum::<f64>()
                        / total_floor_area;
                    let target_use = target_energy_use
                        .values()
                        .map(|fuel_energy_use| fuel_energy_use.get(key).unwrap_or(&0.))
                        .sum::<f64>()
                        / total_floor_area;
                    PerformanceValue {
                        actual: dwelling_use,
                        notional: target_use,
                    }
                })
            })
            .collect::<IndexMap<_, _>>();

        DeliveredEnergyUse { total, by_system }
    }
}

#[derive(Clone, Copy, Serialize)]
pub(super) struct PerformanceValue {
    actual: f64,
    notional: f64,
}

pub(super) struct CalculatedComplianceResult {
    dwelling_final_rates: FinalRates,
    target_final_rates: FinalRates,
    dwelling_fabric_energy_efficiency: f64,
    target_fabric_energy_efficiency: f64,
    delivered_energy_use: DeliveredEnergyUse,
    energy_demand: EnergyDemand,
    energy_use_by_fuel: IndexMap<String, PerformanceValue>,
}

impl FhsComplianceCalculationResult for CalculatedComplianceResult {
    fn dwelling_emission_rate(&self) -> f64 {
        self.dwelling_final_rates.emission_rate
    }

    fn target_emission_rate(&self) -> f64 {
        self.target_final_rates.emission_rate
    }

    fn dwelling_primary_energy_rate(&self) -> f64 {
        self.dwelling_final_rates.primary_energy_rate
    }

    fn target_primary_energy_rate(&self) -> f64 {
        self.target_final_rates.primary_energy_rate
    }

    fn dwelling_fabric_energy_efficiency(&self) -> f64 {
        self.dwelling_fabric_energy_efficiency
    }

    fn target_fabric_energy_efficiency(&self) -> f64 {
        self.target_fabric_energy_efficiency
    }

    fn energy_demand(&self) -> &EnergyDemand {
        &self.energy_demand
    }

    fn delivered_energy_use(&self) -> &DeliveredEnergyUse {
        &self.delivered_energy_use
    }

    fn energy_use_by_fuel(&self) -> &IndexMap<String, PerformanceValue> {
        &self.energy_use_by_fuel
    }
}

impl TryFrom<&HashMap<CalculationKey, CalculationResultsWithContext<'_>>>
    for CalculatedComplianceResult
{
    type Error = anyhow::Error;

    fn try_from(
        results: &HashMap<CalculationKey, CalculationResultsWithContext>,
    ) -> Result<Self, Self::Error> {
        let dwelling_fhs_results = results
            .get(&CalculationKey::Fhs)
            .ok_or_else(|| anyhow!("Results were not available for the FHS calculation key."))?;
        let notional_fhs_results = results.get(&CalculationKey::FhsNotional).ok_or_else(|| {
            anyhow!("Results were not available for the FHS Notional calculation key.")
        })?;
        let dwelling_fhs_fee_results = results.get(&CalculationKey::FhsFee).ok_or_else(|| {
            anyhow!("Results were not available for the FHS FEE calculation key.")
        })?;
        let notional_fhs_fee_results =
            results
                .get(&CalculationKey::FhsNotionalFee)
                .ok_or_else(|| {
                    anyhow!("Results were not available for the FHS Notional FEE calculation key.")
                })?;

        let SummaryData {
            delivered_energy_map: dwelling_energy_use,
            ..
        } = build_summary_data(dwelling_fhs_results.try_into()?);
        let SummaryData {
            delivered_energy_map: target_energy_use,
            ..
        } = build_summary_data(notional_fhs_results.try_into()?);
        let total_floor_area = dwelling_fhs_results.context.corpus.total_floor_area();

        Ok(Self {
            dwelling_final_rates: calc_final_rates(
                dwelling_fhs_results.context.input,
                &dwelling_fhs_results.results.energy_import,
                &dwelling_fhs_results.results.energy_export,
                &dwelling_fhs_results.results.results_end_user,
                dwelling_fhs_results.results.timestep_array.len(),
            )?,
            target_final_rates: calc_final_rates(
                notional_fhs_results.context.input,
                &notional_fhs_results.results.energy_import,
                &notional_fhs_results.results.energy_export,
                &notional_fhs_results.results.results_end_user,
                notional_fhs_results.results.timestep_array.len(),
            )?,
            dwelling_fabric_energy_efficiency: calc_fabric_energy_efficiency(
                dwelling_fhs_fee_results.results.space_heat_demand_total(),
                dwelling_fhs_fee_results.results.space_cool_demand_total(),
                dwelling_fhs_fee_results.context.corpus.total_floor_area,
            ),
            target_fabric_energy_efficiency: calc_fabric_energy_efficiency(
                notional_fhs_fee_results.results.space_heat_demand_total(),
                notional_fhs_fee_results.results.space_cool_demand_total(),
                notional_fhs_fee_results.context.corpus.total_floor_area,
            ),
            delivered_energy_use: (&dwelling_energy_use, &target_energy_use, total_floor_area)
                .into(),
            energy_demand: (dwelling_fhs_results, notional_fhs_results, total_floor_area).into(),
            energy_use_by_fuel: dwelling_energy_use
                .keys()
                .map(|fuel| {
                    (fuel.clone(), {
                        let dwelling_fuel_total = dwelling_energy_use[fuel].values().sum::<f64>();
                        let target_fuel_total = target_energy_use[fuel].values().sum::<f64>();
                        PerformanceValue {
                            actual: dwelling_fuel_total,
                            notional: target_fuel_total,
                        }
                    })
                })
                .collect(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;
    use serde_json::json;
    use std::sync::LazyLock;

    #[fixture]
    fn canned_result() -> impl FhsComplianceCalculationResult {
        CannedResult {}
    }

    #[rstest]
    fn test_serialize_result_to_response(canned_result: impl FhsComplianceCalculationResult) {
        let expected_json = json!({
            "dwelling_emission_rate": 2.115,
            "target_emission_rate": 2.758,
            "emission_rate_compliant": true,
            "dwelling_primary_energy_rate": 61.190,
            "target_primary_energy_rate": 63.019,
            "primary_energy_rate_compliant": true,
            "dwelling_fabric_energy_efficiency": 41.487,
            "target_fabric_energy_efficiency": 43.529,
            "fabric_energy_efficiency_compliant": true,
            "energy_demand": {
                "space_heating": {
                    "actual": 34.962,
                    "notional": 42.550,
                },
                "space_cooling": {
                    "actual": 0.0,
                    "notional": 0.0,
                }
            },
            "delivered_energy_use": {
                "total": {
                    "actual": 73.750,
                    "notional": 54.866,
                },
                "by_system": {
                    "space_heating": {
                        "actual": 17.255,
                        "notional": 10.897,
                    },
                    "auxiliary": {
                        "actual": 4.54,
                        "notional": 1.904,
                    },
                    "water_heating": {
                        "actual": 20.624,
                        "notional": 13.293,
                    },
                    "electric_showers": {
                        "actual": 0.0,
                        "notional": 0.0,
                    },
                    "space_cooling": {
                        "actual": 0.0,
                        "notional": 0.0,
                    },
                    "ventilation": {
                        "actual": 1.336,
                        "notional": 0.0,
                    },
                    "lighting": {
                        "actual": 2.446,
                        "notional": 1.223,
                    },
                    "cooking": {
                        "actual": 4.769,
                        "notional": 4.769,
                    },
                    "appliances": {
                        "actual": 22.78,
                        "notional": 22.78,
                    },
                }
            },
            "energy_use_by_fuel": {
                "mains_electricity": {
                    "actual": 73.750,
                    "notional": 54.866,
                }
            }
        });

        assert_eq!(
            expected_json,
            serde_json::to_value(FhsComplianceResponse::build_from(&canned_result).unwrap())
                .unwrap()
        );
    }

    struct CannedResult;

    impl FhsComplianceCalculationResult for CannedResult {
        fn dwelling_emission_rate(&self) -> f64 {
            2.115
        }

        fn target_emission_rate(&self) -> f64 {
            2.758
        }

        fn dwelling_primary_energy_rate(&self) -> f64 {
            61.190
        }

        fn target_primary_energy_rate(&self) -> f64 {
            63.019
        }

        fn dwelling_fabric_energy_efficiency(&self) -> f64 {
            41.487
        }

        fn target_fabric_energy_efficiency(&self) -> f64 {
            43.529
        }

        fn energy_demand(&self) -> &EnergyDemand {
            &EnergyDemand {
                space_heating: PerformanceValue {
                    actual: 34.962,
                    notional: 42.550,
                },
                space_cooling: PerformanceValue {
                    actual: 0.0,
                    notional: 0.0,
                },
            }
        }

        fn delivered_energy_use(&self) -> &DeliveredEnergyUse {
            &DELIVERED_ENERGY_USE
        }

        fn energy_use_by_fuel(&self) -> &IndexMap<String, PerformanceValue> {
            &ENERGY_USE_BY_FUEL
        }
    }

    static DELIVERED_ENERGY_USE: LazyLock<DeliveredEnergyUse> =
        LazyLock::new(|| DeliveredEnergyUse {
            total: PerformanceValue {
                actual: 73.750,
                notional: 54.866,
            },
            by_system: IndexMap::from(
                [
                    ("space_heating", 17.255, 10.897),
                    ("auxiliary", 4.54, 1.904),
                    ("water_heating", 20.624, 13.293),
                    ("electric_showers", 0.0, 0.0),
                    ("space_cooling", 0.0, 0.0),
                    ("ventilation", 1.336, 0.0),
                    ("lighting", 2.446, 1.223),
                    ("cooking", 4.769, 4.769),
                    ("appliances", 22.78, 22.78),
                ]
                .map(|(k, actual, notional)| (k.into(), PerformanceValue { actual, notional })),
            ),
        });

    static ENERGY_USE_BY_FUEL: LazyLock<IndexMap<String, PerformanceValue>> = LazyLock::new(|| {
        IndexMap::from(
            [("mains_electricity", 73.750, 54.866)]
                .map(|(k, actual, notional)| (k.into(), PerformanceValue { actual, notional })),
        )
    });
}
