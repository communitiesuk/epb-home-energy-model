#![allow(dead_code)]

use indexmap::IndexMap;
use serde::Serialize;

#[derive(Serialize)]
struct FhsComplianceResponse {
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
    fn build_from<T>(result: &T) -> anyhow::Result<Self>
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
            energy_demand: result.energy_demand(),
            delivered_energy_use: result.delivered_energy_use(),
            energy_use_by_fuel: result.energy_use_by_fuel(),
        })
    }
}

/// A trait to define data on the underlying calculation result that an FHS compliance result requires.
trait FhsComplianceCalculationResult {
    fn dwelling_emission_rate(&self) -> f64;
    fn target_emission_rate(&self) -> f64;
    fn dwelling_primary_energy_rate(&self) -> f64;
    fn target_primary_energy_rate(&self) -> f64;
    fn dwelling_fabric_energy_efficiency(&self) -> f64;
    fn target_fabric_energy_efficiency(&self) -> f64;
    fn energy_demand(&self) -> EnergyDemand;
    fn delivered_energy_use(&self) -> DeliveredEnergyUse;
    fn energy_use_by_fuel(&self) -> IndexMap<String, PerformanceValue>;
}

#[derive(Serialize)]
struct EnergyDemand {
    space_heating: PerformanceValue,
    space_cooling: PerformanceValue,
}

#[derive(Serialize)]
struct DeliveredEnergyUse {
    total: PerformanceValue,
    by_system: IndexMap<String, PerformanceValue>,
}

#[derive(Clone, Copy, Serialize)]
struct PerformanceValue {
    actual: f64,
    notional: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;
    use serde_json::json;

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

        fn energy_demand(&self) -> EnergyDemand {
            EnergyDemand {
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

        fn delivered_energy_use(&self) -> DeliveredEnergyUse {
            DeliveredEnergyUse {
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
                    .map(|(k, actual, notional)| {
                        (k.to_owned(), PerformanceValue { actual, notional })
                    }),
                ),
            }
        }

        fn energy_use_by_fuel(&self) -> IndexMap<String, PerformanceValue> {
            IndexMap::from(
                [("mains_electricity", 73.750, 54.866)].map(|(k, actual, notional)| {
                    (k.to_owned(), PerformanceValue { actual, notional })
                }),
            )
        }
    }
}