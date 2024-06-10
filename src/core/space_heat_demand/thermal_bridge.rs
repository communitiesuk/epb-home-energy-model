use indexmap::IndexMap;
use serde_json::Value;

#[derive(Debug, PartialEq)]
pub enum ThermalBridging {
    Number(f64),
    Bridges(IndexMap<String, ThermalBridge>),
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum ThermalBridge {
    Linear {
        linear_thermal_transmittance: f64,
        length: f64,
    },
    Point {
        heat_transfer_coefficient: f64,
    },
}

pub fn heat_transfer_coefficient_for_thermal_bridge(thermal_bridge: &ThermalBridge) -> f64 {
    match *thermal_bridge {
        ThermalBridge::Linear {
            linear_thermal_transmittance: t,
            length: l,
        } => t * l,
        ThermalBridge::Point {
            heat_transfer_coefficient: h,
        } => h,
    }
}

/// Converts the serde_json value for ThermalBridging into structs as have not been able to handle these variants
/// Ultimately this field in the JSON should be divided into different attributes (one for bridging elements, one for a number)
/// and then this would not be necessary (as serde_json would handle it) and could be removed!
pub fn thermal_bridging_from_input(input: Value) -> ThermalBridging {
    match input {
        Value::Object(map) => {
            let mut bridges = IndexMap::new();
            for (name, bridge_object) in map.clone().iter() {
                bridges.insert(
                    (*name).clone(),
                    match bridge_object.get("type") {
                        Some(Value::String(s)) if s == &String::from("ThermalBridgeLinear") => {
                            ThermalBridge::Linear {
                                linear_thermal_transmittance: bridge_object
                                    .get("linear_thermal_transmittance")
                                    .unwrap()
                                    .as_f64()
                                    .unwrap(),
                                length: bridge_object.get("length").unwrap().as_f64().unwrap(),
                            }
                        }
                        Some(Value::String(s)) if s == &String::from("ThermalBridgePoint") => {
                            ThermalBridge::Point {
                                heat_transfer_coefficient: bridge_object
                                    .get("heat_transfer_coeff")
                                    .unwrap()
                                    .as_f64()
                                    .unwrap(),
                            }
                        }
                        _ => panic!("unknown type of thermal bridging sent"),
                    },
                );
            }

            ThermalBridging::Bridges(bridges)
        }
        Value::Number(n) => ThermalBridging::Number(n.as_f64().unwrap()),
        _ => panic!("thermal bridging input was not in an expected format"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_ulps_eq;
    use rstest::*;
    use serde_json::json;

    #[rstest]
    pub fn should_convert_serde_json_value_into_thermal_bridging() {
        assert_eq!(
            thermal_bridging_from_input(json!(2.4)),
            ThermalBridging::Number(2.4)
        );
        assert_eq!(
            thermal_bridging_from_input(json!(
                {"TB1": {
                    "type": "ThermalBridgeLinear",
                    "linear_thermal_transmittance": 1.0,
                    "length": 3.0
                },
                "TB2": {
                    "type": "ThermalBridgeLinear",
                    "linear_thermal_transmittance": 0.1,
                    "length": 2.0
                },
                "TB3": {
                    "type": "ThermalBridgePoint",
                    "heat_transfer_coeff": 2.0
                }}
            )),
            ThermalBridging::Bridges(IndexMap::from([
                (
                    "TB1".to_string(),
                    ThermalBridge::Linear {
                        linear_thermal_transmittance: 1.0,
                        length: 3.0
                    }
                ),
                (
                    "TB2".to_string(),
                    ThermalBridge::Linear {
                        linear_thermal_transmittance: 0.1,
                        length: 2.0
                    }
                ),
                (
                    "TB3".to_string(),
                    ThermalBridge::Point {
                        heat_transfer_coefficient: 2.0
                    }
                )
            ]))
        );
    }

    #[rstest]
    pub fn should_have_correct_heat_transfer_coefficient_for_thermal_bridge() {
        assert_eq!(
            heat_transfer_coefficient_for_thermal_bridge(&ThermalBridge::Linear {
                linear_thermal_transmittance: 1.2,
                length: 4.0,
            }),
            4.8
        );
        assert_eq!(
            heat_transfer_coefficient_for_thermal_bridge(&ThermalBridge::Point {
                heat_transfer_coefficient: 6.7
            }),
            6.7
        );
    }

    #[rstest]
    pub fn test_heat_trans_coeff_for_linear() {
        let tb = ThermalBridge::Linear {
            linear_thermal_transmittance: 0.28,
            length: 5.0,
        };
        assert_ulps_eq!(heat_transfer_coefficient_for_thermal_bridge(&tb), 1.4,);
    }

    #[rstest]
    pub fn test_heat_trans_coeff_for_point() {
        let tb = ThermalBridge::Point {
            heat_transfer_coefficient: 1.4,
        };
        assert_eq!(
            heat_transfer_coefficient_for_thermal_bridge(&tb),
            1.4,
            "incorrect heat transfer coefficient returned"
        );
    }
}
