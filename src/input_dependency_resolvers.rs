use crate::core::common::WaterSupplyBehaviour;
use crate::core::heating_systems::heat_battery_pcm::HeatBatteryChargingSource;
use crate::input::{HeatBatteryPcmChargingSource, PreHeatedWaterSourceDetails};
use anyhow::bail;
use arcstr::ArcStr;
use indexmap::IndexMap;
use petgraph::algo::toposort;
use petgraph::Graph;
use std::collections::HashSet;
use thiserror::Error;

/// Build a dependency graph for PreHeatedWaterSource objects.
pub(crate) fn build_preheated_water_source_dependency_graph(
    preheated_sources_input: &IndexMap<ArcStr, PreHeatedWaterSourceDetails>,
) -> Graph<ArcStr, ArcStr> {
    let mut graph = Graph::<ArcStr, ArcStr>::new();
    let mut nodes: IndexMap<ArcStr, _> = IndexMap::new();

    for name in preheated_sources_input.keys() {
        let node_index = graph.add_node(name.into());
        nodes.insert(name.into(), node_index);
    }

    let mut edges = Vec::new();

    for (source_name, source_details) in preheated_sources_input {
        let cold_source_name = source_details.cold_water_source();
        if preheated_sources_input.contains_key(cold_source_name) {
            edges.push((nodes[cold_source_name], nodes[source_name]));
        }
    }
    graph.extend_with_edges(&edges);

    graph
}

/// Returns a list of names in dependency order (dependencies first).
fn topological_sort<T: Clone>(
    dependency_graph: &Graph<T, T>,
    error_context: ArcStr,
) -> Result<Vec<T>, CircularDependencyError> {
    match toposort(dependency_graph, None) {
        Ok(ordered_nodes) => Ok(ordered_nodes
            .into_iter()
            .map(|node| dependency_graph[node].clone())
            .collect()),
        Err(_) => Err(CircularDependencyError(error_context)),
    }
}

/// Topological sort for PreHeatedWaterSource dependencies.
pub(crate) fn topological_sort_preheated_water_sources<T: Clone>(
    dependency_graph: &Graph<T, T>,
) -> Result<Vec<T>, CircularDependencyError> {
    topological_sort(dependency_graph, "PreHeatedWaterSource".into())
}

/// Build a dependency graph for PCM heat battery hydronic charging.
pub(crate) fn build_heat_battery_charging_dependency_graph<T: WaterSupplyBehaviour>(
    pcm_hydronic_charging_pending: &[(
        ArcStr,
        HeatBatteryChargingSource<T>,
        HeatBatteryPcmChargingSource,
    )],
) -> anyhow::Result<Graph<ArcStr, ArcStr>> {
    let mut graph = Graph::<ArcStr, ArcStr>::new();
    let mut nodes: IndexMap<ArcStr, _> = IndexMap::new();

    let battery_names: HashSet<&ArcStr> = pcm_hydronic_charging_pending
        .iter()
        .map(|(name, _, _)| name)
        .collect();
    for name in &battery_names {
        let node_index = graph.add_node((*name).into());
        nodes.insert((*name).into(), node_index);
    }

    let mut edges = Vec::new();

    for (battery_name, _, src_data) in pcm_hydronic_charging_pending.iter() {
        let charging_source_name = match &src_data {
            HeatBatteryPcmChargingSource::Hydronic { name, .. } => name,
            _ => bail!("Hydronic charging source expected for a PCM heat battery"),
        };
        if battery_names.contains(&charging_source_name) {
            edges.push((nodes[charging_source_name], nodes[battery_name]));
        }
    }

    graph.extend_with_edges(&edges);

    Ok(graph)
}

/// Topological sort for heat battery charging dependencies.
pub(crate) fn topological_sort_heat_battery_charging<T: Clone>(
    dependency_graph: &Graph<T, T>,
) -> Result<Vec<T>, CircularDependencyError> {
    topological_sort(dependency_graph, "heat battery charging".into())
}

#[derive(Debug, Error, PartialEq)]
#[error("Circular dependency detected in {0}")]
pub(crate) struct CircularDependencyError(ArcStr);

#[cfg(test)]
mod tests {
    use super::*;
    use indexmap::indexmap;
    use serde_json::json;

    fn hot_water_source_details(cold_water_source: ArcStr) -> PreHeatedWaterSourceDetails {
        serde_json::from_value(json!(
            {
                "type": "StorageTank",
                "volume": 24.0,
                "daily_losses": 1.55,
                "init_temp": 48.0,
                "ColdWaterSource": cold_water_source,
                "HeatSource": {
                    "{name}_immersion": {
                        "type": "ImmersionHeater",
                        "power": 3.0,
                        "EnergySupply": "mains elec",
                        "Controlmin": "min_temp",
                        "Controlmax": "setpoint_temp_max",
                        "heater_position": 0.3,
                        "thermostat_position": 0.33
                    }
                }
            }
        ))
        .unwrap()
    }

    mod test_build_preheated_water_source_dependency_graph {
        use super::*;

        #[test]
        /// Sources with cold water feeds have no predecessors.
        fn test_no_dependencies() {
            let sources: IndexMap<ArcStr, PreHeatedWaterSourceDetails> = IndexMap::from([
                (
                    "preheat_A".into(),
                    hot_water_source_details("mains_cold".into()),
                ),
                (
                    "preheat_B".into(),
                    hot_water_source_details("mains_cold".into()),
                ),
            ]);

            let graph = build_preheated_water_source_dependency_graph(&sources);

            assert_eq!(graph.edge_count(), 0);
            assert_eq!(graph.node_count(), 2);
        }

        #[test]
        /// Source B feeds into source A — A depends on B.
        fn test_chain_dependency() {
            let sources: IndexMap<ArcStr, PreHeatedWaterSourceDetails> = IndexMap::from([
                (
                    "preheat_A".into(),
                    hot_water_source_details("preheat_B".into()),
                ),
                (
                    "preheat_B".into(),
                    hot_water_source_details("mains_cold".into()),
                ),
            ]);

            let graph = build_preheated_water_source_dependency_graph(&sources);
            let a_idx = graph
                .node_indices()
                .find(|&idx| graph[idx] == "preheat_A")
                .unwrap();
            let b_idx = graph
                .node_indices()
                .find(|&idx| graph[idx] == "preheat_B")
                .unwrap();

            assert_eq!(graph.edge_count(), 1);
            assert_eq!(graph.node_count(), 2);
            // A depends on B (edge flows from B to A)
            assert!(graph.contains_edge(b_idx, a_idx));
            // B doesn't depend on A (no edge from A to B)
            assert!(!graph.contains_edge(a_idx, b_idx));
        }

        #[test]
        /// Empty input returns empty graph.
        fn test_empty_map() {
            let graph = build_preheated_water_source_dependency_graph(&indexmap! {});

            assert_eq!(graph.edge_count(), 0);
            assert_eq!(graph.node_count(), 0);
        }
    }

    mod test_topological_sort_preheated_water_sources {
        use super::*;

        #[test]
        /// B must come before A when A depends on B.
        fn test_chain_order() {
            let sources: IndexMap<ArcStr, PreHeatedWaterSourceDetails> = IndexMap::from([
                (
                    "preheat_A".into(),
                    hot_water_source_details("preheat_B".into()),
                ),
                (
                    "preheat_B".into(),
                    hot_water_source_details("mains_cold".into()),
                ),
            ]);

            let graph = build_preheated_water_source_dependency_graph(&sources);

            let result = topological_sort_preheated_water_sources(&graph).unwrap();

            assert_eq!(result, ["preheat_B", "preheat_A"]);
        }

        #[test]
        /// Circular dependency produces an error
        fn test_circular_dependency_raises() {
            let sources: IndexMap<ArcStr, PreHeatedWaterSourceDetails> = IndexMap::from([
                (
                    "preheat_A".into(),
                    hot_water_source_details("preheat_B".into()),
                ),
                (
                    "preheat_B".into(),
                    hot_water_source_details("preheat_A".into()),
                ),
            ]);

            let graph = build_preheated_water_source_dependency_graph(&sources);

            let result = topological_sort_preheated_water_sources(&graph);

            assert!(result.is_err());
        }

        #[test]
        /// Independent sources are all returned (order doesn't matter).
        fn test_independent_sources() {
            let sources: IndexMap<ArcStr, PreHeatedWaterSourceDetails> = IndexMap::from([
                ("A".into(), hot_water_source_details("mains_cold".into())),
                ("B".into(), hot_water_source_details("mains_cold".into())),
                ("C".into(), hot_water_source_details("mains_cold".into())),
            ]);

            let graph = build_preheated_water_source_dependency_graph(&sources);

            let mut result = topological_sort_preheated_water_sources(&graph).unwrap();
            result.sort();

            assert_eq!(result, ["A", "B", "C"]);
        }
    }
}
