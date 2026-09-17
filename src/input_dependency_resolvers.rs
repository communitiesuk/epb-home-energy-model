use crate::input::PreHeatedWaterSourceDetails;
use arcstr::ArcStr;
use indexmap::IndexMap;
use petgraph::algo::toposort;
use petgraph::Graph;
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
    error_context: String,
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

#[derive(Debug, Error, PartialEq)]
#[error("Circular dependency detected in {0}")]
pub(crate) struct CircularDependencyError(String);

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn test_preheated_water_sources_are_reordered() {
        let hot_water_source_details: PreHeatedWaterSourceDetails = serde_json::from_value(json!(
        {"type": "StorageTank",
        "volume": 24.0,
        "daily_losses": 1.55,
        "init_temp": 48.0,
        "ColdWaterSource": "storagetank2",
        "HeatSource": {
            "{name}_immersion": {
                "type": "ImmersionHeater",
                "power": 3.0,
                "EnergySupply": "mains elec",
                "Controlmin": "min_temp",
                "Controlmax": "setpoint_temp_max",
                "heater_position": 0.3,
                "thermostat_position": 0.33}}
            }))
        .unwrap();
        let hot_water_source_details_2: PreHeatedWaterSourceDetails =
            serde_json::from_value(json!(
            {"type": "StorageTank",
            "volume": 24.0,
            "daily_losses": 1.55,
            "init_temp": 48.0,
            "ColdWaterSource": "mains water",
            "HeatSource": {
                "{name}_immersion": {
                    "type": "ImmersionHeater",
                    "power": 3.0,
                    "EnergySupply": "mains elec",
                    "Controlmin": "min_temp",
                    "Controlmax": "setpoint_temp_max",
                    "heater_position": 0.3,
                    "thermostat_position": 0.33}}
                }))
            .unwrap();
        let preheated_sources_input: IndexMap<ArcStr, PreHeatedWaterSourceDetails> =
            IndexMap::from([
                ("storagetank1".into(), hot_water_source_details),
                ("storagetank2".into(), hot_water_source_details_2),
            ]);

        let graph = build_preheated_water_source_dependency_graph(&preheated_sources_input);
        let result = topological_sort_preheated_water_sources(&graph).unwrap();

        assert_eq!(result, vec!["storagetank2", "storagetank1"]);
    }

    #[test]
    fn test_preheated_water_sources_with_circular_dependency_produces_error() {
        let hot_water_source_details: PreHeatedWaterSourceDetails = serde_json::from_value(json!(
        {"type": "StorageTank",
        "volume": 24.0,
        "daily_losses": 1.55,
        "init_temp": 48.0,
        "ColdWaterSource": "storagetank2",
        "HeatSource": {
            "{name}_immersion": {
                "type": "ImmersionHeater",
                "power": 3.0,
                "EnergySupply": "mains elec",
                "Controlmin": "min_temp",
                "Controlmax": "setpoint_temp_max",
                "heater_position": 0.3,
                "thermostat_position": 0.33}}
            }))
        .unwrap();
        let hot_water_source_details_2: PreHeatedWaterSourceDetails =
            serde_json::from_value(json!(
            {"type": "StorageTank",
            "volume": 24.0,
            "daily_losses": 1.55,
            "init_temp": 48.0,
            "ColdWaterSource": "storagetank1",
            "HeatSource": {
                "{name}_immersion": {
                    "type": "ImmersionHeater",
                    "power": 3.0,
                    "EnergySupply": "mains elec",
                    "Controlmin": "min_temp",
                    "Controlmax": "setpoint_temp_max",
                    "heater_position": 0.3,
                    "thermostat_position": 0.33}}
                }))
            .unwrap();
        let preheated_sources_input: IndexMap<ArcStr, PreHeatedWaterSourceDetails> =
            IndexMap::from([
                ("storagetank1".into(), hot_water_source_details),
                ("storagetank2".into(), hot_water_source_details_2),
            ]);

        let graph = build_preheated_water_source_dependency_graph(&preheated_sources_input);
        let result = topological_sort_preheated_water_sources(&graph);

        assert!(result.is_err());
    }
}
