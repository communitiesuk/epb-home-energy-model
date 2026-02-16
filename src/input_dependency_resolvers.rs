use crate::input::HotWaterSourceDetails;
use indexmap::IndexMap;
use petgraph::Graph;
use smartstring::alias::String;

/// Build a dependency graph for PreHeatedWaterSource objects.
pub(crate) fn build_preheated_water_source_dependency_graph(
    preheated_sources_input: IndexMap<String, HotWaterSourceDetails>,
) -> Graph<String, String> {
    let mut graph = Graph::<String, String>::new();
    let mut nodes = IndexMap::new();

    for name in preheated_sources_input.keys() {
        let node_index = graph.add_node(name.clone());
        nodes.insert(name, node_index);
    }

    let mut edges = Vec::new();

    for (source_name, source_details) in &preheated_sources_input {
        let cold_source_name = source_details.cold_water_source();
        if preheated_sources_input.contains_key(cold_source_name.as_str()) {
            edges.push((nodes[&cold_source_name], nodes[source_name]));
        }
    }
    graph.extend_with_edges(&edges);

    graph
}

/// Returns a list of source names in initialization order (dependencies first)
pub(crate) fn topological_sort_preheated_water_sources(_dependency_graph: Graph<String, String>) {

    // try:
    //     sorter = TopologicalSorter(dependency_graph)
    // return list(sorter.static_order())
    // except CycleError as err:
    //     raise ValueError(f"Circular dependency detected in PreHeatedWaterSource: {err}") from err
}
