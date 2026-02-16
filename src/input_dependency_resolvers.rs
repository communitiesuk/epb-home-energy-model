use crate::input::HotWaterSourceDetails;
use indexmap::IndexMap;
use petgraph::Graph;
use smartstring::alias::String;
use std::borrow::Borrow;
use std::hash::Hash;

/// Build a dependency graph for PreHeatedWaterSource objects.
pub(crate) fn build_preheated_water_source_dependency_graph<T>(
    preheated_sources_input: IndexMap<T, HotWaterSourceDetails>,
) -> Graph<T, T>
where
    T: Hash + Eq + Borrow<str> + AsRef<str> + Clone + Default,
{
    let mut graph = Graph::<T, T>::new();
    let mut nodes = IndexMap::new();

    for name in preheated_sources_input.keys() {
        let node_index = graph.add_node(name.clone());
        nodes.insert(name, node_index);
    }

    let mut edges = Vec::new();

    for (source_name, source_details) in &preheated_sources_input {
        let cold_source_name = source_details.cold_water_source();
        if preheated_sources_input.contains_key(cold_source_name) {
            edges.push((nodes[source_name], nodes[source_name]));
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
