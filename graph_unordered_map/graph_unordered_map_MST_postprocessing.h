#pragma once


#include <graph_unordered_map.h> 
#include <graph_unordered_map_extract_subgraph_for_a_hash_of_vertices.h> 

graph_unordered_map graph_unordered_map_MST_postprocessing(graph_unordered_map& input_graph, graph_unordered_map& theta) {

	/*time complexity O(|subgraph_V|+|adj_v of subgraph_V in input_graph|+|subgraph_E|+|subgraph_V|log|subgraph_V|)*/

	graph_unordered_map postprocessed_theta; // nw and ec will be the same with input_graph, may not theta

	/*time complexity O(|subgraph_V|)*/
	std::unordered_set<int> hash_of_v;
	for (auto it = theta.begin(); it != theta.end(); it++) {
		hash_of_v.insert(it->first);
		graph_unordered_map_add_vertex(postprocessed_theta, it->first, input_graph[it->first].first); // insert vertex
	}

	/*time complexity O(|subgraph_V|+|adj_v of subgraph_V in input_graph|)*/
	graph_unordered_map subgraph_theta = graph_unordered_map_extract_subgraph_for_a_hash_of_vertices(
		input_graph, hash_of_v);

	/*time complexity: O(|subgraph_E|+|subgraph_V|log|subgraph_V|)*/
	std::unordered_map<int, int> predesscors = graph_unordered_map_minimum_spanning_tree(subgraph_theta);
	for (auto it = predesscors.begin(); it != predesscors.end(); it++) {
		int v1 = it->first, v2 = it->second;
		if (v1 != v2) { // predesscor of start_v is start_v
			graph_unordered_map_add_edge(postprocessed_theta, v1, v2, input_graph[v1].second[v2]);
		}
	}

	return postprocessed_theta;

}



graph_unordered_map graph_unordered_map_MST_postprocessing_hash_of_vertices
(graph_unordered_map& input_graph, std::unordered_set<int>& hash_of_v) {

	/*time complexity O(|subgraph_V|+|adj_v of subgraph_V in input_graph|+|subgraph_E|+|subgraph_V|log|subgraph_V|)*/

	graph_unordered_map postprocessed_theta; // nw and ec will be the same with input_graph, may not theta

	/*time complexity O(|subgraph_V)*/
	for (auto it = hash_of_v.begin(); it != hash_of_v.end(); it++) {
		graph_unordered_map_add_vertex(postprocessed_theta, *it, input_graph[*it].first); // insert vertex
	}

	/*time complexity O(|subgraph_V|+|adj_v of subgraph_V in input_graph|)*/
	graph_unordered_map subgraph_theta = graph_unordered_map_extract_subgraph_for_a_hash_of_vertices(
		input_graph, hash_of_v);

	//graph_unordered_map_print(subgraph_theta);

	/*time complexity: O(|subgraph_E|+|subgraph_V|log|subgraph_V|)*/
	std::unordered_map<int, int> predesscors = graph_unordered_map_minimum_spanning_tree(subgraph_theta);
	for (auto it = predesscors.begin(); it != predesscors.end(); it++) {
		int v1 = it->first, v2 = it->second;
		if (v1 != v2) {
			graph_unordered_map_add_edge(postprocessed_theta, v1, v2, input_graph[v1].second[v2]);
		}
	}

	return postprocessed_theta;

}
