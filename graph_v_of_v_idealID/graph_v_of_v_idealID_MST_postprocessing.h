
#pragma once


#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_minimum_spanning_tree.h>
#include <graph_v_of_v_idealID/graph_v_of_v_idealID_extract_subgraph.h>
graph_hash_of_mixed_weighted graph_v_of_v_idealID_MST_postprocessing(graph_v_of_v_idealID& input_graph, unordered_set<int>& V_set) {

	/*time complexity O(|subgraph_V|+|adj_v of subgraph_V in input_graph|+|subgraph_E|+|subgraph_V|log|subgraph_V|); nw and ec in theta will be the same with input_graph*/

	graph_hash_of_mixed_weighted postprocessed_theta; // nw and ec will be the same with input_graph, may not theta

	/*time complexity O(|subgraph_V)*/
	for (auto it = V_set.begin(); it != V_set.end(); it++) {
		graph_hash_of_mixed_weighted_add_vertex(postprocessed_theta, *it, 0); // insert vertex
	}

	/*time complexity O(|subgraph_V|+|adj_v of subgraph_V in input_graph|)*/
	graph_hash_of_mixed_weighted subgraph_theta = graph_v_of_v_idealID_extract_subgraph(input_graph, V_set);

	//graph_hash_of_mixed_weighted_print(subgraph_theta);

	/*time complexity: O(|subgraph_E|+|subgraph_V|log|subgraph_V|)*/
	std::unordered_map<int, int> predesscors = graph_hash_of_mixed_weighted_minimum_spanning_tree(subgraph_theta);
	for (auto it = predesscors.begin(); it != predesscors.end(); it++) {
		int v1 = it->first, v2 = it->second;
		if (v1 != v2) {
			graph_hash_of_mixed_weighted_add_edge(postprocessed_theta, v1, v2, graph_v_of_v_idealID_edge_weight(input_graph, v1, v2));
		}
	}

	return postprocessed_theta;

}






