#pragma once


#include <graph_unordered_map.h> 

void graph_unordered_map_copy_graph_to_another_graph(graph_unordered_map& target_graph, graph_unordered_map& copied_graph) {

	/*time complexity: O(|V_copied_graph| +|E_copied_graph|);
	this is to merge copied_graph into target_graph*/

	for (auto it1 = copied_graph.begin(); it1 != copied_graph.end(); it1++) {
		int i = it1->first;
		double w_i = it1->second.first;
		graph_unordered_map_add_vertex(target_graph, i, w_i);
		for (auto it2 = it1->second.second.begin(); it2 != it1->second.second.end(); it2++) {
			int j = it2->first;
			if (i < j) {
				double c_ij = it2->second;
				graph_unordered_map_add_edge(target_graph, i, j, c_ij);
			}
		}
	}
}