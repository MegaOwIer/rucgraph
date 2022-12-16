#pragma once


#include <graph_unordered_map.h> 

void graph_unordered_map_copy_weights_of_graph1_to_graph2(graph_unordered_map& graph1, graph_unordered_map& graph2) {

	/*we assume that every vertex and edge of graph2 is in graph1*/

	/*time complexity: O(|V_graph2| +|E_graph2|)*/

	for (auto it1 = graph2.begin(); it1 != graph2.end(); it1++) {
		int i = it1->first;
		it1->second.first = graph1[i].first; // copy nw
		for (auto it2 = it1->second.second.begin(); it2 != it1->second.second.end(); it2++) {
			int j = it2->first;
			it2->second = graph1[i].second[j]; // copy ec
		}
	}
}