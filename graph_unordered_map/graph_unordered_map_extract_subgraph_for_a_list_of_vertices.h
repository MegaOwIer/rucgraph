#pragma once

#include <graph_unordered_map/graph_unordered_map.h> 

graph_unordered_map graph_unordered_map_extract_subgraph_for_a_list_of_vertices(
	graph_unordered_map& input_graph, std::list<int>& list_of_v) {

	/*extract a smaller_graph, which contains all the vertices in list_of_v,
	and all the edges between vertices in list_of_v;
	time complexity O(|V_list|+|adj_v of V_list in input_graph|)*/



	/*time complexity O(|V_list|)*/
	std::unordered_set<int> list_unordered_set;
	for (auto i = list_of_v.begin(); i != list_of_v.end(); i++) {
		list_unordered_set.insert(*i);
	}

	/*time complexity O(|V_list|+|adj_v of V_list in input_graph|)*/
	graph_unordered_map smaller_graph;
	for (auto i = list_of_v.begin(); i != list_of_v.end(); i++) {
		int v1 = *i;
		double nw = input_graph[v1].first;
		graph_unordered_map_add_vertex(smaller_graph, v1, nw); // add vertex
		for (auto it = input_graph[v1].second.begin(); it != input_graph[v1].second.end(); it++) {
			int v2 = it->first;
			double ec = it->second;
			if (list_unordered_set.count(v2) > 0 & v2 > v1) { // v2 is in the list and only add edge once
				graph_unordered_map_add_edge(smaller_graph, v1, v2, ec); // add edge
			}

		}
	}

	return smaller_graph;

}