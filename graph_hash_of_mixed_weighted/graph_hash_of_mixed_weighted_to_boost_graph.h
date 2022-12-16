#pragma once

#include <boost_graph/boost_graph.h>
graph graph_hash_of_mixed_weighted_to_boost_graph(graph_hash_of_mixed_weighted& input_graph) {

	/*this function transforms a graph_hash_of_mixed_weighted to a boost graph;
	we assume that vertices in input_graph are 0,...,N-1;
	the time complexity depends on the boost graph structure*/

	int N = input_graph.hash_of_vectors.size();
	graph boost_graph(N);

	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		int v = it->first;
		double nw = it->second.vertex_weight;
		boost::put(boost::vertex_name_t(), boost_graph, v, nw); // put node weight

		auto search = input_graph.hash_of_hashs.find(v);
		if (search != input_graph.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int vertex = it2->first;
				if (vertex < v) {
					double ec = it2->second;
					boost::add_edge(v, vertex, ec, boost_graph);
				}
			}
		}
		else {
			auto search2 = input_graph.hash_of_vectors.find(v);
			for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
				int vertex = it2->first;
				if (vertex < v) {
					double ec = it2->second;
					boost::add_edge(v, vertex, ec, boost_graph);
				}
			}
		}
	}

	return boost_graph;
}