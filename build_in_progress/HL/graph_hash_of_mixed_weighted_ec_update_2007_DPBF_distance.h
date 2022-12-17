#pragma once
#include <cmath> 
void graph_hash_of_mixed_weighted_ec_update_2007_DPBF_distance(graph_hash_of_mixed_weighted& input_graph, int max_N) {

	vector<int> degrees(max_N);
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		int v = it->first;
		degrees[v] = input_graph.degree(v);
	}


	/*time complexity: O(|V|+|E|)*/
	for (auto it1 = input_graph.hash_of_vectors.begin(); it1 != input_graph.hash_of_vectors.end(); it1++) {
		int v1 = it1->first;
		std::vector<int> adj_v = input_graph.adj_v(v1);
		for (int i = adj_v.size() - 1; i >= 0; i--) {
			int v2 = adj_v[i];
			double ec = log2(1 + (degrees[v1] > degrees[v2] ? degrees[v1] : degrees[v2]));
			graph_hash_of_mixed_weighted_add_edge(input_graph, v1, v2, ec);// update_pairwise_jaccard_distance for edge (j,i)
		}
	}

}