
#pragma once

void graph_hash_of_mixed_weighted_distribute_nw_into_ec_evenly(graph_hash_of_mixed_weighted& input_graph) {

	/*time complexity: O(|V|+|E|)*/

	for (auto it1 = input_graph.hash_of_vectors.begin(); it1 != input_graph.hash_of_vectors.end(); it1++) {
		int i = it1->first;
		double w_i = it1->second.vertex_weight;
		double degree = graph_hash_of_mixed_weighted_degree(input_graph, i);
		double divided_w_i = w_i / degree;


		auto search = input_graph.hash_of_hashs.find(i);
		if (search != input_graph.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int j = it2->first;
				it2->second = it2->second + divided_w_i; // evenly distribute nw into ec
			}
		}
		else {
			auto search2 = input_graph.hash_of_vectors.find(i);
			for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
				int j = it2->first;
				it2->second = it2->second + divided_w_i; // evenly distribute nw into ec
			}
		}
	}

}


void graph_hash_of_mixed_weighted_distribute_nw_into_ec_evenly_reverse(graph_hash_of_mixed_weighted& input_graph) {

	/*time complexity: O(|V|+|E|)*/

	for (auto it1 = input_graph.hash_of_vectors.begin(); it1 != input_graph.hash_of_vectors.end(); it1++) {
		int i = it1->first;
		double w_i = it1->second.vertex_weight;
		double degree = graph_hash_of_mixed_weighted_degree(input_graph, i);
		double divided_w_i = w_i / degree;

		auto search = input_graph.hash_of_hashs.find(i);
		if (search != input_graph.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int j = it2->first;
				it2->second = it2->second - divided_w_i; // remove the evenly distribute nw from ec
			}
		}
		else {
			auto search2 = input_graph.hash_of_vectors.find(i);
			for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
				int j = it2->first;
				it2->second = it2->second - divided_w_i; // remove the evenly distribute nw from ec
			}
		}
	}

}

