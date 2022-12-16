#pragma once


#include <graph_unordered_map.h> 

void graph_unordered_map_nw_ec_normalization(graph_unordered_map& input_graph) {

	/*this is "min-max normalization" for bringing all node weights to [0,1];
	https://en.wikipedia.org/wiki/Normalization_(statistics)
	time complexity: O(|V|+|E|)*/

	double min_nw = std::numeric_limits<double>::max();
	double max_nw = std::numeric_limits<double>::min();
	double min_ec = std::numeric_limits<double>::max();
	double max_ec = std::numeric_limits<double>::min();

	for (auto it1 = input_graph.begin(); it1 != input_graph.end(); it1++) {
		int i = it1->first;
		double w_i = it1->second.first;
		if (w_i < min_nw) {
			min_nw = w_i;
		}
		if (w_i > max_nw) {
			max_nw = w_i;
		}
		for (auto it2 = it1->second.second.begin(); it2 != it1->second.second.end(); it2++) {
			int j = it2->first;
			if (i < j) { // edge (i,j)
				double c_ij = it2->second;
				if (c_ij < min_ec) {
					min_ec = c_ij;
				}
				if (c_ij > max_ec) {
					max_ec = c_ij;
				}
			}
		}
	}


	for (auto it1 = input_graph.begin(); it1 != input_graph.end(); it1++) {
		int i = it1->first;
		double w_i = it1->second.first;
		double new_w_i = (w_i - min_nw) / (max_nw - min_nw);
		it1->second.first = new_w_i;
		for (auto it2 = it1->second.second.begin(); it2 != it1->second.second.end(); it2++) {
			int j = it2->first;
			if (i < j) { // edge (i,j)
				double c_ij = it2->second;
				double new_c_ij = (c_ij - min_ec) / (max_ec - min_ec);
				it2->second = new_c_ij;
				input_graph[j].second[i] = new_c_ij;
			}
		}
	}

}