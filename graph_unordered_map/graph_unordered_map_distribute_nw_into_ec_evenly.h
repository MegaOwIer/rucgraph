#pragma once

#include <graph_unordered_map.h>

void graph_unordered_map_distribute_nw_into_ec_evenly(graph_unordered_map& input_graph) {

	/*time complexity: O(|V|+|E|)*/

	for (auto it1 = input_graph.begin(); it1 != input_graph.end(); it1++) {
		int i = it1->first;
		double w_i = it1->second.first;
		double degree = it1->second.second.size();
		double divided_w_i = w_i / degree;
		for (auto it2 = it1->second.second.begin(); it2 != it1->second.second.end(); it2++) {
			int j = it2->first;
			it2->second = it2->second + divided_w_i; // evenly distribute nw into ec
		}
	}

}


void graph_unordered_map_distribute_nw_into_ec_evenly_reverse(graph_unordered_map& input_graph) {

	/*time complexity: O(|V|+|E|)*/

	for (auto it1 = input_graph.begin(); it1 != input_graph.end(); it1++) {
		int i = it1->first;
		double w_i = it1->second.first;
		double degree = it1->second.second.size();
		double divided_w_i = w_i / degree;
		for (auto it2 = it1->second.second.begin(); it2 != it1->second.second.end(); it2++) {
			int j = it2->first;
			it2->second = it2->second - divided_w_i; // remove the evenly distribute nw from ec
		}
	}

}
