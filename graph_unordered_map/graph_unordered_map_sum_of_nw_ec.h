#pragma once


double graph_unordered_map_sum_of_nw_ec(graph_unordered_map& input_graph) {

	double sum = 0;

	for (auto it1 = input_graph.begin(); it1 != input_graph.end(); it1++) {
		int i = it1->first;
		double w_i = it1->second.first;
		sum = sum + w_i;
		for (auto it2 = it1->second.second.begin(); it2 != it1->second.second.end(); it2++) {
			int j = it2->first;
			if (i < j) { // edge (i,j)
				double c_ij = it2->second;
				sum = sum + c_ij;
			}
		}
	}

	return sum;

}




double graph_unordered_map_sum_of_nw(graph_unordered_map& input_graph) {

	double sum = 0;

	for (auto it1 = input_graph.begin(); it1 != input_graph.end(); it1++) {
		int i = it1->first;
		double w_i = it1->second.first;
		sum = sum + w_i;
	}

	return sum;

}



double graph_unordered_map_sum_of_ec(graph_unordered_map& input_graph) {

	double sum = 0;

	for (auto it1 = input_graph.begin(); it1 != input_graph.end(); it1++) {
		int i = it1->first;
		for (auto it2 = it1->second.second.begin(); it2 != it1->second.second.end(); it2++) {
			int j = it2->first;
			if (i < j) { // edge (i,j)
				double c_ij = it2->second;
				sum = sum + c_ij;
			}
		}
	}

	return sum;

}
