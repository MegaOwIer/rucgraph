#pragma once

/* in the generated graph, 0 to V-1 are vertices, and V to V+time_slot_num-1 are time slots */

#include <boost/random.hpp>
graph_hash_of_mixed_weighted graph_hash_of_mixed_weighted_generate_random_query_state_graph(int V, int time_slot_num, int query_state_graph_E, boost::random::mt19937& boost_random_time_seed) {

	/*time complexity: O(time_slot_num*query_state_graph_E)*/

	/*time complexity: O(|V|)*/
	graph_hash_of_mixed_weighted random_graph; // generate vertices
	for (int i = 0; i < V; i++) {
		graph_hash_of_mixed_weighted_add_vertex(random_graph, i, 0);
	}

	/*add edges to random_graph*/
	if (query_state_graph_E == V * time_slot_num) { // complete graphs
		/*time complexity: O(|V|+|E|)*/
		for (int i = 0; i < V; i++) {
			for (int j = 0; j < time_slot_num; j++) {
				graph_hash_of_mixed_weighted_add_edge(random_graph, i, V + j, 0);
			}
		}
	}
	else if (query_state_graph_E > V* time_slot_num) {
		std::cout <<
			"query_state_graph_E > V * time_slot_num in graph_hash_of_mixed_weighted_generate_random_query_state_graph!"
			<< '\n';
		exit(1);
	}
	else { // incomplete graphs

		/*time complexity: O(|V|)*/
		std::vector<int> not_full_vertices; // vertices without a full degree
		for (int i = 0; i < V; i++) {
			not_full_vertices.insert(not_full_vertices.end(), i);
		}


		/*time complexity: O(time_slot_num*query_state_graph_E)*/
		int edge_num = 0;
		while (edge_num < query_state_graph_E) {
			boost::random::uniform_int_distribution<> dist_id
			{ static_cast<int>(0), static_cast<int>(not_full_vertices.size() - 1) };
			int RAND = dist_id(boost_random_time_seed); // generate int random number  0, not_full_vertices.size()-1
			if (graph_hash_of_mixed_weighted_degree(random_graph, not_full_vertices[RAND]) < time_slot_num) { // this is a vertex without a full degree


				/*time complexity: O(time_slot_num)*/
				std::vector<int> unchecked(time_slot_num);
				std::iota(std::begin(unchecked), std::end(unchecked), 0);
				bool added = false;
				while (added == false) {
					boost::random::uniform_int_distribution<> dist_id2
					{ static_cast<int>(0), static_cast<int>(unchecked.size() - 1) };
					int x = dist_id2(boost_random_time_seed);
					int j = unchecked[x] + V;
					if (graph_hash_of_mixed_weighted_contain_edge(random_graph, not_full_vertices[RAND], j) == 0) {
						// This edge does not exist
						graph_hash_of_mixed_weighted_add_edge(random_graph, not_full_vertices[RAND], j, 0); // add a new edge
						edge_num++;
						added = true;
						break; // break after adding one edge
					}
					else {
						unchecked.erase(unchecked.begin() + x);
					}
				}




			}
			else { // this is a vertex with a full degree
				not_full_vertices.erase(not_full_vertices.begin() + RAND);
			}
		}


	}

	return random_graph;
}
