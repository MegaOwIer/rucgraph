#pragma once


#include <graph_unordered_map/graph_unordered_map.h>


void graph_unordered_map_ec_update_pairwise_jaccard_distance(graph_unordered_map& input_graph) {

	/*time complexity: O(|V|+|E|)*/

	/*for introducing this distance: https://en.wikipedia.org/wiki/Jaccard_index ;
	The Jaccard distance, which measures dissimilarity between sample sets, is complementary
	to the Jaccard coefficient and is obtained by subtracting the Jaccard coefficient from 1;
	dist(i,j) = 1 - |V_i \cap V_j| / |V_i \cup V_j|, where V_i is the set of adjacent vertices of i;
	this is different from that in T Lappse's 2009 paper */

	/*time complexity: O(|V|+|E|)*/
	for (auto it1 = input_graph.begin(); it1 != input_graph.end(); it1++) {
		int i = it1->first;
		for (auto it2 = it1->second.second.begin(); it2 != it1->second.second.end(); it2++) {
			int j = it2->first;

			if (i < j) { // edge (i,j)
				/*update_pairwise_Jaccard_distance for edge (i,j)*/
				int V_i_cap_V_j = 0;
				for (auto it3 = it1->second.second.begin(); it3 != it1->second.second.end(); it3++) {
					int k = it3->first; // k is adjacenct to i
					if (input_graph[j].second.count(k) > 0) { // k is also adjacenct to j
						V_i_cap_V_j++;
					}
				}

				int V_i_cup_V_j = input_graph[i].second.size() + input_graph[j].second.size() - V_i_cap_V_j;
				double ec = 1 - (double)V_i_cap_V_j / V_i_cup_V_j;

				graph_unordered_map_add_edge(input_graph, i, j, ec);// update_pairwise_jaccard_distance for edge (j,i)

			}
		}
	}

}