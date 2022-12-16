#pragma once

#include <graph_unordered_map/graph_unordered_map.h> 
#include <queue>

std::list<std::list<int>> graph_unordered_map_connected_components(graph_unordered_map& input_graph) {

	/*this is to find connected_components using deppth first search; time complexity O(|V|+|E|);
	related content: https://www.boost.org/doc/libs/1_68_0/boost/graph/connected_components.hpp
	https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm*/

	std::list<std::list<int>> components;
	int N = input_graph.size();

	std::unordered_map<int, bool> discovered;
	for (auto it = input_graph.begin(); it != input_graph.end(); it++) {
		discovered[it->first] = false;
	}
	//std::vector<bool> discovered(N, false); // this is wrong, as the vertices may not be 1, ..., N-1

	for (auto it = input_graph.begin(); it != input_graph.end(); it++) {
		int i = it->first;
		if (discovered[i] == false) {

			std::list<int> component;
			/*below is a deppth first search*/
			std::queue<int> Q; // Queue is a data structure designed to operate in FIFO (First in First out) context.
			Q.push(i);
			component.push_back(i);
			discovered[i] = true;
			while (Q.size() > 0) {
				int v = Q.front();
				Q.pop(); //Removing that vertex from queue,whose neighbour will be visited now
				for (auto it = input_graph[v].second.begin(); it != input_graph[v].second.end(); it++) {
					int adj_v = it->first;
					if (discovered[adj_v] == false) {
						Q.push(adj_v);
						component.push_back(adj_v);
						discovered[adj_v] = true;
					}
				}
			}

			components.push_back(component);

		}
	}

	return components;

}