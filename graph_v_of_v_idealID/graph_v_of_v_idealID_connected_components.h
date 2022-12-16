#pragma once

#include <queue>


std::list<std::list<int>> graph_v_of_v_idealID_connected_components(graph_v_of_v_idealID& input_graph) {

	/*this is to find connected_components using breadth first search; time complexity O(|V|+|E|);
	related content: https://www.boost.org/doc/libs/1_68_0/boost/graph/connected_components.hpp
	https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm*/

	std::list<std::list<int>> components;

	/*time complexity: O(V)*/
	int N = input_graph.size();
	vector<bool> discovered(N, false);

	for (int i = 0; i < N; i++) {

		if (discovered[i] == false) {

			std::list<int> component;
			/*below is a deppth first search; time complexity O(|V|+|E|)*/
			std::queue<int> Q; // Queue is a data structure designed to operate in FIFO (First in First out) context.
			Q.push(i);
			component.push_back(i);
			discovered[i] = true;
			while (Q.size() > 0) {
				int v = Q.front();
				Q.pop(); //Removing that vertex from queue,whose neighbour will be visited now

				int adj_size = input_graph[v].size();
				for (int j = 0; j < adj_size; j++) {
					int adj_v = input_graph[v][j].first;
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

