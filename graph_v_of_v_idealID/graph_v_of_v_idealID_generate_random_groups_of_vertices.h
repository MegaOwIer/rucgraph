#pragma once


#include <boost/random.hpp>

graph_v_of_v_idealID graph_v_of_v_idealID_generate_random_groups_of_vertices(int G, int g_size_min, int g_size_max, int input_graph_N, std::unordered_set<int>& generated_group_vertices,
	boost::random::mt19937& boost_random_time_seed) {

	/*time complexity: O(|G||V|)*/

	boost::random::uniform_int_distribution<> dist{ g_size_min, g_size_max };

	/*time complexity: O(|V|)*/
	graph_v_of_v_idealID group_graph(G + input_graph_N);


	/*add groups*/
	/*time complexity: O(|G||V|)*/
	int add_group_num = 0;
	while (add_group_num < G) {
		int group_size = dist(boost_random_time_seed); // generate int random number
		std::vector<int> to_be_linked_v(input_graph_N);
		std::iota(std::begin(to_be_linked_v), std::end(to_be_linked_v), 0);
		std::vector<int> linked_vertices;
		while (linked_vertices.size() < group_size) {
			boost::random::uniform_int_distribution<> dist2{ 0, (int)(to_be_linked_v.size() - 1) };
			int randID = dist2(boost_random_time_seed);
			linked_vertices.insert(linked_vertices.end(), to_be_linked_v[randID]);
			to_be_linked_v.erase(to_be_linked_v.begin() + randID);
		}

		// add this group
		int group_vertex = input_graph_N + add_group_num;
		generated_group_vertices.insert(group_vertex);
		for (int j = 0; j < linked_vertices.size(); j++) {
			graph_v_of_v_idealID_add_edge(group_graph, linked_vertices[j], group_vertex, 1);
		}

		add_group_num++;
	}

	return group_graph;

}
