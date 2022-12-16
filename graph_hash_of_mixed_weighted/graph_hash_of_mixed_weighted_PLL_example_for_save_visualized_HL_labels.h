#pragma once


#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_PLL_single_thread.h>


void save_visualized_HL_labels(string save_name, vector<vector<PLL_sorted_label>>& L) {

	std::ofstream outputFile;
	outputFile.precision(5);
	outputFile.setf(std::ios::fixed);
	outputFile.setf(std::ios::showpoint);
	outputFile.open(save_name);

	outputFile << "visualize_HL_labels <vertex,parent_vertex,distance>:" << endl;
	for (int i = 0; i < L.size(); i++) {
		outputFile << "ID " << i << ": ";
		for (auto it = L[i].begin(); it != L[i].end(); it++) {
			outputFile << "<" << it->vertex << "," << it->parent_vertex << "," << it->distance << "> ";
		}
		outputFile << endl;
	}

}

#include <graph_v_of_v_idealID/graph_v_of_v_idealID_save_for_GSTP.h>

void PLL_example_for_save_visualized_HL_labels(graph_hash_of_mixed_weighted& input_graph, int max_N) {

	/*this function is different from (and faster than) PLL_hash_of_vectors_generate_indexes_basic
	in that this function change vertexIDs to ideal IDs,
	and then use the rank and T P tricks in 2013 Japan paper*/

	int N = input_graph.hash_of_vectors.size();

	/*sort vertices by degrees*/
	vector<pair<int, int>> sorted_vertices;
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		sorted_vertices.push_back({ it->first, input_graph.degree(it->first) });
	}
	sort(sorted_vertices.begin(), sorted_vertices.end(), PLL_Creating_Indexes_compare);


	/*graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID*/
	unordered_map<int, int> vertexID_old_to_new;
	vector<int> vertexID_new_to_old(N);
	for (int i = 0; i < N; i++) {
		vertexID_old_to_new[sorted_vertices[i].first] = i;
		vertexID_new_to_old[i] = sorted_vertices[i].first;
	}
	graph_v_of_v_idealID ideal_graph = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(input_graph, vertexID_old_to_new);

	graph_v_of_v_idealID group_graph;
	std::unordered_set<int> group_vertices;
	graph_v_of_v_idealID_save_for_GSTP("example_visualize_PLL_Labels_ideal_graph.txt", ideal_graph, group_graph, group_vertices);


	vector<vector<PLL_sorted_label>> L(N); // Labels for ideal IDs
	vector<double> P(N, std::numeric_limits<double>::max()); // P is distance to v_k, and priorities in Q
	vector<double> T(N, std::numeric_limits<double>::max());
	vector<PLL_handle_t_for_sp> Q_handles(N);

	/*PLL: Pruned Landmark Labeling*/
	for (int v_k = 0; v_k < N; v_k++) {

		/*Pruned Dijkstra from vertex v_k; see Algorithm 1 in 2013 Japan SIGMOD paper*/

		queue<int> P_changed_vertices;

		PLL_node_for_sp node;
		boost::heap::fibonacci_heap<PLL_node_for_sp> Q;
		PLL_sorted_label xx;

		node.vertex = v_k;
		node.parent_vertex = v_k;
		node.priority_value = 0;
		Q_handles[v_k] = Q.push(node);
		P[v_k] = 0;
		P_changed_vertices.push(v_k);


		int L_v_k_size = L[v_k].size();
		for (int i = 0; i < L_v_k_size; i++) {
			T[L[v_k][i].vertex] = L[v_k][i].distance; //allocate T values for L[v_k]
		}


		while (Q.size() > 0) {

			node = Q.top();
			Q.pop();
			int u = node.vertex;
			int u_parent = node.parent_vertex;
			double P_u = node.priority_value;
			double P_u_with_error = P_u + 1e-5; // 1e-5 is error

			/*merge-joinly compute query_v_k_u, which is the QUERY value in Line 7 of Algo. 1 in 2013 Japan paper*/
			double query_v_k_u = std::numeric_limits<double>::max(); // if disconnected, retun this large value
			auto L_u_size = L[u].size(); // a vector<PLL_sorted_label>
			for (int i = 0; i < L_u_size; i++) {
				double dis = L[u][i].distance + T[L[u][i].vertex];
				if (query_v_k_u > dis) {
					query_v_k_u = dis;
				}
			}

			if (P_u_with_error < query_v_k_u) { // this is pruning

				xx.vertex = v_k;
				xx.distance = P_u;
				xx.parent_vertex = u_parent;
				L[u].push_back(xx);

				int u_adj_size = ideal_graph[u].size();
				for (int i = 0; i < u_adj_size; i++) {
					int adj_v = ideal_graph[u][i].first;
					double ec = ideal_graph[u][i].second;
					if (P[adj_v] == std::numeric_limits<double>::max()) {
						node.vertex = adj_v;
						node.parent_vertex = u;
						node.priority_value = P_u + ec;
						Q_handles[adj_v] = Q.push(node);
						P[adj_v] = node.priority_value;
						P_changed_vertices.push(adj_v);
					}
					else {
						if (P[adj_v] > P_u + ec) {
							node.vertex = adj_v;
							node.parent_vertex = u;
							node.priority_value = P_u + ec;
							Q.update(Q_handles[adj_v], node);
							P[adj_v] = node.priority_value;
						}
					}
				}
			}
		}


		for (int i = 0; i < L_v_k_size; i++) {
			T[L[v_k][i].vertex] = std::numeric_limits<double>::max(); // reverse-allocate T values for L[v_k]
		}

		while (P_changed_vertices.size() > 0) {
			P[P_changed_vertices.front()] = std::numeric_limits<double>::max(); // reverse-allocate P values
			P_changed_vertices.pop();
		}

	}


	save_visualized_HL_labels("example_visualize_PLL_Labels_labels.txt", L);

}

void example_visualize_PLL_Labels() {

	/*parameters*/
	int V = 5, E = 8, precision = 1;
	double ec_min = 1, ec_max = 10; // set ec_min=ec_max=1 for testing unweighted PLL

	std::unordered_set<int> generated_group_vertices;
	graph_hash_of_mixed_weighted instance_graph, generated_group_graph;
	instance_graph = graph_hash_of_mixed_weighted_generate_random_graph(V, E, 0, 0, ec_min, ec_max, precision, boost_random_time_seed);

	PLL_example_for_save_visualized_HL_labels(instance_graph, V + 1);
}



/*
an example C++ code:

#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
using namespace std;

// header files in the Boost library: https://www.boost.org/
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };

#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_PLL_example_for_save_visualized_HL_labels.h>

int main()
{
	example_visualize_PLL_Labels();
}


check the saved two files
*/