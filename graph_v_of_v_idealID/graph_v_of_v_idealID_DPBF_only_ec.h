#pragma once

/*this is the DPBF algorithm in Ding, Bolin, et al. "Finding top-k min-cost connected trees in databases." 2007 IEEE 23rd International Conference on Data Engineering. IEEE, 2007.

time complexity: O( 4^|Gamma| + 3^|Gamma||V|+ 2^|Gamma|* (|E| + |V|*(|Gamma| + log V)) )*/

#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <boost/heap/fibonacci_heap.hpp> 
#include <graph_v_of_v_idealID/graph_v_of_v_idealID_PrunedDPPlusPlus.h>

graph_hash_of_mixed_weighted graph_v_of_v_idealID_DPBF_only_ec(
	graph_v_of_v_idealID& input_graph, graph_v_of_v_idealID& group_graph, std::unordered_set<int>& cumpulsory_group_vertices, double& RAM_MB) {

	/*time complexity: O( 4^|Gamma| + 3^|Gamma||V|+ 2^|Gamma|* (|E| + |V|*(|Gamma| + log V)) ) */
	double bit_num = 0;

	if (cumpulsory_group_vertices.size() >= 20) {
		std::cout << "cumpulsory_group_vertices.size() is too large for graph_hash_of_mixed_weighted_DPBF_edge_weighted!" << std::endl;
		exit(1);
	}

	int N = input_graph.size();

	/*initialization; time complexity: O(4^|Gamma|)*/
	double inf = std::numeric_limits<double>::max(); // cost of empty tree is inf
	int group_sets_ID_range = pow(2, cumpulsory_group_vertices.size()) - 1; // the number of group sets: 2^|Gamma|, including empty set;   |Gamma| should be smaller than 31 due to precision
	/*Group		G1	G0	group_set_ID
				0   0   0
				0	1	1
				1	0	2
				1	1	3*/
	vector<vector<int>> non_overlapped_group_sets_IDs = graph_v_of_v_idealID_PrunedDPPlusPlus_non_overlapped_group_sets(group_sets_ID_range);  // time complexity: O(4^|Gamma|)
	for (auto it = non_overlapped_group_sets_IDs.begin(); it != non_overlapped_group_sets_IDs.end(); it++) {
		bit_num += sizeof(vector<int>) + (*it).size() * 4;
	}

	boost::heap::fibonacci_heap<graph_v_of_v_idealID_PrunedDPPlusPlus_min_node> Q_T; // min queues of trees
	std::unordered_map<string, handle_graph_v_of_v_idealID_PrunedDPPlusPlus_min_node> Q_T_handles; // key is string "v_p" ("vertex_ID" + "_" + "group set ID"); Q_T_handles only contain keys that are in Q_T




	/*initialize trees with vertices;

	time complexity: O(2^|Gamma|*|V|)

	every vertex v is associated with T(v,p) only when v covers p, otherwise the cost of T(v,p) is considered as inf;

	every vertex v is associated with at most 2^|Gamma| trees;
	*/
	vector<std::unordered_map<int, graph_v_of_v_idealID_PrunedDPPlusPlus_tree_node>> trees(N); // <v, <p, T(v,p)>>; this cannot be changed to vector of vectors
	int xxm = sizeof(graph_v_of_v_idealID_PrunedDPPlusPlus_min_node) + sizeof(graph_v_of_v_idealID_PrunedDPPlusPlus_tree_node);
	for (int v = 0; v < N; v++) {
		int group_set_ID_v = graph_v_of_v_idealID_PrunedDPPlusPlus_vertex_group_set_ID(v, group_graph, cumpulsory_group_vertices); /*time complexity: O(|Gamma|)*/
		for (int p = 1; p <= group_set_ID_v; p++) { // p is non-empty; time complexity: O(2^|Gamma|)
			if ((p | group_set_ID_v) == group_set_ID_v) { // p represents a non-empty group set inside group_set_ID_v, including group_set_ID_v

				/*T(v,p)*/
				graph_v_of_v_idealID_PrunedDPPlusPlus_tree_node node;
				node.cost = 0;
				node.type = 0;
				trees[v][p] = node;

				/*insert T(v,p) into Q_T*/
				graph_v_of_v_idealID_PrunedDPPlusPlus_min_node x;
				x.v = v;
				x.p = p;
				x.priority_value = 0;
				string handle_ID = to_string(v) + "_" + to_string(p);
				Q_T_handles[handle_ID] = Q_T.push(x);

				bit_num += xxm + 2 * 8;
			}
		}
	}


	double Q_T_max_size = 0;
	/*Big while loop*/
	while (Q_T.size() > 0) { // at most 2^|Gamma|*V loops

		double Q_T_size = Q_T.size();
		if (Q_T_size > Q_T_max_size) {
			Q_T_max_size = Q_T_size;
		}

		graph_v_of_v_idealID_PrunedDPPlusPlus_min_node top_node = Q_T.top();
		int v = top_node.v, p = top_node.p;

		Q_T.pop(); // O(2^|Gamma|*V*(|Gamma| + log V)) throught the loop, since Q_T contains at most 2^|Gamma|*V elements
		string handle_ID = to_string(v) + "_" + to_string(p);
		Q_T_handles.erase(handle_ID); // Q_T_handles only contains handles of elements in Q_T


		if (p == group_sets_ID_range) { // T(v,p) covers all groups 	

			bit_num += Q_T_max_size * (sizeof(graph_v_of_v_idealID_PrunedDPPlusPlus_min_node) + 8 + sizeof(handle_graph_v_of_v_idealID_PrunedDPPlusPlus_min_node) + 8); // Q_T + Q_T_handles
			RAM_MB = bit_num / 1024 / 1024;

			return graph_v_of_v_idealID_PrunedDPPlusPlus_build_tree(v, p, input_graph, trees); // time complexity: O(|V|)
		}


		/*grow*/
		for (auto it = input_graph[v].begin(); it != input_graph[v].end(); it++) {

			/*below: O(2^|Gamma|*E) in all loops, since each v has 2^|Gamma| times*/
			int u = it->first;
			double cost_euv = it->second;
			double grow_tree_cost = trees[v][p].cost + cost_euv;
			double T_up_cost;
			if (trees[u].count(p) == 0) {
				T_up_cost = inf;
			}
			else {
				T_up_cost = trees[u][p].cost;
			}


			if (grow_tree_cost < T_up_cost) {

				/*below: O(2^|Gamma|*V*(|Gamma| + log V)) throught the loop, since each u is checked 2^|Gamma| times, and Q_T contains at most 2^|Gamma|*V elements */

				/*update T(u,p) by grow T(v,p) with (u,v)*/
				trees[u][p].cost = grow_tree_cost;
				trees[u][p].type = 1;
				trees[u][p].u = v;

				/*update T(u,p) in Q_T*/
				graph_v_of_v_idealID_PrunedDPPlusPlus_min_node x;
				x.v = u;
				x.p = p;
				x.priority_value = grow_tree_cost;
				handle_ID = to_string(u) + "_" + to_string(p);
				if (Q_T_handles.count(handle_ID) == 0) { // T(u,p) is not in Q_T
					Q_T_handles[handle_ID] = Q_T.push(x);
				}
				else { // T(u,p) is in Q_T
					Q_T.update(Q_T_handles[handle_ID], x); // O(1) for decrease key
				}
			}
		}


		/*merge*/
		int p1 = p;
		for (auto it = non_overlapped_group_sets_IDs[p1].begin(); it != non_overlapped_group_sets_IDs[p1].end(); it++) {
			int p2 = *it; // p2 is not overlapped with p1
			double cost_Tvp1;
			if (trees[v].count(p1) == 0) {
				cost_Tvp1 = inf;
			}
			else {
				cost_Tvp1 = trees[v][p1].cost;
			}
			double cost_Tvp2;
			if (trees[v].count(p2) == 0) {
				cost_Tvp2 = inf;
			}
			else {
				cost_Tvp2 = trees[v][p2].cost;
			}
			int p1_cup_p2 = p1 + p2;
			double cost_Tvp1_cup_p2;
			if (trees[v].count(p1_cup_p2) == 0) {
				cost_Tvp1_cup_p2 = inf;
			}
			else {
				cost_Tvp1_cup_p2 = trees[v][p1_cup_p2].cost;
			}
			double merged_tree_cost = cost_Tvp1 + cost_Tvp2;
			if (merged_tree_cost < cost_Tvp1_cup_p2) { // O(3^|Gamma||V| comparisons in totel, see the DPBF paper)

				/*update T(v,p1_cup_p2) by merge T(v,p1) with T(v,v2)*/
				trees[v][p1_cup_p2].cost = merged_tree_cost;
				trees[v][p1_cup_p2].type = 2;
				trees[v][p1_cup_p2].p1 = p1;
				trees[v][p1_cup_p2].p2 = p2;

				/*update T(v,p1_cup_p2) in Q_T*/
				graph_v_of_v_idealID_PrunedDPPlusPlus_min_node x;
				x.v = v;
				x.p = p1_cup_p2;
				x.priority_value = merged_tree_cost;
				handle_ID = to_string(v) + "_" + to_string(p1_cup_p2);
				if (Q_T_handles.count(handle_ID) == 0) { // T(v,p1_cup_p2) is not in Q_T
					Q_T_handles[handle_ID] = Q_T.push(x);
				}
				else { // T(v,p1_cup_p2) is in Q_T
					Q_T.update(Q_T_handles[handle_ID], x); // O(1) for decrease key
				}

			}
		}

	}

	std::cout << "graph_v_of_v_idealID_DPBF_only_ec did not find a feasible solution!" << std::endl;
	getchar();
	exit(1);

}


#pragma region
void test_graph_v_of_v_idealID_DPBF_only_ec() {

	/*parameters*/
	int iteration_times = 100;
	int V = 1000, E = 5000, G = 5, g_size_min = 10, g_size_max = 50, precision = 3;
	double ec_min = 0.001, ec_max = 1; // PrunedDP does not allow zero edge weight



	double solution_cost_DPBF_sum = 0, solution_cost_PrunedDPPlusPlus_sum = 0;

	double time_DPBF_avg = 0, time_PrunedDPPlusPlus_avg = 0;

	/*iteration*/
	for (int i = 0; i < iteration_times; i++) {

		/*input and output*/
		int generate_new_graph = 1;
		double lambda = 1;
		std::unordered_set<int> generated_group_vertices;
		graph_hash_of_mixed_weighted instance_graph, generated_group_graph;
		if (generate_new_graph == 1) {
			instance_graph = graph_hash_of_mixed_weighted_generate_random_connected_graph(V, E, 0, 0, ec_min, ec_max, precision);
			graph_hash_of_mixed_weighted_generate_random_groups_of_vertices(G, g_size_min, g_size_max,
				instance_graph, instance_graph.hash_of_vectors.size(), generated_group_vertices, generated_group_graph);
			graph_hash_of_mixed_weighted_save_for_GSTP("simple_iterative_tests.text", instance_graph,
				generated_group_graph, generated_group_vertices, lambda);
		}
		else {
			graph_hash_of_mixed_weighted_read_for_GSTP("simple_iterative_tests.text", instance_graph,
				generated_group_graph, generated_group_vertices, lambda);
		}

		/*graph_hash_of_mixed_weighted_DPBF_vertex_edge_weighted*/
		if (1) {
			auto begin = std::chrono::high_resolution_clock::now();
			graph_hash_of_mixed_weighted solu = graph_hash_of_mixed_weighted_DPBF_edge_weighted(instance_graph, generated_group_graph, generated_group_vertices);
			auto end = std::chrono::high_resolution_clock::now();
			double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
			time_DPBF_avg = time_DPBF_avg + (double)runningtime / iteration_times;

			double cost = graph_hash_of_mixed_weighted_sum_of_ec(solu);

			//graph_hash_of_mixed_weighted_print(solu);
			//cout << "cost=" << cost << endl;

			solution_cost_DPBF_sum = solution_cost_DPBF_sum + cost;

			if (!this_is_a_feasible_solution(solu, generated_group_graph, generated_group_vertices)) {
				cout << "Error: graph_hash_of_mixed_weighted_DPBF_vertex_edge_weighted is not feasible!" << endl;
				graph_hash_of_mixed_weighted_print(solu);
				exit(1);
			}
		}

		unordered_map<int, int> vertexID_old_to_new;
		for (int mm = 0; mm < V; mm++) {
			vertexID_old_to_new[mm] = mm;
		}
		graph_v_of_v_idealID v_instance_graph = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(instance_graph, vertexID_old_to_new);
		for (int mm = V; mm < V + G; mm++) {
			vertexID_old_to_new[mm] = mm;
		}
		graph_v_of_v_idealID v_generated_group_graph = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(generated_group_graph, vertexID_old_to_new);

		/*graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted*/
		if (1) {
			double RAM;
			auto begin = std::chrono::high_resolution_clock::now();
			graph_hash_of_mixed_weighted solu = graph_v_of_v_idealID_DPBF_only_ec(v_instance_graph, v_generated_group_graph, generated_group_vertices, RAM);
			auto end = std::chrono::high_resolution_clock::now();
			double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
			time_PrunedDPPlusPlus_avg = time_PrunedDPPlusPlus_avg + (double)runningtime / iteration_times;

			//graph_hash_of_mixed_weighted_print(solu);

			double cost = graph_hash_of_mixed_weighted_sum_of_ec(solu);
			solution_cost_PrunedDPPlusPlus_sum = solution_cost_PrunedDPPlusPlus_sum + cost;

			if (!this_is_a_feasible_solution(solu, generated_group_graph, generated_group_vertices)) {
				cout << "Error: graph_v_of_v_idealID_DPBF_only_ec is not feasible!" << endl;
				graph_hash_of_mixed_weighted_print(solu);
				exit(1);
			}
		}

		if (solution_cost_DPBF_sum + 1e-8 < solution_cost_PrunedDPPlusPlus_sum) {
			cout << "solution_cost_DPBF_sum=" << solution_cost_DPBF_sum << endl;
			cout << "solution_cost_PrunedDPPlusPlus_sum=" << solution_cost_PrunedDPPlusPlus_sum << endl;
			cout << "wrong answer" << endl;
			getchar();
		}


	}

	cout << "solution_cost_DPBF_sum=" << solution_cost_DPBF_sum << endl;
	cout << "solution_cost_PrunedDPPlusPlus_sum=" << solution_cost_PrunedDPPlusPlus_sum << endl;

	cout << "time_DPBF_avg=" << time_DPBF_avg << "s" << endl;
	cout << "time_PrunedDPPlusPlus_avg=" << time_PrunedDPPlusPlus_avg << "s" << endl;
}
#pragma endregion test_graph_v_of_v_idealID_DPBF_only_ec