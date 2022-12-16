#pragma once


/*  General_pruning: this function is to optimally prune a tree, and returns the maximum weight sub-tree;

	details in Sun, Yahui, et al. "The fast heuristic algorithms and post-processing techniques to design
	large and low-cost communication networks." IEEE/ACM Transactions on Networking 27.1 (2019): 375-388.

	time complexity: O(V)
*/
#include <chrono>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_breadth_first_search_a_set_of_vertices.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_save_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_read_graph_with_weight.h>
#include <vector>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_breadth_first_search_a_set_of_vertices.h>

class General_Pruning_tree_node {
public:
	double nw;
	bool unprocessed;
	int processing_degree;
};

graph_hash_of_mixed_weighted graph_hash_of_mixed_weighted_General_Pruning_tree(graph_hash_of_mixed_weighted& input_trees) {

	/*the input is a tree (may be a single vertex with the largest weight), or empty*/
	if (input_trees.hash_of_vectors.size() == 0) {
		graph_hash_of_mixed_weighted solution_tree;
		return solution_tree;
	}

	/*initialization; time complexity: O(V)*/
	unordered_map<int, General_Pruning_tree_node> nodes;
	General_Pruning_tree_node node;
	vector<int> target_vertex; // vertices which are unprocessed and have a processing_degree of 1
	for (auto it = input_trees.hash_of_vectors.begin(); it != input_trees.hash_of_vectors.end(); it++) {
		int v = it->first;
		node.nw = it->second.vertex_weight; // initial nw values are node weights
		node.unprocessed = true; // all the vertices in trees are unprocessed initially
		node.processing_degree = input_trees.degree(v); /*initialize processing_degree as the degree*/
		nodes[v] = node;
		if (node.processing_degree == 1) {
			target_vertex.push_back(v); // pruning should start from these leaves
		}
	}

	/*update nw values without root; time complexity: O(V)*/
	while (target_vertex.size() > 0) { // end the process until there is no target vertex any more

		/*the last vetex (or vertices for forests) poped out of target_vertex will be the last unprocessed vertex
		(or vertices for forests); the following code will do nothing on this vertex*/

		int v = target_vertex[0]; // processing target_vertex[0]
		auto pointer_v = nodes.find(v);

		auto search = input_trees.hash_of_hashs.find(v);
		if (search != input_trees.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int adj_v = it2->first;
				auto pointer_adj_v = nodes.find(adj_v);
				if (pointer_adj_v->second.unprocessed == true) { // adj_v is unprocessed, so adj_v is v_adj
					double ec = it2->second; // do not remove edges of v in this for loop, otherwise it2 points to outside
					if (ec < pointer_v->second.nw) {
						pointer_adj_v->second.nw = pointer_adj_v->second.nw + pointer_v->second.nw - ec; // update nw[adj_v]
					}
					pointer_v->second.unprocessed = false; // mark v as processed
					pointer_adj_v->second.processing_degree--; // update processing_degree[adj_v]
					if (pointer_adj_v->second.processing_degree == 1) { // adj_v becomes a new target_vertex
						target_vertex.insert(target_vertex.end(), adj_v);
					}
					break; // there is at most one v_adj (finally, target_vertex[0] is the remaining unprocessed vertex)
				}
			}
		}
		else {
			auto search2 = input_trees.hash_of_vectors.find(v);
			for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
				int adj_v = it2->first;
				auto pointer_adj_v = nodes.find(adj_v);
				if (pointer_adj_v->second.unprocessed == true) { // adj_v is unprocessed, so adj_v is v_adj
					double ec = it2->second; // do not remove edges of v in this for loop, otherwise it2 points to outside
					if (ec < pointer_v->second.nw) {
						pointer_adj_v->second.nw = pointer_adj_v->second.nw + pointer_v->second.nw - ec; // update nw[adj_v]
					}
					pointer_v->second.unprocessed = false; // mark v as processed
					pointer_adj_v->second.processing_degree--; // update processing_degree[adj_v]
					if (pointer_adj_v->second.processing_degree == 1) { // adj_v becomes a new target_vertex
						target_vertex.insert(target_vertex.end(), adj_v);
					}
					break; // there is at most one v_adj (finally, target_vertex[0] is the remaining unprocessed vertex)
				}
			}
		}

		target_vertex.erase(target_vertex.begin()); // erase target_vertex[0]
	}

	/* find the vertex with the largest nw as the root */
	int root;
	double root_nw = -INT_MAX;
	for (auto it = input_trees.hash_of_vectors.begin(); it != input_trees.hash_of_vectors.end(); it++) {
		int v = it->first;
		double xx = nodes[v].nw;
		if (root_nw < xx) {
			root_nw = xx;
			root = v;
		}
	}
	//cout << "root = " << root << " root_nw = " << root_nw << endl;

	/*initialization with root; time complexity: O(V)*/
	graph_hash_of_mixed_weighted solution_tree = graph_hash_of_mixed_weighted_copy_graph(input_trees);
	for (auto it = input_trees.hash_of_vectors.begin(); it != input_trees.hash_of_vectors.end(); it++) {
		int v = it->first;
		auto pointer = nodes.find(v);
		pointer->second.nw = input_trees.hash_of_vectors[v].vertex_weight; // initial nw values are node weights
		pointer->second.unprocessed = true; // all the vertices in trees are unprocessed initially
		pointer->second.processing_degree = input_trees.degree(v); /*initialize processing_degree as the degree*/
		if (pointer->second.processing_degree == 1 && v != root) { // v is not root
			target_vertex.push_back(v); // pruning should start from these leaves
		}
	}

	/*pruning with root; time complexity: O(V)*/
	while (target_vertex.size() > 0) { // end the process until there is no target vertex any more

		/*all vertices in target_vertex (without root) should be processed ultimately*/

		int v = target_vertex[0]; // processing target_vertex[0]
		auto pointer_v = nodes.find(v);

		auto search = solution_tree.hash_of_hashs.find(v);
		if (search != solution_tree.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int adj_v = it2->first;
				auto pointer_adj_v = nodes.find(adj_v);
				if (pointer_adj_v->second.unprocessed == true) { // adj_v is unprocessed, so adj_v is v_adj
					double ec = it2->second; // do not remove edges of v in this for loop, otherwise it2 points to outside
					if (ec < pointer_v->second.nw) {
						pointer_adj_v->second.nw = pointer_adj_v->second.nw + pointer_v->second.nw - ec; // update nw[adj_v]
					}
					else {
						graph_hash_of_mixed_weighted_remove_edge_but_not_isolated_vertices(solution_tree, v, adj_v); // remove edge
					}
					pointer_v->second.unprocessed = false; // mark v as processed
					pointer_adj_v->second.processing_degree--; // update processing_degree[adj_v]
					if (pointer_adj_v->second.processing_degree == 1 && adj_v != root) { // adj_v becomes a new target_vertex
						target_vertex.insert(target_vertex.end(), adj_v);
					}
					break; // there is at most one v_adj (finally, target_vertex[0] is the remaining unprocessed vertex)
				}
			}
		}
		else {
			auto search2 = solution_tree.hash_of_vectors.find(v);
			for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
				int adj_v = it2->first;
				auto pointer_adj_v = nodes.find(adj_v);
				if (pointer_adj_v->second.unprocessed == true) { // adj_v is unprocessed, so adj_v is v_adj
					double ec = it2->second; // do not remove edges of v in this for loop, otherwise it2 points to outside
					if (ec < pointer_v->second.nw) {
						pointer_adj_v->second.nw = pointer_adj_v->second.nw + pointer_v->second.nw - ec; // update nw[adj_v]
					}
					else {
						graph_hash_of_mixed_weighted_remove_edge_but_not_isolated_vertices(solution_tree, v, adj_v); // remove edge
					}
					pointer_v->second.unprocessed = false; // mark v as processed
					pointer_adj_v->second.processing_degree--; // update processing_degree[adj_v]
					if (pointer_adj_v->second.processing_degree == 1 && adj_v != root) { // adj_v becomes a new target_vertex
						target_vertex.insert(target_vertex.end(), adj_v);
					}
					break; // there is at most one v_adj (finally, target_vertex[0] is the remaining unprocessed vertex)
				}
			}
		}

		target_vertex.erase(target_vertex.begin()); // erase target_vertex[0]
	}


	/* remove disconnected part from solution_tree */
	unordered_set<int> new_root_component = graph_hash_of_mixed_weighted_breadth_first_search_a_set_of_vertices(solution_tree, root); // v is connected to root; including root
	for (auto it = input_trees.hash_of_vectors.begin(); it != input_trees.hash_of_vectors.end(); it++) {
		int v = it->first;
		if (new_root_component.count(v) == 0) {
			graph_hash_of_mixed_weighted_remove_vertex(solution_tree, v);
		}
	}

	return solution_tree;
}



#pragma region
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_General_Pruning_tree_old.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_generate_random_connected_graph.h>
double PCST_weight1(graph_hash_of_mixed_weighted& solu_tree) {

	double total_weight = 0;

	for (auto it1 = solu_tree.hash_of_vectors.begin(); it1 != solu_tree.hash_of_vectors.end(); it1++) {
		int i = it1->first;
		double nw = it1->second.vertex_weight;
		total_weight = total_weight + nw;

		auto search = solu_tree.hash_of_hashs.find(i);
		if (search != solu_tree.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int j = it2->first;
				if (i < j) {
					double ec = it2->second;
					total_weight = total_weight - ec;
				}
			}
		}
		else {
			auto search2 = solu_tree.hash_of_vectors.find(i);
			for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
				int j = it2->first;
				if (i < j) {
					double ec = it2->second;
					total_weight = total_weight - ec;
				}
			}
		}
	}

	return total_weight;

}

void compare_General_Pruning_tree_old_new() {

	/*parameters*/
	int iteration_times = 100;
	int V = 1000, E = V - 1, precision = 2;
	double alpha, ec_min = 2, ec_max = 10, nw_min = 1, nw_max = 5;

	double solution_old_sum = 0, solution_new_sum = 0;
	double time_old_avg = 0, time_new_avg = 0;


	/*iteration*/
	for (int i = 0; i < iteration_times; i++) {

		/*input and output*/
		int generate_new_graph = 1;
		graph_hash_of_mixed_weighted instance_graph;
		if (generate_new_graph == 1) {
			instance_graph = graph_hash_of_mixed_weighted_generate_random_connected_graph(V, E, nw_min, nw_max, ec_min, ec_max, precision);
			graph_hash_of_mixed_weighted_save_graph_with_weight("instance_graph.txt", instance_graph, alpha);
		}
		else {
			graph_hash_of_mixed_weighted_read_graph_with_weight("instance_graph.txt", instance_graph, alpha);
		}

		if (1) {
			auto begin1 = std::chrono::high_resolution_clock::now();
			graph_hash_of_mixed_weighted solu = graph_hash_of_mixed_weighted_General_Pruning_tree_old(instance_graph);
			auto end1 = std::chrono::high_resolution_clock::now();
			double runningtime1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1).count() / 1e9; // s
			time_old_avg = time_old_avg + (double)runningtime1 / iteration_times;
			solution_old_sum = solution_old_sum + PCST_weight1(solu);
		}

		if (1) {
			auto begin1 = std::chrono::high_resolution_clock::now();
			graph_hash_of_mixed_weighted solu = graph_hash_of_mixed_weighted_General_Pruning_tree(instance_graph);
			auto end1 = std::chrono::high_resolution_clock::now();
			double runningtime1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1).count() / 1e9; // s
			time_new_avg = time_new_avg + (double)runningtime1 / iteration_times;
			solution_new_sum = solution_new_sum + PCST_weight1(solu);
		}

	}


	cout << "solution_old_sum = " << solution_old_sum << endl;
	cout << "solution_new_sum = " << solution_new_sum << endl;

	cout << "time_old_avg = " << time_old_avg << "s" << endl;
	cout << "time_new_avg = " << time_new_avg << "s" << endl;
}
#pragma endregion compare_General_Pruning_tree_old_new