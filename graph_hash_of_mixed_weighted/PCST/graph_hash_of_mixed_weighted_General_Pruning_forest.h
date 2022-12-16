#pragma once


/*  General_pruning: this function is to optimal prune a forest, and returns the maximum weight sub-tree;

	details in Sun, Yahui, et al. "The fast heuristic algorithms and post-processing techniques to design
	large and low-cost communication networks." IEEE/ACM Transactions on Networking 27.1 (2019): 375-388.

	time complexity: O(V)
*/


#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_breadth_first_search_a_set_of_vertices.h>

graph_hash_of_mixed_weighted graph_hash_of_mixed_weighted_General_Pruning_forest(graph_hash_of_mixed_weighted& input_trees) {

	/*the input is a tree or a forest, or empty*/
	if (input_trees.hash_of_vectors.size() == 0) {
		graph_hash_of_mixed_weighted solution_tree;
		return solution_tree;
	}

	/*initialization; time complexity: O(V)*/
	unordered_map<int, double> nw;
	unordered_map<int, bool> unprocessed;
	unordered_map<int, int> processing_degree;
	vector<int> target_vertex; // vertices which are unprocessed and have a processing_degree of 1
	for (auto it = input_trees.hash_of_vectors.begin(); it != input_trees.hash_of_vectors.end(); it++) {
		int v = it->first;
		nw[v] = it->second.vertex_weight; // initial nw values are node weights
		unprocessed[v] = true; // all the vertices in trees are unprocessed initially
		processing_degree[v] = input_trees.degree(v); /*initialize processing_degree as the degree*/
		if (processing_degree[v] == 1) {
			target_vertex.push_back(v); // pruning should start from these leaves
		}
	}

	/*update nw values without root; time complexity: O(V)*/
	while (target_vertex.size() > 0) { // end the process until there is no target vertex any more

		/*the last vetex (or vertices for forests) poped out of target_vertex will be the last unprocessed vertex
		(or vertices for forests); the following code will do nothing on this vertex*/

		int v = target_vertex[0]; // processing target_vertex[0]

		auto search = input_trees.hash_of_hashs.find(v);
		if (search != input_trees.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int adj_v = it2->first;
				if (unprocessed[adj_v] == true) { // adj_v is unprocessed, so adj_v is v_adj
					double ec = it2->second; // do not remove edges of v in this for loop, otherwise it2 points to outside
					if (ec < nw[v]) {
						nw[adj_v] = nw[adj_v] + nw[v] - ec; // update nw[adj_v]
					}
					unprocessed[v] = false; // mark v as processed
					processing_degree[adj_v]--; // update processing_degree[adj_v]
					if (processing_degree[adj_v] == 1) { // adj_v becomes a new target_vertex
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
				if (unprocessed[adj_v] == true) { // adj_v is unprocessed, so adj_v is v_adj
					double ec = it2->second; // do not remove edges of v in this for loop, otherwise it2 points to outside
					if (ec < nw[v]) {
						nw[adj_v] = nw[adj_v] + nw[v] - ec; // update nw[adj_v]
					}
					unprocessed[v] = false; // mark v as processed
					processing_degree[adj_v]--; // update processing_degree[adj_v]
					if (processing_degree[adj_v] == 1) { // adj_v becomes a new target_vertex
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
		if (root_nw < nw[v]) {
			root_nw = nw[v];
			root = v;
		}
	}
	//cout << "root = " << root << " root_nw = " << root_nw << endl;

	/* breadth_first_search all vertices connected to root; time complexity: O(V) */
	unordered_set<int> root_component = graph_hash_of_mixed_weighted_breadth_first_search_a_set_of_vertices(input_trees, root); // v is connected to root; including root

	/*initialization with root; time complexity: O(V)*/
	graph_hash_of_mixed_weighted solution_tree;
	for (auto it = root_component.begin(); it != root_component.end(); it++) {
		int v = *it;
		nw[v] = input_trees.hash_of_vectors[v].vertex_weight; // initial nw values are node weights
		unprocessed[v] = true; // all the vertices in trees are unprocessed initially
		processing_degree[v] = input_trees.degree(v); /*initialize processing_degree as the degree*/
		if (processing_degree[v] == 1 && v != root) { // v is not root
			target_vertex.push_back(v); // pruning should start from these leaves
		}

		graph_hash_of_mixed_weighted_add_vertex(solution_tree, v, nw[v]); // insert vertex
		auto search = input_trees.hash_of_hashs.find(v);
		if (search != input_trees.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int vertex = it2->first;
				if (vertex < v) { // vertex has already been inserted
					double ec = it2->second; // do not remove edges of v in this for loop, otherwise it2 points to outside
					graph_hash_of_mixed_weighted_add_edge(solution_tree, v, vertex, ec); // insert edge
				}
			}
		}
		else {
			auto search2 = input_trees.hash_of_vectors.find(v);
			for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
				int vertex = it2->first;
				if (vertex < v) { // vertex has already been inserted
					double ec = it2->second; // do not remove edges of v in this for loop, otherwise it2 points to outside
					graph_hash_of_mixed_weighted_add_edge(solution_tree, v, vertex, ec); // insert edge
				}
			}
		}
	}

	/*pruning with root; time complexity: O(V)*/
	while (target_vertex.size() > 0) { // end the process until there is no target vertex any more

		/*all vertices in target_vertex (without root) should be processed ultimately*/

		int v = target_vertex[0]; // processing target_vertex[0]

		auto search = solution_tree.hash_of_hashs.find(v);
		if (search != solution_tree.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int adj_v = it2->first;
				if (unprocessed[adj_v] == true) { // adj_v is unprocessed, so adj_v is v_adj
					double ec = it2->second; // do not remove edges of v in this for loop, otherwise it2 points to outside
					if (ec < nw[v]) {
						nw[adj_v] = nw[adj_v] + nw[v] - ec; // update nw[adj_v]
					}
					else {
						graph_hash_of_mixed_weighted_remove_edge_but_not_isolated_vertices(solution_tree, v, adj_v); // remove edge
					}
					unprocessed[v] = false; // mark v as processed
					processing_degree[adj_v]--; // update processing_degree[adj_v]
					if (processing_degree[adj_v] == 1 && adj_v != root) { // adj_v becomes a new target_vertex
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
				if (unprocessed[adj_v] == true) { // adj_v is unprocessed, so adj_v is v_adj
					double ec = it2->second; // do not remove edges of v in this for loop, otherwise it2 points to outside
					if (ec < nw[v]) {
						nw[adj_v] = nw[adj_v] + nw[v] - ec; // update nw[adj_v]
					}
					else {
						graph_hash_of_mixed_weighted_remove_edge_but_not_isolated_vertices(solution_tree, v, adj_v); // remove edge
					}
					unprocessed[v] = false; // mark v as processed
					processing_degree[adj_v]--; // update processing_degree[adj_v]
					if (processing_degree[adj_v] == 1 && adj_v != root) { // adj_v becomes a new target_vertex
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
	for (auto it = root_component.begin(); it != root_component.end(); it++) {
		int v = *it;
		if (new_root_component.count(v) == 0) {
			graph_hash_of_mixed_weighted_remove_vertex(solution_tree, v);
		}
	}

	return solution_tree;
}

