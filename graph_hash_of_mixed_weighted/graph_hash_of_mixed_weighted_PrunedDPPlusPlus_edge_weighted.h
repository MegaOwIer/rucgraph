#pragma once

/*
time complexity: O(  )

cost function: (1 - lambda) * vertex_weights + lambda * edge_weights
*/
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <boost/heap/fibonacci_heap.hpp> 
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_shortest_paths.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_MST_postprocessing.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_sum_of_nw_ec.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_save_for_GSTP.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_read_for_GSTP.h>

#pragma region
pair<std::unordered_map<int, int>, std::unordered_map<int, double>> graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_find_SPs_to_g
(graph_hash_of_mixed_weighted& group_graph, graph_hash_of_mixed_weighted& input_graph, int g_vertex) {

	/*time complexity: O(|E|+|V|log|V|)*/

	/*add dummy vertex and edges; time complexity: O(|V|)*/
	graph_hash_of_mixed_weighted_add_vertex(input_graph, g_vertex, 0); // add dummy vertex
	auto search = group_graph.hash_of_hashs.find(g_vertex);
	if (search != group_graph.hash_of_hashs.end()) {
		for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
			int vertex = it2->first; // vertex is in group g_vertex
			graph_hash_of_mixed_weighted_add_edge(input_graph, g_vertex, vertex, 0); // add dummy edge
		}
	}
	else {
		auto search2 = group_graph.hash_of_vectors.find(g_vertex); // if v is not in g, error is here
		for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
			int vertex = it2->first; // vertex is in group g_vertex
			graph_hash_of_mixed_weighted_add_edge(input_graph, g_vertex, vertex, 0); // add dummy edge
		}
	}



	/*time complexity: O(|E|+|V|log|V|)*/
	std::unordered_map<int, double> distances; // total vertex and edge weights of paths
	std::unordered_map<int, int> predecessors;
	graph_hash_of_mixed_weighted_shortest_paths_source_to_all(input_graph, g_vertex, distances, predecessors);
	graph_hash_of_mixed_weighted_remove_vertex(input_graph, g_vertex);  // all dummy vertex and edges are removed; time complexity: O(|V|)
	distances.erase(g_vertex);
	predecessors.erase(g_vertex);
	for (auto it = predecessors.begin(); it != predecessors.end(); it++) {
		int pre = it->second;
		if (pre == g_vertex) {
			it->second = it->first; // since g_vertex is not in predecessors, it->second points to it->first, i.e., the path ends at it->first.
		}
	}

	return { predecessors, distances };

}

std::unordered_map<int, pair<std::unordered_map<int, int>, std::unordered_map<int, double>>> graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_find_SPs
(graph_hash_of_mixed_weighted& input_graph, graph_hash_of_mixed_weighted& group_graph, std::unordered_set<int>& cumpulsory_group_vertices) {


	/*time complexity: O(|T||E|+|T||V|log|V|)*/
	std::unordered_map<int, pair<std::unordered_map<int, int>, std::unordered_map<int, double>>> SPs_to_groups;
	for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++) {
		int g_vertex = *it;
		SPs_to_groups[g_vertex] = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_find_SPs_to_g(group_graph, input_graph, g_vertex);
	}

	return SPs_to_groups;
}

#pragma endregion graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_find_SPs

#pragma region
std::unordered_map<int, std::unordered_set<int>> graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_non_overlapped_group_sets(int group_sets_ID_range) {

	/*this function calculate the non-empty and non_overlapped_group_sets_IDs of each non-empty group_set ID;

	time complexity: O(4^|Gamma|), since group_sets_ID_range=2^|Gamma|;

	the original DPBF code use the same method in this function, and thus has the same O(4^|Gamma|) complexity;*/

	std::unordered_map<int, std::unordered_set<int>> non_overlapped_group_sets_IDs; // <set_ID, non_overlapped_group_sets_IDs>

	for (int i = 1; i <= group_sets_ID_range; i++) { // i is a nonempty group_set ID
		non_overlapped_group_sets_IDs[i] = {};
		for (int j = 1; j < group_sets_ID_range; j++) { // j is another nonempty group_set ID
			if ((i & j) == 0) { // i and j are non-overlapping group sets
				/* The & (bitwise AND) in C or C++ takes two numbers as operands and does AND on every bit of two numbers. The result of AND for each bit is 1 only if both bits are 1.
				https://www.programiz.com/cpp-programming/bitwise-operators */
				non_overlapped_group_sets_IDs[i].insert(j);
			}
		}
	}

	return non_overlapped_group_sets_IDs;

}
#pragma endregion graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_non_overlapped_group_sets

#pragma region
struct graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node {
	int v;
	int p; // group_set_ID
	//double tree_cost; // cost of T(v,p)
	double priority_value; // this is lb in PrunedDP++, but is tree_cost in others
};
bool operator<(graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node const& x, graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node const& y) {
	return x.priority_value > y.priority_value; // < is the max-heap; > is the mean heap; PriorityQueue is expected to be a max-heap of integer values
}
typedef typename boost::heap::fibonacci_heap<graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node>::handle_type handle_graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node;
#pragma endregion graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted priority queue

#pragma region
class graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_tree_node
{
	/*this is like the tree T(v,p) in the DPBF paper*/

public:

	int type; // =0: this is the single vertex v; =1: this tree is built by grown; =2: built by merge

	double cost; // cost of this tree T(v,p);

	int u; // if this tree is built by grown, then it's built by growing edge (v,u);

	int p1, p2; // if this tree is built by merge, then it's built by merge T(v,p1) and T(v,p2);

};
#pragma endregion graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_tree_node

#pragma region
int graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_vertex_group_set_ID(int vertex, graph_hash_of_mixed_weighted& group_graph,
	std::unordered_set<int>& cumpulsory_group_vertices) {

	/*time complexity: O(|Gamma|); this function returns the maximum group set ID for a single vertex*/

	int ID = 0;
	int pow_num = 0;
	for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++) {
		if (graph_hash_of_mixed_weighted_contain_edge(group_graph, vertex, *it)) { // vertex is in group *it
			ID = ID + pow(2, pow_num);
		}
		pow_num++;
	}

	return ID;

}
#pragma endregion graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_vertex_group_set_ID

#pragma region
int graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_vertex_group_vertex_2_group_set_ID(int group_vertex,
	std::unordered_set<int>& cumpulsory_group_vertices) {

	/*time complexity: O(|Gamma|); this function returns the maximum group set ID for a single vertex*/

	int pow_num = 0;
	for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++) {
		if (*it == group_vertex) {
			return pow(2, pow_num);
		}
		pow_num++;
	}
}
#pragma endregion graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_vertex_group_vertex_2_group_set_ID

#pragma region
void graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_covered_uncovered_groups(int group_sets_ID_range,
	std::unordered_set<int>& cumpulsory_group_vertices, std::unordered_map<int, std::unordered_set<int>>& covered_groups, std::unordered_map<int, std::unordered_set<int>>& uncovered_groups) {

	/*time complexity: O(|Gamma|*2^|Gamma|); for each p \in [1,group_sets_ID_range], this function calculate the groups that have not been coverred by p*/

	for (int p = 1; p <= group_sets_ID_range; p++) {

		std::unordered_set<int> c_groups, unc_groups;

		int pow_num = 0;
		for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++) {
			int id = pow(2, pow_num);
			if ((id | p) != p) { // id is not covered by p
				unc_groups.insert(*it); // *it is a group not covered by p
			}
			else {
				c_groups.insert(*it); // *it is a group covered by p
			}
			pow_num++;
		}
		covered_groups[p] = c_groups;
		uncovered_groups[p] = unc_groups;
	}

}
#pragma endregion graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_covered_uncovered_groups

#pragma region
std::unordered_map<int, std::unordered_set<int>> graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_uncovered_groups(int group_sets_ID_range,
	std::unordered_set<int>& cumpulsory_group_vertices) {

	/*time complexity: O(|Gamma|*2^|Gamma|); for each p \in [1,group_sets_ID_range], this function calculate the groups that have not been coverred by p*/

	std::unordered_map<int, std::unordered_set<int>> uncovered_groups; // <p, <uncovered_groups>>

	for (int p = 1; p <= group_sets_ID_range; p++) {

		std::unordered_set<int> groups;

		int pow_num = 0;
		for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++) {
			int id = pow(2, pow_num);
			if ((id | p) != p) { // id is not covered by p
				groups.insert(*it); // *it is a group not covered by p
			}
			pow_num++;
		}

		uncovered_groups[p] = groups;

	}

	return uncovered_groups;
}
#pragma endregion graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_uncovered_groups

#pragma region
graph_hash_of_mixed_weighted graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_build_tree(int root_v, int root_p, graph_hash_of_mixed_weighted& input_graph,
	std::unordered_map<int, std::unordered_map<int, graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_tree_node>>& trees) {

	/*this function builds tree T(v,p) at a cost of O(|V|)*/

	graph_hash_of_mixed_weighted solution_tree;

	std::queue<pair<int, int>> waited_to_processed_trees; // <v, p>
	waited_to_processed_trees.push({ root_v,root_p });

	while (waited_to_processed_trees.size() > 0) {

		int v = waited_to_processed_trees.front().first, p = waited_to_processed_trees.front().second;
		waited_to_processed_trees.pop();

		/*insert v*/
		double w_v = input_graph.hash_of_vectors[v].vertex_weight;
		graph_hash_of_mixed_weighted_add_vertex(solution_tree, v, w_v);

		int form_type = trees[v][p].type;
		if (form_type == 0) { // T(v,p) is a single vertex		
		}
		else if (form_type == 1) { // T(v,p) is formed by grow
			int u = trees[v][p].u;
			waited_to_processed_trees.push({ u,p });
			/*insert (u,v); no need to insert weight of u here, which will be inserted later for T(u,p)*/
			double c_uv = graph_hash_of_mixed_weighted_edge_weight(input_graph, u, v);
			graph_hash_of_mixed_weighted_add_edge(solution_tree, u, v, c_uv);
		}
		else { // T(v,p) is formed by merge
			int p1 = trees[v][p].p1, p2 = trees[v][p].p2;
			waited_to_processed_trees.push({ v,p1 });
			waited_to_processed_trees.push({ v,p2 });
		}

	}

	return solution_tree;

}
#pragma endregion graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_build_tree

#pragma region
std::unordered_map<int, std::unordered_map<int, double>> virtual_node_distances(graph_hash_of_mixed_weighted& group_graph, std::unordered_set<int>& cumpulsory_group_vertices,
	std::unordered_map<int, pair<std::unordered_map<int, int>, std::unordered_map<int, double>>>& SPs_to_groups) {

	std::unordered_map<int, std::unordered_map<int, double>> virtual_node_distances;

	for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++) {
		int g1 = *it;
		for (auto it2 = cumpulsory_group_vertices.begin(); it2 != cumpulsory_group_vertices.end(); it2++) {
			int g2 = *it2;
			if (g1 <= g2) {
				auto xx = SPs_to_groups.find(g1);
				double distance = INT_MAX;
				std::vector<int> g_2_vertices = group_graph.adj_v(g2);
				for (int i = g_2_vertices.size() - 1; i >= 0; i--) {
					double dis = xx->second.second[g_2_vertices[i]];
					if (dis < distance) {
						distance = dis;
					}
				}
				virtual_node_distances[g1][g2] = distance;
				virtual_node_distances[g2][g1] = distance;
			}
		}
	}

	return virtual_node_distances;
}
#pragma endregion virtual_node_distances

#pragma region
struct AllPaths_min_node {
	int v_i, v_j;
	int Xslash; // group_set_ID
	double priority_value; // W(v_i,v_j,X)
};
bool operator<(AllPaths_min_node const& x, AllPaths_min_node const& y) {
	return x.priority_value > y.priority_value; // < is the max-heap; > is the min heap; PriorityQueue is expected to be a max-heap of integer values
}
typedef typename boost::heap::fibonacci_heap<AllPaths_min_node>::handle_type handle_AllPaths_min_node;


std::unordered_map<string, double> AllPaths(std::unordered_set<int>& cumpulsory_group_vertices, std::unordered_map<int, std::unordered_map<int, double>>& virtual_node_distances) {

	std::unordered_map<string, double> W; // String is ("v_i" + "_" + "v_j" + "_" + "group set ID")

	boost::heap::fibonacci_heap<AllPaths_min_node> Q; // min queue
	std::unordered_map<string, double> Q_priorities;
	std::unordered_map<string, handle_AllPaths_min_node> Q_handles; // key is String is ("v_i" + "_" + "v_j" + "_" + "group set ID")

	/*D records the popped out optimal subtrees; String is ("v_i" + "_" + "v_j" + "_" + "group set ID") */
	std::unordered_set<string> D;

	for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++) {
		int p = *it;
		//int Xslash = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_vertex_group_vertex_2_group_set_ID(p, cumpulsory_group_vertices);
		AllPaths_min_node x;


		//x.v_i = p;
		//x.v_j = p;
		//x.Xslash = Xslash;
		//x.priority_value = 0;
		//string handle_ID = to_string(p) + "_" + to_string(p) + "_" + to_string(Xslash);
		//Q_handles[handle_ID] = Q.push(x);
		//Q_priorities[handle_ID] = 0;

		/*the following code for Xslash=0 is not in 2016 paper, but is necessary for computing every W(v_i,v_j,X)*/
		x.v_i = p;
		x.v_j = p;
		x.Xslash = 0;
		x.priority_value = 0;
		string handle_ID = to_string(p) + "_" + to_string(p) + "_" + to_string(0);
		Q_handles[handle_ID] = Q.push(x);
		Q_priorities[handle_ID] = 0;
	}

	while (Q.size() > 0) {

		AllPaths_min_node top_node = Q.top();
		int v_i = top_node.v_i, v_j = top_node.v_j, Xslash = top_node.Xslash;
		double cost = top_node.priority_value;
		Q.pop();

		string handle_ID = to_string(v_i) + "_" + to_string(v_j) + "_" + to_string(Xslash);
		W[handle_ID] = cost;
		D.insert(handle_ID);
		Q_handles.erase(handle_ID);
		Q_priorities.erase(handle_ID);
		//cout << "Q pop " + handle_ID << " priority " << cost << endl;


		/*the following code is not in 2016 paper, but is necessary for computing every W(v_i,v_j,X)*/
		for (int i = 0; i < Xslash; i++) {
			if ((i | Xslash) == Xslash) { 
				handle_ID = to_string(v_i) + "_" + to_string(v_j) + "_" + to_string(i);
				if (W.count(handle_ID) == 0) {
					W[handle_ID] = cost;
				}
				else {
					if (W[handle_ID] > cost) {
						W[handle_ID] = cost;
					}
				}
			}
		}



		for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++) {
			int p = *it;
			int p_setID = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_vertex_group_vertex_2_group_set_ID(p, cumpulsory_group_vertices);

			if ((p_setID | Xslash) != Xslash) { // p_setID is not covered by Xslash

				int Xwave = Xslash + p_setID;

				double cost_wave = cost + virtual_node_distances[v_j][p];
				handle_ID = to_string(v_i) + "_" + to_string(p) + "_" + to_string(Xwave);

				if (D.count(handle_ID) > 0) {
					continue;
				}

				AllPaths_min_node x;
				x.v_i = v_i;
				x.v_j = p;
				x.Xslash = Xwave;
				x.priority_value = cost_wave;

				if (Q_handles.count(handle_ID) == 0) {
					Q_handles[handle_ID] = Q.push(x);
					Q_priorities[handle_ID] = cost_wave;
					//cout << "Q push " + handle_ID << " priority " << cost_wave << endl;
				}
				else {
					if (cost_wave < Q_priorities[handle_ID]) {
						Q.update(Q_handles[handle_ID], x); // O(1) for decrease key
						Q_priorities[handle_ID] = cost_wave;
						//cout << "Q update " + handle_ID << " priority " << cost_wave << endl;
					}
				}
			}
		}
	}

	//cout << "W.size(): " << W.size() << endl;

	return W;

}
#pragma endregion AllPaths

#pragma region
double LB_procedure(int v, int X, double cost, int group_sets_ID_range, std::unordered_map<int, std::unordered_set<int>>& uncovered_groups,
	std::unordered_map<int, pair<std::unordered_map<int, int>, std::unordered_map<int, double>>>& SPs_to_groups,
	std::unordered_map<string, double>& W, std::unordered_map<string, double>& W2) {

	double lb = 0; // lb should be lower bound cost of a feasible solution contains T(v,X)

	if (group_sets_ID_range != X) {
		int X_slash = group_sets_ID_range - X; // X_slash \cup X equals group_sets_ID_range
		std::unordered_set<int> groups_in_X_slash = uncovered_groups[X];

		double lb1 = INT_MAX, lb2 = -1, lb_one_label = -1;
		for (auto it = groups_in_X_slash.begin(); it != groups_in_X_slash.end(); it++) {
			int i = *it;
			double dis_v_i = SPs_to_groups[i].second[v];
			double xxx = INT_MAX;
			for (auto it2 = groups_in_X_slash.begin(); it2 != groups_in_X_slash.end(); it2++) {
				int j = *it2;
				double dis_v_j = SPs_to_groups[j].second[v];
				if (xxx > dis_v_j) {
					xxx = dis_v_j;
				}
				double lb1_value = (dis_v_i + W[to_string(i) + "_" + to_string(j) + "_" + to_string(X_slash)] + dis_v_j) / 2;
				if (lb1 > lb1_value) {
					lb1 = lb1_value;
				}
			}
			double lb2_value = (dis_v_i + W2[to_string(i) + "_" + to_string(X_slash)] + xxx) / 2;
			if (lb2 < lb2_value) {
				lb2 = lb2_value;
			}
			if (lb_one_label < dis_v_i) {
				lb_one_label = dis_v_i;
			}
		}

		if (lb1 < lb2) {
			lb = lb2;
		}
		else {
			lb = lb1;
		}
		if (lb < lb_one_label) {
			lb = lb_one_label;
		}

		//cout << "lb_one_label=" << lb_one_label << endl;
	}

	return cost + lb;


}
#pragma endregion LB_procedure

graph_hash_of_mixed_weighted graph_hash_of_mixed_weighted_Basic_edge_weighted(graph_hash_of_mixed_weighted& input_graph, graph_hash_of_mixed_weighted& group_graph,
	std::unordered_set<int>& cumpulsory_group_vertices, double maximum_return_app_ratio) {

	/*this function returns the first found feasible solution that has an approximation ratio not larger than maximum_return_app_ratio*/

	double error_safeguard = 1e-4;

	if (cumpulsory_group_vertices.size() >= 20) {
		std::cout << "cumpulsory_group_vertices.size() is too large for graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted!" << std::endl;
		exit(1);
	}


	/*finding lowest-weighted paths from groups to vertices;
	time complexity: O(|T||E|+|T||V|log|V|);
	return {g_ID, { distances, predecessors }} */
	std::unordered_map<int, pair<std::unordered_map<int, int>, std::unordered_map<int, double>>> SPs_to_groups =
		graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_find_SPs(input_graph, group_graph, cumpulsory_group_vertices);

	/*initialize Q*/
	boost::heap::fibonacci_heap<graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node> Q_T; // min queues of trees
	std::unordered_map<string, handle_graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node>
		Q_T_handles; // key is string "v_p" ("vertex_ID" + "_" + "group set ID"); Q_T_handles only contain keys that are in Q_T

	double inf = std::numeric_limits<double>::max(); // cost of empty tree is inf

	/* this is the cost of the best found solution yet */
	graph_hash_of_mixed_weighted best_solu;
	double best_cost = inf;

	/*initialize non_overlapped_group_sets; time complexity: O(4^|Gamma|);
	Group		G1	G0	group_set_ID
				0   0   0
				0	1	1
				1	0	2
				1	1	3*/
	int group_sets_ID_range = pow(2, cumpulsory_group_vertices.size()) - 1; // the number of group sets: 2^|Gamma|, including empty set;   |Gamma| should be smaller than 31 due to precision
	std::unordered_map<int, std::unordered_set<int>> non_overlapped_group_sets_IDs =
		graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_non_overlapped_group_sets(group_sets_ID_range);  // time complexity: O(4^|Gamma|)


	/*initialize uncovered_groups;  <p, <uncovered_groups>>;  time complexity: O(|Gamma|*2^|Gamma|)*/
	std::unordered_map<int, std::unordered_set<int>> uncovered_groups = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_uncovered_groups(
		group_sets_ID_range, cumpulsory_group_vertices);






	/*initialize trees with vertices;
	time complexity: O(2^|Gamma|*|V|);
	every vertex v is associated with T(v,p) only when v covers p, otherwise the cost of T(v,p) is considered as inf;
	every vertex v is associated with at most 2^|Gamma| trees;
	*/
	std::unordered_map<int, std::unordered_map<int, graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_tree_node>> trees; // <v, <p, T(v,p)>>
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		int v = it->first; // a vertex
		trees[v] = {};
		int group_set_ID_v = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_vertex_group_set_ID(v, group_graph, cumpulsory_group_vertices); /*time complexity: O(|Gamma|)*/
		for (int p = 1; p <= group_set_ID_v; p++) { // p is non-empty; time complexity: O(2^|Gamma|)
			if ((p | group_set_ID_v) == group_set_ID_v) { // p represents a non-empty group set inside group_set_ID_v, including group_set_ID_v

				/*T(v,p)*/
				graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_tree_node node;
				node.cost = 0;
				node.type = 0;
				trees[v][p] = node;

				/*insert T(v,p) into Q_T*/
				graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
				x.v = v;
				x.p = p;
				x.priority_value = 0;
				string handle_ID = to_string(v) + "_" + to_string(p);
				Q_T_handles[handle_ID] = Q_T.push(x);
			}
		}
	}








	/*D records the popped out optimal subtrees; String is "v_p" ("vertex_ID" + "_" + "group set ID") */
	std::unordered_set<string> D;



	/*Big while loop*/
	while (Q_T.size() > 0) { // at most 2^|Gamma|*V loops


		graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node top_node = Q_T.top();
		int v = top_node.v, p = top_node.p;

		Q_T.pop(); // O(2^|Gamma|*V*(|Gamma| + log V)) throught the loop, since Q_T contains at most 2^|Gamma|*V elements
		string handle_ID = to_string(v) + "_" + to_string(p);
		Q_T_handles.erase(handle_ID); // Q_T_handles only contains handles of elements in Q_T

		D.insert(handle_ID); // optimal T(v,p) has been found
		//cout << "D.insert " << handle_ID << endl;

		if (p == group_sets_ID_range) { // T(v,p) covers all groups 	
			graph_hash_of_mixed_weighted feasible_solu = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_build_tree(v, p, input_graph, trees); // time complexity: O(|V|)
			return feasible_solu;
		}

		/*build a feasible solution, report app ratio, and update best; O(|T||V| + |V|log|V|)*/
		double feasible_solu_cost = top_node.priority_value;
		for (auto it = uncovered_groups[p].begin(); it != uncovered_groups[p].end(); it++) {
			feasible_solu_cost = feasible_solu_cost + SPs_to_groups[*it].second[v];
		}
		if (feasible_solu_cost < best_cost) { // best_cost is also updated in merge!
			graph_hash_of_mixed_weighted feasible_solu = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_build_tree(v, p, input_graph, trees); // time complexity: O(|V|)
			for (auto it = uncovered_groups[p].begin(); it != uncovered_groups[p].end(); it++) {
				int g_id = *it;
				/*merge LWP(v to g_id) into feasible_solu; O(|V|)*/
				int v_start = v;
				int pre = SPs_to_groups[g_id].first[v_start]; // LWPs_to_groups[g_id].second is the predecessor index
				while (v_start != pre) {
					double w_v_pre = input_graph.hash_of_vectors[pre].vertex_weight;
					double ec = graph_hash_of_mixed_weighted_edge_weight(input_graph, v_start, pre);
					graph_hash_of_mixed_weighted_add_vertex(feasible_solu, pre, w_v_pre);
					graph_hash_of_mixed_weighted_add_edge(feasible_solu, v_start, pre, ec);
					v_start = pre;
					pre = SPs_to_groups[g_id].first[v_start];
				}
			}
			best_solu = graph_hash_of_mixed_weighted_MST_postprocessing_no_whole_graph(feasible_solu);
			best_cost = graph_hash_of_mixed_weighted_sum_of_ec(best_solu);
		}


		double T_v_p_cost = trees[v][p].cost; // since an optimal solution has not been popped out, the optimal cost must be larger than or equal to T_v_p_cost
		if (T_v_p_cost > 0) {
			double ratio = best_cost / T_v_p_cost;
			if (ratio <= maximum_return_app_ratio) { // this feasible solution can be returned
				return best_solu;
			}
		}

		/*grow*/
		std::vector<int> v_adjs = input_graph.adj_v(v);
		for (auto it = v_adjs.begin(); it != v_adjs.end(); it++) {

			/*below: O(2^|Gamma|*E) in all loops, since each v has 2^|Gamma| times*/
			int u = *it;
			double cost_euv = graph_hash_of_mixed_weighted_edge_weight(input_graph, u, v);
			double grow_tree_cost = trees[v][p].cost + cost_euv;

			handle_ID = to_string(u) + "_" + to_string(p);
			//cout << "grow " << handle_ID << endl;
			//cout << "grow_tree_cost " << grow_tree_cost << endl;
			//cout << "best_cost " << best_cost << endl;
			if (D.count(handle_ID) > 0 || grow_tree_cost > best_cost + error_safeguard) { // error_safeguard is error
				//cout << "D.count(handle_ID) " << D.count(handle_ID) << endl;
				//cout << grow_tree_cost - best_cost << endl;
				continue;
			}

			double T_up_cost;
			if (trees[u].count(p) == 0) {
				T_up_cost = inf;
			}
			else {
				T_up_cost = trees[u][p].cost;
			}

			//cout << "T_up_cost " << T_up_cost << endl;

			if (grow_tree_cost < T_up_cost) {

				/*below: O(2^|Gamma|*V*(|Gamma| + log V)) throught the loop, since each u is checked 2^|Gamma| times, and Q_T contains at most 2^|Gamma|*V elements */

				/*update T(u,p) by grow T(v,p) with (u,v)*/
				trees[u][p].cost = grow_tree_cost;
				trees[u][p].type = 1;
				trees[u][p].u = v;

				/*update T(u,p) in Q_T*/
				graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
				x.v = u;
				x.p = p;
				x.priority_value = grow_tree_cost;
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
			handle_ID = to_string(v) + "_" + to_string(p2);
			if (D.count(handle_ID) > 0) { // only merge optimal T(v,p2), this is recorded in the pseudo code of 2016 paper

				int p1_cup_p2 = p1 + p2;
				handle_ID = to_string(v) + "_" + to_string(p1_cup_p2);
				if (D.count(handle_ID) > 0) {
					continue;
				}

				double cost_Tvp1 = trees[v][p1].cost;
				double cost_Tvp2 = trees[v][p2].cost;
				double cost_Tvp1_cup_p2;
				if (trees[v].count(p1_cup_p2) == 0) {
					cost_Tvp1_cup_p2 = inf;
				}
				else {
					cost_Tvp1_cup_p2 = trees[v][p1_cup_p2].cost;
				}
				double merged_tree_cost = cost_Tvp1 + cost_Tvp2;

				if (merged_tree_cost > best_cost + error_safeguard) {
					continue;
				}



				if (merged_tree_cost < cost_Tvp1_cup_p2) { // O(3^|Gamma||V| comparisons in totel, see the DPBF paper)

					/*update T(v,p1_cup_p2) by merge T(v,p1) with T(v,v2)*/
					trees[v][p1_cup_p2].cost = merged_tree_cost;
					trees[v][p1_cup_p2].type = 2;
					trees[v][p1_cup_p2].p1 = p1;
					trees[v][p1_cup_p2].p2 = p2;

					/*update T(v,p1_cup_p2) in Q_T*/
					graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
					x.v = v;
					x.p = p1_cup_p2;
					x.priority_value = merged_tree_cost;
					if (Q_T_handles.count(handle_ID) == 0) { // T(v,p1_cup_p2) is not in Q_T
						Q_T_handles[handle_ID] = Q_T.push(x);
					}
					else { // T(v,p1_cup_p2) is in Q_T
						Q_T.update(Q_T_handles[handle_ID], x); // O(1) for decrease key
					}

				}

				if (p1_cup_p2 == group_sets_ID_range) {
					best_cost = merged_tree_cost;
					best_solu = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_build_tree(v, p1_cup_p2, input_graph, trees);
				}

			}



		}

	}



	std::cout << "graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_vertex_edge_weighted_ProgressiveReturn did not find a feasible solution!" << std::endl;
	getchar();
	exit(1);

}


graph_hash_of_mixed_weighted graph_hash_of_mixed_weighted_PrunedDP_edge_weighted(graph_hash_of_mixed_weighted& input_graph, graph_hash_of_mixed_weighted& group_graph,
	std::unordered_set<int>& cumpulsory_group_vertices, double maximum_return_app_ratio) {

	/*this function returns the first found feasible solution that has an approximation ratio not larger than maximum_return_app_ratio*/

	if (cumpulsory_group_vertices.size() >= 20) {
		std::cout << "cumpulsory_group_vertices.size() is too large for graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted!" << std::endl;
		exit(1);
	}

	double error_safeguard = 1e-4;

	/*finding lowest-weighted paths from groups to vertices; time complexity: O(|T||E|+|T||V|log|V|); return {g_ID, { distances, predecessors }} */
	std::unordered_map<int, pair<std::unordered_map<int, int>, std::unordered_map<int, double>>> SPs_to_groups =
		graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_find_SPs(input_graph, group_graph, cumpulsory_group_vertices);

	/*initialize Q*/
	boost::heap::fibonacci_heap<graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node> Q_T; // min queues of trees
	std::unordered_map<string, handle_graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node> Q_T_handles; // key is string "v_p" ("vertex_ID" + "_" + "group set ID"); Q_T_handles only contain keys that are in Q_T

	double inf = std::numeric_limits<double>::max(); // cost of empty tree is inf

	/* this is the cost of the best found solution yet */
	graph_hash_of_mixed_weighted best_solu;
	double best_cost = inf;

	/*initialize non_overlapped_group_sets; time complexity: O(4^|Gamma|);
	Group		G1	G0	group_set_ID
				0   0   0
				0	1	1
				1	0	2
				1	1	3*/
	int group_sets_ID_range = pow(2, cumpulsory_group_vertices.size()) - 1; // the number of group sets: 2^|Gamma|, including empty set;   |Gamma| should be smaller than 31 due to precision
	std::unordered_map<int, std::unordered_set<int>> non_overlapped_group_sets_IDs =
		graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_non_overlapped_group_sets(group_sets_ID_range);  // time complexity: O(4^|Gamma|)

	/*initialize uncovered_groups;  <p, <uncovered_groups>>;  time complexity: O(|Gamma|*2^|Gamma|)*/
	std::unordered_map<int, std::unordered_set<int>> uncovered_groups = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_uncovered_groups(group_sets_ID_range, cumpulsory_group_vertices);


	/*initialize trees with vertices;
	time complexity: O(2^|Gamma|*|V|);
	every vertex v is associated with T(v,p) only when v covers p, otherwise the cost of T(v,p) is considered as inf;
	every vertex v is associated with at most 2^|Gamma| trees;*/
	std::unordered_map<int, std::unordered_map<int, graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_tree_node>> trees; // <v, <p, T(v,p)>>
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		int v = it->first; // a vertex
		trees[v] = {};
		int group_set_ID_v = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_vertex_group_set_ID(v, group_graph, cumpulsory_group_vertices); /*time complexity: O(|Gamma|)*/
		for (int p = 1; p <= group_set_ID_v; p++) { // p is non-empty; time complexity: O(2^|Gamma|)
			if ((p | group_set_ID_v) == group_set_ID_v) { // p represents a non-empty group set inside group_set_ID_v, including group_set_ID_v

				/*T(v,p)*/
				graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_tree_node node;
				node.cost = 0;
				node.type = 0;
				trees[v][p] = node;

				/*insert T(v,p) into Q_T*/
				graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
				x.v = v;
				x.p = p;
				x.priority_value = 0;
				string handle_ID = to_string(v) + "_" + to_string(p);
				Q_T_handles[handle_ID] = Q_T.push(x);
			}
		}
	}



	/*D records the popped out optimal subtrees; String is "v_p" ("vertex_ID" + "_" + "group set ID") */
	std::unordered_set<string> D;

	/*Big while loop*/
	while (Q_T.size() > 0) { // at most 2^|Gamma|*V loops


		graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node top_node = Q_T.top();
		int v = top_node.v, p = top_node.p;

		Q_T.pop(); // O(2^|Gamma|*V*(|Gamma| + log V)) throught the loop, since Q_T contains at most 2^|Gamma|*V elements
		string handle_ID = to_string(v) + "_" + to_string(p);
		Q_T_handles.erase(handle_ID); // Q_T_handles only contains handles of elements in Q_T

		D.insert(handle_ID); // optimal T(v,p) has been found
		//cout << "D.insert " << handle_ID << endl;

		if (p == group_sets_ID_range) { // T(v,p) covers all groups 	
			graph_hash_of_mixed_weighted feasible_solu = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_build_tree(v, p, input_graph, trees); // time complexity: O(|V|)
			return feasible_solu;
		}

		/*build a feasible solution, report app ratio, and update best; O(|T||V| + |V|log|V|)*/
		double feasible_solu_cost = top_node.priority_value;
		for (auto it = uncovered_groups[p].begin(); it != uncovered_groups[p].end(); it++) {
			feasible_solu_cost = feasible_solu_cost + SPs_to_groups[*it].second[v];
		}
		if (feasible_solu_cost < best_cost) { // best_cost is also updated in merge!
			graph_hash_of_mixed_weighted feasible_solu = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_build_tree(v, p, input_graph, trees); // time complexity: O(|V|)
			for (auto it = uncovered_groups[p].begin(); it != uncovered_groups[p].end(); it++) {
				int g_id = *it;
				/*merge LWP(v to g_id) into feasible_solu; O(|V|)*/
				int v_start = v;
				int pre = SPs_to_groups[g_id].first[v_start]; // LWPs_to_groups[g_id].second is the predecessor index
				while (v_start != pre) {
					double w_v_pre = input_graph.hash_of_vectors[pre].vertex_weight;
					double ec = graph_hash_of_mixed_weighted_edge_weight(input_graph, v_start, pre);
					graph_hash_of_mixed_weighted_add_vertex(feasible_solu, pre, w_v_pre);
					graph_hash_of_mixed_weighted_add_edge(feasible_solu, v_start, pre, ec);
					v_start = pre;
					pre = SPs_to_groups[g_id].first[v_start];
				}
			}
			best_solu = graph_hash_of_mixed_weighted_MST_postprocessing_no_whole_graph(feasible_solu);
			best_cost = graph_hash_of_mixed_weighted_sum_of_ec(best_solu);
		}


		double T_v_p_cost = trees[v][p].cost; // since an optimal solution has not been popped out, the optimal cost must be larger than or equal to T_v_p_cost
		if (T_v_p_cost > 0) {
			double ratio = best_cost / T_v_p_cost;
			if (ratio <= maximum_return_app_ratio) { // this feasible solution can be returned
				return best_solu;
			}
		}


		/*Lines 16-18 in PrunedDP*/
		int p_slash = group_sets_ID_range - p; // p \cup p_slash equals group_sets_ID_range
		handle_ID = to_string(v) + "_" + to_string(p_slash);
		if (D.count(handle_ID) > 0) {
			handle_ID = to_string(v) + "_" + to_string(group_sets_ID_range);
			double cost_Tvp1 = trees[v][p].cost;
			double cost_Tvp2 = trees[v][p_slash].cost;

			double cost_Tvp1_cup_p2;
			if (trees[v].count(group_sets_ID_range) == 0) {
				cost_Tvp1_cup_p2 = inf;
			}
			else {
				cost_Tvp1_cup_p2 = trees[v][group_sets_ID_range].cost;
			}
			double merged_tree_cost = cost_Tvp1 + cost_Tvp2;

			if (merged_tree_cost > best_cost + error_safeguard) { 
				continue;
			}

			if (merged_tree_cost < cost_Tvp1_cup_p2) {

				/*update T(v,p1_cup_p2) by merge T(v,p1) with T(v,v2)*/
				trees[v][group_sets_ID_range].cost = merged_tree_cost;
				trees[v][group_sets_ID_range].type = 2;
				trees[v][group_sets_ID_range].p1 = p;
				trees[v][group_sets_ID_range].p2 = p_slash;

				/*update T(v,p1_cup_p2) in Q_T*/
				graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
				x.v = v;
				x.p = group_sets_ID_range;
				x.priority_value = merged_tree_cost;
				if (Q_T_handles.count(handle_ID) == 0) { // T(v,p1_cup_p2) is not in Q_T
					Q_T_handles[handle_ID] = Q_T.push(x);
				}
				else { // T(v,p1_cup_p2) is in Q_T
					Q_T.update(Q_T_handles[handle_ID], x); // O(1) for decrease key
				}

			}

			best_cost = merged_tree_cost;
			best_solu = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_build_tree(v, group_sets_ID_range, input_graph, trees);

			continue;
		}




		/* Optimal-Tree Decomposition Theorem (Theorem 1 in 2016 SIGMOD paper);
		note the following place of error_safeguard; this if condition should be met even if error_safeguard=inf*/
		if (T_v_p_cost < best_cost / 2 + error_safeguard) {

			/*grow*/
			std::vector<int> v_adjs = input_graph.adj_v(v);
			for (auto it = v_adjs.begin(); it != v_adjs.end(); it++) {

				/*below: O(2^|Gamma|*E) in all loops, since each v has 2^|Gamma| times*/
				int u = *it;
				double cost_euv = graph_hash_of_mixed_weighted_edge_weight(input_graph, u, v);
				double grow_tree_cost = trees[v][p].cost + cost_euv;

				handle_ID = to_string(u) + "_" + to_string(p);
				//cout << "grow " << handle_ID << endl;
				//cout << "grow_tree_cost " << grow_tree_cost << endl;
				//cout << "best_cost " << best_cost << endl;
				if (D.count(handle_ID) > 0 || grow_tree_cost > best_cost + error_safeguard) {
					//cout << "D.count(handle_ID) " << D.count(handle_ID) << endl;
					//cout << grow_tree_cost - best_cost << endl;
					continue;
				}

				double T_up_cost;
				if (trees[u].count(p) == 0) {
					T_up_cost = inf;
				}
				else {
					T_up_cost = trees[u][p].cost;
				}

				//cout << "T_up_cost " << T_up_cost << endl;

				if (grow_tree_cost < T_up_cost) {

					/*below: O(2^|Gamma|*V*(|Gamma| + log V)) throught the loop, since each u is checked 2^|Gamma| times, and Q_T contains at most 2^|Gamma|*V elements */

					/*update T(u,p) by grow T(v,p) with (u,v)*/
					trees[u][p].cost = grow_tree_cost;
					trees[u][p].type = 1;
					trees[u][p].u = v;

					/*update T(u,p) in Q_T*/
					graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
					x.v = u;
					x.p = p;
					x.priority_value = grow_tree_cost;
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
				handle_ID = to_string(v) + "_" + to_string(p2);
				if (D.count(handle_ID) > 0) { // only merge optimal T(v,p2), this is recorded in the pseudo code of 2016 paper

					int p1_cup_p2 = p1 + p2;
					handle_ID = to_string(v) + "_" + to_string(p1_cup_p2);
					if (D.count(handle_ID) > 0) {
						continue;
					}

					double cost_Tvp1 = trees[v][p1].cost;
					double cost_Tvp2 = trees[v][p2].cost;

					double cost_Tvp1_cup_p2;
					if (trees[v].count(p1_cup_p2) == 0) {
						cost_Tvp1_cup_p2 = inf;
					}
					else {
						cost_Tvp1_cup_p2 = trees[v][p1_cup_p2].cost;
					}
					double merged_tree_cost = cost_Tvp1 + cost_Tvp2;

					/*Conditional Tree Merging Theorem (Theorem 2 in 2016 SIGMOD paper)*/
					if (merged_tree_cost > best_cost * 2 / 3 + error_safeguard) { // error_safeguard is error
						continue;
					}

					if (merged_tree_cost < cost_Tvp1_cup_p2) { // O(3^|Gamma||V| comparisons in totel, see the DPBF paper)

						/*update T(v,p1_cup_p2) by merge T(v,p1) with T(v,v2)*/
						trees[v][p1_cup_p2].cost = merged_tree_cost;
						trees[v][p1_cup_p2].type = 2;
						trees[v][p1_cup_p2].p1 = p1;
						trees[v][p1_cup_p2].p2 = p2;

						/*update T(v,p1_cup_p2) in Q_T*/
						graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
						x.v = v;
						x.p = p1_cup_p2;
						x.priority_value = merged_tree_cost;
						if (Q_T_handles.count(handle_ID) == 0) { // T(v,p1_cup_p2) is not in Q_T
							Q_T_handles[handle_ID] = Q_T.push(x);
						}
						else { // T(v,p1_cup_p2) is in Q_T
							Q_T.update(Q_T_handles[handle_ID], x); // O(1) for decrease key
						}

					}

					if (p1_cup_p2 == group_sets_ID_range) {
						best_cost = merged_tree_cost;
						best_solu = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_build_tree(v, p1_cup_p2, input_graph, trees);
					}

				}



			}
		}



	}



	std::cout << "graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_vertex_edge_weighted_ProgressiveReturn did not find a feasible solution!" << std::endl;
	getchar();
	exit(1);

}



graph_hash_of_mixed_weighted graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted(graph_hash_of_mixed_weighted& input_graph, graph_hash_of_mixed_weighted& group_graph,
	std::unordered_set<int>& cumpulsory_group_vertices, double maximum_return_app_ratio) {

	/*this function returns the first found feasible solution that has an approximation ratio not larger than maximum_return_app_ratio*/

	if (cumpulsory_group_vertices.size() >= 15) {
		std::cout << "cumpulsory_group_vertices.size() is too large!" << std::endl;
		exit(1);
	}

	double error_safeguard = 1e-3;

	/*finding lowest-weighted paths from groups to vertices; time complexity: O(|T||E|+|T||V|log|V|); return {g_ID, { distances, predecessors }} */
	std::unordered_map<int, pair<std::unordered_map<int, int>, std::unordered_map<int, double>>> SPs_to_groups =
		graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_find_SPs(input_graph, group_graph, cumpulsory_group_vertices);

	/*initialize Q*/
	boost::heap::fibonacci_heap<graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node> Q_T; // min queues of trees
	std::unordered_map<string, handle_graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node> Q_T_handles; // key is string "v_p" ("vertex_ID" + "_" + "group set ID"); Q_T_handles only contain keys that are in Q_T
	std::unordered_map<string, double> Q_T_priorities; // key is string "v_p" ("vertex_ID" + "_" + "group set ID"); Q_T_priorities only contain keys that are in Q_T


	double inf = std::numeric_limits<double>::max(); // cost of empty tree is inf

	/* this is the cost of the best found solution yet */
	graph_hash_of_mixed_weighted best_solu;
	double best_cost = inf;

	/*initialize non_overlapped_group_sets; time complexity: O(4^|Gamma|);
	Group		G1	G0	group_set_ID
				0   0   0
				0	1	1
				1	0	2
				1	1	3*/
	int group_sets_ID_range = pow(2, cumpulsory_group_vertices.size()) - 1; // the number of group sets: 2^|Gamma|, including empty set;   |Gamma| should be smaller than 31 due to precision
	std::unordered_map<int, std::unordered_set<int>> non_overlapped_group_sets_IDs = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_non_overlapped_group_sets(group_sets_ID_range);  // time complexity: O(4^|Gamma|)

	/*initialize uncovered_groups;  <p, <uncovered_groups>>;  time complexity: O(|Gamma|*2^|Gamma|)*/
	std::unordered_map<int, std::unordered_set<int>> covered_groups, uncovered_groups;
	graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_covered_uncovered_groups(group_sets_ID_range, cumpulsory_group_vertices, covered_groups, uncovered_groups);
	//for (auto it = covered_groups.begin(); it != covered_groups.end(); it++) {
	//	cout << "group_set_ID " << it->first << " covered_groups:" << endl;
	//	print_unordered_set_int(it->second);
	//}

	std::unordered_map<int, std::unordered_map<int, double>> virtual_distances = virtual_node_distances(group_graph, cumpulsory_group_vertices, SPs_to_groups);
	//for (auto it = virtual_distances.begin(); it != virtual_distances.end(); it++) {
	//	cout << "virtual_node_distances from " << it->first << " to " << endl;
	//	for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
	//		cout << it2->first << " dis: " << it2->second << endl;
	//	}
	//}

	std::unordered_map<string, double> W = AllPaths(cumpulsory_group_vertices, virtual_distances); // String is ("v_i" + "_" + "v_j" + "_" + "group set ID")
	std::unordered_map<string, double> W2; // String is ("v_i" + "_" + "group set ID")
	for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++) {
		int v_i = *it;
		for (int Xslash = 1; Xslash <= group_sets_ID_range; Xslash++) {
			double dis = INT_MAX;
			for (auto it2 = cumpulsory_group_vertices.begin(); it2 != cumpulsory_group_vertices.end(); it2++) {
				int v_j = *it2;
				string handle_ID = to_string(v_i) + "_" + to_string(v_j) + "_" + to_string(Xslash);
				if (dis > W[handle_ID]) {
					dis = W[handle_ID];
				}
			}
			string handle_ID = to_string(v_i) + "_" + to_string(Xslash);
			W2[handle_ID] = dis;
		}
	}

	//print_unordered_map_string_double(W);
	//print_unordered_map_string_double(W2);
	//if (W.size() != pow(cumpulsory_group_vertices.size(), 2) * (group_sets_ID_range + 1)) {
	//	cout << "pow(cumpulsory_group_vertices.size(), 2) * (group_sets_ID_range + 1): " << pow(cumpulsory_group_vertices.size(), 2) * (group_sets_ID_range + 1) << endl;
	//	cout << "W.size(): " << W.size() << endl;
	//	getchar();
	//}


	/*initialize trees with vertices; time complexity: O(2^|Gamma|*|V|);
	every vertex v is associated with T(v,p) only when v covers p, otherwise the cost of T(v,p) is considered as inf;
	every vertex v is associated with at most 2^|Gamma| trees*/
	std::unordered_map<int, std::unordered_map<int, graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_tree_node>> trees; // <v, <p, T(v,p)>>
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		int v = it->first; // a vertex
		trees[v] = {};
		int group_set_ID_v = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_vertex_group_set_ID(v, group_graph, cumpulsory_group_vertices); /*time complexity: O(|Gamma|)*/
		for (int p = 1; p <= group_set_ID_v; p++) { // p is non-empty; time complexity: O(2^|Gamma|)
			if ((p | group_set_ID_v) == group_set_ID_v) { // p represents a non-empty group set inside group_set_ID_v, including group_set_ID_v

				/*T(v,p)*/
				graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_tree_node node;
				node.cost = 0;
				node.type = 0;
				trees[v][p] = node;

				/*insert T(v,p) into Q_T*/
				graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
				x.v = v;
				x.p = p;
				x.priority_value = 0;
				string handle_ID = to_string(v) + "_" + to_string(p);
				Q_T_handles[handle_ID] = Q_T.push(x);
				Q_T_priorities[handle_ID] = x.priority_value;
				//cout << "initial Q push " + handle_ID << " priority: " << x.priority_value << endl;
			}
		}
	}



	/*D records the popped out optimal subtrees; String is "v_p" ("vertex_ID" + "_" + "group set ID") */
	std::unordered_set<string> D;

	//cout << "group_sets_ID_range:" << group_sets_ID_range << endl;

	/*Big while loop*/
	while (Q_T.size() > 0) { // at most 2^|Gamma|*V loops


		graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node top_node = Q_T.top();
		int v = top_node.v, X = top_node.p;
		double v_X_tree_cost = trees[v][X].cost;
		Q_T.pop(); // O(2^|Gamma|*V*(|Gamma| + log V)) throught the loop, since Q_T contains at most 2^|Gamma|*V elements

		string handle_ID = to_string(v) + "_" + to_string(X);
		Q_T_handles.erase(handle_ID); // Q_T_handles only contains handles of elements in Q_T
		Q_T_priorities.erase(handle_ID);

		//cout << "X:" << X << endl;
		if (X == group_sets_ID_range) { // T(v,p) covers all groups 	
			graph_hash_of_mixed_weighted feasible_solu = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_build_tree(v, X, input_graph, trees);
			return feasible_solu;
		}

		//cout << "Q pop " + handle_ID << endl;
		D.insert(handle_ID); // optimal T(v,p) has been found

		std::unordered_set<int> X_uncovered_groups = uncovered_groups[X];

		/*build a feasible solution, report app ratio, and update best; O(|T||V| + |V|log|V|)*/
		double feasible_solu_cost = v_X_tree_cost;
		for (auto it = X_uncovered_groups.begin(); it != X_uncovered_groups.end(); it++) {
			feasible_solu_cost = feasible_solu_cost + SPs_to_groups[*it].second[v];
		}
		if (feasible_solu_cost < best_cost) {
			graph_hash_of_mixed_weighted feasible_solu = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_build_tree(v, X, input_graph, trees);
			for (auto it = X_uncovered_groups.begin(); it != X_uncovered_groups.end(); it++) {
				int g_id = *it;
				/*merge LWP(v to g_id) into feasible_solu; O(|V|)*/
				int v_start = v;
				int pre = SPs_to_groups[g_id].first[v_start];
				while (v_start != pre) {
					double w_v_pre = input_graph.hash_of_vectors[pre].vertex_weight;
					double ec = graph_hash_of_mixed_weighted_edge_weight(input_graph, v_start, pre);
					graph_hash_of_mixed_weighted_add_vertex(feasible_solu, pre, w_v_pre);
					graph_hash_of_mixed_weighted_add_edge(feasible_solu, v_start, pre, ec);
					v_start = pre;
					pre = SPs_to_groups[g_id].first[v_start];
				}
			}
			best_solu = graph_hash_of_mixed_weighted_MST_postprocessing_no_whole_graph(feasible_solu);
			best_cost = graph_hash_of_mixed_weighted_sum_of_ec(best_solu);
		}


		/*since an optimal solution has not been popped out, the optimal cost must be larger than or equal to v_X_tree_cost*/
		if (v_X_tree_cost > 0) {
			double ratio = best_cost / v_X_tree_cost;
			//cout << "ratio:" << ratio << " v_X_tree_cost:" << v_X_tree_cost << endl;
			if (ratio <= maximum_return_app_ratio + 1e-5) { // this feasible solution can be returned
				return best_solu;
			}
		}


		/*Lines 17-19 in PrunedDP++*/
		int X_slash = group_sets_ID_range - X; // p \cup X_slash equals group_sets_ID_range
		handle_ID = to_string(v) + "_" + to_string(X_slash);
		if (D.count(handle_ID) > 0) {

			double merged_tree_cost = v_X_tree_cost + trees[v][X_slash].cost;

			double cost_Tvp1_cup_p2 = inf;
			if (trees[v].count(group_sets_ID_range) > 0) {
				cost_Tvp1_cup_p2 = trees[v][group_sets_ID_range].cost;
			}
			if (merged_tree_cost < cost_Tvp1_cup_p2) { // update tree T(v,group_sets_ID_range)
				trees[v][group_sets_ID_range].cost = merged_tree_cost;
				trees[v][group_sets_ID_range].type = 2;
				trees[v][group_sets_ID_range].p1 = X;
				trees[v][group_sets_ID_range].p2 = X_slash;

				double lb = merged_tree_cost;
				//if (v_X_lb_in_Q > lb) {
				//	lb = v_X_lb_in_Q;
				//}

				if (merged_tree_cost > best_cost + error_safeguard) { // cannot have >= here, since best_solu may not be in Q; if >=, then best_solu may never be in Q
					continue;
				}
				if (merged_tree_cost < best_cost) {
					best_cost = merged_tree_cost;
					best_solu = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_build_tree(v, group_sets_ID_range, input_graph, trees);
				}

				handle_ID = to_string(v) + "_" + to_string(group_sets_ID_range);
				if (Q_T_priorities.count(handle_ID) == 0) {
					graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
					x.v = v;
					x.p = group_sets_ID_range;
					x.priority_value = merged_tree_cost;
					Q_T_handles[handle_ID] = Q_T.push(x);
					Q_T_priorities[handle_ID] = x.priority_value;
					//cout << "Q push " + handle_ID << " priority: " << Q_T_priorities[handle_ID] << endl;
				}
				else {
					if (lb < Q_T_priorities[handle_ID]) {
						graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
						x.v = v;
						x.p = group_sets_ID_range;
						x.priority_value = merged_tree_cost;
						Q_T.update(Q_T_handles[handle_ID], x); // O(1) for decrease key
						Q_T_priorities[handle_ID] = x.priority_value;
					}
				}
			}
			continue;
		}


		/* Optimal-Tree Decomposition Theorem (Theorem 1 in 2016 SIGMOD paper); optimal solution may not be found if below is / 2, why?*/
		if (v_X_tree_cost < best_cost / 2 + error_safeguard) {

			/*grow*/
			std::vector<int> v_adjs = input_graph.adj_v(v);
			for (auto it = v_adjs.begin(); it != v_adjs.end(); it++) {

				/*below: O(2^|Gamma|*E) in all loops, since each v has 2^|Gamma| times*/
				int u = *it;
				double cost_euv = graph_hash_of_mixed_weighted_edge_weight(input_graph, u, v);
				double grow_tree_cost = trees[v][X].cost + cost_euv;

				//if (v == 3 && u == 2 && X == 7) {
				//	cout << "here-2 grow u=" << u << " X=" << X << " grow_tree_cost=" << grow_tree_cost << endl;
				//	cout << "v:" << v << " trees[v][X].cost:" << trees[v][X].cost << endl;
				//	cout << "top_node.priority_value:" << top_node.priority_value << endl;
				//	cout << "group_set_ID " << X << " covered_groups:" << endl;
				//	print_unordered_set_int(covered_groups[X]);
				//	getchar();
				//}

				handle_ID = to_string(u) + "_" + to_string(X);
				//cout << "grow " << handle_ID <<  " grow_tree_cost " << grow_tree_cost << " best_cost " << best_cost << endl;
				if (D.count(handle_ID) > 0) {
					continue;
				}

				

				double T_up_cost;
				if (trees[u].count(X) == 0) {
					T_up_cost = inf;
				}
				else {
					T_up_cost = trees[u][X].cost;
				}
				if (grow_tree_cost < T_up_cost) {
					trees[u][X].cost = grow_tree_cost;
					trees[u][X].type = 1;
					trees[u][X].u = v;

					double lb = LB_procedure(u, X, grow_tree_cost, group_sets_ID_range, uncovered_groups, SPs_to_groups, W, W2);



					//vector<int> opt = { 2,0,3,9,5,7,4 }; // vertices in optimal solution
					//bool is_in = false;
					//for (int ss = 0; ss < opt.size(); ss++) {
					//	if (u == opt[ss]) {
					//		is_in = 1;
					//		break;
					//	}
					//}
					//if (is_in && lb > 1.465) { // lb is wrongly larger than optimal cost
					//	if (v == 9 && X == 7) {
					//		cout << "grow u=" << u << " X=" << X << " grow_tree_cost=" << grow_tree_cost <<
					//			" lb=" << lb << " best_cost=" << best_cost << endl;
					//		cout << "v:" << v << " trees[v][X].cost:" << trees[v][X].cost << endl;
					//		cout << "top_node.priority_value:" << top_node.priority_value << endl;
					//		cout << "group_set_ID " << X << " covered_groups:" << endl;
					//		print_unordered_set_int(covered_groups[X]);
					//		getchar();
					//	}
					//}
					//if (v == 3 && X == 7) {
					//	cout << "here-1 grow u=" << u << " X=" << X << " grow_tree_cost=" << grow_tree_cost <<
					//		" lb=" << lb << " best_cost=" << best_cost << endl;
					//	cout << "v:" << v << " trees[v][X].cost:" << trees[v][X].cost << endl;
					//	cout << "top_node.priority_value:" << top_node.priority_value << endl;
					//	cout << "group_set_ID " << X << " covered_groups:" << endl;
					//	print_unordered_set_int(covered_groups[X]);
					//	getchar();
					//}

					//cout << "grow T_up_cost=" << T_up_cost << " lb=" << lb << " best_cost=" << best_cost << endl;

					if (lb > best_cost + error_safeguard) {
						//cout << "grow T_up_cost=" << T_up_cost << " lb=" << lb << " best_cost=" << best_cost << endl;
						continue;
					}

					if (Q_T_priorities.count(handle_ID) == 0) {
						graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
						x.v = u;
						x.p = X;
						x.priority_value = grow_tree_cost;
						Q_T_handles[handle_ID] = Q_T.push(x);
						Q_T_priorities[handle_ID] = x.priority_value;

						

						//cout << "Q push " + handle_ID << " priority: " << Q_T_priorities[handle_ID] << endl;

						//if (handle_ID == "3_7") {
						//	cout << "here-6" << endl;
						//	cout << "v:" << v << " trees[v][X].cost:" << trees[v][X].cost << " grow_tree_cost:" << grow_tree_cost << endl;
						//	cout << "grow u=" << u << " X=" << X << " grow_tree_cost=" << grow_tree_cost <<
						//		" lb=" << lb << " best_cost=" << best_cost << endl;
						//	cout << "top_node.priority_value:" << top_node.priority_value << endl;
						//	cout << "group_set_ID " << X << " covered_groups:" << endl;
						//	print_unordered_set_int(covered_groups[X]);
						//	getchar();
						//}

					}
					else {
						if (grow_tree_cost < Q_T_priorities[handle_ID]) {
							graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
							x.v = u;
							x.p = X;
							x.priority_value = grow_tree_cost;
							Q_T.update(Q_T_handles[handle_ID], x); // O(1) for decrease key
							Q_T_priorities[handle_ID] = x.priority_value;
						}
					}
				}
			}


			/*merge*/
			int p1 = X;
			std::unordered_set<int> non_overlapped_group_sets_p1 = non_overlapped_group_sets_IDs[p1];
			for (auto it = non_overlapped_group_sets_p1.begin(); it != non_overlapped_group_sets_p1.end(); it++) {

				int p2 = *it; // p2 is not overlapped with p1

				if (D.count(to_string(v) + "_" + to_string(p2)) > 0) { // only merge optimal T(v,p2), this is recorded in the pseudo code of 2016 paper

					int p1_cup_p2 = p1 + p2;
					handle_ID = to_string(v) + "_" + to_string(p1_cup_p2);

					if (D.count(handle_ID) > 0) {
						continue;
					}

					double merged_tree_cost = v_X_tree_cost + trees[v][p2].cost;

					//cout << "merge " << handle_ID << " merged_tree_cost " << merged_tree_cost << " best_cost " << best_cost << endl;

					double cost_Tvp1_cup_p2;
					if (trees[v].count(p1_cup_p2) == 0) {
						cost_Tvp1_cup_p2 = inf;
					}
					else {
						cost_Tvp1_cup_p2 = trees[v][p1_cup_p2].cost;
					}

					if (merged_tree_cost < cost_Tvp1_cup_p2) { // O(3^|Gamma||V| comparisons in totel, see the DPBF paper)
						trees[v][p1_cup_p2].cost = merged_tree_cost;
						trees[v][p1_cup_p2].type = 2;
						trees[v][p1_cup_p2].p1 = p1;
						trees[v][p1_cup_p2].p2 = p2;

						

						/*Conditional Tree Merging Theorem (Theorem 2 in 2016 SIGMOD paper);
						note the following place of error_safeguard; this if condition should be met even if error_safeguard=inf*/
						if (merged_tree_cost <= best_cost * 2 / 3 + error_safeguard) { // error_safeguard is error

							double lb = LB_procedure(v, p1_cup_p2, merged_tree_cost, group_sets_ID_range, uncovered_groups, SPs_to_groups, W, W2);
							//cout << "merge cost_Tvp1_cup_p2=" << cost_Tvp1_cup_p2 << " lb_slash=" << lb_slash << " lb=" << lb << " best_cost=" << best_cost << endl;

							//if (v == 3 && p1_cup_p2 == 7) {
							//	cout << "here-5 merged_tree_cost=" << merged_tree_cost << endl;
							//	cout << "v:" << v << " trees[v][X].cost:" << trees[v][X].cost << " best_cost=" << best_cost << " lb=" << lb << endl;
							//	cout << "top_node.priority_value:" << top_node.priority_value << endl;
							//	cout << "handle_ID:" << handle_ID << endl;
							//	cout << "Q_T_priorities.count(handle_ID):" << Q_T_priorities.count(handle_ID) << endl;
							//	cout << "Q_T_priorities[handle_ID]:" << Q_T_priorities[handle_ID] << endl;
							//	cout << "group_set_ID " << X << " covered_groups:" << endl;
							//	print_unordered_set_int(covered_groups[X]);
							//	getchar();
							//}


							if (lb > best_cost + error_safeguard) {
								continue;
							}
							if (p1_cup_p2 == group_sets_ID_range && merged_tree_cost < best_cost) {
								best_cost = merged_tree_cost;
								best_solu = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_build_tree(v, group_sets_ID_range, input_graph, trees);
							}

							if (Q_T_priorities.count(handle_ID) == 0) {
								graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
								x.v = v;
								x.p = p1_cup_p2;
								x.priority_value = merged_tree_cost;
								Q_T_handles[handle_ID] = Q_T.push(x);
								Q_T_priorities[handle_ID] = x.priority_value;
								//cout << "Q push " + handle_ID << " priority: " << Q_T_priorities[handle_ID] << endl;

								//if (x.v == 3) {
								//	cout << "Q push " + handle_ID << " priority: " << Q_T_priorities[handle_ID] << endl;
								//}
								//if (handle_ID == "3_7") {
								//	cout << "here-3 merged_tree_cost=" << merged_tree_cost << endl;
								//	cout << "v:" << v << " trees[v][X].cost:" << trees[v][X].cost << endl;
								//	cout << "top_node.priority_value:" << top_node.priority_value << endl;
								//	cout << "group_set_ID " << X << " covered_groups:" << endl;
								//	print_unordered_set_int(covered_groups[X]);
								//	getchar();
								//}

							}
							else {
								if (merged_tree_cost < Q_T_priorities[handle_ID]) { // if priority is lb values, then here is lv < < Q_T_priorities[handle_ID]
									graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
									x.v = v;
									x.p = p1_cup_p2;
									x.priority_value = merged_tree_cost;
									Q_T.update(Q_T_handles[handle_ID], x); // O(1) for decrease key
									Q_T_priorities[handle_ID] = x.priority_value;
									//cout << "Q update " + handle_ID << " priority: " << Q_T_priorities[handle_ID] << endl;

									//if (handle_ID == "3_7") {
									//	cout << "here-4 merged_tree_cost=" << merged_tree_cost << endl;
									//	cout << "v:" << v << " trees[v][X].cost:" << trees[v][X].cost << endl;
									//	cout << "top_node.priority_value:" << top_node.priority_value << endl;
									//	cout << "group_set_ID " << X << " covered_groups:" << endl;
									//	print_unordered_set_int(covered_groups[X]);
									//	getchar();
									//}
								}
							}
						}
					}
				}
			}
		}
	}


	std::cout << "graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted did not find a feasible solution!" << std::endl;
	graph_hash_of_mixed_weighted_save_for_GSTP("PrunedDPPlusPlus_irregular.txt", input_graph, group_graph, cumpulsory_group_vertices, 1);
	getchar();
	exit(1);

}


/*the following graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_lb_priority is wrong; it may not be right to use lb as priority*/
graph_hash_of_mixed_weighted graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_lb_priority(graph_hash_of_mixed_weighted& input_graph, graph_hash_of_mixed_weighted& group_graph,
	std::unordered_set<int>& cumpulsory_group_vertices, double maximum_return_app_ratio) {

	/*this function returns the first found feasible solution that has an approximation ratio not larger than maximum_return_app_ratio*/

	if (cumpulsory_group_vertices.size() >= 15) {
		std::cout << "cumpulsory_group_vertices.size() is too large!" << std::endl;
		exit(1);
	}

	double error_safeguard = 1e-3;

	/*finding lowest-weighted paths from groups to vertices; time complexity: O(|T||E|+|T||V|log|V|); return {g_ID, { distances, predecessors }} */
	std::unordered_map<int, pair<std::unordered_map<int, int>, std::unordered_map<int, double>>> SPs_to_groups =
		graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_find_SPs(input_graph, group_graph, cumpulsory_group_vertices);

	/*initialize Q*/
	boost::heap::fibonacci_heap<graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node> Q_T; // min queues of trees
	std::unordered_map<string, handle_graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node> Q_T_handles; // key is string "v_p" ("vertex_ID" + "_" + "group set ID"); Q_T_handles only contain keys that are in Q_T
	std::unordered_map<string, double> Q_T_priorities; // key is string "v_p" ("vertex_ID" + "_" + "group set ID"); Q_T_priorities only contain keys that are in Q_T


	double inf = std::numeric_limits<double>::max(); // cost of empty tree is inf

	/* this is the cost of the best found solution yet */
	graph_hash_of_mixed_weighted best_solu;
	double best_cost = inf;

	/*initialize non_overlapped_group_sets; time complexity: O(4^|Gamma|);
	Group		G1	G0	group_set_ID
				0   0   0
				0	1	1
				1	0	2
				1	1	3*/
	int group_sets_ID_range = pow(2, cumpulsory_group_vertices.size()) - 1; // the number of group sets: 2^|Gamma|, including empty set;   |Gamma| should be smaller than 31 due to precision
	std::unordered_map<int, std::unordered_set<int>> non_overlapped_group_sets_IDs = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_non_overlapped_group_sets(group_sets_ID_range);  // time complexity: O(4^|Gamma|)

	/*initialize uncovered_groups;  <p, <uncovered_groups>>;  time complexity: O(|Gamma|*2^|Gamma|)*/
	std::unordered_map<int, std::unordered_set<int>> covered_groups, uncovered_groups;
	graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_covered_uncovered_groups(group_sets_ID_range, cumpulsory_group_vertices, covered_groups, uncovered_groups);
	//for (auto it = covered_groups.begin(); it != covered_groups.end(); it++) {
	//	cout << "group_set_ID " << it->first << " covered_groups:" << endl;
	//	print_unordered_set_int(it->second);
	//}

	std::unordered_map<int, std::unordered_map<int, double>> virtual_distances = virtual_node_distances(group_graph, cumpulsory_group_vertices, SPs_to_groups);
	//for (auto it = virtual_distances.begin(); it != virtual_distances.end(); it++) {
	//	cout << "virtual_node_distances from " << it->first << " to " << endl;
	//	for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
	//		cout << it2->first << " dis: " << it2->second << endl;
	//	}
	//}

	std::unordered_map<string, double> W = AllPaths(cumpulsory_group_vertices, virtual_distances); // String is ("v_i" + "_" + "v_j" + "_" + "group set ID")
	std::unordered_map<string, double> W2; // String is ("v_i" + "_" + "group set ID")
	for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++) {
		int v_i = *it;
		for (int Xslash = 1; Xslash <= group_sets_ID_range; Xslash++) {
			double dis = INT_MAX;
			for (auto it2 = cumpulsory_group_vertices.begin(); it2 != cumpulsory_group_vertices.end(); it2++) {
				int v_j = *it2;
				string handle_ID = to_string(v_i) + "_" + to_string(v_j) + "_" + to_string(Xslash);
				if (dis > W[handle_ID]) {
					dis = W[handle_ID];
				}
			}
			string handle_ID = to_string(v_i) + "_" + to_string(Xslash);
			W2[handle_ID] = dis;
		}
	}

	//print_unordered_map_string_double(W);
	//print_unordered_map_string_double(W2);
	//if (W.size() != pow(cumpulsory_group_vertices.size(), 2) * (group_sets_ID_range + 1)) {
	//	cout << "pow(cumpulsory_group_vertices.size(), 2) * (group_sets_ID_range + 1): " << pow(cumpulsory_group_vertices.size(), 2) * (group_sets_ID_range + 1) << endl;
	//	cout << "W.size(): " << W.size() << endl;
	//	getchar();
	//}


	/*initialize trees with vertices; time complexity: O(2^|Gamma|*|V|);
	every vertex v is associated with T(v,p) only when v covers p, otherwise the cost of T(v,p) is considered as inf;
	every vertex v is associated with at most 2^|Gamma| trees*/
	std::unordered_map<int, std::unordered_map<int, graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_tree_node>> trees; // <v, <p, T(v,p)>>
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		int v = it->first; // a vertex
		trees[v] = {};
		int group_set_ID_v = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_vertex_group_set_ID(v, group_graph, cumpulsory_group_vertices); /*time complexity: O(|Gamma|)*/
		for (int p = 1; p <= group_set_ID_v; p++) { // p is non-empty; time complexity: O(2^|Gamma|)
			if ((p | group_set_ID_v) == group_set_ID_v) { // p represents a non-empty group set inside group_set_ID_v, including group_set_ID_v

				/*T(v,p)*/
				graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_tree_node node;
				node.cost = 0;
				node.type = 0;
				trees[v][p] = node;

				/*insert T(v,p) into Q_T*/
				graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
				x.v = v;
				x.p = p;
				x.priority_value = LB_procedure(v, p, 0, group_sets_ID_range, uncovered_groups, SPs_to_groups, W, W2);
				string handle_ID = to_string(v) + "_" + to_string(p);
				Q_T_handles[handle_ID] = Q_T.push(x);
				Q_T_priorities[handle_ID] = x.priority_value;

				if (x.v == 46 && x.p == 16) {
					cout << "here0" << endl;
					cout << "x.priority_value:" << x.priority_value << endl;
					getchar();
				}

				//cout << "initial Q push " + handle_ID << " priority: " << x.priority_value << endl;
			}
		}
	}



	/*D records the popped out optimal subtrees; String is "v_p" ("vertex_ID" + "_" + "group set ID") */
	std::unordered_set<string> D;

	//cout << "group_sets_ID_range:" << group_sets_ID_range << endl;

	/*Big while loop*/
	while (Q_T.size() > 0) { // at most 2^|Gamma|*V loops


		graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node top_node = Q_T.top();
		int v = top_node.v, X = top_node.p;
		double v_X_tree_cost = trees[v][X].cost, v_X_lb_in_Q = top_node.priority_value;
		Q_T.pop(); // O(2^|Gamma|*V*(|Gamma| + log V)) throught the loop, since Q_T contains at most 2^|Gamma|*V elements

		string handle_ID = to_string(v) + "_" + to_string(X);
		Q_T_handles.erase(handle_ID); // Q_T_handles only contains handles of elements in Q_T
		Q_T_priorities.erase(handle_ID);

		//cout << "X:" << X << endl;
		if (X == group_sets_ID_range) { // T(v,p) covers all groups 	
			graph_hash_of_mixed_weighted feasible_solu = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_build_tree(v, X, input_graph, trees);
			return feasible_solu;
		}

		//cout << "Q pop " + handle_ID << endl;
		D.insert(handle_ID); // optimal T(v,p) has been found

		std::unordered_set<int> X_uncovered_groups = uncovered_groups[X];

		/*build a feasible solution, report app ratio, and update best; O(|T||V| + |V|log|V|)*/
		double feasible_solu_cost = v_X_tree_cost;
		for (auto it = X_uncovered_groups.begin(); it != X_uncovered_groups.end(); it++) {
			feasible_solu_cost = feasible_solu_cost + SPs_to_groups[*it].second[v];
		}
		if (feasible_solu_cost < best_cost) {
			graph_hash_of_mixed_weighted feasible_solu = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_build_tree(v, X, input_graph, trees);
			for (auto it = X_uncovered_groups.begin(); it != X_uncovered_groups.end(); it++) {
				int g_id = *it;
				/*merge LWP(v to g_id) into feasible_solu; O(|V|)*/
				int v_start = v;
				int pre = SPs_to_groups[g_id].first[v_start];
				while (v_start != pre) {
					double w_v_pre = input_graph.hash_of_vectors[pre].vertex_weight;
					double ec = graph_hash_of_mixed_weighted_edge_weight(input_graph, v_start, pre);
					graph_hash_of_mixed_weighted_add_vertex(feasible_solu, pre, w_v_pre);
					graph_hash_of_mixed_weighted_add_edge(feasible_solu, v_start, pre, ec);
					v_start = pre;
					pre = SPs_to_groups[g_id].first[v_start];
				}
			}
			best_solu = graph_hash_of_mixed_weighted_MST_postprocessing_no_whole_graph(feasible_solu);
			best_cost = graph_hash_of_mixed_weighted_sum_of_ec(best_solu);
		}


		/*since an optimal solution has not been popped out, the optimal cost must be larger than or equal to v_X_tree_cost*/ 
		if (v_X_tree_cost > 0) {
			double ratio = best_cost / v_X_tree_cost;
			//cout << "ratio:" << ratio << " v_X_tree_cost:" << v_X_tree_cost << endl;
			if (ratio <= maximum_return_app_ratio + 1e-5) { // this feasible solution can be returned
				return best_solu;
			}
		}


		/*Lines 17-19 in PrunedDP++*/
		int X_slash = group_sets_ID_range - X; // p \cup X_slash equals group_sets_ID_range
		handle_ID = to_string(v) + "_" + to_string(X_slash);
		if (D.count(handle_ID) > 0) {

			double merged_tree_cost = v_X_tree_cost + trees[v][X_slash].cost;

			double cost_Tvp1_cup_p2 = inf;
			if (trees[v].count(group_sets_ID_range) > 0) {
				cost_Tvp1_cup_p2 = trees[v][group_sets_ID_range].cost;
			}
			if (merged_tree_cost < cost_Tvp1_cup_p2) { // update tree T(v,group_sets_ID_range)
				trees[v][group_sets_ID_range].cost = merged_tree_cost;
				trees[v][group_sets_ID_range].type = 2;
				trees[v][group_sets_ID_range].p1 = X;
				trees[v][group_sets_ID_range].p2 = X_slash;

				double lb = merged_tree_cost;
				//if (v_X_lb_in_Q > lb) {
				//	lb = v_X_lb_in_Q;
				//}

				if (merged_tree_cost > best_cost + error_safeguard) { // cannot have >= here, since best_solu may not be in Q; if >=, then best_solu may never be in Q
					continue;
				}
				if (merged_tree_cost < best_cost) {
					best_cost = merged_tree_cost;
					best_solu = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_build_tree(v, group_sets_ID_range, input_graph, trees);
				}

				handle_ID = to_string(v) + "_" + to_string(group_sets_ID_range);
				if (Q_T_priorities.count(handle_ID) == 0) {
					graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
					x.v = v;
					x.p = group_sets_ID_range;
					x.priority_value = lb;
					Q_T_handles[handle_ID] = Q_T.push(x);
					Q_T_priorities[handle_ID] = lb;
					//cout << "Q push " + handle_ID << " priority: " << Q_T_priorities[handle_ID] << endl;
				}
				else {
					if (lb < Q_T_priorities[handle_ID]) {
						graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
						x.v = v;
						x.p = group_sets_ID_range;
						x.priority_value = lb;
						Q_T.update(Q_T_handles[handle_ID], x); // O(1) for decrease key
						Q_T_priorities[handle_ID] = lb;
					}
				}
			}
			continue;
		}


		/* Optimal-Tree Decomposition Theorem (Theorem 1 in 2016 SIGMOD paper); optimal solution may not be found if below is / 2, why?*/
		if (v_X_tree_cost < best_cost / 1 + error_safeguard) {

			/*grow*/
			std::vector<int> v_adjs = input_graph.adj_v(v);
			for (auto it = v_adjs.begin(); it != v_adjs.end(); it++) {

				/*below: O(2^|Gamma|*E) in all loops, since each v has 2^|Gamma| times*/
				int u = *it;
				double cost_euv = graph_hash_of_mixed_weighted_edge_weight(input_graph, u, v);
				double grow_tree_cost = trees[v][X].cost + cost_euv;

				handle_ID = to_string(u) + "_" + to_string(X);
				//cout << "grow " << handle_ID <<  " grow_tree_cost " << grow_tree_cost << " best_cost " << best_cost << endl;
				if (D.count(handle_ID) > 0) {
					continue;
				}

				double T_up_cost;
				if (trees[u].count(X) == 0) {
					T_up_cost = inf;
				}
				else {
					T_up_cost = trees[u][X].cost;
				}
				if (grow_tree_cost < T_up_cost) {
					trees[u][X].cost = grow_tree_cost;
					trees[u][X].type = 1;
					trees[u][X].u = v;

					if (u == 56 && X == 24) {
						cout << "here5" << endl;
						cout << "grow u=" << u << " X=" << X << " grow_tree_cost=" << grow_tree_cost << endl;
						cout << "v:" << v << " trees[v][X].cost:" << trees[v][X].cost << endl;
						getchar();
					}

					double lb = LB_procedure(u, X, grow_tree_cost, group_sets_ID_range, uncovered_groups, SPs_to_groups, W, W2);


					vector<int> opt = { 46, 22, 56, 64, 0, 69, 1, 3, 42, 37, 79, 49, 90, 28 }; // vertices in optimal solution
					bool is_in = false;
					for (int ss = 0; ss < opt.size(); ss++) {
						if (u == opt[ss]) {
							is_in = 1;
							break;
						}
					}
					if (is_in && lb >3.615) { // lb is wrongly larger than optimal cost

						//if (u == 0 && X == 24) {
						//	cout << "grow u=" << u << " X=" << X << " grow_tree_cost=" << grow_tree_cost <<
						//		" lb=" << lb << " best_cost=" << best_cost << endl;
						//	cout << "v:" << v << " trees[v][X].cost:" << trees[v][X].cost << endl;
						//	cout << "group_set_ID " << X << " covered_groups:" << endl;
						//	print_unordered_set_int(covered_groups[X]);




						//	getchar();
						//}

						
					}


					if (v_X_lb_in_Q > lb) {
						lb = v_X_lb_in_Q;
					}
					//cout << "grow T_up_cost=" << T_up_cost << " lb=" << lb << " best_cost=" << best_cost << endl;


					




					if (lb > best_cost + error_safeguard) {
						//cout << "grow T_up_cost=" << T_up_cost << " lb=" << lb << " best_cost=" << best_cost << endl;
						continue;
					}

					if (Q_T_priorities.count(handle_ID) == 0) {
						graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
						x.v = u;
						x.p = X;
						x.priority_value = lb;
						Q_T_handles[handle_ID] = Q_T.push(x);
						Q_T_priorities[handle_ID] = lb;
						//cout << "Q push " + handle_ID << " priority: " << Q_T_priorities[handle_ID] << endl;

						
						if (x.v == 56 && x.p == 16) {
							cout << "here1" << endl;
							cout << "grow u=" << u << " X=" << X << " grow_tree_cost=" << grow_tree_cost 
								<< " x.priority_value=" << x.priority_value << endl;
							cout << "v:" << v << " trees[v][X].cost:" << trees[v][X].cost << endl;
							cout << "D.count(to_string(56) + - + to_string(8)): " << D.count(to_string(56) + "_" + to_string(8)) << endl;
							getchar();
						}

					}
					else {
						if (lb < Q_T_priorities[handle_ID]) {
							graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
							x.v = u;
							x.p = X;
							x.priority_value = lb;
							Q_T.update(Q_T_handles[handle_ID], x); // O(1) for decrease key
							Q_T_priorities[handle_ID] = lb;

							if (x.v == 56 && x.p == 16) {
								cout << "here2" << endl;
								getchar();
							}
						}
					}
				}
			}


			/*merge*/
			int p1 = X;
			std::unordered_set<int> non_overlapped_group_sets_p1 = non_overlapped_group_sets_IDs[p1];
			for (auto it = non_overlapped_group_sets_p1.begin(); it != non_overlapped_group_sets_p1.end(); it++) {

				int p2 = *it; // p2 is not overlapped with p1

				//if (v == 56 && p1+p2==24) {
				//	cout << "here8" << endl;
				//	cout << "v_X_lb_in_Q: " << v_X_lb_in_Q << endl;
				//	cout << "group_set_ID " << p1 << " covered_groups:" << endl;
				//	print_unordered_set_int(covered_groups[p1]);
				//	cout << "group_set_ID " << p2 << " covered_groups:" << endl;
				//	print_unordered_set_int(covered_groups[p2]);
				//	cout << "D.count(to_string(v) + - + to_string(p2)): " << D.count(to_string(v) + "_" + to_string(p2)) << endl;
				//	getchar();
				//}


				if (D.count(to_string(v) + "_" + to_string(p2)) > 0) { // only merge optimal T(v,p2), this is recorded in the pseudo code of 2016 paper

					int p1_cup_p2 = p1 + p2;
					handle_ID = to_string(v) + "_" + to_string(p1_cup_p2);

					if (D.count(handle_ID) > 0) {
						continue;
					}

					double merged_tree_cost = v_X_tree_cost + trees[v][p2].cost;

					if (v == 56) {
						cout << "here7" << endl;
						cout << "merge " << handle_ID << " merged_tree_cost " << merged_tree_cost << " best_cost " << best_cost << endl;
						getchar();
					}
					//cout << "merge " << handle_ID << " merged_tree_cost " << merged_tree_cost << " best_cost " << best_cost << endl;

					double cost_Tvp1_cup_p2;
					if (trees[v].count(p1_cup_p2) == 0) {
						cost_Tvp1_cup_p2 = inf;
					}
					else {
						cost_Tvp1_cup_p2 = trees[v][p1_cup_p2].cost;
					}

					if (merged_tree_cost < cost_Tvp1_cup_p2) { // O(3^|Gamma||V| comparisons in totel, see the DPBF paper)
						trees[v][p1_cup_p2].cost = merged_tree_cost;
						trees[v][p1_cup_p2].type = 2;
						trees[v][p1_cup_p2].p1 = p1;
						trees[v][p1_cup_p2].p2 = p2;

						if (v == 56 && p1_cup_p2 == 24) {
							cout << "here6" << endl;
							getchar();
						}

						/*Conditional Tree Merging Theorem (Theorem 2 in 2016 SIGMOD paper);
						note the following place of error_safeguard; this if condition should be met even if error_safeguard=inf*/
						if (merged_tree_cost <= best_cost * 2 / 3 + error_safeguard) { // error_safeguard is error

							double lb = LB_procedure(v, p1_cup_p2, merged_tree_cost, group_sets_ID_range, uncovered_groups, SPs_to_groups, W, W2);
							if (v_X_lb_in_Q > lb) {
								lb = v_X_lb_in_Q;
							}
							//cout << "merge cost_Tvp1_cup_p2=" << cost_Tvp1_cup_p2 << " lb_slash=" << lb_slash << " lb=" << lb << " best_cost=" << best_cost << endl;

							if (lb > best_cost + error_safeguard) {
								continue;
							}
							if (p1_cup_p2 == group_sets_ID_range && merged_tree_cost < best_cost) {
								best_cost = merged_tree_cost;
								best_solu = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_build_tree(v, group_sets_ID_range, input_graph, trees);
							}

							if (Q_T_priorities.count(handle_ID) == 0) {
								graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
								x.v = v;
								x.p = p1_cup_p2;
								x.priority_value = lb;
								Q_T_handles[handle_ID] = Q_T.push(x);
								Q_T_priorities[handle_ID] = lb;
								//cout << "Q push " + handle_ID << " priority: " << Q_T_priorities[handle_ID] << endl;

								if (x.v == 56 && x.p == 24) {
									cout << "here3" << endl;
									getchar();
								}


							}
							else {
								if (lb < Q_T_priorities[handle_ID]) {
									graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted_min_node x;
									x.v = v;
									x.p = p1_cup_p2;
									x.priority_value = lb;
									Q_T.update(Q_T_handles[handle_ID], x); // O(1) for decrease key
									Q_T_priorities[handle_ID] = lb;
									//cout << "Q update " + handle_ID << " priority: " << Q_T_priorities[handle_ID] << endl;

									if (x.v == 56 && x.p == 24) {
										cout << "here4" << endl;
										getchar();
									}
								}
							}
						}
					}
				}
			}
		}
	}


	std::cout << "graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted did not find a feasible solution!" << std::endl;
	graph_hash_of_mixed_weighted_save_for_GSTP("PrunedDPPlusPlus_irregular.txt", input_graph, group_graph, cumpulsory_group_vertices, 1);
	getchar();
	exit(1);

}





/*debug codes*/

#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_DPBF_edge_weighted.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_generate_random_connected_graph.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_generate_random_groups_of_vertices.h>

#pragma region
bool this_is_a_feasible_solution(graph_hash_of_mixed_weighted& solu, graph_hash_of_mixed_weighted& group_graph,
	std::unordered_set<int>& group_vertices) {

	/*time complexity O(|V_solu|+|E_solu|)*/
	if (graph_hash_of_mixed_weighted_connected_components(solu).size() != 1) { // it's not connected
		cout << "this_is_a_feasible_solution: solu is disconnected!" << endl;
		return false;
	}

	for (auto it = group_vertices.begin(); it != group_vertices.end(); it++) {
		int g = *it;
		bool covered = false;
		for (auto it2 = solu.hash_of_vectors.begin(); it2 != solu.hash_of_vectors.end(); it2++) {
			int v = it2->first;
			if (graph_hash_of_mixed_weighted_contain_edge(group_graph, v, g)) {
				covered = true;
				break;
			}
		}
		if (covered == false) {
			cout << "this_is_a_feasible_solution: a group is not covered!" << endl;
			return false;
		}
	}

	return true;

}
#pragma endregion this_is_a_feasible_solution

#pragma region
void simple_iterative_tests_compare_DPBF_PrunedDPPlusPlus() {

	/*parameters*/
	int iteration_times = 3;
	int V = 10000, E = 100000, G = 6, g_size_min = 10, g_size_max = 100, precision = 3;
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

		/*graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted*/
		if (1) {
			auto begin = std::chrono::high_resolution_clock::now();
			graph_hash_of_mixed_weighted solu = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted(instance_graph, generated_group_graph, generated_group_vertices, 1);
			auto end = std::chrono::high_resolution_clock::now();
			double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
			time_PrunedDPPlusPlus_avg = time_PrunedDPPlusPlus_avg + (double)runningtime / iteration_times;

			//graph_hash_of_mixed_weighted_print(solu);

			double cost = graph_hash_of_mixed_weighted_sum_of_ec(solu);
			solution_cost_PrunedDPPlusPlus_sum = solution_cost_PrunedDPPlusPlus_sum + cost;

			if (!this_is_a_feasible_solution(solu, generated_group_graph, generated_group_vertices)) {
				cout << "Error: graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted is not feasible!" << endl;
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
#pragma endregion simple_iterative_tests_compare_DPBF_PrunedDPPlusPlus

#pragma region
void PrunedDPPlusPlus_dedug_irregular_case() {

	double lambda = 1;
	std::unordered_set<int> generated_group_vertices;
	graph_hash_of_mixed_weighted instance_graph, generated_group_graph;
	graph_hash_of_mixed_weighted_read_for_GSTP("PrunedDPPlusPlus_irregular.txt", instance_graph,
		generated_group_graph, generated_group_vertices, lambda);

	if (1) {
		auto begin = std::chrono::high_resolution_clock::now();
		graph_hash_of_mixed_weighted solu = graph_hash_of_mixed_weighted_DPBF_edge_weighted(instance_graph, generated_group_graph, generated_group_vertices);
		auto end = std::chrono::high_resolution_clock::now();
		double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

		double cost = graph_hash_of_mixed_weighted_sum_of_ec(solu);

		graph_hash_of_mixed_weighted_print(solu);
		cout << "cost=" << cost << endl;
	}

	/*graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted*/
	if (1) {
		auto begin = std::chrono::high_resolution_clock::now();
		graph_hash_of_mixed_weighted solu = graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted(instance_graph, generated_group_graph, generated_group_vertices, 1);
		auto end = std::chrono::high_resolution_clock::now();
		double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

		//graph_hash_of_mixed_weighted_print(solu);

		double cost = graph_hash_of_mixed_weighted_sum_of_ec(solu);

		if (!this_is_a_feasible_solution(solu, generated_group_graph, generated_group_vertices)) {
			cout << "Error: graph_hash_of_mixed_weighted_PrunedDPPlusPlus_edge_weighted is not feasible!" << endl;
			graph_hash_of_mixed_weighted_print(solu);
			exit(1);
		}
	}
}
#pragma endregion PrunedDPPlusPlus_dedug_irregular_case
