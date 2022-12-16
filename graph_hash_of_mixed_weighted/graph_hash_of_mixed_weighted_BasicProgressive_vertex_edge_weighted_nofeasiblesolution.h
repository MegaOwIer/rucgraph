#pragma once

/*
time complexity: O(  )

cost function: (1 - lambda) * vertex_weights + lambda * edge_weights
*/
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <boost/heap/fibonacci_heap.hpp> 

#pragma region
graph_hash_of_mixed_weighted graph_hash_of_mixed_weighted_BasicProgressive_G2(graph_hash_of_mixed_weighted& input_graph, double lambda) {

	/*time complexity: O(|V|+|E|);

	it outputs a new graph, does not change input_graph;

	edge weights in G_2 have been embedded with regulated vertex weights;

	vertex weights in G_2 are the same with those in input_graph;
	*/

	graph_hash_of_mixed_weighted G_2 = graph_hash_of_mixed_weighted_copy_graph(input_graph);

	for (auto it1 = G_2.hash_of_vectors.begin(); it1 != G_2.hash_of_vectors.end(); it1++) {
		int i = it1->first;

		auto search = G_2.hash_of_hashs.find(i);
		if (search != G_2.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int j = it2->first;
				if (i < j) { // edge (i,j)
					double c_ij = it2->second;
					double w_i = it1->second.vertex_weight;
					double w_j = G_2.hash_of_vectors[j].vertex_weight;
					double new_c_ij = lambda * c_ij + (1 - lambda) * (w_i / 2 + w_j / 2);
					graph_hash_of_mixed_weighted_add_edge(G_2, i, j, new_c_ij); // embed nw into ec
				}
			}
		}
		else {
			for (auto it2 = it1->second.adj_vertices.begin(); it2 != it1->second.adj_vertices.end(); it2++) {
				int j = it2->first;
				if (i < j) { // edge (i,j)
					double c_ij = it2->second;
					double w_i = it1->second.vertex_weight;
					double w_j = G_2.hash_of_vectors[j].vertex_weight;
					double new_c_ij = lambda * c_ij + (1 - lambda) * (w_i / 2 + w_j / 2);
					graph_hash_of_mixed_weighted_add_edge(G_2, i, j, new_c_ij); // embed nw into ec
				}
			}
		}
	}

	return G_2;

}

std::unordered_map<int, double> graph_hash_of_mixed_weighted_BasicProgressive_find_LWPs_to_g
(graph_hash_of_mixed_weighted& group_graph, graph_hash_of_mixed_weighted& G_2, int g_vertex, double lambda) {

	/*time complexity: O(|E|+|V|log|V|)*/

	/*add dummy vertex and edges; time complexity: O(|V|)*/
	graph_hash_of_mixed_weighted_add_vertex(G_2, g_vertex, 0); // add dummy vertex
	auto search = group_graph.hash_of_hashs.find(g_vertex);
	if (search != group_graph.hash_of_hashs.end()) {
		for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
			int vertex = it2->first; // vertex is in group g_vertex

			double nw = G_2.hash_of_vectors[vertex].vertex_weight; // vertex weights in G_2 are the same with those in input_graph
			double dummy_ec = (1 - lambda) * nw / 2; // original dummy edge weight is 0, this is the weight embedded with vertex weights
			graph_hash_of_mixed_weighted_add_edge(G_2, g_vertex, vertex, dummy_ec); // add dummy edge

		}
	}
	else {
		auto search2 = group_graph.hash_of_vectors.find(g_vertex); // if v is not in g, error is here
		for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
			int vertex = it2->first; // vertex is in group g_vertex

			double nw = G_2.hash_of_vectors[vertex].vertex_weight; // vertex weights in G_2 are the same with those in input_graph
			double dummy_ec = (1 - lambda) * nw / 2; // original dummy edge weight is 0, this is the weight embedded with vertex weights
			graph_hash_of_mixed_weighted_add_edge(G_2, g_vertex, vertex, dummy_ec); // add dummy edge

		}
	}

	/*time complexity: O(|E|+|V|log|V|)*/
	std::unordered_map<int, double> distances;
	std::unordered_map<int, int> predecessors;
	graph_hash_of_mixed_weighted_shortest_paths_source_to_all(G_2, g_vertex, distances, predecessors); /*time complexity: O(|V|+|E|)*/
	graph_hash_of_mixed_weighted_remove_vertex(G_2, g_vertex);  // all dummy vertex and edges are removed; time complexity: O(|V|)
	distances.erase(g_vertex);
	predecessors.erase(g_vertex);
	for (auto it = distances.begin(); it != distances.end(); it++) {
		int v = it->first;
		double nw = G_2.hash_of_vectors[v].vertex_weight; // vertex weights in G_2 are the same with those in input_graph
		it->second = it->second + (1 - lambda) * nw / 2; // update the distance to be the regulated weight of the LWP from v to g_vertex
	}
	for (auto it = predecessors.begin(); it != predecessors.end(); it++) {
		int pre = it->second;
		if (pre == g_vertex) {
			it->second = it->first; // since g_vertex is not in predecessors, it->second points to it->first, i.e., the path ends at it->first.
		}
	}

	return distances;

}

std::unordered_map<int, std::unordered_map<int, double>> graph_hash_of_mixed_weighted_BasicProgressive_find_LWPs
(graph_hash_of_mixed_weighted& input_graph, graph_hash_of_mixed_weighted& group_graph, std::unordered_set<int>& cumpulsory_group_vertices, double lambda) {

	/*return {g_ID, { distances, predecessors }}*/

	/*time complexity: O(|V|+|E|)*/
	graph_hash_of_mixed_weighted G_2 = graph_hash_of_mixed_weighted_BasicProgressive_G2(input_graph, lambda);

	/*time complexity: O(|T||E|+|T||V|log|V|)*/
	std::unordered_map<int, std::unordered_map<int, double>> LWPs_to_groups;
	for (auto it = cumpulsory_group_vertices.begin(); it != cumpulsory_group_vertices.end(); it++) {
		int g_vertex = *it;
		LWPs_to_groups[g_vertex] = graph_hash_of_mixed_weighted_BasicProgressive_find_LWPs_to_g(group_graph, G_2, g_vertex, lambda);
		//cout << "LWPs_to_groups[" << g_vertex << "]:" << endl;
		//print_unordered_map_int_double(LWPs_to_groups[g_vertex].first);
		//print_unordered_map_int_int(LWPs_to_groups[g_vertex].second);
	}

	return LWPs_to_groups;
}

#pragma endregion graph_hash_of_mixed_weighted_BasicProgressive_find_LWPs

#pragma region
std::unordered_map<int, std::unordered_set<int>> graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_non_overlapped_group_sets(int group_sets_ID_range) {

	/*this function calculate the non-empty and non_overlapped_group_sets_IDs of each non-empty group_set ID (including itself);

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
#pragma endregion graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_non_overlapped_group_sets

#pragma region
struct graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_min_node {
	int v;
	int p; // group_set_ID
	double priority_value; // cost of T(v,p)
};
bool operator<(graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_min_node const& x, graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_min_node const& y) {
	return x.priority_value > y.priority_value; // < is the max-heap; > is the mean heap; PriorityQueue is expected to be a max-heap of integer values
}
typedef typename boost::heap::fibonacci_heap<graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_min_node>::handle_type handle_graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_min_node;
#pragma endregion graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted priority queue

#pragma region
class graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_tree_node
{
	/*this is like the tree T(v,p) in the DPBF paper*/

public:

	int type; // =0: this is the single vertex v; =1: this tree is built by grown; =2: built by merge

	double cost; // cost of this tree T(v,p); cost function: (1 - lambda) * vertex_weights + lambda * edge_weights

	int u; // if this tree is built by grown, then it's built by growing edge (v,u);

	int p1, p2; // if this tree is built by merge, then it's built by merge T(v,p1) and T(v,p2);

};
#pragma endregion graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_tree_node

#pragma region
int graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_vertex_group_set_ID(int vertex, graph_hash_of_mixed_weighted& group_graph,
	std::unordered_set<int>& cumpulsory_group_vertices) {

	/*time complexity: O(|Gamma|);

	this function returns the maximum group set ID for a single vertex*/

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
#pragma endregion graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_vertex_group_set_ID

#pragma region
std::unordered_map<int, std::unordered_set<int>> graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_uncovered_groups(int group_sets_ID_range,
	std::unordered_set<int>& cumpulsory_group_vertices) {

	/*time complexity: O(|Gamma|*2^|Gamma|);

	for each p \in [1,group_sets_ID_range], this function calculate the groups that have not been coverred by p*/

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
#pragma endregion graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_uncovered_groups

#pragma region
graph_hash_of_mixed_weighted graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_build_tree(int root_v, int root_p, graph_hash_of_mixed_weighted& input_graph,
	std::unordered_map<int, std::unordered_map<int, graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_tree_node>>& trees) {

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
#pragma endregion graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_build_tree

#pragma region
double graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_treeCost(graph_hash_of_mixed_weighted& theta, double lambda) {

	/*time complexity: O(|V|)*/

	double lambda_cost = 0;

	for (auto it1 = theta.hash_of_vectors.begin(); it1 != theta.hash_of_vectors.end(); it1++) {
		int i = it1->first;
		double w_i = it1->second.vertex_weight;
		lambda_cost = lambda_cost + (1 - lambda) * w_i;
		auto search = theta.hash_of_hashs.find(i);
		if (search != theta.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int j = it2->first;
				if (i < j) { // edge (i,j)
					double c_ij = it2->second;
					lambda_cost = lambda_cost + lambda * c_ij;
				}
			}
		}
		else {
			for (auto it2 = it1->second.adj_vertices.begin(); it2 != it1->second.adj_vertices.end(); it2++) {
				int j = it2->first;
				if (i < j) { // edge (i,j)
					double c_ij = it2->second;
					lambda_cost = lambda_cost + lambda * c_ij;
				}
			}
		}

	}

	return lambda_cost;

}
#pragma endregion graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_treeCost

graph_hash_of_mixed_weighted graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted(graph_hash_of_mixed_weighted& input_graph, graph_hash_of_mixed_weighted& group_graph,
	std::unordered_set<int>& cumpulsory_group_vertices, double lambda) {

	/*cost function: (1 - lambda) * vertex_weights + lambda * edge_weights */


	if (cumpulsory_group_vertices.size() >= 20) {
		std::cout << "cumpulsory_group_vertices.size() is too large for graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted!" << std::endl;
		exit(1);
	}


	/*finding lowest-weighted paths from groups to vertices;
	time complexity: O(|T||E|+|T||V|log|V|);
	return {g_ID, { distances, predecessors }} */
	std::unordered_map<int, std::unordered_map<int, double>> LWPs_to_groups =
		graph_hash_of_mixed_weighted_BasicProgressive_find_LWPs(input_graph, group_graph, cumpulsory_group_vertices, lambda);

	/*initialize Q*/
	boost::heap::fibonacci_heap<graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_min_node> Q_T; // min queues of trees
	std::unordered_map<string, handle_graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_min_node>
		Q_T_handles; // key is string "v_p" ("vertex_ID" + "_" + "group set ID"); Q_T_handles only contain keys that are in Q_T

	double inf = std::numeric_limits<double>::max(); // cost of empty tree is inf

	/* this is the cost of the best found solution yet */
	double best_cost = inf;

	/*initialize non_overlapped_group_sets; time complexity: O(4^|Gamma|);
	Group		G1	G0	group_set_ID
				0   0   0
				0	1	1
				1	0	2
				1	1	3*/
	int group_sets_ID_range = pow(2, cumpulsory_group_vertices.size()) - 1; // the number of group sets: 2^|Gamma|, including empty set;   |Gamma| should be smaller than 31 due to precision
	std::unordered_map<int, std::unordered_set<int>> non_overlapped_group_sets_IDs =
		graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_non_overlapped_group_sets(group_sets_ID_range);  // time complexity: O(4^|Gamma|)


	/*initialize uncovered_groups;  <p, <uncovered_groups>>;  time complexity: O(|Gamma|*2^|Gamma|)*/
	std::unordered_map<int, std::unordered_set<int>> uncovered_groups = graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_uncovered_groups(
		group_sets_ID_range, cumpulsory_group_vertices);






	/*initialize trees with vertices;
	time complexity: O(2^|Gamma|*|V|);
	every vertex v is associated with T(v,p) only when v covers p, otherwise the cost of T(v,p) is considered as inf;
	every vertex v is associated with at most 2^|Gamma| trees;
	*/
	std::unordered_map<int, std::unordered_map<int, graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_tree_node>> trees; // <v, <p, T(v,p)>>
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		int v = it->first; // a vertex
		double cost_v = (1 - lambda) * it->second.vertex_weight; // lambda is embedded
		trees[v] = {};
		int group_set_ID_v = graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_vertex_group_set_ID(v, group_graph, cumpulsory_group_vertices); /*time complexity: O(|Gamma|)*/
		for (int p = 1; p <= group_set_ID_v; p++) { // p is non-empty; time complexity: O(2^|Gamma|)
			if ((p | group_set_ID_v) == group_set_ID_v) { // p represents a non-empty group set inside group_set_ID_v, including group_set_ID_v

				/*T(v,p)*/
				graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_tree_node node;
				node.cost = cost_v;
				node.type = 0;
				trees[v][p] = node;

				/*insert T(v,p) into Q_T*/
				graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_min_node x;
				x.v = v;
				x.p = p;
				x.priority_value = cost_v;
				string handle_ID = to_string(v) + "_" + to_string(p);
				Q_T_handles[handle_ID] = Q_T.push(x);
			}
		}
	}








	/*D records the popped out optimal subtrees; String is "v_p" ("vertex_ID" + "_" + "group set ID") */
	std::unordered_set<string> D;



	/*Big while loop*/
	while (Q_T.size() > 0) { // at most 2^|Gamma|*V loops


		graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_min_node top_node = Q_T.top();
		int v = top_node.v, p = top_node.p;
		double cost_v = (1 - lambda) * input_graph.hash_of_vectors[v].vertex_weight; // lambda is embedded

		Q_T.pop(); // O(2^|Gamma|*V*(|Gamma| + log V)) throught the loop, since Q_T contains at most 2^|Gamma|*V elements
		string handle_ID = to_string(v) + "_" + to_string(p);
		Q_T_handles.erase(handle_ID); // Q_T_handles only contains handles of elements in Q_T

		D.insert(handle_ID); // optimal T(v,p) has been found
		//cout << "D.insert " << handle_ID << endl;

		if (p == group_sets_ID_range) { // T(v,p) covers all groups 	
			graph_hash_of_mixed_weighted feasible_solu = graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_build_tree(v, p, input_graph, trees); // time complexity: O(|V|)
			return feasible_solu;
		}

		/*build a feasible solution, report app ratio, and update best; O(|T||V| + |V|log|V|)*/
		double feasible_solu_cost = top_node.priority_value;
		for (auto it = uncovered_groups[p].begin(); it != uncovered_groups[p].end(); it++) {
			int g_id = *it;
			/*merge LWP(v to g_id) into feasible_solu; O(|V|)*/
			int v_start = v;
			feasible_solu_cost = feasible_solu_cost + LWPs_to_groups[g_id][v_start] - cost_v;
		}
		if (feasible_solu_cost < best_cost) { // best_cost is also updated in merge!
			best_cost = feasible_solu_cost;
			//double T_v_p_cost = trees[v][p].cost;
			//double ratio = feasible_solu_cost / T_v_p_cost; // best.app_ratio may not always decrease
		}


		/*grow*/
		std::vector<int> v_adjs = input_graph.adj_v(v);
		for (auto it = v_adjs.begin(); it != v_adjs.end(); it++) {

			/*below: O(2^|Gamma|*E) in all loops, since each v has 2^|Gamma| times*/
			int u = *it;
			double cost_u = (1 - lambda) * input_graph.hash_of_vectors[u].vertex_weight; // lambda is embedded
			double cost_euv = lambda * graph_hash_of_mixed_weighted_edge_weight(input_graph, u, v); // lambda is embedded
			double grow_tree_cost = trees[v][p].cost + cost_euv + cost_u;

			handle_ID = to_string(u) + "_" + to_string(p);
			//cout << "grow " << handle_ID << endl;
			//cout << "grow_tree_cost " << grow_tree_cost << endl;
			//cout << "best_cost " << best_cost << endl;
			if (D.count(handle_ID) > 0 || grow_tree_cost > best_cost + 1e-10) { // 1e-10 is error
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
				graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_min_node x;
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
			if (D.count(handle_ID) > 0) { // only merge optimal T(v,p2)

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
				double merged_tree_cost = cost_Tvp1 + cost_Tvp2 - cost_v;

				if (merged_tree_cost > best_cost + 1e-10) { // 1e-10 is error
					continue;
				}

				if (p1_cup_p2 == group_sets_ID_range) {
					best_cost = merged_tree_cost;
				}

				if (merged_tree_cost < cost_Tvp1_cup_p2) { // O(3^|Gamma||V| comparisons in totel, see the DPBF paper)

					/*update T(v,p1_cup_p2) by merge T(v,p1) with T(v,v2)*/
					trees[v][p1_cup_p2].cost = merged_tree_cost;
					trees[v][p1_cup_p2].type = 2;
					trees[v][p1_cup_p2].p1 = p1;
					trees[v][p1_cup_p2].p2 = p2;

					/*update T(v,p1_cup_p2) in Q_T*/
					graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted_min_node x;
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

	}



	std::cout << "graph_hash_of_mixed_weighted_BasicProgressive_vertex_edge_weighted did not find a feasible solution!" << std::endl;
	getchar();
	exit(1);

}

