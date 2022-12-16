#pragma once



#include <unordered_map>
#include <boost/heap/fibonacci_heap.hpp> 

/*
this is the find the lowest weight path, i.e., the sum of vertex and edge weights in this path is minimized,
from a source vertex to the other vertices, using a modified Dijkstra¡¯s algorithm.

In these codes, we find shortest paths in a graph with new edge weights: c_new(v_1,v_2) = c(v_1,v_2) + w(v_1)/2 + w(v_2)/2

time complexity: O(|E|+|V|log|V|)
*/

struct node_for_lwp {
	int index;
	double priority_value;
}; // define the node in the queue
bool operator<(node_for_lwp const& x, node_for_lwp const& y) {
	return x.priority_value > y.priority_value; // < is the max-heap; > is the min heap
}
typedef typename boost::heap::fibonacci_heap<node_for_lwp>::handle_type handle_t_for_lwp;



void graph_hash_of_mixed_weighted_lowest_weight_paths_source_to_all(graph_hash_of_mixed_weighted& input_graph, int source,
	std::unordered_map<int, double>& path_weights, std::unordered_map<int, int>& predecessors) {

	/*Dijkstra¡¯s shortest path algorithm: https://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-greedy-algo-7/
	time complexity: O(|E|+|V|log|V|);*/

	double inf = std::numeric_limits<double>::max();

	node_for_lwp node;
	boost::heap::fibonacci_heap<node_for_lwp> Q;
	std::unordered_map<int, double> Q_keys;
	std::unordered_map<int, handle_t_for_lwp> Q_handles;

	/*time complexity: O(|V|)*/
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		int v = it->first;
		node.index = v;
		node.priority_value = inf;
		Q_keys[v] = inf;
		Q_handles[v] = Q.push(node);
	}

	/*initialize the source*/
	Q_keys[source] = 0;
	node.index = source;
	node.priority_value = Q_keys[source];
	Q.update(Q_handles[source], node);
	//std::cout << "Q.update " << node.index << " key:" << node.priority_value << '\n';
	predecessors[source] = source;

	/*time complexity: O(|E|+|V|log|V|) based on fibonacci_heap, not on pairing_heap, which is
	O((|E|+|V|)log|V|);
	there are at most |V| elements in Q; Q only contains vertices not touched*/
	while (Q.size() > 0) {

		node = Q.top();
		Q.pop();
		int top_v = node.index;
		double top_key = Q_keys[top_v];

		path_weights[top_v] = top_key; // path_weights here are NOT the sum of vertex and edge weights from source to top_v

		double half_top_v_weight = input_graph.hash_of_vectors[top_v].vertex_weight / 2;

		auto search = input_graph.hash_of_hashs.find(top_v);
		if (search != input_graph.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int adj_v = it2->first;
				double ec = it2->second;
				double half_adj_v_weight = input_graph.hash_of_vectors[adj_v].vertex_weight / 2;
				double G_2_new_ec = ec + half_top_v_weight + half_adj_v_weight;
				if (path_weights.count(adj_v) == 0) {
					if (Q_keys[adj_v] > path_weights[top_v] + G_2_new_ec) { // needs to update key
						Q_keys[adj_v] = path_weights[top_v] + G_2_new_ec;
						node.index = adj_v;
						node.priority_value = Q_keys[adj_v];
						Q.update(Q_handles[adj_v], node);
						predecessors[adj_v] = top_v;
						//std::cout << "Q.update " << node.index << " key:" << node.priority_value << '\n';
					}
				}
			}
		}
		else {
			auto search2 = input_graph.hash_of_vectors.find(top_v);
			for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
				int adj_v = it2->first;
				double ec = it2->second;
				double half_adj_v_weight = input_graph.hash_of_vectors[adj_v].vertex_weight / 2;
				double G_2_new_ec = ec + half_top_v_weight + half_adj_v_weight;
				if (path_weights.count(adj_v) == 0) {
					if (Q_keys[adj_v] > path_weights[top_v] + G_2_new_ec) { // needs to update key
						Q_keys[adj_v] = path_weights[top_v] + G_2_new_ec;
						node.index = adj_v;
						node.priority_value = Q_keys[adj_v];
						Q.update(Q_handles[adj_v], node);
						predecessors[adj_v] = top_v;
						//std::cout << "Q.update " << node.index << " key:" << node.priority_value << '\n';
					}
				}
			}
		}
	}

	/*update path_weights*/
	double half_source_weight = input_graph.hash_of_vectors[source].vertex_weight / 2;
	for (auto it = path_weights.begin(); it != path_weights.end(); it++) {
		int terminal = it->first;
		double half_terminal_weight = input_graph.hash_of_vectors[terminal].vertex_weight / 2;
		it->second = it->second + half_source_weight + half_terminal_weight; // path_weights here are the sum of vertex and edge weights from source to top_v (including weights of source and terminal)
	}

}



void graph_hash_of_mixed_weighted_lowest_weight_paths_source_to_all_lambda(graph_hash_of_mixed_weighted& input_graph, int source,
	std::unordered_map<int, double>& path_weights, std::unordered_map<int, int>& predecessors, double lambda) {

	/*Dijkstra¡¯s shortest path algorithm: https://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-greedy-algo-7/
	time complexity: O(|E|+|V|log|V|);

	cost function: (1 - lambda) * vertex_weights + lambda * edge_weights*/

	double inf = std::numeric_limits<double>::max();

	node_for_lwp node;
	boost::heap::fibonacci_heap<node_for_lwp> Q;
	std::unordered_map<int, double> Q_keys;
	std::unordered_map<int, handle_t_for_lwp> Q_handles;

	/*time complexity: O(|V|)*/
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		int v = it->first;
		node.index = v;
		node.priority_value = inf;
		Q_keys[v] = inf;
		Q_handles[v] = Q.push(node);
	}

	/*initialize the source*/
	Q_keys[source] = 0;
	node.index = source;
	node.priority_value = Q_keys[source];
	Q.update(Q_handles[source], node);
	//std::cout << "Q.update " << node.index << " key:" << node.priority_value << '\n';
	predecessors[source] = source;

	/*time complexity: O(|E|+|V|log|V|) based on fibonacci_heap, not on pairing_heap, which is
	O((|E|+|V|)log|V|);
	there are at most |V| elements in Q; Q only contains vertices not touched*/
	while (Q.size() > 0) {

		node = Q.top();
		Q.pop();
		int top_v = node.index;
		double top_key = Q_keys[top_v];

		path_weights[top_v] = top_key; // path_weights here are NOT the sum of vertex and edge weights from source to top_v

		double half_top_v_weight = (1 - lambda) * input_graph.hash_of_vectors[top_v].vertex_weight / 2;

		auto search = input_graph.hash_of_hashs.find(top_v);
		if (search != input_graph.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int adj_v = it2->first;
				double ec = lambda * it2->second;
				double half_adj_v_weight = (1 - lambda) * input_graph.hash_of_vectors[adj_v].vertex_weight / 2;
				double G_2_new_ec = ec + half_top_v_weight + half_adj_v_weight; // edges weights are updated (regulated) edges weights c(i,j)+w(i)/2+w(j)/2
				if (path_weights.count(adj_v) == 0) {
					if (Q_keys[adj_v] > path_weights[top_v] + G_2_new_ec) { // needs to update key
						Q_keys[adj_v] = path_weights[top_v] + G_2_new_ec;
						node.index = adj_v;
						node.priority_value = Q_keys[adj_v];
						Q.update(Q_handles[adj_v], node);
						predecessors[adj_v] = top_v;
						//std::cout << "Q.update " << node.index << " key:" << node.priority_value << '\n';
					}
				}
			}
		}
		else {
			auto search2 = input_graph.hash_of_vectors.find(top_v);
			for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
				int adj_v = it2->first;
				double ec = lambda * it2->second;
				double half_adj_v_weight = (1 - lambda) * input_graph.hash_of_vectors[adj_v].vertex_weight / 2;
				double G_2_new_ec = ec + half_top_v_weight + half_adj_v_weight; // edges weights are updated (regulated) edges weights c(i,j)+w(i)/2+w(j)/2
				if (path_weights.count(adj_v) == 0) {
					if (Q_keys[adj_v] > path_weights[top_v] + G_2_new_ec) { // needs to update key
						Q_keys[adj_v] = path_weights[top_v] + G_2_new_ec;
						node.index = adj_v;
						node.priority_value = Q_keys[adj_v];
						Q.update(Q_handles[adj_v], node);
						predecessors[adj_v] = top_v;
						//std::cout << "Q.update " << node.index << " key:" << node.priority_value << '\n';
					}
				}
			}
		}
	}

	/*update path_weights; in the above queueing process, just imagine that edges weights are updated edges weights (regulated) c(i,j)+w(i)/2+w(j)/2*/
	double half_source_weight = input_graph.hash_of_vectors[source].vertex_weight / 2;
	for (auto it = path_weights.begin(); it != path_weights.end(); it++) {
		int terminal = it->first;
		double half_terminal_weight = (1 - lambda) * input_graph.hash_of_vectors[terminal].vertex_weight / 2;
		it->second = it->second + half_source_weight + half_terminal_weight; // path_weights here are the sum of vertex and edge weights from source to top_v (including weights of source and terminal)
	}

}
