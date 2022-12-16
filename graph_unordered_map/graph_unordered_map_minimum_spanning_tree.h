#pragma once


#include <unordered_map>
#include <boost/heap/pairing_heap.hpp> 

/*pairing_heap is 1/6 fatser than fibonacci_heap*/

/*in a DBLP graph with 150k vertices and 500k edges,
the below code based on pairing_heap uses 0.5s; 
the below code based on fibonacci_heap uses 0.65s;
while the boost code uses 0.12s;
the slowness is likely due to the frequent use of unordered_map*/

struct node_for_mst {
	int index;
	double priority_value;
}; // define the node in the queue
bool operator<(node_for_mst const& x, node_for_mst const& y) {
	return x.priority_value > y.priority_value; // < is the max-heap; > is the min heap
}
typedef typename boost::heap::pairing_heap<node_for_mst>::handle_type handle_t_for_mst;

std::unordered_map<int, int> graph_unordered_map_minimum_spanning_tree(graph_unordered_map& input_graph) {

	/*Prim¡¯s Minimum Spanning Tree: https://www.geeksforgeeks.org/prims-minimum-spanning-tree-mst-greedy-algo-5/
	time complexity: O(|E|+|V|log|V|);*/

	std::unordered_map<int, int> mst; // edges: vertex, and its pair vertex in mst
	double inf = std::numeric_limits<double>::max();
	node_for_mst node;
	boost::heap::pairing_heap<node_for_mst> Q;
	std::unordered_map<int, double> Q_keys;
	std::unordered_map<int, handle_t_for_mst> Q_handles;

	/*time complexity: O(|V|)*/
	for (auto it = input_graph.begin(); it != input_graph.end(); it++) {
		int v = it->first;
		node.index = v;
		node.priority_value = inf;
		Q_keys[v] = inf;
		Q_handles[v] = Q.push(node);
	}

	/*select a start vertex*/
	auto start = input_graph.begin();
	int start_v = start->first;
	Q_keys[start_v] = 0;
	node.index = start_v;
	node.priority_value = Q_keys[start_v];
	Q.update(Q_handles[start_v], node);
	mst[start_v] = start_v; // edge (start_v, start_v) is in mst

	/*time complexity: O(|E|+|V|log|V|) based on fibonacci_heap, not on pairing_heap, which is
	O((|E|+|V|)log|V|);
	there are at most |V| elements in Q; Q only contains vertices not in mst*/
	while (Q.size() > 0) {

		node = Q.top();
		int top_v = node.index;
		if (node.priority_value == inf) {
			break; // this happens when the graph is disconnected
		}
		Q.pop();

		for (auto it = input_graph[top_v].second.begin(); it != input_graph[top_v].second.end(); it++) {
			int adj_v = it->first;
			double ec = it->second;
			if (mst.count(adj_v) > 0) { // adj_v is in mst
				if (ec == Q_keys[top_v]) {
					mst[top_v] = adj_v; // find an mst edge (top_v,adj_v); top_v is now in mst
				}
			}
			else {
				if (Q_keys[adj_v] > ec) { // needs to update key

					Q_keys[adj_v] = ec;
					node.index = adj_v;
					node.priority_value = Q_keys[adj_v];
					Q.update(Q_handles[adj_v], node);
				}
			}
		}

	}

	return mst; // this is not edges! this is predessors, and include {start_v, start_v}

}



void graph_unordered_map_minimum_spanning_tree_print(std::unordered_map<int, int>& mst) {

	std::cout << "graph_unordered_map_minimum_spanning_tree_print: " << std::endl;
	std::cout << "mst.size(): " << mst.size() << std::endl;

	for (auto it = mst.begin(); it != mst.end(); ++it) {
		std::cout << it->first << " - " << it->second << '\n';
	}

	std::cout << "graph_unordered_map_minimum_spanning_tree_print END" << '\n';

}




double graph_unordered_map_minimum_spanning_tree_cost(graph_unordered_map& input_graph,
	std::unordered_map<int, int>& mst) {

	double cost = 0;

	for (auto it = mst.begin(); it != mst.end(); ++it) {

		if (it->first != it->second) {

			double ec = input_graph[it->first].second[it->second];

			cost = cost + ec;
		}
	}

	return cost;

}
