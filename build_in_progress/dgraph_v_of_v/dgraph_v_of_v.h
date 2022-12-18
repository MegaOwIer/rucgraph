#pragma once
#include<vector>
#include<iostream>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_binary_operations.h>


/*
the paired adjacency lists of node v;

if v1 is in INs of v2, then v2 must be in OUTs of v1, and both records should be associated with the same weight.
*/
template <typename weight_type> // weight_type may be int, float, double
class dgraph_v_of_v_adjs {
public:
	std::vector<std::pair<int, weight_type>> INs; // set of adjacency vertices that can reach v (with arc/edge weights)
	std::vector<std::pair<int, weight_type>> OUTs; // set of adjacency vertices that can be reached by v (with arc/edge weights)
};


/*a directed graph with edge weights, but not node weights*/
template <typename weight_type>
class dgraph_v_of_v {
public:
	std::vector<dgraph_v_of_v_adjs<weight_type>> adj_lists; // adj_lists.size()==n is the number of vertices, IDs from 0 to n-1

	/*constructors*/

	dgraph_v_of_v() {}

	dgraph_v_of_v(int n) {
		adj_lists.resize(n); // initialize n vertices
	}




	/*class member functions*/
	inline void add_edge(int, int, weight_type); // this function can change edge weights
	inline void remove_edge(int, int);
	inline weight_type edge_weight(int, int);
	inline bool contain_edge(int, int); // whether there is an edge
	inline long long int edge_number(); // the total number of edges
};





/*class member functions*/

template <typename weight_type>
void dgraph_v_of_v<weight_type>::add_edge(int v1, int v2, weight_type weight) { 

	/*
	edge direction: v1 to v2;
	this function adds v1 into INs of v2, and also adds v2 into OUTs of v1;
	this function can change edge weights;
	*/

	graph_hash_of_mixed_weighted_binary_operations_insert(adj_lists[v1].OUTs, v2, weight);
	graph_hash_of_mixed_weighted_binary_operations_insert(adj_lists[v2].INs, v1, weight);
}


template <typename weight_type>
void dgraph_v_of_v<weight_type>::remove_edge(int v1, int v2) {

	/*edge direction: v1 to v2*/

	graph_hash_of_mixed_weighted_binary_operations_erase(adj_lists[v1].OUTs, v2);
	graph_hash_of_mixed_weighted_binary_operations_erase(adj_lists[v2].INs, v1);
}

template <typename weight_type>
weight_type dgraph_v_of_v<weight_type>::edge_weight(int v1, int v2) {

	/*edge direction: v1 to v2*/

	return graph_hash_of_mixed_weighted_binary_operations_search_weight(adj_lists[v1].OUTs, v2);
}


template <typename weight_type>
bool dgraph_v_of_v<weight_type>::contain_edge(int v1, int v2) {

	/*this function checks two conditions: v1 is in INs of v2; v2 is in OUTs of v1*/

	int x = graph_hash_of_mixed_weighted_binary_operations_search(adj_lists[v1].OUTs, v2) +
		graph_hash_of_mixed_weighted_binary_operations_search(adj_lists[v2].INs, v1);

	if (x == 0) {
		return false;
	}
	else if (x == 2) {
		return true;
	}
	else {
		std::cout << "dgraph_v_of_v check_edge error!" << std::endl;
		if (graph_hash_of_mixed_weighted_binary_operations_search(adj_lists[v1].OUTs, v2)) {
			std::cout << v2 << " is in OUTs of" << v1 << std::endl;
		}
		else {
			std::cout << v2 << " is NOT in OUTs of" << v1 << std::endl;
		}
		if (graph_hash_of_mixed_weighted_binary_operations_search(adj_lists[v2].INs, v1)) {
			std::cout << v1 << " is in INs of" << v2 << std::endl;
		}
		else {
			std::cout << v1 << " is NOT in INs of" << v2 << std::endl;
		}
		exit(1);
	}
}


template <typename weight_type>
long long int dgraph_v_of_v<weight_type>::edge_number() {

	/*only check INs*/
	long long int num = 0;
	for (int i = adj_lists.size() - 1; i >= 0; i--) {
		num += adj_lists[i].INs.size();
	}

	return num;
}







/*
examples
----------------------
#include <build_in_progress/dgraph_v_of_v/dgraph_v_of_v.h>

int main()
{
	example_dgraph_v_of_v();
}
-----------------------
*/

void example_dgraph_v_of_v() {

	dgraph_v_of_v<double> g(10); // initialize a graph with 10 vertices

	g.add_edge(1, 5, 1.02); 
	g.add_edge(5, 1, 1.42);
	g.add_edge(5, 2, 122);
	g.remove_edge(5, 1);

	std::cout << g.edge_weight(1, 5) << std::endl;
	std::cout << g.contain_edge(1, 5) << std::endl;
	std::cout << g.contain_edge(5, 1) << std::endl;
	std::cout << g.edge_number() << std::endl;
	std::cout << g.adj_lists.size() << std::endl;
	

}