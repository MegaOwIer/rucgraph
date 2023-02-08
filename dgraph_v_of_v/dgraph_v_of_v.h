#pragma once

/*
a directed graph with edge weights, but not node weights
*/
template <typename weight_type> // weight_type may be int, float, double
class dgraph_v_of_v {
    std::vector<std::vector<std::pair<int, weight_t>>> adj_lists;

public:
	/*
	INs.size() == OUTs.size() is the number of vertices, IDs from 0 to n-1
	INs: adj_lists: set of adjacency vertices that can reach v (with arc/edge weights)
	OUTs: // adj_lists: set of adjacency vertices that can be reached by v (with arc/edge weights)
	if v1 is in INs of v2, then v2 must be in OUTs of v1, and both records should be associated with the same weight.
	*/
	std::vector<std::vector<std::pair<int, weight_type>>> INs, OUTs;

	/*constructors*/
	dgraph_v_of_v() {}
	dgraph_v_of_v(int n) {
		INs.resize(n); // initialize n vertices
		OUTs.resize(n);
	}

    /**
     * @brief Get number of edges.
     */
    size_t edge_number();

	/*class member functions*/
	inline void add_edge(int, int, weight_type); // this function can change edge weights
	inline void remove_edge(int, int);
	inline weight_type edge_weight(int, int);
	inline bool contain_edge(int, int); // whether there is an edge
	inline long long int edge_number(); // the total number of edges
};

/* class member functions */




/*class member functions*/

template <typename weight_type>
void dgraph_v_of_v<weight_type>::add_edge(int v1, int v2, weight_type weight) { 

	/*
	edge direction: v1 to v2;
	this function adds v1 into INs of v2, and also adds v2 into OUTs of v1;
	this function can change edge weights;
	*/

	graph_hash_of_mixed_weighted_binary_operations_erase(OUTs[v1], v2);
	graph_hash_of_mixed_weighted_binary_operations_erase(INs[v2], v1);
}

template <typename weight_t>
weight_t dgraph_v_of_v<weight_t>::edge_weight(int v1, int v2) {

    /*edge direction: v1 to v2*/

	return graph_hash_of_mixed_weighted_binary_operations_search_weight(OUTs[v1], v2);
}


template <typename weight_type>
bool dgraph_v_of_v<weight_type>::contain_edge(int v1, int v2) {

	/*this function checks two conditions: v1 is in INs of v2; v2 is in OUTs of v1*/

	int x = graph_hash_of_mixed_weighted_binary_operations_search(OUTs[v1], v2) +
		graph_hash_of_mixed_weighted_binary_operations_search(INs[v2], v1);

	if (x == 0) {
		return false;
	}
	else if (x == 2) {
		return true;
	}
	else {
		std::cout << "dgraph_v_of_v check_edge error!" << std::endl;
		if (graph_hash_of_mixed_weighted_binary_operations_search(OUTs[v1], v2)) {
			std::cout << v2 << " is in OUTs of" << v1 << std::endl;
		}
		else {
			std::cout << v2 << " is NOT in OUTs of" << v1 << std::endl;
		}
		if (graph_hash_of_mixed_weighted_binary_operations_search(INs[v2], v1)) {
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
	for (int i = INs.size() - 1; i >= 0; i--) {
		num += INs[i].size();
	}

	return num;
}

/*
examples
----------------------
#include <dgraph_v_of_v/dgraph_v_of_v.h>

int main()
{
    example_dgraph_v_of_v();
}
-----------------------
*/

void example_dgraph_v_of_v() {

    dgraph_v_of_v<double> g(10);  // initialize a graph with 10 vertices

    g.add_edge(1, 5, 1.02);
    g.add_edge(5, 1, 1.42);
    g.add_edge(5, 2, 122);
    g.remove_edge(5, 1);

	std::cout << g.edge_weight(1, 5) << std::endl;
	std::cout << g.contain_edge(1, 5) << std::endl;
	std::cout << g.contain_edge(5, 1) << std::endl;
	std::cout << g.edge_number() << std::endl;
	std::cout << g.INs.size() << std::endl;
	

}
