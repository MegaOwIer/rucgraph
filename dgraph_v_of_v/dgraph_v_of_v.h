#pragma once

#include <iostream>
#include <vector>

#include "graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_binary_operations.h"

/**
 * @brief Directed graph with edge weights. Nodes are labeled from 0.
 *
 * @tparam weight_t
 */
template <typename weight_t>
class dgraph_v_of_v {
    std::vector<std::vector<std::pair<int, weight_t>>> adj_lists;

public:
    dgraph_v_of_v() = delete;

    /**
     * @brief Construct a new dgraph_v_of_v object
     *
     * @param n Number of nodes.
     */
    dgraph_v_of_v(int n) { adj_lists.resize(n); }

    // TODO: dir = 0 -> copy, dir = 1 -> reverse
    dgraph_v_of_v(const dgraph_v_of_v &g, int dir = 0) {}

    /**
     * @brief Get number of vertices.
     */
    size_t getV() const { return adj_lists.size(); }

    /**
     * @brief Get number of edges.
     */
    size_t edge_number();

    /**
     * @brief Get all the edges starting from node u.
     * 
     * @return const std::vector<std::pair<int, weight_t>>& 
     */
    const std::vector<std::pair<int, weight_t>> &operator[](size_t u) const {
        return adj_lists[u];
    }

    /**
     * @brief Add an edge from v1 to v2 with given weight. If such edge exisits, update its weight.
     */
    void add_edge(int v1, int v2, weight_t weight);

    /**
     * @brief Remove the edge from v1 to v2.
     */
    void remove_edge(int, int);

    weight_t edge_weight(int, int);
    bool contain_edge(int, int);  // whether there is an edge
};

/* class member functions */

template <typename weight_t>
void dgraph_v_of_v<weight_t>::add_edge(int v1, int v2, weight_t weight) {
    graph_hash_of_mixed_weighted_binary_operations_insert(adj_lists[v1], v2, weight);
}

template <typename weight_t>
void dgraph_v_of_v<weight_t>::remove_edge(int v1, int v2) {
    graph_hash_of_mixed_weighted_binary_operations_erase(adj_lists[v1], v2);
}

template <typename weight_t>
weight_t dgraph_v_of_v<weight_t>::edge_weight(int v1, int v2) {

    /*edge direction: v1 to v2*/

    return graph_hash_of_mixed_weighted_binary_operations_search_weight(adj_lists[v1], v2);
}

/*this function checks two conditions: v1 is in INs of v2; v2 is in OUTs of
    * v1*/
template <typename weight_t>
bool dgraph_v_of_v<weight_t>::contain_edge(int v1, int v2) {
    return graph_hash_of_mixed_weighted_binary_operations_search(adj_lists[v1], v2);
}

template <typename weight_t>
size_t dgraph_v_of_v<weight_t>::edge_number() {
    size_t num = 0;
    for (auto &cur : adj_lists) {
        num += cur.size();
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
    std::cout << g.getV() << std::endl;
}
