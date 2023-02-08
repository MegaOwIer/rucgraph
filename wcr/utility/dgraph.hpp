#pragma once

#include <algorithm>
#include <numeric>
#include <utility>
#include <vector>

/**
 * @brief Directed graph with edge weights for PSL. Nodes are labeled from 0.
 *
 * @tparam weight_t
 */
template <class weight_t>
class dgraph {
    std::vector<std::vector<std::pair<size_t, weight_t>>> adj_lists;
    std::vector<size_t> sorted_vertex;
    std::vector<size_t> rank;

public:
    dgraph() = delete;
    dgraph(const dgraph &&g) = delete;

    /**
     * @brief Construct a new dgraph object
     *
     * @param n Number of nodes.
     */
    dgraph(size_t n) { adj_lists.resize(n); }

    /**
     * @brief Copy constructor
     *
     * @param g
     * @param reverse whether to reverse direction of edges
     */
    dgraph(const dgraph &g, bool reverse = false);

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
     * @return const std::vector<std::pair<size_t, weight_t>>&
     */
    const std::vector<std::pair<size_t, weight_t>> &operator[](size_t u) const {
        return adj_lists[u];
    }

    /**
     * @brief Add an edge from v1 to v2 with given weight. No check to prevent duplicate edges.
     */
    void add_edge(size_t v1, size_t v2, weight_t weight);

    /**
     * @brief Sort nodes in out degree DESC.
     */
    void sort_nodes();

    /**
     * @brief Return node rank of node u.
     */
    size_t rnk(size_t u) const;
};

template <class weight_t>
dgraph<weight_t>::dgraph(const dgraph &g, bool reverse)
    : dgraph(g.getV()) {

    for (size_t u = 0; u < g.getV(); u++) {
        for (auto [v, w] : g[u]) {
            if (reverse) {
                add_edge(v, u, w);
            } else {
                add_edge(u, v, w);
            }
        }
    }
}

/* class member functions */

template <class weight_t>
size_t dgraph<weight_t>::edge_number() {
    size_t num = 0;
    for (auto &cur : adj_lists) {
        num += cur.size();
    }
    return num;
}

template <class weight_t>
void dgraph<weight_t>::add_edge(size_t v1, size_t v2, weight_t weight) {
    adj_lists[v1].emplace_back(v2, weight);
}

template <class weight_t>
void dgraph<weight_t>::sort_nodes() {
    sorted_vertex.resize(getV());
    std::iota(sorted_vertex.begin(), sorted_vertex.end(), 0);
    std::stable_sort(sorted_vertex.begin(), sorted_vertex.end(),
                     [&](size_t u, size_t v) { return adj_lists[u].size() > adj_lists[v].size(); });

    rank.resize(getV());
    for (size_t i = 0; i < getV(); i++) {
        rank[sorted_vertex[i]] = i;
    }
}

template <class weight_t>
size_t dgraph<weight_t>::rnk(size_t u) const {
    return rank.at(u);
}
