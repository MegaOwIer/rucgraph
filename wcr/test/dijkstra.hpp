#pragma once

#include <functional>
#include <queue>
#include <utility>
#include <vector>

namespace PSL {

template <class weight_t>
struct Dijkstra {
    std::vector<std::vector<std::pair<int, weight_t>>> G;

    Dijkstra(int V) : G(V) {}

    void add_edge(int u, int v, weight_t w) { G[u].emplace_back(v, w); }

    weight_t operator()(int s, int t) {
        std::vector<bool> vis(G.size(), false);
        std::vector<int64_t> dis(G.size(), std::numeric_limits<int64_t>::max());

        using info_t = std::pair<int64_t, int>;
        std::priority_queue<info_t, std::vector<info_t>, std::greater<info_t>> heap;

        heap.emplace(dis[s] = 0, s);
        while (!heap.empty()) {
            auto [d, u] = heap.top();
            heap.pop();
            if (u == t) {
                return d;
            }
            if (vis[u]) {
                continue;
            }

            vis[u] = true;
            for (auto [v, w] : G[u]) {
                if (dis[v] > d + w) {
                    heap.emplace(dis[v] = d + w, v);
                }
            }
        }

        // unreachable
        return std::numeric_limits<int64_t>::max();
    }
};

}  // namespace PSL
