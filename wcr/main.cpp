/**
 * @brief Sample testing codes for PSL for Directed graphs.
 */

#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "PSL.hpp"
#include "test/dijkstra.hpp"
#include "utility/config.h"
#include "utility/dgraph.hpp"

const int MXVAL = 1000000;

void run_small_sample() {
    PSL::dgraph<int> g(5);

    // Sample graph #1
    // g.add_edge(0, 1, 1);
    // g.add_edge(0, 4, 2);
    // g.add_edge(0, 5, 4);
    // g.add_edge(1, 2, 4);
    // g.add_edge(1, 4, 2);
    // g.add_edge(2, 3, 3);
    // g.add_edge(2, 5, 1);
    // g.add_edge(4, 3, 6);
    // g.add_edge(4, 5, 1);
    // g.add_edge(5, 3, 1);

    // Sample graph #2
    // g.add_edge(0, 2, 1);
    // g.add_edge(1, 2, 1);
    // g.add_edge(2, 3, 1);
    // g.add_edge(3, 4, 1);
    // g.add_edge(3, 5, 1);

    // Sample graph #3
    // g.add_edge(0, 1, 1);
    // g.add_edge(1, 2, 2);
    // g.add_edge(2, 3, 4);
    // g.add_edge(3, 0, 8);
    // g.add_edge(0, 0, 16);

    // Sample graph #4
    g.add_edge(4, 3, 219408);
    g.add_edge(0, 1, 948346);
    g.add_edge(2, 0, 672150);
    g.add_edge(2, 2, 128888);
    g.add_edge(1, 2, 49863);
    g.add_edge(3, 2, 996739);
    g.add_edge(3, 3, 223673);
    g.add_edge(2, 4, 504573);

    /* reduction method selection*/
    PSL::runtime_info cfg;
    cfg.use_lms = true;

    PSL::PSL<int> solve(g, 6, cfg);

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            std::cout << solve.query_dist(i, j) << "\t";
        }
        std::cout << "\n";
    }
}

void display(int cur, int tot) {
    std::cout << "\r(" << cur << "/" << tot << ") [";
    int cnt = 50. * cur / tot;
    for (int i = 1; i <= 50; i++) {
        std::cout << (i <= cnt ? ">" : "-");
    }
    std::cout << "]" << std::flush;
}

void run_test(std::string file_name) {
    std::cout << file_name << std::endl;
    int q = 10000;
    std::ifstream fin(file_name);

    // generating graph
    std::mt19937 gen(19260817);
    std::uniform_int_distribution<> distrib(1, MXVAL);

    int cntV = 0;
    std::unordered_map<std::string, int> NodeID;
    std::vector<std::tuple<int, int, int>> Edges;
    std::vector<std::pair<int, int>> Queries;

    {
        std::string u, v;
        while (fin >> u >> v) {
            int uid = NodeID.count(u) ? NodeID[u] : NodeID[u] = cntV++;
            int vid = NodeID.count(v) ? NodeID[v] : NodeID[v] = cntV++;
            if (uid != vid) {
                Edges.emplace_back(uid, vid, distrib(gen));
            }
        }
    }
    for (int i = 1; i <= q; i++) {
        Queries.emplace_back(distrib(gen) % cntV, distrib(gen) % cntV);
    }
    std::cerr << "Load graph finish." << std::endl;

    decltype(std::chrono::high_resolution_clock::now()) begin_time, end_time;

    if(0) {  // naive PSL
        std::cerr << "trivial PSL" << std::flush;

        PSL::dgraph<int64_t> g(cntV);
        for (auto [u, v, w] : Edges) {
            g.add_edge(u, v, w);
        }

        PSL::runtime_info cfg;
        cfg.use_lms = false;

        PSL::PSL<int64_t> solve(g, 48, cfg);

        begin_time = std::chrono::high_resolution_clock::now();
        int64_t ans = 0;
        for (auto [s, t] : Queries) {
            ans ^= solve.query_dist(s, t);
        }
        end_time = std::chrono::high_resolution_clock::now();

        cfg.display();
        std::cout
            << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count()
            << ' ' << ans << std::endl;
        std::cerr << " finish." << std::endl;
    }

    if(0) {  // optimized PSL
        std::cerr << "optimized PSL" << std::flush;

        PSL::dgraph<int64_t> g(cntV);
        for (auto [u, v, w] : Edges) {
            g.add_edge(u, v, w);
        }

        PSL::runtime_info cfg;
        cfg.use_lms = true;

        PSL::PSL<int64_t> solve(g, 48, cfg);

        begin_time = std::chrono::high_resolution_clock::now();
        int64_t ans = 0;
        for (auto [s, t] : Queries) {
            ans ^= solve.query_dist(s, t);
        }
        end_time = std::chrono::high_resolution_clock::now();

        cfg.display();
        std::cout
            << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count()
            << ' ' << ans << std::endl;
        std::cerr << " finish." << std::endl;
    }

    if(1) {  // Dijkstra
        std::cerr << "Dijkstra" << std::flush;
        // Initialize Dijkstra
        PSL::Dijkstra<int64_t> dij_sol(cntV);
        for (auto [u, v, w] : Edges) {
            dij_sol.add_edge(u, v, w);
        }

        begin_time = std::chrono::high_resolution_clock::now();
        int64_t ans = 0;
        for (auto [s, t] : Queries) {
            ans ^= dij_sol(s, t);
        }
        end_time = std::chrono::high_resolution_clock::now();

        std::cout
            << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count()
            << ' ' << ans << std::endl;
        std::cerr << " finish." << std::endl;
    }

    std::cout << std::endl;
}

void run_small_random_sample_correctness_test(int n, int m, int seed = 19260817) {
    // generating graph
    std::mt19937 gen(seed);
    std::uniform_int_distribution<> distrib(1, std::numeric_limits<int>::max());

    std::vector<std::tuple<int, int, int>> Edges;
    for (int i = 1; i <= m; i++) {
        int u = distrib(gen) % n, v = distrib(gen) % n;
        Edges.emplace_back(u, v, distrib(gen) % MXVAL + 1);
    }
    std::cout << "Load graph finish." << std::endl;

    // for (auto [u, v, w] : Edges) {
    //     std::cout << u << " " << v << " " << w << std::endl;
    // }

    // Initialize Dijkstra
    PSL::Dijkstra<int64_t> dij_sol(n);
    for (auto [u, v, w] : Edges) {
        dij_sol.add_edge(u, v, w);
    }

    std::cout << "Initialize Dijkstra finish." << std::endl;

    // PSL build label
    PSL::dgraph<int64_t> g(n);
    for (auto [u, v, w] : Edges) {
        g.add_edge(u, v, w);
    }

    PSL::runtime_info cfg;
    cfg.use_lms = true;

    PSL::PSL<int64_t> solve(g, 48, cfg);
    std::cout << "PSL build label finish." << std::endl;

    // test
    std::cout << "Testing..." << std::endl;
    int TOTAL_ROUND = 1000;
    std::uniform_int_distribution<> rnd(0, n - 1);
    for (int i = 1; i <= TOTAL_ROUND; i++) {
        int u = rnd(gen), v = rnd(gen);
        if (dij_sol(u, v) != solve.query_dist(u, v)) {
            printf("[Error] Query dist(%d, %d) : expect %ld, got %ld.", u, v, dij_sol(u, v),
                   solve.query_dist(u, v));
            exit(0);
        }
        display(i, TOTAL_ROUND);
    }
    std::cout << std::endl;
}

void run_small_random_sample_test(int n, int m, int q) {
    std::cout << n << ' ' << m << ' ' << q << std::endl;
    std::cerr << n << ' ' << m << ' ' << q << std::endl;

    // generating graph
    std::mt19937 gen(19260817);
    std::uniform_int_distribution<> distrib(1, std::numeric_limits<int>::max());

    std::map<std::pair<int, int>, int> Edges;
    std::vector<std::pair<int, int>> Queries;
    for (int i = 1; i <= m; i++) {
        int u = distrib(gen) % n, v = distrib(gen) % n;
        int w = distrib(gen) % MXVAL + 1;
        if (u == v) {
            continue;
        }
        if (Edges.count(std::make_pair(u, v))) {
            Edges[{u, v}] = std::min(Edges[{u, v}], w);
        } else {
            Edges[{u, v}] = w;
        }
    }
    for (int i = 1; i <= q; i++) {
        Queries.emplace_back(distrib(gen) % n, distrib(gen) % n);
    }
    std::cerr << "Load graph finish." << std::endl;

    decltype(std::chrono::high_resolution_clock::now()) begin_time, end_time;

    {  // naive PSL
        std::cerr << "trivial PSL" << std::flush;

        PSL::dgraph<int64_t> g(n);
        for (auto [cur, w] : Edges) {
            g.add_edge(cur.first, cur.second, w);
        }

        PSL::runtime_info cfg;
        cfg.use_lms = false;

        PSL::PSL<int64_t> solve(g, 48, cfg);

        begin_time = std::chrono::high_resolution_clock::now();
        int64_t ans = 0;
        for (auto [s, t] : Queries) {
            ans ^= solve.query_dist(s, t);
        }
        end_time = std::chrono::high_resolution_clock::now();

        cfg.display();
        std::cout
            << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count()
            << ' ' << ans << std::endl;
        std::cerr << " finish." << std::endl;
    }

    {  // optimized PSL
        std::cerr << "optimized PSL" << std::flush;

        PSL::dgraph<int64_t> g(n);
        for (auto [cur, w] : Edges) {
            g.add_edge(cur.first, cur.second, w);
        }

        PSL::runtime_info cfg;
        cfg.use_lms = true;

        PSL::PSL<int64_t> solve(g, 48, cfg);

        begin_time = std::chrono::high_resolution_clock::now();
        int64_t ans = 0;
        for (auto [s, t] : Queries) {
            ans ^= solve.query_dist(s, t);
        }
        end_time = std::chrono::high_resolution_clock::now();

        cfg.display();
        std::cout
            << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count()
            << ' ' << ans << std::endl;
        std::cerr << " finish." << std::endl;
    }

    {  // Dijkstra
        std::cerr << "Dijkstra" << std::flush;
        // Initialize Dijkstra
        PSL::Dijkstra<int64_t> dij_sol(n);
        for (auto [cur, w] : Edges) {
            dij_sol.add_edge(cur.first, cur.second, w);
        }

        begin_time = std::chrono::high_resolution_clock::now();
        int64_t ans = 0;
        for (auto [s, t] : Queries) {
            ans ^= dij_sol(s, t);
        }
        end_time = std::chrono::high_resolution_clock::now();

        std::cout
            << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count()
            << ' ' << ans << std::endl;
        std::cerr << " finish." << std::endl;
    }

    std::cout << std::endl;
}

int main() {
    // run_small_sample();

    // std::mt19937 gen(19260817);
    // for(int i = 1; i <= 100; i++) {
    //     run_small_random_sample_correctness_test(100, 5000, gen());
    // }

    run_small_random_sample_test(1000, 20000, 10000);
    run_small_random_sample_test(5000, 200000, 10000);
    run_small_random_sample_test(10000, 1000000, 10000);
    run_small_random_sample_test(30000, 5000000, 10000);

    // run_test("/home/wangcr/test_data/wiki-Vote/wiki-Vote.txt");
    // run_test("/home/wangcr/test_data/soc-Slashdot/soc-Slashdot0902.txt");
    // run_test("/home/wangcr/test_data/wiki-Talk/wiki-Talk.txt");
    // run_test("/home/wangcr/test_data/cit-HepPh/out.cit-HepPh");
    // run_test("/home/wangcr/test_data/munmun_twitter_social/out.munmun_twitter_social");
    return 0;
}
