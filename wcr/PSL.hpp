/**
 * @file PSL.hpp
 * @brief Implementation of PSL algorithm for directed graph.
 * @version 0.1
 * @ref <paper-bibinfo>
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include "tool_functions/ThreadPool.h"

#include "utility/config.h"
#include "utility/dgraph.hpp"
#include "utility/label.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <memory>
#include <numeric>
#include <queue>
#include <shared_mutex>
#include <unordered_map>
#include <utility>
#include <vector>

namespace PSL {

template <class weight_t>
class PSL {
    using label_set = std::vector<std::vector<label<weight_t>>>;

    // graph info & node ranks
    std::array<std::shared_ptr<dgraph<weight_t>>, 2> G;
    std::vector<size_t> sorted_vertex, rank;

    // results. 0 -> out, 1 -> in
    label_set L[2], L_temp[2];
    std::vector<size_t> endpos1[2];  // begin of L_{d-1}
    std::vector<size_t> endpos2[2];  // begin of L_{d-2}

    // for multi-thread
    std::queue<int> thread_id;
    std::shared_mutex mtx;
    std::vector<std::vector<size_t>> dirt;
    std::vector<std::vector<weight_t>> dmin;

    // Reduction 2 : Local Minimum Set
    bool use_lms = false;
    size_t MG_size;
    std::vector<size_t> is_lms;

    void propagate(int k, size_t u);

    void append(int k, size_t u);

    void shrink(int k, size_t u);

    void build_PSL_labels(int num_of_threads);

    // for debug
    void print_label_set(const label_set &_L);

public:
    /**
     * @brief Building process of PSL algorithm. Implemented according to Algorithm 1 in the paper.
     *
     * @tparam weight_t
     * @param g
     * @param num_of_threads
     * @param cfg
     */
    PSL(const dgraph<weight_t> &g, int num_of_threads, runtime_info &cfg);

    /**
     * @brief Query length of shortest path from u to v.
     */
    weight_t query_dist(size_t u, size_t v);
};

template <class weight_t>
PSL<weight_t>::PSL(const dgraph<weight_t> &g, int num_of_threads, runtime_info &cfg) {

    decltype(std::chrono::high_resolution_clock::now()) begin_time, end_time;

    /* ------- Initialization ------- */

    begin_time = std::chrono::high_resolution_clock::now();

    // Prepare for thread pool
    dmin.resize(num_of_threads);
    dirt.resize(num_of_threads);
    for (int i = 0; i < num_of_threads; i++) {
        dmin[i].resize(g.getV());
        dirt[i].assign(g.getV(), -1U);
        thread_id.push(i);
    }

    // Build graph and rev_graph
    G[0] = std::make_shared<dgraph<weight_t>>(g, false);
    G[1] = std::make_shared<dgraph<weight_t>>(g, true);

    // Sort nodes
    sorted_vertex.resize(G[0]->getV());
    rank.resize(G[0]->getV());

    std::iota(sorted_vertex.begin(), sorted_vertex.end(), 0);
    for (size_t i = 0; i < G[0]->getV(); i++) {
        rank[i] = (*G[0])[i].size() + (*G[1])[i].size();  // sum of in-deg and out-deg
    }
    std::stable_sort(sorted_vertex.begin(), sorted_vertex.end(),
                     [&](size_t u, size_t v) { return rank[u] > rank[v]; });

    for (size_t i = 0; i < G[0]->getV(); i++) {
        rank[sorted_vertex[i]] = i;
    }

    end_time = std::chrono::high_resolution_clock::now();
    cfg.time_initialization =
        std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count();

    /* ------- Reduction ------- */

    /* Reduction 1: Equivalence between Vertices (Not Implemented) */
    if (cfg.use_equiv) {
        // begin_time = std::chrono::high_resolution_clock::now();

        // ThreadPool pool(num_of_threads);
        // std::vector<std::future<int>> results;  // return typename: xxx
        // for (int i = 0; i < g.getV(); i++) {
        //     std::vector<int> *xx = &(cfg.reduction_measures_2019R1);
        //     std::vector<int> *yy = &(cfg.f_2019R1);
        //     results.emplace_back(
        //         pool.enqueue([i, ideal_graph_size, xx,
        //                       yy] {  // pass const type value j to thread; [] can be empty
        //             update_2019R1_condition_PSL_enhancedoriginalR2(i, ideal_graph_size, xx, yy);
        //             return 1;  // return to results; the return type must be the same with
        //             results
        //         }));
        // }
        // for (auto &&result : results)
        //     result.get();  // all threads finish here
        // /* remove edges */
        // cfg.reduce_V_num_2019R1 = 0;
        // for (int i = 0; i < max_N_ID; i++) {
        //     if (cfg.f_2019R1[i] != i) {
        //         if (cfg.f_2019R1[cfg.f_2019R1[i]] != cfg.f_2019R1[i]) {
        //             cout << "f error due to the above parallelly updating f" << endl;
        //             // getchar();
        //         }
        //         cfg.reduce_V_num_2019R1++;
        //         graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_595,
        //                                                        vertexID_old_to_new[i]);
        //     }
        // }

        // end_time = std::chrono::high_resolution_clock::now();
        // cfg.time_2019R1 =
        //     std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
    }

    /* Reduction 2 : Local Minimum Set */
    if (cfg.use_lms) {
        begin_time = std::chrono::high_resolution_clock::now();

        use_lms = true;
        is_lms.assign(G[0]->getV(), true);
        for (int i = 0; i < 2; i++) {
            for (size_t u = 0; u < G[i]->getV(); u++) {
                for (auto &[v, _] : (*G[i])[u]) {
                    if (rank[u] < rank[v]) {
                        is_lms[u] = false;
                    }
                }
            }
        }
        MG_size = std::count(is_lms.begin(), is_lms.end(), true);

        // std::cout << "[Local Minimum]" << "\n";
        // for(int u = 0; u < G[0]->getV(); u++) {
        //     if(is_lms[u]) {
        //         std::cout << u << ' ';
        //     }
        // }
        // std::cout << std::endl;

        end_time = std::chrono::high_resolution_clock::now();
        cfg.time_rdc_lms =
            std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count();
    }
    // if (cfg.use_enhanced2019R2) {  // e.g., 2019R2enhance_10
    //     cfg.MG_num = 0;
    //     for (int x = 0; x < N; x++) {
    //         if (ideal_graph_595[x].size() > 0) {  // Prevent memory overflow
    //             if (x > ideal_graph_595[x][ideal_graph_595[x].size() - 1]
    //                         .first) {  // Here, only one comparison is needed. A little trick.
    //                 cfg.MG_num++;
    //                 cfg.reduction_measures_2019R2[vertexID_new_to_old_595[x]] = 2;
    //                 // cout << "reduce " << vertexID_new_to_old_595[x] << endl;
    //             }
    //         }
    //     }
    //     int bound = cfg.max_degree_MG_enhanced2019R2;
    //     for (int x = N - 1; x >= 0; x--) {  // from low ranking to high ranking
    //         if (cfg.reduction_measures_2019R2[vertexID_new_to_old_595[x]] == 0 &&
    //             ideal_graph_595[x].size() <= bound) {  // bound is the max degree for reduction
    //             bool no_adj_MG_vertices = true;
    //             for (auto it = ideal_graph_595[x].begin(); it != ideal_graph_595[x].end(); it++)
    //             {
    //                 if (cfg.reduction_measures_2019R2[vertexID_new_to_old_595[it->first]]
    //                 ==
    //                     2) {
    //                     no_adj_MG_vertices = false;
    //                     break;
    //                 }
    //             }
    //             if (no_adj_MG_vertices) {
    //                 cfg.MG_num++;
    //                 cfg.reduction_measures_2019R2[vertexID_new_to_old_595[x]] =
    //                     2;  // new reduction
    //                         // cout << "new reduce " << vertexID_new_to_old_595[x] << endl;
    //             }
    //         }
    //     }
    // } else if (cfg.use_non_adj_reduc_degree) {  // e.g., 2019R2enhance_10
    //     cfg.MG_num = 0;
    //     int bound = cfg.max_degree_MG_enhanced2019R2;
    //     for (int x = N - 1; x >= 0; x--) {  // from low ranking to high ranking
    //         if (cfg.reduction_measures_2019R2[vertexID_new_to_old_595[x]] == 0 &&
    //             ideal_graph_595[x].size() <= bound) {  // bound is the max degree for reduction
    //             bool no_adj_MG_vertices = true;
    //             for (auto it = ideal_graph_595[x].begin(); it != ideal_graph_595[x].end(); it++)
    //             {
    //                 if (cfg.reduction_measures_2019R2[vertexID_new_to_old_595[it->first]]
    //                 ==
    //                     2) {
    //                     no_adj_MG_vertices = false;
    //                     break;
    //                 }
    //             }
    //             if (no_adj_MG_vertices) {
    //                 cfg.MG_num++;
    //                 cfg.reduction_measures_2019R2[vertexID_new_to_old_595[x]] =
    //                     2;  // new reduction
    //                         // cout << "new reduce " << vertexID_new_to_old_595[x] << endl;
    //             }
    //         }
    //     }
    // }

    /* ------- Generate Labels ------- */
    begin_time = std::chrono::high_resolution_clock::now();

    // build labels
    build_PSL_labels(num_of_threads);

    // free heaps
    L_temp[0].clear();
    L_temp[1].clear();
    endpos1[0].clear();
    endpos1[1].clear();
    endpos2[0].clear();
    endpos2[1].clear();
    dirt.clear();
    dmin.clear();

    end_time = std::chrono::high_resolution_clock::now();
    cfg.time_generate_labels =
        std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count();

    //----------------------------------------------- step 4: canonical_repair
    //---------------------------------------------------------------

    // cout << "print_L_temp_595:" << endl;
    // for (int i = 0; i < L_temp_595.size(); i++) {
    //	cout << "L[" << i << "]=";
    //	for (int j = 0; j < L_temp_595[i].size(); j++) {
    //		cout << "{" << L_temp_595[i][j].vertex << "," << L_temp_595[i][j].distance << "," <<
    // L_temp_595[i][j].parent_vertex << "}";
    //	}
    //	cout << endl;
    // }

    /*canonical_repair based on the sorted new ID order, not the original ID order!*/
    // if (cfg.use_canonical_repair) {
    //     begin = std::chrono::high_resolution_clock::now();
    //     reduction_measures_2019R1_new_ID.resize(max_N_ID, 0);
    //     reduction_measures_2019R2_new_ID.resize(max_N_ID, 0);
    //     f_2019R1_new_ID.resize(max_N_ID, 0);
    //     for (int i = 0; i < N; i++) {
    //         reduction_measures_2019R1_new_ID[vertexID_old_to_new[i]] =
    //             cfg.reduction_measures_2019R1[i];
    //         reduction_measures_2019R2_new_ID[vertexID_old_to_new[i]] =
    //             cfg.reduction_measures_2019R2[i];
    //         f_2019R1_new_ID[vertexID_old_to_new[i]] = vertexID_old_to_new[cfg.f_2019R1[i]];
    //         sort(L_temp_595[i].begin(), L_temp_595[i].end(),
    //              compare_two_hop_label_small_to_large);  // sort is necessary
    //     }
    //     graph_hash_of_mixed_weighted new_ID_g =
    //         graph_hash_of_mixed_weighted_update_vertexIDs(g, vertexID_old_to_new);
    //     std::vector<std::vector<pair<int, double>>>().swap(adjs_new_IDs);
    //     adjs_new_IDs.resize(max_N_ID);
    //     std::vector<pair<int, double>>().swap(min_adjs_new_IDs);
    //     min_adjs_new_IDs.resize(max_N_ID);
    //     for (auto it = new_ID_g.hash_of_std::vectors.begin();
    //          it != new_ID_g.hash_of_std::vectors.end(); it++) {
    //         adjs_new_IDs[it->first] = new_ID_g.adj_v_and_ec(it->first);
    //         min_adjs_new_IDs[it->first] = new_ID_g.min_adj(it->first);
    //     }
    //     end = std::chrono::high_resolution_clock::now();
    //     cfg.time_canonical_repair1 =
    //         std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;  //
    //         s

    //     begin = std::chrono::high_resolution_clock::now();
    //     canonical_repair_multi_threads(new_ID_g, cfg.label_size_before_canonical_repair,
    //                                    cfg.label_size_after_canonical_repair,
    //                                    cfg.canonical_repair_remove_label_ratio,
    //                                    num_of_threads);
    //     end = std::chrono::high_resolution_clock::now();
    //     cfg.time_canonical_repair2 =
    //         std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;  //
    //         s
    // }

    // cout << "print_L_temp_595:" << endl;
    // for (int i = 0; i < L_temp_595.size(); i++) {
    //	cout << "L[" << i << "]=";
    //	for (int j = 0; j < L_temp_595[i].size(); j++) {
    //		cout << "{" << L_temp_595[i][j].vertex << "," << L_temp_595[i][j].distance << "," <<
    // L_temp_595[i][j].parent_vertex << "}";
    //	}
    //	cout << endl;
    // }

    //---------------------------------------------------------------------------------------------------------------------------------------

    //----------------------------------------------- step 5: update_old_IDs_in_labels
    //---------------------------------------------------------------
    // begin = std::chrono::high_resolution_clock::now();

    // /*return unordered_map_L for old_IDs*/
    // cfg.L = graph_hash_of_mixed_weighted_PSL_v1_transform_labels_to_old_vertex_IDs(
    //     N, max_N_ID, num_of_threads);

    // end = std::chrono::high_resolution_clock::now();
    // cfg.time_update_old_IDs_in_labels =
    //     std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count() /
    //     1e9;
    //---------------------------------------------------------------------------------------------------------------------------------------

    // graph_hash_of_mixed_weighted_two_hop_clear_global_values();
}

template <class weight_t>
weight_t PSL<weight_t>::query_dist(size_t u, size_t v) {
    std::unordered_map<size_t, weight_t> H[2];
    for (int k = 0; k < 2; k++) {
        size_t x = k ? v : u;
        if (use_lms && is_lms[x]) {
            for (auto &[w, d] : (*G[k])[x]) {
                for (const label<weight_t> &cur : L[k][w]) {
                    if (H[k].count(cur.vertex)) {
                        H[k][cur.vertex] = std::min(H[k][cur.vertex], d + cur.dist);
                    } else {
                        H[k][cur.vertex] = d + cur.dist;
                    }
                }
            }
            H[k][x] = 0;
        } else {
            for (const label<weight_t> &cur : L[k][x]) {
                if (H[k].count(cur.vertex)) {
                    H[k][cur.vertex] = std::min(H[k][cur.vertex], cur.dist);
                } else {
                    H[k][cur.vertex] = cur.dist;
                }
            }
        }
    }

    weight_t ans = std::numeric_limits<weight_t>::max();
    int x = H[0].size() > H[1].size();
    for (const auto &[w, d] : H[x]) {
        if (H[x ^ 1].count(w)) {
            ans = std::min(ans, d + H[x ^ 1][w]);
        }
    }
    return ans;
}

template <class weight_t>
void PSL<weight_t>::propagate(int k, size_t u) {
    mtx.lock();
    int current_tid = thread_id.front();
    thread_id.pop();
    mtx.unlock();

    // Line 5
    for (const label<weight_t> &cur : L[k][u]) {
        size_t v = cur.vertex;
        weight_t d = cur.dist;
        if (dirt[current_tid][v] == -1U) {
            dirt[current_tid][v] = 0;
            dmin[current_tid][v] = d;
        } else {
            dmin[current_tid][v] = std::min(dmin[current_tid][v], d);
        }
    }

    std::function<bool(size_t, weight_t)> is_redundant = [&](size_t x, weight_t d_temp) {
        for (size_t j = endpos1[k ^ 1][x]; j < L[k ^ 1][x].size(); j++) {
            size_t y = L[k ^ 1][x][j].vertex;
            if (!dirt[current_tid][y] && dmin[current_tid][y] + L[k ^ 1][x][j].dist <= d_temp) {
                return true;
            }
        }
        return false;
    };

    for (auto [v, dd] : (*G[k])[u]) {
        if (use_lms && is_lms[v]) {
            // Redction 2: if v in M(G), enumerate neighbor of neighbor
            for (auto [z, dd2] : (*G[k])[v]) {
                for (size_t i = endpos2[k][z]; i < endpos1[k][z]; i++) {
                    size_t x = L[k][z][i].vertex;  // x is a hub in L_{d-2}(z)
                    if (rank[x] > rank[z]) {
                        continue;
                    }

                    weight_t d_temp = dd + dd2 + L[k][z][i].dist;  // u -> v -> z --> x
                    if (!is_redundant(x, d_temp)) {
                        L_temp[k][u].emplace_back(x, d_temp, v);
                    }
                }
            }
        } else {
            for (size_t i = endpos1[k][v]; i < L[k][v].size(); i++) {
                size_t x = L[k][v][i].vertex;  // x is a hub in L_{d-1}(v)
                if (rank[x] > rank[u]) {       // higher rank means lower value of `rank`
                    continue;
                }

                weight_t d_temp = dd + L[k][v][i].dist;  // u -> v --> x

                // Line 12-17
                if (!is_redundant(x, d_temp)) {
                    L_temp[k][u].emplace_back(x, d_temp, v);
                }
            }
        }
    }

    for (const label<weight_t> &cur : L[k][u]) {
        size_t v = cur.vertex;
        dirt[current_tid][v] = -1U;
    }

    mtx.lock();
    thread_id.push(current_tid);
    mtx.unlock();
}

template <class weight_t>
void PSL<weight_t>::append(int k, size_t u) {
    L[k][u].insert(L[k][u].end(), L_temp[k][u].begin(), L_temp[k][u].end());
    L[k][u].shrink_to_fit();
}

template <class weight_t>
void PSL<weight_t>::shrink(int k, size_t u) {
    mtx.lock();
    int current_tid = thread_id.front();
    thread_id.pop();
    mtx.unlock();

    for (const label<weight_t> &cur : L[k][u]) {
        size_t v = cur.vertex;
        weight_t d = cur.dist;
        if (dirt[current_tid][v] == -1U) {
            dmin[current_tid][v] = d;
            dirt[current_tid][v] = cur.prev;
        } else if (d < dmin[current_tid][v]) {
            dmin[current_tid][v] = d;
            dirt[current_tid][v] = cur.prev;
        }
    }

    L[k][u].clear();
    for (size_t i = 0; i < dirt[current_tid].size(); i++) {
        if (dirt[current_tid][i] == -1U) {
            continue;
        }
        L[k][u].emplace_back(i, dmin[current_tid][i], dirt[current_tid][i]);
        dirt[current_tid][i] = -1U;
    }

    mtx.lock();
    thread_id.push(current_tid);
    mtx.unlock();
}

template <class weight_t>
void PSL<weight_t>::build_PSL_labels(int num_of_threads) {
    // Line 1-2
    for (int k = 0; k < 2; k++) {
        L[k].resize(G[k]->getV());
        L_temp[k].resize(G[k]->getV());
        endpos1[k].resize(G[k]->getV());
        for (size_t i = 0; i < G[k]->getV(); i++) {
            L[k][i].emplace_back(i, 0, i);
        }
    }

    // when reduction 2 is enabled, generate labels for d = 1 manually.
    if (use_lms) {
        for (int k = 0; k < 2; k++) {
            endpos2[k].resize(G[k]->getV());

            for (size_t u = 0; u < G[k]->getV(); u++) {
                if (is_lms[u]) {
                    L[k][u].clear();
                    continue;
                }
                endpos1[k][u] = L[k][u].size();
                for (auto &[v, w] : (*G[k])[u]) {
                    if (rank[u] > rank[v]) {
                        L[k][u].emplace_back(v, w, v);
                    }
                }
            }
        }
    }

    ThreadPool pool(num_of_threads);
    std::vector<std::future<int>> results;
    // int num_of_threads_per_push = num_of_threads * 1000;
    // 每次push进去 num_of_threads_per_push
    // 线程，如果没有异常，继续push进去num_of_threads_per_push线程；如果全都一起push进去必须全部线程都结束才能catch异常

    // When reduction 2 is enabled, L_{d} is generated by both L_{d-1} and L_{d-2}, hence we can
    // stop only when both of them are empty
    bool empty_2 = false;

    while (true) {
        std::cout << "[Iter - Out]" << std::endl;
        print_label_set(L[0]);
        std::cout << "[Iter - In]" << std::endl;
        print_label_set(L[1]);
        std::cout << std::endl;

        // Line 4-20
        bool is_empty[2] = {true, true};

        for (int k = 0; k < 2; k++) {
            for (size_t u = 0; u < G[k]->getV(); u++) {
                if (use_lms && is_lms[u]) {
                    continue;
                }

                // new labels add into L_temp[u], but also read L in the process
#ifndef DISABLE_MULTITHREAD
                results.emplace_back(pool.enqueue([k, u, this] {
#endif
                    this->propagate(k, u);
#ifndef DISABLE_MULTITHREAD
                    return 1;
                }));
#endif
                // push_num++;
                // if (push_num % num_of_threads_per_push == 0) {
                //     for (auto &&result : results)
                //         result.get();  // all threads finish here
                //     results.clear();
                // }
            }

#ifndef DISABLE_MULTITHREAD
            for (auto &&result : results) {
                result.get();
            }
            results.clear();
#endif

            for (size_t u = 0; u < G[k]->getV(); u++) {
                if (use_lms && is_lms[u]) {
                    continue;
                }

                // new labels in L_temp[u] add into L[u], to avoid locking L in propagate process
#ifndef DISABLE_MULTITHREAD
                results.emplace_back(pool.enqueue([k, u, this] {
#endif
                    this->append(k, u);
#ifndef DISABLE_MULTITHREAD
                    return 1;
                }));
#endif
            }

#ifndef DISABLE_MULTITHREAD
            for (auto &&result : results) {
                result.get();
            }
            results.clear();
#endif

            if (use_lms) {
                std::copy(endpos1[k].begin(), endpos1[k].end(), endpos2[k].begin());
            }
            for (size_t u = 0; u < G[k]->getV(); u++) {
                if (!L_temp[k][u].empty()) {
                    is_empty[k] = false;
                }
                endpos1[k][u] = L[k][u].size() - L_temp[k][u].size();
                L_temp[k][u].clear();
            }
        }  // for loop to choose graph

        if (is_empty[0] && is_empty[1] && (!use_lms || empty_2)) {
            break;
        }
        empty_2 = is_empty[0] && is_empty[1];
    }  // while loop

    // remove redundant labels in L
    for (int k = 0; k < 2; k++) {
        for (size_t u = 0; u < G[k]->getV(); u++) {
#ifndef DISABLE_MULTITHREAD
            results.emplace_back(pool.enqueue([k, u, this] {
#endif
                this->shrink(k, u);
#ifndef DISABLE_MULTITHREAD
                return 1;
            }));
#endif
        }
    }

    std::cout << "[Shrinked - Out]" << std::endl;
    print_label_set(L[0]);
    std::cout << "[Shrinked - In]" << std::endl;
    print_label_set(L[1]);
    std::cout << std::endl;

#ifndef DISABLE_MULTITHREAD
    for (auto &&result : results) {
        result.get();
    }
    results.clear();
#endif
}

template <class weight_t>
void PSL<weight_t>::print_label_set(const label_set &_L) {
    for (size_t i = 0; i < _L.size(); i++) {
        std::cout << "[" << i << "] ";
        for (const label<weight_t> &cur : _L[i]) {
            std::cout << cur << " ";
        }
        std::cout << "\n";
    }
}

}  // namespace PSL
