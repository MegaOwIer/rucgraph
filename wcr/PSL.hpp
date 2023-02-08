/**
 * @file PSL.hpp
 * @brief Implementation of PSL algorithm.
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
#include <numeric>
#include <queue>
#include <shared_mutex>
#include <unordered_map>
#include <utility>
#include <vector>

template <class weight_t>
class PSL {
    using label_set = std::vector<std::vector<label<weight_t>>>;

    // node ranks
    std::vector<size_t> sorted_vertex;
    std::vector<size_t> rank;

    // results. 0 -> out, 1 -> in
    label_set L[2], L_temp[2];
    std::vector<size_t> endpos1[2];

    // for multi-thread
    std::queue<int> thread_id;
    std::shared_mutex mtx;
    std::vector<std::vector<size_t>> dirt;
    std::vector<std::vector<weight_t>> dmin;

    void propagate(const std::array<dgraph<weight_t>, 2> &G, int k, size_t u);

    void append(int k, size_t u);

    void shrink(int k, size_t u);

    void build_PSL_labels(const std::array<dgraph<weight_t>, 2> &G, int num_of_threads);

    // for debug
    void print_label_set(const std::vector<std::vector<label<weight_t>>> &_L);

public:
    /**
     * @brief Building process of PSL algorithm. Implemented according to Algorithm 1 in the paper.
     *
     * @tparam weight_t
     * @param g
     * @param num_of_threads
     * @param case_info
     */
    PSL(const dgraph<weight_t> &g, int num_of_threads, PSL_runtime_info &case_info);

    /**
     * @brief Query length of shortest path from u to v.
     */
    weight_t query_dist(size_t u, size_t v);
};

template <class weight_t>
PSL<weight_t>::PSL(const dgraph<weight_t> &g, int num_of_threads, PSL_runtime_info &case_info) {

    decltype(std::chrono::high_resolution_clock::now()) begin_time, end_time;

    /* ------- Initialization ------- */

    begin_time = std::chrono::high_resolution_clock::now();

    // std::vector<std::vector<pair<int, double>>>().swap(adjs);
    // adjs.resize(max_N_ID);
    // std::vector<pair<int, double>>().swap(min_adjs);
    // min_adjs.resize(max_N_ID);
    // std::vector<std::pair<int, int>> sorted_vertices;
    // for (auto it = g.hash_of_std::vectors.begin(); it != g.hash_of_std::vectors.end(); it++) {
    //     sorted_vertices.push_back({it->first, g.degree(it->first)});
    //     adjs[it->first] = g.adj_v_and_ec(it->first);
    //     min_adjs[it->first] = g.min_adj(it->first);
    // }
    // sort(sorted_vertices.begin(), sorted_vertices.end(), compare_pair_second_large_to_small);

    /*graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID*/
    // unordered_map<int, int> vertexID_old_to_new;
    // vertexID_new_to_old_595.resize(N);
    // for (int i = 0; i < N; i++) {
    //     vertexID_old_to_new[sorted_vertices[i].first] = i;
    //     vertexID_new_to_old_595[i] = sorted_vertices[i].first;
    // }
    // std::vector<pair<int, int>>().swap(sorted_vertices);
    // ideal_graph_595 = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(g,
    // vertexID_old_to_new);

    // /*redcution: add and remove certain edges*/
    // case_info.reduction_measures_2019R2.clear();  // for using this function multiple times
    // case_info.reduction_measures_2019R2.resize(max_N_ID);
    // /*clear case_info.f_2019R1*/
    // case_info.reduction_measures_2019R1.clear();  // for using this function multiple times
    // case_info.reduction_measures_2019R1.resize(max_N_ID);
    // case_info.f_2019R1.resize(max_N_ID);
    // std::iota(std::begin(case_info.f_2019R1), std::end(case_info.f_2019R1),
    //           0);  // Fill with 0, 1, ...

    end_time = std::chrono::high_resolution_clock::now();
    case_info.time_initialization =
        std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count();

    /* ------- Reduction ------- */

    /* redcution1: equivalence between vertices */
    if (case_info.use_2019R1) {
        // begin_time = std::chrono::high_resolution_clock::now();

        // ThreadPool pool(num_of_threads);
        // std::vector<std::future<int>> results;  // return typename: xxx
        // for (int i = 0; i < g.getV(); i++) {
        //     std::vector<int> *xx = &(case_info.reduction_measures_2019R1);
        //     std::vector<int> *yy = &(case_info.f_2019R1);
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
        // case_info.reduce_V_num_2019R1 = 0;
        // for (int i = 0; i < max_N_ID; i++) {
        //     if (case_info.f_2019R1[i] != i) {
        //         if (case_info.f_2019R1[case_info.f_2019R1[i]] != case_info.f_2019R1[i]) {
        //             cout << "f error due to the above parallelly updating f" << endl;
        //             // getchar();
        //         }
        //         case_info.reduce_V_num_2019R1++;
        //         graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_595,
        //                                                        vertexID_old_to_new[i]);
        //     }
        // }

        // end_time = std::chrono::high_resolution_clock::now();
        // case_info.time_2019R1 =
        //     std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;
    }

    // begin = std::chrono::high_resolution_clock::now();
    // if (case_info.use_2019R2) {  // no edge is removed here
    //     case_info.MG_num = 0;
    //     for (int x = 0; x < N; x++) {
    //         if (ideal_graph_595[x].size() > 0) {  // Prevent memory overflow
    //             if (x > ideal_graph_595[x][ideal_graph_595[x].size() - 1]
    //                         .first) {  // Here, only one comparison is needed. A little trick.
    //                 case_info.MG_num++;
    //                 case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] = 2;
    //                 // cout << "reduce " << vertexID_new_to_old_595[x] << endl;
    //             }
    //         }
    //     }
    // } else if (case_info.use_enhanced2019R2) {  // e.g., 2019R2enhance_10
    //     case_info.MG_num = 0;
    //     for (int x = 0; x < N; x++) {
    //         if (ideal_graph_595[x].size() > 0) {  // Prevent memory overflow
    //             if (x > ideal_graph_595[x][ideal_graph_595[x].size() - 1]
    //                         .first) {  // Here, only one comparison is needed. A little trick.
    //                 case_info.MG_num++;
    //                 case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] = 2;
    //                 // cout << "reduce " << vertexID_new_to_old_595[x] << endl;
    //             }
    //         }
    //     }
    //     int bound = case_info.max_degree_MG_enhanced2019R2;
    //     for (int x = N - 1; x >= 0; x--) {  // from low ranking to high ranking
    //         if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] == 0 &&
    //             ideal_graph_595[x].size() <= bound) {  // bound is the max degree for reduction
    //             bool no_adj_MG_vertices = true;
    //             for (auto it = ideal_graph_595[x].begin(); it != ideal_graph_595[x].end(); it++)
    //             {
    //                 if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[it->first]]
    //                 ==
    //                     2) {
    //                     no_adj_MG_vertices = false;
    //                     break;
    //                 }
    //             }
    //             if (no_adj_MG_vertices) {
    //                 case_info.MG_num++;
    //                 case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] =
    //                     2;  // new reduction
    //                         // cout << "new reduce " << vertexID_new_to_old_595[x] << endl;
    //             }
    //         }
    //     }
    // } else if (case_info.use_non_adj_reduc_degree) {  // e.g., 2019R2enhance_10
    //     case_info.MG_num = 0;
    //     int bound = case_info.max_degree_MG_enhanced2019R2;
    //     for (int x = N - 1; x >= 0; x--) {  // from low ranking to high ranking
    //         if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] == 0 &&
    //             ideal_graph_595[x].size() <= bound) {  // bound is the max degree for reduction
    //             bool no_adj_MG_vertices = true;
    //             for (auto it = ideal_graph_595[x].begin(); it != ideal_graph_595[x].end(); it++)
    //             {
    //                 if (case_info.reduction_measures_2019R2[vertexID_new_to_old_595[it->first]]
    //                 ==
    //                     2) {
    //                     no_adj_MG_vertices = false;
    //                     break;
    //                 }
    //             }
    //             if (no_adj_MG_vertices) {
    //                 case_info.MG_num++;
    //                 case_info.reduction_measures_2019R2[vertexID_new_to_old_595[x]] =
    //                     2;  // new reduction
    //                         // cout << "new reduce " << vertexID_new_to_old_595[x] << endl;
    //             }
    //         }
    //     }
    // }
    // end = std::chrono::high_resolution_clock::now();
    // case_info.time_2019R2_or_enhanced_pre =
    //     std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count() /
    //     1e9;

    /* ------- Generate Labels ------- */
    begin_time = std::chrono::high_resolution_clock::now();

    // prepare for thread pool
    dmin.resize(num_of_threads);
    dirt.resize(num_of_threads);
    for (int i = 0; i < num_of_threads; i++) {
        dmin[i].resize(g.getV());
        dirt[i].assign(g.getV(), -1U);
        thread_id.push(i);
    }

    // build labels
    std::array<dgraph<weight_t>, 2> G {dgraph<weight_t>(g), dgraph<weight_t>(g, true)};
    G[0].sort_nodes();
    G[1].sort_nodes();
    build_PSL_labels(G, num_of_threads);

    // free heaps
    L_temp[0].clear();
    L_temp[1].clear();
    endpos1[0].clear();
    endpos1[1].clear();
    dirt.clear();
    dmin.clear();

    end_time = std::chrono::high_resolution_clock::now();
    case_info.time_generate_labels =
        std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count();
    //---------------------------------------------------------------------------------------------------------------------------------------

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
    // if (case_info.use_canonical_repair) {
    //     begin = std::chrono::high_resolution_clock::now();
    //     reduction_measures_2019R1_new_ID.resize(max_N_ID, 0);
    //     reduction_measures_2019R2_new_ID.resize(max_N_ID, 0);
    //     f_2019R1_new_ID.resize(max_N_ID, 0);
    //     for (int i = 0; i < N; i++) {
    //         reduction_measures_2019R1_new_ID[vertexID_old_to_new[i]] =
    //             case_info.reduction_measures_2019R1[i];
    //         reduction_measures_2019R2_new_ID[vertexID_old_to_new[i]] =
    //             case_info.reduction_measures_2019R2[i];
    //         f_2019R1_new_ID[vertexID_old_to_new[i]] = vertexID_old_to_new[case_info.f_2019R1[i]];
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
    //     case_info.time_canonical_repair1 =
    //         std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;  //
    //         s

    //     begin = std::chrono::high_resolution_clock::now();
    //     canonical_repair_multi_threads(new_ID_g, case_info.label_size_before_canonical_repair,
    //                                    case_info.label_size_after_canonical_repair,
    //                                    case_info.canonical_repair_remove_label_ratio,
    //                                    num_of_threads);
    //     end = std::chrono::high_resolution_clock::now();
    //     case_info.time_canonical_repair2 =
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
    // case_info.L = graph_hash_of_mixed_weighted_PSL_v1_transform_labels_to_old_vertex_IDs(
    //     N, max_N_ID, num_of_threads);

    // end = std::chrono::high_resolution_clock::now();
    // case_info.time_update_old_IDs_in_labels =
    //     std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count() /
    //     1e9;
    //---------------------------------------------------------------------------------------------------------------------------------------

    // graph_hash_of_mixed_weighted_two_hop_clear_global_values();
}

template <class weight_t>
weight_t PSL<weight_t>::query_dist(size_t u, size_t v) {
    auto i = L[0][u].begin(), j = L[1][v].begin();
    weight_t ans = std::numeric_limits<weight_t>::max();
    while (i < L[0][u].end() && j < L[1][v].end()) {
        if (i->vertex < j->vertex) {
            i++;
        } else if (i->vertex > j->vertex) {
            j++;
        } else {
            ans = std::min(ans, i->dist + j->dist);
            i++, j++;
        }
    }
    return ans;
}

template <class weight_t>
void PSL<weight_t>::propagate(const std::array<dgraph<weight_t>, 2> &G, int k, size_t u) {
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

    for (auto [v, dd] : G[k][u]) {
        for (size_t i = endpos1[k][v]; i < L[k][v].size(); i++) {
            size_t x = L[k][v][i].vertex;
            if (G[k].rnk(x) < G[k].rnk(u)) {  // higher rank means lower value of `rank`
                continue;
            }
            weight_t d_temp = dd + L[k][v][i].dist;

            // Line 12-16
            bool flag = false;
            for (size_t j = endpos1[k ^ 1][x]; j < L[k ^ 1][x].size(); j++) {
                size_t y = L[k ^ 1][x][j].vertex;
                if (!dirt[current_tid][y] && dmin[current_tid][y] + L[k ^ 1][x][j].dist <= d_temp) {
                    flag = true;
                }
            }
            if (flag) {
                continue;
            }

            // Line 17
            L_temp[k][u].emplace_back(x, d_temp, v);
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
void PSL<weight_t>::build_PSL_labels(const std::array<dgraph<weight_t>, 2> &G,
                                     int num_of_threads) {
    // Line 1-2
    for (int k = 0; k < 2; k++) {
        L[k].resize(G[k].getV());
        L_temp[k].resize(G[k].getV());
        endpos1[k].resize(G[k].getV());
        for (size_t i = 0; i < G[k].getV(); i++) {
            L[k][i].emplace_back(i, 0, i);
        }
    }

    ThreadPool pool(num_of_threads);
    std::vector<std::future<int>> results;
    // int num_of_threads_per_push = num_of_threads * 1000;
    // ÿ��push��ȥ num_of_threads_per_push
    // �̣߳����û���쳣������push��ȥnum_of_threads_per_push�̣߳����ȫ��һ��push��ȥ����ȫ���̶߳���������catch�쳣

    // since R2 skip some vertices, some new labels can only be
    // generated when d increases 2, not 1, thus terminate the
    // loop only when if_continue_595==false twice
    while (true) {
        std::cout << "[Iter - Out]" << std::endl;
        print_label_set(L[0]);
        std::cout << "[Iter - In]" << std::endl;
        print_label_set(L[1]);
        std::cout << std::endl;

        // Line 4-20
        bool is_empty[2] = {true, true};

        for (int k = 0; k < 2; k++) {
            for (size_t u = 0; u < G[k].getV(); u++) {
                // if (ideal_graph_595[u].size() == 0 ||
                //     case_info.reduction_measures_2019R2[vertexID_new_to_old_595[u]] == 2)
                //     continue;  // do not search isolated vertices
                // auto *xx = &case_info;

                // new labels add into L_temp[u], but also read L in the process
#ifndef DISABLE_MULTITHREAD
                results.emplace_back(pool.enqueue([k, u, this, &G] {
#endif
                    this->propagate(G, k, u);
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

            for (size_t u = 0; u < G[k].getV(); u++) {
                // if (ideal_graph_595[u].size() == 0 ||
                //     case_info.reduction_measures_2019R2[vertexID_new_to_old_595[u]] == 2)
                //     continue;                            // do not search isolated vertices

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

            for (size_t u = 0; u < G[k].getV(); u++) {
                if (!L_temp[k][u].empty()) {
                    is_empty[k] = false;
                }
                endpos1[k][u] = L[k][u].size() - L_temp[k][u].size();
                L_temp[k][u].clear();
            }
        }  // for loop to choose graph

        if (is_empty[0] && is_empty[1]) {
            break;
        }
    }  // while loop

    // remove redundant labels in L
    for (int k = 0; k < 2; k++) {
        for (size_t u = 0; u < G[k].getV(); u++) {
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
void PSL<weight_t>::print_label_set(const std::vector<std::vector<label<weight_t>>> &_L) {
    for (size_t i = 0; i < _L.size(); i++) {
        std::cout << "[" << i << "] ";
        for (const label<weight_t> &cur : _L[i]) {
            std::cout << cur << " ";
        }
        std::cout << "\n";
    }
}
