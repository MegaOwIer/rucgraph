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

#include "dgraph_v_of_v/dgraph_v_of_v.h"
#include "tool_functions/ThreadPool.h"

#include "config.h"
#include "utility/label.h"

#include <algorithm>
#include <chrono>
#include <numeric>
#include <queue>
#include <shared_mutex>
#include <unordered_map>
#include <utility>
#include <vector>

template <class weight_t>
class PSL {
    std::vector<int> sorted_vertex;
    std::vector<size_t> rank;

    std::vector<std::vector<label<weight_t>>> L;
    std::vector<std::vector<label<weight_t>>> L_temp;
    std::vector<size_t> endpos1;

    std::queue<int> thread_id;
    std::shared_mutex mtx;
    std::vector<std::vector<int>> dirt;
    std::vector<std::vector<weight_t>> dmin;

    void propagate(const dgraph_v_of_v<weight_t> &g, int u);

    void append(int u);

    void shrink(int u);

public:
    PSL(const dgraph_v_of_v<weight_t> &g, int num_of_threads, PSL_runtime_info &case_info);

    weight_t query_dist(int u, int v) {}
};

/**
 * @brief Building process of PSL algorithm. Implemented according to Algorithm 1 in the paper.
 *
 * @tparam weight_t
 * @param g
 * @param num_of_threads
 * @param case_info
 */
template <class weight_t>
PSL<weight_t>::PSL(const dgraph_v_of_v<weight_t> &g, int num_of_threads,
                   PSL_runtime_info &case_info) {

    size_t Vnum = g.getV();
    std::chrono::_V2::system_clock::time_point begin_time, end_time;

    /* ------- Initialization ------- */

    begin_time = std::chrono::high_resolution_clock::now();

    // sort vertices by degrees
    sorted_vertex.resize(Vnum);
    std::iota(sorted_vertex.begin(), sorted_vertex.end(), 0);
    std::stable_sort(sorted_vertex.begin(), sorted_vertex.end(),
                     [&](int u, int v) { return g[u].size() > g[v].size(); });

    // node rank
    rank.resize(Vnum);
    for (size_t i = 0; i < Vnum; i++) {
        rank[sorted_vertex[i]] = i;
    }

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
        (double) std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count() / 1e9;

    /* ------- Reduction ------- */

    /* redcution1: equivalence between vertices */
    if (case_info.use_2019R1) {
        // begin_time = std::chrono::high_resolution_clock::now();

        // ThreadPool pool(num_of_threads);
        // std::vector<std::future<int>> results;  // return typename: xxx
        // for (int i = 0; i < Vnum; i++) {
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

    L.resize(Vnum);
    L_temp.resize(Vnum);
    endpos1.resize(Vnum);

    // Line 1-2.
    for (size_t i = 0; i < Vnum; i++) {
        L[i].emplace_back(i, 0, i);
    }

    dmin.resize(num_of_threads);
    dirt.resize(num_of_threads);
    for (int i = 0; i < num_of_threads; i++) {
        dmin[i].resize(Vnum);
        dirt[i].assign(Vnum, -1);
        thread_id.push(i);
    }
    // if_continue_595 = true;
    // int if_continue_595_false_time = 0;

    ThreadPool pool(num_of_threads);
    std::vector<std::future<int>> results;
    // int num_of_threads_per_push = num_of_threads * 1000;
    // ÿ��push��ȥ num_of_threads_per_push
    // �̣߳����û���쳣������push��ȥnum_of_threads_per_push�̣߳����ȫ��һ��push��ȥ����ȫ���̶߳���������catch�쳣

    // since R2 skip some vertices, some new labels can only be
    // generated when d increases 2, not 1, thus terminate the
    // loop only when if_continue_595==false twice

    while (true) {
        // cout << "here" << endl;

        // if (if_continue_595 == true) {
        //     if_continue_595_false_time = 0;
        // }
        // if_continue_595 = false;
        // int push_num = 0;

        // Line 4-20
        for (size_t u = 0; u < Vnum; u++) {
            // if (ideal_graph_595[u].size() == 0 ||
            //     case_info.reduction_measures_2019R2[vertexID_new_to_old_595[u]] == 2)
            //     continue;  // do not search isolated vertices
            // auto *xx = &case_info;

            // new labels add into L_temp[u], but also read L in the process
            results.emplace_back(pool.enqueue([u, this, &g] {
                this->propagate(g, (int)u);
                return 1;
            }));
            // push_num++;
            // if (push_num % num_of_threads_per_push == 0) {
            //     for (auto &&result : results)
            //         result.get();  // all threads finish here
            //     results.clear();
            // }
        }

        for (auto &&result : results)
            result.get();
        results.clear();

        for (size_t u = 0; u < Vnum; u++) {
            // if (ideal_graph_595[u].size() == 0 ||
            //     case_info.reduction_measures_2019R2[vertexID_new_to_old_595[u]] == 2)
            //     continue;                            // do not search isolated vertices

            // new labels in L_temp[u] add into L[u], to avoid locking L in propagate process
            results.emplace_back(pool.enqueue([u, this] {
                this->append((int)u);
                return 1;
            }));
        }

        for (auto &&result : results)
            result.get();
        results.clear();  // 如果本轮没有异常则继续

        bool is_empty = true;
        for (size_t u = 0; u < Vnum; u++) {
            if (!L_temp[u].empty()) {
                is_empty = false;
            }
            endpos1[u] = L[u].size() - L_temp[u].size();
            L_temp[u].clear();
        }
        if (is_empty) {
            break;
        }
        // if (if_continue_595 == false) {
        //     if_continue_595_false_time++;  // if if_continue_595_false_time==2, then
        //                                    // if_continue_595==false twice
        // }
    }

    // remove redundant labels in L
    for (size_t u = 0; u < Vnum; u++) {
        results.emplace_back(pool.enqueue([u, this] {
            this->shrink((int)u);
            return 1;
        }));
    }
    for (auto &&result : results)
        result.get();
    results.clear();

    end_time = std::chrono::high_resolution_clock::now();
    case_info.time_generate_labels =
        (double)std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - begin_time).count() / 1e9;
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
    puts("ok");
}

template <class weight_t>
void PSL<weight_t>::propagate(const dgraph_v_of_v<weight_t> &g, int u) {
    mtx.lock();
    int current_tid = thread_id.front();
    thread_id.pop();
    mtx.unlock();

    // Line 5
    for (const label<weight_t> &cur : L[u]) {
        int v = cur.vertex;
        weight_t d = cur.dist;
        if (dirt[current_tid][v] == -1) {
            dirt[current_tid][v] = 0;
            dmin[current_tid][v] = d;
        } else {
            dmin[current_tid][v] = std::min(dmin[current_tid][v], d);
        }
    }

    for (auto [v, dd] : g[u]) {
        for (size_t i = endpos1[v]; i < L[v].size(); i++) {
            int x = L[v][i].vertex;
            if (rank[x] < rank[u]) {
                continue;
            }
            weight_t d_temp = dd + L[v][i].dist;

            // Line 12-16
            bool flag = false;
            for (size_t j = endpos1[u]; j < L[u].size(); j++) {
                int y = L[u][j].vertex;
                if (!dirt[current_tid][y] && dmin[current_tid][y] + L[u][j].dist <= d_temp) {
                    flag = true;
                }
            }
            for (size_t j = endpos1[x]; j < L[x].size(); j++) {
                int y = L[x][j].vertex;
                if (!dirt[current_tid][y] && dmin[current_tid][y] + L[x][j].dist <= d_temp) {
                    flag = true;
                }
            }
            if (flag) {
                continue;
            }

            // Line 17
            L_temp[u].emplace_back(x, d_temp, v);
        }
    }

    for (const label<weight_t> &cur : L[u]) {
        int v = cur.vertex;
        dirt[current_tid][v] = -1;
    }

    mtx.lock();
    thread_id.push(current_tid);
    mtx.unlock();
}

template <class weight_t>
void PSL<weight_t>::append(int u) {
    L[u].insert(L[u].end(), L_temp[u].begin(), L_temp[u].end());
    L[u].shrink_to_fit();
}

template <class weight_t>
void PSL<weight_t>::shrink(int u) {
    mtx.lock();
    int current_tid = thread_id.front();
    thread_id.pop();
    mtx.unlock();

    for (const label<weight_t> &cur : L[u]) {
        int v = cur.vertex;
        weight_t d = cur.dist;
        if (dirt[current_tid][v] == -1) {
            dmin[current_tid][v] = d;
            dirt[current_tid][v] = cur.prev;
        } else if (d < dmin[current_tid][v]) {
            dmin[current_tid][v] = d;
            dirt[current_tid][v] = cur.prev;
        }
    }

    L[u].clear();
    for (size_t i = 0; i < dirt[current_tid].size(); i++) {
        if (dirt[current_tid][i] == -1) {
            continue;
        }
        L[u].emplace_back(i, dmin[current_tid][i], dirt[current_tid][i]);
        dirt[current_tid][i] = -1;
    }

    mtx.lock();
    thread_id.push(current_tid);
    mtx.unlock();
}
