/**
 * @file config.h
 * @version 0.1
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <iostream>
#include <vector>

namespace PSL {

struct runtime_info {
    /* Running time */
    long long time_initialization = 0;
    long long time_generate_labels = 0;
    // long long time_canonical_repair1 = 0;
    // long long time_canonical_repair2 = 0;

    /* Reduction 1: Equivalence between Vertices */
    bool use_equiv = false;
    long long time_rdc_equiv = 0;

    /* Reduction 2 : Local Minimum Set */
    bool use_lms = false;  // Local minimum set elimination
    long long time_rdc_lms = 0;

    // /*running limits*/
    // long long int max_labal_size = 1e12;  // 2-hop-label num
    // double max_run_time_seconds = 1e12;   // s

    /*labels*/
    // std::vector<int>
    //     reduction_measures_2019R1;  // for 2019 R1;  11 means equivalent_1 relation (no edge
    //                                 // between), 12 means equivalent_2 relation (edge between)
    // std::vector<int> f_2019R1;      // for 2019 R1
    // std::vector<std::vector<two_hop_label_v1>> L;

    // /*canonical_repair info*/
    // bool use_canonical_repair = false;
    // long long int label_size_before_canonical_repair = 0;
    // long long int label_size_after_canonical_repair = 0;
    // double canonical_repair_remove_label_ratio = 0;

    // /*compute label size; this should equal label_size_after_canonical_repair when
    //  * use_canonical_repair==true*/
    // long long int compute_label_bit_size() {
    //     long long int size = 0;
    //     size = size + reduction_measures_2019R2.size() * 4;
    //     size = size + reduction_measures_2019R1.size() * 4;
    //     size = size + f_2019R1.size() * 4;
    //     for (auto it = L.begin(); it != L.end(); it++) {
    //         size = size + (*it).size() * sizeof(two_hop_label_v1);  // 12 bit per
    //         two_hop_label_v1
    //     }
    //     return size;
    // }

    // /*clear labels*/
    // void clear_labels() {
    //     std::vector<int>().swap(reduction_measures_2019R2);
    //     std::vector<int>().swap(reduction_measures_2019R1);
    //     std::vector<int>().swap(f_2019R1);
    //     std::vector<std::vector<two_hop_label_v1>>().swap(L);
    // }

    // /*printing*/
    // void print_L() {
    //     cout << "print_L:" << endl;
    //     for (int i = 0; i < L.size(); i++) {
    //         cout << "L[" << i << "]=";
    //         for (int j = 0; j < L[i].size(); j++) {
    //             cout << "{" << L[i][j].vertex << "," << L[i][j].distance << ","
    //                  << L[i][j].parent_vertex << "}";
    //         }
    //         cout << endl;
    //     }
    // }
    // void print_reduction_measures_2019R1() {
    //     cout << "print_reduction_measures_2019R1:" << endl;
    //     for (int i = 0; i < reduction_measures_2019R1.size(); i++) {
    //         cout << "reduction_measures_2019R1[" << i << "]=" << reduction_measures_2019R1[i]
    //              << endl;
    //     }
    // }
    // void print_reduction_measures_2019R2() {
    //     cout << "print_reduction_measures_2019R2:" << endl;
    //     for (int i = 0; i < reduction_measures_2019R2.size(); i++) {
    //         cout << "reduction_measures_2019R2[" << i << "]=" << reduction_measures_2019R2[i]
    //              << endl;
    //     }
    // }
    // void print_f_2019R1() {
    //     cout << "print_f_2019R1:" << endl;
    //     for (int i = 0; i < f_2019R1.size(); i++) {
    //         cout << "f_2019R1[" << i << "]=" << f_2019R1[i] << endl;
    //     }
    // }
};

}  // namespace PSL
