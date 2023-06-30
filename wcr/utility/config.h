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
    long long time_total = 0;
    long long time_initialization = 0;
    long long time_generate_labels = 0;
    // long long time_canonical_repair1 = 0;
    // long long time_canonical_repair2 = 0;

    /* Memory Use */
    size_t label_size_before = 0;
    size_t label_size_after = 0;

    /* Reduction 1: Equivalence between Vertices */
    bool use_equiv = false;
    long long time_rdc_equiv = 0;

    /* Reduction 2 : Local Minimum Set */
    bool use_lms = false;  // Local minimum set elimination
    long long time_rdc_lms = 0;

    // /*running limits*/
    // long long int max_labal_size = 1e12;  // 2-hop-label num
    // double max_run_time_seconds = 1e12;   // s

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

    void display(std::ostream &os = std::cout);
};

void runtime_info::display(std::ostream &os) {
    os << "time: " << time_total << "(" << time_initialization << " + " << time_generate_labels
       << ")" << std::endl;
    if (use_lms) {
        os << "time_lms: " << time_rdc_lms << std::endl;
    }
    os << "size: " << label_size_before << " -> " << label_size_after << std::endl;
}

}  // namespace PSL
