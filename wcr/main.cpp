/**
 * @brief Sample testing codes for PSL for Directed graphs.
 */

#include <iostream>
#include <string>

#include "PSL.hpp"
#include "utility/config.h"
#include "utility/dgraph.hpp"

void run_small_sample() {

    dgraph<int> g(6);

    g.add_edge(0, 1, 1);
    g.add_edge(0, 4, 2);
    g.add_edge(0, 5, 4);
    g.add_edge(1, 2, 4);
    g.add_edge(1, 4, 2);
    g.add_edge(2, 3, 3);
    g.add_edge(2, 5, 1);
    g.add_edge(4, 3, 6);
    g.add_edge(4, 5, 1);
    g.add_edge(5, 3, 1);

    /* reduction method selection*/
    PSL_runtime_info mm;
    // mm.use_2019R1 = 0;
    // mm.use_2019R2 = 0;
    // mm.use_enhanced2019R2 = 0;
    // mm.use_non_adj_reduc_degree = 0;
    // mm.max_degree_MG_enhanced2019R2 = 100;
    // mm.max_labal_size = 6e9;
    // mm.max_run_time_seconds = 1e9;
    // mm.use_canonical_repair = true;

    PSL<int> solve(g, 1, mm);

    for(int i = 0; i < 6; i++) {
        for(int j = 0; j < 6; j++) {
            std::cout << solve.query_dist(i, j) << "\t";
        }
        std::cout << "\n";
    }
}

int main() {
    run_small_sample();
    return 0;
}
