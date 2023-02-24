/**
 * @brief Sample testing codes for PSL for Directed graphs.
 */

#include <iostream>
#include <string>

#include "PSL.hpp"
#include "utility/config.h"
#include "utility/dgraph.hpp"

void run_small_sample() {

    PSL::dgraph<int> g(6);

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
    PSL::runtime_info cfg;
    cfg.use_lms = false;

    PSL::PSL<int> solve(g, 6, cfg);

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            std::cout << solve.query_dist(i, j) << "\t";
        }
        std::cout << "\n";
    }
}

int main() {
    run_small_sample();
    return 0;
}
