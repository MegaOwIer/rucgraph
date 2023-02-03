#pragma once

template <typename weight_type>
class d_2hoplabel {
public:
	int vertex; // std::vector<PLL_with_non_adj_reduction_sorted_label> should be sorted by vertex
	weight_type distance;
};

template <typename weight_type>
class d_2hoplabels {
public:
	std::vector<std::vector<d_2hoplabel<weight_type>>> INlabels, OUTlabels;
};

/*
example:
------------------------------------

#include <build_in_progress/HL/dgraph/dHL_base.h>

int main()
{
	d_2hoplabels<float> L; // ≥ı ºªØ
}

----------------------
*/















