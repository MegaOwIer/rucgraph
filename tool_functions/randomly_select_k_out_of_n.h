#pragma once



#include <copy_items.h>

vector<int> randomly_select_k_out_of_n_vector_int(vector<int>& input_vector, int target_num) {

	vector<int> candidates = copy_vector_int(input_vector);

	vector<int> output_vector;

	std::time_t now = std::time(0);
	boost::random::mt19937 gen{ static_cast<std::uint32_t>(now) };

	while (output_vector.size() < target_num) {
		boost::random::uniform_int_distribution<> dist{ 0, int(candidates.size() - 1) }; // generating random number in [0, candidates.size()-1]
		int rand = dist(gen);
		output_vector.push_back(candidates[rand]);
		candidates.erase(candidates.begin() + rand);
	}

	return output_vector;
}