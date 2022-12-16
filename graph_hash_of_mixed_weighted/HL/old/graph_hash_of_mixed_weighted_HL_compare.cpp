#include <iostream>
#include <fstream>
using namespace std;

// header files in the Boost library: https://www.boost.org/
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };

#include <graph_hash_of_mixed_weighted/HL/PSL_weighted/graph_hash_of_mixed_weighted_HL_PSL_weighted_multiple_threads.h>
#include <graph_hash_of_mixed_weighted/HL/PLL/graph_hash_of_mixed_weighted_HL_PLL_multiple_threads.h>
#include <graph_hash_of_mixed_weighted/HL/PLL_with_non_adj_reduction/graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_multiple_threads.h>


void graph_hash_of_mixed_weighted_HL_PLL_VS_PLL_weighted() {

	/*parameters*/
	int iteration_graph_times = 10;
	int V = 1e4, E = 5e4, precision = 1;
	double ec_min = 1, ec_max = 10; // set ec_min=ec_max=1 for testing unweighted PLL

	vector<int> thread_nums = { 1, 5, 10, 30, 60 };

	for (auto its = thread_nums.begin(); its != thread_nums.end(); its++) {

		double avg_index_size_per_v_PLL = 0, avg_index_size_per_v_PSL = 0;
		double time_per_graph_PLL = 0, time_per_graph_PSL = 0;

		int thread_num = *its;

		/*iteration*/
		for (int i = 0; i < iteration_graph_times; i++) {

			if (1) {
				/*input and output; below is for generating random new graph, or read saved graph*/
				int generate_new_graph = 1;
				std::unordered_set<int> generated_group_vertices;
				graph_hash_of_mixed_weighted instance_graph, generated_group_graph;
				if (generate_new_graph == 1) {
					instance_graph = graph_hash_of_mixed_weighted_generate_random_graph(V, E, 0, 0, 1, 1, 1, boost_random_time_seed);
					graph_hash_of_mixed_weighted_save_graph_with_weight("simple_iterative_tests.txt", instance_graph, 0);
				}
				else {
					double lambda;
					graph_hash_of_mixed_weighted_read_graph_with_weight("simple_iterative_tests.txt", instance_graph, lambda);
				}

				if (1) {
					auto begin = std::chrono::high_resolution_clock::now();
					vector<vector<graph_hash_of_mixed_weighted_HL_PLL_label>> L = graph_hash_of_mixed_weighted_HL_PLL_generate_indexes_multiple_threads(instance_graph, V + 1, true, thread_num);
					auto end = std::chrono::high_resolution_clock::now();
					double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
					time_per_graph_PLL = time_per_graph_PLL + runningtime / iteration_graph_times;

					long long int index_size = 0;
					for (auto it = L.begin(); it != L.end(); it++) {
						index_size = index_size + (*it).size();
					}
					avg_index_size_per_v_PLL = avg_index_size_per_v_PLL + (double)index_size / V / iteration_graph_times;
				}

				if (1) {
					auto begin = std::chrono::high_resolution_clock::now();
					vector<vector<graph_hash_of_mixed_weighted_HL_PSL_weighted_label>> L = graph_hash_of_mixed_weighted_HL_PSL_weighted_generate_indexes_multiple_threads(instance_graph, V + 1, thread_num);
					auto end = std::chrono::high_resolution_clock::now();
					double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
					time_per_graph_PSL = time_per_graph_PSL + runningtime / iteration_graph_times;

					long long int index_size = 0;
					for (auto it = L.begin(); it != L.end(); it++) {
						index_size = index_size + (*it).size();
					}
					avg_index_size_per_v_PSL = avg_index_size_per_v_PSL + (double)index_size / V / iteration_graph_times;
				}
			}
		}

		cout << "thread_num: " << thread_num << endl;
		cout << "time_per_graph_PLL: " << time_per_graph_PLL << " avg_index_size_per_v_PLL: " << avg_index_size_per_v_PLL << endl;
		cout << "time_per_graph_PSL: " << time_per_graph_PSL << " avg_index_size_per_v_PSL: " << avg_index_size_per_v_PSL << endl << endl;

	}

}


void graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_compare() {

	/*parameters*/
	int iteration_graph_times = 10;
	int V = 1e4, E = 5e4, precision = 1;
	double ec_min = 1, ec_max = 10; // set ec_min=ec_max=1 for testing unweighted PLL
	int thread_num = 20;


	double avg_index_size_per_v_PLL_2019R2 = 0, avg_index_size_per_v_PLL_2019R2enhance_15 = 0;
	double time_per_graph_PLL_2019R2 = 0, time_per_graph_PLL_2019R2enhance_15 = 0;

	/*iteration*/
	for (int i = 0; i < iteration_graph_times; i++) {

		if (1) {
			/*input and output; below is for generating random new graph, or read saved graph*/
			int generate_new_graph = 1;
			std::unordered_set<int> generated_group_vertices;
			graph_hash_of_mixed_weighted instance_graph, generated_group_graph;
			if (generate_new_graph == 1) {
				instance_graph = graph_hash_of_mixed_weighted_generate_random_graph(V, E, 0, 0, 1, 1, 1, boost_random_time_seed);
				graph_hash_of_mixed_weighted_save_graph_with_weight("simple_iterative_tests.txt", instance_graph, 0);
			}
			else {
				double lambda;
				graph_hash_of_mixed_weighted_read_graph_with_weight("simple_iterative_tests.txt", instance_graph, lambda);
			}

			if (1) {
				auto begin = std::chrono::high_resolution_clock::now();
				vector<vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>> L =
					graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_generate_indexes_multiple_threads(instance_graph, V + 1, true, thread_num, "2019R2");
				auto end = std::chrono::high_resolution_clock::now();
				double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
				time_per_graph_PLL_2019R2 = time_per_graph_PLL_2019R2 + runningtime / iteration_graph_times;

				long long int index_size = 0;
				for (auto it = L.begin(); it != L.end(); it++) {
					index_size = index_size + (*it).size();
				}
				avg_index_size_per_v_PLL_2019R2 = avg_index_size_per_v_PLL_2019R2 + (double)index_size / V / iteration_graph_times;
			}

			if (1) {
				auto begin = std::chrono::high_resolution_clock::now();
				vector<vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>> L =
					graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_generate_indexes_multiple_threads(instance_graph, V + 1, true, thread_num, "2019R2enhance_15");
				auto end = std::chrono::high_resolution_clock::now();
				double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
				time_per_graph_PLL_2019R2enhance_15 = time_per_graph_PLL_2019R2enhance_15 + runningtime / iteration_graph_times;

				long long int index_size = 0;
				for (auto it = L.begin(); it != L.end(); it++) {
					index_size = index_size + (*it).size();
				}
				avg_index_size_per_v_PLL_2019R2enhance_15 = avg_index_size_per_v_PLL_2019R2enhance_15 + (double)index_size / V / iteration_graph_times;
			}
		}
	}

	cout << "thread_num: " << thread_num << endl;
	cout << "time_per_graph_PLL_2019R2: " << time_per_graph_PLL_2019R2 << " avg_index_size_per_v_PLL_2019R2: " << avg_index_size_per_v_PLL_2019R2 << endl;
	cout << "time_per_graph_PLL_2019R2enhance_15: " << time_per_graph_PLL_2019R2enhance_15 << " avg_index_size_per_v_PLL_2019R2enhance_15: " << avg_index_size_per_v_PLL_2019R2enhance_15 << endl << endl;



}


int main()
{
	graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_compare();
}