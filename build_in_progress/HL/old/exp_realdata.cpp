#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
using namespace std;

// header files in the Boost library: https://www.boost.org/
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };

#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_ec_update_pairwise_jaccard_distance.h>
#include <graph_hash_of_mixed_weighted/HL/graph_hash_of_mixed_weighted_ec_update_2007_DPBF_distance.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_binary_save_read.h>
#include <graph_hash_of_mixed_weighted/HL/PLL_with_non_adj_reduction/graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_multiple_threads.h>
#include <graph_hash_of_mixed_weighted/HL/PSL_with_non_adj_reduction/graph_hash_of_mixed_weighted_HL_PSL_with_non_adj_reduction_multiple_threads.h>
#include <graph_hash_of_mixed_weighted/HL/PSL_enhancedoriginalR2/graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_multiple_threads.h>
#include <ThreadPool.h>

/*read txt file to binary graph data*/
void generate_binary_graph() {

	vector<pair<string, string>> files_names;
	files_names.push_back({ "/home/dengs/HLData/TOPC/wiki-topcats.txt",  "TOPC" });
	files_names.push_back({ "/home/dengs/HLData/DBLP/dblp_coauthor/dblp_coauthor.txt",  "DBLP" });
	files_names.push_back({ "/home/dengs/HLData/FLIC/flickr-growth/flickr_growth.txt",  "FLIC" });
	files_names.push_back({ "/home/dengs/HLData/Douban/douban/douban.txt",  "Douban" });

	int max_N = 1e7;

	srand((int)time(0));
	for (int i = files_names.size() - 1; i >= 0; i--) {
		ifstream myfile(files_names[i].first);
		string name = files_names[i].second;
		if (!myfile.is_open())
		{
			cout << "can't open this file" << endl;
			exit(1);
		}
		graph_hash_of_mixed_weighted input_graph;

		int v1, v2;
		while (1)
		{
			myfile >> v1 >> v2;
			if (v1 == -1 && v2 == -1)
				break;
			float w = (rand() / float(RAND_MAX)) + 1e-4;
			graph_hash_of_mixed_weighted_add_edge(input_graph, v1, v2, w);

			v1 = -1;
			v2 = -1;
		}

		cout << files_names[i].second << endl;

		graph_hash_of_mixed_weighted_print_size(input_graph);

		graph_hash_of_mixed_weighted_binary_save(input_graph, files_names[i].second + "_random.bin");

		graph_hash_of_mixed_weighted_ec_update_pairwise_jaccard_distance_fast(input_graph, max_N);

		graph_hash_of_mixed_weighted_binary_save(input_graph, files_names[i].second + "_jaccard.bin");

		graph_hash_of_mixed_weighted_ec_update_2007_DPBF_distance(input_graph, max_N);

		graph_hash_of_mixed_weighted_binary_save(input_graph, files_names[i].second + "_2007_DPBF.bin");

	}
}


template <typename T>
long long int HL_label_num(vector<vector<T>>& L) {

	long long int num = 0;
	for (int i = L.size() - 1; i >= 0; i--) {
		num = num + L[i].size();
	}

	return num;
}

/*query dis path*/

pair<double, double> compute_PLL_query_dis_path_average_time_parallel_elelment(vector<vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>>* L, graph_hash_of_mixed_weighted* instance_graph, int s, int t) {

	pair<double, double> result;

	if (1) {
		auto begin = std::chrono::high_resolution_clock::now();
		graph_hash_of_mixed_weighted_HL_PLL_with_2019R1_extract_distance
		(*L, graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1,
			graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_2, graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_f_2019R1, *instance_graph, s, t);
		auto end = std::chrono::high_resolution_clock::now();
		result.first = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	}
	if (1) {
		auto begin = std::chrono::high_resolution_clock::now();
		graph_hash_of_mixed_weighted_HL_PLL_with_2019R1_extract_shortest_path
		(*L, graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1,
			graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_2, graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_f_2019R1, *instance_graph, s, t);
		auto end = std::chrono::high_resolution_clock::now();
		result.second = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	}

	return result;

}

pair<double, double> compute_PLL_query_dis_path_average_time(vector<vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>>& L,
	int query_times, int max_V_ID, graph_hash_of_mixed_weighted& instance_graph, int num_of_threads) {

	vector<pair<int, int>> s_t;
	boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(max_V_ID) };
	for (int x = 0; x < query_times; x++) {
		s_t.push_back({ dist(boost_random_time_seed), dist(boost_random_time_seed) });
	}

	vector<vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>>* p1 = &L;
	graph_hash_of_mixed_weighted* p2 = &instance_graph;

	ThreadPool pool(num_of_threads);
	std::vector< std::future<pair<double, double>> > results; // return typename: xxx
	for (int x = 0; x < query_times; x++) {
		int s = s_t[x].first, t = s_t[x].second;
		results.emplace_back(
			pool.enqueue([p1, p2, s, t] { // pass const type value j to thread; [] can be empty
			return compute_PLL_query_dis_path_average_time_parallel_elelment(p1, p2, s, t);
		})
		);
	}

	pair<double, double> final_result = { 0, 0 };
	for (auto&& result : results) {
		pair<double, double> x = result.get(); // all threads finish here
		final_result.first = final_result.first + x.first / query_times; // s
		final_result.second = final_result.second + x.second / query_times; // s
	}

	return final_result;

}


pair<double, double> compute_PSL_query_dis_path_average_time_parallel_elelment(vector<vector<graph_hash_of_mixed_weighted_HL_PSL_with_non_adj_reduction_label>>* L, graph_hash_of_mixed_weighted* instance_graph, int s, int t) {

	pair<double, double> result;

	if (1) {
		auto begin = std::chrono::high_resolution_clock::now();
		graph_hash_of_mixed_weighted_HL_PSL_with_2019R1_extract_distance
		(*L, graph_hash_of_mixed_weighted_HL_PSL_with_non_adj_reduction_reduction_measures_1,
			graph_hash_of_mixed_weighted_HL_PSL_with_non_adj_reduction_reduction_measures_2, graph_hash_of_mixed_weighted_HL_PSL_with_non_adj_reduction_f_2019R1, *instance_graph, s, t);
		auto end = std::chrono::high_resolution_clock::now();
		result.first = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	}
	if (1) {
		auto begin = std::chrono::high_resolution_clock::now();
		graph_hash_of_mixed_weighted_HL_PSL_with_2019R1_extract_shortest_path
		(*L, graph_hash_of_mixed_weighted_HL_PSL_with_non_adj_reduction_reduction_measures_1,
			graph_hash_of_mixed_weighted_HL_PSL_with_non_adj_reduction_reduction_measures_2, graph_hash_of_mixed_weighted_HL_PSL_with_non_adj_reduction_f_2019R1, *instance_graph, s, t);
		auto end = std::chrono::high_resolution_clock::now();
		result.second = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	}

	return result;

}

pair<double, double> compute_PSL_query_dis_path_average_time(vector<vector<graph_hash_of_mixed_weighted_HL_PSL_with_non_adj_reduction_label>>& L,
	int query_times, int max_V_ID, graph_hash_of_mixed_weighted& instance_graph, int num_of_threads) {

	vector<pair<int, int>> s_t;
	boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(max_V_ID) };
	for (int x = 0; x < query_times; x++) {
		s_t.push_back({ dist(boost_random_time_seed), dist(boost_random_time_seed) });
	}

	vector<vector<graph_hash_of_mixed_weighted_HL_PSL_with_non_adj_reduction_label>>* p1 = &L;
	graph_hash_of_mixed_weighted* p2 = &instance_graph;

	ThreadPool pool(num_of_threads);
	std::vector< std::future<pair<double, double>> > results; // return typename: xxx
	for (int x = 0; x < query_times; x++) {
		int s = s_t[x].first, t = s_t[x].second;
		results.emplace_back(
			pool.enqueue([p1, p2, s, t] { // pass const type value j to thread; [] can be empty
			return compute_PSL_query_dis_path_average_time_parallel_elelment(p1, p2, s, t);
		})
		);
	}

	pair<double, double> final_result = { 0, 0 };
	for (auto&& result : results) {
		pair<double, double> x = result.get(); // all threads finish here
		final_result.first = final_result.first + x.first / query_times; // s
		final_result.second = final_result.second + x.second / query_times; // s
	}

	return final_result;

}


pair<double, double> compute_PSL_enhancedoriginalR2_query_dis_path_average_time_parallel_elelment(vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>>* L, graph_hash_of_mixed_weighted* instance_graph, int s, int t) {

	pair<double, double> result;

	if (1) {
		auto begin = std::chrono::high_resolution_clock::now();
		graph_hash_of_mixed_weighted_HL_PSL_with_2019R1_extract_distance
		(*L, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1,
			graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, *instance_graph, s, t);
		auto end = std::chrono::high_resolution_clock::now();
		result.first = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	}
	if (1) {
		auto begin = std::chrono::high_resolution_clock::now();
		graph_hash_of_mixed_weighted_HL_PSL_with_2019R1_extract_shortest_path
		(*L, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1,
			graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, *instance_graph, s, t);
		auto end = std::chrono::high_resolution_clock::now();
		result.second = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	}

	return result;

}


pair<double, double> compute_PSL_enhancedoriginalR2_query_dis_path_average_time(vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>>& L,
	int query_times, int max_V_ID, graph_hash_of_mixed_weighted& instance_graph, int num_of_threads) {

	vector<pair<int, int>> s_t;
	boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(max_V_ID) };
	for (int x = 0; x < query_times; x++) {
		s_t.push_back({ dist(boost_random_time_seed), dist(boost_random_time_seed) });
	}

	vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>>* p1 = &L;
	graph_hash_of_mixed_weighted* p2 = &instance_graph;

	ThreadPool pool(num_of_threads);
	std::vector< std::future<pair<double, double>> > results; // return typename: xxx
	for (int x = 0; x < query_times; x++) {
		int s = s_t[x].first, t = s_t[x].second;
		results.emplace_back(
			pool.enqueue([p1, p2, s, t] { // pass const type value j to thread; [] can be empty
			return compute_PSL_enhancedoriginalR2_query_dis_path_average_time_parallel_elelment(p1, p2, s, t);
		})
		);
	}

	pair<double, double> final_result = { 0, 0 };
	for (auto&& result : results) {
		pair<double, double> x = result.get(); // all threads finish here
		final_result.first = final_result.first + x.first / query_times; // s
		final_result.second = final_result.second + x.second / query_times; // s
	}

	return final_result;

}



/*experiments*/
void experiment_element_PLL_element(ofstream& outputFile, string& file_path, graph_hash_of_mixed_weighted& instance_graph,
	int thread_num, int max_degree_MG_enhanced2019R2, int query_times, graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_Use_Reduction_info& mm, string algorithm_name) {

	bool catch_error = false;
	vector<vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>> L;
	auto begin = std::chrono::high_resolution_clock::now();
	try {
		L = graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_generate_indexes_multiple_threads(instance_graph, instance_graph.hash_of_vectors.size(), true, thread_num, mm);
	}
	catch (string s) {
		cout << s << endl;
		terminate_procedures_592();
		catch_error = true;
	}
	auto end = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	if (catch_error == false) {
		pair<double, double> query_t = compute_PLL_query_dis_path_average_time(L, query_times, instance_graph.hash_of_vectors.size() - 1, instance_graph, thread_num);
		outputFile << file_path << "," + algorithm_name + "," << thread_num << "," << max_degree_MG_enhanced2019R2 << "," << runningtime << "," << mm.time_2019R1 << "," << mm.time_2019R2_or_enhanced_pre << ","
			<< mm.time_2019R2_or_enhanced_fixlabels << "," << HL_label_num(L) << "," << query_t.first << "," << query_t.second << "," << mm.reduce_V_num_2019R1
			<< "," << mm.MG_num << endl;
	}
	else {
		outputFile << file_path << "," + algorithm_name + "," << thread_num << "," << max_degree_MG_enhanced2019R2 << "," << -1 << "," << -1 << "," << -1 << ","
			<< -1 << "," << -1 << "," << -1 << "," << -1 << "," << -1 << "," << -1 << endl;
	}

	graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_clear_labels(L);

}

void experiment_element_PSL_New_edges_element(ofstream& outputFile, string& file_path, graph_hash_of_mixed_weighted& instance_graph,
	int thread_num, int max_degree_MG_enhanced2019R2, int query_times, graph_hash_of_mixed_weighted_HL_PSL_with_non_adj_reduction_Use_Reduction_info& mm, string algorithm_name) {

	bool catch_error = false;
	vector<vector<graph_hash_of_mixed_weighted_HL_PSL_with_non_adj_reduction_label>> L;
	auto begin = std::chrono::high_resolution_clock::now();
	try {
		L = graph_hash_of_mixed_weighted_HL_PSL_with_non_adj_reduction_generate_indexes_multiple_threads(instance_graph, instance_graph.hash_of_vectors.size(), thread_num, mm);
	}
	catch (string s) {
		cout << s << endl;
		terminate_procedures_173();
		catch_error = true;
	}
	auto end = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	if (catch_error == false) {
		pair<double, double> query_t = compute_PSL_query_dis_path_average_time(L, query_times, instance_graph.hash_of_vectors.size() - 1, instance_graph, thread_num);
		outputFile << file_path << "," + algorithm_name + "," << thread_num << "," << max_degree_MG_enhanced2019R2 << "," << runningtime << "," << mm.time_2019R1 << "," << mm.time_2019R2_or_enhanced_pre << ","
			<< mm.time_2019R2_or_enhanced_fixlabels << "," << HL_label_num(L) << "," << query_t.first << "," << query_t.second << "," << mm.reduce_V_num_2019R1
			<< "," << mm.MG_num << endl;
	}
	else {
		outputFile << file_path << "," + algorithm_name + "," << thread_num << "," << max_degree_MG_enhanced2019R2 << "," << -1 << "," << -1 << "," << -1 << ","
			<< -1 << "," << -1 << "," << -1 << "," << -1 << "," << -1 << "," << -1 << endl;
	}

	graph_hash_of_mixed_weighted_HL_PSL_with_non_adj_reduction_clear_labels(L);

}


void experiment_element_PSL_enhancedoriginalR2_element(ofstream& outputFile, string& file_path, graph_hash_of_mixed_weighted& instance_graph,
	int thread_num, int max_degree_MG_enhanced2019R2, int query_times, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_Use_Reduction_info& mm, string algorithm_name) {

	bool catch_error = false;
	vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>> L;
	auto begin = std::chrono::high_resolution_clock::now();
	try {
		L = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_generate_indexes_multiple_threads(instance_graph, instance_graph.hash_of_vectors.size(), thread_num, mm);
	}
	catch (string s) {
		cout << s << endl;
		terminate_procedures_316();
		catch_error = true;
	}
	auto end = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	if (catch_error == false) {
		pair<double, double> query_t = compute_PSL_enhancedoriginalR2_query_dis_path_average_time(L, query_times, instance_graph.hash_of_vectors.size() - 1, instance_graph, thread_num);
		outputFile << file_path << "," + algorithm_name + "," << thread_num << "," << max_degree_MG_enhanced2019R2 << "," << runningtime << "," << mm.time_2019R1 << "," << mm.time_2019R2_or_enhanced_pre << ","
			<< 0 << "," << HL_label_num(L) << "," << query_t.first << "," << query_t.second << "," << mm.reduce_V_num_2019R1
			<< "," << mm.MG_num << endl;
	}
	else {
		outputFile << file_path << "," + algorithm_name + "," << thread_num << "," << max_degree_MG_enhanced2019R2 << "," << -1 << "," << -1 << "," << -1 << ","
			<< -1 << "," << -1 << "," << -1 << "," << -1 << "," << -1 << "," << -1 << endl;
	}

	graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_clear_labels(L);

}

void experiment_element(string data_name, string save_file_name) {

	int thread_num = 20;
	int max_degree_MG_enhanced2019R2 = 100;
	int query_times = 1e4; // this cannot be too large
	long long int max_labal_size = ((long long int)100 * 1024 * 1024 * 1024) / 12; // 100 GB max; 12 bit per label
	double max_run_time_seconds = 3600 * 12.0; // 12h

	/*output*/
	ofstream outputFile;
	outputFile.precision(8);
	outputFile.setf(ios::fixed);
	outputFile.setf(ios::showpoint);
	outputFile.open(save_file_name);
	outputFile << "graph_name,method_name,thread_num,max_degree_MG_enhanced2019R2,index_time(s),R1_time(s),time_2019R2_or_enhanced_pre(s),time_2019R2_or_enhanced_fixlabels(s)," <<
		"index_size,query_dis_average_time(s),query_path_average_time(s),R1_removeV_num,R2MG_removeV_num" << endl;

	/*default experiments*/
	for (int j = 0; j < 3; j++) {
		string file_path = "/home/data/HL_data/" + data_name;
		if (j == 0) {
			file_path = file_path + "_jaccard.bin";
		}
		else if (j == 1) {
			//continue;
			file_path = file_path + "_2007_DPBF.bin";
		}
		else if (j == 2) {
			file_path = file_path + "_random.bin";
		}

		cout << "indexing " + file_path + "..." << endl;

		graph_hash_of_mixed_weighted instance_graph = graph_hash_of_mixed_weighted_binary_read(file_path);

		/*PLL*/
		if (1) {
			graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_Use_Reduction_info mm;
			mm.use_2019R1 = 0;
			mm.use_2019R2 = 0;
			mm.use_enhanced2019R2 = 0;
			mm.max_degree_MG_enhanced2019R2 = max_degree_MG_enhanced2019R2;
			mm.max_labal_size = max_labal_size;
			mm.max_run_time_seconds = max_run_time_seconds;

			experiment_element_PLL_element(outputFile, file_path, instance_graph, thread_num, max_degree_MG_enhanced2019R2, query_times, mm, "PLL");
		}

		/*PLL with 2019 reductions*/
		if (1) {
			graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_Use_Reduction_info mm;
			mm.use_2019R1 = 1;
			mm.use_2019R2 = 1;
			mm.use_enhanced2019R2 = 0;
			mm.use_non_adj_reduc_degree = 0;
			mm.max_degree_MG_enhanced2019R2 = max_degree_MG_enhanced2019R2;
			mm.max_labal_size = max_labal_size;
			mm.max_run_time_seconds = max_run_time_seconds;

			experiment_element_PLL_element(outputFile, file_path, instance_graph, thread_num, max_degree_MG_enhanced2019R2, query_times, mm, "PLL_withR1R2");
		}

		/*PLL with enhanced reductions*/
		if (1) {
			graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_Use_Reduction_info mm;
			mm.use_2019R1 = 1;
			mm.use_2019R2 = 0;
			mm.use_enhanced2019R2 = 1;
			mm.use_non_adj_reduc_degree = 0;		
			mm.max_degree_MG_enhanced2019R2 = max_degree_MG_enhanced2019R2;
			mm.max_labal_size = max_labal_size;
			mm.max_run_time_seconds = max_run_time_seconds;

			experiment_element_PLL_element(outputFile, file_path, instance_graph, thread_num, max_degree_MG_enhanced2019R2, query_times, mm, "PLL_withR1EnhancedR2");
		}

		/*PLL with use_non_adj_reduc_degree*/
		if (1) {
			graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_Use_Reduction_info mm;
			mm.use_2019R1 = 1;
			mm.use_2019R2 = 0;
			mm.use_enhanced2019R2 = 0;
			mm.use_non_adj_reduc_degree = 1;
			mm.max_degree_MG_enhanced2019R2 = max_degree_MG_enhanced2019R2;
			mm.max_labal_size = max_labal_size;
			mm.max_run_time_seconds = max_run_time_seconds;

			experiment_element_PLL_element(outputFile, file_path, instance_graph, thread_num, max_degree_MG_enhanced2019R2, query_times, mm, "PLL_withR1use_non_adj_reduc_degree");
		}


		/*PSL with 2019 reductions with new edges*/
		if (1) {
			graph_hash_of_mixed_weighted_HL_PSL_with_non_adj_reduction_Use_Reduction_info mm;
			mm.use_2019R1 = 1;
			mm.use_2019R2 = 1;
			mm.use_enhanced2019R2 = 0;
			mm.max_degree_MG_enhanced2019R2 = max_degree_MG_enhanced2019R2;
			mm.max_labal_size = max_labal_size;
			mm.max_run_time_seconds = max_run_time_seconds;

			experiment_element_PSL_New_edges_element(outputFile, file_path, instance_graph, thread_num, max_degree_MG_enhanced2019R2, query_times, mm, "PSL_withR1R2_withNewEdges");
		}


		/*PSL with 2019 reductions*/
		if (1) {
			graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_Use_Reduction_info mm;
			mm.use_2019R1 = 1;
			mm.use_2019R2 = 1;
			mm.use_enhanced2019R2 = 0;
			mm.use_non_adj_reduc_degree = 0;
			mm.max_degree_MG_enhanced2019R2 = max_degree_MG_enhanced2019R2;
			mm.max_labal_size = max_labal_size;
			mm.max_run_time_seconds = max_run_time_seconds;

			experiment_element_PSL_enhancedoriginalR2_element(outputFile, file_path, instance_graph, thread_num, max_degree_MG_enhanced2019R2, query_times, mm, "PSL_withR1R2");
		}

		/*PSL with enhanced reductions*/
		if (1) {
			graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_Use_Reduction_info mm;
			mm.use_2019R1 = 1;
			mm.use_2019R2 = 0;
			mm.use_enhanced2019R2 = 1;
			mm.use_non_adj_reduc_degree = 0;
			mm.max_degree_MG_enhanced2019R2 = max_degree_MG_enhanced2019R2;
			mm.max_labal_size = max_labal_size;
			mm.max_run_time_seconds = max_run_time_seconds;

			experiment_element_PSL_enhancedoriginalR2_element(outputFile, file_path, instance_graph, thread_num, max_degree_MG_enhanced2019R2, query_times, mm, "PSL_withR1EnhancedR2");
		}

		/*PSL with use_non_adj_reduc_degree*/
		if (1) {
			graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_Use_Reduction_info mm;
			mm.use_2019R1 = 1;
			mm.use_2019R2 = 0;
			mm.use_enhanced2019R2 = 0;
			mm.use_non_adj_reduc_degree = 1;
			mm.max_degree_MG_enhanced2019R2 = max_degree_MG_enhanced2019R2;
			mm.max_labal_size = max_labal_size;
			mm.max_run_time_seconds = max_run_time_seconds;

			experiment_element_PSL_enhancedoriginalR2_element(outputFile, file_path, instance_graph, thread_num, max_degree_MG_enhanced2019R2, query_times, mm, "PSL_withR1use_non_adj_reduc_degree");
		}

	}





}

void experiments() {

	vector<string> data_names = { "Douban" ,"DBLP", "FLIC", "TOPC" };

	int i = 0;
	experiment_element(data_names[i], "experiment_result_" + data_names[i] + ".csv");  // 里面的HL算法不能并行，只能多次编译独立运行
}



int main()
{
	experiments();
	return 0;
}







