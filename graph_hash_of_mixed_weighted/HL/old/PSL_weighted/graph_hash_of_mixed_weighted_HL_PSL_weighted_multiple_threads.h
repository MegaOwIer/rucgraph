/*
 Notes
 1. We can try to figure out a better way to handle the query process.
 2. PSL applied on weighted graphs causes redundancy in labels,
	the same hub may even have different distances in the labels with the same vertex.
 3. The key to the reduction of labels may be finding a more suitable shortest path algorithm when generating the indexes.
 4. The altenative way is reserved in the explanatory notes. Do not delete them. I'm not sure whether it would perform better.
 5. The function "PSL_generate_indexes_vectorFORMAT_para" is called in the file "graph_hash_of_mixed_weighted_PSL_weighted.h" by the function "PSL_weighted_generate_and_save_indexes_vectorFORMAT_binary", which is called by "test_PSL_weighted_vectorFORMAT".
	 And I call the function "test_PSL_weighted_vectorFORMAT" in the file "ConsoleApplication1.cpp" to go through some random tests.
 */

#pragma once
#include <iostream>
#include <ThreadPool.h>
#include <mutex>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID.h>
#include <graph_hash_of_mixed_weighted/HL/PSL_weighted/graph_hash_of_mixed_weighted_HL_PSL_weighted_labels.h>

 /*unique code for this file: 579*/

static std::mutex mtx_579;
static graph_v_of_v_idealID ideal_graph_579;
static vector<vector<graph_hash_of_mixed_weighted_HL_PSL_weighted_label>> L_579;
static vector<vector<graph_hash_of_mixed_weighted_HL_PSL_weighted_label>> L_temp_579;
static vector<int> pos_579;
static vector<int> increment_579;
static vector<vector<float>> T_579;
static vector<vector<bool>> dirty_tag_579;
static queue<int> Qid_579; // IDs of available elements for T
static bool if_continue_579;



bool graph_hash_of_mixed_weighted_HL_PSL_weighted_compare(const pair<int, int>& i, pair<int, int>& j)
{
	/*< is nearly 10 times slower than >*/
	return i.second > j.second;  // < is from small to big; > is from big to small.  sort by the second item of pair<int, int>
}


void graph_hash_of_mixed_weighted_HL_PSL_weighted_thread_function_dij(int u)
{
	/*the id of T and dirty_flag uesd here */
	int used_id;
	mtx_579.lock();
	used_id = Qid_579.front();
	Qid_579.pop();
	mtx_579.unlock();

	/* save labels of u in T */
	int u_label_size = L_579[u].size();
	for (int i = 0; i < u_label_size; i++)
	{
		int w = L_579[u][i].vertex; float dis = L_579[u][i].distance;
		if (dirty_tag_579[used_id][w]) // the same hub may have redundancy, record the shortest distance
		{
			dirty_tag_579[used_id][w] = false;
			T_579[used_id][w] = dis;
		}
		else if (dis < T_579[used_id][w]) T_579[used_id][w] = dis;
	}


	int u_adj_size = ideal_graph_579[u].size();
	for (int i = 0; i < u_adj_size; i++)
	{

		int v = ideal_graph_579[u][i].first; // visit each neighbour of u
		float ec = ideal_graph_579[u][i].second; // length of edge(u,v)
		int v_label_size = pos_579[v] + increment_579[v]; //label size in the last iteration
		for (int j = pos_579[v]; j < v_label_size; j++)
		{
			int w = L_579[v][j].vertex; // visit each hub of v
			if (w > u) continue; // cannot be the hub of u

			/* query pruning */
			float dis = L_579[v][j].distance + ec;
			int w_label_size = L_579[w].size();
			bool flag = false;
			for (int k = 0; k < w_label_size; k++)
				if (!dirty_tag_579[used_id][L_579[w][k].vertex])
				{
					float query_dis = T_579[used_id][L_579[w][k].vertex] + L_579[w][k].distance;
					if (query_dis - 1e-6 <= dis) { flag = true; break; }
				}
			if (flag) continue;

			/*add a new label*/
			graph_hash_of_mixed_weighted_HL_PSL_weighted_label xx;
			xx.vertex = w;
			xx.distance = dis;
			xx.parent_vertex = v;
			L_temp_579[u].push_back(xx);
			if (!if_continue_579) if_continue_579 = true;
		}
	}

	for (int i = 0; i < u_label_size; i++) dirty_tag_579[used_id][L_579[u][i].vertex] = true; // recover dirty tag

	mtx_579.lock();
	Qid_579.push(used_id);
	mtx_579.unlock();

}


void graph_hash_of_mixed_weighted_HL_PSL_weighted_thread_function_add_new_labels(int u)
{
	pos_579[u] += increment_579[u]; increment_579[u] = 0;
	for (vector<graph_hash_of_mixed_weighted_HL_PSL_weighted_label>::iterator it = L_temp_579[u].begin(); it != L_temp_579[u].end();)
	{
		L_579[u].push_back(*it);
		L_temp_579[u].erase(it);
		increment_579[u]++;
	}
}


vector<vector<graph_hash_of_mixed_weighted_HL_PSL_weighted_label>> graph_hash_of_mixed_weighted_HL_PSL_weighted_generate_indexes_multiple_threads(graph_hash_of_mixed_weighted& input_graph, int max_N, int num_of_threads) {

	int N = input_graph.hash_of_vectors.size();

	/*sort vertices by degrees*/
	vector<pair<int, int>> sorted_vertices;
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		sorted_vertices.push_back({ it->first, input_graph.degree(it->first) });
	}
	sort(sorted_vertices.begin(), sorted_vertices.end(), graph_hash_of_mixed_weighted_HL_PSL_weighted_compare);

	/*graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID*/
	unordered_map<int, int> vertexID_old_to_new;
	vector<int> vertexID_new_to_old(N);
	for (int i = 0; i < N; i++) {
		vertexID_old_to_new[sorted_vertices[i].first] = i;
		vertexID_new_to_old[i] = sorted_vertices[i].first;
	}
	ideal_graph_579 = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(input_graph, vertexID_old_to_new);

	//graph_v_of_v_idealID_print(ideal_graph_579);


	L_579.resize(N); pos_579.resize(N); increment_579.resize(N);
	for (int v_k = 0; v_k < N; v_k++)  //initialization
	{
		pos_579[v_k] = 0; // record the position of the begining of the label in last iteration
		increment_579[v_k] = 1;
		graph_hash_of_mixed_weighted_HL_PSL_weighted_label xx;
		xx.vertex = v_k;
		xx.distance = 0;
		xx.parent_vertex = v_k;
		L_579[v_k].push_back(xx);
	}


	T_579.resize(num_of_threads);
	dirty_tag_579.resize(num_of_threads);
	for (int i = 0; i < num_of_threads; i++)
	{
		T_579[i].resize(N);
		dirty_tag_579[i].resize(N);
		for (int j = 0; j < N; j++) dirty_tag_579[i][j] = true;
		Qid_579.push(i);
	}

	L_temp_579.resize(N);
	if_continue_579 = true;

	ThreadPool pool(num_of_threads);
	std::vector< std::future<int> > results;

	while (if_continue_579)
	{
		if_continue_579 = false;	
		for (int u = 0; u < N; u++)
		{
			results.emplace_back(
				pool.enqueue([u] { // pass const type value j to thread; [] can be empty
				graph_hash_of_mixed_weighted_HL_PSL_weighted_thread_function_dij(u);
				return 1; // return to results; the return type must be the same with results
			})
			);
		}

		for (auto&& result : results)
			result.get();
		results.clear();

		for (int u = 0; u < N; u++)
		{
			results.emplace_back(
				pool.enqueue([u] { // pass const type value j to thread; [] can be empty
				graph_hash_of_mixed_weighted_HL_PSL_weighted_thread_function_add_new_labels(u);
				return 1; // return to results; the return type must be the same with results
			})
			);
		}

		for (auto&& result : results)
			result.get();
		results.clear();
	}

	/* clean the labels */
	vector<int> pre(N);
	for (int u = 0; u < N; u++)
	{
		/* save the smallest label for each hub of u in T */
		int u_label_size = L_579[u].size();
		for (int i = 0; i < u_label_size; i++)
		{
			int w = L_579[u][i].vertex; float dis = L_579[u][i].distance;
			if (dirty_tag_579[0][w]) // the same hub may have redundancy, record the shortest distance
			{
				dirty_tag_579[0][w] = false;
				T_579[0][w] = dis;
				pre[w] = L_579[u][i].parent_vertex;
			}
			else if (dis < T_579[0][w])
			{
				T_579[0][w] = dis;
				pre[w] = L_579[u][i].parent_vertex;
			}
		}
		for (int i = 0; i < u_label_size; i++)
		{
			int w = L_579[u][i].vertex;
			if (!dirty_tag_579[0][w])
			{
				dirty_tag_579[0][w] = true;
				graph_hash_of_mixed_weighted_HL_PSL_weighted_label xx;
				xx.vertex = w;
				xx.distance = T_579[0][w];
				xx.parent_vertex = pre[w];
				L_temp_579[u].push_back(xx);
			}
		}
		L_579[u].clear(); // do not have two L in RAM simultaneously
	}

	/*return unordered_map_L for old_IDs*/
	vector<vector<graph_hash_of_mixed_weighted_HL_PSL_weighted_label>> output_L(max_N);
	for (int v_k = 0; v_k < N; v_k++) {
		vector<graph_hash_of_mixed_weighted_HL_PSL_weighted_label> xx_vector = {};
		int L_v_k_size = L_temp_579[v_k].size();
		//  std :: cerr << "labels of " << vertexID_new_to_old[v_k] << std :: endl;
		for (int i = 0; i < L_v_k_size; i++)
		{
			auto it = &L_temp_579[v_k][i];
			graph_hash_of_mixed_weighted_HL_PSL_weighted_label_binary_insert(xx_vector, vertexID_new_to_old[it->vertex], vertexID_new_to_old[it->parent_vertex], it->distance);
		} //binary insert here due to mapping, there's no need to sort labels beforehand.
		L_temp_579[v_k].clear();
		output_L[vertexID_new_to_old[v_k]] = xx_vector;
	}

	/*clear global values*/
	ideal_graph_579.clear();
	L_579.clear();
	L_temp_579.clear();
	pos_579.clear();
	increment_579.clear();
	T_579.clear();
	dirty_tag_579.clear();
	Qid_579 = queue<int>(); // IDs of available elements for T


	return output_L;


}


#include <text mining/binary_save_read_vector_of_vectors.h>
void graph_hash_of_mixed_weighted_HL_PSL_weighted_multiple_threads(graph_hash_of_mixed_weighted& input_graph, string index_file_name, int max_N, int num_of_threads) {

	auto begin1 = std::chrono::high_resolution_clock::now();

	vector<vector<graph_hash_of_mixed_weighted_HL_PSL_weighted_label>> L = graph_hash_of_mixed_weighted_HL_PSL_weighted_generate_indexes_multiple_threads(input_graph, max_N, num_of_threads);

	auto end1 = std::chrono::high_resolution_clock::now();
	float runningtime1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1).count() / 1e9; // s

	cout << index_file_name + " PSL_weighted runningtime1: " << runningtime1 << "s" << std::endl;


	auto begin2 = std::chrono::high_resolution_clock::now();

	long long int index_size = 0;
	for (auto it = L.begin(); it != L.end(); it++) {
		index_size = index_size + (*it).size();
	}

	/*save indexes*/
	std::ofstream outputFile;
	outputFile.precision(6);
	outputFile.setf(std::ios::fixed);
	outputFile.setf(std::ios::showpoint);
	outputFile.open("readme_" + index_file_name);
	outputFile << "SECTION Comments" << std::endl;
	outputFile << "Creator Yahui SUN" << std::endl;
	outputFile << "Graph Size: |V|=" << input_graph.hash_of_vectors.size() << " |E|=" << graph_hash_of_mixed_weighted_num_edges(input_graph) << std::endl;
	outputFile << "Total Number of Indexes: " << index_size << " Consumed Time: " << runningtime1 << "s" << std::endl;
	outputFile << "Number of Indexes Per Vertex: " << (float)index_size / input_graph.hash_of_vectors.size() << std::endl;
	outputFile << "Comments END at Line 6" << std::endl;

	binary_save_vector_of_vectors(index_file_name, L);

	auto end2 = std::chrono::high_resolution_clock::now();
	float runningtime2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2).count() / 1e9; // s

	cout << index_file_name + " Saving runningtime2: " << runningtime2 << "s" << std::endl;

}












/*the following codes are for testing

---------------------------------------------------
a cpp file (try.cpp) for running the following test code:
----------------------------------------
#include <iostream>
#include <fstream>
using namespace std;

// header files in the Boost library: https://www.boost.org/
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };

#include <graph_hash_of_mixed_weighted/HL/PSL_weighted/graph_hash_of_mixed_weighted_HL_PSL_weighted_multiple_threads.h>


int main()
{
	test_graph_hash_of_mixed_weighted_HL_PSL_weighted_multiple_threads();
}
------------------------------------------------------------------------------------------
Commends for running the above cpp file on Linux:

g++ -std=c++17 -I/home/boost_1_75_0 -I/root/ysgraph try.cpp -lpthread -Ofast -o A
./A
rm A

(optional to put the above commends in run.sh, and then use the comment: sh run.sh)

*/
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_generate_random_graph.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_read_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_save_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_shortest_paths.h>

void graph_hash_of_mixed_weighted_HL_PSL_weighted_check_correctness_of_PSL_weighted_labels(vector<vector<graph_hash_of_mixed_weighted_HL_PSL_weighted_label>>& L, graph_hash_of_mixed_weighted& instance_graph, int iteration_source_times, int iteration_terminal_times, boost::random::mt19937& boost_random_time_seed) {

	/*below is for checking whether the above labels are right (by randomly computing shortest paths)

	this function can only be used when 0 to n-1 is in the graph, i.e., the graph is an ideal graph

	*/

	boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(instance_graph.hash_of_vectors.size() - 1) };

	for (int yy = 0; yy < iteration_source_times; yy++) {
		int source = dist(boost_random_time_seed);
		std::unordered_map<int, double> distances;
		std::unordered_map<int, int> predecessors;
		graph_hash_of_mixed_weighted_shortest_paths_source_to_all(instance_graph, source, distances, predecessors);

		for (int xx = 0; xx < iteration_terminal_times; xx++) {
			int terminal = dist(boost_random_time_seed);
			float dis = graph_hash_of_mixed_weighted_HL_PSL_weighted_extract_distance(L, source, terminal);
			if (abs(dis - distances[terminal]) > 1e-5 && (dis < std::numeric_limits<float>::max() || distances[terminal] < std::numeric_limits<float>::max())) {
				cout << "source = " << source << endl;
				cout << "terminal = " << terminal << endl;
				cout << "source vector:" << endl;
				for (auto it = L_579[source].begin(); it != L_579[source].end(); it++) {
					cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << ">";
				}
				cout << endl;
				cout << "terminal vector:" << endl;
				for (auto it = L_579[terminal].begin(); it != L_579[terminal].end(); it++) {
					cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << ">";
				}
				cout << endl;

				cout << "dis = " << dis << endl;
				cout << "distances[terminal] = " << distances[terminal] << endl;
				cout << "abs(dis - distances[terminal]) > 1e-5!" << endl;
				getchar();
			}
			vector<pair<int, int>> path = graph_hash_of_mixed_weighted_HL_PSL_weighted_extract_shortest_path(L, source, terminal);
			float path_dis = 0;
			if (path.size() == 0) {
				if (source != terminal) { // disconnected
					path_dis = std::numeric_limits<float>::max();
				}
			}
			else {
				for (auto it = path.begin(); it != path.end(); it++) {
					path_dis = path_dis + graph_hash_of_mixed_weighted_edge_weight(instance_graph, it->first, it->second);
					if (path_dis > std::numeric_limits<float>::max()) {
						path_dis = std::numeric_limits<float>::max();
					}
				}
			}
			if (abs(dis - path_dis) > 1e-5) {
				cout << "source = " << source << endl;
				cout << "terminal = " << terminal << endl;

				cout << "source vector:" << endl;
				for (auto it = L_579[source].begin(); it != L_579[source].end(); it++) {
					cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << ">";
				}
				cout << endl;
				cout << "terminal vector:" << endl;
				for (auto it = L_579[terminal].begin(); it != L_579[terminal].end(); it++) {
					cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << ">";
				}
				cout << endl;

				print_vector_pair_int(path);
				cout << "dis = " << dis << endl;
				cout << "path_dis = " << path_dis << endl;
				cout << "abs(dis - path_dis) > 1e-5!" << endl;
				getchar();
			}

		}

	}

}

void test_graph_hash_of_mixed_weighted_HL_PSL_weighted_multiple_threads() {

	/*parameters*/
	int iteration_graph_times = 100, iteration_source_times = 10, iteration_terminal_times = 10;
	int V = 1e2, E = 5e2, precision = 1, thread_num = 5;
	float ec_min = 0, ec_max = 10; // set ec_min=ec_max=1 for testing unweighted PSL_weighted
	float avg_index_size_per_v = 0;

	bool weighted = true;
	if (ec_min == 1 && ec_max == 1) {
		weighted = false;
	}

	/*iteration*/
	for (int i = 0; i < iteration_graph_times; i++) {

		/*input and output; below is for generating random new graph, or read saved graph*/
		int generate_new_graph = 1;
		std::unordered_set<int> generated_group_vertices;
		graph_hash_of_mixed_weighted instance_graph, generated_group_graph;
		if (generate_new_graph == 1) {
			instance_graph = graph_hash_of_mixed_weighted_generate_random_graph(V, E, 0, 0, ec_min, ec_max, precision, boost_random_time_seed);
			graph_hash_of_mixed_weighted_save_graph_with_weight("simple_iterative_tests.txt", instance_graph, 0);
		}
		else {
			double lambda;
			graph_hash_of_mixed_weighted_read_graph_with_weight("simple_iterative_tests.txt", instance_graph, lambda);
		}
		//graph_hash_of_mixed_weighted_print(instance_graph);

		vector<vector<graph_hash_of_mixed_weighted_HL_PSL_weighted_label>> L = graph_hash_of_mixed_weighted_HL_PSL_weighted_generate_indexes_multiple_threads(instance_graph, V + 1, thread_num);

		graph_hash_of_mixed_weighted_HL_PSL_weighted_check_correctness_of_PSL_weighted_labels(L, instance_graph, iteration_source_times, iteration_terminal_times, boost_random_time_seed);

		long long int index_size = 0;
		for (auto it = L.begin(); it != L.end(); it++) {
			index_size = index_size + (*it).size();
		}
		avg_index_size_per_v = avg_index_size_per_v + (float)index_size / V / iteration_graph_times;

		//cout << "here" << endl;

	}

	cout << "avg_index_size_per_v: " << avg_index_size_per_v << endl;
}
