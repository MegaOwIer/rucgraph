#pragma once
#include <iostream>
#include <ThreadPool.h>
#include <mutex>
#include <stdexcept>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_shortest_paths.h>
#include <boost/heap/fibonacci_heap.hpp>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_generate_random_connected_graph.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_read_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_save_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_PLL_labels.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_PLL_single_thread.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID.h>

static graph_v_of_v_idealID ideal_graph;
static vector<vector<PLL_sorted_label>> L;

static bool this_parallel_PLL_is_running = false;
static int max_N_for_mtx = 3e7;
vector<std::mutex> mtx(max_N_for_mtx);

static queue<int> Qid; // IDs of available elements of P T
static vector<vector<double>> P_dij;
static vector<vector<double>> T_dij;
static vector<vector<int>> P_bfs;
static vector<vector<int>> T_bfs;

void thread_function_dij_mixed(int v_k, int N)
{
	/*Pruned Dijkstra from vertex v_k; see Algorithm 1 in 2013 Japan SIGMOD paper*/

	mtx[max_N_for_mtx-1].lock();
	int used_id = Qid.front();
	Qid.pop();
	mtx[max_N_for_mtx-1].unlock();

	queue<int> P_changed_vertices, T_changed_vertices;
	vector<PLL_handle_t_for_sp> Q_handles(N);

	PLL_node_for_sp node;
	boost::heap::fibonacci_heap<PLL_node_for_sp> Q;
	PLL_sorted_label xx;

	node.vertex = v_k;
	node.parent_vertex = v_k;
	node.priority_value = 0;
	Q_handles[v_k] = Q.push(node);
	P_dij[used_id][v_k] = 0;
	P_changed_vertices.push(v_k);

	mtx[v_k].lock();
	int L_v_k_size = L[v_k].size();
	for (int i = 0; i < L_v_k_size; i++) {
		int L_v_k_i_vertex = L[v_k][i].vertex;
		T_dij[used_id][L_v_k_i_vertex] = L[v_k][i].distance; //allocate T values for L[v_k]
		T_changed_vertices.push(L_v_k_i_vertex);
	}
	mtx[v_k].unlock();
	//因为v-k的标签在从自己出发的过程中不会发生改变，并且在求query的过程中每次都会用到，所以可以提前取出来放在T数组，节省后面查找的时间

	while (Q.size() > 0) {

		node = Q.top();
		Q.pop();
		int u = node.vertex;

		if (v_k <= u) { // this condition is not in 2013 paper, but in 2019 paper (Lemma 3.16)
			int u_parent = node.parent_vertex;
			double P_u = node.priority_value;
			double P_u_with_error = P_u + 1e-5;
			double query_v_k_u = std::numeric_limits<double>::max();

#ifdef _WIN32
			mtx[u].lock();
			auto L_u_size = L[u].size(); // a vector<PLL_sorted_label>
			mtx[u].unlock();
			for (int i = 0; i < L_u_size; i++) {
				mtx[u].lock();      // put lock in for loop is very slow, but it may be the only way under Windows
				double dis = L[u][i].distance + T_dij[used_id][L[u][i].vertex];
				mtx[u].unlock();
				if (query_v_k_u > dis) { query_v_k_u = dis; }
			} //求query的值		
#else
			mtx[u].lock();
			auto L_u_size = L[u].size(); // a vector<PLL_sorted_label>
			for (int i = 0; i < L_u_size; i++) {
				double dis = L[u][i].distance + T_dij[used_id][L[u][i].vertex];   // dont know why this code does not work under Windows
				if (query_v_k_u > dis) { query_v_k_u = dis; }
			} //求query的值
			mtx[u].unlock();
#endif

			if (P_u_with_error < query_v_k_u) { // this is pruning

				xx.vertex = v_k;
				xx.distance = P_u;
				xx.parent_vertex = u_parent;

				mtx[u].lock();
				L[u].push_back(xx); //新增标签
				mtx[u].unlock();
				//下面是dij更新邻接点的过程，同时更新优先队列和距离
				int u_adj_size = ideal_graph[u].size();
				for (int i = 0; i < u_adj_size; i++) {
					int adj_v = ideal_graph[u][i].first; // this needs to be locked
					double ec = ideal_graph[u][i].second;
					if (P_dij[used_id][adj_v] == std::numeric_limits<double>::max()) { //尚未到达的点
						node.vertex = adj_v;
						node.parent_vertex = u;
						node.priority_value = P_u + ec;
						Q_handles[adj_v] = Q.push(node);
						P_dij[used_id][adj_v] = node.priority_value;
						P_changed_vertices.push(adj_v);
					}
					else {
						if (P_dij[used_id][adj_v] > P_u + ec) {
							node.vertex = adj_v;
							node.parent_vertex = u;
							node.priority_value = P_u + ec;
							Q.update(Q_handles[adj_v], node);
							P_dij[used_id][adj_v] = node.priority_value;
						}
					}
				}

				
			}


		}
	}

	while (P_changed_vertices.size() > 0) {
		P_dij[used_id][P_changed_vertices.front()] = std::numeric_limits<double>::max(); // reverse-allocate P values
		P_changed_vertices.pop();
	}
	while (T_changed_vertices.size() > 0) {
		T_dij[used_id][T_changed_vertices.front()] = std::numeric_limits<double>::max(); // reverse-allocate T values
		T_changed_vertices.pop();
	}


	mtx[max_N_for_mtx - 1].lock();
	Qid.push(used_id);
	mtx[max_N_for_mtx - 1].unlock();

}

void thread_function_bfs_mixed(int v_k, int N)
{
	mtx[max_N_for_mtx-1].lock();
	int used_id = Qid.front(); // which P T elements does this thread use
	Qid.pop();
	mtx[max_N_for_mtx-1].unlock();

	queue<int> P_changed_vertices, T_changed_vertices;

	queue<PLL_node_for_sp> Q;

	PLL_node_for_sp node;
	PLL_sorted_label xx;

	node.vertex = v_k;
	node.parent_vertex = v_k;
	node.priority_value = 0;
	Q.push(node);
	P_bfs[used_id][v_k] = 0;
	P_changed_vertices.push(v_k);

	mtx[v_k].lock();
	int L_v_k_size = L[v_k].size();
	for (int i = 0; i < L_v_k_size; i++) {
		int L_v_k_i_vertex = L[v_k][i].vertex;
		T_bfs[used_id][L_v_k_i_vertex] = L[v_k][i].distance; //allocate T values for L[v_k]
		T_changed_vertices.push(L_v_k_i_vertex);
	}
	mtx[v_k].unlock();

	while (Q.size() > 0) {

		node = Q.front();
		Q.pop();
		int u = node.vertex;

		if (v_k <= u) { // this condition is not in 2013 paper, but in 2019 paper (Lemma 3.16)

			int u_parent = node.parent_vertex;
			int P_u = node.priority_value;

			double query_v_k_u = std::numeric_limits<double>::max();

			mtx[u].lock();
			int L_u_size = L[u].size();
			for (int i = 0; i < L_u_size; i++) {
				double dis = L[u][i].distance + T_bfs[used_id][L[u][i].vertex];		 // cannot lock mtx[u] in this for loop, since the locking time is very large		
				if (query_v_k_u > dis) { query_v_k_u = dis; }
			}
			mtx[u].unlock();

			if (P_u < query_v_k_u) { // this is pruning

				xx.vertex = v_k;
				xx.distance = P_u;
				xx.parent_vertex = u_parent;

				mtx[u].lock();
				L[u].push_back(xx); //新增标签
				mtx[u].unlock();
				auto u_adj_size = ideal_graph[u].size();
				for (int i = 0; i < u_adj_size; i++) {
					int adj_v = ideal_graph[u][i].first;
					if (P_bfs[used_id][adj_v] == INT_MAX) {
						node.vertex = adj_v;
						node.parent_vertex = u;
						node.priority_value = P_u + 1;
						Q.push(node);
						P_bfs[used_id][adj_v] = node.priority_value;
						P_changed_vertices.push(adj_v);
					}
				}

				


			}
		}
	}

	while (P_changed_vertices.size() > 0) {
		P_bfs[used_id][P_changed_vertices.front()] = INT_MAX; // reverse-allocate P values
		P_changed_vertices.pop();
	}
	while (T_changed_vertices.size() > 0) {
		T_bfs[used_id][T_changed_vertices.front()] = INT_MAX; // reverse-allocate T values
		T_changed_vertices.pop();
	}

	mtx[max_N_for_mtx - 1].lock();
	Qid.push(used_id);
	mtx[max_N_for_mtx - 1].unlock();
}


/*the following parallel PLL code cannot be run parallelly, due to the above (static) globel values*/

vector<vector<PLL_sorted_label>> PLL_generate_indexes_multiple_threads(graph_hash_of_mixed_weighted& input_graph, int max_N, bool weighted, int num_of_threads)
{

	if (max_N > max_N_for_mtx) {
		cout << "max_N > max_N_for_mtx; max_N_for_mtx is too small! location: PLL_generate_indexes_unweighted_multiple_threads_onlyforLinux" << endl;
		exit(1);
	}

	mtx[max_N_for_mtx - 1].lock();
	if (this_parallel_PLL_is_running == true) {
		cout << "the following parallel PLL code cannot be run parallelly, due to the above (static) globel values" << endl;
		exit(1);
	}
	this_parallel_PLL_is_running = true;
	mtx[max_N_for_mtx - 1].unlock();

	L.resize(max_N);

	int N = input_graph.hash_of_vectors.size();

	vector<pair<int, int>> sorted_vertices;
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		sorted_vertices.push_back({ it->first, input_graph.degree(it->first) });
	}
	sort(sorted_vertices.begin(), sorted_vertices.end(), PLL_Creating_Indexes_compare);


	unordered_map<int, int> vertexID_old_to_new;
	vector<int> vertexID_new_to_old(N);
	for (int i = 0; i < N; i++) {
		vertexID_old_to_new[sorted_vertices[i].first] = i;
		vertexID_new_to_old[i] = sorted_vertices[i].first;
	}
	ideal_graph = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(input_graph, vertexID_old_to_new);


	ThreadPool pool(num_of_threads);
	std::vector< std::future<int> > results; // return typename: xxx
	if (weighted) {
		P_dij.resize(num_of_threads);
		T_dij.resize(num_of_threads);
		for (int i = 0; i < num_of_threads; i++)
		{
			P_dij[i].resize(N);
			T_dij[i].resize(N);
			for (int j = 0; j < N; j++)
			{
				P_dij[i][j] = std::numeric_limits<double>::max();
				T_dij[i][j] = std::numeric_limits<double>::max();
			}
			Qid.push(i);
		}
		for (int v_k = 0; v_k < N; v_k++) {
			results.emplace_back(
				pool.enqueue([v_k, N] { // pass const type value j to thread; [] can be empty
				thread_function_dij_mixed(v_k, N);
				return 1; // return to results; the return type must be the same with results
			})
			);
		}
	}
	else {
		P_bfs.resize(num_of_threads);
		T_bfs.resize(num_of_threads);
		for (int i = 0; i < num_of_threads; i++)
		{
			P_bfs[i].resize(N);
			T_bfs[i].resize(N);
			for (int j = 0; j < N; j++)
			{
				P_bfs[i][j] = INT_MAX;
				T_bfs[i][j] = INT_MAX;
			}
			Qid.push(i);

		}
		for (int v_k = 0; v_k < N; v_k++) {
			results.emplace_back(
				pool.enqueue([v_k, N] { // pass const type value j to thread; [] can be empty
				thread_function_bfs_mixed(v_k, N);
				return 1; // return to results; the return type must be the same with results
			})
			);
		}
	}
	for (auto&& result : results)
		result.get(); // all threads finish here

	/*return L for old_IDs*/
	vector<vector<PLL_sorted_label>> output_L(max_N);
	for (int v_k = 0; v_k < N; v_k++) {
		vector<PLL_sorted_label> xx_vector = {};
		int L_v_k_size = L[v_k].size();
		for (int i = 0; i < L_v_k_size; i++) {
			PLL_sorted_label_binary_insert(xx_vector, vertexID_new_to_old[L[v_k][i].vertex], vertexID_new_to_old[L[v_k][i].parent_vertex], L[v_k][i].distance);
		}
		output_L[vertexID_new_to_old[v_k]] = xx_vector;
	}

	/*clear static global*/
	graph_v_of_v_idealID().swap(ideal_graph);
	vector<vector<PLL_sorted_label>>().swap(L);
	vector<vector<double>>().swap(P_dij);
	vector<vector<double>>().swap(T_dij);
	vector<vector<int>>().swap(P_bfs);
	vector<vector<int>>().swap(T_bfs);
	Qid = queue<int>();
	this_parallel_PLL_is_running = false;


	return output_L;

}

void PLL_generate_and_save_indexes_multiple_threads(graph_hash_of_mixed_weighted& input_graph, string index_file_name, int max_N, bool weighted, int num_of_threads) {

	auto begin1 = std::chrono::high_resolution_clock::now();

	vector<vector<PLL_sorted_label>> L = PLL_generate_indexes_multiple_threads(input_graph, max_N, weighted, num_of_threads);

	auto end1 = std::chrono::high_resolution_clock::now();
	double runningtime1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1).count() / 1e9; // s

	cout << index_file_name + " PLL runningtime1: " << runningtime1 << "s" << std::endl;


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
	outputFile << "Total Number of Indexes: " << index_size << " Consumed Time: " << runningtime1 << "s"
		<< " Consumed RAM: " << Current_Memory_Consumption_of_This_Process() << "MB" << std::endl;
	outputFile << "Number of Indexes Per Vertex: " << (double)index_size / input_graph.hash_of_vectors.size() << std::endl;
	outputFile << "Comments END at Line 6" << std::endl;

	binary_save_vector_of_vectors(index_file_name, L);

	auto end2 = std::chrono::high_resolution_clock::now();
	double runningtime2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2).count() / 1e9; // s

	cout << index_file_name + " Saving runningtime2: " << runningtime2 << "s" << std::endl;

}












/*the following codes are for testing

---------------------------------------------------
a cpp file (try.cpp) for running the following test code:
----------------------------------------
#include <iostream>
#include <string>
#include <chrono>
#include <list>
#include <fstream>
using namespace std;

// header files in the Boost library: https://www.boost.org/
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };

#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_PLL_multiple_threads.h>

int main()
{
	PLL_generate_indexes_single_thread_VS_multiple_threads();
}
------------------------------------------------------------------------------------------
Commends for running the above cpp file on Linux:

g++ -std=c++17 -I/home/boost_1_75_0 -I/root/ysgraph try.cpp -lpthread -Ofast -o A
./A
rm A

(optional to put the above commends in run.sh, and then use the comment: sh run.sh)

*/

void test_PLL_generate_indexes_multiple_threads() {

	/*parameters*/
	int iteration_graph_times = 100, iteration_source_times = 10, iteration_terminal_times = 10;
	int V = 1e2, E = 2e2, precision = 1, thread_num = 5;
	double ec_min = 0, ec_max = 10; // set ec_min=ec_max=1 for testing unweighted PLL
	double avg_index_size_per_v = 0;

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

		///*below is for generating PLL labels*/
		PLL_generate_and_save_indexes_multiple_threads(instance_graph, "simple_iterative_tests_PLL_Indexes.txt", V + 1, weighted, thread_num);
		///*below is for reading PLL labels*/
		//vector<vector<PLL_sorted_label>> L = PLL_read_indexes_vectorFORMAT_binary("simple_iterative_tests_PLL_Indexes.txt");

		
		vector<vector<PLL_sorted_label>> L = PLL_generate_indexes_multiple_threads(instance_graph, V + 1, weighted, thread_num);

		check_correctness_of_PLL_labels(L, instance_graph, iteration_source_times, iteration_terminal_times, boost_random_time_seed);

		long long int index_size = 0;
		for (auto it = L.begin(); it != L.end(); it++) {
			index_size = index_size + (*it).size();
		}
		avg_index_size_per_v = avg_index_size_per_v + (double)index_size / V / iteration_graph_times;

	}

	cout << "avg_index_size_per_v: " << avg_index_size_per_v << endl;
}


void PLL_generate_indexes_single_thread_VS_multiple_threads() {

	/*parameters*/
	int iteration_graph_times = 2;
	int V = 1e3, E = 1e4, precision = 1, thread_num = 5;
	double ec_min = 1, ec_max = 10; // set ec_min=ec_max=1 for testing unweighted PLL

	double avg_index_size_per_v_unweighted_single_thread = 0, avg_index_size_per_v_weighted_single_thread = 0, avg_index_size_per_v_unweighted_multiple_threads = 0, avg_index_size_per_v_weighted_multiple_threads = 0;
	double PLL_time_per_graph_unweighted_single_thread = 0, PLL_time_per_graph_weighted_single_thread = 0, PLL_time_per_graph_unweighted_multiple_threads = 0, PLL_time_per_graph_weighted_multiple_threads = 0;

	/*iteration*/
	for (int i = 0; i < iteration_graph_times; i++) {

		/*unweighted*/
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
				vector<vector<PLL_sorted_label>> L = PLL_generate_indexes_unweighted_single_thread(instance_graph, V + 1);
				auto end = std::chrono::high_resolution_clock::now();
				double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
				PLL_time_per_graph_unweighted_single_thread = PLL_time_per_graph_unweighted_single_thread + runningtime / iteration_graph_times;

				long long int index_size = 0;
				for (auto it = L.begin(); it != L.end(); it++) {
					index_size = index_size + (*it).size();
				}
				avg_index_size_per_v_unweighted_single_thread = avg_index_size_per_v_unweighted_single_thread + (double)index_size / V / iteration_graph_times;
			}

			if (1) {
				auto begin = std::chrono::high_resolution_clock::now();
				vector<vector<PLL_sorted_label>> L = PLL_generate_indexes_multiple_threads(instance_graph, V + 1, false, thread_num);
				auto end = std::chrono::high_resolution_clock::now();
				double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
				PLL_time_per_graph_unweighted_multiple_threads = PLL_time_per_graph_unweighted_multiple_threads + runningtime / iteration_graph_times;

				long long int index_size = 0;
				for (auto it = L.begin(); it != L.end(); it++) {
					index_size = index_size + (*it).size();
				}
				avg_index_size_per_v_unweighted_multiple_threads = avg_index_size_per_v_unweighted_multiple_threads + (double)index_size / V / iteration_graph_times;
			}
		}


		/*weighted*/
		if (1) {
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

			if (1) {
				auto begin = std::chrono::high_resolution_clock::now();
				vector<vector<PLL_sorted_label>> L = PLL_generate_indexes_weighted_single_thread(instance_graph, V + 1);
				auto end = std::chrono::high_resolution_clock::now();
				double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
				PLL_time_per_graph_weighted_single_thread = PLL_time_per_graph_weighted_single_thread + runningtime / iteration_graph_times;

				long long int index_size = 0;
				for (auto it = L.begin(); it != L.end(); it++) {
					index_size = index_size + (*it).size();
				}
				avg_index_size_per_v_weighted_single_thread = avg_index_size_per_v_weighted_single_thread + (double)index_size / V / iteration_graph_times;
			}

			if (1) {
				auto begin = std::chrono::high_resolution_clock::now();
				vector<vector<PLL_sorted_label>> L = PLL_generate_indexes_multiple_threads(instance_graph, V + 1, true, thread_num);
				auto end = std::chrono::high_resolution_clock::now();
				double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
				PLL_time_per_graph_weighted_multiple_threads = PLL_time_per_graph_weighted_multiple_threads + runningtime / iteration_graph_times;

				long long int index_size = 0;
				for (auto it = L.begin(); it != L.end(); it++) {
					index_size = index_size + (*it).size();
				}
				avg_index_size_per_v_weighted_multiple_threads = avg_index_size_per_v_weighted_multiple_threads + (double)index_size / V / iteration_graph_times;
			}
		}


	}

	cout << "thread_num: " << thread_num << endl;
	cout << "avg_index_size_per_v_unweighted_single_thread: " << avg_index_size_per_v_unweighted_single_thread << endl;
	cout << "avg_index_size_per_v_unweighted_multiple_threads: " << avg_index_size_per_v_unweighted_multiple_threads << endl;
	cout << "avg_index_size_per_v_unweighted_multiple_threads / avg_index_size_per_v_unweighted_single_thread: " << avg_index_size_per_v_unweighted_multiple_threads / avg_index_size_per_v_unweighted_single_thread << endl;
	cout << endl;
	cout << "avg_index_size_per_v_weighted_single_thread: " << avg_index_size_per_v_weighted_single_thread << endl;
	cout << "avg_index_size_per_v_weighted_multiple_threads: " << avg_index_size_per_v_weighted_multiple_threads << endl;
	cout << "avg_index_size_per_v_weighted_multiple_threads / avg_index_size_per_v_weighted_single_thread: " << avg_index_size_per_v_weighted_multiple_threads / avg_index_size_per_v_weighted_single_thread << endl;
	cout << endl;
	cout << "PLL_time_per_graph_unweighted_single_thread: " << PLL_time_per_graph_unweighted_single_thread << endl;
	cout << "PLL_time_per_graph_unweighted_multiple_threads: " << PLL_time_per_graph_unweighted_multiple_threads << endl;
	cout << "PLL_time_per_graph_unweighted_single_thread / PLL_time_per_graph_unweighted_multiple_threads: " << PLL_time_per_graph_unweighted_single_thread / PLL_time_per_graph_unweighted_multiple_threads << endl;
	cout << endl;
	cout << "PLL_time_per_graph_weighted_single_thread: " << PLL_time_per_graph_weighted_single_thread << endl;
	cout << "PLL_time_per_graph_weighted_multiple_threads: " << PLL_time_per_graph_weighted_multiple_threads << endl;
	cout << "PLL_time_per_graph_weighted_single_thread / PLL_time_per_graph_weighted_multiple_threads: " << PLL_time_per_graph_weighted_single_thread / PLL_time_per_graph_weighted_multiple_threads << endl;

	/*
	
	for V = 1e3, E = 1e4
	thread_num: 50
avg_index_size_per_v_unweighted_single_thread: 1008.63
avg_index_size_per_v_unweighted_multiple_threads: 1012.19
avg_index_size_per_v_unweighted_multiple_threads / avg_index_size_per_v_unweighted_single_thread: 1.00354

avg_index_size_per_v_weighted_single_thread: 1099.09
avg_index_size_per_v_weighted_multiple_threads: 1101.4
avg_index_size_per_v_weighted_multiple_threads / avg_index_size_per_v_weighted_single_thread: 1.0021

PLL_time_per_graph_unweighted_single_thread: 54.1492
PLL_time_per_graph_unweighted_multiple_threads: 3.19735
PLL_time_per_graph_unweighted_single_thread / PLL_time_per_graph_unweighted_multiple_threads: 16.9357

PLL_time_per_graph_weighted_single_thread: 70.7582
PLL_time_per_graph_weighted_multiple_threads: 4.13053
PLL_time_per_graph_weighted_single_thread / PLL_time_per_graph_weighted_multiple_threads: 17.1305
	
	*/

}
