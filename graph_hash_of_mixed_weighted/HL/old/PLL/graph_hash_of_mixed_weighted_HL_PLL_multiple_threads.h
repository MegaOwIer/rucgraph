#pragma once
#include <iostream>
#include <ThreadPool.h>
#include <mutex>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>
#include <boost/heap/fibonacci_heap.hpp>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID.h>
#include <graph_hash_of_mixed_weighted/HL/PLL/graph_hash_of_mixed_weighted_HL_PLL_labels.h>

/*unique code for this file: 987*/


/*static global values; these values cannot be passed to multiple threads as parameters*/
static graph_v_of_v_idealID ideal_graph_987;
static vector<vector<graph_hash_of_mixed_weighted_HL_PLL_label>> L_987;
static bool this_parallel_PLL_is_running_987 = false;
static int max_N_for_mtx_987 = 1e7;  // this is the max N to run
static vector<std::mutex> mtx_987(max_N_for_mtx_987);  // std::mutex has no copy or move constructor, while std::vector::resize() requires that; you cannot resize mtx;    moreover, do not change mtx to a pointer and then points to local values, it is very slow!!
static queue<int> Qid_987; // IDs of available elements of P T
static vector<vector<float>> P_dij_987;
static vector<vector<float>> T_dij_987;
static vector<vector<int>> P_bfs_987;
static vector<vector<int>> T_bfs_987;

bool graph_hash_of_mixed_weighted_HL_PLL_compare(const pair<int, int>& i, pair<int, int>& j)
{
	/*< is nearly 10 times slower than >*/
	return i.second > j.second;  // < is from small to big; > is from big to small.  sort by the second item of pair<int, int>
}

struct graph_hash_of_mixed_weighted_HL_PLL_node_for_sp {
public:
	int vertex, parent_vertex;
	float priority_value;
}; // define the node in the queue
bool operator<(graph_hash_of_mixed_weighted_HL_PLL_node_for_sp const& x, graph_hash_of_mixed_weighted_HL_PLL_node_for_sp const& y) {
	return x.priority_value > y.priority_value; // < is the max-heap; > is the min heap
}
typedef typename boost::heap::fibonacci_heap<graph_hash_of_mixed_weighted_HL_PLL_node_for_sp>::handle_type graph_hash_of_mixed_weighted_HL_PLL_handle_t_for_sp;


void graph_hash_of_mixed_weighted_HL_PLL_thread_function_dij_mixed(int v_k, int N)
{
	/*Pruned Dijkstra from vertex v_k; see Algorithm 1 in 2013 Japan SIGMOD paper*/

	mtx_987[max_N_for_mtx_987 -1].lock();
	int used_id = Qid_987.front();
	Qid_987.pop();
	mtx_987[max_N_for_mtx_987 -1].unlock();

	queue<int> P_changed_vertices, T_changed_vertices;
	vector<graph_hash_of_mixed_weighted_HL_PLL_handle_t_for_sp> Q_handles(N);

	graph_hash_of_mixed_weighted_HL_PLL_node_for_sp node;
	boost::heap::fibonacci_heap<graph_hash_of_mixed_weighted_HL_PLL_node_for_sp> Q;
	graph_hash_of_mixed_weighted_HL_PLL_label xx;

	node.vertex = v_k;
	node.parent_vertex = v_k;
	node.priority_value = 0;
	Q_handles[v_k] = Q.push(node);
	P_dij_987[used_id][v_k] = 0;
	P_changed_vertices.push(v_k);

	mtx_987[v_k].lock();
	int L_v_k_size = L_987[v_k].size();
	for (int i = 0; i < L_v_k_size; i++) {
		int L_v_k_i_vertex = L_987[v_k][i].vertex;
		T_dij_987[used_id][L_v_k_i_vertex] = L_987[v_k][i].distance; //allocate T values for L_987[v_k]
		T_changed_vertices.push(L_v_k_i_vertex);
	}
	mtx_987[v_k].unlock();
	//因为v-k的标签在从自己出发的过程中不会发生改变，并且在求query的过程中每次都会用到，所以可以提前取出来放在T数组，节省后面查找的时间

	while (Q.size() > 0) {

		node = Q.top();
		Q.pop();
		int u = node.vertex;

		if (v_k <= u) { // this condition is not in 2013 paper, but in 2019 paper (Lemma 3.16)
			int u_parent = node.parent_vertex;
			float P_u = node.priority_value;
			float P_u_with_error = P_u + 1e-5;
			float query_v_k_u = std::numeric_limits<float>::max();

#ifdef _WIN32
			mtx_987[u].lock();
			auto L_u_size = L_987[u].size(); // a vector<PLL_sorted_label>
			mtx_987[u].unlock();
			for (int i = 0; i < L_u_size; i++) {
				mtx_987[u].lock();      // put lock in for loop is very slow, but it may be the only way under Windows
				float dis = L_987[u][i].distance + T_dij_987[used_id][L_987[u][i].vertex];
				mtx_987[u].unlock();
				if (query_v_k_u > dis) { query_v_k_u = dis; }
			} //求query的值		
#else
			mtx_987[u].lock();
			auto L_u_size = L_987[u].size(); // a vector<PLL_sorted_label>
			for (int i = 0; i < L_u_size; i++) {
				float dis = L_987[u][i].distance + T_dij_987[used_id][L_987[u][i].vertex];   // dont know why this code does not work under Windows
				if (query_v_k_u > dis) { query_v_k_u = dis; }
			} //求query的值
			mtx_987[u].unlock();
#endif

			if (P_u_with_error < query_v_k_u) { // this is pruning

				xx.vertex = v_k;
				xx.distance = P_u;
				xx.parent_vertex = u_parent;

				mtx_987[u].lock();
				L_987[u].push_back(xx); //新增标签
				mtx_987[u].unlock();
				//下面是dij更新邻接点的过程，同时更新优先队列和距离
				int u_adj_size = ideal_graph_987[u].size();
				for (int i = 0; i < u_adj_size; i++) {
					int adj_v = ideal_graph_987[u][i].first; // this needs to be locked
					float ec = ideal_graph_987[u][i].second;
					if (P_dij_987[used_id][adj_v] == std::numeric_limits<float>::max()) { //尚未到达的点
						node.vertex = adj_v;
						node.parent_vertex = u;
						node.priority_value = P_u + ec;
						Q_handles[adj_v] = Q.push(node);
						P_dij_987[used_id][adj_v] = node.priority_value;
						P_changed_vertices.push(adj_v);
					}
					else {
						if (P_dij_987[used_id][adj_v] > P_u + ec) {
							node.vertex = adj_v;
							node.parent_vertex = u;
							node.priority_value = P_u + ec;
							Q.update(Q_handles[adj_v], node);
							P_dij_987[used_id][adj_v] = node.priority_value;
						}
					}
				}

				
			}


		}
	}

	while (P_changed_vertices.size() > 0) {
		P_dij_987[used_id][P_changed_vertices.front()] = std::numeric_limits<float>::max(); // reverse-allocate P values
		P_changed_vertices.pop();
	}
	while (T_changed_vertices.size() > 0) {
		T_dij_987[used_id][T_changed_vertices.front()] = std::numeric_limits<float>::max(); // reverse-allocate T values
		T_changed_vertices.pop();
	}


	mtx_987[max_N_for_mtx_987 - 1].lock();
	Qid_987.push(used_id);
	mtx_987[max_N_for_mtx_987 - 1].unlock();

}

void graph_hash_of_mixed_weighted_HL_PLL_thread_function_bfs_mixed(int v_k, int N)
{
	mtx_987[max_N_for_mtx_987 -1].lock();
	int used_id = Qid_987.front(); // which P T elements does this thread use
	Qid_987.pop();
	mtx_987[max_N_for_mtx_987 -1].unlock();

	queue<int> P_changed_vertices, T_changed_vertices;

	queue<graph_hash_of_mixed_weighted_HL_PLL_node_for_sp> Q;

	graph_hash_of_mixed_weighted_HL_PLL_node_for_sp node;
	graph_hash_of_mixed_weighted_HL_PLL_label xx;

	node.vertex = v_k;
	node.parent_vertex = v_k;
	node.priority_value = 0;
	Q.push(node);
	P_bfs_987[used_id][v_k] = 0;
	P_changed_vertices.push(v_k);

	mtx_987[v_k].lock();
	int L_v_k_size = L_987[v_k].size();
	for (int i = 0; i < L_v_k_size; i++) {
		int L_v_k_i_vertex = L_987[v_k][i].vertex;
		T_bfs_987[used_id][L_v_k_i_vertex] = L_987[v_k][i].distance; //allocate T values for L_987[v_k]
		T_changed_vertices.push(L_v_k_i_vertex);
	}
	mtx_987[v_k].unlock();

	while (Q.size() > 0) {

		node = Q.front();
		Q.pop();
		int u = node.vertex;

		if (v_k <= u) { // this condition is not in 2013 paper, but in 2019 paper (Lemma 3.16)

			int u_parent = node.parent_vertex;
			int P_u = node.priority_value;

			float query_v_k_u = std::numeric_limits<float>::max();

			mtx_987[u].lock();
			int L_u_size = L_987[u].size();
			for (int i = 0; i < L_u_size; i++) {
				float dis = L_987[u][i].distance + T_bfs_987[used_id][L_987[u][i].vertex];		 // cannot lock mtx_987[u] in this for loop, since the locking time is very large		
				if (query_v_k_u > dis) { query_v_k_u = dis; }
			}
			mtx_987[u].unlock();

			if (P_u < query_v_k_u) { // this is pruning

				xx.vertex = v_k;
				xx.distance = P_u;
				xx.parent_vertex = u_parent;

				mtx_987[u].lock();
				L_987[u].push_back(xx); //新增标签
				mtx_987[u].unlock();
				auto u_adj_size = ideal_graph_987[u].size();
				for (int i = 0; i < u_adj_size; i++) {
					int adj_v = ideal_graph_987[u][i].first;
					if (P_bfs_987[used_id][adj_v] == INT_MAX) {
						node.vertex = adj_v;
						node.parent_vertex = u;
						node.priority_value = P_u + 1;
						Q.push(node);
						P_bfs_987[used_id][adj_v] = node.priority_value;
						P_changed_vertices.push(adj_v);
					}
				}

				


			}
		}
	}

	while (P_changed_vertices.size() > 0) {
		P_bfs_987[used_id][P_changed_vertices.front()] = INT_MAX; // reverse-allocate P values
		P_changed_vertices.pop();
	}
	while (T_changed_vertices.size() > 0) {
		T_bfs_987[used_id][T_changed_vertices.front()] = INT_MAX; // reverse-allocate T values
		T_changed_vertices.pop();
	}

	mtx_987[max_N_for_mtx_987 - 1].lock();
	Qid_987.push(used_id);
	mtx_987[max_N_for_mtx_987 - 1].unlock();
}


/*the following parallel PLL code cannot be run parallelly, due to the above (static) globel values*/

vector<vector<graph_hash_of_mixed_weighted_HL_PLL_label>> graph_hash_of_mixed_weighted_HL_PLL_generate_indexes_multiple_threads(graph_hash_of_mixed_weighted& input_graph, int max_N, bool weighted, int num_of_threads)
{

	if (max_N > max_N_for_mtx_987) {
		cout << "max_N > max_N_for_mtx; max_N_for_mtx is too small!" << endl;
		exit(1);
	}

	mtx_987[max_N_for_mtx_987 - 1].lock();
	if (this_parallel_PLL_is_running_987 == true) {
		cout << "the following parallel PLL code cannot be run parallelly, due to the above (static) globel values" << endl;
		exit(1);
	}
	this_parallel_PLL_is_running_987 = true;
	mtx_987[max_N_for_mtx_987 - 1].unlock();

	L_987.resize(max_N);

	int N = input_graph.hash_of_vectors.size();

	vector<pair<int, int>> sorted_vertices;
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		sorted_vertices.push_back({ it->first, input_graph.degree(it->first) });
	}
	sort(sorted_vertices.begin(), sorted_vertices.end(), graph_hash_of_mixed_weighted_HL_PLL_compare);


	unordered_map<int, int> vertexID_old_to_new;
	vector<int> vertexID_new_to_old(N);
	for (int i = 0; i < N; i++) {
		vertexID_old_to_new[sorted_vertices[i].first] = i;
		vertexID_new_to_old[i] = sorted_vertices[i].first;
	}
	ideal_graph_987 = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(input_graph, vertexID_old_to_new);


	ThreadPool pool(num_of_threads);
	std::vector< std::future<int> > results; // return typename: xxx
	if (weighted) {
		P_dij_987.resize(num_of_threads);
		T_dij_987.resize(num_of_threads);
		for (int i = 0; i < num_of_threads; i++)
		{
			P_dij_987[i].resize(N);
			T_dij_987[i].resize(N);
			for (int j = 0; j < N; j++)
			{
				P_dij_987[i][j] = std::numeric_limits<float>::max();
				T_dij_987[i][j] = std::numeric_limits<float>::max();
			}
			Qid_987.push(i);
		}
		for (int v_k = 0; v_k < N; v_k++) {
			results.emplace_back(
				pool.enqueue([v_k, N] { // pass const type value j to thread; [] can be empty
				graph_hash_of_mixed_weighted_HL_PLL_thread_function_dij_mixed(v_k, N);
				return 1; // return to results; the return type must be the same with results
			})
			);
		}
	}
	else {
		P_bfs_987.resize(num_of_threads);
		T_bfs_987.resize(num_of_threads);
		for (int i = 0; i < num_of_threads; i++)
		{
			P_bfs_987[i].resize(N);
			T_bfs_987[i].resize(N);
			for (int j = 0; j < N; j++)
			{
				P_bfs_987[i][j] = INT_MAX;
				T_bfs_987[i][j] = INT_MAX;
			}
			Qid_987.push(i);

		}
		for (int v_k = 0; v_k < N; v_k++) {
			results.emplace_back(
				pool.enqueue([v_k, N] { // pass const type value j to thread; [] can be empty
				graph_hash_of_mixed_weighted_HL_PLL_thread_function_bfs_mixed(v_k, N);
				return 1; // return to results; the return type must be the same with results
			})
			);
		}
	}
	for (auto&& result : results)
		result.get(); // all threads finish here

	/*return L for old_IDs*/
	vector<vector<graph_hash_of_mixed_weighted_HL_PLL_label>> output_L(max_N);
	for (int v_k = 0; v_k < N; v_k++) {
		vector<graph_hash_of_mixed_weighted_HL_PLL_label> xx_vector = {};
		int L_v_k_size = L_987[v_k].size();
		for (int i = 0; i < L_v_k_size; i++) {
			auto it = &L_987[v_k][i];
			graph_hash_of_mixed_weighted_HL_PLL_label_binary_insert(xx_vector, vertexID_new_to_old[it->vertex], vertexID_new_to_old[it->parent_vertex], it->distance);
		}
		L_987[v_k].clear();
		output_L[vertexID_new_to_old[v_k]] = xx_vector;
	}

	/*clear static global*/
	ideal_graph_987.clear();
	L_987.clear();
	P_dij_987.clear();
	T_dij_987.clear();
	P_bfs_987.clear();
	T_bfs_987.clear();
	Qid_987 = queue<int>();
	this_parallel_PLL_is_running_987 = false;


	return output_L;

}


#include <text mining/binary_save_read_vector_of_vectors.h>
void graph_hash_of_mixed_weighted_HL_PLL_multiple_threads(graph_hash_of_mixed_weighted& input_graph, string index_file_name, int max_N, bool weighted, int num_of_threads) {

	auto begin1 = std::chrono::high_resolution_clock::now();

	vector<vector<graph_hash_of_mixed_weighted_HL_PLL_label>> L = graph_hash_of_mixed_weighted_HL_PLL_generate_indexes_multiple_threads(input_graph, max_N, weighted, num_of_threads);

	auto end1 = std::chrono::high_resolution_clock::now();
	float runningtime1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1).count() / 1e9; // s

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

#include <graph_hash_of_mixed_weighted/HL/PLL/graph_hash_of_mixed_weighted_HL_PLL_multiple_threads.h>


int main()
{
	test_graph_hash_of_mixed_weighted_HL_PLL_multiple_threads();
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

void graph_hash_of_mixed_weighted_HL_PLL_check_correctness_of_PLL_labels(vector<vector<graph_hash_of_mixed_weighted_HL_PLL_label>>& L, graph_hash_of_mixed_weighted& instance_graph, int iteration_source_times, int iteration_terminal_times, boost::random::mt19937& boost_random_time_seed) {

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
			float dis = graph_hash_of_mixed_weighted_HL_PLL_extract_distance(L, source, terminal);
			if (abs(dis - distances[terminal]) > 1e-5 && (dis< std::numeric_limits<float>::max() || distances[terminal] < std::numeric_limits<float>::max()) ) {
				cout << "source = " << source << endl;
				cout << "terminal = " << terminal << endl;
				cout << "source vector:" << endl;
				for (auto it = L[source].begin(); it != L[source].end(); it++) {
					cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << ">";
				}
				cout << endl;
				cout << "terminal vector:" << endl;
				for (auto it = L[terminal].begin(); it != L[terminal].end(); it++) {
					cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << ">";
				}
				cout << endl;

				cout << "dis = " << dis << endl;
				cout << "distances[terminal] = " << distances[terminal] << endl;
				cout << "abs(dis - distances[terminal]) > 1e-5!" << endl;
				getchar();
			}
			vector<pair<int, int>> path = graph_hash_of_mixed_weighted_HL_PLL_extract_shortest_path(L, source, terminal);
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
				for (auto it = L[source].begin(); it != L[source].end(); it++) {
					cout << "<" << it->vertex << "," << it->distance << "," << it->parent_vertex << ">";
				}
				cout << endl;
				cout << "terminal vector:" << endl;
				for (auto it = L[terminal].begin(); it != L[terminal].end(); it++) {
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

void test_graph_hash_of_mixed_weighted_HL_PLL_multiple_threads() {

	/*parameters*/
	int iteration_graph_times = 100, iteration_source_times = 10, iteration_terminal_times = 10;
	int V = 1e2, E = 2e2, precision = 1, thread_num = 5;
	float ec_min = 0, ec_max = 10; // set ec_min=ec_max=1 for testing unweighted PLL
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

		///*below is for generating PLL labels*/
		//graph_hash_of_mixed_weighted_HL_PLL_multiple_threads(instance_graph, "simple_iterative_tests_PLL_Indexes.txt", V + 1, weighted, thread_num);
		///*below is for reading PLL labels*/
		//vector<vector<PLL_sorted_label>> L = PLL_read_indexes_vectorFORMAT_binary("simple_iterative_tests_PLL_Indexes.txt");

		
		vector<vector<graph_hash_of_mixed_weighted_HL_PLL_label>> L = graph_hash_of_mixed_weighted_HL_PLL_generate_indexes_multiple_threads(instance_graph, V + 1, weighted, thread_num);

		graph_hash_of_mixed_weighted_HL_PLL_check_correctness_of_PLL_labels(L, instance_graph, iteration_source_times, iteration_terminal_times, boost_random_time_seed);

		long long int index_size = 0;
		for (auto it = L.begin(); it != L.end(); it++) {
			index_size = index_size + (*it).size();
		}
		avg_index_size_per_v = avg_index_size_per_v + (float)index_size / V / iteration_graph_times;

	}

	cout << "avg_index_size_per_v: " << avg_index_size_per_v << endl;
}

