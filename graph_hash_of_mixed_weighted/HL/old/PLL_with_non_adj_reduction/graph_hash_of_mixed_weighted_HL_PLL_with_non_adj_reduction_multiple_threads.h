#pragma once
#include <iostream>
#include <ThreadPool.h>
#include <shared_mutex>
#include <chrono>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>
#include <boost/heap/fibonacci_heap.hpp>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID.h>
#include <graph_hash_of_mixed_weighted/HL/PLL_with_non_adj_reduction/graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_labels.h>

/*unique code for this file: 592*/


/*static global values that will be cleared automatically after use; these values cannot be passed to multiple threads as parameters*/
graph_v_of_v_idealID ideal_graph_592;
bool this_parallel_PLL_with_non_adj_reduction_is_running_592 = false;
vector<vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>> L_592;
int max_N_for_mtx_592 = 1e7;  // this is the max N to run
vector<std::shared_mutex> mtx_592(max_N_for_mtx_592);  // std::mutex has no copy or move constructor, while std::vector::resize() requires that; you cannot resize mtx;    moreover, do not change mtx to a pointer and then points to local values, it is very slow!!
queue<int> Qid_592; // IDs of available elements of P T
vector<vector<float>> P_dij_592;
vector<vector<float>> T_dij_592;
vector<vector<int>> P_bfs_592;
vector<vector<int>> T_bfs_592;
vector<int> vertexID_new_to_old_592;
long long int max_labal_size_592;
long long int labal_size_592;
auto begin_time_592 = std::chrono::high_resolution_clock::now();
double max_run_time_nanoseconds_592;


bool graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_compare(const pair<int, int>& i, pair<int, int>& j)
{
	/*< is nearly 10 times slower than >*/
	return i.second > j.second;  // < is from small to big; > is from big to small.  sort by the second item of pair<int, int>
}

struct graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_node_for_sp {
public:
	int vertex, parent_vertex;
	float priority_value;
}; // define the node in the queue
bool operator<(graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_node_for_sp const& x, graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_node_for_sp const& y) {
	return x.priority_value > y.priority_value; // < is the max-heap; > is the min heap
}
typedef typename boost::heap::fibonacci_heap<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_node_for_sp>::handle_type graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_handle_t_for_sp;


void update_2019R1_condition_PLL_with_non_adj_reduction(int v1, int ideal_graph_size) {
	/*here, we assume v1 and v2 have the same number of adjs*/

	for (int v2 = v1 + 1; v2 < ideal_graph_size; v2++)
	{
		/*here is a little trick. it's certain that i has adjs no less than j*/
		if (ideal_graph_592[v1].size() > ideal_graph_592[v2].size())
			break; // no need to j++ any more

		int condition;

		if (graph_v_of_v_idealID_contain_edge(ideal_graph_592, v1, v2)) { // may be equivalent_2
			bool is_equivalent_2 = true;
			int size = ideal_graph_592[v1].size();
			auto it1 = ideal_graph_592[v1].begin();
			auto it2 = ideal_graph_592[v2].begin();
			while (it1 != ideal_graph_592[v1].end() && it2 != ideal_graph_592[v2].end()) {
				if (it1->first == v2) {
					it1++;
					continue;
				}
				if (it2->first == v1) {
					it2++;
					continue;
				}
				if (it1->first == it2->first) {
					if (abs(it1->second - it2->second) > 1e-5) {
						is_equivalent_2 = false;
						break;
					}
				}
				else {
					is_equivalent_2 = false;
					break;
				}
				it1++;
				it2++;
			}
			if (is_equivalent_2)
				condition = 2;
			else
				condition = 0;
		}
		else { // may be equivalent_1
			bool is_equivalent_1 = true;
			int size = ideal_graph_592[v1].size();
			auto it1 = ideal_graph_592[v1].begin();
			auto it2 = ideal_graph_592[v2].begin();
			while (it1 != ideal_graph_592[v1].end()) {
				if (it1->first == it2->first) {
					if (abs(it1->second - it2->second) > 1e-5) {
						is_equivalent_1 = false;
						break;
					}
				}
				else {
					is_equivalent_1 = false;
					break;
				}
				it1++;
				it2++;
			}
			if (is_equivalent_1)
				condition = 1;
			else
				condition = 0;
		}

		if (condition == 1)
		{
			graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_2[vertexID_new_to_old_592[v1]] = 11;
			graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_2[vertexID_new_to_old_592[v2]] = 11;
			graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_f_2019R1[vertexID_new_to_old_592[v2]] = graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_f_2019R1[vertexID_new_to_old_592[v1]];
		}
		else if (condition == 2)
		{
			graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_2[vertexID_new_to_old_592[v1]] = 12;
			graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_2[vertexID_new_to_old_592[v2]] = 12;
			graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_f_2019R1[vertexID_new_to_old_592[v2]] = graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_f_2019R1[vertexID_new_to_old_592[v1]];
		}
	}

}


void terminate_procedures_592() {
	/*clear static global*/

	graph_v_of_v_idealID().swap(ideal_graph_592);
	vector<vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>>().swap(L_592);
	queue<int>().swap(Qid_592);
	vector<vector<float>>().swap(P_dij_592);
	vector<vector<float>>().swap(T_dij_592);
	vector<vector<int>>().swap(P_bfs_592);
	vector<vector<int>>().swap(T_bfs_592);
	vector<int>().swap(vertexID_new_to_old_592);
	this_parallel_PLL_with_non_adj_reduction_is_running_592 = false;
}

void graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_thread_function_dij_mixed(int v_k, int N)
{
	/*Pruned Dijkstra from vertex v_k; see Algorithm 1 in 2013 Japan SIGMOD paper*/

	mtx_592[max_N_for_mtx_592 -1].lock();
	int used_id = Qid_592.front();
	Qid_592.pop();
	mtx_592[max_N_for_mtx_592 -1].unlock();

	queue<int> P_changed_vertices, T_changed_vertices;
	vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_handle_t_for_sp> Q_handles(N);

	graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_node_for_sp node;
	boost::heap::fibonacci_heap<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_node_for_sp> Q;
	graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label xx;

	node.vertex = v_k;
	node.parent_vertex = v_k;
	node.priority_value = 0;
	Q_handles[v_k] = Q.push(node);
	P_dij_592[used_id][v_k] = 0;
	P_changed_vertices.push(v_k);

	mtx_592[v_k].lock();
	int L_v_k_size = L_592[v_k].size();
	for (int i = 0; i < L_v_k_size; i++) {
		int L_v_k_i_vertex = L_592[v_k][i].vertex;
		T_dij_592[used_id][L_v_k_i_vertex] = L_592[v_k][i].distance; //allocate T values for L_592[v_k]
		T_changed_vertices.push(L_v_k_i_vertex);
	}
	mtx_592[v_k].unlock();
	//因为v-k的标签在从自己出发的过程中不会发生改变，并且在求query的过程中每次都会用到，所以可以提前取出来放在T数组，节省后面查找的时间

	long long int new_label_num = 0;

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
			mtx_592[u].lock();
			auto L_u_size = L_592[u].size(); // a vector<PLL_with_non_adj_reduction_sorted_label>
			mtx_592[u].unlock();
			for (int i = 0; i < L_u_size; i++) {
				mtx_592[u].lock();      // put lock in for loop is very slow, but it may be the only way under Windows
				float dis = L_592[u][i].distance + T_dij_592[used_id][L_592[u][i].vertex];
				mtx_592[u].unlock();
				if (query_v_k_u > dis) { query_v_k_u = dis; }
			} //求query的值		
#else
			mtx_592[u].lock();
			auto L_u_size1 = L_592[u].size(); // a vector<PLL_with_non_adj_reduction_sorted_label>
			for (int i = 0; i < L_u_size1; i++) {
				float dis = L_592[u][i].distance + T_dij_592[used_id][L_592[u][i].vertex];   // dont know why this code does not work under Windows
				if (query_v_k_u > dis) { query_v_k_u = dis; }
			} //求query的值
			mtx_592[u].unlock();
#endif

			if (P_u_with_error < query_v_k_u) { // this is pruning

				xx.vertex = v_k;
				xx.distance = P_u;
				xx.parent_vertex = u_parent;

				mtx_592[u].lock();
				L_592[u].push_back(xx); //新增标签，并行时L_592[u]里面的标签不一定是按照vertex ID排好序的，但是因为什么query时用了T_dij_592的trick，没必要让L_592[u]里面的标签排好序
				mtx_592[u].unlock();
				new_label_num++;
				//下面是dij更新邻接点的过程，同时更新优先队列和距离
				int u_adj_size = ideal_graph_592[u].size();
				for (int i = 0; i < u_adj_size; i++) {
					int adj_v = ideal_graph_592[u][i].first; // this needs to be locked
					float ec = ideal_graph_592[u][i].second;
					if (P_dij_592[used_id][adj_v] == std::numeric_limits<float>::max()) { //尚未到达的点
						node.vertex = adj_v;
						node.parent_vertex = u;
						node.priority_value = P_u + ec;
						Q_handles[adj_v] = Q.push(node);
						P_dij_592[used_id][adj_v] = node.priority_value;
						P_changed_vertices.push(adj_v);
					}
					else {
						if (P_dij_592[used_id][adj_v] > P_u + ec) {
							node.vertex = adj_v;
							node.parent_vertex = u;
							node.priority_value = P_u + ec;
							Q.update(Q_handles[adj_v], node);
							P_dij_592[used_id][adj_v] = node.priority_value;
						}
					}
				}

				
			}


		}
	}

	while (P_changed_vertices.size() > 0) {
		P_dij_592[used_id][P_changed_vertices.front()] = std::numeric_limits<float>::max(); // reverse-allocate P values
		P_changed_vertices.pop();
	}
	while (T_changed_vertices.size() > 0) {
		T_dij_592[used_id][T_changed_vertices.front()] = std::numeric_limits<float>::max(); // reverse-allocate T values
		T_changed_vertices.pop();
	}

	mtx_592[v_k].lock();
	vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>(L_592[v_k]).swap(L_592[v_k]); // swap释放vector中多余空间： https://blog.csdn.net/qq_41929943/article/details/103190891 
	mtx_592[v_k].unlock();

	mtx_592[max_N_for_mtx_592 - 1].lock();
	Qid_592.push(used_id);
	labal_size_592 = labal_size_592 + new_label_num;
	mtx_592[max_N_for_mtx_592 - 1].unlock();

	if (labal_size_592 > max_labal_size_592) {
		throw string("labal size limit error");  // after catching error, must call terminate_procedures_592(), otherwise this PLL cannot be reused
	}

	if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time_592).count() > max_run_time_nanoseconds_592) {
		throw string("run time limit error");  // after catching error, must call terminate_procedures_592(), otherwise this PLL cannot be reused
	}
}

void graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_thread_function_bfs_mixed(int v_k, int N)
{
	mtx_592[max_N_for_mtx_592 -1].lock();
	int used_id = Qid_592.front(); // which P T elements does this thread use
	Qid_592.pop();
	mtx_592[max_N_for_mtx_592 -1].unlock();

	queue<int> P_changed_vertices, T_changed_vertices;

	queue<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_node_for_sp> Q;

	graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_node_for_sp node;
	graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label xx;

	node.vertex = v_k;
	node.parent_vertex = v_k;
	node.priority_value = 0;
	Q.push(node);
	P_bfs_592[used_id][v_k] = 0;
	P_changed_vertices.push(v_k);

	mtx_592[v_k].lock();
	int L_v_k_size = L_592[v_k].size();
	for (int i = 0; i < L_v_k_size; i++) {
		int L_v_k_i_vertex = L_592[v_k][i].vertex;
		T_bfs_592[used_id][L_v_k_i_vertex] = L_592[v_k][i].distance; //allocate T values for L_592[v_k]
		T_changed_vertices.push(L_v_k_i_vertex);
	}
	mtx_592[v_k].unlock();

	long long int new_label_num = 0;

	while (Q.size() > 0) {

		node = Q.front();
		Q.pop();
		int u = node.vertex;

		if (v_k <= u) { // this condition is not in 2013 paper, but in 2019 paper (Lemma 3.16)

			int u_parent = node.parent_vertex;
			int P_u = node.priority_value;

			float query_v_k_u = std::numeric_limits<float>::max();

			mtx_592[u].lock();
			int L_u_size = L_592[u].size();
			for (int i = 0; i < L_u_size; i++) {
				float dis = L_592[u][i].distance + T_bfs_592[used_id][L_592[u][i].vertex];		 // cannot lock mtx_592[u] in this for loop, since the locking time is very large		
				if (query_v_k_u > dis) { query_v_k_u = dis; }
			}
			mtx_592[u].unlock();

			if (P_u < query_v_k_u) { // this is pruning

				xx.vertex = v_k;
				xx.distance = P_u;
				xx.parent_vertex = u_parent;

				mtx_592[u].lock();
				L_592[u].push_back(xx); //新增标签，并行时L_592[u]里面的标签不一定是按照vertex ID排好序的，但是因为什么query时用了T_bfs_592的trick，没必要让L_592[u]里面的标签排好序
				mtx_592[u].unlock();
				new_label_num++;
				auto u_adj_size = ideal_graph_592[u].size();
				for (int i = 0; i < u_adj_size; i++) {
					int adj_v = ideal_graph_592[u][i].first;
					if (P_bfs_592[used_id][adj_v] == INT_MAX) {
						node.vertex = adj_v;
						node.parent_vertex = u;
						node.priority_value = P_u + 1;
						Q.push(node);
						P_bfs_592[used_id][adj_v] = node.priority_value;
						P_changed_vertices.push(adj_v);
					}
				}

				


			}
		}
	}

	while (P_changed_vertices.size() > 0) {
		P_bfs_592[used_id][P_changed_vertices.front()] = INT_MAX; // reverse-allocate P values
		P_changed_vertices.pop();
	}
	while (T_changed_vertices.size() > 0) {
		T_bfs_592[used_id][T_changed_vertices.front()] = INT_MAX; // reverse-allocate T values
		T_changed_vertices.pop();
	}

	mtx_592[v_k].lock();
	vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>(L_592[v_k]).swap(L_592[v_k]); // swap释放vector中多余空间： https://blog.csdn.net/qq_41929943/article/details/103190891 
	mtx_592[v_k].unlock();

	mtx_592[max_N_for_mtx_592 - 1].lock();
	Qid_592.push(used_id);
	mtx_592[max_N_for_mtx_592 - 1].unlock();

	labal_size_592 = labal_size_592 + new_label_num;
	if (labal_size_592 > max_labal_size_592) {
		throw string("labal size limit error");  // after catching error, must call terminate_procedures_592(), otherwise this PLL cannot be reused
	}

	if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time_592).count() > max_run_time_nanoseconds_592) {
		throw string("run time limit error");  // after catching error, must call terminate_procedures_592(), otherwise this PLL cannot be reused
	}
}



bool graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_transform_labels_to_old_vertex_IDs_compare
(graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label& i, graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label& j)
{
	return i.vertex < j.vertex;  // < is from small to big; > is from big to small
}

void graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_transform_labels_to_old_vertex_IDs_element(vector<vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>>* output_L, int v_k) {

	int L_v_k_size = L_592[v_k].size();
	for (int i = 0; i < L_v_k_size; i++) {
		auto it = &L_592[v_k][i];
		it->vertex = vertexID_new_to_old_592[it->vertex];
		it->parent_vertex = vertexID_new_to_old_592[it->parent_vertex];
	}
	sort(L_592[v_k].begin(), L_592[v_k].end(), graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_transform_labels_to_old_vertex_IDs_compare);

	(*output_L)[vertexID_new_to_old_592[v_k]] = L_592[v_k];
	vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>().swap(L_592[v_k]); // clear new labels for RAM efficiency

}


vector<vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>> graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_transform_labels_to_old_vertex_IDs
(int N, int max_N, int num_of_threads) {

	/*time complexity: O(V*L*logL), where L is average number of labels per vertex*/

	vector<vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>> output_L(max_N);

	/*time complexity: O(V*L*logL), where L is average number of labels per vertex*/
	ThreadPool pool(num_of_threads);
	std::vector< std::future<int> > results; // return typename: xxx
	vector<vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>>* p = &output_L;
	for (int v_k = 0; v_k < N; v_k++) {
		results.emplace_back(
			pool.enqueue([p, v_k] { // pass const type value j to thread; [] can be empty
			graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_transform_labels_to_old_vertex_IDs_element(p, v_k);
			return 1; // return to results; the return type must be the same with results
		})
		);
	}
	for (auto&& result : results)
		result.get(); // all threads finish here

	return output_L;
}



/*the following parallel PLL_with_non_adj_reduction code cannot be run parallelly, due to the above (static) globel values*/

vector<vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>> graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_generate_indexes_multiple_threads
(graph_hash_of_mixed_weighted& input_graph, int max_N, bool weighted, int num_of_threads, graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_Use_Reduction_info& reduction_method)
{
	begin_time_592 = std::chrono::high_resolution_clock::now();
	max_run_time_nanoseconds_592 = reduction_method.max_run_time_seconds * 1e9;
	labal_size_592 = 0;
	max_labal_size_592 = reduction_method.max_labal_size;

	if (max_N > max_N_for_mtx_592) {
		cout << "max_N > max_N_for_mtx_592; max_N_for_mtx_592 is too small!" << endl;
		exit(1);
	}

	mtx_592[max_N_for_mtx_592 - 1].lock();
	if (this_parallel_PLL_with_non_adj_reduction_is_running_592 == true) {
		cout << "the following parallel PLL_with_non_adj_reduction code cannot be run parallelly, due to the above (static) globel values" << endl;
		exit(1);
	}
	this_parallel_PLL_with_non_adj_reduction_is_running_592 = true;
	mtx_592[max_N_for_mtx_592 - 1].unlock();

	L_592.resize(max_N);
	int N = input_graph.hash_of_vectors.size();


	/*sorting vertices*/
	vector<pair<int, int>> sorted_vertices;
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		sorted_vertices.push_back({ it->first, input_graph.degree(it->first) });
	}
	sort(sorted_vertices.begin(), sorted_vertices.end(), graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_compare);
	unordered_map<int, int> vertexID_old_to_new;
	vertexID_new_to_old_592.resize(N);
	for (int i = 0; i < N; i++) {
		vertexID_old_to_new[sorted_vertices[i].first] = i;
		vertexID_new_to_old_592[i] = sorted_vertices[i].first;
	}
	vector<pair<int, int>>().swap(sorted_vertices);
	ideal_graph_592 = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(input_graph, vertexID_old_to_new);



	vector<pair<pair<int, int>, int>> new_edges_with_middle_v;	//Record newly added edges
	/*redcution: add and remove certain edges*/
	graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1.clear(); // for using this function multiple times
	graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1.resize(N);
	/*clear graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_f_2019R1*/
	graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_2.clear(); // for using this function multiple times
	graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_2.resize(N);
	graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_f_2019R1.resize(N); 
	std::iota(std::begin(graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_f_2019R1), std::end(graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_f_2019R1), 0); // Fill with 0, 1, ...


	/* redcution1: equivalence between vertices; we assume that vIDs in ideal_graph_592 are sorted by degree from large to small*/
	if (reduction_method.use_2019R1)
	{
		auto begin = std::chrono::high_resolution_clock::now();
		int ideal_graph_size = ideal_graph_592.size();

		ThreadPool pool(num_of_threads);
		std::vector< std::future<int> > results; // return typename: xxx
		for (int i = 0; i < ideal_graph_size; i++)
		{
			results.emplace_back(
				pool.enqueue([i, ideal_graph_size] { // pass const type value j to thread; [] can be empty
				update_2019R1_condition_PLL_with_non_adj_reduction(i, ideal_graph_size);
				return 1; // return to results; the return type must be the same with results
			})
			);
		}
		for (auto&& result : results)
			result.get(); // all threads finish here
		/* remove edges */
		reduction_method.reduce_V_num_2019R1 = 0;
		for (int i = 0; i < ideal_graph_size; i++)
		{
			if (graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_f_2019R1[i] != i)
			{
				reduction_method.reduce_V_num_2019R1++;
				graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_592, vertexID_old_to_new[i]);
			}
		}
		auto end = std::chrono::high_resolution_clock::now();
		reduction_method.time_2019R1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	}


    /*reduction 2; 用了2019 R2 enhance之后的图就是weighted，不能使用Unweighted bfs了！*/
	if (weighted == 0 && reduction_method.use_2019R2 + reduction_method.use_enhanced2019R2 > 0) {
		cout << "weighted = 1; // 用了2019 R2 enhance之后的图就是weighted，不能使用Unweighted bfs了！" << endl;
		weighted = 1;
	}
	auto begin = std::chrono::high_resolution_clock::now();
	if (reduction_method.use_2019R2) {
		reduction_method.MG_num = 0;
		for (int x = 0; x < N; x++) {
			if (ideal_graph_592[x].size() > 0) {										//Prevent memory overflow
				if (x > ideal_graph_592[x][ideal_graph_592[x].size() - 1].first) {		//Here, only one comparison is needed. A little trick.
					graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1[vertexID_new_to_old_592[x]] = 2;
					//cout << "reduce " << vertexID_new_to_old_592[x] << endl;
				}
			}
		}
		for (int x = N - 1; x >= 0; x--) {
			if (graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1[vertexID_new_to_old_592[x]] == 2) {
				/*add edge*/
				auto it1 = ideal_graph_592[x].begin();
				for (int m = ideal_graph_592[x].size() - 1; m > 0; m--)
				{
					for (int n = m - 1; n >= 0; n--)
					{
						double s_vPLUSv_t = (it1 + m)->second + (it1 + n)->second;
						int e1 = (it1 + m)->first;
						int e2 = (it1 + n)->first;
						if (s_vPLUSv_t < graph_v_of_v_idealID_edge_weight(ideal_graph_592, e1, e2)) {
							// (s,v)+(v,t) is the shorter path or there is no edge between s and t
							graph_v_of_v_idealID_add_edge(ideal_graph_592, e1, e2, s_vPLUSv_t);
							new_edges_with_middle_v.push_back({ {e1, e2}, x });
						}
					}
				}
				/*remove edge*/
				reduction_method.MG_num++;
				graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_592, x);
			}
		}
	}
	else if (reduction_method.use_enhanced2019R2) { // e.g., 2019R2enhance_10
		reduction_method.MG_num = 0;
		for (int x = 0; x < N; x++) {
			if (ideal_graph_592[x].size() > 0) {										//Prevent memory overflow
				if (x > ideal_graph_592[x][ideal_graph_592[x].size() - 1].first) {		//Here, only one comparison is needed. A little trick.
					graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1[vertexID_new_to_old_592[x]] = 2;
					//cout << "reduce " << vertexID_new_to_old_592[x] << endl;
				}
			}
		}
		int bound = reduction_method.max_degree_MG_enhanced2019R2;
		for (int x = N - 1; x >= 0; x--) { // from low ranking to high ranking
			if (graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1[vertexID_new_to_old_592[x]] == 0 && ideal_graph_592[x].size() <= bound) { // bound is the max degree for reduction
				bool no_adj_MG_vertices = true;
				for (auto it = ideal_graph_592[x].begin(); it != ideal_graph_592[x].end(); it++) {
					if (graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1[vertexID_new_to_old_592[it->first]] == 2) {
						no_adj_MG_vertices = false;
						break;
					}
				}
				if (no_adj_MG_vertices) {
					graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1[vertexID_new_to_old_592[x]] = 2; // new reduction
					//cout << "new reduce " << vertexID_new_to_old_592[x] << endl;
				}
			}
		}
		for (int x = N - 1; x >= 0; x--) {
			if (graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1[vertexID_new_to_old_592[x]] == 2) {
				/*add edge*/
				auto it1 = ideal_graph_592[x].begin();
				for (int m = ideal_graph_592[x].size() - 1; m > 0; m--)
				{
					for (int n = m - 1; n >= 0; n--)
					{
						double s_vPLUSv_t = (it1 + m)->second + (it1 + n)->second;
						int e1 = (it1 + m)->first;
						int e2 = (it1 + n)->first;
						if (s_vPLUSv_t < graph_v_of_v_idealID_edge_weight(ideal_graph_592, e1, e2)) {
							// (s,v)+(v,t) is the shorter path or there is no edge between s and t
							graph_v_of_v_idealID_add_edge(ideal_graph_592, e1, e2, s_vPLUSv_t);
							new_edges_with_middle_v.push_back({ {e1, e2}, x });
						}
					}
				}
				/*remove edge*/
				reduction_method.MG_num++;
				graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_592, x);
			}
		}
	}
	else if (reduction_method.use_non_adj_reduc_degree) { 
		reduction_method.MG_num = 0;
		int bound = reduction_method.max_degree_MG_enhanced2019R2;
		for (int x = N - 1; x >= 0; x--) { // from low ranking to high ranking
			if (graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1[vertexID_new_to_old_592[x]] == 0 && ideal_graph_592[x].size() <= bound) { // bound is the max degree for reduction
				bool no_adj_MG_vertices = true;
				for (auto it = ideal_graph_592[x].begin(); it != ideal_graph_592[x].end(); it++) {
					if (graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1[vertexID_new_to_old_592[it->first]] == 2) {
						no_adj_MG_vertices = false;
						break;
					}
				}
				if (no_adj_MG_vertices) {
					graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1[vertexID_new_to_old_592[x]] = 2; // new reduction
					//cout << "new reduce " << vertexID_new_to_old_592[x] << endl;
				}
			}
		}
		for (int x = N - 1; x >= 0; x--) {
			if (graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1[vertexID_new_to_old_592[x]] == 2) {
				/*add edge*/
				auto it1 = ideal_graph_592[x].begin();
				for (int m = ideal_graph_592[x].size() - 1; m > 0; m--)
				{
					for (int n = m - 1; n >= 0; n--)
					{
						double s_vPLUSv_t = (it1 + m)->second + (it1 + n)->second;
						int e1 = (it1 + m)->first;
						int e2 = (it1 + n)->first;
						if (s_vPLUSv_t < graph_v_of_v_idealID_edge_weight(ideal_graph_592, e1, e2)) {
							// (s,v)+(v,t) is the shorter path or there is no edge between s and t
							graph_v_of_v_idealID_add_edge(ideal_graph_592, e1, e2, s_vPLUSv_t);
							new_edges_with_middle_v.push_back({ {e1, e2}, x });
						}
					}
				}
				/*remove edge*/
				reduction_method.MG_num++;
				graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_592, x);
			}
		}
	}
	auto end = std::chrono::high_resolution_clock::now();
	reduction_method.time_2019R2_or_enhanced_pre = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	
	unordered_map<int, int>().swap(vertexID_old_to_new);

	/*seaching shortest paths*/
	ThreadPool pool(num_of_threads);
	std::vector< std::future<int> > results; // return typename: xxx
	int num_of_threads_per_push = num_of_threads * 100; // 每次push进去 num_of_threads_per_push 线程，如果没有异常，继续push进去num_of_threads_per_push线程；如果全都一起push进去必须全部线程都结束才能catch异常
	if (weighted) {
		P_dij_592.resize(num_of_threads);
		T_dij_592.resize(num_of_threads);
		for (int i = 0; i < num_of_threads; i++)
		{
			P_dij_592[i].resize(N);
			T_dij_592[i].resize(N);
			for (int j = 0; j < N; j++)
			{
				P_dij_592[i][j] = std::numeric_limits<float>::max();
				T_dij_592[i][j] = std::numeric_limits<float>::max();
			}
			Qid_592.push(i);
		}
		int push_num = 0;
		for (int v_k = 0; v_k < N; v_k++) {
			if (ideal_graph_592[v_k].size() > 0) {  // not from isolated vertices
				results.emplace_back(
					pool.enqueue([v_k, N] { // pass const type value j to thread; [] can be empty
					graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_thread_function_dij_mixed(v_k, N);
					return 1; // return to results; the return type must be the same with results
				})
				);
				push_num++;
			}
			if (push_num % num_of_threads_per_push == 0) { 
				for (auto&& result : results)
					result.get(); //all threads finish here
				results.clear();  
			}
		}
	}
	else {
		P_bfs_592.resize(num_of_threads);
		T_bfs_592.resize(num_of_threads);
		for (int i = 0; i < num_of_threads; i++)
		{
			P_bfs_592[i].resize(N);
			T_bfs_592[i].resize(N);
			for (int j = 0; j < N; j++)
			{
				P_bfs_592[i][j] = INT_MAX;
				T_bfs_592[i][j] = INT_MAX;
			}
			Qid_592.push(i);

		}
		int push_num = 0;
		for (int v_k = 0; v_k < N; v_k++) {
			if (ideal_graph_592[v_k].size() > 0) {  // not from isolated vertices
				results.emplace_back(
					pool.enqueue([v_k, N] { // pass const type value j to thread; [] can be empty
					graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_thread_function_bfs_mixed(v_k, N);
					return 1; // return to results; the return type must be the same with results
				})
				);
				push_num++;
				if (push_num % num_of_threads_per_push == 0) {
					for (auto&& result : results)
						result.get(); // all threads finish here
					results.clear();
				}
			}
		}
	}

	for (auto&& result : results)
			result.get(); //all threads finish here


	/*
		update predecessors for this non_adj_reduction,
		this update is for correct recursive direction.
	*/
	begin = std::chrono::high_resolution_clock::now();
	for (int i = new_edges_with_middle_v.size() - 1; i >= 0; i--) {
		int e1 = new_edges_with_middle_v[i].first.first;
		int e2 = new_edges_with_middle_v[i].first.second;
		int middle_k = new_edges_with_middle_v[i].second;
		/*
			why just change the labels of e1 and e2 ?
			because 'parent_vertex' stands for 'next vertex on the shortest path', so it can only be shown in e1 and e2's labels
		*/
		for (int j = L_592[e1].size() - 1; j >= 0; j--) {
			if (L_592[e1][j].parent_vertex == e2) {
				L_592[e1][j].parent_vertex = middle_k;
			}
		}
		for (int j = L_592[e2].size() - 1; j >= 0; j--) {
			if (L_592[e2][j].parent_vertex == e1) {
				L_592[e2][j].parent_vertex = middle_k;
			}
		}
	}
	end = std::chrono::high_resolution_clock::now();
	reduction_method.time_2019R2_or_enhanced_fixlabels = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	/*return L for old_IDs*/
	vector<vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>> output_L = graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_transform_labels_to_old_vertex_IDs(N, max_N, num_of_threads);


	terminate_procedures_592();

	return output_L;

}


#include <text mining/binary_save_read_vector_of_vectors.h>
void graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_multiple_threads(graph_hash_of_mixed_weighted& input_graph, string index_file_name, int max_N, bool weighted, int num_of_threads, graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_Use_Reduction_info& reduction_method) {

	auto begin1 = std::chrono::high_resolution_clock::now();

	vector<vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>> L = graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_generate_indexes_multiple_threads(input_graph, max_N, weighted, num_of_threads, reduction_method);

	auto end1 = std::chrono::high_resolution_clock::now();
	float runningtime1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1).count() / 1e9; // s

	cout << index_file_name + " PLL_with_non_adj_reduction runningtime1: " << runningtime1 << "s" << std::endl;


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

#include <graph_hash_of_mixed_weighted/HL/PLL_with_non_adj_reduction/graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_multiple_threads.h>


int main()
{
	test_graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_multiple_threads();
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
#include <Current_Memory_Consumption_of_This_Process.h>

void graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_check_correctness_of_PLL_with_non_adj_reduction_labels(vector<vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>>& L, graph_hash_of_mixed_weighted& instance_graph, int iteration_source_times, int iteration_terminal_times, boost::random::mt19937& boost_random_time_seed) {

	/*below is for checking whether the above labels are right (by randomly computing shortest paths)

	this function can only be used when 0 to n-1 is in the graph, i.e., the graph is an ideal graph

	*/

	boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(instance_graph.hash_of_vectors.size() - 1) };

	//graph_hash_of_mixed_weighted_print(instance_graph);

	for (int yy = 0; yy < iteration_source_times; yy++) {
		int source = dist(boost_random_time_seed);
		std::unordered_map<int, double> distances;
		std::unordered_map<int, int> predecessors;

		//source = 2;

		graph_hash_of_mixed_weighted_shortest_paths_source_to_all(instance_graph, source, distances, predecessors);

		for (int xx = 0; xx < iteration_terminal_times; xx++) {

			int terminal = dist(boost_random_time_seed);
		
			//terminal = 3;

			float dis = graph_hash_of_mixed_weighted_HL_PLL_with_2019R1_extract_distance
			(L, graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1, graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_2, graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_f_2019R1, instance_graph, source, terminal);
			if (abs(dis - distances[terminal]) > 1e-4 && (dis< std::numeric_limits<float>::max() || distances[terminal] < std::numeric_limits<float>::max()) ) {
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

			//cout << 0 << endl;
			//cout << source << " " << terminal << endl;
			//getchar();

			vector<pair<int, int>> path = graph_hash_of_mixed_weighted_HL_PLL_with_2019R1_extract_shortest_path(L, graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1, graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_2, graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_f_2019R1, instance_graph, source, terminal);

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
			if (abs(dis - path_dis) > 1e-4 && (dis < std::numeric_limits<float>::max() || distances[terminal] < std::numeric_limits<float>::max())) {
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

void test_graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_multiple_threads() {

	/*parameters*/
	int iteration_graph_times = 1e2, iteration_source_times = 10, iteration_terminal_times = 10;
	int V = 1e3, E = 5e3, precision = 1, thread_num = 5;
	float ec_min = 0.1, ec_max = 1; // set ec_min=ec_max=1 for testing unweighted PLL_with_non_adj_reduction

	float avg_index_time = 0, avg_index_size_per_v = 0, avg_reduce_V_num_2019R1 = 0, avg_MG_num = 0;

	bool weighted = true;
	if (ec_min == 1 && ec_max == 1) {
		weighted = false;
	}

	/*reduction method selection*/
	graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_Use_Reduction_info mm;
	mm.use_2019R1 = 1;
	mm.use_2019R2 = 1;
	mm.use_enhanced2019R2 = 0;
	mm.use_non_adj_reduc_degree = 0;
	mm.max_degree_MG_enhanced2019R2 = 100;
	mm.max_labal_size = 6e9;
	mm.max_run_time_seconds = 1e9;

	/*iteration*/
	for (int i = 0; i < iteration_graph_times; i++) {
		//cout << i << endl;

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


		auto begin = std::chrono::high_resolution_clock::now();
		vector<vector<graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_label>> L;
		try {
			L = graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_generate_indexes_multiple_threads(instance_graph, V + 1, weighted, thread_num, mm);
		}
		catch (string s) {
			cout << s << endl;
			terminate_procedures_592();
			continue;
		}
		auto end = std::chrono::high_resolution_clock::now();
		double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
		avg_index_time = avg_index_time + runningtime / iteration_graph_times;

		avg_reduce_V_num_2019R1 = avg_reduce_V_num_2019R1 + (float)mm.reduce_V_num_2019R1 / iteration_graph_times;
		avg_MG_num = avg_MG_num + (float)mm.MG_num / iteration_graph_times;

		//for (int xx = 0;  xx < graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1.size(); xx++) {
		//	cout << xx << ":" << graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_1[xx] << endl;
		//}
		//for (int xx = 0; xx < graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_2.size(); xx++) {
		//	cout << xx << ":" << graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_reduction_measures_2[xx] << endl;
		//}
		//for (int xx = 0; xx < graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_f_2019R1.size(); xx++) {
		//	cout << xx << ":" << graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_f_2019R1[xx] << endl;
		//}
		//save_visualized_HL_labels("a.txt", L);
		//getchar();


		graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_check_correctness_of_PLL_with_non_adj_reduction_labels(L, instance_graph, iteration_source_times, iteration_terminal_times, boost_random_time_seed);

		long long int index_size = 0;
		for (auto it = L.begin(); it != L.end(); it++) {
			index_size = index_size + (*it).size();
		}
		avg_index_size_per_v = avg_index_size_per_v + (float)index_size / V / iteration_graph_times;

		graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_clear_labels(L);

		//cout << Current_Memory_Consumption_of_This_Process() << endl;

	}

	cout << "avg_index_time: " << avg_index_time << endl;
	cout << "avg_index_size_per_v: " << avg_index_size_per_v << endl;
	cout << "avg_reduce_V_num_2019R1: " << avg_reduce_V_num_2019R1 << endl;
	cout << "avg_MG_num: " << avg_MG_num << endl;
}



