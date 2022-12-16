#pragma once
#include <iostream>
#include <ThreadPool.h>
#include <shared_mutex>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID.h>
#include <graph_hash_of_mixed_weighted/HL/PSL_enhancedoriginalR2/graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_labels.h>

 /*unique code for this file: 316*/

static std::shared_mutex mtx_316;
static graph_v_of_v_idealID ideal_graph_316;
static vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>> L_316;
static vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>> L_temp_316; // d轮中d-1轮的label
static vector<int> pos_316;
static vector<int> pos_2_316;
static vector<int> increment_316;
static vector<vector<float>> T_316;
static vector<vector<bool>> dirty_tag_316;
static queue<int> Qid_316; // IDs of available elements for T
static bool if_continue_316;
static vector<int> vertexID_new_to_old_316;
static bool this_parallel_PSL_enhancedoriginalR2_is_running_316 = false;
long long int max_labal_size_316;
long long int labal_size_316;
auto begin_time_316 = std::chrono::high_resolution_clock::now();
double max_run_time_nanoseconds_316;

void update_2019R1_condition_PSL_enhancedoriginalR2(int v1, int ideal_graph_size) {
	/*here, we assume v1 and v2 have the same number of adjs*/


	for (int v2 = v1 + 1; v2 < ideal_graph_size; v2++)
	{
		/*here is a little trick. it's certain that i has adjs no less than j*/
		if (ideal_graph_316[v1].size() > ideal_graph_316[v2].size())
			break; // no need to j++ any more


		int condition;

		if (graph_v_of_v_idealID_contain_edge(ideal_graph_316, v1, v2)) { // may be equivalent_2
			bool is_equivalent_2 = true;
			int size = ideal_graph_316[v1].size();
			auto it1 = ideal_graph_316[v1].begin();
			auto it2 = ideal_graph_316[v2].begin();
			while (it1 != ideal_graph_316[v1].end() && it2 != ideal_graph_316[v2].end()) {
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
			int size = ideal_graph_316[v1].size();
			auto it1 = ideal_graph_316[v1].begin();
			auto it2 = ideal_graph_316[v2].begin();
			while (it1 != ideal_graph_316[v1].end()) {
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
			graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2[vertexID_new_to_old_316[v1]] = 11;
			graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2[vertexID_new_to_old_316[v2]] = 11;
			graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[vertexID_new_to_old_316[v2]] = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[vertexID_new_to_old_316[v1]];
		}
		else if (condition == 2)
		{
			graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2[vertexID_new_to_old_316[v1]] = 12;
			graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2[vertexID_new_to_old_316[v2]] = 12;
			graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[vertexID_new_to_old_316[v2]] = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[vertexID_new_to_old_316[v1]];
		}


	}

}

void terminate_procedures_316() {
	/*clear global values*/
	graph_v_of_v_idealID().swap(ideal_graph_316);
	vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>>().swap(L_316);
	vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>>().swap(L_temp_316);
	vector<int>().swap(pos_316);
	vector<int>().swap(pos_2_316);
	vector<int>().swap(increment_316);
	vector<vector<float>>().swap(T_316);
	vector<vector<bool>>().swap(dirty_tag_316);
	queue<int>().swap(Qid_316); // IDs of available elements for T
	vector<int>().swap(vertexID_new_to_old_316);
	this_parallel_PSL_enhancedoriginalR2_is_running_316 = false;
}


bool graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_compare(const pair<int, int>& i, pair<int, int>& j)
{
	/*< is nearly 10 times slower than >*/
	return i.second > j.second;  // < is from small to big; > is from big to small.  sort by the second item of pair<int, int>
}

void graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_thread_function_dij(int u)
{
	/*the id of T and dirty_flag uesd here */
	int used_id;
	mtx_316.lock();
	used_id = Qid_316.front();
	Qid_316.pop();
	mtx_316.unlock();

	/* save labels of u in T */ //有了这个T之后,后面的query就会快一点,这个应该不用改
	int u_label_size = L_316[u].size();
	for (int i = 0; i < u_label_size; i++)
	{
		int w = L_316[u][i].vertex;
		float dis = L_316[u][i].distance;
		if (dirty_tag_316[used_id][w]) // the same hub may have redundancy, record the shortest distance
		{
			dirty_tag_316[used_id][w] = false;
			T_316[used_id][w] = dis;
		}
		else if (dis < T_316[used_id][w])
			T_316[used_id][w] = dis;
	}

	/*
	if (reduction_measures[vertexID_new_to_old[u]] == 2) //如果u是lms中的点,那么相当于它被删除掉了,不需要再给它做label了
	{
		for (int i = 0; i < u_label_size; i++)
			dirty_tag_316[used_id][L_316[u][i].vertex] = true; // recover dirty tag
	}
	*/

	int u_adj_size = ideal_graph_316[u].size(); // u的邻居的数量
	for (int i = 0; i < u_adj_size; i++)        // 遍历u的邻居v
	{
		if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1[vertexID_new_to_old_316[ideal_graph_316[u][i].first]] == 2) //如果点v是lms集合中的点,那么进入这个循环,找到v除了u以外的邻接点(N2部分)
		{
			int v = ideal_graph_316[u][i].first;
			float ec_1 = ideal_graph_316[u][i].second;        //(u,v)距离
			int v_label_size = pos_316[v] + increment_316[v]; // label size in the last iteration
			int v_adj_size = ideal_graph_316[v].size();
			for (int j = 0; j < v_adj_size; j++)
			{
				if (ideal_graph_316[v][j].first == u) //需要排除是u点的情况
					continue;
				int v_1 = ideal_graph_316[v][j].first;     //这是点v1的index
				float ec_2 = ideal_graph_316[v][j].second; // v到v1的边的长度
				//遍历v_1的d-2级别label,用上pos_2
				for (int k = pos_2_316[v_1]; k < pos_316[v_1]; k++)
				{
					int w = L_316[v_1][k].vertex;

					if (w > u)
						continue; // cannot be the hub of u

					/* query pruning */
					float dis = L_316[v_1][k].distance + ec_1 + ec_2;
					int w_label_size = L_316[w].size();
					bool flag = false;
					for (int k = 0; k < w_label_size; k++)
						if (!dirty_tag_316[used_id][L_316[w][k].vertex])
						{
							float query_dis = T_316[used_id][L_316[w][k].vertex] + L_316[w][k].distance;
							if (query_dis - 1e-6 <= dis)
							{
								flag = true;
								break;
							}
						}
					if (flag)
						continue;

					/*add a new label*/
					graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label xx;
					// graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label xx;
					xx.vertex = w;
					xx.distance = dis;
					xx.parent_vertex = v;
					L_temp_316[u].push_back(xx);
					if (!if_continue_316)
						if_continue_316 = true;
				}
			}
		}
		else // v不是lms集合中的点,那么进入这个循环,找到v的d-1label(N1部分)
		{
			int v = ideal_graph_316[u][i].first;// visit each neighbour of u
			float ec = ideal_graph_316[u][i].second;// length of edge(u,v)
			int v_label_size = pos_316[v] + increment_316[v]; // label size in the last iteration.这个要关注一下
			for (int j = pos_316[v]; j < v_label_size; j++)
			{
				int w = L_316[v][j].vertex; // visit each hub of v
				if (w > u)
					continue; // cannot be the hub of u

				/* query pruning */
				float dis = L_316[v][j].distance + ec;
				int w_label_size = L_316[w].size();
				bool flag = false;
				for (int k = 0; k < w_label_size; k++)
					if (!dirty_tag_316[used_id][L_316[w][k].vertex])
					{
						float query_dis = T_316[used_id][L_316[w][k].vertex] + L_316[w][k].distance;
						if (query_dis - 1e-6 <= dis)
						{
							flag = true;
							break;
						}
					}
				if (flag)
					continue;

				/*add a new label*/
				graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label xx;
				// graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label xx;
				xx.vertex = w;
				xx.distance = dis;
				xx.parent_vertex = v;
				L_temp_316[u].push_back(xx);
				if (!if_continue_316)
					if_continue_316 = true;
			}
		}
	}

	for (int i = 0; i < u_label_size; i++)
		dirty_tag_316[used_id][L_316[u][i].vertex] = true; // recover dirty tag



	mtx_316.lock();
	Qid_316.push(used_id);
	labal_size_316 = labal_size_316 + L_temp_316[u].size(); // notably, there are redundent labels that may be removed in the end, so max_labal_size_316 is a limit for middle process, not only for the end
	mtx_316.unlock();

	if (labal_size_316 > max_labal_size_316) {
		throw string("labal size limit error");  // after catching error, must call terminate_procedures_316(), otherwise this PSL cannot be reused
	}
	if (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin_time_316).count() > max_run_time_nanoseconds_316) {
		throw string("run time limit error");  // after catching error, must call terminate_procedures_316(), otherwise this PSL cannot be reused
	}
}

void graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_thread_function_add_new_labels(int u)
{
	pos_2_316[u] = pos_316[u];
	pos_316[u] += increment_316[u]; increment_316[u] = L_temp_316[u].size();
	L_316[u].insert(L_316[u].end(), L_temp_316[u].begin(), L_temp_316[u].end());

	vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>().swap(L_temp_316[u]);
	vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>(L_316[u]).swap(L_316[u]);  // swap释放vector中多余空间： https://blog.csdn.net/qq_41929943/article/details/103190891 
}

bool graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_transform_labels_to_old_vertex_IDs_compare
(graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label& i, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label& j)
{
	return i.vertex < j.vertex;  // < is from small to big; > is from big to small
}

void graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_transform_labels_to_old_vertex_IDs_element(vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>>* output_L, int v_k) {

	int L_v_k_size = L_temp_316[v_k].size();
	for (int i = 0; i < L_v_k_size; i++) {
		auto it = &L_temp_316[v_k][i];
		it->vertex = vertexID_new_to_old_316[it->vertex];
		it->parent_vertex = vertexID_new_to_old_316[it->parent_vertex];
	}
	sort(L_temp_316[v_k].begin(), L_temp_316[v_k].end(), graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_transform_labels_to_old_vertex_IDs_compare);

	(*output_L)[vertexID_new_to_old_316[v_k]] = L_temp_316[v_k];
	vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>().swap(L_temp_316[v_k]); // clear new labels for RAM efficiency

}

vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>> graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_transform_labels_to_old_vertex_IDs
(int N, int max_N, int num_of_threads) {

	/*time complexity: O(V*L*logL), where L is average number of labels per vertex*/

	vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>> output_L(max_N);

	/*time complexity: O(V*L*logL), where L is average number of labels per vertex*/
	ThreadPool pool(num_of_threads);
	std::vector< std::future<int> > results; // return typename: xxx
	vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>>* p = &output_L;
	for (int v_k = 0; v_k < N; v_k++) {
		results.emplace_back(
			pool.enqueue([p, v_k] { // pass const type value j to thread; [] can be empty
			graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_transform_labels_to_old_vertex_IDs_element(p, v_k);
			return 1; // return to results; the return type must be the same with results
		})
		);
	}
	for (auto&& result : results)
		result.get(); // all threads finish here

	return output_L;
}


vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>> 
graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_generate_indexes_multiple_threads(graph_hash_of_mixed_weighted& input_graph, int max_N, int num_of_threads, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_Use_Reduction_info& reduction_method) {

	begin_time_316 = std::chrono::high_resolution_clock::now();
	max_run_time_nanoseconds_316 = reduction_method.max_run_time_seconds * 1e9;
	labal_size_316 = 0;
	max_labal_size_316 = reduction_method.max_labal_size;


	mtx_316.lock();
	if (this_parallel_PSL_enhancedoriginalR2_is_running_316 == true) {
		cout << "the following parallel PSL_enhancedoriginalR2 code cannot be run parallelly, due to the above (static) globel values" << endl;
		exit(1);
	}
	this_parallel_PSL_enhancedoriginalR2_is_running_316 = true;
	mtx_316.unlock();


	
	
	int N = input_graph.hash_of_vectors.size();

	/*sort vertices by degrees*/
	vector<pair<int, int>> sorted_vertices;
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		sorted_vertices.push_back({ it->first, input_graph.degree(it->first) });
	}
	sort(sorted_vertices.begin(), sorted_vertices.end(), graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_compare);

	/*graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID*/
	unordered_map<int, int> vertexID_old_to_new;
	vertexID_new_to_old_316.resize(N);
	for (int i = 0; i < N; i++) {
		vertexID_old_to_new[sorted_vertices[i].first] = i;
		vertexID_new_to_old_316[i] = sorted_vertices[i].first;
	}
	vector<pair<int, int>>().swap(sorted_vertices);
	ideal_graph_316 = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(input_graph, vertexID_old_to_new);



	/*redcution: add and remove certain edges*/
	graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1.clear(); // for using this function multiple times
	graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1.resize(N);
	/*clear graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1*/
	graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2.clear(); // for using this function multiple times
	graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2.resize(N);
	graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1.resize(N);
	std::iota(std::begin(graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1), std::end(graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1), 0); // Fill with 0, 1, ...


	/* redcution1: equivalence between vertices */
	if (reduction_method.use_2019R1)
	{
		auto begin = std::chrono::high_resolution_clock::now();
		int ideal_graph_size = ideal_graph_316.size();
		ThreadPool pool(num_of_threads);
		std::vector< std::future<int> > results; // return typename: xxx
		for (int i = 0; i < ideal_graph_size; i++)
		{
			results.emplace_back(
				pool.enqueue([i, ideal_graph_size] { // pass const type value j to thread; [] can be empty
				update_2019R1_condition_PSL_enhancedoriginalR2(i, ideal_graph_size);
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
			if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[i] != i)
			{
				reduction_method.reduce_V_num_2019R1++;
				graph_v_of_v_idealID_remove_all_adjacent_edges(ideal_graph_316, vertexID_old_to_new[i]);
			}
		}
		auto end = std::chrono::high_resolution_clock::now();
		reduction_method.time_2019R1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
	}

	auto begin = std::chrono::high_resolution_clock::now();
	if (reduction_method.use_2019R2) {
		reduction_method.MG_num = 0;
		for (int x = 0; x < N; x++) {
			if (ideal_graph_316[x].size() > 0) {										//Prevent memory overflow
				if (x > ideal_graph_316[x][ideal_graph_316[x].size() - 1].first) {		//Here, only one comparison is needed. A little trick.
					graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1[vertexID_new_to_old_316[x]] = 2;
					//cout << "reduce " << vertexID_new_to_old_316[x] << endl;
				}
			}
		}
	}
	else if (reduction_method.use_enhanced2019R2) { // e.g., 2019R2enhance_10
		reduction_method.MG_num = 0;
		for (int x = 0; x < N; x++) {
			if (ideal_graph_316[x].size() > 0) {										//Prevent memory overflow
				if (x > ideal_graph_316[x][ideal_graph_316[x].size() - 1].first) {		//Here, only one comparison is needed. A little trick.
					graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1[vertexID_new_to_old_316[x]] = 2;
					//cout << "reduce " << vertexID_new_to_old_316[x] << endl;
				}
			}
		}
		int bound = reduction_method.max_degree_MG_enhanced2019R2;
		for (int x = N - 1; x >= 0; x--) { // from low ranking to high ranking
			if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1[vertexID_new_to_old_316[x]] == 0 && ideal_graph_316[x].size() <= bound) { // bound is the max degree for reduction
				bool no_adj_MG_vertices = true;
				for (auto it = ideal_graph_316[x].begin(); it != ideal_graph_316[x].end(); it++) {
					if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1[vertexID_new_to_old_316[it->first]] == 2) {
						no_adj_MG_vertices = false;
						break;
					}
				}
				if (no_adj_MG_vertices) {
					graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1[vertexID_new_to_old_316[x]] = 2; // new reduction
					//cout << "new reduce " << vertexID_new_to_old_316[x] << endl;
				}
			}
		}
	}
	auto end = std::chrono::high_resolution_clock::now();
	reduction_method.time_2019R2_or_enhanced_pre = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s


	unordered_map<int, int>().swap(vertexID_old_to_new);


	L_316.resize(N); pos_316.resize(N); pos_2_316.resize(N); increment_316.resize(N);
	for (int v_k = 0; v_k < N; v_k++)  //initialization
	{
		if (ideal_graph_316[v_k].size() == 0) continue; // do not search isolated vertices
		pos_316[v_k] = 0; // record the position of the begining of the label in last iteration
		pos_2_316[v_k] = 0;
		increment_316[v_k] = 1;
		graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label xx;
		xx.vertex = v_k;
		xx.distance = 0;
		xx.parent_vertex = v_k;
		L_316[v_k].push_back(xx);
	}


	T_316.resize(num_of_threads);
	dirty_tag_316.resize(num_of_threads);
	for (int i = 0; i < num_of_threads; i++)
	{
		T_316[i].resize(N);
		dirty_tag_316[i].resize(N);
		for (int j = 0; j < N; j++) dirty_tag_316[i][j] = true;
		Qid_316.push(i);
	}

	L_temp_316.resize(N);
	if_continue_316 = true;
	int if_continue_316_false_time = 0;


	ThreadPool pool(num_of_threads);
	std::vector< std::future<int> > results;
	int num_of_threads_per_push = num_of_threads * 100; // 每次push进去 num_of_threads_per_push 线程，如果没有异常，继续push进去num_of_threads_per_push线程；如果全都一起push进去必须全部线程都结束才能catch异常
	while (if_continue_316 || if_continue_316_false_time < 2) // since R2 skip some vertices, some new labels can only be generated when d increases 2, not 1, thus terminals the loop only when if_continue_316==false twice
	{
		if (if_continue_316 == true) {
			if_continue_316_false_time = 0;
		}
		if_continue_316 = false;

		int push_num = 0;
		for (int u = 0; u < N; u++)
		{
			if (ideal_graph_316[u].size() == 0 || graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1[vertexID_new_to_old_316[u]] == 2) continue; // do not search isolated vertices
			results.emplace_back(
				pool.enqueue([u] { // pass const type value j to thread; [] can be empty
				graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_thread_function_dij(u);
				return 1; // return to results; the return type must be the same with results
			})
			);
			push_num++;
			if (push_num % num_of_threads_per_push == 0) {
				for (auto&& result : results)
					result.get(); //all threads finish here
				results.clear();
			}
		}

		for (auto&& result : results)
			result.get();
		results.clear();

		for (int u = 0; u < N; u++)
		{
			if (ideal_graph_316[u].size() == 0 || graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1[vertexID_new_to_old_316[u]] == 2) continue; // do not search isolated vertices
			results.emplace_back(
				pool.enqueue([u] { // pass const type value j to thread; [] can be empty
				graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_thread_function_add_new_labels(u);
				return 1; // return to results; the return type must be the same with results
			})
			);
		}

		for (auto&& result : results)
			result.get();
		results.clear(); // 如果本轮没有异常则继续

		if (if_continue_316 == false) {
			if_continue_316_false_time++; // if if_continue_316_false_time==2, then if_continue_316==false twice
		}
	}

	/* clean the labels */
	vector<int> pre(N);
	for (int u = 0; u < N; u++)
	{
		/* save the smallest label for each hub of u in T */
		int u_label_size = L_316[u].size();
		for (int i = 0; i < u_label_size; i++)
		{
			int w = L_316[u][i].vertex; float dis = L_316[u][i].distance;
			if (dirty_tag_316[0][w]) // the same hub may have redundancy, record the shortest distance
			{
				dirty_tag_316[0][w] = false;
				T_316[0][w] = dis;
				pre[w] = L_316[u][i].parent_vertex;
			}
			else if (dis < T_316[0][w])
			{
				T_316[0][w] = dis;
				pre[w] = L_316[u][i].parent_vertex;
			}
		}
		for (int i = 0; i < u_label_size; i++)
		{
			int w = L_316[u][i].vertex;
			if (!dirty_tag_316[0][w])
			{
				dirty_tag_316[0][w] = true;
				graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label xx;
				xx.vertex = w;
				xx.distance = T_316[0][w];
				xx.parent_vertex = pre[w];
				L_temp_316[u].push_back(xx);
			}
		}
		vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>().swap(L_316[u]); // do not have two L in RAM simultaneously
	}

	/*return unordered_map_L for old_IDs*/
	vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>> output_L = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_transform_labels_to_old_vertex_IDs(N, max_N, num_of_threads);


	terminate_procedures_316();

	return output_L;


}


#include <text mining/binary_save_read_vector_of_vectors.h>
void graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_multiple_threads
(graph_hash_of_mixed_weighted& input_graph, string index_file_name, int max_N, int num_of_threads, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_Use_Reduction_info& reduction_method) {

	auto begin1 = std::chrono::high_resolution_clock::now();

	vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>> L = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_generate_indexes_multiple_threads(input_graph, max_N, num_of_threads, reduction_method);

	auto end1 = std::chrono::high_resolution_clock::now();
	float runningtime1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1).count() / 1e9; // s

	cout << index_file_name + " PSL_enhancedoriginalR2 runningtime1: " << runningtime1 << "s" << std::endl;


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

#include <graph_hash_of_mixed_weighted/HL/PSL_enhancedoriginalR2/graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_multiple_threads.h>


int main()
{
	test_graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_multiple_threads();
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
//#include <graph_hash_of_mixed_weighted/HL/PLL_with_non_adj_reduction/graph_hash_of_mixed_weighted_HL_PLL_with_non_adj_reduction_labels.h>

void graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_check_correctness_of_PSL_enhancedoriginalR2_labels(vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>>& L, graph_hash_of_mixed_weighted& instance_graph, int iteration_source_times, int iteration_terminal_times, boost::random::mt19937& boost_random_time_seed) {

	/*below is for checking whether the above labels are right (by randomly computing shortest paths)

	this function can only be used when 0 to n-1 is in the graph, i.e., the graph is an ideal graph

	*/

	boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(instance_graph.hash_of_vectors.size() - 1) };

	//graph_hash_of_mixed_weighted_print(instance_graph);

	for (int yy = 0; yy < iteration_source_times; yy++) {
		int source = dist(boost_random_time_seed);
		std::unordered_map<int, double> distances;
		std::unordered_map<int, int> predecessors;

		//source = 3;

		graph_hash_of_mixed_weighted_shortest_paths_source_to_all(instance_graph, source, distances, predecessors);

		for (int xx = 0; xx < iteration_terminal_times; xx++) {

			int terminal = dist(boost_random_time_seed);

			//terminal = 0;

			float dis = graph_hash_of_mixed_weighted_HL_PSL_with_2019R1_extract_distance(L, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, source, terminal);
			if (abs(dis - distances[terminal]) > 1e-4 && (dis < std::numeric_limits<float>::max() || distances[terminal] < std::numeric_limits<float>::max())) {
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

			vector<pair<int, int>> path = graph_hash_of_mixed_weighted_HL_PSL_with_2019R1_extract_shortest_path(L, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, source, terminal);

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

void test_graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_multiple_threads() {

	/*parameters*/
	int iteration_graph_times = 5e1, iteration_source_times = 10, iteration_terminal_times = 10;
	int V = 1e3, E = 5e3, precision = 1, thread_num = 5;
	float ec_min = 0.1, ec_max = 1; // set ec_min=ec_max=1 for testing unweighted PSL_enhancedoriginalR2

	float avg_index_time = 0, avg_index_size_per_v = 0, avg_reduce_V_num_2019R1 = 0, avg_MG_num = 0;


	/*reduction method selection*/
	graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_Use_Reduction_info mm;
	mm.use_2019R1 = 1;
	mm.use_2019R2 = 0;
	mm.use_enhanced2019R2 = 1;
	mm.max_degree_MG_enhanced2019R2 = 20;
	mm.max_labal_size = 1e9;
	mm.max_run_time_seconds = 1e9;


	/*iteration*/
	for (int i = 0; i < iteration_graph_times; i++) {
		//cout << i << endl;

		/*input and output; below is for generating random new graph, or read saved graph*/
		int generate_new_graph = 0;
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
		vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>> L;
		try {
			L = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_generate_indexes_multiple_threads(instance_graph, V + 1, thread_num, mm);
		}
		catch (string s) {
			cout << s << endl;
			terminate_procedures_316();
			continue;
		}
		auto end = std::chrono::high_resolution_clock::now();
		double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
		avg_index_time = avg_index_time + runningtime / iteration_graph_times;

		avg_reduce_V_num_2019R1 = avg_reduce_V_num_2019R1 + (float)mm.reduce_V_num_2019R1 / iteration_graph_times;
		avg_MG_num = avg_MG_num + (float)mm.MG_num / iteration_graph_times;

		//for (int xx = 0;  xx < graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1.size(); xx++) {
		//	cout << xx << ":" << graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1[xx] << endl;
		//}
		//for (int xx = 0; xx < graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2.size(); xx++) {
		//	cout << xx << ":" << graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2[xx] << endl;
		//}
		//for (int xx = 0; xx < graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1.size(); xx++) {
		//	cout << xx << ":" << graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[xx] << endl;
		//}
		//save_visualized_HL_labels("a.txt", L);
		//getchar();

		graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_check_correctness_of_PSL_enhancedoriginalR2_labels(L, instance_graph, iteration_source_times, iteration_terminal_times, boost_random_time_seed);

		long long int index_size = 0;
		for (auto it = L.begin(); it != L.end(); it++) {
			index_size = index_size + (*it).size();
		}
		avg_index_size_per_v = avg_index_size_per_v + (float)index_size / V / iteration_graph_times;

		graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_clear_labels(L);

	}

	cout << "avg_index_time: " << avg_index_time << endl;
	cout << "avg_index_size_per_v: " << avg_index_size_per_v << endl;
	cout << "avg_reduce_V_num_2019R1: " << avg_reduce_V_num_2019R1 << endl;
	cout << "avg_MG_num: " << avg_MG_num << endl;
}
