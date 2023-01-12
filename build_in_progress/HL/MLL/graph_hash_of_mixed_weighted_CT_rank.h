#pragma once
/*this CT records ranks*/
#include <build_in_progress/HL/MLL/graph_hash_of_mixed_weighted_PLL_rank.h>
#include <build_in_progress/HL/two_hop_v2/graph_hash_of_mixed_weighted_CT_v2.h>

vector<int> global_CT_rank;			// for mll



void CT_rank(graph_hash_of_mixed_weighted& input_graph, int max_N_ID, graph_hash_of_mixed_weighted_CT_v2_case_info& case_info) {


	//--------------------------------- step 1: initialization ---------------------------
	//cout << "step 1: initialization" << endl;
	auto begin1 = std::chrono::high_resolution_clock::now();

	auto& Bags = case_info.Bags;
	auto& isIntree = case_info.isIntree;
	auto& root = case_info.root;
	auto& tree_st = case_info.tree_st;
	auto& tree_st_r = case_info.tree_st_r;
	auto& first_pos = case_info.first_pos;
	auto& lg = case_info.lg;
	auto& dep = case_info.dep;

	int N = input_graph.hash_of_vectors.size();

	// used to extract the path
	vector<vector<int>> midnode_bag(max_N_ID);
	graph_hash_of_mixed_weighted_to_ideal_graph_of_CT(input_graph, max_N_ID);
	isIntree.resize(max_N_ID, 0); // whether it is in the CT-tree

	/*
	priority_queue for maintaining the degrees of vertices  (we do not update degrees in q, so everytime you pop out a degree in q, you check whether it is the right one, ignore it if wrong)

	We need to maintain the degree of vertexes in graph dynamically, otherwise it would increase the cost in time
	*/
	priority_queue<node_degree> q;
	for (int i = 0; i < N; i++) // change N to max_N_ID cause bugs, why?
	{
		node_degree nd;
		nd.degree = global_ideal_graph_CT[i].size();
		nd.vertex = i;
		q.push(nd);
	}

	auto end1 = std::chrono::high_resolution_clock::now();
	case_info.time1_initialization = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1).count() / 1e9;
	//---------------------------------------------------------------------------------------------------




	//-------------------------------------------------- step 2: MDE-based tree decomposition ------------------------------------------------------------
	//cout << "step 2: MDE-based tree decomposition" << endl;
	auto begin2 = std::chrono::high_resolution_clock::now();


	/*MDE-based tree decomposition; generating bags*/
	int bound_lambda = N;
	Bags.resize(N);
	vector<int> right_temp(N);
	vector<int> node_order(N + 1); // merging ID to original ID

	ThreadPool pool(case_info.thread_num);
	std::vector< std::future<int> > results; // return typename: xxx

	for (int i = 1; i <= N; i++)  // repeat
	{
		node_degree nd;
		while (1)
		{
			nd = q.top();
			q.pop();
			if (!isIntree[nd.vertex] && global_ideal_graph_CT[nd.vertex].size() == nd.degree) break; // nd.vertex is the lowest degree vertex not in tree
		}
		int v_x = nd.vertex;  // the node with the minimum degree in G 
		if (nd.degree >= case_info.d) // reach the boudary
		{
			bound_lambda = i - 1;
			q.push(nd);
			break;  // until |Ni| >= d
		}

		isIntree[v_x] = 1; // add to CT-tree
		node_order[i] = v_x;

		vector<pair<int, double>>& adj_temp = global_ideal_graph_CT[v_x]; // global_ideal_graph_CT is G_i-1
		int v_adj_size = global_ideal_graph_CT[v_x].size();
		for (int j = 0; j < v_adj_size; j++)
		{
			auto it = &global_ideal_graph_CT[v_x][j];
			Bags[v_x].push_back({ it->first,it->second }); // Bags[v_x] stores adj vertices and weights of v_x
			auto it2 = &global_midnode_graph[v_x][j]; // it2->left is predecessor of v_x to global_ideal_graph_CT[v_x][j].first
			right_temp[j] = it2->right;
			midnode_bag[v_x].push_back(it2->left); // midnode_bag[v_x] stores predecessors v_x for vertices in Bags[v_x]
		}

		/*add new edge*/
		for (int j = 0; j < v_adj_size; j++) {
			int adj_j = adj_temp[j].first;
			int right_j = right_temp[j];
			for (int k = j + 1; k < v_adj_size; k++) {
				int adj_k = adj_temp[k].first;
				double new_ec = adj_temp[j].second + adj_temp[k].second;
				int right_k = right_temp[k];
				results.emplace_back(
					pool.enqueue([adj_j, adj_k, new_ec, right_j, right_k] { // pass const type value j to thread; [] can be empty
						substitute_parallel(adj_j, adj_k, new_ec, right_j, right_k);
						return 1;
						})
				);
			}
		}
		for (auto&& result : results) {
			result.get();
		}
		results.clear();


		// delete edge related to v_x and update degree;
		for (int j = 0; j < v_adj_size; j++) {
			//update degree
			nd.vertex = adj_temp[j].first;
			nd.degree = global_ideal_graph_CT[nd.vertex].size() - 1; // v_x will be removed from global_ideal_graph_CT[nd.vertex] below
			q.push(nd);
			// remove v_x
			int m = global_ideal_graph_CT[v_x][j].first;
			results.emplace_back(
				pool.enqueue([m, v_x] { // pass const type value j to thread; [] can be empty
					remove_parallel(m, v_x);
					return 1;
					})
			);
		}
		for (auto&& result : results) {
			result.get();
		}
		results.clear();
		// delete v_x from ideal graph directly
		vector<pair<int, double>>().swap(global_ideal_graph_CT[v_x]);
		vector<midnode>().swap(global_midnode_graph[v_x]);
	}


	auto end2 = std::chrono::high_resolution_clock::now();
	case_info.time2_tree_decomposition = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2).count() / 1e9;
	//--------------------------------------------------------------------------------------------------------

	//cout << "case_info.time2_tree_decomposition = " << case_info.time2_tree_decomposition << endl;


	//---------------------------------------------- step 3: generate CT-tree indexs ------------------------------------------------
	//cout << "step 3: generate CT-tree indexs" << endl;
	auto begin3 = std::chrono::high_resolution_clock::now();

	/* generate CT-tree indexs */
	vector<vector<two_hop_label_v1>> L1(N); // Labels for CT-tree index merge later, otherwise may increase query time on PLL
	vector<int> fa(N);
	root.resize(N);
	vector<double> temp_dis(N);
	vector<int> temp_parent(N);
	vector<int> order_mapping(N + 1); // original ID to merging ID
	for (int i = 1; i <= bound_lambda; i++)  order_mapping[node_order[i]] = i;
	vector<bool> popped_isIncore_node(N, 0);
	for (int i = bound_lambda + 1; i <= N; i++)  // take advantage of the priority queue
	{
		struct node_degree nd;
		while (1)
		{
			nd = q.top();
			q.pop();
			if (!isIntree[nd.vertex] && !popped_isIncore_node[nd.vertex] && global_ideal_graph_CT[nd.vertex].size() == nd.degree) break;
		}
		node_order[i] = nd.vertex; // merging ID to original ID
		popped_isIncore_node[nd.vertex] = 1;
		order_mapping[nd.vertex] = i; // original ID to merging ID
	}

	dep.resize(N);
	vector<int> islabel(N, 0);
	first_pos.resize(N);
	vector<vector<int>> son(N);
	vector<int> index_node(N);

	vector <int> isneighour(N);
	int neighbournum = 0;
	vector<double> T_temp_n(N);
	vector<int> T_temp_p(N);


	int labelnum = 0;

	for (int i = bound_lambda; i >= 1; i--)
	{
		int v_x = node_order[i]; //  node_order(N + 1); // merging ID to original ID

		fa[v_x] = INT_MAX; // merging ID
		int v_adj_size = Bags[v_x].size();
		for (int j = 0; j < v_adj_size; j++)
			if (order_mapping[Bags[v_x][j].first] < fa[v_x]) fa[v_x] = order_mapping[Bags[v_x][j].first]; // renew fa[v_x] to be the smallest merging_ID (the lowest ancestor bag)

		if (fa[v_x] > bound_lambda || v_adj_size == 0)  // a root in the forest (a bag with interfaces)
		{
			root[v_x] = v_x;  // original ID
			fa[v_x] = -1; dep[v_x] = 0;// multi_fa[v_x][0] = -1;

			/*below is the tree indexes of a root;
			Lines 24-25 of CT in 2020 paper;

			tree indexes of a root only contain interfaces, but no ancetors*/
			two_hop_label_v1 xx;
			for (int j = 0; j < v_adj_size; j++) // (bound_lambda - 1) - local distance to the interface
			{
				xx.vertex = Bags[v_x][j].first;
				xx.distance = Bags[v_x][j].second;
				xx.parent_vertex = midnode_bag[v_x][j];
				L1[v_x].push_back(xx); // interpfact parts of tree indexes of root v_x
			}
		}
		else //a non-root in the forest
		{
			fa[v_x] = node_order[fa[v_x]]; //  node_order(N + 1); // merging ID to original ID
			root[v_x] = root[fa[v_x]]; // int i = bound_lambda; i >= 1; i--, already has the right root[fa[v_x]];  root[v_x] = root[fa[v_x]], from high to low to get roots
			dep[v_x] = dep[fa[v_x]] + 1;// multi_fa[v_x][0] = fa[v_x]; for LCA
			son[fa[v_x]].push_back(v_x); // for LCA

			int index_node_num = 0;
			labelnum++;  // to ensure the order, we can not use "push_back"

			int root_adj_size = Bags[root[v_x]].size(); // iterfaces, already added above


			/*add interface node*/
			for (int j = 0; j < root_adj_size; j++)   // put the labels to interface in the beginnig position
			{  // add interface
				islabel[Bags[root[v_x]][j].first] = labelnum; // labelnum means that a hub vertex is added into bag v_x
				index_node_num++;
				index_node[index_node_num] = Bags[root[v_x]][j].first;
				temp_dis[Bags[root[v_x]][j].first] = std::numeric_limits<double>::max(); // initial dis
			}

			/*add ancestor node: representation nodes of all ancetor bags*/
			int v_y = v_x;
			while (fa[v_y] != -1)
			{
				// add ancestor
				if (islabel[fa[v_y]] != labelnum) // fa[v_y] is not in bag v_x yet
				{
					index_node_num++;
					index_node[index_node_num] = fa[v_y];
					islabel[fa[v_y]] = labelnum; // add fa[v_y] into bag v_x
					temp_dis[fa[v_y]] = std::numeric_limits<double>::max(); // initial dis
				}
				v_y = fa[v_y];
			}



			/*corresponds to Line 30 of CT: the first value after min: delta_u = Bags[v_x][j].second*/
			for (int j = 0; j < v_adj_size; j++)
			{
				// add neighbours
				if (islabel[Bags[v_x][j].first] != labelnum)
				{
					islabel[Bags[v_x][j].first] = labelnum;
					index_node.push_back(Bags[v_x][j].first);
					temp_dis[Bags[v_x][j].first] = Bags[v_x][j].second;
					temp_parent[Bags[v_x][j].first] = midnode_bag[v_x][j];
				}
				else
				{
					temp_dis[Bags[v_x][j].first] = Bags[v_x][j].second;
					temp_parent[Bags[v_x][j].first] = midnode_bag[v_x][j];
				}
			}
			// query (bound_lambda - 1)_local_distance to ancestor or the interface through neighbours


			/*corresponds to Line 30 of CT: the second min value: dis_vj + L1[vj][k].distance*/
			for (int j = 0; j < v_adj_size; j++)
				if (isIntree[Bags[v_x][j].first]) // isneighour and isintree --> can be used as an intermediate node to update labels
				{
					double dis_vj = Bags[v_x][j].second;
					int vj = Bags[v_x][j].first;
					int Lj_size = L1[vj].size();
					for (int k = 0; k < Lj_size; k++)  // update the (bound_lambda-1)_local_distance
					{
						if (islabel[L1[vj][k].vertex] == labelnum && dis_vj + L1[vj][k].distance < temp_dis[L1[vj][k].vertex])
						{
							temp_dis[L1[vj][k].vertex] = dis_vj + L1[vj][k].distance;
							temp_parent[L1[vj][k].vertex] = midnode_bag[v_x][j];
						}
					}
				}

			/*add correct indexes of Lines 29-30 of CT into L1, and possibly wrong distances for Lines 31-32 into L1*/
			// add labels to L1
			L1[v_x].resize(index_node_num);
			for (int j = 1; j <= index_node_num; j++)
			{
				two_hop_label_v1 xx;
				xx.vertex = index_node[j];
				xx.distance = temp_dis[index_node[j]];
				xx.parent_vertex = temp_parent[index_node[j]];
				L1[v_x][j - 1] = xx;

			}

			/*Lines 31-32 of CT; update possibly wrong distances for Lines 31-32 in L1*/
			// update conversely
			neighbournum++;
			for (int j = 0; j < v_adj_size; j++)
			{
				isneighour[Bags[v_x][j].first] = neighbournum;
				T_temp_n[Bags[v_x][j].first] = Bags[v_x][j].second;
				T_temp_p[Bags[v_x][j].first] = midnode_bag[v_x][j];
			}
			for (int j = 1; j <= index_node_num; j++)
			{
				int vj = index_node[j];
				int Lj_size = L1[vj].size();
				for (int k = 0; k < Lj_size; k++)
				{
					int vk = L1[vj][k].vertex;
					if ((isneighour[vk] == neighbournum) && (T_temp_n[vk] + L1[vj][k].distance < temp_dis[vj]))
					{
						temp_dis[vj] = T_temp_n[vk] + L1[vj][k].distance;
						temp_parent[vj] = T_temp_p[vk];
					}
				}
			}
			for (int j = 1; j <= index_node_num; j++)
				if (temp_dis[index_node[j]] < L1[v_x][j - 1].distance)
				{
					L1[v_x][j - 1].distance = temp_dis[index_node[j]];
					L1[v_x][j - 1].parent_vertex = temp_parent[index_node[j]];
				}

		}
	}


	/*add distance-0 labels to tree nodes; this is needed in querying functions*/
	two_hop_label_v1 node;
	for (int i = 0; i < N; i++) {
		if (isIntree[i]) {
			node.vertex = i;
			node.distance = 0;
			node.parent_vertex = i;
			L1[i].push_back(node);
		}
	}



	auto end3 = std::chrono::high_resolution_clock::now();
	case_info.time3_tree_indexs = std::chrono::duration_cast<std::chrono::nanoseconds>(end3 - begin3).count() / 1e9; // s
	//-------------------------------------------------------------------------------------------------------






	//------------------------------------------------ step 4: LCA --------------------------------------------------------------
	//cout << "step 4: LCA" << endl;
	auto begin4 = std::chrono::high_resolution_clock::now();

	/* LCA code; already get the root, the father and the depth, here is the preprocessing of querying LCA */
	int total = 0;
	vector<int> dfn(2 * N + 5);
	for (int i = 1; i <= bound_lambda; i++)
	{
		int v_x = node_order[i];
		if (root[v_x] == v_x) dfs(total, first_pos, v_x, son, dfn);
	}

	if (total > 0)
	{
		int multi_step = ceil(log(total) / log(2)) + 2;

		tree_st.resize(total + 5);
		tree_st_r.resize(total + 5);
		for (int i = 1; i <= total; i++)
		{
			tree_st[i].resize(multi_step + 2);
			tree_st_r[i].resize(multi_step + 2);
		}

		vector<int> pow_2(multi_step);

		pow_2[0] = 1;
		for (int i = 1; i < multi_step; i++)
			pow_2[i] = pow_2[i - 1] << 1;


		for (int i = 1; i <= total; i++)
		{
			tree_st[i][0] = dfn[i];
			tree_st_r[i][0] = dfn[i];
		}

		for (int j = 1; j < multi_step; j++)
			for (int i = 1; i <= total; i++)
			{
				int k = i + pow_2[j - 1];
				if (k > total || dep[tree_st[i][j - 1]] <= dep[tree_st[k][j - 1]])
					tree_st[i][j] = tree_st[i][j - 1];
				else tree_st[i][j] = tree_st[k][j - 1];
				k = i - pow_2[j - 1];
				if (k <= 0 || dep[tree_st_r[i][j - 1]] <= dep[tree_st_r[k][j - 1]])
					tree_st_r[i][j] = tree_st_r[i][j - 1];
				else  tree_st_r[i][j] = tree_st_r[k][j - 1];
			}

	}

	lg.resize(total + 1);
	for (int i = 1; i <= total; i++)
		lg[i] = floor(log(i) / log(2));



	/*clear variables not used below*/
	global_CT_rank.swap(order_mapping);
	vector<vector<int>>().swap(midnode_bag);
	priority_queue<node_degree>().swap(q);
	vector<int>().swap(node_order);
	vector<int>().swap(fa);
	vector<double>().swap(temp_dis);
	vector<int>().swap(temp_parent);
	vector<int>().swap(order_mapping);
	vector<int>().swap(islabel);
	vector<vector<int>>().swap(son);
	vector<int>().swap(index_node);
	vector <int>().swap(isneighour);
	vector<double>().swap(T_temp_n);
	vector<int>().swap(T_temp_p);
	vector<int>().swap(dfn);
	vector<int>().swap(right_temp);

	auto end4 = std::chrono::high_resolution_clock::now();
	case_info.time4_lca = std::chrono::duration_cast<std::chrono::nanoseconds>(end4 - begin4).count() / 1e9; // s
	//--------------------------------------------------------------------------------------------------------------------------




	//----------------------------------------------- step 5: 2-hop labeling -------------------------------------------
	//cout << "step 5: 2-hop labeling" << endl;
	auto begin5 = std::chrono::high_resolution_clock::now();

	/*update limits*/
	double to_date_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end4 - begin1).count() / 1e9;
	case_info.two_hop_case_info.max_run_time_seconds = case_info.max_run_time_seconds - to_date_time;
	if (case_info.two_hop_case_info.max_run_time_seconds < 0) {
		throw reach_limit_error_string_time;
	}
	long long int to_date_bit_size = case_info.compute_label_bit_size();
	case_info.two_hop_case_info.max_labal_size = (case_info.max_bit_size - to_date_bit_size) / sizeof(two_hop_label_v1); // this is slightly inaccurate, since reduction measures of R1 R2 are not counted
	if (case_info.two_hop_case_info.max_labal_size < 0) {
		throw reach_limit_error_string_MB;
	}


	/* construct 2-hop labels on core */
	/* graph_v_of_v_idealID to graph_hash_of_mixed_weighted */
	graph_hash_of_mixed_weighted hash_g;
	int size = global_ideal_graph_CT.size();
	for (int i = 0; i < size; i++) {
		int v_size = global_ideal_graph_CT[i].size();
		for (int x = 0; x < v_size; x++) {
			int j = global_ideal_graph_CT[i][x].first;
			if (i < j) {
				double ec = global_ideal_graph_CT[i][x].second;
				graph_hash_of_mixed_weighted_add_edge(hash_g, i, j, ec);
			}
		}
	}
	case_info.core_graph = hash_g;
	if (case_info.use_PLL) {
		graph_hash_of_mixed_weighted_PLL_rank(hash_g, max_N_ID, 1, case_info.thread_num, case_info.two_hop_case_info, global_CT_rank, bound_lambda);
	}
	else {
		graph_hash_of_mixed_weighted_PSL_v1(hash_g, max_N_ID, case_info.thread_num, case_info.two_hop_case_info);
	}
	//case_info.two_hop_case_info.print_L();
	auto end5 = std::chrono::high_resolution_clock::now();
	case_info.time5_core_indexs = std::chrono::duration_cast<std::chrono::nanoseconds>(end5 - begin5).count() / 1e9;
	//--------------------------------------------------------------------------------------------------------------------




	//-------------------------------------------------- step 6: postprocessing -------------------------------------------------------------------
	//cout << "step 6: postprocessing" << endl;
	auto begin6 = std::chrono::high_resolution_clock::now();

	/* merge tree_index: L1 into case_info.two_hop_case_info.L */
	for (int v_k = 0; v_k < N; v_k++)
	{
		if (L1[v_k].size() > 0) {
			vector<two_hop_label_v1>(L1[v_k]).swap(L1[v_k]);
			case_info.two_hop_case_info.L[v_k] = L1[v_k];
			vector<two_hop_label_v1>().swap(L1[v_k]);
		}
	}

	/*build merging_predecessors*/
	std::vector<std::vector<pair<int, midnode>>>& merging_predecessors = case_info.merging_predecessors;
	merging_predecessors.resize(max_N_ID);
	for (int i = 0; i < max_N_ID; i++) {
		int size = global_ideal_graph_CT[i].size();
		auto it1 = global_ideal_graph_CT[i].begin();
		auto it2 = global_midnode_graph[i].begin();
		for (int x = 0; x < size; x++, it1++, it2++) {
			int j = it1->first;
			midnode node = *it2;
			graph_hash_of_mixed_weighted_binary_operations_insert(merging_predecessors[i], j, node);
		}
	}

	auto end6 = std::chrono::high_resolution_clock::now();
	case_info.time6_post = std::chrono::duration_cast<std::chrono::nanoseconds>(end6 - begin6).count() / 1e9;
	//---------------------------------------------------------------------------------------------------------------------------------

	clear_gloval_values_CT();

	case_info.time_total = std::chrono::duration_cast<std::chrono::nanoseconds>(end6 - begin1).count() / 1e9;
}
