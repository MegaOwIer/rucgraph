#pragma once
#include <build_in_progress/HL/MLL/graph_hash_of_mixed_weighted_CT_rank.h>
#include <cfloat>

/*global values*/
vector<vector<pair<int, int>>> MLL; // for mll
double MLL_time_check_rank_and_monotonicity = 0; 
double MLL_time_query_dist_in_monotonicity = 0; 
double MLL_time_update_part = 0; 
double MLL_time_update_part_1 = 0; 
double MLL_time_later_update = 0;

void print_mll_time() {
	cout << "check_rank_and_monotonicity: \t" << MLL_time_check_rank_and_monotonicity << "s" << endl;
	cout << "query_dist_in_monotonicity: \t" << MLL_time_query_dist_in_monotonicity << "s" << endl;
	cout << "update_part: \t" << MLL_time_update_part << "s" << endl;
	cout << "later_update: \t" << MLL_time_later_update << "s" << endl;
}

long long int compute_mll_label_bit_size()
{
	long long int size = 0;

	for (auto it1 = MLL.begin(); it1 != MLL.end(); it1++)
		for (auto it2 = it1->begin(); it2 != it1->end(); it2++)
		{
			if (it2->second == -1)
				size += sizeof(int);
			else
				size += sizeof(pair<int, int>);
		}
	return size;
}




struct dij_node
{
public:
	int vertex, parent_vertex;
	int highest_rank;
	double distance;
}; // define the node in the queue
bool operator<(dij_node const &x, dij_node const &y)
{
	return x.distance > y.distance; // < is the max-heap; > is the min heap
}
typedef typename boost::heap::fibonacci_heap<dij_node>::handle_type dij_pointer;




/**
 * Time Complexity:	O(1) 
 *
 */
int max_ranked_vertex(vector<int> &rank, int u, int hu, int hw)
{
	int r_hu, r_hw, r_u;
	if (hu < 0)
		r_hu = -1;
	else
		r_hu = rank[hu];
	if (hw < 0)
		r_hw = -1;
	else
		r_hu = rank[hu];
	r_u = rank[u];

	if (r_u > r_hu)
	{
		if (r_u > r_hw)
			return u;
		else
			return hw;
	}
	else // r_hu is bigger
	{
		if (r_hu > r_hw)
			return hu;
		else
			return hw;
	}
}


double query_dist_in_monotonicity = 0;


/**
 * Function: 		check if the shortest path between v and i is a monotonous path
 * Time Complexity:	O( delta*(ct) )
 *					delta:	maximum in label L(v)
 *					(ct):	the time comlexity of CT_extract_distance
 */
bool check_monotonicity(graph_hash_of_mixed_weighted_CT_v2_case_info &case_info,
						int v, int i, double dist,
						vector<vector<double>> &Dist_CT,
						int print) // rank[v] < rank[i]
{
	/* push_back is a stupid operation to the extreme */
	if (v == i)
		return false;

	auto& L_v = case_info.two_hop_case_info.L[v];
	auto& is_tree = case_info.isIntree;
	auto& Bag_v = case_info.Bags[v];
	auto& rank = global_CT_rank;

	int f = 0; /* f means if i is in v */
	vector<int> Lu;
	if (!is_tree[v])
	{
		// if (print)
		// 	cout << "v is in core" << endl;
		// for (int i = L[v].size(); i > 0; i--)
		// {
		// 	Lu.push_back(L[v][i - 1].vertex);
		// }
		if (L_v.size() > 0)
		{
			for (auto it2 = L_v.begin(); it2 != L_v.end(); it2++)
			{
				if (it2->vertex == i)
				{
					f = 1;
					continue;
				}
				if (it2->vertex == v)
					continue;

				auto begin2 = std::chrono::high_resolution_clock::now();

				double dist1, dist2;
				if (Dist_CT[it2->vertex][v] == DBL_MAX)
					Dist_CT[it2->vertex][v] = Dist_CT[v][it2->vertex] = dist1 = CT_extract_distance(case_info, it2->vertex, v);
				else 
					dist1 = Dist_CT[it2->vertex][v];
				if (Dist_CT[it2->vertex][i] == DBL_MAX)
					Dist_CT[it2->vertex][i] = Dist_CT[i][it2->vertex] = dist2 = CT_extract_distance(case_info, it2->vertex, i);
				else
					dist2 = Dist_CT[it2->vertex][i];
				// double dist1 = CT_extract_distance(case_info, it2->vertex, v); // dist(u,w)
				// double dist2 = CT_extract_distance(case_info, it2->vertex, i); // dist(v,w)

				auto end2 = std::chrono::high_resolution_clock::now();
				double runningtime2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2).count() / 1e9;
				query_dist_in_monotonicity += runningtime2;

				double rlt = abs(dist1 + dist2 - dist);
				/* so far, we can make sure w(*it2) is in L[v], but we can not promise w is in L[i] */
				if (rlt < 1e-4)
				{
					// if (if_w_in_Lv(L, i, *it2))
					if (print)
						cout << it2->vertex << endl;
					return false;
				}
			}
		}
	}
	else
	{
		if (Bag_v.size() > 0)
		{
			for (auto it2 = Bag_v.begin(); it2 != Bag_v.end(); it2++)
			{
				if (it2->first == i)
				{
					f = 1;
					continue;
				}
				if (it2->first == v)
					continue;

				auto begin2 = std::chrono::high_resolution_clock::now();
				
				double dist1, dist2;
				if (Dist_CT[it2->first][v] == DBL_MAX)
					Dist_CT[it2->first][v] = Dist_CT[v][it2->first] = dist1 = CT_extract_distance(case_info, it2->first, v);
				else 
					dist1 = Dist_CT[it2->first][v];
				if (Dist_CT[it2->first][i] == DBL_MAX)
					Dist_CT[it2->first][i] = Dist_CT[i][it2->first] = dist2 = CT_extract_distance(case_info, it2->first, i);
				else
					dist2 = Dist_CT[it2->first][i];
				// double dist1 = CT_extract_distance(case_info, it2->first, v); // dist(u,w)
				// double dist2 = CT_extract_distance(case_info, it2->first, i); // dist(v,w)

				auto end2 = std::chrono::high_resolution_clock::now();
				double runningtime2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2).count() / 1e9;
				query_dist_in_monotonicity += runningtime2;

				double rlt = abs(dist1 + dist2 - dist);
				/* so far, we can make sure w(*it2) is in L[v], but we can not promise w is in L[i] */
				if (rlt < 1e-4)
				{
					return false;
				}
			}
		}
	}

	if (!f)
		return false;
	else
		return true;
}

/**
 * compared to the VLDBJ paper Algorithm3
 * i stands for v
 * v stands for u
 * vv stands for w
 *
 * Time Complexity: O( nlgn + m + n*delta*(ct) )
 */
void ct_mll_index_construction(graph_hash_of_mixed_weighted &input_graph, int max_N_ID,
							   graph_hash_of_mixed_weighted_CT_v2_case_info &case_info)
{
	int N = max_N_ID - 1;
	auto& rank = global_CT_rank;
	MLL.resize(N);

	int debug = 0;
	pair<int, int> iv = {9, 4};
	
	double check_rank_and_monotonicity = 0;
	query_dist_in_monotonicity = 0;
	double update_part = 0;
	double update_part_1 = 0;
	double later_update = 0;

	vector<vector<double>> Dist_CT(N, vector<double>(N, DBL_MAX));

	/*--- n times dijkstra ---*/
	for (int i = 0; i < N; i++)
	{
		boost::heap::fibonacci_heap<dij_node> Q;
		vector<dij_pointer> Q_pointer(N);
		vector<int> visited(N, 0);
		vector<int> InQ(N, 0);
		vector<double> dij_dist(N, std::numeric_limits<double>::max());
		vector<int> dij_parent(N, -1);
		vector<int> parent_changed_vertex;
		vector<int> h(N, -1); // h(x) records the max ranked inner vertex between i and x
		dij_node node;
		node.vertex = i;
		node.parent_vertex = i;
		node.distance = 0;
		node.highest_rank = rank[node.vertex];
		Q_pointer[i] = Q.push(node);
		dij_dist[i] = 0;
		InQ[i] = 1;

		int pp = 0;
		int print = 0;

		while (Q.size() > 0)
		{
			auto begin1 = std::chrono::high_resolution_clock::now();

			node = Q.top();		/* extract-min */
			Q.pop();			/* delete */
			int v = node.vertex;
			int v_higest_rank = node.highest_rank;
			visited[v] = 1;

			/* check the rank */
			/* check the monotonicity */
			/* i is the source, v is the traversed vertex */
			if (rank[i] >= rank[v])
			{
				if (check_monotonicity(case_info, v, i, dij_dist[v], Dist_CT, print))
				{
					MLL[v].push_back({i, h[v]});
				}
			}

			auto end1 = std::chrono::high_resolution_clock::now();
			double runningtime1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1).count() / 1e9;
			check_rank_and_monotonicity += runningtime1;

			/* update dij_dist and hw */
			auto begin2 = std::chrono::high_resolution_clock::now();

			vector<pair<int, double>> v_adj = input_graph.adj_v_and_ec(v);
			int v_adj_size = v_adj.size();
			bool update_dist;
			bool update_parent;
			double dist_tmp;
			double tmp_123;
			int tmp_124, tmp_125, tmp_126;
			for (int j = 0; j < v_adj_size; j++)
			{
				int vv = v_adj[j].first; /* vv is the neighbor of v */
				update_dist = 0;
				update_parent = 0;
				dist_tmp = dij_dist[v] + v_adj[j].second;
				tmp_123 = dist_tmp - dij_dist[vv];

				if (visited[vv] == 0)
				{
					auto begin4 = std::chrono::high_resolution_clock::now();

					if (tmp_123 < -1e-4)
					{ /* just update the dij_dist */
						update_dist = 1;
						dij_dist[vv] = dist_tmp;
						dij_parent[vv] = v;
					}
					else if (abs(tmp_123) < 1e-4)
					{ /* the same sp, but we need to find the max-ranking parent vertex */
						if (dij_parent[vv] >= 0)
						{
							tmp_124 = (rank[h[vv]] - rank[h[v]]);
							tmp_125 = (rank[h[vv]] - rank[v]);
							if (tmp_124 < 0 || tmp_125 < 0)
							{ /* v is not included in h[v], so we need to add v */
								update_parent = 1;
								dij_parent[vv] = v;
							}
						}
					}

					node.vertex = vv;
					node.distance = dij_dist[vv];
					node.parent_vertex = v;

					if (InQ[vv] == 0)
					{
						Q_pointer[vv] = Q.push(node);
						InQ[vv] = 1;
					} /* the 'else' means vv is already in Q */
					else if (update_dist)
					{
						Q.update(Q_pointer[vv], node);
					}

					if ((update_dist || update_parent) && i != v) // update_dist==1 means v is on the sp of i to vv
					{
						h[vv] = max_ranked_vertex(rank, v, h[v], h[vv]);
					}

					auto end4 = std::chrono::high_resolution_clock::now();
					double runningtime4 = std::chrono::duration_cast<std::chrono::nanoseconds>(end4 - begin4).count() / 1e9;
					update_part_1 += runningtime4;
				}
				else /* vv is already visited, but we need to check whether update parent */
				{
					if (abs(tmp_123) < 1e-4)
					{
						if (dij_parent[vv] >= 0)
						{ /* h() means the max-rank inner vertex */
							tmp_126 = rank[h[dij_parent[vv]]] - rank[h[v]];
							if (tmp_126 < 0)
							{ /* since vv is already visited, so the changed parent can't be changed in MLL, we must later update */
								update_parent = 1;
								dij_parent[vv] = v;
								parent_changed_vertex.push_back(vv);
							}
						}
					}

					if (update_parent && i != v) // update_dist==1 means v is on the sp of i to vv
					{
						h[vv] = max_ranked_vertex(rank, v, h[v], h[vv]);
					}
				}
			}

			auto end2 = std::chrono::high_resolution_clock::now();
			double runningtime2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2).count() / 1e9;
			update_part += runningtime2;
		} // while end
		int size = parent_changed_vertex.size();

		auto begin3 = std::chrono::high_resolution_clock::now();

		for (int k = 0; k < size; k++)
		{
			int v = parent_changed_vertex[k];
			if (rank[i] >= rank[v])
			{
				if (check_monotonicity(case_info, v, i, dij_dist[v], Dist_CT, 0))
				{
					int MLL_v_size = MLL[v].size();
					for (int j = 0; j < MLL_v_size; j++)
					{
						if (MLL[v][j].first == i)
						{
							MLL[v][j].second = h[v];
							break;
						}
					}
				}
			}
		}

		auto end3 = std::chrono::high_resolution_clock::now();
		double runningtime3 = std::chrono::duration_cast<std::chrono::nanoseconds>(end3 - begin3).count() / 1e9;
		later_update += runningtime3;

		if (0)
		{
			int dij_size = dij_dist.size();
			for (int j = 0; j < dij_size; j++)
			{
				double real = CT_extract_distance(case_info, i, j);
				double rslt = abs(real - dij_dist[j]);

				if (rslt > 1e-4)
				{
					cout << i << " " << j << ":";
					cout << dij_dist[j] << " " << real;
					cout << "|error|" << endl;
				}
			}
			// cout<<endl;
		}
	}


	MLL_time_check_rank_and_monotonicity = check_rank_and_monotonicity;
	MLL_time_query_dist_in_monotonicity = query_dist_in_monotonicity;
	MLL_time_update_part = update_part;
	MLL_time_update_part_1 = update_part_1;
	MLL_time_later_update = later_update;

	if (0)
	{
		for (int i = MLL.size(); i > 0; i--)
		{
			cout << i - 1 << ":\t";
			for (int j = MLL[i - 1].size(); j > 0; j--)
				cout << MLL[i - 1][j - 1].first << "," << MLL[i - 1][j - 1].second << " ";
			cout << endl;
		}
	}
}


/*
	x is the highest-rank inner vertex between u and v,
	this func is to unfold the path between u and v.
	We've already promised that u is not a adj of v.
 */
void unfold(graph_hash_of_mixed_weighted &input_graph,
			graph_hash_of_mixed_weighted_CT_v2_case_info &case_info,
			int u, int v, int x,
			vector<pair<int, int>> &path)
{

	double edge_ux = graph_hash_of_mixed_weighted_edge_weight(input_graph, u, x);
	double edge_vx = graph_hash_of_mixed_weighted_edge_weight(input_graph, v, x);
	double dist_ux = CT_extract_distance(case_info, u, x);
	double dist_vx = CT_extract_distance(case_info, v, x);
	double edge_reduc_dist_ux = abs(edge_ux - dist_ux);
	double edge_reduc_dist_vx = abs(edge_vx - dist_vx);
	int MLL_x_size = MLL[x].size();

	if (edge_reduc_dist_ux < 1e-4)
		path.push_back({u, x});
	else /* unfold path(u,x), u is a hub of x, so find hx */
	{
		int hx;
		for (int i = 0; i < MLL_x_size; i++)
		{
			if (MLL[x][i].first == u)
			{
				hx = MLL[x][i].second;
				break;
			}
		}
		unfold(input_graph, case_info, u, x, hx, path);
	}

	if (edge_reduc_dist_vx < 1e-4)
		path.push_back({v, x});
	else
	{
		int hx;
		for (int i = 0; i < MLL_x_size; i++)
		{
			if (MLL[x][i].first == v)
			{
				hx = MLL[x][i].second;
				break;
			}
		}
		unfold(input_graph, case_info, v, x, hx, path);
	}
}

void MLL_extract_path(graph_hash_of_mixed_weighted &input_graph,
					  graph_hash_of_mixed_weighted_CT_v2_case_info &case_info,
					  int s, int t, vector<pair<int, int>> &path)
{
	if (s == t)
		return;

	auto &rank = global_CT_rank;

	/* make sure r(s) < r(t), so we can always check MLL[s] */
	if (rank[s] > rank[t])
	{
		int tmp = s;
		s = t;
		t = tmp;
	}
	// cout << "s,t = " << s << " " << t << endl;

	/* check if (s,t) is connected */
	double dist_st = CT_extract_distance(case_info, s, t);
	double is_disconnected = abs(dist_st - std::numeric_limits<double>::max());
	if (is_disconnected < 1e-4)
	{
		path.push_back({INT_MAX, INT_MAX});
		return;
	}
	double edge_st = graph_hash_of_mixed_weighted_edge_weight(input_graph, s, t);
	double edge_reduc_dist_st = abs(dist_st - edge_st);
	if (edge_reduc_dist_st < 1e-4)
	{
		path.push_back({s, t});
		return;
	}

	/* use MLL to seek the path */
	int MLL_s_size = MLL[s].size();
	int w, hs;
	double dist_sw, dist_tw;
	for (int i = 0; i < MLL_s_size; i++)
	{
		w = MLL[s][i].first;
		hs = MLL[s][i].second;

		dist_sw = CT_extract_distance(case_info, s, w);
		dist_tw = CT_extract_distance(case_info, t, w);
		double rlt = abs(dist_st - dist_sw - dist_tw);
		if (rlt < 1e-4)
			break;
	}
	// cout << "w=" << w << endl;
	// cout << "dist_sw=" << dist_sw << endl;

	double edge_sw = graph_hash_of_mixed_weighted_edge_weight(input_graph, s, w);
	double edge_reduc_dist_sw = abs(edge_sw - dist_sw);
	if (edge_reduc_dist_sw > 1e-4)
	{	/* (s,w) is not sp */
		// cout << "unfold:" << s << " " << w << " " << hs << endl;
		unfold(input_graph, case_info, s, w, hs, path);
	}
	else
	{	/* (s,w) is sp */
		// cout << "(s,w) is sp" << endl;
		path.push_back({s, w});
	}

	MLL_extract_path(input_graph, case_info, w, t, path);
}

