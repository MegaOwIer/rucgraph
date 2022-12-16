#pragma once

class graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label {
public:
	int vertex, parent_vertex; // std::vector<PSL_enhancedoriginalR2_sorted_label> should be sorted by vertex
	float distance;
};

class graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_Use_Reduction_info {
public:
	bool use_2019R1 = true, use_2019R2 = false, use_enhanced2019R2 = false, use_non_adj_reduc_degree = false;
	int max_degree_MG_enhanced2019R2 = 100;
	int reduce_V_num_2019R1 = 0, MG_num = 0;
	double time_2019R1 = 0, time_2019R2_or_enhanced_pre = 0; // s
	long long int max_labal_size = 1e12;
	double max_run_time_seconds = 1e12; // s
};

/*static global values that will not be cleared automatically after use*/
static vector<int> graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1; // for 2019 R2
static vector<int> graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2; // for 2019 R1
static vector<int> graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1; // for 2019 R1





/*the following function is for clearing labels and reduction info*/
void graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_clear_labels(vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>>& L) {
	vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>>().swap(L);
	vector<int>().swap(graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_1);
	vector<int>().swap(graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2);
	vector<int>().swap(graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1);
}



#include <text mining/binary_save_read_vector_of_vectors.h>
vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>> graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_read_indexes_binary(string index_file_name) {

	vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>> L;

	binary_read_vector_of_vectors(index_file_name, L);

	return L;

}


/* -------non_adj_reduction------- */
/*query dis and path after 2019R1*/

float graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance_no_reduc(vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>>& L, int source, int terminal) {

	/*return std::numeric_limits<float>::max() is not connected*/

	if (source == terminal) {
		return 0;
	}

	float distance = std::numeric_limits<float>::max(); // if disconnected, return this large value
	auto vector1_check_pointer = L[source].begin();
	auto vector2_check_pointer = L[terminal].begin();
	auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
	while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end) {
		if (vector1_check_pointer->vertex == vector2_check_pointer->vertex) {
			float dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
			if (distance > dis) {
				distance = dis;
			}
			vector1_check_pointer++;
		}
		else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex) {
			vector2_check_pointer++;
		}
		else {
			vector1_check_pointer++;
		}
	}

	return distance;

}

float graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance
(vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>>& L, vector<int>& reduction_measures, vector<int>& graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, graph_hash_of_mixed_weighted& instance_graph, int source, int terminal)
{

	/* we assume that source and terminal are not reduced by 2019R1*/

	if (source == terminal) {
		return 0;
	}

	vector<float> selected_distance = { std::numeric_limits<float>::max() }; // Store all path lengths to be selected; disconnected = std::numeric_limits<float>::max()
	vector<pair<int, double>> s_adj = instance_graph.adj_v_and_ec(source);
	vector<pair<int, double>> t_adj = instance_graph.adj_v_and_ec(terminal);
	if (reduction_measures[source] == 2)
	{
		if (reduction_measures[terminal] == 2)
		{
			/*"Both reduced"*/
			for (auto it1 = s_adj.begin(); it1 != s_adj.end(); it1++)
			{
				for (auto it2 = t_adj.begin(); it2 != t_adj.end(); it2++)
				{
					if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[it1->first] == it1->first && graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[it2->first] == it2->first) {
						float x = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance_no_reduc(L, it1->first, it2->first);
						if (x == std::numeric_limits<float>::max()) {
							return x;
						}
						else {
							selected_distance.push_back(x + float(it1->second) + float(it2->second));
						}
					}
				}
			}
		}
		else
		{
			/*"Only source reduced"*/
			for (auto it1 = s_adj.begin(); it1 != s_adj.end(); it1++)
			{
				if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[it1->first] == it1->first) { 
					float x = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance_no_reduc(L, it1->first, terminal);
					if (x == std::numeric_limits<float>::max()) {
						return x;
					}
					else {
						selected_distance.push_back(x + float(it1->second));
					}
				}
			}
		}
	}
	else
	{
		if (reduction_measures[terminal] == 2)
		{
			/*"Only terminal reduced"*/
			for (auto it2 = t_adj.begin(); it2 != t_adj.end(); it2++)
			{
				if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[it2->first] == it2->first) { 
					float x = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance_no_reduc(L, source, it2->first);
					if (x == std::numeric_limits<float>::max()) {
						return x;
					}
					else {
						selected_distance.push_back(x + float(it2->second));
					}
				}
			}
		}
		else
		{
			/*"Nothing happened"*/
			selected_distance.push_back(graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance_no_reduc(L, source, terminal));
		}
	}

	float dis = *min_element(selected_distance.begin(), selected_distance.end());

	return dis;
}

vector<pair<int, int>> graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_shortest_path
(vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>>& L, vector<int>& reduction_measures, vector<int>& graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, graph_hash_of_mixed_weighted& instance_graph, int source, int terminal)
{

	/* we assume that source and terminal are not reduced by 2019R1*/

	vector<pair<int, int>> paths;
	if (source == terminal)
	{
		return paths;
	}

	float min_dis = std::numeric_limits<float>::max();
	vector<pair<int, int>> partial_edges(2);

	vector<pair<int, double>> s_adj = instance_graph.adj_v_and_ec(source);
	vector<pair<int, double>> t_adj = instance_graph.adj_v_and_ec(terminal);

	if (reduction_measures[source] == 2)
	{
		if (reduction_measures[terminal] == 2)
		{
			/*"Both reduced"*/
			for (auto it1 = s_adj.begin(); it1 != s_adj.end(); it1++)
			{
				for (auto it2 = t_adj.begin(); it2 != t_adj.end(); it2++)
				{
					if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[it1->first] == it1->first && graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[it2->first] == it2->first) { 
						float x = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance_no_reduc(L, it1->first, it2->first) + float(it1->second) + float(it2->second);
						/*After removing the two edges, it becomes the fourth case: nothing happened*/
						if (min_dis > x)
						{
							min_dis = x;
							partial_edges[0] = { it1->first, source };
							partial_edges[1] = { it2->first, terminal };
						}
					}				
				}
			}

			if (min_dis == std::numeric_limits<float>::max())
			{ // disconnected
				return paths;
			}

			// add partial_edges into paths
			paths.push_back(partial_edges[0]);
			paths.push_back(partial_edges[1]);

			// find new edges
			vector<pair<int, int>> new_edges = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_shortest_path(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, partial_edges[0].first, partial_edges[1].first);
			if (new_edges.size() > 0)
			{
				for (int i = new_edges.size() - 1; i >= 0; i--)
				{
					paths.push_back(new_edges[i]);
				}
			}
		}
		/*"Only source reduced"*/
		else
		{
			for (auto it1 = s_adj.begin(); it1 != s_adj.end(); it1++)
			{
				if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[it1->first] == it1->first) { 
					float x = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance_no_reduc(L, it1->first, terminal) + float(it1->second);
					if (min_dis > x)
					{
						min_dis = x;
						partial_edges[0] = { it1->first, source };
					}
				}
			}
			if (min_dis == std::numeric_limits<float>::max())
			{ // disconnected
				return paths;
			}

			// add partial_edges into paths
			paths.push_back(partial_edges[0]);

			// find new edges
			vector<pair<int, int>> new_edges = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_shortest_path(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, partial_edges[0].first, terminal);
			if (new_edges.size() > 0)
			{
				for (int i = new_edges.size() - 1; i >= 0; i--)
				{
					paths.push_back(new_edges[i]);
				}
			}
		}
	}
	else
	{
		/*"Only terminal reduced"*/
		if (reduction_measures[terminal] == 2)
		{
			for (auto it2 = t_adj.begin(); it2 != t_adj.end(); it2++)
			{
				if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[it2->first] == it2->first) { 
					float x = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance_no_reduc(L, source, it2->first) + float(it2->second);
					if (min_dis > x)
					{
						min_dis = x;
						partial_edges[0] = { it2->first, terminal };
					}
				}
			}

			if (min_dis == std::numeric_limits<float>::max())
			{ // disconnected
				return paths;
			}

			// add partial_edges into paths
			paths.push_back(partial_edges[0]);

			// find new edges
			vector<pair<int, int>> new_edges = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_shortest_path(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, source, partial_edges[0].first);
			if (new_edges.size() > 0)
			{
				for (int i = new_edges.size() - 1; i >= 0; i--)
				{
					paths.push_back(new_edges[i]);
				}
			}
		}
		/*
		"Nothing happened"
		In this case, the problem that the removed vertices appear in the path needs to be solved
		*/
		else
		{
			int vector1_capped_v_parent = 0, vector2_capped_v_parent = 0;
			float distance = std::numeric_limits<float>::max(); // if disconnected, retun this large value
			bool connected = false;
			auto vector1_check_pointer = L[source].begin();
			auto vector2_check_pointer = L[terminal].begin();
			auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
			{
				if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
				{
					connected = true;
					float dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
					if (distance > dis)
					{
						distance = dis;
						vector1_capped_v_parent = vector1_check_pointer->parent_vertex;
						vector2_capped_v_parent = vector2_check_pointer->parent_vertex;
					}
					vector1_check_pointer++;
				}
				else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex)
				{
					vector2_check_pointer++;
				}
				else
				{
					vector1_check_pointer++;
				}
			}

			if (connected)
			{
				/*the following code will not induce redundent edge, since for each vertex v_k, there is a label (v_k,0,v_k) in L[v_k];
				you must ascending from both source and terminal, otherwise you may not be able to extract SP */
				if (source != vector1_capped_v_parent)
				{
					paths.push_back({ source, vector1_capped_v_parent });
					source = vector1_capped_v_parent; // ascending from source
				}
				if (terminal != vector2_capped_v_parent)
				{
					paths.push_back({ terminal, vector2_capped_v_parent });
					terminal = vector2_capped_v_parent; // ascending from terminal
				}
			}
			else
			{
				return paths;
			}

			// find new edges
			vector<pair<int, int>> new_edges = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_shortest_path(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, source, terminal);

			if (new_edges.size() > 0)
			{
				for (int i = new_edges.size() - 1; i >= 0; i--)
				{
					paths.push_back(new_edges[i]);
				}
			}
		}
	}

	return paths;
}



/* -------2019R1------- 此处函数与PSL_with_non_adj_reduction里面的重名函数构成重载关系*/

/*query dis*/
float graph_hash_of_mixed_weighted_HL_PSL_with_2019R1_extract_distance
(vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>>& L, vector<int>& reduction_measures, vector<int>& graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592, vector<int>& graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, graph_hash_of_mixed_weighted& instance_graph, int source, int terminal)
{

	if (source == terminal) {
		return 0;
	}

	if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592[source] == 11)
	{
		/* case 2 */
		if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592[terminal] == 11)
		{
			if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source] == graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal])
			{
				pair<int, double> s_min_adj = instance_graph.min_adj(source);
				return float(s_min_adj.second * 2);
			}
			else
			{
				return graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal]);
			}
		}
		/* case 3 */
		else if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592[terminal] == 12)
		{
			return graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal]);
		}
		else { /* case 1 */
			return graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], terminal);
		}	
	}
	else if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592[source] == 12)
	{
		/* case 5 */
		if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592[terminal] == 11)
		{
			return graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal]);
		}
		/* case 6 -- same with case 3 */
		else if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592[terminal] == 12)
		{
			if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source] == graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal])
			{
				pair<int, double> s_min_adj = instance_graph.min_adj(source);
				float s_t_weight = float(graph_hash_of_mixed_weighted_edge_weight(instance_graph, source, terminal));
				return float(s_min_adj.second * 2) > s_t_weight ? s_t_weight : float(s_min_adj.second * 2);
			}
			else
			{
				return graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal]);
			}
		}
		/* case 4 */
		else {
			return graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], terminal);
		}	
	}
	else {
		/* case 8 -- same with case 1 */
		if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592[terminal] == 11)
		{
			return graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, source, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal]);
		}
		/* case 9 -- same with case 4 */
		else if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592[terminal] == 12)
		{
			return graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, source, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal]);
		}
		else {
			return graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_distance(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, source, terminal);
		}	
	}

}

/*query paths*/
void vector_replace_on_paths_PSL_enhancedoriginalR2(vector<pair<int, int>>& paths, int old_, int new_)
{
	/*
		this function can change specific vector on paths, used for shortesd paths query
		take care, we assume that the 'old_' vector only appear once
	*/
	for (auto it = paths.begin(); it != paths.end(); it++)
	{
		if (it->first == old_)
		{
			it->first = new_;
			break;
		}
		if (it->second == old_)
		{
			it->second = new_;
			break;
		}
	}
}

vector<pair<int, int>> graph_hash_of_mixed_weighted_HL_PSL_with_2019R1_extract_shortest_path
(vector<vector<graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_label>>& L, vector<int>& reduction_measures, 
	vector<int>& graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592, vector<int>& graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, graph_hash_of_mixed_weighted& instance_graph, int source, int terminal)
{

	vector<pair<int, int>> paths;
	if (source == terminal)
	{
		return paths;
	}
	
	if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592[source] == 11)
	{
		/* case 2 */
		if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592[terminal] == 11)
		{
			if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source] == graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal])
			{
				pair<int, double> s_min_adj = instance_graph.min_adj(source);
				paths.push_back({ source, s_min_adj.first });
				paths.push_back({ terminal, s_min_adj.first });
				return paths;
			}
			else
			{
				paths = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_shortest_path(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal]);
				vector_replace_on_paths_PSL_enhancedoriginalR2(paths, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], source);
				vector_replace_on_paths_PSL_enhancedoriginalR2(paths, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal], terminal);
				return paths;
			}
		}
		/* case 3 */
		else if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592[terminal] == 12)
		{
			paths = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_shortest_path(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal]);
			vector_replace_on_paths_PSL_enhancedoriginalR2(paths, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], source);
			vector_replace_on_paths_PSL_enhancedoriginalR2(paths, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal], terminal);
			return paths;
		}
		/* case 1 */
		else {
			paths = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_shortest_path(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], terminal);
			vector_replace_on_paths_PSL_enhancedoriginalR2(paths, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], source);
			return paths;
		}	
	}
	else if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592[source] == 12)
	{
		/* case 5 -- same with case 3 */
		if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592[terminal] == 11)
		{
			paths = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_shortest_path(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal]);
			vector_replace_on_paths_PSL_enhancedoriginalR2(paths, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], source);
			vector_replace_on_paths_PSL_enhancedoriginalR2(paths, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal], terminal);
			return paths;
		}
		/* case 6 */
		else if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592[terminal] == 12)
		{
			if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source] == graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal])
			{
				pair<int, double> s_min_adj = instance_graph.min_adj(source);
				float s_t_weight = float(graph_hash_of_mixed_weighted_edge_weight(instance_graph, source, terminal));
				if (float(s_min_adj.second * 2) < s_t_weight)
				{
					paths.push_back({ source, s_min_adj.first });
					paths.push_back({ terminal, s_min_adj.first });
					return paths;
				}
				else
				{
					paths.push_back({ source, terminal });
					return paths;
				}
			}
			else // graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source] != graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal]
			{
				paths = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_shortest_path(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal]);
				vector_replace_on_paths_PSL_enhancedoriginalR2(paths, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], source);
				vector_replace_on_paths_PSL_enhancedoriginalR2(paths, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal], terminal);
				return paths;
			}
		}
		/* case 4 */
		else {
			paths = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_shortest_path(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], terminal);
			vector_replace_on_paths_PSL_enhancedoriginalR2(paths, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[source], source);
			return paths;
		}
		
	}
	else {
		/* case 8 -- same with case 1 */
		if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592[terminal] == 11)
		{
			paths = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_shortest_path(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, source, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal]);
			vector_replace_on_paths_PSL_enhancedoriginalR2(paths, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal], terminal);
			return paths;
		}
		/* case 9 -- same with case 4 */
		else if (graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_reduction_measures_2019R1_592[terminal] == 12)
		{
			paths = graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_shortest_path(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, source, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal]);
			vector_replace_on_paths_PSL_enhancedoriginalR2(paths, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1[terminal], terminal);
			return paths;
		}
		/* case 7 */
		else {
			return graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_extract_shortest_path(L, reduction_measures, graph_hash_of_mixed_weighted_HL_PSL_enhancedoriginalR2_f_2019R1, instance_graph, source, terminal);
		}
		
	}

}