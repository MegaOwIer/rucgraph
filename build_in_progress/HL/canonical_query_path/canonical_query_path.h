#pragma once
#include <vector>
#include <tool_functions/ThreadPool.h>
#include <shared_mutex>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>




/*codes for querying distances or paths*/

/*query dis and path after 2019R1*/

double canonical_query_path_binary_operations_search_weight(vector<two_hop_label_v1>& input_vector, int key) {

	/*return std::numeric_limits<T>::max() if key is not in vector; time complexity O(log n)*/

	int left = 0, right = input_vector.size() - 1;

	while (left <= right) {
		int mid = left + ((right - left) / 2); // mid is between left and right (may be equal); 
		if (input_vector[mid].vertex == key) {
			return input_vector[mid].distance;
		}
		else if (input_vector[mid].vertex > key) {
			right = mid - 1;
		}
		else {
			left = mid + 1;
		}
	}

	return std::numeric_limits<double>::max();

}

double canonical_query_path_extract_distance_no_reduc(vector<vector<two_hop_label_v1>>& L, int source, int terminal, pair<bool, int>& common_hub_checked) {

	/*return std::numeric_limits<double>::max() is not connected*/

	if (source == terminal) {
		return 0;
	}

	if (!common_hub_checked.first) {
		double distance = std::numeric_limits<double>::max(); // if disconnected, return this large value
		auto vector1_check_pointer = L[source].begin();
		auto vector2_check_pointer = L[terminal].begin();
		auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
		while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end) {
			if (vector1_check_pointer->vertex == vector2_check_pointer->vertex) {
				double dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
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
	else { // already have the common hub, assume that labels are sorted from small to large

	}

	

}

double canonical_query_path_extract_distance_st_no_R1 /* we assume that source and terminal are not reduced by 2019R1 here*/
(vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, vector<int>& f_2019R1, int source, int terminal)
{

	/* we assume that source and terminal are not reduced by 2019R1*/

	//cout << "source=" << source << endl;
	//cout << "terminal=" << terminal << endl;
	//cout << "reduction_measures_2019R2[source] " << reduction_measures_2019R2[source] << endl;
	//cout << "reduction_measures_2019R2[terminal] " << reduction_measures_2019R2[terminal] << endl;


	if (source == terminal) {
		return 0;
	}

	vector<double> selected_distance = { std::numeric_limits<double>::max() }; // Store all path lengths to be selected; disconnected = std::numeric_limits<double>::max()
	auto s_adj_begin = adjs[source].begin();
	auto s_adj_end = adjs[source].end();
	auto t_adj_begin = adjs[terminal].begin();
	auto t_adj_end = adjs[terminal].end();
	if (reduction_measures_2019R2[source] == 2)
	{
		if (reduction_measures_2019R2[terminal] == 2)
		{
			/*"Both reduced"*/
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
				for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
				{
					if (f_2019R1[it1->first] == it1->first && f_2019R1[it2->first] == it2->first) {
						double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, it2->first);
						if (x == std::numeric_limits<double>::max() && full_two_hop_labels) { // if not full_two_hop_labels, then we may query a max dis while connected
							return x;
						}
						else {
							selected_distance.push_back(x + double(it1->second) + double(it2->second));
						}
					}
				}
			}
		}
		else
		{
			//cout << "here1" << endl;
			/*"Only source reduced"*/
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
				//cout << "it1->first: " << it1->first << " it1->second: " << it1->second << endl;
				if (f_2019R1[it1->first] == it1->first) { 
					double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, terminal);
					if (x == std::numeric_limits<double>::max() && full_two_hop_labels) {
						return x;
					}
					else {
						selected_distance.push_back(x + double(it1->second));
					}
				}
			}
		}
	}
	else
	{
		if (reduction_measures_2019R2[terminal] == 2)
		{
			//cout << "here1" << endl;
			/*"Only terminal reduced"*/
			for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
			{
				//cout << "it2->first: " << it2->first << " it2->second: " << it2->second << endl;
				if (f_2019R1[it2->first] == it2->first) { 
					double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, source, it2->first);
					//cout << "x " << x << endl;
					if (x == std::numeric_limits<double>::max() && full_two_hop_labels) {
						return x;
					}
					else {
						selected_distance.push_back(x + double(it2->second));
					}
				}
			}
		}
		else
		{
			/*"Nothing happened"*/
			selected_distance.push_back(graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, source, terminal));
		}
	}

	double dis = *min_element(selected_distance.begin(), selected_distance.end());

	return dis;
}

vector<pair<int, int>> canonical_query_path_extract_shortest_path_st_no_R1 /* we assume that source and terminal are not reduced by 2019R1*/
(vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, vector<int>& f_2019R1, graph_hash_of_mixed_weighted& instance_graph, int source, int terminal)
{

	/* we assume that source and terminal are not reduced by 2019R1*/

	vector<pair<int, int>> paths;
	if (source == terminal)
	{
		return paths;
	}

	double min_dis = std::numeric_limits<double>::max();
	vector<pair<int, int>> partial_edges(2);

	auto s_adj_begin = adjs[source].begin();
	auto s_adj_end = adjs[source].end();
	auto t_adj_begin = adjs[terminal].begin();
	auto t_adj_end = adjs[terminal].end();

	if (reduction_measures_2019R2[source] == 2)
	{
		if (reduction_measures_2019R2[terminal] == 2)
		{
			/*"Both reduced"*/
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
				for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
				{
					if (f_2019R1[it1->first] == it1->first && f_2019R1[it2->first] == it2->first) { 
						double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, it2->first) + double(it1->second) + double(it2->second);
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

			if (min_dis == std::numeric_limits<double>::max())
			{ // disconnected
				return paths;
			}

			// add partial_edges into paths
			paths.push_back(partial_edges[0]);
			paths.push_back(partial_edges[1]);

			// find new edges
			vector<pair<int, int>> new_edges = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, partial_edges[0].first, partial_edges[1].first);
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
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
				if (f_2019R1[it1->first] == it1->first) { 
					double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, terminal) + double(it1->second);
					if (min_dis > x)
					{
						min_dis = x;
						partial_edges[0] = { it1->first, source };
					}
				}
			}
			if (min_dis == std::numeric_limits<double>::max())
			{ // disconnected
				return paths;
			}

			// add partial_edges into paths
			paths.push_back(partial_edges[0]);

			// find new edges
			vector<pair<int, int>> new_edges = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, partial_edges[0].first, terminal);
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
		if (reduction_measures_2019R2[terminal] == 2)
		{
			for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
			{
				if (f_2019R1[it2->first] == it2->first) { 
					double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, source, it2->first) + double(it2->second);
					if (min_dis > x)
					{
						min_dis = x;
						partial_edges[0] = { it2->first, terminal };
					}
				}
			}
			if (min_dis == std::numeric_limits<double>::max())
			{ // disconnected
				return paths;
			}
			// add partial_edges into paths
			paths.push_back(partial_edges[0]);
			// find new edges
			vector<pair<int, int>> new_edges = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, partial_edges[0].first);
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
			double distance = std::numeric_limits<double>::max(); // if disconnected, retun this large value
			bool connected = false;
			auto vector1_check_pointer = L[source].begin();
			auto vector2_check_pointer = L[terminal].begin();
			auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
			{
				if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
				{
					connected = true;
					double dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
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
			vector<pair<int, int>> new_edges = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, terminal);

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

/*query dis (final)*/
double canonical_query_path_extract_distance
(vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, vector<int>& reduction_measures_2019R1, 
	vector<int>& f_2019R1, graph_hash_of_mixed_weighted& instance_graph, int source, int terminal)
{

	if (source == terminal) {
		return 0;
	}

	if (reduction_measures_2019R1[source] == 11)
	{
		/* case 2 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			if (f_2019R1[source] == f_2019R1[terminal])
			{
				pair<int, double> s_min_adj = min_adjs[source];
				return double(s_min_adj.second * 2);
			}
			else
			{
				return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, f_2019R1[source], f_2019R1[terminal]);
			}
		}
		/* case 3 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, f_2019R1[source], f_2019R1[terminal]);
		}
		else { /* case 1 */
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, f_2019R1[source], terminal);
		}	
	}
	else if (reduction_measures_2019R1[source] == 12)
	{
		/* case 5 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, f_2019R1[source], f_2019R1[terminal]);
		}
		/* case 6 -- same with case 3 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			if (f_2019R1[source] == f_2019R1[terminal])
			{
				pair<int, double> s_min_adj = min_adjs[source];
				double s_t_weight = double(graph_hash_of_mixed_weighted_edge_weight(instance_graph, source, terminal));
				return double(s_min_adj.second * 2) > s_t_weight ? s_t_weight : double(s_min_adj.second * 2);
			}
			else
			{
				return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, f_2019R1[source], f_2019R1[terminal]);
			}
		}
		/* case 4 */
		else {
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, f_2019R1[source], terminal);
		}	
	}
	else {
		/* case 8 -- same with case 1 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, source, f_2019R1[terminal]);
		}
		/* case 9 -- same with case 4 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			//cout << "here" << endl;
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, source, f_2019R1[terminal]);
		}
		else {
			//cout << "here" << endl;
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_st_no_R1(L, reduction_measures_2019R2, f_2019R1, source, terminal);
		}	
	}

}

/*query paths (final)*/
void canonical_query_path_vector_replace_on_paths(vector<pair<int, int>>& paths, int old_, int new_)
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

vector<pair<int, int>> canonical_query_path_extract_shortest_path
(vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, vector<int>& reduction_measures_2019R1, 
	vector<int>& f_2019R1, graph_hash_of_mixed_weighted& instance_graph, int source, int terminal)
{

	vector<pair<int, int>> paths;
	if (source == terminal)
	{
		return paths;
	}
	
	if (reduction_measures_2019R1[source] == 11)
	{
		/* case 2 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			if (f_2019R1[source] == f_2019R1[terminal])
			{
				pair<int, double> s_min_adj = min_adjs[source];
				paths.push_back({ source, s_min_adj.first });
				paths.push_back({ terminal, s_min_adj.first });
				return paths;
			}
			else
			{
				paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
				graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
				graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
				return paths;
			}
		}
		/* case 3 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
			return paths;
		}
		/* case 1 */
		else {
			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], terminal);
			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
			return paths;
		}	
	}
	else if (reduction_measures_2019R1[source] == 12)
	{
		/* case 5 -- same with case 3 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
			return paths;
		}
		/* case 6 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			if (f_2019R1[source] == f_2019R1[terminal])
			{
				pair<int, double> s_min_adj = min_adjs[source];
				double s_t_weight = double(graph_hash_of_mixed_weighted_edge_weight(instance_graph, source, terminal));
				if (double(s_min_adj.second * 2) < s_t_weight)
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
			else // f_2019R1[source] != f_2019R1[terminal]
			{
				paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
				graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
				graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
				return paths;
			}
		}
		/* case 4 */
		else {
			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], terminal);
			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[source], source);
			return paths;
		}
		
	}
	else {
		/* case 8 -- same with case 1 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, f_2019R1[terminal]);
			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
			return paths;
		}
		/* case 9 -- same with case 4 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			paths = graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, f_2019R1[terminal]);
			graph_hash_of_mixed_weighted_two_hop_v1_vector_replace_on_paths(paths, f_2019R1[terminal], terminal);
			return paths;
		}
		/* case 7 */
		else {
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_shortest_path_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, terminal);
		}
		
	}

}







/*query predecessors in instance_graph*/

pair<int, int> canonical_query_path_extract_two_predecessors_st_no_R1 /* we assume that source and terminal are not reduced by 2019R1*/
(vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, vector<int>& f_2019R1, graph_hash_of_mixed_weighted& instance_graph, int source, int terminal)
{

	/* we assume that source and terminal are not reduced by 2019R1*/
	
	if (source == terminal)
	{
		return { source, terminal }; // source_predecessor and terminal_predecessor;
	}

	double min_dis = std::numeric_limits<double>::max();
	vector<pair<int, int>> partial_edges(2);

	auto s_adj_begin = adjs[source].begin();
	auto s_adj_end = adjs[source].end();
	auto t_adj_begin = adjs[terminal].begin();
	auto t_adj_end = adjs[terminal].end();

	if (reduction_measures_2019R2[source] == 2)
	{
		if (reduction_measures_2019R2[terminal] == 2)
		{
			/*"Both reduced"*/
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
				for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
				{
					if (f_2019R1[it1->first] == it1->first && f_2019R1[it2->first] == it2->first) {
						double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, it2->first) + double(it1->second) + double(it2->second);
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

			if (min_dis == std::numeric_limits<double>::max())
			{ // disconnected
				return { source, terminal };
			}

			return { partial_edges[0].first, partial_edges[1].first };
		}
		/*"Only source reduced"*/
		else
		{
			for (auto it1 = s_adj_begin; it1 != s_adj_end; it1++)
			{
				if (f_2019R1[it1->first] == it1->first) {
					double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, it1->first, terminal) + double(it1->second);
					if (min_dis > x)
					{
						min_dis = x;
						partial_edges[0] = { it1->first, source };
					}
				}
			}
			if (min_dis == std::numeric_limits<double>::max())
			{ // disconnected
				return { source, terminal };
			}

			pair<int, int> recursive_two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, partial_edges[0].first, terminal);

			return { partial_edges[0].first , recursive_two_predecessors.second };
		}
	}
	else
	{
		/*"Only terminal reduced"*/
		if (reduction_measures_2019R2[terminal] == 2)
		{
			for (auto it2 = t_adj_begin; it2 != t_adj_end; it2++)
			{
				if (f_2019R1[it2->first] == it2->first) {
					double x = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance_no_reduc(L, source, it2->first) + double(it2->second);
					if (min_dis > x)
					{
						min_dis = x;
						partial_edges[0] = { it2->first, terminal };
					}
				}
			}

			if (min_dis == std::numeric_limits<double>::max())
			{ // disconnected
				return { source, terminal };
			}

			pair<int, int> recursive_two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, partial_edges[0].first);

			return { recursive_two_predecessors.first, partial_edges[0].first };
		}
		/*
		"Nothing happened"
		In this case, the problem that the removed vertices appear in the path needs to be solved
		*/
		else
		{
			int vector1_capped_v_parent = 0, vector2_capped_v_parent = 0;
			double distance = std::numeric_limits<double>::max(); // if disconnected, return this large value
			bool connected = false;
			auto vector1_check_pointer = L[source].begin();
			auto vector2_check_pointer = L[terminal].begin();
			auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
			while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end)
			{
				if (vector1_check_pointer->vertex == vector2_check_pointer->vertex)
				{
					connected = true;
					double dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
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
				pair<int, int> recursive_two_predecessors;
				if (source != vector1_capped_v_parent) // ascending from source
				{
					recursive_two_predecessors.first = vector1_capped_v_parent;
				}
				else {
					recursive_two_predecessors.first = source;
				}
				if (terminal != vector2_capped_v_parent) // ascending from terminal
				{
					recursive_two_predecessors.second = vector2_capped_v_parent;
				}
				else {
					recursive_two_predecessors.second = terminal;
				}
				return recursive_two_predecessors;
			}
			else
			{
				return { source, terminal };
			}
		}
	}
}

pair<int, int> canonical_query_path_extract_two_predecessors
(vector<vector<two_hop_label_v1>>& L, vector<int>& reduction_measures_2019R2, vector<int>& reduction_measures_2019R1,
	vector<int>& f_2019R1, graph_hash_of_mixed_weighted& instance_graph, int source, int terminal)
{

	pair<int, int> two_predecessors; // source_predecessor and terminal_predecessor;  
	/*
	if (source_predecessor == source && terminal_predecessor == terminal){
	no edge is in path
	}

	if (source_predecessor == terminal && terminal_predecessor == source){
	one edge is in path
	}

	*/

	if (source == terminal)
	{
		return { source, terminal };
	}

	if (reduction_measures_2019R1[source] == 11)
	{
		/* case 2 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			if (f_2019R1[source] == f_2019R1[terminal])
			{
				pair<int, double> s_min_adj = min_adjs[source];
				return { s_min_adj.first , s_min_adj.first };
			}
			else
			{
				two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
				if (two_predecessors.first == f_2019R1[source]) {
					two_predecessors.first = source;
				}
				if (two_predecessors.second == f_2019R1[source]) {
					two_predecessors.second = source;
				}
				if (two_predecessors.first == f_2019R1[terminal]) {
					two_predecessors.first = terminal;
				}
				if (two_predecessors.second == f_2019R1[terminal]) {
					two_predecessors.second = terminal;
				}
				return two_predecessors;
			}
		}
		/* case 3 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
			if (two_predecessors.first == f_2019R1[source]) {
				two_predecessors.first = source;
			}
			if (two_predecessors.second == f_2019R1[source]) {
				two_predecessors.second = source;
			}
			if (two_predecessors.first == f_2019R1[terminal]) {
				two_predecessors.first = terminal;
			}
			if (two_predecessors.second == f_2019R1[terminal]) {
				two_predecessors.second = terminal;
			}
			return two_predecessors;
		}
		/* case 1 */
		else {
			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], terminal);
			if (two_predecessors.first == f_2019R1[source]) {
				two_predecessors.first = source;
			}
			if (two_predecessors.second == f_2019R1[source]) {
				two_predecessors.second = source;
			}
			return two_predecessors;
		}
	}
	else if (reduction_measures_2019R1[source] == 12)
	{
		/* case 5 -- same with case 3 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
			if (two_predecessors.first == f_2019R1[source]) {
				two_predecessors.first = source;
			}
			if (two_predecessors.second == f_2019R1[source]) {
				two_predecessors.second = source;
			}
			if (two_predecessors.first == f_2019R1[terminal]) {
				two_predecessors.first = terminal;
			}
			if (two_predecessors.second == f_2019R1[terminal]) {
				two_predecessors.second = terminal;
			}
			return two_predecessors;
		}
		/* case 6 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			if (f_2019R1[source] == f_2019R1[terminal])
			{
				pair<int, double> s_min_adj = min_adjs[source];
				double s_t_weight = double(graph_hash_of_mixed_weighted_edge_weight(instance_graph, source, terminal));
				if (double(s_min_adj.second * 2) < s_t_weight)
				{
					return { s_min_adj.first , s_min_adj.first };
				}
				else
				{
					return { terminal, source };
				}
			}
			else // f_2019R1[source] != f_2019R1[terminal]
			{
				two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], f_2019R1[terminal]);
				if (two_predecessors.first == f_2019R1[source]) {
					two_predecessors.first = source;
				}
				if (two_predecessors.second == f_2019R1[source]) {
					two_predecessors.second = source;
				}
				if (two_predecessors.first == f_2019R1[terminal]) {
					two_predecessors.first = terminal;
				}
				if (two_predecessors.second == f_2019R1[terminal]) {
					two_predecessors.second = terminal;
				}
				return two_predecessors;
			}
		}
		/* case 4 */
		else {
			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, f_2019R1[source], terminal);
			if (two_predecessors.first == f_2019R1[source]) {
				two_predecessors.first = source;
			}
			if (two_predecessors.second == f_2019R1[source]) {
				two_predecessors.second = source;
			}
			return two_predecessors;
		}
	}
	else {
		/* case 8 -- same with case 1 */
		if (reduction_measures_2019R1[terminal] == 11)
		{
			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, f_2019R1[terminal]);
			if (two_predecessors.first == f_2019R1[terminal]) {
				two_predecessors.first = terminal;
			}
			if (two_predecessors.second == f_2019R1[terminal]) {
				two_predecessors.second = terminal;
			}
			return two_predecessors;
		}
		/* case 9 -- same with case 4 */
		else if (reduction_measures_2019R1[terminal] == 12)
		{
			two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, f_2019R1[terminal]);
			if (two_predecessors.first == f_2019R1[terminal]) {
				two_predecessors.first = terminal;
			}
			if (two_predecessors.second == f_2019R1[terminal]) {
				two_predecessors.second = terminal;
			}
			return two_predecessors;
		}
		/* case 7 */
		else {
			return graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors_st_no_R1(L, reduction_measures_2019R2, f_2019R1, instance_graph, source, terminal);
		}

	}

}




void CT_canonical_query_path(graph_hash_of_mixed_weighted_CT_v2_case_info& case_info, int source, int terminal, vector<pair<int, int>>& path) {

	pair<bool, int> common_hub_checked = { false, 0 }; // if true, then second is the common hub for SP

	/*may return INT_MAX, INT_MAX*/

	auto& root = case_info.root;
	auto& isIntree = case_info.isIntree;
	auto& Bags = case_info.Bags;
	auto& L = case_info.two_hop_case_info.L;
	auto& merging_predecessors = case_info.merging_predecessors;

	/*if source and terminal are disconnected, then return empty vector; note that, if source==terminal, then also return empty vector*/
	if (source == terminal) return;

	int vector1_capped_v_parent = -1, vector2_capped_v_parent = -1;

	double distance = std::numeric_limits<double>::max();  // record the shortest distance

	if (!isIntree[source] && !isIntree[terminal]) // both on core, use PLL indexs only
	{
		//cout << "a" << endl;

		pair<int, int> added_edge = { INT_MAX, INT_MAX };
		pair<int, int> two_predecessors = graph_hash_of_mixed_weighted_two_hop_v1_extract_two_predecessors
		(L, case_info.two_hop_case_info.reduction_measures_2019R2, case_info.two_hop_case_info.reduction_measures_2019R1, case_info.two_hop_case_info.f_2019R1, case_info.core_graph, source, terminal);
		if (two_predecessors.first == source && two_predecessors.second == terminal) { // disconnected
			path.push_back({ INT_MAX , INT_MAX });
			return;
		}

		//cout << "two_predecessors " << two_predecessors.first << " " << two_predecessors.second << endl;

		if (source != two_predecessors.first)
		{
			midnode node = search_midnode(case_info.merging_predecessors, source, two_predecessors.first);
			if (two_predecessors.first != node.left) {
				two_predecessors.first = node.left; // node.left should be the neighbour of source in the SP from source to two_predecessors.first
			}
			path.push_back({ source , two_predecessors.first });
			added_edge = { source , two_predecessors.first };
			source = two_predecessors.first;
		}
		if (terminal != two_predecessors.second)
		{
			midnode node = search_midnode(case_info.merging_predecessors, terminal, two_predecessors.second);
			if (two_predecessors.second != node.left) {
				two_predecessors.second = node.left; // node.left should be the neighbour of terminal in the SP from terminal to two_predecessors.second
			}
			if (!(added_edge.first == two_predecessors.second && added_edge.second == terminal)) {
				path.push_back({ two_predecessors.second , terminal });
				terminal = two_predecessors.second; // else, two_predecessors.second == source, terminal should not change to source, since source has changed to terminal above
			}
			else {
				return;
			}
		}

		//cout<< "resursive CT_extract_path " << source << " " << terminal << endl;

		CT_extract_path(case_info, source, terminal, path);
		return;
	}
	else if ((!isIntree[source] && isIntree[terminal]) || (isIntree[source] && !isIntree[terminal])) // source on core and terminal on tree
	{
		//cout << "b" << endl;
		// 3-hop, use the interface
		if (isIntree[source])
		{
			int t = source; source = terminal; terminal = t;
		}
		/*the following: source is in core, terminal is in tree*/
		int r = root[terminal];
		int r_size = Bags[r].size();
		int terminal_predecessor = INT_MAX;
		//cout << "source=" << source << endl;
		//cout << "terminal=" << terminal << endl;
		//cout << "r=" << r << endl;
		for (int i = 0; i < r_size; i++)
		{
			if (L[terminal][i].distance > distance)  continue; // already exceed the present minimum distance

			int x = L[terminal][i].vertex;
			double d_dis = L[terminal][i].distance;

			double dis = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance
			(L, case_info.two_hop_case_info.reduction_measures_2019R2, case_info.two_hop_case_info.reduction_measures_2019R1, case_info.two_hop_case_info.f_2019R1, case_info.core_graph, source, x);
			//cout << "x=" << x << endl;
			//cout << "d_dis=" << d_dis << endl;
			//cout << "dis=" << dis << endl;
			if (distance > dis + d_dis) {
				distance = dis + d_dis;
				terminal_predecessor = L[terminal][i].parent_vertex; // terminal is in the tree, so this is the predecessor in the original graph
			}
		}
		if (terminal_predecessor == INT_MAX) { // disconnected
			path.push_back({ INT_MAX , INT_MAX });
			return;
		}
		path.push_back({ terminal_predecessor, terminal });
		terminal = terminal_predecessor;
		CT_extract_path(case_info, source, terminal, path);

		return;
	}
	else if (root[source] != root[terminal])  // locate in different trees
	{
		//cout << "c" << endl;
		int r_s = root[source];
		int r_s_size = Bags[r_s].size();
		int r_t = root[terminal];
		int r_t_size = Bags[r_t].size();
		int source_predecessor = INT_MAX, terminal_predecessor = INT_MAX;
		for (int i = 0; i < r_s_size; i++) {
			int u = L[source][i].vertex;
			double d_s_u = L[source][i].distance;
			for (int j = 0; j < r_t_size; j++) {
				int w = L[terminal][j].vertex;
				double d_t_w = L[terminal][j].distance;
				double d_u_w = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance
				(L, case_info.two_hop_case_info.reduction_measures_2019R2, case_info.two_hop_case_info.reduction_measures_2019R1, case_info.two_hop_case_info.f_2019R1, case_info.core_graph, u, w);
				double dis = d_s_u + d_t_w + d_u_w;
				if (dis < distance) {
					distance = dis;
					source_predecessor = L[source][i].parent_vertex;
					terminal_predecessor = L[terminal][j].parent_vertex;
				}
			}
		}
		if (source_predecessor == INT_MAX || terminal_predecessor == INT_MAX) { // disconnected
			path.push_back({ INT_MAX , INT_MAX });
			return;
		}
		pair<int, int> added_edge = { INT_MAX, INT_MAX };
		if (source != source_predecessor)
		{
			path.push_back({ source , source_predecessor });
			added_edge = { source , source_predecessor };
			source = source_predecessor;
		}
		if (terminal != terminal_predecessor)
		{
			if (!(added_edge.first == terminal_predecessor && added_edge.second == terminal)) {
				path.push_back({ terminal_predecessor , terminal });
				terminal = terminal_predecessor;
			}
			else {
				return;
			}
		}
		CT_extract_path(case_info, source, terminal, path);

		return;
	}
	else
	{
		//cout << "d" << endl;
		int source_predecessor = INT_MAX, terminal_predecessor = INT_MAX;

		int grand = lca(case_info, source, terminal);

		/*compute d2*/
		unordered_set<int> Bk = { grand };
		for (int i = Bags[grand].size() - 1; i >= 0; i--) {
			Bk.insert(Bags[grand][i].first);
		}
		unordered_map<int, pair<double, int>> Bk_source_dis;
		for (int i = L[source].size() - 1; i >= 0; i--) {
			int v = L[source][i].vertex;
			if (Bk.count(v) > 0) {
				Bk_source_dis[v] = { L[source][i].distance , L[source][i].parent_vertex };
			}
		}
		for (int i = L[terminal].size() - 1; i >= 0; i--) {
			int v = L[terminal][i].vertex;
			if (Bk_source_dis.count(v) > 0) {
				double d2 = Bk_source_dis[v].first + L[terminal][i].distance;
				if (distance > d2) {
					distance = d2;
					source_predecessor = Bk_source_dis[v].second;
					terminal_predecessor = L[terminal][i].parent_vertex;
				}
			}
		}

		/*compute d4*/
		int r = root[source]; // source and terminal have the same root
		int r_size = Bags[r].size();
		for (int i = 0; i < r_size; i++)
		{
			int u = L[source][i].vertex;
			double d_s_u = L[source][i].distance;
			for (int j = 0; j < r_size; j++)
			{
				int w = L[terminal][j].vertex;
				double d_t_w = L[terminal][j].distance;
				double d_u_w = graph_hash_of_mixed_weighted_two_hop_v1_extract_distance
				(L, case_info.two_hop_case_info.reduction_measures_2019R2, case_info.two_hop_case_info.reduction_measures_2019R1, case_info.two_hop_case_info.f_2019R1, case_info.core_graph, u, w);
				double d4 = d_s_u + d_u_w + d_t_w;
				if (distance > d4) {
					distance = d4;
					source_predecessor = L[source][i].parent_vertex;
					terminal_predecessor = L[terminal][j].parent_vertex;
				}
			}
		}

		if (source_predecessor == INT_MAX || terminal_predecessor == INT_MAX) { // disconnected
			path.push_back({ INT_MAX , INT_MAX });
			return;
		}
		pair<int, int> added_edge = { INT_MAX, INT_MAX };
		if (source != source_predecessor)
		{
			path.push_back({ source , source_predecessor });
			added_edge = { source , source_predecessor };
			source = source_predecessor;
		}
		if (terminal != terminal_predecessor)
		{
			if (!(added_edge.first == terminal_predecessor && added_edge.second == terminal)) {
				path.push_back({ terminal_predecessor , terminal });
				terminal = terminal_predecessor;
			}
			else {
				return;
			}
		}
		CT_extract_path(case_info, source, terminal, path);


		return;

	}

}