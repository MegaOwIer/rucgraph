#pragma once

/*PLL;

the follow codes are much slower than PLL in 2013 Japan paper;
the major reason is that the 2013 Japan paper does not consider edge costs, and use BFS O(|E|+|V|), and their vertex IDs are from 0 to V-1?,
while the follow codes consider edge costs, and use Dijkstra's algorithm O(|E|+|V|log|V|), and the follow vertex IDs can be any int value*/

bool graph_v_of_v_idealID_PLL_Creating_Indexes_compare(const pair<int, int>& i, pair<int, int>& j)
{
	return i.second > j.second;  // < is from small to big; > is from big to small.  sort by the second item of pair<int, int>
}

struct graph_v_of_v_idealID_PLL_node_for_sp {
	int vertex, parent_vertex;
	double priority_value;
}; // define the node in the queue
bool operator<(graph_v_of_v_idealID_PLL_node_for_sp const& x, graph_v_of_v_idealID_PLL_node_for_sp const& y) {
	return x.priority_value > y.priority_value; // < is the max-heap; > is the min heap
}
typedef typename boost::heap::fibonacci_heap<graph_v_of_v_idealID_PLL_node_for_sp>::handle_type graph_v_of_v_idealID_PLL_handle_t_for_sp;


#include<graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_PLL_labels.h>
#include <graph_v_of_v_idealID/graph_v_of_v_idealID_change_new_vertexIDs.h>

vector<vector<PLL_sorted_label>> graph_v_of_v_idealID_PLL_generate_indexes(graph_v_of_v_idealID& input_graph) {

	/*this function is different from (and faster than) PLL_hash_of_vectors_generate_indexes_basic
	in that this function change vertexIDs to ideal IDs,
	and then use the rank and T P tricks in 2013 Japan paper*/

	int N = input_graph.size();

	/*sort vertices by degrees*/
	vector<pair<int, int>> sorted_vertices;
	for (int i = 0; i < N; i++) {
		sorted_vertices.push_back({ i, input_graph[i].size() });
	}
	sort(sorted_vertices.begin(), sorted_vertices.end(), graph_v_of_v_idealID_PLL_Creating_Indexes_compare);

	vector<int> vertexIDs_new_to_old(N);
	vector<int> vertexIDs_old_to_new(N);
	for (int i = 0; i < N; i++) {
		vertexIDs_new_to_old[i] = sorted_vertices[i].first;
		vertexIDs_old_to_new[sorted_vertices[i].first] = i;
	}
	graph_v_of_v_idealID ideal_graph = graph_v_of_v_idealID_change_new_vertexIDs(input_graph, vertexIDs_old_to_new);



	/*the following variables use rank IDs*/
	vector<vector<PLL_sorted_label>> L(N); // Labels for rank IDs
	vector<double> P(N, std::numeric_limits<double>::max()); // P is distance to v_k, and priorities in Q
	vector<double> T(N, std::numeric_limits<double>::max());
	vector<graph_v_of_v_idealID_PLL_handle_t_for_sp> Q_handles(N);

	/*PLL: Pruned Landmark Labeling*/
	for (int v_k = 0; v_k < N; v_k++) {

		/*Pruned Dijkstra from vertex v_k; see Algorithm 1 in 2013 Japan SIGMOD paper*/

		queue<int> P_changed_vertices;

		graph_v_of_v_idealID_PLL_node_for_sp node;
		boost::heap::fibonacci_heap<graph_v_of_v_idealID_PLL_node_for_sp> Q;
		PLL_sorted_label xx;

		node.vertex = v_k;
		node.parent_vertex = v_k;
		node.priority_value = 0;
		Q_handles[v_k] = Q.push(node);
		P[v_k] = 0;
		P_changed_vertices.push(v_k);


		int L_v_k_size = L[v_k].size();
		for (int i = 0; i < L_v_k_size; i++) {
			T[L[v_k][i].vertex] = L[v_k][i].distance; //allocate T values for L[v_k]
		}


		while (Q.size() > 0) {

			node = Q.top();
			Q.pop();
			int u = node.vertex;
			int u_parent = node.parent_vertex;
			double P_u = node.priority_value;
			double P_u_with_error = P_u + 1e-5; // 1e-5 is error

			/*merge-joinly compute query_v_k_u, which is the QUERY value in Line 7 of Algo. 1 in 2013 Japan paper*/
			double query_v_k_u = std::numeric_limits<double>::max(); // if disconnected, retun this large value
			auto L_u_size = L[u].size(); // a vector<PLL_sorted_label>
			for (int i = 0; i < L_u_size; i++) {
				double dis = L[u][i].distance + T[L[u][i].vertex];
				if (query_v_k_u > dis) {
					query_v_k_u = dis;
				}
			}

			if (P_u_with_error < query_v_k_u) { // this is pruning

				xx.vertex = v_k;
				xx.distance = P_u;
				xx.parent_vertex = u_parent;
				L[u].push_back(xx);

				int u_adj_size = ideal_graph[u].size();
				for (int i = 0; i < u_adj_size; i++) {
					int adj_v = ideal_graph[u][i].first;
					double ec = ideal_graph[u][i].second;
					if (P[adj_v] == std::numeric_limits<double>::max()) {
						node.vertex = adj_v;
						node.parent_vertex = u;
						node.priority_value = P_u + ec;
						Q_handles[adj_v] = Q.push(node);
						P[adj_v] = node.priority_value;
						P_changed_vertices.push(adj_v);
					}
					else {
						if (P[adj_v] > P_u + ec) {
							node.vertex = adj_v;
							node.parent_vertex = u;
							node.priority_value = P_u + ec;
							Q.update(Q_handles[adj_v], node);
							P[adj_v] = node.priority_value;
						}
					}
				}
			}
		}


		for (int i = 0; i < L_v_k_size; i++) {
			T[L[v_k][i].vertex] = std::numeric_limits<double>::max(); // reverse-allocate T values for L[v_k]
		}

		while (P_changed_vertices.size() > 0) {
			P[P_changed_vertices.front()] = std::numeric_limits<double>::max(); // reverse-allocate P values
			P_changed_vertices.pop();
		}

	}


	/*return unordered_map_L for old_IDs*/
	vector<vector<PLL_sorted_label>> L_for_oldIDs(N);
	for (int v_k = 0; v_k < N; v_k++) {
		vector<PLL_sorted_label> xx_vector = {};
		int L_v_k_size = L[v_k].size();
		for (int i = 0; i < L_v_k_size; i++) {
			PLL_sorted_label_binary_insert(xx_vector, vertexIDs_new_to_old[L[v_k][i].vertex], vertexIDs_new_to_old[L[v_k][i].parent_vertex], L[v_k][i].distance);
		}
		L_for_oldIDs[vertexIDs_new_to_old[v_k]] = xx_vector;
	}
	return L_for_oldIDs;

}

void graph_v_of_v_idealID_PLL_generate_and_save_indexes(graph_v_of_v_idealID& input_graph, string index_file_name) {

	auto begin1 = std::chrono::high_resolution_clock::now();

	vector<vector<PLL_sorted_label>> L = graph_v_of_v_idealID_PLL_generate_indexes(input_graph);

	auto end1 = std::chrono::high_resolution_clock::now();
	double runningtime1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1).count() / 1e9; // s


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
	outputFile.open(index_file_name);
	outputFile << "SECTION Comments" << std::endl;
	outputFile << "Creator Yahui SUN ËïÑÇ»Ô" << std::endl;
	outputFile << "Graph_Size: |V|= " << input_graph.size() << " |E|= " << graph_v_of_v_idealID_total_edge_num(input_graph) << std::endl;
	outputFile << "Total Number of Indexes: " << index_size << " Consumed Time: " << runningtime1 << "s"
		<< " Consumed RAM: " << Current_Memory_Consumption_of_This_Process() << "MB" << std::endl;
	outputFile << "Number of Indexes Per Vertex: " << (double)index_size / input_graph.size() << std::endl;
	outputFile << "Comments END at Line 6" << std::endl;
	int L_size = L.size();
	for (int i = 0; i < L_size; i++) {
		outputFile << i << ":"; // v_1:
		for (auto it2 = L[i].begin(); it2 != L[i].end(); it2++) {
			outputFile << "<" << it2->vertex << "," << it2->distance << "," << it2->parent_vertex << ">"; // <v_x,dis,pre>
		}
		outputFile << std::endl;
	}

	auto end2 = std::chrono::high_resolution_clock::now();
	double runningtime2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2).count() / 1e9; // s

	/*runningtime2 is around 1/5 runningtime1*/
	cout << index_file_name + " PLL runningtime1: " << runningtime1 << "s" << std::endl;
	cout << index_file_name + " Saving runningtime2: " << runningtime2 << "s" << std::endl;

}

vector<vector<PLL_sorted_label>> graph_v_of_v_idealID_PLL_generate_indexes_one_edge_weight(graph_v_of_v_idealID& input_graph) {

	/*this function is different from (and faster than) PLL_hash_of_vectors_generate_indexes_basic
	in that this function change vertexIDs to ideal IDs,
	and then use the rank and T P tricks in 2013 Japan paper*/

	int N = input_graph.size();

	/*sort vertices by degrees*/
	vector<pair<int, int>> sorted_vertices;
	for (int i = 0; i < N; i++) {
		sorted_vertices.push_back({ i, input_graph[i].size() });
	}
	sort(sorted_vertices.begin(), sorted_vertices.end(), graph_v_of_v_idealID_PLL_Creating_Indexes_compare);

	vector<int> vertexIDs_new_to_old(N);
	vector<int> vertexIDs_old_to_new(N);
	for (int i = 0; i < N; i++) {
		vertexIDs_new_to_old[i] = sorted_vertices[i].first;
		vertexIDs_old_to_new[sorted_vertices[i].first] = i;
	}
	graph_v_of_v_idealID ideal_graph = graph_v_of_v_idealID_change_new_vertexIDs(input_graph, vertexIDs_old_to_new);



	/*the following variables use rank IDs*/
	vector<vector<PLL_sorted_label>> L(N); // Labels for rank IDs
	vector<int> P(N, INT_MAX); // P is distance to v_k, and priorities in Q
	vector<int> T(N, INT_MAX);

	/*PLL: Pruned Landmark Labeling*/
	for (int v_k = 0; v_k < N; v_k++) {

		/*Pruned Dijkstra from vertex v_k; see Algorithm 1 in 2013 Japan SIGMOD paper*/

		queue<int> P_changed_vertices;
		queue<graph_v_of_v_idealID_PLL_node_for_sp> Q;

		graph_v_of_v_idealID_PLL_node_for_sp node;
		PLL_sorted_label xx;

		node.vertex = v_k;
		node.parent_vertex = v_k;
		node.priority_value = 0;
		Q.push(node);
		P[v_k] = 0;
		P_changed_vertices.push(v_k);


		int L_v_k_size = L[v_k].size();
		for (int i = 0; i < L_v_k_size; i++) {
			T[L[v_k][i].vertex] = L[v_k][i].distance; //allocate T values for L[v_k]
		}


		while (Q.size() > 0) {

			node = Q.front();
			Q.pop();
			int u = node.vertex;
			int u_parent = node.parent_vertex;
			int P_u = node.priority_value;

			/*merge-joinly compute query_v_k_u, which is the QUERY value in Line 7 of Algo. 1 in 2013 Japan paper*/
			double query_v_k_u = std::numeric_limits<double>::max(); // if disconnected, retun this large value
			auto L_u_size = L[u].size(); // a vector<PLL_sorted_label>
			for (int i = 0; i < L_u_size; i++) {
				double dis = L[u][i].distance + T[L[u][i].vertex];
				if (query_v_k_u > dis) {
					query_v_k_u = dis;
				}
			}

			if (P_u < query_v_k_u) { // this is pruning

				xx.vertex = v_k;
				xx.distance = P_u;
				xx.parent_vertex = u_parent;
				L[u].push_back(xx);

				int u_adj_size = ideal_graph[u].size();
				for (int i = 0; i < u_adj_size; i++) {
					int adj_v = ideal_graph[u][i].first;
					if (P[adj_v] == INT_MAX) {
						node.vertex = adj_v;
						node.parent_vertex = u;
						node.priority_value = P_u + 1;
						Q.push(node);
						P[adj_v] = node.priority_value;
						P_changed_vertices.push(adj_v);
					}
				}
			}
		}


		for (int i = 0; i < L_v_k_size; i++) {
			T[L[v_k][i].vertex] = INT_MAX; // reverse-allocate T values for L[v_k]
		}

		while (P_changed_vertices.size() > 0) {
			P[P_changed_vertices.front()] = INT_MAX; // reverse-allocate P values
			P_changed_vertices.pop();
		}

	}


	/*return unordered_map_L for old_IDs*/
	vector<vector<PLL_sorted_label>> L_for_oldIDs(N);
	for (int v_k = 0; v_k < N; v_k++) {
		vector<PLL_sorted_label> xx_vector = {};
		int L_v_k_size = L[v_k].size();
		for (int i = 0; i < L_v_k_size; i++) {
			PLL_sorted_label_binary_insert(xx_vector, vertexIDs_new_to_old[L[v_k][i].vertex], vertexIDs_new_to_old[L[v_k][i].parent_vertex], L[v_k][i].distance);
		}
		L_for_oldIDs[vertexIDs_new_to_old[v_k]] = xx_vector;
	}
	return L_for_oldIDs;

}

void graph_v_of_v_idealID_PLL_generate_and_save_indexes_one_edge_weight(graph_v_of_v_idealID& input_graph, string index_file_name) {

	auto begin1 = std::chrono::high_resolution_clock::now();

	vector<vector<PLL_sorted_label>> L = graph_v_of_v_idealID_PLL_generate_indexes_one_edge_weight(input_graph);

	auto end1 = std::chrono::high_resolution_clock::now();
	double runningtime1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1).count() / 1e9; // s


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
	outputFile.open(index_file_name);
	outputFile << "SECTION Comments" << std::endl;
	outputFile << "Creator Yahui SUN ËïÑÇ»Ô" << std::endl;
	outputFile << "Graph_Size: |V|= " << input_graph.size() << " |E|= " << graph_v_of_v_idealID_total_edge_num(input_graph) << std::endl;
	outputFile << "Total Number of Indexes: " << index_size << " Consumed Time: " << runningtime1 << "s"
		<< " Consumed RAM: " << Current_Memory_Consumption_of_This_Process() << "MB" << std::endl;
	outputFile << "Number of Indexes Per Vertex: " << (double)index_size / input_graph.size() << std::endl;
	outputFile << "Comments END at Line 6" << std::endl;
	int L_size = L.size();
	for (int i = 0; i < L_size; i++) {
		outputFile << i << ":"; // v_1:
		for (auto it2 = L[i].begin(); it2 != L[i].end(); it2++) {
			outputFile << "<" << it2->vertex << "," << it2->distance << "," << it2->parent_vertex << ">"; // <v_x,dis,pre>
		}
		outputFile << std::endl;
	}

	auto end2 = std::chrono::high_resolution_clock::now();
	double runningtime2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2).count() / 1e9; // s

	/*runningtime2 is around 1/5 runningtime1*/
	cout << index_file_name + " PLL runningtime1: " << runningtime1 << "s" << std::endl;
	cout << index_file_name + " Saving runningtime2: " << runningtime2 << "s" << std::endl;

}




#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_save_for_GSTP.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_read_for_GSTP.h>

void test_graph_v_of_v_idealID_PLL() {

	/*parameters*/
	int iteration_times = 1;
	int V = 5000, E = 20000, precision = 3;
	double ec_min = 1, ec_max = 10;

	/*iteration*/
	std::time_t now = std::time(0);
	boost::random::mt19937 gen{ static_cast<std::uint32_t>(now) };
	for (int i = 0; i < iteration_times; i++) {

		/*input and output*/
		int generate_new_graph = 1;

		graph_hash_of_mixed_weighted old_hash_graph, old_hash_generated_group_graph;;
		graph_v_of_v_idealID instance_graph;
		if (generate_new_graph == 1) {
			instance_graph = graph_v_of_v_idealID_generate_random_connected_graph(V, E, ec_min, ec_max, precision, boost_random_time_seed);
			std::unordered_set<int> generated_group_vertices;
			graph_v_of_v_idealID generated_group_graph;
			graph_v_of_v_idealID_save_for_GSTP("simple_iterative_tests.txt", instance_graph, generated_group_graph, generated_group_vertices);
			double lambda;
			graph_hash_of_mixed_weighted_read_for_GSTP("simple_iterative_tests.txt", old_hash_graph,
				old_hash_generated_group_graph, generated_group_vertices, lambda);
		}
		else {
			std::unordered_set<int> generated_group_vertices;
			graph_v_of_v_idealID generated_group_graph;
			graph_v_of_v_idealID_read_for_GSTP("simple_iterative_tests.txt", instance_graph, generated_group_graph, generated_group_vertices);
			double lambda;
			graph_hash_of_mixed_weighted_read_for_GSTP("simple_iterative_tests.txt", old_hash_graph,
				old_hash_generated_group_graph, generated_group_vertices, lambda);
		}

		graph_v_of_v_idealID_PLL_generate_and_save_indexes(instance_graph, "simple_iterative_tests_PLL_Indexes.txt");

		vector<vector<PLL_sorted_label>> L = PLL_vector_of_vectors_read_indexes("simple_iterative_tests_PLL_Indexes.txt", V+1);

		boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(V - 1) };
		int source = dist(gen);

		std::unordered_map<int, double> distances;
		std::unordered_map<int, int> predecessors;
		graph_hash_of_mixed_weighted_shortest_paths_source_to_all(old_hash_graph, source, distances, predecessors);

		for (int xx = 0; xx < 10; xx++) {
			int terminal = dist(gen);
			double dis = PLL_extract_distance_vectorFORMAT(L, source, terminal);
			if (abs(dis - distances[terminal]) > 1e-5) {
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
			vector<pair<int, int>> path = PLL_extract_shortest_path_vectorFORMAT(L, source, terminal);
			double path_dis = 0;
			for (auto it = path.begin(); it != path.end(); it++) {
				path_dis = path_dis + graph_v_of_v_idealID_edge_weight(instance_graph, it->first, it->second);
			}
			if (abs(dis - path_dis) > 1e-5) {
				cout << "source = " << source << endl;
				cout << "terminal = " << terminal << endl;
				print_vector_pair_int(path);
				cout << "dis = " << dis << endl;
				cout << "path_dis = " << path_dis << endl;
				cout << "abs(dis - path_dis) > 1e-5!" << endl;
				getchar();
			}

		}
	}
}

void test_graph_v_of_v_idealID_PLL_one_edge_weight() {

	/*parameters*/
	int iteration_times = 100;
	int V = 1000, E = 5000, precision = 1;
	double ec_min = 1, ec_max = 1;

	/*iteration*/
	std::time_t now = std::time(0);
	boost::random::mt19937 gen{ static_cast<std::uint32_t>(now) };
	for (int i = 0; i < iteration_times; i++) {

		/*input and output*/
		int generate_new_graph = 1;

		graph_hash_of_mixed_weighted old_hash_graph, old_hash_generated_group_graph;;
		graph_v_of_v_idealID instance_graph;
		if (generate_new_graph == 1) {
			instance_graph = graph_v_of_v_idealID_generate_random_connected_graph(V, E, ec_min, ec_max, precision, boost_random_time_seed);
			std::unordered_set<int> generated_group_vertices;
			graph_v_of_v_idealID generated_group_graph;
			graph_v_of_v_idealID_save_for_GSTP("simple_iterative_tests.txt", instance_graph, generated_group_graph, generated_group_vertices);
			double lambda;
			graph_hash_of_mixed_weighted_read_for_GSTP("simple_iterative_tests.txt", old_hash_graph,
				old_hash_generated_group_graph, generated_group_vertices, lambda);
		}
		else {
			std::unordered_set<int> generated_group_vertices;
			graph_v_of_v_idealID generated_group_graph;
			graph_v_of_v_idealID_read_for_GSTP("simple_iterative_tests.txt", instance_graph, generated_group_graph, generated_group_vertices);
			double lambda;
			graph_hash_of_mixed_weighted_read_for_GSTP("simple_iterative_tests.txt", old_hash_graph,
				old_hash_generated_group_graph, generated_group_vertices, lambda);
		}

		graph_v_of_v_idealID_PLL_generate_and_save_indexes_one_edge_weight(instance_graph, "simple_iterative_tests_PLL_Indexes.txt");

		vector<vector<PLL_sorted_label>> L = PLL_vector_of_vectors_read_indexes("simple_iterative_tests_PLL_Indexes.txt", V + 1);

		boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(V - 1) };
		int source = dist(gen);

		std::unordered_map<int, double> distances;
		std::unordered_map<int, int> predecessors;
		graph_hash_of_mixed_weighted_shortest_paths_source_to_all(old_hash_graph, source, distances, predecessors);

		for (int xx = 0; xx < 10; xx++) {
			int terminal = dist(gen);
			double dis = PLL_extract_distance_vectorFORMAT(L, source, terminal);
			if (abs(dis - distances[terminal]) > 1e-5) {
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
			vector<pair<int, int>> path = PLL_extract_shortest_path_vectorFORMAT(L, source, terminal);
			double path_dis = 0;
			for (auto it = path.begin(); it != path.end(); it++) {
				path_dis = path_dis + graph_v_of_v_idealID_edge_weight(instance_graph, it->first, it->second);
			}
			if (abs(dis - path_dis) > 1e-5) {
				cout << "source = " << source << endl;
				cout << "terminal = " << terminal << endl;
				print_vector_pair_int(path);
				cout << "dis = " << dis << endl;
				cout << "path_dis = " << path_dis << endl;
				cout << "abs(dis - path_dis) > 1e-5!" << endl;
				getchar();
			}

		}
	}
}





#include <Current_Memory_Consumption_of_This_Process.h>

void for_yang_shuang_20211009() {

	vector<vector<PLL_sorted_label>> L = PLL_vector_of_vectors_read_indexes("simple_iterative_tests_PLL_Indexes.txt", 5001);

	double ram = Current_Memory_Consumption_of_This_Process();
	cout << "ram = " << ram << " MB" << endl;


	int total = 0;
	for (int i = 0; i < L.size(); i++) {
		total = total + L[i].size();
	}
	cout << total << endl;


	L.clear();

	ram = Current_Memory_Consumption_of_This_Process();
	cout << "ram = " << ram << " MB" << endl;


}