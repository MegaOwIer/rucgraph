#pragma once

/*PLL;

the follow codes are much slower than PLL in 2013 Japan paper;
the major reason is that the 2013 Japan paper does not consider edge costs, and use BFS O(|E|+|V|),
while the follow codes consider edge costs, and use Dijkstra's algorithm O(|E|+|V|log|V|).

Moreover, PLL are faster when all edge costs are equal, and have smaller indexes. For example,
in the following "void test_PLL_vectorFORMAT_time()", when 
    int iteration_times = 100;
	int V = 500, E = 1000, precision = 1;
	double ec_min = 1, ec_max = 1;
then
    time_avg = 0.00989188s
    L_size_avg = 15477.8
when 
    int iteration_times = 100;
	int V = 500, E = 1000, precision = 1;
	double ec_min = 1, ec_max = 10;
then
    time_avg = 0.01335s
    L_size_avg = 18968.3


the following PLL cannot be paralleled, due to that mutiple threads writing a global L simultaneously is not thread-safe;
however, you can process each connected component parallelly.

*/

#include <Current_Memory_Consumption_of_This_Process.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_shortest_paths.h>
#include <boost/heap/fibonacci_heap.hpp> 

bool PLL_Creating_Indexes_compare(const pair<int, int>& i, pair<int, int>& j)
{
	/*< is nearly 10 times slower than >*/
	return i.second > j.second;  // < is from small to big; > is from big to small.  sort by the second item of pair<int, int>
}

struct PLL_node_for_sp {
	int vertex, parent_vertex;
	double priority_value; // double float do not change the PLL running time
}; // define the node in the queue
bool operator<(PLL_node_for_sp const& x, PLL_node_for_sp const& y) {
	return x.priority_value > y.priority_value; // < is the max-heap; > is the min heap
}
typedef typename boost::heap::fibonacci_heap<PLL_node_for_sp>::handle_type PLL_handle_t_for_sp;


#pragma region
void graph_hash_of_mixed_weighted_PLL_Create_Indexes(graph_hash_of_mixed_weighted& input_graph, string index_file_name) {

	/*indexes; L(v_1) = { v_x, <distance to v_x, predecessor of v_x in the SP between v_1 and v_x>}*/
	unordered_map<int, unordered_map<int, pair<double, int>>> L;
	long long int index_size = 0;

	auto begin = std::chrono::high_resolution_clock::now();

	/*sort vertices by degrees*/
	vector<pair<int, int>> sorted_vertices;
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		sorted_vertices.push_back({ it->first, input_graph.degree(it->first) });
		L[it->first] = {};
	}
	sort(sorted_vertices.begin(), sorted_vertices.end(), PLL_Creating_Indexes_compare);



	/*PLL: Pruned Landmark Labeling*/
	for (auto it = sorted_vertices.begin(); it != sorted_vertices.end(); it++) {

		/*Pruned Dijkstra from vertex v_k = it->first; see Algorithm 1 in 2013 Japan SIGMOD paper*/
		int v_k = it->first;

		PLL_node_for_sp node;
		boost::heap::fibonacci_heap<PLL_node_for_sp> Q;
		std::unordered_map<int, double> P; // distance to v_k; priorities in Q
		std::unordered_map<int, PLL_handle_t_for_sp> Q_handles;

		node.vertex = v_k;
		node.parent_vertex = v_k;
		node.priority_value = 0;
		Q_handles[v_k] = Q.push(node);
		P[v_k] = 0;


		while (Q.size() > 0) {

			node = Q.top();
			Q.pop();
			int u = node.vertex;
			int u_parent = node.parent_vertex;
			double P_u = node.priority_value;
			double P_u_with_error = P_u + 1e-5; // 1e-5 is error

			double query_v_k_u = std::numeric_limits<double>::max(); // if disconnected, retun this large value
			auto mm = L.find(u);
			auto nn = L.find(v_k);
			for (auto xx = mm->second.begin(); xx != mm->second.end(); xx++) {
				int v = xx->first;
				if (nn->second.count(v) > 0) {
					double dis = xx->second.first + nn->second[v].first;
					if (query_v_k_u > dis) {
						query_v_k_u = dis;
					}
				}
			}

			if (P_u_with_error < query_v_k_u) {

				L[u][v_k] = { P_u, u_parent };
				index_size++;

				auto search = input_graph.hash_of_hashs.find(u);
				if (search != input_graph.hash_of_hashs.end()) {
					for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
						int adj_v = it2->first;
						double ec = it2->second;
						if (P.count(adj_v) == 0) {
							node.vertex = adj_v;
							node.parent_vertex = u;
							node.priority_value = P_u + ec;
							Q_handles[adj_v] = Q.push(node);
							P[adj_v] = node.priority_value;
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
				else {
					auto search2 = input_graph.hash_of_vectors.find(u);
					for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
						int adj_v = it2->first;
						double ec = it2->second;
						if (P.count(adj_v) == 0) {
							node.vertex = adj_v;
							node.parent_vertex = u;
							node.priority_value = P_u + ec;
							Q_handles[adj_v] = Q.push(node);
							P[adj_v] = node.priority_value;
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
		}
	}


	auto end = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	/*save indexes*/
	std::ofstream outputFile;
	outputFile.precision(6);
	outputFile.setf(std::ios::fixed);
	outputFile.setf(std::ios::showpoint);
	outputFile.open(index_file_name);
	outputFile << "SECTION Comments" << std::endl;
	outputFile << "Creator Yahui SUN ËïÑÇ»Ô" << std::endl;
	outputFile << "Graph Size: |V|=" << input_graph.hash_of_vectors.size() << " |E|=" << graph_hash_of_mixed_weighted_num_edges(input_graph) << std::endl;
	outputFile << "Total Number of Indexes: " << index_size << " Consumed Time: " << runningtime << "s"
		<< " Consumed RAM: " << Current_Memory_Consumption_of_This_Process() << "MB" << std::endl;
	outputFile << "Number of Indexes Per Vertex: " << (double)index_size / input_graph.hash_of_vectors.size() << std::endl;
	outputFile << "Comments END at Line 6" << std::endl;
	for (auto it = L.begin(); it != L.end(); it++) {
		outputFile << it->first << ":"; // v_1:
		for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
			outputFile << "<" << it2->first << "," << it2->second.first << "," << it2->second.second << ">"; // <v_x,dis,pre>
		}
		outputFile << std::endl;
	}

}


unordered_map<int, unordered_map<int, pair<double, int>>> graph_hash_of_mixed_weighted_PLL_Generate_Indexes_no_save(graph_hash_of_mixed_weighted& input_graph) {

	/*indexes; L(v_1) = { v_x, <distance to v_x, predecessor of v_x in the SP between v_1 and v_x>}*/
	unordered_map<int, unordered_map<int, pair<double, int>>> L;
	long long int index_size = 0;

	auto begin = std::chrono::high_resolution_clock::now();

	/*sort vertices by degrees*/
	vector<pair<int, int>> sorted_vertices;
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		sorted_vertices.push_back({ it->first, input_graph.degree(it->first) });
		L[it->first] = {};
	}
	sort(sorted_vertices.begin(), sorted_vertices.end(), PLL_Creating_Indexes_compare);



	/*PLL: Pruned Landmark Labeling*/
	for (auto it = sorted_vertices.begin(); it != sorted_vertices.end(); it++) {

		/*Pruned Dijkstra from vertex v_k = it->first; see Algorithm 1 in 2013 Japan SIGMOD paper*/
		int v_k = it->first;

		PLL_node_for_sp node;
		boost::heap::fibonacci_heap<PLL_node_for_sp> Q;
		std::unordered_map<int, double> P; // distance to v_k; priorities in Q
		std::unordered_map<int, PLL_handle_t_for_sp> Q_handles;

		node.vertex = v_k;
		node.parent_vertex = v_k;
		node.priority_value = 0;
		Q_handles[v_k] = Q.push(node);
		P[v_k] = 0;


		while (Q.size() > 0) {

			node = Q.top();
			Q.pop();
			int u = node.vertex;
			int u_parent = node.parent_vertex;
			double P_u = node.priority_value;


			double query_v_k_u = std::numeric_limits<double>::max(); // if disconnected, retun this large value
			auto mm = L.find(u);
			auto nn = L.find(v_k);
			for (auto xx = mm->second.begin(); xx != mm->second.end(); xx++) {
				int v = xx->first;
				if (nn->second.count(v) > 0) {
					double dis = xx->second.first + nn->second[v].first;
					if (query_v_k_u > dis) {
						query_v_k_u = dis;
					}
				}
			}

			if (P_u < query_v_k_u) {

				L[u][v_k] = { P_u, u_parent };
				index_size++;

				auto search = input_graph.hash_of_hashs.find(u);
				if (search != input_graph.hash_of_hashs.end()) {
					for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
						int adj_v = it2->first;
						double ec = it2->second;
						if (P.count(adj_v) == 0) {
							node.vertex = adj_v;
							node.parent_vertex = u;
							node.priority_value = P_u + ec;
							Q_handles[adj_v] = Q.push(node);
							P[adj_v] = node.priority_value;
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
				else {
					auto search2 = input_graph.hash_of_vectors.find(u);
					for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
						int adj_v = it2->first;
						double ec = it2->second;
						if (P.count(adj_v) == 0) {
							node.vertex = adj_v;
							node.parent_vertex = u;
							node.priority_value = P_u + ec;
							Q_handles[adj_v] = Q.push(node);
							P[adj_v] = node.priority_value;
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
		}
	}


	auto end = std::chrono::high_resolution_clock::now();
	double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s

	return L;

}


#include <text mining/parse_substring_between_pairs_of_delimiters.h>

unordered_map<int, unordered_map<int, pair<double, int>>> graph_hash_of_mixed_weighted_PLL_Read_Indexes(string index_file_name) {

	/*indexes; L(v_1) = { v_x, <distance to v_x, predecessor of v_x in the SP between v_1 and v_x>}*/
	unordered_map<int, unordered_map<int, pair<double, int>>> L;

	std::string line_content;
	std::ifstream myfile(index_file_name); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		long long int line_id = 0;

		while (getline(myfile, line_content)) // read file line by line
		{
			line_id++;
			if (line_id > 6) {
				std::vector<std::string> Parsed_content = parse_string(line_content, ":");
				if (Parsed_content.size() == 2) {
					int v = stoi(Parsed_content[0]);
					L[v] = {};
					std::vector<string> index_strings = parse_substring_between_pairs_of_delimiters(Parsed_content[1], "<", ">");
					for (auto it = index_strings.begin(); it != index_strings.end(); it++) {
						std::vector<string> an_index = parse_string(*it, ",");
						L[v][stoi(an_index[0])] = { stod(an_index[1]), stoi(an_index[2]) };
						//cout << stoi(an_index[0]) << "," << stod(an_index[1]) << "," << stoi(an_index[2]) << endl;
						//getchar();
					}
				}
			}
		}
		myfile.close(); //close the file

		return L;
	}
	else
	{
		std::cout << "Unable to open file " << index_file_name << std::endl
			<< "Please check the file location or file name." << std::endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}



}


double graph_hash_of_mixed_weighted_PLL_extract_distance(unordered_map<int, unordered_map<int, pair<double, int>>>& L, int source, int terminal) {

	/*return std::numeric_limits<double>::max() is not connected*/

	double distance = INT_MAX;
	auto mm = L.find(source);
	auto nn = L.find(terminal);
	for (auto xx = mm->second.begin(); xx != mm->second.end(); xx++) {
		int v = xx->first;
		if (nn->second.count(v) > 0) {
			double dis = xx->second.first + nn->second[v].first;
			if (distance > dis) {
				distance = dis;
			}
		}
	}
	return distance;

}


vector<pair<int, int>> graph_hash_of_mixed_weighted_PLL_extract_shortest_path(unordered_map<int, unordered_map<int, pair<double, int>>>& L, int source, int terminal) {

	vector<pair<int, int>> path; // edges

	while (source != terminal) {

		int capped_v;

		double distance = INT_MAX;
		bool connected = false;
		auto mm = L.find(source);
		auto nn = L.find(terminal);
		for (auto xx = mm->second.begin(); xx != mm->second.end(); xx++) {
			int v = xx->first;
			if (nn->second.count(v) > 0) {
				connected = true;
				double dis = xx->second.first + nn->second[v].first;
				if (distance > dis) {
					distance = dis;
					capped_v = v;
				}
			}
		}

		if (connected) {
			if (source != mm->second[capped_v].second) {
				path.push_back({ source , mm->second[capped_v].second });
				source = mm->second[capped_v].second; // ascending from source
			}

			if (terminal != nn->second[capped_v].second) {
				path.push_back({ terminal , nn->second[capped_v].second });
				terminal = nn->second[capped_v].second; // ascending from terminal
			}
		}
		else {
			break;
		}

	}

	return path;
}


#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_generate_random_connected_graph.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_read_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_save_graph_with_weight.h>


void simple_iterative_tests_forPLL() {

	/*parameters*/
	int iteration_times = 1;
	int V = 1000, E = 10000, precision = 3;
	double ec_min = 1, ec_max = 10;

	/*iteration*/
	std::time_t now = std::time(0);
	boost::random::mt19937 gen{ static_cast<std::uint32_t>(now) };
	for (int i = 0; i < iteration_times; i++) {

		/*input and output*/
		int generate_new_graph = 1;

		std::unordered_set<int> generated_group_vertices;
		graph_hash_of_mixed_weighted instance_graph, generated_group_graph;
		if (generate_new_graph == 1) {
			instance_graph = graph_hash_of_mixed_weighted_generate_random_connected_graph(V, E, 0, 0, ec_min, ec_max, precision);
			graph_hash_of_mixed_weighted_save_graph_with_weight("simple_iterative_tests.txt", instance_graph, 0);
		}
		else {
			double lambda;
			graph_hash_of_mixed_weighted_read_graph_with_weight("simple_iterative_tests.txt", instance_graph, lambda);
		}

		graph_hash_of_mixed_weighted_PLL_Create_Indexes(instance_graph, "simple_iterative_tests_PLL_Indexes.txt");

		unordered_map<int, unordered_map<int, pair<double, int>>> L = graph_hash_of_mixed_weighted_PLL_Read_Indexes("simple_iterative_tests_PLL_Indexes.txt");

		boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(V - 1) };
		int source = dist(gen);

		std::unordered_map<int, double> distances;
		std::unordered_map<int, int> predecessors;
		graph_hash_of_mixed_weighted_shortest_paths_source_to_all(instance_graph, source, distances, predecessors);

		for (int xx = 0; xx < 10; xx++) {
			int terminal = dist(gen);
			double dis = graph_hash_of_mixed_weighted_PLL_extract_distance(L, source, terminal);
			if (abs(dis - distances[terminal]) > 1e-5) {
				cout << "dis = " << dis << endl;
				cout << "distances[terminal] = " << distances[terminal] << endl;
				cout << "abs(dis - distances[terminal]) > 1e-5!" << endl;
				getchar();
			}
			vector<pair<int, int>> path = graph_hash_of_mixed_weighted_PLL_extract_shortest_path(L, source, terminal);
			double path_dis = 0;
			if (path.size() == 0) {
				if (source != terminal) { // disconnected
					path_dis = std::numeric_limits<double>::max();
				}
			}
			else {
				for (auto it = path.begin(); it != path.end(); it++) {
					path_dis = path_dis + graph_hash_of_mixed_weighted_edge_weight(instance_graph, it->first, it->second);
					if (path_dis > std::numeric_limits<double>::max()) {
						path_dis = std::numeric_limits<double>::max();
					}
				}
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

#pragma endregion slow_hash_version


#include<graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_PLL_labels.h>

#pragma region
/*unordered_map<int, vector<PLL_sorted_label>> is faster than unordered_map<int, unordered_map<int, pair<double, int>>> in large and dense graphs*/
#include<graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID.h>

/*this is the fastest version*/
unordered_map<int, vector<PLL_sorted_label>> PLL_hash_of_vectors_generate_indexes(graph_hash_of_mixed_weighted& input_graph) {

	/*this function is different from (and faster than) PLL_hash_of_vectors_generate_indexes_basic 
	in that this function change vertexIDs to ideal IDs,
	and then use the rank and T P tricks in 2013 Japan paper*/

	int N = input_graph.hash_of_vectors.size();

	/*sort vertices by degrees*/
	vector<pair<int, int>> sorted_vertices;
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		sorted_vertices.push_back({ it->first, input_graph.degree(it->first) });
	}
	sort(sorted_vertices.begin(), sorted_vertices.end(), PLL_Creating_Indexes_compare);


	/*graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID*/
	unordered_map<int, int> vertexID_old_to_new;
	vector<int> vertexID_new_to_old(N);
	for (int i = 0; i < N; i++) {
		vertexID_old_to_new[sorted_vertices[i].first] = i;
		vertexID_new_to_old[i] = sorted_vertices[i].first;
	}
	graph_v_of_v_idealID ideal_graph = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(input_graph, vertexID_old_to_new);



	vector<vector<PLL_sorted_label>> L(N); // Labels for ideal IDs
	vector<double> P(N, std::numeric_limits<double>::max()); // P is distance to v_k, and priorities in Q
	vector<double> T(N, std::numeric_limits<double>::max());
	vector<PLL_handle_t_for_sp> Q_handles(N);

	/*PLL: Pruned Landmark Labeling*/
	for (int v_k = 0; v_k < N; v_k++) {

		/*Pruned Dijkstra from vertex v_k; see Algorithm 1 in 2013 Japan SIGMOD paper*/

		queue<int> P_changed_vertices;

		PLL_node_for_sp node;
		boost::heap::fibonacci_heap<PLL_node_for_sp> Q;
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
	unordered_map<int, vector<PLL_sorted_label>> unordered_map_L;
	for (int v_k = 0; v_k < N; v_k++) {
		vector<PLL_sorted_label> xx_vector = {};
		int L_v_k_size = L[v_k].size();
		for (int i = 0; i < L_v_k_size; i++) {
			PLL_sorted_label_binary_insert(xx_vector, vertexID_new_to_old[L[v_k][i].vertex], vertexID_new_to_old[L[v_k][i].parent_vertex], L[v_k][i].distance);
		}
		unordered_map_L[vertexID_new_to_old[v_k]] = xx_vector;
	}
	return unordered_map_L;

}

unordered_map<int, vector<PLL_sorted_label>> PLL_hash_of_vectors_generate_indexes_basic_slow(graph_hash_of_mixed_weighted& input_graph) {

	unordered_map<int, vector<PLL_sorted_label>> L;

	/*sort vertices by degrees*/
	vector<pair<int, int>> sorted_vertices;
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		sorted_vertices.push_back({ it->first, input_graph.degree(it->first) });
		L[it->first] = {};
	}
	sort(sorted_vertices.begin(), sorted_vertices.end(), PLL_Creating_Indexes_compare);


	/*PLL: Pruned Landmark Labeling*/
	for (auto it = sorted_vertices.begin(); it != sorted_vertices.end(); it++) {

		/*Pruned Dijkstra from vertex v_k = it->first; see Algorithm 1 in 2013 Japan SIGMOD paper*/
		int v_k = it->first;

		PLL_node_for_sp node;
		boost::heap::fibonacci_heap<PLL_node_for_sp> Q;
		std::unordered_map<int, double> P; // distance to v_k; priorities in Q
		std::unordered_map<int, PLL_handle_t_for_sp> Q_handles;

		node.vertex = v_k;
		node.parent_vertex = v_k;
		node.priority_value = 0;
		Q_handles[v_k] = Q.push(node);
		P[v_k] = 0;


		while (Q.size() > 0) {

			node = Q.top();
			Q.pop();
			int u = node.vertex;
			int u_parent = node.parent_vertex;
			double P_u = node.priority_value;


			/*merge-joinly compute query_v_k_u, which is the QUERY value in Line 7 of Algo. 1 in 2013 Japan paper*/
			double query_v_k_u = std::numeric_limits<double>::max(); // if disconnected, retun this large value
			auto vector1_pointer = L.find(u)->second; // a vector<PLL_sorted_label>
			auto vector2_pointer = L.find(v_k)->second; // another vector<PLL_sorted_label>
			auto vector1_check_pointer = vector1_pointer.begin();
			auto vector2_check_pointer = vector2_pointer.begin();
			while (vector1_check_pointer != vector1_pointer.end() && vector2_check_pointer != vector2_pointer.end()) {
				if (vector1_check_pointer->vertex == vector2_check_pointer->vertex) {
					double dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
					if (query_v_k_u > dis) {
						query_v_k_u = dis;
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

			if (P_u < query_v_k_u) { // this is pruning

				PLL_sorted_label_binary_insert(L[u], v_k, u_parent, P_u); // this cannot be PLL_sorted_label_binary_insert(vector1_pointer, v_k, u_parent, P_u), otherwise L[u] will not change any more

				auto search = input_graph.hash_of_hashs.find(u);
				if (search != input_graph.hash_of_hashs.end()) {
					for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
						int adj_v = it2->first;
						double ec = it2->second;
						if (P.count(adj_v) == 0) {
							node.vertex = adj_v;
							node.parent_vertex = u;
							node.priority_value = P_u + ec;
							Q_handles[adj_v] = Q.push(node);
							P[adj_v] = node.priority_value;
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
				else {
					auto search2 = input_graph.hash_of_vectors.find(u);
					for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
						int adj_v = it2->first;
						double ec = it2->second;
						if (P.count(adj_v) == 0) {
							node.vertex = adj_v;
							node.parent_vertex = u;
							node.priority_value = P_u + ec;
							Q_handles[adj_v] = Q.push(node);
							P[adj_v] = node.priority_value;
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
		}
	}

	return L;

}

/*this is the fastest version*/
void PLL_hash_of_vectors_generate_and_save_indexes(graph_hash_of_mixed_weighted& input_graph, string index_file_name) {

	auto begin1 = std::chrono::high_resolution_clock::now();

	unordered_map<int, vector<PLL_sorted_label>> L = PLL_hash_of_vectors_generate_indexes(input_graph);

	auto end1 = std::chrono::high_resolution_clock::now();
	double runningtime1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - begin1).count() / 1e9; // s


	

	auto begin2 = std::chrono::high_resolution_clock::now();

	long long int index_size = 0;
	for (auto it = L.begin(); it != L.end(); it++) {
		index_size = index_size + it->second.size();
	}

	/*save indexes*/
	std::ofstream outputFile;
	outputFile.precision(6);
	outputFile.setf(std::ios::fixed);
	outputFile.setf(std::ios::showpoint);
	outputFile.open(index_file_name);
	outputFile << "SECTION Comments" << std::endl;
	outputFile << "Creator Yahui SUN ËïÑÇ»Ô" << std::endl;
	outputFile << "Graph Size: |V|=" << input_graph.hash_of_vectors.size() << " |E|=" << graph_hash_of_mixed_weighted_num_edges(input_graph) << std::endl;
	outputFile << "Total Number of Indexes: " << index_size << " Consumed Time: " << runningtime1 << "s"
		<< " Consumed RAM: " << Current_Memory_Consumption_of_This_Process() << "MB" << std::endl;
	outputFile << "Number of Indexes Per Vertex: " << (double)index_size / input_graph.hash_of_vectors.size() << std::endl;
	outputFile << "Comments END at Line 6" << std::endl;
	for (auto it = L.begin(); it != L.end(); it++) {
		outputFile << it->first << ":"; // v_1:
		for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
			outputFile << "<" << it2->vertex << "," << it2->distance << "," << it2->parent_vertex << ">"; // <v_x,dis,pre>
		}
		outputFile << std::endl;
	}

	auto end2 = std::chrono::high_resolution_clock::now();
	double runningtime2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2).count() / 1e9; // s

	/*runningtime2 is around 1/5 runningtime1;
	
	for Amazon data on Linux: 
PLL_amazon.txt PLL runningtime1: 2882.16s
PLL_amazon.txt Saving runningtime2: 271.276s
*/	
	cout << index_file_name + " PLL runningtime1: " << runningtime1 << "s" << std::endl;
	cout << index_file_name + " Saving runningtime2: " << runningtime2 << "s" << std::endl;

}

unordered_map<int, vector<PLL_sorted_label>> PLL_hash_of_vectors_read_indexes(string index_file_name) {

	/*indexes; L(v_1) = { v_x, <distance to v_x, predecessor of v_x in the SP between v_1 and v_x>}*/
	unordered_map<int, vector<PLL_sorted_label>> L;

	std::string line_content;
	std::ifstream myfile(index_file_name); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		long long int line_id = 0;

		while (getline(myfile, line_content)) // read file line by line
		{
			line_id++;
			if (line_id > 6) {
				std::vector<std::string> Parsed_content = parse_string(line_content, ":");
				if (Parsed_content.size() == 2) {
					int v = stoi(Parsed_content[0]);
					L[v] = {};
					std::vector<string> index_strings = parse_substring_between_pairs_of_delimiters(Parsed_content[1], "<", ">");
					for (auto it = index_strings.begin(); it != index_strings.end(); it++) {
						std::vector<string> an_index = parse_string(*it, ",");
						PLL_sorted_label_binary_insert(L[v], stoi(an_index[0]), stoi(an_index[2]), stod(an_index[1]));
						//cout << stoi(an_index[0]) << "," << stod(an_index[1]) << "," << stoi(an_index[2]) << endl;
						//getchar();
					}
				}
			}
		}
		myfile.close(); //close the file

		return L;
	}
	else
	{
		std::cout << "Unable to open file " << index_file_name << std::endl
			<< "Please check the file location or file name." << std::endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}



}

double PLL_hash_of_vectors_extract_distance(unordered_map<int, vector<PLL_sorted_label>>& L, int source, int terminal) {

	/*return std::numeric_limits<double>::max() is not connected*/

	double distance = std::numeric_limits<double>::max(); // if disconnected, retun this large value
	auto vector1_pointer = L.find(source)->second; // a vector<PLL_sorted_label>
	auto vector2_pointer = L.find(terminal)->second; // another vector<PLL_sorted_label>
	auto vector1_check_pointer = vector1_pointer.begin();
	auto vector2_check_pointer = vector2_pointer.begin();
	while (vector1_check_pointer != vector1_pointer.end() && vector2_check_pointer != vector2_pointer.end()) {
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

vector<pair<int, int>> PLL_hash_of_vectors_extract_shortest_path(unordered_map<int, vector<PLL_sorted_label>>& L, int source, int terminal) {

	vector<pair<int, int>> path; // edges

	while (source != terminal) {

		double distance = std::numeric_limits<double>::max(); // if disconnected, retun this large value
		bool connected = false;
		auto vector1_pointer = L.find(source)->second; // a vector<PLL_sorted_label>
		auto vector2_pointer = L.find(terminal)->second; // another vector<PLL_sorted_label>
		auto vector1_check_pointer = vector1_pointer.begin();
		auto vector2_check_pointer = vector2_pointer.begin();
		int vector1_capped_v_parent = 0, vector2_capped_v_parent = 0;
		while (vector1_check_pointer != vector1_pointer.end() && vector2_check_pointer != vector2_pointer.end()) {
			if (vector1_check_pointer->vertex == vector2_check_pointer->vertex) {
				connected = true;
				double dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
				if (distance > dis) {
					distance = dis;
					vector1_capped_v_parent = vector1_check_pointer->parent_vertex;
					vector2_capped_v_parent = vector2_check_pointer->parent_vertex;
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


		if (connected) {
			/*the following code will not induce redundent edge, since for each vertex v_k, there is a label (v_k,0,v_k) in L[v_k]*/
			if (source != vector1_capped_v_parent) {
				path.push_back({ source , vector1_capped_v_parent });
				source = vector1_capped_v_parent; // ascending from source
			}
			if (terminal != vector2_capped_v_parent) {
				path.push_back({ terminal , vector2_capped_v_parent });
				terminal = vector2_capped_v_parent; // ascending from terminal
			}
		}
		else {
			break;
		}

	}

	return path;
}

#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_generate_random_graph.h>

void test_PLL_hash_of_vectors() {

	/*parameters*/
	int iteration_times = 100;
	int V = 500, E = 100, precision = 3;
	double ec_min = 1, ec_max = 10;

	/*iteration*/
	std::time_t now = std::time(0);
	boost::random::mt19937 gen{ static_cast<std::uint32_t>(now) };
	for (int i = 0; i < iteration_times; i++) {

		/*input and output*/
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

		PLL_hash_of_vectors_generate_and_save_indexes(instance_graph, "simple_iterative_tests_PLL_Indexes.txt");

		unordered_map<int, vector<PLL_sorted_label>> L = PLL_hash_of_vectors_read_indexes("simple_iterative_tests_PLL_Indexes.txt");

		boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(V - 1) };
		int source = dist(gen);

		std::unordered_map<int, double> distances;
		std::unordered_map<int, int> predecessors;
		graph_hash_of_mixed_weighted_shortest_paths_source_to_all(instance_graph, source, distances, predecessors);

		for (int xx = 0; xx < 10; xx++) {
			int terminal = dist(gen);
			double dis = PLL_hash_of_vectors_extract_distance(L, source, terminal);
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
			vector<pair<int, int>> path = PLL_hash_of_vectors_extract_shortest_path(L, source, terminal);
			double path_dis = 0;
			if (path.size() == 0) {
				if (source != terminal) { // disconnected
					path_dis = std::numeric_limits<double>::max();
				}
			}
			else {
				for (auto it = path.begin(); it != path.end(); it++) {
					path_dis = path_dis + graph_hash_of_mixed_weighted_edge_weight(instance_graph, it->first, it->second);
					if (path_dis > std::numeric_limits<double>::max()) {
						path_dis = std::numeric_limits<double>::max();
					}
				}
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


void test_PLL_hash_of_vectors_time() {

	/*parameters*/
	int iteration_times = 1;
	int V = 5000, E = 100000, precision = 3;
	double ec_min = 1, ec_max = 10;

	/*iteration*/
	std::time_t now = std::time(0);
	boost::random::mt19937 gen{ static_cast<std::uint32_t>(now) };
	for (int i = 0; i < iteration_times; i++) {

		/*input and output*/
		int generate_new_graph = 1;

		std::unordered_set<int> generated_group_vertices;
		graph_hash_of_mixed_weighted instance_graph, generated_group_graph;
		if (generate_new_graph == 1) {
			instance_graph = graph_hash_of_mixed_weighted_generate_random_connected_graph(V, E, 0, 0, ec_min, ec_max, precision);
			graph_hash_of_mixed_weighted_save_graph_with_weight("simple_iterative_tests.txt", instance_graph, 0);
		}
		else {
			double lambda;
			graph_hash_of_mixed_weighted_read_graph_with_weight("simple_iterative_tests.txt", instance_graph, lambda);
		}



		PLL_hash_of_vectors_generate_and_save_indexes(instance_graph, "simple_iterative_tests_PLL_Indexes.txt");



		unordered_map<int, vector<PLL_sorted_label>> L = PLL_hash_of_vectors_read_indexes("simple_iterative_tests_PLL_Indexes.txt");


	}
}

#pragma endregion hash_of_vector_version











#pragma region
/*vectorFORMAT; when vertex IDs are within a constant value max_N*/

vector<vector<PLL_sorted_label>> PLL_generate_indexes_vectorFORMAT(graph_hash_of_mixed_weighted& input_graph, int max_N) {

	/*this function is different from (and faster than) PLL_hash_of_vectors_generate_indexes_basic
	in that this function change vertexIDs to ideal IDs,
	and then use the rank and T P tricks in 2013 Japan paper*/

	int N = input_graph.hash_of_vectors.size();

	/*sort vertices by degrees*/
	vector<pair<int, int>> sorted_vertices;
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		sorted_vertices.push_back({ it->first, input_graph.degree(it->first) });
	}
	sort(sorted_vertices.begin(), sorted_vertices.end(), PLL_Creating_Indexes_compare);


	/*graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID*/
	unordered_map<int, int> vertexID_old_to_new;
	vector<int> vertexID_new_to_old(N);
	for (int i = 0; i < N; i++) {
		vertexID_old_to_new[sorted_vertices[i].first] = i;
		vertexID_new_to_old[i] = sorted_vertices[i].first;
	}
	graph_v_of_v_idealID ideal_graph = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(input_graph, vertexID_old_to_new);



	vector<vector<PLL_sorted_label>> L(N); // Labels for ideal IDs
	vector<double> P(N, std::numeric_limits<double>::max()); // P is distance to v_k, and priorities in Q
	vector<double> T(N, std::numeric_limits<double>::max());
	vector<PLL_handle_t_for_sp> Q_handles(N);

	/*PLL: Pruned Landmark Labeling*/
	for (int v_k = 0; v_k < N; v_k++) {

		/*Pruned Dijkstra from vertex v_k; see Algorithm 1 in 2013 Japan SIGMOD paper*/

		queue<int> P_changed_vertices;

		PLL_node_for_sp node;
		boost::heap::fibonacci_heap<PLL_node_for_sp> Q;
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
	vector<vector<PLL_sorted_label>> output_L(max_N);
	for (int v_k = 0; v_k < N; v_k++) {
		vector<PLL_sorted_label> xx_vector = {};
		int L_v_k_size = L[v_k].size();
		for (int i = 0; i < L_v_k_size; i++) {
			PLL_sorted_label_binary_insert(xx_vector, vertexID_new_to_old[L[v_k][i].vertex], vertexID_new_to_old[L[v_k][i].parent_vertex], L[v_k][i].distance);
		}
		output_L[vertexID_new_to_old[v_k]] = xx_vector;
	}
	return output_L;

}

#include <text mining/binary_save_read_vector_of_vectors.h>

void PLL_generate_and_save_indexes_vectorFORMAT_binary(graph_hash_of_mixed_weighted& input_graph, string index_file_name, int max_N) {

	auto begin1 = std::chrono::high_resolution_clock::now();

	vector<vector<PLL_sorted_label>> L = PLL_generate_indexes_vectorFORMAT(input_graph, max_N);

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
	outputFile << "Creator Yahui SUN ËïÑÇ»Ô" << std::endl;
	outputFile << "Graph Size: |V|=" << input_graph.hash_of_vectors.size() << " |E|=" << graph_hash_of_mixed_weighted_num_edges(input_graph) << std::endl;
	outputFile << "Total Number of Indexes: " << index_size << " Consumed Time: " << runningtime1 << "s"
		<< " Consumed RAM: " << Current_Memory_Consumption_of_This_Process() << "MB" << std::endl;
	outputFile << "Number of Indexes Per Vertex: " << (double)index_size / input_graph.hash_of_vectors.size() << std::endl;
	outputFile << "Comments END at Line 6" << std::endl;

	binary_save_vector_of_vectors(index_file_name, L);

	auto end2 = std::chrono::high_resolution_clock::now();
	double runningtime2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2).count() / 1e9; // s

	/*runningtime2 is around 1/5 runningtime1;

	for Amazon data on Linux:
PLL_amazon.txt PLL runningtime1: 2882.16s
PLL_amazon.txt Saving runningtime2: 271.276s
*/
	
	cout << index_file_name + " Saving runningtime2: " << runningtime2 << "s" << std::endl;

}

vector<vector<PLL_sorted_label>> PLL_generate_indexes_vectorFORMAT_one_edge_weight(graph_hash_of_mixed_weighted& input_graph, int max_N) {

	/*this function is different from (and faster than) PLL_hash_of_vectors_generate_indexes_basic
	in that this function change vertexIDs to ideal IDs,
	and then use the rank and T P tricks in 2013 Japan paper*/

	int N = input_graph.hash_of_vectors.size();

	/*sort vertices by degrees*/
	vector<pair<int, int>> sorted_vertices;
	for (auto it = input_graph.hash_of_vectors.begin(); it != input_graph.hash_of_vectors.end(); it++) {
		sorted_vertices.push_back({ it->first, input_graph.degree(it->first) });
	}
	sort(sorted_vertices.begin(), sorted_vertices.end(), PLL_Creating_Indexes_compare);


	/*graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID*/
	unordered_map<int, int> vertexID_old_to_new;
	vector<int> vertexID_new_to_old(N);
	for (int i = 0; i < N; i++) {
		vertexID_old_to_new[sorted_vertices[i].first] = i;
		vertexID_new_to_old[i] = sorted_vertices[i].first;
	}
	graph_v_of_v_idealID ideal_graph = graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID(input_graph, vertexID_old_to_new);



	vector<vector<PLL_sorted_label>> L(N); // Labels for ideal IDs
	vector<int> P(N, INT_MAX); // P is distance to v_k, and priorities in Q
	vector<int> T(N, INT_MAX);

	/*PLL: Pruned Landmark Labeling*/
	for (int v_k = 0; v_k < N; v_k++) {

		/*Pruned Dijkstra from vertex v_k; see Algorithm 1 in 2013 Japan SIGMOD paper*/

		queue<int> P_changed_vertices;
		queue<PLL_node_for_sp> Q;

		PLL_node_for_sp node;
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
	vector<vector<PLL_sorted_label>> output_L(max_N);
	for (int v_k = 0; v_k < N; v_k++) {
		vector<PLL_sorted_label> xx_vector = {};
		int L_v_k_size = L[v_k].size();
		for (int i = 0; i < L_v_k_size; i++) {
			PLL_sorted_label_binary_insert(xx_vector, vertexID_new_to_old[L[v_k][i].vertex], vertexID_new_to_old[L[v_k][i].parent_vertex], L[v_k][i].distance);
		}
		output_L[vertexID_new_to_old[v_k]] = xx_vector;
	}
	return output_L;

}

void PLL_generate_and_save_indexes_vectorFORMAT_binary_one_edge_weight(graph_hash_of_mixed_weighted& input_graph, string index_file_name, int max_N) {

	auto begin1 = std::chrono::high_resolution_clock::now();

	vector<vector<PLL_sorted_label>> L = PLL_generate_indexes_vectorFORMAT_one_edge_weight(input_graph, max_N);

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
	outputFile << "Creator Yahui SUN ËïÑÇ»Ô" << std::endl;
	outputFile << "Graph Size: |V|=" << input_graph.hash_of_vectors.size() << " |E|=" << graph_hash_of_mixed_weighted_num_edges(input_graph) << std::endl;
	outputFile << "Total Number of Indexes: " << index_size << " Consumed Time: " << runningtime1 << "s"
		<< " Consumed RAM: " << Current_Memory_Consumption_of_This_Process() << "MB" << std::endl;
	outputFile << "Number of Indexes Per Vertex: " << (double)index_size / input_graph.hash_of_vectors.size() << std::endl;
	outputFile << "Comments END at Line 6" << std::endl;

	binary_save_vector_of_vectors(index_file_name, L);

	auto end2 = std::chrono::high_resolution_clock::now();
	double runningtime2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - begin2).count() / 1e9; // s

	/*runningtime2 is around 1/5 runningtime1;

	for Amazon data on Linux:
PLL_amazon.txt PLL runningtime1: 2882.16s
PLL_amazon.txt Saving runningtime2: 271.276s
*/

	cout << index_file_name + " Saving runningtime2: " << runningtime2 << "s" << std::endl;

}



void test_PLL_vectorFORMAT() {

	/*parameters*/
	int iteration_graph_times = 100, iteration_source_times = 10, iteration_terminal_times = 10;
	int V = 200, E = 5000, precision = 3;
	double ec_min = 0.001, ec_max = 1;

	double avg_index_size_per_v = 0;

	/*iteration*/
	std::time_t now = std::time(0);
	boost::random::mt19937 gen{ static_cast<std::uint32_t>(now) };
	for (int i = 0; i < iteration_graph_times; i++) {

		/*input and output*/
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

		PLL_generate_and_save_indexes_vectorFORMAT_binary(instance_graph, "simple_iterative_tests_PLL_Indexes.txt", V + 1);
		vector<vector<PLL_sorted_label>> L = PLL_read_indexes_vectorFORMAT_binary("simple_iterative_tests_PLL_Indexes.txt");

		//vector<vector<PLL_sorted_label>> L = PLL_generate_indexes_vectorFORMAT(instance_graph, V + 1);

		long long int index_size = 0;
		for (auto it = L.begin(); it != L.end(); it++) {
			index_size = index_size + (*it).size();
		}
		avg_index_size_per_v = avg_index_size_per_v + (double)index_size / V / iteration_graph_times;

		boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(V - 1) };

		for (int yy = 0; yy < iteration_source_times; yy++) {

			int source = dist(gen);

			std::unordered_map<int, double> distances;
			std::unordered_map<int, int> predecessors;
			graph_hash_of_mixed_weighted_shortest_paths_source_to_all(instance_graph, source, distances, predecessors);

			for (int xx = 0; xx < iteration_terminal_times; xx++) {
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
				if (path.size() == 0) {
					if (source != terminal) { // disconnected
						path_dis = std::numeric_limits<double>::max();
					}
				}
				else {
					for (auto it = path.begin(); it != path.end(); it++) {
						path_dis = path_dis + graph_hash_of_mixed_weighted_edge_weight(instance_graph, it->first, it->second);
						if (path_dis > std::numeric_limits<double>::max()) {
							path_dis = std::numeric_limits<double>::max();
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

	cout << "avg_index_size_per_v: " << avg_index_size_per_v << endl;
}

void test_PLL_vectorFORMAT_time() {

	/*parameters*/
	int iteration_times = 100;
	int V = 500, E = 1000, precision = 1;
	double ec_min = 1, ec_max = 10;

	double time_avg = 0, L_size_avg = 0;

	/*iteration*/
	std::time_t now = std::time(0);
	boost::random::mt19937 gen{ static_cast<std::uint32_t>(now) };
	for (int i = 0; i < iteration_times; i++) {

		/*input and output*/
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



		//PLL_generate_and_save_indexes_vectorFORMAT_binary(instance_graph, "simple_iterative_tests_PLL_Indexes.txt", V + 1);


		auto begin = std::chrono::high_resolution_clock::now();
		vector<vector<PLL_sorted_label>> L = PLL_generate_indexes_vectorFORMAT(instance_graph, V+1);
		auto end = std::chrono::high_resolution_clock::now();
		double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
		time_avg = time_avg + (double)runningtime / iteration_times;

		long long int index_size = 0;
		for (auto it = L.begin(); it != L.end(); it++) {
			index_size = index_size + (*it).size();
		}
		L_size_avg = L_size_avg + (double)index_size / iteration_times;
	}

	cout << "time_avg = " << time_avg << "s" << endl;
	cout << "L_size_avg = " << L_size_avg << endl;
}











#include "para.h"


void test_PLL_vectorFORMAT_one_edge_weight() {

	/*parameters*/
	int iteration_graph_times = 1, iteration_source_times = 10, iteration_terminal_times = 10;
	int V = 1e4, E = 1e5, precision = 1;
	double ec_min = 1, ec_max = 1;

	double avg_index_size_per_v = 0;

	/*iteration*/
	std::time_t now = std::time(0);
	boost::random::mt19937 gen{ static_cast<std::uint32_t>(now) };
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
		//PLL_generate_and_save_indexes_vectorFORMAT_binary_one_edge_weight(instance_graph, "simple_iterative_tests_PLL_Indexes.txt", V + 1);
		///*below is for reading PLL labels*/
		//vector<vector<PLL_sorted_label>> L = PLL_read_indexes_vectorFORMAT_binary("simple_iterative_tests_PLL_Indexes.txt");


		vector<vector<PLL_sorted_label>> L = PLL_generate_indexes_vectorFORMAT_one_edge_weight(instance_graph, V + 1);

		//graph_hash_of_mixed_weighted_print(instance_graph);

		/*below is for checking whether the above labels are right (by randomly computing shortest paths)*/
		long long int index_size = 0;
		for (auto it = L.begin(); it != L.end(); it++) {
			index_size = index_size + (*it).size();
		}
		avg_index_size_per_v = avg_index_size_per_v + (double)index_size / V / iteration_graph_times;

		boost::random::uniform_int_distribution<> dist{ static_cast<int>(0), static_cast<int>(V - 1) };

		for (int yy = 0; yy < iteration_source_times; yy++) {
			int source = dist(gen);
			std::unordered_map<int, double> distances;
			std::unordered_map<int, int> predecessors;
			graph_hash_of_mixed_weighted_shortest_paths_source_to_all(instance_graph, source, distances, predecessors);

			for (int xx = 0; xx < iteration_terminal_times; xx++) {
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
				if (path.size() == 0) {
					if (source != terminal) { // disconnected
						path_dis = std::numeric_limits<double>::max();
					}
				}
				else {
					for (auto it = path.begin(); it != path.end(); it++) {
						path_dis = path_dis + graph_hash_of_mixed_weighted_edge_weight(instance_graph, it->first, it->second);
						if (path_dis > std::numeric_limits<double>::max()) {
							path_dis = std::numeric_limits<double>::max();
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

	cout << "avg_index_size_per_v: " << avg_index_size_per_v << endl;
}




#pragma endregion vectorFORMAT and binary




void compare_PLL_hash_of_hashes_vs_hash_of_vectors_new() {

	/*parameters*/
	int iteration_times = 2;
	int V = 2000, E = 10000, precision = 3;
	double ec_min = 1, ec_max = 1;

	double time_old_PLL_avg = 0, time_new_PLL_avg = 0, time_new2_PLL_avg = 0;

	/*iteration*/
	std::time_t now = std::time(0);
	boost::random::mt19937 gen{ static_cast<std::uint32_t>(now) };
	for (int i = 0; i < iteration_times; i++) {

		/*input and output*/
		int generate_new_graph = 1;

		std::unordered_set<int> generated_group_vertices;
		graph_hash_of_mixed_weighted instance_graph, generated_group_graph;
		if (generate_new_graph == 1) {
			instance_graph = graph_hash_of_mixed_weighted_generate_random_connected_graph(V, E, 0, 0, ec_min, ec_max, precision);
			graph_hash_of_mixed_weighted_save_graph_with_weight("simple_iterative_tests.txt", instance_graph, 0);
		}
		else {
			double lambda;
			graph_hash_of_mixed_weighted_read_graph_with_weight("simple_iterative_tests.txt", instance_graph, lambda);
		}

		/*old PLL code*/
		if (1) {
			auto begin = std::chrono::high_resolution_clock::now();
			unordered_map<int, unordered_map<int, pair<double, int>>> L = graph_hash_of_mixed_weighted_PLL_Generate_Indexes_no_save(instance_graph);
			auto end = std::chrono::high_resolution_clock::now();
			double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
			time_old_PLL_avg = time_old_PLL_avg + (double)runningtime / iteration_times;
		}

		/*new PLL code*/
		if (1) {
			auto begin = std::chrono::high_resolution_clock::now();
			unordered_map<int, vector<PLL_sorted_label>> L = PLL_hash_of_vectors_generate_indexes(instance_graph);
			auto end = std::chrono::high_resolution_clock::now();
			double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
			time_new_PLL_avg = time_new_PLL_avg + (double)runningtime / iteration_times;
		}

	}

	cout << "time_old_PLL_avg = " << time_old_PLL_avg << "s" << endl;
	cout << "time_new_PLL_avg = " << time_new_PLL_avg << "s" << endl;

}



