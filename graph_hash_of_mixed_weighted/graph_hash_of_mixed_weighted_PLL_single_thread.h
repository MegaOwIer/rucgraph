#pragma once

#include <Current_Memory_Consumption_of_This_Process.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_shortest_paths.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_PLL_labels.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_to_graph_v_of_v_idealID.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_generate_random_graph.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_generate_random_connected_graph.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_read_graph_with_weight.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_save_graph_with_weight.h>
#include <boost/heap/fibonacci_heap.hpp> 

bool PLL_Creating_Indexes_compare(const pair<int, int>& i, pair<int, int>& j)
{
	/*< is nearly 10 times slower than >*/
	return i.second > j.second;  // < is from small to big; > is from big to small.  sort by the second item of pair<int, int>
}

struct PLL_node_for_sp {
public:
	int vertex, parent_vertex;
	double priority_value; // double float do not change the PLL running time
}; // define the node in the queue
bool operator<(PLL_node_for_sp const& x, PLL_node_for_sp const& y) {
	return x.priority_value > y.priority_value; // < is the max-heap; > is the min heap
}
typedef typename boost::heap::fibonacci_heap<PLL_node_for_sp>::handle_type PLL_handle_t_for_sp;



#pragma region

/*the following two files are for generating PLL labels using single thread*/
vector<vector<PLL_sorted_label>> PLL_generate_indexes_weighted_single_thread(graph_hash_of_mixed_weighted& input_graph, int max_N) {

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

vector<vector<PLL_sorted_label>> PLL_generate_indexes_unweighted_single_thread(graph_hash_of_mixed_weighted& input_graph, int max_N) {

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


/*the following two files are for generating and saving PLL labels using binary files*/
#include <text mining/binary_save_read_vector_of_vectors.h>
void PLL_generate_and_save_indexes_weighted_single_thread(graph_hash_of_mixed_weighted& input_graph, string index_file_name, int max_N) {

	auto begin1 = std::chrono::high_resolution_clock::now();

	vector<vector<PLL_sorted_label>> L = PLL_generate_indexes_weighted_single_thread(input_graph, max_N);

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

void PLL_generate_and_save_indexes_unweighted_single_thread(graph_hash_of_mixed_weighted& input_graph, string index_file_name, int max_N) {

	auto begin1 = std::chrono::high_resolution_clock::now();

	vector<vector<PLL_sorted_label>> L = PLL_generate_indexes_unweighted_single_thread(input_graph, max_N);

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

#pragma endregion single thread PLL codes



#pragma region
/*
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
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_PLL_single_thread.h>

int main()
{
	test_PLL_generate_indexes_single_thread();
}
------------------------------------------------------------------------------------------
Commends for running the above cpp file on Linux:

g++ -std=c++17 -I/home/boost_1_75_0 -I/root/ysgraph try.cpp -lpthread -Ofast -o A
./A
rm A

(optional to put the above commends in run.sh, and then use the comment: sh run.sh)
*/


void check_correctness_of_PLL_labels(vector<vector<PLL_sorted_label>>& L, graph_hash_of_mixed_weighted& instance_graph, int iteration_source_times, int iteration_terminal_times, boost::random::mt19937& boost_random_time_seed) {

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


void test_PLL_generate_indexes_single_thread() {

	/*parameters*/
	int iteration_graph_times = 1, iteration_source_times = 10, iteration_terminal_times = 10;
	int V = 10, E = 15, precision = 1;
	double ec_min = 1, ec_max = 1; // set ec_min=ec_max=1 for testing unweighted PLL
	double avg_index_size_per_v = 0;

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
		//PLL_generate_and_save_indexes_weighted_single_thread(instance_graph, "simple_iterative_tests_PLL_Indexes.txt", V + 1);
		///*below is for reading PLL labels*/
		//vector<vector<PLL_sorted_label>> L = PLL_read_indexes_vectorFORMAT_binary("simple_iterative_tests_PLL_Indexes.txt");


		vector<vector<PLL_sorted_label>> L = PLL_generate_indexes_unweighted_single_thread(instance_graph, V + 1);

		check_correctness_of_PLL_labels(L, instance_graph, iteration_source_times, iteration_terminal_times, boost_random_time_seed);

		long long int index_size = 0;
		for (auto it = L.begin(); it != L.end(); it++) {
			index_size = index_size + (*it).size();
		}
		avg_index_size_per_v = avg_index_size_per_v + (double)index_size / V / iteration_graph_times;

	}

	cout << "avg_index_size_per_v: " << avg_index_size_per_v << endl;
}





#pragma endregion test code