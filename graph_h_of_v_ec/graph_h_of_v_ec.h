#pragma once

#include <unordered_map>
#include <print_items.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_binary_operations.h> // this works for sorted vectors of pair<int, double>


typedef std::unordered_map<int, std::vector<pair<int, double>>> graph_h_of_v_ec; // hash of vectors for all vertices (only edge weight, no vertex weight)



void graph_h_of_v_ec_add_vertex(graph_h_of_v_ec& input_graph, int vertex) {

	/*time complexity O(1);

	since unordered_map containers do not allow for duplicate keys, all the vertices are unique */

	if (input_graph.count(vertex) == 0) { // only add vertex if it is not in graph
		input_graph[vertex] = {};
	}
}



void graph_h_of_v_ec_remove_vertex(graph_h_of_v_ec& input_graph, int vertex) {

	/*this function removes a vertex and its adjacent edges;

	time complexity is close to O(vertex_degree * average degree of adjacenct vertices of vertex);

	worse-case time complexity O(|V|^2), use this as less as possible*/

	auto search1 = input_graph.find(vertex);
	if (search1 != input_graph.end()) { // vertex is in input_graph
		for (int i = search1->second.size() - 1; i >= 0; i--) {
			graph_hash_of_mixed_weighted_binary_operations_erase(input_graph[search1->second[i].first], vertex); // remove vertex from adj_list of search1->second[i].first (search1->second[i].first is adj_v)
		}
		input_graph.erase(vertex); // remove vertex
	}
}


bool graph_h_of_v_ec_contain_vertex(graph_h_of_v_ec& input_graph, int vertex) {

	/*this function checks whether a vertex/key is in input_graph or not;
	time complexity O(1)*/

	if (input_graph.count(vertex) > 0) {
		return true;
	}
	return false;
}


void graph_h_of_v_ec_add_edge(graph_h_of_v_ec& input_graph, int e1, int e2, double ec) {

	/*this function adds a weighted edge via binary search,
	and may add e1 and e2 into input_graph if they are new vertices;

	this function can update edge weights;

	time complexity O(largest_degree of e1 and e2);  worst-case time complexity: O(|V|)
	*/

	/*add e2 to adj of e1*/
	auto search = input_graph.find(e1);
	if (search == input_graph.end()) { // e1 is a new vertex
		graph_h_of_v_ec_add_vertex(input_graph, e1); // add e1
		search = input_graph.find(e1);
	}
	graph_hash_of_mixed_weighted_binary_operations_insert(search->second, e2, ec);

	/*add e1 to adj of e2*/
	search = input_graph.find(e2);
	if (search == input_graph.end()) { // e2 is a new vertex
		graph_h_of_v_ec_add_vertex(input_graph, e2); // add e2; initial weight is 0
		search = input_graph.find(e2);
	}
	graph_hash_of_mixed_weighted_binary_operations_insert(search->second, e1, ec);
}


double graph_h_of_v_ec_edge_weight(graph_h_of_v_ec& input_graph, int v1, int v2) {

	/*this function return an edge weight;

	we assume that this edge is in the graph, otherwise large value is returned;

	time complexity O(log (largest degree of v1 and v2))*/

	auto search3 = input_graph.find(v1);
	if (search3 != input_graph.end()) { // v1 is in input_graph
		int pos = graph_hash_of_mixed_weighted_binary_operations_search_position(search3->second, v2);
		if (pos != -1) {
			return search3->second[pos].second;
		}
		else {
			return INT_MAX;
		}
	}
	else {
		return INT_MAX;
	}
}


bool graph_h_of_v_ec_contain_edge(graph_h_of_v_ec& input_graph, int v1, int v2) {

	/*this function checks whether an edge is in input_graph or not;

	 time complexity O(log (largest degree of v1 and v2))*/

	if (input_graph.count(v1) == 0 || input_graph.count(v2) == 0) {
		return false;
	}

	if (graph_hash_of_mixed_weighted_binary_operations_search(input_graph[v1], v2)) {
		return true;
	}
	return false;
}


void graph_h_of_v_ec_remove_edge_and_isolated_vertices(graph_h_of_v_ec& input_graph, int e1, int e2) {

	/*this function removes an edge, and may remove e1 and e2 if they are isolated;

	time complexity O(log (largest degree of v1 and v2))*/


	auto search2 = input_graph.find(e1);
	if (search2 != input_graph.end()) { // e1 is in hash_of_vectors
		graph_hash_of_mixed_weighted_binary_operations_erase(search2->second, e2);  // remove e2 from adj of e1
	}
	if (search2->second.size() == 0) { // e1 is isolated
		input_graph.erase(e1); // remove e1
	}

	search2 = input_graph.find(e2);
	if (search2 != input_graph.end()) { // e2 is in hash_of_vectors
		graph_hash_of_mixed_weighted_binary_operations_erase(search2->second, e1);  // remove e1 from adj of e2
	}
	if (search2->second.size() == 0) { // e2 is isolated
		input_graph.erase(e2); // remove e2
	}
}


void graph_h_of_v_ec_remove_edge_but_not_isolated_vertices(graph_h_of_v_ec& input_graph, int e1, int e2) {

	/*this function removes an edge, and does not remove e1 and e2 if they are isolated;

	time complexity O(log (largest degree of v1 and v2))*/


	auto search2 = input_graph.find(e1);
	if (search2 != input_graph.end()) { // e1 is in hash_of_vectors
		graph_hash_of_mixed_weighted_binary_operations_erase(search2->second, e2);  // remove e2 from adj of e1
	}

	search2 = input_graph.find(e2);
	if (search2 != input_graph.end()) { // e2 is in hash_of_vectors
		graph_hash_of_mixed_weighted_binary_operations_erase(search2->second, e1);  // remove e1 from adj of e2
	}

}


long long int graph_h_of_v_ec_num_edges(graph_h_of_v_ec& input_graph) {

	/*time complexity O(|V|)*/

	long long int num = 0;
	for (auto it = input_graph.begin(); it != input_graph.end(); ++it) {
		int v = it->first;
		num = num + it->second.size();
	}
	return num / 2;
}


void graph_h_of_v_ec_print(graph_h_of_v_ec& input_graph) {

	/*this function prints the whole graph_h_of_v_ec (with bounded size)*/

	int print_relex_num = 100;
	int num = 0;
	//std::cout << '\n';
	std::cout << "graph_h_of_v_ec_print: |V|=" << input_graph.size() << " |E|=" << graph_h_of_v_ec_num_edges(input_graph) << std::endl;

	for (auto it = input_graph.begin(); it != input_graph.end(); ++it) {
		std::cout << "Vertex: " << it->first << " adj_vertices: ";
		for (int i = 0; i < it->second.size(); i++) {
			std::cout << "<" << it->second[i].first << "," << it->second[i].second << "> ";
		}
		std::cout << '\n';

		num++;
		if (num % print_relex_num == 0) {
			getchar();
		}
	}

	std::cout << "Print END" << '\n';
	std::cout << '\n';
}


void graph_h_of_v_ec_print_size(graph_h_of_v_ec& input_graph) {

	int num = 0;
	//std::cout << '\n';
	std::cout << "graph_h_of_v_ec_print_size: |V|=" << input_graph.size() << " |E|=" << graph_h_of_v_ec_num_edges(input_graph) << std::endl;
}





/*the following are test codes*/

void test_graph_h_of_v_ec() {

	graph_h_of_v_ec g;

	graph_h_of_v_ec_add_vertex(g, 1);
	graph_h_of_v_ec_add_vertex(g, 2);
	graph_h_of_v_ec_add_vertex(g, 3);
	graph_h_of_v_ec_add_vertex(g, 4);

	graph_h_of_v_ec_add_edge(g, 1, 12, 0.5);
	graph_h_of_v_ec_add_edge(g, 1, 2, 0.5);
	graph_h_of_v_ec_add_edge(g, 4, 2, 0.5);
	graph_h_of_v_ec_add_edge(g, 12, 1, 0.5);
	graph_h_of_v_ec_add_edge(g, 1, 12, 0.5);
	graph_h_of_v_ec_add_edge(g, 1, 4, 0.5);
	graph_h_of_v_ec_add_edge(g, 1, 3, 0.5);

	graph_h_of_v_ec_print(g);

	cout << "graph_h_of_v_ec_contain_vertex(g,2): " << graph_h_of_v_ec_contain_vertex(g, 2) << endl;
	graph_h_of_v_ec_remove_vertex(g, 2);
	cout << "graph_h_of_v_ec_contain_vertex(g,2): " << graph_h_of_v_ec_contain_vertex(g, 2) << endl;

	cout << "graph_h_of_v_ec_edge_weight(g,1,2): " << graph_h_of_v_ec_edge_weight(g, 1, 2) << endl;
	cout << "graph_h_of_v_ec_edge_weight(g,1,3): " << graph_h_of_v_ec_edge_weight(g, 1, 3) << endl;

	cout << "graph_h_of_v_ec_contain_edge(g,1,2): " << graph_h_of_v_ec_contain_edge(g, 1, 2) << endl;
	cout << "graph_h_of_v_ec_contain_edge(g,1,3): " << graph_h_of_v_ec_contain_edge(g, 1, 3) << endl;

	graph_h_of_v_ec_print(g);

	graph_h_of_v_ec_remove_edge_but_not_isolated_vertices(g, 1, 4);
	graph_h_of_v_ec_remove_edge_but_not_isolated_vertices(g, 1, 3);
	graph_h_of_v_ec_remove_edge_but_not_isolated_vertices(g, 1, 3);
	graph_h_of_v_ec_print(g);

	graph_h_of_v_ec_remove_edge_and_isolated_vertices(g, 1, 4);
	graph_h_of_v_ec_remove_edge_and_isolated_vertices(g, 1, 3);
	graph_h_of_v_ec_remove_edge_and_isolated_vertices(g, 1, 4);


	graph_h_of_v_ec_print(g);

}