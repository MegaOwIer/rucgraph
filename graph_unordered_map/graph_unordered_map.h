#pragma once

/*this is an adjacency list built using hashes*/

#include <unordered_map>
#include <unordered_set>

/* define subgraph*/
typedef std::unordered_map<int, double> graph_node_adj_hash;
// key: an adjacent node; value: weight of this adjacent node
typedef std::unordered_map<int, pair<double, graph_node_adj_hash>> graph_unordered_map;
// key: a node_ID; value: pair <node_weight, the hash of the adjacent nodes>



void graph_unordered_map_add_vertex(graph_unordered_map& input_graph, int vertex, double weight) {

	/*time complexity O(1)*/
	/*since unordered_map containers do not allow for duplicate keys,
	all the vertices in graph_unordered_map are unique*/

	/*this can be used to update vertex weight*/

	if (input_graph.count(vertex) > 0) {
		input_graph[vertex].first = weight;
	}
	else {
		pair<double, graph_node_adj_hash> x;
		x.first = weight;
		input_graph[vertex] = x;
	}
	
}


void graph_unordered_map_remove_vertex(graph_unordered_map& input_graph, int vertex) {

	/*this function removes a vertex and its adjacent edges;
	time complexity O(degree of vertex), which is O(|V|) in the worst case*/

	/*remove vertex from hashes of adj vetices*/
	for (auto it = input_graph[vertex].second.begin(); it != input_graph[vertex].second.end(); ++it) {
		input_graph[it->first].second.erase(vertex);
	}
	input_graph.erase(vertex); // remove vertex from input_graph

}


void graph_unordered_map_remove_leaf(graph_unordered_map& input_graph, int vertex) {

	/*this function removes a leaf and its adjacent edge;
	time complexity O(1)*/

	int adj_v = input_graph[vertex].second.begin()->first; // the adjacency vertex
	input_graph.erase(vertex); // remove this vertex
	input_graph[adj_v].second.erase(vertex);  // remove the adjacency edge

}


bool graph_unordered_map_contain_vertex(graph_unordered_map& input_graph, int vertex) {

	/*this function checks whether a vertex/key is in input_graph or not;
	time complexity O(1)*/

	if (input_graph.count(vertex) > 0) {
		return true;
	}
	else {
		return false;
	}

}


void graph_unordered_map_add_edge(graph_unordered_map& input_graph, int e1, int e2, double ec) {

	/*this function adds a weighted edge, and may add e1 and e2 into input_graph if they are new vertices;
	time complexity O(1)*/

	/*this can be used to update edge weight*/

	if (!graph_unordered_map_contain_vertex(input_graph, e1)) { // e1 is a new vertex
		graph_unordered_map_add_vertex(input_graph, e1, 0); // add e1; initial weight is 0
	}
	if (!graph_unordered_map_contain_vertex(input_graph, e2)) { // e2 is a new vertex
		graph_unordered_map_add_vertex(input_graph, e2, 0); // add e2; initial weight is 0
	}

	input_graph[e1].second[e2] = ec; // add e2 to adj_hash of e1
	input_graph[e2].second[e1] = ec; // add e1 to adj_hash of e2

}


bool graph_unordered_map_contain_edge(graph_unordered_map& input_graph, int v1, int v2) {

	/*this function checks whether an edge is in input_graph or not;
	time complexity O(1)*/

	if (input_graph.count(v1) == 0 || input_graph.count(v2) == 0) {
		return false;
	}


	if (input_graph[v1].second.count(v2) > 0) {
		return true;
	}
	else {
		return false;
	}

}


void graph_unordered_map_remove_edge_and_isolated_vertices(graph_unordered_map& input_graph, int e1, int e2) {

	/*this function removes an edge, and may remove e1 and e2 if they are isolated;
	time complexity O(1)*/

	auto search = input_graph.find(e1);
	search->second.second.erase(e2);
	if (input_graph[e1].second.size() == 0) {
		input_graph.erase(e1);
	}

	search = input_graph.find(e2);
	search->second.second.erase(e1);
	if (input_graph[e2].second.size() == 0) {
		input_graph.erase(e2);
	}

}


void graph_unordered_map_remove_edge_but_not_isolated_vertices(graph_unordered_map& input_graph, int e1, int e2) {

	/*this function removes an edge, but not remove e1 and e2 if they are isolated;
	time complexity O(1)*/

	auto search = input_graph.find(e1);
	search->second.second.erase(e2);
	search = input_graph.find(e2);
	search->second.second.erase(e1);


	/*input_graph[e1].second.erase[e2];
	input_graph[e2].second.erase[e1];
	don't know why the above codes do not work;
	the reason is that you use count[], not count(), stupid!*/

}

graph_unordered_map graph_unordered_map_copy_graph(graph_unordered_map& input_graph) {

	return input_graph;

}



int graph_unordered_map_num_vertices(graph_unordered_map& input_graph) {

	/*time complexity O(1)*/

	return input_graph.size();

}


int graph_unordered_map_num_edges(graph_unordered_map& input_graph) {

	/*time complexity O(|V|)*/

	int num = 0;
	for (auto it = input_graph.begin(); it != input_graph.end(); ++it) {
		num = num + it->second.second.size();
	}

	return num / 2;

}


void graph_unordered_map_print(graph_unordered_map& input_graph) {

	/*time complexity O(|V|+|E|)*/

	int print_relex_num = 100;
	int num = 0;
	//std::cout << '\n';
	std::cout << "graph_unordered_map_print: " << std::endl;
	std::cout << "|V|=" << graph_unordered_map_num_vertices(input_graph) << 
		" |E|=" << graph_unordered_map_num_edges(input_graph) << std::endl;

	for (auto it = input_graph.begin(); it != input_graph.end(); ++it) {
		std::cout << "Vertex: " << it->first << " Weight: " << it->second.first << '\n';
		num++;
		if (num%print_relex_num == 0) {
			getchar();
		}
	}


	for (auto it = input_graph.begin(); it != input_graph.end(); ++it) {
		for (auto it2 = it->second.second.begin(); it2 != it->second.second.end(); ++it2) {
			if (it->first <= it2->first) {
				std::cout << "Edge: (" << it->first << "," << it2->first<<") Cost: " << it2->second << '\n';
				num++;
				if (num%print_relex_num == 0) {
					getchar();
				}
			}
		}
	}


	std::cout << "graph_unordered_map_print END" << '\n';
	std::cout << '\n';
}


std::list<int> graph_unordered_map_adjacent_vertices(graph_unordered_map& input_graph, int vertex) {

	/*this function returns the adjacent vertices of the input vertex;
	time complexity O(degree of vertex), which is O(|V|) in the worst case*/

	std::list<int> adj_list;

	auto search = input_graph.find(vertex);
	//std::cout << "vertex: " << vertex << std::endl;
	//std::cout << "search->second.size(): " << search->second.size() << std::endl;
	//subgraph_print(input_graph);
	if (search != input_graph.end()) {
		if (search->second.second.size() > 0) {
			for (auto it = search->second.second.begin(); it != search->second.second.end(); ++it) {
				adj_list.push_back(it->first);
			}
		}
	}

	return adj_list;

}


int graph_unordered_map_degree(graph_unordered_map& input_graph, int vertex) {

	/*time complexity O(1)*/

	if (input_graph.count(vertex) > 0) {
		return input_graph[vertex].second.size();
	}

	return 0;

}


int graph_unordered_map_vertex_list_ID(graph_unordered_map& input_graph, int vertex) {

	/*this function returns the ID of the hash that corresponds to the vertex;
	time complexity O(|V|)*/
	int i = 0;
	for (auto it = input_graph.begin(); it != input_graph.end(); it++) {
		if ((*it).first == vertex) {
			return i;
		}
		i++;
	}
	return INT_MAX;
}



