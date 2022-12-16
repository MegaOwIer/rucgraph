#pragma once

/*this is an adjacency vector; it is less efficient than a list for Insertion and Deletion*/



/* define subgraph: pair<vertex, adjacent vertices> */
typedef std::pair<int, std::vector<int>> subgraph_vector_node;
typedef std::vector<std::pair<int, std::vector<int>>> subgraph_vector; // you initially just call it subgraph



#pragma region
void subgraph_vector_add_vertex(subgraph_vector& input_graph, int vertex) {

	/*this function only adds a unique vertex; time complexity O(|V|)*/

	/*time complexity O(|V|)*/
	bool unique = true;
	int i = 0;
	while (i < input_graph.size()) {
		if (input_graph[i].first == vertex) { // this vertex is already inside
			unique = false;
			break;
		}
		i++;
	}

	if (unique == true) {
		subgraph_vector_node x;
		x.first = vertex;
		input_graph.insert(input_graph.end(), x); // insert a subgraph_node
	}

}
#pragma endregion subgraph_vector_add_vertex


#pragma region
void subgraph_vector_add_nonunique_vertex(subgraph_vector& input_graph, int vertex) {

	/*this function may add a nonunique vertex; time complexity O(1)*/

	subgraph_vector_node x;
	x.first = vertex;
	input_graph.insert(input_graph.end(), x); // insert a subgraph_node

}
#pragma endregion subgraph_vector_add_nonunique_vertex



#pragma region
void subgraph_vector_remove_vertex(subgraph_vector& input_graph, int vertex) {

	/*time complexity O(|V|+|E|)*/

	int i = 0;
	while (i < input_graph.size()) {
		if (input_graph[i].first == vertex) { // i is the list of this vertex
			input_graph.erase(input_graph.begin() + i); // remove this list
			i--;
		}
		else {
			/*remove vertex from input_graph[i].second*/
			int j = 0;
			while (j < input_graph[i].second.size()) {
				if (input_graph[i].second[j] == vertex) {
					input_graph[i].second.erase(input_graph[i].second.begin() + j); // remove vertex from input_graph[i].second
					break; // move on from this list
				}
				j++;
			}
		}
		i++;
	}

}
#pragma endregion subgraph_vector_remove_vertex



#pragma region
void subgraph_vector_add_edge(subgraph_vector& input_graph, int e1, int e2) {

	/*this function only adds a unique edge, and may add e1 and e2 into input_graph if they are new vertices;
	time complexity O(|V|+|E|)*/

	/*time complexity O(|V|+|E|)*/
	bool unique_edge = true;
	bool e1_is_new_vertex = true;
	int i = 0;
	while (i < input_graph.size()) {
		if (input_graph[i].first == e1) {
			e1_is_new_vertex = false;
			for (int j = 0; j < input_graph[i].second.size(); j++) {
				if (input_graph[i].second[j] == e2) {
					unique_edge == false;
					break;
				}
			}
			if (unique_edge == true) {
				input_graph[i].second.insert(input_graph[i].second.end(), e2); // add e2 into adj_list of e1
			}
			break;
		}
		i++;
	}
	if (e1_is_new_vertex == true) {
		subgraph_vector_node x;
		x.first = e1;
		x.second = { e2 };
		input_graph.insert(input_graph.end(), x);  // add e1 as a new vertex, and add e2 into adj_list of e1
	}


	if (unique_edge == true) { // only continue when the edge is unique
		bool e2_is_new_vertex = true;
		int i = 0;
		while (i < input_graph.size()) {
			if (input_graph[i].first == e2) {
				e2_is_new_vertex = false;
				input_graph[i].second.insert(input_graph[i].second.end(), e1); // add e1 into adj_list of e2
				break;
			}
			i++;
		}
		if (e2_is_new_vertex == true) {
			subgraph_vector_node x;
			x.first = e2;
			x.second = { e1 };
			input_graph.insert(input_graph.end(), x);  // add e2 as a new vertex, and add e1 into adj_list of e2
		}
	}

}
#pragma endregion subgraph_vector_add_edge



#pragma region
subgraph_vector subgraph_vector_copy_subgraph_vector(subgraph_vector& input_subgraph) {

	return input_subgraph;

}
#pragma endregion subgraph_vector_copy_subgraph_vector




subgraph_vector subgraph_vector_merge_subgraph_vector(
	subgraph_vector& vector1, subgraph_vector& vector2) {

	subgraph_vector merged_vector = subgraph_vector_copy_subgraph_vector(vector1);

	for (int i = 0; i < vector2.size(); i++) {
		int v = vector2[i].first;
		std::vector<int> adj_v_vector = vector2[i].second;
		subgraph_vector_add_vertex(merged_vector, v); // copy vertex
		for (int j = 0; j < adj_v_vector.size(); j++) {
			int adj_v = adj_v_vector[j];
			if (adj_v > v) {
				subgraph_vector_add_edge(merged_vector, v, adj_v); // copy edge
			}
		}
	}

	return merged_vector;

}



#include <graph_unordered_map/graph_unordered_map.h>

graph_unordered_map subgraph_vector_2_graph_unordered_map(graph_unordered_map& input_graph,
	subgraph_vector& input_subgraph_vector) {

	/* this function changes a subgraph_vector to a graph_unordered_map;
	time complexity: O(|V_subgraph|+|E_subgraph|)*/

	graph_unordered_map output_graph;

	for (int i = 0; i < input_subgraph_vector.size(); i++) {
		int v = input_subgraph_vector[i].first;
		double v_w = input_graph[v].first;
		graph_unordered_map_add_vertex(output_graph, v, v_w); // add vertex

		std::vector<int> adj_v_vector = input_subgraph_vector[i].second;
		for (int j = 0; j < adj_v_vector.size(); j++) {
			int adj_v = adj_v_vector[j];
			if (adj_v > v) {
				double ec = input_graph[v].second[adj_v];
				graph_unordered_map_add_edge(output_graph, v, adj_v, ec); // add edge
			}
		}
	}

	return output_graph;

}




#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>

graph_hash_of_mixed_weighted subgraph_vector_2_graph_hash_of_mixed_weighted(graph_hash_of_mixed_weighted& input_graph,
	subgraph_vector& input_subgraph_vector) {

	/* this function changes a subgraph_vector to a graph_unordered_map;
	time complexity: O(|V_subgraph|+|E_subgraph|)*/

	graph_hash_of_mixed_weighted output_graph;

	for (int i = 0; i < input_subgraph_vector.size(); i++) {
		int v = input_subgraph_vector[i].first;
		double v_w = input_graph.hash_of_vectors[v].vertex_weight;
		graph_hash_of_mixed_weighted_add_vertex(output_graph, v, v_w); // add vertex

		for (int j = 0; j < input_subgraph_vector[i].second.size(); j++) {
			int adj_v = input_subgraph_vector[i].second[j];
			if (adj_v > v) {
				double ec = graph_hash_of_mixed_weighted_edge_weight(input_graph, v, adj_v);
				graph_hash_of_mixed_weighted_add_edge(output_graph, v, adj_v, ec); // add edge
			}
		}
	}

	return output_graph;

}