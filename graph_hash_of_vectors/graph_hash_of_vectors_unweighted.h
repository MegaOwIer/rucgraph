#pragma once


/* define graph: a hash of vectors */
typedef std::unordered_map<int, std::vector<int>> graph_hash_of_vectors_unweighted;


void graph_hash_of_vectors_unweighted_add_vertex(graph_hash_of_vectors_unweighted& input_graph, int vertex) {

	/*time complexity O(1)*/
	/*since unordered_map containers do not allow for duplicate keys,
	all the vertices in subgraph are unique*/

	if (input_graph.count(vertex) == 0) { // vertex is not in input_graph
		input_graph[vertex] = {};
	}

}


void graph_hash_of_vectors_unweighted_add_edge(graph_hash_of_vectors_unweighted& input_graph, int e1, int e2) {

	/*this function adds an unweighted edge, and may add e1 and e2 into input_graph if they are new vertices;
	time complexity O(sum of degrees of e1 and e2), which is O(|V|) is the worst case*/

	/*may add e1 or e2*/
	if (input_graph.count(e1) == 0) { // e1 is a new vertex
		graph_hash_of_vectors_unweighted_add_vertex(input_graph, e1); // add e1; initial weight is 0
	}
	if (input_graph.count(e2) == 0) { // e2 is a new vertex
		graph_hash_of_vectors_unweighted_add_vertex(input_graph, e2); // add e2; initial weight is 0
	}

	/*add e2 to adj of e1*/
	auto search = input_graph.find(e1);
	bool added = false;
	for (int i = 0; i < search->second.size(); i++) {
		int adj_v = search->second[i];
		if (adj_v == e2) {
			added = true;
			break; // e2 is already in adj of e1
		}
		else if (adj_v > e2) {
			search->second.insert(search->second.begin() + i, e2);
			// add e2 to adj of e1; make sure the adj_vector is ordered from small to large
			added = true;
			break;
		}
	}
	if (added == false) { // this happens when all adj_v are smaller than e2, i.e., e2 should be pushed back
		search->second.push_back(e2); // add e2 to adj of e1
	}

	/*add e1 to adj of e2*/
	search = input_graph.find(e2);
	added = false;
	for (int i = 0; i < search->second.size(); i++) {
		int adj_v = search->second[i];
		if (adj_v == e1) {
			added = true;
			break; // e1 is already in adj of e2
		}
		else if (adj_v > e1) {
			search->second.insert(search->second.begin() + i, e1);
			// add e1 to adj of e2; make sure the adj_vector is ordered from small to large
			added = true;
			break;
		}
	}
	if (added == false) { // this happens when all adj_v are smaller than e1, i.e., e1 should be pushed back
		search->second.push_back(e1); // add e1 to adj of e2
	}

}