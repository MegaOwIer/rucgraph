#pragma once


/* define graph: a hash of vectors */

/* define graph: a hash of vectors */
class graph_hash_of_vectors_weighted_vertex_content {
public:
	double vertex_weight; // weight of this vertex
	std::vector<pair<int, double>> adj_vertices; // adjacenct vertices and weights of edges; ordered from small to large
};
typedef std::unordered_map<int, graph_hash_of_vectors_weighted_vertex_content> graph_hash_of_vectors_weighted;
/*this is an undirected, static graph*/


void graph_hash_of_vectors_weighted_add_vertex(graph_hash_of_vectors_weighted& input_graph, int vertex, double weight) {

	/*time complexity O(1)*/
	/*since unordered_map containers do not allow for duplicate keys, all the vertices are unique;
	this function can also be used to update vertex weight*/

	input_graph[vertex].vertex_weight = weight;

}


void graph_hash_of_vectors_weighted_remove_vertex(graph_hash_of_vectors_weighted& input_graph, int vertex) {

	/*this function removes a vertex and its adjacent edges;
	time complexity O(sum of degrees of adj_vertices), which is O(|E|) in the worst case*/

	/*remove vertex from adj_vertices_vectors of adj vetices*/
	auto search1 = input_graph.find(vertex);
	for (int i = 0; i < search1->second.adj_vertices.size(); i++) {
		int adj_v = search1->second.adj_vertices[i].first;
		auto search2 = input_graph.find(adj_v);
		for (int j = 0; j < search2->second.adj_vertices.size(); j++) {
			if (search2->second.adj_vertices[j].first == vertex) {
				search2->second.adj_vertices.erase(search2->second.adj_vertices.begin() + j); // remove vertex
				break; // go to the next adj_v
			}
		}
	}

	input_graph.erase(vertex); // remove vertex from input_graph

}

bool graph_hash_of_vectors_weighted_contain_vertex(graph_hash_of_vectors_weighted& input_graph, int vertex) {

	/*this function checks whether a vertex/key is in input_graph or not;
	time complexity O(1)*/

	if (input_graph.count(vertex) > 0) {
		return true;
	}
	else {
		return false;
	}

}

void graph_hash_of_vectors_weighted_add_edge(graph_hash_of_vectors_weighted& input_graph, int e1, int e2, double ec) {

	/*this function adds a weighted edge, and may add e1 and e2 into input_graph if they are new vertices;
	time complexity O(sum of degrees of e1 and e2), which is O(|V|) is the worst case;
	this can be used to update edge weight*/

	/*may add e1 or e2*/
	if (input_graph.count(e1) == 0) { // e1 is a new vertex
		graph_hash_of_vectors_weighted_add_vertex(input_graph, e1, 0); // add e1; initial weight is 0
	}
	if (input_graph.count(e2) == 0) { // e2 is a new vertex
		graph_hash_of_vectors_weighted_add_vertex(input_graph, e2, 0); // add e2; initial weight is 0
	}

	/*add e2 to adj of e1*/
	auto search = input_graph.find(e1);
	bool added = false;
	for (int i = 0; i < search->second.adj_vertices.size(); i++) {
		int adj_v = search->second.adj_vertices[i].first;
		if (adj_v == e2) {
			search->second.adj_vertices[i].second = ec; // update edge cost
			added = true;
			break; // e2 is already in adj of e1
		}
		else if (adj_v > e2) {
			search->second.adj_vertices.insert(search->second.adj_vertices.begin() + i, { e2,ec });
			// add e2 to adj of e1; make sure the adj_vector is ordered from small to large
			added = true;
			break;
		}
	}
	if (added == false) { // this happens when all adj_v are smaller than e2, i.e., e2 should be pushed back
		search->second.adj_vertices.push_back({ e2, ec }); // add e2 to adj of e1
	}

	/*add e1 to adj of e2*/
	search = input_graph.find(e2);
	added = false;
	for (int i = 0; i < search->second.adj_vertices.size(); i++) {
		int adj_v = search->second.adj_vertices[i].first;
		if (adj_v == e1) {
			search->second.adj_vertices[i].second = ec;
			added = true;
			break; // e1 is already in adj of e2
		}
		else if (adj_v > e1) {
			search->second.adj_vertices.insert(search->second.adj_vertices.begin() + i, { e1,ec });
			// add e1 to adj of e2; make sure the adj_vector is ordered from small to large
			added = true;
			break;
		}
	}
	if (added == false) { // this happens when all adj_v are smaller than e1, i.e., e1 should be pushed back
		search->second.adj_vertices.push_back({ e1, ec }); // add e1 to adj of e2
	}

}


double graph_hash_of_vectors_weighted_edge_weight(graph_hash_of_vectors_weighted& input_graph, int v1, int v2) {

	/*this function return an edge weight; we assume that this edge is in the graph, otherwise large value is returned;
	time complexity O(the smallest degree of v1 and v2), which is O(|V|) is the worst case*/

	auto search1 = input_graph.find(v1);
	auto search2 = input_graph.find(v2);

	if (search1->second.adj_vertices.size() < search2->second.adj_vertices.size()) {
		for (int i = 0; i < search1->second.adj_vertices.size(); i++) {
			int adj_v = search1->second.adj_vertices[i].first;
			if (adj_v == v2) {
				return search1->second.adj_vertices[i].second;
			}
			else if (adj_v > v2) {
				return INT_MAX;
			}
		}
	}
	else {
		for (int i = 0; i < search2->second.adj_vertices.size(); i++) {
			int adj_v = search2->second.adj_vertices[i].first;
			if (adj_v == v1) {
				return search2->second.adj_vertices[i].second;
			}
			else if (adj_v > v1) {
				return INT_MAX;
			}
		}
	}

	return INT_MAX;

}



bool graph_hash_of_vectors_weighted_contain_edge(graph_hash_of_vectors_weighted& input_graph, int v1, int v2) {

	/*this function checks whether an edge is in input_graph or not;
	time complexity O(the smallest degree of v1 and v2), which is O(|V|) is the worst case*/

	if (input_graph.count(v1) == 0 || input_graph.count(v2) == 0) {
		return false;
	}

	auto search1 = input_graph.find(v1);
	auto search2 = input_graph.find(v2);

	if (search1->second.adj_vertices.size() < search2->second.adj_vertices.size()) {
		for (int i = 0; i < search1->second.adj_vertices.size(); i++) {
			int adj_v = search1->second.adj_vertices[i].first;
			if (adj_v == v2) {
				return true;
			}
			else if (adj_v > v2) {
				return false;
			}
		}
	}
	else {
		for (int i = 0; i < search2->second.adj_vertices.size(); i++) {
			int adj_v = search2->second.adj_vertices[i].first;
			if (adj_v == v1) {
				return true;
			}
			else if (adj_v > v1) {
				return false;
			}
		}
	}

	return false;

}

void graph_hash_of_vectors_weighted_remove_edge_and_isolated_vertices(graph_hash_of_vectors_weighted& input_graph, int e1, int e2) {

	/*this function removes an edge, and may remove e1 and e2 if they are isolated;
	time complexity O(sum of degrees of e1 and e2), which is O(|V|) is the worst case*/

	auto search = input_graph.find(e1);
	for (int i = 0; i < search->second.adj_vertices.size(); i++) {
		if (search->second.adj_vertices[i].first == e2) {
			search->second.adj_vertices.erase(search->second.adj_vertices.begin() + i); // remove e2 from adj of e1
			break;
		}
	}
	if (search->second.adj_vertices.size() == 0) { // e1 is isolated
		graph_hash_of_vectors_weighted_remove_vertex(input_graph, e1); // remove e1
	}

	search = input_graph.find(e2);
	for (int i = 0; i < search->second.adj_vertices.size(); i++) {
		if (search->second.adj_vertices[i].first == e1) {
			search->second.adj_vertices.erase(search->second.adj_vertices.begin() + i); // remove e1 from adj of e2
			break;
		}
	}
	if (search->second.adj_vertices.size() == 0) { // e2 is isolated
		graph_hash_of_vectors_weighted_remove_vertex(input_graph, e2); // remove e2
	}

}

void graph_hash_of_vectors_weighted_remove_edge_but_not_isolated_vertices(graph_hash_of_vectors_weighted& input_graph, int e1, int e2) {

	/*this function removes an edge, but not remove e1 and e2 if they are isolated;
	time complexity O(sum of degrees of e1 and e2), which is O(|V|) is the worst case*/

	auto search = input_graph.find(e1);
	for (int i = 0; i < search->second.adj_vertices.size(); i++) {
		if (search->second.adj_vertices[i].first == e2) {
			search->second.adj_vertices.erase(search->second.adj_vertices.begin() + i); // remove e2 from adj of e1
			break;
		}
	}

	search = input_graph.find(e2);
	for (int i = 0; i < search->second.adj_vertices.size(); i++) {
		if (search->second.adj_vertices[i].first == e1) {
			search->second.adj_vertices.erase(search->second.adj_vertices.begin() + i); // remove e1 from adj of e2
			break;
		}
	}

}


graph_hash_of_vectors_weighted graph_hash_of_vectors_weighted_copy_graph(graph_hash_of_vectors_weighted& input_graph) {

	return input_graph;

}


int graph_hash_of_vectors_weighted_num_vertices(graph_hash_of_vectors_weighted& input_graph) {

	/*time complexity O(1)*/

	return input_graph.size();

}

long int graph_hash_of_vectors_weighted_num_edges(graph_hash_of_vectors_weighted& input_graph) {

	/*time complexity O(|V|)*/

	long int num = 0;
	for (auto it = input_graph.begin(); it != input_graph.end(); ++it) {
		num = num + it->second.adj_vertices.size();
	}

	return num / 2;

}

void graph_hash_of_vectors_weighted_print(graph_hash_of_vectors_weighted& input_graph) {

	/*time complexity O(|V|+|E|)*/

	int print_relex_num = 100;
	int num = 0;
	//std::cout << '\n';
	std::cout << "graph_hash_of_vectors_weighted_print: " << std::endl;
	std::cout << "|V|=" << graph_hash_of_vectors_weighted_num_vertices(input_graph) <<
		" |E|=" << graph_hash_of_vectors_weighted_num_edges(input_graph) << std::endl;

	for (auto it = input_graph.begin(); it != input_graph.end(); ++it) {
		std::cout << "Vertex: " << it->first << '\n';
		std::cout << "vertex_weight: " << it->second.vertex_weight << '\n';
		std::cout << "adj_vertices: ";
		for (int i = 0; i < it->second.adj_vertices.size(); i++) {
			std::cout << "<" << it->second.adj_vertices[i].first << "," << it->second.adj_vertices[i].second << "> ";
		}
		std::cout << '\n';

		num++;
		if (num%print_relex_num == 0) {
			getchar();
		}
	}

	std::cout << "graph_hash_of_vectors_weighted_print END" << '\n';
	std::cout << '\n';
}




void graph_hash_of_vectors_weighted_test() {

	graph_hash_of_vectors_weighted g;

	/*random*/
	//while (graph_hash_of_vectors_weighted_num_vertices(g) < 5) {
	//	graph_hash_of_vectors_weighted_add_vertex(g, rand() % 10, rand() % 100);
	//}
	//while (graph_hash_of_vectors_weighted_num_vertices(g) < 10 || graph_hash_of_vectors_weighted_num_edges(g) < 20) {
	//	graph_hash_of_vectors_weighted_add_edge(g, rand() % 10, rand() % 10, rand() % 100);
	//}
	//graph_hash_of_vectors_weighted_print(g);


	g.clear();
	graph_hash_of_vectors_weighted_add_vertex(g, 0, 1);
	graph_hash_of_vectors_weighted_add_vertex(g, 1, 1);
	graph_hash_of_vectors_weighted_add_vertex(g, 2, 1);

	graph_hash_of_vectors_weighted_remove_vertex(g, 0);
	cout << "graph_hash_of_vectors_weighted_contain_vertex(g,0): " << graph_hash_of_vectors_weighted_contain_vertex(g, 0) << endl;
	cout << "graph_hash_of_vectors_weighted_contain_vertex(g,1): " << graph_hash_of_vectors_weighted_contain_vertex(g, 1) << endl;



	graph_hash_of_vectors_weighted_add_edge(g, 0, 2, 1);
	graph_hash_of_vectors_weighted_add_edge(g, 1, 2, 1);
	graph_hash_of_vectors_weighted_add_edge(g, 1, 2, 2);

	cout << "graph_hash_of_vectors_weighted_contain_edge(g,1,2): " << graph_hash_of_vectors_weighted_contain_edge(g, 1, 2) << endl;
	cout << "graph_hash_of_vectors_weighted_contain_edge(g,1,0): " << graph_hash_of_vectors_weighted_contain_edge(g, 1, 0) << endl;

	cout << "graph_hash_of_vectors_weighted_edge_weight(g,1,2): " << graph_hash_of_vectors_weighted_edge_weight(g, 1, 2) << endl;
	cout << "graph_hash_of_vectors_weighted_edge_weight(g,1,0): " << graph_hash_of_vectors_weighted_edge_weight(g, 1, 0) << endl;


	graph_hash_of_vectors_weighted_remove_edge_and_isolated_vertices(g, 1, 2);
	graph_hash_of_vectors_weighted_remove_edge_but_not_isolated_vertices(g, 0, 2);



	graph_hash_of_vectors_weighted_print(g);





}