#pragma once


#include <graph_unordered_map/graph_unordered_map.h> 
#include <boost_graph/boost_graph.h> 

void graph_unordered_map_to_boost_graph(graph_unordered_map& input_graph, graph& boost_graph,
	std::unordered_map<int, int>& VID_graph_unordered_map_to_boost_graph) {

	boost_graph.clear();
	VID_graph_unordered_map_to_boost_graph.clear();

	typedef graph::edge_descriptor Edge;

	int N = input_graph.size();


	int id = 0;
	for (auto it = input_graph.begin(); it != input_graph.end(); it++) {
		VID_graph_unordered_map_to_boost_graph[it->first] = id;

		//cout << it->first << " to " << id << '\n';

		/*insert vertex*/
		//cout << "insert vertex " << id << '\n';
		boost::add_vertex(id, boost_graph);

		
		id++;
	}

	for (auto it1 = input_graph.begin(); it1 != input_graph.end(); it1++) {
		int i = it1->first;
		double i_nw = it1->second.first;

		boost::put(boost::vertex_name_t(), boost_graph, VID_graph_unordered_map_to_boost_graph[i], i_nw);

		for (auto it2 = it1->second.second.begin(); it2 != it1->second.second.end(); it2++) {
			int j = it2->first;
			double e_ij = it2->second;

			/*insert edge*/
			pair<Edge, bool> ed = boost::edge(VID_graph_unordered_map_to_boost_graph[i],
				VID_graph_unordered_map_to_boost_graph[j], boost_graph);
			if (!ed.second) { // this edge does not exit
				//cout << "insert edge " << VID_graph_unordered_map_to_boost_graph[i] <<
				//	"," << VID_graph_unordered_map_to_boost_graph[j] << '\n';
				boost::add_edge(VID_graph_unordered_map_to_boost_graph[i], VID_graph_unordered_map_to_boost_graph[j],
					e_ij, boost_graph);
			}

		}
	}
}
