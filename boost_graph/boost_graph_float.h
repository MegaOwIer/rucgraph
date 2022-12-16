#pragma once


#include <boost/graph/adjacency_list.hpp>

#pragma region 
// define an adjacency list with edge weights
typedef boost::property<boost::edge_weight_t, float> EdgeWeightProperty_float; // define edge weight property
typedef boost::property<boost::vertex_name_t, float> VertexWeightProperty_float; // define node weight property; note that: vertex_index_t is not mutable
typedef boost::adjacency_list<boost::setS, boost::vecS,
	boost::undirectedS, VertexWeightProperty_float, EdgeWeightProperty_float> graph_float; // define all the graph properties
typedef boost::graph_traits<graph_float>::adjacency_iterator AdjacencyIterator_float;
#pragma endregion define graph_float_weight property 


#pragma region
graph_float copy_graph_float(graph_float& input_graph) {

	return input_graph; // this input_graph is just a copy of the outside input_graph

}
#pragma endregion copy_graph 2018_5_9_17:39:01
