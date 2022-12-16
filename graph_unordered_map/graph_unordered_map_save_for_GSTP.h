#pragma once


#include <graph_unordered_map.h> 

void graph_unordered_map_save_for_GSTP(std::string instance_name,
	graph_unordered_map& input_graph, graph_unordered_map& group_graph, 
	std::unordered_set<int>& group_vertices, double& lambda) {

	std::ofstream outputFile;
	outputFile.precision(10);
	outputFile.setf(std::ios::fixed);
	outputFile.setf(std::ios::showpoint);
	outputFile.open(instance_name);

	// comments
	outputFile << "SECTION Comments" << std::endl;
	outputFile << "Name \"" << instance_name << "\"" << std::endl;
	outputFile << "Creator \"graph_unordered_map_save_for_GSTP\"" << std::endl;
	outputFile << "END" << std::endl;
	outputFile << std::endl;

	// lambda
	outputFile << "lambda " << lambda << '\n';
	outputFile << std::endl;

	// input_graph
	outputFile << "input_graph |V|=" << graph_unordered_map_num_vertices(input_graph) <<
		" |E|=" << graph_unordered_map_num_edges(input_graph) << std::endl;
	for (auto it = input_graph.begin(); it != input_graph.end(); ++it) {
		outputFile << "input_graph Vertex " << it->first << " " << it->second.first << '\n';
	}
	outputFile << std::endl;
	for (auto it = input_graph.begin(); it != input_graph.end(); ++it) {
		for (auto it2 = it->second.second.begin(); it2 != it->second.second.end(); ++it2) {
			if (it->first <= it2->first) {
				outputFile << "input_graph Edge " << it->first << " " << it2->first << " " << it2->second << '\n';
			}
		}
	}
	outputFile << std::endl;

	// group_graph
	outputFile << "group_graph |V|=" << graph_unordered_map_num_vertices(group_graph) <<
		" |E|=" << graph_unordered_map_num_edges(group_graph) << std::endl;
	for (auto it = group_graph.begin(); it != group_graph.end(); ++it) {
		outputFile << "group_graph Vertex " << it->first << " " << it->second.first << '\n';
	}
	outputFile << std::endl;
	for (auto it = group_graph.begin(); it != group_graph.end(); ++it) {
		for (auto it2 = it->second.second.begin(); it2 != it->second.second.end(); ++it2) {
			if (it->first <= it2->first) {
				outputFile << "group_graph Edge " << it->first << " " << it2->first << " " << it2->second << '\n';
			}
		}
	}
	outputFile << std::endl;

	// group_vertices
	outputFile << "|T|=" << group_vertices.size() << std::endl;
	for (auto it = group_vertices.begin(); it != group_vertices.end(); ++it) {
		outputFile << "group_vertices " << *it << '\n';
	}
	outputFile << std::endl;

	outputFile << "EOF" << std::endl;

}