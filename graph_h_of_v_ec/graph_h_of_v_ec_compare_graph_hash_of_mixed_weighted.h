#pragma once


/*

1. in dense graphs, graph_h_of_v_ec consumes less memory than graph_hash_of_mixed_weighted (but speed is similar)
2. in sparse graphs, graph_h_of_v_ec consumes similar memory than graph_hash_of_mixed_weighted (speed is similar)

e.g.:

graph_h_of_v_ec_print_size: |V|=50000 |E|=500000
graph_h_of_v_ec_generate_random_graph    runningtime: 9.82479s
graph_h_of_v_ec_generate_random_graph    RAM: 35.9727 MB

graph_hash_of_mixed_weighted_print_size: |V|=50000 |E|=500000
graph_hash_of_mixed_weighted_generate_random_graph    runningtime: 9.80442s
graph_hash_of_mixed_weighted_generate_random_graph    RAM: 36 MB

graph_h_of_v_ec_print_size: |V|=50000 |E|=5000000
graph_h_of_v_ec_generate_random_graph    runningtime: 100.899s
graph_h_of_v_ec_generate_random_graph    RAM: 231.965 MB

graph_hash_of_mixed_weighted_print_size: |V|=50000 |E|=5000000
graph_hash_of_mixed_weighted_generate_random_graph    runningtime: 108.39s
graph_hash_of_mixed_weighted_generate_random_graph    RAM: 1059.09 MB








*/


#include <graph_h_of_v_ec/graph_h_of_v_ec.h> 
#include <graph_h_of_v_ec/graph_h_of_v_ec_generate_random_graph.h> 
#include <Current_Memory_Consumption_of_This_Process.h> 
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h> 
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_generate_random_graph.h> 

void graph_h_of_v_ec_compare_graph_hash_of_mixed_weighted() {

	int V = 5e4, E = 5e5;

	/*graph_h_of_v_ec*/
	if (1) {
		auto begin = std::chrono::high_resolution_clock::now();
		graph_h_of_v_ec g = graph_h_of_v_ec_generate_random_graph(V, E, 1, 10, 2);
		auto end = std::chrono::high_resolution_clock::now();
		double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
		graph_h_of_v_ec_print_size(g);
		cout << "graph_h_of_v_ec_generate_random_graph    runningtime: " << runningtime << "s" << endl;
		cout << "graph_h_of_v_ec_generate_random_graph    RAM: " << Current_Memory_Consumption_of_This_Process() << " MB" << endl;
		cout << endl;
	}

	/*graph_hash_of_mixed_weighted*/
	if (1) {
		auto begin = std::chrono::high_resolution_clock::now();
		graph_hash_of_mixed_weighted g = graph_hash_of_mixed_weighted_generate_random_graph(V, E, 1, 10, 1, 10, 2);
		auto end = std::chrono::high_resolution_clock::now();
		double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
		graph_hash_of_mixed_weighted_print_size(g);
		cout << "graph_hash_of_mixed_weighted_generate_random_graph    runningtime: " << runningtime << "s" << endl;
		cout << "graph_hash_of_mixed_weighted_generate_random_graph    RAM: " << Current_Memory_Consumption_of_This_Process() << " MB" << endl;
		cout << endl;
	}


}

