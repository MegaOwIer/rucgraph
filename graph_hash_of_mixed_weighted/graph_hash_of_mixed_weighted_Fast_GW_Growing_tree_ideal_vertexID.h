#pragma once

/*  This is the fast implementation of the GW growing process for solving the prize-collecting Steiner
	tree problem. Details in Sun, Yahui, et al. "The fast heuristic algorithms and post-processing techniques to design
	large and low-cost communication networks." IEEE/ACM Transactions on Networking 27.1 (2019): 375-388.

	Time complexity: O(V + V log V + d E log V); In 2014 C hedge paper, it's claimed that the time compleity
	of this growing process is O(d E log V), which is incorrect, as all vertices may have positive weights, i.e.,
	we need to initialize V clusters in general cases.

	We assume that vertex IDs are 0,1,...,N-1 in this function; otherwise you need to
	change vertex ID before calling this function; without this assumption, the growing process may consume
	much more RAM, as hashs must be used.

	There is no root in this function;

	The number of vertices is limited to INT_MAX in this function;

	If return_content = "tree", then this function only returns the last active cluster/tree, but not inactive clusters;
	if return_content = "forest", then this function returns both active and inactive clusters;

	The input graph may be disconnected.

	We assume that all weight values are smaller than std::numeric_limits<float>::max(), i.e., 3.40282e+38, otherwise
	there may be errors. The reason is that an active cluster C in a disconnected graph may be isolated, in which case
	both two heaps C_event_PQ and E_event_PQ contain C, and E_event_PQ[C] is empty, and C_ep_PQ[C].get_top()
	returns <std::numeric_limits<float>::max(),std::numeric_limits<float>::max()>, which makes
	E_event_PQ[C] = <std::numeric_limits<float>::max(),C> (see our #include <PairingHeapYS.h> and
	#include <PairingHeapYS_with_offset.h>)

	if there is no vertex with positive weight, then this function returns empty

*/

#pragma region
#include <Union_Find.h>
#include <PairingHeapYS.h>
typedef PairingHeapYS_min<double, int> PairingHeap_3; // define data types in heap
#include <PairingHeapYS_with_offset.h>
typedef PairingHeapYS_with_offset<double, int> PairingHeap_4; // define data types in heap

graph_hash_of_mixed_weighted graph_hash_of_mixed_weighted_Fast_GW_Growing_tree_ideal_vertexID
(graph_hash_of_mixed_weighted& input_graph, string return_content) {

	if (return_content != "tree" && return_content != "forest") {
		std::cout << "Error in graph_hash_of_mixed_weighted_Fast_GW_Growing_idealvertexID: " <<
			"return_content != tree && return_content != forest" << std::endl;
		exit(1);
	}

	/*we assume that distribution_ratio = 2, i.e., the slacks of edges are divide evenly at the initial stage;
	details in Sun, Yahui, et al. "The fast heuristic algorithms and post-processing techniques to design
	large and low-cost communication networks." IEEE/ACM Transactions on Networking 27.1 (2019): 375-388.*/
	double distribution_ratio = 2;


	/*initialization; time complexity O(|V|) (to calculate number of edges)*/
	graph_hash_of_mixed_weighted solution_tree;
	double Global_time = 0; // global time
	int Active_C_num = 0; // the number of active clusters
	int N = input_graph.hash_of_vectors.size(); // number of vertices
	long long int ep_num = graph_hash_of_mixed_weighted_num_edges(input_graph) * 2; // total number of edge parts
	int ep_order = 0;

	/*Clusters: the number of clusters is always N; time complexity O(|V|)*/
	vector<bool> C_activity(N, false); // activity value of each C; false means inactive
	vector<double> C_event_time(N); // the event time for each C

	/*the use of Union_Find is suggested by C Hedge in 2014 DIMACS paper;
	vertex i is in cluster Union_Find_find(i, vertex_unions, parent);
	initially Union_Find_find(i, vertex_unions, parent) = i*/
	vector<int> vertex_unions(N), parent(N), set_size(N); int sets_num;
	Union_Find_initialization(vertex_unions, parent, set_size, sets_num);

	/*edge parts: PQ and their handles; time complexity O(|E|)*/
	std::vector<PairingHeap_4::ItemHandle> handle_ep(ep_num); // store the handle for each edge part
	vector<PairingHeap_4> C_ep_PQ(N);
	vector<int> ep_v1_list(ep_num);
	vector<int> ep_v2_list(ep_num);
	vector<double> ep_EventTime_list(ep_num);

	/*the event PQ and their handles; time complexity O(|V|)*/
	std::vector<PairingHeap_3::ItemHandle> handle_Cevent(N), handle_Eevent(N);
	PairingHeap_3 C_event_PQ, E_event_PQ;

	/*initialize the clusters; time complexity O(|V|+|E|)*/
	for (int i = 0; i < N; i++)
	{
		/*add edge parts into C*/
		auto search = input_graph.hash_of_hashs.find(i);
		if (search != input_graph.hash_of_hashs.end()) {
			for (auto it2 = search->second.begin(); it2 != search->second.end(); it2++) {
				int j = it2->first; // the adjacent vetex to i
				if (j > i) { // don't overcheck an edge
				// the first ep
					ep_v1_list[ep_order] = i;
					ep_v2_list[ep_order] = j;
					ep_EventTime_list[ep_order] = graph_hash_of_mixed_weighted_edge_weight(input_graph, i, j)
						/ distribution_ratio; // halve the edge cost
					handle_ep[ep_order] = C_ep_PQ[i].insert(ep_EventTime_list[ep_order], ep_order); // put this ep into cluster i
					ep_order++;
					// the second ep
					ep_v1_list[ep_order] = j;
					ep_v2_list[ep_order] = i;
					ep_EventTime_list[ep_order] = ep_EventTime_list[ep_order - 1]
						* (distribution_ratio - 1); // halve the edge cost
					handle_ep[ep_order] = C_ep_PQ[j].insert(ep_EventTime_list[ep_order], ep_order); // put this ep into cluster j
					ep_order++;
				}
			}
		}
		else {
			auto search2 = input_graph.hash_of_vectors.find(i);
			for (auto it2 = search2->second.adj_vertices.begin(); it2 != search2->second.adj_vertices.end(); it2++) {
				int j = it2->first; // the adjacent vetex to i
				if (j > i) { // don't overcheck an edge
				// the first ep
					ep_v1_list[ep_order] = i;
					ep_v2_list[ep_order] = j;
					ep_EventTime_list[ep_order] = graph_hash_of_mixed_weighted_edge_weight(input_graph, i, j)
						/ distribution_ratio; // halve the edge cost
					handle_ep[ep_order] = C_ep_PQ[i].insert(ep_EventTime_list[ep_order], ep_order); // put this ep into cluster i
					ep_order++;
					// the second ep
					ep_v1_list[ep_order] = j;
					ep_v2_list[ep_order] = i;
					ep_EventTime_list[ep_order] = ep_EventTime_list[ep_order - 1]
						* (distribution_ratio - 1); // halve the edge cost
					handle_ep[ep_order] = C_ep_PQ[j].insert(ep_EventTime_list[ep_order], ep_order); // put this ep into cluster j
					ep_order++;
				}
			}
		}

		/*for active cluster*/
		double nw_i = input_graph.hash_of_vectors[i].vertex_weight;
		if (nw_i > 0) {
			Active_C_num++; // the number of active clusters
			C_activity[i] = true; // this cluster is active
			C_event_time[i] = nw_i; // the event time is the node weight
			// push node into C_event_PQ
			handle_Cevent[i] = C_event_PQ.insert(C_event_time[i], i); // C_event_PQ only contain active clusters
			// all the ep for cluster i have been inserted into C_ep_PQ[i]; 
			//Note that, it's only true when i starts from 0 and j>i above
			// push node into E_event_PQ
			handle_Eevent[i] = E_event_PQ.insert(C_ep_PQ[i].get_top().first, i); // E_event_PQ only contain active clusters
		}
	}

	//cout << "Fast growth starts!" << endl;

	/*Fast growth starts!; time complexity O(V log V + d E log V)*/
	int C0, C1, C2, ep1, ep2;
	double Tc, Te, r, tiny_value = 1e-7;  // tiny_value is like the pricision of weights: d
	vector<pair<int, int>> merged_edges;
	pair<int, int> edge;
	while (Active_C_num > 1) // stop the loop when there is only 1 active clusters left
	{
		// find the closest event
		Tc = C_event_PQ.get_top().first; // this cluster event time
		Te = E_event_PQ.get_top().first; // this edge event time


		//cout << "E_event_PQ.size(): " << E_event_PQ.size() << endl;

		//cout << "Tc = " << Tc << " Te = " << Te << endl;

		if (Tc >= Te) {
			/*edge event; there are at most O(d E) edge events; time complexity O(d E log V);

			Tc >= Te makes it possible to get such a clusting tree:
			v1 (nw=4) e(v1,v2)=8 v2(nw=4) e(v2,v3)=8 v3(nw=4)
			Tc > Te will not produce such a result; however, it is preferable
			to set Tc >= Te to get as many edge events as possible,
			as in these cases the next general prunning result may be better*/

			//cout << "edge event starts!" << endl;

			/*time complexity: O(alpha E) for Union_Find_find*/
			C1 = E_event_PQ.get_top().second; // the cluster C1 for this edge event
			ep1 = C_ep_PQ[C1].get_top().second; // the top ep in C1
			C2 = Union_Find_find(ep_v2_list[ep1], vertex_unions, parent); // the cluster C2 for this edge event

			//cout << "C1 = " << C1 << " C2 = " << C2 << endl;

			if (C1 == C2) { /*the inside ep is triggered; C_ep_PQ[C1].size at most V; time complexity: O(log V)*/
				C_ep_PQ[C1].pop_top(); // pop out the inside ep
				// decrease E_event_PQ for the change of event_C1
				E_event_PQ.change_key(handle_Eevent[C1], C_ep_PQ[C1].get_top().first);
				/* in disconnected graphs, C_ep_PQ[C1] may become empty, in which case
				C_ep_PQ[C1].get_top().first is std::numeric_limits<double>::max()*/
			}
			else { /*the outside ep is triggered; time complexity: O(logV)*/
				Global_time = Te;
				/*details in the initialize the clusters (you initial ep_order by pairs from 0 to 2E-1)*/
				if (ep1 % 2 == 0) {
					ep2 = ep1 + 1;
				}
				else {
					ep2 = ep1 - 1;
				}
				if (C_activity[C2] == true) { // C2 is active; time complexity: O(logV)
					r = ep_EventTime_list[ep2] - Global_time; // the slack of the responsible edge
					if (r > tiny_value) { /*r is big; time complexity: O(logV)*/
						// change two ep event time
						ep_EventTime_list[ep1] = Global_time + r / 2;
						ep_EventTime_list[ep2] = Global_time + r / 2;
						// update C_ep_PQ in C1
						C_ep_PQ[C1].change_key(handle_ep[ep1], Global_time, ep_EventTime_list[ep1]);
						// update C_ep_PQ in C2
						C_ep_PQ[C2].change_key(handle_ep[ep2], Global_time + r, ep_EventTime_list[ep2]);
						// update E_event_PQ for the change of C1
						E_event_PQ.change_key(handle_Eevent[C1], C_ep_PQ[C1].get_top().first);
						// update E_event_PQ for the change of C2
						if (C_ep_PQ[C2].get_top().second == ep2) { // the new smaller ep_EventTime_list[ep2] is top
							E_event_PQ.change_key(handle_Eevent[C2], C_ep_PQ[C2].get_top().first);
						}
					}
					else { /*r is small; merge event*/

						//cout << "merge event between two active clusters" << endl;

						merged_edges.push_back({ ep_v1_list[ep1] ,ep_v1_list[ep2] }); // merge edge 
						// merge V_list of C2 into C1; smaller set merge into larger set; time complexity: O(1)
						bool C1_merge_into_c2, C2_merge_into_c1;
						if (set_size[Union_Find_find(C1, vertex_unions, parent)]
							< set_size[Union_Find_find(C2, vertex_unions, parent)]) {
							C1_merge_into_c2 = true; // see Union_Find_Union source code
							C2_merge_into_c1 = false;
						}
						else {
							C1_merge_into_c2 = false; // see Union_Find_Union source code
							C2_merge_into_c1 = true;
						}
						Union_Find_Union(C1, C2, vertex_unions, parent, set_size, sets_num);


						//cout << "C1_merge_into_c2 = " << C1_merge_into_c2 << endl;
						//cout << "Union_Find_find(C1, vertex_unions, parent) = "
						//	<< Union_Find_find(C1, vertex_unions, parent) << endl;
						//cout << "Union_Find_find(C2, vertex_unions, parent) = "
						//	<< Union_Find_find(C2, vertex_unions, parent) << endl;

						C_ep_PQ[C1].pop_top(); // pop out the responsible ep; time complexity: O(logV)
						Active_C_num--; // minus one active cluster
						if (C1_merge_into_c2 == true) {
							/*update event time*/
							C_event_time[C2] = C_event_time[C1] + C_event_time[C2] - Global_time;
							C_activity[C1] = false;
							C_ep_PQ[C2].meld(C_ep_PQ[C1]); // merge two C_ep_PQ
							// update C2 in C_event_time and E_event_time
							C_event_PQ.change_key(handle_Cevent[C2], C_event_time[C2]);
							E_event_PQ.change_key(handle_Eevent[C2], C_ep_PQ[C2].get_top().first);
							// remove C1 from C_event_time and E_event_time
							C_event_PQ.erase(handle_Cevent[C1]);
							E_event_PQ.erase(handle_Eevent[C1]);
						}
						else {
							/*update event time*/
							C_event_time[C1] = C_event_time[C1] + C_event_time[C2] - Global_time;
							C_activity[C2] = false;
							C_ep_PQ[C1].meld(C_ep_PQ[C2]); // merge two C_ep_PQ

							//cout << "here 1" << endl;

							// update C1 in C_event_time and E_event_time
							C_event_PQ.change_key(handle_Cevent[C1], C_event_time[C1]);

							//cout << "here 2" << endl;


							E_event_PQ.change_key(handle_Eevent[C1], C_ep_PQ[C1].get_top().first);

							//cout << "here 3" << endl;

							// remove C2 from C_event_time and E_event_time
							C_event_PQ.erase(handle_Cevent[C2]);

							E_event_PQ.erase(handle_Eevent[C2]);
						}
					}
				}
				else { // C2 is inactive; time complexity: O(logV)
					r = ep_EventTime_list[ep2] - C_event_time[C2]; // the slack of the responsible edge
					if (r > tiny_value) { // r is big; d is not used in this coding
						// change two ep event time
						ep_EventTime_list[ep1] = Global_time + r;
						ep_EventTime_list[ep2] = C_event_time[C2];
						// update C_ep_PQ in C1
						C_ep_PQ[C1].change_key(handle_ep[ep1], Global_time, ep_EventTime_list[ep1]);
						// update C_ep_PQ in C2
						C_ep_PQ[C2].change_key(handle_ep[ep2], ep_EventTime_list[ep2] + r, ep_EventTime_list[ep2]);
						// update E_event_PQ for the change of C1
						E_event_PQ.change_key(handle_Eevent[C1], C_ep_PQ[C1].get_top().first);
					}
					else { // r is small; merge event
						merged_edges.push_back({ ep_v1_list[ep1] ,ep_v1_list[ep2] }); // merge edge 

						// merge V_list of C2 into C1; time complexity: O(1)
						bool C1_merge_into_c2, C2_merge_into_c1;
						if (set_size[Union_Find_find(C1, vertex_unions, parent)]
							< set_size[Union_Find_find(C2, vertex_unions, parent)]) {
							C1_merge_into_c2 = true; // see Union_Find_Union source code
							C2_merge_into_c1 = false;
						}
						else {
							C1_merge_into_c2 = false; // see Union_Find_Union source code
							C2_merge_into_c1 = true;
						}
						Union_Find_Union(C1, C2, vertex_unions, parent, set_size, sets_num);

						C_ep_PQ[C1].pop_top(); // pop out the responsible ep; time complexity: O(logV)	

						/*merge two C_ep_PQ; 2014 C Hedge: worst case time complexity: O(logV):
						The requirements for the priority queue maintaining edge parts events
						are the standard operations of a mergeable heap data structure, combined with
						an operation that adds a constant offset to all elements in a heap in O(log n)
						amortized time. We can build such a data structure by augmenting a pairing
						heap [FSST86] with an offset value at each node.
						details in pairing_heap.h at
						https://github.com/YahuiSun/pcst-fast/tree/master/src */
						C_ep_PQ[C2].add_to_heap(Global_time - C_event_time[C2]);

						if (C1_merge_into_c2 == true) {
							C_activity[C1] = false;
							C_activity[C2] = true;
							C_event_time[C2] = C_event_time[C1];
							C_ep_PQ[C2].meld(C_ep_PQ[C1]);
							// update C2 in C_event_time and E_event_time
							delete(handle_Cevent[C2]); // delete nodes after using their handles to insert new node, to avoid memory leak
							delete(handle_Eevent[C2]);
							handle_Cevent[C2] = C_event_PQ.insert(C_event_time[C2], C2);
							handle_Eevent[C2] = E_event_PQ.insert(C_ep_PQ[C2].get_top().first, C2);
							// remove C1 from C_event_time and E_event_time
							C_event_PQ.erase(handle_Cevent[C1]);
							E_event_PQ.erase(handle_Eevent[C1]);
						}
						else {
							C_ep_PQ[C1].meld(C_ep_PQ[C2]);
							E_event_PQ.change_key(handle_Eevent[C1], C_ep_PQ[C1].get_top().first);
						}
					}
				}
			}
		}
		else { // cluster event; time complexity: O(V logV)

		//cout << "cluster event starts!" << endl;

			Global_time = Tc; // decrease time
			C0 = C_event_PQ.get_top().second; // the cluster for this cluster event

			//cout << "C0 = " << C0 << endl;

			Active_C_num--; // minus one active cluster
			C_event_PQ.pop_top(); // remove the cluster from C_event_PQ

			//cout << "C_event_PQ.size() = " << C_event_PQ.size() << endl;
			//E_event_PQ.print(E_event_PQ.root);

			E_event_PQ.erase(handle_Eevent[C0]); // remove the cluster from E_event_PQ

			//cout << "E_event_PQ.size() = " << E_event_PQ.size() << endl;

			C_activity[C0] = false; // deactivate it
		}
	}


	/*if there is no vertex with positive weight, then this function returns empty*/
	if (return_content == "tree") {
		// output solution_tree; time complexity O(V)
		for (int i = 0; i < N; i++) {
			if (C_activity[i] == true) {
				double nw1 = input_graph.hash_of_vectors[i].vertex_weight;
				graph_hash_of_mixed_weighted_add_vertex(solution_tree, i, nw1); // solution may only be this single vertex
			}
		}
		for (int i = 0; i < merged_edges.size(); i++) {
			if (C_activity[Union_Find_find(merged_edges[i].first, vertex_unions, parent)] == true) { // an edge in active cluster
				int v1 = merged_edges[i].first;
				int v2 = merged_edges[i].second;
				double nw1 = input_graph.hash_of_vectors[v1].vertex_weight;
				double nw2 = input_graph.hash_of_vectors[v2].vertex_weight;
				double ec = graph_hash_of_mixed_weighted_edge_weight(input_graph, v1, v2);
				graph_hash_of_mixed_weighted_add_vertex(solution_tree, v1, nw1);
				graph_hash_of_mixed_weighted_add_vertex(solution_tree, v2, nw2);
				graph_hash_of_mixed_weighted_add_edge(solution_tree, v1, v2, ec);
			}
		}
	}
	else {
		// output solution_tree; time complexity O(V)
		for (int i = 0; i < N; i++) {
			double nw1 = input_graph.hash_of_vectors[i].vertex_weight;
			if (nw1 > 0) {
				graph_hash_of_mixed_weighted_add_vertex(solution_tree, i, nw1); // insert all vertices with positive weights
			}
		}
		for (int i = 0; i < merged_edges.size(); i++) {
			int v1 = merged_edges[i].first;
			int v2 = merged_edges[i].second;
			double nw1 = input_graph.hash_of_vectors[v1].vertex_weight;
			double nw2 = input_graph.hash_of_vectors[v2].vertex_weight;
			double ec = graph_hash_of_mixed_weighted_edge_weight(input_graph, v1, v2);
			graph_hash_of_mixed_weighted_add_vertex(solution_tree, v1, nw1);
			graph_hash_of_mixed_weighted_add_vertex(solution_tree, v2, nw2);
			graph_hash_of_mixed_weighted_add_edge(solution_tree, v1, v2, ec);
		}
	}

	/*release memory of priority queues*/
	for (int i = handle_Cevent.size() - 1; i >= 0; i--) {
		delete(handle_Cevent[i]);
	}
	for (int i = handle_Eevent.size() - 1; i >= 0; i--) {
		delete(handle_Eevent[i]);
	}
	for (int i = handle_ep.size() - 1; i >= 0; i--) {
		delete(handle_ep[i]);
	}



	return solution_tree;
}
#pragma endregion graph_hash_of_mixed_weighted_Fast_GW_Growing_tree_ideal_vertexID 

