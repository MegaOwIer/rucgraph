#pragma once

#include <limits.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stack>
#include <list>
#include <map>
#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <utility>
#include <fstream>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_General_Pruning_forest.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_minimum_spanning_tree.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_Fast_GW_Growing_tree.h>
#include <boost/random.hpp>
using namespace std;

/*___________________TGraph____________________*/
// Temporal graph
class TGraph {
public:
    int VNum;
    int ENum;
    int T;
    //std::vector <int> e; //store all edges(edgeindex->endnode) in gragh

    std::vector <std::vector <int>> w; //edge weights

    std::unordered_map <int, std::vector<std::pair<int, int>>> hash;//hash of vectors store the index of edge

public:

    TGraph() { VNum = 0; ENum = 0; T = 0; };

    void insert_edge(int u, int v);// insert a edge, u is the start node, v is the end node, input_w is the edge weight sequence varies by time
    void insert_time(std::vector <int>& edgeweights);// insert edge weight in times

    void edge_sort();

    int vertex_degree(int u);
    int getMaxOutDegree(); // return maxdegree of all edges in e
    bool contain_edge(int u, int v);
    void getEW_inT(std::vector<int>& ew, int t); // return edge weight in time t
    void getSumEW_in_lr(std::vector<int>& ew, int l, int r); // return sum edgeweight in [l. r]
    void getMaxEW_in_lr(std::vector<int>& ew, int l, int k1, int k2); // return max edgeweight in {i, k1, k2} 

    //void graph_to_file(std::string map_file, std::string time_file);
    //void build_from_file(std::string map_file, std::string time_file, int begin_time_l, int end_time_r); //read input from file named file_name, and build temporal graph from file, input format as follow:
    void generate_tgraph_from_two_vectors(std::vector <std::pair <int, int>>& graph_data, std::vector <std::vector <int>>& temporal_data, int begin_time_l, int end_time_r);
    //void print_edge();//print all edge in graph
};

// generate tgraph from graph_data and temporal_data
void TGraph::generate_tgraph_from_two_vectors(std::vector <std::pair <int, int>>& graph_data, std::vector <std::vector <int>>& temporal_data, int begin_time_l, int end_time_r) {
    // find the max vertex number
    int MaxVnum = 0;
    for (auto iter = graph_data.begin(); iter != graph_data.end(); iter++) {
        if (iter->first > MaxVnum) MaxVnum = iter->first;
        if (iter->second > MaxVnum) MaxVnum = iter->second;
        insert_edge(iter->first, iter->second);
    }
    // for vertex from 0-MaxVnum, if its degree=0(i.e. vertex not include in the graph_data)
    for (int i = 0; i < MaxVnum + 1; i++) {
        if (hash.find(i) == hash.end()) hash[i].resize(0);
    }
    VNum = hash.size();
    // T = temporal_data.size();
    T = end_time_r - begin_time_l + 1;
    // compute temporal weight of each edge
    for (int i = 0; i < ENum; i++) {
        std::vector<int> temp_vector(T, -1);
        w.emplace_back(temp_vector);
    }

    // input temporal data
    for (auto iter1 = temporal_data.begin() + begin_time_l; iter1 != temporal_data.end() && iter1 - temporal_data.begin() <= end_time_r; iter1++) {
        int temp_time = iter1 - temporal_data.begin() - begin_time_l;
        for (auto iter2 = iter1->begin(); iter2 != iter1->end(); iter2++) {
            w[*iter2][temp_time] = 1; // *iter2 超过w的size的话说明出入的边里面有重边；因为下面TGraph::insert_edge保证w里面没有重边
        }
    }

    // compute sumprex
    for (auto iter1 = w.begin(); iter1 != w.end(); iter1++) {
        for (auto iter2 = iter1->begin() + 1; iter2 != iter1->end(); iter2++) {
            *iter2 += *(iter2 - 1);
        }
    }
    edge_sort();
}


void TGraph::insert_edge(int u, int v)
{
    if (contain_edge(u, v)) return;
    // for undirected map, (u, v) and (v, u) share the same edge index
    hash[u].push_back(std::make_pair(v, ENum));
    hash[v].push_back(std::make_pair(u, ENum));
    ENum++;
}

// function used for random graph
void TGraph::insert_time(std::vector <int>& edgeweights) {
    if (edgeweights.size() == 0) w.push_back(edgeweights);
    else {
        w.push_back(edgeweights);
        int T_size = edgeweights.size();
        int edge_index = w.size() - 1;
        for (int temp_t = 1; temp_t < T_size; temp_t++) {
            // int edge_index = it - w[w_size].begin();
            w[edge_index][temp_t] += w[edge_index][temp_t - 1];
        }
    }
}


int TGraph::vertex_degree(int u) {
    return hash[u].size();
}

int TGraph::getMaxOutDegree() {
    int u = 0;
    int nodenum = hash.size();
    for (int i = 0; i < nodenum; i++) {
        if (hash[i].size() > hash[u].size()) {
            u = i;
        }
    }
    return u;
}

bool TGraph::contain_edge(int u, int v) {
    bool bl = false;
    if (hash.find(u) != hash.end()) {
        for (int i = 0; i < hash[u].size(); i++) {
            if (hash[u][i].first == v) {
                bl = true;
                break;
            }
        }
    }
    return bl;
}

void TGraph::getEW_inT(std::vector<int>& ew, int t) {

    if (t == 0) {
        for (int i = 0; i < ENum; i++) {
            ew[i] = w[i][t];
        }
    }
    else if (t > 0){
        for (int i = 0; i < ENum; i++) {
            ew[i] = w[i][t] - w[i][t - 1];
        }
    }
}

void TGraph::getSumEW_in_lr(std::vector<int>& ew, int l, int r) {
    if (l == 0) {
        for (int i = 0; i < ENum; i++) {
            ew[i] = w[i][r];
        }
    }
    else {
        for (int i = 0; i < ENum; i++) {
            ew[i] = w[i][r] - w[i][l - 1];
        }
    }
}

void TGraph::getMaxEW_in_lr(std::vector<int>& ew, int l, int k1, int k2) {;
    for (int i = 0; i < ENum; i++) {
        for (int j = k1; j <= k2; j++) {
            if (w[i][j] > ew[i]) ew[i] = w[i][j];
        }
    }
    if (l != 0) {
        for (int i = 0; i < ENum; i++) {
            ew[i] = ew[i] - w[i][l - 1];
        }
    }
}

void TGraph::edge_sort()
{
    for (std::unordered_map <int, std::vector<std::pair<int, int>>>::iterator it = hash.begin(); it != hash.end(); it++) {
        std::sort(it->second.begin(), it->second.end());
    }
}

//void TGraph::graph_to_file(std::string map_file, std::string time_file) {
//    FILE* fp;
//    fopen_s(&fp, map_file.c_str(), "w");
//    fprintf_s(fp, "%d %d\n", VNum, ENum);
//    for (auto it = hash.begin(); it != hash.end(); it++) {
//        for (int i = 0; i < it->second.size(); i++) {
//            fprintf_s(fp, "%d %d\n", it->first, it->second[i].first);
//        }
//    }
//    fclose(fp);
//
//    FILE* fpt;
//    fopen_s(&fpt, time_file.c_str(), "w");
//    fprintf_s(fp, "%d\n", T);
//    for (auto it = w.begin(); it != w.end(); it++) {
//        fprintf_s(fpt, "%d ", it->size());
//        for (auto it2 = it->begin(); it2 != it->end(); it2++) {
//            fprintf_s(fpt, "%d ", *it2);
//        }
//        fprintf_s(fpt, "\n");
//    }
//    fclose(fpt);
//    std::cout << "map & time files complete" << std::endl;
//}
//
//void TGraph::build_from_file(std::string map_file, std::string time_file, int begin_time_l, int end_time_r)
//{
//    hash.clear(); w.clear();
//
//    FILE* fp;
//    if (fopen_s(&fp, map_file.c_str(), "r")) {
//        std::cout << "File doesn't exist!\n";
//        exit(1);
//    }
//    int vertex_num, edge_num;
//    std::vector <std::pair <int, int>> graph_data;
//    std::vector <std::vector <int>> temporal_data;
//    fscanf_s(fp, "%d %d", &vertex_num, &edge_num);
//    for (int i = 0; i < edge_num; i++) {
//        int u, v;
//        fscanf_s(fp, "%d %d", &u, &v);
//        graph_data.push_back(std::pair<int, int>{u, v});        
//    }
//    fclose(fp);
//
//    FILE* fpt;
//    if (fopen_s(&fpt, time_file.c_str(), "r")) {
//        std::cout << "File doesn't exist!\n";
//        exit(1);
//    }
//    int timelen;
//    fscanf_s(fpt, "%d", &timelen);
//    for (int i = 0; i < timelen; i++) {
//        int t;
//        fscanf_s(fpt, "%d", &t);
//        std::vector<int> tempew;
//        for (int j = 0; j < t; j++) {
//            int u;
//            fscanf_s(fpt, "%d", &u);
//            tempew.push_back(u);
//        }
//        temporal_data.push_back(tempew);
//    }
//    fclose(fpt);
//    
//    generate_tgraph_from_two_vectors(graph_data, temporal_data, begin_time_l, end_time_r);
//    std::cout << "build complete" << std::endl;
//}
//
//void TGraph::print_edge()
//{
//    std::cout << "VNum: " << VNum << "\t ENum: " << ENum << "\t T: " << T << "\n";
//    for (auto it = hash.begin(); it != hash.end(); it++) {
//        std::cout << "vertex " << it->first << " has edge:" << std::endl;
//        for (int i = 0; i < it->second.size(); i++) {
//            int temp_e = it->second[i].second;
//            std::cout << "edge index:" << temp_e << "\tendpoint at: " << it->second[i].first /*<< " has value: "*/;
//            /*for (int j = 0; j < w[temp_e].size(); j++) {
//                std::cout << w[temp_e][j] << " ";
//            }*/
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << "Edge Weight of Times: " << std::endl;
//    for (auto it = w.begin(); it != w.end(); it++) {
//        std::cout << "edge " << it - w.begin() << ": ";
//        for (auto it2 = it->begin(); it2 != it->end(); it2++) {
//            std::cout << *it2 << " ";
//        }
//        std::cout << std::endl;
//    }
//}

/*_______________________________________*/

/*______________EGraph___________________*/

class EGraph {
// edge-weighted graph
public:
    int VNum;
    int ENum;

    std::vector <int> ew; //edge weights

    std::unordered_map <int, std::vector<std::pair<int, int>>>* hashgraph;//hash of vectors store the index of edge

public:
    EGraph() { VNum = 0; ENum = 0; };
    EGraph(TGraph* tgraph) {
        VNum = tgraph->VNum;
        ENum = tgraph->ENum;
        hashgraph = &(tgraph->hash);
        // ew.resize(ENum, 0);
    };
    /*~EGraph() {
        ew.~vector();
    }*/

    void addEdge(int src, int dest, int weight) {
        int edge_cnt = ew.size();
        (*hashgraph)[src].push_back(std::make_pair(dest, edge_cnt));
        (*hashgraph)[dest].push_back(std::make_pair(src, edge_cnt));
        ew.push_back(weight);
        ENum++;
    };
    int computeUBsop(); // compute the sum of positive edges in graph
    int computeUBstr(); // compute the UBstr in paper
    int uf_find(std::vector<int>& parent, int x); // union find - find
    void uf_union(std::vector<int>& parent, std::vector<int>& rank, int x, int y); // union find - union
    int merge(std::vector<int>& parent, int x);
    int Connected_Components(std::vector<int>& nodes); // find connected components
    int Connected_Pos_Components(std::vector<int>& nodes); // find positive connected components

    void dfs_PosConn(int beginnode, std::vector<int>& subgraph, std::vector<int>& tag);

    // set edge weights
    void setew(std::vector<int>& _ew) {
        ew.clear();
        ew = _ew;
    }
    void setMGraph(TGraph* tgraph, int l, int r) {
        ew.clear();
        ew.resize(ENum, 0);
        tgraph->getSumEW_in_lr(ew, l, r);
    } 
    void setDGraph(TGraph* tgraph, int l, int k1, int k2) {
        ew.clear();
        ew.resize(ENum, INT_MIN);
        tgraph->getMaxEW_in_lr(ew, l, k1, k2);
    }

    EGraph operator + (const EGraph& egraph) const {
        EGraph sum_egraph;
        sum_egraph.VNum = egraph.VNum;
        sum_egraph.ENum = egraph.ENum;
        sum_egraph.hashgraph = hashgraph;
        sum_egraph.ew = ew;
        for (int i = 0; i < ENum; i++) {
            sum_egraph.ew[i] += egraph.ew[i];
        }
        return sum_egraph;
    }
};

int EGraph::computeUBsop() {
    // compute UBSOP: sum all positive edges
    int UBsop = 0;
    for (int i = 0; i < ENum; i++) {
        if (ew[i] > 0) UBsop += ew[i];
    }
    return UBsop;
}

int EGraph::computeUBstr() {
    // compute UBSTR: 
    // detail of this process is demonstrated in 2011 paper
    int UBstr = 0;
    int maxPN = 0;
    int minN = 0;
    std::vector<int> nodes;
    int ans = Connected_Pos_Components(nodes);
    std::vector<std::pair<int, int>> pairC(ans, { 0, 0 });  // {P, N}

    for (auto iter1 = (*hashgraph).begin(); iter1 != (*hashgraph).end(); iter1++) {
        int u = nodes[iter1->first];
        for (auto iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++) {
            int v = nodes[iter2->first];
            int edge_weight = ew[iter2->second];
            if (edge_weight > 0 && u == v) {
                pairC[u].first += edge_weight;
            }
            else if (edge_weight < 0 && u != v) {
                if (pairC[u].second > -edge_weight) {
                    pairC[u].second = -edge_weight; // lowest score among neg edges i.e. highest score among -neg edges
                }
            }
        }
    }

    for (auto it = pairC.begin(); it != pairC.end(); it++) {
        if ((double)it->first / 2 - it->second > 0) {
            maxPN += (double)it->first / 2 - it->second;
        }
        if (it->second > minN) minN = it->second;
    }
    UBstr = maxPN + minN;
    return UBstr;
}

int EGraph::uf_find(std::vector<int>& parent, int x) {
    //int k, j, r;
    //r = x;
    //while (r != parent[r])
    //    r = parent[r];
    //k = x;
    //while (k != r) // non-recursive
    //{
    //    j = parent[k];
    //    parent[k] = r;
    //    k = j;
    //}
    //return r;
    if (parent[x] != x) // recursive
        parent[x] = uf_find(parent, parent[x]);
    return parent[x];
}

void EGraph::uf_union(std::vector<int>& parent, std::vector<int>& rank, int x, int y) {
    x = uf_find(parent, x);
    y = uf_find(parent, y);
    if (x == y) return;
    if (rank[x] < rank[y]) {
        parent[x] = y;
    }
    else {
        parent[y] = x;
        if (rank[x] == rank[y]) rank[x]++;
    }
}

int EGraph::merge(std::vector<int>& parent, int x) {
    if (parent[x] == x)
        return x;
    return merge(parent, parent[x]);
}

int EGraph::Connected_Components(std::vector<int>& nodes) {
    // standard UnionFind algorithm
    nodes.clear();
    nodes.resize(VNum, 0);
    std::vector<int> rank(VNum, 0);
    for (auto iter = nodes.begin(); iter != nodes.end(); iter++) {
        *iter = iter - nodes.begin();
    }
    for (auto x = (*hashgraph).begin(); x != (*hashgraph).end(); x++) {
        for (auto y : x->second) {
            uf_union(nodes, rank, x->first, y.first);
            // nodes[merge(nodes, x->first)] = merge(nodes, y.first);
        }
    }
    int ans = 0;
    std::map<int, std::list<int>> m;
    for (int i = 0; i < VNum; i++) {
        ans += (nodes[i] == i);
        nodes[i] = merge(nodes, nodes[i]);
        m[nodes[i]].push_back(i);
    }
    /*for (int i = 0; i < VNum; i++) {
        nodes[i] = merge(nodes, nodes[i]);
    }*/
    /*for (int i = 0; i < VNum; i++) {
        m[nodes[i]].push_back(i);
    }*/
    int i = 0;
    for (auto it = m.begin(); it != m.end(); it++) {
        std::list<int> l = it->second;
        for (auto x : l) {
            nodes[x] = i;
        }
        i++;
    }

    return ans;
}

int EGraph::Connected_Pos_Components(std::vector<int>& nodes) {
    // same as Connected_Components
    nodes.clear();
    nodes.resize(VNum, 0);
    std::vector<int> rank(VNum, 0);
    for (auto iter = nodes.begin(); iter != nodes.end(); iter++) {
        *iter = iter - nodes.begin();
    }
    for (auto x = (*hashgraph).begin(); x != (*hashgraph).end(); x++) {
        for (auto y : x->second) {
            if(ew[y.second] >= 0) // only change to find pos components
                uf_union(nodes, rank, x->first, y.first);
            // nodes[merge(nodes, x->first)] = merge(nodes, y.first);
        }
    }
    int ans = 0;
    std::map<int, std::list<int>> m;
    for (int i = 0; i < VNum; i++) {
        ans += (nodes[i] == i);
        nodes[i] = merge(nodes, nodes[i]);
        m[nodes[i]].push_back(i);
    }
    /*for (int i = 0; i < VNum; i++) {
        nodes[i] = merge(nodes, nodes[i]);
    }*/
    /*for (int i = 0; i < VNum; i++) {
        m[nodes[i]].push_back(i);
    }*/
    int i = 0;
    for (auto it = m.begin(); it != m.end(); it++) {
        std::list<int> l = it->second;
        for (auto x : l) {
            nodes[x] = i;
        }
        i++;
    }
    return ans;
}

void EGraph::dfs_PosConn(int beginnode, std::vector<int>& subgraph, std::vector<int>& tag) {
    for (auto iter = (*hashgraph)[beginnode].begin(); iter != (*hashgraph)[beginnode].end(); ++iter) {
        std::pair<int, int> edge = *iter;
        if (tag[edge.first] == 0 && ew[edge.second] >= 0) {
            subgraph.push_back(edge.first);
            tag[edge.first] = 1;
            dfs_PosConn(edge.first, subgraph, tag);
        }
    }
}

/*_________________________________________*/

/*______________MinHeap_____________________*/
// Minary Heap https://www.geeksforgeeks.org/prims-mst-for-adjacency-list-representation-greedy-algo-6/
// Structure to represent a min heap node
struct MinHeapNode {
    int v;
    int key;
};

// Structure to represent a min heap
struct MinHeap {
    int size; // Number of heap nodes present currently
    int capacity; // Capacity of min heap
    int* pos; // This is needed for decreaseKey()
    struct MinHeapNode** array;
};

// A utility function to create a new Min Heap Node
struct MinHeapNode* newMinHeapNode(int v, int key)
{
    struct MinHeapNode* minHeapNode = (struct MinHeapNode*)malloc(sizeof(struct MinHeapNode));
    minHeapNode->v = v;
    minHeapNode->key = key;
    return minHeapNode;
}

// A utilit function to create a Min Heap
struct MinHeap* createMinHeap(int capacity)
{
    struct MinHeap* minHeap = (struct MinHeap*)malloc(sizeof(struct MinHeap));
    minHeap->pos = (int*)malloc(capacity * sizeof(int));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array = (struct MinHeapNode**)malloc(capacity * sizeof(struct MinHeapNode*));
    return minHeap;
}

// A utilit function to delete a Min Heap
void deleteMinHeap(MinHeap* minHeap)
{
    int capacity = minHeap->size;
    for (int i = 0; i < capacity; i++) {
        free(minHeap->array[i]);
    }
    free(minHeap->array);
    free(minHeap->pos);
    free(minHeap);
}

// A utility function to swap two nodes of min heap. Needed for min heapify
void swapMinHeapNode(struct MinHeapNode** a, struct MinHeapNode** b)
{
    struct MinHeapNode* t = *a;
    *a = *b;
    *b = t;
}

// A standard function to heapify at given idx
// This function also updates position of nodes when they are swapped.
// Position is needed for decreaseKey()
void minHeapify(struct MinHeap* minHeap, int idx)
{
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;

    if (left < minHeap->size && minHeap->array[left]->key > minHeap->array[smallest]->key) // changed
        smallest = left;

    if (right < minHeap->size && minHeap->array[right]->key > minHeap->array[smallest]->key) // changed
        smallest = right;

    if (smallest != idx) {
        // The nodes to be swapped in min heap
        MinHeapNode* smallestNode = minHeap->array[smallest];
        MinHeapNode* idxNode = minHeap->array[idx];

        // Swap positions
        minHeap->pos[smallestNode->v] = idx;
        minHeap->pos[idxNode->v] = smallest;

        // Swap nodes
        swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);

        minHeapify(minHeap, smallest);
    }
}

// A utility function to check if the given minHeap is ampty or not
int isEmpty(struct MinHeap* minHeap)
{
    return minHeap->size == 0;
}

// Standard function to extract minimum node from heap
struct MinHeapNode* extractMin(struct MinHeap* minHeap)
{
    if (isEmpty(minHeap))
        return NULL;

    // Store the root node
    struct MinHeapNode* root = minHeap->array[0];

    // Replace root node with last node
    struct MinHeapNode* lastNode = minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;

    // Update position of last node
    minHeap->pos[root->v] = minHeap->size - 1;
    minHeap->pos[lastNode->v] = 0;

    // Reduce heap size and heapify root
    --minHeap->size;
    minHeapify(minHeap, 0);

    return root;
}

// Function to decrease key value of a given vertex v. This function
// uses pos[] of min heap to get the current index of node in min heap
void decreaseKey(struct MinHeap* minHeap, int v, int key)
{
    // Get the index of v in  heap array
    int i = minHeap->pos[v];

    // Get the node and update its key value
    minHeap->array[i]->key = key;

    // Travel up while the complete tree is not hepified.
    // This is a O(Logn) loop
    while (i && minHeap->array[i]->key > minHeap->array[(i - 1) / 2]->key) { // changed
        // Swap this node with its parent
        minHeap->pos[minHeap->array[i]->v] = (i - 1) / 2;
        minHeap->pos[minHeap->array[(i - 1) / 2]->v] = i;
        swapMinHeapNode(&minHeap->array[i], &minHeap->array[(i - 1) / 2]);

        // move to parent index
        i = (i - 1) / 2;
    }
}

// A utility function to check if a given vertex
// 'v' is in min heap or not
bool isInMinHeap(struct MinHeap* minHeap, int v)
{
    if (minHeap->pos[v] < minHeap->size)
        return true;
    return false;
}

/*_________________________________________*/

/*___________________STree_________________*/

// Spanning Tree for the TopDown algorithm
class STree {
public:
    EGraph graph;
    // std::vector<std::pair<int,int>> Tree; // parents of index
    int VNum;

    struct NWTnode { 
        int nw; // node weitht
        std::vector<std::pair<int, int>> children;  //pair of child, edgeindex
        std::vector<int> nodes;  //pair of component node
        NWTnode() : nw(0) {};
    };
    std::vector<NWTnode> NWTree; // NW-PCST-Tree
    int maxroot; // find maxroot
    int maxnw; 
    int resnw;

    std::vector<bool> resNodes; // bit_map for result_nodes
    std::vector<bool> resEdges; // bit_map for result_edges
    graph_hash_of_mixed_weighted gwgraph;

public:
    STree() {};
    STree(EGraph _graph) { graph = _graph; maxroot = 0; maxnw = INT_MIN; resNodes.resize(graph.VNum, 0); resEdges.resize(graph.ENum, 0); resnw = 0; };
    STree(const STree& st) {
        resnw = st.resnw;
        resNodes = st.resNodes;
        resEdges = st.resEdges;
    };
    /*~STree() {
        resNodes.~vector(); resEdges.~vector(); NWTree.~vector();
    }*/

    void PrimMST(std::vector<int> nodes); // use Prim and binary heap to find maximum spanning tree 
    void printMST();
    void printMST(int u);
    void mergeNWTree(int u); // merge nodes, convert MST to NW-PCST-Tree
    void SPrA(int u); // Strong Pruning process in PCST-algorithm
    void getResNodes(int u); // get result nodes
    void addPosEdges(); // spread NW-PCST-Tree
    void TopDown(); // TopDown Algorithm
    std::unordered_map <int, std::vector<std::pair<int, int>>> convert_to_map(); // convert resEdges to unordered_map

};

// The main function that constructs Minimum Spanning Tree (MST)
// using Prim's algorithm
void STree::PrimMST(std::vector<int> nodes)
{
    VNum = graph.VNum; // Get the number of vertices in graph
    std::vector<std::pair<int, int>> Tree(VNum, { -1, 0 }); // Array to store constructed MST
    std::vector<int> key(VNum, INT_MIN); // Key values used to pick maxmum weight edge in cut

    // init NWTree
    NWTree.clear();
    NWTree.resize(VNum);
    NWTree[0].nw = 0;
    // NWTree[0].nodes.push_back({0,0});

    // minHeap represents set E
    struct MinHeap* minHeap = createMinHeap(VNum);

    // Initialize min heap with all vertices. Key value of
    // all vertices (except 0th vertex) is initially infinite
    for (int v = 1; v < VNum; v++) {
        minHeap->array[v] = newMinHeapNode(v, key[v]);
        minHeap->pos[v] = v;
        NWTree[v].nw = 0;
        // NWTree[v].nodes.push_back({ v ,0});
    }

    // Make key value of 0th vertex as 0 so that it
    // is extracted first (begin nodes)
    for (auto it = nodes.begin(); it != nodes.end(); it++) {
        key[*it] = 0;
        minHeap->array[*it] = newMinHeapNode(*it, key[*it]);
        minHeap->pos[*it] = 0;
    }

    // Initially size of min heap is equal to V
    minHeap->size = VNum;


    // In the following loop, min heap contains all nodes
    // not yet added to MST.
    while (!isEmpty(minHeap)) {
        // Extract the vertex with minimum key value
        struct MinHeapNode* minHeapNode = extractMin(minHeap);
        int u = minHeapNode->v; // Store the extracted vertex number

        // Traverse through all adjacent vertices of u (the extracted
        // vertex) and update their key values
        int pCrawl = (*graph.hashgraph)[u].size();
        for (int i = 0; i < pCrawl; i++) {
            int v = (*graph.hashgraph)[u][i].first;
            int edgeindex = (*graph.hashgraph)[u][i].second;
            // If v is not yet included in MST and weight of u-v is
            // less than key value of v, then update key value and
            // parent of v
            int vweight = graph.ew[edgeindex];
            if (isInMinHeap(minHeap, v) && vweight > key[v]) {
                key[v] = vweight;
                Tree[v] = { u , edgeindex };
                decreaseKey(minHeap, v, key[v]);
            }
        }
    }
    for (int v = 1; v < VNum; v++) {
        if (Tree[v].first != -1)
            NWTree[Tree[v].first].children.push_back({ v ,Tree[v].second });
    }
    deleteMinHeap(minHeap);
}

// A utility function used to print the constructed MST
//void STree::printMST()
//{
//    //for (int i = 1; i < VNum; i++)
//    //    printf("%d <- %d - %d\n", Tree[i].first,Tree[i].second, i);
//    //printf("\n");
//    for (int i = 0; i < VNum; i++) {
//        printf("node: %d\t nw: %d\t nodes: ", i, NWTree[i].nw);
//        for (int j = 0; j < NWTree[i].nodes.size(); j++) {
//            printf("%d\t", NWTree[i].nodes[j]);
//        }
//        printf("children: ");
//        for (int j = 0; j < NWTree[i].children.size(); j++) {
//            printf("%d %d\t", NWTree[i].children[j].first, NWTree[i].children[j].second);
//        }
//        printf("\n");
//    }
//}
//
//void STree::printMST(int u)
//{
//    printf("node: %d\t nw: %d\t nodes: ", u, NWTree[u].nw);
//    for (int j = 0; j < NWTree[u].nodes.size(); j++) {
//        printf("%d\t", NWTree[u].nodes[j]);
//    }
//    printf("children: ");
//    for (int j = 0; j < NWTree[u].children.size(); j++) {
//        printf("%d %d\t", NWTree[u].children[j].first, NWTree[u].children[j].second);
//    }
//    printf("\n");
//    for (auto iter = NWTree[u].children.begin(); iter != NWTree[u].children.end(); iter++) {
//        int v = (*iter).first;
//        int edgeindex = (*iter).second;
//        int edgeweight = graph.ew[edgeindex];
//        printMST(v);
//    }
//}

void STree::mergeNWTree(int u) {
    // begin from u
    for (auto iter = NWTree[u].children.begin(); iter != NWTree[u].children.end();) {
        int v = iter->first;
        int edgeindex = iter->second;
        int edgeweight = graph.ew[edgeindex];
         
        mergeNWTree(v);
        if (edgeweight >= 0) {
            NWTree[u].nodes.push_back(v);
            NWTree[u].nodes.insert(NWTree[u].nodes.end(), NWTree[v].nodes.begin(), NWTree[v].nodes.end());
            NWTree[v].nodes.clear();
            NWTree[u].nw += NWTree[v].nw;
            NWTree[u].nw += edgeweight;
            iter = NWTree[u].children.erase(iter);
            iter = NWTree[u].children.insert(iter, NWTree[v].children.begin(), NWTree[v].children.end()) + NWTree[v].children.size();
            NWTree[v].children.clear();
        }
        else {
            iter++;
        }
    }
}

// Strong Pruncing Approuch of NW-PCST Tree
void STree::SPrA(int u) {
    for (auto iter = NWTree[u].children.begin(); iter != NWTree[u].children.end(); ) {
        int v = (*iter).first;
        int edgeindex = (*iter).second;
        int edgeweight = graph.ew[edgeindex];
        SPrA(v);
        if (-edgeweight >= NWTree[v].nw) {
            iter = NWTree[u].children.erase(iter);
            if (NWTree[v].nw > maxnw) {
                maxnw = NWTree[v].nw;
                maxroot = v;
            }
        }
        else {
            NWTree[u].nw = NWTree[u].nw + NWTree[v].nw + edgeweight;
            iter++;
        }
    }
    if (NWTree[u].nw > maxnw) {
        maxnw = NWTree[u].nw;
        maxroot = u;
    }
}

void STree::getResNodes(int u) {
    resNodes[u] = 1;
    for (auto iter = NWTree[u].nodes.begin(); iter != NWTree[u].nodes.end(); iter++) {
        int v = (*iter);
        getResNodes(v);
    }
    for (auto iter = NWTree[u].children.begin(); iter != NWTree[u].children.end(); iter++) {
        int v = (*iter).first;
        int edgeindex = (*iter).second;
        resEdges[edgeindex] = 1;
        getResNodes(v);
    }
}

void STree::addPosEdges() {
    resnw = 0;
    for (std::unordered_map <int, std::vector<std::pair<int, int>>>::iterator it = (*graph.hashgraph).begin(); it != (*graph.hashgraph).end(); it++) {
        int temp_u = it->first;
        graph_hash_of_mixed_weighted_add_vertex(gwgraph, temp_u, 0);
        for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
            int temp_v = it2->first;
            int temp_edgeindex = it2->second;
            int temp_edgeweight = graph.ew[temp_edgeindex];
            if (temp_edgeweight > 0 && (resNodes[temp_u] == 1 || resNodes[temp_v] == 1)) {
                resEdges[temp_edgeindex] = 1;
                graph_hash_of_mixed_weighted_add_vertex(gwgraph, temp_v, 0);
                graph_hash_of_mixed_weighted_add_edge(gwgraph, temp_u, temp_v, graph.ew[temp_edgeindex]);
            }
        }
    }
    for (auto it = resEdges.begin(); it != resEdges.end(); it++) {
        if (*it) {
            resnw += graph.ew[(it - resEdges.begin())];
        }
    }
}

void STree::TopDown() {
    std::vector<bool> tmp_bitmap_nodes(graph.VNum, 0);
    int tmp_max = INT_MIN;
    // find Connected_Components
    std::vector<int> nodes;
    std::vector<int> ans_nodes(1, 0);
    int ans = graph.Connected_Components(nodes);
    for (auto it = nodes.begin(); it != nodes.end(); it++) {
        if (tmp_bitmap_nodes[*it] == 0) {
            tmp_bitmap_nodes[*it] = 1;
            ans_nodes.push_back(*it);
        } 
    }
    nodes.clear();
    PrimMST(ans_nodes);
    // compute all components
    for (auto it = ans_nodes.begin(); it != ans_nodes.end(); it++) {
        mergeNWTree(*it);
        SPrA(*it);
    }
    // find root with maximum weight
    // mergeNWTree(0);
    // SPrA(0);
    getResNodes(maxroot);
    addPosEdges();
    // maxroot = 0;
    NWTree.clear();
}

std::unordered_map <int, std::vector<std::pair<int, int>>> STree::convert_to_map() {
    // convert resEdges to edges
    std::unordered_map <int, std::vector<std::pair<int, int>>> resmap;
    for (auto it = (*graph.hashgraph).begin(); it != (*graph.hashgraph).end(); it++) {
        for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
            if (resEdges[it2->second]) {
                int u = it->first;
                int v = it2->first;
                resmap[u].push_back({ v, it2->second });
            }
        }
    }
    return resmap;
}

/*_______________________________________*/

/*_________________PCST__________________*/

class PCST {
public:
    graph_hash_of_mixed_weighted wgraph;
    int resnw;

public:
    PCST() {};
    PCST(EGraph& egraph); // convert EGraph to graph_hash_of_mixed_weighted

    void convert_to_wgraph(EGraph& egraph);  // convert egraph to PCST graph
    void computePCST();  // solve PCST problem, including growing and pruning
    void print_PCST_Result();

    std::unordered_map <int, std::vector<std::pair<int, int>>> convert_to_map(); // convert to unordered_map
};

void PCST::computePCST() {
    graph_hash_of_mixed_weighted tempgraph = graph_hash_of_mixed_weighted_Fast_GW_Growing_tree(wgraph, "forest"); // tree
    graph_hash_of_mixed_weighted resgraph = graph_hash_of_mixed_weighted_General_Pruning_forest(tempgraph); //Pruning Tree

    resnw = 0;
    for (auto it = resgraph.hash_of_vectors.begin(); it != resgraph.hash_of_vectors.end(); it++) {
        resnw += it->second.vertex_weight * 2;
        for (auto it2 = it->second.adj_vertices.begin(); it2 != it->second.adj_vertices.end(); it2++) {
            resnw -= it2->second;
        }
    }
    resnw = resnw / 2;
    wgraph = resgraph;
}

PCST::PCST(EGraph& egraph) {
    convert_to_wgraph(egraph);
}

void PCST::convert_to_wgraph(EGraph& egraph) {
    std::vector<int> nodes;
    int ans = egraph.Connected_Pos_Components(nodes);
    std::vector<int> new_weight(egraph.VNum, 0);

    // add edge
    for (auto x = (*egraph.hashgraph).begin(); x != (*egraph.hashgraph).end(); x++) {
        int u = nodes[x->first];
        for (auto y : x->second) {
            int v = nodes[y.first];
            int ew = egraph.ew[y.second];
            // in the same component
            if (ew > 0 && u == v) {
                new_weight[u] += ew;
            }
            // in different components
            else if (ew < 0 && u != v) {
                if (graph_hash_of_mixed_weighted_contain_edge(wgraph, u, v)) {
                    if (graph_hash_of_mixed_weighted_edge_weight(wgraph, u, v) > -ew) {
                        graph_hash_of_mixed_weighted_add_edge(wgraph, u, v, -ew);
                    }
                }
                else {
                    graph_hash_of_mixed_weighted_add_edge(wgraph, u, v, -ew);
                }
            }
        }
    }
    // add vertex weight
    for (int i = 0; i < egraph.VNum; i++) {
        if(new_weight[i] != 0)
            graph_hash_of_mixed_weighted_add_vertex(wgraph, i, new_weight[i] / 2);
    }

}

void PCST::print_PCST_Result() {
    std::cout << resnw << std::endl;
    graph_hash_of_mixed_weighted_print(wgraph);
}

std::unordered_map <int, std::vector<std::pair<int, int>>> PCST::convert_to_map() {
    std::unordered_map <int, std::vector<std::pair<int, int>>> resmap;
    int temp_edge_index = 0;
    for (auto it = wgraph.hash_of_vectors.begin(); it != wgraph.hash_of_vectors.end(); it++) {
        if (it->second.adj_vertices.size() == 0) {
            resmap[it->first].resize(0);
        }
        else {
            for (auto it2 = it->second.adj_vertices.begin(); it2 != it->second.adj_vertices.end(); it2++) {
                resmap[it->first].push_back({ it2->first, temp_edge_index++ });
            }
        }
    }
    return resmap;
}
/*_______________________________________*/

/*____________GroupsofGraph______________*/

class GroupS {
public:
    TGraph* tgraph;
    int T;
    double alpha;
    int ENum;

    struct Tnode {
        std::vector<int> mw; // merge weight
        std::vector<int> dw; // dominate weight
        std::pair<int, int> lr; // time interval [l, r]
        Tnode* left, * right; // left and right child
        Tnode() : left(NULL), right(NULL) {};
    };
    Tnode* Index_Root;  // Root of G index (detail in 2011 paper)

    struct group {
        std::tuple<int, int, int> index; // S^i={S(i,k1,k2)}
        // std::vector<int> dw; // dominate weight
        int UB;  // Upper Bound
    };
    std::vector<group> groups; // bar{G}(i, (k1, k2))

    std::vector<std::pair<int, int>> prunedres; // prune results of intervals
    int maxUB; // maximum Upper Bound
    int LB;  // LowerBound of HDS


    int resvalue; // result value
    std::pair<int, int> reslr; // result time interval [l, r]
    graph_hash_of_mixed_weighted resmap; // result vector of edges

public:
    GroupS(TGraph* _tgraph, double _alpha) {
        Index_Root = NULL; maxUB = INT_MIN; resvalue = 0;
        alpha = _alpha;
        tgraph = _tgraph;
        T = tgraph->T;
        ENum = tgraph->ENum;

        // Generate all groups
        for (int i = 0; i < T; i++) {
            group tgb;
            tgb.index = { i, i ,i };
            groups.push_back(tgb);
            int j = 1;
            int k1 = i + (int)(1.0 / pow(alpha, j - 1));
            int k2 = i + (int)(1.0 / pow(alpha, j)) - 1;
            while (k2 < T - 1) {
                group tgt;
                tgt.index = { i, k1, k2 };
                groups.push_back(tgt);
                j++;
                k1 = i + (int)(1 / pow(alpha, j - 1));
                k2 = i + (int)(1 / pow(alpha, j)) - 1;
            }
            if (k1 < T && k2 >= T - 1) {
                group tgl;
                tgl.index = { i, k1, T - 1 };
                groups.push_back(tgl);
            }
        }
        LB = -1;
    }

    Tnode* buildIndex(int l, int r); // build Graph Index Tree for search
    void getmw(Tnode* node, int i, int k1, std::vector<int>& mw);
    void getdw(Tnode* node, int k1, int k2, std::vector<int>& dw);  // get dominate weight in Graph Index Tree
    void computeDGraph(std::string ub);  // compute UB of DGraph
    void estimateLB(std::string alg); // estimate lower bound
    void pruning(std::string ub); // prune groups and intervals

    void MEDEN(std::string alg, std::string ub); // MEDEN 
    void print_TopDown_Result(STree& resTree);
    void print_PCST_Result(PCST& pcst);
    void print_Result();
};

GroupS::Tnode* GroupS::buildIndex(int l, int r) {
    if (l == r) {
        Tnode* leaf = new Tnode();
        leaf->dw.resize(ENum, 0);
        tgraph->getEW_inT(leaf->dw, l);
        leaf->mw = leaf->dw;
        leaf->lr = { l,r };
        leaf->left = NULL; leaf->right = NULL;
        return leaf;
    }
    int mid = (l + r) / 2;
    Tnode* leaf = new Tnode();
    leaf->left = buildIndex(l, mid);
    leaf->right = buildIndex(mid + 1, r);
    std::vector<int> mw(ENum, 0);
    std::vector<int> dw;
    for (int i = 0; i < ENum; i++) {
        mw[i] = leaf->left->mw[i] + leaf->right->mw[i];
        int tempright = leaf->right->dw[i] + leaf->left->mw[i];
        int temp_max = leaf->left->dw[i] > tempright ? leaf->left->dw[i] : tempright;
        dw.push_back(temp_max);
    }
    leaf->mw = mw;
    leaf->dw = dw;
    leaf->lr = { l,r };
    return leaf;
}

void GroupS::getmw(Tnode* node, int i, int k1, std::vector<int>& mw) {
    int l = node->lr.first;
    int r = node->lr.second;
    int mid = (l + r) / 2;
    if (i <= l && r <= k1) {
        for (int a = 0; a < ENum; a++) {
            mw[a] += node->mw[a];
        }
        return;
    }
    if (node->left != NULL) {
        if (i <= mid ) {
            getmw(node->left, i, k1, mw);
        }
    }
    if (node->right != NULL) {
        if (k1 > mid) {
            getmw(node->right, i, k1, mw);
        }
    }
}

void GroupS::getdw(Tnode* node, int k1, int k2, std::vector<int>& dw) {
    std::stack<Tnode*> s;
    Tnode *p = node;
    int l, r, mid;
    // std::vector<int> mw(ENum, 0);
    // do Traversal, from right nodes (deep)
    while (p || !s.empty()) {
        while (p) {
            s.push(p);
            l = p->lr.first;
            r = p->lr.second;
            mid = (l + r) / 2;
            if (k1 <= l && r <= k2) break; // find node
            if (k2 <= mid) break; // not in right nodes
            p = p->right;
        }
        p = s.top();
        s.pop();
        l = p->lr.first;
        r = p->lr.second;
        mid = (l + r) / 2;
        if (k1 <= l && r <= k2) {
            for (int j = 0; j < ENum; j++) {
                int temp_dw = p->mw[j] + dw[j];
                dw[j] = p->dw[j] > temp_dw ? p->dw[j] : temp_dw;
            }
            p = NULL;
            continue;
        }
        if (k1 <= mid) {
            p = p->left;
        }
        else { p = NULL; }
    }
    //int l = node->lr.first;
    //int r = node->lr.second;
    //int mid = (l + r) / 2;
    //if (k1 <= l && r <= k2) {
    //    mw = node->mw;
    //    dw = node->dw;
    //    /*for (int a = 0; a < ENum; a++) {
    //        mw[a] += node->mw[a];
    //    }*/
    //    return;
    //}
    //int tag = 0;
    //std::vector<int> temp_mwr(ENum, 0);
    //std::vector<int> temp_dwr(ENum, INT_MIN);
    //// std::vector<int> temp_dwr(ENum, INT_MIN);
    //if (node->right != NULL) {
    //    if (k2 > mid) {
    //        getdw(node->right, k1, k2, temp_mwr, temp_dwr);
    //        tag = 1;
    //    }
    //}
    //std::vector<int> temp_mwl(ENum, 0);
    //std::vector<int> temp_dwl(ENum, INT_MIN);
    //if (node->left != NULL) {
    //    if (k1 <= mid) {
    //        getdw(node->left, k1, k2, temp_mwl, temp_dwl);
    //        if (tag == 1) tag = 2;
    //        else tag = 3;
    //    }
    //}
    //if (tag == 1) {
    //    mw = temp_mwr;
    //    dw = temp_dwr;
    //}
    //else if (tag == 2) {
    //    mw = temp_mwl;
    //    dw = temp_dwl;
    //}
    //else if (tag == 3) {
    //    for (int j = 0; j < ENum; j++) {
    //        mw[j] = temp_mwl[j] + temp_mwr[j];
    //        int dm_sum = temp_mwl[j] + temp_dwr[j];
    //        dw[j] = temp_mwl[j] > dm_sum ? temp_mwl[j] : dm_sum;
    //    }
    //}
}

// compute UB and dw in groups
void GroupS::computeDGraph(std::string ub) {
    maxUB = -1;
    for (auto iter = groups.begin(); iter != groups.end(); iter++) {
        std::tuple<int, int, int> tempS = iter->index;
        int i = std::get<0>(tempS);
        int k1 = std::get<1>(tempS);
        int k2 = std::get<2>(tempS);

        std::vector<int> mw(ENum, 0);
        getmw(Index_Root, i, k1, mw);

        std::vector<int> dw(ENum, 0);
        if (k1 + 1 > k2) {  }
        else {
            getdw(Index_Root, k1 + 1, k2, dw);
        }
        for (int j = 0; j < ENum; j++) { 
            int dm_sum = dw[j] + mw[j];
            dw[j] = mw[j] > dm_sum ? mw[j] : dm_sum;
        }
        // test
        /*std::vector<int> test_dw(ENum, INT_MIN);
        tgraph->getMaxEW_in_lr(test_dw, i, k1, k2);
        for (int j = 0; j < tgraph->ENum; j++) {
            if (dw[j] != test_dw[j]) cout << "get_dgraph error" << endl;
        }
        cout << "test: " << i << " "<< k1 << " " << k2 << endl;*/
        // iter->dw = dw;
        EGraph dgraph = EGraph(tgraph);
        dgraph.setew(dw);
        if (ub == "UBsop")
            iter->UB = dgraph.computeUBsop();
        else if (ub == "UBstr")
           iter->UB = dgraph.computeUBstr();
        if (iter->UB > maxUB) {
            maxUB = iter->UB;
        }
        // dgraph.~EGraph();
        // (*iter).UB = dgraph.computeUBsop();
    }
}

// estimate Lower Bound
void GroupS::estimateLB(std::string alg) {
    // compute max indexs
    LB = 0;
    int temp_max_index = INT_MIN;
    std::pair<int, int> max_index;
    for (auto iter = groups.begin(); iter != groups.end(); iter++) {
        if (iter->UB == maxUB) {
            std::tuple<int, int, int> maxS = (*iter).index;
            int i = std::get<0>(maxS);
            int k1 = std::get<1>(maxS);
            int k2 = std::get<2>(maxS);

            // estimateLB
            for (int j = k1; j <= k2; j++) {
                EGraph temp_graph = EGraph(tgraph);
                temp_graph.setMGraph(tgraph, i, j);
                // temp_graph.ew.resize(ENum, 0);
                // getmw(Index_Root, i, j, temp_graph.ew);
                int UBtemp = temp_graph.computeUBsop();
                /*if (alg == "TopDown") {
                    STree mst = STree(temp_graph);
                    mst.TopDown();
                    LBtemp = mst.resnw;
                }
                else if (alg == "PCST") {
                    PCST npcst(temp_graph);
                    npcst.computePCST();
                    LBtemp = npcst.resnw;
                }*/
                if (UBtemp > temp_max_index) {
                    temp_max_index = UBtemp;
                    max_index = {i, j};
                }
                // temp_graph.~EGraph();
            }
        }
    }
    EGraph temp_graph = EGraph(tgraph);
    temp_graph.setMGraph(tgraph, max_index.first, max_index.second);
    // temp_graph.ew.resize(ENum, 0);
    // getmw(Index_Root, max_index.first, max_index.second, temp_graph.ew);
    if (alg == "TopDown") {
        STree mst = STree(temp_graph);
        mst.TopDown();
        LB = mst.resnw;
    }
    else if (alg == "PCST") {
        PCST npcst(temp_graph);
        npcst.computePCST();
        LB = npcst.resnw;
    }
    // temp_graph.~EGraph();
}

// Pruning intervals
void GroupS::pruning(std::string ub) {
    for (auto iter = groups.begin(); iter != groups.end(); iter++) {
        if (iter->UB >= LB) {
            std::tuple<int, int, int> index = (*iter).index;
            int i = std::get<0>(index);
            int k1 = std::get<1>(index);
            int k2 = std::get<2>(index);
            // cout << i <<" "<< k1<<" "<< k2 << " " << endl;
            // for all intervals, compute UB and compared with LB
            for (int j = k1; j <= k2; j++) {
                EGraph temp_graph = EGraph(tgraph);
                // temp_graph.ew.resize(ENum, 0);
                // getmw(Index_Root, i, j, temp_graph.ew);
                temp_graph.setMGraph(tgraph, i, j);
                int UB_temp = 0;
                if (ub == "UBsop")
                    UB_temp = temp_graph.computeUBsop();
                else if (ub == "UBstr")
                    UB_temp = temp_graph.computeUBstr();
                if (UB_temp >= LB) prunedres.push_back({ i, j });
                // temp_graph.~EGraph();
            }
        }
    }
}

// MEDEN algorithm, return STree as result
void GroupS::MEDEN(std::string alg, std::string ub) {

    Index_Root = buildIndex(0, T - 1);
    computeDGraph(ub);
    estimateLB(alg);
    pruning(ub);

    if (alg == "TopDown") {
        STree resTree;
        int maxRes = -1;
        for (auto iter = prunedres.begin(); iter != prunedres.end(); iter++) {
            int l = iter->first;
            int r = iter->second;
            EGraph temp_graph(tgraph);
            // temp_graph.ew.resize(ENum, 0);
            // getmw(Index_Root, l, r, temp_graph.ew);
            temp_graph.setMGraph(tgraph, l, r);
            STree mst = STree(temp_graph);
            mst.TopDown();
            if (mst.resnw > maxRes) {
                maxRes = mst.resnw;
                resTree = mst;
                reslr = *iter;
            }
        }
        resvalue = resTree.resnw;
        resmap = resTree.gwgraph;
        // print_TopDown_Result(resTree);
    }
    else if (alg == "PCST") {
        PCST respcst;
        int maxRes = -1;
        for (auto iter = prunedres.begin(); iter != prunedres.end(); iter++) {
            int l = iter->first;
            int r = iter->second;
            EGraph temp_graph(tgraph);
            // temp_graph.ew.resize(ENum, 0);
            // getmw(Index_Root, l, r, temp_graph.ew);
            temp_graph.setMGraph(tgraph, l, r);
            PCST mst(temp_graph);
            mst.computePCST();
            if (mst.resnw > maxRes) {
                maxRes = mst.resnw;
                respcst = mst;
                reslr = *iter;
            }
        }
        resvalue = respcst.resnw;
        resmap = respcst.wgraph;

        // print_PCST_Result(respcst);
    }
}

//void GroupS::print_TopDown_Result(STree& resTree) {
//    std::cout << "\nCompute Result: " << resTree.resnw << "\n";
//    std::cout << "\nTime [l, r]: [" << reslr.first << ", " << reslr.second << "]\n";
//    std::cout << "nodes: ";
//    for (auto iter = resTree.resNodes.begin(); iter != resTree.resNodes.end(); iter++) {
//        if ((*iter) == 1)
//            std::cout << iter - resTree.resNodes.begin() << "\t";
//    }
//    std::cout << "\nedges: ";
//    for (auto iter = resTree.resEdges.begin(); iter != resTree.resEdges.end(); iter++) {
//        if ((*iter) == 1) std::cout << iter - resTree.resEdges.begin() << "\t";
//    }
//}
//
//void GroupS::print_PCST_Result(PCST& pcst) {
//    std::cout << "\nCompute Result: " << pcst.resnw << "\n";
//    std::cout << "\nTime [l, r]: [" << reslr.first << ", " << reslr.second << "]\n";
//    graph_hash_of_mixed_weighted_print(pcst.wgraph);
//}
//
//void GroupS::print_Result() {
//    std::cout << "\nCompute Result: " << resvalue << "\n";
//    std::cout << "\nTime [l, r]: [" << reslr.first << ", " << reslr.second << "]\n";
//    graph_hash_of_mixed_weighted_print(resmap);
//    /*std::cout << "Result_Nodes: \n";
//    for (auto it = resmap.begin(); it != resmap.end(); it++) {
//        printf_s("%d\t", it->first);
//    }
//    printf_s("\n");
//    std::cout << "Result_Edges: \n";
//    for (auto it = resmap.begin(); it != resmap.end(); it++) {
//        for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
//            printf_s("%d %d\n", it->first, it2->first);
//        }
//    }*/
//}

/*_____________________________________*/

/*____________RandomGraph______________*/

/*this function generates a random graph with vertex and edge weights, and this graph
may not be connected.*/
void generate_random_graph(TGraph* random_graph, int V, int E, int T
    , boost::random::mt19937& boost_random_time_seed) {

    /*time complexity: O(|V||E|)*/

    //double precision = std::pow(10, input_precision);
    /*add edges to random_graph*/
    int max_E = V * (V - 1) / 2; // must use long long int for large V
    if (E == max_E) { // complete graphs
        /*time complexity: O(|E|)*/
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < i; j++) {
                //std::vector <int> w;
                //for(int k = 0; k < T; k++){
                //    int new_weight = dist_tw(boost_random_time_seed); // generate random weights of timestamp
                //    w.push_back(new_weight);
                //}
                random_graph->insert_edge(i, j);
            }
        }
    }
    else if (E > max_E) {
        std::cout << "E: " << E << std::endl;
        std::cout << "V * (V - 1) / 2: " << max_E << std::endl;
        std::cout << "E > V * (V - 1) / 2 in graph_hash_of_mixed_weighted_generate_random_graph!" << '\n';
        exit(1);
    }
    else { // incomplete graphs

        /*time complexity: O(|V|)*/
        std::vector<int> not_full_vertices; // vertices without a full degree
        for (int i = 0; i < V; i++) {
            not_full_vertices.insert(not_full_vertices.end(), i);
        }

        /*time complexity: O(|V||E|)*/
        int edge_num = 0;
        while (edge_num < E) {
            boost::random::uniform_int_distribution<> dist_id
            { static_cast<int>(0), static_cast<int>(not_full_vertices.size() - 1) };
            int RAND = dist_id(boost_random_time_seed); // generate int random number  0, V-1
            if ((*random_graph).vertex_degree(not_full_vertices[RAND]) < V - 1) { // this is a vertex without a full degree

                /*time complexity: O(|V|)*/
                std::vector<int> unchecked(V);
                std::iota(std::begin(unchecked), std::end(unchecked), 0);
                bool added = false;
                while (added == false) {
                    boost::random::uniform_int_distribution<> dist_id2
                    { static_cast<int>(0), static_cast<int>(unchecked.size() - 1) };
                    int x = dist_id2(boost_random_time_seed);
                    int j = unchecked[x];
                    if (not_full_vertices[RAND] != j &&
                        random_graph->contain_edge(not_full_vertices[RAND], j) == 0) {
                        // This edge does not exist
                        //std::vector <int> w;
                        //for (int k = 0; k < T; k++) {
                        //	int new_weight = dist_tw(boost_random_time_seed); // generate random weights of timestamp
                        //	w.push_back(new_weight);
                        //}
                        random_graph->insert_edge(not_full_vertices[RAND], j); // add a new edge
                        edge_num++;
                        added = true;
                        break; // break after adding one edge
                    }
                    else {
                        unchecked.erase(unchecked.begin() + x);
                    }
                }
            }
            else { // this is a vertex with a full degree
                not_full_vertices.erase(not_full_vertices.begin() + RAND);
            }
        }



    }
    random_graph->VNum = V;
    random_graph->edge_sort();
    //boost::random::uniform_int_distribution<> dist_ec{ static_cast<int>(ec_min* precision), static_cast<int>(ec_max* precision) };
    int edge_num = random_graph->ENum;
    boost::random::uniform_int_distribution<> dist_tw{ static_cast<int>(0), static_cast<int>(1) };

    for (int i = 0; i < edge_num; i++) {
        std::vector<int> tempew;
        boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };
        for (int j = 0; j < T; j++) {
            int temp = dist_tw(boost_random_time_seed);
            if (temp == 0) temp = -1;
            // cout << temp << " ";
            tempew.push_back(temp);
        }
        random_graph->insert_time(tempew);
    }
    random_graph->T = T;
}

void generate_random_temporal_graph(TGraph* random_graph, int V, int E, int T) {

    /*parameters*/

    boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };
    generate_random_graph(random_graph, V, E, T, boost_random_time_seed);

}

/*_______________________________________*/

/*_________________MEDEN Algotithm____________________*/

/* MEDEN Algorithm:
   Input: graph_data; temporal data; algorithm optin "TopDown" or "PCST"; upper bound option "UBsop" or "UBstr"; alpha in [0, 1]
   Output: result value and vector of edges(result subgraph)
   */
std::pair<int, std::pair<std::pair<int, int>, graph_hash_of_mixed_weighted > > MEDEN_interface(std::vector <std::pair <int, int>>& graph_data, std::vector <std::vector <int>>& temporal_data,
    int begin_time_l, int end_time_r, 
    std::string alg_option, std::string ub_option, double alpha) {

    TGraph* tgraph = new TGraph();
    tgraph->generate_tgraph_from_two_vectors(graph_data, temporal_data, begin_time_l, end_time_r);
    // alg_option: TopDown or PCST
    // ub_option: UBsop or UBstr
    GroupS groups = GroupS(tgraph, alpha);
    groups.MEDEN(alg_option, ub_option);
    return { groups.resvalue, {{groups.reslr.first + begin_time_l, groups.reslr.second + begin_time_l}, groups.resmap } };
}


// test MEDEN in random temporal graph

void test_MEDEN(int VNum, int ENum, int T) {
    TGraph* tgraph = new TGraph();
    // tgraph->build_from_file("map_file.txt", "time_file.txt", 0 , 4);
    generate_random_temporal_graph(tgraph, VNum, ENum, T);

    time_t start, end;
    double cost;
    time(&start);
    GroupS groups = GroupS(tgraph, 0.5);
    groups.MEDEN("TopDown", "UBsop");
    time(&end);
    std::cout << groups.resvalue << endl;
    std::cout << "[" << groups.reslr.first << ", " << groups.reslr.second << "]" << endl;
    cost = difftime(end, start);
    printf("%f\n", cost);
    // groups.print_Result();
    time(&start);
    GroupS groups1 = GroupS(tgraph, 0.5);
    groups1.MEDEN("PCST", "UBsop");
    time(&end);
    std::cout << groups1.resvalue <<endl;
    std::cout << "[" << groups1.reslr.first << ", " << groups1.reslr.second << "]" << endl;
    cost = difftime(end, start);
    printf("%f\n", cost);
    if (groups.resvalue == 0 || groups1.resvalue == 0);
    printf("________________\n");
    // groups.~GroupS();
    // groups1.~GroupS();
    // groups1.print_Result();
}

/*_________________________________________*/