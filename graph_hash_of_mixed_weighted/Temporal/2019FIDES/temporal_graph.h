#pragma once

#include <queue>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <utility>
#include <cmath>

#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>

//disjoint_set_union is a class for disjoint set union data structure
class disjoint_set_union{
private:
    int n, com_num; //size of set and number of components in graph
    std::unordered_map <int, int> fa; //father node of each node

public:
    disjoint_set_union() 
    {
        n = 0;
        com_num = 0;
    } 

    int size();
    int component_num();
    int find(int x);
    int _union(int x, int y); //union two nodes, return 1 if union succeed, return 0 if x and y are already in one set
    int union_check(int x,int y);//return 1 if x and y are in one set, else return 0   
};

class temporal_graph {
private:

    std::vector <int> all_edge_sum; //all_edge_sum[i] = \sum all_edge_weight_in_graph at time i
    std::vector <std::vector <int>> edge_sum;//edge_sum[i][j] = \sum (t = 0 to j) edge_has_index_i's weight at time t
    std::map <std::pair <int, int>, int> neg_edge_map;//two map used for record edge convert from origin graph to aggregate graph
    std::unordered_map <int, int> pos_edge_map;

public:

    int max_ts, edge_num, vertex_num;//max_ts: max_time_stamp

    temporal_graph() 
    {
    }

    temporal_graph(std::vector <std::pair <int, int>> &g, int l, int r) 
    {
        max_ts = r - l + 1;
        edge_num = g.size();
        vertex_num = -1;
    }

    void prefix_sum(std::vector <std::vector <int>> &temporal_data, int l, int r);//init prefix sum, call prefix_sum before using the next two func

    int cohesive_density(int l, int r);//calc cohesive density of the whole graph in time interval [l, r]
    int positive_cohesive_density(int l, int r);//calc positive cohesive density of the whole graph in time interval [l, r]
    int get_maxts();
    int get_vnum();

    void static_graph(graph_hash_of_mixed_weighted &g, std::vector <std::pair <int, int>> &e, int l, int r);//convert temporal graph in time interval[l,r] to static aggregate graph
    void origin_graph(std::unordered_map<int, int> &g, graph_hash_of_mixed_weighted &origin_g, std::vector <std::pair <int, int>> &e);//convert aggregate graph to origin tempoal graph
};

int disjoint_set_union::size()
{
    return n;
}

int disjoint_set_union::component_num()
{
    return com_num;
}

int disjoint_set_union::find(int x)
{
    int rt = x;
    while(fa.find(rt) != fa.end() && fa[rt] != rt) {
        rt = fa[rt];
    }//find root node of x


    if (fa.find(rt) == fa.end()) {
        fa[rt] = rt;
        com_num++;
        n++;
    }

    int v = x;
    int rt_v;   
    while(v != rt) {
        rt_v = fa[v];
        fa[v] = rt;
        v = rt_v;
    }//path merge

    return rt;
}

int disjoint_set_union::_union(int x, int y)
{
    int rt_x = find(x), rt_y = find(y);
    if (rt_x != rt_y) {
        fa[rt_y] = rt_x;
        com_num --;
        return 1;
    }
    return 0;
}

int disjoint_set_union::union_check(int x, int y)
{
    int rt_x = find(x), rt_y = find(y);
    if (rt_x != rt_y)
        return 0;
    return 1;
}

//this func calc prefix sum in time interval [l, r], result sum store in [0, r-l+1], 0 correspond to l and r correspond to r-l+1
void temporal_graph::prefix_sum(std::vector <std::vector <int>> &temporal_data, int l, int r)
{
    std::vector <std::vector <int>> w;

    std::vector <int> empty_vec;

    for(int i = 0; i < edge_num; i++) {
        w.emplace_back(empty_vec);
    }

    int t_len = r - l + 1;

    //calc everytime's edge weight

    for (int i = l; i <= r; i++) {
        int upper_index = temporal_data[i].size();

        for (int j = 0; j < upper_index; j++) {
            int edge_index = temporal_data[i][j];
            int cur_size = w[edge_index].size();

            for (int k = cur_size; k < i - l; k++) {
                w[edge_index].push_back(-1);
            }
            w[edge_index].push_back(1);
        }
    }

    for (int i = 0; i < edge_num; i++) {
        int cur_size = w[i].size();
        for (int j = cur_size; j <= r - l; j++) {
            w[i].push_back(-1);
        }
    }

    //calc sum

    for (int i = 0; i < t_len; i++) {
        int sum = 0;
        for (int j = 0; j < edge_num; j++) {
            sum += w[j][i];
        }
        if (i > 0) {
            sum += all_edge_sum[i-1];
        }
        all_edge_sum.push_back(sum);
    }

    for (int i = 0; i < edge_num; i++) {
        std::vector <int> temp_sum;
        temp_sum.push_back(w[i][0]);
        for (int j = 1; j < t_len; j++) {
            temp_sum.push_back(temp_sum[j-1] + w[i][j]);
        }
        edge_sum.push_back(temp_sum);
    }
}

int temporal_graph::get_maxts()
{
    return max_ts;
}

int temporal_graph::get_vnum()
{
    return vertex_num;
}

//conhesive_density = all edge weight sum from time l to r
int temporal_graph::cohesive_density(int l, int r)
{
    if (l == 0)
        return all_edge_sum[r];
    return all_edge_sum[r] - all_edge_sum[l-1];
}

//positive_cohesive_density = all positive edge weight sum from time l to r
int temporal_graph::positive_cohesive_density(int l, int r)
{
    int psum = 0;
    for (int i = 0; i < edge_num; i++) {
        int sum = 0;
        if (l == 0) {
            sum = edge_sum[i][r];
        }   
        else {
            if (l - 1 < 0) {
                cout << l << " " << l - 1 << " here" << endl;
            }
            sum = edge_sum[i][r] - edge_sum[i][((int)l) - 1];
        }
        if (sum > 0)
            psum += sum;
    }

    return psum;
}

//convert temporal graph to aggregate graph
void temporal_graph::static_graph(graph_hash_of_mixed_weighted &g, std::vector <std::pair <int, int>> &e, int l, int r)
{
    std::vector <int> aggr_w;
    neg_edge_map.clear();
    pos_edge_map.clear();
    disjoint_set_union dsu = disjoint_set_union();

    //aggr_w[i] store edge i's aggregate weight from time l to time r
    for (int i = 0; i < edge_num; i++) {
        int sum_temp = edge_sum[i][r];
        if (l > 0)
            sum_temp -= edge_sum[i][l-1];
        aggr_w.push_back(sum_temp);
    }

    //positive edge connect two components
    for (int i = 0; i < edge_num; i++) {
        if(aggr_w[i] >= 0) {
            dsu._union(e[i].first, e[i].second);
        }
    }

    std::unordered_map <int, int> v_w;
    //every neg edge corresponds an edge between two vertex in new graph
    for (int i = 0; i < edge_num; i++) {
        if (aggr_w[i] >= 0) {
            int com_no = dsu.find(e[i].first);
            if (v_w.find(com_no) == v_w.end())
                v_w[com_no] = aggr_w[i];
            else
                v_w[com_no] = v_w[com_no] + aggr_w[i];
            
            pos_edge_map[i] = com_no;
        }
    }

    //add new graph's vertexs

    for (int i = 0; i < edge_num; i++) {
        int com_no1 = dsu.find(e[i].first), com_no2 = dsu.find(e[i].second);
        if (v_w.find(com_no1) == v_w.end())
            v_w[com_no1] = 0;
         if (v_w.find(com_no2) == v_w.end())
            v_w[com_no2] = 0;
    }

    for (std::unordered_map <int, int>::iterator it = v_w.begin(); it != v_w.end(); it++) {
        graph_hash_of_mixed_weighted_add_vertex(g, it->first, it->second);
        vertex_num = std::max(vertex_num, it->first + 1);
    }

    //add new graph's edge(neg edge in origin graph)

    for(int i = 0; i < edge_num; i++) {
        if (aggr_w[i] < 0) {
            int com_no1 = dsu.find(e[i].first), com_no2 = dsu.find(e[i].second);
            if (com_no1 == com_no2)
                continue;
            int old_weight = (int)graph_hash_of_mixed_weighted_edge_weight(g, com_no1, com_no2), new_weight = abs(aggr_w[i]);
            if (!graph_hash_of_mixed_weighted_contain_edge(g, com_no1, com_no2)) {
                graph_hash_of_mixed_weighted_add_edge(g, com_no1, com_no2, new_weight);
                neg_edge_map[std::make_pair(com_no1, com_no2)] = i;
                neg_edge_map[std::make_pair(com_no2, com_no1)] = i;
            }
            else {
                if (old_weight > new_weight) {
                    graph_hash_of_mixed_weighted_add_edge(g, com_no1, com_no2, new_weight);
                    neg_edge_map[std::make_pair(com_no1, com_no2)] = i;
                    neg_edge_map[std::make_pair(com_no2, com_no1)] = i;
                }
            }
        }
    }
}

//use infor record in neg_edge_map and pos_edge_map to convert aggregate graph to orgin graph
 void temporal_graph::origin_graph(std::unordered_map<int, int> &g, graph_hash_of_mixed_weighted &origin_g, std::vector <std::pair <int, int>> &e)
 {
     std::set <int> contain_vertex, origin_vertex;

     std::vector <std::pair<int, int>> origin_g_vec;

     for (auto it = g.begin(); it != g.end(); it++) {
        contain_vertex.insert(it->first);
        contain_vertex.insert(it->second);
        if (it->first == it->second)
            continue;
        
        int edge_index = neg_edge_map[std::make_pair(it->first, it->second)];
        origin_g_vec.push_back(std::make_pair(e[edge_index].first , e[edge_index].second));
        origin_vertex.insert(e[edge_index].first);
        origin_vertex.insert(e[edge_index].second);
     }

    for (int i = 0; i < edge_num; i++) {
        int com_no =  pos_edge_map[i];
        if (contain_vertex.find(com_no) != contain_vertex.end()) {
            origin_g_vec.push_back(std::make_pair(e[i].first, e[i].second));
            origin_vertex.insert(e[i].first);
            origin_vertex.insert(e[i].second);
        }
    }

    for (auto it = origin_vertex.begin(); it != origin_vertex.end(); it++) {
        graph_hash_of_mixed_weighted_add_vertex(origin_g, *it, 0);
    }

    for (int i = 0; i < origin_g_vec.size(); i++) {
        graph_hash_of_mixed_weighted_add_edge(origin_g, origin_g_vec[i].first, origin_g_vec[i].second, 0);
    }
 }