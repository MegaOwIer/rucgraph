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
#include <ctime>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_Fast_GW_Growing_tree.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_General_Pruning_tree.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_minimum_spanning_tree.h>
#include <graph_hash_of_mixed_weighted/graph_hash_of_mixed_weighted_extract_subgraph_for_a_hash_of_vertices.h>
#include "temporal_graph.h"


#define NINF -2e9

//#define PCST_DEBUG

struct maxima
{
    int x, l, u;
};

struct time_interval
{
    int l, r;
    bool operator < (const time_interval & b) const {
        if (l == b.l)
            return r < b.r;
        return l < b.l;
    }
};

struct simple_path
{
    std::vector <int> vertex_list;
    int networth;
    bool operator < (const simple_path & b) const {
        return networth < b.networth;
    }
};

struct bfs_node
{
    int vertex, deep,nw;
    std::vector <int> v_list;
};

//calc everytime's cohesive density and smooth factor, smooth factor equal to alpha * average change of conhesive density 

void cohesive_density_curve(temporal_graph &g, std::vector <double> &f, int maxt, double alpha, double &smooth_fac)
{
    double cden_sum = 0;
    for (int i = 0; i < maxt; i++) {
        int cd = g.cohesive_density(i, i);
        f.push_back(cd);
        if (i > 0) {
            cden_sum += cd - f[i-1];
        }
    }

    if (maxt > 1)
        cden_sum /= (maxt - 1);

    smooth_fac = cden_sum * alpha;

    f.push_back(NINF);
}

//find local max and their corresponding interval
void find_local_max(std::vector <double> &f, std::vector <time_interval> &lm, int maxt, double smooth_fac)
{
    


    //if max timestamp <= 10, all time interval is no more than 100, return all time interval
    if (maxt <= 10) {
        for (int i = 0; i < maxt; i++) {
            time_interval t;
            t.l = i;
            t.r = i;
            lm.push_back(t);
        }
        return;
    }

    std::vector <maxima> local_max;
    //find all local max
    for (int i = 0; i < maxt; i++) {
        if (i > 0 && f[i] < f[i - 1])
            continue;
        if (f[i] < f[i + 1])
            continue;
        maxima x;
        x.l = 0;
        x.u = 0;
        x.x = i;
        local_max.push_back(x);
    }

    //find local max's interval
    int h = local_max.size();
    for (int i = 0; i < h; i++) {

        
        for (int j = local_max[i].x + 1; ; j++) {
            if (i != h - 1) {
                if (j >= local_max[i+1].x) {
                    local_max[i].u = std::min(j, maxt - 1);
                    break;
                }
            }
            if (j >= maxt - 1) {
                local_max[i].u = std::min(j, maxt - 1);
                break;
            }
            if (fabs(f[j] - f[local_max[i].x]) >= smooth_fac) {
                local_max[i].u = j;
                break;
            }
        }
        for (int j = std::max(local_max[i].x - 1, 0); ;j--) {
            if (j <= 0 && i == 0) {
                local_max[i].l = j;
                break;
            }
            //merge two overlap intervals
            if (i > 0) {
                if (j <= 0 || j <= local_max[i - 1].u) {
                    local_max[i].l = local_max[i - 1].l;
                    double val_1 = f[local_max[i - 1].x], val_2 = f[local_max[i].x];
                    if (val_1 > val_2)
                        local_max[i].x = local_max[i - 1].x;
                    break;
                }
            }
            if (fabs(f[j] - f[local_max[i].x]) >= smooth_fac) {
                local_max[i].l = j;
                break;
            }
        }

    }

    for (int i = 0; i < h - 1; i++) {
        if (local_max[i].l == local_max[i+1].l)
            continue;
        time_interval x1;
        x1.l = local_max[i].l;
        x1.r = local_max[i].u;
        lm.push_back(x1);
    }
    time_interval x2;
    x2.l = local_max[h-1].l;
    x2.r = local_max[h-1].u;
    lm.push_back(x2);

}

//find local min and their corresponding interval
void find_local_min(std::vector <double> &f, std::vector <time_interval> &lm, int maxt, double smooth_fac)
{

    if (maxt <= 10) {
        for (int i = 0; i < maxt; i++) {
            time_interval t;
            t.l = i;
            t.r = i;
            lm.push_back(t);
        }
        return;
    }

    std::vector <maxima> local_min;

    for (int i = 0; i < maxt; i++) {
        if (i > 0 && f[i] > f[i - 1])
            continue;
        if (f[i] > f[i + 1] && i != maxt - 1)
            continue;
        maxima x;
        x.l = 0;
        x.u = 0;
        x.x = i;
        local_min.push_back(x);
    }

    int h = local_min.size();
    for (int i = 0; i < h; i++) {

        for (int j = local_min[i].x + 1; ; j++) {
            if (i != h - 1) {
                if (j >= local_min[i+1].x) {
                    local_min[i].u = std::min(j, maxt - 1);
                    break;
                }
            }
            if (j >= maxt - 1) {
                local_min[i].u = std::min(j, maxt - 1);
                break;
            }
            if (fabs(f[j] - f[local_min[i].x]) >= smooth_fac) {
                local_min[i].u = j;
                break;
            }
        }
        for (int j = std::max(local_min[i].x - 1, 0); ;j--) {
            if (j <= 0 && i == 0) {
                local_min[i].l = j;
                break;
            }
            if (i > 0) {
                if (j <= 0 || j <= local_min[i - 1].u) {
                    local_min[i].l = local_min[i - 1].l;
                    double val_1 = f[local_min[i - 1].x], val_2 = f[local_min[i].x];
                    if (val_1 > val_2)
                        local_min[i].x = local_min[i - 1].x;
                    break;
                }
            }
            if (fabs(f[j] - f[local_min[i].x]) >= smooth_fac) {
                local_min[i].l = j;
                break;
            }
        }
    }

    for (int i = 0; i < h - 1; i++) {
        if (local_min[i].l == local_min[i+1].l)
            continue;
        time_interval x1;
        x1.l = local_min[i].l;
        x1.r = local_min[i].u;
        lm.push_back(x1);
    }
    time_interval x2;
    x2.l = local_min[h-1].l;
    x2.r = local_min[h-1].u;
    lm.push_back(x2);

}

//generate k/2 interval has top k pcd with local max and local min data
void calc_s_k2(temporal_graph &g, std::vector <time_interval> &lm, int k, std::vector <time_interval> &s2k)
{   
    std::vector <std::pair<int, time_interval>> s;
    int h = lm.size();
    for (int i = 0; i < h ; i++) {
        for (int j = i; j < h; j++) {
            int pcd = g.positive_cohesive_density(lm[i].l, lm[j].r);
            time_interval x;
            x.l = lm[i].l;
            x.r = lm[j].r;
            s.push_back(std::make_pair(pcd, x));
        }
    }
    
    h = s.size();
    for (int i = 0; i < k; i++) {
        int max_val = -1, max_x;
        for (int j = 0; j < h; j++) {
            if (s[j].first > max_val) {
                max_val = s[j].first;
                max_x = j;
            }
        }
        if (max_val != -1) {
            s2k.push_back(s[max_x].second);
            s[max_x].first = -1;
        }
    }
}

//combine interval generate from local max and local min
void combine_k2(std::vector <time_interval> &lm, std::vector <time_interval> &lm_1, std::vector <time_interval> &lm_2) {
    for (int i = 0; i < lm_1.size(); i++) {
        lm.push_back(lm_1[i]);
    }

    for (int i = 0; i < lm_2.size(); i++) {
        lm.push_back(lm_2[i]);
    }
}

//try to optimize interval use adding-doubling algorithm
void time_interval_tune(temporal_graph &g, std::vector <time_interval> &s2k, int k, int maxt)
{
    for (int i = 0; i < s2k.size(); i++) {
        int base = 1;
        int ctime = 0;
        int pcd = g.positive_cohesive_density(s2k[i].l, s2k[i].r);

        while (1) {
            if (base > maxt)
                break;
            bool flag =false;

            int newl = std::max(0, s2k[i].l - base);
            int new_pcd = g.positive_cohesive_density(newl, s2k[i].r);

            if (new_pcd > pcd) {
                s2k[i].l = newl;
                pcd = new_pcd;
                flag = true;
            }

            int newr = std::min(maxt - 1, s2k[i].r + base);
            new_pcd = g.positive_cohesive_density(s2k[i].l, newr);

            if (new_pcd > pcd) {
                s2k[i].r = newr;
                pcd = new_pcd;
                flag = true;
            }

            //double base after every four steps
            if (!flag) break;
            ctime++;
            if (ctime == 4) {
                base *= 2;
                ctime = 0;
            }
        }
    }

}

bool sp_cmp(const simple_path & a, const simple_path & b) {
    return a.networth > b.networth;
}

bool sort_vertex_cmp(std::pair <int, int> &a, std::pair <int, int> &b){
    return a.second > b.second;
}

void bounded_probing(graph_hash_of_mixed_weighted &st, graph_hash_of_mixed_weighted &h, int vertex_num)
{
    int max_dep = 3;

    std::vector <bool> in_solution;
    
    for (int i = 0; i < vertex_num; i++) {
        if(graph_hash_of_mixed_weighted_contain_vertex(st, i)) {
            in_solution.push_back(true);
        }
        else {
            in_solution.push_back(false);
        }
    }

    std::vector <bool> previs, vis;
    std::vector <std::vector<std::pair <int, int>>> prepath, path, bestpath;
    std::vector <int> pregain, gain, bestgain;

    for (int i = 0; i < vertex_num; i++) {
        previs.push_back(false);
        vis.push_back(false);
        gain.push_back(-1);
        pregain.push_back(-1);
        bestgain.push_back(-1);
        std::vector<std::pair <int, int>> temp_vec;
        for (int k = 0; k < max_dep; k++) {
            temp_vec.push_back(std::make_pair(-1, -1));
        }
        prepath.push_back(temp_vec);
        path.push_back(temp_vec);
        bestpath.push_back(temp_vec);
    }
    for (int ii = 0 ; ii < max_dep; ii++) {
        for (int i = 0; i < vertex_num; i++) {
            bestgain[i] = -1;
            previs[i] = false;
            if (in_solution[i]) {
                previs[i] = true;
                pregain[i] = 0;
            }
        }

        for (int l = 0; l < max_dep; l++) {
            for(int i = 0; i < vertex_num; i++) {
                gain[i] = -1;
                vis[i] = false;
            }

            for (int i = 0; i < vertex_num; i++) {
                if (!previs[i])
                    continue;
                std::vector <int> adj_list;
                adj_list = h.adj_v(i);    
                for (int j = 0; j < adj_list.size(); j++) {
                    int v = adj_list[j];
                    if (in_solution[v])
                        continue;

                    bool overlap_flag = false;
                    for (int k = 0; k < l; k++) {
                        if(prepath[i][k].second == v)
                            overlap_flag = true;
                    }
                    if (overlap_flag)
                        continue;
                    
                    int new_val  = pregain[i] + h.hash_of_vectors[v].vertex_weight - (int)graph_hash_of_mixed_weighted_edge_weight(h, i, v);
                    if (new_val > gain[v]) {
                        vis[v] = true;
                        gain[v] = new_val;
                        for (int k = 0; k < l; k++) {
                            path[v][k] = prepath[i][k];
                        }
                        path[v][l] = std::make_pair(i, v);
                        if(gain[v] > bestgain[v]) {
                            bestgain[v] = gain[v];
                            for (int k = 0; k <= l; k++) {
                                bestpath[v][k] = path[v][k];
                            }
                            for (int k = l+1; k < max_dep; k++) {
                                bestpath[v][k] = std::make_pair(-1, -1);
                            }
                        }
                    }
                }
            }
            std::vector <bool> tempvis;
            for (int i = 0; i < vertex_num; i++) {
                tempvis.push_back(previs[i]);
                previs[i] = vis[i];
                vis[i] = tempvis[i];
            }

            for (int i = 0; i < vertex_num; i++) {
                std::vector <std::pair<int, int>> temp_vec;
                for (int k = 0; k < max_dep; k++) {
                    temp_vec.push_back(prepath[i][k]);
                    prepath[i][k] = path[i][k];
                    path[i][k] = temp_vec[k];
                }
            }

            std::vector <int> tempgain;
            for (int i = 0; i < vertex_num; i++) {
                tempgain.push_back(pregain[i]);
                pregain[i] = gain[i];
                gain[i] = tempgain[i];
            }

        }

        std::vector <std::pair<int, int>> sorted_vertex;
        for (int i = 0; i < vertex_num; i++) {
            if (!in_solution[i] && bestgain[i] > 0 && h.hash_of_vectors[i].vertex_weight > 0) {
                sorted_vertex.push_back(std::make_pair(i, bestgain[i]));
            }
        }

        std::sort(sorted_vertex.begin(), sorted_vertex.end(), sort_vertex_cmp);

        bool change = false;
        for (int i = 0; i < sorted_vertex.size(); i++) {
            bool overlap = false;
            int t = sorted_vertex[i].first;
            for (int k = 0; k < max_dep; k++) {
                if (bestpath[t][k].first == -1)
                    break;
                if (in_solution[bestpath[t][k].second]) {
                    overlap = true;
                    break;
                }
            }

            if (!overlap) {
                for (int k = 0; k < max_dep; k++) {
                    if (bestpath[t][k].first == -1)
                        break;
                    in_solution[bestpath[t][k].second] = true;
                    int v  = bestpath[t][k].second;
                    graph_hash_of_mixed_weighted_add_vertex(st, v, h.hash_of_vectors[v].vertex_weight);
                    graph_hash_of_mixed_weighted_add_edge(st, bestpath[t][k].first, bestpath[t][k].second, graph_hash_of_mixed_weighted_edge_weight(h, bestpath[t][k].first, bestpath[t][k].second));
                    
                }
                change = true;

            }
        }

        if (change == false)
            break;
    }

}

int calc_pcst_result(graph_hash_of_mixed_weighted &g, std::unordered_map <int, int> mst) {
    int sum = 0;
    for (auto it = mst.begin(); it != mst.end(); it++) {
        if (it->first == it->second)
            continue;
        //std::cout << it->first<<" " <<it->second << " "<<graph_hash_of_mixed_weighted_edge_weight(g, it->first, it->second) <<std::endl;
        sum -= graph_hash_of_mixed_weighted_edge_weight(g, it->first, it->second);
    }

    //std::cout << "edge sum:" <<sum <<std::endl;

    std::set <int> vertex_set;
    for (auto it = mst.begin(); it != mst.end(); it++) {
        if (vertex_set.find(it->first) == vertex_set.end()) {
            sum += (int)g.hash_of_vectors[it->first].vertex_weight;
            vertex_set.insert(it->first);
        }
        if (vertex_set.find(it->second) == vertex_set.end()) {
            sum += (int)g.hash_of_vectors[it->second].vertex_weight;
            vertex_set.insert(it->second);
        }
        
    }
    //std::cout << "all sum:" <<sum <<std::endl;
    return sum;
}

void map_to_graph(graph_hash_of_mixed_weighted &origin_g, graph_hash_of_mixed_weighted &g, std::unordered_map <int, int> &mp)
{
    std::set <int> vertex_set;
    std::unordered_map <int, int> vertex_weight;
    for (auto it = mp.begin(); it != mp.end(); it++) {
        vertex_set.insert(it->first);
        vertex_weight[it->first] = origin_g.hash_of_vectors[it->first].vertex_weight;
        vertex_set.insert(it->second);
        vertex_weight[it->second] = origin_g.hash_of_vectors[it->second].vertex_weight;
    }

    for (auto it = vertex_set.begin(); it != vertex_set.end(); it++) {
        graph_hash_of_mixed_weighted_add_vertex(g, *it, vertex_weight[*it]);
    }

    for (auto it = mp.begin(); it != mp.end(); it++) {
        if (it->first == it->second)
            continue;
        int w = graph_hash_of_mixed_weighted_edge_weight(origin_g, it->first, it->second);
        graph_hash_of_mixed_weighted_add_edge(g, it->first, it->second, w);
    }
    
}

void debug_print_all_pcd_sorted(temporal_graph &g)
{
    int maxt = g.get_maxts();
    std::vector <std::pair <int, time_interval>> pcd;
    for (int i = 0; i < maxt; i++) {
        for (int j = i; j < maxt; j++) {
            time_interval x;
            x.l = i;
            x.r = j;
            pcd.push_back(std::make_pair(g.positive_cohesive_density(i, j), x) );
        }
    }
    std::sort(pcd.begin(), pcd.end());
    for (int i = 0; i < pcd.size(); i++) {
        std::cout << "l:" << pcd[i].second.l << " r:" << pcd[i].second.r << " pcd:" << pcd[i].first<< "\n";
    }
    std::cout << std::endl;
}

void time_interval_selector (temporal_graph &g, int k, double alpha, std::vector <time_interval> & sk)
{
    std::vector <double> f;
    std::vector <time_interval> lmax, sk2_max, sk2_min, lmin;


    double smooth_fac;

    int maxt = g.get_maxts();

    cohesive_density_curve(g, f, maxt, alpha, smooth_fac);
    find_local_max(f, lmax, maxt, smooth_fac);
    find_local_min(f, lmin, maxt, smooth_fac);
    calc_s_k2(g, lmax, (k+1)/2, sk2_max);
    calc_s_k2(g, lmin, k/2, sk2_min);

    combine_k2(sk, sk2_max, sk2_min);
    time_interval_tune(g, sk, k, maxt);
}

std::pair < int, std::pair <std::pair <int, int>, graph_hash_of_mixed_weighted>> fides(std::vector <std::pair <int, int>> &graph_data,  std::vector <std::vector <int>> &temporal_data, int time_interval_l, int time_interval_u,  
                                                    int k, double alpha, int option)
{

    clock_t start_time = clock(), cur_time;

    temporal_graph g = temporal_graph(graph_data, time_interval_l, time_interval_u);
    g.prefix_sum(temporal_data, time_interval_l, time_interval_u);

    std::vector <time_interval> sk2;

    time_interval_selector(g, k, alpha, sk2);

    #ifdef PRINT_RUN_TIME
        cur_time = clock();
        std::cout << "time_interval use: " << cur_time - start_time << std::endl;
        start_time = cur_time;
    #endif

    int pcst_result = NINF;
    std::unordered_map<int, int> result_graph;
    std::pair <int, int> result_interval;
    
    for (int i = 0; i < sk2.size(); i++) {
        //std::cout<<sk2[i].l<<" "<<sk2[i].r<<std::endl;
        graph_hash_of_mixed_weighted aggr_g;
        g.static_graph(aggr_g, graph_data, sk2[i].l, sk2[i].r);
        #ifdef PRINT_RUN_TIME
            cur_time = clock();
            std::cout << "static graph use: " << cur_time - start_time << std::endl;
            start_time = cur_time;
        #endif
        //run pcst with aggr_g here
        graph_hash_of_mixed_weighted mst_g, pruning_tree;
        if (option == 1) {
            mst_g = graph_hash_of_mixed_weighted_Fast_GW_Growing_tree(aggr_g, "tree");;
        }
        else {
            std::unordered_map <int, int> mst_g_map = graph_hash_of_mixed_weighted_minimum_spanning_tree(aggr_g);
            map_to_graph(aggr_g, mst_g, mst_g_map);
        }

        #ifdef PRINT_RUN_TIME
            cur_time = clock();
            std::cout << "spanning tree use: " << cur_time - start_time << std::endl;
            start_time = cur_time;
        #endif 
        pruning_tree = graph_hash_of_mixed_weighted_General_Pruning_tree(mst_g);
        #ifdef PRINT_RUN_TIME
            cur_time = clock();
            std::cout << "pruning tree use: " << cur_time - start_time << std::endl;
            start_time = cur_time;
        #endif 
        bounded_probing(pruning_tree, aggr_g, g.get_vnum());
        #ifdef PRINT_RUN_TIME
            cur_time = clock();
            std::cout << "bounded_probing use: " << cur_time - start_time << std::endl;
            start_time = cur_time;
        #endif
        std::unordered_map<int, int> final_mst;
        final_mst = graph_hash_of_mixed_weighted_minimum_spanning_tree(pruning_tree);
        int cur_pcst_result = calc_pcst_result(pruning_tree, final_mst);
        if (cur_pcst_result > pcst_result) {
            pcst_result = cur_pcst_result;
            result_graph = final_mst;
            result_interval = std::make_pair(sk2[i].l, sk2[i].r);
        }
    }
    graph_hash_of_mixed_weighted origin_result_graph;

    g.origin_graph(result_graph, origin_result_graph, graph_data);


    //std::cout << "result interval:" << result_interval.first << "," << result_interval.second << std::endl;

    return std::make_pair(pcst_result, std::make_pair(result_interval, origin_result_graph));
    //debug_print_all_pcd_sorted(g);
}