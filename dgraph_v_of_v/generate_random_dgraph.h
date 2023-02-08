#pragma once
#include "dgraph_v_of_v.h"
#include <dgraph_v_of_v/dgraph_v_of_v.h>
#include <boost/random.hpp>

/**
 * @brief this function generates a random directed graph , and this graph may not be connected.
 *
 * @param V : the number of vertices
 * @param E : the number of edges
 * @param ec_min : minimum value of edge weight
 * @param ec_max : maximum value of edge weight
 * @param input_precision : decimal
 * @param boost_random_time_seed : random seed
 * @return dgraph_v_of_v<double>
 */

dgraph_v_of_v<double> generate_random_dgraph(long long int V, long long int E, double ec_min, double ec_max,
                                             int input_precision, boost::random::mt19937 &boost_random_time_seed)
{
    /*time complexity: O(|E|)*/

    double precision = std::pow(10, input_precision);
    boost::random::uniform_int_distribution<> dist_ec{static_cast<int>(ec_min * precision),
                                                      static_cast<int>(ec_max * precision)};

    dgraph_v_of_v<double> random_graph(V);

    long long int max_E = V * (V - 1); // edge number of dgraph should double

    if (E == max_E)
    { // complete graphs
        /*time complexity: O(|V|^2)*/
        for (int i = 0; i < V; i++)
        {
            for (int j = 0; j < V; j++)
            {
                if (i != j)
                {
                    double new_cost = (double)dist_ec(boost_random_time_seed) / precision;
                    random_graph.add_edge(i, j, new_cost);
                }
            }
        }
    }
    else if (E > max_E)
    {
        std::cout << "E: " << E << std::endl;
        std::cout << "V * (V - 1): " << max_E << std::endl;
        std::cout << "E > V * (V - 1) in generate_random_dgraph!" << '\n';
        exit(1);
    }
    else
    { // incomplete graphs

        /*time complexity: O(|V|)*/
        std::vector<int> not_full_vertices; // vertices without a full degree
        for (int i = 0; i < V; i++)
        {
            not_full_vertices.insert(not_full_vertices.end(), i);
        }

        /*time complexity: O(|V||E|)*/
        int edge_num = 0;
        while (edge_num < E)
        {
            boost::random::uniform_int_distribution<> dist_id{static_cast<int>(0),
                                                              static_cast<int>(not_full_vertices.size() - 1)};
            int RAND = dist_id(boost_random_time_seed); // generate int random number  0, not_full_vertices.size()-1
            if (random_graph.degree_out(not_full_vertices[RAND]) < V - 1)
            { // here, we only need to check the out degree; later the 'contain_edge' will check the in degree
                /*time complexity: O(|V|)*/
                std::vector<int> unchecked(V);
                std::iota(std::begin(unchecked), std::end(unchecked), 0);
                bool added = false;
                while (added == false)
                {
                    boost::random::uniform_int_distribution<> dist_id2{static_cast<int>(0),
                                                                       static_cast<int>(unchecked.size() - 1)};
                    int x = dist_id2(boost_random_time_seed);
                    int j = unchecked[x];
                    if (not_full_vertices[RAND] != j && random_graph.contain_edge(not_full_vertices[RAND], j) == 0)
                    { // This edge does not exist
                        double new_cost = (double)dist_ec(boost_random_time_seed) / precision;
                        random_graph.add_edge(not_full_vertices[RAND], j, new_cost);
                        edge_num++;
                        added = true;
                        break;
                    }
                    else
                    {
                        unchecked.erase(unchecked.begin() + x);
                    }
                }
            }
            else
            { // this is a vertex with a full degree
                not_full_vertices.erase(not_full_vertices.begin() + RAND);
            }
        }
    }

    return random_graph;
}



void example_generate_random_dgraph()
{
    int V = 7, E = 20;
    double ec_min = 0.1, ec_max = 1;
    int precision = 1;
    
    dgraph_v_of_v<double> instance_graph = generate_random_dgraph(V, E, ec_min, ec_max, precision, boost_random_time_seed);
}