#pragma once

/*the following codes are for testing

---------------------------------------------------
a cpp file (try.cpp) for running the following test code:
----------------------------------------

#include <iostream>
#include <fstream>
using namespace std;
#include <boost/random.hpp>
boost::random::mt19937 boost_random_time_seed{ static_cast<std::uint32_t>(std::time(0)) };

#include <dgraph_v_of_v/test_dgraph.h>

int main()
{
    test_dgraph();
}

------------------------------------------------------------------------------------------
Commends for running the above cpp file on Linux:

g++ -std=c++17 -I/home/boost_1_75_0 -I/home/dengs/dgraph try.cpp -lpthread -Ofast -o A
./A
rm A

(optional to put the above commends in run.sh, and then use the comment: sh run.sh)


*/

#include <dgraph_v_of_v/dgraph_v_of_v.h>
#include <dgraph_v_of_v/generate_random_dgraph.h>
#include <dgraph_v_of_v/save_dgraph_with_weight.h>
#include <dgraph_v_of_v/read_dgraph_with_weight.h>
void test_dgraph()
{
    int V = 7, E = 20;
    int thread_num = 5;
    double ec_min = 0.1, ec_max = 1;
    int precision = 1;
    
    dgraph_v_of_v<double> instance_graph = generate_random_dgraph(V, E, ec_min, ec_max, precision, boost_random_time_seed);
    save_dgraph_with_weight("random_dgraph_test.txt", instance_graph);

    read_dgraph_with_weight("random_dgraph_test.txt", instance_graph);
    save_dgraph_with_weight("random_dgraph_test_copy.txt", instance_graph);
}