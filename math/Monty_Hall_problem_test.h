#pragma once

/* https://en.wikipedia.org/wiki/Monty_Hall_problem */
#include <boost/random.hpp>
void Monty_Hall_problem_test() {

    std::time_t now = std::time(0);
    boost::random::mt19937 gen{ static_cast<std::uint32_t>(now) };

    int example_num = 1e6;

    int not_switch_win_times = 0;
    int switch_win_times = 0;

    for (int i = 0; i < example_num; i++) {

        vector<bool> doors(3, false); // false is goat
        boost::random::uniform_int_distribution<> dist_1st_0_2{ 0, 2 };
        doors[dist_1st_0_2(gen)] = true;

        boost::random::uniform_int_distribution<> dist_2nd_0_2{ 0, 2 };
        int first_choice = dist_2nd_0_2(gen);

        vector<int> revealed_goats_candicates;
        for (int j = 0; j < 3; j++) {
            if (j != first_choice && doors[j] == false) {
                revealed_goats_candicates.push_back(j);
            }
        }
        boost::random::uniform_int_distribution<> dist{ 0, (int)revealed_goats_candicates.size() - 1 };
        int reveal_id = revealed_goats_candicates[dist(gen)];

        int switch_id;
        for (int j = 0; j < 3; j++) {
            if (j != first_choice && j != reveal_id) {
                switch_id = j;
                break;
            }
        }

        if (doors[first_choice]) { // always 0.33 probability
            not_switch_win_times++;
        }
        else { // 0.67 probability
            switch_win_times++;
        }

    }

    cout << "not_switch_win_probability = " << (double)not_switch_win_times / example_num << endl;
    cout << "switch_win_probability = " << (double)switch_win_times / example_num << endl;
}