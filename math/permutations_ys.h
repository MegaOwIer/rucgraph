#pragma once


#include <iostream>
#include <vector>
#include <numeric>
#include <cstdint>
#include <cassert>
#include <math/Combinations_Permutations.h>


// functor called for each permutation
class permutations_ys_function
{
    unsigned len;
    std::uint64_t count;
    vector<vector<int>> permutation_elements;
public:

    explicit permutations_ys_function(unsigned l) : len(l), count(0) {}

    template <class It>
    bool operator()(It first, It last)  // called for each permutation
    {
        // count the number of times this is called
        ++count;

        vector<int> this_permutation_elements;
        for (auto it = first; it != last; it++) {
            this_permutation_elements.push_back(*it);
        }
        permutation_elements.push_back(this_permutation_elements);

        return false;  // Don't break out of the loop
    }

    operator std::uint64_t() const { return count; }

    vector<vector<int>> GetVect()
    {
        return { permutation_elements };
    }
};

void example_permutations_ys() {

    const int r = 2;
    const int n = 3;

    std::vector<int> v = { 4,2,6 };

    /*for_each_reversible_circular_permutation(v.begin(), v.begin() + r, v.end(), f(v.size())) return a f class;

    see #include <Combinations_Permutations.h>:

    template <class BidirIter, class Function>
Function
for_each_reversible_circular_permutation(BidirIter first,
    BidirIter mid,
    BidirIter last, Function f)
{
    for_each_combination(first, mid, last, detail::reverse_circular_permutation<Function&,
        BidirIter>(f, std::distance(first, mid)));
    return f;
}*/
    std::vector<vector<int>> permutation_elements = for_each_reversible_circular_permutation(v.begin(), v.begin() + r, v.end(), permutations_ys_function(v.size())).GetVect();

    cout << "permutation elements: " << endl; // from first to last
    for (int i = 0; i < permutation_elements.size(); i++) {
        for (int j = 0; j < permutation_elements[i].size(); j++) {
            cout << permutation_elements[i][j];
            if (j != permutation_elements[i].size() - 1) {
                cout << ", ";
            }
        }
        cout << endl;
    }
    std::cout << "Found " << permutation_elements.size() << " permutations of " << v.size() << " objects taken " << r << " at a time.\n";
}