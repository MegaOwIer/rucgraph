#pragma once

/*these codes are for a disjoint-set data structure (also called a union¨Cfind data structure or merge¨Cfind set)
is a data structure that tracks a set of elements partitioned into a number of disjoint (non-overlapping) subsets.:
https://en.wikipedia.org/wiki/Disjoint-set_data_structure

e.g., 0,1,2,...,N are initially all disjoint-sets;
1 may be merged into set 0, etc. this function keeps track of which set is a number in, how many sets,
and also merge sets.

time complexity for all operations below: O(alpha(n)) approx O(1):
https://en.wikipedia.org/wiki/Disjoint-set_data_structure
Using both path compression, splitting, or halving and union by rank or size ensures that
the amortized time per operation is only {\displaystyle O(\alpha (n))}O(\alpha (n)),[4][5]
which is optimal,[6] where {\displaystyle \alpha (n)}\alpha (n) is the inverse Ackermann
function. This function has a value {\displaystyle \alpha (n)<5}{\displaystyle \alpha (n)<5} for any value of n
that can be written in this physical universe, so the disjoint-set operations take
place in essentially constant time. (\alpha (n) is smaller than log n)

These codes are copied from
https://stackoverflow.com/questions/8300125/union-find-data-structure
https://github.com/kartikkukreja/blog-codes/blob/master/src/Union%20Find%20%28Disjoint%20Set%29%20Data%20Structure.cpp#L46
*/

#include<vector>

void Union_Find_initialization(vector<int>& parent, vector<int>& set_size) { // pass arrays

	/*we assume that parent, and set_size are arrays with the same size*/

	int N = parent.size();

	for (int i = 0; i < N; i++) {
		parent[i] = i; // Create an empty union find data structure with N isolated sets.
		set_size[i] = 1;
	}
}

int Union_Find_find(int v, vector<int>& parent) {

	/*this function finds the set_ID of v via Path compression:
	https://en.wikipedia.org/wiki/Disjoint-set_data_structure */

	int root = v;
	while (root != parent[root])
		root = parent[root];
	while (v != root) {
		int newv = parent[v];
		parent[v] = root;
		v = newv;
	}
	return root;
}

void Union_Find_Union(int x, int y, vector<int>& parent, vector<int>& set_size, int& sets_num) {

	/*this function Replace sets containing x and y with their union via the Union by size method:
	https://en.wikipedia.org/wiki/Disjoint-set_data_structure */

	int i = Union_Find_find(x, parent);
	int j = Union_Find_find(y, parent);
	if (i == j) return;

	// make smaller root point to larger one
	if (set_size[i] < set_size[j]) {
		parent[i] = j; // smaller set merge to larger set; the merge method is used in Fast_GW_Growing_idealvertexID
		set_size[j] += set_size[i];
	}
	else {
		parent[j] = i;
		set_size[i] += set_size[j];
	}
	sets_num--;
}


bool Union_Find_in_the_same_Union(int x, int y, vector<int>& parent) {
	return Union_Find_find(x, parent) == Union_Find_find(y, parent);
}



/*below is an example*/

#include<iostream>
void example_Union_Find() {

	int sets_num = 5;
	vector<int> parent(sets_num), set_size(sets_num);

	Union_Find_initialization(parent, set_size);

	Union_Find_Union(1, 2, parent, set_size, sets_num);

	cout << "Union_Find_find(1, input_array, parent): " << Union_Find_find(1, parent) << endl;
	cout << "Union_Find_find(2, input_array, parent): " << Union_Find_find(2, parent) << endl;
	cout << "Union_Find_in_the_same_Union(1, 2, input_array, parent): "
		<< Union_Find_in_the_same_Union(1, 2, parent) << endl;

}

