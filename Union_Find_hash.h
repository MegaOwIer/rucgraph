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

class Union_Find_hash_element_info
{
public:
	int parent;
	int set_size;

};


int Union_Find_hash_find(int v, unordered_map<int, Union_Find_hash_element_info>& Union_Find_hash_elements) {

	/*this function finds the set_ID of v via Path compression:
	https://en.wikipedia.org/wiki/Disjoint-set_data_structure */

	int root = v;
	while (root != Union_Find_hash_elements[root].parent)
		root = Union_Find_hash_elements[root].parent;
	while (v != root) {
		int newv = Union_Find_hash_elements[v].parent;
		Union_Find_hash_elements[v].parent = root;
		v = newv;
	}
	return root;
}

void Union_Find_hash_Union(int x, int y, unordered_map<int, Union_Find_hash_element_info>& Union_Find_hash_elements, int& sets_num) {

	/*this function Replace sets containing x and y with their union via the Union by size method:
	https://en.wikipedia.org/wiki/Disjoint-set_data_structure */

	int i = Union_Find_hash_find(x, Union_Find_hash_elements);
	int j = Union_Find_hash_find(y, Union_Find_hash_elements);
	if (i == j) return;

	// make smaller root point to larger one
	if (Union_Find_hash_elements[i].set_size < Union_Find_hash_elements[j].set_size) {
		Union_Find_hash_elements[i].parent = j; // smaller set merge to larger set; the merge method is used in Fast_GW_Growing_idealvertexID
		Union_Find_hash_elements[j].set_size += Union_Find_hash_elements[i].set_size;
	}
	else {
		Union_Find_hash_elements[j].parent = i;
		Union_Find_hash_elements[i].set_size += Union_Find_hash_elements[j].set_size;
	}
	sets_num--;
}


bool Union_Find_hash_in_the_same_Union(int x, int y, unordered_map<int, Union_Find_hash_element_info>& Union_Find_hash_elements) {
	return Union_Find_hash_find(x, Union_Find_hash_elements) == Union_Find_hash_find(y, Union_Find_hash_elements);
}


void try_Union_hash_Find() {

	/*below is an example of the initialization of all elements*/
	unordered_map<int, Union_Find_hash_element_info> Union_Find_hash_elements;
	for (int i = 0; i < 5; i++) {
		Union_Find_hash_element_info x;
		x.parent = i;
		x.set_size = 1;
		Union_Find_hash_elements[i] = x;
	}
	int sets_num = Union_Find_hash_elements.size();

	Union_Find_hash_Union(1, 3, Union_Find_hash_elements, sets_num);

	cout << "Union_Find_hash_find(1, Union_Find_hash_elements): " << Union_Find_hash_find(1, Union_Find_hash_elements) << endl;
	cout << "Union_Find_hash_find(2, Union_Find_hash_elements): " << Union_Find_hash_find(2, Union_Find_hash_elements) << endl;
	cout << "Union_Find_hash_in_the_same_Union(1, 2, Union_Find_hash_elements): " << Union_Find_hash_in_the_same_Union(1, 2, Union_Find_hash_elements) << endl;

}

