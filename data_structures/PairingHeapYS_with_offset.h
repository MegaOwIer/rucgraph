#pragma once
#include <vector>
#include<iostream>

/* PairingHeapYS_with_offset is specially designed the fast GW growing algorithm. In this heap,
there is an offset value for every node. Using these values, we can change the
key values of all nodes in time O(1) (only for conditions in the fast GW growing algorithm)

the offset values induce a lot of calculations in the operations of this heap;
as a result, the precision of key values in this heap is roughly around 1e-5
when n is large (>1e6) and ValueType is double


size of pointer: On 32-bit machine sizeof pointer is 32 bits ( 4 bytes), while on 64 bit machine it's 8 byte.
	most of the memory consumption of this heap is done by pointers
*/

template <typename ValueType, typename PayloadType> // input two data types
class PairingHeapYS_with_offset {

private:
	struct Node {
		Node* sibling; // pointer to other node; 
		Node* child;
		Node* left_up; // the adj node from this node to the root, i.e., the adj_root
		ValueType value; // key value with offset
		//ValueType true_key_value_without_offset; // it's not possible to efficiently store this value if you want add_to_heap() in O(1)
		ValueType child_offset;
		/*offset value for direct child (initial key value of grandchild is value
		of grandchild plus child_offset of child, plus child_offset of left_up, etc, to child_offset of root )*/
		PayloadType payload; // item
	};
	typedef pair<ValueType, PayloadType> pair_content;
	size_t node_num; // size_t has the largest range

	bool heap_type_compare(ValueType i, ValueType j)
	{
		return i < j; // this is a min heap
	}



public: // all below functions can be accessed outside
	typedef Node* ItemHandle; // ItemHandle is required to make the node pointer data type public
	Node* root = NULL;

private:

	/*the link function:  this is to insert new node or melt
	two heaps by two roots, see https://www.geeksforgeeks.org/pairing-heap/ :
	To join the two heap, first, we compare the root node of the heap;
	if the root node of the first heap is smaller than the root node of the second
	heap, then root node of the second heap becomes a left child of the
	root node of the first heap.
	time complexity O(1)*/
	Node* link(Node* node1, Node* node2) { // two nodes are two roots of two trees; roots have no leftup or sibling
		if (node1 == NULL) { // the linked root/heap may be empty
			return node2;
		}
		if (node2 == NULL) {
			return node1;
		}
		Node* top_node = node2;
		Node* below_node = node1;
		if (heap_type_compare(node1->value, node2->value)) {
			top_node = node1;
			below_node = node2;
		}
		below_node->sibling = top_node->child;
		below_node->left_up = top_node;
		top_node->child = below_node;
		top_node->sibling = NULL; // root->sibling = NULL
		top_node->left_up = NULL; // root->left_up = NULL
		below_node->value -= top_node->child_offset; // all below_node childs value should minus smaller_node child_offset
		below_node->child_offset -= top_node->child_offset; // update below_node->child_offset
		return top_node; // return the new root pointer
	}

	/*the two pass method link a vector of nodes (usually in the same layer) by pairs; time complexity O(log n)*/
	void two_pass_method(std::vector<ItemHandle>& buffer) {
		size_t merged_children = 0;
		while (merged_children + 2 <= buffer.size()) {
			buffer[merged_children / 2] = link(
				buffer[merged_children], buffer[merged_children + 1]); 
			merged_children += 2;
		}
		if (merged_children != buffer.size()) {
			buffer[merged_children / 2] = buffer[merged_children]; 
			buffer.resize(merged_children / 2 + 1); 
		}
		else {
			buffer.resize(merged_children / 2);
		}
	}


public:


	/* some knowledge:

	root->sibling = NULL;
	root->left_up = NULL;

	root = NULL, when heap is empty;

	sibling, child and left_up values are correctly and constantly maintained for all nodes

	sibling values are directed from left to right

	----------------------
	A (child_offset 1)
	|
	B (child_offset 2) - C (child_offset 3) - D (child_offset 4)
	|                    |
	E - F                G - H
	---------------------------
	An operation of puting B as child of A: 1. B->left_up = A; 2. A->child = B.

	*/

	PairingHeapYS_with_offset() { // initialize a heap with heap type; time complexity O(1)
		node_num = 0; // initialize node_num to 0
	}

	bool is_empty() { // time complexity O(1)
		return root == NULL; // PQ is empty
	}

	size_t size() { // time complexity O(1)
		return node_num;
	}

	pair_content get_top() { // time complexity O(1)
		if (root != NULL) {
			return { root->value , root->payload };
		}
		else {
			return { std::numeric_limits<float>::max(), std::numeric_limits<float>::max() }; // nonsense value is returned
		}
	}


	/* pair<ItemHandle, ValueType> record the node pointer and node key value (without offset) */
	ItemHandle insert(ValueType value, PayloadType payload) { // time complexity O(1)
		Node* new_node = new Node(); // new_node is a node pointer
		new_node->sibling = NULL;
		new_node->child = NULL;
		new_node->left_up = NULL;
		new_node->value = value;
		new_node->payload = payload;
		new_node->child_offset = 0; // initial offset is 0
		root = link(root, new_node);
		node_num++;
		return new_node; // return pointer to a new node
	}

	void add_to_heap(ValueType value) { // add a value to key values of all items; time complexity O(1)
		if (root != NULL) {
			root->value += value;
			root->child_offset += value;
		}
	}

	/* merge/insert target_heap into this heap; time complexity O(1)*/
	void meld(PairingHeapYS_with_offset& target_heap) {
		root = link(root, target_heap.root); // just merge two roots
		node_num += target_heap.size();
		target_heap.node_num = 0; // target_heap is marked empty
	}


	/* First delete links between root, left child and all the siblings of the left child.
	Then Merge tree subtrees that are obtained by detaching the left child and all siblings
	by the two pass method and delete the root node. Merge the detached subtrees from left to right
	in one pass and then merge the subtrees from right to left to form the new heap without
	violation of conditions of min-heap. This process takes O(log n) time where n is the number of
	nodes. https://www.geeksforgeeks.org/pairing-heap/
	this function will release momory of the popped out node */
	pair_content pop_top() { // time complexity O(log n)
		if (root == NULL) {
			return { std::numeric_limits<float>::max(), std::numeric_limits<float>::max() }; // return nonsense value
		}
		pair_content result = { root->value , root->payload }; // the popped out node
		std::vector<ItemHandle> buffer; // pointers to all children of root
		Node* cur_child = root->child; // the direct child
		Node* next_child = NULL;
		while (cur_child != NULL) { // put pointers to all children of root into buffer
			buffer.push_back(cur_child);
			next_child = cur_child->sibling;
			cur_child->value += root->child_offset; // cur_child->value is updated by offset value of root
			cur_child->child_offset += root->child_offset; // cur_child->child_offset is updated by offset value of root
			cur_child = next_child;
		}
		node_num--; // size minus 1
		two_pass_method(buffer);
		if (buffer.size() > 0) {
			root = buffer[buffer.size() - 1];// buffer[2]; this kind of removed the new root
			root->sibling = NULL; // root->sibling = NULL
			root->left_up = NULL; // root->left_up = NULL
			for (int ii = buffer.size() - 2; ii >= 0; --ii) {
				root = link(root, buffer[ii]); // first link buffer[1] and buffer[2], then link buffer[0] and the previous link result (root)
			}
		}
		else { // this heppens when the old root has no child
			root = NULL;
		}
		return result;
	}

	/*this function changes the key value of a node from from_value to to_value (these two values are true
	key values without offset);
	time complexity O(log n)*/
	void change_key(ItemHandle node, ValueType from_value, ValueType to_value) {
		//cout << "change_key " << node->payload << " from_value:" << from_value << " to_value:" << to_value << endl;
		ValueType additional_offset = from_value - node->value; // all offsets above node
		if (heap_type_compare(to_value, from_value)) { // to_value is more prioritized				
			if (node->left_up != NULL) { // node is not root, otherwise it's done
				std::vector<ItemHandle> buffer; // pointers to all children of node->left_up
				Node* cur_child = node->left_up->child; // the direct child
				node->left_up->child = NULL; // cut the link from node->left_up to all below nodes
				node->value += to_value - from_value;
				Node* next_child = NULL;
				while (cur_child != NULL) { // put pointers to all children of node->left_up into buffer
					buffer.push_back(cur_child);
					next_child = cur_child->sibling;
					cur_child->value += additional_offset;
					cur_child->child_offset += additional_offset; 
					cur_child = next_child;
				}				
				two_pass_method(buffer);
				for (int ii = buffer.size() - 1; ii >= 0; --ii) {
					root = link(root, buffer[ii]); // link all children of node->left_up to the root
				}
			}
			else {
				node->value += to_value - from_value;
			}
		}
		else { // to_value is less prioritized
			/*this is different from above in that all children of node should also be in buffer,
			and even when node is root, there are actions*/
			if (node->left_up != NULL) { // node is not root
				pair_content x = get_top();
				ValueType new_value_to_top = -std::numeric_limits<double>::max();
				change_key(node, from_value, new_value_to_top); // change node to top
				x = pop_top(); // pop node out
				node_num++; // pop_top() decreases node_num
				node->sibling = NULL;
				node->child = NULL;
				node->left_up = NULL;
				node->value = to_value;  // update node->value
				node->payload = x.second;
				node->child_offset = 0; // initial offset is 0
				root = link(root, node); // link the isolated node to the remaining part
			}
			else { // node is root
				pair_content x = pop_top(); // pop out node
				node_num++; // pop_top() decreases node_num
				node->sibling = NULL;
				node->child = NULL;
				node->left_up = NULL;
				node->value = to_value;  // update node->value
				node->payload = x.second;
				node->child_offset = 0; // initial offset is 0
				root = link(root, node); // link the isolated node to the remaining part
			}
		}
	}

	/* the input is root; the root is printed first */
	void print(ItemHandle node) {
		if (node != NULL) {
			std::cout << "<" << node->value << "," << node->payload << "> child_offset = " << node->child_offset << endl;
		}
		if (node->child != NULL) {
			print(node->child);
		}
		if (node->sibling != NULL) {
			print(node->sibling);
		}
	}

};



/* 
examples are as follows. 

-------------an example cpp main file-------------------------------

using namespace std;

#include <data_structures/PairingHeapYS_with_offset.h>

int main()
{
	example_PairingHeapYS_with_offset();
}

--------------------------------------------

*/


typedef PairingHeapYS_with_offset<double, int> PairingHeapType_example2; // define data types in heap

void example_PairingHeapYS_with_offset() {

	cout << "the following are the printing results of running an example of PairingHeapYS_with_offset" << endl;

	std::vector<PairingHeapType_example2::ItemHandle> pairing_heap_buffer; // define handles/pointers for items in heap
	PairingHeapType_example2 heap; // define heap with handles

	pairing_heap_buffer.push_back(heap.insert(0.1, 1));
	pairing_heap_buffer.push_back(heap.insert(0.2, 2));
	pairing_heap_buffer.push_back(heap.insert(0.3, 3));


	pair<double, int> old_top = heap.pop_top();
	pair<double, int> new_top = heap.get_top();
	cout << "old_top key: " << old_top.first << endl;
	cout << "new_top key: " << new_top.first << endl;

	cout << "heap.size(): " << heap.size() << endl;


	PairingHeapType_example2 heap2; // define heap with handles
	pairing_heap_buffer.push_back(heap2.insert(0.4, 4));
	pairing_heap_buffer.push_back(heap2.insert(0.5, 5));

	heap2.add_to_heap(-1);

	heap.meld(heap2);



	pair<double, int> meld_top = heap.get_top();
	cout << "meld_top key: " << meld_top.first << endl;
	cout << "heap.size(): " << heap.size() << endl;
	cout << "heap2.size(): " << heap2.size() << endl;

	heap.change_key(pairing_heap_buffer[1], 0.2, 1.1);

	heap.change_key(pairing_heap_buffer[3], -0.6, 1);


	while (heap.size() > 0) {
		pair<double, int> top = heap.pop_top();
		cout << "key: " << top.first << " item: " << top.second << endl;
	}

}

#include<chrono>
#include <math.h>
void test_PairingHeapYS_with_offset() {

	vector<double> sizes = { 1e3, 1e4, 1e5, 1e6, 1e7 };
	for (int ii = sizes.size() - 1; ii >= 0; ii--) {
		cout << "heap_YS n = " << sizes[ii] << endl;
		vector<double> keys(sizes[ii] * 2);
		auto beginx = std::chrono::high_resolution_clock::now();
		std::vector<PairingHeapType_example2::ItemHandle> pairing_heap_buffer(sizes[ii] * 2);
		PairingHeapType_example2 heap, heap2; // define heap with handles
		for (int j = sizes[ii] - 1; j >= 0; j--) {
			double x = 1 / (double)(rand() % (j + 1) + 1);
			pairing_heap_buffer[j] = heap.insert(x, j);
			//cout << "insert <" << x << "," << j << ">" << endl;
			keys[j] = x;
			x = 1 / (double)(rand() % (j + 1) + 1);
			pairing_heap_buffer[j + sizes[ii]] = heap2.insert(x, j + sizes[ii]);
			//cout << "insert <" << x << "," << j + sizes[ii] << ">" << endl;
			keys[j + sizes[ii]] = x;
		}
		heap.meld(heap2);
		//heap.print(heap.root);
		for (int j = sizes[ii] * 2 - 1; j >= 0; j--) {
			int ran = rand() % ((int)sizes[ii] * 2);
			double x = 1 / (double)(rand() % (j + 1) + 1);
			heap.change_key(pairing_heap_buffer[ran], keys[ran], x);
			//cout << "change_key <" << x << "," << ran << ">" << endl;
			keys[ran] = x;
		}
		//heap.print(heap.root);
		heap.add_to_heap(0.005);
		//heap.print(heap.root);
		for (int j = sizes[ii] * 2 - 1; j >= 0; j--) {
			keys[j] = keys[j] + 0.005;
		}
		int x = 0;
		pair<double, int> new_top = heap.get_top();
		double old_top_value = new_top.first;
		while (heap.size() > 0) {
			pair<double, int> top = heap.pop_top();
			//cout << "pop_top <" << top.first << "," << top.second << ">" << endl;
			if (abs(keys[top.second] - top.first) > 1e-4) {
				cout << "keys[top.second] != top.first error!" << endl;
				cout << "pop_top <" << top.first << "," << top.second << ">" << endl;
				cout << "keys[top.second] = " << keys[top.second] << endl;
				cout << "top.first = " << top.first << endl;
				exit(1);
			}
			if (old_top_value > top.first + 1e-5) {
				cout << "heap priority error!" << endl;
				cout << "pop_top <" << top.first << "," << top.second << ">" << endl;
				cout << "old_top_value = " << old_top_value << endl;
				exit(1);
			}
			old_top_value = top.first;
			x++;
		}
		if (x != sizes[ii] * 2) {
			cout << "heap size error! x = " << x << " sizes[ii] * 2 = "
				<< sizes[ii] * 2 << endl;
			exit(1);
		}
		auto endx = std::chrono::high_resolution_clock::now();
		double runningtimex = std::chrono::duration_cast<std::chrono::nanoseconds>(endx - beginx).count() / 1e9; // s
		cout << "time = " << runningtimex << "s; time / nlog n = " << runningtimex / x / log(x) << "s" << endl;
		/*release memory*/
		for (int j = pairing_heap_buffer.size() - 1; j >= 0; j--) {
			delete(pairing_heap_buffer[j]);
		}
	}
}

