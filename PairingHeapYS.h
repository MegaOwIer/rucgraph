#pragma once

#include<chrono>
/* this header file is for a pairing heap;
see how a pairing heap works: https://www.geeksforgeeks.org/pairing-heap/ 

A heap is merely a collection of subheaps associated with a value with a well-
   defined order, so there is no distinction between a heap and its nodes.

we assume that key values are smaller than std::numeric_limits<float>::max(), i.e., 3.40282e+38
(see the following get_top and pop_top functions),
otherwise there will be errors 

size of pointer: On 32-bit machine sizeof pointer is 32 bits ( 4 bytes),
	while on 64 bit machine it's 8 byte.
	most of the memory consumption of this heap is done by pointers */

template <typename ValueType, typename PayloadType> // input two data types
class PairingHeapYS {

private:
	struct Node {
		Node* sibling; // pointer to other node; 
		Node* child;
		Node* left_up; // the adj node from this node to the root, i.e., the adj_root
		ValueType value; // key value
		PayloadType payload; // item
	};
	typedef pair<ValueType, PayloadType> pair_content;
	size_t node_num; // size_t has the largest range
	string heap_type;

	bool heap_type_compare(ValueType i, ValueType j)
	{
		if (heap_type == "min") {
			return i < j;  // < is from small to big; > is from big to small
		}
		else {
			return i > j;
		}
	}

public: 
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
		//if (top_node->child == below_node) { // this should never happen
		//	std::cout << "PairingHeapYS top_node->child == below_node" << endl;
		//	getchar();
		//}
		below_node->sibling = top_node->child;
		below_node->left_up = top_node;
		top_node->child = below_node;
		top_node->sibling = NULL; // root->sibling = NULL
		top_node->left_up = NULL; // root->left_up = NULL
		return top_node; // return the new root pointer
	}

	/*the two pass method link a vector of nodes (usually in the same layer) by pairs; time complexity O(log n);
	
	https://en.wikipedia.org/wiki/Pairing_heap :
		The only non-trivial fundamental operation is the deletion of the minimum element from the heap.
		This requires performing repeated melds of its children until only one tree remains. The standard strategy
		first melds the subheaps in pairs (this is the step that gave this data structure its name) from left to right
		and then melds the resulting list of heaps from right to left.

		https://www.geeksforgeeks.org/pairing-heap/  : Merge tree subtrees that are
		obtained by detaching the left child and
		all siblings by the two pass method and delete the root node. (this is the step
		that gave this datastructure its name)

		this two pass method enables O(logn)
		https://stackoverflow.com/questions/22478773/why-does-pairing-heap-need-that-special-two-passes-when-delete-min :
		With this combining rule, the amortized complexity of delete-min is O(log n). 
		With the strict left-to-right rule, it's O(n).

		Example: buffer = {0 1 2 3 4 }
	*/
	void two_pass_method(std::vector<ItemHandle>& buffer) {
		size_t merged_children = 0;
		while (merged_children + 2 <= buffer.size()) {
			buffer[merged_children / 2] = link(
				buffer[merged_children], buffer[merged_children + 1]); // buffer[0]={0,1}, buffer[1]={2,3}, ..., buffer[4] = {4}
			merged_children += 2;
		}
		if (merged_children != buffer.size()) {
			buffer[merged_children / 2] = buffer[merged_children]; // buffer[2] = {4}
			buffer.resize(merged_children / 2 + 1); // buffer.resize(3); when buffer.size()=1, this resize buffer.size()=1
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
	A
	|
	B - C - D
	|       |
	E - F   G - H
	---------------------------
	An operation of puting B as child of A: 1. B->left_up = A; 2. A->child = B.

	*/

	PairingHeapYS(string input_heap_type) { // initialize a heap with heap type; time complexity O(1)
		if (input_heap_type != "min" && input_heap_type != "max") {
			std::cout << "Error: heap_type should be either min or max in PairingHeap" << std::endl;
			exit(1);
		}
		heap_type = input_heap_type;
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
			PayloadType u;
			return { std::numeric_limits<float>::max(), u }; // nonsense value is returned
			/*notably these values are smaller than std::numeric_limits<double>::max()
			in the following change_key and erase functions*/
		}
	}

	/* record all ItemHandles of inserted nodes; free memory using these handles after using this heap;
	 you must record the returned ItemHandle for memory free*/
	ItemHandle insert(ValueType value, PayloadType payload) { // time complexity O(1)
		Node* new_node = new Node(); // new_node is a node pointer
		new_node->sibling = NULL;
		new_node->child = NULL;
		new_node->left_up = NULL;
		new_node->value = value;
		new_node->payload = payload;
		root = link(root, new_node);
		node_num++;
		return new_node; // return pointer to a new node
	}

	/* merge/insert target_heap into this heap; time complexity O(1);
	 this function just point the content of target_heap into this heap, and it
	 does not change memory consumption */
	void meld(PairingHeapYS& target_heap) {
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
			PayloadType u;
			return { std::numeric_limits<float>::max(), u }; // return nonsense value
		}
		pair_content result = { root->value , root->payload };
		std::vector<ItemHandle> buffer; // pointers to all children of root
		Node* cur_child = root->child;
		Node* next_child = NULL;
		while (cur_child != NULL) { // put pointers to all children of root into buffer
			buffer.push_back(cur_child);
			next_child = cur_child->sibling;
			cur_child = next_child;
		}
		node_num--; // size minus 1
		two_pass_method(buffer);
		if (buffer.size() > 0) {
			root = buffer[buffer.size() - 1];// buffer[2]; this kind of removed the new root
			root->sibling = NULL;
			root->left_up = NULL;
			for (int ii = buffer.size() - 2; ii >= 0; --ii) {
				root = link(root, buffer[ii]); // first link buffer[1] and buffer[2], then link buffer[0] and the previous link result (root)
			}
		}
		else { // this heppens when the old root has no child
			root = NULL;
		}
		return result;
	}


	/*this function changes the key value of a node; time complexity O(log n)*/
	void change_key(ItemHandle node, ValueType to_value) {
		if (heap_type_compare(to_value, node->value)) { // to_value is more prioritized		
			node->value = to_value; // update node->value
			if (node->left_up != NULL) { // node is not root, otherwise it's done			
				std::vector<ItemHandle> buffer; // pointers to all children of node->left_up
				Node* cur_child = node->left_up->child; // the direct child
				node->left_up->child = NULL; // cut the link from node->left_up to all below nodes
				Node* next_child = NULL;
				while (cur_child != NULL) { // put pointers to all children of node->left_up into buffer
					buffer.push_back(cur_child);
					next_child = cur_child->sibling;
					cur_child = next_child;
				}				
				two_pass_method(buffer);
				for (int ii = buffer.size() - 1; ii >= 0; ii--) {
					root = link(root, buffer[ii]); // link all children of node->left_up to the root
				}
			}
		}
		else { // to_value is less prioritized
			/*this is different from above in that all children of node should also be in buffer,
			and even when node is root, there are actions*/
			if (node->left_up != NULL) { // node is not root
				pair_content x = get_top();
				ValueType new_value_to_top;
				if (heap_type == "min") {
					new_value_to_top = -std::numeric_limits<double>::max();
				}
				else {
					new_value_to_top = std::numeric_limits<double>::max();
				}
				change_key(node, new_value_to_top); // change node to top
				x = pop_top(); // pop node out
				node_num++; // pop_top() decreases node_num
				node->sibling = NULL;
				node->child = NULL;
				node->left_up = NULL;
				node->value = to_value;  // update node->value
				node->payload = x.second;
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
				root = link(root, node); // link the isolated node to the remaining part
			}
		}
	}

	/*this function erase an item; time complexity O(log n)*/
	void erase(ItemHandle node) {
		pair_content x = get_top();
		ValueType new_value_to_top;
		if (heap_type == "min") {
			new_value_to_top = -std::numeric_limits<double>::max();
		}
		else {
			new_value_to_top = std::numeric_limits<double>::max();
		}
		change_key(node, new_value_to_top); // change node to top 
		x = pop_top(); // pop node out
	}

	/* the input is root; the root is printed first */
	void print(ItemHandle node) {
		if (node != NULL) {
			std::cout << "<" << node->value << "," << node->payload << ">" << endl;
		}
		if (node->child != NULL) {
			print(node->child);
		}
		if (node->sibling != NULL) {
			print(node->sibling);
		}
	}

};



template <typename ValueType, typename PayloadType> // input two data types
class PairingHeapYS_min {

private:
	struct Node {
		Node* sibling; // pointer to other node; 
		Node* child;
		Node* left_up; // the adj node from this node to the root, i.e., the adj_root
		ValueType value; // key value
		PayloadType payload; // item
	};
	typedef pair<ValueType, PayloadType> pair_content;
	size_t node_num; // size_t has the largest range
	string heap_type;

	bool heap_type_compare(ValueType i, ValueType j)
	{
		return i < j;  // < is from small to big; > is from big to small
	}

public:
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
		//if (top_node->child == below_node) { // this should never happen
		//	std::cout << "PairingHeapYS top_node->child == below_node" << endl;
		//	getchar();
		//}
		below_node->sibling = top_node->child;
		below_node->left_up = top_node;
		top_node->child = below_node;
		top_node->sibling = NULL; // root->sibling = NULL
		top_node->left_up = NULL; // root->left_up = NULL
		return top_node; // return the new root pointer
	}

	/*the two pass method link a vector of nodes (usually in the same layer) by pairs; time complexity O(log n);

	https://en.wikipedia.org/wiki/Pairing_heap :
		The only non-trivial fundamental operation is the deletion of the minimum element from the heap.
		This requires performing repeated melds of its children until only one tree remains. The standard strategy
		first melds the subheaps in pairs (this is the step that gave this data structure its name) from left to right
		and then melds the resulting list of heaps from right to left.

		https://www.geeksforgeeks.org/pairing-heap/  : Merge tree subtrees that are
		obtained by detaching the left child and
		all siblings by the two pass method and delete the root node. (this is the step
		that gave this datastructure its name)

		this two pass method enables O(logn)
		https://stackoverflow.com/questions/22478773/why-does-pairing-heap-need-that-special-two-passes-when-delete-min :
		With this combining rule, the amortized complexity of delete-min is O(log n).
		With the strict left-to-right rule, it's O(n).

		Example: buffer = {0 1 2 3 4 }
	*/
	void two_pass_method(std::vector<ItemHandle>& buffer) {
		size_t merged_children = 0;
		while (merged_children + 2 <= buffer.size()) {
			buffer[merged_children / 2] = link(
				buffer[merged_children], buffer[merged_children + 1]); // buffer[0]={0,1}, buffer[1]={2,3}, ..., buffer[4] = {4}
			merged_children += 2;
		}
		if (merged_children != buffer.size()) {
			buffer[merged_children / 2] = buffer[merged_children]; // buffer[2] = {4}
			buffer.resize(merged_children / 2 + 1); // buffer.resize(3); when buffer.size()=1, this resize buffer.size()=1
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
	A
	|
	B - C - D
	|       |
	E - F   G - H
	---------------------------
	An operation of puting B as child of A: 1. B->left_up = A; 2. A->child = B.

	*/

	PairingHeapYS_min() { // initialize a heap with heap type; time complexity O(1)
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
			PayloadType u;
			return { std::numeric_limits<float>::max(), u }; // nonsense value is returned
			/*notably these values are smaller than std::numeric_limits<double>::max()
			in the following change_key and erase functions*/
		}
	}

	/* record all ItemHandles of inserted nodes; free memory using these handles after using this heap;
	 you must record the returned ItemHandle for memory free*/
	ItemHandle insert(ValueType value, PayloadType payload) { // time complexity O(1)
		Node* new_node = new Node(); // new_node is a node pointer
		new_node->sibling = NULL;
		new_node->child = NULL;
		new_node->left_up = NULL;
		new_node->value = value;
		new_node->payload = payload;
		root = link(root, new_node);
		node_num++;
		return new_node; // return pointer to a new node
	}

	/* merge/insert target_heap into this heap; time complexity O(1);
	 this function just point the content of target_heap into this heap, and it
	 does not change memory consumption */
	void meld(PairingHeapYS_min& target_heap) {
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
			PayloadType u;
			return { std::numeric_limits<float>::max(), u }; // return nonsense value
		}
		pair_content result = { root->value , root->payload };
		std::vector<ItemHandle> buffer; // pointers to all children of root
		Node* cur_child = root->child;
		Node* next_child = NULL;
		while (cur_child != NULL) { // put pointers to all children of root into buffer
			buffer.push_back(cur_child);
			next_child = cur_child->sibling;
			cur_child = next_child;
		}
		node_num--; // size minus 1
		two_pass_method(buffer);
		if (buffer.size() > 0) {
			root = buffer[buffer.size() - 1];// buffer[2]; this kind of removed the new root
			root->sibling = NULL;
			root->left_up = NULL;
			for (int ii = buffer.size() - 2; ii >= 0; --ii) {
				root = link(root, buffer[ii]); // first link buffer[1] and buffer[2], then link buffer[0] and the previous link result (root)
			}
		}
		else { // this heppens when the old root has no child
			root = NULL;
		}
		return result;
	}


	/*this function changes the key value of a node; time complexity O(log n)*/
	void change_key(ItemHandle node, ValueType to_value) {
		if (heap_type_compare(to_value, node->value)) { // to_value is more prioritized		
			node->value = to_value; // update node->value
			if (node->left_up != NULL) { // node is not root, otherwise it's done			
				std::vector<ItemHandle> buffer; // pointers to all children of node->left_up
				Node* cur_child = node->left_up->child; // the direct child
				node->left_up->child = NULL; // cut the link from node->left_up to all below nodes
				Node* next_child = NULL;
				while (cur_child != NULL) { // put pointers to all children of node->left_up into buffer
					buffer.push_back(cur_child);
					next_child = cur_child->sibling;
					cur_child = next_child;
				}
				two_pass_method(buffer);
				for (int ii = buffer.size() - 1; ii >= 0; ii--) {
					root = link(root, buffer[ii]); // link all children of node->left_up to the root
				}
			}
		}
		else { // to_value is less prioritized
			/*this is different from above in that all children of node should also be in buffer,
			and even when node is root, there are actions*/
			if (node->left_up != NULL) { // node is not root
				pair_content x = get_top();
				ValueType new_value_to_top;
				new_value_to_top = -std::numeric_limits<double>::max();
				change_key(node, new_value_to_top); // change node to top
				x = pop_top(); // pop node out
				node_num++; // pop_top() decreases node_num
				node->sibling = NULL;
				node->child = NULL;
				node->left_up = NULL;
				node->value = to_value;  // update node->value
				node->payload = x.second;
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
				root = link(root, node); // link the isolated node to the remaining part
			}
		}
	}

	/*this function erase an item; time complexity O(log n)*/
	void erase(ItemHandle node) {
		pair_content x = get_top();
		ValueType new_value_to_top;
		new_value_to_top = -std::numeric_limits<double>::max();
		change_key(node, new_value_to_top); // change node to top 
		x = pop_top(); // pop node out
	}

	/* the input is root; the root is printed first */
	void print(ItemHandle node) {
		if (node != NULL) {
			std::cout << "<" << node->value << "," << node->payload << ">" << endl;
		}
		if (node->child != NULL) {
			print(node->child);
		}
		if (node->sibling != NULL) {
			print(node->sibling);
		}
	}

};




template <typename ValueType, typename PayloadType> // input two data types
class PairingHeapYS_max {

private:
	struct Node {
		Node* sibling; // pointer to other node; 
		Node* child;
		Node* left_up; // the adj node from this node to the root, i.e., the adj_root
		ValueType value; // key value
		PayloadType payload; // item
	};
	typedef pair<ValueType, PayloadType> pair_content;
	size_t node_num; // size_t has the largest range
	string heap_type;

	bool heap_type_compare(ValueType i, ValueType j)
	{
		return i > j;  // < is from small to big; > is from big to small
	}

public:
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
		//if (top_node->child == below_node) { // this should never happen
		//	std::cout << "PairingHeapYS top_node->child == below_node" << endl;
		//	getchar();
		//}
		below_node->sibling = top_node->child;
		below_node->left_up = top_node;
		top_node->child = below_node;
		top_node->sibling = NULL; // root->sibling = NULL
		top_node->left_up = NULL; // root->left_up = NULL
		return top_node; // return the new root pointer
	}

	/*the two pass method link a vector of nodes (usually in the same layer) by pairs; time complexity O(log n);

	https://en.wikipedia.org/wiki/Pairing_heap :
		The only non-trivial fundamental operation is the deletion of the minimum element from the heap.
		This requires performing repeated melds of its children until only one tree remains. The standard strategy
		first melds the subheaps in pairs (this is the step that gave this data structure its name) from left to right
		and then melds the resulting list of heaps from right to left.

		https://www.geeksforgeeks.org/pairing-heap/  : Merge tree subtrees that are
		obtained by detaching the left child and
		all siblings by the two pass method and delete the root node. (this is the step
		that gave this datastructure its name)

		this two pass method enables O(logn)
		https://stackoverflow.com/questions/22478773/why-does-pairing-heap-need-that-special-two-passes-when-delete-min :
		With this combining rule, the amortized complexity of delete-min is O(log n).
		With the strict left-to-right rule, it's O(n).

		Example: buffer = {0 1 2 3 4 }
	*/
	void two_pass_method(std::vector<ItemHandle>& buffer) {
		size_t merged_children = 0;
		while (merged_children + 2 <= buffer.size()) {
			buffer[merged_children / 2] = link(
				buffer[merged_children], buffer[merged_children + 1]); // buffer[0]={0,1}, buffer[1]={2,3}, ..., buffer[4] = {4}
			merged_children += 2;
		}
		if (merged_children != buffer.size()) {
			buffer[merged_children / 2] = buffer[merged_children]; // buffer[2] = {4}
			buffer.resize(merged_children / 2 + 1); // buffer.resize(3); when buffer.size()=1, this resize buffer.size()=1
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
	A
	|
	B - C - D
	|       |
	E - F   G - H
	---------------------------
	An operation of puting B as child of A: 1. B->left_up = A; 2. A->child = B.

	*/

	PairingHeapYS_max() { // initialize a heap with heap type; time complexity O(1)
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
			PayloadType u;
			return { std::numeric_limits<float>::max(), u }; // nonsense value is returned
			/*notably these values are smaller than std::numeric_limits<double>::max()
			in the following change_key and erase functions*/
		}
	}

	/* record all ItemHandles of inserted nodes; free memory using these handles after using this heap;
	 you must record the returned ItemHandle for memory free*/
	ItemHandle insert(ValueType value, PayloadType payload) { // time complexity O(1)
		Node* new_node = new Node(); // new_node is a node pointer
		new_node->sibling = NULL;
		new_node->child = NULL;
		new_node->left_up = NULL;
		new_node->value = value;
		new_node->payload = payload;
		root = link(root, new_node);
		node_num++;
		return new_node; // return pointer to a new node
	}

	/* merge/insert target_heap into this heap; time complexity O(1);
	 this function just point the content of target_heap into this heap, and it
	 does not change memory consumption */
	void meld(PairingHeapYS_max& target_heap) {
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
			PayloadType u;
			return { std::numeric_limits<float>::max(), u }; // return nonsense value
		}
		pair_content result = { root->value , root->payload };
		std::vector<ItemHandle> buffer; // pointers to all children of root
		Node* cur_child = root->child;
		Node* next_child = NULL;
		while (cur_child != NULL) { // put pointers to all children of root into buffer
			buffer.push_back(cur_child);
			next_child = cur_child->sibling;
			cur_child = next_child;
		}
		node_num--; // size minus 1
		two_pass_method(buffer);
		if (buffer.size() > 0) {
			root = buffer[buffer.size() - 1];// buffer[2]; this kind of removed the new root
			root->sibling = NULL;
			root->left_up = NULL;
			for (int ii = buffer.size() - 2; ii >= 0; --ii) {
				root = link(root, buffer[ii]); // first link buffer[1] and buffer[2], then link buffer[0] and the previous link result (root)
			}
		}
		else { // this heppens when the old root has no child
			root = NULL;
		}
		return result;
	}


	/*this function changes the key value of a node; time complexity O(log n)*/
	void change_key(ItemHandle node, ValueType to_value) {
		if (heap_type_compare(to_value, node->value)) { // to_value is more prioritized		
			node->value = to_value; // update node->value
			if (node->left_up != NULL) { // node is not root, otherwise it's done			
				std::vector<ItemHandle> buffer; // pointers to all children of node->left_up
				Node* cur_child = node->left_up->child; // the direct child
				node->left_up->child = NULL; // cut the link from node->left_up to all below nodes
				Node* next_child = NULL;
				while (cur_child != NULL) { // put pointers to all children of node->left_up into buffer
					buffer.push_back(cur_child);
					next_child = cur_child->sibling;
					cur_child = next_child;
				}
				two_pass_method(buffer);
				for (int ii = buffer.size() - 1; ii >= 0; ii--) {
					root = link(root, buffer[ii]); // link all children of node->left_up to the root
				}
			}
		}
		else { // to_value is less prioritized
			/*this is different from above in that all children of node should also be in buffer,
			and even when node is root, there are actions*/
			if (node->left_up != NULL) { // node is not root
				pair_content x = get_top();
				ValueType new_value_to_top;
				new_value_to_top = std::numeric_limits<double>::max();
				change_key(node, new_value_to_top); // change node to top
				x = pop_top(); // pop node out
				node_num++; // pop_top() decreases node_num
				node->sibling = NULL;
				node->child = NULL;
				node->left_up = NULL;
				node->value = to_value;  // update node->value
				node->payload = x.second;
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
				root = link(root, node); // link the isolated node to the remaining part
			}
		}
	}

	/*this function erase an item; time complexity O(log n)*/
	void erase(ItemHandle node) {
		pair_content x = get_top();
		ValueType new_value_to_top;
		new_value_to_top = std::numeric_limits<double>::max();
		change_key(node, new_value_to_top); // change node to top 
		x = pop_top(); // pop node out
	}

	/* the input is root; the root is printed first */
	void print(ItemHandle node) {
		if (node != NULL) {
			std::cout << "<" << node->value << "," << node->payload << ">" << endl;
		}
		if (node->child != NULL) {
			print(node->child);
		}
		if (node->sibling != NULL) {
			print(node->sibling);
		}
	}

};



/* examples: */

typedef PairingHeapYS<double, int> PairingHeapType_example; // define data types in heap

void example_PairingHeapYS() {

	std::vector<PairingHeapType_example::ItemHandle> pairing_heap_buffer; // define handles/pointers for items in heap
	PairingHeapType_example heap("min"); // define heap with handles

	pairing_heap_buffer.push_back(heap.insert(1.1, 1));
	pairing_heap_buffer.push_back(heap.insert(1.5, 2));
	pairing_heap_buffer.push_back(heap.insert(0.5, 3));


	pair<double, int> old_top = heap.pop_top();
	pair<double, int> new_top = heap.get_top();
	cout << "old_top key: " << old_top.first << endl;
	cout << "new_top key: " << new_top.first << endl;

	cout << "heap.size(): " << heap.size() << endl;


	PairingHeapType_example heap2("min"); // define heap with handles
	pairing_heap_buffer.push_back(heap2.insert(0.1, 4));
	pairing_heap_buffer.push_back(heap2.insert(0.2, 5));

	heap.meld(heap2);

	pair<double, int> meld_top = heap.get_top();
	cout << "meld_top key: " << meld_top.first << endl;
	cout << "heap.size(): " << heap.size() << endl;
	cout << "heap2.size(): " << heap2.size() << endl;

	heap.change_key(pairing_heap_buffer[0], 1.1);
	heap.change_key(pairing_heap_buffer[0], 1.2);

	while (heap.size() > 0) {
		pair<double, int> top = heap.pop_top();
		cout << "key: " << top.first << " item: " << top.second << endl;
	}

	if (1 > 0) {
		int sizex = 1e5;
		for (int i = 0; i < sizex; i++) {
			pairing_heap_buffer.push_back(heap.insert(1.1, 1));
		}
	}
}


#include<cmath>
void test_PairingHeapYS() {
	vector<string> types = { "min", "max" };
	for (int mm = 0; mm < 2; mm++) {
		string type = types[mm];
		vector<double> sizes = { 1e3, 1e4, 1e5 };
		for (int ii = sizes.size() - 1; ii >= 0; ii--) {
			cout << type << " heap_YS n = " << sizes[ii] << endl;
			vector<double> keys(sizes[ii] * 2);
			auto beginx = std::chrono::high_resolution_clock::now();
			std::vector<PairingHeapType_example::ItemHandle> pairing_heap_buffer(sizes[ii] * 2);
			PairingHeapType_example heap(type), heap2(type); // define heap with handles
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
			for (int j = sizes[ii] * 2 - 1; j >= 0; j--) {
				int ran = rand() % ((int)sizes[ii] * 2);
				double x = 1 / (double)(rand() % (j + 1) + 1);
				heap.change_key(pairing_heap_buffer[ran], x);
				//cout << "change_key <" << x << "," << ran << ">" << endl;
				keys[ran] = x;
			}
			int x = 0;
			pair<double, int> new_top = heap.get_top();
			double value = new_top.first;
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
				if (type == "min") {
					if (top.first < value) {
						cout << "heap priority error!" << endl;
						exit(1);
					}
				}
				else {
					if (top.first > value) {
						cout << "heap priority error!" << endl;
						exit(1);
					}
				}
				value = top.first;
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
}


