#pragma once

/*vector*/

std::vector<int> copy_vector_int(std::vector<int>& input_vector) {

	return input_vector;

}






/*list*/

std::list<int> copy_list_int(std::list<int>& input_list) {

	return input_list;

}





/*unordered_set*/

std::unordered_set<int> copy_unordered_set_int(std::unordered_set<int>& input_set) {

	return input_set;

}


std::vector<int> copy_unordered_set_int_to_vector_int(std::unordered_set<int>& input_set) {

	std::vector<int> output_vector;

	for (auto it = input_set.begin(); it != input_set.end(); it++) {
		output_vector.insert(output_vector.end(), *it);
	}

	return output_vector;

}