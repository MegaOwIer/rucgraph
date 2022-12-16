#pragma once

class PLL_sorted_label {
public:
	int vertex, parent_vertex; // std::vector<PLL_sorted_label> should be sorted by vertex
	float distance;
}; 

void PLL_sorted_label_binary_insert(std::vector<PLL_sorted_label>& input_vector, int vertex, int parent_vertex, double distance) {

	/*insert <key, load> into vector, if key is already inside, then load is updated; time complexity O(log n + size()-position ), which is O(n) in the worst case, as
	the time complexity of inserting an element into a vector is the number of elements behind this element*/

	int left = 0, right = input_vector.size() - 1;

	while (left <= right) // it will be skept when input_vector.size() == 0
	{
		int mid = left + ((right - left) / 2); // mid is between left and right (may be equal); 
		if (input_vector[mid].vertex == vertex) { // vertex is already in input_vector, this cannot happen in PLL
			return;
		}
		else if (input_vector[mid].vertex > vertex) {
			right = mid - 1; // the elements after right are always either empty, or have larger keys than input key
		}
		else {
			left = mid + 1; // the elements before left are always either empty, or have smaller keys than input key
		}
	}

	/*the following code is used when key is not in vector, i.e., left > right, specifically, left = right + 1;
	the elements before left are always either empty, or have smaller keys than input key;
	the elements after right are always either empty, or have larger keys than input key;
	so, the input key should be insert between right and left at this moment*/
	PLL_sorted_label x;
	x.vertex = vertex;
	x.parent_vertex = parent_vertex;
	x.distance = distance;
	input_vector.insert(input_vector.begin() + left, x);

}



#include <text mining/binary_save_read_vector_of_vectors.h>
vector<vector<PLL_sorted_label>> PLL_read_indexes_vectorFORMAT_binary(string index_file_name) {

	vector<vector<PLL_sorted_label>> L;

	binary_read_vector_of_vectors(index_file_name, L);

	return L;

}

#include <text mining/parse_substring_between_pairs_of_delimiters.h>
vector<vector<PLL_sorted_label>> PLL_vector_of_vectors_read_indexes(string index_file_name, int max_N) {

	/*indexes; L(v_1) = { v_x, <distance to v_x, predecessor of v_x in the SP between v_1 and v_x>}*/
	vector<vector<PLL_sorted_label>> L(max_N);

	std::string line_content;
	std::ifstream myfile(index_file_name); // open the file
	if (myfile.is_open()) // if the file is opened successfully
	{
		long long int line_id = 0;

		while (getline(myfile, line_content)) // read file line by line
		{
			line_id++;
			if (line_id > 6) {
				std::vector<std::string> Parsed_content = parse_string(line_content, ":");
				if (Parsed_content.size() == 2) {
					int v = stoi(Parsed_content[0]);
					if (v > max_N - 1) {
						cout << "v=" << v << " max_N=" << max_N << endl;
					}
					std::vector<string> index_strings = parse_substring_between_pairs_of_delimiters(Parsed_content[1], "<", ">");
					for (auto it = index_strings.begin(); it != index_strings.end(); it++) {
						std::vector<string> an_index = parse_string(*it, ",");
						PLL_sorted_label_binary_insert(L[v], stoi(an_index[0]), stoi(an_index[2]), stod(an_index[1]));
						//cout << stoi(an_index[0]) << "," << stod(an_index[1]) << "," << stoi(an_index[2]) << endl;
						//getchar();
					}
				}
			}
		}
		myfile.close(); //close the file

		return L;
	}
	else
	{
		std::cout << "Unable to open file " << index_file_name << std::endl
			<< "Please check the file location or file name." << std::endl; // throw an error message
		getchar(); // keep the console window
		exit(1); // end the program
	}



}


double PLL_extract_distance_vectorFORMAT(vector<vector<PLL_sorted_label>>& L, int source, int terminal) {

	/*return std::numeric_limits<double>::max() is not connected*/

	double distance = std::numeric_limits<double>::max(); // if disconnected, retun this large value
	auto vector1_check_pointer = L[source].begin();
	auto vector2_check_pointer = L[terminal].begin();
	auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
	while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end) {
		if (vector1_check_pointer->vertex == vector2_check_pointer->vertex) {
			double dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
			if (distance > dis) {
				distance = dis;
			}
			vector1_check_pointer++;
		}
		else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex) {
			vector2_check_pointer++;
		}
		else {
			vector1_check_pointer++;
		}
	}

	return distance;

}

vector<pair<int, int>> PLL_extract_shortest_path_vectorFORMAT(vector<vector<PLL_sorted_label>>& L, int source, int terminal) {

	/*if source and terminal are disconnected, then return empty vector; note that, if source==terminal, then also return empty vector*/

	vector<pair<int, int>> path; // edges

	while (source != terminal) {

		int vector1_capped_v_parent = 0, vector2_capped_v_parent = 0;
		double distance = std::numeric_limits<double>::max(); // if disconnected, retun this large value
		bool connected = false;
		auto vector1_check_pointer = L[source].begin();
		auto vector2_check_pointer = L[terminal].begin();
		auto pointer_L_s_end = L[source].end(), pointer_L_t_end = L[terminal].end();
		while (vector1_check_pointer != pointer_L_s_end && vector2_check_pointer != pointer_L_t_end) {
			if (vector1_check_pointer->vertex == vector2_check_pointer->vertex) {
				connected = true;
				double dis = vector1_check_pointer->distance + vector2_check_pointer->distance;
				if (distance > dis) {
					distance = dis;
					vector1_capped_v_parent = vector1_check_pointer->parent_vertex;
					vector2_capped_v_parent = vector2_check_pointer->parent_vertex;
				}
				vector1_check_pointer++;
			}
			else if (vector1_check_pointer->vertex > vector2_check_pointer->vertex) {
				vector2_check_pointer++;
			}
			else {
				vector1_check_pointer++;
			}
		}

		if (connected) {
			/*the following code will not induce redundent edge, since for each vertex v_k, there is a label (v_k,0,v_k) in L[v_k];
			you must ascending from both source and terminal, otherwise you may not be able to extract SP */
			if (source != vector1_capped_v_parent) {
				path.push_back({ source , vector1_capped_v_parent });
				source = vector1_capped_v_parent; // ascending from source
			}
			if (terminal != vector2_capped_v_parent) {
				path.push_back({ terminal , vector2_capped_v_parent });
				terminal = vector2_capped_v_parent; // ascending from terminal
			}
		}
		else {
			break;
		}

	}

	return path;
}



