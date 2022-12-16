#pragma once


#include <text mining/parse_string.h> 

std::vector<string> parse_substring_between_pairs_of_delimiters
(std::string& parse_target, std::string delimiter1, std::string delimiter2) {

	/*this pair of delimiters must be different*/

	std::vector<string> results;

	std::vector<string> Parsed_content1 = parse_string(parse_target, delimiter1);
	for (int i = 0; i < Parsed_content1.size(); i++) {
		if (Parsed_content1[i].find(delimiter2) != std::string::npos)
		{
			std::vector<string> Parsed_content2 = parse_string(Parsed_content1[i], delimiter2);
			results.push_back(Parsed_content2[0]);
		}
	}

	return results;

}

