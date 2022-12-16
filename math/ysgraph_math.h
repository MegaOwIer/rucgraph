#pragma once


double min_double(std::list<double> numbers) {

	/*if you use std::list<double>& above, then below codes cause errors*/

	double min = *numbers.begin();

	for (auto it = numbers.begin(); it != numbers.end(); it++) {
		if (min > *it) {
			min = *it;
		}
	}

	return min;

}




#pragma region
bool top_percent_threshold_compare(const double&i, const double&j)
{
	return i > j;  // < is from small to big; > is from big to small
}

double top_percent_threshold_change_vector(std::vector<double>& numbers, double top_p) {

	/*this function returns the threshold value such that values in numbers
	that are larger than or equal to this value are within the top top_p %
	of numbers from large to small;


	time complexity: O(log(numbers.size()) * numbers.size())

	you can do this function by maintaining a sorted top_vector of the top numbers.size() * top_p / 100 numbers;
	you check all the numbers, each time you find the place of a number possibly in top_vector
	by iteratively divide the search space in top_vector; in this way, if a number should be in top_vector,
	you can find its place in top_vector within logarithm time log_2(top_vector.size()).
	In this way, the time complexity of this function is O(log_2(top_vector.size()) * numbers.size()) */


	sort(numbers.begin(), numbers.end(), top_percent_threshold_compare); // from large to small

	return numbers[(int)numbers.size() * top_p / 100];

}
#pragma endregion top_percent_threshold