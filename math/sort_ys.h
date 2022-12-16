#pragma once

template <typename T1, typename T2> // T2 is priority value, int double etc
bool sort_pair_compare_large_to_small(const pair<T1, T2>& i, const pair<T1, T2>& j)
{
	return i.second > j.second;  // < is from small to big; > is from big to small.  sort by the second item of pair
}

template <typename T1, typename T2> // T2 is priority value, int double etc
bool sort_pair_compare_small_to_large(const pair<T1, T2>& i, const pair<T1, T2>& j)
{
	return i.second < j.second;  // < is from small to big; > is from big to small.  sort by the second item of pair
}

/*example:

sort(vector_pair.begin(), vector_pair.end(), sort_pair_compare_large_to_small);

why it doesn't work???

*/