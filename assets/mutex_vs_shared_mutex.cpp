/*
mutex is slower than shared_mutex on both Linux and Windows.
For num = 1e8 in the following codes, it's 2.4s vs 1.8s on Windows,
and is 12e-8s vs 4.6e-8s on Linux.

-----------
The reason why ThreadPool using mutex is faster than ThreadPool2 using shared_mutex 
may be that condition_variable in ThreadPool is faster than condition_variable_any in ThreadPool2
*/

#include <iostream>
#include <mutex>
#include <shared_mutex>
using namespace std;

int main()
{
	int num = 1e8;
	int x = 1;
	{
		auto begin = std::chrono::high_resolution_clock::now();
		mutex l;
		for (int i = 0; i < num; i++) {
			l.lock();
			x++;
			l.unlock();
		}
		auto end = std::chrono::high_resolution_clock::now();
		double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
		cout << "mutex: " << runningtime << "s" << endl;
	}

	{
		auto begin = std::chrono::high_resolution_clock::now();
		shared_mutex l;
		for (int i = 0; i < num; i++) {
			l.lock();
			x++;
			l.unlock();
		}
		auto end = std::chrono::high_resolution_clock::now();
		double runningtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9; // s
		cout << "shared_mutex: " << runningtime << "s" << endl;
	}
	cout << x << endl;
}