#include <Windows.h>        // For Sleep()
#include <iostream>         // For cout
#include <chrono>

extern "C" int TestASMFunc();
using namespace std;

int main()
{
	// Record start time
	auto start = std::chrono::high_resolution_clock::now();

	// Portion of code to be timed
	::Sleep(5000);

	// Record end time
	auto finish = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed = finish - start;

	std::cout << "Elapsed time: " << elapsed.count() << std::endl;
	cout << TestASMFunc() << endl;
	system("pause");
	return 0;
}