#include "PDE.h"
#include "FDM.h"
#include <iostream>
#include <chrono>

int main()
{
	double x_dom = 2.0;     
	long J = 10;
	double y_dom = 1.0;
	long K = 10;
	double t_dom = 0.075;
	long N = 10;

	Heat2D* heat2d = new Heat2D;
	ADI adi(x_dom, J, y_dom, K, t_dom, N, heat2d);

	auto start = std::chrono::high_resolution_clock::now();
	adi.stepMarch();
	auto finish = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed = finish - start;
	std::cout << " elapsed time: " << elapsed.count() << std::endl;

	delete heat2d;

	system("pause");
	
}