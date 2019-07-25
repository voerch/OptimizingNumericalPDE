#include "PDE2D.h"
#include "FDM.h"
#include <iostream>
#include <chrono>

int main()
{
	auto start = std::chrono::high_resolution_clock::now();

	double x_dom = 1.0;     
	long J = 10;
	double y_dom = 1.0;
	long K = 10;
	double t_dom = 0.075;
	long N = 10;

	HeatEqn* heat_pde = new HeatEqn;

	ADI adi(x_dom, J, y_dom, K, t_dom, N, heat_pde);
	adi.stepMarch();

	delete heat_pde;

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << " elapsed time: " << elapsed.count() << std::endl;
	
	system("pause");
	return 0;

}