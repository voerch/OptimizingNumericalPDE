#include "PDE.h"
#include "FDM.h"
#include "Option.h"
#include <iostream>
#include <chrono>
#include <fstream>

int main()
{
	//double x_dom = 1.0;     
	//long J = 20;
	//double t_dom = 0.075;
	//long N = 20;

	//HeatEqn* heat_pde = new HeatEqn;

	//ExplicitMethod euler(x_dom, J, t_dom, N, heat_pde);
	//euler.stepMarch();

	//CrankNicholson cn(x_dom, J, t_dom, N, heat_pde);
	//cn.stepMarch();

	//delete heat_pde;

	std::ofstream timing("BSTiming.csv");
	timing << "Explicit,CN" << std::endl;
	for (int timeCounter = 0; timeCounter < 1000; timeCounter++)
	{

		double K = 1.0;  // Strike price
		double r = 0.05;   // Risk-free rate (5%)
		double v = 0.2;    // Volatility of the underlying (20%)
		double T = 2.00;    // One year until expiry

		double x_dom = 2.0;       // Spot goes from [0.0, 2.0]
		unsigned long J = 50;
		double t_dom = T;         // Time period as for the option
		unsigned long N = 50;

		VanillaOption* callOption = new EurCall(K, r, T, v);
		BlackScholesPDE* bsPDE = new BlackScholesPDE(callOption);

		auto start = std::chrono::high_resolution_clock::now();
		ExplicitMethod euler(x_dom, J, t_dom, N, bsPDE);
		euler.stepMarch();
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;

		auto begin = std::chrono::high_resolution_clock::now();
		CrankNicholson cn(x_dom, J, t_dom, N, bsPDE);
		cn.stepMarch();
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed1 = end - begin;

		timing << elapsed.count() << "," << elapsed1.count() << std::endl;

		delete bsPDE;

	}

	timing.close();
	system("pause");
	return 0;

}