#include "PDE.h"
#include "FDM.h"
#include "Option.h"
#include <iostream>
#include <chrono>
#include <fstream>

int main()
{
	std::ofstream timing("BSTiming.csv");
	timing << "Explicit,CN" << std::endl;
	for (int timeCounter = 0; timeCounter < 1000; timeCounter++)
	{

		double K = 1.0;  // Strike price
		double r = 0.05;   // Risk-free rate (5%)
		double v = 0.2;    // Volatility of the underlying (20%)
		double T = 2.00;    // One year until expiry

		double x_dom = 2.0;       // Spot goes from [0.0, 2.0]
		unsigned long J = 64;
		double t_dom = T;         // Time period as for the option
		unsigned long N = 64;

		VanillaOption* callOption = new EurCall(K, r, T, v);
		BlackScholesPDE* bsPDE = new BlackScholesPDE(callOption);

		auto start = std::chrono::high_resolution_clock::now();
		ExplicitMethod euler(x_dom, J, t_dom, N, bsPDE);
		euler.timeMarch();
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;

		auto begin = std::chrono::high_resolution_clock::now();
		CrankNicholson cn(x_dom, J, t_dom, N, bsPDE);
		cn.timeMarch();
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed1 = end - begin;

		timing << elapsed.count() << "," << elapsed1.count() << std::endl;

		delete bsPDE;

	}
	timing.close();

	/*std::ofstream heattiming("HeatTiming.csv");
	heattiming << "Explicit,CN,ADI" << std::endl;
	for (int timeCounter = 0; timeCounter < 1000; timeCounter++)
	{
		double x_dom = 1.0;     
		long J = 50;
		double t_dom = 0.06;
		long N = 50;

		HeatEqn* heat_pde = new HeatEqn;

		auto ExpStart = std::chrono::high_resolution_clock::now();
		ExplicitMethod euler(x_dom, J, t_dom, N, heat_pde);
		euler.timeMarch();
		auto ExpFinish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> ExpTime = ExpFinish - ExpStart;

		auto CNStart = std::chrono::high_resolution_clock::now();
		CrankNicholson cn(x_dom, J, t_dom, N, heat_pde);
		cn.timeMarch();
		auto CNFinish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> CNTime = CNFinish - CNStart;

		delete heat_pde;

		long K = 20;
		J = 20;
		N = 20;
		double y_dom = 1.0;

		Heat2D* heat2d = new Heat2D;

		auto ADIStart = std::chrono::high_resolution_clock::now();
		ADI adi(x_dom, J, y_dom, K, t_dom, N, heat2d);
		adi.TimeMarch();
		auto ADIFinish = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> ADITime = ADIFinish - ADIStart;


		delete heat2d;

		heattiming << ExpTime.count() << "," << CNTime.count() << "," << ADITime.count() << std::endl;

	}
	heattiming.close();*/

}