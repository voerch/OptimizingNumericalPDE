/*
Author: Mustafa Berke Erdis, August 2019

The main to time the runtime of the solutions under different conditions. 
Timings for Black-Scholes equation and heat equation are stored in "BSTiming.csv", "HeatTiming.csv" files.
1000 trials are outputted to the CSV files to ensure accurate timing.

*/

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

		double K = 1.0;			// Strike price
		double r = 0.05;		// Risk-free rate
		double v = 0.2;			// Volatility
		double T = 2.00;		// Maturity

		double x_dom = 2.0;		// Spot goes from [0.0, 2.0]
		unsigned long J = 64;	// x grid size
		double t_dom = T;
		unsigned long N = 64;	// t grid size

		VanillaOption* callOption = new EurCall(K, r, T, v);
		BlackScholesPDE* bsPDE = new BlackScholesPDE(callOption);

		// Timing explicit method for Black - Scholes equation.
		auto start = std::chrono::high_resolution_clock::now();
		ExplicitMethod euler(x_dom, J, t_dom, N, bsPDE);
		euler.timeMarch();
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;

		// Timing Crank - Nicolson method for Black - Scholes equation.
		auto begin = std::chrono::high_resolution_clock::now();
		CrankNicolson cn(x_dom, J, t_dom, N, bsPDE);
		cn.timeMarch();
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed1 = end - begin;

		timing << elapsed.count() << "," << elapsed1.count() << std::endl;

		delete bsPDE;

	}
	timing.close();

	std::ofstream heattiming("HeatTiming.csv");
	heattiming << "Explicit,CN,ADI" << std::endl;
	for (int timeCounter = 0; timeCounter < 1000; timeCounter++)
	{
		double x_dom = 1.0;     // Maximum value for x. [0, 1.0]
		long J = 64;			// x grid size.
		double t_dom = 0.06;	// Maximum value for t. [0, 0.06]
		long N = 64;			// t grid size.

		HeatEqn* heat_pde = new HeatEqn;

		// Timing explicit method for heat equation.
		auto ExpStart = std::chrono::high_resolution_clock::now();
		ExplicitMethod euler(x_dom, J, t_dom, N, heat_pde);
		euler.timeMarch();
		auto ExpFinish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> ExpTime = ExpFinish - ExpStart;

		// Timing Crank - Nicolson method for heat equation.
		auto CNStart = std::chrono::high_resolution_clock::now();
		CrankNicolson cn(x_dom, J, t_dom, N, heat_pde);
		cn.timeMarch();
		auto CNFinish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> CNTime = CNFinish - CNStart;

		delete heat_pde;

		long K = 32;			// x grid size.
		J = 32;					// y grid size.
		N = 32;					// t grid size.
		double y_dom = 1.0;		// Maximum value for y. [0, 1.0]
		Heat2D* heat2d = new Heat2D;

		// Timing alternating direction implicit method for two-dimensional heat equation.
		auto ADIStart = std::chrono::high_resolution_clock::now();
		ADI adi(x_dom, J, y_dom, K, t_dom, N, heat2d);
		adi.TimeMarch();
		auto ADIFinish = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double> ADITime = ADIFinish - ADIStart;

		delete heat2d;
		heattiming << ExpTime.count() << "," << CNTime.count() << "," << ADITime.count() << std::endl;
	}
	heattiming.close();

}