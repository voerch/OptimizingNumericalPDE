//#include "PDE.h"
//#include "FDM.h"
//#include "Option.h"
//#include <iostream>
//#include <chrono>
//#include <fstream>
//
//int main()
//{
//	long xStep = 50;
//	long tStep = 50;
//
//	///////////////////////////////////////
//
//	double x_dom = 1.0;     
//	double t_dom = 0.06;
//
//	HeatEqn* heat_pde = new HeatEqn;
//
//	//ExplicitMethod euler(x_dom, xStep, t_dom, tStep, heat_pde);
//	//euler.stepMarch();
//
//	auto CNStart = std::chrono::high_resolution_clock::now();
//	CrankNicholson cn(x_dom, xStep, t_dom, tStep, heat_pde);
//	cn.stepMarch();
//	auto CNFinish = std::chrono::high_resolution_clock::now();
//	std::chrono::duration<double> CNTime = CNFinish - CNStart;
//
//	std::cout << CNTime.count() << std::endl;
//
//	delete heat_pde;
//
//	////////////////////////////////////////
//
//	//double K = 1.0;  // Strike price
//	//double r = 0.05;   // Risk-free rate (5%)
//	//double v = 0.2;    // Volatility of the underlying (20%)
//	//double T = 2.00;    // One year until expiry
//
//	//double x_dom = 2.0;       // Spot goes from [0.0, 1.0]
//	//double t_dom = T;         // Time period as for the option
//
//
//	//VanillaOption* callOption = new EurCall(K, r, T, v);
//	//BlackScholesPDE* bsPDE = new BlackScholesPDE(callOption);
//
//	////LogSpotBlackScholesPDE* bsPDE = new LogSpotBlackScholesPDE(callOption);
//
//	//ExplicitMethod euler(x_dom, xStep, t_dom, tStep, bsPDE);
//	//euler.stepMarch();
//
//	//CrankNicholson cn(x_dom, xStep, t_dom, tStep, bsPDE);
//	//cn.stepMarch();
//
//	//delete bsPDE;
//	//
//
//	system("pause");
//}