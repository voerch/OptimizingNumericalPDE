//#include "PDE.h"
//#include "FDM.h"
//#include "Option.h"
//#include <iostream>
//#include <chrono>
//#include <fstream>
//#include <vector>
//
//using namespace std;
//
//int main()
//{
//	std::vector<double> xValues;
//	std::vector<double> tValues;

//	double x_dom = 1.0;     
//	double t_dom = 0.075;
//		
//	HeatEqn* heat_pde = new HeatEqn;
//
//	double AnalyticalSol;
//	std::ofstream optimalExplicit("heatOptimalExplicit.csv");
//	optimalExplicit << "GridSize,x,Error" << std::endl;
//
//	//std::ofstream optimalCN("heatOptimalCN.csv");
//	//optimalCN << "step,x,y,Error" << std::endl;
//
//	for (unsigned long j = 25; j < 101; j++)
//	{
//		xValues.resize(j, 0.0);
//		for (long xCounter = 0; xCounter < j; xCounter++)
//		{
//			xValues[xCounter] = static_cast<double>(xCounter) *  x_dom / static_cast<double>(j-1);
//		}
//
//		long J = j;
//		long N = j;
//		ExplicitMethod euler(x_dom, J, t_dom, N, heat_pde);
//		euler.stepMarch();
//
//		std::cout << heat_pde->HeatAnalyticalSolution(t_dom, xValues[5]) << std::endl << euler.oldResult[5];
//
//		for (int iCounter = 0; iCounter < static_cast<int>(j); iCounter++)
//		{
//			AnalyticalSol = heat_pde->HeatAnalyticalSolution(t_dom, xValues[iCounter]);
//			optimalExplicit << j << "," << xValues[iCounter] << "," << 100 * abs(euler.oldResult[iCounter] - AnalyticalSol) / AnalyticalSol << std::endl;
//		}
//		
//
//		
///*
//		CrankNicholson cn(x_dom, J, t_dom, N, heat_pde);
//		cn.stepMarch();
//		optimalCN << j << "," << j << "," << 100 * abs(cn.newResult[j - 1] - AnalyticalSol) / AnalyticalSol << std::endl;*/
//
//	}
//
//	//optimalCN.close();
//	optimalExplicit.close();
//	delete heat_pde;
//	
//	double K = 0.5;  // Strike price
//	double r = 0.05;   // Risk-free rate (5%)
//	double v = 0.2;    // Volatility of the underlying (20%)
//	double T = 1.00;    // One year until expiry
//	double x_dom = 1.0;       // Spot goes from [0.0, 1.0]
//	double t_dom = T;         // Time period as for the option
//	
//	VanillaOption* callOption = new EurCall(K, r, T, v);
//	BlackScholesPDE* bsPDE = new BlackScholesPDE(callOption);
//
//	double AnalyticalSol;
//	std::ofstream optimalExplicit("optimalExplicit.csv");
//	std::ofstream optimalCN("optimalCN.csv");
//
//	for (unsigned long xSize = 4; xSize <= 50; xSize++)
//	{
//		unsigned long J = xSize;
//		xValues.resize(xSize, 0.0);
//		for (long xCounter = 0; xCounter < xSize; xCounter++)
//		{
//			xValues[xCounter] = static_cast<double>(xCounter) *  x_dom / static_cast<double>(xSize - 1);
//		}
//
//		AnalyticalSol = callOption->PriceByBS(xValues[(xSize/2)]);
//
//		for (unsigned long tSize = 4; tSize <= 100; tSize++)
//		{
//			unsigned long N = tSize;
//			auto start = std::chrono::high_resolution_clock::now();
//			ExplicitMethod euler(x_dom, J, t_dom, N, bsPDE);
//			euler.stepMarch();
//			auto finish = std::chrono::high_resolution_clock::now();
//			chrono::duration<double> elapsed = finish - start;
//			optimalExplicit << xSize << "," << tSize << "," << abs(euler.newResult[(xSize/2)] - AnalyticalSol) / AnalyticalSol << "," << elapsed.count() << endl;
//
//			auto begin = std::chrono::high_resolution_clock::now();
//			CrankNicholson cn(x_dom, J, t_dom, N, bsPDE);
//			cn.stepMarch();
//			auto end = std::chrono::high_resolution_clock::now();
//			chrono::duration<double> sure = end - begin;
//			optimalCN << xSize << "," << tSize << "," << abs(cn.newResult[(xSize / 2)] - AnalyticalSol) / AnalyticalSol << "," << sure.count() << endl;
//		}
//	}
//	optimalCN.close();
//	optimalExplicit.close();
//	
//	delete bsPDE;
//
//	system("pause");
//	return 0;
//
//}