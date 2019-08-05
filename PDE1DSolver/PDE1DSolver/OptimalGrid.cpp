#include "PDE.h"
#include "FDM.h"
#include "Option.h"
#include <iostream>
#include <chrono>
#include <fstream>
#include <vector>

using namespace std;

int main()
{
	std::vector<double> xValues;
	unsigned long N, J;
	double AnalyticalSol;
	double Error;
	unsigned long gridMax = 50;

	/////////////////////////////////////////////////////////////

	/*HeatEqn* heatpde = new HeatEqn;

	double x_dom = 1.0;
	double t_dom = 0.05;

	std::ofstream optimalExplicit("HeatOptimalExplicit.csv");
	optimalExplicit << "xSize,tSize,Solution,Analytic,Error,Time" << endl;

	for (unsigned long xSize = 10; xSize <= gridMax; xSize++)
	{
		J = xSize;
		xValues.resize(xSize, 0.0);
		for (unsigned long xCounter = 0; xCounter < xSize; xCounter++)
		{
			xValues[xCounter] = static_cast<double>(xCounter) *  x_dom / static_cast<double>(xSize - 1);
		}

		AnalyticalSol = heatpde->AnalyticSol(t_dom, xValues[(xSize / 2)]);
			
		for (unsigned long tSize = 10; tSize <= gridMax; tSize++)
		{
			N = tSize;

			auto start = std::chrono::high_resolution_clock::now();
			ExplicitMethod euler(x_dom, J, t_dom, N, heatpde);
			euler.stepMarch();
			auto finish = std::chrono::high_resolution_clock::now();
			chrono::duration<double> elapsed = finish - start;
			optimalExplicit << xSize << "," << tSize << "," << euler.newResult[(xSize / 2)] << "," << AnalyticalSol << "," << abs(euler.newResult[(xSize / 2)] - AnalyticalSol) / AnalyticalSol << "," << elapsed.count() << endl;
		}
	}
	optimalExplicit.close();

	std::ofstream optimalCN("HeatOptimalCN.csv");
	optimalCN << "xSize,tSize,Solution,Analytic,Error,Time" << endl;
	for (unsigned long xSize = 10; xSize <= gridMax; xSize++)
	{
		J = xSize;
		xValues.resize(xSize, 0.0);
		for (unsigned long xCounter = 0; xCounter < xSize; xCounter++)
		{
			xValues[xCounter] = static_cast<double>(xCounter) *  x_dom / static_cast<double>(xSize - 1);
		}

		AnalyticalSol = heatpde->AnalyticSol(t_dom, xValues[(xSize / 2)]);
		for (unsigned long tSize = 10; tSize <= gridMax; tSize++)
		{
			N = tSize;

			auto begin = std::chrono::high_resolution_clock::now();
			CrankNicholson cn(x_dom, J, t_dom, N, heatpde);
			cn.stepMarch();
			auto end = std::chrono::high_resolution_clock::now();
			chrono::duration<double> sure = end - begin;
			Error = abs(cn.newResult[(xSize / 2)] - AnalyticalSol) / AnalyticalSol;
			if (Error > 1.0)
				Error = 0.0;
			optimalCN << xSize << "," << tSize << "," << cn.newResult[(xSize / 2)] << "," << AnalyticalSol << "," << Error << "," << sure.count() << endl;
		}
	}
	optimalCN.close();


	delete heatpde;*/
	
	///////////////////////////////////////////////////////////////////

	double K = 1;  // Strike price
	double r = 0.05;   // Risk-free rate (5%)
	double v = 0.2;    // Volatility of the underlying (20%)
	double T = 2.00;    // 2 year until expiry
	double x_dom = 2.00;       // Spot goes from [0.0, 1.0]
	double t_dom = T;         // Time period as for the option

	VanillaOption* callOption = new EurCall(K, r, T, v);
	BlackScholesPDE* bsPDE = new BlackScholesPDE(callOption);
	//LogSpotBlackScholesPDE* bsPDE = new LogSpotBlackScholesPDE(callOption);

	std::ofstream optimalExplicit("BSoptimalExplicit.csv");

	optimalExplicit << "xSize,tSize,Solution,Analytic,Error,Time" << endl;

	for (unsigned long xSize = 10; xSize <= gridMax ; xSize++)
	{
		J = xSize;
		xValues.resize(xSize, 0.0);
		for (unsigned long xCounter = 0; xCounter < xSize; xCounter++)
		{
			xValues[xCounter] = static_cast<double>(xCounter) *  x_dom / static_cast<double>(xSize - 1);
		}

		AnalyticalSol = callOption->PriceByBS(xValues[(xSize/2)]);
		for (unsigned long tSize = 10; tSize <= gridMax; tSize++)
		{
			 N = tSize;

			auto start = std::chrono::high_resolution_clock::now();
			ExplicitMethod euler(x_dom, J, t_dom, N, bsPDE);
			euler.stepMarch();
			auto finish = std::chrono::high_resolution_clock::now();
			chrono::duration<double> elapsed = finish - start;
			optimalExplicit << xSize << "," << tSize << "," << euler.newResult[(xSize / 2)] << "," << AnalyticalSol << "," << abs(euler.newResult[(xSize/2)] - AnalyticalSol) / AnalyticalSol << "," << elapsed.count() << endl;
		}
	}
	optimalExplicit.close();

	std::ofstream optimalCN("BSoptimalCN.csv");
	optimalCN << "xSize,tSize,Solution,Analytic,Error,Time" << endl;
	for (unsigned long xSize = 10; xSize <= gridMax; xSize++)
	{
		J = xSize;
		xValues.resize(xSize, 0.0);
		for (unsigned long xCounter = 0; xCounter < xSize; xCounter++)
		{
			xValues[xCounter] = static_cast<double>(xCounter) *  x_dom / static_cast<double>(xSize - 1);
		}

		AnalyticalSol = callOption->PriceByBS(xValues[(xSize / 2)]);
		for (unsigned long tSize = 10; tSize <= gridMax; tSize++)
		{
			N = tSize;

			auto begin = std::chrono::high_resolution_clock::now();
			CrankNicholson cn(x_dom, J, t_dom, N, bsPDE);
			cn.stepMarch();
			auto end = std::chrono::high_resolution_clock::now();
			chrono::duration<double> sure = end - begin;
			optimalCN << xSize << "," << tSize << "," << cn.newResult[(xSize / 2)] << "," << AnalyticalSol << "," << abs(cn.newResult[(xSize / 2)] - AnalyticalSol) / AnalyticalSol << "," << sure.count() << endl;
		}
	}
	optimalCN.close();

	delete bsPDE;

}