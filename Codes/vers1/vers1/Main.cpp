#include <iostream>
#include <chrono>
#include "BSModel01.h"
#include "Option.h"

#include "BSEq.h"
#include "HeatEq.h"

#include "ExplicitMethod.h"
#include "CNMethod.h"

using namespace std;

int main()
{


	double S0 = 100.0, r = 0.05, sigma = 0.2;
	BSModel Model(S0, r, sigma);

	double T = 1. / 12., K = 100.0, zl = 1.0, zu = 2.0*S0;
	Put EuropeanPut(K, T, zl, zu);

	BSEq BSPDE(&Model,&EuropeanPut);

	int imax=3000, jmax=1000;

	auto start = std::chrono::high_resolution_clock::now();

	ExplicitMethod Method1(&BSPDE, imax, jmax);
	Method1.SolvePDE();
	cout << "Explicit Price = " << Method1.v(0.0, S0) << endl;

	auto finish = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed = finish - start;
	cout << "Elapsed time: " << elapsed.count() << endl;

	auto begin = std::chrono::high_resolution_clock::now();

	CNMethod Method2(&BSPDE, imax, jmax);
	Method2.SolvePDE();
	cout << "CN Price = " << Method2.v(0.0, S0) << endl;
	
	auto end = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed2 = end - begin;
	cout << "Elapsed time: " << elapsed2.count() << endl;


	imax = 200;
	jmax = 2000;
	auto go = std::chrono::high_resolution_clock::now();

	HeatEq HeatPDE(&Model, &EuropeanPut);

	CNMethod Method3(&HeatPDE, imax, jmax);
	Method3.SolvePDE();

	double t = 0.0;
	double z = S0;

	double x = HeatPDE.X(t, z);
	cout << "Heat Equation Price = " << HeatPDE.U(t, Method3.v(t, x)) << endl;
	auto stop = std::chrono::high_resolution_clock::now();
	chrono::duration<double> elapsed3 = stop - go;
	cout << "Elapsed time: " << elapsed3.count() << endl;



	system("pause");
	return 0;
}