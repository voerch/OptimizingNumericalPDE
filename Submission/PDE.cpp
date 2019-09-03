/*
Author: Mustafa Berke Erdis, August 2019

This file consists of Black - Scholes equation, one-dimensional and two dimensional heat equation implementations 

*/

#include "PDE.h"
#include "Option.h"
#include <cmath>

//Implementation of the heat equation
double HeatEqn::DiffusionCoeff(double t, double x) const
{
	return 1;
}
double HeatEqn::ConvectionCoeff(double t, double x) const
{
	return 0;
}
double HeatEqn::ZeroCoeff(double t, double x) const
{
	return 0;
}
double HeatEqn::SourceCoeff(double t, double x) const
{
	return 0;
}
double HeatEqn::BoundaryLeft(double t, double x) const
{
	return 0;
}
double HeatEqn::BoundaryRight(double t, double x) const
{
	return 0;
}
double HeatEqn::InitCond(double x) const
{
	return sin(PI * x);
}

double HeatEqn::AnalyticSol(double t, double x) const
{
	return exp(-t * pow(PI, 2)) * sin(PI * x);
}

// Black - Scholes partial differential equation implementation.
double BlackScholesPDE::DiffusionCoeff(double t, double x) const
{
	return 0.5 * option->sigma * option->sigma * x * x;
}
double BlackScholesPDE::ConvectionCoeff(double t, double x) const
{
	return option->interestRate * x;
}
double BlackScholesPDE::ZeroCoeff(double t, double x) const
{
	return -option->interestRate;
}
double BlackScholesPDE::SourceCoeff(double t, double x) const
{
	return 0.0;
}
double BlackScholesPDE::BoundaryLeft(double t, double x) const
{
	return 0.0;
}
double BlackScholesPDE::BoundaryRight(double t, double x) const
{
	return (x - (option->strike)*exp(-(option->interestRate)*((option->timeToExpiry) - t)));
}
double BlackScholesPDE::InitCond(double x) const
{
	return option->Payoff(x);
}

double BlackScholesPDE::AnalyticSol(double t, double x) const
{
	return option->PriceByBS(x);
}

// 2D Heat equation implementation
double Heat2D::DiffusionCoeff2D(double t, double x, double y) const
{
	return 1;
}
double Heat2D::SecondDiffusionCoeff2D(double t, double x, double y) const
{
	return 1;
}

double Heat2D::BoundaryLeft2D(double t, double x, double y) const
{
	return 0;
}
double Heat2D::BoundaryRight2D(double t, double x, double y) const
{
	return 0;
}
double Heat2D::InitCond2D(double x, double y) const
{
	//return sin(3.14159265358979323846 * x);
	return 1.0;
}
double Heat2D::AnalyticSol(double x, double y) const
{
	return  0;
		//(32 / (pow(PI, 3) * sinh(PI / 2))) * sin((PI * x) / 2)  * sinh((PI * y) / 2);
}