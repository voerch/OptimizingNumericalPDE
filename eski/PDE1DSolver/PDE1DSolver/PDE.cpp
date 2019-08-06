#include "PDE.h"
#include "Option.h"
#include <cmath>

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

double LogSpotBlackScholesPDE::DiffusionCoeff(double t, double x) const
{
	return 0.5 * option->sigma * option->sigma;
}
double LogSpotBlackScholesPDE::ConvectionCoeff(double t, double x) const
{
	return option->interestRate - (0.5 * option->sigma * option->sigma);
}
double LogSpotBlackScholesPDE::ZeroCoeff(double t, double x) const
{
	return -option->interestRate;
}
double LogSpotBlackScholesPDE::SourceCoeff(double t, double x) const
{
	return 0.0;
}
double LogSpotBlackScholesPDE::BoundaryLeft(double t, double x) const
{
	return 0.0;
}
double LogSpotBlackScholesPDE::BoundaryRight(double t, double x) const
{
	return (exp(x) - (option->strike)*exp(-(option->interestRate)*((option->timeToExpiry) - t)));
}
double LogSpotBlackScholesPDE::InitCond(double x) const
{
	return option->Payoff(exp(x));
}
double LogSpotBlackScholesPDE::AnalyticSol(double t, double x) const
{
	return option->PriceByBS(exp(x));
}
