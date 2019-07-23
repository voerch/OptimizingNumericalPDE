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
	return sin(3.14159265358979323846 * x);
	//return 1;
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