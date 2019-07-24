#include "PDE2D.h"
#include <cmath>

double HeatEqn::DiffusionCoeff(double t, double x, double y) const
{
	return 1;
}
double HeatEqn::SecondDiffusionCoeff(double t, double x, double y) const
{
	return 1;
}

double HeatEqn::BoundaryLeft(double t, double x, double y) const
{
	return 0;
}
double HeatEqn::BoundaryRight(double t, double x, double y) const
{
	return 0;
}
double HeatEqn::InitCond(double x, double y) const
{
	//return sin(3.14159265358979323846 * x);
	return 1.0;
}
