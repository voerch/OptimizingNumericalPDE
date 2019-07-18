#include "PDE.h"
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
}