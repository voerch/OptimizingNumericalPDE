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
	//return 1.0;
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


// 2D
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