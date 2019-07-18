#pragma once

// Second-order Parabolic PDE
class ParabolicPDE{
public:
	// PDE Coefficients 
	virtual double DiffusionCoeff(double t, double x) const = 0;
	virtual double ConvectionCoeff(double t, double x) const = 0;
	virtual double ZeroCoeff(double t, double x) const = 0;
	virtual double SourceCoeff(double t, double x) const = 0;

	// Boundary and initial conditions
	virtual double BoundaryLeft(double t, double x) const = 0;
	virtual double BoundaryRight(double t, double x) const = 0;
	virtual double InitCond(double x) const = 0;
};

// Implementing Heat Equation
class HeatEqn : public ParabolicPDE
{
	double DiffusionCoeff(double t, double x) const;
	double ConvectionCoeff(double t, double x) const;
	double ZeroCoeff(double t, double x) const;
	double SourceCoeff(double t, double x) const;

	double BoundaryLeft(double t, double x) const;
	double BoundaryRight(double t, double x) const;
	double InitCond(double x) const;
};

// Implementing Black Scholes PDE
//class BlackScholes : public ParabolicPDE
//{
//	double DiffusionCoeff(double t, double x) const;
//	double ConvectionCoeff(double t, double x) const;
//	double ZeroCoeff(double t, double x) const;
//	double SourceCoeff(double t, double x) const;
//
//	double BoundaryLeft(double t, double x) const;
//	double BoundaryRight(double t, double x) const;
//	double InitCond(double x) const;
//};