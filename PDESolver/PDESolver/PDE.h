/*
Author: Mustafa Berke Erdis, August 2019

Header file for storing parabolic partial differential equations efficiently.

*/
#pragma once
#include "Option.h"

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
	virtual double AnalyticSol(double t, double x) const = 0;
};

// Implementing 1D Heat Equation
class HeatEqn : public ParabolicPDE
{
public:
	double DiffusionCoeff(double t, double x) const;
	double ConvectionCoeff(double t, double x) const;
	double ZeroCoeff(double t, double x) const;
	double SourceCoeff(double t, double x) const;

	double BoundaryLeft(double t, double x) const;
	double BoundaryRight(double t, double x) const;
	double InitCond(double x) const;
	double AnalyticSol(double t, double x) const;
};

// Implementing Black Scholes PDE
class BlackScholesPDE : public ParabolicPDE
{
public:
	VanillaOption* option;
	BlackScholesPDE(VanillaOption* option_) : option(option_) {};

	double DiffusionCoeff(double t, double x) const;
	double ConvectionCoeff(double t, double x) const;
	double ZeroCoeff(double t, double x) const;
	double SourceCoeff(double t, double x) const;

	double BoundaryLeft(double t, double x) const;
	double BoundaryRight(double t, double x) const;
	double InitCond(double x) const;
	double AnalyticSol (double t, double x) const;
};

// Second-order Parabolic PDE
class ParabolicPDE2D {
public:
	// PDE Coefficients 
	virtual double DiffusionCoeff2D(double t, double x, double y) const = 0;
	virtual double SecondDiffusionCoeff2D(double t, double x, double y) const = 0;

	virtual double BoundaryLeft2D(double t, double x, double y) const = 0;
	virtual double BoundaryRight2D(double t, double x, double y) const = 0;
	virtual double InitCond2D(double x, double y) const = 0;

	virtual double AnalyticSol(double x, double y) const = 0;
};

// Implementing 2D Heat Equation
class Heat2D : public ParabolicPDE2D
{
	double DiffusionCoeff2D(double t, double x, double y) const;
	double SecondDiffusionCoeff2D(double t, double x, double y) const;

	double BoundaryLeft2D(double t, double x, double y) const;
	double BoundaryRight2D(double t, double x, double y) const;
	double InitCond2D(double x, double y) const;
	double AnalyticSol(double x, double y) const;
};
