#pragma once

// Second-order Parabolic PDE
class ParabolicPDE2D{
public:
	// PDE Coefficients 
	virtual double DiffusionCoeff(double t, double x, double y) const = 0;
	virtual double SecondDiffusionCoeff(double t, double x, double y) const = 0;
		
	//virtual double ConvectionCoeff(double t, double x, double y) const = 0;
	//virtual double ZeroCoeff(double t, double x, double y) const = 0;
	//virtual double SourceCoeff(double t, double x, double y) const = 0;

	// Boundary and initial conditions
	virtual double BoundaryLeft(double t, double x, double y) const = 0;
	virtual double BoundaryRight(double t, double x, double y) const = 0;
	virtual double InitCond(double x, double y) const = 0;
};

// Implementing Heat Equation
class HeatEqn : public ParabolicPDE2D
{
	double DiffusionCoeff(double t, double x, double y) const;
	double SecondDiffusionCoeff(double t, double x, double y) const;


	//double ConvectionCoeff(double t, double x, double y) const;
	//double ZeroCoeff(double t, double x, double y) const;
	//double SourceCoeff(double t, double x, double y) const;

	double BoundaryLeft(double t, double x, double y) const;
	double BoundaryRight(double t, double x, double y) const;
	double InitCond(double x, double y) const;

};