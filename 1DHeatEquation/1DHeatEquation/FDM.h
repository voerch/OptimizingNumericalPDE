#pragma once

#include "PDE.h"
#include <vector>

class FDM
{
	protected:
		ParabolicPDE* PDE;

		// Discretisating space
		double xDomain; //Spatial extent [0, xDomain]
		long xNumberSteps; // Number of steps for x
		double xStepSize; // Calculated by xDomain/xNumberSteps
		std::vector<double> xValues; //Stores coordinates of x

		// Discretisating time
		double tDomain; // [0, tDomain]
		long tNumberSteps; // Number of steps for t
		double tStepSize; // Calculated by tDomain/tNumberSteps

		double tPrevious, tCurrent;

		double r; // tStepSize / (xStepSize * xStepSize)
		// Differencing coeffs
		double alpha, beta, gamma;

		std::vector<double> newResult; // N+1
		std::vector<double> oldResult; // N

		// Constructor
		FDM(double xDomain_, long xNumberSteps_, double tDomain_, long tNumberSteps_, ParabolicPDE* PDE_) : 
			xDomain(xDomain_), xNumberSteps(xNumberSteps_), tDomain(tDomain_), tNumberSteps(tNumberSteps_), PDE(PDE_) {};

		// Override these virtual methods in derived classes for 
		// specific FDM techniques, such as explicit Euler, Crank-Nicolson, etc.
		virtual void calculateStepSize() = 0;
		virtual void setInitialConditions() = 0;
		virtual void calculateBoundaryConditions() = 0;
		virtual void calculateInnerDomain() = 0;

		void ThomasAlgorithm(const std::vector<double>& a,	const std::vector<double>& b, const std::vector<double>& c,	const std::vector<double>& d, std::vector<double>& f);

	public:
		// Carry out the actual time-stepping
		virtual void stepMarch() = 0;
};

class ExplicitMethod : public FDM
{
	protected:
		void calculateStepSize();
		void setInitialConditions();
		void calculateBoundaryConditions();
		void calculateInnerDomain();
	
	public:
		ExplicitMethod(double xDomain_, long xNumberSteps_,	double tDomain_, long tNumberSteps_, ParabolicPDE* PDE_) 
			: FDM(xDomain_, xNumberSteps_,	tDomain_, tNumberSteps_, PDE_) 
		{
			calculateStepSize();
			setInitialConditions();
		}
		
		void stepMarch();
};


class CrankNicholson : public FDM
{
protected:
	void calculateStepSize();
	void setInitialConditions();
	void calculateBoundaryConditions();
	void calculateInnerDomain();


public:
	CrankNicholson(double xDomain_, long xNumberSteps_, double tDomain_, long tNumberSteps_, ParabolicPDE* PDE_)
		: FDM(xDomain_, xNumberSteps_, tDomain_, tNumberSteps_, PDE_)
	{
		calculateStepSize();
		setInitialConditions();


		

	}

	void stepMarch();
};