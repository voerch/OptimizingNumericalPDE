#pragma once

#include "PDE2D.h"
#include <vector>

class FDM
{
	protected:
		ParabolicPDE2D* PDE;

		// Discretisating x
		double xDomain; //Spatial extent [0, xDomain]
		long xNumberSteps; // Number of steps for x
		double xStepSize; // Calculated by xDomain/xNumberSteps
		std::vector<double> xValues; //Stores coordinates of x

		// Discretisating y
		double yDomain; //Spatial eytent [0, yDomain]
		long yNumberSteps; // Number of steps for y
		double yStepSize; // Calculated by yDomain/yNumberSteps
		std::vector<double> yValues; //Stores coordinates of y


		// Discretisating time
		double tDomain; // [0, tDomain]
		long tNumberSteps; // Number of steps for t
		double tStepSize; // Calculated by tDomain/tNumberSteps

		double tPrevious, tCurrent;

		double rX; // tStepSize / (xStepSize * xStepSize)
		double rY; // tStepSize / (xStepSize * xStepSize)
		// Differencing coeffs
		double alpha, beta, gamma;

		std::vector<double> newResult; // N+1
		std::vector<double> oldResult; // N

		// Constructor
		FDM(double xDomain_, long xNumberSteps_, double yDomain_, long yNumberSteps_, double tDomain_, long tNumberSteps_, ParabolicPDE2D* PDE_) :
			xDomain(xDomain_), xNumberSteps(xNumberSteps_), yDomain(xDomain_), yNumberSteps(xNumberSteps_), tDomain(tDomain_), tNumberSteps(tNumberSteps_), PDE(PDE_) {};

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



class ADI : public FDM
{
protected:
	void calculateStepSize();
	void setInitialConditions();
	void calculateBoundaryConditions();
	void calculateInnerDomain();

	std::vector<double> LowerDiag;
	std::vector<double> Diag;
	std::vector<double> UpperDiag;

	std::vector<std::vector<double>> HalfStep;
	std::vector<std::vector<double>> FullStep;

public:
	ADI(double xDomain_, long xNumberSteps_, double yDomain_, long yNumberSteps_, double tDomain_, long tNumberSteps_, ParabolicPDE2D* PDE_)
		: FDM(xDomain_, xNumberSteps_, yDomain_, yNumberSteps_, tDomain_, tNumberSteps_, PDE_)
	{
		calculateStepSize();
		setInitialConditions();

	}

	void stepMarch();
};