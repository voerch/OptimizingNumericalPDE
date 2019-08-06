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
	

		// Constructor
		FDM(double xDomain_, long xNumberSteps_, double tDomain_, long tNumberSteps_, ParabolicPDE* PDE_) : 
			xDomain(xDomain_), xNumberSteps(xNumberSteps_), tDomain(tDomain_), tNumberSteps(tNumberSteps_), PDE(PDE_) {};

		// Override these virtual methods in derived classes for 
		// specific FDM techniques, such as explicit Euler, Crank-Nicolson, etc.
		virtual void calculateStepSize() = 0;
		virtual void setInitialConditions() = 0;
		virtual void calculateBoundaryConditions() = 0;
		virtual void calculateInnerDomain() = 0;

		std::vector<double> oldResult; // N

	public:
		std::vector<double> newResult; // N+1
		
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

	std::vector<double> LowerDiag;
	std::vector<double> Diag;
	std::vector<double> UpperDiag;


public:
	CrankNicholson(double xDomain_, long xNumberSteps_, double tDomain_, long tNumberSteps_, ParabolicPDE* PDE_)
		: FDM(xDomain_, xNumberSteps_, tDomain_, tNumberSteps_, PDE_)
	{
		calculateStepSize();
		setInitialConditions();
	}

	void stepMarch();
};


class FDM2D
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
	FDM2D(double xDomain_, long xNumberSteps_, double yDomain_, long yNumberSteps_, double tDomain_, long tNumberSteps_, ParabolicPDE2D* PDE_) :
		xDomain(xDomain_), xNumberSteps(xNumberSteps_), yDomain(xDomain_), yNumberSteps(xNumberSteps_), tDomain(tDomain_), tNumberSteps(tNumberSteps_), PDE(PDE_) {};

	// Override these virtual methods in derived classes for 
	// specific FDM techniques, such as explicit Euler, Crank-Nicolson, etc.
	virtual void StepSize() = 0;
	virtual void InitialConditions() = 0;
	virtual void BoundaryConditions() = 0;
	virtual void InnerDomain() = 0;


public:
	// Carry out the actual time-stepping
	virtual void stepMarch() = 0;
};



class ADI : public FDM2D
{
protected:
	void StepSize();
	void InitialConditions();
	void BoundaryConditions();
	void InnerDomain();

	std::vector<double> LowerDiag;
	std::vector<double> Diag;
	std::vector<double> UpperDiag;

	std::vector<std::vector<double>> HalfStep;
	std::vector<std::vector<double>> FullStep;

public:
	ADI(double xDomain_, long xNumberSteps_, double yDomain_, long yNumberSteps_, double tDomain_, long tNumberSteps_, ParabolicPDE2D* PDE_)
		: FDM2D(xDomain_, xNumberSteps_, yDomain_, yNumberSteps_, tDomain_, tNumberSteps_, PDE_)
	{
		StepSize();
		InitialConditions();

	}

	void stepMarch();
};