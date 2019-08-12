#pragma once
#include "PDE.h"
#include <vector>

class FDM
{
	protected:
		ParabolicPDE* PDE;

		// Discretisating space
		double xDomain; //Spatial extent [0, x_max]
		long xNumberSteps; // Number of steps for x
		double xStepSize; // Calculated by xDomain/xNumberSteps
		std::vector<double> xValues; //Stores coordinates of x

		// Discretisating time
		double tDomain; // [0, T]
		long tNumberSteps; // Number of steps for t
		double tStepSize; // Calculated by tDomain/tNumberSteps

		double tPrevious, tCurrent;

		double r; // tStepSize / (xStepSize * xStepSize)
		
		// Differencing coefficients
		double alpha, beta, gamma;
	
		// Constructor
		FDM(double xDomain_, long xNumberSteps_, double tDomain_, long tNumberSteps_, ParabolicPDE* PDE_) : xDomain(xDomain_), xNumberSteps(xNumberSteps_), tDomain(tDomain_), tNumberSteps(tNumberSteps_), PDE(PDE_) {};

		// These virtual methods will be overriden in derived classes for explicit Euler scheme and Crank-Nicolson method.
		virtual void stepSize() = 0;
		virtual void initialConditions() = 0;
		virtual void boundaryConditions() = 0;

		// Loops through x values.
		virtual void innerDomain() = 0;

		std::vector<double> oldResult; // t = N

	public:
		std::vector<double> newResult; // t = N+1
		

		// Loops through time.
		virtual void timeMarch() = 0;
};

class ExplicitMethod : public FDM
{
	protected:
		void stepSize();
		void initialConditions();
		void boundaryConditions();
		void innerDomain();
	
	public:
		ExplicitMethod(double xDomain_, long xNumberSteps_,	double tDomain_, long tNumberSteps_, ParabolicPDE* PDE_) 
			: FDM(xDomain_, xNumberSteps_,	tDomain_, tNumberSteps_, PDE_) 
		{
			stepSize();
			initialConditions();
		}
		
		void timeMarch();
};


class CrankNicholson : public FDM
{
	protected:
		void stepSize();
		void initialConditions();
		void boundaryConditions();
		void innerDomain();

		std::vector<double> LowerDiag;
		std::vector<double> Diag;
		std::vector<double> UpperDiag;


	public:
		CrankNicholson(double xDomain_, long xNumberSteps_, double tDomain_, long tNumberSteps_, ParabolicPDE* PDE_)
			: FDM(xDomain_, xNumberSteps_, tDomain_, tNumberSteps_, PDE_)
		{
			stepSize();
			initialConditions();
		}

		void timeMarch();
};

// Numerical methods for 2D PDEs are more complicated hence FDM2D class was needed.
class FDM2D
{
	protected:
		ParabolicPDE2D* PDE;

		// Discretisating x
		double xDomain; //Spatial extent [0, x_max]
		long xNumberSteps; // Number of steps for x
		double xStepSize; // Calculated by xDomain/xNumberSteps
		std::vector<double> xValues; //Stores coordinates of x

									 // Discretisating y
		double yDomain; //Spatial eytent [0, y_max]
		long yNumberSteps; // Number of steps for y
		double yStepSize; // Calculated by yDomain/yNumberSteps
		std::vector<double> yValues; //Stores coordinates of y

									 // Discretisating time
		double tDomain; // [0, T]
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
		FDM2D(double xDomain_, long xNumberSteps_, double yDomain_, long yNumberSteps_, double tDomain_, long tNumberSteps_, ParabolicPDE2D* PDE_) : xDomain(xDomain_), xNumberSteps(xNumberSteps_), yDomain(xDomain_), yNumberSteps(xNumberSteps_), tDomain(tDomain_), tNumberSteps(tNumberSteps_), PDE(PDE_) {};

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
		ADI(double xDomain_, long xNumberSteps_, double yDomain_, long yNumberSteps_, double tDomain_, long tNumberSteps_, ParabolicPDE2D* PDE_) : FDM2D(xDomain_, xNumberSteps_, yDomain_, yNumberSteps_, tDomain_, tNumberSteps_, PDE_)
		{
			StepSize();
			InitialConditions();

		}
		void TimeMarch();
};