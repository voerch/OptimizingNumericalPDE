#include "FDM.h"
#include <fstream>

void ExplicitMethod::calculateStepSize() 
{
	xStepSize = xDomain / static_cast<double>(xNumberSteps-1);
	tStepSize = tDomain / static_cast<double>(tNumberSteps-1);
}

void ExplicitMethod::setInitialConditions() 
{
	double currentX = 0;

	oldResult.resize(xNumberSteps, 0.0);
	newResult.resize(xNumberSteps, 0.0);
	xValues.resize(xNumberSteps, 0.0);

	for (long xCounter = 0; xCounter < xNumberSteps; xCounter++)
	{
		currentX = static_cast<double>(xCounter) * xStepSize;
		oldResult[xCounter] = PDE->InitCond(currentX);
		xValues[xCounter] = currentX;
	}

	tPrevious = 0;
	tCurrent = 0;

}

void ExplicitMethod::calculateBoundaryConditions() 
{
	newResult[0] = PDE->BoundaryLeft(tPrevious, xValues[0]);
	newResult[xNumberSteps-1] = PDE->BoundaryRight(tPrevious, xValues[xNumberSteps-1]);
}

// Loops through x values on a given time.
void ExplicitMethod::calculateInnerDomain() 
{
	for (long xCounter = 1; xCounter < xNumberSteps - 1; xCounter++)
	{
		double TempVar = tStepSize * (PDE->DiffusionCoeff(tPrevious, xValues[xCounter]));
		double TempVarTwo = 0.5 * tStepSize * xStepSize * (PDE->ConvectionCoeff(tPrevious, xValues[xCounter]));
	
		alpha = TempVar - TempVarTwo;
		beta = (xStepSize * xStepSize) - (2.0 * TempVar) + (tStepSize * xStepSize * xStepSize * PDE->ZeroCoeff(tPrevious, xValues[xCounter]));
		gamma = TempVar + TempVarTwo;

		newResult[xCounter] = ((alpha * oldResult[xCounter - 1]) + (beta * oldResult[xCounter]) + (gamma * oldResult[xCounter + 1])) / (xStepSize * xStepSize) - (tStepSize * PDE->SourceCoeff(tPrevious, xValues[xCounter]));

	}
}

// Loops through time.
void ExplicitMethod::stepMarch()
{
	std::ofstream grid("ExplicitGrid.csv");
	//grid << "xValues,tValues,Solution" << std::endl;
	while (tCurrent < tDomain)
	{
		tCurrent = tPrevious + tStepSize;
		calculateBoundaryConditions();
		calculateInnerDomain();

		for (int outputCounter = 0; outputCounter < xNumberSteps; outputCounter++)
		{
			grid << xValues[outputCounter] << "," << tPrevious << "," << newResult[outputCounter] << std::endl;
		}

		oldResult = newResult;
		tPrevious = tCurrent;
	}

	grid.close();
}

void FDM::ThomasAlgorithm(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d, std::vector<double>& f)
{
	size_t N = d.size();

	// Create the temporary vectors                                                                                                                                                                                    
	std::vector<double> c_star(N, 0.0);
	std::vector<double> d_star(N, 0.0);

	// This updates the coefficients in the first row                                                                                                                                                                  
	c_star[0] = c[0] / b[0];
	d_star[0] = d[0] / b[0];

	// Create the c_star and d_star coefficients in the forward sweep                                                                                                                                                  
	for (int i = 1; i<N-1; i++) {
		double m = 1.0 / (b[i] - a[i] * c_star[i - 1]);
		c_star[i] = c[i] * m;
		d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m;
	}

	// This is the reverse sweep, used to update the solution vector f                                                                                                                                                 
	for (int i = N - 1; i-- > 0; ) {
		f[i] = d_star[i] - c_star[i] * d[i + 1];
	}
}


//Defining functions for CN Method
void CrankNicholson::calculateStepSize()
{
	xStepSize = xDomain / static_cast<double>(xNumberSteps - 1);
	tStepSize = tDomain / static_cast<double>(tNumberSteps - 1);
	r = tStepSize / (xStepSize * xStepSize);
}

void CrankNicholson::setInitialConditions()
{
	double currentX = 0;

	//a.resize(xNumberSteps - 1, -r / 2.0);
	//b.resize(xNumberSteps, 1.0 + r);
	//c.resize(xNumberSteps - 1, -r / 2.0);

	LowerDiag.resize(xNumberSteps - 1, 0.0);
	Diag.resize(xNumberSteps, 0.0);
	UpperDiag.resize(xNumberSteps - 1, 0.0);

	oldResult.resize(xNumberSteps, 0.0);
	newResult.resize(xNumberSteps, 0.0);
	xValues.resize(xNumberSteps, 0.0);

	for (long xCounter = 0; xCounter < xNumberSteps; xCounter++)
	{
		currentX = static_cast<double>(xCounter) * xStepSize;
		newResult[xCounter] = PDE->InitCond(currentX);
		xValues[xCounter] = currentX;
	}
	
	tPrevious = 0;
	tCurrent = 0;

}

void CrankNicholson::calculateBoundaryConditions()
{
	newResult[0] = PDE->BoundaryLeft(tPrevious, xValues[0]);
	newResult[xNumberSteps - 1] = PDE->BoundaryRight(tPrevious, xValues[xNumberSteps - 1]);

}

// Loops through x values on a given time.
void CrankNicholson::calculateInnerDomain()
{
	alpha = PDE->DiffusionCoeff(tPrevious, xValues[0]) * tStepSize / (2 * xStepSize * xStepSize);
	beta = PDE->ConvectionCoeff(tPrevious, xValues[0]) * tStepSize / (4.0 * xStepSize);
	gamma = PDE->ZeroCoeff(tPrevious, xValues[0]) * tStepSize * 0.5;

	LowerDiag[0] = -alpha - beta;
	Diag[0] = (1.0 + (2 * alpha) - gamma);
	UpperDiag[0] = -alpha + beta;

	for (long iCounter = 1; iCounter < xNumberSteps - 1; iCounter++)
	{
		alpha = PDE->DiffusionCoeff(tPrevious, xValues[iCounter]) * tStepSize / (2 * xStepSize * xStepSize);
		beta = PDE->ConvectionCoeff(tPrevious, xValues[iCounter]) * tStepSize / (4.0 * xStepSize);
		gamma = PDE->ZeroCoeff(tPrevious, xValues[iCounter]) * tStepSize * 0.5;

		LowerDiag[iCounter] = - alpha - beta;
		Diag[iCounter] = (1.0 + (2 * alpha) - gamma);
		UpperDiag[iCounter] = - alpha + beta ;
		
		oldResult[iCounter] = (alpha + beta) * newResult[iCounter + 1] + (1.0 - (2 * alpha) + gamma) * newResult[iCounter] + (alpha - beta) * newResult[iCounter - 1];
	}

	ThomasAlgorithm(LowerDiag, Diag, UpperDiag, oldResult, newResult);

}

// Loops through time.
void CrankNicholson::stepMarch()
{
	std::ofstream grid("CNGrid.csv");
	//grid << "xValues,tValues,Solution" << std::endl;
	while (tCurrent < tDomain)
	{
		tCurrent = tPrevious + tStepSize;
		calculateBoundaryConditions();
		
		for (int outputCounter = 0; outputCounter < xNumberSteps; outputCounter++)
		{
			grid << xValues[outputCounter] << "," << tPrevious << "," << newResult[outputCounter] << std::endl;
		}
		calculateInnerDomain();
		tPrevious = tCurrent;
	}

	grid.close();
}