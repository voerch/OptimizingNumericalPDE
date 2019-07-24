#include "FDM.h"
#include <fstream>


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
void ADI::calculateStepSize()
{
	xStepSize = xDomain / static_cast<double>(xNumberSteps - 1);
	yStepSize = yDomain / static_cast<double>(yNumberSteps - 1);
	tStepSize = tDomain / static_cast<double>(tNumberSteps - 1);
	rX = tStepSize / (2.0 * xStepSize * xStepSize);
	rY = tStepSize / (2.0 * yStepSize * yStepSize);
}

void ADI::setInitialConditions()
{
	double currentX = 0;
	double currentY = 0;

	oldResult.resize(xNumberSteps, 0.0);
	newResult.resize(xNumberSteps, 0.0);
	xValues.resize(xNumberSteps, 0.0);
	yValues.resize(yNumberSteps, 0.0);
	
	for (long xCounter = 0; xCounter < xNumberSteps; xCounter++)
	{
		currentX = static_cast<double>(xCounter) * xStepSize;
		currentY = static_cast<double>(xCounter) * yStepSize;
		xValues[xCounter] = currentX;
		yValues[xCounter] = currentY;
		newResult[xCounter] = PDE->InitCond(currentX, currentY);
	}
	
	tPrevious = 0;
	tCurrent = 0;

}

void ADI::calculateBoundaryConditions()
{
	newResult[0] = PDE->BoundaryLeft(tPrevious, xValues[0], yValues[0]);
	newResult[xNumberSteps - 1] = PDE->BoundaryRight(tPrevious, xValues[xNumberSteps - 1], yValues[yNumberSteps - 1]);

}

// Loops through x values on a given time.
void ADI::calculateInnerDomain()
{
	alpha = -rX;
	beta = 1.0 + (2.0 * rX);
	gamma = -rX;

	LowerDiag.resize(xNumberSteps - 1, alpha);
	Diag.resize(xNumberSteps, beta);
	UpperDiag.resize(xNumberSteps - 1, gamma);
	
	//new result becomes n + 1/2
	for (long iCounter = 1; iCounter < xNumberSteps - 1; iCounter++)
	{		
		oldResult[iCounter] = rY * newResult[iCounter + 1] + (1 - 2.0 *rY) * newResult[iCounter] + (rY) * newResult[iCounter - 1];
	}
	ThomasAlgorithm(LowerDiag, Diag, UpperDiag, oldResult, newResult);

	alpha = -rY;
	beta = 1.0 + (2.0 * rY);
	gamma = -rY;

	LowerDiag.assign(yNumberSteps - 1, alpha);
	Diag.assign(yNumberSteps, beta);
	UpperDiag.assign(yNumberSteps - 1, gamma);
	//new result becomes n + 1
	for (long iCounter = 1; iCounter < xNumberSteps - 1; iCounter++)
	{
		oldResult[iCounter] = rX * newResult[iCounter + 1] + (1 - 2.0 *rX) * newResult[iCounter] + (rX)* newResult[iCounter - 1];
	}
	ThomasAlgorithm(LowerDiag, Diag, UpperDiag, oldResult, newResult);
}

// Loops through time.
void ADI::stepMarch()
{
	std::ofstream grid("ADIGrid.csv");
	//grid << "xValues,tValues,Solution" << std::endl;
	while (tCurrent < tDomain)
	{
		tCurrent = tPrevious + tStepSize;
		calculateBoundaryConditions();
		
		for (int xCounter = 0; xCounter < xNumberSteps; xCounter++)
		{
			for (int yCounter = 0; yCounter < yNumberSteps; yCounter++)
			{
				grid << xValues[xCounter] << "," << yValues[yCounter] << "," << tPrevious << "," << newResult[yNumberSteps - 1] << std::endl;
			}
			calculateInnerDomain();		
		}

		tPrevious = tCurrent;
	}

	grid.close();
}