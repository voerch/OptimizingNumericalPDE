#include "FDM2D.h"
#include <fstream>
using namespace std;

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
	
	HalfStep.resize(xNumberSteps, vector<double>(yNumberSteps, 0));
	FullStep.resize(xNumberSteps, vector<double>(yNumberSteps, 0));


	for (long xCounter = 0; xCounter < xNumberSteps; xCounter++)
	{
		currentX = static_cast<double>(xCounter) * xStepSize;
		xValues[xCounter] = currentX;
		for (long yCounter = 0; yCounter < yNumberSteps; yCounter++)
		{
			currentY = static_cast<double>(yCounter) * yStepSize;
			yValues[yCounter] = currentY;
			HalfStep[xCounter][yCounter] = 1;
			FullStep[xCounter][yCounter] = 1;
		}
	}
	
	tPrevious = 0;
	tCurrent = 0;

}

void ADI::calculateBoundaryConditions()
{
	for (long xCounter = 1; xCounter < xNumberSteps; xCounter++)
	{
		HalfStep[0][xCounter] = 0;
		HalfStep[xNumberSteps - 1][xCounter] = 0;
		HalfStep[xCounter][0] = 0;
		HalfStep[xCounter][xNumberSteps - 1] = 0;

		FullStep[0][xCounter] = 0;
		FullStep[xNumberSteps - 1][xCounter] = 0;
		FullStep[xCounter][0] = 0;
		FullStep[xCounter][xNumberSteps - 1] = 0;
	}
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

	//new result becomes n + 1/2 and fills half step.
	for (long yCounter = 1; yCounter < yNumberSteps - 1; yCounter++)
	{
		for (long xCounter = 1; xCounter < xNumberSteps - 1; xCounter++)
		{
			oldResult[xCounter] = rY * FullStep[xCounter][yCounter + 1] + (1 - 2.0 *rY) * FullStep[xCounter][yCounter] + (rY)* FullStep[xCounter][yCounter - 1];
		}

		ThomasAlgorithm(LowerDiag, Diag, UpperDiag, oldResult, newResult);

		for (long xCounter = 1; xCounter < xNumberSteps - 1; xCounter++)
		{
			HalfStep[xCounter][yCounter] = newResult[xCounter];
		}
	}

	alpha = -rY;
	beta = 1.0 + (2.0 * rY);
	gamma = -rY;

	LowerDiag.assign(yNumberSteps - 1, alpha);
	Diag.assign(yNumberSteps, beta);
	UpperDiag.assign(yNumberSteps - 1, gamma);
	//new result becomes n + 1
	for (long xCounter = 1; xCounter < xNumberSteps - 1; xCounter++)
	{
		for (long yCounter = 1; yCounter < yNumberSteps - 1; yCounter++)
		{
			oldResult[yCounter] = rX * FullStep[xCounter][yCounter + 1] + (1 - 2.0 *rX) * FullStep[xCounter][yCounter] + (rX)* FullStep[xCounter][yCounter - 1];
		}

		ThomasAlgorithm(LowerDiag, Diag, UpperDiag, oldResult, newResult);

		for (long yCounter = 1; yCounter < yNumberSteps - 1; yCounter++)
		{
			FullStep[xCounter][yCounter] = newResult[xCounter];
		}
	}
}

// Loops through time.
void ADI::stepMarch()
{
	std::ofstream grid("ADIGrid.csv");
	//grid << "xValues,tValues,Solution" << std::endl;
	while (tCurrent <= tDomain)
	{
		for (int xCounter = 0; xCounter < xNumberSteps; xCounter++)
		{
			for (int yCounter = 0; yCounter < yNumberSteps; yCounter++)
			{
				grid << xValues[xCounter] << "," << yValues[yCounter] << "," << tCurrent << "," << FullStep[xCounter][yCounter] << std::endl;
			}
		}
		calculateBoundaryConditions();
		calculateInnerDomain();
		tCurrent += tStepSize;
	}

	grid.close();
}