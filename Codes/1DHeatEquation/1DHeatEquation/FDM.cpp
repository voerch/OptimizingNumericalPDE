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


//Defining functions for CN Method
void CrankNicholson::calculateStepSize()
{
	xStepSize = xDomain / static_cast<double>(xNumberSteps - 1);
	tStepSize = tDomain / static_cast<double>(tNumberSteps - 1);
}

void CrankNicholson::setInitialConditions()
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

void CrankNicholson::calculateBoundaryConditions()
{
	newResult[0] = PDE->BoundaryLeft(tPrevious, xValues[0]);
	newResult[xNumberSteps - 1] = PDE->BoundaryRight(tPrevious, xValues[xNumberSteps - 1]);
}

// Loops through x values on a given time.
void CrankNicholson::calculateInnerDomain()
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
void CrankNicholson::stepMarch()
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