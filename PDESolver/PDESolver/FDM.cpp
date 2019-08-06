#include "FDM.h"
#include <fstream>
#include <iostream>
#include <random>
#include "TridiagonalSolvers.h"

using namespace std;

void ExplicitMethod::calculateStepSize() 
{
	xStepSize = xDomain / static_cast<double>(xNumberSteps - 1);
	tStepSize = tDomain / static_cast<double>(tNumberSteps-1);

	mt19937 RandomEngine(0);
	uniform_real_distribution<double> RandomSmallNumber(0.0, 0.000000000000000000001);
	xStepSize += RandomSmallNumber(RandomEngine);
	tStepSize += RandomSmallNumber(RandomEngine);
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
	grid << "xValues,tValues,Solution" << std::endl;

	calculateBoundaryConditions();
	for (int outputCounter = 0; outputCounter < xNumberSteps; outputCounter++)
	{
		grid << xValues[outputCounter] << "," << tCurrent << "," << newResult[outputCounter] << std::endl;
	}

	for (tCurrent = tStepSize; tCurrent < tDomain + tStepSize; tCurrent += tStepSize)
	{
		calculateBoundaryConditions();
		calculateInnerDomain();
		for (int outputCounter = 0; outputCounter < xNumberSteps; outputCounter++)
		{
			grid << xValues[outputCounter] << "," << tCurrent << "," << newResult[outputCounter] << std::endl;
		}

	}

	grid.close();
}



//Defining functions for CN Method
void CrankNicholson::calculateStepSize()
{
	xStepSize = xDomain / static_cast<double>(xNumberSteps - 1);
	tStepSize = tDomain / static_cast<double>(tNumberSteps - 1);

	mt19937 RandomEngine(0);
	uniform_real_distribution<double> RandomSmallNumber(0.0, 0.000000000000000000001);
	xStepSize += RandomSmallNumber(RandomEngine);
	tStepSize += RandomSmallNumber(RandomEngine);

	r = tStepSize / (xStepSize * xStepSize);
}

void CrankNicholson::setInitialConditions()
{
	double currentX = 0;

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
	grid << "xValues,tValues,Solution" << std::endl;

	calculateBoundaryConditions();
	for (int outputCounter = 0; outputCounter < xNumberSteps; outputCounter++)
	{
		grid << xValues[outputCounter] << "," << tCurrent << "," << newResult[outputCounter] << std::endl;
	}

	for (tCurrent = tStepSize; tCurrent < tDomain + tStepSize; tCurrent += tStepSize)
	{
		calculateBoundaryConditions();
		calculateInnerDomain();
		for (int outputCounter = 0; outputCounter < xNumberSteps; outputCounter++)
		{
			grid << xValues[outputCounter] << "," << tCurrent << "," << newResult[outputCounter] << std::endl;
		}

	}

	grid.close();
}

//Defining functions for CN Method
void ADI::StepSize()
{
	xStepSize = xDomain / static_cast<double>(xNumberSteps - 1);
	yStepSize = yDomain / static_cast<double>(yNumberSteps - 1);
	tStepSize = tDomain / static_cast<double>(tNumberSteps - 1);
	rX = tStepSize / (2.0 * xStepSize * xStepSize);
	rY = tStepSize / (2.0 * yStepSize * yStepSize);
}

void ADI::InitialConditions()
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

void ADI::BoundaryConditions()
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
void ADI::InnerDomain()
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
		BoundaryConditions();
		InnerDomain();
		tCurrent += tStepSize;
	}

	grid.close();
}