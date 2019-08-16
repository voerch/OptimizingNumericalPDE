#pragma once
#include <vector>
#include <cmath>
#include <omp.h>
#include <mkl.h>
#include <mkl_scalapack.h>


void ThomasAlgorithm(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d, std::vector<double>& f)
{
	int N = d.size();

	// Create the temporary vectors                                                                                                                                                                                    
	std::vector<double> c_star(N, 0.0);
	std::vector<double> d_star(N, 0.0);

	// This updates the coefficients in the first row                                                                                                                                                                  
	c_star[0] = c[0] / b[0];
	d_star[0] = d[0] / b[0];

	// Create the c_star and d_star coefficients in the forward sweep                                                                                                                                                  
	for (int i = 1; i<N - 1; i++) {
		double m = 1.0 / (b[i] - a[i] * c_star[i - 1]);
		c_star[i] = c[i] * m;
		d_star[i] = (d[i] - a[i] * d_star[i - 1]) * m;
	}

	// This is the reverse sweep, used to update the solution vector f                                                                                                                                                 
	for (int i = N - 1; i-- > 0; ) {
		f[i] = d_star[i] - c_star[i] * d[i + 1];
	}
}

// Documentation for Intel Math Kernel Library: 
// https://software.intel.com/en-us/mkl-developer-reference-c-gtsv
// https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/

void IntelSolver(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d, std::vector<double>& f)
{
	int size = a.size();

	//Pointers to vectors as the function requires double arrays.
	double * LowerDiag = &a[0];
	double * Diag = &b[0];
	double * UpperDiag = &c[0];
	double * oldResult = &d[0];

	//Specifies whether matrix storage layout is row major
	int matrix_layout = 102;

	//The order of A.
	lapack_int n = f.size();

	//The number of right-hand sides, the number of columns in B
	lapack_int nrhs = 1;

	//The leading dimension of b
	lapack_int ldb = f.size();

	LAPACKE_dgtsv(matrix_layout, n, nrhs, LowerDiag, Diag, UpperDiag, oldResult, ldb);

}

////Implementation of cyclic reduction from the package ScaLAPACK
//// https://software.intel.com/en-us/mkl-developer-reference-c-scalapack-routines
//// https://software.intel.com/en-us/mkl-developer-reference-c-p-dbtrs
//// https://software.intel.com/en-us/mkl-developer-reference-fortran-p-dbtrf
//void IntelParallelSolver(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d, std::vector<double>& f)
//{
//	//Pointers to vectors as the function requires double arrays.
//	double * LowerDiag = &a[0];
//	double * Diag = &b[0];
//	double * UpperDiag = &c[0];
//	double * oldResult = &d[0];
//	
//	char trans = 'N';
//	MKL_INT n = a.size();
//
//
//	//Number of upper and lower diagonals of A.
//	MKL_INT bwl = 1;
//	MKL_INT bwu = 1;
//
//	//The number of right-hand sides, the number of columns in B
//	MKL_INT nrhs = 1;
//
//	
//	double * ab;
//
//	MKL_INT ja;
//	MKL_INT desca;
//
//	MKL_INT ib = 0;
//	MKL_INT descb = d.size();
//
//	double af;
//	MKL_INT laf;
//	
//	double * work = &f[0];
//	MKL_INT lwork = f.size();
//
//	MKL_INT info;
//
//	pddbtrs(&trans, &n, &bwl, &bwu, &nrhs, ab, &ja, &desca, oldResult, &ib, &descb, &af, &laf, work, &lwork, &info);
//
//}

void CyclicReduction(std::vector<double> LDiag, std::vector<double> Diag, std::vector<double> UDiag, std::vector<double>  oldResult, std::vector<double>& newResult)
{
	double Temp1, Temp2;
	int iReducePlus, iReduceMinus, offset, backSubPlus, backSubMinus;
	int size = newResult.size();
	int iReduce, Step;
	LDiag.insert(LDiag.begin(), 0.0);
	UDiag.push_back(0.0);
	
//Forward Reduction phase
#pragma omp parallel for collapse(2)
	for (Step = 1; Step <= int(log2(size)); Step++)
	{
		for (int iReduce = pow(2, Step) - 1; iReduce < Diag.size(); iReduce += pow(2, Step))
		{
			offset = pow(2, Step - 1);
			iReduceMinus = iReduce - offset;
			iReducePlus = iReduce + offset;

			//if (iReduce == Diag.size() - 1)
			//{
			//	Temp1 = LDiag[iReduce] / Diag[iReduceMinus];
			//	LDiag[iReduce] = -LDiag[iReduceMinus] * Temp1;
			//	UDiag[iReduce] = 0;
			//	Diag[iReduce] = Diag[iReduce] - (UDiag[iReduceMinus] * Temp1);

			//	oldResult[iReduce] = oldResult[iReduce] - oldResult[iReduceMinus] * Temp1;

			//	//newResult[iReduce] = (oldResult[iReduce] - (LDiag[iReduce] * newResult[iReduce - 1])) / Diag[iReduce];
			//}
			//else
			//{
			//	Temp1 = LDiag[iReduce] / Diag[iReduceMinus];
			//	Temp2 = UDiag[iReduce] / Diag[iReducePlus];

			//	LDiag[iReduce] = -LDiag[iReduceMinus] * Temp1;
			//	UDiag[iReduce] = -UDiag[iReducePlus] * Temp2;
			//	Diag[iReduce] = Diag[iReduce] - (LDiag[iReducePlus] * Temp2) - (UDiag[iReduceMinus] * Temp1);
			//	oldResult[iReduce] = oldResult[iReduce] - (oldResult[iReduceMinus] * Temp1) - (oldResult[iReducePlus] * Temp2);


			//}

			//newResult[iReduce] = (oldResult[iReduce] - (LDiag[iReduce] * newResult[iReduceMinus]) - (UDiag[iReduce] * newResult[iReducePlus])) / Diag[iReduce];


			if (iReduceMinus <= 0)
			{
				Temp2 = UDiag[iReduce] / Diag[iReducePlus];
				LDiag[iReduce] = 0;
				UDiag[iReduce] = -UDiag[iReducePlus] * Temp2;
				Diag[iReduce] = Diag[iReduce] - (LDiag[iReducePlus] * Temp2);
				oldResult[iReduce] = oldResult[iReduce] - (oldResult[iReducePlus] * Temp2);
			}
			else if (iReducePlus >= size)
			{
				Temp1 = LDiag[iReduce] / Diag[iReduceMinus];
				LDiag[iReduce] = -LDiag[iReduceMinus] * Temp1;
				UDiag[iReduce] = 0;
				Diag[iReduce] = Diag[iReduce] - (UDiag[iReduceMinus] * Temp1);
				oldResult[iReduce] = oldResult[iReduce] - oldResult[iReduceMinus] * Temp1;
			}
			else
			{
				Temp1 = LDiag[iReduce] / Diag[iReduceMinus];
				Temp2 = UDiag[iReduce] / Diag[iReducePlus];
				LDiag[iReduce] = -LDiag[iReduceMinus] * Temp1;
				UDiag[iReduce] = -UDiag[iReducePlus] * Temp2;
				Diag[iReduce] = Diag[iReduce] - (LDiag[iReducePlus] * Temp2) - (UDiag[iReduceMinus] * Temp1);
				oldResult[iReduce] = oldResult[iReduce] - (oldResult[iReduceMinus] * Temp1) - (oldResult[iReducePlus] * Temp2);
			}


		}


	}


//Backward Substitution phase


//newResult.back() = (oldResult.back()) / Diag.back();
//newResult[pow(2, Step - 1) - 1] = (oldResult[pow(2, Step - 1) - 1] - (UDiag[(2, Step - 1) - 1] * newResult.back())) / Diag.back();

#pragma omp parallel for collapse(2)
	for (int Step = log2(size); Step >= 0; Step--)
	{
		for (int backSub = pow(2, Step + 1) - 1; backSub < size; backSub += pow(2, Step))
		{
			offset = pow(2, Step - 1);
			backSubMinus = backSub - offset;
			backSubPlus = backSub + offset;

			if (backSubMinus <= 0)
			{
				newResult[backSub] = (oldResult[backSub] - (UDiag[backSub] * newResult[backSubPlus])) / Diag[backSub];
			}
			else if (backSubPlus >= size)
			{
				newResult[backSub] = (oldResult[backSub] - (newResult[backSubMinus] * LDiag[backSub])) / Diag[backSub];

			}
			else
			{
				newResult[backSub] = (oldResult[backSub] - (newResult[backSubMinus] * LDiag[backSub]) - (UDiag[backSub] * newResult[backSubPlus])) / Diag[backSub];
			}
		}
	}
	//LDiag.erase(LDiag.begin());
	//UDiag.pop_back();

}