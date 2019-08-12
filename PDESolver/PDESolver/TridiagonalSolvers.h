#pragma once
#include <vector>
#include <cmath>

#include <mkl.h>
#include <mkl_scalapack.h>


#include <omp.h>

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

//Solves A x = B using cyclic reduction from the package ScaLAPACK
// https://software.intel.com/en-us/mkl-developer-reference-c-scalapack-routines
// https://software.intel.com/en-us/mkl-developer-reference-c-p-dbtrs
//
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
//	double a;
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
//	pddbtrs(&trans, &n, &bwl, &bwu, &nrhs, &a, &ja, &desca, oldResult, &ib, &descb, &af, &laf, work, &lwork, &info);
//
//}
