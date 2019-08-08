#pragma once
#include <vector>
#include "mkl.h"
#include <fftw_mpi.h>


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

#include <cmath>

void cpu_cr(float *a, float *b, float *c, float *F, float *x, int size) {
	int i, j, index1, index2, offset;
	float k1, k2;


	/*Part 1 - Forward Reduction */
	for (i = 0; i<log2(size + 1) - 1; i++) {
		for (j = pow(2.0, i + 1) - 1; j<size; j = j + pow(2.0, i + 1)) {
			offset = pow(2.0, i);
			index1 = j - offset;
			index2 = j + offset;

			k1 = a[j] / b[index1];
			k2 = c[j] / b[index2];

			if (j == size - 1) {
				k1 = a[j] / b[j - offset];
				b[j] = b[j] - c[j - offset] * k1;
				F[j] = F[j] - F[j - offset] * k1;
				a[j] = -a[j - offset] * k1;
				c[j] = 0.0;
			}
			else {
				k1 = a[j] / b[j - offset];
				k2 = c[j] / b[j + offset];
				b[j] = b[j] - c[j - offset] * k1 - a[j + offset] * k2;
				F[j] = F[j] - F[j - offset] * k1 - F[j + offset] * k2;
				a[j] = -a[j - offset] * k1;
				c[j] = -c[j + offset] * k2;
			}
		}
	}



	/*part 2 - find the middle  */
	int index = (size - 1) / 2;
	x[index] = F[index] / b[index];

	/*part 3 - back substitution */
	for (i = log2(size + 1) - 2; i >= 0; i--) {
		for (j = pow(2.0, i + 1) - 1; j<size; j = j + pow(2.0, i + 1)) {
			offset = pow(2.0, i);
			index1 = j - offset;
			index2 = j + offset;


			if (j != index1) {
				if (index1 - offset < 0) x[index1] = (F[index1] - c[index1] * x[index1 + offset]) / b[index1];
				else x[index1] = (F[index1] - a[index1] * x[index1 - offset] - c[index1] * x[index1 + offset]) / b[index1];
			}if (j != index2) {
				if (index2 + offset >= size) x[index2] = (F[index2] - a[index2] * x[index2 - offset]) / b[index2];
				else x[index2] = (F[index2] - a[index2] * x[index2 - offset] - c[index2] * x[index2 + offset]) / b[index2];
			}

		}
	}
}
