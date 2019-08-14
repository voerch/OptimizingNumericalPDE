#pragma once
#include <vector>
#include <cmath>

//#include <mkl.h>
//#include <mkl_scalapack.h>

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
//
//// Documentation for Intel Math Kernel Library: 
//// https://software.intel.com/en-us/mkl-developer-reference-c-gtsv
//// https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/
//
//void IntelSolver(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d, std::vector<double>& f)
//{
//	int size = a.size();
//
//	//Pointers to vectors as the function requires double arrays.
//	double * LowerDiag = &a[0];
//	double * Diag = &b[0];
//	double * UpperDiag = &c[0];
//	double * oldResult = &d[0];
//
//	//Specifies whether matrix storage layout is row major
//	int matrix_layout = 102;
//
//	//The order of A.
//	lapack_int n = f.size();
//
//	//The number of right-hand sides, the number of columns in B
//	lapack_int nrhs = 1;
//
//	//The leading dimension of b
//	lapack_int ldb = f.size();
//
//	LAPACKE_dgtsv(matrix_layout, n, nrhs, LowerDiag, Diag, UpperDiag, oldResult, ldb);
//
//	delete[] LowerDiag;
//	delete[] Diag;
//	delete[] UpperDiag;
//	delete[] oldResult;
//}

//
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

//
//#define OMP_SCHEDULE dynamic
//#define CHUNK 32
//
//// Implementation of cyclic reduction using OpenMP
//std::vector<double> CyclicReduction(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d, std::vector<double>& f)
//{
//
//	int equation_count = b.size();
//	int reductions_count = static_cast<int>(floor(log(equation_count) /	log(2)));
//	int coeff_array_size = 2 * (equation_count + 1) - reductions_count - 3;
//	
//	double* a_vecs = &a[0];
//	double* b_vecs = &b[0];
//	double* c_vecs = &c[0];
//	double* rhs_vecs = &d[0];
//
//	int base_idx = equation_count;
//
//	while (equation_count > 1)
//	{
//		int old_base_idx = base_idx - equation_count;
//		equation_count /= 2;
//
//		double* old_a = a_vecs + old_base_idx;
//		double* old_b = b_vecs + old_base_idx;
//		double* old_c = c_vecs + old_base_idx;
//		double* old_rhs = rhs_vecs + old_base_idx;
//		double* new_a = a_vecs + base_idx;
//		double* new_b = b_vecs + base_idx;
//		double* new_c = c_vecs + base_idx;
//		double* new_rhs = rhs_vecs + base_idx;
//
//#pragma omp parallel for if (equation_count > 16)
//		for (int i = 0; i < equation_count; i++)
//		{
//			int old_i = 2 * i + 1;
//			double oldA_next = old_a[old_i + 1];
//			double oldA_prev = old_a[old_i - 1];
//
//			double CA = old_c[old_i] / oldA_next;
//			double BA = old_b[old_i] / oldA_prev;
//
//			new_a[i] = old_a[old_i] - BA * old_c[old_i - 1] - CA * old_b[old_i + 1];
//			new_b[i] = -old_b[old_i - 1] * BA;
//			new_c[i] = -CA * old_c[old_i + 1];
//			new_rhs[i] = old_rhs[old_i] - old_rhs[old_i - 1] * BA -
//			old_rhs[old_i + 1] * CA;
//		}
//
//		base_idx += equation_count;
//	}
//
//	// result is padded with one zero one both sides to emulate
//	// implicit variables x_0 and x_N
//	std::vector<double> result(b.size() + 2, 0.0);
//	// factor maps indices of equations in full system extended with x_0=0
//	// and x_N=0 onto indices of equations in reduced extended system
//	int factor = b.size() / 2 + 1;
//	for (unsigned int i = reductions_count + 1; i-- > 0;)
//	{
//		base_idx -= equation_count;
//		double* as = a_vecs + base_idx;
//		double* bs = b_vecs + base_idx;
//		double* cs = c_vecs + base_idx;
//		double* rhs = rhs_vecs + base_idx;
//
//		// j is an index of an equation in the reduced system extended with
//		// x_0=0 and x_N=0
//#pragma omp parallel for if (equation_count > 16)
//		for (int j = 1; j < equation_count + 1; j += 2)
//		{
//			// index in coefficents arrays is less by 1 since they don't
//			// include x_0=0 and x_N=0 equations
//			int coeff_j = j - 1;
//			// mapping j to the index in the full extended system
//			int result_j = j * factor;
//			int result_j_minus1 = result_j - factor;
//			int result_j_plus1 = result_j + factor;
//
//			result[result_j] =
//				(rhs[coeff_j] - bs[coeff_j] * result[result_j_minus1] -
//					cs[coeff_j] * result[result_j_plus1]) / as[coeff_j];
//		}
//
//		equation_count = equation_count * 2 + 1;
//		factor /= 2;
//	}
//
//	// drop the padding
//	result.pop_back();
//	result.erase(result.begin());
//
//	delete[] a_vecs;
//	delete[] b_vecs;
//	delete[] c_vecs;
//	delete[] rhs_vecs;
//
//	return result;
////	int EquationCount = b.size();
////	int iReduce, Step, k; 
////	int index1, index2, offset;
////	double alpha, beta, gamma;
////	//Forward Reduction.
////
////	std::vector<double> LowerDiag;
////	std::vector<double> Diag;
////	std::vector<double> UpperDiag;
////
////	LowerDiag.resize(alpha);
////	Diag.resize(beta);
////	UpperDiag.resize(gamma);
////
////	for (Step = 1; Step < log2(EquationCount) ; Step++) {
////#pragma omp parallel for 
////		for (iReduce = pow(2, Step) - 1; iReduce < EquationCount; iReduce = iReduce + pow(2, Step)) {
////			f[iReduce] = d[iReduce - 1] + d[iReduce] + d[iReduce + 1];
////
////		}
////	}
//
//
//}

//#define OMP_SCHEDULE dynamic
//#define CHUNK 32

void CyclicReduction(std::vector<double> LDiag, std::vector<double> Diag, std::vector<double> UDiag, std::vector<double> oldResult, std::vector<double>& newResult)
{
	double Temp1;
	double Temp2;
	int iReducePlus, iReduceMinus, offset;

	LDiag.insert(LDiag.begin(), 0.0);
	UDiag.push_back(0.0);
	//Forward Reduction
	//#pragma omp parallel for
	for (int Step = 1; Step < log2(Diag.size() + 1); Step++)
	{
		//#pragma omp parallel for
		for (int iReduce = pow(2, Step) - 1; iReduce < Diag.size(); iReduce += pow(2, Step))
		{
			offset = pow(2.0, Step - 1);
			iReduceMinus = iReduce - offset;
			iReducePlus = iReduce + offset;

			if (iReduce == Diag.size() - 1)
			{
				Temp1 = LDiag[iReduce] / Diag[iReduceMinus];
				LDiag[iReduce] = -LDiag[iReduceMinus] * Temp1;
				UDiag[iReduce] = 0;
				Diag[iReduce] = Diag[iReduce] - (UDiag[iReduceMinus] * Temp1);

				oldResult[iReduce] = oldResult[iReduce] - oldResult[iReduceMinus] * Temp1;

				//newResult[iReduce] = (oldResult[iReduce] - (LDiag[iReduce] * newResult[iReduce - 1])) / Diag[iReduce];
			}
			else
			{
				Temp1 = LDiag[iReduce] / Diag[iReduceMinus];
				Temp2 = UDiag[iReduce] / Diag[iReducePlus];

				LDiag[iReduce] = -LDiag[iReduceMinus] * Temp1;
				UDiag[iReduce] = -UDiag[iReducePlus] * Temp2;
				Diag[iReduce] = Diag[iReduce] - (LDiag[iReducePlus] * Temp2) - (UDiag[iReduceMinus] * Temp1);
				oldResult[iReduce] = oldResult[iReduce] - (oldResult[iReduceMinus] * Temp1) - (oldResult[iReducePlus] * Temp2);

				//newResult[iReduce] = (oldResult[iReduce] - (LDiag[iReduce] * newResult[iReduce - 1]) - (UDiag[iReduce] * newResult[iReduce + 1])) / Diag[iReduce];
			}
			

		}

	}

	int size = newResult.size();
	int backSubPlus, backSubMinus;
	//Backward Substitution
	newResult[size - 1] = (oldResult[size - 1] - (LDiag[size - 1] * newResult[size - 2])) / Diag[size - 1];

	for (int Step = log2(size + 1) - 2; Step >= 0; Step--)
	{
		//#pragma omp parallel for
		for (int backSub = pow(2.0, Step + 1) - 1; backSub < size; backSub += pow(2.0, Step))
		{
			offset = pow(2.0, Step - 1);
			backSubMinus = backSub - offset;
			backSubPlus = backSub + offset;

			//newResult[backSub] = (oldResult[backSub] - (newResult[backSubMinus] * LDiag[backSub]) - (UDiag[backSub] * newResult[backSubPlus])) / Diag[backSub];
			if (backSub != backSubMinus) {
				if (backSubMinus - offset < 0) newResult[backSubMinus] = (oldResult[backSubMinus] - UDiag[backSubMinus] * newResult[backSubMinus + offset]) / Diag[backSubMinus];
				else newResult[backSubMinus] = (oldResult[backSubMinus] - LDiag[backSubMinus] * newResult[backSubMinus - offset] - UDiag[backSubMinus] * newResult[backSubMinus + offset]) / Diag[backSubMinus];
			}if (backSub != backSubPlus) {
				if (backSubPlus + offset >= size) newResult[backSubPlus] = (oldResult[backSubPlus] - LDiag[backSubPlus] * newResult[backSubPlus - offset]) / Diag[backSubPlus];
				else newResult[backSubPlus] = (oldResult[backSubPlus] - LDiag[backSubPlus] * newResult[backSubPlus - offset] - UDiag[backSubPlus] * newResult[backSubPlus + offset]) / Diag[backSubPlus];
			}
		}
	}
	//newResult[0] = oldResult[0] - (UDiag[0] * newResult[1]) / Diag[0];
	LDiag.erase(LDiag.begin());
	UDiag.pop_back();

}


// https://github.com/minavouronikou/cyclic_reduction/blob/master/cyclic%20reduction%20(L)/kernel.cu
// https://github.com/cudpp/cudpp/blob/master/src/cudpp/kernel/tridiagonal_kernel.cuh