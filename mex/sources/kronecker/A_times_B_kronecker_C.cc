/*
 * Copyright (C) 2007-2011 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * This mex file computes A*kron(B,C) or A*kron(B,B) without explicitely building kron(B,C) or kron(B,B), so that
 * one can consider large matrices B and/or C.
 */

#include <string.h>

#include <dynmex.h>
#include <dynblas.h>

#ifdef USE_OMP
# include <omp.h>
#endif

#define DEBUG_OMP 0

void
full_A_times_kronecker_B_C(double *A, double *B, double *C, double *D,
                           blas_int mA, blas_int nA, blas_int mB, blas_int nB, blas_int mC, blas_int nC, int number_of_threads)
{
#if USE_OMP
# pragma omp parallel for num_threads(number_of_threads)
  for (blas_int colD = 0; colD < nB*nC; colD++)
    {
# if DEBUG_OMP
      mexPrintf("%d thread number is %d (%d).\n", colD, omp_get_thread_num(), omp_get_num_threads());
# endif
      blas_int colB = colD/nC;
      blas_int colC = colD%nC;
      for (blas_int colA = 0; colA < nA; colA++)
        {
          blas_int rowB = colA/mC;
          blas_int rowC = colA%mC;
          blas_int idxA = colA*mA;
          blas_int idxD = colD*mA;
          double BC = B[colB*mB+rowB]*C[colC*mC+rowC];
          for (blas_int rowD = 0; rowD < mA; rowD++)
            {
              D[idxD+rowD] += A[idxA+rowD]*BC;
            }
        }
    }
#else
  const blas_int shiftA = mA*mC;
  const blas_int shiftD = mA*nC;
  blas_int kd = 0, ka = 0;
  char transpose[2] = "N";
  double one = 1.0;
  for (blas_int col = 0; col < nB; col++)
    {
      ka = 0;
      for (blas_int row = 0; row < mB; row++)
        {
          dgemm(transpose, transpose, &mA, &nC, &mC, &B[mB*col+row], &A[ka], &mA, &C[0], &mC, &one, &D[kd], &mA);
          ka += shiftA;
        }
      kd += shiftD;
    }
#endif
}

void
full_A_times_kronecker_B_B(double *A, double *B, double *D, blas_int mA, blas_int nA, blas_int mB, blas_int nB, int number_of_threads)
{
#if USE_OMP
# pragma omp parallel for num_threads(number_of_threads)
  for (blas_int colD = 0; colD < nB*nB; colD++)
    {
# if DEBUG_OMP
      mexPrintf("%d thread number is %d (%d).\n", colD, omp_get_thread_num(), omp_get_num_threads());
# endif
      blas_int j1B = colD/nB;
      blas_int j2B = colD%nB;
      for (blas_int colA = 0; colA < nA; colA++)
        {
          blas_int i1B = colA/mB;
          blas_int i2B = colA%mB;
          blas_int idxA = colA*mA;
          blas_int idxD = colD*mA;
          double BB = B[j1B*mB+i1B]*B[j2B*mB+i2B];
          for (blas_int rowD = 0; rowD < mA; rowD++)
            {
              D[idxD+rowD] += A[idxA+rowD]*BB;
            }
        }
    }
#else
  const blas_int shiftA = mA*mB;
  const blas_int shiftD = mA*nB;
  blas_int kd = 0, ka = 0;
  char transpose[2] = "N";
  double one = 1.0;
  for (blas_int col = 0; col < nB; col++)
    {
      ka = 0;
      for (blas_int row = 0; row < mB; row++)
        {
          dgemm(transpose, transpose, &mA, &nB, &mB, &B[mB*col+row], &A[ka], &mA, &B[0], &mB, &one, &D[kd], &mA);
          ka += shiftA;
        }
      kd += shiftD;
    }
#endif
}

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Check input and output:
  if (nrhs > 4 || nrhs < 3 || nlhs != 2)
    DYN_MEX_FUNC_ERR_MSG_TXT("A_times_B_kronecker_C takes 3 or 4 input arguments and provides exactly 2 output arguments.");

  // Get & Check dimensions (columns and rows):
  mwSize mA, nA, mB, nB, mC, nC;
  mA = mxGetM(prhs[0]);
  nA = mxGetN(prhs[0]);
  mB = mxGetM(prhs[1]);
  nB = mxGetN(prhs[1]);
  if (nrhs == 4) // A*kron(B,C) is to be computed.
    {
      mC = mxGetM(prhs[2]);
      nC = mxGetN(prhs[2]);
      if (mB*mC != nA)
        DYN_MEX_FUNC_ERR_MSG_TXT("Input dimension error!");
    }
  else // A*kron(B,B) is to be computed.
    {
      if (mB*mB != nA)
        DYN_MEX_FUNC_ERR_MSG_TXT("Input dimension error!");
    }
  // Get input matrices:
  double *B, *C, *A;
  int numthreads;
  A = mxGetPr(prhs[0]);
  B = mxGetPr(prhs[1]);
  if (nrhs == 4)
    {
      C = mxGetPr(prhs[2]);
      numthreads = (int) mxGetScalar(prhs[3]);
    }
  else
    numthreads = (int) mxGetScalar(prhs[2]);

  // Initialization of the ouput:
  double *D;
  if (nrhs == 4)
    {
      plhs[1] = mxCreateDoubleMatrix(mA, nB*nC, mxREAL);
    }
  else
    {
      plhs[1] = mxCreateDoubleMatrix(mA, nB*nB, mxREAL);
    }
  D = mxGetPr(plhs[1]);
  // Computational part:
  if (nrhs == 3)
    {
      full_A_times_kronecker_B_B(A, B, &D[0], mA, nA, mB, nB, numthreads);
    }
  else
    {
      full_A_times_kronecker_B_C(A, B, C, &D[0], mA, nA, mB, nB, mC, nC, numthreads);
    }
  plhs[0] = mxCreateDoubleScalar(0);
}
