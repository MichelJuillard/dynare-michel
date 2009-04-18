/*
 * Copyright (C) 2007-2009 Dynare Team
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
#include "mex.h"
#include "../matlab_versions_compatibility.h"

#ifdef NO_OPENMP
  #define USE_OMP 0
  #define DEBUG_OMP 0
  #ifdef NO_BLAS_H
    #if defined(__linux__) || defined(OCTAVE)
      #define dgemm dgemm_
    #endif
    extern "C" {
    lapack_int dgemm(char*, char*, lapack_int*, lapack_int*, lapack_int*, double*, double*, lapack_int*, double*, lapack_int*, double*, double*, lapack_int*);
    }
  #else /* NO_BLAS_H */
    #include "blas.h"
  #endif /* NO_BLAS_H */
#else
  #define USE_OMP 1
  #define DEBUG_OMP 0
  #include <omp.h>
#endif

#ifdef MWTYPES_NOT_DEFINED
typedef int mwIndex;
typedef int mwSize;
#endif

void full_A_times_kronecker_B_C(double *A, double *B, double *C, double *D,
				lapack_int mA, lapack_int nA, lapack_int mB, lapack_int nB, lapack_int mC, lapack_int nC)
{
  #if USE_OMP
    #pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
    for(long int colD = 0; colD<nB*nC; colD++)
      {
        #if DEBUG_OMP
          mexPrintf("%d thread number is %d (%d).\n",colD,omp_get_thread_num(),omp_get_num_threads());
        #endif
        unsigned int colB = colD/nC;
        unsigned int colC = colD%nC;
        for(int  colA = 0; colA<nA; colA++)
          {
            unsigned int rowB = colA/mC;
            unsigned int rowC = colA%mC;
            unsigned int idxA = colA*mA;
            unsigned int idxD = colD*mA;
            double BC = B[colB*mB+rowB]*C[colC*mC+rowC];
            for(int rowD = 0; rowD<mA; rowD++)
              {
                D[idxD+rowD] += A[idxA+rowD]*BC;
              }
          }
      }
  #else 
    const unsigned long shiftA = mA*mC ;
    const unsigned long shiftD = mA*nC ;
    unsigned long int kd = 0, ka = 0 ;
    char transpose[2] = "N";
    double one = 1.0 ;
    for(unsigned long int col=0; col<nB; col++)
      {
        ka = 0 ;
        for(unsigned long int row=0; row<mB; row++)
          {
            dgemm(transpose, transpose, &mA, &nC, &mC, &B[mB*col+row], &A[ka], &mA, &C[0], &mC, &one, &D[kd], &mA);
            ka += shiftA;
          }
        kd += shiftD;
      }
  #endif
}


void full_A_times_kronecker_B_B(double *A, double *B, double *D, lapack_int mA, lapack_int nA, lapack_int mB, lapack_int nB)
{
  #if USE_OMP
    #pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
    for(long int colD = 0; colD<nB*nB; colD++)
      {
        #if DEBUG_OMP
          mexPrintf("%d thread number is %d (%d).\n",colD,omp_get_thread_num(),omp_get_num_threads());
        #endif
        unsigned int j1B = colD/nB;
        unsigned int j2B = colD%nB;
        for(int  colA = 0; colA<nA; colA++)
          {
            unsigned int i1B = colA/mB;
            unsigned int i2B = colA%mB;
            unsigned int idxA = colA*mA;
            unsigned int idxD = colD*mA;
            double BB = B[j1B*mB+i1B]*B[j2B*mB+i2B];
            for(int rowD = 0; rowD<mA; rowD++)
              {
                D[idxD+rowD] += A[idxA+rowD]*BB;
              }
          }
      }
  #else
    const unsigned long int shiftA = mA*mB ;
    const unsigned long int shiftD = mA*nB ;
    unsigned long int kd = 0, ka = 0 ;
    char transpose[2] = "N";
    double one = 1.0;
    for(unsigned long int col=0; col<nB; col++)
      {
        ka = 0;
        for(unsigned long int row=0; row<mB; row++)
          {
            dgemm(transpose, transpose, &mA, &nB, &mB, &B[mB*col+row], &A[ka], &mA, &B[0], &mB, &one, &D[kd], &mA);
            ka += shiftA;
          }
        kd += shiftD;
      }
   #endif
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  // Check input and output:
  if ( (nrhs > 3) || (nrhs <2) )
    {
      mexErrMsgTxt("Two or Three input arguments required.");
    }
  if (nlhs>1)
    {
      mexErrMsgTxt("Too many output arguments.");
    }
  // Get & Check dimensions (columns and rows):
  mwSize mA, nA, mB, nB, mC, nC;
  mA = mxGetM(prhs[0]);
  nA = mxGetN(prhs[0]);
  mB = mxGetM(prhs[1]);
  nB = mxGetN(prhs[1]);
  if (nrhs == 3)// A*kron(B,C) is to be computed.
    {
      mC = mxGetM(prhs[2]);
      nC = mxGetN(prhs[2]);
      if (mB*mC != nA)
	{
	  mexErrMsgTxt("Input dimension error!");
	}
    }
  else// A*kron(B,B) is to be computed.
    {
      if (mB*mB != nA)
	{
	  mexErrMsgTxt("Input dimension error!");
	}
    }
  // Get input matrices:
  double *B, *C, *A;
  A = mxGetPr(prhs[0]);
  B = mxGetPr(prhs[1]);
  if (nrhs == 3)
    {
      C = mxGetPr(prhs[2]);
    }
  // Initialization of the ouput:
  double *D;
  if (nrhs == 3)
    {
      plhs[0] = mxCreateDoubleMatrix(mA,nB*nC,mxREAL);
    }
  else
    {
      plhs[0] = mxCreateDoubleMatrix(mA,nB*nB,mxREAL);
    }
  D = mxGetPr(plhs[0]);
  // Computational part:
  if (nrhs == 2)
    {
      full_A_times_kronecker_B_B(A, B, &D[0], (int) mA, (int) nA, (int) mB, (int) nB);
    }
  else
    {
      full_A_times_kronecker_B_C(A, B, C, &D[0], (int) mA, (int) nA, (int) mB, (int) nB, (int) mC, (int) nC);
    }
}
