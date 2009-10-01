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
 * This mex file computes A*kron(B,C) or A*kron(B,B) without explicitly building kron(B,C) or kron(B,B), so that
 * one can consider large matrices A, B and/or C, and assuming that A is a the hessian of a dsge model 
 * (dynare format). This mex file should not be used outside dr1.m.
 */

#include <string.h>

#include <dynmex.h>

#ifdef USE_OMP
# include <omp.h>
#endif

void sparse_hessian_times_B_kronecker_B(mwIndex *isparseA, mwIndex *jsparseA, double *vsparseA, 
					double *B, double *D, mwSize mA, mwSize nA, mwSize mB, mwSize nB)
{
  /* 
  **   Loop over the columns of kron(B,B) (or of the result matrix D).
  **   This loop is splitted into two nested loops because we use the
  **   symmetric pattern of the hessian matrix.     
  */
  #if USE_OMP
    #pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  #endif
  for(int j1B=0; j1B<nB; j1B++)
    {
      #if DEBUG_OMP
        mexPrintf("%d thread number is %d (%d).\n",j1B,omp_get_thread_num(),omp_get_num_threads());
      #endif
      for(unsigned int j2B=j1B; j2B<nB; j2B++)
	{
	  unsigned long int jj = j1B*nB+j2B;// column of kron(B,B) index.
	  unsigned long int iv =0;
          unsigned int nz_in_column_ii_of_A = 0;
	  unsigned int k1 = 0;
          unsigned int k2 = 0;
	  /*
	  ** Loop over the rows of kron(B,B) (column jj).
	  */
	  for(unsigned long int ii=0; ii<nA; ii++)
	    {
	      k1 = jsparseA[ii];
	      k2 = jsparseA[ii+1];
	      if (k1 < k2)// otherwise column ii of A does not have non zero elements (and there is nothing to compute).
		{
		  ++nz_in_column_ii_of_A;
		  unsigned int i1B = (ii/mB);
		  unsigned int i2B = (ii%mB);
		  double bb  = B[j1B*mB+i1B]*B[j2B*mB+i2B];
		  /*
		  ** Loop over the non zero entries of A(:,ii).
		  */
		  for(unsigned int k=k1; k<k2; k++)
		    {
		      unsigned int kk = isparseA[k];
		      D[jj*mA+kk] = D[jj*mA+kk] + bb*vsparseA[iv];
		      iv++;
		    }
		}
	    }
	  if (nz_in_column_ii_of_A>0)
	    {
	      memcpy(&D[(j2B*nB+j1B)*mA],&D[jj*mA],mA*sizeof(double));
	    }
	}
    }
}

void sparse_hessian_times_B_kronecker_C(mwIndex *isparseA, mwIndex *jsparseA, double *vsparseA, 
					double *B, double *C, double *D, 
					mwSize mA, mwSize nA, mwSize mB, mwSize nB, mwSize mC, mwSize nC)
{
  /* 
  **   Loop over the columns of kron(B,B) (or of the result matrix D).
  */
  #if USE_OMP
    #pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  #endif
  for(long int jj=0; jj<nB*nC; jj++)// column of kron(B,C) index.
    {
      // Uncomment the following line to check if all processors are used.
      #if DEBUG_OMP
        mexPrintf("%d thread number is %d (%d).\n",jj,omp_get_thread_num(),omp_get_num_threads());
      #endif
      unsigned int jB = jj/nC;
      unsigned int jC = jj%nC;
      unsigned int k1 = 0;
      unsigned int k2 = 0;
      unsigned long int iv = 0;
      unsigned int nz_in_column_ii_of_A = 0;
      /*
      ** Loop over the rows of kron(B,C) (column jj).
      */
      for(unsigned long int ii=0; ii<nA; ii++)
	{
	  k1 = jsparseA[ii];
	  k2 = jsparseA[ii+1];
	  if (k1 < k2)// otherwise column ii of A does not have non zero elements (and there is nothing to compute).
	    {
	      ++nz_in_column_ii_of_A;
	      unsigned int iC = (ii%mB);
	      unsigned int iB = (ii/mB);
	      double cb = C[jC*mC+iC]*B[jB*mB+iB];
	      /*
	      ** Loop over the non zero entries of A(:,ii).
	      */
	      for(unsigned int k=k1; k<k2; k++)
		{
		  unsigned int kk = isparseA[k];
		  D[jj*mA+kk] = D[jj*mA+kk] + cb*vsparseA[iv];
		  iv++;
		}
	    }
	}
    }
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
  if (!mxIsSparse(prhs[0]))
    {
      mexErrMsgTxt("First input must be a sparse (dynare) hessian matrix.");
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
  double *B, *C;
  B = mxGetPr(prhs[1]);
  if (nrhs == 3)
    {
      C = mxGetPr(prhs[2]);
    }
  // Sparse (dynare) hessian matrix.
  mwIndex *isparseA = (mwIndex*)mxGetIr(prhs[0]);
  mwIndex *jsparseA = (mwIndex*)mxGetJc(prhs[0]);
  double  *vsparseA = mxGetPr(prhs[0]);
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
      sparse_hessian_times_B_kronecker_B(isparseA, jsparseA, vsparseA, B, D, mA, nA, mB, nB);
    }
  else
    {
      sparse_hessian_times_B_kronecker_C(isparseA, jsparseA, vsparseA, B, C, D, mA, nA, mB, nB, mC, nC);
    }
}
