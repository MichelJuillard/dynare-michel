/* This mex file computes A*kron(B,C) or A*kron(B,B) without explicitly building kron(B,C) or kron(B,B), so that
** one can consider large matrices A, B and/or C, and assuming that A is a the hessian of a dsge model 
** (dynare format). This mex file should not be used outside dr1.m.
**
** (linux)SYNTAX:
** mex sparse_hessian_times_B_kronecker_B.cc
**
** stephane.adjemian@ens.fr
** Dynare Team, 2007.
*/
#include <string.h>
#include "matrix.h"
#include "mex.h"

#ifdef MWTYPES_NOT_DEFINED
typedef unsigned int mwIndex;
typedef unsigned int mwSize;
#endif

void sparse_hessian_times_B_kronecker_B(mwIndex *isparseA, mwIndex *jsparseA, double *vsparseA, 
					double *B, double *D, mwSize mA, mwSize nA, mwSize mB, mwSize nB)
{
  /* 
  **   Loop over the columns of kron(B,B) (or of the result matrix D).
  **   This loop is splitted into two nested loops because we use the
  **   symmetric pattern of the hessian matrix.     
  */
  unsigned long int jj, ii, iv;
  unsigned int i1B, i2B, j1B, j2B, k1, k2, kk, k;
  unsigned int nz_in_column_ii_of_A;
  double bb;
  for(j1B=0; j1B<nB; j1B++)
    {
      for(j2B=j1B; j2B<nB; j2B++)
	{
	  jj = j1B*nB+j2B;// column of kron(B,B) index.
	  nz_in_column_ii_of_A = 0;
	  k1 = k2 = iv = 0;
	  /*
	  ** Loop over the rows of kron(B,B) (column jj).
	  */
	  for(ii=0; ii<nA; ii++)
	    {
	      k1 = jsparseA[ii];
	      k2 = jsparseA[ii+1];
	      if (k1 < k2)// otherwise column ii of A does not have non zero elements (and there is nothing to compute).
		{
		  ++nz_in_column_ii_of_A;
		  i1B = (ii/mB);
		  i2B = (ii%mB);
		  bb = B[j1B*mB+i1B]*B[j2B*mB+i2B];
		  /*
		  ** Loop over the non zero entries of A(:,ii).
		  */
		  for(k=k1; k<k2; k++)
		    {
		      kk = isparseA[k];
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
  unsigned long int jj, ii, iv;
  unsigned int iB, iC, jB, jC, k1, k2, kk, k;
  unsigned int nz_in_column_ii_of_A;
  double cb;
  for(jj=0; jj<nB*nC; jj++)// column of kron(B,B) index.
    {
      jB = jj/nC;
      jC = jj%nC;
      iv = k1 = k2 = 0;
      nz_in_column_ii_of_A = 0;
      /*
      ** Loop over the rows of kron(B,B) (column jj).
      */
      for(ii=0; ii<nA; ii++)
	{
	  k1 = jsparseA[ii];
	  k2 = jsparseA[ii+1];
	  if (k1 < k2)// otherwise column ii of A does not have non zero elements (and there is nothing to compute).
	    {
	      ++nz_in_column_ii_of_A;
	      iC = (ii%mB);
	      iB = (ii/mB);
	      cb = C[jC*mC+iC]*B[jB*mB+iB];
	      /*
	      ** Loop over the non zero entries of A(:,ii).
	      */
	      for(k=k1; k<k2; k++)
		{
		  kk = isparseA[k];
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
