/* This mex file computes A*kron(B,C) or A*kron(B,B) without explicitely building kron(B,C) or kron(B,B), so that 
** one can consider large matrices B and/or C. 
**
** (linux)SYNTAX:
** mex AkronBC.cc /opt/matlab2007b/bin/glnx86/mkl.so
**
** stephane.adjemian@ens.fr [15-11-2007]
** Dynare Team, 2007.
*/
#include <string.h>
#include "mex.h"
#include "matrix.h"

#ifdef MWTYPES_NOT_DEFINED
typedef int mwIndex;
typedef int mwSize;
#endif

#ifdef NO_BLAS_H
#  if defined(__linux__)
#    define DGEMM dgemm_
#  else
#    define DGEMM dgemm
#  endif
extern "C"{
  int DGEMM(char*, char*, mwSize*, mwSize*, mwSize*, double*, double*, mwSize*, double*, mwSize*, double*, double*, mwSize*);
}
#else
#  include "blas.h"
#  define DGEMM dgemm
#endif

void full_A_times_kronecker_B_C(double *A, double *B, double *C, double *D,
			   mwSize mA, mwSize nA, mwSize mB, mwSize nB, mwSize mC, mwSize nC)
{
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
	  DGEMM(transpose, transpose, &mA, &nC, &mC, &B[mB*col+row], &A[ka], &mA, &C[0], &mC, &one, &D[kd], &mA);
	  ka += shiftA;
	}
      kd += shiftD;
    }
}


void full_A_times_kronecker_B_B(double *A, double *B, double *D, mwSize mA, mwSize nA, mwSize mB, mwSize nB)
{
  const unsigned long int shiftA = mA*mB ;
  const unsigned long int shiftD = mA*nB ;
  unsigned long int kd = 0, ka = 0 ;
  char transpose[2] = "N";
  double one = 1.0;
  for(unsigned long int col=0; col<nB; col++)
    {
      ka = 0 ;
      for(unsigned long int row=0; row<mB; row++)
	{
	  DGEMM(transpose, transpose, &mA, &nB, &mB, &B[mB*col+row], &A[ka], &mA, &B[0], &mB, &one, &D[kd], &mA);
	  ka += shiftA;
	}
      kd += shiftD;
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
      full_A_times_kronecker_B_B(A, B, &D[0], mA, nA, mB, nB);
    }
  else
    {
      full_A_times_kronecker_B_C(A, B, C, &D[0], mA, nA, mB, nB, mC, nC);
    }
}
