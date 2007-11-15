// This mex file computes A*kron(B,C) or A*kron(B,B) without explicitely building kron(B,C) or kron(B,B), so that 
// one can consider large matrices B and/or C. 
//
// (linux)SYNTAX:
// mex AkronBC.cc /opt/matlab2007b/bin/glnx86/mkl.so
//
// stephane.adjemian@ens.fr [15-11-2007]
// Dynare Team, 2007.
#include "mex.h"
#include "blas.h"


void A_times_kronecker_B_B(double *A, double *B, double *D, int mA, int nA, int mB, int nB)
{
  const int shiftA = mA*mB ;
  const int shiftD = mA*nB ;
  char transpose[2] = "N";
  double one = 1.0;
  int kd = 0 ; 
  for(int col=0; col<nB; col++)
    {
      int ka = 0 ;
      for(int row=0; row<mB; row++)
	{
	  dgemm(transpose, transpose, &mA, &nB, &mB, &B[mB*col+row], &A[ka], &mA, &B[0], &mB, &one, &D[kd], &mA);
	  ka += shiftA;
	}
      kd += shiftD;
    }
}


void A_times_kronecker_B_C(double *A, double *B, double *C, double *D,
			   int mA, int nA, int mB, int nB, int mC, int nC)
{
  const int shiftA = mA*mC ;
  const int shiftD = mA*nC ;
  double one = 1.0 ;
  char transpose[2] = "N";
  int kd = 0 ; 
  for(int col=0; col<nB; col++)
    {
      int ka = 0 ;
      for(int row=0; row<mB; row++)
	{
	  dgemm(transpose, transpose, &mA, &nC, &mC, &B[mB*col+row], &A[ka], &mA, &C[0], &mC, &one, &D[kd], &mA);
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
  int mA, nA, mB, nB, mC, nC;
  mA = (int)mxGetM(prhs[0]);
  nA = (int)mxGetN(prhs[0]);
  mB = (int)mxGetM(prhs[1]);
  nB = (int)mxGetN(prhs[1]);
  if (nrhs == 3)// A*kron(B,C) is to be computed.
    {
      mC = (int)mxGetM(prhs[2]);
      nC = (int)mxGetN(prhs[2]);
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
  double *A, *B, *C;
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
      A_times_kronecker_B_B(A, B, &D[0], mA, nA, mB, nB);
    }
  else
    {
      A_times_kronecker_B_C(A, B, C, &D[0], mA, nA, mB, nB, mC, nC);
    }
}
