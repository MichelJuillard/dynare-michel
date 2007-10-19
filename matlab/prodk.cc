#include "mex.h"
#include "blas.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  int inc = 1;
  char tA[2] = {'N','\0'};
  double beta = 1.0;


  // Check input and output:
  if (nrhs != 3)
    {
      mexErrMsgTxt("Three input arguments required.");
    }
  else if (nlhs>1)
    {
      mexErrMsgTxt("Too many output arguments.");
    }

  // Get & Check dimensions:
  int mA, nA, mb, nb, mc, nc;
  mA = (int)mxGetM(prhs[0]);
  nA = (int)mxGetN(prhs[0]);
  mb = (int)mxGetM(prhs[1]);
  nb = (int)mxGetN(prhs[1]);
  if (nb != 1)
    {
      mexErrMsgTxt("Second argument must be a column vector.");
    }
  mc = (int)mxGetM(prhs[2]);
  nc = (int)mxGetN(prhs[1]);
  if (nc != 1)
    {
      mexErrMsgTxt("Third argument must be a column vector.");
    }
  if (mb*mc != nA)
    {
      mexErrMsgTxt("Input dimension error.");
    }
  
  // Get input matrices:
  double *A, *b, *c ;
  A = mxGetPr(prhs[0]);
  b = mxGetPr(prhs[1]);
  c = mxGetPr(prhs[2]);

  // Initialization of the ouput:
  double *d;
  plhs[0] = mxCreateDoubleMatrix(mA,1,mxREAL);
  d = mxGetPr(plhs[0]);  

  // Computational part:
  int k = 0;
  for (int i=0; i<mb ; i++)
    {
      dgemv(tA, &mA, &mc, &b[i], &A[k], &mA, &c[0], &inc, &beta, &d[0], &inc);
      k += mc*mA;
    }
}
