/* kalman_steady_state.cc
**
** Copyright (C) 2009-2013 Dynare Team.
**
** This file is part of Dynare.
**
** Dynare is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** Dynare is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
**
** This mex file calls fortran routines from the Slicot library.
*/

/*
  ++    INPUTS
  ++    ======
  ++
  ++
  ++      [0]  T       (double)   n-by-n transition matrix.
  ++
  ++      [1]  QQ      (double)   n-by-n matrix (=R*Q*R', where Q is the covariance matrix of the structural innovations).
  ++
  ++      [2]  Z       (double)   n-by-p selection matrix.
  ++
  ++      [3]  H       (double)   p-by-p covariance matrix of the measurement errors.
  ++
  ++
  ++
  ++
  ++    OUTPUTS
  ++    =======
  ++
  ++
  ++      [0]  P       (double)   n-by-n covariance matrix of the state vector.
  ++
  ++
  ++    NOTES
  ++    =====
  ++
  ++    [1] T = transpose(dynare transition matrix) and Z = transpose(dynare selection matrix).
*/

#include <string.h>
#include <stdlib.h>
#include <dynmex.h>
#include <dynlapack.h>

#if !defined(MATLAB_MEX_FILE) || !defined(_WIN32)
# define sb02od sb02od_
#endif

extern "C"
{
  int sb02od(char *, char *, char *, char *, char *, char *, mwSize *, mwSize *, mwSize *, double *, mwSize *, double *, mwSize *, double *, mwSize *, double *, mwSize *, double *, mwSize *, double *, double *, mwSize *, double *, double *, double *, double *, mwSize *, double *, mwSize *, double *, mwSize *, double *, lapack_int *, double *, mwSize *, lapack_int *, lapack_int *);
}

template <typename T>
T
max(T x, T y)
{
  return x < y ? y : x;
}

template <typename T>
T
max(T x, T y, T z)
{
  return max(x, max(y, z));
}

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Check the number of arguments and set some flags.
  int measurement_error_flag = 1;
  if (nrhs < 3 || 4 < nrhs)
    DYN_MEX_FUNC_ERR_MSG_TXT("kalman_steady_state accepts either 3 or 4 input arguments!");

  if (nlhs < 1 || 2 < nlhs)
    DYN_MEX_FUNC_ERR_MSG_TXT("kalman_steady_state requires at least 1, but no more than 2, output arguments!");

  if (nrhs == 3)
    measurement_error_flag = 0;

  // Check the type of the input arguments and get the size of the matrices.
  mwSize n = mxGetM(prhs[0]);
  if ((size_t) n != mxGetN(prhs[0]))
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("kalman_steady_state: The first input argument (T) must be a square matrix!");
    }
  if ((mxIsNumeric(prhs[0]) == 0) || (mxIsComplex(prhs[0]) == 1))
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("kalman_steady_state: The first input argument (T) must be a real matrix!");
    }
  mwSize q = mxGetM(prhs[1]);
  if ((size_t) q != mxGetN(prhs[1]))
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("kalman_steady_state: The second input argument (QQ) must be a square matrix!");
    }
  if ((mxIsNumeric(prhs[1]) == 0) || (mxIsComplex(prhs[1]) == 1))
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("kalman_steady_state: The second input argument (QQ) must be a real matrix!");
    }
  if (q != n)
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("kalman_steady_state: The size of the second input argument (QQ) must match the size of the first argument (T)!");
    }
  mwSize p = mxGetN(prhs[2]);
  if (mxGetM(prhs[2]) != (size_t) n)
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("kalman_steady_state: The number of rows of the third argument (Z) must match the number of rows of the first argument (T)!");
    }
  if ((mxIsNumeric(prhs[2]) == 0) || (mxIsComplex(prhs[2]) == 1))
    {
      DYN_MEX_FUNC_ERR_MSG_TXT("kalman_steady_state: The third input argument (Z) must be a real matrix!");
    }
  if (measurement_error_flag)
    {
      if (mxGetM(prhs[3]) != mxGetN(prhs[3]))
        {
          DYN_MEX_FUNC_ERR_MSG_TXT("kalman_steady_state: The fourth input argument (H) must be a square matrix!");
        }
      if (mxGetM(prhs[3]) != (size_t) p)
        {
          DYN_MEX_FUNC_ERR_MSG_TXT("kalman_steady_state: The number of rows of the fourth input argument (H) must match the number of rows of the third input argument!");
        }
      if ((mxIsNumeric(prhs[3]) == 0) || (mxIsComplex(prhs[3]) == 1))
        {
          DYN_MEX_FUNC_ERR_MSG_TXT("kalman_steady_state: The fifth input argument (H) must be a real matrix!");
        }
    }
  // Get input matrices.
  double *T, *QQ, *Z, *H, *L; // Remark. L will not be used.
  T = (double *) mxCalloc(n*n, sizeof(double));
  memcpy(T, mxGetPr(prhs[0]), n*n*sizeof(double));
  QQ = (double *) mxCalloc(n*n, sizeof(double));
  memcpy(QQ, mxGetPr(prhs[1]), n*n*sizeof(double));
  Z = (double *) mxCalloc(n*p, sizeof(double));
  memcpy(Z, mxGetPr(prhs[2]), n*p*sizeof(double));
  H = (double *) mxCalloc(p*p, sizeof(double));
  if (measurement_error_flag)
    {
      memcpy(H, mxGetPr(prhs[3]), p*p*sizeof(double));
    }
  L = (double *) mxCalloc(n*p, sizeof(double));
  char *DICO, *JOBB, *FACT, *UPLO, *JOBL, *SORT;
  DICO = (char *) mxCalloc(2, sizeof(char));
  memcpy(DICO, "D", 2*sizeof(char)); // We want to solve a discrete Riccati equation.
  JOBB = (char *) mxCalloc(2, sizeof(char));
  memcpy(JOBB, "B", 2*sizeof(char)); // Matrices Z and H are given.
  FACT = (char *) mxCalloc(2, sizeof(char));
  memcpy(FACT, "N", 2*sizeof(char)); // Given matrices H and QQ are not factored.
  UPLO = (char *) mxCalloc(2, sizeof(char));
  memcpy(UPLO, "U", 2*sizeof(char)); // Upper triangle of matrix H is stored.
  JOBL = (char *) mxCalloc(2, sizeof(char));
  memcpy(JOBL, "Z", 2*sizeof(char)); // L matrix is zero.
  SORT = (char *) mxCalloc(2, sizeof(char));
  memcpy(SORT, "S", 2*sizeof(char)); // Stable eigenvalues come first.
  mwSize nn = 2*n;
  mwSize LDA = max((mwSize) 1, n);
  mwSize LDQ = LDA;
  mwSize LDU = max((mwSize) 1, nn);
  mwSize LDS = max((mwSize) 1, nn+p);
  mwSize LIWORK = max((mwSize) 1, p, nn);
  mwSize LDR = max((mwSize) 1, p);
  mwSize LDB = LDA;
  mwSize LDL = LDA;
  mwSize LDT = LDS;
  mwSize LDX = LDA;
  mwSize LDWORK = max((mwSize) 7*((mwSize) 2*n + (mwSize) 1) + (mwSize) 16, (mwSize) 16*n);
  LDWORK = max(LDWORK, (mwSize) 2*n + p, (mwSize) 3*p);
  double tolerance = .0000000000000001;
  lapack_int INFO;
  // Outputs of subroutine sb02OD
  double rcond;
  double *WR, *WI, *BETA, *S, *TT, *UU;
  WR = (double *) mxCalloc(nn, sizeof(double));
  WI = (double *) mxCalloc(nn, sizeof(double));
  BETA = (double *) mxCalloc(nn, sizeof(double));
  S = (double *) mxCalloc(LDS*(nn+p), sizeof(double));
  TT = (double *) mxCalloc(LDT*nn, sizeof(double));
  UU = (double *) mxCalloc(LDU*nn, sizeof(double));
  // Working arrays
  lapack_int *IWORK;
  IWORK = (lapack_int *) mxCalloc(LIWORK, sizeof(lapack_int));
  double *DWORK;
  DWORK = (double *) mxCalloc(LDWORK, sizeof(double));
  lapack_int *BWORK;
  BWORK = (lapack_int *) mxCalloc(nn, sizeof(lapack_int));
  // Initialize the output of the mex file
  double *P;
  plhs[1] = mxCreateDoubleMatrix(n, n, mxREAL);
  P = mxGetPr(plhs[1]);
  // Call the slicot routine
  sb02od(DICO, JOBB, FACT, UPLO, JOBL, SORT, &n, &p, &p, &T[0], &LDA, &Z[0], &LDB, &QQ[0], &LDQ, &H[0], &LDR, &L[0], &LDL, &rcond, &P[0], &LDX, &WR[0], &WI[0], &BETA[0], &S[0], &LDS, &TT[0], &LDT, &UU[0], &LDU, &tolerance, &IWORK[0], &DWORK[0], &LDWORK, &BWORK[0], &INFO);
  mxFree(T);
  mxFree(QQ);
  mxFree(Z);
  mxFree(H);
  mxFree(L);
  mxFree(DICO);
  mxFree(JOBB);
  mxFree(FACT);
  mxFree(UPLO);
  mxFree(JOBL);
  mxFree(SORT);
  mxFree(WR);
  mxFree(WI);
  mxFree(BETA);
  mxFree(S);
  mxFree(TT);
  mxFree(UU);
  mxFree(IWORK);
  mxFree(DWORK);
  mxFree(BWORK);
  if (INFO != 0)
    {
      switch (INFO)
        {
        case 1:
          {
            DYN_MEX_FUNC_ERR_MSG_TXT("The computed extended matrix pencil is singular, possibly due to rounding errors.\n");
            break;
          }
        case 2:
          {
            DYN_MEX_FUNC_ERR_MSG_TXT("The QZ (or QR) algorithm failed!\n");
            break;
          }
        case 3:
          {
            DYN_MEX_FUNC_ERR_MSG_TXT("The reordering of the (generalized) eigenvalues failed!\n");
            break;
          }
        case 4:
          {
            DYN_MEX_FUNC_ERR_MSG_TXT("After reordering, roundoff changed values of some complex eigenvalues so that leading eigenvalues\n in the (generalized) Schur form no longer satisfy the stability condition; this could also be caused due to scaling.");
            break;
          }
        case 5:
          {
            DYN_MEX_FUNC_ERR_MSG_TXT("The computed dimension of the solution does not equal n!\n");
            break;
          }
        case 6:
          {
            DYN_MEX_FUNC_ERR_MSG_TXT("A singular matrix was encountered during the computation of the solution matrix P!\n");
            break;
          }
        default:
          {
            DYN_MEX_FUNC_ERR_MSG_TXT("Unknown problem!\n");
          }
        }
    }
  plhs[0] = mxCreateDoubleScalar(0);
}
