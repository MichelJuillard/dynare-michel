
#ifndef __BMATRIX__
#define __BMATRIX__

#include "prcsn.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* Unary Operators */
int bNegative(PRECISION *x, PRECISION *y, int m);
int bAbs(PRECISION *x, PRECISION *y, int m);
int bTranspose(PRECISION *x, PRECISION *y, int m, int n, int t);
int bTransposeInPlace(PRECISION *x, int m);

/* Addition */
int bAdd(PRECISION *x, PRECISION *y, PRECISION *z, int m);
int bSubtract(PRECISION *x, PRECISION *y, PRECISION *z, int m);
int bLinearUpdateScalar(PRECISION *x, PRECISION *y, PRECISION a, int m);
int bLinearCombination(PRECISION *x, PRECISION a, PRECISION *y, PRECISION b, PRECISION *z, int m);
int bMatrixAdd(PRECISION *x, PRECISION *y, PRECISION *z, int m, int n, int xt, int yt, int zt);
int bMatrixSubtract(PRECISION *x, PRECISION *y, PRECISION *z, int m, int n, int xt, int yt, int zt);
int bMatrixLinearCombination(PRECISION *x, PRECISION a, PRECISION *y, PRECISION b,PRECISION *z, int m, int n, int xt, int yt, int zt);

/* Multiplication */
int bMultiply(PRECISION *x, PRECISION *y, PRECISION s, int m);
int bMatrixMultiply(PRECISION *x, PRECISION *y, PRECISION *z, int m, int n, int p, int xt, int yt, int zt);

/* LU Decomposition */
int bLU(int *p, PRECISION *x, int m, int n, int xt);
int bSolveTriangular(PRECISION *x, PRECISION *b, int m, int n,int u,  int xt, int bt);
int bSolveUnitTriangular(PRECISION *x, PRECISION *b, int m, int n, int u, int xt, int bt);

/* QR Decompositions */
int bQR(PRECISION *Q, PRECISION *R, PRECISION *X, int m, int n, int q, int qt, int rt, int xt);

/* Singular Value Decomposition */
int bSVD(PRECISION *U, PRECISION *d, PRECISION *V, PRECISION *A, int m, int n, int ut, int vt, int at);
int bSVD_new(PRECISION *U, PRECISION *d, PRECISION *V, PRECISION *A, int m, int n, int ut, int vt, int at, int compact);

/* Generalize Schur Decomposition */
int bQZ_real(PRECISION *Q, PRECISION *Z, PRECISION *S, PRECISION *T, PRECISION *A, PRECISION *B, int n, int qt, int zt, int st, int tt, int at, int bt,
         PRECISION *alpha_r, PRECISION *alpha_i, PRECISION *beta);

/* Cholesky Decompositions */
int bCholesky(PRECISION *X, int m, int u, int t);

/* Permutation Routines */
int bPermutationMultiply(int *p, PRECISION *y, int m, int n, int q, int pt, int yt);
int bPermutation(PRECISION *x, int *p, int m, int q, int t);

/* Tensor Calculus */
int bMatrixTensor(PRECISION *x, PRECISION *y, PRECISION *z, int m, int n, int r, int s, int xt, int yt, int zt);
int bVectorTensor(PRECISION *x, PRECISION *y, PRECISION *z, int m, int n);

//int bQRPivot_R(PRECISION *R, int *p, int m, int n);
//int bQRPivot_QR(PRECISION *Q, PRECISION *R, int *p, int m, int n);

#ifdef __cplusplus
}
#endif

#endif


