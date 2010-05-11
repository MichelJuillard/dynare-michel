
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
#include <dynmex.h>
#include <dynblas.h>
#include <dynlapack.h>
#else

#ifndef __BLAS_LAPACK__
#define __BLAS_LAPACK__

typedef int lapack_int;
typedef int blas_int;

#ifdef __cplusplus
extern "C"
{
#endif

/* Linux defines */
#define sscal     sscal_
#define saxpy     saxpy_
#define sgemm     sgemm_
#define sgetrf    sgetrf_
#define sgesdd    sgesdd_
#define sgesvd    sgesvd_
#define sgetrf    sgetrf_
#define sorgqr    sorgqr_
#define sgelqf    sgelqf_
#define sorglq    sorglq_
#define sgges     sgges_
#define stgsen    stgsen_
#define stgexc    stgexc_

#define dscal     dscal_     // Blas scalar times vector
#define daxpy     daxpy_     // Blas vector plus scalar times vector
#define dgemm     dgemm_     // Blas matrix multiplication
#define dgetrf    dgetrf_
#define dgesdd    dgesdd_    // SVD decomposition (divide and conquer)
#define dgesvd    dgesvd_    // SVD decomposition (QR)
#define dgetrf    dgetrf_    // LU decomposition
#define dgeqrf    dgeqrf_    // QR decomposition
#define dorgqr    dorgqr_    // Forms orthogonal matrix from Housholder matrices created by dgeqrf
#define dgelqf    dgelqf_    // LQ decompostion
#define dorglq    dorglq_    // Forms orthogonal matrix from Housholder matrices created by dgeqrf
#define dgges     dgges_     // Generalized Schur decomposition
#define dtgsen    dtgsen_    // Reorders generalized Schur decomposition
#define dtgexc    dtgexc_    // Reorders generalized Schur decomposition

#define dsyev     dsyev_
#define dgeev     dgeev_
#define dpotrf    dpotrf_
#define dpotri    dpotri_
#define dtrtri    dtrtri_
#define dgetri    dgetri_
#define dgeqp3    dgeqp3_
#define dormqr    dormqr_
#define dgesv     dgesv_
/*******************************************************************************/


/* cblas defines *
#define cblas_daxpy    daxpy
/*******************************************************************************/


/* Window defines */
void sscal(int *n, float *alpha, float *x, int *incx);
void saxpy(int *n, float *alpha, float *x, int *incx, float *y, int *incy);
void sgemm(char *transa, char *transb, int *m, int *n, int *k, float *alpha, float *a, int *lda, float *b, int *ldb, float *beta, float *c, int *ldc);
void sgetrf(int *M, int *N, float *A, int *LDA, int *IPIV, int *INFO);
void sgesdd(char *jobz, int *m, int *n, float *a, int *lda, float *s, float *u, int *ldu, float *vt, int *ldvt, float *work, int *lwork, int *iwork, int *info);
void sgesvd(char *jobu, char *jobvt, int *m, int *n, float *a, int *lda, float *s, float *u, int *ldu, float *vt, int *ldvt, float *work, int *lwork, int *info);
void sgeqrf(int *M, int *N, float *A, int *LDA, float *TAU, float *WORK, int *LWORK, int *INFO);
void sorgqr(int *M, int *N, int *K, float *A, int *LDA, float *TAU, float *WORK, int *LWORK, int *INFO);
void sgelqf(int *M, int *N, float *A, int *LDA, float *TAU, float *WORK, int *LWORK, int *INFO);
void sorglq(int *M, int *N, int *K, float *A, int *LDA, float *TAU, float *WORK, int *LWORK, int *INFO);
void sgges(char *jobvsl, char *jobvsr, char *sort, void *selctg, int *n, float *a, int *lda, float *b, int *ldb, int *sdim, float *alphar, float *alphai, float *beta, float *vsl, int *ldvsl, float *vsr, int *ldvsr, float *work, int *lwork, void *bwork, int *info);
void stgsen(int *ijob, void *wantq, void *wantz, void *select, int *n, float *a, int *lda, float *b, int *ldb, float *alphar, float *alphai, float *beta, float *q, int *ldq, float *z, int *ldz, int *m, float *pl, float *pr, float *dif, float *work, int *lwork, int *iwork, int *liwork, int *info);
void stgexc(void *wantq, void *wantz, int *n, float *a, int *lda, float *b, int *ldb, float *q, int *ldq, float *z, int *ldz, int *ifst, int *ilst, float *work, int *lwork, int *info);

void dscal(int*,double*,double*,int*);
void daxpy(int*,double*,double*,int*,double*,int*);
void dgemm(char*,char*,int*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);
void dgetrf(int*,int*,double*,int*,int*,int*);
void dgesdd(char*,int*,int*,double*,int*,double*,double*,int*,double*,int*,double*,int*,int*,int*);
void dgesvd(char*,char*,int*,int*,double*,int*,double*,double*,int*,double*,int*,double*,int*,int*);
void dgeqrf(int*,int*,double*,int*,double*,double*,int*,int*);
void dorgqr(int*,int*,int*,double*,int*,double*,double*,int*,int*);
void dgelqf(int*,int*,double*,int*,double*,double*,int*,int*);
void dorglq(int*,int*,int*,double*,int*,double*,double*,int*,int*);
void dgges(char*,char*,char*,void*,int*,double*,int*,double*,int*,int*,double*,double*,double*,double*,int*,double*,int*,double*,int*,void*,int*);
void dtgsen(int*,void*,void*,void*,int*,double*,int*,double*,int*,double*,double*,double*,double*,int*,double*,int*,int*,double*,double*,double*,double*,int*,int*,int*,int*);
void dtgexc(void*,void*,int*,double*,int*,double*,int*,double*,int*,double*,int*,int*,int*,double*,int*,int*);
void dsyev(char*,char*,int*,double*,int*,double*,double*,int*,int*);
void dgeev(char*,char*,int*,double*,int*,double*,double*,double*,int*,double*,int*,double*,int*,int*);
void dpotrf(char*,int*,double*,int*,int*);
void dpotri(char*,int*,double*,int*,int*);
void dgeqp3(int*,int*,double*,int*,int*,double*,double*,int*,int*);
void dtrtri(char*,char*,int*,double*,int*,int*);
void dgetri(int*,double*,int*,int*,double*,int*,int*);
void dormqr(char*,char*,int*,int*,int*,double*,int*,double*,double*,int*,double*,int*,int*);
void dgesv(int*,int*,double*,int*,int*,double*,int*,int*);
/*******************************************************************************/

#ifdef __cplusplus
}
#endif

#endif

#endif
