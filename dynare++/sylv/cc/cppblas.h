/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/cppblas.h,v 1.2 2004/11/24 20:42:52 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef CPPBLAS_H
#define CPPBLAS_H

#if defined(MATLAB) && !defined(__linux__)
#define BLAS_dgemm dgemm
#define BLAS_dgemv dgemv
#define BLAS_dtrsv dtrsv
#define BLAS_dtrmv dtrmv
#define BLAS_daxpy daxpy
#define BLAS_dcopy dcopy
#define BLAS_zaxpy zaxpy
#define BLAS_dscal dscal
#define BLAS_dtrsm dtrsm
#define BLAS_ddot  ddot
#else /* defined(MATLAB) && !defined(__linux__) */
#define BLAS_dgemm dgemm_
#define BLAS_dgemv dgemv_
#define BLAS_dtrsv dtrsv_
#define BLAS_dtrmv dtrmv_
#define BLAS_daxpy daxpy_
#define BLAS_dcopy dcopy_
#define BLAS_zaxpy zaxpy_
#define BLAS_dscal dscal_
#define BLAS_dtrsm dtrsm_
#define BLAS_ddot  ddot_
#endif /* defined(MATLAB) && !defined(__linux__) */

#if defined(MATLAB)
#include ../../../mex/sources/matlab_versions_compatibility.h
#endif

#define BLCHAR const char*
#define CONST_BLINT const blas_int*
#define CONST_BLDOU const double*
#define BLDOU double*

extern "C" {
	void BLAS_dgemm(BLCHAR transa, BLCHAR transb, CONST_BLINT m, CONST_BLINT n,
					CONST_BLINT k, CONST_BLDOU alpha, CONST_BLDOU a, CONST_BLINT lda,
					CONST_BLDOU b, CONST_BLINT ldb, CONST_BLDOU beta,
					BLDOU c, CONST_BLINT ldc);
	void BLAS_dgemv(BLCHAR trans, CONST_BLINT m, CONST_BLINT n, CONST_BLDOU alpha,
					CONST_BLDOU a, CONST_BLINT lda, CONST_BLDOU x, CONST_BLINT incx,
					CONST_BLDOU beta, BLDOU y, CONST_BLINT incy);
	void BLAS_dtrsv(BLCHAR uplo, BLCHAR trans, BLCHAR diag, CONST_BLINT n,
					CONST_BLDOU a, CONST_BLINT lda, BLDOU x, CONST_BLINT incx);
	void BLAS_dtrmv(BLCHAR uplo, BLCHAR trans, BLCHAR diag, CONST_BLINT n,
					CONST_BLDOU a, CONST_BLINT lda, BLDOU x, CONST_BLINT incx);
	void BLAS_daxpy(CONST_BLINT n, CONST_BLDOU a, CONST_BLDOU x, CONST_BLINT incx,
					BLDOU y, CONST_BLINT incy);
	void BLAS_dcopy(CONST_BLINT n, CONST_BLDOU x, CONST_BLINT incx,
					BLDOU y, CONST_BLINT incy);
	void BLAS_zaxpy(CONST_BLINT n, CONST_BLDOU a, CONST_BLDOU x, CONST_BLINT incx,
					BLDOU y, CONST_BLINT incy);
	void BLAS_dscal(CONST_BLINT n, CONST_BLDOU a, BLDOU x, CONST_BLINT incx);
	void BLAS_dtrsm(BLCHAR side, BLCHAR uplo, BLCHAR transa, BLCHAR diag, CONST_BLINT m,
					CONST_BLINT n, CONST_BLDOU alpha, CONST_BLDOU a, CONST_BLINT lda,
					BLDOU b, CONST_BLINT ldb);
	double BLAS_ddot(CONST_BLINT n, CONST_BLDOU x, CONST_BLINT incx, CONST_BLDOU y,
					 CONST_BLINT incy);
};


#endif /* CPPBLAS_H */

// Local Variables:
// mode:C++
// End:
