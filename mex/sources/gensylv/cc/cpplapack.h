/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/cpplapack.h,v 1.3 2004/11/24 20:43:10 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef CPPLAPACK_H
#define CPPLAPACK_H

#if defined(MATLAB) && !defined(__linux__) && !defined(OCTAVE)
#define LAPACK_dgetrs dgetrs
#define LAPACK_dgetrf dgetrf
#define LAPACK_dgees  dgees
#define LAPACK_dgecon dgecon
#define LAPACK_dtrexc dtrexc
#define LAPACK_dtrsyl dtrsyl
#define LAPACK_dpotrf dpotrf
#define LAPACK_dgges  dgges
#define LAPACK_dsyev  dsyev
#else
#define LAPACK_dgetrs dgetrs_
#define LAPACK_dgetrf dgetrf_
#define LAPACK_dgees  dgees_
#define LAPACK_dgecon dgecon_
#define LAPACK_dtrexc dtrexc_
#define LAPACK_dtrsyl dtrsyl_
#define LAPACK_dpotrf dpotrf_
#define LAPACK_dgges  dgges_
#define LAPACK_dsyev  dsyev_
#endif

#define LACHAR const char*
#define CONST_LAINT const int*
#define LAINT int*
#define CONST_LADOU const double*
#define LADOU double*
typedef int (*DGGESCRIT)(const double*, const double*, const double*);

extern "C" {
	void LAPACK_dgetrs(LACHAR trans, CONST_LAINT n, CONST_LAINT nrhs,
					   CONST_LADOU a, CONST_LAINT lda, CONST_LAINT ipiv,
					   LADOU b, CONST_LAINT ldb, LAINT info);
	void LAPACK_dgetrf(CONST_LAINT m, CONST_LAINT n, LADOU a,
					   CONST_LAINT lda, LAINT ipiv, LAINT info);
	void  LAPACK_dgees(LACHAR jobvs, LACHAR sort, const void* select,
					   CONST_LAINT n, LADOU a, CONST_LAINT lda, LAINT sdim,
					   LADOU wr, LADOU wi, LADOU vs, CONST_LAINT ldvs,
					   LADOU work, CONST_LAINT lwork, const void* bwork, LAINT info);
	void LAPACK_dgecon(LACHAR norm, CONST_LAINT n, CONST_LADOU a, CONST_LAINT lda,
					   CONST_LADOU anorm, LADOU rnorm, LADOU work, LAINT iwork,
					   LAINT info);
	void LAPACK_dtrexc(LACHAR compq, CONST_LAINT n, LADOU t, CONST_LAINT ldt,
					   LADOU q, CONST_LAINT ldq, LAINT ifst, LAINT ilst, LADOU work,
					   LAINT info);
	void LAPACK_dtrsyl(LACHAR trana, LACHAR tranb, CONST_LAINT isgn, CONST_LAINT m,
					   CONST_LAINT n, CONST_LADOU a, CONST_LAINT lda, CONST_LADOU b,
					   CONST_LAINT ldb, LADOU c, CONST_LAINT ldc, LADOU scale,
					   LAINT info);
	void LAPACK_dpotrf(LACHAR uplo, CONST_LAINT n, LADOU a, CONST_LAINT lda,
					   LAINT info);
	void LAPACK_dgges(LACHAR jobvsl, LACHAR jobvsr, LACHAR sort, DGGESCRIT delztg,
					  CONST_LAINT n, LADOU a, CONST_LAINT lda, LADOU b, CONST_LAINT ldb,
					  LAINT sdim, LADOU alphar, LADOU alphai, LADOU beta,
					  LADOU vsl, CONST_LAINT ldvsl, LADOU vsr, CONST_LAINT ldvsr,
					  LADOU work, CONST_LAINT lwork, LAINT bwork, LAINT info);
	void LAPACK_dsyev(LACHAR jobz, LACHAR uplo, CONST_LAINT n, LADOU a, CONST_LAINT lda,
					  LADOU w, LADOU work, CONST_LAINT lwork, LAINT info); 
};


#endif /* CPPLAPACK_H */


// Local Variables:
// mode:C++
// End:

