/*
 * Defines the prototypes for LAPACK Fortran functions.
 *
 * Also defines a typedef lapack_int to be used for all integers passed to LAPACK
 * functions.
 *
 * When used in the context of a MATLAB MEX file, you must define MATLAB_MEX_FILE
 * and MATLAB_VERSION (for version 7.4, define it to 0x0704).
 *
 *
 * Copyright (C) 2009-2010 Dynare Team
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

#ifndef _DYNLAPACK_H
#define _DYNLAPACK_H

/* Starting from version 7.8, MATLAB LAPACK expects ptrdiff_t arguments for integers */
#if defined(MATLAB_MEX_FILE) && MATLAB_VERSION >= 0x0708
# ifdef __cplusplus
#  include <cstdlib>
# else
#  include <stdlib.h>
# endif
typedef ptrdiff_t lapack_int;
#else
typedef int lapack_int;
#endif

#if defined(MATLAB_MEX_FILE) && defined(_WIN32)
# define FORTRAN_WRAPPER(x) x
#else
# define FORTRAN_WRAPPER(x) x ## _
#endif

#ifdef __cplusplus
extern "C" {
#endif

  typedef const char *LACHAR;
  typedef const lapack_int *CONST_LAINT;
  typedef lapack_int *LAINT;
  typedef const double *CONST_LADOU;
  typedef double *LADOU;
  typedef lapack_int (*DGGESCRIT)(const double *, const double *, const double *);

#define dgetrs FORTRAN_WRAPPER(dgetrs)
  void dgetrs(LACHAR trans, CONST_LAINT n, CONST_LAINT nrhs, CONST_LADOU a, CONST_LAINT lda, CONST_LAINT ipiv,
              LADOU b, CONST_LAINT ldb, LAINT info);

#define dgetrf FORTRAN_WRAPPER(dgetrf)
  void dgetrf(CONST_LAINT m, CONST_LAINT n, LADOU a,
              CONST_LAINT lda, LAINT ipiv, LAINT info);

#define dgees FORTRAN_WRAPPER(dgees)
  void dgees(LACHAR jobvs, LACHAR sort, const void *select,
             CONST_LAINT n, LADOU a, CONST_LAINT lda, LAINT sdim,
             LADOU wr, LADOU wi, LADOU vs, CONST_LAINT ldvs,
             LADOU work, CONST_LAINT lwork, const void *bwork, LAINT info);

#define dgecon FORTRAN_WRAPPER(dgecon)
  void dgecon(LACHAR norm, CONST_LAINT n, CONST_LADOU a, CONST_LAINT lda,
              CONST_LADOU anorm, LADOU rnorm, LADOU work, LAINT iwork,
              LAINT info);

#define dtrexc FORTRAN_WRAPPER(dtrexc)
  void dtrexc(LACHAR compq, CONST_LAINT n, LADOU t, CONST_LAINT ldt,
              LADOU q, CONST_LAINT ldq, LAINT ifst, LAINT ilst, LADOU work,
              LAINT info);

#define dtrsyl FORTRAN_WRAPPER(dtrsyl)
  void dtrsyl(LACHAR trana, LACHAR tranb, CONST_LAINT isgn, CONST_LAINT m,
              CONST_LAINT n, CONST_LADOU a, CONST_LAINT lda, CONST_LADOU b,
              CONST_LAINT ldb, LADOU c, CONST_LAINT ldc, LADOU scale,
              LAINT info);

#define dpotrf FORTRAN_WRAPPER(dpotrf)
  void dpotrf(LACHAR uplo, CONST_LAINT n, LADOU a, CONST_LAINT lda,
              LAINT info);

#define dgges FORTRAN_WRAPPER(dgges)
  void dgges(LACHAR jobvsl, LACHAR jobvsr, LACHAR sort, DGGESCRIT delztg,
             CONST_LAINT n, LADOU a, CONST_LAINT lda, LADOU b, CONST_LAINT ldb,
             LAINT sdim, LADOU alphar, LADOU alphai, LADOU beta,
             LADOU vsl, CONST_LAINT ldvsl, LADOU vsr, CONST_LAINT ldvsr,
             LADOU work, CONST_LAINT lwork, LAINT bwork, LAINT info);

#define dsyev FORTRAN_WRAPPER(dsyev)
  void dsyev(LACHAR jobz, LACHAR uplo, CONST_LAINT n, LADOU a, CONST_LAINT lda,
             LADOU w, LADOU work, CONST_LAINT lwork, LAINT info);

#define dsyevr FORTRAN_WRAPPER(dsyevr)
  void dsyevr(LACHAR jobz, LACHAR range, LACHAR uplo, CONST_LAINT n, LADOU a,
              CONST_LAINT lda, LADOU lv, LADOU vu, CONST_LAINT il, CONST_LAINT iu,
              CONST_LADOU abstol, LAINT m, LADOU w, LADOU z, CONST_LAINT ldz,
              LAINT isuppz, LADOU work, CONST_LAINT lwork, LAINT iwork, CONST_LAINT liwork,
              LAINT info);

#define dgeqrf FORTRAN_WRAPPER(dgeqrf)
  void dgeqrf(CONST_LAINT m, CONST_LAINT n, LADOU a, CONST_LAINT lda,
              LADOU tau, LADOU work, CONST_LAINT lwork, LAINT info);

#define dormqr FORTRAN_WRAPPER(dormqr)
  void dormqr(LACHAR side, LACHAR trans, CONST_LAINT m, CONST_LAINT n, CONST_LAINT k,
              CONST_LADOU a, CONST_LAINT lda, CONST_LADOU tau, LADOU c, CONST_LAINT ldc,
              LADOU work, CONST_LAINT lwork, LAINT info);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _DYNLAPACK_H */
