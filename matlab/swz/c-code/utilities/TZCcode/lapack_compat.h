// Created: August 4,2009

#ifndef __LAPACKCOMPAT__
#define __LAPACKCOMPAT__

#define USE_LAPACK

#if defined(USE_LAPACK)
   #include "blas_lapack.h"
#endif

#define dgetrf dgetrf_
#define dgesv dgesv_
#define dpotrf dpotrf_
#define dsyev dsyev_
#define dgeev dgeev_
#define dpotri dpotri_
#define vdDiv vdDiv_
#define vdInv vdInv_
#define vdSqrt vdSqrt_
#define vdLn vdLn_
#define vdExp vdExp_
#endif

