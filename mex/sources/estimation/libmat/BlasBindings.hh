/*
 * Copyright (C) 2010 Dynare Team
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

#ifndef _BLAS_BINDINGS_HH
#define _BLAS_BINDINGS_HH

#include <dynblas.h>

#include "Vector.hh"
#include "Matrix.hh"

namespace blas
{
  /* Level 2 */
  
  //! Symmetric rank 1 operation: A = alpha*X*X' + A
  template<class Mat, class Vec>
  inline void
  syr(const char *uplo, double alpha, Vec X, Mat A)
  {
    assert(X.getSize() == A.getCols() && A.getCols() == A.getRows());
    blas_int n = X.getSize();
    blas_int incx = X.getStride();
    blas_int lda = A.getLd();
    dsyr(uplo, &n, &alpha, X.getData(), &incx, A.getData(), &lda);
  }

  /* Level 3 */

  //! General matrix multiplication
  template<class Mat1, class Mat2, class Mat3>
  inline void
  gemm(const char *transa, const char *transb,
       double alpha, const Mat1 &A, const Mat2 &B,
       double beta, Mat3 &C)
  {
    blas_int m = A.getRows(), n = B.getCols(), k = A.getCols();
    if (*transa == 'N')
      {
        if (*transb == 'N')
          {
            assert(A.getRows() == C.getRows());
            assert(A.getCols() == B.getRows());
            assert(B.getCols() == C.getCols());
          }
        else if (*transb == 'T')
          {
            assert(A.getRows() == C.getRows());
            assert(A.getCols() == B.getCols());
            assert(B.getRows() == C.getCols());
            n = B.getRows();
          }
      }
    else if (*transa == 'T')
      {
        m = A.getCols();
        k = A.getRows();
        if (*transb == 'N')
          {
            assert(A.getCols() == C.getRows());
            assert(A.getRows() == B.getRows());
            assert(B.getCols() == C.getCols());
          }
        else if (*transb == 'T')
          {
            assert(A.getCols() == C.getRows());
            assert(A.getRows() == B.getCols());
            assert(B.getRows() == C.getCols());
            n = B.getRows();
          }
      }
    blas_int lda = A.getLd(), ldb = B.getLd(), ldc = C.getLd();
    dgemm(transa, transb, &m, &n, &k, &alpha, A.getData(), &lda,
          B.getData(), &ldb, &beta, C.getData(), &ldc);
  }

  //! Symmetric matrix multiplication
  template<class Mat1, class Mat2, class Mat3>
  inline void
  symm(const char *side, const char *uplo,
       double alpha, const Mat1 &A, const Mat2 &B,
       double beta, Mat3 &C)
  {
    assert(A.getRows() == A.getCols());
    assert(A.getRows() == C.getRows());
    assert(A.getCols() == B.getRows());
    assert(B.getCols() == C.getCols());
    blas_int m = A.getRows(), n = B.getCols();
    blas_int lda = A.getLd(), ldb = B.getLd(), ldc = C.getLd();
    dsymm(side, uplo, &m, &n, &alpha, A.getData(), &lda,
          B.getData(), &ldb, &beta, C.getData(), &ldc);
  }
} // End of namespace

#endif
