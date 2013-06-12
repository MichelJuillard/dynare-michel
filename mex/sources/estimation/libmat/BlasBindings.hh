/*
 * Copyright (C) 2010-2013 Dynare Team
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
  /* Level 1 */

  //! dot product of two vectors
  template<class Vec1, class Vec2>
  inline double
  dot(const Vec1 &A, Vec2 &B)
  {
    assert(A.getSize() == B.getSize());
    blas_int n = A.getSize();
    blas_int lda = A.getStride(), ldb = B.getStride();
    return ddot(&n, A.getData(), &lda, B.getData(), &ldb);
  }

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

  //! General matrix * vector multiplication
  //  c = alpha*A*b + beta*c,   or   c := alpha*A'*b + beta*c,
  // where alpha and beta are scalars, b and c are vectors and A is an
  // m by n matrix.
  template<class Mat1, class Vec2, class Vec3>
  inline void
  gemv(const char *transa, double alpha, const Mat1 &A,
       const Vec2 &B, double beta, Vec3 &C)
  {
    blas_int m = A.getRows(), n = A.getCols();
    if (*transa == 'T')
      {
        assert(C.getSize() == A.getCols());
        assert(B.getSize() == A.getRows());
      }
    else
      {
        assert(C.getSize() == A.getRows());
        assert(B.getSize() == A.getCols());
      }
    blas_int lda = A.getLd(), ldb = B.getStride(), ldc = C.getStride();
    dgemv(transa, &m, &n, &alpha, A.getData(), &lda,
          B.getData(), &ldb, &beta, C.getData(), &ldc);
  }

  //! Symmetric matrix * vector multiplication
  //  c = alpha*A*b + beta*c,
  // where alpha and beta are scalars, b and c are vectors and A is a
  // m by m symmetric matrix.
  template<class Mat1, class Vec2, class Vec3>
  inline void
  symv(const char *uplo, double alpha, const Mat1 &A,
       const Vec2 &B, double beta, Vec3 &C)
  {
    assert(A.getRows() == A.getCols());
    blas_int n = A.getRows();
    assert(A.getRows() == B.getSize());
    assert(A.getRows() == C.getSize());
    blas_int lda = A.getLd(), ldb = B.getStride(), ldc = C.getStride();
    dsymv(uplo, &n,  &alpha, A.getData(), &lda,
          B.getData(), &ldb, &beta, C.getData(), &ldc);
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

  //! Symmetric matrix A * (poss. rectangular) matrix B multiplication
  template<class Mat1, class Mat2, class Mat3>
  inline void
  symm(const char *side, const char *uplo,
       double alpha, const Mat1 &A, const Mat2 &B,
       double beta, Mat3 &C)
  {
    assert(A.getRows() == A.getCols());
    assert(B.getRows() == C.getRows());
    assert(B.getCols() == C.getCols());
    if (*side == 'L' || *side == 'l')
      assert(A.getCols() == B.getRows());
    else if (*side == 'R' || *side == 'r')
      assert(A.getRows() == B.getCols());

    blas_int m = B.getRows(), n = B.getCols();
    blas_int lda = A.getLd(), ldb = B.getLd(), ldc = C.getLd();
    dsymm(side, uplo, &m, &n, &alpha, A.getData(), &lda,
          B.getData(), &ldb, &beta, C.getData(), &ldc);
  }

} // End of namespace

#endif
