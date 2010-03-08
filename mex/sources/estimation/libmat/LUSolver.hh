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

#include <cstdlib>
#include <cassert>

#include <dynlapack.h>

class LUSolver
{
private:
  const size_t dim;
  lapack_int *ipiv;
public:
  class LUException
  {
  public:
    const lapack_int info;
    LUException(lapack_int info_arg) : info(info_arg) {};
  };
  LUSolver(size_t dim_arg);
  virtual ~LUSolver();
  /*!
    Computes A^(-1)*B (possibly transposing A).
    The output is stored in B.
    The input matrix A is modified.
  */
  template<class Mat1, class Mat2>
  void invMult(const char *trans, Mat1 &A, Mat2 &B) throw (LUException);
};

template<class Mat1, class Mat2>
void
LUSolver::invMult(const char *trans, Mat1 &A, Mat2 &B) throw (LUException)
{
  assert(A.getRows() == dim && A.getCols() == dim);
  assert(B.getRows() == dim);
  lapack_int n = dim, lda = A.getLd(), info;
  dgetrf(&n, &n, A.getData(), &lda, ipiv, &info);

  if (info != 0)
    throw LUException(info);

  lapack_int nrhs = B.getCols(), ldb = B.getLd();
  dgetrs(trans, &n, &nrhs, A.getData(), &lda, ipiv, B.getData(), &ldb, &info);
  assert(info == 0);
}
