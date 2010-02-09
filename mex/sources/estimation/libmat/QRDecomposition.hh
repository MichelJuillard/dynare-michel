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

#include <algorithm> // For std::min()

#include <dynlapack.h>

#include "Vector.hh"
#include "Matrix.hh"
#include "BlasBindings.hh"

class QRDecomposition
{
private:
  const size_t rows, cols, mind;
  lapack_int lwork;
  double *work, *tau;
  Matrix H, Q2;
  Vector v;
public:
  /*!
    \todo Replace heuristic choice for workspace size by a query to determine the optimal size
   */
  QRDecomposition(size_t rows_arg, size_t cols_arg);
  ~QRDecomposition();
  //! Performs the QR decomposition of a matrix, and multiplies another matrix by Q
  /*!
    \param[in,out] A On input, the matrix to be decomposed. On output, equals to the output of dgeqrf
    \param[in] side Specifies on which side is Q in the multiplication, either "L" or "R"
    \param[in] trans Specifies whether Q should be transposed before the multiplication, either "T" or "N"
    \param[in,out] C The matrix to be multiplied by Q, modified in place
  */
  template<class Mat1, class Mat2>
  void computeAndMultByQ(Mat1 &A, const char *side, const char *trans, Mat2 &C);
};

template<class Mat1, class Mat2>
void
QRDecomposition::computeAndMultByQ(Mat1 &A, const char *side, const char *trans, Mat2 &C)
{
  assert(A.getRows() == rows && A.getCols() == cols);
  assert(C.getRows() == rows);

  lapack_int m = rows, n = cols, lda = A.getLd();
  lapack_int info;
  dgeqrf(&m, &n, A.getData(), &lda, tau, work, &lwork, &info);
  assert(info == 0);

  n = C.getCols();
  lapack_int k = mind, ldc = C.getLd();
  dormqr(side, trans, &m, &n, &k, A.getData(), &lda, tau, C.getData(), &ldc,
         work, &lwork, &info);
  assert(info == 0);
}
