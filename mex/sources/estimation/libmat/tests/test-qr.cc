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

#include <iostream>

#include "BlasBindings.hh"
#include "QRDecomposition.hh"

int
main(int argc, char **argv)
{
  size_t m = 4, n = 3;
  Matrix S(m, n), Q(m), A(m, n), B(m), S2(m, n);
  QRDecomposition QRD(m, n, m);

  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < n; j++)
      S(i, j) = i*n + j + 1;

  std::cout << "Matrix to be decomposed:" << std::endl << S << std::endl;

  mat::set_identity(Q);

  S2 = S;
  QRD.computeAndLeftMultByQ(S2, "N", Q);

  std::cout << "Q =" << std::endl << Q << std::endl;

  blas::gemm("T", "N", 1.0, Q, Q, 0.0, B);

  std::cout << "Q'*Q =" << std::endl << B << std::endl;

  for (size_t j = 0; j < n; j++)
    mat::col_set(S2, j, j+1, m-j-1, 0);

  std::cout << "R =" << std::endl << S2 << std::endl;

  blas::gemm("N", "N", 1.0, Q, S2, 0.0, A);

  std::cout << "Q*R =" << std::endl << A << std::endl;

  // Check values
  Matrix B2(m);
  mat::set_identity(B2);
  mat::sub(B2, B);
  assert(mat::nrminf(B2) < 1e-4);

  mat::sub(A, S);
  assert(mat::nrminf(A) < 1e-4);
}
