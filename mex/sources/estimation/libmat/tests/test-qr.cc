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
  Matrix S(m, n), Q(m), A(m, n), B(m);
  QRDecomposition QRD(m, n);

  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < n; j++)
      S(i, j) = i*n + j + 1;

  std::cout << "Matrix to be decomposed:" << std::endl << S << std::endl;

  mat::set_identity(Q);

  QRD.computeAndMultByQ(S, "L", "N", Q);

  std::cout << "Matrix Q:" << std::endl << Q << std::endl;

  blas::gemm("T", "N", 1.0, Q, Q, 0.0, B);

  std::cout << "Matrix Q'*Q:" << std::endl << B << std::endl;

  for (size_t j = 0; j < n; j++)
    mat::col_set(S, j, j+1, m-j-1, 0);

  std::cout << "Matrix R:" << std::endl << S << std::endl;

  blas::gemm("N", "N", 1.0, Q, S, 0.0, A);

  std::cout << "Product Q*R:" << std::endl << A << std::endl;
}
