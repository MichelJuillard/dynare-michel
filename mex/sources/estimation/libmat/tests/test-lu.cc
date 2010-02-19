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
#include "LUSolver.hh"

int
main(int argc, char **argv)
{
  size_t m = 4, n = 3;

  double A_data[] = { -1, 2, 3,
                      4, -5, 6,
                      7, 8, -9 };
  double B_data[] = { 1, -3, 4, 5,
                      -7, 9, 1, 7,
                      -3, 4, 0, -2 };
  MatrixView A(A_data, n, n, n), B_prime(B_data, m, n, m);
  Matrix B(n, m);

  mat::transpose(A);
  mat::transpose(B, B_prime);

  std::cout << "A =" << std::endl << A << std::endl
            << "B =" << std::endl << B << std::endl;

  LUSolver LU(n);

  LU.invMult("N", A, B);

  std::cout << "A\\B =" << std::endl << B;

  // Check this is the right result
  double C_data[] = { -1.0500, 1.3750, 0.0833, 0.6250,
                      0.3000, -0.7500, 0.8333, 0.7500,
                      -0.2167, -0.0417, 0.8056, 1.3750 };
  MatrixView C_prime(C_data, m, n, m);
  Matrix C(n, m);
  mat::transpose(C, C_prime);
  mat::sub(B, C);
  assert(mat::nrminf(B) < 1e-4);
}
