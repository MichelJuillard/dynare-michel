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

#include "GeneralizedSchurDecomposition.hh"

int
main(int argc, char **argv)
{
  size_t n = 3;
  double D_data[] = { 1, 2, 3,
                      4, 5, 6,
                      7, 8, 9 };
  double E_data[] = { 1, -3, 4,
                      -7, 9, 1,
                      -3, 4, 0 };
  MatrixView D(D_data, n, n, n), E(E_data, n, n, n);

  // Need to transpose because internally matrices are in column-major order
  mat::transpose(D);
  mat::transpose(E);

  std::cout << "D =" << std::endl << D << std::endl;
  std::cout << "E =" << std::endl << E << std::endl;

  GeneralizedSchurDecomposition GSD(n, 1.00001);

  Matrix S(n), T(n), Z(n);
  size_t sdim;

  GSD.compute(D, E, S, T, Z, sdim);

  std::cout << "S =" << std::endl << S << std::endl;
  std::cout << "T =" << std::endl << T << std::endl;
  std::cout << "Z =" << std::endl << Z << std::endl;

  Vector eig_real(n), eig_cmplx(n);
  GSD.getGeneralizedEigenvalues(eig_real, eig_cmplx);

  std::cout << "Real part of generalized eigenvalues: " << std::endl << eig_real << std::endl;
  std::cout << "Complex part of generalized eigenvalues: " << std::endl << eig_cmplx << std::endl;

  return 0;
}
