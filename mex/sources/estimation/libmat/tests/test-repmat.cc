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
#include "Matrix.hh"

int
main(int argc, char **argv)
{
  size_t m = 4, n = 3;
  Matrix S(m, n), A(5*m, 3*n), A2(5*m, 3*n);

  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < n; j++)
      S(i, j) = i*n + j + 1;

  std::cout << "Matrix to be tiled:" << std::endl << S << std::endl;

  Matrix   S1(m, 2 *n), S2(2*m, n);
  std::vector<size_t>  toC(n), toR(m);
  for (int i = 0; i < n; i++)
    toC[i] = 2*i;

  for (int j = 0; j < m; j++)
    toR[j] = m-j-1;

  mat::assignByVectors(S1, mat::nullVec, toC, S, mat::nullVec, mat::nullVec);
  std::cout << "Matrix assigned and col reorder by vectors:" << std::endl << S1 << std::endl;

  mat::assignByVectors(S2, toR, mat::nullVec, S, mat::nullVec, mat::nullVec);
  std::cout << "Matrix assigned and row ereorder by vectors:" << std::endl << S2 << std::endl;

  mat::assignByVectors(S1, toR, toC, S, mat::nullVec, mat::nullVec);
  std::cout << "Matrix assigned by vectors:" << std::endl << S1 << std::endl;

  mat::repmat(S, 5, 3, A);
  std::cout << "Tiled=" << std::endl << A << std::endl;

  A2 = A;
  bool b = mat::isDiff(A, A2, 0.0000001);
  std::cout << "No Diff=" << b << std::endl;
  mat::add(A, 0.000001);
  b = mat::isDiff(A, A2, 0.0000001);
  std::cout << "Yes Diff=" << b << std::endl;
  Matrix A4(5), A5(5);
  A4.setAll(1.0);
  A5.setAll(1.0);
  b = mat::isDiffSym(A4, A5, 0.0000001);
  std::cout << "No DiffSym=" << b << std::endl;
  mat::sub(A4, 0.0001);
  b = mat::isDiffSym(A4, A5, 0.0000001);
  std::cout << "Yes DiffSym=" << b << std::endl;

}
