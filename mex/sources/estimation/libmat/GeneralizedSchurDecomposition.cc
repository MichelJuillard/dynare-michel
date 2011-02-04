/*
 * Copyright (C) 2010-2011 Dynare Team
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

#include "GeneralizedSchurDecomposition.hh"

#include <cassert>
#include <cstdlib>

double GeneralizedSchurDecomposition::criterium_static;

GeneralizedSchurDecomposition::GeneralizedSchurDecomposition(size_t n_arg, double criterium_arg) :
  n(n_arg), criterium(criterium_arg)
{
  alphar = new double[n];
  alphai = new double[n];
  beta = new double[n];

  lwork = 16*n+16; // Same heuristic choice than mjdgges
  work = new double[(int) lwork];
  vsl = new double[n*n];

  bwork = new lapack_int[n];
}

GeneralizedSchurDecomposition::~GeneralizedSchurDecomposition()
{
  delete[] alphar;
  delete[] alphai;
  delete[] beta;
  delete[] work;
  delete[] vsl;
  delete[] bwork;
}

lapack_int
GeneralizedSchurDecomposition::selctg(const double *alphar, const double *alphai, const double *beta)
{
  return ((*alphar * *alphar + *alphai * *alphai) < criterium_static * *beta * *beta);
}

std::ostream &
operator<<(std::ostream &out, const GeneralizedSchurDecomposition::GSDException &e)
{
  out << "DGGES return code " << e.info << ": ";
  if (e.info < 0)
    out << "argument " << -e.info << " has an illegal value";
  else if (e.info <= e.n)
    out << "the QZ iteration failed";
  else if (e.info == e.n + 1)
    out << "other than QZ iteration failed in DHGEQZ";
  else if (e.info == e.n + 2)
    out << "after reordering, roundoff changed values of some complex eigenvalues so that leading eigenvalues in the Generalized Schur form no longer satisfy SELCTG=TRUE. This could also be caused due to scaling";
  else if (e.info == e.n + 3)
    out << "reordering failed in DTGSEN";
  return out;
}
