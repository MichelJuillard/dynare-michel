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

#include "VDVEigDecomposition.hh"

VDVEigDecomposition::VDVEigDecomposition(const Matrix &m) throw(VDVEigException) :
  lda(m.getLd()), n(m.getCols()), lwork(-1),
  info(0),   converged(false), V(m), D(n)
{
  if (m.getRows() != m.getCols())
    throw(VDVEigException(info, "Matrix is not square in VDVEigDecomposition constructor"));

  double tmpwork;
  dsyev("V", "U", &n, V.getData(), &lda, D.getData(), &tmpwork, &lwork, &info);
  lwork = (lapack_int) tmpwork;
  work = new double[lwork];
  if (info < 0)
    throw(VDVEigException(info, "Internal error in VDVEigDecomposition constructor"));
}

VDVEigDecomposition::VDVEigDecomposition(size_t inn) throw(VDVEigException) :
  lda(inn), n(inn), lwork(3*inn-1),
  info(0),   converged(false), V(inn), D(inn)
{
  work = new double[lwork];
};

std::ostream &
operator<<(std::ostream &out, const VDVEigDecomposition::VDVEigException &e)
{
  out << " info " << e.info << ": " << e.message  << "\n";
  return out;
}
