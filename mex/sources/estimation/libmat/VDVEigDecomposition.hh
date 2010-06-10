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

/**
 *  Calculates vector of Eigen values D and matrix of Eigen vectors V
 *  of input matrix m
 */

#if !defined(VDVEigDec_INCLUDE)
#define VDVEigDec_INCLUDE

#include <cstdlib>
#include "Vector.hh"
#include "Matrix.hh"
#include <dynlapack.h>

class VDVEigDecomposition
{
  lapack_int lda, n;
  lapack_int lwork, info;
  double tmpwork;
  double *work;
  bool converged;
  Matrix V;
  Vector D;
public:

public:
  class VDVEigException
  {
  public:
    const lapack_int info;
    std::string message;
    VDVEigException(lapack_int info_arg, std::string message_arg) :
      info(info_arg), message(message_arg) {
    };
  };

/**
 *  This constructor only creates optimal workspace using
 *  input matrix m
 */
  VDVEigDecomposition(const Matrix &m) throw(VDVEigException);

/**
 *  This constructoro only crates workspace using the size of
 *  the input matrix m
 */
  VDVEigDecomposition(size_t n) throw(VDVEigException);

  virtual ~VDVEigDecomposition()
  {
    delete[] work;
  };
  template <class Mat>
  void calculate(const Mat &H) throw(VDVEigException);
  // get eigenvalues
  Vector &
  getD()
  {
    return D;
  };
  // get eigen vector
  Matrix &
  getV()
  {
    return V;
  };
  // check if converged
  bool
  hasConverged()
  {
    return converged;
  };
};

template <class Mat>
void
VDVEigDecomposition::calculate(const Mat &m) throw(VDVEigException)
{
  info = 0;
  if (m.getRows() != m.getCols())
    throw(VDVEigException(info, "Matrix is not square in VDVEigDecomposition calculate"));

  if (m.getCols() != (size_t) n  || m.getLd() != (size_t) lda)
    throw(VDVEigException(info, "Matrix not matching VDVEigDecomposition class"));
  lapack_int tmplwork = -1;
  V = m;
  dsyev("V", "U", &n, V.getData(), &lda, D.getData(), &tmpwork, &tmplwork, &info);
  if (lwork < tmpwork)
    {
      lwork = tmpwork;
      delete[] work;
      work = new double[lwork];
    }
  dsyev("V", "U", &n, V.getData(), &lda, D.getData(), work, &lwork, &info);

  if (info < 0)
    throw(VDVEigException(info, "Internal error in VDVEigDecomposition calculation"));
  converged = true;
  if (info)
    converged = false;
}

std::ostream &operator<<(std::ostream &out, const VDVEigDecomposition::VDVEigException &e);

#endif
