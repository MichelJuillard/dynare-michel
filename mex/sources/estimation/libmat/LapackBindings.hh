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

#ifndef _LAPACK_BINDINGS_HH
#define _LAPACK_BINDINGS_HH

#include <dynlapack.h>

#include "Vector.hh"
#include "Matrix.hh"

namespace lapack
{
  // calc Cholesky Decomposition (Mat A, char "U"pper/"L"ower)
  template<class Mat>
  inline int
  choleskyDecomp(Mat &A, const char *UL)
  {
    assert(A.getCols() == A.getRows());
    lapack_int lpinfo = 0;
    lapack_int lrows = A.getRows();
    lapack_int ldl = A.getLd();
    dpotrf(UL, &lrows, A.getData(), &ldl, &lpinfo);
    int info = (int) lpinfo;
    return info;
  }

  // calc Cholesky Decomposition based solution X to A*X=B
  // for A pos. def. and symmetric supplied as uppper/lower triangle
  // packed in a vector if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
  // solutino X is returned in B and if B_in=I then B_out=X=inv(A)
  // A_out contains uupper or lower Cholesky decomposition
  template<class VecA, class MatB>
  inline int
  choleskySolver(VecA &A, MatB &B, const char *UL)
  {
    //assert(A.getCols() == A.getRows());
    lapack_int lpinfo = 0;
    lapack_int size = B.getRows();
    lapack_int bcols = B.getCols();
    lapack_int ldl = B.getLd();
    dppsv(UL, &size, &bcols, A.getData(), B.getData(), &ldl, &lpinfo);
    int info = (int) lpinfo;
    return info;
  }

} // End of namespace

#endif
