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

#include <dynlapack.h>

#include "Vector.hh"
#include "Matrix.hh"

class GeneralizedSchurDecomposition
{
private:
  const size_t n;
  const double criterium;
  lapack_int lwork;
  double *alphar, *alphai, *beta, *vsl, *work;
  lapack_int *bwork;
  static double criterium_static;
  static lapack_int selctg(const double *alphar, const double *alphai, const double *beta);
public:
  class GSDException
  {
  public:
    const lapack_int info, n;
    GSDException(lapack_int info_arg, lapack_int n_arg) : info(info_arg), n(n_arg) {};
  };
  //! \todo Replace heuristic choice for workspace size by a query to determine the optimal size
  GeneralizedSchurDecomposition(size_t n_arg, double criterium_arg);
  ~GeneralizedSchurDecomposition();
  //! \todo Add a lock around the modification of criterium_static for making it thread-safe
  template<class Mat1, class Mat2, class Mat3>
  void compute(Mat1 &S, Mat2 &T, Mat3 &Z, size_t &sdim) throw (GSDException);
  template<class Mat1, class Mat2, class Mat3, class Mat4, class Mat5>
  /*!
    \param[out] sdim Number of non-explosive generalized eigenvalues
  */
  void compute(const Mat1 &D, const Mat2 &E, Mat3 &S, Mat4 &T, Mat5 &Z, size_t &sdim) throw (GSDException);
  template<class Vec1, class Vec2>
  void getGeneralizedEigenvalues(Vec1 &eig_real, Vec2 &eig_cmplx);
};

std::ostream &operator<<(std::ostream &out, const GeneralizedSchurDecomposition::GSDException &e);

template<class Mat1, class Mat2, class Mat3>
void
GeneralizedSchurDecomposition::compute(Mat1 &S, Mat2 &T, Mat3 &Z, size_t &sdim) throw (GSDException)
{
  assert(S.getRows() == n && S.getCols() == n
         && T.getRows() == n && T.getCols() == n
         && Z.getRows() == n && Z.getCols() == n);

  lapack_int n2 = n;
  lapack_int info, sdim2;
  lapack_int lds = S.getLd(), ldt = T.getLd(), ldz = Z.getLd();

  criterium_static = criterium;
  // Here we are forced to give space for left Schur vectors, even if we don't use them, because of a bug in dgges()
  dgges("N", "V", "S", &selctg, &n2, S.getData(), &lds, T.getData(), &ldt,
        &sdim2, alphar, alphai, beta, vsl, &n2, Z.getData(), &ldz,
        work, &lwork, bwork, &info);

  if (info != 0)
    throw GSDException(info, n2);

  sdim = sdim2;
}

template<class Mat1, class Mat2, class Mat3, class Mat4, class Mat5>
void
GeneralizedSchurDecomposition::compute(const Mat1 &D, const Mat2 &E,
                                       Mat3 &S, Mat4 &T, Mat5 &Z, size_t &sdim) throw (GSDException)
{
  assert(D.getRows() == n && D.getCols() == n
         && E.getRows() == n && E.getCols() == n);
  S = D;
  T = E;
  compute(S, T, Z, sdim);
}

template<class Vec1, class Vec2>
void
GeneralizedSchurDecomposition::getGeneralizedEigenvalues(Vec1 &eig_real, Vec2 &eig_cmplx)
{
  assert(eig_real.getSize() == n && eig_cmplx.getSize() == n);

  double *par = alphar, *pai = alphai, *pb = beta,
    *per = eig_real.getData(), *pei = eig_cmplx.getData();
  while (par < alphar + n)
    {
      *per = *par++ / *pb;
      *pei = *pai++ / *pb++;
      per += eig_real.getStride();
      pei += eig_cmplx.getStride();
    }
}
