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

/****************************************************************
   % a wrapper class for function X=disclyap_fast(G,V,ch)
   %
   % Solve the discrete Lyapunov Equation
   % X=G*X*G'+V
   % Using the Doubling Algorithm
   %
   % If ch is defined then the code will check if the resulting X
   % is positive definite and generate an error message if it is not
   %
   % based on work of Joe Pearlman and Alejandro Justiniano
   % 3/5/2005
   % C++ version 28/07/09 by Dynare team
****************************************************************/

#if !defined(DiscLyapFast_INCLUDE)
#define DiscLyapFast_INCLUDE

#include "dynlapack.h"
#include "Matrix.hh"
#include "BlasBindings.hh"

class DiscLyapFast
{
  Matrix A0, A1, Ptmp, P0, P1, I;

public:
  class DLPException
  {
  public:
    const int info;
    std::string message;
    DLPException(int info_arg, std::string message_arg) :
      info(info_arg), message(message_arg) {
    };
  };

  DiscLyapFast(size_t n) :
    A0(n), A1(n), Ptmp(n), P0(n), P1(n), I(n)
  {
    mat::set_identity(I);
  };
  virtual ~DiscLyapFast(){};
  template <class MatG, class MatV, class MatX >
  void solve_lyap(const MatG &G, const MatV &V, MatX &X, double tol, size_t flag_ch) throw(DLPException);

};

template <class MatG, class MatV, class MatX >
void
DiscLyapFast::solve_lyap(const MatG &G, const MatV &V, MatX &X, double tol = 1e-16, size_t flag_ch = 0) throw(DLPException)
{
  P0 = V;
  P1 = V;
  A0 = G;

  const double alpha = 1.0;
  const double half = 0.5;
  const double omega = 0.0;

  bool matd = true;
  while (matd) // matrix diff > tol
    {
      //P1=P0+A0*P0*A0';
      // first step Ptmp=P0*A0';
      // DGEMM: C := alpha*op( A )*op( B ) + beta*C,
      blas::gemm("N", "T", alpha, P0, A0, omega, Ptmp);
      // P1=P0+A0*Ptmp;
      blas::gemm("N", "N", alpha, A0, Ptmp, alpha, P1);
      // A1=A0*A0;
      blas::gemm("N", "N", alpha, A0, A0, omega, A1);

      // ensure symmetry of P1=(P1+P1')/2;
      Ptmp = P1;
      blas::gemm("T", "N", half, Ptmp, I, half, P1);

      // check if max( max( abs( P1 - P0 ) ) )>tol
      matd = mat::isDiffSym(P1, P0, tol);
      P0 = P1;
      A0 = A1;
    } //end while

  // ensure symmetry of X=P0=(P0+P0')/2;
  blas::gemm("T", "N", half, P1, I, half, P0);
  X = P0;

  // Check that X is positive definite
  if (flag_ch == 1) // calc NormCholesky (P0)
    {
      lapack_int lpinfo = 0;
      lapack_int lrows = P0.getRows();
      lapack_int ldl = P0.getLd();
      dpotrf("L", &lrows, P0.getData(), &ldl, &lpinfo);
      if (lpinfo < 0)
        throw DLPException((int) lpinfo, std::string("DiscLyapFast:Internal error in NormCholesky calculator"));
      else if (lpinfo > 0)
        throw DLPException((int) lpinfo, std::string("DiscLyapFast:The matrix is not positive definite in NormCholesky calculator"));

    }
}

#endif //if !defined(DiscLyapFast_INCLUDE)
