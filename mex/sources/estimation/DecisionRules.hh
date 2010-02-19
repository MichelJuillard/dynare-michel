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

#include <cstdlib>
#include <vector>

#include "Vector.hh"
#include "Matrix.hh"
#include "QRDecomposition.hh"
#include "GeneralizedSchurDecomposition.hh"
#include "LUSolver.hh"

class DecisionRules
{
private:
  const size_t n, p;
  const std::vector<size_t> zeta_fwrd, zeta_back, zeta_mixed, zeta_static;
  const size_t n_fwrd, n_back, n_mixed, n_static, n_back_mixed, n_fwrd_mixed, n_dynamic;
  std::vector<size_t> zeta_fwrd_mixed, zeta_back_mixed, zeta_dynamic,
    beta_back, beta_fwrd, pi_back, pi_fwrd;
  Matrix S, A, D, E, Z_prime;
  QRDecomposition QR;
  GeneralizedSchurDecomposition GSD;
  LUSolver LU1, LU2, LU3;
  Matrix Z21, g_y_back, g_y_back_tmp;
  Matrix g_y_static, A0s, A0d, g_y_dynamic, g_y_static_tmp;
public:
  class BlanchardKahnException
  {
  public:
    //! True if the model fails the order condition. False if it fails the rank condition.
    const bool order;
    const int n_fwrd_vars, n_explosive_eigenvals;
    BlanchardKahnException(bool order_arg, int n_fwrd_vars_arg, int n_explosive_eigenvals_arg) : order(order_arg), n_fwrd_vars(n_fwrd_vars_arg), n_explosive_eigenvals(n_explosive_eigenvals_arg) {};
  };
  /*!
    The zetas are supposed to follow C convention (first vector index is zero).
  */
  DecisionRules(size_t n_arg, size_t p_arg, const std::vector<size_t> &zeta_fwrd_arg,
                const std::vector<size_t> &zeta_back_arg, const std::vector<size_t> &zeta_mixed_arg,
                const std::vector<size_t> &zeta_static_arg, double qz_criterium);
  /*!
    \param jacobian First columns are backetermined vars at t-1 (in the order of zeta_back_mixed), then all vars at t (in the orig order), then forward vars at t+1 (in the order of zeta_fwrd_mixed), then exogenous vars.
  */
  void compute(const Matrix &jacobian, Matrix &g_y, Matrix &g_u) throw (BlanchardKahnException, GeneralizedSchurDecomposition::GSDException);
  template<class Vec1, class Vec2>
  void getGeneralizedEigenvalues(Vec1 &eig_real, Vec2 &eig_cmplx);
};

std::ostream &operator<<(std::ostream &out, const DecisionRules::BlanchardKahnException &e);

template<class Vec1, class Vec2>
void
DecisionRules::getGeneralizedEigenvalues(Vec1 &eig_real, Vec2 &eig_cmplx)
{
  GSD.getGeneralizedEigenvalues(eig_real, eig_cmplx);
}
