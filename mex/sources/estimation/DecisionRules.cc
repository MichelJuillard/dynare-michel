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

#include <cassert>

#include <algorithm>

#include "DecisionRules.hh"

DecisionRules::DecisionRules(size_t n_arg, size_t p_arg,
                             const std::vector<size_t> &zeta_fwrd_arg,
                             const std::vector<size_t> &zeta_back_arg,
                             const std::vector<size_t> &zeta_mixed_arg,
                             const std::vector<size_t> &zeta_static_arg,
                             double qz_criterium) :
  n(n_arg), p(p_arg), zeta_fwrd(zeta_fwrd_arg), zeta_back(zeta_back_arg),
  zeta_mixed(zeta_mixed_arg), zeta_static(zeta_static_arg),
  n_fwrd(zeta_fwrd.size()), n_back(zeta_back.size()),
  n_mixed(zeta_mixed.size()), n_static(zeta_static.size()),
  n_back_mixed(n_back+n_mixed), n_fwrd_mixed(n_fwrd+n_mixed),
  n_dynamic(n-n_static),
  S(n, n_static),
  A(n, n_back_mixed + n + n_fwrd_mixed),
  D(n_fwrd + n_back + 2*n_mixed),
  E(n_fwrd + n_back + 2*n_mixed),
  Z_prime(n_fwrd + n_back + 2*n_mixed),
  QR(n, n_static, n_back_mixed + n + n_fwrd_mixed),
  GSD(n_fwrd + n_back + 2*n_mixed, qz_criterium),
  LU1(n_fwrd_mixed),
  LU2(n_back_mixed),
  LU3(n_static),
  Z21(n_fwrd_mixed, n_back_mixed),
  g_y_back(n_back_mixed),
  g_y_back_tmp(n_back_mixed),
  g_y_static(n_static, n_back_mixed),
  A0s(n_static),
  A0d(n_static, n_dynamic),
  g_y_dynamic(n_dynamic, n_back_mixed),
  g_y_static_tmp(n_fwrd_mixed, n_back_mixed),
  g_u_tmp1(n, n_back_mixed),
  g_u_tmp2(n),
  LU4(n)
{
  assert(n == n_back + n_fwrd + n_mixed + n_static);

  set_union(zeta_fwrd.begin(), zeta_fwrd.end(),
            zeta_mixed.begin(), zeta_mixed.end(),
            back_inserter(zeta_fwrd_mixed));
  set_union(zeta_back.begin(), zeta_back.end(),
            zeta_mixed.begin(), zeta_mixed.end(),
            back_inserter(zeta_back_mixed));
  set_union(zeta_back_mixed.begin(), zeta_back_mixed.end(),
            zeta_fwrd.begin(), zeta_fwrd.end(),
            back_inserter(zeta_dynamic));

  // Compute beta_back and pi_back
  for(size_t i = 0; i < n_back_mixed; i++)
    if (find(zeta_mixed.begin(), zeta_mixed.end(), zeta_back_mixed[i])
        == zeta_mixed.end())
      pi_back.push_back(i);
    else
      beta_back.push_back(i);

  // Compute beta_fwrd and pi_fwrd
  for(size_t i = 0; i < n_fwrd_mixed; i++)
    if (find(zeta_mixed.begin(), zeta_mixed.end(), zeta_fwrd_mixed[i])
        == zeta_mixed.end())
      pi_fwrd.push_back(i);
    else
      beta_fwrd.push_back(i);

  // Construct the fixed part of D and E
  D.setAll(0.0);
  for (size_t i = 0; i < n_mixed; i++)
    D(n - n_static + i, beta_back[i]) = 1.0;
  E.setAll(0.0);
  for (size_t i = 0; i < n_mixed; i++)
    E(n - n_static + i, n_back_mixed + beta_fwrd[i]) = 1.0;
}

void
DecisionRules::compute(const Matrix &jacobian, Matrix &g_y, Matrix &g_u) throw (BlanchardKahnException, GeneralizedSchurDecomposition::GSDException)
{
  assert(jacobian.getRows() == n
         && jacobian.getCols() == (n_back_mixed + n + n_fwrd_mixed + p));
  assert(g_y.getRows() == n && g_y.getCols() == n_back_mixed);
  assert(g_u.getRows() == n && g_u.getCols() == p);

  // Construct S, perform QR decomposition and get A = Q*jacobian
  for(size_t i = 0; i < n_static; i++)
    mat::col_copy(jacobian, n_back_mixed+zeta_static[i], S, i);

  A = MatrixConstView(jacobian, 0, 0, n, n_back_mixed + n + n_fwrd_mixed);
  QR.computeAndLeftMultByQ(S, "T", A);

  // Construct matrix D
  for (size_t j = 0; j < n_back_mixed; j++)
    mat::col_copy(A, n_back_mixed + zeta_back_mixed[j], n_static, n - n_static,
                  D, j, 0);
  MatrixView(D, 0, n_back_mixed, n - n_static, n_fwrd_mixed) = MatrixView(A, n_static, n_back_mixed + n, n - n_static, n_fwrd_mixed);

  // Construct matrix E
  MatrixView(E, 0, 0, n - n_static, n_back_mixed) = MatrixView(A, n_static, 0, n - n_static, n_back_mixed);
  for (size_t j = 0; j < n_fwrd; j++)
    mat::col_copy(A, n_back_mixed + zeta_fwrd_mixed[pi_fwrd[j]], n_static, n - n_static,
                  E, n_back_mixed + pi_fwrd[j], 0);
  MatrixView E_tmp(E, 0, 0, n - n_static, n_fwrd + n_back + 2*n_mixed);
  mat::negate(E_tmp); // Here we take the opposite of some of the zeros initialized in the constructor, but it is not a problem

  // Perform the generalized Schur
  size_t sdim;
  GSD.compute(E, D, Z_prime, sdim);

  if (n_back_mixed != sdim)
    throw BlanchardKahnException(true, n_fwrd_mixed, n_fwrd + n_back + 2*n_mixed - sdim);

  // Compute DR for forward variables w.r. to endogenous
  MatrixView Z21_prime(Z_prime, 0, n_back_mixed, n_back_mixed, n_fwrd_mixed),
    Z22_prime(Z_prime, n_back_mixed, n_back_mixed, n_fwrd_mixed, n_fwrd_mixed);

  mat::transpose(Z21, Z21_prime);

  try
    {
      LU1.invMult("T", Z22_prime, Z21);
    }
  catch (LUSolver::LUException &e)
    {
      throw BlanchardKahnException(false, n_fwrd_mixed, n_fwrd + n_back + 2*n_mixed - sdim);
    }
  mat::negate(Z21);
  const Matrix &g_y_fwrd = Z21;

  for (size_t i = 0; i < n_fwrd_mixed; i++)
    mat::row_copy(g_y_fwrd, i, g_y, zeta_fwrd_mixed[i]);

  // Compute DR for backward variables w.r. to endogenous
  MatrixView Z11_prime(Z_prime, 0, 0, n_back_mixed, n_back_mixed),
    T11(D, 0, 0, n_back_mixed, n_back_mixed),
    S11(E, 0, 0, n_back_mixed, n_back_mixed);
  mat::set_identity(g_y_back);
  g_y_back_tmp = Z11_prime;
  LU2.invMult("N", g_y_back_tmp, g_y_back);
  g_y_back_tmp = g_y_back;
  blas::gemm("N", "N", 1.0, S11, g_y_back_tmp, 0.0, g_y_back);
  LU2.invMult("N", T11, g_y_back);
  g_y_back_tmp = g_y_back;
  blas::gemm("N", "N", 1.0, Z11_prime, g_y_back_tmp, 0.0, g_y_back);

  // TODO: avoid to copy mixed variables again, rather test it...
  for (size_t i = 0; i < n_back_mixed; i++)
    mat::row_copy(g_y_back, i, g_y, zeta_back_mixed[i]);

  // Compute DR for static variables w.r. to endogenous
  g_y_static = MatrixView(A, 0, 0, n_static, n_back_mixed);
  for (size_t i = 0; i < n_dynamic; i++)
    {
      mat::row_copy(g_y, zeta_dynamic[i], g_y_dynamic, i);
      mat::col_copy(A, n_back_mixed + zeta_dynamic[i], 0, n_static, A0d, i, 0);
    }
  blas::gemm("N", "N", 1.0, A0d, g_y_dynamic, 1.0, g_y_static);
  blas::gemm("N", "N", 1.0, g_y_fwrd, g_y_back, 0.0, g_y_static_tmp);
  blas::gemm("N", "N", 1.0, MatrixView(A, 0, n_back_mixed + n, n_static, n_fwrd_mixed),
             g_y_static_tmp, 1.0, g_y_static);
  for (size_t i = 0; i < n_static; i++)
    mat::col_copy(A, n_back_mixed + zeta_static[i], 0, n_static, A0s, i, 0);
  LU3.invMult("N", A0s, g_y_static);
  mat::negate(g_y_static);

  for (size_t i = 0; i < n_static; i++)
    mat::row_copy(g_y_static, i, g_y, zeta_static[i]);

  // Compute DR for all endogenous w.r. to shocks
  blas::gemm("N", "N", 1.0, MatrixConstView(jacobian, 0, n_back_mixed + n, n, n_fwrd_mixed), g_y_fwrd, 0.0, g_u_tmp1);
  g_u_tmp2 = MatrixConstView(jacobian, 0, n_back_mixed, n, n);
  for (size_t i = 0; i < n_back_mixed; i++)
    {
      VectorView c1 = mat::get_col(g_u_tmp2, zeta_back_mixed[i]),
        c2 = mat::get_col(g_u_tmp1, i);
      vec::add(c1, c2);
    }
  g_u = MatrixConstView(jacobian, 0, n_back_mixed + n + n_fwrd_mixed, n, p);
  LU4.invMult("N", g_u_tmp2, g_u);
  mat::negate(g_u);
}

std::ostream &
operator<<(std::ostream &out, const DecisionRules::BlanchardKahnException &e)
{
  if (e.order)
    out << "The Blanchard-Kahn order condition is not satisfied: you have " << e.n_fwrd_vars << " forward variables for " << e.n_explosive_eigenvals << " explosive eigenvalues";
  else
    out << "The Blanchard Kahn rank condition is not satisfied";
  return out;
}

