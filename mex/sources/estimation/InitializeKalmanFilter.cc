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

///////////////////////////////////////////////////////////
//  InitializeKalmanFilter.cpp
//  Implementation of the Class InitializeKalmanFilter
//  Created on:      02-Feb-2010 12:25:28
///////////////////////////////////////////////////////////

#include "InitializeKalmanFilter.hh"
#include "BlasBindings.hh"

InitializeKalmanFilter::~InitializeKalmanFilter()
{
}

InitializeKalmanFilter::InitializeKalmanFilter(const std::string &modName, size_t n_endo_arg, size_t n_exo_arg,
                                               const std::vector<size_t> &zeta_fwrd_arg, const std::vector<size_t> &zeta_back_arg,
                                               const std::vector<size_t> &zeta_mixed_arg, const std::vector<size_t> &zeta_static_arg,
                                               const Matrix &llincidence, double qz_criterium,  const std::vector<size_t> &order_var_arg, const std::vector<size_t> &inv_order_var_arg,
                                               const std::vector<size_t> &rivIN, const std::vector<size_t> &ricIN, double lyapunov_tolIN, int &info) :
  riv(rivIN), inv_ric(ricIN.size()), order_var(order_var_arg), lyapunov_tol(lyapunov_tolIN),
  detrendData(false), modelSolution(modName, n_endo_arg, n_exo_arg, zeta_fwrd_arg, zeta_back_arg,
                                    zeta_mixed_arg, zeta_static_arg, llincidence, qz_criterium), discLyapFast(riv.size()),
  ghx(n_endo_arg, zeta_back_arg.size() + zeta_mixed_arg.size()),
  ghx_raw(n_endo_arg, zeta_back_arg.size() + zeta_mixed_arg.size()),
  ghu(n_endo_arg, n_exo_arg),  ghu_raw(n_endo_arg, n_exo_arg),
  Rt(n_exo_arg, riv.size()), RQ(riv.size(), n_exo_arg) 
{
  size_t n_pred = ghx.getCols();
  size_t n_static = zeta_static_arg.size();
  size_t j = 0;
  for (size_t i = 0; i < n_endo_arg; ++i)
    if (inv_order_var_arg[i]-n_static < n_pred && inv_order_var_arg[i]-n_static >= 0)
      inv_ric[j++] = ricIN[inv_order_var_arg[i]-n_static];
}

void
InitializeKalmanFilter::initialize(Vector &steadyState, const Vector &deepParams, Matrix &R, const Matrix &Z,  
                                   const Matrix &Q, Matrix &RQRt, Matrix &T, Matrix &Pstar, Matrix &Pinf,
                                   double &penalty, const MatrixView &dataView, Matrix &Y, int &info)
{
  modelSolution.compute(steadyState, deepParams, ghx_raw, ghu_raw);
  detrendData.detrend(steadyState, dataView, Y);

  mat::reorderRowsByVectors(ghx, mat::nullVec, ghx_raw, order_var);
  mat::reorderRowsByVectors(ghu, mat::nullVec, ghu_raw, order_var);

  setT(T, info);
  setR(R, info);
  setPstar(Pstar, Pinf, T, R, Q, RQRt, info);
}

void
InitializeKalmanFilter::setT(Matrix &T, int &info)
{
  // here is the content of [T R]=[A,B] = kalman_transition_matrix(oo_.dr,iv,ic,aux,M_.exo_nbr);

  T.setAll(0.0);

  //T(i_n_iv,ic) = dr.ghx(iv,:);
  mat::assignByVectors(T, mat::nullVec, inv_ric, ghx, riv, mat::nullVec);
}

void
InitializeKalmanFilter::setR(Matrix &R, int &info)
{
  R.setAll(0.0);
  //B(i_n_iv,:) = dr.ghu(iv,:); and R=B;
  mat::assignByVectors(R, mat::nullVec, mat::nullVec, ghu, riv, mat::nullVec);
}

void
InitializeKalmanFilter::setPstar(Matrix &Pstar, Matrix &Pinf, Matrix &T, Matrix &R, const Matrix &Q, Matrix &RQRt, int &info)
{

  //  Matrix RQRt=R*Q*R'
  RQ.setAll(0.0);
  blas::gemm("N", "N", 1.0, R, Q, 0.0, RQ); // R*Q
  RQRt.setAll(0.0);
  //mat::transpose(Rt, R);
  blas::gemm("N", "T", 1.0, RQ, R, 0.0, RQRt); // R*Q*R'
  //mat::transpose(RQR);

  try
  {
    // disclyap_fast(T, RQR, Pstar, lyapunov_tol, 0 or 1 to check chol)
    discLyapFast.solve_lyap(T, RQRt, Pstar, lyapunov_tol, 0); 

    Pinf.setAll(0.0);
  }
  catch(const DiscLyapFast::DLPException &e)
  {
    if (e.info > 0) // The matrix is not positive definite in NormCholesky calculator
      {
        printf(e.message.c_str());
        info = -1; //likelihood = penalty;
        return;
      }
    else if (e.info < 0)
      {
        printf("Caugth unhandled TS exception with Pstar matrix: ");
        printf(e.message.c_str());
        info = -1; //likelihood = penalty;
        throw;
      }
  }
}

