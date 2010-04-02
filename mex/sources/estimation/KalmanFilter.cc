/*
 * Copyright (C) 2009-2010 Dynare Team
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
//  KalmanFilter.cpp
//  Implementation of the Class KalmanFilter
//  Created on:      02-Feb-2010 12:44:41
///////////////////////////////////////////////////////////

#include "KalmanFilter.hh"

KalmanFilter::~KalmanFilter()
{

}

KalmanFilter::KalmanFilter(const std::string &modName, size_t n_endo, size_t n_exo,
                           const std::vector<size_t> &zeta_fwrd_arg, const std::vector<size_t> &zeta_back_arg,
                           const std::vector<size_t> &zeta_mixed_arg, const std::vector<size_t> &zeta_static_arg,
                           const Matrix &ll_incidence, double qz_criterium,  const std::vector<size_t> &order_var,
                           const std::vector<size_t> &inv_order_var, const std::vector<size_t> &varobs_arg,
                           const std::vector<size_t> &riv, const std::vector<size_t> &ric,
                           double riccati_tol_in, double lyapunov_tol, int &info) :
  Z(varobs_arg.size(), riv.size()), T(riv.size()), R(riv.size(), n_exo),
  Pstar(riv.size(), riv.size()), Pinf(riv.size(), riv.size()), RQRt(riv.size(), riv.size()),
  Pold(riv.size(), riv.size()), Ptmp(riv.size(), riv.size()), F(varobs_arg.size(), varobs_arg.size()),
  Finv(varobs_arg.size(), varobs_arg.size()), K(riv.size(), varobs_arg.size()),
  KFinv(riv.size(), varobs_arg.size()), Finverter(varobs_arg.size()), a_init(riv.size(), 1),
  a_new(riv.size(), 1), vt(varobs_arg.size(), 1), vtFinv(1, varobs_arg.size()), vtFinvVt(1), riccati_tol(riccati_tol_in),
  initKalmanFilter(modName, n_endo, n_exo, zeta_fwrd_arg, zeta_back_arg, zeta_mixed_arg,
                   zeta_static_arg, ll_incidence, qz_criterium, order_var, inv_order_var, riv, ric, lyapunov_tol, info)
{
  a_init.setAll(0.0);
  Z.setAll(0.0);
  for (size_t i = 0; i < varobs_arg.size(); ++i)
    Z(i, varobs_arg[i]) = 1.0;

}

double
KalmanFilter::compute(const MatrixView &dataView, Vector &steadyState,
                      const Matrix &Q, const Matrix &H, const Vector &deepParams,
                      Vector &vll, size_t start, double &penalty, int &info)
{
  double ll;
  Matrix Y(dataView.getRows(), dataView.getCols());    // data

  initKalmanFilter.initialize(steadyState, deepParams, R, Z, Q, RQRt, T, Pstar, Pinf,
                              penalty, dataView, Y, info);

  return ll = filter(Y, H, vll, start, info);

};

/**
 * 30:*
 */
double
KalmanFilter::filter(const Matrix &dataView,  const Matrix &H, Vector &vll, size_t start, int &info)
{
  double loglik, ll, logFdet;
  int p = Finv.getRows();

  bool nonstationary = true;
  for (size_t t = 0; t < dataView.getCols(); ++t)
    {
      if (nonstationary)
        {
          // K=PZ'
          K.setAll(0.0);
          blas::gemm("N", "T", 1.0, Pstar, Z, 1.0, K);

          //F=ZPZ' +H = ZK+H
          F = H;
          blas::gemm("N", "N", 1.0, Z, K, 1.0, F);
          // logFdet=log|F|

          // Finv=inv(F)
          mat::set_identity(Finv);
          Finverter.invMult("N", F, Finv); // F now contains its LU decomposition!
          // KFinv
          KFinv.setAll(0.0);
          blas::gemm("N", "N", 1.0, K, Finv, 1.0, KFinv);
          // deteminant of F:
          logFdet = 1;
          for (int d = 0; d < p; ++d)
            logFdet *= (-F(d, d));
          Pold = Pstar;
          Ptmp = Pstar;
          // Pt+1= T(Pt - KFinvK')T' +RQR'
          // 1) Ptmp= Pt - K*FinvK'
          blas::gemm("N", "T", -1.0, KFinv, K, 1.0, Ptmp);
          // 2) Ptmp= T*Ptmp
          Pstar = Ptmp;
          Ptmp.setAll(0.0);
          blas::gemm("N", "N", 1.0, T, Pstar, 1.0, Ptmp);
          // 3) Pt+1= Ptmp*T' +RQR'
          Pstar = RQRt;
          blas::gemm("N", "T", 1.0, Ptmp, T, 1.0, Pstar);
          nonstationary = mat::isDiffSym(Pstar, Pold, riccati_tol);
        }
      else
        {
        }

      // err= Yt - Za
      //VectorView yt(dataView.getData()+t*dataView.getRows(),dataView.getRows(),1); // current observation vector
      MatrixConstView yt(dataView, 0, t, dataView.getRows(), 1); // current observation vector
      vt = yt;
      blas::gemm("N", "N", -1.0, Z, a_init, 1.0, vt);
      // at+1= T(at+ KFinv *err)
      blas::gemm("N", "N", 1.0, KFinv, vt, 1.0, a_init);
      blas::gemm("N", "N", 1.0, T, a_init, 0.0, a_new);
      a_init = a_new;

      /*****************
         Here we calc likelihood and store results.
      *****************/
      blas::gemm("T", "N", 1.0, vt, Finv, 0.0, vtFinv);
      blas::gemm("N", "N", 1.0, vtFinv, vt, 0.0, vtFinvVt);

      ll = -0.5*(p*log(2*M_PI)+logFdet+*(vtFinvVt.getData()));

      vll(t) = ll;
      if (t > start) loglik += ll;

    }

  return loglik;
}

/**
 * 89:*
   double KalmanFilter::calcStepLogLik(const PLUFact &Finv, const Vector &v){

   }
 */

