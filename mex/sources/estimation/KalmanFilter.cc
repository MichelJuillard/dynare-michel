/*
 * Copyright (C) 2009-2012 Dynare Team
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
#include "LapackBindings.hh"

KalmanFilter::~KalmanFilter()
{

}

KalmanFilter::KalmanFilter(const std::string &dynamicDllFile, size_t n_endo, size_t n_exo,
                           const std::vector<size_t> &zeta_fwrd_arg, const std::vector<size_t> &zeta_back_arg,
                           const std::vector<size_t> &zeta_mixed_arg, const std::vector<size_t> &zeta_static_arg,
                           double qz_criterium_arg, const std::vector<size_t> &varobs_arg,
                           double riccati_tol_arg, double lyapunov_tol_arg, int &info) :
  zeta_varobs_back_mixed(compute_zeta_varobs_back_mixed(zeta_back_arg, zeta_mixed_arg, varobs_arg)),
  Z(varobs_arg.size(), zeta_varobs_back_mixed.size()), Zt(Z.getCols(), Z.getRows()), T(zeta_varobs_back_mixed.size()), R(zeta_varobs_back_mixed.size(), n_exo),
  Pstar(zeta_varobs_back_mixed.size(), zeta_varobs_back_mixed.size()), Pinf(zeta_varobs_back_mixed.size(), zeta_varobs_back_mixed.size()),
  RQRt(zeta_varobs_back_mixed.size(), zeta_varobs_back_mixed.size()), Ptmp(zeta_varobs_back_mixed.size(), zeta_varobs_back_mixed.size()), F(varobs_arg.size(), varobs_arg.size()),
  Finv(varobs_arg.size(), varobs_arg.size()), K(zeta_varobs_back_mixed.size(), varobs_arg.size()), KFinv(zeta_varobs_back_mixed.size(), varobs_arg.size()),
  oldKFinv(zeta_varobs_back_mixed.size(), varobs_arg.size()), a_init(zeta_varobs_back_mixed.size()),
  a_new(zeta_varobs_back_mixed.size()), vt(varobs_arg.size()), vtFinv(varobs_arg.size()), riccati_tol(riccati_tol_arg),
  initKalmanFilter(dynamicDllFile, n_endo, n_exo, zeta_fwrd_arg, zeta_back_arg, zeta_mixed_arg,
                   zeta_static_arg, zeta_varobs_back_mixed, qz_criterium_arg, lyapunov_tol_arg, info),
  FUTP(varobs_arg.size()*(varobs_arg.size()+1)/2)
{
  Z.setAll(0.0);
  Zt.setAll(0.0);
  for (size_t i = 0; i < varobs_arg.size(); ++i)
    {
      size_t j = find(zeta_varobs_back_mixed.begin(), zeta_varobs_back_mixed.end(),
                      varobs_arg[i]) - zeta_varobs_back_mixed.begin();
      Z(i, j) = 1.0;
      Zt(j, i) = 1.0;
    }
}

std::vector<size_t>
KalmanFilter::compute_zeta_varobs_back_mixed(const std::vector<size_t> &zeta_back_arg, const std::vector<size_t> &zeta_mixed_arg, const std::vector<size_t> &varobs_arg)
{
  std::vector<size_t> varobs_sorted = varobs_arg;
  sort(varobs_sorted.begin(), varobs_sorted.end());
  std::vector<size_t> zeta_back_mixed, zeta_varobs_back_mixed;
  set_union(zeta_back_arg.begin(), zeta_back_arg.end(),
            zeta_mixed_arg.begin(), zeta_mixed_arg.end(),
            back_inserter(zeta_back_mixed));
  set_union(zeta_back_mixed.begin(), zeta_back_mixed.end(),
            varobs_sorted.begin(), varobs_sorted.end(),
            back_inserter(zeta_varobs_back_mixed));
  return zeta_varobs_back_mixed;
}


/**
 * Multi-variate standard Kalman Filter
 */
double
KalmanFilter::filter(const MatrixView &detrendedDataView,  const Matrix &H, VectorView &vll, size_t start, int &info)
{
  double loglik = 0.0, ll, logFdet = 0.0, Fdet, dvtFinvVt;
  size_t p = Finv.getRows();
  bool nonstationary = true;
  a_init.setAll(0.0);
  for (size_t t = 0; t < detrendedDataView.getCols(); ++t)
    {
      if (nonstationary)
        {
          // K=PZ'
          //blas::gemm("N", "T", 1.0, Pstar, Z, 0.0, K);
          blas::symm("L", "U", 1.0, Pstar, Zt, 0.0, K);

          //F=ZPZ' +H = ZK+H
          F = H;
          blas::gemm("N", "N", 1.0, Z, K, 1.0, F);
          // logFdet=log|F|

          // Finv=inv(F)
          mat::set_identity(Finv);
          // Pack F upper trinagle as vector
          for (size_t i = 1; i <= p; ++i)
            for (size_t j = i; j <= p; ++j)
              FUTP(i + (j-1)*j/2 -1) = F(i-1, j-1);

          info = lapack::choleskySolver(FUTP, Finv, "U"); // F now contains its Chol decomposition!
          if (info < 0)
            {
              std::cout << "Info:" << info << std::endl;
              std::cout << "t:" << t << std::endl;
              return 0;
            }
          if (info > 0)
            {
              //enforce Pstar symmetry with P=(P+P')/2=0.5P+0.5P'
              mat::transpose(Ptmp, Pstar);
              mat::add(Pstar, Ptmp);
              for (size_t i = 0; i < Pstar.getCols(); ++i)
                for (size_t j = 0; j < Pstar.getCols(); ++j)
                  Pstar(i, j) *= 0.5;

              // K=PZ'
              //blas::gemm("N", "T", 1.0, Pstar, Z, 0.0, K);
              blas::symm("L", "U", 1.0, Pstar, Zt, 0.0, K);

              //F=ZPZ' +H = ZK+H
              F = H;
              blas::gemm("N", "N", 1.0, Z, K, 1.0, F);

              // Finv=inv(F)
              mat::set_identity(Finv);
              // Pack F upper trinagle as vector
              for (size_t i = 1; i <= p; ++i)
                for (size_t j = i; j <= p; ++j)
                  FUTP(i + (j-1)*j/2 -1) = F(i-1, j-1);

              info = lapack::choleskySolver(FUTP, Finv, "U"); // F now contains its Chol decomposition!
              if (info != 0)
                {
                  return 0;
                }
            }
          // KFinv gain matrix
          blas::symm("R", "U", 1.0, Finv, K, 0.0, KFinv);
          // deteminant of F:
          Fdet = 1;
          for (size_t d = 1; d <= p; ++d)
            Fdet *= FUTP(d + (d-1)*d/2 -1);
          Fdet *= Fdet;

          logFdet = log(fabs(Fdet));

          Ptmp = Pstar;
          // Pt+1= T(Pt - KFinvK')T' +RQR'
          // 1) Ptmp= Pt - K*FinvK'
          blas::gemm("N", "T", -1.0, KFinv, K, 1.0, Ptmp);
          // 2) Ptmp= T*Ptmp
          Pstar = Ptmp;
          //blas::gemm("N", "N", 1.0, T, Pstar, 0.0, Ptmp);
          blas::symm("R", "U", 1.0, Pstar, T, 0.0, Ptmp);
          // 3) Pt+1= Ptmp*T' +RQR'
          Pstar = RQRt;
          blas::gemm("N", "T", 1.0, Ptmp, T, 1.0, Pstar);

          if (t > 0)
            nonstationary = mat::isDiff(KFinv, oldKFinv, riccati_tol);
          oldKFinv = KFinv;
        }

      // err= Yt - Za
      VectorConstView yt = mat::get_col(detrendedDataView, t);
      vt = yt;
      blas::gemv("N", -1.0, Z, a_init, 1.0, vt);

      // at+1= T(at+ KFinv *err)
      blas::gemv("N", 1.0, KFinv, vt, 1.0, a_init);
      blas::gemv("N", 1.0, T, a_init, 0.0, a_new);
      a_init = a_new;

      /*****************
         Here we calc likelihood and store results.
      *****************/
      blas::symv("U", 1.0, Finv, vt, 0.0, vtFinv);
      dvtFinvVt = blas::dot(vtFinv, vt);

      ll = -0.5*(p*log(2*M_PI)+logFdet+dvtFinvVt);

      vll(t) = ll;
      if (t >= start)
        loglik += ll;

    }

  return loglik;
}

