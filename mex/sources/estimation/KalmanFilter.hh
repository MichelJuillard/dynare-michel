/*
 * Copyright (C) 2009-2013 Dynare Team
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
//  KalmanFilter.h
//  Implementation of the Class KalmanFilter
//  Created on:      02-Feb-2010 12:44:41
///////////////////////////////////////////////////////////

#if !defined(KF_213B0417_532B_4027_9EDF_36C004CB4CD1__INCLUDED_)
#define KF_213B0417_532B_4027_9EDF_36C004CB4CD1__INCLUDED_

#include "InitializeKalmanFilter.hh"

/**
 * Vanilla Kalman filter without constant and with measurement error (use scalar
 * 0 when no measurement error).
 * If multivariate filter is faster, do as in Matlab: start with multivariate
 * filter and switch to univariate filter only in case of singularity
 *
 * mamber functions: compute() and filter()
 * OUTPUT
 *    LIK:    likelihood
 *
 * REFERENCES
 *   See "Filtering and Smoothing of State Vector for Diffuse State Space
 *   Models", S.J. Koopman and J. Durbin (2003, in Journal of Time Series
 *   Analysis, vol. 24(1), pp. 85-98).
 */

class KalmanFilter
{

public:
  virtual ~KalmanFilter();
  KalmanFilter(const std::string &basename, size_t n_endo, size_t n_exo, const std::vector<size_t> &zeta_fwrd_arg,
               const std::vector<size_t> &zeta_back_arg, const std::vector<size_t> &zeta_mixed_arg, const std::vector<size_t> &zeta_static_arg,
               double qz_criterium_arg, const std::vector<size_t> &varobs_arg,
               double riccati_tol_arg, double lyapunov_tol_arg,
               bool noconstant_arg);

  template <class Vec1, class Vec2, class Mat1>
  double compute(const MatrixConstView &dataView, Vec1 &steadyState,
                 const Mat1 &Q, const Matrix &H, const Vec2 &deepParams,
                 VectorView &vll, MatrixView &detrendedDataView, size_t start, size_t period)
  {
	if (period == 0) // initialise all KF matrices
	  initKalmanFilter.initialize(steadyState, deepParams, R, Q, RQRt, T, Pstar, Pinf,
                                dataView, detrendedDataView);
	else             // initialise parameter dependent KF matrices only but not Ps
	  initKalmanFilter.initialize(steadyState, deepParams, R, Q, RQRt, T,
                                dataView, detrendedDataView);

  return filter(detrendedDataView, H, vll, start);
  }

private:
  const std::vector<size_t> zeta_varobs_back_mixed;
  static std::vector<size_t> compute_zeta_varobs_back_mixed(const std::vector<size_t> &zeta_back_arg, const std::vector<size_t> &zeta_mixed_arg, const std::vector<size_t> &varobs_arg);
  Matrix Z, Zt;   //nob*mm matrix mapping endogeneous variables and observations and its transpose
  Matrix T;   //mm*mm transition matrix of the state equation.
  Matrix R;   //mm*rr matrix, mapping structural innovations to state variables.
  Matrix Pstar; //mm*mm variance-covariance matrix of stationary variables
  Matrix Pinf;  //mm*mm variance-covariance matrix of diffuse variables
  // allocate space for intermediary matrices
  Matrix RQRt, Ptmp;  //mm*mm variance-covariance matrix of variable disturbances
  Matrix F, Finv;  // nob*nob F=ZPZt +H an inv(F)
  Matrix K,  KFinv, oldKFinv; // mm*nobs K=PZt and K*Finv gain matrices
  Vector a_init, a_new; // state vector
  Vector vt; // current observation error vectors
  Vector vtFinv; // intermediate observation error *Finv vector
  double riccati_tol;
  InitializeKalmanFilter initKalmanFilter; //Initialise KF matrices
  Vector FUTP; // F upper triangle packed as vector FUTP(i + (j-1)*j/2) = F(i,j) for 1<=i<=j;

  // Method
  double filter(const MatrixView &detrendedDataView,  const Matrix &H, VectorView &vll, size_t start);

};

#endif // !defined(213B0417_532B_4027_9EDF_36C004CB4CD1__INCLUDED_)
