/*
 * Copyright (C) 2010-2012 Dynare Team
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
//  InitializeKalmanFilter.h
//  Implementation of the Class InitializeKalmanFilter
//  Created on:      02-Feb-2010 12:25:28
///////////////////////////////////////////////////////////

#if !defined(C3D996B8_22AB_4b77_B693_BA4777AFB091__INCLUDED_)
#define C3D996B8_22AB_4b77_B693_BA4777AFB091__INCLUDED_

#include "DetrendData.hh"
#include "ModelSolution.hh"
#include "DiscLyapFast.hh"
#include <string>

/**
 * if model is declared stationary ?compute covariance matrix of endogenous
 * variables () by doubling algorithm
 *
 */
class InitializeKalmanFilter
{

public:
  /*!
    \param[in] zeta_varobs_back_mixed_arg The union of indices of observed, backward and mixed variables
  */
  InitializeKalmanFilter(const std::string &dynamicDllFile, size_t n_endo, size_t n_exo, const std::vector<size_t> &zeta_fwrd_arg,
                         const std::vector<size_t> &zeta_back_arg, const std::vector<size_t> &zeta_mixed_arg, const std::vector<size_t> &zeta_static_arg,
                         const std::vector<size_t> &zeta_varobs_back_mixed_arg,
                         double qz_criterium_arg, double lyapunov_tol_arg, int &info);
  virtual ~InitializeKalmanFilter();
  // initialise parameter dependent KF matrices only but not Ps
  template <class VEC>
  void initialize(VEC &steadyState, const VectorView &deepParams, Matrix &R,
				     const MatrixView &Q, Matrix &RQRt, Matrix &T,
				     double &penalty, const MatrixConstView &dataView,
				     MatrixView &detrendedDataView, int &info)
  {
    modelSolution.compute(steadyState, deepParams, g_x, g_u);
    detrendData.detrend(steadyState, dataView, detrendedDataView);

    setT(T, info);
    setRQR(R, Q, RQRt, info);
  }

  // initialise all KF matrices
  template <class VEC>
  void initialize(VEC &steadyState, const VectorView &deepParams, Matrix &R,
				     const MatrixView &Q, Matrix &RQRt, Matrix &T, Matrix &Pstar, Matrix &Pinf,
				     double &penalty, const MatrixConstView &dataView,
				     MatrixView &detrendedDataView, int &info)
  {
    initialize(steadyState, deepParams, R, Q, RQRt, T, penalty, dataView, detrendedDataView, info);
    setPstar(Pstar, Pinf, T, RQRt, info);
  }

private:
  const double lyapunov_tol;
  const std::vector<size_t> zeta_varobs_back_mixed;
  //! Indices of back+mixed zetas inside varobs+back+mixed zetas
  std::vector<size_t> pi_bm_vbm;

  DetrendData detrendData;
  ModelSolution modelSolution;
  DiscLyapFast discLyapFast; //Lyapunov solver
  Matrix g_x;
  Matrix g_u;
  Matrix Rt, RQ;
  void setT(Matrix &T, int &info);
  void setRQR(Matrix &R, const MatrixView &Q, Matrix &RQRt, int &info);
  void setPstar(Matrix &Pstar, Matrix &Pinf, const Matrix &T, const Matrix &RQRt, int &info);

};

#endif // !defined(C3D996B8_22AB_4b77_B693_BA4777AFB091__INCLUDED_)
