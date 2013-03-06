/*
 * Copyright (C) 2010-2013 Dynare Team
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
  InitializeKalmanFilter(const std::string &basename, size_t n_endo, size_t n_exo, const std::vector<size_t> &zeta_fwrd_arg,
                         const std::vector<size_t> &zeta_back_arg, const std::vector<size_t> &zeta_mixed_arg, const std::vector<size_t> &zeta_static_arg,
                         const std::vector<size_t> &zeta_varobs_back_mixed_arg,
                         const std::vector<size_t> &varobs_arg,
                         double qz_criterium_arg, double lyapunov_tol_arg,
                         bool noconstant_arg);
  virtual ~InitializeKalmanFilter();
  // initialise parameter dependent KF matrices only but not Ps
  template <class Vec1, class Vec2, class Mat1, class Mat2>
  void initialize(Vec1 &steadyState, const Vec2 &deepParams, Mat1 &R,
				     const Mat2 &Q, Matrix &RQRt, Matrix &T,
				     const MatrixConstView &dataView,
				     MatrixView &detrendedDataView)
  {
    modelSolution.compute(steadyState, deepParams, g_x, g_u);
    detrendData.detrend(steadyState, dataView, detrendedDataView);

    setT(T);
    setRQR(R, Q, RQRt);
  }

  // initialise all KF matrices
  template <class Vec1, class Vec2, class Mat1, class Mat2>
  void initialize(Vec1 &steadyState, const Vec2 &deepParams, Mat1 &R,
				     const Mat2 &Q, Matrix &RQRt, Matrix &T, Matrix &Pstar, Matrix &Pinf,
				     const MatrixConstView &dataView,
				     MatrixView &detrendedDataView)
  {
    initialize(steadyState, deepParams, R, Q, RQRt, T, dataView, detrendedDataView);
    setPstar(Pstar, Pinf, T, RQRt);
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
  void setT(Matrix &T);

  template <class Mat1, class Mat2>
  void setRQR(Mat1 &R, const Mat2 &Q, Matrix &RQRt)
  {
    mat::assignByVectors(R, mat::nullVec, mat::nullVec, g_u, zeta_varobs_back_mixed, mat::nullVec);

    //  Matrix RQRt=R*Q*R'
    blas::gemm("N", "N", 1.0, R, Q, 0.0, RQ); // R*Q
    blas::gemm("N", "T", 1.0, RQ, R, 0.0, RQRt); // R*Q*R'
  }
  void setPstar(Matrix &Pstar, Matrix &Pinf, const Matrix &T, const Matrix &RQRt) throw (DiscLyapFast::DLPException);

};

#endif // !defined(C3D996B8_22AB_4b77_B693_BA4777AFB091__INCLUDED_)
