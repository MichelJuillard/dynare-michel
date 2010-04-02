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
  InitializeKalmanFilter(const std::string& modName, size_t n_endo, size_t n_exo, const std::vector<size_t> &zeta_fwrd_arg,
                         const std::vector<size_t> &zeta_back_arg, const std::vector<size_t> &zeta_mixed_arg, const std::vector<size_t> &zeta_static_arg,
                         const Matrix &llincidence, double qz_criterium, const std::vector<size_t> &order_var_arg, const std::vector<size_t> &inv_order_var_arg,
                         const std::vector<size_t> &riv, const std::vector<size_t> &ric, double lyapunov_tol, int &info);
  virtual ~InitializeKalmanFilter();
  void initialize(Vector &steadyState, const Vector &deepParams, Matrix &R, const Matrix &Z, const Matrix &Q, Matrix &RQRt,
                  Matrix &T, Matrix &Pstar, Matrix &Pinf, double &penalty, const MatrixView &dataView, Matrix &Y, int &info);

private:
  const std::vector<size_t> riv; // restrict_var_list
  std::vector<size_t> inv_ric; // inverse restrict_columns
  const std::vector<size_t> order_var;
  const double lyapunov_tol;

  DetrendData detrendData;
  ModelSolution modelSolution;
  DiscLyapFast discLyapFast; //Lyapunov solver
  Matrix ghx, ghx_raw;
  Matrix ghu, ghu_raw;
  Matrix Rt, RQ;
  void setT(Matrix &T, int &info);
  void setR(Matrix &R, int &info);
  void setPstar(Matrix &Pstar, Matrix &Pinf, Matrix &T, Matrix &R,  const Matrix &Q, Matrix &RQRt, int &info);

};

#endif // !defined(C3D996B8_22AB_4b77_B693_BA4777AFB091__INCLUDED_)
