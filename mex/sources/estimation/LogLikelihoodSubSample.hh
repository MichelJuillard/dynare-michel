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
//  LogLikelihoodSubSample.h
//  Implementation of the Class LogLikelihoodSubSample
//  Created on:      14-Jan-2010 22:39:14
///////////////////////////////////////////////////////////

#if !defined(DF8B7AF5_8169_4587_9037_2CD2C82E2DDF__INCLUDED_)
#define DF8B7AF5_8169_4587_9037_2CD2C82E2DDF__INCLUDED_

#include "KalmanFilter.hh"

class LogLikelihoodSubSample {

public:
  LogLikelihoodSubSample(const std::string &modName,  size_t n_endo, size_t n_exo,
                         const std::vector<size_t> &zeta_fwrd_arg, const std::vector<size_t> &zeta_back_arg, const std::vector<size_t> &zeta_mixed_arg,
                         const std::vector<size_t> &zeta_static_arg, const Matrix &ll_incidence, const double qz_criterium,  const std::vector<size_t> &order_var,
                         const std::vector<size_t> &inv_order_var, const std::vector<size_t> &varobs, const std::vector<size_t> &riv,
                         const std::vector<size_t> &ric, double riccati_tol_in, double lyapunov_tol, int &info);

  double compute(Vector &steadyState, const MatrixView &dataView, const Vector &deepParams, //const estPeriod &estPeriod,
                 VectorView &vll, int &info,  size_t start, size_t period, const Matrix &Q, const Matrix &H);
  virtual ~LogLikelihoodSubSample();

public:
  KalmanFilter kalmanFilter;

private:
  double penalty;
  double  logLikelihood;

};

#endif // !defined(DF8B7AF5_8169_4587_9037_2CD2C82E2DDF__INCLUDED_)
