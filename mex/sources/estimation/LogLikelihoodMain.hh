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
//  LogLikelihoodMain.h
//  Implementation of the Class LogLikelihoodMain
//  Created on:      02-Feb-2010 12:57:09
///////////////////////////////////////////////////////////

#if !defined(E126AEF5_AC28_400a_821A_3BCFD1BC4C22__INCLUDED_)
#define E126AEF5_AC28_400a_821A_3BCFD1BC4C22__INCLUDED_

//#include "EstimatedParametersDescription.hh"
#include "LogLikelihoodSubSample.hh"

class LogLikelihoodMain {
private:
  double logLikelihood;
  std::vector<EstimationSubsample> &estSubsamples; // reference to member of EstimatedParametersDescription
  LogLikelihoodSubSample logLikelihoodSubSample;
  Vector deepParams;
  Vector &vll;  // vector of all KF step likelihoods
  Matrix &data; // input data
  Matrix &steadyState;
  //GeneralParams& estimOptions;
  int presampleStart;
  Matrix Q;
  Matrix H;

public:
  virtual ~LogLikelihoodMain();
  LogLikelihoodMain(const Matrix &data, //GeneralParams& estimOptions,
                    const std::string &modName, EstimatedParametersDescription &estiParDesc, size_t n_endo, size_t n_exo,
                    const std::vector<size_t> &zeta_fwrd_arg, const std::vector<size_t> &zeta_back_arg, const std::vector<size_t> &zeta_mixed_arg,
                    const std::vector<size_t> &zeta_static_arg, const Matrix &ll_incidence, const double qz_criterium,  const std::vector<size_t> &order_var,
                    const std::vector<size_t> &inv_order_var, const std::vector<size_t> &varobs, const std::vector<size_t> &riv,
                    const std::vector<size_t> &ric, double riccati_tol_in, double lyapunov_tol, int &info);

  double compute(Matrix &steadyState, Vector &estParams, int &info); // for calls from estimation and to set Steady State
  double compute(Vector &estParams);  // for calls in loop from optimiser

  Vector &
  getVll()
  {
    return vll;
  };
  double
  getLogLikelihood()
  {
    return logLikelihood;
  };

};

#endif // !defined(E126AEF5_AC28_400a_821A_3BCFD1BC4C22__INCLUDED_)
