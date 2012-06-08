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
//  LogLikelihoodMain.h
//  Implementation of the Class LogLikelihoodMain
//  Created on:      02-Feb-2010 12:57:09
///////////////////////////////////////////////////////////

#if !defined(E126AEF5_AC28_400a_821A_3BCFD1BC4C22__INCLUDED_)
#define E126AEF5_AC28_400a_821A_3BCFD1BC4C22__INCLUDED_

#include "LogLikelihoodSubSample.hh"

class LogLikelihoodMain
{
private:
  std::vector<EstimationSubsample> &estSubsamples; // reference to member of EstimatedParametersDescription
  LogLikelihoodSubSample logLikelihoodSubSample;
  Vector vll;  // vector of all KF step likelihoods
  Matrix detrendedData;

public:
  virtual ~LogLikelihoodMain();
  LogLikelihoodMain(const std::string &dynamicDllFile, EstimatedParametersDescription &estiParDesc, size_t n_endo, size_t n_exo,
                    const std::vector<size_t> &zeta_fwrd_arg, const std::vector<size_t> &zeta_back_arg, const std::vector<size_t> &zeta_mixed_arg,
                    const std::vector<size_t> &zeta_static_arg, const double qz_criterium_arg, const std::vector<size_t> &varobs_arg,
                    double riccati_tol_arg, double lyapunov_tol_arg, int &info);

  /**
   * Compute method Inputs:
   * Matrix &steadyState; Matrix of sub-sample periods column-vectors of steady states, one column vectro for each sub-sample period
   * vectors of deep deepParams and estimated estParams
   * Matrix &data input data reference
   * Q and H KF matrices of shock and measurement error varinaces and covariances
   * KF logLikelihood calculation start period.
   */

  template <class VEC1, class VEC2>
  double compute(VEC1 &steadyState, VEC2 &estParams, VectorView &deepParams, const MatrixConstView &data, 
		 MatrixView &Q, Matrix &H, size_t start, int &info)
  {
    double logLikelihood = 0;
    for (size_t i = 0; i < estSubsamples.size(); ++i)
      {
	VectorView vSteadyState (steadyState,0,steadyState.getSize());

	MatrixConstView dataView(data, 0, estSubsamples[i].startPeriod,
				 data.getRows(), estSubsamples[i].endPeriod-estSubsamples[i].startPeriod+1);
	MatrixView detrendedDataView(detrendedData, 0, estSubsamples[i].startPeriod,
				     data.getRows(), estSubsamples[i].endPeriod-estSubsamples[i].startPeriod+1);

	VectorView vllView(vll, estSubsamples[i].startPeriod, estSubsamples[i].endPeriod-estSubsamples[i].startPeriod+1);
	logLikelihood += logLikelihoodSubSample.compute(vSteadyState, dataView, estParams, deepParams,
							Q, H, vllView, detrendedDataView, info, start, i);
      }
    return logLikelihood;
  };

  Vector &getVll() { return vll; };
};

#endif // !defined(E126AEF5_AC28_400a_821A_3BCFD1BC4C22__INCLUDED_)
