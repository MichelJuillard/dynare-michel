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
//  LogPosteriorDensity.hh
//  Implementation of the Class LogPosteriorDensity
//  Created on:      10-Feb-2010 20:54:18
///////////////////////////////////////////////////////////

#if !defined(LPD_052A31B5_53BF_4904_AD80_863B52827973__INCLUDED_)
#define LPD_052A31B5_53BF_4904_AD80_863B52827973__INCLUDED_

#include "EstimatedParametersDescription.hh"
#include "LogPriorDensity.hh"
#include "LogLikelihoodMain.hh"

/**
 * Class that calculates Log Posterior Density using kalman, based on Dynare
 * DsgeLikelihood.m
 */
class LogPosteriorDensity
{

private:
  LogPriorDensity logPriorDensity;
  LogLikelihoodMain logLikelihoodMain;

public:
  virtual ~LogPosteriorDensity();

  LogPosteriorDensity(const std::string &modName, EstimatedParametersDescription &estParamsDesc, size_t n_endo, size_t n_exo,
                      const std::vector<size_t> &zeta_fwrd_arg, const std::vector<size_t> &zeta_back_arg, const std::vector<size_t> &zeta_mixed_arg,
                      const std::vector<size_t> &zeta_static_arg, const double qz_criterium_arg, const std::vector<size_t> &varobs_arg,
                      double riccati_tol_arg, double lyapunov_tol_arg, int &info_arg);

  template <class VEC1, class VEC2>
  double
  compute(VEC1 &steadyState, VEC2 &estParams, VectorView &deepParams, const MatrixConstView &data, MatrixView &Q, Matrix &H, size_t presampleStart, int &info)
  {
    return -logLikelihoodMain.compute(steadyState, estParams, deepParams, data, Q, H, presampleStart, info)
      -logPriorDensity.compute(estParams);
  }

  Vector&getLikVector();

};

#endif // !defined(052A31B5_53BF_4904_AD80_863B52827973__INCLUDED_)
