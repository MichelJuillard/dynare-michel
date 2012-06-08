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
//  LogPosteriorDensity.cpp
//  Implementation of the Class LogPosteriorDensity
//  Created on:      10-Feb-2010 20:54:18
///////////////////////////////////////////////////////////

#include "LogPosteriorDensity.hh"

LogPosteriorDensity::~LogPosteriorDensity()
{
}

LogPosteriorDensity::LogPosteriorDensity(const std::string &modName, EstimatedParametersDescription &estParamsDesc, size_t n_endo, size_t n_exo,
                                         const std::vector<size_t> &zeta_fwrd_arg, const std::vector<size_t> &zeta_back_arg, const std::vector<size_t> &zeta_mixed_arg,
                                         const std::vector<size_t> &zeta_static_arg, const double qz_criterium_arg, const std::vector<size_t> &varobs_arg,
                                         double riccati_tol_arg, double lyapunov_tol_arg, int &info_arg) :
  logPriorDensity(estParamsDesc),
  logLikelihoodMain(modName, estParamsDesc, n_endo, n_exo, zeta_fwrd_arg, zeta_back_arg, zeta_mixed_arg,
                    zeta_static_arg, qz_criterium_arg, varobs_arg, riccati_tol_arg, lyapunov_tol_arg, info_arg)
{

}


/**
 * vector of log likelihoods for each Kalman step
 */
Vector &
LogPosteriorDensity::getLikVector()
{
  return logLikelihoodMain.getVll();
}

