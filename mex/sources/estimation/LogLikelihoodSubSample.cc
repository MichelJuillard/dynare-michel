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
//  LogLikelihoodSubSample.cpp
//  Implementation of the Class LogLikelihoodSubSample
//  Created on:      14-Jan-2010 22:39:14
///////////////////////////////////////////////////////////

//#include "LogLikelihoodSubSample.hh"
#include "LogLikelihoodMain.hh" // use ...Main.hh for testing only

LogLikelihoodSubSample::~LogLikelihoodSubSample()
{
};

LogLikelihoodSubSample::LogLikelihoodSubSample(const std::string &dynamicDllFile, EstimatedParametersDescription &INestiParDesc, size_t n_endo, size_t n_exo,
                                               const std::vector<size_t> &zeta_fwrd_arg, const std::vector<size_t> &zeta_back_arg,
                                               const std::vector<size_t> &zeta_mixed_arg, const std::vector<size_t> &zeta_static_arg, const double qz_criterium,
                                               const std::vector<size_t> &varobs, double riccati_tol, double lyapunov_tol, int &INinfo) :
  startPenalty(-1e8), estiParDesc(INestiParDesc),
  kalmanFilter(dynamicDllFile, n_endo, n_exo, zeta_fwrd_arg, zeta_back_arg, zeta_mixed_arg, zeta_static_arg, qz_criterium,
               varobs, riccati_tol, lyapunov_tol, INinfo), eigQ(n_exo), eigH(varobs.size()), info(INinfo)
{
};

