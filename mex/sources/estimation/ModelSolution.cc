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
//  ModelSolution.cpp
//  Implementation of the Class ModelSolution
//  Created on:      02-Feb-2010 13:06:35
///////////////////////////////////////////////////////////

#include <string>

#include "ModelSolution.hh"

/**
 * compute the steady state (2nd stage), and computes first order approximation
 */
ModelSolution::ModelSolution(const std::string &dynamicDllFile,  size_t n_endo_arg, size_t n_exo_arg, const std::vector<size_t> &zeta_fwrd_arg,
                             const std::vector<size_t> &zeta_back_arg, const std::vector<size_t> &zeta_mixed_arg,
                             const std::vector<size_t> &zeta_static_arg, double INqz_criterium) :
  n_endo(n_endo_arg), n_exo(n_exo_arg),  // n_jcols = Num of Jacobian columns = nStat+2*nPred+3*nBoth+2*nForw+nExog
  n_jcols(n_exo+n_endo+ zeta_back_arg.size() /*nsPred*/ + zeta_fwrd_arg.size() /*nsForw*/ +2*zeta_mixed_arg.size()),
  jacobian(n_endo, n_jcols), residual(n_endo), Mx(1, n_exo),
  decisionRules(n_endo_arg, n_exo_arg, zeta_fwrd_arg, zeta_back_arg, zeta_mixed_arg, zeta_static_arg, INqz_criterium),
  dynamicDLLp(dynamicDllFile, n_exo),
  llXsteadyState(n_jcols-n_exo)
{
  Mx.setAll(0.0);
  jacobian.setAll(0.0);

  set_union(zeta_fwrd_arg.begin(), zeta_fwrd_arg.end(),
            zeta_mixed_arg.begin(), zeta_mixed_arg.end(),
            back_inserter(zeta_fwrd_mixed));
  set_union(zeta_back_arg.begin(), zeta_back_arg.end(),
            zeta_mixed_arg.begin(), zeta_mixed_arg.end(),
            back_inserter(zeta_back_mixed));
}


