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
//  ModelSolution.cpp
//  Implementation of the Class ModelSolution
//  Created on:      02-Feb-2010 13:06:35
///////////////////////////////////////////////////////////

#include <string>

#include "ModelSolution.hh"

/**
* compute the steady state (2nd stage), and computes first order approximation
*/
ModelSolution::ModelSolution(const std::string& modName,  size_t n_endo_arg, size_t n_exo_arg, const std::vector<size_t>& zeta_fwrd_arg, 
                             const std::vector<size_t>& zeta_back_arg, const std::vector<size_t>& zeta_mixed_arg, 
                             const std::vector<size_t>& zeta_static_arg, const Matrix& llincidence, double INqz_criterium)
                             : n_endo(n_endo_arg), n_exo(n_exo_arg),  // n_jcols = Num of Jacobian columns = nStat+2*nPred+3*nBoth+2*nForw+nExog
                             n_jcols (n_exo+n_endo+ zeta_back_arg.size() /*nsPred*/ + zeta_fwrd_arg.size() /*nsForw*/ +2*zeta_mixed_arg.size()), 
                             ll_incidence(llincidence), jacobian (n_endo,n_jcols), residual(n_endo), Mx(2,n_exo),
                               decisionRules ( n_endo_arg, n_exo_arg, zeta_fwrd_arg, zeta_back_arg, zeta_mixed_arg, zeta_static_arg, INqz_criterium),
                               dynamicDLLp(modName, n_endo,  n_jcols,  /* nMax_lag= */ 1,  n_exo, std::string(""))
{
  Mx.setAll(0.0);
  jacobian.setAll(0.0);
}

void 
ModelSolution::compute(Vector& steadyState, const Vector& deepParams, Matrix& ghx, Matrix& ghu) throw (DecisionRules::BlanchardKahnException, GeneralizedSchurDecomposition::GSDException)
{
  // compute Steady State
  ComputeSteadyState(steadyState, deepParams);

  // then get jacobian and 

  ComputeModelSolution( steadyState,  deepParams, ghx, ghu);

}

void 
ModelSolution::ComputeModelSolution(Vector& steadyState, const Vector& deepParams, Matrix& ghx, Matrix& ghu) throw (DecisionRules::BlanchardKahnException, GeneralizedSchurDecomposition::GSDException)
{
  // set extended Steady State

  Vector llXsteadyState(n_jcols-n_exo);

  for (int ll_row = 0; ll_row < ll_incidence.getRows(); ll_row++)
    {
      // populate (non-sparse) vector with ysteady values
      for (int i = 0; i < n_endo; i++)
      {
        if (ll_incidence(ll_row, i))
          llXsteadyState(((int) ll_incidence(ll_row, i))-1) = steadyState(i);
      }
    }
#ifdef DEBUG
    //    std::cout << "Vector llXsteadyState: " << std::endl << llXsteadyState << std::endl;
    mexPrintf(" get jacobian \n");
#endif
    //get jacobian 
    dynamicDLLp.eval(llXsteadyState, Mx, &deepParams,  1,  residual, &jacobian, NULL, NULL);

    std::cout << "jacobian: " << std::endl << jacobian << std::endl;
#ifdef DEBUG
    mexPrintf(" compute rules \n");
#endif
    //compute rules
    decisionRules.compute(jacobian,ghx, ghu);
}
void 
ModelSolution::ComputeSteadyState(Vector& steadyState, const Vector& deepParams)
{
  // does nothig for time being.
}

