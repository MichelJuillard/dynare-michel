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
//  ComputeModelSolution.h
//  Implementation of the Class ModelSolution
//  Created on:      15-Jan-2010 07:37:47
///////////////////////////////////////////////////////////

#if !defined(ModelSolution_5ADFF920_9C74_46f5_9FE9_88AD4D4BBF19__INCLUDED_)
#define ModelSolution_5ADFF920_9C74_46f5_9FE9_88AD4D4BBF19__INCLUDED_

#include "DecisionRules.hh"
#include "dynamic_dll.hh"

/**
 * compute the steady state (2nd stage), and
 * computes first order approximation
 *
 */
class ModelSolution
{

public:
  ModelSolution(const std::string &dynamicDllFile,  size_t n_endo, size_t n_exo, const std::vector<size_t> &zeta_fwrd_arg,
                const std::vector<size_t> &zeta_back_arg, const std::vector<size_t> &zeta_mixed_arg,
                const std::vector<size_t> &zeta_static_arg, double qz_criterium);
  virtual ~ModelSolution() {};
  template <class VEC>
  void compute(VEC &steadyState, const VectorView &deepParams,      Matrix &ghx, Matrix &ghu) throw (DecisionRules::BlanchardKahnException, GeneralizedSchurDecomposition::GSDException)
  {
    // compute Steady State
    ComputeSteadyState(steadyState, deepParams);

    // then get jacobian and

    ComputeModelSolution(steadyState, deepParams, ghx, ghu);

  }

private:
  const size_t n_endo;
  const size_t n_exo;
  const size_t n_jcols; // Num of Jacobian columns
  std::vector<size_t> zeta_fwrd_mixed, zeta_back_mixed;
  Matrix jacobian;
  Vector residual;
  Matrix Mx;
  DecisionRules decisionRules;
  DynamicModelDLL dynamicDLLp;
  Vector llXsteadyState;
  //Matrix jacobian;
  template <class VEC>
  void ComputeModelSolution(VEC &steadyState, const VectorView &deepParams,         
			    Matrix &ghx, Matrix &ghu) 
    throw (DecisionRules::BlanchardKahnException, GeneralizedSchurDecomposition::GSDException)
  {
    // set extended Steady State

    for (size_t i = 0; i < zeta_back_mixed.size(); i++)
      llXsteadyState(i) = steadyState(zeta_back_mixed[i]);

    for (size_t i = 0; i < n_endo; i++)
      llXsteadyState(zeta_back_mixed.size() + i) = steadyState(i);

    for (size_t i = 0; i < zeta_fwrd_mixed.size(); i++)
      llXsteadyState(zeta_back_mixed.size() + n_endo + i) = steadyState(zeta_fwrd_mixed[i]);

    //get jacobian
    dynamicDLLp.eval(llXsteadyState, Mx, deepParams, steadyState, residual, &jacobian, NULL, NULL);

    //compute rules
    decisionRules.compute(jacobian, ghx, ghu);
  }
  template <class VEC>
  void ComputeSteadyState(VEC &steadyState, const VectorView &deepParams)
  {
    // does nothig for time being.
  }


};

#endif // !defined(5ADFF920_9C74_46f5_9FE9_88AD4D4BBF19__INCLUDED_)
