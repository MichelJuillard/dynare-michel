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
class ModelSolution{

public:
  ModelSolution(const std::string& dynamicDllFile,  size_t n_endo, size_t n_exo, const std::vector<size_t>& zeta_fwrd_arg, 
    const std::vector<size_t>& zeta_back_arg, const std::vector<size_t>& zeta_mixed_arg, 
    const std::vector<size_t>& zeta_static_arg, double qz_criterium);
  virtual ~ModelSolution(){};
  void compute( VectorView& steadyState, const Vector& deepParams, 	Matrix& ghx, Matrix& ghu ) throw (DecisionRules::BlanchardKahnException, GeneralizedSchurDecomposition::GSDException);

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
  void ComputeModelSolution( VectorView& steadyState, const Vector& deepParams, 	Matrix& ghx, Matrix& ghu ) throw (DecisionRules::BlanchardKahnException, GeneralizedSchurDecomposition::GSDException);
  void ComputeSteadyState( VectorView& steadyState, const Vector& deepParams);

};




#endif // !defined(5ADFF920_9C74_46f5_9FE9_88AD4D4BBF19__INCLUDED_)
