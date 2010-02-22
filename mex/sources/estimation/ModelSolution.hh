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
  class ModelSolutionException: public TSException
  {
  public:
    const lapack_int info, n;
    ModelSolutionException(lapack_int info_arg, lapack_int n_arg, std::string& mes) 
      : TSException(__FILE__, __LINE__, mes), info(info_arg), n(n_arg) {};
  };

  virtual ~ModelSolution();
  ModelSolution(const std::string& modName,  size_t n_endo, size_t n_exo, const std::vector<size_t>& zeta_fwrd_arg, 
    const std::vector<size_t>& zeta_back_arg, const std::vector<size_t>& zeta_mixed_arg, 
    const std::vector<size_t>& zeta_static_arg, const Matrix& llincidence, double qz_criterium);
  void compute( Vector& steadyState, const Vector& deepParams, 	Matrix& ghx, Matrix& ghu );

private:
  const int n_endo;
  const int n_exo;
  const int n_jcols; // Num of Jacobian columns
  const Matrix ll_incidence; // leads and lags indices
  Matrix jacobian;
  Vector residual;
  Matrix Mx;
  DecisionRules decisionRules;
  DynamicModelDLL*  dynamicDLLp;
  //Matrix jacobian;
  void ComputeModelSolution( Vector& steadyState, const Vector& deepParams, 	Matrix& ghx, Matrix& ghu );
  void ComputeSteadyState( Vector& steadyState, const Vector& deepParams);

};




#endif // !defined(5ADFF920_9C74_46f5_9FE9_88AD4D4BBF19__INCLUDED_)
