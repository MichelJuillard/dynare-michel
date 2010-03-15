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

// Test  for ModelSolution
// Uses fs2000k2.mod and its ..._dynamic.mexw32 

#include "ModelSolution.hh"

int
main (int argc, char** argv)
{
  if (argc < 2)
    {
      std::cerr << argv[0] << ": please provide as argument the name of the dynamic DLL generated from fs2000k2.mod (typically fs2000k2_dynamic.mex*)" << std::endl;
      exit(EXIT_FAILURE);
    }

  std::string modName = argv[1];
  const int npar = 7;
  const size_t n_endo=15, n_exo=2; 
  std::vector<size_t> zeta_fwrd_arg; 
  std::vector<size_t> zeta_back_arg; 
  std::vector<size_t> zeta_mixed_arg;
  std::vector<size_t> zeta_static_arg;
  double qz_criterium=1.0+1.0e-9;  
  Vector steadyState(n_endo), deepParams(npar);

  double dYSparams [] = {
    1.0110,  2.2582,  0.4477,  1.0000,
    4.5959,  1.0212,  5.8012,  0.8494,
    0.1872,  0.8604,  1.0030,  1.0080,
    0.5808,  1.0030,  2.2093
  };
  double vcov[] = {
    0.1960e-3, 0.0,
    0.0, 0.0250e-3
  };
  int nVCVpar = 2;
  MatrixView vCovVW(vcov,nVCVpar,nVCVpar,nVCVpar);
  Matrix vCov (nVCVpar, nVCVpar);
  vCov = vCovVW;

  Matrix ll_incidence(3,n_endo); // leads and lags indices
  double inllincidence[]={
    1,   5,  0,
    2,   6,  20,
    0,   7,  21,
    0,   8,   0,
    0,   9,   0,
    0,  10,   0,
    3,  11,   0,
    0,  12,   0,
    0,  13,   0,
    0,  14,   0,
    0,  15,   0,
    0,  16,   0,
    4,  17,   0,
    0,  18,   0,
    0,  19,  22,
  };
  MatrixView llincidence(inllincidence, 3, n_endo,3); // leads and lags indices
  ll_incidence= llincidence;

  double dparams[] = { 0.3300,
    0.9900,
    0.0030,
    1.0110,
    0.7000,
    0.7870,
    0.0200};

  VectorView modParamsVW (dparams, npar,1);
  deepParams=modParamsVW;
  VectorView steadyStateVW(dYSparams,n_endo,1);
  steadyState=steadyStateVW;
  std::cout << "Vector deepParams: " << std::endl << deepParams << std::endl;
  std::cout << "Matrix vCov: " << std::endl << vCov << std::endl;
  std::cout << "MatrixVw llincidence: " << std::endl << llincidence << std::endl;
  std::cout << "Matrix ll_incidence: " << std::endl << ll_incidence << std::endl;
  std::cout << "Vector steadyState: " << std::endl << steadyState << std::endl;

  // Set zeta vectors [0:(n-1)] from Matlab indices [1:n]
  //order_var = [ stat_var(:); pred_var(:); both_var(:); fwrd_var(:)];
  size_t statc[]={ 4, 5, 6, 8, 9, 10, 11, 12, 14};
  size_t back[]={1,7,13};
  size_t both[]={2};
  size_t fwd[]={ 3,15};
  for (int i=0;i<9;++i)
    zeta_static_arg.push_back(statc[i]-1);
  for (int i=0;i<3;++i)
    zeta_back_arg.push_back(back[i]-1);
  for (int i=0;i<1;++i)
    zeta_mixed_arg.push_back(both[i]-1);
  for (int i=0;i<2;++i)
    zeta_fwrd_arg.push_back(fwd[i]-1);

  Matrix ghx(n_endo, zeta_back_arg.size() + zeta_mixed_arg.size());
  Matrix ghu(n_endo,n_exo);

  ModelSolution modelSolution( modName, n_endo, n_exo
    , zeta_fwrd_arg, zeta_back_arg, zeta_mixed_arg, zeta_static_arg, ll_incidence, qz_criterium);

  modelSolution.compute(steadyState, deepParams, ghx,  ghu);

  std::cout << "Matrix ghx: " << std::endl << ghx << std::endl;
  std::cout << "Matrix ghu: " << std::endl << ghu << std::endl;
}
