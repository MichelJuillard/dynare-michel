/*
 * Copyright (C) 2010-2013 Dynare Team
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

#include "InitializeKalmanFilter.hh"

int
main(int argc, char **argv)
{
  if (argc < 2)
    {
      std::cerr << argv[0] << ": please provide as argument the name of the dynamic DLL generated from fs2000k2e.mod (typically fs2000k2e_dynamic.mex*)" << std::endl;
      exit(EXIT_FAILURE);
    }

  std::string modName = argv[1];
  const int npar = 7; //(int)mxGetM(mxFldp);
  const size_t n_endo = 15, n_exo = 2;
  std::vector<size_t> zeta_fwrd_arg;
  std::vector<size_t> zeta_back_arg;
  std::vector<size_t> zeta_mixed_arg;
  std::vector<size_t> zeta_static_arg;
  //std::vector<size_t>
  double qz_criterium = 1.000001; //1.0+1.0e-9;
  Vector steadyState(n_endo), deepParams(npar);

  double dYSparams [] = {
    1.000199998312523,
    0.993250551764778,
    1.006996670195112,
    1,
    2.718562165733039,
    1.007250753636589,
    18.982191739915155,
    0.860847884886309,
    0.316729149714572,
    0.861047883198832,
    1.00853622757204,
    0.991734328394345,
    1.355876776121869,
    1.00853622757204,
    0.992853374047708
  };

  double vcov[] = {
    0.001256631601,     0.0,
    0.0,        0.000078535044
  };

  double dparams[] = {
    0.3560,
    0.9930,
    0.0085,
    1.0002,
    0.1290,
    0.6500,
    0.0100
  };

  VectorView modParamsVW(dparams, npar, 1);
  deepParams = modParamsVW;
  VectorView steadyStateVW(dYSparams, n_endo, 1);
  steadyState = steadyStateVW;
  std::cout << "Vector deepParams: " << std::endl << deepParams << std::endl;
  std::cout << "Vector steadyState: " << std::endl << steadyState << std::endl;

  // Set zeta vectors [0:(n-1)] from Matlab indices [1:n] so that:
  // order_var = [ stat_var(:); pred_var(:); both_var(:); fwrd_var(:)];
  size_t statc[] = { 4, 5, 6, 8, 9, 10, 11, 12, 14};
  size_t back[] = {1, 7, 13};
  size_t both[] = {2};
  size_t fwd[] = { 3, 15};
  for (int i = 0; i < 9; ++i)
    zeta_static_arg.push_back(statc[i]-1);
  for (int i = 0; i < 3; ++i)
    zeta_back_arg.push_back(back[i]-1);
  for (int i = 0; i < 1; ++i)
    zeta_mixed_arg.push_back(both[i]-1);
  for (int i = 0; i < 2; ++i)
    zeta_fwrd_arg.push_back(fwd[i]-1);

  size_t varobs[] = {12, 11};
  std::vector<size_t> varobs_arg;
  for (size_t i = 0; i < 2; ++i)
    varobs_arg.push_back(varobs[i]-1);

  // Compute zeta_varobs_back_mixed
  sort(varobs_arg.begin(), varobs_arg.end());
  std::vector<size_t> zeta_back_mixed, zeta_varobs_back_mixed;
  set_union(zeta_back_arg.begin(), zeta_back_arg.end(),
            zeta_mixed_arg.begin(), zeta_mixed_arg.end(),
            back_inserter(zeta_back_mixed));
  set_union(zeta_back_mixed.begin(), zeta_back_mixed.end(),
            varobs_arg.begin(), varobs_arg.end(),
            back_inserter(zeta_varobs_back_mixed));

  size_t n_vbm = zeta_varobs_back_mixed.size();

  Matrix T(n_vbm, n_vbm), R(n_vbm, n_exo),
    RQRt(n_vbm, n_vbm), Pstar(n_vbm, n_vbm),
    Pinf(n_vbm, n_vbm), Q(n_exo);

  MatrixView vCovVW(vcov, n_exo, n_exo, n_exo);
  Q = vCovVW;

  double lyapunov_tol = 1e-16;
  int info = 0;
  size_t nobs = 2;
  Matrix yView(nobs, 192); // dummy
  yView.setAll(0.2);
  const MatrixConstView dataView(yView, 0,  0, nobs, yView.getCols()); // dummy
  Matrix yDetrendView(nobs, yView.getCols()); // dummy
  MatrixView dataDetrendView(yDetrendView, 0,  0, nobs, yDetrendView.getCols()); // dummy

  const Vector xparams1(0); // dummy
  double penalty = 1e8;

  InitializeKalmanFilter initializeKalmanFilter(modName, n_endo, n_exo,
                                                zeta_fwrd_arg, zeta_back_arg, zeta_mixed_arg, zeta_static_arg,
                                                zeta_varobs_back_mixed, qz_criterium,
                                                lyapunov_tol, info);

  std::cout << "Initialize KF with Q: " << std::endl << Q << std::endl;

  initializeKalmanFilter.initialize(steadyStateVW, deepParams, R, Q, RQRt, T, Pstar, Pinf,
                                    penalty, dataView, dataDetrendView, info);

  std::cout << "Matrix T: " << std::endl << T << std::endl;
  std::cout << "Matrix R: " << std::endl << R << std::endl;
  std::cout << "Matrix RQRt: " << std::endl << RQRt << std::endl;
  std::cout << "Matrix Pstar: " << std::endl << Pstar << std::endl;
  std::cout << "Matrix Pinf: " << std::endl << Pinf << std::endl;

}

