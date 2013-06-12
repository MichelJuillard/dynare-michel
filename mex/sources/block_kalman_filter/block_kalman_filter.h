/*
 * Copyright (C) 2007-2012 Dynare Team
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

#ifndef BLOCK_KALMAN_FILTER
#define BLOCK_KALMAN_FILTER

#ifndef DEBUG_EX
# include <dynmex.h>
#else
# include "mex_interface.hh"
#endif

#include <dynblas.h>
#include <dynlapack.h>
using namespace std;

class BlockKalmanFilter
{
  public:
  mxArray *pT , *pR, *pQ, *pH, *pP, *pY, *pQQ, *pv, *pF, *piF, *p_P_t_t1, *pK, *p_K_P;
  double *T , *R, *Q , *H, *Y, *mfd, *QQ, *v, *F, *iF;
  int start, pure_obs, smpl, n, n_state, n_shocks, H_size;
  double kalman_tol, riccati_tol, dF, LIK, Inf, pi;
  lapack_int pp, lw, info;

  double* nz_state_var;
  int *i_nz_state_var, *mf;
  int n_diag, t;
  mxArray *M_;
  mxArray* pa, *p_tmp, *p_tmp1, *plik;
  double *tmp_a, *tmp1, *lik, *v_n, *w, *oldK;
  bool notsteady, F_singular, missing_observations;
  lapack_int *iw, *ipiv;
  double anorm, rcond;
  lapack_int size_d_index;
  int no_more_missing_observations, number_of_observations;
  const mxArray* pdata_index;
  vector<int> d_index;
  const mxArray* pd_index;
  double* dd_index;

  public:
  BlockKalmanFilter(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], double *P_mf[], double *v_pp[], double *K[], double *v_n[], double *a[], double *K_P[], double *P_t_t1[], double *tmp[], double *P[]);
  bool block_kalman_filter(int nlhs, mxArray *plhs[], double *P_mf, double *v_pp, double *K, double *v_n, double *a, double *K_P, double *P_t_t1, double *tmp, double *P);
  void block_kalman_filter_ss(double *P_mf, double *v_pp, double *K, double *v_n, double *a, double *K_P, double *P_t_t1, double *tmp, double *P);
  void return_results_and_clean(int nlhs, mxArray *plhs[], double *P_mf[], double *v_pp[], double *K[], double *v_n[], double *a[], double *K_P[], double *P_t_t1[], double *tmp[], double *P[]);
};
#endif
