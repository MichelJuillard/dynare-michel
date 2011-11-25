/*
 * Copyright (C) 2007-2011 Dynare Team
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
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
//#define BLAS

#define DIRECT

#ifndef DEBUG_EX
# include <dynmex.h>
#else
# include "mex_interface.hh"
#endif

#include <dynblas.h>
#include <dynlapack.h>

#define BLOCK

void
mexDisp(mxArray* P)
{
  unsigned int n = mxGetN(P);
  unsigned int m = mxGetM(P);
  double *M = mxGetPr(P);
  mexPrintf("%d x %d\n", m, n);
  mexEvalString("drawnow;");
  for (unsigned int i = 0; i < m; i++)
    {
      for (unsigned int j = 0; j < n; j++)
        mexPrintf(" %9.4f",M[i+ j * m]);
      mexPrintf("\n");
    }
}

/*
(T,R,Q,H,P,Y,start,mf,kalman_tol,riccati_tol, block)
% Computes the likelihood of a stationnary state space model.
%
% INPUTS
%    T                      [double]    n*n transition matrix of the state equation.
%    R                      [double]    n*rr matrix, mapping structural innovations to state variables.
%    Q                      [double]    rr*rr covariance matrix of the structural innovations.
%    H                      [double]    pp*pp (or 1*1 =0 if no measurement error) covariance matrix of the measurement errors. 
%    P                      [double]    n*n variance-covariance matrix with stationary variables
%    Y                      [double]    pp*smpl matrix of (detrended) data, where pp is the maximum number of observed variables.
%    start                  [integer]   scalar, likelihood evaluation starts at 'start'.
%    mf                     [integer]   pp*1 vector of indices.
%    kalman_tol             [double]    scalar, tolerance parameter (rcond).
%    riccati_tol            [double]    scalar, tolerance parameter (riccati iteration).
%
% OUTPUTS
%    LIK        [double]    scalar, MINUS loglikelihood
%    lik        [double]    vector, density of observations in each period.
%
% REFERENCES
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003, in Journal of Time Series
%   Analysis, vol. 24(1), pp. 85-98).
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.

% Copyright (C) 2004-2011 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
%global options_
smpl = size(Y,2);                               % Sample size.
n   = size(T,2);                               % Number of state variables.
pp   = size(Y,1);                               % Maximum number of observed variables.
a    = zeros(n,1);                             % State vector.
dF   = 1;                                       % det(F).
QQ   = R*Q*transpose(R);                        % Variance of R times the vector of structural innovations.
t    = 0;                                       % Initialization of the time index.
lik  = zeros(smpl,1);                           % Initialization of the vector gathering the densities.
LIK  = Inf;                                     % Default value of the log likelihood.
oldK = Inf;
notsteady   = 1;                                % Steady state flag.
F_singular  = 1;
if block
    %nz_state_var = M_.nz_state_var;
    while notsteady && t<smpl
        t  = t+1;
        v  = Y(:,t)-a(mf);
        F  = P(mf,mf) + H;
        if rcond(F) < kalman_tol
            if ~all(abs(F(:))<kalman_tol)
                return
            else
                a = T*a;
                P = T*P*transpose(T)+QQ;
            end
        else
            F_singular = 0;
            dF     = det(F);
            iF     = inv(F);
            lik(t) = log(dF)+transpose(v)*iF*v;
            K      = P(:,mf)*iF;
            a      = T*(a+K*v);
            P = block_pred_Vcov_KF(mf, P, K, T, QQ);
            %P      = T*(P-K*P(mf,:))*transpose(T)+QQ;
            notsteady = max(abs(K(:)-oldK)) > riccati_tol;
            oldK = K(:);
        end
    end;
else
    while notsteady && t<smpl
        t  = t+1;
        v  = Y(:,t)-a(mf);
        F  = P(mf,mf) + H;
        if rcond(F) < kalman_tol
            if ~all(abs(F(:))<kalman_tol)
                return
            else
                a = T*a;
                P = T*P*transpose(T)+QQ;
            end
        else
            F_singular = 0;
            dF     = det(F);
            iF     = inv(F);
            lik(t) = log(dF)+transpose(v)*iF*v;
            K      = P(:,mf)*iF;
            a      = T*(a+K*v);
            P      = T*(P-K*P(mf,:))*transpose(T)+QQ;
            notsteady = max(abs(K(:)-oldK)) > riccati_tol;
            oldK = K(:);
        
        end
    end
end

if F_singular
    error('The variance of the forecast error remains singular until the end of the sample')
end

if t < smpl
    t0 = t+1;
    while t < smpl
        t = t+1;
        v = Y(:,t)-a(mf);
        a = T*(a+K*v);
        lik(t) = transpose(v)*iF*v;
    end
    lik(t0:smpl) = lik(t0:smpl) + log(dF);
end    

% adding log-likelihhod constants
lik = (lik + pp*log(2*pi))/2;

LIK = sum(lik(start:end)); % Minus the log-likelihood.*/



bool
not_all_abs_F_bellow_crit(double* F, int size, double crit)
{
   int i = 0;
   while (i < size && abs(F[i])<crit)
     {
       i++;
     }
   if (i < size)
     return false;
   else
     return true;
}


double
det(double* F, int dim)
{
  double det = 1.0;
  for (int i = 0; i < dim; i++)
    det *= F[i * (1 + dim)];
  return det;
}

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nlhs > 3)
    DYN_MEX_FUNC_ERR_MSG_TXT("kalman_filter provides at most 3 output argument.");
  if (nrhs != 13)
    DYN_MEX_FUNC_ERR_MSG_TXT("kalman_filter requires exactly 13 input arguments.");
//(T,R,Q,H,P,Y,start,mf,kalman_tol,riccati_tol, block)
  mxArray *pT = mxDuplicateArray(prhs[0]);
  mxArray *pR = mxDuplicateArray(prhs[1]);
  mxArray *pQ = mxDuplicateArray(prhs[2]);
  mxArray *pH = mxDuplicateArray(prhs[3]);
  mxArray *pP = mxDuplicateArray(prhs[4]);
  mxArray *pY = mxDuplicateArray(prhs[5]);
  
  double *T = mxGetPr(pT);
  double *R = mxGetPr(pR);
  double *Q = mxGetPr(pQ);
  double *H = mxGetPr(pH);
  double *P = mxGetPr(pP);
  double *Y = mxGetPr(pY);
  int start = mxGetScalar(prhs[6]);
  double* mfd = (double*)mxGetData(prhs[7]);
  
    
  double kalman_tol = mxGetScalar(prhs[8]);
  double riccati_tol = mxGetScalar(prhs[9]);
  int pure_obs = mxGetScalar(prhs[12]);
  
  /* Reading the sparse structure of matrix in (§A) */
  typedef struct 
  {
      int indx_1;
      int indx_2;
      int indx_3;
  } t_Var;

  /*Defining the initials values*/
  int smpl = mxGetN(pY);                                 // Sample size.          ;
  int n   = mxGetN(pT);                                  // Number of state variables.
  lapack_int pp   = mxGetM(pY);                          // Maximum number of observed variables.
  int n_state = n - pure_obs;
  int n_shocks = mxGetM(pQ);


#ifdef DIRECT  
  double* nz_state_var = (double*)mxGetData(prhs[10]);
  int *i_nz_state_var = (int*)mxMalloc(n*sizeof(int));
  for (int i = 0; i < n; i++)
    i_nz_state_var[i] = nz_state_var[i];

  int n_diag = mxGetScalar(prhs[11]);
#else
  mxArray *M_;
  M_ = mexGetVariable("global", "M_");
  if (M_ == NULL)
    mexErrMsgTxt(" in main, global variable not found: M_\n");
  mxArray *mxa = mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "fname"));
  int buflen = mxGetM(mxa) * mxGetN(mxa) + 1;
  char *fname;
  fname = (char *) mxCalloc(buflen+1, sizeof(char));
  int status = mxGetString(mxa, fname, buflen);
  fname[buflen] = ' ';
  if (status != 0)
    mexWarnMsgTxt("Not enough space. Filename is truncated.");
  string file_name = fname;
  file_name.append(".kfi");
  
  ifstream KF_index_file;
  KF_index_file.open(file_name.c_str(), ios::in | ios::binary);
  int n_diag;
  KF_index_file.read(reinterpret_cast<char *>(&n_diag), sizeof(n_diag));
  /*mexPrintf("n_diag=%d\n",n_diag);*/
  
  int size_index_KF;
  KF_index_file.read(reinterpret_cast<char *>(&size_index_KF), sizeof(size_index_KF));
  //mexPrintf("size_index_KF=%d\n", size_index_KF);
  t_Var *Var = (t_Var*)mxMalloc(size_index_KF * sizeof(t_Var));
  KF_index_file.read(reinterpret_cast<char *>(Var), size_index_KF * sizeof(t_Var));
  
  int size_index_KF_2;
  KF_index_file.read(reinterpret_cast<char *>(&size_index_KF_2), sizeof(size_index_KF_2));
  //mexPrintf("size_index_KF_2=%d\n", size_index_KF_2);
  t_Var *Var_2 = (t_Var*)mxMalloc(size_index_KF_2 * sizeof(t_Var));
  KF_index_file.read(reinterpret_cast<char *>(Var_2), size_index_KF_2 * sizeof(t_Var));
  KF_index_file.close();
#endif

  mxArray* pa = mxCreateDoubleMatrix(n, 1, mxREAL);         // State vector.
  double* a = mxGetPr(pa);
  double* tmp_a = (double*)mxMalloc(n * sizeof(double));
  double dF = 0.0;                                            // det(F).
  mxArray* p_tmp = mxCreateDoubleMatrix(n, n_state, mxREAL);
  double *tmp = mxGetPr(p_tmp);
  mxArray* p_tmp1 = mxCreateDoubleMatrix(n, n_shocks, mxREAL);
  double *tmp1 = mxGetPr(p_tmp1);
  int t = 0;                                               // Initialization of the time index.
  mxArray* plik  = mxCreateDoubleMatrix(smpl, 1, mxREAL);
  double* lik = mxGetPr(plik);
  double Inf =  mxGetInf();                                   
  double LIK  = 0.0;                                          // Default value of the log likelihood.

  bool notsteady   = true;                                 // Steady state flag.
  bool F_singular  = true;
  double* v_pp = (double*)mxMalloc(pp * sizeof(double));
  double* v_n = (double*)mxMalloc(n * sizeof(double));
  int* mf = (int*)mxMalloc(pp * sizeof(int));
  for (int i = 0; i < pp; i++)
    mf[i] = mfd[i] - 1;
  double pi = atan2((double)0.0,(double)-1.0);
  
  /*compute QQ = R*Q*transpose(R)*/                        // Variance of R times the vector of structural innovations.;
  // tmp = R * Q;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n_shocks; j++)
      {
        double res = 0.0;
        for (int k = 0; k < n_shocks; k++)
          res += R[i + k * n] * Q[j * n_shocks + k];
        tmp1[i + j * n] = res;
      }

  // QQ = tmp * transpose(R)
  mxArray* pQQ = mxCreateDoubleMatrix(n, n, mxREAL);
  double* QQ = mxGetPr(pQQ);
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      {
        double res = 0.0;
        for (int k = 0; k < n_shocks; k++)
          res += tmp1[i + k * n] * R[k * n + j];
        QQ[i + j * n] = QQ[j + i * n] = res;
      }
  mxDestroyArray(p_tmp1);

  mxArray* pv = mxCreateDoubleMatrix(pp, 1, mxREAL);
  double* v = mxGetPr(pv);
  mxArray* pF =  mxCreateDoubleMatrix(pp, pp, mxREAL);
  double* F = mxGetPr(pF);
  mxArray* piF =  mxCreateDoubleMatrix(pp, pp, mxREAL);
  double* iF = mxGetPr(piF);
  lapack_int lw = pp * 4;
  double* w = (double*)mxMalloc(lw * sizeof(double));
  lapack_int* iw = (lapack_int*)mxMalloc(pp * sizeof(lapack_int));
  lapack_int* ipiv = (lapack_int*)mxMalloc(pp * sizeof(lapack_int));
  lapack_int info = 0;
  double anorm, rcond;
  #ifdef BLAS
  mxArray* p_P_t_t1 = mxCreateDoubleMatrix(n, n, mxREAL);
  #else
  mxArray* p_P_t_t1 = mxCreateDoubleMatrix(n_state, n_state, mxREAL);
  #endif
  double* P_t_t1 = mxGetPr(p_P_t_t1);
  mxArray* pK = mxCreateDoubleMatrix(n, pp, mxREAL);
  double* K = mxGetPr(pK);
  #ifdef BLAS
  mxArray* p_K_P = mxCreateDoubleMatrix(n, n, mxREAL);
  #else
  mxArray* p_K_P = mxCreateDoubleMatrix(n_state, n_state, mxREAL);
  #endif
  double* K_P = mxGetPr(p_K_P);
  double* oldK  = (double*)mxMalloc(n * pp * sizeof(double));
  double* P_mf = (double*)mxMalloc(n * pp * sizeof(double));
  for (int i = 0; i < pp; i++)
    oldK[i] = Inf;
  
  while (notsteady && t < smpl)
    {
      //v = Y(:,t) - a(mf)
      for (int i = 0; i < pp; i++)
        v[i]  = Y[i + t * pp] - a[mf[i]];
      
      //F  = P(mf,mf) + H;
      for (int i = 0; i < pp; i++)
        for (int j = 0; j < pp; j++)
          iF[i + j * pp] = F[i + j * pp] = P[mf[i] + mf[j] * n] + H[i + j * pp];
      
      /* Computes the norm of x */ 
      double anorm = dlange("1", &pp, &pp, iF, &pp, w); 

      /* Modifies F in place with a LU decomposition */ 
      dgetrf(&pp, &pp, iF, &pp, ipiv, &info); 
      if (info != 0) fprintf(stderr, "dgetrf failure with error %d\n", (int) info); 
 
      /* Computes the reciprocal norm */ 
      dgecon("1", &pp, iF, &pp, &anorm, &rcond, w, iw, &info); 
      if (info != 0) fprintf(stderr, "dgecon failure with error %d\n", (int) info); 
      
      if (rcond < kalman_tol)
        if (not_all_abs_F_bellow_crit(F, pp * pp, kalman_tol))   //~all(abs(F(:))<kalman_tol)
          {
            mexPrintf("error: F singular\n");
            double LIK  = Inf;
            if (nlhs >= 1)
              {
                plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
                double* pind = mxGetPr(plhs[0]);
                pind[0] = LIK;
              }
            if (nlhs == 2)
              {
                for (int i = t; i < smpl; i++)
                  lik[i] = Inf;
                plhs[1] = plik;
              }
            return;
          }
        else
          {
            mexPrintf("F singular\n");
            //a = T*a;
            for (int i = 0; i < n; i++)
              {
                double res = 0.0;
                for (int j = pure_obs; j < n; j++)
                  res += T[i + j *n] * a[j];
                tmp_a[i] = res;
              }
           memcpy(a, tmp_a, n * sizeof(double));

           //P = T*P*transpose(T)+QQ;
           memset(tmp, 0, n * n_state * sizeof(double));
#ifdef DIRECT
           for (int i = 0; i < n; i++)
            for (int j = pure_obs; j < n; j++)
              {
                int j1 = j - pure_obs;
                int j1_n_state = j1 * n_state - pure_obs;
                //if ((i < pp) || (i >= n_diag + pp) || (j1 >= n_diag))
                  for (int k = pure_obs; k < i_nz_state_var[i]; k++)
                    {
                      tmp[i + j1 * n ] += T[i + k * n] * P[k + j1_n_state];
                    }
              }
          /*for (int i = pp; i < pp + n_diag; i++)
             for (int j = pp; j < pp + n_diag; j++)
               tmp[i + (j - pp) * n] = T[i + i * n] * P_t_t1[j - pp + (i - pp) * n_state];*/
          memset(P, 0, n * n * sizeof(double));
          int n_n_obs = n * pure_obs;
          for (int i = 0; i < n; i++)
            for (int j = i; j < n; j++)
              {
                //if ((i < pp) || (i >= n_diag + pp) || (j < pp) || (j >= n_diag + pp))
                  for (int k = pure_obs; k < i_nz_state_var[j]; k++)
                    {
                      int k_n = k * n;
                      P[i * n + j] += tmp[i + k_n - n_n_obs] * T[j + k_n];
                    }
              }
#else
           #pragma omp parallel for shared(T, P_t_t1, tmp, Var) num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
           //for (int i = 0; i < n; i++)
             for (int j = 0; j < size_index_KF; j++)
               {
                 t_Var *s_Var = &Var[j];
                 tmp[s_Var->indx_1 ] += T[s_Var->indx_2] * P[s_Var->indx_3];
               }
           for (int i = pp; i < pp + n_diag; i++)
             for (int j = pp; j < pp + n_diag; j++)
               tmp[i + (j - pp) * n] = T[i + i * n] * P[j + i * n];
           memset(P, 0, n * n * sizeof(double));
           #pragma omp parallel for shared(T, tmp, P, Var_2) num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
           //for (int i = 0; i < n; i++)
             for (int j = 0; j < size_index_KF_2; j++)
               {
                 t_Var *s_Var = &Var_2[j];
                 P[s_Var->indx_1 /*+ i*/] += tmp[s_Var->indx_2 /*+ i*/] * T[s_Var->indx_3];
               }
#endif
          /*for (int i = pp; i < pp + n_diag; i++)
            for (int j = i; j < pp + n_diag; j++)
              P[j + i * n] = T[i + i * n] * T[j + j * n] * P[j + i * n];*/
              //P[j + i * n] = tmp[i +]
          
          for ( int i = 0; i < n; i++)
            {
              for ( int j = i ; j < n; j++)
                P[j + i * n] += QQ[j + i * n];
              for ( int j = i + 1; j < n; j++)
                P[i + j * n] = P[j + i * n];
            }
         }
      else
        {
          F_singular = false;
          
          //dF     = det(F);
          dF     = abs(det(iF, pp));

          //iF     = inv(F);
          //int lwork = 4/*2*/* pp;
          dgetri(&pp, iF, &pp, ipiv, w, &lw, &info);
          if (info != 0) fprintf(stderr, "dgetri failure with error %d\n", (int) info); 

          //lik(t) = log(dF)+transpose(v)*iF*v;
          for (int i = 0; i < pp; i++)
            {
              double res = 0.0;
              for (int j = 0; j < pp; j++)
                res += v[j] * iF[j  * pp + i];
              v_pp[i] = res;
            }
          double res = 0.0;
          for (int i = 0; i < pp; i++)
            res += v_pp[i] * v[i];

          lik[t] = (log(dF) + res + pp * log(2.0*pi))/2;
          if (t + 1 >= start)
            LIK += lik[t];

          //K      = P(:,mf)*iF;
          for (int i = 0; i < n; i++)
            for (int j = 0; j < pp; j++)
              P_mf[i + j * n] = P[i + mf[j] * n];
          for (int i = 0; i < n; i++)
            for (int j = 0; j < pp; j++)
              {
                double res = 0.0;
                int j_pp = j * pp;
                for (int k = 0; k < pp; k++)
                  res += P_mf[i + k * n] * iF[j_pp + k];
                K[i + j * n] = res;
              }

          //a      = T*(a+K*v);
          for (int i = pure_obs; i < n; i++)
            {
              double res = 0.0;
              for (int j = 0; j < pp; j++)
                res += K[j  * n + i] * v[j];
              v_n[i] = res + a[i];
            }
          
          for (int i = 0; i < n; i++)
            {
              double res = 0.0;
              for (int j = pure_obs; j < n; j++)
                res += T[j  * n + i] * v_n[j];
              a[i] = res;
            }

          //P      = T*(P-K*P(mf,:))*transpose(T)+QQ;
          for (int i = 0; i < pp; i++)
            for (int j = pure_obs; j < n; j++)
              P_mf[i + j * pp] = P[mf[i] + j * n];
#ifdef BLAS
          #pragma omp parallel for shared(K, P, K_P) num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
          for (int i = 0; i < n; i++)
            for (int j = i; j < n; j++)
              {
                double res = 0.0;
                //int j_pp = j * pp;
                for (int k = 0; k < pp; k++)
                  res += K[i + k * n] * P_mf[k + j * pp];
                K_P[i * n + j] = K_P[j * n + i] = res;
              }
          //#pragma omp parallel for shared(P, K_P, P_t_t1)
          for (int i = pp; i < n; i++)
            for (int j = i; j < n; j++)
              {
                unsigned int k = i * n + j;
                P_t_t1[j * n + i] = P_t_t1[k] = P[k] - K_P[k];
              }
          double one  = 1.0;
          double zero = 0.0;
          memcpy(P, QQ, n * n *sizeof(double));
          dsymm("R", "U", &n, &n,
                &one, P_t_t1, &n,
                T, &n, &zero,
                tmp, &n);
          dgemm("N", "T", &n, &n,
                &n, &one, tmp, &n,
                T, &n, &one,
                P, &n);
          mexPrintf("P\n");
          mexDisp(pP);
#else
          //#pragma omp parallel for shared(K, P, K_P) num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
          
          for (int i = pure_obs; i < n; i++)
            {
              unsigned int i1 = i - pure_obs;
              for (int j = i ; j < n; j++)
                {
                  unsigned int j1 = j - pure_obs;
                  double res = 0.0;
                  int j_pp = j * pp;
                  for (int k = 0; k < pp; k++)
                    res += K[i + k * n] * P_mf[k + j_pp];
                  K_P[i1 * n_state + j1] = K_P[j1 * n_state + i1] = res;
                }
            }
          /*for (int j = pp; j < n; j++)
            {
              unsigned int j1 = j - pp;
              double res = P[mf[j] + j * n]
              for (int i = 0; i < n; i++)
                {
                  unsigned int i1 = i - pp;
                  K_P[i1 * n_state + j1] = res * K[i + j * n];
                }
              for (int i = pp; i < n; i++)
                {
                  unsigned int i1 = i - pp;
                  double res = 0.0;
                  for (int k = 0; k < pp; k++)
                    res += K[i + k * n] * P[mf[k] + j * n];
                  K_P[i1 * n_state + j1] = K_P[j1 * n_state + i1] = res;
                }
            }*/
            
          //#pragma omp parallel for shared(P, K_P, P_t_t1)
          for (int i = pure_obs; i < n; i++)
            {
              unsigned int i1 = i - pure_obs;
              for (int j = i; j < n; j++)
                {
                  unsigned int j1 = j - pure_obs;
                  unsigned int k1 = i1 * n_state + j1;
                  P_t_t1[j1 * n_state + i1] = P_t_t1[k1] = P[i * n + j] - K_P[k1];
                }
            }
          
          memset(tmp, 0, n * n_state * sizeof(double));
#ifdef DIRECT
          for (int i = 0; i < n; i++)
            {
              int max_k = i_nz_state_var[i];
              for (int j = pure_obs; j < n; j++)
                {
                  int j1 = j - pure_obs;
                  int j1_n_state = j1 * n_state - pure_obs;
                  int indx_tmp = i + j1 * n ;
                  //if ((i < pp) || (i >= n_diag + pp) || (j1 >= n_diag))
                  for (int k = pp; k < max_k; k++)
                    tmp[indx_tmp] += T[i + k * n] * P_t_t1[k + j1_n_state];
                }
            }
          memset(P, 0, n * n * sizeof(double));
          
          for (int i = 0; i < n; i++)
            {
              int n_n_obs = - n * pure_obs;
              for (int j = i; j < n; j++)
                {
                  int max_k = i_nz_state_var[j];
                  int P_indx = i * n + j;
                  //if ((i < pp) || (i >= n_diag + pp) || (j < pp) || (j >= n_diag + pp))
                  for (int k = pure_obs; k < max_k; k++)
                    {
                      int k_n = k * n;
                      P[P_indx] += tmp[i + k_n + n_n_obs] * T[j + k_n];
                    }
                }
            }
#else
          //#pragma omp parallel for shared(T, P_t_t1, tmp, Var) num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
          for (int j = 0; j < size_index_KF; j++)
            {
              t_Var *s_Var = &Var[j];
              tmp[s_Var->indx_1 ] += T[s_Var->indx_2 ] * P_t_t1[s_Var->indx_3];
            }
         for (int i = pp; i < pp + n_diag; i++)
           for (int j = pp; j < pp + n_diag; j++)
             tmp[i + (j - pp) * n] = T[i + i * n] * P_t_t1[j - pp + (i - pp) * n_state];
         memset(P, 0, n * n * sizeof(double));
         //#pragma omp parallel for shared(T, tmp, P, Var_2) num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
         for (int j = 0; j < size_index_KF_2; j++)
           {
             t_Var *s_Var = &Var_2[j];
             P[s_Var->indx_1 /*+ i*/] += tmp[s_Var->indx_2 /*+ i*/] * T[s_Var->indx_3];
           }
#endif
          for ( int i = 0; i < n; i++)
            {
              for ( int j = i ; j < n; j++)
                P[j + i * n] += QQ[j + i * n];
              for ( int j = i + 1; j < n; j++)
                P[i + j * n] = P[j + i * n];
            }
#endif
          //notsteady = max(abs(K(:)-oldK)) > riccati_tol;
          double max_abs = 0.0;
          for (int i = 0; i < n * pp; i++)
            {
              double res = abs(K[i] - oldK[i]);
              if (res > max_abs)
                max_abs = res;
            }
          notsteady = max_abs > riccati_tol;

          //oldK = K(:);
          memcpy(oldK, K, n * pp * sizeof(double));
        }
      t++;
    }

  if (F_singular)
    mexErrMsgTxt("The variance of the forecast error remains singular until the end of the sample\n");

  if (t+1 < smpl)
    {
      while (t < smpl)
        {
          //v = Y(:,t)-a(mf);
          for (int i = 0; i < pp; i++)
            v[i]  = Y[i + t * pp] - a[mf[i]];
          
          //a = T*(a+K*v);
          for (int i = pure_obs; i < n; i++)
            {
              double res = 0.0;
              for (int j = 0; j < pp; j++)
                res += K[j  * n + i] * v[j];
              v_n[i] = res + a[i];
            }
          for (int i = 0; i < n; i++)
            {
              double res = 0.0;
              for (int j = pure_obs; j < n; j++)
                res += T[j  * n + i] * v_n[j];
              a[i] = res;
            }
        
          //lik(t) = transpose(v)*iF*v;
          for (int i = 0; i < pp; i++)
            {
              double res = 0.0;
              for (int j = 0; j < pp; j++)
                res += v[j] * iF[j * pp + i];
              v_pp[i] = res;
            }
          double res = 0.0;
          for (int i = 0; i < pp; i++)
            res += v_pp[i] * v[i];

          lik[t] = (log(dF) + res + pp * log(2.0*pi))/2;
          if (t + 1 > start)
            LIK += lik[t];

          t++;
        }
    }

  // info = 0
  plhs[0] = mxCreateDoubleScalar(0);
  
  if (nlhs >= 2)
    {
      plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
      double* pind = mxGetPr(plhs[1]);
      pind[0] = LIK;
    }

  if (nlhs == 3)
    plhs[2] = plik;
  else
    mxDestroyArray(plik);
  
  mxFree(w);
#ifdef DIRECT
  mxFree(i_nz_state_var);
#else
  mxFree(Var);
  mxFree(Var_2);
#endif
  mxFree(tmp_a);
  mxFree(v_pp);
  mxFree(v_n);
  mxFree(mf);
  mxFree(w);
  mxFree(iw);
  mxFree(ipiv);
  mxFree(oldK);
  mxFree(P_mf);
  mxDestroyArray(pa);
  mxDestroyArray(p_tmp);
  mxDestroyArray(pQQ);
  mxDestroyArray(pv);
  mxDestroyArray(pF);
  mxDestroyArray(piF);
  mxDestroyArray(p_P_t_t1);
  mxDestroyArray(pK);
  mxDestroyArray(p_K_P);
}

