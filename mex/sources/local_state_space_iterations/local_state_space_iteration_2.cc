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

/*
 * This mex file computes particles at time t+1 given particles and innovations at time t,
 * using a second order approximation of the nonlinear state space model.
 */

#include <iostream>
#include <cstring>
#include <vector>
#include <dynmex.h>
#include <dynblas.h>

#ifdef USE_OMP
#include <omp.h>
#endif

using namespace std;

#define FIRST_ORDER_LOOP 1// Comment out this line to use mkl-blas instead of loops when computing ghx*yhat and ghu*epsilon

void set_vector_of_indices(const int n, const int r, vector<int> &v1, vector<int> &v2, vector<int> &v3)
{
  const int m = n*(n+1)/2;
  v1.resize(m,0);
  v2.resize(m,0);
  v3.resize(m,0);
  for(int i=0, index=0, jndex=0;i<n; i++)
    {
      jndex+=i;
      for(int j=i; j<n; j++, index++, jndex++)
        {
          v1[index] = i;
          v2[index] = j;
          v3[index] = jndex*r;
        }
    }
}

void ss2Iteration_pruning(double* y2, double* y1, const double* yhat2, const double* yhat1, const double *epsilon,
                  double* ghx, double* ghu,
                  const double* constant, const double* ghxx, const double* ghuu, const double* ghxu, const double* ss,
			  const blas_int m, const blas_int n, const blas_int q, const blas_int s, const int number_of_threads)
{
  #ifndef FIRST_ORDER_LOOP
      const char transpose[2] = "N";
      const double one = 1.0;
      const blas_int ONE = 1;
  #endif
  vector<int> ii1, ii2, ii3;// vector indices for ghxx
  vector<int> jj1, jj2, jj3;// vector indices for ghuu
  set_vector_of_indices(n, m, ii1, ii2, ii3);
  set_vector_of_indices(q, m, jj1, jj2, jj3);
  #ifdef USE_OMP
  #pragma omp parallel for num_threads(number_of_threads)
  #endif
  for (int particle = 0; particle<s; particle++)
  {
    int particle_ = particle*m;
    int particle__ = particle*n;
    int particle___ = particle*q;
    memcpy(&y2[particle_],&constant[0],m*sizeof(double));
    memcpy(&y1[particle_],&ss[0],m*sizeof(double));
    #ifndef FIRST_ORDER_LOOP
        dgemv(transpose,&m,&n,&one,&ghx[0],&m,&yhat2[particle__],&ONE,&one,&y2[particle_],&ONE);
	dgemv(transpose,&m,&q,&one,&ghu[0],&m,&epsilon[particle___],&ONE,&one,&y2[particle_],&ONE);
    #endif
    for (int variable = 0; variable<m; variable++)
      {
        int variable_ = variable + particle_;
        // +ghx*yhat2+ghu*u
        #ifdef FIRST_ORDER_LOOP
            for (int column = 0, column_=0; column<q; column++, column_ += m)
	      {
		int i1 = variable+column_;
		int i2 = column+particle__;
		int i3 = column+particle___;
		y2[variable_] += ghx[i1]*yhat2[i2];
		y2[variable_] += ghu[i1]*epsilon[i3];
	      }
	    for (int column = q, column_=q*m; column<n; column++, column_ += m)
	      {
		y2[variable_] += ghx[variable+column_]*yhat2[column+particle__];
	      }
        #endif
        // +ghxx*kron(yhat1,yhat1)
        for(int i=0; i<n*(n+1)/2; i++)
          {
            int i1 = particle__+ii1[i];
            int i2 = particle__+ii2[i];
            if(i1==i2)
              {
                y2[variable_] += .5*ghxx[variable+ii3[i]]*yhat1[i1]*yhat1[i1];
              }
            else
              {
                y2[variable_] += ghxx[variable+ii3[i]]*yhat1[i1]*yhat1[i2];
              }
          }
        // +ghuu*kron(u,u)
        for(int j=0; j<q*(q+1)/2; j++)
          {
            int j1 = particle___+jj1[j];
            int j2 = particle___+jj2[j];
            if(j1==j2)
              {
                y2[variable_] += .5*ghuu[variable+jj3[j]]*epsilon[j1]*epsilon[j1];
              }
            else
              {
                y2[variable_] += ghuu[variable+jj3[j]]*epsilon[j1]*epsilon[j2];
              }
          }
        // +ghxu*kron(yhat1,u)
        for (int v = particle__, i = 0; v<particle__+n; v++)
          for (int s = particle___; s<particle___+q; s++, i += m)
            y2[variable_] += ghxu[variable+i]*epsilon[s]*yhat2[v];
            #ifdef FIRST_ORDER_LOOP
                for (int column = 0, column_=0; column<q; column++, column_ += m)
		  {
		    int i1 = variable+column_;
		    int i2 = column+particle__;
		    int i3 = column+particle___;
		    y1[variable_] += ghx[i1]*yhat1[i2];
		    y1[variable_] += ghu[i1]*epsilon[i3];
		  }
		for (int column = q, column_=q*m; column<n; column++, column_ += m)
		  {
		    y1[variable_] += ghx[variable+column_]*yhat1[column+particle__];
		  }
            #endif
      }
      #ifndef FIRST_ORDER_LOOP
          dgemv(transpose,&m,&n,&one,&ghx[0],&m,&yhat1[particle__],&ONE,&one,&y1[particle_],&ONE);
	  dgemv(transpose,&m,&q,&one,&ghu[0],&m,&epsilon[particle___],&ONE,&one,&y1[particle_],&ONE);
      #endif
  }
}

void ss2Iteration(double* y, const double* yhat, const double *epsilon,
                  double* ghx, double* ghu,
                  const double* constant, const double* ghxx, const double* ghuu, const double* ghxu,
                  const blas_int m, const blas_int n, const blas_int q, const blas_int s, const int number_of_threads)
{
  #ifndef FIRST_ORDER_LOOP
      const char transpose[2] = "N";
      const double one = 1.0;
      const blas_int ONE = 1;
  #endif
  vector<int> ii1, ii2, ii3;// vector indices for ghxx
  vector<int> jj1, jj2, jj3;// vector indices for ghuu
  set_vector_of_indices(n, m, ii1, ii2, ii3);
  set_vector_of_indices(q, m, jj1, jj2, jj3);
  #ifdef USE_OMP
  #pragma omp parallel for num_threads(number_of_threads)
  #endif
  for (int particle = 0; particle<s; particle++)
  {
    int particle_ = particle*m;
    int particle__ = particle*n;
    int particle___ = particle*q;
    memcpy(&y[particle_],&constant[0],m*sizeof(double));
    #ifndef FIRST_ORDER_LOOP
        dgemv(transpose,&m,&n,&one,&ghx[0],&m,&yhat[particle__],&ONE,&one,&y[particle_],&ONE);
        dgemv(transpose,&m,&q,&one,&ghu[0],&m,&epsilon[particle___],&ONE,&one,&y[particle_],&ONE);
    #endif
    for (int variable = 0; variable<m; variable++)
      {
        int variable_ = variable + particle_;
        // +ghx*yhat+ghu*u
        #ifdef FIRST_ORDER_LOOP
            for (int column = 0, column_=0; column<q; column++, column_ += m)
	      {
		int i1 = variable+column_;
		int i2 = column+particle__;
		int i3 = column+particle___;
		y[variable_] += ghx[i1]*yhat[i2];
		y[variable_] += ghu[i1]*epsilon[i3];
	      }
	    for (int column = q, column_=q*m; column<n; column++, column_ += m)
	      {
		y[variable_] += ghx[variable+column_]*yhat[column+particle__];
	      }
        #endif
        // +ghxx*kron(yhat,yhat)
	for(int i=0; i<n*(n+1)/2; i++)
          {
            int i1 = particle__+ii1[i];
            int i2 = particle__+ii2[i];
            if(i1==i2)
              {
                y[variable_] += .5*ghxx[variable+ii3[i]]*yhat[i1]*yhat[i1];
              }
            else
              {
                y[variable_] += ghxx[variable+ii3[i]]*yhat[i1]*yhat[i2];
              }
          }
        // +ghuu*kron(u,u)
        for(int j=0; j<q*(q+1)/2; j++)
          {
            int j1 = particle___+jj1[j];
            int j2 = particle___+jj2[j];
            if(j1==j2)
              {
                y[variable_] += .5*ghuu[variable+jj3[j]]*epsilon[j1]*epsilon[j1];
              }
            else
              {
                y[variable_] += ghuu[variable+jj3[j]]*epsilon[j1]*epsilon[j2];
              }
          }
        // +ghxu*kron(yhat,u)
        for (int v = particle__, i = 0; v<particle__+n; v++)
          for (int s = particle___; s<particle___+q; s++, i += m)
            y[variable_] += ghxu[variable+i]*epsilon[s]*yhat[v];
      }
  }
}




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*
  ** prhs[0] yhat          [double]  n*s array, time t particles.
  ** prhs[1] epsilon       [double]  q*s array, time t innovations.
  ** prhs[2] ghx           [double]  m*n array, first order reduced form.
  ** prhs[3] ghu           [double]  m*q array, first order reduced form.
  ** prhs[4] constant      [double]  m*1 array, deterministic steady state + second order correction for the union of the states and observed variables.
  ** prhs[5] ghxx          [double]  m*n^2 array, second order reduced form.
  ** prhs[6] ghuu          [double]  m*q^2 array, second order reduced form.
  ** prhs[7] ghxu          [double]  m*nq array, second order reduced form.
  ** prhs[8] yhat_         [double]  [OPTIONAL] n*s array, time t particles (pruning additional latent variables).
  ** prhs[9] ss            [double]  [OPTIONAL] m*1 array, steady state for the union of the states and the observed variables (needed for the pruning mode).
  **
  ** plhs[0] y             [double]  n*s array, time t+1 particles.
  ** plhs[1] y_            [double]  n*s array, time t+1 particles for the pruning latent variables.
  **
  */

  // Check the number of input and output.
  if ((nrhs != 9) && (nrhs != 11))
    {
      mexErrMsgTxt("Eight or ten input arguments are required.");
    }
  if (nlhs > 2)
    {
      mexErrMsgTxt("Too many output arguments.");
    }
  // Get dimensions.
  size_t n = mxGetM(prhs[0]);// Number of states.
  size_t s = mxGetN(prhs[0]);// Number of particles.
  size_t q = mxGetM(prhs[1]);// Number of innovations.
  size_t m = mxGetM(prhs[2]);// Number of elements in the union of states and observed variables.
  //mexPrintf("\n s (the number of column of yhat) is equal to %d.", s);
  //mexPrintf("\n The number of column of epsilon is %d.", mxGetN(prhs[1]));
  // Check the dimensions.
  if (
      (s != mxGetN(prhs[1]))   || // Number of columns for epsilon
      (n != mxGetN(prhs[2]))   || // Number of columns for ghx
      (m != mxGetM(prhs[3]))   || // Number of rows for ghu
      (q != mxGetN(prhs[3]))   || // Number of columns for ghu
      (m != mxGetM(prhs[4]))   || // Number of rows for 2nd order constant correction + deterministic steady state
      (m != mxGetM(prhs[5]))   || // Number of rows for ghxx
      (n*n != mxGetN(prhs[5])) || // Number of columns for ghxx
      (m != mxGetM(prhs[6]))   || // Number of rows for ghuu
      (q*q != mxGetN(prhs[6])) || // Number of columns for ghuu
      (m != mxGetM(prhs[7]))   || // Number of rows for ghxu
      (n*q != mxGetN(prhs[7]))    // Number of rows for ghxu
      )
    {
      mexErrMsgTxt("Input dimension mismatch!.");
    }
  if (nrhs>9)
    {
      if (
          (n != mxGetM(prhs[8]))   || // Number of rows for yhat_
          (s != mxGetN(prhs[8]))   || // Number of columns for yhat_
          (m != mxGetM(prhs[9]))      // Number of rows for ss
          )
        {
          mexErrMsgTxt("Input dimension mismatch!.");
        }
    }
  // Get Input arrays.
  double *yhat = mxGetPr(prhs[0]);
  double *epsilon = mxGetPr(prhs[1]);
  double *ghx = mxGetPr(prhs[2]);
  double *ghu = mxGetPr(prhs[3]);
  double *constant = mxGetPr(prhs[4]);
  double *ghxx = mxGetPr(prhs[5]);
  double *ghuu = mxGetPr(prhs[6]);
  double *ghxu = mxGetPr(prhs[7]);
  double *yhat_ = NULL, *ss = NULL;
  if (nrhs>9)
    {
      yhat_ = mxGetPr(prhs[8]);
      ss = mxGetPr(prhs[9]);
    }
  if (nrhs==9)
    {
      int numthreads = (int) mxGetScalar(prhs[8]);
      double *y;
      plhs[0] = mxCreateDoubleMatrix(m, s, mxREAL);
      y = mxGetPr(plhs[0]);
      ss2Iteration(y, yhat, epsilon, ghx, ghu, constant, ghxx, ghuu, ghxu, (int) m, (int) n, (int) q, (int) s, numthreads);
    }
  else
    {
      int numthreads = (int) mxGetScalar(prhs[10]);
      double *y, *y_;
      plhs[0] = mxCreateDoubleMatrix(m, s, mxREAL);
      plhs[1] = mxCreateDoubleMatrix(m, s, mxREAL);
      y = mxGetPr(plhs[0]);
      y_ = mxGetPr(plhs[1]);
      ss2Iteration_pruning(y, y_, yhat, yhat_, epsilon, ghx, ghu, constant, ghxx, ghuu, ghxu, ss, (int) m, (int) n, (int) q, (int) s, numthreads);
    }
}
