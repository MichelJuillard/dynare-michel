// $Id: kalman_filter.cpp 532 2005-11-30 13:51:33Z kamenik $

/*
* Copyright (C) 2008-2009 Dynare Team
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

/* derived from c++kalman_filter library by O. Kamenik */

// This provides an interface to KalmanTask::filter.

/****************************************************** 
% kalman_filter.cpp : Defines the entry point for 
% Computing the likelihood of a stationnary state space model.
% It is called from Dynare DsgeLikelihood.m, 
%
% function [loglik per d vll] = kalman_filter_dll(T,R,Q,H,Y,start,a, Z, P. [Pinf | u/U flag]
%
% INPUTS
%    T                      [double]    mm*mm transition matrix of the state equation.
%    R                      [double]    mm*rr matrix, mapping structural innovations to state variables.
%    Q                      [double]    rr*rr covariance matrix of the structural innovations.
%    H                      [double]    pp*pp (or 1*1 =0 if no measurement error) covariance matrix of the measurement errors. 
%    Y                      [double]    pp*smpl matrix of (detrended) data, where pp is the maximum number of observed variables.
%    start                  [integer]   scalar, likelihood evaluation starts at 'start'.
%    Z                      [double]    pp*mm matrix mapping state to pp observations 
%    a                      [vector]    mm vector of initial state, usually of 0s
%    P                      [double]    mm*mm variance-covariance matrix with stationary variables
%    Pinf   [optional]      [double]    mm*mm variance-covariance matrix with stationary variables
%    u/U    [optional]      [char]      u/U univariate
% OUTPUTS
%    loglik                 [double]    scalar, total likelihood
%    per                    [int]       number of succesfully filtered periods; if no error
%                           [int]       then per equals to the number of columns of Y
%    d                                  number of initial periods for which the state is
%                                       still diffuse (d is always 0 for non-diffuse case)
%    vll                    [double]    vector, density of observations in each period.
%
% REFERENCES
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003, in Journal of Time Series
%   Analysis, vol. 24(1), pp. 85-98).
%
% NOTES
%   The vector "vll" is used to evaluate the jacobian of the likelihood.
**********************************************************/

#include "kalman.h"
#include "ts_exception.h"

#include "GeneralMatrix.h"
#include "Vector.h"
#include "SylvException.h"

#include "mex.h"

extern "C" {
  void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[])
    {

    if (nrhs < 9 || nrhs > 11)
      mexErrMsgTxt("Must have  9, 10 or 11 input parameters.\n");
    if (nlhs < 1 || nlhs > 4)
      mexErrMsgTxt("Must have 1, 2, 3 or 4 output parameters.\n");
    //int start = 1; // default start of likelihood calculation
    // test for univariate case
    bool uni = false;
    const mxArray* const last = prhs[nrhs-1];
    if (mxIsChar(last)
      && ((*mxGetChars(last)) == 'u' || (*mxGetChars(last)) == 'U'))
      uni = true;
    
    // test for diffuse case
    bool diffuse = false;
    if ((mxIsChar(last) && nrhs == 11) ||
      (!mxIsChar(last) && nrhs == 10))
      diffuse = true;
    
    try {
      // make input matrices
      GeneralMatrix T(mxGetPr(prhs[0]), mxGetM(prhs[0]), mxGetN(prhs[0]));
      GeneralMatrix R(mxGetPr(prhs[1]), mxGetM(prhs[1]), mxGetN(prhs[1]));
      GeneralMatrix Q(mxGetPr(prhs[2]), mxGetM(prhs[2]), mxGetN(prhs[2]));
      GeneralMatrix H(mxGetPr(prhs[3]), mxGetM(prhs[3]), mxGetN(prhs[3]));
      GeneralMatrix Y(mxGetPr(prhs[4]), mxGetM(prhs[4]), mxGetN(prhs[4]));
      int start = (int)mxGetScalar(prhs[5]);
      GeneralMatrix Z(mxGetPr(prhs[6]), mxGetM(prhs[6]), mxGetN(prhs[6]));
      GeneralMatrix a(mxGetPr(prhs[7]), mxGetM(prhs[7]), mxGetN(prhs[7]));
      GeneralMatrix P(mxGetPr(prhs[8]), mxGetM(prhs[8]), mxGetN(prhs[8]));

      int nper=Y.numCols();
#ifdef DEBUG		
    mexPrintf("kalman_filter: periods=%d start=%d, a.length=%d, uni=%d diffuse=%d\n", nper, start,a.numRows(), uni, diffuse);
#endif		

      // make storage for output
      double loglik;
      int per;
      int d;
      // create state init
      StateInit* init = NULL;
      std::vector<double>* vll=new std::vector<double> (nper);
      if (diffuse||uni) 
        {
        if (diffuse) 
          {
          GeneralMatrix Pinf(mxGetPr(prhs[9]), mxGetM(prhs[9]), mxGetN(prhs[9]));
          init = new StateInit(P, Pinf, a.getData());
          } 
        else 
          {
          init = new StateInit(P, a.getData());
          }
        // fork, create objects and do filtering
        KalmanTask kt(Y, Z, H, T, R, Q, *init);
        if (uni) 
          {
          KalmanUniTask kut(kt);
          loglik = kut.filter(per, d, (start-1), vll);
          per = per / Y.numRows();
          d = d / Y.numRows();
          } 
        else 
          {
#ifdef TIMING_LOOP
  for (int tt=0;tt<1000;++tt)
    {
#endif
          loglik = kt.filter(per, d, (start-1), vll);
#ifdef TIMING_LOOP
    }
    mexPrintf("kalman_filter: finished 1000 loops");
#endif
          }
        }
      else // basic Kalman
        {
        init = new StateInit(P, a.getData());
        BasicKalmanTask bkt(Y, Z, H, T, R, Q, *init);
#ifdef TIMING_LOOP
  for (int tt=0;tt<1000;++tt)
    {
#endif
        loglik = bkt.filter( per, d, (start-1), vll);
#ifdef DEBUG		
    mexPrintf("Basickalman_filter: loglik=%d \n", loglik);
#endif		
#ifdef TIMING_LOOP
    }
    mexPrintf("Basickalman_filter: finished 1000 loops");
#endif

        }
        // destroy init
      delete init;
      
      // create output and upload output data
      if (nlhs >= 1)
        plhs[0] = mxCreateDoubleScalar(loglik);
      if (nlhs >= 2) {
        plhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
        (*((int*)mxGetData(plhs[1]))) = per;
        }
      if (nlhs >= 3) {
        plhs[2] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
        (*((int*)mxGetData(plhs[2]))) = d;
        }
      if (nlhs >= 4)
        {
        // output full log-likelihood array
        /* Set the output pointer to the  array of log likelihood. */
        plhs[3] = mxCreateDoubleMatrix(nper,1, mxREAL);
        double * mxll= mxGetPr(plhs[3]);
        // allocate likelihood array
        for (int j=0;j<nper;++j)
          mxll[j]=(*vll)[j];
        }
      delete vll;
      
      } 
    catch (const TSException& e) 
      {
      mexErrMsgTxt(e.getMessage());
      } 
    catch (SylvException& e) 
      {
      char mes[300];
      e.printMessage(mes, 299);
      mexErrMsgTxt(mes);
      }

    } // mexFunction
  }; // extern 'C'
