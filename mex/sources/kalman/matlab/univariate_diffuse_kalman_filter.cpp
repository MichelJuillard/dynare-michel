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
% function [LIK, lik] = kalman_filter_dll(T,R,Q,H,P,Y,start,Z/mf[, kalman_tol,riccati_tol])
%
% INPUTS
%    T                      [double]    mm*mm transition matrix of the state equation.
%    R                      [double]    mm*rr matrix, mapping structural innovations to state variables.
%    Q                      [double]    rr*rr covariance matrix of the structural innovations.
%    H                      [double]    pp*pp (or 1*1 =0 if no measurement error) covariance matrix of the measurement errors. 
%    P                      [double]    mm*mm variance-covariance matrix with stationary variables
%    Y                      [double]    pp*smpl matrix of (detrended) data, where pp is the maximum number of observed variables.
%    start                  [integer]   scalar, likelihood evaluation starts at 'start'.
%    Z or mf:               [double]   Z: pp*mm matrix mapping state to pp observations 
%  Alternative parameters
%    mf                     [integer]   pp*1 vector of indices - alternative to Z matrix.
%  Additional optional parameters
%    kalman_tol             [double]    scalar, tolerance parameter (rcond).
%    riccati_tol            [double]    scalar, tolerance parameter (riccati iteration).
%
% OUTPUTS
%    LIK        [double]    scalar, likelihood
%    lik        [double]    vector, density of observations in each period.
%
% REFERENCES
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003, in Journal of Time Series
%   Analysis, vol. 24(1), pp. 85-98).
%
% NOTES
%   The vector "lik" is used to evaluate the jacobian of the likelihood.
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
    if (nrhs < 8 || nrhs > 10)
      mexErrMsgTxt("Must have 8, 9, or 10 input parameters.\n");
    if (nlhs < 1 || nlhs > 3)
      mexErrMsgTxt("Must have 1, 2, or 3 output parameters.\n");
    
    int start = 1; // default start of likelihood calculation
    try 
      {
      // make input matrices
      GeneralMatrix T(mxGetPr(prhs[0]), mxGetM(prhs[0]), mxGetN(prhs[0]));
      GeneralMatrix R(mxGetPr(prhs[1]), mxGetM(prhs[1]), mxGetN(prhs[1]));
      GeneralMatrix Q(mxGetPr(prhs[2]), mxGetM(prhs[2]), mxGetN(prhs[2]));
      GeneralMatrix H(mxGetPr(prhs[3]), mxGetM(prhs[3]), mxGetN(prhs[3]));
      GeneralMatrix Pinf(mxGetPr(prhs[4]), mxGetM(prhs[4]), mxGetN(prhs[4]));
      GeneralMatrix P(mxGetPr(prhs[5]), mxGetM(prhs[5]), mxGetN(prhs[5]));
      GeneralMatrix Y(mxGetPr(prhs[6]), mxGetM(prhs[6]), mxGetN(prhs[6]));
      if (nrhs>6) start = (int)mxGetScalar(prhs[7]);
      GeneralMatrix Z(mxGetPr(prhs[8]), mxGetM(prhs[8]), mxGetN(prhs[8]));
      int nper = mxGetN(prhs[5]); // no of periods
      GeneralMatrix a( mxGetN(prhs[0]), 1);// initiate inital state to 0s
      a.zeros();
#ifdef DEBUG		
      mexPrintf("kalman_filter: periods = %d ", nper);
#endif		
     
      // make storage for output
      double loglik;
      int per;
      int d;
      // output for full log-likelihood array
      std::vector<double>* vll=new std::vector<double> (nper);
      // create state init
      StateInit* init = NULL;
      init = new StateInit(P, Pinf, a.getData());
      // fork, create objects and do filtering
      KalmanTask kt(Y, Z, H, T, R, Q, *init);
      // developement of the output.
#ifdef DEBUG		
      mexPrintf("kalman_filter: running and filling outputs.\n");
#endif			
      KalmanUniTask kut(kt);
      loglik = kut.filter(per, d, (start-1), vll);
      per = per / Y.numRows();
      d = d / Y.numRows();
      // destroy init
      delete init;
      
      // create output and upload output data
      if (nlhs >= 1)
        plhs[0] = mxCreateDoubleScalar(loglik);
      if (nlhs >= 2)
        {
        // output full log-likelihood array
        /* Set the output pointer to the  array of log likelihood. */
        plhs[1] = mxCreateDoubleMatrix(nper,1, mxREAL);
        double * mxll= mxGetPr(plhs[1]);
        // allocate likelihood array
        for (int j=0;j<nper;++j)
          mxll[j]=(*vll)[j];
        }
      if (nlhs >= 3) 
        {
        plhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
        (*((int*)mxGetData(plhs[2]))) = per;
        }
      if (nlhs == 4) 
        {
        plhs[2] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
        (*((int*)mxGetData(plhs[3]))) = d;
        }
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
    }
  };
