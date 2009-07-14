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

#include <iostream>
using namespace std;

#include "kalman.h"
#include "ts_exception.h"

#include "GeneralMatrix.h"
#include "Vector.h"
#include "SylvException.h"


#include "mex.h"

/*************************************
* This main() is for testing kalman DLL entry point by linking to 
* the kalman library statically and passing its hard-coded data: 
* parameters, covar, 
***************************************/


int main(int argc, char* argv[])
//main (int nrhs, mxArray* plhs[],
    {
    int nrhs=9; 
    int nlhs=4;
    if (nrhs < 9 || nrhs > 11)
      mexErrMsgTxt("Must have  9, 10 or 11 input parameters.\n");
    if (nlhs < 1 || nlhs > 4)
      mexErrMsgTxt("Must have 1, 2, 3 or 4 output parameters.\n");
    //int start = 1; // default start of likelihood calculation
    // test for univariate case
    bool uni = false;
//    const mxArray* const last = prhs[nrhs-1];
//    if (mxIsChar(last)
//      && ((*mxGetChars(last)) == 'u' || (*mxGetChars(last)) == 'U'))
//      uni = true;
    
    // test for diffuse case
    bool diffuse = false;
//    if ((mxIsChar(last) && nrhs == 11) ||
//      (!mxIsChar(last) && nrhs == 10))
//      diffuse = true;
    
    double Tmat[]={// need to pass transposed matrices!!??
      
         0,         0,         0,         0,         0,         0,         0,         0,
   -0.0013,    0.5000,    0.0000,   -0.0000,    0.0188,   -0.0013,    0.1182,   -0.0017,
    0.2158,    0.0000,    0.9502,   -0.0000,    0.0127,    0.2158,    0.0438,   -0.0088,
    0.0273,   -0.0000,   -0.0000,    0.8522,   -0.1260,   -0.8249,   -0.4720,    0.0356,
   -0.0716,   -0.0000,    0.0000,    0.0000,    0.5491,   -0.0716,   -0.9573,   -0.0935,
   -0.0000,   -0.0000,    0.0000,   -0.0000,   -0.0000,   -0.0000,    0.0000,   -0.0000,
         0,         0,         0,         0,         0,         0,         0,         0,
    0.6464,    0.0000,   -0.0000,   -0.0000,    0.0573,    0.6464,    0.2126,    0.8441
      
      
 //     0	,-0.001294119891461,	0.21578807493606	,0.027263201686985,	-0.071633450625617,	-0,	0,	0.646379181371765,
 //     0,	0.5,	0,	-0,	-0,	-0,	0,	0,
 //     0,	0,	0.9502,	-0,	0,	0,	0,	-0,
 //     0,	-0,	-0,	0.8522,	0,	-0,	0,	-0,
 //     0,	0.018758765757513,	0.012692095232426,	-0.126035674083997,	0.549074256326045,	-0,	0,	0.05730910985981,
 //     0,	-0.001294119891461,	0.21578807493606,	-0.824936798313015,	-0.071633450625617,	-0,	0,	0.646379181371766,
 //     0,	0.118192240459753,	0.04380066554165,	-0.471963836695487,	-0.957255289691476,	0,	0,	0.212592467520726,
 //     0,	-0.00168993250228,	-0.008835241183444,	0.035601779209991,	-0.093542875943306,	-0,	0,	0.844077271823789
      };

      double  Rmat[]={// need to pass transposed matrices!!??
      0.2271,         0,    1.0000,         0,    0.0134,    0.2271,    0.0461,   -0.0093,
      0.0320,         0,         0,    1.0000,   -0.1479,   -0.9680,   -0.5538,    0.0418,
     -0.0026,    1.0000,         0,         0,    0.0375,   -0.0026,    0.2364,   -0.0034,
     -0.0895,         0,         0,         0,    0.6863,   -0.0895,   -1.1966,   -0.1169

//            0.2271,    0.0320,   -0.0026,   -0.0895,
//                 0,         0,    1.0000,         0,
//            1.0000,         0,         0,         0,
//                 0,    1.0000,         0,         0,
//            0.0134,   -0.1479,    0.0375,    0.6863,
//            0.2271,   -0.9680,   -0.0026,   -0.0895,
//            0.0461,   -0.5538,    0.2364,   -1.1966,
//           -0.0093,    0.0418,   -0.0034,   -0.1169
        };
      double  Qmat[]={
         0.0931,    0,         0,         0,
         0,    0.1849,         0,         0,
         0,         0,    0.0931,         0,
         0,         0,         0,    0.0100
        };

      double  Zmat[]={ // need to pass transposed matrices!!??
          0,     0,     1,     0,
         0,     0,     0,     0,
         0,     0,     0,     0,
         0,     0,     0,     0,
         0,     0,     0 ,    1,
         0,     0,     0,     0,
         1,     0,     0,     0,
         0,     1,     0,     0
//         0,     0,     0,     0,     0,     0,     1,     0,
//         0,     0,     0,     0,     0,     0,     0,     1,
//         1,     0,     0,     0,     0,     0,     0,     0,
//         0,     0,     0,     0,     1,     0,     0,     0
        };
      double Ymat[]={
       -0.4073,	0.2674,	0.2896,	0.0669,	0.1166,	-0.1699,	-0.2518,	-0.0562,	-0.3269,-0.0703,-0.1046,	-0.4888	,-0.3524,	-0.2485	,-0.587,	-0.4546,	-0.397,	-0.2353,	-0.0352	-0.2171,	-0.3754,	-0.4322,	-0.4572,	-0.4903,	-0.4518,	-0.6435,	-0.6304	,-0.4148,	-0.2892,	-0.4318,	-0.601,	-0.4148,	-0.4315,	-0.3531,	-0.8053,	-0.468,	-0.4263,
        3.1739,	3.738 ,	3.8285,	3.3342,	3.7447,	3.783,	3.1039,	2.8413,	3.0338,	0.3669,	0.0847	,0.0104,	0.2115,	-0.6649,	-0.9625,	-0.733,	-0.8664,	-1.4441,	-1.0179,	-1.2729	,-1.9539,	-1.4427,	-2.0371,	-1.9764,	-2.5654,	-2.857,	-2.5842,	-3.0427,	-2.8312,	-2.332,	-2.2768,	-2.1816,	-2.1043,	-1.8969,	-2.2388,	-2.1679,	-2.1172,
        3.2174,	3.1903,	3.3396,	3.1358,	2.8625,	3.3546,	2.4609,	1.9534,	0.9962,	-0.7904,-1.1672,	-1.2586,	-1.3593,	-1.3443	,-0.9413,	-0.6023,	-0.4516,	-0.5129,	-0.8741,	-1.0784,	-1.4091,	-1.3627,	-1.5731,	-1.6037	-1.8814,	-2.1482	,-1.3597,	-1.1855,	-1.1122,	-0.8424,	-0.9747,	-1.1385,	-1.4548,	-1.4284,	-1.4633,	-1.0621,	-0.7871,
        0.8635,	0.9058,	0.7656,	0.7936,	0.8631,	0.9074,	0.9547,	1.2045,	1.085,	0.9178,	0.5242,	0.3178	,0.1472,	0.0227,	-0.0799,	-0.0611,	-0.014,	0.1132,	0.1774,	0.0782,	0.0436,	-0.1596,	-0.2691,	-0.2895,	-0.3791,	-0.402,	-0.4166	,-0.4037,	-0.3636,	-0.4075,	-0.4311,	-0.447,	-0.5111,	-0.6274,	-0.7261,	-0.6974,	-0.5012

        };

    try {
      // make input matrices
      GeneralMatrix T(Tmat, 8, 8);
      GeneralMatrix R(Rmat, 8, 4);
      GeneralMatrix Q(Qmat, 4, 4);
      GeneralMatrix H(4, 4);
      H.zeros();
/*********use simlated data for time being *********/
      GeneralMatrix Y( 4, 109);
      Y.zeros();
      for (int i=0;i<4;++i)
        {
        for (int j=0;j<109;++j)
          {
          Y.get(i,j)=  ((double) ( rand() % 10 -5.0))/2.0;
#ifdef DEBUG		
    mexPrintf("Y [%d %d] =%f, \n", i, j,Y.get(i,j));
#endif
   
          }
        }
/***********
      GeneralMatrix Y(Ymat, 4, 109);
      for (int i=0;i<4;++i)
        {
        for (int j=0;j<109;++j)
          {
#ifdef DEBUG		
    mexPrintf("Y [%d %d] =%f, \n", i, j,Y.get(i,j));
#endif
          }
        }
***********/   
      double riccatiTol=0.000001;
      int start = 1;
      GeneralMatrix Z(Zmat, 4, 8);
      GeneralMatrix a(8, 1);
      a.zeros();
      GeneralMatrix P( 8, 8);
      P.zeros();
      for (int i=0;i<8;++i)
        P.get(i,i)=10.0;

      int nper=Y.numCols();
#ifdef DEBUG		
    mexPrintf("kalman_filter: periods=%d start=%d, a.length=%d, uni=%d diffuse=%d\n", nper, start,a.numRows(), uni, diffuse);
#endif		

      // make storage for output
      double loglik=-1.1111;
      int per;
      int d;
      // create state init
      StateInit* init = NULL;
      std::vector<double>* vll=new std::vector<double> (nper);

      if (diffuse||uni) 
        {
        if (diffuse) 
          {
          GeneralMatrix Pinf(P.numRows(),P.numCols());
          Pinf.zeros();
          init = new StateInit(P, Pinf, a.getData());
          } 
        else 
          {
          init = new StateInit(P, a.getData());
          }
        // fork, create objects and do filtering
  #ifdef TIMING_LOOP
    for (int tt=0;tt<10000;++tt)
      {
  #endif
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
          loglik = kt.filter(per, d, (start-1), vll);
          }
  #ifdef TIMING_LOOP
    //    mexPrintf("kalman_filter: finished %d loops", tt);
      }
      mexPrintf("kalman_filter: finished 10,000 loops");
  #endif

        }
      else // basic Kalman
        {
        init = new StateInit(P, a.getData());
        BasicKalmanTask bkt(Y, Z, H, T, R, Q, *init, riccatiTol);
#ifdef TIMING_LOOP
  for (int tt=0;tt<10000;++tt)
    {
#endif 
        loglik = bkt.filter( per, d, (start-1), vll);
#ifdef DEBUG		
//    mexPrintf("Basickalman_filter: loglik=%f \n", loglik);
//    cout << "loglik " << loglik << "\n";
#endif		
#ifdef TIMING_LOOP
    }
    mexPrintf("Basickalman_filter: finished 10,000 loops");
#endif

        }


        // destroy init
      delete init;
      mexPrintf("logLik = %f \n", loglik);
      delete vll;
      // create output and upload output data
/************
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
        // Set the output pointer to the  array of log likelihood. 
        plhs[3] = mxCreateDoubleMatrix(nper,1, mxREAL);
        double * mxll= mxGetPr(plhs[3]);
        // allocate likelihood array
        for (int j=0;j<nper;++j)
          mxll[j]=(*vll)[j];
        }
******************************/      
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

  //  } // mexFunction
  }; // main  extern 'C'
