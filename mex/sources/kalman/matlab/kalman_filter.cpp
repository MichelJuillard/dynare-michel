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

#include "kalman.h"
#include "ts_exception.h"

#include "GeneralMatrix.h"
#include "Vector.h"
#include "SylvException.h"

#include "mex.h"

extern "C" {
  void mexFunction(int nhls, mxArray* plhs[],
    int nhrs, const mxArray* prhs[])
    {
    if (nhrs < 8 || nhrs > 10)
      mexErrMsgTxt("Must have 8, 9, or 10 input parameters.\n");
    if (nhls < 1 || nhls > 3)
      mexErrMsgTxt("Must have 1, 2, or 3 output parameters.\n");
    
    // test for univariate case
    bool uni = false;
    const mxArray* const last = prhs[nhrs-1];
    if (mxIsChar(last)
      && ((*mxGetChars(last)) == 'u' || (*mxGetChars(last)) == 'U'))
      uni = true;
    
    // test for diffuse case
    bool diffuse = false;
    if ((mxIsChar(last) && nhrs == 10) ||
      (!mxIsChar(last) && nhrs == 9))
      diffuse = true;
    
    try {
      // make input matrices
      GeneralMatrix Z(mxGetPr(prhs[0]), mxGetM(prhs[0]), mxGetN(prhs[0]));
      GeneralMatrix H(mxGetPr(prhs[1]), mxGetM(prhs[1]), mxGetN(prhs[1]));
      GeneralMatrix T(mxGetPr(prhs[2]), mxGetM(prhs[2]), mxGetN(prhs[2]));
      GeneralMatrix R(mxGetPr(prhs[3]), mxGetM(prhs[3]), mxGetN(prhs[3]));
      GeneralMatrix Q(mxGetPr(prhs[4]), mxGetM(prhs[4]), mxGetN(prhs[4]));
      GeneralMatrix Y(mxGetPr(prhs[5]), mxGetM(prhs[5]), mxGetN(prhs[5]));
      GeneralMatrix a(mxGetPr(prhs[6]), mxGetM(prhs[6]), mxGetN(prhs[6]));
      GeneralMatrix P(mxGetPr(prhs[7]), mxGetM(prhs[7]), mxGetN(prhs[7]));
      
      // make storage for output
      double loglik;
      int per;
      int d;
      // create state init
      StateInit* init = NULL;
      if (diffuse) {
        GeneralMatrix Pinf(mxGetPr(prhs[8]), mxGetM(prhs[8]), mxGetN(prhs[8]));
        init = new StateInit(P, Pinf, a.getData());
        } else {
        init = new StateInit(P, a.getData());
          }
        // fork, create objects and do filtering
        KalmanTask kt(Y, Z, H, T, R, Q, *init);
        if (uni) {
          KalmanUniTask kut(kt);
          loglik = kut.filter(per, d);
          per = per / Y.numRows();
          d = d / Y.numRows();
          } else {
          loglik = kt.filter(per, d);
            }
          // destroy init
          delete init;
          
          // create output and upload output data
          if (nhls >= 1)
            plhs[0] = mxCreateDoubleScalar(loglik);
          if (nhls >= 2) {
            plhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
            (*((int*)mxGetData(plhs[1]))) = per;
            }
          if (nhls == 3) {
            plhs[2] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
            (*((int*)mxGetData(plhs[2]))) = d;
            }
      } catch (const TSException& e) {
      mexErrMsgTxt(e.getMessage());
        } catch (SylvException& e) {
        char mes[300];
        e.printMessage(mes, 299);
        mexErrMsgTxt(mes);
          }
    }
  };
