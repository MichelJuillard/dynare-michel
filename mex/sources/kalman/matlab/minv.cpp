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

// This provides a matrix inversion

/****************************************************** 
% function [Tinv] = minv(T)
%
% INPUTS
%    T                      [double]    mm*mm transition matrix of the state equation.

% OUTPUTS
%    Tinverse               [double]    mm*mm transition matrix of the state equation.
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

    if (nrhs < 1 )
      mexErrMsgTxt("Must have  min 1 input parameters.\n");
    if (nlhs < 1 )
      mexErrMsgTxt("Must have  min 1 output parameters.\n")
      
      ;
    //int start = 1; // default start of likelihood calculation
    // test for univariate case
    try 
      {
      // make input matrices
      ConstGeneralMatrix T(mxGetPr(prhs[0]), mxGetM(prhs[0]), mxGetN(prhs[0]));
      // create output and upload output data
      plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), mxGetM(prhs[0]), mxREAL);
//      double * mxinv= mxGetPr(plhs[0]);
      GeneralMatrix Tinv(mxGetPr(plhs[0]), mxGetM(prhs[0]), mxGetN(prhs[0]));
//    Tinv.unit();
//    Tinv.zeros();


      // make storage for output

#ifdef TIMING_LOOP
      int loops=1000;
      if (nrhs >1 )
         loops = (int)mxGetScalar(prhs[1]);
      for (int tt=0;tt<loops;++tt)
      {
#endif
        Tinv.unit();
        T.multInvLeft(Tinv);
        //Tinv.print();

#ifdef TIMING_LOOP
      }
      mexPrintf("minv: finished: %d loops\n",loops);
#endif
      // create output and upload output data
/*      if (nlhs >= 1)
        {
        plhs[0] = mxCreateNumericMatrix(mxGetM(prhs[0]), mxGetM(prhs[0]), mxINT32_CLASS, mxREAL);
        double * mxinv= mxGetPr(plhs[0]);
        // allocate likelihood array
        for (int j=0;j<nper;++j)
          mxinv[j]=(*vll)[j];
        }
*/      
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
