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
/****************************************************************
% entry Matlab DLL for function X=disclyap_fast(G,V,ch)
% 
% which solve the discrete Lyapunov Equation 
% X=G*X*G'+V 
% Using the Doubling Algorithm 
%
% If ch is defined then the code will check if the resulting X 
% is positive definite and generate an error message if it is not 
% 
****************************************************************/
#include "ts_exception.h"

#include "GeneralMatrix.h"
#include "SylvException.h"
#include "mex.h"
#include "disclyap_fast.h"

//void disclyap_fast(GeneralMatrix &G,GeneralMatrix & V, double tol= 1e-16, int ch=0);


extern "C" {
  void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[])
    {

    if (nrhs < 2 || nrhs > 4)
      mexErrMsgTxt("Must have  2, 3 or 4 input parameters.\n");
    if (nlhs != 1 )
      mexErrMsgTxt("Must have 1 output parameters.\n");
    int cholCheck = 0;
    double LyapTol=1e-06;
    try 
      {
      // make input matrices
      int s = mxGetM(prhs[0]);
      GeneralMatrix G(mxGetPr(prhs[0]), mxGetM(prhs[0]), mxGetN(prhs[0]));
      GeneralMatrix V(mxGetPr(prhs[1]), mxGetM(prhs[1]), mxGetN(prhs[1]));

      // create output
      plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]),mxGetN(prhs[0]), mxREAL);
      GeneralMatrix X(mxGetPr(plhs[0]),mxGetM(plhs[0]),mxGetN(plhs[0]));
      if (nrhs > 2)
        LyapTol = (double)mxGetScalar(prhs[2]);
      if (nrhs > 3)
        cholCheck = (int)mxGetScalar(prhs[3]);

#ifdef TIMING_LOOP
  for (int tt=0;tt<1000;++tt)
    {
#endif
      disclyap_fast(G, V, X, LyapTol, cholCheck);
#ifdef TIMING_LOOP
    }
    mexPrintf("disclyap_fast: finished 1000 loops");
#endif
      
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
