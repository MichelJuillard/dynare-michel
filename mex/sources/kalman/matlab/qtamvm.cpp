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

/****************************************************** 
%
% This provides an interface to QT f90 library by Andrea Pagano 
% to multiply Quasi trinagular matrix (T) with a vector a
%
% function [a] = qtmvm(QT,a)
%
%  1. T1 = QT2T(QT;n) and Ld = QT2Ld(QT;n);
%  2. Ta = LdV(Ld;a;n)+TV(T1;a;n).
%
% INPUTS
%    T                      [double]    mm*mm transition matrix of the state equation.
%    a                      [double]    mm state vector.
%
% OUTPUTS
%    Tinverse               [double]    mm*mm transition matrix of the state equation.
**********************************************************/

#include "qt.h"
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
    
    if (nrhs < 2 )
      mexErrMsgTxt("Must have  min 2 input parameters.\n");
    if (nlhs < 1 )
      mexErrMsgTxt("Must have  min 1 output parameters.\n")
      
      ;
    //int start = 1; // default start of likelihood calculation
    // test for univariate case
    try 
      {
      // make input matrices
      int n=mxGetM(prhs[0]);

      double *T1, *Ld, *TV ;     
      double * QT=mxGetPr(prhs[0]);
      //      ConstGeneralMatrix a (mxGetPr(prhs[1]), mxGetM(prhs[1]), mxGetN(prhs[1]));
      double *a =(mxGetPr(prhs[1]));
      
      // create output and upload output data
      plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[1]),1, mxREAL);// mxGetM(prhs[1]), mxREAL);
      //      double * mxinv= mxGetPr(plhs[0]);
      //      GeneralMatrix Ta(mxGetPr(plhs[0]), mxGetM(prhs[1]), mxGetN(prhs[1]));
      double * Ta =mxGetPr(plhs[0]);
      //    Tinv.unit();
      //      Ta.zeros();

      T1=(double *)mxCalloc(n*n, sizeof(double));
      Ld=(double *)mxCalloc(n*n,sizeof(double));
      TV=(double *)mxCalloc(n, sizeof(double));

      
#ifdef TIMING_LOOP
      int loops=1;//000;
      if (nrhs >2 )
        loops = (int)mxGetScalar(prhs[2]);
      for (int tt=0;tt<loops;++tt)
        {
#endif
#ifdef DEBUG	
      //  QT.print();
#endif
        // 1. T1 = QT2T(QT;n) and Ld = QT2Ld(QT;n);
     //   double *T1, *Ld, *dTa;//, dT1=-7.77;
//  mexPrintf("start dT1 = %f\n", dT1);
//        dT1 = qt2t_(QT.base() ,&n) ;
// mexPrintf("end dT1 = %f\n", dT1);
        //T1=&dT1;

        qt2t_(T1, QT ,&n) ;
//        T1=*T1p;
//        GeneralMatrix T1gm(T1,n,n);
//        Ld = qt2ld_(QT.base(),&n);
        qt2ld_(Ld , QT,&n);
        // 2. Ta = LdV(Ld;a;n)+TV(T1;a;n).
//        dTa = ldv_(Ld,a.base() ,&n);
        //Vector Ta( n);
        ldv_(Ta, Ld,a ,&n);
//        Ta2 = tv_(T1,a.base(),&n);
        tv_(TV, T1 ,a,&n);
//        Ta.add(1.0,ConstVector(Ta2.base(), n));
//        Ta.add(1.0,TV);
          for (int j=0; j<n;++j)
            Ta[j]+=TV[j];


#ifdef TIMING_LOOP
        }
      mexPrintf("QT array mvm: finished: %d loops\n",loops);
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
      mxFree(T1);
      mxFree(Ld);
      mxFree(TV);
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
