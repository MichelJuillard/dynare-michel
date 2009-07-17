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
% function [a] = qtamvm(QT,a)
%
% use:
%     qtamvm_exe QTt_file a_file size [loops - if enabled]
%
% NOTE: due to fortran matrix orientation, input matrices need to be passed 
% as transposed so QTt instead QT
%
%  1. T1 = QT2T(QT;n) and Ld = QT2Ld(QT;n);
%  2. Ta = LdV(Ld;a;n)+TV(T1;a;n).
%
% INPUTS
%    T                      [double]    mm*mm transition matrix of the state equation.
%    a                      [double]    mm state vector.
%
% OUTPUTS
%    a update               [double]    mm vector of the state equation.
%                                       as file: a_file_out
**********************************************************/

#include "qt.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ascii_array.h"
#include <iostream>
#include <stdexcept>
#include <malloc.h>

int main(int argc, char* argv[])
    {
    
    if (argc < 3 )
      {
      printf("Must have  min 2 input parameters.\n");
      exit(1);
      }
  
    try 
      {
      // make input matrices
      int n=atoi(argv[3]);
      double *T1, *Ld, *TV, * Ta ;
      AsciiNumberArray QT, a;
      QT.GetMX(argv[1],n,n);
      a.GetMX(argv[2],n,1);
      T1=(double *)calloc(n*n, sizeof(double));
      Ld=(double *)calloc(n*n,sizeof(double));
      TV=(double *)calloc(n, sizeof(double));
      // create output and upload output data
      Ta=(double *)calloc(n, sizeof(double));

      
#ifdef TIMING_LOOP
      int loops=1;//000;
      if (argc >3 )
        loops = atoi(argv[4]);
      for (int tt=0;tt<loops;++tt)
        {
#endif
#ifdef DEBUG	
      //  QT.print();
#endif
        // 1. T1 = QT2T(QT;n) and Ld = QT2Ld(QT;n);
        qt2t_(T1, QT.data ,&n) ;
        qt2ld_(Ld , QT.data,&n);
        // 2. Ta = LdV(Ld;a;n)+TV(T1;a;n).
        ldv_(Ta, Ld,a.data ,&n);
        tv_(TV, T1 ,a.data,&n);
//        Ta.add(1.0,TV);
          for (int j=0; j<n;++j)
            Ta[j]+=TV[j];


#ifdef TIMING_LOOP
        }
      printf("QT array mvm: finished: %d loops\n",loops);
#endif
      // create output and upload output data
      WriteMX(argv[2], Ta,n,1);
      free(T1);
      free(Ld);
      free(TV);
      free(Ta);
      } 
    catch (std::exception e) 
      {
      std::cout <<"Error" << std::endl;
      }
    
  }; //main
