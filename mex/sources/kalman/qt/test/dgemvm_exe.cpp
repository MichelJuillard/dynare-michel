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
% This provides an interface to BLAS dgemv f90 library function  
% to multiply Quasi trinagular matrix (T) with a vector a
%
% function [a] = dgevm(QT,a)
%
% use:
%     dgemvm_exe QTt_file a_file size [loops - if enabled]
%
% NOTE: due to fortran matrix orientation, input matrices need to be passed 
% as transposed so QTt instead QT
%
%  1. Ta = dgemv(Ld;a;n)+TV(T1;a;n).
%
% INPUTS
%    T                      [double]    mm*mm transition matrix of the state equation.
%    a                      [double]    mm state vector.
%
% OUTPUTS
%    a update               [double]    mm vector of the state equation.
%                                       as file: a_file_out
**********************************************************/

#include "cppblas.h"
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
      double *Ta ;
      AsciiNumberArray QT, a;
      QT.GetMX(argv[1],n,n);
      a.GetMX(argv[2],n,1);
      const double alpha=1.0;
      const double beta=0.0;
      int inc =1;
      //alpha=1.0;beta=0.0;
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
//	void BLAS_dgemv(BLCHAR trans, CONST_BLINT m, CONST_BLINT n, CONST_BLDOU alpha,
//					CONST_BLDOU a, CONST_BLINT lda, CONST_BLDOU x, CONST_BLINT incx,
//					CONST_BLDOU beta, BLDOU y, CONST_BLINT incy);

	  BLAS_dgemv("N",  &n,  &n, &alpha, QT.data, &n, a.data, &inc, &beta, Ta, &inc);


#ifdef TIMING_LOOP
        }
      printf("QT array mvm: finished: %d loops\n",loops);
#endif
      // create output and upload output data
      WriteMX(argv[2], Ta,n,1);
      free(Ta);
      } 
    catch (std::exception e) 
      {
      std::cout <<"Error" << std::endl;
      }
    
  }; //main
