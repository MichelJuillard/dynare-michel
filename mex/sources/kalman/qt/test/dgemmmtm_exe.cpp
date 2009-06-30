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
%
% use:
%     dgemmmtm_exe QTt_file SS_file size [loops - if enabled]
%
% NOTE: due to fortran matrix orientation, input matrices need to be passed 
% as transposed so QTt instead QT
%
%
% INPUTS
%    QT                      [double]    mm*mm transition matrix of the state equation.
%    SS                      [double]    mm*mm state cov matrix.
%
% OUTPUTS
%    TSTt update             [double]    mm*mm state cov matrix updated.
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
      double *TSTt, *TS ;
      AsciiNumberArray QT, SS;
      QT.GetMX(argv[1],n,n);
      SS.GetMX(argv[2],n,n);
      const double alpha=1.0;
      const double beta=0.0;
      // create output and upload output data
      TS=(double *)calloc(n*n, sizeof(double));
      TSTt=(double *)calloc(n*n, sizeof(double));

      
#ifdef TIMING_LOOP
      int loops=1;//000;
      if (argc >3 )
        loops = atoi(argv[4]);
      for (int tt=0;tt<loops;++tt)
        {
#endif

/* DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
* void BLAS_dgemm(BLCHAR transa, BLCHAR transb, CONST_BLINT m, CONST_BLINT n,
*					CONST_BLINT k, CONST_BLDOU alpha, CONST_BLDOU a, CONST_BLINT lda,
*					CONST_BLDOU b, CONST_BLINT ldb, CONST_BLDOU beta,
*					BLDOU c, CONST_BLINT ldc);
*
* SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*/
	  BLAS_dgemm("N", "N",  &n,  &n, &n, &alpha, QT.data, &n, SS.data, &n, &beta, TS, &n);

	  BLAS_dgemm("N", "T", &n,  &n, &n, &alpha, TS, &n, QT.data, &n, &beta, TSTt, &n);

#ifdef TIMING_LOOP
        }
      printf("QT array mvm: finished: %d loops\n",loops);
#endif
      // create output and upload output data
      WriteMX(argv[2], TSTt,n,n);
      free(TSTt);
      free(TS);
      } 
    catch (std::exception e) 
      {
      std::cout <<"Error" << std::endl;
      }
    
  }; //main
