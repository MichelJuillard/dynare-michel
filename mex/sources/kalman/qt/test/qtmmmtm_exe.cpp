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
% This provides an interface to BLAS dgemm  library function  
% to multiply Quasi trinagular matrix (T) with another matrix
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

#include "qt.h"
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
      double *TSTt ;
      AsciiNumberArray QT, SS;
      QT.GetMX(argv[1],n,n);
      SS.GetMX(argv[2],n,n);
      // create output and upload output data
      TSTt=(double *)calloc(n*n, sizeof(double));

      
#ifdef TIMING_LOOP
      int loops=1;//000;
      if (argc >3 )
        loops = atoi(argv[4]);
      for (int tt=0;tt<loops;++tt)
        {
#endif

/* qtsqtt_  performs one of the matrix-matrix operations
*
*     C := QT*S*QT'
*/

      qtsqtt_(TSTt, QT.data, SS.data, &n);


#ifdef TIMING_LOOP
        }
      printf("QT array mvm: finished: %d loops\n",loops);
#endif
      // create output and upload output data
      WriteMX(argv[2], TSTt,n,n);
      free(TSTt);
      } 
    catch (std::exception e) 
      {
      std::cout <<"Error" << std::endl;
      }
    
  }; //main
