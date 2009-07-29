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
% function X=disclyap_fast(G,V,tol,ch)
% 
% Solve the discrete Lyapunov Equation 
% X=G*X*G'+V 
% Using the Doubling Algorithm 
%
% INPUT: 
%   G and V - square General matrices of same size
%   tol     - double tollerance level
%   flag_ch - integer flag: if 1 check if the result is positive 
%             definite and generate an error message if it is not 
% OUTPUT:
%   on completion V - square General matrice contains solution
% 
% based on work of Joe Pearlman and Alejandro Justiniano 
% 3/5/2005 
% C++ version 28/07/09 by Dynare team
****************************************************************/
#include "ts_exception.h"
#include "cppblas.h"
#include "GeneralMatrix.h"
//#include "Vector.h"
#include "SylvException.h"
#include "utils.h"
#include "mex.h"

void disclyap_fast(const GeneralMatrix &G, const GeneralMatrix & V, GeneralMatrix &X, double tol = 1e-16, int flag_ch=0)
  {
  /**
  if nargin == 2 | isempty( ch ) == 1 
  flag_ch = 0; 
  else 
  flag_ch = 1; 
  end 
  **/
  //P0=V; 
  GeneralMatrix P0(V);
  //A0=G; 
  GeneralMatrix A0(G);
  
  //n=size(G,1); 
  int n= A0.numCols();
  const double alpha=1.0;
  const double half=0.5;
  const double neg_alpha=-1.0;
  const double omega=0.0;
 
  GeneralMatrix A1(n,n);
  GeneralMatrix Ptmp(n,n);
  GeneralMatrix P1(P0);
  GeneralMatrix I(n,n);
  I.unit();  
  bool matd=true; 
  while (matd ) // matrix diff > tol
    {
    //P1=P0+A0*P0*A0'; 
    // first step Ptmp=P0*A0'; 
    // DGEMM: C := alpha*op( A )*op( B ) + beta*C,
		BLAS_dgemm("N", "T", &n, &n, &n, &alpha, P0.base(), &n,
				   A0.base(), &n, &omega, Ptmp.base(), &n); 
    // P1=P0+A0*Ptmp; 
		BLAS_dgemm("N", "N", &n, &n, &n, &alpha, A0.base(), &n,
				   Ptmp.base(), &n, &alpha, P1.base(), &n); 
    // A1=A0*A0;  
    // A0=A1 (=A0*A0);  
    // A0.multRight(A0);
		BLAS_dgemm("N", "N", &n, &n, &n, &alpha, A0.base(), &n,
				   A0.base(), &n, &omega, A1.base(), &n); 

    // check if max( max( abs( P1 - P0 ) ) )>tol 
    matd=P0.isDiffSym(P1, tol);
    P0=P1; 
    A0=A1;
    }//end while
  
  //  X=P0=(P0+P0')/2; 
	BLAS_dgemm("T", "N", &n, &n, &n, &half, P1.base(), &n,
			 I.base(), &n, &half, P0.base(), &n); 
  X=P0;
  // Check that X is positive definite 
  if (flag_ch==1) 
    NormCholesky chol(P0);    
  }
