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
% function X=disclyap_fast(G,V,ch)
% 
% Solve the discrete Lyapunov Equation 
% X=G*X*G'+V 
% Using the Doubling Algorithm 
%
% If ch is defined then the code will check if the resulting X 
% is positive definite and generate an error message if it is not 
% 
% based on work of Joe Pearlman and Alejandro Justiniano 
% 3/5/2005 
% C++ version 28/07/09 by Dynare team
****************************************************************/
#include "GeneralMatrix.h"
void disclyap_fast(const GeneralMatrix &G, const GeneralMatrix & V, GeneralMatrix &X, double tol, int ch);
