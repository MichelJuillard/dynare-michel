function y = local_state_space_iteration_2(yhat,epsilon,ghx,ghu,constant,half_ghxx,half_ghuu,ghxu)
% Given an initial condition (y) and an innovation (epsilon), this
% routines computes the next value of the endogenous variables if the
% model is approximated by an order two taylor expansion around the
% deterministic steady state.
%
% INPUTS
%    yhat          [double]     n*1 vector, initial condition, where n is the number of state variables.
%    epsilon       [double]     q*1 vector, structural innovations.
%    ys            [double]     m*1 vector, steady state (the variables are ordered as in ghx) where m 
%                               is the number of elements in the union of the states and observed 
%                               variables.
%    ghx           [double]     m*n matrix, is a subset of dr.ghx we only consider the lines corresponding
%                               to the states and the observed variables.
%    ghu           [double]     m*q matrix, is a subset of dr.ghu 
%    constant      [double]     m*1 vector (steady state + second order correction).
%    half_ghxx     [double]     m*n² matrix, subset of .5*dr.ghxx. 
%    half_ghuu     [double]     m*q² matrix, subset of .5*dr.ghuu.
%    ghxu          [double]     m*nq matrix, subset of dr.ghxu.
%    index         [integer]    
%
% OUTPUTS
%    y             [double]     stochastic simulations results
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2009 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
       
    y    = constant + ghx*yhat + ghu*epsilon ...
           + A_times_B_kronecker_C(half_ghxx,yhat) ...
           + A_times_B_kronecker_C(half_ghuu,epsilon) ...
           + A_times_B_kronecker_C(ghxu,yhat,epsilon);