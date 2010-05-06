function r = ep_residuals(x, y, ix, iy, s)
% Inversion of the extended path simulation approach. This routine computes the innovations needed to
% reproduce the time path of a subset of endogenous variables.    
%    
% INPUTS
%  o x    [double]   n*1 vector, time t innovations.    
%  o y    [double]   n*1 vector, time t restricted endogenous variables.
%  o ix   [integer]  index of control innovations in the full vector of innovations.
%  o iy   [integer]  index of controlled variables in the full vector of endogenous variables.    
%  o s    [double]   m*1 vector, endogenous variables at time t-1. 
%    
% 
% OUTPUTS
%  o r    [double]  n*1 vector of residuals.
%    
% ALGORITHM
%  
% SPECIAL REQUIREMENTS

% Copyright (C) 2010 Dynare Team.
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

global oo_

weight = 1.0;
tdx = M_.maximum_lag+1;

x = exp(transpose(x));

oo_.exo_simul(tdx,ix) = x;
exogenous_variables = zeros(size(oo_.exo_simul));
exogenous_variables(tdx,ix) = x;
initial_path = simult_(oo_.steady_state,dr,exogenous_variables,1);
oo_.endo_simul = weight*initial_path(:,1:end-1) + (1-weight)*oo_.endo_simul;

    
    