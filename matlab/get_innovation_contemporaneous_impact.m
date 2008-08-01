function B = get_innovation_contemporaneous_impact(type,info)

% function F = get_innovation_contemporaneous_impact('type')
% The approximated reduced form model is 
% 
%   Y^*_t = Z Y_t                   [Measure]
%   Y_t   = A*Y_{t-1} + B*E_t       [State]
%
% where Z is an p*m selection matrix (p<=m), Y^* is the p*1 vector of
% observable endogenous variables, Y is an m*1 vector of endogeneous
% variables, A is an m*m matrix, B is an m*r matrix (r<=m) and E an r*1 
% vector of structural innovations.   
% 
% The contemporaneous is return impact (on the observables) of an innovation is      
% given by F = Z*B. Matrix F is returned by this function.      
%
% INPUTS 
%   o type   [string]      "mode","mean"
%   o info   [integer]     if equal to 1, matrix B is saved in a mat file.
%  
% OUTPUTS 
%   o F (F is also saved in a file)
%
% SPECIAL REQUIREMENTS
%   This function needs to be run after the estimation of a model.

% Copyright (C) 2006-2008 Dynare Team
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

global oo_ M_ bayestopt_ options_

if nargin == 0
    type = 'mode';
end

if nargin == 1
    info = 0;
end

get_posterior_parameters(type);

[dr,info,M_,options_,oo_]=dr1(oo_.dr,0,M_,options_,oo_);

B(dr.order_var,M_.exo_names_orig_ord) = dr.ghu*sqrt(M_.Sigma_e);
B = B(bayestopt_.mfys,:);

save([M_.fname '_InnovImpact.mat'],'B');