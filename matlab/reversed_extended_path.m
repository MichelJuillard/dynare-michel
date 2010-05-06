function innovation_paths = reversed_extended_path(controlled_time_series, controlled_variable_names, control_innovation_names)
% Inversion of the extended path simulation approach. This routine computes the innovations needed to
% reproduce the time path of a subset of endogenous variables. The initial condition is teh deterministic
% steady state.   
%    
% INPUTS
%  o controlled_time_series           [double]  n*T matrix.    
%  o controlled_variable_names   [string]    n*1 matlab's cell. 
% o control_innovation_names     [string]    n*1 matlab's cell.  
%
% OUTPUTS
%  o innovations                        [double]  n*T matrix.
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

global M_ oo_ options_   

%% Initialization

% Compute the deterministic steady state.
steady_;

%  Compute the first order perturbation reduced form.
oldopt = options_;
options_.order = 1;
[dr,info]=resol(oo_.steady_state,0);
oo_.dr = dr; options_ = oldopt;

% Set-up oo_.exo_simul.
make_ex_; 

% Set-up oo_.endo_simul.
make_y_;

% Get indices of the controlled endogenous variables
n  = length(controlled_variable_names);
iy = NaN(n,1);
for k=1:n
    iy(k) = strmatch(controlled_variable_names{k},M_.endo_names,'exact');
end

% Get indices of the control innovations.
ix = NaN(n,1);
for k=1:n
    ix(k) = strmatch(control_innovation_names{k},M_.exo_names,'exact');
end

% Get the length of the sample.
T = size(controlled_time_series,2);

% Output initialization.
innovation_paths = zeros(n,T);

% Initialization of the perfect foresight model solver.
perfect_foresight_simulation();

