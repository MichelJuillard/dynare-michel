function make_y_()
% function make_y_
% forms oo_.endo_simul as guess values for deterministic simulations
%  
% INPUTS
%   ...
% OUTPUTS
%   ...
% ALGORITHM
%   ...
% SPECIAL REQUIREMENTS
%   none
%  

% Copyright (C) 1996-2012 Dynare Team
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

global M_ options_ oo_ ys0_ 

if options_.steadystate_flag
    [oo_.steady_state,M_.params,check] = ...
        evaluate_steady_state_file(oo_.steady_state,oo_.exo_steady_state,M_, ...
                                   options_);
end

if isempty(oo_.steady_state)
    oo_.steady_state = zeros(M_.endo_nbr,1);
end

if isempty(M_.endo_histval)
    if isempty(ys0_)
        oo_.endo_simul = [oo_.steady_state*ones(1,M_.maximum_lag+options_.periods+M_.maximum_lead)];
    else
        oo_.endo_simul = [ys0_*ones(1,M_.maximum_lag) oo_.steady_state*ones(1,options_.periods+M_.maximum_lead)];
    end
else
    if ~isempty(ys0_)
        error('histval and endval cannot be used simultaneously')
    end
    oo_.endo_simul = [M_.endo_histval ...
                      oo_.steady_state*ones(1,options_.periods+M_.maximum_lead)];
end
