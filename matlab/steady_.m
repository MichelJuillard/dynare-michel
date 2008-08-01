
function steady_()

% function steady_()
% computes the steady state 
%  
% INPUTS
%   none
%  
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2007 Dynare Team
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

  global M_ oo_ it_ options_
  
  
  if options_.steadystate_flag
    [oo_.steady_state,check] = feval([M_.fname '_steadystate'],...
				     oo_.steady_state,...
				     [oo_.exo_steady_state; ...
		                      oo_.exo_det_steady_state]);
    % Check if the steady state obtained from the _steadystate file is a 
    % steady state.
    check1 = 0;
    if isfield(options_,'unit_root_vars')
      if isempty(options_.unit_root_vars)
	check1 = max(abs(feval([M_.fname '_static'],...
			       oo_.steady_state,...
			       [oo_.exo_steady_state; ...
		    oo_.exo_det_steady_state], M_.params))) > options_.dynatol ;
	if check1
	  error(['The seadystate values returned by ' M_.fname ...
		'_steadystate.m don''t solve the static model!' ])
	end
      end
    end
    if ~isempty(options_.steadystate_partial)
      ssvar = options_.steadystate_partial.ssvar;
      nov   = length(ssvar);
      indv  = zeros(nov,1);
      for i = 1:nov
	indv(i) = strmatch(ssvar(i),M_.endo_names,'exact');
      end
      [oo_.steady_state,check] = dynare_solve('restricted_steadystate',...
					      oo_.steady_state(indv),...
					      options_.jacobian_flag, ...	    
					      [oo_.exo_steady_state;oo_.exo_det_steady_state],indv);
    end
  else
    [oo_.steady_state,check] = dynare_solve([M_.fname '_static'],...
				     oo_.steady_state,...
				     options_.jacobian_flag, ...	    
			             [oo_.exo_steady_state; ...
                    oo_.exo_det_steady_state], M_.params);
  end

  if check ~= 0
    error('STEADY: convergence problems')
  end