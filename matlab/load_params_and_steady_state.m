function load_params_and_steady_state(filename)
% function load_params_and_steady_state(filename)
%
% For all parameters, endogenous and exogenous variables, loads
% their value from a file created with save_params_and_steady_state.
% * for parameters, their value will be initialized as if they
%   had been calibrated in the .mod file
% * for endogenous and exogenous, their value will be initialized
%   as they would have been from an initval block
%
% INPUTS
%   filename:   where to load from the saved values
%  
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2008 Dynare Team
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

  global M_ oo_

  load(filename);
  
  if ~exist('stored_values')
    error('LOAD_PARAMS_AND_INITVAL: filename provided was probably not created by save_params_and_initval')
  end
  
  names = fieldnames(stored_values);
  
  for i = 1:size(names,1)
    field = names{i};
    j = strmatch(field, M_.param_names, 'exact');
    if ~isempty(j)
      M_.params(j) = stored_values.(field);
    else
      j = strmatch(field, M_.endo_names, 'exact');
      if ~isempty(j)
        oo_.steady_state(j) = stored_values.(field);
      else
        j = strmatch(field, M_.exo_names, 'exact');
        if ~isempty(j)
          oo_.exo_steady_state(j) = stored_values.(field);
        else
          warning(['LOAD_PARAMS_AND_INITVAL: Unknown symbol name: ', field])
        end
      end
    end
  end
