function print_info(info)

% function print_info(info)
% Prints error messages
%
% INPUTS
%    info=1:         the model doesn't determine the current variables '...' uniquely
%    info=2:         MJDGGES returns the following error code'
%    info=3:         Blanchard Kahn conditions are not satisfied: no stable '...' equilibrium
%    info=4:         Blanchard Kahn conditions are not satisfied:'...' indeterminacy
%    info=5:         Blanchard Kahn conditions are not satisfied:'...' indeterminacy due to rank failure
%    info=20:        Impossible to find the steady state. Either the model' ...' doesn't have 
%                    a unique steady state of the guess values' ...' are too far from the solution
%    info=30:        Variance can't be computed
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2005-2007 Dynare Team
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

  global options_

  options_ = set_default_option(options_,'noprint',0);
  
  if ~options_.noprint
    switch info(1)
     case 1
      error(['The model doesn''t determine the current variables' ...
	     ' uniquely'])
     case 2
      error(['MJDGGES returns the following error code' ...
	     int2str(info(2))])
     case 3
      error(['Blanchard Kahn conditions are not satisfied: no stable' ...
	     ' equilibrium'])
     case 4
      error(['Blanchard Kahn conditions are not satisfied:' ...
	     ' indeterminacy'])
     case 5
      error(['Blanchard Kahn conditions are not satisfied:' ...
	     ' indeterminacy due to rank failure'])
     case 20
      error(['Impossible to find the steady state. Either the model' ...
	     ' doesn''t have a unique steady state of the guess values' ...
	     ' are too far from the solution']) 
     case 30
      error('Variance can''t be computed')
     otherwise
      error('This case shouldn''t happen. Contact the authors of Dynare')
    end
  end
  