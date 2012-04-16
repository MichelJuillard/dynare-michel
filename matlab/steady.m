function steady()
% function steady()
% computes and prints the steady state calculations
%  
% INPUTS
%   none
%  
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2010 Dynare Team
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

global M_ oo_ options_ ys0_ 

test_for_deep_parameters_calibration(M_);

if options_.steadystate_flag && options_.homotopy_mode
    error('STEADY: Can''t use homotopy when providing a steady state external file');
end


% Keep of a copy of M_.Sigma_e
Sigma_e = M_.Sigma_e;

% Set M_.Sigma_e=0 (we compute the *deterministic* steady state)
M_.Sigma_e = zeros(size(Sigma_e));


switch options_.homotopy_mode
  case 1
    homotopy1(options_.homotopy_values, options_.homotopy_steps);
  case 2
    homotopy2(options_.homotopy_values, options_.homotopy_steps);
  case 3
    homotopy3(options_.homotopy_values, options_.homotopy_steps);
end

[oo_.steady_state,M_.params,info] = steady_(M_,options_,oo_);

if info(1) && options_.steady.stop_on_error
    print_info(info,options_.noprint);
end

disp_steady_state(M_,oo_);

M_.Sigma_e = Sigma_e;
