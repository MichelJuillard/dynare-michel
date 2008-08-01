function prior_analysis(var_list)

% function prior_analysis(var_list)
% performs stochastic simulations for value of parameters drawn from
% the prior
%
% INPUTS:
%   var_list: list of variable names for which results are requested
%
% OUTPUTS:
%   none
%
% SPECIAL REQUIREMENTS
%   none.
%    

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

  global M_ options_ oo_ estim_params_ bayestopt_

  old_options = options_;
  if options_.replic < 100
  warning('Prior analysis requires at least 100 replications, preferably many more! options replic reset to 100')
    options_.replic = 100;
  end
  
  options_.order = 1;
  if options_.forecast
    forcst_unc(oo_.endo_simul(:,1:M_.maximum_lag),var_list);
  end
  
  options_ = old_options;
