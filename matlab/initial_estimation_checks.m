function initial_estimation_checks(xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations)

% function initial_estimation_checks(xparam1,gend,data)
% Checks data (complex values, ML evaluation, initial values, BK conditions,..)
% 
% INPUTS
%    xparam1: vector of parameters to be estimated
%    gend:    scalar specifying the number of observations
%    data:    matrix of data
%    
% OUTPUTS
%    none
%        
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2008 Dynare Team
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

  global dr1_test bayestopt_ estim_params_ options_ oo_ M_

  nv = size(data,1);
  if nv-size(options_.varobs,1)
    disp(' ')
    disp(['Declared number of observed variables = ' int2str(size(options_.varobs,1))])
    disp(['Number of variables in the database   = ' int2str(nv)])
    disp(' ')
    error(['Estimation can''t take place because the declared number of observed' ...
	   'variables doesn''t match the number of variables in the database.'])
  end
  if nv > M_.exo_nbr+estim_params_.nvn
    error(['Estimation can''t take place because there are less shocks than' ...
	   'observed variables'])
  end
  if (number_of_observations==gend*nv)% No missing observations...
      k = find(all(~isnan(data),2));
      r = rank(data(unique(k),:));
      if r < nv
          error(['Estimation can''t take place because the data are perfectly' ...
                 ' correlated']);
      end
  end
  
  if ~isempty(strmatch('dsge_prior_weight',M_.param_names))
    [fval,cost_flag,info] = DsgeVarLikelihood(xparam1,gend);
  else
    [fval,cost_flag,ys,trend_coeff,info] = DsgeLikelihood(xparam1,gend,data,data_index,number_of_observations,no_more_missing_observations);
  end

  % when their is an analytical steadystate, check that the values
  % returned by *_steadystate match with the static model
  if options_.steadystate_flag
    [oo_.steady_state,check] = feval([M_.fname '_steadystate'],...
				     oo_.steady_state,...
				     [oo_.exo_steady_state; ...
		    oo_.exo_det_steady_state]);
    % Check if the steady state obtained from the _steadystate file is a 
    % steady state.
    check1 = 0;
    if isfield(options_,'unit_root_vars') & options_.diffuse_filter == 0
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
  end
  
  if info(1) > 0
      disp('Error in computing likelihood for initial parameter values')
      print_info(info)
  end
  
  if any(abs(oo_.steady_state(bayestopt_.mfys))>1e-9) & (options_.prefilter==1) 
      disp(['You are trying to estimate a model with a non zero steady state for the observed endogenous'])
      disp(['variables using demeaned data!'])
      error('You should change something in your mod file...')
  end
  
  disp(['Initial value of the log posterior (or likelihood): ' num2str(-fval)]);
