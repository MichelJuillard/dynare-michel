function initial_estimation_checks(xparam1,gend,data)

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
%  
% part of DYNARE, copyright Dynare Team (2003-2008)
% Gnu Public License.


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
  r = rank(data);
  if r < nv
    error(['Estimation can''t take place because the data are perfectly' ...
 	   ' correlated']);
  end

  if ~isempty(strmatch('dsge_prior_weight',M_.param_names))
    [fval,cost_flag,info] = DsgeVarLikelihood(xparam1,gend);
  else
    [fval,cost_flag,ys,trend_coeff,info] = DsgeLikelihood(xparam1,gend,data);
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
    if isfield(options_,'unit_root_vars')
      if isempty(options_.unit_root_vars)
	check1 = max(abs(feval([M_.fname '_static'],...
			       oo_.steady_state,...
			       [oo_.exo_steady_state; ...
		    oo_.exo_det_steady_state]))) > options_.dynatol ;
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

  disp(['Initial value of the log posterior (or likelihood): ' num2str(-fval)]);
