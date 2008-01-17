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
% part of DYNARE, copyright Dynare Team (2006-2008)
% Gnu Public License.

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
