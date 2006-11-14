function prior_analysis(var_list)
% function prior_analysis(var_list)
% performs stochastic simulations for value of parameters drawn from
% the prior
% INPUTS:
%   var_list: list of variable names for which results are requested
% OUTPUTS:
%   none
% ALGORITHM
%   uses antithetic draws for the shocks
%
% SPECIAL REQUIREMENTS
%   None.
%    
% part of DYNARE, copyright S. Adjemian, M. Juillard (2006)
% Gnu Public License.

  global M_ options_ oo_ estim_params_ bayestopt_
  
  old_options = options_;
  if options_.replic == 0;
    options_.replic = 1000;
    options_.order = 1;
  end
  
  if options_.forecast
    forcst_unc(oo_.endo_simul(:,1:M_.maximum_lag),var_list);
  end
  
  options_ = old_options;