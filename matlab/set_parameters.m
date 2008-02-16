function set_parameters(xparam1)

% function set_parameters(xparam1)
% Sets parameters value (except measurement errors)
% This is called for computations such as IRF and forecast
% when measurement errors aren't taken into account
% 
% INPUTS
%    xparam1:   vector of parameters to be estimated (initial values)
%    
% OUTPUTS
%    none
%        
% SPECIAL REQUIREMENTS
%    none
%  
% part of DYNARE, copyright Dynare Team (2003-2008)
% Gnu Public License.

  global estim_params_ M_
  
global estim_params_ M_
 
  nvx = estim_params_.nvx;
  ncx = estim_params_.ncx;
  nvn = estim_params_.nvn;
  ncn = estim_params_.ncn;
  np = estim_params_.np;
  Sigma_e = M_.Sigma_e;
  offset = 0;
 
  % stderrs of the exogenous shocks
  if nvx
    var_exo = estim_params_.var_exo;
    for i=1:nvx
      k = var_exo(i,1);
      Sigma_e(k,k) = xparam1(i)^2;
    end
  end
  % and update offset
  offset = offset + nvx + nvn;
 
  % correlations amonx shocks (ncx)
  if ncx
    corrx = estim_params_.corrx;
    for i=1:ncx
      k1 = corrx(i,1);
      k2 = corrx(i,2);
      Sigma_e(k1,k2) = xparam1(i+offset)*sqrt(Sigma_e(k1,k1)*Sigma_e(k2,k2));
      Sigma_e(k2,k1) = Sigma_e(k1,k2);
    end
  end
  % and update offset
  offset = offset + ncx + ncn;

  % structural parameters
  if np
    M_.params(estim_params_.param_vals(:,1)) = xparam1(offset+1:end);
  end
 
  M_.Sigma_e = Sigma_e;