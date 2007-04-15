function xparam = get_posterior_parameters(type)
% Selects (estimated) parameters (posterior mode or posterior mean).
% 
% INPUTS 
%   o type       [char]     = 'mode' or 'mean'.
%  
% OUTPUTS 
%   o xparam  
%
%
% ALGORITHM 
%   None.       
%
% SPECIAL REQUIREMENTS
%   None.
%  
%  
% part of DYNARE, copyright S. Adjemian, M. Juillard (2006)
% Gnu Public License.
global estim_params_ oo_ options_ M_ 
  
nvx = estim_params_.nvx;
nvn = estim_params_.nvn;
ncx = estim_params_.ncx;
ncn = estim_params_.ncn;
np  = estim_params_.np;

xparam = zeros(nvx+nvn+ncx+ncn+np,1);

m = 1;
for i=1:nvx
    k1 = estim_params_.var_exo(i,1);
    name1 = deblank(M_.exo_names(k1,:));
    xparam(m) = eval(['oo_.posterior_' type '.shocks_std.' name1]);
    M_.Sigma_e(k1,k1) = xparam(m)^2;
    m = m+1;
end

for i=1:nvn
    k1 = estim_params_.var_endo(i,1);
    name1 = deblank(options_.varobs(k1,:));
    xparam(m) = eval(['oo_.posterior_' type '.measurement_errors_std.' name1]);
    m = m+1;
end
  
for i=1:ncx
    k1 = estim_params_.corrx(i,1);
    k2 = estim_params_.corrx(i,2);
    name1 = deblank(M_.var_exo(k1,:));
    name2 = deblank(M_.var_exo(k2,:));
    xparam(m) = eval(['oo_.posterior_' type '.shocks_corr.' name1 '_' name2]);
    M_.Sigma_e(k1,k2) = xparam(m);
    M_.Sigma_e(k2,k1) = xparam(m);
    m = m+1;
end
  
for i=1:ncn
    k1 = estim_params_.corrn(i,1);
    k2 = estim_params_.corrn(i,2);
    name1 = deblank(options_.varobs(k1,:));
    name2 = deblank(options_.varobs(k2,:));
    xparam(m) = eval(['oo_.posterior_' type '.measurement_errors_corr.' name1 '_' name2]);
    m = m+1;
end

FirstDeep = m;

for i=1:np
    name1 = deblank(M_.param_names(estim_params_.param_vals(i,1),:));
    xparam(m) = eval(['oo_.posterior_' type '.parameters.' name1]);
    assignin('base',name1,xparam(m));% Useless with version 4 (except maybe for users)
    m = m+1;
end

if np
    M_.params(estim_params_.param_vals(:,1)) = xparam(FirstDeep:end);
end