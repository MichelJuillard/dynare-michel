function [m0,s0] = compute_mh_covariance_matrix()

% function [m0,s0] = compute_mh_covariance_matrix()
% Estimation of the posterior covariance matrix and expectation. 
% 
% INPUTS 
%   None.
%  
% OUTPUTS
%   o  m0  [double]  (n*1) vector, posterior expectation of the parameters.
%   o  s0  [double]  (n*n) matrix, posterior covariance of the parameters 
%                    (computed from previous metropolis hastings).
%
% SPECIAL REQUIREMENTS
%   None.
%  
% part of DYNARE, copyright Dynare Team (2006-2008)
% Gnu Public License.

global M_ options_ estim_params_

n = estim_params_.np + ...
       estim_params_.nvn+ ...
       estim_params_.ncx+ ...
       estim_params_.ncn+ ...
       estim_params_.nvx;
nblck = options_.mh_nblck;

MhDirectoryName = CheckPath('metropolis');
load([ MhDirectoryName '/'  M_.fname '_mh_history'])

FirstMhFile = record.KeepedDraws.FirstMhFile;
FirstLine = record.KeepedDraws.FirstLine;
TotalNumberOfMhFiles = sum(record.MhDraws(:,2));

params = zeros(1,n); 
oldlogpo2 = -Inf;
offset = 0;
m0 = zeros(n,1); 
s0 = zeros(n,n);

for n = FirstMhFile:TotalNumberOfMhFiles
  for b = 1:nblck
    load([ MhDirectoryName '/' M_.fname '_mh' int2str(n) '_blck' int2str(b)],'x2','logpo2'); 
    [tmp,idx] = max(logpo2);
    if tmp>oldlogpo2
        oldlogpo2 = tmp;
        params = x2(idx,:);
    end
    [m0,s0,offset] = recursive_moments(m0,s0,x2(FirstLine,:),offset);
  end
  FirstLine = 1;
end

xparam1 = params';
hh = inv(s0);
fval = oldlogpo2;

save([M_.fname '_mh_mode'],'xparam1','hh','fval');