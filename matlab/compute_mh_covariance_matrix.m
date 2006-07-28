function covar = compute_mh_covariance_matrix()
% Estimation of the posterior covariance matrix. 
% 
% INPUTS 
%   None.
%  
% OUTPUTS 
%   o covar   [double] p*p matrix, posterior covariance of the estimated
%   parameters (computed from previous metropolis hastings).  
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
global M_ options_ estim_params_ oo_

npar = estim_params_.np+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.nvx;
nblck = options_.mh_nblck;

MhDirectoryName = CheckPath('metropolis');
load([ MhDirectoryName '/'  M_.fname '_mh_history'])

FirstMhFile = record.KeepedDraws.FirstMhFile;
FirstLine = record.KeepedDraws.FirstLine; ifil = FirstLine;
TotalNumberOfMhFiles = sum(record.MhDraws(:,2));
TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);
TODROP = floor(options_.mh_drop*TotalNumberOfMhDraws);

MU = zeros(1,npar);
SIGMA = zeros(npar,npar);

fprintf('MH: I''m computing the posterior mean... ');
for n = FirstMhFile:TotalNumberOfMhFiles
  for b = 1:nblck
    load([ MhDirectoryName '/' M_.fname '_mh' int2str(n) '_blck' int2str(b)],'x2','logpo2'); 
    MU = MU + sum(x2(ifil:end,:));
  end
  ifil = 1;
end
MU = MU/((TotalNumberOfMhDraws-TODROP)*nblck);
MU1 = repmat(MU,MAX_nruns,1);	
fprintf(' Done!\n');
fprintf('MH: I''m computing the posterior covariance matrix... ');
ifil = FirstLine;
for n = FirstMhFile:TotalNumberOfMhFiles
  for b = 1:nblck
    load([MhDirectoryName '/' M_.fname '_mh' int2str(n) '_blck' int2str(b)],'x2');
    x = x2(ifil:end,:)-MU1(1:size(x2(ifil:end,:),1),:);
    SIGMA = SIGMA + x'*x;
  end
  ifil = 1;
end
SIGMA =  SIGMA/((TotalNumberOfMhDraws-TODROP)*nblck);%<=== Variance of the parameters (ok!)
