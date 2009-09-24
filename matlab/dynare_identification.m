function [pdraws, idemodel, idemoments] = dynare_identification(iload)

% main 
%
% Copyright (C) 2008 Dynare Team
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

global M_ options_ oo_ bayestopt_ estim_params_

if nargin==0 | isempty(iload),
  iload=0;
end

options_ = set_default_option(options_,'datafile',[]);
options_.mode_compute = 0;
[data,rawdata]=dynare_estimation_init([],1);
% computes a first linear solution to set up various variables
dynare_resolve;

options_.prior_mc=2000;

SampleSize = options_.prior_mc;

% results = prior_sampler(0,M_,bayestopt_,options_,oo_);

prior_draw(1,bayestopt_);
IdentifDirectoryName = CheckPath('identification');

indx = estim_params_.param_vals(:,1);
indexo=[];
if ~isempty(estim_params_.var_exo)
  indexo = estim_params_.var_exo(:,1);
end
useautocorr = 1;
nlags = 3;
nparam = length(bayestopt_.name);

if iload ==0, 
  
iteration = 0;
loop_indx = 0;

h = waitbar(0,'Monte Carlo identification checks ...');

while iteration < SampleSize,
  loop_indx = loop_indx+1;
  params = prior_draw();
  set_all_parameters(params);
  [A,B,ys,info]=dynare_resolve;

  
  if info(1)==0,
    iteration = iteration + 1;
    tau=[vec(A); vech(B*M_.Sigma_e*B')];
    [JJ, H, GAM] = getJJ(A, B, M_,oo_,options_,0,indx,indexo,bayestopt_.mf2,nlags,useautocorr);  
    siJ = abs(JJ(find(GAM),:).*(1./GAM(find(GAM))*params));
    siH = abs(H(find(abs(tau)>1.e-10),:).*(1./tau(find(abs(tau)>1.e-10))*params));
    stock_params(iteration,:) = params;
    if iteration ==1,
      siJmean = siJ./SampleSize;
      siHmean = siH./SampleSize;
    else
      siJmean = siJ./SampleSize+siJmean;
      siHmean = siH./SampleSize+siHmean;
    end
    pdraws(iteration,:) = params';
    [idemodel.Mco(:,iteration), idemoments.Mco(:,iteration), ...
      idemodel.Pco(:,:,iteration), idemoments.Pco(:,:,iteration), ...
      idemodel.cond(iteration), idemoments.cond(iteration), ...
      idemodel.ee(:,iteration), idemoments.ee(:,iteration), ...
      idemodel.ind(:,iteration), idemoments.ind(:,iteration), ...
      idemodel.indno{iteration}, idemoments.indno{iteration}] = ...
      identification_checks(H,JJ, bayestopt_);   
    
      waitbar(iteration/SampleSize,h)
  end
end
siHmean = siHmean./(max(siHmean')'*ones(size(params)));
siJmean = siJmean./(max(siJmean')'*ones(size(params)));

close(h)


save([IdentifDirectoryName '/' M_.fname '_identif'], 'pdraws', 'idemodel', 'idemoments', ...
  'siHmean', 'siJmean', 'stock_params')
else
load([IdentifDirectoryName '/' M_.fname '_identif'], 'pdraws', 'idemodel', 'idemoments', ...
  'siHmean', 'siJmean', 'stock_params')
end  

disp_identification(pdraws, idemodel, idemoments)

figure,
myboxplot(siHmean)
set(gca,'ylim',[0 1])
set(gca,'xticklabel','')
for ip=1:nparam,
  text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Sensitivity in the model')

figure,
myboxplot(siJmean)
set(gca,'ylim',[0 1])
set(gca,'xticklabel','')
for ip=1:nparam,
  text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Sensitivity in the moments')

figure,
myboxplot(idemodel.Mco')
set(gca,'ylim',[0 1])
set(gca,'xticklabel','')
for ip=1:nparam,
  text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Multicollinearity in the model')

figure,
myboxplot(idemoments.Mco')
set(gca,'ylim',[0 1])
set(gca,'xticklabel','')
for ip=1:nparam,
  text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Multicollinearity in the moments')
