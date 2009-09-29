function [pdraws, TAU, GAM0, H, JJ] = dynare_identification(iload, pdraws0)

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


if nargin==2,
options_.prior_mc=size(pdraws0,1);
else
options_.prior_mc=2000;
end

SampleSize = options_.prior_mc;

% results = prior_sampler(0,M_,bayestopt_,options_,oo_);

prior_draw(1,bayestopt_);
if ~(exist('sylvester3mr','file')==2),

dynareroot = strrep(which('dynare'),'dynare.m','');
addpath([dynareroot 'gensylv'])
end

IdentifDirectoryName = CheckPath('identification');

indx = estim_params_.param_vals(:,1);
indexo=[];
if ~isempty(estim_params_.var_exo)
  indexo = estim_params_.var_exo(:,1);
end
useautocorr = 1;
nlags = 3;
nparam = length(bayestopt_.name);
    options_.ar=nlags;

MaxNumberOfBytes=options_.MaxNumberOfBytes;
           
             
if iload <=0, 
  
iteration = 0;
burnin_iteration = 0;
loop_indx = 0;
file_index = 0;
run_index = 0;

h = waitbar(0,'Monte Carlo identification checks ...');

while iteration < SampleSize,
  loop_indx = loop_indx+1;
  if nargin==2 & burnin_iteration>=50,
    params = pdraws0(iteration+1,:);
  else
    params = prior_draw();
  end
  set_all_parameters(params);
  [A,B,ys,info]=dynare_resolve;

  
  if info(1)==0,
    oo0=oo_;
%     [Aa,Bb] = kalman_transition_matrix(oo0.dr, ...
%        bayestopt_.restrict_var_list, ...
%        bayestopt_.restrict_columns, ...
%        bayestopt_.restrict_aux, M_.exo_nbr);
%     tau=[vec(Aa); vech(Bb*M_.Sigma_e*Bb')];
    tau=[vec(A); vech(B*M_.Sigma_e*B')];
    if burnin_iteration<50,
      burnin_iteration = burnin_iteration + 1;
      TAU(:,burnin_iteration)=tau;
      [gam,stationary_vars] = th_autocovariances(oo0.dr,bayestopt_.mfys,M_,options_);
      sdy = sqrt(diag(gam{1}));
      sy = sdy*sdy';
      if useautocorr,
        sy=sy-diag(diag(sy))+eye(length(sy));
        gam{1}=gam{1}./sy;
      else
        for j=1:nlags,
          gam{j+1}=gam{j+1}.*sy;
        end
      end
      dum = vech(gam{1});
      for j=1:nlags,
        dum = [dum; vec(gam{j+1})];
      end
      GAM(:,burnin_iteration)=dum;
    else
      iteration = iteration + 1;
      run_index = run_index + 1;
      if iteration==1,
        indJJ = (find(std(GAM')>1.e-10));
        indH = (find(std(TAU')>1.e-10));
        TAU = zeros(length(indH),SampleSize);
        GAM = zeros(length(indJJ),SampleSize);
        MAX_tau   = min(SampleSize,ceil(MaxNumberOfBytes/(length(indH)*nparam)/8));
        MAX_gam   = min(SampleSize,ceil(MaxNumberOfBytes/(length(indJJ)*nparam)/8));
        stoH = zeros([length(indH),nparam,MAX_tau]);
        stoJJ = zeros([length(indJJ),nparam,MAX_tau]);
        delete([IdentifDirectoryName '/' M_.fname '_identif_*.mat'])
      end
    end

    if iteration,
      TAU(:,iteration)=tau(indH);
      [JJ, H, gam] = getJJ(A, B, M_,oo0,options_,0,indx,indexo,bayestopt_.mf2,nlags,useautocorr);
      GAM(:,iteration)=gam(indJJ);
      stoH(:,:,run_index) = H(indH,:);
      stoJJ(:,:,run_index) = JJ(indJJ,:);
      % use relative changes
      siJ = abs(JJ(indJJ,:).*(1./gam(indJJ)*params));
      siH = abs(H(indH,:).*(1./tau(indH)*params));
      % use prior uncertainty
      siJ = abs(JJ(indJJ,:));
      siH = abs(H(indH,:));
%       siJ = abs(JJ(indJJ,:).*(ones(length(indJJ),1)*bayestopt_.p2'));
%       siH = abs(H(indH,:).*(ones(length(indH),1)*bayestopt_.p2'));
%       siJ = abs(JJ(indJJ,:).*(1./mGAM'*bayestopt_.p2'));
%       siH = abs(H(indH,:).*(1./mTAU'*bayestopt_.p2'));

      if iteration ==1,
        siJmean = siJ./SampleSize;
        siHmean = siH./SampleSize;
      else
        siJmean = siJ./SampleSize+siJmean;
        siHmean = siH./SampleSize+siHmean;
      end
      pdraws(iteration,:) = params;
      [idemodel.Mco(:,iteration), idemoments.Mco(:,iteration), ...
        idemodel.Pco(:,:,iteration), idemoments.Pco(:,:,iteration), ...
        idemodel.cond(iteration), idemoments.cond(iteration), ...
        idemodel.ee(:,:,iteration), idemoments.ee(:,:,iteration), ...
        idemodel.ind(:,iteration), idemoments.ind(:,iteration), ...
        idemodel.indno{iteration}, idemoments.indno{iteration}] = ...
        identification_checks(H(indH,:),JJ(indJJ,:), bayestopt_);
      if run_index==MAX_tau | iteration==SampleSize,
        file_index = file_index + 1;
        if run_index<MAX_tau,
      stoH = stoH(:,:,1:run_index);
      stoJJ = stoJJ(:,:,1:run_index);
        end          
      save([IdentifDirectoryName '/' M_.fname '_identif_' int2str(file_index)], 'stoH', 'stoJJ')
      run_index = 0;
        
      end

      waitbar(iteration/SampleSize,h)
    end
  end
end

siJmean = siJmean.*(ones(length(indJJ),1)*std(pdraws));
siHmean = siHmean.*(ones(length(indH),1)*std(pdraws));

siHmean = siHmean./(max(siHmean')'*ones(size(params)));
siJmean = siJmean./(max(siJmean')'*ones(size(params)));

close(h)


save([IdentifDirectoryName '/' M_.fname '_identif'], 'pdraws', 'idemodel', 'idemoments', ...
  'siHmean', 'siJmean', 'TAU', 'GAM')
else
load([IdentifDirectoryName '/' M_.fname '_identif'], 'pdraws', 'idemodel', 'idemoments', ...
  'siHmean', 'siJmean', 'TAU', 'GAM')
options_.prior_mc=size(pdraws,1);
SampleSize = options_.prior_mc;

end  

if nargout>3 & iload,
  filnam = dir([IdentifDirectoryName '/' M_.fname '_identif_*.mat']);
  H=[];
  JJ = [];
  for j=1:length(filnam),
    load([IdentifDirectoryName '/' M_.fname '_identif_',int2str(j),'.mat']);
    H = cat(3,H, stoH(:,abs(iload),:));
    JJ = cat(3,JJ, stoJJ(:,abs(iload),:));

  end
end

mTAU = mean(TAU');
mGAM = mean(GAM');
sTAU = std(TAU');
sGAM = std(GAM');
if nargout>=3,
  GAM0=GAM;
end
if useautocorr,
  idiag = find(vech(eye(size(options_.varobs,1))));
  GAM(idiag,:) = GAM(idiag,:)./(sGAM(idiag)'*ones(1,SampleSize));
%   siJmean(idiag,:) = siJmean(idiag,:)./(sGAM(idiag)'*ones(1,nparam));
%   siJmean = siJmean./(max(siJmean')'*ones(size(params)));
end

[pcc, dd] = eig(cov(GAM'));
[latent, isort] = sort(-diag(dd));
latent = -latent;
pcc=pcc(:,isort);
siPCA = (siJmean'*abs(pcc')).^2';
siPCA = siPCA./(max(siPCA')'*ones(1,nparam)).*(latent*ones(1,nparam));
siPCA = sum(siPCA,1);
siPCA = siPCA./max(siPCA);

[pcc, dd] = eig(corrcoef(GAM'));
[latent, isort] = sort(-diag(dd));
latent = -latent;
pcc=pcc(:,isort);
siPCA2 = (siJmean'*abs(pcc')).^2';
siPCA2 = siPCA2./(max(siPCA2')'*ones(1,nparam)).*(latent*ones(1,nparam));
siPCA2 = sum(siPCA2,1);
siPCA2 = siPCA2./max(siPCA2);

[pcc, dd] = eig(cov(TAU'));
[latent, isort] = sort(-diag(dd));
latent = -latent;
pcc=pcc(:,isort);
siHPCA = (siHmean'*abs(pcc')).^2';
siHPCA = siHPCA./(max(siHPCA')'*ones(1,nparam)).*(latent*ones(1,nparam));
siHPCA = sum(siHPCA,1);
siHPCA = siHPCA./max(siHPCA);

[pcc, dd] = eig(corrcoef(TAU'));
[latent, isort] = sort(-diag(dd));
latent = -latent;
pcc=pcc(:,isort);
siHPCA2 = (siHmean'*abs(pcc')).^2';
siHPCA2 = siHPCA2./(max(siHPCA2')'*ones(1,nparam)).*(latent*ones(1,nparam));
siHPCA2 = sum(siHPCA2,1);
siHPCA2 = siHPCA2./max(siHPCA2);


disp_identification(pdraws, idemodel, idemoments)

% figure,
% % myboxplot(siPCA(1:(max(find(cumsum(latent)./length(indJJ)<0.99))+1),:))
% subplot(221)
% bar(siHPCA)
% % set(gca,'ylim',[0 1])
% set(gca,'xticklabel','')
% set(gca,'xlim',[0.5 nparam+0.5])
% for ip=1:nparam,
%   text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
% end
% title('Sensitivity in TAU''s PCA')
% 
% subplot(222)
% % myboxplot(siPCA(1:(max(find(cumsum(latent)./length(indJJ)<0.99))+1),:))
% bar(siHPCA2)
% % set(gca,'ylim',[0 1])
% set(gca,'xticklabel','')
% set(gca,'xlim',[0.5 nparam+0.5])
% for ip=1:nparam,
%   text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
% end
% title('Sensitivity in standardized TAU''s PCA')
% 
% 
% subplot(223)
% % myboxplot(siPCA(1:(max(find(cumsum(latent)./length(indJJ)<0.99))+1),:))
% bar(siPCA)
% % set(gca,'ylim',[0 1])
% set(gca,'xticklabel','')
% set(gca,'xlim',[0.5 nparam+0.5])
% for ip=1:nparam,
%   text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
% end
% title('Sensitivity in moments'' PCA')
% 
% subplot(224)
% % myboxplot(siPCA(1:(max(find(cumsum(latent)./length(indJJ)<0.99))+1),:))
% bar(siPCA2)
% % set(gca,'ylim',[0 1])
% set(gca,'xticklabel','')
% set(gca,'xlim',[0.5 nparam+0.5])
% for ip=1:nparam,
%   text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
% end
% title('Sensitivity in standardized moments'' PCA')

figure,
subplot(221)
myboxplot(siHmean)
set(gca,'ylim',[0 1])
set(gca,'xticklabel','')
for ip=1:nparam,
  text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Sensitivity in the model')

subplot(222)
myboxplot(siJmean)
set(gca,'ylim',[0 1])
set(gca,'xticklabel','')
for ip=1:nparam,
  text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Sensitivity in the moments')

subplot(223)
myboxplot(idemodel.Mco')
set(gca,'ylim',[0 1])
set(gca,'xticklabel','')
for ip=1:nparam,
  text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Multicollinearity in the model')

subplot(224)
myboxplot(idemoments.Mco')
set(gca,'ylim',[0 1])
set(gca,'xticklabel','')
for ip=1:nparam,
  text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Multicollinearity in the moments')
