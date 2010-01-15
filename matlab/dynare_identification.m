function [pdraws, TAU, GAM, LRE, gp, H, JJ] = dynare_identification(options_ident, pdraws0)

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

options_ident = set_default_option(options_ident,'load_ident_files',0);
options_ident = set_default_option(options_ident,'useautocorr',1);
options_ident = set_default_option(options_ident,'ar',3);
options_ident = set_default_option(options_ident,'prior_mc',2000);
if nargin==2,
    options_ident.prior_mc=size(pdraws0,1);
end

iload = options_ident.load_ident_files;
nlags = options_ident.ar;
useautocorr = options_ident.useautocorr;
options_.ar=nlags;
options_.prior_mc = options_ident.prior_mc;
options_.options_ident = options_ident;


options_ = set_default_option(options_,'datafile',[]);
options_.mode_compute = 0;
[data,rawdata]=dynare_estimation_init([],1);
% computes a first linear solution to set up various variables



SampleSize = options_ident.prior_mc;

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

nparam = length(bayestopt_.name);

MaxNumberOfBytes=options_.MaxNumberOfBytes;


if iload <=0, 
    
    iteration = 0;
    burnin_iteration = 0;
    loop_indx = 0;
    file_index = 0;
    run_index = 0;

    h = waitbar(0,'Monte Carlo identification checks ...');
    [I,J]=find(M_.lead_lag_incidence');

    while iteration < SampleSize,
        loop_indx = loop_indx+1;
        if nargin==2,
            if burnin_iteration>=50,
                params = pdraws0(iteration+1,:);
            else
                params = pdraws0(burnin_iteration+1,:);
            end
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
            tau=[oo_.dr.ys(oo_.dr.order_var); vec(A); vech(B*M_.Sigma_e*B')];
            yy0=oo_.dr.ys(I);    
            [residual, g1 ] = feval([M_.fname,'_dynamic'],yy0, oo_.exo_steady_state', M_.params,1);    

            if burnin_iteration<50,
                burnin_iteration = burnin_iteration + 1;
                pdraws(burnin_iteration,:) = params;
                TAU(:,burnin_iteration)=tau;
                LRE(:,burnin_iteration)=vec(g1);      
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
                GAM(:,burnin_iteration)=[oo_.dr.ys(bayestopt_.mfys); dum];
            else
                iteration = iteration + 1;
                run_index = run_index + 1;
                if iteration==1,
                    indJJ = (find(std(GAM')>1.e-8));
                    indH = (find(std(TAU')>1.e-8));
                    indLRE = (find(std(LRE')>1.e-8));
                    TAU = zeros(length(indH),SampleSize);
                    GAM = zeros(length(indJJ),SampleSize);
                    LRE = zeros(length(indLRE),SampleSize);
                    MAX_tau   = min(SampleSize,ceil(MaxNumberOfBytes/(length(indH)*nparam)/8));
                    MAX_gam   = min(SampleSize,ceil(MaxNumberOfBytes/(length(indJJ)*nparam)/8));
                    stoH = zeros([length(indH),nparam,MAX_tau]);
                    stoJJ = zeros([length(indJJ),nparam,MAX_tau]);
                    delete([IdentifDirectoryName '/' M_.fname '_identif_*.mat'])
                end
            end

            if iteration,
                TAU(:,iteration)=tau(indH);
                vg1 = vec(g1);
                LRE(:,iteration)=vg1(indLRE);
                [JJ, H, gam, gp] = getJJ(A, B, M_,oo0,options_,0,indx,indexo,bayestopt_.mf2,nlags,useautocorr);      
                GAM(:,iteration)=gam(indJJ);
                stoLRE(:,:,run_index) = gp(indLRE,:);
                stoH(:,:,run_index) = H(indH,:);
                stoJJ(:,:,run_index) = JJ(indJJ,:);
                % use relative changes
                %       siJ = abs(JJ(indJJ,:).*(1./gam(indJJ)*params));
                %       siH = abs(H(indH,:).*(1./tau(indH)*params));
                % use prior uncertainty
                siJ = (JJ(indJJ,:));
                siH = (H(indH,:));
                siLRE = (gp(indLRE,:));
                %       siJ = abs(JJ(indJJ,:).*(ones(length(indJJ),1)*bayestopt_.p2'));
                %       siH = abs(H(indH,:).*(ones(length(indH),1)*bayestopt_.p2'));
                %       siJ = abs(JJ(indJJ,:).*(1./mGAM'*bayestopt_.p2'));
                %       siH = abs(H(indH,:).*(1./mTAU'*bayestopt_.p2'));

                if iteration ==1,
                    siJmean = abs(siJ)./SampleSize;        
                    siHmean = abs(siH)./SampleSize;        
                    siLREmean = abs(siLRE)./SampleSize;        
                    derJmean = (siJ)./SampleSize;        
                    derHmean = (siH)./SampleSize;        
                    derLREmean = (siLRE)./SampleSize;      
                else
                    siJmean = abs(siJ)./SampleSize+siJmean;        
                    siHmean = abs(siH)./SampleSize+siHmean;        
                    siLREmean = abs(siLRE)./SampleSize+siLREmean;        
                    derJmean = (siJ)./SampleSize+derJmean;        
                    derHmean = (siH)./SampleSize+derHmean;        
                    derLREmean = (siLRE)./SampleSize+derLREmean;      
                end
                pdraws(iteration,:) = params;
                [idemodel.Mco(:,iteration), idemoments.Mco(:,iteration), idelre.Mco(:,iteration), ...        
                 idemodel.Pco(:,:,iteration), idemoments.Pco(:,:,iteration), idelre.Pco(:,:,iteration), ...        
                 idemodel.cond(iteration), idemoments.cond(iteration), idelre.cond(iteration), ...        
                 idemodel.ee(:,:,iteration), idemoments.ee(:,:,iteration), idelre.ee(:,:,iteration), ...        
                 idemodel.ind(:,iteration), idemoments.ind(:,iteration), ...
                 idemodel.indno{iteration}, idemoments.indno{iteration}] = ...
                 identification_checks(H(indH,:),JJ(indJJ,:), gp(indLRE,:), bayestopt_);      
                if run_index==MAX_tau | iteration==SampleSize,
                    file_index = file_index + 1;
                    if run_index<MAX_tau,
                        stoH = stoH(:,:,1:run_index);
                        stoJJ = stoJJ(:,:,1:run_index);
                        stoLRE = stoLRE(:,:,1:run_index);        
                    end          
                    save([IdentifDirectoryName '/' M_.fname '_identif_' int2str(file_index)], 'stoH', 'stoJJ', 'stoLRE')      
                    run_index = 0;
                    
                end

                waitbar(iteration/SampleSize,h)
            end
        end
    end



    close(h)


    save([IdentifDirectoryName '/' M_.fname '_identif'], 'pdraws', 'idemodel', 'idemoments', 'idelre', 'indJJ', 'indH', 'indLRE', ...  
        'siHmean', 'siJmean', 'siLREmean', 'derHmean', 'derJmean', 'derLREmean', 'TAU', 'GAM', 'LRE')
else
    load([IdentifDirectoryName '/' M_.fname '_identif'], 'pdraws', 'idemodel', 'idemoments', 'idelre', 'indJJ', 'indH', 'indLRE', ...  
        'siHmean', 'siJmean', 'siLREmean', 'derHmean', 'derJmean', 'derLREmean', 'TAU', 'GAM', 'LRE')
        options_ident.prior_mc=size(pdraws,1);
        SampleSize = options_ident.prior_mc;  
        options_.options_ident = options_ident;
end  

offset = estim_params_.nvx;
offset = offset + estim_params_.nvn;
offset = offset + estim_params_.ncx;
offset = offset + estim_params_.ncn;

siJmean = siJmean.*(ones(length(indJJ),1)*std(pdraws));
siHmean = siHmean.*(ones(length(indH),1)*std(pdraws));
siLREmean = siLREmean.*(ones(length(indLRE),1)*std(pdraws(:, offset+1:end )));

derJmean = derJmean.*(ones(length(indJJ),1)*std(pdraws));
derHmean = derHmean.*(ones(length(indH),1)*std(pdraws));
derLREmean = derLREmean.*(ones(length(indLRE),1)*std(pdraws(:, offset+1:end )));

derHmean = abs(derHmean./(max(siHmean')'*ones(1,size(pdraws,2))));
derJmean = abs(derJmean./(max(siJmean')'*ones(1,size(pdraws,2))));
derLREmean = abs(derLREmean./(max(siLREmean')'*ones(1,estim_params_.np)));

siHmean = siHmean./(max(siHmean')'*ones(1,size(pdraws,2)));
siJmean = siJmean./(max(siJmean')'*ones(1,size(pdraws,2)));
siLREmean = siLREmean./(max(siLREmean')'*ones(1,estim_params_.np));

tstJmean = derJmean*0;
tstHmean = derHmean*0;
tstLREmean = derLREmean*0;
for j=1:nparam,
  indd = 1:length(siJmean(:,j));
tstJmean(indd,j) = abs(derJmean(indd,j))./siJmean(indd,j);
  indd = 1:length(siHmean(:,j));
tstHmean(indd,j) = abs(derHmean(indd,j))./siHmean(indd,j);
if j>offset
  indd = 1:length(siLREmean(:,j-offset));
tstLREmean(indd,j-offset) = abs(derLREmean(indd,j-offset))./siLREmean(indd,j-offset);
end  
end

if nargout>3 & iload,
    filnam = dir([IdentifDirectoryName '/' M_.fname '_identif_*.mat']);
    H=[];
    JJ = [];
    gp = [];  
    for j=1:length(filnam),
        load([IdentifDirectoryName '/' M_.fname '_identif_',int2str(j),'.mat']);
        H = cat(3,H, stoH(:,abs(iload),:));
        JJ = cat(3,JJ, stoJJ(:,abs(iload),:));
        gp = cat(3,gp, stoLRE(:,abs(iload),:));

    end
end

% mTAU = mean(TAU');
% mGAM = mean(GAM');
% sTAU = std(TAU');
% sGAM = std(GAM');
% if nargout>=3,
%   GAM0=GAM;
% end
% if useautocorr,
%   idiag = find(vech(eye(size(options_.varobs,1))));
%   GAM(idiag,:) = GAM(idiag,:)./(sGAM(idiag)'*ones(1,SampleSize));
% %   siJmean(idiag,:) = siJmean(idiag,:)./(sGAM(idiag)'*ones(1,nparam));
% %   siJmean = siJmean./(max(siJmean')'*ones(size(params)));
% end
% 
% [pcc, dd] = eig(cov(GAM'));
% [latent, isort] = sort(-diag(dd));
% latent = -latent;
% pcc=pcc(:,isort);
% siPCA = (siJmean'*abs(pcc')).^2';
% siPCA = siPCA./(max(siPCA')'*ones(1,nparam)).*(latent*ones(1,nparam));
% siPCA = sum(siPCA,1);
% siPCA = siPCA./max(siPCA);
% 
% [pcc, dd] = eig(corrcoef(GAM'));
% [latent, isort] = sort(-diag(dd));
% latent = -latent;
% pcc=pcc(:,isort);
% siPCA2 = (siJmean'*abs(pcc')).^2';
% siPCA2 = siPCA2./(max(siPCA2')'*ones(1,nparam)).*(latent*ones(1,nparam));
% siPCA2 = sum(siPCA2,1);
% siPCA2 = siPCA2./max(siPCA2);
% 
% [pcc, dd] = eig(cov(TAU'));
% [latent, isort] = sort(-diag(dd));
% latent = -latent;
% pcc=pcc(:,isort);
% siHPCA = (siHmean'*abs(pcc')).^2';
% siHPCA = siHPCA./(max(siHPCA')'*ones(1,nparam)).*(latent*ones(1,nparam));
% siHPCA = sum(siHPCA,1);
% siHPCA = siHPCA./max(siHPCA);
% 
% [pcc, dd] = eig(corrcoef(TAU'));
% [latent, isort] = sort(-diag(dd));
% latent = -latent;
% pcc=pcc(:,isort);
% siHPCA2 = (siHmean'*abs(pcc')).^2';
% siHPCA2 = siHPCA2./(max(siHPCA2')'*ones(1,nparam)).*(latent*ones(1,nparam));
% siHPCA2 = sum(siHPCA2,1);
% siHPCA2 = siHPCA2./max(siHPCA2);


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
subplot(231)
myboxplot(siLREmean)
set(gca,'ylim',[0 1.05])
set(gca,'xticklabel','')
for ip=1:estim_params_.np,
  text(ip,-0.02,deblank(M_.param_names(estim_params_.param_vals(ip,1),:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Sensitivity in the LRE model')

subplot(232)
myboxplot(siHmean)
set(gca,'ylim',[0 1.05])
set(gca,'xticklabel','')
for ip=1:nparam,
    text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Sensitivity in the model')

subplot(233)
myboxplot(siJmean)
set(gca,'ylim',[0 1.05])
set(gca,'xticklabel','')
for ip=1:nparam,
    text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Sensitivity in the moments')

subplot(234)
myboxplot(idelre.Mco')
set(gca,'ylim',[0 1])
set(gca,'xticklabel','')
for ip=1:estim_params_.np,
  text(ip,-0.02,deblank(M_.param_names(estim_params_.param_vals(ip,1),:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Multicollinearity in the LRE model')

subplot(235)
myboxplot(idemodel.Mco')
set(gca,'ylim',[0 1])
set(gca,'xticklabel','')
for ip=1:nparam,
    text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Multicollinearity in the model')

subplot(236)
myboxplot(idemoments.Mco')
set(gca,'ylim',[0 1])
set(gca,'xticklabel','')
for ip=1:nparam,
    text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Multicollinearity in the moments')
saveas(gcf,[IdentifDirectoryName,'/',M_.fname,'_ident'])
eval(['print -depsc2 ' IdentifDirectoryName '/' M_.fname '_ident']);
eval(['print -dpdf ' IdentifDirectoryName '/' M_.fname '_ident']);


figure,
subplot(221)
hist(log10(idemodel.cond))
title('log10 of Condition number in the model')
subplot(222)
hist(log10(idemoments.cond))
title('log10 of Condition number in the moments')
subplot(223)
hist(log10(idelre.cond))
title('log10 of Condition number in the LRE model')
saveas(gcf,[IdentifDirectoryName,'/',M_.fname,'_ident_COND'])
eval(['print -depsc2 ' IdentifDirectoryName '/' M_.fname '_ident_COND']);
eval(['print -dpdf ' IdentifDirectoryName '/' M_.fname '_ident_COND']);
ifig=0;
nbox = min(estim_params_.np-1,12);
for j=1:estim_params_.np,  
    if mod(j,12)==1,    
        ifig = ifig+1;    
        figure('name','Partial correlations in the LRE model'),    
        iplo=0;  
    end
    iplo=iplo+1;    
    mmm = mean(squeeze(idelre.Pco(:,j,:))');    
    [sss, immm] = sort(-mmm);    
    subplot(3,4,iplo), 
    myboxplot(squeeze(idelre.Pco(immm(2:nbox+1),j,:))'),
    set(gca,'ylim',[0 1])
    set(gca,'xticklabel','')
    for ip=1:nbox, %estim_params_.np,  
        text(ip,-0.02,deblank(M_.param_names(estim_params_.param_vals(immm(ip+1),1),:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
    end
    title(deblank(M_.param_names(estim_params_.param_vals(j,1),:))),
    if j==estim_params_.np | mod(j,12)==0  
        saveas(gcf,[IdentifDirectoryName,'/',M_.fname,'_ident_PCORR_LRE',int2str(ifig)])
        eval(['print -depsc2 ' IdentifDirectoryName '/' M_.fname '_ident_PCORR_LRE',int2str(ifig)]);
        eval(['print -dpdf ' IdentifDirectoryName '/' M_.fname '_ident_PCORR_LRE',int2str(ifig)]);
    end
end
