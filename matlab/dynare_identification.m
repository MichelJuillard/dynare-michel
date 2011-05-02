function [pdraws, TAU, GAM, LRE, gp, H, JJ] = dynare_identification(options_ident, pdraws0)
%function [pdraws, TAU, GAM, LRE, gp, H, JJ] = dynare_identification(options_ident, pdraws0)
%
% INPUTS
%    o options_ident    [structure] identification options
%    o pdraws0          [matrix] optional: matrix of MC sample of model params. 
%    
% OUTPUTS
%    o pdraws           [matrix] matrix of MC sample of model params used
%    o TAU,             [matrix] MC sample of entries in the model solution (stacked vertically)
%    o GAM,             [matrix] MC sample of entries in the moments (stacked vertically)
%    o LRE,             [matrix] MC sample of entries in LRE model (stacked vertically)
%    o gp,              [matrix] derivatives of the Jacobian (LRE model)
%    o H,               [matrix] derivatives of the model solution
%    o JJ               [matrix] derivatives of the  moments
%    
% SPECIAL REQUIREMENTS
%    None

% main 
%
% Copyright (C) 2010-2011 Dynare Team
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

if exist('OCTAVE_VERSION')
    warning('off'),
else
    warning off,
end

fname_ = M_.fname;
options_ident = set_default_option(options_ident,'load_ident_files',0);
options_ident = set_default_option(options_ident,'useautocorr',0);
options_ident = set_default_option(options_ident,'ar',3);
options_ident = set_default_option(options_ident,'prior_mc',1);
options_ident = set_default_option(options_ident,'prior_range',0);
options_ident = set_default_option(options_ident,'periods',300);
options_ident = set_default_option(options_ident,'replic',100);
options_ident = set_default_option(options_ident,'advanced',0);
if nargin==2,
    options_ident.prior_mc=size(pdraws0,1);
end
if isempty(estim_params_),
    options_ident.prior_mc=1;
    options_ident.prior_range=0;
    prior_exist=0;
else
    prior_exist=1;
end

% options_ident.load_ident_files=1;
iload = options_ident.load_ident_files;
%options_ident.advanced=1;
advanced = options_ident.advanced;
nlags = options_ident.ar;
periods = options_ident.periods;
replic = options_ident.replic;
useautocorr = options_ident.useautocorr;
options_.ar=nlags;
options_.prior_mc = options_ident.prior_mc;
options_.options_ident = options_ident;
options_.Schur_vec_tol = 1.e-8;

options_ = set_default_option(options_,'datafile',[]);
options_.mode_compute = 0;
options_.plot_priors = 0;
[data,rawdata,xparam1,data_info]=dynare_estimation_init(M_.endo_names,fname_,1);
if isempty(data_info),
    data_info.gend = periods;
    data_info.data = [];
    data_info.data_index = [];
    data_info.number_of_observations = periods*size(options_.varobs,1);
    data_info.no_more_missing_observations = 0;
    data_info.missing_value = 0;
end


SampleSize = options_ident.prior_mc;

% results = prior_sampler(0,M_,bayestopt_,options_,oo_);

if prior_exist
    if options_ident.prior_range
        prior_draw(1,1);
    else
        prior_draw(1);
    end
end

if ~(exist('sylvester3mr','file')==2),

    dynareroot = strrep(which('dynare'),'dynare.m','');
    addpath([dynareroot 'gensylv'])
end

IdentifDirectoryName = CheckPath('identification');
if prior_exist,

    indx = [];
    if ~isempty(estim_params_.param_vals),
        indx = estim_params_.param_vals(:,1);
    end
    indexo=[];
    if ~isempty(estim_params_.var_exo)
        indexo = estim_params_.var_exo(:,1);
    end

    nparam = length(bayestopt_.name);
    np = estim_params_.np;
    name = bayestopt_.name;
    name_tex = char(M_.exo_names_tex(indexo,:),M_.param_names_tex(indx,:));

    offset = estim_params_.nvx;
    offset = offset + estim_params_.nvn;
    offset = offset + estim_params_.ncx;
    offset = offset + estim_params_.ncn;
else
    indx = [1:M_.param_nbr];
    indexo = [1:M_.exo_nbr];
    offset = M_.exo_nbr;
    np = M_.param_nbr;
    nparam = np+offset;
    name = [cellstr(M_.exo_names); cellstr(M_.param_names)];
    name_tex = [cellstr(M_.exo_names_tex); cellstr(M_.param_names_tex)];
end

options_ident = set_default_option(options_ident,'max_dim_cova_group',min([2,nparam-1]));
max_dim_cova_group = options_ident.max_dim_cova_group;


MaxNumberOfBytes=options_.MaxNumberOfBytes;

disp(' ')
disp(['==== Identification analysis ====' ]),
disp(' ')

if iload <=0,
    
    [I,J]=find(M_.lead_lag_incidence');
    if prior_exist,
        if exist([fname_,'_mean.mat'],'file'),
%             disp('Testing posterior mean')
            load([fname_,'_mean'],'xparam1')
            pmean = xparam1';
            clear xparam1
        end
        if exist([fname_,'_mode.mat'],'file'),
%             disp('Testing posterior mode')
            load([fname_,'_mode'],'xparam1')
            pmode = xparam1';
            clear xparam1
        end
        disp('Testing prior mean')
        params = set_prior(estim_params_,M_,options_)';
    else
        params = [sqrt(diag(M_.Sigma_e))', M_.params'];
    end
    [idehess_prior, idemoments_prior, idemodel_prior, idelre_prior, derivatives_info_prior] = ...
        identification_analysis(params,indx,indexo,options_ident,data_info, prior_exist, name_tex,1);
    idehess_prior.params=params;
%     siH = idemodel_prior.siH;
%     siJ = idemoments_prior.siJ;
%     siLRE = idelre_prior.siLRE;
%     normH = max(abs(siH)')';
%     normJ = max(abs(siJ)')';
%     normLRE = max(abs(siLRE)')';
    save([IdentifDirectoryName '/' M_.fname '_identif.mat'], 'idehess_prior', 'idemoments_prior','idemodel_prior', 'idelre_prior')
    disp_identification(params, idemodel_prior, idemoments_prior, name);
    plot_identification(params,idemoments_prior,idehess_prior,idemodel_prior,idelre_prior,advanced,'Prior mean - ',name);


    if exist('pmode','var'),
         disp('Testing posterior mode')
         [idehess_mode, idemoments_mode, idemodel_mode, idelre_mode, derivatives_info_mode] = ...
             identification_analysis(pmode,indx,indexo,options_ident,data_info, prior_exist, name_tex,0);
         idehess_mode.params=pmode;
         save([IdentifDirectoryName '/' M_.fname '_identif.mat'], 'idehess_mode', 'idemoments_mode','idemodel_mode', 'idelre_mode','-append')
         disp_identification(params, idemodel_mode, idemoments_mode, name);
         plot_identification(params,idemoments_mode,idehess_mode,idemodel_mode,idelre_mode,advanced,'Posterior mode - ',name);
    end
    if exist('pmean','var'),
         disp('Testing posterior mean')
         [idehess_mean, idemoments_mean, idemodel_mean, idelre_mean, derivatives_info_mean] = ...
             identification_analysis(pmean,indx,indexo,options_ident,data_info, prior_exist, name_tex,0);
         idehess_mean.params=pmean;
         save([IdentifDirectoryName '/' M_.fname '_identif.mat'], 'idehess_mean', 'idemoments_mean','idemodel_mean', 'idelre_mean','-append')
         disp_identification(params, idemodel_mean, idemoments_mean, name);
         plot_identification(params,idemoments_mean,idehess_mean,idemodel_mean,idelre_mean,advanced,'Posterior mean - ',name);
    end

    if SampleSize > 1,
        disp(' ')
        disp('Monte Carlo Testing')
        h = waitbar(0,'Monte Carlo identification checks ...');
        iteration = 0;
        loop_indx = 0;
        file_index = 0;
        run_index = 0;
        options_MC=options_ident;
        options_MC.advanced=0;
    else
        iteration = 1;
        pdraws = [];
    end
    while iteration < SampleSize,
        loop_indx = loop_indx+1;
        if nargin==2,
            params = pdraws0(iteration+1,:);
        else
            params = prior_draw();
        end
        [~, ideJ, ideH, ideGP, ~ , info] = ...
            identification_analysis(params,indx,indexo,options_MC,data_info, prior_exist, name_tex,0);
        if iteration==0,
            MAX_tau   = min(SampleSize,ceil(MaxNumberOfBytes/(size(ideH.siH,1)*nparam)/8));
            stoH = zeros([size(ideH.siH,1),nparam,MAX_tau]);
            stoJJ = zeros([size(ideJ.siJ,1),nparam,MAX_tau]);
            stoLRE = zeros([size(ideGP.siLRE,1),np,MAX_tau]);
            TAU = zeros(size(ideH.siH,1),SampleSize);
            GAM = zeros(size(ideJ.siJ,1),SampleSize);
            LRE = zeros(size(ideGP.siLRE,1),SampleSize);
            pdraws = zeros(SampleSize,nparam);
            idemoments.indJJ = ideJ.indJJ;
            idemodel.indH = ideH.indH;
            idelre.indLRE = ideGP.indLRE;
            idemoments.ind0 = zeros(SampleSize,nparam);
            idemodel.ind0 = zeros(SampleSize,nparam);
            idelre.ind0 = zeros(SampleSize,np);
            idemoments.jweak = zeros(SampleSize,nparam);
            idemodel.jweak = zeros(SampleSize,nparam);
            idelre.jweak = zeros(SampleSize,np);
            idemoments.jweak_pair = zeros(SampleSize,nparam*(nparam+1)/2);
            idemodel.jweak_pair = zeros(SampleSize,nparam*(nparam+1)/2);
            idelre.jweak_pair = zeros(SampleSize,np*(np+1)/2);
            idemoments.cond = zeros(SampleSize,1);
            idemodel.cond = zeros(SampleSize,1);
            idelre.cond = zeros(SampleSize,1);
            idemoments.Mco = zeros(SampleSize,nparam);
            idemodel.Mco = zeros(SampleSize,nparam);
            idelre.Mco = zeros(SampleSize,np);
            delete([IdentifDirectoryName '/' M_.fname '_identif_*.mat'])
        end
        if info(1)==0,
            iteration = iteration + 1;
            run_index = run_index + 1;
            TAU(:,iteration)=ideH.TAU;
            LRE(:,iteration)=ideGP.LRE;
            GAM(:,iteration)=ideJ.GAM;
            idemoments.cond(iteration,1)=ideJ.cond;
            idemodel.cond(iteration,1)=ideH.cond;
            idelre.cond(iteration,1)=ideGP.cond;
            idemoments.ino(iteration,1)=ideJ.ino;
            idemodel.ino(iteration,1)=ideH.ino;
            idelre.ino(iteration,1)=ideGP.ino;
            idemoments.ind0(iteration,:)=ideJ.ind0;
            idemodel.ind0(iteration,:)=ideH.ind0;
            idelre.ind0(iteration,:)=ideGP.ind0;
            idemoments.jweak(iteration,:)=ideJ.jweak;
            idemodel.jweak(iteration,:)=ideH.jweak;
            idelre.jweak(iteration,:)=ideGP.jweak;
            idemoments.jweak_pair(iteration,:)=ideJ.jweak_pair;
            idemodel.jweak_pair(iteration,:)=ideH.jweak_pair;
            idelre.jweak_pair(iteration,:)=ideGP.jweak_pair;
            idemoments.Mco(iteration,:)=ideJ.Mco;
            idemodel.Mco(iteration,:)=ideH.Mco;
            idelre.Mco(iteration,:)=ideGP.Mco;
            stoLRE(:,:,run_index) = ideGP.siLRE;
            stoH(:,:,run_index) = ideH.siH;
            stoJJ(:,:,run_index) = ideJ.siJ;
            pdraws(iteration,:) = params;
            if run_index==MAX_tau || iteration==SampleSize,
                file_index = file_index + 1;
                if run_index<MAX_tau,
                    stoH = stoH(:,:,1:run_index);
                    stoJJ = stoJJ(:,:,1:run_index);
                    stoLRE = stoLRE(:,:,1:run_index);
                end
                save([IdentifDirectoryName '/' M_.fname '_identif_' int2str(file_index) '.mat'], 'stoH', 'stoJJ', 'stoLRE')
                run_index = 0;
                stoH = zeros(size(stoH));
                stoJJ = zeros(size(stoJJ));
                stoLRE = zeros(size(stoLRE));
                
            end
            
            if SampleSize > 1,
                waitbar(iteration/SampleSize,h,['MC identification checks ',int2str(iteration),'/',int2str(SampleSize)])
            end
        end
        
    end
    
    
    if SampleSize > 1,
        close(h),
        normTAU=std(TAU')';
        normLRE=std(LRE')';
        normGAM=std(GAM')';
        normaliz1=std(pdraws);
        iter=0;
        for ifile_index = 1:file_index,
            load([IdentifDirectoryName '/' M_.fname '_identif_' int2str(ifile_index) '.mat'], 'stoH', 'stoJJ', 'stoLRE')
            for irun=1:size(stoH,3),
                iter=iter+1;
                siJnorm(iter,:) = vnorm(stoJJ(:,:,irun)./repmat(normGAM,1,nparam)).*normaliz1;
                siHnorm(iter,:) = vnorm(stoH(:,:,irun)./repmat(normTAU,1,nparam)).*normaliz1;
                siLREnorm(iter,:) = vnorm(stoLRE(:,:,irun)./repmat(normLRE,1,nparam-offset)).*normaliz1(offset+1:end);
            end
            
        end
        idemoments.siJnorm = siJnorm;
        idemodel.siHnorm = siHnorm;
        idelre.siLREnorm = siLREnorm;
        save([IdentifDirectoryName '/' M_.fname '_identif.mat'], 'pdraws', 'idemodel', 'idemoments', 'idelre', ... %'indJJ', 'indH', 'indLRE', ...
            'TAU', 'GAM', 'LRE','-append')
    else
        siJnorm = idemoments_prior.siJnorm;
        siHnorm = idemodel_prior.siHnorm;
        siLREnorm = idelre_prior.siLREnorm;
    end
    
else
    load([IdentifDirectoryName '/' M_.fname '_identif'])
    identFiles = dir([IdentifDirectoryName '/' M_.fname '_identif_*']);
    options_ident.prior_mc=size(pdraws,1);
    SampleSize = options_ident.prior_mc;
    options_.options_ident = options_ident;
end  

if nargout>3 && iload,
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

if iload,
    disp('Testing prior mean')
    disp_identification(idehess_prior.params, idemodel_prior, idemoments_prior, name);
    if advanced, disp('Press ENTER to display advanced diagnostics'), pause, end
    plot_identification(idehess_prior.params,idemoments_prior,idehess_prior,idemodel_prior,idelre_prior,advanced,'Prior mean - ',name);
end
if SampleSize > 1,
    fprintf('\n\n')
    disp('Testing MC sample')
    disp_identification(pdraws, idemodel, idemoments, name);
    plot_identification(pdraws,idemoments,idehess_prior,idemodel,idelre,advanced,'MC sample - ',name, IdentifDirectoryName, 1);
    if advanced,
        jcrit=find(idemoments.ino);
        for j=1:length(jcrit),
            if ~iload,
                [idehess_(j), idemoments_(j), idemodel_(j), idelre_(j), derivatives_info_(j)] = ...
                    identification_analysis(pdraws(jcrit(j),:),indx,indexo,options_ident,data_info, prior_exist, name_tex,1);
            end
            tittxt = ['Critical draw n. ',int2str(j),' - '];
            fprintf('\n\n')
            disp(['Testing ',tittxt, 'Press ENTER']), pause,
            plot_identification(pdraws(jcrit(j),:),idemoments_(j),idehess_(j),idemodel_(j),idelre_(j),1,tittxt,name);
%             disp_identification(pdraws(jcrit(j),:), idemodel_(j), idemoments_(j), name, advanced);
        end        
        if ~iload,
            save([IdentifDirectoryName '/' M_.fname '_identif.mat'], 'idehess_', 'idemoments_','idemodel_', 'idelre_', 'jcrit', '-append');
        end
    end
end

if exist('OCTAVE_VERSION')
    warning('on'),
else
    warning on,
end

disp(' ')
disp(['==== Identification analysis completed ====' ]),
disp(' ')
disp(' ')
