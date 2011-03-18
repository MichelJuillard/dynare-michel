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

iload = options_ident.load_ident_files;
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
[data,rawdata,xparam1,data_info]=dynare_estimation_init([],fname_,1);
if isempty(data_info),
    data_info.gend = periods;
    data_info.data = [];
    data_info.data_index = [];
    data_info.number_of_observations = periods*length(options_.varobs);
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

MaxNumberOfBytes=options_.MaxNumberOfBytes;

disp(' ')
disp(['==== Identification analysis ====' ]),
disp(' ')

if iload <=0,
    
    iteration = 0;
    burnin_iteration = 0;
    if SampleSize==1,
        BurninSampleSize=0;
    else
        BurninSampleSize=50;
    end
    loop_indx = 0;
    file_index = 0;
    run_index = 0;

    if SampleSize > 1,
        h = waitbar(0,'Monte Carlo identification checks ...');
    end
    [I,J]=find(M_.lead_lag_incidence');

    while iteration < SampleSize,
        loop_indx = loop_indx+1;
        if prior_exist,
            if SampleSize==1,
                if exist([fname_,'_mean.mat'],'file'),
                    disp('Testing posterior mean')
                    load([fname_,'_mean'],'xparam1')
                    params = xparam1';
                    clear xparam1
                elseif exist([fname_,'_mode.mat'],'file'),
                    disp('Testing posterior mode')
                    load([fname_,'_mode'],'xparam1')
                    params = xparam1';
                    clear xparam1
                else
                    disp('Testing prior mean')
                    params = set_prior(estim_params_,M_,options_)';
                end
            else
                if nargin==2,
                    if burnin_iteration>=BurninSampleSize,
                        params = pdraws0(iteration+1,:);
                    else
                        params = pdraws0(burnin_iteration+1,:);
                    end
                else
                    params = prior_draw();
                end
            end
            set_all_parameters(params);
        else
            params = [sqrt(diag(M_.Sigma_e))', M_.params'];
        end
        [A,B,ys,info]=dynare_resolve;

        
        if info(1)==0,
            oo0=oo_;
            tau=[oo_.dr.ys(oo_.dr.order_var); vec(A); dyn_vech(B*M_.Sigma_e*B')];
            yy0=oo_.dr.ys(I);    
            [residual, g1 ] = feval([M_.fname,'_dynamic'],yy0, ...
                                    oo_.exo_steady_state', M_.params, ...
                                    oo_.dr.ys, 1);    

            if burnin_iteration<BurninSampleSize,
                burnin_iteration = burnin_iteration + 1;
                pdraws(burnin_iteration,:) = params;
                TAU(:,burnin_iteration)=tau;
                LRE(:,burnin_iteration)=[oo_.dr.ys(oo_.dr.order_var); vec(g1)];
                [gam,stationary_vars] = th_autocovariances(oo0.dr,bayestopt_.mfys,M_,options_);
                if exist('OCTAVE_VERSION')
                    warning('off')
                else
                    warning off,
                end
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
                dum = dyn_vech(gam{1});
                for j=1:nlags,
                    dum = [dum; vec(gam{j+1})];
                end
                GAM(:,burnin_iteration)=[oo_.dr.ys(bayestopt_.mfys); dum];
                %                 warning warning_old_state;
            else
                iteration = iteration + 1;
                run_index = run_index + 1;
                if iteration==1 && BurninSampleSize,
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
                [JJ, H, gam, gp, dA, dOm, dYss] = getJJ(A, B, M_,oo0,options_,0,indx,indexo,bayestopt_.mf2,nlags,useautocorr);
                derivatives_info.dA=dA;
                derivatives_info.dOm=dOm;
                derivatives_info.dYss=dYss;
                if BurninSampleSize == 0,
                    indJJ = (find(max(abs(JJ'))>1.e-8));
                    indH = (find(max(abs(H'))>1.e-8));
                    indLRE = (find(max(abs(gp'))>1.e-8));
                    TAU = zeros(length(indH),SampleSize);
                    GAM = zeros(length(indJJ),SampleSize);
                    LRE = zeros(length(indLRE),SampleSize);
                    MAX_tau   = min(SampleSize,ceil(MaxNumberOfBytes/(length(indH)*nparam)/8));
                    MAX_gam   = min(SampleSize,ceil(MaxNumberOfBytes/(length(indJJ)*nparam)/8));
                    stoH = zeros([length(indH),nparam,MAX_tau]);
                    stoJJ = zeros([length(indJJ),nparam,MAX_tau]);
                    delete([IdentifDirectoryName '/' M_.fname '_identif_*.mat'])
                end
                TAU(:,iteration)=tau(indH);
                vg1 = [oo_.dr.ys(oo_.dr.order_var); vec(g1)];
                LRE(:,iteration)=vg1(indLRE);
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
                
                siJnorm(iteration,:) = vnorm(siJ./repmat(GAM(:,iteration),1,nparam)).*params;        
                siHnorm(iteration,:) = vnorm(siH./repmat(TAU(:,iteration),1,nparam)).*params;        
                siLREnorm(iteration,:) = vnorm(siLRE./repmat(LRE(:,iteration),1,nparam-offset)).*params(offset+1:end);        
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
                normH = max(abs(stoH(:,:,run_index))')';
                normJ = max(abs(stoJJ(:,:,run_index))')';
                normLRE = max(abs(stoLRE(:,:,run_index))')';
                %               normH = TAU(:,iteration);
                %               normJ = GAM(:,iteration);
                %               normLRE = LRE(:,iteration);
                [idemodel.Mco(:,iteration), idemoments.Mco(:,iteration), idelre.Mco(:,iteration), ...
                    idemodel.Pco(:,:,iteration), idemoments.Pco(:,:,iteration), idelre.Pco(:,:,iteration), ...
                    idemodel.cond(iteration), idemoments.cond(iteration), idelre.cond(iteration), ...
                    idemodel.ee(:,:,iteration), idemoments.ee(:,:,iteration), idelre.ee(:,:,iteration), ...
                    idemodel.ind(:,iteration), idemoments.ind(:,iteration), ...
                    idemodel.indno{iteration}, idemoments.indno{iteration}, ...
                    idemodel.ino(iteration), idemoments.ino(iteration)] = ...
                    identification_checks(H(indH,:)./normH(:,ones(nparam,1)),JJ(indJJ,:)./normJ(:,ones(nparam,1)), gp(indLRE,:)./normLRE(:,ones(size(gp,2),1)));
                %                  identification_checks(H(indH,:),JJ(indJJ,:), gp(indLRE,:), bayestopt_);
                indok = find(max(idemoments.indno{iteration},[],1)==0);
                if iteration ==1 && ~isempty(indok),
                    ide_strength_J=NaN(1,nparam);
                    ide_strength_J_prior=NaN(1,nparam);
                    if advanced,
                        [pars, cosnJ] = ident_bruteforce(JJ(indJJ,:)./normJ(:,ones(nparam,1)),2,1,name_tex);
                    end
                    normaliz = abs(params);
                    if prior_exist,
                        if ~isempty(estim_params_.var_exo),
                            normaliz1 = estim_params_.var_exo(:,7); % normalize with prior standard deviation
                        else
                            normaliz1=[];
                        end
                        if ~isempty(estim_params_.param_vals),
                            normaliz1 = [normaliz1; estim_params_.param_vals(:,7)]'; % normalize with prior standard deviation
                        end
%                         normaliz = max([normaliz; normaliz1]);
                    else
                        normaliz1 = ones(1,nparam);
                    end
                    replic = max([replic, length(indJJ)*3]);
                    try,
                        options_.irf = 0;
                        options_.noprint = 1;
                        options_.order = 1;
                        options_.periods = data_info.gend+100;
                        info = stoch_simul(options_.varobs);
                        datax=oo_.endo_simul(options_.varobs_id,100+1:end);
%                         datax=data;
                        [fval,cost_flag,ys,trend_coeff,info,DLIK,AHess] = DsgeLikelihood(params',data_info.gend,datax,data_info.data_index,data_info.number_of_observations,derivatives_info);
                        cparam = inv(-AHess);
                        normaliz(indok) = sqrt(diag(cparam))';
                        cmm = siJ*((-AHess)\siJ');
                    catch,
                        cmm = simulated_moment_uncertainty(indJJ, periods, replic);
                        %                 Jinv=(siJ(:,indok)'*siJ(:,indok))\siJ(:,indok)';
                        %                 MIM=inv(Jinv*cmm*Jinv');
                        MIM=siJ(:,indok)'*(cmm\siJ(:,indok));
                        deltaM = sqrt(diag(MIM));
                        tildaM = MIM./((deltaM)*(deltaM'));
                        rhoM=sqrt(1-1./diag(inv(tildaM)));
                        deltaM = deltaM.*normaliz(indok)';
                        normaliz(indok) = sqrt(diag(inv(MIM)))';
                    end
                    ide_strength_J(indok) = (1./(normaliz(indok)'./abs(params(indok)')));
                    ide_strength_J_prior(indok) = (1./(normaliz(indok)'./normaliz1(indok)'));
                    %                 indok = find(max(idemodel.indno{iteration},[],1)==0);
                    %                 ide_strength_H(iteration,:)=zeros(1,nparam);
                    %                 mim=inv(siH(:,indok)'*siH(:,indok))*siH(:,indok)';
                    % %                 mim=mim*diag(GAM(:,iteration))*mim';
                    % %                 MIM=inv(mim);
                    %                 mim=mim.*repmat(TAU(:,iteration),1,length(indok))';
                    %                 MIM=inv(mim*mim');
                    %                 deltaM = sqrt(diag(MIM));
                    %                 tildaM = MIM./((deltaM)*(deltaM'));
                    %                 rhoM=sqrt(1-1./diag(inv(tildaM)));
                    %                 deltaM = deltaM.*params(indok)';
                    %                 ide_strength_H(iteration,indok) = (1./[sqrt(diag(inv(MIM)))./params(indok)']);
                    %                 inok = find((abs(GAM(:,iteration))==0));
                    %                 isok = find((abs(GAM(:,iteration))));
                    %                 quant(isok,:) = siJ(isok,:)./repmat(GAM(isok,iteration),1,nparam);
                    %                 quant(inok,:) = siJ(inok,:)./repmat(mean(abs(GAM(:,iteration))),length(inok),nparam);
                    quant = siJ./repmat(sqrt(diag(cmm)),1,nparam);
                    siJnorm(iteration,:) = vnorm(quant).*normaliz;
                    %                 siJnorm(iteration,:) = vnorm(siJ(inok,:)).*normaliz;
                    quant=[];
                    inok = find((abs(TAU(:,iteration))==0));
                    isok = find((abs(TAU(:,iteration))));
                    quant(isok,:) = siH(isok,:)./repmat(TAU(isok,iteration),1,nparam);
                    quant(inok,:) = siH(inok,:)./repmat(mean(abs(TAU(:,iteration))),length(inok),nparam);
                    siHnorm(iteration,:) = vnorm(quant).*normaliz;
                    %                 siHnorm(iteration,:) = vnorm(siH./repmat(TAU(:,iteration),1,nparam)).*normaliz;
                    quant=[];
                    inok = find((abs(LRE(:,iteration))==0));
                    isok = find((abs(LRE(:,iteration))));
                    quant(isok,:) = siLRE(isok,:)./repmat(LRE(isok,iteration),1,np);
                    quant(inok,:) = siLRE(inok,:)./repmat(mean(abs(LRE(:,iteration))),length(inok),np);
                    siLREnorm(iteration,:) = vnorm(quant).*normaliz(offset+1:end);
                    %                 siLREnorm(iteration,:) = vnorm(siLRE./repmat(LRE(:,iteration),1,nparam-offset)).*normaliz(offset+1:end);
                end,
                if run_index==MAX_tau || iteration==SampleSize,
                    file_index = file_index + 1;
                    if run_index<MAX_tau,
                        stoH = stoH(:,:,1:run_index);
                        stoJJ = stoJJ(:,:,1:run_index);
                        stoLRE = stoLRE(:,:,1:run_index);
                    end
                    save([IdentifDirectoryName '/' M_.fname '_identif_' int2str(file_index) '.mat'], 'stoH', 'stoJJ', 'stoLRE')
                    run_index = 0;
                    
                end
                
                if SampleSize > 1,
                    waitbar(iteration/SampleSize,h,['MC Identification checks ',int2str(iteration),'/',int2str(SampleSize)])
                end
            end
        end
    end
    
    
    if SampleSize > 1,
        close(h)
    end
    
    
    save([IdentifDirectoryName '/' M_.fname '_identif.mat'], 'pdraws', 'idemodel', 'idemoments', 'idelre', 'indJJ', 'indH', 'indLRE', ...
        'siHmean', 'siJmean', 'siLREmean', 'derHmean', 'derJmean', 'derLREmean', 'TAU', 'GAM', 'LRE')
else
    load([IdentifDirectoryName '/' M_.fname '_identif'], 'pdraws', 'idemodel', 'idemoments', 'idelre', 'indJJ', 'indH', 'indLRE', ...
        'siHmean', 'siJmean', 'siLREmean', 'derHmean', 'derJmean', 'derLREmean', 'TAU', 'GAM', 'LRE')
    identFiles = dir([IdentifDirectoryName '/' M_.fname '_identif_*']);
    options_ident.prior_mc=size(pdraws,1);
    SampleSize = options_ident.prior_mc;
    options_.options_ident = options_ident;
    if iload>1,
        idemodel0=idemodel;
        idemoments0=idemoments;
        idelre0 = idelre;
        iteration = 0;
        h = waitbar(0,'Monte Carlo identification checks ...');
        for file_index=1:length(identFiles)
            load([IdentifDirectoryName '/' M_.fname '_identif_' int2str(file_index)], 'stoH', 'stoJJ', 'stoLRE')
            for index=1:size(stoH,3),
                iteration = iteration+1;
                normH = max(abs(stoH(:,:,index))')';
                normJ = max(abs(stoJJ(:,:,index))')';
                normLRE = max(abs(stoLRE(:,:,index))')';
                %             normH = TAU(:,iteration);
                %             normJ = GAM(:,iteration);
                %             normLRE = LRE(:,iteration);
                [idemodel.Mco(:,iteration), idemoments.Mco(:,iteration), idelre.Mco(:,iteration), ...
                    idemodel.Pco(:,:,iteration), idemoments.Pco(:,:,iteration), idelre.Pco(:,:,iteration), ...
                    idemodel.cond(iteration), idemoments.cond(iteration), idelre.cond(iteration), ...
                    idemodel.ee(:,:,iteration), idemoments.ee(:,:,iteration), idelre.ee(:,:,iteration), ...
                    idemodel.ind(:,iteration), idemoments.ind(:,iteration), ...
                    idemodel.indno{iteration}, idemoments.indno{iteration}, ...
                    idemodel.ino(iteration), idemoments.ino(iteration)] = ...
                    identification_checks(stoH(:,:,index)./normH(:,ones(nparam,1)), ...
                    stoJJ(:,:,index)./normJ(:,ones(nparam,1)), ...
                    stoLRE(:,:,index)./normLRE(:,ones(size(stoLRE,2),1)));
                waitbar(iteration/SampleSize,h,['MC Identification checks ',int2str(iteration),'/',int2str(SampleSize)])
            end
        end
        close(h);
        save([IdentifDirectoryName '/' M_.fname '_identif.mat'], 'idemodel', 'idemoments', 'idelre', '-append')
    end
    iteration = 0;
    h = waitbar(0,'Monte Carlo identification checks ...');
    for file_index=1:length(identFiles)
        load([IdentifDirectoryName '/' M_.fname '_identif_' int2str(file_index)], 'stoH', 'stoJJ', 'stoLRE')
        for index=1:size(stoH,3),
            iteration = iteration+1;
            fobj(iteration)=sum(((GAM(:,iteration)-GAM(:,1))).^2);
            fobjH(iteration)=sum(((TAU(:,iteration)-TAU(:,1))).^2);
            fobjR(iteration)=sum(((GAM(:,iteration)./GAM(:,1)-1)).^2);
            fobjHR(iteration)=sum(((TAU(:,iteration)./TAU(:,1)-1)).^2);
            FOC(iteration,:) = (GAM(:,iteration)-GAM(:,1))'*stoJJ(:,:,index);
            FOCH(iteration,:) = (TAU(:,iteration)-TAU(:,1))'*stoH(:,:,index);
            FOCR(iteration,:) = ((GAM(:,iteration)./GAM(:,1)-1)./GAM(:,1))'*stoJJ(:,:,index);
            FOCHR(iteration,:) = ((TAU(:,iteration)./TAU(:,1)-1)./TAU(:,1))'*stoH(:,:,index);
            
            waitbar(iteration/SampleSize,h,['MC Identification checks ',int2str(iteration),'/',int2str(SampleSize)])
        end
    end
    close(h);
end  

if SampleSize>1,
    siJmean = siJmean.*(ones(length(indJJ),1)*std(pdraws));
    siHmean = siHmean.*(ones(length(indH),1)*std(pdraws));
    siLREmean = siLREmean.*(ones(length(indLRE),1)*std(pdraws(:, offset+1:end )));

    derJmean = derJmean.*(ones(length(indJJ),1)*std(pdraws));
    derHmean = derHmean.*(ones(length(indH),1)*std(pdraws));
    derLREmean = derLREmean.*(ones(length(indLRE),1)*std(pdraws(:, offset+1:end )));

    derHmean = abs(derHmean./(max(siHmean')'*ones(1,size(pdraws,2))));
    derJmean = abs(derJmean./(max(siJmean')'*ones(1,size(pdraws,2))));
    derLREmean = abs(derLREmean./(max(siLREmean')'*ones(1,np)));

    siHmean = siHmean./(max(siHmean')'*ones(1,size(pdraws,2)));
    siJmean = siJmean./(max(siJmean')'*ones(1,size(pdraws,2)));
    siLREmean = siLREmean./(max(siLREmean')'*ones(1,np));

    tstJmean = derJmean*0;
    tstHmean = derHmean*0;
    tstLREmean = derLREmean*0;

    if exist('OCTAVE_VERSION')
        warning('off'),
    else
        warning off,
    end

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

disp_identification(pdraws, idemodel, idemoments, name, advanced)


if advanced,
    figure('Name','Identification LRE model form'),
    subplot(211)
    if SampleSize > 1,
        mmm = mean(siLREnorm);
    else
        mmm = (siLREnorm);
    end
    [ss, is] = sort(mmm);
    if SampleSize ==1,
        bar(siLREnorm(:,is))
    else
        myboxplot(log(siLREnorm(:,is)))
    end
    % mmm = mean(siLREmean);
    % [ss, is] = sort(mmm);
    % myboxplot(siLREmean(:,is))
    % set(gca,'ylim',[0 1.05])
    set(gca,'xticklabel','')
    for ip=1:np,
        text(ip,-0.02,deblank(M_.param_names(indx(is(ip)),:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
    end
    title('Sensitivity in the LRE model')
    
    subplot(212)
    if SampleSize>1,
        mmm = mean(-idelre.Mco');
    else
        mmm = (-idelre.Mco');
    end
    [ss, is] = sort(mmm);
    if SampleSize ==1,
        bar(idelre.Mco(is,:)')
    else
        myboxplot(idelre.Mco(is,:)')
    end
    set(gca,'ylim',[0 1])
    set(gca,'xticklabel','')
    for ip=1:np,
        text(ip,-0.02,deblank(M_.param_names(indx(is(ip)),:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
    end
    title('Multicollinearity in the LRE model')
    saveas(gcf,[IdentifDirectoryName,'/',M_.fname,'_ident_LRE'])
    eval(['print -depsc2 ' IdentifDirectoryName '/' M_.fname '_ident_LRE']);
    eval(['print -dpdf ' IdentifDirectoryName '/' M_.fname '_ident_LRE']);
    if options_.nograph, close(gcf); end
    
    figure('Name','Identification in the model'),
    % subplot(311)
    %
    % if SampleSize>1,
    % mmm = mean(ide_strength_H);
    % else
    % mmm = (ide_strength_H);
    % end
    % [ss, is] = sort(mmm);
    % if SampleSize>1,
    % myboxplot(ide_strength_H(:,is))
    % else
    % bar(ide_strength_H(:,is))
    % end
    % % set(gca,'ylim',[0 1.05])
    % set(gca,'xticklabel','')
    % dy = get(gca,'ylim');
    % % dy=dy(2)-dy(1);
    % for ip=1:nparam,
    %     text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
    % end
    % title('Identification strength in the model')
    
    subplot(211)
    % mmm = mean(siHmean);
    % [ss, is] = sort(mmm);
    % myboxplot(siHmean(:,is))
    if SampleSize>1,
        mmm = mean(siHnorm);
    else
        mmm = (siHnorm);
    end
    [ss, is] = sort(mmm);
    if SampleSize>1,
        myboxplot(log(siHnorm(:,is)))
    else
        bar(siHnorm(:,is))
    end
    % set(gca,'ylim',[0 1.05])
    set(gca,'xticklabel','')
    dy = get(gca,'ylim');
    % dy=dy(2)-dy(1);
    for ip=1:nparam,
        text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
    end
    title('Sensitivity in the model')
    
    subplot(212)
    if SampleSize>1,
        mmm = mean(-idemodel.Mco');
    else
        mmm = (-idemodel.Mco');
    end
    % [ss, is] = sort(mmm);
    if SampleSize>1,
        myboxplot(idemodel.Mco(is,:)')
    else
        bar(idemodel.Mco(is,:)')
    end
    set(gca,'ylim',[0 1])
    set(gca,'xticklabel','')
    for ip=1:nparam,
        text(ip,-0.02,name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
    end
    title('Multicollinearity in the model')
    saveas(gcf,[IdentifDirectoryName,'/',M_.fname,'_ident_model'])
    eval(['print -depsc2 ' IdentifDirectoryName '/' M_.fname '_ident_model']);
    eval(['print -dpdf ' IdentifDirectoryName '/' M_.fname '_ident_model']);
    if options_.nograph, close(gcf); end
end

figure('Name','Identification in the moments'),
subplot(211)
% mmm = mean(siHmean);
% [ss, is] = sort(mmm);
% myboxplot(siHmean(:,is))
if SampleSize>1,
    mmm = mean(ide_strength_J);
else
    mmm = (ide_strength_J);
end
[ss, is] = sort(mmm);
if SampleSize>1,
    myboxplot(log(ide_strength_J(:,is)))
else
    bar(log([ide_strength_J(:,is)' ide_strength_J_prior(:,is)']))
end
% set(gca,'ylim',[0 1.05])
set(gca,'xlim',[0 nparam+1])
set(gca,'xticklabel','')
dy = get(gca,'ylim');
% dy=dy(2)-dy(1);
for ip=1:nparam,
    text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
legend('relative to param value','relative to prior std','Location','Best')
title('Identification strength in the moments (log-scale)')

subplot(212)
% mmm = mean(siJmean);
% [ss, is] = sort(mmm);
% myboxplot(siJmean(:,is))
if SampleSize > 1,
    mmm = mean(siJnorm);
else
    mmm = (siJnorm);
end
% [ss, is] = sort(mmm);
if SampleSize > 1,
    myboxplot(log(siJnorm(:,is)))
else
    bar(siJnorm(:,is))
end
% set(gca,'ylim',[0 1.05])
set(gca,'xlim',[0 nparam+1])
set(gca,'xticklabel','')
dy = get(gca,'ylim');
% dy=dy(2)-dy(1);
for ip=1:nparam,
    text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Sensitivity in the moments')

% subplot(313)
% if SampleSize>1,
% mmm = mean(-idemoments.Mco');
% else
% mmm = (-idemoments.Mco');
% end
% % [ss, is] = sort(mmm);
% if SampleSize>1,
% myboxplot(idemoments.Mco(is,:)')
% else
% bar(idemoments.Mco(is,:)')
% end
% set(gca,'ylim',[0 1])
% set(gca,'xticklabel','')
% for ip=1:nparam,
%     text(ip,-0.02,name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
% end
% title('Multicollinearity in the moments')
% saveas(gcf,[IdentifDirectoryName,'/',M_.fname,'_ident_moments'])
% eval(['print -depsc2 ' IdentifDirectoryName '/' M_.fname '_ident_moments']);
% eval(['print -dpdf ' IdentifDirectoryName '/' M_.fname '_ident_moments']);
% if options_.nograph, close(gcf); end

if SampleSize==1 && advanced,
    % identificaton patterns
    for  j=1:size(cosnJ,2),
        fprintf('\n\n')
        disp(['Collinearity patterns with ', int2str(j) ,' parameter(s)'])
        fprintf('%-15s [%-*s] %10s\n','Parameter',(15+1)*j,' Expl. params ','cosn')
        for i=1:nparam,
            namx='';
            for in=1:j,
                namx=[namx ' ' sprintf('%-15s',name{pars{i,j}(in)})];
            end
            fprintf('%-15s [%s] %10.3f\n',name{i},namx,cosnJ(i,j))
        end
    end
    disp('')
    [U,S,V]=svd(siJ./normJ(:,ones(nparam,1)),0);
    if nparam<5,
      f1 = figure('name','Identification patterns (moments)');
    else
      f1 = figure('name','Identification patterns (moments): SMALLEST SV');
      f2 = figure('name','Identification patterns (moments): HIGHEST SV');
    end
    for j=1:min(nparam,8),
        if j<5,
            figure(f1),
            jj=j;
        else
            figure(f2),
            jj=j-4;
        end
        subplot(4,1,jj),
        if j<5
            bar(abs(V(:,end-j+1))),
            Stit = S(end-j+1,end-j+1);
%             if j==4 || j==nparam, 
%                 xlabel('SMALLEST singular values'), 
%             end,
        else
            bar(abs(V(:,j))),
            Stit = S(j,j);
%             if j==8 || j==nparam, 
%                 xlabel('LARGEST singular values'), 
%             end,
        end
        set(gca,'xticklabel','')
        if j==4 || j==nparam || j==8, 
        for ip=1:nparam,
            text(ip,-0.02,name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
        end
        end
        title(['Singular value ',num2str(Stit)])
    end
    figure(f1);
    saveas(f1,[IdentifDirectoryName,'/',M_.fname,'_ident_pattern_1'])
    eval(['print -depsc2 ' IdentifDirectoryName '/' M_.fname '_ident_pattern_1']);
    eval(['print -dpdf ' IdentifDirectoryName '/' M_.fname '_ident_pattern_1']);
    if nparam>4,
    figure(f2),
    saveas(f2,[IdentifDirectoryName,'/',M_.fname,'_ident_pattern_2'])
    eval(['print -depsc2 ' IdentifDirectoryName '/' M_.fname '_ident_pattern_2']);
    eval(['print -dpdf ' IdentifDirectoryName '/' M_.fname '_ident_pattern_2']);
    end
end

if advanced,
    if SampleSize>1
        figure('Name','Condition Number'),
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
        if options_.nograph, close(gcf); end
    end
    
    
    pco = NaN(np,np);
    for j=1:np,
        if SampleSize>1
            pco(j+1:end,j) = mean(abs(squeeze(idelre.Pco(j+1:end,j,:))'));
        else
            pco(j+1:end,j) = abs(idelre.Pco(j+1:end,j))';
        end
    end
    figure('name','Pairwise correlations in the LRE model'),
    imagesc(pco',[0 1]);
    set(gca,'xticklabel','')
    set(gca,'yticklabel','')
    for ip=1:nparam,
        text(ip,(0.5),name{ip},'rotation',90,'HorizontalAlignment','left','interpreter','none')
        text(0.5,ip,name{ip},'rotation',0,'HorizontalAlignment','right','interpreter','none')
    end
    colorbar;
    saveas(gcf,[IdentifDirectoryName,'/',M_.fname,'_ident_PCORR_LRE'])
    eval(['print -depsc2 ' IdentifDirectoryName '/' M_.fname '_ident_PCORR_LRE']);
    eval(['print -dpdf ' IdentifDirectoryName '/' M_.fname '_ident_PCORR_LRE']);
    if options_.nograph, close(gcf); end
    
    pco = NaN(nparam,nparam);
    for j=1:nparam,
        if SampleSize>1
            pco(j+1:end,j) = mean(abs(squeeze(idemodel.Pco(j+1:end,j,:))'));
        else
            pco(j+1:end,j) = abs(idemodel.Pco(j+1:end,j))';
        end
    end
    figure('name','Pairwise correlations in the model'),
    imagesc(pco',[0 1]);
    set(gca,'xticklabel','')
    set(gca,'yticklabel','')
    for ip=1:nparam,
        text(ip,(0.5),name{ip},'rotation',90,'HorizontalAlignment','left','interpreter','none')
        text(0.5,ip,name{ip},'rotation',0,'HorizontalAlignment','right','interpreter','none')
    end
    colorbar;
    saveas(gcf,[IdentifDirectoryName,'/',M_.fname,'_ident_PCORR_model'])
    eval(['print -depsc2 ' IdentifDirectoryName '/' M_.fname '_ident_PCORR_model']);
    eval(['print -dpdf ' IdentifDirectoryName '/' M_.fname '_ident_PCORR_model']);
    if options_.nograph, close(gcf); end
    
    for j=1:nparam,
        if SampleSize>1
            pco(j+1:end,j) = mean(abs(squeeze(idemoments.Pco(j+1:end,j,:))'));
        else
            pco(j+1:end,j) = abs(idemoments.Pco(j+1:end,j))';
        end
    end
    figure('name','Pairwise correlations in the moments'),
    imagesc(pco',[0 1]);
    set(gca,'xticklabel','')
    set(gca,'yticklabel','')
    for ip=1:nparam,
        text(ip,(0.5),name{ip},'rotation',90,'HorizontalAlignment','left','interpreter','none')
        text(0.5,ip,name{ip},'rotation',0,'HorizontalAlignment','right','interpreter','none')
    end
    colorbar;
    saveas(gcf,[IdentifDirectoryName,'/',M_.fname,'_ident_PCORR_moments'])
    eval(['print -depsc2 ' IdentifDirectoryName '/' M_.fname '_ident_PCORR_moments']);
    eval(['print -dpdf ' IdentifDirectoryName '/' M_.fname '_ident_PCORR_moments']);
    if options_.nograph, close(gcf); end
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
