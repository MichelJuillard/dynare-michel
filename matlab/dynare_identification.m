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

options_ident = set_default_option(options_ident,'max_ord_bruteforce',min([2,nparam-1]));
max_ord_bruteforce = options_ident.max_ord_bruteforce;


MaxNumberOfBytes=options_.MaxNumberOfBytes;

disp(' ')
disp(['==== Identification analysis ====' ]),
disp(' ')

if iload <=0,
    
    iteration = 0;
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
                    params = pdraws0(iteration+1,:);
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

            iteration = iteration + 1;
            run_index = run_index + 1;

            if iteration,
                [JJ, H, gam, gp, dA, dOm, dYss] = getJJ(A, B, M_,oo0,options_,0,indx,indexo,bayestopt_.mf2,nlags,useautocorr);
                derivatives_info.DT=dA;
                derivatives_info.DOm=dOm;
                derivatives_info.DYss=dYss;
                if iteration==1,
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
                % use prior uncertainty
                siJ = (JJ(indJJ,:));
                siH = (H(indH,:));

                siLRE = (gp(indLRE,:));
                pdraws(iteration,:) = params;
                normH = max(abs(stoH(:,:,run_index))')';
                normJ = max(abs(stoJJ(:,:,run_index))')';
                normLRE = max(abs(stoLRE(:,:,run_index))')';
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
                if iteration ==1,
                    ide_strength_J=NaN(1,nparam);
                    ide_strength_J_prior=NaN(1,nparam);
                end
                if iteration ==1 && advanced,
                    [pars, cosnJ] = ident_bruteforce(JJ(indJJ,:)./normJ(:,ones(nparam,1)),max_ord_bruteforce,options_.TeX,name_tex);
                end
                if iteration ==1 && ~isempty(indok),
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
                        normaliz1(isinf(normaliz1)) = 1;

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
                        derivatives_info.no_DLIK=1;
                        [fval,cost_flag,ys,trend_coeff,info,DLIK,AHess] = DsgeLikelihood(params',data_info.gend,datax,data_info.data_index,data_info.number_of_observations,data_info.no_more_missing_observations,derivatives_info);
                        cparam = inv(-AHess);
                        normaliz(indok) = sqrt(diag(cparam))';
                        cmm = siJ*((-AHess)\siJ');
                        flag_score=1;
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
                        flag_score=0;
                    end
                    ide_strength_J(indok) = (1./(normaliz(indok)'./abs(params(indok)')));
                    ide_strength_J_prior(indok) = (1./(normaliz(indok)'./normaliz1(indok)'));
                    ide_strength_J(params==0)=ide_strength_J_prior(params==0);
                    quant = siJ./repmat(sqrt(diag(cmm)),1,nparam);
                    siJnorm(iteration,:) = vnorm(quant).*normaliz;
                    %                 siJnorm(iteration,:) = vnorm(siJ(inok,:)).*normaliz;
                    quant=[];
                    inok = find((abs(TAU(:,iteration))<1.e-8));
                    isok = find((abs(TAU(:,iteration))));
                    quant(isok,:) = siH(isok,:)./repmat(TAU(isok,iteration),1,nparam);
                    quant(inok,:) = siH(inok,:)./repmat(mean(abs(TAU(:,iteration))),length(inok),nparam);
                    siHnorm(iteration,:) = vnorm(quant).*normaliz;
                    %                 siHnorm(iteration,:) = vnorm(siH./repmat(TAU(:,iteration),1,nparam)).*normaliz;
                    quant=[];
                    inok = find((abs(LRE(:,iteration))<1.e-8));
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
    end
    
    save([IdentifDirectoryName '/' M_.fname '_identif.mat'], 'pdraws', 'idemodel', 'idemoments', 'idelre', 'indJJ', 'indH', 'indLRE', ...
        'siJnorm', 'siHnorm', 'siLREnorm', 'TAU', 'GAM', 'LRE')
else
    load([IdentifDirectoryName '/' M_.fname '_identif'], 'pdraws', 'idemodel', 'idemoments', 'idelre', 'indJJ', 'indH', 'indLRE', ...
        'siJnorm', 'siHnorm', 'siLREnorm', 'TAU', 'GAM', 'LRE')
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

disp_identification(pdraws, idemodel, idemoments, name, advanced)

figure('Name','Identification in the moments'),
subplot(211)
mmm = (ide_strength_J);
[ss, is] = sort(mmm);
    bar(log([ide_strength_J(:,is)' ide_strength_J_prior(:,is)']))
set(gca,'xlim',[0 nparam+1])
set(gca,'xticklabel','')
dy = get(gca,'ylim');
for ip=1:nparam,
    text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
legend('relative to param value','relative to prior std','Location','Best')
if  flag_score,
    title('Identification strength using asymptotic Information matrix (log-scale)')
else
    title('Identification strength in the moments (log-scale)')
end
subplot(212)

if SampleSize > 1,
    mmm = mean(siJnorm);
else
    mmm = (siJnorm);
end
bar(mmm(is))
set(gca,'xlim',[0 nparam+1])
set(gca,'xticklabel','')
dy = get(gca,'ylim');
for ip=1:nparam,
    text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Sensitivity in the moments')

if advanced,
    figure('Name','Identification in the model'),
    subplot(211)
    if SampleSize > 1,
        mmm = mean(siLREnorm);
    else
        mmm = (siLREnorm);
    end
    mmm = [NaN(1, offset), mmm];
    bar(mmm(is))
    set(gca,'xticklabel','')
    for ip=1:nparam,
        text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
    end
    title('Sensitivity in the LRE model')
        
    subplot(212)
    if SampleSize>1,
        mmm = mean(siHnorm);
    else
        mmm = (siHnorm);
    end
    bar(mmm(is))
    set(gca,'xticklabel','')
    dy = get(gca,'ylim');
    for ip=1:nparam,
        text(ip,dy(1),name{is(ip)},'rotation',90,'HorizontalAlignment','right','interpreter','none')
    end
    title('Sensitivity in the model')

    if options_.nograph, close(gcf); end

    % identificaton patterns
    for  j=1:size(cosnJ,2),
        pax=NaN(nparam,nparam);
        fprintf('\n\n')
        disp(['Collinearity patterns with ', int2str(j) ,' parameter(s)'])
        fprintf('%-15s [%-*s] %10s\n','Parameter',(15+1)*j,' Expl. params ','cosn')
        for i=1:nparam,
            namx='';
            for in=1:j,
                namx=[namx ' ' sprintf('%-15s',name{pars{i,j}(in)})];
            end
            fprintf('%-15s [%s] %10.3f\n',name{i},namx,cosnJ(i,j))
            pax(i,pars{i,j})=cosnJ(i,j);
        end
        figure('name',['Collinearity patterns with ', int2str(j) ,' parameter(s)']), 
        imagesc(pax,[0 1]);
        set(gca,'xticklabel','')
        set(gca,'yticklabel','')
        for ip=1:nparam,
            text(ip,(0.5),name{ip},'rotation',90,'HorizontalAlignment','left','interpreter','none')
            text(0.5,ip,name{ip},'rotation',0,'HorizontalAlignment','right','interpreter','none')
        end
        colorbar;
        ax=colormap;
        ax(1,:)=[0.9 0.9 0.9];
        colormap(ax);
        if nparam>10,
            set(gca,'xtick',(5:5:nparam))
            set(gca,'ytick',(5:5:nparam))
        end
        set(gca,'xgrid','on')
        set(gca,'ygrid','on')
        saveas(gcf,[IdentifDirectoryName,'/',M_.fname,'_ident_collinearity_', int2str(j)])
        eval(['print -depsc2 ' IdentifDirectoryName '/' M_.fname '_ident_collinearity_', int2str(j)]);
        eval(['print -dpdf ' IdentifDirectoryName '/' M_.fname '_ident_collinearity_', int2str(j)]);
        if options_.nograph, close(gcf); end
    end
    disp('')
    if flag_score,
        [U,S,V]=svd(AHess,0);
        if nparam<5,
            f1 = figure('name','Identification patterns (Information matrix)');
        else
            f1 = figure('name','Identification patterns (Information matrix): SMALLEST SV');
            f2 = figure('name','Identification patterns (Information matrix): HIGHEST SV');
        end
    else
        [U,S,V]=svd(siJ./normJ(:,ones(nparam,1)),0);
        if nparam<5,
            f1 = figure('name','Identification patterns (moments)');
        else
            f1 = figure('name','Identification patterns (moments): SMALLEST SV');
            f2 = figure('name','Identification patterns (moments): HIGHEST SV');
        end
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
        else
            bar(abs(V(:,jj))),
            Stit = S(jj,jj);
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
    
    if SampleSize>1,
        options_.nograph=1;
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
        ncut=floor(SampleSize/10*9);
        [~,is]=sort(idelre.cond);
        [proba, dproba] = stab_map_1(pdraws, is(1:ncut), is(ncut+1:end), 'HighestCondNumberLRE', 1, [], IdentifDirectoryName);
        [~,is]=sort(idemodel.cond);
        [proba, dproba] = stab_map_1(pdraws, is(1:ncut), is(ncut+1:end), 'HighestCondNumberModel', 1, [], IdentifDirectoryName);
        [~,is]=sort(idemoments.cond);
        [proba, dproba] = stab_map_1(pdraws, is(1:ncut), is(ncut+1:end), 'HighestCondNumberMoments', 1, [], IdentifDirectoryName);
%         [proba, dproba] = stab_map_1(idemoments.Mco', is(1:ncut), is(ncut+1:end), 'HighestCondNumberMoments_vs_Mco', 1, [], IdentifDirectoryName);
        for j=1:nparam,
%             ibeh=find(idemoments.Mco(j,:)<0.9);
%             inonbeh=find(idemoments.Mco(j,:)>=0.9);
%             if ~isempty(ibeh) && ~isempty(inonbeh)
%                 [proba, dproba] = stab_map_1(pdraws, ibeh, inonbeh, ['HighestMultiCollinearity_',name{j}], 1, [], IdentifDirectoryName);
%             end
            [~,is]=sort(idemoments.Mco(j,:));
            [proba, dproba] = stab_map_1(pdraws, is(1:ncut), is(ncut+1:end), ['HighestMultiCollinearity_',name{j}], 1, [], IdentifDirectoryName);
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
