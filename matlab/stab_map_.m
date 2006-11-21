function x0 = stab_map_(Nsam, fload, ksstat, alpha2, prepSA, pprior, ilptau, OutputDirectoryName)
%
% function x0 = stab_map_(Nsam, fload, alpha2, prepSA, pprior)
%
% Mapping of stability regions in the prior ranges applying
% Monte Carlo filtering techniques.
%
% M. Ratto, Global Sensitivity Analysis for Macroeconomic models
% I. Mapping stability, MIMEO, 2005.
%
% INPUTS
% Nsam = MC sample size
% fload = 0 to run new MC; 1 to load prevoiusly generated analysis
% alpha2 =  significance level for bivariate sensitivity analysis
% [abs(corrcoef) > alpha2]
% prepSA = 1: save transition matrices for mapping reduced form
%        = 0: no transition matrix saved (default)
% pprior = 1: sample from prior ranges (default): sample saved in
%            _prior.mat   file
%        = 0: sample from posterior ranges: sample saved in
%            _mc.mat file
% OUTPUT: 
% x0: one parameter vector for which the model is stable.
%
% GRAPHS
% 1) Pdf's of marginal distributions under the stability (dotted
%     lines) and unstability (solid lines) regions
% 2) Cumulative distributions of: 
%   - stable subset (dotted lines) 
%   - unacceptable subset (solid lines)
% 3) Bivariate plots of significant correlation patterns 
%  ( abs(corrcoef) > alpha2) under the stable and unacceptable subsets
%
% USES lptauSEQ, 
%      stab_map_1, stab_map_2
%
% Copyright (C) 2005 Marco Ratto
% THIS PROGRAM WAS WRITTEN FOR MATLAB BY
% Marco Ratto,
% Unit of Econometrics and Statistics AF
% (http://www.jrc.cec.eu.int/uasa/),
% IPSC, Joint Research Centre
% The European Commission,
% TP 361, 21020 ISPRA(VA), ITALY
% marco.ratto@jrc.it 
%
% ALL COPIES MUST BE PROVIDED FREE OF CHARGE AND MUST INCLUDE THIS COPYRIGHT
% NOTICE.
%

%global bayestopt_ estim_params_ dr_ options_ ys_ fname_
global bayestopt_ estim_params_ options_ oo_ M_

opt_gsa=options_.opt_gsa;

nliv = opt_gsa.morris_nliv;    
ntra = opt_gsa.morris_ntra;

dr_ = oo_.dr;
if isfield(dr_,'ghx'),
  ys_ = oo_.dr.ys;
  nspred = size(dr_.ghx,2);
  nboth = dr_.nboth;
  nfwrd = dr_.nfwrd;
end
fname_ = M_.fname;

nshock = estim_params_.nvx;
nshock = nshock + estim_params_.nvn;
nshock = nshock + estim_params_.ncx;
nshock = nshock + estim_params_.ncn;
lpmat0=[];

if nargin==0,
  Nsam=2000; %2^13; %256;
end
if nargin<2,
  fload=0;
end
if nargin<3,
  ksstat=0.1;
end
if nargin<4,
  alpha2=0.3;
end
if nargin<5,
  prepSA=0;
end
if nargin<6,
  pprior=1;
end
if nargin<7,
  ilptau=1;
end
if nargin<8,
  OutputDirectoryName='';
end

options_.periods=0;
options_.nomoments=1;
options_.irf=0;
options_.noprint=1;
options_.simul=0;

if fload==0 | nargin<2 | isempty(fload),
  if prepSA
    T=zeros(size(dr_.ghx,1),size(dr_.ghx,2)+size(dr_.ghu,2),Nsam/2);
  end
  
  if isfield(dr_,'ghx'),
    egg=zeros(length(dr_.eigval),Nsam);
  end
  yys=zeros(length(dr_.ys),Nsam);
  
  if opt_gsa.morris
    if opt_gsa.morris == 1
      [lpmat, OutFact] = Sampling_Function_2(nliv, estim_params_.np, ntra, ones(estim_params_.np, 1), zeros(estim_params_.np,1), []);
      lpmat = lpmat.*(nliv-1)/nliv+1/nliv/2;
      Nsam=size(lpmat,1);
    elseif opt_gsa.morris==2
      lpmat = prep_ide(Nsam,estim_params_.np,5);
      Nsam=size(lpmat,1);
    end
  else
  if estim_params_.np<52 & ilptau>0,
    [lpmat] = lptauSEQ(Nsam,estim_params_.np); % lptau
    if estim_params_.np>30 | ilptau==2, % scrambled lptau
      for j=1:estim_params_.np,
        lpmat(:,j)=lpmat(randperm(Nsam),j);
      end
    end
  else ilptau==0
    %[lpmat] = rand(Nsam,estim_params_.np);
    for j=1:estim_params_.np,
      lpmat(:,j) = randperm(Nsam)'./(Nsam+1); %latin hypercube
    end

  end
  end
  prior_draw_gsa(1);
  if pprior,
    for j=1:nshock,
      lpmat0(:,j) = randperm(Nsam)'./(Nsam+1); %latin hypercube
      %lpmat0(:,j)=lpmat0(:,j).*(bayestopt_.ub(j)-bayestopt_.lb(j))+bayestopt_.lb(j);
    end
    %for j=1:estim_params_.np,
    %  lpmat(:,j)=lpmat(:,j).*(bayestopt_.ub(j+nshock)-bayestopt_.lb(j+nshock))+bayestopt_.lb(j+nshock);
    %end
    xx=prior_draw_gsa(0,[lpmat0 lpmat]);
    lpmat0=xx(:,1:nshock);
    lpmat=xx(:,nshock+1:end);
    clear xx;
  else
    %         for j=1:nshock,
    %             xparam1(j) = oo_.posterior_mode.shocks_std.(bayestopt_.name{j});
    %             sd(j) = oo_.posterior_std.shocks_std.(bayestopt_.name{j});
    %             lpmat0(:,j) = randperm(Nsam)'./(Nsam+1); %latin hypercube
    %             lb = max(bayestopt_.lb(j), xparam1(j)-2*sd(j));
    %             ub1=xparam1(j)+(xparam1(j) - lb); % define symmetric range around the mode!
    %             ub = min(bayestopt_.ub(j),ub1);
    %             if ub<ub1,
    %                 lb=xparam1(j)-(ub-xparam1(j)); % define symmetric range around the mode!
    %             end
    %             lpmat0(:,j) = lpmat0(:,j).*(ub-lb)+lb;
    %         end
    %         % 
    %         for j=1:estim_params_.np,
    %             xparam1(j+nshock) = oo_.posterior_mode.parameters.(bayestopt_.name{j+nshock});
    %             sd(j+nshock) = oo_.posterior_std.parameters.(bayestopt_.name{j+nshock});
    %             lb = max(bayestopt_.lb(j+nshock),xparam1(j+nshock)-2*sd(j+nshock));
    %             ub1=xparam1(j+nshock)+(xparam1(j+nshock) - lb); % define symmetric range around the mode!
    %             ub = min(bayestopt_.ub(j+nshock),ub1);
    %             if ub<ub1,
    %                 lb=xparam1(j+nshock)-(ub-xparam1(j+nshock)); % define symmetric range around the mode!
    %             end
    %             %ub = min(bayestopt_.ub(j+nshock),xparam1(j+nshock)+2*sd(j+nshock));
    %             if estim_params_.np>30 & estim_params_.np<52
    %                 lpmat(:,j) = lpmat(randperm(Nsam),j).*(ub-lb)+lb;
    %             else
    %                 lpmat(:,j) = lpmat(:,j).*(ub-lb)+lb;
    %             end
    %         end
    %load([fname_,'_mode'])  
    eval(['load ' options_.mode_file ';']');
    d = chol(inv(hh));
    lp=randn(Nsam*2,nshock+estim_params_.np)*d+kron(ones(Nsam*2,1),xparam1');
    for j=1:Nsam*2,
        lnprior(j) = any(lp(j,:)'<=bayestopt_.lb | lp(j,:)'>=bayestopt_.ub);
    end
    ireal=[1:2*Nsam]; 
    ireal=ireal(find(lnprior==0));
    lp=lp(ireal,:);
    Nsam=min(Nsam, length(ireal));
    lpmat0=lp(1:Nsam,1:nshock);
    lpmat=lp(1:Nsam,nshock+1:end);
    clear lp lnprior ireal;
  end
  % 
  h = waitbar(0,'Please wait...');
  istable=[1:Nsam];
  jstab=0;
  iunstable=[1:Nsam];
  iindeterm=zeros(1,Nsam);
  iwrong=zeros(1,Nsam);
  for j=1:Nsam,
    M_.params(estim_params_.param_vals(:,1)) = lpmat(j,:)';
    stoch_simul([]);
    dr_ = oo_.dr;
    if isfield(dr_,'ghx'),
      egg(:,j) = sort(dr_.eigval);
      iunstable(j)=0;
      if prepSA
        jstab=jstab+1;
        T(:,:,jstab) = [dr_.ghx dr_.ghu];
      end            
      if ~exist('nspred'),
        nspred = size(dr_.ghx,2);
        nboth = dr_.nboth;
        nfwrd = dr_.nfwrd;
      end
    else
      istable(j)=0;
      if isfield(dr_,'eigval')
        egg(:,j) = sort(dr_.eigval);
        if exist('nspred')
        if any(isnan(egg(1:nspred,j)))
          iwrong(j)=j;
        else
          if (nboth | nfwrd) & abs(egg(nspred+1,j))<=options_.qz_criterium,
            iindeterm(j)=j;
          end                                      
        end  
        end
      else
        if exist('egg'),
        egg(:,j)=ones(size(egg,1),1).*1.1;
        end
        iwrong(j)=j;
      end
    end
    ys_=real(dr_.ys);
    yys(:,j) = ys_;
    ys_=yys(:,1);
    waitbar(j/Nsam,h,['MC iteration ',int2str(j),'/',int2str(Nsam)])
  end
  close(h)
  if prepSA,
    T=T(:,:,1:jstab);
  end
  istable=istable(find(istable));  % stable params
  iunstable=iunstable(find(iunstable));   % unstable params
  iindeterm=iindeterm(find(iindeterm));  % indeterminacy
  iwrong=iwrong(find(iwrong));  % dynare could not find solution
  
  %     % map stable samples
  %     istable=[1:Nsam];
  %     for j=1:Nsam,
  %         if any(isnan(egg(1:nspred,j)))
  %             istable(j)=0;
  %         else
  %             if abs(egg(nspred,j))>=options_.qz_criterium; %(1-(options_.qz_criterium-1)); %1-1.e-5;
  %                 istable(j)=0;
  %                 %elseif (dr_.nboth | dr_.nfwrd) & abs(egg(nspred+1,j))<=options_.qz_criterium; %1+1.e-5;
  %             elseif (nboth | nfwrd) & abs(egg(nspred+1,j))<=options_.qz_criterium; %1+1.e-5;
  %                 istable(j)=0;
  %             end
  %         end
  %     end
  %     istable=istable(find(istable));  % stable params
  %     
  %     % map unstable samples
  %     iunstable=[1:Nsam];
  %     for j=1:Nsam,
  %         %if abs(egg(dr_.npred+1,j))>1+1.e-5 & abs(egg(dr_.npred,j))<1-1.e-5;
  %         %if (dr_.nboth | dr_.nfwrd),
  %         if ~any(isnan(egg(1:5,j)))
  %             if (nboth | nfwrd),
  %                 if abs(egg(nspred+1,j))>options_.qz_criterium & abs(egg(nspred,j))<options_.qz_criterium; %(1-(options_.qz_criterium-1));
  %                     iunstable(j)=0;
  %                 end
  %             else
  %                 if abs(egg(nspred,j))<options_.qz_criterium; %(1-(options_.qz_criterium-1));
  %                     iunstable(j)=0;
  %                 end
  %             end
  %         end
  %     end
  %     iunstable=iunstable(find(iunstable));   % unstable params
  if pprior,
    if ~prepSA
      save([OutputDirectoryName '\' fname_ '_prior'],'lpmat','lpmat0','iunstable','istable','iindeterm','iwrong','egg','yys','nspred','nboth','nfwrd')
    else
      save([OutputDirectoryName '\' fname_ '_prior'],'lpmat','lpmat0','iunstable','istable','iindeterm','iwrong','egg','yys','T','nspred','nboth','nfwrd')
    end
    
  else
    if ~prepSA
      save([OutputDirectoryName '\' fname_ '_mc'], ...
        'lpmat','lpmat0','iunstable','istable','iindeterm','iwrong','egg','yys','nspred','nboth','nfwrd')
    else
      save([OutputDirectoryName '\' fname_ '_mc'], ...
        'lpmat','lpmat0','iunstable','istable','iindeterm','iwrong','egg','yys','T','nspred','nboth','nfwrd')
    end
  end
else
  if pprior,
    filetoload=[OutputDirectoryName '\' fname_ '_prior'];
  else
    filetoload=[OutputDirectoryName '\' fname_ '_mc'];
  end
  load(filetoload,'lpmat','lpmat0','iunstable','istable','iindeterm','iwrong','egg','yys','nspred','nboth','nfwrd')
  Nsam = size(lpmat,1);    
  
  if prepSA & isempty(strmatch('T',who('-file', filetoload),'exact')),
    h = waitbar(0,'Please wait...');
    options_.periods=0;
    options_.nomoments=1;
    options_.irf=0;
    options_.noprint=1;
    stoch_simul([]);
    T=zeros(size(dr_.ghx,1),size(dr_.ghx,2)+size(dr_.ghu,2),length(istable));
    ntrans=length(istable);
    for j=1:ntrans,
      M_.params(estim_params_.param_vals(:,1)) = lpmat(istable(j),:)';
      stoch_simul([]);
      dr_ = oo_.dr;
      T(:,:,j) = [dr_.ghx dr_.ghu];
      if ~exist('nspred')
        nspred = size(dr_.ghx,2);
        nboth = dr_.nboth;
        nfwrd = dr_.nfwrd;
      end
      ys_=real(dr_.ys);
      yys(:,j) = ys_;
      ys_=yys(:,1);
      waitbar(j/ntrans,h,['MC iteration ',int2str(j),'/',int2str(ntrans)])
    end
    close(h)
    save(filetoload,'T','-append')    
  else
    load(filetoload,'T')        
  end
end

if pprior
  aname='prior_stab';
  auname='prior_unacceptable';
  aunstname='prior_unstable';
  aindname='prior_indeterm';
  asname='prior_stable';
else
  aname='mc_stab';
  auname='mc_unacceptable';
  aunstname='mc_unstable';
  aindname='mc_indeterm';
  asname='mc_stable';
end
delete([OutputDirectoryName,'\',fname_,'_',aname,'_*.*']);
%delete([OutputDirectoryName,'\',fname_,'_',aname,'_SA_*.*']);
delete([OutputDirectoryName,'\',fname_,'_',asname,'_corr_*.*']);
delete([OutputDirectoryName,'\',fname_,'_',auname,'_corr_*.*']);
delete([OutputDirectoryName,'\',fname_,'_',aunstname,'_corr_*.*']);
delete([OutputDirectoryName,'\',fname_,'_',aindname,'_corr_*.*']);

if length(iunstable)>0 & length(iunstable)<Nsam,
  disp([num2str(length(istable)/Nsam*100),'\% of the prior support is stable.'])
  disp([num2str( (length(iunstable)-length(iwrong)-length(iindeterm) )/Nsam*100),'\% of the prior support is unstable.'])
  if ~isempty(iindeterm),
    disp([num2str(length(iindeterm)/Nsam*100),'\% of the prior support gives indeterminacy.'])
  end
  if ~isempty(iwrong),
    disp(' ');
    disp(['For ',num2str(length(iwrong)/Nsam*100),'\% of the prior support dynare could not find a solution.'])      
  end
  % Blanchard Kahn
  [proba, dproba] = stab_map_1(lpmat, istable, iunstable, aname,0);
  indstab=find(dproba>ksstat);
  disp('The following parameters mostly drive acceptable behaviour')
  disp(M_.param_names(estim_params_.param_vals(indstab,1),:))  
  stab_map_1(lpmat, istable, iunstable, aname, 1, indstab, OutputDirectoryName);
  if ~isempty(iindeterm),
    ixun=iunstable(find(~ismember(iunstable,[iindeterm,iwrong])));
    [proba, dproba] = stab_map_1(lpmat, [1:Nsam], iindeterm, [aname, '_indet'],0);
    indindet=find(dproba>ksstat);
    disp('The following parameters mostly drive indeterminacy')
    disp(M_.param_names(estim_params_.param_vals(indindet,1),:))  
    stab_map_1(lpmat, [1:Nsam], iindeterm, [aname, '_indet'], 1, indindet, OutputDirectoryName);
    if ~isempty(ixun),
      [proba, dproba] = stab_map_1(lpmat, [1:Nsam], ixun, [aname, '_unst'],0);
      indunst=find(dproba>ksstat);
      disp('The following parameters mostly drive instability')
      disp(M_.param_names(estim_params_.param_vals(indunst,1),:))  
      stab_map_1(lpmat, [1:Nsam], ixun, [aname, '_unst'], 1, indunst, OutputDirectoryName);
    end
  end
  
  disp(' ')
  disp('Starting bivariate analysis:')
  
  c0=corrcoef(lpmat(istable,:));
  c00=tril(c0,-1);
  
  stab_map_2(lpmat(istable,:),alpha2, asname, OutputDirectoryName);
  stab_map_2(lpmat(iunstable,:),alpha2, auname, OutputDirectoryName);
  if ~isempty(iindeterm),
    stab_map_2(lpmat(iindeterm,:),alpha2, aindname, OutputDirectoryName);
    if ~isempty(ixun),
      stab_map_2(lpmat(ixun,:),alpha2, aunstname, OutputDirectoryName);
    end
  end
  
  x0=0.5.*(bayestopt_.ub(1:nshock)-bayestopt_.lb(1:nshock))+bayestopt_.lb(1:nshock);
  x0 = [x0; lpmat(istable(1),:)'];
  if istable(end)~=Nsam
    M_.params(estim_params_.param_vals(:,1)) = lpmat(istable(1),:)';
    stoch_simul([]);        
  end
else
  if length(iunstable)==0,
    disp('All parameter values in the specified ranges are stable!')
    x0=0.5.*(bayestopt_.ub(1:nshock)-bayestopt_.lb(1:nshock))+bayestopt_.lb(1:nshock);
    x0 = [x0; lpmat(istable(1),:)'];
  else
    disp('All parameter values in the specified ranges are not acceptable!')        
    x0=[];
  end
  
end

if opt_gsa.redform,
  [nr,nc,nn]=size(T);
  j0=0;
  for j=1:nr,
    for i=1:nc,
      y0=squeeze(T(j,i,:));
      if max(y0)-min(y0)>1.e-10,
        j0=j0+1;
        y1=ones(size(lpmat,1),1)*NaN;
        y1(istable,1)=y0;
        yt(:,j0)=y1;
      end
    end
  end
  
end

if opt_gsa.morris==1 & opt_gsa.redform,
  OutputDir = CheckPath('GSA');  

  for j=1:j0,
    [SAmeas, SAMorris(:,:,j)] = Morris_Measure_Groups(estim_params_.np, lpmat, yt(:,j),nliv);
  end

  SAM = squeeze(SAMorris(:,1,:));
  for j=1:j0
    SAnorm(:,j)=SAM(:,j)./max(SAM(:,j));
    irex(j)=length(find(SAnorm(:,j)>0.01));
  end
  [dum, irel]=sort(irex);

  figure, bar(SAnorm(:,irel))
  set(gca,'xtick',[1:estim_params_.np])
  set(gca,'xlim',[0.5 estim_params_.np+0.5])
  title('Elementary effects parameters')
  saveas(gcf,[OutputDir,'\',fname_,'_morris_par'])
  eval(['print -depsc2 ' OutputDir '\' fname_ '_morris_par']);
  eval(['print -dpdf ' OutputDir '\' fname_ '_morris_par']);
  figure, bar(SAnorm(:,irel)')
  set(gca,'xtick',[1:j0])
  set(gca,'xlim',[0.5 j0+0.5])
  title('Elementary effects relationships')
  saveas(gcf,[OutputDir,'\',fname_,'_morris_redform'])
  eval(['print -depsc2 ' OutputDir '\' fname_ '_morris_redform']);
  eval(['print -dpdf ' OutputDir '\' fname_ '_morris_redform']);
elseif opt_gsa.morris==2 & opt_gsa.redform,
  np=estim_params_.np;
  na=(4*np+1)*opt_gsa.Nsam;
  for j=1:j0,
    [idex(j,:), yd(j,:)] = spop_ide(lpmat, yt(:,j), opt_gsa.Nsam, 5-1);
  end
  iok=find(~isnan(yt(1:opt_gsa.Nsam,1)));
  yr=NaN*ones(size(lpmat,1),j0);
  for j=1:j0, 
    ys(j,:)=yd(j,:)./max(yd(j,:));
    [dum, is]=sort(yt(iok,j));
    yr(iok(is),j)=[1:length(iok)]'./length(iok);
    yr(istable(length(iok)+1:end),j) = interp1(yt(iok,j),yr(iok,j),yt(istable(length(iok)+1:end),j),'','extrap');
    ineg=find(yr(:,j)<0);
    if any(ineg),
      [dum, is]=sort(yr(ineg,j));
      yr(ineg(is),j)=-[length(ineg):-1:1]./length(iok);
      
    end
    [idex_r(j,:), yd_r(j,:)] = spop_ide(lpmat, yr(:,j), opt_gsa.Nsam, 5-1);
    ys_r(j,:)=yd_r(j,:)./max(yd_r(j,:));
      
  end,
  figure, bar((idex.*ys)./opt_gsa.Nsam), title('Relationships')
  figure, bar((idex.*ys)'./opt_gsa.Nsam), title('Parameters')
  figure, bar((idex_r.*ys_r)./opt_gsa.Nsam), title('Relationships rank')
  figure, bar((idex_r.*ys_r)'./opt_gsa.Nsam), title('Parameters rank')
  [v0,d0]=eig(corrcoef(yt(iok,:)));
  ee=diag(d0);
  ee=ee([end:-1:1])./j0;
  i0=length(find(ee>0.01));
  v0=v0(:,[end:-1:1]);
  for j=1:i0,
    [idex_pc(j,:), yd_pc(j,:)] = spop_ide(lpmat, yt*v0(:,j), opt_gsa.Nsam, 5-1);
  end
  for j=1:i0, 
    ys_pc(j,:)=yd_pc(j,:)./max(yd_pc(j,:));
  end,
  figure, bar((idex_pc.*ys_pc)./opt_gsa.Nsam), title('Relationships PCA')
  figure, bar((idex_pc.*ys_pc)'./opt_gsa.Nsam), title('Parameters PCA')
  
  [vr,dr]=eig(corrcoef(yr(iok,:)));
  er=diag(dr);
  er=er([end:-1:1])./j0;
  ir0=length(find(er>0.01));
  vr=vr(:,[end:-1:1]);
  for j=1:ir0,
    [idex_pcr(j,:), yd_pcr(j,:)] = spop_ide(lpmat, yr*vr(:,j), opt_gsa.Nsam, 5-1);
  end
  for j=1:ir0, 
    ys_pcr(j,:)=yd_pcr(j,:)./max(yd_pcr(j,:));
  end,
  figure, bar((idex_pcr.*ys_pcr)./opt_gsa.Nsam), title('Relationships rank PCA')
  figure, bar((idex_pcr.*ys_pcr)'./opt_gsa.Nsam), title('Parameters rank PCA')
elseif opt_gsa.redform
  yr=[];
  for j=1:j0,
    [dum, is]=sort(yt(istable,j));
    yr(is,j)=[1:length(istable)]'./length(istable);
  end
  gsa_flag=0; %-2 for GSA estimation
  RividDir = CheckPath('GSA/rivid');
  for j=1:j0,
    gsa_(j) = gsa_sdp_fn(yt(istable,j), lpmat(istable,:), ...
      bayestopt_.pshape, [bayestopt_.p1 bayestopt_.p2 bayestopt_.p3 bayestopt_.p4], ...
      gsa_flag, [RividDir,'/map_T',int2str(j)], M_.param_names);
    gsa_r(j) = gsa_sdp_fn(yr(:,j), lpmat(istable,:), ...
      bayestopt_.pshape, [bayestopt_.p1 bayestopt_.p2 bayestopt_.p3 bayestopt_.p4], ...
      gsa_flag, [RividDir,'/map_T_rank',int2str(j)], M_.param_names);
    S(j,:)=gsa_(j).multivariate.si; 
    Sr(j,:)=gsa_r(j).multivariate.si; 
  end
  figure, bar(S')
  set(gca,'xtick',[1:estim_params_.np])
  set(gca,'xlim',[0.5 estim_params_.np+0.5])
  title('Main effects parameters')
  saveas(gcf,[RividDir,'\',fname_,'_SA_par'])
  eval(['print -depsc2 ' RividDir '\' fname_ '_SA_par']);
  eval(['print -dpdf ' RividDir '\' fname_ '_SA_par']);
  figure, bar(S)
  set(gca,'xtick',[1:j0])
  set(gca,'xlim',[0.5 j0+0.5])
  title('Main effects relationships')
  saveas(gcf,[RividDir,'\',fname_,'_SA_redform'])
  eval(['print -depsc2 ' RividDir '\' fname_ '_SA_redform']);
  eval(['print -dpdf ' RividDir '\' fname_ '_SA_redform']);
  
  figure, bar(Sr')
  set(gca,'xtick',[1:estim_params_.np])
  set(gca,'xlim',[0.5 estim_params_.np+0.5])
  title('Main effects parameters rank')
  saveas(gcf,[RividDir,'\',fname_,'_SA_par'])
  eval(['print -depsc2 ' RividDir '\' fname_ '_SA_par_rank']);
  eval(['print -dpdf ' RividDir '\' fname_ '_SA_par']);
  figure, bar(Sr)
  set(gca,'xtick',[1:j0])
  set(gca,'xlim',[0.5 j0+0.5])
  title('Main effects relationships rank')
  saveas(gcf,[RividDir,'\',fname_,'_SA_redform_rank'])
  eval(['print -depsc2 ' RividDir '\' fname_ '_SA_redform_rank']);
  eval(['print -dpdf ' RividDir '\' fname_ '_SA_redform_rank']);

  [v0,d0]=eig(corrcoef(yt(istable,:)));
  ee=diag(d0);
  ee=ee(end:-1:1);
  v0=v0(:,end:-1:1);
  %id=find((er./length(er))>0.01);
  id=find(cumsum(ee)./j0<0.99);
  jpc=length(id)+1;
  gsa_flag=-2;
  for j=1:jpc,
    gsa_pca(j) = gsa_sdp_fn(yt(istable,:)*v0(:,j), lpmat(istable,:), ...
      bayestopt_.pshape, [bayestopt_.p1 bayestopt_.p2 bayestopt_.p3 bayestopt_.p4], ...
      gsa_flag, [RividDir,'/map_T_pca',int2str(j)], M_.param_names);
    S_pca(j,:)=gsa_pca(j).multivariate.si;
  end
  figure, bar(S_pca')
  set(gca,'xtick',[1:estim_params_.np])
  set(gca,'xlim',[0.5 estim_params_.np+0.5])
  title('Main effects parameters PCA')
  saveas(gcf,[RividDir,'\',fname_,'_SA_par_PCA'])
  eval(['print -depsc2 ' RividDir '\' fname_ '_SA_par_PCA']);
  eval(['print -dpdf ' RividDir '\' fname_ '_SA_par_PCA']);
  figure, bar(S_pca)
  set(gca,'xtick',[1:j0])
  set(gca,'xlim',[0.5 j0+0.5])
  title('Main effects relationships PCA')
  saveas(gcf,[RividDir,'\',fname_,'_SA_redform_PCA'])
  eval(['print -depsc2 ' RividDir '\' fname_ '_SA_redform_PCA']);
  eval(['print -dpdf ' RividDir '\' fname_ '_SA_redform_PCA']);
  
end
