function x0 = stab_map_(Nsam, fload, alpha2, prepSA, pprior, ilptau)
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


if nargin==0,
  Nsam=2000; %2^13; %256;
end
if nargin<2,
  fload=0;
end
if nargin<3,
  alpha2=0.3;
end
if nargin<4,
  prepSA=0;
end
if nargin<5,
  pprior=1;
end
if nargin<6,
  ilptau=1;
end

options_.periods=0;
options_.nomoments=1;
options_.irf=0;
options_.noprint=1;

if fload==0 | nargin<2 | isempty(fload),
  if estim_params_.np<52 & ilptau,
    [lpmat] = lptauSEQ(Nsam,estim_params_.np);
    if estim_params_.np>30
      for j=1:estim_params_.np,
        lpmat(:,j)=lpmat(randperm(Nsam),j).*(bayestopt_.ub(j+nshock)-bayestopt_.lb(j+nshock))+bayestopt_.lb(j+nshock);
      end
    end
  else
    %[lpmat] = rand(Nsam,estim_params_.np);
    for j=1:estim_params_.np,
      lpmat(:,j) = randperm(Nsam)'./(Nsam+1); %latin hypercube
    end
  end
  
  if pprior,
    for j=1:nshock,
      lpmat0(:,j) = randperm(Nsam)'./(Nsam+1); %latin hypercube
      lpmat0(:,j)=lpmat0(:,j).*(bayestopt_.ub(j)-bayestopt_.lb(j))+bayestopt_.lb(j);
    end
    for j=1:estim_params_.np,
      lpmat(:,j)=lpmat(:,j).*(bayestopt_.ub(j+nshock)-bayestopt_.lb(j+nshock))+bayestopt_.lb(j+nshock);
    end
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
    lp=randn(Nsam,nshock+estim_params_.np)*d+kron(ones(Nsam,1),xparam1');
    lpmat0=lp(:,1:nshock);
    lpmat=lp(:,nshock+1:end);
  end
  % 
  h = waitbar(0,'Please wait...');
  istable=[1:Nsam];
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
        T(:,:,j) = [dr_.ghx dr_.ghu];
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
        if any(isnan(egg(1:nspred,j)))
          iwrong(j)=j;
        else
          if (nboth | nfwrd) & abs(egg(nspred+1,j))<=options_.qz_criterium,
            iindeterm(j)=j;
          end                                      
        end  
      else
        egg(:,j)=ones(size(egg,1),1).*1.1;
        iwrong(j)=j;
      end
    end
    ys_=real(dr_.ys);
    yys(:,j) = ys_;
    ys_=yys(:,1);
    waitbar(j/Nsam,h,['MC iteration ',int2str(j),'/',int2str(Nsam)])
  end
  close(h)
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
      save([fname_,'_prior'],'lpmat','lpmat0','iunstable','istable','iindeterm','iwrong','egg','yys','nspred','nboth','nfwrd')
    else
      save([fname_,'_prior'],'lpmat','lpmat0','iunstable','istable','iindeterm','iwrong','egg','yys','T','nspred','nboth','nfwrd')
    end
    
  else
    if ~prepSA
      save([fname_,'_mc'],'lpmat','lpmat0','iunstable','istable','iindeterm','iwrong','egg','yys','nspred','nboth','nfwrd')
    else
      save([fname_,'_mc'],'lpmat','lpmat0','iunstable','istable','iindeterm','iwrong','egg','yys','T','nspred','nboth','nfwrd')
    end
  end
else
  if pprior,
    load([fname_,'_prior'])
  else
    load([fname_,'_mc'])
  end
  Nsam = size(lpmat,1);    
end

if prepSA & ~exist('T'),
  h = waitbar(0,'Please wait...');
  options_.periods=0;
  options_.nomoments=1;
  options_.irf=0;
  options_.noprint=1;
  stoch_simul([]);
  T=zeros(size(dr_.ghx,1),size(dr_.ghx,2)+size(dr_.ghu,2),length(istable));
  
  for j=1:length(istable),
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
    waitbar(j/Nsam,h,['MC iteration ',int2str(j),'/',int2str(Nsam)])
  end
  close(h)
  if pprior
    save([fname_,'_prior'],'T','-append')    
  else
    save([fname_,'_mc'],'T','-append')    
  end
end

if pprior
  aname='prior_stab';
  auname='prior_unacceptable';
  asname='prior_stable';
else
  aname='mc_stab';
  auname='mc_unacceptable';
  asname='mc_stable';
end
delete([fname_,'_',aname,'_*.*']);
delete([fname_,'_',aname,'_SA_*.*']);
delete([fname_,'_',asname,'_corr_*.*']);
delete([fname_,'_',auname,'_corr_*.*']);

if length(iunstable)>0 & length(iunstable)<Nsam,
  disp([num2str(length(istable)/Nsam*100),'\% of the prior support is stable.'])
  if ~isempty(iwrong),
    disp(['For ',num2str(length(iwrong)/Nsam*100),'\% of the prior support dynare could not find a solution.'])      
  end
  if ~isempty(iindeterm),
    disp([num2str(length(iindeterm)/Nsam*100),'\% of the prior support gives indeterminacy.'])
  end
  % Blanchard Kahn
  proba = stab_map_1(lpmat, istable, iunstable, aname);
  disp(' ')
  disp(' ')
  disp('Starting bivariate analysis:')
  
  c0=corrcoef(lpmat(istable,:));
  c00=tril(c0,-1);
  
  stab_map_2(lpmat(istable,:),alpha2, asname);
  stab_map_2(lpmat(iunstable,:),alpha2, auname);
  
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
