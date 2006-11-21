function [fit,fitse,par,parse,zs,pars,parses,rsq,nvr,y0]=...
    sdp0(y,z,x,TVP,nvr,opts,P0,x0,nvr0,tab)
% SDP  Non-parametric state dependent regression modelling (backfitting)
%
% [fit,fitse,par,parse,xs,pars,parses,rsq,nvre,y0]=...
%                    sdp(y,z,x,TVP,nvr,opts,P0,x0,nvr0,tab)
%                        1 2 3  4   5   6   7  8   9
%
% y: Time series (*)
% z: Regressors (*)
% x: States on which the parameters depend, column for each regressor (z)
% TVP: Model type for each TVP (0-RW, 1-IRW) (0)
% nvr: NVR hyper-parameters (0)
%        if an nvr is a negative integer (-n) it is
%        optimised over the first n backfit steps
% opts: estimation options [iter con meth sm ALG plotopt], use -1 for defaults
%         iter: number of backfitting iterations (10)
%         con: backfitting convergence threshold (0.001)
%         meth: Optimisation estimation method (0)
%                 0: Maximum Likelihood
%                 Integer:Sum of squares of the #-step-ahead forecasting errors
%         sm: Smoothing on (1-default) or off (0-saves memory)
%         ALG: Smoothing algorithm P (0) or Q (1-default)
%         plotopt: Plot results during estimation on (1) or off (0-default)
% P0: Initial P matrix diagonal (1e5)
% x0: Initial state vector (0)
% nvr0: Initial NVRs specified by user (0.0001)
% tab:  controls progress display, 0 - no display, 2 - continuous display
%
% fit: Model fit
% fitse: Standard error of the fit
% par: Parameter estimates
% parse: Standard errors of parameters
% xs: Sorted states
% pars: Sorted parameter estimates
% parses: Sorted standard errors of parameters
% rsq: R squared value
% nvre: Estimated NVRs
% y0: Interpolated data
%
% Example: sdr(y,[u1 u2],[x1 x2],[0 1],[-2 -1])
%   regression type model y = c1(x1)*u1 + c2(x1)*u2, with an RW model for
%   c1 where the dependent state is x1 and an IRW model for c2 where the
%   dependent state is x2; the NVR for the first SDP (c1) is optimised at
%   the first iteration and at the first two iterations for the second SDP
%
% See also FCAST, STAND, DLR, DHR
%

% Copyright (c) 2004 by CRES, Lancaster University, United Kingdom
% Authors : Peter Young, Wlodek Tych, Diego Pedregal, James Taylor
% Additional author: Paul McKenna

% Captain Tbx. functions used: dlr01, ksirw2c, dlropt01, dlrfun01

% Revision history :
%   28/01/2005, WT, conversion into SDR solver, code cosmetics
%   24/09/2003, PM, corrected for x0 input with IRW
%   23/09/2003, PM, corrected for x0 input
%   05/03/2003, JT, interpolated data output argument
%   22/10/2002, PM, adapted for convergence with missing values in output
%   10/07/2002, JT, default NVR changed to 1 in line nvr(j)=min(nvr2,1)
%   30/05/2002, JT, included in toolbox, P0 now passed to DLROPT as a vector
%   18/04/2002, PM, nvr input made more flexible for optimisation
%   17/04/2002, PM, optimisation set to ALG=0 to avoid problems with for Matlab 4
%   17/04/2002, PM, nvr input padded with end values if too short
%   15/04/2002, PM, opts input argument created to make the function call more compact

% Argument sizes
Ndat=size(z,1);nr=size(z,2);

% ******* Default argument settings
% NVR0
if nargin<10,tab=[]; end
if nargin<9;nvr0=[];end
if isempty(nvr0);nvr0=0.0001;end
if length(nvr0)==1 & nr~=1;nvr0=nvr0*ones(1,nr);end
lnvr0=length(nvr0);nvr0=nvr0(:)';nvr0=[nvr0 ones(1,nr-lnvr0)*nvr0(lnvr0)];
% x0: initial conditions (state)
if nargin<8,x0=[];end
% P0: initial conditions on P-matrix
if nargin<7,P0=1e5;end
if isempty(P0);P0=1e5;end
% options
if nargin<6,opts=[];end
if isempty(opts),opts=[10 0.001 0 1 1 0];end
opts0=[10 0.001 0 1 1 0];
if length(opts)~=6;
    opts=[opts opts0(length(opts)+1:6)];
end
if isempty(tab), tab=2; end

% ******* Options setting
opts(find(opts<0))=opts0(find(opts<0));
bf_it=opts(1);
con=opts(2);
sm=opts(4);
ALG=opts(5);
plotopt=opts(6);
if opts(3)==0;
    meth='ml';
else
    meth=['f' int2str(opts(3))];
end

%  ********
% Default NVR settings
if nargin<5,nvr=0;end
if isempty(nvr);nvr=0;end
if length(nvr)==1 & nr~=1;nvr=nvr*ones(1,nr);end
lnvr=length(nvr);nvr=nvr(:)';nvr=[nvr ones(1,nr-lnvr)*nvr(lnvr)];

% Default TVP settings (type of RW model)
if nargin<4,TVP=0;end
if isempty(TVP);TVP=0;end
if length(TVP)==1 & nr~=1;TVP=TVP*ones(1,nr);end
lTVP=length(TVP);TVP=TVP(:)';TVP=[TVP ones(1,nr-lTVP)*TVP(lTVP)];
if ~isempty(P0),lP0=length(P0);P0=P0(:)';P0=[P0 ones(1,sum(TVP)+nr-lP0)*P0(lP0)];end

% Default settings for state and regressors
if nargin <3, x=z;end
if isempty(x);x=z;end

%marking the state locations that relate to nvrs (not alpha)
ipp=[];
for i=1:length(TVP);ipp=[ipp 1 zeros(1,TVP(i))];end
ipp=find(ipp);

optim=zeros(size(nvr));
par=zeros(Ndat,nr);
parse=zeros(Ndat,nr);fit=zeros(Ndat,1);
ccv=zeros(bf_it,1); % convergence criteria vector
for i=1:length(nvr);
    if nvr(i)<0;
        optim(i)=-nvr(i);
    end
end

% Estimate time invariant parameter estimates
par0=(z\y)';
par=par0(ones(length(y),1),:);

% Optional graphics output
if plotopt;
    clf;
    for i=1:nr;
        subplot(nr,1,i);
        %plot(z(:,i),par(:,i),'b-');drawnow;hold on
        ll=plot(x(:,i),par(:,i),'b.');set(ll,'markersize',1);drawnow;hold on
    end
end

% Optimisation options - variable recycling
opts=foptions;
opts(2)=0.01;%termination tolerance
i=0;
for i=1:bf_it;
    for j=1:nr;
        par0=par;
        % Generate partial residual
        npar=[1:nr];npar=[npar(npar<j) npar(npar>j)];
        if nr==1;
            yc=y;
        else
            if nr>2;
                yc=y-sum(par0(:,npar)'.*z(:,npar)')';
            else
                yc=y-(par0(:,npar).*z(:,npar));
            end
        end
        % Sort data
        [ys,I]=sort(x(:,j));		%sort w.r.t. j-th regressor
        ysc=yc(I);			    %partial residual, sorted and corrected
        xr=z(:,j);			    %regressor
        xs=xr(I);			        %regressor, sorted
        %optimise
        if i<=optim(j);
            if i==1;
                [nvr2,alpha2,opts2,parse2]=dlropt01(ysc,xs,TVP(j),meth,-2,1,nvr0(j),1,opts,0,tab,P0(ipp(j):ipp(j)+TVP(j)));% ALG=0 for Matlab4
            else
                [nvr2,alpha2,opts2,parse2]=dlropt01(ysc,xs,TVP(j),meth,-2,1,nvr(j),1,opts,0,tab,P0(ipp(j):ipp(j)+TVP(j)));% ALG=0 for Matlab4
            end
            nvr(j)=min(nvr2,1);  %  default=1
        end

        % Estimate
        if i==1
            if isempty(x0);
                [fit,fitse,par1,parse1,comp,e]=dlr01(ysc,xs,TVP(j),nvr(j),1,P0(ipp(j):ipp(j)+TVP(j)),x0,1);
            else
                [fit,fitse,par1,parse1,comp,e]=dlr01(ysc,xs,TVP(j),nvr(j),1,P0(ipp(j):ipp(j)+TVP(j)),[x0(j);zeros(TVP(j),1)],1);
            end
        else
            [fit,fitse,par1,parse1,comp,e]=dlr01(ysc,xs,TVP(j),nvr(j),1,P0(ipp(j):ipp(j)+TVP(j)),[par(I(1),j);zeros(TVP(j),1)],1);
        end

        par(I,j)=par1;
        parse(I,j)=parse1;
        if plotopt;
            subplot(nr,1,j);
            if i>1;eval(['set(l' int2str(j) ',''color'',[0.5 0.5 0.5]);']);end
            eval(['l' int2str(j) '=plot(x(:,j),par(:,j),''r.'');']);eval(['set(l' int2str(j) ',''markersize'',1);']);drawnow
        end

    end

    if i<=max(optim);
        disp(['Optimised NVRs (after ' int2str(i) ' of ' int2str(bf_it) ' backfitting iterations)'])
        disp(nvr)
        if all(optim<=i) & j==nr;
            disp('Retaining these values for subsequent iterations')
        end
    end

    NotMiss=find(~isnan(y) & ~sum(isnan(par)')'  & ~sum(isnan(z)')' );% Flag "not missing" values

    if nr>1;
        disp([i 1-(cov(y(NotMiss)-sum(par(NotMiss,:)'.*z(NotMiss,:)')')/cov(y(NotMiss)))])
    else
        disp([i 1-(cov(y(NotMiss)-(par(NotMiss,:).*z(NotMiss,:)))/cov(y(NotMiss)))])
    end
    if nr>1;
        ccv(i)=1-(cov(y(NotMiss)-sum(par(NotMiss,:)'.*z(NotMiss,:)')')/cov(y(NotMiss)));
    else
        ccv(i)=1-(cov(y(NotMiss)-(par(NotMiss,:).*z(NotMiss,:)))/cov(y(NotMiss)));
    end
    if i>3;
        if abs(ccv(i)-ccv(i-1))<con | abs(ccv(i)-ccv(i-2))<con;
            disp(['Converged on backfit iteration number ' int2str(i)]);break
        else
            if i==bf_it;
                disp('Maximum number of backfit iterations reached, convergence criterion not achieved')
            end
        end
    end
end
[zs,zI]=sort(x);
pars=zeros(size(par));
parses=zeros(size(parse));
for i=1:nr;
    pars(:,i)=par(zI(:,i),i);
    parses(:,i)=parse(zI(:,i),i);
end
if nr>1;
    fit=sum(par'.*z')';
else
    fit=par.*z;
end
nvre=nvr;
rsq=1-(cov(y(NotMiss)-fit(NotMiss))/cov(y(NotMiss)));

% interpolated data
ii=find(isnan(y));
y0=y;
y0(ii)=fit(ii);

% end of m-file