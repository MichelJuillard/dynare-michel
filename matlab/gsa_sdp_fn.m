function gsa_ = gsa_sdp_fn(y, x, pshape, p1, p2, tvp, fname, pnames, ifig, ipred, nvr, ialias, nvrT, yi_anal)

% copyright Marco Ratto 2006
%
% function gsa_ =
% gsa_sdp_fn(y, x, pshape, p1, p2, tvp, k, pnames, ifig, ipred, nvr, ialias, nvrT, yi_anal)
%
% INPUTS
% y, the output
% x, the inputs
% pshape, input distributions
% p1, p2, parameters of input distributions
% tvp: tipe of estimation
%      0 = load previous estimation
%      1 = RW for all parameters
%      2 = IRW for all parameters
%      -1 = RW for all parameters, fixed nvr's (no estimation)
%      -2 = IRW for all parameters, fixed nvr's (no estimation)
%      an array of length(x) specifies the process for each parameter
% fname, file name to save the analysis
% pnames: names of the params
% ifig: 1 plot figures, 0 no plots
% ipred: use polynomial interpolation to predict model output
% nvr
% yi_anal, analytical value of the first order HDMR terms (optional)
%
% OUTPUT
%   gsa_.univariate.f
%   gsa_.univariate.fs
%   gsa_.univariate.fses
%   gsa_.univariate.si
%   gsa_.univariate.si_std
%   gsa_.multivariate.f
%   gsa_.multivariate.fs
%   gsa_.multivariate.fses
%   gsa_.multivariate.si
%   gsa_.multivariate.si_std
%   gsa_.multivariate.stat
%   gsa_.x0
%   gsa_.xx
%   gsa_.y
% 
% f, function estimates
% fs, sorted function estimates
% fses, sorted standard error of function estimates
% si, sensitivity indices
% si_std, st. error of sensitivity indices
% stat, euristic tstat for significance of f
% xx, transformed inputs used (rank-transformed)
% x0, original inputs
% y, output


[kn, np]=size(x);
if nargin<6 | isempty(tvp),
  tvp=-2;
  TVP=ones(1,np).*(abs(tvp)-1);
elseif length(tvp)==1,
  if tvp>0
    TVP=ones(1,np).*(tvp-1);
  else
    TVP=ones(1,np).*(abs(tvp)-1);
  end
else
  TVP=tvp-1;    
end
if nargin<7 | isempty(fname),
  fname='gsa_sdp';
end
if nargin<9 | isempty(ifig),
  ifig=0;
end
if nargin<10 | isempty(ipred),
  ipred=0;
end
if nargin<12 | isempty(ialias),
  ialias=0;
end
if nargin<13 | isempty(nvrT),
  nvrT=16;
end
if nargin<14 | isempty(yi_anal),
  yi_anal=[];
end

yin=y;
[y,my,sy]=stand_(y);
x0=x;
% x=stand_(x);
if tvp~=0,
  [xcum, paracdf] = cumdens(x,pshape,p1,p2);
  [dum, ix]=sort(x);
  [dum, iy]=sort(ix);
  xx=kron([1:kn]',ones(1,np));
  xx=xx(iy)./kn;
  x=xx;
  x=x+1;
  nvr1 = ( 2*sin( 0.5*(2*pi/(kn/2)) ) )^2;
  nvr2 = ( 2*sin( 0.5*(2*pi/(kn/2)) ) )^4;
  %nvr2 = ( 2*sin( 0.5*(2*pi/(3*(kn/4))) ) )^4;
  
  
  opts=[40 0.001 0 1 1 0]; %[iter con meth sm ALG plotopt]
  for j=1:np,
    
    lastwarn('')
    if tvp<0
      if tvp==-1
        [fit(:,j),fitse(:,j),par(:,j),parse(:,j),zs(:,j),pars(:,j),parses(:,j),rsq(j),nvre(j)]= ...
          sdr(y,x(:,j),x(:,j), 0, nvr1,[1 -1 -1 -1 -1 -1]);
      else
        [fit(:,j),fitse(:,j),par(:,j),parse(:,j),zs(:,j),pars(:,j),parses(:,j),rsq(j),nvre(j)]= ...
          sdr(y,x(:,j),x(:,j), 1, nvr2,[1 -1 -1 -1 -1 -1]);
      end
    else
      [fit(:,j),fitse(:,j),par(:,j),parse(:,j),zs(:,j),pars(:,j),parses(:,j),rsq(j),nvre(j)]= ...
        sdr(y,x(:,j),x(:,j), TVP(j), -1,[1 -1 -1 -1 -1 -1],[],[],[],0);
    end
    dwarn=lastwarn;
    if length(dwarn)>0,
      nvre(j)=0;
    end  
    m0=mean(zs(:,j).*pars(:,j));
    vi(1, j) = var(zs(:,j).*pars(:,j));
    vi(2, j) = max(abs([(par(:,j)-mean(par(:,j)))./parse(:,j)])) ;         
    rT(1, j) = rsq(j);
    fu(:,j) =(x(:,j).*par(:,j)-m0).*sy;
    fus(:,j)=zs(:,j).*pars(:,j).*sy; 
    fuses(:,j)=zs(:,j).*parses(:,j).*sy;
    lmod(j) = {'RW'};
    if TVP(j), lmod(j)={'IRW'}; end,
  end
  
  %imult = find(rT(1,:)>0.005);
  [dum, imult]=sort(-vi(2,:));
  %imult=[11:15];
  %imult=1:np;
  if isempty(imult), 
    mrT=max(rT(1,:));
    imult = find(rT(1,:)>(0.05*mrT));
  end
  %[fit,fitse,par,parse,zs,pars,parses,rsq,nvre]=sdr(y,x,x,ones(1,np), ones(1,np).*(-4));
  opts=[40 0.001 0 1 1 0]; %[iter con meth sm ALG plotopt]
  %nvrP=ones(1,np).*(-2).*(nvre~=0);
  y0=y;
  if tvp<0,
    if tvp==-1
      nvrP=ones(1,np).*nvr1;
    else
      nvrP=ones(1,np).*nvr2;
    end
  else
    nvrP=ones(1,np).*(-2);
    if ialias,
      nvr0=nvre.*(nvre~=0)+ones(1,np).*(1.e-10).*(nvre==0);    
      nvrP(ialias)=nvre(ialias);
      nvr0(ialias)=nvre(ialias);
      y0 = y0 - par(:,ialias).*x(:,ialias);
      imult=imult(find(imult~=ialias));
    else
      nvr0=nvre.*(nvre~=0)+ones(1,np).*(1.e-10).*(nvre==0);    
    end
  end
  parM=par;
  parseM=parse;
  zsM=zs;
  parsM=pars;
  parsesM=parses;
  nvreM=nvre;
  if tvp<0
    [fit2,fitse2,par2,parse2,zs2,pars2,parses2,rsq2,nvre2]= ...
      sdr(y0,x(:,imult),x(:,imult),abs(tvp)-1,nvrP(imult));
  else
    [fit2,fitse2,par2,parse2,zs2,pars2,parses2,rsq2,nvre2]= ...
      sdr(y0,x(:,imult),x(:,imult),TVP(imult),nvrP(imult),opts,[],[],nvr0(imult),0);
  end
  if any(nvre2>nvrT),
    ik=find(nvre2>nvrT);
    nvre2(ik)=ones(1,length(ik)).*nvrT;
    [fit2,fitse2,par2,parse2,zs2,pars2,parses2,rsq2,nvre2]= ...
      sdr(y0,x(:,imult),x(:,imult),TVP(imult),nvre2);
  end
  parM(:,imult)=par2;
  parseM(:,imult)=parse2;
  zsM(:,imult)=zs2;
  parsM(:,imult)=pars2;
  parsesM(:,imult)=parses2;
  nvreM(imult)=nvre2;
  for j=1:np,
    m0=mean(zsM(:,j).*parsM(:,j));
    v_di_e(1, j) = var(zsM(:,j).*parsM(:,j));
    v_di_e(2,j) = mean((zsM(:,j).*parsesM(:,j)).^2);
    stat(1, j) = max(abs([(parM(:,j)-mean(parM(:,j)))./parseM(:,j)])) ;         
    f(:,j)=(x(:,j).*parM(:,j)-m0).*sy;
    fs(:,j)=(zsM(:,j).*parsM(:,j)-m0).*sy;
    fses(:,j)=zsM(:,j).*parsesM(:,j).*sy;
  end
  pcoef = ChRegressGSA(xcum,8,1,f);
  for ij=1:np,
    scoef{ij}=['[',num2str(pcoef(ij,:),'%20.16e '),']'];
    gg{ij}=inline(['polyval(',scoef{ij},',',char(paracdf{ij}),')'],'X');
  end
  gsa_.univariate.f0=my;
  gsa_.univariate.f=fu;
  gsa_.univariate.fs=fus;
  gsa_.univariate.fses=fuses;
  gsa_.univariate.nvr=nvre;
  gsa_.univariate.r2=rsq2;
  gsa_.univariate.si=vi(1,:);
  gsa_.univariate.si_std=vi(2,:);
  gsa_.multivariate.f0=my;
  gsa_.multivariate.f=f;
  gsa_.multivariate.fs=fs;
  gsa_.multivariate.fses=fses;
  gsa_.multivariate.nvr=nvreM;
  gsa_.multivariate.r2=rsq2;
  gsa_.multivariate.si=v_di_e(1,:);
  gsa_.multivariate.si_std=v_di_e(2,:);
  gsa_.multivariate.stat=stat;
  gsa_.polynomial.gpred=gg;
  gsa_.x0=x0;
  gsa_.xx=xx;
  gsa_.y=yin;
  save(fname,'gsa_')
else
  load(fname)
%   xx=gsa_.xx;
%   yin=gsa_.y;
%   [xcum, paracdf] = cumdens(x,pshape,p1,p2);
%   pcoef = ChRegressGSA(xcum,8,1,gsa_.multivariate.f);
%   for ij=1:np,
%     scoef{ij}=['[',num2str(pcoef(ij,:),'%20.16e '),']'];
%     gg{ij}=inline(['polyval(',scoef{ij},',',char(paracdf{ij}),')'],'X');
%   end
%   gsa_.univariate.f0=my;
%   gsa_.multivariate.f0=my;
%   gsa_.multivariate.gpred=gg;
%   save(fname,'gsa_')
end  

if ifig
  nplo = ceil(sqrt(np));
  if nplo*(nplo-2)>=np,
    nplo2=nplo-2;
  elseif nplo*(nplo-1)>=np,
    nplo2=nplo-1;    
  elseif (nplo+1)*(nplo-1)>=np,
    nplo2=nplo-1;    
    nplo=nplo+1;    
  else
    nplo2=nplo;
  end
  
  if nargin<8 | isempty(pnames),
    pnames=['X_{',num2str(1),'}'];
    for j=2:np,
      pnames=str2mat(pnames,['X_{',num2str(j),'}']);
    end
    
  end
  
  figure('name',['Nrun=',num2str(kn),', ',fname,', SDP multivariate']),
  f=gsa_.multivariate.f;
  fs=gsa_.multivariate.fs;
  fses=gsa_.multivariate.fses;
  stat=gsa_.multivariate.stat;
  v_di_e(1,:)=gsa_.multivariate.si;
  for j=1:np,
    %subplot(3,4,j)
    subplot(nplo2,nplo,j)
    xp=sort(gsa_.x0(:,j));
    %xp=zsM(:,jj)-1;
    if exist('yi_anal') & ~isempty(yi_anal),
      %plot(xp, yi_anal(:,j)./sy,'r:'), hold on, 
      plot(xp, yi_anal(:,j),'r:'), hold on, 
    end
    plot(xp, fs(:,j),'k', xp, 3.09.*fses(:,j),':k', xp, -3.09.*fses(:,j),':k')
    set(gca,'xlim',[min(xp) max(xp)])        
    lmod(j) = {'I1'};
    if TVP(j), lmod(j)={'I2'}; end,
    title([deblank(pnames(j,:)),', Si=',num2str(v_di_e(1,j),'%5.2f')],'fontsize',9,'interpreter','none')
  end
  saveas(gcf,[fname,'_SE_multiv'])
  eval(['print -depsc2 ' fname '_SE_multiv']);
  eval(['print -dpdf ' fname '_SA_multiv']);
  close(gcf)
  
%   figure('name',['Nrun=',num2str(kn),', ',fname,', SDP univariate']),    
%   f=gsa_.univariate.f;
%   fs=gsa_.univariate.fs;
%   fses=gsa_.univariate.fses;
%   vi(1,:)=gsa_.univariate.si;
%   vi(2,:)=gsa_.univariate.si_std;
%   for j=1:np,
%     subplot(nplo2,nplo,j)
%     xp=sort(x0(:,j));
%     if exist('yi_anal') & ~isempty(yi_anal),
%       plot(xp, yi_anal(:,j),'r'), hold on, 
%     end
%     plot(xp, [fs(:,j) 3.09.*fses(:,j) -3.09.*fses(:,j)])
%     set(gca,'xlim',[min(xp) max(xp)])        
%     lmod(j) = {'I1'};
%     if TVP(j), lmod(j)={'I2'}; end,
%     title([deblank(pnames(j,:)),', ',lmod{j},', tstat=',num2str(vi(2,j),'%5.2f'), ...
%         ', S=',num2str(vi(1,j),'%5.2f')])
%   end
%   %saveas(gcf,[fname,'_SE_univ'])
%   close(gcf)
end
%end

if ipred,
  gg=gsa_.polynomial.gpred;
  for j=1:np; 
    ff(:,j)=gg{j}(x0(:,j)); 
  end
  gsa_.polynomial.ypred=sum(ff,2)+gsa_.multivariate.f0;  
  gsa_.polynomial.r2=1-cov(gsa_.polynomial.ypred-yin)/cov(yin);
  for j=1:12; subplot(3,4,j), plot(x0(:,j),ff(:,j),'.',gsa_.x0(:,j),f(:,j),'b.'), end
  gsa_.xpred=x0;  
%   nvreM=gsa_.multivariate.nvr;
%   x=[gsa_.x0; x0];
%   xcum = cumdens(x,pshape,p1,p2);
%   [yin,my,sy]=stand_(gsa_.y);
%   yf=[yin; NaN.*ones(size(x0,1),1)];    
%   yfit2= sdp(yf,xcum,xcum,TVP,nvreM);
%   gsa_.ypred_sdr = yfit2(length(yin)+1:end).*sy+my;  
  
  save(fname,'gsa_')
end