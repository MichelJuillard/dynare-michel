function redform_map(dirname)
%function redform_map(dirname)
% inputs (from opt_gsa structure
% anamendo    = options_gsa_.namendo;
% anamlagendo = options_gsa_.namlagendo;
% anamexo     = options_gsa_.namexo;
% iload       = options_gsa_.load_redform;
% pprior      = options_gsa_.pprior;
% ilog        = options_gsa_.logtrans_redform;
% threshold   = options_gsa_.threshold_redform;
% ksstat      = options_gsa_.ksstat_redform;
% alpha2      = options_gsa_.alpha2_redform;
%
% Part of the Sensitivity Analysis Toolbox for DYNARE
%
% Written by Marco Ratto, 2006
% Joint Research Centre, The European Commission,
% (http://eemc.jrc.ec.europa.eu/),
% marco.ratto@jrc.it 
%
% Disclaimer: This software is not subject to copyright protection and is in the public domain. 
% It is an experimental system. The Joint Research Centre of European Commission 
% assumes no responsibility whatsoever for its use by other parties
% and makes no guarantees, expressed or implied, about its quality, reliability, or any other
% characteristic. We would appreciate acknowledgement if the software is used.
% Reference:
% M. Ratto, Global Sensitivity Analysis for Macroeconomic models, MIMEO, 2006.
%

global M_ oo_ estim_params_ options_ bayestopt_

options_gsa_ = options_.opt_gsa;

anamendo = options_gsa_.namendo;
anamlagendo = options_gsa_.namlagendo;
anamexo = options_gsa_.namexo;
iload = options_gsa_.load_redform;
pprior = options_gsa_.pprior;
ilog = options_gsa_.logtrans_redform;
threshold = options_gsa_.threshold_redform;
ksstat = options_gsa_.ksstat_redform;
alpha2 = options_gsa_.alpha2_redform;

pnames = M_.param_names(estim_params_.param_vals(:,1),:);
if nargin==0,
  dirname='';
end

if pprior
  load([dirname,'\',M_.fname,'_prior']);
  adir=[dirname '\redform_stab'];
else
  load([dirname,'\',M_.fname,'_mc']);
  adir=[dirname '\redform_mc'];
end
if isempty(dir(adir))
  mkdir(adir)
end
adir0=pwd;
%cd(adir)

nspred=size(T,2)-M_.exo_nbr;
x0=lpmat(istable,:);
[kn, np]=size(x0);
offset = length(bayestopt_.pshape)-np;
pshape = bayestopt_.pshape(offset+1:end);
pd =  [bayestopt_.p1(offset+1:end) bayestopt_.p2(offset+1:end) bayestopt_.p3(offset+1:end) bayestopt_.p4(offset+1:end)];

nsok = length(find(M_.lead_lag_incidence(M_.maximum_lag,:)));
clear lpmat lpmat0 egg iunstable yys
js=0;
for j=1:size(anamendo,1)
  namendo=deblank(anamendo(j,:));
  iendo=strmatch(namendo,M_.endo_names(oo_.dr.order_var,:),'exact');
  ifig=0;
  iplo=0;
  for jx=1:size(anamexo,1)
    namexo=deblank(anamexo(jx,:));
    iexo=strmatch(namexo,M_.exo_names,'exact');

    if ~isempty(iexo),
      %y0=squeeze(T(iendo,iexo+nspred,istable));
      y0=squeeze(T(iendo,iexo+nspred,:));
      if (max(y0)-min(y0))>1.e-10,
        if mod(iplo,9)==0,
          ifig=ifig+1;
          figure('name',[namendo,' vs. shocks ',int2str(ifig)]),
          iplo=0;
        end
        iplo=iplo+1;
        js=js+1;
        xdir0 = [adir,'\',namendo,'_vs_', namexo];
        if ilog==0,
          if isempty(threshold)
            si(:,js) = redform_private(x0, y0, pshape, pd, iload, pnames, namendo, namexo, xdir0);
          else
            iy=find( (y0>threshold(1)) & (y0<threshold(2)));
            iyc=find( (y0<=threshold(1)) | (y0>=threshold(2)));
            xdir = [xdir0,'_cut'];
            si(:,js) = redform_private(x0(iy,:), y0(iy), pshape, pd, iload, pnames, namendo, namexo, xdir);
            delete([xdir, '\*cut*.*'])
            [proba, dproba] = stab_map_1(x0, iy, iyc, 'cut',0);
            indsmirnov = find(dproba>ksstat);
            stab_map_1(x0, iy, iyc, 'cut',1,indsmirnov,xdir);
            stab_map_2(x0(iy,:),alpha2,'cut',xdir)
            stab_map_2(x0(iyc,:),alpha2,'trim',xdir)
          end
        else
          [yy, xdir] = log_trans_(y0,xdir0);
          silog(:,js) = redform_private(x0, yy, pshape, pd, iload, pnames, namendo, namexo, xdir);
        end

        subplot(3,3,iplo),
        if ilog,
          [saso, iso] = sort(-silog(:,js));
          bar([silog(iso(1:10),js)])
          logflag='log';
        else
          [saso, iso] = sort(-si(:,js));
          bar(si(iso(1:10),js))
          logflag='';
        end
        %set(gca,'xticklabel',pnames(iso(1:10),:),'fontsize',8)
        set(gca,'xticklabel',' ','fontsize',10)
        set(gca,'xlim',[0.5 10.5])
        for ip=1:10,
          text(ip,-0.02,deblank(pnames(iso(ip),:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
        end
        title([logflag,' ',namendo,' vs. ',namexo],'interpreter','none')
        if iplo==9,
          saveas(gcf,[dirname,'\',M_.fname,'_redform_', namendo,'_vs_shocks_',logflag,num2str(ifig)])
          eval(['print -depsc2 ' dirname,'\',M_.fname,'_redform_', namendo,'_vs_shocks_',logflag,num2str(ifig)]);
          eval(['print -dpdf ' dirname,'\',M_.fname,'_redform_', namendo,'_vs_shocks_',logflag,num2str(ifig)]);
          close(gcf)
        end
      
      end
    end
  end
  if iplo<9 & iplo>0 & ifig,
    saveas(gcf,[dirname,'\',M_.fname,'_redform_', namendo,'_vs_shocks_',logflag,num2str(ifig)])
    eval(['print -depsc2 ' dirname,'\',M_.fname,'_redform_', namendo,'_vs_shocks_',logflag,num2str(ifig)]);
    eval(['print -dpdf ' dirname,'\',M_.fname,'_redform_', namendo,'_vs_shocks_',logflag,num2str(ifig)]);
    close(gcf)
  end
  ifig=0;
  iplo=0;
  for je=1:size(anamlagendo,1)
    namlagendo=deblank(anamlagendo(je,:));
    ilagendo=strmatch(namlagendo,M_.endo_names(oo_.dr.order_var(oo_.dr.nstatic+1:oo_.dr.nstatic+nsok),:),'exact');

    if ~isempty(ilagendo),
      %y0=squeeze(T(iendo,ilagendo,istable));
      y0=squeeze(T(iendo,ilagendo,:));
      if (max(y0)-min(y0))>1.e-10,
        if mod(iplo,9)==0,
          ifig=ifig+1;
          figure('name',[namendo,' vs. lags ',int2str(ifig)]),
          iplo=0;
        end
        iplo=iplo+1;
        js=js+1;
        xdir0 = [adir,'\',namendo,'_vs_', namlagendo];
        if ilog==0,
        if isempty(threshold)
          si(:,js) = redform_private(x0, y0, pshape, pd, iload, pnames, namendo, namlagendo, xdir0);
        else
          iy=find( (y0>threshold(1)) & (y0<threshold(2)));
          iyc=find( (y0<=threshold(1)) | (y0>=threshold(2)));
          xdir = [xdir0,'_cut'];
          si(:,js) = redform_private(x0(iy,:), y0(iy), pshape, pd, iload, pnames, namendo, namlagendo, xdir);
          delete([xdir, '\*cut*.*'])
          [proba, dproba] = stab_map_1(x0, iy, iyc, 'cut',0);
          indsmirnov = find(dproba>ksstat);
          stab_map_1(x0, iy, iyc, 'cut',1,indsmirnov,xdir);
          stab_map_2(x0(iy,:),alpha2,'cut',xdir)
          stab_map_2(x0(iyc,:),alpha2,'trim',xdir)
        end
        else
          [yy, xdir] = log_trans_(y0,xdir0);
          silog(:,js) = redform_private(x0, yy, pshape, pd, iload, pnames, namendo, namlagendo, xdir);
        end

        subplot(3,3,iplo),
        if ilog,
          [saso, iso] = sort(-silog(:,js));
          bar([silog(iso(1:10),js)])
          logflag='log';
        else
          [saso, iso] = sort(-si(:,js));
          bar(si(iso(1:10),js))
          logflag='';
        end
        %set(gca,'xticklabel',pnames(iso(1:10),:),'fontsize',8)
        set(gca,'xticklabel',' ','fontsize',10)
        set(gca,'xlim',[0.5 10.5])
        for ip=1:10,
          text(ip,-0.02,deblank(pnames(iso(ip),:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
        end
        title([logflag,' ',namendo,' vs. ',namlagendo,'(-1)'],'interpreter','none')
        if iplo==9,
          saveas(gcf,[dirname,'\',M_.fname,'_redform_', namendo,'_vs_lags_',logflag,num2str(ifig)])
          eval(['print -depsc2 ' dirname,'\',M_.fname,'_redform_', namendo,'_vs_lags_',logflag,num2str(ifig)]);
          eval(['print -dpdf ' dirname,'\',M_.fname,'_redform_', namendo,'_vs_lags_',logflag,num2str(ifig)]);
          close(gcf)
        end
      
      end
    end
  end
  if iplo<9 & iplo>0 & ifig,
    saveas(gcf,[dirname,'\',M_.fname,'_redform_', namendo,'_vs_lags_',logflag,num2str(ifig)])
    eval(['print -depsc2 ' dirname,'\',M_.fname,'_redform_', namendo,'_vs_lags_',logflag,num2str(ifig)]);
    eval(['print -dpdf ' dirname,'\',M_.fname,'_redform_', namendo,'_vs_lags_',logflag,num2str(ifig)]);
    close(gcf)
  end
end

if ilog==0,
figure, bar(si)
set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
set(gca,'xlim',[0.5 np+0.5])
set(gca,'position',[0.13 0.2 0.775 0.7])
for ip=1:np,
  text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Reduced form GSA')

saveas(gcf,[dirname,'\',M_.fname,'_redform_gsa'])
eval(['print -depsc2 ' dirname,'\',M_.fname,'_redform_gsa']);
eval(['print -dpdf ' dirname,'\',M_.fname,'_redform_gsa']);

else
figure, bar(silog)
set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
set(gca,'xlim',[0.5 np+0.5])
set(gca,'position',[0.13 0.2 0.775 0.7])
for ip=1:np,
  text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
title('Reduced form GSA - Log-transformed elements')

saveas(gcf,[dirname,'\',M_.fname,'_redform_gsa_log'])
eval(['print -depsc2 ' dirname,'\',M_.fname,'_redform_gsa_log']);
eval(['print -dpdf ' dirname,'\',M_.fname,'_redform_gsa_log']);

end
function si  = redform_private(x0, y0, pshape, pd, iload, pnames, namy, namx, xdir)
global bayestopt_

fname=[xdir,'\map'];
if iload==0,
  figure, hist(y0,30), title([namy,' vs. ', namx])
  if isempty(dir(xdir))
    mkdir(xdir)
  end
  saveas(gcf,[xdir,'\', namy,'_vs_', namx])
  eval(['print -depsc2 ' xdir,'\', namy,'_vs_', namx]);
  eval(['print -dpdf ' xdir,'\', namy,'_vs_', namx]);
  close(gcf)
  gsa_ = gsa_sdp_dyn(y0, x0, -2, [],[],[],1,fname, pnames);
else
  gsa_ = gsa_sdp_dyn(y0, x0, 0, [],[],[],0,fname, pnames);
end
si = gsa_.multivariate.si;

function [yy, xdir]=log_trans_(y0,xdir0)

f=inline('skewness(log(y+lam))','lam','y');
if ~(max(y0)<0 | min(y0)>0)
  isig=1;
  if skewness(y0)<0,
    isig=-1;
    y0=-y0;
  end
  n=hist(y0,10);
  if n(1)>20*n(end),
    try lam=fzero(f,[-min(y0)+10*eps -min(y0)+abs(median(y0))],[],y0);
    catch
      yl(1)=f(-min(y0)+10*eps,y0);
      yl(2)=f(-min(y0)+abs(median(y0)),y0);
      if abs(yl(1))<abs(yl(2))
        lam=-min(y0)+eps;
      else
        lam = -min(y0)+abs(median(y0)); %abs(100*(1+min(y0)));
      end
    end
    yy = log(y0+lam);
    xdir=[xdir0,'_logskew'];
  else
    yy = log(y0.^2);
    xdir=[xdir0,'_logsquared'];
  end
else
  if max(y0)<0
    y0=-y0;
    %yy=log(-y0);
    xdir=[xdir0,'_minuslog'];
  elseif min(y0)>0
    %yy=log(y0);
    xdir=[xdir0,'_log'];
  end
  try lam=fzero(f,[-min(y0)+10*eps -min(y0)+median(y0)],[],y0);
  catch
    yl(1)=f(-min(y0)+10*eps,y0);
      yl(2)=f(-min(y0)+abs(median(y0)),y0);
    if abs(yl(1))<abs(yl(2))
      lam=-min(y0)+eps;
    else
        lam = -min(y0)+abs(median(y0)); %abs(100*(1+min(y0)));
    end
  end
  yy = log(y0+lam);
end
