function redform_screen(dirname)
%function redform_map(dirname)
% inputs (from opt_gsa structure
% anamendo    = options_gsa_.namendo;
% anamlagendo = options_gsa_.namlagendo;
% anamexo     = options_gsa_.namexo;
% iload       = options_gsa_.load_redform;
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
nliv = options_gsa_.morris_nliv;

pnames = M_.param_names(estim_params_.param_vals(:,1),:);
if nargin==0,
  dirname='';
end

load([dirname,'\',M_.fname,'_prior'],'lpmat','lpmat0','istable','T');

nspred=oo_.dr.nspred;

[kn, np]=size(lpmat);
nshock = length(bayestopt_.pshape)-np;

nsok = length(find(M_.lead_lag_incidence(M_.maximum_lag,:)));

js=0;
for j=1:size(anamendo,1),
  namendo=deblank(anamendo(j,:));
  iendo=strmatch(namendo,M_.endo_names(oo_.dr.order_var,:),'exact');

  iplo=0;
  ifig=0;
  for jx=1:size(anamexo,1)
    namexo=deblank(anamexo(jx,:));
    iexo=strmatch(namexo,M_.exo_names,'exact');

    if ~isempty(iexo),
      y0=teff(T(iendo,iexo+nspred,:),kn,istable);
      if ~isempty(y0),
        if mod(iplo,9)==0,
          ifig=ifig+1;
          figure('name',[namendo,' vs. shocks ',int2str(ifig)]),
          iplo=0;
        end
        iplo=iplo+1;
        js=js+1;
        subplot(3,3,iplo),
        [SAmeas, SAMorris] = Morris_Measure_Groups(np+nshock, [lpmat0 lpmat], y0,nliv);
        SAM = squeeze(SAMorris(nshock+1:end,1));
        SA(:,js)=SAM./(max(SAM)+eps);
        [saso, iso] = sort(-SA(:,js));
        bar(SA(iso(1:min(np,10)),js))
        %set(gca,'xticklabel',pnames(iso(1:min(np,10)),:),'fontsize',8)
        set(gca,'xticklabel',' ','fontsize',10)
        set(gca,'xlim',[0.5 10.5])
        for ip=1:min(np,10),
          text(ip,-0.02,deblank(pnames(iso(ip),:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
        end
        title([namendo,' vs. ',namexo],'interpreter','none')
        if iplo==9,
          saveas(gcf,[dirname,'\',M_.fname,'_', namendo,'_vs_shock_',num2str(ifig)])
          eval(['print -depsc2 ' dirname,'\',M_.fname,'_', namendo,'_vs_shock_',num2str(ifig)]);
          eval(['print -dpdf ' dirname,'\',M_.fname,'_', namendo,'_vs_shock_',num2str(ifig)]);
          close(gcf)
        end

      end
    end
  end
  if iplo<9 & iplo>0 & ifig,
    saveas(gcf,[dirname,'\',M_.fname,'_', namendo,'_vs_shocks_',num2str(ifig)])
    eval(['print -depsc2 ' dirname,'\',M_.fname,'_', namendo,'_vs_shocks_',num2str(ifig)]);
    eval(['print -dpdf ' dirname,'\',M_.fname,'_', namendo,'_vs_shocks_',num2str(ifig)]);
    close(gcf)
  end

  iplo=0;
  ifig=0;
  for je=1:size(anamlagendo,1)
    namlagendo=deblank(anamlagendo(je,:));
    ilagendo=strmatch(namlagendo,M_.endo_names(oo_.dr.order_var(oo_.dr.nstatic+1:oo_.dr.nstatic+nsok),:),'exact');

    if ~isempty(ilagendo),
      y0=teff(T(iendo,ilagendo,:),kn,istable);
      if ~isempty(y0),
        if mod(iplo,9)==0,
          ifig=ifig+1;
          figure('name',[namendo,' vs. lagged endogenous ',int2str(ifig)]),
          iplo=0;
        end
        iplo=iplo+1;
        js=js+1;
        subplot(3,3,iplo),
        [SAmeas, SAMorris] = Morris_Measure_Groups(np+nshock, [lpmat0 lpmat], y0,nliv);
        SAM = squeeze(SAMorris(nshock+1:end,1));
        SA(:,js)=SAM./(max(SAM)+eps);
        [saso, iso] = sort(-SA(:,js));
        bar(SA(iso(1:min(np,10)),js))
        %set(gca,'xticklabel',pnames(iso(1:min(np,10)),:),'fontsize',8)
        set(gca,'xticklabel',' ','fontsize',10)
        set(gca,'xlim',[0.5 10.5])
        for ip=1:min(np,10),
          text(ip,-0.02,deblank(pnames(iso(ip),:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
        end

        title([namendo,' vs. ',namlagendo,'(-1)'],'interpreter','none')
        if iplo==9,
          saveas(gcf,[dirname,'\',M_.fname,'_', namendo,'_vs_lags_',num2str(ifig)])
          eval(['print -depsc2 ' dirname,'\',M_.fname,'_', namendo,'_vs_lags_',num2str(ifig)]);
          eval(['print -dpdf ' dirname,'\',M_.fname,'_', namendo,'_vs_lags_',num2str(ifig)]);
          close(gcf)
        end
      end
    end
  end
  if iplo<9 & iplo>0 & ifig,
    saveas(gcf,[dirname,'\',M_.fname,'_', namendo,'_vs_lags_',num2str(ifig)])
    eval(['print -depsc2 ' dirname,'\',M_.fname,'_', namendo,'_vs_lags_',num2str(ifig)]);
    eval(['print -dpdf ' dirname,'\',M_.fname,'_', namendo,'_vs_lags_',num2str(ifig)]);
    close(gcf)
  end
end

figure, 
%bar(SA)
% boxplot(SA','whis',10,'symbol','r.')
myboxplot(SA',[],'.',[],10)
set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
set(gca,'xlim',[0.5 np+0.5])
set(gca,'ylim',[0 1])
set(gca,'position',[0.13 0.2 0.775 0.7])
for ip=1:np,
  text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
end
xlabel(' ')
ylabel('Elementary Effects')
title('Reduced form screening')

saveas(gcf,[dirname,'\',M_.fname,'_redform_screen'])
eval(['print -depsc2 ' dirname,'\',M_.fname,'_redform_screen']);
eval(['print -dpdf ' dirname,'\',M_.fname,'_redform_screen']);

