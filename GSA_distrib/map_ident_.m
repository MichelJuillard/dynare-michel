function map_ident_(OutputDirectoryName)
global bayestopt_ M_ options_ estim_params_ oo_

opt_gsa = options_.opt_gsa;
fname_ = M_.fname;
nliv   = opt_gsa.morris_nliv;
ntra   = opt_gsa.morris_ntra;
itrans = opt_gsa.trans_ident;

np = estim_params_.np;
if opt_gsa.load_ident,
  gsa_flag=0;
else
  gsa_flag=-2;
end

pnames = M_.param_names(estim_params_.param_vals(:,1),:);

filetoload=[OutputDirectoryName '\' fname_ '_prior'];
load(filetoload,'lpmat','lpmat0','istable','T','nspred','nboth','nfwrd')
Nsam = size(lpmat,1);
nshock = size(lpmat0,2);
npT = np+nshock;

fname_ = M_.fname;

% th moments
[vdec, cc, ac] = mc_moments(T, lpmat0, oo_.dr);
if opt_gsa.morris==0,
  ifig=0;
  for j=1:M_.exo_nbr,
    if mod(j,6)==1
      figure('name',['Variance decomposition shocks']);
      ifig=ifig+1;
      iplo=0;
    end
    iplo=iplo+1;
    subplot(2,3,iplo)
    boxplot(squeeze(vdec(:,j,:))','whis',10,'symbol','.r')
    set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:size(options_.varobs,1)])
    set(gca,'xlim',[0.5 size(options_.varobs,1)+0.5])
    set(gca,'ylim',[-2 102])
    for ip=1:size(options_.varobs,1),
      text(ip,-4,deblank(options_.varobs(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
    end
    xlabel(' ')
    ylabel(' ')
    title(M_.exo_names(j,:),'interpreter','none')
    if mod(j,6)==0 | j==M_.exo_nbr,
      saveas(gcf,[OutputDirectoryName,'\',fname_,'_vdec_exo_',int2str(ifig)])
      eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_vdec_exo_',int2str(ifig)]);
      eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_vdec_exo_',int2str(ifig)]);
    end
  end
end
[vdec, j0, ir_vdec, ic_vdec] = teff(vdec,Nsam,istable);
[cc, j0, ir_cc, ic_cc] = teff(cc,Nsam,istable);
[ac, j0, ir_ac, ic_ac] = teff(ac,Nsam,istable);

[Aa,Bb] = ghx2transition(squeeze(T(:,:,1)), ...
  bayestopt_.restrict_var_list, ...
  bayestopt_.restrict_columns, ...
  bayestopt_.restrict_aux);
A = zeros(size(Aa,1),size(Aa,2)+size(Bb,2),length(istable));
A(:,:,1)=[Aa, Bb];
for j=2:length(istable),
  [Aa,Bb] = ghx2transition(squeeze(T(:,:,j)), ...
    bayestopt_.restrict_var_list, ...
    bayestopt_.restrict_columns, ...
    bayestopt_.restrict_aux);
  A(:,:,j)=[Aa, Bb];
end
clear T;

[nr,nc,nn]=size(A);
io=bayestopt_.mf1;
T1=A(io,1:nr,:);
ino=find(~ismember([1:nr],io));
T2=A(ino,1:nr,:);
R=A(:,nr+1:nc,:);
[tadj, iff] = speed(A(1:nr,1:nr,:),R,io,0.5);
[tadj, j0, ir_tadj, ic_tadj] = teff(tadj,Nsam,istable);
[iff, j0, ir_if, ic_if] = teff(iff,Nsam,istable);


[yt, j0]=teff(A,Nsam,istable);
[yt1, j01]=teff(T1,Nsam,istable);
[yt2, j02]=teff(T2,Nsam,istable);
[ytr, j0r]=teff(R,Nsam,istable);

yt=[yt1 yt2 ytr];

%   for j=1:nr,
%     for i=1:nc,
%       y0=squeeze(A(j,i,:));
%       if max(y0)-min(y0)>1.e-10,
%         j0=j0+1;
%         y1=ones(size(lpmat,1),1)*NaN;
%         y1(istable,1)=y0;
%         yt(:,j0)=y1;
%       end
%     end
%   end
%   yt = yt(:,j0);

if opt_gsa.morris==1,
  %OutputDir = CheckPath('GSA\SCREEN');
  SAMorris = [];
  for i=1:size(vdec,2),
    [SAmeas, SAMorris(:,:,i)] = Morris_Measure_Groups(npT, [lpmat0 lpmat], vdec(:,i),nliv);
  end
  SAvdec = squeeze(SAMorris(:,1,:))';
  figure,
  boxplot(SAvdec,'whis',10,'symbol','r.')
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
  set(gca,'xlim',[0.5 npT+0.5])
  ydum = get(gca,'ylim');
  set(gca,'ylim',[0 ydum(2)])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  for ip=1:npT,
    text(ip,-2,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  title('All variance decomposition')
  saveas(gcf,[OutputDirectoryName,'\',fname_,'_morris_vdec'])
  eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_morris_vdec']);
  eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_morris_vdec']);

  ifig = 0;
  for j=1:size(options_.varobs,1)
    if mod(j,6)==1
      figure('name',['EET variance decomposition observed variables']);
      ifig=ifig+1;
      iplo=0;
    end
    iplo=iplo+1;
    subplot(3,2,iplo)
    iv = find(ir_vdec==j);
    if ~isempty(iv)
      if length(iv)>1
        boxplot(SAvdec(iv,:),'whis',3,'symbol','r.');
      else
        plot(SAvdec(iv,:),'r.');
      end
      set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
      set(gca,'xlim',[0.5 npT+0.5])
      ydum = get(gca,'ylim');
      set(gca,'ylim',[0 ydum(2)])
      for ip=1:npT,
        text(ip,-2,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
      end
      xlabel(' ')
    end
    title(options_.varobs(j,:),'interpreter','none')
    if mod(j,6)==0 | j==size(options_.varobs,1)
      saveas(gcf,[OutputDirectoryName,'\',fname_,'_morris_vdec_varobs_',int2str(ifig)])
      eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_morris_vdec_varobs_',int2str(ifig)]);
      eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_morris_vdec_varobs_',int2str(ifig)]);
    end
  end

  ifig = 0;
  for j=1:M_.exo_nbr,
    if mod(j,6)==1
      figure('name',['EET variance decomposition shocks']);
      ifig=ifig+1;
      iplo=0;
    end
    iplo=iplo+1;
    subplot(3,2,iplo)
    iv = find(ic_vdec==j);
    if ~isempty(iv)
      if length(iv)>1
        boxplot(SAvdec(iv,:),'whis',3,'symbol','r.');
      else
        plot(SAvdec(iv,:),'r.');
      end
      set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
      set(gca,'xlim',[0.5 npT+0.5])
      ydum = get(gca,'ylim');
      set(gca,'ylim',[0 ydum(2)])
      for ip=1:npT,
        text(ip,-2,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
      end
      xlabel(' ')
    end
    title(M_.exo_names(j,:),'interpreter','none')
    if mod(j,6)==0 | j==M_.exo_nbr,
      saveas(gcf,[OutputDirectoryName,'\',fname_,'_morris_vdec_exo_',int2str(ifig)])
      eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_morris_vdec_exo_',int2str(ifig)]);
      eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_morris_vdec_exo_',int2str(ifig)]);
    end
  end


  SAMorris = [];
  for i=1:size(cc,2),
    [SAmeas, SAMorris(:,:,i)] = Morris_Measure_Groups(npT, [lpmat0 lpmat], cc(:,i),nliv);
  end
  SAcc = squeeze(SAMorris(:,1,:))';
  figure,
  boxplot(SAcc,'whis',10,'symbol','r.')
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
  set(gca,'xlim',[0.5 npT+0.5])
  ydum = get(gca,'ylim');
  set(gca,'ylim',[0 ydum(2)])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  for ip=1:npT,
    text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  title('EET All cross-correlation matrix')
  saveas(gcf,[OutputDirectoryName,'\',fname_,'_morris_cc'])
  eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_morris_cc']);
  eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_morris_cc']);

  ifig = 0;
  for j=1:size(options_.varobs,1)
    if mod(j,6)==1
      figure('name',['EET cross-correlations']);
      ifig=ifig+1;
      iplo=0;
    end
    iplo=iplo+1;
    subplot(3,2,iplo)
    iv = find(ir_cc==j);
    iv = [iv; find(ic_cc==j)];
    if ~isempty(iv)
      if length(iv)>1
        boxplot(SAcc(iv,:),'whis',3,'symbol','r.');
      else
        plot(SAcc(iv,:),'r.');
      end
      set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
      set(gca,'xlim',[0.5 npT+0.5])
      ydum = get(gca,'ylim');
      set(gca,'ylim',[0 ydum(2)])
      for ip=1:npT,
        text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
      end
      xlabel(' ')
    end
    title(options_.varobs(j,:),'interpreter','none')
    if mod(j,6)==0 | j==size(options_.varobs,1)
      saveas(gcf,[OutputDirectoryName,'\',fname_,'_morris_cc_',int2str(ifig)])
      eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_morris_cc_',int2str(ifig)]);
      eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_morris_cc_',int2str(ifig)]);
    end
  end


  SAMorris = [];
  for i=1:size(ac,2),
    [SAmeas, SAMorris(:,:,i)] = Morris_Measure_Groups(npT, [lpmat0 lpmat], ac(:,i),nliv);
  end
  %end
  SAac = squeeze(SAMorris(:,1,:))';
  figure,
  boxplot(SAac,'whis',10,'symbol','r.')
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
  set(gca,'xlim',[0.5 npT+0.5])
  ydum = get(gca,'ylim');
  set(gca,'ylim',[0 ydum(2)])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  for ip=1:npT,
    text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  title('EET All auto-correlations')
  saveas(gcf,[OutputDirectoryName,'\',fname_,'_morris_ac'])
  eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_morris_ac']);
  eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_morris_ac']);

  ifig = 0;
  for j=1:size(options_.varobs,1)
    if mod(j,6)==1
      figure('name',['EET auto-correlations']);
      ifig=ifig+1;
      iplo=0;
    end
    iplo=iplo+1;
    subplot(3,2,iplo)
    iv = find(ir_ac==j);
    if ~isempty(iv)
      if length(iv)>1
        boxplot(SAac(iv,:),'whis',3,'symbol','r.');
      else
        plot(SAac(iv,:),'r.');
      end
      set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
      set(gca,'xlim',[0.5 npT+0.5])
      ydum = get(gca,'ylim');
      set(gca,'ylim',[0 ydum(2)])
      for ip=1:npT,
        text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
      end
      xlabel(' ')
    end
    title(options_.varobs(j,:),'interpreter','none')
    if mod(j,6)==0 | j==size(options_.varobs,1)
      saveas(gcf,[OutputDirectoryName,'\',fname_,'_morris_ac_',int2str(ifig)])
      eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_morris_ac_',int2str(ifig)]);
      eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_morris_ac_',int2str(ifig)]);
    end
  end

  js=0;
  %for j=1:size(tadj,1),
  SAMorris = [];
  for i=1:size(tadj,2),
    js=js+1;
    [SAmeas, SAMorris(:,:,js)] = Morris_Measure_Groups(npT, [lpmat0 lpmat], tadj(:,i),nliv);
  end
  %end
  SAM = squeeze(SAMorris(nshock+1:end,1,:));
  for j=1:js,
    SAtadj(:,j)=SAM(:,j)./(max(SAM(:,j))+eps);
  end
  SAtadj = SAtadj';
  js=0;
  SAMorris = [];
  for i=1:size(iff,2),
    js=js+1;
    [SAmeas, SAMorriss(:,:,js)] = Morris_Measure_Groups(npT, [lpmat0 lpmat], iff(:,i),nliv);
  end
  SAM = squeeze(SAMorriss(nshock+1:end,1,:));
  for j=1:js,
    SAIF(:,j)=SAM(:,j)./(max(SAM(:,j))+eps);
  end
  SAIF = SAIF';

  figure,
  %bar(SAtadj),
  boxplot(SAtadj,'whis',10,'symbol','r.')
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
  set(gca,'xlim',[0.5 np+0.5])
  set(gca,'ylim',[0 1])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  for ip=1:np,
    text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  title('All half-life')
  saveas(gcf,[OutputDirectoryName,'\',fname_,'_morris_tadj'])
  eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_morris_tadj']);
  eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_morris_tadj']);

  ifig = 0;
  for j=1:size(options_.varobs,1)
    if mod(j,6)==1
      figure('name',['EET speed of adjustment observed variables']);
      ifig=ifig+1;
      iplo=0;
    end
    iplo=iplo+1;
    subplot(3,2,iplo)
    iv = find(ir_tadj==j);
    if ~isempty(iv)
      if length(iv)>1
        boxplot(SAtadj(iv,:),'whis',3,'symbol','r.');
      else
        plot(SAtadj(iv,:),'r.');
      end
      set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
      set(gca,'xlim',[0.5 np+0.5])
      ydum = get(gca,'ylim');
      set(gca,'ylim',[0 ydum(2)])
      for ip=1:np,
        text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
      end
      xlabel(' ')
    end
    title(options_.varobs(j,:),'interpreter','none')
    if mod(j,6)==0 | j==size(options_.varobs,1)
      saveas(gcf,[OutputDirectoryName,'\',fname_,'_morris_tadj_varobs_',int2str(ifig)])
      eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_morris_tadj_varobs_',int2str(ifig)]);
      eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_morris_tadj_varobs_',int2str(ifig)]);
    end
  end

  ifig = 0;
  for j=1:M_.exo_nbr,
    if mod(j,6)==1
      figure('name',['EET speed of adjustment shocks']);
      ifig=ifig+1;
      iplo=0;
    end
    iplo=iplo+1;
    subplot(3,2,iplo)
    iv = find(ic_tadj==j);
    if ~isempty(iv)
      if length(iv)>1
        boxplot(SAtadj(iv,:),'whis',3,'symbol','r.');
      else
        plot(SAtadj(iv,:),'r.');
      end
      set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
      set(gca,'xlim',[0.5 np+0.5])
      ydum = get(gca,'ylim');
      set(gca,'ylim',[0 ydum(2)])
      for ip=1:np,
        text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
      end
      xlabel(' ')
    end
    title(M_.exo_names(j,:),'interpreter','none')
    if mod(j,6)==0 | j==M_.exo_nbr,
      saveas(gcf,[OutputDirectoryName,'\',fname_,'_morris_tadj_exo_',int2str(ifig)])
      eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_morris_tadj_exo_',int2str(ifig)]);
      eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_morris_tadj_exo_',int2str(ifig)]);
    end
  end

  figure,
  %bar(SAIF),
  boxplot(SAIF,'whis',10,'symbol','r.')
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
  set(gca,'xlim',[0.5 np+0.5])
  set(gca,'ylim',[0 1])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  for ip=1:np,
    text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  ylabel('Elementary Effects')
  title('Steady state gains (impact factors)')
  saveas(gcf,[OutputDirectoryName,'\',fname_,'_morris_gain'])
  eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_morris_gain']);
  eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_morris_gain']);
  %figure, bar(SAIF'), title('All Gain Relationships')
  ifig = 0;
  for j=1:size(options_.varobs,1)
    if mod(j,6)==1
      figure('name',['EET steady state gain observed series']);
      ifig=ifig+1;
      iplo=0;
    end
    iplo=iplo+1;
    subplot(3,2,iplo)
    iv = find(ir_if==j);
    if ~isempty(iv)
      if length(iv)>1
        boxplot(SAIF(iv,:),'whis',10,'symbol','r.');
      else
        plot(SAIF(iv,:),'r.');
      end
      set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
      set(gca,'xlim',[0.5 np+0.5])
      ydum = get(gca,'ylim');
      set(gca,'ylim',[0 ydum(2)])
      for ip=1:np,
        text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
      end
      xlabel(' ')
    end
    title(options_.varobs(j,:),'interpreter','none')
    if mod(j,6)==0 | j==size(options_.varobs,1)
      saveas(gcf,[OutputDirectoryName,'\',fname_,'_morris_gain_varobs_',int2str(ifig)])
      eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_morris_gain_varobs_',int2str(ifig)]);
      eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_morris_gain_varobs_',int2str(ifig)]);
    end
  end

  ifig = 0;
  for j=1:M_.exo_nbr,
    if mod(j,6)==1
      figure('name',['EET steady state gain shocks']);
      ifig=ifig+1;
      iplo=0;
    end
    iplo=iplo+1;
    subplot(3,2,iplo)
    iv = find(ic_if==j);
    if ~isempty(iv)
      if length(iv)>1
        boxplot(SAIF(iv,:),'whis',3,'symbol','r.');
      else
        plot(SAIF(iv,:),'r.');
      end
      set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
      set(gca,'xlim',[0.5 np+0.5])
      ydum = get(gca,'ylim');
      set(gca,'ylim',[0 ydum(2)])
      for ip=1:np,
        text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
      end
      xlabel(' ')
    end
    title(M_.exo_names(j,:),'interpreter','none')
    if mod(j,6)==0 | j==M_.exo_nbr,
      saveas(gcf,[OutputDirectoryName,'\',fname_,'_morris_gain_exo_',int2str(ifig)])
      eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_morris_gain_exo_',int2str(ifig)]);
      eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_morris_gain_exo_',int2str(ifig)]);
    end
  end


  SAMorris = [];
  for j=1:j0,
    [SAmeas, SAMorris(:,:,j)] = Morris_Measure_Groups(npT, [lpmat0 lpmat], yt(:,j),nliv);
  end

  SAM = squeeze(SAMorris(nshock+1:end,1,:));
  for j=1:j0
    SAnorm(:,j)=SAM(:,j)./max(SAM(:,j));
    irex(j)=length(find(SAnorm(:,j)>0.01));
  end
  [dum, irel]=sort(irex);

  SAMmu = squeeze(SAMorris(:,2,:));
  for j=1:j0
    SAmunorm(:,j)=SAMmu(:,j)./max(SAM(:,j));  % normalised w.r.t. mu*
  end
  SAMsig = squeeze(SAMorris(:,3,:));
  for j=1:j0
    SAsignorm(:,j)=SAMsig(:,j)./max(SAMsig(:,j));
  end

  figure, %bar(SAnorm(:,irel))
  boxplot(SAnorm','whis',10,'symbol','r.')
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
  set(gca,'xlim',[0.5 np+0.5])
  set(gca,'ylim',[0 1])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  xlabel(' ')
  for ip=1:np,
    text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  title('Elementary effects parameters')
  saveas(gcf,[OutputDirectoryName,'\',fname_,'_morris_par'])
  eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_morris_par']);
  eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_morris_par']);

  figure, %bar(SAmunorm(:,irel))
  boxplot(SAmunorm','whis',10,'symbol','r.')
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
  set(gca,'xlim',[0.5 np+0.5])
  set(gca,'ylim',[-1 1])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  xlabel(' ')
  for ip=1:np,
    text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  title('\mu parameters')
  saveas(gcf,[OutputDirectoryName,'\',fname_,'_morrismu_par'])
  eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_morrismu_par']);
  eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_morrismu_par']);
  figure, %bar(SAsignorm(:,irel))
  boxplot(SAsignorm','whis',10,'symbol','r.')
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
  set(gca,'xlim',[0.5 np+0.5])
  set(gca,'ylim',[0 1])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  xlabel(' ')
  for ip=1:np,
    text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  title('\sigma parameters')
  saveas(gcf,[OutputDirectoryName,'\',fname_,'_morrissig_par'])
  eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_morrissig_par']);
  eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_morrissig_par']);

  %     figure, bar(SAnorm(:,irel)')
  %     set(gca,'xtick',[1:j0])
  %     set(gca,'xlim',[0.5 j0+0.5])
  %     title('Elementary effects relationships')
  %     saveas(gcf,[OutputDirectoryName,'\',fname_,'_morris_redform'])
  %     eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_morris_redform']);
  %     eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_morris_redform']);

elseif opt_gsa.morris==2,
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
else,

  if itrans==0,
    fsuffix = '';
  elseif itrans==1,
    fsuffix = '_log';
  else
    fsuffix = '_rank';
  end

  nrun=length(istable);
  nest=min(250,nrun);
  nfit=min(1000,nrun);

  for j=1:size(vdec,2),
    if itrans==0,
      y0 = vdec(istable,j);
    elseif itrans==1,
      y0 = log_trans_(vdec(istable,j));
    else
      y0 = trank(vdec(istable,j));
    end
  gsa_(j) = gsa_sdp(y0(1:nest), lpmat(istable(1:nest),:), ...
      2, [],[],[],0,[OutputDirectoryName,'/map_vdec',fsuffix,int2str(j)], pnames);
  if nfit>nest,
    gsa_(j) = gsa_sdp(y0(1:nest), lpmat(istable(1:nest),:), ...
        -2, gsa_(j).nvr*nest^3/nfit^3,[],[],0,[OutputDirectoryName,'/map_vdec',fsuffix,int2str(j)], pnames);
  end
    
    SAvdec(j,:)=gsa_(j).si;

  end

  figure,
  boxplot(SAvdec,'whis',10,'symbol','r.')
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
  set(gca,'xlim',[0.5 np+0.5])
  ydum = get(gca,'ylim');
  set(gca,'ylim',[0 ydum(2)])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  for ip=1:np,
    text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  title(['Main effects variance decomposition ',fsuffix],'interpreter','none')
  saveas(gcf,[OutputDirectoryName,'\',fname_,'_map_vdec',fsuffix])
  eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_map_vdec',fsuffix]);
  eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_map_vdec',fsuffix]);

  ifig = 0;
  for j=1:size(options_.varobs,1)
    if mod(j,6)==1
      figure('name',['Main effects observed variance decomposition ',fsuffix]);
      ifig=ifig+1;
      iplo=0;
    end
    iplo=iplo+1;
    subplot(3,2,iplo)
    iv = find(ir_vdec==j);
    if ~isempty(iv)
      if length(iv)>1
        boxplot(SAvdec(iv,:),'whis',10,'symbol','r.');
      else
        plot(SAvdec(iv,:),'r.');
      end
      set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
      set(gca,'xlim',[0.5 np+0.5])
      ydum = get(gca,'ylim');
      set(gca,'ylim',[0 ydum(2)])
      for ip=1:np,
        text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
      end
      xlabel(' ')
    end
    title(options_.varobs(j,:),'interpreter','none')
    if mod(j,6)==0 | j==size(options_.varobs,1)
      saveas(gcf,[OutputDirectoryName,'\',fname_,'_map_vdec',fsuffix,'_varobs_',int2str(ifig)])
      eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_map_vdec',fsuffix,'_varobs_',int2str(ifig)]);
      eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_map_vdec',fsuffix,'_varobs_',int2str(ifig)]);
    end
  end

  ifig = 0;
  for j=1:M_.exo_nbr,
    if mod(j,6)==1
      figure('name',['Main effects shocks variance decomposition ',fsuffix]);
      ifig=ifig+1;
      iplo=0;
    end
    iplo=iplo+1;
    subplot(3,2,iplo)
    iv = find(ic_vdec==j);
    if ~isempty(iv)
      if length(iv)>1
        boxplot(SAvdec(iv,:),'whis',3,'symbol','r.');
      else
        plot(SAvdec(iv,:),'r.');
      end
      set(gca,'xticklabel',' ','fontsize',3,'xtick',[1:np])
      set(gca,'xlim',[0.5 np+0.5])
      ydum = get(gca,'ylim');
      set(gca,'ylim',[0 ydum(2)])
      for ip=1:np,
        text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
      end
      xlabel(' ')
    end
    title(M_.exo_names(j,:),'interpreter','none')
    if mod(j,6)==0 | j==M_.exo_nbr
      saveas(gcf,[OutputDirectoryName,'\',fname_,'_map_vdec',fsuffix,'_exo_',int2str(ifig)])
      eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_map_vdec',fsuffix,'_exo_',int2str(ifig)]);
      eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_map_vdec',fsuffix,'_exo_',int2str(ifig)]);
    end
  end

  for j=1:size(cc,2),
    if itrans==0,
      y0 = cc(istable,j);
    elseif itrans==1,
      y0 = log_trans_(cc(istable,j));
    else
      y0 = trank(cc(istable,j));
    end
  gsa_(j) = gsa_sdp(y0(1:nest), lpmat(istable(1:nest),:), ...
      2, [],[],[],0,[OutputDirectoryName,'/map_cc',fsuffix,int2str(j)], pnames);
  if nfit>nest,
    gsa_(j) = gsa_sdp(y0(1:nest), lpmat(istable(1:nest),:), ...
        -2, gsa_(j).nvr*nest^3/nfit^3,[],[],0,[OutputDirectoryName,'/map_cc',fsuffix,int2str(j)], pnames);
  end
    SAcc(j,:)=gsa_(j).si;

  end

  figure,
  boxplot(SAcc,'whis',10,'symbol','r.')
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
  set(gca,'xlim',[0.5 np+0.5])
  ydum = get(gca,'ylim');
  set(gca,'ylim',[0 ydum(2)])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  for ip=1:np,
    text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  ylabel(' ')
  title(['Main effects cross-covariances ',fsuffix],'interpreter','none')
  saveas(gcf,[OutputDirectoryName,'\',fname_,'_map_cc',fsuffix])
  eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_map_cc',fsuffix]);
  eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_map_cc',fsuffix]);

  ifig = 0;
  for j=1:size(options_.varobs,1)
    if mod(j,6)==1
      figure('name',['Main effects cross-covariances ',fsuffix]);
      ifig=ifig+1;
      iplo=0;
    end
    iplo=iplo+1;
    subplot(3,2,iplo)
    iv = find(ir_cc==j);
    iv = [iv; find(ic_cc==j)];
    if ~isempty(iv)
      if length(iv)>1
        boxplot(SAcc(iv,:),'whis',10,'symbol','r.');
      else
        plot(SAcc(iv,:),'r.');
      end
      set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
      set(gca,'xlim',[0.5 np+0.5])
      ydum = get(gca,'ylim');
      set(gca,'ylim',[0 ydum(2)])
      for ip=1:np,
        text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
      end
      xlabel(' ')
    end
    title(options_.varobs(j,:),'interpreter','none')
    if mod(j,6)==0 | j==size(options_.varobs,1)
      saveas(gcf,[OutputDirectoryName,'\',fname_,'_map_cc',fsuffix,'_',int2str(ifig)])
      eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_map_cc',fsuffix,'_',int2str(ifig)]);
      eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_map_cc',fsuffix,'_',int2str(ifig)]);
    end
  end

  for j=1:size(ac,2),
    if itrans==0,
      y0 = ac(istable,j);
    elseif itrans==1,
      y0 = log_trans_(ac(istable,j));
    else
      y0 = trank(ac(istable,j));
    end
%     gsa_(j) = gsa_sdp_dyn( y0, lpmat(istable,:), ...
%       gsa_flag, [],[],[],0,[OutputDirectoryName,'/map_ac',fsuffix,int2str(j)], pnames);
  gsa_(j) = gsa_sdp(y0(1:nest), lpmat(istable(1:nest),:), ...
      2, [],[],[],0,[OutputDirectoryName,'/map_ac',fsuffix,int2str(j)], pnames);
  if nfit>nest,
    gsa_(j) = gsa_sdp(y0(1:nest), lpmat(istable(1:nest),:), ...
        -2, gsa_(j).nvr*nest^3/nfit^3,[],[],0,[OutputDirectoryName,'/map_ac',fsuffix,int2str(j)], pnames);
  end
    SAac(j,:)=gsa_(j).si;

  end

  figure,
  boxplot(SAac,'whis',10,'symbol','r.')
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
  set(gca,'xlim',[0.5 np+0.5])
  ydum = get(gca,'ylim');
  set(gca,'ylim',[0 ydum(2)])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  for ip=1:np,
    text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  title(['Main effects 1 lag auto-covariances ',fsuffix],'interpreter','none')
  saveas(gcf,[OutputDirectoryName,'\',fname_,'_map_ac',fsuffix])
  eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_map_ac',fsuffix]);
  eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_map_ac',fsuffix]);

  ifig = 0;
  for j=1:size(options_.varobs,1)
    if mod(j,6)==1
      figure('name',['Main effects auto-covariances ',fsuffix]);
      ifig=ifig+1;
      iplo = 0;
    end
    iplo=iplo+1;
    subplot(3,2,iplo)
    iv = find(ir_ac==j);
    %iv = [iv; find(ic_ac==j)];
    if ~isempty(iv)
      if length(iv)>1
        boxplot(SAac(iv,:),'whis',10,'symbol','r.');
      else
        plot(SAac(iv,:),'r.');
      end
      set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
      set(gca,'xlim',[0.5 np+0.5])
      ydum = get(gca,'ylim');
      set(gca,'ylim',[0 ydum(2)])
      for ip=1:np,
        text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
      end
      xlabel(' ')
    end
    title(options_.varobs(j,:),'interpreter','none')
    if mod(j,6)==0 | j==size(options_.varobs,1)
      saveas(gcf,[OutputDirectoryName,'\',fname_,'_map_ac',fsuffix,'_',int2str(ifig)])
      eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_map_ac',fsuffix,'_',int2str(ifig)]);
      eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_map_ac',fsuffix,'_',int2str(ifig)]);
    end
  end

  for j=1:size(tadj,2),
    if itrans==0,
      y0 = tadj(istable,j);
    elseif itrans==1,
      y0 = log_trans_(tadj(istable,j));
    else
      y0 = trank(tadj(istable,j));
    end
%     gsa_(j) = gsa_sdp_dyn( y0, lpmat(istable,:), ...
%       gsa_flag, [],[],[],0,[OutputDirectoryName,'/map_tadj',fsuffix,int2str(j)], pnames);
  gsa_(j) = gsa_sdp(y0(1:nest), lpmat(istable(1:nest),:), ...
      2, [],[],[],0,[OutputDirectoryName,'/map_tadj',fsuffix,int2str(j)], pnames);
  if nfit>nest,
    gsa_(j) = gsa_sdp(y0(1:nest), lpmat(istable(1:nest),:), ...
        -2, gsa_(j).nvr*nest^3/nfit^3,[],[],0,[OutputDirectoryName,'/map_tadj',fsuffix,int2str(j)], pnames);
  end
    SAtadj(j,:)=gsa_(j).si;

  end

  figure,
  boxplot(SAtadj,'whis',10,'symbol','r.')
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
  set(gca,'xlim',[0.5 np+0.5])
  ydum = get(gca,'ylim');
  set(gca,'ylim',[0 ydum(2)])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  for ip=1:np,
    text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  title(['Main effects speed of adjustment ',fsuffix],'interpreter','none')
  saveas(gcf,[OutputDirectoryName,'\',fname_,'_map_tadj',fsuffix])
  eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_map_tadj',fsuffix]);
  eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_map_tadj',fsuffix]);

  ifig = 0;
  for j=1:size(options_.varobs,1)
    if mod(j,6)==1
      figure('name',['Main effects observed speed adjustment ',fsuffix]);
      ifig=ifig+1;
      iplo = 0;
    end
    iplo=iplo+1;
    subplot(3,2,iplo)
    iv = find(ir_tadj==j);
    if ~isempty(iv)
      if length(iv)>1
        boxplot(SAtadj(iv,:),'whis',3,'symbol','r.');
      else
        plot(SAtadj(iv,:),'r.');
      end
      set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
      set(gca,'xlim',[0.5 np+0.5])
      ydum = get(gca,'ylim');
      set(gca,'ylim',[0 ydum(2)])
      for ip=1:np,
        text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
      end
      xlabel(' ')
    end
    title(options_.varobs(j,:),'interpreter','none')
    if mod(j,6)==0 | j==size(options_.varobs,1)
      saveas(gcf,[OutputDirectoryName,'\',fname_,'_map_tadj',fsuffix,'_varobs_',int2str(ifig)])
      eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_map_tadj',fsuffix,'_varobs_',int2str(ifig)]);
      eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_map_tadj',fsuffix,'_varobs_',int2str(ifig)]);
    end
  end

  ifig = 0;
  for j=1:M_.exo_nbr,
    if mod(j,6)==1
      figure('name',['Main effects shocks speed of adjustment ',fsuffix]);
      ifig=ifig+1;
      iplo=0;
    end
    iplo=iplo+1;
    subplot(3,2,iplo)
    iv = find(ic_tadj==j);
    if ~isempty(iv)
      if length(iv)>1
        boxplot(SAtadj(iv,:),'whis',3,'symbol','r.');
      else
        plot(SAtadj(iv,:),'r.');
      end
      set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
      set(gca,'xlim',[0.5 np+0.5])
      ydum = get(gca,'ylim');
      set(gca,'ylim',[0 ydum(2)])
      for ip=1:np,
        text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
      end
      xlabel(' ')
    end
    title(M_.exo_names(j,:),'interpreter','none')
    if mod(j,6)==0 | j==M_.exo_nbr,
      saveas(gcf,[OutputDirectoryName,'\',fname_,'_map_tadj',fsuffix,'_exo_',int2str(ifig)])
      eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_map_tadj',fsuffix,'_exo_',int2str(ifig)]);
      eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_map_tadj',fsuffix,'_exo_',int2str(ifig)]);
    end
  end


  for j=1:size(iff,2),
    if itrans==0,
      y0 = iff(istable,j);
    elseif itrans==1,
      y0 = log_trans_(iff(istable,j));
    else
      y0 = trank(iff(istable,j));
    end
%     gsa_(j) = gsa_sdp_dyn( y0, lpmat(istable,:), ...
%       gsa_flag, [],[],[],0,[OutputDirectoryName,'/map_if',fsuffix,int2str(j)], pnames);
  gsa_(j) = gsa_sdp(y0(1:nest), lpmat(istable(1:nest),:), ...
      2, [],[],[],0,[OutputDirectoryName,'/map_if',fsuffix,int2str(j)], pnames);
  if nfit>nest,
    gsa_(j) = gsa_sdp(y0(1:nest), lpmat(istable(1:nest),:), ...
        -2, gsa_(j).nvr*nest^3/nfit^3,[],[],0,[OutputDirectoryName,'/map_if',fsuffix,int2str(j)], pnames);
  end
    SAif(j,:)=gsa_(j).si;

  end

  figure,
  boxplot(SAif,'whis',10,'symbol','r.')
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
  set(gca,'xlim',[0.5 np+0.5])
  ydum = get(gca,'ylim');
  set(gca,'ylim',[0 ydum(2)])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  for ip=1:np,
    text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  title(['Main effects impact factors ',fsuffix],'interpreter','none')
  saveas(gcf,[OutputDirectoryName,'\',fname_,'_map_if',fsuffix])
  eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_map_if',fsuffix]);
  eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_map_if',fsuffix]);

  ifig = 0;
  for j=1:size(options_.varobs,1)
    if mod(j,6)==1
      figure('name',['Main effects observed impact factors ',fsuffix]);
      ifig=ifig+1;
      iplo = 0;
    end
    iplo=iplo+1;
    subplot(3,2,iplo)
    iv = find(ir_if==j);
    if ~isempty(iv)
      if length(iv)>1
        boxplot(SAif(iv,:),'whis',3,'symbol','r.');
      else
        plot(SAif(iv,:),'r.');
      end
      set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
      set(gca,'xlim',[0.5 np+0.5])
      ydum = get(gca,'ylim');
      set(gca,'ylim',[0 ydum(2)])
      for ip=1:np,
        text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
      end
      xlabel(' ')
    end
    title(options_.varobs(j,:),'interpreter','none')
    if mod(j,6)==0 | j==size(options_.varobs,1)
      saveas(gcf,[OutputDirectoryName,'\',fname_,'_map_if',fsuffix,'_varobs_',int2str(ifig)])
      eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_map_if',fsuffix,'_varobs_',int2str(ifig)]);
      eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_map_if',fsuffix,'_varobs_',int2str(ifig)]);
    end
  end

  ifig = 0;
  for j=1:M_.exo_nbr,
    if mod(j,6)==1
      figure('name',['Main effects shocks impact factors ',fsuffix]);
      ifig=ifig+1;
      iplo=0;
    end
    iplo=iplo+1;
    subplot(3,2,iplo)
    iv = find(ic_if==j);
    if ~isempty(iv)
      if length(iv)>1
        boxplot(SAif(iv,:),'whis',3,'symbol','r.');
      else
        plot(SAif(iv,:),'r.');
      end
      set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
      set(gca,'xlim',[0.5 np+0.5])
      ydum = get(gca,'ylim');
      set(gca,'ylim',[0 ydum(2)])
      for ip=1:np,
        text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
      end
      xlabel(' ')
    end
    title(M_.exo_names(j,:),'interpreter','none')
    if mod(j,6)==0 | j==M_.exo_nbr
      saveas(gcf,[OutputDirectoryName,'\',fname_,'_map_if',fsuffix,'_exo_',int2str(ifig)])
      eval(['print -depsc2 ' OutputDirectoryName '\' fname_ '_map_if',fsuffix,'_exo_',int2str(ifig)]);
      eval(['print -dpdf ' OutputDirectoryName '\' fname_ '_map_if',fsuffix,'_exo_',int2str(ifig)]);
    end
  end

end
