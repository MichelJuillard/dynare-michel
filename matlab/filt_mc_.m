function [rmse_MC, ixx] = filt_mc_(vvarvecm, loadSA, pfilt, alpha, alpha2, istart, alphaPC)
% copyright Marco Ratto 2006
global bayestopt_ estim_params_ M_ options_ oo_

if nargin<1 | isempty(vvarvecm),
  vvarvecm = options_.varobs;
end
if nargin<2,
  loadSA=0;
end
if nargin<3 | isempty(pfilt),
  pfilt=0.1;  % cut the best 10% of runs
end
if nargin<4 | isempty(alpha),
  alpha=0.002;
end
if nargin<5 | isempty(alpha2),
  alpha2=0.5;
end
if nargin<6 | isempty(istart),
  istart=1;
end
if nargin<7,
  alphaPC=0.5;
end

fname_ = M_.fname;
lgy_ = M_.endo_names;
dr_ = oo_.dr;

disp(' ')
disp(' ')
disp('Starting sensitivity analysis')
disp('for the fit of EACH observed series ...')
disp(' ')
disp('Deleting old SA figures...')
a=dir('*.fig');
if options_.opt_gsa.pprior
  tmp=['_SA_fit_prior'];
else
  tmp=['_SA_fit_mc'];
end
for j=1:length(a), 
  if strmatch([fname_,tmp],a(j).name), 
    disp(a(j).name)
    delete(a(j).name)
  end, 
end
disp('done !')


nshock=estim_params_.nvx + estim_params_.nvn + estim_params_.ncx + estim_params_.ncn;
npar=estim_params_.np;
for j=1:npar+nshock,
  if j>nshock
    xparam1(j)=oo_.posterior_mode.parameters.(bayestopt_.name{j});
  else
    xparam1(j)=oo_.posterior_mode.shocks_std.(bayestopt_.name{j});
  end
end

if options_.opt_gsa.pprior
  fnamtmp=[fname_,'_prior'];
else
  fnamtmp=[fname_,'_mc'];      
end

if ~loadSA,
  set_all_parameters(xparam1);
  steady_;
  eval(options_.datafile)
  load(fnamtmp);
  nruns=size(x,1);
  nfilt=floor(pfilt*nruns);
  disp(' ')
  disp('Computing RMSE''s...')
  fobs = options_.first_obs;
  nobs=options_.nobs;
  pp= corrcoef(x);
  for i=1:size(vvarvecm,1),
    vj=deblank(vvarvecm(i,:));
    if options_.prefilter == 1
      eval([vj,'=',vj,'-bayestopt_.mean_varobs(i);'])
    end
    
    jxj = strmatch(vj,lgy_(dr_.order_var,:),'exact');
    js = strmatch(vj,lgy_,'exact');
    eval(['rmse_mode(i) = sqrt(mean((',vj,'(fobs-1+istart:fobs-1+nobs)-oo_.steady_state(js)-oo_.FilteredVariables.',vj,'(istart:end-1)).^2));'])
    y0=zeros(nobs+1,nruns);
    nb = size(stock_filter,3);
    y0 = squeeze(stock_filter(:,jxj,:)) + ...
      kron(stock_ys(js,:),ones(size(stock_filter,1),1));
    y0M=mean(y0,2);
    for j=1:nruns,
      eval(['rmse_MC(j,i) = sqrt(mean((',vj,'(fobs-1+istart:fobs-1+nobs)-y0(istart:end-1,j)).^2));'])
    end
    eval(['rmse_pmean(j,i) = sqrt(mean((',vj,'(fobs-1+istart:fobs-1+nobs)-y0M(istart:end-1)).^2));'])
  end
  disp('... done!')
  
  save(fnamtmp, 'rmse_MC', 'rmse_mode','-append')    
else
  tmp=load(fnamtmp);
  x=tmp.x;
  logpo2=tmp.logpo2;
  rmse_MC=tmp.rmse_MC;
  rmse_mode=tmp.rmse_mode;
  clear tmp
  nruns=size(x,1);
  nfilt=floor(pfilt*nruns);
end
% smirnov tests
[dum, ipost]=sort(logpo2);
for i=1:size(vvarvecm,1),
  [dum, ixx(:,i)]=sort(rmse_MC(:,i));
  for j=1:npar+nshock,
    [H,P,KSSTAT] = smirnov(x(ixx(nfilt+1:end,i),j),x(ixx(1:nfilt,i),j), alpha);
    [H1,P1,KSSTAT1] = smirnov(x(ixx(nfilt+1:end,i),j),x(ixx(1:nfilt,i),j),alpha,1);
    [H2,P2,KSSTAT2] = smirnov(x(ixx(nfilt+1:end,i),j),x(ixx(1:nfilt,i),j),alpha,-1);
    if H1 & H2==0,
      SS(j,i)=1;      
    elseif H1==0,
      SS(j,i)=-1;     
    else
      SS(j,i)=0;     
    end
    PP(j,i)=P;
  end
end

param_names='';
for j=1:npar+nshock,
  param_names=str2mat(param_names, bayestopt_.name{j});
end
param_names=param_names(2:end,:);

disp(' ')
disp('RMSE over the MC sample:')
disp('            min yr RMSE    max yr RMSE')
for j=1:size(vvarvecm,1),
  disp([vvarvecm(j,:), sprintf('%15.5g',[(min(rmse_MC(:,j))) [(max(rmse_MC(:,j)))]])])
end
invar = find( std(rmse_MC)./mean(rmse_MC)<=0.0001 );
if ~isempty(invar)
  disp(' ')
  disp(' ')
  disp('RMSE is not varying significantly over the MC sample for the following variables:')
  disp(vvarvecm(invar,:))
  disp('These variables are excluded from SA')
  disp('[Unless you treat these series as exogenous, there is something wrong in your estimation !]')
end
ivar = find( std(rmse_MC)./mean(rmse_MC)>0.0001 );
vvarvecm=vvarvecm(ivar,:);
rmse_MC=rmse_MC(:,ivar);

disp(' ')
disp(['Sample filtered the ',num2str(pfilt*100),'% best RMSE''s for each observed series ...' ])

% figure, boxplot(rmse_MC)
% set(gca,'xticklabel',vvarvecm)
% saveas(gcf,[fname_,'_SA_RMSE'])

disp(' ')
disp(' ')
disp('RMSE ranges after filtering:')
disp(['             best ',num2str(pfilt*100),'% filtered             remaining 90%'])
disp(['             min            max            min            max            posterior mode'])
for j=1:size(vvarvecm,1),
  disp([vvarvecm(j,:), sprintf('%15.5g',[min(rmse_MC(ixx(1:nfilt,j),j)) ...
        max(rmse_MC(ixx(1:nfilt,j),j))  ...
        min(rmse_MC(ixx(nfilt+1:end,j),j)) ...
        max(rmse_MC(ixx(nfilt+1:end,j),j)) ...
        rmse_mode(j)])])
end

%stab_map_1(x, ipost(1:nfilt), ipost(nfilt+1:end), 'SA_post', 1);
%stab_map_2(x(ipost(1:nfilt),:),alpha2,'SA_post', 1);
% for i=1:size(vvarvecm,1),
%     aname=['SA_fit_ALL_',deblank(vvarvecm(i,:))];        
%     stab_map_1(x, ixx(1:nfilt,i), ixx(nfilt+1:end,i), aname, 1);
%     close all
% end

SP=zeros(npar+nshock,size(vvarvecm,1));
for j=1:size(vvarvecm,1),
  ns=find(PP(:,j)<alpha);
  SP(ns,j)=ones(size(ns));
  SS(:,j)=SS(:,j).*SP(:,j);
end

for j=1:npar+nshock, %estim_params_.np,
  nsp(j)=length(find(SP(j,:)));
end
snam0=param_names(find(nsp==0),:);
snam1=param_names(find(nsp==1),:);
snam2=param_names(find(nsp>1),:);
snam=param_names(find(nsp>0),:);
% snam0=bayestopt_.name(find(nsp==0));
% snam1=bayestopt_.name(find(nsp==1));
% snam2=bayestopt_.name(find(nsp>1));
% snam=bayestopt_.name(find(nsp>0));
nsnam=(find(nsp>1));

disp(' ')
disp(' ')
disp('These parameters do not affect significantly the fit of ANY observed series:')
disp(snam0)
disp(' ')
disp('These parameters affect ONE single observed series:')
disp(snam1)
disp(' ')
disp('These parameters affect MORE THAN ONE observed series: trade off exists!')
disp(snam2)


%pnam=bayestopt_.name(end-estim_params_.np+1:end);
pnam=bayestopt_.name;

% plot trade-offs
a00=jet(size(vvarvecm,1));
for ix=1:ceil(length(nsnam)/6),
  figure,
  for j=1+6*(ix-1):min(size(snam2,1),6*ix),
    subplot(2,3,j-6*(ix-1))
    %h0=cumplot(x(:,nsnam(j)+nshock));
    h0=cumplot(x(:,nsnam(j)));
    %set(h0,'color',[1 1 1])
    hold on,
    np=find(SP(nsnam(j),:));
    %a0=jet(nsp(nsnam(j)));
    a0=a00(np,:);
    for i=1:nsp(nsnam(j)), %size(vvarvecm,1),
      %h0=cumplot(x(ixx(1:nfilt,np(i)),nsnam(j)+nshock));
      h0=cumplot(x(ixx(1:nfilt,np(i)),nsnam(j)));
      set(h0,'color',a0(i,:))
    end
    ydum=get(gca,'ylim');
    %xdum=xparam1(nshock+nsnam(j));
    xdum=xparam1(nsnam(j));
    h1=plot([xdum xdum],ydum);
    set(h1,'color',[0.85 0.85 0.85],'linewidth',2)
    h0=legend(str2mat('base',vvarvecm(np,:)),0); 
    %h0=legend({'base',vnam{np}}',0); 
    xlabel('')
    set(findobj(get(h0,'children'),'type','text'),'interpreter','none')
    title([pnam{nsnam(j)}],'interpreter','none')
    if options_.opt_gsa.pprior
      saveas(gcf,[fname_,'_SA_fit_prior_',num2str(ix)])
    else
      saveas(gcf,[fname_,'_SA_fit_mc_',num2str(ix)])
    end
  end
end
close all

for j=1:size(SP,2),
  nsx(j)=length(find(SP(:,j)));
end

number_of_grid_points = 2^9;      % 2^9 = 512 !... Must be a power of two.
bandwidth = 0;                    % Rule of thumb optimal bandwidth parameter.
kernel_function = 'gaussian';     % Gaussian kernel for Fast Fourrier Transform approximaton.  
%kernel_function = 'uniform';     % Gaussian kernel for Fast Fourrier Transform approximaton.  

% for j=1:size(SP,2),
%     nfig=0;
%     np=find(SP(:,j));
%     for i=1:nsx(j), %size(vvarvecm,1),
%         if mod(i,12)==1,
%             nfig=nfig+1;
%             %figure('name',['Sensitivity of fit of ',vnam{j}]),
%             figure('name',['Sensitivity of fit of ',deblank(vvarvecm(j,:)),' ',num2str(nfig)]),
%         end
%         
%         subplot(3,4,i-12*(nfig-1))
%         optimal_bandwidth = mh_optimal_bandwidth(x(ixx(1:nfilt,j),np(i)),nfilt,bandwidth,kernel_function); 
%         [x1,f1] = kernel_density_estimate(x(ixx(1:nfilt,j),np(i)),number_of_grid_points,...
%             optimal_bandwidth,kernel_function);
%         plot(x1, f1,':k','linewidth',2)
%         optimal_bandwidth = mh_optimal_bandwidth(x(ixx(nfilt+1:end,j),np(i)),nruns-nfilt,bandwidth,kernel_function); 
%         [x1,f1] = kernel_density_estimate(x(ixx(nfilt+1:end,j),np(i)),number_of_grid_points,...
%             optimal_bandwidth,kernel_function);
%         hold on, plot(x1, f1,'k','linewidth',2)
%         ydum=get(gca,'ylim');
%         %xdum=xparam1(nshock+np(i));
%         xdum=xparam1(np(i));
%         h1=plot([xdum xdum],ydum);
%         set(h1,'color',[0.85 0.85 0.85],'linewidth',2)
%         %xdum1=mean(x(ixx(1:nfilt,j),np(i)+nshock));
%         xdum1=mean(x(ixx(1:nfilt,j),np(i)));
%         h2=plot([xdum1 xdum1],ydum);
%         set(h2,'color',[0 1 0],'linewidth',2)
%         %         h0=cumplot(x(nfilt+1:end,np(i)+nshock));
%         %         set(h0,'color',[1 1 1])
%         %         hold on,
%         %         h0=cumplot(x(ixx(1:nfilt,j),np(i)+nshock));
%         %         set(h0,'linestyle',':','color',[1 1 1])
%         %title([pnam{np(i)}])
%         title([pnam{np(i)},'. K-S prob ', num2str(PP(np(i),j))],'interpreter','none')
%         xlabel('')
%         if mod(i,12)==0 | i==nsx(j),
%             saveas(gcf,[fname_,'_SA_fit_',deblank(vvarvecm(j,:)),'_',int2str(nfig)])
%             close(gcf)
%         end
%     end
% end


disp(' ')
disp(' ')
disp('Sensitivity table (significance and direction):')
vav=char(zeros(1, size(param_names,2)+3 ));
ibl = 12-size(vvarvecm,2);
for j=1:size(vvarvecm,1), 
  vav = [vav, char(zeros(1,ibl)),vvarvecm(j,:)];
end
disp(vav)
for j=1:npar+nshock, %estim_params_.np,
  %disp([param_names(j,:), sprintf('%8.5g',SP(j,:))])    
  disp([param_names(j,:),'   ', sprintf('%12.3g',PP(j,:))])    
  disp([char(zeros(1, size(param_names,2)+3 )),sprintf('    (%6g)',SS(j,:))])    
end


disp(' ')
disp(' ')
disp('Starting bivariate analysis:')

for i=1:size(vvarvecm,1)
  if options_.opt_gsa.pprior
    fnam = ['SA_fit_prior_',deblank(vvarvecm(i,:))];
  else
    fnam = ['SA_fit_mc_',deblank(vvarvecm(i,:))];
  end
  stab_map_2(x(ixx(1:nfilt,i),:),alpha2,fnam, 1);
  
  %     [pc,latent,explained] = pcacov(c0);
  %     %figure, bar([explained cumsum(explained)])
  %     ifig=0;
  %     j2=0;
  %     for j=1:npar+nshock,
  %         i2=find(abs(pc(:,j))>alphaPC);
  %         if ~isempty(i2),
  %             j2=j2+1;
  %             if mod(j2,12)==1,
  %                 ifig=ifig+1;
  %                 figure('name',['PCA of the filtered sample ',deblank(vvarvecm(i,:)),' ',num2str(ifig)]),
  %             end
  %             subplot(3,4,j2-(ifig-1)*12)
  %             bar(pc(i2,j)), 
  %             set(gca,'xticklabel',bayestopt_.name(i2)), 
  %             set(gca,'xtick',[1:length(i2)])
  %             title(['PC ',num2str(j),'. Explained ',num2str(explained(j)),'%'])
  %         end
  %         if (mod(j2,12)==0 | j==(npar+nshock)) & j2,
  %             saveas(gcf,[fname_,'_SA_PCA_',deblank(vvarvecm(i,:)),'_',int2str(ifig)])
  %         end
  %     end
  %     close all
end

