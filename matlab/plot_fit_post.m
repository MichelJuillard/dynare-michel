function plot_fit_post(varargin);
%function [rmse, lgobs_] = plot_fit(varargin);
global options_ bayestopt_ M_ oo_

dirname = CheckPath('Output');


ys_ = oo_.steady_state;
lgy_ = M_.endo_names;
fname_ = M_.fname;

if exist(options_.datafile)
  instr = options_.datafile;
else
  instr = ['load ' options_.datafile];
end
eval(instr);
if ~exist('T','var'), 
  temp = eval(deblank(options_.varobs(1,:)));
  T=[1:length(temp)]'; 
  clear temp;
end
istart = options_.presample+1;

if isempty(varargin),
  lgobs_= options_.varobs;
  mfys=bayestopt_.mfys;
else
  mfys=[];
  for j=1:length(varargin),
    dum = strmatch(varargin{j},lgy_,'exact');
    mfys = [mfys dum];
    if j==1,
      lgobs_ = varargin{j};      
    else
      lgobs_ = str2mat(lgobs_,varargin{j});
    end
  end
end

nobs = size(lgobs_,1);
fobs = options_.first_obs;
if options_.loglinear == 1
  constant = log(ys_(mfys));
else
  constant = ys_(mfys);
end

trend = constant*ones(1,options_.nobs);

%disp(' ')
%disp(' ')
%disp('            RMSE')
y00=getFilterRuns(lgobs_);
if isempty(y00), return, end,
x = T(fobs+istart-1:fobs+options_.nobs-1);
for ifig = 1:ceil(size(lgobs_,1)/9),
  h1 = figure('Name',['Bayesian DSGE 1-step ahead prediction']);
  h2 = figure('Name',['Bayesian DSGE Innovations']);
  for j=1+9*(ifig-1):min(9*ifig,size(lgobs_,1)),
    figure(h1)
    subplot(3,3,j-9*(ifig-1)),
    vj = deblank(lgobs_(j,:));
    i = strmatch(vj,lgy_,'exact');
    for k=1:options_.nobs,
    [y0(k,1),Median(k,1),Var(k,1),HPD(k,:),Distrib(k,:)] = ...
          posterior_moments(y00(j,k,:),0);
    end
    yobs=eval(vj);
%     patch([x, x(end:-1:1)], [HPD(istart:end ,1)', HPD(end:-1:istart,2)'],[0.75 0.75 0.75]);
%     hold on,
    plot(x, [y0(istart:end) yobs(fobs+istart-1:fobs+options_.nobs-1)])
    if options_.TeX,
      tnam=deblank(M_.endo_names_tex(i,:));
      if strmatch('\log',tnam),
        tnam=tnam(2:end);
      end
      title(tnam,'interpreter','tex')
    else
      title(vj,'interpreter','none')
    end
    vv(:,j) = y0(istart:end)-yobs(fobs+istart-1:fobs+options_.nobs-1);
    hh=get(gca,'children');
    set(hh(1),'linestyle',':','color',[0 0 0])
    set(hh(2),'linestyle','-','color',[0 0 0])
    set(gca,'xlim',[T(fobs)-1 T(fobs+options_.nobs-1)+1])
    hold off,
    %disp([lgobs_(j,:), sprintf('%15.5g',[rmse(j)'])])
    figure(h2)
    subplot(3,3,j-9*(ifig-1)),
    plot(x, vv(:,j),'k')
    if options_.TeX,
      title(tnam,'interpreter','tex')
    else
      title(vj,'interpreter','none')
    end
    set(gca,'xlim',[T(fobs)-1 T(fobs+options_.nobs-1)+1])
  end
  saveas(h1,[dirname,'\',fname_,'_Bayesian_Fit',int2str(ifig)])
  eval(['print -depsc2 ' dirname,'\',fname_,'_Bayesian_Fit',int2str(ifig)]);
  eval(['print -dpdf ' dirname,'\',fname_,'_Bayesian_Fit',int2str(ifig)]);
  saveas(h2,[dirname,'\',fname_,'_Bayesian_Innovations',int2str(ifig)])
  if options_.nograph
    close(h1)
    close(h2)
  end
end
