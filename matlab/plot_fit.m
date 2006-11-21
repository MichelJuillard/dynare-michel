function [rmse, lgobs_, vv] = plot_fit(varargin);
%function [rmse, lgobs_] = plot_fit(varargin);
global options_ bayestopt_ M_ oo_

ys_ = oo_.steady_state;
lgy_ = M_.endo_names;
fname_ = M_.fname;
texname = M_.endo_names_tex;

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
for ifig = 1:ceil(size(lgobs_,1)/9),
  h1 = figure('Name',['DSGE 1-step ahead prediction']);
  h2 = figure('Name',['DSGE Innovations']);
  for j=1+9*(ifig-1):min(9*ifig,size(lgobs_,1)),
    figure(h1)
    subplot(3,3,j-9*(ifig-1)),
    vj = deblank(lgobs_(j,:));
    i = strmatch(vj,lgy_,'exact');
    if ~isempty(texname),
    if strmatch('\log',texname(i,:));
      texname(i,1:end-1)=texname(i,2:end);
    end
    end
    eval(['plot(T(fobs+istart-1:fobs+options_.nobs-1), [oo_.FilteredVariables.', ...
        vj,'(istart:end-1)+trend(j,istart:end)'' oo_.SmoothedVariables.',...
        vj,'(istart:end)+trend(j,istart:end)''])'])
    if options_.TeX,
      title(deblank(texname(i,:)),'interpreter','tex')
    else  
      title(vj,'interpreter','none')
    end
    eval(['rmse(j) = sqrt(mean((oo_.SmoothedVariables.',vj,'(istart:end)-oo_.FilteredVariables.',vj,'(istart:end-1)).^2));'])
    eval(['vv(:,j) = (oo_.SmoothedVariables.',vj,'(istart:end)-oo_.FilteredVariables.',vj,'(istart:end-1));'])
    hh=get(gca,'children');
    set(hh(1),'linestyle',':','color',[0 0 0])
    set(hh(2),'linestyle','-','color',[0 0 0])
    set(gca,'xlim',[T(fobs)-1 T(fobs+options_.nobs-1)+1])
    %disp([lgobs_(j,:), sprintf('%15.5g',[rmse(j)'])])
    figure(h2)
    subplot(3,3,j-9*(ifig-1)),
    plot(T(fobs+istart-1:fobs+options_.nobs-1), vv(:,j),'k')
    if options_.TeX,
      title(deblank(texname(i,:)),'interpreter','tex')
    else  
      title(vj,'interpreter','none')
    end
    set(gca,'xlim',[T(fobs)-1 T(fobs+options_.nobs-1)+1])
  end
  saveas(h1,[fname_,'_Fit',int2str(ifig)])
  eval(['print -depsc2 ' fname_,'_Fit',int2str(ifig)]);
  eval(['print -dpdf ' fname_,'_Fit',int2str(ifig)]);
  saveas(h2,[fname_,'_Innovations',int2str(ifig)])
  eval(['print -depsc2 ' fname_,'_Innovations',int2str(ifig)]);
  eval(['print -dpdf ' fname_,'_Innovations',int2str(ifig)]);
  if options_.nograph
    close(h1)
    close(h2)
  end
end
