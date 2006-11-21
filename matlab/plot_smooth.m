function plot_smooth(varargin)
global options_ M_ oo_

ys_ = oo_.steady_state;
lgy_ = M_.endo_names;
fname_ = M_.fname;
texname=M_.endo_names_tex;

ifig0=0;
if strmatch('-',varargin{end})
    optn = varargin{end};
    varargin=varargin(1:end-1);
    if strmatch('-new',optn)
        ifig0=length(dir([fname_,'_SmoothedUnobserved*.fig']));
    end
else
    afig=dir([fname_,'_SmoothedUnobserved*.fig']);
    for j=1:length(afig),
        delete(afig(j).name);
    end
end
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

fobs = options_.first_obs;
nv=length(varargin);
for ifig = 1:ceil(nv/9),
    figure('Name','Smoothed unobserved variables'),
    for j=1+9*(ifig-1):min(9*ifig,nv),
        subplot(3,3,j-9*(ifig-1)),
        i = strmatch(varargin{j},lgy_,'exact');
        if ~isempty(texname)
        if strmatch('\log',texname(i,:));
          texname(i,1:end-1)=texname(i,2:end);
        end
        end
        eval(['plot(T(fobs:fobs+options_.nobs-1),[oo_.SmoothedVariables.',varargin{j},'(1:end)+ys_(i)],''-k'')']) 
        if exist(varargin{j},'var')
            hold on,
            eval(['plot(T(fobs:fobs+options_.nobs-1),',varargin{j},'(fobs:fobs+options_.nobs-1),'':k'')'])         
        end
        
        if options_.TeX,
          title(deblank(texname(i,:)),'interpreter','tex')
        else  
          title(varargin{j},'interpreter','none')
        end
        set(gca,'xlim',[T(fobs)-1 T(fobs+options_.nobs-1)+1])
    end
    saveas(gcf,[fname_,'_SmoothedUnobserved',int2str(ifig+ifig0)])
    eval(['print -depsc2 ' fname_,'_SmoothedUnobserved',int2str(ifig+ifig0)]);
    eval(['print -dpdf ' fname_,'_SmoothedUnobserved',int2str(ifig+ifig0)]);
    if options_.nograph
        close(gcf)
    end
end
