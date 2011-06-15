function gsa_plotmatrix(type,varargin)
global bayestopt_ options_ M_

RootDirectoryName = CheckPath('GSA');

if options_.opt_gsa.pprior
    load([ RootDirectoryName filesep  M_.fname '_prior.mat'],'lpmat0','lpmat','istable','iunstable','iindeterm','iwrong')
else
    load([ RootDirectoryName filesep  M_.fname '_mc.mat'],'lpmat0','lpmat','istable','iunstable','iindeterm','iwrong')
    eval(['load ' options_.mode_file ' xparam1;']');
end

iexplosive = iunstable(~ismember(iunstable,[iindeterm;iwrong]));

switch type
    case 'all'
        x=[lpmat0 lpmat];
        NumberOfDraws=size(x,1);
        B=NumberOfDraws;
    case 'stable'
        x=[lpmat0(istable,:) lpmat(istable,:)];
        NumberOfDraws=size(x,1);
        B=NumberOfDraws;
    case 'nosolution'
        x=[lpmat0(iunstable,:) lpmat(iunstable,:)];
        NumberOfDraws=size(x,1);
        B=NumberOfDraws;
    case 'unstable'
        x=[lpmat0(iexplosive,:) lpmat(iexplosive,:)];
        NumberOfDraws=size(x,1);
        B=NumberOfDraws;
    case 'indeterm'
        x=[lpmat0(iindeterm,:) lpmat(iindeterm,:)];
        NumberOfDraws=size(x,1);
        B=NumberOfDraws;
    case 'wrong'
        x=[lpmat0(iwrong,:) lpmat(iwrong,:)];
        NumberOfDraws=size(x,1);
        B=NumberOfDraws;
        
end

if isempty(x),
    disp('Empty parameter set!')
    return
end

for j=1:length(varargin),
    jcol(j)=strmatch(varargin{j},bayestopt_.name,'exact');
end

[H,AX,BigA,P,PAx]=plotmatrix(x(:,jcol));

for j=1:length(varargin),
    %      axes(AX(1,j)), title(varargin{j})
    %      axes(AX(j,1)), ylabel(varargin{j})
    %      set(AX(1,j),'title',varargin{j}),
    set(get(AX(j,1),'ylabel'),'string',varargin{j})
    set(get(AX(end,j),'xlabel'),'string',varargin{j})
end

if options_.opt_gsa.pprior==0,
    xparam1=xparam1(jcol);
    for j=1:length(varargin),
        for i=1:j-1,
            axes(AX(j,i)),
            hold on, plot(xparam1(i),xparam1(j),'*r')
        end
        for i=j+1:length(varargin),
            axes(AX(j,i)),
            hold on, plot(xparam1(i),xparam1(j),'*r')
        end
    end
end
