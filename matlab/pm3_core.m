function myoutput=pm3_core(myinputs,fpar,nvar,whoiam, ThisMatlab)


if nargin<4,
    whoiam=0;
end
struct2local(myinputs);


global options_ M_ oo_


if whoiam
      waitbarString = ['Parallel plots pm3 ...'];
      if Parallel(ThisMatlab).Local,
        waitbarTitle=['Local '];
      else
        waitbarTitle=[Parallel(ThisMatlab).PcName];
      end        
        fMessageStatus(0,whoiam,waitbarString, waitbarTitle, Parallel(ThisMatlab), MasterName, DyMo);   
 end



figunumber = 0;
subplotnum = 0;
hh = figure('Name',[tit1 ' ' int2str(figunumber+1)]);
RemoteFlag = 0;
if whoiam,
    if Parallel(ThisMatlab).Local ==0
        RemoteFlag=1;
    end
end

OutputFileName = {};

for i=fpar:nvar
    NAMES = [];
    if max(abs(Mean(:,i))) > 10^(-6)
        subplotnum = subplotnum+1;
        set(0,'CurrentFigure',hh)
        subplot(nn,nn,subplotnum);
        plot([1 n2],[0 0],'-r','linewidth',0.5);
        hold on
        for k = 1:9
            plot(1:n2,squeeze(Distrib(k,:,i)),'-g','linewidth',0.5)
        end
        plot(1:n2,Mean(:,i),'-k','linewidth',1)
        xlim([1 n2]);
        hold off
        name = deblank(varlist(i,:));
        NAMES = strvcat(NAMES,name);
        title(name,'Interpreter','none')
    end
    if subplotnum == MaxNumberOfPlotsPerFigure | i == nvar
        eval(['print -depsc2 ' M_.dname '/Output/'  M_.fname '_' name3 '_' deblank(tit3(i,:)) '.eps' ]);
        if ~exist('OCTAVE_VERSION')
            eval(['print -dpdf ' M_.dname '/Output/' M_.fname  '_' name3 '_' deblank(tit3(i,:))]);
            saveas(hh,[M_.dname '/Output/' M_.fname '_' name3 '_' deblank(tit3(i,:)) '.fig']);
        end
        if RemoteFlag==1,
            OutputFileName = [OutputFileName; {[M_.dname, filesep, 'Output',filesep], [M_.fname '_' name3 '_' deblank(tit3(i,:)) '.*']}];
        end
        if options_.nograph, close(hh), end
        subplotnum = 0;
        figunumber = figunumber+1;
        if (i ~= nvar)
            hh = figure('Name',[name3 ' ' int2str(figunumber+1)]);
        end
    end
    
    if whoiam,
        waitbarString = [ 'Variable ' int2str(i) '/' int2str(nvar) ' done.'];
        fMessageStatus((i-fpar+1)/(nvar-fpar+1),whoiam,waitbarString, waitbarTitle, Parallel(ThisMatlab), MasterName, DyMo)
    end
    
    
end


myoutput.OutputFileName=OutputFileName;
