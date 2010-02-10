function myoutput=PosteriorIRF_core2(myinputs,fpar,npar,whoiam, ThisMatlab)


if nargin<4,
    whoiam=0;
end
struct2local(myinputs);

% To save the figures where the function is computed!

DirectoryName = CheckPath('Output');

RemoteFlag = 0;
if whoiam
      waitbarString = ['Please wait... PosteriorIRF Plots (exog. shocks ' int2str(fpar) 'of' int2str(npar) ')...'];
      if Parallel(ThisMatlab).Local,
        waitbarTitle=['Local '];
      else
        waitbarTitle=[Parallel(ThisMatlab).PcName];
        RemoteFlag = 1;
      end        
        fMessageStatus(0,whoiam,waitbarString, waitbarTitle, Parallel(ThisMatlab), MasterName, DyMo);   
 end

OutputFileName={};

subplotnum = 0;
for i=fpar:npar
  NAMES = [];
  figunumber = 0;
  
  for j=1:nvar
    if max(abs(MeanIRF(:,j,i))) > 10^(-6)
      subplotnum = subplotnum+1;
      if options_.nograph
        if subplotnum == 1 & options_.relative_irf
          hh = figure('Name',['Relative response to orthogonalized shock to ' tit(i,:)],'Visible','off');
        elseif subplotnum == 1 & ~options_.relative_irf
          hh = figure('Name',['Orthogonalized shock to ' tit(i,:)],'Visible','off');
        end
      else
        if subplotnum == 1 & options_.relative_irf
          hh = figure('Name',['Relative response to orthogonalized shock to ' tit(i,:)]);
        elseif subplotnum == 1 & ~options_.relative_irf
          hh = figure('Name',['Orthogonalized shock to ' tit(i,:)]);
        end
      end
      set(0,'CurrentFigure',hh)
      subplot(nn,nn,subplotnum);
      if ~MAX_nirfs_dsgevar
          h1 = area(1:options_.irf,HPDIRF(:,2,j,i));
          set(h1,'FaceColor',[.9 .9 .9]);
          set(h1,'BaseValue',min(HPDIRF(:,1,j,i)));
          hold on
          h2 = area(1:options_.irf,HPDIRF(:,1,j,i),'FaceColor',[1 1 1],'BaseValue',min(HPDIRF(:,1,j,i)));
          set(h2,'FaceColor',[1 1 1]);
          set(h2,'BaseValue',min(HPDIRF(:,1,j,i)));
          plot(1:options_.irf,MeanIRF(:,j,i),'-k','linewidth',3)
          % plot([1 options_.irf],[0 0],'-r','linewidth',0.5);          
          box on
          axis tight
          xlim([1 options_.irf]);
          hold off
      else    
          h1 = area(1:options_.irf,HPDIRF(:,2,j,i));
          set(h1,'FaceColor',[.9 .9 .9]);
          set(h1,'BaseValue',min([min(HPDIRF(:,1,j,i)),min(HPDIRFdsgevar(:,1,j,i))]));
          hold on;
          h2 = area(1:options_.irf,HPDIRF(:,1,j,i));
          set(h2,'FaceColor',[1 1 1]);
          set(h2,'BaseValue',min([min(HPDIRF(:,1,j,i)),min(HPDIRFdsgevar(:,1,j,i))]));
          plot(1:options_.irf,MeanIRF(:,j,i),'-k','linewidth',3)
          % plot([1 options_.irf],[0 0],'-r','linewidth',0.5);
          plot(1:options_.irf,MeanIRFdsgevar(:,j,i),'--k','linewidth',2)
          plot(1:options_.irf,HPDIRFdsgevar(:,1,j,i),'--k','linewidth',1)
          plot(1:options_.irf,HPDIRFdsgevar(:,2,j,i),'--k','linewidth',1)
          box on
          axis tight
          xlim([1 options_.irf]);
          hold off
      end
      name = deblank(varlist(j,:));
      NAMES = strvcat(NAMES,name);
      title(name,'Interpreter','none')
    end
    
    
%      TempPath=DirectoryName;
%      DirectoryNamePar='C:\dynare_calcs\ModelTest\ls2003\metropolis';
%      DirectoryName=DirectoryNamePar;

    
    if subplotnum == MaxNumberOfPlotPerFigure | (j == nvar  & subplotnum> 0)
      figunumber = figunumber+1;
      set(hh,'visible','on')
      eval(['print -depsc2 ' DirectoryName '/'  M_.fname '_Bayesian_IRF_' deblank(tit(i,:)) '_' int2str(figunumber) '.eps']);
      if ~exist('OCTAVE_VERSION')
          eval(['print -dpdf ' DirectoryName '/' M_.fname  '_Bayesian_IRF_' deblank(tit(i,:)) '_' int2str(figunumber)]);
          saveas(hh,[DirectoryName '/' M_.fname  '_Bayesian_IRF_' deblank(tit(i,:))  '_' int2str(figunumber) '.fig']);
      end
      if RemoteFlag==1,
      OutputFileName = [OutputFileName; {[DirectoryName,filesep], [M_.fname '_Bayesian_IRF_' deblank(tit(i,:)) '_' int2str(figunumber) '.*']}];
      end
      set(hh,'visible','off')
      if options_.nograph, close(hh), end
      subplotnum = 0;
    end
  end% loop over selected endo_var
    if whoiam,
           fprintf('Done! \n');
           waitbarString = [ 'Exog. shocks ' int2str(i) '/' int2str(npar) ' done.'];
           fMessageStatus((i-fpar+1)/(npar-fpar+1),whoiam,waitbarString, waitbarTitle, Parallel(ThisMatlab), MasterName, DyMo)
    end
end% loop over exo_var  
 

% DirectoryName=TempPath;
 
myoutput.OutputFileName = OutputFileName;

