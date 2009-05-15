
function fparallel(fblck,nblck,whoiam,ThisMatlab,fname);

global funcName

funcName=fname;

warning off;
diary off;

delete( [fname,'_',int2str(whoiam),'.log']);

diary( [fname,'_',int2str(whoiam),'.log']);

    
% configure dynare environment
dynareroot = dynare_config();

load( [fname,'_input']) 
% keyboard;
if exist('fGlobalVar'),
  globalVars = fieldnames(fGlobalVar);
  for j=1:length(globalVars),
    eval(['global ',globalVars{j},';'])
  end
  struct2local(fGlobalVar);
end

% lounch the routine to be run in parallel
tic,
[fOutputVar, OutputFileName] = feval(fname, fInputVar ,fblck, nblck, whoiam, ThisMatlab);
toc,

%%% Sincronismo "Esterno" %%%%%%%%%%%%%
%%% Ogni Processo quando ha finito lo notifica cancellando un file ... 
% keyboard;
if(whoiam)
  save([ fname,'_output_',int2str(whoiam)],'fOutputVar' )
  
  
  if options_.parallel(ThisMatlab).Local
    delete(['P_',fname,'_',int2str(whoiam),'End.txt']);

  else
    if isunix,
%     for j=1:size(OutputFileName,1),
%       copyfile([OutputFileName{j,1},OutputFileName{j,2}],['\\',fInputVar.MasterName,'\',fInputVar.DyMo(1),'$\',fInputVar.DyMo(4:end),'\',OutputFileName{j,1}])
%     end
%     copyfile([fname,'_output_',int2str(whoiam),'.mat'],['\\',fInputVar.MasterName,'\',fInputVar.DyMo(1),'$\',fInputVar.DyMo(4:end)]);
%     delete(['\\',fInputVar.MasterName,'\',fInputVar.DyMo(1),'$\',fInputVar.DyMo(4:end),'\P_',fname,'_',int2str(whoiam),'End.txt']);
    else
      for j=1:size(OutputFileName,1),
        copyfile([OutputFileName{j,1},OutputFileName{j,2}],['\\',fInputVar.MasterName,'\',fInputVar.DyMo(1),'$\',fInputVar.DyMo(4:end),'\',OutputFileName{j,1}])
      end
      copyfile([fname,'_output_',int2str(whoiam),'.mat'],['\\',fInputVar.MasterName,'\',fInputVar.DyMo(1),'$\',fInputVar.DyMo(4:end)]);
      delete(['\\',fInputVar.MasterName,'\',fInputVar.DyMo(1),'$\',fInputVar.DyMo(4:end),'\P_',fname,'_',int2str(whoiam),'End.txt']);
    end
  end
end

disp(['fParallel ',int2str(whoiam),' completed.'])
diary off;

exit;
