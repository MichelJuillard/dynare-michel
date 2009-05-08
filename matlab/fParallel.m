
function fparallel(fblck,nblck,whoiam,ThisMatlab,fname);

    
% configure dynare environment
dynareroot = dynare_config();

load( [fname,'_input']) 
% keyboard;
global fname
if exist('fGlobalVar'),
  globalVars = fieldnames(fGlobalVar);
  for j=1:length(globalVars),
    eval(['global ',globalVars{j},';'])
  end
  struct2local(fGlobalVar);
end

% lounch the routine to be run in parallel
fOutputVar = feval(fname, fInputVar ,fblck, nblck, whoiam, ThisMatlab);


%%% Sincronismo "Esterno" %%%%%%%%%%%%%
%%% Ogni Processo quando ha finito lo notifica cancellando un file ... 
% keyboard;
if(whoiam)
  save([ fname,'_output_',int2str(whoiam)],'fOutputVar' )
  
  
  if options_.parallel(ThisMatlab).Local
    delete(['P_',fname,'_',int2str(whoiam),'End.txt']);
  else
    delete(['\\',fInputVar.MasterName,'\',fInputVar.DyMo(1),'$\',fInputVar.DyMo(4:end),'\P_',fname,'_',int2str(whoiam),'End.txt']);
   
  end
end



exit;
