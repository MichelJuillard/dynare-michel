function slaveParallel(whoiam,ThisMatlab)
% PARALLEL CONTEXT
% In a parallelization context, this function is launched on slave
% machines, to initialize MATLAB and DYNARE environment and waits for
% instructions sent by the Master. 
% This function is invoked by masterParallel only when the strategy (1),
% i.e. always open, is actived.
%
%
% INPUTS
%  o whoiam [int]         index number of this CPU among all CPUs in the
%                       cluster.
%  o ThisMatlab [int]     index number of this slave machine in the cluster.
%
% OUTPUTS 
%   None  

% Copyright (C) 2006-2008,2010 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global funcName

warning off;
diary off;

delete( ['slaveParallel_',int2str(whoiam),'.log']);

diary( ['slaveParallel_',int2str(whoiam),'.log']);

% Load input data.
load( ['slaveParallel_input',int2str(whoiam)])

% configure dynare environment
if ~empty(Parallel(ThisMatlab).DynarePath)
    addpath(Parallel(ThisMatlab).DynarePath)
end
dynareroot = dynare_config();

%Loads fGlobalVar Parallel.
if exist('fGlobalVar'),
    globalVars = fieldnames(fGlobalVar);
    for j=1:length(globalVars),
        eval(['global ',globalVars{j},';'])
        evalin('base',['global ',globalVars{j},';'])
    end
    struct2local(fGlobalVar);
    clear fGlobalVar
    % create global variables in the base workspace as well
    evalin('base',['load( [''slaveParallel_input',int2str(whoiam),'''],''fGlobalVar'')']) 
    evalin('base','struct2local(fGlobalVar)');
    evalin('base','clear fGlobalVar');
end

t0=clock;
fslave = dir( ['slaveParallel_input',int2str(whoiam),'.mat']);

while (etime(clock,t0)<1200 && ~isempty(fslave)) || ~isempty(dir(['stayalive',int2str(whoiam),'.txt'])),
    if ~isempty(dir(['stayalive',int2str(whoiam),'.txt'])),
        t0=clock;
        delete(['stayalive',int2str(whoiam),'.txt']);
    end
    % I wait for 20 min or while mater asks to exit (i.e. it cancels fslave file)
    pause(1);
    
    % -> Da Sistemare!!!!!!!!!!!!!!!!
    % Con testing su reti vere e con core reali!!!!
    
    fjob = dir(['slaveJob',int2str(whoiam),'.mat']);
    
    if ~isempty(fjob),
        clear fGlobalVar fInputVar fblck nblck fname
        
        while(1)
           Go=0;
           
           Go=fopen(['slaveJob',int2str(whoiam),'.mat']);
          
           if Go>0    
               fclose(Go);
               pause(1);
               load(['slaveJob',int2str(whoiam),'.mat']);
               break
           else
               % Only for testing, will be remouved!
               
               if isunx
                 E1=fopen('/home/ivano/Works/Errore-slaveParallel.txt','w+');
                 fclose(E1);
               else            
                 E1=fopen('c:\dynare_calcs\Errore-slaveParallel.txt','w+');
                 fclose(E1);
               end
                       
           end
         end
               
        funcName=fname;  % Update global job name.

        if exist('fGlobalVar') && ~isempty (fGlobalVar)
            globalVars = fieldnames(fGlobalVar);
            for j=1:length(globalVars),
                info_whos = whos(globalVars{j});
                if isempty(info_whos) || ~info_whos.global,
                    eval(['global ',globalVars{j},';'])
                    evalin('base',['global ',globalVars{j},';'])
                end
            end
            struct2local(fGlobalVar);
            evalin('base',['load( [''slaveJob',int2str(whoiam),'''],''fGlobalVar'')']) 
            evalin('base','struct2local(fGlobalVar)');
            evalin('base','clear fGlobalVar');
        end
        delete(['slaveJob',int2str(whoiam),'.mat']);
        fInputVar.Parallel = Parallel;
        
        % Launch the routine to be run in parallel.
        tic,
        fOutputVar = feval(fname, fInputVar ,fblck, nblck, whoiam, ThisMatlab);
        toc,
        if isfield(fOutputVar,'OutputFileName'),
            OutputFileName = fOutputVar.OutputFileName;
        else
            OutputFileName = '';
        end

        if(whoiam)

            % Save the output result.
            save([ fname,'_output_',int2str(whoiam),'.mat'],'fOutputVar' )

            % Inform the master that the job is finished, and transfer the output data
            if Parallel(ThisMatlab).Local
                delete(['P_',fname,'_',int2str(whoiam),'End.txt']);
            else
                if isunix || (~matlab_ver_less_than('7.4') && ismac),
                    for j=1:size(OutputFileName,1),
                        system(['scp ',OutputFileName{j,1},OutputFileName{j,2},' ',Parallel(ThisMatlab).user,'@',fInputVar.MasterName,':',fInputVar.DyMo,'/',OutputFileName{j,1}]);
                    end
                    system(['scp ',fname,'_output_',int2str(whoiam),'.mat ',Parallel(ThisMatlab).user,'@',fInputVar.MasterName,':',fInputVar.DyMo]);
                    system(['ssh ',Parallel(ThisMatlab).user,'@',fInputVar.MasterName,' rm -f ',fInputVar.DyMo,'/P_',fname,'_',int2str(whoiam),'End.txt']);
                    %                 system(['fusermount -u ~/MasterRemoteMirror_',fname,'_',int2str(whoiam)]);
                    %                 system(['rm -r ~/MasterRemoteMirror_',fname,'_',int2str(whoiam)]);
                else
                    for j=1:size(OutputFileName,1),
                        copyfile([OutputFileName{j,1},OutputFileName{j,2}],['\\',fInputVar.MasterName,'\',fInputVar.DyMo(1),'$\',fInputVar.DyMo(4:end),'\',OutputFileName{j,1}])
                    end
                    copyfile([fname,'_output_',int2str(whoiam),'.mat'],['\\',fInputVar.MasterName,'\',fInputVar.DyMo(1),'$\',fInputVar.DyMo(4:end)]);
                    delete(['\\',fInputVar.MasterName,'\',fInputVar.DyMo(1),'$\',fInputVar.DyMo(4:end),'\P_',fname,'_',int2str(whoiam),'End.txt']);
                end
            end
        end

        disp(['Job ',fname,' on CPU ',int2str(whoiam),' completed.'])
        t0 =clock; % re-set waiting time of 20 mins
    end
    fslave = dir( ['slaveParallel_input',int2str(whoiam),'.mat']); % check if Master asks to exit
end


if Parallel(ThisMatlab).Local
    delete(['P_slave_',int2str(whoiam),'End.txt']);
else
    if isunix || (~matlab_ver_less_than('7.4') && ismac),
        system(['ssh ',Parallel(ThisMatlab).user,'@',fInputVar.MasterName,' rm -f ',fInputVar.DyMo,'/P_slave_',int2str(whoiam),'End.txt']);
    else
        delete(['\\',fInputVar.MasterName,'\',fInputVar.DyMo(1),'$\',fInputVar.DyMo(4:end),'\P_slave_',int2str(whoiam),'End.txt']);
    end
end
disp(['slaveParallel on CPU ',int2str(whoiam),' completed.'])
diary off;

exit;
