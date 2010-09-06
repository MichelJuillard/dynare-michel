function fParallel(fblck,nblck,whoiam,ThisMatlab,fname)
% PARALLEL CONTEXT
% In a parallel context, this function is launched on slave
% machines, and acts as a wrapper around the function containing the
% computing task itself.
%
% INPUTS
%  o fblck [int]          index number of the first thread to run in this
%                       MATLAB instance
%  o nblck [int]          number of threads to run in this
%                       MATLAB instance
%  o whoiam [int]         index number of this CPU among all CPUs in the
%                       cluster
%  o ThisMatlab [int]     index number of this slave machine in the cluster
%                       (entry in options_.parallel)
%  o fname [string]       function to be run, containing the computing task
%
% OUTPUTS 
%   None 
%
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

funcName=fname;

warning off;
diary off;

delete( [fname,'_',int2str(whoiam),'.log']);

diary( [fname,'_',int2str(whoiam),'.log']);

% Load input data
load( [fname,'_input']) 

% configure dynare environment
if ~empty(Parallel(ThisMatlab).DynarePath)
    addpath(Parallel(ThisMatlab).DynarePath)
end
dynareroot = dynare_config();

if exist('fGlobalVar') && ~isempty (fGlobalVar)
    globalVars = fieldnames(fGlobalVar);
    for j=1:length(globalVars),
        eval(['global ',globalVars{j},';'])
        evalin('base',['global ',globalVars{j},';'])
    end
    struct2local(fGlobalVar);
    % create global variables in the base workspace as well
    evalin('base',['load( [''',fname,'_input''],''fGlobalVar'')']) 
    evalin('base','struct2local(fGlobalVar)');
end

% On UNIX, mount the master working directory through SSH FS
% if (isunix | (~matlab_ver_less_than('7.4') & ismac) ) & Parallel(ThisMatlab).Local==0,
%     system(['mkdir ~/MasterRemoteMirror_',fname,'_',int2str(whoiam)]);
%     system(['sshfs ',Parallel(ThisMatlab).user,'@',fInputVar.MasterName,':/',fInputVar.DyMo,' ~/MasterRemoteMirror_',fname,'_',int2str(whoiam)]);
% end

% Special hack for MH directory: possibly no longer needed ?????? now all
% files are managed locally and sent backwards later on in this routine!
% if isfield(fInputVar,'MhDirectoryName') & Parallel(ThisMatlab).Local==0,
%     if isunix | (~matlab_ver_less_than('7.4') & ismac),
%         fInputVar.MhDirectoryName = ['~/MasterRemoteMirror_',fname,'_',int2str(whoiam),'/',fInputVar.MhDirectoryName];
%     else
%         fInputVar.MhDirectoryName = ['\\',fInputVar.MasterName,'\',fInputVar.DyMo(1),'$\',fInputVar.DyMo(4:end),'\',fInputVar.MhDirectoryName];
%     end
% end

fInputVar.Parallel = Parallel;
% Launch the routine to be run in parallel
tic,
% keyboard;
fOutputVar = feval(fname, fInputVar ,fblck, nblck, whoiam, ThisMatlab);
toc,
if isfield(fOutputVar,'OutputFileName'),
    OutputFileName = fOutputVar.OutputFileName;
else
    OutputFileName = '';
end

if(whoiam)
    
    % Save the output result
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
%             system(['fusermount -u ~/MasterRemoteMirror_',fname,'_',int2str(whoiam)]);
%             system(['rm -r ~/MasterRemoteMirror_',fname,'_',int2str(whoiam)]);      
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
